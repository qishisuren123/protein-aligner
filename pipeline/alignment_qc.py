"""
模块4: Alignment & Quality Control - 对齐与质量控制

功能：
- 验证原子模型与密度图的坐标系对齐
- 使用 gemmi 进行坐标插值（确保坐标系正确）
- 使用三线性插值放置原子生成高精度模拟密度图
- 计算 map-model 拟合质量指标（CC_mask 等）
- 根据质量阈值过滤低质量数据
"""
import os
import json
import logging
import numpy as np
import gemmi
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)

# 元素散射权重（归一化到碳=1.0）
ELEMENT_WEIGHT = {
    'C': 1.0, 'N': 1.17, 'O': 1.33, 'S': 2.67, 'P': 2.5,
    'H': 0.17, 'FE': 4.33, 'ZN': 5.0, 'MG': 2.0, 'CA': 3.33,
    'NA': 1.83, 'CL': 2.83, 'MN': 4.17, 'CO': 4.5, 'CU': 4.83,
}


class AlignmentQC:
    """对齐验证与质量控制"""

    def __init__(self, config):
        self.config = config
        self.sim_sigma = config['alignment_qc'].get('sim_sigma', 0.8)
        self.mask_threshold = config['alignment_qc'].get('mask_threshold', 0.05)
        self.cc_threshold = config['alignment_qc']['cc_threshold']

    def load_structure(self, cif_path):
        """加载 mmCIF/PDB 原子模型"""
        st = gemmi.read_structure(cif_path)
        st.setup_entities()
        return st

    def get_atom_coords(self, structure):
        """从结构中提取所有原子坐标和信息"""
        atoms = []
        model = structure[0]
        for chain in model:
            for residue in chain:
                if residue.name == "HOH":
                    continue
                for atom in residue:
                    atoms.append({
                        'pos': atom.pos,
                        'coord': np.array([atom.pos.x, atom.pos.y, atom.pos.z]),
                        'name': atom.name,
                        'element': atom.element.name,
                        'residue_name': residue.name,
                        'residue_seq': residue.seqid.num,
                        'chain': chain.name,
                        'b_factor': atom.b_iso,
                    })
        return atoms

    def compute_cc_at_atoms(self, grid, structure):
        """
        在 CA 原子位置采样密度并统计诊断信息

        参数:
            grid: gemmi Grid 对象（归一化后的密度图）
            structure: gemmi Structure 对象
        返回:
            quality_score (正值密度比例), 诊断信息
        """
        model = structure[0]
        exp_values = []

        for chain in model:
            for residue in chain:
                if residue.name == "HOH":
                    continue
                ca = residue.find_atom("CA", '*')
                if ca:
                    val = grid.interpolate_value(ca.pos)
                    exp_values.append(val)

        exp_values = np.array(exp_values)

        if len(exp_values) == 0:
            return 0.0, {"n_ca": 0}

        info = {
            "n_ca": len(exp_values),
            "density_mean_at_ca": float(exp_values.mean()),
            "density_median_at_ca": float(np.median(exp_values)),
            "density_positive_ratio": float((exp_values > 0).mean()),
            "density_above_0.1_ratio": float((exp_values > 0.1).mean()),
        }

        logger.info(f"  CA 原子处密度: mean={info['density_mean_at_ca']:.4f}, "
                    f"median={info['density_median_at_ca']:.4f}, "
                    f">0 比例={info['density_positive_ratio']*100:.1f}%")

        quality_score = float((exp_values > 0).mean())
        info["quality_score"] = quality_score

        return quality_score, info

    def _place_atoms_trilinear(self, model, cell, nu, nv, nw):
        """
        三线性插值放置原子到 grid 上（保留亚体素精度）

        将每个原子的权重分配到最近的 8 个体素，
        按距离加权，避免四舍五入到最近体素的精度损失
        """
        sim = np.zeros((nu, nv, nw), dtype=np.float32)

        for chain in model:
            for residue in chain:
                if residue.name == "HOH":
                    continue
                for atom in residue:
                    frac = cell.fractionalize(atom.pos)
                    # 连续坐标
                    fi = frac.x * nu
                    fj = frac.y * nv
                    fk = frac.z * nw
                    # 下界索引
                    i0 = int(np.floor(fi)) % nu
                    j0 = int(np.floor(fj)) % nv
                    k0 = int(np.floor(fk)) % nw
                    i1 = (i0 + 1) % nu
                    j1 = (j0 + 1) % nv
                    k1 = (k0 + 1) % nw
                    # 小数部分
                    di = fi - np.floor(fi)
                    dj = fj - np.floor(fj)
                    dk = fk - np.floor(fk)
                    # 元素权重
                    w = ELEMENT_WEIGHT.get(atom.element.name, 1.0)
                    # 分配到 8 个相邻体素
                    sim[i0, j0, k0] += w * (1-di) * (1-dj) * (1-dk)
                    sim[i1, j0, k0] += w * di     * (1-dj) * (1-dk)
                    sim[i0, j1, k0] += w * (1-di) * dj     * (1-dk)
                    sim[i0, j0, k1] += w * (1-di) * (1-dj) * dk
                    sim[i1, j1, k0] += w * di     * dj     * (1-dk)
                    sim[i1, j0, k1] += w * di     * (1-dj) * dk
                    sim[i0, j1, k1] += w * (1-di) * dj     * dk
                    sim[i1, j1, k1] += w * di     * dj     * dk

        return sim

    def generate_simulated_map(self, grid, structure, sigma_A=None):
        """
        从原子模型生成模拟密度图

        使用三线性插值放置原子 + 高斯模糊模拟电子密度
        """
        if sigma_A is None:
            sigma_A = self.sim_sigma

        cell = grid.unit_cell
        nu, nv, nw = grid.nu, grid.nv, grid.nw
        model = structure[0]

        # 三线性插值放置原子
        sim_data = self._place_atoms_trilinear(model, cell, nu, nv, nw)
        n_atoms = sum(
            1 for chain in model for res in chain
            if res.name != "HOH" for _ in res
        )
        logger.info(f"  放置了 {n_atoms} 个原子（三线性插值）")

        # 高斯模糊
        sigma_voxels = sigma_A / grid.spacing[0]
        sim_data = gaussian_filter(sim_data, sigma=sigma_voxels)

        if sim_data.max() > 0:
            sim_data /= sim_data.max()

        # 包装为 gemmi grid
        sim_grid = gemmi.FloatGrid(nu, nv, nw)
        sim_grid.set_unit_cell(cell)
        sim_grid.spacegroup = grid.spacegroup
        arr = np.array(sim_grid, copy=False)
        arr[:] = sim_data

        return sim_grid

    def compute_cc_mask(self, exp_grid, sim_grid, mask_threshold=None):
        """
        计算 masked CC

        mask 为模拟密度图中信号强度超过阈值的区域
        """
        if mask_threshold is None:
            mask_threshold = self.mask_threshold

        exp_data = np.array(exp_grid, copy=True)
        sim_data = np.array(sim_grid, copy=True)

        mask = sim_data > mask_threshold
        if mask.sum() == 0:
            return 0.0

        exp_masked = exp_data[mask]
        sim_masked = sim_data[mask]

        if exp_masked.std() == 0 or sim_masked.std() == 0:
            return 0.0

        cc, _ = pearsonr(exp_masked.flatten(), sim_masked.flatten())
        return cc

    def evaluate_entry(self, entry_dir):
        """评估单个条目的对齐质量"""
        map_path = os.path.join(entry_dir, "map_normalized.mrc")
        model_path = os.path.join(entry_dir, "model.cif")

        if not os.path.exists(map_path) or not os.path.exists(model_path):
            logger.warning(f"跳过 {entry_dir}: 缺少文件")
            return None

        # 加载密度图
        ccp4 = gemmi.read_ccp4_map(map_path)
        ccp4.setup(float('nan'))
        grid = ccp4.grid

        logger.info(f"  Grid: ({grid.nu},{grid.nv},{grid.nw}), "
                    f"spacing: ({grid.spacing[0]:.3f},{grid.spacing[1]:.3f},{grid.spacing[2]:.3f})")

        # 加载结构
        structure = self.load_structure(model_path)
        atoms = self.get_atom_coords(structure)
        logger.info(f"  加载了 {len(atoms)} 个原子")

        # 在原子位置采样密度
        quality_score, diag_info = self.compute_cc_at_atoms(grid, structure)
        logger.info(f"  质量分数 (正密度比例): {quality_score:.4f}")

        # 生成模拟密度图（三线性 + 高斯模糊）
        sim_grid = self.generate_simulated_map(grid, structure)

        # 计算 CC_mask
        cc_mask = self.compute_cc_mask(grid, sim_grid)
        logger.info(f"  CC_mask = {cc_mask:.4f}")

        # 判断质量
        passed = quality_score >= 0.5 or cc_mask >= self.cc_threshold
        quality = "pass" if passed else "fail"
        logger.info(f"  质量判定: {quality}")

        metrics = {
            "cc_mask": float(cc_mask),
            "quality_score": float(quality_score),
            "n_atoms": len(atoms),
            "grid_shape": [grid.nu, grid.nv, grid.nw],
            "spacing": [float(s) for s in grid.spacing],
            "passed": bool(passed),
            "quality": quality,
            "diagnostics": diag_info,
        }

        # 保存模拟密度图
        sim_path = os.path.join(entry_dir, "sim_map.mrc")
        sim_ccp4 = gemmi.Ccp4Map()
        sim_ccp4.grid = sim_grid
        sim_ccp4.update_ccp4_header()
        sim_ccp4.write_ccp4_map(sim_path)

        # 保存质量指标
        qc_path = os.path.join(entry_dir, "qc_metrics.json")
        with open(qc_path, 'w') as f:
            json.dump(metrics, f, indent=2)

        return metrics

    def run(self, entry_dirs):
        """对所有条目进行对齐验证和质量控制"""
        results = []
        passed = []
        failed = []

        for entry_dir in entry_dirs:
            logger.info(f"质量控制: {entry_dir}")
            try:
                metrics = self.evaluate_entry(entry_dir)
                if metrics is None:
                    continue
                results.append({"entry_dir": entry_dir, "metrics": metrics})
                if metrics["passed"]:
                    passed.append(entry_dir)
                else:
                    failed.append(entry_dir)
            except Exception as e:
                logger.error(f"  QC 失败: {e}", exc_info=True)

        logger.info(f"质量控制完成: {len(passed)} 通过, {len(failed)} 未通过")
        return results, passed
