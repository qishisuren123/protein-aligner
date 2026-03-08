"""
模块4: Alignment & Quality Control - 对齐与质量控制

功能：
- 验证原子模型与密度图的坐标系对齐
- 使用 gemmi 进行坐标插值（确保坐标系正确）
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


class AlignmentQC:
    """对齐验证与质量控制"""

    def __init__(self, config):
        self.config = config
        self.sim_resolution = config['alignment_qc']['sim_resolution']
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

    def compute_cc_at_atoms(self, grid, structure, radius=2.0):
        """
        计算原子位置处的密度相关性（使用 gemmi 插值）

        更可靠的 CC 计算方法：直接在原子位置采样密度值，
        与理论值（1.0 对于有原子的位置）比较

        参数:
            grid: gemmi Grid 对象（归一化后的密度图）
            structure: gemmi Structure 对象
            radius: 采样半径 (Å)
        返回:
            cc_value, 诊断信息
        """
        model = structure[0]
        exp_values = []  # 实验密度值（在原子位置）
        bg_values = []   # 背景密度值

        # 在 CA 原子位置采样实验密度
        ca_coords = []
        for chain in model:
            for residue in chain:
                if residue.name == "HOH":
                    continue
                ca = residue.find_atom("CA", '*')
                if ca:
                    val = grid.interpolate_value(ca.pos)
                    exp_values.append(val)
                    ca_coords.append(ca.pos)

        exp_values = np.array(exp_values)

        if len(exp_values) == 0:
            return 0.0, {"n_ca": 0}

        # 统计诊断信息
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

        # 生成模拟密度图并计算 CC
        # 方法：在同一组位置采样模拟密度和实验密度
        # 模拟密度：CA 位置 = 1.0，远离原子 = 0.0
        sim_values = np.ones(len(exp_values))  # CA 位置的理论密度

        # CC 只在有信号的区域计算
        if exp_values.std() == 0:
            return 0.0, info

        cc, _ = pearsonr(exp_values, sim_values)
        # 这个 CC 不太有意义因为 sim_values 全是 1
        # 改用: 正值密度比例作为质量指标
        quality_score = (exp_values > 0).mean()
        info["quality_score"] = float(quality_score)

        return quality_score, info

    def generate_simulated_map(self, grid, structure, resolution=None):
        """
        从原子模型生成模拟密度图（使用 gemmi 坐标系）

        参数:
            grid: 参考 gemmi Grid（用于确定网格大小和坐标系）
            structure: gemmi Structure
            resolution: 模拟分辨率 (Å)
        返回:
            sim_grid: gemmi FloatGrid（模拟密度图）
        """
        if resolution is None:
            resolution = self.sim_resolution

        # 创建与输入 grid 同形状的模拟 grid
        sim_grid = gemmi.FloatGrid(grid.nu, grid.nv, grid.nw)
        sim_grid.set_unit_cell(grid.unit_cell)
        sim_grid.spacegroup = grid.spacegroup
        sim_data = np.array(sim_grid, copy=False)

        model = structure[0]
        n_atoms = 0
        for chain in model:
            for residue in chain:
                if residue.name == "HOH":
                    continue
                for atom in residue:
                    # 将原子坐标转换为 grid 索引
                    frac = grid.unit_cell.fractionalize(atom.pos)
                    i = int(round(frac.x * grid.nu)) % grid.nu
                    j = int(round(frac.y * grid.nv)) % grid.nv
                    k = int(round(frac.z * grid.nw)) % grid.nw
                    sim_data[i, j, k] += 1.0
                    n_atoms += 1

        logger.info(f"  放置了 {n_atoms} 个原子到模拟密度图")

        # 高斯模糊
        sigma_A = resolution / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        sigma_voxels = sigma_A / grid.spacing[0]
        sim_data_blurred = gaussian_filter(sim_data, sigma=sigma_voxels)

        if sim_data_blurred.max() > 0:
            sim_data_blurred /= sim_data_blurred.max()

        sim_data[:] = sim_data_blurred
        return sim_grid

    def compute_cc_mask(self, exp_grid, sim_grid):
        """计算 masked CC (在模拟密度>0的区域)"""
        exp_data = np.array(exp_grid, copy=True)
        sim_data = np.array(sim_grid, copy=True)

        mask = sim_data > 0.01
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

        # 使用 gemmi 加载密度图
        ccp4 = gemmi.read_ccp4_map(map_path)
        ccp4.setup(float('nan'))
        grid = ccp4.grid

        logger.info(f"  Grid: ({grid.nu},{grid.nv},{grid.nw}), "
                    f"spacing: ({grid.spacing[0]:.3f},{grid.spacing[1]:.3f},{grid.spacing[2]:.3f})")

        # 加载结构
        structure = self.load_structure(model_path)
        atoms = self.get_atom_coords(structure)
        logger.info(f"  加载了 {len(atoms)} 个原子")

        # 在原子位置采样密度并评估
        quality_score, diag_info = self.compute_cc_at_atoms(grid, structure)
        logger.info(f"  质量分数 (正密度比例): {quality_score:.4f}")

        # 生成模拟密度图
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
            "spacing": [float(grid.spacing[0]), float(grid.spacing[1]), float(grid.spacing[2])],
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

                result = {"entry_dir": entry_dir, "metrics": metrics}
                results.append(result)

                if metrics["passed"]:
                    passed.append(entry_dir)
                else:
                    failed.append(entry_dir)
            except Exception as e:
                logger.error(f"  QC 失败: {e}", exc_info=True)

        logger.info(f"质量控制完成: {len(passed)} 通过, {len(failed)} 未通过")
        return results, passed
