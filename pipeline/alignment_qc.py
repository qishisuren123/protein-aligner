"""
模块4: Alignment & Quality Control - 对齐与质量控制

功能：
- 验证原子模型与密度图的坐标系对齐
- 按 ChimeraX molmap 规范生成模拟密度图
- 支持生物学组装体展开（对称蛋白完整覆盖密度）
- 计算 CC_box、CC_mask、CC_volume、CC_peaks 四项指标
  （Afonine et al., Acta Cryst. D74, 531–544, 2018）
- 根据质量阈值过滤低质量数据

CC 指标定义 (Afonine 2018):
- CC_box:    整个 box 内的 Pearson CC
- CC_mask:   分子掩膜内的 Pearson CC（掩膜由原子模型定义）
- CC_volume: 模拟密度最高 N 个体素的 Pearson CC（N = 掩膜内体素数）
- CC_peaks:  模拟 + 实验密度各自最高 N 个体素的并集的 Pearson CC
"""
import os
import json
import logging
import numpy as np
import gemmi
from scipy.spatial import cKDTree

from pipeline.bio_assembly import expand_to_assembly
from pipeline.qscore import compute_qscore_per_atom, compute_qscore_summary

logger = logging.getLogger(__name__)


def _pearson_cc(x, y):
    """
    计算 Pearson 相关系数（先减均值）

    等价于 scipy.stats.pearsonr 但更轻量，避免额外的 p-value 计算。
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if len(x) < 2:
        return 0.0
    xm = x - x.mean()
    ym = y - y.mean()
    denom = np.sqrt(np.sum(xm ** 2) * np.sum(ym ** 2))
    if denom < 1e-15:
        return 0.0
    return float(np.sum(xm * ym) / denom)


class AlignmentQC:
    """对齐验证与质量控制"""

    def __init__(self, config):
        self.config = config
        qc_cfg = config['alignment_qc']
        self.cc_threshold = qc_cfg.get('cc_threshold', 0.3)
        self.density_calculator = qc_cfg.get('density_calculator', 'gemmi_dce')
        self.default_resolution = qc_cfg.get('default_resolution', 2.0)
        self.dc_blur = qc_cfg.get('dc_blur', 0)
        # 分子掩膜半径（Å），默认根据分辨率自动设定
        self.mask_radius = qc_cfg.get('mask_radius', None)
        # 生物组装体配置
        bio_cfg = config.get('bio_assembly', {})
        self.bio_assembly_enabled = bio_cfg.get('enabled', True)
        self.assembly_id = bio_cfg.get('assembly_id', '1')
        self.max_chains = bio_cfg.get('max_chains', 200)

    def load_structure(self, cif_path):
        """加载 mmCIF/PDB 原子模型"""
        st = gemmi.read_structure(cif_path)
        st.setup_entities()
        return st

    def _get_resolution(self, entry_dir):
        """从 metadata.json 读取分辨率，无则使用默认值"""
        meta_path = os.path.join(entry_dir, "metadata.json")
        if os.path.exists(meta_path):
            with open(meta_path, 'r') as f:
                meta = json.load(f)
            res = meta.get("resolution")
            if res is not None:
                return float(res)
        logger.info(f"  未找到分辨率信息，使用默认值 {self.default_resolution}Å")
        return self.default_resolution

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

    def generate_simulated_map(self, grid, structure, resolution=None):
        """
        按 ChimeraX molmap 规范生成模拟密度图

        使用统一的 molmap 参数（sigma = resolution/(pi*sqrt(2))），
        保证 sim_map 和 mol_map 输出一致。

        参考: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/molmap.html
        """
        if resolution is None:
            resolution = self.default_resolution

        from pipeline.molmap import generate_molmap
        sim_data = generate_molmap(grid, structure, resolution)

        # 包装为 gemmi grid
        cell = grid.unit_cell
        sim_grid = gemmi.FloatGrid(grid.nu, grid.nv, grid.nw)
        sim_grid.set_unit_cell(cell)
        sim_grid.spacegroup = grid.spacegroup
        arr = np.array(sim_grid, copy=False)
        arr[:] = sim_data

        return sim_grid

    def _build_molecular_mask(self, grid, structure, resolution):
        """
        构建分子掩膜（Jiang & Brünger, 1994）

        以原子模型定义 mask：每个网格点到最近原子的距离 <= mask_radius
        则属于分子内部。

        mask_radius 默认 = max(3.0, resolution * 0.7) Å，
        这是 phenix/CCP4 常用的 mask 生成标准。

        参数:
            grid: gemmi.FloatGrid 参考网格
            structure: gemmi.Structure 原子模型（已展开组装体）
            resolution: float 分辨率 (Å)
        返回:
            mask: np.ndarray(bool), shape=(nu, nv, nw)
        """
        # 确定 mask 半径
        if self.mask_radius is not None:
            mask_r = self.mask_radius
        else:
            mask_r = max(3.0, resolution * 0.7)

        # 收集原子坐标
        model = structure[0]
        positions = []
        for chain in model:
            for res in chain:
                if res.name == "HOH":
                    continue
                for atom in res:
                    positions.append([atom.pos.x, atom.pos.y, atom.pos.z])

        if not positions:
            nu, nv, nw = grid.nu, grid.nv, grid.nw
            return np.zeros((nu, nv, nw), dtype=bool)

        positions = np.array(positions)
        atom_tree = cKDTree(positions)

        cell = grid.unit_cell
        nu, nv, nw = grid.nu, grid.nv, grid.nw

        # 生成网格点的笛卡尔坐标
        # 分数坐标 → 笛卡尔坐标
        orth_mat = np.array(cell.orth.mat.tolist())

        # 分块处理避免内存溢出（对大网格）
        chunk_size = 500000  # 每块 50 万点
        total_points = nu * nv * nw
        mask_flat = np.zeros(total_points, dtype=bool)

        # 生成分数坐标索引
        ii = np.arange(nu, dtype=np.float64) / nu
        jj = np.arange(nv, dtype=np.float64) / nv
        kk = np.arange(nw, dtype=np.float64) / nw

        for i_start in range(0, nu, max(1, chunk_size // (nv * nw))):
            i_end = min(i_start + max(1, chunk_size // (nv * nw)), nu)
            fi = ii[i_start:i_end]
            fi_grid, fj_grid, fk_grid = np.meshgrid(fi, jj, kk, indexing='ij')
            frac = np.stack([fi_grid.ravel(), fj_grid.ravel(), fk_grid.ravel()], axis=1)
            cart = frac @ orth_mat.T

            # 查询最近原子距离
            dists, _ = atom_tree.query(cart, k=1)

            # 写入 mask
            start_idx = i_start * nv * nw
            end_idx = start_idx + len(dists)
            mask_flat[start_idx:end_idx] = (dists <= mask_r)

        mask = mask_flat.reshape(nu, nv, nw)
        n_masked = mask.sum()
        logger.info(f"  分子掩膜: radius={mask_r:.1f}Å, "
                    f"掩膜内体素={n_masked} ({100*n_masked/total_points:.1f}%)")

        return mask

    def compute_all_cc(self, exp_grid, sim_grid, molecular_mask):
        """
        计算四项 CC 指标 (Afonine et al., Acta Cryst. D74, 2018)

        所有 CC 均为 Pearson 相关系数（先减均值再计算）。

        - CC_box:    整个 box 内所有体素
        - CC_mask:   分子掩膜内的体素（由原子模型定义）
        - CC_volume: 模拟密度最高 N 个体素（N = 掩膜内体素数）
        - CC_peaks:  模拟 + 实验各自最高 N 个体素的并集

        参数:
            exp_grid: gemmi.FloatGrid 实验密度
            sim_grid: gemmi.FloatGrid 模拟密度
            molecular_mask: np.ndarray(bool) 分子掩膜
        返回:
            dict: cc_box, cc_mask, cc_volume, cc_peaks
        """
        exp_data = np.array(exp_grid, copy=True).astype(np.float64)
        sim_data = np.array(sim_grid, copy=True).astype(np.float64)

        results = {}
        N = int(molecular_mask.sum())  # 掩膜内体素数

        # CC_box: 整个 box
        results["cc_box"] = _pearson_cc(exp_data.ravel(), sim_data.ravel())

        # CC_mask: 分子掩膜内
        if N > 0:
            results["cc_mask"] = _pearson_cc(exp_data[molecular_mask],
                                              sim_data[molecular_mask])
        else:
            results["cc_mask"] = 0.0

        # CC_volume: 模拟密度最高 N 个体素
        if N > 0:
            sim_flat = sim_data.ravel()
            # argpartition 比 argsort 更快（O(n) vs O(n log n)）
            top_n_indices = np.argpartition(sim_flat, -N)[-N:]
            results["cc_volume"] = _pearson_cc(exp_data.ravel()[top_n_indices],
                                                sim_flat[top_n_indices])
        else:
            results["cc_volume"] = 0.0

        # CC_peaks: 模拟 + 实验各自最高 N 个体素的并集
        if N > 0:
            sim_flat = sim_data.ravel()
            exp_flat = exp_data.ravel()
            top_sim = set(np.argpartition(sim_flat, -N)[-N:].tolist())
            top_exp = set(np.argpartition(exp_flat, -N)[-N:].tolist())
            peaks_indices = np.array(sorted(top_sim | top_exp))
            results["cc_peaks"] = _pearson_cc(exp_flat[peaks_indices],
                                               sim_flat[peaks_indices])
        else:
            results["cc_peaks"] = 0.0

        return results

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
        n_orig_atoms = sum(
            1 for ch in structure[0] for res in ch
            if res.name != "HOH" for _ in res
        )

        # 生物组装体展开
        if self.bio_assembly_enabled:
            structure = expand_to_assembly(
                structure,
                assembly_id=self.assembly_id,
                max_chains=self.max_chains
            )

        atoms = self.get_atom_coords(structure)
        logger.info(f"  原子数: {n_orig_atoms} -> {len(atoms)} (展开后)")

        # 在原子位置采样密度
        quality_score, diag_info = self.compute_cc_at_atoms(grid, structure)
        logger.info(f"  质量分数 (正密度比例): {quality_score:.4f}")

        # 获取分辨率
        resolution = self._get_resolution(entry_dir)

        # 生成模拟密度图
        sim_grid = self.generate_simulated_map(grid, structure, resolution)

        # 构建分子掩膜（Afonine 2018 标准）
        molecular_mask = self._build_molecular_mask(grid, structure, resolution)

        # 计算全部 CC 指标（Afonine 2018）
        cc_results = self.compute_all_cc(grid, sim_grid, molecular_mask)
        logger.info(f"  CC_box={cc_results['cc_box']:.4f}, "
                    f"CC_mask={cc_results['cc_mask']:.4f}, "
                    f"CC_volume={cc_results['cc_volume']:.4f}, "
                    f"CC_peaks={cc_results['cc_peaks']:.4f}")

        # 综合判断：CC_mask 或 CC_volume 有一项达标即通过
        cc_mask = cc_results['cc_mask']
        cc_volume = cc_results['cc_volume']
        passed = (quality_score >= 0.5 or
                  cc_mask >= self.cc_threshold or
                  cc_volume >= self.cc_threshold)
        quality = "pass" if passed else "fail"
        logger.info(f"  质量判定: {quality}")

        metrics = {
            "cc_box": cc_results['cc_box'],
            "cc_mask": cc_mask,
            "cc_volume": cc_results['cc_volume'],
            "cc_peaks": cc_results['cc_peaks'],
            "quality_score": float(quality_score),
            "resolution": resolution,
            "n_atoms_orig": n_orig_atoms,
            "n_atoms_expanded": len(atoms),
            "n_chains": len(set(a['chain'] for a in atoms)),
            "bio_assembly_expanded": self.bio_assembly_enabled,
            "density_calculator": self.density_calculator,
            "dc_blur": self.dc_blur,
            "grid_shape": [grid.nu, grid.nv, grid.nw],
            "spacing": [float(s) for s in grid.spacing],
            "mask_voxels": int(molecular_mask.sum()),
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

        # 计算 Q-score（原子级密度拟合置信度，Pintilie 2020）
        try:
            qscore_cfg = self.config.get('qscore', {})
            qscores, _ = compute_qscore_per_atom(
                grid, structure,
                sigma=qscore_cfg.get('sigma', 0.6),
                n_directions=qscore_cfg.get('n_directions', 64),
                max_radius=qscore_cfg.get('max_radius', 2.0),
                step_size=qscore_cfg.get('step_size', 0.1),
                n_points_per_shell=qscore_cfg.get('n_points_per_shell', 8),
            )
            qscore_summary = compute_qscore_summary(qscores)
            metrics.update(qscore_summary)
            logger.info(f"  Q-score: mean={qscore_summary['q_score_mean']:.3f}, "
                        f"median={qscore_summary['q_score_median']:.3f}")

            # 持久化逐原子 Q-score，供 correspondence 步骤生成 label_qscore.mrc
            # key 格式: "{chain}_{resseq}_{atomname}" → float
            qscore_dict = {}
            for (chain, resseq, atomname), val in qscores.items():
                key = f"{chain}_{resseq}_{atomname}"
                qscore_dict[key] = round(val, 6)
            qscore_path = os.path.join(entry_dir, "qscores.json")
            with open(qscore_path, 'w') as f:
                json.dump(qscore_dict, f)
            logger.info(f"  Q-score 逐原子数据已保存: {qscore_path} ({len(qscore_dict)} 原子)")

        except Exception as e:
            logger.warning(f"  Q-score 计算失败: {e}")

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
