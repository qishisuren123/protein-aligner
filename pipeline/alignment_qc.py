"""
模块4: Alignment & Quality Control - 对齐与质量控制

功能：
- 验证原子模型与密度图的坐标系对齐
- 使用 gemmi.DensityCalculatorX 生成物理级模拟密度图（原子散射因子）
- 支持生物学组装体展开（对称蛋白完整覆盖密度）
- 计算 CC_mask、CC_volume、CC_overall 三项指标
- 根据质量阈值过滤低质量数据
"""
import os
import json
import logging
import numpy as np
import gemmi
from scipy.stats import pearsonr
from scipy.ndimage import zoom

from pipeline.bio_assembly import expand_to_assembly
from pipeline.qscore import compute_qscore_per_atom, compute_qscore_summary

logger = logging.getLogger(__name__)


class AlignmentQC:
    """对齐验证与质量控制"""

    def __init__(self, config):
        self.config = config
        qc_cfg = config['alignment_qc']
        self.mask_threshold = qc_cfg.get('mask_threshold', 0.05)
        self.cc_threshold = qc_cfg.get('cc_threshold', 0.3)
        self.density_calculator = qc_cfg.get('density_calculator', 'gemmi_dcx')
        self.default_resolution = qc_cfg.get('default_resolution', 2.0)
        self.dc_blur = qc_cfg.get('dc_blur', 0)
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
        使用 gemmi.DensityCalculatorX 生成物理级模拟密度图

        DensityCalculatorX 使用 5 高斯近似的原子散射因子，
        物理上更准确，CC 显著高于手写三线性+高斯模糊
        """
        if resolution is None:
            resolution = self.default_resolution

        cell = grid.unit_cell
        model = structure[0]

        n_atoms = sum(
            1 for chain in model for res in chain
            if res.name != "HOH" for _ in res
        )
        logger.info(f"  使用 DensityCalculatorX (blur={self.dc_blur}) 生成模拟密度图")
        logger.info(f"  原子数: {n_atoms}, 分辨率: {resolution:.2f}Å")

        # 创建 DensityCalculatorX 并设置参数
        dc = gemmi.DensityCalculatorX()
        dc.d_min = resolution
        dc.blur = self.dc_blur
        dc.grid.set_unit_cell(cell)
        dc.grid.set_size(grid.nu, grid.nv, grid.nw)
        dc.grid.spacegroup = grid.spacegroup

        # 放置原子（使用物理散射因子）
        dc.put_model_density_on_grid(model)

        # 提取数据（DensityCalculatorX 可能自动调整 grid 尺寸）
        sim_data = np.array(dc.grid, copy=True)

        # 清除 NaN（DCX 偶尔产生极少量 NaN，会在 zoom 插值时传播）
        nan_count = np.isnan(sim_data).sum()
        if nan_count > 0:
            logger.info(f"  清除 {nan_count} 个 NaN 值")
            sim_data = np.nan_to_num(sim_data, nan=0.0)

        # 如果 DC 自动调整了 grid 尺寸，需要重采样到目标尺寸
        target_shape = (grid.nu, grid.nv, grid.nw)
        if sim_data.shape != target_shape:
            logger.info(f"  DC grid {sim_data.shape} -> 重采样到 {target_shape}")
            zoom_factors = tuple(t / s for t, s in zip(target_shape, sim_data.shape))
            sim_data = zoom(sim_data, zoom_factors, order=3)

        # 归一化到 [0, 1]
        if sim_data.max() > 0:
            sim_data /= sim_data.max()

        # 包装为 gemmi grid
        sim_grid = gemmi.FloatGrid(grid.nu, grid.nv, grid.nw)
        sim_grid.set_unit_cell(cell)
        sim_grid.spacegroup = grid.spacegroup
        arr = np.array(sim_grid, copy=False)
        arr[:] = sim_data.astype(np.float32)

        return sim_grid

    def compute_cc_mask(self, exp_grid, sim_grid, mask_threshold=None):
        """
        计算 CC_mask（模拟密度 mask 区域内的 Pearson 相关系数）
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
        return float(cc)

    def compute_all_cc(self, exp_grid, sim_grid, mask_threshold=None):
        """
        计算三项 CC 指标：CC_mask, CC_volume, CC_overall

        - CC_mask: 模拟密度 mask 区域的 Pearson CC
        - CC_volume: 实验密度 mask 区域的 Pearson CC
        - CC_overall: 全体素 Pearson CC
        """
        if mask_threshold is None:
            mask_threshold = self.mask_threshold

        exp_data = np.array(exp_grid, copy=True)
        sim_data = np.array(sim_grid, copy=True)

        results = {}

        # CC_mask: 模拟密度 mask 内
        sim_mask = sim_data > mask_threshold
        if sim_mask.sum() > 0:
            e, s = exp_data[sim_mask], sim_data[sim_mask]
            if e.std() > 0 and s.std() > 0:
                results["cc_mask"] = float(pearsonr(e, s)[0])
            else:
                results["cc_mask"] = 0.0
        else:
            results["cc_mask"] = 0.0

        # CC_volume: 实验密度 mask 内
        exp_mask = exp_data > mask_threshold
        if exp_mask.sum() > 0:
            e, s = exp_data[exp_mask], sim_data[exp_mask]
            if e.std() > 0 and s.std() > 0:
                results["cc_volume"] = float(pearsonr(e, s)[0])
            else:
                results["cc_volume"] = 0.0
        else:
            results["cc_volume"] = 0.0

        # CC_overall: 全体素（有信号区域）
        combined_mask = sim_mask | exp_mask
        if combined_mask.sum() > 0:
            e, s = exp_data[combined_mask], sim_data[combined_mask]
            if e.std() > 0 and s.std() > 0:
                results["cc_overall"] = float(pearsonr(e, s)[0])
            else:
                results["cc_overall"] = 0.0
        else:
            results["cc_overall"] = 0.0

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

        # 生成模拟密度图（DensityCalculatorX）
        sim_grid = self.generate_simulated_map(grid, structure, resolution)

        # 计算全部 CC 指标
        cc_results = self.compute_all_cc(grid, sim_grid)
        logger.info(f"  CC_mask={cc_results['cc_mask']:.4f}, "
                    f"CC_volume={cc_results['cc_volume']:.4f}, "
                    f"CC_overall={cc_results['cc_overall']:.4f}")

        # 综合判断：CC_mask 或 CC_volume 有一项达标即通过
        cc_mask = cc_results['cc_mask']
        cc_volume = cc_results['cc_volume']
        passed = (quality_score >= 0.5 or
                  cc_mask >= self.cc_threshold or
                  cc_volume >= self.cc_threshold)
        quality = "pass" if passed else "fail"
        logger.info(f"  质量判定: {quality}")

        metrics = {
            "cc_mask": cc_mask,
            "cc_volume": cc_results['cc_volume'],
            "cc_overall": cc_results['cc_overall'],
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

        # 计算 Q-score（原子级密度拟合置信度）
        try:
            qscore_cfg = self.config.get('qscore', {})
            qscores, _ = compute_qscore_per_atom(
                grid, structure,
                sigma=qscore_cfg.get('sigma', 0.6),
                n_directions=qscore_cfg.get('n_directions', 40),
                max_radius=qscore_cfg.get('max_radius', 2.0),
                step_size=qscore_cfg.get('step_size', 0.1),
            )
            qscore_summary = compute_qscore_summary(qscores)
            metrics.update(qscore_summary)
            logger.info(f"  Q-score: mean={qscore_summary['q_score_mean']:.3f}, "
                        f"median={qscore_summary['q_score_median']:.3f}")
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
