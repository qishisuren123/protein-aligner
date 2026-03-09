"""
模块6: Enhancement Labeling - 增强标注

功能：
- 从原子模型生成模拟密度图（低噪声），作为去噪训练标签
- 使用 gemmi.DensityCalculatorX 替代手写三线性插值
- 支持生物学组装体展开
- 从 metadata 读取分辨率，用 b_factor 控制模糊程度
"""
import os
import json
import logging
import numpy as np
import gemmi
from scipy.ndimage import zoom

from pipeline.bio_assembly import expand_to_assembly

logger = logging.getLogger(__name__)


class EnhancementLabeler:
    """从原子模型生成模拟密度图"""

    def __init__(self, config):
        self.config = config
        enh_cfg = config['enhancement']
        self.b_factor = enh_cfg.get('sim_b_factor', 100.0)
        self.resolution = enh_cfg.get('sim_resolution', 2.0)
        self.adaptive_blur = enh_cfg.get('adaptive_blur', False)
        # 生物组装体配置
        bio_cfg = config.get('bio_assembly', {})
        self.bio_assembly_enabled = bio_cfg.get('enabled', True)
        self.assembly_id = bio_cfg.get('assembly_id', '1')
        self.max_chains = bio_cfg.get('max_chains', 200)

    def _get_resolution(self, entry_dir):
        """从 metadata.json 读取分辨率"""
        meta_path = os.path.join(entry_dir, "metadata.json")
        if os.path.exists(meta_path):
            with open(meta_path, 'r') as f:
                meta = json.load(f)
            res = meta.get("resolution")
            if res is not None:
                return float(res)
        return self.resolution

    def generate_mol_map(self, entry_dir):
        """
        使用 DensityCalculatorX 生成模拟密度图 mol_map.mrc

        使用物理散射因子，通过 b_factor 控制平滑程度
        """
        map_path = os.path.join(entry_dir, "map_normalized.mrc")
        model_path = os.path.join(entry_dir, "model.cif")

        if not os.path.exists(map_path) or not os.path.exists(model_path):
            logger.warning(f"跳过 {entry_dir}: 缺少文件")
            return None

        # 读取参考 grid
        ccp4 = gemmi.read_ccp4_map(map_path)
        ccp4.setup(float('nan'))
        ref_grid = ccp4.grid
        cell = ref_grid.unit_cell
        nu, nv, nw = ref_grid.nu, ref_grid.nv, ref_grid.nw

        # 加载结构
        st = gemmi.read_structure(model_path)
        st.setup_entities()

        # 生物组装体展开
        if self.bio_assembly_enabled:
            st = expand_to_assembly(
                st,
                assembly_id=self.assembly_id,
                max_chains=self.max_chains
            )

        model = st[0]
        n_atoms = sum(
            1 for chain in model for res in chain
            if res.name != "HOH" for _ in res
        )

        # 获取分辨率
        resolution = self._get_resolution(entry_dir)

        # 使用 DensityCalculatorX 生成密度
        # blur = b_factor 来控制模拟密度的平滑度
        # 自适应模式：按分辨率缩放，高分辨率保留更多细节
        if self.adaptive_blur:
            blur = self.b_factor * (resolution / 2.0)
            logger.info(f"  自适应 blur: base={self.b_factor}, "
                        f"resolution={resolution:.2f}Å, blur={blur:.1f}")
        else:
            blur = self.b_factor

        dc = gemmi.DensityCalculatorX()
        dc.d_min = resolution
        dc.blur = blur
        dc.grid.set_unit_cell(cell)
        dc.grid.set_size(nu, nv, nw)
        dc.grid.spacegroup = ref_grid.spacegroup

        dc.put_model_density_on_grid(model)
        mol_data = np.array(dc.grid, copy=True)

        # DensityCalculatorX 可能自动调整 grid 尺寸，需要重采样
        target_shape = (nu, nv, nw)
        if mol_data.shape != target_shape:
            logger.info(f"  DC grid {mol_data.shape} -> 重采样到 {target_shape}")
            zoom_factors = tuple(t / s for t, s in zip(target_shape, mol_data.shape))
            mol_data = zoom(mol_data, zoom_factors, order=3)

        logger.info(f"  DensityCalculatorX: {n_atoms} 原子, "
                    f"resolution={resolution:.2f}Å, blur={blur:.1f}")

        # 归一化
        if mol_data.max() > 0:
            mol_data /= mol_data.max()

        # 保存
        output_path = os.path.join(entry_dir, "mol_map.mrc")
        new_grid = gemmi.FloatGrid(nu, nv, nw)
        new_grid.set_unit_cell(cell)
        new_grid.spacegroup = ref_grid.spacegroup
        arr = np.array(new_grid, copy=False)
        arr[:] = mol_data.astype(np.float32)

        new_ccp4 = gemmi.Ccp4Map()
        new_ccp4.grid = new_grid
        new_ccp4.update_ccp4_header()
        new_ccp4.write_ccp4_map(output_path)

        logger.info(f"  模拟密度图已保存: {output_path}")
        return output_path

    def run(self, entry_dirs):
        """对所有条目生成模拟密度图"""
        results = []
        for entry_dir in entry_dirs:
            logger.info(f"增强标注: {entry_dir}")
            try:
                output = self.generate_mol_map(entry_dir)
                if output:
                    results.append({"entry_dir": entry_dir, "output": output})
            except Exception as e:
                logger.error(f"  增强标注失败: {e}", exc_info=True)

        logger.info(f"增强标注完成: {len(results)} 个条目")
        return results
