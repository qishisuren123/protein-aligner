"""
模块6: Enhancement Labeling - 增强标注

功能：
- 从原子模型生成模拟密度图（低噪声），作为去噪训练标签
- 按 ChimeraX molmap 规范生成，与 alignment_qc sim_map 参数一致
- 支持生物学组装体展开
- 参考: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/molmap.html
"""
import os
import json
import logging
import numpy as np
import gemmi

from pipeline.bio_assembly import expand_to_assembly

logger = logging.getLogger(__name__)


class EnhancementLabeler:
    """从原子模型生成模拟密度图"""

    def __init__(self, config):
        self.config = config
        enh_cfg = config['enhancement']
        self.resolution = enh_cfg.get('sim_resolution', 2.0)
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
        按 ChimeraX molmap 规范生成模拟密度图 mol_map.mrc

        使用与 alignment_qc sim_map 相同的 molmap 参数，保证输出一致。
        参考: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/molmap.html
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

        # 获取分辨率
        resolution = self._get_resolution(entry_dir)

        # 按 molmap 规范生成模拟密度
        from pipeline.molmap import generate_molmap
        mol_data = generate_molmap(ref_grid, st, resolution)

        # 保存
        output_path = os.path.join(entry_dir, "mol_map.mrc")
        new_grid = gemmi.FloatGrid(nu, nv, nw)
        new_grid.set_unit_cell(cell)
        new_grid.spacegroup = ref_grid.spacegroup
        arr = np.array(new_grid, copy=False)
        arr[:] = mol_data

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
