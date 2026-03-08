"""
模块6: Enhancement Labeling - 增强标注

从原子模型生成模拟密度图（低噪声），作为去噪训练标签
使用 gemmi 坐标系
"""
import os
import logging
import numpy as np
import gemmi
from scipy.ndimage import gaussian_filter

logger = logging.getLogger(__name__)


class EnhancementLabeler:
    """从原子模型生成模拟密度图"""

    def __init__(self, config):
        self.config = config
        self.b_factor = config['enhancement']['sim_b_factor']
        self.resolution = config['enhancement']['sim_resolution']

    def generate_mol_map(self, entry_dir):
        """生成模拟密度图 mol_map.mrc"""
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
        model = st[0]

        # 元素权重（归一化散射因子的近似）
        element_weight = {
            'C': 1.0, 'N': 1.17, 'O': 1.33, 'S': 2.67, 'P': 2.5,
            'H': 0.17, 'FE': 4.33, 'ZN': 5.0, 'MG': 2.0,
        }

        mol_data = np.zeros((nu, nv, nw), dtype=np.float32)
        n_atoms = 0

        # 三线性插值放置原子（保留亚体素精度）
        for chain in model:
            for residue in chain:
                if residue.name == "HOH":
                    continue
                for atom in residue:
                    frac = cell.fractionalize(atom.pos)
                    fi, fj, fk = frac.x * nu, frac.y * nv, frac.z * nw
                    i0 = int(np.floor(fi)) % nu
                    j0 = int(np.floor(fj)) % nv
                    k0 = int(np.floor(fk)) % nw
                    i1, j1, k1 = (i0+1)%nu, (j0+1)%nv, (k0+1)%nw
                    di = fi - np.floor(fi)
                    dj = fj - np.floor(fj)
                    dk = fk - np.floor(fk)
                    w = element_weight.get(atom.element.name, 1.0)
                    mol_data[i0,j0,k0] += w*(1-di)*(1-dj)*(1-dk)
                    mol_data[i1,j0,k0] += w*di*(1-dj)*(1-dk)
                    mol_data[i0,j1,k0] += w*(1-di)*dj*(1-dk)
                    mol_data[i0,j0,k1] += w*(1-di)*(1-dj)*dk
                    mol_data[i1,j1,k0] += w*di*dj*(1-dk)
                    mol_data[i1,j0,k1] += w*di*(1-dj)*dk
                    mol_data[i0,j1,k1] += w*(1-di)*dj*dk
                    mol_data[i1,j1,k1] += w*di*dj*dk
                    n_atoms += 1

        logger.info(f"  放置了 {n_atoms} 个原子（三线性插值）")

        # 高斯模糊
        sigma_A = np.sqrt(self.b_factor / (8 * np.pi**2))
        sigma_voxel = sigma_A / ref_grid.spacing[0]
        mol_data = gaussian_filter(mol_data, sigma=sigma_voxel)

        if mol_data.max() > 0:
            mol_data /= mol_data.max()

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
                logger.error(f"  增强标注失败: {e}")

        logger.info(f"增强标注完成: {len(results)} 个条目")
        return results
