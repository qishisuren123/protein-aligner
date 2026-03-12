"""
模块: Molmap - ChimeraX 规范模拟密度图生成

按照 ChimeraX molmap 命令的默认参数生成模拟密度图：
- sigma = resolution / (pi * sqrt(2)) ≈ 0.225 * resolution
- B_iso = 8 * pi^2 * sigma^2 = 4 * resolution^2
- 振幅 = 原子电子散射因子（gemmi DensityCalculatorE）
- cutoff = 5 * sigma

参考: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/molmap.html
"""
import math
import logging
import numpy as np
import gemmi
from scipy.ndimage import zoom

logger = logging.getLogger(__name__)


def generate_molmap(grid, structure, resolution):
    """
    按 ChimeraX molmap 规范生成模拟密度图

    参数:
        grid: gemmi.FloatGrid 参考 grid（提供尺寸和坐标系）
        structure: gemmi.Structure 原子模型（已展开组装体）
        resolution: float 分辨率 (Å)

    返回:
        mol_data: numpy.ndarray 归一化后的模拟密度 [0, 1], shape=(nu, nv, nw)
    """
    # ChimeraX molmap 默认: sigma = resolution / (pi * sqrt(2))
    sigma_factor = 1.0 / (math.pi * math.sqrt(2))
    sigma = sigma_factor * resolution
    # B = 8 * pi^2 * sigma^2 = 4 * resolution^2
    B_iso = 8.0 * math.pi ** 2 * sigma ** 2

    logger.info(f"  Molmap 参数: resolution={resolution:.2f}Å, "
                f"sigma={sigma:.3f}Å, B_iso={B_iso:.1f}")

    cell = grid.unit_cell
    nu, nv, nw = grid.nu, grid.nv, grid.nw

    # 克隆结构，清零原子 B-factor（匹配 molmap 均匀高斯行为）
    st_copy = structure.clone()
    model = st_copy[0]
    n_atoms = 0
    for chain in model:
        for res in chain:
            if res.name == "HOH":
                continue
            for atom in res:
                atom.b_iso = 0.0
                n_atoms += 1

    logger.info(f"  原子数: {n_atoms}")

    # 使用 DensityCalculatorE（电子散射因子）
    dc = gemmi.DensityCalculatorE()
    dc.d_min = resolution
    dc.blur = B_iso
    dc.grid.set_unit_cell(cell)
    dc.grid.set_size(nu, nv, nw)
    dc.grid.spacegroup = grid.spacegroup

    dc.put_model_density_on_grid(model)
    mol_data = np.array(dc.grid, copy=True)

    # 清除 NaN
    nan_count = np.isnan(mol_data).sum()
    if nan_count > 0:
        logger.info(f"  清除 {nan_count} 个 NaN 值")
        mol_data = np.nan_to_num(mol_data, nan=0.0)

    # DC 可能自动调整 grid 尺寸，需要重采样
    target_shape = (nu, nv, nw)
    if mol_data.shape != target_shape:
        logger.info(f"  DC grid {mol_data.shape} -> 重采样到 {target_shape}")
        zoom_factors = tuple(t / s for t, s in zip(target_shape, mol_data.shape))
        mol_data = zoom(mol_data, zoom_factors, order=3)

    # 归一化到 [0, 1]
    if mol_data.max() > 0:
        mol_data = mol_data / mol_data.max()

    return mol_data.astype(np.float32)
