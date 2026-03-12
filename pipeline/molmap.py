"""
模块: Molmap - 按 ChimeraX molmap 规范生成模拟密度图

ChimeraX molmap 算法核心：
- 每个原子用一个高斯 blob 表示
- 振幅 = 原子序数 Z（非电子散射因子）
- sigma = resolution / (pi * sqrt(2)) ≈ 0.225 * resolution
- 内部网格步长 = resolution / 3（保证高斯被充分采样）
- 输出重采样到参考网格尺寸

实现方式：三线性插值放置原子到细网格 → 3D 高斯模糊 → 重采样
等价于 ChimeraX 的逐原子高斯求和，但更高效

参考: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/molmap.html
"""
import math
import logging
import numpy as np
import gemmi
from scipy.ndimage import gaussian_filter, zoom

logger = logging.getLogger(__name__)


def generate_molmap(grid, structure, resolution):
    """
    按 ChimeraX molmap 规范生成模拟密度图

    算法：
    1. 创建内部细网格（步长 = resolution/3）
    2. 将每个原子的原子序数 Z 通过三线性插值放置到细网格上
    3. 用 sigma = resolution / (pi * sqrt(2)) 的高斯核卷积
    4. 重采样到参考 grid 尺寸
    5. 归一化到 [0, 1]

    参数:
        grid: gemmi.FloatGrid 参考 grid（提供尺寸和坐标系）
        structure: gemmi.Structure 原子模型（已展开组装体）
        resolution: float 分辨率 (Å)

    返回:
        mol_data: numpy.ndarray 归一化后的模拟密度 [0, 1], shape=(nu, nv, nw)
    """
    # ChimeraX molmap 默认参数
    sigma = resolution / (math.pi * math.sqrt(2))

    cell = grid.unit_cell
    nu, nv, nw = grid.nu, grid.nv, grid.nw

    # 参考网格步长 (Å)
    ref_spacing = np.array([cell.a / nu, cell.b / nv, cell.c / nw])

    # ChimeraX 内部细网格：步长 = resolution / 3
    fine_step = resolution / 3.0
    # 确保细网格步长不超过参考网格步长
    fine_step = min(fine_step, ref_spacing.min())
    fine_nu = max(int(round(cell.a / fine_step)), nu)
    fine_nv = max(int(round(cell.b / fine_step)), nv)
    fine_nw = max(int(round(cell.c / fine_step)), nw)

    fine_spacing = np.array([cell.a / fine_nu, cell.b / fine_nv, cell.c / fine_nw])
    sigma_grid = sigma / fine_spacing

    logger.info(f"  Molmap 参数: resolution={resolution:.2f}Å, sigma={sigma:.3f}Å")
    logger.info(f"  参考 grid: {nu}x{nv}x{nw}, 步长={ref_spacing}")
    logger.info(f"  细网格: {fine_nu}x{fine_nv}x{fine_nw}, 步长={fine_spacing}")
    logger.info(f"  sigma_grid = ({sigma_grid[0]:.2f}, "
                f"{sigma_grid[1]:.2f}, {sigma_grid[2]:.2f})")

    # 收集原子坐标和原子序数
    positions = []
    atomic_nums = []
    model = structure[0]
    for chain in model:
        for res in chain:
            if res.name == "HOH":
                continue
            for atom in res:
                z = atom.element.atomic_number
                if z == 0:
                    continue
                positions.append([atom.pos.x, atom.pos.y, atom.pos.z])
                atomic_nums.append(z)

    n_atoms = len(positions)
    logger.info(f"  原子数: {n_atoms}")

    if n_atoms == 0:
        return np.zeros((nu, nv, nw), dtype=np.float32)

    positions = np.array(positions, dtype=np.float64)
    atomic_nums = np.array(atomic_nums, dtype=np.float64)

    # 将笛卡尔坐标转换为分数坐标
    frac_mat = np.array(cell.frac.mat.tolist())
    frac_coords = positions @ frac_mat.T

    # 分数坐标 -> 细网格索引（浮点）
    grid_coords = frac_coords * np.array([fine_nu, fine_nv, fine_nw])

    # 三线性插值放置原子到密度网格
    density = np.zeros((fine_nu, fine_nv, fine_nw), dtype=np.float64)

    gi0 = np.floor(grid_coords).astype(int)
    gf = grid_coords - gi0

    for di in range(2):
        for dj in range(2):
            for dk in range(2):
                wi = gf[:, 0] if di == 1 else (1.0 - gf[:, 0])
                wj = gf[:, 1] if dj == 1 else (1.0 - gf[:, 1])
                wk = gf[:, 2] if dk == 1 else (1.0 - gf[:, 2])
                w = wi * wj * wk

                ii = (gi0[:, 0] + di) % fine_nu
                jj = (gi0[:, 1] + dj) % fine_nv
                kk = (gi0[:, 2] + dk) % fine_nw

                np.add.at(density, (ii, jj, kk), atomic_nums * w)

    # 高斯模糊
    density = gaussian_filter(density, sigma=sigma_grid)

    # 重采样到参考 grid 尺寸
    if (fine_nu, fine_nv, fine_nw) != (nu, nv, nw):
        zoom_factors = (nu / fine_nu, nv / fine_nv, nw / fine_nw)
        logger.info(f"  重采样: {fine_nu}x{fine_nv}x{fine_nw} -> {nu}x{nv}x{nw}")
        density = zoom(density, zoom_factors, order=3)

    # 确保尺寸精确匹配（zoom 可能有 ±1 偏差）
    if density.shape != (nu, nv, nw):
        result = np.zeros((nu, nv, nw), dtype=np.float64)
        sn = min(density.shape[0], nu)
        sv = min(density.shape[1], nv)
        sw = min(density.shape[2], nw)
        result[:sn, :sv, :sw] = density[:sn, :sv, :sw]
        density = result

    # 清除 NaN
    nan_count = np.isnan(density).sum()
    if nan_count > 0:
        logger.info(f"  清除 {nan_count} 个 NaN 值")
        density = np.nan_to_num(density, nan=0.0)

    # 归一化到 [0, 1]
    if density.max() > 0:
        density = density / density.max()

    return density.astype(np.float32)
