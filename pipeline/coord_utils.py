"""
坐标工具 - 处理 MRC 密度图与原子模型之间的坐标转换

使用 gemmi 的 CCP4 map 读取确保正确处理坐标系
"""
import numpy as np
import gemmi
import mrcfile


def load_map_gemmi(map_path):
    """
    使用 gemmi 加载 MRC/CCP4 密度图，正确处理坐标系

    返回:
        grid: gemmi Grid 对象（支持坐标插值）
        data: numpy 数组
        metadata: 元数据字典
    """
    ccp4 = gemmi.read_ccp4_map(map_path)
    ccp4.setup(float('nan'))
    grid = ccp4.grid

    data = np.array(grid, copy=True)

    # 提取关键元数据
    unit_cell = grid.unit_cell
    metadata = {
        'shape': data.shape,
        'spacing': (grid.spacing[0], grid.spacing[1], grid.spacing[2]),
        'unit_cell': (unit_cell.a, unit_cell.b, unit_cell.c),
    }

    return grid, data, metadata


def save_map_with_ref(data, ref_map_path, output_path):
    """
    使用参考 map 的坐标系保存新的密度图数据

    参数:
        data: numpy 数组
        ref_map_path: 参考 MRC 文件（用于复制 header）
        output_path: 输出路径
    """
    # 用 gemmi 读取参考 map 获取坐标信息
    ref = gemmi.read_ccp4_map(ref_map_path)
    ref.setup(float('nan'))
    ref_grid = ref.grid

    with mrcfile.new(output_path, overwrite=True) as mrc:
        mrc.set_data(data.astype(np.float32))
        # 从参考 map 复制坐标系信息
        with mrcfile.open(ref_map_path, mode='r') as ref_mrc:
            mrc.header.origin = ref_mrc.header.origin
            mrc.header.nxstart = ref_mrc.header.nxstart
            mrc.header.nystart = ref_mrc.header.nystart
            mrc.header.nzstart = ref_mrc.header.nzstart
            mrc.header.cella = ref_mrc.header.cella
            mrc.header.mapc = ref_mrc.header.mapc
            mrc.header.mapr = ref_mrc.header.mapr
            mrc.header.maps = ref_mrc.header.maps
        mrc.update_header_stats()


def atomic_to_fractional(coord, grid):
    """
    将原子坐标 (Å) 转换为 grid 的分数坐标

    参数:
        coord: [x, y, z] 实际坐标 (Å)
        grid: gemmi Grid 对象
    返回:
        分数坐标 [u, v, w]
    """
    pos = gemmi.Position(coord[0], coord[1], coord[2])
    frac = grid.unit_cell.fractionalize(pos)
    return np.array([frac.x, frac.y, frac.z])


def atomic_to_grid_index(coord, grid):
    """
    将原子坐标 (Å) 转换为 grid 整数索引

    参数:
        coord: [x, y, z] 实际坐标 (Å)
        grid: gemmi Grid 对象
    返回:
        grid 索引 [i, j, k]（对应 data[i, j, k]）
    """
    frac = atomic_to_fractional(coord, grid)
    # 分数坐标乘以 grid 维度得到索引
    idx = np.array([
        frac[0] * grid.nu,
        frac[1] * grid.nv,
        frac[2] * grid.nw
    ])
    return np.round(idx).astype(int)


def interpolate_at_position(grid, x, y, z):
    """
    在指定坐标处插值密度值

    参数:
        grid: gemmi Grid 对象
        x, y, z: 实际坐标 (Å)
    返回:
        插值后的密度值
    """
    return grid.interpolate_value(gemmi.Position(x, y, z))
