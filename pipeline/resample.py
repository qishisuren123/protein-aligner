"""
模块2: Resample - 密度图重采样

功能：
- 将不同体素大小的密度图重采样到统一的目标体素大小
- 使用 scipy.ndimage.zoom 做快速重采样
- 使用 gemmi 正确读取/保存 CCP4/MRC 格式的坐标信息
"""
import os
import logging
import numpy as np
import gemmi
import mrcfile
from scipy.ndimage import zoom

logger = logging.getLogger(__name__)


class Resampler:
    """将密度图重采样到统一体素大小"""

    def __init__(self, config):
        self.config = config
        self.target_voxel_sizes = config['resample']['target_voxel_sizes']
        self.interp_order = config['resample']['interpolation_order']

    def resample_map(self, input_path, output_path, target_voxel_size):
        """
        将 MRC 密度图重采样到目标体素大小

        参数:
            input_path: 输入 MRC/MAP 文件路径
            output_path: 输出 MRC 文件路径
            target_voxel_size: 目标体素大小 (Å)
        """
        # 使用 gemmi 读取（正确处理坐标系）
        ccp4 = gemmi.read_ccp4_map(input_path)
        ccp4.setup(float('nan'))
        grid = ccp4.grid
        cell = grid.unit_cell
        data = np.array(grid, copy=True)

        orig_spacing = grid.spacing
        logger.info(f"  原始体素大小: ({orig_spacing[0]:.3f}, {orig_spacing[1]:.3f}, {orig_spacing[2]:.3f})")
        logger.info(f"  原始 Grid 维度: ({grid.nu}, {grid.nv}, {grid.nw})")
        logger.info(f"  Unit cell: ({cell.a:.1f}, {cell.b:.1f}, {cell.c:.1f})")
        logger.info(f"  目标体素大小: {target_voxel_size}Å")

        # 计算缩放因子
        zoom_factors = np.array([
            orig_spacing[0] / target_voxel_size,
            orig_spacing[1] / target_voxel_size,
            orig_spacing[2] / target_voxel_size,
        ])

        if np.allclose(zoom_factors, 1.0, atol=0.01):
            logger.info("  体素大小已接近目标值，直接复制")
            resampled = data.copy()
        else:
            logger.info(f"  缩放因子: {zoom_factors}")
            resampled = zoom(data, zoom_factors, order=self.interp_order)

        logger.info(f"  重采样后形状: {resampled.shape}")

        # 使用 gemmi 保存，保持正确坐标系
        new_grid = gemmi.FloatGrid(resampled.shape[0], resampled.shape[1], resampled.shape[2])
        new_grid.set_unit_cell(cell)  # unit cell 不变（物理空间一样大）
        new_grid.spacegroup = grid.spacegroup

        # 将 numpy 数据复制到 gemmi grid
        arr = np.array(new_grid, copy=False)
        arr[:] = resampled

        new_ccp4 = gemmi.Ccp4Map()
        new_ccp4.grid = new_grid
        new_ccp4.update_ccp4_header()
        new_ccp4.write_ccp4_map(output_path)

        logger.info(f"  重采样完成: {output_path}")
        return output_path

    def run(self, entry_dirs):
        """对所有条目的密度图进行重采样"""
        results = []
        for entry_dir in entry_dirs:
            raw_map = os.path.join(entry_dir, "raw_map.map")
            if not os.path.exists(raw_map):
                logger.warning(f"跳过 {entry_dir}: 未找到 raw_map.map")
                continue

            logger.info(f"处理: {entry_dir}")
            for vs in self.target_voxel_sizes:
                output_path = os.path.join(entry_dir, f"map_std_{vs:.1f}A.mrc")
                try:
                    self.resample_map(raw_map, output_path, vs)
                    results.append({
                        "entry_dir": entry_dir,
                        "voxel_size": vs,
                        "output": output_path
                    })
                except Exception as e:
                    logger.error(f"  重采样失败 ({vs}Å): {e}", exc_info=True)

        return results
