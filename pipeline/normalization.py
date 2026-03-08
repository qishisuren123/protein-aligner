"""
模块3: Normalization - 密度图归一化

功能：
- 将密度图值归一化到 [0, 1] 范围
- 使用第 95 百分位数进行归一化
- 使用 gemmi 保持坐标系正确
"""
import os
import logging
import numpy as np
import gemmi

logger = logging.getLogger(__name__)


class Normalizer:
    """密度图值归一化"""

    def __init__(self, config):
        self.config = config
        self.percentile = config['normalization']['percentile']
        self.clip_min = config['normalization']['clip_min']
        self.clip_max = config['normalization']['clip_max']

    def normalize_map(self, input_path, output_path):
        """归一化密度图: new_value = clip(raw_value / 95th_percentile, 0, 1)"""
        ccp4 = gemmi.read_ccp4_map(input_path)
        ccp4.setup(float('nan'))
        grid = ccp4.grid
        data = np.array(grid, copy=True)

        logger.info(f"  原始值范围: [{data.min():.4f}, {data.max():.4f}]")
        logger.info(f"  原始均值: {data.mean():.4f}, 标准差: {data.std():.4f}")

        # 计算正值区域的第 N 百分位数
        positive_values = data[data > 0]
        if len(positive_values) == 0:
            logger.warning("  没有正值体素，使用全局百分位数")
            p_value = np.percentile(data, self.percentile)
        else:
            p_value = np.percentile(positive_values, self.percentile)

        logger.info(f"  第{self.percentile}百分位数: {p_value:.4f}")

        if p_value <= 0:
            p_value = data.max()

        # 归一化
        normalized = np.clip(data / p_value, self.clip_min, self.clip_max)
        logger.info(f"  归一化后值范围: [{normalized.min():.4f}, {normalized.max():.4f}]")

        # 使用 gemmi 保存（保持坐标系）
        new_grid = gemmi.FloatGrid(normalized.shape[0], normalized.shape[1], normalized.shape[2])
        new_grid.set_unit_cell(grid.unit_cell)
        new_grid.spacegroup = grid.spacegroup
        arr = np.array(new_grid, copy=False)
        arr[:] = normalized.astype(np.float32)

        new_ccp4 = gemmi.Ccp4Map()
        new_ccp4.grid = new_grid
        new_ccp4.update_ccp4_header()
        new_ccp4.write_ccp4_map(output_path)

        return output_path

    def run(self, entry_dirs, voxel_size=1.0):
        """对所有条目的重采样后密度图进行归一化"""
        results = []
        for entry_dir in entry_dirs:
            input_map = os.path.join(entry_dir, f"map_std_{voxel_size:.1f}A.mrc")
            if not os.path.exists(input_map):
                logger.warning(f"跳过 {entry_dir}: 未找到 {input_map}")
                continue

            output_path = os.path.join(entry_dir, "map_normalized.mrc")
            logger.info(f"归一化: {entry_dir}")

            try:
                self.normalize_map(input_map, output_path)
                results.append({"entry_dir": entry_dir, "output": output_path})
            except Exception as e:
                logger.error(f"  归一化失败: {e}")

        return results
