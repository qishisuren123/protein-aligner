"""
模块3: Normalization - 密度图归一化

功能：
- 支持多种归一化方法：robust_zscore（默认）、zscore、percentile
- robust z-score 基于中位数和 MAD，保留密度动态范围
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
        norm_cfg = config['normalization']
        self.method = norm_cfg.get('method', 'robust_zscore')
        self.percentile = norm_cfg.get('percentile', 95)
        self.clip_min = norm_cfg.get('clip_min', 0.0)
        self.clip_max = norm_cfg.get('clip_max', 1.0)

    def _robust_zscore_normalize(self, data):
        """
        Robust z-score 归一化：基于中位数和 MAD

        公式: normalized = (x - median) / (1.4826 * MAD)
        然后将正值区域映射到 [0, 1]

        优点：对异常值不敏感，保留密度动态范围
        """
        # 只在有信号的区域计算统计量
        positive = data[data > 0]
        if len(positive) == 0:
            logger.warning("  没有正值体素")
            return np.zeros_like(data)

        median = np.median(positive)
        mad = np.median(np.abs(positive - median))

        # MAD 为 0 时回退到标准差
        if mad < 1e-10:
            logger.warning("  MAD 接近 0，回退到 z-score")
            return self._zscore_normalize(data)

        # 1.4826 使 MAD 与正态分布的标准差一致
        scale = 1.4826 * mad
        robust_z = (data - median) / scale

        # 将信号区间映射到 [0, 1]
        # 使正值中位数对应 ~0.5，用正值范围的合理上界做归一化
        p99 = np.percentile(positive, 99)
        upper_z = (p99 - median) / scale
        lower_z = -median / scale  # 0 对应的 z-score

        # 线性映射 [lower_z, upper_z] -> [0, 1]
        normalized = (robust_z - lower_z) / (upper_z - lower_z)
        normalized = np.clip(normalized, self.clip_min, self.clip_max)

        logger.info(f"  Robust z-score: median={median:.4f}, MAD={mad:.4f}, "
                    f"scale={scale:.4f}")
        return normalized.astype(np.float32)

    def _zscore_normalize(self, data):
        """
        标准 z-score 归一化

        公式: normalized = (x - mean) / std，然后映射到 [0, 1]
        """
        positive = data[data > 0]
        if len(positive) == 0:
            return np.zeros_like(data)

        mean = positive.mean()
        std = positive.std()
        if std < 1e-10:
            logger.warning("  标准差接近 0")
            return np.zeros_like(data)

        z = (data - mean) / std
        # 映射到 [0, 1]: 0 -> clip_min, mean+3*std -> 1
        p99 = np.percentile(positive, 99)
        upper_z = (p99 - mean) / std
        lower_z = -mean / std

        normalized = (z - lower_z) / (upper_z - lower_z)
        normalized = np.clip(normalized, self.clip_min, self.clip_max)

        logger.info(f"  Z-score: mean={mean:.4f}, std={std:.4f}")
        return normalized.astype(np.float32)

    def _percentile_normalize(self, data):
        """
        百分位数归一化（原有方法）

        new_value = clip(raw_value / 95th_percentile, 0, 1)
        """
        positive_values = data[data > 0]
        if len(positive_values) == 0:
            logger.warning("  没有正值体素，使用全局百分位数")
            p_value = np.percentile(data, self.percentile)
        else:
            p_value = np.percentile(positive_values, self.percentile)

        logger.info(f"  第{self.percentile}百分位数: {p_value:.4f}")

        if p_value <= 0:
            p_value = data.max()

        normalized = np.clip(data / p_value, self.clip_min, self.clip_max)
        return normalized.astype(np.float32)

    def normalize_map(self, input_path, output_path):
        """归一化密度图"""
        ccp4 = gemmi.read_ccp4_map(input_path)
        ccp4.setup(float('nan'))
        grid = ccp4.grid
        data = np.array(grid, copy=True)

        logger.info(f"  原始值范围: [{data.min():.4f}, {data.max():.4f}]")
        logger.info(f"  原始均值: {data.mean():.4f}, 标准差: {data.std():.4f}")
        logger.info(f"  归一化方法: {self.method}")

        # 根据配置选择归一化方法
        if self.method == 'robust_zscore':
            normalized = self._robust_zscore_normalize(data)
        elif self.method == 'zscore':
            normalized = self._zscore_normalize(data)
        elif self.method == 'percentile':
            normalized = self._percentile_normalize(data)
        else:
            logger.warning(f"  未知方法 '{self.method}'，使用 robust_zscore")
            normalized = self._robust_zscore_normalize(data)

        logger.info(f"  归一化后值范围: [{normalized.min():.4f}, {normalized.max():.4f}]")

        # 使用 gemmi 保存（保持坐标系）
        new_grid = gemmi.FloatGrid(normalized.shape[0], normalized.shape[1], normalized.shape[2])
        new_grid.set_unit_cell(grid.unit_cell)
        new_grid.spacegroup = grid.spacegroup
        arr = np.array(new_grid, copy=False)
        arr[:] = normalized

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
