"""
模块5b: Interface Detection - 蛋白-蛋白界面检测

功能：
- 基于 CA-CA 跨链距离检测蛋白-蛋白界面残基
- 对每条链的 CA 原子，批量查询其他链的 KDTree
- 若最近跨链 CA 距离 < 阈值 → interface(1)，否则 non-interface(0)
"""
import logging
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict

logger = logging.getLogger(__name__)


def detect_interface(ca_atoms, threshold=5.0):
    """
    检测蛋白-蛋白界面残基

    算法：
    1. 按链分组 CA 原子坐标（真实 Å 空间）
    2. 每条链建 cKDTree
    3. 对每条链的所有 CA 原子，批量查询其他链的 KDTree
    4. 若最近跨链 CA 距离 < threshold → interface(1)

    参数:
        ca_atoms: list of dict，每个包含 'coord' (Å 坐标) 和 'chain'
        threshold: CA-CA 跨链距离阈值 (Å)，默认 5.0

    返回:
        interface_labels: np.array of int，与 ca_atoms 等长，1=界面，0=非界面
    """
    n_atoms = len(ca_atoms)
    interface_labels = np.zeros(n_atoms, dtype=np.int16)

    if n_atoms == 0:
        return interface_labels

    # 按链分组（记录全局索引）
    chain_groups = defaultdict(list)
    for i, atom in enumerate(ca_atoms):
        chain_groups[atom['chain']].append(i)

    chain_names = sorted(chain_groups.keys())

    # 只有单链，无界面
    if len(chain_names) <= 1:
        logger.info(f"  界面检测: 仅 {len(chain_names)} 条链，无跨链界面")
        return interface_labels

    # 每条链构建 KDTree（真实 Å 坐标）
    chain_trees = {}
    chain_coords = {}
    for cname in chain_names:
        indices = chain_groups[cname]
        coords = np.array([ca_atoms[i]['coord'] for i in indices])
        chain_coords[cname] = coords
        chain_trees[cname] = cKDTree(coords)

    # 对每条链的 CA 原子，查询所有其他链的最近距离
    n_interface = 0
    for cname in chain_names:
        indices = chain_groups[cname]
        coords = chain_coords[cname]
        other_chains = [c for c in chain_names if c != cname]

        if not other_chains:
            continue

        # 对当前链的每个 CA 原子，找跨链最近距离
        min_cross_dist = np.full(len(indices), np.inf)

        for other_name in other_chains:
            other_tree = chain_trees[other_name]
            # 批量查询：一次查完当前链所有 CA 到另一条链的最近距离
            dists, _ = other_tree.query(coords, k=1)
            min_cross_dist = np.minimum(min_cross_dist, dists)

        # 标记界面残基
        is_interface = min_cross_dist < threshold
        for local_idx, global_idx in enumerate(indices):
            if is_interface[local_idx]:
                interface_labels[global_idx] = 1
                n_interface += 1

    logger.info(f"  界面检测: {n_interface}/{n_atoms} CA 原子为界面 "
                f"(阈值={threshold}Å, {len(chain_names)} 条链)")
    return interface_labels
