"""
Q-score 物理置信度计算 (Pintilie et al., 2020)

Q-score 是 cryo-EM 领域标准的原子级密度拟合质量指标。
对每个原子，沿径向方向采样实验密度，与参考高斯做 Pearson 相关。
返回值范围 [-1, 1]，越接近 1 表示局部密度与原子位置匹配越好。

参考: Pintilie et al., Nature Methods 17, 328–334 (2020)
"""
import logging
import numpy as np
import gemmi

logger = logging.getLogger(__name__)


def _fibonacci_sphere(n_points):
    """
    生成 Fibonacci 球面均匀分布的方向向量

    参数:
        n_points: 采样方向数
    返回:
        directions: (n_points, 3) 的单位方向向量数组
    """
    indices = np.arange(n_points, dtype=np.float64)
    # 黄金角
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    theta = golden_angle * indices
    # z 坐标均匀分布在 [-1, 1]
    z = 1.0 - (2.0 * indices / (n_points - 1)) if n_points > 1 else np.array([0.0])
    radius = np.sqrt(1.0 - z * z)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    return np.column_stack([x, y, z])


def _reference_gaussian(distances, sigma):
    """
    参考高斯衰减 profile

    参数:
        distances: 径向距离数组 (Å)
        sigma: 高斯宽度参数 (Å)
    返回:
        高斯值数组
    """
    return np.exp(-distances ** 2 / (2.0 * sigma ** 2))


def compute_qscore_per_atom(grid, structure, sigma=0.6, n_directions=40,
                            max_radius=2.0, step_size=0.1):
    """
    计算结构中每个原子的 Q-score

    参数:
        grid: gemmi.FloatGrid，实验密度图
        structure: gemmi.Structure，原子模型
        sigma: 参考高斯宽度 (Å)，经验最优 0.6
        n_directions: 径向采样方向数（Fibonacci sphere）
        max_radius: 最大采样半径 (Å)
        step_size: 径向距离步长 (Å)
    返回:
        qscores: dict[(chain_name, res_seq, atom_name) → float]
        residue_qscores: dict[(chain_name, res_seq) → float] (残基中位数)
    """
    # 生成采样方向和距离
    directions = _fibonacci_sphere(n_directions)
    distances = np.arange(0, max_radius + step_size * 0.5, step_size)
    n_dist = len(distances)

    # 参考 profile
    ref_profile = _reference_gaussian(distances, sigma)

    # 预计算所有 (direction, distance) 偏移向量: (n_directions * n_dist, 3)
    # 对每个方向，生成该方向上所有距离的偏移
    offsets = directions[:, np.newaxis, :] * distances[np.newaxis, :, np.newaxis]
    # offsets shape: (n_directions, n_dist, 3)

    model = structure[0]
    qscores = {}
    residue_atoms = {}  # (chain, resseq) -> [qscore, ...]

    n_atoms_total = 0
    n_atoms_valid = 0

    for chain in model:
        for residue in chain:
            if residue.name == "HOH":
                continue
            for atom in residue:
                n_atoms_total += 1
                pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])

                # 对每个方向采样密度 profile
                profiles = np.zeros(n_dist)
                for d_idx in range(n_directions):
                    for r_idx in range(n_dist):
                        sample_pos = pos + offsets[d_idx, r_idx]
                        gpos = gemmi.Position(*sample_pos)
                        val = grid.interpolate_value(gpos)
                        profiles[r_idx] += val

                # 方向平均
                profiles /= n_directions

                # 计算 Pearson 相关
                if np.std(profiles) > 1e-10:
                    # 标准化
                    p_mean = profiles.mean()
                    r_mean = ref_profile.mean()
                    p_std = profiles.std()
                    r_std = ref_profile.std()
                    corr = np.sum((profiles - p_mean) * (ref_profile - r_mean)) / (
                        n_dist * p_std * r_std
                    )
                    corr = float(np.clip(corr, -1.0, 1.0))
                    n_atoms_valid += 1
                else:
                    corr = -1.0

                key = (chain.name, residue.seqid.num, atom.name)
                qscores[key] = corr

                res_key = (chain.name, residue.seqid.num)
                if res_key not in residue_atoms:
                    residue_atoms[res_key] = []
                residue_atoms[res_key].append(corr)

    # 残基级 Q-score（中位数）
    residue_qscores = {}
    for res_key, values in residue_atoms.items():
        residue_qscores[res_key] = float(np.median(values))

    logger.info(f"  Q-score: {n_atoms_valid}/{n_atoms_total} 原子有效, "
                f"采样方向={n_directions}, 步长={step_size}Å, σ={sigma}Å")

    return qscores, residue_qscores


def compute_qscore_summary(qscores):
    """
    计算 Q-score 汇总统计

    参数:
        qscores: compute_qscore_per_atom 返回的原子级 Q-score dict
    返回:
        summary: dict 包含 mean, median, std, q25, q75
    """
    if not qscores:
        return {"q_score_mean": 0.0, "q_score_median": 0.0,
                "q_score_std": 0.0, "q_score_q25": 0.0, "q_score_q75": 0.0,
                "q_score_n_atoms": 0}

    values = np.array(list(qscores.values()))
    return {
        "q_score_mean": float(np.mean(values)),
        "q_score_median": float(np.median(values)),
        "q_score_std": float(np.std(values)),
        "q_score_q25": float(np.percentile(values, 25)),
        "q_score_q75": float(np.percentile(values, 75)),
        "q_score_n_atoms": len(values),
    }
