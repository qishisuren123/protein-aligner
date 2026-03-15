"""
Q-score 物理置信度计算 (Pintilie et al., 2020, Nature Methods)

Q-score 是 cryo-EM 领域标准的原子级密度拟合质量指标。
对每个原子，沿径向方向采样实验密度，与参考高斯做 Pearson 相关。
返回值范围 [-1, 1]，越接近 1 表示局部密度与原子位置匹配越好。

关键改进（V3.4）：
- 每个半径 r 上只使用"比任何其他原子都更接近当前原子"的球面采样点
- 每个半径保留 N=8 个有效点
- r=0 使用原子中心密度值（重复 N 次）
- 残基聚合使用 mean（非 median）

参考: Pintilie et al., Nature Methods 17, 328–334 (2020)
"""
import logging
import numpy as np
import gemmi
from scipy.spatial import cKDTree

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
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    theta = golden_angle * indices
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


def compute_qscore_per_atom(grid, structure, sigma=0.6, n_directions=64,
                            max_radius=2.0, step_size=0.1, n_points_per_shell=8):
    """
    计算结构中每个原子的 Q-score（Pintilie et al. 2020）

    核心改进：每个半径 r 上只使用离当前原子最近（比任何其他原子都近）的
    球面采样点，避免邻近原子密度对径向 profile 的污染。

    参数:
        grid: gemmi.FloatGrid，实验密度图
        structure: gemmi.Structure，原子模型
        sigma: 参考高斯宽度 (Å)，经验最优 0.6
        n_directions: 候选采样方向数（Fibonacci sphere），需 > n_points_per_shell
        max_radius: 最大采样半径 (Å)
        step_size: 径向距离步长 (Å)
        n_points_per_shell: 每个半径壳层保留的有效采样点数 (N=8)
    返回:
        qscores: dict[(chain_name, res_seq, atom_name) → float]
        residue_qscores: dict[(chain_name, res_seq) → float] (残基均值)
    """
    # 生成采样方向和距离
    directions = _fibonacci_sphere(n_directions)
    distances = np.arange(0, max_radius + step_size * 0.5, step_size)
    n_dist = len(distances)
    N = n_points_per_shell

    # 参考 profile
    ref_profile = _reference_gaussian(distances, sigma)

    # 预计算偏移向量: (n_directions, n_dist, 3)
    offsets = directions[:, np.newaxis, :] * distances[np.newaxis, :, np.newaxis]

    # 收集所有原子坐标，构建 KDTree
    model = structure[0]
    all_atoms = []  # [(chain, resseq, atomname, pos_array), ...]
    all_positions = []

    for chain in model:
        for residue in chain:
            if residue.name == "HOH":
                continue
            for atom in residue:
                pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                all_atoms.append((chain.name, residue.seqid.num, atom.name, pos))
                all_positions.append(pos)

    if not all_positions:
        return {}, {}

    all_positions = np.array(all_positions)
    atom_tree = cKDTree(all_positions)

    qscores = {}
    residue_atoms = {}  # (chain, resseq) -> [qscore, ...]

    n_atoms_total = len(all_atoms)
    n_atoms_valid = 0

    for atom_idx, (chain_name, resseq, atomname, atom_pos) in enumerate(all_atoms):
        # 构建径向 profile（只用 closest-to-this-atom 的采样点）
        profile = np.zeros(n_dist)

        for r_idx in range(n_dist):
            r = distances[r_idx]

            if r == 0:
                # r=0: 原子中心密度值，重复 N 次取均值 = 中心密度值
                gpos = gemmi.Position(*atom_pos)
                profile[r_idx] = grid.interpolate_value(gpos)
                continue

            # 生成球面采样点: atom_pos + r * direction
            sample_points = atom_pos + offsets[:, r_idx, :]  # (n_directions, 3)

            # 查询每个采样点的最近原子
            _, nearest_indices = atom_tree.query(sample_points, k=1)

            # 只保留最近原子是当前原子的采样点
            valid_mask = (nearest_indices == atom_idx)
            valid_points = sample_points[valid_mask]

            if len(valid_points) == 0:
                # 没有有效点，此半径的密度设为 0（被其他原子完全遮挡）
                profile[r_idx] = 0.0
                continue

            # 取 N 个有效点（如果多于 N，取前 N 个）
            if len(valid_points) > N:
                valid_points = valid_points[:N]

            # 采样密度并平均
            density_vals = np.zeros(len(valid_points))
            for p_idx, pt in enumerate(valid_points):
                gpos = gemmi.Position(*pt)
                density_vals[p_idx] = grid.interpolate_value(gpos)

            profile[r_idx] = density_vals.mean()

        # 计算 Pearson 相关（Q-score）
        if np.std(profile) > 1e-10:
            p_mean = profile.mean()
            r_mean = ref_profile.mean()
            p_std = profile.std()
            r_std = ref_profile.std()
            corr = np.sum((profile - p_mean) * (ref_profile - r_mean)) / (
                n_dist * p_std * r_std
            )
            corr = float(np.clip(corr, -1.0, 1.0))
            n_atoms_valid += 1
        else:
            corr = -1.0

        key = (chain_name, resseq, atomname)
        qscores[key] = corr

        res_key = (chain_name, resseq)
        if res_key not in residue_atoms:
            residue_atoms[res_key] = []
        residue_atoms[res_key].append(corr)

    # 残基级 Q-score（均值，非中位数）
    residue_qscores = {}
    for res_key, values in residue_atoms.items():
        residue_qscores[res_key] = float(np.mean(values))

    logger.info(f"  Q-score: {n_atoms_valid}/{n_atoms_total} 原子有效, "
                f"方向={n_directions}, N/shell={N}, 步长={step_size}Å, σ={sigma}Å")

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
