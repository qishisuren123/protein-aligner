"""
标签可视化脚本

生成密度图与各层标注的切片对比图，用于目视检查标注质量
支持输出：
  1. 综合面板：密度图 + 5层标签在同一切面的对比
  2. 多切面概览：沿3个轴的多个切面
  3. 标注覆盖率统计
  4. Voronoi vs Hard 标注对比
  5. Before/After 覆盖率对比
"""
import os
import sys
import argparse
import numpy as np
import gemmi
import json
import matplotlib
matplotlib.use('Agg')  # 无头模式
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.gridspec as gridspec

# 标签颜色映射
ATOM_NAMES = {0: 'bg', 1: 'CA', 2: 'N', 3: 'C', 4: 'O', 5: 'CB', 6: 'other'}
AA_NAMES = {
    0: 'bg', 1: 'ALA', 2: 'ARG', 3: 'ASN', 4: 'ASP', 5: 'CYS',
    6: 'GLN', 7: 'GLU', 8: 'GLY', 9: 'HIS', 10: 'ILE',
    11: 'LEU', 12: 'LYS', 13: 'MET', 14: 'PHE', 15: 'PRO',
    16: 'SER', 17: 'THR', 18: 'TRP', 19: 'TYR', 20: 'VAL'
}
SS_NAMES = {0: 'bg', 1: 'coil', 2: 'helix', 3: 'strand'}


def load_mrc(path):
    """读取 MRC 文件为 numpy 数组"""
    ccp4 = gemmi.read_ccp4_map(path)
    ccp4.setup(float('nan'))
    grid = ccp4.grid
    data = np.array(grid, copy=True)
    return data, grid


def find_signal_center(density):
    """找到密度信号的质心位置，用于选择最有内容的切面"""
    threshold = np.percentile(density[density > 0], 50) if (density > 0).any() else 0
    mask = density > threshold
    if not mask.any():
        return [s // 2 for s in density.shape]
    coords = np.argwhere(mask)
    center = coords.mean(axis=0).astype(int)
    return center


def make_discrete_cmap(n, name='tab20'):
    """创建离散颜色映射，0为透明背景"""
    base = plt.colormaps.get_cmap(name).resampled(max(n, 2))
    colors = [base(i) for i in range(n)]
    # 背景设为透明黑色
    colors[0] = (0, 0, 0, 0)
    return ListedColormap(colors)


def plot_comprehensive_panel(entry_dir, output_path=None):
    """
    综合面板：密度图 + 5层标签 + mol_map 在3个正交切面的对比
    """
    # 加载所有数据
    density, grid = load_mrc(os.path.join(entry_dir, "map_normalized.mrc"))
    label_atom, _ = load_mrc(os.path.join(entry_dir, "label_atom.mrc"))
    label_aa, _ = load_mrc(os.path.join(entry_dir, "label_aa.mrc"))
    label_ss, _ = load_mrc(os.path.join(entry_dir, "label_ss.mrc"))
    label_chain, _ = load_mrc(os.path.join(entry_dir, "label_chain.mrc"))
    label_conf, _ = load_mrc(os.path.join(entry_dir, "label_confidence.mrc"))

    mol_map_path = os.path.join(entry_dir, "mol_map.mrc")
    has_mol_map = os.path.exists(mol_map_path)
    if has_mol_map:
        mol_map, _ = load_mrc(mol_map_path)

    sim_map_path = os.path.join(entry_dir, "sim_map.mrc")
    has_sim_map = os.path.exists(sim_map_path)
    if has_sim_map:
        sim_map, _ = load_mrc(sim_map_path)

    # 找到信号中心
    center = find_signal_center(density)
    entry_name = os.path.basename(entry_dir)

    # 标签取整
    label_atom = np.round(label_atom).astype(int)
    label_aa = np.round(label_aa).astype(int)
    label_ss = np.round(label_ss).astype(int)
    label_chain = np.round(label_chain).astype(int)

    # === 图1：XY 切面综合面板 ===
    n_extra = int(has_mol_map) + int(has_sim_map)
    ncols = 6 + n_extra
    fig, axes = plt.subplots(3, ncols, figsize=(3.2 * ncols, 10))
    fig.suptitle(f'{entry_name} - Label Visualization\n'
                 f'Grid: {density.shape}, Spacing: {grid.spacing[0]:.3f}A',
                 fontsize=14, fontweight='bold')

    slice_labels = ['XY (z={z})', 'XZ (y={y})', 'YZ (x={x})']
    slices_density = [
        density[:, :, center[2]],
        density[:, center[1], :],
        density[center[0], :, :],
    ]
    slices_atom = [
        label_atom[:, :, center[2]],
        label_atom[:, center[1], :],
        label_atom[center[0], :, :],
    ]
    slices_aa = [
        label_aa[:, :, center[2]],
        label_aa[:, center[1], :],
        label_aa[center[0], :, :],
    ]
    slices_ss = [
        label_ss[:, :, center[2]],
        label_ss[:, center[1], :],
        label_ss[center[0], :, :],
    ]
    slices_chain = [
        label_chain[:, :, center[2]],
        label_chain[:, center[1], :],
        label_chain[center[0], :, :],
    ]
    slices_conf = [
        label_conf[:, :, center[2]],
        label_conf[:, center[1], :],
        label_conf[center[0], :, :],
    ]

    if has_mol_map:
        slices_mol = [
            mol_map[:, :, center[2]],
            mol_map[:, center[1], :],
            mol_map[center[0], :, :],
        ]
    if has_sim_map:
        slices_sim = [
            sim_map[:, :, center[2]],
            sim_map[:, center[1], :],
            sim_map[center[0], :, :],
        ]

    atom_cmap = make_discrete_cmap(7, 'tab10')
    ss_cmap = make_discrete_cmap(4, 'Set1')
    n_chains = int(label_chain.max()) + 1
    chain_cmap = make_discrete_cmap(max(n_chains, 2), 'tab10')
    n_aa = int(label_aa.max()) + 1
    aa_cmap = make_discrete_cmap(max(n_aa, 2), 'tab20')

    for row in range(3):
        col = 0

        # 密度图
        im = axes[row, col].imshow(slices_density[row].T, cmap='gray',
                                    origin='lower', aspect='equal')
        if row == 0:
            axes[row, col].set_title('Density', fontsize=10)
        axes[row, col].set_ylabel(slice_labels[row].format(
            x=center[0], y=center[1], z=center[2]), fontsize=9)
        col += 1

        # 原子类型
        im = axes[row, col].imshow(slices_density[row].T, cmap='gray',
                                    origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_atom[row].T, cmap=atom_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=6)
        if row == 0:
            axes[row, col].set_title('Atom Type', fontsize=10)
        col += 1

        # 氨基酸类型
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_aa[row].T, cmap=aa_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=20)
        if row == 0:
            axes[row, col].set_title('Amino Acid', fontsize=10)
        col += 1

        # 二级结构
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_ss[row].T, cmap=ss_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=3)
        if row == 0:
            axes[row, col].set_title('Sec. Struct.', fontsize=10)
        col += 1

        # 链
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_chain[row].T, cmap=chain_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=max(n_chains - 1, 1))
        if row == 0:
            axes[row, col].set_title('Chain', fontsize=10)
        col += 1

        # 置信度
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.3)
        axes[row, col].imshow(slices_conf[row].T, cmap='hot',
                               origin='lower', aspect='equal',
                               alpha=0.8, vmin=0, vmax=1)
        if row == 0:
            axes[row, col].set_title('Confidence', fontsize=10)
        col += 1

        # mol_map（如果有）
        if has_mol_map:
            axes[row, col].imshow(slices_mol[row].T, cmap='inferno',
                                   origin='lower', aspect='equal')
            if row == 0:
                axes[row, col].set_title('mol_map\n(denoise label)', fontsize=10)
            col += 1

        # sim_map（如果有）
        if has_sim_map:
            axes[row, col].imshow(slices_sim[row].T, cmap='inferno',
                                   origin='lower', aspect='equal')
            if row == 0:
                axes[row, col].set_title('sim_map\n(QC simulated)', fontsize=10)
            col += 1

    # 去掉坐标轴刻度
    for ax in axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()

    if output_path is None:
        output_path = os.path.join(entry_dir, "vis_labels_panel.png")
    fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Panel saved: {output_path}")

    # === 图2：多切面概览 ===
    plot_multi_slices(density, label_atom, label_ss, label_conf, center,
                      entry_name, entry_dir)

    # === 图3：统计图 ===
    plot_label_stats(entry_dir, label_atom, label_aa, label_ss,
                     label_chain, label_conf, density, entry_name)

    # === 图4：Voronoi vs Hard 标注对比 ===
    plot_voronoi_comparison(density, label_atom, label_conf, center,
                            entry_name, entry_dir)


def plot_multi_slices(density, label_atom, label_ss, label_conf,
                      center, entry_name, entry_dir):
    """沿 Z 轴等间距取多个切面，展示标注分布"""
    nz = density.shape[2]
    n_slices = 8
    # 选择有信号的范围
    z_signal = np.where(density.max(axis=(0, 1)) > 0.1)[0]
    if len(z_signal) > 0:
        z_start, z_end = z_signal[0], z_signal[-1]
    else:
        z_start, z_end = 0, nz - 1
    z_indices = np.linspace(z_start, z_end, n_slices, dtype=int)

    atom_cmap = make_discrete_cmap(7, 'tab10')

    fig, axes = plt.subplots(3, n_slices, figsize=(2.5 * n_slices, 8))
    fig.suptitle(f'{entry_name} - Z-axis Multi-slice (Density / Atom Label / Confidence)',
                 fontsize=13, fontweight='bold')

    for col, z in enumerate(z_indices):
        # 密度
        axes[0, col].imshow(density[:, :, z].T, cmap='gray',
                             origin='lower', aspect='equal')
        axes[0, col].set_title(f'z={z}', fontsize=9)

        # 原子标签 overlay
        axes[1, col].imshow(density[:, :, z].T, cmap='gray',
                             origin='lower', aspect='equal', alpha=0.4)
        axes[1, col].imshow(label_atom[:, :, z].T, cmap=atom_cmap,
                             origin='lower', aspect='equal',
                             alpha=0.8, vmin=0, vmax=6)

        # 置信度
        axes[2, col].imshow(label_conf[:, :, z].T, cmap='hot',
                             origin='lower', aspect='equal', vmin=0, vmax=1)

    for ax in axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])

    axes[0, 0].set_ylabel('Density', fontsize=10)
    axes[1, 0].set_ylabel('Atom Label', fontsize=10)
    axes[2, 0].set_ylabel('Confidence', fontsize=10)

    plt.tight_layout()
    out = os.path.join(entry_dir, "vis_multi_slices.png")
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Multi-slice overview saved: {out}")


def plot_label_stats(entry_dir, label_atom, label_aa, label_ss,
                     label_chain, label_conf, density, entry_name):
    """标注统计图"""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(f'{entry_name} - Label Statistics', fontsize=14, fontweight='bold')

    # 1. 原子类型分布
    atom_counts = {}
    for val, name in ATOM_NAMES.items():
        if val == 0:
            continue
        cnt = (label_atom == val).sum()
        if cnt > 0:
            atom_counts[name] = cnt
    if atom_counts:
        axes[0, 0].bar(atom_counts.keys(), atom_counts.values(), color='steelblue')
        axes[0, 0].set_title('Atom Type Distribution')
        axes[0, 0].set_ylabel('Voxel Count')
        for i, (k, v) in enumerate(atom_counts.items()):
            axes[0, 0].text(i, v, str(v), ha='center', va='bottom', fontsize=7)

    # 2. 氨基酸分布
    aa_counts = {}
    for val, name in AA_NAMES.items():
        if val == 0:
            continue
        cnt = (label_aa == val).sum()
        if cnt > 0:
            aa_counts[name] = cnt
    if aa_counts:
        axes[0, 1].bar(range(len(aa_counts)), aa_counts.values(),
                        color='coral', tick_label=list(aa_counts.keys()))
        axes[0, 1].set_title('Amino Acid Distribution')
        axes[0, 1].set_ylabel('Voxel Count')
        axes[0, 1].tick_params(axis='x', rotation=45, labelsize=7)

    # 3. 二级结构分布
    ss_counts = {}
    for val, name in SS_NAMES.items():
        if val == 0:
            continue
        cnt = (label_ss == val).sum()
        if cnt > 0:
            ss_counts[name] = cnt
    if ss_counts:
        colors_ss = {'coil': '#2196F3', 'helix': '#F44336', 'strand': '#4CAF50'}
        axes[0, 2].bar(ss_counts.keys(), ss_counts.values(),
                        color=[colors_ss.get(k, 'gray') for k in ss_counts.keys()])
        axes[0, 2].set_title('Secondary Structure Distribution')
        axes[0, 2].set_ylabel('Voxel Count')
        for i, (k, v) in enumerate(ss_counts.items()):
            axes[0, 2].text(i, v, str(v), ha='center', va='bottom', fontsize=8)

    # 4. 置信度直方图（区分 Hard vs Voronoi）
    conf_nonzero = label_conf[label_conf > 0]
    if len(conf_nonzero) > 0:
        # 高置信度区域（Hard标注）和低置信度区域（Voronoi标注）
        conf_hard = conf_nonzero[conf_nonzero > 0.3]
        conf_voronoi = conf_nonzero[conf_nonzero <= 0.3]

        axes[1, 0].hist(conf_hard.flatten(), bins=30, color='coral', alpha=0.8,
                         label=f'Hard (n={len(conf_hard)})')
        if len(conf_voronoi) > 0:
            axes[1, 0].hist(conf_voronoi.flatten(), bins=30, color='steelblue', alpha=0.6,
                             label=f'Voronoi (n={len(conf_voronoi)})')
        axes[1, 0].set_title(f'Confidence Distribution')
        axes[1, 0].set_xlabel('Confidence')
        axes[1, 0].set_ylabel('Voxel Count')
        axes[1, 0].legend(fontsize=8)

    # 5. 标注覆盖率对比
    total = density.size
    signal = (density > 0.05).sum()
    labeled = (label_atom > 0).sum()
    high_conf = (label_conf > 0.3).sum()
    low_conf = ((label_conf > 0) & (label_conf <= 0.3)).sum()

    labels_bar = ['Total', 'Signal\n(>0.05)', 'All Labeled', 'Hard\n(conf>0.3)', 'Voronoi\n(conf<=0.3)']
    values = [total, signal, labeled, high_conf, low_conf]
    colors = ['gray', 'steelblue', 'coral', '#FF6B35', '#4ECDC4']
    bars = axes[1, 1].bar(labels_bar, values, color=colors)
    axes[1, 1].set_title('Coverage Breakdown')
    axes[1, 1].set_ylabel('Voxel Count')
    for bar, val in zip(bars, values):
        axes[1, 1].text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                         f'{val}\n({100*val/total:.1f}%)',
                         ha='center', va='bottom', fontsize=7)

    # 6. 链标签分布（显示多链信息）
    chain_vals = label_chain[label_chain > 0]
    if len(chain_vals) > 0:
        unique_chains, chain_counts = np.unique(chain_vals, return_counts=True)
        if len(unique_chains) > 30:
            # 太多链时只显示前 30
            top_idx = np.argsort(-chain_counts)[:30]
            unique_chains = unique_chains[top_idx]
            chain_counts = chain_counts[top_idx]
        axes[1, 2].bar(range(len(unique_chains)), chain_counts, color='mediumpurple')
        axes[1, 2].set_title(f'Chain Distribution ({len(unique_chains)} chains)')
        axes[1, 2].set_xlabel('Chain ID')
        axes[1, 2].set_ylabel('Voxel Count')
        if len(unique_chains) <= 30:
            axes[1, 2].set_xticks(range(len(unique_chains)))
            axes[1, 2].set_xticklabels([str(c) for c in unique_chains], fontsize=6, rotation=45)
    else:
        axes[1, 2].text(0.5, 0.5, 'No chain labels', ha='center', va='center',
                          transform=axes[1, 2].transAxes)

    plt.tight_layout()
    out = os.path.join(entry_dir, "vis_label_stats.png")
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Label stats saved: {out}")


def plot_voronoi_comparison(density, label_atom, label_conf, center,
                             entry_name, entry_dir):
    """
    Voronoi vs Hard 标注对比图

    在同一切面上对比展示：
    - 仅 Hard 标注区域（高置信度）
    - 全部标注区域（Hard + Voronoi）
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'{entry_name} - Hard vs Voronoi Labeling Comparison',
                 fontsize=14, fontweight='bold')

    atom_cmap = make_discrete_cmap(7, 'tab10')

    # 三个切面
    slices_density = [
        density[:, :, center[2]],
        density[:, center[1], :],
        density[center[0], :, :],
    ]
    slices_atom = [
        label_atom[:, :, center[2]],
        label_atom[:, center[1], :],
        label_atom[center[0], :, :],
    ]
    slices_conf = [
        label_conf[:, :, center[2]],
        label_conf[:, center[1], :],
        label_conf[center[0], :, :],
    ]

    slice_names = [f'XY (z={center[2]})', f'XZ (y={center[1]})', f'YZ (x={center[0]})']

    for col in range(3):
        # 第一行：仅 Hard 标注（置信度 > 0.3）
        hard_mask = (slices_conf[col] > 0.3).astype(float)
        hard_atoms = slices_atom[col] * (hard_mask > 0).astype(int)

        axes[0, col].imshow(slices_density[col].T, cmap='gray',
                             origin='lower', aspect='equal', alpha=0.5)
        axes[0, col].imshow(hard_atoms.T, cmap=atom_cmap,
                             origin='lower', aspect='equal',
                             alpha=0.8, vmin=0, vmax=6)
        if col == 0:
            axes[0, col].set_ylabel('Hard Only\n(conf>0.3)', fontsize=11)
        axes[0, col].set_title(slice_names[col], fontsize=10)

        # 第二行：全部标注（Hard + Voronoi）
        all_atoms = slices_atom[col]
        axes[1, col].imshow(slices_density[col].T, cmap='gray',
                             origin='lower', aspect='equal', alpha=0.5)
        axes[1, col].imshow(all_atoms.T, cmap=atom_cmap,
                             origin='lower', aspect='equal',
                             alpha=0.8, vmin=0, vmax=6)
        if col == 0:
            axes[1, col].set_ylabel('Hard + Voronoi\n(all labeled)', fontsize=11)

    for ax in axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    out = os.path.join(entry_dir, "vis_voronoi_comparison.png")
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Voronoi comparison saved: {out}")


def main():
    parser = argparse.ArgumentParser(description='标签可视化')
    parser.add_argument('entry_dir', nargs='?',
                        help='条目目录路径，不指定则处理所有')
    parser.add_argument('--data-dir', default='data/raw',
                        help='数据根目录 (默认: data/raw)')
    args = parser.parse_args()

    if args.entry_dir:
        dirs = [args.entry_dir]
    else:
        base = args.data_dir
        dirs = [os.path.join(base, d) for d in sorted(os.listdir(base))
                if os.path.isdir(os.path.join(base, d))]

    for entry_dir in dirs:
        if not os.path.exists(os.path.join(entry_dir, "map_normalized.mrc")):
            print(f"跳过 {entry_dir}: 无 map_normalized.mrc")
            continue
        if not os.path.exists(os.path.join(entry_dir, "label_atom.mrc")):
            print(f"跳过 {entry_dir}: 无标签文件")
            continue

        print(f"\n{'='*60}")
        print(f"可视化: {entry_dir}")
        print(f"{'='*60}")
        plot_comprehensive_panel(entry_dir)


if __name__ == "__main__":
    main()
