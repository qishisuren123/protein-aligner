"""
标签可视化脚本 V3

生成密度图与 V3 各层标注的切片对比图，用于目视检查标注质量
V3 标签通道：Density / Segment / Atom / AA / SS / Q-score / Chain / Domain / Interface
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

# V3 标签名映射
ATOM_NAMES = {0: 'bg', 1: 'CA', 2: 'C', 3: 'N', 4: 'O', 5: 'CB', 6: 'other'}
AA_NAMES = {
    0: 'bg', 1: 'ALA', 2: 'ARG', 3: 'ASN', 4: 'ASP', 5: 'CYS',
    6: 'GLN', 7: 'GLU', 8: 'GLY', 9: 'HIS', 10: 'ILE',
    11: 'LEU', 12: 'LYS', 13: 'MET', 14: 'PHE', 15: 'PRO',
    16: 'SER', 17: 'THR', 18: 'TRP', 19: 'TYR', 20: 'VAL'
}
# V3 二级结构重新编号
SS_NAMES = {0: 'bg', 1: 'helix', 2: 'strand', 3: 'coil'}


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
    V3 综合面板：密度图 + 8层标签在3个正交切面的对比
    通道：Density / Segment / Atom / AA / SS / Q-score / Chain / Domain / Interface
    """
    # 加载所有数据
    density, grid = load_mrc(os.path.join(entry_dir, "map_normalized.mrc"))

    # V3 标签通道
    label_segment, _ = load_mrc(os.path.join(entry_dir, "label_segment.mrc"))
    label_atom, _ = load_mrc(os.path.join(entry_dir, "label_atom.mrc"))
    label_aa, _ = load_mrc(os.path.join(entry_dir, "label_aa.mrc"))
    label_ss, _ = load_mrc(os.path.join(entry_dir, "label_ss.mrc"))
    label_qscore, _ = load_mrc(os.path.join(entry_dir, "label_qscore.mrc"))
    label_chain, _ = load_mrc(os.path.join(entry_dir, "label_chain.mrc"))

    # domain/interface 可选（可能全零）
    dom_path = os.path.join(entry_dir, "label_domain.mrc")
    label_domain, _ = load_mrc(dom_path) if os.path.exists(dom_path) else (np.zeros_like(density), None)

    iface_path = os.path.join(entry_dir, "label_interface.mrc")
    label_interface, _ = load_mrc(iface_path) if os.path.exists(iface_path) else (np.zeros_like(density), None)

    mol_map_path = os.path.join(entry_dir, "mol_map.mrc")
    has_mol_map = os.path.exists(mol_map_path)
    if has_mol_map:
        mol_map, _ = load_mrc(mol_map_path)

    # 找到信号中心
    center = find_signal_center(density)
    entry_name = os.path.basename(entry_dir)

    # 标签取整
    label_segment = np.round(label_segment).astype(int)
    label_atom = np.round(label_atom).astype(int)
    label_aa = np.round(label_aa).astype(int)
    label_ss = np.round(label_ss).astype(int)
    label_chain = np.round(label_chain).astype(int)
    label_domain = np.round(label_domain).astype(int)
    label_interface = np.round(label_interface).astype(int)

    # V3 面板：9 列基础 + mol_map 可选
    ncols = 9 + int(has_mol_map)
    fig, axes = plt.subplots(3, ncols, figsize=(3.0 * ncols, 10))
    fig.suptitle(f'{entry_name} - V3 Label Visualization\n'
                 f'Grid: {density.shape}, Spacing: {grid.spacing[0]:.3f}A',
                 fontsize=14, fontweight='bold')

    # 切片数据
    def get_slices(vol):
        return [vol[:, :, center[2]], vol[:, center[1], :], vol[center[0], :, :]]

    slices_density = get_slices(density)
    slices_segment = get_slices(label_segment)
    slices_atom = get_slices(label_atom)
    slices_aa = get_slices(label_aa)
    slices_ss = get_slices(label_ss)
    slices_qscore = get_slices(label_qscore)
    slices_chain = get_slices(label_chain)
    slices_domain = get_slices(label_domain)
    slices_interface = get_slices(label_interface)
    if has_mol_map:
        slices_mol = get_slices(mol_map)

    slice_labels = ['XY (z={z})', 'XZ (y={y})', 'YZ (x={x})']

    # 颜色映射
    seg_cmap = make_discrete_cmap(2, 'Set1')
    atom_cmap = make_discrete_cmap(7, 'tab10')
    n_chains = int(label_chain.max()) + 1
    chain_cmap = make_discrete_cmap(max(n_chains, 2), 'tab10')
    n_aa = int(label_aa.max()) + 1
    aa_cmap = make_discrete_cmap(max(n_aa, 2), 'tab20')
    ss_cmap = make_discrete_cmap(4, 'Set1')
    n_domains = int(label_domain.max()) + 1
    domain_cmap = make_discrete_cmap(max(n_domains, 2), 'Set3')
    iface_cmap = make_discrete_cmap(2, 'Set1')

    for row in range(3):
        col = 0

        # 1. 密度图
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal')
        if row == 0:
            axes[row, col].set_title('Density', fontsize=9)
        axes[row, col].set_ylabel(slice_labels[row].format(
            x=center[0], y=center[1], z=center[2]), fontsize=8)
        col += 1

        # 2. Segment（蛋白/背景）
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_segment[row].T, cmap=seg_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=1)
        if row == 0:
            axes[row, col].set_title('Segment', fontsize=9)
        col += 1

        # 3. 原子类型
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_atom[row].T, cmap=atom_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=6)
        if row == 0:
            axes[row, col].set_title('Atom Type', fontsize=9)
        col += 1

        # 4. 氨基酸类型
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_aa[row].T, cmap=aa_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=20)
        if row == 0:
            axes[row, col].set_title('Amino Acid', fontsize=9)
        col += 1

        # 5. 二级结构
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_ss[row].T, cmap=ss_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=3)
        if row == 0:
            axes[row, col].set_title('Sec. Struct.', fontsize=9)
        col += 1

        # 6. Q-score
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.3)
        axes[row, col].imshow(slices_qscore[row].T, cmap='RdYlGn',
                               origin='lower', aspect='equal',
                               alpha=0.8, vmin=-1, vmax=1)
        if row == 0:
            axes[row, col].set_title('Q-score', fontsize=9)
        col += 1

        # 7. 链
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_chain[row].T, cmap=chain_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=max(n_chains - 1, 1))
        if row == 0:
            axes[row, col].set_title('Chain', fontsize=9)
        col += 1

        # 8. 结构域
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_domain[row].T, cmap=domain_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=max(n_domains - 1, 1))
        if row == 0:
            axes[row, col].set_title('Domain', fontsize=9)
        col += 1

        # 9. 界面
        axes[row, col].imshow(slices_density[row].T, cmap='gray',
                               origin='lower', aspect='equal', alpha=0.5)
        axes[row, col].imshow(slices_interface[row].T, cmap=iface_cmap,
                               origin='lower', aspect='equal',
                               alpha=0.7, vmin=0, vmax=1)
        if row == 0:
            axes[row, col].set_title('Interface', fontsize=9)
        col += 1

        # 10. mol_map（如果有）
        if has_mol_map:
            axes[row, col].imshow(slices_mol[row].T, cmap='inferno',
                                   origin='lower', aspect='equal')
            if row == 0:
                axes[row, col].set_title('mol_map', fontsize=9)
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
    print(f"V3 Panel saved: {output_path}")

    # === 统计图 ===
    plot_label_stats(entry_dir, label_segment, label_atom, label_aa, label_ss,
                     label_qscore, label_chain, label_domain, label_interface,
                     density, entry_name)


def plot_label_stats(entry_dir, label_segment, label_atom, label_aa, label_ss,
                     label_qscore, label_chain, label_domain, label_interface,
                     density, entry_name):
    """V3 标注统计图"""
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    fig.suptitle(f'{entry_name} - V3 Label Statistics', fontsize=14, fontweight='bold')

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

    # 3. 二级结构分布（V3 编号：helix=1, strand=2, coil=3）
    ss_counts = {}
    for val, name in SS_NAMES.items():
        if val == 0:
            continue
        cnt = (label_ss == val).sum()
        if cnt > 0:
            ss_counts[name] = cnt
    if ss_counts:
        colors_ss = {'helix': '#F44336', 'strand': '#4CAF50', 'coil': '#2196F3'}
        axes[0, 2].bar(ss_counts.keys(), ss_counts.values(),
                        color=[colors_ss.get(k, 'gray') for k in ss_counts.keys()])
        axes[0, 2].set_title('Secondary Structure Distribution')
        axes[0, 2].set_ylabel('Voxel Count')
        for i, (k, v) in enumerate(ss_counts.items()):
            axes[0, 2].text(i, v, str(v), ha='center', va='bottom', fontsize=8)

    # 4. Q-score 直方图
    qs_nonzero = label_qscore[label_qscore != 0]
    if len(qs_nonzero) > 0:
        axes[0, 3].hist(qs_nonzero.flatten(), bins=50, color='teal', alpha=0.8)
        axes[0, 3].set_title(f'Q-score Distribution (n={len(qs_nonzero)})')
        axes[0, 3].set_xlabel('Q-score')
        axes[0, 3].set_ylabel('Voxel Count')
        axes[0, 3].axvline(np.median(qs_nonzero), color='red', linestyle='--',
                            label=f'median={np.median(qs_nonzero):.2f}')
        axes[0, 3].legend()
    else:
        axes[0, 3].text(0.5, 0.5, 'No Q-score data', ha='center', va='center',
                          transform=axes[0, 3].transAxes)

    # 5. 标注覆盖率对比
    total = density.size
    signal = (density > 0.05).sum()
    segment_labeled = (label_segment > 0).sum()
    ca_labeled = (label_aa > 0).sum()

    labels_bar = ['Total', 'Signal\n(>0.05)', 'Segment\n(protein)', 'CA-only\n(aa/ss/...)']
    values = [total, signal, segment_labeled, ca_labeled]
    colors = ['gray', 'steelblue', 'coral', '#FF6B35']
    bars = axes[1, 0].bar(labels_bar, values, color=colors)
    axes[1, 0].set_title('Coverage Breakdown')
    axes[1, 0].set_ylabel('Voxel Count')
    for bar, val in zip(bars, values):
        axes[1, 0].text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                         f'{val}\n({100*val/total:.1f}%)',
                         ha='center', va='bottom', fontsize=7)

    # 6. 链标签分布
    chain_vals = label_chain[label_chain > 0]
    if len(chain_vals) > 0:
        unique_chains, chain_counts = np.unique(chain_vals, return_counts=True)
        if len(unique_chains) > 30:
            top_idx = np.argsort(-chain_counts)[:30]
            unique_chains = unique_chains[top_idx]
            chain_counts = chain_counts[top_idx]
        axes[1, 1].bar(range(len(unique_chains)), chain_counts, color='mediumpurple')
        axes[1, 1].set_title(f'Chain Distribution ({len(unique_chains)} chains)')
        axes[1, 1].set_xlabel('Chain ID')
        axes[1, 1].set_ylabel('Voxel Count')
    else:
        axes[1, 1].text(0.5, 0.5, 'No chain labels', ha='center', va='center',
                          transform=axes[1, 1].transAxes)

    # 7. Domain 分布
    dom_vals = label_domain[label_domain > 0]
    if len(dom_vals) > 0:
        unique_doms, dom_counts = np.unique(dom_vals, return_counts=True)
        n_doms = len(unique_doms)
        if n_doms > 30:
            # domain 太多时只显示前30个（按体素数排序）
            top_idx = np.argsort(-dom_counts)[:30]
            unique_doms = unique_doms[top_idx]
            dom_counts = dom_counts[top_idx]
            # 按 domain ID 重新排序
            sort_idx = np.argsort(unique_doms)
            unique_doms = unique_doms[sort_idx]
            dom_counts = dom_counts[sort_idx]
            title_suffix = f' (top 30 of {n_doms})'
        else:
            title_suffix = ''
        axes[1, 2].bar(range(len(unique_doms)), dom_counts, color='seagreen')
        # 只在 domain 较少时显示 tick label
        if len(unique_doms) <= 20:
            axes[1, 2].set_xticks(range(len(unique_doms)))
            axes[1, 2].set_xticklabels([str(int(d)) for d in unique_doms], fontsize=7)
        axes[1, 2].set_title(f'Domain Distribution ({n_doms} domains){title_suffix}')
        axes[1, 2].set_xlabel('Domain ID')
        axes[1, 2].set_ylabel('Voxel Count')
    else:
        axes[1, 2].text(0.5, 0.5, 'No domain labels\n(Merizo unavailable?)',
                          ha='center', va='center', transform=axes[1, 2].transAxes)

    # 8. Interface 统计
    n_iface = int((label_interface > 0).sum())
    n_non_iface = int((label_interface == 0).sum()) - int((label_segment == 0).sum())
    n_non_iface = max(0, n_non_iface)
    if n_iface > 0 or n_non_iface > 0:
        labels_iface = ['Non-interface', 'Interface']
        vals_iface = [n_non_iface, n_iface]
        colors_iface = ['#4C72B0', '#C44E52']
        axes[1, 3].bar(labels_iface, vals_iface, color=colors_iface)
        axes[1, 3].set_title(f'Interface Detection')
        axes[1, 3].set_ylabel('Voxel Count')
        for i, v in enumerate(vals_iface):
            axes[1, 3].text(i, v, str(v), ha='center', va='bottom', fontsize=8)
    else:
        axes[1, 3].text(0.5, 0.5, 'No interface data\n(single chain?)',
                          ha='center', va='center', transform=axes[1, 3].transAxes)

    plt.tight_layout()
    out = os.path.join(entry_dir, "vis_label_stats.png")
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"V3 Label stats saved: {out}")


def main():
    parser = argparse.ArgumentParser(description='V3 标签可视化')
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
        if not os.path.exists(os.path.join(entry_dir, "label_segment.mrc")):
            print(f"跳过 {entry_dir}: 无 V3 标签文件")
            continue

        print(f"\n{'='*60}")
        print(f"V3 可视化: {entry_dir}")
        print(f"{'='*60}")
        plot_comprehensive_panel(entry_dir)


if __name__ == "__main__":
    main()
