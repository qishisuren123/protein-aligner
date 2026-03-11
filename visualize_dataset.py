#!/usr/bin/env python
"""
数据集统计可视化 V3 — 生成发表级数据集概览图

V3 变更：
- 移除分子类型饼图（蛋白专注）
- 新增：Domain 数量分布
- 新增：Interface 比例分布
- 更新覆盖率为 segment/CA-only

生成 6 个子图：
1. 分辨率分布直方图
2. CC_mask vs 分辨率散点图（按 tier 着色）
3. Domain 数量分布
4. 覆盖率分布（Segment vs CA-only）
5. Q-score 分布
6. Interface 比例分布

输出: data/output/dataset_overview.png
"""
import os
import sys
import json
import glob
import argparse
import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
except ImportError:
    print("需要安装 matplotlib: pip install matplotlib")
    sys.exit(1)


def collect_data(raw_dir):
    """从所有条目收集统计数据"""
    entries = []
    pattern = os.path.join(raw_dir, "EMD-*", "qc_metrics.json")
    for qc_path in sorted(glob.glob(pattern)):
        entry_dir = os.path.dirname(qc_path)
        entry_name = os.path.basename(entry_dir)

        entry = {"name": entry_name, "dir": entry_dir}

        # QC 指标
        with open(qc_path) as f:
            qc = json.load(f)
        entry["resolution"] = qc.get("resolution")
        entry["cc_mask"] = qc.get("cc_mask")
        entry["cc_volume"] = qc.get("cc_volume")
        entry["q_score_mean"] = qc.get("q_score_mean")
        entry["q_score_median"] = qc.get("q_score_median")

        # V3 标注统计
        label_path = os.path.join(entry_dir, "labeling_stats.json")
        if os.path.exists(label_path):
            with open(label_path) as f:
                stats = json.load(f)
            entry["coverage_segment_pct"] = stats.get("coverage_segment_pct", 0)
            entry["coverage_ca_pct"] = stats.get("coverage_ca_pct", 0)
            entry["n_interface_voxels"] = stats.get("n_interface_voxels", 0)
            entry["n_unique_domains"] = stats.get("n_unique_domains", 0)
            entry["n_voxels_segment"] = stats.get("n_voxels_segment", 0)
            entry["n_voxels_total"] = stats.get("n_voxels_total", 1)

        # Domain 信息
        dom_path = os.path.join(entry_dir, "domain_assignment.json")
        if os.path.exists(dom_path):
            with open(dom_path) as f:
                dom = json.load(f)
            entry["n_domains"] = dom.get("n_domains", 0)
            entry["domain_method"] = dom.get("method", "none")

        entries.append(entry)

    return entries


def classify_tier(entry):
    """根据质量指标分层（V3.1: 4 层）"""
    res = entry.get("resolution", 99)
    cc = entry.get("cc_mask", 0)
    qs = entry.get("q_score_mean", 0)

    if res is not None and cc is not None:
        if res < 2.5 and cc > 0.80 and (qs or 0) > 0.60:
            return "Gold"
        elif res < 4.0 and cc > 0.70 and (qs or 0) > 0.40:
            return "Silver"
        elif cc > 0.60 and (qs or 0) > 0.20:
            return "Copper"
    return "Hard"


def plot_dataset_overview(entries, output_path):
    """生成 V3 发表级数据集概览图"""
    if not entries:
        print("没有数据可以绘制")
        return

    # 设置发表级样式
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.dpi': 150,
    })

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # 颜色方案
    tier_colors = {"Gold": "#FFD700", "Silver": "#C0C0C0", "Copper": "#B87333", "Hard": "#CD7F32"}

    # === 1. 分辨率分布直方图 ===
    ax1 = fig.add_subplot(gs[0, 0])
    resolutions = [e["resolution"] for e in entries if e.get("resolution") is not None]
    if resolutions:
        ax1.hist(resolutions, bins=max(3, len(resolutions) // 2),
                 color="#4C72B0", edgecolor="white", alpha=0.8)
        ax1.set_xlabel("Resolution (Å)")
        ax1.set_ylabel("Count")
        ax1.set_title("Resolution Distribution")
        ax1.axvline(np.median(resolutions), color='red', linestyle='--',
                     label=f"Median={np.median(resolutions):.1f}Å")
        ax1.legend()

    # === 2. CC_mask vs 分辨率（按 tier 着色） ===
    ax2 = fig.add_subplot(gs[0, 1])
    for entry in entries:
        tier = classify_tier(entry)
        res = entry.get("resolution")
        cc = entry.get("cc_mask")
        if res is not None and cc is not None:
            ax2.scatter(res, cc, c=tier_colors[tier], s=80, edgecolors='black',
                        linewidth=0.5, label=tier, zorder=3)
    # 去重图例
    handles_labels = ax2.get_legend_handles_labels()
    if handles_labels[0]:
        unique = {}
        for h, l in zip(*handles_labels):
            if l not in unique:
                unique[l] = h
        ax2.legend(unique.values(), unique.keys(), title="Tier")
    ax2.set_xlabel("Resolution (Å)")
    ax2.set_ylabel("CC_mask")
    ax2.set_title("CC_mask vs Resolution")
    ax2.grid(True, alpha=0.3)

    # === 3. Domain 数量分布（V3 新增，替代分子类型饼图） ===
    ax3 = fig.add_subplot(gs[0, 2])
    n_domains_list = [e.get("n_domains", 0) for e in entries]
    names_short = [e["name"].replace("EMD-", "").replace("_", "\n") for e in entries]
    if names_short:
        colors_dom = ['seagreen' if d > 0 else 'lightgray' for d in n_domains_list]
        ax3.bar(range(len(names_short)), n_domains_list, color=colors_dom, edgecolor='white')
        ax3.set_xlabel("Entry")
        ax3.set_ylabel("# Domains")
        ax3.set_title("Domain Count per Entry (Merizo)")
        ax3.set_xticks(range(len(names_short)))
        ax3.set_xticklabels(names_short, fontsize=7, rotation=45, ha='right')
        for i, v in enumerate(n_domains_list):
            ax3.text(i, v + 0.1, str(v), ha='center', va='bottom', fontsize=8)

    # === 4. 覆盖率分布（V3: Segment vs CA-only） ===
    ax4 = fig.add_subplot(gs[1, 0])
    seg_cov = [e.get("coverage_segment_pct", 0) for e in entries]
    ca_cov = [e.get("coverage_ca_pct", 0) for e in entries]
    names_cov = [e["name"].replace("EMD-", "").replace("_", "\n") for e in entries]
    if names_cov:
        x = np.arange(len(names_cov))
        width = 0.35
        ax4.bar(x - width/2, seg_cov, width, label="Segment (all atom)", color="#4C72B0", alpha=0.8)
        ax4.bar(x + width/2, ca_cov, width, label="CA-only", color="#DD8452", alpha=0.8)
        ax4.set_xlabel("Entry")
        ax4.set_ylabel("Coverage (%)")
        ax4.set_title("V3 Labeling Coverage")
        ax4.set_xticks(x)
        ax4.set_xticklabels(names_cov, fontsize=7, rotation=45, ha='right')
        ax4.legend()

    # === 5. Q-score 分布 ===
    ax5 = fig.add_subplot(gs[1, 1])
    q_means = [e.get("q_score_mean") for e in entries
               if e.get("q_score_mean") is not None]
    q_names = [e["name"].replace("EMD-", "").replace("_", "\n")
               for e in entries if e.get("q_score_mean") is not None]
    if q_means:
        colors_q = []
        for qm in q_means:
            if qm > 0.5:
                colors_q.append("#55A868")
            elif qm > 0.3:
                colors_q.append("#FFD700")
            else:
                colors_q.append("#C44E52")
        ax5.bar(range(len(q_means)), q_means, color=colors_q, edgecolor='white')
        ax5.set_xlabel("Entry")
        ax5.set_ylabel("Q-score (mean)")
        ax5.set_title("Q-score per Entry")
        ax5.set_xticks(range(len(q_names)))
        ax5.set_xticklabels(q_names, fontsize=7, rotation=45, ha='right')
        ax5.axhline(0.5, color='green', linestyle='--', alpha=0.5, label="Good (>0.5)")
        ax5.axhline(0.3, color='orange', linestyle='--', alpha=0.5, label="Fair (>0.3)")
        ax5.legend()
        ax5.set_ylim(0, 1)

    # === 6. Interface 比例分布（V3 新增） ===
    ax6 = fig.add_subplot(gs[1, 2])
    iface_pcts = []
    iface_names = []
    for e in entries:
        n_iface = e.get("n_interface_voxels", 0)
        n_seg = e.get("n_voxels_segment", 0)
        if n_seg > 0:
            iface_pcts.append(100.0 * n_iface / n_seg)
        else:
            iface_pcts.append(0)
        iface_names.append(e["name"].replace("EMD-", "").replace("_", "\n"))
    if iface_names:
        colors_if = ['#C44E52' if p > 0 else 'lightgray' for p in iface_pcts]
        ax6.bar(range(len(iface_names)), iface_pcts, color=colors_if, edgecolor='white')
        ax6.set_xlabel("Entry")
        ax6.set_ylabel("Interface Voxels (%)")
        ax6.set_title("Protein-Protein Interface Ratio")
        ax6.set_xticks(range(len(iface_names)))
        ax6.set_xticklabels(iface_names, fontsize=7, rotation=45, ha='right')
        for i, v in enumerate(iface_pcts):
            if v > 0:
                ax6.text(i, v + 0.5, f'{v:.1f}%', ha='center', va='bottom', fontsize=7)

    # 总标题
    fig.suptitle("Cryo-EM Dataset Overview (V3)", fontsize=14, fontweight='bold', y=0.98)

    # 保存
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"V3 数据集概览图已保存: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="V3 数据集统计概览图")
    parser.add_argument("--raw-dir", default=os.path.join(PROJECT_ROOT, "data/raw"),
                        help="原始数据目录")
    parser.add_argument("--output", default=os.path.join(PROJECT_ROOT, "data/output/dataset_overview.png"),
                        help="输出图片路径")
    args = parser.parse_args()

    entries = collect_data(args.raw_dir)
    print(f"收集到 {len(entries)} 个条目")
    if entries:
        plot_dataset_overview(entries, args.output)
    else:
        print("未找到任何条目数据")


if __name__ == '__main__':
    main()
