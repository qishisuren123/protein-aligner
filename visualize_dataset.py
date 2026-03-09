#!/usr/bin/env python
"""
数据集统计可视化 — 生成发表级数据集概览图

生成 5 个子图：
1. 分辨率分布直方图
2. CC_mask vs 分辨率散点图（按 tier 着色）
3. 分子类型饼图
4. 覆盖率分布（Hard vs Voronoi）
5. Q-score 分布

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

        # 标注统计
        label_path = os.path.join(entry_dir, "labeling_stats.json")
        if os.path.exists(label_path):
            with open(label_path) as f:
                stats = json.load(f)
            entry["coverage_hard_pct"] = stats.get("coverage_hard_pct", 0)
            entry["coverage_total_pct"] = stats.get("coverage_total_pct", 0)
            entry["moltype_residue_counts"] = stats.get("moltype_residue_counts", {})
            entry["moltype_voxel_counts"] = stats.get("moltype_voxel_counts", {})

        entries.append(entry)

    return entries


def classify_tier(entry):
    """根据质量指标分层"""
    res = entry.get("resolution", 99)
    cc = entry.get("cc_mask", 0)
    qs = entry.get("q_score_mean", 0)

    if res is not None and cc is not None:
        if res < 2.5 and cc > 0.7 and (qs or 0) > 0.5:
            return "Gold"
        elif res < 3.5 and cc > 0.5:
            return "Silver"
    return "Hard"


def classify_moltype(entry):
    """确定条目的主要分子类型"""
    counts = entry.get("moltype_residue_counts", {})
    has_protein = counts.get("protein", 0) > 0
    has_rna = counts.get("RNA", 0) > 0
    has_dna = counts.get("DNA", 0) > 0
    has_sugar = counts.get("sugar", 0) > 0

    if has_rna and has_protein:
        return "Protein-RNA"
    if has_dna and has_protein:
        return "Protein-DNA"
    if has_sugar and has_protein:
        return "Glycoprotein"
    if has_protein:
        return "Protein"
    if has_rna:
        return "RNA"
    return "Other"


def plot_dataset_overview(entries, output_path):
    """生成发表级数据集概览图"""
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
    tier_colors = {"Gold": "#FFD700", "Silver": "#C0C0C0", "Hard": "#CD7F32"}
    moltype_colors = {
        "Protein": "#4C72B0", "Protein-RNA": "#DD8452",
        "Protein-DNA": "#55A868", "Glycoprotein": "#C44E52",
        "RNA": "#8172B3", "Other": "#937860"
    }

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

    # === 3. 分子类型饼图 ===
    ax3 = fig.add_subplot(gs[0, 2])
    mol_counts = {}
    for entry in entries:
        mt = classify_moltype(entry)
        mol_counts[mt] = mol_counts.get(mt, 0) + 1
    if mol_counts:
        labels = list(mol_counts.keys())
        sizes = list(mol_counts.values())
        colors = [moltype_colors.get(l, "#999999") for l in labels]
        wedges, texts, autotexts = ax3.pie(
            sizes, labels=labels, colors=colors, autopct='%1.0f%%',
            startangle=90, pctdistance=0.85
        )
        for t in autotexts:
            t.set_fontsize(9)
        ax3.set_title("Molecule Types")

    # === 4. 覆盖率分布（Hard vs Total） ===
    ax4 = fig.add_subplot(gs[1, 0])
    hard_cov = [e.get("coverage_hard_pct", 0) for e in entries
                if e.get("coverage_hard_pct") is not None]
    total_cov = [e.get("coverage_total_pct", 0) for e in entries
                 if e.get("coverage_total_pct") is not None]
    names = [e["name"].replace("EMD-", "").replace("_", "\n")
             for e in entries if e.get("coverage_total_pct") is not None]
    if names:
        x = np.arange(len(names))
        width = 0.35
        ax4.bar(x - width/2, hard_cov, width, label="Hard", color="#4C72B0", alpha=0.8)
        ax4.bar(x + width/2, total_cov, width, label="Total (+ Voronoi)",
                color="#DD8452", alpha=0.8)
        ax4.set_xlabel("Entry")
        ax4.set_ylabel("Coverage (%)")
        ax4.set_title("Labeling Coverage")
        ax4.set_xticks(x)
        ax4.set_xticklabels(names, fontsize=7, rotation=45, ha='right')
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

    # === 6. 数据分层统计 ===
    ax6 = fig.add_subplot(gs[1, 2])
    tier_counts = {"Gold": 0, "Silver": 0, "Hard": 0}
    for entry in entries:
        tier = classify_tier(entry)
        tier_counts[tier] += 1
    labels_t = list(tier_counts.keys())
    sizes_t = list(tier_counts.values())
    colors_t = [tier_colors[l] for l in labels_t]
    bars = ax6.bar(labels_t, sizes_t, color=colors_t, edgecolor='black', linewidth=0.5)
    for bar, val in zip(bars, sizes_t):
        ax6.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                 str(val), ha='center', va='bottom', fontweight='bold')
    ax6.set_ylabel("Count")
    ax6.set_title("Data Quality Tiers")

    # 总标题
    fig.suptitle("Cryo-EM Dataset Overview", fontsize=14, fontweight='bold', y=0.98)

    # 保存
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"数据集概览图已保存: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="生成数据集统计概览图")
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
