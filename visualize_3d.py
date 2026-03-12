#!/usr/bin/env python
"""
3D 等值面可视化脚本 (V3.1)

使用 marching cubes 提取等值面 + plotly 交互式渲染：
- 离散标签（segment, atom, aa, ss, chain, domain, interface）：逐类别提取等值面，不同颜色
- 连续标签（qscore, mol_map）：等值面 + colorscale 着色

输出：
- vis_3d_labels.html — 交互式 3D 查看器
- vis_3d_labels.png — 静态截图（需要 kaleido）
"""
import os
import sys
import argparse
import numpy as np
import json

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

try:
    import gemmi
except ImportError:
    print("需要安装 gemmi: pip install gemmi")
    sys.exit(1)

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    print("需要安装 plotly: pip install plotly")
    sys.exit(1)

try:
    from skimage.measure import marching_cubes
except ImportError:
    print("需要安装 scikit-image: pip install scikit-image")
    sys.exit(1)

from scipy.ndimage import gaussian_filter


# 标签名称映射
ATOM_NAMES = {1: 'CA', 2: 'C', 3: 'N', 4: 'O', 5: 'CB', 6: 'other'}
AA_NAMES = {
    1: 'ALA', 2: 'ARG', 3: 'ASN', 4: 'ASP', 5: 'CYS',
    6: 'GLN', 7: 'GLU', 8: 'GLY', 9: 'HIS', 10: 'ILE',
    11: 'LEU', 12: 'LYS', 13: 'MET', 14: 'PHE', 15: 'PRO',
    16: 'SER', 17: 'THR', 18: 'TRP', 19: 'TYR', 20: 'VAL',
}
SS_NAMES = {1: 'helix', 2: 'strand', 3: 'coil'}

# 颜色方案（plotly 格式）
DISCRETE_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
    '#393b79', '#5254a3', '#6b6ecf', '#9c9ede',
]

SS_COLORS = {'helix': '#F44336', 'strand': '#4CAF50', 'coil': '#2196F3'}


def load_mrc(path):
    """读取 MRC 文件为 numpy 数组"""
    ccp4 = gemmi.read_ccp4_map(path)
    ccp4.setup(float('nan'))
    grid = ccp4.grid
    data = np.array(grid, copy=True)
    return data, grid


def downsample_volume(data, max_size=128):
    """
    下采样体数据到合理尺寸，避免 marching cubes 过慢
    对离散标签使用最近邻，对连续数据使用线性
    """
    shape = np.array(data.shape)
    if shape.max() <= max_size:
        return data, np.ones(3)

    scale = max_size / shape.max()
    from scipy.ndimage import zoom
    new_data = zoom(data, scale, order=0)  # 最近邻
    actual_scale = np.array(new_data.shape) / np.array(data.shape)
    return new_data, actual_scale


def mesh_from_discrete_label(label_data, label_val, sigma=1.0, level=0.5):
    """
    从离散标签中提取单个类别的等值面

    步骤：二值化 → 高斯模糊平滑 → marching cubes
    """
    binary = (label_data == label_val).astype(np.float32)
    if binary.sum() < 10:
        return None

    smoothed = gaussian_filter(binary, sigma=sigma)
    try:
        verts, faces, normals, _ = marching_cubes(smoothed, level=level)
        return verts, faces
    except (ValueError, RuntimeError):
        return None


def mesh_from_continuous(data, level, sigma=0.5):
    """
    从连续数据中提取等值面
    """
    if sigma > 0:
        data = gaussian_filter(data, sigma=sigma)
    try:
        verts, faces, normals, values = marching_cubes(data, level=level)
        return verts, faces, values
    except (ValueError, RuntimeError):
        return None


def make_mesh3d(verts, faces, color, name, opacity=0.6):
    """创建 plotly Mesh3d trace"""
    return go.Mesh3d(
        x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        color=color,
        opacity=opacity,
        name=name,
        showlegend=True,
        flatshading=True,
    )


def make_mesh3d_continuous(verts, faces, values, name, colorscale='Viridis', opacity=0.7):
    """创建带 colorscale 的 plotly Mesh3d trace（连续数据）"""
    return go.Mesh3d(
        x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        intensity=values,
        colorscale=colorscale,
        opacity=opacity,
        name=name,
        showlegend=True,
        flatshading=True,
        colorbar=dict(title=name, len=0.5),
    )


def visualize_discrete_label(label_data, label_names, title, colors=None, max_labels=24):
    """
    可视化离散标签的 3D 等值面

    返回 plotly Figure
    """
    data_ds, scale = downsample_volume(label_data)
    unique_vals = sorted(set(np.unique(data_ds).astype(int)) - {0})

    if len(unique_vals) > max_labels:
        # 取体素数最多的前 max_labels 个
        counts = {v: (data_ds == v).sum() for v in unique_vals}
        unique_vals = sorted(counts.keys(), key=lambda v: counts[v], reverse=True)[:max_labels]

    fig = go.Figure()
    for i, val in enumerate(unique_vals):
        result = mesh_from_discrete_label(data_ds, val, sigma=1.0, level=0.5)
        if result is None:
            continue
        verts, faces = result
        name = label_names.get(val, str(val))
        color = colors.get(val, DISCRETE_COLORS[i % len(DISCRETE_COLORS)]) if colors else DISCRETE_COLORS[i % len(DISCRETE_COLORS)]
        fig.add_trace(make_mesh3d(verts, faces, color, name))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X', yaxis_title='Y', zaxis_title='Z',
            aspectmode='data',
        ),
        legend=dict(itemsizing='constant'),
        margin=dict(l=0, r=0, t=40, b=0),
    )
    return fig


def visualize_continuous_label(data, title, threshold_pct=50, colorscale='RdYlGn'):
    """
    可视化连续标签的 3D 等值面

    返回 plotly Figure
    """
    data_ds, scale = downsample_volume(data.astype(np.float32), max_size=128)

    # 计算等值面阈值（非零值的百分位数）
    nonzero = data_ds[data_ds != 0]
    if len(nonzero) < 10:
        return go.Figure().update_layout(title=f"{title} (no data)")

    level = np.percentile(nonzero, threshold_pct)
    result = mesh_from_continuous(data_ds, level=level, sigma=0.5)
    if result is None:
        return go.Figure().update_layout(title=f"{title} (marching cubes failed)")

    verts, faces, values = result
    fig = go.Figure()
    fig.add_trace(make_mesh3d_continuous(verts, faces, values, title, colorscale=colorscale))

    fig.update_layout(
        title=f"{title} (level={level:.3f})",
        scene=dict(
            xaxis_title='X', yaxis_title='Y', zaxis_title='Z',
            aspectmode='data',
        ),
        margin=dict(l=0, r=0, t=40, b=0),
    )
    return fig


def create_combined_figure(entry_dir):
    """
    创建合并的 3D 可视化页面：多个标签的 3D 等值面

    返回组合 HTML 字符串
    """
    entry_name = os.path.basename(entry_dir)
    figures = {}

    # === 密度图（连续） ===
    density_path = os.path.join(entry_dir, "map_normalized.mrc")
    if os.path.exists(density_path):
        density, _ = load_mrc(density_path)
        figures['density'] = visualize_continuous_label(
            density, f"Density Map", threshold_pct=70, colorscale='Greys')

    # === Segment（离散，binary） ===
    seg_path = os.path.join(entry_dir, "label_segment.mrc")
    if os.path.exists(seg_path):
        seg_data, _ = load_mrc(seg_path)
        seg_data = np.round(seg_data).astype(int)
        figures['segment'] = visualize_discrete_label(
            seg_data, {1: 'protein'}, "Segment (Protein Mask)",
            colors={1: '#4C72B0'})

    # === Atom Type（离散） ===
    atom_path = os.path.join(entry_dir, "label_atom.mrc")
    if os.path.exists(atom_path):
        atom_data, _ = load_mrc(atom_path)
        atom_data = np.round(atom_data).astype(int)
        atom_colors = {1: '#e41a1c', 2: '#377eb8', 3: '#4daf4a', 4: '#984ea3', 5: '#ff7f00', 6: '#999999'}
        figures['atom'] = visualize_discrete_label(
            atom_data, ATOM_NAMES, "Atom Type", colors=atom_colors)

    # === Secondary Structure（离散） ===
    ss_path = os.path.join(entry_dir, "label_ss.mrc")
    if os.path.exists(ss_path):
        ss_data, _ = load_mrc(ss_path)
        ss_data = np.round(ss_data).astype(int)
        ss_colors = {1: '#F44336', 2: '#4CAF50', 3: '#2196F3'}
        figures['ss'] = visualize_discrete_label(
            ss_data, SS_NAMES, "Secondary Structure", colors=ss_colors)

    # === Chain（离散） ===
    chain_path = os.path.join(entry_dir, "label_chain.mrc")
    if os.path.exists(chain_path):
        chain_data, _ = load_mrc(chain_path)
        chain_data = np.round(chain_data).astype(int)
        n_chains = int(chain_data.max())
        chain_names = {i: f'Chain {i}' for i in range(1, n_chains + 1)}
        figures['chain'] = visualize_discrete_label(
            chain_data, chain_names, f"Chain ID ({n_chains} chains)")

    # === Domain（离散） ===
    dom_path = os.path.join(entry_dir, "label_domain.mrc")
    if os.path.exists(dom_path):
        dom_data, _ = load_mrc(dom_path)
        dom_data = np.round(dom_data).astype(int)
        if dom_data.max() > 0:
            n_doms = int(dom_data.max())
            dom_names = {i: f'Domain {i}' for i in range(1, n_doms + 1)}
            figures['domain'] = visualize_discrete_label(
                dom_data, dom_names, f"Domain ({n_doms} domains)")

    # === Interface（离散，binary） ===
    iface_path = os.path.join(entry_dir, "label_interface.mrc")
    if os.path.exists(iface_path):
        iface_data, _ = load_mrc(iface_path)
        iface_data = np.round(iface_data).astype(int)
        if iface_data.max() > 0:
            figures['interface'] = visualize_discrete_label(
                iface_data, {1: 'interface'}, "Protein-Protein Interface",
                colors={1: '#C44E52'})

    # === Q-score（连续） ===
    qs_path = os.path.join(entry_dir, "label_qscore.mrc")
    if os.path.exists(qs_path):
        qs_data, _ = load_mrc(qs_path)
        if (qs_data != 0).sum() > 0:
            figures['qscore'] = visualize_continuous_label(
                qs_data, "Q-score", threshold_pct=30, colorscale='RdYlGn')

    # === mol_map（连续） ===
    mol_path = os.path.join(entry_dir, "mol_map.mrc")
    if os.path.exists(mol_path):
        mol_data, _ = load_mrc(mol_path)
        figures['mol_map'] = visualize_continuous_label(
            mol_data, "Simulated Density (mol_map)", threshold_pct=50, colorscale='Inferno')

    return figures, entry_name


def build_combined_html(figures, entry_name):
    """
    将多个 plotly figure 组合成一个带选项卡的 HTML 页面

    使用懒渲染策略：只在标签页首次可见时调用 Plotly.newPlot，
    避免在 display:none 的 div 上渲染导致 0 尺寸问题。
    """
    html_parts = []
    html_parts.append(f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{entry_name} - 3D Label Visualization</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        h1 {{ color: #333; text-align: center; }}
        .tab-container {{ display: flex; flex-wrap: wrap; gap: 5px; justify-content: center; margin: 20px 0; }}
        .tab-btn {{ padding: 8px 16px; border: 1px solid #ddd; background: #fff; cursor: pointer;
                    border-radius: 4px; font-size: 14px; }}
        .tab-btn.active {{ background: #4C72B0; color: white; border-color: #4C72B0; }}
        .plot-container {{ display: none; width: 100%; height: 700px; background: white;
                          border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .plot-container.active {{ display: block; }}
        .info {{ text-align: center; color: #666; margin: 10px 0; font-size: 13px; }}
    </style>
</head>
<body>
    <h1>{entry_name} - 3D Label Visualization (V3.2)</h1>
    <p class="info">Click tabs to switch between label channels. Use mouse to rotate/zoom/pan.</p>
    <div class="tab-container">
""")

    tab_names = list(figures.keys())
    for i, name in enumerate(tab_names):
        active = ' active' if i == 0 else ''
        html_parts.append(f'        <button class="tab-btn{active}" onclick="showTab(\'{name}\')">{name}</button>\n')

    html_parts.append('    </div>\n')

    for i, (name, fig) in enumerate(figures.items()):
        active = ' active' if i == 0 else ''
        div_id = f'plot_{name}'
        html_parts.append(f'    <div id="{div_id}" class="plot-container{active}"></div>\n')

    # 懒渲染：存储 figure JSON，仅在标签页首次可见时渲染
    html_parts.append("""
    <script>
        var figureData = {};
        var rendered = {};

        function showTab(name) {
            document.querySelectorAll('.plot-container').forEach(function(el) { el.classList.remove('active'); });
            document.querySelectorAll('.tab-btn').forEach(function(el) { el.classList.remove('active'); });
            document.getElementById('plot_' + name).classList.add('active');
            event.target.classList.add('active');

            // 懒渲染：首次显示时创建 plotly 图
            if (!rendered[name]) {
                var data = figureData[name];
                Plotly.newPlot('plot_' + name, data.data, data.layout, {responsive: true});
                rendered[name] = true;
            } else {
                Plotly.Plots.resize(document.getElementById('plot_' + name));
            }
        }
""")

    # 嵌入每个 figure 的 JSON 数据（不立即渲染）
    for name, fig in figures.items():
        fig_json = fig.to_json()
        html_parts.append(f"""
        figureData['{name}'] = {fig_json};
""")

    # 立即渲染第一个标签页（它是可见的）
    if tab_names:
        first = tab_names[0]
        html_parts.append(f"""
        // 渲染首个可见标签页
        Plotly.newPlot('plot_{first}', figureData['{first}'].data, figureData['{first}'].layout, {{responsive: true}});
        rendered['{first}'] = true;
""")

    html_parts.append("""
    </script>
</body>
</html>""")

    return ''.join(html_parts)


def save_static_png(figures, output_path):
    """尝试导出静态 PNG（需要 kaleido）"""
    try:
        import kaleido  # noqa: F401
    except ImportError:
        print("  kaleido 未安装，跳过 PNG 导出 (pip install kaleido)")
        return False

    # 选一个最有代表性的图：segment 或第一个
    fig = figures.get('segment') or figures.get('density') or list(figures.values())[0]
    try:
        fig.write_image(output_path, width=1200, height=800, scale=2)
        print(f"  静态截图已保存: {output_path}")
        return True
    except Exception as e:
        print(f"  PNG 导出失败: {e}")
        return False


def process_entry(entry_dir, output_html=None, output_png=None):
    """处理单个条目的 3D 可视化"""
    print(f"\n3D 可视化: {entry_dir}")

    # 检查必要文件
    if not os.path.exists(os.path.join(entry_dir, "map_normalized.mrc")):
        print(f"  跳过: 无 map_normalized.mrc")
        return
    if not os.path.exists(os.path.join(entry_dir, "label_segment.mrc")):
        print(f"  跳过: 无 V3 标签文件")
        return

    figures, entry_name = create_combined_figure(entry_dir)

    if not figures:
        print(f"  无可视化数据")
        return

    print(f"  生成 {len(figures)} 个 3D 视图: {', '.join(figures.keys())}")

    # 生成 HTML
    if output_html is None:
        output_html = os.path.join(entry_dir, "vis_3d_labels.html")
    html_content = build_combined_html(figures, entry_name)
    with open(output_html, 'w') as f:
        f.write(html_content)
    print(f"  交互式 HTML: {output_html}")

    # 尝试导出 PNG
    if output_png is None:
        output_png = os.path.join(entry_dir, "vis_3d_labels.png")
    save_static_png(figures, output_png)


def main():
    parser = argparse.ArgumentParser(description='V3.1 3D 等值面可视化')
    parser.add_argument('entry_dir', nargs='?',
                        help='条目目录路径，不指定则处理所有')
    parser.add_argument('--data-dir', default='data/raw',
                        help='数据根目录 (默认: data/raw)')
    parser.add_argument('--output-html', default=None,
                        help='输出 HTML 路径（仅单条目模式）')
    parser.add_argument('--output-png', default=None,
                        help='输出 PNG 路径（仅单条目模式）')
    args = parser.parse_args()

    if args.entry_dir:
        dirs = [args.entry_dir]
    else:
        base = args.data_dir
        if not os.path.exists(base):
            print(f"数据目录不存在: {base}")
            sys.exit(1)
        dirs = [os.path.join(base, d) for d in sorted(os.listdir(base))
                if os.path.isdir(os.path.join(base, d))]

    for entry_dir in dirs:
        process_entry(
            entry_dir,
            output_html=args.output_html if len(dirs) == 1 else None,
            output_png=args.output_png if len(dirs) == 1 else None,
        )

    print(f"\n3D 可视化完成: {len(dirs)} 个条目")


if __name__ == '__main__':
    main()
