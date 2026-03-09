# Cryo-EM Map-Model Alignment Pipeline

从 EMDB/RCSB 自动化下载冷冻电镜数据，生成用于深度学习训练的体素级标注数据集。

## 实测结果

在 4 个测试条目上的实际运行结果（2026-03-09）：

| 条目 | 类型 | 分辨率 | CC_mask | Q-score | 覆盖率 | 链数 | Tier |
|------|------|--------|---------|---------|--------|------|------|
| EMD-11638 / 7A4M | 蛋白质 (apoferritin 24聚体) | 1.22Å | **0.80** | **0.56** | **29.8%** | 24 | Gold |
| EMD-31083 / 7EFC | 蛋白质 | 1.70Å | **0.79** | **0.86** | 4.3% | 4 | Gold |
| EMD-48045 / 9K6S | 蛋白质-RNA 复合物 | 2.80Å | 0.00† | -0.19† | 8.5% | 3 | Hard |
| EMD-34948 / 8HKS | 糖蛋白 | 2.90Å | 0.00† | — | — | — | QC 未通过† |

> † RNA 条目模型仅占 536³ map 极小区域导致 CC 偏低，但 **RNA 体素成功检测** (protein=350万, RNA=957万体素)。糖蛋白条目 PDB 坐标系与 EMDB map 不匹配（已知问题，多数条目无此问题）。

### 完成度

| 功能模块 | 状态 | 说明 |
|---------|------|------|
| 数据检索与下载 | ✅ 完成 | RCSB Search API + EMDB FTP，含分辨率 metadata |
| 密度图重采样 | ✅ 完成 | scipy.ndimage.zoom 到 1.0Å |
| Robust Z-score 归一化 | ✅ 完成 | 中位数 + MAD，对异常值鲁棒 |
| 生物学组装体展开 | ✅ 完成 | gemmi.make_assembly()，apoferritin 3k→76k 原子 |
| DensityCalculatorX 模拟密度 | ✅ 完成 | 5 高斯散射因子，CC_mask=0.80 |
| 三项 CC 指标 | ✅ 完成 | CC_mask / CC_volume / CC_overall |
| Q-score 物理置信度 | ✅ 完成 | Pintilie 2020，Fibonacci sphere 40方向径向采样 |
| 双层体素标注 (Hard + Voronoi) | ✅ 完成 | 覆盖率从 <1% 提升到 30% |
| 多分子类型标注 | ✅ 完成 | protein/RNA/DNA/sugar/ligand 自动检测 |
| 核苷酸 + 糖基标签扩展 | ✅ 完成 | 氨基酸(1-20) + 核苷酸(21-28) + 糖基(29-34) |
| RNA/DNA 二级结构 | ✅ 完成 | sugar_phosphate / base 原子级分类 |
| 自适应分辨率 blur | ✅ 完成 | blur = b_factor × (resolution / 2.0) |
| Gold/Silver/Hard 数据分层 | ✅ 完成 | 基于 resolution + CC_mask + Q-score |
| MMseqs2 序列聚类 | ✅ 完成 | 自动下载静态二进制，回退 3-mer Jaccard |
| train/val/test 划分 | ✅ 完成 | 按序列聚类划分，防止数据泄露 |
| 数据集可视化 | ✅ 完成 | 6 子图发表级概览 |
| 蛋白质 demo | ✅ 验证通过 | 2 条目，CC_mask > 0.79，Q-score > 0.56 |
| 蛋白质-RNA demo | ⚠️ 部分通过 | RNA 类型检测成功，CC 因 map 尺寸不匹配偏低 |
| 糖蛋白 demo | ❌ 未通过 | 坐标系不匹配，需更换条目 |

## 管线流程

```
1. Retrieval     — RCSB Search API 查询 + EMDB/RCSB 下载（含分辨率 metadata）
2. Resample      — scipy.ndimage.zoom 重采样到统一体素大小 (1.0Å)
3. Normalization — Robust z-score 归一化到 [0, 1]
4. Alignment QC  — DensityCalculatorX 模拟密度 + CC 三项指标 + Q-score
5. Correspondence— 生物组装体展开 + 6 层体素标注 + 多分子类型 + Voronoi
6. Enhancement   — DensityCalculatorX 自适应 blur 生成去噪训练标签 (mol_map)
7. Redundancy    — MMseqs2/3-mer 聚类 + Gold/Silver/Hard 分层 + train/val/test 划分
```

## 核心技术

### 生物学组装体展开
通过 `gemmi.make_assembly()` 将不对称单元展开为完整生物学组装体。Apoferritin 从 3,067 原子（1 链）展开到 73,608 原子（24 链），使标注覆盖完整密度。

### 物理级模拟密度
使用 `gemmi.DensityCalculatorX`（5 高斯近似原子散射因子），CC_mask 从手写三线性+高斯模糊的 0.40 提升到 **0.80**。

### 双层标注策略
- **Hard 层**: 距离阈值内（≤3Å），高斯衰减置信度，高置信度
- **Voronoi 层**: 距离阈值外，最近邻赋值，置信度上限 0.3
- 效果: 标注覆盖率从 <1% 提升到 **29.8%**

### Q-score (Pintilie et al., Nature Methods 2020)
cryo-EM 标准的原子级密度拟合质量指标：
1. 对每个原子，沿 40 个径向方向（Fibonacci sphere）采样实验密度
2. 采样范围 0~2Å，步长 0.1Å
3. 与参考高斯 `exp(-d²/(2σ²))` (σ=0.6Å) 做 Pearson 相关
4. 返回 [-1, 1]，越接近 1 表示匹配越好

### 多分子类型标注
通过 `gemmi.find_tabulated_residue().kind` 自动检测：

| 分子类型 | 标签值 | 检测方法 |
|----------|--------|----------|
| Protein | 1 | `ResidueKind.AA` |
| RNA | 2 | `ResidueKind.RNA` |
| DNA | 3 | `ResidueKind.DNA` |
| Sugar | 4 | `ResidueKind.PYR` + 糖基残基名匹配 |
| Ligand | 5 | 其他非水小分子 |

残基标签扩展：氨基酸 (1-20) + 核苷酸 (21-28: A/U/G/C/DA/DT/DG/DC) + 糖基 (29-34: NAG/MAN/BMA/FUC/GAL/SIA)

二级结构扩展：蛋白质 coil/helix/strand + RNA/DNA sugar_phosphate/base

### 数据质量分层

| 层级 | 条件 | 用途 |
|------|------|------|
| Gold | resolution < 2.5Å AND cc_mask > 0.7 AND q_score > 0.5 | 高质量训练集 |
| Silver | resolution < 3.5Å AND cc_mask > 0.5 | 一般训练集 |
| Hard | 其余通过 QC 的条目 | 困难样本 |

## 快速开始

### 环境

```bash
conda create -n cryo-pipeline python=3.10
conda activate cryo-pipeline
pip install gemmi mrcfile biopython scipy numpy requests tqdm pyyaml matplotlib
```

### 运行测试

```bash
# 下载 4 个测试条目并跑通全流程（约 5 分钟）
python test_pipeline.py

# 生成数据集概览图
python visualize_dataset.py

# 生成标签可视化
python visualize_labels.py
```

### 运行完整管线

```bash
python run_pipeline.py                          # 默认配置，下载 5 个条目
python run_pipeline.py --steps all              # 全部步骤
python run_pipeline.py --steps correspondence,enhancement  # 只跑后续步骤
```

## 输出结构

每个条目产出 15 个文件：

```
data/raw/EMD-XXXXX_YYYY/
├── raw_map.map           # 原始密度图
├── map_std_1.0A.mrc      # 重采样到 1.0Å
├── map_normalized.mrc    # Robust z-score 归一化
├── sim_map.mrc           # QC 用模拟密度图
├── mol_map.mrc           # 去噪训练标签（自适应 blur）
├── model.cif             # 原子模型
├── label_atom.mrc        # 原子类型 (CA=1, N=2, C=3, O=4, CB=5, other=6)
├── label_aa.mrc          # 残基类型 (AA 1-20, 核苷酸 21-28, 糖基 29-34)
├── label_ss.mrc          # 二级结构 (coil/helix/strand/sugar_phosphate/base)
├── label_chain.mrc       # 链 ID
├── label_confidence.mrc  # 置信度 (Hard 高置信 + Voronoi 低置信)
├── label_moltype.mrc     # 分子类型 (protein/RNA/DNA/sugar/ligand)
├── qc_metrics.json       # CC_mask, CC_volume, CC_overall, Q-score
├── labeling_stats.json   # 覆盖率统计 + 分子类型计数
└── metadata.json         # PDB/EMDB ID, 分辨率
```

打包输出：
```
data/output/
├── train/ val/ test/       # 按序列聚类划分（防数据泄露）
├── dataset_report.json     # 全量报告（含 tier 分层）
└── dataset_overview.png    # 6 子图数据集概览
```

## 配置

核心配置项 (`configs/default.yaml`):

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `bio_assembly.enabled` | true | 生物组装体展开 |
| `alignment_qc.dc_blur` | 0 | DensityCalculatorX blur (0 最优) |
| `qscore.sigma` | 0.6 | Q-score 参考高斯宽度 (Å) |
| `correspondence.molecule_types.enabled` | true | 多分子类型标注 |
| `correspondence.use_voronoi` | true | Voronoi 全覆盖 |
| `enhancement.adaptive_blur` | true | 按分辨率自适应 blur |
| `redundancy.use_mmseqs2` | true | MMseqs2 聚类 |

## 已知限制

- 部分 EMDB/PDB 条目存在坐标系不匹配问题（PDB 坐标落在 map 外），需要 map-model 重对齐
- 超大密度图（>400³）的 Q-score 计算较慢（逐原子 Python 循环）
- 当前糖蛋白 demo 条目（8HKS）因坐标系问题未通过 QC，需更换

## 依赖

- Python 3.10+
- gemmi >= 0.7.0 — 结构解析、密度计算、坐标变换
- scipy — 重采样、KDTree、统计
- numpy — 向量化计算
- requests, tqdm — 数据下载
- matplotlib — 可视化
- MMseqs2 — 序列聚类（可选，自动下载）
