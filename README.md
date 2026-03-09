# Cryo-EM Map-Model Alignment Pipeline

从 EMDB/RCSB 自动化下载冷冻电镜数据，生成用于深度学习训练的体素级标注数据集。

## 核心特性

- **生物学组装体展开** — 通过 `gemmi.make_assembly()` 将不对称单元展开为完整生物学组装体（如 apoferritin 24聚体从 3k 原子展开到 76k 原子）
- **物理级模拟密度** — 使用 `gemmi.DensityCalculatorX`（5高斯原子散射因子），CC_mask 显著优于手写三线性+高斯模糊
- **双层标注策略** — 距离阈值内高置信度 + Voronoi 全覆盖低置信度，标注覆盖率从 <1% 提升到 25-35%
- **多分子类型标注** — 支持蛋白质、RNA、DNA、糖基、配体的自动检测与分类标注
- **Q-score 物理置信度** — Pintilie 2020 标准原子级密度拟合指标，替代启发式置信度
- **数据质量分层** — Gold/Silver/Hard 三级自动分类，基于分辨率、CC_mask、Q-score
- **自适应分辨率 blur** — mol_map 生成按条目分辨率自动缩放平滑度
- **Robust Z-score 归一化** — 基于中位数和 MAD，对异常值不敏感

## 管线流程

```
1. Retrieval     — RCSB Search API 查询 + EMDB/RCSB 下载（含分辨率 metadata）
2. Resample      — scipy.ndimage.zoom 重采样到统一体素大小 (1.0Å)
3. Normalization — Robust z-score 归一化到 [0, 1]
4. Alignment QC  — DensityCalculatorX 模拟密度 + CC 三项指标 + Q-score
5. Correspondence— 生物组装体展开 + 6层体素标注 + 多分子类型 + Voronoi
6. Enhancement   — DensityCalculatorX 自适应 blur 生成去噪训练标签 (mol_map)
7. Redundancy    — MMseqs2/3-mer 聚类 + Gold/Silver/Hard 分层 + train/val/test 划分
```

## 多分子类型标注

自动检测 mmCIF 文件中的分子类型并分类标注：

| 分子类型 | 标签值 | 检测方法 |
|----------|--------|----------|
| Protein | 1 | `gemmi.find_tabulated_residue().is_amino_acid()` |
| RNA | 2 | `is_nucleotide()` + 含 O2' 核糖 |
| DNA | 3 | `is_nucleotide()` + 脱氧核糖 |
| Sugar | 4 | 残基名匹配 (NAG, MAN, BMA, FUC, GAL, SIA 等) |
| Ligand | 5 | 其他非水小分子 |

**残基标签扩展：**
- 氨基酸 (1-20): ALA...VAL
- 核苷酸 (21-28): A, U, G, C, DA, DT, DG, DC
- 糖基 (29-34): NAG, MAN, BMA, FUC, GAL, SIA

**二级结构扩展：**
- 蛋白质: coil=1, helix=2, strand=3
- RNA/DNA: sugar_phosphate=4, base=5

## Q-score 原理

Q-score (Pintilie et al., Nature Methods 2020) 是 cryo-EM 领域标准的原子级密度拟合质量指标：

1. 对每个原子，沿 40 个径向方向（Fibonacci sphere）采样实验密度
2. 采样范围 0~2Å，步长 0.1Å
3. 将径向平均密度 profile 与参考高斯 `exp(-d²/(2σ²))` (σ=0.6Å) 做 Pearson 相关
4. 返回 [-1, 1]，越接近 1 表示局部密度与原子位置匹配越好

## 数据质量分层

| 层级 | 条件 | 适用场景 |
|------|------|---------|
| Gold | resolution < 2.5Å AND cc_mask > 0.7 AND q_score > 0.5 | 高质量训练数据 |
| Silver | resolution < 3.5Å AND cc_mask > 0.5 | 一般训练数据 |
| Hard | 其余通过 QC 的条目 | 困难样本/增广数据 |

## 快速开始

### 环境

```bash
conda create -n cryo-pipeline python=3.10
conda activate cryo-pipeline
pip install gemmi mrcfile biopython scipy numpy requests tqdm pyyaml matplotlib
```

### 运行测试

```bash
# 下载 4 个测试条目（蛋白+RNA+糖蛋白）并跑通全流程
python test_pipeline.py

# 生成数据集概览图
python visualize_dataset.py

# 生成标签可视化
python visualize_labels.py
```

### 运行完整管线

```bash
# 使用默认配置（下载 5 个条目）
python run_pipeline.py

# 指定配置和步骤
python run_pipeline.py --config configs/default.yaml --steps all

# 只跑后续步骤（已下载数据）
python run_pipeline.py --steps normalization,alignment_qc,correspondence,enhancement
```

## 输出结构

```
data/raw/EMD-XXXXX_YYYY/
├── raw_map.map           # 原始密度图
├── map_std_1.0A.mrc      # 重采样后
├── map_normalized.mrc    # 归一化后
├── sim_map.mrc           # 模拟密度图 (QC)
├── mol_map.mrc           # 模拟密度图 (去噪标签, 自适应blur)
├── model.cif             # 原子模型
├── label_atom.mrc        # 原子类型 (CA/N/C/O/CB/other)
├── label_aa.mrc          # 残基类型 (氨基酸1-20/核苷酸21-28/糖基29-34)
├── label_ss.mrc          # 二级结构 (coil/helix/strand/sugar_phosphate/base)
├── label_chain.mrc       # 链 ID
├── label_confidence.mrc  # 置信度 (Hard高置信 + Voronoi低置信)
├── label_moltype.mrc     # 分子类型 (protein/RNA/DNA/sugar/ligand)
├── qc_metrics.json       # CC_mask, CC_volume, CC_overall, Q-score
├── labeling_stats.json   # 覆盖率统计 + 分子类型计数
└── metadata.json         # PDB/EMDB ID, 分辨率
```

## 配置说明

核心配置项 (`configs/default.yaml`):

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `bio_assembly.enabled` | true | 启用生物组装体展开 |
| `bio_assembly.max_chains` | 200 | 最大链数安全限制 |
| `normalization.method` | robust_zscore | 归一化方法 |
| `alignment_qc.dc_blur` | 0 | DensityCalculatorX blur (0最优) |
| `qscore.sigma` | 0.6 | Q-score 参考高斯宽度 (Å) |
| `qscore.n_directions` | 40 | Q-score 采样方向数 |
| `correspondence.use_voronoi` | true | Voronoi 全覆盖标注 |
| `correspondence.molecule_types.enabled` | true | 多分子类型标注 |
| `enhancement.adaptive_blur` | true | 按分辨率自适应 blur |
| `redundancy.use_mmseqs2` | true | MMseqs2 聚类（自动回退） |
| `redundancy.tiers.gold.max_resolution` | 2.5 | Gold 层最大分辨率 |

## 依赖

- Python 3.10+
- gemmi >= 0.7.0
- mrcfile, biopython, scipy, numpy, requests, tqdm, pyyaml
- matplotlib (可视化)
- MMseqs2 (可选，自动下载静态二进制)
