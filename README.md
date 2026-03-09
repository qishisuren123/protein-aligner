# Cryo-EM Map-Model Alignment Pipeline

从 EMDB/RCSB 自动化下载冷冻电镜数据，生成用于深度学习训练的体素级标注数据集。

## 核心特性

- **生物学组装体展开** — 通过 `gemmi.make_assembly()` 将不对称单元展开为完整生物学组装体（如 apoferritin 24聚体从 3k 原子展开到 76k 原子）
- **物理级模拟密度** — 使用 `gemmi.DensityCalculatorX`（5高斯原子散射因子），CC_mask 显著优于手写三线性+高斯模糊
- **双层标注策略** — 距离阈值内高置信度 + Voronoi 全覆盖低置信度，标注覆盖率从 <1% 提升到 25-35%
- **Robust Z-score 归一化** — 基于中位数和 MAD，对异常值不敏感

## 管线流程

```
1. Retrieval     — RCSB Search API 查询 + EMDB/RCSB 下载（含分辨率 metadata）
2. Resample      — scipy.ndimage.zoom 重采样到统一体素大小 (1.0Å)
3. Normalization — Robust z-score 归一化到 [0, 1]
4. Alignment QC  — DensityCalculatorX 模拟密度 + CC_mask/CC_volume/CC_overall
5. Correspondence— 生物组装体展开 + 5层体素标注 + Voronoi 全覆盖
6. Enhancement   — DensityCalculatorX 生成去噪训练标签 (mol_map)
7. Redundancy    — 序列聚类 + train/val/test 划分
```

## 快速开始

### 环境

```bash
conda create -n cryo-pipeline python=3.10
conda activate cryo-pipeline
pip install gemmi mrcfile biopython scipy numpy requests tqdm pyyaml matplotlib
```

### 运行测试

```bash
# 下载 2 个测试条目并跑通全流程
python test_pipeline.py

# 生成可视化
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
├── mol_map.mrc           # 模拟密度图 (去噪标签)
├── model.cif             # 原子模型
├── label_atom.mrc        # 原子类型 (CA/N/C/O/CB/other)
├── label_aa.mrc          # 氨基酸类型 (20类)
├── label_ss.mrc          # 二级结构 (coil/helix/strand)
├── label_chain.mrc       # 链 ID
├── label_confidence.mrc  # 置信度 (Hard高置信 + Voronoi低置信)
├── qc_metrics.json       # CC_mask, CC_volume, CC_overall
├── labeling_stats.json   # 覆盖率统计
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
| `correspondence.use_voronoi` | true | Voronoi 全覆盖标注 |
| `correspondence.voronoi_max_confidence` | 0.3 | Voronoi 最大置信度 |

## 依赖

- Python 3.10+
- gemmi >= 0.7.0
- mrcfile, biopython, scipy, numpy, requests, tqdm, pyyaml
- matplotlib (可视化)
