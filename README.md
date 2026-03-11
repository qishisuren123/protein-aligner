# Cryo-EM Map-Model Alignment Pipeline V3

面向 Cryo-EM Model Building 的多层级标注数据集生成管线。

## V3 设计理念

V3 针对 model building 任务进行了根本性重构：

- **CA-only 标注范式**: 语义标签（氨基酸、二级结构、链、结构域、界面）只在 CA 原子位置标注，减少噪声，聚焦残基级预测
- **三组 KDTree 架构**: all-atom (segment/atom) + CA-only (aa/ss/chain/domain/interface) + backbone (qscore)
- **蛋白质专注**: 移除 RNA/DNA/糖基/配体标注，专注 model building 核心需求
- **物理级 Q-score**: 替代启发式高斯衰减置信度，使用 Pintilie 2020 原子级密度拟合指标
- **Domain 标签**: Merizo 深度学习工具进行结构域分割
- **Interface 标签**: CA-CA 跨链距离检测蛋白-蛋白界面
- **Pfam Fold/Family Split**: HMMER+Pfam 做更严格的数据划分防止结构同源泄露

## V3 标签通道 (8 + 1)

| # | 文件 | 内容 | 标注范围 | 值域 |
|---|------|------|---------|------|
| 1 | `label_segment.mrc` | 蛋白/背景 | 所有信号体素 | 0=背景, 1=蛋白 |
| 2 | `label_atom.mrc` | 原子类型 | 所有原子位置 | CA=1,C=2,N=3,O=4,CB=5,Other=6 |
| 3 | `label_aa.mrc` | 氨基酸类型 | CA-only | 1-20 (标准 AA) |
| 4 | `label_ss.mrc` | 二级结构 | CA-only | Helix=1,Strand=2,Coil=3 |
| 5 | `label_qscore.mrc` | Q-score 置信度 | 骨架原子 | float [-1, 1] |
| 6 | `label_chain.mrc` | 链 ID | CA-only | 1-N |
| 7 | `label_domain.mrc` | 结构域 ID | CA-only | 1-M (Merizo) |
| 8 | `label_interface.mrc` | 界面标签 | CA-only | 0=非界面, 1=界面 |
| 9 | `mol_map.mrc` | 去噪模拟密度 | 所有体素 | float [0, 1] |

### 相比 V2 的变化

| V2 | V3 | 变化原因 |
|----|----|----|
| `label_confidence.mrc` (高斯衰减) | `label_qscore.mrc` (物理 Q-score) | 物理级原子置信度替代启发式 |
| `label_moltype.mrc` (5类分子类型) | 删除 | 蛋白专注，无需多分子类型 |
| 无 | `label_segment.mrc` | 新增蛋白/背景二值分割标签 |
| 无 | `label_domain.mrc` | 新增结构域分割 |
| 无 | `label_interface.mrc` | 新增蛋白-蛋白界面检测 |
| 所有原子同一 KDTree | 三组 KDTree | CA-only 减少噪声，语义标签更准确 |
| Voronoi 默认开启 | Voronoi 默认关闭 | CA-only 标签无需扩展覆盖 |

## 管线流程 (8 步)

```
1. Retrieval        — RCSB/EMDB 下载
2. Resample         — 重采样到 1.0Å
3. Normalization    — Robust z-score
4. Alignment QC     — CC 指标 + Q-score 持久化
5. Domain Seg.      — Merizo 结构域分割（新增）
6. Correspondence   — 三 KDTree + CA-only + 8 层体素标注
7. Enhancement      — 自适应 blur mol_map
8. Redundancy       — MMseqs2 + Pfam Fold/Family Split
```

## 核心技术

### 三组 KDTree 架构

| KDTree | 源原子 | 产出标签 |
|--------|--------|---------|
| `tree_all` | 所有蛋白原子 | `label_segment`(距离内=1), `label_atom`(原子类型) |
| `tree_ca` | 仅 CA 原子 | `label_aa`, `label_ss`, `label_chain`, `label_domain`, `label_interface` |
| `tree_backbone` | CA/C/N/O/CB | `label_qscore`(从 qscores.json 读取 Q-score 值) |

### Q-score (Pintilie et al., Nature Methods 2020)

- 对每个原子，沿 40 个径向方向 (Fibonacci sphere) 采样实验密度
- 采样范围 0~2Å，步长 0.1Å
- 与参考高斯 `exp(-d²/(2σ²))` (σ=0.6Å) 做 Pearson 相关
- 返回 [-1, 1]，越接近 1 匹配越好
- Step 4 保存逐原子 Q-score 到 `qscores.json`，Step 6 读取生成 `label_qscore.mrc`

### Domain Segmentation (Merizo)

- 基于深度学习的蛋白质结构域分割
- 在原始结构上运行，通过残基序号复制到展开后的对称链
- 回退策略：Merizo 不可用时 `label_domain.mrc` 全零，管线不中断

### Interface Detection

- 基于 CA-CA 跨链距离 (阈值 5.0Å)
- 每条链建 cKDTree，批量查询其他链的最近距离
- 高效处理多链蛋白（如 apoferritin 24 聚体）

### Pfam Fold/Family Split

- HMMER `hmmscan` 搜索 Pfam-A 数据库
- 共享 Pfam family 的条目归入同一组
- 按 Pfam 分组做 train/val/test 划分，防止结构同源泄露
- 回退：Pfam 不可用时自动使用 MMseqs2 序列聚类

## 快速开始

### 环境

```bash
conda create -n cryo-pipeline python=3.10
conda activate cryo-pipeline
pip install gemmi mrcfile biopython scipy numpy requests tqdm pyyaml matplotlib
```

### 安装外部工具（可选）

```bash
# 安装 Merizo + HMMER + Pfam-A（需要代理时先设置）
bash scripts/install_tools.sh
```

### 运行测试

```bash
# 下载 4 个测试条目并跑通全流程
python test_pipeline.py

# 生成数据集概览图
python visualize_dataset.py

# 生成标签可视化
python visualize_labels.py
```

### 运行完整管线

```bash
python run_pipeline.py                          # 默认配置，下载 5 个条目
python run_pipeline.py --steps all              # 全部 8 步
python run_pipeline.py --steps domain_segmentation,correspondence,enhancement  # 只跑部分步骤
```

## 输出结构

每个条目产出 18 个文件：

```
data/raw/EMD-XXXXX_YYYY/
├── raw_map.map           # 原始密度图
├── map_std_1.0A.mrc      # 重采样到 1.0Å
├── map_normalized.mrc    # Robust z-score 归一化
├── sim_map.mrc           # QC 用模拟密度图
├── mol_map.mrc           # 去噪训练标签（自适应 blur）
├── model.cif             # 原子模型
├── label_segment.mrc     # 蛋白/背景 (0/1)
├── label_atom.mrc        # 原子类型 (CA=1, C=2, N=3, O=4, CB=5, Other=6)
├── label_aa.mrc          # 氨基酸类型 (1-20)，CA-only
├── label_ss.mrc          # 二级结构 (Helix=1, Strand=2, Coil=3)，CA-only
├── label_qscore.mrc      # Q-score 置信度 (float)，骨架原子
├── label_chain.mrc       # 链 ID (1-N)，CA-only
├── label_domain.mrc      # 结构域 ID (1-M)，CA-only
├── label_interface.mrc   # 界面标签 (0/1)，CA-only
├── qscores.json          # 逐原子 Q-score
├── domain_assignment.json # 结构域分割结果
├── qc_metrics.json       # CC + Q-score 指标
├── labeling_stats.json   # V3 标注统计
└── metadata.json         # PDB/EMDB ID, 分辨率
```

打包输出：
```
data/output/
├── train/ val/ test/       # 按 Pfam/序列聚类划分（防同源泄露）
├── dataset_report.json     # 全量报告（含 tier 分层 + Pfam 注释）
└── dataset_overview.png    # V3 数据集概览
```

## 配置

核心配置项 (`configs/default.yaml`):

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `domain_segmentation.enabled` | true | Merizo 结构域分割 |
| `correspondence.use_voronoi` | false | V3 默认关闭 Voronoi |
| `correspondence.interface.distance_threshold` | 5.0 | CA-CA 跨链界面阈值 (Å) |
| `redundancy.pfam.enabled` | true | Pfam Fold/Family Split |

## 数据质量分层

| 层级 | 条件 | 用途 |
|------|------|------|
| Gold | resolution < 2.5Å AND cc_mask > 0.7 AND q_score > 0.5 | 高质量训练集 |
| Silver | resolution < 3.5Å AND cc_mask > 0.5 | 一般训练集 |
| Hard | 其余通过 QC 的条目 | 困难样本 |

## 风险与回退

| 风险 | 回退方案 |
|------|---------|
| Merizo 安装失败 | `label_domain.mrc` 全零，管线继续 |
| Pfam-A 下载失败 | 回退到 MMseqs2 聚类 split |
| PyTorch CPU 安装失败 | Merizo 不可用，domain 全零 |
| hmmscan 不可用 | 回退到 MMseqs2 |

## 依赖

- Python 3.10+
- gemmi >= 0.7.0 — 结构解析、密度计算
- scipy — KDTree、重采样
- numpy — 向量化计算
- requests, tqdm — 数据下载
- matplotlib — 可视化
- MMseqs2 — 序列聚类（可选，自动下载）
- PyTorch CPU — Merizo 依赖（可选）
- HMMER — Pfam 搜索（可选）
