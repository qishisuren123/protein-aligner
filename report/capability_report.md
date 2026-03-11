# Cryo-EM Model Building 标注数据集管线 — 能力汇报

> 项目路径：`/mnt/shared-storage-user/renyiming/protein-aligner`
> 迭代历程：V1 (2026-03-09) → V2 (2026-03-09) → V3 (2026-03-11)
> 环境：conda `cryo-pipeline` (Python 3.10, gemmi 0.7.5)，CPU only

---

## 一、项目定位

本管线从 EMDB/RCSB 公共数据库自动化下载冷冻电镜密度图与原子模型，经过标准化预处理后，生成面向 **cryo-EM model building** 深度学习训练的 **体素级多层级标注数据集**。最终输出可直接用于训练蛋白质主链追踪、氨基酸类型预测、二级结构识别、结构域分割、界面检测等下游模型。

---

## 二、迭代历程

### V1：全流程搭建（2026-03-09）

从零搭建 7 步管线，跑通端到端流程：

- RCSB Search API 查询 + EMDB/RCSB 自动下载
- scipy.ndimage.zoom 重采样到统一体素大小（1.0Å）
- Robust z-score 归一化（中位数 + MAD，对异常值鲁棒）
- 手写三线性插值 + 高斯模糊模拟密度，CC 评估对齐质量
- 单 KDTree 最近邻体素标注（原子类型、氨基酸、二级结构、链 ID、高斯衰减置信度）
- 3-mer Jaccard 序列聚类 + train/val/test 划分

**问题**：apoferritin 只标注 1/24 密度（未展开组装体），CC_mask 仅 0.40（手写模拟密度粗糙），标注覆盖率 <1%。

### V1.5：对齐质量提升（2026-03-09）

针对 V1 核心痛点做了 5 项关键改进：

| 改进 | 效果 |
|------|------|
| `gemmi.make_assembly()` 生物组装体展开 | apoferritin: 3,067 → 73,608 原子，1 → 24 链 |
| `DensityCalculatorX` 物理级模拟密度 | CC_mask 从 0.40 → **0.80**（5 高斯原子散射因子）|
| Voronoi 双层标注（Hard + 低置信度扩展）| 标注覆盖率从 <1% → **29.8%** |
| Robust z-score 归一化 | 保留密度动态范围，消除异常值影响 |
| RCSB API 分辨率 metadata | 每条条目记录实际分辨率，用于后续分层 |

### V2：发表级质量（2026-03-09）

在对齐质量已达标的基础上，丰富标注维度和数据管理能力：

| 新增能力 | 说明 |
|---------|------|
| Q-score 物理置信度 | Pintilie 2020 标准，Fibonacci sphere 40 方向径向采样 |
| 多分子类型标注 | protein/RNA/DNA/sugar/ligand 自动检测 |
| 核苷酸 + 糖基标签扩展 | AA(1-20) + 核苷酸(21-28) + 糖基(29-34) |
| RNA/DNA 二级结构 | sugar_phosphate / base 原子级分类 |
| 自适应分辨率 blur | `blur = b_factor × (resolution / 2.0)` |
| Gold/Silver/Hard 数据分层 | 基于 resolution + CC_mask + Q-score |
| MMseqs2 序列聚类 | 自动下载静态二进制，回退 3-mer Jaccard |
| 数据集概览可视化 | 6 子图发表级概览 |

测试覆盖扩展至 4 条目：蛋白质 × 2、蛋白-RNA 复合物、糖蛋白。

### V3：Model Building 专注重构（2026-03-11）

针对 model building 任务需求做根本性架构重构：

| 核心变化 | 动机 |
|---------|------|
| CA-only 标注范式 | 语义标签聚焦残基级预测，减少原子级噪声 |
| 三组 KDTree 架构 | 不同粒度标签使用不同原子子集，更精准 |
| 蛋白专注 | 移除 RNA/DNA/糖基/配体，聚焦 model building 核心 |
| Q-score 替代高斯衰减 | 物理级原子置信度替代启发式 |
| Domain 标签 (Merizo) | 结构域分割，支持多 domain 蛋白分析 |
| Interface 标签 | 蛋白-蛋白界面检测，支持复合物分析 |
| Pfam Fold/Family Split | 更严格防泄露，替代纯序列聚类 |

---

## 三、当前完整能力（V3）

### 3.1 管线流程（8 步）

```
Step 1  Retrieval           RCSB Search API 查询 + EMDB/RCSB 自动下载
Step 2  Resample            scipy.ndimage.zoom 重采样到 1.0Å 统一体素
Step 3  Normalization       Robust z-score 归一化到 [0, 1]
Step 4  Alignment QC        DensityCalculatorX 模拟密度 + 三项 CC + Q-score 持久化
Step 5  Domain Seg.         Merizo 深度学习结构域分割（V3 新增）
Step 6  Correspondence      三 KDTree + CA-only + 8 层体素标注
Step 7  Enhancement         自适应 blur 生成去噪训练标签 (mol_map)
Step 8  Redundancy          Pfam Fold/Family Split + MMseqs2 聚类 + 分层 + 打包
```

### 3.2 标签通道（8 + 1）

| # | 文件 | 内容 | 标注范围 | 值域 | 下游用途 |
|---|------|------|---------|------|---------|
| 1 | `label_segment.mrc` | 蛋白/背景 二值分割 | 所有信号体素 | {0, 1} | 蛋白区域检测 |
| 2 | `label_atom.mrc` | 原子类型 | 所有原子位置 | CA=1,C=2,N=3,O=4,CB=5,Other=6 | 原子级定位 |
| 3 | `label_aa.mrc` | 20 种氨基酸类型 | CA-only | 1-20 | 残基类型预测 |
| 4 | `label_ss.mrc` | 二级结构 | CA-only | Helix=1,Strand=2,Coil=3 | 骨架追踪 |
| 5 | `label_qscore.mrc` | Q-score 原子置信度 | 骨架原子 | float [-1, 1] | 置信度加权 |
| 6 | `label_chain.mrc` | 链 ID | CA-only | 1-N | 多链分割 |
| 7 | `label_domain.mrc` | 结构域 ID | CA-only | 1-M | 结构域分割 |
| 8 | `label_interface.mrc` | 蛋白-蛋白界面 | CA-only | {0, 1} | 界面预测 |
| 9 | `mol_map.mrc` | 去噪模拟密度 | 所有体素 | float [0, 1] | 密度去噪训练 |

### 3.3 核心技术能力

#### 数据获取
- RCSB Search API v2 查询，支持分辨率 / 实验方法筛选
- EMDB FTP 密度图 + RCSB mmCIF 原子模型自动下载
- 分辨率 metadata 自动回填

#### 密度图预处理
- **重采样**：scipy.ndimage.zoom，支持任意目标体素大小（默认 1.0Å），cubic 插值
- **归一化**：Robust z-score（中位数 + MAD），对异常值鲁棒；备选 zscore / percentile
- **坐标系**：gemmi CCP4 标准读写，保证原子坐标与体素精确对齐

#### 对齐与质量控制
- **DensityCalculatorX**：5 高斯近似原子散射因子，物理级模拟密度
- **三项 CC 指标**：CC_mask（模拟 mask 内）、CC_volume（实验 mask 内）、CC_overall（全体素）
- **Q-score**：Pintilie 2020 标准，Fibonacci sphere 40 方向径向采样，逐原子保存到 JSON
- **生物组装体展开**：`gemmi.make_assembly()` 支持对称蛋白（apoferritin 24 聚体等）
- **QC 判定**：quality_score ≥ 0.5 OR CC_mask ≥ 0.3 OR CC_volume ≥ 0.3

#### 体素标注（V3 核心）
- **三组 KDTree**：all-atom / CA-only / backbone 三种粒度独立查询
- **CA-only 范式**：语义标签仅在 CA 原子位置标注，减少噪声
- **向量化赋值**：scipy.cKDTree 批量最近邻查询，高效处理大体素网格
- **Voronoi 扩展**（可选）：距离阈值外低置信度赋值，可配置开关

#### 结构域分割
- **Merizo** 深度学习工具封装，支持 CPU 推理
- mmCIF → PDB 自动转换
- 在原始结构上运行，残基序号复制到展开后的对称链
- 回退策略：不可用时 `label_domain.mrc` 全零，管线不中断

#### 界面检测
- CA-CA 跨链距离检测（默认阈值 5.0Å）
- 每条链建 cKDTree，批量查询所有其他链
- 高效处理大型多链复合物

#### 去噪标签
- DensityCalculatorX 自适应 blur 生成模拟密度
- `blur = b_factor × (resolution / 2.0)`，高分辨率保留细节，低分辨率增加平滑

#### 去冗余与数据划分
- **Pfam Fold/Family Split**（V3 首选）：HMMER hmmscan 搜索 Pfam-A，Union-Find 按 family 分组
- **MMseqs2 精确聚类**（回退方案）：自动下载静态二进制，identity=0.3
- **3-mer Jaccard**（兜底方案）：无外部依赖
- **Gold/Silver/Hard 数据分层**：基于 resolution + CC_mask + Q-score
- **train/val/test 按聚类划分**：同一聚类/family 不跨集，防止数据泄露
- **符号链接打包**：节省磁盘空间

#### 可视化
- **标签面板**：9 列正交切面（Density / Segment / Atom / AA / SS / Q-score / Chain / Domain / Interface）
- **统计图**：原子类型 / 氨基酸 / 二级结构 / Q-score / 覆盖率 / 链 / Domain / Interface 共 8 子图
- **数据集概览**：分辨率分布 / CC vs Resolution / Domain 数量 / 覆盖率 / Q-score / Interface 比例

### 3.4 每条目产出文件（18 个）

```
密度图（5）: raw_map.map / map_std_1.0A.mrc / map_normalized.mrc / sim_map.mrc / mol_map.mrc
标签（8）:   label_segment / label_atom / label_aa / label_ss / label_qscore / label_chain / label_domain / label_interface
元数据（5）: metadata.json / qc_metrics.json / labeling_stats.json / qscores.json / domain_assignment.json
```

### 3.5 回退与鲁棒性

管线设计了多层回退机制，确保在外部工具不可用时仍可正常运行：

| 组件 | 不可用时的行为 |
|------|--------------|
| Merizo (结构域分割) | `label_domain.mrc` 全零，管线继续 |
| HMMER + Pfam-A | 回退到 MMseqs2 聚类 split |
| MMseqs2 | 回退到 3-mer Jaccard 近似聚类 |
| PyTorch CPU | Merizo 不可用，domain 标签全零 |
| 生物组装体信息缺失 | 使用原始不对称单元 |
| Q-score 计算失败 | `label_qscore.mrc` 全零 |

---

## 四、实测数据

在 4 个覆盖不同类型/分辨率的测试条目上验证：

| 条目 | 类型 | 分辨率 | CC_mask | Q-score | 链数 | Segment 覆盖 | Interface 体素 | Tier |
|------|------|--------|---------|---------|------|-------------|---------------|------|
| EMD-11638 / 7A4M | apoferritin 24聚体 | 1.22Å | **0.80** | **0.56** | 24 | 10.6% | **4,938** | Gold |
| EMD-31083 / 7EFC | 蛋白质 | 1.70Å | **0.79** | **0.87** | 4 | 3.2% | **1,478** | Gold |
| EMD-62134 / 9K6S | 蛋白-RNA 复合物 | 2.80Å | 0.29 | 0.46 | 3 | 0.8% | 0 | Hard |
| EMD-38560 / 8XPS | 糖蛋白 | 3.22Å | **0.80** | **0.65** | 7 | 0.7% | **346** | Silver |

> V3 的 Segment 覆盖率相比 V2 的总覆盖率(含 Voronoi)较低，这是设计预期——V3 关闭了 Voronoi 扩展，只标注距离阈值内的高质量体素。V3 的 CA-only 标签专注残基级精度而非体素覆盖率。

---

## 五、代码结构

```
protein-aligner/
├── run_pipeline.py               # 主管线脚本（8 步）
├── test_pipeline.py              # 4 条目测试脚本
├── visualize_labels.py           # V3 标签可视化（9 列面板）
├── visualize_dataset.py          # V3 数据集概览（6 子图）
├── configs/
│   └── default.yaml              # V3 配置文件
├── scripts/
│   └── install_tools.sh          # 外部工具安装（Merizo/HMMER/Pfam）
├── pipeline/
│   ├── retrieval.py              # Step 1: EMDB/RCSB 数据下载
│   ├── resample.py               # Step 2: 密度图重采样
│   ├── normalization.py          # Step 3: Robust z-score 归一化
│   ├── alignment_qc.py           # Step 4: CC + Q-score 质量控制
│   ├── domain.py                 # Step 5: Merizo 结构域分割（V3 新增）
│   ├── correspondence.py         # Step 6: 三 KDTree + 8 通道标注（V3 核心）
│   ├── enhancement.py            # Step 7: 自适应 blur mol_map
│   ├── redundancy.py             # Step 8: Pfam/MMseqs2 聚类 + 打包
│   ├── pfam.py                   # HMMER + Pfam 注释（V3 新增）
│   ├── interface.py              # CA-CA 界面检测（V3 新增）
│   ├── bio_assembly.py           # 生物组装体展开工具
│   ├── qscore.py                 # Q-score 计算（Pintilie 2020）
│   ├── coord_utils.py            # 坐标转换工具
│   └── config.py                 # 配置加载
├── tools/                        # 外部工具（MMseqs2、Merizo、Pfam）
├── data/
│   ├── raw/                      # 处理中的条目数据
│   └── output/                   # 打包后的 train/val/test 数据集
└── report/                       # 管线报告
```

共 **15 个 Python 模块**，约 **2,800 行**核心代码。

---

## 六、关键指标演进

| 指标 | V1 | V1.5 | V2 | V3 |
|------|-----|------|-----|-----|
| 管线步骤数 | 7 | 7 | 7 | **8** |
| 标签通道数 | 5 | 5 | 6 | **8 + 1** |
| CC_mask (apoferritin) | 0.40 | **0.80** | 0.80 | 0.80 |
| 原子覆盖 (apoferritin) | 3,067 | **73,608** | 73,608 | 73,608 |
| 链数 (apoferritin) | 1 | **24** | 24 | 24 |
| Q-score 支持 | 无 | 无 | **原子级** | **体素级 MRC** |
| 结构域分割 | 无 | 无 | 无 | **Merizo DL** |
| 界面检测 | 无 | 无 | 无 | **CA-CA 跨链** |
| 标注范式 | 全原子 | 全原子+Voronoi | 全原子+Voronoi | **CA-only 三 KDTree** |
| 分子类型 | 蛋白only | 蛋白only | 蛋白/RNA/DNA/糖/配体 | **蛋白专注** |
| 数据划分 | 3-mer Jaccard | 3-mer Jaccard | MMseqs2 | **Pfam Fold/Family** |
| 数据分层 | 无 | 无 | Gold/Silver/Hard | Gold/Silver/Hard |
| 测试条目 | 2 | 2 | 4 | 4 |
| 回退机制 | 无 | 无 | MMseqs2→3-mer | **多层级回退** |

---

## 七、下一步可选方向

1. **安装 Merizo + HMMER + Pfam-A** → 激活 Domain 标签和 Pfam Split
2. **大规模数据生成** → 提交 rlaunch 任务跑 100+ 条目
3. **下游模型训练** → 用 V3 标签训练 model building 网络
4. **Voronoi 对比实验** → 打开/关闭 Voronoi 对下游模型性能的影响
5. **Q-score 加权损失** → 利用 `label_qscore.mrc` 做逐体素样本加权
