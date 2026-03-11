# 项目文件说明

> 项目路径：`/mnt/shared-storage-user/renyiming/protein-aligner`
> 版本：V3 (2026-03-11)
> 总代码量：~4,750 行（15 个 Python 模块 + 配置 + 脚本）

---

## 目录结构总览

```
protein-aligner/
├── run_pipeline.py               # 管线主入口（8步流程）
├── test_pipeline.py              # 4条目端到端测试
├── visualize_labels.py           # 单条目标签可视化
├── visualize_dataset.py          # 数据集统计概览图
├── configs/
│   └── default.yaml              # 全局配置
├── scripts/
│   └── install_tools.sh          # 外部工具安装
├── pipeline/                     # 核心管线模块（13个）
│   ├── __init__.py
│   ├── config.py                 # 配置加载
│   ├── retrieval.py              # Step 1: 数据下载
│   ├── resample.py               # Step 2: 密度图重采样
│   ├── normalization.py          # Step 3: 归一化
│   ├── alignment_qc.py           # Step 4: 对齐质量控制
│   ├── domain.py                 # Step 5: 结构域分割
│   ├── correspondence.py         # Step 6: 体素标注（核心）
│   ├── enhancement.py            # Step 7: 去噪标签
│   ├── redundancy.py             # Step 8: 去冗余+打包
│   ├── bio_assembly.py           # 工具：生物组装体展开
│   ├── qscore.py                 # 工具：Q-score 计算
│   ├── interface.py              # 工具：界面检测
│   ├── pfam.py                   # 工具：Pfam 注释
│   └── coord_utils.py            # 工具：坐标转换
├── data/
│   ├── raw/                      # 每条目原始+处理数据
│   └── output/                   # 打包后的数据集
├── tools/                        # 外部工具（MMseqs2等）
└── report/                       # 文档与可视化
```

---

## 一、顶层脚本

### `run_pipeline.py`（224 行）

管线主入口。按顺序执行 8 个步骤，支持按步骤选择运行。

```bash
python run_pipeline.py                          # 全部 8 步
python run_pipeline.py --steps all              # 同上
python run_pipeline.py --steps correspondence,enhancement  # 只跑部分步骤
python run_pipeline.py --max-entries 10         # 覆盖最大条目数
```

- 自动将配置中的相对路径转为绝对路径
- 跳过 retrieval 时自动从 `data/raw/` 加载已有条目
- 日志同时输出到终端和 `pipeline.log`

### `test_pipeline.py`（324 行）

端到端测试脚本。在 4 个手选条目上跑完整 V3 管线并验证输出。

测试条目：
| 条目 | 类型 | 分辨率 | 特点 |
|------|------|--------|------|
| EMD-11638 / 7A4M | apoferritin | 1.22Å | 24 聚体，测试组装体展开+界面检测 |
| EMD-31083 / 7EFC | 蛋白质 | 1.70Å | 4 链，常规蛋白质 |
| EMD-62134 / 9K6S | 蛋白-RNA | 2.80Å | 含 RNA（V3 只标蛋白部分） |
| EMD-38560 / 8XPS | 糖蛋白 | 3.22Å | 含糖基修饰，7 链 |

V3 断言检查：
- `label_segment.mrc` 值域 = {0, 1}
- `label_qscore.mrc` 有非零浮点值
- apoferritin 的 `label_interface.mrc` 有 interface 体素
- 全部 18 个预期文件存在

### `visualize_labels.py`（425 行）

单条目可视化。生成两张图：

1. **标签面板**（`vis_labels_panel.png`）：9 列 × 3 行正交切面
   - 列：Density / Segment / Atom / AA / SS / Q-score / Chain / Domain / Interface（+ 可选 mol_map）
   - 行：XY / XZ / YZ 三个正交平面
   - 切面位置：自动选择密度信号质心

2. **统计图**（`vis_label_stats.png`）：2×4 = 8 个子图
   - 原子类型分布、氨基酸分布、二级结构分布、Q-score 直方图
   - 覆盖率分解、链分布、Domain 分布、Interface 统计

```bash
python visualize_labels.py                       # 处理 data/raw 下所有条目
python visualize_labels.py data/raw/EMD-11638_7A4M  # 指定条目
```

### `visualize_dataset.py`（266 行）

数据集级统计概览。从所有条目的 JSON 元数据聚合生成 6 子图：

1. 分辨率分布直方图
2. CC_mask vs 分辨率散点图（按 Gold/Silver/Hard 着色）
3. Domain 数量分布
4. 覆盖率对比（Segment vs CA-only）
5. Q-score 分布
6. Interface 比例分布

输出：`data/output/dataset_overview.png`

---

## 二、管线模块（`pipeline/`）

### `config.py`（15 行）

加载 YAML 配置文件。默认读取 `configs/default.yaml`。

### `retrieval.py`（321 行）— Step 1

从 RCSB/EMDB 下载冷冻电镜数据。

- **RCSB Search API v2** 查询：按分辨率 + 实验方法（EM）筛选
- **EMDB FTP** 下载密度图（`.map.gz`）
- **RCSB** 下载原子模型（`.cif`）
- `update_existing_metadata()`：从 RCSB API 获取分辨率回填到 `metadata.json`

每条目输出：`raw_map.map` + `model.cif` + `metadata.json`

### `resample.py`（104 行）— Step 2

将密度图重采样到统一体素大小（默认 1.0Å）。

- 使用 `gemmi.read_ccp4_map()` 读取保证坐标系正确
- `scipy.ndimage.zoom()` 做三次插值重采样
- 自动计算 zoom 因子：`原始spacing / 目标spacing`
- 更新 MRC header 的 cell 参数

输出：`map_std_1.0A.mrc`

### `normalization.py`（174 行）— Step 3

密度值归一化到 [0, 1]。

三种模式（`config.normalization.method`）：
- **robust_zscore**（默认）：`(x - median) / (1.4826 * MAD)`，对异常值鲁棒
- **zscore**：`(x - mean) / std`
- **percentile**：按百分位裁剪后线性缩放

最后 clip 到 [0, 1] 并保存。

输出：`map_normalized.mrc`

### `alignment_qc.py`（402 行）— Step 4

对齐质量评估 + Q-score 计算。这是最复杂的步骤之一。

流程：
1. 加载归一化密度图和原子模型
2. `bio_assembly.py` 展开生物组装体
3. 检查 CA 原子处的密度值分布（快速诊断对齐质量）
4. `DensityCalculatorX`（gemmi）生成模拟密度图（5 高斯原子散射因子）
5. 计算三项 CC 指标：CC_mask、CC_volume、CC_overall
6. 调用 `qscore.py` 计算逐原子 Q-score
7. **V3 新增**：保存逐原子 Q-score 到 `qscores.json`

QC 判定：`quality_score >= 0.5 OR cc_mask >= 0.3 OR cc_volume >= 0.3`

输出：`sim_map.mrc` + `qc_metrics.json` + `qscores.json`

### `domain.py`（309 行）— Step 5（V3 新增）

Merizo 深度学习结构域分割的封装。

流程：
1. 检查 Merizo 是否可用（`tools/merizo/predict.py`）
2. mmCIF → PDB 转换（Merizo 需要 PDB 输入）
3. subprocess 调用 `predict.py -d cpu`
4. 解析输出文件，构建残基→domain 映射
5. **关键设计**：在原始结构（未展开）上运行，通过残基序号复制到对称链

回退：Merizo 不可用时生成空映射，`label_domain.mrc` 全零，管线不中断。

输出：`domain_assignment.json`

### `correspondence.py`（512 行）— Step 6（V3 核心）

体素级标注，V3 管线的核心模块。

**三组 KDTree 架构**：

| KDTree | 源原子 | 产出标签 |
|--------|--------|---------|
| `tree_all` | 所有蛋白原子 | `label_segment`(距离内=1), `label_atom`(原子类型) |
| `tree_ca` | 仅 CA 原子 | `label_aa`, `label_ss`, `label_chain`, `label_domain`, `label_interface` |
| `tree_backbone` | CA/C/N/O/CB | `label_qscore`(Q-score 浮点值) |

**CA-only 标注范式**：语义标签（氨基酸、二级结构、链、结构域、界面）只在 CA 原子位置标注，减少噪声。

流程：
1. `_parse_structure()` 解析蛋白质原子（只保留 `gemmi.ResidueKind.AA`），返回三组原子列表
2. 加载 `qscores.json` 和 `domain_assignment.json`
3. 调用 `interface.py` 检测跨链界面
4. 建三组 KDTree，查询所有信号体素（`density > 0.05`）的最近邻
5. 按距离阈值（默认 3.0Å）赋值 8 个标签通道
6. 保存 8 个 MRC 文件 + `labeling_stats.json`

标签值域：
- `label_segment`: {0, 1}
- `label_atom`: CA=1, C=2, N=3, O=4, CB=5, Other=6
- `label_aa`: 1-20（标准氨基酸）
- `label_ss`: Helix=1, Strand=2, Coil=3
- `label_qscore`: float [-1, 1]
- `label_chain`: 1-N
- `label_domain`: 1-M
- `label_interface`: {0, 1}

### `enhancement.py`（158 行）— Step 7

生成去噪训练标签 `mol_map.mrc`。

- 使用 `DensityCalculatorX` 从原子模型生成模拟密度
- **自适应 blur**：`blur = base_b_factor × (resolution / 2.0)`
  - 高分辨率（<2Å）保留更多细节
  - 低分辨率（>3Å）增加平滑

输出：`mol_map.mrc`

### `redundancy.py`（533 行）— Step 8

去冗余、数据划分和打包。

1. **序列聚类**（三级回退）：
   - Pfam Fold/Family Split（`pfam.py`，首选）
   - MMseqs2 精确聚类（identity=0.3，自动下载二进制）
   - 3-mer Jaccard 近似聚类（兜底）

2. **数据分层**：
   - Gold：resolution < 2.5Å AND cc_mask > 0.7 AND q_score > 0.5
   - Silver：resolution < 3.5Å AND cc_mask > 0.5
   - Hard：其余通过 QC 的条目

3. **train/val/test 划分**：同一聚类/family 不跨集

4. **打包**：符号链接 → `data/output/{train,val,test}/`

输出：`data/output/dataset_report.json` + 符号链接目录

### `bio_assembly.py`（98 行）— 工具模块

生物组装体展开。

- `gemmi.make_assembly()` 将不对称单元展开为完整生物组装体
- apoferritin：1 链 → 24 链（3,067 → 73,608 原子）
- 安全限制：max_chains=200 防止内存溢出

### `qscore.py`（163 行）— 工具模块

Q-score 原子级密度拟合指标（Pintilie et al., Nature Methods 2020）。

算法：
1. 对每个原子，用 Fibonacci sphere 生成 40 个均匀分布的径向方向
2. 沿每个方向以 0.1Å 步长采样 0~2Å 的实验密度值
3. 与参考高斯 `exp(-d²/(2σ²))`（σ=0.6Å）做 Pearson 相关
4. 40 个方向取平均，返回 [-1, 1]

### `interface.py`（89 行）— 工具模块（V3 新增）

蛋白-蛋白界面检测。

- 按链分组 CA 原子，每条链建 cKDTree
- 批量查询跨链最近 CA-CA 距离
- 距离 < 5.0Å → interface(1)
- 高效处理多链蛋白（如 apoferritin 24 聚体）

### `pfam.py`（266 行）— 工具模块（V3 新增）

HMMER + Pfam 注释。

- 从原子模型提取序列 → FASTA
- `hmmscan --tblout ... Pfam-A.hmm seqs.fasta`
- 解析 tblout 提取 Pfam family
- Union-Find 按共享 family 分组条目
- 回退：HMMER/Pfam 不可用时跳过

### `coord_utils.py`（111 行）— 工具模块

密度图坐标系处理。

- `get_grid_info()`：从 gemmi Grid 提取 origin、spacing、shape
- `fractional_to_grid()` / `grid_to_cartesian()`：坐标转换
- `cartesian_to_grid_index()`：原子坐标 → 体素索引（核心，correspondence 依赖）

---

## 三、配置

### `configs/default.yaml`（137 行）

全局配置文件，控制管线所有参数。关键配置项：

```yaml
retrieval:
  resolution_max: 3.5             # 最大分辨率（Å）
  max_entries: 5                  # 最大下载条目数

resample:
  target_spacing: 1.0             # 目标体素大小（Å）

normalization:
  method: robust_zscore           # 归一化方法

alignment_qc:
  qc_threshold: 0.5               # QC 通过阈值

domain_segmentation:
  enabled: true                   # Merizo 结构域分割开关
  merizo_path: "tools/merizo"
  device: "cpu"

correspondence:
  distance_threshold: 3.0         # 标注距离阈值（Å）
  use_voronoi: false              # V3 默认关闭 Voronoi
  interface:
    distance_threshold: 5.0       # CA-CA 界面阈值（Å）

enhancement:
  adaptive_blur: true             # 自适应 blur

redundancy:
  pfam:
    enabled: true                 # Pfam split 开关
```

---

## 四、脚本

### `scripts/install_tools.sh`（112 行）

外部工具一键安装。安装内容：
1. **PyTorch CPU**：Merizo 的运行依赖
2. **Merizo**：从 GitHub 克隆 + Python 依赖
3. **HMMER**：尝试 conda 安装，失败则手动下载二进制
4. **Pfam-A**：从 EBI FTP 下载数据库（~2GB）+ hmmpress

每个组件独立安装，失败不影响其他组件。

```bash
bash scripts/install_tools.sh
```

---

## 五、数据目录

### `data/raw/EMD-XXXXX_YYYY/`

每个条目包含 18 个文件（V3）：

| 文件 | 类型 | 说明 | 产生步骤 |
|------|------|------|----------|
| `raw_map.map` | 密度图 | EMDB 原始密度图 | Step 1 |
| `model.cif` | 结构 | RCSB mmCIF 原子模型 | Step 1 |
| `metadata.json` | 元数据 | PDB/EMDB ID、分辨率、类型 | Step 1 |
| `map_std_1.0A.mrc` | 密度图 | 重采样到 1.0Å | Step 2 |
| `map_normalized.mrc` | 密度图 | 归一化到 [0,1] | Step 3 |
| `sim_map.mrc` | 密度图 | QC 用模拟密度 | Step 4 |
| `qc_metrics.json` | 元数据 | CC 三项指标 + Q-score 统计 | Step 4 |
| `qscores.json` | 元数据 | 逐原子 Q-score（key=chain_resseq_atom） | Step 4 |
| `domain_assignment.json` | 元数据 | 残基→domain 映射 | Step 5 |
| `label_segment.mrc` | 标签 | 蛋白/背景 {0,1} | Step 6 |
| `label_atom.mrc` | 标签 | 原子类型 CA=1,C=2,N=3,O=4,CB=5,Other=6 | Step 6 |
| `label_aa.mrc` | 标签 | 氨基酸类型 1-20（CA-only） | Step 6 |
| `label_ss.mrc` | 标签 | 二级结构 Helix=1,Strand=2,Coil=3（CA-only） | Step 6 |
| `label_qscore.mrc` | 标签 | Q-score 浮点值（骨架原子） | Step 6 |
| `label_chain.mrc` | 标签 | 链 ID 1-N（CA-only） | Step 6 |
| `label_domain.mrc` | 标签 | 结构域 ID 1-M（CA-only） | Step 6 |
| `label_interface.mrc` | 标签 | 界面 {0,1}（CA-only） | Step 6 |
| `mol_map.mrc` | 标签 | 去噪模拟密度 float [0,1] | Step 7 |
| `labeling_stats.json` | 元数据 | 覆盖率 + 标注统计 | Step 6 |

可视化文件（由 visualize_labels.py 生成）：
| `vis_labels_panel.png` | 图片 | 9 列标签面板 |
| `vis_label_stats.png` | 图片 | 8 子图统计 |

### `data/output/`

打包后的数据集（符号链接）。

```
data/output/
├── train/                         # 训练集
│   ├── EMD-11638_7A4M/ → 符号链接
│   ├── EMD-31083_7EFC/
│   └── EMD-62134_9K6S/
├── val/                           # 验证集
│   └── EMD-38560_8XPS/
├── test/                          # 测试集（当前为空）
├── dataset_report.json            # 完整报告（QC + 标注统计 + 分层 + 划分）
└── dataset_overview.png           # 数据集概览图
```

---

## 六、报告目录

### `report/`

| 文件 | 说明 |
|------|------|
| `capability_report.md` | V1→V3 全能力汇报 |
| `pipeline_report.md` | 管线功能说明（V2 版，待更新） |
| `file_reference.md` | 本文件 |
| `figures/` | 可视化图片 |

### `report/figures/`

| 图片 | 说明 |
|------|------|
| `00_dataset_overview.png` | V3 数据集统计概览（6 子图） |
| `*_labels_panel.png` | V3 标签面板（9 列 × 3 正交切面） |
| `*_label_stats.png` | V3 标注统计（8 子图） |

---

## 七、V3 管线测试结果

2026-03-11 在全部 4 个条目上完整运行 Step 4-8，总耗时 3.1 分钟。

| 条目 | 类型 | 分辨率 | CC_mask | Q-score | 链数 | Segment 覆盖 | Interface 体素 | Tier |
|------|------|--------|---------|---------|------|-------------|---------------|------|
| EMD-11638 / 7A4M | apoferritin 24聚体 | 1.22Å | **0.80** | **0.56** | 24 | 10.5% | **4,938** | Gold |
| EMD-31083 / 7EFC | 蛋白质 | 1.70Å | **0.79** | **0.87** | 4 | 3.2% | **1,478** | Gold |
| EMD-38560 / 8XPS | 糖蛋白 | 3.22Å | **0.80** | **0.65** | 7 | 0.7% | **346** | Silver |
| EMD-62134 / 9K6S | 蛋白-RNA 复合物 | 2.80Å | 0.29 | 0.46 | 3 | 0.8% | 0 | Hard |

文件完整性：4/4 条目全部 V3 完整（18 个文件）。

回退状态：
- Merizo 未安装 → `label_domain.mrc` 全零（预期行为）
- HMMER/Pfam 未安装 → 回退到 MMseqs2 聚类（预期行为）

---

## 八、外部依赖

### 必需

| 包 | 版本 | 用途 |
|----|------|------|
| Python | 3.10+ | |
| gemmi | >= 0.7.0 | 结构解析、密度计算、坐标系 |
| scipy | | KDTree、重采样 |
| numpy | | 向量化计算 |
| mrcfile | | MRC 文件写入 |
| biopython | | 序列提取 |
| requests | | API 调用、文件下载 |
| tqdm | | 进度条 |
| pyyaml | | 配置文件 |
| matplotlib | | 可视化 |

### 可选

| 工具 | 用途 | 不可用时 |
|------|------|---------|
| MMseqs2 | 序列聚类 | 回退到 3-mer Jaccard |
| PyTorch CPU | Merizo 依赖 | Domain 标签全零 |
| Merizo | 结构域分割 | Domain 标签全零 |
| HMMER | Pfam 搜索 | 回退到 MMseqs2 |
| Pfam-A | 数据库 | 回退到 MMseqs2 |

Conda 环境：`cryo-pipeline`
```bash
conda activate cryo-pipeline
```
