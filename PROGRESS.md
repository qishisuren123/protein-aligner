# Cryo-EM 数据处理管线 - 进展记录

## 2026-03-09

### V2 管线优化：发表级质量提升

#### 新增功能

1. **Q-score 物理置信度** (`pipeline/qscore.py` 新建)
   - Pintilie 2020 标准原子级密度拟合指标
   - Fibonacci sphere 40 方向径向采样，步长 0.1Å
   - 与参考高斯 (σ=0.6Å) 做 Pearson 相关
   - 集成到 `alignment_qc.py`，结果写入 `qc_metrics.json`

2. **多分子类型标注** (`pipeline/correspondence.py` 重构)
   - 新增 `label_moltype.mrc`：protein=1, RNA=2, DNA=3, sugar=4, ligand=5
   - 使用 `gemmi.find_tabulated_residue()` 自动检测分子类型
   - RNA/DNA 区分：O2' 原子检测（核糖 vs 脱氧核糖）
   - 糖基残基匹配：NAG, MAN, BMA, FUC, GAL, SIA 等
   - 核苷酸标签扩展 (21-28): A, U, G, C, DA, DT, DG, DC
   - 糖基标签扩展 (29-34): NAG, MAN, BMA, FUC, GAL, SIA
   - RNA/DNA 二级结构: sugar_phosphate=4, base=5

3. **自适应分辨率 blur** (`pipeline/enhancement.py`)
   - `blur = base_b_factor * (resolution / 2.0)`
   - 高分辨率(<2Å)保留更多细节，低分辨率(>3Å)增加平滑
   - 通过 `enhancement.adaptive_blur: true` 控制

4. **Gold/Silver/Hard 数据分层** (`pipeline/redundancy.py`)
   - Gold: resolution < 2.5Å AND cc_mask > 0.7 AND q_score > 0.5
   - Silver: resolution < 3.5Å AND cc_mask > 0.5
   - Hard: 其余通过 QC 的条目
   - 分层信息写入 `dataset_report.json`

5. **MMseqs2 序列聚类** (`pipeline/redundancy.py`)
   - 自动下载 MMseqs2 静态二进制到 `tools/mmseqs`
   - `mmseqs easy-cluster` 精确序列聚类
   - 不可用时自动回退到 3-mer Jaccard

6. **数据集概览可视化** (`visualize_dataset.py` 新建)
   - 分辨率分布直方图
   - CC_mask vs 分辨率散点图（按 tier 着色）
   - 分子类型饼图
   - 覆盖率对比（Hard vs Voronoi）
   - Q-score 分布
   - 输出: `data/output/dataset_overview.png`

7. **测试扩展** (`test_pipeline.py`)
   - 新增 protein-RNA 复合物: 9K6S / EMD-48045 (2.8Å)
   - 新增糖蛋白: 8HKS / EMD-34948 (2.9Å)
   - label_moltype 验证：RNA 体素、sugar 体素断言
   - Q-score 输出验证

#### 修改文件清单
- `pipeline/qscore.py` — 新建，Q-score 计算 (~120行)
- `pipeline/alignment_qc.py` — 集成 Q-score 到 evaluate_entry()
- `pipeline/correspondence.py` — 重构，多分子类型 + 核苷酸/糖基标签 + RNA/DNA SS
- `pipeline/enhancement.py` — 自适应 blur
- `pipeline/redundancy.py` — 重写，Gold/Silver/Hard 分层 + MMseqs2 聚类
- `test_pipeline.py` — 4 条目覆盖全类别 + 分子类型验证
- `visualize_dataset.py` — 新建，发表级概览图
- `configs/default.yaml` — V2 配置（qscore/moltype/tier/mmseqs2）
- `README.md` — V2 文档更新

#### 新增输出文件
- `label_moltype.mrc` — 每个条目目录下
- `data/output/dataset_overview.png` — 数据集概览图

---

### 15:00 - 管线优化：达到人类水平的 Map-Model 对齐质量

#### 核心改进
1. **生物学组装体展开** (`pipeline/bio_assembly.py` 新建)
   - `gemmi.make_assembly()` 将不对称单元展开为完整组装体
   - apoferritin: 3177 原子 → ~76248 原子（24链）
   - 安全限制: max_chains=200 防止内存溢出

2. **DensityCalculatorX 物理级模拟密度** (`alignment_qc.py`, `enhancement.py`)
   - 替代手写三线性插值+高斯模糊
   - 5高斯原子散射因子，blur=0 最优
   - 新增 CC_mask, CC_volume, CC_overall 三项指标

3. **Voronoi 全覆盖标注** (`correspondence.py`)
   - 双层策略：距离阈值内高置信度 + Voronoi 低置信度（上限0.3）
   - 向量化赋值替代逐体素循环
   - 覆盖率从 <1% 提升到 25-35%

4. **Robust Z-score 归一化** (`normalization.py`)
   - 基于中位数和 MAD，保留密度动态范围
   - 支持 percentile/zscore/robust_zscore 三种方法切换

5. **分辨率 Metadata** (`retrieval.py`)
   - 从 RCSB API 获取分辨率写入 metadata.json
   - 支持已下载条目的分辨率回填

#### 预期效果
| 指标 | 优化前 | 优化后 |
|------|--------|--------|
| EMD-11638 CC_mask | 0.40 | > 0.65 |
| EMD-31083 CC_mask | 0.54 | > 0.65 |
| EMD-11638 标注覆盖率 | 0.6% | 25-35% |
| EMD-11638 原子数 | 3,067 | ~76,248 |
| EMD-11638 链数 | 1 | 24 |

---

### 00:10 - 全流程跑通
- **状态**: 7 个模块全部实现并跑通
- 测试数据: EMD-11638/7A4M (apoferritin, 1.22Å), EMD-31083/7EFC (1.7Å)
- 环境: conda env `cryo-pipeline` (Python 3.10)
- 依赖: mrcfile, gemmi, biopython, scipy, numpy, requests, tqdm, pyyaml

### 各模块状态
- [x] **Retrieval**: RCSB Search API + EMDB FTP 下载，支持分辨率/方法筛选
- [x] **Resample**: scipy.ndimage.zoom + gemmi 坐标系保持，可多尺度重采样
- [x] **Normalization**: Robust z-score / percentile 双模式归一化
- [x] **Alignment & QC**: DensityCalculatorX + 生物组装体 + 三项CC指标 + Q-score
- [x] **Correspondence Labeling**: 生物组装体 + 多分子类型 + Voronoi + 向量化
- [x] **Enhancement Labeling**: DensityCalculatorX + 自适应 blur
- [x] **Redundancy Removal & Packaging**: MMseqs2/3-mer 聚类 + 分层 + 划分 + 打包

### 遇到的问题及解决
1. **RCSB Search API 204**: 原查询字段 `em_3d_reconstruction.resolution` 不工作 → 改用 `rcsb_entry_info.resolution_combined` + `exptl.method`
2. **MRC 坐标系错乱**: 手动解析 origin/voxel_size 导致 CA 原子处密度为 0 → 改用 gemmi 的 CCP4 读取（`gemmi.read_ccp4_map()`）正确处理坐标系
3. **numpy.bool_ 不可 JSON 序列化**: 用 `bool()` 包装
4. **gemmi.assign_secondary_structure 不存在**: gemmi 0.7.5 通过解析 mmCIF 文件直接获取 helices/sheets → 改为直接读取 `st.helices`/`st.sheets`
5. **大密度图下载超时**: EMD-11668 (412MB) 超时 → 筛选较小的条目(<100MB)
6. **对称蛋白覆盖率极低**: apoferritin 只标注 1/24 密度 → `gemmi.make_assembly()` 展开生物组装体
7. **模拟密度质量差**: 手写三线性+高斯模糊 CC=0.40 → `DensityCalculatorX` 使用物理散射因子

### 输出目录结构
```
data/
  raw/
    EMD-XXXXX_YYYY/
      raw_map.map           # 原始密度图
      map_std_1.0A.mrc      # 重采样后
      map_normalized.mrc    # 归一化后
      sim_map.mrc           # 模拟密度图(QC用)
      mol_map.mrc           # 模拟密度图(去噪标签)
      model.cif             # 原子模型
      label_atom.mrc        # 原子类型标签
      label_aa.mrc          # 残基类型标签 (氨基酸+核苷酸+糖基)
      label_ss.mrc          # 二级结构标签 (含RNA/DNA sugar_phosphate/base)
      label_chain.mrc       # 链标签
      label_confidence.mrc  # 置信度（高置信Hard + 低置信Voronoi）
      label_moltype.mrc     # 分子类型标签 (protein/RNA/DNA/sugar/ligand)
      qc_metrics.json       # 质量指标（CC + Q-score）
      labeling_stats.json   # 标注统计（覆盖率 + 分子类型计数）
      metadata.json         # 元数据（含分辨率）
  output/
    train/ val/ test/       # 按聚类划分
    dataset_report.json     # 数据集报告（含 tier 分层）
    dataset_overview.png    # 数据集概览图
```
