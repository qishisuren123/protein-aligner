# Cryo-EM Map-Model Alignment Pipeline V3 — 功能说明

> 面向 cryo-EM model building 的多层级标注数据集生成管线。
> 版本：V3 (2026-03-11)，已在 4 个测试条目上完整验证。

---

## 一、管线总览（8 步）

```
原始密度图 + 原子模型
     │
     ├─ 1. Retrieval          ─ RCSB Search API + EMDB FTP 下载
     ├─ 2. Resample            ─ 重采样到 1.0Å 统一体素
     ├─ 3. Normalization       ─ Robust z-score 归一化 [0, 1]
     ├─ 4. Alignment QC        ─ DensityCalculatorX + CC 三项 + Q-score 持久化
     ├─ 5. Domain Seg.         ─ Merizo 结构域分割（V3 新增）
     ├─ 6. Correspondence      ─ 三 KDTree + CA-only + 8 层体素标注
     ├─ 7. Enhancement         ─ 自适应 blur 生成去噪训练标签
     └─ 8. Redundancy          ─ Pfam Fold/Family Split + MMseqs2 + 打包
```

---

## 二、V3 标签通道（8 + 1）

| # | 文件 | 内容 | 标注范围 | 值域 | 下游用途 |
|---|------|------|---------|------|---------|
| 1 | `label_segment.mrc` | 蛋白/背景 | 所有信号体素 | {0, 1} | 蛋白区域检测 |
| 2 | `label_atom.mrc` | 原子类型 | 所有原子位置 | CA=1,C=2,N=3,O=4,CB=5,Other=6 | 原子级定位 |
| 3 | `label_aa.mrc` | 氨基酸类型 | CA-only | 1-20 | 残基类型预测 |
| 4 | `label_ss.mrc` | 二级结构 | CA-only | Helix=1,Strand=2,Coil=3 | 骨架追踪 |
| 5 | `label_qscore.mrc` | Q-score 置信度 | 骨架原子 | float [-1, 1] | 置信度加权 |
| 6 | `label_chain.mrc` | 链 ID | CA-only | 1-N | 多链分割 |
| 7 | `label_domain.mrc` | 结构域 ID | CA-only | 1-M | 结构域分割 |
| 8 | `label_interface.mrc` | 蛋白-蛋白界面 | CA-only | {0, 1} | 界面预测 |
| 9 | `mol_map.mrc` | 去噪模拟密度 | 所有体素 | float [0, 1] | 密度去噪训练 |

### V3 vs V2 变化

| V2 | V3 | 原因 |
|----|----|------|
| `label_confidence.mrc`（高斯衰减） | `label_qscore.mrc`（物理 Q-score） | 物理级原子置信度替代启发式 |
| `label_moltype.mrc`（5类分子类型） | 删除 | 蛋白专注，无需多分子类型 |
| 无 | `label_segment.mrc` | 新增蛋白/背景二值分割 |
| 无 | `label_domain.mrc` | 新增结构域分割 |
| 无 | `label_interface.mrc` | 新增界面检测 |
| 全原子同一 KDTree | 三组 KDTree | CA-only 减少噪声 |
| Voronoi 默认开启 | 默认关闭 | CA-only 标签无需扩展 |

---

## 三、核心技术

### 3.1 三组 KDTree 架构

| KDTree | 源原子 | 产出标签 |
|--------|--------|---------|
| `tree_all` | 所有蛋白原子 | `label_segment`（距离阈值内=1）, `label_atom`（原子类型） |
| `tree_ca` | 仅 CA 原子 | `label_aa`, `label_ss`, `label_chain`, `label_domain`, `label_interface` |
| `tree_backbone` | CA/C/N/O/CB | `label_qscore`（从 qscores.json 读取浮点值） |

CA-only 标注范式将语义标签聚焦在残基级，减少全原子标注的体素噪声。

### 3.2 生物组装体展开

`gemmi.make_assembly()` 将不对称单元展开为完整生物组装体。

| 条目 | 展开前 | 展开后 |
|------|--------|--------|
| Apoferritin (7A4M) | 3,067 原子 / 1 链 | 73,608 原子 / 24 链 |
| 7EFC | 962 原子 / 1 链 | 3,848 原子 / 4 链 |

### 3.3 DensityCalculatorX 物理级模拟密度

gemmi 内置 5 高斯近似原子散射因子，CC_mask 从 V1 的 0.40 提升到 **0.80**。

### 3.4 Q-score（Pintilie et al., Nature Methods 2020）

- 每个原子 40 个径向方向（Fibonacci sphere）采样实验密度
- 采样范围 0~2Å，步长 0.1Å
- 与参考高斯 `exp(-d²/(2×0.6²))` 做 Pearson 相关
- Step 4 保存逐原子 Q-score 到 `qscores.json`，Step 6 读取生成 `label_qscore.mrc`

### 3.5 Domain Segmentation（Merizo）

- 深度学习蛋白质结构域分割工具
- 在原始结构上运行，通过残基序号复制到展开后的对称链
- 回退：Merizo 不可用时全零

### 3.6 Interface Detection

- CA-CA 跨链距离检测（阈值 5.0Å）
- 每条链建 cKDTree 批量查询
- 高效处理多链蛋白（apoferritin 24 聚体：168/4536 CA 为界面）

### 3.7 数据划分（三级回退）

1. **Pfam Fold/Family Split**（首选）：HMMER hmmscan → Pfam family → Union-Find 分组
2. **MMseqs2 聚类**（回退）：identity=0.3，自动下载二进制
3. **3-mer Jaccard**（兜底）：无外部依赖

### 3.8 数据分层

| 层级 | 条件 | 用途 |
|------|------|------|
| Gold | resolution < 2.5Å AND cc_mask > 0.7 AND q_score > 0.5 | 高质量训练 |
| Silver | resolution < 3.5Å AND cc_mask > 0.5 | 一般训练 |
| Hard | 其余通过 QC | 困难样本 |

---

## 四、测试结果

2026-03-11 在 4 个条目上完整运行 V3 管线（Steps 4-8），总耗时 3.1 分钟。

| 条目 | 类型 | 分辨率 | CC_mask | Q-score | 链数 | Segment覆盖 | Interface体素 | Tier |
|------|------|--------|---------|---------|------|-------------|---------------|------|
| EMD-11638 / 7A4M | apoferritin 24聚体 | 1.22Å | **0.80** | **0.56** | 24 | 10.5% | **4,938** | Gold |
| EMD-31083 / 7EFC | 蛋白质 | 1.70Å | **0.79** | **0.87** | 4 | 3.2% | **1,478** | Gold |
| EMD-38560 / 8XPS | 糖蛋白 | 3.22Å | **0.80** | **0.65** | 7 | 0.7% | **346** | Silver |
| EMD-62134 / 9K6S | 蛋白-RNA | 2.80Å | 0.29 | 0.46 | 3 | 0.8% | 0 | Hard |

文件完整性：4/4 条目全部 V3 完整（18 个文件）。

回退状态：
- Merizo 未安装 → `label_domain.mrc` 全零（预期）
- HMMER/Pfam 未安装 → 回退到 MMseqs2 聚类（预期）

---

## 五、可视化

所有图片位于 `report/figures/`（共 9 张 V3 图）。

| 文件 | 内容 |
|------|------|
| `00_dataset_overview.png` | V3 数据集概览（6 子图） |
| `*_labels_panel.png` × 4 | 9 列标签面板（Density/Segment/Atom/AA/SS/Q-score/Chain/Domain/Interface） |
| `*_label_stats.png` × 4 | 8 子图统计 |

---

## 六、每条目输出文件（18 个）

```
data/raw/EMD-XXXXX_YYYY/
├── raw_map.map              # 原始 EMDB 密度图
├── map_std_1.0A.mrc         # 重采样到 1.0Å
├── map_normalized.mrc       # Robust z-score 归一化
├── sim_map.mrc              # QC 用模拟密度图
├── mol_map.mrc              # 去噪训练标签（自适应 blur）
├── model.cif                # 原子模型 (mmCIF)
├── label_segment.mrc        # 蛋白/背景 (0/1)
├── label_atom.mrc           # 原子类型 (CA=1,C=2,N=3,O=4,CB=5,Other=6)
├── label_aa.mrc             # 氨基酸类型 (1-20)，CA-only
├── label_ss.mrc             # 二级结构 (Helix=1,Strand=2,Coil=3)，CA-only
├── label_qscore.mrc         # Q-score 置信度 (float)，骨架原子
├── label_chain.mrc          # 链 ID (1-N)，CA-only
├── label_domain.mrc         # 结构域 ID (1-M)，CA-only
├── label_interface.mrc      # 界面标签 (0/1)，CA-only
├── qscores.json             # 逐原子 Q-score
├── domain_assignment.json   # 结构域分割结果
├── qc_metrics.json          # CC + Q-score 指标
├── labeling_stats.json      # V3 标注统计
└── metadata.json            # PDB/EMDB ID, 分辨率
```

---

## 七、回退机制

| 组件 | 不可用时的行为 |
|------|--------------|
| Merizo | `label_domain.mrc` 全零，管线继续 |
| HMMER + Pfam-A | 回退到 MMseqs2 聚类 split |
| MMseqs2 | 回退到 3-mer Jaccard 近似聚类 |
| PyTorch CPU | Merizo 不可用，domain 标签全零 |
| 生物组装体信息缺失 | 使用原始不对称单元 |
| Q-score 计算失败 | `label_qscore.mrc` 全零 |
