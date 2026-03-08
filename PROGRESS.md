# Cryo-EM 数据处理管线 - 进展记录

## 2026-03-09

### 00:10 - 全流程跑通
- **状态**: 7 个模块全部实现并跑通
- 测试数据: EMD-11638/7A4M (apoferritin, 1.22Å), EMD-31083/7EFC (1.7Å)
- 环境: conda env `cryo-pipeline` (Python 3.10)
- 依赖: mrcfile, gemmi, biopython, scipy, numpy, requests, tqdm, pyyaml

### 各模块状态
- [x] **Retrieval**: RCSB Search API + EMDB FTP 下载，支持分辨率/方法筛选
- [x] **Resample**: scipy.ndimage.zoom + gemmi 坐标系保持，可多尺度重采样
- [x] **Normalization**: 第95百分位数归一化到 [0,1]
- [x] **Alignment & QC**: gemmi 坐标插值评估 map-model 拟合，通过率 100%
- [x] **Correspondence Labeling**: 原子/氨基酸/二级结构/链/置信度 5 层标注
- [x] **Enhancement Labeling**: 从原子模型生成模拟密度图（去噪标签）
- [x] **Redundancy Removal & Packaging**: 序列聚类 + train/val/test 划分 + 打包

### QC 结果
| 条目 | CA 密度均值 | CA>0 比例 | CC_mask | 状态 |
|------|-----------|----------|---------|------|
| EMD-11638/7A4M | 0.6691 | 100% | 0.1443 | pass |
| EMD-31083/7EFC | 0.9780 | 100% | 0.2832 | pass |

### 遇到的问题及解决
1. **RCSB Search API 204**: 原查询字段 `em_3d_reconstruction.resolution` 不工作 → 改用 `rcsb_entry_info.resolution_combined` + `exptl.method`
2. **MRC 坐标系错乱**: 手动解析 origin/voxel_size 导致 CA 原子处密度为 0 → 改用 gemmi 的 CCP4 读取（`gemmi.read_ccp4_map()`）正确处理坐标系
3. **numpy.bool_ 不可 JSON 序列化**: 用 `bool()` 包装
4. **gemmi.assign_secondary_structure 不存在**: gemmi 0.7.5 通过解析 mmCIF 文件直接获取 helices/sheets → 改为直接读取 `st.helices`/`st.sheets`
5. **大密度图下载超时**: EMD-11668 (412MB) 超时 → 筛选较小的条目(<100MB)

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
      label_aa.mrc          # 氨基酸类型标签
      label_ss.mrc          # 二级结构标签
      label_chain.mrc       # 链标签
      label_confidence.mrc  # 置信度
      qc_metrics.json       # 质量指标
      labeling_stats.json   # 标注统计
      metadata.json         # 元数据
  output/
    train/ val/ test/       # 按聚类划分
    dataset_report.json     # 数据集报告
```

### 待优化
- CC_mask 数值偏低（0.14-0.28），可能需要优化模拟密度的分辨率匹配
- 序列聚类是简化版（3-mer Jaccard），生产环境应用 MMseqs2
- 未安装 pydssp（备选的二级结构计算工具）
