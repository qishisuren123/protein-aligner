#!/usr/bin/env python
"""
测试脚本 V3 - 下载小样本数据并跑通完整管线

覆盖全类别：
- 蛋白质 (apoferritin 24聚体, 另一个蛋白)
- 蛋白质-RNA 复合物
- 糖蛋白 (含 NAG 修饰)

V3 断言：
- label_segment.mrc 只有 0/1
- label_qscore.mrc 有非零浮点值
- label_interface.mrc 在多链条目中有 interface 体素
- label_domain.mrc 生成（可能全零若 Merizo 不可用）
- 8 个标签通道 + qscores.json + domain_assignment.json
"""
import os
import sys
import json
import logging
import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

from pipeline.config import load_config
from pipeline.retrieval import Retriever
from pipeline.resample import Resampler
from pipeline.normalization import Normalizer
from pipeline.alignment_qc import AlignmentQC
from pipeline.domain import DomainSegmenter
from pipeline.correspondence import CorrespondenceLabeler
from pipeline.enhancement import EnhancementLabeler
from pipeline.redundancy import RedundancyRemover

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(PROJECT_ROOT, 'test_pipeline.log'))
    ]
)
logger = logging.getLogger('test')

# 手动选择的 EMDB 条目，覆盖多种分子类型
TEST_ENTRIES = [
    # 蛋白质
    {"pdb_id": "7A4M", "emdb_id": "EMD-11638", "resolution": 1.22,
     "type": "protein", "desc": "apoferritin 24聚体"},
    {"pdb_id": "7EFC", "emdb_id": "EMD-31083", "resolution": 1.70,
     "type": "protein", "desc": "蛋白质"},
    # 蛋白质-RNA 复合物
    {"pdb_id": "9K6S", "emdb_id": "EMD-62134", "resolution": 2.80,
     "type": "protein-RNA", "desc": "蛋白-RNA 复合物"},
    # 糖蛋白（含 NAG 等糖基修饰）
    {"pdb_id": "8XPS", "emdb_id": "EMD-38560", "resolution": 3.22,
     "type": "glycoprotein", "desc": "糖蛋白"},
]


def _load_mrc_data(path):
    """加载 MRC 文件数据为 numpy 数组"""
    import gemmi
    ccp4 = gemmi.read_ccp4_map(path)
    ccp4.setup(float('nan'))
    return np.array(ccp4.grid, copy=True)


def main():
    config = load_config()

    # 绝对路径
    for key in ['raw_dir', 'processed_dir', 'output_dir']:
        if not os.path.isabs(config['paths'][key]):
            config['paths'][key] = os.path.join(PROJECT_ROOT, config['paths'][key])
    os.makedirs(config['paths']['raw_dir'], exist_ok=True)
    os.makedirs(config['paths']['output_dir'], exist_ok=True)

    # ===== Step 1: 下载测试数据 =====
    logger.info("=" * 60)
    logger.info("Step 1: 下载测试数据 (4 条目)")
    logger.info("=" * 60)

    retriever = Retriever(config)
    entry_dirs = []

    for entry in TEST_ENTRIES:
        logger.info(f"下载: {entry['emdb_id']} / {entry['pdb_id']} "
                    f"({entry['resolution']}Å, {entry['type']})")
        entry_dir = retriever.download_entry(entry, output_dir=config['paths']['raw_dir'])
        if entry_dir:
            entry_dirs.append(entry_dir)

    if not entry_dirs:
        logger.error("没有成功下载的条目")
        return

    logger.info(f"成功下载 {len(entry_dirs)} 个条目")

    # 回填分辨率到 metadata
    retriever.update_existing_metadata(entry_dirs)

    # ===== Step 2: Resample =====
    logger.info("\n" + "=" * 60)
    logger.info("Step 2: Resample")
    logger.info("=" * 60)

    resampler = Resampler(config)
    resampler.run(entry_dirs)

    # ===== Step 3: Normalization =====
    logger.info("\n" + "=" * 60)
    logger.info("Step 3: Normalization (robust z-score)")
    logger.info("=" * 60)

    normalizer = Normalizer(config)
    normalizer.run(entry_dirs)

    # ===== Step 4: Alignment & QC =====
    logger.info("\n" + "=" * 60)
    logger.info("Step 4: Alignment & QC (DensityCalculatorX + Bio Assembly + Q-score)")
    logger.info("=" * 60)

    aligner = AlignmentQC(config)
    results, passed_dirs = aligner.run(entry_dirs)

    # 质量断言
    logger.info("\n--- 质量断言 ---")
    for r in results:
        entry_name = os.path.basename(r['entry_dir'])
        cc_mask = r['metrics']['cc_mask']
        cc_volume = r['metrics']['cc_volume']
        n_atoms = r['metrics']['n_atoms_expanded']
        n_chains = r['metrics']['n_chains']
        q_mean = r['metrics'].get('q_score_mean', 'N/A')
        q_median = r['metrics'].get('q_score_median', 'N/A')
        logger.info(f"  {entry_name}: CC_mask={cc_mask:.4f}, CC_volume={cc_volume:.4f}, "
                    f"atoms={n_atoms}, chains={n_chains}, "
                    f"Q-score mean={q_mean}, median={q_median}")

        # apoferritin 应展开为 24 条链
        if "11638" in entry_name:
            if n_chains >= 20:
                logger.info(f"    [PASS] 链数 {n_chains} >= 20 (apoferritin 24聚体)")
            else:
                logger.warning(f"    [WARN] 链数 {n_chains} < 20 (预期 ~24)")

        # V3: 检查 qscores.json 已生成
        qscore_file = os.path.join(r['entry_dir'], "qscores.json")
        if os.path.exists(qscore_file):
            with open(qscore_file) as f:
                qs_data = json.load(f)
            logger.info(f"    [PASS] qscores.json: {len(qs_data)} 原子")
        else:
            logger.warning(f"    [WARN] qscores.json 未生成")

    # 即使 QC 不通过也继续处理（测试目的）
    if not passed_dirs:
        logger.warning("没有通过 QC 的条目，使用所有条目继续测试")
        passed_dirs = entry_dirs

    # ===== Step 5: Domain Segmentation (V3 新增) =====
    logger.info("\n" + "=" * 60)
    logger.info("Step 5: Domain Segmentation (Merizo)")
    logger.info("=" * 60)

    segmenter = DomainSegmenter(config)
    domain_results = segmenter.run(passed_dirs)

    for dr in domain_results:
        entry_name = os.path.basename(dr['entry_dir'])
        n_domains = dr['data'].get('n_domains', 0)
        method = dr['data'].get('method', 'unknown')
        logger.info(f"  {entry_name}: {n_domains} domains (method={method})")

    # ===== Step 6: Correspondence Labeling (V3) =====
    logger.info("\n" + "=" * 60)
    logger.info("Step 6: Correspondence Labeling V3 (三 KDTree + CA-only + 8 通道)")
    logger.info("=" * 60)

    labeler = CorrespondenceLabeler(config)
    label_results = labeler.run(passed_dirs)

    # V3 覆盖率断言
    logger.info("\n--- V3 标签断言 ---")
    for lr in label_results:
        entry_name = os.path.basename(lr['entry_dir'])
        stats = lr['stats']
        coverage_seg = stats.get('coverage_segment_pct', 0)
        coverage_ca = stats.get('coverage_ca_pct', 0)
        n_interface = stats.get('n_interface_voxels', 0)
        n_domain = stats.get('n_domain_voxels', 0)
        logger.info(f"  {entry_name}: segment覆盖={coverage_seg:.1f}%, "
                    f"CA覆盖={coverage_ca:.1f}%, "
                    f"interface体素={n_interface}, domain体素={n_domain}")

        # 1. label_segment.mrc 只有 0/1
        seg_path = os.path.join(lr['entry_dir'], "label_segment.mrc")
        if os.path.exists(seg_path):
            seg_data = _load_mrc_data(seg_path)
            unique_vals = set(np.unique(np.round(seg_data).astype(int)))
            is_binary = unique_vals <= {0, 1}
            logger.info(f"    [{'PASS' if is_binary else 'FAIL'}] "
                        f"label_segment 二值: {unique_vals}")

        # 2. label_qscore.mrc 有非零浮点值
        qs_path = os.path.join(lr['entry_dir'], "label_qscore.mrc")
        if os.path.exists(qs_path):
            qs_data = _load_mrc_data(qs_path)
            has_nonzero = np.any(qs_data != 0)
            qs_range = (float(qs_data.min()), float(qs_data.max()))
            logger.info(f"    [{'PASS' if has_nonzero else 'WARN'}] "
                        f"label_qscore 非零: {has_nonzero}, 范围: {qs_range}")

        # 3. apoferritin（24链）的 label_interface.mrc 有 interface 体素
        if "11638" in entry_name:
            iface_path = os.path.join(lr['entry_dir'], "label_interface.mrc")
            if os.path.exists(iface_path):
                iface_data = _load_mrc_data(iface_path)
                has_interface = np.any(iface_data > 0)
                n_iface = int((iface_data > 0).sum())
                logger.info(f"    [{'PASS' if has_interface else 'FAIL'}] "
                            f"label_interface 有界面: {has_interface} ({n_iface} 体素)")

        # 4. label_domain.mrc 存在
        dom_path = os.path.join(lr['entry_dir'], "label_domain.mrc")
        if os.path.exists(dom_path):
            dom_data = _load_mrc_data(dom_path)
            n_dom_voxels = int((dom_data > 0).sum())
            logger.info(f"    [INFO] label_domain: {n_dom_voxels} 体素标注")

    # ===== Step 7: Enhancement Labeling =====
    logger.info("\n" + "=" * 60)
    logger.info("Step 7: Enhancement Labeling (自适应 blur)")
    logger.info("=" * 60)

    enhancer = EnhancementLabeler(config)
    enhancer.run(passed_dirs)

    # ===== Step 8: Redundancy Removal & Packaging =====
    logger.info("\n" + "=" * 60)
    logger.info("Step 8: Redundancy Removal & Packaging (含 Pfam)")
    logger.info("=" * 60)

    remover = RedundancyRemover(config)
    remover.run(passed_dirs, config['paths']['output_dir'])

    # ===== 验证输出 =====
    logger.info("\n" + "=" * 60)
    logger.info("验证输出文件 (V3)")
    logger.info("=" * 60)

    for entry_dir in passed_dirs:
        dirname = os.path.basename(entry_dir)
        logger.info(f"\n条目: {dirname}")
        # V3 文件列表
        expected_files = [
            "raw_map.map", "map_std_1.0A.mrc", "map_normalized.mrc",
            "sim_map.mrc", "mol_map.mrc",
            # V3: 8 个标签通道
            "label_segment.mrc", "label_atom.mrc", "label_aa.mrc",
            "label_ss.mrc", "label_qscore.mrc", "label_chain.mrc",
            "label_domain.mrc", "label_interface.mrc",
            # 元数据
            "qc_metrics.json", "labeling_stats.json", "metadata.json",
            "qscores.json", "domain_assignment.json",
        ]
        for f in expected_files:
            path = os.path.join(entry_dir, f)
            if os.path.exists(path):
                size = os.path.getsize(path)
                logger.info(f"  [OK] {f} ({size/1024/1024:.1f}MB)")
            else:
                logger.warning(f"  [MISSING] {f}")

    # ===== 验证 dataset_report.json =====
    report_path = os.path.join(config['paths']['output_dir'], "dataset_report.json")
    if os.path.exists(report_path):
        with open(report_path) as f:
            report = json.load(f)
        tiers = report.get("tiers", {})
        logger.info(f"\n数据分层: {tiers}")
        pfam = report.get("pfam_annotations", {})
        if pfam:
            logger.info(f"Pfam 注释: {len(pfam)} 条目有 Pfam hit")
        for name, info in report.get("entries", {}).items():
            logger.info(f"  {name}: tier={info.get('tier', 'N/A')}")

    # ===== 最终报告 =====
    logger.info("\n" + "=" * 60)
    logger.info("V3 最终质量报告")
    logger.info("=" * 60)

    for entry_dir in passed_dirs:
        dirname = os.path.basename(entry_dir)
        qc_path = os.path.join(entry_dir, "qc_metrics.json")
        stats_path = os.path.join(entry_dir, "labeling_stats.json")

        if os.path.exists(qc_path):
            with open(qc_path) as f:
                qc = json.load(f)
            logger.info(f"\n{dirname} QC:")
            logger.info(f"  CC_mask:      {qc.get('cc_mask', 'N/A')}")
            logger.info(f"  CC_volume:    {qc.get('cc_volume', 'N/A')}")
            logger.info(f"  CC_overall:   {qc.get('cc_overall', 'N/A')}")
            logger.info(f"  Q-score mean: {qc.get('q_score_mean', 'N/A')}")
            logger.info(f"  Q-score med:  {qc.get('q_score_median', 'N/A')}")
            logger.info(f"  原子数:       {qc.get('n_atoms_orig', '?')} -> {qc.get('n_atoms_expanded', '?')}")
            logger.info(f"  链数:         {qc.get('n_chains', '?')}")

        if os.path.exists(stats_path):
            with open(stats_path) as f:
                stats = json.load(f)
            logger.info(f"  segment覆盖:  {stats.get('coverage_segment_pct', 0):.1f}%")
            logger.info(f"  CA覆盖:       {stats.get('coverage_ca_pct', 0):.1f}%")
            logger.info(f"  interface体素: {stats.get('n_interface_voxels', 0)}")
            logger.info(f"  domain体素:   {stats.get('n_domain_voxels', 0)}")

    logger.info("\nV3 测试完成!")


if __name__ == '__main__':
    main()
