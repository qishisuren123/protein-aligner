#!/usr/bin/env python
"""
测试脚本 - 下载小样本数据并跑通完整管线

选择几个分辨率高、密度图较小的 EMDB 条目进行测试
"""
import os
import sys
import logging

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

from pipeline.config import load_config
from pipeline.retrieval import Retriever
from pipeline.resample import Resampler
from pipeline.normalization import Normalizer
from pipeline.alignment_qc import AlignmentQC
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

# 手动选择的小型 EMDB 条目（密度图较小，适合快速测试）
TEST_ENTRIES = [
    {"pdb_id": "7A4M", "emdb_id": "EMD-11638"},  # apoferritin, ~1.22Å, 1668 atoms
    {"pdb_id": "7EFC", "emdb_id": "EMD-31083"},  # ~1.7Å, 1001 atoms
]


def main():
    config = load_config()

    # 绝对路径
    for key in ['raw_dir', 'processed_dir', 'output_dir']:
        if not os.path.isabs(config['paths'][key]):
            config['paths'][key] = os.path.join(PROJECT_ROOT, config['paths'][key])
    os.makedirs(config['paths']['raw_dir'], exist_ok=True)
    os.makedirs(config['paths']['output_dir'], exist_ok=True)

    # ===== Step 1: 下载测试数据 =====
    logger.info("=" * 50)
    logger.info("Step 1: 下载测试数据")
    logger.info("=" * 50)

    retriever = Retriever(config)
    entry_dirs = []

    for entry in TEST_ENTRIES:
        logger.info(f"下载: {entry['emdb_id']} / {entry['pdb_id']}")
        entry_dir = retriever.download_entry(entry, output_dir=config['paths']['raw_dir'])
        if entry_dir:
            entry_dirs.append(entry_dir)

    if not entry_dirs:
        logger.error("没有成功下载的条目")
        return

    logger.info(f"成功下载 {len(entry_dirs)} 个条目")

    # ===== Step 2: Resample =====
    logger.info("\n" + "=" * 50)
    logger.info("Step 2: Resample")
    logger.info("=" * 50)

    resampler = Resampler(config)
    resampler.run(entry_dirs)

    # ===== Step 3: Normalization =====
    logger.info("\n" + "=" * 50)
    logger.info("Step 3: Normalization")
    logger.info("=" * 50)

    normalizer = Normalizer(config)
    normalizer.run(entry_dirs)

    # ===== Step 4: Alignment & QC =====
    logger.info("\n" + "=" * 50)
    logger.info("Step 4: Alignment & QC")
    logger.info("=" * 50)

    aligner = AlignmentQC(config)
    results, passed_dirs = aligner.run(entry_dirs)

    # 即使 QC 不通过也继续处理（测试目的）
    if not passed_dirs:
        logger.warning("没有通过 QC 的条目，使用所有条目继续测试")
        passed_dirs = entry_dirs

    # ===== Step 5: Correspondence Labeling =====
    logger.info("\n" + "=" * 50)
    logger.info("Step 5: Correspondence Labeling")
    logger.info("=" * 50)

    labeler = CorrespondenceLabeler(config)
    labeler.run(passed_dirs)

    # ===== Step 6: Enhancement Labeling =====
    logger.info("\n" + "=" * 50)
    logger.info("Step 6: Enhancement Labeling")
    logger.info("=" * 50)

    enhancer = EnhancementLabeler(config)
    enhancer.run(passed_dirs)

    # ===== Step 7: Redundancy Removal & Packaging =====
    logger.info("\n" + "=" * 50)
    logger.info("Step 7: Redundancy Removal & Packaging")
    logger.info("=" * 50)

    remover = RedundancyRemover(config)
    remover.run(passed_dirs, config['paths']['output_dir'])

    # ===== 验证输出 =====
    logger.info("\n" + "=" * 50)
    logger.info("验证输出文件")
    logger.info("=" * 50)

    for entry_dir in passed_dirs:
        dirname = os.path.basename(entry_dir)
        logger.info(f"\n条目: {dirname}")
        expected_files = [
            "raw_map.map", "map_std_1.0A.mrc", "map_normalized.mrc",
            "sim_map.mrc", "mol_map.mrc",
            "label_atom.mrc", "label_aa.mrc", "label_ss.mrc",
            "label_chain.mrc", "label_confidence.mrc",
            "qc_metrics.json", "labeling_stats.json", "metadata.json"
        ]
        for f in expected_files:
            path = os.path.join(entry_dir, f)
            if os.path.exists(path):
                size = os.path.getsize(path)
                logger.info(f"  [OK] {f} ({size/1024/1024:.1f}MB)")
            else:
                logger.warning(f"  [MISSING] {f}")

    logger.info("\n测试完成!")


if __name__ == '__main__':
    main()
