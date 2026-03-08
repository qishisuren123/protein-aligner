#!/usr/bin/env python
"""
Cryo-EM 数据处理管线 - 主运行脚本

完整流程:
1. Retrieval: 从 EMDB/RCSB 检索下载数据
2. Resample: 重采样密度图到统一体素大小
3. Normalization: 归一化密度值
4. Alignment & QC: 对齐验证和质量控制
5. Correspondence Labeling: 体素标注（原子/氨基酸/二级结构）
6. Enhancement Labeling: 生成模拟密度图（去噪标签）
7. Redundancy Removal & Packaging: 去冗余、划分数据集、打包

用法:
    python run_pipeline.py [--config configs/default.yaml] [--steps all]
"""
import os
import sys
import time
import json
import logging
import argparse

# 项目根目录
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

# 日志配置
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(PROJECT_ROOT, 'pipeline.log'))
    ]
)
logger = logging.getLogger('pipeline')


def run_pipeline(config, steps='all'):
    """执行完整管线"""
    start_time = time.time()
    logger.info("=" * 60)
    logger.info("Cryo-EM 数据处理管线启动")
    logger.info("=" * 60)

    # 转为绝对路径
    for key in ['raw_dir', 'processed_dir', 'output_dir']:
        if not os.path.isabs(config['paths'][key]):
            config['paths'][key] = os.path.join(PROJECT_ROOT, config['paths'][key])
    os.makedirs(config['paths']['raw_dir'], exist_ok=True)
    os.makedirs(config['paths']['processed_dir'], exist_ok=True)
    os.makedirs(config['paths']['output_dir'], exist_ok=True)

    all_steps = ['retrieval', 'resample', 'normalization', 'alignment_qc',
                 'correspondence', 'enhancement', 'redundancy']

    if steps == 'all':
        run_steps = all_steps
    else:
        run_steps = [s.strip() for s in steps.split(',')]

    entry_dirs = []
    passed_dirs = []

    # ========== Step 1: Retrieval ==========
    if 'retrieval' in run_steps:
        logger.info("\n" + "=" * 50)
        logger.info("Step 1: Retrieval - 数据检索与下载")
        logger.info("=" * 50)
        t0 = time.time()

        retriever = Retriever(config)
        results = retriever.run()

        entry_dirs = [r['dir'] for r in results]
        logger.info(f"Step 1 完成，耗时 {time.time()-t0:.1f}s，下载 {len(entry_dirs)} 个条目")
    else:
        # 从已有数据目录加载
        raw_dir = config['paths']['raw_dir']
        if os.path.exists(raw_dir):
            for d in sorted(os.listdir(raw_dir)):
                full = os.path.join(raw_dir, d)
                if os.path.isdir(full):
                    entry_dirs.append(full)
        logger.info(f"跳过 Retrieval，加载已有 {len(entry_dirs)} 个条目")

    if not entry_dirs:
        logger.error("没有可处理的条目，退出")
        return

    # ========== Step 2: Resample ==========
    if 'resample' in run_steps:
        logger.info("\n" + "=" * 50)
        logger.info("Step 2: Resample - 密度图重采样")
        logger.info("=" * 50)
        t0 = time.time()

        resampler = Resampler(config)
        resampler.run(entry_dirs)
        logger.info(f"Step 2 完成，耗时 {time.time()-t0:.1f}s")

    # ========== Step 3: Normalization ==========
    if 'normalization' in run_steps:
        logger.info("\n" + "=" * 50)
        logger.info("Step 3: Normalization - 密度图归一化")
        logger.info("=" * 50)
        t0 = time.time()

        normalizer = Normalizer(config)
        normalizer.run(entry_dirs)
        logger.info(f"Step 3 完成，耗时 {time.time()-t0:.1f}s")

    # ========== Step 4: Alignment & QC ==========
    if 'alignment_qc' in run_steps:
        logger.info("\n" + "=" * 50)
        logger.info("Step 4: Alignment & QC - 对齐与质量控制")
        logger.info("=" * 50)
        t0 = time.time()

        aligner = AlignmentQC(config)
        results, passed_dirs = aligner.run(entry_dirs)
        logger.info(f"Step 4 完成，耗时 {time.time()-t0:.1f}s")
    else:
        passed_dirs = entry_dirs

    # ========== Step 5: Correspondence Labeling ==========
    if 'correspondence' in run_steps:
        logger.info("\n" + "=" * 50)
        logger.info("Step 5: Correspondence Labeling - 体素标注")
        logger.info("=" * 50)
        t0 = time.time()

        labeler = CorrespondenceLabeler(config)
        labeler.run(passed_dirs)
        logger.info(f"Step 5 完成，耗时 {time.time()-t0:.1f}s")

    # ========== Step 6: Enhancement Labeling ==========
    if 'enhancement' in run_steps:
        logger.info("\n" + "=" * 50)
        logger.info("Step 6: Enhancement Labeling - 增强标注")
        logger.info("=" * 50)
        t0 = time.time()

        enhancer = EnhancementLabeler(config)
        enhancer.run(passed_dirs)
        logger.info(f"Step 6 完成，耗时 {time.time()-t0:.1f}s")

    # ========== Step 7: Redundancy Removal & Packaging ==========
    if 'redundancy' in run_steps:
        logger.info("\n" + "=" * 50)
        logger.info("Step 7: Redundancy Removal & Packaging - 去冗余与打包")
        logger.info("=" * 50)
        t0 = time.time()

        remover = RedundancyRemover(config)
        remover.run(passed_dirs, config['paths']['output_dir'])
        logger.info(f"Step 7 完成，耗时 {time.time()-t0:.1f}s")

    # ========== 总结 ==========
    total_time = time.time() - start_time
    logger.info("\n" + "=" * 60)
    logger.info(f"管线执行完成！总耗时: {total_time:.1f}s ({total_time/60:.1f}min)")
    logger.info(f"处理条目数: {len(entry_dirs)}")
    logger.info(f"通过 QC 条目数: {len(passed_dirs)}")
    logger.info(f"输出目录: {config['paths']['output_dir']}")
    logger.info("=" * 60)


def main():
    parser = argparse.ArgumentParser(description='Cryo-EM 数据处理管线')
    parser.add_argument('--config', type=str, default=None,
                        help='配置文件路径 (默认: configs/default.yaml)')
    parser.add_argument('--steps', type=str, default='all',
                        help='要执行的步骤 (逗号分隔或 "all")')
    parser.add_argument('--max-entries', type=int, default=None,
                        help='覆盖最大下载条目数')
    parser.add_argument('--resolution-max', type=float, default=None,
                        help='覆盖最大分辨率')
    args = parser.parse_args()

    config = load_config(args.config)

    # 命令行覆盖配置
    if args.max_entries is not None:
        config['retrieval']['max_entries'] = args.max_entries
    if args.resolution_max is not None:
        config['retrieval']['resolution_max'] = args.resolution_max

    run_pipeline(config, steps=args.steps)


if __name__ == '__main__':
    main()
