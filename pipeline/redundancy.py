"""
模块7: Remove Redundancy & Packaging - 去冗余与打包

功能：
- 从蛋白质序列角度进行去冗余（基于序列相似度聚类）
- 划分 train/val/test 集（按聚类划分防止数据泄露）
- 生成数据集统计报告
- 打包成可分发的数据集格式
"""
import os
import json
import logging
import shutil
import hashlib
from collections import defaultdict

import numpy as np
import gemmi

logger = logging.getLogger(__name__)


class RedundancyRemover:
    """去冗余与数据集打包"""

    def __init__(self, config):
        self.config = config
        self.seq_identity = config['redundancy']['sequence_identity']

    def extract_sequences(self, entry_dirs):
        """
        从所有条目中提取蛋白质序列

        返回:
            sequences: {entry_dir: {chain_id: sequence}}
        """
        sequences = {}
        for entry_dir in entry_dirs:
            model_path = os.path.join(entry_dir, "model.cif")
            if not os.path.exists(model_path):
                continue

            try:
                st = gemmi.read_structure(model_path)
                model = st[0]
                chain_seqs = {}

                for chain in model:
                    seq = []
                    for residue in chain:
                        if residue.name == "HOH":
                            continue
                        # 三字母转单字母
                        one_letter = gemmi.find_tabulated_residue(residue.name).one_letter_code
                        if one_letter != '?':
                            seq.append(one_letter)
                    if seq:
                        chain_seqs[chain.name] = ''.join(seq)

                sequences[entry_dir] = chain_seqs
            except Exception as e:
                logger.warning(f"提取序列失败 {entry_dir}: {e}")

        return sequences

    def simple_sequence_clustering(self, sequences):
        """
        简化的序列聚类（基于序列相似度）

        由于没有 MMseqs2，使用简化的基于序列长度和组成的聚类
        实际生产环境应使用 MMseqs2 进行精确聚类

        返回:
            clusters: {cluster_id: [entry_dir, ...]}
        """
        # 提取每个条目的代表序列（最长链的序列）
        representatives = {}
        for entry_dir, chain_seqs in sequences.items():
            if chain_seqs:
                # 取最长的链作为代表
                rep_seq = max(chain_seqs.values(), key=len)
                representatives[entry_dir] = rep_seq

        # 简化聚类：根据序列长度和简单的 k-mer 相似度
        entries = list(representatives.keys())
        visited = set()
        clusters = {}
        cluster_id = 0

        for i, entry_i in enumerate(entries):
            if entry_i in visited:
                continue

            cluster = [entry_i]
            visited.add(entry_i)
            seq_i = representatives[entry_i]

            for j in range(i + 1, len(entries)):
                entry_j = entries[j]
                if entry_j in visited:
                    continue

                seq_j = representatives[entry_j]

                # 简单的序列相似度估计
                sim = self._estimate_similarity(seq_i, seq_j)
                if sim >= self.seq_identity:
                    cluster.append(entry_j)
                    visited.add(entry_j)

            clusters[cluster_id] = cluster
            cluster_id += 1

        logger.info(f"序列聚类: {len(entries)} 条目 -> {len(clusters)} 簇")
        return clusters

    def _estimate_similarity(self, seq1, seq2):
        """
        估计两条序列的相似度（简化版）
        使用 3-mer 组成的 Jaccard 相似系数
        """
        if len(seq1) < 3 or len(seq2) < 3:
            return 0.0

        # 长度差异太大的直接排除
        len_ratio = min(len(seq1), len(seq2)) / max(len(seq1), len(seq2))
        if len_ratio < 0.5:
            return 0.0

        # 3-mer 集合
        kmers1 = set(seq1[i:i+3] for i in range(len(seq1) - 2))
        kmers2 = set(seq2[i:i+3] for i in range(len(seq2) - 2))

        if not kmers1 or not kmers2:
            return 0.0

        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)
        return intersection / union

    def split_dataset(self, clusters, train_ratio=0.8, val_ratio=0.1, test_ratio=0.1):
        """
        按聚类划分数据集（同一聚类的条目不会跨集划分）

        参数:
            clusters: 聚类结果
            train_ratio/val_ratio/test_ratio: 划分比例
        返回:
            splits: {"train": [...], "val": [...], "test": [...]}
        """
        cluster_ids = sorted(clusters.keys())
        np.random.seed(42)
        np.random.shuffle(cluster_ids)

        n_clusters = len(cluster_ids)
        n_train = max(1, int(n_clusters * train_ratio))
        n_val = max(1, int(n_clusters * val_ratio))

        train_clusters = cluster_ids[:n_train]
        val_clusters = cluster_ids[n_train:n_train + n_val]
        test_clusters = cluster_ids[n_train + n_val:]

        splits = {"train": [], "val": [], "test": []}
        for cid in train_clusters:
            splits["train"].extend(clusters[cid])
        for cid in val_clusters:
            splits["val"].extend(clusters[cid])
        for cid in test_clusters:
            splits["test"].extend(clusters[cid])

        logger.info(f"数据集划分: train={len(splits['train'])}, "
                    f"val={len(splits['val'])}, test={len(splits['test'])}")
        return splits

    def generate_report(self, entry_dirs, clusters, splits, output_dir):
        """生成数据集统计报告"""
        report = {
            "total_entries": len(entry_dirs),
            "n_clusters": len(clusters),
            "splits": {k: len(v) for k, v in splits.items()},
            "sequence_identity_threshold": self.seq_identity,
            "entries": {}
        }

        for entry_dir in entry_dirs:
            # 读取已有的 QC 指标
            qc_path = os.path.join(entry_dir, "qc_metrics.json")
            label_path = os.path.join(entry_dir, "labeling_stats.json")
            meta_path = os.path.join(entry_dir, "metadata.json")

            entry_info = {"dir": entry_dir}
            for path, key in [(qc_path, "qc"), (label_path, "labeling"), (meta_path, "metadata")]:
                if os.path.exists(path):
                    with open(path) as f:
                        entry_info[key] = json.load(f)

            # 确定所属 split
            for split_name, dirs in splits.items():
                if entry_dir in dirs:
                    entry_info["split"] = split_name
                    break

            # 确定所属 cluster
            for cid, dirs in clusters.items():
                if entry_dir in dirs:
                    entry_info["cluster_id"] = cid
                    break

            dirname = os.path.basename(entry_dir)
            report["entries"][dirname] = entry_info

        # 保存报告
        report_path = os.path.join(output_dir, "dataset_report.json")
        os.makedirs(output_dir, exist_ok=True)
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)

        logger.info(f"数据集报告已保存: {report_path}")
        return report

    def package_dataset(self, splits, output_dir):
        """
        将数据集打包成标准格式

        目录结构:
        output_dir/
          train/
            EMD-XXXXX_YYYY/
              map_normalized.mrc
              model.cif
              label_*.mrc
              mol_map.mrc
              ...
          val/
          test/
          dataset_report.json
        """
        os.makedirs(output_dir, exist_ok=True)

        # 需要打包的文件列表
        files_to_copy = [
            "map_normalized.mrc", "model.cif", "raw_map.map",
            "sim_map.mrc", "mol_map.mrc",
            "label_atom.mrc", "label_aa.mrc", "label_ss.mrc",
            "label_chain.mrc", "label_confidence.mrc",
            "qc_metrics.json", "labeling_stats.json", "metadata.json"
        ]

        for split_name, entry_dirs in splits.items():
            split_dir = os.path.join(output_dir, split_name)
            os.makedirs(split_dir, exist_ok=True)

            for entry_dir in entry_dirs:
                dirname = os.path.basename(entry_dir)
                dest_dir = os.path.join(split_dir, dirname)
                os.makedirs(dest_dir, exist_ok=True)

                for fname in files_to_copy:
                    src = os.path.join(entry_dir, fname)
                    dst = os.path.join(dest_dir, fname)
                    if os.path.exists(src):
                        # 使用符号链接节省空间
                        if os.path.exists(dst):
                            os.remove(dst)
                        os.symlink(os.path.abspath(src), dst)

                logger.info(f"  打包: {dirname} -> {split_name}")

        logger.info(f"数据集打包完成: {output_dir}")

    def run(self, entry_dirs, output_dir=None):
        """执行完整的去冗余和打包流程"""
        if output_dir is None:
            output_dir = self.config['paths'].get('output_dir', 'data/output')

        # 1. 提取序列
        sequences = self.extract_sequences(entry_dirs)

        # 2. 序列聚类
        clusters = self.simple_sequence_clustering(sequences)

        # 3. 划分数据集
        splits = self.split_dataset(clusters)

        # 4. 生成报告
        self.generate_report(entry_dirs, clusters, splits, output_dir)

        # 5. 打包
        self.package_dataset(splits, output_dir)

        return clusters, splits
