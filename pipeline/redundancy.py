"""
模块7: Remove Redundancy & Packaging - 去冗余与打包 (V3)

功能：
- 从蛋白质序列角度进行去冗余（支持 MMseqs2 精确聚类或 3-mer Jaccard 近似聚类）
- V3 新增：Pfam Fold/Family Split 防止结构同源泄露
- Gold/Silver/Hard 数据质量分层
- 划分 train/val/test 集（按聚类划分防止数据泄露）
- 生成数据集统计报告
- 打包成可分发的数据集格式（V3: 8 个标签通道）
"""
import os
import json
import logging
import shutil
import subprocess
import tempfile
import stat
import urllib.request
from collections import defaultdict

import numpy as np
import gemmi

from pipeline.pfam import PfamAnnotator

logger = logging.getLogger(__name__)


class RedundancyRemover:
    """去冗余与数据集打包"""

    def __init__(self, config):
        self.config = config
        red_cfg = config['redundancy']
        self.seq_identity = red_cfg.get('sequence_identity', 0.7)
        self.use_mmseqs2 = red_cfg.get('use_mmseqs2', False)
        self.mmseqs2_identity = red_cfg.get('mmseqs2_identity', 0.3)
        self.mmseqs2_binary = red_cfg.get('mmseqs2_binary', 'tools/mmseqs')
        # 数据分层配置
        self.tier_cfg = red_cfg.get('tiers', {})
        # V3: Pfam split 配置
        pfam_cfg = red_cfg.get('pfam', {})
        self.pfam_enabled = pfam_cfg.get('enabled', False)

    # ===== 序列提取 =====

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

    # ===== 3-mer Jaccard 聚类（原有方法） =====

    def simple_sequence_clustering(self, sequences):
        """
        简化的序列聚类（基于序列相似度）

        使用 3-mer Jaccard 相似系数做近似聚类。
        MMseqs2 不可用时的回退方案。

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

        logger.info(f"序列聚类 (3-mer Jaccard): {len(entries)} 条目 -> {len(clusters)} 簇")
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

    # ===== MMseqs2 精确聚类 =====

    def _download_mmseqs2(self):
        """
        下载 MMseqs2 静态二进制（Linux x86_64）

        返回:
            binary_path: 下载后的可执行文件路径，失败返回 None
        """
        binary_path = os.path.abspath(self.mmseqs2_binary)
        if os.path.exists(binary_path) and os.access(binary_path, os.X_OK):
            logger.info(f"  MMseqs2 已存在: {binary_path}")
            return binary_path

        url = ("https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz")
        os.makedirs(os.path.dirname(binary_path) or '.', exist_ok=True)

        try:
            logger.info(f"  下载 MMseqs2 静态二进制...")
            with tempfile.TemporaryDirectory() as tmpdir:
                tarball = os.path.join(tmpdir, "mmseqs.tar.gz")
                urllib.request.urlretrieve(url, tarball)
                # 解压
                import tarfile
                with tarfile.open(tarball, 'r:gz') as tf:
                    tf.extractall(tmpdir)
                # 找到 mmseqs 可执行文件
                extracted = os.path.join(tmpdir, "mmseqs", "bin", "mmseqs")
                if not os.path.exists(extracted):
                    logger.warning(f"  解压后未找到 mmseqs 二进制")
                    return None
                shutil.copy2(extracted, binary_path)
                os.chmod(binary_path, os.stat(binary_path).st_mode | stat.S_IEXEC)
            logger.info(f"  MMseqs2 下载完成: {binary_path}")
            return binary_path
        except Exception as e:
            logger.warning(f"  MMseqs2 下载失败: {e}")
            return None

    def _run_mmseqs2(self, sequences):
        """
        使用 MMseqs2 进行序列聚类

        参数:
            sequences: {entry_dir: {chain_id: sequence}}
        返回:
            clusters: {cluster_id: [entry_dir, ...]}，失败返回 None
        """
        binary_path = self._download_mmseqs2()
        if binary_path is None:
            return None

        # 准备代表序列（取每个条目最长链）
        representatives = {}
        for entry_dir, chain_seqs in sequences.items():
            if chain_seqs:
                rep_seq = max(chain_seqs.values(), key=len)
                representatives[entry_dir] = rep_seq

        if len(representatives) < 2:
            # 只有一个条目，直接返回单簇
            return {0: list(representatives.keys())}

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                # 写入 FASTA
                fasta_path = os.path.join(tmpdir, "seqs.fasta")
                entry_order = []
                with open(fasta_path, 'w') as f:
                    for i, (entry_dir, seq) in enumerate(representatives.items()):
                        name = os.path.basename(entry_dir)
                        f.write(f">{name}\n{seq}\n")
                        entry_order.append(entry_dir)

                # 运行 mmseqs easy-cluster
                result_prefix = os.path.join(tmpdir, "result")
                mmseqs_tmp = os.path.join(tmpdir, "tmp")
                cmd = [
                    binary_path, "easy-cluster",
                    fasta_path, result_prefix, mmseqs_tmp,
                    "--min-seq-id", str(self.mmseqs2_identity),
                    "-c", "0.8",
                    "--cov-mode", "0",
                ]
                logger.info(f"  运行 MMseqs2: identity={self.mmseqs2_identity}")
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                if proc.returncode != 0:
                    logger.warning(f"  MMseqs2 运行失败: {proc.stderr[:500]}")
                    return None

                # 解析聚类结果: result_prefix_cluster.tsv
                cluster_tsv = result_prefix + "_cluster.tsv"
                if not os.path.exists(cluster_tsv):
                    logger.warning(f"  MMseqs2 聚类结果文件不存在")
                    return None

                # TSV 格式: representative_name \t member_name
                name_to_dir = {os.path.basename(d): d for d in entry_order}
                cluster_map = {}  # rep_name -> [entry_dirs]
                with open(cluster_tsv) as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) < 2:
                            continue
                        rep, member = parts[0], parts[1]
                        if rep not in cluster_map:
                            cluster_map[rep] = []
                        member_dir = name_to_dir.get(member)
                        if member_dir:
                            cluster_map[rep].append(member_dir)

                # 转换为 {cluster_id: [entry_dirs]}
                clusters = {}
                for i, (rep, members) in enumerate(cluster_map.items()):
                    clusters[i] = members

                logger.info(f"序列聚类 (MMseqs2): {len(representatives)} 条目 -> {len(clusters)} 簇")
                return clusters

        except Exception as e:
            logger.warning(f"  MMseqs2 聚类失败: {e}")
            return None

    def mmseqs2_clustering(self, sequences):
        """
        MMseqs2 聚类入口，失败时自动回退到 3-mer Jaccard

        返回:
            clusters: {cluster_id: [entry_dir, ...]}
        """
        clusters = self._run_mmseqs2(sequences)
        if clusters is not None:
            return clusters
        logger.info("  MMseqs2 不可用，回退到 3-mer Jaccard 聚类")
        return self.simple_sequence_clustering(sequences)

    # ===== 数据分层 =====

    def classify_tier(self, entry_dir):
        """
        根据数据质量将条目分为 Gold/Silver/Copper/Hard 四级（V3.1 更新）

        Gold:   resolution < 2.5Å AND cc_mask > 0.80 AND q_score_mean > 0.60
        Silver: resolution < 4.0Å AND cc_mask > 0.70 AND q_score_mean > 0.40
        Copper: cc_mask > 0.60 AND q_score_mean > 0.20
        Hard:   其余通过 QC 的条目
        """
        qc_path = os.path.join(entry_dir, "qc_metrics.json")
        if not os.path.exists(qc_path):
            return "hard"

        with open(qc_path) as f:
            qc = json.load(f)

        resolution = qc.get("resolution", 99.0)
        cc_mask = qc.get("cc_mask", 0.0)
        q_score_mean = qc.get("q_score_mean", 0.0)

        # 可配置阈值
        gold_cfg = self.tier_cfg.get('gold', {})
        silver_cfg = self.tier_cfg.get('silver', {})
        copper_cfg = self.tier_cfg.get('copper', {})

        gold_res = gold_cfg.get('max_resolution', 2.5)
        gold_cc = gold_cfg.get('min_cc_mask', 0.80)
        gold_qs = gold_cfg.get('min_qscore', 0.60)

        silver_res = silver_cfg.get('max_resolution', 4.0)
        silver_cc = silver_cfg.get('min_cc_mask', 0.70)
        silver_qs = silver_cfg.get('min_qscore', 0.40)

        copper_cc = copper_cfg.get('min_cc_mask', 0.60)
        copper_qs = copper_cfg.get('min_qscore', 0.20)

        if resolution < gold_res and cc_mask > gold_cc and q_score_mean > gold_qs:
            return "gold"
        elif resolution < silver_res and cc_mask > silver_cc and q_score_mean > silver_qs:
            return "silver"
        elif cc_mask > copper_cc and q_score_mean > copper_qs:
            return "copper"
        else:
            return "hard"

    # ===== 数据集划分 =====

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

    # ===== 报告生成 =====

    def generate_report(self, entry_dirs, clusters, splits, output_dir):
        """生成数据集统计报告（含 tier 分层信息）"""
        # 数据分层统计
        tier_counts = {"gold": 0, "silver": 0, "copper": 0, "hard": 0}
        tier_map = {}

        for entry_dir in entry_dirs:
            tier = self.classify_tier(entry_dir)
            tier_map[entry_dir] = tier
            tier_counts[tier] += 1

        logger.info(f"数据分层: gold={tier_counts['gold']}, "
                    f"silver={tier_counts['silver']}, copper={tier_counts['copper']}, "
                    f"hard={tier_counts['hard']}")

        report = {
            "total_entries": len(entry_dirs),
            "n_clusters": len(clusters),
            "splits": {k: len(v) for k, v in splits.items()},
            "sequence_identity_threshold": self.seq_identity,
            "tiers": tier_counts,
            "entries": {}
        }

        for entry_dir in entry_dirs:
            # 读取已有的 QC 指标
            qc_path = os.path.join(entry_dir, "qc_metrics.json")
            label_path = os.path.join(entry_dir, "labeling_stats.json")
            meta_path = os.path.join(entry_dir, "metadata.json")

            entry_info = {"dir": entry_dir, "tier": tier_map.get(entry_dir, "hard")}
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

    # ===== 打包 =====

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

        # V3 文件列表（8 个标签通道 + 辅助文件）
        files_to_copy = [
            "map_normalized.mrc", "model.cif", "raw_map.map",
            "sim_map.mrc", "mol_map.mrc",
            # V3: 8 个标签通道
            "label_segment.mrc", "label_atom.mrc", "label_aa.mrc",
            "label_ss.mrc", "label_qscore.mrc", "label_chain.mrc",
            "label_domain.mrc", "label_interface.mrc",
            # 元数据
            "qc_metrics.json", "labeling_stats.json", "metadata.json",
            "qscores.json", "domain_assignment.json",
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

    # ===== 主入口 =====

    def run(self, entry_dirs, output_dir=None):
        """执行完整的去冗余和打包流程"""
        if output_dir is None:
            output_dir = self.config['paths'].get('output_dir', 'data/output')

        # 1. 提取序列
        sequences = self.extract_sequences(entry_dirs)

        # 2. V3: 优先使用 Pfam Fold/Family Split
        pfam_groups = None
        entry_pfam_map = {}
        if self.pfam_enabled:
            try:
                pfam_annotator = PfamAnnotator(self.config)
                pfam_groups, entry_pfam_map = pfam_annotator.annotate_entries(entry_dirs)
            except Exception as e:
                logger.warning(f"  Pfam 注释失败: {e}，回退到序列聚类")

        # 3. 序列聚类（Pfam 不可用时使用）
        if pfam_groups is not None and len(pfam_groups) > 0:
            clusters = pfam_groups
            logger.info(f"使用 Pfam Fold/Family Split: {len(clusters)} 组")
        elif self.use_mmseqs2:
            clusters = self.mmseqs2_clustering(sequences)
        else:
            clusters = self.simple_sequence_clustering(sequences)

        # 4. 划分数据集
        splits = self.split_dataset(clusters)

        # 5. 生成报告（含数据分层 + Pfam 注释）
        report = self.generate_report(entry_dirs, clusters, splits, output_dir)
        # 附加 Pfam 信息到报告
        if entry_pfam_map:
            report["pfam_annotations"] = entry_pfam_map
            report_path = os.path.join(output_dir, "dataset_report.json")
            with open(report_path, 'w') as f:
                json.dump(report, f, indent=2)

        # 6. 打包
        self.package_dataset(splits, output_dir)

        return clusters, splits
