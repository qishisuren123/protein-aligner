"""
模块8a: Pfam Annotation - Pfam 家族/折叠注释

功能：
- 使用 HMMER (hmmscan) 搜索 Pfam-A 数据库
- 解析 hmmscan 输出，提取每条序列的 Pfam family 集合
- 按 Pfam family 分组条目（用于 Fold/Family Split 防止同源泄露）
- 回退策略：HMMER/Pfam 不可用时返回空注释
"""
import os
import json
import logging
import subprocess
import tempfile
from collections import defaultdict

import gemmi

logger = logging.getLogger(__name__)


class PfamAnnotator:
    """HMMER + Pfam-A 蛋白家族注释"""

    def __init__(self, config):
        red_cfg = config.get('redundancy', {})
        pfam_cfg = red_cfg.get('pfam', {})
        self.enabled = pfam_cfg.get('enabled', True)
        self.hmmer_binary = pfam_cfg.get('hmmer_binary', 'hmmscan')
        self.pfam_db = pfam_cfg.get('pfam_db', 'tools/pfam/Pfam-A.hmm')
        self.evalue_threshold = pfam_cfg.get('evalue_threshold', 1e-5)

    def _check_available(self):
        """检查 HMMER 和 Pfam-A 数据库是否可用"""
        # 检查 hmmscan
        try:
            proc = subprocess.run(
                [self.hmmer_binary, '-h'],
                capture_output=True, text=True, timeout=10
            )
            if proc.returncode != 0:
                logger.warning(f"  hmmscan 不可用: returncode={proc.returncode}")
                return False
        except (FileNotFoundError, subprocess.TimeoutExpired):
            logger.warning(f"  hmmscan 未找到或超时: {self.hmmer_binary}")
            return False

        # 检查 Pfam-A 数据库
        db_path = self.pfam_db
        if not os.path.exists(db_path):
            logger.warning(f"  Pfam-A 数据库不存在: {db_path}")
            return False
        # 检查索引文件
        if not os.path.exists(db_path + '.h3m'):
            logger.warning(f"  Pfam-A 索引不存在: {db_path}.h3m (需要运行 hmmpress)")
            return False

        return True

    def _extract_sequences(self, entry_dirs):
        """
        从条目中提取蛋白质序列

        返回:
            sequences: {entry_name: sequence_string}
        """
        sequences = {}
        for entry_dir in entry_dirs:
            model_path = os.path.join(entry_dir, "model.cif")
            if not os.path.exists(model_path):
                continue
            try:
                st = gemmi.read_structure(model_path)
                model = st[0]
                # 取最长链的序列作为代表
                best_seq = ""
                for chain in model:
                    seq = []
                    for residue in chain:
                        if residue.name == "HOH":
                            continue
                        one_letter = gemmi.find_tabulated_residue(residue.name).one_letter_code
                        if one_letter != '?':
                            seq.append(one_letter)
                    chain_seq = ''.join(seq)
                    if len(chain_seq) > len(best_seq):
                        best_seq = chain_seq
                if best_seq:
                    entry_name = os.path.basename(entry_dir)
                    sequences[entry_name] = best_seq
            except Exception as e:
                logger.warning(f"  序列提取失败 {entry_dir}: {e}")
        return sequences

    def run_hmmscan(self, sequences):
        """
        运行 hmmscan 搜索 Pfam-A

        参数:
            sequences: {entry_name: sequence_string}
        返回:
            tblout_path: hmmscan tblout 输出文件路径，失败返回 None
        """
        if not sequences:
            return None

        try:
            tmpdir = tempfile.mkdtemp(prefix="pfam_")

            # 写 FASTA
            fasta_path = os.path.join(tmpdir, "seqs.fasta")
            with open(fasta_path, 'w') as f:
                for name, seq in sequences.items():
                    f.write(f">{name}\n{seq}\n")

            # 运行 hmmscan
            tblout_path = os.path.join(tmpdir, "pfam_tblout.txt")
            cmd = [
                self.hmmer_binary,
                '--tblout', tblout_path,
                '--noali',
                '-E', str(self.evalue_threshold),
                '--cpu', '4',
                self.pfam_db,
                fasta_path,
            ]

            logger.info(f"  运行 hmmscan: {len(sequences)} 条序列")
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=1800
            )

            if proc.returncode != 0:
                logger.warning(f"  hmmscan 失败: {proc.stderr[:500]}")
                return None

            logger.info(f"  hmmscan 完成")
            return tblout_path

        except subprocess.TimeoutExpired:
            logger.warning("  hmmscan 超时 (>30min)")
            return None
        except Exception as e:
            logger.warning(f"  hmmscan 异常: {e}")
            return None

    def parse_hmmscan_output(self, tblout_path):
        """
        解析 hmmscan tblout 输出

        返回:
            entry_pfam_map: {entry_name: set of pfam_family_ids}
        """
        entry_pfam_map = defaultdict(set)

        if not tblout_path or not os.path.exists(tblout_path):
            return dict(entry_pfam_map)

        with open(tblout_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) < 4:
                    continue
                # tblout 格式: target_name accession query_name ...
                pfam_name = fields[0]       # e.g. PF00210.25
                pfam_acc = fields[1]        # e.g. PF00210
                query_name = fields[2]      # entry_name
                entry_pfam_map[query_name].add(pfam_acc)

        return dict(entry_pfam_map)

    def group_by_pfam(self, entry_pfam_map, entry_dirs):
        """
        按 Pfam family 分组条目（同 family 的条目归入同一组）

        使用 Union-Find 将共享 Pfam family 的条目合并为一组。

        参数:
            entry_pfam_map: {entry_name: set of pfam_acc}
            entry_dirs: list of entry_dir paths
        返回:
            groups: {group_id: [entry_dir, ...]}
        """
        name_to_dir = {os.path.basename(d): d for d in entry_dirs}
        entries = list(name_to_dir.keys())

        # Union-Find
        parent = {e: e for e in entries}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry

        # 按 Pfam family 合并（共享同一 family 的条目归入同组）
        pfam_to_entries = defaultdict(list)
        for entry_name, pfam_set in entry_pfam_map.items():
            if entry_name in parent:
                for pfam_acc in pfam_set:
                    pfam_to_entries[pfam_acc].append(entry_name)

        for pfam_acc, entry_names in pfam_to_entries.items():
            for i in range(1, len(entry_names)):
                union(entry_names[0], entry_names[i])

        # 构建分组
        group_map = defaultdict(list)
        for entry_name in entries:
            root = find(entry_name)
            if entry_name in name_to_dir:
                group_map[root].append(name_to_dir[entry_name])

        # 重新编号
        groups = {}
        for i, (root, dirs) in enumerate(group_map.items()):
            groups[i] = dirs

        logger.info(f"  Pfam 分组: {len(entries)} 条目 → {len(groups)} 组")
        return groups

    def annotate_entries(self, entry_dirs):
        """
        主入口：提取序列 → hmmscan → 解析 → 分组

        返回:
            pfam_groups: {group_id: [entry_dir, ...]}，失败返回 None
            entry_pfam_map: {entry_name: [pfam_acc, ...]}
        """
        if not self.enabled:
            logger.info("  Pfam 注释已禁用")
            return None, {}

        if not self._check_available():
            logger.info("  HMMER/Pfam 不可用，跳过 Pfam 注释")
            return None, {}

        # 1. 提取序列
        sequences = self._extract_sequences(entry_dirs)
        if not sequences:
            logger.warning("  无序列可注释")
            return None, {}

        # 2. 运行 hmmscan
        tblout_path = self.run_hmmscan(sequences)
        if tblout_path is None:
            return None, {}

        # 3. 解析输出
        entry_pfam_map = self.parse_hmmscan_output(tblout_path)
        logger.info(f"  Pfam 注释: {len(entry_pfam_map)} 条目有 Pfam hit")

        # 转为可序列化格式
        serializable_map = {k: list(v) for k, v in entry_pfam_map.items()}

        # 4. 分组
        pfam_groups = self.group_by_pfam(entry_pfam_map, entry_dirs)

        return pfam_groups, serializable_map
