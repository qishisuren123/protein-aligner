"""
模块5a: Domain Segmentation - 结构域分割

功能：
- 封装 Merizo 深度学习工具进行蛋白质结构域分割
- 将 mmCIF 转换为 PDB（Merizo 需要 PDB 输入）
- 解析 Merizo 输出，构建残基→结构域 ID 映射
- 在原始结构上运行，通过残基序号复制到展开后的对称链
- 回退策略：Merizo 不可用时返回空映射
"""
import os
import json
import logging
import subprocess
import tempfile

logger = logging.getLogger(__name__)


class DomainSegmenter:
    """Merizo 结构域分割封装"""

    def __init__(self, config):
        self.config = config
        dom_cfg = config.get('domain_segmentation', {})
        self.enabled = dom_cfg.get('enabled', True)
        self.merizo_path = dom_cfg.get('merizo_path', 'tools/merizo')
        self.device = dom_cfg.get('device', 'cpu')

    def _check_merizo(self):
        """检查 Merizo 是否可用"""
        predict_script = os.path.abspath(os.path.join(self.merizo_path, 'predict.py'))
        if not os.path.exists(predict_script):
            logger.warning(f"  Merizo 不可用: {predict_script} 不存在")
            return False
        try:
            import torch
            return True
        except ImportError:
            logger.warning("  Merizo 不可用: PyTorch 未安装")
            return False

    def _cif_to_pdb(self, cif_path, pdb_path):
        """
        使用 gemmi 将 mmCIF 转换为 PDB 格式（Merizo 需要 PDB 输入）
        """
        import gemmi
        st = gemmi.read_structure(cif_path)
        st.setup_entities()
        st.write_pdb(pdb_path)
        logger.info(f"  mmCIF → PDB 转换完成: {pdb_path}")

    def _run_merizo(self, pdb_path, output_dir):
        """
        运行 Merizo 结构域分割

        Merizo 输出在 {pdb_path_without_ext}_merizo_v2/ 目录下
        使用 --save_domains 保存 domain PDB，--return_indices 保存残基索引

        返回:
            output_file: Merizo 输出目录路径，失败返回 None
        """
        predict_script = os.path.abspath(os.path.join(self.merizo_path, 'predict.py'))
        merizo_dir = os.path.abspath(self.merizo_path)

        cmd = [
            'python', predict_script,
            '-d', self.device,
            '-i', pdb_path,
            '--save_domains',
            '--return_indices',
        ]

        logger.info(f"  运行 Merizo: {' '.join(cmd)}")
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=600,
                cwd=merizo_dir,
            )
            if proc.returncode != 0:
                logger.warning(f"  Merizo 运行失败 (returncode={proc.returncode})")
                if proc.stderr:
                    logger.warning(f"  stderr: {proc.stderr[:500]}")
                return None

            # Merizo 输出目录: {pdb_path_without_ext}_merizo_v2/
            pdb_bn = os.path.splitext(pdb_path)[0]
            merizo_out = pdb_bn + "_merizo_v2"

            # 解析 stdout 获取 chopping 信息
            stdout_text = proc.stdout.strip() if proc.stdout else ""
            logger.info(f"  Merizo stdout: {stdout_text[:200]}")

            return merizo_out

        except subprocess.TimeoutExpired:
            logger.warning("  Merizo 运行超时 (>600s)")
            return None
        except Exception as e:
            logger.warning(f"  Merizo 运行异常: {e}")
            return None

    def _parse_merizo_output(self, output_dir, pdb_name):
        """
        解析 Merizo 输出文件

        优先解析 .idx 文件（格式: residue_num:domain_id,...）
        回退到 domain PDB 文件

        返回:
            domain_map: {(chain, resseq): domain_id} 映射
            n_domains: 结构域总数
            domain_boundaries: 结构域边界信息
        """
        domain_map = {}
        domain_boundaries = {}

        # 方式1：解析 .idx 文件（最可靠）
        idx_path = output_dir + ".idx"
        if os.path.exists(idx_path):
            try:
                with open(idx_path) as f:
                    content = f.read().strip()
                if content:
                    # 格式: "resnum:domid,resnum:domid,..."
                    # resnum 是残基序号，domid 是 domain ID（0=未分配）
                    dom_ranges = {}  # dom_id -> [resseq, ...]
                    for pair in content.split(','):
                        pair = pair.strip()
                        if ':' not in pair:
                            continue
                        resnum_str, domid_str = pair.split(':')
                        resseq = int(float(resnum_str))
                        domid = int(float(domid_str))
                        if domid == 0:
                            continue  # 未分配的残基
                        # Merizo 只处理单链（默认 chain A），使用 "A" 作为链名
                        domain_map[("A", resseq)] = domid
                        if domid not in dom_ranges:
                            dom_ranges[domid] = []
                        dom_ranges[domid].append(resseq)

                    # 构建边界信息
                    for dom_id, residues in dom_ranges.items():
                        domain_boundaries[str(dom_id)] = {
                            "chain": "A",
                            "start": min(residues),
                            "end": max(residues),
                            "n_residues": len(residues),
                        }
                    n_domains = len(set(domain_map.values())) if domain_map else 0
                    logger.info(f"  Merizo .idx 解析: {n_domains} 个结构域, {len(domain_map)} 个残基")
                    return domain_map, n_domains, domain_boundaries
            except Exception as e:
                logger.warning(f"  .idx 文件解析失败: {e}")

        # 方式2：解析 domain PDB 文件
        import glob as globmod
        if os.path.isdir(output_dir):
            domain_pdbs = sorted(globmod.glob(os.path.join(output_dir, f"*domain*.pdb")))
            if not domain_pdbs:
                domain_pdbs = sorted(globmod.glob(os.path.join(output_dir, "*.pdb")))
                domain_pdbs = [p for p in domain_pdbs if pdb_name not in os.path.basename(p)]

            if domain_pdbs:
                import gemmi
                for dom_idx, pdb_file in enumerate(domain_pdbs, start=1):
                    try:
                        st = gemmi.read_structure(pdb_file)
                        model = st[0]
                        res_start = None
                        res_end = None
                        chain_name = None
                        for chain in model:
                            chain_name = chain.name
                            for residue in chain:
                                resseq = residue.seqid.num
                                domain_map[(chain_name, resseq)] = dom_idx
                                if res_start is None or resseq < res_start:
                                    res_start = resseq
                                if res_end is None or resseq > res_end:
                                    res_end = resseq

                        if chain_name and res_start is not None:
                            domain_boundaries[str(dom_idx)] = {
                                "chain": chain_name,
                                "start": res_start,
                                "end": res_end,
                            }
                    except Exception as e:
                        logger.warning(f"  解析 domain PDB 失败 {pdb_file}: {e}")

        # 方式3：解析文本输出（chopping 格式）
        if not domain_map:
            import glob as globmod
            txt_files = globmod.glob(os.path.join(output_dir, "*.txt")) if os.path.isdir(output_dir) else []
            for txt_file in txt_files:
                try:
                    domain_map, domain_boundaries = self._parse_chopping_file(txt_file)
                    if domain_map:
                        break
                except Exception:
                    continue

        n_domains = len(set(domain_map.values())) if domain_map else 0
        logger.info(f"  Merizo 结果: {n_domains} 个结构域, {len(domain_map)} 个残基")
        return domain_map, n_domains, domain_boundaries

    def _parse_chopping_file(self, txt_path):
        """
        解析 Merizo chopping 格式文本文件

        典型格式: A:1-50,A:60-120 表示两个 domain
        """
        domain_map = {}
        domain_boundaries = {}

        with open(txt_path) as f:
            content = f.read().strip()

        if not content:
            return domain_map, domain_boundaries

        # 尝试按行解析
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # 格式: domain_segments 用逗号分隔
            segments = line.split(',')
            for dom_idx, seg in enumerate(segments, start=1):
                seg = seg.strip()
                # 格式: A:1-50 或 1-50
                if ':' in seg:
                    chain_part, range_part = seg.split(':', 1)
                else:
                    chain_part = 'A'
                    range_part = seg

                if '-' in range_part:
                    try:
                        start, end = range_part.split('-')
                        start, end = int(start), int(end)
                        for resseq in range(start, end + 1):
                            domain_map[(chain_part, resseq)] = dom_idx
                        domain_boundaries[str(dom_idx)] = {
                            "chain": chain_part,
                            "start": start,
                            "end": end,
                        }
                    except ValueError:
                        continue

        return domain_map, domain_boundaries

    def segment_entry(self, entry_dir):
        """
        对单个条目进行结构域分割

        策略：在原始结构（未展开）上运行 Merizo，
        然后通过残基序号复制到展开后的对称链（与 SS 映射策略一致）。

        返回:
            domain_data: dict 包含 n_domains, domains 映射等
        """
        model_path = os.path.join(entry_dir, "model.cif")
        output_path = os.path.join(entry_dir, "domain_assignment.json")

        # 空结果模板
        empty_result = {
            "n_domains": 0,
            "domains": {},
            "domain_boundaries": {},
            "method": "none",
        }

        if not self.enabled:
            logger.info(f"  结构域分割已禁用")
            with open(output_path, 'w') as f:
                json.dump(empty_result, f, indent=2)
            return empty_result

        if not os.path.exists(model_path):
            logger.warning(f"  跳过 {entry_dir}: 缺少 model.cif")
            with open(output_path, 'w') as f:
                json.dump(empty_result, f, indent=2)
            return empty_result

        if not self._check_merizo():
            logger.info(f"  Merizo 不可用，跳过结构域分割")
            empty_result["method"] = "fallback_unavailable"
            with open(output_path, 'w') as f:
                json.dump(empty_result, f, indent=2)
            return empty_result

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                # 1. mmCIF → PDB
                pdb_path = os.path.join(tmpdir, "input.pdb")
                self._cif_to_pdb(model_path, pdb_path)

                # 2. 运行 Merizo（输出路径由 Merizo 自动生成在 pdb 文件旁边）
                merizo_out = self._run_merizo(pdb_path, tmpdir)

                if merizo_out is None:
                    logger.warning(f"  Merizo 运行失败，使用空映射")
                    empty_result["method"] = "fallback_merizo_failed"
                    with open(output_path, 'w') as f:
                        json.dump(empty_result, f, indent=2)
                    return empty_result

                # 3. 解析输出（merizo_out = {pdb_without_ext}_merizo_v2）
                domain_map, n_domains, boundaries = self._parse_merizo_output(
                    merizo_out, "input"
                )

            # 4. 构建输出（key 格式: "{chain}_{resseq}"）
            domains_serializable = {}
            for (chain, resseq), dom_id in domain_map.items():
                key = f"{chain}_{resseq}"
                domains_serializable[key] = dom_id

            domain_data = {
                "n_domains": n_domains,
                "domains": domains_serializable,
                "domain_boundaries": boundaries,
                "method": "merizo",
            }

            with open(output_path, 'w') as f:
                json.dump(domain_data, f, indent=2)
            logger.info(f"  结构域分割完成: {n_domains} domains → {output_path}")
            return domain_data

        except Exception as e:
            logger.error(f"  结构域分割异常: {e}", exc_info=True)
            empty_result["method"] = "fallback_error"
            with open(output_path, 'w') as f:
                json.dump(empty_result, f, indent=2)
            return empty_result

    def run(self, entry_dirs):
        """对所有条目进行结构域分割"""
        results = []
        for entry_dir in entry_dirs:
            logger.info(f"结构域分割: {entry_dir}")
            try:
                data = self.segment_entry(entry_dir)
                results.append({"entry_dir": entry_dir, "data": data})
            except Exception as e:
                logger.error(f"  结构域分割失败: {e}", exc_info=True)

        logger.info(f"结构域分割完成: {len(results)} 个条目")
        return results
