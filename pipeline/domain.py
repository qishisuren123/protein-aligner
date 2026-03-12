"""
模块5a: Domain Segmentation - 结构域分割

功能：
- 封装 Merizo 深度学习工具进行蛋白质结构域分割
- 将 mmCIF 转换为 PDB（Merizo 需要 PDB 输入）
- 对 ASU 中每条蛋白链分别运行 Merizo
- 通过残基序号将 domain 结果复制到展开后的对称链，ID 递增
- 回退策略：Merizo 不可用时返回空映射
"""
import os
import json
import logging
import shutil
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
        """使用 gemmi 将 mmCIF 转换为 PDB 格式（Merizo 需要 PDB 输入）"""
        import gemmi
        st = gemmi.read_structure(cif_path)
        st.setup_entities()
        st.write_pdb(pdb_path)
        logger.info(f"  mmCIF → PDB 转换完成: {pdb_path}")

    def _run_merizo(self, pdb_path, output_dir, pdb_chain='A'):
        """
        运行 Merizo 结构域分割

        返回:
            Merizo 输出目录路径（{pdb_without_ext}_merizo_v2），失败返回 None
        """
        predict_script = os.path.abspath(os.path.join(self.merizo_path, 'predict.py'))
        merizo_dir = os.path.abspath(self.merizo_path)

        cmd = [
            'python', predict_script,
            '-d', self.device,
            '-i', pdb_path,
            '--save_domains',
            '--return_indices',
            '--pdb_chain', pdb_chain,
        ]

        logger.info(f"  运行 Merizo (chain {pdb_chain}): {' '.join(cmd)}")
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

            pdb_bn = os.path.splitext(pdb_path)[0]
            merizo_out = pdb_bn + "_merizo_v2"

            stdout_text = proc.stdout.strip() if proc.stdout else ""
            logger.info(f"  Merizo stdout: {stdout_text[:200]}")

            return merizo_out

        except subprocess.TimeoutExpired:
            logger.warning("  Merizo 运行超时 (>600s)")
            return None
        except Exception as e:
            logger.warning(f"  Merizo 运行异常: {e}")
            return None

    def _parse_merizo_output(self, output_dir, pdb_name, target_chain="A"):
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
                    dom_ranges = {}
                    for pair in content.split(','):
                        pair = pair.strip()
                        if ':' not in pair:
                            continue
                        resnum_str, domid_str = pair.split(':')
                        resseq = int(float(resnum_str))
                        domid = int(float(domid_str))
                        if domid == 0:
                            continue
                        domain_map[(target_chain, resseq)] = domid
                        if domid not in dom_ranges:
                            dom_ranges[domid] = []
                        dom_ranges[domid].append(resseq)

                    for dom_id, residues in dom_ranges.items():
                        domain_boundaries[str(dom_id)] = {
                            "chain": target_chain,
                            "start": min(residues),
                            "end": max(residues),
                            "n_residues": len(residues),
                        }
                    n_domains = len(set(domain_map.values())) if domain_map else 0
                    logger.info(f"  Merizo .idx 解析 (chain {target_chain}): "
                                f"{n_domains} 个结构域, {len(domain_map)} 个残基")
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
                        for chain in model:
                            for residue in chain:
                                resseq = residue.seqid.num
                                domain_map[(target_chain, resseq)] = dom_idx
                                if res_start is None or resseq < res_start:
                                    res_start = resseq
                                if res_end is None or resseq > res_end:
                                    res_end = resseq

                        if res_start is not None:
                            domain_boundaries[str(dom_idx)] = {
                                "chain": target_chain,
                                "start": res_start,
                                "end": res_end,
                            }
                    except Exception as e:
                        logger.warning(f"  解析 domain PDB 失败 {pdb_file}: {e}")

        # 方式3：解析文本输出（chopping 格式）
        if not domain_map:
            import glob as globmod
            txt_files = globmod.glob(os.path.join(output_dir, "*.txt")) \
                if os.path.isdir(output_dir) else []
            for txt_file in txt_files:
                try:
                    domain_map, domain_boundaries = self._parse_chopping_file(txt_file)
                    if domain_map:
                        break
                except Exception:
                    continue

        n_domains = len(set(domain_map.values())) if domain_map else 0
        logger.info(f"  Merizo 结果 (chain {target_chain}): "
                    f"{n_domains} 个结构域, {len(domain_map)} 个残基")
        return domain_map, n_domains, domain_boundaries

    def _parse_chopping_file(self, txt_path):
        """解析 Merizo chopping 格式文本文件"""
        domain_map = {}
        domain_boundaries = {}

        with open(txt_path) as f:
            content = f.read().strip()

        if not content:
            return domain_map, domain_boundaries

        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            segments = line.split(',')
            for dom_idx, seg in enumerate(segments, start=1):
                seg = seg.strip()
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

    def _expand_domains_to_assembly(self, model_path, asu_domain_map, n_asu_domains):
        """
        将 ASU 的 domain 结果展开到整个组装体

        每个 ASU copy 的 domain ID 递增：
        若 ASU 有 N 个 domain，则 copy 0 保持 1..N，copy 1 为 N+1..2N，依此类推

        返回:
            expanded_map: {(chain, resseq): dom_id} 展开后的映射
            total_domains: 总 domain 数
        """
        import gemmi
        from pipeline.bio_assembly import expand_to_assembly

        bio_cfg = self.config.get('bio_assembly', {})
        if not bio_cfg.get('enabled', True):
            return asu_domain_map, n_asu_domains

        st = gemmi.read_structure(model_path)
        st.setup_entities()

        # ASU 各链的蛋白残基序号集合
        orig_chain_resseqs = {}
        for chain in st[0]:
            resseqs = set()
            for res in chain:
                info = gemmi.find_tabulated_residue(res.name)
                if info.kind == gemmi.ResidueKind.AA:
                    resseqs.add(res.seqid.num)
            if resseqs:
                orig_chain_resseqs[chain.name] = resseqs

        # ASU domain 按原始链分组: {chain_name: {resseq: dom_id}}
        asu_by_chain = {}
        for (chain, resseq), dom_id in asu_domain_map.items():
            if chain not in asu_by_chain:
                asu_by_chain[chain] = {}
            asu_by_chain[chain][resseq] = dom_id

        # 展开组装体
        expanded_st = expand_to_assembly(
            st,
            assembly_id=bio_cfg.get('assembly_id', '1'),
            max_chains=bio_cfg.get('max_chains', 200),
        )

        # 对展开后的每条链，匹配原始链并复制 domain ID
        expanded_map = {}
        copy_counter = {}

        for chain in expanded_st[0]:
            exp_resseqs = set()
            for res in chain:
                info = gemmi.find_tabulated_residue(res.name)
                if info.kind == gemmi.ResidueKind.AA:
                    exp_resseqs.add(res.seqid.num)

            if not exp_resseqs:
                continue

            # 找匹配度最高的原始链
            best_match = None
            best_overlap = 0
            for orig_name, orig_resseqs in orig_chain_resseqs.items():
                overlap = len(exp_resseqs & orig_resseqs)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_match = orig_name

            if best_match is None or best_match not in asu_by_chain:
                continue

            # copy index
            if best_match not in copy_counter:
                copy_counter[best_match] = 0
            copy_idx = copy_counter[best_match]
            copy_counter[best_match] += 1

            # 复制 domain IDs，偏移量 = copy_idx * n_asu_domains
            offset = copy_idx * n_asu_domains
            for resseq, dom_id in asu_by_chain[best_match].items():
                expanded_map[(chain.name, resseq)] = dom_id + offset

        total_domains = max(expanded_map.values()) if expanded_map else 0
        n_copies = max(copy_counter.values()) if copy_counter else 0
        logger.info(f"  Domain 展开: ASU {n_asu_domains} domains × "
                    f"{n_copies} copies → {total_domains} total domains")

        return expanded_map, total_domains

    def segment_entry(self, entry_dir):
        """
        对单个条目进行结构域分割

        策略：
        1. 对 ASU 中每条蛋白链分别运行 Merizo
        2. 合并结果（domain ID 全局唯一）
        3. 通过残基序号复制到展开后的对称链，每个 copy 的 ID 递增
        """
        model_path = os.path.join(entry_dir, "model.cif")
        output_path = os.path.join(entry_dir, "domain_assignment.json")

        empty_result = {
            "n_domains": 0,
            "n_asu_domains": 0,
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
            # Phase 1: 在 ASU 上运行 Merizo（每条蛋白链分别处理）
            asu_domain_map = {}
            next_dom_id = 1

            with tempfile.TemporaryDirectory() as tmpdir:
                # mmCIF → PDB
                pdb_path = os.path.join(tmpdir, "input.pdb")
                self._cif_to_pdb(model_path, pdb_path)

                # 枚举 ASU 中的蛋白链
                import gemmi
                st = gemmi.read_structure(pdb_path)
                asu_chains = []
                for chain in st[0]:
                    has_protein = any(
                        gemmi.find_tabulated_residue(res.name).kind == gemmi.ResidueKind.AA
                        for res in chain
                    )
                    if has_protein:
                        asu_chains.append(chain.name)

                logger.info(f"  ASU 蛋白链: {asu_chains}")

                # 对每条链运行 Merizo
                for chain_name in asu_chains:
                    # 使用不同文件名避免 Merizo 输出目录冲突
                    chain_pdb = os.path.join(tmpdir, f"chain_{chain_name}.pdb")
                    shutil.copy2(pdb_path, chain_pdb)

                    merizo_out = self._run_merizo(
                        chain_pdb, tmpdir, pdb_chain=chain_name
                    )
                    if merizo_out is None:
                        continue

                    chain_map, _, _ = self._parse_merizo_output(
                        merizo_out, f"chain_{chain_name}", target_chain=chain_name
                    )
                    if not chain_map:
                        continue

                    # 重映射 domain IDs 为全局唯一
                    remap = {}
                    for (ch, resseq), dom_id in chain_map.items():
                        if dom_id not in remap:
                            remap[dom_id] = next_dom_id
                            next_dom_id += 1
                        asu_domain_map[(chain_name, resseq)] = remap[dom_id]

                    logger.info(f"    Chain {chain_name}: {len(remap)} domains, "
                                f"{len(chain_map)} residues")

            n_asu_domains = next_dom_id - 1

            if n_asu_domains == 0:
                logger.warning(f"  ASU 无 domain 结果")
                empty_result["method"] = "fallback_no_domains"
                with open(output_path, 'w') as f:
                    json.dump(empty_result, f, indent=2)
                return empty_result

            # Phase 2: 展开到组装体（domain ID 按 copy 递增）
            expanded_map, total_domains = self._expand_domains_to_assembly(
                model_path, asu_domain_map, n_asu_domains
            )

            # Phase 3: 序列化输出
            domains_serializable = {}
            for (chain, resseq), dom_id in expanded_map.items():
                key = f"{chain}_{resseq}"
                domains_serializable[key] = dom_id

            # 构建边界信息
            boundaries = {}
            dom_residues = {}
            for (chain, resseq), dom_id in expanded_map.items():
                if dom_id not in dom_residues:
                    dom_residues[dom_id] = {"chain": chain, "residues": []}
                dom_residues[dom_id]["residues"].append(resseq)
            for dom_id, info in dom_residues.items():
                boundaries[str(dom_id)] = {
                    "chain": info["chain"],
                    "start": min(info["residues"]),
                    "end": max(info["residues"]),
                    "n_residues": len(info["residues"]),
                }

            domain_data = {
                "n_domains": total_domains,
                "n_asu_domains": n_asu_domains,
                "domains": domains_serializable,
                "domain_boundaries": boundaries,
                "method": "merizo",
            }

            with open(output_path, 'w') as f:
                json.dump(domain_data, f, indent=2)
            logger.info(f"  结构域分割完成: ASU {n_asu_domains} domains → "
                        f"组装体 {total_domains} domains → {output_path}")
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
