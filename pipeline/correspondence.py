"""
模块5: Correspondence Labeling V3 - 体素对应标注

V3 核心变化：
- CA-only 标注范式：大多数语义标签只在 CA 原子位置标注，减少噪声
- 三组 KDTree：
  - tree_all (所有蛋白原子) → label_segment, label_atom
  - tree_ca (仅 CA 原子) → label_aa, label_ss, label_chain, label_domain, label_interface
  - tree_backbone (CA/C/N/O/CB) → label_qscore
- 8 个标签通道 + mol_map
- 蛋白专注：移除 RNA/DNA/糖基/配体标注
- Voronoi 默认关闭
- Q-score 替代高斯衰减置信度
- Domain/Interface 新标签
"""
import os
import json
import logging
import numpy as np
import gemmi
from scipy.spatial import cKDTree

from pipeline.bio_assembly import expand_to_assembly
from pipeline.interface import detect_interface

logger = logging.getLogger(__name__)

# === V3 氨基酸标签 (1-20) ===
STANDARD_AA = {
    'ALA': 1, 'ARG': 2, 'ASN': 3, 'ASP': 4, 'CYS': 5,
    'GLN': 6, 'GLU': 7, 'GLY': 8, 'HIS': 9, 'ILE': 10,
    'LEU': 11, 'LYS': 12, 'MET': 13, 'PHE': 14, 'PRO': 15,
    'SER': 16, 'THR': 17, 'TRP': 18, 'TYR': 19, 'VAL': 20
}

# === V3 原子类型标签（重新排序：CA 优先）===
ATOM_TYPES = {'CA': 1, 'C': 2, 'N': 3, 'O': 4, 'CB': 5}

# === V3 二级结构标签（重新编号）===
SS_LABELS = {'helix': 1, 'strand': 2, 'coil': 3}

# === 骨架原子集合（用于 tree_backbone）===
BACKBONE_ATOMS = {'CA', 'C', 'N', 'O', 'CB'}


class CorrespondenceLabeler:
    """V3 体素对应标注 — 三 KDTree + CA-only 范式"""

    def __init__(self, config):
        self.config = config
        corr_cfg = config['correspondence']
        self.distance_threshold = corr_cfg['distance_threshold']
        self.use_voronoi = corr_cfg.get('use_voronoi', False)  # V3 默认关闭
        # 界面检测配置
        iface_cfg = corr_cfg.get('interface', {})
        self.interface_threshold = iface_cfg.get('distance_threshold', 5.0)
        # 生物组装体配置
        bio_cfg = config.get('bio_assembly', {})
        self.bio_assembly_enabled = bio_cfg.get('enabled', True)
        self.assembly_id = bio_cfg.get('assembly_id', '1')
        self.max_chains = bio_cfg.get('max_chains', 200)

    def _build_ss_map_expanded(self, orig_structure, expanded_structure):
        """
        为展开后的结构构建 SS 映射

        策略：先从原始结构读取 SS 信息（按残基序号），
        然后对展开后的每条链，通过残基序号匹配原始链的 SS 类型
        """
        # 从原始结构获取 SS（按残基序号）
        orig_ss_by_seq = {}
        try:
            for helix in orig_structure.helices:
                start_n = helix.start.res_id.seqid.num
                end_n = helix.end.res_id.seqid.num
                for n in range(start_n, end_n + 1):
                    orig_ss_by_seq[n] = 'helix'

            for sheet in orig_structure.sheets:
                for strand in sheet.strands:
                    start_n = strand.start.res_id.seqid.num
                    end_n = strand.end.res_id.seqid.num
                    for n in range(start_n, end_n + 1):
                        orig_ss_by_seq[n] = 'strand'
        except Exception as e:
            logger.warning(f"  原始结构 SS 读取失败: {e}")

        # 对展开后的所有链应用相同的 SS 映射
        ss_map = {}
        model = expanded_structure[0]
        for chain in model:
            for residue in chain:
                seq_num = residue.seqid.num
                ss_type = orig_ss_by_seq.get(seq_num, 'coil')
                ss_map[(chain.name, seq_num)] = ss_type

        n_helix = sum(1 for v in ss_map.values() if v == 'helix')
        n_strand = sum(1 for v in ss_map.values() if v == 'strand')
        n_coil = sum(1 for v in ss_map.values() if v == 'coil')
        logger.info(f"  二级结构: {n_helix} helix, {n_strand} strand, {n_coil} coil")

        return ss_map

    def _build_ss_map(self, structure):
        """构建二级结构映射表（非展开场景）"""
        ss_map = {}
        model = structure[0]

        # 默认 coil
        for chain in model:
            for residue in chain:
                ss_map[(chain.name, residue.seqid.num)] = 'coil'

        try:
            for helix in structure.helices:
                cname = helix.start.chain_name
                start_n = helix.start.res_id.seqid.num
                end_n = helix.end.res_id.seqid.num
                for n in range(start_n, end_n + 1):
                    ss_map[(cname, n)] = 'helix'

            for sheet in structure.sheets:
                for strand in sheet.strands:
                    cname = strand.start.chain_name
                    start_n = strand.start.res_id.seqid.num
                    end_n = strand.end.res_id.seqid.num
                    for n in range(start_n, end_n + 1):
                        ss_map[(cname, n)] = 'strand'

            n_helix = sum(1 for v in ss_map.values() if v == 'helix')
            n_strand = sum(1 for v in ss_map.values() if v == 'strand')
            n_coil = sum(1 for v in ss_map.values() if v == 'coil')
            logger.info(f"  二级结构: {n_helix} helix, {n_strand} strand, {n_coil} coil")
        except Exception as e:
            logger.warning(f"  二级结构读取失败: {e}")

        return ss_map

    def _parse_structure(self, cif_path):
        """
        解析结构文件，提取三组原子列表（仅蛋白质残基）

        返回:
            all_atoms: 所有蛋白原子 → tree_all (segment, atom)
            ca_atoms: 仅 CA 原子 → tree_ca (aa, ss, chain, domain, interface)
            backbone_atoms: CA/C/N/O/CB → tree_backbone (qscore)
            ss_map: 二级结构映射
            chain_ids: 排序后的链 ID 列表
        """
        orig_st = gemmi.read_structure(cif_path)
        orig_st.setup_entities()

        # 生物组装体展开
        if self.bio_assembly_enabled:
            expanded_st = expand_to_assembly(
                orig_st,
                assembly_id=self.assembly_id,
                max_chains=self.max_chains
            )
            ss_map = self._build_ss_map_expanded(orig_st, expanded_st)
        else:
            expanded_st = orig_st.clone()
            ss_map = self._build_ss_map(expanded_st)

        model = expanded_st[0]
        all_atoms = []
        ca_atoms = []
        backbone_atoms = []
        chain_ids = set()

        for chain in model:
            chain_ids.add(chain.name)
            for residue in chain:
                if residue.name == "HOH":
                    continue

                # V3: 只处理蛋白质残基
                info = gemmi.find_tabulated_residue(residue.name)
                if info.kind != gemmi.ResidueKind.AA:
                    continue

                aa_type = STANDARD_AA.get(residue.name, 0)

                for atom in residue:
                    coord = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                    atom_type = ATOM_TYPES.get(atom.name, 6)  # 6=Other

                    atom_dict = {
                        'coord': coord,
                        'atom_name': atom.name,
                        'atom_type': atom_type,
                        'aa_type': aa_type,
                        'residue_name': residue.name,
                        'residue_seq': residue.seqid.num,
                        'chain': chain.name,
                    }

                    # 所有蛋白原子 → tree_all
                    all_atoms.append(atom_dict)

                    # CA 原子 → tree_ca
                    if atom.name == 'CA':
                        ca_atoms.append(atom_dict)

                    # 骨架原子 → tree_backbone
                    if atom.name in BACKBONE_ATOMS:
                        backbone_atoms.append(atom_dict)

        logger.info(f"  蛋白原子: all={len(all_atoms)}, CA={len(ca_atoms)}, "
                    f"backbone={len(backbone_atoms)}, 链={len(chain_ids)}")

        return all_atoms, ca_atoms, backbone_atoms, ss_map, sorted(chain_ids)

    def _load_qscores(self, entry_dir):
        """
        加载逐原子 Q-score（Phase 1 产出）

        返回:
            qscore_data: dict, key="{chain}_{resseq}_{atomname}" → float
        """
        qscore_path = os.path.join(entry_dir, "qscores.json")
        if not os.path.exists(qscore_path):
            logger.warning(f"  qscores.json 不存在，label_qscore 将全零")
            return {}
        with open(qscore_path) as f:
            return json.load(f)

    def _load_domain_assignment(self, entry_dir):
        """
        加载结构域分割结果（Phase 2 产出）

        返回:
            domain_data: dict, key="{chain}_{resseq}" → domain_id
        """
        domain_path = os.path.join(entry_dir, "domain_assignment.json")
        if not os.path.exists(domain_path):
            logger.warning(f"  domain_assignment.json 不存在，label_domain 将全零")
            return {}
        with open(domain_path) as f:
            data = json.load(f)
        return data.get("domains", {})

    def _save_label_grid(self, data, ref_grid, output_path):
        """保存标签数据为 MRC 格式（保持坐标系）"""
        new_grid = gemmi.FloatGrid(data.shape[0], data.shape[1], data.shape[2])
        new_grid.set_unit_cell(ref_grid.unit_cell)
        new_grid.spacegroup = ref_grid.spacegroup
        arr = np.array(new_grid, copy=False)
        arr[:] = data.astype(np.float32)

        ccp4 = gemmi.Ccp4Map()
        ccp4.grid = new_grid
        ccp4.update_ccp4_header()
        ccp4.write_ccp4_map(output_path)

    def label_voxels(self, entry_dir):
        """
        V3 体素标注 — 三 KDTree + CA-only 范式

        输出 8 个标签通道:
        1. label_segment.mrc: 蛋白/背景 (0/1)
        2. label_atom.mrc: 原子类型 (CA=1,C=2,N=3,O=4,CB=5,Other=6)
        3. label_aa.mrc: 氨基酸类型 (1-20), CA-only
        4. label_ss.mrc: 二级结构 (Helix=1,Sheet=2,Coil=3), CA-only
        5. label_qscore.mrc: Q-score 置信度 (float), 骨架原子
        6. label_chain.mrc: 链 ID (1-N), CA-only
        7. label_domain.mrc: 结构域 ID (1-M), CA-only
        8. label_interface.mrc: 界面标签 (0/1), CA-only
        """
        map_path = os.path.join(entry_dir, "map_normalized.mrc")
        model_path = os.path.join(entry_dir, "model.cif")

        if not os.path.exists(map_path) or not os.path.exists(model_path):
            logger.warning(f"跳过 {entry_dir}: 缺少文件")
            return None

        # 读取密度图
        ccp4 = gemmi.read_ccp4_map(map_path)
        ccp4.setup(float('nan'))
        grid = ccp4.grid
        data = np.array(grid, copy=True)
        nu, nv, nw = grid.nu, grid.nv, grid.nw
        cell = grid.unit_cell

        logger.info(f"  Grid: ({nu},{nv},{nw}), spacing: ({grid.spacing[0]:.3f}Å)")

        # 解析结构 → 三组原子列表
        all_atoms, ca_atoms, backbone_atoms, ss_map, chain_ids = \
            self._parse_structure(model_path)

        if not all_atoms:
            logger.warning(f"  无蛋白原子信息")
            return None

        chain_to_num = {c: i + 1 for i, c in enumerate(chain_ids)}

        # 加载依赖数据
        qscore_data = self._load_qscores(entry_dir)
        domain_data = self._load_domain_assignment(entry_dir)

        # 界面检测（在真实 Å 坐标空间进行）
        interface_labels = detect_interface(ca_atoms, threshold=self.interface_threshold)

        # === 坐标转换：真实 Å → grid 索引 ===
        def coords_to_grid(atoms_list):
            coords_real = np.array([a['coord'] for a in atoms_list])
            grid_idx = np.zeros_like(coords_real)
            for i in range(len(atoms_list)):
                pos = gemmi.Position(*coords_real[i])
                frac = cell.fractionalize(pos)
                grid_idx[i] = [frac.x * nu, frac.y * nv, frac.z * nw]
            return grid_idx

        all_grid_idx = coords_to_grid(all_atoms)
        ca_grid_idx = coords_to_grid(ca_atoms) if ca_atoms else np.zeros((0, 3))
        bb_grid_idx = coords_to_grid(backbone_atoms) if backbone_atoms else np.zeros((0, 3))

        # === 建 KDTree ===
        tree_all = cKDTree(all_grid_idx)
        tree_ca = cKDTree(ca_grid_idx) if len(ca_atoms) > 0 else None
        tree_backbone = cKDTree(bb_grid_idx) if len(backbone_atoms) > 0 else None

        # === 初始化 8 个标签通道 ===
        label_segment = np.zeros((nu, nv, nw), dtype=np.int16)
        label_atom = np.zeros((nu, nv, nw), dtype=np.int16)
        label_aa = np.zeros((nu, nv, nw), dtype=np.int16)
        label_ss = np.zeros((nu, nv, nw), dtype=np.int16)
        label_qscore = np.zeros((nu, nv, nw), dtype=np.float32)
        label_chain = np.zeros((nu, nv, nw), dtype=np.int16)
        label_domain = np.zeros((nu, nv, nw), dtype=np.int16)
        label_interface = np.zeros((nu, nv, nw), dtype=np.int16)

        # === 预计算属性数组 ===
        # tree_all 属性
        atom_types_all = np.array([a['atom_type'] for a in all_atoms], dtype=np.int16)

        # tree_ca 属性
        if ca_atoms:
            aa_types_ca = np.array([a['aa_type'] for a in ca_atoms], dtype=np.int16)
            ss_types_ca = np.array([
                SS_LABELS.get(ss_map.get((a['chain'], a['residue_seq']), 'coil'), 3)
                for a in ca_atoms
            ], dtype=np.int16)
            chain_nums_ca = np.array([
                chain_to_num.get(a['chain'], 0) for a in ca_atoms
            ], dtype=np.int16)
            # Domain ID：从 domain_assignment.json 读取
            domain_ids_ca = np.array([
                domain_data.get(f"{a['chain']}_{a['residue_seq']}", 0)
                for a in ca_atoms
            ], dtype=np.int16)
            # Interface：已计算
            interface_ca = interface_labels

        # tree_backbone 属性：Q-score
        if backbone_atoms:
            qscore_values_bb = np.array([
                qscore_data.get(f"{a['chain']}_{a['residue_seq']}_{a['atom_name']}", 0.0)
                for a in backbone_atoms
            ], dtype=np.float32)

        # === 查询信号体素 ===
        threshold = 0.05
        signal_mask = data > threshold
        signal_coords = np.argwhere(signal_mask)
        logger.info(f"  有信号体素: {len(signal_coords)}/{data.size} "
                    f"({100 * len(signal_coords) / data.size:.1f}%)")

        if len(signal_coords) == 0:
            logger.warning("  无信号体素")
            return None

        # 距离阈值（grid 索引单位）
        spacing_avg = (grid.spacing[0] + grid.spacing[1] + grid.spacing[2]) / 3.0
        dist_thresh_grid = self.distance_threshold / spacing_avg

        query_coords = signal_coords.astype(np.float64)

        # === 1. tree_all → label_segment + label_atom ===
        dist_all, idx_all = tree_all.query(query_coords, k=1)
        within_all = dist_all <= dist_thresh_grid
        within_idx_all = np.where(within_all)[0]

        vi = signal_coords[within_idx_all, 0]
        vj = signal_coords[within_idx_all, 1]
        vk = signal_coords[within_idx_all, 2]
        ai = idx_all[within_idx_all]

        label_segment[vi, vj, vk] = 1
        label_atom[vi, vj, vk] = atom_types_all[ai]

        n_segment = len(within_idx_all)
        logger.info(f"  label_segment/atom: {n_segment} 体素 "
                    f"({100 * n_segment / data.size:.2f}%)")

        # Voronoi 扩展（可选，V3 默认关闭）
        if self.use_voronoi:
            beyond_all = ~within_all
            beyond_idx_all = np.where(beyond_all)[0]
            if len(beyond_idx_all) > 0:
                vi_v = signal_coords[beyond_idx_all, 0]
                vj_v = signal_coords[beyond_idx_all, 1]
                vk_v = signal_coords[beyond_idx_all, 2]
                ai_v = idx_all[beyond_idx_all]
                label_segment[vi_v, vj_v, vk_v] = 1
                label_atom[vi_v, vj_v, vk_v] = atom_types_all[ai_v]
                logger.info(f"  Voronoi 扩展: +{len(beyond_idx_all)} 体素")

        # === 2. tree_ca → label_aa + label_ss + label_chain + label_domain + label_interface ===
        if tree_ca is not None:
            dist_ca, idx_ca = tree_ca.query(query_coords, k=1)
            within_ca = dist_ca <= dist_thresh_grid
            within_idx_ca = np.where(within_ca)[0]

            vi = signal_coords[within_idx_ca, 0]
            vj = signal_coords[within_idx_ca, 1]
            vk = signal_coords[within_idx_ca, 2]
            ci = idx_ca[within_idx_ca]

            label_aa[vi, vj, vk] = aa_types_ca[ci]
            label_ss[vi, vj, vk] = ss_types_ca[ci]
            label_chain[vi, vj, vk] = chain_nums_ca[ci]
            label_domain[vi, vj, vk] = domain_ids_ca[ci]
            label_interface[vi, vj, vk] = interface_ca[ci]

            n_ca_labeled = len(within_idx_ca)
            logger.info(f"  CA-only 标签 (aa/ss/chain/domain/interface): "
                        f"{n_ca_labeled} 体素")
        else:
            n_ca_labeled = 0

        # === 3. tree_backbone → label_qscore ===
        if tree_backbone is not None:
            dist_bb, idx_bb = tree_backbone.query(query_coords, k=1)
            within_bb = dist_bb <= dist_thresh_grid
            within_idx_bb = np.where(within_bb)[0]

            vi = signal_coords[within_idx_bb, 0]
            vj = signal_coords[within_idx_bb, 1]
            vk = signal_coords[within_idx_bb, 2]
            bi = idx_bb[within_idx_bb]

            label_qscore[vi, vj, vk] = qscore_values_bb[bi]

            n_qscore = len(within_idx_bb)
            logger.info(f"  label_qscore (backbone): {n_qscore} 体素")
        else:
            n_qscore = 0

        # === 保存 8 个标签文件 ===
        label_files = {
            'label_segment.mrc': label_segment,
            'label_atom.mrc': label_atom,
            'label_aa.mrc': label_aa,
            'label_ss.mrc': label_ss,
            'label_qscore.mrc': label_qscore,
            'label_chain.mrc': label_chain,
            'label_domain.mrc': label_domain,
            'label_interface.mrc': label_interface,
        }

        for fname, ldata in label_files.items():
            self._save_label_grid(ldata, grid, os.path.join(entry_dir, fname))

        # === 统计 ===
        n_interface_voxels = int((label_interface > 0).sum())
        n_domain_voxels = int((label_domain > 0).sum())
        n_unique_domains = len(set(label_domain[label_domain > 0].flatten().tolist())) \
            if n_domain_voxels > 0 else 0

        stats = {
            "n_atoms_total": len(all_atoms),
            "n_atoms_ca": len(ca_atoms),
            "n_atoms_backbone": len(backbone_atoms),
            "n_chains": len(chain_ids),
            "n_voxels_segment": int(n_segment),
            "n_voxels_ca_labeled": int(n_ca_labeled),
            "n_voxels_qscore": int(n_qscore),
            "n_voxels_signal": int(len(signal_coords)),
            "n_voxels_total": int(data.size),
            "coverage_segment_pct": float(100 * n_segment / data.size),
            "coverage_ca_pct": float(100 * n_ca_labeled / data.size),
            "n_interface_voxels": n_interface_voxels,
            "n_domain_voxels": n_domain_voxels,
            "n_unique_domains": n_unique_domains,
            "chain_map": chain_to_num,
            "label_files": list(label_files.keys()),
            "distance_threshold_A": self.distance_threshold,
            "interface_threshold_A": self.interface_threshold,
            "use_voronoi": self.use_voronoi,
            "bio_assembly_expanded": self.bio_assembly_enabled,
            "version": "V3",
        }
        with open(os.path.join(entry_dir, "labeling_stats.json"), 'w') as f:
            json.dump(stats, f, indent=2)

        return stats

    def run(self, entry_dirs):
        """对所有条目进行体素标注"""
        results = []
        for entry_dir in entry_dirs:
            logger.info(f"体素标注 (V3): {entry_dir}")
            try:
                stats = self.label_voxels(entry_dir)
                if stats:
                    results.append({"entry_dir": entry_dir, "stats": stats})
            except Exception as e:
                logger.error(f"  标注失败: {e}", exc_info=True)

        logger.info(f"体素标注完成: {len(results)} 个条目")
        return results
