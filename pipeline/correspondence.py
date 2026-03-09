"""
模块5: Correspondence Labeling - 体素对应标注

功能：
- 使用生物学组装体展开实现完整密度覆盖
- 双层标注策略：
  - 距离阈值内: 高置信度（高斯衰减）
  - 距离阈值外: Voronoi 低置信度（最近邻赋值，上限 0.3）
- 多分子类型标注（蛋白质、RNA、DNA、糖基、配体）
- 核苷酸/糖基标签扩展
- RNA/DNA 二级结构区分（sugar_phosphate/base）
- 向量化赋值提升性能
- 二级结构注释复制到所有扩展链
"""
import os
import json
import logging
import numpy as np
import gemmi
from scipy.spatial import cKDTree

from pipeline.bio_assembly import expand_to_assembly

logger = logging.getLogger(__name__)

# === 氨基酸标签 (1-20) ===
STANDARD_AA = {
    'ALA': 1, 'ARG': 2, 'ASN': 3, 'ASP': 4, 'CYS': 5,
    'GLN': 6, 'GLU': 7, 'GLY': 8, 'HIS': 9, 'ILE': 10,
    'LEU': 11, 'LYS': 12, 'MET': 13, 'PHE': 14, 'PRO': 15,
    'SER': 16, 'THR': 17, 'TRP': 18, 'TYR': 19, 'VAL': 20
}

# === 核苷酸标签 (21-28) ===
NUCLEOTIDE_LABELS = {
    'A': 21, 'U': 22, 'G': 23, 'C': 24,     # RNA
    'DA': 25, 'DT': 26, 'DG': 27, 'DC': 28,  # DNA
}

# === 糖基标签 (29-34) ===
SUGAR_LABELS = {
    'NAG': 29, 'MAN': 30, 'BMA': 31, 'FUC': 32, 'GAL': 33, 'SIA': 34
}

# === 统一残基/单体标签字典 ===
RESIDUE_LABELS = {}
RESIDUE_LABELS.update(STANDARD_AA)
RESIDUE_LABELS.update(NUCLEOTIDE_LABELS)
RESIDUE_LABELS.update(SUGAR_LABELS)

# === 原子类型标签 ===
ATOM_TYPES = {'CA': 1, 'N': 2, 'C': 3, 'O': 4, 'CB': 5}

# === 二级结构标签（扩展） ===
SS_LABELS = {
    'coil': 1, 'helix': 2, 'strand': 3,
    'sugar_phosphate': 4, 'base': 5  # RNA/DNA 专用
}

# === 分子类型标签 ===
MOLTYPE_LABELS = {
    'protein': 1, 'RNA': 2, 'DNA': 3, 'sugar': 4, 'ligand': 5
}

# 常见糖基残基集合（用于分子类型检测）
DEFAULT_SUGAR_RESIDUES = {'NAG', 'MAN', 'BMA', 'FUC', 'GAL', 'SIA',
                          'BGC', 'GLC', 'NDG', 'XYS'}

# RNA/DNA 骨架原子名（sugar-phosphate backbone）
_BACKBONE_ATOMS = {
    "P", "OP1", "OP2", "OP3",
    "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
    # 部分 PDB 用无引号命名
    "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*",
}


def _detect_moltype(residue, sugar_residues=None):
    """
    检测残基的分子类型

    使用 gemmi.find_tabulated_residue() 获取残基信息，
    通过 info.kind (ResidueKind) 判断分子类型。

    返回:
        moltype: 'protein', 'RNA', 'DNA', 'sugar', 'ligand', 或 None（水）
    """
    if sugar_residues is None:
        sugar_residues = DEFAULT_SUGAR_RESIDUES

    # 水分子不标注
    if residue.name == 'HOH':
        return None

    # 糖基残基优先检测（常见修饰糖）
    if residue.name in sugar_residues:
        return 'sugar'

    info = gemmi.find_tabulated_residue(residue.name)
    kind = info.kind

    # 蛋白质
    if kind == gemmi.ResidueKind.AA:
        return 'protein'
    # RNA
    if kind == gemmi.ResidueKind.RNA:
        return 'RNA'
    # DNA
    if kind == gemmi.ResidueKind.DNA:
        return 'DNA'
    # 糖（pyranose，gemmi 内置分类）
    if kind == gemmi.ResidueKind.PYR:
        return 'sugar'

    # 其他非水小分子（ELS=配体, BUF=缓冲液 等）
    return 'ligand'


def _get_nucleotide_ss(atom_name):
    """
    RNA/DNA 原子的二级结构分类

    返回 'sugar_phosphate' 或 'base'
    """
    if atom_name in _BACKBONE_ATOMS:
        return 'sugar_phosphate'
    return 'base'


class CorrespondenceLabeler:
    """体素对应标注"""

    def __init__(self, config):
        self.config = config
        corr_cfg = config['correspondence']
        self.distance_threshold = corr_cfg['distance_threshold']
        self.use_voronoi = corr_cfg.get('use_voronoi', True)
        self.voronoi_max_confidence = corr_cfg.get('voronoi_max_confidence', 0.3)
        # 多分子类型配置
        mol_cfg = corr_cfg.get('molecule_types', {})
        self.moltype_enabled = mol_cfg.get('enabled', True)
        self.sugar_residues = set(mol_cfg.get('sugar_residues',
                                              list(DEFAULT_SUGAR_RESIDUES)))
        # 生物组装体配置
        bio_cfg = config.get('bio_assembly', {})
        self.bio_assembly_enabled = bio_cfg.get('enabled', True)
        self.assembly_id = bio_cfg.get('assembly_id', '1')
        self.max_chains = bio_cfg.get('max_chains', 200)

    def _build_ss_map(self, structure):
        """
        构建二级结构映射表

        从原始结构的 helices/sheets 信息构建 (chain_name, residue_seq) -> ss_type 映射。
        对于展开后的链，通过残基序号匹配原始链的 SS 注释。
        """
        ss_map = {}
        model = structure[0]

        # 默认 coil
        for chain in model:
            for residue in chain:
                ss_map[(chain.name, residue.seqid.num)] = 'coil'

        try:
            # 标记 helix 区域
            for helix in structure.helices:
                cname = helix.start.chain_name
                start_n = helix.start.res_id.seqid.num
                end_n = helix.end.res_id.seqid.num
                for n in range(start_n, end_n + 1):
                    ss_map[(cname, n)] = 'helix'

            # 标记 sheet/strand 区域
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

    def _build_ss_map_expanded(self, orig_structure, expanded_structure):
        """
        为展开后的结构构建 SS 映射

        策略：先从原始结构读取 SS 信息，然后对展开后的每条链，
        通过残基序号匹配对应的 SS 类型
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

        return ss_map

    def _parse_structure(self, cif_path):
        """解析结构文件，提取原子信息、二级结构和分子类型"""
        orig_st = gemmi.read_structure(cif_path)
        orig_st.setup_entities()

        # 生物组装体展开
        if self.bio_assembly_enabled:
            expanded_st = expand_to_assembly(
                orig_st,
                assembly_id=self.assembly_id,
                max_chains=self.max_chains
            )
            # 展开后的 SS 映射
            ss_map = self._build_ss_map_expanded(orig_st, expanded_st)
        else:
            expanded_st = orig_st.clone()
            ss_map = self._build_ss_map(expanded_st)

        model = expanded_st[0]
        atoms_info = []
        chain_ids = set()
        # 分子类型统计
        moltype_counts = {'protein': 0, 'RNA': 0, 'DNA': 0, 'sugar': 0, 'ligand': 0}

        for chain in model:
            chain_ids.add(chain.name)
            for residue in chain:
                if residue.name == "HOH":
                    continue

                # 检测分子类型
                moltype = _detect_moltype(residue, self.sugar_residues)
                if moltype is None:
                    continue

                # 残基标签：统一查 RESIDUE_LABELS
                aa_type = RESIDUE_LABELS.get(residue.name, 0)

                for atom in residue:
                    # RNA/DNA 的二级结构特殊处理
                    if moltype in ('RNA', 'DNA'):
                        atom_ss = _get_nucleotide_ss(atom.name)
                    else:
                        atom_ss = None  # 用 ss_map

                    atoms_info.append({
                        'coord': np.array([atom.pos.x, atom.pos.y, atom.pos.z]),
                        'atom_name': atom.name,
                        'atom_type': ATOM_TYPES.get(atom.name, 6),
                        'aa_type': aa_type,
                        'residue_name': residue.name,
                        'residue_seq': residue.seqid.num,
                        'chain': chain.name,
                        'moltype': moltype,
                        'nucleotide_ss': atom_ss,
                    })

                moltype_counts[moltype] += 1

        # 日志输出分子类型统计（残基级）
        nonzero = {k: v for k, v in moltype_counts.items() if v > 0}
        logger.info(f"  分子类型（残基数）: {nonzero}")

        return atoms_info, ss_map, sorted(chain_ids), moltype_counts

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
        为密度图中的每个体素生成标签

        双层策略：
        1. 距离阈值内: 按最近原子赋标签，高置信度（高斯衰减）
        2. 距离阈值外（Voronoi）: 按最近原子赋标签，低置信度（上限 voronoi_max_confidence）

        标签通道：
        - label_atom.mrc: 原子类型
        - label_aa.mrc: 氨基酸/核苷酸/糖基类型 (1-34)
        - label_ss.mrc: 二级结构 (1-5，含 RNA/DNA sugar_phosphate/base)
        - label_chain.mrc: 链 ID
        - label_confidence.mrc: 置信度
        - label_moltype.mrc: 分子类型 (protein=1, RNA=2, DNA=3, sugar=4, ligand=5)
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

        # 解析结构（包含生物组装体展开 + 分子类型检测）
        atoms_info, ss_map, chain_ids, moltype_counts = self._parse_structure(model_path)
        if not atoms_info:
            logger.warning(f"  无原子信息")
            return None

        logger.info(f"  原子数: {len(atoms_info)}, 链数: {len(chain_ids)}")
        chain_to_num = {c: i+1 for i, c in enumerate(chain_ids)}

        # 将原子坐标转换为 grid 索引空间（向量化）
        all_coords_real = np.array([a['coord'] for a in atoms_info])

        # 向量化坐标转换
        all_grid_idx = np.zeros_like(all_coords_real)
        for ai in range(len(atoms_info)):
            pos = gemmi.Position(*all_coords_real[ai])
            frac = cell.fractionalize(pos)
            all_grid_idx[ai] = [frac.x * nu, frac.y * nv, frac.z * nw]

        tree = cKDTree(all_grid_idx)

        # 初始化标签（6 个通道）
        label_atom = np.zeros((nu, nv, nw), dtype=np.int16)
        label_aa = np.zeros((nu, nv, nw), dtype=np.int16)
        label_ss = np.zeros((nu, nv, nw), dtype=np.int16)
        label_chain = np.zeros((nu, nv, nw), dtype=np.int16)
        label_confidence = np.zeros((nu, nv, nw), dtype=np.float32)
        label_moltype = np.zeros((nu, nv, nw), dtype=np.int16)

        # 预计算原子属性数组（向量化赋值用）
        atom_types = np.array([a['atom_type'] for a in atoms_info], dtype=np.int16)
        aa_types = np.array([a['aa_type'] for a in atoms_info], dtype=np.int16)

        # 二级结构：RNA/DNA 用 nucleotide_ss，蛋白质用 ss_map
        ss_types = np.array([
            SS_LABELS.get(a['nucleotide_ss'], 0) if a['nucleotide_ss'] is not None
            else SS_LABELS.get(ss_map.get((a['chain'], a['residue_seq']), 'coil'), 1)
            for a in atoms_info
        ], dtype=np.int16)

        chain_nums = np.array([
            chain_to_num.get(a['chain'], 0) for a in atoms_info
        ], dtype=np.int16)

        # 分子类型数组
        moltype_arr = np.array([
            MOLTYPE_LABELS.get(a['moltype'], 0) for a in atoms_info
        ], dtype=np.int16)

        # 只标注有信号的体素
        threshold = 0.05
        signal_mask = data > threshold
        signal_coords = np.argwhere(signal_mask)
        logger.info(f"  有信号体素: {len(signal_coords)}/{data.size} "
                    f"({100*len(signal_coords)/data.size:.1f}%)")

        if len(signal_coords) == 0:
            logger.warning("  无信号体素")
            return None

        # 距离阈值（grid 索引单位）
        spacing_avg = (grid.spacing[0] + grid.spacing[1] + grid.spacing[2]) / 3.0
        dist_thresh_grid = self.distance_threshold / spacing_avg
        sigma = dist_thresh_grid / 2.0

        # 批量最近邻查询
        query_coords = signal_coords.astype(np.float64)
        distances, indices = tree.query(query_coords, k=1)

        # === 第一层：距离阈值内 - 高置信度 ===
        within_mask = distances <= dist_thresh_grid
        within_idx = np.where(within_mask)[0]

        # 向量化赋值（阈值内）
        voxel_i = signal_coords[within_idx, 0]
        voxel_j = signal_coords[within_idx, 1]
        voxel_k = signal_coords[within_idx, 2]
        atom_idx = indices[within_idx]

        label_atom[voxel_i, voxel_j, voxel_k] = atom_types[atom_idx]
        label_aa[voxel_i, voxel_j, voxel_k] = aa_types[atom_idx]
        label_ss[voxel_i, voxel_j, voxel_k] = ss_types[atom_idx]
        label_chain[voxel_i, voxel_j, voxel_k] = chain_nums[atom_idx]
        label_moltype[voxel_i, voxel_j, voxel_k] = moltype_arr[atom_idx]
        label_confidence[voxel_i, voxel_j, voxel_k] = np.exp(
            -distances[within_idx]**2 / (2 * sigma**2)
        ).astype(np.float32)

        n_hard = len(within_idx)
        logger.info(f"  高置信度标注: {n_hard} 体素 ({100*n_hard/data.size:.2f}%)")

        # === 第二层：Voronoi - 低置信度 ===
        n_voronoi = 0
        if self.use_voronoi:
            beyond_mask = ~within_mask
            beyond_idx = np.where(beyond_mask)[0]

            if len(beyond_idx) > 0:
                voxel_i = signal_coords[beyond_idx, 0]
                voxel_j = signal_coords[beyond_idx, 1]
                voxel_k = signal_coords[beyond_idx, 2]
                atom_idx = indices[beyond_idx]

                label_atom[voxel_i, voxel_j, voxel_k] = atom_types[atom_idx]
                label_aa[voxel_i, voxel_j, voxel_k] = aa_types[atom_idx]
                label_ss[voxel_i, voxel_j, voxel_k] = ss_types[atom_idx]
                label_chain[voxel_i, voxel_j, voxel_k] = chain_nums[atom_idx]
                label_moltype[voxel_i, voxel_j, voxel_k] = moltype_arr[atom_idx]

                # Voronoi 置信度：随距离衰减，上限为 voronoi_max_confidence
                voronoi_conf = self.voronoi_max_confidence * np.exp(
                    -distances[beyond_idx]**2 / (2 * (dist_thresh_grid * 3)**2)
                )
                label_confidence[voxel_i, voxel_j, voxel_k] = voronoi_conf.astype(np.float32)
                n_voronoi = len(beyond_idx)

            logger.info(f"  Voronoi 标注: {n_voronoi} 体素 ({100*n_voronoi/data.size:.2f}%)")

        total_labeled = n_hard + n_voronoi
        logger.info(f"  总标注体素: {total_labeled} ({100*total_labeled/data.size:.2f}%)")

        # 保存标签
        label_files = {
            'label_atom.mrc': label_atom,
            'label_aa.mrc': label_aa,
            'label_ss.mrc': label_ss,
            'label_chain.mrc': label_chain,
            'label_confidence.mrc': label_confidence,
            'label_moltype.mrc': label_moltype,
        }

        for fname, ldata in label_files.items():
            self._save_label_grid(ldata, grid, os.path.join(entry_dir, fname))

        # 分子类型体素统计
        moltype_voxel_counts = {}
        for name, val in MOLTYPE_LABELS.items():
            count = int((label_moltype == val).sum())
            if count > 0:
                moltype_voxel_counts[name] = count

        stats = {
            "n_atoms_total": len(atoms_info),
            "n_chains": len(chain_ids),
            "n_voxels_hard_labeled": int(n_hard),
            "n_voxels_voronoi_labeled": int(n_voronoi),
            "n_voxels_total_labeled": int(total_labeled),
            "n_voxels_signal": int(len(signal_coords)),
            "n_voxels_total": int(data.size),
            "coverage_hard_pct": float(100 * n_hard / data.size),
            "coverage_total_pct": float(100 * total_labeled / data.size),
            "coverage_of_signal_pct": float(100 * total_labeled / max(len(signal_coords), 1)),
            "chain_map": chain_to_num,
            "label_files": list(label_files.keys()),
            "distance_threshold_A": self.distance_threshold,
            "use_voronoi": self.use_voronoi,
            "voronoi_max_confidence": self.voronoi_max_confidence,
            "bio_assembly_expanded": self.bio_assembly_enabled,
            "moltype_residue_counts": {k: v for k, v in moltype_counts.items() if v > 0},
            "moltype_voxel_counts": moltype_voxel_counts,
        }
        with open(os.path.join(entry_dir, "labeling_stats.json"), 'w') as f:
            json.dump(stats, f, indent=2)

        return stats

    def run(self, entry_dirs):
        """对所有条目进行体素标注"""
        results = []
        for entry_dir in entry_dirs:
            logger.info(f"体素标注: {entry_dir}")
            try:
                stats = self.label_voxels(entry_dir)
                if stats:
                    results.append({"entry_dir": entry_dir, "stats": stats})
            except Exception as e:
                logger.error(f"  标注失败: {e}", exc_info=True)

        logger.info(f"体素标注完成: {len(results)} 个条目")
        return results
