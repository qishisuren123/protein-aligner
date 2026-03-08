"""
模块5: Correspondence Labeling - 体素对应标注

功能：
- 使用 gemmi 坐标系正确处理原子到体素的映射
- 原子标签、氨基酸标签、二级结构标签、链标签、置信度
"""
import os
import json
import logging
import numpy as np
import gemmi
from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)

STANDARD_AA = {
    'ALA': 1, 'ARG': 2, 'ASN': 3, 'ASP': 4, 'CYS': 5,
    'GLN': 6, 'GLU': 7, 'GLY': 8, 'HIS': 9, 'ILE': 10,
    'LEU': 11, 'LYS': 12, 'MET': 13, 'PHE': 14, 'PRO': 15,
    'SER': 16, 'THR': 17, 'TRP': 18, 'TYR': 19, 'VAL': 20
}

ATOM_TYPES = {'CA': 1, 'N': 2, 'C': 3, 'O': 4, 'CB': 5}
SS_LABELS = {'coil': 1, 'helix': 2, 'strand': 3}


class CorrespondenceLabeler:
    """体素对应标注"""

    def __init__(self, config):
        self.config = config
        self.distance_threshold = config['correspondence']['distance_threshold']

    def _parse_structure(self, cif_path):
        """解析结构文件，提取原子信息和二级结构"""
        st = gemmi.read_structure(cif_path)
        st.setup_entities()
        model = st[0]

        atoms_info = []
        chain_ids = set()

        for chain in model:
            chain_ids.add(chain.name)
            for residue in chain:
                if residue.name == "HOH":
                    continue
                for atom in residue:
                    atoms_info.append({
                        'coord': np.array([atom.pos.x, atom.pos.y, atom.pos.z]),
                        'atom_name': atom.name,
                        'atom_type': ATOM_TYPES.get(atom.name, 6),
                        'aa_type': STANDARD_AA.get(residue.name, 0),
                        'residue_name': residue.name,
                        'residue_seq': residue.seqid.num,
                        'chain': chain.name,
                    })

        # 从 mmCIF/PDB 文件中读取二级结构信息
        # gemmi 在解析文件时会自动读取 HELIX/SHEET 记录
        ss_map = {}
        try:
            # 默认 coil
            for chain in model:
                for residue in chain:
                    ss_map[(chain.name, residue.seqid.num)] = 'coil'

            # 标记 helix 区域
            for helix in st.helices:
                cname = helix.start.chain_name
                start_n = helix.start.res_id.seqid.num
                end_n = helix.end.res_id.seqid.num
                for n in range(start_n, end_n + 1):
                    ss_map[(cname, n)] = 'helix'

            # 标记 sheet/strand 区域
            for sheet in st.sheets:
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

        return atoms_info, ss_map, sorted(chain_ids)

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
        """为密度图中的每个体素生成标签"""
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

        # 解析结构
        atoms_info, ss_map, chain_ids = self._parse_structure(model_path)
        if not atoms_info:
            logger.warning(f"  无原子信息")
            return None

        logger.info(f"  原子数: {len(atoms_info)}, 链数: {len(chain_ids)}")
        chain_to_num = {c: i+1 for i, c in enumerate(chain_ids)}

        # 将原子坐标转换为 grid 索引空间
        # gemmi grid 索引: frac = cell.fractionalize(pos), idx = frac * nu/nv/nw
        all_coords_real = np.array([a['coord'] for a in atoms_info])  # 实际坐标 (Å)

        # 转换为 grid 索引坐标
        all_grid_idx = np.zeros_like(all_coords_real)
        for ai in range(len(atoms_info)):
            pos = gemmi.Position(*all_coords_real[ai])
            frac = cell.fractionalize(pos)
            all_grid_idx[ai] = [frac.x * nu, frac.y * nv, frac.z * nw]

        tree = cKDTree(all_grid_idx)

        # 初始化标签
        label_atom = np.zeros((nu, nv, nw), dtype=np.int16)
        label_aa = np.zeros((nu, nv, nw), dtype=np.int16)
        label_ss = np.zeros((nu, nv, nw), dtype=np.int16)
        label_chain = np.zeros((nu, nv, nw), dtype=np.int16)
        label_confidence = np.zeros((nu, nv, nw), dtype=np.float32)

        # 只标注有信号的体素
        threshold = 0.05
        signal_mask = data > threshold
        signal_coords = np.argwhere(signal_mask)  # [i, j, k] in grid space
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

        n_labeled = 0
        for idx_i in range(len(signal_coords)):
            i, j, k = signal_coords[idx_i]
            dist = distances[idx_i]
            atom_idx = indices[idx_i]

            if dist > dist_thresh_grid:
                continue

            atom = atoms_info[atom_idx]
            label_atom[i, j, k] = atom['atom_type']
            label_aa[i, j, k] = atom['aa_type']

            ss_key = (atom['chain'], atom['residue_seq'])
            ss_type = ss_map.get(ss_key, 'coil')
            label_ss[i, j, k] = SS_LABELS.get(ss_type, 1)

            label_chain[i, j, k] = chain_to_num.get(atom['chain'], 0)
            label_confidence[i, j, k] = np.exp(-dist**2 / (2 * sigma**2))
            n_labeled += 1

        logger.info(f"  标注体素数: {n_labeled} ({100*n_labeled/data.size:.2f}%)")

        # 保存标签
        label_files = {
            'label_atom.mrc': label_atom,
            'label_aa.mrc': label_aa,
            'label_ss.mrc': label_ss,
            'label_chain.mrc': label_chain,
            'label_confidence.mrc': label_confidence,
        }

        for fname, ldata in label_files.items():
            self._save_label_grid(ldata, grid, os.path.join(entry_dir, fname))

        stats = {
            "n_atoms_total": len(atoms_info),
            "n_voxels_labeled": int(n_labeled),
            "n_voxels_signal": int(len(signal_coords)),
            "n_chains": len(chain_ids),
            "chain_map": chain_to_num,
            "label_files": list(label_files.keys()),
            "distance_threshold_A": self.distance_threshold,
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
