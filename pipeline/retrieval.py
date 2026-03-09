"""
模块1: Retrieval - 从 EMDB/RCSB 检索和下载冷冻电镜数据

功能：
- 根据筛选条件（分辨率、方法等）查询 EMDB 数据库
- 获取关联的 PDB ID
- 下载密度图 (.map/.mrc) 和原子模型 (.cif)
"""
import os
import json
import gzip
import shutil
import logging
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)


class Retriever:
    """从 EMDB 和 RCSB 检索并下载冷冻电镜数据"""

    def __init__(self, config):
        self.config = config
        self.emdb_api = config['retrieval']['emdb_api_base']
        self.emdb_ftp = config['retrieval']['emdb_ftp_base']
        self.rcsb_api = config['retrieval']['rcsb_api_base']
        self.rcsb_download = config['retrieval']['rcsb_download_base']
        self.raw_dir = config['paths']['raw_dir']
        # 使用 requests session 复用连接
        self.session = requests.Session()

    def search_emdb(self, resolution_min=None, resolution_max=None,
                    method=None, max_entries=None):
        """
        通过 RCSB Search API 查询同时有 EMDB map 和 PDB model 的条目
        返回 (emdb_id, pdb_id) 对的列表
        """
        if resolution_min is None:
            resolution_min = self.config['retrieval']['resolution_min']
        if resolution_max is None:
            resolution_max = self.config['retrieval']['resolution_max']
        if max_entries is None:
            max_entries = self.config['retrieval']['max_entries']

        logger.info(f"搜索 EMDB 条目: 分辨率 {resolution_min}-{resolution_max}Å, 最大 {max_entries} 条")

        # 使用 RCSB Search API v2 查询 EM 方法 + 分辨率范围的 PDB 条目
        search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        # 多请求一些来过滤（有些可能没有 EMDB 关联）
        query_rows = max_entries * 3
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "exptl.method",
                            "operator": "exact_match",
                            "value": "ELECTRON MICROSCOPY"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal",
                            "value": resolution_max
                        }
                    }
                ]
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": query_rows
                },
                "sort": [
                    {
                        "sort_by": "rcsb_entry_info.resolution_combined",
                        "direction": "asc"
                    }
                ]
            }
        }

        try:
            resp = self.session.post(search_url, json=query, timeout=30)
            if resp.status_code == 204:
                logger.warning("RCSB 搜索返回 204 No Content（无匹配条目）")
                return []
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            logger.error(f"RCSB 搜索失败: {e}")
            return []

        pdb_ids = []
        if "result_set" in data:
            for item in data["result_set"]:
                pdb_ids.append(item["identifier"].upper())

        logger.info(f"找到 {len(pdb_ids)} 个 PDB 条目")

        # 对每个 PDB 条目获取关联的 EMDB ID
        entries = []
        for pdb_id in pdb_ids:
            emdb_id = self._get_emdb_id_for_pdb(pdb_id)
            if emdb_id:
                entries.append({"pdb_id": pdb_id, "emdb_id": emdb_id})
                logger.info(f"  {pdb_id} -> {emdb_id}")
            if len(entries) >= max_entries:
                break

        return entries

    def _get_emdb_id_for_pdb(self, pdb_id):
        """通过 RCSB API 获取 PDB 条目关联的 EMDB ID"""
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        try:
            resp = self.session.get(url, timeout=15)
            resp.raise_for_status()
            data = resp.json()
            # 从 pdbx_database_related 中找 EMDB 关联
            related = data.get("pdbx_database_related", [])
            for item in related:
                if item.get("db_name") == "EMDB":
                    return item.get("db_id")
            # 备用：从 em_db_mapping 获取
            em_info = data.get("rcsb_entry_container_identifiers", {})
            emdb_ids = em_info.get("emdb_ids", [])
            if emdb_ids:
                return emdb_ids[0]
        except Exception as e:
            logger.warning(f"获取 {pdb_id} 的 EMDB ID 失败: {e}")
        return None

    def _get_resolution(self, pdb_id, emdb_num):
        """
        从 RCSB API 获取分辨率信息

        优先从 RCSB entry 获取 resolution_combined，
        备用从 EMDB API 获取
        """
        # 方法1: RCSB entry API
        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            resp = self.session.get(url, timeout=15)
            resp.raise_for_status()
            data = resp.json()
            res = data.get("rcsb_entry_info", {}).get("resolution_combined", [None])
            if res and res[0] is not None:
                return float(res[0])
        except Exception as e:
            logger.debug(f"从 RCSB 获取分辨率失败: {e}")

        # 方法2: EMDB API
        try:
            url = f"{self.emdb_api}/entry/EMD-{emdb_num}"
            resp = self.session.get(url, timeout=15)
            resp.raise_for_status()
            data = resp.json()
            res = data.get("map", {}).get("resolution", {}).get("value")
            if res is not None:
                return float(res)
        except Exception as e:
            logger.debug(f"从 EMDB 获取分辨率失败: {e}")

        return None

    def download_entry(self, entry, output_dir=None):
        """
        下载一个条目的密度图和原子模型
        entry: {"pdb_id": "XXXX", "emdb_id": "EMD-XXXX"}
        """
        if output_dir is None:
            output_dir = self.raw_dir

        pdb_id = entry["pdb_id"]
        emdb_id = entry["emdb_id"]
        # 标准化 EMDB ID 格式
        emdb_num = emdb_id.replace("EMD-", "").replace("emd-", "")
        emdb_id_std = f"EMD-{emdb_num}"

        entry_dir = os.path.join(output_dir, f"{emdb_id_std}_{pdb_id}")
        os.makedirs(entry_dir, exist_ok=True)

        # 获取分辨率
        resolution = entry.get("resolution")
        if resolution is None:
            resolution = self._get_resolution(pdb_id, emdb_num)
            if resolution is not None:
                logger.info(f"  分辨率: {resolution:.2f}Å")

        # 保存元数据（包含分辨率）
        metadata = dict(entry)
        if resolution is not None:
            metadata["resolution"] = resolution
        meta_path = os.path.join(entry_dir, "metadata.json")
        with open(meta_path, 'w') as f:
            json.dump(metadata, f, indent=2)

        # 下载密度图
        map_ok = self._download_map(emdb_num, entry_dir)
        # 下载原子模型
        model_ok = self._download_model(pdb_id, entry_dir)

        if map_ok and model_ok:
            logger.info(f"✓ {emdb_id_std}/{pdb_id} 下载完成")
            return entry_dir
        else:
            logger.warning(f"✗ {emdb_id_std}/{pdb_id} 下载不完整 (map={map_ok}, model={model_ok})")
            return None

    def _download_map(self, emdb_num, entry_dir):
        """下载 EMDB 密度图"""
        map_path = os.path.join(entry_dir, "raw_map.map.gz")
        map_final = os.path.join(entry_dir, "raw_map.map")

        if os.path.exists(map_final):
            logger.info(f"  密度图已存在，跳过: {map_final}")
            return True

        # EMDB FTP 下载路径格式
        url = f"{self.emdb_ftp}/EMD-{emdb_num}/map/emd_{emdb_num}.map.gz"
        logger.info(f"  下载密度图: {url}")

        try:
            resp = self.session.get(url, stream=True, timeout=120)
            resp.raise_for_status()
            total = int(resp.headers.get('content-length', 0))
            with open(map_path, 'wb') as f:
                with tqdm(total=total, unit='B', unit_scale=True, desc="Map") as pbar:
                    for chunk in resp.iter_content(chunk_size=8192):
                        f.write(chunk)
                        pbar.update(len(chunk))
            # 解压
            with gzip.open(map_path, 'rb') as f_in:
                with open(map_final, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(map_path)
            return True
        except Exception as e:
            logger.error(f"  密度图下载失败: {e}")
            # 清理不完整的文件
            for p in [map_path, map_final]:
                if os.path.exists(p):
                    os.remove(p)
            return False

    def _download_model(self, pdb_id, entry_dir):
        """下载 RCSB 原子模型 (mmCIF 格式)"""
        model_path = os.path.join(entry_dir, "model.cif")

        if os.path.exists(model_path):
            logger.info(f"  原子模型已存在，跳过: {model_path}")
            return True

        url = f"{self.rcsb_download}/{pdb_id.lower()}.cif"
        logger.info(f"  下载原子模型: {url}")

        try:
            resp = self.session.get(url, timeout=60)
            resp.raise_for_status()
            with open(model_path, 'wb') as f:
                f.write(resp.content)
            return True
        except Exception as e:
            logger.error(f"  原子模型下载失败: {e}")
            if os.path.exists(model_path):
                os.remove(model_path)
            return False

    def update_existing_metadata(self, entry_dirs):
        """
        回填已下载条目的分辨率信息

        对已有 metadata.json 但缺少 resolution 字段的条目，从 API 获取并写入
        """
        updated = 0
        for entry_dir in entry_dirs:
            meta_path = os.path.join(entry_dir, "metadata.json")
            if not os.path.exists(meta_path):
                continue
            with open(meta_path, 'r') as f:
                meta = json.load(f)
            if "resolution" in meta and meta["resolution"] is not None:
                continue

            pdb_id = meta.get("pdb_id", "")
            emdb_id = meta.get("emdb_id", "")
            emdb_num = emdb_id.replace("EMD-", "").replace("emd-", "")
            resolution = self._get_resolution(pdb_id, emdb_num)
            if resolution is not None:
                meta["resolution"] = resolution
                with open(meta_path, 'w') as f:
                    json.dump(meta, f, indent=2)
                logger.info(f"  更新 {os.path.basename(entry_dir)} 分辨率: {resolution:.2f}Å")
                updated += 1

        logger.info(f"回填分辨率完成: {updated} 个条目已更新")
        return updated

    def run(self, **kwargs):
        """执行完整的检索和下载流程"""
        entries = self.search_emdb(**kwargs)
        results = []
        for entry in entries:
            entry_dir = self.download_entry(entry)
            if entry_dir:
                results.append({
                    "entry": entry,
                    "dir": entry_dir
                })
        logger.info(f"成功下载 {len(results)}/{len(entries)} 个条目")
        return results
