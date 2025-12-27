#!/usr/bin/env python3
"""
Single-Cell Database Metadata Collector - 完整版
采集主流单细胞数据库的元数据，用于规模-引用-用途分析

作者: AI Assistant
版本: 1.0.0
日期: 2025-11-25

支持的数据库:
1. Broad SCP
2. EMBL-EBI SCXA
3. CZI CELLxGENE
4. CNGBdb
5. Heidelberg Reuse-Trends
6. NCBI SRA
7. NCBI GEO
8. HCA Data Portal
9. Figshare
10. Dryad
11. HTAN
12. PsychAN

使用示例:
    # 采集所有数据库
    python sc_database_collector.py --all --output ./output --parallel --workers 8
    
    # 采集特定数据库
    python sc_database_collector.py --databases cellxgene hca --output ./output
    
    # 仅合并已有数据
    python sc_database_collector.py --merge-only --input ./raw_data --output ./output
"""

import argparse
import logging
import sys
import json
import csv
import hashlib
import time
import re
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Any
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlencode, quote

try:
    import requests
except ImportError:
    print("错误: 需要安装 requests 库")
    print("请运行: pip install requests")
    sys.exit(1)


# ============================================================================
# 工具函数模块
# ============================================================================

def setup_logging(level='INFO', log_file=None):
    """
    设置日志配置
    
    Args:
        level: 日志级别 (DEBUG, INFO, WARNING, ERROR)
        log_file: 日志文件路径
    """
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    handlers = [logging.StreamHandler(sys.stdout)]
    
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding='utf-8'))
    
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format=log_format,
        datefmt=date_format,
        handlers=handlers,
        force=True
    )


def validate_output_dir(output_dir: Path):
    """验证并创建输出目录"""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / 'raw').mkdir(exist_ok=True)


def generate_hash_key(dataset: Dict, keys: List[str]) -> str:
    """
    生成数据集哈希键用于去重
    
    Args:
        dataset: 数据集字典
        keys: 用于生成哈希的字段列表
    
    Returns:
        str: MD5哈希键
    """
    values = []
    for key in keys:
        value = dataset.get(key)
        if value is not None:
            # 转换为字符串并标准化
            if isinstance(value, (int, float)):
                values.append(str(value))
            elif isinstance(value, str):
                values.append(value.strip().lower())
            else:
                values.append(str(value))
        else:
            values.append('')
    
    combined = '|'.join(values)
    return hashlib.md5(combined.encode('utf-8')).hexdigest()


def deduplicate_datasets(datasets: List[Dict], keys: List[str]) -> List[Dict]:
    """
    去重数据集
    
    Args:
        datasets: 数据集列表
        keys: 去重键字段列表
    
    Returns:
        List[Dict]: 去重后的数据集列表
    """
    logger = logging.getLogger(__name__)
    seen_hashes = {}
    deduplicated = []
    duplicate_count = 0
    
    for idx, dataset in enumerate(datasets):
        hash_key = generate_hash_key(dataset, keys)
        
        if hash_key not in seen_hashes:
            seen_hashes[hash_key] = idx
            deduplicated.append(dataset)
        else:
            duplicate_count += 1
            # 合并引用计数等信息
            original_idx = seen_hashes[hash_key]
            original = deduplicated[original_idx - duplicate_count]
            
            # 如果新记录有引用计数而旧记录没有，更新之
            if not original.get('citation_count') and dataset.get('citation_count'):
                original['citation_count'] = dataset['citation_count']
    
    logger.info(f"去重完成: 原始 {len(datasets)} 条, 去重后 {len(deduplicated)} 条, 重复 {duplicate_count} 条")
    
    return deduplicated


def export_to_csv(datasets: List[Dict], output_file: Path):
    """
    导出数据集到CSV
    
    Args:
        datasets: 数据集列表
        output_file: 输出CSV文件路径
    """
    logger = logging.getLogger(__name__)
    
    if not datasets:
        logger.warning("没有数据可导出")
        return
    
    # 确定所有字段
    all_fields = set()
    for dataset in datasets:
        all_fields.update(dataset.keys())
    
    # 移除raw_metadata字段(太大，不适合CSV)
    all_fields.discard('raw_metadata')
    
    # 标准字段顺序
    standard_fields = [
        'database', 'dataset_id', 'title', 'n_cells', 'data_size_gb',
        'doi', 'accession', 'lineage_annotation_rate', 'citation_count',
        'tissue', 'disease', 'organism', 'technology', 'study_type',
        'keywords', 'publication_date', 'last_update', 'authors',
        'url', 'download_url'
    ]
    
    # 标准字段在前，其他字段在后
    ordered_fields = []
    for field in standard_fields:
        if field in all_fields:
            ordered_fields.append(field)
            all_fields.remove(field)
    
    # 添加剩余字段
    ordered_fields.extend(sorted(all_fields))
    
    # 写入CSV
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        with open(output_file, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.DictWriter(f, fieldnames=ordered_fields, extrasaction='ignore')
            writer.writeheader()
            
            for dataset in datasets:
                # 清理数据
                clean_row = {}
                for field in ordered_fields:
                    value = dataset.get(field, '')
                    # 处理列表和字典
                    if isinstance(value, (list, dict)):
                        value = json.dumps(value, ensure_ascii=False)
                    clean_row[field] = value
                
                writer.writerow(clean_row)
        
        logger.info(f"CSV已导出: {output_file}, {len(datasets)} 条记录, {len(ordered_fields)} 个字段")
        
    except Exception as e:
        logger.error(f"导出CSV失败: {e}", exc_info=True)
        raise


def safe_get(data: Dict, *keys, default=None) -> Any:
    """
    安全获取嵌套字典的值
    
    Args:
        data: 字典数据
        *keys: 嵌套键路径
        default: 默认值
    
    Returns:
        获取到的值或默认值
    """
    for key in keys:
        if isinstance(data, dict):
            data = data.get(key, {})
        else:
            return default
    return data if data != {} else default


# ============================================================================
# 基础收集器抽象类
# ============================================================================

class BaseCollector(ABC):
    """数据库收集器基类"""
    
    DATABASE_NAME = "base"  # 子类需要重写
    
    # 标准字段映射
    STANDARD_FIELDS = [
        'database',           # 数据库名称
        'dataset_id',         # 数据集ID
        'title',              # 标题
        'n_cells',            # 细胞数
        'data_size_gb',       # 数据大小(GB)
        'doi',                # DOI
        'accession',          # 登录号
        'lineage_annotation_rate',  # 血统注释率
        'citation_count',     # 引用量
        'tissue',             # 组织
        'disease',            # 疾病
        'organism',           # 物种
        'technology',         # 测序技术
        'study_type',         # 研究类型
        'keywords',           # 关键词
        'publication_date',   # 发布日期
        'last_update',        # 最后更新
        'authors',            # 作者
        'url',                # 数据集URL
        'download_url',       # 下载URL
        'raw_metadata'        # 原始元数据(JSON)
    ]
    
    def __init__(self, max_records=None, timeout=300, retry=3, 
                 skip_citation=False, use_cache=False):
        """
        初始化收集器
        
        Args:
            max_records: 最大采集记录数
            timeout: 超时时间(秒)
            retry: 重试次数
            skip_citation: 是否跳过引用量查询
            use_cache: 是否使用缓存
        """
        self.max_records = max_records
        self.timeout = timeout
        self.retry = retry
        self.skip_citation = skip_citation
        self.use_cache = use_cache
        self.logger = logging.getLogger(f"collector.{self.DATABASE_NAME}")
        
        # 创建HTTP会话
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'SingleCellMetadataCollector/1.0 (Research Tool; mailto:research@example.com)'
        })
    
    @abstractmethod
    def fetch_datasets(self) -> List[Dict]:
        """
        获取数据集列表(需子类实现)
        
        Returns:
            List[Dict]: 原始数据集列表
        """
        pass
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """
        规范化数据集字段(子类可重写以自定义映射)
        
        Args:
            raw_data: 原始数据
        
        Returns:
            Dict: 规范化后的数据
        """
        # 默认实现，子类应重写此方法
        return {
            'database': self.DATABASE_NAME,
            'raw_metadata': raw_data
        }
    
    def fetch_citation_count(self, doi: str) -> Optional[int]:
        """
        查询DOI引用量
        
        Args:
            doi: DOI标识符
        
        Returns:
            Optional[int]: 引用量
        """
        if self.skip_citation or not doi:
            return None
        
        # 清理DOI
        doi = doi.strip()
        if doi.startswith('http'):
            doi = doi.split('doi.org/')[-1]
        
        try:
            # 使用Crossref API查询引用量
            url = f"https://api.crossref.org/works/{quote(doi)}"
            response = self.session.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                count = safe_get(data, 'message', 'is-referenced-by-count', default=0)
                return int(count) if count else 0
            elif response.status_code == 404:
                self.logger.debug(f"DOI未找到: {doi}")
                return 0
        except Exception as e:
            self.logger.warning(f"查询DOI {doi} 引用量失败: {e}")
        
        return None
    
    def safe_request(self, url: str, method='GET', **kwargs) -> Optional[requests.Response]:
        """
        安全的HTTP请求，带重试机制
        
        Args:
            url: 请求URL
            method: 请求方法
            **kwargs: 其他请求参数
        
        Returns:
            Optional[requests.Response]: 响应对象
        """
        for attempt in range(self.retry):
            try:
                if method.upper() == 'GET':
                    response = self.session.get(url, timeout=self.timeout, **kwargs)
                elif method.upper() == 'POST':
                    response = self.session.post(url, timeout=self.timeout, **kwargs)
                else:
                    raise ValueError(f"不支持的HTTP方法: {method}")
                
                response.raise_for_status()
                return response
                
            except requests.RequestException as e:
                self.logger.warning(f"请求失败 (尝试 {attempt + 1}/{self.retry}): {e}")
                if attempt < self.retry - 1:
                    time.sleep(min(2 ** attempt, 10))  # 指数退避，最大10秒
                else:
                    self.logger.error(f"请求最终失败: {url}")
                    return None
        
        return None
    
    def collect(self) -> List[Dict]:
        """
        执行采集流程
        
        Returns:
            List[Dict]: 规范化后的数据集列表
        """
        self.logger.info(f"开始从 {self.DATABASE_NAME} 采集数据")
        start_time = time.time()
        
        try:
            # 获取原始数据
            raw_datasets = self.fetch_datasets()
            self.logger.info(f"获取到 {len(raw_datasets)} 条原始记录")
            
            # 限制记录数
            if self.max_records and len(raw_datasets) > self.max_records:
                raw_datasets = raw_datasets[:self.max_records]
                self.logger.info(f"限制为 {len(raw_datasets)} 条记录")
            
            # 规范化数据
            normalized_datasets = []
            for i, raw_data in enumerate(raw_datasets, 1):
                try:
                    normalized = self.normalize_dataset(raw_data)
                    
                    # 查询引用量
                    if not self.skip_citation and normalized.get('doi'):
                        citation_count = self.fetch_citation_count(normalized['doi'])
                        if citation_count is not None:
                            normalized['citation_count'] = citation_count
                    
                    normalized_datasets.append(normalized)
                    
                    if i % 100 == 0:
                        self.logger.info(f"已处理 {i}/{len(raw_datasets)} 条记录")
                        
                except Exception as e:
                    self.logger.error(f"规范化第 {i} 条记录失败: {e}")
                    continue
            
            elapsed = time.time() - start_time
            self.logger.info(f"采集完成: {len(normalized_datasets)} 条记录, 耗时 {elapsed:.2f} 秒")
            
            return normalized_datasets
            
        except Exception as e:
            self.logger.error(f"采集过程发生错误: {e}", exc_info=True)
            raise


# ============================================================================
# 具体数据库收集器实现
# ============================================================================

class CELLxGENECollector(BaseCollector):
    """CZI CELLxGENE 数据库收集器"""
    
    DATABASE_NAME = "cellxgene"
    API_BASE_URL = "https://api.cellxgene.cziscience.com"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取CELLxGENE数据集"""
        datasets = []
        
        # 获取所有集合
        collections_url = f"{self.API_BASE_URL}/curation/v1/collections"
        response = self.safe_request(collections_url)
        
        if not response:
            self.logger.error("无法获取collections列表")
            return datasets
        
        collections = response.json()
        self.logger.info(f"获取到 {len(collections)} 个collections")
        
        # 遍历每个collection
        for i, collection in enumerate(collections, 1):
            collection_id = collection.get('collection_id')
            
            # 获取collection详情
            detail_url = f"{self.API_BASE_URL}/curation/v1/collections/{collection_id}"
            detail_response = self.safe_request(detail_url)
            
            if detail_response:
                detail_data = detail_response.json()
                datasets.append(detail_data)
                
                if i % 50 == 0:
                    self.logger.info(f"已获取 {i}/{len(collections)} 个collections详情")
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化CELLxGENE数据"""
        # 计算总细胞数和数据大小
        n_cells = 0
        total_size_bytes = 0
        
        for dataset in raw_data.get('datasets', []):
            n_cells += dataset.get('cell_count', 0)
            # 获取数据集资源大小
            for asset in dataset.get('dataset_assets', []):
                total_size_bytes += asset.get('filesize', 0)
        
        # 提取组织和疾病
        tissues = set()
        diseases = set()
        technologies = set()
        organisms = set()
        
        for dataset in raw_data.get('datasets', []):
            # 组织
            for tissue in dataset.get('tissue', []):
                if tissue.get('label'):
                    tissues.add(tissue['label'])
            
            # 疾病
            for disease in dataset.get('disease', []):
                if disease.get('label'):
                    diseases.add(disease['label'])
            
            # 技术
            for assay in dataset.get('assay', []):
                if assay.get('label'):
                    technologies.add(assay['label'])
            
            # 物种
            organism = safe_get(dataset, 'organism', 'label')
            if organism:
                organisms.add(organism)
        
        # 提取作者
        authors = []
        for contact in raw_data.get('contact_name', '').split(','):
            contact = contact.strip()
            if contact:
                authors.append(contact)
        
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': raw_data.get('collection_id'),
            'title': raw_data.get('name'),
            'n_cells': n_cells,
            'data_size_gb': round(total_size_bytes / (1024**3), 2) if total_size_bytes > 0 else None,
            'doi': raw_data.get('doi'),
            'tissue': ', '.join(sorted(filter(None, tissues))),
            'disease': ', '.join(sorted(filter(None, diseases))),
            'organism': ', '.join(sorted(organisms)),
            'technology': ', '.join(sorted(technologies)),
            'publication_date': raw_data.get('published_at'),
            'last_update': raw_data.get('revised_at'),
            'authors': ', '.join(authors) if authors else None,
            'url': raw_data.get('collection_url'),
            'raw_metadata': raw_data
        }
        
        return normalized


class HCACollector(BaseCollector):
    """Human Cell Atlas (HCA) 数据库收集器"""
    
    DATABASE_NAME = "hca"
    API_BASE_URL = "https://service.azul.data.humancellatlas.org"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取HCA数据集"""
        datasets = []
        
        # 搜索所有项目
        search_url = f"{self.API_BASE_URL}/index/projects"
        
        params = {
            'catalog': 'dcp2',
            'size': 100,
            'sort': 'projectTitle',
            'order': 'asc'
        }
        
        while True:
            response = self.safe_request(search_url, params=params)
            
            if not response:
                break
            
            data = response.json()
            hits = data.get('hits', [])
            
            if not hits:
                break
            
            for hit in hits:
                project = hit.get('projects', [{}])[0]
                datasets.append(project)
            
            self.logger.info(f"已获取 {len(datasets)} 个项目")
            
            # 检查是否有下一页
            pagination = data.get('pagination', {})
            if pagination.get('next'):
                params['search_after'] = pagination.get('search_after')
                params['search_after_uid'] = pagination.get('search_after_uid')
            else:
                break
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化HCA数据"""
        # 提取贡献者
        contributors = raw_data.get('contributors', [])
        authors = []
        for contrib in contributors[:5]:  # 限制前5个作者
            name = contrib.get('contactName', '')
            if name:
                authors.append(name)
        
        # 提取组织和疾病
        tissues = set()
        diseases = set()
        
        for specimen in raw_data.get('specimens', []):
            organ = safe_get(specimen, 'organ', 'text')
            if organ:
                tissues.add(organ)
            
            disease_info = safe_get(specimen, 'disease', 'text')
            if disease_info:
                diseases.add(disease_info)
        
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': raw_data.get('projectId'),
            'title': raw_data.get('projectTitle'),
            'n_cells': raw_data.get('estimatedCellCount'),
            'accession': raw_data.get('accessions', [{}])[0].get('accession') if raw_data.get('accessions') else None,
            'tissue': ', '.join(sorted(tissues)),
            'disease': ', '.join(sorted(diseases)),
            'organism': ', '.join([org.get('text', '') for org in raw_data.get('genusSpecies', [])]),
            'publication_date': raw_data.get('publicationDate'),
            'last_update': raw_data.get('lastModifiedDate'),
            'authors': ', '.join(authors) if authors else None,
            'url': f"https://data.humancellatlas.org/explore/projects/{raw_data.get('projectId')}",
            'raw_metadata': raw_data
        }
        
        return normalized


class NCBIGEOCollector(BaseCollector):
    """NCBI GEO 数据库收集器"""
    
    DATABASE_NAME = "ncbi_geo"
    API_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取NCBI GEO数据集"""
        datasets = []
        
        # 搜索单细胞相关数据集
        search_terms = [
            '"single cell"[Title/Abstract] AND "RNA-seq"[DataSet Type]',
            '"scRNA-seq"[Title/Abstract]',
            '"single cell transcriptome"[Title/Abstract]'
        ]
        
        seen_ids = set()
        
        for term in search_terms:
            # 搜索GEO数据集
            search_url = f"{self.API_BASE_URL}/esearch.fcgi"
            params = {
                'db': 'gds',
                'term': term,
                'retmax': 500,
                'retmode': 'json'
            }
            
            response = self.safe_request(search_url, params=params)
            if not response:
                continue
            
            search_data = response.json()
            id_list = safe_get(search_data, 'esearchresult', 'idlist', default=[])
            
            # 过滤已见过的ID
            new_ids = [id for id in id_list if id not in seen_ids]
            seen_ids.update(new_ids)
            
            if not new_ids:
                continue
            
            # 批量获取详情
            for i in range(0, len(new_ids), 50):
                batch_ids = new_ids[i:i+50]
                
                summary_url = f"{self.API_BASE_URL}/esummary.fcgi"
                params = {
                    'db': 'gds',
                    'id': ','.join(batch_ids),
                    'retmode': 'json'
                }
                
                response = self.safe_request(summary_url, params=params)
                if response:
                    summary_data = response.json()
                    result = summary_data.get('result', {})
                    
                    for geo_id in batch_ids:
                        if geo_id in result:
                            datasets.append(result[geo_id])
                
                self.logger.info(f"已获取 {len(datasets)} 个GEO数据集")
                time.sleep(0.5)  # 避免请求过快
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化NCBI GEO数据"""
        accession = raw_data.get('accession', '')
        
        # 从summary中提取信息
        summary = raw_data.get('summary', '')
        
        # 尝试从summary中提取细胞数
        n_cells = None
        cell_patterns = [
            r'(\d+[\d,]*)\s*cells',
            r'(\d+[\d,]*)\s*single[\s-]cells',
        ]
        for pattern in cell_patterns:
            match = re.search(pattern, summary, re.IGNORECASE)
            if match:
                n_cells_str = match.group(1).replace(',', '')
                try:
                    n_cells = int(n_cells_str)
                    break
                except ValueError:
                    pass
        
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': raw_data.get('uid'),
            'accession': accession,
            'title': raw_data.get('title'),
            'n_cells': n_cells,
            'organism': raw_data.get('taxon'),
            'technology': raw_data.get('gdstype'),
            'publication_date': raw_data.get('pdat'),
            'url': f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
            'raw_metadata': raw_data
        }
        
        return normalized


class BroadSCPCollector(BaseCollector):
    """Broad Single Cell Portal 收集器"""
    
    DATABASE_NAME = "broad_scp"
    API_BASE_URL = "https://singlecell.broadinstitute.org/single_cell/api/v1"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取Broad SCP数据集"""
        datasets = []
        
        # 获取所有公开研究
        studies_url = f"{self.API_BASE_URL}/studies"
        params = {
            'page': 1,
            'per_page': 100
        }
        
        while True:
            response = self.safe_request(studies_url, params=params)
            
            if not response:
                break
            
            data = response.json()
            studies = data.get('studies', [])
            
            if not studies:
                break
            
            datasets.extend(studies)
            self.logger.info(f"已获取 {len(datasets)} 个研究")
            
            # 检查是否还有更多页
            if len(studies) < params['per_page']:
                break
            
            params['page'] += 1
            time.sleep(0.5)
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化Broad SCP数据"""
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': raw_data.get('accession'),
            'accession': raw_data.get('accession'),
            'title': raw_data.get('name'),
            'n_cells': raw_data.get('cell_count'),
            'organism': raw_data.get('species'),
            'publication_date': raw_data.get('created_at'),
            'url': f"https://singlecell.broadinstitute.org/single_cell/study/{raw_data.get('accession')}",
            'raw_metadata': raw_data
        }
        
        return normalized


class SCXACollector(BaseCollector):
    """EMBL-EBI Single Cell Expression Atlas 收集器"""
    
    DATABASE_NAME = "scxa"
    API_BASE_URL = "https://www.ebi.ac.uk/gxa/sc/json/experiments"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取SCXA数据集"""
        datasets = []
        
        response = self.safe_request(self.API_BASE_URL)
        
        if response:
            data = response.json()
            experiments = data.get('experiments', [])
            datasets.extend(experiments)
            self.logger.info(f"获取到 {len(datasets)} 个实验")
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化SCXA数据"""
        accession = raw_data.get('experimentAccession', '')
        
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': accession,
            'accession': accession,
            'title': raw_data.get('experimentDescription'),
            'organism': raw_data.get('species'),
            'technology': raw_data.get('experimentType'),
            'last_update': raw_data.get('lastUpdate'),
            'url': f"https://www.ebi.ac.uk/gxa/sc/experiments/{accession}",
            'raw_metadata': raw_data
        }
        
        return normalized


class CNGBdbCollector(BaseCollector):
    """CNGBdb (China National GeneBank) 收集器"""
    
    DATABASE_NAME = "cngbdb"
    API_BASE_URL = "https://ftp.cngb.org/pub/SciRAID"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取CNGBdb数据集"""
        # CNGBdb 需要特殊的API访问权限或网页爬取
        # 这里提供基础框架，实际实现需要根据API文档调整
        self.logger.warning("CNGBdb收集器需要进一步实现具体API调用")
        return []
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化CNGBdb数据"""
        return {
            'database': self.DATABASE_NAME,
            'raw_metadata': raw_data
        }


class HeidelbergCollector(BaseCollector):
    """Heidelberg Reuse-Trends 收集器"""
    
    DATABASE_NAME = "heidelberg"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取Heidelberg数据集"""
        # Heidelberg是一个分析工具而非数据库，这里提供占位符
        self.logger.warning("Heidelberg是分析工具，不是数据存储库")
        return []
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化Heidelberg数据"""
        return {
            'database': self.DATABASE_NAME,
            'raw_metadata': raw_data
        }


class NCBISRACollector(BaseCollector):
    """NCBI SRA 收集器"""
    
    DATABASE_NAME = "ncbi_sra"
    API_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取NCBI SRA数据集"""
        datasets = []
        
        # 搜索单细胞RNA-seq数据
        search_url = f"{self.API_BASE_URL}/esearch.fcgi"
        params = {
            'db': 'sra',
            'term': '("single cell"[Title] OR "scRNA-seq"[Title]) AND "rna seq"[Strategy]',
            'retmax': 500,
            'retmode': 'json',
            'usehistory': 'y'
        }
        
        response = self.safe_request(search_url, params=params)
        if not response:
            return datasets
        
        search_data = response.json()
        id_list = safe_get(search_data, 'esearchresult', 'idlist', default=[])
        
        # 批量获取详情
        for i in range(0, len(id_list), 50):
            batch_ids = id_list[i:i+50]
            
            summary_url = f"{self.API_BASE_URL}/esummary.fcgi"
            params = {
                'db': 'sra',
                'id': ','.join(batch_ids),
                'retmode': 'json'
            }
            
            response = self.safe_request(summary_url, params=params)
            if response:
                summary_data = response.json()
                result = summary_data.get('result', {})
                
                for sra_id in batch_ids:
                    if sra_id in result:
                        datasets.append(result[sra_id])
            
            self.logger.info(f"已获取 {len(datasets)} 个SRA数据集")
            time.sleep(0.5)
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化NCBI SRA数据"""
        # 从expxml中提取信息
        expxml = raw_data.get('expxml', '')
        
        # 解析Run信息
        runs = raw_data.get('runs', '').split(',') if raw_data.get('runs') else []
        accession = runs[0] if runs else raw_data.get('uid')
        
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': raw_data.get('uid'),
            'accession': accession,
            'title': raw_data.get('title'),
            'organism': raw_data.get('organism'),
            'technology': raw_data.get('platform'),
            'publication_date': raw_data.get('createdate'),
            'last_update': raw_data.get('updatedate'),
            'url': f"https://www.ncbi.nlm.nih.gov/sra/{accession}",
            'raw_metadata': raw_data
        }
        
        return normalized


class FigshareCollector(BaseCollector):
    """Figshare 收集器"""
    
    DATABASE_NAME = "figshare"
    API_BASE_URL = "https://api.figshare.com/v2"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取Figshare数据集"""
        datasets = []
        
        # 搜索单细胞相关数据
        search_url = f"{self.API_BASE_URL}/articles/search"
        
        for page in range(1, 11):  # 限制前10页
            params = {
                'search_for': 'single cell RNA-seq',
                'page': page,
                'page_size': 100
            }
            
            response = self.safe_request(search_url, method='POST', json=params)
            
            if not response:
                break
            
            articles = response.json()
            
            if not articles:
                break
            
            datasets.extend(articles)
            self.logger.info(f"已获取 {len(datasets)} 个Figshare数据集")
            time.sleep(1)
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化Figshare数据"""
        # 计算总大小
        total_size_bytes = sum(file.get('size', 0) for file in raw_data.get('files', []))
        
        # 提取作者
        authors = [author.get('full_name', '') for author in raw_data.get('authors', [])]
        
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': raw_data.get('id'),
            'title': raw_data.get('title'),
            'data_size_gb': round(total_size_bytes / (1024**3), 2) if total_size_bytes > 0 else None,
            'doi': raw_data.get('doi'),
            'publication_date': raw_data.get('published_date'),
            'authors': ', '.join(authors) if authors else None,
            'url': raw_data.get('url'),
            'download_url': raw_data.get('url_public_html'),
            'raw_metadata': raw_data
        }
        
        return normalized


class DryadCollector(BaseCollector):
    """Dryad 收集器"""
    
    DATABASE_NAME = "dryad"
    API_BASE_URL = "https://datadryad.org/api/v2"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取Dryad数据集"""
        datasets = []
        
        # 搜索单细胞数据
        search_url = f"{self.API_BASE_URL}/search"
        
        for page in range(1, 11):
            params = {
                'q': 'single cell RNA-seq',
                'page': page,
                'per_page': 100
            }
            
            response = self.safe_request(search_url, params=params)
            
            if not response:
                break
            
            data = response.json()
            items = data.get('_embedded', {}).get('stash:datasets', [])
            
            if not items:
                break
            
            datasets.extend(items)
            self.logger.info(f"已获取 {len(datasets)} 个Dryad数据集")
            time.sleep(1)
        
        return datasets
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化Dryad数据"""
        # 提取作者
        authors = [author.get('author', '') for author in raw_data.get('authors', [])]
        
        # 计算总大小
        total_size_bytes = raw_data.get('storageSize', 0)
        
        normalized = {
            'database': self.DATABASE_NAME,
            'dataset_id': raw_data.get('id'),
            'title': raw_data.get('title'),
            'data_size_gb': round(total_size_bytes / (1024**3), 2) if total_size_bytes > 0 else None,
            'doi': safe_get(raw_data, 'identifier'),
            'publication_date': raw_data.get('publicationDate'),
            'last_update': raw_data.get('lastModificationDate'),
            'authors': ', '.join(authors) if authors else None,
            'url': safe_get(raw_data, '_links', 'stash:dataset', 'href'),
            'raw_metadata': raw_data
        }
        
        return normalized


class HTANCollector(BaseCollector):
    """HTAN (Human Tumor Atlas Network) 收集器"""
    
    DATABASE_NAME = "htan"
    API_BASE_URL = "https://data.humantumoratlas.org/api"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取HTAN数据集"""
        # HTAN的API可能需要特殊认证，这里提供基础框架
        self.logger.warning("HTAN收集器需要进一步实现具体API调用")
        return []
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化HTAN数据"""
        return {
            'database': self.DATABASE_NAME,
            'raw_metadata': raw_data
        }


class PsychANCollector(BaseCollector):
    """PsychAN (Psychencode Atlas) 收集器"""
    
    DATABASE_NAME = "psychan"
    
    def fetch_datasets(self) -> List[Dict]:
        """获取PsychAN数据集"""
        # PsychAN可能需要特殊访问，这里提供基础框架
        self.logger.warning("PsychAN收集器需要进一步实现具体API调用")
        return []
    
    def normalize_dataset(self, raw_data: Dict) -> Dict:
        """规范化PsychAN数据"""
        return {
            'database': self.DATABASE_NAME,
            'raw_metadata': raw_data
        }


# ============================================================================
# 主程序
# ============================================================================

# 数据库收集器映射
DATABASE_COLLECTORS = {
    'broad_scp': BroadSCPCollector,
    'scxa': SCXACollector,
    'cellxgene': CELLxGENECollector,
    'cngbdb': CNGBdbCollector,
    'heidelberg': HeidelbergCollector,
    'ncbi_sra': NCBISRACollector,
    'ncbi_geo': NCBIGEOCollector,
    'hca': HCACollector,
    'figshare': FigshareCollector,
    'dryad': DryadCollector,
    'htan': HTANCollector,
    'psychan': PsychANCollector
}


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='单细胞数据库元数据采集工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 采集所有数据库
  python sc_database_collector.py --all --output ./output
  
  # 采集特定数据库
  python sc_database_collector.py --databases cellxgene hca --output ./output
  
  # 并行采集，使用8个线程
  python sc_database_collector.py --all --parallel --workers 8
  
  # 仅去重和合并已有数据
  python sc_database_collector.py --merge-only --input ./raw_data --output ./output
  
  # 跳过引用量查询以加快速度
  python sc_database_collector.py --all --skip-citation --output ./output
        """
    )
    
    # 数据库选择参数
    db_group = parser.add_mutually_exclusive_group(required=True)
    db_group.add_argument(
        '--all',
        action='store_true',
        help='采集所有支持的数据库'
    )
    db_group.add_argument(
        '--databases',
        nargs='+',
        choices=list(DATABASE_COLLECTORS.keys()),
        help='指定要采集的数据库'
    )
    db_group.add_argument(
        '--merge-only',
        action='store_true',
        help='仅合并和去重已有数据，不进行采集'
    )
    
    # 输出参数
    parser.add_argument(
        '--output',
        type=str,
        default='./output',
        help='输出目录路径 (默认: ./output)'
    )
    parser.add_argument(
        '--input',
        type=str,
        default='./raw_data',
        help='原始数据输入目录 (用于merge-only模式, 默认: ./raw_data)'
    )
    
    # 并行参数
    parser.add_argument(
        '--parallel',
        action='store_true',
        help='启用并行采集'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=4,
        help='并行工作线程数 (默认: 4)'
    )
    
    # 采集参数
    parser.add_argument(
        '--max-records',
        type=int,
        default=None,
        help='每个数据库最大采集记录数 (默认: 无限制)'
    )
    parser.add_argument(
        '--timeout',
        type=int,
        default=300,
        help='单个请求超时时间(秒) (默认: 300)'
    )
    parser.add_argument(
        '--retry',
        type=int,
        default=3,
        help='失败重试次数 (默认: 3)'
    )
    
    # 去重参数
    parser.add_argument(
        '--dedup-keys',
        nargs='+',
        default=['doi', 'n_cells'],
        help='去重使用的字段组合 (默认: doi n_cells)'
    )
    
    # 其他参数
    parser.add_argument(
        '--skip-citation',
        action='store_true',
        help='跳过引用量查询(加快采集速度)'
    )
    parser.add_argument(
        '--cache',
        action='store_true',
        help='启用本地缓存(暂未实现)'
    )
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='日志级别 (默认: INFO)'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        default=None,
        help='日志文件路径 (默认: 仅输出到控制台)'
    )
    
    return parser.parse_args()


def collect_database(collector_class, output_dir, config):
    """
    采集单个数据库
    
    Args:
        collector_class: 收集器类
        output_dir: 输出目录
        config: 配置字典
    
    Returns:
        dict: 采集结果统计
    """
    logger = logging.getLogger(__name__)
    db_name = collector_class.DATABASE_NAME
    
    try:
        logger.info(f"开始采集数据库: {db_name}")
        
        # 初始化收集器
        collector = collector_class(
            max_records=config.get('max_records'),
            timeout=config.get('timeout'),
            retry=config.get('retry'),
            skip_citation=config.get('skip_citation'),
            use_cache=config.get('cache')
        )
        
        # 采集数据
        datasets = collector.collect()
        
        # 保存原始数据
        raw_file = output_dir / 'raw' / f'{db_name}.json'
        raw_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(raw_file, 'w', encoding='utf-8') as f:
            json.dump(datasets, f, ensure_ascii=False, indent=2)
        
        logger.info(f"数据库 {db_name} 采集完成: {len(datasets)} 条记录")
        
        return {
            'database': db_name,
            'status': 'success',
            'records': len(datasets),
            'file': str(raw_file)
        }
        
    except Exception as e:
        logger.error(f"数据库 {db_name} 采集失败: {str(e)}", exc_info=True)
        return {
            'database': db_name,
            'status': 'failed',
            'records': 0,
            'error': str(e)
        }


def main():
    """主函数"""
    args = parse_arguments()
    
    # 设置日志
    setup_logging(args.log_level, args.log_file)
    logger = logging.getLogger(__name__)
    
    logger.info("=" * 80)
    logger.info("单细胞数据库元数据采集系统 v1.0.0")
    logger.info("=" * 80)
    
    # 验证输出目录
    output_dir = Path(args.output)
    validate_output_dir(output_dir)
    
    # 如果是仅合并模式
    if args.merge_only:
        logger.info("仅合并模式: 合并和去重已有数据")
        input_dir = Path(args.input) / 'raw'
        
        if not input_dir.exists():
            logger.error(f"输入目录不存在: {input_dir}")
            sys.exit(1)
        
        # 加载所有JSON文件
        all_datasets = []
        for json_file in sorted(input_dir.glob('*.json')):
            try:
                with open(json_file, 'r', encoding='utf-8') as f:
                    datasets = json.load(f)
                    all_datasets.extend(datasets)
                logger.info(f"✓ 加载文件: {json_file.name}, {len(datasets)} 条记录")
            except Exception as e:
                logger.error(f"✗ 加载文件 {json_file} 失败: {e}")
        
        logger.info(f"\n合并前总记录数: {len(all_datasets)}")
        
        # 去重
        deduplicated = deduplicate_datasets(all_datasets, args.dedup_keys)
        logger.info(f"去重后记录数: {len(deduplicated)}")
        
        # 导出CSV
        csv_file = output_dir / f'sc_metadata_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'
        export_to_csv(deduplicated, csv_file)
        
        logger.info(f"\n统计表格已保存: {csv_file}")
        logger.info("=" * 80)
        logger.info("合并完成!")
        logger.info("=" * 80)
        return
    
    # 确定要采集的数据库
    if args.all:
        databases_to_collect = list(DATABASE_COLLECTORS.keys())
    else:
        databases_to_collect = args.databases
    
    logger.info(f"\n将采集以下数据库 ({len(databases_to_collect)}个):")
    for i, db in enumerate(databases_to_collect, 1):
        logger.info(f"  {i}. {db}")
    
    # 配置参数
    config = {
        'max_records': args.max_records,
        'timeout': args.timeout,
        'retry': args.retry,
        'skip_citation': args.skip_citation,
        'cache': args.cache
    }
    
    logger.info("\n" + "=" * 80)
    logger.info("开始采集数据...")
    logger.info("=" * 80)
    
    # 采集数据
    results = []
    start_time = time.time()
    
    if args.parallel:
        # 并行采集
        logger.info(f"\n使用并行模式采集, 工作线程数: {args.workers}\n")
        
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    collect_database,
                    DATABASE_COLLECTORS[db_name],
                    output_dir,
                    config
                ): db_name
                for db_name in databases_to_collect
            }
            
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                
                # 实时显示进度
                completed = len(results)
                total = len(databases_to_collect)
                logger.info(f"进度: {completed}/{total} 完成")
    else:
        # 串行采集
        logger.info("\n使用串行模式采集\n")
        for i, db_name in enumerate(databases_to_collect, 1):
            logger.info(f"\n[{i}/{len(databases_to_collect)}] 正在采集 {db_name}...")
            result = collect_database(
                DATABASE_COLLECTORS[db_name],
                output_dir,
                config
            )
            results.append(result)
    
    total_time = time.time() - start_time
    
    # 统计采集结果
    logger.info("\n" + "=" * 80)
    logger.info("采集结果统计")
    logger.info("=" * 80)
    
    total_records = 0
    success_count = 0
    failed_databases = []
    
    for result in sorted(results, key=lambda x: x['database']):
        if result['status'] == 'success':
            success_count += 1
            total_records += result['records']
            logger.info(f"✓ {result['database']:20s}: {result['records']:6d} 条记录")
        else:
            failed_databases.append(result['database'])
            logger.error(f"✗ {result['database']:20s}: {result.get('error', 'Unknown error')}")
    
    logger.info("-" * 80)
    logger.info(f"成功: {success_count}/{len(results)} 个数据库")
    logger.info(f"总记录数: {total_records:,}")
    logger.info(f"总耗时: {total_time:.2f} 秒")
    
    if failed_databases:
        logger.warning(f"\n失败的数据库: {', '.join(failed_databases)}")
    
    # 合并和去重
    logger.info("\n" + "=" * 80)
    logger.info("开始合并和去重数据...")
    logger.info("=" * 80)
    
    all_datasets = []
    raw_dir = output_dir / 'raw'
    
    for json_file in sorted(raw_dir.glob('*.json')):
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                datasets = json.load(f)
                all_datasets.extend(datasets)
        except Exception as e:
            logger.error(f"加载文件 {json_file} 失败: {e}")
    
    logger.info(f"合并前总记录数: {len(all_datasets):,}")
    
    # 去重
    deduplicated = deduplicate_datasets(all_datasets, args.dedup_keys)
    dedup_rate = (1 - len(deduplicated)/len(all_datasets))*100 if all_datasets else 0
    logger.info(f"去重后记录数: {len(deduplicated):,}")
    logger.info(f"去重率: {dedup_rate:.2f}%")
    
    # 导出CSV
    csv_file = output_dir / f'sc_metadata_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'
    export_to_csv(deduplicated, csv_file)
    
    # 保存采集报告
    report = {
        'timestamp': datetime.now().isoformat(),
        'databases_collected': len(databases_to_collect),
        'successful': success_count,
        'failed': len(failed_databases),
        'failed_databases': failed_databases,
        'total_records_raw': len(all_datasets),
        'total_records_dedup': len(deduplicated),
        'dedup_rate': f"{dedup_rate:.2f}%",
        'total_time_seconds': round(total_time, 2),
        'config': config,
        'results': results
    }
    
    report_file = output_dir / f'report_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json'
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    
    # 最终总结
    logger.info("\n" + "=" * 80)
    logger.info("采集任务完成!")
    logger.info("=" * 80)
    logger.info(f"统计表格: {csv_file}")
    logger.info(f"采集报告: {report_file}")
    logger.info(f"原始数据: {raw_dir}")
    logger.info("=" * 80)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n用户中断操作")
        sys.exit(1)
    except Exception as e:
        logging.error(f"程序异常退出: {e}", exc_info=True)
        sys.exit(1)