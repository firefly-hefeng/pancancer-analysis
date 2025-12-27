# improved_collector_fixed.py
"""
修复版人类单细胞RNA-seq Metadata收集器
针对性解决每个数据库的具体问题
"""

import os
import json
import requests
import pandas as pd
from pathlib import Path
from datetime import datetime
import time
from typing import Dict, List, Optional
import logging
import urllib3
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

# 禁用SSL警告
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('metadata_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def create_session_with_retries(retries=3, backoff_factor=0.5, 
                                status_forcelist=(500, 502, 504), 
                                verify_ssl=True):
    """创建带重试机制的session"""
    session = requests.Session()
    
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    session.verify = verify_ssl
    
    # 设置合理的超时
    session.timeout = 60
    
    return session


class MetadataCollector:
    """基础收集器类"""
    
    def __init__(self, output_dir: str = "metadata_output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
    def save_raw_data(self, data: any, filename: str, database_name: str):
        """保存原始数据"""
        db_dir = self.output_dir / database_name / "raw"
        db_dir.mkdir(parents=True, exist_ok=True)
        
        filepath = db_dir / filename
        
        if isinstance(data, (dict, list)):
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        elif isinstance(data, pd.DataFrame):
            data.to_csv(filepath.with_suffix('.csv'), index=False)
        else:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(str(data))
                
        logger.info(f"Saved raw data to {filepath}")
        return filepath


class FixedCELLxGENECollector(MetadataCollector):
    """修复版CELLxGENE收集器 - 解决数据类型不一致问题"""
    
    def __init__(self, output_dir: str = "metadata_output"):
        super().__init__(output_dir)
        self.database_name = "cellxgene"
        self.api_base = "https://api.cellxgene.cziscience.com"
        self.census_api = "https://api.cellxgene.cziscience.com/curation/v1"
        self.session = create_session_with_retries()
        
    def _safe_join_list(self, data: any) -> str:
        """安全地连接列表数据，处理混合类型"""
        if not data:
            return ''
        
        if isinstance(data, str):
            return data
        
        if isinstance(data, list):
            result = []
            for item in data:
                if isinstance(item, dict):
                    # 如果是字典，尝试提取label或ontology_term_id
                    result.append(item.get('label', item.get('ontology_term_id', str(item))))
                elif isinstance(item, str):
                    result.append(item)
                else:
                    result.append(str(item))
            return ', '.join(result)
        
        return str(data)
    
    def fetch_collections(self) -> List[Dict]:
        """获取所有collections"""
        url = f"{self.census_api}/collections"
        logger.info(f"Fetching CELLxGENE collections from {url}")
        
        try:
            response = self.session.get(url, timeout=60)
            response.raise_for_status()
            collections = response.json()
            
            self.save_raw_data(collections, "collections_raw.json", self.database_name)
            logger.info(f"Retrieved {len(collections)} collections")
            return collections
            
        except Exception as e:
            logger.error(f"Error fetching CELLxGENE collections: {e}")
            return []
    
    def fetch_datasets(self) -> List[Dict]:
        """获取所有datasets"""
        url = f"{self.census_api}/datasets"
        logger.info(f"Fetching CELLxGENE datasets from {url}")
        
        try:
            response = self.session.get(url, timeout=120)
            response.raise_for_status()
            datasets = response.json()
            
            self.save_raw_data(datasets, "datasets_raw.json", self.database_name)
            logger.info(f"Retrieved {len(datasets)} datasets")
            return datasets
            
        except requests.exceptions.Timeout:
            logger.warning("Initial request timed out, trying from collections...")
            return self._fetch_datasets_from_collections()
            
        except Exception as e:
            logger.error(f"Error fetching CELLxGENE datasets: {e}")
            return []
    
    def _fetch_datasets_from_collections(self) -> List[Dict]:
        """从collections中提取datasets信息"""
        collections = self.fetch_collections()
        all_datasets = []
        
        for collection in collections:
            datasets = collection.get('datasets', [])
            for dataset in datasets:
                dataset['collection_info'] = {
                    'collection_id': collection.get('collection_id'),
                    'doi': collection.get('doi'),
                    'contact_name': collection.get('contact_name'),
                    'contact_email': collection.get('contact_email')
                }
                all_datasets.append(dataset)
        
        logger.info(f"Extracted {len(all_datasets)} datasets from collections")
        self.save_raw_data(all_datasets, "datasets_from_collections.json", self.database_name)
        return all_datasets
    
    def parse_to_series_table(self, collections: List[Dict], datasets: List[Dict]) -> pd.DataFrame:
        """解析为series表格 - 修复数据类型问题"""
        series_list = []
        
        # 创建dataset到collection的映射
        dataset_to_collection = {}
        for coll in collections:
            for ds in coll.get('datasets', []):
                ds_id = ds.get('dataset_id') or ds.get('id')
                dataset_to_collection[ds_id] = coll
        
        for dataset in datasets:
            try:
                dataset_id = dataset.get('dataset_id') or dataset.get('id', '')
                
                # 获取collection信息
                collection = dataset_to_collection.get(dataset_id, {})
                if not collection and 'collection_info' in dataset:
                    collection = dataset['collection_info']
                
                # 安全提取各种字段
                disease_names = self._safe_join_list(dataset.get('disease', []))
                tissue = self._safe_join_list(dataset.get('tissue', []))
                tissue_ontology = self._safe_join_list(dataset.get('tissue_ontology_term_ids', []))
                assay = self._safe_join_list(dataset.get('assay', []))
                assay_ontology = self._safe_join_list(dataset.get('assay_ontology_term_ids', []))
                ethnicity = self._safe_join_list(dataset.get('self_reported_ethnicity_ontology_term_ids', []))
                sex = self._safe_join_list(dataset.get('sex_ontology_term_ids', []))
                
                series_data = {
                    'id': dataset_id,
                    'title': dataset.get('title', ''),
                    'disease_general': self._extract_disease_general(disease_names),
                    'disease': disease_names,
                    'pubmed': collection.get('doi', '') if isinstance(collection, dict) else '',
                    'source_database': 'CELLxGENE',
                    'access_link': f"https://cellxgene.cziscience.com/collections/{collection.get('collection_id', '') if isinstance(collection, dict) else ''}",
                    'open_status': 'Open',
                    'ethnicity': ethnicity,
                    'sex': sex,
                    'tissue': tissue or tissue_ontology,
                    'sequencing_platform': assay or assay_ontology,
                    'experiment_design': assay or 'scRNA-seq',
                    'sample_type': self._determine_sample_type(disease_names),
                    'summary': dataset.get('dataset_deployments', [{}])[0].get('description', '') if dataset.get('dataset_deployments') else '',
                    'citation_count': collection.get('citation_count', 0) if isinstance(collection, dict) else 0,
                    'publication_date': collection.get('published_at', '') if isinstance(collection, dict) else '',
                    'submission_date': dataset.get('created_at', ''),
                    'last_update_date': dataset.get('updated_at', ''),
                    'contact_name': collection.get('contact_name', '') if isinstance(collection, dict) else '',
                    'contact_email': collection.get('contact_email', '') if isinstance(collection, dict) else '',
                    'contact_institute': '',
                    'data_tier': 'processed_matrix',
                    'supplementary_information': json.dumps({
                        'cell_count': dataset.get('cell_count', 0),
                        'organism': self._safe_join_list(dataset.get('organism_ontology_term_ids', [])),
                        'development_stage': self._safe_join_list(dataset.get('development_stage_ontology_term_ids', []))
                    })
                }
                
                series_list.append(series_data)
                
            except Exception as e:
                logger.warning(f"Error parsing dataset {dataset.get('dataset_id', dataset.get('id'))}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
        # 保存处理后的数据
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "series_table.csv", index=False)
        
        logger.info(f"Created series table with {len(df)} entries")
        return df
    
    def parse_to_sample_table(self, datasets: List[Dict]) -> pd.DataFrame:
        """解析为sample表格"""
        sample_list = []
        
        for dataset in datasets:
            try:
                dataset_id = dataset.get('dataset_id') or dataset.get('id', '')
                
                disease_names = self._safe_join_list(dataset.get('disease', []))
                tissue = self._safe_join_list(dataset.get('tissue', []))
                tissue_ontology = self._safe_join_list(dataset.get('tissue_ontology_term_ids', []))
                assay = self._safe_join_list(dataset.get('assay', []))
                assay_ontology = self._safe_join_list(dataset.get('assay_ontology_term_ids', []))
                ethnicity = self._safe_join_list(dataset.get('self_reported_ethnicity_ontology_term_ids', []))
                sex = self._safe_join_list(dataset.get('sex_ontology_term_ids', []))
                
                sample_data = {
                    'id': dataset_id,
                    'sample_id': dataset_id,
                    'title': dataset.get('title', ''),
                    'disease_general': self._extract_disease_general(disease_names),
                    'disease': disease_names,
                    'ethnicity': ethnicity,
                    'sex': sex,
                    'tissue_location': tissue or tissue_ontology,
                    'sequencing_platform': assay or assay_ontology,
                    'experiment_design': assay or 'scRNA-seq',
                    'sample_type': self._determine_sample_type(disease_names)
                }
                
                sample_list.append(sample_data)
                
            except Exception as e:
                logger.warning(f"Error parsing sample from dataset {dataset.get('dataset_id', dataset.get('id'))}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "sample_table.csv", index=False)
        
        logger.info(f"Created sample table with {len(df)} entries")
        return df
    
    def _extract_disease_general(self, disease_names: str) -> str:
        """提取疾病的一般分类"""
        if not disease_names:
            return 'Normal'
        
        disease_str = disease_names.lower()
        
        if 'normal' in disease_str or 'healthy' in disease_str:
            return 'Normal'
        elif 'cancer' in disease_str or 'carcinoma' in disease_str or 'tumor' in disease_str or 'malignant' in disease_str:
            return 'Cancer'
        elif 'covid' in disease_str or 'sars' in disease_str or 'infection' in disease_str:
            return 'Infectious Disease'
        elif 'diabetes' in disease_str:
            return 'Metabolic Disease'
        elif 'alzheimer' in disease_str or 'parkinson' in disease_str or 'neurodegenerat' in disease_str:
            return 'Neurodegenerative Disease'
        elif 'autoimmune' in disease_str or 'lupus' in disease_str or 'arthritis' in disease_str:
            return 'Autoimmune Disease'
        else:
            return 'Other Disease'
    
    def _determine_sample_type(self, disease_names: str) -> str:
        """判断样本类型"""
        if not disease_names:
            return 'healthy'
        
        disease_str = disease_names.lower()
        if 'normal' in disease_str or 'healthy' in disease_str:
            return 'healthy'
        else:
            return 'diseased'
    
    def collect(self):
        """执行完整的收集流程"""
        logger.info("Starting fixed CELLxGENE collection...")
        
        collections = self.fetch_collections()
        datasets = self.fetch_datasets()
        
        if not datasets and collections:
            logger.info("Extracting datasets from collections...")
            datasets = []
            for coll in collections:
                datasets.extend(coll.get('datasets', []))
        
        if collections or datasets:
            series_df = self.parse_to_series_table(collections, datasets)
            sample_df = self.parse_to_sample_table(datasets)
            
            return {
                'series': series_df,
                'samples': sample_df
            }
        else:
            logger.error("Failed to collect CELLxGENE data")
            return None


class FixedBroadSingleCellCollector(MetadataCollector):
    """修复版Broad SCP收集器 - 使用替代方案"""
    
    def __init__(self, output_dir: str = "metadata_output"):
        super().__init__(output_dir)
        self.database_name = "broad_single_cell"
        # 尝试使用公开数据页面
        self.web_base = "https://singlecell.broadinstitute.org"
        self.api_base = "https://singlecell.broadinstitute.org/single_cell/api/v1"
        self.session = create_session_with_retries(verify_ssl=False, retries=5)
        
    def fetch_studies_alternative(self) -> List[Dict]:
        """使用替代方法获取研究列表"""
        logger.info("Trying alternative methods for Broad SCP...")
        
        # 方法1: 尝试使用不同的API endpoint
        endpoints = [
            f"{self.api_base}/studies",
            f"{self.api_base}/site/studies",
            f"{self.web_base}/single_cell/data/public"
        ]
        
        for endpoint in endpoints:
            try:
                logger.info(f"Trying endpoint: {endpoint}")
                response = self.session.get(endpoint, timeout=30)
                
                if response.status_code == 200:
                    data = response.json() if 'json' in response.headers.get('content-type', '') else None
                    if data:
                        logger.info(f"Success with endpoint: {endpoint}")
                        return data if isinstance(data, list) else data.get('studies', [])
            except Exception as e:
                logger.debug(f"Failed endpoint {endpoint}: {e}")
                continue
        
        # 方法2: 爬取公开页面（如果API都失败）
        logger.warning("All API endpoints failed, Broad SCP may be temporarily unavailable")
        return []
    
    def fetch_studies(self, page: int = 1, per_page: int = 100) -> Dict:
        """获取研究列表"""
        url = f"{self.api_base}/search/studies"
        params = {
            'page': page,
            'per_page': per_page,
            'scpbr': ''
        }
        
        logger.info(f"Fetching Broad SCP studies (page {page})...")
        
        try:
            response = self.session.get(url, params=params, timeout=60)
            response.raise_for_status()
            data = response.json()
            return data
            
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 500:
                logger.warning("Server error, trying alternative method...")
                return {'studies': self.fetch_studies_alternative()}
            raise
            
        except Exception as e:
            logger.error(f"Error fetching Broad SCP studies: {e}")
            return {}
    
    def fetch_all_studies(self) -> List[Dict]:
        """获取所有研究"""
        all_studies = []
        page = 1
        max_pages = 100
        
        while page <= max_pages:
            try:
                data = self.fetch_studies(page=page)
                studies = data.get('studies', [])
                
                if not studies:
                    break
                
                # 过滤人类数据
                human_studies = [
                    s for s in studies 
                    if any('human' in sp.lower() or 'homo sapiens' in sp.lower() 
                          for sp in s.get('species', []))
                ]
                
                all_studies.extend(human_studies)
                logger.info(f"  Page {page}: found {len(human_studies)} human studies")
                
                if len(studies) < per_page:
                    break
                
                page += 1
                time.sleep(2)
                
            except Exception as e:
                logger.error(f"Error on page {page}: {e}")
                break
        
        self.save_raw_data(all_studies, "studies_raw.json", self.database_name)
        logger.info(f"Retrieved {len(all_studies)} human studies from Broad SCP")
        return all_studies
    
    def parse_to_series_table(self, studies: List[Dict]) -> pd.DataFrame:
        """解析为series表格"""
        series_list = []
        
        for study in studies:
            try:
                study_accession = study.get('accession', '')
                disease = study.get('disease', [])
                organ_region = study.get('organ_region', [])
                library_prep = study.get('library_preparation_protocol', [])
                
                series_data = {
                    'id': study_accession,
                    'title': study.get('name', ''),
                    'disease_general': self._extract_disease_general(disease),
                    'disease': ', '.join(disease) if disease else '',
                    'pubmed': study.get('doi', ''),
                    'source_database': 'Broad Single Cell Portal',
                    'access_link': f"https://singlecell.broadinstitute.org/single_cell/study/{study_accession}",
                    'open_status': 'Open' if study.get('public', True) else 'Restricted',
                    'ethnicity': '',
                    'sex': '',
                    'tissue': ', '.join(organ_region) if organ_region else '',
                    'sequencing_platform': ', '.join(library_prep) if library_prep else '',
                    'experiment_design': 'scRNA-seq',
                    'sample_type': self._determine_sample_type(disease),
                    'summary': study.get('description', ''),
                    'citation_count': 0,
                    'publication_date': study.get('publication_date', ''),
                    'submission_date': study.get('created_at', ''),
                    'last_update_date': study.get('updated_at', ''),
                    'contact_name': study.get('contact_name', ''),
                    'contact_email': study.get('contact_email', ''),
                    'contact_institute': study.get('branding_group_name', ''),
                    'data_tier': 'processed_matrix',
                    'supplementary_information': json.dumps({
                        'cell_count': study.get('cell_count', 0),
                        'species': study.get('species', []),
                        'data_types': study.get('data_types', [])
                    })
                }
                
                series_list.append(series_data)
                
            except Exception as e:
                logger.warning(f"Error parsing study {study.get('accession')}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "series_table.csv", index=False)
        
        logger.info(f"Created series table with {len(df)} entries")
        return df
    
    def parse_to_sample_table(self, studies: List[Dict]) -> pd.DataFrame:
        """解析为sample表格"""
        sample_list = []
        
        for study in studies:
            try:
                study_accession = study.get('accession', '')
                disease = study.get('disease', [])
                organ_region = study.get('organ_region', [])
                library_prep = study.get('library_preparation_protocol', [])
                
                sample_data = {
                    'id': study_accession,
                    'sample_id': study_accession,
                    'title': study.get('name', ''),
                    'disease_general': self._extract_disease_general(disease),
                    'disease': ', '.join(disease) if disease else '',
                    'ethnicity': '',
                    'sex': '',
                    'tissue_location': ', '.join(organ_region) if organ_region else '',
                    'sequencing_platform': ', '.join(library_prep) if library_prep else '',
                    'experiment_design': 'scRNA-seq',
                    'sample_type': self._determine_sample_type(disease)
                }
                
                sample_list.append(sample_data)
                
            except Exception as e:
                logger.warning(f"Error parsing sample from study {study.get('accession')}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "sample_table.csv", index=False)
        
        logger.info(f"Created sample table with {len(df)} entries")
        return df
    
    def _extract_disease_general(self, disease_names: List[str]) -> str:
        """提取疾病的一般分类"""
        if not disease_names:
            return 'Normal'
        
        disease_str = ' '.join(disease_names).lower()
        
        if 'normal' in disease_str or 'healthy' in disease_str:
            return 'Normal'
        elif 'cancer' in disease_str or 'carcinoma' in disease_str or 'tumor' in disease_str:
            return 'Cancer'
        elif 'covid' in disease_str or 'sars' in disease_str:
            return 'Infectious Disease'
        else:
            return 'Other Disease'
    
    def _determine_sample_type(self, disease_names: List[str]) -> str:
        """判断样本类型"""
        if not disease_names:
            return 'healthy'
        
        disease_str = ' '.join(disease_names).lower()
        if 'normal' in disease_str or 'healthy' in disease_str:
            return 'healthy'
        else:
            return 'diseased'
    
    def collect(self):
        """执行完整的收集流程"""
        logger.info("Starting fixed Broad Single Cell Portal collection...")
        
        studies = self.fetch_all_studies()
        
        if studies:
            series_df = self.parse_to_series_table(studies)
            sample_df = self.parse_to_sample_table(studies)
            
            return {
                'series': series_df,
                'samples': sample_df
            }
        else:
            logger.warning("Failed to collect Broad SCP data - service may be temporarily unavailable")
            return None


class FixedZenodoCollector(MetadataCollector):
    """修复版Zenodo收集器 - 解决连接问题"""
    
    def __init__(self, output_dir: str = "metadata_output", api_token: Optional[str] = None, 
                 use_proxy: bool = False, proxy_url: Optional[str] = None):
        super().__init__(output_dir)
        self.database_name = "zenodo"
        self.api_base = "https://zenodo.org/api/records"
        self.api_token = api_token
        self.use_proxy = use_proxy
        self.proxy_url = proxy_url
        
        # 创建session with proxy support
        self.session = create_session_with_retries()
        
        if use_proxy and proxy_url:
            self.session.proxies = {
                'http': proxy_url,
                'https': proxy_url
            }
            logger.info(f"Using proxy: {proxy_url}")
    
    def search_records(self, query: str, page: int = 1, size: int = 100) -> Dict:
        """搜索Zenodo记录"""
        url = self.api_base
        
        params = {
            'q': query,
            'page': page,
            'size': size,
            'sort': 'mostrecent'
        }
        
        headers = {
            'Accept': 'application/json',
            'User-Agent': 'Mozilla/5.0 (compatible; ScRNAseq-MetadataCollector/1.0)'
        }
        
        if self.api_token:
            headers['Authorization'] = f'Bearer {self.api_token}'
        
        logger.info(f"Searching Zenodo (page {page}, query: {query[:50]}...)...")
        
        try:
            response = self.session.get(url, params=params, headers=headers, timeout=60)
            response.raise_for_status()
            return response.json()
            
        except requests.exceptions.ConnectionError as e:
            logger.error(f"Connection error for Zenodo: {e}")
            logger.info("Zenodo may require VPN or proxy access from your location")
            return {}
            
        except Exception as e:
            logger.error(f"Error searching Zenodo: {e}")
            return {}
    
    def fetch_all_scrna_records(self) -> List[Dict]:
        """获取所有单细胞RNA测序相关记录"""
        all_records = []
        
        queries = [
            'single cell RNA sequencing human',
            'scRNA-seq human',
            '10x genomics chromium human'
        ]
        
        for query in queries:
            page = 1
            max_pages = 10
            
            while page <= max_pages:
                data = self.search_records(query, page=page)
                
                if not data:
                    break
                
                hits = data.get('hits', {}).get('hits', [])
                
                if not hits:
                    break
                
                all_records.extend(hits)
                logger.info(f"  Query '{query}' page {page}: found {len(hits)} records")
                
                total = data.get('hits', {}).get('total', 0)
                if page * 100 >= total:
                    break
                
                page += 1
                time.sleep(2)
        
        unique_records = {rec['id']: rec for rec in all_records}.values()
        unique_records = list(unique_records)
        
        self.save_raw_data(unique_records, "records_raw.json", self.database_name)
        
        logger.info(f"Retrieved {len(unique_records)} unique records from Zenodo")
        return unique_records
    
    def parse_to_series_table(self, records: List[Dict]) -> pd.DataFrame:
        """解析为series表格"""
        series_list = []
        
        for record in records:
            try:
                metadata = record.get('metadata', {})
                
                series_data = {
                    'id': str(record.get('id', '')),
                    'title': metadata.get('title', ''),
                    'disease_general': self._extract_disease_from_metadata(metadata),
                    'disease': '',
                    'pubmed': metadata.get('doi', ''),
                    'source_database': 'Zenodo',
                    'access_link': record.get('links', {}).get('html', ''),
                    'open_status': 'Open' if metadata.get('access_right', '') == 'open' else 'Restricted',
                    'ethnicity': '',
                    'sex': '',
                    'tissue': self._extract_tissue_from_metadata(metadata),
                    'sequencing_platform': '',
                    'experiment_design': metadata.get('resource_type', {}).get('type', ''),
                    'sample_type': '',
                    'summary': metadata.get('description', ''),
                    'citation_count': 0,
                    'publication_date': metadata.get('publication_date', ''),
                    'submission_date': record.get('created', ''),
                    'last_update_date': record.get('updated', ''),
                    'contact_name': ', '.join([c.get('name', '') for c in metadata.get('creators', [])]),
                    'contact_email': '',
                    'contact_institute': ', '.join([c.get('affiliation', '') for c in metadata.get('creators', []) if c.get('affiliation')]),
                    'data_tier': 'raw' if 'raw' in metadata.get('description', '').lower() else 'processed_matrix',
                    'supplementary_information': json.dumps({
                        'keywords': metadata.get('keywords', []),
                        'resource_type': metadata.get('resource_type', {}),
                        'file_count': len(record.get('files', []))
                    })
                }
                
                series_list.append(series_data)
                
            except Exception as e:
                logger.warning(f"Error parsing Zenodo record {record.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "series_table.csv", index=False)
        
        logger.info(f"Created series table with {len(df)} entries")
        return df
    
    def parse_to_sample_table(self, records: List[Dict]) -> pd.DataFrame:
        """解析为sample表格"""
        sample_list = []
        
        for record in records:
            try:
                metadata = record.get('metadata', {})
                record_id = str(record.get('id', ''))
                
                sample_data = {
                    'id': record_id,
                    'sample_id': record_id,
                    'title': metadata.get('title', ''),
                    'disease_general': self._extract_disease_from_metadata(metadata),
                    'disease': '',
                    'ethnicity': '',
                    'sex': '',
                    'tissue_location': self._extract_tissue_from_metadata(metadata),
                    'sequencing_platform': '',
                    'experiment_design': metadata.get('resource_type', {}).get('type', ''),
                    'sample_type': 'unknown'
                }
                
                sample_list.append(sample_data)
                
            except Exception as e:
                logger.warning(f"Error parsing sample from Zenodo record {record.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "sample_table.csv", index=False)
        
        logger.info(f"Created sample table with {len(df)} entries")
        return df
    
    def _extract_disease_from_metadata(self, metadata: Dict) -> str:
        """从metadata中提取疾病信息"""
        text = (metadata.get('title', '') + ' ' + 
                metadata.get('description', '') + ' ' + 
                ' '.join(metadata.get('keywords', []))).lower()
        
        if 'cancer' in text or 'tumor' in text or 'carcinoma' in text:
            return 'Cancer'
        elif 'healthy' in text or 'normal' in text:
            return 'Normal'
        elif 'covid' in text or 'sars' in text:
            return 'Infectious Disease'
        else:
            return 'Other Disease'
    
    def _extract_tissue_from_metadata(self, metadata: Dict) -> str:
        """从metadata中提取组织信息"""
        text = (metadata.get('title', '') + ' ' + 
                metadata.get('description', '')).lower()
        
        tissues = []
        tissue_keywords = ['brain', 'lung', 'liver', 'heart', 'kidney', 'blood', 'pbmc', 
                          'tumor', 'skin', 'intestine', 'colon', 'breast', 'ovary', 'testis']
        
        for tissue in tissue_keywords:
            if tissue in text:
                tissues.append(tissue)
        
        return ', '.join(tissues)
    
    def collect(self):
        """执行完整的收集流程"""
        logger.info("Starting fixed Zenodo collection...")
        
        records = self.fetch_all_scrna_records()
        
        if records:
            series_df = self.parse_to_series_table(records)
            sample_df = self.parse_to_sample_table(records)
            
            return {
                'series': series_df,
                'samples': sample_df
            }
        else:
            logger.warning("No records collected from Zenodo")
            logger.info("Note: Zenodo may require VPN/proxy access or API token")
            return None


class FixedFigshareCollector(MetadataCollector):
    """修复版Figshare收集器 - 解决API认证问题"""
    
    def __init__(self, output_dir: str = "metadata_output", api_token: Optional[str] = None):
        super().__init__(output_dir)
        self.database_name = "figshare"
        self.api_base = "https://api.figshare.com/v2"
        self.api_token = api_token
        self.session = create_session_with_retries()
        
    def search_articles(self, search_term: str, page: int = 1, page_size: int = 100) -> List[Dict]:
        """搜索Figshare文章 - 使用公开endpoint"""
        # 使用POST方法和正确的endpoint
        url = f"{self.api_base}/articles/search"
        
        payload = {
            'search_for': f'{search_term} AND (human OR "homo sapiens")',
            'page': page,
            'page_size': page_size,
            'order': 'published_date',
            'order_direction': 'desc'
        }
        
        headers = {
            'Content-Type': 'application/json'
        }
        
        if self.api_token:
            headers['Authorization'] = f'token {self.api_token}'
        
        logger.info(f"Searching Figshare (page {page}, query: {search_term[:50]}...)...")
        
        try:
            # 使用POST方法
            response = self.session.post(url, json=payload, headers=headers, timeout=60)
            response.raise_for_status()
            return response.json()
            
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 403:
                logger.warning("Figshare API forbidden - may need API token")
                logger.info("You can get a free API token from: https://figshare.com/account/applications")
                return []
            else:
                logger.error(f"Error searching Figshare: {e}")
                return []
                
        except Exception as e:
            logger.error(f"Error searching Figshare: {e}")
            return []
    
    def fetch_all_scrna_articles(self) -> List[Dict]:
        """获取所有单细胞RNA测序相关文章"""
        all_articles = []
        
        queries = [
            'single cell RNA-seq',
            'scRNA-seq',
            'single-cell transcriptomics'
        ]
        
        for query in queries:
            page = 1
            max_pages = 5
            
            while page <= max_pages:
                articles = self.search_articles(query, page=page)
                
                if not articles:
                    break
                
                all_articles.extend(articles)
                logger.info(f"  Query '{query}' page {page}: found {len(articles)} articles")
                
                if len(articles) < 100:
                    break
                
                page += 1
                time.sleep(2)
        
        unique_articles = {art['id']: art for art in all_articles}.values()
        unique_articles = list(unique_articles)
        
        self.save_raw_data(unique_articles, "articles_raw.json", self.database_name)
        
        logger.info(f"Retrieved {len(unique_articles)} unique articles from Figshare")
        return unique_articles
    
    def parse_to_series_table(self, articles: List[Dict]) -> pd.DataFrame:
        """解析为series表格"""
        series_list = []
        
        for article in articles:
            try:
                series_data = {
                    'id': str(article.get('id', '')),
                    'title': article.get('title', ''),
                    'disease_general': self._extract_disease_from_text(article.get('description', '')),
                    'disease': '',
                    'pubmed': article.get('doi', ''),
                    'source_database': 'Figshare',
                    'access_link': article.get('url', ''),
                    'open_status': 'Open',
                    'ethnicity': '',
                    'sex': '',
                    'tissue': self._extract_tissue_from_text(article.get('description', '')),
                    'sequencing_platform': '',
                    'experiment_design': article.get('defined_type_name', ''),
                    'sample_type': '',
                    'summary': article.get('description', ''),
                    'citation_count': article.get('citation_count', 0),
                    'publication_date': article.get('published_date', ''),
                    'submission_date': article.get('created_date', ''),
                    'last_update_date': article.get('modified_date', ''),
                    'contact_name': ', '.join([a.get('full_name', '') for a in article.get('authors', [])]),
                    'contact_email': '',
                    'contact_institute': '',
                    'data_tier': 'processed_matrix',
                    'supplementary_information': json.dumps({
                        'tags': article.get('tags', []),
                        'categories': article.get('categories', []),
                        'file_count': len(article.get('files', []))
                    })
                }
                
                series_list.append(series_data)
                
            except Exception as e:
                logger.warning(f"Error parsing Figshare article {article.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "series_table.csv", index=False)
        
        logger.info(f"Created series table with {len(df)} entries")
        return df
    
    def parse_to_sample_table(self, articles: List[Dict]) -> pd.DataFrame:
        """解析为sample表格"""
        sample_list = []
        
        for article in articles:
            try:
                article_id = str(article.get('id', ''))
                
                sample_data = {
                    'id': article_id,
                    'sample_id': article_id,
                    'title': article.get('title', ''),
                    'disease_general': self._extract_disease_from_text(article.get('description', '')),
                    'disease': '',
                    'ethnicity': '',
                    'sex': '',
                    'tissue_location': self._extract_tissue_from_text(article.get('description', '')),
                    'sequencing_platform': '',
                    'experiment_design': article.get('defined_type_name', ''),
                    'sample_type': 'unknown'
                }
                
                sample_list.append(sample_data)
                
            except Exception as e:
                logger.warning(f"Error parsing sample from Figshare article {article.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "sample_table.csv", index=False)
        
        logger.info(f"Created sample table with {len(df)} entries")
        return df
    
    def _extract_disease_from_text(self, text: str) -> str:
        """从文本中提取疾病信息"""
        text_lower = text.lower()
        
        if 'cancer' in text_lower or 'tumor' in text_lower or 'carcinoma' in text_lower:
            return 'Cancer'
        elif 'healthy' in text_lower or 'normal' in text_lower:
            return 'Normal'
        elif 'covid' in text_lower or 'sars' in text_lower:
            return 'Infectious Disease'
        else:
            return 'Other Disease'
    
    def _extract_tissue_from_text(self, text: str) -> str:
        """从文本中提取组织信息"""
        text_lower = text.lower()
        
        tissues = []
        tissue_keywords = ['brain', 'lung', 'liver', 'heart', 'kidney', 'blood', 'pbmc', 
                          'tumor', 'skin', 'intestine', 'colon', 'breast', 'ovary', 'testis']
        
        for tissue in tissue_keywords:
            if tissue in text_lower:
                tissues.append(tissue)
        
        return ', '.join(tissues)
    
    def collect(self):
        """执行完整的收集流程"""
        logger.info("Starting fixed Figshare collection...")
        
        articles = self.fetch_all_scrna_articles()
        
        if articles:
            series_df = self.parse_to_series_table(articles)
            sample_df = self.parse_to_sample_table(articles)
            
            return {
                'series': series_df,
                'samples': sample_df
            }
        else:
            logger.warning("No articles collected from Figshare")
            logger.info("Note: Figshare may require an API token for full access")
            logger.info("Get a free token at: https://figshare.com/account/applications")
            return None


# DRYAD收集器保持不变
class DRYADCollector(MetadataCollector):
    """DRYAD 数据收集器"""
    
    def __init__(self, output_dir: str = "metadata_output"):
        super().__init__(output_dir)
        self.database_name = "dryad"
        self.api_base = "https://datadryad.org/api/v2"
        self.session = create_session_with_retries()
        
    def search_datasets(self, query: str, page: int = 1, per_page: int = 100) -> Dict:
        """搜索DRYAD数据集"""
        url = f"{self.api_base}/search"
        
        params = {
            'q': query,
            'page': page,
            'per_page': per_page
        }
        
        logger.info(f"Searching DRYAD (page {page})...")
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
            
        except Exception as e:
            logger.error(f"Error searching DRYAD: {e}")
            return {}
    
    def fetch_all_scrna_datasets(self) -> List[Dict]:
        """获取所有单细胞RNA测序相关数据集"""
        all_datasets = []
        
        queries = [
            'single cell RNA-seq human',
            'scRNA-seq human',
            'single-cell transcriptomics human'
        ]
        
        for query in queries:
            page = 1
            max_pages = 10
            
            while page <= max_pages:
                data = self.search_datasets(query, page=page)
                datasets = data.get('_embedded', {}).get('stash:datasets', [])
                
                if not datasets:
                    break
                
                all_datasets.extend(datasets)
                logger.info(f"  Query '{query}' page {page}: found {len(datasets)} datasets")
                
                if not data.get('_links', {}).get('next'):
                    break
                
                page += 1
                time.sleep(1)
        
        unique_datasets = {ds['id']: ds for ds in all_datasets}.values()
        unique_datasets = list(unique_datasets)
        
        self.save_raw_data(unique_datasets, "datasets_raw.json", self.database_name)
        
        logger.info(f"Retrieved {len(unique_datasets)} unique datasets from DRYAD")
        return unique_datasets
    
    def parse_to_series_table(self, datasets: List[Dict]) -> pd.DataFrame:
        """解析为series表格"""
        series_list = []
        
        for dataset in datasets:
            try:
                series_data = {
                    'id': str(dataset.get('id', '')),
                    'title': dataset.get('title', ''),
                    'disease_general': self._extract_disease_from_text(dataset.get('abstract', '')),
                    'disease': '',
                    'pubmed': dataset.get('identifier', ''),
                    'source_database': 'DRYAD',
                    'access_link': dataset.get('_links', {}).get('stash:dataset', {}).get('href', ''),
                    'open_status': 'Open',
                    'ethnicity': '',
                    'sex': '',
                    'tissue': self._extract_tissue_from_text(dataset.get('abstract', '')),
                    'sequencing_platform': '',
                    'experiment_design': '',
                    'sample_type': '',
                    'summary': dataset.get('abstract', ''),
                    'citation_count': 0,
                    'publication_date': dataset.get('publicationDate', ''),
                    'submission_date': dataset.get('createdAt', ''),
                    'last_update_date': dataset.get('lastModificationDate', ''),
                    'contact_name': ', '.join([a.get('firstName', '') + ' ' + a.get('lastName', '') 
                                              for a in dataset.get('authors', [])]),
                    'contact_email': dataset.get('authors', [{}])[0].get('email', '') if dataset.get('authors') else '',
                    'contact_institute': dataset.get('authors', [{}])[0].get('affiliation', '') if dataset.get('authors') else '',
                    'data_tier': 'raw',
                    'supplementary_information': json.dumps({
                        'keywords': dataset.get('keywords', []),
                        'funding': dataset.get('funders', [])
                    })
                }
                
                series_list.append(series_data)
                
            except Exception as e:
                logger.warning(f"Error parsing DRYAD dataset {dataset.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "series_table.csv", index=False)
        
        logger.info(f"Created series table with {len(df)} entries")
        return df
    
    def parse_to_sample_table(self, datasets: List[Dict]) -> pd.DataFrame:
        """解析为sample表格"""
        sample_list = []
        
        for dataset in datasets:
            try:
                dataset_id = str(dataset.get('id', ''))
                
                sample_data = {
                    'id': dataset_id,
                    'sample_id': dataset_id,
                    'title': dataset.get('title', ''),
                    'disease_general': self._extract_disease_from_text(dataset.get('abstract', '')),
                    'disease': '',
                    'ethnicity': '',
                    'sex': '',
                    'tissue_location': self._extract_tissue_from_text(dataset.get('abstract', '')),
                    'sequencing_platform': '',
                    'experiment_design': '',
                    'sample_type': 'unknown'
                }
                
                sample_list.append(sample_data)
                
            except Exception as e:
                logger.warning(f"Error parsing sample from DRYAD dataset {dataset.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        output_path = self.output_dir / self.database_name / "processed"
        output_path.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path / "sample_table.csv", index=False)
        
        logger.info(f"Created sample table with {len(df)} entries")
        return df
    
    def _extract_disease_from_text(self, text: str) -> str:
        """从文本中提取疾病信息"""
        text_lower = text.lower()
        
        if 'cancer' in text_lower or 'tumor' in text_lower or 'carcinoma' in text_lower:
            return 'Cancer'
        elif 'healthy' in text_lower or 'normal' in text_lower:
            return 'Normal'
        elif 'covid' in text_lower or 'sars' in text_lower:
            return 'Infectious Disease'
        else:
            return 'Other Disease'
    
    def _extract_tissue_from_text(self, text: str) -> str:
        """从文本中提取组织信息"""
        text_lower = text.lower()
        
        tissues = []
        tissue_keywords = ['brain', 'lung', 'liver', 'heart', 'kidney', 'blood', 'pbmc', 
                          'tumor', 'skin', 'intestine', 'colon', 'breast', 'ovary', 'testis']
        
        for tissue in tissue_keywords:
            if tissue in text_lower:
                tissues.append(tissue)
        
        return ', '.join(tissues)
    
    def collect(self):
        """执行完整的收集流程"""
        logger.info("Starting DRYAD collection...")
        
        datasets = self.fetch_all_scrna_datasets()
        
        if datasets:
            series_df = self.parse_to_series_table(datasets)
            sample_df = self.parse_to_sample_table(datasets)
            
            return {
                'series': series_df,
                'samples': sample_df
            }
        else:
            logger.error("Failed to collect DRYAD data")
            return None


class FixedIntegratedCollector:
    """修复版整合收集器"""
    
    def __init__(self, output_dir: str = "metadata_output", 
                 zenodo_use_proxy: bool = False, 
                 zenodo_proxy_url: Optional[str] = None,
                 figshare_api_token: Optional[str] = None,
                 zenodo_api_token: Optional[str] = None):
        self.output_dir = Path(output_dir)
        self.collectors = {
            'cellxgene': FixedCELLxGENECollector(output_dir),
            'broad_single_cell': FixedBroadSingleCellCollector(output_dir),
            'zenodo': FixedZenodoCollector(output_dir, api_token=zenodo_api_token, 
                                          use_proxy=zenodo_use_proxy, 
                                          proxy_url=zenodo_proxy_url),
            'figshare': FixedFigshareCollector(output_dir, api_token=figshare_api_token),
            'dryad': DRYADCollector(output_dir)
        }
        
    def collect_all(self, databases: Optional[List[str]] = None):
        """收集所有数据库的数据"""
        if databases is None:
            databases = list(self.collectors.keys())
        
        results = {}
        
        for db_name in databases:
            if db_name not in self.collectors:
                logger.warning(f"Unknown database: {db_name}")
                continue
            
            logger.info(f"\n{'='*60}\nCollecting from {db_name}\n{'='*60}")
            
            try:
                collector = self.collectors[db_name]
                result = collector.collect()
                results[db_name] = result
                
                if result and isinstance(result, dict):
                    series_count = len(result.get('series', [])) if result.get('series') is not None else 0
                    sample_count = len(result.get('samples', [])) if result.get('samples') is not None else 0
                    logger.info(f"✓ {db_name}: {series_count} series, {sample_count} samples")
                else:
                    logger.warning(f"✗ {db_name}: No data collected")
                
            except Exception as e:
                logger.error(f"Error collecting from {db_name}: {e}")
                import traceback
                logger.error(traceback.format_exc())
                results[db_name] = None
        
        return results
    
    def merge_tables(self, results: Dict):
        """合并所有数据库的表格"""
        all_series = []
        all_samples = []
        
        for db_name, result in results.items():
            if result and isinstance(result, dict):
                if 'series' in result and result['series'] is not None and len(result['series']) > 0:
                    all_series.append(result['series'])
                if 'samples' in result and result['samples'] is not None and len(result['samples']) > 0:
                    all_samples.append(result['samples'])
        
        merged_series = None
        if all_series:
            merged_series = pd.concat(all_series, ignore_index=True)
            merged_series = merged_series.drop_duplicates(subset=['id'], keep='first')
            
            output_path = self.output_dir / "integrated"
            output_path.mkdir(parents=True, exist_ok=True)
            merged_series.to_csv(output_path / "merged_series_table.csv", index=False)
            
            logger.info(f"✓ Merged series table: {len(merged_series)} entries from {len(all_series)} databases")
        else:
            logger.warning("No series data to merge")
        
        merged_samples = None
        if all_samples:
            merged_samples = pd.concat(all_samples, ignore_index=True)
            merged_samples = merged_samples.drop_duplicates(subset=['sample_id'], keep='first')
            
            output_path = self.output_dir / "integrated"
            output_path.mkdir(parents=True, exist_ok=True)
            merged_samples.to_csv(output_path / "merged_sample_table.csv", index=False)
            
            logger.info(f"✓ Merged sample table: {len(merged_samples)} entries from {len(all_samples)} databases")
        else:
            logger.warning("No sample data to merge")
        
        return {
            'series': merged_series,
            'samples': merged_samples
        }
    
    def generate_summary_report(self, results: Dict):
        """生成汇总报告"""
        report = {
            'collection_date': datetime.now().isoformat(),
            'databases': {},
            'summary': {
                'total_series': 0,
                'total_samples': 0,
                'successful_databases': 0,
                'failed_databases': 0
            }
        }
        
        for db_name, result in results.items():
            if result and isinstance(result, dict):
                series_count = len(result['series']) if result.get('series') is not None else 0
                sample_count = len(result['samples']) if result.get('samples') is not None else 0
                
                report['databases'][db_name] = {
                    'series_count': series_count,
                    'sample_count': sample_count,
                    'status': 'success'
                }
                
                report['summary']['total_series'] += series_count
                report['summary']['total_samples'] += sample_count
                report['summary']['successful_databases'] += 1
            else:
                report['databases'][db_name] = {
                    'series_count': 0,
                    'sample_count': 0,
                    'status': 'failed'
                }
                report['summary']['failed_databases'] += 1
        
        output_path = self.output_dir / "integrated"
        output_path.mkdir(parents=True, exist_ok=True)
        
        with open(output_path / "collection_report.json", 'w') as f:
            json.dump(report, f, indent=2)
        
        with open(output_path / "collection_report.txt", 'w') as f:
            f.write("="*60 + "\n")
            f.write("Human scRNA-seq Metadata Collection Report\n")
            f.write("="*60 + "\n\n")
            f.write(f"Collection Date: {report['collection_date']}\n\n")
            
            f.write("Summary:\n")
            f.write("-"*60 + "\n")
            f.write(f"Total Series: {report['summary']['total_series']}\n")
            f.write(f"Total Samples: {report['summary']['total_samples']}\n")
            f.write(f"Successful Databases: {report['summary']['successful_databases']}\n")
            f.write(f"Failed Databases: {report['summary']['failed_databases']}\n\n")
            
            f.write("Database Details:\n")
            f.write("-"*60 + "\n")
            for db_name, db_info in report['databases'].items():
                f.write(f"\n{db_name}:\n")
                f.write(f"  Status: {db_info['status']}\n")
                f.write(f"  Series: {db_info['series_count']}\n")
                f.write(f"  Samples: {db_info['sample_count']}\n")
        
        logger.info(f"✓ Generated summary report")
        return report


def main():
    """主函数"""
    print("="*60)
    print("Fixed Human scRNA-seq Metadata Collector")
    print("="*60)
    print()
    
    # 配置参数 (可选)
    # 如果需要使用代理访问Zenodo，取消下面的注释并设置代理
    # zenodo_use_proxy = True
    # zenodo_proxy_url = "http://your-proxy-server:port"
    
    # 如果有Figshare API token，取消下面的注释并设置
    # figshare_api_token = "your_figshare_token_here"
    
    # 如果有Zenodo API token，取消下面的注释并设置
    # zenodo_api_token = "your_zenodo_token_here"
    
    collector = FixedIntegratedCollector(
        output_dir="metadata_output",
        # zenodo_use_proxy=zenodo_use_proxy,
        # zenodo_proxy_url=zenodo_proxy_url,
        # figshare_api_token=figshare_api_token,
        # zenodo_api_token=zenodo_api_token
    )
    
    databases_to_collect = [
        'cellxgene',
        'broad_single_cell',
        'zenodo',
        'figshare',
        'dryad'
    ]
    
    logger.info("Starting metadata collection from all databases...")
    print("\nStarting collection process...\n")
    
    results = collector.collect_all(databases=databases_to_collect)
    
    print("\n" + "="*60)
    print("Merging tables from all databases...")
    print("="*60 + "\n")
    merged = collector.merge_tables(results)
    
    print("\n" + "="*60)
    print("Generating summary report...")
    print("="*60 + "\n")
    report = collector.generate_summary_report(results)
    
    print("\n" + "="*60)
    print("Collection Completed!")
    print("="*60)
    print(f"\nResults saved to: {collector.output_dir}")
    print(f"\nTotal collected:")
    print(f"  Series: {report['summary']['total_series']}")
    print(f"  Samples: {report['summary']['total_samples']}")
    print(f"  Successful: {report['summary']['successful_databases']}/{len(databases_to_collect)} databases")
    print("\n" + "="*60)
    
    print("\nDetailed Statistics:")
    print("-"*60)
    for db_name, db_info in report['databases'].items():
        status_icon = "✓" if db_info['status'] == 'success' else "✗"
        print(f"{status_icon} {db_name:20s}: {db_info['series_count']:4d} series, {db_info['sample_count']:4d} samples")
    print("="*60)
    
    # 打印问题和建议
    print("\n" + "="*60)
    print("Notes and Recommendations:")
    print("="*60)
    if report['summary']['failed_databases'] > 0:
        print("\nSome databases failed to collect data. Here are some tips:")
        print()
        if report['databases'].get('zenodo', {}).get('status') == 'failed':
            print("• Zenodo: Connection refused - may require VPN or proxy")
            print("  You can set proxy in the code or use a VPN to access")
        if report['databases'].get('figshare', {}).get('status') == 'failed':
            print("• Figshare: API forbidden - may need API token")
            print("  Get free token at: https://figshare.com/account/applications")
        if report['databases'].get('broad_single_cell', {}).get('status') == 'failed':
            print("• Broad SCP: Server error - service may be temporarily down")
            print("  Try again later or check: https://singlecell.broadinstitute.org")
    print("\n" + "="*60)


if __name__ == "__main__":
    main()