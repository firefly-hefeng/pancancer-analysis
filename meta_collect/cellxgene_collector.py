# cellxgene_collector.py
"""
CELLxGENE 单细胞数据收集器
独立版本 - 修复了数据类型不一致问题
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

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('cellxgene_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def create_session_with_retries(retries=3, backoff_factor=0.5):
    """创建带重试机制的session"""
    session = requests.Session()
    
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=(500, 502, 504),
    )
    
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    
    return session


class CELLxGENECollector:
    """CELLxGENE 数据收集器"""
    
    def __init__(self, output_dir: str = "cellxgene_output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.database_name = "cellxgene"
        self.api_base = "https://api.cellxgene.cziscience.com"
        self.census_api = "https://api.cellxgene.cziscience.com/curation/v1"
        self.session = create_session_with_retries()
        
    def save_raw_data(self, data: any, filename: str):
        """保存原始数据"""
        raw_dir = self.output_dir / "raw"
        raw_dir.mkdir(exist_ok=True)
        
        filepath = raw_dir / filename
        
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
    
    def _safe_join_list(self, data: any) -> str:
        """安全地连接列表数据，处理混合类型（字符串/字典）"""
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
            
            self.save_raw_data(collections, "collections_raw.json")
            logger.info(f"✓ Retrieved {len(collections)} collections")
            return collections
            
        except Exception as e:
            logger.error(f"✗ Error fetching CELLxGENE collections: {e}")
            return []
    
    def fetch_datasets(self) -> List[Dict]:
        """获取所有datasets"""
        url = f"{self.census_api}/datasets"
        logger.info(f"Fetching CELLxGENE datasets from {url}")
        
        try:
            response = self.session.get(url, timeout=120)
            response.raise_for_status()
            datasets = response.json()
            
            self.save_raw_data(datasets, "datasets_raw.json")
            logger.info(f"✓ Retrieved {len(datasets)} datasets")
            return datasets
            
        except requests.exceptions.Timeout:
            logger.warning("Request timed out, trying to extract from collections...")
            return self._fetch_datasets_from_collections()
            
        except Exception as e:
            logger.error(f"✗ Error fetching CELLxGENE datasets: {e}")
            return []
    
    def _fetch_datasets_from_collections(self) -> List[Dict]:
        """从collections中提取datasets信息（备用方法）"""
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
        
        logger.info(f"✓ Extracted {len(all_datasets)} datasets from collections")
        self.save_raw_data(all_datasets, "datasets_from_collections.json")
        return all_datasets
    
    def parse_to_series_table(self, collections: List[Dict], datasets: List[Dict]) -> pd.DataFrame:
        """解析为series表格"""
        series_list = []
        
        # 创建dataset到collection的映射
        dataset_to_collection = {}
        for coll in collections:
            for ds in coll.get('datasets', []):
                ds_id = ds.get('dataset_id') or ds.get('id')
                dataset_to_collection[ds_id] = coll
        
        success_count = 0
        error_count = 0
        
        for dataset in datasets:
            try:
                dataset_id = dataset.get('dataset_id') or dataset.get('id', '')
                
                # 获取collection信息
                collection = dataset_to_collection.get(dataset_id, {})
                if not collection and 'collection_info' in dataset:
                    collection = dataset['collection_info']
                
                # 安全提取各种字段（处理混合类型）
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
                success_count += 1
                
            except Exception as e:
                error_count += 1
                logger.warning(f"Error parsing dataset {dataset.get('dataset_id', dataset.get('id'))}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
        # 保存处理后的数据
        processed_dir = self.output_dir / "processed"
        processed_dir.mkdir(exist_ok=True)
        df.to_csv(processed_dir / "series_table.csv", index=False)
        
        logger.info(f"✓ Created series table: {len(df)} entries (success: {success_count}, errors: {error_count})")
        return df
    
    def parse_to_sample_table(self, datasets: List[Dict]) -> pd.DataFrame:
        """解析为sample表格"""
        sample_list = []
        
        success_count = 0
        error_count = 0
        
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
                success_count += 1
                
            except Exception as e:
                error_count += 1
                logger.warning(f"Error parsing sample from dataset {dataset.get('dataset_id', dataset.get('id'))}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        processed_dir = self.output_dir / "processed"
        processed_dir.mkdir(exist_ok=True)
        df.to_csv(processed_dir / "sample_table.csv", index=False)
        
        logger.info(f"✓ Created sample table: {len(df)} entries (success: {success_count}, errors: {error_count})")
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
        logger.info("="*60)
        logger.info("Starting CELLxGENE collection...")
        logger.info("="*60)
        
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
            
            logger.info("="*60)
            logger.info("CELLxGENE Collection Summary:")
            logger.info(f"  Collections: {len(collections)}")
            logger.info(f"  Datasets: {len(datasets)}")
            logger.info(f"  Series entries: {len(series_df)}")
            logger.info(f"  Sample entries: {len(sample_df)}")
            logger.info(f"  Output directory: {self.output_dir}")
            logger.info("="*60)
            
            return {
                'series': series_df,
                'samples': sample_df,
                'collections': collections,
                'datasets': datasets
            }
        else:
            logger.error("✗ Failed to collect CELLxGENE data")
            return None


def main():
    """测试主函数"""
    print("\n" + "="*60)
    print("CELLxGENE Data Collector - Standalone Version")
    print("="*60 + "\n")
    
    collector = CELLxGENECollector(output_dir="cellxgene_output")
    result = collector.collect()
    
    if result:
        print("\n✓ Collection successful!")
        print(f"\nFiles saved to: {collector.output_dir}")
        print("\nGenerated files:")
        print("  - raw/collections_raw.json")
        print("  - raw/datasets_raw.json")
        print("  - processed/series_table.csv")
        print("  - processed/sample_table.csv")
    else:
        print("\n✗ Collection failed!")
    
    print("\n" + "="*60 + "\n")


if __name__ == "__main__":
    main()