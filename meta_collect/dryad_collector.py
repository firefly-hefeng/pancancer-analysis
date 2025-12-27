# dryad_collector.py
"""
DRYAD 数据收集器
独立版本
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
        logging.FileHandler('dryad_collection.log'),
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


class DRYADCollector:
    """DRYAD 数据收集器"""
    
    def __init__(self, output_dir: str = "dryad_output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.database_name = "dryad"
        self.api_base = "https://datadryad.org/api/v2"
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
    
    def search_datasets(self, query: str, page: int = 1, per_page: int = 100) -> Dict:
        """搜索DRYAD数据集"""
        url = f"{self.api_base}/search"
        
        params = {
            'q': query,
            'page': page,
            'per_page': per_page
        }
        
        logger.info(f"Searching DRYAD (page {page}, query: '{query[:50]}...')")
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
            
        except Exception as e:
            logger.error(f"  ✗ Error searching DRYAD: {e}")
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
                    logger.info(f"  No more datasets for query '{query}' at page {page}")
                    break
                
                all_datasets.extend(datasets)
                logger.info(f"  ✓ Query '{query}' page {page}: found {len(datasets)} datasets")
                
                # 检查是否还有更多页
                if not data.get('_links', {}).get('next'):
                    break
                
                page += 1
                time.sleep(1)
        
        # 去重
        unique_datasets = {ds['id']: ds for ds in all_datasets}.values()
        unique_datasets = list(unique_datasets)
        
        self.save_raw_data(unique_datasets, "datasets_raw.json")
        
        logger.info(f"✓ Retrieved {len(unique_datasets)} unique datasets from DRYAD")
        return unique_datasets
    
    def parse_to_series_table(self, datasets: List[Dict]) -> pd.DataFrame:
        """解析为series表格"""
        series_list = []
        
        success_count = 0
        error_count = 0
        
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
                success_count += 1
                
            except Exception as e:
                error_count += 1
                logger.warning(f"Error parsing DRYAD dataset {dataset.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
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
                success_count += 1
                
            except Exception as e:
                error_count += 1
                logger.warning(f"Error parsing sample from DRYAD dataset {dataset.get('id')}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        processed_dir = self.output_dir / "processed"
        processed_dir.mkdir(exist_ok=True)
        df.to_csv(processed_dir / "sample_table.csv", index=False)
        
        logger.info(f"✓ Created sample table: {len(df)} entries (success: {success_count}, errors: {error_count})")
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
        logger.info("="*60)
        logger.info("Starting DRYAD collection...")
        logger.info("="*60)
        
        datasets = self.fetch_all_scrna_datasets()
        
        if datasets:
            series_df = self.parse_to_series_table(datasets)
            sample_df = self.parse_to_sample_table(datasets)
            
            logger.info("="*60)
            logger.info("DRYAD Collection Summary:")
            logger.info(f"  Datasets: {len(datasets)}")
            logger.info(f"  Series entries: {len(series_df)}")
            logger.info(f"  Sample entries: {len(sample_df)}")
            logger.info(f"  Output directory: {self.output_dir}")
            logger.info("="*60)
            
            return {
                'series': series_df,
                'samples': sample_df,
                'datasets': datasets
            }
        else:
            logger.error("✗ Failed to collect DRYAD data")
            return None


def main():
    """测试主函数"""
    print("\n" + "="*60)
    print("DRYAD Data Collector - Standalone Version")
    print("="*60 + "\n")
    
    collector = DRYADCollector(output_dir="dryad_output")
    result = collector.collect()
    
    if result:
        print("\n✓ Collection successful!")
        print(f"\nFiles saved to: {collector.output_dir}")
        print("\nGenerated files:")
        print("  - raw/datasets_raw.json")
        print("  - processed/series_table.csv")
        print("  - processed/sample_table.csv")
    else:
        print("\n✗ Collection failed - see log for details")
    
    print("\n" + "="*60 + "\n")


if __name__ == "__main__":
    main()