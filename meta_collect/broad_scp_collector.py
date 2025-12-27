# broad_scp_collector.py
"""
Broad Single Cell Portal 数据收集器
独立版本 - 包含多种备用访问方法
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
        logging.FileHandler('broad_scp_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def create_session_with_retries(retries=5, backoff_factor=1.0, verify_ssl=False):
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
    session.verify = verify_ssl
    
    return session


class BroadSCPCollector:
    """Broad Single Cell Portal 数据收集器"""
    
    def __init__(self, output_dir: str = "broad_scp_output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.database_name = "broad_single_cell"
        self.web_base = "https://singlecell.broadinstitute.org"
        self.api_base = "https://singlecell.broadinstitute.org/single_cell/api/v1"
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
    
    def fetch_studies_alternative(self) -> List[Dict]:
        """使用替代方法获取研究列表"""
        logger.info("Trying alternative methods for Broad SCP...")
        
        # 尝试多个可能的endpoint
        endpoints = [
            f"{self.api_base}/studies",
            f"{self.api_base}/site/studies",
            f"{self.api_base}/search",
            f"{self.web_base}/single_cell/data/public"
        ]
        
        for endpoint in endpoints:
            try:
                logger.info(f"  Trying endpoint: {endpoint}")
                response = self.session.get(endpoint, timeout=30)
                
                if response.status_code == 200:
                    content_type = response.headers.get('content-type', '')
                    
                    if 'json' in content_type:
                        data = response.json()
                        if data:
                            studies = data if isinstance(data, list) else data.get('studies', [])
                            if studies:
                                logger.info(f"  ✓ Success with endpoint: {endpoint}")
                                return studies
                    else:
                        logger.debug(f"  Non-JSON response from {endpoint}")
                        
            except Exception as e:
                logger.debug(f"  ✗ Failed endpoint {endpoint}: {e}")
                continue
        
        logger.warning("All alternative endpoints failed")
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
                logger.warning("  Server error (500), trying alternative method...")
                alt_studies = self.fetch_studies_alternative()
                return {'studies': alt_studies} if alt_studies else {}
            else:
                logger.error(f"  HTTP error: {e}")
                return {}
            
        except requests.exceptions.Timeout:
            logger.warning("  Request timeout, trying alternative method...")
            alt_studies = self.fetch_studies_alternative()
            return {'studies': alt_studies} if alt_studies else {}
            
        except Exception as e:
            logger.error(f"  Error fetching Broad SCP studies: {e}")
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
                    logger.info(f"  No more studies at page {page}")
                    break
                
                # 过滤人类数据
                human_studies = [
                    s for s in studies 
                    if any('human' in sp.lower() or 'homo sapiens' in sp.lower() 
                          for sp in s.get('species', []))
                ]
                
                all_studies.extend(human_studies)
                logger.info(f"  ✓ Page {page}: found {len(human_studies)} human studies (total: {len(studies)})")
                
                if len(studies) < 100:
                    break
                
                page += 1
                time.sleep(2)  # 避免请求过快
                
            except Exception as e:
                logger.error(f"  ✗ Error on page {page}: {e}")
                break
        
        self.save_raw_data(all_studies, "studies_raw.json")
        logger.info(f"✓ Retrieved {len(all_studies)} human studies from Broad SCP")
        return all_studies
    
    def parse_to_series_table(self, studies: List[Dict]) -> pd.DataFrame:
        """解析为series表格"""
        series_list = []
        
        success_count = 0
        error_count = 0
        
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
                success_count += 1
                
            except Exception as e:
                error_count += 1
                logger.warning(f"Error parsing study {study.get('accession')}: {e}")
                continue
        
        df = pd.DataFrame(series_list)
        
        processed_dir = self.output_dir / "processed"
        processed_dir.mkdir(exist_ok=True)
        df.to_csv(processed_dir / "series_table.csv", index=False)
        
        logger.info(f"✓ Created series table: {len(df)} entries (success: {success_count}, errors: {error_count})")
        return df
    
    def parse_to_sample_table(self, studies: List[Dict]) -> pd.DataFrame:
        """解析为sample表格"""
        sample_list = []
        
        success_count = 0
        error_count = 0
        
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
                success_count += 1
                
            except Exception as e:
                error_count += 1
                logger.warning(f"Error parsing sample from study {study.get('accession')}: {e}")
                continue
        
        df = pd.DataFrame(sample_list)
        
        processed_dir = self.output_dir / "processed"
        processed_dir.mkdir(exist_ok=True)
        df.to_csv(processed_dir / "sample_table.csv", index=False)
        
        logger.info(f"✓ Created sample table: {len(df)} entries (success: {success_count}, errors: {error_count})")
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
        logger.info("="*60)
        logger.info("Starting Broad Single Cell Portal collection...")
        logger.info("="*60)
        
        studies = self.fetch_all_studies()
        
        if studies:
            series_df = self.parse_to_series_table(studies)
            sample_df = self.parse_to_sample_table(studies)
            
            logger.info("="*60)
            logger.info("Broad SCP Collection Summary:")
            logger.info(f"  Studies: {len(studies)}")
            logger.info(f"  Series entries: {len(series_df)}")
            logger.info(f"  Sample entries: {len(sample_df)}")
            logger.info(f"  Output directory: {self.output_dir}")
            logger.info("="*60)
            
            return {
                'series': series_df,
                'samples': sample_df,
                'studies': studies
            }
        else:
            logger.warning("="*60)
            logger.warning("Broad SCP Collection Failed")
            logger.warning("Possible reasons:")
            logger.warning("  1. Service temporarily unavailable (500 error)")
            logger.warning("  2. API endpoint changed")
            logger.warning("  3. Network connectivity issues")
            logger.warning("Recommendation: Try again later or check service status at:")
            logger.warning("  https://singlecell.broadinstitute.org")
            logger.warning("="*60)
            return None


def main():
    """测试主函数"""
    print("\n" + "="*60)
    print("Broad Single Cell Portal Collector - Standalone Version")
    print("="*60 + "\n")
    
    collector = BroadSCPCollector(output_dir="broad_scp_output")
    result = collector.collect()
    
    if result:
        print("\n✓ Collection successful!")
        print(f"\nFiles saved to: {collector.output_dir}")
        print("\nGenerated files:")
        print("  - raw/studies_raw.json")
        print("  - processed/series_table.csv")
        print("  - processed/sample_table.csv")
    else:
        print("\n✗ Collection failed - see log for details")
    
    print("\n" + "="*60 + "\n")


if __name__ == "__main__":
    main()