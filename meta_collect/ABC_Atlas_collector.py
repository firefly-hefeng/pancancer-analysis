"""
Allen Brain Cell Atlas (ABC Atlas) Metadata 全面收集工具 - 修复版
支持多种数据源和访问方式的系统化收集
"""

import requests
import pandas as pd
import numpy as np
from datetime import datetime
import json
import time
from typing import Dict, List, Optional, Tuple
import xml.etree.ElementTree as ET
from urllib.parse import urlencode, quote
import logging
from pathlib import Path
import re
from bs4 import BeautifulSoup
import warnings
warnings.filterwarnings('ignore')

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('abc_atlas_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class ABCAtlasCollector:
    """Allen Brain Cell Atlas 数据收集器"""
    
    def __init__(self, output_dir: str = "./abc_atlas_metadata"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Allen Brain Map API endpoints
        self.api_base = "https://api.brain-map.org/api/v2"
        self.allen_portal = "https://portal.brain-map.org"
        self.knowledge_base = "https://knowledge.brain-map.org"
        
        # CELLxGENE API (Allen 数据也发布在这里)
        self.cellxgene_api = "https://api.cellxgene.cziscience.com"
        
        # 数据存储
        self.raw_data = {}
        self.processed_data = []
        
        logger.info(f"输出目录: {self.output_dir}")
    
    def save_raw_data(self, data: any, filename: str, format: str = 'json'):
        """保存原始数据"""
        if format == 'json':
            filepath = self.output_dir / f"raw_{filename}.json"
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        elif format == 'csv':
            filepath = self.output_dir / f"raw_{filename}.csv"
            if isinstance(data, pd.DataFrame):
                data.to_csv(filepath, index=False, encoding='utf-8-sig')
            else:
                pd.DataFrame(data).to_csv(filepath, index=False, encoding='utf-8-sig')
        
        logger.info(f"✓ 原始数据已保存: {filepath}")
        return filepath
    
    def save_processed_data(self, df: pd.DataFrame, filename: str):
        """保存处理后的数据"""
        filepath = self.output_dir / f"processed_{filename}.csv"
        df.to_csv(filepath, index=False, encoding='utf-8-sig')
        logger.info(f"✓ 处理后数据已保存: {filepath}")
        return filepath
    
    # ========== 方法 1: Allen Brain Map API ==========
    
    def fetch_allen_datasets_via_api(self) -> List[Dict]:
        """通过 Allen Brain Map API 获取数据集"""
        logger.info("\n" + "="*80)
        logger.info("方法 1: 通过 Allen Brain Map API 收集数据")
        logger.info("="*80)
        
        all_datasets = []
        
        # 1. 获取所有项目
        projects = self._fetch_allen_projects()
        if projects:
            all_datasets.extend([{'type': 'project', 'data': p} for p in projects])
        
        # 2. 获取数据集列表
        datasets = self._fetch_allen_data_products()
        if datasets:
            all_datasets.extend([{'type': 'data_product', 'data': d} for d in datasets])
        
        # 3. 获取实验数据
        experiments = self._fetch_allen_experiments()
        if experiments:
            all_datasets.extend([{'type': 'experiment', 'data': e} for e in experiments])
        
        self.raw_data['allen_api'] = all_datasets
        self.save_raw_data(all_datasets, 'allen_api_all_data')
        
        logger.info(f"✓ Allen API: 收集到 {len(all_datasets)} 条记录")
        return all_datasets
    
    def _fetch_allen_projects(self) -> List[Dict]:
        """获取 Allen 项目列表"""
        logger.info("→ 获取 Allen 项目列表...")
        
        url = f"{self.api_base}/data/query.json"
        criteria = "model::Project,rma::criteria,[name$il'*single*cell*'],rma::options[num_rows$eq2000]"
        params = {'criteria': criteria}
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            result = response.json()
            projects = result.get('msg', [])
            logger.info(f"  找到 {len(projects)} 个项目")
            return projects
        except Exception as e:
            logger.error(f"  获取项目失败: {e}")
            return []
    
    def _fetch_allen_data_products(self) -> List[Dict]:
        """获取 Allen 数据产品"""
        logger.info("→ 获取 Allen 数据产品...")
        
        url = f"{self.api_base}/data/query.json"
        criteria = "model::WellKnownFile,rma::criteria,well_known_file_type[name$il'*RNAseq*'],rma::options[num_rows$eq5000]"
        params = {'criteria': criteria}
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            result = response.json()
            files = result.get('msg', [])
            logger.info(f"  找到 {len(files)} 个数据文件")
            return files
        except Exception as e:
            logger.error(f"  获取数据产品失败: {e}")
            return []
    
    def _fetch_allen_experiments(self) -> List[Dict]:
        """获取实验数据"""
        logger.info("→ 获取实验数据...")
        
        url = f"{self.api_base}/data/query.json"
        criteria = "model::RnaSeqExperiment,rma::options[num_rows$eq5000]"
        params = {'criteria': criteria}
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            result = response.json()
            experiments = result.get('msg', [])
            logger.info(f"  找到 {len(experiments)} 个实验")
            return experiments
        except Exception as e:
            logger.error(f"  获取实验数据失败: {e}")
            return []
    
    # ========== 方法 2: CELLxGENE (Allen 数据) ==========
    
    def fetch_allen_from_cellxgene(self) -> List[Dict]:
        """从 CELLxGENE 获取 Allen Institute 发布的数据"""
        logger.info("\n" + "="*80)
        logger.info("方法 2: 从 CELLxGENE 收集 Allen Institute 数据")
        logger.info("="*80)
        
        url = f"{self.cellxgene_api}/curation/v1/collections"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            all_collections = response.json()
            
            # 筛选 Allen Institute 的数据
            allen_collections = []
            for collection in all_collections:
                contact_info = str(collection.get('contact_name', '')) + str(collection.get('contact_email', ''))
                consortia = collection.get('consortia', [])
                publisher = collection.get('publisher', '')
                
                if any([
                    'allen' in contact_info.lower(),
                    'allen' in publisher.lower(),
                    any('allen' in str(c).lower() for c in consortia),
                    'brain' in str(collection.get('name', '')).lower() and 'cell' in str(collection.get('name', '')).lower()
                ]):
                    allen_collections.append(collection)
                    
                    # 获取该 collection 的详细信息
                    collection_id = collection.get('id')
                    if collection_id:
                        detail = self._fetch_cellxgene_collection_detail(collection_id)
                        if detail:
                            allen_collections[-1]['datasets'] = detail.get('datasets', [])
            
            self.raw_data['cellxgene_allen'] = allen_collections
            self.save_raw_data(allen_collections, 'cellxgene_allen_collections')
            
            logger.info(f"✓ CELLxGENE: 找到 {len(allen_collections)} 个 Allen 相关的 collections")
            
            # 展开所有 datasets
            all_datasets = []
            for coll in allen_collections:
                for ds in coll.get('datasets', []):
                    ds['collection_info'] = {k: v for k, v in coll.items() if k != 'datasets'}
                    all_datasets.append(ds)
            
            logger.info(f"✓ CELLxGENE: 展开为 {len(all_datasets)} 个 datasets")
            
            return all_datasets
            
        except Exception as e:
            logger.error(f"从 CELLxGENE 获取数据失败: {e}")
            return []
    
    def _fetch_cellxgene_collection_detail(self, collection_id: str) -> Optional[Dict]:
        """获取 CELLxGENE collection 详细信息"""
        url = f"{self.cellxgene_api}/curation/v1/collections/{collection_id}"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error(f"  获取 collection {collection_id} 详情失败: {e}")
            return None
    
    # ========== 方法 3: ABC Atlas 专用数据 ==========
    
    def fetch_abc_atlas_data(self) -> List[Dict]:
        """获取 ABC Atlas 特定数据"""
        logger.info("\n" + "="*80)
        logger.info("方法 3: 收集 ABC Atlas 专用数据")
        logger.info("="*80)
        
        abc_data = []
        
        # ABC Atlas 主要数据集（手动整理的已知数据集）
        known_datasets = self._get_known_abc_datasets()
        abc_data.extend(known_datasets)
        
        # 尝试从 Knowledge Base 抓取
        kb_data = self._scrape_knowledge_base()
        abc_data.extend(kb_data)
        
        self.raw_data['abc_atlas'] = abc_data
        self.save_raw_data(abc_data, 'abc_atlas_specific_data')
        
        logger.info(f"✓ ABC Atlas: 收集到 {len(abc_data)} 条记录")
        
        return abc_data
    
    def _get_known_abc_datasets(self) -> List[Dict]:
        """获取已知的 ABC Atlas 核心数据集"""
        logger.info("→ 加载已知的 ABC Atlas 数据集...")
        
        known_datasets = [
            {
                'id': 'ABC_Atlas_Human_MTG',
                'title': 'A multimodal cell census and atlas of the mammalian primary motor cortex',
                'tissue': 'primary motor cortex',
                'organism': 'Homo sapiens',
                'pubmed': '34795446',
                'doi': '10.1038/s41586-021-03950-0',
                'publication_date': '2021-10',
                'modality': 'scRNA-seq, snRNA-seq, MERFISH, Patch-seq',
                'cell_count': '~450,000',
                'data_link': 'https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x',
                'description': 'Comprehensive cellular census of human and mouse primary motor cortex'
            },
            {
                'id': 'ABC_Atlas_Human_Brain_Full',
                'title': 'A comparative atlas of single-cell chromatin accessibility in the human brain',
                'tissue': 'multiple brain regions',
                'organism': 'Homo sapiens',
                'pubmed': '38092916',
                'doi': '10.1126/science.adf7044',
                'publication_date': '2023-12',
                'modality': 'snATAC-seq, snRNA-seq',
                'cell_count': '~3,000,000',
                'data_link': 'https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas',
                'description': 'Comprehensive human brain cell atlas across development and adulthood'
            },
            {
                'id': 'ABC_Atlas_WHB_10Xv3',
                'title': 'Whole Human Brain scRNA-seq (10X Chromium v3)',
                'tissue': 'whole brain',
                'organism': 'Homo sapiens',
                'publication_date': '2023',
                'modality': 'scRNA-seq (10X Chromium v3)',
                'cell_count': '>3,000,000',
                'data_link': 'https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html',
                'description': 'Comprehensive single-cell transcriptomic atlas of the entire human brain'
            },
            {
                'id': 'ABC_Atlas_WHB_10Xv3_SMART',
                'title': 'Whole Human Brain scRNA-seq (10X v3 + SMART-Seq)',
                'tissue': 'whole brain - multiple regions',
                'organism': 'Homo sapiens',
                'publication_date': '2023',
                'modality': 'scRNA-seq (10X v3), SMART-Seq v4',
                'cell_count': '>3,000,000',
                'data_link': 'https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas',
                'description': 'Multi-platform transcriptomic profiling of human brain'
            },
            {
                'id': 'ABC_Atlas_Human_Development',
                'title': 'Developmental atlas of the human brain',
                'tissue': 'brain - developmental stages',
                'organism': 'Homo sapiens',
                'publication_date': '2022-2023',
                'modality': 'snRNA-seq, snATAC-seq',
                'development_stage': 'fetal to adult',
                'data_link': 'https://knowledge.brain-map.org/data',
                'description': 'Temporal dynamics of cell types during human brain development'
            }
        ]
        
        logger.info(f"  加载了 {len(known_datasets)} 个已知数据集")
        
        return known_datasets
    
    def _scrape_knowledge_base(self) -> List[Dict]:
        """从 Allen Knowledge Base 抓取数据信息"""
        logger.info("→ 从 Knowledge Base 抓取数据...")
        
        kb_data = []
        
        try:
            urls_to_check = [
                'https://knowledge.brain-map.org/data',
                'https://knowledge.brain-map.org/data/HT29V88EW95CW1FJ6B2/summary',
            ]
            
            for url in urls_to_check:
                try:
                    response = requests.get(url, timeout=30)
                    response.raise_for_status()
                    
                    soup = BeautifulSoup(response.text, 'html.parser')
                    titles = soup.find_all(['h1', 'h2', 'h3'])
                    
                    logger.info(f"  从 {url} 找到 {len(titles)} 个潜在标题")
                    
                except Exception as e:
                    logger.warning(f"  无法抓取 {url}: {e}")
                    continue
            
        except Exception as e:
            logger.error(f"  Knowledge Base 抓取失败: {e}")
        
        return kb_data
    
    # ========== 方法 4: S3 数据目录 ==========
    
    def fetch_allen_s3_manifest(self) -> List[Dict]:
        """获取 Allen S3 存储的数据清单"""
        logger.info("\n" + "="*80)
        logger.info("方法 4: 收集 Allen S3 数据目录")
        logger.info("="*80)
        
        s3_base = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
        s3_data = []
        
        try:
            manifest_urls = [
                f"{s3_base}/index.html",
                f"{s3_base}/metadata/manifest.json",
                f"{s3_base}/METADATA.json",
            ]
            
            for url in manifest_urls:
                try:
                    response = requests.get(url, timeout=30)
                    if response.status_code == 200:
                        logger.info(f"  成功访问: {url}")
                        
                        try:
                            data = response.json()
                            s3_data.append({'source': url, 'data': data})
                        except:
                            s3_data.append({
                                'source': url,
                                'content_type': response.headers.get('content-type'),
                                'size': len(response.content)
                            })
                except:
                    continue
            
            self.raw_data['s3_manifest'] = s3_data
            self.save_raw_data(s3_data, 's3_manifest')
            
            logger.info(f"✓ S3: 收集到 {len(s3_data)} 条清单信息")
            
        except Exception as e:
            logger.error(f"S3 清单收集失败: {e}")
        
        return s3_data
    
    # ========== 数据整合与标准化 ==========
    
    def integrate_all_sources(self) -> pd.DataFrame:
        """整合所有数据源"""
        logger.info("\n" + "="*80)
        logger.info("整合所有数据源")
        logger.info("="*80)
        
        all_records = []
        
        # 1. 处理 Allen API 数据
        if 'allen_api' in self.raw_data:
            logger.info("→ 处理 Allen API 数据...")
            for item in self.raw_data['allen_api']:
                # 修复: 检查 item 是否为字典，并提取实际数据
                if isinstance(item, dict):
                    # 如果有 'data' 字段，提取它
                    actual_data = item.get('data', item)
                    record_type = item.get('type', 'unknown')
                    
                    if isinstance(actual_data, dict):
                        record = self._transform_allen_api_to_standard(actual_data, record_type)
                        if record:
                            all_records.append(record)
        
        # 2. 处理 CELLxGENE 数据
        if 'cellxgene_allen' in self.raw_data:
            logger.info("→ 处理 CELLxGENE 数据...")
            for coll in self.raw_data['cellxgene_allen']:
                if isinstance(coll, dict):
                    for ds in coll.get('datasets', []):
                        if isinstance(ds, dict):
                            record = self._transform_cellxgene_to_standard(ds, coll)
                            if record:
                                all_records.append(record)
        
        # 3. 处理 ABC Atlas 专用数据
        if 'abc_atlas' in self.raw_data:
            logger.info("→ 处理 ABC Atlas 数据...")
            for item in self.raw_data['abc_atlas']:
                if isinstance(item, dict):
                    record = self._transform_abc_atlas_to_standard(item)
                    if record:
                        all_records.append(record)
        
        # 创建 DataFrame
        df = pd.DataFrame(all_records)
        
        # 去重
        if not df.empty:
            initial_count = len(df)
            df = df.drop_duplicates(subset=['id'], keep='first')
            logger.info(f"  去重: {initial_count} → {len(df)} 条记录")
        
        # 保存整合数据
        self.save_processed_data(df, 'integrated_all_sources')
        
        logger.info(f"✓ 整合完成: 共 {len(df)} 条记录")
        
        return df
    
    def _transform_allen_api_to_standard(self, item: Dict, record_type: str = 'unknown') -> Optional[Dict]:
        """转换 Allen API 数据为标准格式"""
        if not item or not isinstance(item, dict):
            return None
        
        record = {
            'id': f"ALLEN_{record_type}_{item.get('id', '')}",
            'title': item.get('name', item.get('description', '')),
            'disease_general': 'normal',
            'disease': 'healthy control',
            'pubmed': '',
            'source_database': f'Allen Brain Map API ({record_type})',
            'access_link': f"{self.allen_portal}/explore/rnaseq/{item.get('id', '')}",
            'open_status': 'public',
            'ethnicity': '',
            'sex': '',
            'tissue': 'brain',
            'sequencing_platform': '',
            'experiment_design': item.get('type', record_type),
            'sample_type': 'single cell',
            'summary': item.get('description', ''),
            'citation_count': '',
            'publication_date': '',
            'submission_date': '',
            'last_update_date': '',
            'contact_name': 'Allen Institute for Brain Science',
            'contact_email': 'support@brain-map.org',
            'contact_institute': 'Allen Institute',
            'data_tier': 'processed_matrix',
            'tissue_location': 'brain',
            'supplementary_information': json.dumps(item, ensure_ascii=False),
        }
        
        return record
    
    def _transform_cellxgene_to_standard(self, dataset: Dict, collection: Dict) -> Optional[Dict]:
        """转换 CELLxGENE 数据为标准格式"""
        if not isinstance(dataset, dict) or not isinstance(collection, dict):
            return None
        
        record = {
            'id': dataset.get('id', ''),
            'title': dataset.get('title', collection.get('name', '')),
            'disease_general': self._extract_disease_general_cellxgene(dataset),
            'disease': self._extract_diseases_cellxgene(dataset),
            'pubmed': self._extract_pubmed_cellxgene(collection),
            'source_database': 'CELLxGENE (Allen Institute)',
            'access_link': f"https://cellxgene.cziscience.com/e/{dataset.get('id', '')}.cxg/",
            'open_status': 'public',
            'ethnicity': self._extract_field_cellxgene(dataset, 'self_reported_ethnicity'),
            'sex': self._extract_field_cellxgene(dataset, 'sex'),
            'tissue': self._extract_field_cellxgene(dataset, 'tissue'),
            'sequencing_platform': '',
            'experiment_design': self._extract_field_cellxgene(dataset, 'assay'),
            'sample_type': self._extract_field_cellxgene(dataset, 'suspension_type'),
            'summary': collection.get('description', ''),
            'citation_count': '',
            'publication_date': collection.get('published_at', ''),
            'submission_date': collection.get('created_at', ''),
            'last_update_date': collection.get('revised_at', ''),
            'contact_name': collection.get('contact_name', ''),
            'contact_email': collection.get('contact_email', ''),
            'contact_institute': 'Allen Institute for Brain Science',
            'data_tier': 'processed_matrix',
            'tissue_location': 'brain',
            'supplementary_information': f"Cell count: {dataset.get('cell_count', '')}",
            'cell_count': dataset.get('cell_count', ''),
            'organism': self._extract_field_cellxgene(dataset, 'organism'),
        }
        
        return record
    
    def _transform_abc_atlas_to_standard(self, item: Dict) -> Optional[Dict]:
        """转换 ABC Atlas 数据为标准格式"""
        if not isinstance(item, dict):
            return None
        
        record = {
            'id': item.get('id', ''),
            'title': item.get('title', ''),
            'disease_general': 'normal',
            'disease': 'healthy control',
            'pubmed': item.get('pubmed', ''),
            'source_database': 'ABC Atlas',
            'access_link': item.get('data_link', ''),
            'open_status': 'public',
            'ethnicity': '',
            'sex': '',
            'tissue': item.get('tissue', 'brain'),
            'sequencing_platform': self._extract_platform_from_modality(item.get('modality', '')),
            'experiment_design': item.get('modality', ''),
            'sample_type': 'single cell/nucleus',
            'summary': item.get('description', ''),
            'citation_count': '',
            'publication_date': item.get('publication_date', ''),
            'submission_date': '',
            'last_update_date': '',
            'contact_name': 'Allen Institute for Brain Science',
            'contact_email': 'celltypes@alleninstitute.org',
            'contact_institute': 'Allen Institute',
            'data_tier': 'processed_matrix',
            'tissue_location': 'brain',
            'supplementary_information': f"DOI: {item.get('doi', '')}; Cell count: {item.get('cell_count', '')}",
            'cell_count': item.get('cell_count', ''),
            'doi': item.get('doi', ''),
            'organism': item.get('organism', 'Homo sapiens'),
            'development_stage': item.get('development_stage', ''),
        }
        
        return record
    
    # ========== 辅助函数 ==========
    
    def _extract_disease_general_cellxgene(self, ds: Dict) -> str:
        """提取疾病大类"""
        diseases = ds.get('disease', [])
        if not diseases:
            return 'normal'
        
        disease_labels = [d.get('label', '').lower() for d in diseases if isinstance(d, dict)]
        
        if any('normal' in d or 'healthy' in d for d in disease_labels):
            return 'normal'
        elif any('alzheimer' in d or 'parkinson' in d for d in disease_labels):
            return 'neurodegenerative'
        else:
            return 'other_disease'
    
    def _extract_diseases_cellxgene(self, ds: Dict) -> str:
        """提取具体疾病"""
        diseases = ds.get('disease', [])
        if isinstance(diseases, list):
            disease_labels = [d.get('label', '') for d in diseases if isinstance(d, dict)]
            return '; '.join(disease_labels)
        return ''
    
    def _extract_pubmed_cellxgene(self, coll: Dict) -> str:
        """提取 PubMed ID"""
        links = coll.get('links', [])
        if isinstance(links, list):
            for link in links:
                if isinstance(link, dict):
                    url = link.get('link_url', '')
                    if 'pubmed' in url.lower():
                        parts = url.split('/')
                        if parts:
                            return parts[-1]
        return ''
    
    def _extract_field_cellxgene(self, ds: Dict, field: str) -> str:
        """提取 CELLxGENE 字段"""
        data = ds.get(field, [])
        if isinstance(data, list):
            labels = [d.get('label', '') for d in data if isinstance(d, dict)]
            return '; '.join(labels)
        return ''
    
    def _extract_platform_from_modality(self, modality: str) -> str:
        """从技术平台提取测序平台"""
        if '10X' in modality or '10x' in modality:
            return '10X Genomics'
        elif 'SMART' in modality:
            return 'SMART-Seq'
        elif 'MERFISH' in modality:
            return 'MERFISH'
        return ''
    
    # ========== 数据质量检查 ==========
    
    def quality_check(self, df: pd.DataFrame) -> pd.DataFrame:
        """数据质量检查和统计"""
        logger.info("\n" + "="*80)
        logger.info("数据质量检查")
        logger.info("="*80)
        
        # 1. 完整性检查
        logger.info("\n→ 字段完整性:")
        completeness = {}
        for col in df.columns:
            non_empty = df[col].notna().sum()
            completeness[col] = f"{non_empty}/{len(df)} ({non_empty/len(df)*100:.1f}%)"
        
        completeness_df = pd.DataFrame([completeness]).T
        completeness_df.columns = ['完整性']
        print(completeness_df)
        
        # 2. 数据统计
        logger.info("\n→ 数据来源分布:")
        if 'source_database' in df.columns:
            print(df['source_database'].value_counts())
        
        logger.info("\n→ 组织分布:")
        if 'tissue' in df.columns:
            tissue_dist = df['tissue'].value_counts().head(10)
            print(tissue_dist)
        
        logger.info("\n→ 疾病类型分布:")
        if 'disease_general' in df.columns:
            print(df['disease_general'].value_counts())
        
        logger.info("\n→ 实验设计分布:")
        if 'experiment_design' in df.columns:
            design_dist = df['experiment_design'].value_counts().head(10)
            print(design_dist)
        
        # 保存质量报告
        quality_report = {
            'total_records': len(df),
            'completeness': completeness,
            'source_distribution': df['source_database'].value_counts().to_dict() if 'source_database' in df.columns else {},
            'tissue_distribution': df['tissue'].value_counts().head(20).to_dict() if 'tissue' in df.columns else {},
        }
        
        self.save_raw_data(quality_report, 'quality_report')
        
        return df
    
    # ========== 主收集流程 ==========
    
    def collect_all_metadata(self) -> pd.DataFrame:
        """执行完整的元数据收集流程"""
        
        logger.info("\n" + "="*80)
        logger.info("Allen Brain Cell Atlas - 完整元数据收集")
        logger.info("="*80)
        
        # 1. Allen Brain Map API
        try:
            allen_api_data = self.fetch_allen_datasets_via_api()
        except Exception as e:
            logger.error(f"Allen API 收集失败: {e}")
            allen_api_data = []
        
        # 2. CELLxGENE (Allen 数据)
        try:
            cellxgene_data = self.fetch_allen_from_cellxgene()
        except Exception as e:
            logger.error(f"CELLxGENE 收集失败: {e}")
            cellxgene_data = []
        
        # 3. ABC Atlas 专用数据
        try:
            abc_data = self.fetch_abc_atlas_data()
        except Exception as e:
            logger.error(f"ABC Atlas 收集失败: {e}")
            abc_data = []
        
        # 4. S3 清单
        try:
            s3_data = self.fetch_allen_s3_manifest()
        except Exception as e:
            logger.error(f"S3 清单收集失败: {e}")
            s3_data = []
        
        # 5. 整合所有数据源
        integrated_df = self.integrate_all_sources()
        
        # 6. 数据质量检查
        final_df = self.quality_check(integrated_df)
        
        # 7. 保存最终结果
        final_path = self.save_processed_data(final_df, 'FINAL_abc_atlas_metadata')
        
        logger.info("\n" + "="*80)
        logger.info(f"✓ 收集完成!")
        logger.info(f"✓ 总计: {len(final_df)} 条记录")
        logger.info(f"✓ 最终文件: {final_path}")
        logger.info("="*80)
        
        return final_df


# ========== 补充: 引用计数获取 ==========

class CitationEnricher:
    """为数据集添加引用计数"""
    
    def __init__(self):
        self.semantic_scholar_api = "https://api.semanticscholar.org/graph/v1"
        
    def get_citation_count(self, pubmed_id: str = None, doi: str = None) -> int:
        """获取引用次数"""
        if not pubmed_id and not doi:
            return 0
        
        try:
            if pubmed_id:
                url = f"{self.semantic_scholar_api}/paper/PMID:{pubmed_id}"
            elif doi:
                url = f"{self.semantic_scholar_api}/paper/DOI:{doi}"
            else:
                return 0
            
            params = {'fields': 'citationCount'}
            response = requests.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                return data.get('citationCount', 0)
            
        except Exception as e:
            logger.debug(f"获取引用失败: {e}")
        
        return 0
    
    def enrich_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """为 DataFrame 添加引用计数"""
        logger.info("添加引用计数信息...")
        
        citation_counts = []
        
        for idx, row in df.iterrows():
            if idx % 10 == 0:
                logger.info(f"  处理进度: {idx+1}/{len(df)}")
            
            pubmed = row.get('pubmed', '')
            doi = row.get('doi', '')
            
            count = self.get_citation_count(pubmed_id=pubmed, doi=doi)
            citation_counts.append(count)
            
            time.sleep(0.3)  # 避免 API 限流
        
        df['citation_count'] = citation_counts
        
        logger.info(f"✓ 引用计数添加完成")
        
        return df


# ========== 主程序 ==========

def main():
    """主程序入口"""
    
    print("\n" + "="*80)
    print("Allen Brain Cell Atlas (ABC Atlas) - 元数据全面收集工具")
    print("="*80)
    print("\n本工具将从以下数据源收集 ABC Atlas 的元数据:")
    print("  1. Allen Brain Map API")
    print("  2. CELLxGENE (Allen Institute 数据)")
    print("  3. ABC Atlas 专用数据库")
    print("  4. Allen S3 存储清单")
    print("\n" + "="*80 + "\n")
    
    # 创建收集器
    collector = ABCAtlasCollector(output_dir="./abc_atlas_metadata")
    
    # 执行收集
    df = collector.collect_all_metadata()
    
    # 可选: 添加引用计数
    add_citations = input("\n是否添加引用计数? (这将花费较长时间) [y/N]: ").strip().lower()
    
    if add_citations == 'y':
        enricher = CitationEnricher()
        df = enricher.enrich_dataframe(df)
        
        # 保存更新后的数据
        final_path = collector.save_processed_data(df, 'FINAL_abc_atlas_with_citations')
        logger.info(f"✓ 包含引用计数的最终文件: {final_path}")
    
    # 生成摘要报告
    print("\n" + "="*80)
    print("收集摘要")
    print("="*80)
    print(f"\n总数据集数量: {len(df)}")
    print(f"输出目录: {collector.output_dir}")
    print(f"\n主要输出文件:")
    print(f"  - 原始数据: raw_*.json")
    print(f"  - 处理数据: processed_*.csv")
    print(f"  - 最终整合: processed_FINAL_abc_atlas_metadata.csv")
    
    if not df.empty:
        print(f"\n数据覆盖范围:")
        if 'tissue' in df.columns:
            print(f"  - 组织类型: {df['tissue'].nunique()} 种")
        if 'experiment_design' in df.columns:
            print(f"  - 实验设计: {df['experiment_design'].nunique()} 种")
        if 'cell_count' in df.columns:
            total_cells = df['cell_count'].apply(lambda x: 
                int(re.sub(r'[^0-9]', '', str(x))) if pd.notna(x) and str(x).strip() else 0
            ).sum()
            print(f"  - 总细胞数: ~{total_cells:,}")
    
    print("\n" + "="*80)
    print("收集完成!")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()