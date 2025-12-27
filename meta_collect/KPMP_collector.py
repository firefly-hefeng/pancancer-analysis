"""
KPMP Single-Cell RNA-seq Metadata Collection Script
完整采集KPMP数据库中所有scRNA-seq研究的元数据
"""

import requests
import pandas as pd
import json
from datetime import datetime
import time
from bs4 import BeautifulSoup
import re
from typing import Dict, List, Optional
import logging
from pathlib import Path

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class KPMPDataCollector:
    """KPMP数据采集器"""
    
    def __init__(self, output_dir: str = "./kpmp_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # KPMP API endpoints
        self.base_urls = {
            'atlas': 'https://atlas.kpmp.org/api',
            'repository': 'https://repository.kpmp.org/api',
            'dataexplorer': 'https://dataexplorer.kpmp.org/api'
        }
        
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Research Data Collection)'
        })
        
        # 存储所有原始metadata
        self.raw_metadata = []
        
    def collect_atlas_data(self) -> List[Dict]:
        """从KPMP Atlas采集数据"""
        logger.info("开始采集KPMP Atlas数据...")
        atlas_data = []
        
        try:
            # Atlas API endpoints
            endpoints = [
                '/v1/datasets',
                '/v1/participants',
                '/v1/summary',
                '/v1/spatial'
            ]
            
            for endpoint in endpoints:
                url = f"{self.base_urls['atlas']}{endpoint}"
                try:
                    response = self.session.get(url, timeout=30)
                    if response.status_code == 200:
                        data = response.json()
                        atlas_data.append({
                            'source': 'atlas',
                            'endpoint': endpoint,
                            'data': data,
                            'timestamp': datetime.now().isoformat()
                        })
                        logger.info(f"成功采集: {endpoint}")
                    else:
                        logger.warning(f"无法访问 {endpoint}: {response.status_code}")
                except Exception as e:
                    logger.error(f"采集 {endpoint} 时出错: {str(e)}")
                
                time.sleep(1)  # 避免请求过快
                
        except Exception as e:
            logger.error(f"Atlas数据采集失败: {str(e)}")
        
        return atlas_data
    
    def collect_repository_data(self) -> List[Dict]:
        """从KPMP Repository采集数据"""
        logger.info("开始采集KPMP Repository数据...")
        repo_data = []
        
        try:
            # Repository包含文件级别的metadata
            endpoints = [
                '/v1/file',
                '/v1/package',
                '/v1/participant'
            ]
            
            for endpoint in endpoints:
                url = f"{self.base_urls['repository']}{endpoint}"
                try:
                    response = self.session.get(url, timeout=30)
                    if response.status_code == 200:
                        data = response.json()
                        repo_data.append({
                            'source': 'repository',
                            'endpoint': endpoint,
                            'data': data,
                            'timestamp': datetime.now().isoformat()
                        })
                        logger.info(f"成功采集: {endpoint}")
                except Exception as e:
                    logger.error(f"采集 {endpoint} 时出错: {str(e)}")
                
                time.sleep(1)
                
        except Exception as e:
            logger.error(f"Repository数据采集失败: {str(e)}")
        
        return repo_data
    
    def collect_dataexplorer_metadata(self) -> List[Dict]:
        """从Data Explorer采集metadata"""
        logger.info("开始采集KPMP Data Explorer数据...")
        explorer_data = []
        
        try:
            # Data Explorer提供更详细的实验metadata
            url = f"{self.base_urls['dataexplorer']}/metadata"
            
            # 尝试获取不同类型的数据
            data_types = ['scrna', 'snrna', 'spatial', 'bulk']
            
            for dtype in data_types:
                try:
                    params = {'dataType': dtype}
                    response = self.session.get(url, params=params, timeout=30)
                    
                    if response.status_code == 200:
                        data = response.json()
                        explorer_data.append({
                            'source': 'dataexplorer',
                            'data_type': dtype,
                            'data': data,
                            'timestamp': datetime.now().isoformat()
                        })
                        logger.info(f"成功采集 {dtype} 数据")
                except Exception as e:
                    logger.error(f"采集 {dtype} 时出错: {str(e)}")
                
                time.sleep(1)
                
        except Exception as e:
            logger.error(f"Data Explorer数据采集失败: {str(e)}")
        
        return explorer_data
    
    def scrape_web_metadata(self) -> List[Dict]:
        """从KPMP网站爬取额外的metadata"""
        logger.info("开始从KPMP网站爬取数据...")
        web_data = []
        
        urls = [
            'https://atlas.kpmp.org/explorer',
            'https://atlas.kpmp.org/repository',
            'https://www.kpmp.org/data-release'
        ]
        
        for url in urls:
            try:
                response = self.session.get(url, timeout=30)
                if response.status_code == 200:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # 提取页面中的JSON-LD或其他结构化数据
                    scripts = soup.find_all('script', type='application/ld+json')
                    for script in scripts:
                        try:
                            data = json.loads(script.string)
                            web_data.append({
                                'source': 'web_scraping',
                                'url': url,
                                'data': data,
                                'timestamp': datetime.now().isoformat()
                            })
                        except:
                            pass
                    
                    # 提取表格数据
                    tables = soup.find_all('table')
                    for i, table in enumerate(tables):
                        df = pd.read_html(str(table))[0]
                        web_data.append({
                            'source': 'web_table',
                            'url': url,
                            'table_index': i,
                            'data': df.to_dict('records'),
                            'timestamp': datetime.now().isoformat()
                        })
                
                logger.info(f"成功爬取: {url}")
                time.sleep(2)
                
            except Exception as e:
                logger.error(f"爬取 {url} 时出错: {str(e)}")
        
        return web_data
    
    def query_geo_crossref(self) -> List[Dict]:
        """查询GEO/SRA中的KPMP相关数据"""
        logger.info("查询GEO/SRA数据库...")
        geo_data = []
        
        # 使用NCBI E-utilities查询
        search_terms = [
            'KPMP[All Fields]',
            'Kidney Precision Medicine Project[All Fields]',
            'KPMP AND single cell[All Fields]'
        ]
        
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
        
        for term in search_terms:
            try:
                # 搜索GEO
                search_url = f"{base_url}/esearch.fcgi"
                params = {
                    'db': 'gds',
                    'term': term,
                    'retmode': 'json',
                    'retmax': 100
                }
                
                response = self.session.get(search_url, params=params, timeout=30)
                if response.status_code == 200:
                    search_results = response.json()
                    
                    if 'esearchresult' in search_results and 'idlist' in search_results['esearchresult']:
                        ids = search_results['esearchresult']['idlist']
                        
                        # 获取详细信息
                        if ids:
                            fetch_url = f"{base_url}/esummary.fcgi"
                            fetch_params = {
                                'db': 'gds',
                                'id': ','.join(ids),
                                'retmode': 'json'
                            }
                            
                            fetch_response = self.session.get(fetch_url, params=fetch_params, timeout=30)
                            if fetch_response.status_code == 200:
                                geo_data.append({
                                    'source': 'geo',
                                    'search_term': term,
                                    'data': fetch_response.json(),
                                    'timestamp': datetime.now().isoformat()
                                })
                
                logger.info(f"GEO查询完成: {term}")
                time.sleep(1)
                
            except Exception as e:
                logger.error(f"GEO查询出错 ({term}): {str(e)}")
        
        return geo_data
    
    def query_publications(self) -> List[Dict]:
        """查询PubMed中的KPMP相关文献"""
        logger.info("查询PubMed文献...")
        pubmed_data = []
        
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
        
        search_terms = [
            'KPMP[Affiliation] AND single cell sequencing',
            'Kidney Precision Medicine Project AND transcriptomics',
            'KPMP AND (scRNA-seq OR single-cell RNA-seq)'
        ]
        
        for term in search_terms:
            try:
                # 搜索PubMed
                search_url = f"{base_url}/esearch.fcgi"
                params = {
                    'db': 'pubmed',
                    'term': term,
                    'retmode': 'json',
                    'retmax': 200,
                    'sort': 'relevance'
                }
                
                response = self.session.get(search_url, params=params, timeout=30)
                if response.status_code == 200:
                    search_results = response.json()
                    
                    if 'esearchresult' in search_results and 'idlist' in search_results['esearchresult']:
                        pmids = search_results['esearchresult']['idlist']
                        
                        if pmids:
                            # 批量获取文献详细信息
                            fetch_url = f"{base_url}/efetch.fcgi"
                            fetch_params = {
                                'db': 'pubmed',
                                'id': ','.join(pmids),
                                'retmode': 'xml'
                            }
                            
                            fetch_response = self.session.get(fetch_url, params=fetch_params, timeout=30)
                            if fetch_response.status_code == 200:
                                pubmed_data.append({
                                    'source': 'pubmed',
                                    'search_term': term,
                                    'pmids': pmids,
                                    'data': fetch_response.text,
                                    'timestamp': datetime.now().isoformat()
                                })
                
                logger.info(f"PubMed查询完成: {term} (找到 {len(pmids) if pmids else 0} 篇)")
                time.sleep(1)
                
            except Exception as e:
                logger.error(f"PubMed查询出错 ({term}): {str(e)}")
        
        return pubmed_data
    
    def collect_all_raw_metadata(self):
        """采集所有原始metadata"""
        logger.info("=" * 50)
        logger.info("开始全面采集KPMP数据库metadata")
        logger.info("=" * 50)
        
        # 1. KPMP Atlas数据
        atlas_data = self.collect_atlas_data()
        self.save_raw_data(atlas_data, 'atlas_raw_metadata.json')
        
        # 2. KPMP Repository数据
        repo_data = self.collect_repository_data()
        self.save_raw_data(repo_data, 'repository_raw_metadata.json')
        
        # 3. Data Explorer数据
        explorer_data = self.collect_dataexplorer_metadata()
        self.save_raw_data(explorer_data, 'dataexplorer_raw_metadata.json')
        
        # 4. 网站爬取数据
        web_data = self.scrape_web_metadata()
        self.save_raw_data(web_data, 'web_raw_metadata.json')
        
        # 5. GEO交叉引用
        geo_data = self.query_geo_crossref()
        self.save_raw_data(geo_data, 'geo_crossref_metadata.json')
        
        # 6. PubMed文献
        pubmed_data = self.query_publications()
        self.save_raw_data(pubmed_data, 'pubmed_metadata.json')
        
        # 合并所有原始数据
        self.raw_metadata = {
            'atlas': atlas_data,
            'repository': repo_data,
            'dataexplorer': explorer_data,
            'web': web_data,
            'geo': geo_data,
            'pubmed': pubmed_data,
            'collection_timestamp': datetime.now().isoformat()
        }
        
        self.save_raw_data(self.raw_metadata, 'all_raw_metadata.json')
        logger.info("所有原始metadata采集完成!")
        
        return self.raw_metadata
    
    def save_raw_data(self, data: any, filename: str):
        """保存原始数据到JSON文件"""
        filepath = self.output_dir / filename
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        logger.info(f"数据已保存到: {filepath}")
    
    def parse_to_series_table(self, raw_metadata: Dict) -> pd.DataFrame:
        """将原始metadata解析为标准series表格"""
        logger.info("开始解析metadata到series表格...")
        
        series_records = []
        
        # 解析Atlas数据
        if 'atlas' in raw_metadata:
            for item in raw_metadata['atlas']:
                if 'data' in item and isinstance(item['data'], dict):
                    records = self._parse_atlas_data(item['data'])
                    series_records.extend(records)
        
        # 解析Repository数据
        if 'repository' in raw_metadata:
            for item in raw_metadata['repository']:
                if 'data' in item:
                    records = self._parse_repository_data(item['data'])
                    series_records.extend(records)
        
        # 解析Data Explorer数据
        if 'dataexplorer' in raw_metadata:
            for item in raw_metadata['dataexplorer']:
                if 'data' in item:
                    records = self._parse_explorer_data(item['data'])
                    series_records.extend(records)
        
        # 解析GEO数据
        if 'geo' in raw_metadata:
            for item in raw_metadata['geo']:
                if 'data' in item:
                    records = self._parse_geo_data(item['data'])
                    series_records.extend(records)
        
        # 创建DataFrame
        df = pd.DataFrame(series_records)
        
        # 确保包含所有必需字段
        required_columns = [
            'id', 'title', 'disease_general', 'disease', 'pubmed',
            'source_database', 'access_link', 'open_status', 'ethnicity',
            'sex', 'tissue', 'sequencing_platform', 'experiment_design',
            'sample_type', 'summary', 'citation_count', 'publication_date',
            'submission_date', 'last_update_date', 'contact_name',
            'contact_email', 'contact_institute', 'data_tier',
            'tissue_location', 'supplementary_information'
        ]
        
        for col in required_columns:
            if col not in df.columns:
                df[col] = None
        
        # 去重并排序
        df = df.drop_duplicates(subset=['id'], keep='first')
        df = df.sort_values('submission_date', ascending=False)
        
        return df[required_columns]
    
    def _parse_atlas_data(self, data: Dict) -> List[Dict]:
        """解析Atlas数据"""
        records = []
        
        # Atlas通常包含datasets信息
        if isinstance(data, list):
            for dataset in data:
                record = self._create_base_record()
                
                record['id'] = dataset.get('id', f"KPMP_ATLAS_{len(records)}")
                record['title'] = dataset.get('title') or dataset.get('name')
                record['source_database'] = 'KPMP_Atlas'
                record['access_link'] = f"https://atlas.kpmp.org/explorer/{dataset.get('id', '')}"
                
                # 提取组织信息
                record['tissue'] = self._extract_tissue(dataset)
                record['tissue_location'] = dataset.get('tissue_location') or dataset.get('anatomical_site')
                
                # 提取疾病信息
                disease_info = self._extract_disease(dataset)
                record['disease'] = disease_info.get('specific')
                record['disease_general'] = disease_info.get('general')
                
                # 提取实验信息
                record['sequencing_platform'] = dataset.get('platform') or dataset.get('technology')
                record['experiment_design'] = dataset.get('experiment_type')
                record['sample_type'] = dataset.get('sample_type') or 'single_cell'
                
                # 提取人群信息
                record['sex'] = dataset.get('sex') or dataset.get('gender')
                record['ethnicity'] = dataset.get('race') or dataset.get('ethnicity')
                
                # 数据层级
                record['data_tier'] = self._determine_data_tier(dataset)
                
                # 状态
                record['open_status'] = dataset.get('access_type', 'controlled')
                
                # 日期
                record['submission_date'] = dataset.get('created_date') or dataset.get('upload_date')
                record['last_update_date'] = dataset.get('modified_date') or dataset.get('update_date')
                
                # 摘要
                record['summary'] = dataset.get('description') or dataset.get('summary')
                
                # 补充信息
                record['supplementary_information'] = json.dumps({
                    'n_cells': dataset.get('cell_count'),
                    'n_samples': dataset.get('sample_count'),
                    'protocols': dataset.get('protocols'),
                    'original_data': dataset
                }, ensure_ascii=False)
                
                records.append(record)
        
        return records
    
    def _parse_repository_data(self, data: Dict) -> List[Dict]:
        """解析Repository数据"""
        records = []
        
        if isinstance(data, dict) and 'packages' in data:
            for package in data['packages']:
                record = self._create_base_record()
                
                record['id'] = package.get('packageId', f"KPMP_REPO_{len(records)}")
                record['title'] = package.get('packageName')
                record['source_database'] = 'KPMP_Repository'
                record['access_link'] = f"https://repository.kpmp.org/package/{package.get('packageId', '')}"
                
                # 从files中提取技术信息
                files = package.get('files', [])
                for file in files:
                    file_type = file.get('fileType', '').lower()
                    if 'scrna' in file_type or 'single' in file_type:
                        record['experiment_design'] = 'scRNA-seq'
                        record['sample_type'] = 'single_cell'
                        break
                
                # 数据层级
                record['data_tier'] = 'raw' if any('raw' in f.get('fileType', '').lower() for f in files) else 'processed'
                
                # 提取participant信息
                participants = package.get('participants', [])
                if participants:
                    record['sex'] = participants[0].get('sex')
                    record['ethnicity'] = participants[0].get('race')
                
                # 组织和疾病
                record['tissue'] = package.get('tissueType') or 'kidney'
                record['disease'] = package.get('protocol')
                
                # 日期
                record['submission_date'] = package.get('createdAt')
                record['last_update_date'] = package.get('updatedAt')
                
                record['supplementary_information'] = json.dumps(package, ensure_ascii=False)
                
                records.append(record)
        
        return records
    
    def _parse_explorer_data(self, data: Dict) -> List[Dict]:
        """解析Data Explorer数据"""
        records = []
        
        if isinstance(data, dict):
            for key, dataset in data.items():
                if isinstance(dataset, dict):
                    record = self._create_base_record()
                    
                    record['id'] = dataset.get('datasetId', f"KPMP_EXPLORER_{key}")
                    record['title'] = dataset.get('title') or dataset.get('datasetName')
                    record['source_database'] = 'KPMP_DataExplorer'
                    
                    # 技术平台
                    record['sequencing_platform'] = dataset.get('platform') or dataset.get('instrument')
                    record['experiment_design'] = dataset.get('assayType')
                    
                    # 样本信息
                    record['tissue'] = dataset.get('tissue', 'kidney')
                    record['sex'] = dataset.get('sex')
                    
                    record['supplementary_information'] = json.dumps(dataset, ensure_ascii=False)
                    
                    records.append(record)
        
        return records
    
    def _parse_geo_data(self, data: Dict) -> List[Dict]:
        """解析GEO数据"""
        records = []
        
        if isinstance(data, dict) and 'result' in data:
            for gse_id, gse_data in data['result'].items():
                if gse_id == 'uids':
                    continue
                
                record = self._create_base_record()
                
                record['id'] = gse_data.get('accession', gse_id)
                record['title'] = gse_data.get('title')
                record['source_database'] = 'GEO_KPMP'
                record['access_link'] = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={record['id']}"
                
                # PubMed链接
                pubmed_ids = gse_data.get('pubmedids', [])
                if pubmed_ids:
                    record['pubmed'] = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_ids[0]}/"
                
                # 平台信息
                record['sequencing_platform'] = gse_data.get('gpl')
                
                # 日期
                record['submission_date'] = gse_data.get('pdat')
                record['publication_date'] = gse_data.get('pubdate')
                
                # 摘要
                record['summary'] = gse_data.get('summary')
                
                # 样本数
                n_samples = gse_data.get('n_samples', 0)
                
                record['supplementary_information'] = json.dumps({
                    'n_samples': n_samples,
                    'organism': gse_data.get('taxon'),
                    'gse_type': gse_data.get('gdstype'),
                    'original_data': gse_data
                }, ensure_ascii=False)
                
                records.append(record)
        
        return records
    
    def _create_base_record(self) -> Dict:
        """创建基础记录模板"""
        return {
            'id': None,
            'title': None,
            'disease_general': 'kidney_disease',
            'disease': None,
            'pubmed': None,
            'source_database': 'KPMP',
            'access_link': None,
            'open_status': 'controlled',
            'ethnicity': None,
            'sex': None,
            'tissue': 'kidney',
            'sequencing_platform': None,
            'experiment_design': None,
            'sample_type': 'single_cell',
            'summary': None,
            'citation_count': None,
            'publication_date': None,
            'submission_date': None,
            'last_update_date': None,
            'contact_name': None,
            'contact_email': None,
            'contact_institute': 'KPMP Consortium',
            'data_tier': None,
            'tissue_location': None,
            'supplementary_information': None
        }
    
    def _extract_tissue(self, data: Dict) -> str:
        """提取组织信息"""
        tissue_fields = ['tissue', 'organ', 'tissue_type', 'sample_type']
        for field in tissue_fields:
            if field in data and data[field]:
                return data[field]
        return 'kidney'
    
    def _extract_disease(self, data: Dict) -> Dict:
        """提取疾病信息"""
        disease_info = {
            'general': 'kidney_disease',
            'specific': None
        }
        
        disease_fields = ['disease', 'diagnosis', 'condition', 'protocol']
        for field in disease_fields:
            if field in data and data[field]:
                disease_info['specific'] = data[field]
                break
        
        # 根据具体疾病判断一般分类
        if disease_info['specific']:
            disease_lower = disease_info['specific'].lower()
            if 'ckd' in disease_lower or 'chronic' in disease_lower:
                disease_info['general'] = 'chronic_kidney_disease'
            elif 'aki' in disease_lower or 'acute' in disease_lower:
                disease_info['general'] = 'acute_kidney_injury'
            elif 'diabetic' in disease_lower:
                disease_info['general'] = 'diabetic_kidney_disease'
            elif 'glomerulo' in disease_lower:
                disease_info['general'] = 'glomerulonephritis'
        
        return disease_info
    
    def _determine_data_tier(self, data: Dict) -> str:
        """确定数据层级"""
        if 'data_tier' in data:
            return data['data_tier']
        
        # 根据可用文件判断
        files = data.get('files', [])
        has_raw = any('raw' in str(f).lower() for f in files)
        has_processed = any('processed' in str(f).lower() or 'matrix' in str(f).lower() for f in files)
        
        if has_raw and has_processed:
            return 'raw/processed_matrix'
        elif has_raw:
            return 'raw'
        elif has_processed:
            return 'processed_matrix'
        else:
            return 'unknown'
    
    def run_full_collection(self):
        """运行完整的数据采集流程"""
        logger.info("启动KPMP数据库完整采集流程...")
        
        # 第一步: 采集所有原始metadata
        raw_metadata = self.collect_all_raw_metadata()
        
        # 第二步: 解析为标准series表格
        series_df = self.parse_to_series_table(raw_metadata)
        
        # 保存series表格
        series_file = self.output_dir / 'kpmp_series_metadata.csv'
        series_df.to_csv(series_file, index=False, encoding='utf-8-sig')
        logger.info(f"Series表格已保存到: {series_file}")
        
        # 保存Excel格式(带格式)
        excel_file = self.output_dir / 'kpmp_series_metadata.xlsx'
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            series_df.to_excel(writer, sheet_name='KPMP_Series', index=False)
        logger.info(f"Excel文件已保存到: {excel_file}")
        
        # 生成统计报告
        self.generate_summary_report(series_df)
        
        logger.info("=" * 50)
        logger.info("KPMP数据采集完成!")
        logger.info(f"共采集 {len(series_df)} 条记录")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info("=" * 50)
        
        return series_df
    
    def generate_summary_report(self, df: pd.DataFrame):
        """生成采集总结报告"""
        report_file = self.output_dir / 'collection_summary_report.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 70 + "\n")
            f.write("KPMP单细胞测序数据采集总结报告\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"采集时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总记录数: {len(df)}\n\n")
            
            # 数据来源统计
            f.write("数据来源分布:\n")
            if 'source_database' in df.columns:
                source_counts = df['source_database'].value_counts()
                for source, count in source_counts.items():
                    f.write(f"  {source}: {count}\n")
            f.write("\n")
            
            # 实验设计统计
            f.write("实验设计类型:\n")
            if 'experiment_design' in df.columns:
                exp_counts = df['experiment_design'].value_counts()
                for exp, count in exp_counts.items():
                    f.write(f"  {exp}: {count}\n")
            f.write("\n")
            
            # 测序平台统计
            f.write("测序平台分布:\n")
            if 'sequencing_platform' in df.columns:
                platform_counts = df['sequencing_platform'].value_counts()
                for platform, count in platform_counts.items():
                    f.write(f"  {platform}: {count}\n")
            f.write("\n")
            
            # 疾病类型统计
            f.write("疾病类型分布:\n")
            if 'disease_general' in df.columns:
                disease_counts = df['disease_general'].value_counts()
                for disease, count in disease_counts.items():
                    f.write(f"  {disease}: {count}\n")
            f.write("\n")
            
            # 数据完整性统计
            f.write("数据完整性:\n")
            for col in df.columns:
                non_null = df[col].notna().sum()
                completeness = (non_null / len(df)) * 100
                f.write(f"  {col}: {completeness:.1f}% ({non_null}/{len(df)})\n")
            
            f.write("\n" + "=" * 70 + "\n")
        
        logger.info(f"总结报告已保存到: {report_file}")


def main():
    """主函数"""
    # 创建采集器实例
    collector = KPMPDataCollector(output_dir="./kpmp_metadata_collection")
    
    # 运行完整采集
    series_df = collector.run_full_collection()
    
    # 显示前几行结果
    print("\n采集结果预览:")
    print(series_df.head())
    print(f"\n数据维度: {series_df.shape}")
    
    return series_df


if __name__ == "__main__":
    df = main()