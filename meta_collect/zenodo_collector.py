import requests
import pandas as pd
import time
import json
from datetime import datetime
import re
from typing import Dict, List, Optional
import logging
from pathlib import Path
import pickle
import urllib3
import os

# 禁用SSL警告（如果需要）
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('zenodo_scraping.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class ZenodoScRNASeqCollector:
    """
    Zenodo单细胞测序数据收集器
    使用Zenodo REST API v1进行数据检索和元数据收集
    """
    
    def __init__(self, access_token: Optional[str] = None, use_proxy: bool = True, 
                 proxy_url: Optional[str] = None, verify_ssl: bool = True):
        """
        初始化收集器
        
        Parameters:
        -----------
        access_token : str, optional
            Zenodo API访问令牌
        use_proxy : bool
            是否使用代理
        proxy_url : str, optional
            代理URL，格式如 'http://proxy.example.com:7890'
        verify_ssl : bool
            是否验证SSL证书
        """
        self.base_url = "https://zenodo.org/api"
        self.access_token = access_token
        self.verify_ssl = verify_ssl
        
        self.headers = {
            'Accept': 'application/json',
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        }
        
        if access_token:
            self.headers['Authorization'] = f'Bearer {access_token}'
        
        # 设置代理 - 使用您配置的代理服务器
        self.proxies = None
        if use_proxy:
            if proxy_url:
                self.proxies = {
                    'http': proxy_url,
                    'https': proxy_url
                }
            else:
                # 默认使用您配置的代理
                default_proxy = 'http://proxy.mornai.cn:7890'
                self.proxies = {
                    'http': default_proxy,
                    'https': default_proxy
                }
            logger.info(f"使用代理: {self.proxies['http']}")
        
        # 创建session以复用连接
        self.session = requests.Session()
        self.session.headers.update(self.headers)
        if self.proxies:
            self.session.proxies.update(self.proxies)
        
        # 定义搜索关键词
        self.search_keywords = [
            'single cell RNA-seq',
            'scRNA-seq',
            'single-cell transcriptomics',
            'single cell sequencing',
            '10x Genomics',
            'Drop-seq',
            'Smart-seq',
            'CITE-seq',
            'single nucleus RNA-seq',
            'snRNA-seq',
            'single cell gene expression',
            'sc-RNA-seq',
            'single-cell RNA sequencing'
        ]
        
        self.all_records = []
        self.raw_metadata_dir = Path('zenodo_raw_metadata')
        self.raw_metadata_dir.mkdir(exist_ok=True)
    
    def test_connection(self) -> bool:
        """
        测试与Zenodo的连接
        
        Returns:
        --------
        bool : 连接是否成功
        """
        logger.info("=" * 80)
        logger.info("测试与Zenodo的连接...")
        logger.info("=" * 80)
        
        try:
            response = self.session.get(
                f"{self.base_url}/records",
                params={'size': 1},
                timeout=30,
                verify=self.verify_ssl
            )
            
            if response.status_code == 200:
                logger.info("✓ 连接成功！")
                logger.info(f"✓ 状态码: {response.status_code}")
                data = response.json()
                total = data.get('hits', {}).get('total', 0)
                logger.info(f"✓ Zenodo总记录数: {total:,}")
                logger.info("=" * 80)
                return True
            else:
                logger.error(f"✗ 连接失败，状态码: {response.status_code}")
                logger.error(f"✗ 响应内容: {response.text[:200]}")
                logger.info("=" * 80)
                return False
                
        except requests.exceptions.ProxyError as e:
            logger.error(f"✗ 代理错误: {str(e)}")
            logger.info("提示: 请检查代理设置")
            logger.info(f"当前代理配置: {self.proxies}")
            logger.info("=" * 80)
            return False
            
        except requests.exceptions.SSLError as e:
            logger.error(f"✗ SSL错误: {str(e)}")
            logger.info("提示: SSL证书验证失败，可以尝试设置 verify_ssl=False")
            logger.info("=" * 80)
            return False
            
        except requests.exceptions.ConnectionError as e:
            logger.error(f"✗ 连接错误: {str(e)}")
            logger.info("提示: 无法建立连接，请检查:")
            logger.info("  1. 代理服务器是否可访问")
            logger.info("  2. 网络连接是否正常")
            logger.info("  3. 防火墙设置")
            logger.info("=" * 80)
            return False
            
        except requests.exceptions.Timeout as e:
            logger.error(f"✗ 连接超时: {str(e)}")
            logger.info("提示: 请求超时，可能是网络较慢或代理响应慢")
            logger.info("=" * 80)
            return False
            
        except Exception as e:
            logger.error(f"✗ 未知错误: {str(e)}")
            logger.error(f"错误类型: {type(e).__name__}")
            logger.info("=" * 80)
            return False
    
    def search_records(self, max_records: int = 10000) -> List[Dict]:
        """
        搜索Zenodo中的单细胞测序相关记录
        
        Parameters:
        -----------
        max_records : int
            最大检索记录数
            
        Returns:
        --------
        List[Dict] : 检索到的所有记录
        """
        all_results = []
        seen_ids = set()
        
        for keyword in self.search_keywords:
            logger.info(f"\n{'='*60}")
            logger.info(f"正在搜索关键词: {keyword}")
            logger.info('='*60)
            
            # 构建查询 - 搜索人类单细胞数据
            query = f'("{keyword}") AND (human OR "homo sapiens")'
            page = 1
            size = 100  # 每页100条记录
            keyword_count = 0
            
            while len(all_results) < max_records:
                params = {
                    'q': query,
                    'size': size,
                    'page': page,
                    'sort': 'mostrecent',
                    'type': 'dataset'
                }
                
                try:
                    logger.info(f"  请求第 {page} 页...")
                    response = self.session.get(
                        f"{self.base_url}/records",
                        params=params,
                        timeout=60,
                        verify=self.verify_ssl
                    )
                    
                    if response.status_code == 200:
                        data = response.json()
                        hits = data.get('hits', {}).get('hits', [])
                        total = data.get('hits', {}).get('total', 0)
                        
                        if not hits:
                            logger.info(f"  ✓ 关键词 '{keyword}' 搜索完成，共找到 {keyword_count} 条记录")
                            break
                        
                        # 去重
                        new_records = 0
                        for hit in hits:
                            hit_id = hit.get('id')
                            if hit_id not in seen_ids:
                                seen_ids.add(hit_id)
                                all_results.append(hit)
                                new_records += 1
                                keyword_count += 1
                        
                        logger.info(f"  ✓ 第 {page} 页: 获取 {len(hits)} 条，新增 {new_records} 条")
                        logger.info(f"  ✓ 当前总计: {len(all_results)} 条唯一记录")
                        
                        # 保存原始JSON
                        self._save_raw_json(hits, f"{keyword.replace(' ', '_')}_page_{page}")
                        
                        # 如果已经获取了所有结果，跳出
                        if page * size >= total:
                            logger.info(f"  ✓ 已获取该关键词的所有结果")
                            break
                        
                        page += 1
                        time.sleep(1)  # 避免请求过快
                        
                    elif response.status_code == 429:
                        logger.warning("  ⚠ 请求过快，等待60秒...")
                        time.sleep(60)
                        continue
                        
                    else:
                        logger.error(f"  ✗ 请求失败: HTTP {response.status_code}")
                        logger.debug(f"  响应内容: {response.text[:500]}")
                        break
                        
                except requests.exceptions.Timeout:
                    logger.warning(f"  ⚠ 请求超时，跳过剩余页面")
                    break
                    
                except Exception as e:
                    logger.error(f"  ✗ 搜索出错: {str(e)}")
                    break
                
                # 检查是否达到最大记录数
                if len(all_results) >= max_records:
                    logger.info(f"\n✓ 已达到最大记录数限制 ({max_records})，停止搜索")
                    break
        
        self.all_records = all_results
        
        logger.info("\n" + "=" * 80)
        logger.info(f"搜索完成！总共收集到 {len(self.all_records)} 条唯一记录")
        logger.info("=" * 80)
        return self.all_records
    
    def _save_raw_json(self, data: List[Dict], filename: str):
        """保存原始JSON数据"""
        filepath = self.raw_metadata_dir / f"{filename}.json"
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            logger.warning(f"保存JSON文件失败: {str(e)}")
    
    def get_record_details(self, record_id: str) -> Optional[Dict]:
        """获取单个记录的详细信息"""
        try:
            response = self.session.get(
                f"{self.base_url}/records/{record_id}",
                timeout=30,
                verify=self.verify_ssl
            )
            
            if response.status_code == 200:
                return response.json()
            else:
                logger.warning(f"无法获取记录 {record_id}: HTTP {response.status_code}")
                return None
                
        except Exception as e:
            logger.error(f"获取记录详情出错 {record_id}: {str(e)}")
            return None
    
    def enrich_metadata(self):
        """为所有收集的记录获取详细信息"""
        logger.info("\n" + "=" * 80)
        logger.info("开始丰富元数据（获取详细记录信息）...")
        logger.info("=" * 80)
        
        enriched_records = []
        total = len(self.all_records)
        
        for i, record in enumerate(self.all_records):
            record_id = record['id']
            
            if (i + 1) % 10 == 0 or i == 0:
                logger.info(f"处理进度: {i+1}/{total} ({(i+1)/total*100:.1f}%)")
            
            # 如果记录已经包含完整信息，直接使用
            if 'files' in record and 'metadata' in record and record.get('files'):
                enriched_records.append(record)
                continue
            
            # 否则获取详细信息
            detailed_record = self.get_record_details(record_id)
            if detailed_record:
                enriched_records.append(detailed_record)
                self._save_raw_json([detailed_record], f"detail_{record_id}")
            else:
                enriched_records.append(record)
            
            time.sleep(0.5)  # 避免请求过快
        
        self.all_records = enriched_records
        logger.info(f"✓ 元数据丰富完成，共 {len(enriched_records)} 条记录")
        logger.info("=" * 80)
    
    def extract_pubmed_id(self, text: str) -> Optional[str]:
        """从文本中提取PubMed ID"""
        if not text:
            return None
        
        patterns = [
            r'PMID:\s*(\d+)',
            r'PubMed:\s*(\d+)',
            r'pubmed/(\d+)',
            r'ncbi\.nlm\.nih\.gov/pubmed/(\d+)',
            r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)'
        ]
        
        for pattern in patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                return match.group(1)
        
        return None
    
    def extract_doi(self, record: Dict) -> Optional[str]:
        """提取DOI"""
        metadata = record.get('metadata', {})
        
        # 首先查看metadata中的doi字段
        doi = metadata.get('doi')
        if doi:
            return doi
        
        # 其次查看related_identifiers
        for identifier in metadata.get('related_identifiers', []):
            if identifier.get('scheme') == 'doi':
                return identifier.get('identifier')
        
        return None
    
    def classify_disease(self, text: str) -> tuple:
        """分类疾病类型"""
        if not text:
            return (None, None)
        
        text_lower = text.lower()
        
        disease_categories = {
            'cancer': ['cancer', 'tumor', 'carcinoma', 'leukemia', 'lymphoma', 'melanoma', 
                      'glioma', 'sarcoma', 'adenocarcinoma', 'glioblastoma', 'neuroblastoma',
                      'hepatocellular', 'metastasis', 'malignant', 'oncology'],
            'neurological': ['alzheimer', 'parkinson', 'dementia', 'brain', 'neural', 'neuron',
                           'neurological', 'neurodegenerative', 'epilepsy', 'stroke', 'autism',
                           'multiple sclerosis', 'huntington'],
            'immunological': ['autoimmune', 'inflammation', 'immune', 'arthritis', 'lupus',
                            'immunology', 'immunodeficiency', 'allergy', 'psoriasis'],
            'cardiovascular': ['heart', 'cardiac', 'cardiovascular', 'vascular', 'atherosclerosis',
                             'hypertension', 'coronary', 'myocardial'],
            'metabolic': ['diabetes', 'obesity', 'metabolic', 'insulin', 'glucose'],
            'infectious': ['covid', 'virus', 'bacterial', 'infection', 'sepsis', 'viral',
                         'influenza', 'hepatitis', 'hiv', 'tuberculosis', 'malaria'],
            'respiratory': ['lung', 'pulmonary', 'respiratory', 'asthma', 'copd', 'pneumonia'],
            'renal': ['kidney', 'renal', 'nephropathy'],
            'healthy': ['normal', 'healthy', 'control', 'wild-type', 'wild type', 'non-diseased']
        }
        
        disease_general = None
        disease_specific = []
        
        for category, keywords in disease_categories.items():
            for keyword in keywords:
                if keyword in text_lower:
                    if not disease_general:
                        disease_general = category
                    disease_specific.append(keyword)
        
        disease_specific_str = '; '.join(set(disease_specific)) if disease_specific else None
        
        return (disease_general, disease_specific_str)
    
    def extract_ethnicity(self, text: str) -> Optional[str]:
        """提取种族/人群信息"""
        if not text:
            return None
        
        text_lower = text.lower()
        
        ethnicities = {
            'Asian': ['asian', 'chinese', 'japanese', 'korean', 'indian', 'south asian', 'east asian'],
            'European': ['european', 'caucasian', 'white'],
            'African': ['african', 'black', 'african american'],
            'Hispanic': ['hispanic', 'latino', 'latina'],
            'Middle Eastern': ['middle eastern', 'arab'],
            'Mixed': ['mixed', 'diverse', 'multi-ethnic']
        }
        
        found = []
        for ethnicity, keywords in ethnicities.items():
            for keyword in keywords:
                if keyword in text_lower:
                    found.append(ethnicity)
                    break
        
        return '; '.join(found) if found else None
    
    def extract_sex(self, text: str) -> Optional[str]:
        """提取性别信息"""
        if not text:
            return None
        
        text_lower = text.lower()
        
        # 检查男性关键词
        has_male = any(word in text_lower for word in ['male', 'men', 'man', 'males'])
        # 检查女性关键词
        has_female = any(word in text_lower for word in ['female', 'women', 'woman', 'females'])
        
        if has_male and has_female:
            return 'Mixed'
        elif has_male:
            return 'Male'
        elif has_female:
            return 'Female'
        
        return None
    
    def extract_tissue(self, text: str) -> Optional[str]:
        """提取组织类型"""
        if not text:
            return None
        
        text_lower = text.lower()
        
        tissues = [
            'blood', 'pbmc', 'peripheral blood', 'bone marrow',
            'brain', 'cortex', 'hippocampus', 'cerebellum',
            'liver', 'lung', 'heart', 'kidney',
            'skin', 'muscle', 'bone', 'pancreas', 'spleen', 'thymus',
            'intestine', 'colon', 'stomach', 'esophagus',
            'breast', 'ovary', 'testis', 'uterus',
            'prostate', 'bladder',
            'retina', 'cornea', 'eye',
            'tumor', 'organoid', 'spheroid',
            'lymph node', 'tonsil',
            'adipose', 'fat',
            'placenta', 'umbilical cord'
        ]
        
        found = [tissue for tissue in tissues if tissue in text_lower]
        
        return '; '.join(found) if found else None
    
    def extract_platform(self, text: str) -> Optional[str]:
        """提取测序平台"""
        if not text:
            return None
        
        text_lower = text.lower()
        
        platforms = {
            '10x Genomics': ['10x', '10x genomics', 'chromium', '10x chromium'],
            'Drop-seq': ['drop-seq', 'dropseq'],
            'Smart-seq': ['smart-seq', 'smartseq', 'smart-seq2', 'smart-seq3', 'smartseq2'],
            'CITE-seq': ['cite-seq', 'citeseq'],
            'Seq-Well': ['seq-well', 'seqwell'],
            'inDrop': ['indrop'],
            'BD Rhapsody': ['bd rhapsody', 'rhapsody'],
            'STRT-seq': ['strt-seq'],
            'Quartz-Seq': ['quartz-seq'],
            'CEL-seq': ['cel-seq', 'celseq'],
            'MARS-seq': ['mars-seq'],
            'sci-RNA-seq': ['sci-rna-seq']
        }
        
        found = []
        for platform, keywords in platforms.items():
            for keyword in keywords:
                if keyword in text_lower:
                    found.append(platform)
                    break
        
        return '; '.join(found) if found else None
    
    def extract_data_tier(self, files: List[Dict]) -> str:
        """判断数据层级"""
        if not files:
            return 'Unknown'
        
        tiers = set()
        
        for file in files:
            filename = file.get('key', '').lower()
            
            # 原始数据
            if any(ext in filename for ext in ['.fastq', '.fq', '.bam', '.sam', '.fastq.gz', '.fq.gz']):
                tiers.add('raw')
            
            # 处理后的数据
            if any(ext in filename for ext in ['.h5ad', '.rds', '.h5', '.loom', '.h5ad']):
                tiers.add('processed')
            
            # 矩阵数据
            if any(word in filename for word in ['matrix', 'count', 'expression']):
                if any(ext in filename for ext in ['.mtx', '.csv', '.tsv', '.txt']):
                    tiers.add('processed_matrix')
        
        if not tiers:
            return 'Unknown'
        
        return '; '.join(sorted(tiers))
    
    def parse_single_record(self, record: Dict) -> Dict:
        """解析单条记录为目标格式"""
        metadata = record.get('metadata', {})
        
        # 基本信息
        record_id = str(record.get('id', ''))
        title = metadata.get('title', '')
        description = metadata.get('description', '')
        
        # 合并所有文本用于信息提取
        combined_text = f"{title} {description}"
        for creator in metadata.get('creators', []):
            combined_text += f" {creator.get('name', '')}"
        combined_text += ' ' + ' '.join([kw for kw in metadata.get('keywords', [])])
        
        # 提取各类信息
        disease_general, disease_specific = self.classify_disease(combined_text)
        pubmed_id = self.extract_pubmed_id(combined_text)
        doi = self.extract_doi(record)
        
        # 文件信息
        files = record.get('files', [])
        data_tier = self.extract_data_tier(files)
        
        # 访问权限
        access_right = metadata.get('access_right', 'unknown')
        if access_right == 'open':
            open_status = 'Open'
        elif access_right == 'restricted':
            open_status = 'Restricted'
        elif access_right == 'embargoed':
            open_status = 'Embargoed'
        elif access_right == 'closed':
            open_status = 'Closed'
        else:
            open_status = 'Unknown'
        
        # 联系信息
        creators = metadata.get('creators', [])
        contact_name = creators[0].get('name') if creators else None
        contact_affiliation = creators[0].get('affiliation') if creators else None
        
        # 日期信息
        publication_date = metadata.get('publication_date')
        created = record.get('created')
        updated = record.get('updated')
        
        # 访问链接
        links = record.get('links', {})
        access_link = links.get('html', links.get('self', f"https://zenodo.org/record/{record_id}"))
        
        # 构建解析后的数据
        parsed_data = {
            'id': record_id,
            'sample_id': record_id,
            'title': title,
            'disease_general': disease_general,
            'disease': disease_specific,
            'pubmed': pubmed_id,
            'doi': doi,
            'source_database': 'Zenodo',
            'access_link': access_link,
            'open_status': open_status,
            'ethnicity': self.extract_ethnicity(combined_text),
            'sex': self.extract_sex(combined_text),
            'tissue': self.extract_tissue(combined_text),
            'sequencing_platform': self.extract_platform(combined_text),
            'experiment_design': metadata.get('resource_type', {}).get('title'),
            'sample_type': None,
            'summary': description[:500] if description else None,
            'citation_count': None,
            'publication_date': publication_date,
            'submission_date': created,
            'last_update_date': updated,
            'contact_name': contact_name,
            'contact_email': None,
            'contact_institute': contact_affiliation,
            'data_tier': data_tier,
            'tissue_location': None,
            'supplementary_information': json.dumps({
                'keywords': metadata.get('keywords', []),
                'license': metadata.get('license', {}).get('id'),
                'version': metadata.get('version'),
                'file_count': len(files),
                'total_size': sum(f.get('size', 0) for f in files),
                'total_size_mb': round(sum(f.get('size', 0) for f in files) / 1024 / 1024, 2),
                'communities': [c.get('id') for c in metadata.get('communities', [])],
                'related_identifiers': metadata.get('related_identifiers', []),
                'grants': metadata.get('grants', [])
            }, ensure_ascii=False)
        }
        
        return parsed_data
    
    def create_dataframe(self) -> pd.DataFrame:
        """将所有记录转换为DataFrame"""
        logger.info("\n" + "=" * 80)
        logger.info("开始解析记录为结构化数据...")
        logger.info("=" * 80)
        
        parsed_records = []
        total = len(self.all_records)
        
        for i, record in enumerate(self.all_records):
            try:
                parsed = self.parse_single_record(record)
                parsed_records.append(parsed)
                
                if (i + 1) % 100 == 0:
                    logger.info(f"已解析 {i + 1}/{total} 条记录 ({(i+1)/total*100:.1f}%)")
            except Exception as e:
                logger.error(f"解析记录出错 {record.get('id')}: {str(e)}")
                continue
        
        df = pd.DataFrame(parsed_records)
        logger.info(f"✓ 成功创建DataFrame，共 {len(df)} 条记录")
        logger.info("=" * 80)
        
        return df
    
    def filter_scrna_records(self, df: pd.DataFrame) -> pd.DataFrame:
        """进一步过滤，确保是单细胞RNA测序数据"""
        logger.info("\n" + "=" * 80)
        logger.info("开始过滤单细胞RNA测序记录...")
        logger.info("=" * 80)
        
        # 单细胞相关关键词
        scrna_keywords = [
            'single cell', 'single-cell', 'scRNA', 'sc-RNA', 'sc RNA',
            'scrna-seq', 'single cell RNA', 'single nucleus', 'sn-RNA',
            'snrna-seq', 'single-nucleus', '10x', 'drop-seq', 'smart-seq'
        ]
        
        def contains_scrna_keyword(row):
            text = f"{row['title']} {row['summary']}".lower()
            return any(keyword.lower() in text for keyword in scrna_keywords)
        
        # 排除关键词
        exclusion_keywords = [
            'review article', 'review only', 'protocol paper', 
            'method paper', 'software paper', 'protocol only',
            'methodology', 'software only', 'tool only'
        ]
        
        def not_excluded(row):
            text = f"{row['title']} {row['summary']}".lower()
            return not any(keyword in text for keyword in exclusion_keywords)
        
        # 应用过滤
        mask = df.apply(contains_scrna_keyword, axis=1)
        logger.info(f"  包含单细胞关键词: {mask.sum()}/{len(df)} 条记录")
        
        mask = mask & df.apply(not_excluded, axis=1)
        logger.info(f"  排除非数据记录后: {mask.sum()}/{len(df)} 条记录")
        
        filtered_df = df[mask].copy()
        
        logger.info(f"✓ 过滤完成，保留 {len(filtered_df)}/{len(df)} 条记录")
        logger.info("=" * 80)
        
        return filtered_df
    
    def save_results(self, df: pd.DataFrame, prefix: str = 'zenodo_scrna'):
        """保存结果到多种格式"""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # 保存CSV
        csv_file = f"{prefix}_{timestamp}.csv"
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        logger.info(f"✓ 已保存CSV文件: {csv_file}")
        
        # 保存Excel
        excel_file = f"{prefix}_{timestamp}.xlsx"
        try:
            df.to_excel(excel_file, index=False, engine='openpyxl')
            logger.info(f"✓ 已保存Excel文件: {excel_file}")
        except Exception as e:
            logger.warning(f"保存Excel文件失败: {str(e)}")
        
        # 保存JSON
        json_file = f"{prefix}_{timestamp}.json"
        df.to_json(json_file, orient='records', indent=2, force_ascii=False)
        logger.info(f"✓ 已保存JSON文件: {json_file}")
        
        # 保存原始记录pickle
        pickle_file = f"{prefix}_raw_records_{timestamp}.pkl"
        with open(pickle_file, 'wb') as f:
            pickle.dump(self.all_records, f)
        logger.info(f"✓ 已保存原始记录pickle文件: {pickle_file}")
        
        # 保存统计信息
        stats_file = f"{prefix}_statistics_{timestamp}.txt"
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("Zenodo单细胞RNA测序数据收集统计\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"收集时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总记录数: {len(df)}\n\n")
            
            if len(df) > 0:
                f.write("=" * 80 + "\n")
                f.write("疾病类型分布:\n")
                f.write("-" * 80 + "\n")
                f.write(df['disease_general'].value_counts().to_string())
                f.write("\n\n")
                
                f.write("=" * 80 + "\n")
                f.write("测序平台分布:\n")
                f.write("-" * 80 + "\n")
                platform_counts = df['sequencing_platform'].value_counts()
                f.write(platform_counts.to_string())
                f.write("\n\n")
                
                f.write("=" * 80 + "\n")
                f.write("组织类型分布 (Top 20):\n")
                f.write("-" * 80 + "\n")
                tissue_counts = df['tissue'].value_counts().head(20)
                f.write(tissue_counts.to_string())
                f.write("\n\n")
                
                f.write("=" * 80 + "\n")
                f.write("开放状态分布:\n")
                f.write("-" * 80 + "\n")
                f.write(df['open_status'].value_counts().to_string())
                f.write("\n\n")
                
                f.write("=" * 80 + "\n")
                f.write("数据层级分布:\n")
                f.write("-" * 80 + "\n")
                f.write(df['data_tier'].value_counts().to_string())
                f.write("\n\n")
                
                f.write("=" * 80 + "\n")
                f.write("其他统计:\n")
                f.write("-" * 80 + "\n")
                f.write(f"包含PubMed ID的记录: {df['pubmed'].notna().sum()}\n")
                f.write(f"包含DOI的记录: {df['doi'].notna().sum()}\n")
                f.write(f"开放访问记录: {(df['open_status'] == 'Open').sum()}\n")
                f.write(f"包含种族信息: {df['ethnicity'].notna().sum()}\n")
                f.write(f"包含性别信息: {df['sex'].notna().sum()}\n")
        
        logger.info(f"✓ 已保存统计文件: {stats_file}")


def main():
    """主函数：执行完整的数据收集流程"""
    
    print("\n" + "=" * 80)
    print(" " * 20 + "Zenodo单细胞RNA测序数据收集工具")
    print("=" * 80 + "\n")
    
    # ==================== 配置区域 ====================
    # 使用您配置的代理服务器
    USE_PROXY = True
    PROXY_URL = "http://proxy.mornai.cn:7890"
    VERIFY_SSL = True  # 通常应该保持True以确保安全
    ACCESS_TOKEN = None  # 如果有Zenodo API token可以填写
    MAX_RECORDS = 10000  # 最大收集记录数
    # ==================================================
    
    # 设置环境变量（作为备份方案）
    os.environ['HTTP_PROXY'] = PROXY_URL
    os.environ['HTTPS_PROXY'] = PROXY_URL
    
    logger.info(f"配置信息:")
    logger.info(f"  代理服务器: {PROXY_URL}")
    logger.info(f"  SSL验证: {VERIFY_SSL}")
    logger.info(f"  最大记录数: {MAX_RECORDS}")
    logger.info("")
    
    # 初始化收集器
    collector = ZenodoScRNASeqCollector(
        access_token=ACCESS_TOKEN,
        use_proxy=USE_PROXY,
        proxy_url=PROXY_URL,
        verify_ssl=VERIFY_SSL
    )
    
    # 测试连接
    if not collector.test_connection():
        logger.error("\n" + "=" * 80)
        logger.error("❌ 连接测试失败！")
        logger.error("=" * 80)
        logger.error("可能的原因:")
        logger.error("  1. 代理服务器 proxy.mornai.cn:7890 无法访问")
        logger.error("  2. 代理服务器需要认证")
        logger.error("  3. 防火墙阻止了连接")
        logger.error("  4. Zenodo服务暂时不可用")
        logger.error("")
        logger.error("建议:")
        logger.error("  1. 测试代理: curl -x http://proxy.mornai.cn:7890 https://zenodo.org")
        logger.error("  2. 检查代理配置是否正确")
        logger.error("  3. 尝试不使用代理（设置 USE_PROXY=False）")
        logger.error("=" * 80 + "\n")
        
        # 创建空DataFrame以避免错误
        df_empty = pd.DataFrame(columns=[
            'id', 'sample_id', 'title', 'disease_general', 'disease', 'pubmed', 
            'doi', 'source_database', 'access_link', 'open_status', 'ethnicity', 
            'sex', 'tissue', 'sequencing_platform', 'experiment_design', 
            'sample_type', 'summary', 'citation_count', 'publication_date', 
            'submission_date', 'last_update_date', 'contact_name', 'contact_email', 
            'contact_institute', 'data_tier', 'tissue_location', 'supplementary_information'
        ])
        collector.save_results(df_empty, prefix='zenodo_scrna_empty')
        return df_empty, df_empty
    
    try:
        # 步骤1: 搜索记录
        logger.info("\n" + "=" * 80)
        logger.info("步骤1/4: 搜索Zenodo单细胞测序记录")
        logger.info("=" * 80)
        records = collector.search_records(max_records=MAX_RECORDS)
        
        if not records:
            logger.warning("未找到任何记录")
            df_empty = pd.DataFrame(columns=[
                'id', 'sample_id', 'title', 'disease_general', 'disease', 'pubmed', 
                'doi', 'source_database', 'access_link', 'open_status', 'ethnicity', 
                'sex', 'tissue', 'sequencing_platform', 'experiment_design', 
                'sample_type', 'summary', 'citation_count', 'publication_date', 
                'submission_date', 'last_update_date', 'contact_name', 'contact_email', 
                'contact_institute', 'data_tier', 'tissue_location', 'supplementary_information'
            ])
            collector.save_results(df_empty, prefix='zenodo_scrna_no_results')
            return df_empty, df_empty
        
        # 步骤2: 丰富元数据
        logger.info("\n" + "=" * 80)
        logger.info("步骤2/4: 获取详细元数据")
        logger.info("=" * 80)
        collector.enrich_metadata()
        
        # 步骤3: 创建DataFrame
        logger.info("\n" + "=" * 80)
        logger.info("步骤3/4: 解析记录并创建数据框")
        logger.info("=" * 80)
        df_raw = collector.create_dataframe()
        
        # 保存原始数据框
        collector.save_results(df_raw, prefix='zenodo_scrna_raw')
        
        # 步骤4: 过滤单细胞数据
        logger.info("\n" + "=" * 80)
        logger.info("步骤4/4: 过滤单细胞RNA测序数据")
        logger.info("=" * 80)
        df_filtered = collector.filter_scrna_records(df_raw)
        
        # 保存过滤后的数据
        collector.save_results(df_filtered, prefix='zenodo_scrna_filtered')
        
        # 打印最终统计
        print("\n" + "=" * 80)
        print(" " * 30 + "收集完成!")
        print("=" * 80)
        print(f"原始记录数: {len(df_raw)}")
        print(f"过滤后记录数: {len(df_filtered)}")
        
        if len(df_filtered) > 0:
            print(f"包含PubMed ID: {df_filtered['pubmed'].notna().sum()}")
            print(f"包含DOI: {df_filtered['doi'].notna().sum()}")
            print(f"开放访问: {(df_filtered['open_status'] == 'Open').sum()}")
            print(f"包含组织信息: {df_filtered['tissue'].notna().sum()}")
            print(f"包含平台信息: {df_filtered['sequencing_platform'].notna().sum()}")
        
        print("=" * 80 + "\n")
        
        return df_raw, df_filtered
        
    except KeyboardInterrupt:
        logger.info("\n\n用户中断程序")
        return None, None
        
    except Exception as e:
        logger.error(f"\n程序执行出错: {str(e)}", exc_info=True)
        return None, None


if __name__ == "__main__":
    df_raw, df_filtered = main()
    
    if df_filtered is not None and len(df_filtered) > 0:
        print("\n" + "=" * 80)
        print("数据预览 (前5条):")
        print("=" * 80)
        print(df_filtered[['id', 'title', 'disease_general', 'tissue', 'sequencing_platform', 'open_status']].head())
        
        print("\n" + "=" * 80)
        print("DataFrame信息:")
        print("=" * 80)
        print(df_filtered.info())
        
        print("\n" + "=" * 80)
        print("数据已保存到文件，请查看当前目录")
        print("=" * 80)
    else:
        print("\n未能收集到数据")