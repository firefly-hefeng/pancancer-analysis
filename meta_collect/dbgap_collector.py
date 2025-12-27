#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dbGaP 单细胞测序数据检索与整合工具 v2.0
改进版：使用多种策略从NCBI获取dbGaP单细胞数据
"""

import requests
import xml.etree.ElementTree as ET
import pandas as pd
import time
import os
import json
from datetime import datetime
from pathlib import Path
import re
from typing import Dict, List, Optional, Tuple
import logging
from urllib.parse import quote
from collections import defaultdict

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('dbgap_scraper_v2.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class dbGaPScraperV2:
    """dbGaP 数据库爬取器 - 改进版"""
    
    def __init__(self, output_dir="dbgap_metadata_v2", email="your_email@example.com"):
        """
        初始化爬取器
        
        Args:
            output_dir: 输出目录路径
            email: 必需的联系邮箱
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # NCBI API配置
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.email = email
        self.api_key = None  # 如有API key请填入
        
        # 请求间隔（秒）
        self.request_delay = 0.34 if not self.api_key else 0.1
        
        logger.info(f"初始化完成，输出目录: {self.output_dir}")
        logger.info(f"联系邮箱: {self.email}")
    
    def strategy_1_search_sra_for_dbgap(self) -> List[Dict]:
        """
        策略1: 通过SRA数据库搜索关联到dbGaP的单细胞研究
        
        Returns:
            研究信息列表
        """
        logger.info("=" * 70)
        logger.info("策略1: 通过SRA数据库搜索dbGaP相关的单细胞研究")
        logger.info("=" * 70)
        
        # 构建SRA搜索查询
        search_terms = [
            '("single cell RNA"[Strategy] OR "single-cell RNA"[Strategy])',
            'AND "Homo sapiens"[Organism]',
            'AND "dbgap"[Filter]'  # 限定有dbGaP链接的数据
        ]
        
        query = ' '.join(search_terms)
        logger.info(f"SRA搜索查询: {query}")
        
        # Step 1: 搜索SRA
        study_ids = self._search_ncbi_db('sra', query, max_results=1000)
        logger.info(f"✓ 在SRA中找到 {len(study_ids)} 个相关研究")
        
        if not study_ids:
            logger.warning("SRA中未找到符合条件的研究")
            return []
        
        # Step 2: 获取SRA详情并提取dbGaP链接
        studies = []
        for idx, sra_id in enumerate(study_ids[:100], 1):  # 先处理前100个
            logger.info(f"[{idx}/{min(len(study_ids), 100)}] 处理SRA记录: {sra_id}")
            
            study_info = self._fetch_sra_details(sra_id)
            if study_info and study_info.get('dbgap_accession'):
                studies.append(study_info)
                logger.info(f"  ✓ 找到dbGaP研究: {study_info['dbgap_accession']}")
            
            if idx % 10 == 0:
                self._save_json(studies, "strategy1_intermediate.json")
            
            time.sleep(self.request_delay)
        
        logger.info(f"✓ 策略1完成，找到 {len(studies)} 个dbGaP研究")
        self._save_json(studies, "strategy1_results.json")
        
        return studies
    
    def strategy_2_search_bioproject(self) -> List[Dict]:
        """
        策略2: 通过BioProject搜索dbGaP项目
        
        Returns:
            研究信息列表
        """
        logger.info("=" * 70)
        logger.info("策略2: 通过BioProject搜索dbGaP项目")
        logger.info("=" * 70)
        
        # 搜索BioProject
        search_terms = [
            '(single cell[Title] OR scRNA[Title] OR "single-cell"[Title])',
            'AND "Homo sapiens"[Organism]',
            'AND dbgap[Properties]'
        ]
        
        query = ' '.join(search_terms)
        logger.info(f"BioProject搜索查询: {query}")
        
        project_ids = self._search_ncbi_db('bioproject', query, max_results=500)
        logger.info(f"✓ 找到 {len(project_ids)} 个BioProject")
        
        studies = []
        for idx, project_id in enumerate(project_ids, 1):
            logger.info(f"[{idx}/{len(project_ids)}] 处理BioProject: {project_id}")
            
            study_info = self._fetch_bioproject_details(project_id)
            if study_info:
                studies.append(study_info)
                logger.info(f"  ✓ 成功解析项目")
            
            if idx % 10 == 0:
                self._save_json(studies, "strategy2_intermediate.json")
            
            time.sleep(self.request_delay)
        
        logger.info(f"✓ 策略2完成，找到 {len(studies)} 个研究")
        self._save_json(studies, "strategy2_results.json")
        
        return studies
    
    def strategy_3_direct_dbgap_portal(self) -> List[Dict]:
        """
        策略3: 直接访问dbGaP网页API
        
        Returns:
            研究信息列表
        """
        logger.info("=" * 70)
        logger.info("策略3: 直接访问dbGaP网页API")
        logger.info("=" * 70)
        
        # dbGaP提供了一个JSON格式的研究列表
        dbgap_api_url = "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetStudyListJson.cgi"
        
        try:
            logger.info(f"请求dbGaP API: {dbgap_api_url}")
            response = requests.get(dbgap_api_url, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            all_studies = data.get('studies', [])
            
            logger.info(f"✓ 获取到 {len(all_studies)} 个dbGaP研究")
            
            # 筛选单细胞相关研究
            sc_studies = []
            sc_keywords = [
                'single cell', 'single-cell', 'scRNA', 'scRNA-seq',
                '10x', 'droplet', 'drop-seq', 'chromium'
            ]
            
            for study in all_studies:
                study_text = f"{study.get('study_name', '')} {study.get('disease', '')}".lower()
                
                if any(keyword.lower() in study_text for keyword in sc_keywords):
                    # 标准化研究信息
                    study_info = self._parse_dbgap_json_study(study)
                    sc_studies.append(study_info)
                    logger.info(f"  ✓ 找到单细胞研究: {study.get('study_accession')}")
            
            logger.info(f"✓ 策略3完成，找到 {len(sc_studies)} 个单细胞研究")
            self._save_json(sc_studies, "strategy3_results.json")
            
            return sc_studies
            
        except Exception as e:
            logger.error(f"✗ 策略3失败: {str(e)}")
            return []
    
    def strategy_4_pubmed_linked_search(self) -> List[Dict]:
        """
        策略4: 通过PubMed文献搜索关联的dbGaP数据
        
        Returns:
            研究信息列表
        """
        logger.info("=" * 70)
        logger.info("策略4: 通过PubMed文献搜索dbGaP数据")
        logger.info("=" * 70)
        
        # 搜索相关文献
        search_terms = [
            '("single cell RNA sequencing"[Title/Abstract] OR "scRNA-seq"[Title/Abstract])',
            'AND "Homo sapiens"[Organism]',
            'AND "dbgap"[Text Word]'
        ]
        
        query = ' '.join(search_terms)
        logger.info(f"PubMed搜索查询: {query}")
        
        pmids = self._search_ncbi_db('pubmed', query, max_results=200)
        logger.info(f"✓ 找到 {len(pmids)} 篇相关文献")
        
        studies = []
        dbgap_accessions = set()
        
        for idx, pmid in enumerate(pmids, 1):
            logger.info(f"[{idx}/{len(pmids)}] 处理文献 PMID: {pmid}")
            
            # 获取文献详情
            pub_info = self._fetch_pubmed_details(pmid)
            
            # 提取dbGaP accession
            accessions = self._extract_dbgap_accessions(pub_info)
            
            for accession in accessions:
                if accession not in dbgap_accessions:
                    dbgap_accessions.add(accession)
                    
                    study_info = {
                        'dbgap_accession': accession,
                        'source_pmid': pmid,
                        'source_strategy': 'pubmed_linked',
                        'publication': pub_info
                    }
                    
                    studies.append(study_info)
                    logger.info(f"  ✓ 提取到dbGaP accession: {accession}")
            
            if idx % 10 == 0:
                self._save_json(studies, "strategy4_intermediate.json")
            
            time.sleep(self.request_delay)
        
        logger.info(f"✓ 策略4完成，找到 {len(studies)} 个dbGaP研究")
        self._save_json(studies, "strategy4_results.json")
        
        return studies
    
    def _search_ncbi_db(self, db: str, query: str, max_results: int = 500) -> List[str]:
        """通用NCBI数据库搜索"""
        params = {
            'db': db,
            'term': query,
            'retmax': max_results,
            'retmode': 'json',
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        try:
            response = requests.get(f"{self.base_url}/esearch.fcgi", params=params, timeout=30)
            response.raise_for_status()
            
            result = response.json()
            id_list = result.get('esearchresult', {}).get('idlist', [])
            
            return id_list
            
        except Exception as e:
            logger.error(f"搜索{db}失败: {str(e)}")
            return []
    
    def _fetch_sra_details(self, sra_id: str) -> Optional[Dict]:
        """获取SRA详情并提取dbGaP信息"""
        try:
            params = {
                'db': 'sra',
                'id': sra_id,
                'retmode': 'xml',
                'email': self.email
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(f"{self.base_url}/efetch.fcgi", params=params, timeout=30)
            response.raise_for_status()
            
            root = ET.fromstring(response.content)
            
            # 提取关键信息
            study_info = {
                'sra_id': sra_id,
                'source_strategy': 'sra_linked',
                'dbgap_accession': None,
                'bioproject': None,
                'biosample': None,
                'title': '',
                'organism': '',
                'strategy': '',
                'platform': '',
                'raw_xml': ET.tostring(root, encoding='unicode')[:5000]  # 保留部分XML
            }
            
            # 提取标题
            title_elem = root.find('.//TITLE')
            if title_elem is not None:
                study_info['title'] = title_elem.text
            
            # 提取BioProject
            bioproject_elem = root.find('.//EXTERNAL_ID[@namespace="BioProject"]')
            if bioproject_elem is not None:
                study_info['bioproject'] = bioproject_elem.text
            
            # 提取dbGaP accession (通常在XREF_LINK中)
            for xref in root.findall('.//XREF_LINK'):
                db_elem = xref.find('DB')
                id_elem = xref.find('ID')
                if db_elem is not None and id_elem is not None:
                    if db_elem.text == 'dbGaP':
                        study_info['dbgap_accession'] = id_elem.text
            
            # 从描述中提取
            if not study_info['dbgap_accession']:
                description = root.find('.//STUDY_DESCRIPTION')
                if description is not None and description.text:
                    accessions = re.findall(r'phs\d{6}\.v\d+\.p\d+', description.text)
                    if accessions:
                        study_info['dbgap_accession'] = accessions[0]
            
            # 提取测序策略
            strategy_elem = root.find('.//LIBRARY_STRATEGY')
            if strategy_elem is not None:
                study_info['strategy'] = strategy_elem.text
            
            # 提取平台
            platform_elem = root.find('.//PLATFORM')
            if platform_elem is not None:
                for child in platform_elem:
                    study_info['platform'] = child.tag
                    break
            
            return study_info
            
        except Exception as e:
            logger.error(f"获取SRA {sra_id} 详情失败: {str(e)}")
            return None
    
    def _fetch_bioproject_details(self, project_id: str) -> Optional[Dict]:
        """获取BioProject详情"""
        try:
            params = {
                'db': 'bioproject',
                'id': project_id,
                'retmode': 'xml',
                'email': self.email
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(f"{self.base_url}/efetch.fcgi", params=params, timeout=30)
            response.raise_for_status()
            
            root = ET.fromstring(response.content)
            
            study_info = {
                'bioproject_id': project_id,
                'source_strategy': 'bioproject',
                'dbgap_accession': None,
                'title': '',
                'description': '',
                'organism': '',
                'data_type': []
            }
            
            # 提取项目信息
            project = root.find('.//Project')
            if project is not None:
                # 标题
                title_elem = project.find('.//ProjectDescr/Title')
                if title_elem is not None:
                    study_info['title'] = title_elem.text
                
                # 描述
                desc_elem = project.find('.//ProjectDescr/Description')
                if desc_elem is not None:
                    study_info['description'] = desc_elem.text
                
                # 提取dbGaP链接
                for ext_link in project.findall('.//ExternalLink'):
                    dbname = ext_link.find('dbname')
                    if dbname is not None and dbname.text == 'dbGaP':
                        label = ext_link.find('label')
                        if label is not None:
                            # 从label提取accession
                            accessions = re.findall(r'phs\d{6}', label.text)
                            if accessions:
                                study_info['dbgap_accession'] = accessions[0]
                
                # 从描述中提取
                if not study_info['dbgap_accession'] and study_info['description']:
                    accessions = re.findall(r'phs\d{6}\.v\d+\.p\d+', study_info['description'])
                    if accessions:
                        study_info['dbgap_accession'] = accessions[0]
            
            return study_info
            
        except Exception as e:
            logger.error(f"获取BioProject {project_id} 详情失败: {str(e)}")
            return None
    
    def _parse_dbgap_json_study(self, study_data: Dict) -> Dict:
        """解析dbGaP JSON格式的研究数据"""
        return {
            'dbgap_accession': study_data.get('study_accession', ''),
            'source_strategy': 'dbgap_api',
            'title': study_data.get('study_name', ''),
            'disease': study_data.get('disease', ''),
            'study_type': study_data.get('study_type', ''),
            'participants': study_data.get('number_of_participants', ''),
            'registration_date': study_data.get('registration_date', ''),
            'raw_data': study_data
        }
    
    def _fetch_pubmed_details(self, pmid: str) -> Dict:
        """获取PubMed文献详情"""
        try:
            params = {
                'db': 'pubmed',
                'id': pmid,
                'retmode': 'xml',
                'email': self.email
            }
            
            if self.api_key:
                params['api_key'] = self.api_key
            
            response = requests.get(f"{self.base_url}/efetch.fcgi", params=params, timeout=30)
            response.raise_for_status()
            
            root = ET.fromstring(response.content)
            
            article = root.find('.//Article')
            if article is None:
                return {}
            
            pub_info = {
                'pmid': pmid,
                'title': self._get_xml_text(article, './/ArticleTitle'),
                'journal': self._get_xml_text(article, './/Journal/Title'),
                'pub_date': self._get_xml_text(article, './/PubDate/Year'),
                'abstract': self._get_xml_text(article, './/Abstract/AbstractText'),
                'raw_xml': ET.tostring(article, encoding='unicode')[:5000]
            }
            
            return pub_info
            
        except Exception as e:
            logger.error(f"获取PMID {pmid} 失败: {str(e)}")
            return {}
    
    def _extract_dbgap_accessions(self, pub_info: Dict) -> List[str]:
        """从文献信息中提取dbGaP accession号"""
        accessions = []
        
        text = f"{pub_info.get('title', '')} {pub_info.get('abstract', '')} {pub_info.get('raw_xml', '')}"
        
        # 匹配 phs000000.v1.p1 格式
        matches = re.findall(r'phs\d{6}\.v\d+\.p\d+', text)
        accessions.extend(matches)
        
        # 匹配简化格式 phs000000
        matches = re.findall(r'phs\d{6}', text)
        accessions.extend(matches)
        
        return list(set(accessions))
    
    def _get_xml_text(self, parent: ET.Element, xpath: str) -> str:
        """从XML元素中提取文本"""
        elem = parent.find(xpath)
        return elem.text if elem is not None and elem.text else ""
    
    def enrich_dbgap_studies(self, studies: List[Dict]) -> List[Dict]:
        """
        丰富dbGaP研究信息
        通过accession号获取完整的研究详情
        """
        logger.info("=" * 70)
        logger.info(f"开始丰富 {len(studies)} 个研究的详细信息...")
        logger.info("=" * 70)
        
        enriched_studies = []
        unique_accessions = {}
        
        # 去重
        for study in studies:
            acc = study.get('dbgap_accession')
            if acc and acc not in unique_accessions:
                unique_accessions[acc] = study
        
        logger.info(f"去重后剩余 {len(unique_accessions)} 个唯一研究")
        
        for idx, (accession, study) in enumerate(unique_accessions.items(), 1):
            logger.info(f"[{idx}/{len(unique_accessions)}] 丰富研究: {accession}")
            
            try:
                # 通过dbGaP网页获取详细信息
                details = self._fetch_dbgap_web_details(accession)
                
                if details:
                    # 合并信息
                    enriched = {**study, **details}
                    enriched_studies.append(enriched)
                    logger.info(f"  ✓ 成功获取详细信息")
                else:
                    enriched_studies.append(study)
                    logger.info(f"  ⚠ 使用原始信息")
                
                time.sleep(self.request_delay)
                
            except Exception as e:
                logger.error(f"  ✗ 失败: {str(e)}")
                enriched_studies.append(study)
        
        logger.info(f"✓ 信息丰富完成")
        self._save_json(enriched_studies, "enriched_studies.json")
        
        return enriched_studies
    
    def _fetch_dbgap_web_details(self, accession: str) -> Optional[Dict]:
        """
        从dbGaP网页获取研究详情
        """
        # 清理accession（去除版本号）
        base_accession = accession.split('.')[0]
        
        # dbGaP study页面URL
        url = f"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={base_accession}"
        
        try:
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            }
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            
            html = response.text
            
            details = {
                'web_url': url,
                'html_fetched': True
            }
            
            # 简单的HTML解析提取关键信息（可以使用BeautifulSoup进一步优化）
            # 这里提供基本的正则提取
            
            # 提取PI信息
            pi_match = re.search(r'Principal Investigator.*?<td[^>]*>(.*?)</td>', html, re.DOTALL)
            if pi_match:
                details['pi'] = re.sub(r'<[^>]+>', '', pi_match.group(1)).strip()
            
            # 提取机构
            inst_match = re.search(r'Institution.*?<td[^>]*>(.*?)</td>', html, re.DOTALL)
            if inst_match:
                details['institution'] = re.sub(r'<[^>]+>', '', inst_match.group(1)).strip()
            
            return details
            
        except Exception as e:
            logger.error(f"获取{accession}网页详情失败: {str(e)}")
            return None
    
    def create_integrated_tables(self, studies: List[Dict]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        创建整合的series和sample表
        """
        logger.info("=" * 70)
        logger.info("创建整合数据表...")
        logger.info("=" * 70)
        
        series_records = []
        
        for study in studies:
            # 解析疾病信息
            disease_text = study.get('disease', '')
            disease_general = self._classify_disease(disease_text)
            
            # 获取第一篇相关文献
            pmid = study.get('source_pmid', '')
            if not pmid and study.get('publication'):
                pmid = study.get('publication', {}).get('pmid', '')
            
            series_record = {
                'id': study.get('dbgap_accession', ''),
                'title': study.get('title', '')[:200],
                'disease_general': disease_general,
                'disease': disease_text,
                'pubmed': pmid,
                'source_database': 'dbGaP',
                'access_link': study.get('web_url', f"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={study.get('dbgap_accession', '').split('.')[0]}"),
                'open_status': 'Controlled Access',
                'ethnicity': self._extract_field(study, ['ethnicity', 'race', 'population']),
                'sex': self._extract_field(study, ['sex', 'gender']),
                'tissue': self._extract_field(study, ['tissue', 'sample_type', 'biosample']),
                'sequencing_platform': study.get('platform', self._extract_platform(study)),
                'experiment_design': 'single-cell RNA-seq',
                'sample_type': self._determine_sample_type(study),
                'summary': study.get('description', '')[:500],
                'citation_count': '',
                'publication_date': study.get('publication', {}).get('pub_date', ''),
                'submission_date': study.get('registration_date', ''),
                'last_update_date': '',
                'contact_name': study.get('pi', ''),
                'contact_email': '',
                'contact_institute': study.get('institution', ''),
                'data_tier': 'raw',
                'supplementary_information': json.dumps({
                    'source_strategy': study.get('source_strategy', ''),
                    'bioproject': study.get('bioproject', ''),
                    'sra_id': study.get('sra_id', ''),
                    'participants': study.get('participants', '')
                }, ensure_ascii=False)
            }
            
            series_records.append(series_record)
        
        series_df = pd.DataFrame(series_records)
        
        # 保存Series表
        series_file = self.output_dir / "series_table.tsv"
        series_df.to_csv(series_file, sep='\t', index=False)
        logger.info(f"✓ Series表已保存: {series_file}")
        logger.info(f"  包含 {len(series_df)} 条记录")
        
        # 创建样本表骨架（实际样本信息需要授权访问）
        sample_df = self._create_sample_template(series_df)
        
        sample_file = self.output_dir / "sample_table_template.tsv"
        sample_df.to_csv(sample_file, sep='\t', index=False)
        logger.info(f"✓ 样本表模板已保存: {sample_file}")
        
        return series_df, sample_df
    
    def _classify_disease(self, disease_text: str) -> str:
        """疾病分类"""
        if not disease_text:
            return 'unknown'
        
        disease_lower = disease_text.lower()
        
        categories = {
            'cancer': ['cancer', 'tumor', 'carcinoma', 'leukemia', 'lymphoma', 'melanoma'],
            'neurodegenerative': ['alzheimer', 'parkinson', 'dementia', 'als', 'huntington'],
            'autoimmune': ['lupus', 'arthritis', 'diabetes', 'sclerosis', 'crohn'],
            'cardiovascular': ['heart', 'cardiovascular', 'coronary', 'cardiac'],
            'infectious': ['covid', 'virus', 'infection', 'bacterial', 'sepsis'],
            'metabolic': ['diabetes', 'obesity', 'metabolic'],
            'respiratory': ['asthma', 'copd', 'lung', 'pulmonary']
        }
        
        for category, keywords in categories.items():
            if any(kw in disease_lower for kw in keywords):
                return category
        
        return 'other'
    
    def _extract_field(self, study: Dict, field_names: List[str]) -> str:
        """从研究数据中提取字段"""
        for field in field_names:
            if field in study and study[field]:
                return str(study[field])
        
        # 从描述中提取
        text = f"{study.get('title', '')} {study.get('description', '')}".lower()
        
        # 根据字段类型进行模式匹配
        if 'ethnicity' in field_names or 'race' in field_names:
            patterns = ['european', 'asian', 'african', 'hispanic', 'caucasian']
            for pattern in patterns:
                if pattern in text:
                    return pattern
        
        return 'unknown'
    
    def _extract_platform(self, study: Dict) -> str:
        """提取测序平台"""
        text = f"{study.get('title', '')} {study.get('description', '')} {study.get('platform', '')}".lower()
        
        platforms = {
            '10x Genomics': ['10x', 'chromium', '10x genomics'],
            'Illumina': ['illumina', 'nextseq', 'hiseq', 'novaseq'],
            'Smart-seq2': ['smart-seq', 'smartseq'],
            'Drop-seq': ['drop-seq', 'dropseq'],
            'inDrop': ['indrop']
        }
        
        for platform, keywords in platforms.items():
            if any(kw in text for kw in keywords):
                return platform
        
        return 'unknown'
    
    def _determine_sample_type(self, study: Dict) -> str:
        """判断样本类型"""
        text = f"{study.get('title', '')} {study.get('disease', '')}".lower()
        
        if any(word in text for word in ['healthy', 'normal', 'control']):
            if any(word in text for word in ['disease', 'patient', 'tumor', 'cancer']):
                return 'mixed'
            return 'healthy'
        elif any(word in text for word in ['disease', 'patient', 'tumor', 'cancer']):
            return 'disease'
        
        return 'unknown'
    
    def _create_sample_template(self, series_df: pd.DataFrame) -> pd.DataFrame:
        """创建样本表模板"""
        sample_records = []
        
        for _, row in series_df.iterrows():
            # 为每个研究创建一个示例样本记录
            sample_record = {
                'id': f"{row['id']}_sample_template",
                'sample_id': 'Requires authorized access',
                'title': row['title'],
                'disease_general': row['disease_general'],
                'disease': row['disease'],
                'ethnicity': row['ethnicity'],
                'sex': row['sex'],
                'tissue_location': row['tissue'],
                'sequencing_platform': row['sequencing_platform'],
                'experiment_design': row['experiment_design'],
                'sample_type': row['sample_type']
            }
            
            sample_records.append(sample_record)
        
        return pd.DataFrame(sample_records)
    
    def _save_json(self, data: any, filename: str):
        """保存JSON文件"""
        filepath = self.output_dir / filename
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        logger.info(f"  💾 已保存: {filepath}")
    
    def generate_summary_report(self, series_df: pd.DataFrame):
        """生成汇总报告"""
        logger.info("=" * 70)
        logger.info("生成汇总报告...")
        logger.info("=" * 70)
        
        report = {
            'generation_time': datetime.now().isoformat(),
            'total_studies': len(series_df),
            'statistics': {
                'disease_distribution': series_df['disease_general'].value_counts().to_dict(),
                'platform_distribution': series_df['sequencing_platform'].value_counts().to_dict(),
                'sample_type_distribution': series_df['sample_type'].value_counts().to_dict(),
                'access_status': series_df['open_status'].value_counts().to_dict()
            },
            'data_sources': {
                'databases': ['dbGaP', 'SRA', 'BioProject', 'PubMed'],
                'strategies_used': [
                    'SRA-linked search',
                    'BioProject search',
                    'dbGaP API',
                    'PubMed literature mining'
                ]
            }
        }
        
        report_file = self.output_dir / "summary_report.json"
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        logger.info(f"✓ 汇总报告已保存: {report_file}")
        logger.info("\n" + "=" * 70)
        logger.info("📊 数据统计:")
        logger.info(f"  总研究数: {report['total_studies']}")
        logger.info(f"  疾病分布: {report['statistics']['disease_distribution']}")
        logger.info(f"  平台分布: {report['statistics']['platform_distribution']}")
        logger.info("=" * 70)


def main():
    """主函数 - 运行所有策略"""
    print("=" * 80)
    print(" " * 15 + "dbGaP 单细胞测序数据检索工具 v2.0")
    print(" " * 20 + "多策略整合版 - 2025-12-02")
    print("=" * 80)
    print()
    
    # 获取用户邮箱
    email = input("请输入您的邮箱地址（NCBI要求）: ").strip()
    if not email:
        email = "user@example.com"
        print(f"使用默认邮箱: {email}")
    
    api_key = input("如有NCBI API key请输入（可选，直接回车跳过）: ").strip()
    
    print()
    
    # 初始化爬取器
    scraper = dbGaPScraperV2(output_dir="dbgap_metadata_v2", email=email)
    if api_key:
        scraper.api_key = api_key
        scraper.request_delay = 0.1
        logger.info("已设置API key，使用更快的请求速率")
    
    all_studies = []
    
    # 运行策略1: SRA搜索
    print("\n" + "🔍 运行策略1: 通过SRA搜索...")
    try:
        studies_1 = scraper.strategy_1_search_sra_for_dbgap()
        all_studies.extend(studies_1)
        print(f"✅ 策略1完成，找到 {len(studies_1)} 个研究\n")
    except Exception as e:
        print(f"❌ 策略1失败: {str(e)}\n")
    
    # 运行策略2: BioProject搜索
    print("🔍 运行策略2: 通过BioProject搜索...")
    try:
        studies_2 = scraper.strategy_2_search_bioproject()
        all_studies.extend(studies_2)
        print(f"✅ 策略2完成，找到 {len(studies_2)} 个研究\n")
    except Exception as e:
        print(f"❌ 策略2失败: {str(e)}\n")
    
    # 运行策略3: dbGaP API
    print("🔍 运行策略3: 访问dbGaP API...")
    try:
        studies_3 = scraper.strategy_3_direct_dbgap_portal()
        all_studies.extend(studies_3)
        print(f"✅ 策略3完成，找到 {len(studies_3)} 个研究\n")
    except Exception as e:
        print(f"❌ 策略3失败: {str(e)}\n")
    
    # 运行策略4: PubMed搜索
    print("🔍 运行策略4: 通过PubMed文献搜索...")
    try:
        studies_4 = scraper.strategy_4_pubmed_linked_search()
        all_studies.extend(studies_4)
        print(f"✅ 策略4完成，找到 {len(studies_4)} 个研究\n")
    except Exception as e:
        print(f"❌ 策略4失败: {str(e)}\n")
    
    if not all_studies:
        print("❌ 所有策略均未找到研究，程序退出")
        return
    
    print(f"\n✅ 合计找到 {len(all_studies)} 个研究记录（包含重复）\n")
    
    # 丰富研究信息
    print("📝 丰富研究详细信息...")
    enriched_studies = scraper.enrich_dbgap_studies(all_studies)
    print(f"✅ 信息丰富完成\n")
    
    # 创建数据表
    print("📊 创建标准化数据表...")
    series_df, sample_df = scraper.create_integrated_tables(enriched_studies)
    print(f"✅ 数据表创建完成\n")
    
    # 生成报告
    scraper.generate_summary_report(series_df)
    
    print("\n" + "=" * 80)
    print("🎉 所有任务完成！")
    print(f"📁 输出目录: {scraper.output_dir.absolute()}")
    print("=" * 80)
    
    print("\n📋 输出文件清单:")
    print("  策略结果:")
    print("    ├─ strategy1_results.json      (SRA搜索结果)")
    print("    ├─ strategy2_results.json      (BioProject搜索结果)")
    print("    ├─ strategy3_results.json      (dbGaP API结果)")
    print("    └─ strategy4_results.json      (PubMed搜索结果)")
    print("  整合数据:")
    print("    ├─ enriched_studies.json       (丰富后的研究信息)")
    print("    ├─ series_table.tsv            (Series数据表)")
    print("    ├─ sample_table_template.tsv   (样本表模板)")
    print("    ├─ summary_report.json         (汇总报告)")
    print("    └─ dbgap_scraper_v2.log        (运行日志)")
    print()
    
    print("⚠️  重要提示:")
    print("   1. 样本级别详细数据需要申请dbGaP授权访问")
    print("   2. 访问 https://dbgap.ncbi.nlm.nih.gov 申请DAR")
    print("   3. 使用 sra-toolkit 配合授权下载原始数据")
    print("   4. 当前表格提供研究级别的元数据概览")
    print()


if __name__ == "__main__":
    main()