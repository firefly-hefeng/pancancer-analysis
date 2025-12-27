#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NCBI BioProject/BioSample/SRA + ENA 人源单细胞数据收集系统
增强版：包含FASTQ下载链接、PubMed详情、ENA数据库检索
"""

import requests
import pandas as pd
import xml.etree.ElementTree as ET
from datetime import datetime
import time
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import sqlite3
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import warnings
warnings.filterwarnings('ignore')

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('ncbi_bioproject_sra_collection.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)

class NCBIBioProjectSRACollector:
    """NCBI BioProject/BioSample/SRA + ENA 人源单细胞数据收集器"""
    
    def __init__(self, output_dir='ncbi_bioproject_sra_data', email='your_email@example.com'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        self.email = email
        self.ncbi_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.ena_base = "https://www.ebi.ac.uk/ena/portal/api"
        
        # 初始化数据库
        self.db_path = self.output_dir / 'bioproject_sra_metadata.db'
        self.init_database()
        
        # 疾病关键词映射
        self.disease_keywords = {
            'cancer': ['cancer', 'tumor', 'carcinoma', 'oncology', 'malignant', 'neoplasm', 
                      'leukemia', 'lymphoma', 'melanoma', 'glioma', 'sarcoma'],
            'covid-19': ['covid', 'sars-cov-2', 'coronavirus', 'covid-19'],
            'diabetes': ['diabetes', 'diabetic', 'type 1 diabetes', 'type 2 diabetes'],
            'alzheimer': ['alzheimer', 'dementia', 'neurodegenerative'],
            'cardiovascular': ['cardiovascular', 'heart', 'cardiac', 'myocardial', 'coronary'],
            'autoimmune': ['autoimmune', 'lupus', 'rheumatoid', 'arthritis', 'sclerosis'],
            'infectious': ['infection', 'viral', 'bacterial', 'sepsis', 'influenza'],
            'neurological': ['neurological', 'brain', 'neural', 'nervous system', 'parkinson'],
            'immunological': ['immune', 'immunology', 'lymphocyte', 't cell', 'b cell'],
            'developmental': ['development', 'developmental', 'embryo', 'fetal'],
            'healthy': ['healthy', 'normal', 'control', 'wild-type', 'unaffected']
        }
        
        # 组织关键词
        self.tissue_keywords = [
            'pbmc', 'blood', 'brain', 'lung', 'liver', 'kidney', 'heart', 'skin',
            'bone marrow', 'lymph node', 'spleen', 'thymus', 'pancreas', 'intestine',
            'colon', 'stomach', 'muscle', 'adipose', 'breast', 'prostate', 'ovary',
            'testis', 'uterus', 'placenta', 'cord blood', 'tumor', 'cancer'
        ]
        
        logging.info(f"初始化NCBI BioProject/SRA + ENA收集器")
        logging.info(f"输出目录: {self.output_dir}")
        logging.info(f"Email: {self.email}")
    
    def init_database(self):
        """初始化SQLite数据库"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # BioProject表
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS bioprojects (
                bioproject_id TEXT PRIMARY KEY,
                accession TEXT,
                title TEXT,
                description TEXT,
                organism TEXT,
                project_type TEXT,
                data_type TEXT,
                submission_date TEXT,
                last_update_date TEXT,
                publications TEXT,
                submitter_organization TEXT,
                submitter_name TEXT,
                submitter_email TEXT,
                raw_xml TEXT,
                collection_date TEXT
            )
        ''')
        
        # BioSample表
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS biosamples (
                biosample_id TEXT PRIMARY KEY,
                accession TEXT,
                bioproject_id TEXT,
                title TEXT,
                organism TEXT,
                tissue TEXT,
                cell_type TEXT,
                disease TEXT,
                age TEXT,
                sex TEXT,
                ethnicity TEXT,
                development_stage TEXT,
                attributes TEXT,
                raw_xml TEXT,
                collection_date TEXT
            )
        ''')
        
        # SRA Study表
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS sra_studies (
                study_accession TEXT PRIMARY KEY,
                bioproject_accession TEXT,
                title TEXT,
                abstract TEXT,
                study_type TEXT,
                center_name TEXT,
                submission_date TEXT,
                publication_date TEXT,
                last_update_date TEXT,
                raw_xml TEXT,
                collection_date TEXT
            )
        ''')
        
        # SRA Experiment表
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS sra_experiments (
                experiment_accession TEXT PRIMARY KEY,
                study_accession TEXT,
                biosample_accession TEXT,
                title TEXT,
                library_name TEXT,
                library_strategy TEXT,
                library_source TEXT,
                library_selection TEXT,
                library_layout TEXT,
                platform TEXT,
                instrument_model TEXT,
                design_description TEXT,
                raw_xml TEXT,
                collection_date TEXT
            )
        ''')
        
        # SRA Run表 - 增强版，包含FASTQ下载链接
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS sra_runs (
                run_accession TEXT PRIMARY KEY,
                experiment_accession TEXT,
                biosample_accession TEXT,
                total_spots INTEGER,
                total_bases INTEGER,
                size_mb REAL,
                published_date TEXT,
                load_done_date TEXT,
                download_path TEXT,
                fastq_ftp TEXT,
                fastq_aspera TEXT,
                fastq_galaxy TEXT,
                fastq_md5 TEXT,
                ena_fastq_files TEXT,
                raw_xml TEXT,
                collection_date TEXT
            )
        ''')
        
        # PubMed文献表 - 新增
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS pubmed_articles (
                pmid TEXT PRIMARY KEY,
                title TEXT,
                abstract TEXT,
                authors TEXT,
                journal TEXT,
                publication_date TEXT,
                doi TEXT,
                pmc_id TEXT,
                citation_count INTEGER,
                mesh_terms TEXT,
                keywords TEXT,
                raw_xml TEXT,
                collection_date TEXT
            )
        ''')
        
        # ENA数据表 - 新增
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS ena_studies (
                study_accession TEXT PRIMARY KEY,
                secondary_accession TEXT,
                title TEXT,
                description TEXT,
                center_name TEXT,
                first_public TEXT,
                last_updated TEXT,
                study_type TEXT,
                sample_count INTEGER,
                experiment_count INTEGER,
                run_count INTEGER,
                tax_id INTEGER,
                scientific_name TEXT,
                fastq_files TEXT,
                collection_date TEXT
            )
        ''')
        
        # 最终整合表 - 增强版
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS final_series (
                id TEXT PRIMARY KEY,
                source_database TEXT,
                title TEXT,
                disease_general TEXT,
                disease TEXT,
                pubmed TEXT,
                pubmed_titles TEXT,
                pubmed_abstracts TEXT,
                access_link TEXT,
                open_status TEXT,
                ethnicity TEXT,
                sex TEXT,
                tissue TEXT,
                sequencing_platform TEXT,
                experiment_design TEXT,
                sample_type TEXT,
                summary TEXT,
                citation_count INTEGER,
                publication_date TEXT,
                submission_date TEXT,
                last_update_date TEXT,
                contact_name TEXT,
                contact_email TEXT,
                contact_institute TEXT,
                data_tier TEXT,
                tissue_location TEXT,
                supplementary_information TEXT,
                bioproject TEXT,
                biosample_list TEXT,
                run_count INTEGER,
                total_bases INTEGER,
                total_spots INTEGER,
                fastq_download_links TEXT,
                ena_accession TEXT,
                ena_fastq_links TEXT,
                metadata_completeness REAL,
                collection_date TEXT
            )
        ''')
        
        conn.commit()
        conn.close()
        logging.info("数据库初始化完成")
    
    def search_bioprojects(self) -> List[str]:
        """搜索所有人类单细胞BioProject"""
        logging.info("="*70)
        logging.info("搜索BioProject")
        logging.info("="*70)
        
        all_bp_ids = set()
        
        # 多种搜索策略
        search_terms = [
            # 核心搜索词
            'Homo sapiens[Organism] AND (single cell[All Fields] OR single-cell[All Fields])',
            'Homo sapiens[Organism] AND scRNA-seq[All Fields]',
            'Homo sapiens[Organism] AND single cell RNA sequencing[All Fields]',
            'Homo sapiens[Organism] AND single-cell RNA-seq[All Fields]',
            
            # 技术平台
            'Homo sapiens[Organism] AND 10x genomics[All Fields]',
            'Homo sapiens[Organism] AND 10X Chromium[All Fields]',
            'Homo sapiens[Organism] AND drop-seq[All Fields]',
            'Homo sapiens[Organism] AND Smart-seq[All Fields]',
            'Homo sapiens[Organism] AND SMART-seq2[All Fields]',
            'Homo sapiens[Organism] AND CEL-seq[All Fields]',
            'Homo sapiens[Organism] AND MARS-seq[All Fields]',
            'Homo sapiens[Organism] AND inDrop[All Fields]',
            
            # 单核测序
            'Homo sapiens[Organism] AND single nucleus[All Fields]',
            'Homo sapiens[Organism] AND snRNA-seq[All Fields]',
            'Homo sapiens[Organism] AND single-nucleus RNA-seq[All Fields]',
            
            # 多组学
            'Homo sapiens[Organism] AND CITE-seq[All Fields]',
            'Homo sapiens[Organism] AND multiome[All Fields]',
            'Homo sapiens[Organism] AND single cell multiomics[All Fields]',
            
            # 空间转录组
            'Homo sapiens[Organism] AND spatial transcriptomics[All Fields]',
            'Homo sapiens[Organism] AND Visium[All Fields]',
            'Homo sapiens[Organism] AND spatial RNA-seq[All Fields]',
        ]
        
        for i, term in enumerate(search_terms, 1):
            logging.info(f"\n策略 {i}/{len(search_terms)}")
            logging.info(f"搜索词: {term[:80]}...")
            
            try:
                search_url = f"{self.ncbi_base}/esearch.fcgi"
                params = {
                    'db': 'bioproject',
                    'term': term,
                    'retmax': 100000,
                    'retmode': 'json',
                    'email': self.email,
                    'usehistory': 'y'
                }
                
                response = requests.get(search_url, params=params, timeout=60)
                response.raise_for_status()
                data = response.json()
                
                ids = data.get('esearchresult', {}).get('idlist', [])
                count = data.get('esearchresult', {}).get('count', '0')
                
                logging.info(f"  找到 {count} 个结果，获取到 {len(ids)} 个ID")
                all_bp_ids.update(ids)
                
                time.sleep(0.5)
                
            except Exception as e:
                logging.error(f"  搜索失败: {e}")
                continue
        
        logging.info(f"\n总计找到 {len(all_bp_ids)} 个唯一BioProject")
        return list(all_bp_ids)
    
    def fetch_bioproject_details(self, bp_id: str) -> Optional[Dict]:
        """获取BioProject详细信息"""
        try:
            fetch_url = f"{self.ncbi_base}/efetch.fcgi"
            params = {
                'db': 'bioproject',
                'id': bp_id,
                'retmode': 'xml',
                'email': self.email
            }
            
            response = requests.get(fetch_url, params=params, timeout=30)
            response.raise_for_status()
            root = ET.fromstring(response.content)
            
            project = root.find('.//Project')
            if project is None:
                return None
            
            details = {
                'bioproject_id': bp_id,
                'accession': '',
                'title': '',
                'description': '',
                'organism': 'Homo sapiens',
                'project_type': '',
                'data_type': '',
                'submission_date': '',
                'last_update_date': '',
                'publications': [],
                'submitter_organization': '',
                'submitter_name': '',
                'submitter_email': '',
                'raw_xml': ET.tostring(root, encoding='unicode')
            }
            
            # Accession
            for acc in project.findall('.//ProjectID/ArchiveID'):
                details['accession'] = acc.get('accession', '')
                break
            
            # 标题和描述
            details['title'] = self.get_xml_text(project, './/Title') or \
                              self.get_xml_text(project, './/Name')
            details['description'] = self.get_xml_text(project, './/Description') or \
                                    self.get_xml_text(project, './/ProjectDescr/Description')
            
            # 项目类型
            proj_type = project.find('.//ProjectType')
            if proj_type is not None:
                for child in proj_type:
                    details['project_type'] = child.tag
                    break
            
            # 数据类型
            for data_type in project.findall('.//ProjectDataTypeSet/DataType'):
                if data_type.text:
                    details['data_type'] = data_type.text
                    break
            
            # 提交日期
            submission = project.find('.//Submission')
            if submission is not None:
                details['submission_date'] = submission.get('submitted', '') or \
                                           submission.get('submission_date', '')
                details['last_update_date'] = submission.get('last_update', '')
            
            # PubMed
            for pub in project.findall('.//Publication'):
                pub_id = pub.get('id', '')
                if pub_id:
                    details['publications'].append(pub_id)
                # 也尝试从子元素获取
                pubmed_id = self.get_xml_text(pub, './/DbType[@id="pubmed"]/../Id')
                if pubmed_id:
                    details['publications'].append(pubmed_id)
            
            # 提交者信息
            submitter = project.find('.//Submission/Organization')
            if submitter is not None:
                details['submitter_organization'] = self.get_xml_text(submitter, './/Name')
                
                contact = submitter.find('.//Contact')
                if contact is not None:
                    name_elem = contact.find('.//Name')
                    if name_elem is not None:
                        first_name = self.get_xml_text(name_elem, './/First')
                        last_name = self.get_xml_text(name_elem, './/Last')
                        details['submitter_name'] = f"{first_name} {last_name}".strip()
                    
                    details['submitter_email'] = self.get_xml_text(contact, './/Email')
            
            if details['accession']:
                self.save_bioproject(details)
                return details
            
            return None
            
        except Exception as e:
            logging.error(f"获取BioProject {bp_id} 失败: {e}")
            return None
    
    def save_bioproject(self, details: Dict):
        """保存BioProject到数据库"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO bioprojects
                (bioproject_id, accession, title, description, organism,
                 project_type, data_type, submission_date, last_update_date,
                 publications, submitter_organization, submitter_name,
                 submitter_email, raw_xml, collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                details['bioproject_id'],
                details['accession'],
                details['title'],
                details['description'],
                details['organism'],
                details['project_type'],
                details['data_type'],
                details['submission_date'],
                details['last_update_date'],
                json.dumps(details['publications']),
                details['submitter_organization'],
                details['submitter_name'],
                details['submitter_email'],
                details['raw_xml'],
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存BioProject失败: {e}")
    
    def fetch_linked_biosamples(self, bp_accession: str) -> List[str]:
        """获取关联的BioSample ID"""
        try:
            # 先通过accession搜索BioProject的ID
            search_url = f"{self.ncbi_base}/esearch.fcgi"
            params = {
                'db': 'bioproject',
                'term': f'{bp_accession}[Accession]',
                'retmode': 'json',
                'email': self.email
            }
            
            response = requests.get(search_url, params=params, timeout=30)
            data = response.json()
            
            bp_ids = data.get('esearchresult', {}).get('idlist', [])
            if not bp_ids:
                return []
            
            bp_id = bp_ids[0]
            
            # 获取链接的BioSample
            link_url = f"{self.ncbi_base}/elink.fcgi"
            params = {
                'dbfrom': 'bioproject',
                'db': 'biosample',
                'id': bp_id,
                'retmode': 'json',
                'email': self.email
            }
            
            response = requests.get(link_url, params=params, timeout=30)
            data = response.json()
            
            biosample_ids = []
            linksets = data.get('linksets', [])
            if linksets:
                linksetdbs = linksets[0].get('linksetdbs', [])
                for linksetdb in linksetdbs:
                    links = linksetdb.get('links', [])
                    biosample_ids.extend(links)
            
            return biosample_ids
            
        except Exception as e:
            logging.error(f"获取BioSample链接失败: {e}")
            return []
    
    def fetch_biosample_details(self, bs_id: str) -> Optional[Dict]:
        """获取BioSample详细信息"""
        try:
            fetch_url = f"{self.ncbi_base}/efetch.fcgi"
            params = {
                'db': 'biosample',
                'id': bs_id,
                'retmode': 'xml',
                'email': self.email
            }
            
            response = requests.get(fetch_url, params=params, timeout=30)
            response.raise_for_status()
            root = ET.fromstring(response.content)
            
            biosample = root.find('.//BioSample')
            if biosample is None:
                return None
            
            details = {
                'biosample_id': bs_id,
                'accession': biosample.get('accession', ''),
                'bioproject_id': '',
                'title': '',
                'organism': 'Homo sapiens',
                'tissue': '',
                'cell_type': '',
                'disease': '',
                'age': '',
                'sex': '',
                'ethnicity': '',
                'development_stage': '',
                'attributes': {},
                'raw_xml': ET.tostring(root, encoding='unicode')
            }
            
            # 标题
            description = biosample.find('.//Description')
            if description is not None:
                title = description.find('.//Title')
                if title is not None:
                    details['title'] = title.text or ''
            
            # 提取所有属性
            for attr in biosample.findall('.//Attribute'):
                attr_name = attr.get('attribute_name', attr.get('harmonized_name', ''))
                attr_value = attr.text or ''
                
                if attr_name and attr_value:
                    details['attributes'][attr_name] = attr_value
                    
                    # 映射到标准字段
                    attr_lower = attr_name.lower()
                    
                    if any(x in attr_lower for x in ['tissue', 'source']):
                        if not details['tissue']:
                            details['tissue'] = attr_value
                    
                    elif 'cell type' in attr_lower or 'cell_type' in attr_lower:
                        details['cell_type'] = attr_value
                    
                    elif any(x in attr_lower for x in ['disease', 'condition', 'phenotype']):
                        details['disease'] = attr_value
                    
                    elif 'age' in attr_lower:
                        details['age'] = attr_value
                    
                    elif 'sex' in attr_lower or 'gender' in attr_lower:
                        details['sex'] = attr_value
                    
                    elif any(x in attr_lower for x in ['race', 'ethnicity', 'ancestry']):
                        details['ethnicity'] = attr_value
                    
                    elif 'development' in attr_lower or 'developmental' in attr_lower:
                        details['development_stage'] = attr_value
            
            # BioProject链接
            for id_elem in biosample.findall('.//Id'):
                if id_elem.get('db') == 'BioProject':
                    details['bioproject_id'] = id_elem.text
                    break
            
            # 也从Links获取
            for link in biosample.findall('.//Links/Link'):
                if link.get('type') == 'entrez':
                    target = link.get('target', '')
                    if 'bioproject' in target.lower():
                        label = link.get('label', '')
                        if label:
                            details['bioproject_id'] = label
                            break
            
            if details['accession']:
                self.save_biosample(details)
                return details
            
            return None
            
        except Exception as e:
            logging.error(f"获取BioSample {bs_id} 失败: {e}")
            return None
    
    def save_biosample(self, details: Dict):
        """保存BioSample到数据库"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO biosamples
                (biosample_id, accession, bioproject_id, title, organism,
                 tissue, cell_type, disease, age, sex, ethnicity,
                 development_stage, attributes, raw_xml, collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                details['biosample_id'],
                details['accession'],
                details['bioproject_id'],
                details['title'],
                details['organism'],
                details['tissue'],
                details['cell_type'],
                details['disease'],
                details['age'],
                details['sex'],
                details['ethnicity'],
                details['development_stage'],
                json.dumps(details['attributes'], ensure_ascii=False),
                details['raw_xml'],
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存BioSample失败: {e}")
    
    def fetch_linked_sra_studies(self, bp_accession: str) -> List[str]:
        """获取关联的SRA Study"""
        try:
            # 搜索SRA
            search_url = f"{self.ncbi_base}/esearch.fcgi"
            params = {
                'db': 'sra',
                'term': f'{bp_accession}[BioProject]',
                'retmax': 10000,
                'retmode': 'json',
                'email': self.email
            }
            
            response = requests.get(search_url, params=params, timeout=30)
            data = response.json()
            
            sra_ids = data.get('esearchresult', {}).get('idlist', [])
            return sra_ids
            
        except Exception as e:
            logging.error(f"获取SRA链接失败: {e}")
            return []
    
    def fetch_sra_details_batch(self, sra_ids: List[str]) -> Dict:
        """批量获取SRA详细信息（增强版：包含FASTQ下载链接）"""
        if not sra_ids:
            return {'studies': [], 'experiments': [], 'runs': []}
        
        try:
            # 批量获取
            fetch_url = f"{self.ncbi_base}/efetch.fcgi"
            params = {
                'db': 'sra',
                'id': ','.join(sra_ids),
                'retmode': 'xml',
                'email': self.email
            }
            
            response = requests.get(fetch_url, params=params, timeout=120)
            response.raise_for_status()
            root = ET.fromstring(response.content)
            
            studies = []
            experiments = []
            runs = []
            
            # 解析每个EXPERIMENT_PACKAGE
            for exp_package in root.findall('.//EXPERIMENT_PACKAGE'):
                # Study信息
                study = exp_package.find('.//STUDY')
                if study is not None:
                    study_acc = study.get('accession', '')
                    if study_acc:
                        study_info = {
                            'study_accession': study_acc,
                            'bioproject_accession': '',
                            'title': self.get_xml_text(study, './/STUDY_TITLE'),
                            'abstract': self.get_xml_text(study, './/STUDY_ABSTRACT'),
                            'study_type': self.get_xml_text(study, './/STUDY_TYPE'),
                            'center_name': study.get('center_name', ''),
                            'submission_date': '',
                            'publication_date': '',
                            'last_update_date': '',
                            'raw_xml': ET.tostring(study, encoding='unicode')
                        }
                        
                        # 获取BioProject链接
                        for link in study.findall('.//STUDY_LINK'):
                            xref_link = link.find('.//XREF_LINK')
                            if xref_link is not None:
                                db = self.get_xml_text(xref_link, './/DB')
                                if db == 'bioproject':
                                    study_info['bioproject_accession'] = self.get_xml_text(xref_link, './/ID')
                        
                        studies.append(study_info)
                        self.save_sra_study(study_info)
                
                # Experiment信息
                experiment = exp_package.find('.//EXPERIMENT')
                if experiment is not None:
                    exp_acc = experiment.get('accession', '')
                    if exp_acc:
                        exp_info = {
                            'experiment_accession': exp_acc,
                            'study_accession': study_acc if study is not None else '',
                            'biosample_accession': '',
                            'title': self.get_xml_text(experiment, './/TITLE'),
                            'library_name': self.get_xml_text(experiment, './/LIBRARY_NAME'),
                            'library_strategy': self.get_xml_text(experiment, './/LIBRARY_STRATEGY'),
                            'library_source': self.get_xml_text(experiment, './/LIBRARY_SOURCE'),
                            'library_selection': self.get_xml_text(experiment, './/LIBRARY_SELECTION'),
                            'library_layout': '',
                            'platform': '',
                            'instrument_model': '',
                            'design_description': self.get_xml_text(experiment, './/DESIGN_DESCRIPTION'),
                            'raw_xml': ET.tostring(experiment, encoding='unicode')
                        }
                        
                        # Library layout
                        layout = experiment.find('.//LIBRARY_LAYOUT')
                        if layout is not None:
                            for child in layout:
                                exp_info['library_layout'] = child.tag
                                break
                        
                        # Platform
                        platform = experiment.find('.//PLATFORM')
                        if platform is not None:
                            for instrument in platform:
                                exp_info['platform'] = instrument.tag
                                exp_info['instrument_model'] = self.get_xml_text(instrument, './/INSTRUMENT_MODEL')
                                break
                        
                        # BioSample
                        sample = exp_package.find('.//SAMPLE')
                        if sample is not None:
                            exp_info['biosample_accession'] = sample.get('accession', '')
                        
                        experiments.append(exp_info)
                        self.save_sra_experiment(exp_info)
                
                # Run信息 - 增强版
                for run in exp_package.findall('.//RUN'):
                    run_acc = run.get('accession', '')
                    if run_acc:
                        run_info = {
                            'run_accession': run_acc,
                            'experiment_accession': exp_acc if experiment is not None else '',
                            'biosample_accession': '',
                            'total_spots': 0,
                            'total_bases': 0,
                            'size_mb': 0.0,
                            'published_date': run.get('published', ''),
                            'load_done_date': run.get('load_done', ''),
                            'download_path': '',
                            'fastq_ftp': '',
                            'fastq_aspera': '',
                            'fastq_galaxy': '',
                            'fastq_md5': '',
                            'ena_fastq_files': '',
                            'raw_xml': ET.tostring(run, encoding='unicode')
                        }
                        
                        # 样本信息
                        sample = exp_package.find('.//SAMPLE')
                        if sample is not None:
                            run_info['biosample_accession'] = sample.get('accession', '')
                        
                        # 数据量统计
                        try:
                            run_info['total_spots'] = int(run.get('total_spots', 0))
                            run_info['total_bases'] = int(run.get('total_bases', 0))
                        except:
                            pass
                        
                        # 文件大小
                        try:
                            run_info['size_mb'] = float(run.get('size', 0)) / (1024 * 1024)
                        except:
                            pass
                        
                        # 下载路径
                        sra_files = run.find('.//SRAFiles')
                        if sra_files is not None:
                            for sra_file in sra_files.findall('.//SRAFile'):
                                url = sra_file.get('url', '')
                                if url and 'sra' in url:
                                    run_info['download_path'] = url
                                    break
                        
                        # 从ENA获取FASTQ下载链接
                        ena_fastq = self.fetch_ena_fastq_links(run_acc)
                        if ena_fastq:
                            run_info.update(ena_fastq)
                        
                        runs.append(run_info)
                        self.save_sra_run(run_info)
            
            return {
                'studies': studies,
                'experiments': experiments,
                'runs': runs
            }
            
        except Exception as e:
            logging.error(f"批量获取SRA详情失败: {e}")
            return {'studies': [], 'experiments': [], 'runs': []}
    
    def fetch_ena_fastq_links(self, run_accession: str) -> Dict:
        """从ENA获取FASTQ下载链接"""
        try:
            # ENA Portal API
            url = f"{self.ena_base}/filereport"
            params = {
                'accession': run_accession,
                'result': 'read_run',
                'fields': 'run_accession,fastq_ftp,fastq_aspera,fastq_galaxy,fastq_md5',
                'format': 'json'
            }
            
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0:
                    record = data[0]
                    
                    fastq_info = {
                        'fastq_ftp': record.get('fastq_ftp', ''),
                        'fastq_aspera': record.get('fastq_aspera', ''),
                        'fastq_galaxy': record.get('fastq_galaxy', ''),
                        'fastq_md5': record.get('fastq_md5', ''),
                    }
                    
                    # 格式化FASTQ文件列表
                    fastq_files = []
                    if fastq_info['fastq_ftp']:
                        ftp_files = fastq_info['fastq_ftp'].split(';')
                        md5_list = fastq_info['fastq_md5'].split(';') if fastq_info['fastq_md5'] else []
                        
                        for i, ftp_file in enumerate(ftp_files):
                            file_info = {
                                'filename': ftp_file.split('/')[-1],
                                'ftp': f"ftp://{ftp_file}",
                                'http': f"http://{ftp_file}",
                                'md5': md5_list[i] if i < len(md5_list) else ''
                            }
                            fastq_files.append(file_info)
                    
                    fastq_info['ena_fastq_files'] = json.dumps(fastq_files, ensure_ascii=False)
                    
                    return fastq_info
            
            return {}
            
        except Exception as e:
            logging.debug(f"获取ENA FASTQ链接失败 {run_accession}: {e}")
            return {}
    
    def save_sra_study(self, study_info: Dict):
        """保存SRA Study"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO sra_studies
                (study_accession, bioproject_accession, title, abstract,
                 study_type, center_name, submission_date, publication_date,
                 last_update_date, raw_xml, collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                study_info['study_accession'],
                study_info['bioproject_accession'],
                study_info['title'],
                study_info['abstract'],
                study_info['study_type'],
                study_info['center_name'],
                study_info['submission_date'],
                study_info['publication_date'],
                study_info['last_update_date'],
                study_info['raw_xml'],
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存SRA Study失败: {e}")
    
    def save_sra_experiment(self, exp_info: Dict):
        """保存SRA Experiment"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO sra_experiments
                (experiment_accession, study_accession, biosample_accession,
                 title, library_name, library_strategy, library_source,
                 library_selection, library_layout, platform, instrument_model,
                 design_description, raw_xml, collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                exp_info['experiment_accession'],
                exp_info['study_accession'],
                exp_info['biosample_accession'],
                exp_info['title'],
                exp_info['library_name'],
                exp_info['library_strategy'],
                exp_info['library_source'],
                exp_info['library_selection'],
                exp_info['library_layout'],
                exp_info['platform'],
                exp_info['instrument_model'],
                exp_info['design_description'],
                exp_info['raw_xml'],
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存SRA Experiment失败: {e}")
    
    def save_sra_run(self, run_info: Dict):
        """保存SRA Run（增强版）"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO sra_runs
                (run_accession, experiment_accession, biosample_accession,
                 total_spots, total_bases, size_mb, published_date,
                 load_done_date, download_path, fastq_ftp, fastq_aspera,
                 fastq_galaxy, fastq_md5, ena_fastq_files, raw_xml, collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                run_info['run_accession'],
                run_info['experiment_accession'],
                run_info['biosample_accession'],
                run_info['total_spots'],
                run_info['total_bases'],
                run_info['size_mb'],
                run_info['published_date'],
                run_info['load_done_date'],
                run_info['download_path'],
                run_info.get('fastq_ftp', ''),
                run_info.get('fastq_aspera', ''),
                run_info.get('fastq_galaxy', ''),
                run_info.get('fastq_md5', ''),
                run_info.get('ena_fastq_files', ''),
                run_info['raw_xml'],
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存SRA Run失败: {e}")
    
    def fetch_pubmed_details(self, pmid_list: List[str]) -> Dict[str, Dict]:
        """批量获取PubMed文献详细信息"""
        if not pmid_list:
            return {}
        
        pubmed_details = {}
        
        # 去重
        unique_pmids = list(set(pmid_list))
        
        # 批量处理，每次最多200个
        batch_size = 200
        for i in range(0, len(unique_pmids), batch_size):
            batch = unique_pmids[i:i+batch_size]
            
            try:
                fetch_url = f"{self.ncbi_base}/efetch.fcgi"
                params = {
                    'db': 'pubmed',
                    'id': ','.join(batch),
                    'retmode': 'xml',
                    'email': self.email
                }
                
                response = requests.get(fetch_url, params=params, timeout=60)
                response.raise_for_status()
                root = ET.fromstring(response.content)
                
                for article in root.findall('.//PubmedArticle'):
                    pmid_elem = article.find('.//PMID')
                    if pmid_elem is None:
                        continue
                    
                    pmid = pmid_elem.text
                    
                    details = {
                        'pmid': pmid,
                        'title': '',
                        'abstract': '',
                        'authors': [],
                        'journal': '',
                        'publication_date': '',
                        'doi': '',
                        'pmc_id': '',
                        'citation_count': 0,
                        'mesh_terms': [],
                        'keywords': [],
                        'raw_xml': ET.tostring(article, encoding='unicode')
                    }
                    
                    # 标题
                    title_elem = article.find('.//ArticleTitle')
                    if title_elem is not None:
                        details['title'] = title_elem.text or ''
                    
                    # 摘要
                    abstract_texts = []
                    for abstract in article.findall('.//Abstract/AbstractText'):
                        label = abstract.get('Label', '')
                        text = abstract.text or ''
                        if label:
                            abstract_texts.append(f"{label}: {text}")
                        else:
                            abstract_texts.append(text)
                    details['abstract'] = ' '.join(abstract_texts)
                    
                    # 作者
                    for author in article.findall('.//Author'):
                        last_name = self.get_xml_text(author, './/LastName')
                        fore_name = self.get_xml_text(author, './/ForeName')
                        if last_name:
                            author_name = f"{fore_name} {last_name}".strip()
                            details['authors'].append(author_name)
                    
                    # 期刊
                    journal_elem = article.find('.//Journal/Title')
                    if journal_elem is not None:
                        details['journal'] = journal_elem.text or ''
                    
                    # 发表日期
                    pub_date = article.find('.//PubDate')
                    if pub_date is not None:
                        year = self.get_xml_text(pub_date, './/Year')
                        month = self.get_xml_text(pub_date, './/Month')
                        day = self.get_xml_text(pub_date, './/Day')
                        date_parts = [p for p in [year, month, day] if p]
                        details['publication_date'] = '-'.join(date_parts)
                    
                    # DOI
                    for article_id in article.findall('.//ArticleId'):
                        id_type = article_id.get('IdType', '')
                        if id_type == 'doi':
                            details['doi'] = article_id.text or ''
                        elif id_type == 'pmc':
                            details['pmc_id'] = article_id.text or ''
                    
                    # MeSH Terms
                    for mesh in article.findall('.//MeshHeading/DescriptorName'):
                        if mesh.text:
                            details['mesh_terms'].append(mesh.text)
                    
                    # Keywords
                    for keyword in article.findall('.//Keyword'):
                        if keyword.text:
                            details['keywords'].append(keyword.text)
                    
                    pubmed_details[pmid] = details
                    self.save_pubmed_article(details)
                
                time.sleep(0.5)
                
            except Exception as e:
                logging.error(f"获取PubMed详情失败: {e}")
                continue
        
        logging.info(f"  获取到 {len(pubmed_details)} 篇PubMed文献详情")
        return pubmed_details
    
    def save_pubmed_article(self, details: Dict):
        """保存PubMed文献"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO pubmed_articles
                (pmid, title, abstract, authors, journal, publication_date,
                 doi, pmc_id, citation_count, mesh_terms, keywords,
                 raw_xml, collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                details['pmid'],
                details['title'],
                details['abstract'],
                '; '.join(details['authors']),
                details['journal'],
                details['publication_date'],
                details['doi'],
                details['pmc_id'],
                details['citation_count'],
                '; '.join(details['mesh_terms']),
                '; '.join(details['keywords']),
                details['raw_xml'],
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存PubMed文献失败: {e}")
    
    def search_ena_studies(self) -> List[Dict]:
        """从ENA搜索人类单细胞研究"""
        logging.info("\n" + "="*70)
        logging.info("搜索ENA数据库")
        logging.info("="*70)
        
        all_studies = []
        
        # 多种搜索策略
        search_queries = [
            'tax_eq(9606) AND (library_strategy="RNA-Seq" OR library_strategy="scRNA-Seq")',
            'tax_eq(9606) AND description="single cell*"',
            'tax_eq(9606) AND description="single-cell*"',
            'tax_eq(9606) AND description="scRNA-seq*"',
            'tax_eq(9606) AND description="10x*"',
            'tax_eq(9606) AND description="single nucleus*"',
        ]
        
        for i, query in enumerate(search_queries, 1):
            logging.info(f"\n搜索策略 {i}/{len(search_queries)}")
            logging.info(f"查询: {query[:80]}...")
            
            try:
                url = f"{self.ena_base}/search"
                params = {
                    'result': 'study',
                    'query': query,
                    'fields': 'study_accession,secondary_study_accession,study_title,study_description,center_name,first_public,last_updated,study_type',
                    'format': 'json',
                    'limit': 0  # 获取所有结果
                }
                
                response = requests.get(url, params=params, timeout=60)
                
                if response.status_code == 200:
                    studies = response.json()
                    logging.info(f"  找到 {len(studies)} 个研究")
                    all_studies.extend(studies)
                else:
                    logging.warning(f"  搜索失败: {response.status_code}")
                
                time.sleep(1)
                
            except Exception as e:
                logging.error(f"  ENA搜索失败: {e}")
                continue
        
        # 去重
        unique_studies = {}
        for study in all_studies:
            acc = study.get('study_accession', '')
            if acc and acc not in unique_studies:
                unique_studies[acc] = study
        
        logging.info(f"\n总计找到 {len(unique_studies)} 个唯一ENA研究")
        
        # 保存到数据库
        for study in unique_studies.values():
            self.save_ena_study(study)
        
        return list(unique_studies.values())
    
    def save_ena_study(self, study: Dict):
        """保存ENA研究"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO ena_studies
                (study_accession, secondary_accession, title, description,
                 center_name, first_public, last_updated, study_type,
                 sample_count, experiment_count, run_count, tax_id,
                 scientific_name, fastq_files, collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                study.get('study_accession', ''),
                study.get('secondary_study_accession', ''),
                study.get('study_title', ''),
                study.get('study_description', ''),
                study.get('center_name', ''),
                study.get('first_public', ''),
                study.get('last_updated', ''),
                study.get('study_type', ''),
                0,  # sample_count - 需要额外查询
                0,  # experiment_count
                0,  # run_count
                9606,  # tax_id for Homo sapiens
                'Homo sapiens',
                '',
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存ENA研究失败: {e}")
    
    def extract_disease_info(self, text: str, biosample_attrs: Dict = None) -> tuple:
        """从文本和BioSample属性中提取疾病信息"""
        text_lower = text.lower() if text else ''
        
        detected_diseases = []
        specific_disease = ''
        
        # 从文本提取
        for disease, keywords in self.disease_keywords.items():
            if any(kw in text_lower for kw in keywords):
                detected_diseases.append(disease)
        
        # 从BioSample属性提取
        if biosample_attrs:
            for key, value in biosample_attrs.items():
                key_lower = key.lower()
                value_lower = str(value).lower()
                
                if any(x in key_lower for x in ['disease', 'condition', 'phenotype', 'diagnosis']):
                    if value_lower not in ['normal', 'healthy', 'control', 'none', 'na', 'n/a']:
                        specific_disease = value
                        
                        # 根据具体疾病添加到大类
                        for disease, keywords in self.disease_keywords.items():
                            if any(kw in value_lower for kw in keywords):
                                if disease not in detected_diseases:
                                    detected_diseases.append(disease)
        
        # 从文本中提取具体疾病名称
        if not specific_disease:
            disease_patterns = [
                r'([\w\s-]+cancer)',
                r'([\w\s-]+carcinoma)',
                r'(covid[-\s]?19)',
                r'(type\s+[12]\s+diabetes)',
                r'(alzheimer\'?s?\s+disease)',
                r'([\w\s-]+leukemia)',
                r'([\w\s-]+lymphoma)',
            ]
            
            for pattern in disease_patterns:
                match = re.search(pattern, text_lower)
                if match:
                    specific_disease = match.group(1).strip()
                    break
        
        disease_general = ';'.join(detected_diseases) if detected_diseases else 'unknown'
        return disease_general, specific_disease
    
    def extract_tissue_info(self, text: str, biosample_attrs: Dict = None) -> str:
        """提取组织信息"""
        tissues = []
        
        # 从BioSample属性提取
        if biosample_attrs:
            tissue_fields = ['tissue', 'cell type', 'source name', 'source_name',
                           'cell_type', 'tissue_type', 'organ', 'anatomical site']
            
            for key, value in biosample_attrs.items():
                key_lower = key.lower()
                if any(field in key_lower for field in tissue_fields):
                    if value and str(value).lower() not in ['na', 'n/a', 'none', 'unknown']:
                        tissues.append(str(value))
        
        # 从文本中提取
        if text:
            text_lower = text.lower()
            for tissue_kw in self.tissue_keywords:
                if tissue_kw in text_lower:
                    tissues.append(tissue_kw)
        
        # 去重并返回
        unique_tissues = []
        for t in tissues:
            if t not in unique_tissues:
                unique_tissues.append(t)
        
        return '; '.join(unique_tissues[:3]) if unique_tissues else ''
    
    def extract_sex_info(self, biosample_attrs: Dict = None) -> str:
        """提取性别信息"""
        if not biosample_attrs:
            return ''
        
        sex_fields = ['sex', 'gender']
        
        for key, value in biosample_attrs.items():
            key_lower = key.lower()
            if any(field in key_lower for field in sex_fields):
                if value:
                    value_lower = str(value).lower()
                    if 'male' in value_lower and 'female' not in value_lower:
                        return 'male'
                    elif 'female' in value_lower:
                        return 'female'
                    elif any(x in value_lower for x in ['mixed', 'both', 'pooled']):
                        return 'mixed'
        
        return ''
    
    def extract_ethnicity_info(self, biosample_attrs: Dict = None) -> str:
        """提取种族信息"""
        if not biosample_attrs:
            return ''
        
        ethnicity_fields = ['race', 'ethnicity', 'ancestry', 'population']
        
        for key, value in biosample_attrs.items():
            key_lower = key.lower()
            if any(field in key_lower for field in ethnicity_fields):
                if value and str(value).lower() not in ['na', 'n/a', 'none', 'unknown', 'not provided']:
                    return str(value)
        
        return ''
    
    def format_final_record(self, bp_details: Dict, biosamples: List[Dict], 
                           sra_data: Dict, pubmed_data: Dict = None) -> Dict:
        """格式化最终记录（增强版）"""
        # 合并所有样本属性
        all_attributes = {}
        for bs in biosamples:
            all_attributes.update(bs.get('attributes', {}))
        
        # 提取信息
        combined_text = f"{bp_details.get('title', '')} {bp_details.get('description', '')}"
        
        # 疾病信息
        disease_general, specific_disease = self.extract_disease_info(combined_text, all_attributes)
        
        # 组织信息
        tissue_from_samples = '; '.join([bs.get('tissue', '') for bs in biosamples if bs.get('tissue')])
        tissue = self.extract_tissue_info(combined_text, all_attributes)
        if not tissue and tissue_from_samples:
            tissue = tissue_from_samples
        
        # 性别、种族
        sex = self.extract_sex_info(all_attributes)
        ethnicity = self.extract_ethnicity_info(all_attributes)
        
        # 测序平台
        platforms = set()
        for exp in sra_data.get('experiments', []):
            platform = exp.get('instrument_model', '')
            if platform:
                platforms.add(platform)
        
        # 样本类型
        sample_type = 'single-cell'
        library_strategies = set()
        for exp in sra_data.get('experiments', []):
            strategy = exp.get('library_strategy', '').lower()
            library_strategies.add(strategy)
            
            if 'single nucleus' in combined_text.lower() or 'sn-rna' in strategy:
                sample_type = 'single-nucleus'
            elif 'spatial' in strategy or 'spatial' in combined_text.lower():
                sample_type = 'spatial'
            elif 'cite-seq' in combined_text.lower():
                sample_type = 'multiome'
        
        # 数据量统计
        runs = sra_data.get('runs', [])
        total_bases = sum(run.get('total_bases', 0) for run in runs)
        total_spots = sum(run.get('total_spots', 0) for run in runs)
        
        # 发表日期（从runs中获取最早的）
        publication_dates = [run.get('published_date', '') for run in runs if run.get('published_date')]
        publication_date = min(publication_dates) if publication_dates else ''
        
        # PubMed信息
        pmid_list = bp_details.get('publications', [])
        pubmed_titles = []
        pubmed_abstracts = []
        
        if pubmed_data:
            for pmid in pmid_list:
                if pmid in pubmed_data:
                    pm = pubmed_data[pmid]
                    if pm.get('title'):
                        pubmed_titles.append(f"[PMID:{pmid}] {pm['title']}")
                    if pm.get('abstract'):
                        pubmed_abstracts.append(f"[PMID:{pmid}] {pm['abstract'][:500]}...")
        
        # FASTQ下载链接
        fastq_links = []
        ena_fastq_all = []
        
        for run in runs:
            run_acc = run.get('run_accession', '')
            
            # ENA FASTQ files
            ena_files_str = run.get('ena_fastq_files', '')
            if ena_files_str:
                try:
                    ena_files = json.loads(ena_files_str)
                    for file_info in ena_files:
                        fastq_links.append(f"{run_acc}: {file_info['http']}")
                        ena_fastq_all.append(file_info)
                except:
                    pass
            
            # FTP links
            if run.get('fastq_ftp'):
                for ftp in run['fastq_ftp'].split(';'):
                    fastq_links.append(f"{run_acc}: ftp://{ftp}")
        
        record = {
            'id': bp_details['accession'],
            'source_database': 'BioProject/SRA + ENA',
            'title': bp_details['title'],
            'disease_general': disease_general,
            'disease': specific_disease,
            'pubmed': ','.join(pmid_list),
            'pubmed_titles': ' | '.join(pubmed_titles),
            'pubmed_abstracts': ' | '.join(pubmed_abstracts),
            'access_link': f"https://www.ncbi.nlm.nih.gov/bioproject/{bp_details['accession']}",
            'open_status': 'public',
            'ethnicity': ethnicity,
            'sex': sex,
            'tissue': tissue,
            'sequencing_platform': '; '.join(platforms) if platforms else '',
            'experiment_design': sample_type,
            'sample_type': sample_type,
            'summary': bp_details.get('description', ''),
            'citation_count': len(pmid_list),
            'publication_date': publication_date,
            'submission_date': bp_details.get('submission_date', ''),
            'last_update_date': bp_details.get('last_update_date', ''),
            'contact_name': bp_details.get('submitter_name', ''),
            'contact_email': bp_details.get('submitter_email', ''),
            'contact_institute': bp_details.get('submitter_organization', ''),
            'data_tier': 'raw',
            'tissue_location': '',
            'supplementary_information': f"Library strategies: {'; '.join(library_strategies)}",
            'bioproject': bp_details['accession'],
            'biosample_list': '; '.join([bs['accession'] for bs in biosamples if bs.get('accession')]),
            'run_count': len(runs),
            'total_bases': total_bases,
            'total_spots': total_spots,
            'fastq_download_links': '\n'.join(fastq_links[:10]),  # 限制前10个链接
            'ena_accession': '',
            'ena_fastq_links': json.dumps(ena_fastq_all[:20], ensure_ascii=False),  # 限制前20个文件
            'metadata_completeness': 0.0
        }
        
        # 计算元数据完整度
        important_fields = [
            'title', 'summary', 'disease', 'tissue', 'sex', 'ethnicity',
            'sequencing_platform', 'sample_type', 'biosample_list', 
            'contact_name', 'run_count', 'pubmed', 'fastq_download_links'
        ]
        
        filled = sum(1 for field in important_fields
                    if record.get(field) and str(record[field]).strip() 
                    and str(record[field]) not in ['0', 'unknown'])
        
        record['metadata_completeness'] = round(filled / len(important_fields), 2)
        
        return record
    
    def save_final_record(self, record: Dict):
        """保存最终记录到数据库"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO final_series
                (id, source_database, title, disease_general, disease, pubmed,
                 pubmed_titles, pubmed_abstracts, access_link, open_status,
                 ethnicity, sex, tissue, sequencing_platform, experiment_design,
                 sample_type, summary, citation_count, publication_date,
                 submission_date, last_update_date, contact_name, contact_email,
                 contact_institute, data_tier, tissue_location,
                 supplementary_information, bioproject, biosample_list,
                 run_count, total_bases, total_spots, fastq_download_links,
                 ena_accession, ena_fastq_links, metadata_completeness,
                 collection_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                record['id'], record['source_database'], record['title'],
                record['disease_general'], record['disease'], record['pubmed'],
                record['pubmed_titles'], record['pubmed_abstracts'],
                record['access_link'], record['open_status'], record['ethnicity'],
                record['sex'], record['tissue'], record['sequencing_platform'],
                record['experiment_design'], record['sample_type'], record['summary'],
                record['citation_count'], record['publication_date'],
                record['submission_date'], record['last_update_date'],
                record['contact_name'], record['contact_email'], record['contact_institute'],
                record['data_tier'], record['tissue_location'],
                record['supplementary_information'], record['bioproject'],
                record['biosample_list'], record['run_count'], record['total_bases'],
                record['total_spots'], record['fastq_download_links'],
                record['ena_accession'], record['ena_fastq_links'],
                record['metadata_completeness'], datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logging.error(f"保存最终记录失败: {e}")
    
    def get_xml_text(self, element, path: str) -> str:
        """安全获取XML文本"""
        try:
            found = element.find(path)
            return found.text if found is not None and found.text else ''
        except:
            return ''
    
    def collect_all_data(self, max_projects: Optional[int] = None, 
                        batch_size: int = 10, include_ena: bool = True):
        """收集所有BioProject/BioSample/SRA数据（增强版）"""
        logging.info("\n" + "="*70)
        logging.info("开始收集BioProject/BioSample/SRA + ENA数据")
        logging.info("="*70)
        
        # 1. 搜索所有BioProject
        bp_ids = self.search_bioprojects()
        
        # 2. 可选：搜索ENA
        if include_ena:
            ena_studies = self.search_ena_studies()
        
        if max_projects:
            bp_ids = bp_ids[:max_projects]
            logging.info(f"限制处理 {max_projects} 个Projects")
        
        all_records = []
        total = len(bp_ids)
        
        # 3. 处理每个BioProject
        for i, bp_id in enumerate(bp_ids, 1):
            logging.info(f"\n{'='*70}")
            logging.info(f"处理进度: {i}/{total} ({i/total*100:.1f}%)")
            logging.info(f"BioProject ID: {bp_id}")
            logging.info(f"{'='*70}")
            
            try:
                # 3.1 获取BioProject详情
                logging.info("  [1/5] 获取BioProject详情...")
                bp_details = self.fetch_bioproject_details(bp_id)
                
                if not bp_details or not bp_details.get('accession'):
                    logging.warning("  ✗ 无法获取BioProject详情，跳过")
                    continue
                
                logging.info(f"  ✓ {bp_details['accession']}: {bp_details['title'][:60]}...")
                
                # 3.2 获取关联的BioSamples
                logging.info("  [2/5] 获取关联的BioSamples...")
                biosample_ids = self.fetch_linked_biosamples(bp_details['accession'])
                logging.info(f"  找到 {len(biosample_ids)} 个BioSamples")
                
                biosamples = []
                for j, bs_id in enumerate(biosample_ids[:100], 1):
                    if j % 10 == 0:
                        logging.info(f"    处理BioSample: {j}/{min(len(biosample_ids), 100)}")
                    
                    bs_details = self.fetch_biosample_details(bs_id)
                    if bs_details:
                        biosamples.append(bs_details)
                    time.sleep(0.2)
                
                logging.info(f"  ✓ 成功获取 {len(biosamples)} 个BioSamples")
                
                # 3.3 获取关联的SRA数据
                logging.info("  [3/5] 获取关联的SRA数据...")
                sra_ids = self.fetch_linked_sra_studies(bp_details['accession'])
                logging.info(f"  找到 {len(sra_ids)} 个SRA记录")
                
                sra_data = {'studies': [], 'experiments': [], 'runs': []}
                
                if sra_ids:
                    for j in range(0, len(sra_ids), 100):
                        batch = sra_ids[j:j+100]
                        logging.info(f"    处理SRA批次: {j//100 + 1}/{(len(sra_ids)-1)//100 + 1}")
                        
                        batch_data = self.fetch_sra_details_batch(batch)
                        
                        sra_data['studies'].extend(batch_data['studies'])
                        sra_data['experiments'].extend(batch_data['experiments'])
                        sra_data['runs'].extend(batch_data['runs'])
                        
                        time.sleep(1)
                
                logging.info(f"  ✓ Studies: {len(sra_data['studies'])}, "
                           f"Experiments: {len(sra_data['experiments'])}, "
                           f"Runs: {len(sra_data['runs'])}")
                
                # 3.4 获取PubMed文献详情
                logging.info("  [4/5] 获取PubMed文献详情...")
                pubmed_data = {}
                pmid_list = bp_details.get('publications', [])
                
                if pmid_list:
                    pubmed_data = self.fetch_pubmed_details(pmid_list)
                    logging.info(f"  ✓ 获取到 {len(pubmed_data)} 篇文献")
                else:
                    logging.info("  - 无关联文献")
                
                # 3.5 格式化并保存最终记录
                logging.info("  [5/5] 格式化并保存记录...")
                final_record = self.format_final_record(bp_details, biosamples, 
                                                       sra_data, pubmed_data)
                
                self.save_final_record(final_record)
                all_records.append(final_record)
                
                logging.info(f"  ✓ 完成! 元数据完整度: {final_record['metadata_completeness']:.0%}")
                
                # 进度汇总
                if i % 10 == 0:
                    logging.info(f"\n{'='*70}")
                    logging.info(f"阶段性统计 ({i}/{total}):")
                    logging.info(f"  已处理: {len(all_records)} 个项目")
                    avg_completeness = sum(r['metadata_completeness'] for r in all_records) / len(all_records)
                    logging.info(f"  平均完整度: {avg_completeness:.1%}")
                    logging.info(f"{'='*70}\n")
                
                # 限速
                if i % 5 == 0:
                    time.sleep(2)
                else:
                    time.sleep(0.5)
                
            except KeyboardInterrupt:
                logging.warning("\n收集被用户中断")
                break
                
            except Exception as e:
                logging.error(f"  ✗ 处理失败: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        logging.info(f"\n{'='*70}")
        logging.info("数据收集完成!")
        logging.info(f"总计收集: {len(all_records)} 个项目")
        logging.info(f"{'='*70}")
        
        return all_records
    
    def export_final_data(self):
        """导出最终数据"""
        logging.info("\n" + "="*70)
        logging.info("导出最终数据")
        logging.info("="*70)
        
        conn = sqlite3.connect(self.db_path)
        
        # 导出最终整合数据
        final_df = pd.read_sql_query("SELECT * FROM final_series", conn)
        
        if len(final_df) > 0:
            # 去重
            final_df = final_df.sort_values('metadata_completeness', ascending=False)
            final_df = final_df.drop_duplicates(subset=['id'], keep='first')
            
            # 保存CSV
            csv_output = self.output_dir / 'bioproject_sra_metadata_enhanced.csv'
            final_df.to_csv(csv_output, index=False, encoding='utf-8-sig')
            logging.info(f"✓ CSV导出: {csv_output}")
            
            # 保存Excel
            excel_output = self.output_dir / 'bioproject_sra_metadata_enhanced.xlsx'
            final_df.to_excel(excel_output, index=False, engine='openpyxl')
            logging.info(f"✓ Excel导出: {excel_output}")
            
            # 生成统计报告
            self.generate_statistics_report(final_df)
            
            # 导出原始数据表
            self.export_raw_tables(conn)
            
            # 导出FASTQ下载链接单独文件
            self.export_fastq_links(final_df)
            
        else:
            logging.warning("没有数据可导出")
        
        conn.close()
        
        return final_df
    
    def export_fastq_links(self, df: pd.DataFrame):
        """导出FASTQ下载链接到单独文件"""
        try:
            fastq_data = []
            
            for _, row in df.iterrows():
                bioproject = row.get('bioproject', '')
                links_str = row.get('fastq_download_links', '')
                
                if links_str:
                    for link in links_str.split('\n'):
                        if link.strip():
                            parts = link.split(': ', 1)
                            if len(parts) == 2:
                                run_acc, url = parts
                                fastq_data.append({
                                    'BioProject': bioproject,
                                    'Run_Accession': run_acc,
                                    'FASTQ_URL': url,
                                    'Type': 'HTTP' if url.startswith('http') else 'FTP'
                                })
            
            if fastq_data:
                fastq_df = pd.DataFrame(fastq_data)
                output_path = self.output_dir / 'fastq_download_links.csv'
                fastq_df.to_csv(output_path, index=False, encoding='utf-8-sig')
                logging.info(f"✓ FASTQ链接导出: {output_path} ({len(fastq_data)} 个文件)")
                
        except Exception as e:
            logging.error(f"导出FASTQ链接失败: {e}")
    
    def generate_statistics_report(self, df: pd.DataFrame):
        """生成统计报告（增强版）"""
        report_path = self.output_dir / 'collection_statistics_enhanced.txt'
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("NCBI BioProject/SRA + ENA 人源单细胞数据收集统计报告\n")
            f.write("="*70 + "\n\n")
            f.write(f"收集时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总记录数: {len(df)}\n\n")
            
            # 元数据完整度分析
            if 'metadata_completeness' in df.columns:
                avg_completeness = df['metadata_completeness'].mean()
                f.write(f"平均元数据完整度: {avg_completeness:.2%}\n")
                f.write(f"完整度 >= 80%: {(df['metadata_completeness'] >= 0.8).sum()} 条\n")
                f.write(f"完整度 >= 60%: {(df['metadata_completeness'] >= 0.6).sum()} 条\n")
                f.write(f"完整度 < 40%: {(df['metadata_completeness'] < 0.4).sum()} 条\n\n")
            
            # FASTQ可用性
            if 'fastq_download_links' in df.columns:
                has_fastq = df['fastq_download_links'].notna() & (df['fastq_download_links'] != '')
                f.write(f"FASTQ链接可用性:\n")
                f.write(f"  有FASTQ链接: {has_fastq.sum()} ({has_fastq.sum()/len(df)*100:.1f}%)\n")
                f.write(f"  无FASTQ链接: {(~has_fastq).sum()} ({(~has_fastq).sum()/len(df)*100:.1f}%)\n\n")
            
            # PubMed覆盖率
            if 'pubmed' in df.columns:
                has_pubmed = df['pubmed'].notna() & (df['pubmed'] != '')
                f.write(f"PubMed文献覆盖率:\n")
                f.write(f"  有文献: {has_pubmed.sum()} ({has_pubmed.sum()/len(df)*100:.1f}%)\n")
                f.write(f"  无文献: {(~has_pubmed).sum()} ({(~has_pubmed).sum()/len(df)*100:.1f}%)\n\n")
            
            # 疾病分布
            if 'disease_general' in df.columns:
                f.write("疾病类型分布 (Top 20):\n")
                diseases = df['disease_general'].str.split(';').explode()
                disease_counts = diseases.value_counts()
                for disease, count in disease_counts.head(20).items():
                    if disease and disease != 'unknown':
                        f.write(f"  {disease}: {count} ({count/len(df)*100:.1f}%)\n")
                f.write("\n")
            
            # 组织分布
            if 'tissue' in df.columns:
                f.write("组织类型分布 (Top 20):\n")
                tissues = df['tissue'].str.split(';').explode()
                tissue_counts = tissues.value_counts()
                for tissue, count in tissue_counts.head(20).items():
                    if tissue and tissue.strip():
                        f.write(f"  {tissue.strip()}: {count}\n")
                f.write("\n")
            
            # 测序平台
            if 'sequencing_platform' in df.columns:
                f.write("测序平台分布 (Top 15):\n")
                platforms = df['sequencing_platform'].str.split(';').explode()
                for platform, count in platforms.value_counts().head(15).items():
                    if platform and platform.strip():
                        f.write(f"  {platform.strip()}: {count}\n")
                f.write("\n")
            
            # 性别分布
            if 'sex' in df.columns:
                f.write("性别分布:\n")
                for sex, count in df['sex'].value_counts().items():
                    if sex:
                        f.write(f"  {sex}: {count} ({count/len(df)*100:.1f}%)\n")
                f.write("\n")
            
            # 样本类型
            if 'sample_type' in df.columns:
                f.write("样本类型分布:\n")
                for stype, count in df['sample_type'].value_counts().items():
                    if stype:
                        f.write(f"  {stype}: {count} ({count/len(df)*100:.1f}%)\n")
                f.write("\n")
            
            # 数据量统计
            if 'total_bases' in df.columns:
                total_tb = df['total_bases'].sum() / (1024**4)
                avg_gb = df['total_bases'].mean() / (1024**3)
                f.write(f"数据量统计:\n")
                f.write(f"  总数据量: {total_tb:.2f} TB\n")
                f.write(f"  平均每个项目: {avg_gb:.2f} GB\n\n")
            
            # Run数量统计
            if 'run_count' in df.columns:
                total_runs = df['run_count'].sum()
                avg_runs = df['run_count'].mean()
                f.write(f"SRA Run统计:\n")
                f.write(f"  总Run数: {total_runs}\n")
                f.write(f"  平均每个项目: {avg_runs:.1f} runs\n\n")
            
            # 时间分布
            if 'submission_date' in df.columns:
                f.write("提交时间分布 (按年):\n")
                df_temp = df[df['submission_date'].notna()].copy()
                if len(df_temp) > 0:
                    df_temp['year'] = pd.to_datetime(df_temp['submission_date'], 
                                                     errors='coerce').dt.year
                    for year, count in df_temp['year'].value_counts().sort_index().items():
                        if pd.notna(year):
                            f.write(f"  {int(year)}: {count}\n")
        
        logging.info(f"✓ 统计报告: {report_path}")
    
    def export_raw_tables(self, conn):
        """导出原始数据表（增强版）"""
        tables = {
            'bioprojects': 'raw_bioprojects.csv',
            'biosamples': 'raw_biosamples.csv',
            'sra_studies': 'raw_sra_studies.csv',
            'sra_experiments': 'raw_sra_experiments.csv',
            'sra_runs': 'raw_sra_runs.csv',
            'pubmed_articles': 'raw_pubmed_articles.csv',
            'ena_studies': 'raw_ena_studies.csv'
        }
        
        for table_name, filename in tables.items():
            try:
                df = pd.read_sql_query(f"SELECT * FROM {table_name}", conn)
                if len(df) > 0:
                    output_path = self.output_dir / filename
                    df.to_csv(output_path, index=False, encoding='utf-8-sig')
                    logging.info(f"✓ {table_name}: {len(df)} 条记录 -> {output_path}")
            except Exception as e:
                logging.debug(f"导出 {table_name} 失败: {e}")
    
    def run_full_collection(self, max_projects: Optional[int] = None,
                           include_ena: bool = True):
        """运行完整收集流程（增强版）"""
        logging.info("="*70)
        logging.info("NCBI BioProject/BioSample/SRA + ENA 人源单细胞数据收集")
        logging.info("增强版：包含FASTQ链接、PubMed详情、ENA检索")
        logging.info("="*70)
        
        start_time = datetime.now()
        
        # 1. 收集数据
        records = self.collect_all_data(max_projects=max_projects, 
                                       include_ena=include_ena)
        
        # 2. 导出数据
        final_df = self.export_final_data()
        
        end_time = datetime.now()
        duration = end_time - start_time
        
        logging.info("\n" + "="*70)
        logging.info("收集完成!")
        logging.info("="*70)
        logging.info(f"总耗时: {duration}")
        logging.info(f"收集记录: {len(records)}")
        logging.info(f"最终记录: {len(final_df) if final_df is not None else 0}")
        logging.info(f"\n所有文件保存在: {self.output_dir}")
        
        return final_df


def main():
    """主函数"""
    print("""
    ╔════════════════════════════════════════════════════════════════════╗
    ║  NCBI BioProject/BioSample/SRA + ENA 人源单细胞数据收集系统          ║
    ║                        增强版 v2.0                                 ║
    ║                                                                    ║
    ║  新增功能:                                                          ║
    ║  ✓ FASTQ文件直接下载链接 (ENA FTP/HTTP)                            ║
    ║  ✓ PubMed文献完整信息 (标题、摘要、作者)                            ║
    ║  ✓ ENA数据库独立检索                                               ║
    ║  ✓ 增强的元数据完整度评分                                          ║
    ║                                                                    ║
    ║  收集范围:                                                          ║
    ║  • BioProject - 项目信息                                           ║
    ║  • BioSample - 样本元数据                                          ║
    ║  • SRA Study/Experiment/Run - 测序数据详情                         ║
    ║  • PubMed - 文献详细信息                                           ║
    ║  • ENA - 欧洲核酸数据库                                            ║
    ║                                                                    ║
    ║  输出内容:                                                          ║
    ║  • 标准化的元数据表格 (CSV + Excel)                                ║
    ║  • FASTQ下载链接清单                                               ║
    ║  • PubMed文献详情                                                  ║
    ║  • 完整的统计分析报告                                               ║
    ║  • SQLite数据库                                                    ║
    ╚════════════════════════════════════════════════════════════════════╝
    """)
    
    # 配置
    
    email = 'fenghe13254@gmail'
    
    

    ena_choice = '2'
    include_ena = ena_choice != '2'
    
    
    
    max_projects = None
    
    
    
    # 初始化收集器
    collector = NCBIBioProjectSRACollector(
        output_dir='ncbi_bioproject_sra_data',
        email=email
    )
    
    # 执行收集
    try:
        print("\n" + "="*70)
        print("开始收集数据...")
        print("="*70)
        
        final_df = collector.run_full_collection(
            max_projects=max_projects,
            include_ena=include_ena
        )
        
        if final_df is not None and len(final_df) > 0:
            print("\n" + "="*70)
            print("✓ 数据收集成功!")
            print("="*70)
            print(f"\n总记录数: {len(final_df)}")
            print(f"数据目录: {collector.output_dir}")
            
            print("\n主要输出文件:")
            print("  1. bioproject_sra_metadata_enhanced.csv - 最终数据集(CSV)")
            print("  2. bioproject_sra_metadata_enhanced.xlsx - 最终数据集(Excel)")
            print("  3. fastq_download_links.csv - FASTQ下载链接")
            print("  4. collection_statistics_enhanced.txt - 详细统计报告")
            print("  5. bioproject_sra_metadata.db - SQLite数据库")
            print("  6. raw_*.csv - 各表原始数据")
            
            print("\n数据预览:")
            preview_cols = ['id', 'title', 'disease_general', 'tissue', 
                          'sample_type', 'run_count', 'pubmed']
            available_cols = [col for col in preview_cols if col in final_df.columns]
            print(final_df[available_cols].head(10))
            
            print("\n元数据完整度分布:")
            completeness_dist = final_df['metadata_completeness'].describe()
            print(completeness_dist)
            
            # FASTQ可用性
            if 'fastq_download_links' in final_df.columns:
                has_fastq = final_df['fastq_download_links'].notna() & \
                           (final_df['fastq_download_links'] != '')
                print(f"\nFASTQ链接可用性: {has_fastq.sum()}/{len(final_df)} " +
                      f"({has_fastq.sum()/len(final_df)*100:.1f}%)")
            
            # PubMed覆盖率
            if 'pubmed' in final_df.columns:
                has_pubmed = final_df['pubmed'].notna() & (final_df['pubmed'] != '')
                print(f"PubMed文献覆盖率: {has_pubmed.sum()}/{len(final_df)} " +
                      f"({has_pubmed.sum()/len(final_df)*100:.1f}%)")
            
            print("\n疾病类型统计:")
            diseases = final_df['disease_general'].str.split(';').explode()
            print(diseases.value_counts().head(10))
            
        else:
            print("\n⚠ 未收集到数据")
            print("请查看日志: ncbi_bioproject_sra_collection.log")
    
    except KeyboardInterrupt:
        print("\n\n收集被用户中断")
        print("已收集的数据已保存到数据库")
        print(f"可以查看: {collector.output_dir}")
    
    except Exception as e:
        print(f"\n❌ 收集失败: {e}")
        import traceback
        traceback.print_exc()
        print("\n请查看日志文件获取详细错误信息")


if __name__ == "__main__":
    main()