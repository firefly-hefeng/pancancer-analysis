import requests
import json
import pandas as pd
import time
import os
import logging
from tqdm import tqdm
from collections import defaultdict
from urllib.parse import quote
import re
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ===== 配置部分 =====
BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1"
BIOSTUDIES_SEARCH = f"{BIOSTUDIES_API}/search"
BIOSTUDIES_STUDY = f"{BIOSTUDIES_API}/studies"

SCEA_API = "https://www.ebi.ac.uk/gxa/sc/json/experiments"
AE_FILES_BASE = "https://www.ebi.ac.uk/biostudies/files"
EUROPEPMC_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
ENA_API_BASE = "https://www.ebi.ac.uk/ena/portal/api"

# 输出文件
RAW_BS_JSON = "biostudies_human_scrna_raw.json"
RAW_SCEA_JSON = "scea_human_scrna_raw.json"

EXPERIMENT_METADATA_CSV = "1_experiment_metadata.csv"
SAMPLE_METADATA_CSV = "2_sample_metadata.csv"
FILE_MANIFEST_CSV = "3_file_manifest.csv"
SAMPLE_FILE_MAPPING_CSV = "4_sample_file_mapping.csv"
SEQUENCING_RUNS_CSV = "5_sequencing_runs.csv"
PUBLICATIONS_CSV = "6_publications.csv"
STATISTICS_REPORT = "7_collection_statistics.txt"


# ===== 工具函数：标准化和清洗 =====
def clean_value(value):
    """清洗单个值"""
    if pd.isna(value) or value == '' or value == 'not applicable' or value == 'N/A':
        return None
    return str(value).strip()


def standardize_sex(sex_value):
    """标准化性别值"""
    if not sex_value:
        return None
    sex_lower = str(sex_value).lower()
    if 'female' in sex_lower or sex_lower == 'f':
        return 'Female'
    elif 'male' in sex_lower or sex_lower == 'm':
        return 'Male'
    elif 'unknown' in sex_lower or 'not' in sex_lower:
        return 'Unknown'
    else:
        return sex_value


def parse_age(age_value):
    """解析年龄，提取数值和单位"""
    if not age_value:
        return None, None
    
    age_str = str(age_value).lower()
    numbers = re.findall(r'\d+\.?\d*', age_str)
    if not numbers:
        return age_value, None
    
    age_num = float(numbers[0])
    
    if 'year' in age_str or 'yr' in age_str:
        return age_num, 'years'
    elif 'month' in age_str or 'mo' in age_str:
        return age_num, 'months'
    elif 'week' in age_str or 'wk' in age_str:
        return age_num, 'weeks'
    elif 'day' in age_str:
        return age_num, 'days'
    elif 'gestation' in age_str or 'pcw' in age_str:
        return age_num, 'gestational_weeks'
    else:
        return age_num, 'unknown_unit'


def standardize_disease(disease_value):
    """标准化疾病名称"""
    if not disease_value:
        return 'Normal', 'Healthy'
    
    disease_lower = str(disease_value).lower()
    
    if any(term in disease_lower for term in ['normal', 'healthy', 'control', 'non-disease']):
        return 'Normal', 'Healthy'
    
    cancer_keywords = {
        'carcinoma': 'Cancer', 'cancer': 'Cancer', 'tumor': 'Cancer',
        'tumour': 'Cancer', 'melanoma': 'Cancer', 'leukemia': 'Cancer',
        'lymphoma': 'Cancer', 'sarcoma': 'Cancer', 'glioblastoma': 'Cancer',
        'neuroblastoma': 'Cancer'
    }
    
    for keyword, category in cancer_keywords.items():
        if keyword in disease_lower:
            return category, disease_value
    
    if 'diabetes' in disease_lower:
        return 'Metabolic Disease', disease_value
    elif 'covid' in disease_lower or 'sars-cov-2' in disease_lower:
        return 'Infectious Disease', disease_value
    elif any(term in disease_lower for term in ['alzheimer', 'parkinson', 'neurolog']):
        return 'Neurological Disease', disease_value
    elif any(term in disease_lower for term in ['arthritis', 'lupus', 'autoimmune']):
        return 'Autoimmune Disease', disease_value
    elif 'fibrosis' in disease_lower:
        return 'Fibrotic Disease', disease_value
    else:
        return 'Other Disease', disease_value


def standardize_sequencing_platform(platform_value):
    """标准化测序平台名称"""
    if not platform_value:
        return None, None
    
    platform_lower = str(platform_value).lower()
    
    illumina_models = {
        'nextseq 500': ('Illumina', 'NextSeq 500'),
        'nextseq 550': ('Illumina', 'NextSeq 550'),
        'nextseq 2000': ('Illumina', 'NextSeq 2000'),
        'hiseq 2000': ('Illumina', 'HiSeq 2000'),
        'hiseq 2500': ('Illumina', 'HiSeq 2500'),
        'hiseq 3000': ('Illumina', 'HiSeq 3000'),
        'hiseq 4000': ('Illumina', 'HiSeq 4000'),
        'hiseq x': ('Illumina', 'HiSeq X'),
        'novaseq 6000': ('Illumina', 'NovaSeq 6000'),
        'novaseq': ('Illumina', 'NovaSeq'),
        'miseq': ('Illumina', 'MiSeq'),
        'iseq': ('Illumina', 'iSeq'),
        'genome analyzer': ('Illumina', 'Genome Analyzer'),
    }
    
    for keyword, (vendor, model) in illumina_models.items():
        if keyword in platform_lower:
            return vendor, model
    
    if 'pacbio' in platform_lower or 'pacific biosciences' in platform_lower:
        return 'PacBio', platform_value
    elif 'oxford nanopore' in platform_lower or 'minion' in platform_lower or 'promethion' in platform_lower:
        return 'Oxford Nanopore', platform_value
    elif 'bgiseq' in platform_lower or 'dnbseq' in platform_lower:
        return 'BGI', platform_value
    elif 'ion torrent' in platform_lower:
        return 'Ion Torrent', platform_value
    elif '454' in platform_lower:
        return 'Roche', '454'
    elif 'solid' in platform_lower:
        return 'ABI SOLiD', platform_value
    elif 'illumina' in platform_lower:
        return 'Illumina', platform_value
    else:
        return 'Unknown', platform_value


# ===== 第一步：使用BioStudies API收集数据 =====
def collect_biostudies_metadata():
    """从BioStudies API收集单细胞数据"""
    if os.path.exists(RAW_BS_JSON):
        logging.info(f"'{RAW_BS_JSON}' 已存在，跳过收集。")
        with open(RAW_BS_JSON, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    logging.info("开始从BioStudies API收集元数据...")
    all_studies = []
    all_accessions = set()
    
    # *** 使用经过验证的搜索策略 ***
    search_strategies = [
        # 策略1: 单细胞关键词 + 人类 (最精准)
        {
            'query': 'organism:"Homo sapiens" AND (title:"single cell" OR description:"single cell")',
            'pageSize': 100,
            'name': '单细胞关键词搜索'
        },
        # 策略2: RNA-seq + 人类 + ArrayExpress (获取所有ArrayExpress的RNA-seq)
        {
            'query': 'organism:"Homo sapiens" AND study_type:"RNA-seq" AND accession:E-*',
            'pageSize': 100,
            'name': 'ArrayExpress RNA-seq'
        },
        # 策略3: 单细胞技术关键词
        {
            'query': 'organism:"Homo sapiens" AND (title:"scRNA-seq" OR description:"scRNA-seq")',
            'pageSize': 100,
            'name': 'scRNA-seq关键词'
        },
        {
            'query': 'organism:"Homo sapiens" AND (title:"single-cell" OR description:"single-cell")',
            'pageSize': 100,
            'name': 'single-cell关键词'
        },
        {
            'query': 'organism:"Homo sapiens" AND (title:"10x genomics" OR description:"10x genomics")',
            'pageSize': 100,
            'name': '10x Genomics技术'
        },
        {
            'query': 'organism:"Homo sapiens" AND (title:"drop-seq" OR description:"drop-seq")',
            'pageSize': 100,
            'name': 'Drop-seq技术'
        },
        {
            'query': 'organism:"Homo sapiens" AND (title:"single nucleus" OR description:"single nucleus")',
            'pageSize': 100,
            'name': '单核测序'
        },
        # 策略4: 只限定ArrayExpress人类数据 (最宽泛，作为补充)
        {
            'query': 'organism:"Homo sapiens" AND accession:E-*',
            'pageSize': 100,
            'name': 'ArrayExpress所有人类数据'
        }
    ]
    
    for strategy in search_strategies:
        page = 0
        strategy_count = 0
        
        while True:
            try:
                params = {
                    'query': strategy['query'],
                    'pageSize': strategy['pageSize'],
                    'page': page + 1  # BioStudies API使用1-based页码
                }
                
                logging.info(f"[{strategy['name']}] 第{page+1}页...")
                response = requests.get(BIOSTUDIES_SEARCH, params=params, timeout=30)
                response.raise_for_status()
                
                data = response.json()
                hits = data.get('hits', [])
                
                if not hits:
                    logging.info(f"  [{strategy['name']}] 第{page+1}页无结果，该策略完成")
                    break
                
                new_count = 0
                for hit in hits:
                    accession = hit.get('accession')
                    
                    # 只收集ArrayExpress数据（E-开头）
                    if accession and accession.startswith('E-') and accession not in all_accessions:
                        all_accessions.add(accession)
                        all_studies.append(hit)
                        new_count += 1
                        strategy_count += 1
                
                logging.info(f"  本页新增 {new_count} 个，该策略累计 {strategy_count} 个，总计 {len(all_accessions)} 个")
                
                # 检查是否还有更多页
                total_hits = data.get('totalHits', 0)
                if (page + 1) * strategy['pageSize'] >= total_hits:
                    logging.info(f"  [{strategy['name']}] 已获取所有结果")
                    break
                
                page += 1
                time.sleep(0.5)  # 避免请求过快
                
            except Exception as e:
                logging.error(f"[{strategy['name']}] 第{page+1}页失败: {e}")
                break
        
        logging.info(f"✓ [{strategy['name']}] 完成，贡献 {strategy_count} 个实验\n")
        time.sleep(1)  # 策略之间稍微延迟
    
    logging.info(f"\n{'='*60}")
    logging.info(f"BioStudies收集完成，共 {len(all_studies)} 个唯一ArrayExpress实验")
    logging.info(f"{'='*60}\n")
    
    with open(RAW_BS_JSON, 'w', encoding='utf-8') as f:
        json.dump(all_studies, f, ensure_ascii=False, indent=2)
    
    return all_studies


# ===== 第二步：收集SCEA数据 =====
def collect_scea_metadata():
    """从SCEA收集策展后的单细胞数据"""
    if os.path.exists(RAW_SCEA_JSON):
        logging.info(f"'{RAW_SCEA_JSON}' 已存在，跳过收集。")
        with open(RAW_SCEA_JSON, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    logging.info("开始从Single Cell Expression Atlas收集元数据...")
    
    try:
        response = requests.get(SCEA_API, params={'species': 'Homo sapiens'}, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        experiments = data.get('experiments', [])
        
        logging.info(f"SCEA收集完成，共 {len(experiments)} 个策展实验")
        
        with open(RAW_SCEA_JSON, 'w', encoding='utf-8') as f:
            json.dump(experiments, f, ensure_ascii=False, indent=2)
        
        return experiments
        
    except Exception as e:
        logging.error(f"SCEA收集失败: {e}")
        return []


# ===== 第三步：获取单个研究的详细信息 =====
def fetch_study_details(accession):
    """
    从BioStudies获取单个研究的完整详细信息
    返回: (study_info, samples_list, files_list)
    """
    try:
        study_url = f"{BIOSTUDIES_STUDY}/{accession}"
        response = requests.get(study_url, timeout=30)
        response.raise_for_status()
        
        study_data = response.json()
        
        if not study_data:
            logging.warning(f"研究 {accession} 无详细信息")
            return None, [], []
        
        # *** 修复：处理返回数据可能是列表或字典的情况 ***
        if isinstance(study_data, list):
            if len(study_data) == 0:
                logging.warning(f"研究 {accession} 返回空列表")
                return None, [], []
            study_data = study_data[0]  # 取第一个元素
        
        if not isinstance(study_data, dict):
            logging.warning(f"研究 {accession} 返回了意外的数据类型: {type(study_data)}")
            return None, [], []
        
        # 1. 提取研究级信息
        section = study_data.get('section', {})
        
        # 如果section也可能是列表
        if isinstance(section, list):
            if len(section) == 0:
                logging.warning(f"研究 {accession} section为空列表")
                return None, [], []
            section = section[0]
        
        attributes = {attr.get('name'): attr.get('value') 
                     for attr in section.get('attributes', []) 
                     if isinstance(attr, dict)}
        
        study_info = {
            'accession': accession,
            'title': section.get('title', ''),
            'description': attributes.get('Description', ''),
            'releasedate': attributes.get('ReleaseDate', ''),
            'attachto': section.get('accno', '')
        }
        
        # 提取实验类型
        study_type = attributes.get('Study type', attributes.get('Experiment type', ''))
        study_info['experiment_types'] = study_type
        
        # 提取发表信息
        links = section.get('links', [])
        pubmed_ids = []
        dois = []
        
        for link in links:
            if not isinstance(link, dict):
                continue
            
            link_url = link.get('url', '')
            if 'pubmed' in link_url.lower() or link.get('type') == 'pmid':
                pmid = re.search(r'(\d+)', link_url)
                if pmid:
                    pubmed_ids.append(f"PMID:{pmid.group(1)}")
            elif 'doi' in link_url.lower() or link.get('type') == 'doi':
                doi = link_url.replace('https://doi.org/', '')
                dois.append(f"DOI:{doi}")
        
        study_info['pubmed_ids'] = '; '.join(pubmed_ids)
        study_info['doi'] = '; '.join(dois)
        
        # 提取联系人信息
        authors = section.get('authors', [])
        contact_names = []
        contact_emails = []
        contact_institutes = []
        
        for author in authors:
            if not isinstance(author, dict):
                continue
            if author.get('name'):
                contact_names.append(author['name'])
            if author.get('email'):
                contact_emails.append(author['email'])
            if author.get('affiliation'):
                contact_institutes.append(author['affiliation'])
        
        study_info['contact_names'] = '; '.join(contact_names)
        study_info['contact_emails'] = '; '.join(contact_emails)
        study_info['contact_institutes'] = '; '.join(contact_institutes)
        study_info['contact_roles'] = ''
        
        # 提取protocol信息
        protocols = []
        protocol_types = []
        
        for subsection in section.get('subsections', []):
            if not isinstance(subsection, dict):
                continue
            if subsection.get('type', '').lower() == 'protocols':
                for attr in subsection.get('attributes', []):
                    if not isinstance(attr, dict):
                        continue
                    if attr.get('name') == 'Protocol Type':
                        protocol_types.append(attr.get('value', ''))
                    elif attr.get('name') == 'Protocol Description':
                        protocols.append(attr.get('value', '')[:200])
        
        study_info['protocol_types'] = '; '.join(protocol_types)
        study_info['protocols'] = ' | '.join(protocols)
        
        # 2. 提取样本信息
        samples_list = []
        
        for subsection in section.get('subsections', []):
            if not isinstance(subsection, dict):
                continue
                
            subsection_type = subsection.get('type', '').lower()
            
            if subsection_type in ['samples', 'sample', 'study samples']:
                for table in subsection.get('tables', []):
                    if not isinstance(table, dict):
                        continue
                        
                    for row in table.get('rows', []):
                        if not isinstance(row, dict):
                            continue
                            
                        sample_info = {
                            'experiment_accession': accession,
                            'sample_accession': '',
                            'sample_name': '',
                            'source_name': ''
                        }
                        
                        row_data = {}
                        for cell in row.get('cells', []):
                            if not isinstance(cell, dict):
                                continue
                            col_name = cell.get('name', '')
                            col_value = cell.get('value', '')
                            if col_name and col_value:
                                row_data[col_name.lower()] = col_value
                        
                        # 如果没有提取到任何行数据，跳过
                        if not row_data:
                            continue
                        
                        sample_info['sample_accession'] = row_data.get('sample name', row_data.get('source name', ''))
                        sample_info['sample_name'] = row_data.get('sample name', '')
                        sample_info['source_name'] = row_data.get('source name', '')
                        
                        sample_info['organism'] = row_data.get('organism', row_data.get('species', ''))
                        sample_info['organism_part'] = row_data.get('organism part', row_data.get('tissue', ''))
                        sample_info['cell_type'] = row_data.get('cell type', '')
                        sample_info['cell_line'] = row_data.get('cell line', '')
                        
                        disease_raw = row_data.get('disease', row_data.get('disease state', ''))
                        disease_category, disease_specific = standardize_disease(disease_raw)
                        sample_info['disease_category'] = disease_category
                        sample_info['disease_specific'] = disease_specific
                        sample_info['disease_raw'] = disease_raw
                        
                        sample_info['sex'] = standardize_sex(row_data.get('sex', ''))
                        sample_info['age_raw'] = row_data.get('age', '')
                        age_num, age_unit = parse_age(sample_info['age_raw'])
                        sample_info['age_value'] = age_num
                        sample_info['age_unit'] = age_unit
                        
                        sample_info['developmental_stage'] = row_data.get('developmental stage', '')
                        sample_info['ethnicity'] = row_data.get('ethnicity', row_data.get('ethnic group', ''))
                        
                        sample_info['clinical_stage'] = row_data.get('clinical stage', row_data.get('tumor stage', ''))
                        sample_info['treatment'] = row_data.get('treatment', row_data.get('compound', ''))
                        sample_info['dose'] = row_data.get('dose', '')
                        sample_info['time_point'] = row_data.get('time point', row_data.get('time', ''))
                        
                        sample_info['sampling_site'] = row_data.get('sampling site', '')
                        sample_info['specimen_type'] = row_data.get('specimen type', '')
                        sample_info['growth_condition'] = row_data.get('growth condition', '')
                        
                        instrument_raw = row_data.get('instrument model', row_data.get('platform', ''))
                        vendor, model = standardize_sequencing_platform(instrument_raw)
                        sample_info['instrument_vendor'] = vendor
                        sample_info['instrument_model'] = model
                        sample_info['instrument_model_raw'] = instrument_raw
                        
                        sample_info['library_strategy'] = row_data.get('library strategy', '')
                        sample_info['library_source'] = row_data.get('library source', '')
                        sample_info['library_selection'] = row_data.get('library selection', '')
                        
                        sample_info['individual_id'] = row_data.get('individual', row_data.get('donor', ''))
                        
                        other_chars = []
                        processed_keys = {
                            'sample name', 'source name', 'organism', 'species', 'organism part',
                            'tissue', 'cell type', 'cell line', 'disease', 'disease state',
                            'sex', 'age', 'developmental stage', 'ethnicity', 'ethnic group',
                            'clinical stage', 'tumor stage', 'treatment', 'compound', 'dose',
                            'time point', 'time', 'sampling site', 'specimen type',
                            'growth condition', 'instrument model', 'platform', 'library strategy',
                            'library source', 'library selection', 'individual', 'donor'
                        }
                        
                        for key, value in row_data.items():
                            if key not in processed_keys and value:
                                other_chars.append(f"{key}: {value}")
                        
                        sample_info['other_characteristics'] = ' | '.join(other_chars)
                        
                        samples_list.append(sample_info)
        
        # 3. 提取文件信息
        files_list = []
        files_section = section.get('files', [])
        
        for file_item in files_section:
            if not isinstance(file_item, dict):
                continue
                
            file_size_raw = file_item.get('size', 0)
            file_size_mb = int(file_size_raw) / (1024 * 1024) if file_size_raw else 0
            
            file_path = file_item.get('path', '')
            file_name = os.path.basename(file_path)
            
            file_info = {
                'experiment_accession': accession,
                'file_name': file_name,
                'file_path': file_path,
                'file_type': file_item.get('type', ''),
                'file_kind': '',
                'file_size_bytes': int(file_size_raw) if file_size_raw else 0,
                'file_size_mb': round(file_size_mb, 2),
                'file_url': f"{AE_FILES_BASE}/{accession}/{file_path}",
                'is_downloadable': 'Yes'
            }
            
            file_info['file_category'] = classify_file_type(file_name)
            
            # 提取文件描述
            file_info['file_description'] = ''
            for attr in file_item.get('attributes', []):
                if not isinstance(attr, dict):
                    continue
                if attr.get('name') == 'Description':
                    file_info['file_description'] = attr.get('value', '')
            
            files_list.append(file_info)
        
        return study_info, samples_list, files_list
        
    except requests.exceptions.RequestException as e:
        logging.error(f"网络请求失败 {accession}: {e}")
        return None, [], []
    except Exception as e:
        logging.error(f"获取研究详情失败 {accession}: {e}")
        import traceback
        logging.debug(traceback.format_exc())
        return None, [], []


# ===== 第四步：获取ENA测序Run信息 =====
def fetch_ena_runs(accession):
    """从ENA获取测序Run的详细信息"""
    runs_list = []
    
    try:
        params = {
            'accession': accession,
            'result': 'read_run',
            'fields': 'run_accession,experiment_accession,sample_accession,instrument_platform,instrument_model,library_name,library_strategy,library_source,library_selection,read_count,base_count,fastq_ftp,fastq_aspera,fastq_md5,submitted_ftp,submitted_aspera,submitted_md5',
            'format': 'json',
            'limit': 0
        }
        
        response = requests.get(f"{ENA_API_BASE}/filereport", params=params, timeout=30)
        
        if response.status_code == 200 and response.text.strip():
            runs_data = response.json()
            
            for run in runs_data:
                vendor, model = standardize_sequencing_platform(run.get('instrument_model', ''))
                
                fastq_ftp = run.get('fastq_ftp', '')
                fastq_aspera = run.get('fastq_aspera', '')
                submitted_ftp = run.get('submitted_ftp', '')
                
                has_fastq = bool(fastq_ftp)
                has_submitted = bool(submitted_ftp)
                publicly_available = has_fastq or has_submitted
                
                fastq_urls = []
                if fastq_ftp:
                    fastq_files = fastq_ftp.split(';')
                    fastq_urls = [f"http://{f}" if not f.startswith('ftp://') else f for f in fastq_files]
                
                run_info = {
                    'experiment_accession': accession,
                    'run_accession': run.get('run_accession', ''),
                    'sample_accession': run.get('sample_accession', ''),
                    'instrument_platform': run.get('instrument_platform', ''),
                    'instrument_vendor': vendor,
                    'instrument_model': model,
                    'instrument_model_raw': run.get('instrument_model', ''),
                    'library_name': run.get('library_name', ''),
                    'library_strategy': run.get('library_strategy', ''),
                    'library_source': run.get('library_source', ''),
                    'library_selection': run.get('library_selection', ''),
                    'read_count': run.get('read_count', 0),
                    'base_count': run.get('base_count', 0),
                    'publicly_available': 'Yes' if publicly_available else 'No',
                    'has_fastq': 'Yes' if has_fastq else 'No',
                    'has_submitted_files': 'Yes' if has_submitted else 'No',
                    'fastq_ftp_urls': '; '.join(fastq_urls) if fastq_urls else '',
                    'fastq_aspera': fastq_aspera,
                    'fastq_md5': run.get('fastq_md5', ''),
                    'submitted_ftp': submitted_ftp,
                    'submitted_aspera': run.get('submitted_aspera', ''),
                    'submitted_md5': run.get('submitted_md5', ''),
                    'ena_browser_url': f"https://www.ebi.ac.uk/ena/browser/view/{run.get('run_accession', '')}"
                }
                
                runs_list.append(run_info)
        
        return runs_list
        
    except Exception as e:
        logging.warning(f"获取ENA runs失败 {accession}: {e}")
        return []


# ===== 第五步：获取出版物详细信息 =====
def fetch_publication_details(pubmed_ids, dois):
    """从PubMed和Europe PMC获取出版物详细信息"""
    publications = []
    
    if pubmed_ids:
        pmid_list = [p.replace('PMID:', '').strip() for p in pubmed_ids.split(';') if p.strip()]
        
        for pmid in pmid_list:
            try:
                response = requests.get(
                    EUROPEPMC_API,
                    params={'query': f'ext_id:{pmid}', 'format': 'json'},
                    timeout=15
                )
                
                if response.status_code == 200:
                    data = response.json()
                    results = data.get('resultList', {}).get('result', [])
                    
                    if results:
                        pub = results[0]
                        
                        authors = pub.get('authorString', '')
                        first_author = authors.split(',')[0] if authors else ''
                        
                        pub_info = {
                            'pubmed_id': pmid,
                            'doi': pub.get('doi', ''),
                            'title': pub.get('title', ''),
                            'authors': authors,
                            'first_author': first_author,
                            'journal': pub.get('journalTitle', ''),
                            'journal_abbr': pub.get('journalAbbreviation', ''),
                            'publication_year': pub.get('pubYear', ''),
                            'publication_date': pub.get('firstPublicationDate', ''),
                            'abstract': pub.get('abstractText', ''),
                            'citation_count': pub.get('citedByCount', 0),
                            'is_open_access': pub.get('isOpenAccess', 'N'),
                            'pub_type': '; '.join(pub.get('pubTypeList', {}).get('pubType', [])),
                            'pmcid': pub.get('pmcid', ''),
                            'pubmed_url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                            'pmc_url': f"https://europepmc.org/article/MED/{pmid}",
                            'full_text_url': pub.get('fullTextUrlList', {}).get('fullTextUrl', [{}])[0].get('url', '') if pub.get('fullTextUrlList') else ''
                        }
                        
                        publications.append(pub_info)
                
                time.sleep(0.3)
                
            except Exception as e:
                logging.warning(f"获取PubMed详情失败 {pmid}: {e}")
    
    if dois and not publications:
        doi_list = [d.replace('DOI:', '').strip() for d in dois.split(';') if d.strip()]
        
        for doi in doi_list:
            try:
                response = requests.get(
                    EUROPEPMC_API,
                    params={'query': f'doi:{doi}', 'format': 'json'},
                    timeout=15
                )
                
                if response.status_code == 200:
                    data = response.json()
                    results = data.get('resultList', {}).get('result', [])
                    
                    if results:
                        pub = results[0]
                        
                        pub_info = {
                            'pubmed_id': pub.get('pmid', ''),
                            'doi': doi,
                            'title': pub.get('title', ''),
                            'authors': pub.get('authorString', ''),
                            'first_author': pub.get('authorString', '').split(',')[0] if pub.get('authorString') else '',
                            'journal': pub.get('journalTitle', ''),
                            'journal_abbr': pub.get('journalAbbreviation', ''),
                            'publication_year': pub.get('pubYear', ''),
                            'publication_date': pub.get('firstPublicationDate', ''),
                            'abstract': pub.get('abstractText', ''),
                            'citation_count': pub.get('citedByCount', 0),
                            'is_open_access': pub.get('isOpenAccess', 'N'),
                            'pub_type': '; '.join(pub.get('pubTypeList', {}).get('pubType', [])),
                            'pmcid': pub.get('pmcid', ''),
                            'pubmed_url': f"https://pubmed.ncbi.nlm.nih.gov/{pub.get('pmid', '')}/",
                            'pmc_url': f"https://europepmc.org/article/MED/{pub.get('pmid', '')}",
                            'full_text_url': ''
                        }
                        
                        publications.append(pub_info)
                
                time.sleep(0.3)
                
            except Exception as e:
                logging.warning(f"获取DOI详情失败 {doi}: {e}")
    
    return publications


# ===== 文件分类 =====
def classify_file_type(filename):
    """详细的文件类型分类"""
    filename_lower = filename.lower()
    
    if any(ext in filename_lower for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']):
        return 'raw_fastq'
    if any(ext in filename_lower for ext in ['.bam', '.cram', '.sam']):
        return 'aligned_reads'
    if any(ext in filename_lower for ext in ['.h5ad', '_h5ad.gz']):
        return 'processed_h5ad'
    if '.h5' in filename_lower and 'filtered' in filename_lower:
        return 'processed_10x_h5_filtered'
    if '.h5' in filename_lower and 'raw' in filename_lower:
        return 'processed_10x_h5_raw'
    if '.h5' in filename_lower:
        return 'processed_h5_generic'
    if '.loom' in filename_lower:
        return 'processed_loom'
    if '.rds' in filename_lower or '.rdata' in filename_lower:
        return 'processed_r_object'
    if 'matrix.mtx' in filename_lower:
        return 'processed_matrix_market'
    if 'barcodes.tsv' in filename_lower:
        return 'processed_10x_barcodes'
    if 'features.tsv' in filename_lower or 'genes.tsv' in filename_lower:
        return 'processed_10x_features'
    if any(pattern in filename_lower for pattern in ['counts.csv', 'counts.tsv', 'expression.txt', 'counts.txt']):
        return 'processed_counts_table'
    if any(pattern in filename_lower for pattern in ['tpm.csv', 'fpkm.csv', 'normalized']):
        return 'processed_normalized'
    if '.sdrf.txt' in filename_lower:
        return 'metadata_sdrf'
    if '.idf.txt' in filename_lower:
        return 'metadata_idf'
    if 'metadata' in filename_lower or 'annotation' in filename_lower:
        return 'metadata_general'
    if any(pattern in filename_lower for pattern in ['cluster', 'tsne', 'umap', 'pca', 'embedding']):
        return 'analysis_coordinates'
    if 'marker' in filename_lower or 'deg' in filename_lower or 'differential' in filename_lower:
        return 'analysis_markers'
    
    return 'other'


# ===== 处理和整合所有元数据 (其余代码保持不变) =====
def process_all_metadata(bs_studies, scea_experiments):
    """整合BioStudies和SCEA的数据，生成6层表格"""
    logging.info("开始处理和整合元数据...")
    
    scea_dict = {exp.get('experimentAccession'): exp for exp in scea_experiments}
    
    experiment_records = []
    sample_records = []
    file_records = []
    sample_file_mappings = []
    sequencing_runs = []
    all_publications = []
    
    stats = {
        'total_experiments': 0,
        'experiments_with_samples': 0,
        'experiments_with_publications': 0,
        'experiments_with_sequencing_runs': 0,
        'total_samples': 0,
        'total_files': 0,
        'total_runs': 0,
        'total_publications': 0,
        'in_scea': 0,
        'disease_distribution': defaultdict(int),
        'tissue_distribution': defaultdict(int),
        'sex_distribution': defaultdict(int),
        'platform_distribution': defaultdict(int),
        'file_type_distribution': defaultdict(int),
        'publicly_available_runs': 0
    }
    
    processed_publications = set()
    
    for study in tqdm(bs_studies, desc="Processing Studies"):
        accession = study.get('accession')
        if not accession or not accession.startswith('E-'):
            continue
        
        stats['total_experiments'] += 1
        
        study_info, samples_list, files_list = fetch_study_details(accession)
        
        if not study_info:
            continue
        
        in_scea = accession in scea_dict
        scea_info = scea_dict.get(accession, {})
        
        if in_scea:
            stats['in_scea'] += 1
        
        runs_list = fetch_ena_runs(accession)
        if runs_list:
            stats['experiments_with_sequencing_runs'] += 1
            stats['total_runs'] += len(runs_list)
            
            for run in runs_list:
                sequencing_runs.append(run)
                
                if run.get('instrument_vendor'):
                    stats['platform_distribution'][run['instrument_vendor']] += 1
                
                if run.get('publicly_available') == 'Yes':
                    stats['publicly_available_runs'] += 1
        
        pubmed_ids = study_info.get('pubmed_ids', '')
        dois = study_info.get('doi', '')
        
        if pubmed_ids or dois:
            stats['experiments_with_publications'] += 1
            publications = fetch_publication_details(pubmed_ids, dois)
            
            for pub in publications:
                pub_key = pub.get('pubmed_id') or pub.get('doi')
                if pub_key and pub_key not in processed_publications:
                    processed_publications.add(pub_key)
                    pub['experiment_accession'] = accession
                    all_publications.append(pub)
                    stats['total_publications'] += 1
        
        sample_count = len(samples_list)
        if sample_count > 0:
            stats['experiments_with_samples'] += 1
        
        all_diseases = set()
        all_tissues = set()
        all_sexes = set()
        all_ages = []
        all_platforms = set()
        
        for sample in samples_list:
            if sample.get('disease_category'):
                all_diseases.add(sample['disease_category'])
            if sample.get('organism_part'):
                all_tissues.add(sample['organism_part'])
            if sample.get('sex'):
                all_sexes.add(sample['sex'])
            if sample.get('age_value'):
                all_ages.append(sample['age_value'])
            if sample.get('instrument_vendor'):
                all_platforms.add(sample['instrument_vendor'])
        
        for run in runs_list:
            if run.get('instrument_vendor'):
                all_platforms.add(run['instrument_vendor'])
        
        file_count = len(files_list)
        processed_files = [f for f in files_list if f['file_category'].startswith('processed_')]
        has_h5ad = any(f['file_category'] == 'processed_h5ad' for f in files_list)
        has_10x = any('10x' in f['file_category'] for f in files_list)
        has_matrix = any('matrix' in f['file_category'] or 'counts' in f['file_category'] for f in files_list)
        
        total_file_size_gb = sum(f.get('file_size_mb', 0) for f in files_list) / 1024
        
        exp_record = {
            'experiment_accession': accession,
            'experiment_title': study_info.get('title', ''),
            'experiment_description': study_info.get('description', ''),
            'experiment_types': study_info.get('experiment_types', ''),
            
            'total_samples': sample_count,
            'disease_categories': '; '.join(sorted(all_diseases)),
            'tissues': '; '.join(sorted(all_tissues)),
            'sexes': '; '.join(sorted(all_sexes)),
            'age_range': f"{min(all_ages)}-{max(all_ages)}" if all_ages else '',
            
            'sequencing_platforms': '; '.join(sorted(all_platforms)),
            'total_sequencing_runs': len(runs_list),
            'has_public_raw_data': 'Yes' if any(r.get('publicly_available') == 'Yes' for r in runs_list) else 'No',
            
            'total_files': file_count,
            'total_file_size_gb': round(total_file_size_gb, 2),
            'processed_file_count': len(processed_files),
            'has_h5ad': 'Yes' if has_h5ad else 'No',
            'has_10x_format': 'Yes' if has_10x else 'No',
            'has_matrix': 'Yes' if has_matrix else 'No',
            
            'in_scea': 'Yes' if in_scea else 'No',
            'scea_cell_count': scea_info.get('numberOfCells', 0) if in_scea else 0,
            'scea_url': f"https://www.ebi.ac.uk/gxa/sc/experiments/{accession}" if in_scea else '',
            
            'has_publication': 'Yes' if (pubmed_ids or dois) else 'No',
            'pubmed_ids': pubmed_ids,
            'doi': dois,
            'publication_count': len(publications) if pubmed_ids or dois else 0,
            
            'contact_names': study_info.get('contact_names', ''),
            'contact_emails': study_info.get('contact_emails', ''),
            'contact_institutes': study_info.get('contact_institutes', ''),
            'contact_roles': study_info.get('contact_roles', ''),
            
            'release_date': study_info.get('releasedate', ''),
            'submission_date': '',
            'last_update_date': '',
            
            'biostudies_url': f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}",
            'ena_url': f"https://www.ebi.ac.uk/ena/browser/view/{accession}",
            
            'protocol_types': study_info.get('protocol_types', ''),
            'protocols_summary': study_info.get('protocols', '')[:500]
        }
        
        experiment_records.append(exp_record)
        
        for sample in samples_list:
            stats['total_samples'] += 1
            
            if sample.get('disease_category'):
                stats['disease_distribution'][sample['disease_category']] += 1
            if sample.get('organism_part'):
                stats['tissue_distribution'][sample['organism_part']] += 1
            if sample.get('sex'):
                stats['sex_distribution'][sample['sex']] += 1
            
            sample_records.append(sample)
        
        for file_info in files_list:
            stats['total_files'] += 1
            stats['file_type_distribution'][file_info['file_category']] += 1
            
            file_records.append(file_info)
            
            sample_acc_in_filename = None
            for sample in samples_list:
                sample_acc = sample.get('sample_accession', '')
                if sample_acc and sample_acc in file_info['file_name']:
                    sample_acc_in_filename = sample_acc
                    break
            
            if sample_acc_in_filename:
                sample_file_mappings.append({
                    'experiment_accession': accession,
                    'sample_accession': sample_acc_in_filename,
                    'file_name': file_info['file_name'],
                    'file_category': file_info['file_category']
                })
        
        time.sleep(0.3)
    
    logging.info("\n正在保存数据表...")
    
    df_experiments = pd.DataFrame(experiment_records)
    df_experiments.to_csv(EXPERIMENT_METADATA_CSV, index=False, encoding='utf-8-sig')
    logging.info(f"✓ 实验元数据已保存: {EXPERIMENT_METADATA_CSV}")
    
    df_samples = pd.DataFrame(sample_records)
    df_samples = df_samples.applymap(lambda x: None if x == '' else x)
    df_samples.to_csv(SAMPLE_METADATA_CSV, index=False, encoding='utf-8-sig')
    logging.info(f"✓ 样本元数据已保存: {SAMPLE_METADATA_CSV} (包含 {len(df_samples)} 个样本)")
    
    df_files = pd.DataFrame(file_records)
    df_files.to_csv(FILE_MANIFEST_CSV, index=False, encoding='utf-8-sig')
    logging.info(f"✓ 文件清单已保存: {FILE_MANIFEST_CSV}")
    
    if sample_file_mappings:
        df_mapping = pd.DataFrame(sample_file_mappings)
        df_mapping.to_csv(SAMPLE_FILE_MAPPING_CSV, index=False, encoding='utf-8-sig')
        logging.info(f"✓ 样本-文件关联已保存: {SAMPLE_FILE_MAPPING_CSV}")
    
    if sequencing_runs:
        df_runs = pd.DataFrame(sequencing_runs)
        df_runs.to_csv(SEQUENCING_RUNS_CSV, index=False, encoding='utf-8-sig')
        logging.info(f"✓ 测序Run详情已保存: {SEQUENCING_RUNS_CSV} (包含 {len(df_runs)} 个runs)")
    
    if all_publications:
        df_pubs = pd.DataFrame(all_publications)
        df_pubs.to_csv(PUBLICATIONS_CSV, index=False, encoding='utf-8-sig')
        logging.info(f"✓ 出版物详情已保存: {PUBLICATIONS_CSV} (包含 {len(df_pubs)} 篇文献)")
    
    generate_statistics_report(stats, df_experiments, df_samples, df_files, 
                              sequencing_runs if sequencing_runs else [], 
                              all_publications if all_publications else [])


def generate_statistics_report(stats, df_experiments, df_samples, df_files, sequencing_runs, publications):
    """生成详细的统计报告"""
    report_lines = []
    
    report_lines.append("=" * 80)
    report_lines.append("BioStudies (ArrayExpress) 人源单细胞 RNA-seq 数据收集统计报告")
    report_lines.append("=" * 80)
    report_lines.append("")
    
    report_lines.append("【总体统计】")
    report_lines.append(f"  总实验数: {stats['total_experiments']}")
    report_lines.append(f"  有样本信息的实验: {stats['experiments_with_samples']}")
    report_lines.append(f"  有出版物的实验: {stats['experiments_with_publications']} ({stats['experiments_with_publications']/max(stats['total_experiments'], 1)*100:.1f}%)")
    report_lines.append(f"  有测序Run的实验: {stats['experiments_with_sequencing_runs']} ({stats['experiments_with_sequencing_runs']/max(stats['total_experiments'], 1)*100:.1f}%)")
    report_lines.append(f"  在SCEA中的实验: {stats['in_scea']} ({stats['in_scea']/max(stats['total_experiments'], 1)*100:.1f}%)")
    report_lines.append(f"  总样本数: {stats['total_samples']}")
    report_lines.append(f"  总文件数: {stats['total_files']}")
    report_lines.append(f"  总测序Runs: {stats['total_runs']}")
    report_lines.append(f"  总出版物: {stats['total_publications']}")
    report_lines.append("")
    
    report_lines.append("【测序平台分布】")
    platform_sorted = sorted(stats['platform_distribution'].items(), key=lambda x: x[1], reverse=True)
    for platform, count in platform_sorted:
        report_lines.append(f"  {platform}: {count} runs")
    report_lines.append(f"  公开可下载的Runs: {stats['publicly_available_runs']} ({stats['publicly_available_runs']/max(stats['total_runs'], 1)*100:.1f}%)")
    report_lines.append("")
    
    report_lines.append("【数据可用性】")
    if not df_experiments.empty:
        report_lines.append(f"  有原始数据(FASTQ): {(df_experiments['has_public_raw_data'] == 'Yes').sum()} 个实验")
        report_lines.append(f"  有h5ad格式: {(df_experiments['has_h5ad'] == 'Yes').sum()} 个实验")
        report_lines.append(f"  有10x格式: {(df_experiments['has_10x_format'] == 'Yes').sum()} 个实验")
        report_lines.append(f"  有矩阵文件: {(df_experiments['has_matrix'] == 'Yes').sum()} 个实验")
        
        if 'total_file_size_gb' in df_experiments.columns:
            total_size = df_experiments['total_file_size_gb'].sum()
            avg_size = df_experiments['total_file_size_gb'].mean()
            report_lines.append(f"  总文件大小: {total_size:.2f} GB")
            report_lines.append(f"  平均实验大小: {avg_size:.2f} GB")
    report_lines.append("")
    
    if publications:
        report_lines.append("【出版物统计】")
        report_lines.append(f"  总出版物数: {len(publications)}")
        
        df_pubs = pd.DataFrame(publications)
        if 'publication_year' in df_pubs.columns:
            year_counts = df_pubs['publication_year'].value_counts().sort_index(ascending=False).head(10)
            report_lines.append("  发表年份分布 (Top 10):")
            for year, count in year_counts.items():
                report_lines.append(f"    {year}: {count} 篇")
        
        if 'is_open_access' in df_pubs.columns:
            oa_count = (df_pubs['is_open_access'] == 'Y').sum()
            report_lines.append(f"  开放获取文献: {oa_count} ({oa_count/len(df_pubs)*100:.1f}%)")
        
        if 'citation_count' in df_pubs.columns:
            avg_citations = df_pubs['citation_count'].mean()
            max_citations = df_pubs['citation_count'].max()
            report_lines.append(f"  平均引用数: {avg_citations:.1f}")
            report_lines.append(f"  最高引用数: {int(max_citations)}")
        
        if 'journal_abbr' in df_pubs.columns:
            top_journals = df_pubs['journal_abbr'].value_counts().head(5)
            report_lines.append("  Top 5 期刊:")
            for journal, count in top_journals.items():
                report_lines.append(f"    {journal}: {count} 篇")
        
        report_lines.append("")
    
    report_lines.append("【疾病分布 Top 10】")
    disease_sorted = sorted(stats['disease_distribution'].items(), key=lambda x: x[1], reverse=True)[:10]
    for disease, count in disease_sorted:
        report_lines.append(f"  {disease}: {count} 样本")
    report_lines.append("")
    
    report_lines.append("【组织分布 Top 10】")
    tissue_sorted = sorted(stats['tissue_distribution'].items(), key=lambda x: x[1], reverse=True)[:10]
    for tissue, count in tissue_sorted:
        report_lines.append(f"  {tissue}: {count} 样本")
    report_lines.append("")
    
    report_lines.append("【性别分布】")
    for sex, count in stats['sex_distribution'].items():
        report_lines.append(f"  {sex}: {count} 样本")
    report_lines.append("")
    
    if not df_samples.empty and 'age_value' in df_samples.columns:
        age_data = df_samples['age_value'].dropna()
        if len(age_data) > 0:
            report_lines.append("【年龄统计】")
            report_lines.append(f"  有年龄信息的样本: {len(age_data)}")
            report_lines.append(f"  年龄范围: {age_data.min():.1f} - {age_data.max():.1f}")
            report_lines.append(f"  平均年龄: {age_data.mean():.1f}")
            report_lines.append(f"  中位年龄: {age_data.median():.1f}")
            report_lines.append("")
    
    report_lines.append("【文件类型分布 Top 15】")
    file_type_sorted = sorted(stats['file_type_distribution'].items(), key=lambda x: x[1], reverse=True)[:15]
    for ftype, count in file_type_sorted:
        report_lines.append(f"  {ftype}: {count} 文件")
    report_lines.append("")
    
    if not df_samples.empty:
        report_lines.append("【样本信息完整性】")
        total_samples = len(df_samples)
        report_lines.append(f"  有疾病信息: {df_samples['disease_category'].notna().sum()} ({df_samples['disease_category'].notna().sum()/total_samples*100:.1f}%)")
        report_lines.append(f"  有组织信息: {df_samples['organism_part'].notna().sum()} ({df_samples['organism_part'].notna().sum()/total_samples*100:.1f}%)")
        report_lines.append(f"  有性别信息: {df_samples['sex'].notna().sum()} ({df_samples['sex'].notna().sum()/total_samples*100:.1f}%)")
        report_lines.append(f"  有年龄信息: {df_samples['age_value'].notna().sum()} ({df_samples['age_value'].notna().sum()/total_samples*100:.1f}%)")
        report_lines.append(f"  有细胞类型: {df_samples['cell_type'].notna().sum()} ({df_samples['cell_type'].notna().sum()/total_samples*100:.1f}%)")
        report_lines.append(f"  有平台信息: {df_samples['instrument_vendor'].notna().sum()} ({df_samples['instrument_vendor'].notna().sum()/total_samples*100:.1f}%)")
        report_lines.append("")
    
    report_lines.append("=" * 80)
    report_lines.append(f"报告生成时间: {pd.Timestamp.now()}")
    report_lines.append("=" * 80)
    
    report_text = '\n'.join(report_lines)
    
    with open(STATISTICS_REPORT, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    logging.info("\n" + report_text)


# ===== 主函数 =====
def main():
    """主流程"""
    logging.info("=" * 80)
    logging.info("开始BioStudies (ArrayExpress) 人源单细胞RNA-seq数据收集")
    logging.info("=" * 80)
    
    bs_studies = collect_biostudies_metadata()
    
    scea_experiments = collect_scea_metadata()
    
    if bs_studies:
        process_all_metadata(bs_studies, scea_experiments)
    else:
        logging.error("没有收集到任何数据")
    
    logging.info("\n" + "=" * 80)
    logging.info("数据收集完成！")
    logging.info("=" * 80)


if __name__ == "__main__":
    main()