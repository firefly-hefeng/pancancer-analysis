# data_cleaner.py
"""
人类单细胞RNA-seq元数据清洗和重组系统
解决数据质量问题并重组为统一格式
"""

import os
import json
import re
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import logging
from datetime import datetime
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from dataclasses import dataclass
from collections import defaultdict

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('data_cleaning.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


@dataclass
class KimiAPIConfig:
    """Kimi API配置"""
    api_keys: List[str]
    api_url: str = "https://api.moonshot.cn/v1/chat/completions"
    model: str = "moonshot-v1-8k"
    current_key_index: int = 0
    
    def get_next_key(self) -> str:
        """获取下一个可用的API key"""
        key = self.api_keys[self.current_key_index]
        self.current_key_index = (self.current_key_index + 1) % len(self.api_keys)
        return key


class KimiAssistant:
    """Kimi AI助手，用于复杂数据识别"""
    
    def __init__(self, api_keys: List[str]):
        self.config = KimiAPIConfig(api_keys=api_keys)
        self.session = requests.Session()
        
    def query(self, prompt: str, max_retries: int = 3) -> Optional[str]:
        """查询Kimi API"""
        for attempt in range(max_retries):
            try:
                api_key = self.config.get_next_key()
                
                headers = {
                    "Content-Type": "application/json",
                    "Authorization": f"Bearer {api_key}"
                }
                
                payload = {
                    "model": self.config.model,
                    "messages": [
                        {
                            "role": "system",
                            "content": "你是一个生物信息学专家,专门识别单细胞RNA测序数据和物种信息。请简洁准确地回答问题。"
                        },
                        {
                            "role": "user",
                            "content": prompt
                        }
                    ],
                    "temperature": 0.3
                }
                
                response = self.session.post(
                    self.config.api_url,
                    headers=headers,
                    json=payload,
                    timeout=30
                )
                
                if response.status_code == 200:
                    result = response.json()
                    return result['choices'][0]['message']['content']
                else:
                    logger.warning(f"Kimi API error (attempt {attempt + 1}): {response.status_code}")
                    
            except Exception as e:
                logger.error(f"Error querying Kimi (attempt {attempt + 1}): {e}")
                
            time.sleep(1)
        
        return None
    
    def is_human_data(self, title: str, summary: str, keywords: str) -> Tuple[bool, str]:
        """判断是否为人类数据"""
        prompt = f"""
请判断以下研究是否为人类(Homo sapiens)的数据:

标题: {title[:200]}
摘要: {summary[:500]}
关键词: {keywords[:200]}

请只回答以下格式:
答案: YES/NO/UNCERTAIN
理由: 简短说明(一句话)
"""
        
        response = self.query(prompt)
        if not response:
            return None, "API调用失败"
        
        if "YES" in response.upper():
            return True, response.split("理由:")[-1].strip() if "理由:" in response else ""
        elif "NO" in response.upper():
            return False, response.split("理由:")[-1].strip() if "理由:" in response else ""
        else:
            return None, response.split("理由:")[-1].strip() if "理由:" in response else ""
    
    def is_single_cell_rna_seq(self, title: str, summary: str, platform: str, 
                                experiment_design: str) -> Tuple[bool, str]:
        """判断是否为单细胞RNA测序数据"""
        prompt = f"""
请判断以下研究是否为单细胞RNA测序(scRNA-seq)数据,排除bulk RNA-seq、ATAC-seq、ChIP-seq等其他技术:

标题: {title[:200]}
摘要: {summary[:500]}
测序平台: {platform[:200]}
实验设计: {experiment_design[:200]}

请只回答以下格式:
答案: YES/NO/UNCERTAIN
技术类型: 识别出的具体技术(如scRNA-seq, bulk RNA-seq, ATAC-seq等)
理由: 简短说明(一句话)
"""
        
        response = self.query(prompt)
        if not response:
            return None, "API调用失败"
        
        if "YES" in response.upper():
            return True, response
        elif "NO" in response.upper():
            return False, response
        else:
            return None, response


class DataQualityChecker:
    """数据质量检查器"""
    
    def __init__(self, kimi_assistant: Optional[KimiAssistant] = None):
        self.kimi = kimi_assistant
        
        # 人类物种关键词
        self.human_keywords = {
            'positive': ['human', 'homo sapiens', 'hsa', 'sapiens', 'patient', 
                        'clinical', 'pbmc', 'peripheral blood mononuclear'],
            'negative': ['mouse', 'mice', 'murine', 'mus musculus', 'rat', 
                        'drosophila', 'zebrafish', 'yeast', 'arabidopsis',
                        'c. elegans', 'xenopus', 'monkey', 'macaque']
        }
        
        # 单细胞RNA-seq关键词
        self.scrna_keywords = {
            'positive': ['single cell', 'single-cell', 'scrna-seq', 'scrna',
                        'drop-seq', 'droplet', '10x', 'chromium', 'smart-seq',
                        'cel-seq', 'mars-seq', 'indrops', 'sci-rna-seq'],
            'negative': ['bulk', 'atac-seq', 'chip-seq', 'whole genome',
                        'wgs', 'wes', 'exome', 'metagenom', 'dnase-seq',
                        'hi-c', 'cut&run', 'spatial transcriptomics',
                        'proteomics', 'metabolomics']
        }
        
    def check_human_species(self, row: pd.Series) -> Tuple[str, float, str]:
        """
        检查是否为人类数据
        返回: (判断结果, 置信度, 原因)
        """
        text_fields = [
            str(row.get('title', '')),
            str(row.get('summary', '')),
            str(row.get('supplementary_information', '')),
            str(row.get('tissue', '')),
            str(row.get('sequencing_platform', ''))
        ]
        
        combined_text = ' '.join(text_fields).lower()
        
        # 计算正向和负向匹配分数
        positive_score = sum(1 for kw in self.human_keywords['positive'] 
                           if kw in combined_text)
        negative_score = sum(1 for kw in self.human_keywords['negative'] 
                           if kw in combined_text)
        
        # 检查supplementary_information中的organism字段
        try:
            supp_info = json.loads(str(row.get('supplementary_information', '{}')))
            organism = str(supp_info.get('organism', '')).lower()
            if 'human' in organism or 'homo sapiens' in organism or 'hsa' in organism:
                positive_score += 5
            elif any(neg in organism for neg in self.human_keywords['negative']):
                negative_score += 5
        except:
            pass
        
        # 判断逻辑
        if negative_score > 0 and positive_score == 0:
            return 'NON_HUMAN', 0.9, f'发现非人类物种关键词 (负向分数: {negative_score})'
        
        if positive_score >= 2 and negative_score == 0:
            return 'HUMAN', 0.9, f'明确的人类数据标识 (正向分数: {positive_score})'
        
        if positive_score > negative_score * 2:
            return 'HUMAN', 0.7, f'人类关键词占优 (正向: {positive_score}, 负向: {negative_score})'
        
        if negative_score > positive_score * 2:
            return 'NON_HUMAN', 0.7, f'非人类关键词占优 (正向: {positive_score}, 负向: {negative_score})'
        
        # 需要进一步判断
        if positive_score > 0 or negative_score > 0:
            return 'UNCERTAIN', 0.5, f'混合信号 (正向: {positive_score}, 负向: {negative_score})'
        
        return 'UNCERTAIN', 0.3, '缺乏明确的物种信息'
    
    def check_scrna_seq_type(self, row: pd.Series) -> Tuple[str, float, str]:
        """
        检查是否为单细胞RNA-seq
        返回: (判断结果, 置信度, 原因)
        """
        text_fields = [
            str(row.get('title', '')),
            str(row.get('summary', '')),
            str(row.get('sequencing_platform', '')),
            str(row.get('experiment_design', '')),
            str(row.get('supplementary_information', ''))
        ]
        
        combined_text = ' '.join(text_fields).lower()
        
        # 计算分数
        positive_score = sum(1 for kw in self.scrna_keywords['positive'] 
                           if kw in combined_text)
        negative_score = sum(1 for kw in self.scrna_keywords['negative'] 
                           if kw in combined_text)
        
        # 特殊检查:平台信息
        platform_text = str(row.get('sequencing_platform', '')).lower()
        if any(kw in platform_text for kw in ['10x', 'chromium', 'drop-seq', 'smart-seq']):
            positive_score += 3
        
        if 'bulk' in platform_text or 'illumina hiseq' in combined_text:
            negative_score += 2
        
        # 判断逻辑
        if negative_score > 0 and positive_score == 0:
            return 'NOT_SCRNA', 0.9, f'发现非scRNA-seq技术关键词 (负向分数: {negative_score})'
        
        if positive_score >= 2 and negative_score == 0:
            return 'SCRNA', 0.9, f'明确的scRNA-seq标识 (正向分数: {positive_score})'
        
        if positive_score > negative_score * 2:
            return 'SCRNA', 0.7, f'scRNA-seq关键词占优 (正向: {positive_score}, 负向: {negative_score})'
        
        if negative_score > positive_score:
            return 'NOT_SCRNA', 0.7, f'其他技术关键词占优 (正向: {positive_score}, 负向: {negative_score})'
        
        if positive_score > 0 or negative_score > 0:
            return 'UNCERTAIN', 0.5, f'混合信号 (正向: {positive_score}, 负向: {negative_score})'
        
        return 'UNCERTAIN', 0.3, '缺乏明确的技术类型信息'
    
    def quality_check_with_ai(self, row: pd.Series, check_type: str = 'both') -> Dict:
        """
        使用AI辅助进行质量检查
        check_type: 'species', 'technology', 'both'
        """
        if not self.kimi:
            return {}
        
        result = {}
        
        if check_type in ['species', 'both']:
            is_human, reason = self.kimi.is_human_data(
                str(row.get('title', '')),
                str(row.get('summary', '')),
                str(row.get('supplementary_information', ''))
            )
            
            if is_human is not None:
                result['ai_species_check'] = 'HUMAN' if is_human else 'NON_HUMAN'
                result['ai_species_reason'] = reason
            else:
                result['ai_species_check'] = 'UNCERTAIN'
                result['ai_species_reason'] = reason
        
        if check_type in ['technology', 'both']:
            is_scrna, reason = self.kimi.is_single_cell_rna_seq(
                str(row.get('title', '')),
                str(row.get('summary', '')),
                str(row.get('sequencing_platform', '')),
                str(row.get('experiment_design', ''))
            )
            
            if is_scrna is not None:
                result['ai_technology_check'] = 'SCRNA' if is_scrna else 'NOT_SCRNA'
                result['ai_technology_reason'] = reason
            else:
                result['ai_technology_check'] = 'UNCERTAIN'
                result['ai_technology_reason'] = reason
        
        return result


class DataReorganizer:
    """数据重组器"""
    
    def __init__(self, output_dir: str = "cleaned_metadata"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
    def extract_project_ids(self, row: pd.Series) -> Dict[str, str]:
        """提取项目ID"""
        ids = {
            'Study/Project_id_1': '',
            'Study/Project_id_2': '',
            'Study/Project_id_3': ''
        }
        
        # 从不同字段提取ID
        potential_ids = []
        
        # 主ID
        main_id = str(row.get('id', '')).strip()
        if main_id:
            potential_ids.append(main_id)
        
        # 从链接提取ID
        access_link = str(row.get('access_link', ''))
        if access_link:
            # 提取各种格式的ID
            patterns = [
                r'/([A-Z]+\d+)',  # GEO/SRA格式
                r'accession=([A-Z0-9]+)',
                r'/studies/([A-Z0-9]+)',
                r'/collections/([a-f0-9-]+)',
                r'/dataset/([a-f0-9-]+)',
                r'doi\.org/(.+)$'
            ]
            
            for pattern in patterns:
                matches = re.findall(pattern, access_link)
                potential_ids.extend(matches)
        
        # 从pubmed/DOI提取
        pubmed = str(row.get('pubmed', '')).strip()
        if pubmed and pubmed not in potential_ids:
            potential_ids.append(pubmed)
        
        # 去重并填充
        unique_ids = list(dict.fromkeys([id for id in potential_ids if id]))
        
        for i, id_val in enumerate(unique_ids[:3]):
            ids[f'Study/Project_id_{i+1}'] = id_val
        
        return ids
    
    def determine_data_availability(self, row: pd.Series) -> Dict[str, any]:
        """判断数据可用性"""
        availability = {
            'raw_sample_id': '',
            'matrix_sample_id': '',
            'raw_exist': False,
            'raw_open': False,
            'matrix_exist': False,
            'matrix_open': False,
            'file_type': ''
        }
        
        # 根据数据库类型判断
        source_db = str(row.get('source_database', '')).lower()
        data_tier = str(row.get('data_tier', '')).lower()
        open_status = str(row.get('open_status', '')).lower()
        
        # 判断数据类型
        if 'cellxgene' in source_db or 'broad' in source_db:
            # 这些数据库通常提供处理后的矩阵
            availability['matrix_exist'] = True
            availability['matrix_sample_id'] = str(row.get('id', ''))
            availability['matrix_open'] = 'open' in open_status
            availability['file_type'] = 'h5ad, h5, csv'
            
        elif 'geo' in source_db or 'sra' in source_db:
            # GEO/SRA通常有原始数据
            availability['raw_exist'] = True
            availability['raw_sample_id'] = str(row.get('id', ''))
            availability['raw_open'] = 'open' in open_status
            
            # 可能也有处理后的数据
            if 'processed' in data_tier or 'matrix' in data_tier:
                availability['matrix_exist'] = True
                availability['matrix_sample_id'] = str(row.get('id', ''))
                availability['matrix_open'] = 'open' in open_status
            
            availability['file_type'] = 'fastq, bam, h5'
            
        elif any(db in source_db for db in ['zenodo', 'figshare', 'dryad']):
            # 这些数据库可能有各种类型
            if 'raw' in data_tier:
                availability['raw_exist'] = True
                availability['raw_sample_id'] = str(row.get('id', ''))
                availability['raw_open'] = 'open' in open_status
                availability['file_type'] = 'fastq, bam'
            else:
                availability['matrix_exist'] = True
                availability['matrix_sample_id'] = str(row.get('id', ''))
                availability['matrix_open'] = 'open' in open_status
                availability['file_type'] = 'csv, tsv, h5, h5ad'
        
        return availability
    
    def reorganize_row(self, row: pd.Series) -> Dict[str, any]:
        """重组单行数据"""
        # 提取项目ID
        project_ids = self.extract_project_ids(row)
        
        # 判断数据可用性
        availability = self.determine_data_availability(row)
        
        # 组装新的数据行
        new_row = {
            **project_ids,
            **availability,
            'title': str(row.get('title', '')),
            'disease_general': str(row.get('disease_general', '')),
            'disease': str(row.get('disease', '')),
            'pubmed': str(row.get('pubmed', '')),
            'source_database': str(row.get('source_database', '')),
            'access_link': str(row.get('access_link', '')),
            'open_status': str(row.get('open_status', '')),
            'ethnicity': str(row.get('ethnicity', '')),
            'sex': str(row.get('sex', '')),
            'tissue_location': str(row.get('tissue', '')),
            'sequencing_platform': str(row.get('sequencing_platform', '')),
            'experiment_design': str(row.get('experiment_design', '')),
            'sample_type': self._normalize_sample_type(str(row.get('sample_type', ''))),
            'summary': str(row.get('summary', '')),
            'citation_count': row.get('citation_count', 0),
            'publication_date': str(row.get('publication_date', '')),
            'submission_date': str(row.get('submission_date', '')),
            'last_update_date': str(row.get('last_update_date', '')),
            'contact_name': str(row.get('contact_name', '')),
            'contact_email': str(row.get('contact_email', '')),
            'contact_institute': str(row.get('contact_institute', '')),
            'supplementary_information': str(row.get('supplementary_information', ''))
        }
        
        return new_row
    
    def _normalize_sample_type(self, sample_type: str) -> str:
        """标准化样本类型"""
        sample_type_lower = sample_type.lower()
        
        if any(kw in sample_type_lower for kw in ['healthy', 'normal', 'control']):
            return 'healthy'
        elif any(kw in sample_type_lower for kw in ['tumor', 'cancer', 'malignant']):
            return 'tumor'
        elif any(kw in sample_type_lower for kw in ['disease', 'patient', 'illness']):
            return 'other illness'
        else:
            return sample_type if sample_type else 'unknown'


class MetadataCleaningPipeline:
    """完整的元数据清洗流程"""
    
    def __init__(self, input_dir: str = "metadata_output", 
                 output_dir: str = "cleaned_metadata",
                 kimi_api_keys: Optional[List[str]] = None,
                 use_ai_assistance: bool = True):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # 初始化组件
        self.kimi = KimiAssistant(kimi_api_keys) if kimi_api_keys and use_ai_assistance else None
        self.quality_checker = DataQualityChecker(self.kimi)
        self.reorganizer = DataReorganizer(output_dir)
        
        # 统计信息
        self.stats = {
            'total_records': 0,
            'human_records': 0,
            'non_human_records': 0,
            'uncertain_species': 0,
            'scrna_records': 0,
            'non_scrna_records': 0,
            'uncertain_technology': 0,
            'final_clean_records': 0,
            'flagged_for_review': 0
        }
        
    def load_merged_data(self) -> pd.DataFrame:
        """加载合并后的数据"""
        series_file = self.input_dir / "integrated" / "merged_series_table.csv"
        
        if not series_file.exists():
            raise FileNotFoundError(f"未找到合并后的数据文件: {series_file}")
        
        df = pd.read_csv(series_file)
        logger.info(f"加载了 {len(df)} 条记录")
        self.stats['total_records'] = len(df)
        
        return df
    
    def perform_quality_checks(self, df: pd.DataFrame, 
                              use_ai_for_uncertain: bool = True) -> pd.DataFrame:
        """执行质量检查"""
        logger.info("开始质量检查...")
        
        results = []
        uncertain_indices = []
        
        for idx, row in df.iterrows():
            # 基础检查
            species_result, species_conf, species_reason = self.quality_checker.check_human_species(row)
            tech_result, tech_conf, tech_reason = self.quality_checker.check_scrna_seq_type(row)
            
            result = {
                'index': idx,
                'species_check': species_result,
                'species_confidence': species_conf,
                'species_reason': species_reason,
                'technology_check': tech_result,
                'technology_confidence': tech_conf,
                'technology_reason': tech_reason,
                'needs_ai_review': False
            }
            
            # 判断是否需要AI辅助
            if (species_result == 'UNCERTAIN' or tech_result == 'UNCERTAIN' or
                species_conf < 0.7 or tech_conf < 0.7):
                result['needs_ai_review'] = True
                uncertain_indices.append(idx)
            
            results.append(result)
            
            # 更新统计
            if species_result == 'HUMAN':
                self.stats['human_records'] += 1
            elif species_result == 'NON_HUMAN':
                self.stats['non_human_records'] += 1
            else:
                self.stats['uncertain_species'] += 1
            
            if tech_result == 'SCRNA':
                self.stats['scrna_records'] += 1
            elif tech_result == 'NOT_SCRNA':
                self.stats['non_scrna_records'] += 1
            else:
                self.stats['uncertain_technology'] += 1
        
        # 添加检查结果到数据框
        results_df = pd.DataFrame(results)
        df = df.join(results_df.set_index('index'), how='left')
        
        # AI辅助检查不确定的记录
        if use_ai_for_uncertain and self.kimi and uncertain_indices:
            logger.info(f"使用AI辅助检查 {len(uncertain_indices)} 条不确定的记录...")
            self._ai_assisted_review(df, uncertain_indices)
        
        logger.info(f"质量检查完成: 人类数据={self.stats['human_records']}, "
                   f"非人类={self.stats['non_human_records']}, "
                   f"不确定={self.stats['uncertain_species']}")
        
        return df
    
    def _ai_assisted_review(self, df: pd.DataFrame, uncertain_indices: List[int]):
        """AI辅助审查不确定的记录"""
        batch_size = 10
        reviewed = 0
        
        for i in range(0, len(uncertain_indices), batch_size):
            batch_indices = uncertain_indices[i:i + batch_size]
            
            for idx in batch_indices:
                row = df.loc[idx]
                
                # 判断需要检查的类型
                check_type = 'both'
                if row['species_check'] != 'UNCERTAIN' and row['species_confidence'] >= 0.7:
                    check_type = 'technology'
                elif row['technology_check'] != 'UNCERTAIN' and row['technology_confidence'] >= 0.7:
                    check_type = 'species'
                
                # AI检查
                ai_result = self.quality_checker.quality_check_with_ai(row, check_type)
                
                if ai_result:
                    # 更新结果
                    if 'ai_species_check' in ai_result:
                        df.at[idx, 'ai_species_check'] = ai_result['ai_species_check']
                        df.at[idx, 'ai_species_reason'] = ai_result['ai_species_reason']
                    
                    if 'ai_technology_check' in ai_result:
                        df.at[idx, 'ai_technology_check'] = ai_result['ai_technology_check']
                        df.at[idx, 'ai_technology_reason'] = ai_result['ai_technology_reason']
                
                reviewed += 1
                if reviewed % 10 == 0:
                    logger.info(f"AI审查进度: {reviewed}/{len(uncertain_indices)}")
                
                time.sleep(1)  # 避免API限流
        
        logger.info(f"AI辅助审查完成: {reviewed} 条记录")
    
    def filter_and_flag(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        过滤和标记数据
        返回: (清洁数据, 需要人工复核的数据)
        """
        logger.info("开始过滤和标记数据...")
        
        # 创建最终判断列
        df['final_species'] = df.apply(self._determine_final_species, axis=1)
        df['final_technology'] = df.apply(self._determine_final_technology, axis=1)
        df['flag_for_review'] = False
        df['review_reason'] = ''
        
        # 标记需要复核的记录
        conditions_for_review = [
            (df['final_species'] == 'UNCERTAIN', 'Species uncertain'),
            (df['final_technology'] == 'UNCERTAIN', 'Technology uncertain'),
            ((df['species_confidence'] < 0.7) & (df['final_species'] == 'HUMAN'), 
             'Low confidence human classification'),
            ((df['technology_confidence'] < 0.7) & (df['final_technology'] == 'SCRNA'),
             'Low confidence scRNA-seq classification'),
            (df['final_species'] == 'NON_HUMAN', 'Non-human species'),
            (df['final_technology'] == 'NOT_SCRNA', 'Not scRNA-seq technology')
        ]
        
        for condition, reason in conditions_for_review:
            mask = condition & ~df['flag_for_review']
            df.loc[mask, 'flag_for_review'] = True
            df.loc[mask, 'review_reason'] = df.loc[mask, 'review_reason'].apply(
                lambda x: f"{x}; {reason}" if x else reason
            )
        
        # 分离清洁数据和需要复核的数据
        clean_df = df[~df['flag_for_review']].copy()
        review_df = df[df['flag_for_review']].copy()
        
        self.stats['final_clean_records'] = len(clean_df)
        self.stats['flagged_for_review'] = len(review_df)
        
        logger.info(f"过滤完成: 清洁数据={len(clean_df)}, 需要复核={len(review_df)}")
        
        return clean_df, review_df
    
    def _determine_final_species(self, row: pd.Series) -> str:
        """确定最终的物种判断"""
        # 优先使用AI判断
        if 'ai_species_check' in row and pd.notna(row['ai_species_check']):
            return row['ai_species_check']
        
        # 使用基础判断
        if row['species_confidence'] >= 0.7:
            return row['species_check']
        
        return 'UNCERTAIN'
    
    def _determine_final_technology(self, row: pd.Series) -> str:
        """确定最终的技术类型判断"""
        # 优先使用AI判断
        if 'ai_technology_check' in row and pd.notna(row['ai_technology_check']):
            return row['ai_technology_check']
        
        # 使用基础判断
        if row['technology_confidence'] >= 0.7:
            return row['technology_check']
        
        return 'UNCERTAIN'
    
    def reorganize_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """重组数据为新格式"""
        logger.info("开始重组数据...")
        
        reorganized_rows = []
        
        for idx, row in df.iterrows():
            new_row = self.reorganizer.reorganize_row(row)
            
            # 添加质量检查信息到supplementary_information
            try:
                supp_info = json.loads(new_row['supplementary_information']) if new_row['supplementary_information'] else {}
            except:
                supp_info = {}
            
            supp_info['quality_check'] = {
                'species_check': row.get('final_species', ''),
                'species_confidence': float(row.get('species_confidence', 0)),
                'technology_check': row.get('final_technology', ''),
                'technology_confidence': float(row.get('technology_confidence', 0))
            }
            
            new_row['supplementary_information'] = json.dumps(supp_info, ensure_ascii=False)
            
            reorganized_rows.append(new_row)
        
        reorganized_df = pd.DataFrame(reorganized_rows)
        
        # 确保列顺序正确
        column_order = [
            'Study/Project_id_1', 'Study/Project_id_2', 'Study/Project_id_3',
            'raw_sample_id', 'matrix_sample_id', 
            'raw_exist', 'raw_open', 'matrix_exist', 'matrix_open',
            'file_type', 'title', 'disease_general', 'disease', 'pubmed',
            'source_database', 'access_link', 'open_status',
            'ethnicity', 'sex', 'tissue_location', 'sequencing_platform',
            'experiment_design', 'sample_type', 'summary', 'citation_count',
            'publication_date', 'submission_date', 'last_update_date',
            'contact_name', 'contact_email', 'contact_institute',
            'supplementary_information'
        ]
        
        reorganized_df = reorganized_df[column_order]
        
        logger.info(f"数据重组完成: {len(reorganized_df)} 条记录")
        
        return reorganized_df
    
    def save_results(self, clean_df: pd.DataFrame, review_df: pd.DataFrame,
                    clean_reorganized: pd.DataFrame, review_reorganized: pd.DataFrame):
        """保存结果"""
        logger.info("保存结果...")
        
        # 保存清洁数据
        clean_reorganized.to_csv(
            self.output_dir / "cleaned_metadata.csv", 
            index=False, encoding='utf-8-sig'
        )
        
        # 保存需要复核的数据
        review_reorganized.to_csv(
            self.output_dir / "flagged_for_review.csv",
            index=False, encoding='utf-8-sig'
        )
        
        # 保存详细的质量检查报告
        review_report = review_df[[
            'id', 'title', 'source_database',
            'species_check', 'species_confidence', 'species_reason',
            'technology_check', 'technology_confidence', 'technology_reason',
            'ai_species_check', 'ai_species_reason',
            'ai_technology_check', 'ai_technology_reason',
            'flag_for_review', 'review_reason'
        ]].copy()
        
        review_report.to_csv(
            self.output_dir / "quality_check_report.csv",
            index=False, encoding='utf-8-sig'
        )
        
        # 保存统计报告
        stats_report = {
            'cleaning_date': datetime.now().isoformat(),
            'statistics': self.stats,
            'quality_distribution': {
                'species': {
                    'human': self.stats['human_records'],
                    'non_human': self.stats['non_human_records'],
                    'uncertain': self.stats['uncertain_species']
                },
                'technology': {
                    'scrna_seq': self.stats['scrna_records'],
                    'not_scrna_seq': self.stats['non_scrna_records'],
                    'uncertain': self.stats['uncertain_technology']
                }
            },
            'final_results': {
                'clean_records': self.stats['final_clean_records'],
                'flagged_for_review': self.stats['flagged_for_review'],
                'rejection_rate': f"{(self.stats['flagged_for_review'] / self.stats['total_records'] * 100):.2f}%"
            }
        }
        
        with open(self.output_dir / "cleaning_statistics.json", 'w', encoding='utf-8') as f:
            json.dump(stats_report, f, indent=2, ensure_ascii=False)
        
        # 生成文本报告
        self._generate_text_report(stats_report)
        
        logger.info("所有结果已保存")
    
    def _generate_text_report(self, stats_report: Dict):
        """生成文本格式的报告"""
        report_path = self.output_dir / "cleaning_report.txt"
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("人类单细胞RNA-seq元数据清洗报告\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"清洗时间: {stats_report['cleaning_date']}\n\n")
            
            f.write("一、总体统计\n")
            f.write("-" * 80 + "\n")
            f.write(f"总记录数: {self.stats['total_records']}\n")
            f.write(f"清洁记录数: {self.stats['final_clean_records']}\n")
            f.write(f"需要复核记录数: {self.stats['flagged_for_review']}\n")
            f.write(f"拒绝率: {stats_report['final_results']['rejection_rate']}\n\n")
            
            f.write("二、物种分类\n")
            f.write("-" * 80 + "\n")
            f.write(f"人类数据: {self.stats['human_records']} "
                   f"({self.stats['human_records']/self.stats['total_records']*100:.2f}%)\n")
            f.write(f"非人类数据: {self.stats['non_human_records']} "
                   f"({self.stats['non_human_records']/self.stats['total_records']*100:.2f}%)\n")
            f.write(f"不确定: {self.stats['uncertain_species']} "
                   f"({self.stats['uncertain_species']/self.stats['total_records']*100:.2f}%)\n\n")
            
            f.write("三、技术类型分类\n")
            f.write("-" * 80 + "\n")
            f.write(f"单细胞RNA-seq: {self.stats['scrna_records']} "
                   f"({self.stats['scrna_records']/self.stats['total_records']*100:.2f}%)\n")
            f.write(f"非单细胞RNA-seq: {self.stats['non_scrna_records']} "
                   f"({self.stats['non_scrna_records']/self.stats['total_records']*100:.2f}%)\n")
            f.write(f"不确定: {self.stats['uncertain_technology']} "
                   f"({self.stats['uncertain_technology']/self.stats['total_records']*100:.2f}%)\n\n")
            
            f.write("四、输出文件说明\n")
            f.write("-" * 80 + "\n")
            f.write("1. cleaned_metadata.csv - 清洁的元数据 (可直接使用)\n")
            f.write("2. flagged_for_review.csv - 需要人工复核的记录\n")
            f.write("3. quality_check_report.csv - 详细的质量检查报告\n")
            f.write("4. cleaning_statistics.json - 统计信息 (JSON格式)\n")
            f.write("5. cleaning_report.txt - 本报告\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("报告结束\n")
            f.write("=" * 80 + "\n")
    
    def run(self, use_ai_for_uncertain: bool = True):
        """运行完整的清洗流程"""
        logger.info("="*80)
        logger.info("开始元数据清洗流程")
        logger.info("="*80)
        
        # 1. 加载数据
        df = self.load_merged_data()
        
        # 2. 质量检查
        df = self.perform_quality_checks(df, use_ai_for_uncertain)
        
        # 3. 过滤和标记
        clean_df, review_df = self.filter_and_flag(df)
        
        # 4. 重组数据
        clean_reorganized = self.reorganize_data(clean_df)
        review_reorganized = self.reorganize_data(review_df)
        
        # 5. 保存结果
        self.save_results(clean_df, review_df, clean_reorganized, review_reorganized)
        
        logger.info("="*80)
        logger.info("元数据清洗完成!")
        logger.info("="*80)
        logger.info(f"清洁记录: {self.stats['final_clean_records']}")
        logger.info(f"需要复核: {self.stats['flagged_for_review']}")
        logger.info(f"结果保存在: {self.output_dir}")
        logger.info("="*80)
        
        return clean_reorganized, review_reorganized


def main():
    """主函数"""
    print("="*80)
    print("人类单细胞RNA-seq元数据清洗系统")
    print("="*80)
    print()
    
    # Kimi API keys
    kimi_api_keys = [
        "sk-R8WDgcyv5aQuW4OFKQg4v8YQZe0bTK90TbNXzYXFxa9ibE8P",
        "sk-zpkf8mQSlm0hVOtx6KEFtI5wun36NLJ8JVrglwvpEKgLThJW"
    ]
    
    # 创建清洗流程
    pipeline = MetadataCleaningPipeline(
        input_dir="metadata_output",
        output_dir="cleaned_metadata",
        kimi_api_keys=kimi_api_keys,
        use_ai_assistance=True  # 设置为True以使用AI辅助
    )
    
    # 运行清洗流程
    try:
        clean_data, review_data = pipeline.run(use_ai_for_uncertain=True)
        
        print("\n" + "="*80)
        print("清洗完成!")
        print("="*80)
        print(f"\n清洁数据: {len(clean_data)} 条记录")
        print(f"需要复核: {len(review_data)} 条记录")
        print(f"\n结果保存在: cleaned_metadata/")
        print("\n输出文件:")
        print("  - cleaned_metadata.csv (清洁的数据)")
        print("  - flagged_for_review.csv (需要复核的数据)")
        print("  - quality_check_report.csv (质量检查详情)")
        print("  - cleaning_statistics.json (统计信息)")
        print("  - cleaning_report.txt (文本报告)")
        print("="*80)
        
    except Exception as e:
        logger.error(f"清洗流程出错: {e}")
        import traceback
        logger.error(traceback.format_exc())


if __name__ == "__main__":
    main()