#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEO Human scRNA-seq Data Cleaning and Standardization Pipeline
数据清洗和标准化流程 - 移除非人源、非单细胞、非RNA-seq数据
Version 1.0
"""

import os
import re
import time
import pandas as pd
import requests
from datetime import datetime
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import anthropic
from typing import Dict, List, Tuple, Optional
import numpy as np

class KimiAPIClient:
    """Kimi API客户端 - 用于辅助判断复杂样本"""
    
    def __init__(self, api_keys: List[str]):
        self.api_keys = api_keys
        self.current_key_index = 0
        self.base_url = "https://api.moonshot.cn/v1/chat/completions"
        self.max_workers = min(len(api_keys), 5)
    
    def _get_next_key(self) -> str:
        """轮换使用API key"""
        key = self.api_keys[self.current_key_index]
        self.current_key_index = (self.current_key_index + 1) % len(self.api_keys)
        return key
    
    def check_data_quality(self, title: str, summary: str, platform: str, 
                          library_strategy: str, sample_info: str) -> Dict:
        """
        使用Kimi API检查数据质量
        
        Returns:
        --------
        dict: {
            'is_human': bool,
            'is_single_cell': bool,
            'is_rna_seq': bool,
            'confidence': float,
            'reason': str
        }
        """
        prompt = f"""请分析以下GEO数据集信息，判断它是否符合以下三个标准：
1. 是否为人类(Homo sapiens)数据
2. 是否为单细胞精度数据(single-cell resolution)
3. 是否为RNA测序数据(RNA-seq)

数据信息：
标题: {title}
摘要: {summary}
平台: {platform}
测序策略: {library_strategy}
样本信息: {sample_info}

请以JSON格式返回结果：
{{
    "is_human": true/false,
    "is_single_cell": true/false,
    "is_rna_seq": true/false,
    "confidence": 0.0-1.0,
    "reason": "详细原因说明"
}}

只返回JSON，不要其他内容。"""

        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self._get_next_key()}"
        }
        
        data = {
            "model": "moonshot-v1-8k",
            "messages": [
                {"role": "user", "content": prompt}
            ],
            "temperature": 0.3
        }
        
        try:
            response = requests.post(self.base_url, headers=headers, json=data, timeout=30)
            response.raise_for_status()
            result = response.json()
            
            content = result['choices'][0]['message']['content']
            # 提取JSON部分
            json_match = re.search(r'\{[^{}]*\}', content, re.DOTALL)
            if json_match:
                return json.loads(json_match.group())
            else:
                return {
                    'is_human': False,
                    'is_single_cell': False,
                    'is_rna_seq': False,
                    'confidence': 0.0,
                    'reason': 'API返回格式错误'
                }
        except Exception as e:
            print(f"  Kimi API调用失败: {e}")
            return {
                'is_human': None,
                'is_single_cell': None,
                'is_rna_seq': None,
                'confidence': 0.0,
                'reason': f'API错误: {str(e)}'
            }


class DataCleaner:
    """数据清洗器 - 识别和移除非符合条件的数据"""
    
    def __init__(self, series_file: str, samples_file: str, use_kimi: bool = True):
        """
        Parameters:
        -----------
        series_file : str
            Series metadata文件路径
        samples_file : str
            Samples metadata文件路径
        use_kimi : bool
            是否使用Kimi API辅助判断
        """
        self.series_df = pd.read_csv(series_file) if series_file.endswith('.csv') else pd.read_excel(series_file)
        self.samples_df = pd.read_csv(samples_file) if samples_file.endswith('.csv') else pd.read_excel(samples_file)
        
        # Kimi API客户端
        self.kimi_client = None
        if use_kimi:
            api_keys = [
                "sk-UL5YodR7ZL4S9dytpfMWJgmPTXJjkeNSd7Ktq9bbEhElzDfX",
                "sk-WL1hTKhW3sYKuJjBSc1k8wxL0r8ZLcuM7YiYxFjGgVHqAXhU",
                "sk-RAkA28HIT5tiMEfKtXAgbZ9nZKweq5Bnw0WbSwwBdNX7nbi1"
            ]
            self.kimi_client = KimiAPIClient(api_keys)
        
        # 人类标识符
        self.human_identifiers = [
            'homo sapiens', 'human', '9606', 'hg19', 'hg38', 'grch37', 'grch38'
        ]
        
        # 单细胞关键词
        self.single_cell_keywords = [
            'single cell', 'single-cell', 'singlecell', 'scRNA-seq', 'scRNA',
            '10x', 'drop-seq', 'smart-seq', 'cel-seq', 'scrna', 'sc-rna'
        ]
        
        # RNA-seq关键词
        self.rna_seq_keywords = [
            'rna-seq', 'rna seq', 'transcriptome', 'gene expression',
            'mrna-seq', 'mrna seq'
        ]
        
        # 排除关键词
        self.exclude_keywords = {
            'bulk': ['bulk rna', 'bulk rna-seq', 'bulk transcriptome'],
            'atac': ['atac-seq', 'atac seq', 'chromatin accessibility'],
            'chip': ['chip-seq', 'chip seq', 'chromatin immunoprecipitation'],
            'wgs': ['whole genome sequencing', 'wgs', 'genome sequencing'],
            'wes': ['whole exome', 'exome sequencing', 'wes'],
            'methylation': ['methylation', 'bisulfite', 'bs-seq'],
            'microarray': ['microarray', 'array'],
            'spatial': ['spatial transcriptom', 'visium', 'slide-seq']  # 某些spatial可能是单细胞，需要仔细判断
        }
        
        # 结果存储
        self.clean_data = []
        self.flagged_data = []
        self.removed_data = []
    
    def _safe_str(self, value) -> str:
        """安全字符串转换"""
        if pd.isna(value) or value is None:
            return ''
        return str(value)
    
    def _safe_lower(self, value) -> str:
        """安全小写转换"""
        return self._safe_str(value).lower()
    
    def _is_human(self, text: str) -> bool:
        """判断是否为人类数据"""
        text_lower = self._safe_lower(text)
        return any(identifier in text_lower for identifier in self.human_identifiers)
    
    def _is_single_cell(self, text: str, library_strategy: str = '') -> Tuple[bool, float]:
        """
        判断是否为单细胞数据
        
        Returns:
        --------
        tuple: (is_single_cell, confidence)
        """
        text_lower = self._safe_lower(text)
        strategy_lower = self._safe_lower(library_strategy)
        
        # 强排除bulk
        for bulk_kw in self.exclude_keywords['bulk']:
            if bulk_kw in text_lower:
                return False, 1.0
        
        # 检查单细胞关键词
        sc_score = sum(1 for kw in self.single_cell_keywords if kw in text_lower)
        strategy_score = sum(1 for kw in self.single_cell_keywords if kw in strategy_lower)
        
        total_score = sc_score + strategy_score
        
        if total_score >= 2:
            return True, 0.9
        elif total_score == 1:
            return True, 0.6
        else:
            return False, 0.3
    
    def _is_rna_seq(self, text: str, library_strategy: str = '') -> Tuple[bool, float]:
        """
        判断是否为RNA-seq数据
        
        Returns:
        --------
        tuple: (is_rna_seq, confidence)
        """
        text_lower = self._safe_lower(text)
        strategy_lower = self._safe_lower(library_strategy)
        
        # 检查排除类型
        for exclude_type, keywords in self.exclude_keywords.items():
            if exclude_type in ['atac', 'chip', 'wgs', 'wes', 'methylation']:
                for kw in keywords:
                    if kw in text_lower or kw in strategy_lower:
                        return False, 1.0
        
        # 检查RNA-seq关键词
        rna_score = sum(1 for kw in self.rna_seq_keywords if kw in text_lower)
        strategy_score = 1 if 'rna' in strategy_lower else 0
        
        total_score = rna_score + strategy_score
        
        if total_score >= 2:
            return True, 0.9
        elif total_score == 1:
            return True, 0.6
        else:
            return False, 0.3
    
    def check_sample_quality(self, series_id: str, sample_data: pd.Series) -> Dict:
        """
        检查单个样本的质量
        
        Returns:
        --------
        dict: {
            'series_id': str,
            'sample_id': str,
            'is_human': bool,
            'is_single_cell': bool,
            'is_rna_seq': bool,
            'human_confidence': float,
            'sc_confidence': float,
            'rna_confidence': float,
            'overall_confidence': float,
            'needs_review': bool,
            'reason': str
        }
        """
        # 获取Series信息
        series_info = self.series_df[self.series_df['Series_id'] == series_id]
        if len(series_info) == 0:
            return {
                'series_id': series_id,
                'sample_id': self._safe_str(sample_data.get('Sample_id', '')),
                'is_human': False,
                'is_single_cell': False,
                'is_rna_seq': False,
                'human_confidence': 0.0,
                'sc_confidence': 0.0,
                'rna_confidence': 0.0,
                'overall_confidence': 0.0,
                'needs_review': True,
                'reason': 'Series信息缺失'
            }
        
        series_info = series_info.iloc[0]
        
        # 组合文本信息
        title = self._safe_str(series_info.get('Title', ''))
        summary = self._safe_str(series_info.get('Summary', ''))
        platform = self._safe_str(series_info.get('Platform', ''))
        sample_title = self._safe_str(sample_data.get('Title', ''))
        organism = self._safe_str(sample_data.get('Organism', ''))
        library_strategy = self._safe_str(sample_data.get('Library_strategy', ''))
        library_source = self._safe_str(sample_data.get('Library_source', ''))
        
        full_text = f"{title} {summary} {sample_title} {organism} {platform}"
        
        # 1. 检查人类
        is_human = self._is_human(full_text)
        human_confidence = 0.9 if is_human else 0.1
        
        # 2. 检查单细胞
        is_single_cell, sc_confidence = self._is_single_cell(full_text, library_strategy)
        
        # 3. 检查RNA-seq
        is_rna_seq, rna_confidence = self._is_rna_seq(full_text, library_strategy)
        
        # 计算总体置信度
        overall_confidence = (human_confidence + sc_confidence + rna_confidence) / 3.0
        
        # 判断是否需要人工复核
        needs_review = False
        reason_parts = []
        
        if not is_human:
            reason_parts.append(f"非人类数据(organism: {organism})")
        
        if not is_single_cell:
            reason_parts.append("非单细胞数据")
        elif sc_confidence < 0.7:
            needs_review = True
            reason_parts.append(f"单细胞置信度低({sc_confidence:.2f})")
        
        if not is_rna_seq:
            reason_parts.append("非RNA-seq数据")
        elif rna_confidence < 0.7:
            needs_review = True
            reason_parts.append(f"RNA-seq置信度低({rna_confidence:.2f})")
        
        # 特殊情况标记
        text_lower = full_text.lower()
        if 'spatial' in text_lower:
            needs_review = True
            reason_parts.append("包含spatial关键词，需确认是否为单细胞")
        
        if 'multiome' in text_lower or 'cite-seq' in text_lower:
            needs_review = True
            reason_parts.append("多组学数据，需确认RNA部分")
        
        reason = '; '.join(reason_parts) if reason_parts else 'OK'
        
        result = {
            'series_id': series_id,
            'sample_id': self._safe_str(sample_data.get('Sample_id', '')),
            'is_human': is_human,
            'is_single_cell': is_single_cell,
            'is_rna_seq': is_rna_seq,
            'human_confidence': human_confidence,
            'sc_confidence': sc_confidence,
            'rna_confidence': rna_confidence,
            'overall_confidence': overall_confidence,
            'needs_review': needs_review,
            'reason': reason
        }
        
        # 如果需要复核且启用了Kimi API，使用API进一步判断
        if needs_review and self.kimi_client and overall_confidence > 0.3:
            print(f"  使用Kimi API复核: {series_id} - {result['sample_id']}")
            sample_info_text = f"Title: {sample_title}, Organism: {organism}, Library: {library_strategy}"
            
            kimi_result = self.kimi_client.check_data_quality(
                title, summary, platform, library_strategy, sample_info_text
            )
            
            if kimi_result['confidence'] > 0.7:
                result['is_human'] = kimi_result['is_human']
                result['is_single_cell'] = kimi_result['is_single_cell']
                result['is_rna_seq'] = kimi_result['is_rna_seq']
                result['overall_confidence'] = kimi_result['confidence']
                result['reason'] += f" | Kimi: {kimi_result['reason']}"
        
        return result
    
    def clean_all_data(self) -> Tuple[List[Dict], List[Dict], List[Dict]]:
        """
        清洗所有数据
        
        Returns:
        --------
        tuple: (clean_data, flagged_data, removed_data)
        """
        print("\n开始数据清洗...")
        print(f"总Series数: {len(self.series_df)}")
        print(f"总Samples数: {len(self.samples_df)}")
        
        total_samples = len(self.samples_df)
        
        for idx, (_, sample) in enumerate(self.samples_df.iterrows(), 1):
            series_id = self._safe_str(sample.get('Series_id', ''))
            
            if idx % 10 == 0 or idx == total_samples:
                print(f"\r处理进度: {idx}/{total_samples} ({idx/total_samples*100:.1f}%)", end='')
            
            quality_check = self.check_sample_quality(series_id, sample)
            
            # 分类
            if quality_check['is_human'] and quality_check['is_single_cell'] and quality_check['is_rna_seq']:
                if quality_check['needs_review']:
                    self.flagged_data.append({**sample.to_dict(), **quality_check})
                else:
                    self.clean_data.append({**sample.to_dict(), **quality_check})
            else:
                self.removed_data.append({**sample.to_dict(), **quality_check})
        
        print(f"\n\n清洗完成:")
        print(f"  ✓ 清洁数据: {len(self.clean_data)} 个样本")
        print(f"  ⚠ 待复核数据: {len(self.flagged_data)} 个样本")
        print(f"  ✗ 移除数据: {len(self.removed_data)} 个样本")
        
        return self.clean_data, self.flagged_data, self.removed_data
    
    def save_cleaning_report(self, output_dir: str = "cleaning_report"):
        """保存清洗报告"""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # 保存清洁数据
        if self.clean_data:
            clean_df = pd.DataFrame(self.clean_data)
            clean_file = os.path.join(output_dir, f'clean_data_{timestamp}.csv')
            clean_df.to_csv(clean_file, index=False, encoding='utf-8-sig')
            print(f"\n清洁数据已保存: {clean_file}")
        
        # 保存待复核数据
        if self.flagged_data:
            flagged_df = pd.DataFrame(self.flagged_data)
            flagged_file = os.path.join(output_dir, f'flagged_for_review_{timestamp}.csv')
            flagged_df.to_csv(flagged_file, index=False, encoding='utf-8-sig')
            print(f"待复核数据已保存: {flagged_file}")
        
        # 保存移除数据
        if self.removed_data:
            removed_df = pd.DataFrame(self.removed_data)
            removed_file = os.path.join(output_dir, f'removed_data_{timestamp}.csv')
            removed_df.to_csv(removed_file, index=False, encoding='utf-8-sig')
            print(f"移除数据已保存: {removed_file}")
        
        # 保存统计摘要
        summary = {
            'total_samples': len(self.clean_data) + len(self.flagged_data) + len(self.removed_data),
            'clean_samples': len(self.clean_data),
            'flagged_samples': len(self.flagged_data),
            'removed_samples': len(self.removed_data),
            'clean_rate': len(self.clean_data) / (len(self.clean_data) + len(self.flagged_data) + len(self.removed_data)) * 100
        }
        
        summary_file = os.path.join(output_dir, f'cleaning_summary_{timestamp}.json')
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        print(f"清洗摘要已保存: {summary_file}")
        
        return clean_file if self.clean_data else None


class DataStandardizer:
    """数据标准化器 - 转换为目标格式"""
    
    def __init__(self, series_file: str, samples_file: str, clean_samples_file: Optional[str] = None):
        """
        Parameters:
        -----------
        series_file : str
            Series metadata文件
        samples_file : str
            Samples metadata文件
        clean_samples_file : str, optional
            清洗后的样本文件（如果有）
        """
        self.series_df = pd.read_csv(series_file) if series_file.endswith('.csv') else pd.read_excel(series_file)
        
        # 如果提供了清洗后的文件，使用它
        if clean_samples_file:
            self.samples_df = pd.read_csv(clean_samples_file)
        else:
            self.samples_df = pd.read_csv(samples_file) if samples_file.endswith('.csv') else pd.read_excel(samples_file)
        
        # 目标字段
        self.target_fields = [
            'Study/Project_id_1', 'Study/Project_id_2', 'Study/Project_id_3',
            'raw_sample_id', 'matrix_sample_id', 'raw_exist', 'raw_open',
            'matrix_exist', 'matrix_open', 'file_type', 'title',
            'disease_general', 'disease', 'pubmed', 'source_database',
            'access_link', 'open_status', 'ethnicity', 'sex',
            'tissue_location', 'sequencing_platform', 'experiment_design',
            'sample_type', 'summary', 'citation_count', 'publication_date',
            'submission_date', 'last_update_date', 'contact_name',
            'contact_email', 'contact_institute', 'supplementary_information'
        ]
    
    def _safe_str(self, value) -> str:
        """安全字符串转换"""
        if pd.isna(value) or value is None:
            return ''
        return str(value)
    
    def _safe_lower(self, value) -> str:
        """安全小写转换"""
        return self._safe_str(value).lower()
    
    def _extract_sra_id(self, sample_data: pd.Series) -> str:
        """提取SRA ID"""
        # 从supplementary files中查找SRA相关ID
        supp_files = self._safe_str(sample_data.get('Supplementary_files', ''))
        
        # 查找SRR, SRX等ID
        sra_pattern = r'(SRR\d+|SRX\d+|SRP\d+)'
        match = re.search(sra_pattern, supp_files)
        if match:
            return match.group(1)
        
        return ''
    
    def _determine_file_type(self, sample_data: pd.Series) -> str:
        """判断文件类型"""
        supp_files = self._safe_lower(sample_data.get('Supplementary_files', ''))
        
        file_types = []
        if '.fastq' in supp_files or '.fq' in supp_files:
            file_types.append('FASTQ')
        if '.bam' in supp_files or '.sam' in supp_files:
            file_types.append('BAM')
        if '.h5' in supp_files or '.h5ad' in supp_files:
            file_types.append('H5')
        if '.mtx' in supp_files:
            file_types.append('MTX')
        if '.rds' in supp_files:
            file_types.append('RDS')
        if '.csv' in supp_files or '.tsv' in supp_files or '.txt' in supp_files:
            file_types.append('TXT/CSV')
        
        return '; '.join(file_types) if file_types else 'Unknown'
    
    def _check_data_availability(self, sample_data: pd.Series, series_data: pd.Series) -> Tuple[str, str, str, str]:
        """
        检查数据可用性
        
        Returns:
        --------
        tuple: (raw_exist, raw_open, matrix_exist, matrix_open)
        """
        raw_files = self._safe_str(sample_data.get('Raw_data_files', '[]'))
        processed_files = self._safe_str(sample_data.get('Processed_data_files', '[]'))
        ftp_link = self._safe_str(series_data.get('FTP_Link', ''))
        sra_link = self._safe_str(series_data.get('SRA_Link', ''))
        
        # 原始数据
        raw_exist = 'Y' if (raw_files and raw_files != '[]') or sra_link else 'N'
        raw_open = 'Y' if raw_exist == 'Y' else 'N'  # GEO数据默认公开
        
        # 矩阵数据
        matrix_exist = 'Y' if (processed_files and processed_files != '[]') or ftp_link else 'N'
        matrix_open = 'Y' if matrix_exist == 'Y' else 'N'
        
        return raw_exist, raw_open, matrix_exist, matrix_open
    
    def _extract_disease_info(self, series_data: pd.Series, sample_data: pd.Series) -> Tuple[str, str]:
        """
        提取疾病信息
        
        Returns:
        --------
        tuple: (disease_general, disease_specific)
        """
        title = self._safe_lower(series_data.get('Title', ''))
        summary = self._safe_lower(series_data.get('Summary', ''))
        sample_title = self._safe_lower(sample_data.get('Title', ''))
        
        full_text = f"{title} {summary} {sample_title}"
        
        # 疾病映射
        disease_mapping = {
            'breast cancer': ('BC', 'BRCA'),
            'breast carcinoma': ('BC', 'BRCA'),
            'lung cancer': ('LC', 'LUAD'),
            'lung adenocarcinoma': ('LC', 'LUAD'),
            'lung squamous': ('LC', 'LUSC'),
            'colorectal cancer': ('CRC', 'COAD'),
            'colon cancer': ('CRC', 'COAD'),
            'rectal cancer': ('CRC', 'READ'),
            'gastric cancer': ('GC', 'STAD'),
            'stomach cancer': ('GC', 'STAD'),
            'liver cancer': ('HCC', 'LIHC'),
            'hepatocellular carcinoma': ('HCC', 'LIHC'),
            'pancreatic cancer': ('PDAC', 'PAAD'),
            'pancreatic adenocarcinoma': ('PDAC', 'PAAD'),
            'prostate cancer': ('PRAD', 'PRAD'),
            'melanoma': ('MEL', 'SKCM'),
            'glioblastoma': ('GBM', 'GBM'),
            'glioma': ('GBM', 'LGG'),
            'ovarian cancer': ('OV', 'OV'),
            'renal': ('RCC', 'KIRC'),
            'kidney': ('RCC', 'KIRC'),
        }
        
        for disease_name, (general, specific) in disease_mapping.items():
            if disease_name in full_text:
                return general, specific
        
        # 如果没有找到特定癌症，检查是否为健康/正常样本
        if any(term in full_text for term in ['healthy', 'normal', 'control', 'non-cancer']):
            return 'Healthy', 'Normal'
        
        return '', ''
    
    def _extract_sample_type(self, sample_data: pd.Series) -> str:
        """
        提取样本类型
        
        Returns:
        --------
        str: 'Tumor', 'Healthy', 'Other'
        """
        characteristics = self._parse_characteristics(sample_data.get('Characteristics', '{}'))
        title = self._safe_lower(sample_data.get('Title', ''))
        source = self._safe_lower(sample_data.get('Source_name', ''))
        
        combined = f"{title} {source}"
        
        # 检查characteristics
        for key, value in characteristics.items():
            value_lower = self._safe_lower(value)
            if any(term in value_lower for term in ['tumor', 'cancer', 'malignant', 'carcinoma']):
                return 'Tumor'
            if any(term in value_lower for term in ['normal', 'healthy', 'control']):
                return 'Healthy'
        
        # 检查标题和来源
        if any(term in combined for term in ['tumor', 'cancer', 'malignant']):
            return 'Tumor'
        if any(term in combined for term in ['normal', 'healthy', 'control']):
            return 'Healthy'
        
        return 'Other'
    
    def _parse_characteristics(self, characteristics_str) -> Dict:
        """解析characteristics"""
        if pd.isna(characteristics_str) or characteristics_str == '{}':
            return {}
        
        try:
            if isinstance(characteristics_str, str):
                return eval(characteristics_str)
            elif isinstance(characteristics_str, dict):
                return characteristics_str
            return {}
        except:
            return {}
    
    def _extract_sex(self, sample_data: pd.Series) -> str:
        """提取性别"""
        characteristics = self._parse_characteristics(sample_data.get('Characteristics', '{}'))
        
        gender_keys = ['gender', 'sex']
        for key in gender_keys:
            if key in characteristics:
                gender = self._safe_lower(characteristics[key])
                if gender in ['male', 'm']:
                    return 'Male'
                elif gender in ['female', 'f']:
                    return 'Female'
                return characteristics[key]
        
        return ''
    
    def _extract_ethnicity(self, series_data: pd.Series) -> str:
        """提取种族/国籍"""
        institute = self._safe_str(series_data.get('Contact_Institute', ''))
        
        country_keywords = {
            'USA': ['USA', 'United States', 'America', 'American'],
            'China': ['China', 'Chinese', 'Beijing', 'Shanghai'],
            'UK': ['United Kingdom', 'UK', 'England', 'British'],
            'Germany': ['Germany', 'German'],
            'France': ['France', 'French'],
            'Japan': ['Japan', 'Japanese'],
            'Korea': ['Korea', 'Korean'],
            'Canada': ['Canada', 'Canadian'],
            'Australia': ['Australia', 'Australian'],
        }
        
        for country, keywords in country_keywords.items():
            if any(kw in institute for kw in keywords):
                return country
        
        return ''
    
    def _extract_tissue_location(self, sample_data: pd.Series) -> str:
        """提取组织位置"""
        characteristics = self._parse_characteristics(sample_data.get('Characteristics', '{}'))
        
        tissue_keys = ['tissue', 'tissue_type', 'organ', 'anatomic_site', 'location']
        
        for key in tissue_keys:
            if key in characteristics:
                return self._safe_str(characteristics[key])
        
        source = self._safe_str(sample_data.get('Source_name', ''))
        if source:
            return source
        
        return ''
    
    def _extract_platform(self, series_data: pd.Series) -> str:
        """提取测序平台"""
        platform = series_data.get('Platform', [])
        
        if isinstance(platform, str):
            try:
                platform = eval(platform)
            except:
                platform = [platform] if platform else []
        
        if isinstance(platform, list):
            return '; '.join(platform)
        
        return self._safe_str(platform)
    
    def _get_experiment_design(self, series_data: pd.Series, sample_data: pd.Series) -> str:
        """获取实验设计"""
        overall_design = self._safe_str(series_data.get('Overall_Design', ''))
        library_strategy = self._safe_str(sample_data.get('Library_strategy', ''))
        instrument = self._safe_str(sample_data.get('Instrument_model', ''))
        
        design_parts = []
        if library_strategy:
            design_parts.append(library_strategy)
        if instrument:
            design_parts.append(instrument)
        if overall_design and len(design_parts) == 0:
            design_parts.append(overall_design[:200])  # 限制长度
        
        return '; '.join(design_parts)
    
    def standardize_data(self) -> pd.DataFrame:
        """
        将数据标准化为目标格式
        
        Returns:
        --------
        pd.DataFrame: 标准化后的数据
        """
        print("\n开始数据标准化...")
        
        standardized_data = []
        total_samples = len(self.samples_df)
        
        for idx, (_, sample) in enumerate(self.samples_df.iterrows(), 1):
            series_id = self._safe_str(sample.get('Series_id', ''))
            
            if idx % 10 == 0 or idx == total_samples:
                print(f"\r处理进度: {idx}/{total_samples} ({idx/total_samples*100:.1f}%)", end='')
            
            # 获取Series信息
            series_info = self.series_df[self.series_df['Series_id'] == series_id]
            if len(series_info) == 0:
                continue
            series_info = series_info.iloc[0]
            
            # 创建标准化记录
            record = {}
            
            # Study/Project IDs
            record['Study/Project_id_1'] = series_id  # GEO ID
            record['Study/Project_id_2'] = self._extract_sra_id(sample)  # SRA ID if exists
            record['Study/Project_id_3'] = ''  # 预留给其他数据库ID
            
            # Sample IDs
            sample_id = self._safe_str(sample.get('Sample_id', ''))
            record['raw_sample_id'] = sample_id
            record['matrix_sample_id'] = sample_id  # 通常相同，除非特别注明
            
            # 数据可用性
            raw_exist, raw_open, matrix_exist, matrix_open = self._check_data_availability(sample, series_info)
            record['raw_exist'] = raw_exist
            record['raw_open'] = raw_open
            record['matrix_exist'] = matrix_exist
            record['matrix_open'] = matrix_open
            
            # 文件类型
            record['file_type'] = self._determine_file_type(sample)
            
            # 标题
            record['title'] = self._safe_str(series_info.get('Title', ''))
            
            # 疾病信息
            disease_general, disease_specific = self._extract_disease_info(series_info, sample)
            record['disease_general'] = disease_general
            record['disease'] = disease_specific
            
            # 出版信息
            record['pubmed'] = self._safe_str(series_info.get('PubMed_ID', ''))
            record['source_database'] = 'GEO'
            record['access_link'] = self._safe_str(series_info.get('GEO_Link', ''))
            record['open_status'] = 'Public'  # GEO数据默认公开
            
            # 样本特征
            record['ethnicity'] = self._extract_ethnicity(series_info)
            record['sex'] = self._extract_sex(sample)
            record['tissue_location'] = self._extract_tissue_location(sample)
            record['sequencing_platform'] = self._extract_platform(series_info)
            record['experiment_design'] = self._get_experiment_design(series_info, sample)
            record['sample_type'] = self._extract_sample_type(sample)
            
            # 描述信息
            record['summary'] = self._safe_str(series_info.get('Summary', ''))
            record['citation_count'] = ''  # 需要从外部API获取，这里留空
            
            # 日期信息
            record['publication_date'] = self._safe_str(series_info.get('Publication_Date', ''))
            record['submission_date'] = self._safe_str(series_info.get('Submission_Date', ''))
            record['last_update_date'] = self._safe_str(series_info.get('Last_Update_Date', ''))
            
            # 联系信息
            record['contact_name'] = self._safe_str(series_info.get('Contact_Name', ''))
            record['contact_email'] = self._safe_str(series_info.get('Contact_Email', ''))
            record['contact_institute'] = self._safe_str(series_info.get('Contact_Institute', ''))
            
            # 补充信息
            record['supplementary_information'] = self._safe_str(series_info.get('Overall_Design', ''))
            
            standardized_data.append(record)
        
        print(f"\n\n标准化完成: {len(standardized_data)} 条记录")
        
        # 创建DataFrame
        df = pd.DataFrame(standardized_data, columns=self.target_fields)
        
        return df
    
    def save_standardized_data(self, df: pd.DataFrame, output_dir: str = "standardized_data"):
        """保存标准化数据"""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # 保存为Excel和CSV
        excel_file = os.path.join(output_dir, f'standardized_metadata_{timestamp}.xlsx')
        csv_file = os.path.join(output_dir, f'standardized_metadata_{timestamp}.csv')
        
        df.to_excel(excel_file, index=False, engine='openpyxl')
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        
        print(f"\n标准化数据已保存:")
        print(f"  Excel: {excel_file}")
        print(f"  CSV: {csv_file}")
        
        return excel_file, csv_file


def main():
    """主函数 - 完整的清洗和标准化流程"""
    print("=" * 80)
    print("GEO Human scRNA-seq Data Cleaning and Standardization Pipeline")
    print("数据清洗和标准化流程")
    print("=" * 80)
    
    # ===== 输入文件 =====
    # 请修改为您实际的文件路径
    series_file = "/mnt/public7/pancancercol/meta_analysis/geo_metadata/comprehensive_metadata/complete_human_series_20251125_034124.csv"
    samples_file = "/mnt/public7/pancancercol/meta_analysis/geo_metadata/comprehensive_metadata/complete_human_samples_20251125_034124.csv"
    
    # 检查文件是否存在
    if not os.path.exists(series_file) or not os.path.exists(samples_file):
        print("\n错误: 输入文件不存在，请检查文件路径")
        print(f"Series文件: {series_file}")
        print(f"Samples文件: {samples_file}")
        return
    
    # ===== 第一步：数据清洗 =====
    print("\n" + "=" * 80)
    print("第一步: 数据质量检查和清洗")
    print("=" * 80)
    
    # 创建清洗器（use_kimi=True 启用API辅助判断）
    cleaner = DataCleaner(series_file, samples_file, use_kimi=True)
    
    # 执行清洗
    clean_data, flagged_data, removed_data = cleaner.clean_all_data()
    
    # 保存清洗报告
    clean_file = cleaner.save_cleaning_report()
    
    # ===== 第二步：数据标准化 =====
    print("\n" + "=" * 80)
    print("第二步: 数据标准化")
    print("=" * 80)
    
    if clean_file:
        # 使用清洗后的数据进行标准化
        standardizer = DataStandardizer(series_file, samples_file, clean_file)
    else:
        print("\n警告: 没有清洁数据，使用原始样本数据")
        standardizer = DataStandardizer(series_file, samples_file)
    
    # 执行标准化
    standardized_df = standardizer.standardize_data()
    
    # 保存标准化数据
    excel_file, csv_file = standardizer.save_standardized_data(standardized_df)
    
    # ===== 数据质量报告 =====
    print("\n" + "=" * 80)
    print("数据质量最终报告")
    print("=" * 80)
    
    print(f"\n总体统计:")
    print(f"  原始样本数: {len(cleaner.clean_data) + len(cleaner.flagged_data) + len(cleaner.removed_data)}")
    print(f"  清洁样本数: {len(cleaner.clean_data)}")
    print(f"  待复核样本数: {len(cleaner.flagged_data)}")
    print(f"  移除样本数: {len(cleaner.removed_data)}")
    print(f"  标准化记录数: {len(standardized_df)}")
    
    if len(standardized_df) > 0:
        print(f"\n字段完整性:")
        for col in standardized_df.columns:
            non_empty = (standardized_df[col].notna() & (standardized_df[col] != '')).sum()
            percentage = (non_empty / len(standardized_df)) * 100
            print(f"  {col}: {non_empty}/{len(standardized_df)} ({percentage:.1f}%)")
    
    print("\n" + "=" * 80)
    print("所有步骤完成!")
    print("=" * 80)
    print("\n输出文件:")
    print(f"  1. 清洁数据: cleaning_report/clean_data_*.csv")
    print(f"  2. 待复核数据: cleaning_report/flagged_for_review_*.csv")
    print(f"  3. 移除数据: cleaning_report/removed_data_*.csv")
    print(f"  4. 标准化数据: {excel_file}")
    print("\n请手工复核 'flagged_for_review' 文件中的数据!")


if __name__ == "__main__":
    main()