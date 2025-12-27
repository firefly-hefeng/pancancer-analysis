#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEO Data Smart Re-verification Module
智能复核模块 - 优先级分级，减少API调用
Version 2.0 - Cost Optimized
"""

import os
import time
import pandas as pd
import json
from datetime import datetime
from typing import Dict, List, Tuple
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed  # 添加这一行
class DataReVerifier:
    """数据复核器 - 使用Kimi API对可疑数据进行二次确认"""
    
    def __init__(self, kimi_client, series_df: pd.DataFrame, 
                 confidence_threshold_high: float = 0.85,
                 confidence_threshold_low: float = 0.4):
        """
        Parameters:
        -----------
        kimi_client : KimiAPIClient
            Kimi API客户端
        series_df : pd.DataFrame
            Series元数据
        confidence_threshold_high : float
            高置信度阈值 - 超过此值直接判定
        confidence_threshold_low : float
            低置信度阈值 - 低于此值直接排除
        """
        self.kimi_client = kimi_client
        self.series_df = series_df
        self.confidence_threshold_high = confidence_threshold_high
        self.confidence_threshold_low = confidence_threshold_low
        
        # 结果存储
        self.confirmed_clean = []      # 确认为清洁数据
        self.confirmed_removed = []    # 确认应移除
        self.need_manual_review = []   # 需要人工复审
    
    def _safe_str(self, value) -> str:
        """安全字符串转换"""
        if pd.isna(value) or value is None:
            return ''
        return str(value)
    
    def _create_detailed_prompt(self, sample_data: Dict, series_data: pd.Series) -> str:
        """创建详细的复核提示词"""
        
        # 提取关键信息
        title = self._safe_str(series_data.get('Title', ''))
        summary = self._safe_str(series_data.get('Summary', ''))
        overall_design = self._safe_str(series_data.get('Overall_Design', ''))
        platform = self._safe_str(series_data.get('Platform', ''))
        
        sample_id = self._safe_str(sample_data.get('Sample_id', ''))
        sample_title = self._safe_str(sample_data.get('Title', ''))
        organism = self._safe_str(sample_data.get('Organism', ''))
        source = self._safe_str(sample_data.get('Source_name', ''))
        library_strategy = self._safe_str(sample_data.get('Library_strategy', ''))
        library_source = self._safe_str(sample_data.get('Library_source', ''))
        instrument = self._safe_str(sample_data.get('Instrument_model', ''))
        
        # 提取characteristics
        characteristics = self._safe_str(sample_data.get('Characteristics', '{}'))
        
        prompt = f"""你是一个专业的生物信息学数据质量审核专家。请仔细分析以下GEO数据集的样本，判断它是否同时满足以下三个严格标准：

标准1: 必须是人类(Homo sapiens)数据
标准2: 必须是单细胞分辨率(single-cell resolution)的数据
标准3: 必须是RNA测序(RNA-seq)数据

【研究项目信息】
标题: {title}
摘要: {summary[:800]}
总体设计: {overall_design[:400]}
测序平台: {platform}

【样本详细信息】
样本ID: {sample_id}
样本标题: {sample_title}
物种: {organism}
样本来源: {source}
测序策略: {library_strategy}
文库来源: {library_source}
测序仪器: {instrument}
样本特征: {characteristics[:400]}

【重要判断依据】
1. 人类数据判断:
   - 物种必须明确为 "Homo sapiens" 或 "human"
   - 如果是其他物种(小鼠、大鼠等)则不符合
   
2. 单细胞判断:
   - 必须包含 "single cell", "single-cell", "scRNA-seq" 等关键词
   - 如果是 "bulk RNA-seq" 则不符合
   - 如果只是 "spatial transcriptomics" 需要谨慎判断
   
3. RNA-seq判断:
   - 必须是转录组测序，测序策略应包含 "RNA-Seq" 或相关术语
   - 如果是 ATAC-seq, ChIP-seq, WGS, WES 等则不符合
   - 如果是 multiome 数据，需要确认包含RNA部分

【输出要求】
请以严格的JSON格式返回结果（只返回JSON，不要任何其他文字）：

{{
    "is_human": true/false,
    "is_single_cell": true/false, 
    "is_rna_seq": true/false,
    "confidence": 0.0-1.0,
    "detailed_analysis": {{
        "organism_check": "对物种的详细分析",
        "cell_type_check": "对单细胞/bulk的详细分析",
        "seq_type_check": "对测序类型的详细分析",
        "key_evidence": "支持判断的关键证据",
        "concerns": "需要注意的疑点（如果有）"
    }},
    "final_decision": "ACCEPT/REJECT/UNCERTAIN",
    "reason": "综合判断的简要说明"
}}

注意：
- confidence是你对判断的整体置信度(0.0-1.0)
- final_decision中: ACCEPT表示符合所有三个标准，REJECT表示明确不符合，UNCERTAIN表示无法确定
- 如果三个标准有任何一个不满足，final_decision应为REJECT
"""
        
        return prompt
    
    def reverify_sample(self, sample_data: Dict, max_retries: int = 3) -> Dict:
        """
        复核单个样本
        
        Returns:
        --------
        dict: 复核结果
        """
        series_id = self._safe_str(sample_data.get('Series_id', ''))
        
        # 获取Series信息
        series_info = self.series_df[self.series_df['Series_id'] == series_id]
        if len(series_info) == 0:
            return {
                'sample_id': self._safe_str(sample_data.get('Sample_id', '')),
                'series_id': series_id,
                'reverify_result': 'ERROR',
                'confidence': 0.0,
                'reason': 'Series信息缺失',
                'detailed_analysis': {}
            }
        
        series_info = series_info.iloc[0]
        
        # 创建详细提示词
        prompt = self._create_detailed_prompt(sample_data, series_info)
        
        for attempt in range(max_retries):
            try:
                # 调用API
                self.kimi_client._rate_limit()
                
                headers = {
                    "Content-Type": "application/json",
                    "Authorization": f"Bearer {self.kimi_client._get_next_key()}"
                }
                
                data = {
                    "model": "moonshot-v1-32k",  # 使用更大的模型以处理详细信息
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.2  # 降低温度以获得更确定的结果
                }
                
                import requests
                response = requests.post(
                    self.kimi_client.base_url, 
                    headers=headers, 
                    json=data, 
                    timeout=45
                )
                response.raise_for_status()
                result = response.json()
                
                content = result['choices'][0]['message']['content']
                
                # 提取JSON
                import re
                json_match = re.search(r'\{.*\}', content, re.DOTALL)
                if json_match:
                    parsed_result = json.loads(json_match.group())
                    
                    # 添加样本标识
                    parsed_result['sample_id'] = self._safe_str(sample_data.get('Sample_id', ''))
                    parsed_result['series_id'] = series_id
                    
                    return parsed_result
                else:
                    if attempt == max_retries - 1:
                        return self._get_error_result(sample_data, 'API返回格式错误')
                    time.sleep(1)
                    
            except requests.exceptions.Timeout:
                if attempt == max_retries - 1:
                    return self._get_error_result(sample_data, 'API请求超时')
                time.sleep(2)
                
            except Exception as e:
                if attempt == max_retries - 1:
                    return self._get_error_result(sample_data, f'错误: {str(e)}')
                time.sleep(1)
        
        return self._get_error_result(sample_data, '达到最大重试次数')
    
    def _get_error_result(self, sample_data: Dict, reason: str) -> Dict:
        """返回错误结果"""
        return {
            'sample_id': self._safe_str(sample_data.get('Sample_id', '')),
            'series_id': self._safe_str(sample_data.get('Series_id', '')),
            'reverify_result': 'ERROR',
            'confidence': 0.0,
            'reason': reason,
            'detailed_analysis': {}
        }
    
    def batch_reverify(self, samples_data: List[Dict], batch_size: int = 15) -> List[Dict]:
        """
        批量复核样本
        
        Parameters:
        -----------
        samples_data : List[Dict]
            待复核的样本数据列表
        batch_size : int
            批处理大小
        
        Returns:
        --------
        List[Dict]: 复核结果列表
        """
        print(f"\n开始批量复核 {len(samples_data)} 个样本...")
        
        results = []
        total = len(samples_data)
        
        # 分批处理
        for i in range(0, total, batch_size):
            batch = samples_data[i:i+batch_size]
            batch_num = i // batch_size + 1
            total_batches = (total + batch_size - 1) // batch_size
            
            print(f"\n处理批次 {batch_num}/{total_batches} ({len(batch)} 个样本)")
            
            # 使用线程池并发处理
            with ThreadPoolExecutor(max_workers=self.kimi_client.max_workers) as executor:
                futures = {
                    executor.submit(self.reverify_sample, sample): idx
                    for idx, sample in enumerate(batch)
                }
                
                batch_results = [None] * len(batch)
                completed = 0
                
                for future in as_completed(futures):
                    idx = futures[future]
                    try:
                        batch_results[idx] = future.result()
                        completed += 1
                        print(f"\r  进度: {completed}/{len(batch)}", end='')
                    except Exception as e:
                        batch_results[idx] = self._get_error_result(
                            batch[idx], 
                            f'并发处理错误: {str(e)}'
                        )
                
                results.extend(batch_results)
            
            # 批次间等待
            if i + batch_size < total:
                print(f"\n  等待 2 秒后继续...")
                time.sleep(2)
        
        print(f"\n\n批量复核完成!")
        return results
    
    def classify_reverify_results(self, reverify_results: List[Dict], 
                                  original_samples: List[Dict]) -> Tuple[List[Dict], List[Dict], List[Dict]]:
        """
        根据复核结果分类样本
        
        Parameters:
        -----------
        reverify_results : List[Dict]
            复核结果
        original_samples : List[Dict]
            原始样本数据
        
        Returns:
        --------
        tuple: (confirmed_clean, need_manual_review, confirmed_removed)
        """
        print("\n分类复核结果...")
        
        # 创建sample_id到原始数据的映射
        sample_map = {
            self._safe_str(s.get('Sample_id', '')): s 
            for s in original_samples
        }
        
        for result in reverify_results:
            sample_id = result.get('sample_id', '')
            original_data = sample_map.get(sample_id, {})
            
            # 合并原始数据和复核结果
            merged_data = {**original_data, **result}
            
            final_decision = result.get('final_decision', 'UNCERTAIN')
            confidence = result.get('confidence', 0.0)
            
            # 分类逻辑
            if final_decision == 'ACCEPT':
                if confidence >= self.confidence_threshold_high:
                    # 高置信度接受
                    merged_data['reverify_status'] = 'Confirmed Clean (High Confidence)'
                    self.confirmed_clean.append(merged_data)
                elif confidence >= self.confidence_threshold_low:
                    # 中等置信度，需要人工复审
                    merged_data['reverify_status'] = 'Need Manual Review (Medium Confidence Accept)'
                    self.need_manual_review.append(merged_data)
                else:
                    # 低置信度，需要人工复审
                    merged_data['reverify_status'] = 'Need Manual Review (Low Confidence Accept)'
                    self.need_manual_review.append(merged_data)
            
            elif final_decision == 'REJECT':
                if confidence >= self.confidence_threshold_high:
                    # 高置信度拒绝
                    merged_data['reverify_status'] = 'Confirmed Removed (High Confidence)'
                    self.confirmed_removed.append(merged_data)
                elif confidence >= self.confidence_threshold_low:
                    # 中等置信度，需要人工复审
                    merged_data['reverify_status'] = 'Need Manual Review (Medium Confidence Reject)'
                    self.need_manual_review.append(merged_data)
                else:
                    # 低置信度，需要人工复审
                    merged_data['reverify_status'] = 'Need Manual Review (Low Confidence Reject)'
                    self.need_manual_review.append(merged_data)
            
            else:  # UNCERTAIN or ERROR
                merged_data['reverify_status'] = 'Need Manual Review (Uncertain)'
                self.need_manual_review.append(merged_data)
        
        print(f"\n分类完成:")
        print(f"  ✓ 确认为清洁数据: {len(self.confirmed_clean)} 个样本")
        print(f"  ⚠ 需要人工复审: {len(self.need_manual_review)} 个样本")
        print(f"  ✗ 确认应移除: {len(self.confirmed_removed)} 个样本")
        
        return self.confirmed_clean, self.need_manual_review, self.confirmed_removed
    
    def save_reverify_results(self, output_dir: str = "reverification_report"):
        """保存复核结果"""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        print(f"\n保存复核结果到: {output_dir}")
        
        # 1. 确认为清洁数据
        if self.confirmed_clean:
            df_clean = pd.DataFrame(self.confirmed_clean)
            clean_file = os.path.join(output_dir, f'reverify_confirmed_clean_{timestamp}.csv')
            df_clean.to_csv(clean_file, index=False, encoding='utf-8-sig')
            print(f"  ✓ 确认清洁: {clean_file}")
        
        # 2. 需要人工复审
        if self.need_manual_review:
            df_review = pd.DataFrame(self.need_manual_review)
            review_file = os.path.join(output_dir, f'reverify_need_manual_review_{timestamp}.csv')
            df_review.to_csv(review_file, index=False, encoding='utf-8-sig')
            print(f"  ⚠ 需要复审: {review_file}")
            
            # 创建简化版供人工审核
            review_cols = [
                'Sample_id', 'Series_id', 'Title', 'Organism', 
                'Library_strategy', 'reverify_result', 'confidence',
                'reason', 'reverify_status'
            ]
            available_cols = [col for col in review_cols if col in df_review.columns]
            df_review_simple = df_review[available_cols]
            review_simple_file = os.path.join(output_dir, f'reverify_review_simplified_{timestamp}.csv')
            df_review_simple.to_csv(review_simple_file, index=False, encoding='utf-8-sig')
            print(f"  ⚠ 简化版复审表: {review_simple_file}")
        
        # 3. 确认应移除
        if self.confirmed_removed:
            df_removed = pd.DataFrame(self.confirmed_removed)
            removed_file = os.path.join(output_dir, f'reverify_confirmed_removed_{timestamp}.csv')
            df_removed.to_csv(removed_file, index=False, encoding='utf-8-sig')
            print(f"  ✗ 确认移除: {removed_file}")
        
        # 4. 保存详细分析报告
        analysis_report = {
            'total_reverified': len(self.confirmed_clean) + len(self.need_manual_review) + len(self.confirmed_removed),
            'confirmed_clean_count': len(self.confirmed_clean),
            'need_review_count': len(self.need_manual_review),
            'confirmed_removed_count': len(self.confirmed_removed),
            'confidence_threshold_high': self.confidence_threshold_high,
            'confidence_threshold_low': self.confidence_threshold_low,
            'timestamp': timestamp
        }
        
        report_file = os.path.join(output_dir, f'reverify_summary_{timestamp}.json')
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(analysis_report, f, indent=2, ensure_ascii=False)
        
        print(f"  📊 统计报告: {report_file}")
        
        return {
            'clean_file': clean_file if self.confirmed_clean else None,
            'review_file': review_file if self.need_manual_review else None,
            'removed_file': removed_file if self.confirmed_removed else None,
            'report_file': report_file
        }

class SmartReVerifier:
    """智能复核器 - 使用规则+采样+API的混合策略"""
    
    def __init__(self, kimi_client, series_df: pd.DataFrame):
        self.kimi_client = kimi_client
        self.series_df = series_df
        
        # 结果存储
        self.confirmed_clean = []       # 规则直接确认清洁
        self.confirmed_removed = []     # 规则直接确认移除
        self.need_api_check = []        # 需要API检查
        self.need_manual_review = []    # 最终需要人工复审
        
        # 强规则模式
        self.strict_human_keywords = ['homo sapiens', 'human']
        self.strict_exclude_organisms = ['mus musculus', 'mouse', 'rattus', 'rat', 
                                         'drosophila', 'zebrafish', 'danio rerio']
        
        self.strict_sc_keywords = ['10x genomics', 'chromium', 'drop-seq', 
                                   'smart-seq2', 'cel-seq', 'scrna-seq', 'scrnaseq']
        self.strict_bulk_keywords = ['bulk rna-seq', 'bulk transcriptome', 'bulk rna']
        
        self.strict_rna_keywords = ['rna-seq', 'rna seq', 'transcriptome sequencing']
        self.strict_non_rna = ['atac-seq', 'chip-seq', 'dnase-seq', 'wgs', 'wes', 
                               'whole genome sequencing', 'whole exome']
    
    def _safe_str(self, value) -> str:
        if pd.isna(value) or value is None:
            return ''
        return str(value)
    
    def _safe_lower(self, value) -> str:
        return self._safe_str(value).lower()
    
    def calculate_priority_score(self, sample_data: Dict) -> Tuple[float, Dict]:
        """
        计算样本的优先级分数，决定是否需要API检查
        
        Returns:
        --------
        tuple: (priority_score, confidence_indicators)
            priority_score: 0-100，分数越高越需要API检查
            confidence_indicators: 各项指标的详细信息
        """
        score = 0
        indicators = {
            'organism_clarity': 0,     # 物种明确度 (0=不明确, 1=明确人类, -1=明确非人类)
            'sc_clarity': 0,           # 单细胞明确度
            'rna_clarity': 0,          # RNA-seq明确度
            'text_quality': 0,         # 文本信息质量
            'conflict_signals': 0      # 冲突信号数量
        }
        
        # 获取文本信息
        title = self._safe_lower(sample_data.get('Title', ''))
        organism = self._safe_lower(sample_data.get('Organism', ''))
        library_strategy = self._safe_lower(sample_data.get('Library_strategy', ''))
        source = self._safe_lower(sample_data.get('Source_name', ''))
        
        series_id = self._safe_str(sample_data.get('Series_id', ''))
        series_info = self.series_df[self.series_df['Series_id'] == series_id]
        if len(series_info) > 0:
            series_title = self._safe_lower(series_info.iloc[0].get('Title', ''))
            summary = self._safe_lower(series_info.iloc[0].get('Summary', ''))
        else:
            series_title = ''
            summary = ''
        
        full_text = f"{title} {series_title} {summary[:300]} {organism} {source}"
        
        # 1. 物种明确度检查
        if any(kw in organism for kw in self.strict_human_keywords):
            indicators['organism_clarity'] = 1  # 明确人类
        elif any(kw in organism for kw in self.strict_exclude_organisms):
            indicators['organism_clarity'] = -1  # 明确非人类
            score += 0  # 不需要API检查，直接排除
        else:
            indicators['organism_clarity'] = 0  # 不明确
            score += 40  # 需要API检查
        
        # 2. 单细胞明确度
        sc_count = sum(1 for kw in self.strict_sc_keywords if kw in full_text)
        bulk_count = sum(1 for kw in self.strict_bulk_keywords if kw in full_text)
        
        if sc_count >= 2:
            indicators['sc_clarity'] = 1  # 明确单细胞
        elif bulk_count >= 1:
            indicators['sc_clarity'] = -1  # 明确bulk
            score += 0  # 不需要API
        elif sc_count == 1:
            indicators['sc_clarity'] = 0.5  # 可能单细胞
            score += 30
        else:
            indicators['sc_clarity'] = 0  # 不明确
            score += 50  # 高优先级API检查
        
        # 3. RNA-seq明确度
        rna_count = sum(1 for kw in self.strict_rna_keywords if kw in library_strategy or kw in full_text)
        non_rna_count = sum(1 for kw in self.strict_non_rna if kw in full_text)
        
        if rna_count >= 1 and non_rna_count == 0:
            indicators['rna_clarity'] = 1  # 明确RNA-seq
        elif non_rna_count >= 1:
            indicators['rna_clarity'] = -1  # 明确非RNA-seq
            score += 0  # 不需要API
        else:
            indicators['rna_clarity'] = 0  # 不明确
            score += 35
        
        # 4. 文本质量
        text_length = len(summary) + len(title)
        if text_length > 200:
            indicators['text_quality'] = 1
        elif text_length > 50:
            indicators['text_quality'] = 0.5
            score += 10
        else:
            indicators['text_quality'] = 0
            score += 20  # 文本少，难判断，需要API
        
        # 5. 冲突信号检查
        conflicts = 0
        if 'single cell' in full_text and 'bulk' in full_text:
            conflicts += 1
        if 'spatial' in full_text and 'single cell' in full_text:
            conflicts += 1
        if 'multiome' in full_text or 'cite-seq' in full_text:
            conflicts += 1
        
        indicators['conflict_signals'] = conflicts
        score += conflicts * 25  # 有冲突信号，提高优先级
        
        return min(score, 100), indicators
    
    def rule_based_classification(self, sample_data: Dict) -> str:
        """
        基于强规则的分类
        
        Returns:
        --------
        str: 'ACCEPT', 'REJECT', 'UNCERTAIN'
        """
        organism = self._safe_lower(sample_data.get('Organism', ''))
        library_strategy = self._safe_lower(sample_data.get('Library_strategy', ''))
        
        series_id = self._safe_str(sample_data.get('Series_id', ''))
        series_info = self.series_df[self.series_df['Series_id'] == series_id]
        if len(series_info) > 0:
            title = self._safe_lower(series_info.iloc[0].get('Title', ''))
            summary = self._safe_lower(series_info.iloc[0].get('Summary', ''))
        else:
            title = ''
            summary = ''
        
        full_text = f"{title} {summary} {organism} {library_strategy}"
        
        # 强排除规则
        # 1. 非人类物种
        if any(kw in organism for kw in self.strict_exclude_organisms):
            return 'REJECT'
        
        # 2. 明确的bulk数据
        if any(kw in full_text for kw in self.strict_bulk_keywords):
            return 'REJECT'
        
        # 3. 明确的非RNA测序
        if any(kw in full_text for kw in self.strict_non_rna):
            return 'REJECT'
        
        # 强接受规则（需要同时满足）
        is_human = any(kw in organism for kw in self.strict_human_keywords)
        is_sc = any(kw in full_text for kw in self.strict_sc_keywords)
        is_rna = any(kw in full_text for kw in self.strict_rna_keywords)
        
        if is_human and is_sc and is_rna:
            # 检查是否有冲突
            has_conflict = ('spatial' in full_text or 'multiome' in full_text or 
                          'cite-seq' in full_text)
            if not has_conflict:
                return 'ACCEPT'
        
        return 'UNCERTAIN'
    
    def smart_sampling_strategy(self, uncertain_samples: List[Dict], 
                                max_api_calls: int = 1000) -> Tuple[List[Dict], List[Dict]]:
        """
        智能采样策略 - 选择最需要API检查的样本
        
        Parameters:
        -----------
        uncertain_samples : List[Dict]
            不确定的样本列表
        max_api_calls : int
            最大API调用次数
        
        Returns:
        --------
        tuple: (samples_for_api, samples_for_manual)
        """
        print(f"\n开始智能采样 ({len(uncertain_samples)} 个不确定样本)...")
        
        # 计算每个样本的优先级
        samples_with_priority = []
        for sample in uncertain_samples:
            priority_score, indicators = self.calculate_priority_score(sample)
            samples_with_priority.append({
                'sample': sample,
                'priority_score': priority_score,
                'indicators': indicators
            })
        
        # 按优先级排序
        samples_with_priority.sort(key=lambda x: x['priority_score'], reverse=True)
        
        # 选择高优先级样本进行API检查
        samples_for_api = []
        samples_for_manual = []
        
        for item in samples_with_priority:
            if len(samples_for_api) < max_api_calls and item['priority_score'] > 30:
                samples_for_api.append(item['sample'])
            else:
                # 低优先级的直接标记为需要人工复审
                samples_for_manual.append(item['sample'])
        
        print(f"  选择 {len(samples_for_api)} 个高优先级样本进行API检查")
        print(f"  {len(samples_for_manual)} 个低优先级样本直接标记为人工复审")
        
        # 显示优先级分布
        if samples_with_priority:
            scores = [x['priority_score'] for x in samples_with_priority]
            print(f"\n  优先级分数分布:")
            print(f"    最高分: {max(scores):.1f}")
            print(f"    最低分: {min(scores):.1f}")
            print(f"    平均分: {np.mean(scores):.1f}")
            print(f"    中位数: {np.median(scores):.1f}")
        
        return samples_for_api, samples_for_manual
    
    def batch_reverify_with_api(self, samples: List[Dict], batch_size: int = 15) -> List[Dict]:
        
        reverifier = DataReVerifier(
            kimi_client=self.kimi_client,
            series_df=self.series_df
        )
        
        return reverifier.batch_reverify(samples, batch_size)
    
    def classify_all_samples(self, flagged_samples: List[Dict], 
                            removed_samples: List[Dict],
                            max_api_calls: int = 1000) -> Dict:
        """
        智能分类所有样本
        
        Parameters:
        -----------
        flagged_samples : List[Dict]
            待复核样本
        removed_samples : List[Dict]
            已移除样本
        max_api_calls : int
            最大API调用次数预算
        
        Returns:
        --------
        dict: 分类结果统计
        """
        print("\n" + "=" * 80)
        print("智能复核流程 - 三阶段处理")
        print("=" * 80)
        
        all_samples = flagged_samples + removed_samples
        total_samples = len(all_samples)
        
        print(f"\n总样本数: {total_samples}")
        print(f"  - 待复核样本: {len(flagged_samples)}")
        print(f"  - 已移除样本: {len(removed_samples)}")
        print(f"API调用预算: {max_api_calls}")
        
        # ===== 阶段1: 规则分类 =====
        print("\n" + "-" * 80)
        print("阶段1: 基于强规则的快速分类")
        print("-" * 80)
        
        rule_accept = []
        rule_reject = []
        rule_uncertain = []
        
        for idx, sample in enumerate(all_samples, 1):
            if idx % 1000 == 0:
                print(f"\r  进度: {idx}/{total_samples}", end='')
            
            decision = self.rule_based_classification(sample)
            
            if decision == 'ACCEPT':
                sample['classification_method'] = 'Rule-based ACCEPT'
                sample['confidence'] = 0.9
                rule_accept.append(sample)
            elif decision == 'REJECT':
                sample['classification_method'] = 'Rule-based REJECT'
                sample['confidence'] = 0.9
                rule_reject.append(sample)
            else:
                rule_uncertain.append(sample)
        
        print(f"\n\n阶段1完成:")
        print(f"  ✓ 规则确认清洁: {len(rule_accept)}")
        print(f"  ✗ 规则确认移除: {len(rule_reject)}")
        print(f"  ? 需要进一步判断: {len(rule_uncertain)}")
        
        # ===== 阶段2: 智能采样 + API检查 =====
        print("\n" + "-" * 80)
        print("阶段2: 智能采样 + API深度检查")
        print("-" * 80)
        
        api_samples, direct_manual = self.smart_sampling_strategy(
            rule_uncertain, 
            max_api_calls
        )
        
        api_results = []
        if api_samples:
            print(f"\n执行API检查 ({len(api_samples)} 个样本)...")
            api_results = self.batch_reverify_with_api(api_samples, batch_size=15)
        
        # 分类API结果
        api_accept = []
        api_reject = []
        api_uncertain = []
        
        for result in api_results:
            confidence = result.get('confidence', 0.0)
            decision = result.get('final_decision', 'UNCERTAIN')
            
            if decision == 'ACCEPT' and confidence >= 0.85:
                result['classification_method'] = 'API ACCEPT (High Confidence)'
                api_accept.append(result)
            elif decision == 'REJECT' and confidence >= 0.85:
                result['classification_method'] = 'API REJECT (High Confidence)'
                api_reject.append(result)
            else:
                result['classification_method'] = 'API UNCERTAIN'
                api_uncertain.append(result)
        
        print(f"\n阶段2完成:")
        print(f"  ✓ API确认清洁: {len(api_accept)}")
        print(f"  ✗ API确认移除: {len(api_reject)}")
        print(f"  ? API不确定: {len(api_uncertain)}")
        
        # ===== 阶段3: 汇总结果 =====
        print("\n" + "-" * 80)
        print("阶段3: 结果汇总")
        print("-" * 80)
        
        self.confirmed_clean = rule_accept + api_accept
        self.confirmed_removed = rule_reject + api_reject
        self.need_manual_review = api_uncertain + direct_manual
        
        # 标记人工复审的原因
        for sample in direct_manual:
            sample['manual_review_reason'] = 'Low priority - skipped API check'
            sample['classification_method'] = 'Manual review needed (low priority)'
        
        print(f"\n最终分类结果:")
        print(f"  ✓ 确认清洁数据: {len(self.confirmed_clean)} ({len(self.confirmed_clean)/total_samples*100:.1f}%)")
        print(f"  ✗ 确认移除数据: {len(self.confirmed_removed)} ({len(self.confirmed_removed)/total_samples*100:.1f}%)")
        print(f"  ⚠ 需人工复审: {len(self.need_manual_review)} ({len(self.need_manual_review)/total_samples*100:.1f}%)")
        
        print(f"\nAPI调用统计:")
        print(f"  预算: {max_api_calls}")
        print(f"  实际使用: {len(api_samples)}")
        print(f"  节省: {max_api_calls - len(api_samples)} ({(1 - len(api_samples)/max_api_calls)*100:.1f}%)")
        
        return {
            'total': total_samples,
            'confirmed_clean': len(self.confirmed_clean),
            'confirmed_removed': len(self.confirmed_removed),
            'need_manual_review': len(self.need_manual_review),
            'api_calls_used': len(api_samples),
            'api_calls_budget': max_api_calls
        }
    
    def save_results(self, output_dir: str = "smart_reverification_report"):
        """保存结果"""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        print(f"\n保存结果到: {output_dir}")
        
        files = {}
        
        # 1. 确认清洁数据
        if self.confirmed_clean:
            df = pd.DataFrame(self.confirmed_clean)
            file_path = os.path.join(output_dir, f'confirmed_clean_{timestamp}.csv')
            df.to_csv(file_path, index=False, encoding='utf-8-sig')
            files['clean'] = file_path
            print(f"  ✓ 确认清洁: {file_path}")
        
        # 2. 确认移除数据
        if self.confirmed_removed:
            df = pd.DataFrame(self.confirmed_removed)
            file_path = os.path.join(output_dir, f'confirmed_removed_{timestamp}.csv')
            df.to_csv(file_path, index=False, encoding='utf-8-sig')
            files['removed'] = file_path
            print(f"  ✗ 确认移除: {file_path}")
        
        # 3. 需人工复审 - 完整版
        if self.need_manual_review:
            df = pd.DataFrame(self.need_manual_review)
            file_path = os.path.join(output_dir, f'need_manual_review_{timestamp}.csv')
            df.to_csv(file_path, index=False, encoding='utf-8-sig')
            files['review_full'] = file_path
            print(f"  ⚠ 需人工复审(完整): {file_path}")
            
            # 简化版
            review_cols = ['Sample_id', 'Series_id', 'Title', 'Organism', 
                          'Library_strategy', 'classification_method', 
                          'manual_review_reason']
            available_cols = [col for col in review_cols if col in df.columns]
            df_simple = df[available_cols]
            file_path_simple = os.path.join(output_dir, f'need_manual_review_simple_{timestamp}.csv')
            df_simple.to_csv(file_path_simple, index=False, encoding='utf-8-sig')
            files['review_simple'] = file_path_simple
            print(f"  ⚠ 需人工复审(简化): {file_path_simple}")
        
        # 4. 统计报告
        report = {
            'timestamp': timestamp,
            'total_samples': len(self.confirmed_clean) + len(self.confirmed_removed) + len(self.need_manual_review),
            'confirmed_clean': len(self.confirmed_clean),
            'confirmed_removed': len(self.confirmed_removed),
            'need_manual_review': len(self.need_manual_review),
            'clean_rate': len(self.confirmed_clean) / (len(self.confirmed_clean) + len(self.confirmed_removed) + len(self.need_manual_review)) * 100
        }
        
        report_file = os.path.join(output_dir, f'summary_{timestamp}.json')
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        files['report'] = report_file
        print(f"  📊 统计报告: {report_file}")
        
        return files


def main_smart_reverification():
    """主函数 - 智能复核"""
    print("=" * 80)
    print("GEO Smart Re-verification Module v2.0")
    print("智能复核模块 - 成本优化版")
    print("=" * 80)
    
    # 配置
    series_file = "/mnt/public7/pancancercol/meta_analysis/geo_metadata/comprehensive_metadata/complete_human_series_20251125_034124.csv"
    
    # 这里需要替换为实际的文件路径
    flagged_file = "/mnt/public7/pancancercol/meta_analysis/geo_metadata/clean_metadata/cleaning_report/flagged_for_review_20251217_012206.csv"
    removed_file = "/mnt/public7/pancancercol/meta_analysis/geo_metadata/clean_metadata/cleaning_report/removed_data_20251217_012206.csv"
    
    # 读取数据
    print(f"\n加载数据...")
    series_df = pd.read_csv(series_file)
    
    flagged_samples = []
    removed_samples = []
    
    if os.path.exists(flagged_file):
        df = pd.read_csv(flagged_file)
        flagged_samples = df.to_dict('records')
        print(f"  待复核样本: {len(flagged_samples)}")
    
    if os.path.exists(removed_file):
        df = pd.read_csv(removed_file)
        removed_samples = df.to_dict('records')
        print(f"  已移除样本: {len(removed_samples)}")
    
    # 初始化
    from geo_cleaner import KimiAPIClient
    api_keys = [
        "sk-UL5YodR7ZL4S9dytpfMWJgmPTXJjkeNSd7Ktq9bbEhElzDfX",
        "sk-WL1hTKhW3sYKuJjBSc1k8wxL0r8ZLcuM7YiYxFjGgVHqAXhU",
        "sk-RAkA28HIT5tiMEfKtXAgbZ9nZKweq5Bnw0WbSwwBdNX7nbi1"
    ]
    kimi_client = KimiAPIClient(api_keys)
    
    # 创建智能复核器
    reverifier = SmartReVerifier(kimi_client, series_df)
    
    # 执行智能分类
    start_time = time.time()
    
    stats = reverifier.classify_all_samples(
        flagged_samples, 
        removed_samples,
        max_api_calls=1000  # 可根据预算调整
    )
    
    # 保存结果
    files = reverifier.save_results()
    
    elapsed = time.time() - start_time
    
    # 最终报告
    print("\n" + "=" * 80)
    print("智能复核完成!")
    print("=" * 80)
    print(f"\n总耗时: {elapsed:.2f} 秒 ({elapsed/60:.2f} 分钟)")
    print(f"\n成本效益:")
    print(f"  如果全部使用API: {stats['total']} 次调用")
    print(f"  实际API调用: {stats['api_calls_used']} 次")
    print(f"  节省: {stats['total'] - stats['api_calls_used']} 次 ({(1-stats['api_calls_used']/stats['total'])*100:.1f}%)")


if __name__ == "__main__":
    main_smart_reverification()