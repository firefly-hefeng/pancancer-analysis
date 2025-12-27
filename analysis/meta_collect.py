#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEO scRNA-seq Complete Metadata Collector - Human Only
专门收集人类单细胞测序数据的metadata
Version 2.3 - 修复ID转换问题
"""

import os
import re
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from Bio import Entrez
import json
from datetime import datetime
import xml.etree.ElementTree as ET
from collections import defaultdict
import gzip
from io import BytesIO

# 设置你的邮箱（NCBI要求）
Entrez.email = "your_email@example.com"

class ComprehensiveGEOCollector:
    """全面的GEO数据收集器 - 仅收集人类数据"""
    
    def __init__(self, output_dir="geo_comprehensive_data"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # 定义scRNA-seq相关的搜索策略
        self.scrna_keywords = [
            "single cell RNA seq",
            "single-cell RNA seq", 
            "scRNA-seq",
            "scRNAseq",
            "single cell transcriptome",
            "single cell gene expression"
        ]
        
        # 人类物种识别
        self.human_identifiers = [
            "Homo sapiens",
            "homo sapiens", 
            "human",
            "Human",
            "9606"  # NCBI Taxonomy ID
        ]
    
    def convert_gds_id_to_gse(self, gds_ids):
        """
        将GDS数据库的内部ID转换为GSE accession号
        
        Parameters:
        -----------
        gds_ids : list
            GDS数据库ID列表
            
        Returns:
        --------
        list : GSE accession号列表
        """
        gse_list = []
        batch_size = 200
        
        print(f"\n正在将 {len(gds_ids)} 个内部ID转换为GSE号...")
        
        for i in range(0, len(gds_ids), batch_size):
            batch = gds_ids[i:i+batch_size]
            
            try:
                # 使用esummary批量获取摘要信息
                handle = Entrez.esummary(db="gds", id=','.join(batch))
                summaries = Entrez.read(handle)
                handle.close()
                
                for summary in summaries:
                    # 从summary中提取GSE号
                    accession = summary.get('Accession', '')
                    if accession.startswith('GSE'):
                        gse_list.append(accession)
                    elif accession.startswith('GDS'):
                        # GDS需要进一步获取对应的GSE
                        gse = self._get_gse_from_gds(accession)
                        if gse:
                            gse_list.append(gse)
                
                print(f"  已转换 {len(gse_list)} 个GSE号", end='\r')
                time.sleep(0.5)
                
            except Exception as e:
                print(f"\n  转换批次 {i}-{i+batch_size} 失败: {e}")
                continue
        
        print(f"\n转换完成: {len(set(gse_list))} 个唯一GSE号")
        return list(set(gse_list))
    
    def _get_gse_from_gds(self, gds_id):
        """从GDS获取对应的GSE"""
        try:
            handle = Entrez.esearch(db="gds", term=f"{gds_id}[Accession]")
            result = Entrez.read(handle)
            handle.close()
            
            if result["IdList"]:
                handle = Entrez.esummary(db="gds", id=result["IdList"][0])
                summary = Entrez.read(handle)
                handle.close()
                
                # 从relations中找GSE
                if summary and len(summary) > 0:
                    relations = summary[0].get('Relations', [])
                    for relation in relations:
                        if relation.get('RelationType') == 'SuperSeries of':
                            target = relation.get('TargetObject', '')
                            if target.startswith('GSE'):
                                return target
            return None
        except:
            return None
        
    def search_all_scrna_datasets(self, start_date=None, end_date=None, batch_size=500, max_results=None):
        """
        全面搜索GEO中所有人类scRNA-seq数据集
        
        Parameters:
        -----------
        start_date : str
            起始日期 (格式: "YYYY/MM/DD")
        end_date : str
            结束日期 (格式: "YYYY/MM/DD")
        batch_size : int
            每批次获取的数量
        max_results : int
            最多获取多少个结果（None表示全部）
            
        Returns:
        --------
        list : 所有GSE ID列表
        """
        all_gds_ids = set()
        
        # 构建多种搜索策略，添加人类物种过滤
        search_strategies = [
            f'"single cell RNA seq"[All Fields] AND "Homo sapiens"[Organism]',
            f'"single-cell RNA-seq"[All Fields] AND "Homo sapiens"[Organism]',
            f'"scRNA-seq"[All Fields] AND "Homo sapiens"[Organism]',
            f'"10x genomics"[All Fields] AND "single cell"[All Fields] AND "Homo sapiens"[Organism]',
            f'"drop-seq"[All Fields] AND "Homo sapiens"[Organism]',
            f'"smart-seq"[All Fields] AND "Homo sapiens"[Organism]',
            f'"single cell transcriptome"[All Fields] AND "Homo sapiens"[Organism]',
        ]
        
        for strategy in search_strategies:
            print(f"\n使用搜索策略: {strategy}")
            
            # 添加日期限制
            if start_date and end_date:
                date_filter = f' AND "{start_date}"[Publication Date] : "{end_date}"[Publication Date]'
                full_query = strategy + date_filter
            else:
                full_query = strategy
            
            try:
                # 先获取总数
                handle = Entrez.esearch(
                    db="gds",
                    term=full_query,
                    retmax=0
                )
                record = Entrez.read(handle)
                handle.close()
                total_count = int(record["Count"])
                print(f"  找到 {total_count} 个结果")
                
                # 确定实际要获取的数量
                if max_results:
                    total_count = min(total_count, max_results)
                
                # 分批获取所有ID
                for start in range(0, total_count, batch_size):
                    handle = Entrez.esearch(
                        db="gds",
                        term=full_query,
                        retstart=start,
                        retmax=batch_size,
                        sort="relevance"
                    )
                    record = Entrez.read(handle)
                    handle.close()
                    
                    all_gds_ids.update(record["IdList"])
                    print(f"  已收集 {len(all_gds_ids)} 个唯一GDS ID", end='\r')
                    time.sleep(0.5)  # 避免请求过快
                    
            except Exception as e:
                print(f"  搜索出错: {e}")
                continue
        
        print(f"\n\n总计找到 {len(all_gds_ids)} 个唯一的GDS ID")
        
        # 转换GDS ID为GSE号
        gse_list = self.convert_gds_id_to_gse(list(all_gds_ids))
        
        print(f"\n最终得到 {len(gse_list)} 个唯一的人类scRNA-seq数据集（GSE）")
        return sorted(gse_list)
    
    def _is_human(self, text):
        """判断文本中是否包含人类标识"""
        if not text:
            return False
        text_str = str(text).lower()
        return any(identifier.lower() in text_str for identifier in self.human_identifiers)
    
    def fetch_series_metadata(self, gse_id):
        """
        获取Series级别的完整metadata
        
        Parameters:
        -----------
        gse_id : str
            GSE accession号
            
        Returns:
        --------
        dict : Series完整metadata，如果不是人类数据返回None
        """
        # 确保是GSE格式
        if not gse_id.startswith('GSE'):
            gse_id = f'GSE{gse_id}'
        
        metadata = {
            'Series_id': gse_id,
            'Title': '',
            'Summary': '',
            'Overall_Design': '',
            'Platform': [],
            'Sample_Count': 0,
            'Publication_Date': '',
            'Submission_Date': '',
            'Last_Update_Date': '',
            'Contact_Name': '',
            'Contact_Email': '',
            'Contact_Institute': '',
            'PubMed_ID': '',
            'Organism': [],
            'Data_Processing': '',
            'GEO_Link': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}',
            'FTP_Link': '',
            'SRA_Link': '',
            'Supplementary_Files': [],
            'Raw_Metadata': {}
        }
        
        is_human = False  # 假定为非人类，需要验证
        
        try:
            # 方法1: 使用Entrez获取摘要信息
            try:
                handle = Entrez.esearch(db="gds", term=gse_id)
                search_result = Entrez.read(handle)
                handle.close()
                
                if search_result["IdList"]:
                    geo_id = search_result["IdList"][0]
                    handle = Entrez.esummary(db="gds", id=geo_id)
                    summary = Entrez.read(handle)
                    handle.close()
                    
                    if summary:
                        record = summary[0]
                        
                        # 检查物种 - 多种方式
                        taxon = str(record.get('taxon', ''))
                        if taxon == '9606':  # 人类的taxon ID
                            is_human = True
                        
                        # 检查Organism字段
                        organism_list = record.get('Organism', [])
                        if organism_list:
                            metadata['Organism'] = organism_list
                            if self._is_human(' '.join(organism_list)):
                                is_human = True
                        
                        # 保存其他信息
                        metadata['Title'] = record.get('title', '')
                        metadata['Summary'] = record.get('summary', '')
                        metadata['Publication_Date'] = record.get('PDAT', '')
                        metadata['Sample_Count'] = int(record.get('n_samples', 0))
                        metadata['Platform'] = record.get('GPL', '').split(';') if record.get('GPL') else []
                        metadata['PubMed_ID'] = record.get('PubMedIds', [''])[0] if record.get('PubMedIds') else ''
                        metadata['Raw_Metadata']['Entrez_Summary'] = record
            except Exception as e:
                print(f"  Entrez获取失败: {e}")
            
            # 方法2: 解析SOFT格式文件获取详细信息（更可靠）
            soft_data = self._parse_soft_file(gse_id)
            if soft_data:
                # 从SOFT文件检查物种
                if soft_data.get('Organism'):
                    organisms = soft_data['Organism']
                    metadata['Organism'] = organisms
                    if self._is_human(' '.join(organisms)):
                        is_human = True
                
                # 从样本中检查物种
                if soft_data.get('Sample_Organisms'):
                    if self._is_human(' '.join(soft_data['Sample_Organisms'])):
                        is_human = True
                
                metadata.update(soft_data)
            
            # 如果还是不确定，从标题和摘要判断
            if not is_human:
                combined_text = f"{metadata['Title']} {metadata['Summary']}"
                if self._is_human(combined_text):
                    is_human = True
                    print(f"  从标题/摘要判断为人类数据")
            
            # 最终判断
            if not is_human:
                print(f"  {gse_id} 验证为非人类数据，跳过")
                print(f"    Organism: {metadata['Organism']}")
                return None
            
            # 方法3: 获取FTP链接
            metadata['FTP_Link'] = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/'
            
            # 方法4: 检查是否有SRA数据
            sra_link = self._get_sra_link(gse_id)
            if sra_link:
                metadata['SRA_Link'] = sra_link
            
            return metadata
            
        except Exception as e:
            print(f"获取 {gse_id} metadata失败: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def _parse_soft_file(self, gse_id):
        """解析SOFT文件获取详细metadata"""
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&targ=self&form=text&view=full"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            content = response.text
            
            soft_data = {
                'Overall_Design': '',
                'Contact_Name': '',
                'Contact_Email': '',
                'Contact_Institute': '',
                'Data_Processing': '',
                'Submission_Date': '',
                'Last_Update_Date': '',
                'Supplementary_Files': [],
                'Sample_Details': [],
                'Organism': [],
                'Sample_Organisms': []  # 收集样本中的物种信息
            }
            
            current_section = None
            
            for line in content.split('\n'):
                line = line.strip()
                
                if not line or line.startswith('#'):
                    continue
                
                if line.startswith('^SERIES'):
                    current_section = 'SERIES'
                elif line.startswith('^SAMPLE'):
                    current_section = 'SAMPLE'
                    
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    if current_section == 'SERIES':
                        if key == '!Series_overall_design':
                            soft_data['Overall_Design'] += value + ' '
                        elif key == '!Series_contact_name':
                            soft_data['Contact_Name'] = value
                        elif key == '!Series_contact_email':
                            soft_data['Contact_Email'] = value
                        elif key == '!Series_contact_institute':
                            soft_data['Contact_Institute'] = value
                        elif key == '!Series_submission_date':
                            soft_data['Submission_Date'] = value
                        elif key == '!Series_last_update_date':
                            soft_data['Last_Update_Date'] = value
                        elif key == '!Series_supplementary_file':
                            soft_data['Supplementary_Files'].append(value)
                        elif key == '!Series_summary':
                            if 'Summary' not in soft_data or not soft_data.get('Summary'):
                                soft_data['Summary'] = value
                        elif key == '!Series_sample_organism':
                            if value not in soft_data['Organism']:
                                soft_data['Organism'].append(value)
                    
                    elif current_section == 'SAMPLE':
                        if key == '!Sample_organism_ch1':
                            if value not in soft_data['Sample_Organisms']:
                                soft_data['Sample_Organisms'].append(value)
            
            soft_data['Overall_Design'] = soft_data['Overall_Design'].strip()
            return soft_data
            
        except Exception as e:
            print(f"  解析SOFT文件失败: {e}")
            return {}
    
    def _get_sra_link(self, gse_id):
        """获取SRA链接"""
        try:
            handle = Entrez.esearch(db="sra", term=gse_id, retmax=1)
            result = Entrez.read(handle)
            handle.close()
            
            if result["IdList"]:
                return f"https://www.ncbi.nlm.nih.gov/sra?term={gse_id}"
            return ""
            
        except:
            return ""
    
    def fetch_sample_metadata(self, gse_id):
        """
        获取所有样本级别的详细metadata
        
        Parameters:
        -----------
        gse_id : str
            GSE accession号
            
        Returns:
        --------
        list : 样本metadata列表（仅人类样本）
        """
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&targ=gsm&form=text&view=full"
        
        try:
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            content = response.text
            
            samples = []
            current_sample = None
            
            for line in content.split('\n'):
                line = line.strip()
                
                if line.startswith('^SAMPLE'):
                    if current_sample:
                        # 验证是否为人类样本
                        if self._is_human_sample(current_sample):
                            samples.append(current_sample)
                    
                    sample_id = line.split('=')[1].strip()
                    current_sample = {
                        'Series_id': gse_id,
                        'Sample_id': sample_id,
                        'Title': '',
                        'Source_name': '',
                        'Organism': '',
                        'Characteristics': {},
                        'Treatment_protocol': '',
                        'Growth_protocol': '',
                        'Extract_protocol': '',
                        'Library_strategy': '',
                        'Library_source': '',
                        'Library_selection': '',
                        'Instrument_model': '',
                        'Data_processing': '',
                        'Platform_id': '',
                        'Contact_name': '',
                        'Supplementary_files': [],
                        'Raw_data_files': [],
                        'Processed_data_files': []
                    }
                
                elif current_sample and '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip().replace('!Sample_', '')
                    value = value.strip()
                    
                    if key == 'title':
                        current_sample['Title'] = value
                    elif key == 'source_name_ch1':
                        current_sample['Source_name'] = value
                    elif key == 'organism_ch1':
                        current_sample['Organism'] = value
                    elif key == 'characteristics_ch1':
                        # 解析characteristics
                        if ':' in value:
                            char_key, char_value = value.split(':', 1)
                            char_key = char_key.strip().lower().replace(' ', '_')
                            current_sample['Characteristics'][char_key] = char_value.strip()
                    elif key == 'treatment_protocol_ch1':
                        current_sample['Treatment_protocol'] += value + ' '
                    elif key == 'growth_protocol_ch1':
                        current_sample['Growth_protocol'] += value + ' '
                    elif key == 'extract_protocol_ch1':
                        current_sample['Extract_protocol'] += value + ' '
                    elif key == 'library_strategy':
                        current_sample['Library_strategy'] = value
                    elif key == 'library_source':
                        current_sample['Library_source'] = value
                    elif key == 'library_selection':
                        current_sample['Library_selection'] = value
                    elif key == 'instrument_model':
                        current_sample['Instrument_model'] = value
                    elif key == 'data_processing':
                        current_sample['Data_processing'] += value + ' '
                    elif key == 'platform_id':
                        current_sample['Platform_id'] = value
                    elif key == 'contact_name':
                        current_sample['Contact_name'] = value
                    elif key == 'supplementary_file':
                        current_sample['Supplementary_files'].append(value)
                        
                        # 区分原始数据和处理后数据
                        if any(ext in value.lower() for ext in ['.fastq', '.fq', '.bam', '.sam']):
                            current_sample['Raw_data_files'].append(value)
                        elif any(ext in value.lower() for ext in ['.txt', '.csv', '.tsv', '.h5', '.rds', '.mtx']):
                            current_sample['Processed_data_files'].append(value)
            
            # 添加最后一个样本
            if current_sample and self._is_human_sample(current_sample):
                samples.append(current_sample)
            
            return samples
            
        except Exception as e:
            print(f"  获取样本信息失败: {e}")
            return []
    
    def _is_human_sample(self, sample):
        """验证样本是否为人类"""
        organism = sample.get('Organism', '')
        return self._is_human(organism)
    
    def collect_complete_metadata(self, gse_list, delay=2):
        """
        收集完整的原始metadata（Series + Samples），仅人类数据
        
        Parameters:
        -----------
        gse_list : list
            GSE列表
        delay : int
            请求间隔（秒）
            
        Returns:
        --------
        tuple : (series_df, samples_df)
        """
        series_data = []
        all_samples_data = []
        skipped_count = 0
        
        for i, gse_id in enumerate(gse_list, 1):
            print(f"\n[{i}/{len(gse_list)}] 处理: {gse_id}")
            
            # 获取Series metadata
            series_meta = self.fetch_series_metadata(gse_id)
            
            if series_meta is None:
                print(f"  跳过非人类数据")
                skipped_count += 1
                time.sleep(delay)
                continue
            
            series_data.append(series_meta)
            print(f"  Series信息: ✓ (人类数据)")
            
            # 获取Samples metadata
            samples = self.fetch_sample_metadata(gse_id)
            all_samples_data.extend(samples)
            print(f"  样本信息: {len(samples)} 个人类样本 ✓")
            
            # 保存中间结果（每10个数据集）
            if i % 10 == 0:
                self._save_intermediate_results(series_data, all_samples_data, i)
            
            # 延迟
            if i < len(gse_list):
                time.sleep(delay)
        
        print(f"\n总计: {len(series_data)} 个人类数据集, 跳过 {skipped_count} 个非人类数据集")
        
        # 转换为DataFrame
        series_df = pd.DataFrame(series_data)
        samples_df = pd.DataFrame(all_samples_data)
        
        return series_df, samples_df
    
    def _save_intermediate_results(self, series_data, samples_data, count):
        """保存中间结果"""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        if series_data:
            series_df = pd.DataFrame(series_data)
            series_file = os.path.join(self.output_dir, f'series_intermediate_{count}_{timestamp}.csv')
            series_df.to_csv(series_file, index=False, encoding='utf-8-sig')
        
        if samples_data:
            samples_df = pd.DataFrame(samples_data)
            samples_file = os.path.join(self.output_dir, f'samples_intermediate_{count}_{timestamp}.csv')
            samples_df.to_csv(samples_file, index=False, encoding='utf-8-sig')
    
    def save_complete_metadata(self, series_df, samples_df, prefix="complete_human"):
        """保存完整的原始metadata"""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # 如果没有数据，创建空DataFrame但保留列名
        if len(series_df) == 0:
            print("\n警告: 没有收集到人类数据，将创建空文件")
            series_df = pd.DataFrame(columns=['Series_id', 'Title', 'Summary'])
        
        if len(samples_df) == 0:
            samples_df = pd.DataFrame(columns=['Series_id', 'Sample_id', 'Title'])
        
        # 保存Series数据
        series_excel = os.path.join(self.output_dir, f'{prefix}_series_{timestamp}.xlsx')
        series_csv = os.path.join(self.output_dir, f'{prefix}_series_{timestamp}.csv')
        
        series_df.to_excel(series_excel, index=False, engine='openpyxl')
        series_df.to_csv(series_csv, index=False, encoding='utf-8-sig')
        
        print(f"\nSeries数据已保存:")
        print(f"  Excel: {series_excel}")
        print(f"  CSV: {series_csv}")
        
        # 保存Samples数据
        samples_excel = os.path.join(self.output_dir, f'{prefix}_samples_{timestamp}.xlsx')
        samples_csv = os.path.join(self.output_dir, f'{prefix}_samples_{timestamp}.csv')
        
        samples_df.to_excel(samples_excel, index=False, engine='openpyxl')
        samples_df.to_csv(samples_csv, index=False, encoding='utf-8-sig')
        
        print(f"\nSamples数据已保存:")
        print(f"  Excel: {samples_excel}")
        print(f"  CSV: {samples_csv}")
        
        return series_csv, samples_csv  # 返回CSV文件路径


class MetadataTransformer:
    """将原始metadata转换为目标格式 - 增强错误处理"""
    
    def __init__(self, series_file, samples_file):
        """
        Parameters:
        -----------
        series_file : str
            Series metadata文件路径
        samples_file : str
            Samples metadata文件路径
        """
        # 检查文件是否为空
        if os.path.getsize(series_file) == 0:
            print("警告: Series文件为空，创建空DataFrame")
            self.series_df = pd.DataFrame()
        else:
            self.series_df = pd.read_csv(series_file) if series_file.endswith('.csv') else pd.read_excel(series_file)
        
        if os.path.getsize(samples_file) == 0:
            print("警告: Samples文件为空，创建空DataFrame")
            self.samples_df = pd.DataFrame()
        else:
            self.samples_df = pd.read_csv(samples_file) if samples_file.endswith('.csv') else pd.read_excel(samples_file)
        
        # 目标字段定义
        self.target_fields = [
            'Cancer_general', 'Cancer', 'Supplementary_information',
            'Molecular_subtype', 'ID', 'Sample_ID', 'Tissue',
            'TumorNormal', 'Filter', 'Sorting_criteria', 'Metastasis',
            'Gender', 'Nationality', 'Age', 'Stage', 'Location',
            'TNM', 'Treatment_phase', 'TCGA_name', 'Availability',
            'Series_id', 'Sample_id', 'Samples_Count', 'Data_Size',
            'Description', 'Platform', 'Publication_Date', 'GEO_Link', 'FTP_Link'
        ]
    
    def _safe_str(self, value):
        """安全地将值转换为字符串"""
        if pd.isna(value) or value is None:
            return ''
        if isinstance(value, (int, float)):
            if pd.isna(value):
                return ''
            return str(value)
        return str(value)
    
    def _safe_lower(self, value):
        """安全地转换为小写"""
        value_str = self._safe_str(value)
        return value_str.lower() if value_str else ''
    
    def transform_to_target_format(self):
        """
        将原始metadata转换为目标格式
        
        Returns:
        --------
        pd.DataFrame : 目标格式的数据
        """
        # 检查是否有数据
        if len(self.series_df) == 0 or len(self.samples_df) == 0:
            print("\n警告: 没有数据可转换，返回空DataFrame")
            return pd.DataFrame(columns=self.target_fields)
        
        print("\n开始转换metadata到目标格式...")
        
        transformed_data = []
        
        # 按Series分组处理
        for series_id in self.samples_df['Series_id'].unique():
            series_info = self.series_df[self.series_df['Series_id'] == series_id]
            
            if len(series_info) == 0:
                print(f"\n警告: {series_id} 没有找到对应的Series信息，跳过")
                continue
            
            series_info = series_info.iloc[0]
            samples = self.samples_df[self.samples_df['Series_id'] == series_id]
            
            print(f"\n处理 {series_id}: {len(samples)} 个样本")
            
            for _, sample in samples.iterrows():
                record = {}
                
                # ===== 直接映射字段 =====
                record['Series_id'] = series_id
                record['ID'] = series_id
                record['Sample_id'] = self._safe_str(sample.get('Sample_id', ''))
                record['Sample_ID'] = self._safe_str(sample.get('Sample_id', ''))
                record['Samples_Count'] = int(series_info.get('Sample_Count', 0)) if pd.notna(series_info.get('Sample_Count')) else 0
                record['Description'] = self._safe_str(series_info.get('Summary', ''))
                
                # Platform处理
                platform_value = series_info.get('Platform', [])
                if isinstance(platform_value, str):
                    try:
                        platform_value = eval(platform_value)
                    except:
                        platform_value = [platform_value] if platform_value else []
                record['Platform'] = ', '.join(platform_value) if isinstance(platform_value, list) else self._safe_str(platform_value)
                
                record['Publication_Date'] = self._safe_str(series_info.get('Publication_Date', ''))
                record['GEO_Link'] = self._safe_str(series_info.get('GEO_Link', ''))
                record['FTP_Link'] = self._safe_str(series_info.get('FTP_Link', ''))
                
                # ===== 从样本Characteristics中提取字段 =====
                characteristics = self._parse_characteristics(sample.get('Characteristics', '{}'))
                
                # 组织类型
                record['Tissue'] = self._extract_tissue(characteristics, sample)
                
                # 肿瘤/正常
                record['TumorNormal'] = self._extract_tumor_normal(characteristics, sample)
                
                # 性别
                record['Gender'] = self._extract_gender(characteristics)
                
                # 年龄
                record['Age'] = self._extract_age(characteristics)
                
                # 分期
                record['Stage'] = self._extract_stage(characteristics)
                
                # TNM
                record['TNM'] = self._extract_tnm(characteristics)
                
                # 转移
                record['Metastasis'] = self._extract_metastasis(characteristics)
                
                # 分子亚型
                record['Molecular_subtype'] = self._extract_molecular_subtype(characteristics)
                
                # 治疗阶段 - 使用安全方法
                record['Treatment_phase'] = self._extract_treatment_phase(characteristics, sample)
                
                # 位置
                record['Location'] = self._extract_location(characteristics)
                
                # ===== 需要手动补充或推断的字段 =====
                record['Cancer_general'] = self._infer_cancer_general(series_info, sample)
                record['Cancer'] = self._infer_cancer_type(series_info, sample)
                record['TCGA_name'] = self._map_to_tcga(record['Cancer'])
                record['Supplementary_information'] = self._safe_str(series_info.get('Overall_Design', ''))
                record['Filter'] = self._extract_filter_criteria(sample)
                record['Sorting_criteria'] = self._extract_sorting_criteria(sample)
                record['Nationality'] = self._extract_nationality(series_info)
                record['Availability'] = 'Public'
                record['Data_Size'] = self._estimate_data_size(sample)
                
                transformed_data.append(record)
        
        # 创建DataFrame
        df = pd.DataFrame(transformed_data, columns=self.target_fields)
        
        print(f"\n转换完成: {len(df)} 条记录")
        return df
    
    def _parse_characteristics(self, characteristics_str):
        """解析characteristics字符串为字典"""
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
    
    def _extract_tissue(self, characteristics, sample):
        """提取组织类型"""
        tissue_keys = ['tissue', 'tissue_type', 'source_name', 'cell_type', 'organ']
        
        for key in tissue_keys:
            if key in characteristics:
                return self._safe_str(characteristics[key])
        
        # 从source_name提取
        source = self._safe_str(sample.get('Source_name', ''))
        if source:
            return source
        
        return ''
    
    def _extract_tumor_normal(self, characteristics, sample):
        """判断肿瘤/正常"""
        tumor_keys = ['sample_type', 'tissue_type', 'disease_state', 'condition']
        
        for key in tumor_keys:
            if key in characteristics:
                value = self._safe_lower(characteristics[key])
                if any(term in value for term in ['tumor', 'tumour', 'cancer', 'malignant', 'carcinoma']):
                    return 'T'
                elif any(term in value for term in ['normal', 'healthy', 'control', 'non-tumor']):
                    return 'N'
        
        # 从标题推断
        title = self._safe_lower(sample.get('Title', ''))
        if any(term in title for term in ['tumor', 'cancer', 'malignant']):
            return 'T'
        elif any(term in title for term in ['normal', 'healthy', 'control']):
            return 'N'
        
        return ''
    
    def _extract_gender(self, characteristics):
        """提取性别"""
        gender_keys = ['gender', 'sex']
        
        for key in gender_keys:
            if key in characteristics:
                gender = self._safe_lower(characteristics[key])
                if gender in ['male', 'm']:
                    return 'M'
                elif gender in ['female', 'f']:
                    return 'F'
                return self._safe_str(characteristics[key])
        
        return ''
    
    def _extract_age(self, characteristics):
        """提取年龄"""
        age_keys = ['age', 'age_at_diagnosis', 'patient_age']
        
        for key in age_keys:
            if key in characteristics:
                age_str = self._safe_str(characteristics[key])
                return self._parse_age_range(age_str)
        
        return ''
    
    def _parse_age_range(self, age_str):
        """解析年龄范围"""
        if not age_str:
            return ''
        
        # 提取数字
        numbers = re.findall(r'\d+', str(age_str))
        if numbers:
            age = int(numbers[0])
            if age < 20:
                return '<20'
            elif age <= 30:
                return '20-30'
            elif age <= 40:
                return '31-40'
            elif age <= 50:
                return '41-50'
            elif age <= 60:
                return '51-60'
            elif age <= 70:
                return '61-70'
            else:
                return '>70'
        
        return age_str
    
    def _extract_stage(self, characteristics):
        """提取分期"""
        stage_keys = ['stage', 'tumor_stage', 'pathologic_stage', 'clinical_stage']
        
        for key in stage_keys:
            if key in characteristics:
                return self._safe_str(characteristics[key])
        
        return ''
    
    def _extract_tnm(self, characteristics):
        """提取TNM分期"""
        tnm_keys = ['tnm', 'tnm_stage', 'tnm_classification']
        
        for key in tnm_keys:
            if key in characteristics:
                return self._safe_str(characteristics[key])
        
        # 尝试组合T、N、M
        t = self._safe_str(characteristics.get('t_stage', ''))
        n = self._safe_str(characteristics.get('n_stage', ''))
        m = self._safe_str(characteristics.get('m_stage', ''))
        
        if t or n or m:
            return f"{t}{n}{m}".strip()
        
        return ''
    
    def _extract_metastasis(self, characteristics):
        """提取转移状态"""
        metastasis_keys = ['metastasis', 'metastatic', 'm_stage']
        
        for key in metastasis_keys:
            if key in characteristics:
                value = self._safe_lower(characteristics[key])
                if value in ['yes', 'true', '1', 'positive', 'm1']:
                    return 'T'
                elif value in ['no', 'false', '0', 'negative', 'm0']:
                    return 'F'
                return value
        
        return ''
    
    def _extract_molecular_subtype(self, characteristics):
        """提取分子亚型"""
        subtype_keys = ['molecular_subtype', 'subtype', 'pam50', 'receptor_status']
        
        for key in subtype_keys:
            if key in characteristics:
                return self._safe_str(characteristics[key])
        
        # 对于乳腺癌，尝试组合ER/PR/HER2状态
        er = self._safe_str(characteristics.get('er_status', characteristics.get('er', '')))
        pr = self._safe_str(characteristics.get('pr_status', characteristics.get('pr', '')))
        her2 = self._safe_str(characteristics.get('her2_status', characteristics.get('her2', '')))
        
        if er or pr or her2:
            subtype_parts = []
            if er:
                subtype_parts.append(f"ER{'+' if 'pos' in er.lower() else '-'}")
            if pr:
                subtype_parts.append(f"PR{'+' if 'pos' in pr.lower() else '-'}")
            if her2:
                subtype_parts.append(f"HER2{'+' if 'pos' in her2.lower() else '-'}")
            
            return '/'.join(subtype_parts)
        
        return ''
    
    def _extract_treatment_phase(self, characteristics, sample):
        """提取治疗阶段 - 使用安全方法"""
        treatment_keys = ['treatment', 'treatment_status', 'therapy', 'timepoint']
        
        for key in treatment_keys:
            if key in characteristics:
                value = self._safe_lower(characteristics[key])
                if 'pre' in value or 'before' in value or 'naive' in value:
                    return 'Pre'
                elif 'post' in value or 'after' in value:
                    return 'Post'
                elif 'during' in value or 'on' in value:
                    return 'During'
                return self._safe_str(characteristics[key])
        
        # 从protocol提取 - 使用安全方法
        treatment_protocol = sample.get('Treatment_protocol', '')
        treatment_protocol_str = self._safe_lower(treatment_protocol)
        
        if treatment_protocol_str:
            if 'untreated' in treatment_protocol_str or 'naive' in treatment_protocol_str:
                return 'Pre'
        
        return ''
    
    def _extract_location(self, characteristics):
        """提取肿瘤位置"""
        location_keys = ['location', 'site', 'anatomic_site', 'tumor_location']
        
        for key in location_keys:
            if key in characteristics:
                return self._safe_str(characteristics[key])
        
        return ''
    
    def _infer_cancer_general(self, series_info, sample):
        """推断癌症大类"""
        title = self._safe_lower(series_info.get('Title', ''))
        summary = self._safe_lower(series_info.get('Summary', ''))
        
        cancer_mapping = {
            'breast': 'BC',
            'lung': 'LC',
            'colorectal': 'CRC',
            'gastric': 'GC',
            'liver': 'HCC',
            'pancreatic': 'PDAC',
            'prostate': 'PRAD',
            'melanoma': 'MEL',
            'glioma': 'GBM',
            'ovarian': 'OV',
        }
        
        full_text = title + ' ' + summary
        
        for keyword, abbr in cancer_mapping.items():
            if keyword in full_text:
                return abbr
        
        return ''
    
    def _infer_cancer_type(self, series_info, sample):
        """推断详细癌症类型（TCGA格式）"""
        title = self._safe_lower(series_info.get('Title', ''))
        summary = self._safe_lower(series_info.get('Summary', ''))
        
        tcga_mapping = {
            'breast': 'BRCA',
            'lung adenocarcinoma': 'LUAD',
            'lung squamous': 'LUSC',
            'colorectal': 'COAD',
            'colon': 'COAD',
            'rectum': 'READ',
            'gastric': 'STAD',
            'liver': 'LIHC',
            'hepatocellular': 'LIHC',
            'pancreatic': 'PAAD',
            'prostate': 'PRAD',
            'melanoma': 'SKCM',
            'glioblastoma': 'GBM',
            'ovarian': 'OV',
        }
        
        full_text = title + ' ' + summary
        
        for keyword, tcga in tcga_mapping.items():
            if keyword in full_text:
                return tcga
        
        return ''
    
    def _map_to_tcga(self, cancer_type):
        """映射到TCGA名称"""
        return cancer_type  # 已经是TCGA格式
    
    def _extract_filter_criteria(self, sample):
        """提取过滤标准"""
        library_selection = self._safe_str(sample.get('Library_selection', ''))
        if library_selection:
            return library_selection
        return ''
    
    def _extract_sorting_criteria(self, sample):
        """提取分选标准"""
        # 从protocol中提取
        extract_protocol = self._safe_lower(sample.get('Extract_protocol', ''))
        
        if 'facs' in extract_protocol:
            return 'FACS'
        elif 'magnetic' in extract_protocol or 'macs' in extract_protocol:
            return 'MACS'
        elif 'cd45' in extract_protocol:
            return 'CD45+'
        
        return ''
    
    def _extract_nationality(self, series_info):
        """提取国籍"""
        contact_institute = self._safe_str(series_info.get('Contact_Institute', ''))
        
        country_keywords = {
            'USA': ['USA', 'United States', 'America'],
            'China': ['China', 'Chinese'],
            'UK': ['United Kingdom', 'UK', 'England'],
            'Germany': ['Germany', 'German'],
            'France': ['France', 'French'],
            'Japan': ['Japan', 'Japanese'],
        }
        
        for country, keywords in country_keywords.items():
            if any(kw in contact_institute for kw in keywords):
                return country
        
        return ''
    
    def _estimate_data_size(self, sample):
        """估算数据大小"""
        # 基于文件列表估算
        raw_files_str = self._safe_str(sample.get('Raw_data_files', '[]'))
        processed_files_str = self._safe_str(sample.get('Processed_data_files', '[]'))
        
        try:
            raw_files = eval(raw_files_str) if raw_files_str and raw_files_str != '[]' else []
            processed_files = eval(processed_files_str) if processed_files_str and processed_files_str != '[]' else []
        except:
            raw_files = []
            processed_files = []
        
        total_files = len(raw_files) + len(processed_files)
        
        if total_files == 0:
            return ''
        elif total_files < 5:
            return '<1GB'
        elif total_files < 10:
            return '1-5GB'
        else:
            return '>5GB'
    
    def save_transformed_data(self, df, output_dir="transformed_metadata"):
        """保存转换后的数据"""
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        excel_file = os.path.join(output_dir, f'transformed_human_metadata_{timestamp}.xlsx')
        csv_file = os.path.join(output_dir, f'transformed_human_metadata_{timestamp}.csv')
        
        df.to_excel(excel_file, index=False, engine='openpyxl')
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        
        print(f"\n转换后的数据已保存:")
        print(f"  Excel: {excel_file}")
        print(f"  CSV: {csv_file}")
        
        return excel_file


def main():
    """主函数 - 完整工作流（仅人类数据）"""
    print("=" * 80)
    print("GEO Human scRNA-seq Complete Metadata Collection Pipeline")
    print("仅收集人类（Homo sapiens）单细胞测序数据")
    print("=" * 80)
    
    # ===== 第一步：收集完整的原始metadata（仅人类） =====
    print("\n" + "=" * 80)
    print("第一步：收集GEO数据库中所有人类scRNA-seq数据的完整metadata")
    print("=" * 80)
    
    collector = ComprehensiveGEOCollector()
    
    # 选项1: 全面搜索（可能需要很长时间）- 限制数量用于测试
    gse_list = collector.search_all_scrna_datasets()
    
    
    print(f"\n准备收集 {len(gse_list)} 个数据集（将验证是否为人类数据）")
    
    # 收集完整metadata
    series_df, samples_df = collector.collect_complete_metadata(gse_list, delay=2)
    
    # 保存完整的原始metadata - 返回CSV路径
    series_file, samples_file = collector.save_complete_metadata(series_df, samples_df)
    
    print(f"\n第一步完成！")
    print(f"  人类Series数据: {len(series_df)} 个数据集")
    print(f"  人类Samples数据: {len(samples_df)} 个样本")
    
    # 如果没有收集到数据，提前退出
    if len(series_df) == 0 or len(samples_df) == 0:
        print("\n警告: 没有收集到人类数据，程序结束")
        return
    
    # ===== 第二步：转换为目标格式 =====
    print("\n" + "=" * 80)
    print("第二步：将原始metadata转换为目标格式")
    print("=" * 80)
    
    transformer = MetadataTransformer(series_file, samples_file)
    
    # 转换数据
    transformed_df = transformer.transform_to_target_format()
    
    # 保存转换后的数据
    output_file = transformer.save_transformed_data(transformed_df)
    
    print(f"\n第二步完成！")
    print(f"  转换后数据: {len(transformed_df)} 条记录")
    
    # ===== 数据质量报告 =====
    if len(transformed_df) > 0:
        print("\n" + "=" * 80)
        print("数据质量报告")
        print("=" * 80)
        
        print(f"\n字段完整性:")
        for col in transformed_df.columns:
            non_empty = transformed_df[col].notna().sum()
            non_empty_non_blank = (transformed_df[col].notna() & (transformed_df[col] != '')).sum()
            percentage = (non_empty_non_blank / len(transformed_df)) * 100
            print(f"  {col}: {non_empty_non_blank}/{len(transformed_df)} ({percentage:.1f}%)")
    
    print("\n" + "=" * 80)
    print("所有步骤完成！")
    print("=" * 80)


if __name__ == "__main__":
    main()