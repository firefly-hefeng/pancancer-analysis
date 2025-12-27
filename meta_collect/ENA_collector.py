"""
ENA_collector.py
European Nucleotide Archive (ENA) 数据收集器
包含来自DDBJ的镜像数据
"""

import os
import json
import time
import requests
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
import logging
from tqdm import tqdm

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('ena_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ENACollector:
    """ENA数据收集器"""
    
    def __init__(self, output_dir="ENA_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.base_url = "https://www.ebi.ac.uk/ena/portal/api/search"
        
    def collect_studies(self, max_results=5000):
        """收集ENA研究数据"""
        logger.info("开始收集ENA研究数据...")
        
        # 构建查询 - 人类单细胞RNA测序
        query = (
            'tax_eq(9606) AND '  # 人类
            '(experiment_title="*single*cell*" OR experiment_title="*scRNA*" OR '
            'study_title="*single*cell*" OR study_title="*scRNA*" OR '
            'library_name="*single*cell*") AND '
            'library_strategy="RNA-Seq"'
        )
        
        try:
            params = {
                'result': 'study',
                'query': query,
                'fields': 'study_accession,secondary_study_accession,study_title,study_alias,study_description,center_name,first_public,last_updated,tax_id,scientific_name',
                'format': 'json',
                'limit': max_results
            }
            
            logger.info(f"查询ENA Portal API...")
            response = requests.get(self.base_url, params=params, timeout=60)
            
            if response.status_code == 200:
                studies = response.json()
                logger.info(f"找到 {len(studies)} 个研究")
                
                # 保存原始JSON
                json_file = self.output_dir / 'ena_studies_raw.json'
                with open(json_file, 'w', encoding='utf-8') as f:
                    json.dump(studies, f, indent=2, ensure_ascii=False)
                logger.info(f"原始数据已保存: {json_file}")
                
                # 保存为CSV
                if studies:
                    df = pd.DataFrame(studies)
                    csv_file = self.output_dir / 'ena_studies.csv'
                    df.to_csv(csv_file, index=False, encoding='utf-8-sig')
                    logger.info(f"CSV已保存: {csv_file}")
                
                return studies
            else:
                logger.error(f"请求失败: {response.status_code}")
                logger.error(f"响应内容: {response.text}")
                return []
                
        except Exception as e:
            logger.error(f"收集研究数据时出错: {str(e)}", exc_info=True)
            return []
    
    def collect_experiments(self, studies: List[Dict], max_per_study=1000):
        """收集实验和样本数据"""
        logger.info(f"开始收集 {len(studies)} 个研究的实验数据...")
        
        all_experiments = []
        failed_studies = []
        
        for study in tqdm(studies, desc="收集实验"):
            study_acc = study.get('study_accession')
            if not study_acc:
                continue
            
            try:
                params = {
                    'result': 'read_experiment',
                    'query': f'study_accession="{study_acc}"',
                    'fields': 'experiment_accession,experiment_title,experiment_alias,library_name,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,sample_accession,sample_title,sample_alias,scientific_name',
                    'format': 'json',
                    'limit': max_per_study
                }
                
                response = requests.get(self.base_url, params=params, timeout=30)
                
                if response.status_code == 200:
                    experiments = response.json()
                    
                    # 添加研究信息到每个实验
                    for exp in experiments:
                        exp['parent_study_accession'] = study_acc
                        exp['parent_study_title'] = study.get('study_title', '')
                    
                    all_experiments.extend(experiments)
                else:
                    failed_studies.append(study_acc)
                    logger.warning(f"研究 {study_acc} 请求失败: {response.status_code}")
                
                time.sleep(0.3)  # 避免请求过快
                
            except Exception as e:
                failed_studies.append(study_acc)
                logger.error(f"获取研究 {study_acc} 的实验时出错: {str(e)}")
        
        logger.info(f"共收集到 {len(all_experiments)} 个实验")
        if failed_studies:
            logger.warning(f"失败的研究数: {len(failed_studies)}")
        
        # 保存实验数据
        if all_experiments:
            json_file = self.output_dir / 'ena_experiments_raw.json'
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(all_experiments, f, indent=2, ensure_ascii=False)
            
            df = pd.DataFrame(all_experiments)
            csv_file = self.output_dir / 'ena_experiments.csv'
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            
            logger.info(f"实验数据已保存: {csv_file}")
        
        # 保存失败列表
        if failed_studies:
            failed_file = self.output_dir / 'failed_studies.txt'
            with open(failed_file, 'w') as f:
                f.write('\n'.join(failed_studies))
        
        return all_experiments
    
    def collect_runs(self, experiments: List[Dict], max_per_experiment=1000):
        """收集测序运行数据"""
        logger.info(f"开始收集 {len(experiments)} 个实验的运行数据...")
        
        # 获取唯一的实验登录号
        experiment_accs = list(set([exp.get('experiment_accession') for exp in experiments if exp.get('experiment_accession')]))
        
        all_runs = []
        
        # 分批查询
        batch_size = 50
        for i in tqdm(range(0, len(experiment_accs), batch_size), desc="收集运行"):
            batch = experiment_accs[i:i+batch_size]
            
            try:
                # 构建查询
                query = ' OR '.join([f'experiment_accession="{acc}"' for acc in batch])
                
                params = {
                    'result': 'read_run',
                    'query': query,
                    'fields': 'run_accession,experiment_accession,sample_accession,study_accession,instrument_platform,instrument_model,library_layout,library_strategy,library_source,read_count,base_count,fastq_ftp,fastq_md5,submitted_ftp',
                    'format': 'json',
                    'limit': max_per_experiment * len(batch)
                }
                
                response = requests.get(self.base_url, params=params, timeout=30)
                
                if response.status_code == 200:
                    runs = response.json()
                    all_runs.extend(runs)
                
                time.sleep(0.5)
                
            except Exception as e:
                logger.error(f"获取批次运行数据时出错: {str(e)}")
        
        logger.info(f"共收集到 {len(all_runs)} 个运行")
        
        # 保存运行数据
        if all_runs:
            json_file = self.output_dir / 'ena_runs_raw.json'
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(all_runs, f, indent=2, ensure_ascii=False)
            
            df = pd.DataFrame(all_runs)
            csv_file = self.output_dir / 'ena_runs.csv'
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            
            logger.info(f"运行数据已保存: {csv_file}")
        
        return all_runs
    
    def generate_summary(self):
        """生成数据摘要"""
        summary_file = self.output_dir / 'collection_summary.txt'
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("ENA数据收集摘要\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"收集时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # 统计研究数
            studies_file = self.output_dir / 'ena_studies.csv'
            if studies_file.exists():
                df_studies = pd.read_csv(studies_file)
                f.write(f"研究总数: {len(df_studies)}\n")
                f.write(f"研究中心数: {df_studies['center_name'].nunique()}\n\n")
            
            # 统计实验数
            exp_file = self.output_dir / 'ena_experiments.csv'
            if exp_file.exists():
                df_exp = pd.read_csv(exp_file)
                f.write(f"实验总数: {len(df_exp)}\n")
                f.write(f"样本总数: {df_exp['sample_accession'].nunique()}\n")
                
                if 'instrument_platform' in df_exp.columns:
                    f.write("\n测序平台分布:\n")
                    platform_counts = df_exp['instrument_platform'].value_counts()
                    for platform, count in platform_counts.items():
                        f.write(f"  {platform}: {count}\n")
                
                if 'library_layout' in df_exp.columns:
                    f.write("\n文库布局分布:\n")
                    layout_counts = df_exp['library_layout'].value_counts()
                    for layout, count in layout_counts.items():
                        f.write(f"  {layout}: {count}\n")
            
            # 统计运行数
            runs_file = self.output_dir / 'ena_runs.csv'
            if runs_file.exists():
                df_runs = pd.read_csv(runs_file)
                f.write(f"\n运行总数: {len(df_runs)}\n")
                
                if 'read_count' in df_runs.columns:
                    total_reads = df_runs['read_count'].sum()
                    f.write(f"总读数: {total_reads:,}\n")
                
                if 'base_count' in df_runs.columns:
                    total_bases = df_runs['base_count'].sum()
                    f.write(f"总碱基数: {total_bases:,}\n")
        
        logger.info(f"摘要已生成: {summary_file}")
    
    def run_full_collection(self):
        """运行完整的收集流程"""
        logger.info("=" * 80)
        logger.info("开始ENA完整数据收集")
        logger.info("=" * 80)
        
        # 1. 收集研究
        studies = self.collect_studies()
        
        if not studies:
            logger.error("未找到研究数据,终止收集")
            return
        
        # 2. 收集实验
        experiments = self.collect_experiments(studies)
        
        # 3. 收集运行(可选,数据量可能很大)
        logger.info("\n是否收集运行数据? 运行数据包含FASTQ文件链接,数据量较大")
        # 这里可以设置是否收集运行数据
        collect_runs = True  # 改为False可跳过
        
        if collect_runs and experiments:
            runs = self.collect_runs(experiments)
        
        # 4. 生成摘要
        self.generate_summary()
        
        logger.info("=" * 80)
        logger.info("ENA数据收集完成!")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info("=" * 80)


def main():
    """主函数"""
    print("=" * 80)
    print("ENA/DDBJ 单细胞RNA测序数据收集器")
    print("=" * 80)
    print("\n此脚本将收集:")
    print("  1. 研究(Studies)元数据")
    print("  2. 实验(Experiments)元数据")
    print("  3. 运行(Runs)元数据和数据文件链接")
    print("\n数据来源:")
    print("  - ENA (European Nucleotide Archive)")
    print("  - DDBJ (镜像数据)")
    print("=" * 80)
    
    collector = ENACollector()
    collector.run_full_collection()


if __name__ == "__main__":
    main()