import pandas as pd
import cellxgene_census
import json
from typing import List, Dict
import time
from tqdm import tqdm
import warnings
import os
warnings.filterwarnings('ignore')

class CellxGeneSampleExtractor:
    """
    基于已有的Collection和Dataset信息提取Sample层面数据
    专门针对人源单细胞RNA-seq数据
    """
    
    def __init__(self, collections_file: str, datasets_file: str):
        """
        初始化
        
        Parameters:
        -----------
        collections_file : str
            Collections表格文件路径 (CSV)
        datasets_file : str
            Datasets表格文件路径 (CSV)
        """
        print("="*80)
        print("CellxGene Sample信息提取工具 - 人源单细胞RNA-seq数据")
        print("="*80 + "\n")
        
        # 加载数据
        self.collections_df = pd.read_csv(collections_file)
        self.datasets_df = pd.read_csv(datasets_file)
        
        print(f"✓ 加载Collections数据: {len(self.collections_df)} 条记录")
        print(f"  文件路径: {collections_file}\n")
        
        print(f"✓ 加载Datasets数据: {len(self.datasets_df)} 条记录")
        print(f"  文件路径: {datasets_file}\n")
        
        # 验证关键列是否存在
        required_cols_collections = ['collection_id']
        required_cols_datasets = ['collection_id', 'dataset_id', 'organisms', 'assays']
        
        missing_cols = []
        for col in required_cols_collections:
            if col not in self.collections_df.columns:
                missing_cols.append(f"collections表缺少列: {col}")
        
        for col in required_cols_datasets:
            if col not in self.datasets_df.columns:
                missing_cols.append(f"datasets表缺少列: {col}")
        
        if missing_cols:
            raise ValueError("表格缺少必需的列:\n" + "\n".join(missing_cols))
    
    def filter_human_scrna_datasets(self, output_dir: str = '.') -> pd.DataFrame:
        """
        筛选人源单细胞RNA-seq数据集，并保存中间表格
        
        Parameters:
        -----------
        output_dir : str
            输出目录
        
        Returns:
        --------
        pd.DataFrame: 筛选后的datasets表格
        """
        print("开始筛选数据集...")
        print("="*80)
        
        original_count = len(self.datasets_df)
        
        # 1. 筛选人源数据 (Homo sapiens)
        print("\n步骤 1: 筛选人源数据 (Homo sapiens)")
        print("-" * 40)
        
        # 检查organisms列，支持多种格式
        human_mask = self.datasets_df['organisms'].str.contains(
            'Homo sapiens|NCBITaxon:9606|9606',
            case=False,
            na=False,
            regex=True
        )
        
        filtered_df = self.datasets_df[human_mask].copy()
        print(f"  原始数据集: {original_count}")
        print(f"  人源数据集: {len(filtered_df)}")
        print(f"  过滤掉: {original_count - len(filtered_df)}")
        
        if len(filtered_df) == 0:
            print("\n⚠ 警告: 没有找到人源数据集！")
            print("Organisms列的值示例:")
            print(self.datasets_df['organisms'].head(10))
            return pd.DataFrame()
        
        # 2. 筛选单细胞RNA-seq数据
        print("\n步骤 2: 筛选单细胞RNA-seq数据")
        print("-" * 40)
        
        # 单细胞RNA-seq相关的EFO ID
        scrna_efo_ids = [
            'EFO:0009899',  # 10x 3' v2
            'EFO:0009922',  # 10x 3' v3
            'EFO:0009901',  # 10x 3' v1
            'EFO:0030003',  # 10x 3' transcription profiling
            'EFO:0011025',  # 10x 5' v1
            'EFO:0009900',  # 10x 5' v2
            'EFO:0030004',  # 10x 5' transcription profiling
            'EFO:0008931',  # Smart-seq2
            'EFO:0022488',  # Smart-seq3
            'EFO:0008722',  # Drop-seq
            'EFO:0008919',  # Seq-Well
            'EFO:0008780',  # inDrop
            'EFO:0008679',  # CEL-seq
            'EFO:0010010',  # CEL-seq2
            'EFO:0700010',  # TruDrop
            'EFO:0022845',  # modified STRT-seq
            'EFO:0022490',  # ScaleBio single cell RNA sequencing
            'EFO:0700003',  # BD Rhapsody Whole Transcriptome Analysis
        ]
        
        # 构建筛选模式
        scrna_pattern = '|'.join(scrna_efo_ids)
        
        # 包含单细胞RNA-seq技术
        scrna_mask = filtered_df['assays'].str.contains(
            scrna_pattern,
            case=False,
            na=False,
            regex=True
        )
        
        before_scrna = len(filtered_df)
        filtered_df = filtered_df[scrna_mask].copy()
        
        print(f"  筛选前: {before_scrna}")
        print(f"  包含单细胞RNA-seq技术: {len(filtered_df)}")
        print(f"  过滤掉: {before_scrna - len(filtered_df)}")
        
        # 进一步筛选：排除纯空间转录组数据集
        if len(filtered_df) > 0:
            # 检查是否仅包含空间转录组技术
            spatial_only_mask = (
                filtered_df['assays'].str.contains('EFO:0022857', case=False, na=False) &
                ~filtered_df['assays'].str.contains(
                    '|'.join([id for id in scrna_efo_ids if id != 'EFO:0022857']),
                    case=False,
                    na=False,
                    regex=True
                )
            )
            
            spatial_only_count = spatial_only_mask.sum()
            if spatial_only_count > 0:
                print(f"  排除纯空间转录组: {spatial_only_count} 个数据集")
                filtered_df = filtered_df[~spatial_only_mask].copy()
        
        # 显示检测到的技术类型统计
        if len(filtered_df) > 0:
            print("\n  检测到的单细胞RNA-seq技术:")
            
            # 解析assays列，统计各技术
            tech_counts = {}
            tech_names = {
                'EFO:0009899': '10x 3\' v2',
                'EFO:0009922': '10x 3\' v3',
                'EFO:0009901': '10x 3\' v1',
                'EFO:0011025': '10x 5\' v1',
                'EFO:0009900': '10x 5\' v2',
                'EFO:0008931': 'Smart-seq2',
                'EFO:0022488': 'Smart-seq3',
                'EFO:0008722': 'Drop-seq',
                'EFO:0008919': 'Seq-Well',
                'EFO:0008780': 'inDrop',
                'EFO:0008679': 'CEL-seq',
                'EFO:0010010': 'CEL-seq2',
                'EFO:0700010': 'TruDrop',
                'EFO:0022845': 'modified STRT-seq',
                'EFO:0022490': 'ScaleBio',
                'EFO:0700003': 'BD Rhapsody',
                'EFO:0030059': '10x multiome (含RNA)',
            }
            
            for assay_str in filtered_df['assays'].dropna():
                for efo_id, name in tech_names.items():
                    if efo_id in str(assay_str):
                        tech_counts[name] = tech_counts.get(name, 0) + 1
            
            # 按数量排序显示
            for tech, count in sorted(tech_counts.items(), key=lambda x: -x[1]):
                print(f"    - {tech}: {count} 个数据集")
            
            # 显示suspension_type统计
            if 'suspension_types' in filtered_df.columns:
                print("\n  样本类型分布:")
                susp_counts = {}
                for susp_str in filtered_df['suspension_types'].dropna():
                    types = str(susp_str).split(';')
                    for t in types:
                        t = t.strip()
                        if t and t != 'na':
                            susp_counts[t] = susp_counts.get(t, 0) + 1
                
                for susp, count in sorted(susp_counts.items(), key=lambda x: -x[1]):
                    print(f"    - {susp}: {count} 个数据集")
        else:
            print("\n⚠ 警告: 没有找到单细胞RNA-seq数据集！")
            print("Assays列的值示例:")
            print(self.datasets_df['assays'].head(10))
        
        # 3. 合并Collection信息
        if len(filtered_df) > 0:
            print("\n步骤 3: 合并Collection信息")
            print("-" * 40)
            
            # 合并collections表的信息
            merge_cols = [col for col in ['collection_id', 'name', 'contact_name', 
                                          'contact_email', 'published_year', 'doi',
                                          'description', 'authors', 'consortia']
                         if col in self.collections_df.columns]
            
            filtered_df = filtered_df.merge(
                self.collections_df[merge_cols],
                on='collection_id',
                how='left',
                suffixes=('_dataset', '_collection')
            )
            
            print(f"  合并后的列数: {len(filtered_df.columns)}")
            print(f"  涉及的Collections: {filtered_df['collection_id'].nunique()}")
        
        # 4. 保存中间表格
        if len(filtered_df) > 0:
            print("\n步骤 4: 保存中间表格")
            print("-" * 40)
            
            # 创建输出目录
            os.makedirs(output_dir, exist_ok=True)
            
            # 保存clean_datasets
            clean_datasets_file = os.path.join(output_dir, 'clean_datasets.csv')
            # 只保存datasets相关的列（不包括collection的详细信息）
            dataset_cols = [col for col in filtered_df.columns 
                          if not col.endswith('_collection') and col not in merge_cols[1:]]
            clean_datasets = filtered_df[dataset_cols].copy()
            clean_datasets.to_csv(clean_datasets_file, index=False, encoding='utf-8-sig')
            print(f"  ✓ clean_datasets.csv: {clean_datasets_file}")
            print(f"    - 行数: {len(clean_datasets)}")
            print(f"    - 列数: {len(clean_datasets.columns)}")
            
            # 保存clean_collections
            # 筛选出涉及到的collections
            clean_collection_ids = filtered_df['collection_id'].unique()
            clean_collections = self.collections_df[
                self.collections_df['collection_id'].isin(clean_collection_ids)
            ].copy()
            
            clean_collections_file = os.path.join(output_dir, 'clean_collections.csv')
            clean_collections.to_csv(clean_collections_file, index=False, encoding='utf-8-sig')
            print(f"  ✓ clean_collections.csv: {clean_collections_file}")
            print(f"    - 行数: {len(clean_collections)}")
            print(f"    - 列数: {len(clean_collections.columns)}")
        
        # 5. 最终汇总
        print("\n" + "="*80)
        print("筛选结果汇总:")
        print("="*80)
        print(f"✓ 最终保留数据集: {len(filtered_df)}")
        print(f"  来自 {filtered_df['collection_id'].nunique()} 个Collections")
        print(f"  筛选比例: {100*len(filtered_df)/original_count:.1f}%")
        
        if len(filtered_df) > 0:
            # 显示一些统计信息
            print("\n数据集分布:")
            
            # 按组织统计
            if 'tissues' in filtered_df.columns:
                print("\n  Top 5 组织类型:")
                tissue_counts = filtered_df['tissues'].value_counts().head(5)
                for tissue, count in tissue_counts.items():
                    print(f"    - {str(tissue)[:50]}: {count}")
            
            # 按疾病统计
            if 'diseases' in filtered_df.columns:
                print("\n  Top 5 疾病类型:")
                disease_counts = filtered_df['diseases'].value_counts().head(5)
                for disease, count in disease_counts.items():
                    print(f"    - {str(disease)[:50]}: {count}")
            
            # 按年份统计
            if 'published_year' in filtered_df.columns:
                print("\n  按发表年份:")
                year_counts = filtered_df['published_year'].value_counts().sort_index()
                for year, count in year_counts.items():
                    if pd.notna(year):
                        print(f"    - {int(year)}: {count}")
        
        print("="*80 + "\n")
        
        return filtered_df
    
    def prepare_dataset_info(self, filtered_df: pd.DataFrame = None) -> List[Dict]:
        """
        准备dataset信息列表
        
        Parameters:
        -----------
        filtered_df : pd.DataFrame
            筛选后的数据集，如果为None则自动筛选
        
        Returns:
        --------
        List[Dict]: 包含完整信息的dataset列表
        """
        if filtered_df is None:
            filtered_df = self.filter_human_scrna_datasets()
        
        if len(filtered_df) == 0:
            print("⚠ 没有符合条件的数据集")
            return []
        
        # 构建dataset信息列表
        dataset_info_list = []
        for _, row in filtered_df.iterrows():
            info = {
                'dataset_id': row['dataset_id'],
                'collection_id': row['collection_id'],
                'collection_name': row.get('name') or row.get('collection_name'),
                'organisms': row.get('organisms'),
                'diseases': row.get('diseases'),
                'tissues': row.get('tissues'),
                'assays': row.get('assays'),
                'suspension_types': row.get('suspension_types'),
                'contact_name': row.get('contact_name'),
                'contact_email': row.get('contact_email'),
                'published_year': row.get('published_year'),
                'doi': row.get('doi'),
                'description': row.get('description'),
                'authors': row.get('authors'),
            }
            dataset_info_list.append(info)
        
        print(f"✓ 准备了 {len(dataset_info_list)} 个人源单细胞RNA-seq数据集信息\n")
        return dataset_info_list
    
    def get_samples_from_census(self, dataset_info_list: List[Dict] = None, 
                                batch_size: int = 5,
                                save_intermediate: bool = True) -> pd.DataFrame:
        """
        使用CELLxGENE Census API获取Sample层面信息
        
        Parameters:
        -----------
        dataset_info_list : List[Dict]
            Dataset信息列表，如果为None则自动准备
        batch_size : int
            每批处理的dataset数量（建议5-10）
        save_intermediate : bool
            是否保存中间结果（防止中断后重新开始）
        
        Returns:
        --------
        pd.DataFrame: Sample层面的信息表格
        """
        if dataset_info_list is None:
            filtered_df = self.filter_human_scrna_datasets()
            dataset_info_list = self.prepare_dataset_info(filtered_df)
        
        if not dataset_info_list:
            raise ValueError("没有找到符合条件的dataset信息")
        
        all_samples = []
        failed_datasets = []
        processed_count = 0
        
        print(f"开始从CELLxGENE Census提取Sample信息...")
        print(f"总共需要处理 {len(dataset_info_list)} 个人源单细胞RNA-seq数据集")
        print(f"批次大小: {batch_size}")
        print("="*80 + "\n")
        
        # 检查是否有中间结果
        intermediate_file = 'intermediate_samples_human_scrna.csv'
        if save_intermediate and os.path.exists(intermediate_file):
            print(f"发现中间结果文件: {intermediate_file}")
            response = input("是否从中间结果继续？(y/n): ")
            if response.lower() == 'y':
                all_samples.append(pd.read_csv(intermediate_file))
                print(f"✓ 加载了 {len(all_samples[0])} 条已有记录\n")
        
        with cellxgene_census.open_soma(census_version='2024-05-20') as census:
            human = census["census_data"]["homo_sapiens"]
            
            # 使用进度条
            for i, info in enumerate(tqdm(dataset_info_list, desc="处理Datasets")):
                dataset_id = info['dataset_id']
                
                try:
                    # 读取该dataset的obs数据 - 修复列名
                    obs = human.obs.read(
                        column_names=[
                            "dataset_id",
                            "donor_id",
                            "suspension_type",
                            "tissue",
                            "tissue_general",
                            "tissue_ontology_term_id",
                            "assay",
                            "assay_ontology_term_id",
                            "disease",
                            "disease_ontology_term_id",
                            "sex",
                            "sex_ontology_term_id",
                            "development_stage",
                            "development_stage_ontology_term_id",
                            "self_reported_ethnicity",  # 修复: ethnicity -> self_reported_ethnicity
                            "self_reported_ethnicity_ontology_term_id",
                            "cell_type",
                            "cell_type_ontology_term_id",
                            "is_primary_data"
                            # 移除: organism, organism_ontology_term_id (Census中不存在)
                        ],
                        value_filter=f"dataset_id == '{dataset_id}'"
                    ).concat().to_pandas()
                    
                    if len(obs) > 0:
                        # 按donor_id分组，汇总sample信息
                        sample_summary = obs.groupby('donor_id').agg({
                            'dataset_id': 'first',
                            'suspension_type': 'first',
                            'tissue': lambda x: '|'.join(sorted(set(str(v) for v in x if pd.notna(v)))),
                            'tissue_general': lambda x: '|'.join(sorted(set(str(v) for v in x if pd.notna(v)))),
                            'tissue_ontology_term_id': lambda x: '|'.join(sorted(set(str(v) for v in x if pd.notna(v)))),
                            'assay': lambda x: '|'.join(sorted(set(str(v) for v in x if pd.notna(v)))),
                            'assay_ontology_term_id': lambda x: '|'.join(sorted(set(str(v) for v in x if pd.notna(v)))),
                            'disease': lambda x: '|'.join(sorted(set(str(v) for v in x if pd.notna(v)))),
                            'disease_ontology_term_id': lambda x: '|'.join(sorted(set(str(v) for v in x if pd.notna(v)))),
                            'sex': 'first',
                            'sex_ontology_term_id': 'first',
                            'development_stage': 'first',
                            'development_stage_ontology_term_id': 'first',
                            'self_reported_ethnicity': 'first',  # 修复
                            'self_reported_ethnicity_ontology_term_id': 'first',
                            'is_primary_data': 'first',
                            'cell_type': lambda x: len(set(x)),  # 细胞类型数量
                            'cell_type_ontology_term_id': lambda x: len(set(x))  # 细胞类型本体论ID数量
                        }).reset_index()
                        
                        # 添加细胞计数
                        cell_counts = obs.groupby('donor_id').size().reset_index(name='n_cells')
                        sample_summary = sample_summary.merge(cell_counts, on='donor_id')
                        
                        # 添加organism信息（直接设为Homo sapiens，因为已经在homo_sapiens节点中）
                        sample_summary['organism'] = 'Homo sapiens'
                        sample_summary['organism_ontology_term_id'] = 'NCBITaxon:9606'
                        
                        # 添加collection和dataset信息
                        for key in ['collection_id', 'collection_name', 'organisms', 
                                   'diseases', 'tissues', 'assays', 'suspension_types',
                                   'contact_name', 'contact_email', 'published_year', 
                                   'doi', 'description', 'authors']:
                            sample_summary[key] = info.get(key)
                        
                        # 重命名列
                        sample_summary.rename(columns={
                            'donor_id': 'sample_id',
                            'cell_type': 'n_cell_types'
                        }, inplace=True)
                        
                        all_samples.append(sample_summary)
                        processed_count += 1
                        
                        # 每批次后显示进度和保存中间结果
                        if (i + 1) % batch_size == 0:
                            total_samples = sum(len(s) for s in all_samples)
                            print(f"\n进度: {i+1}/{len(dataset_info_list)} datasets")
                            print(f"已提取: {total_samples} 个samples")
                            print(f"成功率: {processed_count}/{i+1} ({100*processed_count/(i+1):.1f}%)")
                            
                            # 保存中间结果
                            if save_intermediate and all_samples:
                                temp_df = pd.concat(all_samples, ignore_index=True)
                                temp_df.to_csv(intermediate_file, index=False)
                                print(f"✓ 中间结果已保存: {intermediate_file}\n")
                    
                    else:
                        tqdm.write(f"⚠ Dataset {dataset_id[:8]}... 没有找到数据")
                        failed_datasets.append({
                            'dataset_id': dataset_id,
                            'collection_name': info.get('collection_name'),
                            'reason': 'No data found in Census'
                        })
                
                except Exception as e:
                    tqdm.write(f"✗ Dataset {dataset_id[:8]}... 处理失败: {str(e)[:50]}")
                    failed_datasets.append({
                        'dataset_id': dataset_id,
                        'collection_name': info.get('collection_name'),
                        'reason': str(e)[:100]
                    })
                
                # 避免请求过快
                if (i + 1) % batch_size == 0:
                    time.sleep(1)
        
        # 合并所有结果
        if all_samples:
            samples_df = pd.concat(all_samples, ignore_index=True)
            
            # 调整列顺序
            column_order = [
                'collection_id', 'collection_name', 'dataset_id',
                'sample_id', 'n_cells', 'n_cell_types',
                'tissue', 'tissue_general', 'tissue_ontology_term_id',
                'assay', 'assay_ontology_term_id',
                'disease', 'disease_ontology_term_id',
                'sex', 'sex_ontology_term_id',
                'development_stage', 'development_stage_ontology_term_id',
                'self_reported_ethnicity', 'self_reported_ethnicity_ontology_term_id',
                'suspension_type', 
                'organism', 'organism_ontology_term_id',
                'is_primary_data',
                'organisms', 'diseases', 'tissues', 'assays', 'suspension_types',
                'contact_name', 'contact_email', 'published_year', 'doi',
                'description', 'authors'
            ]
            
            # 只保留存在的列
            available_columns = [col for col in column_order if col in samples_df.columns]
            other_columns = [col for col in samples_df.columns if col not in column_order]
            samples_df = samples_df[available_columns + other_columns]
            
            # 删除中间文件
            if save_intermediate and os.path.exists(intermediate_file):
                os.remove(intermediate_file)
                print(f"✓ 已删除中间结果文件")
            
            print("\n" + "="*80)
            print("提取完成!")
            print("="*80)
            print(f"✓ 成功提取了 {len(samples_df)} 个人源单细胞RNA-seq Samples")
            print(f"  来自 {samples_df['dataset_id'].nunique()} 个Datasets")
            print(f"  来自 {samples_df['collection_id'].nunique()} 个Collections")
            print(f"  总细胞数: {samples_df['n_cells'].sum():,}")
            print(f"  平均每个sample细胞数: {samples_df['n_cells'].mean():.0f}")
            print(f"  中位数每个sample细胞数: {samples_df['n_cells'].median():.0f}")
            
            if failed_datasets:
                print(f"\n⚠ {len(failed_datasets)} 个Datasets处理失败")
                failed_df = pd.DataFrame(failed_datasets)
                failed_file = 'failed_datasets_human_scrna.csv'
                failed_df.to_csv(failed_file, index=False)
                print(f"  失败列表已保存到: {failed_file}")
            
            print("="*80 + "\n")
            
            return samples_df
        
        else:
            print("\n✗ 没有提取到任何Sample数据")
            if failed_datasets:
                print(f"所有 {len(failed_datasets)} 个Datasets都失败了")
                pd.DataFrame(failed_datasets).to_csv('failed_datasets_human_scrna.csv', index=False)
            return pd.DataFrame()
    
    def generate_summary_statistics(self, samples_df: pd.DataFrame) -> pd.DataFrame:
        """
        生成汇总统计信息
        
        Parameters:
        -----------
        samples_df : pd.DataFrame
            Sample数据表格
        
        Returns:
        --------
        pd.DataFrame: 统计摘要表格
        """
        summary_data = []
        
        # 总体统计
        summary_data.append({
            'Category': 'Overall',
            'Metric': 'Total Samples',
            'Count': len(samples_df),
            'Percentage': 100.0
        })
        summary_data.append({
            'Category': 'Overall',
            'Metric': 'Total Cells',
            'Count': int(samples_df['n_cells'].sum()),
            'Percentage': 100.0
        })
        summary_data.append({
            'Category': 'Overall',
            'Metric': 'Total Collections',
            'Count': samples_df['collection_id'].nunique(),
            'Percentage': 100.0
        })
        summary_data.append({
            'Category': 'Overall',
            'Metric': 'Total Datasets',
            'Count': samples_df['dataset_id'].nunique(),
            'Percentage': 100.0
        })
        
        # 按各个维度统计
        for col, category in [
            ('tissue_general', 'Tissue'),
            ('disease', 'Disease'),
            ('assay', 'Assay'),
            ('sex', 'Sex'),
            ('development_stage', 'Development Stage'),
            ('suspension_type', 'Suspension Type')
        ]:
            if col in samples_df.columns:
                value_counts = samples_df[col].value_counts()
                for value, count in value_counts.items():
                    summary_data.append({
                        'Category': category,
                        'Metric': str(value)[:50],  # 限制长度
                        'Count': count,
                        'Percentage': 100 * count / len(samples_df)
                    })
        
        summary_df = pd.DataFrame(summary_data)
        return summary_df
    
    def save_results(self, samples_df: pd.DataFrame, output_dir: str = '.', 
                    output_prefix: str = 'cellxgene_human_scrna_samples'):
        """
        保存结果到文件
        
        Parameters:
        -----------
        samples_df : pd.DataFrame
            Sample数据表格
        output_dir : str
            输出目录
        output_prefix : str
            输出文件前缀
        """
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        print("保存结果文件...")
        print("="*80)
        
        # 1. 保存主数据CSV
        csv_file = os.path.join(output_dir, f'{output_prefix}.csv')
        samples_df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        print(f"✓ 主数据CSV: {csv_file}")
        print(f"  行数: {len(samples_df)}")
        print(f"  列数: {len(samples_df.columns)}")
        
        # 2. 保存Excel（多sheet）
        try:
            excel_file = os.path.join(output_dir, f'{output_prefix}.xlsx')
            with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                # Sheet 1: 主数据
                samples_df.to_excel(writer, sheet_name='Samples', index=False)
                
                # Sheet 2: 统计摘要
                summary_df = self.generate_summary_statistics(samples_df)
                summary_df.to_excel(writer, sheet_name='Summary Statistics', index=False)
                
                # Sheet 3: 按Collection汇总
                collection_summary = samples_df.groupby(['collection_id', 'collection_name']).agg({
                    'sample_id': 'count',
                    'n_cells': 'sum',
                    'dataset_id': 'nunique'
                }).reset_index()
                collection_summary.columns = ['Collection ID', 'Collection Name', 
                                             'N Samples', 'Total Cells', 'N Datasets']
                collection_summary = collection_summary.sort_values('N Samples', ascending=False)
                collection_summary.to_excel(writer, sheet_name='By Collection', index=False)
                
                # Sheet 4: 按Dataset汇总
                dataset_summary = samples_df.groupby(['collection_name', 'dataset_id']).agg({
                    'sample_id': 'count',
                    'n_cells': 'sum'
                }).reset_index()
                dataset_summary.columns = ['Collection', 'Dataset ID', 'N Samples', 'Total Cells']
                dataset_summary = dataset_summary.sort_values('N Samples', ascending=False)
                dataset_summary.to_excel(writer, sheet_name='By Dataset', index=False)
                
                # Sheet 5: 按Tissue汇总
                if 'tissue_general' in samples_df.columns:
                    tissue_summary = samples_df.groupby('tissue_general').agg({
                        'sample_id': 'count',
                        'n_cells': 'sum',
                        'collection_id': 'nunique'
                    }).reset_index()
                    tissue_summary.columns = ['Tissue', 'N Samples', 'Total Cells', 'N Collections']
                    tissue_summary = tissue_summary.sort_values('N Samples', ascending=False)
                    tissue_summary.to_excel(writer, sheet_name='By Tissue', index=False)
                
                # Sheet 6: 按Assay汇总
                if 'assay' in samples_df.columns:
                    assay_summary = samples_df.groupby('assay').agg({
                        'sample_id': 'count',
                        'n_cells': 'sum',
                        'collection_id': 'nunique'
                    }).reset_index()
                    assay_summary.columns = ['Assay', 'N Samples', 'Total Cells', 'N Collections']
                    assay_summary = assay_summary.sort_values('N Samples', ascending=False)
                    assay_summary.to_excel(writer, sheet_name='By Assay', index=False)
            
            print(f"✓ Excel文件: {excel_file}")
            print(f"  包含 6 个sheets")
        except Exception as e:
            print(f"⚠ Excel保存失败: {e}")
        
        # 3. 保存JSON
        json_file = os.path.join(output_dir, f'{output_prefix}.json')
        samples_df.to_json(json_file, orient='records', indent=2, force_ascii=False)
        print(f"✓ JSON文件: {json_file}")
        
        # 4. 保存详细统计信息
        stats_file = os.path.join(output_dir, f'{output_prefix}_statistics.txt')
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("CellxGene 人源单细胞RNA-seq Sample数据统计报告\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"总览:\n")
            f.write(f"  总Samples数: {len(samples_df):,}\n")
            f.write(f"  总细胞数: {samples_df['n_cells'].sum():,}\n")
            f.write(f"  总Collections数: {samples_df['collection_id'].nunique()}\n")
            f.write(f"  总Datasets数: {samples_df['dataset_id'].nunique()}\n\n")
            
            f.write(f"细胞统计:\n")
            f.write(f"  平均每sample细胞数: {samples_df['n_cells'].mean():.0f}\n")
            f.write(f"  中位数每sample细胞数: {samples_df['n_cells'].median():.0f}\n")
            f.write(f"  最小细胞数: {samples_df['n_cells'].min()}\n")
            f.write(f"  最大细胞数: {samples_df['n_cells'].max()}\n\n")
            
            # 各维度Top 10
            for col, title in [
                ('tissue_general', 'Top 10 Tissues'),
                ('disease', 'Top 10 Diseases'),
                ('assay', 'Top 10 Assays'),
                ('collection_name', 'Top 10 Collections'),
                ('sex', 'Sex Distribution'),
                ('development_stage', 'Development Stage Distribution')
            ]:
                if col in samples_df.columns:
                    f.write(f"{title}:\n")
                    top_values = samples_df[col].value_counts().head(10)
                    for value, count in top_values.items():
                        f.write(f"  {value}: {count} ({100*count/len(samples_df):.1f}%)\n")
                    f.write("\n")
        
        print(f"✓ 统计报告: {stats_file}")
        print("="*80 + "\n")


def main():
    """
    主函数 - 只处理人源单细胞RNA-seq数据
    """
    # ============================================================
    # 配置部分 - 根据你的路径修改
    # ============================================================
    
    # 输入文件路径
    collections_file = '/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/processed/collections_table.csv'
    datasets_file = '/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/processed/datasets_table.csv'
    
    # 输出目录和文件前缀
    output_dir = '/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/processed'
    output_prefix = 'cellxgene_human_scrna_samples'
    
    # 处理参数
    batch_size = 5  # 每批处理5个dataset，可根据网络情况调整
    
    # ============================================================
    
    try:
        # 1. 初始化提取器
        extractor = CellxGeneSampleExtractor(
            collections_file=collections_file,
            datasets_file=datasets_file
        )
        
        # 2. 筛选人源单细胞RNA-seq数据集（会自动保存clean_collections和clean_datasets）
        filtered_df = extractor.filter_human_scrna_datasets(output_dir=output_dir)
        
        if len(filtered_df) == 0:
            print("\n⚠ 警告: 没有找到符合条件的数据集（人源 + 单细胞RNA-seq）")
            print("请检查你的数据集中是否包含:")
            print("  - Homo sapiens (人源)")
            print("  - 单细胞RNA-seq技术 (EFO本体论ID)")
            return
        
        # 3. 提取Sample信息
        samples_df = extractor.get_samples_from_census(
            batch_size=batch_size,
            save_intermediate=True  # 保存中间结果,防止中断
        )
        
        if len(samples_df) > 0:
            # 4. 显示数据预览
            print("数据预览 (前10条):")
            print("="*80)
            print(samples_df.head(10).to_string())
            print("\n")
            
            # 5. 保存结果
            extractor.save_results(
                samples_df, 
                output_dir=output_dir,
                output_prefix=output_prefix
            )
            
            print("✓ 所有操作完成!")
            print(f"\n输出文件位置: {output_dir}/")
            print(f"\n中间文件:")
            print(f"  - clean_collections.csv (筛选后的Collections表)")
            print(f"  - clean_datasets.csv (筛选后的Datasets表)")
            print(f"\n最终结果文件:")
            print(f"  - {output_prefix}.csv")
            print(f"  - {output_prefix}.xlsx")
            print(f"  - {output_prefix}.json")
            print(f"  - {output_prefix}_statistics.txt")
        
        else:
            print("\n✗ 没有提取到任何数据")
    
    except FileNotFoundError as e:
        print(f"\n✗ 错误: 文件未找到")
        print(f"  {e}")
        print("\n请检查文件路径是否正确:")
        print(f"  Collections: {collections_file}")
        print(f"  Datasets: {datasets_file}")
    
    except Exception as e:
        print(f"\n✗ 发生错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()