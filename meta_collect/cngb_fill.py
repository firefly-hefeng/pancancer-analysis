import pandas as pd
from pathlib import Path
import json
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CNGBComprehensiveProcessor:
    """CNGB metadata全面处理器 - 生成超级总表"""
    
    def __init__(self, metadata_dir, output_dir):
        self.metadata_dir = Path(metadata_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 定义超级总表的字段映射（智能映射所有可能的列名变体）
        self.super_mapping = {
            # ============ 项目标识 ============
            'project_id': ['project_accession', 'project_id'],
            
            # ============ 样本标识 ============
            'sample_id': ['sample_accession', 'sample_id', 'biosample_accession', 'biosample'],
            'sample_name': ['sample_name', 'sample_alias', 'sample_title'],
            
            # ============ 实验标识 ============
            'experiment_id': ['experiment_accession', 'experiment_id', 'exp_id'],
            'run_id': ['run_accession', 'run_id', 'srr', 'run'],
            'experiment_title': ['experiment_title', 'experiment_name', 'exp_title'],
            
            # ============ 单细胞标识 ============
            'singlecell_id': ['single-cell_accession', 'singlecell_accession', 'sc_accession'],
            
            # ============ 生物学信息 ============
            'organism': ['organism', 'species', 'scientific_name'],
            'tax_id': ['tax_id', 'taxid', 'taxonomy_id', 'ncbi_taxid'],
            'strain': ['strain', 'isolate', 'breed'],
            'tissue': ['tissue', 'tissue_type', 'organ', 'source_tissue'],
            'cell_line': ['cell_line', 'cellline', 'cell_culture'],
            'cell_type': ['cell_type', 'celltype', 'cell_category'],
            'cell_subtype': ['cell_subtype', 'cell_subset', 'subtype'],
            
            # ============ 发育和生理信息 ============
            'dev_stage': ['dev_stage', 'developmental_stage', 'development_stage'],
            'age': ['age', 'age_at_collection', 'age_at_diagnosis'],
            'sex': ['sex', 'gender'],
            'ethnicity': ['ethnicity', 'ethnic_group', 'race'],
            'population': ['population', 'race', 'ancestry'],
            
            # ============ 疾病信息 ============
            'disease': ['disease', 'disease_state', 'condition', 'diagnosis'],
            'disease_stage': ['disease_stage', 'stage', 'tumor_stage', 'cancer_stage', 'clinical_stage'],
            'phenotype': ['phenotype', 'phenotypic_feature', 'clinical_phenotype'],
            'health_state': ['health_state', 'health_status', 'disease_status'],
            
            # ============ 样本处理信息 ============
            'treatment': ['treatment', 'drug_treatment', 'therapy'],
            'culture_collection': ['culture_collection', 'culture_condition'],
            'biomaterial_provider': ['biomaterial_provider', 'provider', 'source_institute'],
            'karyotype': ['karyotype', 'chromosome'],
            
            # ============ 测序库信息 ============
            'library_name': ['library_name', 'library_id', 'lib_name'],
            'library_strategy': ['library_strategy', 'assay_type', 'sequencing_strategy'],
            'library_source': ['library_source', 'source_material'],
            'library_selection': ['library_selection', 'selection_method'],
            'library_layout': ['library_layout', 'layout', 'read_type'],
            
            # ============ 测序平台信息 ============
            'platform': ['platform', 'sequencing_platform', 'instrument_platform'],
            'instrument_model': ['instrument_model', 'instrument', 'sequencer', 'machine'],
            
            # ============ 实验设计 ============
            'design_description': ['design_description', 'experimental_design', 'protocol', 'description'],
            
            # ============ 文件信息 ============
            'file_type': ['file_type', 'data_type', 'format'],
            
            # FastQ R1
            'fastq_r1': ['file_name', 'file1_name', 'read1_file', 'fastq_1', 'r1_file', 'read1', 'fq1'],
            'fastq_r1_md5': ['file_md5', 'file1_md5', 'read1_md5', 'md5_1', 'r1_md5'],
            
            # FastQ R2
            'fastq_r2': ['file2_name', 'read2_file', 'fastq_2', 'r2_file', 'read2', 'fq2'],
            'fastq_r2_md5': ['file2_md5', 'read2_md5', 'md5_2', 'r2_md5'],
            
            # FastQ R3 (如果有)
            'fastq_r3': ['file3_name', 'read3_file', 'fastq_3', 'r3_file', 'read3'],
            'fastq_r3_md5': ['file3_md5', 'read3_md5', 'md5_3', 'r3_md5'],
            
            # ============ 单细胞数据文件 ============
            'expression_matrix': ['gene_expression_file_name', 'expression_file', 'matrix_file', 'count_matrix', 'expression_matrix'],
            'expression_matrix_md5': ['gene_expression_file_md5', 'matrix_md5', 'expression_md5'],
            'expression_description': ['gene_expression_file_description', 'expression_description', 'matrix_description'],
            'expression_axis_label': ['expression_axis_label', 'value_type', 'count_type'],
            'expression_file_type': ['gene_expression_file_type', 'matrix_type', 'file_format'],
            
            'features_file': ['genes_file_name', 'features_file', 'gene_file', 'feature_file'],
            'features_file_md5': ['genes_file_MD5', 'genes_file_md5', 'features_md5', 'gene_md5'],
            
            'barcodes_file': ['barcodes_file_name', 'barcode_file', 'cells_file'],
            'barcodes_file_md5': ['barcodes_file_md5', 'barcode_md5', 'cells_md5'],
            
            # ============ 其他补充信息 ============
            'sample_type': ['sample_type', 'material_type', 'sample_category'],
            'sample_title': ['sample_title'],
            'additional_info': ['description', 'comment', 'note', 'remarks'],
        }
        
        # 统计信息
        self.stats = {
            'total_projects': 0,
            'total_records': 0,
            'experiment_count': 0,
            'sample_count': 0,
            'singlecell_count': 0,
            'merged_count': 0,
            'field_usage': defaultdict(int),
            'missing_fields': defaultdict(int),
        }
    
    def parse_file(self, file_path):
        """解析文件"""
        try:
            ext = file_path.suffix.lower()
            
            if ext in ['.tsv', '.txt']:
                for sep in ['\t', ',', '|']:
                    try:
                        df = pd.read_csv(file_path, sep=sep, nrows=5)
                        if len(df.columns) > 1:
                            df_full = pd.read_csv(file_path, sep=sep, dtype=str, keep_default_na=False)
                            # 替换空字符串为None
                            df_full = df_full.replace('', None)
                            return df_full
                    except:
                        continue
            
            elif ext == '.csv':
                df = pd.read_csv(file_path, dtype=str, keep_default_na=False)
                df = df.replace('', None)
                return df
            
            elif ext == '.xlsx':
                df = pd.read_excel(file_path, dtype=str, keep_default_na=False)
                df = df.replace('', None)
                return df
            
            return None
            
        except Exception as e:
            logger.error(f"解析文件失败 {file_path.name}: {str(e)}")
            return None
    
    def smart_map_columns(self, df, source_type):
        """智能映射列名 - 尝试所有可能的映射"""
        if df is None or df.empty:
            return None
        
        mapped_data = {}
        original_columns = df.columns.tolist()
        
        # 为每个目标字段寻找最佳源字段
        for target_field, possible_sources in self.super_mapping.items():
            value = None
            matched_source = None
            
            # 尝试所有可能的源字段（优先级从前到后）
            for source_field in possible_sources:
                if source_field in original_columns:
                    # 找到第一个非空的值
                    temp_values = df[source_field].dropna()
                    if len(temp_values) > 0:
                        value = df[source_field]
                        matched_source = source_field
                        self.stats['field_usage'][f"{target_field} <- {source_field}"] += 1
                        break
            
            if value is not None:
                mapped_data[target_field] = value
            else:
                mapped_data[target_field] = None
                self.stats['missing_fields'][target_field] += 1
        
        # 添加数据来源标记
        mapped_data['data_source'] = source_type
        
        # 保留未映射的原始列（作为补充信息）
        mapped_fields = set()
        for sources in self.super_mapping.values():
            mapped_fields.update(sources)
        
        for col in original_columns:
            if col not in mapped_fields:
                mapped_data[f'extra_{col}'] = df[col]
                logger.debug(f"保留额外字段: extra_{col}")
        
        result_df = pd.DataFrame(mapped_data)
        return result_df
    
    def process_project(self, project_folder):
        """处理单个项目的所有文件"""
        project_id = project_folder.name
        logger.info(f"处理项目: {project_id}")
        
        project_dfs = {
            'experiments': None,
            'samples': None,
            'singlecells': None,
        }
        
        try:
            # 1. 处理experiment文件
            experiment_files = list(project_folder.glob('*experiment*.tsv')) + \
                             list(project_folder.glob('*experiment*.csv'))
            
            if experiment_files:
                df = self.parse_file(experiment_files[0])
                if df is not None and not df.empty:
                    mapped_df = self.smart_map_columns(df, 'experiment')
                    if mapped_df is not None:
                        project_dfs['experiments'] = mapped_df
                        self.stats['experiment_count'] += len(mapped_df)
                        logger.info(f"  ✓ Experiment: {len(mapped_df)} 条记录")
            
            # 2. 处理sample文件
            sample_files = list(project_folder.glob('*sample*.tsv')) + \
                          list(project_folder.glob('*sample*.csv'))
            
            if sample_files:
                df = self.parse_file(sample_files[0])
                if df is not None and not df.empty:
                    mapped_df = self.smart_map_columns(df, 'sample')
                    if mapped_df is not None:
                        project_dfs['samples'] = mapped_df
                        self.stats['sample_count'] += len(mapped_df)
                        logger.info(f"  ✓ Sample: {len(mapped_df)} 条记录")
            
            # 3. 处理singlecell文件
            singlecell_files = list(project_folder.glob('*singlecell*.tsv')) + \
                             list(project_folder.glob('*singlecell*.csv'))
            
            if singlecell_files:
                df = self.parse_file(singlecell_files[0])
                if df is not None and not df.empty:
                    mapped_df = self.smart_map_columns(df, 'singlecell')
                    if mapped_df is not None:
                        project_dfs['singlecells'] = mapped_df
                        self.stats['singlecell_count'] += len(mapped_df)
                        logger.info(f"  ✓ SingleCell: {len(mapped_df)} 条记录")
            
            # 4. 合并项目内的所有数据
            merged_df = self.merge_project_data(project_dfs, project_id)
            
            return merged_df
            
        except Exception as e:
            logger.error(f"处理项目失败 {project_id}: {str(e)}")
            return None
    
    def merge_project_data(self, project_dfs, project_id):
        """合并单个项目内的experiment、sample、singlecell数据"""
        
        # 策略：以experiment为主表，左连接sample和singlecell
        base_df = project_dfs['experiments']
        
        if base_df is None or base_df.empty:
            # 如果没有experiment，看是否有sample或singlecell
            if project_dfs['samples'] is not None:
                base_df = project_dfs['samples']
            elif project_dfs['singlecells'] is not None:
                base_df = project_dfs['singlecells']
            else:
                return None
        
        merged = base_df.copy()
        
        # 合并sample信息（通过sample_id）
        if project_dfs['samples'] is not None and not project_dfs['samples'].empty:
            samples = project_dfs['samples']
            
            # 找到连接键
            if 'sample_id' in merged.columns and 'sample_id' in samples.columns:
                # 去除samples中与merged重复的列（除了连接键）
                overlap_cols = set(merged.columns) & set(samples.columns)
                overlap_cols.discard('sample_id')
                
                # 对重复列进行重命名策略
                samples_to_merge = samples.copy()
                for col in overlap_cols:
                    # 如果experiment中的值为空，用sample中的值填充
                    merged[col] = merged[col].fillna(samples_to_merge[col])
                
                # 添加sample中独有的列
                unique_sample_cols = set(samples.columns) - set(merged.columns)
                for col in unique_sample_cols:
                    merged = merged.merge(
                        samples[['sample_id', col]],
                        on='sample_id',
                        how='left'
                    )
        
        # 合并singlecell信息（通过project_id）
        if project_dfs['singlecells'] is not None and not project_dfs['singlecells'].empty:
            singlecells = project_dfs['singlecells']
            
            # 单细胞数据通常是一个项目一条或几条记录
            # 我们将其信息广播到所有experiment记录上
            if len(singlecells) > 0:
                # 取第一条单细胞记录（大多数情况只有一条）
                sc_row = singlecells.iloc[0]
                
                # 添加单细胞特有的字段
                sc_fields = ['singlecell_id', 'expression_matrix', 'expression_matrix_md5',
                           'expression_description', 'expression_axis_label', 'expression_file_type',
                           'features_file', 'features_file_md5', 'barcodes_file', 'barcodes_file_md5']
                
                for field in sc_fields:
                    if field in sc_row and pd.notna(sc_row[field]):
                        # 如果merged中没有这个字段或为空，则添加
                        if field not in merged.columns:
                            merged[field] = sc_row[field]
                        else:
                            merged[field] = merged[field].fillna(sc_row[field])
        
        # 确保有project_id
        if 'project_id' not in merged.columns or merged['project_id'].isna().all():
            merged['project_id'] = project_id
        else:
            merged['project_id'] = merged['project_id'].fillna(project_id)
        
        logger.info(f"  → 合并后: {len(merged)} 条记录")
        self.stats['merged_count'] += len(merged)
        
        return merged
    
    def create_super_table(self):
        """创建超级总表"""
        logger.info("="*80)
        logger.info("开始创建超级总表...")
        logger.info("="*80)
        
        # 获取所有项目文件夹
        project_folders = sorted([
            d for d in self.metadata_dir.iterdir() 
            if d.is_dir() and d.name.startswith('CNP')
        ])
        
        self.stats['total_projects'] = len(project_folders)
        logger.info(f"发现 {len(project_folders)} 个项目")
        
        # 处理每个项目并收集数据
        all_project_data = []
        
        for i, project_folder in enumerate(project_folders, 1):
            logger.info(f"\n进度: {i}/{len(project_folders)}")
            merged_df = self.process_project(project_folder)
            
            if merged_df is not None and not merged_df.empty:
                all_project_data.append(merged_df)
        
        # 合并所有项目的数据
        if not all_project_data:
            logger.error("没有数据可以合并!")
            return None
        
        logger.info("\n" + "="*80)
        logger.info("合并所有项目数据...")
        
        super_table = pd.concat(all_project_data, ignore_index=True, sort=False)
        self.stats['total_records'] = len(super_table)
        
        logger.info(f"超级总表记录数: {len(super_table)}")
        logger.info(f"超级总表字段数: {len(super_table.columns)}")
        
        # 保存超级总表
        output_file = self.output_dir / "CNGB_super_table.tsv"
        super_table.to_csv(output_file, sep='\t', index=False)
        super_table.to_csv(output_file.with_suffix('.csv'), index=False)
        logger.info(f"✓ 超级总表已保存: {output_file}")
        
        # 生成字段说明
        self.generate_field_description(super_table)
        
        # 生成数据质量报告
        self.generate_quality_report(super_table)
        
        # 生成统计报告
        self.generate_statistics_report(super_table)
        
        return super_table
    
    def generate_field_description(self, super_table):
        """生成字段说明文档"""
        logger.info("生成字段说明...")
        
        field_info = []
        
        for col in super_table.columns:
            non_null_count = super_table[col].notna().sum()
            non_null_pct = non_null_count / len(super_table) * 100
            
            # 获取示例值
            sample_values = super_table[col].dropna().head(3).tolist()
            sample_str = ', '.join([str(v)[:50] for v in sample_values])
            
            # 获取唯一值数量
            unique_count = super_table[col].nunique()
            
            field_info.append({
                'field_name': col,
                'non_null_count': non_null_count,
                'non_null_percentage': f"{non_null_pct:.2f}%",
                'unique_values': unique_count,
                'sample_values': sample_str,
            })
        
        field_df = pd.DataFrame(field_info)
        field_df = field_df.sort_values('non_null_count', ascending=False)
        
        output_file = self.output_dir / "field_description.tsv"
        field_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"✓ 字段说明已保存: {output_file}")
        
        # 打印关键字段信息
        print("\n" + "="*80)
        print("关键字段统计（非空率 > 50%）")
        print("="*80)
        
        important_fields = field_df[field_df['non_null_count'] > len(super_table) * 0.5]
        for _, row in important_fields.iterrows():
            print(f"{row['field_name']:30s} | {row['non_null_percentage']:>8s} | 唯一值: {row['unique_values']:>6d} | 示例: {row['sample_values'][:50]}")
    
    def generate_quality_report(self, super_table):
        """生成数据质量报告"""
        logger.info("生成数据质量报告...")
        
        quality_report = {
            'total_records': len(super_table),
            'total_fields': len(super_table.columns),
            
            # 关键字段完整性
            'key_fields_completeness': {},
            
            # 数据来源统计
            'data_source_distribution': super_table['data_source'].value_counts().to_dict() if 'data_source' in super_table.columns else {},
            
            # 项目统计
            'projects': {
                'total_projects': super_table['project_id'].nunique() if 'project_id' in super_table.columns else 0,
                'project_list': super_table['project_id'].unique().tolist() if 'project_id' in super_table.columns else [],
            },
            
            # 生物体统计
            'organisms': super_table['organism'].value_counts().head(10).to_dict() if 'organism' in super_table.columns else {},
            
            # 组织统计
            'tissues': super_table['tissue'].value_counts().head(20).to_dict() if 'tissue' in super_table.columns else {},
            
            # 疾病统计
            'diseases': super_table['disease'].value_counts().head(20).to_dict() if 'disease' in super_table.columns else {},
            
            # 测序平台统计
            'platforms': super_table['platform'].value_counts().to_dict() if 'platform' in super_table.columns else {},
            
            # 文库策略统计
            'library_strategies': super_table['library_strategy'].value_counts().to_dict() if 'library_strategy' in super_table.columns else {},
        }
        
        # 关键字段完整性检查
        key_fields = ['project_id', 'sample_id', 'experiment_id', 'organism', 'tissue', 
                     'library_strategy', 'platform', 'fastq_r1']
        
        for field in key_fields:
            if field in super_table.columns:
                non_null = super_table[field].notna().sum()
                quality_report['key_fields_completeness'][field] = {
                    'count': int(non_null),
                    'percentage': f"{non_null / len(super_table) * 100:.2f}%"
                }
        
        # 保存报告
        output_file = self.output_dir / "quality_report.json"
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(quality_report, f, indent=2, ensure_ascii=False)
        
        logger.info(f"✓ 质量报告已保存: {output_file}")
        
        # 打印关键统计
        print("\n" + "="*80)
        print("数据质量报告")
        print("="*80)
        print(f"总记录数: {quality_report['total_records']}")
        print(f"总字段数: {quality_report['total_fields']}")
        print(f"涉及项目: {quality_report['projects']['total_projects']}")
        
        if quality_report['organisms']:
            print(f"\n主要生物体:")
            for org, count in list(quality_report['organisms'].items())[:5]:
                print(f"  - {org}: {count}")
        
        if quality_report['platforms']:
            print(f"\n测序平台:")
            for platform, count in quality_report['platforms'].items():
                print(f"  - {platform}: {count}")
    
    def generate_statistics_report(self, super_table):
        """生成处理统计报告"""
        logger.info("生成统计报告...")
        
        stats_report = {
            'processing_stats': {
                'total_projects': self.stats['total_projects'],
                'total_records': self.stats['total_records'],
                'experiment_records': self.stats['experiment_count'],
                'sample_records': self.stats['sample_count'],
                'singlecell_records': self.stats['singlecell_count'],
                'merged_records': self.stats['merged_count'],
            },
            'field_mapping_usage': dict(sorted(self.stats['field_usage'].items(), 
                                              key=lambda x: x[1], reverse=True)),
            'missing_fields_frequency': dict(sorted(self.stats['missing_fields'].items(), 
                                                   key=lambda x: x[1], reverse=True)),
        }
        
        output_file = self.output_dir / "processing_statistics.json"
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(stats_report, f, indent=2, ensure_ascii=False)
        
        logger.info(f"✓ 统计报告已保存: {output_file}")
        
        # 打印字段映射使用情况
        print("\n" + "="*80)
        print("字段映射使用统计 (Top 30)")
        print("="*80)
        for mapping, count in list(stats_report['field_mapping_usage'].items())[:30]:
            print(f"{mapping:60s} : {count:>5d} 次")


def main():
    # 配置路径
    metadata_dir = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/cngb/processed/meta_files"
    output_dir = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/cngb/processed/comprehensive_table"
    
    # 创建处理器并运行
    processor = CNGBComprehensiveProcessor(metadata_dir, output_dir)
    super_table = processor.create_super_table()
    
    
    if super_table is not None:
        logger.info("\n" + "="*80)
        logger.info("✓ 超级总表创建完成!")
        logger.info("="*80)
        logger.info(f"总记录数: {len(super_table)}")
        logger.info(f"总字段数: {len(super_table.columns)}")
        logger.info(f"输出目录: {output_dir}")


if __name__ == "__main__":
    main()