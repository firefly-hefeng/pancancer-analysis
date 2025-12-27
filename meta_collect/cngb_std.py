import pandas as pd
from pathlib import Path
import logging
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CNGBStandardizer:
    """将CNGB超级总表转换为标准格式"""
    
    def __init__(self, super_table_path, output_path):
        self.super_table_path = Path(super_table_path)
        self.output_path = Path(output_path)
        
    def load_super_table(self):
        """加载超级总表"""
        logger.info(f"加载超级总表: {self.super_table_path}")
        
        if self.super_table_path.suffix == '.tsv':
            df = pd.read_csv(self.super_table_path, sep='\t', dtype=str, keep_default_na=False)
        else:
            df = pd.read_csv(self.super_table_path, dtype=str, keep_default_na=False)
        
        df = df.replace('', None)
        logger.info(f"加载了 {len(df)} 条记录，{len(df.columns)} 个字段")
        return df
    
    def infer_disease_general(self, disease):
        """推断疾病大类"""
        if pd.isna(disease) or disease == '':
            return None
        
        disease_lower = str(disease).lower()
        
        # 定义疾病大类映射规则
        disease_categories = {
            'Cancer': ['cancer', 'carcinoma', 'tumor', 'tumour', 'melanoma', 'leukemia', 
                      'lymphoma', 'sarcoma', 'glioma', 'blastoma', 'myeloma', 'malignant',
                      'metastatic', 'neoplasm'],
            'Immune/Inflammatory': ['inflammation', 'inflammatory', 'immune', 'arthritis', 
                                   'lupus', 'allergy', 'asthma', 'colitis'],
            'Neurological': ['alzheimer', 'parkinson', 'dementia', 'neurological', 
                           'neuropathy', 'epilepsy', 'stroke', 'brain'],
            'Cardiovascular': ['cardiovascular', 'cardiac', 'heart', 'hypertension', 
                             'atherosclerosis', 'myocardial'],
            'Metabolic': ['diabetes', 'obesity', 'metabolic', 'syndrome'],
            'Infectious': ['infection', 'viral', 'bacterial', 'covid', 'sars', 'hiv', 
                         'hepatitis', 'tuberculosis', 'sepsis'],
            'Developmental': ['developmental', 'congenital', 'birth defect'],
            'Normal/Healthy': ['normal', 'healthy', 'control', 'wild type', 'wild-type'],
        }
        
        for category, keywords in disease_categories.items():
            for keyword in keywords:
                if keyword in disease_lower:
                    return category
        
        return 'Other'
    
    def determine_file_type(self, row):
        """
        判断文件类型：列出该项目/样本在CNGB中存储的所有可能的数据类型
        CNGB主要存储原始测序数据，不直接提供处理后的矩阵
        """
        file_types = []
        
        # FastQ原始测序数据（最常见）
        if pd.notna(row.get('fastq_r1')) or pd.notna(row.get('fastq_r2')):
            file_types.append('FastQ')
        
        # Run数据（即使没有直接的FastQ文件名，有run_id也说明有原始数据）
        elif pd.notna(row.get('run_id')):
            file_types.append('Run')
        
        # BAM/SAM比对数据（如果CNGB存储了比对结果）
        if pd.notna(row.get('bam_file')) or 'bam' in str(row.get('fastq_r1', '')).lower():
            file_types.append('BAM')
        
        # VCF变异数据（如果存储了变异calling结果）
        if pd.notna(row.get('vcf_file')) or 'vcf' in str(row.get('fastq_r1', '')).lower():
            file_types.append('VCF')
        
        # 10X Genomics格式（features + barcodes + matrix，虽然CNGB主要存原始数据）
        if pd.notna(row.get('features_file')) and pd.notna(row.get('barcodes_file')):
            file_types.append('10X')
        
        # 表达矩阵（如果CNGB恰好存储了处理后的矩阵，但这很少见）
        if pd.notna(row.get('expression_matrix')):
            file_types.append('Matrix')
        
        # 单细胞数据（如果有单细胞ID，说明是单细胞项目）
        if pd.notna(row.get('singlecell_id')):
            if 'Single-cell' not in file_types:
                file_types.append('Single-cell')
        
        # 如果什么都没有，返回None
        if not file_types:
            return None
        
        # 用分号分隔多个文件类型
        return ';'.join(file_types)
    
    def determine_raw_exist(self, row):
        """判断原始数据是否存在（FastQ/Run等）"""
        # 有FastQ文件、Run ID或Experiment ID都表示有原始数据
        if (pd.notna(row.get('fastq_r1')) or 
            pd.notna(row.get('run_id')) or 
            pd.notna(row.get('experiment_id'))):
            return 'Yes'
        return 'No'
    
    def determine_matrix_exist(self, row):
        """
        判断矩阵数据是否存在
        注意：CNGB通常不直接提供处理后的矩阵，所以这里大多数情况是'No'
        """
        # 只有在明确有表达矩阵文件时才返回Yes
        if pd.notna(row.get('expression_matrix')):
            return 'Yes'
        
        # SingleCell ID不一定表示有矩阵，可能只是原始数据
        # 所以不能仅凭singlecell_id判断
        return 'No'
    
    def infer_open_status(self, row):
        """推断开放状态"""
        # CNGB数据通常是开放的，如果有文件信息就认为是开放的
        if (pd.notna(row.get('fastq_r1')) or 
            pd.notna(row.get('run_id')) or 
            pd.notna(row.get('experiment_id')) or
            pd.notna(row.get('expression_matrix'))):
            return 'Open'
        return 'Unknown'
    
    def build_access_link(self, row):
        """构建访问链接"""
        project_id = row.get('project_id')
        sample_id = row.get('sample_id')
        experiment_id = row.get('experiment_id')
        singlecell_id = row.get('singlecell_id')
        run_id = row.get('run_id')
        
        links = []
        
        # 项目链接
        if pd.notna(project_id):
            links.append(f"https://db.cngb.org/search/project/{project_id}/")
        
        # 样本链接
        if pd.notna(sample_id):
            links.append(f"https://db.cngb.org/search/sample/{sample_id}/")
        
        # 实验链接
        if pd.notna(experiment_id):
            links.append(f"https://db.cngb.org/search/experiment/{experiment_id}/")
        
        # Run链接
        if pd.notna(run_id):
            links.append(f"https://db.cngb.org/search/run/{run_id}/")
        
        # 单细胞链接
        if pd.notna(singlecell_id):
            links.append(f"https://db.cngb.org/search/single-cell/{singlecell_id}/")
        
        return ';'.join(links) if links else None
    
    def infer_sample_type(self, row):
        """推断样本类型"""
        # 优先使用已有的sample_type
        if pd.notna(row.get('sample_type')):
            return row.get('sample_type')
        
        # 从组织、细胞系等信息推断
        tissue = str(row.get('tissue', '')).lower()
        cell_line = str(row.get('cell_line', '')).lower()
        
        if 'cell line' in cell_line or 'cell-line' in cell_line or pd.notna(row.get('cell_line')):
            return 'Cell Line'
        elif 'primary' in tissue:
            return 'Primary Tissue'
        elif 'blood' in tissue or 'pbmc' in tissue:
            return 'Blood/PBMC'
        elif pd.notna(row.get('tissue')):
            return 'Tissue'
        
        return None
    
    def standardize(self):
        """转换为标准格式"""
        logger.info("开始标准化处理...")
        
        # 加载数据
        super_df = self.load_super_table()
        
        # 创建标准化数据框
        std_data = []
        
        for idx, row in super_df.iterrows():
            if (idx + 1) % 1000 == 0:
                logger.info(f"处理进度: {idx + 1}/{len(super_df)}")
            
            std_row = {
                # ============ Study/Project IDs ============
                # CNGB的项目可能有多个ID，按优先级排列
                # ID槽1: 主项目ID（Project ID）
                'Study/Project_id_1': row.get('project_id'),
                
                # ID槽2: 实验ID（Experiment ID，属于项目下的子层级）
                'Study/Project_id_2': row.get('experiment_id'),
                
                # ID槽3: 单细胞ID或其他特殊ID（如果是单细胞项目）
                # 如果同一个项目还有其他数据库的ID（比如提交到GEO/SRA的ID），
                # 也可以填在这里，但CNGB本身通常只有这几种ID
                'Study/Project_id_3': row.get('singlecell_id'),
                
                # ============ Sample IDs ============
                # raw_sample_id: 原始测序数据的样本ID（CNGB的Sample ID）
                'raw_sample_id': row.get('sample_id'),
                
                # matrix_sample_id: 处理后矩阵的样本ID
                # CNGB通常不直接提供矩阵，所以这里留空
                # 如果研究者后续提交到其他数据库（如GEO），那个ID应该在另一个脚本中处理
                'matrix_sample_id': None,  # CNGB不提供矩阵
                
                # ============ Data Availability ============
                # 原始数据（FastQ/Run）是否存在
                'raw_exist': self.determine_raw_exist(row),
                'raw_open': self.infer_open_status(row),
                
                # 矩阵数据是否存在（CNGB几乎都是No）
                'matrix_exist': self.determine_matrix_exist(row),
                'matrix_open': self.infer_open_status(row) if self.determine_matrix_exist(row) == 'Yes' else None,
                
                # ============ File Type ============
                # 列出该样本在CNGB中存储的所有数据类型（用分号分隔）
                'file_type': self.determine_file_type(row),
                
                # ============ Study Information ============
                'title': row.get('experiment_title') or row.get('sample_title') or row.get('sample_name'),
                
                # ============ Disease Information ============
                'disease_general': self.infer_disease_general(row.get('disease')),
                'disease': row.get('disease'),
                
                # ============ Publication ============
                # CNGB数据通常没有直接的PubMed ID
                'pubmed': None,
                
                # ============ Source ============
                'source_database': 'CNGB',
                'access_link': self.build_access_link(row),
                'open_status': self.infer_open_status(row),
                
                # ============ Sample Metadata ============
                'ethnicity': row.get('ethnicity') or row.get('population'),
                'sex': row.get('sex'),
                'tissue_location': row.get('tissue'),
                
                # ============ Technical Information ============
                'sequencing_platform': row.get('platform') or row.get('instrument_model'),
                'experiment_design': row.get('library_strategy') or row.get('design_description'),
                'sample_type': self.infer_sample_type(row),
                
                # ============ Description ============
                'summary': row.get('design_description') or row.get('experiment_title'),
                
                # ============ Citation (CNGB通常不提供) ============
                'citation_count': None,
                
                # ============ Dates (CNGB通常不公开这些信息) ============
                'publication_date': None,
                'submission_date': None,
                'last_update_date': None,
                
                # ============ Contact Information (通常不公开) ============
                'contact_name': row.get('biomaterial_provider'),
                'contact_email': None,
                'contact_institute': row.get('biomaterial_provider'),
                
                # ============ Data Quality ============
                'data_tier': self.determine_data_tier(row),
                
                # ============ Supplementary Information ============
                'supplementary_information': self.build_supplementary_info(row),
            }
            
            std_data.append(std_row)
        
        # 创建标准化DataFrame
        std_df = pd.DataFrame(std_data)
        
        # 保存
        logger.info(f"保存标准化数据到: {self.output_path}")
        std_df.to_csv(self.output_path, index=False, encoding='utf-8-sig')
        
        logger.info(f"✓ 标准化完成！共 {len(std_df)} 条记录")
        
        # 生成统计报告
        self.generate_statistics(std_df)
        
        return std_df
    
    def determine_data_tier(self, row):
        """判断数据等级"""
        score = 0
        
        # 有原始数据 +1
        if self.determine_raw_exist(row) == 'Yes':
            score += 1
        
        # 有矩阵数据 +1（CNGB很少有）
        if self.determine_matrix_exist(row) == 'Yes':
            score += 1
        
        # 有疾病信息 +1
        if pd.notna(row.get('disease')):
            score += 1
        
        # 有完整的样本信息 +1
        if pd.notna(row.get('tissue')) and pd.notna(row.get('sex')):
            score += 1
        
        # 有测序平台信息 +0.5
        if pd.notna(row.get('platform')) or pd.notna(row.get('instrument_model')):
            score += 0.5
        
        # 分级
        if score >= 4:
            return 'Tier 1'
        elif score >= 3:
            return 'Tier 2'
        elif score >= 2:
            return 'Tier 3'
        else:
            return 'Tier 4'
    
    def build_supplementary_info(self, row):
        """构建补充信息"""
        info_parts = []
        
        # 添加关键的额外信息
        fields_to_include = [
            ('organism', 'Organism'),
            ('strain', 'Strain'),
            ('dev_stage', 'Dev Stage'),
            ('age', 'Age'),
            ('cell_type', 'Cell Type'),
            ('cell_line', 'Cell Line'),
            ('treatment', 'Treatment'),
            ('disease_stage', 'Disease Stage'),
            ('library_layout', 'Library Layout'),
            ('library_selection', 'Library Selection'),
            ('library_source', 'Library Source'),
            ('karyotype', 'Karyotype'),
            ('biomaterial_provider', 'Provider'),
        ]
        
        for field, label in fields_to_include:
            value = row.get(field)
            if pd.notna(value) and str(value).strip():
                info_parts.append(f"{label}: {value}")
        
        # 添加Run ID（如果有的话）
        if pd.notna(row.get('run_id')):
            info_parts.append(f"Run ID: {row.get('run_id')}")
        
        # 添加文件MD5信息（用于数据验证）
        if pd.notna(row.get('fastq_r1_md5')):
            info_parts.append(f"R1_MD5: {row.get('fastq_r1_md5')[:16]}...")  # 只显示前16位
        if pd.notna(row.get('fastq_r2_md5')):
            info_parts.append(f"R2_MD5: {row.get('fastq_r2_md5')[:16]}...")
        
        return '; '.join(info_parts) if info_parts else None
    
    def generate_statistics(self, std_df):
        """生成统计报告"""
        logger.info("\n" + "="*80)
        logger.info("CNGB 标准化数据统计")
        logger.info("="*80)
        
        # 基本统计
        logger.info(f"\n📊 基本信息:")
        logger.info(f"  总记录数: {len(std_df):,}")
        logger.info(f"  唯一项目数 (Project ID): {std_df['Study/Project_id_1'].nunique():,}")
        logger.info(f"  唯一实验数 (Experiment ID): {std_df['Study/Project_id_2'].nunique():,}")
        logger.info(f"  唯一样本数: {std_df['raw_sample_id'].nunique():,}")
        
        # 数据可用性统计
        logger.info(f"\n💾 数据可用性:")
        raw_yes = (std_df['raw_exist'] == 'Yes').sum()
        matrix_yes = (std_df['matrix_exist'] == 'Yes').sum()
        logger.info(f"  有原始数据 (FastQ/Run): {raw_yes:,} ({raw_yes/len(std_df)*100:.1f}%)")
        logger.info(f"  有矩阵数据: {matrix_yes:,} ({matrix_yes/len(std_df)*100:.1f}%)")
        logger.info(f"  开放数据: {(std_df['open_status'] == 'Open').sum():,} ({(std_df['open_status'] == 'Open').sum()/len(std_df)*100:.1f}%)")
        
        # 文件类型统计
        logger.info(f"\n📁 文件类型分布:")
        file_types_all = []
        for ft in std_df['file_type'].dropna():
            file_types_all.extend(ft.split(';'))
        from collections import Counter
        ft_counts = Counter(file_types_all)
        for ft, count in ft_counts.most_common():
            logger.info(f"  {ft}: {count:,} ({count/len(std_df)*100:.1f}%)")
        
        # 疾病分类统计
        if std_df['disease_general'].notna().any():
            logger.info(f"\n🏥 疾病大类分布:")
            disease_counts = std_df['disease_general'].value_counts()
            for disease, count in disease_counts.items():
                logger.info(f"  {disease}: {count:,} ({count/len(std_df)*100:.1f}%)")
        
        # 数据等级统计
        logger.info(f"\n⭐ 数据等级分布:")
        tier_counts = std_df['data_tier'].value_counts().sort_index()
        for tier, count in tier_counts.items():
            logger.info(f"  {tier}: {count:,} ({count/len(std_df)*100:.1f}%)")
        
        # 组织统计（Top 10）
        if std_df['tissue_location'].notna().any():
            logger.info(f"\n🧬 组织分布 (Top 10):")
            tissue_counts = std_df['tissue_location'].value_counts().head(10)
            for tissue, count in tissue_counts.items():
                logger.info(f"  {tissue}: {count:,}")
        
        # 测序平台统计
        if std_df['sequencing_platform'].notna().any():
            logger.info(f"\n🔬 测序平台分布:")
            platform_counts = std_df['sequencing_platform'].value_counts()
            for platform, count in platform_counts.items():
                logger.info(f"  {platform}: {count:,}")
        
        # 物种统计
        logger.info(f"\n🧫 物种分布:")
        organism_info = std_df['supplementary_information'].str.extract(r'Organism: ([^;]+)', expand=False)
        if organism_info.notna().any():
            organism_counts = organism_info.value_counts().head(5)
            for organism, count in organism_counts.items():
                logger.info(f"  {organism}: {count:,}")
        
        logger.info("\n" + "="*80)


def main():
    # 配置路径
    super_table_path = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/cngb/processed/comprehensive_table/CNGB_super_table.tsv"
    output_path = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/cngb/processed/cngb_std.csv"
    
    # 创建标准化器
    standardizer = CNGBStandardizer(super_table_path, output_path)
    
    # 执行标准化
    std_df = standardizer.standardize()
    
    logger.info("\n" + "="*80)
    logger.info("✅ CNGB数据标准化完成！")
    logger.info(f"📄 输出文件: {output_path}")
    logger.info("="*80)


if __name__ == "__main__":
    main()