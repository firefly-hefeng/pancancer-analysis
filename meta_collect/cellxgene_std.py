import pandas as pd
import re

# 读取原始数据
input_file = '/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/processed/cellxgene_human_scrna_samples.csv'
output_file = '/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/processed/cellxgene_human_scrna_samples_std.csv'

df = pd.read_csv(input_file)

# 创建标准格式的DataFrame
std_df = pd.DataFrame()

# Study/Project ID 映射
std_df['Study/Project_id_1'] = df['collection_id']
std_df['Study/Project_id_2'] = df['dataset_id']
std_df['Study/Project_id_3'] = ''  # CellxGene没有第三层级ID

# Sample ID 映射
std_df['raw_sample_id'] = df['sample_id']
std_df['matrix_sample_id'] = df['sample_id']

# 数据存在性和开放性 (CellxGene的数据默认是开放的)
std_df['raw_exist'] = 'Yes'
std_df['raw_open'] = 'Yes'
std_df['matrix_exist'] = df['is_primary_data'].apply(lambda x: 'Yes' if x == True or str(x).lower() == 'true' else 'No')
std_df['matrix_open'] = 'Yes'

# 文件类型
std_df['file_type'] = 'scRNA-seq'

# 标题
std_df['title'] = df['collection_name']

# 疾病信息
def extract_disease_general(disease_str):
    """从disease字段提取一般疾病类型"""
    if pd.isna(disease_str) or disease_str == '':
        return ''
    # 提取括号前的疾病名称
    diseases = [d.split('(')[0].strip() for d in str(disease_str).split(';')]
    return '; '.join(diseases)

std_df['disease_general'] = df['diseases'].apply(extract_disease_general)
std_df['disease'] = df['disease']

# PubMed信息 (从DOI提取)
def extract_pubmed_from_doi(doi):
    """从DOI提取可能的PubMed信息"""
    if pd.isna(doi) or doi == '':
        return ''
    # DOI通常格式为 10.xxxx/xxxxx，可以通过DOI查询PubMed
    # 这里暂时返回空，需要实际API查询
    return ''

std_df['pubmed'] = df['doi'].apply(extract_pubmed_from_doi)

# 数据库来源
std_df['source_database'] = 'CellxGene'

# 访问链接
std_df['access_link'] = df['collection_id'].apply(
    lambda x: f'https://cellxgene.cziscience.com/collections/{x}' if pd.notna(x) else ''
)

# 开放状态
std_df['open_status'] = 'Open'

# 人群学信息
std_df['ethnicity'] = df['self_reported_ethnicity']
std_df['sex'] = df['sex']

# 组织位置
std_df['tissue_location'] = df['tissue']

# 测序平台
def extract_platform(assay_str):
    """从assay字段提取测序平台"""
    if pd.isna(assay_str) or assay_str == '':
        return ''
    # 提取括号前的技术名称
    return str(assay_str).split('(')[0].strip()

std_df['sequencing_platform'] = df['assay'].apply(extract_platform)

# 实验设计
std_df['experiment_design'] = 'Single-cell RNA-seq'

# 样本类型
std_df['sample_type'] = df['suspension_type']

# 摘要
std_df['summary'] = df['description']

# 引用次数 (需要外部API查询,暂时留空)
std_df['citation_count'] = ''

# 日期信息
std_df['publication_date'] = df['published_year'].apply(lambda x: str(int(x)) if pd.notna(x) else '')
std_df['submission_date'] = ''  # CellxGene没有提供
std_df['last_update_date'] = ''  # CellxGene没有提供

# 联系信息
std_df['contact_name'] = df['contact_name']
std_df['contact_email'] = df['contact_email']
std_df['contact_institute'] = ''  # CellxGene没有直接提供

# 数据层级
def determine_data_tier(n_cells, is_primary):
    """根据细胞数和是否为原始数据判断数据层级"""
    if pd.isna(n_cells):
        return 'Tier 3'
    if n_cells > 0 and (is_primary == True or str(is_primary).lower() == 'true'):
        return 'Tier 1'
    return 'Tier 2'

std_df['data_tier'] = df.apply(
    lambda row: determine_data_tier(row['n_cells'], row['is_primary_data']), 
    axis=1
)

# 补充信息
def create_supplementary_info(row):
    """创建补充信息字段"""
    info = []
    if pd.notna(row['n_cells']) and row['n_cells'] > 0:
        info.append(f"n_cells: {int(row['n_cells'])}")
    if pd.notna(row['n_cell_types']) and row['n_cell_types'] > 0:
        info.append(f"n_cell_types: {int(row['n_cell_types'])}")
    if pd.notna(row['tissue_ontology_term_id']):
        info.append(f"tissue_ontology: {row['tissue_ontology_term_id']}")
    if pd.notna(row['disease_ontology_term_id']):
        info.append(f"disease_ontology: {row['disease_ontology_term_id']}")
    if pd.notna(row['doi']):
        info.append(f"DOI: {row['doi']}")
    return '; '.join(info)

std_df['supplementary_information'] = df.apply(create_supplementary_info, axis=1)

# 保存转换后的数据
std_df.to_csv(output_file, index=False)

print(f"转换完成！")
print(f"原始数据行数: {len(df)}")
print(f"转换后数据行数: {len(std_df)}")
print(f"输出文件: {output_file}")

# 显示前几行预览
print("\n前3行预览:")
print(std_df.head(3).to_string())