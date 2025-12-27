import pandas as pd
import numpy as np

# 读取数据
psy_series = pd.read_csv('cleaned_psy_series.csv')
psy_samples = pd.read_csv('cleaned_psy_samples.csv')
htan_series = pd.read_csv('cleaned_htan_series.csv')
htan_samples = pd.read_csv('cleaned_htan_samples.csv')
ega_series = pd.read_csv('cleaned_ega_series.csv')
ega_samples = pd.read_csv('cleaned_ega_samples.csv')

# 定义数据类型统一函数
def normalize_merge_keys(df):
    """
    统一合并键的数据类型为字符串，并处理空值
    """
    merge_keys = ['Study/Project_id_1', 'Study/Project_id_2', 'Study/Project_id_3']
    
    for key in merge_keys:
        if key in df.columns:
            # 将所有值转换为字符串，NaN 转为空字符串
            df[key] = df[key].fillna('').astype(str)
            # 将 'nan' 字符串转回空字符串
            df[key] = df[key].replace('nan', '')
    
    return df

# 定义合并函数
def merge_series_to_samples(samples_df, series_df):
    """
    将 series 表的信息合并到 samples 表中
    
    参数:
    samples_df: samples 数据框
    series_df: series 数据框
    
    返回:
    合并后的标准化数据框
    """
    
    # 统一数据类型
    samples_df = normalize_merge_keys(samples_df.copy())
    series_df = normalize_merge_keys(series_df.copy())
    
    # 创建合并键
    merge_keys = ['Study/Project_id_1', 'Study/Project_id_2', 'Study/Project_id_3']
    
    # 执行左连接，samples 表为主表
    merged_df = samples_df.merge(
        series_df,
        on=merge_keys,
        how='left',
        suffixes=('_sample', '_series')
    )
    
    # 定义字段优先级处理函数
    def merge_column(row, col_name):
        """优先使用 samples 表的值，如果为空则使用 series 表的值"""
        sample_val = row.get(f'{col_name}_sample')
        series_val = row.get(f'{col_name}_series')
        
        # 处理各种空值情况
        def is_empty(val):
            if pd.isna(val):
                return True
            if val == '' or val == 'nan':
                return True
            if str(val).lower() == 'unknown':
                return True
            return False
        
        if is_empty(sample_val):
            return series_val if not is_empty(series_val) else sample_val
        return sample_val
    
    # 需要合并的字段列表
    merge_fields = [
        'title', 'disease_general', 'disease', 'ethnicity', 'sex', 
        'sequencing_platform', 'experiment_design', 'sample_type'
    ]
    
    # 处理每个需要合并的字段
    for field in merge_fields:
        sample_col = f'{field}_sample'
        series_col = f'{field}_series'
        
        if sample_col in merged_df.columns and series_col in merged_df.columns:
            merged_df[field] = merged_df.apply(lambda row: merge_column(row, field), axis=1)
            # 删除临时列
            merged_df = merged_df.drop([sample_col, series_col], axis=1)
        elif sample_col in merged_df.columns:
            merged_df[field] = merged_df[sample_col]
            merged_df = merged_df.drop([sample_col], axis=1)
        elif series_col in merged_df.columns:
            merged_df[field] = merged_df[series_col]
            merged_df = merged_df.drop([series_col], axis=1)
    
    # 处理 tissue_location（samples表中是tissue_location，series表中是tissue）
    if 'tissue_location' in merged_df.columns and 'tissue' in merged_df.columns:
        merged_df['tissue_location'] = merged_df.apply(
            lambda row: row['tissue_location'] if pd.notna(row['tissue_location']) and row['tissue_location'] != '' and row['tissue_location'] != 'nan'
            else row['tissue'], axis=1
        )
        merged_df = merged_df.drop(['tissue'], axis=1)
    elif 'tissue' in merged_df.columns:
        merged_df['tissue_location'] = merged_df['tissue']
        merged_df = merged_df.drop(['tissue'], axis=1)
    
    # 选择并排序最终需要的列
    final_columns = [
        'Study/Project_id_1', 'Study/Project_id_2', 'Study/Project_id_3',
        'raw_sample_id', 'matrix_sample_id',
        'raw_exist', 'raw_open', 'matrix_exist', 'matrix_open',
        'file_type', 'title', 'disease_general', 'disease',
        'pubmed', 'source_database', 'access_link', 'open_status',
        'ethnicity', 'sex', 'tissue_location',
        'sequencing_platform', 'experiment_design', 'sample_type',
        'summary', 'citation_count', 'publication_date',
        'submission_date', 'last_update_date',
        'contact_name', 'contact_email', 'contact_institute',
        'data_tier', 'supplementary_information'
    ]
    
    # 只保留存在的列
    available_columns = [col for col in final_columns if col in merged_df.columns]
    result_df = merged_df[available_columns].copy()
    
    return result_df

# 处理三个数据集
print("处理 PSY 数据集...")
psy_std = merge_series_to_samples(psy_samples, psy_series)
print(f"PSY 标准化表格: {psy_std.shape[0]} 行, {psy_std.shape[1]} 列")

print("\n处理 HTAN 数据集...")
htan_std = merge_series_to_samples(htan_samples, htan_series)
print(f"HTAN 标准化表格: {htan_std.shape[0]} 行, {htan_std.shape[1]} 列")

print("\n处理 EGA 数据集...")
ega_std = merge_series_to_samples(ega_samples, ega_series)
print(f"EGA 标准化表格: {ega_std.shape[0]} 行, {ega_std.shape[1]} 列")

# 保存结果
psy_std.to_csv('psy_std.csv', index=False, encoding='utf-8-sig')
htan_std.to_csv('htan_std.csv', index=False, encoding='utf-8-sig')
ega_std.to_csv('ega_std.csv', index=False, encoding='utf-8-sig')

print("\n保存完成!")
print("生成文件:")
print("- psy_std.csv")
print("- htan_std.csv")
print("- ega_std.csv")

# 显示每个数据集的前几行和基本信息
print("\n" + "="*80)
print("PSY 数据集预览:")
print(psy_std.head(3))

print("\n" + "="*80)
print("HTAN 数据集预览:")
print(htan_std.head(3))

print("\n" + "="*80)
print("EGA 数据集预览:")
print(ega_std.head(3))