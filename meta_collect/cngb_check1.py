import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter
import json

def analyze_metadata_structure(metadata_dir):
    """分析所有metadata文件的结构"""
    metadata_dir = Path(metadata_dir)
    
    # 统计数据
    stats = {
        'projects': [],
        'column_frequency': Counter(),  # 所有列名出现频率
        'file_types': defaultdict(list),  # 按文件类型分组的列名
        'column_samples': defaultdict(list),  # 每个列名的示例值
    }
    
    print("开始分析metadata文件结构...\n")
    print("="*80)
    
    # 遍历所有项目文件夹
    project_folders = sorted([d for d in metadata_dir.iterdir() if d.is_dir() and d.name.startswith('CNP')])
    
    for project_folder in project_folders[:5]:  # 先分析前5个项目
        project_id = project_folder.name
        print(f"\n项目: {project_id}")
        print("-"*80)
        
        project_info = {
            'project_id': project_id,
            'files': []
        }
        
        # 遍历项目中的所有文件
        for file in sorted(project_folder.iterdir()):
            if file.suffix.lower() not in ['.tsv', '.csv', '.txt', '.xlsx', '.json']:
                continue
            
            try:
                # 解析文件
                df = parse_file(file)
                if df is None or df.empty:
                    continue
                
                file_info = {
                    'filename': file.name,
                    'rows': len(df),
                    'columns': df.columns.tolist(),
                    'column_count': len(df.columns)
                }
                
                # 识别文件类型
                file_type = identify_file_type(df, file.name)
                file_info['type'] = file_type
                
                print(f"\n  文件: {file.name}")
                print(f"  类型: {file_type}")
                print(f"  行数: {len(df)}")
                print(f"  列数: {len(df.columns)}")
                print(f"  列名:")
                for col in df.columns:
                    # 获取示例值(非空的前3个)
                    sample_values = df[col].dropna().head(3).tolist()
                    sample_str = ', '.join([str(v)[:50] for v in sample_values])
                    print(f"    - {col}: {sample_str}")
                    
                    # 统计
                    stats['column_frequency'][col] += 1
                    stats['file_types'][file_type].append(col)
                    if len(stats['column_samples'][col]) < 5:
                        stats['column_samples'][col].extend(sample_values)
                
                project_info['files'].append(file_info)
                
            except Exception as e:
                print(f"  ✗ 解析失败: {file.name} - {str(e)}")
                continue
        
        stats['projects'].append(project_info)
    
    # 生成汇总报告
    print("\n\n" + "="*80)
    print("汇总统计")
    print("="*80)
    
    print(f"\n分析的项目数: {len(stats['projects'])}")
    print(f"发现的唯一列名数: {len(stats['column_frequency'])}")
    
    # 按文件类型统计列名
    print("\n按文件类型统计:")
    for file_type, columns in stats['file_types'].items():
        unique_cols = set(columns)
        print(f"\n  {file_type.upper()} 类型:")
        print(f"    唯一列名数: {len(unique_cols)}")
        print(f"    列名: {', '.join(sorted(unique_cols)[:20])}")
        if len(unique_cols) > 20:
            print(f"    ... 还有 {len(unique_cols)-20} 个")
    
    # 最常见的列名
    print("\n最常见的列名 (Top 50):")
    for col, count in stats['column_frequency'].most_common(50):
        sample = stats['column_samples'][col][:2]
        sample_str = ', '.join([str(v)[:30] for v in sample])
        print(f"  {col:40s} (出现{count}次) - 示例: {sample_str}")
    
    # 保存详细结果到JSON
    output_file = Path(metadata_dir).parent / "metadata_structure_analysis.json"
    
    # 准备可序列化的数据
    serializable_stats = {
        'projects': stats['projects'],
        'column_frequency': dict(stats['column_frequency']),
        'file_types': {k: list(set(v)) for k, v in stats['file_types'].items()},
        'column_samples': {k: v[:5] for k, v in stats['column_samples'].items()}
    }
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(serializable_stats, f, indent=2, ensure_ascii=False)
    
    print(f"\n详细分析结果已保存到: {output_file}")
    
    # 生成字段映射建议
    print("\n" + "="*80)
    print("字段映射建议")
    print("="*80)
    generate_mapping_suggestions(stats)
    
    return stats


def parse_file(file_path):
    """解析文件"""
    try:
        ext = file_path.suffix.lower()
        
        if ext in ['.tsv', '.txt']:
            for sep in ['\t', ',', '|']:
                try:
                    df = pd.read_csv(file_path, sep=sep, nrows=5)
                    if len(df.columns) > 1:
                        return pd.read_csv(file_path, sep=sep)
                except:
                    continue
        
        elif ext == '.csv':
            return pd.read_csv(file_path)
        
        elif ext == '.xlsx':
            return pd.read_excel(file_path)
        
        elif ext == '.json':
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            if isinstance(data, list):
                return pd.DataFrame(data)
            elif isinstance(data, dict):
                return pd.DataFrame([data])
        
    except Exception as e:
        return None


def identify_file_type(df, filename):
    """识别文件类型"""
    columns_lower = [col.lower() for col in df.columns]
    filename_lower = filename.lower()
    
    # 通过文件名判断
    if 'experiment' in filename_lower:
        return 'experiment'
    elif 'sample' in filename_lower:
        return 'sample'
    elif 'singlecell' in filename_lower or 'single_cell' in filename_lower:
        return 'singlecell'
    elif 'project' in filename_lower:
        return 'project'
    
    # 通过列名判断
    if any(keyword in columns_lower for keyword in ['experiment_accession', 'run_accession', 'library_strategy']):
        return 'experiment'
    elif any(keyword in columns_lower for keyword in ['sample_accession', 'sample_name', 'biosample']):
        return 'sample'
    elif any(keyword in columns_lower for keyword in ['expression_file', 'matrix', 'barcodes']):
        return 'singlecell'
    elif any(keyword in columns_lower for keyword in ['project_accession', 'project_title']):
        return 'project'
    
    return 'unknown'


def generate_mapping_suggestions(stats):
    """生成字段映射建议"""
    
    # 目标字段
    target_fields = {
        'sample_id': ['sample_accession', 'sample_id', 'sample_name', 'biosample_accession'],
        'sample_name': ['sample_name', 'sample_title', 'sample_alias'],
        'organism': ['organism', 'tax_id', 'taxid', 'species'],
        'tissue': ['tissue', 'tissue_type', 'organ'],
        'cell_type': ['cell_type', 'celltype', 'cell_line'],
        'disease': ['disease', 'disease_state', 'phenotype'],
        'disease_stage': ['disease_stage', 'stage', 'tumor_stage'],
        'age': ['age', 'age_at_diagnosis'],
        'sex': ['sex', 'gender'],
        'experiment_id': ['experiment_accession', 'run_accession', 'experiment_id'],
        'library_strategy': ['library_strategy', 'library_layout'],
        'library_source': ['library_source'],
        'platform': ['platform', 'instrument_platform'],
        'instrument_model': ['instrument_model', 'instrument', 'sequencer'],
        'fastq_r1': ['file_name', 'file1_name', 'fastq_file', 'read1'],
        'fastq_r2': ['file2_name', 'read2'],
        'fastq_r1_md5': ['file_md5', 'file1_md5', 'md5'],
        'fastq_r2_md5': ['file2_md5'],
        'expression_matrix': ['expression_file', 'matrix_file', 'expression_matrix'],
        'features_file': ['features_file', 'gene_file'],
        'barcodes_file': ['barcodes_file', 'barcode_file'],
    }
    
    print("\n建议的字段映射:")
    
    all_columns = set(stats['column_frequency'].keys())
    
    for target, source_patterns in target_fields.items():
        print(f"\n  {target}:")
        matched = []
        for col in all_columns:
            col_lower = col.lower()
            for pattern in source_patterns:
                if pattern.lower() in col_lower or col_lower in pattern.lower():
                    matched.append(col)
                    break
        
        if matched:
            for col in sorted(set(matched)):
                count = stats['column_frequency'][col]
                print(f"    ← {col} (出现{count}次)")
        else:
            print(f"    (未找到匹配)")


def main():
    metadata_dir = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/cngb/processed/meta_files"
    
    stats = analyze_metadata_structure(metadata_dir)
    
    print("\n分析完成!")


if __name__ == "__main__":
    main()