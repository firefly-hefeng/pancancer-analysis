import json
import pandas as pd
from typing import List, Dict, Any
import os

def load_json_file(file_path: str) -> List[Dict]:
    """加载JSON文件"""
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data

def flatten_ontology_terms(terms: List[Dict]) -> str:
    """将本体术语列表转换为字符串格式"""
    if not terms:
        return ""
    return "; ".join([f"{term.get('label', '')} ({term.get('ontology_term_id', '')})" 
                      for term in terms])

def process_collections(collections_data: List[Dict]) -> pd.DataFrame:
    """处理collections数据"""
    collections_rows = []
    
    for collection in collections_data:
        # 基本信息
        base_info = {
            'collection_id': collection.get('collection_id', ''),
            'collection_url': collection.get('collection_url', ''),
            'collection_version_id': collection.get('collection_version_id', ''),
            'name': collection.get('name', ''),
            'description': collection.get('description', ''),
            'doi': collection.get('doi', ''),
            'visibility': collection.get('visibility', ''),
            'created_at': collection.get('created_at', ''),
            'published_at': collection.get('published_at', ''),
            'revised_at': collection.get('revised_at', ''),
            'contact_name': collection.get('contact_name', ''),
            'contact_email': collection.get('contact_email', ''),
            'curator_name': collection.get('curator_name', ''),
            'consortia': "; ".join(collection.get('consortia', [])),
        }
        
        # 出版信息
        publisher_metadata = collection.get('publisher_metadata', {})
        if publisher_metadata:
            base_info['journal'] = publisher_metadata.get('journal', '')
            base_info['published_year'] = publisher_metadata.get('published_year', '')
            base_info['published_month'] = publisher_metadata.get('published_month', '')
            base_info['published_day'] = publisher_metadata.get('published_day', '')
            base_info['is_preprint'] = publisher_metadata.get('is_preprint', '')
            
            # 作者信息
            authors = publisher_metadata.get('authors', [])
            author_list = []
            for author in authors:
                author_name = f"{author.get('given', '')} {author.get('family', '')}".strip()
                if author_name:
                    author_list.append(author_name)
            base_info['authors'] = "; ".join(author_list)
            base_info['author_count'] = len(author_list)
        
        # 链接信息
        links = collection.get('links', [])
        raw_data_links = []
        other_links = []
        for link in links:
            link_type = link.get('link_type', '')
            link_url = link.get('link_url', '')
            link_name = link.get('link_name', '')
            
            if link_type == 'RAW_DATA':
                raw_data_links.append(f"{link_name}: {link_url}" if link_name else link_url)
            else:
                other_links.append(f"{link_name}: {link_url}" if link_name else link_url)
        
        base_info['raw_data_links'] = "; ".join(raw_data_links)
        base_info['other_links'] = "; ".join(other_links)
        
        # 数据集统计
        datasets = collection.get('datasets', [])
        base_info['dataset_count'] = len(datasets)
        
        # 收集所有数据集中的唯一值
        all_organisms = set()
        all_diseases = set()
        all_tissues = set()
        all_assays = set()
        all_suspension_types = set()
        
        for dataset in datasets:
            # Organism
            for org in dataset.get('organism', []):
                all_organisms.add(f"{org.get('label', '')} ({org.get('ontology_term_id', '')})")
            
            # Disease
            for disease in dataset.get('disease', []):
                all_diseases.add(f"{disease.get('label', '')} ({disease.get('ontology_term_id', '')})")
            
            # Tissue
            for tissue in dataset.get('tissue', []):
                tissue_info = f"{tissue.get('label', '')} ({tissue.get('ontology_term_id', '')})"
                tissue_type = tissue.get('tissue_type', '')
                if tissue_type:
                    tissue_info += f" [{tissue_type}]"
                all_tissues.add(tissue_info)
            
            # Assay
            for assay in dataset.get('assay', []):
                all_assays.add(f"{assay.get('label', '')} ({assay.get('ontology_term_id', '')})")
            
            # Suspension type
            all_suspension_types.update(dataset.get('suspension_type', []))
        
        base_info['organisms'] = "; ".join(sorted(all_organisms))
        base_info['diseases'] = "; ".join(sorted(all_diseases))
        base_info['tissues'] = "; ".join(sorted(all_tissues))
        base_info['assays'] = "; ".join(sorted(all_assays))
        base_info['suspension_types'] = "; ".join(sorted(all_suspension_types))
        
        collections_rows.append(base_info)
    
    return pd.DataFrame(collections_rows)

def process_datasets(collections_data: List[Dict]) -> pd.DataFrame:
    """处理datasets数据 - 展开每个collection中的datasets"""
    dataset_rows = []
    
    for collection in collections_data:
        collection_id = collection.get('collection_id', '')
        collection_name = collection.get('name', '')
        
        for dataset in collection.get('datasets', []):
            dataset_info = {
                'collection_id': collection_id,
                'collection_name': collection_name,
                'dataset_id': dataset.get('dataset_id', ''),
                'dataset_version_id': dataset.get('dataset_version_id', ''),
            }
            
            # Organism
            dataset_info['organisms'] = flatten_ontology_terms(dataset.get('organism', []))
            
            # Disease
            dataset_info['diseases'] = flatten_ontology_terms(dataset.get('disease', []))
            
            # Tissue
            tissues = dataset.get('tissue', [])
            tissue_strs = []
            for tissue in tissues:
                tissue_str = f"{tissue.get('label', '')} ({tissue.get('ontology_term_id', '')})"
                tissue_type = tissue.get('tissue_type', '')
                if tissue_type:
                    tissue_str += f" [{tissue_type}]"
                tissue_strs.append(tissue_str)
            dataset_info['tissues'] = "; ".join(tissue_strs)
            
            # Assay
            dataset_info['assays'] = flatten_ontology_terms(dataset.get('assay', []))
            
            # Suspension type
            dataset_info['suspension_types'] = "; ".join(dataset.get('suspension_type', []))
            
            dataset_rows.append(dataset_info)
    
    return pd.DataFrame(dataset_rows)

def main():
    # 文件路径
    collections_path = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/raw/collections_raw.json"
    datasets_path = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/raw/datasets_raw.json"
    
    # 输出路径
    output_dir = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/metadata_output/cellxgene/processed"
    os.makedirs(output_dir, exist_ok=True)
    
    # 加载数据
    print("加载JSON文件...")
    collections_data = load_json_file(collections_path)
    
    # 处理collections
    print("处理collections数据...")
    collections_df = process_collections(collections_data)
    
    # 处理datasets
    print("处理datasets数据...")
    datasets_df = process_datasets(collections_data)
    
    # 保存为CSV
    collections_output = os.path.join(output_dir, "collections_table.csv")
    datasets_output = os.path.join(output_dir, "datasets_table.csv")
    
    print(f"保存collections表格到: {collections_output}")
    collections_df.to_csv(collections_output, index=False, encoding='utf-8-sig')
    
    print(f"保存datasets表格到: {datasets_output}")
    datasets_df.to_csv(datasets_output, index=False, encoding='utf-8-sig')
    
    # 也保存为Excel格式（可选）
    collections_excel = os.path.join(output_dir, "collections_table.xlsx")
    datasets_excel = os.path.join(output_dir, "datasets_table.xlsx")
    
    print(f"保存collections表格到: {collections_excel}")
    collections_df.to_excel(collections_excel, index=False, engine='openpyxl')
    
    print(f"保存datasets表格到: {datasets_excel}")
    datasets_df.to_excel(datasets_excel, index=False, engine='openpyxl')
    
    # 打印统计信息
    print("\n" + "="*60)
    print("处理完成！统计信息：")
    print("="*60)
    print(f"Collections数量: {len(collections_df)}")
    print(f"Datasets数量: {len(datasets_df)}")
    print(f"\nCollections表格列数: {len(collections_df.columns)}")
    print(f"Datasets表格列数: {len(datasets_df.columns)}")
    
    print("\nCollections表格列名:")
    for col in collections_df.columns:
        print(f"  - {col}")
    
    print("\nDatasets表格列名:")
    for col in datasets_df.columns:
        print(f"  - {col}")
    
    # 显示前几行
    print("\n" + "="*60)
    print("Collections表格预览（前3行）:")
    print("="*60)
    print(collections_df.head(3).to_string())
    
    print("\n" + "="*60)
    print("Datasets表格预览（前3行）:")
    print("="*60)
    print(datasets_df.head(3).to_string())

if __name__ == "__main__":
    main()