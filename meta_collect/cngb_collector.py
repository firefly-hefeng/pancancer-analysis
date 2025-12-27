import pandas as pd
import requests
from bs4 import BeautifulSoup
import os
import time
from pathlib import Path
import json
import re
from urllib.parse import urljoin
import logging
from openai import OpenAI
import numpy as np

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Kimi API配置
KIMI_API_KEYS = [
    "sk-UL5YodR7ZL4S9dytpfMWJgmPTXJjkeNSd7Ktq9bbEhElzDfX",
    "sk-WL1hTKhW3sYKuJjBSc1k8wxL0r8ZLcuM7YiYxFjGgVHqAXhU",
    "sk-RAkA28HIT5tiMEfKtXAgbZ9nZKweq5Bnw0WbSwwBdNX7nbi1"
]
current_key_index = 0

def convert_to_json_serializable(obj):
    """将pandas/numpy类型转换为JSON可序列化的Python原生类型"""
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.Series):
        return obj.to_list()
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict('records')
    elif isinstance(obj, dict):
        return {k: convert_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_json_serializable(item) for item in obj]
    elif pd.isna(obj):
        return None
    else:
        return obj

class CNGBMetadataCollector:
    def __init__(self, input_csv, output_dir="./cngb_metadata"):
        self.input_csv = input_csv
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.meta_files_dir = self.output_dir / "meta_files"
        self.meta_files_dir.mkdir(exist_ok=True)
        
        # 读取输入CSV
        self.df = pd.read_csv(input_csv)
        logger.info(f"读取了 {len(self.df)} 个项目")
        
    def get_ftp_links(self, project_id):
        """获取项目的FTP元数据文件链接"""
        base_url = f"https://ftp2.cngb.org/pub/CNSA/data7/public_info/{project_id}/"
        ftp_links = []
        
        try:
            response = requests.get(base_url, timeout=30)
            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # 查找所有链接
                for link in soup.find_all('a'):
                    href = link.get('href')
                    if href and not href.startswith('?') and href != '../':
                        full_url = urljoin(base_url, href)
                        # 只获取元数据文件（通常是.txt, .tsv, .csv, .xml, .json等）
                        if any(full_url.endswith(ext) for ext in ['.txt', '.tsv', '.csv', '.xml', '.json', '.xlsx']):
                            ftp_links.append(full_url)
                
                logger.info(f"项目 {project_id} 找到 {len(ftp_links)} 个元数据文件")
            else:
                logger.warning(f"项目 {project_id} 无法访问: HTTP {response.status_code}")
                
        except Exception as e:
            logger.error(f"获取项目 {project_id} 的FTP链接时出错: {str(e)}")
            
        return ftp_links
    
    def download_metadata_file(self, url, project_id):
        """下载单个元数据文件"""
        try:
            filename = url.split('/')[-1]
            save_path = self.meta_files_dir / project_id / filename
            save_path.parent.mkdir(exist_ok=True, parents=True)
            
            if save_path.exists():
                logger.info(f"文件已存在，跳过: {filename}")
                return save_path
            
            response = requests.get(url, timeout=60)
            if response.status_code == 200:
                with open(save_path, 'wb') as f:
                    f.write(response.content)
                logger.info(f"下载成功: {filename}")
                return save_path
            else:
                logger.warning(f"下载失败: {filename}, HTTP {response.status_code}")
                return None
                
        except Exception as e:
            logger.error(f"下载文件 {url} 时出错: {str(e)}")
            return None
    
    def download_all_metadata(self):
        """下载所有项目的元数据文件"""
        all_downloads = {}
        
        for idx, row in self.df.iterrows():
            project_id = row['Project ID']
            logger.info(f"处理项目 {idx+1}/{len(self.df)}: {project_id}")
            
            # 获取FTP链接
            ftp_links = self.get_ftp_links(project_id)
            downloaded_files = []
            
            # 下载每个文件
            for link in ftp_links:
                file_path = self.download_metadata_file(link, project_id)
                if file_path:
                    downloaded_files.append(str(file_path))
                time.sleep(1)  # 避免请求过快
            
            all_downloads[project_id] = {
                'ftp_links': ftp_links,
                'downloaded_files': downloaded_files
            }
            
            time.sleep(2)  # 项目间延迟
        
        # 保存下载记录
        download_record = self.output_dir / "download_record.json"
        with open(download_record, 'w', encoding='utf-8') as f:
            json.dump(all_downloads, f, indent=2, ensure_ascii=False)
        
        logger.info(f"下载记录已保存到: {download_record}")
        return all_downloads
    
    def parse_metadata_file(self, file_path):
        """解析单个元数据文件"""
        try:
            file_ext = Path(file_path).suffix.lower()
            
            if file_ext in ['.txt', '.tsv']:
                # 尝试多种分隔符
                for sep in ['\t', ',', '|']:
                    try:
                        df = pd.read_csv(file_path, sep=sep, nrows=5)
                        if len(df.columns) > 1:
                            return pd.read_csv(file_path, sep=sep)
                    except:
                        continue
                        
            elif file_ext == '.csv':
                return pd.read_csv(file_path)
                
            elif file_ext == '.xlsx':
                return pd.read_excel(file_path)
                
            elif file_ext == '.json':
                with open(file_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                if isinstance(data, list):
                    return pd.DataFrame(data)
                elif isinstance(data, dict):
                    return pd.DataFrame([data])
                    
            elif file_ext == '.xml':
                # 简单的XML解析
                import xml.etree.ElementTree as ET
                tree = ET.parse(file_path)
                root = tree.getroot()
                # 这里需要根据实际XML结构调整
                return None
                
        except Exception as e:
            logger.error(f"解析文件 {file_path} 时出错: {str(e)}")
            return None
    
    def call_kimi_for_mapping(self, columns, first_row, target_fields):
        """调用Kimi API进行字段映射"""
        global current_key_index
        
        # 转换first_row为可序列化的格式
        first_row_serializable = convert_to_json_serializable(first_row)
        
        prompt = f"""你是一个数据字段映射专家。我有一些元数据字段需要映射到标准字段。

原始字段列名：
{', '.join(columns)}

第一行数据示例：
{json.dumps(first_row_serializable, ensure_ascii=False, indent=2)}

目标标准字段：
{', '.join(target_fields)}

请分析这些原始字段，并返回一个JSON格式的映射关系，格式如下：
{{
    "映射": {{
        "原始字段名1": "目标字段名1",
        "原始字段名2": "目标字段名2",
        ...
    }},
    "无法映射": ["无法确定映射的原始字段1", "无法确定映射的原始字段2", ...]
}}

注意：
1. 只映射有明确对应关系的字段
2. 如果原始字段可以合并或转换为目标字段，请说明
3. 如果无法确定映射关系，放入"无法映射"列表
"""

        for attempt in range(len(KIMI_API_KEYS)):
            try:
                client = OpenAI(
                    api_key=KIMI_API_KEYS[current_key_index],
                    base_url="https://api.moonshot.cn/v1"
                )
                
                response = client.chat.completions.create(
                    model="moonshot-v1-8k",
                    messages=[
                        {"role": "system", "content": "你是一个专业的数据字段映射助手，擅长理解生物信息学和基因组学数据结构。"},
                        {"role": "user", "content": prompt}
                    ],
                    temperature=0.3,
                )
                
                result = response.choices[0].message.content
                
                # 提取JSON
                json_match = re.search(r'\{.*\}', result, re.DOTALL)
                if json_match:
                    mapping = json.loads(json_match.group())
                    return mapping
                else:
                    return {"映射": {}, "无法映射": columns}
                    
            except Exception as e:
                logger.warning(f"Kimi API调用失败 (key {current_key_index}): {str(e)}")
                current_key_index = (current_key_index + 1) % len(KIMI_API_KEYS)
                if attempt == len(KIMI_API_KEYS) - 1:
                    logger.error("所有API key都失败了")
                    return {"映射": {}, "无法映射": columns}
                time.sleep(2)
    
    def extract_structured_data(self, project_id, metadata_files):
        """提取并结构化元数据"""
        target_fields = [
            'Study/Project_id', 'sample_id', 'Study/Project_title', 
            'disease_general', 'disease', 'pubmed', 'source_database',
            'access_link', 'open_status', 'ethnicity', 'sex', 
            'tissue_location', 'sequencing_platform', 'experiment_design',
            'sample_type', 'summary', 'citation_count', 'publication_date',
            'submission_date', 'last_update_date', 'contact_name',
            'contact_email', 'contact_institute', 'supplementary_information'
        ]
        
        structured_data = {field: None for field in target_fields}
        structured_data['Study/Project_id'] = project_id
        
        all_parsed_data = []
        
        for file_path in metadata_files:
            df = self.parse_metadata_file(file_path)
            if df is not None and not df.empty:
                all_parsed_data.append({
                    'file': file_path,
                    'dataframe': df
                })
                
                # 使用Kimi进行字段映射
                columns = df.columns.tolist()
                first_row = df.iloc[0].to_dict() if len(df) > 0 else {}
                
                logger.info(f"调用Kimi API映射字段: {Path(file_path).name}")
                mapping_result = self.call_kimi_for_mapping(
                    columns, 
                    first_row, 
                    target_fields
                )
                
                # 应用映射
                if '映射' in mapping_result:
                    for original_field, target_field in mapping_result['映射'].items():
                        if target_field in target_fields and original_field in df.columns:
                            # 取第一个非空值并转换为Python原生类型
                            value = df[original_field].dropna().iloc[0] if not df[original_field].dropna().empty else None
                            if value is not None and not structured_data[target_field]:
                                structured_data[target_field] = convert_to_json_serializable(value)
                
                time.sleep(1)  # API调用间隔
        
        # 确保所有字段都是可序列化的
        structured_data = convert_to_json_serializable(structured_data)
        
        return structured_data, all_parsed_data
    
    def process_all_projects(self):
        """处理所有项目"""
        # 加载下载记录
        download_record_path = self.output_dir / "download_record.json"
        if not download_record_path.exists():
            logger.info("未找到下载记录，开始下载元数据...")
            download_record = self.download_all_metadata()
        else:
            logger.info("加载已有下载记录...")
            with open(download_record_path, 'r', encoding='utf-8') as f:
                download_record = json.load(f)
        
        # 提取结构化数据
        all_structured_data = []
        
        for project_id, record in download_record.items():
            logger.info(f"提取项目 {project_id} 的结构化数据...")
            
            try:
                structured_data, parsed_data = self.extract_structured_data(
                    project_id, 
                    record['downloaded_files']
                )
                
                all_structured_data.append(structured_data)
                
                # 保存中间结果 - 确保数据可序列化
                interim_output = self.output_dir / f"{project_id}_structured.json"
                with open(interim_output, 'w', encoding='utf-8') as f:
                    json.dump({
                        'structured': structured_data,
                        'raw_files': [str(p['file']) for p in parsed_data]
                    }, f, indent=2, ensure_ascii=False)
                    
                logger.info(f"项目 {project_id} 处理完成")
                
            except Exception as e:
                logger.error(f"处理项目 {project_id} 时出错: {str(e)}")
                # 添加一个空记录以保持连续性
                all_structured_data.append({
                    'Study/Project_id': project_id,
                    'error': str(e)
                })
        
        # 合并所有结构化数据
        final_df = pd.DataFrame(all_structured_data)
        
        # 保存最终结果
        output_csv = self.output_dir / "cngb_structured_metadata.csv"
        final_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
        logger.info(f"最终结果已保存到: {output_csv}")
        
        output_excel = self.output_dir / "cngb_structured_metadata.xlsx"
        final_df.to_excel(output_excel, index=False)
        logger.info(f"Excel版本已保存到: {output_excel}")
        
        return final_df


def main():
    # 配置路径
    input_csv = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/cngb/cngb_raw.csv"
    output_dir = "/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/cngb/processed"
    
    # 创建收集器实例
    collector = CNGBMetadataCollector(input_csv, output_dir)
    
    # 处理所有项目
    result_df = collector.process_all_projects()
    
    print(f"\n处理完成！共处理 {len(result_df)} 个项目")
    print(f"结果保存在: {output_dir}")
    print(f"\n数据预览：")
    print(result_df.head())
    print(f"\n字段填充情况：")
    print(result_df.count())


if __name__ == "__main__":
    main()