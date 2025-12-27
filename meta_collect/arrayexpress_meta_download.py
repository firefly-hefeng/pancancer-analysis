import json
import pandas as pd
import requests
import os
from pathlib import Path
from tqdm import tqdm
import logging

# 配置
BASE_DIR = Path("/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/meta1/arrayexpress")
BIOSTUDIES_JSON = BASE_DIR / "biostudies_human_scrna_raw.json"
SCEA_JSON = BASE_DIR / "scea_human_scrna_raw.json"
OUTPUT_DIR = BASE_DIR / "processed_data"
METADATA_DIR = BASE_DIR / "metadata_files"

OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
METADATA_DIR.mkdir(exist_ok=True, parents=True)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def process_biostudies_data(data):
    records = []
    for item in data:
        records.append({
            'Accession': item.get('accession', ''),
            'Title': item.get('title', ''),
            'Type': item.get('type', ''),
            'Authors': item.get('author', ''),
            'Release_Date': item.get('release_date', ''),
            'Files': item.get('files', 0),
            'Views': item.get('views', 0),
            'Source': 'BioStudies'
        })
    return pd.DataFrame(records)

def process_scea_data(data):
    records = []
    for item in data:
        records.append({
            'Accession': item.get('experimentAccession', ''),
            'Title': item.get('experimentDescription', ''),
            'Species': item.get('species', ''),
            'Technology': ', '.join(item.get('technologyType', [])),
            'Assays': item.get('numberOfAssays', 0),
            'Source': 'SCEA'
        })
    return pd.DataFrame(records)

def download_metadata_from_api(accession, output_dir):
    """
    通过BioStudies API获取并下载IDF和SDRF元数据文件
    """
    api_url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}"
    exp_dir = output_dir / accession
    exp_dir.mkdir(exist_ok=True, parents=True)
    
    downloaded = []
    
    try:
        # 获取实验信息
        response = requests.get(api_url, timeout=30)
        if response.status_code != 200:
            return False, [], f"API返回 {response.status_code}"
        
        data = response.json()
        
        # 查找MAGE-TAB Files部分
        section = data.get('section', {})
        subsections = section.get('subsections', [])
        
        for subsection_group in subsections:
            if isinstance(subsection_group, list):
                for subsection in subsection_group:
                    # 只处理MAGE-TAB Files
                    if subsection.get('type') == 'MAGE-TAB Files':
                        files = subsection.get('files', [[]])[0]
                        
                        for file_info in files:
                            file_path = file_info.get('path', '')
                            
                            # 只下载IDF和SDRF文件
                            if file_path.endswith('.idf.txt') or file_path.endswith('.sdrf.txt'):
                                # 构建FTP URL
                                # E-MTAB-6149 -> E-MTAB-/149/E-MTAB-6149
                                acc_parts = accession.split('-')
                                if len(acc_parts) >= 3:
                                    prefix = '-'.join(acc_parts[:-1])
                                    number = acc_parts[-1]
                                    # 取后3位作为子目录
                                    suffix = number[-3:] if len(number) >= 3 else number
                                    ftp_path = f"ftp://ftp.ebi.ac.uk/biostudies/fire/{prefix}-/{suffix}/{accession}/Files/{file_path}"
                                    
                                    # 下载文件（使用HTTPS替换FTP）
                                    local_path = exp_dir / file_path
                                    try:
                                        https_url = ftp_path.replace('ftp://', 'https://')
                                        file_resp = requests.get(https_url, timeout=60)
                                        if file_resp.status_code == 200:
                                            with open(local_path, 'wb') as f:
                                                f.write(file_resp.content)
                                            downloaded.append(str(local_path))
                                            logging.info(f"  ✓ {file_path}")
                                        else:
                                            logging.debug(f"  ✗ {file_path} (HTTP {file_resp.status_code})")
                                    except Exception as e:
                                        logging.debug(f"  ✗ {file_path}: {e}")
        
        if downloaded:
            return True, downloaded, None
        else:
            return False, [], "未找到MAGE-TAB文件"
            
    except Exception as e:
        return False, [], str(e)

def main():
    logging.info("="*60)
    logging.info("ArrayExpress 元数据处理")
    logging.info("="*60)
    
    # 1. 加载数据
    logging.info("\n正在加载JSON数据...")
    bio_data = load_json(BIOSTUDIES_JSON)
    scea_data = load_json(SCEA_JSON)
    
    # 2. 处理为DataFrame
    logging.info("正在处理数据...")
    bio_df = process_biostudies_data(bio_data)
    scea_df = process_scea_data(scea_data)
    
    # 3. 保存表格
    logging.info("正在保存表格...")
    bio_df.to_csv(OUTPUT_DIR / "biostudies_summary.csv", index=False)
    bio_df.to_excel(OUTPUT_DIR / "biostudies_summary.xlsx", index=False, engine='openpyxl')
    logging.info(f"  BioStudies: {len(bio_df)} 条记录")
    
    scea_df.to_csv(OUTPUT_DIR / "scea_summary.csv", index=False)
    scea_df.to_excel(OUTPUT_DIR / "scea_summary.xlsx", index=False, engine='openpyxl')
    logging.info(f"  SCEA: {len(scea_df)} 条记录")
    
    # 4. 合并
    merged = pd.merge(bio_df[['Accession', 'Title', 'Release_Date']], 
                      scea_df[['Accession', 'Species', 'Technology', 'Assays']], 
                      on='Accession', how='outer')
    merged.to_csv(OUTPUT_DIR / "merged_summary.csv", index=False)
    merged.to_excel(OUTPUT_DIR / "merged_summary.xlsx", index=False, engine='openpyxl')
    logging.info(f"  合并表格: {len(merged)} 条记录")
    
    # 5. 下载元数据
    all_acc = set(bio_df['Accession'].tolist() + scea_df['Accession'].tolist())
    all_acc = [a for a in all_acc if a]
    
    logging.info(f"\n开始下载 {len(all_acc)} 个实验的元数据文件...")
    logging.info("(只下载 IDF 和 SDRF 文件)\n")
    
    results = {
        'success': [],
        'failed': [],
        'notfound': []
    }
    
    for acc in tqdm(all_acc, desc="下载进度"):
        success, files, error = download_metadata_from_api(acc, METADATA_DIR)
        
        if success:
            results['success'].append({
                'Accession': acc,
                'Files_Downloaded': len(files),
                'File_Names': ', '.join([Path(f).name for f in files])
            })
        elif "未找到MAGE-TAB文件" in str(error):
            results['notfound'].append(acc)
        else:
            results['failed'].append({
                'Accession': acc,
                'Error': error
            })
    
    # 6. 保存下载报告
    logging.info("\n保存下载报告...")
    if results['success']:
        success_df = pd.DataFrame(results['success'])
        success_df.to_csv(OUTPUT_DIR / 'download_success.csv', index=False)
        success_df.to_excel(OUTPUT_DIR / 'download_success.xlsx', index=False, engine='openpyxl')
    
    if results['failed']:
        failed_df = pd.DataFrame(results['failed'])
        failed_df.to_csv(OUTPUT_DIR / 'download_failed.csv', index=False)
        failed_df.to_excel(OUTPUT_DIR / 'download_failed.xlsx', index=False, engine='openpyxl')
    
    if results['notfound']:
        notfound_df = pd.DataFrame({'Accession': results['notfound']})
        notfound_df.to_csv(OUTPUT_DIR / 'download_notfound.csv', index=False)
        notfound_df.to_excel(OUTPUT_DIR / 'download_notfound.xlsx', index=False, engine='openpyxl')
    
    # 7. 输出总结
    logging.info("\n" + "="*60)
    logging.info("下载完成总结")
    logging.info("="*60)
    logging.info(f"✓ 成功下载: {len(results['success'])} 个实验")
    logging.info(f"✗ 下载失败: {len(results['failed'])} 个实验")
    logging.info(f"? 未找到文件: {len(results['notfound'])} 个实验")
    logging.info(f"\n输出目录: {OUTPUT_DIR}")
    logging.info(f"元数据目录: {METADATA_DIR}")
    logging.info("="*60)

if __name__ == "__main__":
    main()