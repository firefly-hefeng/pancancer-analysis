#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from ftplib import FTP, error_perm
import logging
from tqdm import tqdm
import socket
import time

# ==================== 配置区域 ====================
# 输入：merged_summary.csv文件路径
MERGED_CSV_PATH = Path("/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/meta1/arrayexpress/processed_data/merged_summary.csv")

# 输出：元文件存储根目录
OUTPUT_ROOT = Path("/mnt/public7/pancancercol/meta_analysis/hefeng_metadata/meta_data_collect/meta1/arrayexpress/metadata_files_batch_ftp")

# FTP服务器配置
FTP_HOST = "ftp.ebi.ac.uk"
FTP_USER = "anonymous"
FTP_PASS = ""

# 日志设置
LOG_LEVEL = logging.INFO
# =================================================

# 创建输出目录
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

# 配置日志
log_file = OUTPUT_ROOT / "download_batch.log"
logging.basicConfig(
    level=LOG_LEVEL,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file, encoding='utf-8'),
        logging.StreamHandler()
    ]
)

def build_ftp_path(accession: str) -> str:
    """根据Accession构建FTP路径"""
    acc_parts = accession.split('-')
    if len(acc_parts) < 3:
        raise ValueError(f"Invalid accession format: {accession}")
    
    prefix = '-'.join(acc_parts[:-1])  # "E-MTAB"
    number = acc_parts[-1]             # "6149"
    suffix = number[-3:] if len(number) >= 3 else number  # "149"
    
    return f"biostudies/fire/{prefix}-/{suffix}/{accession}/Files"

def download_metadata_ftp(accession: str, local_dir: Path, max_retries: int = 3) -> tuple:
    """
    通过FTP下载指定Accession的IDF和SDRF文件
    
    Returns:
        tuple: (success: bool, downloaded_files: list, error_msg: str)
    """
    local_dir.mkdir(parents=True, exist_ok=True)
    downloaded_files = []
    ftp = None
    
    for attempt in range(1, max_retries + 1):
        try:
            # 1. 建立FTP连接（被动模式）
            ftp = FTP()
            ftp.connect(host=FTP_HOST, port=21, timeout=60)
            ftp.login(user=FTP_USER, passwd=FTP_PASS)
            ftp.set_pasv(True)  # 启用被动模式
            
            # 2. 构建并进入FTP路径
            ftp_path = build_ftp_path(accession)
            try:
                ftp.cwd(ftp_path)
            except error_perm:
                return False, [], f"Directory not found: {ftp_path}"
            
            # 3. 筛选目标文件
            all_files = ftp.nlst()
            target_files = [f for f in all_files if f.endswith(('.idf.txt', '.sdrf.txt'))]
            
            if not target_files:
                return False, [], f"No MAGE-TAB files found in {ftp_path}"
            
            # 4. 下载文件
            for filename in target_files:
                local_path = local_dir / filename
                
                # 检查文件是否已存在（避免重复下载）
                if local_path.exists():
                    logging.debug(f"  → {filename} already exists, skipping")
                    downloaded_files.append(str(local_path))
                    continue
                
                # 使用二进制模式下载
                with open(local_path, 'wb') as f:
                    ftp.retrbinary(f'RETR {filename}', f.write)
                
                file_size = local_path.stat().st_size
                downloaded_files.append(str(local_path))
                logging.info(f"  ✓ {filename} ({file_size:,} bytes)")
            
            # 5. 成功下载并断开连接
            ftp.quit()
            return True, downloaded_files, None
            
        except socket.timeout:
            wait_time = attempt * 10
            logging.warning(f"  Timeout on attempt {attempt}/{max_retries}, waiting {wait_time}s...")
            time.sleep(wait_time)
            
        except error_perm as e:
            error_msg = f"FTP permission error: {e}"
            logging.error(f"  {error_msg}")
            return False, [], error_msg
            
        except Exception as e:
            logging.error(f"  Unexpected error on attempt {attempt}: {e}")
            
        finally:
            if ftp:
                try:
                    ftp.quit()
                except:
                    pass
    
    # 所有重试失败
    return False, [], f"Failed after {max_retries} attempts"

def main():
    logging.info("="*80)
    logging.info("ArrayExpress 批量FTP下载工具")
    logging.info(f"输入文件: {MERGED_CSV_PATH}")
    logging.info(f"输出目录: {OUTPUT_ROOT}")
    logging.info("="*80)
    
    # 1. 读取merged_summary.csv
    if not MERGED_CSV_PATH.exists():
        logging.error(f"❌ 输入文件不存在: {MERGED_CSV_PATH}")
        return
    
    try:
        df = pd.read_csv(MERGED_CSV_PATH)
        accessions = df['Accession'].dropna().unique().tolist()
        logging.info(f"📊 从merged_summary.csv中加载了 {len(accessions)} 个Accession")
    except Exception as e:
        logging.error(f"❌ 读取CSV文件失败: {e}")
        return
    
    # 2. 初始化统计
    results = {
        'success': [],
        'failed': [],
        'notfound': []
    }
    
    # 3. 批量下载
    logging.info("\n🚀 开始批量下载元数据文件...")
    logging.info("(仅下载 .idf.txt 和 .sdrf.txt)\n")
    
    for accession in tqdm(accessions, desc="总体进度", unit="exp"):
        local_dir = OUTPUT_ROOT / accession
        
        success, files, error = download_metadata_ftp(accession, local_dir)
        
        if success:
            results['success'].append({
                'Accession': accession,
                'Files_Downloaded': len(files),
                'File_List': ', '.join([Path(f).name for f in files])
            })
        elif "No MAGE-TAB files found" in error:
            results['notfound'].append(accession)
            logging.warning(f"⚠️  Accession {accession}: {error}")
        else:
            results['failed'].append({
                'Accession': accession,
                'Error': error
            })
            logging.error(f"✗ Accession {accession}: {error}")
    
    # 4. 生成报告
    logging.info("\n" + "="*80)
    logging.info("生成下载报告...")
    
    # 成功报告
    if results['success']:
        success_df = pd.DataFrame(results['success'])
        success_df.to_csv(OUTPUT_ROOT / 'download_success_report.csv', index=False)
        logging.info(f"  ✓ 成功报告: {len(success_df)} 条记录")
    
    # 失败报告
    if results['failed']:
        failed_df = pd.DataFrame(results['failed'])
        failed_df.to_csv(OUTPUT_ROOT / 'download_failed_report.csv', index=False)
        logging.info(f"  ✗ 失败报告: {len(failed_df)} 条记录")
    
    # 未找到报告
    if results['notfound']:
        notfound_df = pd.DataFrame({'Accession': results['notfound']})
        notfound_df.to_csv(OUTPUT_ROOT / 'download_notfound_report.csv', index=False)
        logging.info(f"  ? 未找到报告: {len(notfound_df)} 条记录")
    
    # 5. 总结
    total = len(accessions)
    success_count = len(results['success'])
    failed_count = len(results['failed'])
    notfound_count = len(results['notfound'])
    
    logging.info("\n" + "="*80)
    logging.info("🎯 下载任务完成总结")
    logging.info("="*80)
    logging.info(f"📦 总计处理: {total} 个实验")
    logging.info(f"✅ 下载成功: {success_count} ({success_count/total*100:.1f}%)")
    logging.info(f"❌ 下载失败: {failed_count} ({failed_count/total*100:.1f}%)")
    logging.info(f"⚠️  文件缺失: {notfound_count} ({notfound_count/total*100:.1f}%)")
    logging.info(f"\n📂 所有文件存储在: {OUTPUT_ROOT}")
    logging.info(f"📝 详细日志请查看: {log_file}")
    logging.info("="*80)

if __name__ == "__main__":
    main()