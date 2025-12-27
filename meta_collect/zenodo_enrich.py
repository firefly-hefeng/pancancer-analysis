# enrichment_and_validation.py
"""
元数据增强和验证脚本
用于进一步丰富收集到的元数据，并验证数据质量
"""

import pandas as pd
import requests
from typing import Optional, Dict
import time
import logging
from Bio import Entrez
import re

logger = logging.getLogger(__name__)


class MetadataEnricher:
    """元数据增强器 - 从外部数据库获取额外信息"""
    
    def __init__(self, email: str = "your.email@example.com"):
        """
        初始化
        
        Parameters:
        -----------
        email : str
            用于NCBI API的邮箱（必需）
        """
        Entrez.email = email
        self.email = email
    
    def get_pubmed_info(self, pmid: str) -> Optional[Dict]:
        """
        从PubMed获取文献信息
        
        Parameters:
        -----------
        pmid : str
            PubMed ID
            
        Returns:
        --------
        Dict : 文献信息
        """
        if not pmid or pd.isna(pmid):
            return None
        
        try:
            handle = Entrez.efetch(db="pubmed", id=str(pmid), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            if not records['PubmedArticle']:
                return None
            
            article = records['PubmedArticle'][0]['MedlineCitation']['Article']
            
            # 提取作者
            authors = []
            if 'AuthorList' in article:
                for author in article['AuthorList'][:10]:  # 最多10个作者
                    if 'LastName' in author and 'Initials' in author:
                        authors.append(f"{author['LastName']} {author['Initials']}")
            
            # 提取期刊和年份
            journal = article.get('Journal', {}).get('Title', '')
            pub_date = article.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
            year = pub_date.get('Year', '')
            
            # 提取摘要
            abstract = ''
            if 'Abstract' in article and 'AbstractText' in article['Abstract']:
                abstract_parts = article['Abstract']['AbstractText']
                if isinstance(abstract_parts, list):
                    abstract = ' '.join(str(part) for part in abstract_parts)
                else:
                    abstract = str(abstract_parts)
            
            return {
                'authors': '; '.join(authors),
                'journal': journal,
                'year': year,
                'abstract': abstract,
                'title': article.get('ArticleTitle', '')
            }
            
        except Exception as e:
            logger.error(f"获取PubMed信息失败 {pmid}: {str(e)}")
            return None
    
    def get_citation_count(self, doi: str) -> Optional[int]:
        """
        从Crossref获取引用数
        
        Parameters:
        -----------
        doi : str
            DOI
            
        Returns:
        --------
        int : 引用数
        """
        if not doi or pd.isna(doi):
            return None
        
        try:
            url = f"https://api.crossref.org/works/{doi}"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                return data.get('message', {}).get('is-referenced-by-count', 0)
            
        except Exception as e:
            logger.error(f"获取引用数失败 {doi}: {str(e)}")
        
        return None
    
    def enrich_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        增强整个数据框
        
        Parameters:
        -----------
        df : pd.DataFrame
            原始数据框
            
        Returns:
        --------
        pd.DataFrame : 增强后的数据框
        """
        df_enriched = df.copy()
        
        # 添加新列
        df_enriched['pubmed_authors'] = None
        df_enriched['pubmed_journal'] = None
        df_enriched['pubmed_year'] = None
        df_enriched['pubmed_abstract'] = None
        
        # 遍历每一行
        for idx, row in df_enriched.iterrows():
            # 从PubMed获取信息
            if pd.notna(row['pubmed']):
                logger.info(f"处理 {idx+1}/{len(df_enriched)}: PMID {row['pubmed']}")
                pubmed_info = self.get_pubmed_info(row['pubmed'])
                
                if pubmed_info:
                    df_enriched.at[idx, 'pubmed_authors'] = pubmed_info.get('authors')
                    df_enriched.at[idx, 'pubmed_journal'] = pubmed_info.get('journal')
                    df_enriched.at[idx, 'pubmed_year'] = pubmed_info.get('year')
                    df_enriched.at[idx, 'pubmed_abstract'] = pubmed_info.get('abstract')
                
                time.sleep(0.4)  # NCBI限制每秒3个请求
            
            # 从Crossref获取引用数
            if pd.notna(row['doi']):
                citation_count = self.get_citation_count(row['doi'])
                if citation_count is not None:
                    df_enriched.at[idx, 'citation_count'] = citation_count
                
                time.sleep(0.2)
        
        return df_enriched


class DataValidator:
    """数据验证器 - 检查数据质量"""
    
    @staticmethod
    def validate_dataframe(df: pd.DataFrame) -> Dict:
        """
        验证数据框质量
        
        Returns:
        --------
        Dict : 验证报告
        """
        report = {
            'total_records': len(df),
            'completeness': {},
            'issues': [],
            'warnings': []
        }
        
        # 检查关键字段完整性
        key_fields = ['id', 'title', 'source_database', 'access_link']
        for field in key_fields:
            if field in df.columns:
                missing = df[field].isna().sum()
                report['completeness'][field] = {
                    'total': len(df),
                    'filled': len(df) - missing,
                    'missing': missing,
                    'percentage': (len(df) - missing) / len(df) * 100
                }
                
                if missing > 0:
                    report['warnings'].append(f"{field} 字段有 {missing} 条缺失")
        
        # 检查重要字段完整性
        important_fields = ['tissue', 'disease', 'sequencing_platform', 'pubmed']
        for field in important_fields:
            if field in df.columns:
                missing = df[field].isna().sum()
                report['completeness'][field] = {
                    'total': len(df),
                    'filled': len(df) - missing,
                    'missing': missing,
                    'percentage': (len(df) - missing) / len(df) * 100
                }
        
        # 检查数据一致性
        if 'open_status' in df.columns:
            status_values = df['open_status'].unique()
            expected_values = ['Open', 'Restricted', 'Closed', None]
            unexpected = set(status_values) - set(expected_values)
            if unexpected:
                report['warnings'].append(f"open_status 有意外值: {unexpected}")
        
        # 检查日期格式
        date_fields = ['publication_date', 'submission_date', 'last_update_date']
        for field in date_fields:
            if field in df.columns:
                try:
                    pd.to_datetime(df[field], errors='coerce')
                except:
                    report['issues'].append(f"{field} 日期格式可能有问题")
        
        return report


# 使用示例
if __name__ == "__main__":
    # 读取之前收集的数据
    df = pd.read_csv('zenodo_scrna_filtered_latest.csv')
    
    # 增强元数据
    enricher = MetadataEnricher(email="your.email@example.com")
    df_enriched = enricher.enrich_dataframe(df)
    
    # 保存增强后的数据
    df_enriched.to_csv('zenodo_scrna_enriched.csv', index=False)
    df_enriched.to_excel('zenodo_scrna_enriched.xlsx', index=False)
    
    # 验证数据
    validator = DataValidator()
    validation_report = validator.validate_dataframe(df_enriched)
    
    # 打印验证报告
    print("\n验证报告:")
    print("=" * 80)
    print(f"总记录数: {validation_report['total_records']}")
    print("\n字段完整性:")
    for field, stats in validation_report['completeness'].items():
        print(f"  {field}: {stats['filled']}/{stats['total']} ({stats['percentage']:.1f}%)")
    
    if validation_report['warnings']:
        print("\n警告:")
        for warning in validation_report['warnings']:
            print(f"  - {warning}")
    
    if validation_report['issues']:
        print("\n问题:")
        for issue in validation_report['issues']:
            print(f"  - {issue}")