"""
Figshare scRNA-seq Metadata Collection System
用于全面检索和整理Figshare数据库中的人类单细胞测序数据
"""

import requests
import pandas as pd
import json
import time
from datetime import datetime
from typing import List, Dict, Optional
import re
from collections import defaultdict
import logging
from pathlib import Path

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('figshare_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class FigshareScRNASeqCollector:
    """Figshare单细胞测序数据收集器"""
    
    def __init__(self, output_dir: str = "./figshare_data"):
        """
        初始化收集器
        
        Args:
            output_dir: 输出目录路径
        """
        self.base_url = "https://api.figshare.com/v2"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 定义搜索关键词组合
        self.search_terms = [
            "single cell RNA-seq human",
            "scRNA-seq human",
            "single-cell transcriptome human",
            "single cell sequencing human",
            "10x genomics human",
            "drop-seq human",
            "smart-seq human",
            "single cell gene expression human",
            "sc-RNA-seq human",
            "human single cell atlas",
            "human cell atlas"
        ]
        
        # 疾病相关关键词
        self.disease_terms = [
            "cancer", "tumor", "carcinoma", "alzheimer", "parkinson",
            "diabetes", "covid", "inflammation", "autism", "schizophrenia"
        ]
        
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ScRNASeq-Meta-Analysis/1.0'
        })
        
    def search_articles(self, search_term: str, page: int = 1, 
                       page_size: int = 100) -> List[Dict]:
        """
        搜索Figshare文章
        
        Args:
            search_term: 搜索关键词
            page: 页码
            page_size: 每页数量
            
        Returns:
            文章列表
        """
        url = f"{self.base_url}/articles/search"
        params = {
            'search_for': search_term,
            'page': page,
            'page_size': page_size,
            'order': 'published_date',
            'order_direction': 'desc'
        }
        
        try:
            response = self.session.post(url, json=params)
            response.raise_for_status()
            time.sleep(0.5)  # 避免请求过快
            return response.json()
        except Exception as e:
            logger.error(f"搜索失败 {search_term} (页 {page}): {str(e)}")
            return []
    
    def get_article_details(self, article_id: int) -> Optional[Dict]:
        """
        获取文章详细信息
        
        Args:
            article_id: 文章ID
            
        Returns:
            文章详细信息
        """
        url = f"{self.base_url}/articles/{article_id}"
        
        try:
            response = self.session.get(url)
            response.raise_for_status()
            time.sleep(0.3)
            return response.json()
        except Exception as e:
            logger.error(f"获取文章详情失败 {article_id}: {str(e)}")
            return None
    
    def get_article_versions(self, article_id: int) -> List[Dict]:
        """获取文章所有版本"""
        url = f"{self.base_url}/articles/{article_id}/versions"
        
        try:
            response = self.session.get(url)
            response.raise_for_status()
            time.sleep(0.3)
            return response.json()
        except Exception as e:
            logger.error(f"获取版本信息失败 {article_id}: {str(e)}")
            return []
    
    def collect_all_articles(self) -> List[Dict]:
        """
        收集所有相关文章
        
        Returns:
            所有文章的原始数据列表
        """
        all_articles = {}
        
        for search_term in self.search_terms:
            logger.info(f"搜索关键词: {search_term}")
            page = 1
            
            while True:
                articles = self.search_articles(search_term, page)
                
                if not articles:
                    break
                
                logger.info(f"  第 {page} 页: 找到 {len(articles)} 条记录")
                
                for article in articles:
                    article_id = article.get('id')
                    if article_id and article_id not in all_articles:
                        # 获取详细信息
                        details = self.get_article_details(article_id)
                        if details:
                            all_articles[article_id] = details
                            logger.info(f"    收集文章 {article_id}: {details.get('title', 'N/A')[:50]}")
                
                # 如果返回的记录少于page_size,说明已经到最后一页
                if len(articles) < 100:
                    break
                    
                page += 1
                time.sleep(1)  # 避免请求过快
        
        logger.info(f"总共收集到 {len(all_articles)} 个唯一文章")
        return list(all_articles.values())
    
    def extract_metadata_fields(self, article: Dict) -> Dict:
        """
        从文章中提取所需的元数据字段
        
        Args:
            article: 文章详细信息
            
        Returns:
            提取的元数据字典
        """
        # 基础信息
        article_id = article.get('id', '')
        title = article.get('title', '')
        description = article.get('description', '')
        
        # 日期信息
        publication_date = article.get('published_date', '')
        submission_date = article.get('created_date', '')
        last_update_date = article.get('modified_date', '')
        
        # 作者信息
        authors = article.get('authors', [])
        contact_name = ''
        contact_email = ''
        contact_institute = ''
        
        if authors:
            first_author = authors[0]
            contact_name = first_author.get('full_name', '')
            contact_email = first_author.get('email', '')
            # Figshare可能没有直接的机构信息,尝试从其他字段提取
            
        # 引用信息
        citation = article.get('citation', '')
        citation_count = article.get('citation_count', 0)
        
        # DOI和链接
        doi = article.get('doi', '')
        access_link = article.get('url', '') or article.get('url_public_html', '')
        
        # 文件信息
        files = article.get('files', [])
        data_tier = self._determine_data_tier(files, description)
        
        # 标签和分类
        tags = article.get('tags', [])
        categories = article.get('categories', [])
        
        # 从描述和标题中提取信息
        disease_general, disease = self._extract_disease_info(title, description, tags)
        ethnicity = self._extract_ethnicity(description)
        sex = self._extract_sex(description)
        tissue = self._extract_tissue(title, description, tags)
        sequencing_platform = self._extract_platform(title, description)
        experiment_design = self._extract_experiment_design(description)
        sample_type = self._extract_sample_type(description)
        tissue_location = self._extract_tissue_location(description)
        
        # PubMed信息
        pubmed = self._extract_pubmed(description, citation)
        
        # 开放状态
        open_status = 'open' if article.get('is_public', True) else 'restricted'
        
        # 补充信息
        supplementary_info = {
            'tags': tags,
            'categories': [cat.get('title', '') for cat in categories],
            'file_count': len(files),
            'total_size': sum(f.get('size', 0) for f in files),
            'license': article.get('license', {}).get('name', ''),
            'resource_doi': article.get('resource_doi', ''),
            'resource_title': article.get('resource_title', ''),
            'funding': article.get('funding', ''),
            'references': article.get('references', [])
        }
        
        return {
            'id': article_id,
            'sample_id': f"FIGSHARE_{article_id}",
            'title': title,
            'disease_general': disease_general,
            'disease': disease,
            'pubmed': pubmed,
            'source_database': 'Figshare',
            'access_link': access_link,
            'open_status': open_status,
            'ethnicity': ethnicity,
            'sex': sex,
            'tissue': tissue,
            'sequencing_platform': sequencing_platform,
            'experiment_design': experiment_design,
            'sample_type': sample_type,
            'summary': description[:500],  # 截取前500字符
            'citation_count': citation_count,
            'publication_date': publication_date,
            'submission_date': submission_date,
            'last_update_date': last_update_date,
            'contact_name': contact_name,
            'contact_email': contact_email,
            'contact_institute': contact_institute,
            'data_tier': data_tier,
            'tissue_location': tissue_location,
            'supplementary_information': json.dumps(supplementary_info, ensure_ascii=False)
        }
    
    def _determine_data_tier(self, files: List[Dict], description: str) -> str:
        """判断数据层级(raw/processed/matrix)"""
        tiers = set()
        
        for file_info in files:
            filename = file_info.get('name', '').lower()
            
            if any(ext in filename for ext in ['.fastq', '.bam', '.sam']):
                tiers.add('raw')
            if any(ext in filename for ext in ['.h5', '.h5ad', '.rds', '.loom']):
                tiers.add('processed')
            if 'matrix' in filename or 'count' in filename:
                tiers.add('matrix')
        
        # 从描述中推断
        desc_lower = description.lower()
        if 'raw data' in desc_lower or 'fastq' in desc_lower:
            tiers.add('raw')
        if 'processed' in desc_lower or 'normalized' in desc_lower:
            tiers.add('processed')
        if 'count matrix' in desc_lower or 'expression matrix' in desc_lower:
            tiers.add('matrix')
        
        return '/'.join(sorted(tiers)) if tiers else 'unknown'
    
    def _extract_disease_info(self, title: str, description: str, tags: List[str]) -> tuple:
        """提取疾病信息"""
        text = (title + ' ' + description + ' ' + ' '.join(tags)).lower()
        
        disease_mapping = {
            'cancer': ['cancer', 'carcinoma', 'tumor', 'tumour', 'malignant'],
            'neurological': ['alzheimer', 'parkinson', 'dementia', 'neurodegeneration'],
            'metabolic': ['diabetes', 'obesity', 'metabolic syndrome'],
            'infectious': ['covid', 'virus', 'infection', 'bacterial'],
            'autoimmune': ['lupus', 'arthritis', 'autoimmune'],
            'psychiatric': ['depression', 'schizophrenia', 'autism', 'psychiatric'],
            'cardiovascular': ['heart', 'cardiac', 'cardiovascular'],
            'healthy': ['healthy', 'normal', 'control']
        }
        
        disease_general = []
        disease_specific = []
        
        for category, keywords in disease_mapping.items():
            for keyword in keywords:
                if keyword in text:
                    disease_general.append(category)
                    disease_specific.append(keyword)
                    break
        
        return (
            ','.join(set(disease_general)) if disease_general else 'unknown',
            ','.join(set(disease_specific)) if disease_specific else 'unknown'
        )
    
    def _extract_ethnicity(self, description: str) -> str:
        """提取种族信息"""
        text = description.lower()
        
        ethnicities = {
            'asian': ['asian', 'chinese', 'japanese', 'korean', 'indian'],
            'european': ['european', 'caucasian', 'white'],
            'african': ['african', 'black'],
            'hispanic': ['hispanic', 'latino'],
            'mixed': ['mixed', 'diverse']
        }
        
        found = []
        for eth, keywords in ethnicities.items():
            if any(kw in text for kw in keywords):
                found.append(eth)
        
        return ','.join(found) if found else 'not specified'
    
    def _extract_sex(self, description: str) -> str:
        """提取性别信息"""
        text = description.lower()
        
        if 'male and female' in text or 'both sexes' in text:
            return 'mixed'
        elif 'female' in text and 'male' not in text.replace('female', ''):
            return 'female'
        elif 'male' in text:
            return 'male'
        
        return 'not specified'
    
    def _extract_tissue(self, title: str, description: str, tags: List[str]) -> str:
        """提取组织信息"""
        text = (title + ' ' + description + ' ' + ' '.join(tags)).lower()
        
        tissues = [
            'brain', 'liver', 'kidney', 'heart', 'lung', 'pancreas',
            'blood', 'pbmc', 'bone marrow', 'spleen', 'thymus',
            'skin', 'muscle', 'adipose', 'intestine', 'colon',
            'breast', 'prostate', 'ovary', 'testis', 'uterus',
            'eye', 'retina', 'cortex', 'hippocampus', 'cerebellum'
        ]
        
        found = [tissue for tissue in tissues if tissue in text]
        return ','.join(found) if found else 'not specified'
    
    def _extract_platform(self, title: str, description: str) -> str:
        """提取测序平台"""
        text = (title + ' ' + description).lower()
        
        platforms = {
            '10x Genomics': ['10x', '10x genomics', 'chromium'],
            'Drop-seq': ['drop-seq', 'dropseq'],
            'Smart-seq': ['smart-seq', 'smartseq', 'smart-seq2'],
            'MARS-seq': ['mars-seq', 'marsseq'],
            'CEL-seq': ['cel-seq', 'celseq'],
            'Seq-Well': ['seq-well', 'seqwell'],
            'inDrop': ['indrop'],
            'STRT-seq': ['strt-seq']
        }
        
        found = []
        for platform, keywords in platforms.items():
            if any(kw in text for kw in keywords):
                found.append(platform)
        
        return ','.join(found) if found else 'not specified'
    
    def _extract_experiment_design(self, description: str) -> str:
        """提取实验设计信息"""
        text = description.lower()
        
        designs = []
        
        if 'time course' in text or 'temporal' in text:
            designs.append('time_course')
        if 'case control' in text or 'disease vs control' in text:
            designs.append('case_control')
        if 'treatment' in text or 'drug' in text:
            designs.append('treatment')
        if 'developmental' in text or 'differentiation' in text:
            designs.append('developmental')
        if 'spatial' in text:
            designs.append('spatial')
        
        return ','.join(designs) if designs else 'not specified'
    
    def _extract_sample_type(self, description: str) -> str:
        """提取样本类型"""
        text = description.lower()
        
        if 'cell line' in text or 'culture' in text:
            return 'cell_line'
        elif 'organoid' in text:
            return 'organoid'
        elif 'tissue' in text or 'biopsy' in text:
            return 'primary_tissue'
        elif 'pbmc' in text or 'blood' in text:
            return 'blood'
        
        return 'not specified'
    
    def _extract_tissue_location(self, description: str) -> str:
        """提取组织位置细节"""
        text = description.lower()
        
        locations = []
        
        # 脑区
        brain_regions = ['cortex', 'hippocampus', 'cerebellum', 'striatum', 
                        'amygdala', 'thalamus', 'hypothalamus', 'brainstem']
        locations.extend([r for r in brain_regions if r in text])
        
        # 其他位置
        if 'peripheral' in text:
            locations.append('peripheral')
        if 'central' in text:
            locations.append('central')
        
        return ','.join(locations) if locations else 'not specified'
    
    def _extract_pubmed(self, description: str, citation: str) -> str:
        """提取PubMed ID"""
        text = description + ' ' + citation
        
        # 查找PMID模式
        pmid_patterns = [
            r'PMID[:\s]+(\d+)',
            r'PubMed[:\s]+(\d+)',
            r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)'
        ]
        
        for pattern in pmid_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                return match.group(1)
        
        return 'not available'
    
    def save_raw_metadata(self, articles: List[Dict], filename: str = "raw_metadata.json"):
        """保存原始元数据"""
        filepath = self.output_dir / filename
        
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(articles, f, ensure_ascii=False, indent=2)
        
        logger.info(f"原始元数据已保存至: {filepath}")
    
    def create_structured_table(self, articles: List[Dict]) -> pd.DataFrame:
        """创建结构化数据表"""
        structured_data = []
        
        for article in articles:
            try:
                metadata = self.extract_metadata_fields(article)
                structured_data.append(metadata)
            except Exception as e:
                logger.error(f"处理文章 {article.get('id')} 时出错: {str(e)}")
        
        df = pd.DataFrame(structured_data)
        
        # 确保所有需要的列都存在
        required_columns = [
            'id', 'sample_id', 'title', 'disease_general', 'disease', 'pubmed',
            'source_database', 'access_link', 'open_status', 'ethnicity', 'sex',
            'tissue', 'sequencing_platform', 'experiment_design', 'sample_type',
            'summary', 'citation_count', 'publication_date', 'submission_date',
            'last_update_date', 'contact_name', 'contact_email', 'contact_institute',
            'data_tier', 'tissue_location', 'supplementary_information'
        ]
        
        for col in required_columns:
            if col not in df.columns:
                df[col] = ''
        
        return df[required_columns]
    
    def run_collection(self):
        """运行完整的收集流程"""
        logger.info("=" * 60)
        logger.info("开始收集Figshare单细胞测序数据")
        logger.info("=" * 60)
        
        # 步骤1: 收集所有文章
        logger.info("\n步骤 1/3: 收集所有相关文章...")
        all_articles = self.collect_all_articles()
        
        if not all_articles:
            logger.warning("未找到任何文章!")
            return
        
        # 步骤2: 保存原始元数据
        logger.info("\n步骤 2/3: 保存原始元数据...")
        self.save_raw_metadata(all_articles)
        
        # 步骤3: 创建结构化表格
        logger.info("\n步骤 3/3: 创建结构化数据表...")
        structured_df = self.create_structured_table(all_articles)
        
        # 保存为多种格式
        csv_path = self.output_dir / "figshare_scrna_metadata.csv"
        excel_path = self.output_dir / "figshare_scrna_metadata.xlsx"
        tsv_path = self.output_dir / "figshare_scrna_metadata.tsv"
        
        structured_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
        structured_df.to_excel(excel_path, index=False, engine='openpyxl')
        structured_df.to_csv(tsv_path, index=False, sep='\t', encoding='utf-8')
        
        logger.info(f"\n结构化数据已保存:")
        logger.info(f"  - CSV: {csv_path}")
        logger.info(f"  - Excel: {excel_path}")
        logger.info(f"  - TSV: {tsv_path}")
        
        # 生成统计报告
        self.generate_statistics_report(structured_df)
        
        logger.info("\n" + "=" * 60)
        logger.info("数据收集完成!")
        logger.info("=" * 60)
        
        return structured_df
    
    def generate_statistics_report(self, df: pd.DataFrame):
        """生成统计报告"""
        report_path = self.output_dir / "collection_statistics.txt"
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("Figshare单细胞测序数据收集统计报告\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"总数据集数量: {len(df)}\n\n")
            
            # 按疾病统计
            f.write("疾病类型分布:\n")
            disease_counts = df['disease_general'].str.split(',').explode().value_counts()
            for disease, count in disease_counts.head(10).items():
                f.write(f"  {disease}: {count}\n")
            f.write("\n")
            
            # 按组织统计
            f.write("组织类型分布:\n")
            tissue_counts = df['tissue'].str.split(',').explode().value_counts()
            for tissue, count in tissue_counts.head(10).items():
                f.write(f"  {tissue}: {count}\n")
            f.write("\n")
            
            # 按平台统计
            f.write("测序平台分布:\n")
            platform_counts = df['sequencing_platform'].value_counts()
            for platform, count in platform_counts.head(10).items():
                f.write(f"  {platform}: {count}\n")
            f.write("\n")
            
            # 按年份统计
            f.write("发布年份分布:\n")
            df['year'] = pd.to_datetime(df['publication_date'], errors='coerce').dt.year
            year_counts = df['year'].value_counts().sort_index()
            for year, count in year_counts.items():
                if pd.notna(year):
                    f.write(f"  {int(year)}: {count}\n")
            f.write("\n")
            
            # 数据开放性
            f.write("数据开放状态:\n")
            open_counts = df['open_status'].value_counts()
            for status, count in open_counts.items():
                f.write(f"  {status}: {count}\n")
            f.write("\n")
            
            # 数据层级
            f.write("数据层级分布:\n")
            tier_counts = df['data_tier'].value_counts()
            for tier, count in tier_counts.head(10).items():
                f.write(f"  {tier}: {count}\n")
        
        logger.info(f"统计报告已保存至: {report_path}")


def main():
    """主函数"""
    print("""
    ╔═══════════════════════════════════════════════════════════╗
    ║   Figshare scRNA-seq 数据收集系统                          ║
    ║   Human Single-Cell RNA Sequencing Metadata Collector    ║
    ╚═══════════════════════════════════════════════════════════╝
    """)
    
    # 创建收集器实例
    collector = FigshareScRNASeqCollector(output_dir="./figshare_scrna_data")
    
    # 运行收集流程
    df = collector.run_collection()
    
    if df is not None:
        print("\n数据预览:")
        print(df.head())
        print(f"\n总共收集 {len(df)} 个数据集")
        print(f"\n所有文件已保存至: {collector.output_dir}")


if __name__ == "__main__":
    main()