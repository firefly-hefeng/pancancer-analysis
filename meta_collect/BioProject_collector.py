"""
BioProject_collector.py
NCBI BioProject 数据收集器
"""

import os
import json
import time
import requests
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
import logging
import xml.etree.ElementTree as ET
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('bioproject_collection.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class BioProjectCollector:
    """BioProject数据收集器"""
    
    def __init__(self, output_dir="BioProject_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        
    def search_projects(self, search_terms=None, max_results=5000):
        """搜索BioProject"""
        logger.info("开始搜索BioProject...")
        
        if search_terms is None:
            search_terms = [
                '("Homo sapiens"[Organism]) AND ("single cell"[All Fields] OR "scRNA"[All Fields]) AND "RNA-Seq"[All Fields]',
                '("Homo sapiens"[Organism]) AND "single-cell transcriptom*"[All Fields]',
                '("Homo sapiens"[Organism]) AND "10x genomics"[All Fields] AND "RNA"[All Fields]'
            ]
        
        all_ids = set()
        
        for query in search_terms:
            try:
                logger.info(f"执行查询: {query[:80]}...")
                
                search_url = f"{self.base_url}/esearch.fcgi"
                params = {
                    'db': 'bioproject',
                    'term': query,
                    'retmax': max_results,
                    'retmode': 'json',
                    'usehistory': 'y'
                }
                
                response = requests.get(search_url, params=params, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    search_result = data.get('esearchresult', {})
                    count = int(search_result.get('count', 0))
                    id_list = search_result.get('idlist', [])
                    
                    logger.info(f"  找到 {count} 个项目 (返回 {len(id_list)} 个ID)")
                    all_ids.update(id_list)
                    
                    # 保存WebEnv用于大批量下载
                    if count > 0:
                        webenv = search_result.get('webenv')
                        query_key = search_result.get('querykey')
                        
                        # 保存搜索会话信息
                        session_info = {
                            'query': query,
                            'count': count,
                            'webenv': webenv,
                            'query_key': query_key,
                            'timestamp': datetime.now().isoformat()
                        }
                        
                        session_file = self.output_dir / f'search_session_{len(all_ids)}.json'
                        with open(session_file, 'w') as f:
                            json.dump(session_info, f, indent=2)
                
                time.sleep(0.5)
                
            except Exception as e:
                logger.error(f"搜索时出错: {str(e)}", exc_info=True)
        
        logger.info(f"总共找到 {len(all_ids)} 个唯一项目")
        
        # 保存ID列表
        id_file = self.output_dir / 'project_ids.txt'
        with open(id_file, 'w') as f:
            f.write('\n'.join(sorted(all_ids)))
        
        return list(all_ids)
    
    def fetch_project_details(self, project_ids: List[str], batch_size=100):
        """获取项目详细信息"""
        logger.info(f"开始获取 {len(project_ids)} 个项目的详细信息...")
        
        all_projects = []
        failed_ids = []
        
        for i in tqdm(range(0, len(project_ids), batch_size), desc="获取项目详情"):
            batch_ids = project_ids[i:i+batch_size]
            
            try:
                fetch_url = f"{self.base_url}/efetch.fcgi"
                params = {
                    'db': 'bioproject',
                    'id': ','.join(batch_ids),
                    'retmode': 'xml'
                }
                
                response = requests.get(fetch_url, params=params, timeout=60)
                
                if response.status_code == 200:
                    # 保存原始XML
                    xml_file = self.output_dir / f'batch_{i//batch_size}_raw.xml'
                    with open(xml_file, 'wb') as f:
                        f.write(response.content)
                    
                    # 解析XML
                    try:
                        root = ET.fromstring(response.content)
                        
                        for package in root.findall('.//Package'):
                            project_data = self._parse_bioproject_xml(package)
                            if project_data:
                                all_projects.append(project_data)
                        
                    except ET.ParseError as e:
                        logger.error(f"XML解析错误: {str(e)}")
                        failed_ids.extend(batch_ids)
                else:
                    logger.warning(f"批次请求失败: {response.status_code}")
                    failed_ids.extend(batch_ids)
                
                time.sleep(0.5)
                
            except Exception as e:
                logger.error(f"获取批次数据时出错: {str(e)}")
                failed_ids.extend(batch_ids)
        
        logger.info(f"成功获取 {len(all_projects)} 个项目")
        if failed_ids:
            logger.warning(f"失败的项目数: {len(failed_ids)}")
            
            # 保存失败列表
            failed_file = self.output_dir / 'failed_project_ids.txt'
            with open(failed_file, 'w') as f:
                f.write('\n'.join(failed_ids))
        
        # 保存项目数据
        if all_projects:
            # JSON格式
            json_file = self.output_dir / 'bioproject_metadata.json'
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(all_projects, f, indent=2, ensure_ascii=False)
            
            # CSV格式
            df = pd.DataFrame(all_projects)
            csv_file = self.output_dir / 'bioproject_metadata.csv'
            df.to_csv(csv_file, index=False, encoding='utf-8-sig')
            
            logger.info(f"项目数据已保存: {csv_file}")
        
        return all_projects
    
    def _parse_bioproject_xml(self, package_elem) -> Optional[Dict]:
        """解析BioProject XML元素"""
        try:
            data = {}
            
            # 项目ID
            project_id = package_elem.find('.//Project/ProjectID')
            if project_id is not None:
                archive_id = project_id.find('ArchiveID')
                if archive_id is not None:
                    data['accession'] = archive_id.get('accession', '')
                    data['archive'] = archive_id.get('archive', '')
                    data['id'] = archive_id.get('id', '')
            
            # 项目描述
            project_descr = package_elem.find('.//Project/ProjectDescr')
            if project_descr is not None:
                data['title'] = self._get_text(project_descr, 'Title')
                data['description'] = self._get_text(project_descr, 'Description')
                
                # 发表信息
                publication = project_descr.find('Publication')
                if publication is not None:
                    data['pubmed_id'] = publication.get('id', '')
                    data['publication_status'] = publication.get('status', '')
                    data['publication_date'] = publication.get('date', '')
                
                # 相关性
                relevance = project_descr.find('Relevance')
                if relevance is not None:
                    data['relevance_agricultural'] = relevance.findtext('Agricultural', '')
                    data['relevance_medical'] = relevance.findtext('Medical', '')
                    data['relevance_industrial'] = relevance.findtext('Industrial', '')
                    data['relevance_environmental'] = relevance.findtext('Environmental', '')
                    data['relevance_evolution'] = relevance.findtext('Evolution', '')
                    data['relevance_model'] = relevance.findtext('Model', '')
                    data['relevance_other'] = relevance.findtext('Other', '')
            
            # 项目类型
            project_type = package_elem.find('.//Project/ProjectType')
            if project_type is not None:
                data['project_data_type'] = project_type.get('ProjectDataType', '')
                
                # 目标信息
                target = project_type.find('.//Target')
                if target is not None:
                    data['sample_scope'] = target.get('sample_scope', '')
                    data['material'] = target.get('material', '')
                    data['capture'] = target.get('capture', '')
                    
                    organism = target.find('Organism')
                    if organism is not None:
                        data['organism'] = organism.get('species', '')
                        data['taxon_id'] = organism.get('taxID', '')
                        data['organism_label'] = organism.findtext('OrganismName', '')
                
                # 方法信息
                method = project_type.find('.//Method')
                if method is not None:
                    data['method_type'] = method.get('method_type', '')
                
                # 目标信息
                objectives = project_type.find('.//Objectives/Data')
                if objectives is not None:
                    data['data_type'] = objectives.get('data_type', '')
            
            # 提交信息
            submission = package_elem.find('.//Submission')
            if submission is not None:
                data['submission_id'] = submission.get('submission_id', '')
                data['submitted_date'] = submission.get('submitted', '')
                data['last_update'] = submission.get('last_update', '')
                
                # 组织信息
                org = submission.find('.//Organization')
                if org is not None:
                    data['organization_role'] = org.get('role', '')
                    data['organization_type'] = org.get('type', '')
                    
                    name_elem = org.find('.//Name')
                    if name_elem is not None:
                        data['organization_name'] = name_elem.text
                    
                    # 联系人
                    contact = org.find('.//Contact')
                    if contact is not None:
                        data['contact_email'] = contact.get('email', '')
                        name = contact.find('.//Name')
                        if name is not None:
                            data['contact_first'] = name.findtext('First', '')
                            data['contact_last'] = name.findtext('Last', '')
            
            # 外部链接
            external_links = package_elem.findall('.//ExternalLink')
            links = []
            for link in external_links:
                link_data = {
                    'category': link.get('category', ''),
                    'label': link.get('label', ''),
                    'url': link.findtext('.//URL', '')
                }
                links.append(link_data)
            if links:
                data['external_links'] = json.dumps(links)
            
            return data if data.get('accession') else None
            
        except Exception as e:
            logger.error(f"解析XML时出错: {str(e)}")
            return None
    
    def _get_text(self, element, tag, default=''):
        """安全获取XML文本"""
        try:
            text = element.findtext(tag, default)
            return text.strip() if text else default
        except:
            return default
    
    def fetch_linked_sra_data(self, projects: List[Dict]):
        """获取关联的SRA数据"""
        logger.info("开始获取关联的SRA数据...")
        
        project_accs = [p.get('accession') for p in projects if p.get('accession')]
        
        all_sra_links = []
        
        for acc in tqdm(project_accs[:100], desc="查询SRA链接"):  # 限制数量
            try:
                # 使用elink查找关联的SRA数据
                elink_url = f"{self.base_url}/elink.fcgi"
                params = {
                    'dbfrom': 'bioproject',
                    'db': 'sra',
                    'id': acc,
                    'retmode': 'json'
                }
                
                response = requests.get(elink_url, params=params, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    linksets = data.get('linksets', [])
                    
                    for linkset in linksets:
                        sra_ids = linkset.get('linksetdbs', [{}])[0].get('links', [])
                        if sra_ids:
                            all_sra_links.append({
                                'bioproject_accession': acc,
                                'sra_count': len(sra_ids),
                                'sra_ids': ','.join(sra_ids)
                            })
                
                time.sleep(0.4)
                
            except Exception as e:
                logger.error(f"查询 {acc} 的SRA链接时出错: {str(e)}")
        
        # 保存SRA链接信息
        if all_sra_links:
            df = pd.DataFrame(all_sra_links)
            sra_file = self.output_dir / 'bioproject_sra_links.csv'
            df.to_csv(sra_file, index=False, encoding='utf-8-sig')
            logger.info(f"SRA链接已保存: {sra_file}")
        
        return all_sra_links
    
    def generate_summary(self):
        """生成摘要报告"""
        summary_file = self.output_dir / 'collection_summary.txt'
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("BioProject 数据收集摘要\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"收集时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # 读取项目数据
            csv_file = self.output_dir / 'bioproject_metadata.csv'
            if csv_file.exists():
                df = pd.read_csv(csv_file)
                
                f.write(f"项目总数: {len(df)}\n\n")
                
                # 按数据类型统计
                if 'project_data_type' in df.columns:
                    f.write("按项目数据类型:\n")
                    type_counts = df['project_data_type'].value_counts()
                    for dtype, count in type_counts.items():
                        if pd.notna(dtype):
                            f.write(f"  {dtype}: {count}\n")
                    f.write("\n")
                
                # 按生物体统计
                if 'organism' in df.columns:
                    f.write("按生物体:\n")
                    org_counts = df['organism'].value_counts().head(10)
                    for org, count in org_counts.items():
                        if pd.notna(org):
                            f.write(f"  {org}: {count}\n")
                    f.write("\n")
                
                # 按提交年份统计
                if 'submitted_date' in df.columns:
                    df['year'] = pd.to_datetime(df['submitted_date'], errors='coerce').dt.year
                    year_counts = df['year'].value_counts().sort_index()
                    f.write("按提交年份:\n")
                    for year, count in year_counts.items():
                        if pd.notna(year):
                            f.write(f"  {int(year)}: {count}\n")
                    f.write("\n")
                
                # 发表状态统计
                if 'publication_status' in df.columns:
                    f.write("发表状态:\n")
                    pub_counts = df['publication_status'].value_counts()
                    for status, count in pub_counts.items():
                        if pd.notna(status):
                            f.write(f"  {status}: {count}\n")
        
        logger.info(f"摘要已生成: {summary_file}")
    
    def run_full_collection(self, include_sra_links=False):
        """运行完整收集流程"""
        logger.info("=" * 80)
        logger.info("开始BioProject完整数据收集")
        logger.info("=" * 80)
        
        # 1. 搜索项目
        project_ids = self.search_projects()
        
        if not project_ids:
            logger.error("未找到项目,终止收集")
            return
        
        # 2. 获取项目详情
        projects = self.fetch_project_details(project_ids)
        
        # 3. 可选: 获取SRA链接
        if include_sra_links and projects:
            self.fetch_linked_sra_data(projects)
        
        # 4. 生成摘要
        self.generate_summary()
        
        logger.info("=" * 80)
        logger.info("BioProject数据收集完成!")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info("=" * 80)


def main():
    """主函数"""
    print("=" * 80)
    print("NCBI BioProject 单细胞RNA测序数据收集器")
    print("=" * 80)
    print("\n此脚本将收集:")
    print("  1. BioProject元数据")
    print("  2. 项目描述和注释")
    print("  3. 关联的发表信息")
    print("  4. 可选: 关联的SRA数据链接")
    print("=" * 80)
    
    collector = BioProjectCollector()
    
    # 是否包含SRA链接查询
    include_sra = False  # 改为True可查询SRA链接
    
    collector.run_full_collection(include_sra_links=include_sra)


if __name__ == "__main__":
    main()