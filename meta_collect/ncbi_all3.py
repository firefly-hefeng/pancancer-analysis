# ncbi_all_2_improved.py
import os
import json
import logging
import time
import requests
import pandas as pd
from typing import List, Dict, Optional, Set
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from pathlib import Path
import xml.etree.ElementTree as ET
from collections import defaultdict
from threading import Lock
import itertools

class NCBIBioProjectSRACollector:
    """NCBI BioProject和SRA数据收集器 - 多邮箱并发版本"""
    
    def __init__(self, 
                 output_dir: str = 'ncbi_bioproject_sra_data',
                 emails: List[str] = None,
                 workers_per_email: int = 3,
                 api_key: str = None):
        """
        初始化数据收集器
        
        Args:
            output_dir: 输出目录
            emails: 邮箱列表，用于并发请求
            workers_per_email: 每个邮箱的并发线程数
            api_key: NCBI API Key（可选，有key可以提高到10次/秒）
        """
        self.output_dir = output_dir
        self.emails = emails if emails else self._load_or_generate_emails()
        self.workers_per_email = workers_per_email
        self.total_workers = len(self.emails) * workers_per_email
        self.api_key = api_key
        
        # 为每个邮箱创建独立资源
        self.email_resources = {
            email: {
                'session': self._create_session(),
                'lock': Lock(),
                'last_request_time': 0,
                'request_count': 0
            }
            for email in self.emails
        }
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'bioproject_xml'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'sra_xml'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'run_info'), exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'logs'), exist_ok=True)
        
        # 设置日志
        self.logger = self._setup_logger()
        
        # NCBI API设置
        if api_key:
            self.min_interval = 0.1  # 有API key时，10次/秒
        else:
            self.min_interval = 0.34  # 无API key时，3次/秒
            
        self.ncbi_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        self.logger.info(f"初始化完成: {len(self.emails)} 个邮箱, "
                        f"每邮箱 {workers_per_email} 线程, "
                        f"总并发: {self.total_workers}")
        
    def _load_or_generate_emails(self) -> List[str]:
        """加载或生成虚拟邮箱"""
        email_file = 'virtual_emails.txt'
        
        if os.path.exists(email_file):
            with open(email_file, 'r') as f:
                emails = [line.strip() for line in f if line.strip()]
                if emails:
                    return emails
        
        # 生成默认邮箱
        return self._generate_default_emails(10)
    
    def _generate_default_emails(self, count: int = 10) -> List[str]:
        """生成默认虚拟邮箱"""
        import random
        import string
        
        emails = []
        prefixes = ['bioproject', 'sra', 'ncbi', 'genomics', 'research',
                   'bioinf', 'seq', 'data', 'science', 'lab']
        
        for i in range(count):
            prefix = random.choice(prefixes)
            suffix = ''.join(random.choices(string.ascii_lowercase + string.digits, k=6))
            email = f"{prefix}.{suffix}@research.tempmail.org"
            emails.append(email)
        
        # 保存到文件
        with open('virtual_emails.txt', 'w') as f:
            for email in emails:
                f.write(email + '\n')
        
        return emails
    
    def _setup_logger(self) -> logging.Logger:
        """设置日志记录器"""
        logger = logging.getLogger('NCBICollector')
        logger.setLevel(logging.INFO)
        
        # 文件处理器
        fh = logging.FileHandler(
            os.path.join(self.output_dir, 'logs', f'collection_{datetime.now():%Y%m%d_%H%M%S}.log'),
            encoding='utf-8'
        )
        fh.setLevel(logging.INFO)
        
        # 控制台处理器
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        
        # 格式化器
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        
        logger.addHandler(fh)
        logger.addHandler(ch)
        
        return logger
    
    def _create_session(self) -> requests.Session:
        """创建带重试和连接池的Session"""
        session = requests.Session()
        
        # 配置重试策略
        retry_strategy = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["HEAD", "GET", "OPTIONS", "POST"]
        )
        
        adapter = HTTPAdapter(
            max_retries=retry_strategy,
            pool_connections=100,
            pool_maxsize=100,
            pool_block=False
        )
        
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        
        return session
    
    def _get_email_for_task(self, task_index: int) -> str:
        """为任务分配邮箱"""
        return self.emails[task_index % len(self.emails)]
    
    def _make_request(self, email: str, url: str, params: dict = None, 
                     max_retries: int = 3) -> Optional[requests.Response]:
        """
        使用指定邮箱发起请求
        
        Args:
            email: 使用的邮箱
            url: 请求URL
            params: 请求参数
            max_retries: 最大重试次数
        
        Returns:
            响应对象或None
        """
        resources = self.email_resources[email]
        
        if params is None:
            params = {}
        
        params['email'] = email
        if self.api_key:
            params['api_key'] = self.api_key
        
        for attempt in range(max_retries):
            try:
                # 控制请求频率
                with resources['lock']:
                    current_time = time.time()
                    elapsed = current_time - resources['last_request_time']
                    
                    if elapsed < self.min_interval:
                        time.sleep(self.min_interval - elapsed)
                    
                    response = resources['session'].get(url, params=params, timeout=30)
                    resources['last_request_time'] = time.time()
                    resources['request_count'] += 1
                
                if response.status_code == 200:
                    return response
                elif response.status_code == 429:
                    wait_time = (2 ** attempt) * 2
                    self.logger.warning(f"邮箱 {email} 速率限制，等待 {wait_time}秒")
                    time.sleep(wait_time)
                else:
                    self.logger.error(f"请求失败 [{email}]: {response.status_code}")
                    
            except requests.exceptions.Timeout:
                self.logger.warning(f"请求超时 [{email}] (尝试 {attempt + 1}/{max_retries})")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
            except Exception as e:
                self.logger.error(f"请求异常 [{email}] (尝试 {attempt + 1}/{max_retries}): {str(e)}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
        
        return None
    
    def search_bioproject(self, email: str, query: str, max_results: int = 10000) -> List[str]:
        """
        搜索BioProject ID
        
        Args:
            email: 使用的邮箱
            query: 搜索查询
            max_results: 最大结果数
        
        Returns:
            BioProject ID列表
        """
        self.logger.info(f"搜索BioProject [{email}]: {query}")
        
        # 第一步：搜索获取总数
        search_url = f"{self.ncbi_base_url}esearch.fcgi"
        search_params = {
            'db': 'bioproject',
            'term': query,
            'retmode': 'json',
            'retmax': 0  # 只获取总数
        }
        
        response = self._make_request(email, search_url, search_params)
        if not response:
            self.logger.error("搜索失败")
            return []
        
        try:
            data = response.json()
            total_count = int(data['esearchresult']['count'])
            self.logger.info(f"找到 {total_count} 个BioProject")
            
            if total_count == 0:
                return []
            
            # 限制最大结果数
            total_count = min(total_count, max_results)
            
            # 分批获取ID
            batch_size = 500
            all_ids = []
            
            for start in range(0, total_count, batch_size):
                search_params['retstart'] = start
                search_params['retmax'] = min(batch_size, total_count - start)
                
                response = self._make_request(email, search_url, search_params)
                if response:
                    batch_data = response.json()
                    batch_ids = batch_data['esearchresult']['idlist']
                    all_ids.extend(batch_ids)
                    self.logger.info(f"已获取 {len(all_ids)}/{total_count} 个ID")
            
            return all_ids
            
        except Exception as e:
            self.logger.error(f"解析搜索结果失败: {str(e)}")
            return []
    
    def fetch_bioproject_data(self, email: str, bioproject_id: str) -> Optional[Dict]:
        """
        获取单个BioProject的详细信息
        
        Args:
            email: 使用的邮箱
            bioproject_id: BioProject ID
        
        Returns:
            BioProject数据字典
        """
        try:
            # 获取BioProject XML
            fetch_url = f"{self.ncbi_base_url}efetch.fcgi"
            fetch_params = {
                'db': 'bioproject',
                'id': bioproject_id,
                'retmode': 'xml'
            }
            
            response = self._make_request(email, fetch_url, fetch_params)
            if not response:
                return None
            
            # 保存XML
            xml_file = os.path.join(self.output_dir, 'bioproject_xml', f'{bioproject_id}.xml')
            with open(xml_file, 'w', encoding='utf-8') as f:
                f.write(response.text)
            
            # 解析XML
            root = ET.fromstring(response.text)
            
            # 提取基本信息
            project_data = {
                'bioproject_id': bioproject_id,
                'accession': None,
                'title': None,
                'description': None,
                'organism': None,
                'submission_date': None,
                'publication_date': None,
                'sra_experiments': []
            }
            
            # 提取Accession
            project_acc = root.find('.//Project/ProjectID/ArchiveID')
            if project_acc is not None:
                project_data['accession'] = project_acc.get('accession')
            
            # 提取标题
            title = root.find('.//Project/ProjectDescr/Title')
            if title is not None:
                project_data['title'] = title.text
            
            # 提取描述
            description = root.find('.//Project/ProjectDescr/Description')
            if description is not None:
                project_data['description'] = description.text
            
            # 提取物种
            organism = root.find('.//Project/ProjectType/ProjectTypeSubmission/Target/Organism/OrganismName')
            if organism is not None:
                project_data['organism'] = organism.text
            
            # 提取日期
            submission = root.find('.//Submission')
            if submission is not None:
                project_data['submission_date'] = submission.get('submitted')
                project_data['publication_date'] = submission.get('last_update')
            
            return project_data
            
        except Exception as e:
            self.logger.error(f"获取BioProject {bioproject_id} 失败 [{email}]: {str(e)}")
            return None
    
    def link_bioproject_to_sra(self, email: str, bioproject_id: str) -> List[str]:
        """
        查找BioProject关联的SRA实验
        
        Args:
            email: 使用的邮箱
            bioproject_id: BioProject ID
        
        Returns:
            SRA实验ID列表
        """
        try:
            # 使用elink查找关联
            link_url = f"{self.ncbi_base_url}elink.fcgi"
            link_params = {
                'dbfrom': 'bioproject',
                'db': 'sra',
                'id': bioproject_id,
                'retmode': 'json'
            }
            
            response = self._make_request(email, link_url, link_params)
            if not response:
                return []
            
            data = response.json()
            
            # 提取SRA ID
            sra_ids = []
            if 'linksets' in data and len(data['linksets']) > 0:
                linkset = data['linksets'][0]
                if 'linksetdbs' in linkset:
                    for linksetdb in linkset['linksetdbs']:
                        if linksetdb.get('dbto') == 'sra':
                            sra_ids = linksetdb.get('links', [])
                            break
            
            return sra_ids
            
        except Exception as e:
            self.logger.error(f"查找BioProject {bioproject_id} 的SRA关联失败 [{email}]: {str(e)}")
            return []
    
    def fetch_sra_data(self, email: str, sra_id: str) -> Optional[Dict]:
        """
        获取SRA实验详细信息
        
        Args:
            email: 使用的邮箱
            sra_id: SRA ID
        
        Returns:
            SRA数据字典
        """
        try:
            # 获取SRA XML
            fetch_url = f"{self.ncbi_base_url}efetch.fcgi"
            fetch_params = {
                'db': 'sra',
                'id': sra_id,
                'retmode': 'xml'
            }
            
            response = self._make_request(email, fetch_url, fetch_params)
            if not response:
                return None
            
            # 保存XML
            xml_file = os.path.join(self.output_dir, 'sra_xml', f'{sra_id}.xml')
            with open(xml_file, 'w', encoding='utf-8') as f:
                f.write(response.text)
            
            # 解析XML
            root = ET.fromstring(response.text)
            
            sra_data = {
                'sra_id': sra_id,
                'experiment_accession': None,
                'study_accession': None,
                'sample_accession': None,
                'run_accessions': [],
                'title': None,
                'library_strategy': None,
                'library_source': None,
                'library_selection': None,
                'platform': None,
                'instrument_model': None
            }
            
            # 提取实验信息
            experiment = root.find('.//EXPERIMENT')
            if experiment is not None:
                sra_data['experiment_accession'] = experiment.get('accession')
                
                title = experiment.find('.//TITLE')
                if title is not None:
                    sra_data['title'] = title.text
                
                design = experiment.find('.//DESIGN')
                if design is not None:
                    lib_desc = design.find('.//LIBRARY_DESCRIPTOR')
                    if lib_desc is not None:
                        strategy = lib_desc.find('.//LIBRARY_STRATEGY')
                        if strategy is not None:
                            sra_data['library_strategy'] = strategy.text
                        
                        source = lib_desc.find('.//LIBRARY_SOURCE')
                        if source is not None:
                            sra_data['library_source'] = source.text
                        
                        selection = lib_desc.find('.//LIBRARY_SELECTION')
                        if selection is not None:
                            sra_data['library_selection'] = selection.text
                
                platform = experiment.find('.//PLATFORM')
                if platform is not None:
                    for child in platform:
                        sra_data['platform'] = child.tag
                        model = child.find('.//INSTRUMENT_MODEL')
                        if model is not None:
                            sra_data['instrument_model'] = model.text
            
            # 提取Study信息
            study = root.find('.//STUDY')
            if study is not None:
                sra_data['study_accession'] = study.get('accession')
            
            # 提取Sample信息
            sample = root.find('.//SAMPLE')
            if sample is not None:
                sra_data['sample_accession'] = sample.get('accession')
            
            # 提取Run信息
            runs = root.findall('.//RUN')
            sra_data['run_accessions'] = [run.get('accession') for run in runs]
            
            return sra_data
            
        except Exception as e:
            self.logger.error(f"获取SRA {sra_id} 失败 [{email}]: {str(e)}")
            return None
    
    def fetch_run_info(self, email: str, run_accession: str) -> Optional[Dict]:
        """
        获取Run的详细信息
        
        Args:
            email: 使用的邮箱
            run_accession: Run accession (如 SRR123456)
        
        Returns:
            Run信息字典
        """
        try:
            # 使用SRA Run Selector API
            url = "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new"
            params = {
                'acc': run_accession,
                'format': 'json'
            }
            
            response = self._make_request(email, url, params)
            if not response:
                return None
            
            data = response.json()
            
            # 保存JSON
            json_file = os.path.join(self.output_dir, 'run_info', f'{run_accession}.json')
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2)
            
            return data
            
        except Exception as e:
            self.logger.error(f"获取Run {run_accession} 信息失败 [{email}]: {str(e)}")
            return None
    
    def process_bioproject_batch(self, bioproject_ids: List[str]) -> List[Dict]:
        """
        批量处理BioProject - 多线程并发
        
        Args:
            bioproject_ids: BioProject ID列表
        
        Returns:
            处理结果列表
        """
        self.logger.info(f"开始批量处理 {len(bioproject_ids)} 个BioProject")
        
        results = []
        completed = 0
        failed = 0
        
        def process_single_bioproject(index: int, bioproject_id: str):
            """处理单个BioProject"""
            email = self._get_email_for_task(index)
            
            try:
                # 获取BioProject数据
                bioproject_data = self.fetch_bioproject_data(email, bioproject_id)
                if not bioproject_data:
                    return None
                
                # 查找关联的SRA
                sra_ids = self.link_bioproject_to_sra(email, bioproject_id)
                
                # 获取SRA详细信息
                sra_data_list = []
                for sra_id in sra_ids[:10]:  # 限制每个项目最多10个SRA
                    sra_data = self.fetch_sra_data(email, sra_id)
                    if sra_data:
                        sra_data_list.append(sra_data)
                
                bioproject_data['sra_experiments'] = sra_data_list
                bioproject_data['sra_count'] = len(sra_ids)
                
                return bioproject_data
                
            except Exception as e:
                self.logger.error(f"处理 {bioproject_id} 失败 [{email}]: {str(e)}")
                return None
        
        # 使用线程池并发处理
        with ThreadPoolExecutor(max_workers=self.total_workers) as executor:
            # 提交所有任务
            futures = {
                executor.submit(process_single_bioproject, idx, bp_id): bp_id
                for idx, bp_id in enumerate(bioproject_ids)
            }
            
            # 收集结果
            for future in as_completed(futures):
                bioproject_id = futures[future]
                
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                        completed += 1
                    else:
                        failed += 1
                except Exception as e:
                    self.logger.error(f"任务执行失败 {bioproject_id}: {str(e)}")
                    failed += 1
                
                # 进度报告
                if (completed + failed) % 50 == 0:
                    self.logger.info(f"进度: {completed + failed}/{len(bioproject_ids)} "
                                   f"(成功: {completed}, 失败: {failed})")
        
        self.logger.info(f"批量处理完成: 成功 {completed}, 失败 {failed}")
        
        return results
    
    def save_results(self, results: List[Dict], prefix: str = 'bioproject'):
        """
        保存结果到多种格式
        
        Args:
            results: 结果列表
            prefix: 文件名前缀
        """
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # 保存JSON
        json_file = os.path.join(self.output_dir, f'{prefix}_data_{timestamp}.json')
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        self.logger.info(f"已保存JSON: {json_file}")
        
        # 准备DataFrame数据
        df_data = []
        for item in results:
            row = {
                'bioproject_id': item.get('bioproject_id'),
                'accession': item.get('accession'),
                'title': item.get('title'),
                'description': item.get('description'),
                'organism': item.get('organism'),
                'submission_date': item.get('submission_date'),
                'publication_date': item.get('publication_date'),
                'sra_count': item.get('sra_count', 0)
            }
            
            # 添加第一个SRA实验的信息（如果存在）
            if item.get('sra_experiments'):
                sra = item['sra_experiments'][0]
                row.update({
                    'first_experiment_acc': sra.get('experiment_accession'),
                    'first_study_acc': sra.get('study_accession'),
                    'library_strategy': sra.get('library_strategy'),
                    'platform': sra.get('platform'),
                    'instrument_model': sra.get('instrument_model')
                })
            
            df_data.append(row)
        
        # 保存Excel
        df = pd.DataFrame(df_data)
        excel_file = os.path.join(self.output_dir, f'{prefix}_summary_{timestamp}.xlsx')
        df.to_excel(excel_file, index=False, engine='openpyxl')
        self.logger.info(f"已保存Excel: {excel_file}")
        
        # 保存CSV
        csv_file = os.path.join(self.output_dir, f'{prefix}_summary_{timestamp}.csv')
        df.to_csv(csv_file, index=False, encoding='utf-8-sig')
        self.logger.info(f"已保存CSV: {csv_file}")
        
        return json_file, excel_file, csv_file
    
    def get_statistics(self) -> Dict:
        """获取收集器统计信息"""
        stats = {
            'total_emails': len(self.emails),
            'total_workers': self.total_workers,
            'email_usage': {}
        }
        
        for email, resources in self.email_resources.items():
            stats['email_usage'][email] = resources['request_count']
        
        return stats
    
    def print_statistics(self):
        """打印统计信息"""
        stats = self.get_statistics()
        
        print("\n" + "="*60)
        print("收集器统计信息")
        print("="*60)
        print(f"邮箱总数: {stats['total_emails']}")
        print(f"并发线程数: {stats['total_workers']}")
        print(f"\n各邮箱请求次数:")
        
        for email, count in sorted(stats['email_usage'].items(), 
                                   key=lambda x: x[1], reverse=True):
            print(f"  {email}: {count} 次")
        
        total_requests = sum(stats['email_usage'].values())
        print(f"\n总请求次数: {total_requests}")
        print("="*60 + "\n")


def main():
    """主函数"""
    
    # 读取或生成虚拟邮箱
    email_file = 'virtual_emails.txt'
    
    if os.path.exists(email_file):
        with open(email_file, 'r') as f:
            emails = [line.strip() for line in f if line.strip()]
        print(f"从文件加载了 {len(emails)} 个邮箱")
    else:
        print("首次运行，生成虚拟邮箱...")
        from email_generator import generate_virtual_emails
        emails = generate_virtual_emails(10)
        print(f"生成了 {len(emails)} 个虚拟邮箱")
    
    # 创建收集器
    collector = NCBIBioProjectSRACollector(
        output_dir='ncbi_bioproject_sra_data',
        emails=emails,
        workers_per_email=3,  # 每个邮箱3个线程
        api_key=None  # 如果有NCBI API key，填入这里
    )
    
    # 示例1: 搜索特定物种的BioProject
    query = 'Homo sapiens[Organism] AND biomol rna[Properties] AND strategy rna seq[Properties]'
    
    # 使用第一个邮箱进行搜索
    bioproject_ids = collector.search_bioproject(emails[0], query, max_results=1000)
    
    if bioproject_ids:
        print(f"\n找到 {len(bioproject_ids)} 个BioProject")
        
        # 批量处理（使用多邮箱并发）
        results = collector.process_bioproject_batch(bioproject_ids)
        
        # 保存结果
        if results:
            json_file, excel_file, csv_file = collector.save_results(results)
            print(f"\n结果已保存:")
            print(f"  JSON: {json_file}")
            print(f"  Excel: {excel_file}")
            print(f"  CSV: {csv_file}")
        
        # 打印统计信息
        collector.print_statistics()
    
    else:
        print("未找到BioProject")


if __name__ == "__main__":
    main()

