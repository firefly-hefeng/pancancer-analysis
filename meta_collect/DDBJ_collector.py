#!/usr/bin/env python3
"""
DDBJ Human Single-Cell RNA-seq Metadata Collector
Author: Bioinformatics Backend Engineer
Date: 2025-12-06
Python: 3.9+

完整采集 DDBJ 侧所有人源单细胞 RNA-seq 项目元数据（含 GEA）
包含完整类型注解、日志、进度条、自动重试机制
"""

import argparse
import csv
import json
import logging
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
from urllib.parse import urlencode, quote
import urllib3

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_fixed, retry_if_exception_type
from tqdm import tqdm

# 禁用 SSL 警告
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


@dataclass
class ScRNAProject:
    """单细胞 RNA-seq 项目元数据结构"""
    id: str
    sample_id: str
    title: str
    disease_general: str
    disease: str
    pubmed: str
    source_database: str
    access_link: str
    open_status: str
    ethnicity: str
    sex: str
    tissue: str
    sequencing_platform: str
    experiment_design: str
    sample_type: str
    summary: str
    citation_count: int
    publication_date: str
    submission_date: str
    last_update_date: str
    contact_name: str
    contact_email: str
    contact_institute: str
    data_tier: str
    tissue_location: str
    supplementary_information: str


class DDBJCollector:
    """DDBJ 数据采集器"""
    
    ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"
    DDBJ_API_BASE = "https://ddbj.nig.ac.jp/search"  # 备用 URL
    DDBJ_FTP_BASE = "https://ddbj.nig.ac.jp/public/ddbj_database/gea"
    NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    ENA_DELAY = 0.35  # 3 req/s
    DDBJ_DELAY = 0.55  # 2 req/s
    NCBI_DELAY = 0.35
    
    SC_KEYWORDS = [
        "single cell", "scrna-seq", "single-cell", 
        "single nucleus", "snrna-seq", "sc rna", "sn rna",
        "scrnaseq", "snrnaseq", "10x", "drop-seq", "smart-seq"
    ]
    
    def __init__(self, output_dir: Path, email: str = "user@example.com"):
        self.output_dir = output_dir
        self.email = email
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": f"DDBJ-ScRNA-Collector/1.0 ({email})"
        })
        # 禁用 SSL 验证以解决 DDBJ API 的 SSL 问题
        self.session.verify = False
        
        self._setup_logging()
        
        self.bioproject_cache: Dict[str, Dict] = {}
        self.biosample_cache: Dict[str, Dict] = {}
        self.pubmed_cache: Dict[str, int] = {}
        
    def _setup_logging(self) -> None:
        """配置日志系统"""
        log_file = self.output_dir / "ddbj_collection.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"DDBJ Collector initialized. Output: {self.output_dir}")
    
    def _is_single_cell(self, text: str) -> bool:
        """判断文本是否包含单细胞关键词"""
        text_lower = text.lower()
        return any(kw in text_lower for kw in self.SC_KEYWORDS)
    
    @retry(
        stop=stop_after_attempt(3), 
        wait=wait_fixed(2),
        retry=retry_if_exception_type((requests.exceptions.RequestException, requests.exceptions.HTTPError))
    )
    def _fetch_ena_studies_search(self) -> List[str]:
        """搜索 ENA 中所有 PRJDB 项目ID"""
        params = {
            "result": "study",
            "query": 'tax_eq(9606) AND study_accession="PRJDB*"',
            "fields": "study_accession",
            "format": "tsv",
            "limit": 0  # 获取总数
        }
        
        time.sleep(self.ENA_DELAY)
        response = self.session.get(self.ENA_BASE, params=params, timeout=30)
        
        if response.status_code == 204:
            self.logger.warning("No PRJDB studies found")
            return []
        
        response.raise_for_status()
        
        lines = response.text.strip().split('\n')
        if len(lines) <= 1:
            return []
        
        study_ids = [line.strip() for line in lines[1:] if line.strip()]
        return study_ids
    
    @retry(
        stop=stop_after_attempt(3), 
        wait=wait_fixed(2),
        retry=retry_if_exception_type((requests.exceptions.RequestException,))
    )
    def _fetch_ena_study_detail(self, study_id: str) -> Optional[Dict]:
        """获取单个 ENA study 的详细信息"""
        params = {
            "result": "study",
            "query": f'study_accession="{study_id}"',
            "fields": "study_accession,study_title,study_alias,secondary_study_accession,description,center_name,first_public,last_updated",
            "format": "json"
        }
        
        time.sleep(self.ENA_DELAY)
        response = self.session.get(self.ENA_BASE, params=params, timeout=30)
        
        if response.status_code == 204:
            return None
        
        try:
            response.raise_for_status()
            data = response.json()
            if isinstance(data, list) and len(data) > 0:
                return data[0]
        except (requests.exceptions.HTTPError, json.JSONDecodeError) as e:
            self.logger.warning(f"Failed to fetch detail for {study_id}: {e}")
        
        return None
    
    def _fetch_ddbj_gea_via_browser(self) -> List[Dict]:
        """通过 DDBJ 浏览器页面爬取 GEA 数据（备用方案）"""
        self.logger.info("Using DDBJ browser scraping method as fallback...")
        
        try:
            # 尝试访问 GEA 主页
            url = "https://ddbj.nig.ac.jp/gea/browse.html"
            time.sleep(self.DDBJ_DELAY)
            
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            
            # 这里需要解析 HTML，但由于页面可能是动态加载的，
            # 我们返回空列表，让主程序继续处理 ENA 数据
            self.logger.warning("DDBJ GEA browser scraping not fully implemented")
            return []
            
        except Exception as e:
            self.logger.warning(f"Failed to scrape DDBJ GEA browser: {e}")
            return []
    
    def _fetch_ddbj_gea_from_ena(self, limit: int = 1000) -> List[Dict]:
        """从 ENA 获取 E-GEAD 实验（更可靠的方法）"""
        self.logger.info("Fetching DDBJ GEA experiments from ENA...")
        
        params = {
            "result": "study",
            "query": 'tax_eq(9606) AND study_accession="E-GEAD*"',
            "fields": "study_accession,study_title,study_alias,secondary_study_accession,description,center_name,first_public,last_updated",
            "format": "json",
            "limit": limit
        }
        
        try:
            time.sleep(self.ENA_DELAY)
            response = self.session.get(self.ENA_BASE, params=params, timeout=30)
            
            if response.status_code == 204:
                self.logger.warning("No E-GEAD studies found in ENA")
                return []
            
            response.raise_for_status()
            data = response.json()
            
            if isinstance(data, list):
                # 过滤单细胞实验
                sc_experiments = []
                for exp in data:
                    title = exp.get("study_title", "")
                    description = exp.get("description", "")
                    if self._is_single_cell(title + " " + description):
                        sc_experiments.append(exp)
                
                self.logger.info(f"Found {len(sc_experiments)} single-cell E-GEAD experiments from ENA")
                return sc_experiments
            
            return []
            
        except Exception as e:
            self.logger.warning(f"Failed to fetch E-GEAD from ENA: {e}")
            return []
    
    @retry(
        stop=stop_after_attempt(3), 
        wait=wait_fixed(2),
        retry=retry_if_exception_type((requests.exceptions.RequestException,))
    )
    def _fetch_bioproject_xml(self, project_id: str) -> Optional[Dict]:
        """获取 BioProject XML 并转为 JSON"""
        if project_id in self.bioproject_cache:
            return self.bioproject_cache[project_id]
        
        params = {
            "db": "bioproject",
            "id": project_id,
            "retmode": "xml"
        }
        
        time.sleep(self.NCBI_DELAY)
        url = f"{self.NCBI_EUTILS}/efetch.fcgi"
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            root = ET.fromstring(response.content)
            project_data = self._parse_bioproject_xml(root)
            self.bioproject_cache[project_id] = project_data
            return project_data
        except (requests.exceptions.RequestException, ET.ParseError) as e:
            self.logger.warning(f"Failed to fetch/parse BioProject {project_id}: {e}")
            return None
    
    def _parse_bioproject_xml(self, root: ET.Element) -> Dict:
        """解析 BioProject XML"""
        data: Dict[str, Any] = {
            "description": "",
            "organism": "",
            "publication": [],
            "submission_date": "",
            "relevance": ""
        }
        
        package = root.find(".//Package")
        if package is None:
            return data
        
        project = package.find(".//Project")
        if project is not None:
            desc_elem = project.find(".//Description")
            if desc_elem is not None:
                title = desc_elem.find("Title")
                if title is not None and title.text:
                    data["description"] = title.text
                
                desc_text = desc_elem.find("Description")
                if desc_text is not None and desc_text.text:
                    data["description"] += " " + desc_text.text
            
            relevance = project.find(".//ProjectType/ProjectTypeSubmission/Target/Organism")
            if relevance is not None:
                org_name = relevance.get("species", "")
                data["organism"] = org_name
        
        pub_list = package.findall(".//Publication")
        for pub in pub_list:
            pub_id = pub.get("id", "")
            if pub_id and pub_id.isdigit():
                data["publication"].append(pub_id)
        
        submission = package.find(".//Submission")
        if submission is not None:
            sub_date = submission.get("submitted", "")
            data["submission_date"] = sub_date
        
        return data
    
    @retry(
        stop=stop_after_attempt(3), 
        wait=wait_fixed(2),
        retry=retry_if_exception_type((requests.exceptions.RequestException,))
    )
    def _fetch_biosample_xml(self, sample_id: str) -> Optional[Dict]:
        """获取 BioSample XML 并转为 JSON"""
        if sample_id in self.biosample_cache:
            return self.biosample_cache[sample_id]
        
        params = {
            "db": "biosample",
            "id": sample_id,
            "retmode": "xml"
        }
        
        time.sleep(self.NCBI_DELAY)
        url = f"{self.NCBI_EUTILS}/efetch.fcgi"
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            root = ET.fromstring(response.content)
            sample_data = self._parse_biosample_xml(root)
            self.biosample_cache[sample_id] = sample_data
            return sample_data
        except (requests.exceptions.RequestException, ET.ParseError) as e:
            self.logger.warning(f"Failed to fetch/parse BioSample {sample_id}: {e}")
            return None
    
    def _parse_biosample_xml(self, root: ET.Element) -> Dict:
        """解析 BioSample XML"""
        data: Dict[str, Any] = {
            "attributes": {},
            "organism": "",
            "tissue": "",
            "sex": "",
            "disease": "",
            "ethnicity": ""
        }
        
        biosample = root.find(".//BioSample")
        if biosample is None:
            return data
        
        description = biosample.find(".//Description")
        if description is not None:
            org = description.find(".//Organism")
            if org is not None:
                data["organism"] = org.get("taxonomy_name", "")
        
        attributes = biosample.find(".//Attributes")
        if attributes is not None:
            for attr in attributes.findall("Attribute"):
                attr_name = attr.get("attribute_name", "")
                attr_value = attr.text or ""
                data["attributes"][attr_name] = attr_value
                
                attr_lower = attr_name.lower()
                
                if attr_lower in ["tissue", "tissue_type", "source_name", "cell_type"]:
                    data["tissue"] = attr_value
                elif attr_lower in ["sex", "gender"]:
                    data["sex"] = attr_value
                elif attr_lower in ["disease", "disease_state", "health_state", "disease_status"]:
                    data["disease"] = attr_value
                elif attr_lower in ["ethnicity", "race", "population"]:
                    data["ethnicity"] = attr_value
        
        return data
    
    @retry(
        stop=stop_after_attempt(3), 
        wait=wait_fixed(2),
        retry=retry_if_exception_type((requests.exceptions.RequestException,))
    )
    def _fetch_sra_runs(self, project_id: str, max_runs: int = 50) -> List[Dict]:
        """获取关联的 SRA/DRA Runs"""
        params = {
            "db": "sra",
            "term": f"{project_id}[BioProject]",
            "retmax": max_runs,
            "retmode": "xml"
        }
        
        time.sleep(self.NCBI_DELAY)
        search_url = f"{self.NCBI_EUTILS}/esearch.fcgi"
        
        try:
            search_response = self.session.get(search_url, params=params, timeout=30)
            search_response.raise_for_status()
            
            search_root = ET.fromstring(search_response.content)
            id_list = search_root.findall(".//Id")
            sra_ids = [id_elem.text for id_elem in id_list if id_elem.text]
            
            if not sra_ids:
                return []
            
            fetch_params = {
                "db": "sra",
                "id": ",".join(sra_ids[:max_runs]),
                "retmode": "xml"
            }
            
            time.sleep(self.NCBI_DELAY)
            fetch_url = f"{self.NCBI_EUTILS}/efetch.fcgi"
            fetch_response = self.session.get(fetch_url, params=fetch_params, timeout=30)
            fetch_response.raise_for_status()
            
            fetch_root = ET.fromstring(fetch_response.content)
            runs = []
            
            for exp_pkg in fetch_root.findall(".//EXPERIMENT_PACKAGE"):
                run_elem = exp_pkg.find(".//RUN")
                if run_elem is not None:
                    run_acc = run_elem.get("accession", "")
                    platform = exp_pkg.find(".//PLATFORM")
                    platform_name = ""
                    if platform is not None:
                        for child in platform:
                            platform_name = child.tag
                            break
                    
                    runs.append({
                        "run_accession": run_acc,
                        "platform": platform_name
                    })
            
            return runs[:max_runs]
        
        except (requests.exceptions.RequestException, ET.ParseError) as e:
            self.logger.warning(f"Failed to fetch SRA runs for {project_id}: {e}")
            return []
    
    @retry(
        stop=stop_after_attempt(3), 
        wait=wait_fixed(2),
        retry=retry_if_exception_type((requests.exceptions.RequestException,))
    )
    def _fetch_pubmed_citations(self, pubmed_id: str) -> int:
        """获取 PubMed 文献被引次数"""
        if pubmed_id in self.pubmed_cache:
            return self.pubmed_cache[pubmed_id]
        
        params = {
            "dbfrom": "pubmed",
            "id": pubmed_id,
            "linkname": "pubmed_pubmed_citedin",
            "retmode": "json"
        }
        
        time.sleep(self.NCBI_DELAY)
        url = f"{self.NCBI_EUTILS}/elink.fcgi"
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            linksets = data.get("linksets", [])
            if linksets and len(linksets) > 0:
                linksetdbs = linksets[0].get("linksetdbs", [])
                if linksetdbs and len(linksetdbs) > 0:
                    links = linksetdbs[0].get("links", [])
                    count = len(links)
                    self.pubmed_cache[pubmed_id] = count
                    return count
        except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
            self.logger.warning(f"Failed to fetch PubMed citations for {pubmed_id}: {e}")
        
        self.pubmed_cache[pubmed_id] = 0
        return 0
    
    def collect_ena_projects(self) -> List[Dict]:
        """采集所有 ENA PRJDB 项目"""
        self.logger.info("Starting ENA PRJDB project collection...")
        
        study_ids = self._fetch_ena_studies_search()
        self.logger.info(f"Found {len(study_ids)} PRJDB studies in ENA")
        
        all_projects = []
        
        for study_id in tqdm(study_ids, desc="Fetching ENA study details"):
            study_detail = self._fetch_ena_study_detail(study_id)
            
            if study_detail:
                title = study_detail.get("study_title", "")
                description = study_detail.get("description", "")
                
                if self._is_single_cell(title + " " + description):
                    all_projects.append(study_detail)
        
        self.logger.info(f"Collected {len(all_projects)} single-cell projects from ENA")
        
        if all_projects:
            df = pd.DataFrame(all_projects)
            df.to_csv(self.output_dir / "raw_ddbj_bioproject.csv", index=False, encoding='utf-8-sig')
        
        return all_projects
    
    def collect_ddbj_gea_experiments(self) -> List[Dict]:
        """采集所有 DDBJ E-GEAD 实验（使用 ENA 作为数据源）"""
        self.logger.info("Starting DDBJ GEA experiment collection...")
        
        # 优先从 ENA 获取 E-GEAD 数据
        experiments = self._fetch_ddbj_gea_from_ena(limit=10000)
        
        self.logger.info(f"Collected {len(experiments)} single-cell E-GEAD experiments")
        
        if experiments:
            df = pd.DataFrame(experiments)
            df.to_csv(self.output_dir / "raw_ddbj_gea_experiments.csv", index=False, encoding='utf-8-sig')
        
        return experiments
    
    def enrich_project_metadata(self, projects: List[Dict]) -> List[ScRNAProject]:
        """丰富项目元数据"""
        self.logger.info(f"Enriching metadata for {len(projects)} projects...")
        enriched_projects = []
        
        for project in tqdm(projects, desc="Enriching ENA metadata"):
            project_id = project.get("study_accession", "")
            
            bioproject_data = self._fetch_bioproject_xml(project_id)
            
            runs = self._fetch_sra_runs(project_id, max_runs=50)
            platform = runs[0]["platform"] if runs else ""
            
            if runs:
                df_runs = pd.DataFrame(runs)
                self._append_to_csv(self.output_dir / "raw_ddbj_dra_runs.csv", df_runs)
            
            pubmed_ids = bioproject_data.get("publication", []) if bioproject_data else []
            pubmed_id = pubmed_ids[0] if pubmed_ids else ""
            citation_count = 0
            if pubmed_id:
                citation_count = self._fetch_pubmed_citations(pubmed_id)
            
            summary = project.get("description", "")
            if bioproject_data and bioproject_data.get("description"):
                summary = bioproject_data["description"]
            
            sc_project = ScRNAProject(
                id=project_id,
                sample_id=project.get("secondary_study_accession", ""),
                title=project.get("study_title", ""),
                disease_general="",
                disease="",
                pubmed=pubmed_id,
                source_database="ENA/DDBJ",
                access_link=f"https://www.ebi.ac.uk/ena/browser/view/{project_id}",
                open_status="public" if project.get("first_public") else "restricted",
                ethnicity="",
                sex="",
                tissue="",
                sequencing_platform=platform,
                experiment_design="single-cell RNA-seq",
                sample_type="",
                summary=summary,
                citation_count=citation_count,
                publication_date=project.get("first_public", ""),
                submission_date=bioproject_data.get("submission_date", "") if bioproject_data else "",
                last_update_date=project.get("last_updated", ""),
                contact_name="",
                contact_email="",
                contact_institute=project.get("center_name", ""),
                data_tier="raw",
                tissue_location="",
                supplementary_information=json.dumps({"runs_count": len(runs)}, ensure_ascii=False)
            )
            
            enriched_projects.append(sc_project)
        
        return enriched_projects
    
    def enrich_gea_metadata(self, experiments: List[Dict]) -> List[ScRNAProject]:
        """丰富 GEA 实验元数据"""
        self.logger.info(f"Enriching metadata for {len(experiments)} GEA experiments...")
        enriched_experiments = []
        
        for exp in tqdm(experiments, desc="Enriching GEA metadata"):
            exp_id = exp.get("study_accession", "")
            
            # E-GEAD 实验通常也会有 secondary_study_accession
            secondary_id = exp.get("secondary_study_accession", "")
            
            sc_experiment = ScRNAProject(
                id=exp_id,
                sample_id=secondary_id,
                title=exp.get("study_title", ""),
                disease_general="",
                disease="",
                pubmed="",
                source_database="DDBJ GEA",
                access_link=f"https://ddbj.nig.ac.jp/gea/browse/{exp_id}",
                open_status="public",
                ethnicity="",
                sex="",
                tissue="",
                sequencing_platform="",
                experiment_design="single-cell RNA-seq",
                sample_type="",
                summary=exp.get("description", ""),
                citation_count=0,
                publication_date=exp.get("first_public", ""),
                submission_date="",
                last_update_date=exp.get("last_updated", ""),
                contact_name="",
                contact_email="",
                contact_institute=exp.get("center_name", ""),
                data_tier="processed_matrix",
                tissue_location="",
                supplementary_information=json.dumps(exp.get("study_alias", ""), ensure_ascii=False)
            )
            
            enriched_experiments.append(sc_experiment)
        
        return enriched_experiments
    
    def _append_to_csv(self, filepath: Path, df: pd.DataFrame) -> None:
        """追加数据到 CSV 文件"""
        if filepath.exists():
            df.to_csv(filepath, mode='a', header=False, index=False, encoding='utf-8-sig')
        else:
            df.to_csv(filepath, index=False, encoding='utf-8-sig')
    
    def save_results(self, all_projects: List[ScRNAProject]) -> None:
        """保存最终结果"""
        self.logger.info(f"Saving {len(all_projects)} projects to output files...")
        
        df = pd.DataFrame([asdict(p) for p in all_projects])
        
        csv_path = self.output_dir / "ddbj_human_sc_rna.csv"
        df.to_csv(csv_path, index=False, encoding='utf-8-sig')
        self.logger.info(f"Saved CSV to {csv_path}")
        
        xlsx_path = self.output_dir / "ddbj_human_sc_rna.xlsx"
        df.to_excel(xlsx_path, index=False, engine='openpyxl')
        self.logger.info(f"Saved Excel to {xlsx_path}")
        
        self._generate_statistics(all_projects)
    
    def _generate_statistics(self, projects: List[ScRNAProject]) -> None:
        """生成统计信息"""
        stats_path = self.output_dir / "ddbj_statistics.txt"
        
        total_projects = len(projects)
        databases = defaultdict(int)
        platforms = defaultdict(int)
        years = defaultdict(int)
        tissues = defaultdict(int)
        
        for p in projects:
            databases[p.source_database] += 1
            if p.sequencing_platform:
                platforms[p.sequencing_platform] += 1
            if p.tissue:
                tissues[p.tissue] += 1
            
            if p.submission_date:
                year = p.submission_date[:4]
                if year.isdigit():
                    years[year] += 1
        
        with open(stats_path, 'w', encoding='utf-8') as f:
            f.write(f"DDBJ Human Single-Cell RNA-seq Collection Statistics\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"{'='*60}\n\n")
            
            f.write(f"Total Projects: {total_projects}\n\n")
            
            f.write(f"By Database:\n")
            for db, count in sorted(databases.items(), key=lambda x: x[1], reverse=True):
                f.write(f"  {db}: {count}\n")
            f.write(f"\n")
            
            f.write(f"By Platform (Top 10):\n")
            for platform, count in sorted(platforms.items(), key=lambda x: x[1], reverse=True)[:10]:
                f.write(f"  {platform}: {count}\n")
            f.write(f"\n")
            
            f.write(f"By Tissue (Top 10):\n")
            for tissue, count in sorted(tissues.items(), key=lambda x: x[1], reverse=True)[:10]:
                f.write(f"  {tissue}: {count}\n")
            f.write(f"\n")
            
            f.write(f"By Submission Year:\n")
            for year, count in sorted(years.items()):
                f.write(f"  {year}: {count}\n")
        
        self.logger.info(f"Statistics saved to {stats_path}")
    
    def run(self) -> None:
        """执行完整采集流程"""
        self.logger.info("="*60)
        self.logger.info("Starting DDBJ Human Single-Cell RNA-seq Collection")
        self.logger.info("="*60)
        
        ena_projects = self.collect_ena_projects()
        
        gea_experiments = self.collect_ddbj_gea_experiments()
        
        enriched_ena = self.enrich_project_metadata(ena_projects)
        
        enriched_gea = self.enrich_gea_metadata(gea_experiments)
        
        all_projects = enriched_ena + enriched_gea
        
        self.save_results(all_projects)
        
        self.logger.info("="*60)
        self.logger.info(f"Collection completed! Total projects: {len(all_projects)}")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info("="*60)


def main() -> None:
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Collect DDBJ Human Single-Cell RNA-seq Metadata",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./ddbj_scrna_output",
        help="Output directory path (default: ./ddbj_scrna_output)"
    )
    
    parser.add_argument(
        "--email",
        type=str,
        default="bioinformatics@research.edu",
        help="Contact email for NCBI E-utilities (default: bioinformatics@research.edu)"
    )
    
    args = parser.parse_args()
    
    output_path = Path(args.output_dir)
    
    try:
        collector = DDBJCollector(output_dir=output_path, email=args.email)
        collector.run()
        
    except KeyboardInterrupt:
        print("\n\nCollection interrupted by user. Exiting...")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()