"""
KPMP网站智能爬虫 - 基于页面结构自适应爬取
Author: 针对KPMP Atlas Explorer和Repository设计
Date: 2025-12-04
"""

import json
import time
import logging
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from selenium.webdriver.remote.webelement import WebElement

import pandas as pd

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class KPMPWebScraper:
    """KPMP网站智能爬虫"""
    
    def __init__(self, output_dir: str = "./kpmp_web_scraping"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.driver = None
        self.wait = None
        
    def setup_driver(self, headless: bool = True):
        """配置浏览器"""
        chrome_options = Options()
        if headless:
            chrome_options.add_argument('--headless')
        chrome_options.add_argument('--no-sandbox')
        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument('--disable-gpu')
        chrome_options.add_argument('--window-size=1920,1080')
        
        # 启用网络日志
        chrome_options.set_capability('goog:loggingPrefs', {'performance': 'ALL'})
        
        self.driver = webdriver.Chrome(options=chrome_options)
        self.wait = WebDriverWait(self.driver, 20)
        logger.info("浏览器初始化完成")
        
    def wait_for_page_load(self, timeout: int = 15):
        """等待页面完全加载"""
        try:
            # 等待页面加载完成
            self.wait.until(
                lambda d: d.execute_script('return document.readyState') == 'complete'
            )
            time.sleep(3)  # 额外等待JavaScript渲染
            logger.info("页面加载完成")
            return True
        except TimeoutException:
            logger.warning("页面加载超时")
            return False
    
    def extract_network_requests(self) -> List[Dict]:
        """提取所有网络请求"""
        logs = self.driver.get_log('performance')
        requests = []
        
        for entry in logs:
            try:
                log = json.loads(entry['message'])['message']
                
                if log['method'] == 'Network.requestWillBeSent':
                    request = log['params']['request']
                    url = request.get('url', '')
                    
                    # 只保留API相关请求
                    if any(kw in url.lower() for kw in ['api', 'data', 'json', 'graphql', 'query']):
                        requests.append({
                            'url': url,
                            'method': request.get('method', 'GET'),
                            'headers': request.get('headers', {}),
                            'timestamp': entry['timestamp']
                        })
            except Exception as e:
                pass
        
        logger.info(f"捕获到 {len(requests)} 个API请求")
        return requests
    
    def find_all_links(self) -> List[Dict[str, str]]:
        """查找页面所有链接"""
        links = []
        try:
            elements = self.driver.find_elements(By.TAG_NAME, 'a')
            for elem in elements:
                href = elem.get_attribute('href')
                text = elem.text.strip()
                if href and href.startswith('http'):
                    links.append({
                        'url': href,
                        'text': text,
                        'title': elem.get_attribute('title') or ''
                    })
        except Exception as e:
            logger.error(f"提取链接失败: {str(e)}")
        
        logger.info(f"找到 {len(links)} 个链接")
        return links
    
    def extract_data_tables(self) -> List[pd.DataFrame]:
        """提取所有表格数据"""
        tables = []
        try:
            table_elements = self.driver.find_elements(By.TAG_NAME, 'table')
            logger.info(f"找到 {len(table_elements)} 个表格")
            
            for idx, table in enumerate(table_elements):
                try:
                    # 提取表头
                    headers = []
                    header_elements = table.find_elements(By.TAG_NAME, 'th')
                    if header_elements:
                        headers = [h.text.strip() for h in header_elements]
                    
                    # 提取行数据
                    rows_data = []
                    row_elements = table.find_elements(By.TAG_NAME, 'tr')
                    
                    for row in row_elements:
                        cells = row.find_elements(By.TAG_NAME, 'td')
                        if cells:
                            row_data = [cell.text.strip() for cell in cells]
                            rows_data.append(row_data)
                    
                    if rows_data:
                        if headers and len(headers) == len(rows_data[0]):
                            df = pd.DataFrame(rows_data, columns=headers)
                        else:
                            df = pd.DataFrame(rows_data)
                        
                        tables.append(df)
                        logger.info(f"表格 {idx + 1}: {df.shape[0]} 行 x {df.shape[1]} 列")
                
                except Exception as e:
                    logger.warning(f"提取表格 {idx + 1} 失败: {str(e)}")
        
        except Exception as e:
            logger.error(f"查找表格失败: {str(e)}")
        
        return tables
    
    def extract_dataset_cards(self) -> List[Dict]:
        """提取数据集卡片信息（常见于现代Web UI）"""
        datasets = []
        
        # 尝试多种可能的选择器
        selectors = [
            "[class*='dataset']",
            "[class*='card']",
            "[class*='item']",
            "[class*='study']",
            "[class*='sample']",
            "[data-testid*='dataset']",
            "[data-testid*='card']",
            "article",
            ".MuiCard-root",  # Material-UI
            ".card",  # Bootstrap
        ]
        
        for selector in selectors:
            try:
                elements = self.driver.find_elements(By.CSS_SELECTOR, selector)
                
                if len(elements) > 5:  # 找到了合理数量的元素
                    logger.info(f"使用选择器 '{selector}' 找到 {len(elements)} 个元素")
                    
                    for idx, elem in enumerate(elements[:100]):  # 限制数量
                        try:
                            dataset_info = {
                                'index': idx,
                                'selector': selector,
                                'text': elem.text[:1000],  # 限制长度
                                'html': elem.get_attribute('outerHTML')[:2000],
                                'class': elem.get_attribute('class'),
                                'id': elem.get_attribute('id'),
                            }
                            
                            # 尝试提取具体字段
                            for tag in ['h1', 'h2', 'h3', 'h4', 'h5', 'h6']:
                                try:
                                    heading = elem.find_element(By.TAG_NAME, tag)
                                    dataset_info['title'] = heading.text.strip()
                                    break
                                except:
                                    pass
                            
                            # 查找链接
                            try:
                                link = elem.find_element(By.TAG_NAME, 'a')
                                dataset_info['link'] = link.get_attribute('href')
                            except:
                                pass
                            
                            datasets.append(dataset_info)
                        
                        except Exception as e:
                            logger.debug(f"提取元素 {idx} 失败: {str(e)}")
                    
                    break  # 找到一个有效选择器就停止
            
            except Exception as e:
                logger.debug(f"选择器 '{selector}' 无效: {str(e)}")
        
        logger.info(f"共提取 {len(datasets)} 个数据集卡片")
        return datasets
    
    def extract_page_structure(self) -> Dict:
        """分析页面结构"""
        structure = {
            'title': self.driver.title,
            'url': self.driver.current_url,
            'timestamp': datetime.now().isoformat(),
        }
        
        # 统计各种元素
        try:
            structure['element_counts'] = {
                'tables': len(self.driver.find_elements(By.TAG_NAME, 'table')),
                'forms': len(self.driver.find_elements(By.TAG_NAME, 'form')),
                'buttons': len(self.driver.find_elements(By.TAG_NAME, 'button')),
                'inputs': len(self.driver.find_elements(By.TAG_NAME, 'input')),
                'links': len(self.driver.find_elements(By.TAG_NAME, 'a')),
                'divs': len(self.driver.find_elements(By.TAG_NAME, 'div')),
                'articles': len(self.driver.find_elements(By.TAG_NAME, 'article')),
            }
        except Exception as e:
            logger.error(f"统计元素失败: {str(e)}")
        
        # 检测常见框架
        try:
            scripts = self.driver.find_elements(By.TAG_NAME, 'script')
            script_content = ' '.join([s.get_attribute('src') or '' for s in scripts])
            
            structure['frameworks_detected'] = {
                'react': 'react' in script_content.lower(),
                'vue': 'vue' in script_content.lower(),
                'angular': 'angular' in script_content.lower(),
                'bootstrap': 'bootstrap' in script_content.lower(),
                'material-ui': 'material' in script_content.lower(),
            }
        except Exception as e:
            logger.error(f"检测框架失败: {str(e)}")
        
        return structure
    
    def scroll_and_load_more(self, scroll_times: int = 5):
        """滚动页面加载更多内容"""
        logger.info("开始滚动加载...")
        
        for i in range(scroll_times):
            # 滚动到底部
            self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
            time.sleep(2)
            
            # 查找"加载更多"按钮
            try:
                load_more_selectors = [
                    "//button[contains(text(), 'Load More')]",
                    "//button[contains(text(), 'Show More')]",
                    "//button[contains(text(), '加载更多')]",
                    "//button[contains(text(), '更多')]",
                    "//a[contains(text(), 'Load More')]",
                ]
                
                for selector in load_more_selectors:
                    try:
                        button = self.driver.find_element(By.XPATH, selector)
                        button.click()
                        logger.info(f"点击了'加载更多'按钮 (第{i+1}次)")
                        time.sleep(2)
                        break
                    except:
                        pass
            except:
                pass
            
            logger.info(f"完成第 {i+1}/{scroll_times} 次滚动")
    
    def scrape_explorer_page(self) -> Dict:
        """爬取Explorer页面"""
        url = "https://atlas.kpmp.org/explorer"
        logger.info(f"开始爬取: {url}")
        
        self.driver.get(url)
        self.wait_for_page_load()
        
        # 滚动加载更多内容
        self.scroll_and_load_more(scroll_times=3)
        
        # 提取各类数据
        data = {
            'url': url,
            'timestamp': datetime.now().isoformat(),
            'page_structure': self.extract_page_structure(),
            'network_requests': self.extract_network_requests(),
            'links': self.find_all_links(),
            'tables': [],
            'dataset_cards': self.extract_dataset_cards(),
        }
        
        # 提取表格
        tables = self.extract_data_tables()
        for idx, df in enumerate(tables):
            table_file = self.output_dir / f"explorer_table_{idx}.csv"
            df.to_csv(table_file, index=False)
            data['tables'].append({
                'index': idx,
                'shape': df.shape,
                'columns': df.columns.tolist(),
                'file': str(table_file)
            })
        
        # 尝试提取特定内容
        data['specific_elements'] = self._extract_explorer_specific_elements()
        
        return data
    
    def _extract_explorer_specific_elements(self) -> Dict:
        """提取Explorer页面特定元素"""
        specific = {}
        
        # 尝试查找数据集统计信息
        try:
            stat_selectors = [
                "[class*='stat']",
                "[class*='count']",
                "[class*='summary']",
                "[class*='metric']"
            ]
            
            for selector in stat_selectors:
                elements = self.driver.find_elements(By.CSS_SELECTOR, selector)
                if elements:
                    specific['statistics'] = [
                        {'text': e.text, 'class': e.get_attribute('class')}
                        for e in elements[:20]
                    ]
                    break
        except Exception as e:
            logger.debug(f"提取统计信息失败: {str(e)}")
        
        # 尝试查找过滤器/搜索框
        try:
            filters = self.driver.find_elements(By.CSS_SELECTOR, 
                "[class*='filter'], [class*='search'], input[type='search']")
            if filters:
                specific['filters'] = [
                    {
                        'type': f.get_attribute('type'),
                        'name': f.get_attribute('name'),
                        'placeholder': f.get_attribute('placeholder')
                    }
                    for f in filters
                ]
        except Exception as e:
            logger.debug(f"提取过滤器失败: {str(e)}")
        
        return specific
    
    def scrape_repository_page(self) -> Dict:
        """爬取Repository页面"""
        url = "https://atlas.kpmp.org/repository"
        logger.info(f"开始爬取: {url}")
        
        self.driver.get(url)
        self.wait_for_page_load()
        
        # 滚动加载
        self.scroll_and_load_more(scroll_times=3)
        
        # 提取数据
        data = {
            'url': url,
            'timestamp': datetime.now().isoformat(),
            'page_structure': self.extract_page_structure(),
            'network_requests': self.extract_network_requests(),
            'links': self.find_all_links(),
            'tables': [],
            'dataset_cards': self.extract_dataset_cards(),
        }
        
        # 提取表格
        tables = self.extract_data_tables()
        for idx, df in enumerate(tables):
            table_file = self.output_dir / f"repository_table_{idx}.csv"
            df.to_csv(table_file, index=False)
            data['tables'].append({
                'index': idx,
                'shape': df.shape,
                'columns': df.columns.tolist(),
                'file': str(table_file)
            })
        
        # Repository特定元素
        data['specific_elements'] = self._extract_repository_specific_elements()
        
        return data
    
    def _extract_repository_specific_elements(self) -> Dict:
        """提取Repository页面特定元素"""
        specific = {}
        
        # 尝试查找文件列表
        try:
            file_selectors = [
                "[class*='file']",
                "[class*='download']",
                "a[href*='.gz']",
                "a[href*='.zip']",
                "a[href*='.tar']",
            ]
            
            files = []
            for selector in file_selectors:
                elements = self.driver.find_elements(By.CSS_SELECTOR, selector)
                if elements:
                    files.extend([
                        {
                            'text': e.text,
                            'href': e.get_attribute('href'),
                            'class': e.get_attribute('class')
                        }
                        for e in elements[:50]
                    ])
            
            if files:
                specific['files'] = files
                logger.info(f"找到 {len(files)} 个文件链接")
        
        except Exception as e:
            logger.debug(f"提取文件列表失败: {str(e)}")
        
        return specific
    
    def compare_with_geo_data(self, geo_csv_path: str) -> Dict:
        """对比KPMP网站数据与GEO数据"""
        logger.info("对比KPMP网站数据与GEO数据...")
        
        try:
            # 读取GEO数据
            geo_df = pd.read_csv(geo_csv_path)
            logger.info(f"GEO数据: {len(geo_df)} 条记录")
            
            # 提取GEO中的关键信息
            geo_ids = set(geo_df['id'].astype(str))
            geo_titles = set(geo_df['title'].str.lower())
            
            # 读取爬取的数据
            explorer_file = self.output_dir / "explorer_data.json"
            repo_file = self.output_dir / "repository_data.json"
            
            kpmp_ids = set()
            kpmp_titles = set()
            
            for file in [explorer_file, repo_file]:
                if file.exists():
                    with open(file) as f:
                        data = json.load(f)
                        
                        # 从dataset_cards中提取信息
                        for card in data.get('dataset_cards', []):
                            text = card.get('text', '').lower()
                            
                            # 查找GSE编号
                            import re
                            gse_matches = re.findall(r'gse\d+', text)
                            kpmp_ids.update(gse_matches)
                            
                            # 提取标题
                            title = card.get('title', '').lower()
                            if title:
                                kpmp_titles.add(title)
            
            # 对比分析
            comparison = {
                'geo_total': len(geo_df),
                'kpmp_gse_found': len(kpmp_ids),
                'kpmp_unique_titles': len(kpmp_titles),
                'overlap_ids': list(geo_ids & kpmp_ids),
                'geo_only_ids': list(geo_ids - kpmp_ids)[:20],  # 限制数量
                'kpmp_only_ids': list(kpmp_ids - geo_ids),
                'conclusion': self._generate_comparison_conclusion(
                    len(geo_ids), len(kpmp_ids), len(geo_ids & kpmp_ids)
                )
            }
            
            return comparison
        
        except Exception as e:
            logger.error(f"对比失败: {str(e)}")
            return {'error': str(e)}
    
    def _generate_comparison_conclusion(self, geo_count: int, 
                                       kpmp_count: int, overlap: int) -> str:
        """生成对比结论"""
        if kpmp_count == 0:
            return "❌ KPMP网站未找到任何GSE编号，可能需要人工分析页面内容"
        
        overlap_rate = overlap / geo_count if geo_count > 0 else 0
        
        if overlap_rate > 0.9:
            return f"✅ 高度重叠 ({overlap_rate:.1%})：KPMP网站数据与GEO基本一致，爬取必要性低"
        elif overlap_rate > 0.5:
            return f"⚠️ 部分重叠 ({overlap_rate:.1%})：KPMP可能有额外信息，建议爬取补充"
        else:
            return f"✅ 独立数据源 ({overlap_rate:.1%})：KPMP网站包含大量独家数据，强烈建议爬取"
    
    def run_full_scraping(self, geo_csv_path: str = None):
        """运行完整爬取流程"""
        logger.info("="*60)
        logger.info("启动KPMP网站完整爬取流程")
        logger.info("="*60)
        
        try:
            self.setup_driver(headless=False)  # 调试时使用False查看浏览器
            
            # 1. 爬取Explorer
            logger.info("\n[1/3] 爬取Explorer页面...")
            explorer_data = self.scrape_explorer_page()
            
            explorer_file = self.output_dir / "explorer_data.json"
            with open(explorer_file, 'w', encoding='utf-8') as f:
                json.dump(explorer_data, f, indent=2, ensure_ascii=False, default=str)
            logger.info(f"✓ Explorer数据已保存: {explorer_file}")
            
            # 2. 爬取Repository
            logger.info("\n[2/3] 爬取Repository页面...")
            repo_data = self.scrape_repository_page()
            
            repo_file = self.output_dir / "repository_data.json"
            with open(repo_file, 'w', encoding='utf-8') as f:
                json.dump(repo_data, f, indent=2, ensure_ascii=False, default=str)
            logger.info(f"✓ Repository数据已保存: {repo_file}")
            
            # 3. 对比GEO数据（如果提供）
            comparison = None
            if geo_csv_path and Path(geo_csv_path).exists():
                logger.info("\n[3/3] 对比GEO数据...")
                comparison = self.compare_with_geo_data(geo_csv_path)
                
                comparison_file = self.output_dir / "geo_comparison.json"
                with open(comparison_file, 'w', encoding='utf-8') as f:
                    json.dump(comparison, f, indent=2, ensure_ascii=False)
                logger.info(f"✓ 对比结果已保存: {comparison_file}")
                
                print("\n" + "="*60)
                print("对比结论:")
                print(comparison.get('conclusion', '未知'))
                print("="*60)
            
            # 生成总结报告
            report = self._generate_summary_report(explorer_data, repo_data, comparison)
            
            report_file = self.output_dir / "scraping_report.txt"
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write(report)
            
            print("\n" + report)
            logger.info(f"\n✓ 爬取完成！输出目录: {self.output_dir}")
            
            return {
                'explorer_data': explorer_data,
                'repository_data': repo_data,
                'comparison': comparison,
                'output_dir': str(self.output_dir)
            }
        
        finally:
            if self.driver:
                self.driver.quit()
                logger.info("浏览器已关闭")
    
    def _generate_summary_report(self, explorer_data: Dict, 
                                 repo_data: Dict, comparison: Dict = None) -> str:
        """生成总结报告"""
        report_lines = [
            "="*60,
            "KPMP网站爬取总结报告",
            "="*60,
            f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "【Explorer页面】",
            f"  - 数据集卡片: {len(explorer_data.get('dataset_cards', []))} 个",
            f"  - 表格数量: {len(explorer_data.get('tables', []))} 个",
            f"  - 链接数量: {len(explorer_data.get('links', []))} 个",
            f"  - API请求: {len(explorer_data.get('network_requests', []))} 个",
            "",
            "【Repository页面】",
            f"  - 数据集卡片: {len(repo_data.get('dataset_cards', []))} 个",
            f"  - 表格数量: {len(repo_data.get('tables', []))} 个",
            f"  - 链接数量: {len(repo_data.get('links', []))} 个",
            f"  - API请求: {len(repo_data.get('network_requests', []))} 个",
        ]
        
        if comparison:
            report_lines.extend([
                "",
                "【与GEO数据对比】",
                f"  - GEO总记录数: {comparison.get('geo_total', 0)}",
                f"  - KPMP找到的GSE: {comparison.get('kpmp_gse_found', 0)}",
                f"  - 重叠记录: {len(comparison.get('overlap_ids', []))}",
                f"  - 结论: {comparison.get('conclusion', '未知')}",
            ])
        
        # API端点建议
        all_requests = (explorer_data.get('network_requests', []) + 
                       repo_data.get('network_requests', []))
        unique_apis = set(r['url'] for r in all_requests)
        
        if unique_apis:
            report_lines.extend([
                "",
                "【发现的API端点】",
                f"  共 {len(unique_apis)} 个唯一端点:",
            ])
            for api in list(unique_apis)[:10]:  # 只显示前10个
                report_lines.append(f"    • {api}")
        
        report_lines.append("="*60)
        
        return "\n".join(report_lines)


# ============== 使用示例 ==============

if __name__ == "__main__":
    # 创建爬虫实例
    scraper = KPMPWebScraper(output_dir="./kpmp_web_data")
    
    # 运行完整爬取（可选：提供GEO CSV路径进行对比）
    results = scraper.run_full_scraping(
        geo_csv_path="./kpmp_metadata_collection/kpmp_series_metadata.csv"  # 您之前的GEO数据
    )
    
    print("\n" + "="*60)
    print("爬取任务完成！")
    print(f"输出目录: {results['output_dir']}")
    print("="*60)