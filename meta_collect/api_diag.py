#!/usr/bin/env python3
"""
服务器端网络和API访问完整诊断工具
生成详细的诊断报告，用于分析数据收集问题
"""

import requests
import json
import sys
from datetime import datetime
from pathlib import Path
import time

class ServerDiagnostic:
    def __init__(self):
        self.results = {}
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        })
        
    def print_section(self, title):
        """打印分节标题"""
        print("\n" + "="*70)
        print(f"  {title}")
        print("="*70)
    
    def test_basic_connection(self, url, name, timeout=10):
        """测试基本连接"""
        print(f"\n[{name}]")
        print(f"URL: {url}")
        
        result = {
            'url': url,
            'name': name,
            'timestamp': datetime.now().isoformat(),
            'accessible': False,
            'status_code': None,
            'response_time': None,
            'error': None,
            'headers': None,
            'content_preview': None
        }
        
        try:
            start_time = time.time()
            response = self.session.get(url, timeout=timeout)
            response_time = time.time() - start_time
            
            result['accessible'] = True
            result['status_code'] = response.status_code
            result['response_time'] = round(response_time, 2)
            result['headers'] = dict(response.headers)
            
            # 尝试解析JSON
            try:
                json_data = response.json()
                result['content_type'] = 'json'
                result['content_preview'] = str(json_data)[:500]
                print(f"✓ Status: {response.status_code}")
                print(f"✓ Response Time: {response_time:.2f}s")
                print(f"✓ Content-Type: {response.headers.get('Content-Type', 'Unknown')}")
                print(f"✓ JSON Keys: {list(json_data.keys()) if isinstance(json_data, dict) else 'Array'}")
                
                # 如果是列表，显示长度
                if isinstance(json_data, list):
                    print(f"✓ Array Length: {len(json_data)}")
                    if len(json_data) > 0:
                        print(f"✓ First Item Keys: {list(json_data[0].keys()) if isinstance(json_data[0], dict) else 'Not a dict'}")
                
            except ValueError:
                result['content_type'] = 'html/text'
                result['content_preview'] = response.text[:500]
                print(f"✓ Status: {response.status_code}")
                print(f"✓ Response Time: {response_time:.2f}s")
                print(f"⚠ Content is not JSON")
                print(f"  Preview: {response.text[:200]}...")
            
            return result
            
        except requests.exceptions.ConnectionError as e:
            result['error'] = f"ConnectionError: {str(e)}"
            print(f"✗ Connection Error: {e}")
            return result
            
        except requests.exceptions.Timeout:
            result['error'] = "Timeout"
            print(f"✗ Timeout after {timeout}s")
            return result
            
        except requests.exceptions.HTTPError as e:
            result['error'] = f"HTTPError: {str(e)}"
            result['status_code'] = e.response.status_code if e.response else None
            print(f"✗ HTTP Error: {e}")
            return result
            
        except Exception as e:
            result['error'] = f"Exception: {str(e)}"
            print(f"✗ Unexpected Error: {e}")
            return result
    
    def test_api_endpoint(self, url, name, method='GET', params=None, data=None, headers=None):
        """测试API端点（支持不同HTTP方法）"""
        print(f"\n[{name}]")
        print(f"URL: {url}")
        print(f"Method: {method}")
        if params:
            print(f"Params: {params}")
        
        result = {
            'url': url,
            'name': name,
            'method': method,
            'timestamp': datetime.now().isoformat(),
            'accessible': False,
            'status_code': None,
            'error': None,
            'data_sample': None
        }
        
        try:
            if headers:
                session_headers = self.session.headers.copy()
                session_headers.update(headers)
            else:
                session_headers = self.session.headers
            
            if method == 'GET':
                response = requests.get(url, params=params, headers=session_headers, timeout=30)
            elif method == 'POST':
                response = requests.post(url, json=data, headers=session_headers, timeout=30)
            else:
                print(f"✗ Unsupported method: {method}")
                return result
            
            result['status_code'] = response.status_code
            
            if response.status_code == 200:
                result['accessible'] = True
                try:
                    json_data = response.json()
                    result['data_sample'] = str(json_data)[:500]
                    
                    print(f"✓ Status: 200 OK")
                    print(f"✓ Response Type: JSON")
                    
                    # 分析数据结构
                    if isinstance(json_data, dict):
                        print(f"✓ Response Keys: {list(json_data.keys())}")
                        
                        # 检查常见的分页字段
                        if 'results' in json_data:
                            print(f"✓ Results Count: {len(json_data['results'])}")
                        if 'total' in json_data:
                            print(f"✓ Total Count: {json_data['total']}")
                        if 'hits' in json_data:
                            print(f"✓ Hits: {json_data['hits']}")
                            
                    elif isinstance(json_data, list):
                        print(f"✓ Array Length: {len(json_data)}")
                        if len(json_data) > 0:
                            print(f"✓ First Item: {str(json_data[0])[:200]}...")
                    
                except ValueError:
                    print(f"✓ Status: 200 OK")
                    print(f"⚠ Response is not JSON")
                    print(f"  Content Type: {response.headers.get('Content-Type')}")
                    print(f"  Preview: {response.text[:300]}...")
                    
            else:
                result['error'] = f"HTTP {response.status_code}"
                print(f"✗ Status: {response.status_code}")
                print(f"  Response: {response.text[:300]}...")
                
        except Exception as e:
            result['error'] = str(e)
            print(f"✗ Error: {e}")
        
        return result
    
    def test_cellxgene_detailed(self):
        """详细测试CELLxGENE API"""
        self.print_section("CELLxGENE API 详细测试")
        
        base_url = "https://api.cellxgene.cziscience.com/curation/v1"
        
        # 测试collections端点
        collections_result = self.test_api_endpoint(
            f"{base_url}/collections",
            "CELLxGENE Collections"
        )
        
        # 测试datasets端点
        datasets_result = self.test_api_endpoint(
            f"{base_url}/datasets",
            "CELLxGENE Datasets"
        )
        
        self.results['cellxgene'] = {
            'collections': collections_result,
            'datasets': datasets_result
        }
    
    def test_broad_scp_detailed(self):
        """详细测试Broad SCP"""
        self.print_section("Broad Single Cell Portal 详细测试")
        
        # 测试主页
        main_result = self.test_basic_connection(
            "https://singlecell.broadinstitute.org",
            "Broad SCP Main Page"
        )
        
        # 测试API v1
        api_v1_result = self.test_api_endpoint(
            "https://singlecell.broadinstitute.org/single_cell/api/v1/search/studies",
            "Broad SCP API v1 Studies",
            params={'page': 1, 'per_page': 10}
        )
        
        # 尝试不带参数
        api_v1_simple = self.test_api_endpoint(
            "https://singlecell.broadinstitute.org/single_cell/api/v1/studies",
            "Broad SCP API v1 Studies (Simple)"
        )
        
        self.results['broad_scp'] = {
            'main': main_result,
            'api_v1_search': api_v1_result,
            'api_v1_simple': api_v1_simple
        }
    
    def test_zenodo_detailed(self):
        """详细测试Zenodo"""
        self.print_section("Zenodo 详细测试")
        
        # 测试主页
        main_result = self.test_basic_connection(
            "https://zenodo.org",
            "Zenodo Main Page"
        )
        
        # 测试API - 简单查询
        api_simple = self.test_api_endpoint(
            "https://zenodo.org/api/records",
            "Zenodo API - Simple Query",
            params={'q': 'test', 'size': 5}
        )
        
        # 测试API - scRNA-seq查询
        api_scrna = self.test_api_endpoint(
            "https://zenodo.org/api/records",
            "Zenodo API - scRNA-seq Query",
            params={'q': 'single cell RNA-seq human', 'size': 10, 'sort': 'mostrecent'}
        )
        
        # 测试OAI-PMH接口
        oai_pmh = self.test_basic_connection(
            "https://zenodo.org/oai2d?verb=Identify",
            "Zenodo OAI-PMH Interface"
        )
        
        self.results['zenodo'] = {
            'main': main_result,
            'api_simple': api_simple,
            'api_scrna': api_scrna,
            'oai_pmh': oai_pmh
        }
    
    def test_figshare_detailed(self):
        """详细测试Figshare"""
        self.print_section("Figshare 详细测试")
        
        # 测试主页
        main_result = self.test_basic_connection(
            "https://figshare.com",
            "Figshare Main Page"
        )
        
        # 测试API v2 - articles
        api_articles = self.test_api_endpoint(
            "https://api.figshare.com/v2/articles",
            "Figshare API v2 - Articles",
            params={'page_size': 10}
        )
        
        # 测试API v2 - search
        api_search = self.test_api_endpoint(
            "https://api.figshare.com/v2/articles/search",
            "Figshare API v2 - Search",
            method='POST',
            data={'search_for': 'single cell RNA-seq', 'page_size': 10}
        )
        
        self.results['figshare'] = {
            'main': main_result,
            'api_articles': api_articles,
            'api_search': api_search
        }
    
    def test_dryad_detailed(self):
        """详细测试DRYAD"""
        self.print_section("DRYAD 详细测试")
        
        # 测试主页
        main_result = self.test_basic_connection(
            "https://datadryad.org",
            "DRYAD Main Page"
        )
        
        # 测试API v2
        api_v2 = self.test_api_endpoint(
            "https://datadryad.org/api/v2/datasets",
            "DRYAD API v2 - Datasets",
            params={'per_page': 10}
        )
        
        # 测试搜索
        api_search = self.test_api_endpoint(
            "https://datadryad.org/api/v2/search",
            "DRYAD API v2 - Search",
            params={'q': 'single cell RNA-seq', 'per_page': 10}
        )
        
        self.results['dryad'] = {
            'main': main_result,
            'api_v2': api_v2,
            'api_search': api_search
        }
    
    def test_geo_detailed(self):
        """详细测试GEO"""
        self.print_section("GEO (NCBI) 详细测试")
        
        # 测试E-utilities API
        eutils_search = self.test_api_endpoint(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            "GEO E-utilities Search",
            params={
                'db': 'gds',
                'term': 'single cell RNA-seq human[Organism]',
                'retmode': 'json',
                'retmax': 10
            }
        )
        
        self.results['geo'] = {
            'eutils_search': eutils_search
        }
    
    def generate_report(self):
        """生成诊断报告"""
        self.print_section("诊断报告摘要")
        
        report = {
            'timestamp': datetime.now().isoformat(),
            'summary': {},
            'recommendations': [],
            'detailed_results': self.results
        }
        
        # 分析每个数据库的状态
        for db_name, db_results in self.results.items():
            accessible_count = 0
            total_tests = 0
            
            for test_name, test_result in db_results.items():
                total_tests += 1
                if test_result.get('accessible') or test_result.get('status_code') == 200:
                    accessible_count += 1
            
            status = "✓ OK" if accessible_count > 0 else "✗ FAILED"
            report['summary'][db_name] = {
                'status': status,
                'accessible_tests': accessible_count,
                'total_tests': total_tests
            }
            
            print(f"\n{db_name.upper()}:")
            print(f"  Status: {status}")
            print(f"  Successful Tests: {accessible_count}/{total_tests}")
            
            # 详细信息
            for test_name, test_result in db_results.items():
                test_status = "✓" if (test_result.get('accessible') or test_result.get('status_code') == 200) else "✗"
                print(f"    {test_status} {test_name}: {test_result.get('status_code', 'N/A')}")
                if test_result.get('error'):
                    print(f"       Error: {test_result['error']}")
        
        # 生成建议
        self.print_section("建议和解决方案")
        
        recommendations = []
        
        # CELLxGENE
        if not self.results['cellxgene']['collections'].get('accessible'):
            recommendations.append({
                'database': 'CELLxGENE',
                'issue': 'API不可访问',
                'solutions': [
                    '检查网络连接',
                    '确认API地址是否变更',
                    '等待片刻后重试'
                ]
            })
        
        # Broad SCP
        broad_ok = any(r.get('accessible') or r.get('status_code') == 200 
                       for r in self.results['broad_scp'].values())
        if not broad_ok:
            recommendations.append({
                'database': 'Broad SCP',
                'issue': '所有API端点均无法访问',
                'solutions': [
                    '服务可能暂时不可用，稍后重试',
                    '考虑使用网页抓取方式',
                    '联系Broad Institute获取API访问权限',
                    '查看 https://github.com/broadinstitute/single_cell_portal/wiki/API-Documentation'
                ]
            })
        
        # Zenodo
        zenodo_ok = any(r.get('accessible') or r.get('status_code') == 200 
                        for r in self.results['zenodo'].values())
        if not zenodo_ok:
            recommendations.append({
                'database': 'Zenodo',
                'issue': '连接被拒绝',
                'solutions': [
                    '检查防火墙设置',
                    '配置HTTP代理（如果在受限网络环境）',
                    '尝试使用VPN',
                    '检查/etc/hosts文件是否屏蔽了zenodo.org',
                    '尝试使用OAI-PMH协议',
                    '考虑申请API token: https://zenodo.org/account/settings/applications/'
                ]
            })
        
        # Figshare
        figshare_ok = any(r.get('accessible') or r.get('status_code') == 200 
                          for r in self.results['figshare'].values())
        if not figshare_ok:
            fig_status = self.results['figshare']['api_articles'].get('status_code')
            if fig_status == 403:
                recommendations.append({
                    'database': 'Figshare',
                    'issue': 'API访问被禁止 (403)',
                    'solutions': [
                        '申请Figshare API token: https://figshare.com/account/applications',
                        '在请求头中添加token: Authorization: token YOUR_TOKEN',
                        '使用公开搜索API而非直接访问articles',
                        '降低请求频率，避免触发限流'
                    ]
                })
            else:
                recommendations.append({
                    'database': 'Figshare',
                    'issue': 'API无法访问',
                    'solutions': [
                        '检查网络连接',
                        '申请API token',
                        '尝试使用不同的endpoint'
                    ]
                })
        
        # 打印建议
        for rec in recommendations:
            print(f"\n【{rec['database']}】")
            print(f"问题: {rec['issue']}")
            print("解决方案:")
            for i, solution in enumerate(rec['solutions'], 1):
                print(f"  {i}. {solution}")
        
        report['recommendations'] = recommendations
        
        # 保存报告
        output_dir = Path("diagnostic_reports")
        output_dir.mkdir(exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = output_dir / f"diagnostic_report_{timestamp}.json"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        print(f"\n{'='*70}")
        print(f"完整报告已保存到: {report_file}")
        print(f"{'='*70}")
        
        return report
    
    def run_all_tests(self):
        """运行所有测试"""
        print("="*70)
        print("  服务器端数据库API诊断工具")
        print("  开始时间:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print("="*70)
        
        self.test_cellxgene_detailed()
        self.test_broad_scp_detailed()
        self.test_zenodo_detailed()
        self.test_figshare_detailed()
        self.test_dryad_detailed()
        self.test_geo_detailed()
        
        report = self.generate_report()
        
        return report


def main():
    """主函数"""
    diagnostic = ServerDiagnostic()
    
    try:
        report = diagnostic.run_all_tests()
        
        print("\n" + "="*70)
        print("诊断完成!")
        print("="*70)
        print("\n请将生成的JSON报告文件发送给我，我将根据报告调整代码。")
        print("报告位置: diagnostic_reports/diagnostic_report_*.json")
        
        return 0
        
    except KeyboardInterrupt:
        print("\n\n诊断被用户中断")
        return 1
        
    except Exception as e:
        print(f"\n\n发生错误: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())