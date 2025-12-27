"""
Figshare API 诊断工具
精确判断403错误的真实原因
"""

import requests
import json
from datetime import datetime

def test_figshare_api():
    """全面诊断Figshare API"""
    
    print("=" * 70)
    print("Figshare API 诊断工具")
    print("=" * 70)
    print(f"测试时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"测试主机: node3")
    print("=" * 70)
    
    # 测试1: 官方文档的API endpoint
    print("\n【测试 1/5】官方API endpoint")
    print("-" * 70)
    
    endpoints = [
        "https://api.figshare.com/v2/articles",
        "https://api.figshare.com/v2/account/articles",  # 需要认证
        "https://stats.figshare.com/",  # 统计API
    ]
    
    for url in endpoints:
        try:
            response = requests.get(url, timeout=10)
            print(f"URL: {url}")
            print(f"状态码: {response.status_code}")
            print(f"响应头: {dict(list(response.headers.items())[:5])}")
            
            if response.status_code == 403:
                print(f"403详情: {response.text[:200]}")
            elif response.status_code == 200:
                print(f"成功! 返回数据示例: {response.text[:100]}")
            
            print()
        except Exception as e:
            print(f"请求失败: {str(e)}\n")
    
    # 测试2: 不同的User-Agent
    print("\n【测试 2/5】不同User-Agent")
    print("-" * 70)
    
    user_agents = [
        "python-requests/2.31.0",  # 默认
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",  # 浏览器
        "curl/7.68.0",  # curl
        "",  # 空
    ]
    
    for ua in user_agents:
        try:
            headers = {'User-Agent': ua} if ua else {}
            response = requests.get(
                "https://api.figshare.com/v2/articles",
                headers=headers,
                timeout=10
            )
            print(f"User-Agent: {ua or '(空)'}")
            print(f"状态码: {response.status_code}\n")
        except Exception as e:
            print(f"失败: {str(e)}\n")
    
    # 测试3: 带参数请求
    print("\n【测试 3/5】不同请求参数")
    print("-" * 70)
    
    params_list = [
        {},
        {'page': 1},
        {'page': 1, 'page_size': 10},
        {'search_for': ':title: RNA'},
    ]
    
    for params in params_list:
        try:
            response = requests.get(
                "https://api.figshare.com/v2/articles",
                params=params,
                timeout=10
            )
            print(f"参数: {params}")
            print(f"状态码: {response.status_code}")
            print(f"完整URL: {response.url}\n")
        except Exception as e:
            print(f"失败: {str(e)}\n")
    
    # 测试4: 检查API是否需要认证
    print("\n【测试 4/5】API文档检查")
    print("-" * 70)
    
    try:
        # 尝试访问公开文档
        doc_url = "https://docs.figshare.com/"
        response = requests.get(doc_url, timeout=10)
        print(f"文档访问: {response.status_code}")
        
        # 检查是否有robots.txt或API政策说明
        robots_url = "https://api.figshare.com/robots.txt"
        response = requests.get(robots_url, timeout=10)
        print(f"robots.txt: {response.status_code}")
        if response.status_code == 200:
            print(f"内容:\n{response.text[:300]}")
    except Exception as e:
        print(f"文档检查失败: {str(e)}")
    
    # 测试5: 使用curl命令（系统级测试）
    print("\n【测试 5/5】系统级curl测试")
    print("-" * 70)
    print("请在命令行执行以下命令:")
    print()
    print("curl -v https://api.figshare.com/v2/articles")
    print()
    print("curl -H 'User-Agent: Mozilla/5.0' https://api.figshare.com/v2/articles")
    print()
    
    # 尝试系统调用curl
    import subprocess
    try:
        result = subprocess.run(
            ['curl', '-s', '-w', '\\nHTTP_CODE:%{http_code}', 
             'https://api.figshare.com/v2/articles'],
            capture_output=True,
            text=True,
            timeout=10
        )
        print("curl 输出:")
        print(result.stdout[-500:])
        print(result.stderr[-200:] if result.stderr else "")
    except Exception as e:
        print(f"无法执行curl: {str(e)}")
    
    print("\n" + "=" * 70)
    print("诊断完成")
    print("=" * 70)


if __name__ == "__main__":
    test_figshare_api()