def search_mendeley_data(self, query: str = "single cell RNA sequencing human") -> pd.DataFrame:
    """
    搜索Mendeley Data
    注意: Mendeley Data需要API密钥
    """
    logger.info("开始检索Mendeley Data...")
    
    # Mendeley Data使用Figshare API
    base_url = "https://api.figshare.com/v2/articles"
    
    params = {
        'search_for': query,
        'item_type': 3,  # dataset
        'page_size': 1000
    }
    
    try:
        all_items = []
        page = 1
        
        while True:
            params['page'] = page
            response = requests.get(base_url, params=params, timeout=30)
            response.raise_for_status()
            
            items = response.json()
            if not items:
                break
            
            all_items.extend(items)
            page += 1
            
            if page > 10:  # 限制页数
                break
            
            time.sleep(1)
        
        logger.info(f"从Mendeley Data找到 {len(all_items)} 个数据集")
        
        # 解析数据
        parsed_data = []
        for item in all_items:
            parsed_data.append({
                'id': item.get('id'),
                'title': item.get('title'),
                'doi': item.get('doi'),
                'url': item.get('url'),
                'published_date': item.get('published_date'),
                'description': item.get('description'),
                'source_database': 'Mendeley Data'
            })
        
        df = pd.DataFrame(parsed_data)
        raw_file = self.output_dir / "mendeley_raw_metadata.csv"
        df.to_csv(raw_file, index=False, encoding='utf-8-sig')
        
        return df
        
    except Exception as e:
        logger.error(f"Mendeley Data检索失败: {str(e)}")
        return pd.DataFrame()