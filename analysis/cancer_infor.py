#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1. 读入 complete_human_series_20251125_034124.csv
2. 完全由 AI 判断是否为实体瘤
3. 完全由 AI 智能给出癌种标签：
   - Cancer_general：癌症大类缩写
   - Cancer：癌种英文，若为罕见癌请 AI 自动附加 '-稀有癌种'
结果保存在同目录 complete_human_series_20251125_034124_zero_shot.csv
"""
import os
import json
import time
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

API_KEYS = [
    "sk-UL5YodR7ZL4S9dytpfMWJgmPTXJjkeNSd7Ktq9bbEhElzDfX",
    "sk-WL1hTKhW3sYKuJjBSc1k8wxL0r8ZLcuM7YiYxFjGgVHqAXhU",
    "sk-RAkA28HIT5tiMEfKtXAgbZ9nZKweq5Bnw0WbSwwBdNX7nbi1"
]
API_URL = "https://api.moonshot.cn/v1/chat/completions"
MODEL   = "kimi-k2-thinking"
WORKERS = len(API_KEYS)
RATE_PER_KEY = 10
COOLDOWN     = 1 / RATE_PER_KEY

IN_FILE    = "/mnt/public7/pancancercol/meta_analysis/geo_metadata/comprehensive_metadata/complete_human_series_20251125_034124.csv"
CACHE_FILE = IN_FILE.replace(".csv", "_zero_shot_cache.json")
OUT_FILE   = IN_FILE.replace(".csv", "_zero_shot.csv")

# ---------- 实体瘤 AI 判断 ---------- #
def is_solid_tumor(summary: str) -> bool:
    if pd.isna(summary):
        return False
    prompt = (
        "下面是一段 GEO 数据集的摘要，请判断该样本是否为“实体瘤”（Solid Tumor）。"
        "只回答“是”或“否”，不要解释。\n\n" + str(summary)
    )
    payload = {"model": MODEL, "messages": [{"role": "user", "content": prompt}], "temperature": 0}
    try:
        resp = requests.post(API_URL,
                             headers={"Authorization": f"Bearer {API_KEYS[0]}"},
                             json=payload,
                             timeout=90)
        resp.raise_for_status()
        return resp.json()["choices"][0]["message"]["content"].strip() == "是"
    except Exception as e:
        print(f"[WARN] 实体瘤判断API失败: {e}")
        return False

# ---------- 零先验 AI 癌种标签 ---------- #
def kimi_cancer_tags(summary: str, api_key: str, sid: str) -> tuple:
    if pd.isna(summary) or str(summary).strip() == "":
        return "", ""

    system_prompt = (
        "你是肿瘤生物信息学专家。请根据 GEO 数据集摘要，智能给出两个标签：\n"
        "1. Cancer_general：常见癌症大类英文缩写（如 BC、LC、CRC、GC、HCC、PDAC、PRAD、OV、GBM、SKCM 等）。"
        "若为非肿瘤或无法判断，写 None。\n"
        "2. Cancer：具体癌种英文名称（优先 TCGA 项目的38个标准分类，如 BRCA、LUAD、COAD 等）；"
        "若该癌种不在 TCGA-38 框架内，请主动在名称后追加 '-rare'。"
        "若为非肿瘤或无法判断，写 None。\n"
        '仅返回一行 JSON：\n'
        '{"Cancer_general":"XXX","Cancer":"XXX"}'
    )
    payload = {
        "model": MODEL,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user",   "content": str(summary)}
        ],
        "temperature": 0
    }
    try:
        resp = requests.post(API_URL,
                             headers={"Authorization": f"Bearer {api_key}"},
                             json=payload,
                             timeout=120)
        resp.raise_for_status()
        data = resp.json()["choices"][0]["message"]["content"].strip()
        return json.loads(data).get("Cancer_general", ""), json.loads(data).get("Cancer", "")
    except Exception as e:
        print(f"[WARN] [{sid}] API 失败: {e}")
        return "", ""

# ---------- 主流程 ---------- #
def main():
    df = pd.read_csv(IN_FILE)
    print(f"[{datetime.now():%H:%M:%S}] 读取 {len(df)} 条 Series")

    print("开始 AI 智能判断实体瘤……")
    solid_mask = df["Summary"].apply(is_solid_tumor)
    df = df[solid_mask].copy()
    print(f"[{datetime.now():%H:%M:%S}] 过滤后剩余 {len(df)} 条实体瘤")

    cache = {}
    if os.path.exists(CACHE_FILE):
        cache = json.load(open(CACHE_FILE, encoding="utf-8"))
        print(f"[{datetime.now():%H:%M:%S}] 加载缓存 {len(cache)} 条")

    todo = [(idx, row) for idx, row in df.iterrows() if row["Series_id"] not in cache]
    total = len(todo)
    print(f"[{datetime.now():%H:%M:%S}] 待处理 {total} 条")

    if total == 0:
        print("没有需要处理的数据，直接导出并退出")
        _export(df, cache)
        return

    def worker(payload):
        idx, row = payload
        sid = row["Series_id"]
        summary = row.get("Summary", "")
        key = API_KEYS[idx % len(API_KEYS)]
        thread_id = f"Thread-{int(idx % len(API_KEYS))}"

        print(f"[{datetime.now():%H:%M:%S}] {thread_id} 开始处理  {sid}")
        general, cancer = kimi_cancer_tags(summary, key, sid)
        time.sleep(COOLDOWN)
        return sid, general, cancer

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futures = {ex.submit(worker, item): item for item in todo}
        for done in as_completed(futures):
            sid, general, cancer = done.result()
            cache[sid] = [general, cancer]
            finished = len(cache)
            percent = finished / len(df) * 100
            print(f"[{datetime.now():%H:%M:%S}] "
                  f"[{finished:>4}/{len(df)}  {percent:5.1f}%]  {sid}  ->  {general} | {cancer}")
            json.dump(cache, open(CACHE_FILE, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    _export(df, cache)


def _export(df, cache):
    map_df = pd.DataFrame.from_dict(cache, orient="index",
                                    columns=["Cancer_general", "Cancer"])
    map_df.index.name = "Series_id"
    map_df.reset_index(inplace=True)
    out_df = df.merge(map_df, on="Series_id", how="left")
    out_df.to_csv(OUT_FILE, index=False, encoding="utf-8-sig")
    print(f"[{datetime.now():%H:%M:%S}] 已保存 {len(out_df)} 条 -> {OUT_FILE}")


if __name__ == "__main__":
    main()