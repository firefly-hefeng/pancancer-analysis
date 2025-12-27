#!/usr/bin/env bash
nohup python -u meta_collect.py >> meta_collect.log 2>&1 &
pid=$!
echo "$pid" | tee -a meta_collect.log
tail -f meta_collect.log