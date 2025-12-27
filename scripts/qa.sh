#!/usr/bin/env bash
nohup Rscript quality.r > quality.log 2>&1 &
PID=$!
echo "Rscript PID: $PID" >> quality.log
tail -f quality.log