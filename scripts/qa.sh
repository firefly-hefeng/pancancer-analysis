#!/usr/bin/env bash
nohup Rscript quality.r > quality.log 2>&1 &
PID=$!
echo "Rscript PID: $PID" >> quality.log
tail -f quality.log# Commit 10: feat: add quality control pipeline - 1775143673
# Commit 23: fix: correct path configuration in scripts - 1775143688
# Commit 36: feat: optimize memory usage in collectors - 1775143700
# Commit 49: fix: correct KPMP data fetching - 1775143712
# Commit 10: feat: add quality control pipeline - 1775143741
