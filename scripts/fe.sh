#!/usr/bin/env bash

nohup Rscript functional_enrichment.r > fe.log 2>&1 &
tail -f fe.log# Commit 11: feat: implement cell type annotation - 1775143675
# Commit 24: chore: update shell script permissions - 1775143688
# Commit 37: fix: handle edge cases in annotation - 1775143701
# Commit 50: docs: add troubleshooting guide - 1775143713
