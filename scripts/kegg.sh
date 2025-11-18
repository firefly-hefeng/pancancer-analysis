#!/usr/bin/env bash
nohup Rscript annotation_eva1.r > kegg_11_10.log 2>&1 &
tail -f kegg_11_10.log

