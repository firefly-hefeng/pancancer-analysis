#!/usr/bin/env bash
nohup Rscript annotation_eva11.r > kegg_11_12.log 2>&1 &
tail -f kegg_11_12.log

