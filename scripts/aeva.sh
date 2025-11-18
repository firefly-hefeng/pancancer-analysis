#!/usr/bin/env bash
nohup Rscript annotation_eva_all.r > annotation_all.log 2>&1 &
tail -f annotation_all.log

