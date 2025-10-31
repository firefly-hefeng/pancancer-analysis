#!/usr/bin/env bash
nohup Rscript annotation_eva.r > annotation_evaluation.log 2>&1 &
tail -f annotation_evaluation.log