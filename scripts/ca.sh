#!/usr/bin/env bash

nohup Rscript cellchat_analysis.r > ca.log 2>&1 &
tail -f ca.log