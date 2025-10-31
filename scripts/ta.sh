#!/usr/bin/env bash

nohup Rscript trajectory_analysis.r > ta.log 2>&1 &
tail -f ta.log