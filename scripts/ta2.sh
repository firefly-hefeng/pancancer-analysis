#!/usr/bin/env bash

nohup Rscript trajectory_analysis2.r > ta2.log 2>&1 &
tail -f ta2.log