#!/usr/bin/env bash

nohup Rscript anno_plot.r > ap.log 2>&1 &
tail -f ap.log