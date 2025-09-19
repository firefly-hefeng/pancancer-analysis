#!/usr/bin/env bash

nohup Rscript anno-plot-mgp.r > apmp.log 2>&1 &
tail -f apmp.log