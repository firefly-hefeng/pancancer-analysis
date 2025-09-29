#!/usr/bin/env bash

nohup Rscript conserved_specific_analysis.r > tcs.log 2>&1 &
tail -f tcs.log