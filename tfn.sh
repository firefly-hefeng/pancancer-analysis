#!/usr/bin/env bash

nohup Rscript functional_network_analysis.r > tfn.log 2>&1 &
tail -f tfn.log