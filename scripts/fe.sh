#!/usr/bin/env bash

nohup Rscript functional_enrichment.r > fe.log 2>&1 &
tail -f fe.log