#!/usr/bin/env bash

nohup Rscript run_clustering_annotation.R > ra.log 2>&1 &
tail -f ra.log