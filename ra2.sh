#!/usr/bin/env bash

nohup Rscript ra2.r > ra2.log 2>&1 &
tail -f ra2.log