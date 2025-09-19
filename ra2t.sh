#!/usr/bin/env bash

nohup Rscript ra2t2.r > ra2t2.log 2>&1 &
tail -f ra2t2.log