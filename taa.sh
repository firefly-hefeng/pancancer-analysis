#!/usr/bin/env bash

nohup Rscript trajectory-aa.r > taa.log 2>&1 &
tail -f taa.log