#!/usr/bin/env bash
DIRNAME=`dirname $0`
"$DIRNAME/analyze_sims.py" --path=$1 --pattern="sim.*\d+.csv$"
