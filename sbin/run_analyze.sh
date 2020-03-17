#!/usr/bin/env bash
DIRNAME=`dirname $0`
NCORES=20
MAX_FILES=5
#"$DIRNAME/analyze_sims.py" --path=$1 --pattern="sim.*\d+.csv$" --ncores=$NCORES --max_files=$MAX_FILES
"$DIRNAME/analyze_sims.py" --path=$1 --pattern="sim.*\d+.csv$" --ncores=$NCORES 
