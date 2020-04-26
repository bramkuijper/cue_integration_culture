#!/usr/bin/env bash
DIRNAME=`dirname $0`
NCORES=10
MAX_FILES=15

# begin: for testing purposes
#"$DIRNAME/analyze_sims.py" --path=$1 --pattern="sim.*\d+.csv$" --ncores=$NCORES --max_files=$MAX_FILES

# end: for testing purposes
"$DIRNAME/analyze_sims.py" --path=$1 --pattern="sim.*\d+.csv$" --ncores=$NCORES 
