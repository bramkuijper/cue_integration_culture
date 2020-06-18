#!/usr/bin/env bash

# delete older files
find . -iname "xaa*" -exec rm -rf {} \;

# check command line args
if [ $# -eq 0 ]
then
    echo "No arguments supplied"
    exit
fi


# splits a runfile in jobs to run locally
split -a5 -l$2 $1

# make executable
find . -iname "xaaa*" -exec chmod +x {} \;
