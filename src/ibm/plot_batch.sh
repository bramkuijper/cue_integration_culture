#!/usr/bin/env bash

# linux batch scrapt to plot a subset of files

if [ -z "$1" ]
then
    echo "no arguments given, exiting"
    exit
fi

find . -regextype posix-extended -regex "\.\/$1.*[[:digit:]]$"

find . -regextype posix-extended -regex "\.\/$1.*[[:digit:]]$" -exec ./plot_simulation.py {} \;

