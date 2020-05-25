#!/usr/bin/env bash

# linux batch scrapt to plot a subset of files

if [ -z "$1" ]
then
    echo "no arguments given, exiting"
    exit
fi

find . -regextype posix-extended -regex "\.\/$1.*[[:digit:]].csv$"

find . -regextype posix-extended -regex "\.\/$1.*[[:digit:]].csv$" -exec ./plot_simulation.py {} \;

