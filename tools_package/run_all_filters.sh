#!/bin/bash

# first run program without -f (filtering) option, redirect output to log file and error file
./filterDirs2.sh 1> outputs/filterDirs2.log 2> outputs/filterDirs2.err &

# then run program with -f (filtering) option, redirect output to log file and error file
./filterDirs2.sh -f 1> outputs/filterDirs2_f.log 2> outputs/filterDirs2_f.err &

echo "Done!"
