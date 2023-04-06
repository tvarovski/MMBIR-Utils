#!/bin/bash

./filterDirs2.sh 1> outputs/filterDirs2.log 2> outputs/filterDirs2.err
./filterDirs2.sh -f 1> outputs/filterDirs2_f.log 2> outputs/filterDirs2_f.err

echo "Done!"
