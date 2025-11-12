#!/bin/bash

# cd ./cryptominisat/mytest

blocksize=6

totalTimes=10

for index in $(seq 1 $totalTimes); do

    cnfName="./data2/cnf/n${blocksize}/example${blocksize}_${index}.cnf"
    backboneName="./data2/backbone/n${blocksize}/example${blocksize}_${index}.cnf"
    ./utils/cadiback/cadiback -q $cnfName > $backboneName

done

echo "done"

# ./find_backbone.sh