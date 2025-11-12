#!/bin/bash

# cd ./cryptominisat/mytest

blocksize=30

g++ -Dblock_size=$blocksize equation2xcnf.cpp -o e2x.o

totalTimes=11

for index in $(seq 11 $totalTimes); do

    eqname="./data2/eq/n${blocksize}/eq${blocksize}_${index}.txt"
    fileName="./data2/cnf/n${blocksize}/example${blocksize}_${index}"
    
    ./e2x.o $fileName > $eqname

done
echo "done"

# ./generate_equation.sh

# ./e2x.o > ./data/example27.cnf
# minisat ./data/example27.cnf >./data/log/result.log


