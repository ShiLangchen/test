#!/bin/bash

# cd ./cryptominisat/mytest

blocksize=24

# g++ -Dblock_size=$blocksize equation2xcnf.cpp -o e2x.o


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/shilangchen/cryptominisat/build/lib
g++ -Dblock_size=$blocksize -I../build/include -o solve_org solve_org.cpp ../build/lib/libcryptominisat5.so.5.11
g++ -Dblock_size=$blocksize -I../build/include -o solve_randAssign solve_randAssign.cpp ../build/lib/libcryptominisat5.so.5.11


totalTimes=20
resultfile="./data/compare${blocksize}.txt"
max_jobs=20  # 最大并行作业数


# for index in $(seq 1 $totalTimes); do


seq 1 $totalTimes | xargs -n 1 -P $max_jobs bash -c '
    index=$0

    filename="./data/eq/n'${blocksize}'/eq'${blocksize}'_${index}.txt"
    logname_org="./data/log/n'${blocksize}'/logfile_org'${blocksize}'_${index}.log"
    logname_rand="./data/log/n'${blocksize}'/logfile_rand'${blocksize}'_${index}.log"

    ./solve_randAssign $filename $logname_rand >> $logname_rand
    echo "eq'${blocksize}'_${index}_randAssign: ">> '$resultfile'
    grep "success" $logname_rand >> '$resultfile'
    grep "Total time" $logname_rand >> '$resultfile'

    ./solve_org $filename $logname_org >> $logname_org
    echo "eq'${blocksize}'_${index}_origin: ">> '$resultfile'
    grep "Total time" $logname_org >> '$resultfile'    
'
# done
echo "done"


# ./solve_equation.sh