#!/bin/bash

# cd ./cryptominisat/mytest

blocksize=24

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/shilangchen/cryptominisat/build/lib
g++ -Dblock_size=$blocksize -I../build/include -o solve_randAssign solve_randAssign.cpp ../build/lib/libcryptominisat5.so.5.11


totalTimes=60
resultfile="./data2/all_solutions${blocksize}.txt"


# # 并行
max_jobs=20  # 最大并行作业数

seq 21 $totalTimes | xargs -n 1 -P $max_jobs bash -c '
    index=$0

    filename="./data2/eq/n'${blocksize}'/eq'${blocksize}'_${index}.txt"
    solution="./data2/solution/n'${blocksize}'/solution'${blocksize}'_${index}.txt"
    logname_rand="./data2/log/n'${blocksize}'/logfile_rand'${blocksize}'_${index}.log"

    ./solve_randAssign $filename $logname_rand > $solution
'


# for index in $(seq 41 $totalTimes); do

#     filename="./data2/eq/n${blocksize}/eq${blocksize}_${index}.txt"
#     logname_rand="./data2/log/n${blocksize}/logfile_rand${blocksize}_${index}.log"

#     # ./solve_randAssign $filename $logname_rand >> $logname_rand
#     echo "eq${blocksize}_${index}_randAssign: ">> $resultfile
#     # grep "success" $logname_rand >> $resultfile
#     # grep "Total time" $logname_rand >> $resultfile
#     ./solve_randAssign $filename $logname_rand >> $resultfile
# done

echo "done"

# ./solve_allSolutions.sh