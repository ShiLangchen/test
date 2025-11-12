#!/bin/bash  
  
# ������Ҫ�ȽϵĲ����б�  
# parameters=("" "--polar rnd" "--polar stable" "--restart geom" "--restart glue" "--restart luby" "--rstfirst 50" "--rstfirst 100" 
#             "--rstfirst 150" "--gluehist 40" "--gluehist 60") 
# parameters=("--branchstr rand" "--branchstr vsids" "--branchstr vmtf" "--branchstr vmtf+vsids" "--branchstr vsids+vmtf" "--branchstr vsids+rand" "--branchstr rand+vsids" "--branchstr vmtf+rand" "--branchstr rand+vmtf" "--branchstr vmtf+vsids+rand" ) 

parameters=("--branchstr vmtf+vsids+rand" "--branchstr vmtf+rand+vsids" "--branchstr vsids+rand+vmtf" "--branchstr vsids+vmtf+rand" "--branchstr rand+vmtf+vsids" "--branchstr rand+vsids+vmtf")

filename="test30_1.cnf"
echo $filename >> output.txt
# # ��������ÿ�������Ľű�  
for param in "${parameters[@]}"  
do  
    # �����ӽ��̲����нű�  
    (./cryptominisat5 --verb 1 --maxtime 2000 --autodisablegauss 0 --detachxor 1 $param $filename | grep "Total time" >> output.txt && echo "$param" >> output.txt) &  
done  


# files=("test24_1.cnf" "test24_2.cnf" "test24.cnf")  
# for file in "${files[@]}"  
# do  
#     (./cryptominisat5 --verb 1 --autodisablegauss 0 --detachxor 1 $file | grep "Total time" >> output.txt && echo "$file" >> output.txt) &  
# done  
  
wait    
echo "done"

# ./cryptominisat5 --verb 1 --autodisablegauss 0 --detachxor 1 --branchstr vmtf ./test30_2.cnf > test30_2.log
# ./cryptominisat5 --verb 1 --autodisablegauss 0 --detachxor 1 ./test24_1.cnf > test24_1.log
