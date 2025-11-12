#!/bin/bash  

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/shilangchen/cryptominisat/build/lib
  
# # 编译 equation2xcnf.cpp 并将执行结果保存到 eq.txt 文件中  

filename="eq.txt"
g++ equation2xcnf.cpp -o e2x.o  
./e2x.o > $filename

# # g++ -o voorg varordering_org.cpp ./lib/libcryptominisat5.so.5.11
g++ -I./include -o voorg varordering_org.cpp ./lib/libcryptominisat5.so 
# ./voorg > result24_org.txt

# # #  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/shilangchen/cryptominisat/build/lib  
# g++ -I./include -o vo variable_order.cpp ./lib/libcryptominisat5.so 
# # ./vo > result24.txt


# (./vo | grep "Total time" >> compare2Exhusted.txt && echo "after">> compare2Exhusted.txt) &
# (./voorg | grep "Total time" >> compare2Exhusted.txt && echo "origin">> compare2Exhusted.txt)
# wait    
# echo "done"


g++ -I./include -o voincre variable_incre.cpp ./lib/libcryptominisat5.so 

(./voorg | grep "Total time" >> compare2incre.txt && echo "origin">> compare2incre.txt) &
(./voincre | grep "Total time" >> compare2incre.txt && echo "incremental">> compare2incre.txt)
wait    
echo "done"
