#include <cryptominisat5/cryptominisat.h>
// #include </home/shilangchen/cryptominisat/build/include/cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
#include <cassert> 

#include <iostream>
#include <fstream>  
#include <sstream>  
#include <stdio.h>
#include <ctime>
#include <set>
using std::vector;
using namespace CMSat;

// #define block_size 30

int state_alpha[block_size][block_size + block_size * block_size / 2 + 1] = { 0 };

void readFile(const std::string &fileName){
    std::string line;  
    int num = 0;  
    int row = 0, col = 0;  
    std::ifstream file;  
    file.open(fileName);  
    if (file.is_open()) {  
        while (getline(file, line)) {  
            col = 0;  
            for (int i = 0; i < line.size(); i++) {  
                if (isdigit(line[i])) {  
                    num = num * 10 + (line[i] - '0');  
                } else if (line[i] == ' ') {  
                    state_alpha[row][col] = num;  
                    num = 0;  
                    col++;  
                }  
            }  
            state_alpha[row][col] = num;  
            row++;  
        }  
        file.close();  
    }else {  
        std::cout << "Unable to open file";  
    }  
}

int main()
{
    std::string fileName="eq.txt"; 
    // std::cin >> fileName; 
    readFile(fileName); 

    SATSolver solver;
    vector<Lit> clause;
    vector<uint32_t> xvars;
    //use 4 threads
    // solver.set_num_threads(4);
    solver.log_to_file("logfile.txt");
    solver.set_allow_otf_gauss();//����ʵʱ��˹��Ԫ

    solver.set_max_time(24*60*60);


    solver.set_xor_detach(1);//bool val
    solver.set_verbosity(1); //default is 0, silent
    solver.set_verbosity_detach_warning(1); //default is 0, silent

    int n = block_size + block_size * (block_size - 1) / 2; //variable's number, start from 0.
    // int m=block_size + 3 * block_size * (block_size - 1) / 2; //�Ӿ���
    
    solver.new_vars(n);

	std::set<int> occur;
    for (int i = 0; i < block_size;i++){
        for (int j = 0; j < block_size + block_size * (block_size-1) / 2;j++){
            if(state_alpha[i][j]==1){
                // cout << j + 1 << " ";
                xvars.push_back(j);
				if (j >= block_size && occur.find(j)==occur.end()){
					occur.insert(j);
				}
			}
        }
        if(state_alpha[i][block_size+block_size*block_size/2]==0){
            solver.add_xor_clause(xvars, false);
        }else{
            solver.add_xor_clause(xvars, true);//rhs=true:clause is true
        }
        // cout << "0" << endl;
        xvars.clear();
    }

    int k = block_size;
    for (int i = 0; i < block_size - 1; i++){
        for (int j = i + 1; j < block_size; j++)
        {
			if(occur.find(k)!=occur.end()){
				// cout <<"-"<< k << " " << i << " 0" << endl;
            	// cout <<"-"<< k << " " << j << " 0" << endl;	
            	// cout << k << " -" << i <<" -" << j << " 0" << endl;
                clause.clear();
                clause.push_back(Lit(i,false));
                clause.push_back(Lit(k,true,i,j));
                solver.add_clause(clause);

                clause.clear();
                clause.push_back(Lit(j,false));
                clause.push_back(Lit(k,true,i,j));
                solver.add_clause(clause);

                clause.clear();
                clause.push_back(Lit(k,false,i,j));
                clause.push_back(Lit(i,true));
                clause.push_back(Lit(j,true));
                solver.add_clause(clause);
			}
            k++;
        }
    }

    lbool ret = solver.solve();
    assert(ret == l_True);
    std::cout<<ret<< std::endl;
    for(int i=0;i<n;i++){
        std::cout<<solver.get_model()[i]<<", ";
    }
    std::cout<< std::endl;
    // std::cout<< solver.get_vsids_scores();
    solver.print_stats();

    return 0;
}