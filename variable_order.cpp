// #include <cryptominisat5/cryptominisat.h>
#include </home/shilangchen/cryptominisat/build/include/cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
using std::vector;
using namespace CMSat;

int main()
{
    SATSolver solver;
    vector<Lit> clause;
    vector<uint32_t> xvars;

    //use 4 threads
    // solver.set_num_threads(4);

    int n=3; //variable's number, start from 0.
    solver.log_to_file("logfile.txt");
    solver.set_allow_otf_gauss();//允许实时高斯消元

    solver.new_vars(n);

    // map<int,vecter<int>> order;
    // // multimap<string, vector<int> > mulmap;

    // for(int i=0;i<block_size;i++){
    //     vector<int> orginial_variable;
    //     orginial_variable.insert(i);
    //     map.insert(make_pair(i, orginial_variable));
    // }

    // int k = block_size+1;
    // for (int i = 1; i <= block_size - 1;i++){
    //     for (int j = i + 1;j<=block_size; j++)
    //     {
    //         vector<int> orginial_variable;
    //         orginial_variable.insert(i);
    //         orginial_variable.insert(j);
    //         map.insert(make_pair(k, orginial_variable));
    //         k++;
    //     }
    // }
    
    
    /*
    void set_xor_detach(bool val);
    void set_find_xors(bool do_find_xors);
    std::vector<std::pair<Lit, Lit> > get_all_binary_xors() const; //get all binary XORs that are = 0
    std::vector<std::pair<std::vector<uint32_t>, bool> > get_recovered_xors(bool xor_together_xors) const; //get XORs recovered. If "xor_together_xors" is TRUE, then xors that share a variable (and ONLY they share them) will be XORed together
    std::vector<uint32_t> get_var_incidence();
        std::vector<uint32_t> get_lit_incidence();
        std::vector<double> get_vsids_scores();    
    bool implied_by(
        const std::vector<Lit>& lits, std::vector<Lit>& out_implied);//一组lit隐含的lits
    */
    
    // add "x1 2 0"
    xvars.push_back(0);
    xvars.push_back(1);
    xvars.push_back(2);
    //rhs=true:clause is true
    solver.add_xor_clause(xvars, false);


    //add "2 0"
    clause.clear();
    
    Lit lit;
    std::cout<<lit.isorigin;
    
    clause.push_back(Lit(0,true));
    std::cout<<clause[0].isorigin;
    std::cout<<clause[0].var_map[0];
    // clause.push_back(Lit(2, false,0,1));
    solver.add_clause(clause);

    // clause.clear();
    // clause.push_back(Lit(2, false,0,1));
    // clause.push_back(Lit(1,true));
    // solver.add_clause(clause);

    // clause.clear();
    // clause.push_back(Lit(2, true,0,1));
    // clause.push_back(Lit(0,false));
    // clause.push_back(Lit(1,false));
    // solver.add_clause(clause);

    lbool ret = solver.solve();
    assert(ret == l_True);
    for(int i=0;i<n;i++){
        std::cout<<solver.get_model()[i]<<", ";
    }
    std::cout<< std::endl;

    


    return 0;
}