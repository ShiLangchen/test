#include <cstdlib>
#include <cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>  // 对于 std::put_time

using namespace std;
using std::vector;
using namespace CMSat;

// g++ main_lowerBound.cpp -o mainlb -lcryptominisat5
// ./mainlb
int maxIndex=0;
int R = 10;
int lowW=40;
int upW=41;
int apxMaxVarNum=10000;


void creat_sat_Sequence_encoding_given_original_method(SATSolver* solver, int weight, int fuzhu_pos, int* index, int indexlen)
{

	vector<Lit> clause;

	int fuzhupos = fuzhu_pos;
	//actually 48 sbox to estimate a link  96bool to weight weight_
	int num = indexlen;
	//single state without combination

	  //cons 1
	//加入标标准准的序列约束

	//cons 1
	clause.clear();
	clause.push_back(Lit(index[0], true));
	clause.push_back(Lit(fuzhupos, false));
	(*solver).add_clause(clause);

	//clause.clear();
	//clause.push_back(Lit(321, true));
	//clause.push_back(Lit(fuzhupos + 2, false));
	//(*solver).add_clause(clause);
	//cons 2
	for (size_t i = 1; i < weight; i++)
	{
		clause.clear();
		clause.push_back(Lit(fuzhupos + i, true));//与j是第几轮的有关系0-8对于0123
		(*solver).add_clause(clause);
	}

	int k = weight;
	for (size_t i = 1; i < num - 1; i++)
	{
		//cons 3
		//all condition should take two branch into consideration
		clause.clear();
		clause.push_back(Lit(index[i], true));
		clause.push_back(Lit(fuzhupos + i * k, false));//s盒的数量是round*64 
		(*solver).add_clause(clause);

		//clause.clear();
		//clause.push_back(Lit(index[i][1], true));
		//clause.push_back(Lit(fuzhupos + i * k + 2, false));//s盒的数量是round*64 
		//(*solver).add_clause(clause);

		//cons 4


		clause.clear();
		clause.push_back(Lit(fuzhupos + (i - 1) * k, true));
		clause.push_back(Lit(fuzhupos + (i)*k, false));
		(*solver).add_clause(clause);


		/*clause.clear();
		clause.push_back(Lit(fuzhupos + (i - 1) * k + 2, true));
		clause.push_back(Lit(fuzhupos + (i)*k + 2, false));
		(*solver).add_clause(clause);*/
		//just 2 or 3 is adequate
		for (size_t j = 1; j < k; j++)//这里应该取k，而不是k-1，这会导致最后一行没有递归
		{
			//cons 5
			//need to modify

			clause.clear();
			clause.push_back(Lit(index[i], true));
			clause.push_back(Lit(fuzhupos + (i - 1) * k + j - 1, true));
			clause.push_back(Lit(fuzhupos + (i)*k + j, false));
			(*solver).add_clause(clause);



			//cons 6
			clause.clear();
			clause.push_back(Lit(fuzhupos + (i - 1) * k + j, true));
			clause.push_back(Lit(fuzhupos + (i)*k + j, false));
			(*solver).add_clause(clause);

		}

		//cons 7
		clause.clear();
		clause.push_back(Lit(index[i], true));
		clause.push_back(Lit(fuzhupos + (i - 1) * k + k - 1, true));//s盒的数量是round*64 
		(*solver).add_clause(clause);


        
	}
	//cons 8

	clause.clear();
	clause.push_back(Lit(index[num - 1], true));
	clause.push_back(Lit(fuzhupos + (num - 2) * k + k - 1, true));//s盒的数量是round*64 
    maxIndex = fuzhupos + (num - 2) * k + k;
	(*solver).add_clause(clause);
    
}


void creat_sat_Sequence_encoding_neg(SATSolver* solver,  int weight, int neg_fuzhu_pos, int* index, int indexlen)
{
    // sum(x) > m --> sum(x) >= m + 1 --> sum(-x) <= n - m - 1

	// fuzhu_pos = 8000, index, indexlen = round * 48
    // neg_fuzhu_pos = 200000
	vector<Lit> clause;

    neg_fuzhu_pos = maxIndex + 1;

	int fuzhupos = neg_fuzhu_pos;
	//actually 48 sbox to estimate a link  96bool to weight weight_
	int num = indexlen; // round * 48
	//single state without combination

	// cons 1
    // 第一个布尔变量与第一个辅助变量的关系
    clause.clear();
    clause.push_back(Lit(index[0], false)); // x0 -> ~s0,0
    clause.push_back(Lit(fuzhupos, false));
    (*solver).add_clause(clause);

    // cons 2
    // 对前 weight 个辅助变量添加单文字子句，固定这些变量的值为 0。
    // 取值为 0 表示这些变量不能参与计数。
    for (int i = 1; i <= weight; i++)
    {
        clause.clear();
        clause.push_back(Lit(fuzhupos + i, true)); // ~s[i]
        (*solver).add_clause(clause);
    }

    int k = weight + 1; // 辅助变量的分支大小
    k = num - k - 1;
    for (int i = 1; i < num - 1; i++)
    {
        // cons 3
        // 当前布尔变量与其对应的辅助变量关系
        clause.clear();
        clause.push_back(Lit(index[i], false)); // xi -> ~si,0
        clause.push_back(Lit(fuzhupos + i * k, false));
        (*solver).add_clause(clause);

        // cons 4
        // 辅助变量之间的递归关系
        clause.clear();
        clause.push_back(Lit(fuzhupos + (i - 1) * k, true)); // s(i-1) -> ~si
        clause.push_back(Lit(fuzhupos + i * k, false));
        (*solver).add_clause(clause);

        for (int j = 1; j < k; j++)
        {
            // cons 5
            // xi 与前一个辅助变量的联合约束
            clause.clear();
            clause.push_back(Lit(index[i], false));                   // xi
            clause.push_back(Lit(fuzhupos + (i - 1) * k + j - 1, true)); // s(i-1,j-1)
            clause.push_back(Lit(fuzhupos + i * k + j, false));      // -> ~si,j
            (*solver).add_clause(clause);

            // cons 6
            // 辅助变量的递归关系
            clause.clear();
            clause.push_back(Lit(fuzhupos + (i - 1) * k + j, true)); // s(i-1,j)
            clause.push_back(Lit(fuzhupos + i * k + j, false));      // -> ~si,j
            (*solver).add_clause(clause);
        }

        // cons 7
        // 当前布尔变量否定最后一个辅助变量
        clause.clear();
        clause.push_back(Lit(index[i], false));                       // xi
        clause.push_back(Lit(fuzhupos + (i - 1) * k + k - 1, true)); // -> ~s(i-1,k-1)
        (*solver).add_clause(clause);
    }

    // cons 8
    // 最后一个布尔变量与最后一个辅助变量的关系
    clause.clear();
    clause.push_back(Lit(index[num - 1], false));                      // x(num-1)
    clause.push_back(Lit(fuzhupos + (num - 2) * k + k - 1, true));    // -> ~s(num-2,k-1)
    solver->add_clause(clause);
}

void at_least_m(SATSolver* solver, int k, int m, int fuzhu_pos,  int indexlen){
	vector<Lit> clause;

	int num = indexlen; // round * 48

	// r_n-1,m-1
	clause.push_back(Lit(fuzhu_pos + (num - 1) * k + m -1, false)); 

	(*solver).add_clause(clause);
}

void add_constrain_linear_present(SATSolver* solver, int x_pos)
{
	int lineara[64] = { 0 };

	lineara[0] = 0;
	lineara[63] = 63;
	for (size_t i = 1; i < 63; i++)
	{
		lineara[i] = i * 16 % 63;
	}

	vector<Lit> clause;
	/*clause.clear();
	(*solver).add_clause(clause);*/
	for (size_t i = 0; i < 64; i++)
	{
		clause.clear();
		//赋值字句，需要两条

		clause.push_back(Lit(x_pos +64+ i, true));
		clause.push_back(Lit(x_pos  + lineara[i], false));

		(*solver).add_clause(clause);

		clause.clear();
		//赋值字句，需要两条

		clause.push_back(Lit(x_pos + 64 + i, false));
		clause.push_back(Lit(x_pos + lineara[i], true));

		(*solver).add_clause(clause);
	}
}


void add_constrain_sbox_present_3weight(SATSolver* solver, int x_pos, int outx_pos, int consnum, int** cons8)
{
	vector<Lit> clause;

	//16个s盒
	for (size_t i = 0; i < 16; i++)
	{
		//32个s盒每轮 我们的目的是找到一条解就行了，不需要全部找到

		for (size_t j = 0; j < consnum; j++)//每个s盒包含如此之多的约束
		{
			//1表示有 2表示反
			clause.clear();

			//前16个s盒的结果存ori 以及变量
			for (size_t k = 0; k < 4; k++)
			{
				//一个是向前的bit，一个是向后的bit,前面是原来的，要减去整个16*64，当然大家都会减去一个exchange
				if (cons8[j][k] == 1)
				{
					clause.push_back(Lit(x_pos + i + k * 16, false));//与j是第几轮的有关系0-8对于0123

				}
				else  if (cons8[j][k] == 2)
				{
					clause.push_back(Lit(x_pos + i + k * 16, true));//与j是第几轮的有关系0-8对于0123

				}


				//变换的结果存储在临时变量里面，索引照用就行，实际上这个临时变量只使用一半
				//可以考虑临时变量后半段是ori的变换结果 前半段是multi8的变换结果
				if (cons8[j][k + 4] == 1)
				{

					clause.push_back(Lit(x_pos + i + k * 16 + 64, false));//与j是第几轮的有关系0-8对于0123

				}
				else  if (cons8[j][k + 4] == 2)
				{
					clause.push_back(Lit(x_pos + i + k * 16 + 64, true));//与j是第几轮的有关系0-8对于0123

				}

			}

			for (size_t k = 0; k < 3; k++)
			{
				//差分概率

				if (cons8[j][k + 8] == 1)
				{

					clause.push_back(Lit(outx_pos + i * 3 + k, false));//与j是第几轮的有关系0-8对于0123

				}
				else  if (cons8[j][k + 8] == 2)
				{
					clause.push_back(Lit(outx_pos + i * 3 + k, true));//与j是第几轮的有关系0-8对于0123

				}

			}







			(*solver).add_clause(clause);



		}


	}

}


std::vector<std::vector<int>> parseExpression_aim_lengh(const std::string& filename,int lengh,int mid) {
	std::vector<std::vector<int>> expressions;
	std::ifstream file(filename);
	if (!file) {
		std::cerr << "Failed to open file: " << filename << std::endl;
		return expressions;
	}

	std::string line;
	while (std::getline(file, line)) {
		//    std::vector<int> expression(16, 0); // Initialize with all zeros

			// Find expressions within parentheses
		size_t start = line.find('(');
		size_t end;
		while (start != std::string::npos) {
			end = line.find(')', start);
			if (end == std::string::npos) {
				std::cerr << "Invalid expression: " << line << std::endl;
				return expressions;
			}
			std::vector<int> expression(lengh, 0); // Initialize with all zeros
			std::string expr = line.substr(start + 1, end - start);

			size_t length = expr.length();
			size_t i = 0;
			while (i < length) {
				char c = expr[i];
				/*   if (c == 'x' || c == 'y') {
					   int index = 0;
					   if (i + 2 < length && expr[i + 1] == '\'' && expr[i + 2] == '\'') {
						   index = (c == 'x') ? 8 : 16;
						   i += 3;
					   }
					   else {
						   index = (c == 'x') ? 0 : 8;
						   i++;
					   }
					   if (i < length && expr[i] >= '1' && expr[i] <= '8') {
						   index += expr[i] - '1';
						   expression[index] = 1;
					   }
				   }*/

				if (c == 'x')
				{
					if (expr[i + 2] == '\'')//表示反
					{
						expression[expr[i + 1] - '1'] = 2;
					}
					else
					{
						expression[expr[i + 1] - '1'] = 1;
					}

				}
				if (c == 'y')
				{
					if (expr[i + 2] == '\'')//表示反
					{
						expression[mid + expr[i + 1] - '1'] = 2;
					}
					else
					{
						expression[mid + expr[i + 1] - '1'] = 1;
					}

				}


				i++;
			}



			expressions.push_back(expression);
			start = line.find('(', end);
		}

		//  expressions.push_back(expression);
	}

	return expressions;
}


void read_constrant_from_file(std::string filename,int* num,int** result,int cnf_length,int mid)
{
	//一个一个来做
  // std::string filename = "eq00.txt"; // 你的表达式文本文件的路径

	std::vector<std::vector<int>> expressions = parseExpression_aim_lengh(filename, cnf_length, mid);
	int count = 0;
	// 打印解析后的表达式数组
	for (const auto& expression : expressions) {
		for (int num : expression) {
			//  std::cout << num << " ";
		}
		count++;
		//   printf("\n %d \n", count);
	   //    std::cout << std::endl;
	}
	*num = count;
	//  printf("\n %d \n", count);

	for (int i = 0; i < count; i++)
	{
		for (size_t j = 0; j < cnf_length; j++)
		{
			result[i][j] = expressions[i][j];
		}
	}

	for (int i = 0; i < count; i++)
	{
		//    printf("\n");
		for (size_t j = 0; j < 10; j++)
		{
			//    printf(" %d ", constrain_ascon10[i][j]);
		}
	}
}


bool creat_sat_present_differential_xuliebianma(int round_, int thread_, int out_time, int weight)
{

	int** cnf_sbox;
	cnf_sbox = (int**)malloc(150 * 8);//100条约束
	for (size_t i = 0; i < 150; i++)
	{
		cnf_sbox[i] = (int*)malloc(20 * 8);//一般约束比较短
	}
	int sbox_num_cons;
	
	read_constrant_from_file("present_ddt_3one.txt", &sbox_num_cons, cnf_sbox, 11, 5);
	//基于贪心算法的约减工作可以放在后续一起做

	SATSolver solver;

	solver.set_num_threads(thread_);

	solver.new_vars(round_ * (apxMaxVarNum));

    // TRACE
    string logfilename = "r" +to_string(R)+ "w"+ to_string(upW)+"lb.log";
	solver.log_to_file(logfilename);


	int round = round_;

	vector<Lit> clause;

	clause.clear();
	for (size_t i = 0; i < 64; i++)
	{
		clause.push_back(Lit(i, false));
	}

	solver.add_clause(clause);


	for (size_t i = 0; i < round_; i++)
	{
		add_constrain_sbox_present_3weight(&solver, i * 128, 5000 + i * 48, sbox_num_cons, cnf_sbox);
	
		add_constrain_linear_present(&solver, i * 128 + 64);
	}


	int index[3000] = { 0 };
	for (size_t i = 0; i < round * 3 * 16; i++)
	{
		index[i] = 5000 + i;
	}
	
	creat_sat_Sequence_encoding_given_original_method(&solver, weight, 8000, index, round * 48);

    // creat_sat_Sequence_encoding_neg(&solver, lowW, 20000, index, round * 48);
	// at_least_m(&solver, weight, lowW, 8000, round * 48);

	//	printf("\n 开始求解\n");

	clock_t t1, t2;


	t1 = clock();
	lbool ret = solver.solve();

	// at_least_m(&solver, weight, lowW, 8000, round * 48);
	// ret = solver.solve();

	t2 = clock();
	std::ostringstream oss;
	double elapsed_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
	oss << std::fixed << std::setprecision(2) << elapsed_time;  // 设置小数点后保留两位
	std::string time_string = oss.str();


	if (ret == l_False)
	{
		printf("\n 序列 该子模型无解！ %d %d\n", round_, weight);

		std::string filename = "sequnence_code_presentlb.csv";
		std::ofstream file(filename, std::ios_base::app);

		//轮数 重量，时间 有无解 0无解 1有解
		file << round_ << "," << weight << "," << time_string <<","<< 0 << "\n";  // 以逗号分隔每列，换行结束每行

		file.close();
		//  exit(0);
		//system("pause");
		//存储时间
        solver.print_stats();

		return false;
	}
	else if (ret == l_True)
	{

		printf("\n 序列 该模型有解！ %d %d\n", round_, weight);

		std::string filename = "sequnence_code_presentlb.csv";
		std::ofstream file(filename, std::ios_base::app);

		//轮数 重量，时间 有无解 0无解 1有解
		file << round_ << "," << weight << "," << time_string << "," << 1 << "\n";  // 以逗号分隔每列，换行结束每行
		file.close();

        solver.print_stats();
		return true;
	}
	return true;
}


void main_control_center_present(int round)
{

	bool b1 = false;

	clock_t t1, t2;

	t1 = clock();

	int weight = upW;
	
	b1 = creat_sat_present_differential_xuliebianma(round, 4, 0, weight);
	
	b1 = false;

	 t2 = clock();

	 std::ostringstream oss;
	 double elapsed_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
	 oss << std::fixed << std::setprecision(2) << elapsed_time;  // 设置小数点后保留两位
	 std::string time_string = oss.str();

	std::string filename = "sequnence_all_presentlb.csv";
 	std::ofstream file(filename, std::ios_base::app);

 	//轮数 重量，时间 有无解 0无解 1有解
 	file << round << "," << weight << "," << time_string << "\n";  // 以逗号分隔每列，换行结束每行

 	file.close();

}


int main()
{
	main_control_center_present(R);
}