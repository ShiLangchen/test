#include <cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
#include<time.h>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>  // 对于 std::put_time

using namespace std;
using std::vector;
using namespace CMSat;


void creat_sat_Sequence_encoding_given_original_method(SATSolver* solver, int weight, int fuzhu_pos, int* index, int indexlen)
{
	// fuzhu_pos = 8000, index, indexlen = round * 48

	vector<Lit> clause;

	int fuzhupos = fuzhu_pos;
	//actually 48 sbox to estimate a link  96bool to weight weight_
	int num = indexlen; // round * 48
	//single state without combination

	//cons 1
	//加入标标准准的序列约束
	// index[0] -> fuzhupos
	clause.clear();
	clause.push_back(Lit(index[0], true));
	clause.push_back(Lit(fuzhupos, false));
	(*solver).add_clause(clause);

	//clause.clear();
	//clause.push_back(Lit(321, true));
	//clause.push_back(Lit(fuzhupos + 2, false));
	//(*solver).add_clause(clause);

	//cons 2
	// 对前 weight 个辅助变量添加单文字子句，固定这些变量的值为 1。
	// ~(fuzhupos + i)	
	for (size_t i = 1; i < weight; i++)
	{
		clause.clear();
		clause.push_back(Lit(fuzhupos + i, true)); //与j是第几轮的有关系0-8对于0123
		(*solver).add_clause(clause);
	}

	int k = weight;
	for (size_t i = 1; i < num - 1; i++)
	{
		//cons 3
		//all condition should take two branch into consideration
		// 每个变量 index[i] 与其对应的辅助变量建立关系。
		// index[i] -> fuzhupos + i * k
		clause.clear();
		clause.push_back(Lit(index[i], true));
		clause.push_back(Lit(fuzhupos + i * k, false));//s盒的数量是round*64 
		(*solver).add_clause(clause);

		//clause.clear();
		//clause.push_back(Lit(index[i][1], true));
		//clause.push_back(Lit(fuzhupos + i * k + 2, false));//s盒的数量是round*64 
		//(*solver).add_clause(clause);

		//cons 4
		// 辅助变量之间建立递归关系。
		// fuzhupos + (i-1)*k -> fuzhupos + i*k
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
			// 辅助变量与布尔变量之间建立联合关系，包含多个分支。
			// (index[i] ∨ fuzhupos + (i-1)*k + j-1) -> fuzhupos + i*k + j
			clause.clear();
			clause.push_back(Lit(index[i], true));
			clause.push_back(Lit(fuzhupos + (i - 1) * k + j - 1, true));
			clause.push_back(Lit(fuzhupos + (i)*k + j, false));
			(*solver).add_clause(clause);

			//cons 6
			// 递归辅助变量之间的直接关系。
			// fuzhupos + (i-1)*k + j -> fuzhupos + i*k + j
			clause.clear();
			clause.push_back(Lit(fuzhupos + (i - 1) * k + j, true));
			clause.push_back(Lit(fuzhupos + (i)*k + j, false));
			(*solver).add_clause(clause);

		}

		//cons 7
		// index[i] -> ~(fuzhupos + (i-1)*k + k-1)
		clause.clear();
		clause.push_back(Lit(index[i], true));
		clause.push_back(Lit(fuzhupos + (i - 1) * k + k - 1, true));//s盒的数量是round*64 
		(*solver).add_clause(clause);
	}

	//cons 8
	// ~index[num-1] ∨ ~ fuzhupos + (num-2)*k + k-1
	clause.clear();
	clause.push_back(Lit(index[num - 1], true));
	clause.push_back(Lit(fuzhupos + (num - 2) * k + k - 1, true));//s盒的数量是round*64 
	(*solver).add_clause(clause);
}

void creat_sat_Sequence_encoding_increment_method(SATSolver* solver, int weight, int fuzhu_pos, int* index, int indexlen)
{
	// fuzhu_pos = 8000, index, indexlen = round * 48

	vector<Lit> clause;

	int fuzhupos = fuzhu_pos;
	//actually 48 sbox to estimate a link  96bool to weight weight_
	int num = indexlen; // round * 48
	//single state without combination

	//cons 1
	//加入标标准准的序列约束
	// index[0] -> fuzhupos
	clause.clear();
	clause.push_back(Lit(index[0], true));
	clause.push_back(Lit(fuzhupos, false));
	(*solver).add_clause(clause);

	//clause.clear();
	//clause.push_back(Lit(321, true));
	//clause.push_back(Lit(fuzhupos + 2, false));
	//(*solver).add_clause(clause);

	//cons 2
	// 对前 weight 个辅助变量添加单文字子句，固定这些变量的值为 1。
	// ~(fuzhupos + i)	
	for (size_t i = 1; i < weight; i++)
	{
		clause.clear();
		clause.push_back(Lit(fuzhupos + i, true)); //与j是第几轮的有关系0-8对于0123
		(*solver).add_clause(clause);
	}

	int k = weight;
	for (size_t i = 1; i < num - 1; i++)
	{
		//cons 3
		//all condition should take two branch into consideration
		// 每个变量 index[i] 与其对应的辅助变量建立关系。
		// index[i] -> fuzhupos + i * k
		clause.clear();
		clause.push_back(Lit(index[i], true));
		clause.push_back(Lit(fuzhupos + i * num, false));//s盒的数量是round*64 
		(*solver).add_clause(clause);

		//clause.clear();
		//clause.push_back(Lit(index[i][1], true));
		//clause.push_back(Lit(fuzhupos + i * k + 2, false));//s盒的数量是round*64 
		//(*solver).add_clause(clause);

		//cons 4
		// 辅助变量之间建立递归关系。
		// fuzhupos + (i-1)*k -> fuzhupos + i*k
		clause.clear();
		clause.push_back(Lit(fuzhupos + (i - 1) * num, true));
		clause.push_back(Lit(fuzhupos + (i) * num, false));
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
			// 辅助变量与布尔变量之间建立联合关系，包含多个分支。
			// (index[i] ∨ fuzhupos + (i-1)*k + j-1) -> fuzhupos + i*k + j
			clause.clear();
			clause.push_back(Lit(index[i], true));
			clause.push_back(Lit(fuzhupos + (i - 1) * num + j - 1, true));
			clause.push_back(Lit(fuzhupos + (i)*num + j, false));
			(*solver).add_clause(clause);

			//cons 6
			// 递归辅助变量之间的直接关系。
			// fuzhupos + (i-1)*k + j -> fuzhupos + i*k + j
			clause.clear();
			clause.push_back(Lit(fuzhupos + (i - 1) * num + j, true));
			clause.push_back(Lit(fuzhupos + (i)*num + j, false));
			(*solver).add_clause(clause);

		}

		//cons 7
		// index[i] -> ~(fuzhupos + (i-1)*k + k-1)
		clause.clear();
		clause.push_back(Lit(index[i], true));
		clause.push_back(Lit(fuzhupos + (i - 1) * num + k - 1, true));//s盒的数量是round*64 
		(*solver).add_clause(clause);
	}

	//cons 8
	// ~index[num-1] ∨ ~ fuzhupos + (num-2)*k + k-1
	clause.clear();
	clause.push_back(Lit(index[num - 1], true));
	clause.push_back(Lit(fuzhupos + (num - 2) * num + k - 1, true));//s盒的数量是round*64 
	(*solver).add_clause(clause);
}

void creat_sat_Sequence_encoding_new_weight(SATSolver* solver, int weight_old, int weight_new, int fuzhu_pos, int* index, int indexlen)
{
	// fuzhu_pos = 8000, index, indexlen = round * 48

	vector<Lit> clause;

	int fuzhupos = fuzhu_pos;
	//actually 48 sbox to estimate a link  96bool to weight weight_
	int num = indexlen; // round * 48
	//single state without combination

	
	//cons 2
	// 对前 weight 个辅助变量添加单文字子句，固定这些变量的值为 1。
	// ~(fuzhupos + i)	
	for (size_t i = weight_old; i < weight_new; i++)
	{
		clause.clear();
		clause.push_back(Lit(fuzhupos + i, true)); //与j是第几轮的有关系0-8对于0123
		(*solver).add_clause(clause);
	}

	int k = weight_new;
	for (size_t i = 1; i < num - 1; i++)
	{

		//just 2 or 3 is adequate
		for (size_t j = weight_old; j < k; j++)//这里应该取k，而不是k-1，这会导致最后一行没有递归
		{
			//cons 5
			//need to modify
			// 辅助变量与布尔变量之间建立联合关系，包含多个分支。
			// (index[i] ∨ fuzhupos + (i-1)*k + j-1) -> fuzhupos + i*k + j
			clause.clear();
			clause.push_back(Lit(index[i], true));
			clause.push_back(Lit(fuzhupos + (i - 1) * num + j - 1, true));
			clause.push_back(Lit(fuzhupos + (i)*num + j, false));
			(*solver).add_clause(clause);

			//cons 6
			// 递归辅助变量之间的直接关系。
			// fuzhupos + (i-1)*k + j -> fuzhupos + i*k + j
			clause.clear();
			clause.push_back(Lit(fuzhupos + (i - 1) * num + j, true));
			clause.push_back(Lit(fuzhupos + (i)*num + j, false));
			(*solver).add_clause(clause);

		}

		//cons 7
		// index[i] -> ~(fuzhupos + (i-1)*k + k-1)
		clause.clear();
		clause.push_back(Lit(index[i], true));
		clause.push_back(Lit(fuzhupos + (i - 1) * num + k - 1, true));//s盒的数量是round*64 
		(*solver).add_clause(clause);
	}

	//cons 8
	// ~index[num-1] ∨ ~ fuzhupos + (num-2)*k + k-1
	clause.clear();
	clause.push_back(Lit(index[num - 1], true));
	clause.push_back(Lit(fuzhupos + (num - 2) * num + k - 1, true));//s盒的数量是round*64 
	(*solver).add_clause(clause);
}


void add_constrain_linear_present(SATSolver* solver, int x_pos)
{
	// x_pos = r * 128 + 64
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
		//赋值字句，需要两条
		// 每一位的输出约束由两条子句表示：
		// 子句 1：~y_i ∨ x_j（y_i 是线性层的输出位，x_j 是输入位）
		clause.clear();
		
		clause.push_back(Lit(x_pos + 64 + i, true));
		clause.push_back(Lit(x_pos  + lineara[i], false));

		(*solver).add_clause(clause);

		// 子句 2：y_i ∨ ~x_j
		clause.clear();

		clause.push_back(Lit(x_pos + 64 + i, false));
		clause.push_back(Lit(x_pos + lineara[i], true));

		(*solver).add_clause(clause);
	}
}


void add_constrain_sbox_present_3weight(SATSolver* solver, int x_pos, int outx_pos, int consnum, int** cons8)
{
	// x_pos = r * 128, outx_pos = 5000 + r * 48, consnum = sbox_num_cons, cons8 = cnf_sbox

	vector<Lit> clause;
	/*clause.clear();
	(*solver).add_clause(clause);*/

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
				else if (cons8[j][k] == 2)
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
	// length=11, mid=5
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
			/*     for (size_t i = 0; i < expr.size(); i++) {
					 char c = expr[i];
					 if (c >= 'x' && c <= 'x' + 7) {
						 expression[c - 'x'] = 1;
					 }
					 else if (c >= 'y' && c <= 'y' + 7) {
						 expression[8 + c - 'y'] = 1;
					 }
					 else if (c >= 'x' + 8 && c <= 'x' + 15) {
						 expression[c - 'x' - 8] = 2;
					 }
					 else if (c >= 'y' + 8 && c <= 'y' + 15) {
						 expression[8 + c - 'y' - 8] = 2;
					 }
				 }*/
				 //这一段需要修改
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

				if (c == 'x'){
					if (expr[i + 2] == '\'')//表示反
					{
						expression[expr[i + 1] - '1'] = 2;
					}else{
						expression[expr[i + 1] - '1'] = 1;
					}
				}
				if (c == 'y'){
					if (expr[i + 2] == '\'')//表示反
					{
						expression[mid + expr[i + 1] - '1'] = 2;
					}else{
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
	//  printf("\n  TRACE: %d \n", count); //trace

	for (int i = 0; i < count; i++)
	{
		for (size_t j = 0; j < cnf_length; j++)
		{
			result[i][j] = expressions[i][j];
		}
	}
	//for (int i = 0; i < count; i++)
	//{
	//    for (size_t j = 0; j < 10; j++)
	//    {
	//        constrain_ascon10[i][j] = expressions[i][j];
	//    }
	//}

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
	// 存储约束，假设最大 150 条约束，每条约束最多包含 20 个整数。
	int** cnf_sbox;
	cnf_sbox = (int**)malloc(150 * 8); //100条约束
	for (size_t i = 0; i < 150; i++)
	{
		cnf_sbox[i] = (int*)malloc(20 * 8);//一般约束比较短
	}
	int sbox_num_cons;
	
	read_constrant_from_file("present_ddt_3one.txt", &sbox_num_cons, cnf_sbox, 11, 5);
	//基于贪心算法的约减工作可以放在后续一起做


	SATSolver solver;
	solver.set_num_threads(thread_); 	//4
	solver.new_vars(round_ * (10000));	//每轮分配10000个变量
	// TRACE
	solver.log_to_file("r15.log");


	int round = round_;
	vector<Lit> clause;

	clause.clear();
	// 初始约束：添加一个包含 64 个文字的子句（假设为初始状态的固定值）。
	for (size_t i = 0; i < 64; i++)
	{
		clause.push_back(Lit(i, false));
	}
	solver.add_clause(clause);

	// 添加与 S-Box 、线性层的约束。
	for (size_t i = 0; i < round_; i++)
	{
		add_constrain_sbox_present_3weight(&solver, i * 128, 5000 + i * 48, sbox_num_cons, cnf_sbox);
		add_constrain_linear_present(&solver, i * 128 + 64);
	}

	int index[3000] = { 0 };	//定义序列编码中使用的变量索引
	for (size_t i = 0; i < round * 3 * 16; i++)
	{
		index[i] = 5000 + i;
	}
	creat_sat_Sequence_encoding_increment_method(&solver, weight, 8000, index, round * 48);

	//	printf("\n 开始求解\n");
	clock_t t1, t2;

	int k = weight;

	lbool ret = l_False;

	while( ret == l_False ){
		t1 = clock();
	
		ret = solver.solve();

		t2 = clock();
		std::ostringstream oss;
		double elapsed_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
		oss << std::fixed << std::setprecision(2) << elapsed_time;  // 设置小数点后保留两位
		std::string time_string = oss.str();


		if (ret == l_False)
		{
			printf("\n 序列 该子模型无解！ %d %d\n", round_, weight);
			std::string filename = "present_incremental_r15.csv";
			std::ofstream file(filename, std::ios_base::app);
			//轮数 重量，时间 有无解 0无解 1有解
			file << round_ << "," << weight << "," << time_string <<","<< 0 << "\n";  // 以逗号分隔每列，换行结束每行
			file.close();
			creat_sat_Sequence_encoding_new_weight(&solver, k, k+1, 8000, index, round * 48);
			k++;
		}
		else if (ret == l_True)
		{
			printf("\n 序列 该模型有解！ %d %d\n", round_, weight);
			std::string filename = "present_incremental_r15.csv";
			std::ofstream file(filename, std::ios_base::app);
			//轮数 重量，时间 有无解 0无解 1有解
			file << round_ << "," << weight << "," << time_string << "," << 1 << "\n";  // 以逗号分隔每列，换行结束每行
			file.close();
			solver.print_stats();
			break;
		}		
	}
	return true;
}


void main_control_center_present(int round)
{
	bool b1 = false;
	clock_t t1, t2;
	t1 = clock();
	// int weight = 2;
	// while (!b1)
	// {
	// 	//序列编码
	// 	b1 = creat_sat_present_differential_xuliebianma(round, 4, 0, weight++);
	// 	// StaticTimer256ublock::ctimer();
	// }
	// b1 = false;
	int weight = 10;
	b1 = creat_sat_present_differential_xuliebianma(round, 4, 0, weight);
	

	t2 = clock();

	 std::ostringstream oss;
	 double elapsed_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
	 oss << std::fixed << std::setprecision(2) << elapsed_time;  // 设置小数点后保留两位
	 std::string time_string = oss.str();

	std::string filename = "present_incremental_r15.csv";
 	std::ofstream file(filename, std::ios_base::app);

 	//轮数 重量，时间 有无解 0无解 1有解
 	file << round << "," << weight << "," << time_string << "\n";  // 以逗号分隔每列，换行结束每行

 	file.close();

}


int main()
{
//    for (size_t round = 2; round < 30; round++)
//    {
//        main_control_center_present(round);
//    }
	int round=15;
	main_control_center_present(round);
	
}

// // g++ main.cpp -o main_incre -lcryptominisat5