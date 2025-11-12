#include <iostream>
#include <thread>
#include <vector>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <cstdlib>
#include <cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <time.h>
#include <fstream>
#include <string>
#include <cassert>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>  // 对于 std::put_time
#include <climits>

#include <unordered_set>
#include <csignal>
#include <unordered_map>

using namespace std;
using std::vector;
using namespace CMSat;

#define R 20


atomic<int> globalMinSAT(INT_MAX);
atomic<int> globalMaxUNSAT(INT_MIN);
atomic<bool> stop(false); // 全局停止标志
queue<pair<int, int>> taskQueue;
mutex queueMutex;
condition_variable cv;
vector<thread> threads; 
// unordered_set<int> runningThreads; // 记录运行中的线程 ID
unordered_map<int, int> threadCurrentTask; // 记录线程当前的 m 值

mutex threadSetMutex; // 保护 runningThreads 的互斥锁

mutex clauseMutex;
condition_variable clauseCV;
// vector<vector<int>> globalClauses;  // 存储所有需要传播的 UNSAT 子句
int maxUnsatM;
bool newClauseAvailable = false;    // 标记是否有新的子句可用
// atomic<bool> newClauseAvailable(false);  // 是否有新子句


void creat_sat_Sequence_encoding_given_original_method(SATSolver* solver, int weight, int fuzhu_pos, int* index, int indexlen)
{
    // sum(x) <= w
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

void creat_sat_Sequence_encoding_neg(SATSolver* solver, int weight, int neg_fuzhu_pos, int* index, int indexlen)
{
    // sum(x) > m --> sum(x) >= m + 1 --> sum(-x) <= n - m - 1

	// fuzhu_pos = 8000, index, indexlen = round * 48
    // neg_fuzhu_pos = 200000
	vector<Lit> clause;

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

	// R_n-1,m-1 
    // at-least-m

	clause.push_back(Lit(fuzhu_pos + (num - 1) * k + m -1, false)); 

	(*solver).add_clause(clause);
}

void add_constrain_linear_present(SATSolver* solver, int x_pos){
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

void add_constrain_sbox_present_3weight(SATSolver* solver, int x_pos, int outx_pos, int consnum, int** cons8){
	// x_pos = r * 128, outx_pos = 5000 + r * 48, consnum = sbox_num_cons, cons8 = cnf_sbox
	vector<Lit> clause;
	/*clause.clear();
	(*solver).add_clause(clause);*/

	//16个s盒
	for (size_t i = 0; i < 16; i++){
		//32个s盒每轮 我们的目的是找到一条解就行了，不需要全部找到
		for (size_t j = 0; j < consnum; j++)//每个s盒包含如此之多的约束
		{
			//1表示有 2表示反
			clause.clear();

			//前16个s盒的结果存ori 以及变量
			for (size_t k = 0; k < 4; k++){
				//一个是向前的bit，一个是向后的bit,前面是原来的，要减去整个16*64，当然大家都会减去一个exchange
				if (cons8[j][k] == 1){
					clause.push_back(Lit(x_pos + i + k * 16, false));//与j是第几轮的有关系0-8对于0123
				}else if (cons8[j][k] == 2){
					clause.push_back(Lit(x_pos + i + k * 16, true));//与j是第几轮的有关系0-8对于0123
				}

				//变换的结果存储在临时变量里面，索引照用就行，实际上这个临时变量只使用一半
				//可以考虑临时变量后半段是ori的变换结果 前半段是multi8的变换结果
				if (cons8[j][k + 4] == 1){
					clause.push_back(Lit(x_pos + i + k * 16 + 64, false));//与j是第几轮的有关系0-8对于0123
				}else  if (cons8[j][k + 4] == 2){
					clause.push_back(Lit(x_pos + i + k * 16 + 64, true));//与j是第几轮的有关系0-8对于0123
				}
			}

			for (size_t k = 0; k < 3; k++){
				//差分概率
				if (cons8[j][k + 8] == 1){
					clause.push_back(Lit(outx_pos + i * 3 + k, false));//与j是第几轮的有关系0-8对于0123
				}else  if (cons8[j][k + 8] == 2){
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
			
			size_t length = expr.length();
			size_t i = 0;
			while (i < length) {
				char c = expr[i];
				
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

void read_constrant_from_file(std::string filename,int* num,int** result,int cnf_length,int mid){
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
	}
	*num = count;
	//  printf("\n  TRACE: %d \n", count); //trace

	for (int i = 0; i < count; i++){
		for (size_t j = 0; j < cnf_length; j++){
			result[i][j] = expressions[i][j];
		}
	}

	for (int i = 0; i < count; i++){
		//    printf("\n");
		for (size_t j = 0; j < 10; j++){
			//    printf(" %d ", constrain_ascon10[i][j]);
		}
	}
}


// void clauseListenerFunc() {
//     while (solving.load()) {
//         unique_lock<mutex> lk(clauseMutex);
//         clauseCV.wait(lk, [&]() { return newClauseAvailable || !solving.load(); });

//         if (!solving.load()) break;  // 退出监听线程
//         cout << "[ClauseListener] Detected new clause. Notifying solve() thread." << endl;

//         // 只通知solve()线程，不直接添加子句
//         clauseCV.notify_all();
//     }
// }

bool creat_sat_present_differential_xuliebianma(int round_, int thread_, int out_time, int m, int weight)
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

	SATSolver solver;
	solver.set_num_threads(thread_); 	//4
	solver.new_vars(round_ * (10000));	//每轮分配10000个变量
	// TRACE
	// solver.log_to_file("r20.log");

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
	creat_sat_Sequence_encoding_given_original_method(&solver, weight, 8000, index, round * 48);
    at_least_m(&solver, weight, m, 8000,  round * 48);

	//	printf("\n 开始求解\n");
	clock_t t1, t2;

	t1 = clock();
	
	// lbool ret = solver.solve();
	atomic<bool> solving(true);  // 记录求解器是否正在运行
	// bool hasPreviousClause = false;  // 是否已经添加过子句

	/*
	// **在求解过程中，动态监听新子句**
    thread clauseListener([&]() {
        
        while (solving.load()) {  // 只要求解器还在运行，就持续监听
			unique_lock<mutex> lk(clauseMutex);
            clauseCV.wait(lk, [&]() { return newClauseAvailable || !solving.load(); });  // 等待新子句

			if (!solving.load()) break;  // 如果求解器已经结束，退出监听线程
            
			// creat_sat_Sequence_encoding_neg(&solver, weight, 200000, index, round * 48);
            solver.interrupt_asap();

			// at_least_m(&solver, weight, maxUnsatM, 8000,  round * 48);
            cout << "[ClauseListener] New clause added via at_least_m. maxUnsatM = " << maxUnsatM << endl; //trace

            // print learned_clause

			vector<Lit> assumptions;
			assumptions.push_back(Lit(8000+(round*48-1)*weight+maxUnsatM-1, false));

            

			lbool result = solver.solve(&assumptions);
            cout << "[ClauseListener] Re-solving triggered. Result = " << result << endl; //trace

            newClauseAvailable = false;  // 处理完后重置标志
        }
    });

	*/

    // if (newClauseAvailable.load(std::memory_order_relaxed)) {
    //     std::lock_guard<std::mutex> lock(queueMutex);

    //     at_least_m(&solver, weight, maxUnsatM, 8000,  round * 48);
    //     lbool result = solver.solve();
            
    //     newClauseAvailable.store(false, std::memory_order_relaxed);  // 清除标记
    // }

	lbool ret;
	{
		solving.store(true);  // 设置标志，表示正在求解
        cout << "[Main] Starting main solver.solve()..." << endl;
		ret = solver.solve();
        cout << "[Main] Main solver.solve() finished with result: " << ret << endl;
		solving.store(false); // 结束后设置为 false，通知监听线程退出
	}

	// // 确保监听线程正确退出
	// clauseCV.notify_all();  
	// clauseListener.join();  //join() 已经结束线程

	// 终止监听线程
    // clauseListener.detach();


	t2 = clock();
	std::ostringstream oss;
	double elapsed_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
	oss << std::fixed << std::setprecision(2) << elapsed_time;  // 设置小数点后保留两位
	std::string time_string = oss.str();


	if (ret == l_False)
	{
        cout << "[SAT] No solution found for round " << round_ << ", weight " << weight << ", cost " << time_string << endl;
        
		printf("\n 序列 该子模型无解！ %d %d\n", round_, weight);

		// std::string filename = "sequnence_code_present.csv";
		// std::ofstream file(filename, std::ios_base::app);

		// //轮数 重量，时间 有无解 0无解 1有解
		// file << round_ << "," << weight << "," << time_string <<","<< 0 << "\n";  // 以逗号分隔每列，换行结束每行

		// file.close();

		return false;
	}
	else if (ret == l_True)
	{
        cout << "[SAT] A solution was found for round " << round_ << ", weight " << weight <<", cost " << time_string << endl;       
		printf("\n 序列 该模型有解！ %d %d\n", round_, weight);

		// std::string filename = "sequnence_code_present.csv";
		// std::ofstream file(filename, std::ios_base::app);

		// //轮数 重量，时间 有无解 0无解 1有解
		// file << round_ << "," << weight << "," << time_string << "," << 1 << "\n";  // 以逗号分隔每列，换行结束每行
		// file.close();

		// TRACE
		// for(int i = 0; i < n; i++){
        //     std::cout<<solver.get_model()[i]<<", ";
        // }
        // std::cout<< std::endl;
        // solver.print_stats();

		return true;
	}
	return true;
}

// Function to simulate SAT solving (replace with actual solver call)
bool solveSATInstance(int m, int k) {
    // thread=4
    bool b1 = creat_sat_present_differential_xuliebianma(R, 4, 0, m, k);
    return b1;
}

// Worker thread function
void workerThread(int threadID) {

	// 持续运行直到全局停止标志stop被置为true
    while (!stop.load()) {
        // int m;
		pair<int, int> task;

        // Fetch task from queue
        {
            unique_lock<mutex> lock(queueMutex);
			// 线程休眠直到任务队列非空或收到停止信号
            cv.wait(lock, [&]() { return !taskQueue.empty() || stop.load(); });

            if (stop.load()) break;

            task = taskQueue.front();
            taskQueue.pop();
        }
		int left = task.first, m = task.second;

		{
            lock_guard<mutex> lock(threadSetMutex);
            threadCurrentTask[threadID] = m; // 记录当前线程计算的 m
        }

        cout << "[Worker " << threadID << "] Processing m = " << left <<" Processing k = " << m 
             << " for interval [" << left << ", " << m << "]" << endl;

        // Solve SAT instance
        bool isSAT = solveSATInstance(left, m);

        cout << "[Worker " << threadID << "] Result for m = " << m 
             << " isSAT: " << isSAT << endl;

        if (isSAT) {
            // Update globalMinSAT
            int prevMinSAT = globalMinSAT.load();
            while (m < prevMinSAT && !globalMinSAT.compare_exchange_weak(prevMinSAT, m));

			// todo：剪枝，终止所有 > 当前m的任务线程（因为这些m值不是最优的）
			{
                lock_guard<mutex> lock(threadSetMutex);
				for (auto it = threadCurrentTask.begin(); it != threadCurrentTask.end();) {
					int threadID = it->first;
					int taskM = it->second;

					if (taskM > m) {
                        cout << "[Worker " << threadID << "] Canceling thread " << threadID 
                             << " with m = " << taskM << " (current SAT m = " << m << ")" << endl;
                        
						pthread_cancel(threads[threadID].native_handle()); // 强制终止线程
						it = threadCurrentTask.erase(it); // 从记录中移除
					} else {
						++it;
					}
				}
            }

            // Add new tasks with refined intervals
            {
                lock_guard<mutex> lock(queueMutex);
                // if (m - 1 > globalMaxUNSAT.load()) {
				if(left < m - 1){ //区间未收缩到单值
                    int mid = (m + left) / 2;
                    taskQueue.push({left, mid});
                    cout << "[Worker " << threadID << "] Pushed new task: [" << left << ", " << mid << "]" << endl;
                    
					cv.notify_one();
					taskQueue.push({mid, m-1});
                    cout << "[Worker " << threadID << "] Pushed new task: [" << mid << ", " << m - 1 << "]" << endl;
                    
                    cv.notify_one();	//条件通知	唤醒一个线程处理新任务
                }
            }

			// 添加m-1的线程
			

        } else {
            // Update globalMaxUNSAT
            int prevMaxUNSAT = globalMaxUNSAT.load();
            while (m > prevMaxUNSAT && !globalMaxUNSAT.compare_exchange_weak(prevMaxUNSAT, m));
			
			// todo：剪枝，终止所有 < 当前m的任务线程（因为这些m值已被证明不可行）
            {
                lock_guard<mutex> lock(threadSetMutex);
				for (auto it = threadCurrentTask.begin(); it != threadCurrentTask.end();) {
					int threadID = it->first;
					int taskM = it->second;

					if (taskM < m) {
                        cout << "[Worker " << threadID << "] Canceling thread " << threadID 
                             << " with m = " << taskM << " (current UNSAT m = " << m << ")" << endl;
                        
						pthread_cancel(threads[threadID].native_handle()); // 强制终止线程
						it = threadCurrentTask.erase(it); // 从记录中移除
					} else {
						++it;
						// 传递 w > m的信息
						maxUnsatM = m;
						newClauseAvailable = true;  // 标记新子句可用
                        // newClauseAvailable.store(true, std::memory_order_relaxed);  // 标记有新子句
    					clauseCV.notify_all();  // 通知所有等待的 SAT 线程
					}
				}
            }
        }

		{
            lock_guard<mutex> lock(threadSetMutex);
            threadCurrentTask.erase(threadID); // 任务完成后删除记录
        }
        break;
    }
}                  

int main() {
    // int r=20;
    const int initialMValues[] = {70, 80, 90, 100}; // Initial intervals: example values
    const int numThreads = 8; // Number of threads

    // queue<int,int> taskQueue; // Task queue holding m values
    // mutex queueMutex; // Mutex for task queue 互斥锁
    
    // Initialize task queue
    {
        lock_guard<mutex> lock(queueMutex);
        // 初始化任务队列
        taskQueue.push({1, initialMValues[0]});
        for (int i = 0; i < int(sizeof(initialMValues) / sizeof(initialMValues[0])) - 1; i++) {
            taskQueue.push({initialMValues[i], initialMValues[i+1]});
        }
    }
    // lock_guard<mutex> lock(queueMutex);	//lock_guard确保线程安全的队列操作
	// taskQueue.push({1,initialMValues[0]});
	// for(int i=0; i< sizeof(initialMValues) / sizeof(initialMValues[0])-1; i++){
	// 	taskQueue.push({initialMValues[i], initialMValues[i+1]});
	// }
    
    // Start worker threads
	for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back(workerThread, i);
    }

    // 主线程监控搜索状态，当 globalMinSAT 与 globalMaxUNSAT 收敛时结束
    while (!stop.load()) {
        // 输出当前状态，便于追踪
        // if(globalMinSAT<=initialMValues[3]||globalMaxUNSAT>=initialMValues[0]){
        //     cout << "[Main] Monitoring: globalMinSAT = " << globalMinSAT 
        //         << ", globalMaxUNSAT = " << globalMaxUNSAT << endl;
        // }
        if (globalMinSAT.load() <= globalMaxUNSAT.load() + 1) {
            cout << "[Main] Termination condition met: globalMinSAT <= globalMaxUNSAT + 1" << endl;
            stop.store(true);
            cv.notify_all();  // 通知所有等待的线程退出
            break;
        }
        // 每 100 毫秒检查一次
        this_thread::sleep_for(chrono::milliseconds(1000));
    }

    // Wait for all threads to finish
    for (auto& th : threads) {
        th.join();	// 等待所有线程完成
    }

    // Output the result
    if (globalMinSAT < INT_MAX) {
        cout << "Maximum differential probability (minimum m) is: " << globalMinSAT << endl;
    } else {
        cout << "No SAT solution found." << endl;
    }

    return 0;
}

// g++ parallel.cpp -o parallel -lcryptominisat5 -pthread
// ./parallel