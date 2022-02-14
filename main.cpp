#include <cctype> 
#include <cerrno> 
#include <cfloat> 
#include <ciso646> 
#include <climits> 
#include <clocale> 
#include <cmath> 
#include <ctime> 
  
// C++ 
#include <algorithm> 
#include <bitset> 
#include <complex> 
#include <deque> 
#include <exception> 
#include <fstream> 
#include <functional> 
#include <iomanip> 
#include <ios> 
#include <iosfwd> 
#include <iostream> 
#include <istream> 
#include <iterator> 
#include <limits> 
#include <list> 
#include <locale> 
#include <map> 
#include <memory> 
#include <new> 
#include <numeric> 
#include <ostream> 
#include <queue> 
#include <set> 
#include <sstream> 
#include <stack> 
#include <stdexcept> 
#include <streambuf> 
#include <string> 
#include <typeinfo> 
#include <utility> 
#include <valarray> 
#include <vector> 
  
#if __cplusplus >= 201103L 
#include <array> 
#include <atomic> 
#include <chrono> 
#include <condition_variable> 
#include <forward_list> 
#include <future> 
#include <initializer_list> 
#include <mutex> 
#include <random> 
#include <ratio> 
#include <regex> 
#include <scoped_allocator> 
#include <system_error> 
#include <thread> 
#include <tuple> 
#include <typeindex> 
#include <type_traits> 
#include <unordered_map> 
#include <unordered_set> 
#endif 
#include "system_model.h"
using namespace std;

int NO_OF_CORES = 8;

int NO_OF_TASKS = 50;

double DEADLINE = 0;

double LIFETIME_THRESHOLD = 40;

Multicore* multicore;


//**********calc start**************************************************************************************************88

double taskReliability(task* t, Core* c){
	double ro = c->ro;
	double w = c->hardware_coefficient;
	double Fmin = c->Fmin;
	double Fmax = c->Fmax;

	double wc = t->worst_case_time;
	double freq = t->freq_assigned;

	double rf = ro*(pow(10,(w*(1-freq))/(1-Fmin)));
	double Ri = exp((-1)*(rf*wc)/freq);

	//cout<<" task id "<<t->id<<" Ri : "<<Ri<<" \n";

	return Ri;
}

double taskRecReliability(task* t, Core* c){
	double Ri = taskReliability(t,c);
	double Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);
	return Rreci;
}

double taskReliability_dynamic(task* t, Core* t, double execution_time){
	double ro = c->ro;
	double w = c->hardware_coefficient;
	double Fmin = c->Fmin;
	double Fmax = c->Fmax;

	double freq = t->freq_assigned;

	double rf = ro*(pow(10,(w*(1-freq))/(1-Fmin)));
	double Ri = exp((-1)*(rf*execution_time));

	//cout<<" task id "<<t->id<<" Ri : "<<Ri<<" \n";

	return Ri;
}

// schedule: vector with task pointers in a sorted order as per priority
// hasExecuted : bit mask to represent fail/success for a task
// in: the bit index to be fixed in this function call
// i: the index of task in schedule for which we are calculating the Prec
void PrecHelper(Core* core, int hasExecuted,int in,int i,double* PMil){
	if(in==i){
		//calculate PMil
		double slack = DEADLINE - core->total_time;

		for(int j=0;j<i;j++){
			if((hasExecuted & (1 <<j)) == 0){
				double time_rec = core->tasks[j]->worst_case_time/core->tasks[j]->freq_assigned;
				if(slack >= time_rec) slack-=time_rec;
			}
		}
		double time_rec = core->tasks[i]->worst_case_time/core->tasks[i]->freq_assigned;
		if(slack >= time_rec){

			double prob_of_mode=1;
			for(int j=0;j<i;j++){

				double ri = taskReliability(core->tasks[j],core);
				if((hasExecuted & (1 << j)) == 0){
					prob_of_mode*=(1-ri);
				}
				else prob_of_mode*=(ri);
			}

			*PMil+=prob_of_mode;
			//cout<<"prob_of_mode : "<<prob_of_mode<<" ";
		}

		return;
	}

	//two cases for bit mask : in'th task fails or succeeds:
	PrecHelper(core,hasExecuted,in+1,i,PMil);
	PrecHelper(core,hasExecuted|(1 << in) ,in+1,i,PMil);
}

double taskPreci(Core* c,int i){
	double freq_assigned = c->tasks[i]->freq_assigned;
	c->tasks[i]->freq_assigned = 0.9; //fmax
	double RiFmax = taskReliability(c->tasks[i],c);
	c->tasks[i]->freq_assigned = freq_assigned; //restoring old value;

	double pMil = 0;
	double *PMil = &pMil;

	*PMil = 0;
	PrecHelper(c,0,0,i,PMil);

	// cout<<" Task id "<<c->tasks[i]->id<<" Pmil: "<<*PMil<<" * RiFmax : "<<RiFmax<<"\n";
	return RiFmax*(*PMil);
}

double systemReliability(Multicore* multicore){
	double Rsys = 1;

	for(auto core : multicore->cores){
		int i=0;
		for(auto task : core->tasks){

			double Ri = taskReliability(task,core);
			// double Preci = taskPreci(core,i);
			double Rreci = Ri;
			if(task->recovery_assigned){
				Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);
			}

			Rsys*=Rreci;

			// cout<<"Task id : "<<task->id<<" Ri :"<<Ri<<" Preci :"<<Preci<<" Rreci : "<<Rreci<<"Rsys : "<<Rsys<<" \n";

			i++;
		}
	}

	return Rsys;
	
}

double lifetimeReliability(task* t,Core* c){
	double constant_a = 38.92;
	double freq = t->freq_assigned;
	double MTTF = constant_a / (freq-0.1); //freq in ghz please (10^9) 
	//TO-DO(niyatic): multiply by root Ao.

	cout<<"freq "<<freq<<" mttf: "<<MTTF<<" \n";
	return MTTF;
	
}

double systemLifetimeReliability(Multicore* multicore){
	double MTTFsys = DBL_MAX;
	for(auto core : multicore->cores){
		for(auto task : core->tasks){
			double MTTFi = lifetimeReliability(task,core);
			MTTFsys = min(MTTFsys,MTTFi);
		}
	}
	return MTTFsys;
}

//**********calc end**************************************************************************************************88


//**********freq assigning algos**************************************************************************************************88

void assign_rand_freq(DAG* dag){
	for(task* t : dag->nodes){
		t->freq_assigned = ((double)rand()/RAND_MAX) *(0.9 - 0.3) + 0.3;
	}
}

void assign_min_freq(DAG* dag){
	for(task* t : dag->nodes){
		t->freq_assigned = 0.3;
	}
}

void assign_max_freq(DAG* dag){
	for(task* t : dag->nodes){
		t->freq_assigned = 0.9;
	}
}

void assign_freq(DAG* dag,int mode){
	switch(mode){
		case 0 : assign_min_freq(dag); break;
		case 1 : assign_max_freq(dag); break;
		case 2 : assign_rand_freq(dag); break;
		default : break;
	}
}

//**********freq assigning algos end**************************************************************************************************88

//**********recovery assigning algos**************************************************************************************************88
//****sem8
int numberOfRecoveriesStatic(DAG* dag, int start, int end){
	if(start<=end) return start;
	if(start == end - 1){
		DAG* new_dag = DAG* new_dag = getNewDAG(dag, end);
	    bool can_schedule_new_dag = canScheduleStatic(new_dag);
		if(can_schedule_new_dag) return end;
		else return start;
	}
	int mid = start + (end-start)/2;
	DAG* new_dag = getNewDAG(dag, mid);
	bool can_schedule_new_dag = canScheduleStatic(new_dag);
	if(can_schedule_new_dag){
		numberOfRecoveriesStatic(dag, mid, end);
	}
	else numberOfRecoveriesStatic(dag, start, mid-1);
}

bool compWorstCaseTime(task* i,task* j){
	return i->worst_case_time < j->worst_case_time;
}
bool compEndTime(task* i,task* j){
	return i->end_time < j->end_time;
}
bool compVexitDistance(task* i,task* j){
	return i->vexit_dist < j->vexit_dist;
}

void topologicalSortUtil(DAG* dag,int v, bool visited[],stack<int>& stack)
{
	  // Mark the current node as visited
	  visited[v] = true;
	  task* nodev = dag->nodes[v];
	  // Recur for all the vertices adjacent to this vertex
	  for (auto pred:nodev->predecessors) {
	    if (!visited[pred->id])
	      topologicalSortUtil(dag,pred->id, visited, stack);
	  }

	  // Push current vertex to stack which stores topological
	  // sort
	  stack.push(v);
}

void longestDistVexit(DAG* dag, int vexit_id)
{
    stack<int> Stack;
   	
   	int V = vexit_id+1;
    // Mark all the vertices as not visited
    bool* visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;
   
    // Call the recursive helper function to store Topological
    // Sort starting from all vertices one by one
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            topologicalSortUtil(i, visited, Stack);
   
    // Initialize distances to all vertices as infinite and
    // distance to source as 0
    for (auto node:dag->nodes)
        node->vexit_dist = NINF;
    dag->nodes[vexit_id] = 0;
    // Process vertices in topological order
    while (Stack.empty() == false) {
        // Get the next vertex from topological order
        int u = Stack.top();
        Stack.pop();
        task* nodeu = dag->nodes[u];
   
        if (nodeu->vexit_dist != NINF) {
            for (auto pred: nodeu->predecessors){
             
                if (pred->vexit_dist < nodeu->vexit_dist + 1)
                    pred->vexit_dist = nodeu->vexit_dist + 1;
            }
        }
    }
     
    delete [] visited;
}

//returns a new dag with 
// 1) duplicated nodes for tasks with recovery as per criteria
// 2) added dummy nodes vexit and ventry
// 3) dist from vexit is updated by calling a function
DAG* getNewDAG(DAG* dag, int number_of_recoveries){
	DAG* new_dag = new DAG();
	*new_dag = *dag;
	//TODO(krishnahere): check ^
	sort(new_dag->nodes.begin(), new_dag->nodes.end(),compWorstCaseTime);
	for (int i = 0; i < number_of_recoveries; i++){
		task* node = nodes[i];
		node->recovery_assigned=true;
		task* new_node = new task(); //successor - prede.. done
		*new_node = *node;
		//TODO(krishnahere): check ^
		new_dag->nodes.push_back(new_node);
		for(auto pred : node->predecessors){
			pred->successors.push_back(new_node);
		}
		for(auto suc :  node->sucessors){
			suc->predecessors.push_back(new_node);
		}
	}

	task* ventry = new task();
	ventry->worst_case_time = 0;

	for(auto node:new_dag->nodes){
		if(node->predecessors.empty()){
			ventry->successors.push_back(node);
			node->predecessors.push_back(ventry);
		}
	}

	task* vexit = new task();
	vexit->worst_case_time = 0;
	for(auto node: new_dag->nodes){
		if(node->successors.empty()){
			vexit->predecessors.push_back(node);
			node->successors.push_back(vexit);
		}
	}

	ventry->id = new_dag->nodes.size();
	vexit->id = ventry->id + 1;
	new_dag->nodes.push_back(ventry);
	new_dag->nodes.push_back(vexit);

	longestDistVexit(dag,vexit->id);

	return new_dag;
}

bool canScheduleStatic(DAG* node_graph){
	//TODO(krishnahere): DO all.
	
	vector<Core*> free_cores;
	priority_queue<task*, vector<task*>,compEndTime> processing_queue;
	priority_queue<task*, vector<task*>, compVexitDistance> ready_queue;


}





//****sem8
void randomAlgo(Multicore* multicore, double deadline){
	
	for(auto core : multicore->cores){

		random_shuffle(core->tasks.begin(),core->tasks.end());
		double slack = deadline - core->total_time;

		for(auto l:core->tasks){
			if(l->worst_case_time > slack) break;
			slack-=l->worst_case_time;
			l->recovery_assigned = true;
		}
	}

}

bool comp_greedyReliability(task* i,task* j){
	bool iToj = false, jToi = false;
	for(auto l:i->predecessors){
		if(l==j) jToi = true;
	}
	for(auto l:j->predecessors){
		if(l==i) iToj = true;
	}

	if(jToi) return false;
	if(iToj) return true;

	return taskReliability(i,multicore->cores[i->core_assigned]) > taskReliability(j,multicore->cores[j->core_assigned]);
}

void greedyReliabilityAlgo(Multicore* multicore, double deadline){
	//assign priority as per the dependency and reliability

	for(auto core : multicore->cores){
		//tasks vector to be sorted
		sort(core->tasks.begin(),core->tasks.end(),comp_greedyReliability);
		double slack = deadline - core->total_time;

		cout<<"sorted : ";
		for(auto l:core->tasks){
			cout<<l->id<<"  ";
		}
		cout<<endl;

		for(auto l:core->tasks){
			if(l->worst_case_time > slack) break;
			slack-=l->worst_case_time;
			l->recovery_assigned = true;
		}
	}

}

//paper 4 algo :
bool comp_greedyTime(task* i,task* j){
	bool iToj = false, jToi = false;
	for(auto l:i->predecessors){
		if(l==j) jToi = true;
	}
	for(auto l:j->predecessors){
		if(l==i) iToj = true;
	}

	if(jToi) return false;
	if(iToj) return true;

	return i->worst_case_time/i->freq_assigned < j->worst_case_time/j->freq_assigned;
}

void greedyTimeAlgo(Multicore* multicore, double deadline){
	//assign priority as per the dependency and reliability

	for(auto core : multicore->cores){
		//tasks vector to be sorted
		sort(core->tasks.begin(),core->tasks.end(),comp_greedyTime);
		double slack = deadline - core->total_time;

		cout<<"sorted : ";
		for(auto l:core->tasks){
			cout<<l->id<<"  ";
		}
		cout<<endl;

		for(auto l:core->tasks){
			if(l->worst_case_time > slack) break;
			slack-=l->worst_case_time;
			l->recovery_assigned = true;
		}
	}

}

bool isFeasible(task* i, task* j){
	bool jToi = false;
	bool iToj = false;

	for(auto l:i->predecessors){
		if(l==j) jToi = true;
	}
	for(auto l:j->predecessors){
		if(l==i) iToj = true;
	}

	return iToj;
}

// void simulatedAnnealingAlgo(Multicore* multicore, double deadline){
	
// 	greedyTimeAlgo(multicore,deadline);

// 	for(auto core : multicore->cores){
// 		vector<task*> next_state;
// 		vector<task*> curr_state = core->tasks;
// 		int si = curr_state.size();

// 		int base_index = rand()%(si-1 + 0 - 1) + 0;
// 		int movable_index = rand()%(si-1 + 0 - 1) + 0;

// 		while(base_index == movable_index && !isFeasible(curr_state[base_index],curr_state[movable_index])){
// 			movable_index = rand()%(si-1 + 0 - 1) + 0;
// 		} 

// 		next_state = curr_state;

// 		task* tempTask = next_state[base_index];
// 		next_state[base_index] = next_state[movable_index];
// 		next_state[movable_index] = tempTask;


// 		double w = systemReliability(multicore);
// 		core->task = next_state;
// 		double w_new = systemReliability(multicore);
// 		if(w_new > w){
// 			core->task = next_state;
// 		}
// 	}

// }

int no_of_recoveries(){
	int low = 0;
	int high = dag->nodes.size();
	int mid;
	while(high>low)
	{
		mid = low + (high - low)/2;
		if(static_schedule()){
			low = mid;
		}
		else{
			high = mid - 1;
		}
	}
	return mid;
}

//**********recovery assigning algos end**************************************************************************************************88
 
//**********constraint checks**************************************************************************************************88

bool isWithinDeadline(Multicore* multicore, double deadline){

	for(auto core : multicore->cores){

		double time_taken = 0;

		for(auto task : core->tasks){
			double ti = 0;
			if(task->recovery_assigned) ti = 2*task->worst_case_time;
			else ti = task->worst_case_time;
			time_taken+=ti;
		}

		if(time_taken > deadline) return false;
	}

	return true;
}

bool hasMTTF(Multicore* multicore){

	return systemReliability(multicore) >= LIFETIME_THRESHOLD;
}

//**********constraint checks end**************************************************************************************************88

int main(){

	// cout<<"generating cores\n";
	// //generating cores
	// multicore = new Multicore();
	// multicore->generateCores(NO_OF_CORES);


	// cout<<"generating DAG\n";
	// //create or input a DAG
	// DAG* dag = new DAG();
	// dag->generateDAG(NO_OF_TASKS,NO_OF_CORES);

	// for(int i=0;i<dag->nodes.size();i++){
	// 	int core_id = dag->nodes[i]->core_assigned;
	// 	multicore->cores[core_id]->tasks.push_back(dag->nodes[i]);
	// }


	// //dag->displayDAG();

	// //applying algo to assign frequencies to tasks

	// assign_freq(dag,2);

	double max_time_taken = 0;
	// for(int i=0;i<dag->nodes.size();i++){
	// 	int core_id = dag->nodes[i]->core_assigned;
	// 	multicore->cores[core_id]->total_time+=dag->nodes[i]->worst_case_time/dag->nodes[i]->freq_assigned; //so we add time taken actually.
	// 	max_time_taken = max(max_time_taken,multicore->cores[core_id]->total_time);
	// }
	// dag->deadline = 1.2*(max_time_taken);
	// DEADLINE = dag->deadline;


	// applying algo to assign recoveries :


	// cout<<"applying a randomAlgo to give recoveries to some tasks\n";
	// randomAlgo(multicore,dag->deadline);
	// double rsys = systemReliability(multicore);
	// cout<<"System reliability : "<<rsys<<" \n";
	// cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";
	// cout<<"MTTFsys is as:\n"<<systemLifetimeReliability(multicore);
	// cout<<"****************\n";
	// dag->resetAssignment();


	// cout<<"applying a GREEDY Algo to give recoveries to some tasks\n";
	// greedyReliabilityAlgo(multicore,dag->deadline);
	// rsys = systemReliability(multicore);
	// cout<<"System reliability : "<<rsys<<" \n";
	// cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";
	// // cout<<"MTTFsys is as:\n"<<systemLifetimeReliability(multicore);
	// // cout<<"****************\n";
 //    dag->resetAssignment();

	// cout<<"applying a GREEDY TIME Algo to give recoveries to some tasks\n";
	// greedyTimeAlgo(multicore,dag->deadline);
	// rsys = systemReliability(multicore);
	// cout<<"System reliability : "<<rsys<<" \n";
	// cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";

	//calculating reliability


	//create a list of file paths here
	//for each file we create a DAG
	//run that DAG for each algo and get results
	// print the name and result from each algo

	vector<string> files;
	files.push_back("./BenchmarkDagsTXT/Montage_100.txt");

	for(auto filepath : files){
		dag = new DAG();
		dag->inputDAG(filepath,NO_OF_CORES);

		multicore = new Multicore();
		multicore->generateCores(NO_OF_CORES);

		for(int i=0;i<dag->nodes.size();i++)
		{
			int core_id = dag->nodes[i]->core_assigned;
			multicore->cores[core_id]->tasks.push_back(dag->nodes[i]);
		}

		dag->displayDAG();

		assign_freq(dag,2);

		max_time_taken = 0;
		for(int i=0;i<dag->nodes.size();i++){
			int core_id = dag->nodes[i]->core_assigned;
			multicore->cores[core_id]->total_time+=dag->nodes[i]->worst_case_time/dag->nodes[i]->freq_assigned; //so we add time taken actually.
			max_time_taken = max(max_time_taken,multicore->cores[core_id]->total_time);
		}
		dag->deadline = 1.4*(max_time_taken);
		DEADLINE = dag->deadline;

		//randomAlgo
		cout<<"applying a randomAlgo to give recoveries to some tasks\n";
		randomAlgo(multicore,dag->deadline);
		double rsys = systemReliability(multicore);
		cout<<"System reliability : "<<rsys<<" \n";
		cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";
		dag->resetAssignment();


		cout<<"applying a GREEDY Algo to give recoveries to some tasks\n";
		greedyReliabilityAlgo(multicore,dag->deadline);
		rsys = systemReliability(multicore);
		cout<<"System reliability : "<<rsys<<" \n";
		cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";
	    dag->resetAssignment();

		cout<<"applying a GREEDY TIME Algo to give recoveries to some tasks\n";
		greedyTimeAlgo(multicore,dag->deadline);
		rsys = systemReliability(multicore);
		cout<<"System reliability : "<<rsys<<" \n";
		cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";


	}

	return 0;

}