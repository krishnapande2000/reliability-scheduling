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

int NO_OF_CORES = 4;

int NO_OF_TASKS = 50;

double DEADLINE = 0.0;

double LIFETIME_THRESHOLD = 60;

double WC_ADJ_FACTOR = 1;

bool SHORTEST_JOB_FIRST = true;

double RO = 1e-6;

int NUMBER_OF_RECOVERIES = 0;

int NO_OF_INTERVALS = 1000;

Multicore* multicore;


//**********calc start**************************************************************************************************88

double taskReliability(task* t, Core* c){
	double ro = RO;
	double w = c->hardware_coefficient;
	double Fmin = c->Fmin;
	double Fmax = c->Fmax;

	double wc = t->worst_case_time;
	double freq = t->freq_assigned;

	double rf = ro*(pow(10,(w*(1-freq))/(1-Fmin)));
	double Ri = exp((-1)*(rf*wc));

	// cout<<" task id "<<t->id<<" Ri : "<<Ri<<" \n";

	return Ri;
}

double taskRecReliability(task* t, Core* c){
	double Ri = taskReliability(t,c);
	double Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);
	return Rreci;
}

double taskReliability_dynamic(task* t, Core* c){
	double ro = c->ro;
	double w = c->hardware_coefficient;
	double Fmin = c->Fmin;
	double Fmax = c->Fmax;

	double exec_time = t->execution_time;
	double freq = t->freq_assigned;

	double rf = ro*(pow(10,(w*(1-freq))/(1-Fmin)));
	double Ri = exp((-1)*(rf*exec_time));

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

	return RiFmax*(*PMil);
}

double systemReliability(DAG* dag){

	double Rsys = 1;
	for(int i=1;i<=NO_OF_TASKS;i++){

		task* task = dag->nodes[i];
		Core* core = multicore->cores[task->core_assigned];
		double Ri = taskReliability(task,core);
		double Rreci = Ri;
		if(task->recovery_assigned){
			Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);
		}

		Rsys*=Rreci;
	}
	
	return Rsys;
}

double systemReliability_dynamic(DAG* dag){

	double Rsys = 1;
	for(int i=1;i<=NO_OF_TASKS;i++){

		task* task = dag->nodes[i];
		Core* core = multicore->cores[task->core_assigned];
		double Ri = taskReliability_dynamic(task,core);
		double Rreci = Ri;
		if(task->recovery_assigned){
			Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);
		}

		Rsys*=Rreci;
	}
	
	return Rsys;
}

double lifetime_f_value(double freq, double time){
	// double power_of_e = 
	double value = exp((-1)*3.83*time*time*freq*freq*1e-3);
	return value;
}

//calculating integral using simpsons rule
double lifetimeReliability_approx(task* t){

	double ni = NO_OF_INTERVALS;
	double exec_time = t->execution_time;
	double len_of_interval = exec_time/ni;
	double freq = t->freq_assigned;

	double sum = 1 + lifetime_f_value(freq,exec_time);
	double idle_constant = 0;

	for(int k=1;k<ni;k++){

		double xk = k*(len_of_interval);
		double term;
		if(k%2){
			term = 4*lifetime_f_value(freq,xk);
		}
		else{
			term = 2*lifetime_f_value(freq,xk);
		}
		idle_constant = max(term,idle_constant);
		sum+=term;
	}

	sum +=  idle_constant*(task->core_assigned->idle_time);
	double result = (len_of_interval*sum)/3;
	return result;

}


double systemLifetimeReliability(Multicore* multicore){
	double MTTFsys = DBL_MAX;
	for(auto core : multicore->cores){
		double MTTFi = 0;
		for(auto task : core->tasks){
			if(task->worst_case_time == 0) continue; //to skip vexit and ventry nodes in calc
			MTTFi += lifetimeReliability_approx(task);
			cout<<"task id "<<task->id<<" freq "<<task->freq_assigned<<" time "<<task->worst_case_time<<" mttf val "<<MTTFi<<"\n";
		}
		MTTFsys = min(MTTFsys,MTTFi);
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

struct CompWorstCaseTime{
	bool operator()(task* i,task* j){
	return i->worst_case_time < j->worst_case_time;
	}
};
struct CompEndTime{
	bool operator()(task* i,task* j){
	return i->end_time > j->end_time;
	}
};
struct CompVexitDistance{
	bool operator()(task* i,task* j){
	return i->vexit_dist > j->vexit_dist;
	}
};

bool compWorstCaseTime(task* i,task* j){
	if(SHORTEST_JOB_FIRST) return i->worst_case_time < j->worst_case_time;
	return i->worst_case_time > j->worst_case_time;
}

struct CompDynamicEndTime{
	bool operator()(task* i,task* j){
	return i->dynamic_process_end_time > j->dynamic_process_end_time;
	}
};


//as the tasks are already topologically sorted
//this function uses the order to calculate 
// dist from vexit for every node
void longestDistVexit(DAG* dag)
{
    stack<task*> topo_stack;
   	
   	int vexit_id = dag->vexit_id;
   	int V = vexit_id+1;
    // Mark all the vertices as not visited
    vector<bool> visited(V,false);
   
    // Call the recursive helper function to store Topological
    // Sort starting from all vertices one by one
    // for (int i = 0; i < V; i++){
    //     if (visited[i] == false){
    //         topologicalSortUtil(dag, i, visited, topo_stack);
    //     }
    // }
   
    // Initialize distances to all vertices as infinite and
    // distance to source as 0
    for (auto node:dag->nodes){
        node->vexit_dist = INT_MIN;
        topo_stack.push(node);
    }
    dag->nodes[vexit_id]->vexit_dist = 0;

    // Process vertices in topological order
    while (!topo_stack.empty()) {
        // Get the next vertex from topological order
        task* nodeu = topo_stack.top();
        topo_stack.pop();
   
        if (nodeu->vexit_dist != INT_MIN) {
            for (auto pred: nodeu->predecessors){
             
                if (pred->vexit_dist < nodeu->vexit_dist + 1){
                    pred->vexit_dist = nodeu->vexit_dist + 1;
                }
            }
        }
    }
}

// returns a new dag with 
// 1) duplicated nodes for tasks with recovery as per criteria
// 2) added dummy nodes vexit and ventry
// 3) dist from vexit is updated by calling a function
DAG* getNewDAG(DAG* dag, int no_of_recoveries){
	DAG* new_dag = new DAG();
	new_dag->inputDAG(dag->filesname, NO_OF_CORES);

	new_dag->resetAssignment();
	
	vector<task*> sorted_nodes = new_dag->nodes;
	sort(sorted_nodes.begin(), sorted_nodes.end(),compWorstCaseTime);
	for (int i = 0; i < no_of_recoveries; i++){
		task* node = sorted_nodes[i];
		node->recovery_assigned=true;
		task* new_node = new task(); //successor - prede.. done
		*new_node = *node;
		new_node->id = new_dag->nodes.size() + 1;
		new_dag->nodes.push_back(new_node);
		for(auto pred : node->predecessors){
			pred->successors.push_back(new_node);
		}
		for(auto suc :  node->successors){
			suc->predecessors.push_back(new_node);
		}
		// cout<<"recovery: "<<node->id<<" to "<<new_node->id<<endl;
	}

	task* ventry = new task();
	ventry->worst_case_time = 0;
	ventry->start_time = 0;
	ventry->end_time = 0;

	for(auto node:new_dag->nodes){
		if(node->predecessors.empty()){
			ventry->successors.push_back(node);
			node->predecessors.push_back(ventry);
		}
	}

	task* vexit = new task();
	vexit->worst_case_time = 0;
	vexit->start_time = 0;
	vexit->end_time = 0;

	for(auto node: new_dag->nodes){
		if(node->successors.empty()){
			vexit->predecessors.push_back(node);
			node->successors.push_back(vexit);
		}
	}

	ventry->id = 0;
	vexit->id = new_dag->nodes.size() + 1;
	new_dag->nodes.insert(new_dag->nodes.begin(),ventry);
	new_dag->nodes.push_back(vexit);
	new_dag->ventry_id = ventry->id;
	new_dag->vexit_id = vexit->id;

	// cout<<"this is the new DAG with backups and all";
	// new_dag->displayDAG(); //compile and check

	// cout<<"OLD DAG";
	// dag->displayDAG(); //compile and check

	longestDistVexit(new_dag);

	return new_dag;
}


//checks if given DAG can be scheduled within deadline
// assigns tasks to cores if possible
// also assigns start and end time for each task
bool canScheduleStatic(DAG* node_graph){
	multicore->clear_cores();
	vector<Core*> free_cores;
	for(auto core:multicore->cores){
		free_cores.push_back(core);
		core->free_at = 0;
	}

	//put frequency assignment here
	// try to already give max for each that is mixed levels like 4 and 5 as well
	// so that it doesnt happen that a greater frequency was possible within deadline
	// but we assigned 4 and said it wasnt possible to schedule it at all with all constraints
	// (we can calculate mttf beforehand also cause we have homogenous cores)

	for(int i=multicore->freq_levels.size()-1;i>=0;i--){

		//
		for(task* node : node_graph->nodes){
			node->freq_assigned = freq_levels[i];
		}

	}

	priority_queue<task*, vector<task*>,CompEndTime> processing_queue;
	priority_queue<task*, vector<task*>,CompVexitDistance> ready_queue;

	for(task* node : node_graph->nodes){
		node->remaining_predecessors = node->predecessors.size();
		node->core_assigned = 0;
	}
	// cout<<free_cores.size()<<endl;
	free_cores.back()->tasks.push_back(node_graph->nodes[node_graph->ventry_id]);
	node_graph->nodes[node_graph->ventry_id]->core_assigned = free_cores.back()->id;
	free_cores.pop_back();

	// cout<<node_graph->nodes[node_graph->ventry_id]->core_assigned<<" free core size - "<<free_cores.size()<<endl;
	processing_queue.push(node_graph->nodes[node_graph->ventry_id]);


	while (!processing_queue.empty()){
		task* top_node = processing_queue.top();
		// cout<<top_node->id<<" * ";
		processing_queue.pop();
		double end_time = top_node->end_time;
		// cout<<"end time -- "<<end_time<<" - ";
		while(true){

			Core* core_assigned = multicore->cores[top_node->core_assigned];
			for(task* child : top_node->successors){
				child->remaining_predecessors --;
				if(child->remaining_predecessors == 0){
					child->ready_time = end_time;
					ready_queue.push(child);
					// cout<<"in ready q: "<<child->id<<" ";
				}
			}
			
			core_assigned->free_at = end_time;
			free_cores.push_back(core_assigned);

			// cout<<endl;
			// cout<<" ready q size -- "<<ready_queue.size()<<" processing q size -- "<<processing_queue.size()<<endl;

			// if processing queue has same end time for other nodes
			// remove those nodes too 
			if(!processing_queue.empty() && processing_queue.top()->end_time == end_time){
				top_node = processing_queue.top();
				processing_queue.pop();
			}
			else break;
			
		}

		// now assigning free cores to ready queue tasks
		while(!ready_queue.empty()){
			task* next_node = ready_queue.top();
			for(int i=0; i<free_cores.size(); i++){
				Core* free_core = free_cores[i];
				double prob_start_time = max(next_node->ready_time, free_core->free_at);
				double prob_end_time = prob_start_time + WC_ADJ_FACTOR*(next_node->worst_case_time);
				if(prob_end_time <= DEADLINE) {
					next_node->core_assigned = free_core->id;
					next_node->start_time = prob_start_time;
					next_node->end_time = prob_end_time;
					free_core->tasks.push_back(next_node);
					ready_queue.pop();
					processing_queue.push(next_node);
					free_cores.erase(free_cores.begin()+i);
					break;
				}
				else if(i==free_cores.size()-1) return false;
			}
			if(free_cores.size() == 0) break;
		}
	}
	if(!ready_queue.empty()) return false;
	return true;

}

//returns max no of recoveries possible 
int numberOfRecoveriesStatic(DAG* dag){
	int n = dag->nodes.size();
	int start =0;
	int end = n;

	int mid;

	while(end>=start){
		mid = (start+end)/2;
		DAG* new_dag_mid = getNewDAG(dag, mid);
		bool can_schedule_new_dag_mid = canScheduleStatic(new_dag_mid);
		// cout<<"start: "<<start<<" end: "<<end<<" mid: "<<mid<<endl;
		if(start == end ) return start;
		if(start == end - 1){
			DAG* new_dag_end = getNewDAG(dag, end);
			bool a = canScheduleStatic(new_dag_end);
			if(a) return end;
			else return start;
		}
		if(can_schedule_new_dag_mid) start = mid;
		else end=mid-1;
		
	}

	return mid;
}

double dynamicScaleFrequency(task* node){
	double time_available = node->end_time - node->dynamic_process_start_time;
	double freq_factor = 1.0 ;
	double max_freq = multicore->freq_levels[multicore->freq_levels.size()-1];
	// if number of recoveries is equal to the number of nodes
	if(NUMBER_OF_RECOVERIES == NO_OF_TASKS){
		for(int i=0; i<=multicore->freq_levels.size()-2; i++){
			if(node->worst_case_time * (max_freq/multicore->freq_levels[i]) <= time_available){
				freq_factor = max_freq/multicore->freq_levels[i];
				node->freq_assigned = 0.99 / freq_factor;
				break;
			}
		}
	}
	else{
		for(int i = multicore->freq_levels.size()-2; i>=0; i--){
			if(node->worst_case_time * (max_freq/multicore->freq_levels[i]) <= time_available 
			&& systemLifetimeReliability(multicore) < LIFETIME_THRESHOLD){
				freq_factor = max_freq/multicore->freq_levels[i];
				node->freq_assigned = 0.99 / freq_factor;
			}
			else break;
		}	
	}
	node->freq_assigned = 0.99 / freq_factor;
	return freq_factor;
}

double dynamicScaleFrequency_new(task* node){
	double original_freq_assigned = node->freq_assigned;
	double freq_factor = 1;
	for(int i = multicore->freq_levels.size()-1; i>=0; i--){
			freq_factor = max_freq/multicore->freq_levels[i];
			node->freq_assigned = 0.99 / freq_factor
			if(systemLifetimeReliability(multicore) < LIFETIME_THRESHOLD){
				node->freq_assigned = original_freq_assigned;
			}
			else break;
		}	
	return 1/freq_factor;
}

void executeTask(task* node, double freq_factor){
	double factor = (double)rand() / (double) RAND_MAX;
	if(factor<0.5) factor+=0.5;
	if(factor>0.9) {
		factor -= (double)(rand() % 4) / 10.0 - 0.1;
	}
	double original_exec_time = node->execution_time;
	node->execution_time = node->worst_case_time * factor * freq_factor;
	multicore->cores[node->core_assigned]->idle_time = original_exec_time - node->execution_time;
	node->dynamic_process_end_time = node->dynamic_process_start_time + node->execution_time;
}

// scales the frequency of all tasks to adjust as per MTTF
// dynamically at run time using the schedule created by static
void scheduleDynamic(DAG* node_graph){
	for(Core* core:multicore->cores){
		core->free_at = 0;
		for(task* node: core->tasks){
			core->static_tasks_assigned.push(node);

		}
	}

	for(task* node : node_graph->nodes){
		node->dynamic_process_start_time = -1.0;
		node->dynamic_process_end_time = -1.0;
		node->dynamic_ready_time = -1.0;
		node->remaining_predecessors = node->predecessors.size();
	}
	
	node_graph->nodes[node_graph->ventry_id]->dynamic_process_end_time = 0.0;
	node_graph->nodes[node_graph->ventry_id]->dynamic_process_start_time = 0.0;

	priority_queue<task*, vector<task*>,CompDynamicEndTime> processing_queue;
	processing_queue.push(node_graph->nodes[node_graph->ventry_id]);
	while(!processing_queue.empty()){
		task* node = processing_queue.top();
		processing_queue.pop();
		double curr_time = node->dynamic_process_end_time;

		Core* core_assigned_to_node = multicore->cores[node->core_assigned];
		core_assigned_to_node->free_at = curr_time;
		core_assigned_to_node->static_tasks_assigned.pop();
		// next task in the core already ready to be executed now
		if(!core_assigned_to_node->static_tasks_assigned.empty()){
			task* core_top_node = core_assigned_to_node->static_tasks_assigned.front();
			if(core_top_node->dynamic_ready_time != -1.0){
				core_top_node->dynamic_process_start_time = curr_time;
				processing_queue.push(core_top_node);
				// frequency scale 
				double freq_factor = dynamicScaleFrequency(core_top_node);
				executeTask(core_top_node, freq_factor);
			}
		}

		for(task* child : node->successors){
			child->remaining_predecessors --;
			if(child->remaining_predecessors == 0){
				Core* child_core = multicore->cores[child->core_assigned];
				child->dynamic_ready_time = curr_time;
				if(child_core->static_tasks_assigned.front()==child){
					child->dynamic_process_start_time = curr_time;
					processing_queue.push(child);
					// frequency scale 
					double freq_factor = dynamicScaleFrequency(child);
					executeTask(child, freq_factor);
				}
			}
		}

	}
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

	return systemLifetimeReliability(multicore) >= LIFETIME_THRESHOLD;
}

//**********constraint checks end**************************************************************************************************88

//converts a given line string to a vector of words seperated as per spaces
//extra spaces are skipped
vector<string> line_to_words(string str) 
{ 
    // Used to split string around spaces. 
    vector<string> res;
    istringstream ss(str); 

    string word; // for storing each word 
  
    // Traverse through all words 
    // while loop till we get  
    // strings to store in string word 
    while (ss >> word)  
    { 
    	res.push_back(word);
        //push the read word       
    } 
    return res;
}

int executor(string config_filename){

	ifstream input_file;
	string config_filepath = "./INPUTS/" + config_filename + ".txt";
	input_file.open(config_filepath);

	string output_filename = "./OUTPUTS/output_" + config_filename + ".txt";
	ofstream output_file;
	output_file.open(output_filename);

	string csv_output_filename = "./OUTPUTS/output_" + config_filename + ".csv";
	ofstream csv_output_file;
	csv_output_file.open(csv_output_filename);

	csv_output_file<<"Name of file,No of recoveries,Rsys_static,PoF_static,Rsys_dynamic,PoF_dynamic\n";

	string line="";
	int i=1;

	while(getline(input_file,line,'\n'))
	{
		vector<string> line_w=line_to_words(line);
		NO_OF_CORES = stoi(line_w[0]);
		WC_ADJ_FACTOR = stof(line_w[1]);
		SHORTEST_JOB_FIRST = stoi(line_w[2]);
		RO = stof(line_w[3]);

		output_file<<" CONFIG NUMBER : "<<i<<" \n";
		output_file<<line<<endl;

		csv_output_file<<"Config Number,"<<i<<"\n";
		csv_output_file<<NO_OF_CORES<<","<<WC_ADJ_FACTOR<<","<<SHORTEST_JOB_FIRST<<","<<RO<<"\n";

		string benchmark_folder = "./BenchmarkDagsTXT/";

		vector<string> files;
		files.push_back("Montage_25.txt");
		files.push_back("Montage_50.txt");
		files.push_back("CyberShake_30.txt");
		files.push_back("CyberShake_50.txt");
		files.push_back("Inspiral_30.txt");
		files.push_back("Inspiral_50.txt");
		files.push_back("Sipht_30.txt");
		files.push_back("Sipht_60.txt");

		for(auto filename : files){

			DAG* dag = new DAG();
			string filepath = benchmark_folder+filename;
			dag->inputDAG(filepath,NO_OF_CORES);
			NO_OF_TASKS = dag->nodes.size();

			multicore = new Multicore();
			multicore->generateCores(NO_OF_CORES);

			// assign frequency
			assign_freq(dag,2);

			// create a reasonable deadline
			double total_time = 0;
			for(auto node:dag->nodes){
				total_time+=node->worst_case_time;
			}
			dag->deadline = 1.75*(total_time/(double)NO_OF_CORES);
			DEADLINE = dag->deadline;
			
			//use algo
			int no_of_recoveries = numberOfRecoveriesStatic(dag);
		    NUMBER_OF_RECOVERIES = no_of_recoveries;
			DAG* new_dag = getNewDAG(dag,no_of_recoveries);
		    canScheduleStatic(new_dag);

			assign_freq(new_dag,2);

			double Rsys_static = systemReliability(new_dag);
			double PoF_static = 1 - Rsys_static;

			scheduleDynamic(new_dag);
			double Rsys_dynamic = systemReliability_dynamic(new_dag);
			double PoF_dynamic = 1 - Rsys_dynamic;

			//print reliability
			 output_file<<setw(7)<<filename<<setw(8)<<no_of_recoveries<<"\t"<<Rsys_static<<"\t"<<Rsys_dynamic<<"\t"<<PoF_dynamic<<endl;
			 csv_output_file<<filename<<","<<no_of_recoveries<<","<<Rsys_static<<","<<PoF_static<<","<<Rsys_dynamic<<","<<PoF_dynamic<<"\n";
		}

		i++;
		output_file<<"---------------------------------------------------------\n\n";
		csv_output_file<<"\n\n";
	}

	output_file.close();
	input_file.close();

}


int main(){

	// string config_filename;
	// cin>>config_filename;
	// executor(config_filename);
	
	// applying algo to assign recoveries :

	vector<string> files;
	files.push_back("./BenchmarkDagsTXT/Montage_25.txt");


	cout<<"Name of file \t Rsys \n";

	for(auto filepath : files){
		DAG* dag = new DAG();
		NO_OF_CORES = 10;
		dag->inputDAG(filepath,NO_OF_CORES);
		NO_OF_TASKS = dag->nodes.size();

		multicore = new Multicore();
		multicore->generateCores(NO_OF_CORES);

		// assign frequency
		assign_freq(dag,2);

		// create a reasonable deadline
		double total_time = 0;
		for(auto node:dag->nodes){
			total_time+=node->worst_case_time;
		}
		dag->deadline = 1.75*(total_time/(double)NO_OF_CORES);
		DEADLINE = dag->deadline;
		cout<<DEADLINE<<endl;
		//use algo

		int no_of_recoveries = numberOfRecoveriesStatic(dag);
		cout<<" NO OF RECOVERIES "<<no_of_recoveries;
		NUMBER_OF_RECOVERIES = no_of_recoveries;
		cout<<endl;
		// dag->displayDAG();
		DAG* new_dag = getNewDAG(dag, no_of_recoveries);
		// new_dag->displayDAG();
		assign_freq(new_dag,1);
		bool a = canScheduleStatic(new_dag);
		cout<<endl<<endl;
		cout<<"can schedule true or false -- "<<a<<endl;

		for(task* node:new_dag->nodes){
			cout<<" Node ID:"<<node->id;
			cout<<" Core assigned:"<<node->core_assigned;
			cout<<" Start time:"<<node->start_time;
			cout<<" End time:"<<node->end_time;
			cout<<" Frequency Assigned:"<<node->freq_assigned;
			cout<<endl;
		}

		cout<<"MTTF calcccc \n";
		double MTTFsys = systemLifetimeReliability(multicore);
		cout<<"MTTF sys is "<<MTTFsys<<endl;
		cout<<endl<<endl;

		// for(Core* core: multicore->cores){
		// 	cout<<" Core ID:"<<core->id<<endl;
		// 	for(task*  node: core->tasks){
		// 		cout<<" Node ID:"<<node->id;
		// 		cout<<" Core assigned:"<<node->core_assigned;
		// 		cout<<" Start time:"<<node->dynamic_process_start_time;
		// 		cout<<" End time:"<<node->dynamic_process_end_time;
		// 		cout<<endl;
		// 	}
		// 	cout<<endl;
		// }
		cout<<endl<<endl<<"DYNAMIC----"<<endl<<endl;
		scheduleDynamic(new_dag);

		for(task* node:new_dag->nodes){
			cout<<" Node ID:"<<node->id;
			cout<<" Core assigned:"<<node->core_assigned;
			cout<<" Start time:"<<node->dynamic_process_start_time;
			cout<<" End time:"<<node->dynamic_process_end_time;
			cout<<" Frequency Assigned:"<<node->freq_assigned;
			cout<<endl;
		}

		double Rsys = systemReliability_dynamic(new_dag);


		//print reliability
		cout<<hasMTTF(multicore)<<"this is mttf check\n";
		 cout<<filepath<<"\t"<<Rsys<<"\t"<<endl;


	}

	return 0;

}



/*
CODE DUMP

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
*/