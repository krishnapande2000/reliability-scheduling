#include <bits/stdc++.h>
#include "system_model.h"
using namespace std;

int NO_OF_CORES = 4;

int NO_OF_TASKS = 20;

double taskRealibility(task* t, Core* c){
	double ro = c->ro;
	double w = c->hardware_coefficient;
	double Fmin = c->Fmin;
	double wc = t->worst_case_time;
	double freq = t->freq_assigned;

	double rf = ro*(pow(10,(w*(1.0-freq))/(1.0-Fmin)));

	double Ri = exp((-1)*(rf*wc)/freq);

	return Ri;
}

double taskRecReliability(task* t, Core* c){
	double Ri = taskRealibility(t,c);

	double Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);

	return Rreci;
}

double systemReliability(Multicore* multicore){
	double Rsys = 1;

	for(auto core : multicore->cores){
		for(auto task : core->tasks){

			double Ri = 1;
			if(task->recovery_assigned){
				Ri = taskRecReliability(task,core);
			}
			else{
				Ri = taskRealibility(task,core);
			}

			Rsys = Rsys*Ri;
		}
	}

	return Rsys;
	
}

double lifetimeReliability(Core* c){
	double MTTF =1;
	return MTTF;
	
}

double systemLifetimeReliability(Multicore* multicore){
	
	return 1;
}

void randomAlgo(DAG* dag){
	for(task* t : dag->nodes){
		int r = rand()%2;
		if(r) t->recovery_assigned = true;
		t->freq_assigned = rand()%(500-200) + 200;
	}
}

bool comp(task* i,task* j){
	bool iToj = false, jToi = false;
	for(auto l:i->predecessors){
		if(l==j) jToi = true;
	}
	for(auto l:j->predecessors){
		if(l==i) iToj = true;
	}

	if(jToi) return false;
	if(iToj) return true;

	return taskRealibility(i) > taskRealibility(j);
}

void greedyAlgo(Multicore* multicore){

	//assign priority as per the dependency and reliability

	for(auto core : multicore->cores){
		//tasks vector to be sorted
		sort(core->tasks.begin,core->tasks.end(),comp);
	}

	
}

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
	return true;
}

int main(){

	cout<<"generating cores\n";
	//generating cores
	Multicore* multicore = new Multicore();
	multicore->generateCores(NO_OF_CORES);


	cout<<"generating DAG\n";
	//create or input a DAG
	DAG* dag = new DAG();
	dag->generateDAG(NO_OF_TASKS,NO_OF_CORES);

	for(int i=0;i<dag->nodes.size();i++){
		int core_id = dag->nodes[i]->core_assigned;
		multicore->cores[core_id]->tasks.push_back(dag->nodes[i]);
	}



	int n = 10;

	//applying algo to assign recovery, frequencies to tasks
	cout<<"applying a randomAlgo to give recoveries to some tasks\n";
	randomAlgo(dag);



	//calculating reliability





	dag->displayDAG();

	cout<<"System reliability : "<<systemReliability(multicore)<<" \n";
	cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";

	return 0;

}