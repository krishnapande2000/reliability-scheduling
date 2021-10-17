#include <bits/stdc++.h>
#include "system_model.h"
#include "application_model.h"
using namespace std;

int NO_OF_CORES = 4;

int NO_OF_TASKS = 20;

double taskRealibility(task* t){
	double ro = t->.c.ro;
	double w = t.c.w;
	double Fmin = t.c.Fmin;
	double wc = t.worst_case_time;
	double freq = t.freq_assigned;

	double rf = ro*(power(10,(w*(1-freq))/(1-Fmin)));

	double Ri = power(e,(-1)*(rf*wc)/freq);

	return Ri;
}

double taskRecReliability(task* t){
	double Ri = taskRealibility(t);
	// double Prec = ?

	double Rreci = 1 - (1-Ri)*(1-Ri);

	return Rreci;
}

double systemReliability(vector<task*> tasks){
	double Rsys = 1;

	for(task t: tasks){
		double Ri;
		if(t.recovery_Assigned) Ri = taskRecReliability(task);
		else Ri = taskRealibility(task);
		Rsys = Rsys*Ri;
	}

	return Rsys;
}

double lifetimeReliability(core c){
	double MTTF =1;
	return MTTF;
	
}

void randomAlgo(DAG dag){
	for(task* t : dag.nodes){
		int r = rand()2;
		if(r) t->recovery_assigned = true;
		t->freq_assigned = rand()%(500-200) + 200;
	}
}

bool isWithinDeadline(Multicore multicore){

}

bool hasMTTF(Multicore multicore){

}

void main(){
	//generating cores
	Multicore multicore;
	multicore.generateCores(NO_OF_CORES);

	//create or input a DAG
	DAG dag;
	dag.generateDAG(NO_OF_TASKS,NO_OF_CORES);

	for(int i=0;i<dag.nodes.size();i++){
		int core_id = dag.nodes[i]->core_assigned;
		multicore.cores[i]->tasks.push_back(dag.nodes[i]);
	}



	int n = 10;

	//applying algo to assign recovery, frequencies to tasks
	randomAlgo(dag);



	//calculating reliability



}