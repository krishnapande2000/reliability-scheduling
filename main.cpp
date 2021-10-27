#include <bits/stdc++.h>
#include "system_model.h"
using namespace std;

int NO_OF_CORES = 4;

int NO_OF_TASKS = 20;

Multicore* multicore;

double taskRealibility(task* t, Core* c){
	double ro = c->ro;
	double w = c->hardware_coefficient;
	double Fmin = c->Fmin;
	double wc = t->worst_case_time;
	double freq = t->freq_assigned;

	double Fmax = c->Fmax;

	double rf = ro*(pow(10,(w*(1e9-freq))/(1e9-Fmin)));

	double Ri = exp((-1)*(rf*wc)/freq);

	cout<<" ro :"<<ro<<" freq: "<<freq<<" wc: "<<wc<<" rf: "<<rf<<" Ri : "<<Ri<<" \n";

	return Ri;
}

double taskRecReliability(task* t, Core* c){
	double Ri = taskRealibility(t,c);

	double Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);

	cout<<" Rreci : "<<Rreci<<" \n";

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

double lifetimeReliability(task* t,Core* c){
	double constant_a = 38.92e9;
	double freq = t->freq_assigned;
	double MTTF = constant_a / (freq-1e9); //freq in ghz please (10^9) 
	//TO-DO(niyati): multiply by root Ao.

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

void randomAlgo(DAG* dag){
	for(task* t : dag->nodes){
		int r = rand()%2;
		if(r) t->recovery_assigned = true;
		t->freq_assigned = rand()%(3000000000-1500000000 + 1) + 1.5e9;
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

	return taskRealibility(i,multicore->cores[i->core_assigned]) > taskRealibility(j,multicore->cores[j->core_assigned]);
}

void greedyAlgo(Multicore* multicore, double deadline){

	//assign priority as per the dependency and reliability

	for(auto core : multicore->cores){
		//tasks vector to be sorted
		sort(core->tasks.begin(),core->tasks.end(),comp);
		double slack = deadline - core->total_time;

		for(auto l:core->tasks){
			if(l->worst_case_time > slack) break;
			slack-=l->worst_case_time;
			l->recovery_assigned = true;
		}
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
	multicore = new Multicore();
	multicore->generateCores(NO_OF_CORES);


	cout<<"generating DAG\n";
	//create or input a DAG
	DAG* dag = new DAG();
	dag->generateDAG(NO_OF_TASKS,NO_OF_CORES);

	for(int i=0;i<dag->nodes.size();i++){
		int core_id = dag->nodes[i]->core_assigned;
		multicore->cores[core_id]->tasks.push_back(dag->nodes[i]);
	}



	// int n = 10;
	dag->displayDAG();

	//applying algo to assign recovery, frequencies to tasks
	cout<<"applying a randomAlgo to give recoveries to some tasks\n";
	randomAlgo(dag);

	for(int i=0;i<dag->nodes.size();i++){
		int core_id = dag->nodes[i]->core_assigned;
		multicore->cores[core_id]->total_time+=dag->nodes[i]->worst_case_time/dag->nodes[i]->freq_assigned; //so we add time taken actually.
	}


	cout<<"System reliability : "<<systemReliability(multicore)<<" \n";
	cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";

	cout<<"MTTFsys is as:\n"<<systemLifetimeReliability(multicore);

	cout<<"****************\n";

	cout<<"applying a GREEDY Algo to give recoveries to some tasks\n";

	greedyAlgo(multicore,dag->deadline);
	cout<<"System reliability : "<<systemReliability(multicore)<<" \n";
	cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";


	//calculating reliability





	

	

	return 0;

}