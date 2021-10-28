#include <bits/stdc++.h>
#include "system_model.h"
using namespace std;

int NO_OF_CORES = 2;

int NO_OF_TASKS = 10;

Multicore* multicore;

//**********calc start**************************************************************************************************88

double taskRealibility(task* t, Core* c){
	double ro = c->ro;
	double w = c->hardware_coefficient;
	double Fmin = c->Fmin;
	double wc = t->worst_case_time;
	double freq = t->freq_assigned;

	double Fmax = c->Fmax;

	double rf = ro*(pow(10,(w*(1-freq))/(1-Fmin)));

	double Ri = exp((-1)*(rf*wc)/freq);

	cout<<" task id "<<t->id<<" Ri : "<<Ri<<" \n";

	return Ri;
}

double taskRecReliability(task* t, Core* c){
	double Ri = taskRealibility(t,c);

	double Rreci = 1.0 - (1.0-Ri)*(1.0-Ri);

	// cout<<" Rreci : "<<Rreci<<" \n";

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

	return taskRealibility(i,multicore->cores[i->core_assigned]) > taskRealibility(j,multicore->cores[j->core_assigned]);
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

void paperAlgo(){}

void simulatedAnnealingAlgo(){}

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
	return true;
}

//**********constraint checks end**************************************************************************************************88

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

	//applying algo to assign frequencies to tasks

	assign_freq(dag,2);

	double max_time_taken = 0;
	for(int i=0;i<dag->nodes.size();i++){
		int core_id = dag->nodes[i]->core_assigned;
		multicore->cores[core_id]->total_time+=dag->nodes[i]->worst_case_time/dag->nodes[i]->freq_assigned; //so we add time taken actually.
		max_time_taken = max(max_time_taken,multicore->cores[core_id]->total_time);
	}
	dag->deadline = 1.5*(max_time_taken);



	//assign revoeries :


	// cout<<"applying a randomAlgo to give recoveries to some tasks\n";
	// randomAlgo(dag);


	// double rsys = systemReliability(multicore);
	// cout<<"System reliability : "<<rsys<<" \n";
	// cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";

	// cout<<"MTTFsys is as:\n"<<systemLifetimeReliability(multicore);

	// cout<<"****************\n";

	cout<<"applying a GREEDY Algo to give recoveries to some tasks\n";

	greedyReliabilityAlgo(multicore,dag->deadline);
	cout<<"System reliability : "<<systemReliability(multicore)<<" \n";
	cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";


	//calculating reliability





	

	

	return 0;

}