#include <bits/stdc++.h>
#include "system_model.h"
using namespace std;

int NO_OF_CORES = 2;

int NO_OF_TASKS = 10;

double DEADLINE = 0;

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

// schedule: vector with task pointers in a sorted order as per priority
// hasExecuted : bit mask to represent fail/success for a task
// in: the bit index to be fixed in this function call
// i: the index of task in schedule for which we are calculating the Prec
void PrecHelper(Core* core, int hasExecuted,int in,int i,double* PMil){
	if(in==i){
		//calculate PMil
		double slack = DEADLINE - core->total_time;

		for(int j=0;j<i;j++){
			if(hasExecuted & (1 << (j - 1)) == 0){
				double time_rec = core->tasks[j]->worst_case_time/core->tasks[j]->freq_assigned;
				if(slack >= time_rec) slack-=time_rec;
			}
		}
		double time_rec = core->tasks[i]->worst_case_time/core->tasks[i]->freq_assigned;
		if(slack > time_rec){

			double prob_of_mode=1;
			for(int j=0;j<i;j++){

				double ri = taskReliability(core->tasks[j],core);
				if(hasExecuted & (1 << (j - 1)) == 0){
					prob_of_mode*=(1-ri);
				}
				else prob_of_mode*=(ri);
			}

			*PMil+=prob_of_mode;
		}

		return;
	}

	//two cases for bit mask : in'th task fails or succeeds:
	PrecHelper(core,hasExecuted,in+1,i,PMil);
	PrecHelper(core,hasExecuted|1<<i-1,in+1,i,PMil);
}

double taskPreci(Core* c,int i){
	double freq_assigned = c->tasks[i]->freq_assigned;
	c->tasks[i]->freq_assigned = 0.99; //fmax
	double RiFmax = taskReliability(c->tasks[i],c);
	c->tasks[i]->freq_assigned = freq_assigned; //restoring old value;

	double pMil = 0;
	double *PMil = &pMil;

	*PMil = 0;
	PrecHelper(c,0,0,i,PMil);

	return RiFmax*(*PMil);
}

double systemReliability(Multicore* multicore){
	double Rsys = 1;

	for(auto core : multicore->cores){
		int i=0;
		for(auto task : core->tasks){

			double Ri = taskReliability(task,core);
			double Preci = taskPreci(core,i);
			double Rreci = 1 - (1-Ri)*(1-Preci);

			cout<<"Task id : "<<task->id<<" Ri :"<<Ri<<" Preci :"<<Preci<<" Rreci : "<<Rreci<<" \n";

			Rsys*=Rreci;
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


		for(auto l:core->tasks){
			if(l->worst_case_time > slack) break;
			slack-=l->worst_case_time;
			l->recovery_assigned = true;
		}
	}

}

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
	DEADLINE = dag->deadline;


	//applying algo to assign recoveries :


	// cout<<"applying a randomAlgo to give recoveries to some tasks\n";
	// randomAlgo(dag);
	// double rsys = systemReliability(multicore);
	// cout<<"System reliability : "<<rsys<<" \n";
	// cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";
	// cout<<"MTTFsys is as:\n"<<systemLifetimeReliability(multicore);
	// cout<<"****************\n";
	// dag->resetAssignment();


	cout<<"applying a GREEDY Algo to give recoveries to some tasks\n";
	greedyReliabilityAlgo(multicore,dag->deadline);
	cout<<"System reliability : "<<systemReliability(multicore)<<" \n";
	cout<<"Is within deadline : "<<isWithinDeadline(multicore,dag->deadline)<<" \n";
	// cout<<"MTTFsys is as:\n"<<systemLifetimeReliability(multicore);
	// cout<<"****************\n";
	// dag->resetAssignment();


	//calculating reliability





	

	

	return 0;

}