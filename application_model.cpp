#include <system_model.h>

#include <application_model.h>
#include <reliability_model.cpp>

class DAG 
{
	map<task,vector<task>> edges;
	double deadline;
	// edges[task i]->task j, task k, task l. j,k,l depend on i
	DAG(){

	}

	DAG(map<task,vector<task>> edges, double deadline){
		this.edges = edges;
		this.deadline = deadline;
	}

	void setDeadline(double deadline){
		this.deadline = deadline;
	}

	void addDependency(task i, task j){

	}
	
	void removeDependency(task i,task j){

	}


};