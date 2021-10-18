#include <bits/stdc++.h>
using namespace std;

struct task{
	int id;
	double worst_case_time;
	double priority;
	vector<task*> predecessors;
	vector<task*> successors;
	bool recovery_assigned;
	int core_assigned;
	double freq_assigned;

	task(){
		recovery_assigned = false;
	}
};

struct edge{
	task* src;
	task* dest;
};

class DAG 
{
	//use pointers or ids not the entire struct here:
	//use id and define a global map for id->pointer
	public:
		vector<task*> nodes;
		vector<edge*> edges;

	double deadline;
	// edges[task i]->task j, task k, task l. j,k,l depend on i
	DAG();
	DAG(vector<task*> nodes, vector<edge*> edges, double deadline);

	//sort topologically

	void setDeadline(double deadline);
	void addDependency(int i,int j);
	void removeDependency(int i,int j);
	void generateDAG(int no_of_tasks, int no_of_cores);
	void displayDAG();

};