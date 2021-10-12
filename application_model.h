#include <system_model.h>

struct task{
	int id;
	double worst_case_time;
	double priority;
	bool recovery_assigned;
	core core_Assigned;
	double freq_Assigned;

	task(){
		recovery_assigned = false;
	}
};

class DAG 
{
	map<task,vector<task>> edges;
	double deadline;
	// edges[task i]->task j, task k, task l. j,k,l depend on i
	DAG();
	DAG(map<task,vector<task>> edges, double deadline);

	void setDeadline(double deadline);
	void addDependency(task i, task j);
	void removeDependency(task i,task j);

};