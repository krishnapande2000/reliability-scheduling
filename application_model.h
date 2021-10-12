struct task{
	int id;
	int worst_case_time;
	int priority;
};

class DAG 
{
	map<task,vector<task>> edges;
	double deadline;
	// edges[task i]->task j, task k, task l. j,k,l depend on i

	void addDependency(task i, task j);
	void removeDependency(task i,task j);

};