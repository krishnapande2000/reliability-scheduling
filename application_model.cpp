#include "system_model.h"
#include "application_model.h"

double generateRand(int m, int n) {
    double a = rand()%n;
    double b = rand()%m;
    return a/b;
}

class DAG 
{
	vector<task*> nodes;
	vector<edge*> edges;
	double deadline;
	// edges[task i]->task j, task k, task l. j,k,l depend on i
	DAG(){

	}

	DAG(vector<task*> nodes, vector<edge*> edges, double deadline){
		this.edges = edges;
		this.nodes = nodes;
		this.deadline = deadline;
	}

	void setDeadline(double deadline){
		this.deadline = deadline;
	}

	void addDependency(int i, int j){
		nodes[j]->predecessors.push_back(nodes[i]);
		nodes[i]->successors.push_back(nodes[j]);

		edge* newEdge = new edge();
		newEdge->src = nodes[i];
		newEdge->dest = nodes[j];
	}
	
	void removeDependency(int i,int j){
		for(int i=0;i<edges.size();i++){
			if(edges[i].src == nodes[i] && edges[i].dest == nodes[j]){
				edges.remove(edges.begin()+i);
				break;
			}
		}
	}

	void generateDAG(int no_of_tasks, int no_of_cores){

		this.deadline = rand()%(100-70 + 1) + 70;

		for(int i=0;i<no_of_tasks;i++){
			task* newTask = new task();
			newTask->id = i;
			newTask->worst_case_time = generateRand(10,100);
			newTask->core_assigned = rand()%no_of_cores;

			nodes.push_back(newTask);
		}

		for(int i=0;i<no_of_tasks;i++){
			for(int j=i+1;j<no_of_tasks;j++){
				int r = rand()%2;
				if(r){
					nodes[j]->predecessors.push_back(nodes[i]);
					nodes[i]->successors.push_back(nodes[j]);

					edge* newEdge = new edge();
					newEdge->src = nodes[i];
					newEdge->dest = nodes[j];

				}
			}
		}

	}


};