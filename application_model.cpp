#include <cctype> 
#include <cerrno> 
#include <cfloat> 
#include <ciso646> 
#include <climits> 
#include <clocale> 
#include <cmath> 
#include <ctime> 

  
// C++ 
#include <algorithm> 
#include <bitset> 
#include <complex> 
#include <deque> 
#include <exception> 
#include <fstream> 
#include <functional> 
#include <iomanip> 
#include <ios> 
#include <iosfwd> 
#include <iostream> 
#include <istream> 
#include <iterator> 
#include <limits> 
#include <list> 
#include <locale> 
#include <map> 
#include <memory> 
#include <new> 
#include <numeric> 
#include <ostream> 
#include <queue> 
#include <set> 
#include <sstream> 
#include <stack> 
#include <stdexcept> 
#include <streambuf> 
#include <string> 
#include <typeinfo> 
#include <utility> 
#include <valarray> 
#include <vector> 
  
#if __cplusplus >= 201103L 
#include <array> 
#include <atomic> 
#include <chrono> 
#include <condition_variable> 
#include <forward_list> 
#include <future> 
#include <initializer_list> 
#include <mutex> 
#include <random> 
#include <ratio> 
#include <regex> 
#include <scoped_allocator> 
#include <system_error> 
#include <thread> 
#include <tuple> 
#include <typeindex> 
#include <type_traits> 
#include <unordered_map> 
#include <unordered_set> 
#endif 
#include "application_model.h"

using namespace std;

double generateRand(int m, int n) {
    double a = rand()%n;
    double b = rand()%m;
    return a/b;
}


DAG::DAG(){

}

DAG::DAG(vector<task*> nodes, vector<edge*> edges, double deadline){
		this->edges = edges;
		this->nodes = nodes;
		this->deadline = deadline;
	}

void DAG::setDeadline(double deadline){
		this->deadline = deadline;
	}

void DAG::addDependency(int i, int j){
		nodes[j]->predecessors.push_back(nodes[i]);
		nodes[i]->successors.push_back(nodes[j]);

		edge* newEdge = new edge();
		newEdge->src = nodes[i];
		newEdge->dest = nodes[j];
	}
	
void DAG::removeDependency(int i,int j){
		for(int i=0;i<edges.size();i++){
			if(edges[i]->src == nodes[i] && edges[i]->dest == nodes[j]){
				edges.erase(edges.begin()+i);
				break;
			}
		}
	}

void DAG::generateDAG(int no_of_tasks, int no_of_cores){

		cout<<"generate dag started\n";

		this->deadline = 0;

		for(int i=0;i<no_of_tasks;i++){
			task* newTask = new task();
			newTask->id = i;
			newTask->worst_case_time = rand()%(100-10 + 1) + 10;
			newTask->core_assigned = rand()%no_of_cores;

			nodes.push_back(newTask);
		}

		cout<<"generate dag mid\n";

		for(int i=0;i<no_of_tasks;i++){
			for(int j=i+1;j<no_of_tasks;j++){
				int r = rand()%2;
				if(r){
					nodes[j]->predecessors.push_back(nodes[i]);
					nodes[i]->successors.push_back(nodes[j]);

					edge* newEdge = new edge();
					newEdge->src = nodes[i];
					newEdge->dest = nodes[j];

					edges.push_back(newEdge);
				}
			}
		}

		cout<<"generate dag ended\n";

	}

void DAG::resetAssignment(){
	for(auto l:nodes){
		l->recovery_assigned = false;
	}
}

void DAG::inputDAG(string filepath, int no_of_cores){
	//read from file and create a DAG 
	// have filepath as input

	// cout<<"reading DAG from file "<<filepath<<"\n";
	this->filesname = filepath;
	fstream fin(filepath);
	int no_of_jobs;
	fin>>no_of_jobs;

	for(int i=0;i<no_of_jobs;i++){

		task* newTask = new task();
		newTask->id = i+1;
		newTask->worst_case_time = rand()%(100-10 + 1) + 10;
		newTask->core_assigned = rand()%no_of_cores;

		nodes.push_back(newTask);
	}

	for(int i=0;i<no_of_jobs;i++){
		int id;
		double runtime;
		fin >> id >> runtime;
		nodes[id]->worst_case_time = runtime;
		nodes[id]->execution_time = runtime;
	}

	int src,dst;
	while(fin >> src >> dst){

		nodes[dst]->predecessors.push_back(nodes[src]);
		nodes[src]->successors.push_back(nodes[dst]);

		edge* newEdge = new edge();
		newEdge->src = nodes[src];
		newEdge->dest = nodes[dst];

		edges.push_back(newEdge);
	}


}

void DAG::displayDAG(){
	cout<<"Nodes == Tasks are as follows\n";

	for(task*  l:nodes){
		cout<<"Task id : "<<l->id<<" wc time : "<<l->worst_case_time<<" core_assigned : "<<l->core_assigned<<" \n";
		cout<<"Predecessors : ";
		for(task* node: l->predecessors){
			cout<<node->id<<" ";
		}
		cout<<endl;

		cout<<"Successors : ";
		for(task* node: l->successors){
			cout<<node->id<<" ";
		}
		cout<<endl;
	}

	cout<<"END OF GRAPH\n";
}

void DAG::CopyDAG(DAG* new_dag){
	new_dag->deadline =  this->deadline;
	new_dag->ventry_id = this->ventry_id;
	new_dag->vexit_id = this->vexit_id;
	for(task* node: this->nodes){
		task* new_node = new task();
		*new_node = *node;
		new_dag->nodes.push_back(new_node);
	}

	for(edge* edgee: this->edges){
		edge* new_edge = new edge();
		*new_edge = *edgee;
		new_dag->edges.push_back(new_edge);
	}
}
