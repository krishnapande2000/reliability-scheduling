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
	double end_time;
	double execution_time;
	double dynamic_ready_time;
	double dynamic_process_start_time;

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
	void resetAssignment();
	void displayDAG();
	void inputDAG(string filepath, int no_of_cores);

};