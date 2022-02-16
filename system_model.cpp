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
#include "system_model.h"
using namespace std;

double generateRand1(int n, int m) {
    double a = rand()%n;
    double b = rand()%m;
    return a/b;
}

/*
double aging_rate;
	double weibull_alpha;
	vector<double> freq_levels;
	double Fmax;
	double Fmin;
	double ro;
*/

Multicore::Multicore(){

}

void Multicore::generateCores(int no_of_cores){

	for(int i=0;i<no_of_cores;i++){
		Core* newCore = new Core();
		newCore->hardware_coefficient = 3;
		newCore->aging_rate = generateRand1(45,5);
		newCore->weibull_alpha = generateRand1(60,10);
		newCore->Fmax = 0.99;
		newCore->Fmin = 0.01;
		newCore->ro = 1e-5;
		newCore->id = i;
		newCore->free_at = INT_MAX;

		cores.push_back(newCore);

	}
}

void Multicore::clear_cores(){
	for(auto core:cores){
		core->tasks.clear();
	}
}
