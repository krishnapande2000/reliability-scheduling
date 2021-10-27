#include <bits/stdc++.h>
#include "application_model.h"
using namespace std;

struct Core{
	int id;
	double hardware_coefficient;
	double aging_rate;
	double weibull_alpha;
	vector<double> freq_levels;
	double Fmax;
	double Fmin;
	double ro;

	double total_time;

	vector<task*> tasks;

	Core(){
		total_time=0;
	}
};

class Multicore{
public:
	vector<Core*> cores;

	Multicore();

	void generateCores(int no_of_cores);
};