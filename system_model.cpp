#include <bits/stdc++.h>
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

		cores.push_back(newCore);

	}
}
