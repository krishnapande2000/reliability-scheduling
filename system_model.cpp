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
		newCore->hardware_coefficient = generateRand1(30,10);
		newCore->aging_rate = generateRand1(45,5);
		newCore->weibull_alpha = generateRand1(60,10);
		newCore->Fmax = rand()%(300-100 + 1) + 100;
		newCore->Fmin = rand()%((int)newCore->Fmax-100 + 1) + 100;
		newCore->ro = generateRand1(54,3);
		newCore->id = i;

		cores.push_back(newCore);

	}
}
