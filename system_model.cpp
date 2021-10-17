#include "system_model.h"

double generateRand(int n, int m) {
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
class Multicore{
	vector<core*> cores;

	void Multicore(){

	}

	void generateCores(no_of_cores){

		for(int i=0;i<no_of_cores;i++){
			core* newCore;
			core->aging_rate = generateRand(45,5);
			core->weibull_alpha = generateRand(60,10);
			core->Fmax = generateRand(500,4);
			core->Fmin = generateRand((int)core->Fmax, 3);
			core->ro = generateRand(54,3);
			core->id = i;

			cores.push_back(newCore);

		}
	}
};