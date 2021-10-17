struct core{
	int id;
	double aging_rate;
	double weibull_alpha;
	vector<double> freq_levels;
	double Fmax;
	double Fmin;
	double ro;

	vector<task*> tasks;
};

class Multicore{
	vector<core*> cores;

	void Multicore();

	void generateCores(no_of_cores);
};