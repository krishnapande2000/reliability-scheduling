#include <application_model.h>
#include <system_model.h>

double taskRealibility(task t){
	double ro = t.c.ro;
	double w = t.c.w;
	double Fmin = t.c.Fmin;
	double wc = t.worst_case_time;
	double freq = t.freq_assigned;

	double rf = ro*(power(10,(w*(1-freq))/(1-Fmin)));

	double Ri = power(e,(-1)*(rf*wc)/freq);

	return Ri;
}

double taskRecReliability(task t){
	double Ri = taskRealibility(t);
	// double Prec = ?

	double Rreci = 1 - (1-Ri)*(1-Ri);

	return Rreci;
}

double systemReliability(vector<task> tasks){
	double Rsys = 1;

	for(task t: tasks){
		double Ri;
		if(t.recovery_Assigned) Ri = taskRecReliability(task);
		else Ri = taskRealibility(task);
		Rsys = Rsys*Ri;
	}

	return Rsys;
}

double lifetimeReliability(core c){
	double MTTF =1;
	return MTTF;
	
}