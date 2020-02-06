#include "canonical.h"

using namespace std;

int main(int argc,char *argv[]){
	parameterMap parmap;
	parameter::ReadParsFromFile(parmap,"parameters_Omega100.dat");
	CpartitionFunction pf(parmap);
	int ievent,nevents,i;
	double eta=10.0;
	Cmoments moments;
	moments.acceptance->ACCEPTANCE=1.0;
	CRandom *randy=new CRandom(1234);
	moments.ncalls=1000;
	moments.Qbar=eta*moments.ncalls;
	moments.Q2bar=(eta*eta+eta)*moments.ncalls;
	moments.Q3bar=(eta*eta*eta+3*eta*eta+eta)*moments.ncalls;
	moments.Q4bar=(eta*eta*eta*eta+6*eta*eta*eta+7*eta*eta+eta)*moments.ncalls;
	moments.Bbar=eta*moments.ncalls;
	moments.B2bar=(eta*eta+eta)*moments.ncalls;
	moments.B3bar=(eta*eta*eta+3*eta*eta+eta)*moments.ncalls;
	moments.B4bar=(eta*eta*eta*eta+6*eta*eta*eta+7*eta*eta+eta)*moments.ncalls;
	moments.Summarize(pf);
	return 0;
}
