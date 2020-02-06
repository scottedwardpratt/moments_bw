#include "canonical.h"

using namespace std;

int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: canonical parameter_filename\n");
		exit(1);
	}
	int ievent,nevents,i,b0,q0,s0=0; // initial charge, baryon no. and strangeness
	double Omega,rhoB,rhoQ;
	parameterMap parmap;
	Cacceptance acceptance(parameter::getD(parmap,"ACCEPTANCE",0.25));
	Cmoments moments(&acceptance);
	vector<CResInfo *> resinfo;
	parameter::ReadParsFromFile(parmap,string(argv[1]));
	CpartitionFunction pf(parmap);
	pf.randy->reset(-time(NULL));
	pf.CalcZ();
	//pf.WriteZ();
	//pf.ReadZ();
	
	nevents=parameter::getI(parmap,"NEVENTS",100);
	Omega=parameter::getD(parmap,"OMEGA",100);
	rhoB=parameter::getD(parmap,"RHOB",0.0);
	
	rhoQ=0.5*rhoB;
	pf.SetOmega(Omega);
	for(ievent=0;ievent<nevents;ievent++){
		b0=pf.randy->GetNPoissonian(rhoB*Omega);
		q0=pf.randy->GetNPoissonian(rhoQ*Omega);
		pf.GenEvent(b0,q0,s0,resinfo);
		moments.IncrementMoments(resinfo);
		if((ievent+1)%(nevents/10)==0)
			printf("Finished %g percent\n",100.0*(ievent+1.0)/double(nevents));
	}
	moments.Summarize(Omega,rhoB,rhoQ,0.0,pf.T);
	moments.Clear();
	
	return 0;
}
