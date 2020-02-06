#include "canonical.h"
using namespace std;

int main(int argc,char *argv[]){
	parameterMap parmap;
	if(argc!=2)
		parameter::ReadParsFromFile(parmap,"parameters.dat");
	else
		parameter::ReadParsFromFile(parmap,string(argv[1]));
	int ievent,nevents,i;
	int b0,q0,s0=0; // initial charge, baryon no. and strangeness
	double Omega,rhoB,rhoQ;
	Cmultlist multlist;
	vector<CResInfo *> resinfo;
	CpartitionFunction pf(parmap);
	pf.randy->reset(-time(NULL));
	pf.CalcZ();
	//pf.WriteZ();
	//pf.ReadZ();
	//nevents=100000;
	nevents=parameter::getI(parmap,"NEVENTS",100);
	
	Omega=50.0;
	for(Omega=50.0;Omega<500;Omega*=2){
		pf.SetOmega(Omega);
		for(rhoB=0.01;rhoB<0.035;rhoB+=0.01){
			rhoQ=0.5*rhoB;					
			printf("----- Omega=%g ------ rhoB=%g ------ rhoQ=%g ------\n",Omega,rhoB,rhoQ);
			for(ievent=0;ievent<nevents;ievent++){
				//b0=lrint(rhoB*Omega);
				//q0=lrint(rhoQ*Omega);
				b0=pf.randy->GetNPoissonian(rhoB*Omega);
				q0=pf.randy->GetNPoissonian(rhoQ*Omega);
				pf.GenEvent(b0,q0,s0,resinfo);
				multlist.AddMult(resinfo);
				if((ievent+1)%(nevents/10)==0)
					printf("Finished %g percent\n",100.0*(ievent+1.0)/double(nevents));
			}
			multlist.WriteMults(Omega,rhoB,rhoQ);
		}
	}
	return 0;
}
