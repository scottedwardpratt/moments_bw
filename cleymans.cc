#include "canonical.h"
#include "MuBOverT.cc"
using namespace std;

int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: cleymans parameter_filename\n");
		exit(1);
	}
	int b0,q0,s0,iroots,ib,iq,is; // initial charge, baryon no. and strangeness
	double muBoverT,T,rhoN,rhoH,rhoB,rhoHtarget;
	double roots[7]={7.7,11.5,19.6,27.0,39.0,62.4,200.0};
	parameterMap parmap;
	parameter::ReadParsFromFile(parmap,string(argv[1]));
	CpartitionFunction pf(parmap);
	pf.randy->reset(-time(NULL));
	/*
	pf.T=165.0;
	muBoverT=0.0;
	pf.CalcZ();
	rhoN=0.0;  // nucleon density
	for(q0=-2;q0<=2;q0++){
		for(s0=-3;s0<=3;s0++){
			if(pf.CheckRelevance(1,1,q0,s0)){
				b0=-1;
				pf.Getibiqis(1,1,q0,s0,ib,iq,is);
				rhoN+=pf.Z[1][ib][iq][iq]/pf.Omega0;
				b0=1;
				pf.Getibiqis(1,1,q0,s0,ib,iq,is);
				rhoN+=pf.Z[1][ib][iq][is]/pf.Omega0;
			}
		}
	}

	rhoH=0.0; // hadron density
	for(b0=-1;b0<=1;b0++){
		for(q0=-2;q0<=2;q0++){
			for(s0=-3;s0<=3;s0++){
				if(pf.CheckRelevance(1,b0,q0,s0)){
					pf.Getibiqis(1,b0,q0,s0,ib,iq,is);
					rhoH+=pf.Z[1][ib][iq][is]/pf.Omega0;
				}
			}
		}
	}
	//printf("T=%g: With muB=0, rhoN=%g, rhoH=%g\n",T,rhoN,rhoH);
	rhoH=rhoH+rhoN*(cosh(muBoverT)-1.0);
	rhoB=rhoN*sinh(muBoverT);
	printf("For T=%g, rhoH=%g, rhoB=%g\n",pf.T,rhoH,rhoB);
	exit(1);
	*/
	printf("Enter rhoHtarget: ");
	scanf("%lf",&rhoHtarget);
	T=130.0;
	for(iroots=0;iroots<7;iroots++){
		printf("__________ ROOTS=%g ____________\n",roots[iroots]);
		muBoverT=GetMuBOverT(roots[iroots]);
		do{
			T+=0.01;
			pf.T=T;
			pf.CalcZ();
			rhoN=0.0;  // nucleon density
			for(q0=-2;q0<=2;q0++){
				for(s0=-3;s0<=3;s0++){
					if(pf.CheckRelevance(1,1,q0,s0)){
						b0=-1;
						pf.Getibiqis(1,1,q0,s0,ib,iq,is);
						rhoN+=pf.Z[1][ib][iq][iq]/pf.Omega0;
						b0=1;
						pf.Getibiqis(1,1,q0,s0,ib,iq,is);
						rhoN+=pf.Z[1][ib][iq][is]/pf.Omega0;
					}
				}
			}
	
			rhoH=0.0; // hadron density
			for(b0=-1;b0<=1;b0++){
				for(q0=-2;q0<=2;q0++){
					for(s0=-3;s0<=3;s0++){
						if(pf.CheckRelevance(1,b0,q0,s0)){
							pf.Getibiqis(1,b0,q0,s0,ib,iq,is);
							rhoH+=pf.Z[1][ib][iq][is]/pf.Omega0;
						}
					}
				}
			}
			//printf("T=%g: With muB=0, rhoN=%g, rhoH=%g\n",T,rhoN,rhoH);
			rhoH=rhoH+rhoN*(cosh(muBoverT)-1.0);
			rhoB=rhoN*sinh(muBoverT);
		}while(rhoH<rhoHtarget);
		printf("T=%g, muB/T=%g, pbar/p=%g, rhoB=%g, rhoH=%g\n",T,muBoverT,exp(-2.0*muBoverT),rhoB,rhoH);
	}

	
	
	return 0;
}
