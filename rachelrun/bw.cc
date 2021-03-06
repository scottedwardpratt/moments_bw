#include "canonical.h"
//#include "MuBOverT.cc"
#include <string>
using namespace std;

int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: bw parameter_filename\n");
		exit(1);
	}
	const int NROOTS=7;
	int ievent,nevents,i,b0=0,q0=0,s0=0,iroots; // initial charge, baryon no. and strangeness
	int nruns=10;
	string file;
	string altfile;
	bool fixedcharge;
	double Omega,muBoverT,ff,ffa=1.0;
	parameterMap parmap;
	parameter::ReadParsFromFile(parmap,string(argv[1]));
	//double T=parameter::getD(parmap,"T",150.0);
	//double rhoB=parameter::getD(parmap,"RHOB",0.05);
	double roots[NROOTS]={7.7,11.5,19.6,27.0,39.0,62.4,200.0};
	// These correspond to rhoH=0.333333333
	double T[NROOTS]={143.45,149.53,153.4,154.45,155.08,155.44,155.66};
	double rhoB[NROOTS]={0.143511,0.0942362,0.0525105,0.0373777,0.0255775,0.0159063,0.00497411};
	//

	Cacceptance acceptance(parameter::getD(parmap,"ACCEPTANCE",0.25));
	Cmoments moments(&acceptance);
	vector<CResInfo *> resinfovec;
	vector<Cpart> partvec;
	CpartitionFunction pf(parmap);
	pf.DECAY=false;  // don't decay particles (will be done by blastwave)
	CblastWave blastwave(parmap,pf.randy,pf.reslist);
	//pf.randy->reset(-time(NULL));
	pf.randy->reset(1234);
	//pf.WriteZ();
	//pf.ReadZ();

	nevents=10000; //parameter::getI(parmap,"NEVENTS",100);
	Omega=parameter::getD(parmap,"OMEGA",100);
	string tag=parameter::getS(parmap,"FILE_TAG","_");
	//pf.SetOmega(Omega);
	if(parameter::getB(parmap,"FIXEDCHARGE",true))
            fixedcharge=true;
    else
            fixedcharge=false;
    if(fixedcharge) printf("fixed charge\n");
	for(iroots=0;iroots<NROOTS;iroots++){
		printf("____ roots=%lf ____ pf.T=%g ____ pf.Omega=%g ____\n",roots[iroots],pf.T,pf.Omega);
		pf.CalcZofOmega0(T[iroots]);
		pf.ScaleZ(Omega);
		printf("----------- Z Calculated -----------\n");

		if (fixedcharge){
			file="data/"+tag+".dat";
			altfile="altdata/"+tag+".dat";
		}
		else {
			file="data/"+tag+"poissonQ.dat";
			altfile="altdata/"+tag+"poissonQ.dat";
		}

		ff=log(1.0+ffa*(roots[iroots]-roots[0])/(roots[NROOTS-1]-roots[0]))/log(1.0+ffa); // interpolating weight from zero to 1
		blastwave.Tf=120.0-20*ff;
		blastwave.SetYbeam(roots[iroots]);
		printf("sigma_source=%g\n",blastwave.sigma_source);
		blastwave.uperpx=0.5+(0.74-0.5)*ff;
		blastwave.uperpy=blastwave.uperpx;
		for(int irun=0;irun<nruns;irun++){
			printf("run: %d\n",irun);
			moments.Clear();
			for(ievent=0;ievent<int(nevents/nruns);ievent++){
				if(fixedcharge){
                    b0=Omega*rhoB[iroots]; //parameter::getD(parmap,"RHOB",0.0);
                    q0=b0/2.0;
                    //if(ievent==0 && irun==0) printf("b0=%d\n",b0);
                }
                else{
                    do{
                        b0=pf.randy->GetNPoissonian(rhoB[iroots]*Omega); //(parameter::getD(parmap,"RHOB",0.0)*Omega);
                        q0=pf.randy->GetNPoissonian(0.5*rhoB[iroots]*Omega); //(0.5*parameter::getD(parmap,"RHOB",0.0)*Omega);
                        if(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0)){
                            printf("picked b0=%d or q0=%d out of bounds\n",b0,q0);
                            printf("If this happens often, increase pf.NhadMAX\n");
                        }
                    }while(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0));
                }
				pf.GenEvent(b0,q0,s0,resinfovec);
				blastwave.GenerateParts(resinfovec,partvec);
				moments.IncrementMoments(partvec);
				partvec.clear();
			}
			moments.Summarize(file,altfile,Omega,rhoB[iroots],0.5*rhoB[iroots],roots[iroots],T[iroots]);
			moments.Clear();
		}
	}

	return 0;
}
