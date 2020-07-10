#include "canonical.h"
//#include "MuBOverT.cc"
#include <string>
using namespace std;

//int main(int argc,char *argv[]){
int main(int argc,char *argv[]){
	if(argc!=2){
		printf("Usage: bw parameter_filename\n");
		exit(1);
	}
	const int NROOTS=7,NRUNS=8;
	char filename[100];
	int b0=0,q0=0,s0=0,iroots,irun; // initial charge, baryon no. and strangeness
	int ipatch_eta,ipatch_uperp,ievent,nevents;
	//double muBoverT;
	double ff,ffa=1.0;
	bool fixedcharge;
	parameterMap parmap;
	string parfilename="parameters/parameters_"+string(argv[1])+".dat";
	parameter::ReadParsFromFile(parmap,parfilename);
	double roots[NROOTS]={7.7,11.5,19.6,27.0,39.0,62.4,200.0};
	double Omega[NRUNS]={25.0,50.0,100.0,200.0,400.0,100.0,100.0,100.0};
	double sigmaeta[NRUNS]={0.3,0.3,0.3,0.3,0.3,0.1,0.5,0.7};
	nevents=parameter::getI(parmap,"NEVENTS",1000);
	fixedcharge=parameter::getB(parmap,"FIXEDCHARGE",false);
	if(fixedcharge){
		printf("This program does not work for fixed charge\n");
		exit(1);
	}

	// These correspond to rhoH=0.333333333
	double T[NROOTS]={143.45,149.53,153.4,154.45,155.08,155.44,155.66};
	double rhoB[NROOTS]={0.143511,0.0942362,0.0525105,0.0373777,0.0255775,0.0159063,0.00497411};
	//
	CpartitionFunction pf(parmap);
	pf.DECAY=false;  // don't decay particles (will be done by blastwave)
	pf.randy->reset(time(NULL));
	
	Cacceptance acceptance(parameter::getD(parmap,"ACCEPTANCE",0.25));
	
	CblastWave blastwave(parmap,pf.randy,pf.reslist);
	
	vector<vector<Cmoments*>> moments;
	moments.resize(blastwave.NPATCHES_ETA);
	for(ipatch_eta=0;ipatch_eta<blastwave.NPATCHES_ETA;ipatch_eta++){
		moments[ipatch_eta].resize(blastwave.NPATCHES_UPERP);
		for(ipatch_uperp=0;ipatch_uperp<blastwave.NPATCHES_UPERP;ipatch_uperp++){
			moments[ipatch_eta][ipatch_uperp]=new Cmoments(&acceptance);
		}
	}
	Cmoments moments_sum=(&acceptance);

	vector<CResInfo *> resinfovec;
	vector<Cpart> partvec;

	//pf.WriteZ();
	//pf.ReadZ();

	for(iroots=0;iroots<NROOTS;iroots++){
		printf("____ iroots=%d ____ roots=%g ____ pf.T=%g ____\n",iroots,roots[iroots],T[iroots]);
		printf("NPATCHES_ETA=%d, NPATCHES_UPERP=%d\n",blastwave.NPATCHES_ETA,blastwave.NPATCHES_UPERP);
		blastwave.Tf=120.0-20*ff;
		blastwave.SetYbeam(roots[iroots]);
		ff=log(1.0+ffa*(roots[iroots]-roots[0])/(roots[NROOTS-1]-roots[0]))/log(1.0+ffa); // interpolating weight from zero to 1
		blastwave.uperpx=0.5+(0.74-0.5)*ff;
		blastwave.uperpy=blastwave.uperpx;
		
		pf.CalcZofOmega0(T[iroots]);
		printf("----------- Z Calculated, T=%g -----------\n",T[iroots]);
		
		// Calc for different Omega
		for(irun=0;irun<NRUNS;irun++){
			blastwave.sigma_eta=sigmaeta[irun];
			pf.ScaleZ(Omega[irun]);
			sprintf(filename,"results/sigmaeta%g_Omega%g.dat",sigmaeta[irun],Omega[irun]);
			
			for(ipatch_eta=0;ipatch_eta<blastwave.NPATCHES_ETA;ipatch_eta++){
				for(ipatch_uperp=0;ipatch_uperp<blastwave.NPATCHES_UPERP;ipatch_uperp++){
					for(ievent=0;ievent<nevents;ievent++){
						moments[ipatch_eta][ipatch_uperp]->ResetEvent();
						// poissonian fluctuation of net B and net Q=B/2
						do{
							b0=pf.randy->GetNPoissonian(rhoB[iroots]*Omega[irun]);
							q0=pf.randy->GetNPoissonian(0.5*rhoB[iroots]*Omega[irun]);
							if(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0)){
								printf("picked b0=%d or q0=%d out of bounds\n",b0,q0);
								printf("If this happens often, increase pf.NhadMAX\n");
							}
						}while(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0));
						pf.GenEvent(b0,q0,s0,resinfovec);
						blastwave.GenerateParts(resinfovec,partvec);
						blastwave.BoostParts(partvec,ipatch_eta,ipatch_uperp);
						moments[ipatch_eta][ipatch_uperp]->IncrementMoments(partvec);
						partvec.clear();
					}
					moments[ipatch_eta][ipatch_uperp]->CalcCumulants();
				}
			}
			moments_sum.Sum(moments);
			for(ipatch_eta=0;ipatch_eta<blastwave.NPATCHES_ETA;ipatch_eta++){
				for(ipatch_uperp=0;ipatch_uperp<blastwave.NPATCHES_UPERP;ipatch_uperp++){
					moments[ipatch_eta][ipatch_uperp]->Clear();
				}
			}
			moments_sum.Summarize(string(filename),Omega[irun],rhoB[iroots],0.5*rhoB[iroots],roots[iroots],T[iroots]);
			moments_sum.Clear();
		}	
	
	}
	return 0;
}
