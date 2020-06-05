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
	int ievent,nevents,b0=0,q0=0,s0=0,iroots; // initial charge, baryon no. and strangeness
	int nruns=1;
	string filename,tag,altfilename,parfiletag,parfilename;
	//double muBoverT;
	double ff,ffa=1.0;
	bool fixedcharge;
	parameterMap parmap;
	parfiletag=string(argv[1]);
	parfilename="parameters_"+parfiletag+".dat";
	parameter::ReadParsFromFile(parmap,parfilename);
	double Omega=parameter::getD(parmap,"OMEGA",140.0);
	double roots[NROOTS]={7.7,11.5,19.6,27.0,39.0,62.4,200.0};
	if(parameter::getB(parmap,"FIXEDCHARGE",true))
		fixedcharge=true;
	else
		fixedcharge=false;

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
	pf.randy->reset(time(NULL));
	//pf.WriteZ();
	//pf.ReadZ();

	nevents=parameter::getI(parmap,"NEVENTS",100);
	if(fixedcharge)
		tag=parameter::getS(parmap,"FILE_TAG","fixedcharge_")+parfiletag;
	else
		tag=parameter::getS(parmap,"FILE_TAG","randomcharge_")+parfiletag;
		
	//for(iroots=0;iroots<NROOTS;iroots++){
	for(iroots=0;iroots<NROOTS;iroots++){
		printf("____ iroots=%d ____ roots=%g ____ pf.T=%g ____ pf.Omega=%g ____\n",iroots,roots[iroots],pf.T,pf.Omega);
		pf.CalcZofOmega0(T[iroots]);
		pf.ScaleZ(Omega);
		printf("----------- Z Calculated -----------\n");

		filename="results/roots"+to_string(iroots)+"_"+tag+".dat";
		altfilename="results/alt_roots"+to_string(iroots)+"_"+tag+".dat";
		//strcat(file,const char(iroots));
		//muBoverT=GetMuBOverT(roots[iroots]);
		ff=log(1.0+ffa*(roots[iroots]-roots[0])/(roots[NROOTS-1]-roots[0]))/log(1.0+ffa); // interpolating weight from zero to 1
		blastwave.Tf=120.0-20*ff;
		//blastwave.Tf=120.0;
		blastwave.SetYbeam(roots[iroots]);
		printf("sigma_source=%g\n",blastwave.sigma_source);
		blastwave.uperpx=0.5+(0.74-0.5)*ff;
		blastwave.uperpy=blastwave.uperpx;
		for(int irun=0;irun<nruns;irun++){
			moments.Clear();
			for(ievent=0;ievent<int(nevents/nruns);ievent++){
				// poissonian fluctuation of net B and net Q=B/2
				do{
					if(fixedcharge){
						
					}
					else{
						b0=pf.randy->GetNPoissonian(rhoB[iroots]*Omega);
						q0=pf.randy->GetNPoissonian(0.5*rhoB[iroots]*Omega);
					}
					if(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0)){
						printf("picked b0=%d or q0=%d out of bounds\n",b0,q0);
						printf("If this happens often, increase pf.NhadMAX\n");
					}
				}while(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0));
				
				pf.GenEvent(b0,q0,s0,resinfovec);
				blastwave.GenerateParts(resinfovec,partvec);
				moments.IncrementMoments(partvec);
				partvec.clear();
				if(nevents>10 && (ievent+1)%(nevents/(10*nruns))==0)
					printf("Finished %g percent\n",100.0*(ievent+1.0)/double(nevents));
			}
			moments.Summarize(filename,altfilename,Omega,rhoB[iroots],0.5*rhoB[iroots],roots[iroots],T[iroots]);
			moments.Clear();
		}
	}

	return 0;
}
