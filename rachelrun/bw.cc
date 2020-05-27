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

	nevents=1000000; //parameter::getI(parmap,"NEVENTS",100);
	/*
	double pbound=4000;
	int ndp=50,nparts;
  double dp=pbound/double(ndp);
	double dN_dp_p2[ndp];
	for (int i=0;i<ndp;i++) {
		dN_dp_p2[i]=0;
	}
	double pi=3.1415926536;
	double hbarc=197.3269718;
	double pmag;
	int ibin,code;
	Cpart *ipart;
	*/
	Omega=parameter::getD(parmap,"OMEGA",100);
	string tag=parameter::getS(parmap,"FILE_TAG","_");
	//string file=parameter::getS(parmap,"MOMENTS_OUTPUT_FILE","moments.dat");
	//pf.SetOmega(Omega);

	for(iroots=0;iroots<NROOTS;iroots++){
		printf("____ roots=%lf ____ pf.T=%g ____ pf.Omega=%g ____\n",roots[iroots],pf.T,pf.Omega);
		pf.CalcZofOmega0(T[iroots]);
		pf.ScaleZ(Omega);
		printf("----------- Z Calculated -----------\n");

		file="data/STAR"+tag+".dat";
		altfile="altdata/STAR"+tag+".dat";
		//FILE *fptr=fopen(file.c_str(),"w");
		//strcat(file,const char(roots));
		//muBoverT=GetMuBOverT(roots[iroots]);
		ff=log(1.0+ffa*(roots[iroots]-roots[0])/(roots[NROOTS-1]-roots[0]))/log(1.0+ffa); // interpolating weight from zero to 1
		blastwave.Tf=120.0-20*ff;
		blastwave.uperpx=0.5+(0.74-0.5)*ff;
		blastwave.uperpy=blastwave.uperpx;

		for(int irun=0;irun<nruns;irun++){
			printf("run: %d\n",irun);
			moments.Clear();
			for(ievent=0;ievent<int(nevents/nruns);ievent++){

				do{
					b0=pf.randy->GetNPoissonian(rhoB[iroots]*Omega);
					q0=pf.randy->GetNPoissonian(0.5*rhoB[iroots]*Omega);
					if(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0)){
						printf("picked b0=%d or q0=%d out of bounds\n",b0,q0);
						printf("If this happens often, increase pf.NhadMAX\n");
					}
				}while(!pf.CheckRelevance(pf.NhadMAX/2,b0,q0,s0));

				//b0=0; //rhoB*Omega;
				//q0=0; //0.5*rhoB*Omega;
				pf.GenEvent(b0,q0,s0,resinfovec);
				blastwave.GenerateParts(resinfovec,partvec);
				moments.IncrementMoments(partvec);
				/*
				nparts=partvec.size();
				for (int i=0;i<nparts;i++) {
					ipart=&partvec[i];
					code=ipart->resinfo->code;
					if (code==211 || code==-211 || code==111) {
						pmag=sqrt(ipart->p[1]*ipart->p[1]+ipart->p[2]*ipart->p[2]+ipart->p[3]*ipart->p[3]);
		      	ibin=floorl(pmag/dp);
		      	if (ibin<ndp) dN_dp_p2[ibin]+=(2*pi*pi*hbarc*hbarc*hbarc)/(3*(ibin*dp+dp/2)*(ibin*dp+dp/2)*dp);
					}
				}
				*/
				partvec.clear();
			}
			moments.Summarize(file,altfile,Omega,rhoB[iroots],0.5*rhoB[iroots],roots[iroots],T[iroots]);
			moments.Clear();
		}
		/*
		for (ibin=0;ibin<ndp;ibin++) {
			fprintf(fptr,"%lf %lf\n",ibin*dp+dp/2,dN_dp_p2[ibin]/Omega);
		}
		fclose(fptr);
		*/
	}

	return 0;
}
