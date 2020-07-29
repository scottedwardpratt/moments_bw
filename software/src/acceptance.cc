#include "canonical.h"

using namespace std;

void Cacceptance::Acceptance(Cpart *part,bool &acceptQ,bool &acceptP,bool &acceptK,bool &acceptPi,bool &acceptB){
	bool STAR = true;
	bool EFFICIENCY_CORRECTED=true;
	double r=randy->ran();
	acceptQ=acceptP=acceptK=acceptPi=true;
	double efficiency;

	if (STAR) {
		double y,eta,pt,pmag,ptmin=200.0,ptmax=2000.0;
		
		CResInfo *resinfo=part->resinfo;
		pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
		pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
		eta=atanh(part->p[3]/pmag);
		y=atanh(part->p[3]/part->p[0]);
		int pid=part->pid;
		if(abs(pid)==211 || abs(pid)=311){
			ptmin=200.0;
			ptmax=1600;
		}
		if(abs(pid)==2212){
			ptmin=400.0;
			ptmax=2000.0;
		}

		// For net charge
		ptmin=200.0;
		ptmax=2000.0;
		efficiency=0.715;
		if(EFFICIENCY_CORRECTED)
			efficiency=1.0;
		if(resinfo->charge==0 || fabs(eta)>0.5 || pt<ptmin || pt>ptmax || r>efficiency){
			acceptQ=acceptP=acceptK=acceptPi=false;
		}
		
		// For specific species
		if(abs(pid)==2212){
			ptmin=400.0;
			ptmax=2000.0;
			if(fabs(eta)>0.5 || p)
		}
		
		else{
			acceptP=acceptK=acceptPi=true;
			if(pt>ptmax || pt<ptmin || fabs(eta)>0.5 || r>0.6){
				acceptP=acceptK=acceptPi=false;
			}
			acceptQ=true;
			if(pt>2000.0 || pt<200 || fabs(eta)>0.5 || r>0.8){
				acceptQ=false;
			}
		}
	}
	else {
		if(r>ACCEPTANCE)
			acceptQ=acceptP=acceptK=acceptPi=acceptB=false;
		else acceptQ=acceptP=acceptK=acceptPi=acceptB=true;
	}
}
