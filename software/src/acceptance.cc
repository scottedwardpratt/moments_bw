#include "canonical.h"

using namespace std;

void Cacceptance::Acceptance(Cpart *part,bool &acceptQ,bool &acceptP,bool &acceptK,bool &acceptPi,bool &acceptB){
	bool STAR = true;
	bool EFFICIENCY_CORRECTED=true;
	double r=randy->ran();
	acceptQ=acceptP=acceptK=acceptPi=false;
	double efficiency;
	double pt,pmag,eta,y,ptmin,ptmax;
	CResInfo *resinfo;
	resinfo=part->resinfo;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	y=atanh(part->p[3]/part->p[0]);
	int pid=part->resinfo->code;
	
	if (STAR) {
		// For net charge
		if(resinfo->charge!=0){
			acceptQ=false;
			ptmin=200.0;
			ptmax=2000.0;
			EFFICIENCY_CORRECTED=true;
			efficiency=0.715;
			if(EFFICIENCY_CORRECTED)
				efficiency=1.0;
			if(fabs(eta)<0.5 && pt>ptmin && pt<ptmax && r<efficiency){
				acceptQ=true;
			}
		}
		
		// For specific species
		if(abs(pid)==211){
			acceptPi=false;
			ptmin=200.0;
			ptmax=1600.0;
			efficiency=0.6;
			EFFICIENCY_CORRECTED=true;
			if(EFFICIENCY_CORRECTED)
				efficiency=1.0;
			if(fabs(y)<0.5 && pt>ptmin && pt<ptmax && r<efficiency){
				acceptPi=true;
			}	
		}
		else if(abs(pid)==2212){
			acceptP=false;
			ptmin=400.0;
			ptmax=2000.0;
			efficiency=0.6;
			EFFICIENCY_CORRECTED=true;
			if(EFFICIENCY_CORRECTED)
				efficiency=1.0;
			if(fabs(y)<0.5 && pt>ptmin && pt<ptmax && r<efficiency){
				acceptP=true;
			}	
		}
		else if(abs(pid)==321){
			acceptK=false;
			ptmin=200.0;
			ptmax=1600.0;
			efficiency=0.6;
			EFFICIENCY_CORRECTED=true;
			if(EFFICIENCY_CORRECTED)
				efficiency=1.0;
			if(fabs(y)<0.5 && pt>ptmin && pt<ptmax && r<efficiency){
				acceptK=true;
			}	
		}
	}
	else{
		if(resinfo->charge!=0){
			acceptQ=false;
			ptmin=0.0;
			ptmax=20000000.0;
			EFFICIENCY_CORRECTED=true;
			efficiency=1.0;
			if(EFFICIENCY_CORRECTED)
				efficiency=1.0;
			if(fabs(eta)<1.0 && pt>ptmin && pt<ptmax & r<efficiency){
				acceptQ=acceptP=acceptK=acceptPi=false;
			}
		}	
	}
}
