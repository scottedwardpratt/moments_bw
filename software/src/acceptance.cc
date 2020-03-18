#include "canonical.h"

using namespace std;

void Cacceptance::Acceptance(Cpart *part,bool &acceptQ,bool &acceptP,bool &acceptK,bool &acceptPi){
	double y,eta,pt,pmag;
	CResInfo *resinfo=part->resinfo;
	pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
	pmag=sqrt(pt*pt+part->p[3]*part->p[3]);
	eta=atanh(part->p[3]/pmag);
	y=atanh(part->p[3]/part->p[0]);
	double r=randy->ran();

	if(r>ACCEPTANCE) acceptQ=acceptP=acceptK=acceptPi=false;
	else acceptQ=acceptP=acceptK=acceptPi=true;
	/*
	if(resinfo->charge==0 || fabs(eta)>1.0 || pt<150.0 || r>0.8){
		acceptQ=acceptP=acceptK=acceptPi=false;
	}
	else{
		acceptQ=acceptP=acceptK=acceptPi=true;
		if(pt>1500 || fabs(eta)>0.9 || r>0.6){
			acceptP=acceptK=acceptPi=false;
		}
	}
	*/
}
