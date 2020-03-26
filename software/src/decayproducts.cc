#include "canonical.h"

void CpartitionFunction::GetProducts(CResInfo *resinfo0,int &ndaughters,array<CResInfo *,10> &daughter){
	bool decay;
	//resinfo0->decay=false; //this turns off decays
	decay=resinfo0->decay;
	int idaughtertemp;
	int nmothers,imother,ndaughterstemp;
	array <CResInfo *,5> mother;
	array <CResInfo *,5> daughtertemp;
	/** Decay the i-particles */
	if(resinfo0->decay){
		mother[0]=resinfo0;
		nmothers=1;
		ndaughters=0;
	}
	else{
		ndaughters=1;
		nmothers=0;
		daughter[0]=resinfo0;
	}
	imother=0;
	while(imother<nmothers){
		mother[imother]->DecayGetResInfoPtr(ndaughterstemp,daughtertemp);
		for(idaughtertemp=0;idaughtertemp<ndaughterstemp;idaughtertemp++){
			if(daughtertemp[idaughtertemp]->decay){
				mother[nmothers]=daughtertemp[idaughtertemp];
				nmothers+=1;
			}
			else{
				daughter[ndaughters]=daughtertemp[idaughtertemp];
				ndaughters+=1;
			}
		}
		imother+=1;
	}
}