#include "canonical.h"

using namespace std;

Cmoments::Cmoments(Cacceptance *acceptanceset){
	acceptance=acceptanceset;
	Clear();
}

void Cmoments::Clear(){
	ncalls=0.0;
	Qbar=Q2bar=Q3bar=Q4bar=Pbar=P2bar=P3bar=P4bar=Kbar=K2bar=K3bar=K4bar=Pibar=Pi2bar=Pi3bar=Pi4bar=0.0;
	TotQbar=TotPbar=TotKbar=TotPibar=0.0;
	meanpt_pions=0.0;
	meanpt_kaons=0.0;
	meanpt_protons=0.0;
}

void Cmoments::IncrementMoments(vector<CResInfo *> &resinfovec){
	ncalls+=1;
	int NetQ=0,NetP=0,NetK=0,NetPi=0;
	int NetQ2=0,NetP2=0,NetK2=0,NetPi2=0;
	int NetQ3=0,NetP3=0,NetK3=0,NetPi3=0;
	int NetQ4=0,NetP4=0,NetK4=0,NetPi4=0;
	int TotP=0,TotQ=0,TotK=0,TotPi=0;
	int i,mult=resinfovec.size();
	for(i=0;i<mult;i++){
		if(acceptance->Acceptance(resinfovec[i])){
			if(abs(resinfovec[i]->charge)==1){
				NetQ+=resinfovec[i]->charge;
				TotQ+=1;
			}
			if(abs(resinfovec[i]->code)==2212){
				NetP+=resinfovec[i]->charge;
				TotP+=1;
			}
			if(abs(resinfovec[i]->code)==321){
				NetK+=resinfovec[i]->charge;
				TotK+=1;
			}
			if(abs(resinfovec[i]->code)==211){
				NetPi+=resinfovec[i]->charge;
				TotPi+=1;
			}
		}
	}
	Qbar+=NetQ; Q2bar+=NetQ*NetQ; Q3bar+=pow(NetQ,3); Q4bar+=pow(NetQ,4);
	Pbar+=NetP; P2bar+=NetP*NetP; P3bar+=pow(NetP,3); P4bar+=pow(NetP,4);
	Kbar+=NetK; K2bar+=NetK*NetK; K3bar+=pow(NetK,3); K4bar+=pow(NetK,4);
	Pibar+=NetPi; Pi2bar+=NetPi*NetPi; Pi3bar+=pow(NetPi,3); Pi4bar+=pow(NetPi,4);

	TotQbar+=TotQ; TotPbar+=TotP; TotKbar+=TotK; TotPibar+=TotPi;
}

void Cmoments::IncrementMoments(vector<Cpart> &partvec){
	bool acceptQ,acceptP,acceptK,acceptPi;
	double pt;
	CResInfo *resinfo;
	ncalls+=1;
	int NetQ=0,NetP=0,NetK=0,NetPi=0;
	int NetQ2=0,NetP2=0,NetK2=0,NetPi2=0;
	int NetQ3=0,NetP3=0,NetK3=0,NetPi3=0;
	int NetQ4=0,NetP4=0,NetK4=0,NetPi4=0;
	int TotQ=0,TotP=0,TotK=0,TotPi=0;
	int i;
	for(int ipart=0;ipart<partvec.size();ipart++){
		resinfo=partvec[ipart].resinfo;
		pt=sqrt(partvec[ipart].p[1]*partvec[ipart].p[1]+partvec[ipart].p[2]*partvec[ipart].p[2]);
		acceptance->Acceptance(&partvec[ipart],acceptQ,acceptP,acceptK,acceptPi);
		//acceptQ=acceptP=acceptK=acceptPi=true;
		if(abs(resinfo->charge)==1 && acceptQ){
			NetQ+=resinfo->charge;
			TotQ+=1;
		}
		if(abs(resinfo->code)==2212 && acceptP){
			NetP+=resinfo->charge;
			TotP+=1;
			meanpt_protons+=pt;
		}
		if(abs(resinfo->code)==321 && acceptK){
			NetK+=resinfo->charge;
			TotK+=1;
			meanpt_kaons+=pt;
		}
		if(abs(resinfo->code)==211 && acceptPi){
			NetPi+=resinfo->charge;
			TotPi+=1;
			meanpt_pions+=pt;
		}
	}
	Qbar+=NetQ; Q2bar+=NetQ*NetQ; Q3bar+=pow(NetQ,3); Q4bar+=pow(NetQ,4);
	Pbar+=NetP; P2bar+=NetP*NetP; P3bar+=pow(NetP,3); P4bar+=pow(NetP,4);
	Kbar+=NetK; K2bar+=NetK*NetK; K3bar+=pow(NetK,3); K4bar+=pow(NetK,4);
	Pibar+=NetPi; Pi2bar+=NetPi*NetPi; Pi3bar+=pow(NetPi,3); Pi4bar+=pow(NetPi,4);

	TotQbar+=TotQ; TotPbar+=TotP; TotKbar+=TotK; TotPibar+=TotPi;
}

void Cmoments::Summarize(string file,double Omega,double rhoB,double rhoQ,double roots,double T){
	double qbar,kappaq2,kappaq3,kappaq4;
	double pbar,kappap2,kappap3,kappap4;
	double kbar,kappak2,kappak3,kappak4;
	double pibar,kappapi2,kappapi3,kappapi4;
	double totq,totp,totk,totpi;
	double sigma2q,sigma2p,sigma2k,sigma2pi,Ssigmaq,Ssigmap,Ssigmak,Ssigmapi,Ksigma2q,Ksigma2p,Ksigma2k,Ksigma2pi;
	FILE *fptr;
	printf("<pt> for p,K,pi=%g,%g,%g\n",meanpt_protons/TotPbar,meanpt_kaons/TotKbar,meanpt_pions/TotPibar);
	printf("Qbar=%g Q2bar=%g Q3bar=%g Q4bar=%g\n",
	Qbar/double(ncalls),Q2bar/double(ncalls),Q3bar/double(ncalls),Q4bar/double(ncalls));
	printf("Pbar=%g P2bar=%g P3bar=%g P4bar=%g\n",
	Pbar/double(ncalls),P2bar/double(ncalls),P3bar/double(ncalls),P4bar/double(ncalls));

	pbar= Pbar/double(ncalls);
	kappap2=P2bar/double(ncalls)-pbar*pbar;
	kappap3=P3bar/double(ncalls)-3*kappap2*pbar-pbar*pbar*pbar;
	kappap4=P4bar/double(ncalls)-4*kappap3*pbar-3*kappap2*kappap2-6*kappap2*pbar*pbar-pbar*pbar*pbar*pbar;

	qbar=Qbar/double(ncalls);
	kappaq2=Q2bar/double(ncalls)-qbar*qbar;
	kappaq3=Q3bar/double(ncalls)-3*kappaq2*qbar-qbar*qbar*qbar;
	kappaq4=Q4bar/double(ncalls)-4*kappaq3*qbar-3*kappaq2*kappaq2-6*kappaq2*qbar*qbar-qbar*qbar*qbar*qbar;

	kbar=Kbar/double(ncalls);
	kappak2=K2bar/double(ncalls)-kbar*kbar;
	kappak3=K3bar/double(ncalls)-3*kappak2*kbar-kbar*kbar*kbar;
	kappak4=K4bar/double(ncalls)-4*kappak3*kbar-3*kappak2*kappak2-6*kappak2*kbar*kbar-kbar*kbar*kbar*kbar;

	pibar=Pibar/double(ncalls);
	kappapi2=Pi2bar/double(ncalls)-pibar*pibar;
	kappapi3=Pi3bar/double(ncalls)-3*kappapi2*pibar-pibar*pibar*pibar;
	kappapi4=Pi4bar/double(ncalls)-4*kappapi3*pibar-3*kappapi2*kappapi2-6*kappapi2*pibar*pibar-pibar*pibar*pibar*pibar;

	totq=TotQbar/double(ncalls);
	totp=TotPbar/double(ncalls);
	totk=TotKbar/double(ncalls);
	totpi=TotPibar/double(ncalls);


	sigma2q=kappaq2;
	sigma2p=kappap2;
	sigma2k=kappak2;
	sigma2pi=kappapi2;
	Ssigmaq=kappaq3/kappaq2;
	Ssigmap=kappap3/kappap2;
	Ssigmak=kappak3/kappak2;
	Ssigmapi=kappapi3/kappapi2;
	Ksigma2q=kappaq4/kappaq2;
	Ksigma2p=kappap4/kappap2;
	Ksigma2k=kappak4/kappak2;
	Ksigma2pi=kappapi4/kappapi2;

	printf("Q: <Q>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",qbar,sigma2q,Ssigmaq,Ksigma2q);
	printf("P: <P>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",pbar,sigma2p,Ssigmap,Ksigma2p);
	printf("K: <K>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",kbar,sigma2k,Ssigmak,Ksigma2k);
	printf("Pi: <Pi>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",pibar,sigma2pi,Ssigmapi,Ksigma2pi);

	if((bool)ifstream(file)){
		fptr=fopen(file.c_str(),"a");
	}
	else{
		fptr=fopen(file.c_str(),"w");
		fprintf(fptr,"# Omega - roots -     T    - rhoB  - rhoQ -rhoP - sigma^2(P) - Ssigma(P) - Ksigma^2(P) - rhoQ - sigma^2(Q) - Ssigma(Q) - Ksigma^2(Q) - rhoK - sigma^2(K) - Ssigma(K) - Ksigma^2(K) - rhoPi - sigma^2(Pi) - Ssigma(Pi) - Ksigma^2(Pi) - TotrhoP - TotrhoQ - TotrhoK - TotrhoPi\n");
	}
	fprintf(fptr,"%6.1f %6.1f %6.2f %6.4f %6.4f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
	Omega,roots,T,rhoB,rhoQ,pbar/Omega,sigma2p,Ssigmap,Ksigma2p,qbar/Omega,sigma2q,Ssigmaq,Ksigma2q,
	kbar/Omega,sigma2k,Ssigmak,Ksigma2k,pibar/Omega,sigma2pi,Ssigmapi,Ksigma2pi,totp/Omega,totq/Omega,totk/Omega,totpi/Omega);
	fclose(fptr);
}
