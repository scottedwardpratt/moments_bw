#include "canonical.h"

using namespace std;

Cmoments::Cmoments(Cacceptance *acceptanceset){
	acceptance=acceptanceset;
	nevents=0;
	Clear();
}

void Cmoments::Clear(){
	nevents=0;
	Bbar=B2bar=B3bar=B4bar=Qbar=Q2bar=Q3bar=Q4bar=Pbar=P2bar=P3bar=P4bar=Kbar=K2bar=K3bar=K4bar=Pibar=Pi2bar=Pi3bar=Pi4bar=0.0;
	TotBbar=TotQbar=TotPbar=TotKbar=TotPibar=0.0;
	meanpt_pions=0.0;
	meanpt_kaons=0.0;
	meanpt_protons=0.0;

	altQbar=altQ2bar=altQ3bar=altQ4bar=altPbar=altP2bar=altP3bar=altP4bar=altKbar=altK2bar=altK3bar=altK4bar=altPibar=altPi2bar=altPi3bar=altPi4bar=0.0;
	altTotQbar=altTotPbar=altTotKbar=altTotPibar=0.0;
	altmeanpt_pions=0.0;
	altmeanpt_kaons=0.0;
	altmeanpt_protons=0.0;
}

/*void Cmoments::AltIncrementMoments(vector<Cpart> &partvec){
	bool acceptQ,acceptP,acceptK,acceptPi,acceptB;
	double pt;
	CResInfo *resinfo;
	nevents+=1;
	int NetQ=0,NetP=0,NetK=0,NetPi=0,NetB=0;
	int TotQ=0,TotP=0,TotK=0,TotPi=0,TotB=0;

	int altNetQ=0,altNetP=0,altNetK=0,altNetPi=0;
	int altTotQ=0,altTotP=0,altTotK=0,altTotPi=0;
	for(int ipart=0;ipart<partvec.size();ipart++){
		resinfo=partvec[ipart].resinfo;
		pt=sqrt(partvec[ipart].p[1]*partvec[ipart].p[1]+partvec[ipart].p[2]*partvec[ipart].p[2]);
		if(resinfo->bose_pion==true) {
			int n=abs(resinfo->code%100);
			for (int i=1;i<=n;i++) {
				acceptance->Acceptance(&partvec[ipart],acceptQ,acceptP,acceptK,acceptPi,acceptB);
				if (acceptPi) {
					NetQ+=resinfo->charge/double(n);
					TotQ+=1;
					NetPi+=resinfo->charge/double(n);
					TotPi+=1;
				}
			}
		}
		else {
			acceptance->Acceptance(&partvec[ipart],acceptQ,acceptP,acceptK,acceptPi,acceptB);
			//acceptQ=acceptP=acceptK=acceptPi=acceptB=true;
			if(abs(resinfo->charge)==1 && acceptQ){
				if(false) { //(resinfo->bose_pion==true){
					NetQ+=resinfo->charge;
					TotQ+=abs(resinfo->code%100);
				}
				else {
					NetQ+=resinfo->charge;
					TotQ+=1;
				}
			} else {
				altNetQ+=resinfo->charge;
				altTotQ+=1;
			}
			if(abs(resinfo->code)==2212 && acceptP){
				NetP+=resinfo->charge;
				TotP+=1;
				meanpt_protons+=pt;
			} else {
				altNetP+=resinfo->charge;
				altTotP+=1;
				altmeanpt_protons+=pt;
			}
			if(abs(resinfo->code)==321 && acceptK){
				NetK+=resinfo->charge;
				TotK+=1;
				meanpt_kaons+=pt;
			} else {
				altNetK+=resinfo->charge;
				altTotK+=1;
				altmeanpt_kaons+=pt;
			}
			if((abs(resinfo->code)==211 || resinfo->bose_pion) && acceptPi){
				if(false) { //if(resinfo->bose_pion==true){
					NetPi+=resinfo->charge;
					TotPi+=abs(resinfo->code%100);
				}
				else {
					NetPi+=resinfo->charge;
					TotPi+=1;
				}
				meanpt_pions+=pt;
			} else {
				altNetPi+=resinfo->charge;
				altTotPi+=1;
				altmeanpt_pions+=pt;
			}
			if(abs(resinfo->baryon)==1 && acceptB){
				NetB+=resinfo->baryon;
				TotB+=1;
			}
		}
	}
	Bbar+=NetB; B2bar+=NetB*NetB; B3bar+=pow(NetB,3); B4bar+=pow(NetB,4);
	Qbar+=NetQ; Q2bar+=NetQ*NetQ; Q3bar+=pow(NetQ,3); Q4bar+=pow(NetQ,4);
	Pbar+=NetP; P2bar+=NetP*NetP; P3bar+=pow(NetP,3); P4bar+=pow(NetP,4);
	Kbar+=NetK; K2bar+=NetK*NetK; K3bar+=pow(NetK,3); K4bar+=pow(NetK,4);
	Pibar+=NetPi; Pi2bar+=NetPi*NetPi; Pi3bar+=pow(NetPi,3); Pi4bar+=pow(NetPi,4);

	TotBbar+=TotB; TotQbar+=TotQ; TotPbar+=TotP; TotKbar+=TotK; TotPibar+=TotPi;

	altQbar+=altNetQ; altQ2bar+=altNetQ*altNetQ; altQ3bar+=pow(altNetQ,3); altQ4bar+=pow(altNetQ,4);
	altPbar+=altNetP; altP2bar+=altNetP*altNetP; altP3bar+=pow(altNetP,3); altP4bar+=pow(altNetP,4);
	altKbar+=altNetK; altK2bar+=altNetK*altNetK; altK3bar+=pow(altNetK,3); altK4bar+=pow(altNetK,4);
	altPibar+=altNetPi; altPi2bar+=altNetPi*altNetPi; altPi3bar+=pow(altNetPi,3); altPi4bar+=pow(altNetPi,4);
	//printf("%lld\n",altQbar);

	altTotQbar+=altTotQ; altTotPbar+=altTotP; altTotKbar+=altTotK; altTotPibar+=altTotPi;
}*/

void Cmoments::ResetEvent(){
	NetQ=NetP=NetK=NetPi=NetB=0;
	TotQ=TotP=TotK=TotPi=TotB=0;
}

void Cmoments::IncrementMoments(vector<Cpart> &partvec){
	bool acceptQ,acceptP,acceptK,acceptPi,acceptB;
	double pt;
	CResInfo *resinfo;

	for(int ipart=0;ipart<partvec.size();ipart++){
		resinfo=partvec[ipart].resinfo;
		pt=sqrt(partvec[ipart].p[1]*partvec[ipart].p[1]+partvec[ipart].p[2]*partvec[ipart].p[2]);
		
		acceptance->Acceptance(&partvec[ipart],acceptQ,acceptP,acceptK,acceptPi,acceptB);
		//acceptQ=acceptP=acceptK=acceptPi=acceptB=true;
		if(abs(resinfo->charge)==1 && acceptQ){
			NetQ+=resinfo->charge;
			TotQ+=1;
		}
		
		if(abs(resinfo->code)==2212 && acceptP){
			NetP+=resinfo->charge;
			TotP+=1;
			meanpt_protons+=pt;
		}
		else if(abs(resinfo->code)==321 && acceptK){
			NetK+=resinfo->charge;
			TotK+=1;
			meanpt_kaons+=pt;
		}
		else if((abs(resinfo->code)==211) && acceptPi){	
			NetPi+=resinfo->charge;
			TotPi+=1;
			meanpt_pions+=pt;
		}
		
		if(abs(resinfo->baryon)==1 && acceptB){
			NetB+=resinfo->baryon;
			TotB+=1;
		}
	}
	nevents+=1;
	Bbar+=NetB; B2bar+=NetB*NetB; B3bar+=pow(NetB,3); B4bar+=pow(NetB,4);
	Qbar+=NetQ; Q2bar+=NetQ*NetQ; Q3bar+=pow(NetQ,3); Q4bar+=pow(NetQ,4);
	Pbar+=NetP; P2bar+=NetP*NetP; P3bar+=pow(NetP,3); P4bar+=pow(NetP,4);
	Kbar+=NetK; K2bar+=NetK*NetK; K3bar+=pow(NetK,3); K4bar+=pow(NetK,4);
	Pibar+=NetPi; Pi2bar+=NetPi*NetPi; Pi3bar+=pow(NetPi,3); Pi4bar+=pow(NetPi,4);
	TotBbar+=TotB; TotQbar+=TotQ; TotPbar+=TotP; TotKbar+=TotK; TotPibar+=TotPi;
}

/*void Cmoments::AltSummarize(string file,string altfile,double Omega,double rhoB,double rhoQ,double roots,double T){
	double qbar,kappaq2,kappaq3,kappaq4;
	double bbar,kappab2,kappab3,kappab4;
	double pbar,kappap2,kappap3,kappap4;
	double kbar,kappak2,kappak3,kappak4;
	double pibar,kappapi2,kappapi3,kappapi4;
	double totb,totq,totp,totk,totpi;
	double sigma2b,sigma2q,sigma2p,sigma2k,sigma2pi,Ssigmab,Ssigmaq,Ssigmap,Ssigmak,Ssigmapi,Ksigma2b,Ksigma2q,Ksigma2p,Ksigma2k,Ksigma2pi;

	double altqbar,altkappaq2,altkappaq3,altkappaq4;
	double altpbar,altkappap2,altkappap3,altkappap4;
	double altkbar,altkappak2,altkappak3,altkappak4;
	double altpibar,altkappapi2,altkappapi3,altkappapi4;
	double alttotq,alttotp,alttotk,alttotpi;
	double altsigma2q,altsigma2p,altsigma2k,altsigma2pi,altSsigmaq,altSsigmap,altSsigmak,altSsigmapi,altKsigma2q,altKsigma2p,altKsigma2k,altKsigma2pi;

	FILE *fptr;
	printf("<pt> for p,K,pi=%g,%g,%g\n",meanpt_protons/TotPbar,meanpt_kaons/TotKbar,meanpt_pions/TotPibar);
	printf("Qbar=%g Q2bar=%g Q3bar=%g Q4bar=%g\n",
	Qbar/double(nevents),Q2bar/double(nevents),Q3bar/double(nevents),Q4bar/double(nevents));
	printf("Pbar=%g P2bar=%g P3bar=%g P4bar=%g\n",
	Pbar/double(nevents),P2bar/double(nevents),P3bar/double(nevents),P4bar/double(nevents));

	pbar= Pbar/double(nevents);
	kappap2=P2bar/double(nevents)-pbar*pbar;
	kappap3=P3bar/double(nevents)-3*kappap2*pbar-pbar*pbar*pbar;
	kappap4=P4bar/double(nevents)-4*kappap3*pbar-3*kappap2*kappap2-6*kappap2*pbar*pbar-pbar*pbar*pbar*pbar;

	qbar=Qbar/double(nevents);
	kappaq2=Q2bar/double(nevents)-qbar*qbar;
	kappaq3=Q3bar/double(nevents)-3*kappaq2*qbar-qbar*qbar*qbar;
	kappaq4=Q4bar/double(nevents)-4*kappaq3*qbar-3*kappaq2*kappaq2-6*kappaq2*qbar*qbar-qbar*qbar*qbar*qbar;

	bbar=Bbar/double(nevents);
	kappab2=B2bar/double(nevents)-bbar*bbar;
	kappab3=B3bar/double(nevents)-3*kappab2*bbar-bbar*bbar*bbar;
	kappab4=B4bar/double(nevents)-4*kappab3*bbar-3*kappab2*kappab2-6*kappab2*bbar*bbar-bbar*bbar*bbar*bbar;

	kbar=Kbar/double(nevents);
	kappak2=K2bar/double(nevents)-kbar*kbar;
	kappak3=K3bar/double(nevents)-3*kappak2*kbar-kbar*kbar*kbar;
	kappak4=K4bar/double(nevents)-4*kappak3*kbar-3*kappak2*kappak2-6*kappak2*kbar*kbar-kbar*kbar*kbar*kbar;

	pibar=Pibar/double(nevents);
	kappapi2=Pi2bar/double(nevents)-pibar*pibar;
	kappapi3=Pi3bar/double(nevents)-3*kappapi2*pibar-pibar*pibar*pibar;
	kappapi4=Pi4bar/double(nevents)-4*kappapi3*pibar-3*kappapi2*kappapi2-6*kappapi2*pibar*pibar-pibar*pibar*pibar*pibar;

	totb=TotBbar/double(nevents);
	totq=TotQbar/double(nevents);
	totp=TotPbar/double(nevents);
	totk=TotKbar/double(nevents);
	totpi=TotPibar/double(nevents);

	sigma2b=kappab2;
	sigma2q=kappaq2;
	sigma2p=kappap2;
	sigma2k=kappak2;
	sigma2pi=kappapi2;
	Ssigmab=kappab3/kappab2;
	Ssigmaq=kappaq3/kappaq2;
	Ssigmap=kappap3/kappap2;
	Ssigmak=kappak3/kappak2;
	Ssigmapi=kappapi3/kappapi2;
	Ksigma2b=kappab4/kappab2;
	Ksigma2q=kappaq4/kappaq2;
	Ksigma2p=kappap4/kappap2;
	Ksigma2k=kappak4/kappak2;
	Ksigma2pi=kappapi4/kappapi2;

	printf("Q: <Q>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",qbar,sigma2q,Ssigmaq,Ksigma2q);
	printf("P: <P>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",pbar,sigma2p,Ssigmap,Ksigma2p);
	printf("K: <K>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",kbar,sigma2k,Ssigmak,Ksigma2k);
	printf("Pi: <Pi>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",pibar,sigma2pi,Ssigmapi,Ksigma2pi);

	altpbar=altPbar/double(nevents);
	altkappap2=altP2bar/double(nevents)-altpbar*altpbar;
	altkappap3=altP3bar/double(nevents)-3*altkappap2*altpbar-altpbar*altpbar*altpbar;
	altkappap4=altP4bar/double(nevents)-4*altkappap3*altpbar-3*altkappap2*altkappap2-6*altkappap2*altpbar*altpbar-altpbar*altpbar*altpbar*altpbar;

	altqbar=altQbar/double(nevents);
	altkappaq2=altQ2bar/double(nevents)-altqbar*altqbar;
	altkappaq3=altQ3bar/double(nevents)-3*altkappaq2*altqbar-altqbar*altqbar*altqbar;
	altkappaq4=altQ4bar/double(nevents)-4*altkappaq3*altqbar-3*altkappaq2*altkappaq2-6*altkappaq2*altqbar*altqbar-altqbar*altqbar*altqbar*altqbar;

	altkbar=altKbar/double(nevents);
	altkappak2=altK2bar/double(nevents)-altkbar*altkbar;
	altkappak3=altK3bar/double(nevents)-3*altkappak2*altkbar-altkbar*altkbar*altkbar;
	altkappak4=altK4bar/double(nevents)-4*altkappak3*altkbar-3*altkappak2*altkappak2-6*altkappak2*altkbar*altkbar-altkbar*altkbar*altkbar*altkbar;

	altpibar=altPibar/double(nevents);
	altkappapi2=altPi2bar/double(nevents)-altpibar*altpibar;
	altkappapi3=altPi3bar/double(nevents)-3*altkappapi2*altpibar-altpibar*altpibar*altpibar;
	altkappapi4=altPi4bar/double(nevents)-4*altkappapi3*altpibar-3*altkappapi2*altkappapi2-6*altkappapi2*altpibar*altpibar-altpibar*altpibar*altpibar*altpibar;

	alttotq=altTotQbar/double(nevents);
	alttotp=altTotPbar/double(nevents);
	alttotk=altTotKbar/double(nevents);
	alttotpi=altTotPibar/double(nevents);


	altsigma2q=altkappaq2;
	altsigma2p=altkappap2;
	altsigma2k=altkappak2;
	altsigma2pi=altkappapi2;
	altSsigmaq=altkappaq3/altkappaq2;
	altSsigmap=altkappap3/altkappap2;
	altSsigmak=altkappak3/altkappak2;
	altSsigmapi=altkappapi3/altkappapi2;
	altKsigma2q=altkappaq4/altkappaq2;
	altKsigma2p=altkappap4/altkappap2;
	altKsigma2k=altkappak4/altkappak2;
	altKsigma2pi=altkappapi4/altkappapi2;

	if((bool)ifstream(file)){
		fptr=fopen(file.c_str(),"a");
	}
	else{
		fptr=fopen(file.c_str(),"w");
		fprintf(fptr,"# Omega - roots -     T    - rhoB  - rhoQ -rhoP - sigma^2(P) - Ssigma(P) - Ksigma^2(P) - rhoQ - sigma^2(Q) - Ssigma(Q) - Ksigma^2(Q) - rhoB - sigma^2(B) - Ssigma(B) - Ksigma^2(B) - rhoK - sigma^2(K) - Ssigma(K) - Ksigma^2(K) - rhoPi - sigma^2(Pi) - Ssigma(Pi) - Ksigma^2(Pi) - TotrhoP - TotrhoQ - TotrhoK - TotrhoPi\n");
	}
	fprintf(fptr,"%6.1f %6.1f %6.2f %6.4f %6.4f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
	Omega,roots,T,rhoB,rhoQ,pbar/Omega,sigma2p,Ssigmap,Ksigma2p,qbar/Omega,sigma2q,Ssigmaq,Ksigma2q,bbar/Omega,sigma2b,Ssigmab,Ksigma2b,
	kbar/Omega,sigma2k,Ssigmak,Ksigma2k,pibar/Omega,sigma2pi,Ssigmapi,Ksigma2pi,totp/Omega,totq/Omega,totk/Omega,totpi/Omega);
	fclose(fptr);

	if((bool)ifstream(altfile)){
		fptr=fopen(altfile.c_str(),"a");
	}
	else{
		fptr=fopen(altfile.c_str(),"w");
		fprintf(fptr,"# Omega - roots -     T    - rhoB  - rhoQ -rhoP - sigma^2(P) - Ssigma(P) - Ksigma^2(P) - rhoQ - sigma^2(Q) - Ssigma(Q) - Ksigma^2(Q) - rhoK - sigma^2(K) - Ssigma(K) - Ksigma^2(K) - rhoPi - sigma^2(Pi) - Ssigma(Pi) - Ksigma^2(Pi) - TotrhoP - TotrhoQ - TotrhoK - TotrhoPi\n");
	}
	fprintf(fptr,"%6.1f %6.1f %6.2f %6.4f %6.4f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
	Omega,roots,T,rhoB,rhoQ,altpbar/Omega,altsigma2p,altSsigmap,altKsigma2p,altqbar/Omega,altsigma2q,altSsigmaq,altKsigma2q,
	altkbar/Omega,altsigma2k,altSsigmak,altKsigma2k,altpibar/Omega,altsigma2pi,altSsigmapi,altKsigma2pi,alttotp/Omega,alttotq/Omega,alttotk/Omega,alttotpi/Omega);
	fclose(fptr);
}*/

void Cmoments::CalcCumulants(){

	/* printf("<pt> for p,K,pi=%g,%g,%g\n",meanpt_protons/TotPbar,meanpt_kaons/TotKbar,meanpt_pions/TotPibar);
	printf("Qbar=%g Q2bar=%g Q3bar=%g Q4bar=%g\n",
	Qbar/double(nevents),Q2bar/double(nevents),Q3bar/double(nevents),Q4bar/double(nevents));
	printf("Pbar=%g P2bar=%g P3bar=%g P4bar=%g\n",
	Pbar/double(nevents),P2bar/double(nevents),P3bar/double(nevents),P4bar/double(nevents));*/

	pbar= Pbar/double(nevents);
	kappap2=P2bar/double(nevents)-pbar*pbar;
	kappap3=P3bar/double(nevents)-3*kappap2*pbar-pbar*pbar*pbar;
	kappap4=P4bar/double(nevents)-4*kappap3*pbar-3*kappap2*kappap2-6*kappap2*pbar*pbar-pbar*pbar*pbar*pbar;

	qbar=Qbar/double(nevents);
	kappaq2=Q2bar/double(nevents)-qbar*qbar;
	kappaq3=Q3bar/double(nevents)-3*kappaq2*qbar-qbar*qbar*qbar;
	kappaq4=Q4bar/double(nevents)-4*kappaq3*qbar-3*kappaq2*kappaq2-6*kappaq2*qbar*qbar-qbar*qbar*qbar*qbar;

	bbar=Bbar/double(nevents);
	kappab2=B2bar/double(nevents)-bbar*bbar;
	kappab3=B3bar/double(nevents)-3*kappab2*bbar-bbar*bbar*bbar;
	kappab4=B4bar/double(nevents)-4*kappab3*bbar-3*kappab2*kappab2-6*kappab2*bbar*bbar-bbar*bbar*bbar*bbar;

	kbar=Kbar/double(nevents);
	kappak2=K2bar/double(nevents)-kbar*kbar;
	kappak3=K3bar/double(nevents)-3*kappak2*kbar-kbar*kbar*kbar;
	kappak4=K4bar/double(nevents)-4*kappak3*kbar-3*kappak2*kappak2-6*kappak2*kbar*kbar-kbar*kbar*kbar*kbar;

	pibar=Pibar/double(nevents);
	kappapi2=Pi2bar/double(nevents)-pibar*pibar;
	kappapi3=Pi3bar/double(nevents)-3*kappapi2*pibar-pibar*pibar*pibar;
	kappapi4=Pi4bar/double(nevents)-4*kappapi3*pibar-3*kappapi2*kappapi2-6*kappapi2*pibar*pibar-pibar*pibar*pibar*pibar;

	totb=TotBbar/double(nevents);
	totq=TotQbar/double(nevents);
	totp=TotPbar/double(nevents);
	totk=TotKbar/double(nevents);
	totpi=TotPibar/double(nevents);
}

void Cmoments::Sum(vector<vector<Cmoments *>> &moments){
	int ieta,iuperp;
	int Neta=moments.size();
	int Nuperp=moments[0].size();
	pbar=kappap2=kappap3=kappap4=0.0;
	qbar=kappaq2=kappaq3=kappaq4=0.0;
	bbar=kappab2=kappab3=kappab4=0.0;
	kbar=kappak2=kappak3=kappak4=0.0;
	pibar=kappapi2=kappapi3=kappapi4=0.0;
	totb=totq=totp=totk=totpi=0.0;
	/*
	double wpsum,wqsum,wbsum,wksum,wpisum;
	wpsum=wqsum=wbsum=wksum=wpisum=0.0;
	for(ieta=0;ieta<Neta;ieta++){
		for(iuperp=0;iuperp<Nuperp;iuperp++){
			wpsum+=moments[ieta][iuperp]->kappap2;
			wqsum+=moments[ieta][iuperp]->kappaq2;
			wbsum+=moments[ieta][iuperp]->kappab2;
			wksum+=moments[ieta][iuperp]->kappak2;
			wpisum+=moments[ieta][iuperp]->kappapi2;
		}
	}*/
	
	for(ieta=0;ieta<Neta;ieta++){
		for(iuperp=0;iuperp<Nuperp;iuperp++){
			pbar+=moments[ieta][iuperp]->pbar;
			qbar+=moments[ieta][iuperp]->qbar;
			bbar+=moments[ieta][iuperp]->pbar;
			kbar+=moments[ieta][iuperp]->kbar;
			pibar+=moments[ieta][iuperp]->pibar;

			kappap2+=moments[ieta][iuperp]->kappap2;
			kappaq2+=moments[ieta][iuperp]->kappaq2;
			kappab2+=moments[ieta][iuperp]->kappab2;
			kappak2+=moments[ieta][iuperp]->kappak2;
			kappapi2+=moments[ieta][iuperp]->kappapi2;
			
			kappap3+=moments[ieta][iuperp]->kappap3;
			kappaq3+=moments[ieta][iuperp]->kappaq3;
			kappab3+=moments[ieta][iuperp]->kappab3;
			kappak3+=moments[ieta][iuperp]->kappak3;
			kappapi3+=moments[ieta][iuperp]->kappapi3;
			
			kappap4+=moments[ieta][iuperp]->kappap4;
			kappaq4+=moments[ieta][iuperp]->kappaq4;
			kappab4+=moments[ieta][iuperp]->kappab4;
			kappak4+=moments[ieta][iuperp]->kappak4;
			kappapi4+=moments[ieta][iuperp]->kappapi4;
			
			totp+=moments[ieta][iuperp]->totp;
			totq+=moments[ieta][iuperp]->totq;
			totb+=moments[ieta][iuperp]->totb;
			totk+=moments[ieta][iuperp]->totk;
			totpi+=moments[ieta][iuperp]->totpi;
		}
	}
	
}


void Cmoments::Summarize(string filename,double Omega,double rhoB,double rhoQ,double roots,double T){
	FILE *fptr;
	double Ssigmab,Ssigmaq,Ssigmap,Ssigmak,Ssigmapi,Ksigma2b,Ksigma2q,Ksigma2p,Ksigma2k,Ksigma2pi;
	Ssigmab=kappab3/kappab2;
	Ssigmaq=kappaq3/kappaq2;
	Ssigmap=kappap3/kappap2;
	Ssigmak=kappak3/kappak2;
	Ssigmapi=kappapi3/kappapi2;
	Ksigma2b=kappab4/kappab2;
	Ksigma2q=kappaq4/kappaq2;
	Ksigma2p=kappap4/kappap2;
	Ksigma2k=kappak4/kappak2;
	Ksigma2pi=kappapi4/kappapi2;

	printf("Q: <Q>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",qbar,kappaq2,Ssigmaq,Ksigma2q);
	printf("P: <P>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",pbar,kappap2,Ssigmap,Ksigma2p);
	printf("K: <K>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",kbar,kappak2,Ssigmak,Ksigma2k);
	printf("Pi: <Pi>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",pibar,kappapi2,Ssigmapi,Ksigma2pi);

	if((bool)ifstream(filename)){
		fptr=fopen(filename.c_str(),"a");
	}
	else{
		fptr=fopen(filename.c_str(),"w");
		fprintf(fptr,"# Omega - roots -     T    - rhoB  - rhoQ -rhoP - sigma^2(P) - Ssigma(P) - Ksigma^2(P) - rhoQ - sigma^2(Q) - Ssigma(Q) - Ksigma^2(Q) - rhoB - sigma^2(B) - Ssigma(B) - Ksigma^2(B) - rhoK - sigma^2(K) - Ssigma(K) - Ksigma^2(K) - rhoPi - sigma^2(Pi) - Ssigma(Pi) - Ksigma^2(Pi) - TotrhoP - TotrhoQ - TotrhoK - TotrhoPi\n");
	}
	fprintf(fptr,"%6.1f %6.1f %6.2f %6.4f %6.4f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
	Omega,roots,T,rhoB,rhoQ,
	pbar/Omega,kappap2,Ssigmap,Ksigma2p,
	qbar/Omega,kappaq2,Ssigmaq,Ksigma2q,
	bbar/Omega,kappab2,Ssigmab,Ksigma2b,
	kbar/Omega,kappak2,Ssigmak,Ksigma2k,
	pibar/Omega,kappapi2,Ssigmapi,Ksigma2pi,
	totp/Omega,totq/Omega,totk/Omega,totpi/Omega);

	fclose(fptr);
}
