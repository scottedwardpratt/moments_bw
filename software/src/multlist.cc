#include "canonical.h"

using namespace std;

Cmult::Cmult(){
	Zero();
}
void Cmult::Zero(){
	Mplus=Mminus=Mp=Mpbar=Mpiplus=Mpiminus=MKplus=MKminus=0;
}

void Cmult::Increment(CResInfo *resinfo){
	if(resinfo->charge!=0){
		if(resinfo->charge==1)
			Mplus+=1;
		else
			Mminus+=1;
		if(resinfo->code==211)
			Mpiplus+=1;
		else if(resinfo->code==-211)
			Mpiminus+=1;
		else if(resinfo->code==321)
			MKplus+=1;
		else if(resinfo->code==-321)
			MKminus+=1;
		else if(resinfo->code==2212)
			Mp+=1;
		else if(resinfo->code==-2212)
			Mpbar+=1;
		else{
			printf("What's this particle, pid=%d\n",resinfo->code);
		}
	}
}

Cmultlist::Cmultlist(){
	multlist.clear();
}

void Cmultlist::AddMult(vector<CResInfo *> &resinfo){
	int i,n=resinfo.size();
	Cmult mult;
	mult.Zero();
	for(i=0;i<n;i++){
		mult.Increment(resinfo[i]);
	}
	multlist.push_back(mult);
}

void Cmultlist::WriteMults(double Omega,double rhoB,double rhoQ){
	char fname[120];
	sprintf(fname,"Omega%g_rhoB%g_rhoQ%g.dat",Omega,rhoB,rhoQ);
	string fn="multdata/"+string(fname);
	FILE *fptr=fopen(fn.c_str(),"w");
	long long int i,n=multlist.size();
	for(i=0;i<n;i++){
		fprintf(fptr,"%d %d %d %d %d %d %d %d\n",multlist[i].Mplus,multlist[i].Mminus,
		multlist[i].Mpiplus,multlist[i].Mpiminus,multlist[i].MKplus,multlist[i].MKminus,
		multlist[i].Mp,multlist[i].Mpbar);
	}
	multlist.clear();
	fclose(fptr);
}

double Cmultlist::GetProb(int M1,int M2,double a,int n){
	int ell,ellmax,Mswitch,m;
	double term,sum;
	if(n<0){
		n=abs(n);
		Mswitch=M1;
		M1=M2;
		M2=Mswitch;
	}
	if(n>M1)
		return 0.0;
	term=1.0;
	for(m=1;m<=n;m++){
		term*=double(M1-m+1)/double(m);
	}
	term*=pow(a,n)*pow(1.0-a,M1+M2-n);
	ellmax=M1-n;
	if(M2<ellmax)
		ellmax=M2;
	sum=term;
	for(ell=1;ell<=ellmax;ell++){
		term*=double(M2-ell+1)*double(M1-n-ell+1)/(double(n+ell)*double(ell));
		term*=a*a/((1.0-a)*(1.0-a));
		sum+=term;
	}
	return sum;
}

void Cmultlist::AnalyzeMults(double Omega,double rhoB,double rhoQ,double acceptance){
	char fname[120];
	sprintf(fname,"Omega%g_rhoB%g_rhoQ%g.dat",Omega,rhoB,rhoQ);
	string fn="multdata/"+string(fname);
	FILE *fptr=fopen(fn.c_str(),"r");
	long long int nevents=0;
	int Mplus,Mminus,Mpiplus,Mpiminus,MKplus,MKminus,Mp,Mpbar,n;
	double Qbar=0,Q2bar=0,Q3bar=0,Q4bar=0,Bbar=0,B2bar=0,B3bar=0,B4bar=0;
	double prob;
	double qbar,kappaq2,kappaq3,kappaq4,bbar,kappab2,kappab3,kappab4;
	double sigma2q,sigma2b,Ssigmaq,Ssigmab,Ksigma2q,Ksigma2b;
	
	fscanf(fptr,"%d %d %d %d %d %d %d %d",&Mplus,&Mminus,&Mpiplus,&Mpiminus,&MKplus,&MKminus,&Mp,&Mpbar);
	while(!feof(fptr)){
		nevents+=1;
		fscanf(fptr,"%d %d %d %d %d %d %d %d",&Mplus,&Mminus,&Mpiplus,&Mpiminus,&MKplus,&MKminus,&Mp,&Mpbar);
		for(n=-Mminus;n<=Mplus;n++){
			prob=GetProb(Mplus,Mminus,acceptance,n);
			if(prob>1.0E-17){
				Qbar+=n*prob;
				Q2bar+=n*n*prob;
				Q3bar+=double(n*n*n)*prob;
				Q4bar+=n*n*n*n*prob;
			}
		}
		for(n=-Mpbar;n<=Mp;n++){
			prob= GetProb(Mp,Mpbar,acceptance,n);
			if(prob>1.0E-17){
				Bbar  += n*prob;
				B2bar += n*n*prob;
				B3bar += n*n*n*prob;
				B4bar += n*n*n*n*prob;
			}
		}
	}
	fclose(fptr);
	printf("Qbar=%g, Q2bar=%g, Q3bar=%g, Q4bar=%g\n",Qbar/double(nevents),Q2bar/double(nevents),Q3bar/double(nevents),Q4bar/double(nevents));
	
	qbar=Qbar/double(nevents);
	kappaq2=Q2bar/double(nevents)-qbar*qbar;
	kappaq3=Q3bar/double(nevents)-3*kappaq2*qbar-qbar*qbar*qbar;
	kappaq4=Q4bar/double(nevents)-4*kappaq3*qbar-3*kappaq2*kappaq2-6*kappaq2*qbar*qbar-qbar*qbar*qbar*qbar;
	
	bbar= Bbar/double(nevents);
	kappab2=B2bar/double(nevents)-bbar*bbar;
	kappab3=B3bar/double(nevents)-3*kappab2*bbar-bbar*bbar*bbar;
	kappab4=B4bar/double(nevents)-4*kappab3*bbar-3*kappab2*kappab2-6*kappab2*bbar*bbar-bbar*bbar*bbar*bbar;
	
	sigma2q=kappaq2;
	sigma2b=kappab2;
	Ssigmaq=kappaq3/kappaq2;
	Ssigmab=kappab3/kappab2;
	Ksigma2q=kappaq4/kappaq2;
	Ksigma2b=kappab4/kappab2;

	printf("Q: <Q>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",qbar,sigma2q,Ssigmaq,Ksigma2q);
	printf("B: <B>=%g, sigma^2=%g, Ssigma=%g, Ksigma^2=%g\n",bbar,sigma2b,Ssigmab,Ksigma2b);
	
	fn="momentsdata/"+string(fname);
	if((bool)ifstream(fn)){
		fptr=fopen(fn.c_str(),"a");
	}
	else{
		fptr=fopen(fn.c_str(),"w");
		fprintf(fptr,"# Acceptance - 100*rhoB - sigma^2(B) - Ssigma(B) - Ksigma^2(B) - 100*rhoQ - sigma^2(Q) - Ssigma(Q) - Ksigma^2(Q)\n");
	}
	fprintf(fptr,"%6.4f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	acceptance,bbar/Omega,sigma2b,Ssigmab,Ksigma2b,qbar/Omega,sigma2q,Ssigmaq,Ksigma2q);
	fclose(fptr);
	
	char ffnn[120];
	sprintf(ffnn,"momentsdata/Acceptance%g_rhoB%g_rhoQ%g.dat",acceptance,rhoB,rhoQ);
	fn=string(ffnn);
	if((bool)ifstream(fn)){
		fptr=fopen(fn.c_str(),"a");
	}
	else{
		fptr=fopen(fn.c_str(),"w");
		fprintf(fptr,"# Omega  - 100rhoB - sigma^2(B) - Ssigma(B) - Ksigma^2(B) - 100rhoQ - sigma^2(Q) - Ssigma(Q) - Ksigma^2(Q)\n");
	}
	fprintf(fptr,"%6.1f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	Omega,bbar/Omega,sigma2b,Ssigmab,Ksigma2b,qbar/Omega,sigma2q,Ssigmaq,Ksigma2q);
	fclose(fptr);
}

