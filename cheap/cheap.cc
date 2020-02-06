#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include <list>
#include <unordered_map>
using namespace std;
double GetProb(int M1,int M2,double a,int n){
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


int main(int argc,char *argv[]){
	const int NMAX=200;
	double alpha,Z,P[NMAX+1],PQ[NMAX+1];
	int n;
	printf("Enter alpha: ");
	scanf("%lf",&alpha);
	P[0]=Z=1.0;
	for(n=0;n<NMAX;n++){
		P[n+1]=alpha*P[n]/(double((n+1.0)*(n+1.0)));
		Z+=P[n+1];
	}
	for(n=0;n<=NMAX;n++){
		P[n]=P[n]/Z;
		//printf("P(%d)=%g\n",n,P[n]);
	}
	
	double a,K1=0.0,K2=0.0,K3=0.0,K4=0.0;
	double kappa1,kappa2,kappa3,kappa4;
	double ZQ=0.0;
	int M1,M2,Q;
	printf("Enter acceptance: ");
	scanf("%lf",&a);
	ZQ=0.0;
	for(Q=0;Q<=NMAX;Q++){
		PQ[Q]=0.0;
		for(n=Q;n<=NMAX;n++)
			PQ[Q]+=P[n]*GetProb(n,n,a,Q);
		if(PQ[Q]==PQ[Q]){
			if(Q==0)
				ZQ+=PQ[Q];
			else
				ZQ+=2.0*PQ[Q];
		}
		else
			PQ[Q]=0.0;
	}
	printf("ZQ=%g\n",ZQ);
	K1=K2=K3=K4=0.0;
	for(Q=-NMAX;Q<=NMAX;Q++){
		K1+=PQ[abs(Q)]*Q;
		K2+=PQ[abs(Q)]*Q*Q;
		K3+=PQ[abs(Q)]*Q*Q*Q;
		K4+=PQ[abs(Q)]*Q*Q*Q*Q;
	}
	kappa1=K1;
	kappa2=K2-kappa1*kappa1;
	kappa3=K3-3*kappa2*kappa1-kappa1*kappa1*kappa1;
	kappa4=K4-4*kappa3*kappa1-3*kappa2*kappa2-6*kappa2*kappa1*kappa1-kappa1*kappa1*kappa1*kappa1;
	printf("kappa1=%g, kappa2=%g, kappa3=%g, kappa4=%g, k3/k1=%g, k4/k2=%g\n",
	kappa1,kappa2,kappa3,kappa4,kappa3/kappa1,kappa4/kappa2);
	return 0;
}
