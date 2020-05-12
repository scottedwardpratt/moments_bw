#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include <list>
#include <unordered_map>
using namespace std;

void BinomialCumulants(int N,double alpha,double &C1,double &C2,double &C3,double &C4){
	double factorial[101];
	int n;
	factorial[0]=1.0;
	for(n=1;n<=100;n++)
		factorial[n]=factorial[n-1]*n;
	C1=N*alpha;
	C2=N*alpha*(1.0-alpha);
	C3=N*alpha*(1.0-alpha)*(1.0-2.0*alpha);
	C4=N*alpha*(1.0-alpha)*(1.0-6.0*alpha*(1.0-alpha));
	double C4test=0.0,C2test=0.0,C1test=0.0,C3test=0.0;
	double Nbar=N*alpha,prob;

	for(n=0;n<=N;n++){
		prob=pow(alpha,n)*pow(1.0-alpha,N-n)*factorial[N]/(factorial[N-n]*factorial[n]);
		C1test+=n*prob;
		C2test+=(n-Nbar)*(n-Nbar)*prob;
		C3test+=pow(n-Nbar,3)*prob;
		C4test+=pow(n-Nbar,4)*prob;
	}
	C4test-=3.0*C2test*C2test;
	printf("C1=%g=?%g, C2=%g=?%g, C3=%g=?%g, C4=%g=?%g\n",C1,C1test,C2,C2test,C3,C3test,C4,C4test);
}

int main(int argc,char *argv[]){
	double C1,C2,C3,C4;
	//BinomialCumulants(10,0.5,C1,C2,C3,C4);
	//exit(1);
	bool equilibrated_pairs=true,bombs=false;  // one should be true, one false
	const int AMAX=200;
	vector<double> lnfact;
	int A,q,Apairs,Qtot,Qbar,nmax;
	double sigma2mult,Abar,A2bar,normcheck,Kmult,omega,omegaA,omegaa,abar,a2bar,qbar,q2bar,q3bar,q4bar;
	double z,Ztot;
	FILE *fptr=fopen("bomb.dat","w");
	//printf("Enter z and Qtot: ");
	//scanf("%lf %lf",&z,&Qtot);
	z=10.0;
	Qtot=0;
	int nplus,nminus;
	double alpha=0.5; // efficiency  // now look at charge fluctuations
	double x,sigma2,dPA,K,dPQ,norm,guess;
	vector<double> Z;
	int nbomb=10,a;
	Z.resize(AMAX+1);

	printf("Enter nbomb(even integer): ");
	scanf("%d",&nbomb);
	lnfact.resize(AMAX+Qtot+1);
	lnfact[0]=0.0;
	for(A=1;A<=AMAX;A++)
		lnfact[A]=lnfact[A-1]+log(double(A));
	
	for(z=0.005;z<50;z++){

		Z[0]=1.0;
		Ztot=Z[0];
		Abar=A2bar=0.0;

		for(A=nbomb;A<=AMAX;A+=nbomb){
			Z[A]=Z[A-nbomb]*z/double(A);
			Abar+=Z[A]*A;
			A2bar+=Z[A]*A*A;
			Ztot+=Z[A];
		}
		Abar=Abar/Ztot;
		A2bar=(A2bar/Ztot)-Abar*Abar;
		omega=A2bar/Abar;
		printf("Abar=%g, A2bar=%g, omega=%g\n",Abar,A2bar,omega);
		fprintf(fptr."%8.5f %10.5f %10.5f\n",z,Abar,omega);
	}
	fclose(fptr);
	
	

	return 0;
}
