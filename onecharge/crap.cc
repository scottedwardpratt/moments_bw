#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include <list>
#include <unordered_map>
using namespace std;


int main(int argc,char *argv[]){
	int N=9,nplus,nminus,Nalpha=100;
	double Q4bar=0.0,Q2bar=0.0,q,netprob=0.0,alpha,aa,prob,dalpha=1.0/double(Nalpha);
	double Factorial[101];
	Factorial[0]=1.0;
	printf("Enter N: ");
	scanf("%d",&N);
	for(int i=1;i<=N;i++){
		Factorial[i]=Factorial[i-1]*double(i);
	}
	//printf("Enter alpha: ");
	//scanf("%lf",&alpha);
	for(aa=0.5*dalpha;aa<1.0;aa+=dalpha){
		//alpha=aa;
		alpha=0.5+1.0*(aa-0.5);
		printf("alpha=%g\n",alpha);
		for(nplus=0;nplus<=N;nplus++){
			for(nminus=0;nminus<=N;nminus++){
				prob=dalpha*pow(alpha,nplus)*pow(1.0-alpha,N-nplus)*pow(alpha,nminus)*pow(1.0-alpha,N-nminus);
				prob=prob*Factorial[N]*Factorial[N]/(Factorial[nplus]*Factorial[nminus]*Factorial[N-nplus]*Factorial[N-nminus]);
				q=nplus-nminus;
				Q2bar+=q*q*prob;
				Q4bar+=q*q*q*q*prob;
				netprob+=prob;
			}
		}
	}
	printf("netprob=%g\n",netprob);
	double C2,C4;
	C2=Q2bar;
	C4=Q4bar-3.0*Q2bar*Q2bar;
	printf("C2=%g, C4=%g\n",Q2bar,Q4bar-3.0*Q2bar*Q2bar);
	alpha=0.5;
	double x=alpha*(1.0-alpha);
	//printf("Guesses=%g, %g\n",2.0*N*x,N*(-12.0*x*x+2.0*x));
	printf("C4/C2=%g =? %g\n",C4/C2,1.0+3.0*x*(-2.0));
	
	return 0;
}
