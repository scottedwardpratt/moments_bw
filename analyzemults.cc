#include "canonical.h"
using namespace std;

int main(int argc,char *argv[]){
	Cmultlist multlist;
	double Omega,rhoB,rhoQ;
	//printf("Enter Omega,rhoB,rhoQ:");
	//scanf("%lf %lf %lf",&Omega,&rhoB,&rhoQ);
	for(Omega=50;Omega<500;Omega*=2){
		for(rhoB=0.01;rhoB<0.035;rhoB+=0.01){
			rhoQ=0.5*rhoB;
			for(double a=0.05;a<0.995;a+=0.05){
				printf("---------  Omega=%g, rhoB=%g,  a=%g ----------\n",Omega,rhoB,a);
				multlist.AnalyzeMults(Omega,rhoB,rhoQ,a);
			}
		}
	}
	return 0;
}
