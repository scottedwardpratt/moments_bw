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
	const int AMAX=300;
	vector<double> lnfact;
	int A,Q,Qm,iq,iqm,Apairs,Qtot,Qbar;
	double sigma2mult,Abar,A2bar,normmult,Kmult,omega,omegaA,omegaa,abar,a2bar;
	double z,Ztot;
	printf("Enter z and Qtot: ");
	scanf("%lf %lf",&z,&Qtot);
	int nplus,nminus;
	double alpha=0.5; // efficiency  // now look at charge fluctuations
	double x,sigma2,dPA,K,dPQ,norm,guess;
	double gg=z*z;
	double ProbRate[AMAX+1]={0.0};
	vector<double> Z[AMAX+1];

	printf("Enter efficiency: ");
	scanf("%lf",&alpha);
	x=alpha*(1.0-alpha);
	
	lnfact.resize(AMAX+Qtot+1);
	lnfact[0]=1.0;
	for(int A=0;A<=AMAX;A++){
		ProbRate[A]=0.0;
		Z[A].resize(2*A+1);
		for(Q=-A;Q<=A;Q++){
			iq=Q+A;
			Z[A][iq]=0.0;
		}
	}
	Z[0][0]=1.0;
	lnfact[0]=0.0;
	for(A=1;A<=AMAX;A++){
		lnfact[A]=lnfact[A-1]+log(double(A));
		for(Qm=-A+1;Qm<=A-1;Qm++){
			iqm=Qm+(A-1);
			Q=Qm+1;
			iq=Q+A;
			Z[A][iq]+=(z/double(A))*Z[A-1][iqm];
			Q=Qm-1;
			iq=Q+A;
			Z[A][iq]+=(z/double(A))*Z[A-1][iqm];
		}
	}
	/*
	for(A=0;A<=8;A++){
		for(Q=-A;Q<=A;Q++){
			iqm=Q+A;
			if(pow(-1,iqm)==1){
				guess=pow(z,A)/exp(  lnfact[(A-Q)/2]+lnfact[(A+Q)/2]  );
				printf("Z(A=%d,Q=%d)=%12.4f =? %12.4f\n",A,Q,Z[A][iqm],guess);
			}
		}
	}
	*/
	
	// Now do calc. based on reaction rates for Qtot !=0
	Ztot=0.0;
	for(A=Qtot;A<=AMAX;A++){
		Ztot+=Z[A][A+Qtot];
	}
	
	ProbRate[0]=1.0;
	double ZPR=1.0;
	for(Apairs=2;Apairs<=AMAX;Apairs+=2){
		ProbRate[Apairs]=ProbRate[Apairs-2]*gg/double(0.25*Apairs*(Apairs+2.0*Qtot));
		ZPR+=ProbRate[Apairs];
	}
	double nbarZ=0.0,nbarPR=0.0,normZ=0.0,normPR=0.0;
	double sigma2Z=0.0,sigma2PR=0.0;
	Abar=A2bar=normmult=0.0;
	for(A=Qtot;A<=AMAX;A+=2){
		Apairs=A-Qtot;
		//printf("%3d: %8.6f %8.6f\n",A,Z[A][A+Qtot]/Ztot,ProbRate[Apairs]/ZPR);
		dPA=Z[A][A+Qtot]/Ztot;
		normmult+=dPA;
		Abar+=dPA*A;
		A2bar+=dPA*A*A;
		/* nbarZ+=A*Z[A][A+Qtot]/Ztot;
		nbarPR+=Apairs*ProbRate[Apairs]/PRtot;
		sigma2Z+=A*A*Z[A][iq]/Ztot;
		sigma2PR+=A*A*ProbRate[A]/PRtot;*/
		normPR+=ProbRate[Apairs]/ZPR;
		normZ+=Z[A][A+Qtot]/Ztot;
	}
	Abar=Abar/normmult;
	A2bar=A2bar/normmult;
	A2bar=A2bar-Abar*Abar;
	
	printf("testing, normZ=%g, normPR=%g, omega_M=%g\n",normZ,normPR,A2bar/Abar);

	// Now do moments of charge distributions for Qtot != 0
	Ztot=0.0;
	for(A=Qtot;A<=AMAX;A++){
		Ztot+=Z[A][A+Qtot];
	}
	sigma2=K=norm=0.0;
	normmult=sigma2mult=Abar=Kmult=0.0;
	abar=a2bar=0.0;
	for(A=Qtot;A<=AMAX;A+=2){
		dPA=Z[A][A+Qtot]/Ztot;
		normmult+=dPA;
		sigma2mult+=dPA*A*A;
		Abar+=dPA*A;
		Kmult+=dPA*pow(A-2.0*z,4);
		for(nplus=0;nplus<=A/2;nplus++){
			for(nminus=0;nminus<=A/2;nminus++){
				Q=nplus-nminus;
				dPQ=pow(alpha,nplus)*pow(1.0-alpha,A/2-nplus)*pow(alpha,nminus)*pow(1.0-alpha,A/2-nminus);
				dPQ=dPQ*exp(2*lnfact[A/2]-lnfact[nplus]-lnfact[A/2-nplus]-lnfact[nminus]-lnfact[A/2-nminus]);
				sigma2+=Q*Q*dPA*dPQ;
				K+=Q*Q*Q*Q*dPA*dPQ;
				norm+=dPQ*dPA;
				abar+=(nplus+nminus)*dPQ*dPA;
				a2bar+=(nplus+nminus)*(nplus+nminus)*dPQ*dPA;
			}
		}
	}
	K=K-3.0*sigma2*sigma2;
	K=K/sigma2;
	a2bar=a2bar-abar*abar;
	omegaa=a2bar/abar;
	omegaA=A2bar/Abar;
	printf("omegaA/omegaa=%12.10f\n",omegaA/omegaa);
	printf("abar=%g alpha*Abar=%g\n",abar,alpha*Abar);
	printf("------------------------------------------------------\n");
	printf("CANONICAL ENSEMBLE:\n");
	printf("For net charge: norm=%g, sigma2=%g, Ksigma^2=%g\n",norm,sigma2,K);
	guess=1.0+3.0*alpha*(1.0-alpha)*(omegaA-2.0);
	//guess=1.0-3.0*alpha*(1.0-alpha);
	printf("guess=%g\n",guess);
	sigma2mult=sigma2mult/normmult;
	Abar=Abar/normmult;
	sigma2mult=sigma2mult-Abar*Abar;
	Kmult=Kmult/normmult;
	Kmult=Kmult-3.0*sigma2mult*sigma2mult;
	Kmult=Kmult/sigma2mult;
	printf("For multiplicities: norm=%g, <A>=%g, <delta A^2>=%g, Ksigma^2=%g, omega=%g\n",normmult,Abar,sigma2mult,Kmult,omegaA);
	printf("------------------------------------------------------\n");
	
	
	
	//
	// Now from Poissonian Pair Dist
	// Now do calc. based on reaction rates for Qtot != 0
	sigma2=K=norm=0.0;
	normmult=sigma2mult=Abar=Kmult=0.0;
	Qbar=Qtot*alpha;
	Abar=A2bar=abar=a2bar=0.0;
	for(A=0;A<=AMAX;A+=2){
		dPA=exp(-z)*pow(z,A/2)/exp(lnfact[A/2]);
		normmult+=dPA;
		Abar+=dPA*A;
		A2bar+=dPA*A*A;
		Kmult+=dPA*pow(A-2.0*z,4);
		for(nplus=0;nplus<=A/2+Qtot;nplus++){
			for(nminus=0;nminus<=A/2;nminus++){
				Q=nplus-nminus;
				dPQ=pow(alpha,nplus)*pow(1.0-alpha,A/2+Qtot-nplus)*pow(alpha,nminus)*pow(1.0-alpha,A/2-nminus);
				dPQ=dPQ*exp(lnfact[A/2]+lnfact[A/2+Qtot]-lnfact[nplus]-lnfact[A/2+Qtot-nplus]-lnfact[nminus]-lnfact[A/2-nminus]);
				sigma2+=(Q-Qbar)*(Q-Qbar)*dPA*dPQ;
				abar+=(nplus+nminus)*dPQ*dPA;
				a2bar+=(nplus+nminus)*(nplus+nminus)*dPQ*dPA;
				K+=pow(Q-Qbar,4)*dPA*dPQ;
				norm+=dPQ*dPA;
			}
		}
	}
	A2bar=A2bar-Abar*Abar;
	a2bar=a2bar-abar*abar;
	K=K-3.0*sigma2*sigma2;
	K=K/sigma2;
	printf("POISSONIAN PAIRS:\n");
	printf("For net charge: norm=%g, sigma2=%g, Ksigma^2=%g\n",norm,sigma2,K);
	omegaA=A2bar/Abar;
	omega=a2bar/abar;
	guess=1.0+3.0*alpha*(1.0-alpha)*(omegaA-2.0);
	printf("guess=%g\n",guess);
	sigma2mult=sigma2mult/normmult;
	Abar=Abar/normmult;
	sigma2mult=sigma2mult-Abar*Abar;
	Kmult=Kmult/normmult;
	Kmult=Kmult-3.0*sigma2mult*sigma2mult;
	Kmult=Kmult/sigma2mult;
	printf("For multiplicities: norm=%g, <A>=%g, <delta A^2>=%g, Ksigma^2=%g\n",normmult,Abar,sigma2mult,Kmult);
	return 0;
}
