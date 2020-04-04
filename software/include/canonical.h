#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <armadillo>
#include <vector>
#include <array>
#include <fstream>
#include "resonances.h"
#include "misc.h"
#include <list>
#include <unordered_map>

using namespace std;

class Cpart{
public:
	Cpart(CResInfo *resinfoset){
		resinfo=resinfoset;
	}
	FourVector p;
	CResInfo *resinfo;
};

class CpartitionFunction{
public:
	CpartitionFunction(parameterMap &parmap);
	double T;
	int B0,Q0,S0;  // must be set for MC generation of particles
	vector<vector<vector<vector<double>>>> Z; // Canonical partition function(nhad,b,q,s), multiplied by nhad!
	vector<vector<vector<double>>> Ztot0;
	void ZeroZ(); // Sets Z elements to zero
	void ZeroZtot0(); // Sets Ztot0(Z summed over nhad) to zero
	vector<double> z; //
	bool DECAY; // when generating particles, decay resonances (should be false for blastwave)
	void Calcz();
	CResList *reslist;
	int nres,NhadMAX;
	double Omega,Omega0;  // Omega0 is just for calculating Z before scaling for real volume
	void Getibiqis(int ihad,int b,int q,int s,int &ib,int &iq,int &is);
	int PickNhad(int b0,int q0,int s0);
	void GetNextQ(int ihad0,int b0,int q0,int s0,int &b1,int &q1,int &s1);
	void GetProducts(CResInfo *resinfo0,int &ndaughters,array<CResInfo *,10> &daughter);
	CRandom *randy;
	void CalcZofOmega0(double Tset); //Calculates Z(A,Q) with volume Omega0 (actually stores Z*A!)
	void ScaleZ(double OmegaSet); // Scales Z(A,Q) to get Ztot(Q) with right volume, Set T before setting Omega
	CResInfo* GenHad(int b,int q,int s);
	void GenEvent(int B0,int Q0,int S0,vector<CResInfo *> &resinfo);
	vector<int> bmin;
	vector<int> bmax;
	vector<vector<int>> qmin;
	vector<vector<int>> qmax;
	vector<vector<vector<int>>> smin;
	vector<vector<vector<int>>> smax;
	//vector<vector<vector<list<CResInfo *>>>> resinfobycharge;
	vector<vector<vector<multimap<double,CResInfo *>>>> resinfobycharge;
	void GetQminQmax(int ihad,int b,int &qqmin,int &qqmax);
	void GetSminSmax(int ihad,int b,int q,int &ssmin,int &ssmax);
	bool CheckRelevance(int ihad,int b, int q, int s);
	void CheckQminQmax(int ihad,int b,int q,int s);
	void PrintZ(int ihad);
	void PrintZtot0();
	void WriteZ();
	void ReadZ();
	multimap<double,CResInfo*> massresmap;
};

class Cacceptance{
public:
	CRandom *randy;
	double ACCEPTANCE;
	Cacceptance(double a){
		ACCEPTANCE=a;
		randy=new CRandom(1234);
	}
	bool Acceptance(CResInfo *resinfo){
		if(randy->ran()<ACCEPTANCE)
			return true;
		else
			return false;
	}
	void Acceptance(Cpart *part,bool &acceptQ,bool &acceptP,bool &acceptK,bool &acceptPi,bool &acceptB);
};

class Cmult{
public:
	Cmult();
	int Mplus,Mminus,Mpiplus,Mpiminus,MKplus,MKminus,Mp,Mpbar;
	void Increment(CResInfo *);
	void Zero();
};

class Cmultlist{
public:
	Cmultlist();
	vector<Cmult> multlist;
	void AddMult(vector<CResInfo *> &resinfo);
	void WriteMults(double Omega,double rhoB,double rhoQ);
	void AnalyzeMults(double Omega,double rhoB,double rhoQ,double acceptance);
	double GetProb(int M1,int M2,double a,int n);
	string filename;
};

class CblastWave{
public:
	CRandom *randy;
	double Tf,uperpx,uperpy;
	double sigma_eta,sigma_uperp;
	CResList *reslist;
	CblastWave(parameterMap &parmap,CRandom *randyset,CResList *reslistset);
	void GenerateParts(vector<CResInfo *> &resinfovec,vector<Cpart> &partvec);
	void BoostParts(vector<Cpart> &partvec);
	void GetDecayMomenta(Cpart *mother,int &nbodies,vector<Cpart *> &daughterpartvec);
};

class Cmoments{
public:
	Cmoments(){};
	long long int Qbar,Q2bar,Q3bar,Q4bar,Pbar,P2bar,P3bar,P4bar,Kbar,K2bar,K3bar,K4bar,Pibar,Pi2bar,Pi3bar,Pi4bar,Bbar,B2bar,B3bar,B4bar;
	long long int TotQbar,TotPbar,TotKbar,TotPibar,TotBbar;

	long long int altQbar,altQ2bar,altQ3bar,altQ4bar,altPbar,altP2bar,altP3bar,altP4bar,altKbar,altK2bar,altK3bar,altK4bar,altPibar,altPi2bar,altPi3bar,altPi4bar;
	long long int altTotQbar,altTotPbar,altTotKbar,altTotPibar;

	long long int ncalls;
	double meanpt_pions,meanpt_kaons,meanpt_protons;
	double altmeanpt_pions,altmeanpt_kaons,altmeanpt_protons;

	Cmoments(Cacceptance *acceptanceset);
	void Summarize(string file,string altfile,double Omega,double rhoB,double rhoQ,double roots,double T);
	void IncrementMoments(vector<CResInfo *> &resinfovec);
	void IncrementMoments(vector<Cpart> &partvec);
	void Clear();
	Cacceptance *acceptance;
};
