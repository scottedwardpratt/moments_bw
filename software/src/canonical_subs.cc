#include "canonical.h"

using namespace std;

CpartitionFunction::CpartitionFunction(parameterMap &parmap){
	int ihad,ib,iq,is,b,q,s,ihh;
	int qqmin,qqmax,ssmin,ssmax;
	double m;
	CResInfo *resinfo;
	DECAY=true;
	Omega=parameter::getD(parmap,"OMEGA",100.0); // Volume in fm^3
	Omega0=parameter::getD(parmap,"OMEGA0",20.0); // V in fm^3 for calcZ, scaled by Omega later (don't change usually)
	NhadMAX=parameter::getI(parmap,"NhadMAX",35+lrint(0.5*Omega)); // max no. of hadrons in volume, should be even no.
	use_bose_terms=parameter::getB(parmap,"USE_BOSE_TERMS",false);
	n_bose_terms=parameter::getI(parmap,"N_BOSE_TERMS",1);
	if(NhadMAX%2==1)
		NhadMAX+=1;
	T=parameter::getD(parmap,"T",150.0); // Temperature in MeV
	randy=new CRandom(1234);
	reslist=new CResList(&parmap);
	CResInfo::ranptr=randy;
	nres=reslist->resmap.size();
	CResInfoMap::iterator rpos;
	multimap<double,CResInfo*>::iterator mpos;

	bmin.resize(NhadMAX/2+1);
	bmax.resize(NhadMAX/2+1);
	qmin.resize(NhadMAX/2+1);
	qmax.resize(NhadMAX/2+1);
	smin.resize(NhadMAX/2+1);
	smax.resize(NhadMAX/2+1);
	for(ihad=0;ihad<=NhadMAX/2;ihad++){
		bmin[ihad]=-ihad;
		bmax[ihad]=ihad;
		qmin[ihad].resize(bmax[ihad]-bmin[ihad]+1);
		qmax[ihad].resize(bmax[ihad]-bmin[ihad]+1);
		smin[ihad].resize(bmax[ihad]-bmin[ihad]+1);
		smax[ihad].resize(bmax[ihad]-bmin[ihad]+1);
		if(ihad==0)
			qmin[0][0]=qmax[0][0]=0;

		for(ib=0;ib<=bmax[ihad]-bmin[ihad];ib++){
			b=-ihad+ib;
			GetQminQmax(ihad,b,qqmin,qqmax);
			qmin[ihad][ib]=qqmin; qmax[ihad][ib]=qqmax;
			smin[ihad][ib].resize(qqmax-qqmin+1);
			smax[ihad][ib].resize(qqmax-qqmin+1);
			if(ihad==0){
				smin[0][0][0]=smax[0][0][0]=0;
			}
			for(iq=0;iq<=(qmax[ihad][ib]-qmin[ihad][ib]);iq++){
				q=iq+qmin[ihad][ib];
				GetSminSmax(ihad,b,q,ssmin,ssmax);
				smin[ihad][ib][iq]=ssmin; smax[ihad][ib][iq]=ssmax;
			}
		}
	}

	z.resize(nres);
	Z.resize(NhadMAX+1);
	for(ihad=0;ihad<=NhadMAX;ihad++){
		ihh=ihad;
		if(ihh>NhadMAX/2)
			ihh=NhadMAX-ihad;
		Z[ihad].resize(bmax[ihh]-bmin[ihh]+1);
		for(ib=0;ib<=bmax[ihh]-bmin[ihh];ib++){
			b=bmin[ihh]+ib;
			Z[ihad][ib].resize(qmax[ihh][ib]-qmin[ihh][ib]+1);
			for(iq=0;iq<=qmax[ihh][ib]-qmin[ihh][ib];iq++){
				q=qmin[ihh][ib]+iq;
				Z[ihad][ib][iq].resize(smax[ihh][ib][iq]-smin[ihh][ib][iq]+1);
				for(is=0;is<=(smax[ihh][ib][iq]-smin[ihh][ib][iq]);is++){
					Z[ihad][ib][iq][is]=0.0;
				}
			}
		}
	}
	ihh=NhadMAX/2;
	Ztot0.resize(bmax[ihh]-bmin[ihh]+1);
	for(ib=0;ib<=bmax[ihh]-bmin[ihh];ib++){
		Ztot0[ib].resize(qmax[ihh][ib]-qmin[ihh][ib]+1);
		for(iq=0;iq<=qmax[ihh][ib]-qmin[ihh][ib];iq++){
			Ztot0[ib][iq].resize(smax[ihh][ib][iq]-smin[ihh][ib][iq]+1);
			for(is=0;is<=(smax[ihh][ib][iq]-smin[ihh][ib][iq]);is++){
				Ztot0[ib][iq][is]=0.0;
			}
		}
	}
	resinfobycharge.resize(3);
	for(ib=0;ib<3;ib++){
		b=bmin[1]+ib;
		resinfobycharge[ib].resize(qmax[1][ib]-qmin[1][ib]+1);
		for(iq=0;iq<=qmax[1][ib]-qmin[1][ib];iq++){
			resinfobycharge[ib][iq].resize(smax[1][ib][iq]-smin[1][ib][iq]+1);
			for(is=0;is<=smax[1][ib][iq]-smin[1][ib][iq];is++)
				resinfobycharge[ib][iq][is].clear();
		}
	}
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			b=resinfo->baryon;
			q=resinfo->charge;
			s=resinfo->strange;
			CheckQminQmax(1,b,q,s);
			m=resinfo->mass;
			Getibiqis(1,b,q,s,ib,iq,is);
			resinfobycharge[ib][iq][is].insert(pair<double,CResInfo *>(resinfo->mass,resinfo));
		}
	}
	//Calcz(); // Also calculates Z[1][ib][iq][is]
	printf("Initialized CpartitionFunction NhadMAX=%d\n",NhadMAX);
}

void CpartitionFunction::ZeroZ(){
	int ihad,ib,iq,is,b,q,ihh;
	for(ihad=2;ihad<=NhadMAX;ihad++){
		ihh=ihad;
		if(ihh>NhadMAX/2)
			ihh=NhadMAX-ihad;
		for(ib=0;ib<=bmax[ihh]-bmin[ihh];ib++){
			b=bmin[ihh]+ib;
			for(iq=0;iq<=qmax[ihh][ib]-qmin[ihh][ib];iq++){
				q=qmin[ihh][ib]+iq;
				for(is=0;is<=(smax[ihh][ib][iq]-smin[ihh][ib][iq]);is++){
					Z[ihad][ib][iq][is]=0.0;
				}
			}
		}
	}
}

void CpartitionFunction::ZeroZtot0(){
	int ib,iq,is,ihh;
	ihh=NhadMAX/2;
	for(ib=0;ib<=bmax[ihh]-bmin[ihh];ib++){
		for(iq=0;iq<=qmax[ihh][ib]-qmin[ihh][ib];iq++){
			for(is=0;is<=(smax[ihh][ib][iq]-smin[ihh][ib][iq]);is++){
				Ztot0[ib][iq][is]=0.0;
			}
		}
	}
}

void CpartitionFunction::CheckQminQmax(int ihad,int b,int q,int s){
	int ihh=ihad,ib,iq;
	bool answer=false;
	if(ihh>NhadMAX/2)
		ihh=NhadMAX-ihad;
	if(b>=bmin[ihh] || b<=bmax[ihh]){
		ib=b-bmin[ihh];
		if(q>=qmin[ihh][ib] || q<=qmax[ihh][ib]){
			iq=q-qmin[ihh][ib];
			if(s>=smin[ihh][ib][iq] || s<=smax[ihh][ib][iq]){
				answer=true;
			}
		}
	}
	if(answer==false)
		printf("\n b or q or s out of range\n");
}

void CpartitionFunction::GetQminQmax(int ihad,int b,int &qqmin,int &qqmax){
	int bprime,ibprime,qmin1,qmax1,db,dq,ihcharge;
	const int temp_nhcharge=11;
	int nhcharge=temp_nhcharge;
	if (use_bose_terms) nhcharge+=(n_bose_terms-1)*2;
	int hcharge[nhcharge][2];
	int temp_hcharge[temp_nhcharge][2]={{1,-1},{1,0},{1,1},{1,2},
	{-1,1},{-1,0},{-1,-1},{-1,-2},
	{0,1},{0,0},{0,-1}};
	for (int i=0;i<temp_nhcharge;i++) {
		for (int j=0;j<2;j++) {
			hcharge[i][j]=temp_hcharge[i][j];
		}
	}
	if (use_bose_terms) {
		for (int i=0;i<n_bose_terms-1;i++) {
			hcharge[temp_nhcharge+2*i][0]=0;
			hcharge[temp_nhcharge+2*i][1]=i+2;
			hcharge[temp_nhcharge+2*i+1][0]=0;
			hcharge[temp_nhcharge+2*i+1][1]=-(i+2);
		}
	}
	qmin1=qmax1=qqmin=qqmax=0;
	if(ihad>0){
		for(ihcharge=0;ihcharge<nhcharge;ihcharge++){
			db=hcharge[ihcharge][0];
			dq=hcharge[ihcharge][1];
			bprime=b-db;
			if(bprime>=bmin[ihad-1] && bprime<=bmax[ihad-1]){
				ibprime=bprime-bmin[ihad-1];
				qmax1=qmax[ihad-1][ibprime]+dq;
				qqmax=max(qqmax,qmax1);
				qmin1=qmin[ihad-1][ibprime]+dq;
				qqmin=min(qqmin,qmin1);
			}
		}
	}
}

void CpartitionFunction::GetSminSmax(int ihad,int b,int q,int &ssmin,int &ssmax){
	int bprime,ibprime,qprime,iqprime,smin1,smax1;
	int db,dq,ds,ihcharge;
	const int temp_nhcharge=27;
	int nhcharge=temp_nhcharge;
	if (use_bose_terms) nhcharge+=(n_bose_terms-1)*2;
	int hcharge[nhcharge][3];
	const int temp_hcharge[temp_nhcharge][3]={{1,-1,0},{1,0,0},{1,1,0},{1,2,0},
	{1,-1,-1},{1,0,-1},{1,1,-1},
	{1,-1,-2},{1,0,-2},
	{1,-1,-3},
	{-1,1,0},{-1,0,0},{-1,-1,0},{-1,-2,0},
	{-1,1,1},{-1,0,1},{-1,-1,1},
	{-1,1,2},{-1,0,2},
	{-1,1,3},
	{0,1,1},{0,-1,-1},{0,0,1},{0,0,-1},
	{0,1,0},{0,0,0},{0-1,0}};
	for (int i=0;i<temp_nhcharge;i++) {
		for (int j=0;j<3;j++) {
			hcharge[i][j]=temp_hcharge[i][j];
		}
	}
	if (use_bose_terms) {
		for (int i=0;i<n_bose_terms-1;i++) {
			hcharge[temp_nhcharge+2*i][0]=0;
			hcharge[temp_nhcharge+2*i][1]=i+2;
			hcharge[temp_nhcharge+2*i][2]=0;
			hcharge[temp_nhcharge+2*i+1][0]=0;
			hcharge[temp_nhcharge+2*i+1][1]=-(i+2);
			hcharge[temp_nhcharge+2*i+1][2]=0;
		}
	}

	smin1=smax1=0;
	ssmin=0;
	ssmax=0;

	if(ihad>0){
		for(ihcharge=0;ihcharge<nhcharge;ihcharge++){
			db=hcharge[ihcharge][0];
			dq=hcharge[ihcharge][1];
			ds=hcharge[ihcharge][2];

			bprime=b-db;
			if(bprime>=bmin[ihad-1] && bprime<=bmax[ihad-1]){
				ibprime=bprime-bmin[ihad-1];
				qprime=q-dq;
				if(qprime>=qmin[ihad-1][ibprime] && qprime<=qmax[ihad-1][ibprime]){
					iqprime=qprime-qmin[ihad-1][ibprime];
					smax1=smax[ihad-1][ibprime][iqprime]+ds;
					ssmax=max(ssmax,smax1);
					smin1=smin[ihad-1][ibprime][iqprime]+ds;
					ssmin=min(ssmin,smin1);
				}
			}
		}
	}
}

void CpartitionFunction::Getibiqis(int ihad,int b,int q,int s,int &ib,int &iq,int &is){
	int ihh=ihad;
	if(2*ihh>NhadMAX)
		ihh=NhadMAX-ihh;
	ib=b-bmin[ihh];
	iq=q-qmin[ihh][ib];
	is=s-smin[ihh][ib][iq];
}

void CpartitionFunction::Calcz(){
	int b,q,s,ib,iq,is,ires=0,ihh=1;
	double m,p,e,dens,sigma2,dedt;
	double Ti;
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	double rhoH=0.0;
	Z[0][0][0][0]=1.0;
	for(ib=0;ib<=bmax[ihh]-bmin[ihh];ib++){
		for(iq=0;iq<=qmax[ihh][ib]-qmin[ihh][ib];iq++){
			for(is=0;is<=smax[ihh][ib][iq]-smin[ihh][ib][iq];is++){
				Z[1][ib][iq][is]=0.0;
			}
		}
	}
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		ires=resinfo->ires;
		if(ires>=nres || ires<0){
			printf("ires out of range\n");
			exit(1);
		}
		if(resinfo->code!=22){
			b=resinfo->baryon;
			q=resinfo->charge;
			if (resinfo->bose_pion==true) {
				Ti=T/(abs(resinfo->code)%100);
			}
			else Ti=T;
			s=resinfo->strange;
			CheckQminQmax(1,b,q,s);
			m=resinfo->mass;
			reslist->freegascalc_onespecies(m,Ti,e,p,dens,sigma2,dedt);
			z[ires]=Omega0*dens*(2*resinfo->spin+1.0);
			rhoH+=z[ires]/Omega0;
			Getibiqis(1,b,q,s,ib,iq,is);
			Z[1][ib][iq][is]+=z[ires];
			//resinfobycharge[ib][iq][is].insert(pair<double,CResInfo *>(resinfo->mass,resinfo));
		}
		else{
			z[ires]=0.0;
		}
	}
}

void CpartitionFunction::ScaleZ(double OmegaSet){
	Omega=OmegaSet;
	int b,q,s,ihh,ihad,ib,iq,is,ibprime,iqprime,isprime,ihhZtot0=NhadMAX/2;
	double dZ0;
	ZeroZtot0();

	for(ib=0;ib<=(bmax[ihhZtot0]-bmin[ihhZtot0]);ib++){
		b=bmin[ihhZtot0]+ib;
		for(iq=0;iq<=(qmax[ihhZtot0][ib]-qmin[ihhZtot0][ib]);iq++){
			q=qmin[ihhZtot0][ib]+iq;
			for(is=0;is<=(smax[ihhZtot0][ib][iq]-smin[ihhZtot0][ib][iq]);is++)
				Ztot0[ib][iq][is]=0.0;
		}
	}
	double factor=1.0;
	for(ihad=0;ihad<=NhadMAX;ihad++){
		ihh=ihad;
		if(2*ihad>NhadMAX)
			ihh=NhadMAX-ihh;
		for(ib=0;ib<=(bmax[ihh]-bmin[ihh]);ib++){
			b=bmin[ihh]+ib;
			for(iq=0;iq<=(qmax[ihh][ib]-qmin[ihh][ib]);iq++){
				q=qmin[ihh][ib]+iq;
				for(is=0;is<=(smax[ihh][ib][iq]-smin[ihh][ib][iq]);is++){
					s=smin[ihh][ib][iq]+is;
					if(CheckRelevance(ihhZtot0,b,q,s)){
						Getibiqis(ihhZtot0,b,q,s,ibprime,iqprime,isprime);
						dZ0=Z[ihad][ib][iq][is]*factor;
						Ztot0[ibprime][iqprime][isprime]+=dZ0;
					}
				}
			}
		}
		factor*=(Omega/Omega0)/(ihad+1.0);
	}
	//printf("Z[NhadMAX]/Ztot0=%g,  should be small, otherwise increase NhadMAX\n",dZ0/Ztot0[ibprime][iqprime][isprime]);
}

void CpartitionFunction::PrintPofA(){
	double factor;
	int b=0,q=0,s=0,ihad,ib,iq,is,ibprime,iqprime,isprime,ihhZtot0=NhadMAX/2;
	factor=1.0;
	Getibiqis(ihhZtot0,b,q,s,ibprime,iqprime,isprime);
	for(ihad=0;ihad<=NhadMAX;ihad++){
		Getibiqis(ihad,b,q,s,ib,iq,is);
		printf("%3d %g\n",ihad,Z[ihad][ib][iq][is]*factor/Ztot0[ibprime][iqprime][isprime]);
		factor*=(Omega/Omega0)/(ihad+1.0);
	}
}

void CpartitionFunction::CalcZofOmega0(double TSet){
	T=TSet;
	int ihad,ib,iq,is,b,q,s,ihh;
	int ihad1,ib1,iq1,is1,b1,q1,s1;
	int ihadprime,ibprime,iqprime,isprime,bprime,qprime,sprime;
	ZeroZ();
	Calcz();
	for(ihad=1;ihad<NhadMAX;ihad++){
		ihh=ihad;
		if(2*ihad>NhadMAX)
			ihh=NhadMAX-ihh;
		for(ib=0;ib<=(bmax[ihh]-bmin[ihh]);ib++){
			b=bmin[ihh]+ib;
			for(iq=0;iq<=(qmax[ihh][ib]-qmin[ihh][ib]);iq++){
				q=qmin[ihh][ib]+iq;
				for(is=0;is<=(smax[ihh][ib][iq]-smin[ihh][ib][iq]);is++){
					s=smin[ihh][ib][iq]+is;
					ihad1=1;
					for(ib1=0;ib1<=bmax[ihad1]-bmin[ihad1];ib1++){
						b1=bmin[ihad1]+ib1;
						for(iq1=0;iq1<=qmax[ihad1][ib1]-qmin[ihad1][ib1];iq1++){
							q1=qmin[ihad1][ib1]+iq1;
							for(is1=0;is1<=smax[ihad1][ib1][iq1]-smin[ihad1][ib1][iq1];is1++){
								s1=smin[ihad1][ib1][iq1]+is1;
								ihadprime=ihad+1;
								bprime=b+b1;
								qprime=q+q1;
								sprime=s+s1;
								if(CheckRelevance(ihadprime,bprime,qprime,sprime)){
									Getibiqis(ihadprime,bprime,qprime,sprime,ibprime,iqprime,isprime);
									Z[ihadprime][ibprime][iqprime][isprime]
										+=Z[1][ib1][iq1][is1]*Z[ihad][ib][iq][is];
								}
							}
						}
					}
				}
			}
		}
	}
}

int CpartitionFunction::PickNhad(int B0,int Q0,int S0){
	int nhad=-1,ib,iq,is;
	double rancheck,ransum,nfact=1.0;
	Getibiqis(NhadMAX/2,B0,Q0,S0,ib,iq,is);
	rancheck=Ztot0[ib][iq][is]; //Ztot0[ib][iq][is]*randy->ran();
	//printf("rancheck=%lf\n",rancheck);
	ransum=0.0;
	do{
		nhad+=1;
		if(CheckRelevance(nhad,B0,Q0,S0)){
			Getibiqis(nhad,B0,Q0,S0,ib,iq,is);
			ransum+=pow(Omega/Omega0,nhad)*Z[nhad][ib][iq][is]/nfact;
			//printf("Z[%d] contribution: %lf ransum: %Lf\n",nhad,pow(Omega/Omega0,nhad)*Z[nhad][ib][iq][is]/nfact,ransum);
		}
		nfact*=(nhad+1.0);
	}while(rancheck>ransum);
	//printf("done!\n");
	//exit(1);
	return nhad;
}

CResInfo* CpartitionFunction::GenHad(int b,int q,int s){
	int ib,iq,is,ires;
	bool success=false;
	//list<CResInfo *>::iterator iter;
	multimap<double,CResInfo *>::iterator iter;
	Getibiqis(1,b,q,s,ib,iq,is);
	CResInfo *resinfo;
	iter=resinfobycharge[ib][iq][is].begin();
	double zsum=0.0,rancheck=randy->ran();

	while(iter!=resinfobycharge[ib][iq][is].end() && !success){
		resinfo=iter->second;
		ires=resinfo->ires;
		zsum+=z[ires];
		if(rancheck*Z[1][ib][iq][is]<zsum)
			success=true;
		++iter;
	}
	if(!success){
		printf("oops, we ended without success!!\n");
		exit(1);
	}
	return resinfo;
}

void CpartitionFunction::GetNextQ(int ihad0,int b0,int q0,int s0,int &b1,int &q1,int &s1){
	int ib0,iq0,is0,ib1,iq1,is1,ihad1;
	int bprime,qprime,sprime,ibprime,iqprime,isprime,ihadprime;
	double rancheck=randy->ran();
	double zsum=0.0;
	bool success=false;
	Getibiqis(ihad0,b0,q0,s0,ib0,iq0,is0);
	ihad1=1;
	for(ib1=0;ib1<=bmax[1]-bmin[1];ib1++){
		b1=bmin[1]+ib1;
		for(iq1=0;iq1<=qmax[1][ib1]-qmin[1][ib1];iq1++){
			q1=qmin[1][ib1]+iq1;
			for(is1=0;is1<=smax[1][ib1][iq1]-smin[1][ib1][iq1];is1++){
				s1=smin[1][ib1][iq1]+is1;

				ihadprime=ihad0-1;
				bprime=b0-b1;
				qprime=q0-q1;
				sprime=s0-s1;
				if(CheckRelevance(ihadprime,bprime,qprime,sprime)){
					Getibiqis(ihadprime,bprime,qprime,sprime,ibprime,iqprime,isprime);
					zsum+=Z[ihadprime][ibprime][iqprime][isprime]*Z[1][ib1][iq1][is1];

					if(zsum>rancheck*Z[ihad0][ib0][iq0][is0]){
						success=true;
						goto GNQ_FOUNDIT;
					}
				}
			}
		}
	}
	GNQ_FOUNDIT:
	if(!success){
		printf("FAILURE in GetNextQ()\n");
		exit(1);
	}

}

void CpartitionFunction::GenEvent(int B0,int Q0,int S0,vector<CResInfo *> &resinfo){
	CResInfo *resinfo0;
	int b0,q0,s0;
	resinfo.clear();
	int ihad0,mult=0;
	int ndaughters,idaughter;
	array<CResInfo *,10> daughter;
	b0=B0; q0=Q0; s0=S0;
	int b1,q1,s1;
	ihad0=PickNhad(B0,Q0,S0);
	while(ihad0>0){
		GetNextQ(ihad0,b0,q0,s0,b1,q1,s1);
		ihad0-=1;
		b0=b0-b1;
		q0=q0-q1;
		s0=s0-s1;
		resinfo0=GenHad(b1,q1,s1);
		if(resinfo0->decay && DECAY){
			GetProducts(resinfo0,ndaughters,daughter);
			for(idaughter=0;idaughter<ndaughters;idaughter++){
				resinfo.push_back(daughter[idaughter]);
				mult+=1;
			}
		}
		else{
			resinfo.push_back(resinfo0);
			mult+=1;
		}
	}
}

void CpartitionFunction::PrintZ(int ihad){
	int ib,iq,is,b,q,s;
	printf("---------------- Z(ihad=%d) ----------------\n",ihad);
	for(ib=0;ib<=bmax[ihad]-bmin[ihad];ib++){
		b=bmin[ihad]+ib;
		for(iq=0;iq<=qmax[ihad][ib]-qmin[ihad][ib];iq++){
			q=qmin[ihad][ib]+iq;
			printf("------- b=%d, q=%d -------\n",b,q);
			for(is=0;is<=smax[ihad][ib][iq]-smin[ihad][ib][iq];is++){
				s=smin[ihad][ib][iq]+is;
				printf("%d: %g ",s,Z[ihad][ib][iq][is]);
			}
			printf("\n");
		}
	}
}

void CpartitionFunction::PrintZtot0(){
	int ib,iq,is,b,q,s,ihad=NhadMAX/2;
	printf("---------------- Ztot0 ----------------\n");
	for(ib=0;ib<=bmax[ihad]-bmin[ihad];ib++){
		b=bmin[ihad]+ib;
		for(iq=0;iq<=qmax[ihad][ib]-qmin[ihad][ib];iq++){
			q=qmin[ihad][ib]+iq;
			printf("------- b=%d, q=%d -------\n",b,q);
			for(is=0;is<=smax[ihad][ib][iq]-smin[ihad][ib][iq];is++){
				s=smin[ihad][ib][iq]+is;
				printf("%d: %g ",s,Z[ihad][ib][iq][is]);
			}
			printf("\n");
		}
	}
	printf("finished Printing\n");
}

bool CpartitionFunction::CheckRelevance(int ihad,int b, int q, int s){
	int ihh,ib,iq;
	bool answer=false;
	ihh=ihad;
	if(2*ihh>NhadMAX)
		ihh=NhadMAX-ihad;
	if(b>=bmin[ihh] && b<=bmax[ihh]){
		ib=b-bmin[ihh];
		if(q>=qmin[ihh][ib] && q<=qmax[ihh][ib]){
			iq=q-qmin[ihh][ib];
			if(s>=smin[ihh][ib][iq] && s<=smax[ihh][ib][iq]){
				answer=true;
			}
		}
	}
	return answer;
}
