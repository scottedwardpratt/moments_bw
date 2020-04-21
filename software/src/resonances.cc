#ifndef __RESONANCES_CC__
#define __RESONANCES_CC__

#include "resonances.h"

CResList *CResInfo::reslist=NULL;

CResList::~CResList(){
	CResInfo *resinfo;
	CResInfoMap::iterator rpos=resmap.begin();
	while(rpos!=resmap.end()){
		resinfo=rpos->second;
		resmap.erase(rpos);
		delete resinfo;
		rpos=resmap.begin();
	}
}

CResList::CResList(parameterMap* parmap_in){
	parmap=parmap_in;
	RESWIDTH_ALPHA=parameter::getD(*parmap,"B3D_RESWIDTH_ALPHA",0.5);
	RESONANCE_DECAYS=parameter::getB(*parmap,"B3D_RESONANCE_DECAYS",true);
	CResInfo::reslist=this;
	ReadResInfo();
}

CRandom *CResInfo::ranptr=NULL;

CResInfo::CResInfo(){
	minmass=0.0;
	branchlist.clear();
}

CMerge::CMerge(CResInfo *resinfo_in,double branching_in, int L_in){
	resinfo=resinfo_in;
	branching=branching_in;
	L = L_in;
	next=NULL;
}

void CResInfo::DecayGetResInfoPtr(int &nbodies,array<CResInfo *,5> &daughterresinfo){
	double randy,bsum;
	int ibody,ibranch;
	CBranchInfo *bptr;
	bptr=NULL;
	bsum=0.0;
	randy=ranptr->ran();
	ibranch=0;
	do{
		bptr=branchlist[ibranch];
		bsum+=bptr->branching;
		ibranch++;
		if(bsum>1.00000001){
			printf("FATAL: In DecayGetResInfo: bsum too large, = %g\n",bsum);
			exit(1);
		}
	}while(bsum<randy);
	nbodies=bptr->resinfoptr.size();
	for(ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr->resinfoptr[ibody];
	}
}

bool CResInfo::CheckForDaughters(int codecheck){
	//checks to see if any decay daughters match code, for code=0, checking for charged parts
	int ibody,nbodies,ibranch;
	bool exists=false;
	CResInfo *daughter;
	CBranchInfo *bptr;
	CBranchList::iterator bpos;
	if(codecheck!=0){
		if(code==codecheck){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfoptr.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfoptr[ibody];
					if(daughter->code==codecheck){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(codecheck)){
						exists=true;
						return exists;
					}
				}
				bpos++;
			}while(bpos!=branchlist.end());
		}
	}
	else{
		if(charge!=0){
			exists=true;
			return exists;
		}
		if(decay){
			ibranch=0;
			do{
				bptr=branchlist[ibranch];
				ibranch++;
				nbodies=bptr->resinfoptr.size();
				for(ibody=0;ibody<nbodies;ibody++){
					daughter=bptr->resinfoptr[ibody];
					if(daughter->charge!=0){
						exists=true;
						return exists;
					}
					if(daughter->CheckForDaughters(codecheck)){
						exists=true;
						return exists;
					}
				}
				++bpos;
			}while(bpos!=branchlist.end());
		}
	}
	return exists;
}

void CResInfo::DecayGetResInfoPtr_minmass(int &nbodies,array<CResInfo *,5> &daughterresinfo){
	nbodies=bptr_minmass->resinfoptr.size();
	for(int ibody=0;ibody<nbodies;ibody++){
		daughterresinfo[ibody]=bptr_minmass->resinfoptr[ibody];
	}
}

CBranchInfo::CBranchInfo(){
}

bool CResInfo::CheckForNeutral(){
	bool neutral=true;
	if(charge!=0 || strange!=0 || baryon!=0)
		neutral=false;
	return neutral;
}

double CResInfo::GenerateMass(){
	double m,m1,m2;
	int i=0;
	double alpha=reslist->RESWIDTH_ALPHA;
	if(decay){
		printf("check GM a\n");
		m1=branchlist[0]->resinfoptr[0]->mass;
		m2=0.0;
		printf("m1=%g, m2=%g\n",m1,m2);
		for(int n=1;n<(branchlist[0]->resinfoptr.size());n++){
			m2+=branchlist[0]->resinfoptr[n]->mass;
		}
		double kr=sqrt(pow((mass*mass-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0))/(2.0*mass);
		while (i == 0) {
			double r = ranptr->ran();
			m = ((width/2)*tan(PI*(r - .5))) + mass;
			if ((m < (m1+m2))||(m>2.0*mass)) continue;
			double k=sqrt(pow((m*m-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0))/(2.0*m);
			double gamma=width*pow((2.0*k*k)/(k*k+kr*kr),alpha);
			double rho=(2.0/(width*PI))*(0.25*gamma*gamma)/((0.25*gamma*gamma)+(mass-m)*(mass-m));
			double lor = (width/(2*PI))/(pow(width/2,2.0) + pow(mass-m,2.0));
			double weight = rho/(lor*8.0);
			r = ranptr->ran();
			if (r < weight) i = 1;
		}
	}
	else m=mass;
	return m;
}

double CResInfo::GenerateThermalMass(double maxweight, double T){
	double m,kr;
	double alpha=reslist->RESWIDTH_ALPHA;
	if(decay){
		double m1=branchlist[0]->resinfoptr[0]->mass;
		double m2=0.0;
		for(int n=1;n<(branchlist[0]->resinfoptr.size());n++){
			m2+=branchlist[0]->resinfoptr[n]->mass;
		}
		double k2mr = gsl_sf_bessel_Kn(2,(mass/T)); // K2 for resmass
		kr=pow(mass*mass-m1*m1-m2*m2,2) - (4*m1*m1*m2*m2);
		kr = (1/(2*mass))*sqrt(kr); // k at resonant mass
		int i = 0; // for use in while loop
		while(i==0){
			double r1 = ranptr->ran(); // get random numbers
			double r2 = ranptr->ran(); // between [0, 1]
			m = ((width/2)*tan(PI*(r1 - .5))) + mass;// generate random mass value proportional to the lorentz distribution
			if ((m < minmass) ) continue;
			// throw out values out of range
			double k=sqrt(pow((m*m-m1*m1-m2*m2),2.0)-pow((2.0*m1*m2),2.0))/(2.0*m);
			double gamma=width*pow((2.0*k*k)/(k*k+kr*kr),alpha);
			double rho=(2.0/(width*PI))*(0.25*gamma*gamma)/((0.25*gamma*gamma)+(mass-m)*(mass-m));
			double lor = (width/(2*PI))/(pow(width/2,2.0) + pow(mass-m,2.0));
			double k2 = gsl_sf_bessel_Kn(2,(m/T)); // K2 value
			double weight = rho*k2*m*m/(lor*k2mr*mass*mass*maxweight);
			if (r2 < weight) i = 1; // success
		}
	}
	else m=mass;
	return (m); //returns a random mass proportional to n0*L'
}

void CResInfo::Print(){
	printf("+++++++ ID=%d, M=%g, M_min=%g, %s +++++++++\n",code,mass,minmass,name.c_str());
	printf("Gamma=%g, Spin=%g, Decay=%d\n",width,spin,int(decay));
	printf("Q=%d, B=%d, S=%d, G_parity=%d\n",charge,baryon,strange,G_Parity);
}

/*
void CResList::freegascalc_onespecies_offshell(CResInfo *resinfo,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
double width=resinfo->width;
double mass=resinfo->GenerateThermalMass(T);
double minmass=resinfo->minmass;
double wtot=0.0,w;
double m0,dm=0.25*width;
double dP,depsilon,ddens,dsigma2,ddedt;
epsilon=P=dens=sigma2=dedt=0.0;
if(width/T<0.01){
freegascalc_onespecies(mass,T,epsilon,P,dens,sigma2,dedt);
}
else{
m0=mass-4.0*width;
if(m0<minmass)
m0=minmass;
while(m0<mass+4.0*width){
Misc::Pause();
}
}
}
*/

void CResList::freegascalc_onespecies(double m,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt){
	const double prefactor=1.0/(2.0*PI*PI*pow(HBARC,3));
	double k0,k1,z,k0prime,k1prime,m2,m3,m4,t2,t3,I1,I2,Iomega;
	m2=m*m;
	m3=m2*m;
	m4=m2*m2;
	t2=T*T;
	t3=t2*T;
	z=m/T;
	if(z>1000.0){
		P=epsilon=dens=dedt=0.0;
		printf("z is huge=%g, m=%g, t=%g\n",z,m,T);
	}
	else{
		if(z<0.0){
			printf("___z=%g,m=%g,T=%g ___\n",z,m,T);
			exit(1);
		}
		k0=gsl_sf_bessel_K0(z);
		k1=gsl_sf_bessel_K1(z);
		P=prefactor*(m2*t2*k0+2.0*m*t3*k1);
		epsilon=prefactor*(3.0*m2*t2*k0+(m3*T+6.0*m*t3)*k1);
		dens=P/T;
		k0prime=-k1;
		k1prime=-k0-k1/z;
		dedt=prefactor*(6.0*m2*T*k0+(m3+18.0*m*t2)*k1-3.0*m3*k0prime-((m4/T)+6.0*m2*T)*k1prime);
		Iomega=exp(-z)/(30.0*PI*PI*HBARC*HBARC*HBARC);
		I1=pow(m,1.5)*pow(T,3.5)*7.5*sqrt(2.0*PI);
		I2=24.0*pow(T,5);
		sigma2=Iomega*(I1+I2+0.5*sqrt(I1*I2));  // this is an approximation (+/-2%) to messy integral
	}
}

void CResList::freegascalc_onespecies_finitewidth(double resmass, double m1, double m2, double T, double width,
double minmass,double &epsilon,double &P,double &dens,double &sigma2,double &dedt,double &maxweight){

	double kr,k,E,E0,dE,gamma,rho,percent,dp,closest;
	double sum=0.0,esum=0.0,psum=0.0,dsum=0.0,sigsum=0.0,dedtsum=0.0;
	double n0,resn0,lor,weight;
	double alpha=RESWIDTH_ALPHA;
	dE=1.0;
	percent=0.001;
	dp=0.002;
	if (width<1.0){
		dE=0.1;
		percent=0.0005;
		dp=0.001;
	}
	closest=1000.0;
	E0=m1+m2;
	maxweight=0.0;
	kr=sqrt(pow((resmass*resmass-m1*m1-m2*m2),2.0)-4.0*m1*m1*m2*m2)/(2.0*resmass);

	for(E=(m1+m2+0.5*dE);E<2.0*resmass;E+=dE){

		k=sqrt(pow((E*E-m1*m1-m2*m2),2.0)-(4.0*m1*m1*m2*m2))/(2.0*E);
		gamma=width*pow((2.0*k*k)/(k*k+kr*kr),alpha);
		rho=(2.0/(width*PI))*0.25*gamma*gamma/((0.25*gamma*gamma)+(resmass-E)*(resmass-E));
		sum+=rho*dE;

		n0=gsl_sf_bessel_Kn(2,E/T)*E*E*T/(2*PI*PI*pow(HBARC,3.0));
		resn0=gsl_sf_bessel_Kn(2,resmass/T)*resmass*resmass*T/(2*PI*PI*pow(HBARC,3.0));
		lor=(width/(2.0*PI))/(0.25*width*width+(resmass-E)*(resmass-E));
		weight=n0*rho/(resn0*lor);

		if(weight>maxweight) maxweight=weight;
		if (abs(sum-percent)<closest) closest=abs(sum-percent);
		else{
			freegascalc_onespecies(E,T,epsilon,P,dens,sigma2,dedt);
			esum+=epsilon*rho*(E-E0);
			psum+=P*rho*(E-E0);
			dsum+=dens*rho*(E-E0);
			sigsum+=sigma2*rho*(E-E0);
			dedtsum+=dedt*rho*(E-E0);
			closest=1000.0;
			percent+=dp;
			E0=E;
		}
	}

	epsilon=esum/sum;
	P=psum/sum;
	dens=dsum/sum;
	sigma2=sigsum/sum;
	dedt=dedtsum/sum;

}

void CResList::ReadResInfo(){
	CMerge *merge;
	int mothercode,code,decay,NResonances,j,sgn=1;
	double mothermass,bsum,netm,qR2,bmax;
	int ires,jres,ires1,ires2,iresflip,ichannel,nchannels,ibody,nbodies,LDecay;
	int netq,netb,nets;
	string name, filename;
	CResInfo *resinfoptr=NULL, *bresinfoptr=NULL;
	CBranchInfo *bptr=NULL,*firstbptr=NULL;
	FILE *resinfofile;
	FILE * decayinfofile;
	char dummy[200],cname[200];
	filename=parameter::getS(*parmap,"B3D_RESONANCES_INFO_FILE",string("resinfo/resonances_pdg_weak.dat"));
	use_bose_terms=parameter::getB(*parmap,"USE_BOSE_TERMS",false);
	n_bose_terms=parameter::getI(*parmap,"N_BOSE_TERMS",1);
	printf("will read res info from %s\n",filename.c_str());
	resinfofile=fopen(filename.c_str(),"r");
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fgets(dummy,200,resinfofile);
	fscanf(resinfofile,"%d",&NResonances);
	fgets(dummy,200,resinfofile);
	MergeArray=new CMerge **[NResonances];
	SigmaMaxArray=new double *[NResonances];
	for(ires=0;ires<NResonances;ires++){
		MergeArray[ires]=new CMerge *[NResonances];
		SigmaMaxArray[ires]=new double[NResonances];
		for(jres=0;jres<NResonances;jres++){
			MergeArray[ires][jres]=NULL;
			SigmaMaxArray[ires][jres]=0.0;
		}
	}
	for(ires=0;ires<NResonances;ires++){
		resinfoptr=new CResInfo();
		fscanf(resinfofile,"%d %lf %d %d %d %lf %d %d %lf", &resinfoptr->code,&resinfoptr->mass,&resinfoptr->charge,&resinfoptr->baryon, &resinfoptr->strange,&resinfoptr->spin,&resinfoptr->G_Parity,&decay,&resinfoptr->width);
		fgets(cname,100,resinfofile);
		cname[int(strlen(cname))-1]='\0';
		resinfoptr->name=cname;
		resinfoptr->decay=bool(decay);
		//if(RESONANCE_DECAYS) resinfoptr->decay=false;
		resinfoptr->ires=ires;
		resinfoptr->bose_pion=false;
		resinfoptr->branchlist.clear();
		if(!resinfoptr->decay && decay==1){
			resinfoptr->Print();
		}
		resmap.insert(CResInfoPair(resinfoptr->code,resinfoptr));
		if ((resinfoptr->code==211||resinfoptr->code==-211||resinfoptr->code==111) && use_bose_terms){
			if (resinfoptr->code==211) j=0;
			if (resinfoptr->code==-211) j=1;
			if (resinfoptr->code==111) j=2;
			for (int i=0;i<n_bose_terms;i++){ //make extra imaginary resonances for bose corrections
				bresinfoptr=new CResInfo();
				if (resinfoptr->code<0) sgn=-1;
				bresinfoptr->code=resinfoptr->code*1000+(i+2)*sgn; //no real resonance has this code
				bresinfoptr->mass=resinfoptr->mass;
				bresinfoptr->charge=resinfoptr->charge;
				bresinfoptr->baryon=resinfoptr->baryon;
				bresinfoptr->strange=resinfoptr->strange;
				bresinfoptr->spin=resinfoptr->spin;
				bresinfoptr->G_Parity=resinfoptr->G_Parity;
				bresinfoptr->decay=resinfoptr->decay;
				bresinfoptr->width=resinfoptr->width;
				bresinfoptr->name=resinfoptr->name;

				bresinfoptr->ires=NResonances+j*(n_bose_terms-1)+i; //stick the extra resonances at the end of the list
				bresinfoptr->bose_pion=true;
				resinfoptr->branchlist.clear();
				resmap.insert(CResInfoPair(bresinfoptr->code,bresinfoptr));
			}
		}
	}
	fclose(resinfofile);

	filename=parameter::getS(*parmap,"B3D_RESONANCES_DECAYS_FILE",string("resinfo/decays_pdg_weak.dat"));
	printf("will read decay info from %s\n",filename.c_str());
	decayinfofile=fopen(filename.c_str(),"r");
	while(fscanf(decayinfofile,"%d %lf",&mothercode,&mothermass) && !feof(decayinfofile)){
		fgets(dummy,200,decayinfofile);
		fscanf(decayinfofile,"%d %d",&mothercode,&nchannels);
		resinfoptr=GetResInfoPtr(mothercode);
		resinfoptr->minmass=1.0E10;
		bsum=0.0;
		bmax=0.0;
		for(ichannel=0;ichannel<nchannels;ichannel++){
			bptr=new CBranchInfo();
			bptr->resinfoptr.clear();
			resinfoptr->branchlist.push_back(bptr);
			fscanf(decayinfofile,"%d",&nbodies);
			netq=-resinfoptr->charge;
			netb=-resinfoptr->baryon;
			nets=-resinfoptr->strange;
			netm=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				fscanf(decayinfofile,"%d",&code);
				bptr->resinfoptr.push_back(GetResInfoPtr(code));
				netq+=bptr->resinfoptr[ibody]->charge;
				netb+=bptr->resinfoptr[ibody]->baryon;
				nets+=bptr->resinfoptr[ibody]->strange;
				netm+=bptr->resinfoptr[ibody]->mass;
			}
			//total charge and baryon number should be conserved, and shouldn't be larger than single strangeness
			if(netq!=0 || netb!=0 || abs(nets)>1){
				printf("Charge conservation failure while reading decay info,\nnetq=%d, netb=%d, nets=%d\n",netq,netb,nets);
				printf("MOTHER (ichannel=%d, nbodies=%d):\n",ichannel,nbodies);
				resinfoptr->Print();
				printf("DAUGHTERS:\n");
				for(ibody=0;ibody<nbodies;ibody++)
					bptr->resinfoptr[ibody]->Print();
				if(netq!=0 || netb!=0)
					exit(1);
			}
			fscanf(decayinfofile,"%lf %d",&bptr->branching,&LDecay);
			//store two body decays only
			if(nbodies==2){
				ires1=bptr->resinfoptr[0]->ires;
				ires2=bptr->resinfoptr[1]->ires;
				if(ires1>ires2){
					iresflip=ires1; ires1=ires2; ires2=iresflip;
				}
				merge=MergeArray[ires1][ires2];
				if(merge==NULL){
					MergeArray[ires1][ires2]=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				else{
					while(merge->next!=NULL){
						merge=merge->next;
					}
					merge->next=new CMerge(resinfoptr,bptr->branching, LDecay);
				}
				if(resinfoptr->mass>bptr->resinfoptr[0]->mass+bptr->resinfoptr[1]->mass){
					qR2=Misc::triangle(resinfoptr->mass,
					bptr->resinfoptr[0]->mass,bptr->resinfoptr[1]->mass);
					SigmaMaxArray[ires1][ires2]+=
						(bptr->branching*4.0*PI*HBARC*HBARC)*(2.0*resinfoptr->spin+1.0)/
							((2.0*bptr->resinfoptr[0]->spin+1.0)*(2.0*bptr->resinfoptr[1]->spin+1.0)*qR2);
					if(ires2!=ires1)
						SigmaMaxArray[ires2][ires1]=SigmaMaxArray[ires1][ires2];
				}
			}
			bsum+=bptr->branching;
			//if the total mass is smaller than the minimum required mass, replace it
			if(netm<resinfoptr->minmass){
				resinfoptr->minmass=netm;
				resinfoptr->bptr_minmass=bptr;
			}
			// switch places to make sure first branch has largest
			if(bptr->branching>bmax){
				bmax=bptr->branching>bmax;
				if(ichannel>0){
					firstbptr=resinfoptr->branchlist[0];
					resinfoptr->branchlist[0]=bptr;
					resinfoptr->branchlist[ichannel]=firstbptr;
				}
			}
		}
	}
	fclose(decayinfofile);
}

CResInfo* CResList::GetResInfoPtr(int code){
	CResInfoMap::iterator rpos;
	rpos=resmap.find(code);
	if(rpos!=resmap.end())
		return rpos->second;
	else{
		printf("Warning GetResInfoPtr() can't find match for PID=%d\n",code);
		exit(1);
		return NULL;
	}
}

void CResList::CalcEoS(double T0,double Tf,double delT){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	printf("#_____________________\n#  T       s         P        epsilon\n");
	double T,P=0.0,epsilon=0.0,dedT=0.0,dens=0.0,s,m,m1,m2,degen;
	int n;
	double pi,epsiloni,densi,sigma2i,dedti;
	double minmass,width,maxweighti;
	printf("   T         s         P         epsilon     nh      cs^2\n");
	for(T=T0;T<Tf+0.00000001;T+=delT){
		P=epsilon=s=dens=dedT=0.0;
		for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
			resinfoptr=rpos->second;
			if(resinfoptr->code!=22){
				width=resinfoptr->width;
				minmass=resinfoptr->minmass;
				if(resinfoptr->decay){
					m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
					m2=0.0;
					for(n=1;n<(resinfoptr->branchlist[0]->resinfoptr.size());n++){
						m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
					}
				}
				degen=2.0*resinfoptr->spin+1.0;
				m=resinfoptr->mass;

				if((minmass>0.0)&&(width>0.0)) freegascalc_onespecies_finitewidth(m,m1,m2,T,width,minmass,epsiloni,pi,densi,sigma2i,dedti,maxweighti);
				else freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);

				//freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);

				P+=pi*degen;
				dedT+=dedti*degen;
				epsilon+=epsiloni*degen;
				dens+=densi*degen;
				s+=(pi+epsiloni)*degen/T;
			}
		}
		//printf("m=%g, dens=%g, P/T=%g\n",m,dens,P/T);
		printf("%6.2f %10.4e %10.4e %10.4e %10.4e %10.4e\n",T,s,P,epsilon,dens,s/dedT);
	}
}

void CResList::CalcEoS(double T,double &epsilon,double &P,double &nhadrons,vector<double> &density,vector<double> &boseweight,vector<double> &maxweight){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double m,m1,m2,degen,width,minmass;
	double pi,epsiloni,densi,sigma2i,dedti,maxweighti;
	int ires=0,nres,ibose,nbose,n;
	if(boseweight.size()==0){
		boseweight.resize(2);
	}
	P=epsilon=nhadrons=0.0;
	density.clear();
	maxweight.clear();
	nres=resmap.size();
	if(GetResInfoPtr(22)->code==22)
		nres-=1;
	density.resize(nres);
	maxweight.resize(nres);
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;

			m=resinfoptr->mass;
			width=resinfoptr->width;
			minmass=resinfoptr->minmass;

			if(resinfoptr->decay){
				m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
				m2=0.0;
				for(n=1;n<(resinfoptr->branchlist[0]->resinfoptr.size());n++){
					m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
				}
			}

			nbose=1;
			if(abs(resinfoptr->code)==211 || resinfoptr->code==111)
				nbose=boseweight.size()-1;
			if(nbose<1)
				nbose=1;
			density[ires]=0.0;
			maxweight[ires]=0.0;
			boseweight[0]=0.0;
			for(ibose=1;ibose<=nbose;ibose++){
				if((minmass>0.0)&&(width>0.0))
					freegascalc_onespecies_finitewidth(m,m1,m2,T/double(ibose),width,minmass,epsiloni,pi,densi,sigma2i,dedti,maxweighti);
				else
					freegascalc_onespecies(m,T/double(ibose),epsiloni,pi,densi,sigma2i,dedti);
				if(resinfoptr->code==211)
					boseweight[ibose]=boseweight[ibose-1]+densi*degen;
				P+=pi*degen;
				epsilon+=epsiloni*degen;
				density[ires]+=densi*degen;
				maxweight[ires]+=maxweighti;
				nhadrons+=density[ires];
			}
			ires+=1;
		}
	}
	for(ibose=1;ibose<=boseweight.size();ibose++)
		boseweight[ibose]/=boseweight[nbose];
}

void CResList::CalcEoS(double T,double &epsilon,double &P,double &nhadrons,double &cs2,vector<double> &density){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double s,m,m1,m2,degen,dedt;
	double pi,epsiloni,densi,sigma2i,dedti;
	double width,minmass,maxweighti;
	int ires=0,nres,n;
	P=epsilon=nhadrons=dedt=0.0;
	density.clear();
	nres=resmap.size();
	if(GetResInfoPtr(22)->code==22)
		nres-=1;
	density.resize(nres);
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			width=resinfoptr->width;
			minmass=resinfoptr->minmass;
			if(resinfoptr->decay){
				m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
				m2=0.0;
				for(n=1;n<(resinfoptr->branchlist[0]->resinfoptr.size());n++){
					m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
				}
			}
			if((minmass>0.0)&&(width>0.0))
				freegascalc_onespecies_finitewidth(m,m1,m2,T,width,minmass,epsiloni,pi,densi,sigma2i,dedti,maxweighti);
			else
				freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			density[ires]=densi*degen;
			dedt+=dedti*degen;
			nhadrons+=density[ires];
			ires+=1;
		}
	}
	s=(P+epsilon)/T;
	cs2=s/dedt;
	//printf("T=%g, epsilon=%g, P=%g, nhadrons=%g\n",T,epsilon,P,nhadrons);
}

void CResList::CalcEoSandChi(double T){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double P,epsilon,s,m,m1,m2,degen;
	double width,minmass,maxweighti;
	double Q,S,B,chi[3][3],chiBQS[3][3],q[3],BQS[3];
	double pi,epsiloni,densi,sigma2i,dedti;
	int a,b,n;
	for(a=0;a<3;a++){
		for(b=0;b<3;b++)
			chi[a][b]=chiBQS[a][b]=0.0;
	}
	P=epsilon=s=0.0;
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			width=resinfoptr->width;
			minmass=resinfoptr->minmass;
			if(resinfoptr->decay){
				m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
				m2=0.0;
				for(n=1;n<(resinfoptr->branchlist[0]->resinfoptr.size());n++){
					m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
				}
			}
			if((minmass>0.0)&&(width>0.0))
				freegascalc_onespecies_finitewidth(m,m1,m2,T,width,minmass,epsiloni,pi,densi,sigma2i,dedti,maxweighti);
			else freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			s+=(pi+epsiloni)*degen/T;
			Q=resinfoptr->charge;
			B=resinfoptr->baryon;
			S=resinfoptr->strange;
			q[0]=(B+Q);
			q[1]=(2*B+S-Q);
			q[2]=-S;
			BQS[0]=B;
			BQS[1]=Q;
			BQS[2]=-S;
			//sprintf(dummy,"%s",resinfoptr->name.c_str());
			//printf("%40s: u=%3.1f,  d=%3.1f,  s=%3.1f\n",dummy,q[0],q[1],q[2]);
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi[a][b]+=densi*degen*q[a]*q[b];
					chiBQS[a][b]+=densi*degen*BQS[a]*BQS[b];
				}
			}
		}
	}
	char filename[80];
	printf("-------------------\n");
	printf("T=%6.2f  s/T^3=%10.4e\n",T,s*HBARC*HBARC*HBARC/(T*T*T));
	printf("----- CHI_ab for HADRON GAS --------------------\n");
	sprintf(filename,"chi_uds.dat");
	FILE *fptr=fopen(filename,"a");
	fprintf(fptr,"%7.3f ",T);
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			printf("%7.4f ",chi[a][b]/s);
			fprintf(fptr,"%10.3e ",chi[a][b]/s);
		}
		printf("\n");
	}
	fprintf(fptr,"\n");
	fclose(fptr);
	printf("----- CHIBQS_ab for HADRON GAS --------------------\n");
	sprintf(filename,"chi_BQS.dat");
	fptr=fopen(filename,"a");
	fprintf(fptr,"%7.3f ",T);
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			printf("%7.4f ",chiBQS[a][b]/s);
			fprintf(fptr,"%10.3e ",chiBQS[a][b]/s);
		}
		printf("\n");
	}
	fprintf(fptr,"\n");
	fclose(fptr);
	printf("-------------------\n");
}

void CResList::CalcEoSandChi(double T,double &P,double &epsilon,double &s,vector< vector<double> > &chi,vector< vector<double> > &chiBQS){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	bool usepolemass=true;
	double m,m1,m2,degen;
	double width,minmass,maxweighti;
	double Q,S,B,q[3],BQS[3];
	double pi,epsiloni,densi,sigma2i,dedti;
	int a,b,n;
	for(a=0;a<3;a++){
		for(b=0;b<3;b++)
			chi[a][b]=chiBQS[a][b]=0.0;
	}
	P=epsilon=s=0.0;
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			if(usepolemass)
				freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			else{
				width=resinfoptr->width;
				minmass=resinfoptr->minmass;
				if(resinfoptr->decay){
					m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
					m2=0.0;
					for(n=1;n<(resinfoptr->branchlist[0]->resinfoptr.size());n++){
						m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
					}
				}
				if((minmass>0.0)&&(width>0.0))
					freegascalc_onespecies_finitewidth(m,m1,m2,T,width,minmass,epsiloni,pi,densi,sigma2i,dedti,maxweighti);
				else
					freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			}
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			s+=(pi+epsiloni)*degen/T;
			Q=resinfoptr->charge;
			B=resinfoptr->baryon;
			S=resinfoptr->strange;
			q[0]=(B+Q);
			q[1]=(2*B+S-Q);
			q[2]=-S;
			BQS[0]=B;
			BQS[1]=Q;
			BQS[2]=-S;
			//sprintf(dummy,"%s",resinfoptr->name.c_str());
			//printf("%40s: u=%3.1f,  d=%3.1f,  s=%3.1f\n",dummy,q[0],q[1],q[2]);
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					chi[a][b]+=densi*degen*q[a]*q[b];
					chiBQS[a][b]+=densi*degen*BQS[a]*BQS[b];
				}
			}
		}
	}
	char filename[80];
	printf("-------------------\n");
	printf("T=%6.2f  s/T^3=%10.4e\n",T,s*HBARC*HBARC*HBARC/(T*T*T));
	printf("----- CHI_ab for HADRON GAS --------------------\n");
	sprintf(filename,"chi_uds.dat");
	FILE *fptr=fopen(filename,"a");
	fprintf(fptr,"%7.3f ",T);
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			printf("%7.4f ",chi[a][b]/s);
			fprintf(fptr,"%10.3e ",chi[a][b]/s);
		}
		printf("\n");
	}
	fprintf(fptr,"\n");
	fclose(fptr);
	//printf("----- CHIBQS_ab for HADRON GAS --------------------\n");
	sprintf(filename,"chi_BQS.dat");
	fptr=fopen(filename,"a");
	fprintf(fptr,"%7.3f ",T);
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			//printf("%7.4f ",chiBQS[a][b]/s);
			fprintf(fptr,"%10.3e ",chiBQS[a][b]/s);
		}
		//printf("\n");
	}
	fprintf(fptr,"\n");
	fclose(fptr);
	//printf("-------------------\n");
}

void CResList::CalcEosandKubo(double T,double &epsilon,double &P,double &nhadrons,double &sdens,double &kubo){
	CResInfo *resinfoptr;
	CResInfoMap::iterator rpos;
	double m,m1,m2,degen,dedt;
	double width,minmass,maxweighti;
	double pi,epsiloni,densi,sigma2i,dedti;
	double e,p,delp=5.0;
	int ires=0,nres,n;
	vector<double> density;
	P=epsilon=nhadrons=dedt=kubo=0.0;
	density.clear();
	nres=resmap.size();
	//if(GetResInfoPtr(22)->code==22)
	//	nres-=1;
	density.resize(nres);
	for(rpos=resmap.begin();rpos!=resmap.end();rpos++){
		resinfoptr=rpos->second;
		if(resinfoptr->code!=22){
			degen=2.0*resinfoptr->spin+1.0;
			m=resinfoptr->mass;
			width=resinfoptr->width;
			minmass=resinfoptr->minmass;
			if(resinfoptr->decay){
				m1=resinfoptr->branchlist[0]->resinfoptr[0]->mass;
				m2=0.0;
				for(n=1;n<(resinfoptr->branchlist[0]->resinfoptr.size());n++){
					m2+=resinfoptr->branchlist[0]->resinfoptr[n]->mass;
				}
			}
			if((minmass>0.0)&&(width>0.0)) freegascalc_onespecies_finitewidth(m,m1,m2,T,width,minmass,epsiloni,pi,densi,sigma2i,dedti,maxweighti);
			else freegascalc_onespecies(m,T,epsiloni,pi,densi,sigma2i,dedti);
			P+=pi*degen;
			epsilon+=epsiloni*degen;
			density[ires]=densi*degen;
			nhadrons+=density[ires];
			for(p=0.5*delp;p<40*sqrt(m*T);p+=delp){
				e=sqrt(m*m+p*p);
				kubo+=degen*(pow(p,6)/(e*e))*exp(-e/T)*delp;
			}
			ires+=1;
		}
	}
	kubo=kubo/(30.0*PI*PI*T*pow(HBARC,4));
	sdens=(P+epsilon)/T;
}

#endif
