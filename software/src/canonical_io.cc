#include "canonical.h"

using namespace std;

void CpartitionFunction::WriteZ(){
	FILE *fptr;
	char filename[100];
	int ihad,ihh,ib,iq,is,Amax,Bmin,Bmax,Qmin,Qmax,Smin,Smax;
	for(ihad=0;ihad<=NhadMAX;ihad++){
		ihh=ihad;
		if(ihh>NhadMAX/2)
			ihh=NhadMAX-ihad;
		sprintf(filename,"Zdata/T%gOmega0%g_ihad%d.bin",T,Omega0,ihad);
		fptr=fopen(filename,"wb");
		fprintf(fptr,"%d\n",NhadMAX);
		fprintf(fptr,"%d %d\n",bmin[ihh],bmax[ihh]);
		for(ib=0;ib<=bmax[ihh]-bmin[ihh];ib++){
			fprintf(fptr,"%d %d\n",qmin[ihh][ib],qmax[ihh][ib]);
			for(iq=0;iq<=qmax[ihh][ib]-qmin[ihh][ib];iq++){
				fprintf(fptr,"%d %d\n",smin[ihh][ib][iq],smax[ihh][ib][iq]);
				for(is=0;is<=smax[ihh][ib][iq]-smin[ihh][ib][iq];is++){
					fprintf(fptr,"%g ",Z[ihh][ib][iq][is]);
				}
				fprintf(fptr,"\n");
			}
		}
		fclose(fptr);
	}			
}

void CpartitionFunction::ReadZ(){
	FILE *fptr;
	char filename[100];
	int ihad,ihh,ib,iq,is,Amax,Bmin,Bmax,Qmin,Qmax,Smin,Smax;
	int iib,iiq,iis;
	int bminfile,bmaxfile,qminfile,qmaxfile,sminfile,smaxfile;
	int bfile,qfile,sfile,ihhfile;
	int NhadMAXfile;
	double Zfile;
	for(ihad=0;ihad<=NhadMAX;ihad++){
		printf("Reading for ihad=%d\n",ihad);
		ihhfile=ihad;
		if(ihhfile>NhadMAXfile/2)
			ihhfile=NhadMAXfile-ihad;
		ihh=ihad;
		if(ihh>NhadMAX/2)
			ihh=NhadMAX-ihad;
		sprintf(filename,"Zdata/T%gOmega0%g_ihad%d.bin",T,Omega0,ihad);
		fptr=fopen(filename,"rb");
		fscanf(fptr,"%d",&NhadMAXfile);
		if(NhadMAXfile<NhadMAX){
			printf("ReadZ(): NhadMAX larger than NhadMAXfile\n");
			exit(1);
		}
		fscanf(fptr,"%d %d\n",&bminfile,&bmaxfile);
		for(ib=0;ib<=bmaxfile-bminfile;ib++){
			bfile=bminfile+ib;
			fscanf(fptr,"%d %d\n",&qminfile,&qmaxfile);
			for(iq=0;iq<=qmaxfile-qminfile;iq++){
				qfile=qminfile+iq;
				fscanf(fptr,"%d %d\n",&sminfile,&smaxfile);
				for(is=0;is<=smaxfile-sminfile;is++){
					sfile=sminfile+is;
					fscanf(fptr,"%lf ",&Zfile);
					//printf("Zfile(%d,%d,%d,%d)=%g\n",ihhfile,ib,iq,is,Zfile);
					if(bfile>=bmin[ihh] && bfile<=bmax[ihh]){
						Getibiqis(ihad,bfile,qfile,sfile,iib,iiq,iis);
						if(qfile>=qmin[ihh][iib] && qfile<=qmax[ihh][iib]){
							if(sfile>=smin[ihh][iib][iiq] && sfile<=smax[ihh][iib][iiq]){
								Z[ihad][iib][iiq][iis]=Zfile;
							}
						}
					}
				}
			}
		}
		fclose(fptr);
	}			
}