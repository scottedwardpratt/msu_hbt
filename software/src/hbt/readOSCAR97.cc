#include "msu_hbt/hbt.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>

using namespace std;
using namespace NMSUPratt;

void Chbt_master::ReadOSCAR_1997(){
	Chbt_cell *cell;
	Chbt_part *tmp_particle=new Chbt_part();
	double taucompare=parmap.getD("OSCAR_TAUCOMPRE",25.0);
	string filenamefilenames=parmap.getS("OSCAR_FILENAME_FILENAMES","oscarfilenames_URQMD.txt");
	double BMIN=parmap.getD("OSCAR_BMIN",4.0);
	double BMAX=parmap.getD("OSCAR_BMAX",4.6);

	// opening the files from first -> last:
	FILE *fptr_in;

	double t, x, y, z, mass, p0, px, py, pz, bim, dumbo;
	int pid,pdg;
	long long naccept=0;
	int nrParticlesInEvent,tracknumber=0;
	int nr_event;
	int ifile=0,nfilesmax;
	nfilesmax=parmap.getI("OSCAR_NFILESMAX",99999);

	list<string> oscar_filenames;
	list<string>::iterator fiter;
	string filename;
	char dummy[800];
	FILE  *fptr_filenames=fopen(filenamefilenames.c_str(),"r");
	do{
		fscanf(fptr_filenames,"%s",dummy);
		filename=dummy;
		oscar_filenames.push_back(filename);
		ifile+=1;
	}while(!feof(fptr_filenames) && ifile<nfilesmax);
	fclose(fptr_filenames);
	
	for(fiter=oscar_filenames.begin();fiter!=oscar_filenames.end();++fiter){
		filename=*fiter;
		CLog::Info("Reading "+filename+"\n");
		fptr_in=fopen(filename.c_str(),"r");
		//getting rid of some lines
		for(int i = 0; i < 3; i++){
			fgets(dummy,800,fptr_in);
			//printf("%s\n",dummy);
		}
		do{
			fscanf(fptr_in,"%d %d %lf %lf",&nr_event,&nrParticlesInEvent,&bim,&dumbo);
			//printf("nr_event=%d, nrParticlesInEvent=%d, bim=%g, dumbo=%g\n",nr_event,nrParticlesInEvent,bim,dumbo);
			fgets(dummy,800,fptr_in);
			//printf("will try to read %d particles from this event\n",nrParticlesInEvent);
			if(!feof(fptr_in)){
				for(int i = 1; i <= nrParticlesInEvent; i++){//reading particles in event loop
					fscanf(fptr_in,"%d %d  %lf %lf %lf %lf  %lf  %lf %lf %lf %lf",&tracknumber,&pid,&px,&py,&pz,&p0,&mass,&x,&y,&z,&t);
					fgets(dummy,800,fptr_in);
					if(i!=tracknumber){
						CLog::Info("Warning for file "+filename+", tracknumber suspicious\n");
						CLog::Info("nr_event="+to_string(nr_event)+",tracknumber="+to_string(tracknumber)+", i="+to_string(i)+"\n");
						exit(1);
					}
					if(bim>=BMIN && bim<=BMAX){
						mass*=1000.0; p0*=1000.0; px*=1000.0; py*=1000.0; pz*=1000.0;
						p0=sqrt(mass*mass+px*px+py*py+pz*pz);
						pdg=pid;
						if(ANTIPARTSYMM && pid==-PIDA)
							pdg=PIDA;
						if(ANTIPARTSYMM && pid==-PIDB)
							pdg=PIDB;
						bool accept; double eff;						
						if(pdg == PIDA || pdg==PIDB){
								if(OVERRIDE_GAUSS){
								//printf("m=%g, p=(%g,%g,%g,%g), r=(%g,%g,%g,%g)\n",mass,p0,px,py,pz,t,x,y,z);
	
								double pdotr,psquared,gamma;
								x=GAUSS_Rx*randy->ran_gauss();
								y=GAUSS_Ry*randy->ran_gauss();
								z=GAUSS_Rz*randy->ran_gauss();
								t=0.0;
								psquared=px*px+py*py+pz*pz;
								pdotr=px*x+py*y+pz*z;
								p0=sqrt(psquared+mass*mass);
								gamma=p0/mass;
								x=x-px*pdotr*(1.0-1.0/gamma)/psquared;
								y=y-py*pdotr*(1.0-1.0/gamma)/psquared;
								z=z-pz*pdotr*(1.0-1.0/gamma)/psquared;
							}
							
							tmp_particle->p[0]=p0;
							tmp_particle->p[1]=px;
							tmp_particle->p[2]=py;
							tmp_particle->p[3]=pz;
							tmp_particle->pid=pdg;
							acceptance->Smear(tmp_particle);
							tmp_particle->Setp0();
							tmp_particle->x[1]=x-(px/p0)*(t-taucompare);
							tmp_particle->x[2]=y-(py/p0)*(t-taucompare);
							tmp_particle->x[3]=z-(pz/p0)*(t-taucompare);
							tmp_particle->x[0]=taucompare;				
							
							accept=acceptance->OneParticleAcceptance(tmp_particle, eff);
							if(accept){
								cell_list->FindCell(tmp_particle,cell);
								if(cell!=NULL){
									if(pdg==PIDA)
										cell->partlist_a.push_back(tmp_particle);
									else
										cell->partlist_b.push_back(tmp_particle);
									tmp_particle=new Chbt_part;
									naccept+=1;
								}
							}
						}
					}
				} // end particles loop
				nrParticlesInEvent = 0;
			}
		}while(!feof(fptr_in));
		fclose(fptr_in);
	}//end of loop over files
	CLog::Info("readOSCAR: Naccept="+to_string(naccept)+"\n");
	delete tmp_particle;
}
