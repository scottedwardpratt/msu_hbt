#include "msu_hbt/hbt.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace NMSUPratt;

void Chbt_master::ReadOSCAR_2003(){
	Chbt_cell *cell;
	Chbt_part *tmp_particle=new Chbt_part();
	string directory=parmap.getS("OSCAR_DIRNAME","oscar_files");
	double taucompare=parmap.getD("OSCAR_TAUCOMPRE",20.0);
	string filenamefilenames=parmap.getS("OSCAR_FILENAME_FILENAMES","oscarfilenames_URQMD.txt");

	double BMIN=parmap.getD("OSCAR_BMIN",-0.1);
	double BMAX=parmap.getD("OSCAR_BMAX",100.0);
	
	// opening the files from first -> last:
	ifstream f_in;

	string dust, line; // dummy variable
	double t, x, y, z, mass, p0, px, py, pz, charge, bim;
	int pid,pdg,nparts=0;
	int nrParticlesInEvent;
	
	list<string> oscar_filenames;
	list<string>::iterator fiter;
	string filename;
	char dummy[200];
	FILE  *fptr_filenames=fopen(filenamefilenames.c_str(),"r");
	do{
		fscanf(fptr_filenames,"%s",dummy);
		filename=dummy;
		oscar_filenames.push_back(filename);
	}while(!feof(fptr_filenames));
	fclose(fptr_filenames);
	
	
	
	for(fiter=oscar_filenames.begin();fiter!=oscar_filenames.end();++fiter){
		string infile_name=*fiter;
		// open input file
		CLog::Info("Reading "+infile_name+"\n");
		f_in.open( infile_name );
		if( f_in.fail() )
			CLog::Fatal( "cannot open input file" );
		//getting rid of some lines
		for(int i = 0; i < 3; i++)
			getline( f_in, line );
		while (!f_in.eof()){
			for(int i = 0; i < 4; i++) f_in >> dust;
			f_in >> nrParticlesInEvent;
			for(int i = 0; i < nrParticlesInEvent; i++){//reading particles in event loop
				f_in >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pdg >> pid >> charge;
				mass*=1000.0; p0*=1000.0; px*=1000.0; py*=1000.0; pz*=1000.0;
				bool accept; double eff;
				if(pdg == PIDA){
					tmp_particle->p[0]=p0;
					tmp_particle->p[1]=px; tmp_particle->p[2]=py; tmp_particle->p[3]=pz;
					tmp_particle->mass=mass;
					tmp_particle->pid=pdg;
					acceptance->Smear(tmp_particle);
					tmp_particle->Setp0();
					//tmp_particle->x[0]=t;
					tmp_particle->x[1]=x-(px/p0)*(t-taucompare);
					tmp_particle->x[2]=y-(py/p0)*(t-taucompare);
					tmp_particle->x[3]=z-(pz/p0)*(t-taucompare);
					tmp_particle->x[0]=taucompare;
					
					accept=acceptance->OneParticleAcceptance(pdg,tmp_particle, eff);
					if(accept){
						cell_list->FindCell(tmp_particle,cell);
						if(cell!=NULL){
							cell->partlist_a.push_back(tmp_particle);
							tmp_particle=new Chbt_part;
							nparts+=1;
						}
					}
				}
				else if(pdg == PIDB ){
					accept=acceptance->OneParticleAcceptance(pdg, tmp_particle, eff);
					tmp_particle->p[0]=p0;
					tmp_particle->p[1]=px; tmp_particle->p[2]=py; tmp_particle->p[3]=pz;
					tmp_particle->mass=mass;
					tmp_particle->pid=pdg;
					acceptance->Smear(tmp_particle);
					tmp_particle->Setp0();
					tmp_particle->x[1]=x-(px/p0)*(t-taucompare);
					tmp_particle->x[2]=y-(py/p0)*(t-taucompare);
					tmp_particle->x[3]=z-(pz/p0)*(t-taucompare);
					tmp_particle->x[0]=taucompare;
					
					accept=acceptance->OneParticleAcceptance(pdg,tmp_particle, eff);
					if(accept){
						cell_list->FindCell(tmp_particle,cell);
						if(cell!=NULL){
							cell->partlist_b.push_back(tmp_particle);
							tmp_particle=new Chbt_part;
							nparts+=1;
						}
					}
				}
			} // end particles loop
			for(int i = 0; i < 6; i++) f_in >> dust;
			f_in >> bim;
			if(bim<BMIN || bim>BMAX){
				snprintf(message,CLog::CHARLENGTH,"bim=%g, but BMIN=%g and BMAX=%g\n",bim,BMIN,BMAX);
				CLog::Info(message);
			}
			for(int i = 0; i < 2; i++) f_in >> dust;
		}
		CLog::Info("readOSCAR: read in "+to_string(nparts)+" parts\n");
	}//end of loop over files
	delete tmp_particle;
}
