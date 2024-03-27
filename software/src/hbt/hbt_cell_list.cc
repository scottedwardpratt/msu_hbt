#include "msu_hbt/hbt.h"
Chbt_master *Chbt_cell_list::master=NULL;

using namespace std;
using namespace NMSUPratt;

Chbt_cell_list::Chbt_cell_list(CparameterMap *parmap){
	int ix,iy,iz;
	int inx,iny,inz;
	//QMAX=parmap->getD("QMAX",50.0);
	QMAX=parmap->getD("NQMAX",50)*parmap->getD("DQINV",2.0);
	double ma=master->wf->m1;
	double mb=master->wf->m2;
	double mu=ma*mb/(ma+mb);
	DRAPX=DRAPY=DRAPZ=QMAX/mu;

	/* This is the old way
	double PXMAXa=parmap->getD("PXMAXA",1000.0);
	double PYMAXa=parmap->getD("PYMAXA",1000.0);
	double PZMAXa=parmap->getD("PZMAXA",2000.0);
	double PXMAXb=parmap->getD("PXMAXB",1000.0);
	double PYMAXb=parmap->getD("PYMAXB",1000.0);
	double PZMAXb=parmap->getD("PZMAXB",2000.0);
	if(PXMAXa/ma > PXMAXb/mb)
		rapxmax=asinh(PXMAXb/mb);
	else
		rapxmax=asinh(PXMAXa/ma);
	
	if(PYMAXa/ma > PYMAXb/mb)
		rapymax=asinh(PYMAXb/mb);
	else
		rapymax=asinh(PYMAXa/ma);
	if(PZMAXa/ma > PZMAXb/mb)
		rapzmax=asinh(PZMAXb/mb);
	else
		rapzmax=asinh(PZMAXa/ma);

	NRAPX=2+2*floorl(rapxmax/DRAPX);
	NRAPY=2+2*floorl(rapymax/DRAPY);
	NRAPZ=2+2*floorl(rapzmax/DRAPZ);
	*/
	// This is the new way
	double PXMINa,PXMAXa,PXMINb,PXMAXb,PYMINa,PYMAXa,PYMINb,PYMAXb;
	double acceptance_ptmax_a,acceptance_ptmax_b,acceptance_ptmin_a,acceptance_ptmin_b;
	double acceptance_rapmax,acceptance_rapmin;
	double rapxmin_a,rapxmin_b,rapxmax_a,rapxmax_b,rapymin_a,rapymin_b,rapymax_a,rapymax_b;
	acceptance_rapmax=parmap->getD("ACCEPTANCE_YMAX",1.0);
	acceptance_rapmin=parmap->getD("ACCEPTANCE_YMIN",1.0);
	acceptance_ptmin_a=parmap->getD("ACCEPTANCE_PTMIN_A",0.0);
	acceptance_ptmin_b=parmap->getD("ACCEPTANCE_PTMIN_B",0.0);
	acceptance_ptmax_a=parmap->getD("ACCEPTANCE_PTMAX_A",0.0);
	acceptance_ptmax_b=parmap->getD("ACCEPTANCE_PTMAX_B",0.0);
	
	PXMINa=acceptance_ptmin_a;
	PXMAXa=acceptance_ptmax_a;
	rapxmin_a=asinh(PXMINa/ma);
	rapxmax_a=asinh(PXMAXa/ma);
	PXMINb=acceptance_ptmin_b;
	PXMAXb=acceptance_ptmax_b;
	rapxmin_b=asinh(PXMINb/mb);
	rapxmax_b=asinh(PXMAXb/mb);
	if(rapxmin_a<rapxmin_b)
		rapxmin=rapxmin_a;
	else
		rapxmin=rapxmin_b;
	if(rapxmax_a>rapxmax_b)
		rapxmax=rapxmax_a;
	else
		rapxmax=rapxmax_b;
	NRAPX=ceill((rapxmax-rapxmin)/DRAPX);
	rapxmax=rapxmin+NRAPX*DRAPX;
	
	PYMINa=acceptance_ptmin_a;
	PYMAXa=acceptance_ptmax_a;
	rapymin_a=asinh(PYMINa/ma);
	rapymax_a=asinh(PYMAXa/ma);
	PYMINb=acceptance_ptmin_b;
	PYMAXb=acceptance_ptmax_b;
	rapymin_b=asinh(PYMINb/mb);
	rapymax_b=asinh(PYMAXb/mb);
	if(rapymin_a<rapymin_b)
		rapymin=rapymin_a;
	else
		rapymin=rapymin_b;
	if(rapymax_a>rapymax_b)
		rapymax=rapymax_a;
	else
		rapymax=rapymax_b;
	NRAPY=ceill((rapymax-rapymin)/DRAPY);
	rapymax=rapxmin+NRAPX*DRAPY;
	
	rapzmin=acceptance_rapmin;
	rapzmax=acceptance_rapmax;
	NRAPZ=ceill((rapzmax-rapzmin)/DRAPZ);
	rapzmax=rapzmin+NRAPZ*DRAPZ;

	CLog::Info("NRAPX="+to_string(NRAPX)+", NRAPY="+to_string(NRAPY)+", NRAPZ="+to_string(NRAPZ)+"\n");
	
	cell.resize(NRAPX);
	for(ix=0;ix<NRAPX;ix++){
		cell[ix].resize(NRAPY);
		for(iy=0;iy<NRAPY;iy++){
			cell[ix][iy].resize(NRAPZ);
			for(iz=0;iz<NRAPZ;iz++){
				cell[ix][iy][iz]=new Chbt_cell();
			}
		}
	}
	for(ix=0;ix<NRAPX;ix++){
		for(iy=0;iy<NRAPY;iy++){
			for(iz=0;iz<NRAPZ;iz++){
				cell[ix][iy][iz]->neighbor.resize(3);
				for(inx=0;inx<3;inx++){
					cell[ix][iy][iz]->neighbor[inx].resize(3);
					for(iny=0;iny<3;iny++){
						cell[ix][iy][iz]->neighbor[inx][iny].resize(3);
						for(inz=0;inz<3;inz++){
							if((ix+inx<NRAPX && ix+inx-1>=0) 
								&&(iy+iny<NRAPY && iy+iny-1>=0)
									&&(iz+inz<NRAPZ && iz+inz-1>=0)){
										cell[ix][iy][iz]->neighbor[inx][iny][inz]=cell[ix+inx-1][iy+iny-1][iz+inz-1];
							}
						}
					}
				}
			}
		}
	}
}

void Chbt_cell_list::FindCell(Chbt_part *part,Chbt_cell *&cellptr){
	double px=part->psmear[1],py=part->psmear[2],pz=part->psmear[3];
	double E=part->psmear[0];
	double mass=sqrt(E*E-px*px-py*py-pz*pz);
	double rapx=asinh(px/sqrt(mass*mass+py*py+pz*pz));
	double rapy=asinh(py/sqrt(mass*mass+px*px+pz*pz));
	double rapz=asinh(pz/sqrt(mass*mass+py*py+px*px));
	cellptr=NULL;
	if(fabs(rapx)<rapxmax && fabs(rapy)<rapymax && fabs(rapz)<rapzmax){
		double drapx=rapx+NRAPX*DRAPX*0.5;
		int irapx=lrint(floor(drapx/DRAPX));
		if(irapx>=0 && irapx<NRAPX){
			double drapy=rapy+NRAPY*DRAPY*0.5;
			int irapy=lrint(floor(drapy/DRAPY));
			if(irapy>=0 && irapy<NRAPY){
				double drapz=rapz+NRAPZ*DRAPZ*0.5;
				int irapz=lrint(floor(drapz/DRAPZ));
				if(irapz>=0 && irapz<NRAPZ){
					cellptr=cell[irapx][irapy][irapz];
				}
			}
		}
	}	
}
