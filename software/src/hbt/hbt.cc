#include "msu_hbt/hbt.h"
using namespace NMSUPratt;

Chbt_master::Chbt_master(string parsfilename_prefix_set){
	parsfilename_prefix=parsfilename_prefix_set;
	string parsfilename="parameters/"+parsfilename_prefix+".txt";
	parmap.ReadParsFromFile(parsfilename);
	parmap.PrintPars();
	string logfilename=parmap.getS("LOG_FILENAME","log.txt");
	//string coralpars_filename=parmap.getS("CORALPARS_FILENAME","parameters/coralpars.txt");
	//parmap.ReadParsFromFile(coralpars_filename);
	CLog::Init(logfilename);
	PIDA=parmap.getI("PIDA",211);
	PIDB=parmap.getI("PIDB",211);
	GAUSS=parmap.getB("GAUSS",false);
	if(GAUSS){
		GAUSS_Rx=parmap.getD("GAUSS_RX",2.7);
		GAUSS_Ry=parmap.getD("GAUSS_RY",2.7);
		GAUSS_Rz=parmap.getD("GAUSS_RZ",2.7);
	}
	if(PIDA==PIDB){
		parmap.set("XYZSYM",true);
	}	
	
	if((PIDA==2212 && PIDB==2212) || (PIDA==-2212 && PIDB==-2212)){
		wf=new CWaveFunction_pp_schrod(parsfilename);
	}
	else if((PIDA==211 && PIDB==211) || (PIDA==-211 && PIDB==-211)){
		wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);
	}
	else if((abs(PIDA)==2212 && fabs(PIDB)==22122112)){
		wf=new CWaveFunction_pd_sqwell(parsfilename);
	}
	else{
		CLog::Fatal("Cannot recognize PIDA="+to_string(PIDA)+" or PIDB="+to_string(PIDB)+"\n");
	}
	
	Chbt_cell_list::master=this;
	cell_list=new Chbt_cell_list(&parmap);
	
	Chbt_CFs::master=this;
	cfs=new Chbt_CFs(&parmap);
	
	Chbt_acceptance::master=this;
	string smearstring=parmap.getS("HBT_SMEARSTRING","smear");
	if(smearstring=="smear"){
		acceptance=new Chbt_acceptance_smear(&parmap);
	}
	else if(smearstring=="smear_maria"){
		acceptance=new Chbt_acceptance_smear_maria(&parmap);
	}
	else if(smearstring=="nosmear"){
		acceptance=new Chbt_acceptance_nosmear(&parmap);
	}
	else{
		CLog::Fatal("smearstring="+smearstring+" is not recognized\n");
	}
	acceptance->PIDA=PIDA;
	acceptance->PIDB=PIDB;
	
	randy=new Crandy(-12345);
	acceptance->randy=randy;
		
}

void Chbt_master::CalcCFs(){
	int icx,icy,icz,inx,iny,inz,ia,ib,na,nb;
	int inxmin,inymin,inzmin,ibmin;
	int natot=0;
	nincrement=nsuccess=0;
	inxmin=0;
	if(PIDA==PIDB)
		inxmin=1;
	Chbt_cell *cella,*cellb;
	Chbt_part *parta,*partb;
	for(icx=0;icx<cell_list->NRAPX;icx++){
		CLog::Info("icx="+to_string(icx)+" out of NRAPX="+to_string(cell_list->NRAPX)+"\n");
		for(icy=0;icy<cell_list->NRAPY;icy++){
			for(icz=0;icz<cell_list->NRAPZ;icz++){
				cella=cell_list->cell[icx][icy][icz];
				na=cella->partlist_a.size();
				natot+=na;
				for(ia=0;ia<na;ia++){
					parta=cella->partlist_a[ia];
					for(inx=inxmin;inx<3;inx++){
						inymin=0;
						if(PIDA==PIDB && inx==1)
							inymin=1;
						for(iny=inymin;iny<3;iny++){
							inzmin=0;
							if(PIDA==PIDB && inx==1 && iny==1)
								inzmin=1;
							for(inz=inzmin;inz<3;inz++){
								cellb=cella->neighbor[inx][iny][inz];
								if(cellb!=NULL){
									ibmin=0;
									if(PIDA==PIDB && inx==1 && iny==1 && inz==1){
										ibmin=ia+1;
									}
									if(PIDA==PIDB){
										nb=cellb->partlist_a.size();
									}
									else
										nb=cellb->partlist_b.size();
									for(ib=ibmin;ib<nb;ib++){
										if(PIDA==PIDB)
											partb=cellb->partlist_a[ib];
										else
											partb=cellb->partlist_b[ib];
										nincrement+=1;
										IncrementCFs(parta,partb);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	CLog::Info("nincrement="+to_string(nincrement)+", nincrement/npairs_tot="
		+to_string(2.0*double(nincrement)/(double(natot)*double(natot-1)))+"\n");//, nwf/nincrement="+to_string(double(nsuccess)/double(nincrement));
}

void Chbt_master::CalcCFs_Gaussian(){
	double Rx=parmap.getD("GAUSS_RX",3.0);
	double Ry=parmap.getD("GAUSS_RY",3.0);
	double Rz=parmap.getD("GAUSS_RZ",3.0);
	double x,y,z,qx,qy,qz,q,r,ctheta,weight,root2=sqrt(2.0);
	int imc,NMC=parmap.getI("NMC_GAUSSIAN",1000);
	int iq,iqx,iqy,iqz,isx,isy,isz,nsx=2,nsy=2,nsz=2;
	if(cfs->XSYM)
		nsx=1;
	if(cfs->YSYM)
		nsy=1;
	if(cfs->ZSYM)
		nsz=1;
	for(iqx=0;iqx<cfs->NQ3D;iqx++){
		for(isx=0;isx<nsx;isx++){
			for(iqy=0;iqy<cfs->NQ3D;iqy++){
				for(isy=0;isy<nsy;isy++){
					for(iqz=0;iqz<cfs->NQ3D;iqz++){
						for(isz=0;isz<nsz;isz++){
							for(imc=0;imc<NMC;imc++){
								qx=cfs->DELQ3D*(iqx+randy->ran());
								qy=cfs->DELQ3D*(iqy+randy->ran());
								qz=cfs->DELQ3D*(iqz+randy->ran());
								if(isx>0)
									qx=-qx;
								if(isy>0)
									qy=-qy;
								if(isz>0)
									qz=-qz;
								q=sqrt(qx*qx+qy*qy+qz*qz);
								weight=1.0;
								if(q<cell_list->QMAX){
									if(q<cfs->DQINV*cfs->NQINV){
										x=root2*Rx*randy->ran_gauss();
										y=root2*Ry*randy->ran_gauss();
										z=root2*Rz*randy->ran_gauss();
										r=sqrt(x*x+y*y+z*z);
										ctheta=(qx*x+qy*y+qz*z)/(q*r);
										weight=wf->GetPsiSquared(q,r,ctheta);
										iq=lrint(floor(q/cfs->DQINV));
										cfs->C_of_qinv[iq]+=weight;
										cfs->denom_of_qinv[iq]+=1;
										if(fabs(qx)<cfs->Q3DMAX && fabs(qy)<cfs->Q3DMAX && fabs(qz)<cfs->Q3DMAX){
											nsuccess+=1;
											cfs->threed_num->IncrementElement(qx,qy,qz,weight);
											cfs->threed_den->IncrementElement(qx,qy,qz,1.0);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Chbt_master::IncrementCFs(Chbt_part *parta,Chbt_part *partb){
	double qinv,r,ctheta,efficiency,weight=1.0;
	int iq;
	
	double qout,qlong,qside,deleta,dely,delphi,qinv_smeared;
	
	if(acceptance->TwoParticleAcceptance(parta,partb,efficiency)){
		Misc::outsidelong(parta->psmear,partb->psmear,qinv_smeared,qout,qside,qlong,deleta,dely,delphi);
		wf->getqrctheta(parta->p,parta->x,partb->p,partb->x,qinv,r,ctheta);
		if(r>1.0E-8){
			if(qinv_smeared<cell_list->QMAX){
				nsuccess+=1;
				if(qinv_smeared<cfs->DQINV*cfs->NQINV){
					weight=wf->GetPsiSquared(qinv,r,ctheta);
					if(weight!=weight){
						parta->Print();
						partb->Print();
						CLog::Fatal("weight=Nan\n");
					}
				}
			}
		}
	
		if(fabs(qout)<cfs->Q3DMAX && fabs(qside)<cfs->Q3DMAX && fabs(qlong)<cfs->Q3DMAX){
			cfs->threed_num->IncrementElement(qout,qlong,qside,weight*efficiency);
			cfs->threed_den->IncrementElement(qout,qlong,qside,efficiency);
		}
		iq=lrint(floor(qinv_smeared/cfs->DQINV));
		if(iq<int(cfs->C_of_qinv.size())){
			cfs->C_of_qinv[iq]+=weight*efficiency;
			cfs->denom_of_qinv[iq]+=efficiency;
		}
	}
	
}