#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "msu_hbt/hbt.h"

using namespace std;

int main(int argc,char *argv[]){
	Crandy *randy=new Crandy(1234);
	double R,cf,x,y,z,root2=sqrt(2.0);
	unsigned int imc,NMC=1000000;
	string parsfilename;
	printf("Enter Rinv in fm: ");
	scanf("%lf",&R);
	if(argc!=2){
		CLog::Info("Usage hbt_fromGAUSS parameter_file_name_prefix\n");
		exit(1);
	}
	else{
		parsfilename="parameters/"+string(argv[1])+".txt";
	}
	CWaveFunction_pd_sqwell *wf=new CWaveFunction_pd_sqwell(parsfilename);
	//wf->PrintPhaseShifts();
	
	
	double psisquared,qmag,r,ctheta;
	for(qmag=2;qmag<201;qmag+=2){
		cf=0.0;
		for(imc=0.0;imc<NMC;imc++){
			x=R*root2*randy->ran_gauss();
			y=R*root2*randy->ran_gauss();
			z=R*root2*randy->ran_gauss();
			r=sqrt(x*x+y*y+z*z);
			ctheta=z/r;
			psisquared=wf->GetPsiSquared(qmag,r,ctheta);
			cf+=psisquared;
		}
		cf=cf/double(NMC);
		printf("%5.1f %g\n",qmag,cf);
	}
	
	
	return 0;
}
