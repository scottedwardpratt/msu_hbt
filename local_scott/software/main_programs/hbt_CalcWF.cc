#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "msu_hbt/hbt.h"

using namespace std;

int main(int argc,char *argv[]){
	Crandy *randy=new Crandy(1234);
	double R,cf,unc,x,y,z,root2=sqrt(2.0);
	unsigned int imc,NMC=100000;
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
	CWaveFunction_pp_schrod *wf=new CWaveFunction_pp_schrod(parsfilename);
	double psisquared,qmag,r,ctheta;
	for(qmag=2;qmag<81;qmag+=2){
		cf=unc=0.0;
		for(imc=0.0;imc<NMC;imc++){
			x=R*root2*randy->ran_gauss();
			y=R*root2*randy->ran_gauss();
			z=R*root2*randy->ran_gauss();
			r=sqrt(x*x+y*y+z*z);
			ctheta=z/r;
			psisquared=wf->GetPsiSquared(qmag,r,ctheta);
			cf+=psisquared;
			unc+=psisquared*psisquared;
		}
		cf=cf/double(NMC);
		unc=unc/double(NMC);
		unc=sqrt(unc-cf*cf)/sqrt(double(NMC));
		printf("%5.1f %g +/- %g\n",qmag,cf,unc);
	}
	
	return 0;
}
