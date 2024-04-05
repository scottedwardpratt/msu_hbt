#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "msu_hbt/hbt.h"

using namespace std;

int main(int argc,char *argv[]){
	string parsfilename;
	double qmag,r,KL;
	if(argc!=2){
		CLog::Info("Usage hbt_fromGAUSS parameter_file_name_prefix\n");
		exit(1);
	}
	else{
		parsfilename="parameters/"+string(argv[1])+".txt";
	}
	CWaveFunction_pp_schrod *wf=new CWaveFunction_pp_schrod(parsfilename);
	CKernel *kernel=new CKernel(parsfilename);
	kernel->Calc(wf);
	
	int ell,Lmax=kernel->GetLMAX();
	printf(" q         L=0       L=2       L=4       L=6\n");
	for(qmag=5;qmag<30;qmag+=10){
		printf("----- q=%g ------------\n",qmag);
		for(r=0.5;r<10;r+=1.0){
			printf("%5.2f ",r);
			for(ell=0;ell<=Lmax;ell+=2){
				KL=kernel->GetValue(ell,qmag,r);
				printf("%9.4f ",KL);
			}
			printf("\n");
		}
	}
	
	return 0;
}
