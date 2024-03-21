#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "msu_hbt/hbt.h"

using namespace std;

int main(int argc,char *argv[]){
	string parsfilename;
	if(argc!=2){
		printf("Usage hbt_fromGAUSS parameter_file_name_prefix\n");
		exit(1);
	}
	else{
		parsfilename=argv[1];
	}
	Chbt_master *hbt_master=new Chbt_master(parsfilename);
	hbt_master->CalcCFs_Gaussian();
	hbt_master->cfs->PrintC_of_qinv();
	//hbt_master->cfs->WriteC_of_qinv();
	//msu_hbt_master->cfs->WriteC3D("threed_output_gauss");
	return 0;
}
