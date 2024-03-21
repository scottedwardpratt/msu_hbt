#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "hades_hbt/hades_hbt.h"

using namespace std;

int main(int argc,char *argv[]){
	string parsfilename_prefix;
	if(argc!=2){
		printf("Usage hades_hbt_fromGAUSS parameter_file_prefix\n");
		exit(1);
	}
	else{
		parsfilename_prefix=argv[1];
	}
	Chades_hbt_master *hades_hbt_master=new Chades_hbt_master(parsfilename_prefix);
	hades_hbt_master->ReadOSCAR_1997();
	//hades_hbt_master->ReadOSCAR_2003();
	hades_hbt_master->CalcCFs();
	hades_hbt_master->cfs->PrintC_of_qinv();
	hades_hbt_master->cfs->WriteC_of_qinv();
	hades_hbt_master->cfs->WriteC3D();
	return 0;
}

