#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "msu_hbt/hbt.h"

using namespace std;

int main(int argc,char *argv[]){
	string parsfilename_prefix;
	if(argc!=2){
		printf("Usage hbt_fromGAUSS parameter_file_prefix\n");
		exit(1);
	}
	else{
		parsfilename_prefix=argv[1];
	}
	Chbt_master *hbt_master=new Chbt_master(parsfilename_prefix);
	hbt_master->ReadOSCAR_1997();
	//hbt_master->ReadOSCAR_2003();
	hbt_master->CalcCFs();
	hbt_master->cfs->PrintC_of_qinv();
	hbt_master->cfs->WriteC_of_qinv();
	hbt_master->cfs->WriteC3D();
	return 0;
}

