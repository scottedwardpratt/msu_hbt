#ifndef __HBT_H__
#define __HBT_H__

#include <list>
#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/parametermap.h"
#include "msu_coral/coral.h"
#include "msu_hbt/acceptance.h"

using namespace std;

class Chbt_part;
class Chbt_cell;
class Chbt_cell_list;
class Chbt_part;
class Chbt_resinfo;
class Chbt_CFs;
class Chbt_acceptance;

class Chbt_master{
public:
	string parsfilename_prefix;
	Chbt_master(string parmapfilename_prefix);
	int PIDA,PIDB;
	bool GAUSS;
	double GAUSS_Rx,GAUSS_Ry,GAUSS_Rz;
	CparameterMap parmap;
	void ReadOSCAR_1997();
	void ReadOSCAR_2003();
	CWaveFunction *wf;
	double GetCorrelationWeight(Chbt_part *parta,Chbt_part *partb);
	Chbt_acceptance *acceptance;
	void IncrementCFs(Chbt_part *parta,Chbt_part *partb);
	void CalcCFs();
	void CalcCFs_Gaussian();
	Chbt_cell_list *cell_list;
	Chbt_CFs *cfs;
	int nincrement,nsuccess;
	Crandy *randy;
	char message[CLog::CHARLENGTH];
};

class Chbt_cell_list{
public:
	Chbt_cell_list(CparameterMap *parmap);
	int NRAPX,NRAPY,NRAPZ;
	double DRAPX,DRAPY,DRAPZ;
	double rapxmax,rapymax,rapzmin,rapzmax;
	double QMAX;
	vector<vector<vector<Chbt_cell *> >> cell;
	void FindCell(Chbt_part *part,Chbt_cell *&cell);
	void Add2List(Chbt_part &parta);
	static Chbt_master *master;
	char message[CLog::CHARLENGTH];
};

class Chbt_cell{
public:
	Chbt_cell();
	vector<Chbt_part *> partlist_a;
	vector<Chbt_part *> partlist_b;
	vector<vector<vector<Chbt_cell *>>> neighbor;
	//vector<vector<int>> nlist[2][2][2];
};

class Chbt_resinfo{
public:
	Chbt_resinfo();
	int pida,pidb;
	double massa,massb;
};

class Chbt_CFs{
public:
	Chbt_CFs(CparameterMap *parmap);
	int NQINV;
	double DQINV;
	vector<double> C_of_qinv;
	vector<int> denom_of_qinv;
	void PrintC_of_qinv();
	void WriteC_of_qinv();
	void WriteC3D();
	bool XSYM,YSYM,ZSYM;
	
	C3DArray *threed_num,*threed_den;
	double Q3DMAX,DELQ3D;
	int NQ3D;
	static Chbt_master *master;
	
	char message[CLog::CHARLENGTH];
};

class Chbt_part{
public:
	void Print(){
		double mass=sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
		printf("x=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3]);
		printf("mass=%g, p=(%g,%g,%g,%g)\n",mass,p[0],p[1],p[2],p[3]);
		printf("psmear=(%g,%g,%g,%g)\n",psmear[0],psmear[1],psmear[2],psmear[3]);
	}
	int pid;
	FourVector x;
	FourVector p; // true momentum
	FourVector psmear; // smeared momentum, due to resolution
	double mass;
	void Setp0(){
		if(abs(pid)==2112)
			mass=NeutronMass;
		else if(abs(pid)==2212)
			mass=ProtonMass;
		else if(abs(pid)==211)
			mass=PionMass;
		p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
		psmear[0]=sqrt(mass*mass+psmear[1]*psmear[1]+psmear[2]*psmear[2]+psmear[3]*psmear[3]);
	}
	char message[CLog::CHARLENGTH];
};

#endif
