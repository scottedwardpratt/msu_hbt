#ifndef __HBT_ACCEPTANCE_H__#ifndef __HBT_ACCEPTANCE_H__
#define __HBT_ACCEPTANCE_H__
#include "msu_hbt/hbt.h"
using namespace NMSUPratt;

class Chbt_part;
class Chbt_master;

class Chbt_acceptance{
public:
	double pTmax,pTmin,thetamin,thetamax,ymax,ymin;
	char message[CLog::CHARLENGTH];
	static Chbt_master *master;
	CparameterMap *parmap;
	Crandy *randy;
	FourVector ucm;
	Chbt_acceptance(){
		//
	}
	void Init(CparameterMap *parmap);
	
	Chbt_acceptance(CparameterMap *parmap){
		Init(parmap);
	}
	
	virtual bool OneParticleAcceptance(int pid,Chbt_part *part,double &efficiency){
		// this is a dummy function
		return true;
	}
	
	virtual bool TwoParticleAcceptance(Chbt_part *parta,Chbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency){
		// this is a dummy function
		return true;
	}
	
	
	virtual void Smear(Chbt_part *part){
		// this is a dummy function
	}
	
	
};

class Chbt_acceptance_nosmear : public Chbt_acceptance{
public:
	Chbt_acceptance_nosmear(CparameterMap *parmap){
		Init(parmap);
	}
	bool OneParticleAcceptance(int pid,Chbt_part *part,double &efficiency);
	bool TwoParticleAcceptance(Chbt_part *parta,Chbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency);
	void Smear(Chbt_part *part);
};

class Chbt_acceptance_smear : public Chbt_acceptance{
public:
	Chbt_acceptance_smear(CparameterMap *parmap){
		Init(parmap);
	}
	bool OneParticleAcceptance(int pid,Chbt_part *part,double &efficiency);
	bool TwoParticleAcceptance(Chbt_part *parta,Chbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency);
	void Smear(Chbt_part *part);
};

class Chbt_acceptance_smear_maria : public Chbt_acceptance{
public:
	Chbt_acceptance_smear_maria(CparameterMap *parmap){
		Init(parmap);
	}
	bool OneParticleAcceptance(int pid,Chbt_part *part,double &efficiency);
	bool TwoParticleAcceptance(Chbt_part *parta,Chbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency);
	void Smear(Chbt_part *part);
};

#endif
