#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/sf.h"
#include "msu_commonutils/randy.h"
using namespace std;
using namespace NMSUPratt;

CWaveFunction_pd_sqwell::CWaveFunction_pd_sqwell(string parsfilename) : CWaveFunction(){
	CLog::Info("Beware! The p-d wavefunction was tuned to match phase shifts which were only measured for q<100.\n Also, the p-d system becomes inelastic (deuteron breaks up) above q=52 MeV/c\n So this treatment is pretty questionable for q>50!\n");
	// Interaction fit to phaseshifts from T.C. Black et al., PLB 471, p. 103-107 (1999).
  ParsInit(parsfilename);
	Crandy randy(1234);

  m1=ProtonMass; 
  m2=ProtonMass+NeutronMass-2.224;
  IDENTICAL=0;

  q1q2=1;
  nchannels=4;
	
  ellmax=1;
  InitArrays();
  CLog::Info("Arrays Initialized\n");

  ell[0]=0;
  ell[1]=0;
	ell[2]=1;
	ell[3]=1;
  InitWaves();

  nwells=new int[nchannels];
  nwells[0]=3;
  nwells[1]=1;
	nwells[2]=1;
	nwells[3]=2;

  SquareWell_MakeArrays();
	
  a[0][0]=3.23445; a[0][1]=5.57680; a[0][2]=7.71351; // L=0, S=1/2
  a[1][0]=4.99965; // L=0, S=3/2
	//a[2][0]=4.99965; // L=0, S=3/2
	//a[3][0]= 0.54790; a[3][1]=7.17907; // L=0, S=3/2
	

  V0[0][0]=-28.002; V0[0][1]=37.283; V0[0][2]=-5.829;
  V0[1][0]=108.203;
	
	//V0[2][0]=4.0;
	//V0[3][0]=12.63402; V0[3][1]=-2.06386;
	
	a[2][0]=7.93404; 
	V0[2][0]=0.373212; 
	a[3][0]=0.426206; a[3][1]=7.17018; 
	V0[3][0]=12.109; V0[3][1]=-2.01329; 
	//V0[0][0]=V0[0][1]=V0[0][2]=V0[1][0]=0.0;

  SquareWell_Init();
	
	bool atest;
	double chisquared,dd;
	double Vstep=0.3,astep=0.1;
	int ichannel,iq,iwell;
	double bestchisquared=1.0E99;
	vector<vector<double>> Vbest,abest;
	double expdata_e243[4]={0,0,-7.1,22.0};
	double expdata_e363[4]={0,0,-7.3,24.0};
	double expdata_e675[4]={0,0,-5.8,29.5};
	double expdata_e1082[4]={0,0,-3.1,31.0};
	Vbest.resize(nchannels);
	abest.resize(nchannels);
	for(ichannel=0;ichannel<nchannels;ichannel++){
		Vbest[ichannel].resize(nwells[ichannel]);
		abest[ichannel].resize(nwells[ichannel]);
		for(iwell=0;iwell<nwells[ichannel];iwell++){
			Vbest[ichannel][iwell]=V0[ichannel][iwell];
			abest[ichannel][iwell]=a[ichannel][iwell];
		}
	}
	// uses delq=10
	for(int imc=0;imc<1000;imc++){
		SquareWell_Init();
		chisquared=0.0;
		for(ichannel=3;ichannel<nchannels;ichannel++){
			iq=4;
			dd=delta[ichannel][iq]*180.0/PI;
			chisquared+=pow(dd-expdata_e243[ichannel],2);
			iq=5;
			dd=delta[ichannel][iq]*180.0/PI;
			chisquared+=pow(dd-expdata_e363[ichannel],2);
			iq=7;
			dd=delta[ichannel][iq]*180.0/PI;
			chisquared+=pow(dd-expdata_e675[ichannel],2);
			iq=9;
			dd=delta[ichannel][iq]*180.0/PI;
			chisquared+=pow(dd-expdata_e1082[ichannel],2);	
		}
		//printf("chisquared=%g\n",chisquared);
		if(chisquared<bestchisquared){
			printf("success!!!!\n");
			for(ichannel=3;ichannel<nchannels;ichannel++){
				for(iwell=0;iwell<nwells[ichannel];iwell++){
					Vbest[ichannel][iwell]=V0[ichannel][iwell];
					abest[ichannel][iwell]=a[ichannel][iwell];
				}
			}
			bestchisquared=chisquared;
			printf("bestchisquared=%g\n",chisquared);
		}
		for(ichannel=3;ichannel<nchannels;ichannel++){
			for(iwell=0;iwell<nwells[ichannel];iwell++){
				V0[ichannel][iwell]=Vbest[ichannel][iwell]+Vstep*randy.ran_gauss();
				atest=false;
				do{
					a[ichannel][iwell]=abest[ichannel][iwell]+astep*randy.ran_gauss();
					if(a[ichannel][iwell]<0.0)
						atest=true;
					if(iwell>0 && a[ichannel][iwell]<a[ichannel][iwell-1]){
						atest=true;
					}
				}while(atest);
			}
		}
	}
	for(ichannel=3;ichannel<nchannels;ichannel++){
		for(iwell=0;iwell<nwells[ichannel];iwell++){
			printf("a[%d][%d]=%g; ",ichannel,iwell,abest[ichannel][iwell]);
		}
		printf("\n");
		for(iwell=0;iwell<nwells[ichannel];iwell++){
			printf("V0[%d][%d]=%g; ",ichannel,iwell,Vbest[ichannel][iwell]);
		}
		printf("\n");
		
	}
	for(ichannel=0;ichannel<nchannels;ichannel++){
		for(iwell=0;iwell<nwells[ichannel];iwell++){
			V0[ichannel][iwell]=Vbest[ichannel][iwell];
			a[ichannel][iwell]=abest[ichannel][iwell];
		}
	}
	SquareWell_Init();
	
		
}

CWaveFunction_pd_sqwell::~CWaveFunction_pd_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_pd_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,theta=acos(ctheta),x;
  double q=GetQ(iq);
  complex<double> psi,psia,Xlm00,Xlm10;
  psia=planewave[iq]->planewave(r,ctheta);
  complex<double> DelPhi[4];

  SquareWell_GetDelPhi(iq,r,DelPhi);
	x=q*r/HBARC;
  Xlm00=sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)/x;
	Xlm10=ci*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,1,theta,0.0)/x;
	// for S=1/2
  psi=psia+Xlm00*DelPhi[0]+Xlm10*DelPhi[2];
  psisquared=real(psi*conj(psi))/3.0;
	// for S=3/2
  psi=psia+Xlm00*DelPhi[1]+Xlm10*DelPhi[3];
  psisquared+=2.0*real(psi*conj(psi))/3.0;

	//psisquared*=RelativisticCorrection(r,iq);
  return psisquared;

}
