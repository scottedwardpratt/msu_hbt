#include "msu_coral/wavefunction.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/misc.h"
using namespace NMSUPratt;

CWaveFunction_lambdalambdaantiparspin_phaseshift::CWaveFunction_lambdalambdaantiparspin_phaseshift(string  parsfilename){
	int iq,ichannel;
	double q;

	ParsInit(parsfilename);

	m1=1115.7;
	m2=m1;
	IDENTICAL=1;
	q1q2=0;
	nchannels=1;
	ellmax=0;
	InitArrays();

	ell[0]=0;

	InitWaves();

	channelweight[0]=1.0;

	GetPhaseshifts();

	for(ichannel=0;ichannel<nchannels;ichannel++){
		for(iq=0;iq<nqmax;iq++){
			q=qarray[iq];
			Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
				-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
				+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
			Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
				/(4.0*PI*pow(epsilon,3));
		}
	}
}

double CWaveFunction_lambdalambdaantiparspin_phaseshift::CalcPsiSquared(int iq,double r,double ctheta){
	double psisquared,x,dpsi2;
	complex<double> psisymm,psiantisymm,psia,psib,hstar0;
	int ichannel;
	const double ROOT2=sqrt(2.0);
	double q;

	if(iq>=nqmax){
		psisquared=1.0;
	}
	else{
		q=qarray[iq];
		psia=planewave[iq]->planewave(r,ctheta);
		psib=planewave[iq]->planewave(r,-ctheta);
		psisymm=(psia+psib)/ROOT2;
		psiantisymm=(psia-psib)/ROOT2;
		if(STRONG==1){
			x=q*r/HBARC;
			hstar0=partwave[0][iq]->GetPhiIncoming(r)/x;
			if(r<epsilon){
				psisquared=0.5*real(psisymm*conj(psisymm))
					+0.5*real(psiantisymm*conj(psiantisymm));
				for(ichannel=0;ichannel<nchannels;ichannel++){
					dpsi2=channelweight[ichannel]*2.0*PI*Wepsilon[ichannel][iq]
						*pow(HBARC,3)/(q*q);
					psisquared+=dpsi2;
				}
			}
			else{
				psisymm+=0.5*hstar0*ROOT2*(Misc::ceiphi(-2.0*delta[0][iq])-1.0);
				psisquared=0.5*real(psisymm*conj(psisymm))
					+0.5*real(psiantisymm*conj(psiantisymm));
			}
		}
		else psisquared=0.5*real(psisymm*conj(psisymm))
			+0.5*real(psiantisymm*conj(psiantisymm));
	}
	psisquared*=RelativisticCorrection(r,iq);
	return psisquared;
}

void CWaveFunction_lambdalambdaantiparspin_phaseshift::GetPhaseshifts(){
	double q,q0;
	double tandelta,alength,dtandeltadq;
	double MH0,M,EH0,GammaH0,lambda=500.0;
	int iq;

	CLog::Info("Enter the energy of the H0 above the 2Lambda threshold in MeV : ");
	scanf("%lf",&EH0);
	CLog::Info("Enter the width of the H0 in MeV : ");
	scanf("%lf",&GammaH0);

	// Scattering length and effective range parameters
	CLog::Info("Enter scattering length in fm (Do not include effect of H0): ");
	scanf("%lf",&alength);
	lambda=500.0; // Arbitrary choice, returns delta to zero at large q
	//reff=2*HBARC*HBARC/(lambda*lambda*a); // Effective range (not used)

	MH0=m1+m2+EH0;
	q0=sqrt(0.25*MH0*MH0-m1*m1);
	snprintf(message,CLog::CHARLENGTH,"Resonance occurs at q=%g, mesh goes q=%g\n",
	q0,nqmax*delq);
	CLog::Info(message);

	for(iq=0;iq<nqmax;iq++){
		q=qarray[iq];
		M=2.0*sqrt(q*q+m1*m1);
		tandelta=(alength*q/HBARC)/(1.0+q*q/(lambda*lambda));
		dtandeltadq=(tandelta/q)*(1.0-q*q/(lambda*lambda))
			/(1.0+q*q/(lambda*lambda));
		tandelta+=0.5*(q/q0)*GammaH0/(MH0-M);
		dtandeltadq+=(tandelta/q)*(1+(4.0*q*q/M)/(MH0-M));
		delta[0][iq]=atan(tandelta);
		ddeltadq[0][iq]=dtandeltadq*pow(cos(delta[0][iq]),2);
	}

}
