#include "hades_hbt/hades_hbt.h"
using namespace std;

Chades_hbt_master *Chades_hbt_acceptance::master=NULL;

void Chades_hbt_acceptance::Init(CparameterMap *parmap){
	thetamin=parmap->getD("HADES_ACCEPTANCE_THETAMIN",18.0);
	thetamax=parmap->getD("HADES_ACCEPTANCE_THETAMAX",85.0);
	pTmin=parmap->getD("HADES_ACCEPTANCE_PTMIN",0.0);
	pTmax=parmap->getD("HADES_ACCEPTANCE_PTMAX",2500.0);
	ymin=parmap->getD("HADES_ACCEPTANCE_YMIN",-5.0);
	ymax=parmap->getD("HADES_ACCEPTANCE_YMAX",5.0);
	
	double ebeam=master->parmap.getD("HADES_BEAM_ENERGY_GEV",1.26); // KE of beam in lab frame per nucleon
	int    Abeam=master->parmap.getI("HADES_BEAM_A",197);
	int    Atarget=master->parmap.getI("HADES_TARGET_A",197);
	double m=0.931,Mtarget,Mbeam,roots,Pbeam,Ebeam;   // roughly account for binding energy
	Mbeam=Abeam*m;
	Mtarget=Atarget*m;
	Ebeam=(m+ebeam)*Abeam;
	Pbeam=sqrt(Ebeam*Ebeam-Mbeam*Mbeam);
	roots=sqrt((Ebeam+Mtarget)*(Ebeam+Mtarget)-Pbeam*Pbeam);
	ucm[0]=(Ebeam+Mtarget)/roots;
	ucm[1]=ucm[2]=0.0;
	ucm[3]=sqrt(ucm[0]*ucm[0]-1.0);
}

bool Chades_hbt_acceptance_nosmear::OneParticleAcceptance(int pid, Chades_hbt_part *part, double &efficiency){
	//// === HADES acceptance === ////
	efficiency=0.0;
	FourVector Plab;
	Misc::Boost(ucm,part->psmear,Plab);
  
	//momentum:
	int charge = pid / abs(pid);
	double pmag = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]+Plab[3]*Plab[3]);
	if(charge > 0 && pmag < 50.0)
		return false;
	if(charge < 0 && pmag < 105.0)
		return false;

	//pT:
	double pT = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]);
	if(pT<pTmin && pT>pTmax)
		return false;

	//theta:
	double theta = atan(pT/Plab[3])*(180.0/PI);
	if(theta < thetamin || theta > thetamax)
		return false;

	//rapidity:
	double rapidity = 0.5 * log((Plab[0]+Plab[3])/(Plab[0]-Plab[3]));
	if(rapidity > ymax || rapidity < ymin)
		return false;

	//// === Efficiency === ////
	efficiency=1.0;
	
	return true;
}

void Chades_hbt_acceptance_nosmear::Smear(Chades_hbt_part *part){
	part->psmear[1]=part->p[1];
	part->psmear[2]=part->p[2];
	part->psmear[3]=part->p[3];
	part->Setp0();
}

bool Chades_hbt_acceptance_nosmear::TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
double &efficiency){
    
	//// === HADES acceptance === ////
	bool acc=true;
	efficiency=1.0;

	FourVector PlabA;
	Misc::Boost(ucm,parta->psmear,PlabA);
	FourVector PlabB;
	Misc::Boost(ucm,partb->psmear,PlabB);
    
	double ptA = sqrt(PlabA[1]*PlabA[1] + PlabA[2]*PlabA[2]);
	double ptB = sqrt(PlabB[1]*PlabB[1] + PlabB[2]*PlabB[2]);
	double thetaA = atan(ptA/PlabA[3]);
	double thetaB = atan(ptB/PlabB[3]);
    
	if(fabs(thetaA-thetaB)<0.04 && fabs(delphi)<0.25){ // merging suppression
		efficiency=0.0;
		acc=false;
	}
    
	return acc;
}

bool Chades_hbt_acceptance_smear::OneParticleAcceptance(int pid, Chades_hbt_part *part, double &efficiency){
	//// === HADES acceptance === ////
	efficiency=0.0;
	FourVector Plab;
	Misc::Boost(ucm,part->psmear,Plab);
  
	//momentum:
	int charge = pid / abs(pid);
	double pmag = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]+Plab[3]*Plab[3]);
	if(charge > 0 && pmag < 50.0)
		return false;
	if(charge < 0 && pmag < 105.0)
		return false;

	//pT:
	double pT = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]);
	if(pT<pTmin && pT>pTmax)
		return false;

	//theta:
	double theta = atan(pT/Plab[3])*(180.0/PI);
	if(theta < thetamin || theta > thetamax)
		return false;

	//rapidity:
	double rapidity = 0.5 * log((Plab[0]+Plab[3])/(Plab[0]-Plab[3]));
	if(rapidity > ymax || rapidity < ymin)
		return false;

	//// === Efficiency === ////
	efficiency=1.0;
	
	return true;
}


void Chades_hbt_acceptance_smear::Smear(Chades_hbt_part *part){

	FourVector plab;
	part->Setp0();
	Misc::Boost(ucm,part->p,plab);

	int const nrBins  = 46;
	int const nrBinsZ = 23;
	double momRanges[nrBins+1] ={-1200, -700, -650, -600, -550, -500, -450, -400, -375, -350, -325, -300, -275, -250, -225, -200, -175, -150, -125, -100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 1200};
	double momRangesZ[nrBinsZ+1] ={0, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 1200};
        
        
	int momRangeX, momRangeY, momRangeZ = 0;
	for(int iMomX = 0; iMomX < nrBins-1; iMomX++){
		if(plab[1] >= momRanges[iMomX] && plab[1] < momRanges[iMomX+1]){
			momRangeX = iMomX; break;
		}
		// else if (iMomX == nrBins-1 && plab[1] < momRanges[0]) {momRangeX = 0; break;}
		else if (iMomX == nrBins-1) {momRangeX = nrBins; break;}
		momRangeX = nrBins;
	}
	for(int iMomY = 0; iMomY < nrBins-1; iMomY++){
		if(plab[2] >= momRanges[iMomY] && plab[2] < momRanges[iMomY+1]){
			momRangeY = iMomY; break;
		}
		// else if (iMomY == nrBins-1 && plab[2] < momRanges[0]) {momRangeY = 0; break;}
		else if (iMomY == nrBins-1) {momRangeY = nrBins; break;}
		momRangeY = nrBins;
	}
	for(int iMomZ = 0; iMomZ < nrBinsZ-1; iMomZ++){
		if(plab[3] >= momRangesZ[iMomZ] && plab[3] < momRangesZ[iMomZ+1]){
			momRangeZ = iMomZ; break;
		}
		else if (iMomZ == nrBinsZ-1 && plab[3] < momRangesZ[0]) {momRangeZ = 0; break;}
		else if (iMomZ == nrBinsZ-1) {momRangeZ = nrBinsZ; break;}
		momRangeZ = nrBinsZ;
	}
	cout <<  plab[1] << "  ";
	if(part->pid == 2212){ // PROTON //
		// fit parameters //
		double GaussMeanX[nrBins+1] ={-5.58137,-3.66741,-3.43687,-3.35549,-3.37302,-3.57569,-3.9553,-4.35422,-4.67892,-5.06295,-5.44771,-5.71363,-5.80405,-5.55286,-5.02242,-4.23637,-3.30623,-2.49273,-1.91235,-1.55749,-1.19324,-0.72658,-0.250611,0.266152,0.753763,1.20581,1.59098,1.94292,2.51783,3.28986,4.21042,5.01521,5.60318,5.80978,5.76176,5.45222,5.08275,4.70293,4.36503,3.95268,3.53823,3.32256,3.29365,3.34782,3.55956,5.67388};
		double GaussSigmaX[nrBins+1] ={35.5716,23.34,20.5012,17.9763,15.7492,13.8895,12.4151,11.507,11.04,10.6364,10.3451,10.1051,9.81521,9.39345,8.86145,8.15992,7.28527,6.4546,5.76289,5.24318,4.84999,4.50923,4.3173,4.31246,4.50404,4.85214,5.25137,5.76062,6.44251,7.26837,8.14128,8.87766,9.46295,9.8683,10.1881,10.4038,10.6812,11.0611,11.5383,12.4618,13.9592,15.8313,18.0366,20.6001,23.4831,37.1498};
		double GaussMeanY[nrBins+1] ={-6.62789,-4.56648,-4.30764,-4.18323,-4.09656,-4.16742,-4.38092,-4.58909,-4.74152,-4.88759,-4.95917,-4.84416,-4.54032,-4.11168,-3.71471,-3.28698,-2.78907,-2.35197,-2.02235,-1.81035,-1.63472,-1.18765,1.55367,-1.38068,0.97129,1.44627,1.68396,1.97301,2.33265,2.80078,3.26406,3.68325,4.04015,4.45464,4.7943,4.92248,4.82405,4.67414,4.46591,4.22036,3.89998,3.71735,3.65721,3.71646,3.85323,5.77251};
		double GaussSigmaY[nrBins+1] ={34.9362,22.8639,20.2571,18.0363,16.0448,14.4046,13.076,12.2779,11.8211,11.4367,11.0557,10.5545,9.86689,9.06275,8.31344,7.56973,6.80383,6.10664,5.60062,5.31251,5.39539,6.03971,7.41911,6.80347,5.64116,5.16411,5.2123,5.56418,6.11285,6.81507,7.55823,8.31222,9.04,9.85403,10.5882,11.0937,11.5083,11.9204,12.3848,13.1663,14.4984,16.1334,18.1012,20.3409,22.8523,36.1061};
        
		double GaussMeanZ[nrBinsZ+1] ={-1.72437,3.41971,5.92205,6.08429,6.1867,6.8371,7.5259,7.2993,6.87667,6.81288,7.93318,9.50573,10.4636,10.9037,10.8191,10.2872,9.41909,8.4484,7.64279,6.97906,6.49021,6.09366,5.94862};
		double GaussSigmaZ[nrBinsZ+1] ={6.62536,8.63675,8.82611,8.69724,8.87599,9.51675,10.2278,10.3888,10.5968,11.145,12.4995,13.8957,14.6115,14.9409,14.795,14.4281,14.0735,14.0942,14.5052,15.1473,15.9322,16.8431,34.1973};
		normal_distribution<> distX(GaussMeanX[momRangeX], GaussSigmaX[momRangeX]);
		normal_distribution<> distY(GaussMeanY[momRangeY], GaussSigmaY[momRangeY]);
		normal_distribution<> distZ(GaussMeanZ[momRangeZ], GaussSigmaZ[momRangeZ]);
		default_random_engine generator;
		plab[1]=plab[1]+distX(generator);
		plab[2]=plab[2]+distY(generator);
		plab[3]=plab[3]+distZ(generator);

	} //end of proton
	cout <<  plab[1] << endl;
	if(part->pid == 22122112){ // DEUTERON //
		// fit parameters //
		double GaussMeanX[nrBins+1] ={-14.0718,-9.40807,-9.22066,-9.5316,-10.0634,-10.9013,-11.5509,-11.2561,-10.8413,-10.2003,-9.42484,-8.31359,-7.25679,-6.27827,-5.45726,-4.74262,-3.97032,-3.32045,-2.81058,-2.21519,-1.62055,-0.985304,-0.328576,0.36098,1.00696,1.64805,2.25278,2.8284,3.36377,3.98175,4.66423,5.51927,6.25332,7.25867,8.42195,9.37385,10.2855,11.0323,11.4649,11.6161,11.0543,10.1651,9.37574,9.17779,9.23119,16.4325};
		double GaussSigmaX[nrBins+1] ={37.3403,25.4553,23.6337,22.4174,21.1406,20.2204,19.1998,18.0694,17.0987,16.159,15.1907,14.0028,12.8335,11.7182,10.7632,9.9252,8.93986,8.03828,7.28698,6.54784,5.92624,5.45565,5.20211,5.21846,5.45459,5.92708,6.55124,7.25073,8.01678,8.85382,9.86781,10.807,11.6837,12.8084,14.1341,15.1529,16.2273,17.2235,18.2293,19.352,20.4298,21.3066,22.5223,23.7995,25.5813,42.5095};
		double GaussMeanY[nrBins+1] ={-15.8775,-11.2899,-10.8964,-10.8423,-10.8005,-10.8261,-10.1052,-9.06718,-8.46968,-7.90946,-7.31609,-6.71178,-6.09908,-5.39557,-4.84401,-4.20708,-3.77006,-3.39867,-3.1215,-2.98456,-2.69686,-1.39475,15.4396,-5.65873,0.849569,2.19858,2.5375,2.83881,3.21753,3.68487,4.19949,4.81395,5.46293,6.11075,6.67691,7.20333,7.77589,8.29615,8.96711,9.8964,10.6419,10.6292,10.5682,10.4346,10.637,17.0083};
		double GaussSigmaY[nrBins+1] ={37.158,26.7246,25.0665,23.5833,22.4106,21.2498,19.6,17.9495,16.6653,15.6157,14.4939,13.395,12.325,11.2081,10.2195,9.24927,8.53906,7.92703,7.59662,7.63873,8.13047,9.26068,97.9599,9.60283,8.43902,7.50488,7.16169,7.29393,7.71599,8.42244,9.21733,10.1323,11.2212,12.3087,13.3976,14.4316,15.5195,16.5419,17.8578,19.6092,21.3959,22.5296,23.8398,25.2925,26.8601,42.4746};
        
		double GaussMeanZ[nrBinsZ+1] ={-56.1186,-5.05178,6.90977,10.492,10.987,11.3923,12.3028,12.9635,13.1746,15.3676,17.0261,16.0448,14.808,14.371,14.4388,14.1978,15.0696,18.4155,20.1961,20.7071,20.0648,18.8193,19.294};
		double GaussSigmaZ[nrBinsZ+1] ={73.1058,13.4724,12.3167,12.9678,13.506,14.3122,15.447,16.3879,17.2155,19.1484,20.6899,20.7571,20.8496,21.673,22.7067,23.4047,25.0391,27.7783,29.0681,29.7199,29.5497,29.346,50.9619};
		normal_distribution<> distX(GaussMeanX[momRangeX], GaussSigmaX[momRangeX]);
		normal_distribution<> distY(GaussMeanY[momRangeY], GaussSigmaY[momRangeY]);
		normal_distribution<> distZ(GaussMeanZ[momRangeZ], GaussSigmaZ[momRangeZ]);
		default_random_engine generator;
		plab[1]=plab[1]+distX(generator);
		plab[2]=plab[2]+distY(generator);
		plab[3]=plab[3]+distZ(generator);
        
	}//end of deuteron
    
    
	if(part->pid == 50000) // 3He //
	{
		// fit parameters //
		double GaussMeanX[nrBins+1] ={23.3296,22.4,22.6593,21.252,19.8727,17.9199,15.9621,13.8854,12.1824,9.89427,9.07581,8.73657,6.95355,6.49875,6.17407,5.96857,4.81737,4.44065,3.43794,2.79732,1.77529,0.91033,0.244348,-0.424017,-1.23185,-1.95401,-2.77078,-3.58155,-4.65188,-4.99336,-5.38849,-6.11713,-6.7413,-7.40857,-7.46422,-9.24435,-11.5005,-11.9724,-13.2146,-14.7874,-17.7343,-19.7434,-21.5659,-22.1246,-23.1729,-23.6064};
		double GaussSigmaX[nrBins+1] ={31.9781,24.0158,22.9301,22.1416,21.0567,19.4475,18.1967,17.1121,16.1931,14.6449,14.137,13.5271,12.6073,12.0116,12.0505,11.5443,10.7174,10.4623,10.2433,10.0256,9.51295,9.58531,9.44486,9.49105,9.37789,9.81726,10.1131,10.2497,10.8891,11.242,11.1737,11.9779,12.2968,12.4659,12.5769,13.6582,15.969,16.3279,17.6659,18.191,20.3178,20.7817,22.2142,22.7386,24.3755,37.5942};
		double GaussMeanY[nrBins+1] ={24.4264,19.7617,18.542,17.0257,16.1388,14.4979,12.5565,11.0382,10.1297,9.44867,9.42353,8.67237,8.38649,7.20523,6.69474,6.6783,6.39812,5.41189,5.11699,0.0494039,-24.2319,-83.1667,64.4755,90,121.93,6.15868,-1.57629,-3.72272,-5.49585,-5.60813,-5.56123,-6.40782,-7.36308,-8.07872,-7.92471,-9.07267,-9.39266,-10.4065,-11.0338,-12.0231,-14.1562,-16.3187,-17.2099,-18.4523,-19.4598,-24.0843};
		double GaussSigmaY[nrBins+1] ={32.6177,25.5474,23.82,21.6468,20.5493,18.6861,17.0679,15.9983,15.0017,14.6086,14.2801,13.6563,12.7894,12.5973,12.1669,11.9288,11.9382,11.5554,11.4488,11.5871,149.343,2,1005.47,43.9026,540.409,17.4319,11.0552,11.4194,10.904,10.8519,11.2707,11.5948,12.7201,12.8815,13.3689,13.8363,14.5932,15.3355,15.7542,17.2928,19.2653,20.7281,22.1669,23.7516,24.4304,37.79};
            
		int const nrBinsZ_3He =16;
		double momRangesZ_3He[nrBinsZ_3He+1] ={0, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 1200};
		for(int iMomZ = 0; iMomZ < nrBinsZ_3He-1; iMomZ++){
			if(plab[3] >= momRangesZ_3He[iMomZ] && plab[3] < momRangesZ_3He[iMomZ+1]){
				momRangeZ = iMomZ; break;
			}
			else if (iMomZ == nrBinsZ_3He-1 && plab[3] < momRangesZ_3He[0]) {momRangeZ = 0; break;}
			else if (iMomZ == nrBinsZ_3He-1) {momRangeZ = nrBins; break;}
			momRangeZ = nrBinsZ_3He;
		}
            
            
		double GaussMeanZ[nrBinsZ+1] ={-8.06199,-12.6399,-13.1529,-13.495,-13.7742,-15.8369,-16.1523,-16.8729,-17.5998,-19.3194,-21.4478,-23.9666,-26.957,-26.5938,-27.8527,-42.5801};
		double GaussSigmaZ[nrBinsZ+1] ={15.1634,15.0585,15.2009,16.5073,16.0282,17.1987,17.6616,19.3194,19.4671,19.9868,21.4826,22.8324,24.1461,25.6958,27.6953,56.295};
            
		normal_distribution<> distX(GaussMeanX[momRangeX], GaussSigmaX[momRangeX]);
		normal_distribution<> distY(GaussMeanY[momRangeY], GaussSigmaY[momRangeY]);
		normal_distribution<> distZ(GaussMeanZ[momRangeZ], GaussSigmaZ[momRangeZ]);
		default_random_engine generator;
		plab[1]=plab[1]+distX(generator);
		plab[2]=plab[2]+distY(generator);
		plab[3]=plab[3]+distZ(generator);
            
	}
    
	plab[0]=sqrt(part->mass*part->mass+plab[1]*plab[1]+plab[2]*plab[2]+plab[3]*plab[3]);
	part->Setp0();
	Misc::BoostToCM(ucm,plab,part->psmear);
   
}
/*

void Chades_hbt_acceptance_smear::Smear(Chades_hbt_part *part){
double sigmap1,sigmap2,sigmap3;
FourVector plab;
part->Setp0();
Misc::Boost(ucm,part->p,plab);
    
sigmap1=10.0+fabs(plab[1]*0.02);
sigmap2=10.0+fabs(plab[2]*0.02);
sigmap3=10.0+fabs(plab[3]*0.02);
plab[1]=plab[1]+sigmap1*randy->ran_gauss();
plab[2]=plab[2]+sigmap2*randy->ran_gauss();
plab[3]=plab[3]+sigmap3*randy->ran_gauss();
plab[0]=sqrt(part->mass*part->mass+plab[1]*plab[1]+plab[2]*plab[2]+plab[3]*plab[3]);
    
Misc::BoostToCM(ucm,plab,part->psmear);
part->Setp0();
}

*/
bool Chades_hbt_acceptance_smear::TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
double &efficiency){
    
	//// === HADES acceptance === ////
	bool acc=true;
    
	efficiency=1.0;
   
	FourVector PlabA;
	parta->Setp0();
	Misc::Boost(ucm,parta->psmear,PlabA);
	FourVector PlabB;
	partb->Setp0();
	Misc::Boost(ucm,partb->psmear,PlabB);
    
	double ptA = sqrt(PlabA[1]*PlabA[1] + PlabA[2]*PlabA[2]);
	double ptB = sqrt(PlabB[1]*PlabB[1] + PlabB[2]*PlabB[2]);
	double thetaA = atan(ptA/PlabA[3]);
	double thetaB = atan(ptB/PlabB[3]);
    
	if(fabs(thetaA-thetaB)<0.04 && fabs(delphi)<0.25){ // merging suppression
		//efficiency=0.0;
		acc=false;
	}
    
	// kT bin //
	double kT_min = 350;
	double kT_max = 500;
	double kT_pair =  0.5* sqrt((PlabA[1]+PlabB[1])* (PlabA[1]+PlabB[1])+ (PlabA[2]+PlabB[2])*(PlabA[2]+PlabB[2]));
	//cout << "kT_pair " << kT_pair << endl;
	if(kT_pair <=  kT_min ||kT_pair > kT_max){ // merging suppression
		//cout << "here" << endl;
		// efficiency=0.0;
		acc=false;
	}
	return acc;
}

