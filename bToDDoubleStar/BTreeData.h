#ifndef BTREE_DATA_H
#define BTREE_DATA_H
//for the Short_t datatype
#include "TTree.h"

class BTreeData
{
 public:


  int mcBCharge;


  int foundDPiPi;
  int recBToDlNuPiPi;
  //the two pion channel
  int recDDoubleStar;
  //only one pion
  int recBToDlNuPi;
  float logProb;

  int bestBCharge;


  float deltaETag;
  float mBTag;
  int tagId;
  int tagDecay;
  float tagCorr;
  //is a fixed values (4.5%/4.3%)
  float sTagCorr;

  float B_DecayCorr;
  float sB_DecayCorr;

  float pidCorrection;
  float sPidCorrection;
  float CrossSectionLumiCorrection;
  float sCrossSectionLumiCorrection;
  float CrossSectionCorrection;
  float sCrossSectionCorrection;

  float lumiCorrection;
  float sLumiCorrection;

  float FFDCorrection[200];
  float sFFDCorrection[200];

  float FFDdsCorrection[200];
  float sFFDdsCorrection[200];



  float sRelTrackingCorrection;

  int numTracks;

  float pi0Correction[200];
  float sPi0Correction[200];
  Short_t pi0MomBin[200];
  int numPi0Corr;

  float KsCorrection[200];
  float sKsCorrection[200];
  Short_t KsCorrMomBin[200];
  Short_t KsCorrThetaBin[200];
  int numKsCorr;

  float ChargedCorrection[200];
  float sChargedCorrection[200];
  Short_t ChargedCorrThetaBin[200];
  Short_t ChargedCorrMomBin[200];
  Short_t ChargedCorrSVDBin[200];
  Short_t ChargedCorrMisIdType[200];
  int numChargedCorr;

  Short_t D_pBinFF[200];
  Short_t D_q2BinFF[200];
  Short_t D_TypeBinFF[200];
  int numFFDCorr;

  Short_t Dds_cosTBinFF[200];
  Short_t Dds_wBinFF[200];
  Short_t Dds_TypeBinFF[200];
  int numFFDdsCorr;



  float D_DecayCorr[200];
  float sD_DecayCorr[200];
  Short_t D_decType[200];
  int numDDecCorr;

  //for the correction
  Short_t B_decType;

  Short_t experiment;


  int leptonId;

  int found_2SD;
  int found_2SD_Star;
  int mcDecaySignature;
  int recDecaySignature;

  //more generic, decay signature
  int foundAnyDDoubleStar;
  int sig_numPions;
  int sig_numLeptons;
  int sig_numKaons;
  int sig_numBaryons;
  int sig_numPi0;
  int sig_numD;
  int sig_numDStar;


  int sig_dStar_2S;
  int sig_d_2S;

  int sigResDLNu;
  int sigResDPiLNu;
  int sigResDPiPiLNu;

  int sigResDStarLNu;
  int sigResDStarPiLNu;
  int sigResDStarPiPiLNu;

  int sigDLNu;
  int sigDPiLNu;
  int sigDPiPiLNu;

  int sigDStarLNu;
  int sigDStarPiLNu;
  int sigDStarPiPiLNu;
  //

  float leptonP;
  int foundLepton;
  float piPlusP;
  int foundPiPlus;
  float piMinusP;
  int foundPiMinus;
  float dMesonP;
  int foundDMeson;
  int dPID;

  float tagOverlapFractionCharged;
  float tagOverlapFractionPi0;
  int foundDDoubleStarInMC;
  int daughterDPID;
  int overlapEvent;
  int mcDCharge;
  int mcIsDStar;


  float U[1000];
  float mNu2[1000];
  float mB[1000];
  float mXl[1000];
  float mDnPi[1000];
  //recoil for B->Dlnu
  float w[1000];

  int systemCharge[1000];
  int leptonCharge[1000];
  int dCharge[1000];
  int isDStar[1000];



  int recDType[1000];
  int numRecPions[1000];
  int bestD[1000];

  int dDecay[1000];
  int dType[1000];

  int pi1Rec[1000];
  int pi2Rec[1000];

  float leptonMom[1000];
  float leptonTheta[1000];
  float leptonPhi[1000];
  float pi1Mom[1000];
  float pi1Theta[1000];
  float pi1Phi[1000];
  float pi2Mom[1000];
  float pi2Theta[1000];
  float pi2Phi[1000];


  //there is only one real mc possibility
  float pi1Mom_mc;
  float pi1Theta_mc;
  float pi1Phi_mc;
  float pi2Mom_mc;
  float pi2Theta_mc;
  float pi2Phi_mc;
  float DDStarMass_mc;

  float pi1Found;
  float pi2Found;

  float DDiff[1000];
  float DStarDiff[1000];

  //masses of the system with only one pion combined (cut in case it is too close to D*)
  float hypDMass1[1000];
  float hypDMass2[1000];

  int size;

  //call at the beginning of the event to make sure that we don't have any stale data in the structure
  void initData()
  {
    numTracks=0;

    for(int i=0;i<200;i++)
      {

	pi0Correction[i]=0;
	sPi0Correction[i]=0;
	pi0MomBin[i]=0;
	KsCorrection[i]=0;
	sKsCorrection[i]=0;
	KsCorrMomBin[i]=0;
	KsCorrThetaBin[i]=0;
	ChargedCorrection[i]=0;
	sChargedCorrection[i]=0;
	ChargedCorrThetaBin[i]=0;
	ChargedCorrThetaBin[i]=0;
	ChargedCorrMomBin[i]=0;
	ChargedCorrSVDBin[i]=0;
	ChargedCorrMisIdType[i]=0;
	FFDCorrection[i]=0;
	sFFDCorrection[i]=0;

	FFDdsCorrection[i]=0;
	sFFDdsCorrection[i]=0;
	Dds_cosTBinFF[i]=0;
	Dds_wBinFF[i]=0;
	Dds_TypeBinFF[i]=0;
	D_pBinFF[i]=0;
	D_q2BinFF[i]=0;
	D_TypeBinFF[i]=0;
	D_DecayCorr[i]=0;
	sD_DecayCorr[i]=0;
	D_decType[i]=0;
	//	Dds_DecayCorr[i]=0;
	//	sDds_DecayCorr[i]=0;
	//	Dds_decType[i]=-1;
      }
    numPi0Corr=0;
    numKsCorr=0;
    numChargedCorr=0;

    numFFDCorr=0;
    numFFDdsCorr=0;
    numDDecCorr=0;

    B_decType=-1;
    experiment=0;

    mcBCharge=-2;
    mcIsDStar=-2;
    mcDCharge=-2;

    pi1Found=-1;
    pi2Found=-1;

    pi1Phi_mc=-1;
    pi2Phi_mc=-1;


    pi1Theta_mc=-1;
    pi2Theta_mc=-1;


    pi1Mom_mc=-1;
    pi2Mom_mc=-1;

    foundDPiPi=-1;
    recBToDlNuPiPi=-1;
  //the two pion channel
    recDDoubleStar=-1;
  //only one pion
   recBToDlNuPi=-1;
   logProb=-1;

   bestBCharge=-2;
  //systemCharge=-1;

   deltaETag=-1;
   mBTag=-1;
   tagId=-1;
   tagDecay=-1;
   tagCorr=-1;

   B_DecayCorr=1;
   sB_DecayCorr=0;
   pidCorrection=1;
   sPidCorrection=0;

   CrossSectionCorrection=1;
   sCrossSectionCorrection=0;

   CrossSectionLumiCorrection=1;
   sCrossSectionLumiCorrection=0;

   lumiCorrection=1;
   sLumiCorrection=0;



   leptonId=-1;

   found_2SD=-1;
   found_2SD_Star=-1;
   mcDecaySignature=-1;
   recDecaySignature=-1;

  //more generic, decay signature
   foundAnyDDoubleStar=-1;
   sig_numPions=-1;
   sig_numLeptons=-1;
   sig_numKaons=-1;
   sig_numBaryons=-1;
   sig_numPi0=-1;
   sig_numD=-1;
   sig_numDStar=-1;

   sig_dStar_2S=-1;
   sig_d_2S=-1;

   sigDLNu=-1;
   sigDPiLNu=-1;
   sigDPiPiLNu=-1;

   sigDStarLNu=-1;
   sigDStarPiLNu=-1;
   sigDStarPiPiLNu=-1;

   sigResDLNu=-1;
   sigResDPiLNu=-1;
   sigResDPiPiLNu=-1;

   sigResDStarLNu=-1;
   sigResDStarPiLNu=-1;
   sigResDStarPiPiLNu=-1;

   leptonP=-1;
   foundLepton=-1;
   piPlusP=-1;
   foundPiPlus=-1;
   piMinusP=-1;
   foundPiMinus=-1;
   dMesonP=-1;
   foundDMeson=-1;
   dPID=-1;

   tagOverlapFractionCharged=-2;
   tagOverlapFractionPi0=-1;
   foundDDoubleStarInMC=-1;
   daughterDPID=-1;

   size=0;

  }


};

#endif

