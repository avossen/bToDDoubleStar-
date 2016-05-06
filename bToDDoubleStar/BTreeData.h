#ifndef BTREE_DATA_H
#define BTREE_DATA_H

class BTreeData
{
 public:

  int foundDPiPi;
  int recBToDlNuPiPi;
  //the two pion channel
  int recDDoubleStar;
  //only one pion
  int recBToDlNuPi;
  float logProb;

  int bestBCharge;
  int systemCharge;

  float deltaETag;
  float mBTag;
  int tagId;
  int tagDecay;
  float tagCorr;

  float D_DecayCorr;
  float B_DecayCorr;

  float pidCorrection;
  float CrossSectionLumiCorrection;

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

  float mNu2[1000];
  float mB[1000];
  float mXl[1000];
  float mDnPi[1000];
  int recDType[1000];
  int numRecPions[1000];
  int bestD[1000];
  float leptonMom[1000];
  float leptonTheta[1000];
  float leptonPhi[1000];
  float pi1Mom[1000];
  float pi1Theta[1000];
  float pi1Phi[1000];
  float pi2Mom[1000];
  float pi2Theta[1000];
  float pi2Phi[1000];

  int size;

  //call at the beginning of the event to make sure that we don't have any stale data in the structure
  void initData()
  {

    foundDPiPi=-1;
    recBToDlNuPiPi=-1;
  //the two pion channel
    recDDoubleStar=-1;
  //only one pion
   recBToDlNuPi=-1;
   logProb=-1;

  bestBCharge=-1;
  systemCharge=-1;

   deltaETag=-1;
   mBTag=-1;
 tagId=-1;
  tagDecay=-1;
   tagCorr=-1;

   D_DecayCorr=-1;
   B_DecayCorr=-1;
   pidCorrection=-1;
   CrossSectionLumiCorrection=-1;

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

   leptonP=-1;
   foundLepton=-1;
   piPlusP=-1;
   foundPiPlus=-1;
   piMinusP=-1;
   foundPiMinus=-1;
   dMesonP=-1;
   foundDMeson=-1;
   dPID=-1;

   tagOverlapFractionCharged=-1;
   tagOverlapFractionPi0=-1;
   foundDDoubleStarInMC=-1;
   daughterDPID=-1;

  size=0;

  }


};

#endif
