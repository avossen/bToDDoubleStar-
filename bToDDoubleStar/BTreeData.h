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

  float deltaETag;
  float mBTag;
  int tagId;
  int tagDecay;
  float tagCorr;


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



};

#endif
