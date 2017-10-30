#ifndef SYS_DATA_GEN
#define SYS_DATA_GEN

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "Fit/Fitter.h"
#include <math.h>
#include "TRandom3.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TColor.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include "TFractionFitter.h"
#include <sstream>

using namespace std;

#include "idx.h"

class sysDataGen
{

 public: 
  //reads the data from the trees
  sysDataGen(TTree** trees);
  ~sysDataGen(){};

  void readTrees();
  //get the templates for a specific channel. At this point the channel selection cuts are made
  TH1F** getTemplates(int channelIdx, int leptonId, char* channelString, TH1F** components, bool doSysStudies=false, int sysIndex=-1);


  //store data in mNu2,weight and weight uncertainty tuples
  //(here weight uncertainty is a combination of all uncertainties)
  //we store the data separately for each template
  //then each iteration will only change the weights in this template
  struct dataElement
  {
    short tagId;
    float mNu2;
    float weight;
    float eWeight; 

    //the elements needed to do the channel selection
    int mcBCharge;
    int bestBCharge;
    int dCharge;
    int mcDCharge;
    int dType;
    int mcIsDStar;
    int leptonId;
    short iTree;
    float tagCorr;

    /////---->weight details...:
    float B_DecayCorr;
    float sB_DecayCorr;
    float CrossSectionLumiCorrection;
    float sCrossSectionLumiCorrection;
    float sRelTrackingCorrection;
    int numTracks;
    int B_decayType;
    int experiment;
    int numKsCorrection;
    float KsCorrection[20];
    float sKsCorrection[20];
    int numChargedCorrection;
    float ChargedCorrection[50];
    float sChargedCorrection[50];
    int numPi0Correction;
    float pi0Correction[100];
    float sPi0Correction[100];
    int numD_DecayCorr;
    float D_decayCorr[10];
    float sD_decayCorr[10];
    int numFFDCorrection;
    float ffDCorrection[4];
    float sFfDCorrection[4];
    int numFFDdsCorrection;
    float ffDdsCorrection[4];
    float sFfDdsCorrection[4];
    short ChargedCorrCombBin[50];
    short KsCorrCombBin[20];
    short pi0CorrCombBin[100];
    short DFFCorrCombBin[4];
    short DdsFFCorrCombBin[4];
    short D_decTypeCorrBin[10];

    //    int recDType;


    ///////--->

    //so we can build the components per input file
    int fileIndex;

  };

  dataElement* pContinuum;
  dataElement* pDDStar;
  dataElement* pDDStarPi;
  dataElement* pDDStarPiWrongChannel;
  dataElement* pDDStarPiPi;
  dataElement* pDLNu;
  dataElement* pDPiPiLNu;
  dataElement* pDStarLNu;
  dataElement* pDStarPiPiLNu;
  dataElement* pDDStarPiCrossFeed;
  dataElement* pOtherBB;

  dataElement* pData[20];


  char* histoNames[20];

  TTree** mTrees;

  int currentCounts[20];


 protected:
  static const int maxKs =20;
  static const int maxCharged =50;
  static const int maxPi0 =100;
  static const int maxDDec =10;
  static const int maxFFD =4;
  static const int maxFFDds =4;

  //misIdType, SVDBin, mom, theta
  //9x2x30x12
  float ChargedWeights[6480];
  //mom, theta
  //3,4
  float KsWeights[12];
  //24 mom bins
  float pi0Weights[24];
  //24 p, 12 q2, 2 type
  float D_FFWeights[576];

  //17 w bis, 10 cosT, 4 types: 680
  float Dds_FFWeights[680];
  float D_decTypeWeights[26]; //this is for the D_decay
  //dec type wasn't initialized correctly, but we can check for sB_DecayCorr==0 to weed
  //out the cases were we haven't found anything
  float B_decWeights[12];
  //per experiment
  float lumiWeight[100];
  //id 511 and 521
  float tagIdWeight[4];
  //weight shift for each number of tracks
  //(0.35% rel uncert per track)
  float trackWeights[200];



  void fillWeightDetails(dataElement& element,int iTree);
  float getWeight(int iTree,bool foundDDoubleStar);
  float getWeightWithError(dataElement& element, bool foundDDoubleStar);
  TRandom3 rnd;
  bool isContinuum();
      //iDDStar
  bool isDDStar();
      //DDStarPi
  bool isDDStarPi();


      //DDStarPiWrongChannel
  //  bool isDDStarPiWrongChannel();
      
      //DDStarPiPi
  bool isDDStarPiPi();
      //DLNu
  bool isDLNu();
      //DPiPiLNu
  bool isDPiPiLNu();
      //DStarLNu
  bool isDStarLNu();
      //DStarPiPiLNu
  bool isDStarPiPiLNu();

      //DDStarPiCrossFeed
  //  bool isDDStarPiCrossFeed();

      //OtherBB
  bool isOtherBB();
  //branch addresses
  static const int MaxArraySize=100;
  Int_t mNu2Counter;
  Int_t tagId;
  Float_t mNu2[MaxArraySize];
  Float_t tagCorr;
  Float_t B_DecayCorr;
  Float_t sB_DecayCorr;
  //  Float_t PID_Correction;
  //  Float_t sPID_Correction;
  Float_t CrossSectionLumiCorrection;
  Float_t sCrossSectionLumiCorrection;
  Float_t FFCorrection;
  Float_t PIDCorrection;
  Float_t sFFCorrection;
  Float_t sRelTrackingCorrection;


  //the shifts are generated for the systematic studies for each fit. 
  //since the variation of the weights are correlated event-by-event, we should have one shift and not
  //generate a new weight for every event (that would most likely understimate the systematics, because the
  //effect of the weight uncertainty would cancel out over many events
  float FFCorrectionRelShift;
  float tagCorrRelShift;
  float D_DecayCorrRelShift;
  float B_DecayCorrRelShift;
  float CrossSectionLumiCorrectionRelShift;
  Short_t experiment;
  Int_t bestBCharge;
  Int_t mcBCharge;
  Int_t systemCharge;
  Int_t leptonId[MaxArraySize];
  Int_t foundAnyDDoubleStar;
  Int_t sig_numPions;
  Int_t sig_numKaons;
  Int_t sig_numPi0;
  Int_t sig_numBaryons;
  Int_t sig_numLeptons;
  Int_t sig_DLNu;
  Int_t sig_DPiLNu;
  Int_t sig_DPiPiLNu;
  Int_t sig_DStarLNu;
  Int_t sig_DStarPiLNu;
  Int_t sig_DStarPiPiLNu;
  Int_t numRecPions[MaxArraySize];
  Int_t dType[MaxArraySize];
  Int_t dDecay[MaxArraySize];
  Int_t mcIsDStar;
  Float_t mDnPi[MaxArraySize];
  Float_t pi1Mom[MaxArraySize];
  Float_t pi2Mom[MaxArraySize];
  Int_t bestD[MaxArraySize];
  Int_t dCharge[MaxArraySize];
  Int_t recDType[MaxArraySize];
  Int_t mcDCharge;
  Int_t mcDecaySignature;
  Int_t recDecaySignature;
  Float_t mBTag;
  Float_t logProb;
  Float_t deltaETag;
  ///new for corrections
  Int_t numTracks;

  ///all the weight stuff...:
  Float_t KsCorrection[MaxArraySize];
  Float_t sKsCorrection[MaxArraySize];
  Int_t KsCorrectionCounter;

  Float_t ChargedCorrection[MaxArraySize];
  Float_t sChargedCorrection[MaxArraySize];
  Int_t ChargedCorrectionCounter;


  Float_t Pi0Correction[MaxArraySize];
  Float_t sPi0Correction[MaxArraySize];
  Int_t Pi0CorrectionCounter;

  Float_t D_DecayCorr[MaxArraySize];
  Float_t sD_DecayCorr[MaxArraySize];
  Int_t D_DecayCorrCounter;

  Float_t FFDCorrection[MaxArraySize];
  Float_t sFFDCorrection[MaxArraySize];
  Int_t FFDCorrectionCounter;

  Float_t FFDdsCorrection[MaxArraySize];
  Float_t sFFDdsCorrection[MaxArraySize];
  Int_t FFDdsCorrectionCounter;

  Short_t ChargedCorrThetaBin[MaxArraySize];
  Short_t ChargedCorrMomBin[MaxArraySize];
  Short_t ChargedCorrSVDBin[MaxArraySize];
  Short_t ChargedCorrMisIdType[MaxArraySize];

  Short_t KsCorrMomBin[MaxArraySize];
  Short_t KsCorrThetaBin[MaxArraySize];

  Short_t pi0MomBin[MaxArraySize];
  Short_t D_pBinFF[MaxArraySize];
  Short_t D_q2BinFF[MaxArraySize];
  Short_t B_decayType;
  Short_t D_TypeFF[MaxArraySize];
  Short_t Dds_wBinFF[MaxArraySize];
  Short_t Dds_cosTBinFF[MaxArraySize];

  Short_t Dds_TypeFF[MaxArraySize];
  Short_t D_decType[MaxArraySize];






  float mixedFactor;
  float chargedFactor;

};


#endif
