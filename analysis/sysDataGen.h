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
  TH1D** getTemplates(int channelIdx, int leptonId, char* channelString);


  //store data in mNu2,weight and weight uncertainty tuples
  //(here weight uncertainty is a combination of all uncertainties)
  //we store the data separately for each template
  //then each iteration will only change the weights in this template
  struct dataElement
  {
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

  float getWeight();
  float getWeightUncertainty();

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
  Float_t mNu2;
  Float_t tagCorr;
  Int_t bestBCharge;
  Int_t systemCharge;
  Int_t leptonId;
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
  Int_t numRecPions;
  Int_t dType;
  Int_t dDecay;
  Int_t mcIsDStar;
  Float_t mDnPi;
  Float_t pi1Mom;
  Float_t pi2Mom;
  Int_t bestD;
  Int_t dCharge;
  Int_t recDType;
  Int_t mcDCharge;
  Int_t mcDecaySignature;
  Int_t recDecaySignature;
  Float_t mBTag;
  Float_t logProb;
  Float_t deltaETag;


};


#endif
