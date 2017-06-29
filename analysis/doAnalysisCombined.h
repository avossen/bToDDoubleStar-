#ifndef DO_ANALYSIS_H
#define DO_ANALYSIS_H

//implement the factor ~2.5/2.3 to account for the different lumis between signal and other MC

#define Huschle_Signal_Weighting

const bool combineDPiPi=true;
const int oneIdx=-1;

//#define DPiPi_SEARCH
int maxDataTreeSize;
int SIG_IDX;
#define partialBoxFraction 0.15
//#define P0STRING "( recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.0 "
//Robnin's cuts for x-check
#define P0STRING "( recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.24  && logProb > -3.0 "
//#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.0 && mDnPi < 3.5  "
#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.0 && mDnPi < 3.5  &&mDnPi > 2.05"
////


#define P2STRING_SinglePion "bestD==1 && pi2Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb  > -2.6  && mDnPi < 3.0 && !(dType==1 && dDecay==2) && !(dType==1 && dDecay==3) && !(dType==2) && !overlapEvent" 

#define P2STRING_Search "bestD==1 && pi2Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb  > -2.6  && mDnPi < 3.0 && ((bestBCharge==0 && systemCharge==0) || (bestBCharge== -leptonCharge) ) && mDnPi>1.0  && mDnPi < 3.0&&  hypDMass1 > 2.03 && hypDMass2 > 2.03 && U < 2.0 && U > -2.0"  


//3,5,7 are the D pi pi 


#define SIG_IDX_D_PI 2
#define SIG_IDX_D_PI_PI 3

#ifdef DPiPi_SEARCH
#define P2STRING P2STRING_Search
#else
#define P2STRING P2STRING_SinglePion
#endif



double gl_templateScaleFactor;

//#define P2STRING "recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb  > -3.0  && mDnPi < 3.0 "
//&& !(dType==0 && dDecay==1) && !(dType==1 && dDecay==3) && !(dType==2 && dDecay==1) 

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

char channelStringGlobal[300];
int pCount;

TH1D** sigSignificance;
int glChannelIdx;

void getChannelString(int channel, char* buffer);
void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel=false);
TRandom3* rnd;
int numBinsSB=20;
int numBinsSBTop=40;
int numBinsWS=20;
char* allLegendNames[20];

enum selectionNames{};

//int numBins=40;
//float lowerCut=-0.3;
//float upperCut=0.5;

//int numBins=150;

////int numBins=20;
//int numBins=100;


//we get meaningful results with 50 bins between -1 and 2..
///->allChannels, BToD, B0ToDStar, B0ToD, B0ToDStar
//int numBins[]={50,30,30,30,30};
//int numBins[]={70,200,70,200,70};
///////--> use this, but for partial box might need less...int numBins[]={70,70,70,70,70}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
#ifdef PARTIAL_BOX
int numBins[]={40,10,10,10,10}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma..
#else
/////int numBins[]={70,70,70,70,70}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={30,30,40,40,40}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//has to be the same for the channels we plan to combine
///-->int numBins[]={30,40,40,40,40}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
int numBins[]={20,8,8,8,8}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
#endif
//int numBins[]={200,2000,2000,200,200}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={50,50,5,200,5}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={70,100,70,70,70}; 
///_--> quite good int numBins[]={70,80,70,70,70};
///bias for BToD goes away with range -0.5 to 2.0 and 80 bins

///--->B0ToD seems to be fine with range of -0.5 to 1.5 with 70 bins...
////--->optimalfloat lowerCut[]={-0.65,-0.6,-0.65,-0.55,-0.25};
//float lowerCut[]={-0.65,-0.6,-0.65,-0.55,-0.25};
///-->should be -0.5, test -0.3 for fit stability 
//////float lowerCut[]={-0.5,-0.5,-0.5,-0.5,-0.5};
//needs to be double, otherwise there are problems with add...
Double_t lowerCut[]={-0.5,-0.3,-0.3,-0.3,-0.3};
//Double_t lowerCut[]={-0.1,-0.1,-0.1,-0.1,-0.1};
//Double_t lowerCut[]={-0.0,-0.0,-0.0,-0.0,-0.0};
////////float upperCut[]={2.0,2.0,1.5,1.5,1.5};

///for fit stability set the upper cut a bit lower for the DStar where the counts are lower
Double_t upperCut[]={0.5,0.5,0.4,0.5,0.4};
//Double_t upperCut[]={0.3,0.3,0.3,0.3,0.3};

//float lowerCut[]={-0.5,-0.5,-0.5,-0.5,-0.5};
//float upperCut[]={2.0,2.0,1.5,1.5,1.5};
//float upperCut[]={1.0,1.0,1.0,1.0,1.0};
//float upperCut[]={1.0,0.25,0.2,0.3,0.15};
//float lowerCut=-0.5;

//float lowerCut=-1.0;
//float lowerCut=-0.2;
//float lowerCut=-0.5;
////float lowerCut=-0.3;
//float upperCut=2.0;
//float upperCut=0.4;
//float upperCut=0.25;

//float upperCut=1.0;

bool withPIDCorrection;
bool withLumiCorrection;
bool withBCorrection;
bool withDCorrection;
bool withFFCorrection;
float getFitSignal( TH1F* data,TH1F** templates,int numTemplates,     TH1F* &result, TH1F** mcPredictions, double& fitVal, double& fitErr , int fixThresholdCounts, double* allFitVals,double* allFitErrs, int &numEffective, vector<int>& effectiveComponentsIndices, int& status, int numPions);
float getFitSignal_RooFit( TH1F* data,TH1F** templates,int numTemplates, TH1F* &result, TH1F** mcPredictions, double& fitVal, double& fitErr , int fixThresholdCounts, double* allFitVals,double* allFitErrs, int &numEffective, vector<int>& effectiveComponentsIndices, int& status, int numPions);

void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents, TH1D* pulls,double signalFraction, int fixThresholdCounts,int numPions);
void fillTemplates(TH1F** templates,TH1F** summedComponents,char* templateLegendNames,int numPions);
void addPoissonNoise(TH1F* h);

TColor* glColorTable[20];


//just to split some stuff off
#include "idx.h"


//insert Huschle signal weighting. This assumes that the 'foundDDoubleStarFlag' is correct
void insertSignalMCWeighting(char* buffer,  int fileIndex)
{
#ifdef Huschle_Signal_Weighting
  float mixedFactor=1.0/0.852255195;
  float chargedFactor=1.0/0.9266338302;
#else
  float mixedFactor=1.0;
  float chargedFactor=1.0;
#endif
  //this is to make it compatible with the 1+foundAnyDDoubleStar...
  mixedFactor-=1.0;
  chargedFactor-=1.0;
  char tmpBuffer[1000];
  cout <<"getting buffer: " << buffer <<endl;
  //iF==0 -->mixed, iF==1 --> charged
  //  sprintf(huschleMC_lumi_corr,"(1+foundAnyDDoubleStar*%f)*",huschleLumiFactor);
  //is this one of the DDStar channels?
  //  if(component==iDDStar || component==iDDStarPi || component==iDDStarPiPi || component==iDDStarPiWrongChannel|| component==iDDStarPiCrossFeed) 
    {
      //mixed
      if(fileIndex==0)
	{
	  sprintf(tmpBuffer,"%s",buffer);
	  sprintf(buffer,"(1+foundAnyDDoubleStar*%f)* %s",mixedFactor,tmpBuffer);
	}
      if(fileIndex==1)
	{
	  sprintf(tmpBuffer,"%s",buffer);
	  sprintf(buffer,"(1+foundAnyDDoubleStar*%f)* %s",chargedFactor,tmpBuffer);
	}

    }

}

void addCorrections(char* buffer)
{
  stringstream str;
  //D_DecayCorr*B_DecayCorr*PIDCorrection*CrossSectionLumiCorrection

  if(withFFCorrection)
    {
      str << "FFCorrection*";
    }
  if(withPIDCorrection)
    {
      str << "PIDCorrection*";
    }
  if(withDCorrection)
    {
      str <<"D_DecayCorr*";
    }
  if(withBCorrection)
    {
      str <<"B_DecayCorr*";
    }
  if(withLumiCorrection)
    {
      str <<"CrossSectionLumiCorrection*";
    }
  sprintf(buffer,"%s",str.str().c_str());
};



void getChannelString(int channel, char* buffer)
{
  if(channel<0)
    {
      sprintf(buffer,"_allChannels");
      return;
    }
  switch(channel)
    {
    case 0:
      sprintf(buffer,"_BToD_");
      break;
    case 1:
      sprintf(buffer,"_BToDStar_");
      break;
    case 2:
      sprintf(buffer,"_B0ToD_");
      break;
    case 3:
      sprintf(buffer,"_B0ToDStar_");
      break;
    case 4:
      sprintf(buffer,"_BToD_DStar");
      break;
    case 5:
      sprintf(buffer,"_B0ToD_DStar_");
      break;
    default:
      sprintf(buffer,"_wrongChannel_");
      break;
    }


}


//gets the string, similar to charge and star selection, that selects contributions from the D* channels where we loose a pi0
//for now, we don't care for the D* -> D pi+ channel because we hope that we would notice the ensuing change in charge.
void getFeedDown(char* feedDownSelection, int channel, bool dataAndMC, int numPions)
{
  //no feed down if we use all channels
  if(channel<0)
    {
      //need somestatement here since the string gets appended via &&, so for wrong channel accept none (there is no wrong channel if we take all)
      sprintf(feedDownSelection," 2!=2 ");
      return;
    }
  //  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 0==mcIsDStar");
  if(numPions==1)
    {
      switch(channel)
	{
	case 0:
	  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
	  sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 1==mcIsDStar");
	  break;
	  //B to DStar
	case 1:
	  sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && 0==mcIsDStar");
	  break;
	  //B0to D0
	case 2:
	  sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar==1");
	  break;
	  //B0to DStar0
	case 3:
	  sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && mcIsDStar==0");
	  break;
	default:
	  cout <<" wrong channel!!! (feed down) " << endl;
	  break;
	}
    }
  else
    {
      if(numPions==0)
	{
	  switch(channel)
	    {
	      //for the zero pion case, the b and D have different charges, so case 0 is B+ --> D0 lnu
	    case 0:
	      //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
	      sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && 1==mcIsDStar");
	      break;
	      //B to DStar
	    case 1:
	      sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && 0==mcIsDStar");
	  break;
	  //B0to D0
	    case 2:
	      sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && mcIsDStar==1");
	      break;
	      //B0to DStar0
	    case 3:
	      sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && mcIsDStar");
	  break;
	    default:
	      cout <<" wrong channel!!! (feed down) " << endl;
	      break;
	    }
	}
      else
	{
	  if(numPions==2)
	    {
	      switch(channel)
		{
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(feedDownSelection," 1==2 ");
		  break;
		  //B to DStar
		case 1:
		  //is D* but we only reconstruct D. This also means that the charges are different
		  sprintf(feedDownSelection,"   abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge!=mcDCharge && dType!=2 && abs(mcDCharge)==1  && 0=mcIsDStar");
		  break;
		  //B0to D0
		case 2:
		  sprintf(feedDownSelection,"  1==2");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(feedDownSelection,"   abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge!=mcDCharge && dType!=2 && abs(mcDCharge)==0  && 0=mcIsDStar");
		  break;
		default:
		  cout <<" wrong channel!!!  (feed down)" << endl;
		  break;
		}
	      
	    }
	}
    }
}

void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel)
{
  //means that we want all channels
  if(channel<0)
    {
      //need somestatement here since the string gets appended via &&, so for wrong channel accept none (there is no wrong channel if we take all)
      if(wrongChannel)
	sprintf(chargeAndStarSelection," 1!=1");
      else
	sprintf(chargeAndStarSelection," 1==1");
      return;
    }

  if(wrongChannel)
    {
      if(numPions==0)
	{
	  switch(channel)
	    {
	      //B to D
	    case 0:
	      //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar)  --> since we now have the dStar feeddown extra, we have to make sure that 
	      //we don't have overlap with that, but still get all. Should be fine to just remove the mcIsDStar check and make sure that the crossfeed only picks up events where
	      //everything else is correct
	      //question: dType selection seems stupid since it should already be in the regular channel selection? Same for bestBCharge
	      //	      sprintf(chargeAndStarSelection,"  (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge ||  dType==2 || abs(mcDCharge)!=1  ||  0!=mcIsDStar )");
	      //we don't have MC info on the D meson
	      sprintf(chargeAndStarSelection,"  (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=0)");
	      break;
	      //B to DStar
	    case 1:
	      sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=0 )");
	      break;
	      //B0to D0
	    case 2:
	      sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || abs(dCharge)!=1   )");
	      break;
	      //B0to DStar0
	    case 3:
	      sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || abs(dCharge)!=1  )");
	      break;
	    default:
	      cout <<" wrong channel!!!  " << channel <<endl;
	      break;
	    }
	}
      else
	{
	  if(numPions==1)
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar)  --> since we now have the dStar feeddown extra, we have to make sure that 
		  //we don't have overlap with that, but still get all. Should be fine to just remove the mcIsDStar check and make sure that the crossfeed only picks up events where
		  //everything else is correct
		  //question: dType selection seems stupid since it should already be in the regular channel selection? Same for bestBCharge
		  //	      sprintf(chargeAndStarSelection,"  (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge ||  dType==2 || abs(mcDCharge)!=1  ||  0!=mcIsDStar )");

		  sprintf(chargeAndStarSelection,"  (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge || abs(mcDCharge)!=1  )");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge || abs(mcDCharge)!=1  )");
		  break;
	      //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge  || abs(mcDCharge)!=0  )");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge  || abs(mcDCharge)!=0  )");
	      break;
		default:
		  cout <<" wrong channel!!!  " << channel <<endl;
		  break;
		}
	    }	  
	  else//numPions==2
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection," 1==2 ");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge==mcDCharge || abs(mcDCharge)!=1  )");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection," 1==2 ");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge==mcDCharge  ||abs(mcDCharge)!=0  )");
		  break;
		default:
		  cout <<" wrong channel!!!  " << channel <<endl;
		  break;
		}
	    }
	}
      return;
    }
  else// correct channel selection
    {

      //problem is that we do not write out the charge for the d-meson for this channel or if it is a dStar. But it should be fine just checking the charge of the B meson 
      if(numPions==0)
	{
	  if(dataAndMC)
	    {
	      switch(channel)
		{
		  //B to D0
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==0 && dType!=2");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==0 && dType==2");
		  break;
		  //B0to D
		case 2:
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && abs(dCharge)==1 && dType!=2");
		  break;
		  //B0to DStar
		case 3:
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && abs(dCharge)==1 && dType==2");
		  break;
		default:
		  cout <<" wrong channel!!!  " << endl;
		  break;
		}
	    }
	  else //data cuts only
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection,"  abs(bestBCharge)==1 && abs(dCharge)==0 && dType!=2 ");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection,"    abs(bestBCharge)==1 && abs(dCharge)==0 && dType==2");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection,"   abs(bestBCharge)==0 && abs(dCharge)==1 && dType!=2");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"    abs(bestBCharge)==0 && abs(dCharge)==1 && dType==2");
		  break;
		default:
		  cout <<" wrong channel!!!  " << endl;
		  break;
		}
	    }
	}
      if(numPions==1)
	{
	  if(dataAndMC)
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 0==mcIsDStar");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && mcIsDStar");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar==0");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && mcIsDStar");
		  break;
		default:
		  cout <<" wrong channel!!!  " << endl;
		  break;
		}
	    }
	  else //data cuts only
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection,"  abs(bestBCharge)==1 && abs(dCharge)==1 && dType!=2 ");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection,"    abs(bestBCharge)==1 && abs(dCharge)==1 && dType==2");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection,"   abs(bestBCharge)==0 && dCharge==0 && dType!=2");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"    abs(bestBCharge)==0 && dCharge==0 && dType==2");
		  break;
		default:
		  cout <<" wrong channel!!!  " << endl;
		  break;
		}
	    }
	}
      if(numPions==2)
	{
	  if(dataAndMC)
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection," 1==2 ");
		  break;
		  //B to DStar
		case 1:
		  //is D* but we only reconstruct D. This also means that the charges are different
		  sprintf(chargeAndStarSelection,"   abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge!=mcDCharge && dType!=2 && abs(mcDCharge)==1  && mcIsDStar");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection,"  1==2");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"   abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge!=mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar");
		  break;
		default:
		  cout <<" wrong channel!!!  " << endl;
		  break;
		}
	    }
	  else //data cuts only
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection,"  1==2");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection,"    abs(bestBCharge)==1 && abs(dCharge)==0 && dType!=2");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection,"   1==2");
		  break;
		  //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"    abs(bestBCharge)==0 && abs(dCharge)==1 && dType!=2");
		  break;
		default:
		  cout <<" wrong channel!!!  " << endl;
		  break;
		}
	    }
	}
    }
  
}


/*
get the histogram we want to fit from the tree (given numPions and leptonId), fit with the fractions from 'summedComponents'

 */
void fitFractions(TH1F* data, TTree** trees, TH1F** summedComponents, int numComponents,int numPions,int leptonId, int channel, bool dataTree, bool addNoise, TH1D* pulls)
{
  char channelString[500];
  getChannelString(channel,channelString);
  cout <<"fitting amplitude for channel " << channelString <<endl;
  //since we name one of the channels -1
  int channelIdx=channel+1;
  //these scale
   cout <<"using scale factor " << gl_templateScaleFactor <<endl;
   for(int i=0;i<numComponents;i++)
     {
       cout <<"running over comp " << i <<" integral: " << summedComponents[i]->Integral()<<endl;
       summedComponents[i]->Sumw2();
       (*summedComponents)->SetFillStyle(1001);
       cout <<"after scale "  << summedComponents[i]->Integral()<<endl;
     }

   cout <<"using scale factor: "<< gl_templateScaleFactor <<endl;
  TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
  cout <<"fit fractions with " << numComponents <<" components"<<endl;


  char histoName[2009];
  char drawCommand[2000];
  char buffer[2000];
  char corrBuffer[2000];
  int minCounts=0;
  ///// Used for Analysis note:  int fixThresholdCounts=300;
#ifdef PARTIAL_BOX
  //probably doesn't make sense to have different thresholds for partial_box because these are the counts of the templates
  int fixThresholdCounts=100;
#else
  //  int fixThresholdCounts=300;
  int fixThresholdCounts=3000;
  if(numPions==0)
    fixThresholdCounts=400;

#endif
  //  if(leptonId>0)
  //    fixThresholdCounts=1000;
  char* templateLegendNames[50];
  //pointer given as reference and then allocated in 'getTemplates'
  TH1F** templates;
  //given to 'getTemplates' as reference
  int numMergers=0;

  //test here if the number of counts before and after is teh same
  float templateIntegral=0;
  int numTemplateEntries=0;
  for(int i=0;i<numComponents;i++)
    {
      templateIntegral+=summedComponents[i]->Integral();
      numTemplateEntries+=summedComponents[i]->GetEntries();
    }
  cout <<"template integral before summing combining: "<< templateIntegral<<", entries: "<< numTemplateEntries<<endl;
  getTemplates(summedComponents, templates, templateLegendNames, allLegendNames, numComponents, numMergers,numPions, oneIdx,combineDPiPi);
  templateIntegral=0;
  numTemplateEntries=0;
  for(int i=0;i<numComponents-numMergers;i++)
    {
      templateIntegral+=templates[i]->Integral();
      numTemplateEntries+=summedComponents[i]->GetEntries();
    }
  cout <<"and after combining: "<< templateIntegral<<", entries: "<< numTemplateEntries<< endl;
  cout <<"did set up " << numComponents-numMergers << " templates" <<endl;

  //to save counts so we can fix the components which have too little counts
  vector<int> countsOfComponents;
  vector<int> countsOfComponents2;
  vector<int>  indexOfEffComp;
  ////from the example on the root web pages..

  //SIG_IDX is the index of the template that gives the signal we are after
  cout <<" done making data ..  sig ids: "<< SIG_IDX <<endl;
  cout <<"data integral: " << data->Integral() <<endl;
  cout <<"done with 2nd integral"<<endl;

  //this is a rough initial guess
  ////---
  double signalFraction=gl_templateScaleFactor*templates[SIG_IDX]->Integral()/data->Integral();
  cout <<"for BR, signal for channel " << channel<< "  in MC is : " << templates[SIG_IDX]->Integral();
  double mcSignalIntegral=templates[SIG_IDX]->Integral();
  cout <<"signal Fraction estimated to be : " << signalFraction<<endl;
  cout <<"add noise? " << addNoise <<endl;
  getChannelString(channel,channelStringGlobal);
  glChannelIdx=channelIdx;
  
  double fitVal, fitErr;
  cout<<"trying usual fit function: " << endl;
  double* allFitVals=new double[numComponents-numMergers];
  double* allFitErrs=new double[numComponents-numMergers];
  int numEffective=0;
  //shouldn't this be the main fit? Or is there another fit by hand later on
  vector<int> _effectiveComponentsIndices;
  int _status=0;
  cout <<"data integral first: "<< data->Integral() <<" data entries: "<< data->GetEntries()<<endl;
  TH1F* result;
  TH1F* mcPredictions[100];
  //    TH1F* result = (TH1F*) _fit->GetPlot();
    double S=getFitSignal_RooFit(data,templates,numComponents-numMergers,result,mcPredictions,fitVal, fitErr, fixThresholdCounts,allFitVals, allFitErrs, numEffective,_effectiveComponentsIndices,_status, numPions);  
    //double S=getFitSignal(data,templates,numComponents-numMergers,result,mcPredictions,fitVal, fitErr, fixThresholdCounts,allFitVals, allFitErrs, numEffective,_effectiveComponentsIndices,_status, numPions);  



  //  float signalFraction=templates[SIG_IDX]/data->Integral();

  cout <<"got " << fitVal*data->Integral() << " signal counts, fraction: " << fitVal << " +- " << fitErr  << endl;
  /////-----------------
  //void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents)
  if(addNoise)
    {
      //scale data by 0.2 since otherwise it is the sum of the 5 streams...
      //	data->Scale(0.2);
      //-->do this in the 'performFractionFit... function'
      for(int nIt=0;nIt<400;nIt++)
	{
	  performFractionFitStabilityTest(templates,data,numComponents-numMergers,pulls,signalFraction, fixThresholdCounts, numPions);
	}
    }
  /////////-----------

  TCanvas sampleData;
  data->Draw();
  sampleData.SaveAs("sampleData.png");
  sampleData.SaveAs("sampleData.pdf");
  sampleData.SaveAs("sampleData.eps");
  
  int numEffectiveComponents=0;

  if (_status == 0) {                       // check on fit status
    TCanvas c;
    double templatePredIntegral=result->Integral();
    data->Draw("Ep");
    result->Draw("same");
    sprintf(buffer,"fracFit_numPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelString);
    c.SaveAs(buffer);
    sprintf(buffer,"fracFit_numPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelString);
    c.SaveAs(buffer);
    sprintf(buffer,"fracFit_numPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelString);
    c.SaveAs(buffer);
    //and do this for all parameters:
    sprintf(buffer,"fracFitComp_numPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelString);
    THStack* predComponents=new THStack(buffer,buffer);
    //    cout <<"fraction fitter has " << _fit->GetFitter()->GetNumberFreeParameters() << " free and " << _fit->GetFitter()->GetNumberTotalParameters() <<" overall parameters" <<endl;
    //we want to flip the signal (index =2 ) so that it is later
    int signalIdx=-1;
    
    //has to sum to 1.0
    double totalFraction=0.0;
    double integralResult=result->Integral();
    cout <<"fit integral: "<< integralResult <<" data integral: " << data->Integral() <<endl;
    //double integralRe2=data->Integral();

    double sumOfCompInts=0.0;

    for(int i=0;i<numEffective;i++)
      {
	if(_effectiveComponentsIndices[i]!=SIG_IDX)
	  {
	    // TH1F* mcComp=(TH1F*) _fit->GetMCPrediction(i);

	    TH1F* mcComp=mcPredictions[i];

	    double mcPredInt=mcComp->Integral();
	    sumOfCompInts+=mcPredInt;
	    if(allFitVals[i]>0)
	      {

		/////
		//by using the mcPredInt, which, even for fixed ratios contains apparently poission fluctuations, we are not guaranteed to get fixed scalefactors for fixed components back
		////
		cout <<"xa looking at effective component " << i <<" integralResult: " << integralResult <<" mcPredInt: "<< mcPredInt <<" firstr ratio: "<< integralResult/mcPredInt <<" fit val: "<< allFitVals[i] <<endl;
		float scaleFact=integralResult/mcPredInt*allFitVals[i];
		//		scaleFact/=sumOfFractions;
		cout <<"scaling mccomp2 by : "<<scaleFact<<endl;
		totalFraction+=allFitVals[i];
		mcComp->Scale(scaleFact);
	      }
	    mcComp->SetFillStyle(1001);
	    mcComp->SetFillColor(glColorTable[_effectiveComponentsIndices[i]]->GetNumber());
	    predComponents->Add(mcComp);
	    cout <<"add as : " << templateLegendNames[_effectiveComponentsIndices[i]]<<endl;
	    legend->AddEntry(mcComp,templateLegendNames[_effectiveComponentsIndices[i]],"f" );
	  }
	else
	  {
	    signalIdx=i;
	  }
      }
    //add the signal last...
    if(signalIdx>=0)
      {
	cout <<"dealing with the signal " << endl;
	//	TH1F* mcComp=(TH1F*) _fit->GetMCPrediction(signalIdx);
	TH1F* mcComp=mcPredictions[signalIdx];
	double mcPredInt=mcComp->Integral();
	sumOfCompInts+=mcPredInt;
	float scaleFact=integralResult/mcPredInt*allFitVals[signalIdx];
	//just fill with the rest
	//	float scaleFact=integralResult/mcPredInt*(1-totalFraction);
	//	scaleFact/=sumOfFractions;


	/////
	//by using the mcPredInt, which, even for fixed ratios contains apparently poission fluctuations, we are not guaranteed to get fixed scalefactors for fixed components back
	////
	cout <<"scaling mccomp by : "<<scaleFact<<endl;
	mcComp->Scale(scaleFact);
	mcComp->SetFillStyle(1001);
	mcComp->SetFillColor(glColorTable[SIG_IDX]->GetNumber());

	cout <<"BR ratio to MC for channel " << channelString << " is : "<< mcComp->Integral()/mcSignalIntegral<<endl;
	predComponents->Add(mcComp);
	//	legend->AddEntry(mcComp,templateLegendNames[SIG_IDX],"f" );

	cout <<"mcPred is: "<< mcPredInt <<endl;
	double fitFraction,fitUncert;

	//_fit->GetResult(signalIdx,fitFraction,fitUncert);
	fitFraction=fitVal;
	fitUncert=fitErr;
	double fitFraction2=allFitVals[signalIdx];
	cout <<"compare the two fractions: " << fitFraction <<" to : " << fitFraction2 <<endl;
	cout <<"we still need "<< 1.0-totalFraction <<" of the data and have " << fitFraction<< " so missing " << 1.0-totalFraction-fitFraction <<endl;
	double miss=1.0-totalFraction-fitFraction;
	cout <<"after scale mc pred is: "<< mcComp->Integral()<<endl;
	cout <<"BR ratio to MC is : "<< mcComp->Integral()/mcSignalIntegral<<endl;
	legend->AddEntry(mcComp,templateLegendNames[SIG_IDX],"f" );
      }

    cout <<"sum of all component integrals: "<< sumOfCompInts <<endl;

    ////
    if(numPions==0)
      {
#ifdef PARTIAL_BOX
	if(leptonId==0)
	  {
	    data->GetYaxis()->SetRangeUser(0,3300*gl_templateScaleFactor);
	    if(channel==1)
	      data->GetYaxis()->SetRangeUser(0,60);
	    if(channel==3)
	      data->GetYaxis()->SetRangeUser(0,20);
	  }
	else
	  {
	    data->GetYaxis()->SetRangeUser(0,1500*gl_templateScaleFactor);
	  }
#else
	if(leptonId==0)
	  data->GetYaxis()->SetRangeUser(0,4000*gl_templateScaleFactor);
	else
	  data->GetYaxis()->SetRangeUser(0,1500*gl_templateScaleFactor);
#endif
      }
    if(numPions==1)
      {
#ifdef PARTIAL_BOX
	if(leptonId==0)
	  {
	    data->GetYaxis()->SetRangeUser(0,3300*gl_templateScaleFactor);
	    if(channel==1)
	      data->GetYaxis()->SetRangeUser(0,60);
	    if(channel==3)
	      data->GetYaxis()->SetRangeUser(0,20);
	  }
	else
	  {
	    data->GetYaxis()->SetRangeUser(0,1500*gl_templateScaleFactor);
	  }
#else
	if(leptonId==0)
	  data->GetYaxis()->SetRangeUser(0,800*gl_templateScaleFactor);
	else
	  data->GetYaxis()->SetRangeUser(0,250*gl_templateScaleFactor);
#endif
      }
    if(numPions==2)
      {
	if(leptonId==0)
	  {
	    if(combineDPiPi)
	      data->GetYaxis()->SetRangeUser(0,450*gl_templateScaleFactor);
	    else
	      data->GetYaxis()->SetRangeUser(0,1600*gl_templateScaleFactor);
	  }
	else
	  {
	    if(combineDPiPi)
	      data->GetYaxis()->SetRangeUser(0,250*gl_templateScaleFactor);
	    else
	      data->GetYaxis()->SetRangeUser(0,800*gl_templateScaleFactor);
	  }
      }
    data->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
    data->SetTitle("");
    data->Draw("Ep");
    predComponents->SetTitle("");
    //    predComponents->SetStats(0);
    //    predComponents->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
    predComponents->Draw("hist same");
    data->SetLineWidth(2);
    data->SetStats(0);
    data->Draw("same Ep");
    legend->Draw();
    sprintf(buffer,"predComp_numPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelString);
    c.SaveAs(buffer);
    sprintf(buffer,"predComp_numPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelString);
    c.SaveAs(buffer);
    sprintf(buffer,"predComp_numPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelString);
    c.SaveAs(buffer);
  }
}
  ////





void loadComponents(TFile* file,TH1F** components, TH1F** summedComponents, int numPions,int leptonId, int channel,  int numComponents, int numFiles)
{
  int channelIdx=channel+1;
  char channelString[1000];
  getChannelString(channel,channelString);

  char buffer[2000];
  for(int iF=0;iF<numFiles;iF++)
    {
      for(int b=0;b<numComponents;b++)
	{
	  sprintf(buffer,"histo_If_%d_b_%d_numPions_%d_leptonId_%d_%s",iF,b,numPions,leptonId,channelString);
	  components[iF*11+b]=(TH1F*)file->Get(buffer);
	}
    }

  sprintf(buffer,"continuum_%dLept_%dPions_%s",leptonId,numPions,channelString);
  cout <<" done with histos, trying to load: " << buffer <<endl;
  summedComponents[iContinuum]=(TH1F*)file->Get(buffer);

  cout <<"pointer : "<< summedComponents[0]<<endl;
  sprintf(buffer,"DDStar_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDDStar]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDDStarPi]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPi_WrongChannel_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDDStarPiWrongChannel]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPiPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDDStarPiPi]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DLnu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDLNu]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[5]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDPiPiLNu]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DStarLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDStarLNu]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DStarPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[8]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DStarPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[8]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DStarPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDStarPiPiLNu]=(TH1F*)file->Get(buffer);


  sprintf(buffer,"CrossFeed_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iDDStarPiCrossFeed]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"OtherBB_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[iOtherBB]=(TH1F*)file->Get(buffer);
  cout <<" done, setting fillColor " <<endl;
  for(int i=0;i<11;i++)
    {
      cout <<"pointer of component: " << i << ": " << summedComponents[i] <<endl;
      cout <<"loaded component " << i << " with counts " << summedComponents[i]->GetEntries()<<endl;
      summedComponents[i]->SetFillStyle(1001);
      summedComponents[i]->SetFillColor(glColorTable[i]->GetNumber());
    }
  //for the color, make it consistent to the other stacks

}



/**
Combine D and D* channels (so B0->Dpilnu with B0->D*pilnu and the same with B->D so we can fit the D* signal and the feeddown at the same time.
The returned components are one less, since we combine two and the D* missing mass templates are shifted
summedComponents_out will have 2* summedComponentsD -1 size (everything double, but remove both x-feed fields and combine (DStar x-feed doesn't make sense, only there for symmetry)
 */
void combineChannels(TH1F** summedComponentsD_in,  TH1F** summedComponentsDStar_in, TH1F** summedComponents_out, int numComponents, int DChannelIdx, int numPions)
{
  //question: do we get in trouble with binning, mass range assumptions in later fits? Need to keep that flexible...
  
  //already the correct index because this has been computed from the combined channel (so starts at 0)

  cout <<" using channel idx: " << DChannelIdx <<endl;
  //the corresponding dStar channel should be the next one
  int DStarChannelIdx=DChannelIdx+1;
  Double_t upperCutD=upperCut[DChannelIdx];
  Double_t lowerCutD=lowerCut[DChannelIdx];
  int numBinsD=numBins[DChannelIdx];
  Double_t upperCutDStar=upperCut[DStarChannelIdx];
  Double_t lowerCutDStar=lowerCut[DStarChannelIdx];
  int numBinsDStar=numBins[DStarChannelIdx];
  char buffer[500];
  cout <<"numBinsD: "<< numBinsD <<" upperCutD: " << upperCutD <<" lowerCut: "<< lowerCutD <<", bins DS: " << numBinsDStar << " upper DS: " << upperCutDStar <<" lower DS: "<< lowerCutDStar <<endl;
  //for SIG and x_feed we put it in the same template, separate templates for the others

  //run over components twice, first for D, where we add the upward crossfeed (should be non-existent) to the signal and the D* signal to the x-feed
  //then a second run for the D* where we add everything but the signal and xfeed
  int counter=-1;
  for(int i=0;i<numComponents;i++)
    {
      counter++;
      if(i==SIG_IDX)
	{
	  sprintf(buffer,"combined_%s_%s",summedComponentsD_in[SIG_IDX]->GetName(),summedComponentsDStar_in[iDDStarPiCrossFeed]->GetName());
	  summedComponents_out[counter]=new TH1F(buffer,buffer,numBinsD+numBinsDStar,lowerCutD,upperCutD+(upperCutDStar-lowerCutDStar));
	  summedComponents_out[counter]->Sumw2();
	  addShiftedHistos(summedComponentsD_in[SIG_IDX],summedComponentsDStar_in[iDDStarPiCrossFeed],summedComponents_out[counter],numBinsD,numBinsDStar);
	  summedComponents_out[counter]->SetFillStyle(1001);
	  summedComponents_out[counter]->SetFillColor(glColorTable[i]->GetNumber());
	  continue;
	}
      bool isCrossFeed=(i==iDDStarPiCrossFeed && numPions >0 )|| (numPions==0 && i==iDStarLNu);
      if(isCrossFeed)
	{
	  sprintf(buffer,"combined_%s_%s",summedComponentsD_in[iDDStarPiCrossFeed]->GetName(),summedComponentsDStar_in[SIG_IDX]->GetName());
	  summedComponents_out[counter]=new TH1F(buffer,buffer,numBinsD+numBinsDStar,lowerCutD,upperCutD+(upperCutDStar-lowerCutDStar));
	  summedComponents_out[counter]->Sumw2();
	  if(numPions==0)
	    addShiftedHistos(summedComponentsD_in[iDStarLNu],summedComponentsDStar_in[iDStarLNu],summedComponents_out[counter],numBinsD,numBinsDStar);
	  else
	    addShiftedHistos(summedComponentsD_in[iDDStarPiCrossFeed],summedComponentsDStar_in[SIG_IDX],summedComponents_out[counter],numBinsD,numBinsDStar);
	  summedComponents_out[counter]->SetFillStyle(1001);
	  summedComponents_out[counter]->SetFillColor(glColorTable[i]->GetNumber());
	  continue;  
	}
      //all other cases just copy the D templates in the lower half
      if((!isCrossFeed) && (i != SIG_IDX))
	{
	  sprintf(buffer,"combined_D_%s",summedComponentsD_in[i]->GetName());
	  summedComponents_out[counter]=new TH1F(buffer,buffer,numBinsD+numBinsDStar,lowerCutD,upperCutD+(upperCutDStar-lowerCutDStar));
	  summedComponents_out[counter]->Sumw2();
	  copyHisto(summedComponentsD_in[counter],summedComponents_out[counter],numBinsD,numBinsD+numBinsDStar);
	  summedComponents_out[counter]->SetFillStyle(1001);
	  summedComponents_out[counter]->SetFillColor(glColorTable[i]->GetNumber());
	}
    }
  int ommitted=0;
  //now run over the DStar components, but leave out SIG_IDX and cross-feed because we already added that earlier (and don't forget to shift up)
  for(int i=0;i<numComponents;i++)
      {
	bool isCrossFeed=(i==iDDStarPiCrossFeed && numPions >0 )|| (numPions==0 && i==iDStarLNu);
	if((i==SIG_IDX) || (isCrossFeed))
	  {
	    continue;
	  }
	counter++;
	sprintf(buffer,"combined_DStar_%s",summedComponentsDStar_in[i]->GetName());
	summedComponents_out[counter]=new TH1F(buffer,buffer,numBinsD+numBinsDStar,lowerCutD,upperCutD+(upperCutDStar-lowerCutDStar));
	summedComponents_out[counter]->Sumw2();
	copyHisto(summedComponentsDStar_in[i],summedComponents_out[counter],numBinsDStar,numBinsD+numBinsDStar,numBinsD);
	//take the same colors and jump over signal etc..
	summedComponents_out[counter]->SetFillStyle(1001);
	summedComponents_out[counter]->SetFillColor(glColorTable[i]->GetNumber());
    }
}



//trees: the input trees for the 4 MC files, components: The components for each file, so 4*11, summedComponents: The same, but summed over files
void getMCComponents(TTree** trees, TH1F** components, TH1F** summedComponents, int numPions,int leptonId, int channel)
{
  int channelIdx=channel+1;
  char channelSelectionData[1000];
  char channelSelectionDataAndMC[1000];
  char wrongChannelSelection[1000];
  char crossFeedSelection[1000];

    //void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions)
  getChargeAndStarSelection(channelSelectionData,channel,false,numPions);
  getChargeAndStarSelection(channelSelectionDataAndMC,channel,true,numPions);
  getChargeAndStarSelection(wrongChannelSelection,channel,true,numPions,true);
  getFeedDown(crossFeedSelection,channel,false,numPions);

  char channelString[500];
  getChannelString(channel,channelString);

  cout<<" getting pres for channel " << channelString <<" lowerCut: " << lowerCut[channelIdx] << " upperCut " << upperCut[channelIdx] <<" num bins: " << numBins[channelIdx] <<endl;
  cout <<"getting MC components...(comb) " <<endl;
  char histoName[2009];
  char drawCommand[2000];
  //just getting the MC components, so only looking at 4 files
  int numFiles=4;
  //now 11 components because we split off the d* downfeed
  int numComponents=11;
  cout <<"-7"<<endl;
  char buffer[2000];
  sprintf(buffer,"continuum_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hContinuum=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hContinuum->Sumw2();
  //  hContinuum->Fill(0.0,0.0);
  //		  cout <<"second axis limits: "<< summedHistos[b]->GetXaxis()->GetXmin() <<" to " << summedHistos[b]->GetXaxis()->GetXmax()<<endl;
  cout << setprecision(15);
  cout <<"cont min: "<<   hContinuum->GetXaxis()->GetXmin()<<" max: "<< hContinuum->GetXaxis()->GetXmax()<<endl;

  sprintf(buffer,"DDStar_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStar=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  //  hDDStar->Fill(0.0,0.000001);
  cout <<"DDStar min: "<<   hContinuum->GetXaxis()->GetXmin()<<" max: "<< hContinuum->GetXaxis()->GetXmax()<<endl;
  hDDStar->Sumw2();
  sprintf(buffer,"DDStarPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStarPi=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDDStarPi->Sumw2();
  cout <<"-5"<<endl;
  sprintf(buffer,"DDStarPi_WrongChannel_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStarPiWrongChannel=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDDStarPiWrongChannel->Sumw2();
  //this should be the same channel. For consistency we will also have a cross feed channel for the D* channels, but this should be zero there
  sprintf(buffer,"CrossFeed_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStarPiCrossFeed=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDDStarPiCrossFeed->Sumw2();
  cout <<"-4"<<endl;
  sprintf(buffer,"DDStarPiPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStarPiPi=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDDStarPiPi->Sumw2();
  sprintf(buffer,"DLnu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDlNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDlNu->Sumw2();
  sprintf(buffer,"DPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDPilNu->Sumw2();
  sprintf(buffer,"DPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDPiPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDPiPilNu->Sumw2();
  sprintf(buffer,"DStarLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDStarlNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDStarlNu->Sumw2();
  sprintf(buffer,"DStarPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDStarPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDStarPilNu->Sumw2();
  sprintf(buffer,"DStarPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDStarPiPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hDStarPiPilNu->Sumw2();
  sprintf(buffer,"OtherBB_%dLept_%dPions_%s",leptonId,numPions,channelString);

  TH1F* hOtherBB=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  hOtherBB->Sumw2();
  char* selections[11];
  char bufferDDStar[2000];
  char bufferDDStarPi[2000];
  char bufferDDStarPiWrongChannel[2000];
  char bufferDDStarPiCrossFeed[2000];
  char bufferDDStarPiPi[2000];

  char bufferDlNu[2000];
  char bufferDPilNu[2000];
  char bufferDPiPilNu[2000];

  char bufferDStarlNu[2000];
  char bufferDStarPilNu[2000];
  char bufferDStarPiPilNu[2000];

  char bufferAll[2000];
  char bufferNoSelection[2000];
  char corrBuffer[2000];
  cout <<"-1" <<endl;
  addCorrections(corrBuffer);
  cout <<"1" <<endl;
  if(leptonId!=0)
    {
      if(numPions==0)
	{
	  sprintf(buffer,"%s  tagCorr*" P0STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		//	  sprintf(buffer,"%s  " P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
	}
      if(numPions==1)
	{
	  sprintf(buffer,"%s  tagCorr*" P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		//	  sprintf(buffer,"%s  " P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
	}
      if(numPions==2)
	{
	  sprintf(buffer,"%s  tagCorr*(" P2STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
	}
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d  ",numPions,leptonId);
    }
  else
    {
      if(numPions==0)
	{	
	  sprintf(buffer,"%s tagCorr*" P0STRING "  && bestBCharge==((-1)*systemCharge) && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  //	  sprintf(buffer,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
	}
      if(numPions==1)
	{	
	  sprintf(buffer,"%s tagCorr*" P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  //	  sprintf(buffer,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
	}
      if(numPions==2)
	{
	  sprintf(buffer,"%s tagCorr*(" P2STRING "  && bestBCharge==((-1)*systemCharge) && %s",corrBuffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
	}
      //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0   ",numPions);
    }

  sprintf(bufferDDStar,"%s &&   foundAnyDDoubleStar==1 && sig_numPions==0 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1 && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
  //sprintf(bufferDDStar,"%s && foundAnyDDoubleStar==1)",buffer);
  //this is the signal, so here we put the channel selection data and mc (even though the data part should be redundant)
  if(numPions==1 || numPions==2)
    {
      sprintf(bufferDDStarPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu && %s)",buffer,channelSelectionDataAndMC);
    }
  if(numPions==0)
    {
      //for the x-check with Robin we don't care if we pick up wrong combinations..., so no need for the channel selection
      sprintf(bufferDDStarPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
    }

  if(numPions==1 || numPions==2)
    sprintf(bufferDDStarPiWrongChannel,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu && %s)",buffer,wrongChannelSelection);
  if(numPions==0)
    {
      //use this fo rthe D l nu stuff
      sprintf(bufferDDStarPiWrongChannel,"%s && sig_DLNu &&  %s)",buffer,wrongChannelSelection);
    }


//////  //since for the no-pion case we want to keep the same cross-feed index
//////  if(numPions==0)  
//////    {
//////      sprintf(bufferDDStarPiCrossFeed,"%s  && sig_numPions==0 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0 && sig_numLeptons==1 && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu &&  %s)",buffer,crossFeedSelection);
//////    }
//////  else
  if(numPions==1 || numPions==2)
    {
      sprintf(bufferDDStarPiCrossFeed,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0 && sig_numLeptons==1 && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu &&  %s)",buffer,crossFeedSelection);
    }
  if(numPions==0)
    {
      //put all D** in one plot
      sprintf(bufferDDStarPiCrossFeed,"1!=1");
    }
  sprintf(bufferDDStarPiPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==2 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);

  if(numPions==1 || numPions==2)
    {
      sprintf(bufferDlNu,"%s && sig_DLNu )",buffer);
    }
  if(numPions==0)
    {
      //if this is actually our signal, make sure that the channel selection is correct
      sprintf(bufferDlNu,"%s && sig_DLNu && %s)",buffer,channelSelectionDataAndMC);
    }


  sprintf(bufferDPilNu,"%s && sig_DPiLNu&& !sig_DLNu  )",buffer);
  sprintf(bufferDPiPilNu,"%s && sig_DPiPiLNu && !sig_DLNu && !sig_DPiLNu)",buffer);
  sprintf(bufferDStarlNu,"%s &&  sig_DStarLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu)",buffer);
  sprintf(bufferDStarPilNu,"%s  && sig_DStarPiLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu)",buffer);
  sprintf(bufferDStarPiPilNu,"%s  && sig_DStarPiPiLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu)",buffer);

  //weird, the test for sig_X that come first seem useless.
  sprintf(bufferAll,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu|| sig_DStarPiPiLNu) && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);

  sprintf(bufferNoSelection,"%s)",buffer);
      //sprintf(bufferNoSelection,"%s)",bufferAll);
  //    sprintf(bufferAll,"%s)",buffer);
  cout <<"selection for DDstar: "<< bufferDDStar <<endl<<" DDstarPi: "<< bufferDDStarPi<<endl<< " DDstarPiPi: "<< bufferDDStarPiPi <<endl;
  cout <<" DlNu: " << bufferDlNu<<endl <<" DPilNu:"<< bufferDPilNu << endl << " DPiPilNu: "<< bufferDPiPilNu <<endl;
  cout <<" all: "<< bufferAll <<endl;

  selections[iNoSelection]=bufferNoSelection;
  selections[iDDStar]=bufferDDStar;
  selections[iDDStarPi]=bufferDDStarPi;
  selections[iDDStarPiWrongChannel]=bufferDDStarPiWrongChannel;
  selections[iDDStarPiPi]=bufferDDStarPiPi;
  selections[iDLNu]=bufferDlNu;
  //  selections[5]=bufferDPilNu;
  selections[iDPiPiLNu]=bufferDPiPilNu;

  selections[iDStarLNu]=bufferDStarlNu;
  //  selections[8]=bufferDStarPilNu;
  selections[iDStarPiPiLNu]=bufferDStarPiPilNu;
  selections[iDDStarPiCrossFeed]=bufferDDStarPiCrossFeed;
  selections[iAll]=bufferAll;

  //
  // summedHistos: not differentiated between mixed and charged
  //in the end we use something like 		  summedHistos[b+1]->Add(result from selection[b]) (and continuum is treated extra)
  //
  TH1F**  summedHistos=summedComponents;
  summedHistos[iContinuum]=hContinuum;
  summedHistos[iDDStar]=hDDStar;
  summedHistos[iDDStarPi]=hDDStarPi;
  summedHistos[iDDStarPiWrongChannel]=hDDStarPiWrongChannel;
  summedHistos[iDDStarPiPi]=hDDStarPiPi;
  summedHistos[iDLNu]=hDlNu;
  //  summedHistos[5]=hDPilNu;
  summedHistos[iDPiPiLNu]=hDPiPilNu;
  summedHistos[iDStarLNu]=hDStarlNu;
  //  summedHistos[8]=hDStarPilNu;
  summedHistos[iDStarPiPiLNu]=hDStarPiPilNu;
  summedHistos[iDDStarPiCrossFeed]=hDDStarPiCrossFeed;
  summedHistos[iOtherBB]=hOtherBB;
  
  for(int iF=0;iF<numFiles;iF++)
    {
      cout <<" iF: "<< iF <<" numComps: " << numComponents<< " numFiles: "<< numFiles <<endl;
      int allCounts=0;
      int noSelCounts=0;
      for(int b=0;b<numComponents;b++)
	{
	  char currentSelection[1000];
	  sprintf(currentSelection,"%s",selections[b]);
	  //	  insertSignalMCWeighting(currentSelection, b,iF);
	  insertSignalMCWeighting(currentSelection,iF);
	  cout <<"b: "<< b << endl;
	  cout <<"creating histo independently " << endl;
	  sprintf(histoName,"histo_If_%d_b_%d_numPions_%d_leptonId_%d_%s",iF,b,numPions,leptonId,channelString);
	  TH1F* h=new TH1F(histoName,histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	  h->Sumw2();
	  //	  sprintf(drawCommand,"mNu2 >> %s",histoName);
	  	  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	  cout <<"draw command: " << drawCommand << ", selections: " << (char*) currentSelection<<endl;
	  cout <<"drawing tree " << iF <<endl;
	  ///need to add weighting for signal MC

	  ////
	  int counts=trees[iF]->Draw(drawCommand,(char*)currentSelection);
	  if(b>0)
	    allCounts+=counts;
	  if(b==0)
	    noSelCounts=counts;
	    if(b==10)
	  //	  if(b==9)
	      cout <<"all counts so far: "<< allCounts <<" no selection counts: "<< noSelCounts <<" difference: "<< noSelCounts-allCounts <<endl;
	  cout <<"got " << counts <<" counts selected " <<endl;
	  //do we have to clone this before we can return the histo?
	  //if counts == 0, this leads to an error in the add function (axis limit)

	  TH1F* result=(TH1F*)gDirectory->Get(histoName);
	  result->SetFillStyle(1001);
	  result->SetFillColor(glColorTable[b]->GetNumber());
	  components[iF*11+b]=result;
	  //	  if(b!=10)
	    {
	      //uds or charm, but only use the noSelection
	      if(iF>=2 && b==iNoSelection)
		{
		  cout <<"c adding histo with: "<< result->GetNbinsX() <<", "<< result->GetBinCenter(1) <<" to " << result->GetBinCenter(result->GetNbinsX())<<endl;
		  cout <<" to histo with: "<< hContinuum->GetNbinsX() <<", "<< hContinuum->GetBinCenter(1) <<" to " << hContinuum->GetBinCenter(summedHistos[b]->GetNbinsX())<<endl;
		  // counts==0 leads to problems with Axis ranges
		  if(counts >0)
		    hContinuum->Add(result);
		  cout <<"done " <<endl;
		  hContinuum->SetFillStyle(1001);
		  hContinuum->SetFillColor(glColorTable[0]->GetNumber());
		}
	      //mixed or charged. The b==0 case is the noSelection case, which only makes sense for the continuum
	      if(iF<2 && b>iNoSelection)
		{
		  //this is because the selections and the summedHistos are offset by one
		  //
		  if(result->GetZaxis()->GetXmin() != summedHistos[b]->GetZaxis()->GetXmin())
		    cout <<"z aaxis lower limit different " << endl;
		  if(result->GetZaxis()->GetXmax() != summedHistos[b]->GetZaxis()->GetXmax())
		    cout <<"z aaxis upper limit different " << endl;
		  if(result->GetYaxis()->GetXmin() != summedHistos[b]->GetYaxis()->GetXmin())
		    cout <<"y aaxis lower limit different " << endl;
		  if(result->GetYaxis()->GetXmax() != summedHistos[b]->GetYaxis()->GetXmax())
		    cout <<"y aaxis upper limit different " << endl;
		  if(result->GetXaxis()->GetXmin() != summedHistos[b]->GetXaxis()->GetXmin())
		    cout <<"x aaxis lower limit different " << endl;
		  if(result->GetXaxis()->GetXmax() != summedHistos[b]->GetXaxis()->GetXmax())
		    cout <<"x aaxis upper limit different " << endl;
		  cout <<setprecision(15);
		  cout <<" first zaxis limits: "<< result->GetZaxis()->GetXmin() <<" to " << result->GetZaxis()->GetXmax()<<endl;
		  cout <<"secondy zxis limits: "<< summedHistos[b]->GetZaxis()->GetXmin() <<" to " << summedHistos[b]->GetZaxis()->GetXmax()<<endl;
		  cout <<" first yaxis limits: "<< result->GetYaxis()->GetXmin() <<" to " << result->GetYaxis()->GetXmax()<<endl;
		  cout <<"secondy axis limits: "<< summedHistos[b]->GetYaxis()->GetXmin() <<" to " << summedHistos[b]->GetYaxis()->GetXmax()<<endl;
		  cout <<" first axis limits: "<< result->GetXaxis()->GetXmin() <<" to " << result->GetXaxis()->GetXmax()<<endl;
		  cout <<"second axis limits: "<< summedHistos[b]->GetXaxis()->GetXmin() <<" to " << summedHistos[b]->GetXaxis()->GetXmax()<<endl;
		  cout <<" adding  histo with: "<< result->GetNbinsX() <<", "<< result->GetBinCenter(1) <<" to " << result->GetBinCenter(result->GetNbinsX())<<endl;
		  cout <<" to summed histo with: "<< summedHistos[b]->GetNbinsX() <<", "<< summedHistos[b]->GetBinCenter(1) <<" to " << summedHistos[b]->GetBinCenter(summedHistos[b]->GetNbinsX())<<endl;
		  cout <<setprecision(3);
		  // counts==0 leads to problems with Axis ranges
		  if(counts>0)
		    {
		      cout <<"summed histos has: "<< summedHistos[b]->GetEntries() <<" and result: "<< result->GetEntries() <<endl;
		      summedHistos[b]->Add(result);
		    }
		  cout <<"done " <<endl;


		  //for the color, make it consistent to the other stacks
		  summedHistos[b]->SetFillStyle(1001);
		  summedHistos[b]->SetFillColor(glColorTable[b]->GetNumber());
		}
	    }
	}
    }
  cout <<"done with getMCComponents" <<endl;
};





void saveStack(TH1F** components, TH1F** summedComponents, int numPions, int leptonId, int channel, bool isCombinedChannel)
{
  int channelIdx=channel+1;
  int printOrder[]={0,10,1,3,4,5,9,7,8,6,2};
  //for combined
  if(numPions==1)
    {
      //      printOrder[10]=3;
      //      printOrder[3]=2;
      //      printOrder[9]=4;
      //      printOrder[4]=6;
    }

  cout <<"in save stack.. " <<endl;
  //this is for "all", i.e. not mixed, charged separation
  THStack all;
  char* legendNames[20];

  int numFiles=4;
  char* fileNames[numFiles];
  for(int i=0;i<numFiles;i++)
    {
      fileNames[i]=new char[400];
    }


  for(int i=0;i<20;i++)
    {
      legendNames[i]=new char [500];
      allLegendNames[i]=new char[500];
    }

  char channelString[500];
  getChannelString(channel,channelString);


  sprintf(legendNames[iNoSelection]," no Selection");
  sprintf(legendNames[iDDStar],"B #rightarrow D Double Star X #rightarrow D^{(*)} l #nu (no non-res Dn#pi l#nu)");
  sprintf(legendNames[iDDStarPi],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu (no non-res Dn#pi l#nu), %s",channelString);
  sprintf(legendNames[iDDStarPiWrongChannel],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu (no non-res Dn#pi l#nu)-->Wrong Channel");
  sprintf(legendNames[4],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi #pi l #nu (no non-res Dn#pi l#nu)");
  sprintf(legendNames[5],"D l #nu");
  //  sprintf(legendNames[5],"D #pi l #nu");
  sprintf(legendNames[6],"D #pi #pi l #nu");

  sprintf(legendNames[7],"D* l #nu");
  //  sprintf(legendNames[8],"D* #pi l #nu");
  sprintf(legendNames[8],"D* #pi #pi l #nu");
  sprintf(legendNames[9],"B #rightarrow D Double Star X #rightarrow D^{*} #pi l #nu -->feeddown ");
  sprintf(legendNames[10]," no D(*)n #pi l#nu ");

  sprintf(allLegendNames[0],"Continuum");
  sprintf(allLegendNames[1],"B #rightarrow D Double Star X #rightarrow D^{(*)}  l #nu");
  sprintf(allLegendNames[2],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu, %s", channelString);
  sprintf(allLegendNames[3],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu-->Wrong Channel");
  sprintf(allLegendNames[4],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi #pi l #nu");
  sprintf(allLegendNames[5],"D l #nu");
  //  sprintf(allLegendNames[5],"D #pi l #nu");
  sprintf(allLegendNames[6],"D #pi #pi l #nu");
  sprintf(allLegendNames[7],"D* l #nu");
  //  sprintf(allLegendNames[8],"D* #pi l #nu");
  sprintf(allLegendNames[8],"D* #pi #pi l #nu");
  sprintf(allLegendNames[9],"B #rightarrow D Double Star X #rightarrow D^{*} #pi l #nu-->feedDown");
  sprintf(allLegendNames[10]," other B B ");

  //for the combination, make sure that we fill the other legend names as well
  int counter = -1;
  for(int i=0;i<11;i++)
    {
      if(i==SIG_IDX || i==iDDStarPiCrossFeed)
	continue;
      counter++;
      cout <<"trying to populate legend names: "<< counter+11 <<" with " << i%11 <<endl;
      sprintf(allLegendNames[counter+11], "%s",allLegendNames[i%11]);
    }


  sprintf(fileNames[0],"mixed");
  sprintf(fileNames[1],"charged");
  sprintf(fileNames[2],"uds");
  sprintf(fileNames[3],"charm");

  //
  // the stacks for the charm, uds, mixed, charged differentiated plots
  //
  THStack* stacks[4];
  char histoName[2009];
  char outFileName[2000];

  char buffer[2000];


  
  for(int iF=0;iF<4;iF++)
    {
      cout <<"file/tree nr: " << iF << endl;
      char bufferF[200];
      sprintf(bufferF,"stackF_%d",iF);
      stacks[iF]=new THStack(bufferF,bufferF);
      TLegend* legend=0;

      if(channel==0 || channel==2)
	{
	  legend=new TLegend(0.12,0.6,0.5,0.99);
	  legend->SetNColumns(2);
	}
      else
	{
	  legend =new TLegend(0.6,0.7,0.89,0.99);
	}
      //      int allCounts=0;
      for(int b=0;b<11;b++)
	{
	  int index=printOrder[b];
	  cout <<"b: " << b <<" index: " << index<<endl;
	  sprintf(histoName,"histo_If_%d_index_%d_numPions_%d_leptonId_%d_%s",iF,index,numPions,leptonId,channelString);
	  sprintf(outFileName,"%s.png",histoName);
	  //	  components[iF*9+b]=result;
	  cout <<"grabbing component.." << iF*11+index <<endl;
	  TH1F* result=(TH1F*)components[iF*11+index];
	  cout <<"done that .." <<endl;
	  //  if(counts>0)
	    {
	      //other bb is not 10
	      //used to be 10... but 10 is bufferAll, so all but everything else, I think 0 is the no selection (so all)...
	      if(index==0)
		{
		  //		  stacks[iF]->Add(result,"nostack");
		  //		legend->AddEntry(result,legendNames[b],"f");
		}
	      else
		{
		  cout <<"add to stack " << endl;
		  stacks[iF]->Add(result);
		  cout << "Add entry " << legendNames[index] <<endl;
		  legend->AddEntry(result,legendNames[index],"f");
		}


	    }
	  if(!result)
	    {
	      cout <<"null pointer returned" <<endl;
	    }
	  else
	    {
	      TCanvas c;
	      result->Draw("hist");
	      c.Update();
	      c.SaveAs(outFileName);
	      sprintf(outFileName,"%s.pdf",histoName);
	      c.SaveAs(outFileName);
	      sprintf(outFileName,"%s.eps",histoName);
	      c.SaveAs(outFileName);
	    }
	}
      cout <<"move on..." <<endl;
      char stackName[300];

      TCanvas c2;
      sprintf(stackName,"Stack_%s_%d pions_%d_leptonId_%s",fileNames[iF],numPions,leptonId,channelString);
      //      stacks[iF]->SetTitle(stackName);
      stacks[iF]->SetTitle("");
      //      stacks[iF]->SetStats(0);
      stacks[iF]->Draw("hist");
      legend->Draw();
      c2.Update();
      stacks[iF]->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
      c2.Modified();
      sprintf(stackName,"Stack_%s_%d_pions_leptonId_%d_%s.png",fileNames[iF],numPions,leptonId,channelString);
      c2.SaveAs(stackName);
      sprintf(stackName,"Stack_%s_%d_pions_leptonId_%d_%s.pdf",fileNames[iF],numPions,leptonId,channelString);
      c2.SaveAs(stackName);
      sprintf(stackName,"Stack_%s_%d_pions_leptonId_%d_%s.eps",fileNames[iF],numPions,leptonId,channelString);
      c2.SaveAs(stackName);
    }


  ///the summed, all stuff...
  TH1F**  summedHistos=summedComponents;
  TLegend* legend =0;
  if(channel==0 || channel==2)
    {
      legend =new TLegend(0.01,0.7,0.89,0.99);
      //      legend=new TLegend(0.12,0.6,0.5,0.99);-->this would put it in the upper left which interferes with the peak if we chose a window from -0.5 to 1.5 or 2
      legend->SetNColumns(2);
      //      gStyle->SetLegendFont(22);
      cout <<"current legend text size: " << legend->GetTextSize()<<endl;///evaluates to 0
      //of course this doesn't help if it is zero...
      //      gStyle->SetLegendTextSize(legend->GetTextSize()*2);
      legend->SetTextSize(legend->GetTextSize()*2);
      if(legend->GetTextSize()==0)
	legend->SetTextSize(0.03);
	      //      	gStyle->SetLegendTextSize(0.03);
      all.SetMinimum(0);
      all.SetMaximum(1500);
    }
  else
    {
	  legend =new TLegend(0.6,0.7,0.89,0.99);
    }

  TCanvas t;
  //do the D* non-signal channels first, so they are not on top of the signal
  if(isCombinedChannel)
    {
      int count=-1;
      for(int i=0;i<11;i++)
	{

	  //these were joined
	  if(i==SIG_IDX || i== iDDStarPiCrossFeed)
	    continue;
	  count++;
	  //indices are meaningless here, because we are missing two fields
	  //	  int index=printOrder[i];
	  cout <<"accessing summedHisto " << count+11<<endl;
	  sprintf(buffer,"stackTest_DStarJoin_%d.png",i);
	  summedHistos[count+11]->SetFillStyle(1001);
	  summedHistos[count+11]->Draw();
	  t.SaveAs(buffer);

	  //	  summedHistos[count+11]->Update();
	  all.Add(summedHistos[count+11]);
	  //legends should be the same, not clear if one can use one legend twice
	  legend->AddEntry(summedHistos[count+11],allLegendNames[i],"f");
	}
    }
  for(int i=0;i<11;i++)
    {
      int index=printOrder[i];
      cout <<"adding summed histo with " << summedHistos[index]->GetEntries() << " counts, index: " << index <<" i: " << i <<endl;
      sprintf(buffer,"stackTest_%d.png",index);
      summedHistos[index]->SetFillStyle(1001);
      summedHistos[index]->Draw();
      t.SaveAs(buffer);
      //      summedHistos[index]->Update();
      all.Add(summedHistos[index]);

      legend->AddEntry(summedHistos[index],allLegendNames[index],"f");
    }

  TCanvas c;
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s",numPions,leptonId,channelString);
  //  all.SetTitle(buffer);
  all.SetTitle("");
  //  all.GetYaxis()->SetRangeUser(0,1000);
  //  all.SetStats(0);
  all.Draw("hist");
  legend->Draw();
  c.Update();
  all.GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
  c.Modified();
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.png",numPions,leptonId,channelString);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.root",numPions,leptonId,channelString);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.C",numPions,leptonId,channelString);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.pdf",numPions,leptonId,channelString);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.eps",numPions,leptonId,channelString);
  c.SaveAs(buffer);
}


void doSidebandComparison(TTree* mcTree, TTree* dataTree,int leptonId, int numPions, TH1F** lowerSidebandMC, TH1F** upperSidebandMC, TH1F** lowerSidebandData, TH1F** upperSidebandData, int channel)
{

#ifdef Huschle_Signal_Weighting
  float mixedFactor=1.0/0.852255195;
  float chargedFactor=1.0/0.9266338302;
#else
  float mixedFactor=1.0;
  float chargedFactor=1.0;
#endif
  float huschleLumiFactor=mixedFactor;
  if(channel==-1)//for all, should probably take the mean
    huschleLumiFactor=(mixedFactor+chargedFactor)/2;

  if(channel>1)
    huschleLumiFactor=chargedFactor;
  huschleLumiFactor-=1;


  char huschleMC_lumi_corr[1000];
  sprintf(huschleMC_lumi_corr,"(1+foundAnyDDoubleStar*%f)*",huschleLumiFactor);
  cout <<"using additional lumi factor: "<< huschleMC_lumi_corr<<endl;

  int channelIdx=channel+1;
  //2.0 to 0.5
  //  float upperSidebandTop=3.0;
  //  float upperSidebandBottom=-1.5;

    float upperSidebandTop=2.0;
  //  float upperSidebandTop=3.0;
    float upperSidebandBottom=0.5;
  //  float upperSidebandBottom=-1.5;

  float lowerSidebandTop=-0.5;
  float lowerSidebandBottom=-1.0;

  char upperSBSelection[2000];
  char lowerSBSelection[2000];

  char upperSBSelectionData[2000];
  char lowerSBSelectionData[2000];

  char histoName[200];
  char drawCommand[200];
  char buffer[2000];
  char tmpBuffer[2000];

  char channelString[500];
  char channelSelectionData[1000];

  getChannelString(channel,channelString);
  cout <<"getting charge and star selection for channel: "<< channel <<" pions: " << numPions <<endl;
  getChargeAndStarSelection(channelSelectionData,channel,false,numPions);
  cout <<"string: "<< channelSelectionData <<endl;

  addCorrections(tmpBuffer);

  sprintf(buffer,"%s %s",huschleMC_lumi_corr,tmpBuffer);
  //select sidebands from all (need to redo all components because we select a different range
  if(leptonId!=0)
    {
      if(numPions==0)
	{
	  sprintf(upperSBSelection,"%s tagCorr*" P0STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	  //	sprintf(upperSBSelection,"%s " P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	}
      if(numPions==1)
	{
	  sprintf(upperSBSelection,"%s tagCorr*" P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	  //	sprintf(upperSBSelection,"%s " P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	}
      if(numPions==2)
	{
	sprintf(upperSBSelection,"%s tagCorr*(" P2STRING "&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	}

      if(numPions==0)
	sprintf(upperSBSelectionData,P0STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ", upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
      if(numPions==1)
	sprintf(upperSBSelectionData,P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ", upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
      if(numPions==2)
	sprintf(upperSBSelectionData," (" P2STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ", upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
     
      if(numPions==0)	
	sprintf(lowerSBSelection,"%s tagCorr*" P0STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);	
      if(numPions==1)	
	sprintf(lowerSBSelection,"%s tagCorr*" P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);	
      if(numPions==2)
	sprintf(lowerSBSelection,"%s tagCorr*(" P2STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);

      if(numPions==0)
	sprintf(lowerSBSelectionData,P0STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);
      if(numPions==1)
	sprintf(lowerSBSelectionData,P1STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);
      if(numPions==2)
	sprintf(lowerSBSelectionData," (" P2STRING "&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
      if(numPions==0)
	sprintf(upperSBSelection,"%s tagCorr*" P0STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
      if(numPions==1)
	sprintf(upperSBSelection,"%s tagCorr*" P1STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
      if(numPions==2)
	sprintf(upperSBSelection,"%s tagCorr*(" P2STRING " && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);

      if(numPions==0)
	sprintf(upperSBSelectionData,"%s " P0STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
      if(numPions==1)
	sprintf(upperSBSelectionData,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
      if(numPions==2)
	sprintf(upperSBSelectionData,"%s (" P2STRING " && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
        
      if(numPions==0)
	sprintf(lowerSBSelection,"%s tagCorr*" P0STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);
      if(numPions==1)
	sprintf(lowerSBSelection,"%s tagCorr*" P1STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);
      if(numPions==2)
      sprintf(lowerSBSelection,"%s tagCorr*(" P2STRING " && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);

      if(numPions==0)
	sprintf(lowerSBSelectionData,"%s " P0STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);
      if(numPions==1)
	sprintf(lowerSBSelectionData,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);
      if(numPions==2)
	sprintf(lowerSBSelectionData,"%s (" P2STRING " && bestBCharge==((-1)*systemCharge) && %s ) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);

      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
    }


  cout <<"upper selection: "<< upperSBSelection << " lower selection : "<< lowerSBSelection<<endl;
  sprintf(histoName,"upperSideBandMC_numPions_%d_leptonId_%d_%s",numPions,leptonId, channelString);
  TH1D* histo_USMC=new TH1D(histoName,histoName,numBinsSBTop,upperSidebandBottom,upperSidebandTop);
  sprintf(drawCommand,"mNu2 >> %s",histoName);
  int counts=mcTree->Draw(drawCommand,(char*)upperSBSelection);
  cout <<"got " << counts <<" counts from mc upperSB selected for lepton ID "<< leptonId << " channel: "<< channelString<<endl;
  *upperSidebandMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"lowerSideBandMC_numPions_%d_leptonId_%d_%s",numPions,leptonId, channelString);
  TH1D* histo_LSMC=new TH1D(histoName,histoName,numBinsSB,lowerSidebandBottom,lowerSidebandTop);
  sprintf(drawCommand,"mNu2 >> %s",histoName);
  counts=mcTree->Draw(drawCommand,(char*)lowerSBSelection);
  cout <<"got " << counts <<" counts from mc lowerSB selected " <<endl;
  *lowerSidebandMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"upperSideBandData_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_USData=new TH1D(histoName,histoName,numBinsSBTop,upperSidebandBottom,upperSidebandTop);
  sprintf(drawCommand,"mNu2 >> %s",histoName);
  cout <<"about to get data with string: "<< drawCommand <<endl;
  counts=dataTree->Draw(drawCommand,(char*)upperSBSelectionData);
  cout <<"got " << counts <<" counts from data upperSB selected " <<endl;
  *upperSidebandData=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"lowerSideBandData_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_LSData=new TH1D(histoName,histoName,numBinsSB,lowerSidebandBottom,lowerSidebandTop);
  sprintf(drawCommand,"mNu2 >> %s",histoName);
  //  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSB,lowerSidebandBottom,lowerSidebandTop);
  counts=dataTree->Draw(drawCommand,(char*)lowerSBSelectionData);
  cout <<"got " << counts <<" counts from data lowerSB selected " <<endl;
  *lowerSidebandData=(TH1F*)gDirectory->Get(histoName);
}

void doWrongSignComparison(TTree* mcTree,TTree* dataTree, int leptonId,  int numPions, TH1F** sameChargeMC, TH1F** chargeNeutralMC,  TH1F** sameChargeData, TH1F** chargeNeutralData, int channel)
{

  int channelIdx=channel+1;
  char sameChargeSelection[2000];
  char chargeNeutralSelection[2000];

  //since we put the tagCorr in there, we have to have a string w/o for the data
  char sameChargeSelectionData[2000];
  char chargeNeutralSelectionData[2000];
  char histoName[200];
  char drawCommand[2000];
  char buffer[2000];

  char channelString[500];
  char channelSelectionData[1000];

  getChannelString(channel,channelString);
  getChargeAndStarSelection(channelSelectionData,channel,false,numPions);

  addCorrections(buffer);
  //select sidebands from all (need to redo all components because we select a different range
  if(leptonId!=0)
    {

      if(numPions==0)
	sprintf(sameChargeSelection,"%s tagCorr*" P0STRING " && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==1)
	sprintf(sameChargeSelection,"%s tagCorr*" P1STRING " && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==2)
	sprintf(sameChargeSelection,"%s tagCorr*(" P2STRING "&& bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      
      if(numPions==0)
	sprintf(sameChargeSelection,"%s tagCorr*" P0STRING " && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==1)
        sprintf(chargeNeutralSelection,"%s tagCorr*" P1STRING " && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==2)
        sprintf(chargeNeutralSelection,"%s tagCorr*(" P2STRING "&& bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);

      if(numPions==0)
	sprintf(sameChargeSelectionData," " P0STRING " && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ", upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==1)
	sprintf(sameChargeSelectionData," " P1STRING " && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ", upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==2)
      sprintf(sameChargeSelectionData," (" P2STRING "&& bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ", upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      
      if(numPions==0)
        sprintf(chargeNeutralSelectionData," " P0STRING " && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==1)
        sprintf(chargeNeutralSelectionData," " P1STRING " && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      if(numPions==2)
        sprintf(chargeNeutralSelectionData," (" P2STRING "&& bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
      if(numPions==0)
	sprintf(sameChargeSelection,"%s tagCorr*" P0STRING "  && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==1)
	sprintf(sameChargeSelection,"%s tagCorr*" P1STRING "  && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==2)
	sprintf(sameChargeSelection,"%s tagCorr*(" P2STRING " && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);

      if(numPions==0)
	sprintf(chargeNeutralSelection,"%s tagCorr*" P0STRING "  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==1)
	sprintf(chargeNeutralSelection,"%s tagCorr*" P1STRING "  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==2)
	sprintf(chargeNeutralSelection,"%s tagCorr*(" P2STRING " && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);


      if(numPions==0)
	sprintf(sameChargeSelectionData," " P0STRING "  && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==1)
	sprintf(sameChargeSelectionData," " P1STRING "  && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==2)
	sprintf(sameChargeSelectionData," (" P2STRING " && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);

      if(numPions==0)
	sprintf(chargeNeutralSelectionData," " P0STRING "  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==1)
	sprintf(chargeNeutralSelectionData," " P1STRING "  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      if(numPions==2)
	sprintf(chargeNeutralSelectionData," (" P2STRING " && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
    }
  sprintf(histoName,"sameChargeMC_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_WSMC=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  // sprintf(drawCommand,"mNu2 >> %s",histoName);
   sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  int counts=mcTree->Draw(drawCommand,(char*)sameChargeSelection);
  cout <<"got " << counts <<" counts from mc same charge selected " <<endl;
  *sameChargeMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"neutralChargeMC_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_NCMC=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  //sprintf(drawCommand,"mNu2 >> %s",histoName);
   sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  counts=mcTree->Draw(drawCommand,(char*)chargeNeutralSelection);
  cout <<"got " << counts <<" counts from mc charge neutral selected " <<endl;
  *chargeNeutralMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"sameChargeData_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_SCData=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  //  sprintf(drawCommand,"mNu2 >> %s",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
   sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  counts=dataTree->Draw(drawCommand,(char*)sameChargeSelectionData);
  cout <<"got " << counts <<" counts from data same charge selected " <<endl;
  *sameChargeData=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"neutralChargeData_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_NCData=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  //sprintf(drawCommand,"mNu2 >> %s",histoName);
   sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  counts=dataTree->Draw(drawCommand,(char*)chargeNeutralSelectionData);
  cout <<"got " << counts <<" counts from data charge neutral selected " <<endl;
  *chargeNeutralData=(TH1F*)gDirectory->Get(histoName);

}

void addPoissonNoise(TH1F* h)
{
  for(int i=1;i<=h->GetNbinsX();i++)
    {
      double bc=h->GetBinContent(i);
      //
      double nv=0.0;
      double uncert=0.0;
      if(bc>10)
	{
	  //	  nv=rnd->Gaus(bc,h->GetBinError(i));
	  nv=rnd->Poisson(bc);
	//	  uncert=sqrt(bc);
	}
      else{
	//e.g. if errors were computed with 10000 counts and then scaled by 10, we would have an error of 10 and bin count of 1000
	///---> scale factor 1000/(10*10)=10
	////	float scaleFactor=h->GetBinContent(i)/(h->GetBinError(i)*h->GetBinError(i));
	///	bc=bc*scaleFactor;
	nv=rnd->Poisson(bc);
	//	uncert=bc;
	///	nv=nv/scaleFactor;
      }
      if(nv>=0)
	{
	  uncert=sqrt(nv);
	  h->SetBinContent(i,nv);
	  //  h->SetBinError(i,uncert);
	}
      else
	{
	  h->SetBinContent(i,0);
	  //	  h->SetBinError(i,0.0);
	}
      //      cout <<"new error: "<< h->GetBinError(i) <<endl;;
    }
}

void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents, TH1D* pulls,double signalFraction, int fixThresholdCounts, int numPions)
{
  cout <<"trying stability test.. " << endl;
  TH1F** templates=new TH1F*[numComponents];

  float S=0;
  cout <<"doing the templates .." <<endl;
  char buffer[900];
  int iteration=0;


  for(int i=0;i<numComponents;i++)
    {
      cout <<"looking at template " << i <<endl;
      sprintf(buffer,"%s_Clone_it%d",templatesOrg[i]->GetName(),iteration);
      templates[i]=(TH1F*)templatesOrg[i]->Clone(buffer);
      cout <<"before template " << i <<" has " << templates[i]->Integral() <<" int and " << templates[i]->GetEntries()<<endl;
      if(templates[i]->GetEntries()>0)
      	{
	  //get rid of the sumw2
	  //templates[i]->GetSumw2()->Set(0);
			  //			  			  templates[i]->Scale(5);
	  //	  //make the error essentially zero

	  ///-->check if this has an impact on the pulls...
	  //	  	  if(i==6)
	  //	  	    templates[i]->Scale(10000);

	  addPoissonNoise(templates[i]);
	  //	  templates[i]->Scale(1000);
      	}
            cout <<"template " << i <<" has " << templates[i]->Integral() <<" int and " << templates[i]->GetEntries()<<endl;
    }
  cout <<"SIG idx: " << SIG_IDX <<endl;
  float SData=templates[SIG_IDX]->Integral();
  sprintf(buffer,"%s_Clone_it%d",dataOrg->GetName(),iteration);
  TH1F* data=(TH1F*)dataOrg->Clone(buffer);
  cout <<"doing the data " << endl;
  data->GetSumw2()->Set(0);
  //  data->Scale(0.2);
  addPoissonNoise(data);
  float newSignalFraction=SData/(data->Integral());
  cout <<"data has " <<SData <<" counts"<<endl;
  double fitVal=0;
  double fitErr=0;
  double* allFitVals=new double[numComponents];
  double* allFitErrs=new double[numComponents];
  int numEffective=0;
  if(SData>0)
    {
      TFractionFitter* _fit, *_fit2;
      vector<int> _effectiveComponentIndices;
      Int_t _status=0;
      cout <<"data integral before fit: "<< data->Integral()<< " entries: " << data->GetEntries()<<endl;
      TH1F* result;
      TH1F* mcPredictions[100];
            S=getFitSignal_RooFit(data,templates,numComponents, result,mcPredictions,fitVal, fitErr, fixThresholdCounts,allFitVals,allFitErrs, numEffective,_effectiveComponentIndices, _status, numPions);  
	    //S=getFitSignal(data,templates,numComponents, result,mcPredictions,fitVal, fitErr, fixThresholdCounts,allFitVals,allFitErrs, numEffective,_effectiveComponentIndices, _status, numPions);  
      cout <<"data integral: "<< data->Integral()<< " entries: " << data->GetEntries()<<endl;
      cout <<"fitVal: " << fitVal <<" fitErr: "<< fitErr <<endl;
      cout <<"S is : " << S << " fraction times data: " << data->Integral()*fitVal <<" or " << data->GetEntries()*fitVal<<endl;
      if(S>0&& fitErr>0.0)
	{
	  cout <<"Signal significance: " << fitVal/fitErr <<endl;
	  sigSignificance[glChannelIdx]->Fill(fitVal/fitErr);
	  cout <<"signalFraction: " << signalFraction <<" fitVal: "<< fitVal << " diff: "<< signalFraction-fitVal<<endl;
	  cout <<"pull: " << (fitVal-signalFraction)/fitErr <<endl;
	  cout <<"or " << (fitVal-newSignalFraction)/fitErr <<endl;
	  pulls->Fill((fitVal-signalFraction)/fitErr);
	  //pulls->Fill((fitVal-newSignalFraction)/fitErr);
	}
      float errData=1.0;
      float err=1.0;
      if(SData>0)
	errData=sqrt(SData);
      if(S>0)
	err=sqrt(S);
    }
  if(S > 0 && SData> 0)
    //    cout << "S: " << S << " Sdata: " << SData << " diff " << S/SData <<" pull: " << (S/SData)/ sqrt( (errData/SData)*(errData/SData)+(err/S)*(err/S))*S/SData <<endl;
    //    cout << "S: " << S << " Sdata: " << SData << " diff " << S-SData <<" pull: " << (S-SData)/()<<endl;
    {

    }
  else
    {
      cout <<" S or SD =<0"<<endl;
    }

  for(int i=0;i<numComponents;i++)
    {
      delete templates[i];
    }
  delete templates;
}

float getFitSignal_RooFit( TH1F* data,TH1F** templates,int numTemplates,    TH1F* &result, TH1F** mcPredictions,  double& fitVal, double& fitErr, int fixThresholdCounts, double* allFitVals, double* allFitErrs, int& numEffective,vector<int>& effectiveComponentsIndices, Int_t& status, int numPions)
{

  //  RooDataHist data();
  //need to subtract the templates with small counts 
  float S=0;
  vector<int> countsOfComponents;
  vector<float> integralOfComponents;
  //  vector<int> effectiveComponentsIndices; 
  char buffer[500];
  int numEffectiveComponents=0;
  TObjArray *mc = new TObjArray(numTemplates);        // MC histograms are put in this array
  int signalIndex=-1;
  int otherBBIndex1=-1;
  int otherBBIndex2=-1;
  int continuumIndex1=-1;
  int continuumIndex2=-1;
  //maybe it is better to clone the signal and subtract all the fixed components
  TCanvas c;
  double templateIntegral=0.0;

  pCount++;

  RooDataHist* rooHists[100];
  RooHistPdf* rooHistPdfs[100];
  sprintf(buffer, "rooRealDataHist_numP_%d_%s",numPions,channelStringGlobal);

  RooRealVar rooMnu2("mNu2","mNu2",lowerCut[glChannelIdx],upperCut[glChannelIdx]);
  RooDataHist rooData(buffer,buffer,rooMnu2,RooFit::Import(*data));
  //get the non-zero ones...
  for(int i=0;i<numTemplates;i++)
    {
      //      if(templates[i]->GetEntries()>0)
      if(templates[i]->Integral()>0)
	{
	  //	  templates[i]->Sumw2(0);
	  cout <<"template " << i << " ("<<templates[i]->GetName() <<")  has " << templates[i]->Integral()<<endl;
	  numEffectiveComponents++;
	  sprintf(buffer,"poissonTemplate%d_numPions%d_%s.png",i,numPions,channelStringGlobal);
	  templates[i]->Draw();
	  if(((pCount-1)%100)==0 || pCount< 7)
	    {
	      c.SaveAs(buffer);
	      sprintf(buffer,"poissonTemplate%d_numPions%d_%s.pdf",i,numPions,channelStringGlobal);
	      c.SaveAs(buffer);
	      sprintf(buffer,"poissonTemplate%d_numPions%d_%s.png",i,numPions,channelStringGlobal);
	      c.SaveAs(buffer);
	    }
	  //because we add later
	  int effIndex=effectiveComponentsIndices.size();
	  sprintf(buffer, "rooDataHist_%d_numP_%d_%s",i,numPions,channelStringGlobal);
	  rooHists[effIndex]=new RooDataHist(buffer,buffer,rooMnu2,RooFit::Import(*templates[i]));
	  sprintf(buffer, "rooDataHistPdf_%d_numP_%d_%s",i,numPions,channelStringGlobal);
	  cout <<"creating new pdf at index: "<< effIndex <<endl;
	  rooHistPdfs[effIndex]=new RooHistPdf(buffer,buffer, rooMnu2,*rooHists[effIndex],0);


	  //computation otherwise to complex with template combination and the like
	  if((string(templates[i]->GetName()).find("continuum")!=string::npos) && (string(templates[i]->GetName()).find("BToD_")!=string::npos || string(templates[i]->GetName()).find("B0ToD_")!=string::npos))
	    {
	      continuumIndex1=effectiveComponentsIndices.size();
	    }
	  if((string(templates[i]->GetName()).find("continuum")!=string::npos) && (string(templates[i]->GetName()).find("BToDStar_")!=string::npos || string(templates[i]->GetName()).find("B0ToDStar_")!=string::npos))
	    {
	      continuumIndex2=effectiveComponentsIndices.size();
	    }
	  if((string(templates[i]->GetName()).find("OtherBB")!=string::npos) && (string(templates[i]->GetName()).find("BToD_")!=string::npos || string(templates[i]->GetName()).find("B0ToD_")!=string::npos))
	    {
	      otherBBIndex1=effectiveComponentsIndices.size();
	    }
	  if((string(templates[i]->GetName()).find("OtherBB")!=string::npos) && (string(templates[i]->GetName()).find("BToDStar")!=string::npos|| string(templates[i]->GetName()).find("B0ToDStar")!=string::npos))
	    {
	      otherBBIndex2=effectiveComponentsIndices.size();
	    }
	  if(i==SIG_IDX)
	    {
	      signalIndex=effectiveComponentsIndices.size();
	      cout <<"set signalIndex to " << signalIndex<<endl;
	    }
	  effectiveComponentsIndices.push_back(i);

	  //the same as the test effectiveComponentIndices[i]==2 from the regular fit
	  //	  countsOfComponents.push_back(templates[i]->Integral());
	  double tempIntegral=templates[i]->Integral();
	  cout <<"template " << i << " has " << templates[i]->Integral() <<" and " << templates[i]->GetEntries() << " counts " << endl;
	  integralOfComponents.push_back(tempIntegral);
	  //	  if((i==SIG_IDX && numPions > 0 )|| (i==iDLNu && numPions==0))
	  if(i==SIG_IDX)
	    {
	      templateIntegral+=tempIntegral;
	    }
	  else{
	    templateIntegral+=tempIntegral;
	  }
	}
    }
  cout <<"other bb indices: " << otherBBIndex1 <<" and " << otherBBIndex2 <<endl;
  cout <<"we have " << numEffectiveComponents <<" non-zero components " <<endl;
  //  data->Sumw2(0);
  double dataIntegral=data->Integral();
  cout <<"done with integral: "<< templateIntegral <<", data integral: "<< dataIntegral <<endl;
// to scale roughly to the 5 (or 4) streams used. For the poisson noise add the two integrals should be roughly the same
//this was dataIntegral<3*templateIntegral before... probably wrong
  if(3*dataIntegral<templateIntegral)
    {
      cout <<" rescaling data integral by factor 5 " << endl;
      dataIntegral*=5; 
    }

  data->Draw();
  cout <<"done with data draw.." << endl;
  if(((pCount-1)%100)==0 || pCount < 5)
    {
      sprintf(buffer,"poissonData_%s.png",channelStringGlobal);
      c.SaveAs(buffer);
      sprintf(buffer,"poissonData_%s.png",channelStringGlobal);
      c.SaveAs(buffer);
    }

   //   fit->GetFitter()->SetPrecision(0.001);
  double sumOfFractions=0.0;
  double sumOfIntegrals=0.0;
  double signalFraction=0.0;
  vector<int> fixedComponents;
  vector<double> backgroundFractions;
  double totalBGFraction=0.0;
  for(int i=0;i<numEffectiveComponents;i++)
    {

      //  int fixThresholdCounts=2000;
#ifndef PARTIAL_BOX 
      if((integralOfComponents[i]>fixThresholdCounts|| i==signalIndex) && !(i==otherBBIndex1 || i==otherBBIndex2) && !(i==continuumIndex1|| i==continuumIndex2))
#else
	//the low counts in the partial_box case doesn't allow to differentiate between feed-down and other bb
		if((integralOfComponents[i]>fixThresholdCounts|| i==signalIndex) && !(i==otherBBIndex1 || i==otherBBIndex2) && !(i==continuumIndex1|| i==continuumIndex2))
	//	if((integralOfComponents[i]>fixThresholdCounts|| i==signalIndex) )//  && !(i==otherBBIndex1 || i==otherBBIndex2) && !(i==continuumIndex1|| i==continuumIndex2))
#endif 

	{
	  //anyways good to give decent start value. The only caveat is that the signal is most likely less than what we think in MC. So the signal fraction is proably a bit overestimated and 
	  //the other fractions underestimated. But in principel we can calculate the best start value by taking the pdg value for the signal (i.e. half es much as MC pred)
	  double iThFraction=integralOfComponents[i]/templateIntegral;
	  cout <<"iTh: " << iThFraction <<endl;
	  //aim for roughly constant counts for these backgrounds. So if we see much less data than expected the fraction should rise
	  cout <<"component integral : " << integralOfComponents[i] << ", templateIntegral: " << templateIntegral <<" dataIntegral : " << dataIntegral <<endl;
	  //	  iThFraction*=(templateIntegral/dataIntegral);
	  cout <<"now: "<< iThFraction <<endl;
	  sumOfFractions+=iThFraction;
	  sumOfIntegrals+=integralOfComponents[i];
	  if(i==signalIndex)
	    {
	      signalFraction=iThFraction;
	      //to have some variation...
	      //	      iThFraction-=0.05;
	      //two the factor two mojo...

	      //	      iThFraction=integralOfComponents[i]/(2*templateIntegral);
	      //	      iThFraction=integralOfComponents[i]/(templateIntegral);
	      //	      iThFraction/=2;
	      cout << " setting signal fraction to  " << iThFraction <<endl;
	      
	    }
	  else{
	    backgroundFractions.push_back(iThFraction);
	    totalBGFraction+=iThFraction;
	    cout << " setting non signal fraction " << i <<" ("<<templates[effectiveComponentsIndices[i]]->GetName() <<")  to  " << iThFraction <<endl;
	  }
	  sprintf(buffer,"para%d",i);
	  float upperLimit=1.0;
	  float lowerLimit=0.0;
	  upperLimit=2*iThFraction;
	  if(iThFraction>0.1)
	    {
	      lowerLimit=0.5*iThFraction;
	    }
	  float paraUncert=0.1*iThFraction;
	  paraUncert=0.1; //better than above
	  //	  paraUncert=0.;
	  ////->	  fit->GetFitter()->SetParameter(i,buffer,iThFraction,paraUncert,lowerLimit,upperLimit);

	  //does the set parameter do already the same as constrain?
	  //	  fit->Constrain(i,0.0,1.0);               // constrain fraction i to be between 0 and 1
	}
      else
	{
	  fixedComponents.push_back(i);
	  double iThFraction=integralOfComponents[i]/templateIntegral;
	  backgroundFractions.push_back(iThFraction);
	  totalBGFraction+=iThFraction;
	  sumOfFractions+=iThFraction;
	  sumOfIntegrals+=integralOfComponents[i];
	  cout <<"component integral : " << integralOfComponents[i] << ", templateIntegral: " << templateIntegral <<" dataIntegral : " << dataIntegral <<endl;
	  cout <<"fixing parameter " << i <<" ("<<templates[effectiveComponentsIndices[i]]->GetName() <<") to " << iThFraction <<endl;
	  sprintf(buffer,"para%d",i);

	}
    }
  //now construct everything we need for the roofit...:
  ////----
  int bgCounter=-1;
  RooArgList bgShapes;
  //minus the signal
  RooArgList bgFracList; 
  cout<<"constructing roo bg shapes, total bgFraction: "<< totalBGFraction <<endl;
  RooRealVar* rvArr[100];
  for(int i=0;i<numEffectiveComponents;i++)
    {
      bgCounter++;
      if(i!=signalIndex)
	{
	  bgShapes.add(*rooHistPdfs[i]);
	  sprintf(buffer, "roofracVar_%d_numP_%d_%s",i,numPions,channelStringGlobal);
	  RooRealVar* rv=new RooRealVar(buffer,buffer,backgroundFractions[bgCounter]/totalBGFraction,0,100000);
	  //so we can delete more easily
	  rvArr[i]=rv;
	  //if constant
	  if(find(fixedComponents.begin(),fixedComponents.end(),i)!=fixedComponents.end())
	    {
	      cout <<"fixing component " << i  << endl;
	      rv->setConstant(true);
	    }
	  else
	    {
	      cout <<"component " << i << " is left unfixed " << endl;
	    }
	  bgFracList.add(*rv);
	}
    }
  sprintf(buffer, "bgPdf_numP_%d_%s",numPions,channelStringGlobal);
  RooAddPdf bgPdf(buffer,buffer,bgShapes,bgFracList);

  RooArgList sigBgFractions;
  sprintf(buffer, "totalBgFraction_numP_%d_%s",numPions,channelStringGlobal);
  RooRealVar rooBgFrac(buffer,buffer,totalBGFraction,0,100000);
  sprintf(buffer, "totalSigFraction_numP_%d_%s",numPions,channelStringGlobal);
  RooRealVar rooSigFrac(buffer,buffer,totalBGFraction,0,100000);
  sigBgFractions.add(rooBgFrac);
  sigBgFractions.add(rooSigFrac);
  RooArgList sigBgPdfs;
  sigBgPdfs.add(bgPdf);
  sigBgPdfs.add(*rooHistPdfs[signalIndex]);
  RooAddPdf totalPdf(buffer,buffer,sigBgPdfs,sigBgFractions);

  RooPlot* bFrame1=rooMnu2.frame();
  RooPlot* bFrame2=rooMnu2.frame();
  RooPlot* bFrame3=rooMnu2.frame();
  RooPlot* bFrame4=rooMnu2.frame();
  totalPdf.plotOn(bFrame1);
  bFrame1->Draw();
  sprintf(buffer, "before_totalPdf_numP_%d_%s.png",numPions,channelStringGlobal);
  c.SaveAs(buffer);
  bgPdf.plotOn(bFrame2);
  bFrame2->Draw();
  sprintf(buffer, "before_bgPdf_numP_%d_%s.png",numPions,channelStringGlobal);
  c.SaveAs(buffer);
  sprintf(buffer, "before_sigPdf_numP_%d_%s.png",numPions,channelStringGlobal);
  rooHistPdfs[signalIndex]->plotOn(bFrame3);
  bFrame3->Draw();
  c.SaveAs(buffer);
  sprintf(buffer, "before_rooData_numP_%d_%s.png",numPions,channelStringGlobal);
  rooData.plotOn(bFrame4);
  bFrame4->Draw();
  c.SaveAs(buffer);

  totalPdf.fitTo(rooData,RooFit::Extended());
  //get result, produce result and mcPred plots
  RooPlot* mFrame1=rooMnu2.frame();
  RooPlot* mFrame2=rooMnu2.frame();
  RooPlot* mFrame3=rooMnu2.frame();
  RooPlot* mFrame4=rooMnu2.frame();
  totalPdf.plotOn(mFrame1);
  mFrame1->Draw();
  sprintf(buffer, "totalPdf_numP_%d_%s.png",numPions,channelStringGlobal);
  c.SaveAs(buffer);
  bgPdf.plotOn(mFrame2);
  mFrame2->Draw();
  sprintf(buffer, "bgPdf_numP_%d_%s.png",numPions,channelStringGlobal);
  c.SaveAs(buffer);
  sprintf(buffer, "sigPdf_numP_%d_%s.png",numPions,channelStringGlobal);
  rooHistPdfs[signalIndex]->plotOn(mFrame3);
  mFrame3->Draw();
  c.SaveAs(buffer);

  sprintf(buffer, "rooData_numP_%d_%s.png",numPions,channelStringGlobal);
  rooData.plotOn(mFrame4);
  mFrame4->Draw();
  c.SaveAs(buffer);
  cout <<"fitted sig fraction: "<<   rooSigFrac.getValV() <<" +- " <<rooSigFrac.getError()<<endl;
  cout <<"fitted bg fraction: "<<   rooBgFrac.getValV() <<" +- " <<rooBgFrac.getError()<<endl;
  fitVal=rooSigFrac.getValV();
  fitErr=rooSigFrac.getError();

  ////save output histograms to see that we have the right stuff...




  sprintf(buffer, "result_numP_%d_%s",numPions,channelStringGlobal);

  result=(TH1F*)templates[signalIndex]->Clone(buffer);
  for(int i=0;i<numEffectiveComponents;i++)
    {
      sprintf(buffer, "mcPred_%d_numP_%d_%s",i,numPions,channelStringGlobal);
      mcPredictions[i]=(TH1F*)templates[effectiveComponentsIndices[i]]->Clone(buffer);
      if(i!=signalIndex)
	{
	  //but this is relative to the bgFraction
	  allFitVals[i]=rvArr[i]->getValV()*rooBgFrac.getValV();
	  cout <<" component: " << i <<" fit val: "<< allFitVals[i] << " uncert: "<< rvArr[i]->getError()<<endl;
	  if(rvArr[i]->getValV()>0)
	    {
	      //relative error
	    allFitErrs[i]=sqrt(rvArr[i]->getError()*rvArr[i]->getError()/(rvArr[i]->getValV()*rvArr[i]->getValV())+rooBgFrac.getError()*rooBgFrac.getError()/(rooBgFrac.getValV()*rooBgFrac.getValV()));
	    //total error
	    allFitErrs[i]*=allFitVals[i];
	    }
	  else
	    allFitErrs[i]=0.0;
	}
      else
	{
	  allFitVals[i]=rooSigFrac.getValV();
	  allFitErrs[i]=rooSigFrac.getError();
	}
    }


  //delete everyting
  cout <<" deleting roo objects " <<endl;
  for(int i=0;i<numEffectiveComponents;i++)
    {
      if(i!=signalIndex)
	{
	  delete rvArr[i];
	}

      delete rooHistPdfs[i];
      delete rooHists[i];
    }
  cout <<" done " <<endl;

    //    RooHistPdf pdf("mpdf","mpdf",);

}

//other fit method where we subtract all the fixed fractions from the data and just fit the remaining
float getFitSignal_DiffMethod( TH1F* data,TH1F** templates,int numTemplates,     TH1F* result, double& fitVal, double& fitErr, int fixThresholdCounts, double* allFitVals, int& numEffective,vector<int>& effectiveComponentsIndices, Int_t& status)
{
  char buffer[1000];
  TObjArray *mc2 = new TObjArray(numTemplates);        // MC histograms are put in this array
  int signalIndex2=0;
  int numEffectiveComponents2=0;
  vector<int> effectiveComponentsIndices2; 
  double fitVal2,fitErr2;
  float S=0.0;
  int numEffective2=0;
  double* allFitVals2=new double[20];

  double totalIntegral2=0.0;
  sprintf(buffer,"%s_FitClone",data->GetName());
  TH1F* data2=(TH1F*) data->Clone();
  cout <<"cloned data " << endl;
  for(int i=0;i<numTemplates;i++)
    {
      if((templates[i]->Integral()>fixThresholdCounts) || (i==SIG_IDX))
	{
	  cout <<"adding template " << i <<" with " << templates[i]->Integral() <<" counts " << endl;
	  numEffectiveComponents2++;
	  mc2->Add(templates[i]);
	  totalIntegral2+=templates[i]->Integral();
	  if(i==SIG_IDX)
	    {
	      signalIndex2=effectiveComponentsIndices2.size();
	    }
	  effectiveComponentsIndices2.push_back(i);
	}
      else
	{
	  data2->Add(templates[i],-1);
	}
    }
  TFractionFitter* fit2 = new TFractionFitter(data2, mc2); // initialise
  cout <<"total integral after subtraction: " << totalIntegral2 <<endl;
  for(int i=0;i<numEffectiveComponents2;i++)
    {
      cout <<"calculating fraction of eff comp: "<< effectiveComponentsIndices2[i] <<" has : "<< templates[effectiveComponentsIndices2[i]]->Integral() << " counts " << endl;
      double fraction=templates[effectiveComponentsIndices2[i]]->Integral()/totalIntegral2;
      sprintf(buffer,"para_2_%d",i);
      cout <<"setting " << buffer <<" to " << fraction <<endl;
      ////->      fit2->GetFitter()->SetParameter(i,buffer,fraction,0.1,0.0,1.0);
    }
  Int_t status2 = fit2->Fit();               // perform the fit
  std::cout << "fit status2: " << status2 << std::endl;

  if (status2 == 0) 
    {
      numEffective2=numEffectiveComponents2;
    //we want to flip the signal (index =2 ) so that it is later
    //	if(effectiveComponentsIndices[i]!=2)
      TH1F* mcComp=(TH1F*) fit2->GetMCPrediction(signalIndex2);
      S=mcComp->Integral();
      fit2->GetResult(signalIndex2,fitVal2,fitErr2);
      cout <<"signal fit2: " << fitVal2 <<" +- " << fitErr2 <<endl;
      for(int i=0;i<numEffectiveComponents2;i++)
	{
	  cout <<"effective component2 " << i<< " has " << fit2->GetMCPrediction(i)->Integral() <<" integral and  " << fit2->GetMCPrediction(i)->GetEntries() <<" counts, fit: ";
	  double temp=0.0;
	  double tempErr=0.0;
	  fit2->GetResult(i, temp, tempErr);
	  cout << temp<<endl;
	  allFitVals2[i]=temp;
	  //	  allFitVals[i]=temp;
	}
    }
  return S;
}




float getFitSignal(TH1F* data,TH1F** templates,int numTemplates,     TH1F* &result, TH1F** mcPredictions, double& fitVal, double& fitErr, int fixThresholdCounts, double* allFitVals,double* allFitErrs ,int& numEffective,vector<int>& effectiveComponentsIndices, Int_t& status, int numPions)
{

  TDirectory* oldDir=gDirectory;
  char buffer[500];  
  char channelBuffer[500];
  getChannelString(glChannelIdx-1,channelBuffer);
  sprintf(buffer,"fitInput_%s.root",channelBuffer);
  TFile lFile(buffer,"RECREATE");
  for(int t=0;t<numTemplates;t++)
    {
      templates[t]->SetDirectory((TDirectory*) &lFile);
      templates[t]->Write();
    }
  data->SetDirectory((TDirectory*) &lFile);
  data->Write();
  lFile.Write();
  data->SetDirectory(oldDir);
  for(int t=0;t<numTemplates;t++)
    {
      templates[t]->SetDirectory(oldDir);
    }
  lFile.Close();

  cout <<"getFitSignal: expecting " << numTemplates << " templates " << endl;
  float S=0;
  vector<int> countsOfComponents;
  vector<float> integralOfComponents;
  //  vector<int> effectiveComponentsIndices; 
  int numEffectiveComponents=0;
  TObjArray *mc = new TObjArray(numTemplates);        // MC histograms are put in this array
  int signalIndex=-1;
  int otherBBIndex1=-1;
  int otherBBIndex2=-1;
  int continuumIndex1=-1;
  int continuumIndex2=-1;
  //maybe it is better to clone the signal and subtract all the fixed components
  TCanvas c;
  double templateIntegral=0.0;
  oldDir->cd();

  pCount++;
  //get the non-zero ones...
  for(int i=0;i<numTemplates;i++)
    {
      //      if(templates[i]->GetEntries()>0)
      if(templates[i]->Integral()>0)
	{
	  //	  templates[i]->Sumw2(0);
	  cout <<"template " << i << " ("<<templates[i]->GetName() <<")  has " << templates[i]->Integral()<<endl;
	  numEffectiveComponents++;
	  mc->Add(templates[i]);
	  sprintf(buffer,"poissonTemplate%d_numPions%d_%s.png",i,numPions,channelStringGlobal);
	  templates[i]->Draw();
	  if(((pCount-1)%100)==0 || pCount< 7)
	    {
	      c.SaveAs(buffer);
	      sprintf(buffer,"poissonTemplate%d_numPions%d_%s.pdf",i,numPions,channelStringGlobal);
	      c.SaveAs(buffer);
	      sprintf(buffer,"poissonTemplate%d_numPions%d_%s.png",i,numPions,channelStringGlobal);
	      c.SaveAs(buffer);
	    }
	  //computation otherwise to complex with template combination and the like
	  if((string(templates[i]->GetName()).find("continuum")!=string::npos) && (string(templates[i]->GetName()).find("BToD_")!=string::npos || string(templates[i]->GetName()).find("B0ToD_")!=string::npos))
	    {
	      continuumIndex1=effectiveComponentsIndices.size();
	    }
	  if((string(templates[i]->GetName()).find("continuum")!=string::npos) && (string(templates[i]->GetName()).find("BToDStar_")!=string::npos || string(templates[i]->GetName()).find("B0ToDStar_")!=string::npos))
	    {
	      continuumIndex2=effectiveComponentsIndices.size();
	    }
	  if((string(templates[i]->GetName()).find("OtherBB")!=string::npos) && (string(templates[i]->GetName()).find("BToD_")!=string::npos || string(templates[i]->GetName()).find("B0ToD_")!=string::npos))
	    {
	      otherBBIndex1=effectiveComponentsIndices.size();
	    }
	  if((string(templates[i]->GetName()).find("OtherBB")!=string::npos) && (string(templates[i]->GetName()).find("BToDStar")!=string::npos|| string(templates[i]->GetName()).find("B0ToDStar")!=string::npos))
	    {
	      otherBBIndex2=effectiveComponentsIndices.size();
	    }
	  if(i==SIG_IDX)
	    {
	      signalIndex=effectiveComponentsIndices.size();
	      cout <<"set signalIndex to " << signalIndex<<endl;
	    }
	  effectiveComponentsIndices.push_back(i);

	  //the same as the test effectiveComponentIndices[i]==2 from the regular fit
	  //	  countsOfComponents.push_back(templates[i]->Integral());
	  double tempIntegral=templates[i]->Integral();
	  integralOfComponents.push_back(tempIntegral);
	  //	  if((i==SIG_IDX && numPions > 0 )|| (i==iDLNu && numPions==0))
	  if(i==SIG_IDX)
	    {
	      templateIntegral+=tempIntegral;
	    }
	  else{
	    templateIntegral+=tempIntegral;
	  }
	}
    }
  cout <<"other bb indices: " << otherBBIndex1 <<" and " << otherBBIndex2 <<endl;
  cout <<"we have " << numEffectiveComponents <<" non-zero components " <<endl;
  //  data->Sumw2(0);
  double dataIntegral=data->Integral();
  cout <<"done with integral: "<< templateIntegral <<", data integral: "<< dataIntegral <<endl;
// to scale roughly to the 5 (or 4) streams used. For the poisson noise add the two integrals should be roughly the same
//this was dataIntegral<3*templateIntegral before... probably wrong
  if(3*dataIntegral<templateIntegral)
    {
      cout <<" rescaling data integral by factor 5 " << endl;
      dataIntegral*=5; 
    }

  data->Draw();
  cout <<"done with data draw.." << endl;
  if(((pCount-1)%100)==0 || pCount < 5)
    {
      sprintf(buffer,"poissonData_%s.png",channelStringGlobal);
      c.SaveAs(buffer);
      sprintf(buffer,"poissonData_%s.png",channelStringGlobal);
      c.SaveAs(buffer);
    }
  cout <<"overall integral of templates: " << templateIntegral <<" of data: "<< data->Integral()<< " or " << dataIntegral<<endl;
  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise

   //   fit->GetFitter()->SetPrecision(0.001);
  double sumOfFractions=0.0;
  double sumOfIntegrals=0.0;
  for(int i=0;i<numEffectiveComponents;i++)
    {
      cout <<" i: " << i << " otherBB1: "<< otherBBIndex1 << " BB2: "<< otherBBIndex2 <<endl;
      cout <<"if statement: ";
      if((integralOfComponents[i]>fixThresholdCounts|| i==signalIndex) && !(i==otherBBIndex1 || i==otherBBIndex2) && !(i==continuumIndex1|| i==continuumIndex2))
	cout <<true;
      else
	cout <<false;
      cout <<endl;
      //use the leptonId=0
      //  int fixThresholdCounts=2000;
#ifndef PARTIAL_BOX 
      if((integralOfComponents[i]>fixThresholdCounts|| i==signalIndex) && !(i==otherBBIndex1 || i==otherBBIndex2) && !(i==continuumIndex1|| i==continuumIndex2))
#else
	//the low counts in the partial_box case doesn't allow to differentiate between feed-down and other bb
		if((integralOfComponents[i]>fixThresholdCounts|| i==signalIndex) && !(i==otherBBIndex1 || i==otherBBIndex2) && !(i==continuumIndex1|| i==continuumIndex2))
	//	if((integralOfComponents[i]>fixThresholdCounts|| i==signalIndex) )//  && !(i==otherBBIndex1 || i==otherBBIndex2) && !(i==continuumIndex1|| i==continuumIndex2))
#endif 

	{
	  //anyways good to give decent start value. The only caveat is that the signal is most likely less than what we think in MC. So the signal fraction is proably a bit overestimated and 
	  //the other fractions underestimated. But in principel we can calculate the best start value by taking the pdg value for the signal (i.e. half es much as MC pred)
	  double iThFraction=integralOfComponents[i]/templateIntegral;
	  cout <<"iTh: " << iThFraction <<endl;
	  //aim for roughly constat counts for these backgrounds. So if we see much less data than expected the fraction should rise
	  cout <<"component integral : " << integralOfComponents[i] << ", templateIntegral: " << templateIntegral <<" dataIntegral : " << dataIntegral <<endl;
	  //	  iThFraction*=(templateIntegral/dataIntegral);
	  cout <<"now: "<< iThFraction <<endl;
	  sumOfFractions+=iThFraction;
	  sumOfIntegrals+=integralOfComponents[i];
	  if(i==signalIndex)
	    {
	      //to have some variation...
	      //	      iThFraction-=0.05;
	      //two the factor two mojo...

	      //	      iThFraction=integralOfComponents[i]/(2*templateIntegral);
	      //	      iThFraction=integralOfComponents[i]/(templateIntegral);
	      //	      iThFraction/=2;
	      cout << " setting signal fraction to  " << iThFraction <<endl;
	      
	    }
	  else{
	    cout << " setting non signal fraction " << i <<" ("<<templates[effectiveComponentsIndices[i]]->GetName() <<")  to  " << iThFraction <<endl;
	  }
	  sprintf(buffer,"para%d",i);
	  float upperLimit=1.0;
	  float lowerLimit=0.0;
	  upperLimit=2*iThFraction;
	  if(iThFraction>0.1)
	    {
	      lowerLimit=0.5*iThFraction;
	    }
	  float paraUncert=0.1*iThFraction;
	  paraUncert=0.1; //better than above
	  //	  paraUncert=0.;
	  ////->	  fit->GetFitter()->SetParameter(i,buffer,iThFraction,paraUncert,lowerLimit,upperLimit);
	  fit->GetFitter()->Config().ParSettings(i).SetValue(iThFraction);
	  fit->GetFitter()->Config().ParSettings(i).SetLimits(lowerLimit,upperLimit);
	  //does the set parameter do already the same as constrain?
	  //	  fit->Constrain(i,0.0,1.0);               // constrain fraction i to be between 0 and 1
	}
      else
	{
	  double iThFraction=integralOfComponents[i]/templateIntegral;
	  sumOfFractions+=iThFraction;
	  sumOfIntegrals+=integralOfComponents[i];
	  cout <<"component integral : " << integralOfComponents[i] << ", templateIntegral: " << templateIntegral <<" dataIntegral : " << dataIntegral <<endl;
	  cout <<"fixing parameter " << i <<" ("<<templates[effectiveComponentsIndices[i]]->GetName() <<") to " << iThFraction <<endl;
	  sprintf(buffer,"para%d",i);
	  //	  fit->GetFitter()->SetParameter(i,buffer,iThFraction,0.1,0.0,1.0);
	  ///-->	  fit->GetFitter()->SetParameter(i,buffer,iThFraction,0.0,0.0,1.0);
	  fit->GetFitter()->Config().ParSettings(i).SetValue(iThFraction);
	  fit->GetFitter()->Config().ParSettings(i).Set(buffer,iThFraction);
	  fit->GetFitter()->Config().ParSettings(i).SetLimits(iThFraction,iThFraction);
	  //	  fit->GetFitter()->FixParameter(i);
	  //instead of fix, constrain with same upper and lower bound, hopefully that is the same as fixing it...
	  //	  fit->Constrain(i,iThFraction,iThFraction);
	  //	  fit->GetFitter()->Constrain(i);
	}
    }
  cout <<" sum of fractions of templates: " << sumOfFractions <<" sum of integrals: "<< sumOfIntegrals<<endl;
  status = fit->Fit();               // perform the fit
  std::cout << "fit status: " << status << std::endl;
                       // check on fit status
  if (status == 0) 
    {
     result = (TH1F*) fit->GetPlot();
      numEffective=numEffectiveComponents;
    //we want to flip the signal (index =2 ) so that it is later
    //	if(effectiveComponentsIndices[i]!=2)
      cout <<"getting MC pred " << signalIndex <<endl;
      TH1F* mcComp=(TH1F*) fit->GetMCPrediction(signalIndex);
      S=mcComp->Integral();
      fit->GetResult(signalIndex,fitVal,fitErr);
      double sumFractions=0.0;
      for(int i=0;i<numEffectiveComponents;i++)
	{
	  mcPredictions[i]=(TH1F*) fit->GetMCPrediction(i);
	  cout <<"effective component " << i<< " has " << fit->GetMCPrediction(i)->Integral() <<" integral and  " << fit->GetMCPrediction(i)->GetEntries() <<" counts, fit: ";
	  double temp=0.0;
	  double tempErr=0.0;
	  fit->GetResult(i, temp, tempErr);
	  cout << temp<<endl;
	  sumFractions+=temp;
	  allFitVals[i]=temp;
	  allFitErrs[i]=tempErr;
	}
      cout <<" sum of all fractions of first fit: "<< sumFractions <<endl;
    }


  return S;
} /// end of getFitSignal
void fillTemplates(TH1F** templates,TH1F** summedComponents,char* templateLegendNames,int numPions)
{


};
#endif

