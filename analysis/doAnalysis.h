#ifndef DO_ANALYSIS_H
#define DO_ANALYSIS_H

//#define DPiPi_SEARCH
int maxDataTreeSize;
#define partialBoxFraction 0.15

//#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.0 && mDnPi < 3.5  "
////#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.0 && mDnPi < 3.5  &&mDnPi > 2.05"
#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.04 && deltaETag<0.04 && logProb > -3.0 && mDnPi < 3.5  &&mDnPi > 2.05"


#define P2STRING_SinglePion "bestD==1 && pi2Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb  > -2.6  && mDnPi < 3.0 && !(dType==1 && dDecay==2) && !(dType==1 && dDecay==3) && !(dType==2) && !overlapEvent" 



#define P2STRING_Search "bestD==1 && pi2Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb  > -2.6  && mDnPi < 3.0 && ((bestBCharge==0 && systemCharge==0) || (bestBCharge== -leptonCharge) ) && mDnPi>1.0  && mDnPi < 3.0&&  hypDMass1 > 2.03 && hypDMass2 > 2.03 && U < 2.0 && U > -2.0"  


//3,5,7 are the D pi pi 

#define SIG_IDX_D_PI 2
#define SIG_IDX_D_PI_PI 3

#ifdef DPiPi_SEARCH
#define P2STRING P2STRING_Search
#define SIG_IDX SIG_IDX_D_PI_PI
#else
#define P2STRING P2STRING_SinglePion
#define SIG_IDX SIG_IDX_D_PI
#endif


//#define P2STRING "recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb  > -3.0  && mDnPi < 3.0 "
//&& !(dType==0 && dDecay==1) && !(dType==1 && dDecay==3) && !(dType==2 && dDecay==1) 


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


void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel=false);
TRandom3* rnd;
int numBinsSB=20;
int numBinsSBTop=40;
int numBinsWS=20;
char* allLegendNames[10];

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
int numBins[]={30,30,40,40,40}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
#endif
//int numBins[]={200,2000,2000,200,200}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={50,50,5,200,5}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={70,100,70,70,70}; 
///_--> quite good int numBins[]={70,80,70,70,70};
///bias for BToD goes away with range -0.5 to 2.0 and 80 bins

///--->B0ToD seems to be fine with range of -0.5 to 1.5 with 70 bins...
////--->optimalfloat lowerCut[]={-0.65,-0.6,-0.65,-0.55,-0.25};
//float lowerCut[]={-0.65,-0.6,-0.65,-0.55,-0.25};
float lowerCut[]={-0.5,-0.5,-0.5,-0.5,-0.5};
////////float upperCut[]={2.0,2.0,1.5,1.5,1.5};
float upperCut[]={0.5,0.5,0.5,0.5,0.5};

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
float getFitSignal(TH1F* data,TH1F** templates,int numTemplates, double& fitVal, double& fitErr , int fixThresholdCounts, double* allFitVals, int &numEffective);
void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents, TH1D* pulls,double signalFraction, int fixThresholdCounts);
void fillTemplates(TH1F** templates,TH1F** summedComponents,char* templateLegendNames,int numPions);
void addPoissonNoise(TH1F* h);

TColor* glColorTable[10];
void addCorrections(char* buffer)
{
  stringstream str;
  //D_DecayCorr*B_DecayCorr*PIDCorrection*CrossSectionLumiCorrection

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
    default:
      sprintf(buffer,"_wrongChannel_");
      break;
    }


}

void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel)
{

  if(channel<0)
    {
      //need somestatement here since the string gets appended via &&
      if(wrongChannel)
	sprintf(chargeAndStarSelection," 1!=1");
      else
	sprintf(chargeAndStarSelection," 1==1");
      return;
    }

      if(wrongChannel)
	{
	  if(numPions==1)
	    {
	      switch(channel)
		{
		  //B to D
		case 0:
		  //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
		  sprintf(chargeAndStarSelection,"  (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge ||  dType==2 || abs(mcDCharge)!=1  ||  0!=mcIsDStar )");
		  break;
		  //B to DStar
		case 1:
		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge || dType!=2 || abs(mcDCharge)!=1  || 0==mcIsDStar)");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge || dType==2 || abs(mcDCharge)!=0  || mcIsDStar!=0)");
		  break;
	      //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge || dType!=2 || abs(mcDCharge)!=0  ||  0==mcIsDStar)");
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
		  sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge==mcDCharge || dType==2 || abs(mcDCharge)!=1  || mcIsDStar==0)");
		  break;
		  //B0to D0
		case 2:
		  sprintf(chargeAndStarSelection," 1==2 ");

		  break;
	      //B0to DStar0
		case 3:
		  sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge==mcDCharge || dType==2 ||abs(mcDCharge)!=0  || 0==mcIsDStar)");
		  break;
		default:
		  cout <<" wrong channel!!!  " << channel <<endl;
	      break;
		}

	    }


	  return;
	}
      else
	{
  if(numPions==1)
    {
      if(dataAndMC)
	{
	  switch(channel)
	    {
	      //B to D
	    case 0:
	      //mcIsDStar .. dCharge mcDCharge mcBCharge, bestBCharge dType==2 (dstar) 
	      sprintf(chargeAndStarSelection,"  abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 0==mcIsDStar");
	      break;
	      //B to DStar
	    case 1:
	      sprintf(chargeAndStarSelection,"   abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && mcIsDStar");
	      break;
	      //B0to D0
	    case 2:
	      sprintf(chargeAndStarSelection,"   abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar==0");
	      break;
	      //B0to DStar0
	    case 3:
	      sprintf(chargeAndStarSelection,"   abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && mcIsDStar");
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
void fitFractions(TTree** trees, TH1F** summedComponents, int numComponents,int numPions,int leptonId, int channel, bool dataTree, bool addNoise, TH1D* pulls)
{


  char channelString[500];
  getChannelString(channel,channelString);

  //since we name one of the channels -1
  int channelIdx=channel+1;

  float templateScaleFactor=1.0;
  if(dataTree)
    templateScaleFactor=0.2;

#ifdef PARTIAL_BOX
   templateScaleFactor=partialBoxFraction/5.0;
#endif

#ifdef MC_TEST
   templateScaleFactor=0.25;
#endif
   cout <<"using scale factor " << templateScaleFactor <<endl;
   for(int i=0;i<numComponents;i++)
     {
       summedComponents[i]->Sumw2();
       (*summedComponents)->SetFillStyle(1001);
       summedComponents[i]->Scale(templateScaleFactor);
     }
   cout <<"using scale factor: "<< templateScaleFactor <<endl;
  TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
  cout <<"fit fractions with " << numComponents <<" components"<<endl;

  //should be only 10 (or num components, this is just for safety)
  char* templateLegendNames[20];
  for(int i=0;i<numComponents;i++)
    {
      templateLegendNames[i]=new char[300];
    }
  vector<int> effectiveComponentsIndices;
  char histoName[2009];
  char drawCommand[2000];
  char buffer[2000];
  char corrBuffer[2000];
  int minCounts=0;
  ///// Used for Analysis note:  int fixThresholdCounts=300;
  int fixThresholdCounts=300;
  //  if(leptonId>0)
  //    fixThresholdCounts=1000;


  //set to the pionId 1 to merge components...
  int oneIdx=-1;
  bool combineDPiPi=true;

  ///make new templates that take into account that some contributions look the same, so probably create problems while fittign
  //same shape combination would be 5 with 7 (only one merger) --> alternatively used 4 mergers 5, 6, 7, 8,9 where all the others have small statistics
  int numMergers=0;
  if(numPions==oneIdx)
    numMergers=1;
  if(combineDPiPi)
    numMergers=2;
    //    numMergers=4;

  TH1F** templates=new TH1F*[numComponents-numMergers];

  //  fillTemplates(templates,summedComponents,templateLegendNames,numPions);
  cout <<"filling " << numComponents-numMergers << " templates " <<endl;

  if(numPions==oneIdx)
    {
      for(int i=0;i<5;i++)
	{
	  templates[i]=summedComponents[i];
	  templateLegendNames[i]=allLegendNames[i];
	  cout <<"filled template " << i <<endl;
	}
      sprintf(buffer,"%s_clonedAgain",summedComponents[5]->GetName());
      templates[5]=(TH1F*)summedComponents[5]->Clone(buffer);
      templates[5]->Add(summedComponents[6]);
      templates[5]->Add(summedComponents[7]);
      //           templates[5]->Add(summedComponents[8]);
      //      templates[5]->Add(summedComponents[9]);

      //      sprintf(buffer,"%s plus %s and %s and  %s and %s",allLegendNames[5],allLegendNames[6],allLegendNames[7],allLegendNames[8], allLegendNames[9]);
      sprintf(buffer,"%s plus %s and %s",allLegendNames[5],allLegendNames[6],allLegendNames[7]);
      templates[5]->SetTitle(buffer);
      sprintf(templateLegendNames[5],"%s",buffer);
      //            templates[6]=summedComponents[6];
      //      templateLegendNames[6]=allLegendNames[6];
      //      for(int i=10;i<numComponents;i++)
      for(int i=8;i<numComponents;i++)
	{
	  templates[i-numMergers]=summedComponents[i];
	  templateLegendNames[i-numMergers]=allLegendNames[i];
	  cout <<"filling template " << i-numMergers <<endl;
	}
    }
  else{
    if(combineDPiPi)
      {
       for(int i=0;i<4;i++)
	{
	  templates[i]=summedComponents[i];
	  templateLegendNames[i]=allLegendNames[i];
	  cout <<"filled template " << i <<endl;
	}
      sprintf(buffer,"%s_clonedAgain",summedComponents[4]->GetName());
      templates[4]=(TH1F*)summedComponents[4]->Clone(buffer);
      templates[4]->Add(summedComponents[6]);
      templates[4]->Add(summedComponents[8]);
      //           templates[5]->Add(summedComponents[8]);
      //      templates[5]->Add(summedComponents[9]);

      //      sprintf(buffer,"%s plus %s and %s and  %s and %s",allLegendNames[5],allLegendNames[6],allLegendNames[7],allLegendNames[8], allLegendNames[9]);
      sprintf(buffer,"%s plus %s and %s",allLegendNames[4],allLegendNames[6],allLegendNames[8]);
      templates[4]->SetTitle(buffer);
      sprintf(templateLegendNames[4],"%s",buffer);
      //            templates[6]=summedComponents[6];
      //      templateLegendNames[6]=allLegendNames[6];
      //      for(int i=10;i<numComponents;i++)
      templates[5]=summedComponents[5];
      templateLegendNames[5]=allLegendNames[5];
      templates[6]=summedComponents[7];
      templateLegendNames[6]=allLegendNames[7];
      //already got [8] (added to 4), so next up is 7 <-- 9
      for(int i=9;i<numComponents;i++)
	{
	  templates[i-numMergers]=summedComponents[i];
	  templateLegendNames[i-numMergers]=allLegendNames[i];
	  cout <<"filling template " << i-numMergers <<endl;
	}

      }
  else
    {
      for(int i=0;i<numComponents;i++)
	{
	  templates[i]=summedComponents[i];
	  templateLegendNames[i]=allLegendNames[i];
	}
    }
  }
  //////

  cout <<"did set up " << numComponents-numMergers << " templates" <<endl;

  //to save counts so we can fix the components which have too little counts
  vector<int> countsOfComponents;
  vector<int>  indexOfEffComp;

  ////from the example on the root web pages..

    TH1F *data;                              //data histogram
    //use first 4 trees (the last one is data)
    int treeCount=4;
    addCorrections(corrBuffer);
    if(dataTree)
      {
	treeCount=1;
      }
    cout <<"tree count: " << treeCount <<endl;


    char channelSelectionData[1000];
    char channelSelectionDataAndMC[1000];
    //void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions)
    getChargeAndStarSelection(channelSelectionData,channel,false,numPions);
    getChargeAndStarSelection(channelSelectionDataAndMC,channel,true,numPions);

    if(leptonId!=0)
	  {
		//for the MC_test we still have to do the tag corr....
	    if(dataTree)
	      {
#ifdef MC_TEST
		if(numPions==1)
		  sprintf(buffer,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		else
		  sprintf(buffer,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);

#else
		if(numPions==1)
		  {
		  sprintf(buffer,"%s "P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		  }
		else
		  {
		  sprintf(buffer,"%s ("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		  }

#endif
	      }
	    else
	      {
		if(numPions==1)
		  {
		    sprintf(buffer,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId, channelSelectionData);
				  //		  sprintf(buffer,"%s "P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId, channelSelectionData);
		  }
		else
		  {
		    sprintf(buffer,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
				  //		  sprintf(buffer,"%s ("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut,numPions,leptonId,channelSelectionData);
		  }
	      }
		  //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	  }
	else
	  {
	    if(dataTree)
	      {

		//for the MC_test we still have to do the tag corr....
#ifdef MC_TEST
		if(numPions==1)
		  {
		   		    sprintf(buffer,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
				    //		    sprintf(buffer,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
		else
		  {
		    sprintf(buffer,"%s tagCorr*("P2STRING"  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
#else
		if(numPions==1)
		  {
		  sprintf(buffer,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
		else
		  {
		  sprintf(buffer,"%s ("P2STRING"  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
#endif
	      }
	    else
	      {
		if(numPions==1)
		  {
		    sprintf(buffer,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions, channelSelectionData);
		    //		    sprintf(buffer,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions, channelSelectionData);
		  }
		else
		  {
		    		    sprintf(buffer,"%s tagCorr*("P2STRING"  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
				    //		    sprintf(buffer,"%s ("P2STRING"  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
	      }
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
	  }

    cout <<" treeCount: " << treeCount <<endl;
    for(int tc=0;tc<treeCount;tc++)
      {
	cout << " tc: " << tc<< " numPions: "<< numPions<<" leptId: "<< leptonId <<" channel: " << channelString<<endl;
	sprintf(histoName,"histo_Data_%d_pions_%d_leptonId_treeNr%d_%s",numPions,leptonId,tc, channelString);
	cout <<"draw command..." <<endl;
	TH1D* histo=new TH1D(histoName,histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	//	sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	sprintf(drawCommand,"mNu2 >> %s",histoName);
	cout <<" and the actual draw... "<< drawCommand<<"--> >" <<buffer<<endl;
	int counts=0;
	if(!dataTree)
	  {
	    counts=trees[tc]->Draw(drawCommand,(char*)buffer);
	  }
	else
	  {
	    //grab last tree which should be the data tree...
	    counts=trees[4]->Draw(drawCommand,(char*)buffer,"",maxDataTreeSize);
	  }
	cout <<"got " << counts <<" counts from data selected " <<endl;
	if(tc==0)
	  data=(TH1F*)gDirectory->Get(histoName);
	else
	  data->Add((TH1F*)gDirectory->Get(histoName));

		data->Sumw2();
	data->SetFillColor(glColorTable[0]->GetNumber());	
      }

    cout <<" done making data ..  sig ids: "<< SIG_IDX <<endl;
    cout <<"data integral: " << data->Integral() <<endl;
    double signalFraction=templates[SIG_IDX]->Integral()/data->Integral();
    cout <<"for BR, signal in MC is : " << templates[SIG_IDX]->Integral();
    double mcSignalIntegral=templates[SIG_IDX]->Integral();
    cout <<"signal Fraction estimated to be : " << signalFraction<<endl;
    cout <<"add noise? " << addNoise <<endl;
    getChannelString(channel,channelStringGlobal);
    glChannelIdx=channelIdx;

    double fitVal, fitErr;
    cout<<"trying usual fit function: " << endl;
    double* allFitVals=new double[numComponents-numMergers];
    int numEffective=0;
    double S=getFitSignal(data,templates,numComponents-numMergers, fitVal, fitErr, fixThresholdCounts,allFitVals, numEffective);  
    cout <<"got " << fitVal*data->Integral() << " signal counts, fraction: " << fitVal << " +- " << fitErr  << endl;
    /////-----------------
    //void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents)
    if(addNoise)
      {

	//scale data by 0.2 since otherwise it is the sum of the 5 streams...
	//	data->Scale(0.2);
	//-->do this in the 'performFractionFit... function'
	for(int nIt=0;nIt<1000;nIt++)
	  {
	    performFractionFitStabilityTest(templates,data,numComponents-numMergers,pulls,signalFraction, fixThresholdCounts);
	  }
      }
    /////////-----------

    TCanvas sampleData;
    data->Draw();
    sampleData.SaveAs("sampleData.png");
    sampleData.SaveAs("sampleData.pdf");
    sampleData.SaveAs("sampleData.eps");

  int numEffectiveComponents=0;
  for(int i=0;i<numComponents-numMergers;i++)
    {
      //special case of 'other BB' which looks like the signal for some reason...
      if(templates[i]->GetEntries()>minCounts)
	{
	  numEffectiveComponents++;
	  countsOfComponents.push_back(templates[i]->GetEntries());
	  indexOfEffComp.push_back(i);
	}

      //      if(i==6)
      //	templates[i]->Scale(100000);
    }
  cout <<" we have " << numEffectiveComponents<<" effective components " << endl;
  TObjArray *mc = new TObjArray(numEffectiveComponents);        // MC histograms are put in this array

  //numpions=1 --> comp 5 and 7 look the same
  for(int i=0;i<numComponents-numMergers;i++)
    {
      cout <<"template " << i <<" has " << templates[i]->GetEntries() <<" entries" <<endl;
      if(templates[i]->GetEntries()>minCounts)
	{
	  mc->Add(templates[i]);
	  effectiveComponentsIndices.push_back(i);
	}
      sprintf(buffer,"mcComponent_%d_leptId_%d_numPion_%d_%s.png",i,leptonId,numPions,channelString);

      templates[i]->Draw();
      sampleData.SaveAs(buffer);
      sprintf(buffer,"mcComponent_%d_leptId_%d_numPion_%d_%s.pdf",i,leptonId,numPions,channelString);
      sampleData.SaveAs(buffer);
      sprintf(buffer,"mcComponent_%d_leptId_%d_numPion_%d_%s.eps",i,leptonId,numPions,channelString);
      sampleData.SaveAs(buffer);
    }

  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  //default is 1e-6
  //    fit->GetFitter()->SetPrecision(0.001);
  //    fit->GetFitter()->SetPrecision(1e-8);

  for(int i=0;i<numEffectiveComponents;i++)
    {
      cout <<" looking at component " << i <<" have " << countsOfComponents[i] << " counts " <<endl;
      //i = 9 (after merger)is BB background
      bool partialBox=false;
#ifdef PARTIAL_BOX
      partialBox=true;
#endif

      //      if(!partialBox || channel==-1)
      if(!partialBox || channel==-1)
	{
	  if(countsOfComponents[i]>fixThresholdCounts  && (numPions!=oneIdx || indexOfEffComp[i]!=9))
	    {
	  //	  cout << "not fixed component " << i << " scale fact: "<< templateScaleFactor*fixThresholdCounts << " counts: "<< countsOfComponents[i] <<endl;
	      fit->Constrain(i,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	  //	  fit->UnConstrain(i);
	    }
	  else
	    {

	  //	  cout << "fixed component " << i << " scale fact: "<< templateScaleFactor*fixThresholdCounts << " counts: "<< countsOfComponents[i] <<endl;
	  //fit->Constrain(i,1.0,1.0);               ///not enough counts, so fix this component
	      ////---> this worked for MC part
	      //should do that, otherwise the parameter gets fixed to some random value
	      //	  fit->GetFitter()->SetParameter(i,buffer,iThFraction,0.0,0.0,1.0);
	  fit->GetFitter()->FixParameter(i);
	  fit->GetFitter()->FixParameter(i);
	      //	      ->Scale(1000000);


	      //could also achieve this by scaling this MC template up....
	    }
	}
      else
	{
	  if(effectiveComponentsIndices[i]!=SIG_IDX)
	    {
	      cout<< "eff comp index: "<< effectiveComponentsIndices[i] << " i: " << i << " SIG_IDX: " << SIG_IDX<<endl;
	      if((i==0 || i==2) && channel!=(-1)) //continuum or false channel
		{
		  fit->Constrain(i,0.0,1.0);
		}
	      else
		{
		  fit->Constrain(i,0.0,0.01);
		}
	      if(channel==1 ) // no overlap
		{
		  if(i==0|| i==2)
		    {
		      fit->GetFitter()->FixParameter(i);
		    }
		  else
		    {
		      fit->Constrain(i,0.0,0.01);
		    }
		}
	      if(channel==-1)
		{
		  if(i==0 || i==1)
		    {
		      fit->Constrain(i,0.0,1.0);
		    }
		  else
		    {
		      fit->Constrain(i,0.0,0.01);
		    }
		}
	    }
	  else
	    {
	      cout <<" i: " << i << " seems to be signal, constrain to 0-1" <<endl;
	      fit->Constrain(i,0.0,1.0);
	    }
      }


    }
      //                  fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
  Int_t status = fit->Fit();               // perform the fit
  std::cout << "fit status: " << status << std::endl;
  if (status == 0) {                       // check on fit status
    TCanvas c;
    TH1F* result = (TH1F*) fit->GetPlot();
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
    sprintf(buffer,"fracFitComp2_numPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelString);
    THStack* predComponents2=new THStack(buffer,buffer);
    cout <<"fraction fitter has " << fit->GetFitter()->GetNumberFreeParameters() << " free and " << fit->GetFitter()->GetNumberTotalParameters() <<" overall parameters" <<endl;

    //we want to flip the signal (index =2 ) so that it is later
    int signalIdx=-1;

    //has to sum to 1.0
    double totalFraction=0.0;
    for(int i=0;i<numEffectiveComponents;i++)
      {
	if(effectiveComponentsIndices[i]!=SIG_IDX)
	  {
	    cout <<"component " << i<< " (" << templateLegendNames[effectiveComponentsIndices[i] ] <<") is now : "<< fit->GetFitter()->GetParameter(i)<<" other meth: " << allFitVals[i]<<endl;
	    TH1F* mcComp=(TH1F*) fit->GetMCPrediction(i);
	    sprintf(buffer,"%s_Clone",mcComp->GetName());
	    TH1F* mcComp2=(TH1F*) mcComp->Clone();
	    //check what fraction of the total prediction this thing has now, if not consistent with fit, scale accordingly...
	        double mcPredInt=mcComp->Integral();
	    //shouldn't we take the data ?
	    ///	    double mcPredInt=data->Integral();
	    double fitFraction=fit->GetFitter()->GetParameter(i);
	    totalFraction+=fitFraction;
	    double factualFraction=mcPredInt/templatePredIntegral;
	    cout <<"mc Pred is " << factualFraction << " of total template prediction, should be : " << fitFraction<<", so scale by " << fitFraction/factualFraction<<endl;
	    mcComp->Scale(fitFraction/factualFraction);
	    //and scale so that data and mc integrals fit...
	    cout <<"integral of result: " << templatePredIntegral <<" data Integral: "<< data->Integral() <<endl;
	    cout <<"scaling to match data and mc integral by " << data->Integral()/mcPredInt<<endl;
	    //	    mcComp->Scale(data->Integral()/mcPredInt);
	    if(allFitVals[i]>0)
	      {
		cout<<"other method scale is " << fitFraction/allFitVals[i]<<endl;
		mcComp2->Scale(allFitVals[i]/factualFraction);
		//		mcComp2->Scale(data->Integral()/mcPredInt);
	      }
	    cout <<"comp with hist entries: "<< mcComp->GetEntries() <<" and " << templates[effectiveComponentsIndices[i] ]->GetEntries() <<": " << (float)mcComp->GetEntries()/(float)templates[effectiveComponentsIndices[i] ]->GetEntries()<<endl;
	    cout <<"comp with hist scale: "<< (float)mcComp->GetEntries()/(float)data->GetEntries() <<" and " << (float)templates[effectiveComponentsIndices[i] ]->GetEntries()/(float)data->GetEntries()<<endl;
	    //	    mcComp->SetFillColorAlpha(1.0);
	    mcComp->SetFillStyle(1001);
	    mcComp->SetFillColor(glColorTable[effectiveComponentsIndices[i]]->GetNumber());
	    mcComp2->SetFillStyle(1001);
	    mcComp2->SetFillColor(glColorTable[effectiveComponentsIndices[i]]->GetNumber());

	    predComponents->Add(mcComp);
	    predComponents2->Add(mcComp2);
	    legend->AddEntry(mcComp,templateLegendNames[effectiveComponentsIndices[i]],"f" );
	  }
	else
	  {
	    signalIdx=i;
	  }
      }
    //add the signal last...
    if(signalIdx>=0)
      {
	TH1F* mcComp=(TH1F*) fit->GetMCPrediction(signalIdx);
	sprintf(buffer,"%s_cloned",mcComp->GetName());
	TH1F* mcComp2=(TH1F*) mcComp->Clone(buffer);
	double mcPredInt=mcComp->Integral();
	cout <<"mcPred is: "<< mcPredInt <<endl;
	double fitFraction=fit->GetFitter()->GetParameter(signalIdx);
	cout <<"we still need "<< 1.0-totalFraction <<" of the data and have " << fitFraction<< " so missing " << 1.0-totalFraction-fitFraction <<endl;
	double miss=1.0-totalFraction-fitFraction;
#ifdef PARTIAL_BOX
	fitFraction+=miss/2;
#endif
	double factualFraction=mcPredInt/templatePredIntegral;
	cout <<"signal mc Pred is " << factualFraction << " of total template prediction, should be : " << fitFraction<<" (other method: " << allFitVals[signalIdx]<<"), so scale by " << fitFraction/factualFraction<<endl;
	mcComp->Scale(fitFraction/factualFraction);
	mcComp2->Scale(allFitVals[signalIdx]/factualFraction);
	cout <<"after scale mc pred is: "<< mcComp->Integral()<<endl;
	cout <<"BR ratio to MC is : "<< mcComp->Integral()/mcSignalIntegral<<endl;
	cout <<"BR2 ratio to MC is : "<< mcComp2->Integral()/mcSignalIntegral<<endl;
	mcComp->SetFillStyle(1001);
	mcComp->SetFillColor(glColorTable[SIG_IDX]->GetNumber());
	predComponents->Add(mcComp);
	mcComp2->SetFillStyle(1001);
	mcComp2->SetFillColor(glColorTable[SIG_IDX]->GetNumber());
	predComponents2->Add(mcComp2);
	legend->AddEntry(mcComp,templateLegendNames[SIG_IDX],"f" );
      }


    ////
    if(numPions==1)
      {
#ifdef PARTIAL_BOX
	if(leptonId==0)
	  {
	    data->GetYaxis()->SetRangeUser(0,3300*templateScaleFactor);
	    if(channel==1)
	      data->GetYaxis()->SetRangeUser(0,60);
	    if(channel==3)
	      data->GetYaxis()->SetRangeUser(0,20);
	  }
	else
	  {
	    data->GetYaxis()->SetRangeUser(0,1500*templateScaleFactor);
	  }
#else
	if(leptonId==0)
	  data->GetYaxis()->SetRangeUser(0,800*templateScaleFactor);
	else
	  data->GetYaxis()->SetRangeUser(0,250*templateScaleFactor);
#endif
      }
    else
      {
	if(leptonId==0)
	  {
	    if(combineDPiPi)
	      data->GetYaxis()->SetRangeUser(0,450*templateScaleFactor);
	    else
	      data->GetYaxis()->SetRangeUser(0,1600*templateScaleFactor);
	  }
	else
	  {
	    if(combineDPiPi)
	      data->GetYaxis()->SetRangeUser(0,250*templateScaleFactor);
	    else
	      data->GetYaxis()->SetRangeUser(0,800*templateScaleFactor);
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
    data->Draw("Ep");
    predComponents2->SetTitle("");
    //    predComponents->SetStats(0);
    //    predComponents->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
    predComponents2->Draw("hist same");
    data->SetLineWidth(2);
    data->SetStats(0);
    data->Draw("same Ep");
    legend->Draw();


    sprintf(buffer,"predComp2_numPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelString);
    c.SaveAs(buffer);
    sprintf(buffer,"predComp2_numPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelString);
    c.SaveAs(buffer);
    sprintf(buffer,"predComp2_numPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelString);
    c.SaveAs(buffer);
  }


  ////

}




void loadComponents(TFile* file,TH1F** components, TH1F** summedComponents, int numPions,int leptonId, int channel,  int numComponents, int numFiles)
{
  int channelIdx=channel+1;
  char channelString[500];
  getChannelString(channel,channelString);

  char buffer[2000];
  for(int iF=0;iF<numFiles;iF++)
    {
      for(int b=0;b<numComponents;b++)
	{
	  sprintf(buffer,"histo_If_%d_b_%d_numPions_%d_leptonId_%d_%s",iF,b,numPions,leptonId,channelString);
	  components[iF*10+b]=(TH1F*)file->Get(buffer);
	}
    }


  sprintf(buffer,"continuum_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[0]=(TH1F*)file->Get(buffer);

  cout <<"pointer : "<< summedComponents[0]<<endl;
  sprintf(buffer,"DDStar_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[1]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[2]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPi_WrongChannel_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[3]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPiPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[4]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DLnu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[5]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[5]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[6]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DStarLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[7]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DStarPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[8]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DStarPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[8]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DStarPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[8]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"OtherBB_%dLept_%dPions_%s",leptonId,numPions,channelString);
  summedComponents[9]=(TH1F*)file->Get(buffer);

  for(int i=0;i<10;i++)
    {
      cout <<"loaded component " << i << " with counts " << summedComponents[i]->GetEntries()<<endl;
    }
}

//trees: the input trees for the 4 MC files, components: The components for each file, so 4*11, summedComponents: The same, but summed over files
void getMCComponents(TTree** trees, TH1F** components, TH1F** summedComponents, int numPions,int leptonId, int channel)
{

  int channelIdx=channel+1;
  char channelSelectionData[1000];
  char channelSelectionDataAndMC[1000];
  char wrongChannelSelection[1000];

    //void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions)
  getChargeAndStarSelection(channelSelectionData,channel,false,numPions);
  getChargeAndStarSelection(channelSelectionDataAndMC,channel,true,numPions);
  getChargeAndStarSelection(wrongChannelSelection,channel,true,numPions,true);


  char channelString[500];
  getChannelString(channel,channelString);


  cout<<" getting pres for channel " << channelString <<" lowerCut: " << lowerCut[channelIdx] << " upperCut " << upperCut[channelIdx] <<" num bins: " << numBins[channelIdx] <<endl;
  cout <<"getting MC components..." <<endl;
  char histoName[2009];
  char drawCommand[2000];
  //just getting the MC components, so only looking at 4 files
  int numFiles=4;
  int numComponents=10;

  char buffer[2000];
  sprintf(buffer,"continuum_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hContinuum=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DDStar_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStar=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DDStarPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStarPi=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);

  sprintf(buffer,"DDStarPi_WrongChannel_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStarPiWrongChannel=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);

  sprintf(buffer,"DDStarPiPi_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDDStarPiPi=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DLnu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDlNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDPiPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DStarLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDStarlNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DStarPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDStarPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"DStarPiPiLNu_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hDStarPiPilNu=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(buffer,"OtherBB_%dLept_%dPions_%s",leptonId,numPions,channelString);
  TH1F* hOtherBB=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);


  char* selections[10];

  char bufferDDStar[2000];
  char bufferDDStarPi[2000];
  char bufferDDStarPiWrongChannel[2000];
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
  addCorrections(corrBuffer);

  if(leptonId!=0)
    {
      if(numPions==1)
	{
	  sprintf(buffer,"%s  tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		//	  sprintf(buffer,"%s  "P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
	}
      else
	{
	  sprintf(buffer,"%s  tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
	}
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d  ",numPions,leptonId);
    }
  else
    {
      if(numPions==1)
	{	
	  	  sprintf(buffer,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  //	  sprintf(buffer,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge) && %s ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
	}
      else
	{
	  sprintf(buffer,"%s tagCorr*("P2STRING"  && bestBCharge==((-1)*systemCharge) && %s",corrBuffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
	}
      //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0   ",numPions);
    }

  sprintf(bufferDDStar,"%s &&   foundAnyDDoubleStar==1 && sig_numPions==0 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1 && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
  //sprintf(bufferDDStar,"%s && foundAnyDDoubleStar==1)",buffer);
  sprintf(bufferDDStarPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu && %s)",buffer,channelSelectionDataAndMC);
  sprintf(bufferDDStarPiWrongChannel,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu && %s)",buffer,wrongChannelSelection);

  sprintf(bufferDDStarPiPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==2 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);

  sprintf(bufferDlNu,"%s && sig_DLNu )",buffer);
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

  selections[0]=bufferNoSelection;
  selections[1]=bufferDDStar;
  selections[2]=bufferDDStarPi;
  selections[3]=bufferDDStarPiWrongChannel;
  selections[4]=bufferDDStarPiPi;
  selections[5]=bufferDlNu;
  //  selections[5]=bufferDPilNu;
  selections[6]=bufferDPiPilNu;

  selections[7]=bufferDStarlNu;
  //  selections[8]=bufferDStarPilNu;
  selections[8]=bufferDStarPiPilNu;
  selections[9]=bufferAll;

  //
  // summedHistos: not differentiated between mixed and charged
  //in the end we use something like 		  summedHistos[b+1]->Add(result from selection[b]) (and continuum is treated extra)
  //
  TH1F**  summedHistos=summedComponents;
  summedHistos[0]=hContinuum;
  summedHistos[1]=hDDStar;
  summedHistos[2]=hDDStarPi;
  summedHistos[3]=hDDStarPiWrongChannel;
  summedHistos[4]=hDDStarPiPi;
  summedHistos[5]=hDlNu;
  //  summedHistos[5]=hDPilNu;
  summedHistos[6]=hDPiPilNu;
  summedHistos[7]=hDStarlNu;
  //  summedHistos[8]=hDStarPilNu;
  summedHistos[8]=hDStarPiPilNu;
  summedHistos[9]=hOtherBB;
  
  for(int iF=0;iF<numFiles;iF++)
    {
      cout <<" iF: "<< iF <<" numComps: " << numComponents<< " numFiles: "<< numFiles <<endl;
      int allCounts=0;
      int noSelCounts=0;
      for(int b=0;b<numComponents;b++)
	{
	  cout <<"b: "<< b << endl;
	  cout <<"creating histo idependently " << endl;
	  sprintf(histoName,"histo_If_%d_b_%d_numPions_%d_leptonId_%d_%s",iF,b,numPions,leptonId,channelString);
	  TH1F* h=new TH1F(histoName,histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	  sprintf(drawCommand,"mNu2 >> %s",histoName);
	  //	  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	  cout <<"draw command: " << drawCommand << ", selections: " << (char*) selections[b]<<endl;
	  cout <<"drawing tree " << iF <<endl;
	  int counts=trees[iF]->Draw(drawCommand,(char*)selections[b]);
	  if(b>0)
	    allCounts+=counts;
	  if(b==0)
	    noSelCounts=counts;
	  //  if(b==10)
	  if(b==9)
	    cout <<"all counts so far: "<< allCounts <<" no selection counts: "<< noSelCounts <<" difference: "<< noSelCounts-allCounts <<endl;
	  cout <<"got " << counts <<" counts selected " <<endl;
	  //do we have to clone this before we can return the histo?
	  TH1F* result=(TH1F*)gDirectory->Get(histoName);
	  result->SetFillColor(glColorTable[b]->GetNumber());
	  components[iF*10+b]=result;
	  //	  if(b!=10)
	    {
	      //uds or charm, but only use the noSelection
	      if(iF>=2 && b==0)
		{

		  cout <<"c adding histo with: "<< result->GetNbinsX() <<", "<< result->GetBinCenter(1) <<" to " << result->GetBinCenter(result->GetNbinsX())<<endl;
		  cout <<" to histo with: "<< summedHistos[b]->GetNbinsX() <<", "<< summedHistos[b]->GetBinCenter(1) <<" to " << summedHistos[b]->GetBinCenter(summedHistos[b]->GetNbinsX())<<endl;
		  hContinuum->Add(result);
		  cout <<"done " <<endl;
		  hContinuum->SetFillColor(glColorTable[0]->GetNumber());
		}
	      //mixed or charged. The b==0 case is the noSelection case, which only makes sense for the continuum
	      if(iF<2 && b>0)
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


		  cout <<" first zaxis limits: "<< result->GetZaxis()->GetXmin() <<" to " << result->GetZaxis()->GetXmax()<<endl;
		  cout <<"secondy zxis limits: "<< summedHistos[b]->GetZaxis()->GetXmin() <<" to " << summedHistos[b]->GetZaxis()->GetXmax()<<endl;
		  cout <<" first yaxis limits: "<< result->GetYaxis()->GetXmin() <<" to " << result->GetYaxis()->GetXmax()<<endl;
		  cout <<"secondy axis limits: "<< summedHistos[b]->GetYaxis()->GetXmin() <<" to " << summedHistos[b]->GetYaxis()->GetXmax()<<endl;
		  cout <<" first axis limits: "<< result->GetXaxis()->GetXmin() <<" to " << result->GetXaxis()->GetXmax()<<endl;
		  cout <<"second axis limits: "<< summedHistos[b]->GetXaxis()->GetXmin() <<" to " << summedHistos[b]->GetXaxis()->GetXmax()<<endl;
		  cout <<" adding histo with: "<< result->GetNbinsX() <<", "<< result->GetBinCenter(1) <<" to " << result->GetBinCenter(result->GetNbinsX())<<endl;
		  cout <<" to histo with: "<< summedHistos[b]->GetNbinsX() <<", "<< summedHistos[b]->GetBinCenter(1) <<" to " << summedHistos[b]->GetBinCenter(summedHistos[b]->GetNbinsX())<<endl;
		  summedHistos[b]->Add(result);
		  cout <<"done " <<endl;


		  //for the color, make it consistent to the other stacks
		  summedHistos[b]->SetFillColor(glColorTable[b]->GetNumber());
		}
	    }
	}
    }
  cout <<"done with getMCComponents" <<endl;
};





void saveStack(TH1F** components, TH1F** summedComponents, int numPions, int leptonId, int channel)
{
  int channelIdx=channel+1;
  //  int printOrder[]={0,10,1,3,4,5,9,7,8,6,2};
  int printOrder[]={0,9,1,4,5,8,3,6,7,2};
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
  char* legendNames[10];

  int numFiles=4;
  char* fileNames[numFiles];
  for(int i=0;i<numFiles;i++)
    {
      fileNames[i]=new char[400];
    }


  for(int i=0;i<10;i++)
    {
      legendNames[i]=new char [500];
      allLegendNames[i]=new char[500];
    }

  char channelString[500];
  getChannelString(channel,channelString);


  sprintf(legendNames[0]," no Selection");
  sprintf(legendNames[1],"B #rightarrow D Double Star X #rightarrow D^{(*)} l #nu (no non-res Dn#pi l#nu)");
  sprintf(legendNames[2],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu (no non-res Dn#pi l#nu), %s",channelString);
  sprintf(legendNames[3],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu (no non-res Dn#pi l#nu)-->Wrong Channel");
  sprintf(legendNames[4],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi #pi l #nu (no non-res Dn#pi l#nu)");
  sprintf(legendNames[5],"D l #nu");
  //  sprintf(legendNames[5],"D #pi l #nu");
  sprintf(legendNames[6],"D #pi #pi l #nu");

  sprintf(legendNames[7],"D* l #nu");
  //  sprintf(legendNames[8],"D* #pi l #nu");
  sprintf(legendNames[8],"D* #pi #pi l #nu");
  sprintf(legendNames[9]," no D(*)n #pi l#nu ");

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
  sprintf(allLegendNames[9]," other B B ");


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
      for(int b=0;b<10;b++)
	{
	  int index=printOrder[b];
	  cout <<"b: " << b <<" index: " << index<<endl;
	  sprintf(histoName,"histo_If_%d_index_%d_numPions_%d_leptonId_%d_%s",iF,index,numPions,leptonId,channelString);
	  sprintf(outFileName,"%s.png",histoName);
	  //	  components[iF*9+b]=result;
	  cout <<"grabbing component.." << iF*10+index <<endl;
	  TH1F* result=(TH1F*)components[iF*10+index];
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
	      result->Draw();
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
      stacks[iF]->Draw();
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
  for(int i=0;i<10;i++)
    {
      int index=printOrder[i];
      all.Add(summedHistos[index]);
      legend->AddEntry(summedHistos[index],allLegendNames[index],"f");
    }
  TCanvas c;
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s",numPions,leptonId,channelString);
  //  all.SetTitle(buffer);
  all.SetTitle("");
  //  all.GetYaxis()->SetRangeUser(0,1000);
  //  all.SetStats(0);

  all.Draw();
  legend->Draw();
  c.Update();
  all.GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
  c.Modified();
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.png",numPions,leptonId,channelString);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.root",numPions,leptonId,channelString);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.pdf",numPions,leptonId,channelString);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId_%s.eps",numPions,leptonId,channelString);
  c.SaveAs(buffer);


}


void doSidebandComparison(TTree* mcTree, TTree* dataTree,int leptonId, int numPions, TH1F** lowerSidebandMC, TH1F** upperSidebandMC, TH1F** lowerSidebandData, TH1F** upperSidebandData, int channel)
{
  int channelIdx=channel+1;
  float upperSidebandTop=2.0;
  float upperSidebandBottom=0.5;

  float lowerSidebandTop=-0.5;
  float lowerSidebandBottom=-1.0;

  char upperSBSelection[2000];
  char lowerSBSelection[2000];

  char upperSBSelectionData[2000];
  char lowerSBSelectionData[2000];

  char histoName[200];
  char drawCommand[200];
  char buffer[2000];

  char channelString[500];
  char channelSelectionData[1000];

  getChannelString(channel,channelString);
  getChargeAndStarSelection(channelSelectionData,channel,false,numPions);




  addCorrections(buffer);
  //select sidebands from all (need to redo all components because we select a different range
  if(leptonId!=0)
    {
      if(numPions==1)
	{
	  sprintf(upperSBSelection,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	  //	sprintf(upperSBSelection,"%s "P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	}
      else
	{
	sprintf(upperSBSelection,"%s tagCorr*("P2STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
	}


      if(numPions==1)
	sprintf(upperSBSelectionData,P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ", upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
      else
	sprintf(upperSBSelectionData," ("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ", upperSidebandTop,upperSidebandBottom,numPions,leptonId, channelSelectionData);
     
      if(numPions==1)	
	sprintf(lowerSBSelection,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);	
      else
	sprintf(lowerSBSelection,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);

      if(numPions==1)
	sprintf(lowerSBSelectionData,P1STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);
      else
	sprintf(lowerSBSelectionData," ("P2STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",lowerSidebandTop,lowerSidebandBottom,numPions,leptonId, channelSelectionData);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
      if(numPions==1)
	sprintf(upperSBSelection,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
      else
	sprintf(upperSBSelection,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);

      if(numPions==1)
	sprintf(upperSBSelectionData,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
      else
	sprintf(upperSBSelectionData,"%s ("P2STRING" && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, upperSidebandTop,upperSidebandBottom,numPions, channelSelectionData);
        
      if(numPions==1)
      sprintf(lowerSBSelection,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);
      else
      sprintf(lowerSBSelection,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);

      if(numPions==1)
	sprintf(lowerSBSelectionData,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge)  && %s) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);
      else
	sprintf(lowerSBSelectionData,"%s ("P2STRING" && bestBCharge==((-1)*systemCharge) && %s ) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions, channelSelectionData);

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

      if(numPions==1)
	sprintf(sameChargeSelection,"%s tagCorr*"P1STRING" && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      else
	sprintf(sameChargeSelection,"%s tagCorr*("P2STRING"&& bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      

      if(numPions==1)
        sprintf(chargeNeutralSelection,"%s tagCorr*"P1STRING" && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      else
        sprintf(chargeNeutralSelection,"%s tagCorr*("P2STRING"&& bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);

      if(numPions==1)
      sprintf(sameChargeSelectionData," "P1STRING" && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ", upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      else
      sprintf(sameChargeSelectionData," ("P2STRING"&& bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d && %s)  ", upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      
      if(numPions==1)
        sprintf(chargeNeutralSelectionData," "P1STRING" && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      else
        sprintf(chargeNeutralSelectionData," ("P2STRING"&& bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d && %s)  ",upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
      if(numPions==1)
	sprintf(sameChargeSelection,"%s tagCorr*"P1STRING"  && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      else
	sprintf(sameChargeSelection,"%s tagCorr*("P2STRING" && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);

      if(numPions==1)
	sprintf(chargeNeutralSelection,"%s tagCorr*"P1STRING"  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      else
	sprintf(chargeNeutralSelection,"%s tagCorr*("P2STRING" && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",buffer, upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);


      if(numPions==1)
	sprintf(sameChargeSelectionData," "P1STRING"  && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      else
	sprintf(sameChargeSelectionData," ("P2STRING" && bestBCharge==systemCharge && abs(bestBCharge)==1  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);

      if(numPions==1)
	sprintf(chargeNeutralSelectionData," "P1STRING"  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      else
	sprintf(chargeNeutralSelectionData," ("P2STRING" && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0))  && %s) ",upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
    }
  sprintf(histoName,"sameChargeMC_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_WSMC=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(drawCommand,"mNu2 >> %s",histoName);
  //  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  int counts=mcTree->Draw(drawCommand,(char*)sameChargeSelection);
  cout <<"got " << counts <<" counts from mc same charge selected " <<endl;
  *sameChargeMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"neutralChargeMC_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_NCMC=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(drawCommand,"mNu2 >> %s",histoName);
  //  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  counts=mcTree->Draw(drawCommand,(char*)chargeNeutralSelection);
  cout <<"got " << counts <<" counts from mc charge neutral selected " <<endl;
  *chargeNeutralMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"sameChargeData_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_SCData=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(drawCommand,"mNu2 >> %s",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  //  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  counts=dataTree->Draw(drawCommand,(char*)sameChargeSelectionData);
  cout <<"got " << counts <<" counts from data same charge selected " <<endl;
  *sameChargeData=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"neutralChargeData_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
  TH1D* histo_NCData=new TH1D(histoName,histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
  sprintf(drawCommand,"mNu2 >> %s",histoName);
  //  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut[channelIdx],upperCut[channelIdx]);
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
      double nv=0;
      if(bc>10)
	{
	  nv=rnd->Gaus(bc,h->GetBinError(i));
	}
      else{
	//e.g. if errors were computed with 10000 counts and then scaled by 10, we would have an error of 10 and bin count of 1000
	///---> scale factor 1000/(10*10)=10
	float scaleFactor=h->GetBinContent(i)/(h->GetBinError(i)*h->GetBinError(i));
	bc=bc*scaleFactor;
	nv=rnd->Poisson(bc);
	nv=nv/scaleFactor;
      }

      //      cout <<"for bin  " << i <<" old value : "<< bc <<" new : "<< nv <<endl;
      //      cout <<"old error: " << h->GetBinError(i);
      if(nv>=0)
     	h->SetBinContent(i,nv);
      //      cout <<"new error: "<< h->GetBinError(i) <<endl;;
    }
}

void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents, TH1D* pulls,double signalFraction, int fixThresholdCounts)
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
	  //	  templates[i]->GetSumw2()->Set(0);
	  //	  //make the error essentially zero

	  ///-->check if this has an impact on the pulls...
	  //	  	  if(i==6)
	  //	  	    templates[i]->Scale(10000);
	  addPoissonNoise(templates[i]);
	  //	  templates[i]->Scale(1000);
      	}
      //      cout <<"template " << i <<" has " << templates[i]->Integral() <<" int and " << templates[i]->GetEntries()<<endl;
    }


  float SData=templates[SIG_IDX]->Integral();


  sprintf(buffer,"%s_Clone_it%d",dataOrg->GetName(),iteration);
  TH1F* data=(TH1F*)dataOrg->Clone(buffer);
  cout <<"doing the data " << endl;
  data->GetSumw2()->Set(0);
  data->Scale(0.2);
  addPoissonNoise(data);
  float newSignalFraction=SData/(5*data->Integral());
  cout <<"data has " <<SData <<" counts"<<endl;
  double fitVal=0;
  double fitErr=0;
  double* allFitVals=new double[numComponents];
  int numEffective=0;
  if(SData>0)
    {

      S=getFitSignal(data,templates,numComponents, fitVal, fitErr, fixThresholdCounts,allFitVals, numEffective);  
      cout <<"data integral: "<< data->Integral()<< " entries: " << data->GetEntries()<<endl;
      cout <<"fitVal: " << fitVal <<" fitErr: "<< fitErr <<endl;
      cout <<"S is : " << S << " fraction times data: " << data->Integral()*fitVal <<" or " << data->GetEntries()*fitVal<<endl;
      if(S>0&& fitErr>0.0)
	{
	  cout <<"Signal significance: " << fitVal/fitErr <<endl;
	  sigSignificance[glChannelIdx]->Fill(fitVal/fitErr);
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



float getFitSignal(TH1F* data,TH1F** templates,int numTemplates, double& fitVal, double& fitErr, int fixThresholdCounts, double* allFitVals, int& numEffective)
{
  cout <<"expecting " << numTemplates << " templates " << endl;
  float S=0;
  vector<int> countsOfComponents;
  vector<float> integralOfComponents;
  vector<int> effectiveComponentsIndices; 
  int numEffectiveComponents=0;
  TObjArray *mc = new TObjArray(numTemplates);        // MC histograms are put in this array
  int signalIndex=0;
  char buffer[400];
  TCanvas c;

  double templateIntegral=0.0;
  pCount++;
  for(int i=0;i<numTemplates;i++)
    {
      if(templates[i]->GetEntries()>0)
	{
	  cout <<"template " << i << " ("<<templates[i]->GetName() <<")  has " << templates[i]->GetEntries()<<endl;
	  numEffectiveComponents++;

	  mc->Add(templates[i]);
	  sprintf(buffer,"poissonTemplate%d_%s.png",i,channelStringGlobal);
	  templates[i]->Draw();
	  if((pCount%100)==0)
	    {
	      c.SaveAs(buffer);
	      sprintf(buffer,"poissonTemplate%d_%s.pdf",i,channelStringGlobal);
	      c.SaveAs(buffer);
	      sprintf(buffer,"poissonTemplate%d_%s.eps",i,channelStringGlobal);
	      c.SaveAs(buffer);
	    }
	  if(i==SIG_IDX)
	    {
	      signalIndex=effectiveComponentsIndices.size();
	      cout <<"set signalIndex to " << signalIndex<<endl;
	    }
	  effectiveComponentsIndices.push_back(i);
	  //the same as the test effectiveComponentIndices[i]==2 from the regular fit
	  countsOfComponents.push_back(templates[i]->GetEntries());
	  double tempIntegral=templates[i]->Integral();
	  integralOfComponents.push_back(tempIntegral);
	  if(i==SIG_IDX)
	    {
	      //only add half, since on average the BR for the signal in MC is twice as high as the PDG value
	      //	      templateIntegral+=(tempIntegral/2);
	      templateIntegral+=tempIntegral;
	    }
	  else{
	    templateIntegral+=tempIntegral;
	  }
	}
    }
  double dataIntegral=data->Integral();
  cout <<"done with integral: "<< templateIntegral <<", data integral: "<< dataIntegral <<endl;
// to scale roughly to the 5 (or 4) streams used. For the poisson noise add the two integrals should be roughly the same
  if(dataIntegral<3*templateIntegral)
    dataIntegral*=5; 

  data->Draw();
  cout <<"done with data draw.." << endl;
  if(((pCount-1)%100)==0)
    {
      sprintf(buffer,"poissonData_%s.png",channelStringGlobal);
      c.SaveAs(buffer);
      sprintf(buffer,"poissonData_%s.png",channelStringGlobal);
      c.SaveAs(buffer);
    }
  cout <<"overall integral of templates: " << templateIntegral <<" of data: "<< data->Integral()<<endl;
  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  for(int i=0;i<numEffectiveComponents;i++)
    {
      //use the leptonId=0
      //  int fixThresholdCounts=2000;
      if(countsOfComponents[i]>fixThresholdCounts|| i==signalIndex)
	{
	  //anyways good to give decent start value. The only caveat is that the signal is most likely less than what we think in MC. So the signal fraction is proably a bit overestimated and 
	  //the other fractions underestimated. But in principel we can calculate the best start value by taking the pdg value for the signal (i.e. half es much as MC pred)
	  double iThFraction=integralOfComponents[i]/templateIntegral;
	  //aim for roughly constant counts for these backgrounds. So if we see much less data than expected the fraction should rise
	  iThFraction*templateIntegral/dataIntegral;
	  if(i==signalIndex)
	    {
	      //two the factor two mojo...

	      //	      iThFraction=integralOfComponents[i]/(2*templateIntegral);
	      //	      iThFraction=integralOfComponents[i]/(templateIntegral);
	      cout << " setting signal fraction to  " << iThFraction <<endl;
	    }
	  else{
	    cout << " setting non signal fraction " << i <<" ("<<templates[effectiveComponentsIndices[i]]->GetName() <<")  to  " << iThFraction <<endl;
	  }
	  sprintf(buffer,"para%d",i);
	  fit->GetFitter()->SetParameter(i,buffer,iThFraction,0.0,0.0,1.0);
	  //does the set parameter do already the same as constrain?
	  //	  fit->Constrain(i,0.0,1.0);               // constrain fraction i to be between 0 and 1
	}
      else
	{
	  double iThFraction=integralOfComponents[i]/templateIntegral;
	  
	  cout <<"fixing parameter " << i <<" ("<<templates[effectiveComponentsIndices[i]]->GetName() <<") to " << iThFraction <<endl;
	  sprintf(buffer,"para%d",i);
	  fit->GetFitter()->SetParameter(i,buffer,iThFraction,0.0,0.0,1.0);
	  fit->GetFitter()->FixParameter(i);
	}
    }

  Int_t status = fit->Fit();               // perform the fit
  std::cout << "fit status: " << status << std::endl;
                       // check on fit status

  if (status == 0) 
    {
      numEffective=numEffectiveComponents;
    //we want to flip the signal (index =2 ) so that it is later
    //	if(effectiveComponentsIndices[i]!=2)
      cout <<"getting MC pred " << signalIndex <<endl;
      TH1F* mcComp=(TH1F*) fit->GetMCPrediction(signalIndex);
      S=mcComp->Integral();
      fit->GetResult(signalIndex,fitVal,fitErr);
      for(int i=0;i<numEffectiveComponents;i++)
	{
	  cout <<"effective component " << i<< " has " << fit->GetMCPrediction(i)->Integral() <<" integral and  " << fit->GetMCPrediction(i)->GetEntries() <<" counts, fit: ";
	  double temp=0.0;
	  double tempErr=0.0;
	  fit->GetResult(i, temp, tempErr);
	  cout << temp<<endl;
	  allFitVals[i]=temp;
	}
    }

  return S;
}
void fillTemplates(TH1F** templates,TH1F** summedComponents,char* templateLegendNames,int numPions)
{


};
#endif

