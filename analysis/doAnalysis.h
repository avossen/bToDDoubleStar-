#ifndef DO_ANALYSIS_H
#define DO_ANALYSIS_H

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

using namespace std;

//trees: the input trees for the 4 MC files, components: The components for each file, so 4*11, summedComponents: The same, but summed over files
void getMCComponents(TTree** trees, TH1F** components, TH1F** summedComponents, int numPions,int leptonId)
{
    char histoName[2009];
    char drawCommand[2000];
  int numFiles=4;
  int numComponents=11;



    TColor* colorTable[11];
    colorTable[0]=gROOT->GetColor(kBlue);
    colorTable[1]=gROOT->GetColor(kRed);
    colorTable[2]=gROOT->GetColor(kYellow);
    colorTable[3]=gROOT->GetColor(kMagenta);
    colorTable[4]=gROOT->GetColor(kGreen);
    colorTable[5]=gROOT->GetColor(kBlack);
    colorTable[6]=gROOT->GetColor(kWhite);
    colorTable[7]=gROOT->GetColor(kGray);
    colorTable[8]=gROOT->GetColor(kViolet);
    colorTable[9]=gROOT->GetColor(kCyan);
    colorTable[10]=gROOT->GetColor(kTeal);



  int numBins=100;

  TH1F* hContinuum=new TH1F("continuum","continuum",numBins,-1,1);
  TH1F* hDDStar=new TH1F("DDStar","DDstar",numBins,-1,1);
  TH1F* hDDStarPi=new TH1F("DDStarPi","DDstarPi",numBins,-1,1);
  TH1F* hDDStarPiPi=new TH1F("DDStarPiPi","DDStarPiPi",numBins,-1,1);
  TH1F* hDlNu=new TH1F("DLnu","DlNu",numBins,-1,1);
  TH1F* hDPilNu=new TH1F("DPiLNu","DPiLNu",numBins,-1,1);
  TH1F* hDPiPilNu=new TH1F("DPiPiLNu","DPiPiLNu",numBins,-1,1);
  TH1F* hDStarlNu=new TH1F("DStarLNu","DStarLNu",numBins,-1,1);
  TH1F* hDStarPilNu=new TH1F("DStarPiLNu","DStarPiLNu",numBins,-1,1);
  TH1F* hDStarPiPilNu=new TH1F("DStarPiPiLNu","DStarPiPiLNu",numBins,-1,1);
  TH1F* hOtherBB=new TH1F("OtherBB","OtherBB",numBins,-1,1);

  TH1F** histos[4];
  for(int i=0;i<4;i++)
    {
      histos[i]=new TH1F*[7];
    }

  char* selections[11];

    char buffer[2000];

    char bufferDDStar[2000];
    char bufferDDStarPi[2000];
    char bufferDDStarPiPi[2000];

    char bufferDlNu[2000];
    char bufferDPilNu[2000];
    char bufferDPiPilNu[2000];

    char bufferDStarlNu[2000];
    char bufferDStarPilNu[2000];
    char bufferDStarPiPilNu[2000];

    char bufferAll[2000];
    char bufferNoSelection[2000];


    selections[0]=bufferDDStar;
    selections[1]=bufferDDStarPi;
    selections[2]=bufferDDStarPiPi;
    selections[3]=bufferDlNu;
    selections[4]=bufferDPilNu;
    selections[5]=bufferDPiPilNu;

    selections[6]=bufferDStarlNu;
    selections[7]=bufferDStarPilNu;
    selections[8]=bufferDStarPiPilNu;

    selections[9]=bufferAll;
    selections[10]=bufferNoSelection;


    if(leptonId!=0)
      {
	sprintf(buffer,"tagCorr*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  ",numPions,leptonId);
      }
    else
      {
	sprintf(buffer,"tagCorr*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge==((-1)*systemCharge) ",numPions);
      }

    sprintf(bufferDDStar,"%s && foundAnyDDoubleStar==1 && sig_numPions==0 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1 && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
    //sprintf(bufferDDStar,"%s && foundAnyDDoubleStar==1)",buffer);
    sprintf(bufferDDStarPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
    sprintf(bufferDDStarPiPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==2 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);

    sprintf(bufferDlNu,"%s && sig_DLNu )",buffer);
    sprintf(bufferDPilNu,"%s && sig_DPiLNu&& !sig_DLNu  )",buffer);

    sprintf(bufferDPiPilNu,"%s && sig_DPiPiLNu && !sig_DLNu && !sig_DPiLNu)",buffer);

    sprintf(bufferDStarlNu,"%s &&  sig_DStarLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu)",buffer);
    sprintf(bufferDStarPilNu,"%s  && sig_DStarPiLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu)",buffer);
    sprintf(bufferDStarPiPilNu,"%s  && sig_DStarPiPiLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu)",buffer);

    sprintf(bufferAll,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu|| sig_DStarPiPiLNu) && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);

    sprintf(bufferNoSelection,"%s)",buffer);
    //    sprintf(bufferAll,"%s)",buffer);


    cout <<"selection for DDstar: "<< bufferDDStar <<endl<<" DDstarPi: "<< bufferDDStarPi<<endl<< " DDstarPiPi: "<< bufferDDStarPiPi <<endl;
    cout <<" DlNu: " << bufferDlNu<<endl <<" DPilNu:"<< bufferDPilNu << endl << " DPiPilNu: "<< bufferDPiPilNu <<endl;
    cout <<" all: "<< bufferAll <<endl;

    TH1F** summedHistos=summedComponents;
  summedHistos[0]=hContinuum;
  summedHistos[1]=hDDStar;
  summedHistos[2]=hDDStarPi;
  summedHistos[3]=hDDStarPiPi;
  summedHistos[4]=hDlNu;
  summedHistos[5]=hDPilNu;
  summedHistos[6]=hDPiPilNu;
  summedHistos[7]=hDStarlNu;
  summedHistos[8]=hDStarPilNu;
  summedHistos[9]=hDStarPiPilNu;
  summedHistos[10]=hOtherBB;


    for(int iF=0;iF<numFiles;iF++)
      {
	int allCounts=0;
	for(int b=0;b<numComponents;b++)
	  {
	    sprintf(histoName,"histo_If_%d_b_%d",iF,b);
	    sprintf(drawCommand,"mNu2 >> %s(%d,-1,1)",histoName,numBins);
	    cout <<"draw command: " << drawCommand << ", selections: " << (char*) selections[b]<<endl;
	    cout <<"drawing tree " << iF <<endl;
	    int counts=trees[iF]->Draw(drawCommand,(char*)selections[b]);
	    if(b<10)
	      allCounts+=counts;
	    if(b==10)
	      cout <<"all counts so far: "<< allCounts <<" no selection counts: "<< counts <<" difference: "<< counts-allCounts <<endl;
	    cout <<"got " << counts <<" counts selected " <<endl;
	    //do we have to clone this before we can return the histo?
	    TH1F* result=(TH1F*)gDirectory->Get(histoName);
	    result->SetFillColor(colorTable[b]->GetNumber());
	    components[iF*11+b]=result;

	    if(b!=10)
	      {
		//uds or charm
		if(iF>=2)
		  {
		    hContinuum->Add(result);
		    hContinuum->SetFillColor(colorTable[10]->GetNumber());
		  }
		//mixed or charged
		if(iF<2)
		  {
		    summedHistos[b+1]->Add(result);
		    //for the color, make it consistent to the other stacks
		    summedHistos[b]->SetFillColor(colorTable[b]->GetNumber());
		  }
	      }

	

	  }


      }



};


void doSidebandComparison(TTree** trees, TH1F** sidebands)
{

}


#endif
