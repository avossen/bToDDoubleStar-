#include "TStyle.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TColor.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include <math.h>
#include "TGraph.h"

#include "doAnalysisCombined.h"

extern  void loadBinning();
//number of all 5 streams: 41510082

using namespace std;
int main(int argc, char** argv)
{

  loadBinning();
  if(argc!=8)
    {
    cout <<"argc: "<< argc<<" need six file names and numPions.." <<endl;
    exit(0);
    }
  int nMax[6];

  int redFactor=5;
  //=41510082/20;

  int numPions=atoi(argv[7]);
  //  cout <<"looking at filename; " << filename << " and " << numPions<< " pions " <<endl;
  //FOM S/sqrt(S+B)

  TFile* file[6];
  for(int i=1;i<7;i++)
    {
      cout <<"file " << i-1 <<": " << argv[i] <<endl;
      file[i-1]=new TFile(argv[i]);
    }

  TTree* tree[6];
  for(int i=0;i<6;i++)
    {
      tree[i] = (TTree*)file[i]->Get("DataTree");
      nMax[i]=tree[i]->GetEntries()/redFactor;
      cout <<"setting nMax " << i << " to: " << nMax[i] <<endl;
    }


  int leptonId=0;
  //let's use the boundaries of -0.4 to 1.0 for the n=1 case and -0.4 to 0.4 for n=2

  char corrBuffer[]="PIDCorrection*D_DecayCorr*B_DecayCorr*";
  char histoNameSig[200];
  char histoNameBg[200];
  char buffer[3000];
  char bufferBg[3000];
  char bufferSig[3000];
  char drawCommandSig[3000];
  char drawCommandBg[3000];
  int numBins=50;

  int iF=0;
  int b=0;


  float x[300];
  float yFoM[300];
  float yS[300];
  float yB[300];

  float mixedFactor=1.0/0.852255195;
  float chargedFactor=1.0/0.9266338302;



  //  char huschleMC_lumi_corr[1000];
  //  sprintf(huschleMC_lumi_corr,"(1+foundAnyDDoubleStar*%f)*",huschleLumiFactor);


    for(int iC=-1;iC<4;iC++)
  //  for(int iC=0;iC<2;iC++)
    {


      char channelString[1000];
      char channelSelection[1000];
      char channelSelectionWrongChannel[1000];
      //void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel)
      getChargeAndStarSelection(channelSelection,iC,true,numPions,false);
      getChargeAndStarSelection(channelSelectionWrongChannel,iC,true,numPions,true);
      getChannelString(iC,channelString);
      float lowerCut=-0.4;
      float upperCut=1.0;
      if(numPions==2)
	{
	  upperCut=0.4;
	}

      float upperPCut=-2.0;
      float pCutStep=-0.1;
      int numSteps=30;

      
      for(int i=0; i<numSteps;i++)
	{

	float logPCut=upperPCut+i*pCutStep;


	  double bgCounts=0;
	  double signalCounts=0;
	  for(int iT=0;iT<6;iT++)
	    {
	      sprintf(histoNameSig,"histoSig_Tree_%d_iC_%d_step_%d",iT,iC,i);
	      sprintf(histoNameBg,"histoBg_Tree_%d_iC_%d_step_%d",iT,iC,i);
	      sprintf(drawCommandSig,"mNu2 >> %s(%d,%f,%f)",histoNameSig,numBins,lowerCut,upperCut);
	      sprintf(drawCommandBg,"mNu2 >> %s(%d,%f,%f)",histoNameBg,numBins,lowerCut,upperCut);
	      float huschleLumiFactor=mixedFactor;
	      if(iT>1)
		huschleLumiFactor=0;
	      if(iT==1)
		huschleLumiFactor=chargedFactor;
	      //there should be no signal in the uds, charm and the XlNu files
	      huschleLumiFactor-=1.0;
	      sprintf(buffer,"%s (1+foundAnyDDoubleStar*%f)*tagCorr*(recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > %f && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge))  ",corrBuffer,huschleLumiFactor,upperCut,lowerCut,numPions,logPCut);
	      
	      if(numPions==1)
		{
		  sprintf(bufferBg,"%s && (%s || !foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu)",buffer,channelSelectionWrongChannel);
		  //fix sig_numPions to 1,because that is our main signal
		  sprintf(bufferSig,"%s && %s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu",buffer, channelSelection);
		}
	      else
		{
		  sprintf(bufferBg,"%s && (%s || !foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions != 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu)",buffer, channelSelectionWrongChannel);
		  //fix sig_numPions to 2,because that is our main signal
		  sprintf(bufferSig,"%s && %s && foundAnyDDoubleStar==1 && sig_numPions==2 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu",buffer, channelSelection);
		}
	      cout <<"sig selection: "<< bufferSig <<endl;
	      cout <<"bg selection: "<< bufferBg <<endl;
	      double bgCountsLoc=0;
	      cout <<"using nMax of " << nMax[iT] <<endl;
	      bgCounts+=tree[iT]->Draw(drawCommandBg,(char*)bufferBg, "",nMax[iT]);
	      bgCounts+=bgCountsLoc;
	      TH1F* result=(TH1F*)gDirectory->Get(histoNameBg);
	      bgCountsLoc=result->Integral();
	      double signalCountsLoc=0;
	      tree[iT]->Draw(drawCommandSig,(char*)bufferSig, "",nMax[iT]);
	      result=(TH1F*)gDirectory->Get(histoNameSig);
	      signalCountsLoc=result->Integral();
	      signalCounts+=signalCountsLoc;

	      cout <<"tree " << iT << " got " << bgCountsLoc << " bg counts, " << signalCountsLoc <<" signal counts, FoM: " << signalCountsLoc/(double)(sqrt(signalCountsLoc+bgCountsLoc))<<endl;
	      cout <<"total bg: "<< bgCounts <<" total Sig: " << signalCounts <<endl;
	    }

	float fom=0;
	if(signalCounts+bgCounts>0)
	  fom=(float)signalCounts/(float)sqrt(signalCounts+bgCounts);
	x[i]=logPCut;
	cout <<"fom: " << fom << " x: " << x[i] << " i: " << i <<endl;
	
	yFoM[i]=fom;
	yS[i]=(float)signalCounts;
	yB[i]=(float)bgCounts;
      }
    cout <<"done, doing graphs...num points: " << numSteps<<endl;
    TGraph grFoM(numSteps,x,yFoM);
    cout <<"done wiht graph " <<endl;
    //    TGraph grFoS(numSteps,x);
    //    TGraph grFoB(numSteps,x);
    TCanvas c;
    cout <<"canvas.." <<endl;
    grFoM.GetYaxis()->SetTitle("FoM");
    cout <<"x axis " << endl;
    grFoM.GetXaxis()->SetTitle("min log P");
    cout <<"drawing.." <<endl;
    grFoM.Draw("AP*");
    sprintf(buffer,"fomLogProb_%d_%s.png",numPions,channelString);
    c.SaveAs(buffer);

    }
}
