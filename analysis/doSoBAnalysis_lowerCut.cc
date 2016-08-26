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

using namespace std;
int main(int argc, char** argv)
{

  if(argc!=3)
    {
    cout <<"argc: "<< argc<<" need file name and numPions.." <<endl;
    exit(0);
    }

  char* filename=argv[1];
  int numPions=atoi(argv[2]);
  cout <<"looking at filename; " << filename << " and " << numPions<< " pions " <<endl;

  //FOM S/sqrt(S+B)
  TFile file(filename);
  TTree* tree = (TTree*)file.Get("DataTree");
  int leptonId=0;
  //let's use the boundaries of -0.4 to 1.0 for the n=1 case and -0.4 to 0.4 for n=2

  float lowerCut=-2.0;
  float upperCut=1.0;
  if(numPions==2)
    {
      upperCut=0.4;
    }



  char corrBuffer[]="PIDCorrection*D_DecayCorr*B_DecayCorr*";
  char histoName[200];
  char buffer[3000];
  char bufferBg[3000];
  char bufferSignal[3000];
  char drawCommand[3000];
  int numBins=50;

  int iF=0;
  int b=0;


  float x[300];
  float yFoM[300];
  float yS[300];
  float yB[300];


  float lowerPCut=2.5;
  float cutStep=0.05;
  int numSteps=40;

    sprintf(histoName,"histo_If_%d_b_%d",iF,b);

    for(int i=0; i<numSteps;i++)
      {
	//	float pCut=upperPCut+i*pCutStep;
	//	float probCut=lowerPCut+i*cutStep;
	//	float mDnCut=lowerPCut+i*cutStep;
	//	lowerCut+=cutStep;
	lowerCut+=cutStep;
	cout <<"testing cut " << upperCut <<endl;
	sprintf(buffer,"%s tagCorr*(bestD && pi1Mom>0.24&& recDecaySignature  &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge))  ",corrBuffer,upperCut,lowerCut,numPions);
	cout <<"first buffer: "<< buffer <<endl;
	sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins,lowerCut,upperCut);
	if(numPions==1)
	  {
	    sprintf(bufferBg,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu)",buffer);
	//fix sig_numPions to 1,because that is our main signal
	    sprintf(bufferSignal,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu",buffer);
	  }
	else
	  {
	    sprintf(bufferBg,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions != 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu)",buffer);
	//fix sig_numPions to 2,because that is our main signal
	    sprintf(bufferSignal,"%s && foundAnyDDoubleStar==1 && sig_numPions==2 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu",buffer);
	  }
	
	int bgCounts=0;
	bgCounts=tree->Draw(drawCommand,(char*)bufferBg);
	//do we need the result histo? ...only for debug purposes...
	//	cout <<" bg selection: " << bufferBg<<endl <<" signal Selection: " << bufferSignal <<endl;
	//  TH1F* result=(TH1F*)gDirectory->Get(histoName);
	int signalCounts=tree->Draw(drawCommand,(char*)bufferSignal);
		cout <<"got " << bgCounts << " bg counts, " << signalCounts <<" signal counts, FoM: " << signalCounts/(double)(sqrt(signalCounts+bgCounts))<<endl;
	float fom=0;
	if(signalCounts+bgCounts>0)
	  fom=(float)signalCounts/(float)sqrt(signalCounts+bgCounts);
	x[i]=lowerCut;
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
    grFoM.GetXaxis()->SetTitle("lower m_{#nu}^{2} [GeV] cut  ");
    cout <<"drawing.." <<endl;
    grFoM.Draw("AP*");
    if(numPions==1)
      c.SaveAs("fomLowerCut1.png");
    else
      c.SaveAs("fomLowerCut2.png");

}
