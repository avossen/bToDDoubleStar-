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
  if(argc!=2)
    {
    cout <<"argc: "<< argc<<" need onnly file name.." <<endl;
    exit(0);
    }

  char* filename=argv[1];
  cout <<"looking at filename; " << filename <<endl;
  //FOM S/sqrt(S+B)
  TFile file(filename);
  TTree* tree = (TTree*)file.Get("DataTree");
  int numPions=2;
  int leptonId=0;
  //let's use the boundaries of -0.4 to 1.0 for the n=1 case and -0.4 to 0.4 for n=2

  float lowerCut=-0.4;
  float upperCut=1.0;
  if(numPions==2)
    {
      upperCut=0.4;
    }

  int mNu2Counter;
  float mNu2[100];
  int recDecaySignature;
  float pi1Mom[100];
  float pi2Mom[100];
  int numRecPions[100];
  float mBTag;
  float deltaETag;
  float logProb;
  float mDnPi[100];
  int bestD[100];
  int bestBCharge;
  int systemCharge[100];

  int foundAnyDDoubleStar;
  int sig_numKaons;
  int sig_numPi0;
  int sig_numBaryons;
  int sig_numLeptons;
  int sig_numPions;
  int sig_DLNu;
  int sig_DPiLNu;
  int sig_DStarLNu;
  int sig_DStarPiLNu;

  float PIDCorrection;
  float D_DecayCorr;
  float B_DecayCorr;
  float tagCorr;


  tree->SetBranchAddress("bestD",bestD);
  tree->SetBranchAddress("mNu2",mNu2);
  tree->SetBranchAddress("mNu2Counter",&mNu2Counter);
  tree->SetBranchAddress("recDecaySignature",&recDecaySignature);
  tree->SetBranchAddress("pi1Mom",pi1Mom);
  tree->SetBranchAddress("pi2Mom",pi2Mom);
  tree->SetBranchAddress("numRecPions",numRecPions);
  tree->SetBranchAddress("mBTag",&mBTag);
  tree->SetBranchAddress("deltaETag",&deltaETag);
  tree->SetBranchAddress("logProb",&logProb);
  tree->SetBranchAddress("mDnPi",mDnPi);
  tree->SetBranchAddress("bestBCharge",&bestBCharge);
  tree->SetBranchAddress("systemCharge",systemCharge);
  tree->SetBranchAddress("foundAnyDDoubleStar",&foundAnyDDoubleStar);
  tree->SetBranchAddress("sig_numKaons",&sig_numKaons);
  tree->SetBranchAddress("sig_numPi0",&sig_numPi0);
  tree->SetBranchAddress("sig_numBaryons",&sig_numBaryons);
  tree->SetBranchAddress("sig_numLeptons",&sig_numLeptons);
  tree->SetBranchAddress("sig_numPions",&sig_numPions);
  tree->SetBranchAddress("sig_DLNu",&sig_DLNu);
  tree->SetBranchAddress("sig_DPiLNu",&sig_DPiLNu);
  tree->SetBranchAddress("sig_DStarLNu",&sig_DStarLNu);
  tree->SetBranchAddress("sig_DStarPiLNu",&sig_DStarPiLNu);
  tree->SetBranchAddress("PIDCorrection",&PIDCorrection);
  tree->SetBranchAddress("D_DecayCorr",&D_DecayCorr);
  tree->SetBranchAddress("B_DecayCorr",&B_DecayCorr);
  tree->SetBranchAddress("tagCorr",&tagCorr);

  int nevents=tree->GetEntries();

  int singleDCount=0;
  int multipleDCount=0;
  int bestDCount=0;
  int noBestDCount=0;
  for(long i=0;i<nevents;i++)
    {
      tree->GetEntry(i);
      //general event level cuts
      if(recDecaySignature!=1)
	continue;
      if(mBTag<5.27)
	continue;
      if(fabs(deltaETag)>0.16)
	 continue;
      if(logProb<-3)
	continue;


      int numGoodCombinations=0;
      for(int j=0;j<mNu2Counter;j++)
	{
	  if(bestBCharge!=((-1)*systemCharge[j]))
	    continue;
	  if(mNu2[j]>upperCut || mNu2[j]<lowerCut)
	    continue;
	  if(numRecPions[j]!=1)
	    continue;
	  //	  if(pi2Mom[j]>0.25)
	  //	    continue;
	  numGoodCombinations++;
	  if(bestD[j])
	    bestDCount++;
	  else
	    noBestDCount++;

	}
      if(numGoodCombinations==1)
	{
	  singleDCount++;
	}
      else
	{
	  if(numGoodCombinations>1)
	    multipleDCount++;
	}
      //        sprintf(buffer,"%s tagCorr*(recDecaySignature && pi2Mom<%f &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge))  ",corrBuffer,pCut,upperCut,lowerCut,numPions);
	//	cout <<"first buffer: "<< buffer <<endl;
      //	sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins,lowerCut,upperCut);
      //	sprintf(bufferBg,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu)",buffer);
      	//fix sig_numPions to 1,because that is our main signal
      //	sprintf(bufferSignal,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu",buffer);
    }
  cout <<"singleD Count " << singleDCount<<endl;
  cout <<" multiple D count " << multipleDCount << endl;
  cout <<"best Ds: "<< bestDCount <<" noBest D: " << noBestDCount <<endl;
}

////  char corrBuffer[]="PIDCorrection*D_DecayCorr*B_DecayCorr*";
////  char histoName[200];
////  char buffer[3000];
////  char bufferBg[3000];
////  char bufferSignal[3000];
////  char drawCommand[3000];
////  int numBins=50;
////
////  int iF=0;
////  int b=0;
////
////
////  float x[300];
////  float yFoM[300];
////  float yS[300];
////  float yB[300];
////
////
////  float upperPCut=0.1;
////  float pCutStep=0.01;
////  int numSteps=30;
////
////sprintf(histoName,"histo_If_%d_b_%d",iF,b);
////
////    for(int i=0; i<numSteps;i++)
////      {
////	float pCut=upperPCut+i*pCutStep;
////		cout <<"testing cut " << pCut <<endl;
////        sprintf(buffer,"%s tagCorr*(recDecaySignature && pi2Mom<%f &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge))  ",corrBuffer,pCut,upperCut,lowerCut,numPions);
////	//	cout <<"first buffer: "<< buffer <<endl;
////	sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins,lowerCut,upperCut);
////	sprintf(bufferBg,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu)",buffer);
////	//fix sig_numPions to 1,because that is our main signal
////	sprintf(bufferSignal,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu",buffer);
////	
////	int bgCounts=0;
////	bgCounts=tree->Draw(drawCommand,(char*)bufferBg);
////	//do we need the result histo? ...only for debug purposes...
////	//	cout <<" bg selection: " << bufferBg<<endl <<" signal Selection: " << bufferSignal <<endl;
////	//  TH1F* result=(TH1F*)gDirectory->Get(histoName);
////	int signalCounts=tree->Draw(drawCommand,(char*)bufferSignal);
////	//	cout <<"got " << bgCounts << " bg counts, " << signalCounts <<" signal counts, FoM: " << signalCounts/(double)(sqrt(signalCounts+bgCounts))<<endl;
////	float fom=0;
////	if(signalCounts+bgCounts>0)
////	  fom=(float)signalCounts/(float)sqrt(signalCounts+bgCounts);
////	x[i]=pCut;
////	cout <<"fom: " << fom << " x: " << x[i] << " i: " << i <<endl;
////	
////	yFoM[i]=fom;
////	yS[i]=(float)signalCounts;
////	yB[i]=(float)bgCounts;
////      }
////    cout <<"done, doing graphs...num points: " << numSteps<<endl;
////    TGraph grFoM(numSteps,x,yFoM);
////    cout <<"done wiht graph " <<endl;
////    //    TGraph grFoS(numSteps,x);
////    //    TGraph grFoB(numSteps,x);
////    TCanvas c;
////    cout <<"canvas.." <<endl;
////    grFoM.GetYaxis()->SetTitle("FoM");
////    cout <<"x axis " << endl;
////    grFoM.GetXaxis()->SetTitle("max #pi momentum [GeV]");
////    cout <<"drawing.." <<endl;
////    grFoM.Draw("AP*");
////    c.SaveAs("fom.png");
////
////}
