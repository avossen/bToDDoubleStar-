#ifndef DO_ANALYSIS_H
#define DO_ANALYSIS_H

#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.5 && mDnPi < 3.0  "
#define P2STRING "bestD==1 && pi2Mom < 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb  > -2.6  && mDnPi < 3.0 && !(dType==1 && dDecay==2) && !(dType==1 && dDecay==3) && !(dType==2) " 
//#define P2STRING "recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb  > -3.0  && mDnPi < 3.0 "
//&& !(dType==0 && dDecay==1) && !(dType==1 && dDecay==3) && !(dType==2 && dDecay==1) 


#include <math.h>
#include "TRandom.h"
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
TRandom* rnd;
int numBinsSB=20;
int numBinsSBTop=40;
int numBinsWS=20;
char* allLegendNames[9];

enum selectionNames{};

//int numBins=40;
//float lowerCut=-0.3;
//float upperCut=0.5;

//int numBins=150;
int numBins=50;
//int numBins=100;

//float lowerCut=-0.5;
float lowerCut=-1.0;
float upperCut=2.0;

bool withPIDCorrection;
bool withLumiCorrection;
bool withBCorrection;
bool withDCorrection;
float getFitSignal(TH1F* data,TH1F** templates,int numTemplates, double& fitVal, double& fitErr );
void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents, TH1D* pulls,double signalFraction);
void fillTemplates(TH1F** templates,TH1F** summedComponents,char* templateLegendNames,int numPions);
void addPoissonNoise(TH1F* h);

TColor* glColorTable[9];
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

/*
get the histogram we want to fit from the tree (given numPions and leptonId), fit with the fractions from 'summedComponents'

 */
void fitFractions(TTree** trees, TH1F** summedComponents, int numComponents,int numPions,int leptonId,bool dataTree, bool addNoise, TH1D* pulls)
{

  float templateScaleFactor=1.0;
  if(dataTree)
    templateScaleFactor=0.2;

#ifdef MC_TEST
   templateScaleFactor=0.25;
#endif
   cout <<"using scale factor " << templateScaleFactor <<endl;
   for(int i=0;i<numComponents;i++)
     {
       summedComponents[i]->Sumw2();
       //       (*summedComponents)->SetFillStyle(1001);
       summedComponents[i]->Scale(templateScaleFactor);
     }
   cout <<"using scale factor: "<< templateScaleFactor <<endl;
  TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
  cout <<"fit fractions with " << numComponents <<" components"<<endl;
  char* templateLegendNames[9];
  for(int i=0;i<9;i++)
    {
      templateLegendNames[i]=new char[300];
    }
  vector<int> effectiveComponentsIndices;
  char histoName[2009];
  char drawCommand[2000];
  char buffer[2000];
  char corrBuffer[2000];
  int minCounts=0;
  int fixThresholdCounts=2000;
  if(leptonId>0)
    fixThresholdCounts=1000;
  //set to the pionId 1 to merge components...
  int oneIdx=-1;
  ///make new templates that take into account that some contributions look the same, so probably create problems while fittign
  //same shape combination would be 5 with 7 (only one merger) --> alternatively used 4 mergers 5, 6, 7, 8,9 where all the others have small statistics
  int numMergers=0;
  if(numPions==oneIdx)
    numMergers=1;
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
  //  if(numPions==2)
  else
    {
      for(int i=0;i<numComponents;i++)
	{
	  templates[i]=summedComponents[i];
	  templateLegendNames[i]=allLegendNames[i];
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
    if(leptonId!=0)
	  {
		//for the MC_test we still have to do the tag corr....
	    if(dataTree)
	      {

#ifdef MC_TEST
		if(numPions==1)
		  sprintf(buffer,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);
		else
		  sprintf(buffer,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);

#else
		if(numPions==1)
		  sprintf(buffer,"%s "P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);
		else
		  sprintf(buffer,"%s ("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);

#endif
	      }
	    else
	      {
		if(numPions==1)
		  sprintf(buffer,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);
		else
		  sprintf(buffer,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);
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
		  sprintf(buffer,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) ) ",corrBuffer,upperCut,lowerCut,numPions);
		else
		  sprintf(buffer,"%s tagCorr*("P2STRING"  && bestBCharge==((-1)*systemCharge) ) ",corrBuffer,upperCut,lowerCut,numPions);
#else
		if(numPions==1)
		  sprintf(buffer,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge) ) ",corrBuffer,upperCut,lowerCut,numPions);
		else
		  sprintf(buffer,"%s ("P2STRING"  && bestBCharge==((-1)*systemCharge) ) ",corrBuffer,upperCut,lowerCut,numPions);
#endif
	      }
	    else
	      {
		if(numPions==1)
		  sprintf(buffer,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) ) ",corrBuffer,upperCut,lowerCut,numPions);
		else
		  sprintf(buffer,"%s tagCorr*("P2STRING"  && bestBCharge==((-1)*systemCharge) ) ",corrBuffer,upperCut,lowerCut,numPions);
	      }
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
	  }
    cout <<" treeCount: " << treeCount <<endl;
    for(int tc=0;tc<treeCount;tc++)
      {
	cout << " tc: " << tc<< " numPions: "<< numPions<<" leptId: "<< leptonId <<endl;
	sprintf(histoName,"histo_Data_%d_pions_%d_leptonId_treeNr%d",numPions,leptonId,tc);
	cout <<"draw command..." <<endl;
	sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins,lowerCut,upperCut);
	cout <<" and the actual draw... "<< drawCommand<<"--> >" <<buffer<<endl;
	int counts=0;
	if(!dataTree)
	  {
	    counts=trees[tc]->Draw(drawCommand,(char*)buffer);
	  }
	else
	  {
	    //grab last tree which should be the data tree...
	    counts=trees[4]->Draw(drawCommand,(char*)buffer);
	  }
	cout <<"got " << counts <<" counts from data selected " <<endl;
	if(tc==0)
	  data=(TH1F*)gDirectory->Get(histoName);
	else
	  data->Add((TH1F*)gDirectory->Get(histoName));

	data->SetFillColor(glColorTable[0]->GetNumber());	
      }

    cout <<" done making data ..  " << endl;
    double signalFraction=templates[2]->Integral()/data->Integral();
    cout <<"signal Fraction estimated to be : " << signalFraction<<endl;
    /////-----------------
    //void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents)
    if(addNoise)
      {
	for(int nIt=0;nIt<1000;nIt++)
	  {
	    performFractionFitStabilityTest(templates,data,numComponents-numMergers,pulls,signalFraction);
	  }
      }
    /////////-----------



    TCanvas sampleData;
    data->Draw();
    sampleData.SaveAs("sampleData.png");

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
    }
  cout <<" we have " << numEffectiveComponents<<endl;
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
      sprintf(buffer,"mcComponent_%d_leptId_%d_numPion_%d.png",i,leptonId,numPions);
      templates[i]->Draw();
      sampleData.SaveAs(buffer);
    }

  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  fit->GetFitter()->SetPrecision(0.001);
  for(int i=0;i<numEffectiveComponents;i++)
    {
      cout <<" looking at component " << i <<" have " << countsOfComponents[i] << " counts " <<endl;
      //i = 9 (after merger)is BB background
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
	  fit->GetFitter()->FixParameter(i);
	}
    }
      //                  fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
  Int_t status = fit->Fit();               // perform the fit
  std::cout << "fit status: " << status << std::endl;
  if (status == 0) {                       // check on fit status
    TCanvas c;
    TH1F* result = (TH1F*) fit->GetPlot();
    data->Draw("Ep");
    result->Draw("same");
    sprintf(buffer,"fracFit_numPions_%d_leptonId_%d.png",numPions,leptonId);
    c.SaveAs(buffer);
    //and do this for all parameters:
    sprintf(buffer,"fracFitComp_numPions_%d_leptonId_%d.png",numPions,leptonId);
    THStack* predComponents=new THStack(buffer,buffer);
    cout <<"fraction fitter has " << fit->GetFitter()->GetNumberFreeParameters() << " free and " << fit->GetFitter()->GetNumberTotalParameters() <<" overall parameters" <<endl;

    //we want to flip the signal (index =2 ) so that it is later
    int signalIdx=-1;

    for(int i=0;i<numEffectiveComponents;i++)
      {
	if(effectiveComponentsIndices[i]!=2)
	  {
	    cout <<"component " << i<< " (" << templateLegendNames[effectiveComponentsIndices[i] ] <<") is now : "<< fit->GetFitter()->GetParameter(i)<<endl;
	    TH1F* mcComp=(TH1F*) fit->GetMCPrediction(i);
	    cout <<"comp with hist entries: "<< mcComp->GetEntries() <<" and " << templates[effectiveComponentsIndices[i] ]->GetEntries() <<": " << (float)mcComp->GetEntries()/(float)templates[effectiveComponentsIndices[i] ]->GetEntries()<<endl;
	    cout <<"comp with hist scale: "<< (float)mcComp->GetEntries()/(float)data->GetEntries() <<" and " << (float)templates[effectiveComponentsIndices[i] ]->GetEntries()/(float)data->GetEntries()<<endl;
	    //	    mcComp->SetFillColorAlpha(1.0);
	    mcComp->SetFillStyle(1001);
	    mcComp->SetFillColor(glColorTable[effectiveComponentsIndices[i]]->GetNumber());

	    predComponents->Add(mcComp);
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
	mcComp->SetFillStyle(1001);
	mcComp->SetFillColor(glColorTable[2]->GetNumber());
	predComponents->Add(mcComp);
	legend->AddEntry(mcComp,templateLegendNames[2],"f" );
      }


    ////
    if(numPions==1)
      {
	if(leptonId==0)
	  data->GetYaxis()->SetRangeUser(0,6000*templateScaleFactor);
	else
	  data->GetYaxis()->SetRangeUser(0,3000*templateScaleFactor);
      }
    else
      {
	if(leptonId==0)
	  data->GetYaxis()->SetRangeUser(0,1600*templateScaleFactor);
	else
	  data->GetYaxis()->SetRangeUser(0,800*templateScaleFactor);
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
    sprintf(buffer,"predComp_numPions_%d_leptonId_%d.png",numPions,leptonId);
    c.SaveAs(buffer);
  }


  ////

}




void loadComponents(TFile* file,TH1F** components, TH1F** summedComponents, int numPions,int leptonId, int numComponents, int numFiles)
{

  char buffer[2000];
  for(int iF=0;iF<numFiles;iF++)
    {
      for(int b=0;b<numComponents;b++)
	{
	  sprintf(buffer,"histo_If_%d_b_%d_numPions_%d_leptonId_%d",iF,b,numPions,leptonId);
	  components[iF*9+b]=(TH1F*)file->Get(buffer);
	}
    }


  sprintf(buffer,"continuum_%dLept_%dPions",leptonId,numPions);
  summedComponents[0]=(TH1F*)file->Get(buffer);

  cout <<"pointer : "<< summedComponents[0]<<endl;
  sprintf(buffer,"DDStar_%dLept_%dPions",leptonId,numPions);
  summedComponents[1]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPi_%dLept_%dPions",leptonId,numPions);
  summedComponents[2]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DDStarPiPi_%dLept_%dPions",leptonId,numPions);
  summedComponents[3]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DLnu_%dLept_%dPions",leptonId,numPions);
  summedComponents[4]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[5]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DPiPiLNu_%dLept_%dPions",leptonId,numPions);
  summedComponents[5]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DStarLNu_%dLept_%dPions",leptonId,numPions);
  summedComponents[6]=(TH1F*)file->Get(buffer);

  //  sprintf(buffer,"DStarPiLNu_%dLept_%dPions",leptonId,numPions);
  //  summedComponents[8]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"DStarPiPiLNu_%dLept_%dPions",leptonId,numPions);
  summedComponents[7]=(TH1F*)file->Get(buffer);

  sprintf(buffer,"OtherBB_%dLept_%dPions",leptonId,numPions);
  summedComponents[8]=(TH1F*)file->Get(buffer);

  for(int i=0;i<9;i++)
    {
      cout <<"loaded component " << i << " with counts " << summedComponents[i]->GetEntries()<<endl;
    }
}

//trees: the input trees for the 4 MC files, components: The components for each file, so 4*11, summedComponents: The same, but summed over files
void getMCComponents(TTree** trees, TH1F** components, TH1F** summedComponents, int numPions,int leptonId)
{

  cout <<"getting MC components..." <<endl;
  char histoName[2009];
  char drawCommand[2000];
  //just getting the MC components, so only looking at 4 files
  int numFiles=4;
  int numComponents=9;

  char buffer[2000];
  sprintf(buffer,"continuum_%dLept_%dPions",leptonId,numPions);
  TH1F* hContinuum=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DDStar_%dLept_%dPions",leptonId,numPions);
  TH1F* hDDStar=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DDStarPi_%dLept_%dPions",leptonId,numPions);
  TH1F* hDDStarPi=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DDStarPiPi_%dLept_%dPions",leptonId,numPions);
  TH1F* hDDStarPiPi=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DLnu_%dLept_%dPions",leptonId,numPions);
  TH1F* hDlNu=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DPiLNu_%dLept_%dPions",leptonId,numPions);
  TH1F* hDPilNu=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DPiPiLNu_%dLept_%dPions",leptonId,numPions);
  TH1F* hDPiPilNu=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DStarLNu_%dLept_%dPions",leptonId,numPions);
  TH1F* hDStarlNu=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DStarPiLNu_%dLept_%dPions",leptonId,numPions);
  TH1F* hDStarPilNu=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"DStarPiPiLNu_%dLept_%dPions",leptonId,numPions);
  TH1F* hDStarPiPilNu=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);
  sprintf(buffer,"OtherBB_%dLept_%dPions",leptonId,numPions);
  TH1F* hOtherBB=new TH1F(buffer,buffer,numBins,lowerCut,upperCut);


  char* selections[9];

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
  char corrBuffer[2000];
  addCorrections(corrBuffer);

  if(leptonId!=0)
    {
      if(numPions==1)
	sprintf(buffer,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);
      else
 	sprintf(buffer,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  ",corrBuffer,upperCut,lowerCut,numPions,leptonId);
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d  ",numPions,leptonId);
    }
  else
    {
      if(numPions==1)
	sprintf(buffer,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) ",corrBuffer,upperCut,lowerCut,numPions);
      else
	sprintf(buffer,"%s tagCorr*("P2STRING"  && bestBCharge==((-1)*systemCharge) ",corrBuffer,upperCut,lowerCut,numPions);
      //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0   ",numPions);
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
  selections[3]=bufferDDStarPiPi;
  selections[4]=bufferDlNu;
  //  selections[5]=bufferDPilNu;
  selections[5]=bufferDPiPilNu;

  selections[6]=bufferDStarlNu;
  //  selections[8]=bufferDStarPilNu;
  selections[7]=bufferDStarPiPilNu;
  selections[8]=bufferAll;

  //
  // summedHistos: not differentiated between mixed and charged
  //in the end we use something like 		  summedHistos[b+1]->Add(result from selection[b]) (and continuum is treated extra)
  //
  TH1F**  summedHistos=summedComponents;
  summedHistos[0]=hContinuum;
  summedHistos[1]=hDDStar;
  summedHistos[2]=hDDStarPi;
  summedHistos[3]=hDDStarPiPi;
  summedHistos[4]=hDlNu;
  //  summedHistos[5]=hDPilNu;
  summedHistos[5]=hDPiPilNu;
  summedHistos[6]=hDStarlNu;
  //  summedHistos[8]=hDStarPilNu;
  summedHistos[7]=hDStarPiPilNu;
  summedHistos[8]=hOtherBB;
  
  for(int iF=0;iF<numFiles;iF++)
    {
      cout <<" iF: "<< iF <<" numComps: " << numComponents<< " numFiles: "<< numFiles <<endl;
      int allCounts=0;
      int noSelCounts=0;
      for(int b=0;b<numComponents;b++)
	{
	  cout <<"b: "<< b << endl;
	  sprintf(histoName,"histo_If_%d_b_%d_numPions_%d_leptonId_%d",iF,b,numPions,leptonId);
	  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins,lowerCut,upperCut);
	  cout <<"draw command: " << drawCommand << ", selections: " << (char*) selections[b]<<endl;
	  cout <<"drawing tree " << iF <<endl;
	  int counts=trees[iF]->Draw(drawCommand,(char*)selections[b]);
	  if(b>0)
	    allCounts+=counts;
	  if(b==0)
	    noSelCounts=counts;
	  //  if(b==10)
	  if(b==8)
	    cout <<"all counts so far: "<< allCounts <<" no selection counts: "<< noSelCounts <<" difference: "<< noSelCounts-allCounts <<endl;
	  cout <<"got " << counts <<" counts selected " <<endl;
	  //do we have to clone this before we can return the histo?
	  TH1F* result=(TH1F*)gDirectory->Get(histoName);
	  result->SetFillColor(glColorTable[b]->GetNumber());
	  components[iF*9+b]=result;
	  //	  if(b!=10)
	    {
	      //uds or charm, but only use the noSelection
	      if(iF>=2 && b==0)
		{
		  hContinuum->Add(result);
		  hContinuum->SetFillColor(glColorTable[0]->GetNumber());
		}
	      //mixed or charged. The b==0 case is the noSelection case, which only makes sense for the continuum
	      if(iF<2 && b>0)
		{
		  //this is because the selections and the summedHistos are offset by one
		  //
		  summedHistos[b]->Add(result);
		  //for the color, make it consistent to the other stacks
		  summedHistos[b]->SetFillColor(glColorTable[b]->GetNumber());
		}
	    }
	}
    }
  cout <<"done with getMCComponents" <<endl;
};





void saveStack(TH1F** components, TH1F** summedComponents, int numPions, int leptonId)
{

  //  int printOrder[]={0,10,1,3,4,5,9,7,8,6,2};
  int printOrder[]={0,8,1,3,4,7,6,5,2};
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
  char* legendNames[9];

  int numFiles=4;
  char* fileNames[numFiles];
  for(int i=0;i<numFiles;i++)
    {
      fileNames[i]=new char[200];
    }


  for(int i=0;i<9;i++)
    {
      legendNames[i]=new char [300];
      allLegendNames[i]=new char[300];
    }


  sprintf(legendNames[0]," no Selection");
  sprintf(legendNames[1],"B #rightarrow D Double Star X #rightarrow D^{(*)} l #nu (no non-res Dn#pi l#nu)");
  sprintf(legendNames[2],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu (no non-res Dn#pi l#nu)");
  sprintf(legendNames[3],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi #pi l #nu (no non-res Dn#pi l#nu)");
  sprintf(legendNames[4],"D l #nu");
  //  sprintf(legendNames[5],"D #pi l #nu");
  sprintf(legendNames[5],"D #pi #pi l #nu");

  sprintf(legendNames[6],"D* l #nu");
  //  sprintf(legendNames[8],"D* #pi l #nu");
  sprintf(legendNames[7],"D* #pi #pi l #nu");
  sprintf(legendNames[8]," no D(*)n #pi l#nu ");

  sprintf(allLegendNames[0],"Continuum");
  sprintf(allLegendNames[1],"B #rightarrow D Double Star X #rightarrow D^{(*)}  l #nu");
  sprintf(allLegendNames[2],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi l #nu");
  sprintf(allLegendNames[3],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi #pi l #nu");
  sprintf(allLegendNames[4],"D l #nu");
  //  sprintf(allLegendNames[5],"D #pi l #nu");
  sprintf(allLegendNames[5],"D #pi #pi l #nu");
  sprintf(allLegendNames[6],"D* l #nu");
  //  sprintf(allLegendNames[8],"D* #pi l #nu");
  sprintf(allLegendNames[7],"D* #pi #pi l #nu");
  sprintf(allLegendNames[8]," other B B ");


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
      TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
      //      int allCounts=0;
      for(int b=0;b<9;b++)
	{
	  int index=printOrder[b];
	  cout <<"b: " << b <<" index: " << index<<endl;
	  sprintf(histoName,"histo_If_%d_index_%d_numPions_%d_leptonId_%d",iF,index,numPions,leptonId);
	  sprintf(outFileName,"%s.png",histoName);
	  //	  components[iF*9+b]=result;
	  cout <<"grabbing component.." << iF*9+index <<endl;
	  TH1F* result=(TH1F*)components[iF*9+index];
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
	    }
	}
      cout <<"move on..." <<endl;
      char stackName[200];

      TCanvas c2;
      sprintf(stackName,"Stack_%s_%d pions_%d_leptonId",fileNames[iF],numPions,leptonId);
      //      stacks[iF]->SetTitle(stackName);
      stacks[iF]->SetTitle("");
      //      stacks[iF]->SetStats(0);
      stacks[iF]->Draw();
      legend->Draw();
      c2.Update();
      stacks[iF]->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
      c2.Modified();
      sprintf(stackName,"Stack_%s_%d_pions_leptonId_%d.png",fileNames[iF],numPions,leptonId);
      c2.SaveAs(stackName);
    }


  ///the summed, all stuff...
  TH1F**  summedHistos=summedComponents;
  TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
  for(int i=0;i<9;i++)
    {
      int index=printOrder[i];
      all.Add(summedHistos[index]);
      legend->AddEntry(summedHistos[index],allLegendNames[index],"f");
    }
  TCanvas c;
  sprintf(buffer,"All_%d_pions_%d_leptonId",numPions,leptonId);
  //  all.SetTitle(buffer);
  all.SetTitle("");
  //  all.SetStats(0);
  all.Draw();
  legend->Draw();
  c.Update();
  all.GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
  c.Modified();
  sprintf(buffer,"All_%d_pions_%d_leptonId.png",numPions,leptonId);
  c.SaveAs(buffer);
  sprintf(buffer,"All_%d_pions_%d_leptonId.root",numPions,leptonId);
  c.SaveAs(buffer);



}


void doSidebandComparison(TTree* mcTree, TTree* dataTree,int leptonId, int numPions, TH1F** lowerSidebandMC, TH1F** upperSidebandMC, TH1F** lowerSidebandData, TH1F** upperSidebandData)
{
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
  addCorrections(buffer);
  //select sidebands from all (need to redo all components because we select a different range
  if(leptonId!=0)
    {
      if(numPions==1)
	sprintf(upperSBSelection,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId);
      else
	sprintf(upperSBSelection,"%s tagCorr*("P2STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId);


      if(numPions==1)
	sprintf(upperSBSelectionData,P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ", upperSidebandTop,upperSidebandBottom,numPions,leptonId);
      else
	sprintf(upperSBSelectionData," ("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ", upperSidebandTop,upperSidebandBottom,numPions,leptonId);
     
      if(numPions==1)
	sprintf(lowerSBSelection,"%s tagCorr*"P1STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId);
      else
	sprintf(lowerSBSelection,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId);

      if(numPions==1)
	sprintf(lowerSBSelectionData,P1STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",lowerSidebandTop,lowerSidebandBottom,numPions,leptonId);
      else
	sprintf(lowerSBSelectionData," ("P2STRING"&& bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",lowerSidebandTop,lowerSidebandBottom,numPions,leptonId);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
      if(numPions==1)
	sprintf(upperSBSelection,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) ) ",buffer, upperSidebandTop,upperSidebandBottom,numPions);
      else
	sprintf(upperSBSelection,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) ) ",buffer, upperSidebandTop,upperSidebandBottom,numPions);

      if(numPions==1)
	sprintf(upperSBSelectionData,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge) ) ",buffer, upperSidebandTop,upperSidebandBottom,numPions);
      else
	sprintf(upperSBSelectionData,"%s ("P2STRING" && bestBCharge==((-1)*systemCharge) ) ",buffer, upperSidebandTop,upperSidebandBottom,numPions);
        
      if(numPions==1)
      sprintf(lowerSBSelection,"%s tagCorr*"P1STRING"  && bestBCharge==((-1)*systemCharge) ) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions);
      else
      sprintf(lowerSBSelection,"%s tagCorr*("P2STRING" && bestBCharge==((-1)*systemCharge) ) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions);

      if(numPions==1)
	sprintf(lowerSBSelectionData,"%s "P1STRING"  && bestBCharge==((-1)*systemCharge) ) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions);
      else
	sprintf(lowerSBSelectionData,"%s ("P2STRING" && bestBCharge==((-1)*systemCharge) ) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions);

      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
    }
  cout <<"upper selection: "<< upperSBSelection << " lower selection : "<< lowerSBSelection<<endl;
  sprintf(histoName,"upperSideBandMC_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSBTop,upperSidebandBottom,upperSidebandTop);
  int counts=mcTree->Draw(drawCommand,(char*)upperSBSelection);
  cout <<"got " << counts <<" counts from mc upperSB selected for lepton ID "<< leptonId <<endl;
  *upperSidebandMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"lowerSideBandMC_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSB,lowerSidebandBottom,lowerSidebandTop);
  counts=mcTree->Draw(drawCommand,(char*)lowerSBSelection);
  cout <<"got " << counts <<" counts from mc lowerSB selected " <<endl;
  *lowerSidebandMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"upperSideBandData_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSBTop,upperSidebandBottom,upperSidebandTop);
  cout <<"about to get data with string: "<< drawCommand <<endl;
  counts=dataTree->Draw(drawCommand,(char*)upperSBSelectionData);
  cout <<"got " << counts <<" counts from data upperSB selected " <<endl;
  *upperSidebandData=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"lowerSideBandData_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSB,lowerSidebandBottom,lowerSidebandTop);
  counts=dataTree->Draw(drawCommand,(char*)lowerSBSelectionData);
  cout <<"got " << counts <<" counts from data lowerSB selected " <<endl;
  *lowerSidebandData=(TH1F*)gDirectory->Get(histoName);
}

void doWrongSignComparison(TTree* mcTree,TTree* dataTree, int leptonId, int numPions, TH1F** sameChargeMC, TH1F** chargeNeutralMC,  TH1F** sameChargeData, TH1F** chargeNeutralData)
{
  char sameChargeSelection[2000];
  char chargeNeutralSelection[2000];

  //since we put the tagCorr in there, we have to have a string w/o for the data
  char sameChargeSelectionData[2000];
  char chargeNeutralSelectionData[2000];
  char histoName[200];
  char drawCommand[2000];
  char buffer[2000];
  addCorrections(buffer);
  //select sidebands from all (need to redo all components because we select a different range
  if(leptonId!=0)
    {

      if(numPions==1)
	sprintf(sameChargeSelection,"%s tagCorr*"P1STRING" && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d)  ",buffer, upperCut,lowerCut,numPions,leptonId);
      else
	sprintf(sameChargeSelection,"%s tagCorr*("P2STRING"&& bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d)  ",buffer, upperCut,lowerCut,numPions,leptonId);
      

      if(numPions==1)
        sprintf(chargeNeutralSelection,"%s tagCorr*"P1STRING" && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d)  ",buffer, upperCut,lowerCut,numPions,leptonId);
      else
        sprintf(chargeNeutralSelection,"%s tagCorr*("P2STRING"&& bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d)  ",buffer, upperCut,lowerCut,numPions,leptonId);

      if(numPions==1)
      sprintf(sameChargeSelectionData," "P1STRING" && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d)  ", upperCut,lowerCut,numPions,leptonId);
      else
      sprintf(sameChargeSelectionData," ("P2STRING"&& bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d)  ", upperCut,lowerCut,numPions,leptonId);
      
      if(numPions==1)
        sprintf(chargeNeutralSelectionData," "P1STRING" && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d)  ",upperCut,lowerCut,numPions,leptonId);
      else
        sprintf(chargeNeutralSelectionData," ("P2STRING"&& bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d)  ",upperCut,lowerCut,numPions,leptonId);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
      if(numPions==1)
	sprintf(sameChargeSelection,"%s tagCorr*"P1STRING"  && bestBCharge==systemCharge && abs(bestBCharge)==1 ) ",buffer, upperCut,lowerCut,numPions);
      else
	sprintf(sameChargeSelection,"%s tagCorr*("P2STRING" && bestBCharge==systemCharge && abs(bestBCharge)==1 ) ",buffer, upperCut,lowerCut,numPions);

      if(numPions==1)
	sprintf(chargeNeutralSelection,"%s tagCorr*"P1STRING"  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0)) ) ",buffer, upperCut,lowerCut,numPions);
      else
	sprintf(chargeNeutralSelection,"%s tagCorr*("P2STRING" && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0)) ) ",buffer, upperCut,lowerCut,numPions);


      if(numPions==1)
	sprintf(sameChargeSelectionData," "P1STRING"  && bestBCharge==systemCharge && abs(bestBCharge)==1 ) ",upperCut,lowerCut,numPions);
      else
	sprintf(sameChargeSelectionData," ("P2STRING" && bestBCharge==systemCharge && abs(bestBCharge)==1 ) ",upperCut,lowerCut,numPions);

      if(numPions==1)
	sprintf(chargeNeutralSelectionData," "P1STRING"  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0)) ) ",upperCut,lowerCut,numPions);
      else
	sprintf(chargeNeutralSelectionData," ("P2STRING" && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0)) ) ",upperCut,lowerCut,numPions);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
    }
  sprintf(histoName,"sameChargeMC_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut,upperCut);
  int counts=mcTree->Draw(drawCommand,(char*)sameChargeSelection);
  cout <<"got " << counts <<" counts from mc same charge selected " <<endl;
  *sameChargeMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"neutralChargeMC_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut,upperCut);
  counts=mcTree->Draw(drawCommand,(char*)chargeNeutralSelection);
  cout <<"got " << counts <<" counts from mc charge neutral selected " <<endl;
  *chargeNeutralMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"sameChargeData_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut,upperCut);
  counts=dataTree->Draw(drawCommand,(char*)sameChargeSelectionData);
  cout <<"got " << counts <<" counts from data same charge selected " <<endl;
  *sameChargeData=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"neutralChargeData_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut,upperCut);
  counts=dataTree->Draw(drawCommand,(char*)chargeNeutralSelectionData);
  cout <<"got " << counts <<" counts from data charge neutral selected " <<endl;
  *chargeNeutralData=(TH1F*)gDirectory->Get(histoName);

}

void addPoissonNoise(TH1F* h)
{
  for(int i=1;i<=h->GetNbinsX();i++)
    {
      double bc=h->GetBinContent(i);
      double nv=rnd->Poisson(bc);

      //      cout <<"for bin  " << i <<" old value : "<< bc <<" new : "<< nv <<endl;
      //      cout <<"old error: " << h->GetBinError(i);
      if(nv>0)
	h->SetBinContent(i,nv);
      //      cout <<"new error: "<< h->GetBinError(i) <<endl;;
    }
}

void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents, TH1D* pulls,double signalFraction)
{
  cout <<"trying stability test.. " << endl;
  TH1F** templates=new TH1F*[numComponents];
  float S=0;
  cout <<"doing the templates .." <<endl;
  char buffer[500];
  int iteration=0;

  for(int i=0;i<numComponents;i++)
    {
      sprintf(buffer,"%s_Clone_it%d",templatesOrg[i]->GetName(),iteration);
      templates[i]=(TH1F*)templatesOrg[i]->Clone(buffer);
      cout <<"before template " << i <<" has " << templates[i]->Integral() <<" int and " << templates[i]->GetEntries()<<endl;
      if(templates[i]->GetEntries()>0)
	{
	  addPoissonNoise(templates[i]);
	}
      //      cout <<"template " << i <<" has " << templates[i]->Integral() <<" int and " << templates[i]->GetEntries()<<endl;
    }
  float SData=templates[2]->Integral();


  sprintf(buffer,"%s_Clone_it%d",dataOrg->GetName(),iteration);
  TH1F* data=(TH1F*)dataOrg->Clone(buffer);
  cout <<"doing the data " << endl;

    addPoissonNoise(data);

  cout <<"data has " <<SData <<" counts"<<endl;
  double fitVal=0;
  double fitErr=0;
  S=getFitSignal(data,templates,numComponents, fitVal, fitErr);  
  cout <<"data integral: "<< data->Integral()<< " entries: " << data->GetEntries()<<endl;
  cout <<"fitVal: " << fitVal <<" fitErr: "<< fitErr <<endl;
  cout <<"S is : " << S << " fraction times data: " << data->Integral()*fitVal <<" or " << data->GetEntries()*fitVal<<endl;
  if(S>0&& fitErr>0.0)
    {
      cout <<"pull: " << (fitVal-signalFraction)/fitErr <<endl;
      pulls->Fill((fitVal-signalFraction)/fitErr);
    }
  float errData=1.0;
  float err=1.0;
  if(SData>0)
    errData=sqrt(SData);
  if(S>0)
    err=sqrt(S);
  
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



float getFitSignal(TH1F* data,TH1F** templates,int numTemplates, double& fitVal, double& fitErr)
{
  float S=0;
  vector<int> countsOfComponents;
  vector<float> integralOfComponents;
  vector<int> effectiveComponentsIndices; 
  int numEffectiveComponents=0;
  TObjArray *mc = new TObjArray(numTemplates);        // MC histograms are put in this array
  int signalIndex=0;
  char buffer[300];
  TCanvas c;

  double templateIntegral=0.0;

  for(int i=0;i<numTemplates;i++)
    {
      if(templates[i]->GetEntries()>0)
	{
	  cout <<"template " << i << " has " << templates[i]->GetEntries()<<endl;
	  numEffectiveComponents++;

	  mc->Add(templates[i]);
	  sprintf(buffer,"poissonTemplate%d.png",i);
	  templates[i]->Draw();
	  c.SaveAs(buffer);
	  if(i==2)
	    {
	      signalIndex=effectiveComponentsIndices.size();
	      cout <<"set signalIndex to " << signalIndex<<endl;
	    }
	  effectiveComponentsIndices.push_back(i);
	  //the same as the test effectiveComponentIndices[i]==2 from the regular fit
	  countsOfComponents.push_back(templates[i]->GetEntries());
	  double tempIntegral=templates[i]->Integral();
	  integralOfComponents.push_back(tempIntegral);
	  templateIntegral+=tempIntegral;
	}
    }
  data->Draw();
  c.SaveAs("poissonData.png");
  cout <<"overall integral of templates: " << templateIntegral <<" of data: "<< data->Integral()<<endl;
  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  for(int i=0;i<numEffectiveComponents;i++)
    {
      //use the leptonId=0
      //  int fixThresholdCounts=2000;
      if(countsOfComponents[i]>2000)
	{
	  fit->Constrain(i,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	}
      else
	{
	  double iThFraction=integralOfComponents[i]/templateIntegral;
	  cout <<"fixing parameter " << i <<" to " << iThFraction <<endl;
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
    //we want to flip the signal (index =2 ) so that it is later
    //	if(effectiveComponentsIndices[i]!=2)
      cout <<"getting MC pred " << signalIndex <<endl;
      TH1F* mcComp=(TH1F*) fit->GetMCPrediction(signalIndex);
      S=mcComp->Integral();
      fit->GetResult(signalIndex,fitVal,fitErr);
      for(int i=0;i<numEffectiveComponents;i++)
	{
	  cout <<"effective component " << i<< " has " << fit->GetMCPrediction(i)->Integral() <<" integral and  " << fit->GetMCPrediction(i)->GetEntries() <<" counts" <<endl;
	}
    }
  return S;
}
void fillTemplates(TH1F** templates,TH1F** summedComponents,char* templateLegendNames,int numPions)
{


};
#endif

