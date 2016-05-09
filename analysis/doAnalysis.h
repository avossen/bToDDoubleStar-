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
#include "TFractionFitter.h"


using namespace std;

int numBinsSB=20;
int numBinsWS=20;


//int numBins=40;
//float lowerCut=-0.3;
//float upperCut=0.5;

int numBins=100;
float lowerCut=-1.0;
float upperCut=1.0;

bool withPIDCorrection;
bool withLumiCorrection;
bool withBCorrection;
bool withDCorrection;


TColor* glColorTable[11];
void addCorrections(char* buffer)
{
  sprintf(buffer,"");
  //D_DecayCorr*B_DecayCorr*PIDCorrection*CrossSectionLumiCorrection

  if(withPIDCorrection)
    {
      sprintf(buffer,"PIDCorrection*");
    }
  if(withDCorrection)
    {
      cout <<" adding d correction: " << buffer <<endl;
      sprintf(buffer,"D_DecayCorr* %s ",buffer);
      cout <<"fuffer is now: "<< buffer <<endl;
    }
  if(withBCorrection)
    {
      cout <<" adding B correction: " << buffer <<endl;
      sprintf(buffer,"B_DecayCorr* %s ",buffer);
      cout <<"buffer is now: "<< buffer <<endl;
    }
  if(withLumiCorrection)
    {
      sprintf(buffer,"CrossSectionLumiCorrection* %s",buffer);
    }
};

/*
get the histogram we want to fit from the tree (given numPions and leptonId), fit with the fractions from 'summedComponents'

 */
void fitFractions(TTree** trees, TH1F** summedComponents, int numComponents,int numPions,int leptonId,bool dataTree)
{
  char histoName[2009];
  char drawCommand[2000];
  char buffer[2000];



  int minCounts=10;

  ////from the example on the root web pages..

    TH1F *data;                              //data histogram
    int treeCount=4;
    if(dataTree)
      {
	treeCount=1;
      }
    if(leptonId!=0)
	  {
	    addCorrections(buffer);
	    sprintf(buffer,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",buffer,upperCut,lowerCut,numPions,leptonId);
		  //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
	  }
	else
	  {
	    addCorrections(buffer);
	    sprintf(buffer,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge==((-1)*systemCharge) ) ",buffer,upperCut,lowerCut,numPions);
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
	  }

    for(int tc=0;tc<treeCount;tc++)
      {
	sprintf(histoName,"histo_Data_%d_pions_%d_leptonId_treeNr%d",numPions,leptonId,tc);
	sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins,lowerCut,upperCut);
	int counts=trees[tc]->Draw(drawCommand,(char*)buffer);
	cout <<"got " << counts <<" counts from data selected " <<endl;
	if(tc==0)
	  data=(TH1F*)gDirectory->Get(histoName);
	else
	  data->Add((TH1F*)gDirectory->Get(histoName));

	data->SetFillColor(glColorTable[0]->GetNumber());	
      }
    TCanvas sampleData;
    data->Draw();
    sampleData.SaveAs("sampleData.png");

  int numEffectiveComponents=0;
  for(int i=0;i<numComponents;i++)
    {
      if(summedComponents[i]->GetEntries()>minCounts)
	{
	  numEffectiveComponents++;
	}
    }
  cout <<" we have " << numEffectiveComponents<<endl;
  TObjArray *mc = new TObjArray(numEffectiveComponents);        // MC histograms are put in this array

  for(int i=0;i<numComponents;i++)
    {
      cout <<"summed compoentent " << i <<" has " << summedComponents[i]->GetEntries() <<" entries" <<endl;
      if(summedComponents[i]->GetEntries()>minCounts)
	{
	  mc->Add(summedComponents[i]);
	}
      sprintf(buffer,"mcComponent_%d.png",i);
      summedComponents[i]->Draw();
      sampleData.SaveAs(buffer);
    }

  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  for(int i=0;i<numEffectiveComponents;i++)
    {
            fit->Constrain(i,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
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
  }


  ////




}



//trees: the input trees for the 4 MC files, components: The components for each file, so 4*11, summedComponents: The same, but summed over files
void getMCComponents(TTree** trees, TH1F** components, TH1F** summedComponents, int numPions,int leptonId)
{
  char histoName[2009];
  char drawCommand[2000];
  //just getting the MC components, so only looking at 4 files
  int numFiles=4;
  int numComponents=11;

  

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

  TH1F** histos[4];
  for(int i=0;i<4;i++)
    {
      histos[i]=new TH1F*[7];
    }

  char* selections[11];



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




  if(leptonId!=0)
    {
        addCorrections(buffer);
      sprintf(buffer,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  ",buffer,upperCut,lowerCut,numPions,leptonId);
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d  ",numPions,leptonId);

    }
  else
    {
      addCorrections(buffer);
      sprintf(buffer,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge==((-1)*systemCharge) ",buffer,upperCut,lowerCut,numPions);
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

  sprintf(bufferAll,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu|| sig_DStarPiPiLNu) && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);

  sprintf(bufferNoSelection,"%s)",buffer);
  //    sprintf(bufferAll,"%s)",buffer);


  cout <<"selection for DDstar: "<< bufferDDStar <<endl<<" DDstarPi: "<< bufferDDStarPi<<endl<< " DDstarPiPi: "<< bufferDDStarPiPi <<endl;
  cout <<" DlNu: " << bufferDlNu<<endl <<" DPilNu:"<< bufferDPilNu << endl << " DPiPilNu: "<< bufferDPiPilNu <<endl;
  cout <<" all: "<< bufferAll <<endl;



  selections[0]=bufferNoSelection;
  selections[1]=bufferDDStar;
  selections[2]=bufferDDStarPi;
  selections[3]=bufferDDStarPiPi;
  selections[4]=bufferDlNu;
  selections[5]=bufferDPilNu;
  selections[6]=bufferDPiPilNu;

  selections[7]=bufferDStarlNu;
  selections[8]=bufferDStarPilNu;
  selections[9]=bufferDStarPiPilNu;
  selections[10]=bufferAll;


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
  summedHistos[5]=hDPilNu;
  summedHistos[6]=hDPiPilNu;
  summedHistos[7]=hDStarlNu;
  summedHistos[8]=hDStarPilNu;
  summedHistos[9]=hDStarPiPilNu;
  summedHistos[10]=hOtherBB;

  
  for(int iF=0;iF<numFiles;iF++)
    {
      cout <<" iF: "<< iF <<endl;
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
	  if(b==10)
	    cout <<"all counts so far: "<< allCounts <<" no selection counts: "<< noSelCounts <<" difference: "<< noSelCounts-allCounts <<endl;
	  cout <<"got " << counts <<" counts selected " <<endl;
	  //do we have to clone this before we can return the histo?
	  TH1F* result=(TH1F*)gDirectory->Get(histoName);
	  result->SetFillColor(glColorTable[b]->GetNumber());
	  components[iF*11+b]=result;

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


  int printOrder[]={0,10,1,3,4,5,9,7,8,6,2};
  if(numPions==1)
    {
      printOrder[10]=3;
      printOrder[3]=2;
      printOrder[9]=4;
      printOrder[4]=6;
    }

  cout <<"in save stack.. " <<endl;
  //this is for "all", i.e. not mixed, charged separation
  THStack all;
  char* legendNames[11];
  char* allLegendNames[11];
  int numFiles=4;
  char* fileNames[numFiles];
  for(int i=0;i<numFiles;i++)
    {
      fileNames[i]=new char[200];
    }


  for(int i=0;i<11;i++)
    {
      legendNames[i]=new char [200];
      allLegendNames[i]=new char[200];
    }


  sprintf(legendNames[0]," no Selection");
  sprintf(legendNames[1],"D Double Star (no Dn#pi l#nu)");
  sprintf(legendNames[2],"D Double Star #pi (no Dn#pi l#nu)");
  sprintf(legendNames[3],"D Double Star #pi #pi (no Dn#pi l#nu)");
  sprintf(legendNames[4],"D l #nu");
  sprintf(legendNames[5],"D #pi l #nu");
  sprintf(legendNames[6],"D #pi #pi l #nu");

  sprintf(legendNames[7],"D* l #nu");
  sprintf(legendNames[8],"D* #pi l #nu");
  sprintf(legendNames[9],"D* #pi #pi l #nu");
  sprintf(legendNames[10]," no D(*)n #pi l#nu ");

  sprintf(allLegendNames[0],"Continuum");
  sprintf(allLegendNames[1],"D Double Star");
  sprintf(allLegendNames[2],"D Double Star #pi");
  sprintf(allLegendNames[3],"D Double Star #pi #pi");
  sprintf(allLegendNames[4],"D l #nu");
  sprintf(allLegendNames[5],"D #pi l #nu");
  sprintf(allLegendNames[6],"D #pi #pi l #nu");
  sprintf(allLegendNames[7],"D* l #nu");
  sprintf(allLegendNames[8],"D* #pi l #nu");
  sprintf(allLegendNames[9],"D* #pi #pi l #nu");
  sprintf(allLegendNames[10]," other B B ");


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
      for(int b=0;b<11;b++)
	{
	  int index=printOrder[b];
	  cout <<"b: " << b <<" index: " << index<<endl;
	  sprintf(histoName,"histo_If_%d_index_%d_numPions_%d_leptonId_%d",iF,index,numPions,leptonId);
	  sprintf(outFileName,"%s.png",histoName);
	  //	  components[iF*11+b]=result;
	  cout <<"grabbing component.." << iF*11+index <<endl;
	  TH1F* result=(TH1F*)components[iF*11+index];
	  cout <<"done that .." <<endl;
	  //  if(counts>0)
	    {
	      //other bb is not 10
	      if(index==10)
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
      stacks[iF]->SetTitle(stackName);
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
  for(int i=0;i<11;i++)
    {
      int index=printOrder[i];
      all.Add(summedHistos[index]);
      legend->AddEntry(summedHistos[index],allLegendNames[index],"f");
    }
  TCanvas c;
  sprintf(buffer,"All_%d_pions_%d_leptonId",numPions,leptonId);
  all.SetTitle(buffer);
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
  float upperSidebandTop=1.0;
  float upperSidebandBottom=0.5;

  float lowerSidebandTop=-0.5;
  float lowerSidebandBottom=-1.0;

  char upperSBSelection[2000];
  char lowerSBSelection[2000];

  char histoName[200];
  char drawCommand[200];
  char buffer[200];
  //select sidebands from all (need to redo all components because we select a different range
  if(leptonId!=0)
    {
        addCorrections(buffer);
      sprintf(upperSBSelection,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",buffer, upperSidebandTop,upperSidebandBottom,numPions,leptonId);
        
        addCorrections(buffer);
      sprintf(lowerSBSelection,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d)  ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions,leptonId);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
        addCorrections(buffer);
      sprintf(upperSBSelection,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge==((-1)*systemCharge) ) ",buffer, upperSidebandTop,upperSidebandBottom,numPions);
        
        addCorrections(buffer);
      sprintf(lowerSBSelection,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge==((-1)*systemCharge) ) ",buffer, lowerSidebandTop,lowerSidebandBottom,numPions);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
    }
  sprintf(histoName,"upperSideBandMC_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSB,upperSidebandBottom,upperSidebandTop);
  int counts=mcTree->Draw(drawCommand,(char*)upperSBSelection);
  cout <<"got " << counts <<" counts from mc upperSB selected " <<endl;
  *upperSidebandMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"lowerSideBandMC_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSB,lowerSidebandBottom,lowerSidebandTop);
   counts=mcTree->Draw(drawCommand,(char*)lowerSBSelection);
  cout <<"got " << counts <<" counts from mc lowerSB selected " <<endl;
  *lowerSidebandMC=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"upperSideBandData_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSB,upperSidebandBottom,upperSidebandTop);
  cout <<"about to get data with string: "<< drawCommand <<endl;
   counts=dataTree->Draw(drawCommand,(char*)upperSBSelection);
  cout <<"got " << counts <<" counts from data upperSB selected " <<endl;
  *upperSidebandData=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"lowerSideBandData_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsSB,lowerSidebandBottom,lowerSidebandTop);
   counts=dataTree->Draw(drawCommand,(char*)lowerSBSelection);
  cout <<"got " << counts <<" counts from data lowerSB selected " <<endl;
  *lowerSidebandData=(TH1F*)gDirectory->Get(histoName);

}


void doWrongSignComparison(TTree* mcTree,TTree* dataTree, int leptonId, int numPions, TH1F** sameChargeMC, TH1F** chargeNeutralMC,  TH1F** sameChargeData, TH1F** chargeNeutralData)
{
  char sameChargeSelection[2000];
  char chargeNeutralSelection[2000];
  char histoName[200];
  char drawCommand[2000];
  char buffer[200];
  //select sidebands from all (need to redo all components because we select a different range
  if(leptonId!=0)
    {
        addCorrections(buffer);
      sprintf(sameChargeSelection,"%s tagCorr*CrossSectionLumiCorrection*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==systemCharge  && abs(bestBCharge)==1 && abs(leptonId)==%d)  ",buffer, upperCut,lowerCut,numPions,leptonId);
      
        addCorrections(buffer);
        sprintf(chargeNeutralSelection,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge!=systemCharge  && ((bestBCharge==0) || (systemCharge==0)) && abs(leptonId)==%d)  ",buffer, upperCut,lowerCut,numPions,leptonId);
      //sprintf(buffer,"D_DecayCorr*B_DecayCorr*PIDCorrection*tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
	    
    }
  else
    {
         addCorrections(buffer);
      sprintf(sameChargeSelection,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge==systemCharge && abs(bestBCharge)==1 ) ",buffer, upperCut,lowerCut,numPions);
        
        addCorrections(buffer);
      sprintf(chargeNeutralSelection,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge!=systemCharge && ((bestBCharge==0) || (systemCharge==0)) ) ",buffer, upperCut,lowerCut,numPions);
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
  counts=dataTree->Draw(drawCommand,(char*)sameChargeSelection);
  cout <<"got " << counts <<" counts from data same charge selected " <<endl;
  *sameChargeData=(TH1F*)gDirectory->Get(histoName);

  sprintf(histoName,"neutralChargeData_numPions_%d_leptonId_%d",numPions,leptonId);
  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBinsWS,lowerCut,upperCut);
  counts=dataTree->Draw(drawCommand,(char*)chargeNeutralSelection);
  cout <<"got " << counts <<" counts from data charge neutral selected " <<endl;
  *chargeNeutralData=(TH1F*)gDirectory->Get(histoName);

}


#endif

