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
#include <vector>


//This option is used when we read in 1/4th of the MC data as the data tree
//#define MC_TEST
#define GENERATE_SINGLE_STREAM
//open box by 15%
//#define PARTIAL_BOX
//#define  USE_DATA

//#define onlySB


#include "TFractionFitter.h"
#include "doAnalysisCombined.h"

using namespace std;
void doAnalysis(char* fileNameMixed, char* fileNameCharged, char* fileNameUds, char* fileNameCharm, int numPions=2, int leptonId=0);

//get the components that are used for teh TFractionFitter
void getMCComponents(TTree** trees, TH1F** components, TH1F** summedComponents, int numPions, int leptonId,int channel);

void doSidebandComparison(TTree* mcTree, TTree* dataTree,int leptonId, int numPions, TH1F** lowerSidebandMC, TH1F** upperSidebandMC, TH1F** lowerSidebandData, TH1F** upperSidebandData);

void doWrongSignComparison(TTree* mcTree,TTree* dataTree, int leptonId,int numPions, TH1F** sameChargeMC, TH1F** chargeNeutralMC,  TH1F** sameChargeData, TH1F** chargeNeutralData);

void saveStack(TH1F** components, TH1F** summedComponents, int numPions, int leptonId, int channel, bool isCombinedChannel=false);

void addCorrections(char* buffer);

//void fitFractions(TTree* tree, TH1F** summedComponents, int numPions,int leptonId, int channel,  bool dataTree=false, bool addNoise=false);



//add one component to distinguish cross-feed
const int gl_numComponents=11;
const int gl_numFiles=4;

int main(int argc, char** argv)
{
  for(int i=0;i<7;i++)
    {
      gl_signalFraction[i]=-1;
      gl_signalInt[i]=-1;
      gl_dataInt[i]=-1;
    }
#ifdef DPiPi_SEARCH
  SIG_IDX= SIG_IDX_D_PI_PI;
#else
  SIG_IDX =SIG_IDX_D_PI;
#endif
  bool doSBComp=false;
  //  bool doSBComp=true;
  //remember that partial_box changes bins, so have to rerun in between
  //int loadFromFile=false;
     int loadFromFile=true;

  glChannelIdx=0;
  pCount=0;
  rnd=new TRandom3();
  gStyle->SetOptStat(0);
  withPIDCorrection=true;
  withLumiCorrection=true;
  withBCorrection=true;
  withDCorrection=true;
  withFFCorrection=true;

  gStyle->SetOptFit(1111);
  if(argc!=7)
    {
      cout <<"see only " << argc << " arguments, need 7 though.." <<endl;
      exit(0);
    }

  char* fileNameMixed=argv[1];
  char* fileNameCharged=argv[2];
  char* fileNameUds=argv[3];
  char* fileNameCharm=argv[4];
  char* fileNameData=argv[5];

  int numPions=atoi(argv[6]);
  cout <<"want " << numPions <<endl;
  if(numPions==0)
    SIG_IDX=5;
  //outdated...
  //  doAnalysis(fileNameMixed,fileNameCharged,fileNameUds,fileNameCharm, numPions);


  char* fileNames[5];

  fileNames[0]=new char[400];
  fileNames[1]=new char[400];
  fileNames[2]=new char[400];
  fileNames[3]=new char[400];
  fileNames[4]=new char[400];


  TFile* files[5];
  TTree* trees[5];

  files[0]=new TFile(fileNameMixed);
  files[1]=new TFile(fileNameCharged);
  files[2]=new TFile(fileNameUds);
  files[3]=new TFile(fileNameCharm);
  files[4]=new TFile(fileNameData);


  for(int i=0;i<5;i++)
    {
      trees[i] = (TTree*)files[i]->Get("DataTree");
      if(!trees[i])
	cout <<"tree " << i << " is NULL" <<endl;
    }

  //4th tree is the data tree
  int dataTreeSize=trees[4]->GetEntries();
  maxDataTreeSize=dataTreeSize;
#ifdef PARTIAL_BOX
  maxDataTreeSize=partialBoxFraction*maxDataTreeSize;
#endif
  glColorTable[0]=gROOT->GetColor(kBlue);
  glColorTable[1]=gROOT->GetColor(kRed);
  glColorTable[2]=gROOT->GetColor(kYellow);
  glColorTable[3]=gROOT->GetColor(kBlack);
  glColorTable[4]=gROOT->GetColor(kMagenta);
  glColorTable[5]=gROOT->GetColor(kOrange);

  glColorTable[6]=gROOT->GetColor(kSpring);
  glColorTable[7]=gROOT->GetColor(kGray);
  //  glColorTable[8]=gROOT->GetColor(kPink);
  glColorTable[8]=gROOT->GetColor(kCyan);
  //should be some sort of brown
  glColorTable[9]=gROOT->GetColor(28);
  glColorTable[10]=gROOT->GetColor(kWhite);



  //replicate all, but the signal and cross feed channel for the combined ones...
    int counter=10;
    for(int i=0;i<11;i++)
      {
        if(i==SIG_IDX|| (i==iDDStarPiCrossFeed && numPions==0))
	  {
	    continue;
	  }
        counter++;
        glColorTable[counter]=glColorTable[counter-11];
      }


  TH1F* components[11*4];
  //the allocation of the actual histograms is done in the functions
  TH1F* summedComponents[7][3][3][20];

  TH1F** sameChargeMC=new TH1F*;
  TH1F** chargeNeutralMC=new TH1F*;
  TH1F** sameChargeData=new TH1F*;
  TH1F** chargeNeutralData=new TH1F*;

  TH1F** lowerSidebandMC=new TH1F*;
  TH1F** upperSidebandMC=new TH1F*;
  TH1F** lowerSidebandData=new TH1F*;
  TH1F** upperSidebandData=new TH1F*;

  TChain* mcChain=new TChain("DataTree");
  char buffer[200];

  //void getMCComponents(TTree** trees, TH1F** components, int numPions, int leptonId);
  int leptonId=0; //meaning both
  int leptonIds[]={0,11,13}; //meaning both
  mcChain->Add(fileNameMixed);
  mcChain->Add(fileNameCharged);
  mcChain->Add(fileNameUds);
  mcChain->Add(fileNameCharm);
  int pionIndex=0;
  if(numPions==2)
    pionIndex=1;
  //attach to the end to not disturb anything, even if it is out of order
  if(numPions==0)
    pionIndex=2;


  TFile* mF=0;
  char filename[200];
  sprintf(filename,"myFile_numPions%d.root",numPions);
  if(!loadFromFile)
    {
      mF=new TFile(filename,"RECREATE");
    }
  else
    {
      cout <<"do load .." <<endl;
      //READ should be default
      mF=new TFile(filename,"READ");
      cout <<"loading .." <<endl;
    }
  //more just in case, pointers are cheap
  sigSignificance=new TH1D*[20];
  vector<int> channelsUnderConsideration;
  //      channelsUnderConsideration.push_back(-1);
      channelsUnderConsideration.push_back(0);
        channelsUnderConsideration.push_back(1);
	channelsUnderConsideration.push_back(2);
      channelsUnderConsideration.push_back(3);

  // 4&5 are the combined channels where feeddown and D* are fitted simultaneously
  //since 4 needs 0&1 and 5 needs 2&3, these have to be always included
     channelsUnderConsideration.push_back(4);
      channelsUnderConsideration.push_back(5);
  //channel
  //    for(int iC=-1;iC<4;iC++)
      TH1F* data[7];
  for(vector<int>::iterator it=channelsUnderConsideration.begin();it!=channelsUnderConsideration.end();it++)
    {
      cout <<"looking at channel " << (*it) <<endl;
      int iC=(*it);
      glChannelIdx=iC+1;
      if(numPions==2)
	{
	  //these are the non-DStar channels
	  if(iC==0 || iC ==2)
	    {
	      continue;
	    }
	}
      if(numPions==0)
	{
	  if(iC==1 || iC==3)
	    {
	      SIG_IDX=7;
	      //the template 6 is merged, so that the indices shift
	      if(combineDPiPi)
		{
		  SIG_IDX=6;
		}
	    }
	  else
	    SIG_IDX=5;
	}
      char channelBuffer[500];
      getChannelString(iC,channelBuffer);
      sprintf(buffer,"significanceOf_%s",channelBuffer);
      sigSignificance[glChannelIdx]=new TH1D(buffer,buffer,100,0,40);
		      
      cout <<"1" <<endl;
      
      //      for(int i=0;i<3;i++)
      	      for(int i=0;i<1;i++)      
      //	      for(int i=1;i<3;i++)      
	{
	  leptonId=leptonIds[i];
#ifndef onlySB




	  if(!loadFromFile)
	    {
	      //components give the components for each file separately
	      if(iC<4)
		{
		  getMCComponents(trees,components,summedComponents[glChannelIdx][i][pionIndex], numPions,leptonId,iC);
		}
	      else
		{
		  int DChannelIdx=1;
		  if(glChannelIdx==6)
		    {
		      DChannelIdx=3;
		    }
		  int DStarChannelIdx=DChannelIdx+1;
		  combineChannels(summedComponents[DChannelIdx][i][pionIndex],summedComponents[DStarChannelIdx][i][pionIndex],summedComponents[glChannelIdx][i][pionIndex],gl_numComponents,DChannelIdx, numPions);
		}
	    }
	  else
	    {
	      if(iC<4)
		{
		  cout <<"trying to load channel " << iC <<endl;
		  loadComponents(mF,components,summedComponents[glChannelIdx][i][pionIndex],numPions,leptonId,iC,gl_numComponents,gl_numFiles);
		  cout <<"done" <<endl;
		}
	      else
		{
		  //channel index is one more than iC (i.e. starts with 0 and is equal 5 if iC==4)
		  int DChannelIdx=1;
		  if(glChannelIdx==6)
		    DChannelIdx=3;
		  int DStarChannelIdx=DChannelIdx+1;
		  cout <<"combining channel, DChannelIdx: " << DChannelIdx <<endl;
	  for(int iF=0;iF<4;iF++)
	    {
	      for(int b=0;b<11;b++)
		{
		  TH1F* result=(TH1F*)components[iF*11+b];
		  cout <<" checking component  before combine" << result->GetName() <<" with " << result->GetNbinsX() <<", "<< result->GetBinCenter(1) <<" to " << result->GetBinCenter(result->GetNbinsX())<<endl;
		}
	    }
		  combineChannels(summedComponents[DChannelIdx][i][pionIndex],summedComponents[DStarChannelIdx][i][pionIndex],summedComponents[glChannelIdx][i][pionIndex],gl_numComponents,DChannelIdx, numPions);
		}
	    }
	  
	  for(int iF=0;iF<4;iF++)
	    {
	      for(int b=0;b<11;b++)
		{
		  TH1F* result=(TH1F*)components[iF*11+b];
		  cout <<" checking component " << result->GetName() <<" with " << result->GetNbinsX() <<", "<< result->GetBinCenter(1) <<" to " << result->GetBinCenter(result->GetNbinsX())<<endl;
		}
	    }





	  cout <<"done loading " <<endl;
	  if(!loadFromFile)
	    {
	      for(int j=0;j<11;j++)
		{
		  cout <<"saving summed C" << j << endl;
		  summedComponents[glChannelIdx][i][pionIndex][j]->Write();
		  //so that any scaled histo is not saved...
		  summedComponents[glChannelIdx][i][pionIndex][j]->SetDirectory(0);
		}
	      for(int b=0;b<gl_numComponents;b++)
		{
		  cout <<"writing " << b << endl;
		  for(int iF=0;iF<gl_numFiles;iF++)
		    {
		      components[iF*gl_numComponents+b]->Write();
		    }
		}
	      mF->Write();
	    }

	  bool isCombinedChannel=false;
	  if(iC>=4)
	    {
	      isCombinedChannel=true;
	    }
	  else
	    {
	      //make sure that the scaling depending on the data or partial box opening is only happening once
	      //so not twice for the combined channels
	      gl_templateScaleFactor=1.0;
#ifdef USE_DATA
   gl_templateScaleFactor=0.2;
#endif
#ifdef PARTIAL_BOX
   gl_templateScaleFactor=partialBoxFraction/5.0;
#endif

#ifdef MC_TEST
      gl_templateScaleFactor=0.25;
#endif
#ifdef GENERATE_SINGLE_STREAM
      gl_templateScaleFactor=0.2;
#endif
	
	      for(int k=0;k<11;k++)
		{
		  //should not be necessary for the fit and would indeed lead to overestimated uncertainties in the fit
		  //		  summedComponents[glChannelIdx][i][pionIndex][k]->Scale(gl_templateScaleFactor);
		}
	    }




	  int numComponents=11;
	  if(isCombinedChannel)
	    numComponents=20;
	  cout <<" saving stacks " <<endl;
	  for(int k=0;k<numComponents;k++)
	    {
	      cout <<"saving stack of loaded components " << k << " with counts " << summedComponents[glChannelIdx][i][pionIndex][k]->GetEntries()<<endl;
	    }
	  //  cout <<"calling save stack.." <<endl;
	  //have to call 'save Stack' to set e.g. 'allLegendNames'...
	

	  saveStack(components,summedComponents[glChannelIdx][i][pionIndex],numPions,leptonId,iC,isCombinedChannel);
	  //      for(
	  //the 'other BB doesn't seem to be used...'
	  //      fitFractions(trees,summedComponents,10, numPions,leptonId,false);
	  sprintf(buffer,"pulls_numPions%d_%s",numPions,channelBuffer);
	  TH1D* pulls=new TH1D(buffer,buffer,100,-3,3);
	  //this used to be 10 components... I don't understand why, I guess that meanst that the other BB was missing



#ifdef USE_DATA
	  bool dataTree=true;
	  //	  bool addNoise=false;
	  bool addNoise=true;
#else
#ifdef PARTIAL_BOX
	  bool dataTree=true;
	  bool addNoise=false;
	  //	  bool addNoise=true;
#else
#ifdef MC_TEST
	  bool dataTree=true;
	  bool addNoise=true;
	  if(numPions==0)
	    addNoise=false;
#else
	  bool dataTree=false;
	  	  bool addNoise=true;
		  //bool addNoise=false;
#endif
#endif
#endif
	  if(iC<4)
	    {

	      //should just add up template, but do the adding for now to x-check that we get the same
	      if(!dataTree)
		{
		  //careful!, for MC_TEST need to go into the getData routine since the signal fractions are determined there. But dataTree should be set for that
		    getDataFromMC(data[glChannelIdx], iC, numPions, leptonId,summedComponents[glChannelIdx][i][pionIndex],numComponents);
#ifdef GENERATE_SINGLE_STREAM
		    data[glChannelIdx]->Scale(0.2);
		    data[glChannelIdx]->GetSumw2()->Set(0);
#endif
		    //		    cout <<"data has " << data[glChannelIdx]->Integral() <<" integral " <<endl;
		    //		    getData(data[glChannelIdx], dataTree, trees, iC, numPions, leptonId);		   
		    //		    cout <<"from real data data has " << data[glChannelIdx]->Integral() <<" integral " <<endl;
		}
	      else
		{
		    getData(data[glChannelIdx], dataTree, trees, iC, numPions, leptonId);		   
	      }
	    }
	  else
	    {
	      int DChannelIdx=1;
	      if(glChannelIdx==6)
		DChannelIdx=3;
	      int DStarChannelIdx=DChannelIdx+1; 
	      //need to create histo, since it is not created from a tree
	      sprintf(buffer,"histo_Data_%d_pions_%d_leptonId__%s",numPions,leptonId,channelBuffer);
	      float upperCutD=upperCut[DChannelIdx];
	      float lowerCutD=lowerCut[DChannelIdx];
	      int numBinsD=numBins[DChannelIdx];
	      int numBinsDStar=numBins[DStarChannelIdx];
	      float upperCutDStar=upperCut[DStarChannelIdx];
	      float lowerCutDStar=lowerCut[DStarChannelIdx];
	      data[glChannelIdx]=new TH1F(buffer,buffer,numBinsD+numBinsDStar,lowerCutD,upperCutD+(upperCutDStar-lowerCutDStar));
	      addShiftedHistos(data[DChannelIdx],data[DStarChannelIdx],data[glChannelIdx],numBinsD,numBinsDStar);
	      //data[glChannelIdx]=
	    }

	  ///temporary

	  fitFractions(data[glChannelIdx],trees,summedComponents[glChannelIdx][i][pionIndex],numComponents, numPions,leptonId,iC,dataTree,addNoise,pulls);


	  TCanvas cSig;
	  sigSignificance[glChannelIdx]->Draw();
	  sprintf(buffer,"signalSignificance_%s.png",channelBuffer);
	  cSig.SaveAs(buffer);
	  sprintf(buffer,"signalSignificance_%s.pdf",channelBuffer);
	  cSig.SaveAs(buffer);
	  sprintf(buffer,"signalSignificance_%s.eps",channelBuffer);
	  cSig.SaveAs(buffer);


	  //-->this calls fitFraction with addNoise set to true
	  ///////      fitFractions(trees,summedComponents[i][pionIndex],9,numPions,leptonId,false,true,pulls);

	  



	  TCanvas cpulls;
	  pulls->Draw();
	  pulls->Fit("gaus");
	  pulls->Draw();
	  sprintf(buffer,"pulls_numPions%d_%s.png",numPions,channelBuffer);
	  cpulls.SaveAs(buffer);
	  sprintf(buffer,"pulls_numPions%d_%s.pdf",numPions,channelBuffer);
	  cpulls.SaveAs(buffer);
	  sprintf(buffer,"pulls_numPions%d_%s.eps",numPions,channelBuffer);
	  cpulls.SaveAs(buffer);
#endif	  

	  if(doSBComp && iC<4)
	    {
	      doSidebandComparison(mcChain,trees[4],leptonId,numPions,lowerSidebandMC,upperSidebandMC,lowerSidebandData,upperSidebandData,iC);

	      (*upperSidebandData)->Sumw2();
	      (*upperSidebandData)->SetLineColor(kBlack);
	      (*lowerSidebandData)->Sumw2();
	      (*upperSidebandMC)->Sumw2();
	      (*lowerSidebandMC)->Sumw2();
	      (*lowerSidebandData)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*lowerSidebandMC)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*upperSidebandData)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*upperSidebandMC)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*lowerSidebandData)->SetLineWidth(2);
	      (*lowerSidebandData)->SetLineColor(kBlack);
	      (*lowerSidebandMC)->SetLineWidth(2);
	      (*upperSidebandData)->SetLineWidth(2);
	      (*upperSidebandMC)->SetLineWidth(2);
	      (*lowerSidebandData)->SetStats(0);
	      (*lowerSidebandMC)->SetStats(0);
	      (*upperSidebandData)->SetStats(0);
	      (*upperSidebandMC)->SetStats(0);
	      (*lowerSidebandData)->SetTitle("");
	      (*lowerSidebandMC)->SetTitle("");
	      (*upperSidebandData)->SetTitle("");
	      (*upperSidebandMC)->SetTitle("");
	   
	      (*lowerSidebandData)->SetMinimum(0);
	      (*lowerSidebandMC)->SetMinimum(0);
	      (*upperSidebandData)->SetMinimum(0);
	      (*upperSidebandMC)->SetMinimum(0);


	      double scaleUpper=(*upperSidebandMC)->GetEntries()/((double)((*upperSidebandData)->GetEntries()));
	      double scaleLower=(*lowerSidebandMC)->GetEntries()/((double)((*lowerSidebandData)->GetEntries()));

	      //5 streams vs 15% of one stream (not true for the sidebands)
#ifdef PARTIAL_BOX
	      //	      double scaleFactor=5*1.0/partialBoxFraction;
	      double scaleFactor=5*1.0;
	      scaleUpper=scaleFactor;
	      scaleLower=scaleFactor;
#else

#ifdef MC_TEST
	      scaleUpper=4.0;
	      scaleLower=4.0;
#else
	      scaleUpper=5.0;
	      scaleLower=5.0;
#endif
#endif

	      Double_t* zeroErrors=new Double_t[500];
	      for(int i=0;i<500;i++)
		{
		  zeroErrors[i]=0.0;
		}

	      cout <<"scaleUpper: "<< scaleUpper <<" lower: "<< scaleLower <<endl;
	   
	      (*upperSidebandData)->Scale(scaleUpper);
	      (*lowerSidebandData)->Scale(scaleLower);


	      TCanvas c("c","c",0,0,1000,800);
	      c.Divide(2,2);
	      c.cd(1);
	      //      stacks[iF]->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
	   
	      (*lowerSidebandMC)->Draw();
	      c.Update();
	   
	      (*lowerSidebandData)->Draw("SAME");
	   
	      c.cd(2);
	   
	      (*upperSidebandMC)->Draw();
	      c.Update();
	   
	      (*upperSidebandData)->Draw("SAME");
	      c.cd(3);
	      TH1F* lowerSidebandMCDiv=(TH1F*)(*lowerSidebandMC)->Clone("lowerSBMCDiv");
	      lowerSidebandMCDiv->Sumw2();
	      lowerSidebandMCDiv->SetError(zeroErrors);
	      (*lowerSidebandData)->SetError(zeroErrors);
	      lowerSidebandMCDiv->Divide(*lowerSidebandData);

	      lowerSidebandMCDiv->Draw();
	      lowerSidebandMCDiv->Fit("pol0");
	      c.cd(4);

	      TH1F* upperSidebandMCDiv=(TH1F*)(*upperSidebandMC)->Clone("upperSBMCDiv");
	      upperSidebandMCDiv->Sumw2();
	      upperSidebandMCDiv->SetError(zeroErrors);
	      (*upperSidebandData)->SetError(zeroErrors);
	      upperSidebandMCDiv->Divide(*upperSidebandData);

	      upperSidebandMCDiv->Draw();
	      upperSidebandMCDiv->Fit("pol0");
	   
	      sprintf(buffer,"SidebandComparison_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SidebandComparison_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SidebandComparison_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      c.cd(0);
	      (*lowerSidebandMC)->Draw();
	      c.Update();
	      (*lowerSidebandData)->Draw("SAME");
	      sprintf(buffer,"LowerSidebandComparison_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"LowerSidebandComparison_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"LowerSidebandComparison_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      c.cd(0);
	      (*upperSidebandMC)->Draw();
	      c.Update();
	      (*upperSidebandData)->Draw("SAME");
	      sprintf(buffer,"UpperSidebandComparison_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"UpperSidebandComparison_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"LowerSidebandComparison_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	   
	      doWrongSignComparison(mcChain,trees[4],leptonId,numPions,sameChargeMC,chargeNeutralMC,sameChargeData,chargeNeutralData,iC);
	      (*sameChargeMC)->Sumw2();
	      (*sameChargeData)->Sumw2();
	      (*chargeNeutralMC)->Sumw2();
	      (*chargeNeutralData)->Sumw2();


	      (*sameChargeData)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*sameChargeMC)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*chargeNeutralData)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*chargeNeutralMC)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*sameChargeData)->SetLineWidth(2);
	      (*sameChargeMC)->SetLineWidth(2);
	      (*chargeNeutralData)->SetLineWidth(2);
	      (*chargeNeutralMC)->SetLineWidth(2);
	      (*sameChargeData)->SetLineColor(kBlack);
	      (*sameChargeMC)->SetLineColor(kBlue);
	      (*sameChargeData)->SetStats(0);
	      (*sameChargeMC)->SetStats(0);
	      (*chargeNeutralData)->SetStats(0);
	      (*chargeNeutralMC)->SetStats(0);
	      (*sameChargeData)->SetTitle("");
	      (*sameChargeMC)->SetTitle("");
	      (*chargeNeutralData)->SetTitle("");
	      (*chargeNeutralMC)->SetTitle("");
	   
	      (*sameChargeData)->SetMinimum(0);
	      (*sameChargeMC)->SetMinimum(0);
	      (*chargeNeutralData)->SetMinimum(0);
	      (*chargeNeutralMC)->SetMinimum(0);


	      //      double scaleSameCharge=(*sameChargeMC)->GetEntries()/((double)((*sameChargeData)->GetEntries()));
	      //      double scaleChargeNeutral=(*chargeNeutralMC)->GetEntries()/((double)((*chargeNeutralData)->GetEntries()));

	      double scaleSameCharge=5.0;
	      double scaleChargeNeutral=5.0;

	      cout <<"scalesameCharge: "<< scaleSameCharge <<" chargeNeutral: "<< scaleChargeNeutral <<endl;
	   
	      (*sameChargeData)->Scale(scaleSameCharge);
	      (*chargeNeutralData)->Scale(scaleChargeNeutral);

	      (*sameChargeMC)->Draw();
	      c.Update();
	      (*chargeNeutralMC)->SetLineColor(kBlue);
	      (*chargeNeutralData)->SetLineColor(kBlack);
	      (*sameChargeData)->Draw("SAME");
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	   
	      (*chargeNeutralMC)->Draw();
	      c.Update();
	      (*chargeNeutralData)->Draw("SAME");
	      sprintf(buffer,"chargeNeutral_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"chargeNeutral_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"chargeNeutral_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);


	      c.Divide(2,2);
	      c.cd(1);

	      (*sameChargeMC)->Draw();
	      c.Update();
	   
	      (*sameChargeData)->Draw("SAME");
	      c.cd(2);
	   
	      (*chargeNeutralMC)->Draw();
	      c.Update();
	   
	      (*chargeNeutralData)->Draw("SAME");
	      c.cd(3);
	      TH1F* sameChargeMCDiv=(TH1F*)(*sameChargeMC)->Clone("sameChargeMCDiv");
	      sameChargeMCDiv->Sumw2();
	      sameChargeMCDiv->SetError(zeroErrors);
	      (*sameChargeData)->SetError(zeroErrors);
	      sameChargeMCDiv->Divide(*sameChargeData);

	      sameChargeMCDiv->Draw();
	      sameChargeMCDiv->Fit("pol0");
	      c.cd(4);

	      TH1F* chargeNeutralMCDiv=(TH1F*)(*chargeNeutralMC)->Clone("chargeNeutralMCDiv");
	      chargeNeutralMCDiv->Sumw2();
	      chargeNeutralMCDiv->SetError(zeroErrors);
	      (*chargeNeutralData)->SetError(zeroErrors);
	      chargeNeutralMCDiv->Divide(*chargeNeutralData);

	      chargeNeutralMCDiv->Draw();
	      chargeNeutralMCDiv->Fit("pol0");
	      sprintf(buffer,"WrongCharge_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"WrongCharge_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"WrongCharge_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	   
	      c.cd(0);
	      (*sameChargeMC)->Draw();
	      (*sameChargeMC)->SetLineColor(kBlue);
	      (*sameChargeData)->SetLineColor(kBlack);
	      c.Update();
	      (*sameChargeData)->Draw("SAME");
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	   
	      c.cd(0);
	      (*chargeNeutralMC)->Draw();
	      (*chargeNeutralMC)->SetLineColor(kBlue);
	      (*chargeNeutralData)->SetLineColor(kBlack);
	      c.Update();
	      (*chargeNeutralData)->Draw("SAME");
	      sprintf(buffer,"ChargeNeutral_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"ChargeNeutral_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"ChargeNeutral_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	    }//end if(doSB..)

	}
    }

}



///void doAnalysis(char* fileNameMixed, char* fileNameCharged, char* fileNameUds, char* fileNameCharm, int numPions, int leptonId)
//-->see git repository history
/////{
/////  /*
/////    we have 4 data files, the two B MCs, charm and uds. For each we extract DDouble star +0,1,2 Pions, as well as the Dlnu +0,1,2 pions, 
/////    and the rest
/////
/////    We do this for the 0,1,2 pions in the final state
/////
/////    The selection criteria are 
/////    i) The correct final state has to be reconstructed
/////    ii) The correct 
/////
/////  */
/////  char* fileNames[4];
/////  char* legendNames[10];
/////
/////  char* allLegendNames[10];
/////
/////  fileNames[0]=new char[500];
/////  fileNames[1]=new char[500];
/////  fileNames[2]=new char[500];
/////  fileNames[3]=new char[500];
/////
/////  THStack all;
/////  TH1F* hContinuum=new TH1F("continuum","continuum",numBins,lowerCut,upperCut);
/////  TH1F* hDDStar=new TH1F("DDStar","DDstar",numBins,lowerCut,upperCut);
/////  TH1F* hDDStarPi=new TH1F("DDStarPi","DDstarPi",numBins,lowerCut,upperCut);
/////  TH1F* hDDStarPiPi=new TH1F("DDStarPiPi","DDStarPiPi",numBins,lowerCut,upperCut);
/////  TH1F* hDlNu=new TH1F("DLnu","DlNu",numBins,lowerCut,upperCut);
/////  TH1F* hDPilNu=new TH1F("DPiLNu","DPiLNu",numBins,lowerCut,upperCut);
/////  TH1F* hDPiPilNu=new TH1F("DPiPiLNu","DPiPiLNu",numBins,lowerCut,upperCut);
/////  TH1F* hDStarlNu=new TH1F("DStarLNu","DStarLNu",numBins,lowerCut,upperCut);
/////  TH1F* hDStarPilNu=new TH1F("DStarPiLNu","DStarPiLNu",numBins,lowerCut,upperCut);
/////  TH1F* hDStarPiPilNu=new TH1F("DStarPiPiLNu","DStarPiPiLNu",numBins,lowerCut,upperCut);
/////  TH1F* hOtherBB=new TH1F("OtherBB","OtherBB",numBins,lowerCut,upperCut);
/////
/////  TH1F* summedHistos[9];
/////  summedHistos[0]=hContinuum;
/////  summedHistos[1]=hOtherBB;
/////  summedHistos[2]=hDDStar;
/////
/////  summedHistos[3]=hDDStarPiPi;
/////  summedHistos[4]=hDlNu;
/////  //  summedHistos[5]=hDPilNu;
/////  summedHistos[5]=hDPiPilNu;
/////  summedHistos[6]=hDStarlNu;
/////  //  summedHistos[8]=hDStarPilNu;
/////  summedHistos[7]=hDStarPiPilNu;
/////  summedHistos[8]=hDDStarPi;  
/////
/////
/////
/////  for(int i=0;i<9;i++)
/////    {
/////      legendNames[i]=new char [200];
/////      allLegendNames[i]=new char[200];
/////    }
/////
/////  sprintf(legendNames[0],"B #rightarrow D Double Star l #nu(no Dn#pi l#nu non-res)");
/////  sprintf(legendNames[1],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi (no Dn#pi l#nu nonresonant)");
/////  sprintf(legendNames[2],"B #rightarrow Double Star X #rightarrow D^{(*)} #pi #pi (no Dn#pi l#nu nonresonant)");
/////  sprintf(legendNames[3],"D l #nu");
/////  // sprintf(legendNames[4],"D #pi l #nu");
/////  sprintf(legendNames[4],"D #pi #pi l #nu");
/////
/////  sprintf(legendNames[5],"D* l #nu");
/////  //  sprintf(legendNames[7],"D* #pi l #nu");
/////  sprintf(legendNames[6],"D* #pi #pi l #nu");
/////  sprintf(legendNames[7]," no D(*)n #pi l#nu ");
/////  sprintf(legendNames[8]," no Selection");
/////
/////
/////
/////
/////  sprintf(allLegendNames[0],"Continuum");
/////  sprintf(allLegendNames[1]," other B B ");
/////  sprintf(allLegendNames[2],"B #rightarrow D Double Star X #rightarrow D^{(*)}");
/////
/////  sprintf(allLegendNames[3],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi #pi");
/////  sprintf(allLegendNames[4],"D l #nu");
/////  //sprintf(allLegendNames[4],"D #pi l #nu");
/////  sprintf(allLegendNames[5],"D #pi #pi l #nu");
/////
/////  sprintf(allLegendNames[6],"D* l #nu");
/////  //  sprintf(allLegendNames[8],"D* #pi l #nu");
/////  sprintf(allLegendNames[7],"D* #pi #pi l #nu");
/////  sprintf(allLegendNames[8],"B #rightarrow D Double Star X #rightarrow D^{(*)} #pi");
/////
/////
/////  sprintf(fileNames[0],"mixed");
/////  sprintf(fileNames[1],"charged");
/////  sprintf(fileNames[2],"uds");
/////  sprintf(fileNames[3],"charm");
/////
/////  TFile* files[4];
/////  TTree* trees[4];
/////
/////  files[0]=new TFile(fileNameMixed);
/////  files[1]=new TFile(fileNameCharged);
/////  files[2]=new TFile(fileNameUds);
/////  files[3]=new TFile(fileNameCharm);
/////
/////
/////  TH1F** histos[4];
/////  for(int i=0;i<4;i++)
/////    {
/////      histos[i]=new TH1F*[7];
/////    }
/////
/////
/////  for(int i=0;i<4;i++)
/////    {
/////      trees[i] = (TTree*)files[i]->Get("DataTree");
/////      if(!trees[i])
/////	cout <<"tree " << i << " is NULL" <<endl;
/////    }
/////
/////  //  dataTree->Draw(" >> result")
/////  //  TH1F* result=(TH1F*) gDirectory->Get("result");
/////
/////  //selection for the data
/////
/////  //numRecPions  
/////  //  recDecaySignature==true   //recDecay signature found
/////
/////  // get mNu2 for d double star decays
/////
/////
/////  //this is for all the same, the selection of the data. It can be 0, one or two pions
/////
/////  char* selections[9];
/////
/////  char buffer[2000];
/////
/////  char bufferDDStar[2000];
/////  char bufferDDStarPi[2000];
/////  char bufferDDStarPiPi[2000];
/////
/////  char bufferDlNu[2000];
/////  char bufferDPilNu[2000];
/////  char bufferDPiPilNu[2000];
/////
/////  char bufferDStarlNu[2000];
/////  char bufferDStarPilNu[2000];
/////  char bufferDStarPiPilNu[2000];
/////
/////  char bufferAll[2000];
/////  char bufferNoSelection[2000];
/////
/////  selections[0]=bufferDDStar;
/////  selections[1]=bufferDDStarPi;
/////  selections[2]=bufferDDStarPiPi;
/////  selections[3]=bufferDlNu;
/////  //  selections[4]=bufferDPilNu;
/////  selections[4]=bufferDPiPilNu;
/////  selections[5]=bufferDStarlNu;
/////  //  selections[7]=bufferDStarPilNu;
/////  selections[6]=bufferDStarPiPilNu;
/////  selections[7]=bufferAll;
/////  selections[8]=bufferNoSelection;
/////
/////  if(leptonId!=0)
/////    {
/////      //      sprintf(buffer,"tagCorr*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  ",numPions,leptonId);
/////      addCorrections(buffer);
/////      sprintf(buffer,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d  ",buffer,upperCut,lowerCut,numPions,leptonId);
/////    }
/////  else
/////    {
/////      //      sprintf(buffer,"CrossSectionLumiCorrection*tagCorr*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && bestBCharge==((-1)*systemCharge) ",numPions);
/////      addCorrections(buffer);
/////      sprintf(buffer,"%s tagCorr*(mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0   ",buffer,upperCut,lowerCut,numPions);
/////    }
/////  sprintf(bufferDDStar,"%s && foundAnyDDoubleStar==1 && sig_numPions==0 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1 && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
/////  //sprintf(bufferDDStar,"%s && foundAnyDDoubleStar==1)",buffer);
/////  sprintf(bufferDDStarPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
/////  sprintf(bufferDDStarPiPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==2 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1&& !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
/////
/////  sprintf(bufferDlNu,"%s && sig_DLNu )",buffer);
/////  sprintf(bufferDPilNu,"%s && sig_DPiLNu&& !sig_DLNu  )",buffer);
/////  sprintf(bufferDPiPilNu,"%s && sig_DPiPiLNu && !sig_DLNu && !sig_DPiLNu)",buffer);
/////  sprintf(bufferDStarlNu,"%s &&  sig_DStarLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu)",buffer);
/////  sprintf(bufferDStarPilNu,"%s  && sig_DStarPiLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu)",buffer);
/////  sprintf(bufferDStarPiPilNu,"%s  && sig_DStarPiPiLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu)",buffer);
/////
/////  sprintf(bufferAll,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu|| sig_DStarPiPiLNu) && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);
/////
/////  sprintf(bufferNoSelection,"%s)",buffer);
/////  //    sprintf(bufferAll,"%s)",buffer);
/////
/////
/////  cout <<"selection for DDstar: "<< bufferDDStar <<endl<<" DDstarPi: "<< bufferDDStarPi<<endl<< " DDstarPiPi: "<< bufferDDStarPiPi <<endl;
/////  cout <<" DlNu: " << bufferDlNu<<endl <<" DPilNu:"<< bufferDPilNu << endl << " DPiPilNu: "<< bufferDPiPilNu <<endl;
/////  cout <<" all: "<< bufferAll <<endl;
/////
/////
/////  char histoName[2009];
/////  char drawCommand[2000];
/////  char outFileName[2000];
/////
/////  THStack* stacks[4];
/////
/////
/////
/////  TColor* colorTable[9];
/////  colorTable[0]=gROOT->GetColor(kBlue);
/////  colorTable[1]=gROOT->GetColor(kRed);
/////  colorTable[2]=gROOT->GetColor(kCyan);
/////  colorTable[3]=gROOT->GetColor(kMagenta);
/////  //  colorTable[4]=gROOT->GetColor(kGreen);
/////  colorTable[4]=gROOT->GetColor(kBlack);
/////  colorTable[5]=gROOT->GetColor(kWhite);
/////  colorTable[6]=gROOT->GetColor(kGray);
/////  //  colorTable[8]=gROOT->GetColor(kViolet);
/////  colorTable[7]=gROOT->GetColor(kYellow);
/////  colorTable[8]=gROOT->GetColor(kTeal);
/////
/////
/////  for(int iF=0;iF<4;iF++)
/////    {
/////      char bufferF[200];
/////      sprintf(bufferF,"stackF_%d",iF);
/////      stacks[iF]=new THStack(bufferF,bufferF);
/////      TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
/////      int allCounts=0;
/////      for(int b=0;b<9;b++)
/////	{
/////	  sprintf(histoName,"histo_If_%d_b_%d",iF,b);
/////	  sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins,lowerCut,upperCut);
/////	  sprintf(outFileName,"%s.png",histoName);
/////	  cout <<"draw command: " << drawCommand << ", selections: " << (char*) selections[b]<<endl;
/////	  cout <<"drawing tree " << iF <<endl;
/////	  int counts=trees[iF]->Draw(drawCommand,(char*)selections[b]);
/////	  if(b<10)
/////	    allCounts+=counts;
/////	  if(b==10)
/////	    cout <<"all counts so far: "<< allCounts <<" no selection counts: "<< counts <<" difference: "<< counts-allCounts <<endl;
/////	  cout <<"got " << counts <<" counts selected " <<endl;
/////	  TH1F* result=(TH1F*)gDirectory->Get(histoName);
/////	  result->SetFillColor(colorTable[b]->GetNumber());
/////
/////	  if(b!=10)
/////	    {
/////	      if(iF>=2)
/////		{
/////		  hContinuum->Add(result);
/////		  hContinuum->SetFillColor(colorTable[10]->GetNumber());
/////		}
/////		
/////	      if(iF<2)
/////		{
/////		  summedHistos[b+1]->Add(result);
/////		  //for the color, make it consistent to the other stacks
/////		  summedHistos[b]->SetFillColor(colorTable[b]->GetNumber());
/////		}
/////	    }
/////
/////	  if(counts>0)
/////	    {
/////	      //why is the other B,B not in here?
/////	      if(b==10)
/////		{
/////		  //		  stacks[iF]->Add(result,"nostack");
/////		  //		legend->AddEntry(result,legendNames[b],"f");
/////		}
/////	      else
/////		{
/////		  stacks[iF]->Add(result);
/////		  legend->AddEntry(result,legendNames[b],"f");
/////		}
/////	    }
/////	  if(!result)
/////	    {
/////	      cout <<"null pointer returned" <<endl;
/////	    }
/////	  else
/////	    {
/////	      TCanvas c;
/////	      result->Draw();
/////	      c.Update();
/////	      c.SaveAs(outFileName);
/////	    }
/////	}
/////
/////      char stackName[200];
/////
/////      TCanvas c2;
/////
/////      sprintf(stackName,"Stack_%s_%d pions ",fileNames[iF],numPions);
/////      stacks[iF]->SetTitle(stackName);
/////      stacks[iF]->Draw();
/////      legend->Draw();
/////      c2.Update();
/////      stacks[iF]->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
/////      c2.Modified();
/////      sprintf(stackName,"Stack_%s_%d_pions.png",fileNames[iF],numPions);
/////      c2.SaveAs(stackName);
/////    }
/////
/////  TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
/////  for(int i=0;i<9;i++)
/////    {
/////      all.Add(summedHistos[i]);
/////      legend->AddEntry(summedHistos[i],allLegendNames[i],"f");
/////    }
/////  TCanvas c;
/////  sprintf(buffer,"All_%d_pions",numPions);
/////  all.SetTitle(buffer);
/////  all.Draw();
/////  legend->Draw();
/////  c.Update();
/////  all.GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
/////  c.Modified();
/////  sprintf(buffer,"All_%d_pions.png",numPions);
/////  c.SaveAs(buffer);
/////
/////  //the cases where a D** was found in the MC is treated differently from the regular D(*)pipi in data
/////  //     foundAnyDDoubleStar==true  
/////}



///available fields
//    addFieldF("leptonP");
//    addFieldF("piPlusP");
//    addFieldF("piMinusP");
//    addFieldF("dMesonP");
//    addFieldF("deltaETag");
//    addFieldF("mBTag");
//    addFieldF("logProb");
//    addFieldF("tagCorr");
//
//    addFieldI("found_2SD");
//    addFieldI("found_2SD_Star");
//    addFieldI("foundAnyDDoubleStar");
//    addFieldI("sig_numPions");
//    addFieldI("sig_numKaons");
//    addFieldI("sig_numPi0");
//    addFieldI("sig_numLeptons");
//    addFieldI("sig_numBaryons");
//
//
//    addFieldI("sig_DLNu");
//    addFieldI("sig_DPiLNu");
//    addFieldI("sig_DPiPiLNu");
//
//    addFieldI("sig_DStarLNu");
//    addFieldI("sig_DStarPiLNu");
//    addFieldI("sig_DStarPiPiLNu");
//
//    addFieldF("tagOverlapFractionCharged");
//    addFieldF("tagOverlapFractionPi0");
//
//    addFieldI("sig_dStar_2S");
//    addFieldI("sig_d_2S");
//    addFieldI("mcDecaySignature");
//    addFieldI("recDecaySignature");
//
//
//
//    addFieldI("foundDPiPi");
//    addFieldI("recBToDlNuPiPi");
//    addFieldI("recBToDlNuPi");
//    addFieldI("foundLepton");
//    addFieldI("foundPiPlus");
//    addFieldI("foundPiMinus");
//    addFieldI("foundDMeson");
//    addFieldI("DMeson_PID");
//    addFieldI("D_DaughterPID");
//    addFieldI("recDDoubleStar");
//    addFieldI("tagId");
//
//    addArrayI("recDType");
//    addArrayI("numRecPions");
//    //best mNu of all D candidates in the event
//    addArrayI("bestD");
//    //    addFieldI("numPi0");
//
//
//    //    addArrayF("tagDeltaE");
//    //    addArrayF("tagMass");
//
//    addArrayF("mNu2");
//    addArrayF("mB");
//    addArrayF("mXl");
//    addArrayF("mDnPi");
//
//    addArrayF("leptonMom");
//    addArrayF("pi1Mom");
//    addArrayF("pi2Mom");
//
//    addArrayF("leptonTheta");
//    addArrayF("pi1Theta");
//    addArrayF("pi2Theta");
//
//    addArrayF("leptonPhi");
//    addArrayF("pi1Phi");
//    addArrayF("pi2Phi");



