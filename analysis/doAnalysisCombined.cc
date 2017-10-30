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
//#define GENERATE_SINGLE_STREAM
//open box by 15%
//#define PARTIAL_BOX
#define  USE_DATA


//of course need to set doSBComb as well...
//#define onlySB


#include "TFractionFitter.h"
#include "doAnalysisCombined.h"
#include "fetchData.h"
#include "sysDataGen.h"


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
  //  bool USE_TREES=false;
  bool USE_TREES=true;
  bool doSysStudy=true;
  //initialize external vars
#ifdef PARTIAL_BOX
    //={40,20,30,20,20};
  numBins[0]=40;
  numBins[1]=20;
  numBins[2]=30;
  numBins[3]=20;
  numBins[4]=20;
#else
  //  numBins={140,140,70,140,70};
  numBins[0]=140;
  numBins[1]=140;
  numBins[2]=70;
  numBins[3]=140;
  numBins[4]=70;
#endif
  //  upperCut={2.0,2.0,0.6,2.0,0.6};
  //  lowerCut={-0.5,-0.3,-0.3,-0.3,-0.3};
  cout <<" set numb in 1 " << numBins[1] <<endl;
  upperCut[0]=2.0;
  upperCut[1]=2.0;
  upperCut[2]=0.6;
  upperCut[3]=2.0;
  cout <<" set numb in 1 " << numBins[1] <<endl;
  upperCut[4]=0.6;
  lowerCut[0]=-0.5;
  lowerCut[1]=-0.3;
  lowerCut[2]=-0.3;
  lowerCut[3]=-0.3;
  lowerCut[4]=-0.3;
  cout <<" set numb in 1 " << numBins[1] <<endl;
  for(int i=0;i<7;i++)
    {
      gl_xFeedFraction[i]=-1;
      gl_signalFraction[i]=-1;
      gl_signalInt[i]=-1;
      gl_xFeedInt[i]=-1;
      gl_dataInt[i]=-1;
    }
#ifdef DPiPi_SEARCH
  SIG_IDX= SIG_IDX_D_PI_PI;
#else
  SIG_IDX =SIG_IDX_D_PI;
#endif
  bool doSBComp=false;
	//        bool doSBComp=true;
#ifdef onlySB
	//otherwise doesn't make sense
	doSBComp=true;
#endif

  //remember that partial_box changes bins, so have to rerun in between
	int loadFromFile=false;
	//			   int loadFromFile=true;
  cout <<" set numb in 1 " << numBins[1] <<endl;
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
  cout <<" set numb in 1 2 " << numBins[1] <<endl;
  //4th tree is the data tree
  int dataTreeSize=trees[4]->GetEntries();
  maxDataTreeSize=dataTreeSize;
  totalTreeSize=dataTreeSize;
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
  char filename[300];
  sprintf(filename,"myFile_numPions%d.root",numPions);
  if(!loadFromFile)
    {
      mF=new TFile(filename,"RECREATE");
    }
  else
    {
      //READ should be default
      mF=new TFile(filename,"READ");
    }
  //more just in case, pointers are cheap
  sigSignificance=new TH1D*[20];
  vector<int> channelsUnderConsideration;
  channelsUnderConsideration.push_back(-1);
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

      sysDataGen sys(trees);
      cout <<"constructed sysData " <<endl;
#ifndef onlySB
      sys.readTrees();
#endif

      //test with BToD
      char tmpBuffer[200];
      char tmpBufferC[200];
      char tmpBufferCnvs[200];
      cout <<"let's get the templates " <<endl;
//      for(int c=0;c<5;c++)
//	{
//	  getChannelString(c-1,tmpBufferC);
//	  cout <<"channel string: "<< tmpBufferC <<endl;
//	  TH1F** ret=sys.getTemplates(c,0,tmpBufferC);
//	  cout <<"returned " << endl;
//	  cout <<" got histos, looking at " << ret[1]->GetNbinsX()<<endl;
//
//	  TCanvas* cnvs;
//	  for(int i=0;i<11;i++)
//	    {
//	      sprintf(tmpBufferCnvs,"cnvs_77_%d_%s",i,tmpBufferC);
//	      cnvs=new TCanvas(tmpBufferCnvs,tmpBufferCnvs,10,10,800,600);
//	      cout <<"forming file name " << endl;
//	      sprintf(tmpBuffer,"treeLoaded_%s_%s.png",sys.histoNames[i],tmpBufferC);
//	      cout <<"trying to draw: " << i <<endl;
//	      cout <<tmpBuffer <<endl;
//	      cout <<"prehceck : "<< ret[i]->GetNbinsX()<<endl;
//	      ret[i]->Draw();
//	      cout <<" done drawing " << endl;
//	      cnvs->SaveAs(tmpBuffer);
//	      cout << cnvs <<endl;
//	      cout <<" done " <<endl;
//	      //	      delete cnvs;
//	      cout<<" done delete" <<endl;
//	    }
//	}

      TH1F sbChi2("sb_chi2","sb_chi2",20,0,2);
      TH1F* data[7];
      TH1F* allPulls=new TH1F("data_mcPulls","data_mcPulls",60,-3,3);
      allPulls->GetXaxis()->SetTitle("pulls");

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

      //            for(int i=0;i<3;i++)
	      for(int i=0;i<1;i++)      
      //      	      for(int i=1;i<3;i++)      
	{
	  leptonId=leptonIds[i];
#ifndef onlySB
	  if(!loadFromFile)
	    {
	      //components give the components for each file separately
	      if(iC<4)
		{
		  if(pionIndex!=0 || !USE_TREES)
		    getMCComponents(trees,components,summedComponents[glChannelIdx][i][pionIndex], numPions,leptonId,iC);
		  else
		    {
		      cout <<"check templates again" <<endl;
		      TH1F** ret=sys.getTemplates(glChannelIdx,leptonId,channelBuffer,components);
		      cout <<"done" <<endl;
		      cout <<"let's clone "<<endl;
		      for(int iL=0;iL<gl_numComponents;iL++)
			{
			  summedComponents[glChannelIdx][i][pionIndex][iL]=(TH1F*)ret[iL]->Clone();
			}

		      cout <<"done cloneing .." <<endl;
		    }
		}
	      else
		{
		  cout <<"greater iC" <<endl;
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
	  cout <<"and more " << endl;
	  for(int iF=0;iF<4;iF++)
	    {
	      for(int b=0;b<11;b++)
		{
		  cout <<"trying to access component " << iF*11+b <<", iF: "<< iF <<" b: "<< b <<endl;
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
	    cout <<"calling save stack.." <<endl;
	  //have to call 'save Stack' to set e.g. 'allLegendNames'...
	

	  saveStack(components,summedComponents[glChannelIdx][i][pionIndex],numPions,leptonId,iC,isCombinedChannel);
	  cout <<"saved stack " << endl;
	  //      for(
	  //the 'other BB doesn't seem to be used...'
	  //      fitFractions(trees,summedComponents,10, numPions,leptonId,false);
	  sprintf(buffer,"pulls_numPions%d_%s",numPions,channelBuffer);
	  TH1D* pulls=new TH1D(buffer,buffer,100,-3,3);
	  sprintf(buffer,"pullsFD_numPions%d_%s",numPions,channelBuffer);
	  TH1D* pullsFeedDown=new TH1D(buffer,buffer,100,-3,3);
	  //this used to be 10 components... I don't understand why, I guess that meanst that the other BB was missing




  bool dataTree=false;
  bool addNoise=false;
#ifdef  GENERATE_SINGLE_STREAM
   dataTree=false;
   addNoise=true;
#else

#ifdef USE_DATA
	   dataTree=true;
	  	   addNoise=false;
	  //	  bool addNoise=true;
#else
#ifdef PARTIAL_BOX
	 dataTree=true;
	   addNoise=false;
	  //	  bool addNoise=true;
#else
#ifdef MC_TEST
	   dataTree=true;
	   addNoise=true;
	  if(numPions==0)
	    addNoise=false;
#else
	   dataTree=false;
	  //	  bool addNoise=true;
	   addNoise=false;
#endif
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
	  float BRRatio=0.0;
	  float relUncertBRRatio=0.0;
	  float BRRatioCombinedChannel=0.0;
	  float relUncertBRRatioCombinedChannel=0.0;
	  fitFractions(data[glChannelIdx],trees,summedComponents[glChannelIdx][i][pionIndex],numComponents, numPions,leptonId,iC,dataTree,addNoise,pulls,pullsFeedDown,BRRatio,relUncertBRRatio,BRRatioCombinedChannel, relUncertBRRatioCombinedChannel, allPulls);
	
	  //regenerate templates based on the uncertainties on the weights
	  if(doSysStudy)
	    {
	      TCanvas cSys;
	      float trueBR=0.0;
	      float trueBRFeedDown=0.0;

	      //I think we have 9 separate sources of systematic uncertainty
	      for(int sysIndex=-1;sysIndex<9;sysIndex++)
		{
		  char sysName[200];
		  char sysNameFeedDown[200];
		  getChannelString(glChannelIdx-1,channelBuffer);
		  sprintf(sysName,"sysStudy_sysIndex_%d_channel_%s",sysIndex,channelBuffer);
		  //these are the total differences in the brRatios. Should be small
		  TH1D* hSys=new TH1D(sysName,sysName,50,-0.1,0.1);
		  sprintf(sysNameFeedDown,"sysPullFeedDown_sysIndex_%d_channel_%s",sysIndex,channelBuffer);
		  TH1D* hSysFeedDown=new TH1D(sysNameFeedDown,sysNameFeedDown,50,-0.2,0.2);

		  for(int l=0;l<10;l++)
		    {
		      //only needed for iC >=4 but let's just set it...
		      int DChannelIdx=1;
		      if(glChannelIdx==6)
			{
			  DChannelIdx=3;
			}
		      int DStarChannelIdx=DChannelIdx+1;
		      if(iC<4)
			{
			  TH1F** ret=0;
			  if(l==0)//first iteration, no shift to get 'true' value
			    ret==sys.getTemplates(glChannelIdx,leptonId,channelBuffer,components);
			  else
			    ret==sys.getTemplates(glChannelIdx,leptonId,channelBuffer,components,true,sysIndex);
			  for(int iL=0;iL<gl_numComponents;iL++)
			    {
			      delete summedComponents[glChannelIdx][i][pionIndex][iL];
			      summedComponents[glChannelIdx][i][pionIndex][iL]=(TH1F*)ret[iL]->Clone();
			      delete ret[iL];
			    }
			}
		      else
			{
			  for(int iL=0;iL<gl_numComponents;iL++)
			    {
			      //need to regenerate D and D* channel
			      delete summedComponents[DChannelIdx][i][pionIndex][iL];
			      TH1F** ret=0;
			      if(0==l)
				{
				  ret=sys.getTemplates(DChannelIdx,leptonId,channelBuffer,components);
				}
			      else
				{
				  ret=sys.getTemplates(DChannelIdx,leptonId,channelBuffer,components,true,sysIndex);
				}
			      summedComponents[DChannelIdx][i][pionIndex][iL]=(TH1F*)ret[iL]->Clone();
			      delete ret[iL];
			      delete summedComponents[DStarChannelIdx][i][pionIndex][iL];
			      TH1F** ret2=sys.getTemplates(DStarChannelIdx,leptonId,channelBuffer,components);
			      summedComponents[DStarChannelIdx][i][pionIndex][iL]=(TH1F*)ret2[iL]->Clone();
			      delete ret2[iL];
			      combineChannels(summedComponents[DChannelIdx][i][pionIndex],summedComponents[DStarChannelIdx][i][pionIndex],summedComponents[glChannelIdx][i][pionIndex],gl_numComponents,DChannelIdx, numPions);
			    }
			}
		      float BRRatio=0.0;
		      float relUncertBRRatio=0.0;
		      float BRRatioCombinedChannel=0.0;
		      float relUncertBRRatioCombinedChannel=0.0;
		      fitFractions(data[glChannelIdx],trees,summedComponents[glChannelIdx][i][pionIndex],numComponents, numPions,leptonId,iC,dataTree,addNoise,pulls,pullsFeedDown,BRRatio,relUncertBRRatio,BRRatioCombinedChannel, relUncertBRRatioCombinedChannel, allPulls);
		      if(l==0)
			{
			  trueBR=BRRatio;
			  trueBRFeedDown=BRRatioCombinedChannel;
			}
		      else
			{
			  hSys->Fill(trueBR-BRRatio);
			  if(iC>=4)
			    hSysFeedDown->Fill(trueBRFeedDown-BRRatioCombinedChannel);
			}

		    }//what is l? Iterations? 
		  //do something with these numbers
		  //--->

		  hSys->Fit("gaus");
		  hSys->Draw();
		  sprintf(buffer,"%s.png",sysName);
		  cSys.SaveAs(buffer);
		  if(iC>=4)
		    {
		      hSysFeedDown->Fit("gaus");
		      sprintf(buffer,"%s,png",sysNameFeedDown);
		      hSysFeedDown->Draw();
		      cSys.SaveAs(buffer);
		    }

		}
	      ///save histos with sys pulls
	      //////
	    }///end doSysStudy
	


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

	  pullsFeedDown->Draw();
	  pullsFeedDown->Fit("gaus");
	  pullsFeedDown->Draw();
	  sprintf(buffer,"pullsFeedDown_numPions%d_%s.png",numPions,channelBuffer);
	  cpulls.SaveAs(buffer);
	  sprintf(buffer,"pullsFeedDown_numPions%d_%s.pdf",numPions,channelBuffer);
	  cpulls.SaveAs(buffer);
	  sprintf(buffer,"pullsFeedDown_numPions%d_%s.eps",numPions,channelBuffer);
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

	      char buffer2[200];
	      sprintf(buffer2,"counts/ %.3f GeV^{2}",(lowerSidebandTop-lowerSidebandBottom)/numBinsSB);
	      (*lowerSidebandData)->GetYaxis()->SetTitle(buffer2);
	      cout <<" setting y axis title to " << buffer2<<endl;
	      (*lowerSidebandData)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*lowerSidebandMC)->GetYaxis()->SetTitle(buffer2);
	      (*lowerSidebandMC)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");

	      sprintf(buffer2,"counts/ %.3f GeV^{2}",(upperSidebandTop-upperSidebandBottom)/numBinsSBTop);
	      (*upperSidebandData)->GetYaxis()->SetTitle(buffer2);
	      (*upperSidebandMC)->GetYaxis()->SetTitle(buffer2);

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
	      double scaleFactor=0.2*1.0;
	      scaleUpper=scaleFactor;
	      scaleLower=scaleFactor;
#else

#ifdef MC_TEST
	      scaleUpper=0.25;
	      scaleLower=0.25;
#else
	      scaleUpper=0.2;
	      scaleLower=0.2;
#endif
#endif

	      Double_t* zeroErrors=new Double_t[500];
	      for(int i=0;i<500;i++)
		{
		  zeroErrors[i]=0.0;
		}

	      cout <<"scaleUpper: "<< scaleUpper <<" lower: "<< scaleLower <<endl;
	   
	      (*upperSidebandMC)->Scale(scaleUpper);
	      (*lowerSidebandMC)->Scale(scaleLower);

	      cout <<" upper SB chi2: "<< getSBChi2(*upperSidebandData, *upperSidebandMC) <<endl;
	      cout <<" lower SB chi2: "<< getSBChi2(*lowerSidebandData, *lowerSidebandMC) <<endl;

	      sbChi2.Fill(getSBChi2(*upperSidebandData, *upperSidebandMC));
	      sbChi2.Fill(getSBChi2(*lowerSidebandData, *lowerSidebandMC));

	      TCanvas c("c","c",0,0,1000,800);
	      c.Divide(2,2);
	      c.cd(1);
	      //      stacks[iF]->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
	      (*lowerSidebandData)->Draw("E1A");	   
	      c.Update();
	      (*lowerSidebandMC)->Draw("hist");	   

	   
	      c.cd(2);

	      (*upperSidebandData)->Draw("E1A");
	      c.Update();
	      (*upperSidebandMC)->Draw("hist SAME");	   

	      c.cd(3);
	      TH1F* lowerSidebandMCDiv=(TH1F*)(*lowerSidebandMC)->Clone("lowerSBMCDiv");
	      lowerSidebandMCDiv->Sumw2();
	      lowerSidebandMCDiv->SetError(zeroErrors);
	      //	      (*lowerSidebandData)->SetError(zeroErrors);
	      lowerSidebandMCDiv->Divide(*lowerSidebandData);

	      lowerSidebandMCDiv->Draw();
	      lowerSidebandMCDiv->Fit("pol0");
	      c.cd(4);

	      TH1F* upperSidebandMCDiv=(TH1F*)(*upperSidebandMC)->Clone("upperSBMCDiv");
	      upperSidebandMCDiv->Sumw2();
	      upperSidebandMCDiv->SetError(zeroErrors);
	      //	      (*upperSidebandData)->SetError(zeroErrors);
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
	      (*lowerSidebandData)->Draw("E1");
	      c.Update();
	      (*lowerSidebandMC)->Draw("hist SAME");

	      sprintf(buffer,"LowerSidebandComparison_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"LowerSidebandComparison_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"LowerSidebandComparison_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      c.cd(0);
	      (*upperSidebandData)->Draw("E1");
	      for(int i=1;i<(*upperSidebandData)->GetNbinsX();i++)
		{
		  cout <<"bin " << i <<" sideband data " << (*upperSidebandData)->GetBinContent(i) <<" error: "<< (*upperSidebandData)->GetBinError(i) <<endl;
		}
	      c.Update();
	      (*upperSidebandMC)->Draw("hist SAME");
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


	      sprintf(buffer2,"counts/ %.3f GeV^{2}",(upperCut[iC+1]-lowerCut[iC+1])/numBinsWS);
	      (*sameChargeData)->GetYaxis()->SetTitle(buffer2);
	      (*sameChargeData)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*sameChargeMC)->GetYaxis()->SetTitle(buffer2);
	      (*sameChargeMC)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*chargeNeutralData)->GetYaxis()->SetTitle(buffer2);
	      (*chargeNeutralData)->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}] ");
	      (*chargeNeutralMC)->GetYaxis()->SetTitle(buffer2);
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

	      double scaleSameCharge=0.2;
	      double scaleChargeNeutral=0.2;

	      cout <<"scalesameCharge: "<< scaleSameCharge <<" chargeNeutral: "<< scaleChargeNeutral <<endl;
	   
	      (*sameChargeMC)->Scale(scaleSameCharge);
	      (*chargeNeutralMC)->Scale(scaleChargeNeutral);
	      (*sameChargeData)->Draw("E1");
	      cout <<" same charge chi2: "<< getSBChi2(*sameChargeData, *sameChargeMC) <<endl;
	      cout <<" charge neutral chi2: "<< getSBChi2(*chargeNeutralData, *chargeNeutralMC) <<endl;
	      sbChi2.Fill(getSBChi2(*sameChargeData, *sameChargeMC));
	      sbChi2.Fill(getSBChi2(*chargeNeutralData, *chargeNeutralMC));
	      c.Update();
	      (*sameChargeMC)->Draw("hist SAME");
	      (*chargeNeutralMC)->SetLineColor(kBlue);
	      (*chargeNeutralData)->SetLineColor(kBlack);

	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      (*chargeNeutralData)->Draw("E1");	   

	      c.Update();
	      (*chargeNeutralMC)->Draw("hist SAME");
	      sprintf(buffer,"chargeNeutral_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"chargeNeutral_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"chargeNeutral_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);


	      c.Divide(2,2);
	      c.cd(1);
	      (*sameChargeData)->Draw("E1");

	      c.Update();
	      (*sameChargeMC)->Draw("hist SAME");	   

	      c.cd(2);
	   
	      (*chargeNeutralData)->Draw("E1");
	      c.Update();
	      (*chargeNeutralMC)->Draw("hist SAME");	   

	      c.cd(3);
	      TH1F* sameChargeMCDiv=(TH1F*)(*sameChargeMC)->Clone("sameChargeMCDiv");
	      sameChargeMCDiv->Sumw2();
	      sameChargeMCDiv->SetError(zeroErrors);
	      //	      (*sameChargeData)->SetError(zeroErrors);
	      sameChargeMCDiv->Divide(*sameChargeData);

	      sameChargeMCDiv->Draw();
	      sameChargeMCDiv->Fit("pol0");
	      c.cd(4);

	      TH1F* chargeNeutralMCDiv=(TH1F*)(*chargeNeutralMC)->Clone("chargeNeutralMCDiv");
	      chargeNeutralMCDiv->Sumw2();
	      chargeNeutralMCDiv->SetError(zeroErrors);
	      //	      (*chargeNeutralData)->SetError(zeroErrors);
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
	      (*sameChargeData)->Draw("E1");
	      (*sameChargeMC)->SetLineColor(kBlue);
	      (*sameChargeData)->SetLineColor(kBlack);
	      c.Update();
	      (*sameChargeMC)->Draw("hist SAME");

	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"SameCharge_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	   
	      c.cd(0);
	      (*chargeNeutralData)->Draw("E1");

	      (*chargeNeutralMC)->SetLineColor(kBlue);
	      (*chargeNeutralData)->SetLineColor(kBlack);
	      c.Update();
	      (*chargeNeutralMC)->Draw("hist SAME");
	      sprintf(buffer,"ChargeNeutral_NumPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"ChargeNeutral_NumPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	      sprintf(buffer,"ChargeNeutral_NumPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelBuffer);
	      c.SaveAs(buffer);
	    }//end if(doSB..)

	}
    }




  TCanvas c77;
  allPulls->Fit("gaus");
  allPulls->Draw();
  c77.SaveAs("allMCDataPulls.png");
  c77.SaveAs("allMCDataPulls.pdf");
  sbChi2.Draw();
  cout <<"mean chi2: " << sbChi2.GetMean() <<endl;
  sbChi2.GetXaxis()->SetTitle("reduced #chi^{2}");
  c77.SaveAs("sbChi2.pdf");
  c77.SaveAs("sbChi2.png");
  c77.SaveAs("sbChi2.root");

}

void copyHisto(TH1* first, TH1* out, int numBins, int maxBins,int shift)
{

  for(int iBin=1;iBin<=shift;iBin++)
    {
      out->SetBinContent(iBin,0);
      out->SetBinError(iBin,0);
    }
  for(int iBin=1;iBin<=numBins;iBin++)
    {
      int index=iBin+shift;
      out->SetBinContent(index,first->GetBinContent(iBin));
      out->SetBinError(index,first->GetBinError(iBin));
    }
  for(int iBin=numBins+shift+1;iBin<=maxBins;iBin++)
    {
      out->SetBinContent(iBin,0);
      out->SetBinError(iBin,0);
    }


}
void addShiftedHistos(TH1* first, TH1* second, TH1* out, int numBins1, int numBins2)
{
  for(int iBin=1;iBin<=numBins1;iBin++)
    {
      out->SetBinContent(iBin,first->GetBinContent(iBin));
      out->SetBinError(iBin,first->GetBinError(iBin));
    }
  for(int iBin=numBins1+1;iBin<=(numBins1+numBins2);iBin++)
    {
      //jaja.... should be zero, could be ommitted
      out->SetBinContent(iBin,second->GetBinContent(iBin-numBins1));
      out->SetBinError(iBin,second->GetBinError(iBin-numBins1));
    }
}


void getDataFromMC(TH1F* &data, int channel, int numPions, int leptonId, TH1F** summedComponents, int numComponents)
{
  char channelString[1000];
  char histoName[2009];
  getChannelString(channel,channelString);
  sprintf(histoName,"histo_Data_%d_pions_%d_leptonId_treeNr%d_%s",numPions,leptonId,0, channelString);
  data=(TH1F*) summedComponents[0]->Clone(histoName);
   for(int i=1;i<numComponents;i++)
     {
       data->Add(summedComponents[i]);
     }
}

void getTemplates(TH1F** summedComponents_in, TH1F** &templates, char** templateLegendNames, char** allLegendNames, int numComponents, int& numMergers, int numPions,int oneIdx,bool combineDPiPi)
{
  char buffer[3000];
  //should be only 10 (or num components, this is just for safety)

  for(int i=0;i<numComponents;i++)
    {
      templateLegendNames[i]=new char[300];
    }


  //set to the pionId 1 to merge components...


  ///make new templates that take into account that some contributions look the same, so probably create problems while fitting
  //same shape combination would be 5 with 7 (only one merger) --> alternatively used 4 mergers 5, 6, 7, 8,9 where all the others have small statistics
  numMergers=0;
  if(numPions==oneIdx)
    //    numMergers=1;
    //from the code it looks like two mergers
    numMergers=2;
  if(combineDPiPi)
    {
      numMergers=2;
    }
    //    numMergers=4;

  templates=new TH1F*[numComponents-numMergers];
  //  fillTemplates(templates,summedComponents,templateLegendNames,numPions);
  cout <<"filling " << numComponents-numMergers << " templates " <<endl;

  if(numPions==oneIdx)
    {
      for(int i=0;i<5;i++)
	{
	  templates[i]=summedComponents_in[i];
	  templateLegendNames[i]=allLegendNames[i];
	  cout <<"filled template " << i <<endl;
	}
      sprintf(buffer,"%s_clonedAgain",summedComponents_in[5]->GetName());
      templates[5]=(TH1F*)summedComponents_in[5]->Clone(buffer);
      templates[5]->Add(summedComponents_in[6]);
      templates[5]->Add(summedComponents_in[7]);
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
	  templates[i-numMergers]=summedComponents_in[i];
	  templateLegendNames[i-numMergers]=allLegendNames[i];
	  cout <<"filling template " << i-numMergers <<endl;
	}
    }
  else{
    if(combineDPiPi)
      {
       for(int i=0;i<4;i++)
	{
	  templates[i]=summedComponents_in[i];
	  templateLegendNames[i]=allLegendNames[i];
	  cout <<"filled template " << i <<endl;
	}
      sprintf(buffer,"%s_clonedAgain",summedComponents_in[4]->GetName());
      templates[4]=(TH1F*)summedComponents_in[4]->Clone(buffer);
      templates[4]->Add(summedComponents_in[6]);
      templates[4]->Add(summedComponents_in[8]);
      //           templates[5]->Add(summedComponents[8]);
      //      templates[5]->Add(summedComponents[9]);

      //      sprintf(buffer,"%s plus %s and %s and  %s and %s",allLegendNames[5],allLegendNames[6],allLegendNames[7],allLegendNames[8], allLegendNames[9]);
      sprintf(buffer,"%s plus %s and %s",allLegendNames[4],allLegendNames[6],allLegendNames[8]);
      templates[4]->SetTitle(buffer);
      sprintf(templateLegendNames[4],"%s",buffer);
      //            templates[6]=summedComponents[6];
      //      templateLegendNames[6]=allLegendNames[6];
      //      for(int i=10;i<numComponents;i++)
      templates[5]=summedComponents_in[5];
      templateLegendNames[5]=allLegendNames[5];
      templates[6]=summedComponents_in[7];
      templateLegendNames[6]=allLegendNames[7];
      //already got [8] (added to 4), so next up is 7 <-- 9
      for(int i=9;i<numComponents;i++)
	{
	  templates[i-numMergers]=summedComponents_in[i];
	  //names are the same for the DStar stuff
	  templateLegendNames[i-numMergers]=allLegendNames[i];
	  cout <<"filling template " << i-numMergers <<endl;
	}

      }
  else
    {
      for(int i=0;i<numComponents;i++)
	{
	  templates[i]=summedComponents_in[i];
	  templateLegendNames[i]=allLegendNames[i];
	}
    }
  }
  //////
}





/*
get the histogram we want to fit from the tree (given numPions and leptonId), fit with the fractions from 'summedComponents'

 */
void fitFractions(TH1F* data, TTree** trees, TH1F** summedComponents, int numComponents,int numPions,int leptonId, int channel, bool dataTree, bool addNoise, TH1D* pulls,TH1D* pullsFeedDown, float& BRRatio, float& relUncertaintyBRRatio ,float& BRRatioCombChannel, float& relUncertaintyBRRatioCombChannel, TH1F* allPulls)
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
      int fixThresholdCounts=300;
      //    int fixThresholdCounts=3000;
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
      double crossFeedFraction=0;
      //  gl_templateScaleFactor=templates[iDDStarPiCrossFeed-2]->Integral()/data->Integral();
      
      cout <<"for BR, signal for channel " << channel<< "  in MC is : " << templates[SIG_IDX]->Integral() <<endl;;
      cout <<"we have " << numComponents <<" components, " << numMergers <<" mergers, so overall: "<< numComponents-numMergers << " templates"<<endl;
      cout <<" we think that the cross feed has index " << iDDStarPiCrossFeed-2 << " and that its integral is " <<      templates[iDDStarPiCrossFeed-2]->Integral() <<endl;
      
      double mcSignalIntegral=templates[SIG_IDX]->Integral();
      ///---
      
      double mcSignalIntegralCrossFeed=0;
      if(channel>3)
	{
	  //get templates merges templates, so the index has to be reduced by the number of mergers (should be 2)
	  mcSignalIntegralCrossFeed=templates[iDDStarPiCrossFeed-2]->Integral();
	  cout <<"mc signal int cross feed: "<< mcSignalIntegralCrossFeed <<endl;
	}
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
#ifdef DO_ROO_FIT
      double S=getFitSignal_RooFit(data,templates,numComponents-numMergers,result,mcPredictions,fitVal, fitErr, fixThresholdCounts,allFitVals, allFitErrs, numEffective,_effectiveComponentsIndices,_status, numPions);  
#else
      double S=getFitSignal(data,templates,numComponents-numMergers,result,mcPredictions,fitVal, fitErr, fixThresholdCounts,allFitVals, allFitErrs, numEffective,_effectiveComponentsIndices,_status, numPions);  
#endif
      float fitSum=0.0;
      float tempIntegral=0.0;
      for(int i=0;i<numEffective;i++)
	{
	  fitSum+=allFitVals[i];
	  tempIntegral+=templates[i]->Integral();
	}
      
      cout <<"sum of all fractions: " << fitSum <<", templates: " << tempIntegral <<endl;
      float templateSignalFraction=gl_lastSignalFraction;
      float templateCrossFeedFraction=gl_lastCrossFeedFraction;
      
      //  float signalFraction=templates[SIG_IDX]/data->Integral();
      
      cout <<"got " << fitVal*data->Integral() << " signal counts, fraction: " << fitVal << " +- " << fitErr  << endl;
      /////-----------------
      //void performFractionFitStabilityTest(TH1F** templatesOrg, TH1F* dataOrg, int numComponents)
      if(addNoise)
	{
	  //scale data by 0.2 since otherwise it is the sum of the 5 streams...
	  //	data->Scale(0.2);
	  //-->do this in the 'performFractionFit... function'
	  int maxIterations=4000;
	  if(glChannelIdx<5)
	    {
	      maxIterations=1;
	    }
	  
	  for(int nIt=0;nIt<maxIterations;nIt++)
	    {
	      float locSignalFraction=signalFraction;
	      float locCrossFeedFraction=crossFeedFraction;
	      if(glChannelIdx<5 && gl_signalFraction[glChannelIdx]>0)
		{
		  locSignalFraction=gl_signalFraction[glChannelIdx];
		  locCrossFeedFraction=gl_xFeedFraction[glChannelIdx];
		}
	      if(glChannelIdx==5 && gl_signalInt[1]>0 && gl_dataInt[1]>0 && gl_dataInt[2]>0)
		{
		  locSignalFraction=gl_signalInt[1]/(gl_dataInt[1]+gl_dataInt[2]);
		  //D* channel signal is crossfeed for the combined case
		  locCrossFeedFraction=(gl_xFeedInt[1]+gl_signalInt[2])/(gl_dataInt[1]+gl_dataInt[2]);
		}
	      if(glChannelIdx==5 && gl_signalInt[3]>0 && gl_dataInt[3]>0 && gl_dataInt[4]>0)
		{
		  locSignalFraction=gl_signalInt[3]/(gl_dataInt[3]+gl_dataInt[4]);
		  locCrossFeedFraction=(gl_xFeedInt[3]+gl_signalInt[4])/(gl_dataInt[3]+gl_dataInt[4]);
		}
	      cout <<"locSignalFraction: "<< locSignalFraction <<" template signal fraction: " << templateSignalFraction<<endl;
	      //let's try this:
#ifdef GENERATE_SINGLE_STREAM
	      locSignalFraction=templateSignalFraction;
	      locCrossFeedFraction=templateCrossFeedFraction;
#endif
	      performFractionFitStabilityTest(templates,data,numComponents-numMergers,pulls,pullsFeedDown, locSignalFraction, locCrossFeedFraction,fixThresholdCounts, numPions);
	    }
	}
      /////////-----------done with the fraction stability test (pulls etc)

      TCanvas sampleData;
      data->Draw();
      sampleData.SaveAs("sampleData.png");
      sampleData.SaveAs("sampleData.pdf");
      sampleData.SaveAs("sampleData.eps");
      
      int numEffectiveComponents=0;
      cout << " status: "<< _status <<endl;
      if (_status == 0) {                       // check on fit status
	TCanvas c;
	cout <<"grabbing result .." <<endl;
	double templatePredIntegral=result->Integral();
	cout <<" done " <<endl;
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
	sprintf(buffer,"SummedfracFitComp_numPions_%d_leptonId_%d_%s",numPions,leptonId,channelString);
	TH1F* summedPredComps=new TH1F(buffer,buffer,numBins[glChannelIdx],lowerCut[glChannelIdx],upperCut[glChannelIdx]);


	//    cout <<"fraction fitter has " << _fit->GetFitter()->GetNumberFreeParameters() << " free and " << _fit->GetFitter()->GetNumberTotalParameters() <<" overall parameters" <<endl;
	//we want to flip the signal (index =2 ) so that it is later
	int signalIdx=-1;
	int signalIdxCrossFeed=-1;
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
		summedPredComps->Add(mcComp);
		cout <<"add as : " << templateLegendNames[_effectiveComponentsIndices[i]]<<endl;
		legend->AddEntry(mcComp,templateLegendNames[_effectiveComponentsIndices[i]],"f" );
	      }
	    else
	      {
		signalIdx=i;
	      }
	    //again, some of the templates are merged, so there is a shift by 2 in the indices
	    if(_effectiveComponentsIndices[i]==(iDDStarPiCrossFeed-2))
	      {
		cout << " setting cross feed index to : "<< i <<endl;
		signalIdxCrossFeed=i;
	      }
	  }
	//add the signal last...
	if(signalIdx>=0)
	  {
	    cout <<"dealing with the signal " << endl;
	    //	TH1F* mcComp=(TH1F*) _fit->GetMCPrediction(signalIdx);
	    TH1F* mcComp=mcPredictions[signalIdx];
	    TH1F* mcCompCrossFeed;
	    //again, some of the templates are merged, so there is a shift by 2 in the indices
	    if(channel>3)
	      {
		mcCompCrossFeed=mcPredictions[signalIdxCrossFeed];
		cout <<" x-feed predition is: "<< mcCompCrossFeed->Integral() <<" original template: " << templates[iDDStarPiCrossFeed-numMergers] <<" (mergers: " << numMergers <<")" <<endl;
		cout <<"temp Int : " << tempIntegral <<" so template fraction is " << templates[iDDStarPiCrossFeed-numMergers]->Integral()/tempIntegral << ", ratio to fit: " <<  allFitVals[signalIdxCrossFeed]/(templates[iDDStarPiCrossFeed-numMergers]->Integral()/tempIntegral) <<endl; 
	    
	    
	      }
	    cout <<"got pred " << endl;
	    double mcPredInt=mcComp->Integral();
	    double mcPredIntCrossFeed=0;
	    if(channel>3)
	      {
		mcPredIntCrossFeed=mcCompCrossFeed->Integral();
		cout <<"mc predIntCross feed: "<< mcPredIntCrossFeed <<endl;
	      }
	    sumOfCompInts+=mcPredInt;
	    float scaleFact=integralResult/mcPredInt*allFitVals[signalIdx];
	    float scaleFactCrossFeed=0;
	    if(channel>3)
	      scaleFactCrossFeed=integralResult/mcPredIntCrossFeed*allFitVals[signalIdxCrossFeed];
	    cout <<"got scalefact " << endl;
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
	    /////---->	    
	    cout <<"BR ratio to MC for channel " << channelString << " is : "<< mcComp->Integral()/mcSignalIntegral<<endl;
	    cout <<"fit val: "<< allFitVals[signalIdx] << " uncert: "<< allFitErrs[signalIdx] << "  relative uncert  " << allFitErrs[signalIdx]/allFitVals[signalIdx] <<endl; 
	    if(channel>3)
	      {
		mcComp->Scale(scaleFactCrossFeed);
		cout <<"using index: "<< signalIdxCrossFeed <<endl;
		cout <<"BR ratio to DDStar MC for channel " << channelString << " is : "<< mcCompCrossFeed->Integral()/mcSignalIntegralCrossFeed<<endl;
		cout <<"fit val: "<< allFitVals[signalIdxCrossFeed] << " uncert: "<< allFitErrs[signalIdxCrossFeed] << "  relative uncert  " << allFitErrs[signalIdxCrossFeed]/allFitVals[signalIdxCrossFeed] <<endl; 
	      }

	    BRRatio=mcComp->Integral()/mcSignalIntegral;
	    relUncertaintyBRRatio=allFitErrs[signalIdx]/allFitVals[signalIdx];
	    if(channel>3)
	      {
		BRRatioCombChannel=mcCompCrossFeed->Integral()/mcSignalIntegralCrossFeed;
		relUncertaintyBRRatioCombChannel=allFitErrs[signalIdxCrossFeed]/allFitVals[signalIdxCrossFeed];
	      }
	    
	    predComponents->Add(mcComp);
	    summedPredComps->Add(mcComp);
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

	char buffer2[200];
	sprintf(buffer2,"counts/ %.2f GeV^{2}",(upperCut[glChannelIdx]-lowerCut[glChannelIdx])/(float)numBins[glChannelIdx]);
	data->GetYaxis()->SetTitle(buffer2);
	data->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV^{2}]");
	data->SetTitle("");
	data->Draw("Ep");
	predComponents->SetTitle("");
	//    predComponents->SetStats(0);
    //    predComponents->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
	predComponents->Draw("hist same");
	data->SetLineWidth(2);
	data->SetStats(0);
	data->Draw("same Ep");
	TH1F* pulls=getCompPulls(summedPredComps,data);
	allPulls->Add(pulls);

	legend->Draw();
	sprintf(buffer,"predComp_numPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelString);
	c.SaveAs(buffer);
	sprintf(buffer,"predComp_numPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelString);
	c.SaveAs(buffer);
	sprintf(buffer,"predComp_numPions_%d_leptonId_%d_%s.eps",numPions,leptonId,channelString);
	c.SaveAs(buffer);
	pulls->Fit("gaus");
	pulls->Draw();
	sprintf(buffer,"pullsPredComp_numPions_%d_leptonId_%d_%s.png",numPions,leptonId,channelString);
	c.SaveAs(buffer);
	sprintf(buffer,"pullsPredComp_numPions_%d_leptonId_%d_%s.pdf",numPions,leptonId,channelString);
	c.SaveAs(buffer);

      }
}
  ////

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



