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
void doAnalysis(char* fileNameMixed, char* fileNameCharged, char* fileNameUds, char* fileNameCharm, int numPions=2);
int main(int argc, char** argv)
{

  if(argc!=6)
    {
      cout <<"see only " << argc << " arguments, need 6 though.." <<endl;
      exit(0);
    }
  char* fileNameMixed=argv[1];
  char* fileNameCharged=argv[2];
  char* fileNameUds=argv[3];
  char* fileNameCharm=argv[4];

  int numPions=atoi(argv[5]);
  cout <<"want " << numPions <<endl;
  doAnalysis(fileNameMixed,fileNameCharged,fileNameUds,fileNameCharm, numPions);

}

void doAnalysis(char* fileNameMixed, char* fileNameCharged, char* fileNameUds, char* fileNameCharm, int numPions)
{

  /*
we have 4 data files, the two B MCs, charm and uds. For each we extract DDouble star +0,1,2 Pions, as well as the Dlnu +0,1,2 pions, 
and the rest

We do this for the 0,1,2 pions in the final state

The selection criteria are 
i) The correct final state has to be reconstructed
ii) The correct 

*/
  int numBins=100;
  char* fileNames[4];
  char* legendNames[11];

  char* allLegendNames[11];

  fileNames[0]=new char[200];
  fileNames[1]=new char[200];
  fileNames[2]=new char[200];
  fileNames[3]=new char[200];

  THStack all;
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

  TH1F* summedHistos[11];
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



  for(int i=0;i<11;i++)
    {
      legendNames[i]=new char [200];
      allLegendNames[i]=new char[200];
    }

    sprintf(legendNames[0],"D Double Star (no Dn#pi l#nu)");
    sprintf(legendNames[1],"D Double Star #pi (no Dn#pi l#nu)");
    sprintf(legendNames[2],"D Double Star #pi #pi (no Dn#pi l#nu)");
    sprintf(legendNames[3],"D l #nu");
    sprintf(legendNames[4],"D #pi l #nu");
    sprintf(legendNames[5],"D #pi #pi l #nu");

    sprintf(legendNames[6],"D* l #nu");
    sprintf(legendNames[7],"D* #pi l #nu");
    sprintf(legendNames[8],"D* #pi #pi l #nu");
    sprintf(legendNames[9]," no D(*)n #pi l#nu ");
    sprintf(legendNames[10]," no Selection");




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

  TFile* files[4];
  TTree* trees[4];

  files[0]=new TFile(fileNameMixed);
  files[1]=new TFile(fileNameCharged);
  files[2]=new TFile(fileNameUds);
  files[3]=new TFile(fileNameCharm);


  TH1F** histos[4];
  for(int i=0;i<4;i++)
    {
      histos[i]=new TH1F*[7];
    }


  for(int i=0;i<4;i++)
    {
       trees[i] = (TTree*)files[i]->Get("DataTree");
       if(!trees[i])
	 cout <<"tree " << i << " is NULL" <<endl;
    }

  //  dataTree->Draw(" >> result")
  //  TH1F* result=(TH1F*) gDirectory->Get("result");

    //selection for the data

    //numRecPions  
    //  recDecaySignature==true   //recDecay signature found

    // get mNu2 for d double star decays


    //this is for all the same, the selection of the data. It can be 0, one or two pions

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


    sprintf(buffer,"tagCorr*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0",numPions);

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


    char histoName[2009];
    char drawCommand[2000];
    char outFileName[2000];

    THStack* stacks[4];
    






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


    for(int iF=0;iF<4;iF++)
      {
	char bufferF[200];
	sprintf(bufferF,"stackF_%d",iF);
	stacks[iF]=new THStack(bufferF,bufferF);
	TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
	int allCounts=0;
	for(int b=0;b<11;b++)
	  {
	    sprintf(histoName,"histo_If_%d_b_%d",iF,b);
	    sprintf(drawCommand,"mNu2 >> %s(%d,-1,1)",histoName,numBins);
	    sprintf(outFileName,"%s.png",histoName);
	    cout <<"draw command: " << drawCommand << ", selections: " << (char*) selections[b]<<endl;
	    cout <<"drawing tree " << iF <<endl;
	    int counts=trees[iF]->Draw(drawCommand,(char*)selections[b]);
	    if(b<10)
	      allCounts+=counts;
	    if(b==10)
	      cout <<"all counts so far: "<< allCounts <<" no selection counts: "<< counts <<" difference: "<< counts-allCounts <<endl;
	    cout <<"got " << counts <<" counts selected " <<endl;
	    TH1F* result=(TH1F*)gDirectory->Get(histoName);
	    result->SetFillColor(colorTable[b]->GetNumber());

	    if(b!=10)
	      {
		if(iF>=2)
		  {
		    hContinuum->Add(result);
		    hContinuum->SetFillColor(colorTable[10]->GetNumber());
		  }
		
		if(iF<2)
		  {
		    summedHistos[b+1]->Add(result);
		    //for the color, make it consistent to the other stacks
		    summedHistos[b]->SetFillColor(colorTable[b]->GetNumber());
		  }
	      }

	    if(counts>0)
	      {
		if(b==10)
		  {
		    //		  stacks[iF]->Add(result,"nostack");
		    //		legend->AddEntry(result,legendNames[b],"f");
		  }
		else
		  {
		    stacks[iF]->Add(result);
		    legend->AddEntry(result,legendNames[b],"f");
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

        char stackName[200];

	TCanvas c2;

	sprintf(stackName,"Stack_%s_%d pions ",fileNames[iF],numPions);
	stacks[iF]->SetTitle(stackName);
	stacks[iF]->Draw();
	legend->Draw();
	c2.Update();
	stacks[iF]->GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
	c2.Modified();
	sprintf(stackName,"Stack_%s_%d_pions.png",fileNames[iF],numPions);
	c2.SaveAs(stackName);
      }

    TLegend* legend =new TLegend(0.6,0.7,0.89,0.99);
    for(int i=0;i<11;i++)
      {
	all.Add(summedHistos[i]);
	legend->AddEntry(summedHistos[i],allLegendNames[i],"f");
      }
    TCanvas c;
    sprintf(buffer,"All_%d_pions",numPions);
    all.SetTitle(buffer);
    all.Draw();
    legend->Draw();
    c.Update();
    all.GetXaxis()->SetTitle("m_{#nu}^{2} [GeV]");
    c.Modified();
    sprintf(buffer,"All_%d_pions.png",numPions);
    c.SaveAs(buffer);

    //the cases where a D** was found in the MC is treated differently from the regular D(*)pipi in data
    //     foundAnyDDoubleStar==true  
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
