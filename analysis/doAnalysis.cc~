#include "TCut.h"
#include "TTree.h"


void doAnalysis(char* fileNameMixed, char* fileNameCharged, char* fileNameUds, char* fileNameCharm, int numPions=2)
{

  /*
we have 4 data files, the two B MCs, charm and uds. For each we extract DDouble star +0,1,2 Pions, as well as the Dlnu +0,1,2 pions, 
and the rest

We do this for the 0,1,2 pions in the final state

The selection criteria are 
i) The correct final state has to be reconstructed
ii) The correct 

*/


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
    }

  //  dataTree->Draw(" >> result")
  //  TH1F* result=(TH1F*) gDirectory->Get("result");

    //selection for the data

    //numRecPions  
    //  recDecaySignature==true   //recDecay signature found

    // get mNu2 for d double star decays


    //this is for all the same, the selection of the data. It can be 0, one or two pions

  char** selections[7];

    char buffer[2000];

    char bufferDDStar[2000];
    char bufferDDStarPi[2000];
    char bufferDDStarPiPi[2000];

    char bufferDlNu[2000];
    char bufferDPilNu[2000];
    char bufferDPiPilNu[2000];
    char bufferAll[2000];

    selections[0]=bufferDDStar;
    selections[1]=bufferDDStarPi;
    selections[2]=bufferDDStarPiPi;
    selections[3]=bufferDlNu;
    selections[4]=bufferDPilNu;
    selections[5]=bufferDPiPilNu;
    selections[6]=bufferAll;


    sprintf(buffer,"tagCorr*(mNu2<1.0 && mNu2>1.0 && numRecPions==%d",numPions);
    sprintf(bufferDDStar,"%s && foundAnyDDoubleStar==1 && sig_numPions==0 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1)",buffer);
    sprintf(bufferDDStarPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==1 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1)",buffer);
    sprintf(bufferDDStarPiPi,"%s && foundAnyDDoubleStar==1 && sig_numPions==2 && sig_numKaons==0 && sig_numPi0==0 && sig_numBaryons==0&& sig_numLeptons==1)",buffer);

    sprintf(bufferDlNu,"%s && !foundAnyDDoubleStar && sigDLNu )",buffer);
    sprintf(bufferDPilNu,"%s && !foundAnyDDoubleStar && sigDPiLNu )",buffer);
    sprintf(bufferDPiPilNu,"%s && !foundAnyDDoubleStar && sigDPiPiLNu )",buffer);
    sprintf(bufferAll,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2) && !sigDLNu && !sigDPiLNu && !sigDPiPiLNu)",buffer);


    cout <<"selection for DDstar: "<< bufferDDStar <<endl<<" DDstarPi: "<< bufferDDStarPi<<endl<< " DDstarPiPi: "<< bufferDDStarPiPi <<endl;
    cout <<" DlNu: " << bufferDlNu<<endl <<" DPilNu:"<< bufferDPilNu << endl << " DPiPilNu: "<< bufferDPiPilNu <<endl;
    cout <<" all: "<< bufferAll <<endl;


    char histoName[2009];
    char drawCommand[2000];
    char outFileName[2000];
    for(int iF=0;iF<4;iF++)
      {
	for(int b=0;b<7;b++)
	  {
	    sprintf(histoName,"histo_If_%d_b_%d",iF,b);
	    sprintf(drawCommand,"mNu2 >> %s",histoName);
	    sprintf(outFileName,"%s.png",histoName);
	    cout <<"draw command: " << drawCommand << ", selections: " << (char*) selections[b]<<endl;
	    trees[iF]->Draw(drawCommand,(char*)selections[b]);
	    TH1F* result=(TH1F*)gDirectory->Get(histoName);
	    if(!results)
	      {
		cout <<"null pointer returned" <<endl;
	      }
	    else
	      {
		result->SaveAs(outFileName);
	      }
	  }
      }





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
