
#ifndef PIDCorrections__H
#define PIDCorrections__H
#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include "tables/mctype.h"

using namespace std;

/*
class to read in the PID tables and encapuslate the PID corrections
 */
class PIDCorrections
{

 public:


  struct  pidRet{
    int thetaBin;
    int momBin;
    //this will also encode the run, exp ranges for lepton/hadron tables
    int svdBin;
    //encodes hadron/lepton type as well... 0-3 for KP, the rest for the leptons
    int misIdType;
  };


  PIDCorrections(float kidCut=0.6, float piCut=0.4){
    mPidCut=piCut;
    mKidCut=kidCut;

    kidBin=round(mKidCut/0.1);
    pidBin=round((1.0-mPidCut)/0.1);
    cout << " pidBin : "<< pidBin <<endl;
    cout << " kidBin : "<< kidBin <<endl;
    cout <<" pidCut: "<< mPidCut <<endl;
    cout <<"1-pidCut: "<< 1.0-mPidCut <<" div by 0.1: "<< (1.0-mPidCut)/0.1<<endl;


    cosThetaBoundary.push_back(-0.612);
    cosThetaBoundary.push_back(-0.511);
    cosThetaBoundary.push_back(-0.3);
    cosThetaBoundary.push_back(-0.152);
    cosThetaBoundary.push_back(0.017);
    cosThetaBoundary.push_back(0.209);
    cosThetaBoundary.push_back(0.355);
    cosThetaBoundary.push_back(0.435);
    cosThetaBoundary.push_back(0.542);
    cosThetaBoundary.push_back(0.692);
    cosThetaBoundary.push_back(0.842);
    cosThetaBoundary.push_back(100);

    for(int i=0;(i*0.1+0.5)<=3.0;i++)
      {
	pLabBoundary.push_back(i*0.1+0.5);
      }
    pLabBoundary.push_back(3.2);
    pLabBoundary.push_back(3.4);
    pLabBoundary.push_back(3.6);
    pLabBoundary.push_back(4.0);
    pLabBoundary.push_back(4.5);
    pLabBoundary.push_back(100.00);

    //the first entry was the lower boundary of the first bin
    //    limitsThetaE.push_back(-0.87);
    limitsThetaE.push_back(-0.67);
    limitsThetaE.push_back(-0.57);
    limitsThetaE.push_back(0.5);
    limitsThetaE.push_back(0.77);
    limitsThetaE.push_back(0.82);
    limitsThetaE.push_back(0.91);
    limitsThetaE.push_back(0.95);
    limitsThetaE.push_back(1.0);

    //    limitsThetaMu.push_back(-0.87);
    limitsThetaMu.push_back(-0.82);
    limitsThetaMu.push_back(-0.64);
    limitsThetaMu.push_back(-0.45);
    limitsThetaMu.push_back(0.63);
    limitsThetaMu.push_back(0.8);
    limitsThetaMu.push_back(0.91);
    limitsThetaMu.push_back(0.96);
    limitsThetaMu.push_back(1.0);

    for(int i=1;i*0.5<=5.0;i++)
      {
	limitsLeptonPEff.push_back(i*0.5);
      }
 
    limitsLeptonP.push_back(0.5);
    limitsLeptonP.push_back(0.75);
    limitsLeptonP.push_back(1.0);
    limitsLeptonP.push_back(1.25);
    limitsLeptonP.push_back(1.5);
    limitsLeptonP.push_back(1.75);
    limitsLeptonP.push_back(2.0);
    limitsLeptonP.push_back(2.25);
    limitsLeptonP.push_back(2.5);
    limitsLeptonP.push_back(3.0);
    limitsLeptonP.push_back(4.0);
    limitsLeptonP.push_back(100);


    //svd 1/2
    ratio=new float****[2];
    ratioSys=new float****[2];
    ratioStat=new float****[2];
    for(int iSvd=0;iSvd<2;iSvd++)
      {
	//pion/kaon ? 
	ratio[iSvd]=new float***[2];
	ratioSys[iSvd]=new float***[2];
	ratioStat[iSvd]=new float***[2];
	for(int iPid=0;iPid<2;iPid++)
	  {
	    //eff or fake?
	    ratio[iSvd][iPid]=new float**[2];
	    ratioSys[iSvd][iPid]=new float**[2];
	    ratioStat[iSvd][iPid]=new float**[2];
	    for(int iKind=0;iKind<2;iKind++)
	      {
		ratio[iSvd][iPid][iKind]=new float*[cosThetaBoundary.size()];
		ratioSys[iSvd][iPid][iKind]=new float*[cosThetaBoundary.size()];
		ratioStat[iSvd][iPid][iKind]=new float*[cosThetaBoundary.size()];
		for(int i=0;i<cosThetaBoundary.size();i++)
		  {
		    ratio[iSvd][iPid][iKind][i]=new float[pLabBoundary.size()];
		    ratioStat[iSvd][iPid][iKind][i]=new float[pLabBoundary.size()];
		    ratioSys[iSvd][iPid][iKind][i]=new float[pLabBoundary.size()];
		  }

	      }

	  }
     }

    leptonFakes=new float***[2];
    leptonFakesErr=new float***[2];
    for(int iLepton=0;iLepton<2;iLepton++)
      {
	leptonFakes[iLepton]=new float**[2];
	leptonFakesErr[iLepton]=new float**[2];
	for(int iHadron=0;iHadron<2;iHadron++)
	  {
	    leptonFakes[iLepton][iHadron]=new float*[limitsLeptonP.size()];
	    leptonFakesErr[iLepton][iHadron]=new float*[limitsLeptonP.size()];
	    for(int iLP=0;iLP<limitsLeptonP.size();iLP++)
	      {
		//electron
		if(iLepton==0)
		  {
		    leptonFakes[iLepton][iHadron][iLP]=new float[limitsThetaE.size()];
		    leptonFakesErr[iLepton][iHadron][iLP]=new float[limitsThetaE.size()];
		  }
		else
		  {
		    leptonFakes[iLepton][iHadron][iLP]=new float[limitsThetaMu.size()];
		    leptonFakesErr[iLepton][iHadron][iLP]=new float[limitsThetaMu.size()];
		  }
	      }
	  }
      }


    leptonEff=new float***[2];
    leptonEffStat=new float***[2];
    leptonEffSys1=new float***[2];
    leptonEffSys2=new float***[2];
    //mu or e
    for(int iLepton=0;iLepton<2;iLepton++)
      {
	leptonEff[iLepton]=new float**[3];
	leptonEffStat[iLepton]=new float**[3];
	leptonEffSys1[iLepton]=new float**[3];
	leptonEffSys2[iLepton]=new float**[3];
	//for muons we have three different tables, two for electrons
	for(int iFile=0;iFile<3;iFile++)
	  {
	    leptonEff[iLepton][iFile]=new float*[limitsLeptonPEff.size()];
	    leptonEffStat[iLepton][iFile]=new float*[limitsLeptonPEff.size()];
	    leptonEffSys1[iLepton][iFile]=new float*[limitsLeptonPEff.size()];
	    leptonEffSys2[iLepton][iFile]=new float*[limitsLeptonPEff.size()];

	    for(int iLP=0;iLP<limitsLeptonPEff.size();iLP++)
	      {
		//electron
		if(iLepton==0)
		  {
		    leptonEff[iLepton][iFile][iLP]=new float[limitsThetaE.size()];
		    leptonEffStat[iLepton][iFile][iLP]=new float[limitsThetaE.size()];
		    leptonEffSys1[iLepton][iFile][iLP]=new float[limitsThetaE.size()];
		    leptonEffSys2[iLepton][iFile][iLP]=new float[limitsThetaE.size()];
		  }
		else
		  {
		    leptonEff[iLepton][iFile][iLP]=new float[limitsThetaMu.size()];
		    leptonEffStat[iLepton][iFile][iLP]=new float[limitsThetaMu.size()];
		    leptonEffSys1[iLepton][iFile][iLP]=new float[limitsThetaMu.size()];
		    leptonEffSys2[iLepton][iFile][iLP]=new float[limitsThetaMu.size()];
		  }
		
	      }

	  }
      }


    for(int i=1;i<=12;i++)
      {
	limitsPi0P.push_back(i*0.5);
      }
    limitsPi0P.push_back(1000);

    pi0Eff=new float[limitsPi0P.size()];
    pi0EffStat=new float[limitsPi0P.size()];
    //data/mc
    pi0Eff[0]=1.009;
    pi0Eff[1]=0.9226;
    pi0Eff[2]=0.9348;
    pi0Eff[3]=0.9295;
    pi0Eff[4]=0.9463;
    pi0Eff[5]=0.9455;
    pi0Eff[6]=0.9713;
    pi0Eff[7]=0.9699;
    pi0Eff[8]=0.9772;
    pi0Eff[9]=0.9620;
    pi0Eff[10]=0.9839;
    pi0Eff[11]=0.9567;
    pi0Eff[12]=1.0;

    pi0EffStat[0]=0.023;
    pi0EffStat[1]=0.021;
    pi0EffStat[2]=0.021;
    pi0EffStat[3]=0.021;
    pi0EffStat[4]=0.021;
    pi0EffStat[5]=0.022;
    pi0EffStat[6]=0.021;
    pi0EffStat[7]=0.022;
    pi0EffStat[8]=0.022;
    pi0EffStat[9]=0.022;
    pi0EffStat[10]=0.022;
    pi0EffStat[11]=0.023;
    pi0EffStat[12]=0.0;


    loadKIDTables();
    loadLeptonTables();

    //numbers for the lumi corrections

    rdLumi[ 7 ] = 5.922;
    rdLumi[ 9 ] = 4.436;
    rdLumi[ 11 ] = 8.125;
    rdLumi[ 13 ] = 10.73;
    rdLumi[ 15 ] = 12.52;
    rdLumi[ 17 ] = 11.181;
    rdLumi[ 19 ] = 24.963;
    rdLumi[ 21 ] = 4.375;
    rdLumi[ 23 ] = 6.273;
    rdLumi[ 25 ] = 26.954;
    rdLumi[ 27 ] = 25.43;
    rdLumi[ 31 ] = 17.725;
    rdLumi[ 33 ] = 17.508;
    rdLumi[ 35 ] = 16.691;
    rdLumi[ 37 ] = 60.909;
    rdLumi[ 39 ] = 41.157;
    rdLumi[ 41 ] = 58.752;
    rdLumi[ 43 ] = 56.206;
    rdLumi[ 45 ] = 12.946;
    rdLumi[ 47 ] = 37.205;
    rdLumi[ 49 ] = 27.024;
    rdLumi[ 51 ] = 39.237;
    rdLumi[ 55 ] = 72.088;
    rdLumi[ 61 ] = 34.095;
    rdLumi[ 63 ] = 32.858;
    rdLumi[ 65 ] = 37.751;


    mcLumi[7] = 5.795;
    mcLumi[9] = 4.116;
    mcLumi[11] = 7.696;
    mcLumi[13] = 10.704;
    mcLumi[15] = 12.686;
    mcLumi[17] = 9.228;
    mcLumi[19] = 24.624;
    mcLumi[21] = 4.315;
    mcLumi[23] = 6.326;
    mcLumi[25] = 25.485;
    mcLumi[27] = 25.355;
    mcLumi[31] = 17.025;
    mcLumi[33] = 17.383;
    mcLumi[35] = 16.870;
    mcLumi[37] = 60.884;
    mcLumi[39] = 42.355;
    mcLumi[41] = 57.327;
    mcLumi[43] = 54.612;
    mcLumi[45] = 12.794;
    mcLumi[47] = 36.726;
    mcLumi[49] = 26.671;
    mcLumi[51] = 38.771;
    mcLumi[55] = 71.268;
    mcLumi[61] = 34.174;
    mcLumi[63] = 32.163;
    mcLumi[65] = 37.173;

};

  pair<float,float> getWeight(int mcLund, int dataLund, float mom, float theta, int expNr, int runNr, pidRet& mPidRet);
  pair<float,float> getWeightChargedHadron(int mcLund, int dataLund, float mom, float theta, int expNr, pidRet& mPidRet);
  pair<float,float> getWeightPi0(float mom, pidRet& mPidRet);
  pair<float,float> getWeightKs(float mom, float theta, pidRet& mPidRet);
  pair<float,float> getWeightLepton(int mcLund, int dataLund, float mom, float theta, int expNr, int runNr, pidRet& mPidRet);

  pair<float,float> getLumiCorrection(int exp);
  pair<float,float> getXSectionCorrection(bool isMC, int mcType);

 protected:
    map<int,float> rdLumi;
    map<int,float> mcLumi;
  bool loadKIDTables();
  bool loadLeptonTables();
  vector<float> cosThetaBoundary;
  vector<float> pLabBoundary;
  float *****ratio;
  float *****ratioStat;
  float *****ratioSys;


  float ****leptonFakes;
  float ****leptonFakesErr;


  float ****leptonEff;
  float ****leptonEffStat;
  float ****leptonEffSys1;
  float ****leptonEffSys2;


  float *pi0Eff;
  float *pi0EffStat;

  float mPidCut;
  float mKidCut;

  int pidBin;
  int kidBin;

 private:
  vector< float > limitsThetaE;
  vector< float > limitsThetaMu;
  vector< float > limitsLeptonP;
  vector<float> limitsPi0P;



  //for some reason the efficiency limits seem to be different...
  vector< float > limitsLeptonPEff;

};

#endif
