///Implementation of the class
#include "tables/mctype.h"
#include "bToDDoubleStar/PIDCorrections.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "bToDDoubleStar/HelperFunctions.h"
#include "bToDDoubleStar/LeptonFakes.h"
#include "TMath.h"
#include "event/BelleEvent.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "tuple/BelleTupleManager.h"

 #include "belle.h"
 #include "panther/panther.h"
#include "panther/panther_config.h"
#include MCTYPE_H
#include MDST_H
#include EVTCLS_H

//#include "bToDDoubleStar/bToDDoubleStar.h"

using namespace std;


float PIDCorrections::getXSectionCorrection(bool isMC, int mcType)
{
  if(!isMC)
    return 1.0;

  //  mixed: 1001
  //    charged: 1002
  //    charm: 1003
  //    uds: 1004

 float mcXSectionBB=1.09;
 float mcXSectionUDS=2.09;
  float mcXSectionCharm=1.3;

 float rdXSectionBB=1.1;
  float rdXSectionUDSC=3.3;

  float rdY4S_chargedBXSection=0.514;
  float rdY4S_neutralBXSection=0.486;

  float eRdY4S_chargedBXSection=0.006;
  float eRdY4S_neutralBXSection=0.006;


  float mcX=1.0;
  float dataX=1.0;
  //mixed
  if(mcType==1001)
    {
       mcX=0.5*mcXSectionBB;
       dataX=rdXSectionBB*rdY4S_neutralBXSection;
    }
  //charged
  if(mcType==1002)
    {
      mcX=0.5*mcXSectionBB;
      dataX=rdXSectionBB*rdY4S_chargedBXSection;
    }
  //charm
  if(mcType==1003)
    {
      //mc is udsc x-section
      mcX=mcXSectionUDS+mcXSectionCharm;
      dataX=rdXSectionUDSC;
    }
  //uds
  if(mcType==1004)
    {
      mcX=mcXSectionUDS+mcXSectionCharm;
      dataX=rdXSectionUDSC;
    }
  return dataX/mcX;
}

float PIDCorrections::getLumiCorrection(int exp)
{
  return rdLumi[exp]/mcLumi[exp];

};
float PIDCorrections::getWeight(int mcLund, int dataLund, float mom, float theta,int expNr, int runNr)
{
  int absLund=fabs(dataLund);
  int absMC=fabs(mcLund);
  if((absLund==211 || absLund==321) && (absMC==211 || absMC==321))
    {
      return getWeightChargedHadron(mcLund,dataLund,mom,theta,expNr);
    }
  if((absLund==11 || absLund==13) && (absMC==211 || absMC==321))
    {
      return getWeightLepton(mcLund,dataLund,mom,theta,expNr,runNr);
    }
  //pi0
  if(absLund==111)
    {
      return getWeightPi0(mom);
    }
  //ks
  if(absLund==310)
    return getWeightKs(mom,theta);

  return 1.0;
}
float PIDCorrections::getWeightLepton(int mcLund, int dataLund, float mom, float theta,int expNr, int runNr)
{
  int absLundData=fabs(dataLund);
  int absLundMC=fabs(mcLund);

  int momBin=getBin(limitsLeptonP,mom);
  int thetaBin=-1;


  if(absLundData==11)
    {

      thetaBin=getBin(limitsThetaE,cos(theta));
    }
  else
    {
      thetaBin=getBin(limitsThetaMu,cos(theta));
    }

  if(absLundData==11)
    {
      if(absLundMC==321)
	{
	  return k2e[thetaBin][momBin];
	}
      else
	{
	  return pi2e[thetaBin][momBin];
	}
      if(absLundMC==11)
	{
	  //theta bin is the same as for fakes, but mom bin is different
	  momBin=getBin(limitsLeptonPEff,mom);
	  ///		  leptonEff[1][iTbl-2][momBin][thetaBin]=r;
	  //
	  int tableI=0;
	  if(expNr>=31)
	    tableI=1;
	  return leptonEff[0][tableI][momBin][thetaBin];
	}
    }
  if(absLundData==13)
    {
      if(absLundMC==321)
	{
	  return k2mu[thetaBin][momBin];
	}
      else
	{
	  return pi2mu[thetaBin][momBin];
	}
      if(absLundData==13)
	{
	  ///efficiency momentum bins are different
	  momBin=getBin(limitsLeptonPEff,mom);
	  int tableI=0;
	  if((expNr>=31 && expNr <= 39) || (expNr==45 && runNr<=220))
	     tableI=1;
	  if(expNr>=41 && expNr<=49)
	    tableI=2;
	  if(expNr>49)
	    tableI=3;
	  return leptonEff[1][tableI][momBin][thetaBin];
	}
    }

}

float PIDCorrections::getWeightPi0(float mom)
{
  int momBin=getBin(limitsPi0P,mom);
  return pi0Eff[momBin];
}

float PIDCorrections::getWeightKs(float mom, float theta)
{

  float cosTheta=cos(theta);
  //r=100.06 %, error 3.75 %
  if(mom < 0.5)
    return 1.0006;
  if(mom>0.5 && mom< 1.5)
    {
      if(cosTheta<-0.5)
	{
	  //err 3.52%
	  return 0.9808;
	}
      if(cosTheta>-0.5 && cosTheta<0)
	{
	  //err 0.87%
	  return 0.9716;
	}
      if(cosTheta>0 && cosTheta <0.5)
	{
	  //err 0.75%
	  return 0.9729;
	}
      if(cosTheta>0.5)
	{
	  //err 1.04%
	  return 0.9863;
	}
    }
  if(mom > 1.5)
    {
      if(cosTheta<-0.5)
	{
	  //err 3.84%
	  return 0.9423;
	}
      if(cosTheta>-0.5 && cosTheta<0)
	{
	  //err 2.3%
	  return 0.9657;
	}
      if(cosTheta>0 && cosTheta <0.5)
	{
	  //err 1.18%
	  return 0.9984;
	}
      if(cosTheta>0.5)
	{
	  //err 1.22%
	  return 0.9938;
	}
    }
  return 1.0;

}
float PIDCorrections::getWeightChargedHadron(int mcLund, int dataLund, float mom, float theta,int expNr)
{
  int momBin=getBin(pLabBoundary,mom);
  int cosThetaBin=getBin(cosThetaBoundary,cos(theta));
  int svdBin=0;
  if(expNr>=31)
    svdBin=1;
  //pion
  int iPid=0;
  //eff
  int iKind=0;
  //looking at kaon
  if(fabs(dataLund)==321)
    {
      iPid=1;
    }
      //misId, e.g. saw pion but is kaon, so we should read out the pid K fake rate
  if(fabs(mcLund)!=fabs(dataLund))
    {
      iKind=1;
    }
  //  cout <<" mcLund " << mcLund <<" data lund: "<< dataLund <<endl;
  //  cout <<" looking at mom: " << mom << " theta: "<< theta <<" cos(theta): " << cos(theta) <<" mom bin: " << momBin <<" cosThetaBin: " << cosThetaBin <<" ratio: ";
  //  cout <<ratio[svdBin][iPid][iKind][cosThetaBin][momBin] <<endl;
  
 return ratio[svdBin][iPid][iKind][cosThetaBin][momBin];

}

bool PIDCorrections::loadLeptonTables()
{
  float pidCut;
  float minTheta;
  float maxTheta;
  float minMom;
  float maxMom;
  float r;
  float statErr;
  float sys1;
  float sys2;


  ifstream lidFile1("eid_data-mc_corr_svd1.dat");
  ifstream lidFile2("eid_data-mc_corr_exp31-65-caseB.dat");
  ifstream lidFile3("muid_data-mc_corr_svd1.dat");
  ifstream lidFile4("muid_data-mc_corr_exp31-39-45a-caseB.dat");
  ifstream lidFile5("muid_data-mc_corr_exp41-49-caseB.dat");
  ifstream lidFile6("muid_data-mc_corr_exp51-65-caseB.dat");

  string str;
  ifstream* files[6];
  files[0]=&lidFile1;
  files[1]=&lidFile2;
  files[2]=&lidFile3;
  files[3]=&lidFile4;
  files[4]=&lidFile5;
  files[5]=&lidFile6;
  //e cut: 0.8, mu_cut 0.9
  for(int iTbl=0;iTbl<6;iTbl++)
    {
      while(getline(*(files[iTbl]),str))
	{
	  if(str.size()>0 && str[0]!='#')
	    {
	      stringstream ss(str);
	      ss >> pidCut >> minTheta>> maxTheta>> minMom >> maxMom >> r >> statErr >> sys1>> sys2;
	      //muon table
	      bool isMu=false;
	      int momBin=getBin(limitsLeptonPEff,(maxMom-minMom)/2+minMom);
	      if(iTbl>1)
		{
		  isMu=true;
		  if(pidCut!=0.9)
		    continue;
		  //mu max theta bins are 25, 37,51,117,130,145,150
		  float maxRad=maxTheta*TMath::Pi()/180;
		  int thetaBin=-1;

		  for(int i=0;i<limitsThetaMu.size();i++)
		    {
		      if(fabs(cos(maxRad)-limitsThetaMu[i])<0.01)
			thetaBin=i;
		    }
		  if(thetaBin<0)
		    {
		      cout <<"theta bin negative! mu " <<endl;
		      exit(1);
		    }
		  leptonEff[1][iTbl-2][momBin][thetaBin]=r;
		  leptonEffStat[1][iTbl-2][momBin][thetaBin]=statErr;
		  leptonEffSys1[1][iTbl-2][momBin][thetaBin]=sys1;
		  leptonEffSys2[1][iTbl-2][momBin][thetaBin]=sys2;
		}
	      else
		{
		  //e max theta bins are 25,35,40,60,125,132,151
		  if(pidCut!=0.8)
		    continue;
		  float maxRad=maxTheta*TMath::Pi()/180;
		  int thetaBin=-1;

		  for(int i=0;i<limitsThetaE.size();i++)
		    {
		      if(fabs(cos(maxRad)-limitsThetaE[i])<0.01)
			thetaBin=i;
		    }
		  if(thetaBin<0)
		    {
		      cout <<"theta bin negative! electron " <<endl;
		      exit(1);
		    }
		  //first 0 is for electron
		  leptonEff[0][iTbl][momBin][thetaBin]=r;
		  leptonEffStat[0][iTbl][momBin][thetaBin]=statErr;
		  leptonEffSys1[0][iTbl][momBin][thetaBin]=sys1;
		  leptonEffSys2[0][iTbl][momBin][thetaBin]=sys2;
		}
	    }
	}
    }
  return true;
}

bool PIDCorrections::loadKIDTables()
{


  int kind;
  int pid;
  int map;
  float effData;
  float effDataUncert;
  float effDataSys;
  float effMC;
  float effMCUncert;
  float mRatio;
  float mRatioUncert;
  float mRatioSys;
  int flag;


  string str;


  ifstream kidFile1("kidTable_SVD1.data");
  ifstream kidFile2("kidTable_SVD2.data");

  ifstream* files[2];
  files[0]=&kidFile1;
  files[1]=&kidFile2;

  for(int iSvd=0;iSvd<2;iSvd++)
    {
      while(getline(*(files[iSvd]),str))
	{
	  if(str.size()>0 && str[0]!='#')
	    {
	      stringstream ss(str);
	      ss >> kind >> pid >> map >> effData >> effDataUncert >> effDataSys >> effMC >> effMCUncert >> mRatio >> mRatioUncert >> mRatioSys >> flag;
	      //	  cout <<"read kind: "<< kind << " pid: "<< pid << " effData: "<< effData <<endl;
	      
	      int pLabBin=(map /1000)*10+((map / 100) %10);
	      int thetaBin=((map / 10)%10) *10+(map % 10);
	      //we use indices starting from 0, not 1;
	      pLabBin--;
	      thetaBin--;
	      //	  cout <<"map: " << map << " plabBin: " << plabBin <<" thetaBin: " << thetaBin <<endl;


	      int iPid=0;
	      if(kind==0 || 1==kind)
		{
		  //kaon eff or fake
		  iPid=1;
		}

	      //not the correct kinematic bin
	      if(iPid==0 && pid!=pidBin)
		continue;
	      if(iPid==1 && pid!=kidBin)
		continue;

	      int iKind=0;
	      if(1==kind || 3==kind)
		{
		  //fake rate, not efficiency
		  iKind=1;
		}
	      cout <<"read: " << str <<endl;
	      cout <<" iPid: "<< iPid << " iKind : "<< iKind <<" thetaBin: " << thetaBin <<" plabBin: "<< pLabBin << " ratio: "<< mRatio <<" uncert: "<< mRatioUncert <<" sys: "<< mRatioSys <<endl;
	      cout <<"iSVD: " << iSvd <<endl;
	      //fit not good or no good stats
	      if(flag)
		{
		  mRatio=1.0;
		  mRatioUncert=0.0;
		  mRatioSys=0.0;
		}

	      ratio[iSvd][iPid][iKind][thetaBin][pLabBin]=mRatio;
	      ratioStat[iSvd][iPid][iKind][thetaBin][pLabBin]=mRatioUncert;
	      ratioSys[iSvd][iPid][iKind][thetaBin][pLabBin]=mRatioSys;
	    }
	}
    }

}
