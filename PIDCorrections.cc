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


pair<float,float> PIDCorrections::getXSectionCorrection(bool isMC, int mcType)
{
  if(!isMC)
    return pair<float,float>(1.0,0.0);

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
  float eRatio=0.0;
  //mixed
  if(mcType==1001)
    {
       mcX=0.5*mcXSectionBB;
       dataX=rdXSectionBB*rdY4S_neutralBXSection;
       eRatio=eRdY4S_neutralBXSection*rdXSectionBB;
    }
  //charged
  if(mcType==1002)
    {
      mcX=0.5*mcXSectionBB;
      dataX=rdXSectionBB*rdY4S_chargedBXSection;
      eRatio=eRdY4S_chargedBXSection*rdXSectionBB;
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
  return pair<float,float>(dataX/mcX,eRatio/mcX);
}

pair<float,float> PIDCorrections::getLumiCorrection(int exp)
{
  //apparently no significant systematic uncertainty
  return pair<float,float>(rdLumi[exp]/mcLumi[exp],0.0);

};

//return value and uncertainty
pair<float,float> PIDCorrections::getWeight(int mcLund, int dataLund, float mom, float theta,int expNr, int runNr,pidRet& mRet)
{

  int absLund=fabs(dataLund);
  int absMC=fabs(mcLund);
  if((absLund==211 || absLund==321) && (absMC==211 || absMC==321))
    {
      return getWeightChargedHadron(mcLund,dataLund,mom,theta,expNr, mRet);
    }

  //found lepton which is pion kaon or efficiency weight of correctly identified lepton
  if((absLund==11 || absLund==13) && (absMC==211 || absMC==321 || (absLund==11 && absMC==11) || (absLund==13 && absMC==13)))
    {
      return getWeightLepton(mcLund,dataLund,mom,theta,expNr,runNr, mRet);
    }
  //pi0
  if(absLund==111)
    {
      return getWeightPi0(mom,mRet);
    }
  //ks
  if(absLund==310)
    return getWeightKs(mom,theta,mRet);

  return pair<float,float>(1.0,0.0);
}
pair<float,float> PIDCorrections::getWeightLepton(int mcLund, int dataLund, float mom, float theta,int expNr, int runNr, pidRet& mPidRet)
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


  mPidRet.thetaBin=thetaBin;
  mPidRet.momBin=momBin;
  mPidRet.svdBin=0;
  if(absLundData==11)
    {
      if(absLundMC==321)
	{
	  mPidRet.misIdType=4;
	  return pair<float,float>(k2e[thetaBin][momBin],k2e_err[thetaBin][momBin]);
	}
      else
	{
	  mPidRet.misIdType=5;
	  return pair<float,float>(pi2e[thetaBin][momBin],pi2e_err[thetaBin][momBin]);
	}
      if(absLundMC==11)
	{
	  mPidRet.misIdType=6;
	  //theta bin is the same as for fakes, but mom bin is different
	  momBin=getBin(limitsLeptonPEff,mom);
	  ///		  leptonEff[1][iTbl-2][momBin][thetaBin]=r;
	  //
	  int tableI=0;
	  if(expNr>=31)
	    tableI=1;

	  mPidRet.svdBin=tableI;
	  float err=leptonEffStat[0][tableI][momBin][thetaBin]*leptonEffStat[0][tableI][momBin][thetaBin];
	  err+=leptonEffSys1[0][tableI][momBin][thetaBin]*leptonEffSys1[0][tableI][momBin][thetaBin];
	  err+=leptonEffSys2[0][tableI][momBin][thetaBin]*leptonEffSys2[0][tableI][momBin][thetaBin];

	  return pair<float,float>(leptonEff[0][tableI][momBin][thetaBin],sqrt(err));
	}
    }
  if(absLundData==13)
    {
      if(absLundMC==321)
	{
	  mPidRet.misIdType=7;
	  return pair<float,float>(k2mu[thetaBin][momBin],k2mu_err[thetaBin][momBin]);
	}
      else
	{
	  mPidRet.misIdType=8;
	  return pair<float,float>(pi2mu[thetaBin][momBin],pi2mu_err[thetaBin][momBin]);
	}
      if(absLundData==13)
	{
	  mPidRet.misIdType=9;
	  ///efficiency momentum bins are different
	  momBin=getBin(limitsLeptonPEff,mom);
	  int tableI=0;
	  if((expNr>=31 && expNr <= 39) || (expNr==45 && runNr<=220))
	     tableI=1;
	  if(expNr>=41 && expNr<=49)
	    tableI=2;
	  if(expNr>49)
	    tableI=3;

	  mPidRet.svdBin=tableI;

	  float err=leptonEffStat[1][tableI][momBin][thetaBin]*leptonEffStat[1][tableI][momBin][thetaBin];
	  err+=leptonEffSys1[1][tableI][momBin][thetaBin]*leptonEffSys1[1][tableI][momBin][thetaBin];
	  err+=leptonEffSys2[1][tableI][momBin][thetaBin]*leptonEffSys2[1][tableI][momBin][thetaBin];

	  return pair<float,float>(leptonEff[1][tableI][momBin][thetaBin],sqrt(err));
	}
    }

}

pair<float,float> PIDCorrections::getWeightPi0(float mom, pidRet& mPidRet)
{
  int momBin=getBin(limitsPi0P,mom);
  mPidRet.momBin=momBin;
  return pair<float,float>(pi0Eff[momBin],pi0EffStat[momBin]);
}

pair<float,float> PIDCorrections::getWeightKs(float mom, float theta, pidRet& mPidRet)
{

  float cosTheta=cos(theta);

  float weight=1.0;
  float relErr=0.0;
  //r=100.06 %, error 3.75 %, fine to just return 3.75% because the weight is pretty much equal to one
  mPidRet.momBin=0;
  mPidRet.thetaBin=0;
  if(mom < 0.5)
    {
      mPidRet.momBin=0;
      mPidRet.thetaBin=0;
      return pair<float,float>(1.0006,0.0375);
    }
  if(mom>0.5 && mom< 1.5)
    {
      mPidRet.momBin=1;
      if(cosTheta<-0.5)
	{
	  mPidRet.thetaBin=0;
	  relErr=0.0352;
	  weight=0.9808;
	  return pair<float,float>(weight,relErr*weight);
	}
      if(cosTheta>-0.5 && cosTheta<0)
	{
	  mPidRet.thetaBin=1;
	  relErr=0.0087;
	  weight= 0.9716;
	  return pair<float,float>(weight,relErr*weight);
	}
      if(cosTheta>0 && cosTheta <0.5)
	{
	  mPidRet.thetaBin=2;
	  relErr= 0.0075;
	  weight= 0.9729;
	  return pair<float,float>(weight,relErr*weight);
	}
      if(cosTheta>0.5)
	{
	  mPidRet.thetaBin=3;
	  relErr=0.0104;
	  weight=0.9863;
	  return pair<float,float>(weight,relErr*weight);
	}
    }
  if(mom > 1.5)
    {
      mPidRet.momBin=2;

      if(cosTheta<-0.5)
	{
	  mPidRet.thetaBin=0;
	  relErr=0.0384;
	  weight= 0.9423;
	  return pair<float,float>(weight,relErr*weight);
	}
      if(cosTheta>-0.5 && cosTheta<0)
	{
	  mPidRet.thetaBin=1;
	  relErr=0.023;
	  weight=0.9657;
	  return pair<float,float>(weight,relErr*weight);
	}
      if(cosTheta>0 && cosTheta <0.5)
	{
	  mPidRet.thetaBin=2;
	  relErr= 0.0118;
	  weight= 0.9984;
	  return pair<float,float>(weight,relErr*weight);
	}
      if(cosTheta>0.5)
	{
	  mPidRet.thetaBin=3;
	  relErr=0.0122;
	  weight=0.9938;
	  return pair<float,float>(weight,relErr*weight);
	}
    }
  return pair<float,float>(weight,relErr*weight);
  

}
pair<float,float> PIDCorrections::getWeightChargedHadron(int mcLund, int dataLund, float mom, float theta,int expNr, pidRet& mPidRet)
{
  int momBin=getBin(pLabBoundary,mom);
  int cosThetaBin=getBin(cosThetaBoundary,cos(theta));
  int svdBin=0;
  mPidRet.momBin=momBin;
  mPidRet.thetaBin=cosThetaBin;
  mPidRet.svdBin=svdBin;
  mPidRet.misIdType=0;

  if(expNr>=31)
    svdBin=1;

  mPidRet.svdBin=svdBin;
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
  
  mPidRet.misIdType=iPid*2+iKind;
  float err=ratioStat[svdBin][iPid][iKind][cosThetaBin][momBin]*ratioStat[svdBin][iPid][iKind][cosThetaBin][momBin]+ratioSys[svdBin][iPid][iKind][cosThetaBin][momBin]*ratioSys[svdBin][iPid][iKind][cosThetaBin][momBin];;
  err=sqrt(err);
  return pair<float,float>(ratio[svdBin][iPid][iKind][cosThetaBin][momBin],err);

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


  ifstream lidFile1("/home/belle/vossen/myProjects/bToDDoubleStar-/eid_data-mc_corr_svd1.dat");
  ifstream lidFile2("/home/belle/vossen/myProjects/bToDDoubleStar-/eid_data-mc_corr_exp31-65-caseB.dat");
  ifstream lidFile3("/home/belle/vossen/myProjects/bToDDoubleStar-/muid_data-mc_corr_svd1.dat");
  ifstream lidFile4("/home/belle/vossen/myProjects/bToDDoubleStar-/muid_data-mc_corr_exp31-39-45a-caseB.dat");
  ifstream lidFile5("/home/belle/vossen/myProjects/bToDDoubleStar-/muid_data-mc_corr_exp41-49-caseB.dat");
  ifstream lidFile6("/home/belle/vossen/myProjects/bToDDoubleStar-/muid_data-mc_corr_exp51-65-caseB.dat");

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


  ifstream kidFile1("/home/belle/vossen/myProjects/bToDDoubleStar-/kidTable_SVD1.data");
  ifstream kidFile2("/home/belle/vossen/myProjects/bToDDoubleStar-/kidTable_SVD2.data");

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
	      //	      cout <<"read: " << str <<endl;
	      //	      cout <<" iPid: "<< iPid << " iKind : "<< iKind <<" thetaBin: " << thetaBin <<" plabBin: "<< pLabBin << " ratio: "<< mRatio <<" uncert: "<< mRatioUncert <<" sys: "<< mRatioSys <<endl;
	      //	      cout <<"iSVD: " << iSvd <<endl;
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
