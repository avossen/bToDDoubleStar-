#ifndef TREESAVER_H
#define TREESAVER_H
#include "bToDDoubleStar/mc.h"  //one central place to put the define mc
#include <mdst/mdst.h>
#include MDST_H
#include EVTCLS_H
#include MDST_OBS_H
#include HEPEVT_H
#include TRK_H

#include "BTreeData.h"
#include "belle.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "tuple/BelleTupleManager.h"
#include "bToDDoubleStar/EventInfo.h"
#include "bToDDoubleStar/AnaConsts.h"
#include "bToDDoubleStar/AuxFunc.h"
#include "bToDDoubleStar/mc.h"
#include "bToDDoubleStar/ParticleInfoMass.h"
#include <math.h>
//#include "AnaDefs.h"
#include HEPEVT_H
#include BELLETDF_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
/*class that facilitates the 
saving to root trees
*/
// class HadronQuadruple;

//#define NUM_F_FIELDS 18 //without the data that is only saved for no mc
//#define NUM_F_FIELDS 22 //without the data that is only saved for no mc
//#define NUM_F_FIELDS 22 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
//#define NUM_F_FIELDS 25 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
////-->with the jet stuff, subtract 8
#define NUM_F_FIELDS 15 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
//6-->10 with the genIds..
#define NUM_I_FIELDS 6



class TreeSaver
{
public:
  TreeSaver()
  {
    initialize();
  };



  void doBranching()
  {
    addFieldF("leptonP");
    addFieldF("piPlusP");
    addFieldF("piMinusP");
    addFieldF("dMesonP");
    addFieldF("deltaETag");
    addFieldF("mBTag");
    addFieldF("logProb");
    addFieldF("tagCorr");
    addFieldF("D_DecayCorr");
    addFieldF("B_DecayCorr");
    addFieldF("PIDCorrection");
    addFieldF("CrossSectionLumiCorrection");  



    addFieldF("pi1Mom_mc");
    addFieldF("pi2Mom_mc");
    addFieldF("pi1Theta_mc");
    addFieldF("pi2Theta_mc");
    addFieldF("pi1Phi_mc");
    addFieldF("pi2Phi_mc");

    addFieldI("pi1Found");
    addFieldI("pi2Found");
    addFieldI("overlapEvent");
    addFieldI("found_2SD");
    addFieldI("found_2SD_Star");
    addFieldI("foundAnyDDoubleStar");
    addFieldI("bestBCharge");
    addFieldI("mcBCharge");



    addFieldI("sig_numPions");
    addFieldI("sig_numKaons");
    addFieldI("sig_numPi0");
    addFieldI("sig_numLeptons");
    addFieldI("sig_numBaryons");
    addFieldI("sig_numD");
    addFieldI("sig_numDStar");



    addFieldI("sig_DLNu");
    addFieldI("sig_DPiLNu");
    addFieldI("sig_DPiPiLNu");

    addFieldI("sig_DStarLNu");
    addFieldI("sig_DStarPiLNu");
    addFieldI("sig_DStarPiPiLNu");


    addFieldI("sig_ResDLNu");
    addFieldI("sig_ResDPiLNu");
    addFieldI("sig_ResDPiPiLNu");

    addFieldI("sig_ResDStarLNu");
    addFieldI("sig_ResDStarPiLNu");
    addFieldI("sig_ResDStarPiPiLNu");



    addFieldF("tagOverlapFractionCharged");
    addFieldF("tagOverlapFractionPi0");

    addFieldI("sig_dStar_2S");
    addFieldI("sig_d_2S");
    addFieldI("mcDecaySignature");
    addFieldI("recDecaySignature");


    addFieldI("mcDCharge");
    addFieldI("mcIsDStar");
    addFieldI("foundDPiPi");
    addFieldI("recBToDlNuPiPi");
    addFieldI("recBToDlNuPi");
    addFieldI("foundLepton");
    addFieldI("leptonId");
    addFieldI("foundPiPlus");
    addFieldI("foundPiMinus");
    addFieldI("foundDMeson");
    addFieldI("DMeson_PID");
    addFieldI("D_DaughterPID");
    addFieldI("recDDoubleStar");
    addFieldI("tagId");


    addArrayI("leptonCharge");
    addArrayI("dCharge");
    addArrayI("systemCharge");
    addArrayI("recDType");
    addArrayI("numRecPions");

    //best mNu of all D candidates in the event
    addArrayI("bestD");
    //    addFieldI("numPi0");


    //    addArrayF("tagDeltaE");
    //    addArrayF("tagMass");
    addArrayI("dType");
    addArrayI("dDecay");

    addArrayF("mNu2");
    addArrayF("U");
    addArrayF("mB");
    addArrayF("mXl");
    addArrayF("mDnPi");

    addArrayF("leptonMom");
    addArrayF("pi1Mom");
    addArrayF("pi2Mom");


    addArrayF("leptonTheta");
    addArrayF("pi1Theta");
    addArrayF("pi2Theta");


    addArrayF("leptonPhi");
    addArrayF("pi1Phi");
    addArrayF("pi2Phi");




    addArrayF("DDiff");
    addArrayF("DStarDiff");

    addArrayF("hypDMass1");
    addArrayF("hypDMass2");

  }


    static bool sgn(float f)
      {
	return f>=0;
      }


    void saveData(BTreeData* data)
    {
      int index=0;
      (*(float*)treeData[index++])=data->leptonP;
      (*(float*)treeData[index++])=data->piPlusP;
      (*(float*)treeData[index++])=data->piMinusP;
      (*(float*)treeData[index++])=data->dMesonP;
      (*(float*)treeData[index++])=data->deltaETag;
      (*(float*)treeData[index++])=data->mBTag;
      (*(float*)treeData[index++])=data->logProb;
      (*(float*)treeData[index++])=data->tagCorr;
      (*(float*)treeData[index++])=data->D_DecayCorr;
      (*(float*)treeData[index++])=data->B_DecayCorr;
      //      cout <<"saving d decay corr : "<< data->D_DecayCorr <<" b decay corr: "<< data->B_DecayCorr<<endl;
      (*(float*)treeData[index++])=data->pidCorrection;
      (*(float*)treeData[index++])=data->CrossSectionLumiCorrection;


      (*(float*)treeData[index++])=data->pi1Mom_mc;
      (*(float*)treeData[index++])=data->pi2Mom_mc;
      (*(float*)treeData[index++])=data->pi1Theta_mc;
      (*(float*)treeData[index++])=data->pi2Theta_mc;
      (*(float*)treeData[index++])=data->pi1Phi_mc;
      (*(float*)treeData[index++])=data->pi2Phi_mc;


      //      cout <<"saving corrections: tag: " << data->tagCorr <<", D: "<< data->D_DecayCorr <<" B: " << data->B_DecayCorr<<", pid: "<< data->pidCorrection <<", cross section and lumi: "<< data->CrossSectionLumiCorrection<<endl;
      (*(int*)treeData[index++])=data->pi1Found;
      (*(int*)treeData[index++])=data->pi2Found;
      (*(int*)treeData[index++])=data->overlapEvent;


      (*(int*)treeData[index++])=data->found_2SD;
      (*(int*)treeData[index++])=data->found_2SD_Star;
      (*(int*)treeData[index++])=data->foundAnyDDoubleStar;
      (*(int*)treeData[index++])=data->bestBCharge;
      (*(int*)treeData[index++])=data->mcBCharge;
      //      cout <<"any ddouble star: "<<data->foundAnyDDoubleStar <<endl;


      (*(int*)treeData[index++])=data->sig_numPions;
      (*(int*)treeData[index++])=data->sig_numKaons;
      (*(int*)treeData[index++])=data->sig_numPi0;
      (*(int*)treeData[index++])=data->sig_numLeptons;
      (*(int*)treeData[index++])=data->sig_numBaryons;
      (*(int*)treeData[index++])=data->sig_numD;
      (*(int*)treeData[index++])=data->sig_numDStar;

      (*(int*)treeData[index++])=data->sigDLNu;
      //      cout <<"dlnu : "<<data->sigDLNu <<endl;
      //      cout <<"dPilnu : "<<data->sigDPiLNu <<endl;
      //      cout <<"dPiPilnu : "<<data->sigDPiPiLNu <<endl;

      //      cout <<"dStarlnu : "<<data->sigDStarLNu <<endl;
      //      cout <<"dStarPilnu : "<<data->sigDStarPiLNu <<endl;
      //      cout <<"dStarPiPilnu : "<<data->sigDStarPiPiLNu <<endl;
      (*(int*)treeData[index++])=data->sigDPiLNu;
      (*(int*)treeData[index++])=data->sigDPiPiLNu;

      (*(int*)treeData[index++])=data->sigDStarLNu;
      (*(int*)treeData[index++])=data->sigDStarPiLNu;
      (*(int*)treeData[index++])=data->sigDStarPiPiLNu;


      (*(int*)treeData[index++])=data->sigResDLNu;
      (*(int*)treeData[index++])=data->sigResDPiLNu;
      (*(int*)treeData[index++])=data->sigResDPiPiLNu;

      (*(int*)treeData[index++])=data->sigResDStarLNu;
      (*(int*)treeData[index++])=data->sigResDStarPiLNu;
      (*(int*)treeData[index++])=data->sigResDStarPiPiLNu;


      (*(float*)treeData[index++])=data->tagOverlapFractionCharged;
      (*(float*)treeData[index++])=data->tagOverlapFractionPi0;


      (*(int*)treeData[index++])=data->sig_dStar_2S;
      (*(int*)treeData[index++])=data->sig_d_2S;
      (*(int*)treeData[index++])=data->mcDecaySignature;
      (*(int*)treeData[index++])=data->recDecaySignature;
      (*(int*)treeData[index++])=data->mcDCharge;
      (*(int*)treeData[index++])=data->mcIsDStar;

      (*(int*)treeData[index++])=data->foundDPiPi;
      (*(int*)treeData[index++])=data->recBToDlNuPiPi;
      (*(int*)treeData[index++])=data->recBToDlNuPi;
      (*(int*)treeData[index++])=data->foundLepton;
      (*(int*)treeData[index++])=data->leptonId;
      (*(int*)treeData[index++])=data->foundPiPlus;
      (*(int*)treeData[index++])=data->foundPiMinus;
      (*(int*)treeData[index++])=data->foundDMeson;
      (*(int*)treeData[index++])=data->dPID;
      (*(int*)treeData[index++])=data->daughterDPID;
      (*(int*)treeData[index++])=data->recDDoubleStar;
      (*(int*)treeData[index++])=data->tagId;


      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->leptonCharge[i];
	}
      ++index;


      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->dCharge[i];
	}
      ++index;


      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->systemCharge[i];
	}
      ++index;


      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->recDType[i];
	}
      ++index;


      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->numRecPions[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->bestD[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->dType[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((int*)treeData[index])[i]=data->dDecay[i];
	}
      ++index;


      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->mNu2[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->U[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->mB[i];
	}
      ++index;
      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->mXl[i];
	}
      ++index;
      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->mDnPi[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->leptonMom[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->pi1Mom[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->pi2Mom[i];
	}
      ++index;




      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->leptonTheta[i];
	}
      ++index;
      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->pi1Theta[i];
	}
      ++index;
      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->pi2Theta[i];
	}
      ++index;



      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->leptonPhi[i];
	}
      ++index;
      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->pi1Phi[i];
	}
      ++index;
      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->pi2Phi[i];
	}
      ++index;

      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->DDiff[i];
	}
      ++index;
      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->DStarDiff[i];
	}
      ++index;


      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->hypDMass1[i];
	}
      ++index;



      (*(int*)treeData[index++])=data->size;
      for(int i=0;i<data->size;i++)
	{
	  ((float*)treeData[index])[i]=data->hypDMass2[i];
	}
      ++index;

      pDataTree->Fill();

    }






  //get all h-pairs in mc, w/o det acceptance




 //std: float datatype
  void addFieldF(char* fieldname)
  {
    //construct the memory location from which the tree should read the new data field
    float* memLoc=new float;
    treeData.push_back(memLoc);
    pDataTree->Branch(fieldname, memLoc, (fieldname+string("/F")).c_str());
    fieldNamesF.push_back(fieldname);
  };


  void addFieldI(char* fieldname)
  {
    int* memLoc=new int;
    treeData.push_back(memLoc);
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("/I")).c_str());
    fieldNamesI.push_back(fieldname);
  };

  void addArrayI(char* fieldname)
  {
    //standard lenth, shouldn't be more than that
    int* memLoc=new int[1200];
    int* memLocCounter=new int;
    treeData.push_back(memLocCounter);
    treeData.push_back(memLoc);
    string counterName=string(fieldname)+string("Counter");
    pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/I")).c_str());
  };

  void addArrayF(char* fieldname)
  {
    float* memLoc=new float[1200];
    int* memLocCounter=new int;
    treeData.push_back(memLocCounter);
    treeData.push_back(memLoc);
    string counterName=string(fieldname)+string("Counter");
    pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());
  };


  //my guess as to what is going on here:
  //  dataF, dataI is the event only data (these arrays where deleted before handed to the event-fill function), offset is where the single fields start
  void saveData()
  {
    int dataSize=dataF.size()+dataI.size();
#ifdef MC
    if(dataSize !=treeData.size())
#else
    if(dataSize !=treeData.size())
#endif
      {

	cout << "data size does not match number of branches " <<dataSize  << " = " <<treeData.size() <<endl <<flush;
	//10+108=122
	exit(0);
      }
    //	cout << "data size does  match number of branches " <<dataSize <<"+ " << offset << " = " <<treeData.size() <<endl <<flush;
    //10+112
    //this works because dataF was cleared before event specific things there saved, 
    for(int i=0;i<dataF.size();i++)
      {
	*(float*)treeData[i]=dataF[i];
      }
    for(int i=0;i<dataI.size();i++)
      {

	*(int*)treeData[i+dataF.size()]=dataI[i];
      }
    pDataTree->Fill();

  };




private:

  template<class T> void fillTreeEntry(int& counter,T entry)
  {
    counter++;
    *(T*)treeData[counter]=entry;
  }
  template<class T> void fillTreeArrEntry(int& counter,vector<T>& vec)
  {
    counter++;
    *(int*)treeData[counter]=vec.size();
    counter++;
    //cp data from vec in tree array
    for(int i=0;i<vec.size();i++)
      {
	((T*)treeData[counter])[i]=vec[i];
      }

  }

  void initialize()
  {
    if(initialized)
      return;
    //the first time the class is initialized, construct the tree
    pDataTree=new TTree("DataTree","My Transversity Data Tree");
    initialized=true;

  };


  static TTree* pDataTree;
  static bool initialized;
  static vector<float> dataF;
  static vector<int> dataI;


  //the adresses from which the tree should read its data
  static vector<void*> treeData;
  static vector<string> fieldNamesI;
  static vector<string> fieldNamesF;

};
#if defined(BELLE_NAMESPACE)
}
#endif
#endif
