#include "sysDataGen.h"
#include "idx2.h"

const int maxEntries=30000;
sysDataGen::sysDataGen(TTree** trees)
{
  pContinuum=new dataElement[maxEntries];
  pDDStar=new dataElement[maxEntries];
  pDDStarPi=new dataElement[maxEntries];
  pDDStarPiWrongChannel=new dataElement[maxEntries];
  pDDStarPiPi=new dataElement[maxEntries];
  pDLNu=new dataElement[maxEntries];
  pDPiPiLNu=new dataElement[maxEntries];
  pDStarLNu=new dataElement[maxEntries];
  pDStarPiPiLNu=new dataElement[maxEntries];
  pDDStarPiCrossFeed=new dataElement[maxEntries];
  pOtherBB=new dataElement[maxEntries];
  mTrees=trees;


}


void sysDataGen::readTrees()
{


  for(int i=0;i<4;i++)
    {
      mTrees[i]->SetBranchAddress("mNu2",&mNu2);
      mTrees[i]->SetBranchAddress("tagCorr",&tagCorr);
      mTrees[i]->SetBranchAddress("bestBCharge",&bestBCharge);
      mTrees[i]->SetBranchAddress("systemCharge",&systemCharge);
      mTrees[i]->SetBranchAddress("leptonId",&leptonId);
      mTrees[i]->SetBranchAddress("foundAnyDDoubleStar",&foundAnyDDoubleStar);
      mTrees[i]->SetBranchAddress("sig_numPions",&sig_numPions);
      mTrees[i]->SetBranchAddress("sig_numKaons",&sig_numKaons);
      mTrees[i]->SetBranchAddress("sig_numPi0",&sig_numPi0);
      mTrees[i]->SetBranchAddress("sig_numBaryons",&sig_numBaryons);
      mTrees[i]->SetBranchAddress("sig_numLeptons",&sig_numLeptons);
      mTrees[i]->SetBranchAddress("sig_DLNu",&sig_DLNu);
      mTrees[i]->SetBranchAddress("sig_DPiLNu",&sig_DPiLNu);
      mTrees[i]->SetBranchAddress("sig_DPiPiLNu",&sig_DPiPiLNu);
      mTrees[i]->SetBranchAddress("sig_DStarLNu",&sig_DStarLNu);
      mTrees[i]->SetBranchAddress("sig_DStarPiLNu",&sig_DStarPiLNu);
      mTrees[i]->SetBranchAddress("sig_DStarPiPiLNu",&sig_DStarPiPiLNu);
      mTrees[i]->SetBranchAddress("numRecPions",&numRecPions);
      mTrees[i]->SetBranchAddress("dType",&dType);
      mTrees[i]->SetBranchAddress("dDecay",&dDecay);
      mTrees[i]->SetBranchAddress("mcIsDStar",&mcIsDStar);
      mTrees[i]->SetBranchAddress("mDnPi",&mDnPi);
      mTrees[i]->SetBranchAddress("pi1Mom",&pi1Mom);
      mTrees[i]->SetBranchAddress("pi2Mom",&pi2Mom);
      mTrees[i]->SetBranchAddress("bestD",&bestD);
      mTrees[i]->SetBranchAddress("dCharge",&dCharge);
      mTrees[i]->SetBranchAddress("recDType",&recDType);
      mTrees[i]->SetBranchAddress("mcDecaySignature",&mcDecaySignature);
      mTrees[i]->SetBranchAddress("recDecaySignature",&recDecaySignature);
      mTrees[i]->SetBranchAddress("mBTag",&mBTag);
      mTrees[i]->SetBranchAddress("deltaETag",&deltaETag);
      mTrees[i]->SetBranchAddress("logProb",&logProb);

      Long64_t  nevents=mTrees[i]->GetEntries();
      for(Long64_t i=0;i<nevents;i++)
	{
	  mTrees[i]->GetEntry(i);
	  if(recDecaySignature==0)
	    continue;
	  if(numRecPions!=1)
	    continue;
	  if(bestD==0)
	    continue;
	  if(mBTag<5.27)
	    continue;
	  if(logProb<-3.5)
	    continue;
	  if(fabs(deltaETag)>0.18)
	    continue;
	  if(mDnPi>3.5 || mDnPi<2.05)
	    continue;
	  if(pi1Mom<0.24)
	    continue;
	  if(bestBCharge!=((-1)*systemCharge))
	    continue;
	}


      float weight=getWeight();
      float eWeight=getWeightUncertainty();
      //make sure that there is no overlap
      int numClasses=0;
      ///continuum trees are the two last ones

      if(isContinuum()&& (i>1))
	{
	  numClasses++;
      //      pContinuum[ [currentCounts[iContinuum]++ ] ].weight=1.0;
	  pContinuum[currentCounts[iContinuum]++].mNu2=mNu2;
	  pContinuum[currentCounts[iContinuum]++].weight=weight;
	  pContinuum[currentCounts[iContinuum]++].eWeight=eWeight;
	}
      //iDDStar

      //DDStarPi

      //DDStarPiWrongChannel

      
      //DDStarPiPi

      //DLNu

      //DPiPiLNu

      //DStarLNu

      //DStarPiPiLNu


      //DDStarPiCrossFeed


      //OtherBB

    }

}

bool sysDataGen::isContinuum()
{

}
