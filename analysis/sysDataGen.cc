#include "sysDataGen.h"
#include "idx2.h"

const int maxEntries=1000000;
sysDataGen::sysDataGen(TTree** trees)
{

  for(int i=0;i<20;i++)
    {
      currentCounts[i]=0;
    }
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

  pData[iContinuum]=pContinuum;
  pData[iDDStar]=pDDStar;
  pData[iDDStarPi]=pDDStarPi;
  pData[iDDStarPiWrongChannel]=pDDStarPiWrongChannel;
  pData[iDDStarPiPi]=pDDStarPiPi;
  pData[iDLNu]=pDLNu;
  pData[iDPiPiLNu]=pDPiPiLNu;
  pData[iDStarLNu]=pDStarLNu;
  pData[iDStarPiPiLNu]=pDStarPiPiLNu;
  pData[iDDStarPiCrossFeed]=pDDStarPiCrossFeed;
  pData[iOtherBB]=pOtherBB;

  for(int i=0;i<20;i++)
    {
      histoNames[i]=new char[100];
    }
  sprintf(histoNames[iContinuum],"continuum");
  sprintf(histoNames[iDDStar],"DDStar");
  sprintf(histoNames[iDDStarPi],"DDStarPi");
  sprintf(histoNames[iDDStarPiWrongChannel],"DDStarPi_WrongChannel");
  sprintf(histoNames[iDDStarPiPi],"DDStarPiPi");
  sprintf(histoNames[iDLNu],"DLNu");
  sprintf(histoNames[iDPiPiLNu],"DPiPiLNu");
  sprintf(histoNames[iDStarLNu],"DStarLNu");
  sprintf(histoNames[iDStarPiPiLNu],"DStarPiPiLNu");
  sprintf(histoNames[iDDStarPiCrossFeed],"CrossFeed");
  sprintf(histoNames[iOtherBB],"OtherBB");
}



TH1D** sysDataGen::getTemplates(int channelIdx, int leptonId, char* channelString)
{
  int numPions=1;
  char buffer[200];
  TH1D** ret=new TH1D*[13];
  for(int i=0;i<11;i++)
    {
      sprintf(buffer,"%s_%dLept_%Pions_%s",histoNames[i],leptonId,numPions,channelString);
      ret[i]=new TH1D(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
      //for all channels we don't care

      if(i!=iDDStarPiCrossFeed && i!= iDDStarPiWrongChannel)
	{
	  for(int j=0;j<currentCounts[i];j++)
	    {
	      //for channelIdx==0 always true
	      bool accept=true;
	      if(channelIdx==1)
		{
		  //  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 0==mcIsDStar");
		  if(abs(pData[i][j].mcBCharge)!=1 || abs(pData[i][j].bestBCharge)!=1 ||pData[i][j].dCharge!=pData[i][j].mcDCharge || pData[i][j].dType==2 || abs(pData[i][j].mcDCharge)!=1 || 0!=pData[i][j].mcIsDStar)
		    accept=false;


		}
	      if(channelIdx==2)
		{
		  //		  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && mcIsDStar");

		  if(abs(pData[i][j].mcBCharge)!=1 || abs(pData[i][j].bestBCharge)!=1 ||pData[i][j].dCharge!=pData[i][j].mcDCharge || pData[i][j].dType!=2 || abs(pData[i][j].mcDCharge)!=1 || 1!=pData[i][j].mcIsDStar)
		    accept=false;

		}
	      if(channelIdx==3)
		{
		  //		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar==0");
		  if(abs(pData[i][j].mcBCharge)!=0 || abs(pData[i][j].bestBCharge)!=0 ||pData[i][j].dCharge!=pData[i][j].mcDCharge || pData[i][j].dType==2 || abs(pData[i][j].mcDCharge)!=0 || 0!=pData[i][j].mcIsDStar)
		    accept=false;
		}
	      if(channelIdx==4)
		{
		  //		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && mcIsDStar");
		  if(abs(pData[i][j].mcBCharge)!=0 || abs(pData[i][j].bestBCharge)!=0 ||pData[i][j].dCharge!=pData[i][j].mcDCharge || pData[i][j].dType!=2 || abs(pData[i][j].mcDCharge)!=0 || 1!=pData[i][j].mcIsDStar)
		    accept=false;
		}
	      if(accept)
		{
		  ret[i]->Fill(pData[i][j].mNu2,pData[i][j].weight);
		}
	      
	    }
	}
      else //wrong channel or cross feed-->use data from iDDStarPi (signal)
	{
	  //	  don't fill the wrong channel, feed down for 'all channels'
	  if(channelIdx>0)
	    {
	      for(int j=0;j<currentCounts[iDDStarPi];j++)
		{
		  if(i==iDDStarPiCrossFeed)
		    {
		      bool accept=false;
		      //same selection as signal, just different channel selection
		      if(channelIdx==1)
			{
			  //	  sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 1==mcIsDStar");
			  if(abs(pData[iDDStarPi]->mcBCharge)==1 && abs(pData[iDDStarPi]->bestBCharge)==1 && pData[iDDStarPi]->dCharge==pData[iDDStarPi]->mcDCharge&& pData[iDDStarPi]->dType!=2 && abs(pData[iDDStarPi]->mcDCharge)==1 && pData[iDDStarPi]->mcIsDStar==1)
			    accept=true;
			}
		      
		      //	  sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && 0==mcIsDStar");
		      
		      if(channelIdx==2)
			{
			  if(abs(pData[iDDStarPi]->mcBCharge)==1 && abs(pData[iDDStarPi]->bestBCharge)==1 && pData[iDDStarPi]->dCharge==pData[iDDStarPi]->mcDCharge&& pData[iDDStarPi]->dType==2 && abs(pData[iDDStarPi]->mcDCharge)==1 && pData[iDDStarPi]->mcIsDStar==0)
			    accept=true;
			}
		      //	  sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar==1");
		      if(channelIdx==3)
			{
		      if(abs(pData[iDDStarPi]->mcBCharge)==0 && abs(pData[iDDStarPi]->bestBCharge)==0 && pData[iDDStarPi]->dCharge==pData[iDDStarPi]->mcDCharge&& pData[iDDStarPi]->dType!=2 && abs(pData[iDDStarPi]->mcDCharge)==0 && pData[iDDStarPi]->mcIsDStar==1)
			accept=true;
			}
		      //	  sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && mcIsDStar==0");
		      if(channelIdx==4)
			{
			  if(abs(pData[iDDStarPi]->mcBCharge)==0 && abs(pData[iDDStarPi]->bestBCharge)==0 && pData[iDDStarPi]->dCharge==pData[iDDStarPi]->mcDCharge&& pData[iDDStarPi]->dType==2 && abs(pData[iDDStarPi]->mcDCharge)==0 && pData[iDDStarPi]->mcIsDStar==0)
			    accept=true;
			}
		      if(accept)
			{
			  ret[i]->Fill(pData[iDDStarPi][j].mNu2,pData[iDDStarPi][j].weight);
			}
		      
		    }
		  if(i== iDDStarPiWrongChannel)
		    {
		      bool accept=false;
		      //		  sprintf(chargeAndStarSelection,"  (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge || abs(mcDCharge)!=1  )");
		      if(channelIdx==1)
			{
			  if(abs(pData[iDDStarPi]->mcBCharge)!=1 || abs(pData[iDDStarPi]->bestBCharge)!=1 || pData[iDDStarPi]->dCharge!=pData[iDDStarPi]->mcDCharge||  abs(pData[iDDStarPi]->mcDCharge)!=1)
			    accept=true;
			}
		  //		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge || abs(mcDCharge)!=1  )");
		      if(channelIdx==2)
			{
			  if(abs(pData[iDDStarPi]->mcBCharge)!=1 || abs(pData[iDDStarPi]->bestBCharge)!=1 || pData[iDDStarPi]->dCharge!=pData[iDDStarPi]->mcDCharge|| abs(pData[iDDStarPi]->mcDCharge)!=1)
			    accept=true;
			}
		      //      sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge  || abs(mcDCharge)!=0  )");
		      if(channelIdx==3)
			{
			  if(abs(pData[iDDStarPi]->mcBCharge)!=0 || abs(pData[iDDStarPi]->bestBCharge)!=0 || pData[iDDStarPi]->dCharge!=pData[iDDStarPi]->mcDCharge|| abs(pData[iDDStarPi]->mcDCharge)!=0)
			    accept=true;
			}
		      
		      //		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge  || abs(mcDCharge)!=0  )");
		      if(channelIdx==4)
			{
			  if(abs(pData[iDDStarPi]->mcBCharge)!=0 || abs(pData[iDDStarPi]->bestBCharge)!=0 || pData[iDDStarPi]->dCharge!=pData[iDDStarPi]->mcDCharge|| abs(pData[iDDStarPi]->mcDCharge)!=0)
			    accept=true;
			}
		      if(accept)
			{
			  ret[i]->Fill(pData[iDDStarPi][j].mNu2,pData[iDDStarPi][j].weight);
			}
		    }
		}
	    }
	}
    }
  return ret;
}

void sysDataGen::readTrees()
{
  for(int i=0;i<4;i++)
    {
      cout <<"reading tree " << i <<endl;
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
      cout <<"done branching " <<endl;
      Long64_t  nevents=mTrees[i]->GetEntries();
      cout <<" tree has " << nevents <<" entries " <<endl;
      for(Long64_t j=0;j<nevents;j++)
	{
	  //	  if(j%10000==0)
	    cout <<"looking at event j " << endl;
	  //	  cout <<"reading event " << j <<endl;
	  mTrees[i]->GetEntry(j);
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

      if(isDDStar())
	{
	  numClasses++;
	  pDDStar[currentCounts[iDDStar]++].mNu2=mNu2;
	  pDDStar[currentCounts[iDDStar]++].weight=weight;
	  pDDStar[currentCounts[iDDStar]++].eWeight=eWeight;
	}
      //DDStarPi

      if(isDDStarPi())
	{
	  numClasses++;
	  pDDStarPi[currentCounts[iDDStarPi]++].mNu2=mNu2;
	  pDDStarPi[currentCounts[iDDStarPi]++].weight=weight;
	  pDDStarPi[currentCounts[iDDStarPi]++].eWeight=eWeight;
	}

      /// wrong channel and Cross-feed will be constructed when we build the 
      // channel specific histograms since they are channel dependent
      //DDStarPiWrongChannel
      //DDStarPiCrossFeed
      
      //DDStarPiPi
      if(isDDStarPiPi())
	{
	  numClasses++;
	  pDDStarPiPi[currentCounts[iDDStarPiPi]++].mNu2=mNu2;
	  pDDStarPiPi[currentCounts[iDDStarPiPi]++].weight=weight;
	  pDDStarPiPi[currentCounts[iDDStarPiPi]++].eWeight=eWeight;
	}
      //DLNu
      if(isDLNu())
	{
	  numClasses++;
	  pDLNu[currentCounts[iDLNu]++].mNu2=mNu2;
	  pDLNu[currentCounts[iDLNu]++].weight=weight;
	  pDLNu[currentCounts[iDLNu]++].eWeight=eWeight;
	}
      //DPiPiLNu
      if(isDPiPiLNu())
	{
	  numClasses++;
	  pDPiPiLNu[currentCounts[iDPiPiLNu]++].mNu2=mNu2;
	  pDPiPiLNu[currentCounts[iDPiPiLNu]++].weight=weight;
	  pDPiPiLNu[currentCounts[iDPiPiLNu]++].eWeight=eWeight;
	}
      //DStarLNu
      if(isDStarLNu())
	{
	  numClasses++;
	  pDStarLNu[currentCounts[iDStarLNu]++].mNu2=mNu2;
	  pDStarLNu[currentCounts[iDStarLNu]++].weight=weight;
	  pDStarLNu[currentCounts[iDStarLNu]++].eWeight=eWeight;
	}
      //DStarPiPiLNu
      if(isDStarPiPiLNu())
	{
	  numClasses++;
	  pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]++].mNu2=mNu2;
	  pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]++].weight=weight;
	  pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]++].eWeight=eWeight;
	}

      //OtherBB
      if(isOtherBB())
	{
	  numClasses++;
	  pOtherBB[currentCounts[iOtherBB]++].mNu2=mNu2;
	  pOtherBB[currentCounts[iOtherBB]++].weight=weight;
	  pOtherBB[currentCounts[iOtherBB]++].eWeight=eWeight;
	}
      if(numClasses>1)
	cout<<"event " << j <<"  seems to be in more than one class " <<endl;

    }
}

bool sysDataGen::isContinuum()
{
  //all of the uds and charm tree results are accepted
  return true;
}
      //iDDStar
bool sysDataGen::isDDStar()
{
  if(!foundAnyDDoubleStar)
    return false;
  if(sig_numPions>0)
    return false;
  if(sig_numKaons>0 || sig_numPi0>0 || sig_numBaryons>0)
    return false;
  if(sig_numLeptons!=1)
    return false;
  if(sig_DLNu || sig_DPiLNu || sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu || sig_DStarPiPiLNu)
    return false;

  return true;
  
}
//DDStarPi
bool sysDataGen::isDDStarPi()
{
  if(!foundAnyDDoubleStar)
    return false;
  if(sig_numPions!=1)
    return false;
  if(sig_numKaons>0 || sig_numPi0>0 || sig_numBaryons>0)
    return false;
  if(sig_numLeptons!=1)
    return false;
  if(sig_DLNu || sig_DPiLNu || sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu || sig_DStarPiPiLNu)
    return false;

  return true;
}



///??
//DDStarPiWrongChannel 
// bool sysDataGen::isDDStarPiWrongChannel()
// {
//   if(!foundAnyDDoubleStar)
//     return false;
//   if(sig_numPions!=1)
//     return false;
//   if(sig_numKaons>0 || sig_numPi0>0 || sig_numBaryons>0)
//     return false;
//   if(sig_numLeptons!=1)
//     return false;
//   if(sig_DLNu || sig_DPiLNu || sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu || sig_DStarPiPiLNu)
//     return false;
// 
// 
//   return true;
// }
// 
///
// bool sysDataGen::isDDStarPiCrossFeed()
// {
//   if(!foundAnyDDoubleStar)
//     return false;
//   if(sig_numPions!=1)
//     return false;
//   if(sig_numKaons>0 || sig_numPi0>0 || sig_numBaryons>0)
//     return false;
//   if(sig_numLeptons!=1)
//     return false;
//   if(sig_DLNu || sig_DPiLNu || sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu || sig_DStarPiPiLNu)
//     return false;
// 
// }      
//DDStarPiPi
bool sysDataGen::isDDStarPiPi()
{
  if(!foundAnyDDoubleStar)
    return false;
  if(sig_numPions!=2)
    return false;
  if(sig_numKaons>0 || sig_numPi0>0 || sig_numBaryons>0)
    return false;
  if(sig_numLeptons!=1)
    return false;
  if(sig_DLNu || sig_DPiLNu || sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu || sig_DStarPiPiLNu)
    return false;


  return true;
}
 
//DLNu
bool sysDataGen::isDLNu()
{
  if(sig_DLNu)
    return true;

  return false;

}


bool sysDataGen::isDPiPiLNu()
{
  if(sig_DPiPiLNu && !sig_DLNu && !sig_DPiLNu)
    return true;


  return false;
}

bool sysDataGen::isDStarLNu()
{
  if(sig_DStarLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu)
    return true;

  return false;
}

bool sysDataGen::isDStarPiPiLNu()
{
  if( sig_DStarPiPiLNu && !sig_DLNu && !sig_DPiLNu&& !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu)
    return true;


  return false;

}



bool sysDataGen::isOtherBB()
{
  if((!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu|| sig_DStarPiPiLNu) && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)
    return true;

  return false;

}

float sysDataGen::getWeight()
{

  return 1.0;
}

float sysDataGen::getWeightUncertainty()
{
  return 0.0;
}
