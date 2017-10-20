#include "sysDataGen.h"
#include "idx2.h"
#include "TCanvas.h"
#include "TRandom3.h"

const long int maxEntries=40000;


sysDataGen::sysDataGen(TTree** trees)
{

  cout <<" we have " << numBins[1] << " bins " <<endl;


  mixedFactor=1.0/0.852255195;
  chargedFactor=1.0/0.9266338302;

  cout <<"mixed Factor: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;

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
  cout <<"mixed Factor2: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;
}


TH1F** sysDataGen::getTemplates(int channelIdx, int leptonId, char* channelString, TH1F** components, bool doSysStudies)
{
  //generate the gaussian weights
  if(doSysStudies)
    {
      for(int i=0;i<6480;i++)
	{
	  ChargedWeights[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<12;i++)
	{
	  KsWeights[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<24;i++)
	{
	  pi0Weights[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<576;i++)
	{
	  D_FFWeights[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<680;i++)
	{
	  Dds_FFWeights[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<26;i++)
	{
	  D_decTypeWeights[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<12;i++)
	{
	  B_decWeights[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<100;i++)
	{
	  lumiWeight[i]=rnd.Gaus(0,1.0);
	}
      for(int i=0;i<4;i++)
	{
	  tagIdWeight[i]=rnd.Gaus(0,1.0);
	}

      const float sigmaDelta=0.0035;
      const float effMC=0.8;
      const float sigmaW=sigmaDelta/effMC;
      const float smearPerTrack=sigmaW*rnd.Gaus(0,1.0);
      float newWeight=1.0;

      for(int i=0;i<200;i++)
	{
	  //like Robin assyme 0.8 efficiency and 0.35% rel uncertainty
	  trackWeights[i]=newWeight;
	  newWeight*=(1+smearPerTrack);
	}
    }

  cout <<"mixed Factor3: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;
  cout <<"should we really? " <<endl;
  int numPions=1;
  char buffer[200];
  TH1F** ret=new TH1F*[20];
  TCanvas* c;//=new TCanvas();
  for(int i=0;i<11;i++)
    {
      bool foundDDoubleStar=false;
      if(i==iDDStar || i==iDDStarPi || i==iDDStarPiWrongChannel || i==iDDStarPiPi || i==iDDStarPiCrossFeed)
	foundDDoubleStar=true;

      cout <<"histo name: "<< histoNames[i] <<endl;
      cout <<"channel string: "<< channelString <<endl;
      sprintf(buffer,"%s_%dLept_%dPions_%s",histoNames[i],leptonId,numPions,channelString);
      cout <<"creating " << buffer <<" with : "<< numBins[channelIdx] <<" lower Cut: "<< lowerCut[channelIdx] << " upperCut: "<< upperCut[channelIdx] <<endl;
      ret[i]=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
      ret[i]->Sumw2();
      ret[i]->SetFillStyle(1001);
      ret[i]->SetFillColor(glColorTable[i]->GetNumber());
      ///components are additionally sorted by files (4)
      for(int iF=0;iF<4;iF++)
	{
	  sprintf(buffer,"histo_If_%d_b_%d_numPions_%d_leptonId_%d_%s",iF,i,numPions,leptonId,channelString);
	  components[iF*11+i]=new TH1F(buffer,buffer,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	  components[iF*11+i]->SetFillStyle(1001);
	  components[iF*11+i]->SetFillColor(glColorTable[i]->GetNumber());
	}
      //      cout <<"trying to draw 1 " << i << endl;
      //      sprintf(buffer,"cnvs%d_%s",i,channelString);
      //      c=new TCanvas(buffer,buffer);
      //      delete c;
      //      sprintf(buffer,"cnvsf%d_%s",i,channelString);
      //      c=new TCanvas(buffer,buffer);
      //      cout <<"created canvas " << buffer <<endl;
      //      ret[i]->Draw();
      //      sprintf(buffer,"histoDraw_2_%d.png",i);
      //      cout <<"done " << endl;
      //      c->SaveAs(buffer);
      //      c->Clear();//9850
      cout <<"created histo " << i <<endl;
      //      delete c;
      cout <<"done delete 1 " << endl;
      //for all channels we don't care
      //

      if(i!=iDDStarPiCrossFeed && i!= iDDStarPiWrongChannel)
	{
	  if(currentCounts[i]>=maxEntries)
	    {
	      cout <<"too many entries" <<i << endl;
	      exit(0);
	    }
	  cout <<" running over " << currentCounts[i] <<endl;
	  for(int j=0;j<currentCounts[i];j++)
	    {
	      int iF=pData[i][j].fileIndex;
	      //	      if(pData[i][j].mNu2==0)
	      //	      if(i==1)
	      //		cout <<"looking at data element " << j << " for channel " << i <<" mNu2: "<< pData[i][j].mNu2 <<endl;
	      //	      cout <<" mNu2:" << pData[i][j].mNu2 <<endl;
	      //for channelIdx==0 always true
	      bool accept=true;
	      if(leptonId!=0 && (abs(pData[i][j].leptonId)!=leptonId))
		{
		  accept=false;
		  continue;
		}
	      if(channelIdx==1)
		{
		  //  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 0==mcIsDStar");
		  if(abs(pData[i][j].bestBCharge)!=1 ||abs(pData[i][j].dCharge)!=1|| pData[i][j].dType==2 )
		    accept=false;
		  //for the signal, the true mc has also be fulfilled
		  if(i==iDDStarPi)
		    {
		      if(abs(pData[i][j].mcBCharge)!=1  ||pData[i][j].dCharge!=pData[i][j].mcDCharge || abs(pData[i][j].mcDCharge)!=1 || 0!=pData[i][j].mcIsDStar)
			accept=false;
		    }

		}
	      if(channelIdx==2)
		{
		  //		  sprintf(chargeAndStarSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && mcIsDStar");
		  if( abs(pData[i][j].bestBCharge)!=1 || abs(pData[i][j].dCharge)!=1 ||  pData[i][j].dType!=2)
			accept=false;
		  if(i==iDDStarPi)
		    {
		      if(abs(pData[i][j].mcBCharge)!=1 || abs(pData[i][j].bestBCharge)!=1 || pData[i][j].dCharge!=pData[i][j].mcDCharge || pData[i][j].dType!=2 || abs(pData[i][j].mcDCharge)!=1 || 1!=pData[i][j].mcIsDStar)
		    accept=false;
		    }
		}
	      if(channelIdx==3)
		{
		  if( abs(pData[i][j].bestBCharge)!=0 || abs(pData[i][j].dCharge)!=0 ||  pData[i][j].dType==2)
			accept=false;
		  //		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar==0");
		  if(i==iDDStarPi)
		    {
		      if(abs(pData[i][j].mcBCharge)!=0 ||pData[i][j].dCharge!=pData[i][j].mcDCharge || abs(pData[i][j].mcDCharge)!=0 || 0!=pData[i][j].mcIsDStar)
			accept=false;
		    }
		}
	      if(channelIdx==4)
		{
		  if( abs(pData[i][j].bestBCharge)!=0 || abs(pData[i][j].dCharge)!=0 ||  pData[i][j].dType!=2)
			accept=false;
		  //		  sprintf(chargeAndStarSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && mcIsDStar");
		  if(i==iDDStarPi)
		    {
		      if(abs(pData[i][j].mcBCharge)!=0  ||pData[i][j].dCharge!=pData[i][j].mcDCharge  || abs(pData[i][j].mcDCharge)!=0 || 1!=pData[i][j].mcIsDStar)
			accept=false;
		    }
		}
	      if(accept)
		{
		  if(i==1)
		    {
		      //		      cout <<"filling with " << pData[i][j].mNu2 << " and weight; "<< pData[i][j].weight <<endl;
		    }
		  float weight=pData[i][j].weight;
		  if(doSysStudies)
		    {
		      weight=getWeightWithError(pData[i][j],foundDDoubleStar);
		    }
		  //		  cout <<"filling with weight: "<< weight<<endl;
		  ret[i]->Fill(pData[i][j].mNu2,weight);
		  components[iF*11+i]->Fill(pData[i][j].mNu2,weight);

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
		  int iF=pData[iDDStarPi][j].fileIndex;
		  if(i==iDDStarPiCrossFeed)
		    {
		      cout <<"looking at xfeed " << endl;
		      cout <<" mcBCharge: "<< pData[iDDStarPi][j].mcBCharge;
		      cout <<", bestBcharge: "<< pData[iDDStarPi][j].mcBCharge;
			cout <<", dcharge: "<< pData[iDDStarPi][j].dCharge;
			cout <<", mcdcharge: "<< pData[iDDStarPi][j].mcDCharge;
			cout <<", dTypee: "<< pData[iDDStarPi][j].dType;
			cout <<", mcIsDStar: "<< pData[iDDStarPi][j].mcIsDStar<<endl;
		      bool accept=false;
		      //same selection as signal, just different channel selection
		      if(channelIdx==1)
			{
			  //	  sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==1  && 1==mcIsDStar");
			  if(abs(pData[iDDStarPi][j].mcBCharge)==1 && abs(pData[iDDStarPi][j].bestBCharge)==1 && pData[iDDStarPi][j].dCharge==pData[iDDStarPi][j].mcDCharge&& pData[iDDStarPi][j].dType!=2 && abs(pData[iDDStarPi][j].mcDCharge)==1 && pData[iDDStarPi][j].mcIsDStar==1)
			    accept=true;
			}
		      
		      //	  sprintf(feedDownSelection," abs(mcBCharge)==1 && abs(bestBCharge)==1 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==1  && 0==mcIsDStar");
		      
		      if(channelIdx==2)
			{
			  if(abs(pData[iDDStarPi][j].mcBCharge)==1 && abs(pData[iDDStarPi][j].bestBCharge)==1 && pData[iDDStarPi][j].dCharge==pData[iDDStarPi][j].mcDCharge&& pData[iDDStarPi][j].dType==2 && abs(pData[iDDStarPi][j].mcDCharge)==1 && pData[iDDStarPi][j].mcIsDStar==0)
			    accept=true;
			}
		      //	  sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType!=2 && abs(mcDCharge)==0  && mcIsDStar==1");
		      if(channelIdx==3)
			{
			  if(abs(pData[iDDStarPi][j].mcBCharge)==0 && abs(pData[iDDStarPi][j].bestBCharge)==0 && pData[iDDStarPi][j].dCharge==pData[iDDStarPi][j].mcDCharge&& pData[iDDStarPi][j].dType!=2 && abs(pData[iDDStarPi][j].mcDCharge)==0 && pData[iDDStarPi][j].mcIsDStar==1)
			accept=true;
			}
		      //	  sprintf(feedDownSelection," abs(mcBCharge)==0 && abs(bestBCharge)==0 && dCharge==mcDCharge && dType==2 && abs(mcDCharge)==0  && mcIsDStar==0");
		      if(channelIdx==4)
			{
			  if(abs(pData[iDDStarPi][j].mcBCharge)==0 && abs(pData[iDDStarPi][j].bestBCharge)==0 && pData[iDDStarPi][j].dCharge==pData[iDDStarPi][j].mcDCharge&& pData[iDDStarPi][j].dType==2 && abs(pData[iDDStarPi][j].mcDCharge)==0 && pData[iDDStarPi][j].mcIsDStar==0)
			    accept=true;
			}
		      if(accept)
			{
			  float weight=pData[iDDStarPi][j].weight;
			  if(doSysStudies)
			    {
			      weight=getWeightWithError(pData[iDDStarPi][j],foundDDoubleStar);
			    }
		  cout <<"filling with weight: "<< weight<<endl;
			  ret[i]->Fill(pData[iDDStarPi][j].mNu2,weight);
			  components[iF*11+i]->Fill(pData[iDDStarPi][j].mNu2,weight);
			}
		      
		    }
		  if(i==iDDStarPiWrongChannel)
		    {
		      bool accept=false;
		      //		  sprintf(chargeAndStarSelection,"  (abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge || abs(mcDCharge)!=1  )");
		      if(channelIdx==1)
			{
			  if(abs(pData[iDDStarPi][j].mcBCharge)!=1 || abs(pData[iDDStarPi][j].bestBCharge)!=1 || pData[iDDStarPi][j].dCharge!=pData[iDDStarPi][j].mcDCharge||  abs(pData[iDDStarPi][j].mcDCharge)!=1)
			    accept=true;
			  //but check that we would actually accept this event
			  if(abs(pData[iDDStarPi][j].bestBCharge)!=1 ||abs(pData[iDDStarPi][j].dCharge)!=1|| pData[iDDStarPi][j].dType==2 )
			    accept=false;

			}
		  //		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=1 || abs(bestBCharge)!=1 || dCharge!=mcDCharge || abs(mcDCharge)!=1  )");
		      if(channelIdx==2)
			{
			  if(abs(pData[iDDStarPi][j].mcBCharge)!=1 || abs(pData[iDDStarPi][j].bestBCharge)!=1 || pData[iDDStarPi][j].dCharge!=pData[iDDStarPi][j].mcDCharge|| abs(pData[iDDStarPi][j].mcDCharge)!=1)
			    accept=true;
			  //but check that we would actually accept this event
			  if(abs(pData[iDDStarPi][j].bestBCharge)!=1 ||abs(pData[iDDStarPi][j].dCharge)!=1|| pData[iDDStarPi][j].dType!=2 )
			    accept=false;
			}
		      //      sprintf(chargeAndStarSelection,"   (abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge  || abs(mcDCharge)!=0  )");
		      if(channelIdx==3)
			{
			  if(abs(pData[iDDStarPi][j].mcBCharge)!=0 || abs(pData[iDDStarPi][j].bestBCharge)!=0 || pData[iDDStarPi][j].dCharge!=pData[iDDStarPi][j].mcDCharge|| abs(pData[iDDStarPi][j].mcDCharge)!=0)
			    accept=true;
			  //but check that we would actually accept this event
			  if(abs(pData[iDDStarPi][j].bestBCharge)!=0 ||abs(pData[iDDStarPi][j].dCharge)!=0|| pData[iDDStarPi][j].dType==2)
			    accept=false;
			}
		      
		      //		  sprintf(chargeAndStarSelection,"  ( abs(mcBCharge)!=0 || abs(bestBCharge)!=0 || dCharge!=mcDCharge  || abs(mcDCharge)!=0  )");
		      if(channelIdx==4)
			{
			  if(abs(pData[iDDStarPi][j].mcBCharge)!=0 || abs(pData[iDDStarPi][j].bestBCharge)!=0 || pData[iDDStarPi][j].dCharge!=pData[iDDStarPi][j].mcDCharge|| abs(pData[iDDStarPi][j].mcDCharge)!=0)
			    accept=true;
			  //but check that we would actually accept this event
			  if(abs(pData[iDDStarPi][j].bestBCharge)!=0 ||abs(pData[iDDStarPi][j].dCharge)!=0|| pData[iDDStarPi][j].dType!=2 )
			    accept=false;
			}
		      if(accept)
			{
			  float weight=pData[iDDStarPi][j].weight;
			  if(doSysStudies)
			    {
			      weight=getWeightWithError(pData[iDDStarPi][j],foundDDoubleStar);
			    }
		  cout <<"filling with weight: "<< weight<<endl;
			  ret[i]->Fill(pData[iDDStarPi][j].mNu2,weight);
			  components[iF*11+i]->Fill(pData[iDDStarPi][j].mNu2,weight);
			}
		    }
		}
	    }
	}
      if(i==1)
	{
	  cout <<" trying to save before doing the others.." <<endl;
	  //	  ret[1]->Draw();
	  //	  c.SaveAs("tmpTst.png");
	}
    }
  cout <<"returning histos " << endl;
  cout <<" histo 1 has " << ret[1]->GetNbinsX()<<endl;
  char b[200];

  for(int i=0;i<11;i++)
    {
      //      cout <<"saving in sysData Gen " << i <<endl;
      //      sprintf(b,"cnvs2_%d_%s.png",i,channelString);
      //      c=new TCanvas(b,b);
      //      ret[i]->Draw();
      //      sprintf(b,"sysH_%d.png",i);
      //      c->SaveAs(b);

      //      cout <<"done " <<endl;
      //           delete c;
      //      cout <<"done delete" <<endl;
    }
  cout <<"about to return ... " << endl;
  return ret;
}

void sysDataGen::readTrees()
{
  cout <<"mixed Factor read Trees: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;
  for(int i=0;i<4;i++)
    {
      cout <<"mixed Factor read T: "<<i <<" fact: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;
      cout <<"reading tree " << i <<endl;
      mTrees[i]->SetBranchAddress("mNu2",mNu2);
      mTrees[i]->SetBranchAddress("mNu2Counter",&mNu2Counter);
      mTrees[i]->SetBranchAddress("tagCorr",&tagCorr);
      mTrees[i]->SetBranchAddress("B_DecayCorr",&B_DecayCorr);
      mTrees[i]->SetBranchAddress("sB_DecayCorr",&sB_DecayCorr);
      mTrees[i]->SetBranchAddress("B_decayType",&B_decayType);
      mTrees[i]->SetBranchAddress("experiment",&experiment);
      mTrees[i]->SetBranchAddress("CrossSectionLumiCorrection",&CrossSectionLumiCorrection);
      mTrees[i]->SetBranchAddress("sCrossSectionLumiCorrection",&sCrossSectionLumiCorrection);
      ///for x-check
      //           mTrees[i]->SetBranchAddress("FFCorrection",&FFCorrection);
      //           mTrees[i]->SetBranchAddress("PIDCorrection",&PIDCorrection);
	   //           mTrees[i]->SetBranchAddress("Correction",&FFCorrection);
      //      mTrees[i]->SetBranchAddress("sFFCorrection",&sFFCorrection);
      //      mTrees[i]->SetBranchAddress("sRelTrackingCorrection",&sRelTrackingCorrection);
      mTrees[i]->SetBranchAddress("bestBCharge",&bestBCharge);
      mTrees[i]->SetBranchAddress("mcBCharge",&mcBCharge);
      mTrees[i]->SetBranchAddress("systemCharge",&systemCharge);
      mTrees[i]->SetBranchAddress("leptonId",leptonId);
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
      mTrees[i]->SetBranchAddress("numRecPions",numRecPions);
      mTrees[i]->SetBranchAddress("dType",dType);
      mTrees[i]->SetBranchAddress("dDecay",dDecay);
      mTrees[i]->SetBranchAddress("mcIsDStar",&mcIsDStar);
      mTrees[i]->SetBranchAddress("mDnPi",mDnPi);
      mTrees[i]->SetBranchAddress("pi1Mom",pi1Mom);
      mTrees[i]->SetBranchAddress("pi2Mom",pi2Mom);
      mTrees[i]->SetBranchAddress("bestD",bestD);
      mTrees[i]->SetBranchAddress("dCharge",dCharge);
      mTrees[i]->SetBranchAddress("mcDCharge",&mcDCharge);
      mTrees[i]->SetBranchAddress("recDType",recDType);
      mTrees[i]->SetBranchAddress("mcDecaySignature",&mcDecaySignature);
      mTrees[i]->SetBranchAddress("recDecaySignature",&recDecaySignature);
      mTrees[i]->SetBranchAddress("mBTag",&mBTag);
      mTrees[i]->SetBranchAddress("deltaETag",&deltaETag);
      mTrees[i]->SetBranchAddress("logProb",&logProb);
      mTrees[i]->SetBranchAddress("tagId",&tagId);

      mTrees[i]->SetBranchAddress("KsCorrection",KsCorrection);
      mTrees[i]->SetBranchAddress("sKsCorrection",sKsCorrection);
      mTrees[i]->SetBranchAddress("KsCorrectionCounter",&KsCorrectionCounter);

      mTrees[i]->SetBranchAddress("ChargedCorrection",ChargedCorrection);
      mTrees[i]->SetBranchAddress("sChargedCorrection",sChargedCorrection);
      mTrees[i]->SetBranchAddress("ChargedCorrectionCounter",&ChargedCorrectionCounter);

      mTrees[i]->SetBranchAddress("Pi0Correction",Pi0Correction);
      mTrees[i]->SetBranchAddress("sPi0Correction",sPi0Correction);
      mTrees[i]->SetBranchAddress("Pi0CorrectionCounter",&Pi0CorrectionCounter);


      mTrees[i]->SetBranchAddress("D_DecayCorr",D_DecayCorr);
      mTrees[i]->SetBranchAddress("sD_DecayCorr",sD_DecayCorr);
      mTrees[i]->SetBranchAddress("D_DecayCorrCounter",&D_DecayCorrCounter);


      mTrees[i]->SetBranchAddress("FFDCorrection",FFDCorrection);
      mTrees[i]->SetBranchAddress("sFFDCorrection",sFFDCorrection);
      mTrees[i]->SetBranchAddress("FFDCorrectionCounter",&FFDCorrectionCounter);


      mTrees[i]->SetBranchAddress("FFDdsCorrection",FFDdsCorrection);
      mTrees[i]->SetBranchAddress("sFFDdsCorrection",sFFDdsCorrection);
      mTrees[i]->SetBranchAddress("FFDdsCorrectionCounter",&FFDdsCorrectionCounter);


      mTrees[i]->SetBranchAddress("ChargedCorrThetaBin",ChargedCorrThetaBin);
      mTrees[i]->SetBranchAddress("ChargedCorrMomBin",ChargedCorrMomBin);
      mTrees[i]->SetBranchAddress("ChargedCorrSVDBin",ChargedCorrSVDBin);
      mTrees[i]->SetBranchAddress("ChargedCorrMisIdType",ChargedCorrMisIdType);

      mTrees[i]->SetBranchAddress("KsCorrMomBin",KsCorrMomBin);
      mTrees[i]->SetBranchAddress("KsCorrThetaBin",KsCorrThetaBin);

      mTrees[i]->SetBranchAddress("pi0MomBin",pi0MomBin);
      mTrees[i]->SetBranchAddress("D_pBinFF",D_pBinFF);
      mTrees[i]->SetBranchAddress("D_q2BinFF",D_q2BinFF);
      mTrees[i]->SetBranchAddress("D_TypeFF",D_TypeFF);

      mTrees[i]->SetBranchAddress("Dds_wBinFF",Dds_wBinFF);
      mTrees[i]->SetBranchAddress("Dds_cosTBinFF",Dds_cosTBinFF);
      mTrees[i]->SetBranchAddress("Dds_TypeFF",Dds_TypeFF);
      mTrees[i]->SetBranchAddress("D_decType",D_decType);
      cout <<"done branching " <<endl;
      cout <<"mixed Factor read T: "<<i <<" fact: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;
      cout <<endl;
      Long64_t  nevents=mTrees[i]->GetEntries();
      //            nevents/=100;
      cout <<" tree has " << nevents <<" entries " <<endl;
      cout <<" so ... what now? " << endl;
      cout <<"mixed Factor read T: "<<i <<" fact: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;
      for(Long64_t j=0;j<nevents;j++)
	{
	  for(int iC=0;iC<20;iC++)
	    {
	      if(currentCounts[iC]>=maxEntries)
		{
		  cout <<"got too many entries in " << i <<endl;
		  exit(0);
		}
	    }
	  mTrees[i]->GetEntry(j);
	      //	      cout <<"looking at "<< mNu2Counter <<" candidates " <<endl;
	  if(j%100000==0)
	    //	      if(j%100==0)
	    cout <<"looking at event " << j << endl;
	  for(int iM=0;iM<mNu2Counter;iM++)
	    {

	      //	  cout <<"reading event " << j <<endl;

	      if(recDecaySignature==0)
		continue;
	      if(numRecPions[iM]!=1)
		continue;
	      if(bestD[iM]==0)
		continue;
	      if(mBTag<5.27)
		continue;
	      if(logProb<-3.5)
		continue;
	      if(fabs(deltaETag)>0.18)
		continue;
	      if(mDnPi[iM]>3.5 || mDnPi[iM]<2.05)
		continue;
	      if(pi1Mom[iM]<0.24)
		continue;
	      if(bestBCharge!=((-1)*systemCharge))
		continue;
	      if(mixedFactor<0.4)
		{
		  cout <<"event " << j << " has mixed factor "<< mixedFactor <<endl;
		  cout <<" current Counter cont: "<< currentCounts[iContinuum]<<endl;
		  cout <<" current Counter ddstar: "<< currentCounts[iDDStar]<<endl;
		  cout <<" current Counter ddstarpi: "<< currentCounts[iDDStarPi]<<endl;
		  cout <<" current Counter iddstarpipi: "<< currentCounts[iDDStarPiPi]<<endl;
		  cout <<" current Counter wrong channel: "<< currentCounts[iDDStarPiWrongChannel]<<endl;
		  cout <<" current Counter dlnu: "<< currentCounts[iDLNu]<<endl;
		  cout <<" current Counter dpipilnu: "<< currentCounts[iDPiPiLNu]<<endl;
		  cout <<" current Counter dstarlnu: "<< currentCounts[iDStarLNu]<<endl;
		  cout <<" current Counter dstarpipi: "<< currentCounts[iDStarPiPiLNu]<<endl;
		  cout <<" current Counter cross feed: "<< currentCounts[iDDStarPiCrossFeed]<<endl;
		  cout <<" current Counter other: "<< currentCounts[iOtherBB]<<endl;
		}
	      
	      
	      if(chargedFactor<0.4)
		{
		  cout <<"event " << j << " has charged factor "<< chargedFactor <<endl;
		  cout <<"event " << j << " has mixed factor "<< mixedFactor <<endl;
		  cout <<" current Counter cont: "<< currentCounts[iContinuum]<<endl;
		  cout <<" current Counter ddstar: "<< currentCounts[iDDStar]<<endl;
		  cout <<" current Counter ddstarpi: "<< currentCounts[iDDStarPi]<<endl;
		  cout <<" current Counter iddstarpipi: "<< currentCounts[iDDStarPiPi]<<endl;
		  cout <<" current Counter wrong channel: "<< currentCounts[iDDStarPiWrongChannel]<<endl;
		  cout <<" current Counter dlnu: "<< currentCounts[iDLNu]<<endl;
		  cout <<" current Counter dpipilnu: "<< currentCounts[iDPiPiLNu]<<endl;
		  cout <<" current Counter dstarlnu: "<< currentCounts[iDStarLNu]<<endl;
		  cout <<" current Counter dstarpipi: "<< currentCounts[iDStarPiPiLNu]<<endl;
		  cout <<" current Counter cross feed: "<< currentCounts[iDDStarPiCrossFeed]<<endl;
		  cout <<" current Counter other: "<< currentCounts[iOtherBB]<<endl;
	    }
	      //	  cout <<"mixed Factor read T: "<<i <<" fact: "<< mixedFactor <<" chargedFactor: "<< chargedFactor <<endl;
	      //	      cout <<"calling get Weight .. " << endl;
	      float weight=getWeight(i,foundAnyDDoubleStar);
	      //meaningless now, since we compute on the fly given the bins
	      float eWeight=0.0;

	      //make sure that there is no overlap
	      int numClasses=0;
	      ///continuum trees are the two last ones

	      if(isContinuum()&& (i>1))
		{
		  numClasses++;
		  //      pContinuum[ [currentCounts[iContinuum]++ ] ].weight=1.0;
		  pContinuum[currentCounts[iContinuum]].mNu2=mNu2[iM];
		  pContinuum[currentCounts[iContinuum]].weight=weight;
		  pContinuum[currentCounts[iContinuum]].mcBCharge=mcBCharge;
		  pContinuum[currentCounts[iContinuum]].bestBCharge=bestBCharge;
		  pContinuum[currentCounts[iContinuum]].dCharge=dCharge[iM];
		  pContinuum[currentCounts[iContinuum]].mcDCharge=mcDCharge;
		  pContinuum[currentCounts[iContinuum]].dType=dType[iM];
		  pContinuum[currentCounts[iContinuum]].mcIsDStar=mcIsDStar;
		  pContinuum[currentCounts[iContinuum]].leptonId=leptonId[iM];
		  pContinuum[currentCounts[iContinuum]].eWeight=eWeight;
		  fillWeightDetails(pContinuum[currentCounts[iContinuum]],i);
		  pContinuum[currentCounts[iContinuum]++].fileIndex=i;

		  //	      cout <<"is continuum" <<endl;
		}

	  //iDDStar
	      if(i<2)//the on_res trees
		{
		  if(isDDStar())
		    {
		      numClasses++;
		      pDDStar[currentCounts[iDDStar]].mNu2=mNu2[iM];
		      pDDStar[currentCounts[iDDStar]].weight=weight;
		      pDDStar[currentCounts[iDDStar]].mcBCharge=mcBCharge;
		      pDDStar[currentCounts[iDDStar]].bestBCharge=bestBCharge;
		      pDDStar[currentCounts[iDDStar]].dCharge=dCharge[iM];
		      pDDStar[currentCounts[iDDStar]].mcDCharge=mcDCharge;
		      pDDStar[currentCounts[iDDStar]].dType=dType[iM];
		      pDDStar[currentCounts[iDDStar]].mcIsDStar=mcIsDStar;
		      pDDStar[currentCounts[iDDStar]].leptonId=leptonId[iM];
		      pDDStar[currentCounts[iDDStar]].eWeight=eWeight;
		      fillWeightDetails(pDDStar[currentCounts[iDDStar]],i);
		      pDDStar[currentCounts[iDDStar]++].fileIndex=i;
		      //		  cout <<"count for ddstar: "<< currentCounts[iDDStar] <<" " << iDDStar <<endl;
		      ///		  cout <<"is DDStar" <<endl;
		    }
		  //DDStarPi
		  
		  if(isDDStarPi())
		    {
		      numClasses++;
		      pDDStarPi[currentCounts[iDDStarPi]].mNu2=mNu2[iM];
		      pDDStarPi[currentCounts[iDDStarPi]].weight=weight;
		      pDDStarPi[currentCounts[iDDStarPi]].mcBCharge=mcBCharge;
		      pDDStarPi[currentCounts[iDDStarPi]].bestBCharge=bestBCharge;
		      pDDStarPi[currentCounts[iDDStarPi]].dCharge=dCharge[iM];
		      pDDStarPi[currentCounts[iDDStarPi]].mcDCharge=mcDCharge;
		      pDDStarPi[currentCounts[iDDStarPi]].dType=dType[iM];
		      pDDStarPi[currentCounts[iDDStarPi]].mcIsDStar=mcIsDStar;
		      pDDStarPi[currentCounts[iDDStarPi]].leptonId=leptonId[iM];
		      pDDStarPi[currentCounts[iDDStarPi]].eWeight=eWeight;
		      fillWeightDetails(pDDStarPi[currentCounts[iDDStarPi]],i);
		      pDDStarPi[currentCounts[iDDStarPi]++].fileIndex=i;
		      //		  cout <<"count for ddstar pi: "<< currentCounts[iDDStarPi] <<" " << iDDStarPi <<" event: " << j <<endl;
		      ///		  cout <<"is dd star pi" <<endl;
		    }


	      /// wrong channel and Cross-feed will be constructed when we build the 
	      // channel specific histograms since they are channel dependent
	      //DDStarPiWrongChannel
	      //DDStarPiCrossFeed
	      
	      //DDStarPiPi
		  if(isDDStarPiPi())
		    {
		      numClasses++;
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].mNu2=mNu2[iM];
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].weight=weight;
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].mcBCharge=mcBCharge;
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].bestBCharge=bestBCharge;
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].dCharge=dCharge[iM];
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].mcDCharge=mcDCharge;
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].dType=dType[iM];
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].mcIsDStar=mcIsDStar;
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].leptonId=leptonId[iM];
		      pDDStarPiPi[currentCounts[iDDStarPiPi]].eWeight=eWeight;
		      fillWeightDetails(pDDStarPiPi[currentCounts[iDDStarPiPi]],i);
		      pDDStarPiPi[currentCounts[iDDStarPiPi]++].fileIndex=i;
		  ///		  cout <<"is dd star pi pi" <<endl;
		    }
		  //DLNu
		  
		  
		  if(isDLNu())
		    {
		      numClasses++;
		      pDLNu[currentCounts[iDLNu]].mNu2=mNu2[iM];
		      pDLNu[currentCounts[iDLNu]].weight=weight;
		      pDLNu[currentCounts[iDLNu]].mcBCharge=mcBCharge;
		      pDLNu[currentCounts[iDLNu]].bestBCharge=bestBCharge;
		      pDLNu[currentCounts[iDLNu]].dCharge=dCharge[iM];
		      pDLNu[currentCounts[iDLNu]].mcDCharge=mcDCharge;
		      pDLNu[currentCounts[iDLNu]].dType=dType[iM];
		      pDLNu[currentCounts[iDLNu]].mcIsDStar=mcIsDStar;
		      pDLNu[currentCounts[iDLNu]].leptonId=leptonId[iM];
		      pDLNu[currentCounts[iDLNu]].eWeight=eWeight;
		      fillWeightDetails(pDLNu[currentCounts[iDLNu]],i);
		      pDLNu[currentCounts[iDLNu]++].fileIndex=i;
		      ///		  cout <<"is dlnu" <<endl;
		    }
	      //DPiPiLNu
		  if(isDPiPiLNu())
		    {
		      numClasses++;
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].mNu2=mNu2[iM];
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].weight=weight;
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].mcBCharge=mcBCharge;
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].bestBCharge=bestBCharge;
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].dCharge=dCharge[iM];
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].mcDCharge=mcDCharge;
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].dType=dType[iM];
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].mcIsDStar=mcIsDStar;
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].leptonId=leptonId[iM];
		      pDPiPiLNu[currentCounts[iDPiPiLNu]].eWeight=eWeight;
		      fillWeightDetails(pDPiPiLNu[currentCounts[iDPiPiLNu]],i);
		      pDPiPiLNu[currentCounts[iDPiPiLNu]++].fileIndex=i;
		      ///		  cout <<"is d pi pi " <<endl;
		    }
		  
		  
		  //DStarLNu
		  if(isDStarLNu())
		    {
		      numClasses++;
		      pDStarLNu[currentCounts[iDStarLNu]].mNu2=mNu2[iM];
		      pDStarLNu[currentCounts[iDStarLNu]].weight=weight;
		      pDStarLNu[currentCounts[iDStarLNu]].mcBCharge=mcBCharge;
		      pDStarLNu[currentCounts[iDStarLNu]].bestBCharge=bestBCharge;
		      pDStarLNu[currentCounts[iDStarLNu]].dCharge=dCharge[iM];
		      pDStarLNu[currentCounts[iDStarLNu]].mcDCharge=mcDCharge;
		      pDStarLNu[currentCounts[iDStarLNu]].dType=dType[iM];
		      pDStarLNu[currentCounts[iDStarLNu]].mcIsDStar=mcIsDStar;
		      pDStarLNu[currentCounts[iDStarLNu]].leptonId=leptonId[iM];
		      pDStarLNu[currentCounts[iDStarLNu]].eWeight=eWeight;
		      fillWeightDetails( pDStarLNu[currentCounts[iDStarLNu]],i);
		      pDStarLNu[currentCounts[iDStarLNu]++].fileIndex=i;
		      ///		  cout <<"is d star l nu" <<endl;
		    }
		  //DStarPiPiLNu
		  
		  
		  if(isDStarPiPiLNu())
		    {
		      numClasses++;
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].mNu2=mNu2[iM];
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].weight=weight;
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].mcBCharge=mcBCharge;
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].bestBCharge=bestBCharge;
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].dCharge=dCharge[iM];
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].mcDCharge=mcDCharge;
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].dType=dType[iM];
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].mcIsDStar=mcIsDStar;
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].leptonId=leptonId[iM];
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]].eWeight=eWeight;
		      fillWeightDetails( pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]],i);
		      pDStarPiPiLNu[currentCounts[iDStarPiPiLNu]++].fileIndex=i;
		      ///		  cout <<"is d star pipi l nu" <<endl;
		    }


	      //OtherBB
		  if(isOtherBB())
		    {
		      numClasses++;
		      pOtherBB[currentCounts[iOtherBB]].mNu2=mNu2[iM];
		      pOtherBB[currentCounts[iOtherBB]].weight=weight;
		      pOtherBB[currentCounts[iOtherBB]].mcBCharge=mcBCharge;
		      pOtherBB[currentCounts[iOtherBB]].bestBCharge=bestBCharge;
		      pOtherBB[currentCounts[iOtherBB]].dCharge=dCharge[iM];
		      pOtherBB[currentCounts[iOtherBB]].mcDCharge=mcDCharge;
		      pOtherBB[currentCounts[iOtherBB]].dType=dType[iM];
		      pOtherBB[currentCounts[iOtherBB]].mcIsDStar=mcIsDStar;
		      pOtherBB[currentCounts[iOtherBB]].leptonId=leptonId[iM];
		      pOtherBB[currentCounts[iOtherBB]].eWeight=eWeight;
		      fillWeightDetails( pOtherBB[currentCounts[iOtherBB]],i);
		      pOtherBB[currentCounts[iOtherBB]++].fileIndex=i;
		      ///		  cout <<"other" <<endl;
		    }
		  
		  
		}
	      if(numClasses>1)
		cout<<"event " << j <<"  seems to be in more than one class ("<<numClasses<<") "  <<endl;
	      
	      
	    }
	}
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

float sysDataGen::getWeightWithError(dataElement& element, bool foundDDoubleStar)
{
 double weight=1.0;
 float relTagCorrUncert=0.045;
 if(element.iTree==0 || element.iTree==1)
   {
     if(element.iTree==0)
       {
	 //mixed...
	 relTagCorrUncert=0.042;
       }
     int tagIndex=element.tagId*2+element.iTree;
     if(tagIndex<4)
       weight*=(element.tagCorr+relTagCorrUncert*tagIdWeight[tagIndex]);
   }
 else
   {
     weight*=element.tagCorr;
   }
  for(int i=0;i<element.numFFDCorrection;i++)
    {
      weight*=(element.ffDCorrection[i]+element.sFfDCorrection[i]*D_FFWeights[element.DFFCorrCombBin[i]]);
    }
  for(int i=0;i<element.numFFDdsCorrection;i++)
    {
      weight*=(element.ffDdsCorrection[i]+element.sFfDdsCorrection[i]*Dds_FFWeights[element.DdsFFCorrCombBin[i]]);
    }

  for(int i=0;i<element.numKsCorrection;i++)
    {
      weight*=(element.KsCorrection[i]+element.sKsCorrection[i]*KsWeights[element.KsCorrCombBin[i]]);
    }
  for(int i=0;i<element.numChargedCorrection;i++)
    {
      weight*=(element.ChargedCorrection[i]+element.sChargedCorrection[i]*ChargedWeights[element.ChargedCorrCombBin[i]]);
    }
  for(int i=0;i<element.numPi0Correction;i++)
    {
      weight*=(element.pi0Correction[i]+element.sPi0Correction[i]*pi0Weights[element.pi0CorrCombBin[i]]);
    }
  for(int i=0;i<element.numD_DecayCorr;i++)
    {
      weight*=(element.D_decayCorr[i]+element.sD_decayCorr[i]*D_decTypeWeights[element.D_decTypeCorrBin[i]]);
    }
  if(element.B_decayType<12&&element.sB_DecayCorr>0)
    weight*=(element.B_DecayCorr+element.sB_DecayCorr*B_decWeights[element.B_decayType]);
  //1.4 % relative for each experiment
  weight*=(CrossSectionLumiCorrection+0.014*CrossSectionLumiCorrection*lumiWeight[element.experiment]);
  
  weight*=trackWeights[element.numTracks];


  if(foundDDoubleStar)
    {
      if(element.iTree==0)
	{
	  weight*=mixedFactor;
	}
      if(element.iTree==1)
	{
	  weight*=chargedFactor;
	}

    }
  return weight;
}

float sysDataGen::getWeight(int iTree,bool foundDDoubleStar)
{

  double weight=1.0;
  weight*=(tagCorr);
  //  cout <<"weight after tagCorr: " << weight <<endl;
  //  cout <<"applying " << FFDCorrectionCounter <<" ff weights " <<endl;
  for(int i=0;i<FFDCorrectionCounter;i++)
    {
      weight*=(FFDCorrection[i]);
    }
  //  cout <<"weight after FF correction: "<< weight<<endl;
  for(int i=0;i<FFDdsCorrectionCounter;i++)
    {
      weight*=(FFDdsCorrection[i]);
    }
  //  cout <<"weight after FFDds correction: "<< weight <<endl;
  for(int i=0;i<KsCorrectionCounter;i++)
    {
      weight*=KsCorrection[i];
    }
  //  cout <<"weight after Ks correction: "<< weight<<endl;
  for(int i=0;i<ChargedCorrectionCounter;i++)
    {
      weight*=ChargedCorrection[i];
    }
  //  cout <<"weight after Charged correction: "<< weight<<endl;
  for(int i=0;i<Pi0CorrectionCounter;i++)
    {
      weight*=Pi0Correction[i];
    }
  //  cout <<"weight after pi0 correction: "<< weight<<endl;
  for(int i=0;i<D_DecayCorrCounter;i++)
    {
      weight*=D_DecayCorr[i];
    }
  //  cout <<"weight after D_decay correction: "<< weight<<endl;
  weight*=(B_DecayCorr);
  //  cout <<"weight after bdecay correction: "<< weight<<endl;
  weight*=(CrossSectionLumiCorrection);
  //  cout <<"weight after xsection/lumi correction: "<< weight<<endl;
  if(foundDDoubleStar)
    {
      if(iTree==0)
	{
	  weight*=mixedFactor;
	}
      if(iTree==1)
	{
	  weight*=chargedFactor;
	}

    }

  //  cout <<"returning weiht: "<< weight <<endl;
  return weight;
}

void sysDataGen::fillWeightDetails(dataElement& element, int iTree)
  {
    element.tagCorr=tagCorr;
    if(tagId<514)
      element.tagId=0;
    else
      element.tagId=1;
    element.iTree=iTree;
    element.B_DecayCorr=B_DecayCorr;
    element.sB_DecayCorr=sB_DecayCorr;
    element.CrossSectionLumiCorrection=CrossSectionLumiCorrection;
    element.experiment=experiment;
    element.numTracks=numTracks;
    if(sB_DecayCorr>0)
      element.B_decayType=B_decayType;


    element.numKsCorrection=KsCorrectionCounter;
    if(element.numKsCorrection>maxKs)
      element.numKsCorrection=maxKs;
    for(int i=0;i<element.numKsCorrection;i++)
      {
	element.KsCorrection[i]=KsCorrection[i];
	element.sKsCorrection[i]=sKsCorrection[i];
	element.KsCorrCombBin[i]=4*KsCorrMomBin[i]+KsCorrThetaBin[i];
      }
    element.numChargedCorrection=ChargedCorrectionCounter;
    if(element.numChargedCorrection>maxCharged)
      element.numChargedCorrection=maxCharged;
    for(int i=0;i<element.numChargedCorrection;i++)
      {
	element.ChargedCorrection[i]=ChargedCorrection[i];
	element.sChargedCorrection[i]=sChargedCorrection[i];
	element.ChargedCorrCombBin[i]=(9*30*2)*ChargedCorrThetaBin[i]+(9*2)*ChargedCorrMomBin[i]+9*ChargedCorrSVDBin[i]+ChargedCorrMisIdType[i];
      }


    element.numPi0Correction=Pi0CorrectionCounter;
    if(element.numPi0Correction>maxPi0)
      element.numPi0Correction=maxPi0;
    for(int i=0;i<element.numPi0Correction;i++)
      {
	element.pi0Correction[i]=Pi0Correction[i];
	element.sPi0Correction[i]=sPi0Correction[i];
	element.pi0CorrCombBin[i]=pi0MomBin[i];
      }



    element.numD_DecayCorr=D_DecayCorrCounter;
    if(element.numD_DecayCorr>maxDDec)
      element.numD_DecayCorr=maxDDec;
    for(int i=0;i<element.numD_DecayCorr;i++)
      {
	element.D_decayCorr[i]=D_DecayCorr[i];
	element.sD_decayCorr[i]=sD_DecayCorr[i];
	element.D_decTypeCorrBin[i]=D_decType[i];
      }

    element.numFFDCorrection=FFDCorrectionCounter;
    if(element.numFFDCorrection>maxFFD)
      element.numFFDCorrection=maxFFD;
    for(int i=0;i<element.numFFDCorrection;i++)
      {
	element.ffDCorrection[i]=FFDCorrection[i];
	element.sFfDCorrection[i]=sFFDCorrection[i];
	element.DFFCorrCombBin[i]=(2*24)*D_q2BinFF[i]+2*D_pBinFF[i]+D_TypeFF[i];
      }

    element.numFFDdsCorrection=FFDdsCorrectionCounter;
    if(element.numFFDdsCorrection>maxFFDds)
      element.numFFDdsCorrection=maxFFDds;
    for(int i=0;i<element.numFFDdsCorrection;i++)
      {
	element.ffDdsCorrection[i]=FFDdsCorrection[i];
	element.sFfDdsCorrection[i]=sFFDdsCorrection[i];
	element.DdsFFCorrCombBin[i]=(10*4)*Dds_wBinFF[i]+4*Dds_cosTBinFF[i]+Dds_TypeFF[i];
      }

  }
