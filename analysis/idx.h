#ifndef IDX__H__
#define IDX__H__

#include "doAnalysisCombined.h"

#define iNoSelection 0
#define iContinuum 0
#define iDDStar 1
#define iDDStarPi 2
#define iDDStarPiWrongChannel 3
#define iDDStarPiPi 4
#define iDLNu 5
#define iDPiPiLNu 6
#define iDStarLNu 7
#define iDStarPiPiLNu 8
#define iDDStarPiCrossFeed 9
#define iAll 10
#define iOtherBB 10

using namespace std;
void addCorrections(char* buffer);
//void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel=false);
void getFeedDown(char* feedDownSelection, int channel, bool dataAndMC, int numPions);
void insertSignalMCWeighting(char* buffer,  int fileIndex);
void copyHisto(TH1* first, TH1* out, int numBins, int maxBins,int shift=0)
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

//to be used in the channel combination
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

//get the template and apply possible mergers if necessary
//in principle there is a twist for the combined (so D and D*) channels. But since the D* is pretty clean, we probably get away w/o combination there 
//because all is held constant there anyways.
//since we do the simultaneous fit, everything should be more stable anyways...
void getTemplates(TH1F** summedComponents_in, TH1F** &templates, char** templateLegendNames, char** allLegendNames, int numComponents, int& numMergers, int numPions,int oneIdx=-1,bool combineDPiPi=true)
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

//note that 'maxdatatreesize' dictates how much data we fetch (used for the partial_box)
//get the data histogram by either summing over the four MC tree or (in the real data case) getting it from the last (data) tree
//
void getData(TH1F* &data, bool dataTree, TTree** trees, int channel, int numPions, int leptonId)
{
  char histoName[2009];
  char channelString[500];
  char drawCommand[2000];
  getChannelString(channel,channelString);
  int channelIdx=channel+1;
  //  TH1F *data;                              //data histogram
  //use first 4 trees (the last one is data)
  int treeCount=4;
  char corrBuffer[2000];
  char buffer[3000];
  addCorrections(corrBuffer);
  if(dataTree)
    {
	treeCount=1;
    }
  cout <<"tree count: " << treeCount <<endl;
  
  char channelSelectionData[1000];
  //  char feedDownSelection[1000];
  #ifdef Huschle_Signal_Weighting
  float mixedFactor=1.0/0.852255195;
  float chargedFactor=1.0/0.9266338302;
#else
  float mixedFactor=1.0;
  float chargedFactor=1.0;
#endif
  //the lumi correction should only be applied for the MC_TEST case. Otherwise we should just add the templates so we get the correct weighting (which is dependent on the trees)
  //for now we stay with the separate data addition to have a cross-check for the template method (do we get the same number of counts as if we add up the templates?)
  float huschleLumiFactor=mixedFactor;
  if(channel==-1)
    {
      huschleLumiFactor=(mixedFactor+chargedFactor)/2;
    }
  if(channel>1)
    huschleLumiFactor=chargedFactor;
  huschleLumiFactor-=1;


  char huschleMC_lumi_corr[1000];
  sprintf(huschleMC_lumi_corr,"(1+foundAnyDDoubleStar*%f)*",huschleLumiFactor);



  // is that actually used?
    //    char channelSelectionDataAndMC[1000];
    //void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions)
  getChargeAndStarSelection(channelSelectionData,channel,false,numPions);
  //  getFeedDown(feedDownSelection,channel,false,numPions);
  //    getChargeAndStarSelection(channelSelectionDataAndMC,channel,true,numPions);
  
  if(leptonId!=0)
    {
      //for the MC_test we still have to do the tag corr....
      if(dataTree)
	      {
#ifdef MC_TEST
		if(numPions==0)
		  sprintf(buffer,"%s %s tagCorr*" P0STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",huschleMC_lumi_corr,corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		if(numPions==1)
		  sprintf(buffer,"%s %s tagCorr*" P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",huschleMC_lumi_corr,corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		if(numPions==2)
		  sprintf(buffer,"%s %s tagCorr*(" P2STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",huschleMC_lumi_corr,corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		
#else
		if(numPions==0)
		  {
		    sprintf(buffer,"%s " P0STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		  }
		if(numPions==1)
		  {
		    sprintf(buffer,"%s " P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		  }
		if(numPions==2)
		  {
		  sprintf(buffer,"%s (" P2STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		  }
#endif
	      }
	    else
	      {
		if(numPions==0)
		  {
		    sprintf(buffer,"%s tagCorr*" P0STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
		  }
		if(numPions==1)
		  {
		    sprintf(buffer,"%s tagCorr*" P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId, channelSelectionData);
				  //		  sprintf(buffer,"%s " P1STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut,lowerCut,numPions,leptonId, channelSelectionData);
		  }
		if(numPions==2)
		  {
		    sprintf(buffer,"%s tagCorr*(" P2STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,leptonId,channelSelectionData);
				  //		  sprintf(buffer,"%s (" P2STRING " && bestBCharge==((-1)*systemCharge) && abs(leptonId)==%d  && %s)  ",corrBuffer,upperCut[channelIdx],lowerCut,numPions,leptonId,channelSelectionData);
		  }
	      }
		  //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0  && abs(leptonId)==%d)  ",numPions,leptonId);
    }
  else
    {
      if(dataTree)
	{
	  //for the MC_test we still have to do the tag corr....
#ifdef MC_TEST
		if(numPions==0)
		  {
		    sprintf(buffer,"%s %s tagCorr*" P0STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",huschleMC_lumi_corr,corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
				    //		    sprintf(buffer,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
		if(numPions==1)
		  {
		    sprintf(buffer,"%s %s tagCorr*" P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",huschleMC_lumi_corr,corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
				    //		    sprintf(buffer,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
		if(numPions==2)
		  {
		    sprintf(buffer,"%s %s tagCorr*(" P2STRING "  && bestBCharge==((-1)*systemCharge) && %s) ",huschleMC_lumi_corr,corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
#else
		if(numPions==0)
		  {
		    sprintf(buffer,"%s " P0STRING "  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
		if(numPions==1)
		  {
		    sprintf(buffer,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
		if(numPions==2)
		  {
		  sprintf(buffer,"%s (" P2STRING "  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
#endif
	}
	    else
	      {
		if(numPions==0)
		  {
		    sprintf(buffer,"%s tagCorr*" P0STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions, channelSelectionData);
		    //		    sprintf(buffer,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions, channelSelectionData);
		  }
		if(numPions==1)
		  {
		    sprintf(buffer,"%s tagCorr*" P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions, channelSelectionData);
		    //		    sprintf(buffer,"%s " P1STRING "  && bestBCharge==((-1)*systemCharge) && %s ) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions, channelSelectionData);
		  }
		if(numPions==2)
		  {
		    		    sprintf(buffer,"%s tagCorr*(" P2STRING "  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
				    //		    sprintf(buffer,"%s (" P2STRING "  && bestBCharge==((-1)*systemCharge) && %s) ",corrBuffer,upperCut[channelIdx],lowerCut[channelIdx],numPions,channelSelectionData);
		  }
	      }
	    //sprintf(buffer,"tagCorr*CrossSectionLumiCorrection*(mNu2<1.0 && mNu2>-1.0 && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.05 && deltaETag<0.05 && logProb > -3 && mDnPi < 3.0 )  ",numPions);
    }

  char bufferTmp[500];
    cout <<" treeCount: " << treeCount <<endl;
    for(int tc=0;tc<treeCount;tc++)
      {
	cout << " tc: " << tc<< " numPions: "<< numPions<<" leptId: "<< leptonId <<" channel: " << channelString<<endl;
	sprintf(histoName,"histo_Data_%d_pions_%d_leptonId_treeNr%d_%s",numPions,leptonId,tc, channelString);
	cout <<"draw command..." <<endl;
	TH1D* histo=new TH1D(histoName,histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	histo->Sumw2();
	sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
	//sprintf(drawCommand,"mNu2 >> %s",histoName);
	cout << " and the actual draw... "<< drawCommand<<"--> >" <<buffer<<endl;
	int counts=0;
	if(!dataTree)
	  {

	    //add the huschle factor for the mixed and charged trees
	    if(tc<2)
	      {
		//don't want to insert the fractions permanently
		sprintf(bufferTmp,"%s",buffer);
		insertSignalMCWeighting(bufferTmp,tc);
		cout <<"after inserting, using: " << bufferTmp <<endl;
	      }
	    counts=trees[tc]->Draw(drawCommand,(char*)bufferTmp);
	  }
	else
	  {
	    //grab last tree which should be the data tree...
	    counts=trees[4]->Draw(drawCommand,(char*)buffer,"",maxDataTreeSize);
	  }
	cout <<"got " << counts <<" counts from data selected " <<endl;
	if(tc==0)
	  data=(TH1F*)gDirectory->Get(histoName);
	else
	  data->Add((TH1F*)gDirectory->Get(histoName));

	data->Sumw2();
	data->SetFillColor(glColorTable[0]->GetNumber());	
      }
    cout <<"before return, integral: "<< data->Integral() <<endl;

}

#endif
