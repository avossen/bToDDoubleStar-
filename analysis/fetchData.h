#ifndef FETCH_DATA
#define FETCH_DATA

//note that 'maxdatatreesize' dictates how much data we fetch (used for the partial_box)
//get the data histogram by either summing over the four MC tree or (in the real data case) getting it from the last (data) tree
//
void getData(TH1F* &data, bool dataTree, TTree** trees, int channel, int numPions, int leptonId)
{
  char histoName[2009];
  char xFdHistoName[2009];
  char channelString[500];
  char drawCommand[2000];
  getChannelString(channel,channelString);
  int channelIdx=channel+1;
  //  TH1F *data;                              //data histogram
  //use first 4 trees (the last one is data)
  int treeCount=4;
  char corrBuffer[2000];
  char buffer[3000];
  bool tmpFFCorrection=withFFCorrection;
  if(dataTree)
    {
	treeCount=1;
	//not all my data trees have this field (which is  ==1 anyways for data)
	withFFCorrection=false;
    }
  addCorrections(corrBuffer);
  withFFCorrection=tmpFFCorrection;

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
  getChargeAndStarSelection(channelSelectionData,channel,false,numPions,false);
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
	    counts=trees[4]->Draw(drawCommand,(char*)buffer,"",maxDataTreeSize,totalTreeSize-maxDataTreeSize);
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
#ifdef GENERATE_SINGLE_STREAM
    data->Scale(0.2);
    data->GetSumw2()->Set(0);
    cout <<"generate single stream, scale and remove sumw2 return, integral: "<< data->Integral() <<endl;
#endif
    //need to get signal fraction from MC

#ifdef MC_TEST
    cout <<"getting signal fraction for mc_test  with string: "<< gl_signalSelection << endl;
    //the data tree
    sprintf(histoName,"signalHisto");
    sprintf(xFdHistoName,"xFdHisto");
    //    TH1D* sigHisto=new TH1D(histoName,histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
    sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",histoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
    int sigCounts=trees[4]->Draw(drawCommand,gl_signalSelection,"");
    sprintf(drawCommand,"mNu2 >> %s(%d,%f,%f)",xFdHistoName,numBins[channelIdx],lowerCut[channelIdx],upperCut[channelIdx]);
    int xFdCounts=trees[4]->Draw(drawCommand,gl_xFeedSelection,"");
    TH1F* sigHisto=(TH1F*)gDirectory->Get(histoName);
    TH1F* xFdHisto=(TH1F*)gDirectory->Get(xFdHistoName);
    cout <<"got " << sigCounts <<" signal Counts, integral: "<< sigHisto->Integral() << " signal fraction: "<< sigHisto->Integral()/data->Integral()<<endl;
    gl_signalFraction[channelIdx]=sigHisto->Integral()/data->Integral();
    gl_xFeedFraction[channelIdx]=xFdHisto->Integral()/data->Integral();
    gl_dataInt[channelIdx]=data->Integral();
    gl_signalInt[channelIdx]=sigHisto->Integral();
    gl_xFeedInt[channelIdx]=xFdHisto->Integral();
#endif
}
#endif
