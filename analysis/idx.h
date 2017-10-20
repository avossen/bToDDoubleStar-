#ifndef IDX__H__
#define IDX__H__

#include "idx2.h"

//#include "doAnalysisCombined.h"
///////--> use this, but for partial box might need less...int numBins[]={70,70,70,70,70}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...

//int numBins[]={40,20,30,20,20}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma..



extern int numBins[5];
extern TColor* glColorTable[30];
//int numBins[]={140,140,70,140,70}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...



//int numBins[]={30,30,40,40,40}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//has to be the same for the channels we plan to combine
////-->>int numBins[]={30,40,40,40,40}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
///int numBins[]={30,100,100,100,100}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={30,20,20,20,20}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
///--->> for the high binsint numBins[]={30,100,40,100,40}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
////-->>int numBins[]={30,200,200,100,40}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...



//int numBins[]={20,8,8,8,8}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...

extern double lowerCut[5];
extern double upperCut[5];
extern void getChannelString(int channel, char* buffer);
extern void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel);
extern void getFeedDown(char* feedDownSelection, int channel, bool dataAndMC, int numPions);
//extern void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel);
extern bool withFFCorrection;
//we get meaningful results with 50 bins between -1 and 2..
///->allChannels, BToD, B0ToDStar, B0ToD, B0ToDStar
//int numBins[]={50,30,30,30,30};
//int numBins[]={70,200,70,200,70};

//int numBins[]={200,2000,2000,200,200}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={50,50,5,200,5}; //--> this together with lower -0.5, upper 2.0 leads to a mean exactly within one sigma...
//int numBins[]={70,100,70,70,70}; 
///_--> quite good int numBins[]={70,80,70,70,70};
///bias for BToD goes away with range -0.5 to 2.0 and 80 bins

///--->B0ToD seems to be fine with range of -0.5 to 1.5 with 70 bins...
////--->optimalfloat lowerCut[]={-0.65,-0.6,-0.65,-0.55,-0.25};
//float lowerCut[]={-0.65,-0.6,-0.65,-0.55,-0.25};
///-->should be -0.5, test -0.3 for fit stability 
//////float lowerCut[]={-0.5,-0.5,-0.5,-0.5,-0.5};
//needs to be double, otherwise there are problems with add...
////----->Double_t lowerCut[]={-0.5,-0.3,-0.3,-0.3,-0.3};
//////---Double_t lowerCut[]={-0.5,-0.3,-0.3,-0.3,-0.3};
///////Double_t lowerCut[]={-0.3,-0.2,-0.4,-0.2,-0.4};
//Double_t lowerCut[]={-0.3,-0.0,-0.0,-0.0,-0.0};
//Double_t lowerCut[]={-0.1,-0.1,-0.1,-0.1,-0.1};
//Double_t lowerCut[]={-0.0,-0.0,-0.0,-0.0,-0.0};
////////float upperCut[]={2.0,2.0,1.5,1.5,1.5};

///for fit stability set the upper cut a bit lower for the DStar where the counts are lower
//--->Double_t upperCut[]={0.5,0.5,0.4,0.5,0.4};
/////---Double_t upperCut[]={2.0,2.0,0.6,2.0,0.6};
///---->Double_t upperCut[]={0.5,0.5,0.4,0.5,0.4};
//Double_t upperCut[]={0.2,0.2,0.2,0.2,0.2};
//Double_t upperCut[]={0.3,0.3,0.3,0.3,0.3};

//float lowerCut[]={-0.5,-0.5,-0.5,-0.5,-0.5};
//float upperCut[]={2.0,2.0,1.5,1.5,1.5};
//float upperCut[]={1.0,1.0,1.0,1.0,1.0};
//float upperCut[]={1.0,0.25,0.2,0.3,0.15};
//float lowerCut=-0.5;

//float lowerCut=-1.0;
//float lowerCut=-0.2;
//float lowerCut=-0.5;
////float lowerCut=-0.3;
//float upperCut=2.0;
//float upperCut=0.4;
//float upperCut=0.25;

//float upperCut=1.0;

//int numBins=40;
//float lowerCut=-0.3;
//float upperCut=0.5;

//int numBins=150;

////int numBins=20;
//int numBins=100;


using namespace std;
//void addCorrections(char* buffer);
//void getChargeAndStarSelection(char* chargeAndStarSelection,int channel,bool dataAndMC, int numPions, bool wrongChannel=false);
//void getFeedDown(char* feedDownSelection, int channel, bool dataAndMC, int numPions);
//void insertSignalMCWeighting(char* buffer,  int fileIndex);

//get the template and apply possible mergers if necessary
//in principle there is a twist for the combined (so D and D*) channels. But since the D* is pretty clean, we probably get away w/o combination there 
//because all is held constant there anyways.
//since we do the simultaneous fit, everything should be more stable anyways...


#endif
