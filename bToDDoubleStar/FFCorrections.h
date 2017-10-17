#ifndef FFCorrections__H
#define FCorrections__H
#include <iostream>
#include <vector>
#include <math.h>
#include <map>

using namespace std;

/*
class to do the FF corrections based on BN1335
 */
class FFCorrections
{

 public:

  struct ffRet
  {
    int pBin;
    int q2Bin;
    int dType;
  };
  struct ffRetDds
  {
    int cosTBin;
    int wBin;
    int ddsType;
  };
  


  FFCorrections()
    {

      for(int i=1;i<=12;i++)
	{
	  limitsQ2.push_back(float(i));
	}
      float limit=0.1;
      for(int i=0;i<24;i++)
	{
	  limitsL_p.push_back(limit);
	  limit+=0.1;
	}

      limit=1.030;
      for(int i=0;i<20;i++)
	{
	  limitsW.push_back(limit);
	  limit+=0.03;
	}

      limit=-0.8;
      for(int i=0;i<10;i++)
	{
	  limitsCosT.push_back(limit);
	  limit+=0.2;
	}
    }
  void getWeightD(float lP,float q2, int lType, float& weight, float& stat, float& sys,ffRet& mRet);
  void getWeightDDStar(float w,float cosT, int dssIdx, float& weight, float& stat, float& sys, ffRetDds& mRet);
  int getBin(vector<float>& b1, float value);

 protected:
  static float fWeightsDstar[3744];
  static  float fWeightsD[2592];
  static float fWeightsDs[1296];
  static float fWeightsDsstar[1872];
  static float fWeightsD1[3000];
  static  float fWeightsD2[3000];
  static float fWeightsD0star[3000];
  static float fWeightsD1prime[3000];
  

  vector<float> limitsQ2;
  vector<float> limitsL_p;
  vector<float> limitsW;
  vector<float> limitsCosT;


 

};

#endif
