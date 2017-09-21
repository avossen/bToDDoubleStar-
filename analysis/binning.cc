
//the external variables
int numBins[5];
double upperCut[5];
double lowerCut[5];

//in the main analysis code we assign them ourselves since we have to check for the partial box
//this is for the partial box stuff
void loadBinning()
{
  upperCut[0]=2.0;
  upperCut[1]=2.0;
  upperCut[2]=0.6;
  upperCut[3]=2.0;
  upperCut[4]=0.6;
  lowerCut[0]=-0.5;
  lowerCut[1]=-0.3;
  lowerCut[2]=-0.3;
  lowerCut[3]=-0.3;
  lowerCut[4]=-0.3;
  numBins[0]=140;
  numBins[1]=140;
  numBins[2]=70;
  numBins[3]=140;
  numBins[4]=70;

}
