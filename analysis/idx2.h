#ifndef IDX2__H__
#define IDX2__H__

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


//#define P0STRING "( recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.0 "
//Robnin's cuts for x-check
#define P0STRING "( recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.24  && logProb > -3.0 "
//#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.0 && mDnPi < 3.5  "
#define P1STRING "(bestD==1 && pi1Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb > -3.5 && mDnPi < 3.5  &&mDnPi > 2.05"
////


#define P2STRING_SinglePion "bestD==1 && pi2Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb  > -2.6  && mDnPi < 3.0 && !(dType==1 && dDecay==2) && !(dType==1 && dDecay==3) && !(dType==2) && !overlapEvent" 

#define P2STRING_Search "bestD==1 && pi2Mom > 0.24 && recDecaySignature &&mNu2<%f && mNu2>%f && numRecPions==%d && mBTag> 5.27 && deltaETag>-0.18 && deltaETag<0.18 && logProb  > -2.6  && mDnPi < 3.0 && ((bestBCharge==0 && systemCharge==0) || (bestBCharge== -leptonCharge) ) && mDnPi>1.0  && mDnPi < 3.0&&  hypDMass1 > 2.03 && hypDMass2 > 2.03 && U < 2.0 && U > -2.0"  


//3,5,7 are the D pi pi 


#endif
