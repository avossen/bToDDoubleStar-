#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include <math.h>
#include "particle/Particle.h"
#include "bToDDoubleStar/AnaDefs.h"
#include "bToDDoubleStar/AnaConsts.h"
//#include "fastjet/ClusterSequence.hh"


//hopefully the correct lund codes..
#define lc_piPlus 211
#define lc_kPlus 321
#define lc_pi0 111
#define lc_piMinus -211
#define lc_kMinus -321
#define lc_pPlus 2212
#define lc_pMinus -2212



using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class AuxFunc
  {
  public:
    static int sgn(float f)
      {
	if(f>=0)
	  return 1;
	return -1;
      }



  static float getPhi(const Hep3Vector& axis, const Hep3Vector& input)
    {
      //      cout <<"get phi from vec " <<endl;
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input;

      //      cout <<"input phi: " << m_input.phi()<<endl;
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());
      //           cout <<"input phi new :" << m_input.phi()<<endl;
            float norm1=kinematics::firstElectronCM.vect().cross(axis).mag();
            Hep3Vector zAxis(0,0,1);
            float norm1Alt2=zAxis.cross(axis).mag();
            float norm2=axis.cross(input).mag();
	    //	    if(0==norm1)
	    //	      cout <<"norm1 is zero!!!" <<endl;
	    //	    if(0==norm2)
	      {
		//	      cout <<"norm2 is zero!!!" <<endl;
		//	      cout << " axis: " << axis << " input: " << input <<endl;
	      }
	      //	    if(0==norm1Alt2)
	      //	      cout <<"norm1 Alt2 is zero!!" <<endl;

            float macosAlt2=acos((zAxis.cross(axis)).dot(axis.cross(input))/(norm1Alt2*norm2));
	    //            float macos=acos((kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input))/(norm1*norm2));
	    //            cout <<macos<< " norm1: " << norm1 <<"norm2; " << norm2 << endl;
	    //	    cout <<"macos: "<< macos <<" or do angle: " << (kinematics::firstElectronCM.vect().cross(axis)).angle(axis.cross(input)) <<endl;
	    float macos= (kinematics::firstElectronCM.vect().cross(axis)).angle(axis.cross(input));
            float altPhi2=sgn(axis.dot(zAxis.cross(axis).cross(axis.cross(input))))*macosAlt2;


            float altPhi=sgn(axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input))))*macos;

	    //             cout <<"cos phi alt: " << cos(altPhi) <<" cos new zAxis: " << cos(altPhi2) <<endl;
	    //cout <<"phi: " << m_input.phi() <<" alt phi: " << altPhi <<endl;
	    //	    cout << m_input.phi() <<endl;
	    //      return m_input.phi();
	    //     cout << " compare : " << altPhi <<"  swith " <<altPhi2<<" and " << m_input.phi() << endl;
     //altPhi2 and input.phi are exactly the same...
	    //    return m_input.phi();
      //      return altPhi2;
	    if(isnan(altPhi) )
	      {
		//		cout <<"alt phi is not a number... macos: " << macos <<" other fact: " << axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input)))<<endl;
		//		cout <<"input: "<< input<<" firstECM: " << kinematics::firstElectronCM.vect() <<" axis: " << axis <<endl;
		//		cout <<"is: " <<  kinematics::firstElectronCM.vect().cross(axis) <<" parallel to " << axis.cross(input) << " angle; "<<  kinematics::firstElectronCM.vect().cross(axis).angle(axis.cross(input))<<endl;
		//		cout <<"input phi: "<< m_input.phi()<<endl;
		//		cout <<"angle beam to R: "<< kinematics::firstElectronCM.vect().angle(input) <<endl;
	      }
	            return altPhi;
      //      return altPhi2;
    }
  static float getPhi(const Hep3Vector& axis, const Particle& input)
    {
      //      cout <<"get phi from part " <<endl;
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input.p().vect();
      //      cout <<"inputP phi: " << m_input.phi()<<endl;
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());
      //      cout <<"inputP phi new :" << m_input.phi()<<endl;
      //alternative berechnung...  //second electron is e+ (I think)
      float norm1=kinematics::firstElectronCM.vect().cross(axis).mag();
      float norm2=axis.cross(input.p().vect()).mag();
      //float macos=acos((kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect()))/(norm1*norm2));
      float macos=(kinematics::firstElectronCM.vect().cross(axis)).angle(axis.cross(input.p().vect()));
      float altPhi=sgn(axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input.p().vect()))))*macos;
      //         cout <<"old phi: " << m_input.phi() <<" alt phi: " << altPhi <<endl;  //kommt dasselbe raus!!
	    //	    cout <<" altPhi: "<< altPhi <<" other: "<< m_input.phi()<<endl;
      //return m_input.phi();
      if(isnan(altPhi))
	{
	 cout <<"aux::getPhi is not number nan, axis: " << m_axis <<" input: " << m_input << " norm1: " << norm1 <<" norm2: " << norm2 <<" altPHi: " << altPhi <<endl;
	 cout <<"input phi: "<< m_input.phi()<<endl;
	 cout <<"electron vec: " << kinematics::firstElectronCM.vect() <<", first fact "<< axis.dot(kinematics::firstElectronCM.vect().cross(axis).cross(axis.cross(input.p().vect())));
	 cout <<" second factor: " << acos((kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect()))/(norm1*norm2)) <<endl;
	 cout << " second fact first part: " <<kinematics::firstElectronCM.vect().cross(axis) <<" , " << axis.cross(input.p().vect()) <<" dot: " << (kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect())) <<" n*n2: " << norm1*norm2 <<endl;
	 cout <<" is "  <<kinematics::firstElectronCM.vect().cross(axis) <<" parallel to " << axis.cross(input.p().vect())<< kinematics::firstElectronCM.vect().cross(axis).angle(axis.cross(input.p().vect()))<<endl;
	 cout <<"input to acos: " << (kinematics::firstElectronCM.vect().cross(axis)).dot(axis.cross(input.p().vect()))/(norm1*norm2)<<endl;

	}
      return altPhi;
    }
  static float getTheta(const Hep3Vector& axis, const Particle& input)
    {
      //easier to use scaler product and asin...
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input.p().vect();
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());
      return m_input.theta();
    }
  static float getTheta(const Hep3Vector& axis, const Hep3Vector& input)
    {
      //easier to use scaler product and asin...
      Hep3Vector m_axis=axis;
      Hep3Vector m_input=input;
      m_input.rotateZ(-m_axis.phi());
      m_input.rotateY(-m_axis.theta());
      return m_input.theta();
    }
  static AnaDef::TwoHadCharge getCharge(const Ptype& pt1, const Ptype pt2)
    {
      if(pt1.charge() > 0)
	{
	  if(pt2.charge() > 0)
	    return AnaDef::PP;
	  if(pt2.charge() < 0)
	    {
	    return AnaDef::PN;
	    }

	  return AnaDef::PZ;
	}
      else
	{
	  if(pt1.charge()<0)
	    {
	      if(pt2.charge() > 0)
		{
		  //		  cout <<"forbidden charge combination" << endl;
		  //not relly, in mc possible (of course not wanted)
		  return AnaDef::NP;
		}

	      if(pt2.charge() < 0)
		return AnaDef::NN;


	      cout <<"forbidden charge combination" << endl;
	      return AnaDef::NZ;
	      //	      return AnaDef::NZ;
	    }
	      if(pt2.charge() > 0)
		{
		  cout <<"forbidden charge combination" << endl;
		  return AnaDef::NP;
		}
	      if(pt2.charge() < 0)
		return AnaDef::ZN;

	      return AnaDef::ZZ;
	}
      cout <<"no matching charge!" << endl;
      exit(0);

    }

static  AnaDef::TwoHadPType getPType(const Ptype& pt1,const Ptype& pt2)
    {
      //      cout <<"(getPType)first: " << pt1.name() << " second "<< pt2.name()<<endl;
      //      cout <<"bool 1: " << (cPiZero==pt1) <<", bool 2 : " << (cPiZero==pt1) <<endl;
      if((pt1.lund()==lc_piPlus||pt1.lund()==lc_pi0) && (pt2.lund()==lc_piMinus||pt2.lund()==lc_pi0))
	{
	  //	  cout <<"return right" <<endl;
	  return AnaDef::PiPi;
	}
      if((pt1.lund()==lc_piPlus||pt1.lund()==lc_pi0) && pt2.lund()==lc_kMinus)
	{
	  return AnaDef::PiK;
	}
      if(pt1.lund()==lc_kPlus && (pt2.lund()==lc_piMinus || pt2.lund()==lc_pi0))
	{
	  return AnaDef::KPi;
	}
      if(pt1.lund()==lc_kPlus && pt2.lund()==lc_kMinus)
	{
	  return AnaDef::KK;
	}
      //and pt1 and pt2 switched (guess only relevant where one does not order according to charge, e.g.. leading analysis
      if((pt2.lund()==lc_piPlus||pt2.lund()==lc_pi0) && pt1.lund()==lc_piMinus)
	{
	  return AnaDef::PiPi;
	}
      if((pt2.lund()==lc_piPlus||pt2.lund()==lc_pi0) && pt1.lund()==lc_kMinus)
	{
	  return AnaDef::KPi;
	}
      if(pt2.lund()==lc_kPlus && pt1.lund()==lc_piMinus)
	{
	  return AnaDef::PiK;
	}
      if(pt2.lund()==lc_kPlus && pt1.lund()==lc_kMinus)
	{
	  return AnaDef::KK;
	}

      if(pt1.lund()==lc_pi0 && pt2.lund()==lc_pi0)
	{
	  return AnaDef::PiPi;
	}

      //add also for equal charge comginations
      if((pt1.lund()==lc_pi0 || abs(pt1.lund())==lc_piPlus) && (pt2.lund()==lc_pi0 ||abs(pt2.lund())==lc_piPlus))
	return AnaDef::PiPi;
      if((pt1.lund()==lc_pi0 || abs(pt1.lund())==lc_piPlus) && (abs(pt2.lund())==lc_kPlus))
	return AnaDef::PiK;
      if((pt2.lund()==lc_pi0 || abs(pt2.lund())==lc_piPlus) && (abs(pt1.lund())==lc_kPlus))
	return AnaDef::KPi;
      if((abs(pt1.lund())==lc_kPlus)&& (abs(pt2.lund())==lc_kPlus))
	return AnaDef::KK;
      //if nothing matches
      //      cout <<"getPType: nothing matches... lund1: " <<pt1.lund() << " lund2: " <<pt2.lund() <<endl;;
      //      exit(0);//some protons...
      return AnaDef::UNKNOWN;
    }



 static void getMaxLHMassIndexAlt(float e_id, float mu_id, float atcKPi, float atcPK, float& m_mass, int& mass_index)
   {
     bool isLepton=false;
      if(e_id>0.8&& mu_id<0.9)
	{
	  m_mass=m_e;
	  mass_index=0;
	  isLepton=true;
	}
      if(mu_id>0.9 && e_id<0.8)
	{
	  m_mass=m_muon;
	  mass_index=1;
	  isLepton=true;
	}
      if(mu_id>0.9&& e_id>0.8)
	{
	  m_mass=m_e;
	  isLepton=true;
	}

      if(!isLepton)
	{
	  if(atcKPi>0.6) //kaon
	    {
	      m_mass=m_k;
	      mass_index=3;

	    }
	  else
	    {
	      if(atcPK>0.8)  //in martins arbeit L_Kpr< 0.2 ->proton
		{
		  m_mass=m_pr;
		  mass_index=4;
		}
	      else  //pion
		{
		  if(atcKPi<0.3)
		    {
		      mass_index=2;
		      m_mass=m_pi;
		    }
		}
	    }
	}
   }





 static void getMaxLHMassIndexAmi(float e_id, float mu_id, float atcKPi, float atcPK, float& m_mass, int& mass_index)
   {
     bool isLepton=false;
      if(e_id>0.8&& mu_id<0.9)
	{
	  m_mass=m_e;
	  mass_index=0;
	  isLepton=true;
	}
      if(mu_id>0.9 && e_id<0.8)
	{
	  m_mass=m_muon;
	  mass_index=1;
	  isLepton=true;
	}
      if(mu_id>0.9&& e_id>0.8)
	{
	  m_mass=m_e;
	  isLepton=true;
	}

      if(!isLepton)
	{
	  if(atcKPi>0.6) //kaon
	    {
	      m_mass=m_k;
	      mass_index=3;

	    }
	  else
	    {
	      if(atcPK>0.8)  //in martins arbeit L_Kpr< 0.2 ->proton
		{
		  m_mass=m_pr;
		  mass_index=4;
		}
	      else  //pion
		{
		  if(atcKPi<0.3)
		    {
		      mass_index=2;
		      m_mass=m_pi;
		    }
		}
	    }
	}
   }



 static void getMaxLHMassIndex(float pid_e, float pid_mu,float pid_pi, float pid_K, float pid_p, float& m_mass, int& mass_index)
   {
      if(pid_pi > pid_e)
	{
	  if(pid_pi>pid_mu)
	    {
	      if(pid_pi>pid_K)
		{
		  if(pid_pi<pid_p)
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }//m_pi ist schon gesetzt
		}
	      else // K > pi > mu > e
		{
		  if(pid_K>pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	    }
	  else//mu > pi > e
	    { 
	      if(pid_mu > pid_K)
		{
		  if(pid_mu>pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	      else
		{
		  if(pid_K>pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	    }
	}
      else //e> pion
	{
	  if(pid_e>pid_mu)
	    {
	      if(pid_e>pid_K)
		{
		  if(pid_e>pid_p)
		    {
		      m_mass=m_e;
		      mass_index=0;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		} // K > e > mu > pi
	      else
		{
		  if(pid_K > pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }

		}
	    }
	  else// mu > e > pi
	    {
	      if(pid_mu>pid_K)
		{
		  if(pid_mu>pid_p)
		    {
		      m_mass=m_muon;
		      mass_index=1;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }
		}
	      else//K > mu > e >pi
		{
		  if(pid_K > pid_p)
		    {
		      m_mass=m_k;
		      mass_index=3;
		    }
		  else
		    {
		      m_mass=m_pr;
		      mass_index=4;
		    }

		}
	    }

	}
   }

  };








#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif
