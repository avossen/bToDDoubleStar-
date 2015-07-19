#ifndef DIHADANA_H
#define DIHADANA_H

#include MDST_H
#include EVTCLS_H
#include "TFile.h"
#include "TH1F.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "particle/Particle.h"
#include "BTreeData.h"

//#include "bToDDoubleStar/HadronQuadruple.h"
//#include "bToDDoubleStar/HadronPair.h"
#include "bToDDoubleStar/EventInfo.h"
#include "bToDDoubleStar/DebugHistos.h"
#include "bToDDoubleStar/AuxFunc.h"
#include "bToDDoubleStar/TreeSaver.h"
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace fastjet;
using namespace std;



using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  enum DType{dtype_D0,dtype_DCharged,dtype_DStar,dtype_end};

  typedef std::vector<Gen_hepevt*> genhep_vec;
// Module class
class bToDDoubleStar : public Module 
{
public:  
  // constructor
  bToDDoubleStar( void );

  // destructor
  ~bToDDoubleStar( void );
  BTreeData treeData;

  // initialize
  void init ( int * );
  void printMCList();
  int goodHadronB() const;
  // begin_run function
  void begin_run ( BelleEvent*, int* );
  void findDStar(vector<Hep3Vector>& allPB, vector<int>& allPB_Class, vector<int>& allPB_Charge);
  bool findDecaySignature(const Gen_hepevt &mother,bool& dDoubleStar,int& numLeptons, int& numPions, int& numKaons, int& numPi0, int& numBaryons, bool& dStar_2S, bool& d_2S);
  void disp_stat ( const char* ){}
  void saveHistos( vector<Hep3Vector>& v_allParticlesBoosted, vector<Hep3Vector>& v_allParticlesNonBoosted);
  void saveTree();

  void cleanUp();
  void exitEvent();
  // histogram initialize
  void hist_def ( void );

  // event function
  void event ( BelleEvent*, int* );

  // end_run function
  void end_run ( BelleEvent*, int* );

  //  void other ( int*, BelleEvent*, int* ){}

  // terminate
  void term ( void );
  bool validRun;
  int test;
  int zNums[4];
  double smpl_;
  char rFileName[200];

  bool recursivePrint(const Gen_hepevt gen_it, string s);



static int getBin(vector<float>& b1, float value)
{
  int coo1=-1;
  for(int i=0;i<b1.size();i++)
    {
      if(value<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  return coo1;
}

  static float getPhi(const Hep3Vector& axis, const Hep3Vector& input)
    {
      return AuxFunc::getPhi(axis,input);
    }
  static float getPhi(const Hep3Vector& axis, const Particle& input)
    {
      return AuxFunc::getPhi(axis,input);
    }
  static float getTheta(const Hep3Vector& axis, const Particle& input)
    {
      return AuxFunc::getTheta(axis,input);
    }
  //reference particle types
  Ptype cPiPlus;
  Ptype cPiNeg;
  Ptype cPiZero;
  Ptype cKPlus;
  Ptype cKNeg;
 protected:

  bool sig_FoundDDoubleStar;
  bool sig_FoundD;
  int sig_numPions;
  int sig_numKaons;
  int sig_numLeptons;
  int sig_numPi0;
  int sig_numBaryons;
  bool sig_dStar_2S;
  bool sig_d_2S;

  bool found_2SD;
  bool found_2SD_Star;

    vector<Particle*> chargedPiCandidates;
    vector<Particle*> chargedKCandidates;
    vector<Particle*> pi0Candidates;
    vector<Particle*> KsCandidates;
    vector<Particle*> leptonCandidates;
    vector<Particle*> otherChargedTracks;
    vector<Particle*> D0Candidates;
    vector<Particle*> chargedDCandidates;
    vector<Particle*> DStarCandidates;
    bool m_mc;

    void reconstructD0();
    void reconstructChargedD();
    void reconstructDStar();
    unsigned doKmFit(Particle &p, double& conflevel, int debug, double mass=0);
    unsigned doKmVtxFit(Particle &p, double& conflevel, int debug);
    unsigned doKmVtxFit2(Particle &p, double& conflevel, int debug, double mass=0);



private:
  //compute distance between decay vertices of quark and antiquark

    TTree* mesonMassTree;
    float dMass;
    int dType;
    int dDecay;
    int foundDDoubleStarDecay;
    int allDTracksFound;


    bool foundSinglePionDecay;
  float getDecayDist();
  float getDecayDistK();
  float getDecayDistD();
  float getDecayDistPN();
  bool checkForDPiPi(int& bMesonId, bool& foundSinglePionDecay,bool print=false);
  bool checkDDoubleStarDecay(bool& foundSinglePionDecay,const  Gen_hepevt &m_mother, int daughterId,int& numGDaughters, bool& foundSameChargePions, bool& foundPiPlus, bool& foundPiMinus,  bool& foundD, bool& foundDStar,bool print=false);

  genhep_vec *getDaughters(const Gen_hepevt &mother);
  Hep3Vector getVertex(bool firstHemi);

  bool isD0(const Gen_hepevt& gen_it);
  bool isDStar(const Gen_hepevt& gen_it);
  bool isChargedD(const Gen_hepevt& gen_it);


  TH1D* histoD0CandidateMass;
  TH1D* histoKs;
  TH1D* histoNbOut;
  TH1D* histoPi0SlowMom;
  TH1D* histoB0Tag_M;
  TH1D* histoChargedBTag_M;
  TH1D* histoB0Tag_dE;
  TH1D* histoChargedBTag_dE;
  TH1D* histoD0Spect;
  TH1D* histoDStar;
  TH1D* histoPiSlowMom;

  TH1D* histoRecDSpect;
  TH1D* histoRecD0Spect;
  TH1D* histoRecDStarSpect;

  TH1D* histoRecDStarSpectToDPi0;
  TH1D* histoRecDStarSpectToDPi;
  TH1D* histoRecDStarSpectToD0Pi;
  TH1D* histoRecDStarSpectToD0Pi0;

  bool isNu(int lund)
  {
    if(abs(lund)==12 || abs(lund)==14 || abs(lund)==16)
      return true;
    return false;
  }

  bool detectable(int lund)
  {
    int ml=abs(lund);
    if(ml==211)
      return true;
    if(ml==321)
      return true;
    if(111==ml)
      return true;
    if(13==ml)
      return true;
    if(11==ml)
      return true;
    if(22==ml)
      return true;
    if(2212==ml)
      return true;

    return false;

  }

  bool getDecayIds(Particle& p,set<int>& chrId,set<int>& pi0Id, set<int>& gammaId, bool print=false);
  bool getMassFromDecayParticles(Particle& p, HepLorentzVector& mom);
  void recCheck(const Gen_hepevt& gen_it, vector<int>& ids, vector<int>& pids, int& numNu);

  float getTheta(Particle* p1,Particle* p2)
  {
	    double E1,E2;
	    E1=p1->e();
	    E2=p2->e();
	    Hep3Vector mH1=p1->p().vect();
	    Hep3Vector mH2=p2->p().vect();
	    Hep3Vector P_h=p1->p().vect()+p2->p().vect();
	    mH1.rotateZ(-P_h.phi());
	    mH1.rotateY(-P_h.theta());
	    float beta=(mH1.z()+mH2.z())/(E1+E2);
	    float gamma=1/sqrt(1-beta*beta);
	    mH1.setZ(-gamma*beta*E1+gamma*mH1.z());
	    float decayTheta=mH1.theta();
	    return decayTheta;
  }

  float getZ(Particle* p1,Particle* p2)
  {
	    double E1,E2;
	    E1=p1->e();
	    E2=p2->e();
	    float m_z=2*(E1+E2)/kinematics::Q;
	    return m_z;
  }

  //inits the tree and gives the treesaver the mva variable
  void initMvaTree();
  //tree to keep the vars needed for mva to discern charm/uds




  void getDrDz(Mdst_charged_Manager::iterator, int, double&, double&, double&, double&, double&, float);

  bool foundDPiPi;
  //D decay is in a mode we do not reconstruct
  bool noDRec;
  vector<double> v_vertexR;
  vector<double> v_vertexZ;

  //for histogramms of gamma energy coming from pi0 and not from pi0
  vector<float> v_pi0GammaE;
  vector<float> v_gammaE;
  vector<float> v_asyms;
  //event level info
  EventInfo m_evtInfo;
  TreeSaver* pTreeSaver;
  vector<float> dataF; //float data saved in the tree
  vector<int> dataI; //int data    "

  TFile* m_file;
  DebugHistos m_histos;

  //get phi after the given vector is the z direction

};


extern "C" Module_descr *mdcl_bToDDoubleStar()
{

  cout <<"register module...  " <<endl;
  bToDDoubleStar* module = new bToDDoubleStar;
  Module_descr* dscr = new Module_descr("bToDDoubleStar", module);
  dscr->define_param("smpl", "test parameter", &module->smpl_);
  //hopefully the int is the maximum lenght of the string...
  dscr->define_param("rfname","root file name","S",100,&module->rFileName);
  //registers parameters of ipprofile...
  IpProfile::define_global(dscr);
  return dscr;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
