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
#include <iostream>

#include "bToDDoubleStar/PIDCorrections.h"

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

  //should get the tree data info
  float getBRCorrection();
  float calcBRCorrection();
  float getFFCorrection();

  // begin_run function
  void begin_run ( BelleEvent*, int* );
  void findDStar(vector<Hep3Vector>& allPB, vector<int>& allPB_Class, vector<int>& allPB_Charge);
  bool findDecaySignatureForBBRCorrection(const Gen_hepevt &mother,int& numLeptons, int& numPions, int& numKaons, int& numPi0, int& numBaryons,int& numNu,int* br_sigs);
  bool findDecaySignature(const Gen_hepevt &mother,bool& dDoubleStar,int& numLeptons, int& numPions, int& numKaons, int& numPi0, int& numBaryons,int& numD, int& numDStar,int& numNu, bool& dStar_2S, bool& d_2S);
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


  double bestBPx;
  double bestBPy;
  double bestBPz;

  set<int> chargedIds;
    set<int> pi0Ids;
    set<int> gammaIds;

    int getOverlap(set<int>& s1, set<int>& s2);

  bool recursivePrint(const Gen_hepevt gen_it, string s);
  bool isAnyD(int lund);
  bool getDDecayProducts(const Gen_hepevt, int& Kp, int& Km, int& Ks, int& Pip, int& Pim, int& Pi0, int& other);
  void computeD_BR_CorrectionFactor(double& corrFact,int Kp, int Km, int Ks, int Pip,int Pim, int Pi0,int other);

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

  PIDCorrections pidCorrections;

  int evtNr;
  int runNr;
  int expNr; 

  bool sig_FoundDDoubleStar;
  bool sig_FoundD;
  int sig_numPions;
  int sig_numKaons;
  int sig_numLeptons;
  int sig_numPi0;
  int sig_numBaryons;
  int sig_numD;
  int sig_numDStar;


  ///for corrections
  int br_sig_D0LNu;
  int br_sig_DLNu;

  int br_sig_DStarLNu;
  int br_sig_DStar0LNu;
    
  int br_sig_D1LNu;
  int br_sig_D2LNu;

  int br_sig_D1PrimeLNu;
  int br_sig_D0StarLNu;
  
  int br_sig_D10LNu;
  int br_sig_D20LNu;

  int br_sig_D1Prime0LNu;
  int br_sig_D0Star0LNu;
  
static const int br_sig_D0=0;
static const  int br_sig_D=1;

static const  int br_sig_DStar=2;
static const  int br_sig_DStar0=3;

static const  int br_sig_D1=4;
static const  int br_sig_D2=5;

static const  int br_sig_D1Prime=6;
static const  int br_sig_D0Star=7;

static const  int br_sig_D10=8;
static const  int br_sig_D20=9;
  
static const  int br_sig_D1Prime0=10;
static const  int br_sig_D0Star0=11;


    //

  float overlapFractionCharged;
  float overlapFractionPi0;

  int sigDPiLNu;
  int sigDPiPiLNu;
  int sigDLNu;

  int sigDStarPiLNu;
  int sigDStarPiPiLNu;
  int sigDStarLNu;


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
    float mcFactors[12];
    float dataFactors[12];
    float dataFactorError[12];
    float overallCorrFactors[12];


    float dDecayFactorsMC[25];
    float dDecayFactorsData[25];
    float dDecayFactorsDataErrors[25];

    float  D_BR_CorrectionFactor;
    float B_BR_CorrectionFactor;



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
  void recCheck(const Gen_hepevt& gen_it, vector<int>& ids, vector<int>& pids, int& numNu, bool onlyCorrectDDecays=true);

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
