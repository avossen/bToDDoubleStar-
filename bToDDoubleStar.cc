const bool PRINT=false;
const bool kickOut2S=false;
//const bool PRINT=true;
const bool SAVE_MESON_MASS_DISTRIBUTIONS=false;
//const bool SAVE_MESON_MASS_DISTRIBUTIONS=true;
#include "bToDDoubleStar/ParticleInfoDecay.h"
#include <iomanip>
#include "mdst/findKs.h"
#include "bToDDoubleStar/mc.h"  //one central place to put the define mc
#include "tables/ekpfullrecon_panther.h"
#include "tables/mctype.h"
#include "event/BelleEvent.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include <TROOT.h>
//#include "LorentzVector.h" //clhep
#include "belle.h"
#include <kid/atc_pid.h>
#include <eid/eid.h>
#include <mdst/Muid_mdst.h>
#include "math.h"
#include "TMath.h"
//for neutral particles:
#include "ip/IpProfile.h"
#include <mdst/findKs.h>
#include <mdst/findLambda.h>
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"

#include "hamlet/AnaBrecon.h"
#include BRECON_H
#include "benergy/BeamEnergy.h"
#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/khelix2xyz.h"
#include "kfitter/kfitterparticle.h"

#include "bToDDoubleStar/weighting.h"

//#include "exkvertexfitter... " for more than 1 to 2 decay

#include <iostream>
#define pi0Mass 0.1349766
#define PY_RHO_0 113
#define PY_GAMMA 22
#define PY_E 11
#define PY_Tau 15
#define PY_NuE 12
#define PY_NuMu 14
#define PY_NuTau 16
#define PY_ELECTRON 11
#define PY_MU 13
#define PY_Mu 13
#define PY_PI 211
#define PY_K 321
#define PY_Pi0 111
#define PY_KS0 310
#define PY_B0 511
#define PY_B 521


#define PY_D 411
#define PY_D0 421
#define PY_DStar 413
#define PY_DStar0 423
#define PY_D_1 10413
#define PY_DStar0Plus 10411
#define PY_D_10 10423
#define PY_DStar_0s 10431
#define PY_DStar_00 10421
#define PY_DStar_20 425
#define PY_DStar_2 415
//this also lets the particle class crash...
#define PY_DStar_2S 100423

using namespace std;

#include MDST_H
#include EVTCLS_H

//#define SAVE_HISTOS

#define D0Mass 1.865
#define D0Width 0.6
#define D0Lower  D0Mass-D0Width
#define D0Upper  D0Mass+D0Width

#include <cmath>
//for thrust etc.
//strange: including this file in front of toolbox: thrust is not found anymore... (maybe "using belle namespace" a problem??)

//#include <mdst/Evtcls_hadron_info.h>
#if defined(BELLE_NAMESPACE)

namespace Belle {

#endif


  Hep3Vector& retSelf(Hep3Vector& vec)
  {
    return vec;
  };
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

//#define DEBUG_EVENT 287880//please no output
#define DEBUG_EVENT 1//please no output
//#define DEBUG_EVENT 15859
#define DEBUG_EVENT2 -287880
#include "bToDDoubleStar/AnaConsts.h"
#include "bToDDoubleStar/bToDDoubleStar.h"
#include "bToDDoubleStar/ParticleInfo.h"
#include "bToDDoubleStar/ParticleInfoMass.h"
#include "bToDDoubleStar/DebugHistos.h"
#include "particle/Ptype.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <time.h>
#include <fstream>

//#define XCHECK
//#define W_NEUTRAL

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  using namespace std;
  // Constructor


  bToDDoubleStar::bToDDoubleStar():smpl_(12345.),cPiPlus("PI+"), cPiNeg("PI-"),cPiZero("PI0"),cKPlus("K+"),cKNeg("K-")
  {
    ////---initialaize the meson tree
    mesonMassTree=new TTree("MesonMassTree","MesonMassTree");
    mesonMassTree->Branch("dMass", &dMass, "dMass/F");
    mesonMassTree->Branch("dType", &dType, "dType/I");
    mesonMassTree->Branch("dDecay", &dDecay, "dDecay/I");
    mesonMassTree->Branch("isSignal",&dIsSignal,"isSignal/I");
    mesonMassTree->Branch("foundDDoubleStarDecay",&foundDDoubleStarDecay,"foundDDoubleStarDecay/I");
    mesonMassTree->Branch("allDTracksFound",&allDTracksFound,"allDTracksFound/I");
    ///end 
    calcBRCorrection();

    m_mc=false;
#ifdef MC
    m_mc=true;
#endif
    //    cout <<"constructor..." <<endl;
    strcpy(rFileName,"notInitialized.root");
    //set default to charm correction same as uds
    mc_type=1003;
    excludeSignal=0;
    test=0;
    for(int i=0;i<4;i++)
      {
	//	zVals[i]=0;
      }

    histoRecD0Spect=new TH1D("Recd0Spect","Recd0Spect",100,1.5,2.0);
    histoRecDSpect=new TH1D("RecdSpect","RecdSpect",100,1.5,2.0);
    histoRecDStarSpect=new TH1D("RecdStarSpect","RecdStarSpect",150,1.5,2.2);

    histoRecDStarSpectToDPi0=new TH1D("RecdStarSpectToDPi0","RecdStarSpectToDPi0",150,1.5,2.2);
    histoRecDStarSpectToDPi=new TH1D("RecdStarSpectToDPi","RecdStarSpectToDPi0",150,1.5,2.2);
    histoRecDStarSpectToD0Pi=new TH1D("RecdStarSpectD0Pi","RecdStarSpectToD0Pi",150,1.5,2.2);
    histoRecDStarSpectToD0Pi0=new TH1D("RecdStarSpectD0Pi0","RecdStarSpectD0Pi0",150,1.5,2.2);

    histoD0CandidateMass=new TH1D("m_d0","m_d0",100,1.5,2.0);
    histoKs=new TH1D("m_ks","m_ks",100,0,1.0);
    histoNbOut=new TH1D("nbout","nbout",100,-1,1);

    histoPi0SlowMom=new TH1D("slowPi0",",slowPi0",1000,0,2);
    histoD0Spect=new TH1D("d0spect","d0spect",1000,0,3.0);
    histoDStar=new TH1D("dStarspect","dStarspect",300,1.8,5.0);
    histoPiSlowMom=new TH1D("piSlow","piSlow",100,0,3.0);


    histoChargedBTag_M=new TH1D("chargeB_M","chargedB_M",100,5.2,5.3);
    histoB0Tag_M=new TH1D("B0_M","B0_M",100,5.2,5.3);

    histoChargedBTag_dE=new TH1D("chargeBDeltaE","chargedBDeltaE",100,-0.2,0.2);
    histoB0Tag_dE=new TH1D("B0DeltaE","B0DeltaE",100,-0.2,0.2);

  }  
  ofstream* pXCheck;
  // Destructor
  bToDDoubleStar::~bToDDoubleStar(void)
  {

  }
  // initilization
  void bToDDoubleStar::init(int *status)
  {
        cout <<"init... (bTodDDouleStar)" <<endl;
#ifdef XCHECK
    pXCheck=new ofstream("xcheck");
#endif
    gROOT->SetStyle("Plain");
    m_file=new TFile(rFileName,"recreate");

    const double eler(3.499218);//energies of l, h beam
    const double eher(7.998213);


    //gives -999
    //double eler=BeamEnergy::E_LER();
    //    cout <<" eher: "<< BeamEnergy::E_HER() <<endl;
    //	cout <<"eher: "<< eher <<endl;

    const double theta(0.022);
    validRun=true;

    kinematics::cm=HepLorentzVector(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(eher*sin(theta),0.,eher*cos(theta),eher);
    kinematics::secondElectronCM=HepLorentzVector(0.,0.,-eler,eler);
    //they are not boosted yet...
    kinematics::pBeam=kinematics::firstElectronCM+kinematics::secondElectronCM;
    HepLorentzVector CMBoost2=kinematics::cm;
    //    CMBoost2.boost(-kinematics::CMBoost);  ///?????->because the sign of the electron vectors is reverted in order to construct a boost vector with a positive sign...
    kinematics::Q=CMBoost2.t();
    //  kinematics::Q=10.52; //die e energien die ich da hab sind on resonance. Dass hier ist aber continuum
    //    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    //    kinematics::secondElectronCM.boost(kinematics::CMBoost);

    m_histos.setFilenameStart(rFileName);

    srand(time(NULL));

  }



  int bToDDoubleStar::goodHadronB( ) const 
  {
    // initialize return value
    int b = 0;

    // get manager for evtcls_hadronic_flag
    Evtcls_hadronic_flag_Manager & evtMgr
      = Evtcls_hadronic_flag_Manager::get_manager();

    //   Evtcls_evtcls_flag_Manager & evtMgr2
    //     = Evtcls_evtcls_flag2_Manager::get_manager();

    // get flag for HadronB
    //      hadronic_flag(0) =  10 : old HadronA with R2<0.2
    //                       =  20 : old HadronA with R2>=0.2
    //      hadronic_flag(1) =  10 : new HadronA with R2<0.2
    //                       =  20 : new HadronA with R2>=0.2
    //      hadronic_flag(2) =  10 : HadronB with R2<0.2
    //                       =  20 : HadronB with R2>=0.2
    //      hadronic_flag(3) =  10 : new HadronA with #tracks>=5
    //      hadronic_flag(4) =  10 : HadronB with #tracks>=5
    //      hadronic_flag(5) =  10 : HadronC with R2<0.2
    //                       =  20 : HadronC with R2>=0.2
    Panther_ID id(1);
    Evtcls_hadronic_flag & hadflag = evtMgr(id);

    if ( hadflag.hadronic_flag(2) == 10 ||
	 hadflag.hadronic_flag(2) == 20    ) { b = 1; }

    //hadron J
    //      if ( hadflag.evtcls_flag2(2) >= 0  ) { b = 1; }

    if (b !=1)
      {
	//       printf("bad hadb: %d %d %d %d %d\n",hadflag.hadronic_flag(0),hadflag.hadronic_flag(2),hadflag.hadronic_flag(3),hadflag.hadronic_flag(4),hadflag.hadronic_flag(5));

      }

    // to select tau event
    //  if ( hadflag.hadronic_flag(4) != 0 ) { b += 2; }

    return b;
  }
  // begin_run function
  void bToDDoubleStar::begin_run(BelleEvent* evptr, int* status)
  {


    D_lundIds.push_back(PY_D0);
    D_lundIds.push_back(PY_D);
    D_lundIds.push_back(PY_DStar);
    D_lundIds.push_back(PY_DStar0);
    //D_1
    D_lundIds.push_back(10413);
    D_lundIds.push_back(PY_DStar_2);
    //D1Prime
    D_lundIds.push_back(20413);
    D_lundIds.push_back(PY_DStar0Plus);
    D_lundIds.push_back(PY_D_10);
    D_lundIds.push_back(PY_DStar_20);
    //D1Prime0
    D_lundIds.push_back(20423);
    D_lundIds.push_back(PY_DStar_00);

        cout <<" begin run  (bToDDoubleStar)" << endl;
    IpProfile::begin_run();
    eid::init_data();

    BeamEnergy::begin_run();
    double eler=BeamEnergy::E_LER();
    double eher=BeamEnergy::E_HER();


    HepLorentzVector pBeam=BeamEnergy::p_beam();

    //    cout <<"old beam: "<< kinematics::pBeam.px()<<" " << kinematics::pBeam.py()<< " " << kinematics::pBeam.pz()<< " " << kinematics::pBeam.pz() << kinematics::pBeam.t()<<endl;
    //    cout << "new: "<< pBeam.px() << " " << pBeam.py()<< "  "<< pBeam.pz() << "  " << pBeam.t()<<endl;
    kinematics::pBeam=BeamEnergy::p_beam();
    if(eler <3.0 || eher <7.0 || eler > 5.0 || eher > 9.0)
      {
	validRun=false;
	return;
      }
    else
      {
	validRun=true;
      }
    double theta(0.022);
    kinematics::cm=HepLorentzVector(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(eher*sin(theta),0.,eher*cos(theta),eher);
    kinematics::secondElectronCM=HepLorentzVector(0.,0.,-eler,eler);
    HepLorentzVector CMBoost2=kinematics::cm;
    //    CMBoost2.boost(-kinematics::CMBoost);  ///?????->because the sign of the electron vectors is reverted in order to construct a boost vector with a positive sign...
    kinematics::Q=CMBoost2.t();
    //    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    //    kinematics::secondElectronCM.boost(kinematics::CMBoost);
    //    cout <<"end beginrun " <<endl;
    return;
  }

  // hist_def function
  void bToDDoubleStar::hist_def()
  {
    Particle p;
    if(rFileName!=0)
      cout <<endl<<":::----- rFileName (handedness): " << rFileName <<endl<<endl;
    else
      cout <<endl <<":::--------No File Name specified (handedness)" <<endl<<endl;

    pTreeSaver=new TreeSaver();
    cout <<"do branching!! " <<endl;
    pTreeSaver->doBranching();

    //important that first all arrays are defined, F, than I 
#ifdef MC

#endif


#ifdef MC


#endif


#ifdef MC
    //the fields for the mc w/o acceptance
#endif


    //doesn't work anymore with arrays...
    //  pTreeSver->createNTuple(tm);
#ifdef MC

    //  pTreeSaver->addArrayPi0AsymmetryF("realPi0_gammaAsymmetry");
#endif
  }
  int dDoubleStarId=0;
  // event function
  void bToDDoubleStar::event(BelleEvent* evptr, int* status)
  {

    treeData.initData();
    //0.35% systematics per track, efficiency is 0.8, so overall 0.35%/0.8
    //sum up relative uncertainties for tracks for the signal, 
    //uncertainties on the tagside should be in the tag uncertainty
    float relTrackingSys=0.0035/0.8;
    float sTracking=0.0;
    int numTracks=0;

    float sB_BR=0.0;
    float sD_BR=0.0;

    float sXSecWeight=0.0;
    float xPIDWeight=0.0;
    //seems to be global, so we don't see it set
    dssIdx=-1;
    lType=-1;
    //    cout <<"mc type: "<< mc_type<<endl;
    //    cout <<"reading event (bToDDoubleStar)" << endl;
    foundChargedDecIds.clear();
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    float  visEnergyOnFile=hadronInfo_mgr.begin()->Evis();
    kinematics::E_miss=kinematics::Q-visEnergyOnFile;

    //for data make sure that we have a sensible default value
    treeData.initData();
    //should not be necessary
    treeData.B_DecayCorr=1.0;
    //    treeData.D_DecayCorr=1.0;

    bool bgFlag=false;
    bool foundDDStarFlag=false;
    //for the output of the decay signature query
    sig_FoundDDoubleStar=false;
    sig_numLeptons=0;
    sig_numKaons=0;
    sig_numRhos=0;
    sig_numPions=0;
    sig_numD=0;
    sig_DCharge=0;
    sig_numDStar=0;
    sig_numPi0=0;
    sig_numBaryons=0;
    sig_dStar_2S=false;
    sig_d_2S=false;

    sigDStarLNu=0;
    sigDStarPiLNu=0;
    sigDStarPiPiLNu=0;

    sigDLNu=0;
    sigDPiLNu=0;
    sigDPiPiLNu=0;

    sigResDStarLNu=0;
    sigResDStarPiLNu=0;
    sigResDStarPiPiLNu=0;

    sigResDLNu=0;
    sigResDPiLNu=0;
    sigResDPiPiLNu=0;

    ///add the decays for which we want to do branching ratio and FF corrections
    br_sig_D0LNu=0;
    //    br_sig_D0=0;
    br_sig_DLNu=0;

    br_sig_DStarLNu=0;

    br_sig_DStar0LNu=0;
    
    //10413 -->PY_D_1
    br_sig_D1LNu=0;

    //415 PY_DSTAR_2
    br_sig_D2LNu=0;

    //20413
    br_sig_D1PrimeLNu=0;

    //10411# PY_DStar0Plus 
    br_sig_D0StarLNu=0;

    //10423 PY_D_10 
    br_sig_D10LNu=0;


    //425  PY_DStar_20
    br_sig_D20LNu=0;

    //pythia 20423
    br_sig_D1Prime0LNu=0;

    //10421 PY_DStar_00 
    br_sig_D0Star0LNu=0;
    //---
    int br_sigs[12];

    found_2SD=false;
    found_2SD_Star=false;

    const double m_pi0=0.1349766;
    vector<int> foundDecIds;

    int numNu=0;
    vector<int> decIds;
    vector<int> pids;
    int bMesonId;

    foundSinglePionDecay=false;
    foundDlNu=false;
    noDRec=false;

    if(!validRun)
      {
	return;
      }

    /////for xcheck
    AnaBrecon brecon;
    evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    runNr=Belle_event_Manager::get_manager().begin()->RunNo();
    expNr=Belle_event_Manager::get_manager().begin()->ExpNo();
    //    cout <<" evt: " << evtNr << " runNr " << runNr <<endl;
    kinematics::runNr=runNr;
    kinematics::evtNr=evtNr;

    float XSectionWeight=1.0;
    float sXSectionWeight=0.0;

    float lumiWeight=1.0;
    float sLumiWeight=0.0;

    float ffWeight=1.0;
    float sFFWeight=0.0;

    //pid efficiencies only on the signal side (tag correction is done separately)
    float pidWeight=1.0;
    float sPidWeight=0.0;
    //br correction done independently if signal or tag side

    if(m_mc)
      {
	if(foundUnwantedDDoubleStar())
	  {
	    if(kickOut2S)
	      return;
	  }
	Mctype_Manager &mctype_m = Mctype_Manager::get_manager();

	//this doesn't work with Jan's (Huschle) MC, correction factor is about 2%, let's set it to one for now
	//for the usual one we now set the flag 'excludeSignal', so we can use that to see if we can ask for the mc_type

	if(excludeSignal>0)
	  mc_type = mctype_m.begin()->type();

	pair<float,float> ret(1.0,0.0);
	ret=pidCorrections.getXSectionCorrection(m_mc,mc_type);
	XSectionWeight=ret.first;
	sXSectionWeight=ret.second;
	////--->why??	XSectionWeight=1.0;

	treeData.experiment=expNr;
	ret=pidCorrections.getLumiCorrection(expNr);
	lumiWeight=ret.first;
	sLumiWeight=ret.second;

	//find b meson id in the simulation (also does the determination of the subprocess for the BR etc correction
	foundDPiPi=checkForDPiPi(bMesonId,foundSinglePionDecay, foundDlNu);
	//	if(dDoubleStarId==100423)
	//	  cout <<"foundres  dd double star, --> sigres: "<< sigResDPiLNu<<endl;

	//the ff corrections can only be made after a call to 'foundDPiPi', since this is where we check for the correct decay
	float l_ffWeight, l_ffStat,l_ffSys;
	FFCorrections::ffRet mFfRet;	
	FFCorrections::ffRetDds mFfRetDds;	
	if(dssIdx>-1 && lType>-1)
	  {
	    //	    cout <<"both, D**lnu and D(*)lnu??" <<endl;
	  }
	if(lType>-1)
	  {
	    //	    cout <<"get weight for D " <<endl;
	    ffCorrections.getWeightD(ff_pL,q2,lType,l_ffWeight,l_ffStat,l_ffSys,mFfRet);
	    treeData.FFDCorrection[treeData.numFFDCorr]=l_ffWeight;
	    treeData.sFFDCorrection[treeData.numFFDCorr]=(l_ffStat*l_ffStat+l_ffSys*l_ffSys);
	    treeData.D_q2BinFF[treeData.numFFDCorr]=mFfRet.q2Bin;
	    treeData.D_pBinFF[treeData.numFFDCorr]=mFfRet.pBin;
	    treeData.D_TypeBinFF[treeData.numFFDCorr]=mFfRet.dType;
	    treeData.numFFDCorr++;

	  }
	if(dssIdx>-1)
	  {
	    //	    cout <<"get weight for DStar, idx: "<< dssIdx <<endl;

	    ffCorrections.getWeightDDStar(w,cosTheta,dssIdx,l_ffWeight,l_ffStat,l_ffSys,mFfRetDds);
	    treeData.FFDdsCorrection[treeData.numFFDdsCorr]=l_ffWeight;
	    treeData.sFFDdsCorrection[treeData.numFFDdsCorr]=(l_ffStat*l_ffStat+l_ffSys*l_ffSys);
	    treeData.Dds_cosTBinFF[treeData.numFFDdsCorr]=mFfRetDds.cosTBin;
	    treeData.Dds_wBinFF[treeData.numFFDdsCorr]=mFfRetDds.wBin;
	    treeData.Dds_TypeBinFF[treeData.numFFDdsCorr]=mFfRetDds.ddsType;
	    treeData.numFFDdsCorr++;


	    //	    cout <<"adding " <<((l_ffStat*l_ffStat+l_ffSys*l_ffSys)/(l_ffWeight*l_ffWeight)) <<  " with stat:  "<< l_ffStat <<" sys:  " << l_ffSys <<" weight : "<< l_ffWeight <<endl;
	    //	    cout <<"ffWeight now: "<< ffWeight <<endl;
	  }
	if(foundDPiPi || foundSinglePionDecay || foundDlNu)
	  {
	    //print
	    //	foundDPiPi=checkForDPiPi(bMesonId,true);
	    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();

	    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
	      {
		if(gen_it->get_ID()!=bMesonId)
		  continue;
		
		//call recCheck with the meson that has the desired decay
		recCheck((*gen_it),decIds,pids,numNu);
		//	    cout <<"rec check has " << decIds.size() <<" dec dis and " << pids.size()<<" pids, number of Nu: "<< numNu<<endl;
	      }
	  }
      }
    //    printMCList();



    chargedIds.clear();
    pi0Ids.clear();
    gammaIds.clear();

    //the signal decay candidates (after removing decay products from the best B
    //from Robin...
    //was -3 but put lower so we can explore phase space for fom
    const double _log_nbout_min=-4;
    //    const double _tag_Mbc=5.25;

    //    cout <<" valid run " << validRun<<endl;
    //    kinematics::expNr=expNr;

    treeData.recBToDlNuPiPi=0;

    if(!IpProfile::usable())
      return;
    if(!(test%1000))
      {
	//	      cout << "evt " <<test <<endl;
	//	        cout << "nr " <<evtNr <<endl;
      }
    test++;



    //#ifndef MC
    if(!goodHadronB())
      {
	return;
      }
    //#endif
    char ptypName[200];
    int lundPC=-1;


    Mdst_charged_Manager& mdst_chr_Mgr=Mdst_charged_Manager::get_manager();
    //  Mdst_klong_Manager& mdst_klong_Mgr=Mdst_klong_Manager::get_manager();
    Mdst_trk_Manager& mdst_trk_Mgr=Mdst_trk_Manager::get_manager();
    atc_pid selKPi(3,1,5,3,2);  //K/pi separation
    atc_pid selPK(3,1,5,4,3); //proton kaon separation
    atc_pid selPiP(3,1,5,2,3); //pion proton separation
    atc_pid selKP(3,1,5,3,4);

    ////---
    Ekpfullrecon_Manager &ekpfullrecon_m= Ekpfullrecon_Manager::get_manager();
    Ekpfullrecon bestEKP_B;
    int bestId;
    double bestLogProb=-1000000;
    double bestProb=-1000000;
    double bestDeltaE=-10000;
    double bestTagM=-10000;
    int bestDecay=-1000;
    int bestTag=-1000;


    if(foundDPiPi && ekpfullrecon_m.size()!=0)
      {
	//		cout <<" got B and the dpipi" << endl;
      }
    if(ekpfullrecon_m.size()!=0)
      {

	//	cout <<"ekpfull found something" <<endl;
      }
    else
      {
	//	cout <<"ekp didn't find anything " << endl;
      }
    //    cout <<"looking for ekpfullrecon " <<endl;
    //first bind best B meson and decay products
    for(std::vector<Ekpfullrecon>::iterator it=ekpfullrecon_m.begin();it!=ekpfullrecon_m.end();it++)
      {
	Ekpfullrecon& ekpfullrecon=*it;
	double nbout=ekpfullrecon.NBout();
	double log_nbout=log10(nbout);
	double tab_Mbc=ekpfullrecon.Mbc();
	int decay =ekpfullrecon.decay();
	float DeltaE=ekpfullrecon.DeltaE();
	//gives the pythia code of tag
	int tag=abs(ekpfullrecon.tag_id());
	if(log_nbout< _log_nbout_min)
	  continue;
	if(log_nbout<bestLogProb)
	  continue;
	else
	  {
	    bestLogProb=log_nbout;
	    bestProb=nbout;
	    bestEKP_B=ekpfullrecon;
	    bestDeltaE=DeltaE;
	    bestDecay=decay;
	    bestTagM=tab_Mbc;
	    bestTag=tag;
	  }

	//I guess M_bc is beam constrained mass...
	//	cout<<"ex: nbout "<< nbout << " .. tab_Mbc:" << tab_Mbc<<" decay: " << decay<<", deltaE: " << DeltaE <<" tagCorr: " << tagCorr <<endl;
	Particle & bcand = const_cast<Particle &>(brecon.getParticle( (int)ekpfullrecon.get_ID() ));
	bestId=(int) ekpfullrecon.get_ID();
	//	cout <<"found candidate " << bcand.pType().name() <<" with " << bcand.nChildren()<<", pcode: "<< bcand.pType().lund()<<endl;
	histoNbOut->Fill(nbout);
	histoChargedBTag_M->Fill(tab_Mbc);
	histoChargedBTag_dE->Fill(DeltaE);

	histoB0Tag_M->Fill(tab_Mbc);
	histoB0Tag_dE->Fill(DeltaE);

	int numBChildren = bcand.nChildren();
	for(int i=0;i<numBChildren;i++)
	  {
	    const Particle& p=bcand.child(i);
	    if(p.mdstCharged())
	      {
		//		cout <<"looking at charged.." <<endl;
		//	      cout <<"charged id: "<< p.mdstCharged().get_ID()<<endl;
	      }
	    if(p.nChildren()>0)
	      {
		for(int j=0;j<p.nChildren();j++)
		  {
 		    const Particle& p2=p.child(j);
		    //		    cout << "---> " << p2.pType().name() <<" has " << p2.nChildren() <<" own children " <<endl;
		    for(int k=0;k<p2.nChildren();k++)
		      {
			const Particle& p3=p2.child(k);
			//	cout <<"--->--->"<<p3.pType().name()<<endl;
		      }
		  }
	      }
	  }

	//	cout <<"got " << chargedIds.size() << " charged Ids, " << pi0Ids.size() <<" pi0s, " << gammaIds.size() <<" gammas " << endl;
      }
    //need good b
    if(bestLogProb<_log_nbout_min)
      {
	//	if(dDoubleStarId==100423)
	exitEvent();
	return;
      }
    //According to Christoph, there are no non-resonant D(*)pilnu decays...
    if(m_mc && (sigDPiLNu || sigDStarPiLNu))
      {
	if(dDoubleStarId==100423)
	  {
	    //	    cout <<"non res  d(*)pi" <<endl;
	  }
	exitEvent();
	return;
      }
	if(dDoubleStarId==100423)
	  {
	    if(sigResDPiLNu)
	      {
		//	cout <<" there and signal" <<endl;
	      }
	  }
    
    //    cout <<"best log prob: "<< bestLogProb<<endl;
    Particle & bestBcand = const_cast<Particle &>(brecon.getParticle( (int)bestEKP_B.get_ID() ));
    bestBPx=bestBcand.px();
    bestBPy=bestBcand.py();
    bestBPz=bestBcand.pz();


    //   cout <<"looking at b cand" <<endl;

    //    if(bestBcand.mdstVee2())
    //      csout <<"b has vee2 " <<endl;
    //    if(bestBcand.mdstCharged())
    //      cout <<"b has charged " <<endl;


    //    Gen_hepevt hepEvt=get_hepevt(p->mdst);

    float tagCorr=tagcorr(bestDecay,bestProb);
    treeData.mBTag=bestTagM;
    treeData.deltaETag=bestDeltaE;
    treeData.tagDecay=bestDecay;
    treeData.tagId=bestTag;
    treeData.logProb=bestLogProb;
    treeData.foundDPiPi=foundDPiPi;
    treeData.recBToDlNuPi=foundSinglePionDecay;
    treeData.tagCorr=tagCorr;
    // sig_dStar_2S, sig_d_2S);
    treeData.found_2SD=sig_d_2S;
    treeData.found_2SD_Star=sig_dStar_2S;
    treeData.foundAnyDDoubleStar=sig_FoundDDoubleStar;
    treeData.DDStarMass_mc=mcMassDDStar;
    treeData.sig_numPions=sig_numPions;
    treeData.sig_numD=sig_numD;
    treeData.sig_numDStar=sig_numDStar;
    treeData.sig_numKaons=sig_numKaons;

    //otherwise we will already have set in when we looked for the non-resonant decays
    if(sig_FoundDDoubleStar)
      {
	treeData.mcDCharge=sig_DCharge;
	treeData.mcIsDStar=0;
	if(sig_numDStar==1 && sig_numD==0)
	  {
	    treeData.mcIsDStar=1;
	  }
      }
    treeData.sig_numPi0=sig_numPi0;
    treeData.sig_numLeptons=sig_numLeptons;
    treeData.sig_numBaryons=sig_numBaryons;
    treeData.tagOverlapFractionCharged=overlapFractionCharged;
    treeData.tagOverlapFractionPi0=overlapFractionPi0;
    //    if((sigDStarLNu || sigDStarPiLNu || sigDStarPiPiLNu )&& sig_FoundDDoubleStar)
      {
	//	cout <<"found dlnu and ddouble star!!!!!" <<endl;
      }

    treeData.sigDLNu=sigDLNu;
    treeData.sigDPiLNu=sigDPiLNu;
    treeData.sigDPiPiLNu=sigDPiPiLNu;

    treeData.sigDStarLNu=sigDStarLNu;
    treeData.sigDStarPiLNu=sigDStarPiLNu;
    treeData.sigDStarPiPiLNu=sigDStarPiPiLNu;



    treeData.sigResDLNu=sigResDLNu;
    treeData.sigResDPiLNu=sigResDPiLNu;
    treeData.sigResDPiPiLNu=sigResDPiPiLNu;

    treeData.sigResDStarLNu=sigResDStarLNu;
    treeData.sigResDStarPiLNu=sigResDStarPiLNu;
    treeData.sigResDStarPiPiLNu=sigResDStarPiPiLNu;


    treeData.sig_dStar_2S=sig_dStar_2S;
    treeData.sig_d_2S=sig_d_2S;

    getDecayIds(bestBcand,chargedIds,pi0Ids,gammaIds);
    //pi0s and gamma's seem to be balanced
    //    cout <<"b candidate has " << pi0Ids.size() <<" pi0s and " << gammaIds.size()<<" photons.." <<endl;

    //now look at the other particles...
    //find possible Ks..
    Mdst_vee2_Manager &vee2_m=Mdst_vee2_Manager::get_manager();
    for(vector<Mdst_vee2>::iterator vee_it=vee2_m.begin();vee_it!=vee2_m.end();vee_it++)
      {
	//not a K short...
	if(vee_it->kind()!=1)
	  continue;

	/// also quality checks?
	//the below is to make sure that none of the daughter particles are part of the tag or of a Ks that we already found
	if(chargedIds.find(vee_it->chgd(0).get_ID())!=chargedIds.end()){
	  //	  cout <<"charged daugher of Ks is part of tag " <<endl;
	  continue;
	}
	if(chargedIds.find(vee_it->chgd(1).get_ID())!=chargedIds.end()){
	  //	  cout <<"charged daugher of Ks is part of tag " <<endl;
	  continue;
	}
	//second parameter keeps relation
	Particle* p=new Particle(*vee_it,true);
	FindKs findks;

	findks.candidates(*vee_it,IpProfile::position());
	if(findks.goodKs())
	  {
	    KsCandidates.push_back(p);
	    histoKs->Fill(p->mass());
	    chargedIds.insert(vee_it->chgd(0));
	    chargedIds.insert(vee_it->chgd(1));
	    if(m_mc)
	      {
		pair<float,float> ret;
		PIDCorrections::pidRet mPidRet;
		ret=pidCorrections.getWeight(-1,310,p->ptot(), p->p3().theta(),expNr,runNr, mPidRet);
		treeData.KsCorrection[treeData.numKsCorr]=ret.first;
		treeData.sKsCorrection[treeData.numKsCorr]=ret.second;
		treeData.KsCorrMomBin[treeData.numKsCorr]=mPidRet.momBin;
		treeData.KsCorrThetaBin[treeData.numKsCorr]=mPidRet.thetaBin;
		treeData.numKsCorr++;

		pidWeight*=ret.first;
		sPidWeight+=(ret.second*ret.second)/(ret.first*ret.first);


		//	      if(pidCorrections.getWeight(-1,310,p->ptot(), p->p3().theta(),expNr,runNr)>2.0)
		//		if(pidCorrections.getWeight(-1,310,p->ptot(), p->p3().theta(),expNr,runNr).first<0.0)
		  {
		    //		    cout <<"ks pidweight : " << pidCorrections.getWeight(-1,310,p->ptot(), p->p3().theta(),expNr,runNr) <<endl;
		  }
	      
	      }
	  }
	else
	  {
	    delete p;
	  }

      }


    //    cout <<"there are " << mdst_chr_Mgr.size() << " charge tracks in mdst_chr " <<endl;
    for(Mdst_charged_Manager::iterator chr_it=mdst_chr_Mgr.begin();chr_it!=mdst_chr_Mgr.end();chr_it++)
      {
	//	cout <<" looking at charged.." << endl;
	//to check if the pions from in the signal can be reconstructed...
	if(m_mc)
	  foundChargedDecIds.insert(chr_it->get_ID());

	//part of the tag
	if(chargedIds.find(chr_it->get_ID())!=chargedIds.end()){
	  //	 	  cout <<"charged track is part of tag " <<endl;
	  	  continue;
	}


	//      cout <<"looking at " <<(*chr_it).p(0) <<" " << (*chr_it).p(1) <<" " << (*chr_it).p(2) <<endl;
	/////tracking systematics
	sTracking+=relTrackingSys*relTrackingSys;
	numTracks++;
	double m_mass=m_pi;
	int massHyp=2;
	bool isLepton=false;

	//default pion like Dmitri's cuts
	bool isPionKaon=true;
	//as opposed to default pion...
	bool positivelyIdentified=false;
	strcpy(ptypName,"unknown");
	double charge=(*chr_it).charge();
	//defaults...
	if(charge>0)
	  {
	    strcpy(ptypName,"PI+");
	  }
	else
	  {
	    strcpy(ptypName,"PI-");
	  }
	
	//      HepLorentzVector hepvec;
	//immer daran denken in die richtung -boostvector zu boosten! ... nein ist anscheinend schon in der boost vector definition drin...
	eid sel_e(*chr_it);
	double mu_id=0;
	Muid_mdst muID(*chr_it);
	if(muID.Chi_2()>0)
	  mu_id=muID.Muon_likelihood();

	double atcKPi=selKPi.prob(*chr_it);
	double atcKP=selKP.prob(*chr_it);
	double atcPiP=selPiP.prob(*chr_it);
	float e_cut=0.8;
	float mu_cut=0.9;
	double e_id=sel_e.prob(3,-1,5);


	//	cout <<"atcKPi: " << atcKPi <<", atcKP " << atcKP << " atcPiP: "<< atcPiP <<" e_id: "<< e_id <<" mu: "<< mu_id <<endl;

	if(DEBUG_EVENT==evtNr)
	  {
	    //	    cout <<"pid kpi: " << atcKPiAlt <<" pid KP: " << atcKPAlt << " e_id: " << e_id << " mu_id: " << mu_id <<endl;
	  }

	if(e_id>e_cut&& mu_id<0.9)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//	      cout <<"is electron" <<endl;
	      }
	    
	    m_mass=m_e;
	    massHyp=0;
	    positivelyIdentified=true;
	    isLepton=true;

	    if(charge>0)
	      {
		strcpy(ptypName,"E+");	   
		treeData.leptonId=11;
	      }
	    else
	      {
		treeData.leptonId=-11;
		strcpy(ptypName,"E-");
	      }

	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	  }
	//used to be mu_id > e_cut 
	if(mu_id>mu_cut && e_id<e_cut)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//		cout <<"is muon" <<endl;
	      }
	    m_histos.hPidMuPi->Fill(chr_it->trk().pid_mu(),chr_it->trk().pid_pi());
	    m_mass=m_muon;
	    massHyp=1;
	    if(charge>0)
	      {
		treeData.leptonId=13;
		strcpy(ptypName,"MU+");
	      }
	    else
	      {
		treeData.leptonId=-13;
		strcpy(ptypName,"MU-");
	      }

	    positivelyIdentified=true;
	    isLepton=true;
	    m_histos.hPidMu->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	  }
	if(mu_id>0.9&& e_id>e_cut)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//	      cout <<"is electron" <<endl;
	      }
	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	    //this mass is only used for the dr/dz fits, the mass hypothesis for later is pion
	    m_mass=m_e;
	    isLepton=true;
	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	  }


	if(!isLepton)
	  {
	    if(atcKPi>0.6) //kaon
	      {
		m_mass=m_k;
		massHyp=3;
		positivelyIdentified=true;
		isPionKaon=true;
		if(charge>0)
		  {
		    strcpy(ptypName,"K+");
		  }
		else
		  {
		    strcpy(ptypName,"K-");
		  }
		m_histos.hPidK->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_pi());
	      }
	    else
	      {
		//second condition accoriding to Ami's cuts... default to pion, purity goes down
		if(atcKP<0.2 && atcPiP<0.2)
		  {
		    isPionKaon=false;
		    m_mass=m_pr;
		    massHyp=4;
		    positivelyIdentified=true;
		    m_histos.hPidPr->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_p());
		    m_histos.hPidPrPi->Fill(chr_it->trk().pid_p(),chr_it->trk().pid_pi());
		    if(charge>0)
		      strcpy(ptypName,"P+");
		    else
		      strcpy(ptypName,"AP+");
		  }
		else  //pion
		  {
		    //default mass assignment if nothing is found is pion anywasy....
		    //			if(atcKPi<0.3) //Ami has 0.4...
		    if(atcKPi<0.4)
		      {
			massHyp=2;
			positivelyIdentified=true;
			m_mass=m_pi;
			isPionKaon=true;
			if(charge>0)
			  {
			    strcpy(ptypName,"PI+");
			  }
			else
			  {
			    strcpy(ptypName,"PI-");
			  }
			/*		      m_histos.hPidPi->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_pi());
					      m_histos.hPidPiMu->Filpl(chr_it->trk().pid_pi(),chr_it->trk().pid_mu());
					      m_histos.hPidPiE->Fill(chr_it->trk().pid_pi(),chr_it->trk().pid_e());
					      m_histos.hPidPiPr->Fill(chr_it->trk().pid_pi(),chr_it->trk().pid_p());*/

		      }
		  }
	      }
	  }
	else
	  {
	    //	cout <<"e_id : "<< e_id <<" mu_id: "<< mu_id <<endl;
	  }
      
	//default pion, good enough
	//	if(!positivelyIdentified)
	//	  continue;

	double dr,dz, refitPx, refitPy, refitPz;
	//	Momentum mom=p->momentum();
	getDrDz(chr_it, massHyp,dr,dz, refitPx, refitPy, refitPz, m_mass);

	v_vertexR.push_back(dr);
	v_vertexZ.push_back(dz);

	///
	//      cout <<"looking at " <<(*chr_it).p(0) <<" " << (*chr_it).p(1) <<" " << (*chr_it).p(2) <<endl;
	if ( fabs(dr) > cuts::vertexR )//slides from kibayashi 
	  {
	    continue;
	  }
	if ( fabs(dz) > cuts::vertexZ ) 
	  {
	    continue;//used to be 4
	  }
	//check if electron or muon are in reasonable acceptance
	//e: 17 < theta < 150, mu: 25< theta< 145

	Particle* p=new Particle(*chr_it,string(ptypName));
	float deg=180*(p->p3().theta()/TMath::Pi());
	//	cout <<" converted " << p->p3().theta() <<" to  " << deg <<endl;
	if(fabs(p->pType().lund())==PY_E)
	  {
	    if(deg < 17 || deg> 150)
	      continue;
	    addBremsPhoton(p);
	  }
	if(fabs(p->pType().lund())==PY_MU)
	  {
	    if(deg < 25 || deg >145)
	      continue;
	  }


	//------>
	// do mis pid weighting for charged particles which we id'ed by now
	//note that particles that are part of the tag, and that shouldn't be weighted are not part of this loop (we checked earlier)
	//pi0 and ks misid later...
      

	if(m_mc)
	  {
	    Gen_hepevt hepEvt=get_hepevt(*chr_it);  
	    if((fabs(p->pType().lund()==PY_E) || fabs(p->pType().lund())==PY_MU )&& (p->pType().lund()!=hepEvt.idhep()))
	      {
		//		cout <<"thought we had " << p->pType().lund() <<" but we have " << hepEvt.idhep()<<endl;
	      }
	    PIDCorrections::pidRet mPidRet;
	    mPidRet.misIdType=-1;
	    pair<float,float> ret=pidCorrections.getWeight(hepEvt.idhep(), p->pType().lund(), p->ptot(), p->p3().theta(),expNr,runNr,mPidRet);
	    //correction actually applied
	    if(mPidRet.misIdType>-1)
	      {
		treeData.ChargedCorrection[treeData.numChargedCorr]=ret.first;
		treeData.sChargedCorrection[treeData.numChargedCorr]=ret.second;
		treeData.ChargedCorrMomBin[treeData.numChargedCorr]=mPidRet.momBin;
		treeData.ChargedCorrThetaBin[treeData.numChargedCorr]=mPidRet.thetaBin;
		treeData.ChargedCorrSVDBin[treeData.numChargedCorr]=mPidRet.svdBin;
		treeData.ChargedCorrMisIdType[treeData.numChargedCorr]=mPidRet.misIdType;
		treeData.numChargedCorr++;
	      }

	    pidWeight*=ret.first;
	    sPidWeight+=(ret.second*ret.second)/(ret.first*ret.first);
	    //	  if(pidCorrections.getWeight(hepEvt.idhep(), p->pType().lund(), p->ptot(), p->p3().theta(),expNr,runNr)>2.0)
	    //	    if(pidCorrections.getWeight(hepEvt.idhep(), p->pType().lund(), p->ptot(), p->p3().theta(),expNr,runNr).first<0.0)
	      {
		//	      cout <<"charged pidweight : " << pidCorrections.getWeight(hepEvt.idhep(), p->pType().lund(), p->ptot(), p->p3().theta(),expNr,runNr) << ", "<<hepEvt.idhep() <<" rec: "<< p->pType().lund() <<" p toto : " << p->ptot() <<" theta: " << p->p3().theta() << " exp: "<< expNr <<" run  " << runNr <<endl;
	      }

	  }



	  
	//	cout <<"mom px: "<< mom.p().px()<<" py: " << mom.p().py() <<" pz: "<< mom.p().pz() <<" t: "<< mom.p().t()<<endl;
	//	cout <<"par px: "<< p->p().px()<<" py: " << p->p().py() <<" pz: "<< p->p().pz() <<" t: "<< p->p().t()<<endl;
       	//p->momentum(mom);

	//	mom.p().setPx(refitPx);
	//	p->SetPx(refitPx);

	//check if this is part of the d pipi lnu system
	if(m_mc)
	  {
	    Gen_hepevt hepEvt=get_hepevt(p->mdstCharged());


	    if(foundDPiPi)// && p->lund()==PY_PI)
	      {
		for(int i=0;i<decIds.size();i++)
		  {
		    if(decIds[i]==hepEvt.get_ID() && p->pType().lund()==pids[i])
		      {
			foundDecIds.push_back(hepEvt.get_ID());
		      }
		  }

		//	    cout <<"found charged with hepevt id: " << hepEvt.get_ID()<<" and lund id: "<< p->lund()<<endl;
		//	  cout <<"Px: " << p->px() <<" "<<p->py() <<" " << p->pz()<<endl;
		//	  cout <<"hepvt momenta: " << hepEvt.PX() <<" " << hepEvt.PY() <<" " << hepEvt.PZ()<<endl;
	      }
	  }

	//has to be in parantheses
	p->userInfo(*(new ParticleInfo()));
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(p->userInfo());
	Ptype& m_pt=p->pType();
	if(isLepton){
	  leptonCandidates.push_back(p);
	}
	else{
	  if(fabs(p->pType().lund())==PY_PI)
	    {
	      chargedPiCandidates.push_back(p);
	    }
	  if(fabs(p->pType().lund())==PY_K)
	    {
	      chargedKCandidates.push_back(p);
	    }
	  //probably not a good event...
	  if(!isPionKaon)
	    {
	      otherChargedTracks.push_back(p);
	    }
	}
      }

    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
    Mdst_pi0_Manager &pi0_mgr=Mdst_pi0_Manager::get_manager();
    for(std::vector<Mdst_pi0>::const_iterator i =pi0_mgr.begin();i!=pi0_mgr.end();i++)
      {
	const Mdst_pi0& pi0=*i;
	int id =(int)pi0.get_ID();
	//	cout <<"pi0 children: " << pi0.nChildren()<<endl;

	if(pi0Ids.find(id)!=pi0Ids.end())
	  {
	    //	    cout <<"found pi0 which is part of B candidate .." <<endl;
	    continue;
	  }
	
	if(m_mc)
	  {
	    Particle* pi0part=new Particle(*i,"PI0");
	    Gen_hepevt hepEvt=get_hepevt(pi0part->mdstPi0());
	    delete pi0part;

	    if(foundDPiPi)
	      {
		if(find(decIds.begin(),decIds.end(),hepEvt.get_ID())!=decIds.end())
		  {
		    foundDecIds.push_back(hepEvt.get_ID());
		    //		cout <<"found pi0 decay particle.." <<endl;
		  }
	      }
	  }

	double px=pi0.px();
	double py=pi0.py();
	double pz=pi0.pz();


	Mdst_ecl_aux &aux1 =eclaux_mgr(Panther_ID(pi0.gamma(0).ecl().get_ID()));
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25_1 =aux1.e9oe25();
	Mdst_ecl_aux &aux2 =eclaux_mgr(Panther_ID(pi0.gamma(1).ecl().get_ID()));
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25_2 =aux2.e9oe25();
	double mass=pi0.mass(); //mass before fitting ???
	///	if(mass>0.15 || mass<0.12)
	//almost the same as above...
	if(fabs(mass-pi0Mass)>0.015)
	  continue;
	float pLab=sqrt(px*px+py*py+pz*pz);
	//      cout <<"pi0mass: "<< mass <<endl;≈

	float g1Energy= sqrt(pi0.gamma(0).px()*pi0.gamma(0).px()+pi0.gamma(0).py()*pi0.gamma(0).py()+pi0.gamma(0).pz()*pi0.gamma(0).pz());
	float g2Energy= sqrt(pi0.gamma(1).px()*pi0.gamma(1).px()+pi0.gamma(1).py()*pi0.gamma(1).py()+pi0.gamma(1).pz()*pi0.gamma(1).pz());
	//	cout <<"pi0 gamma1: "<< g1Energy <<" gamma2: "<< g2Energy <<endl;

	//	if(abs((g1Energy-g2Energy)/(g1Energy+g2Energy))>cuts::maxPi0GAsym)
	//	  continue;
	//let's make this 100 MeV rather than 50 (see Charlotte's study)
	//	if(g1Energy < 0.1 || g2Energy < 0.1)



	if(!checkGammaEnergy(pi0.gamma(0).px(), pi0.gamma(0).py(), pi0.gamma(0).pz(),g1Energy))
	  continue;
	if(!checkGammaEnergy(pi0.gamma(1).px(), pi0.gamma(1).py(), pi0.gamma(1).pz(),g2Energy))
	  continue;
	Particle* p=new Particle(pi0);

	double confLevel;
	//       cout <<" pi0 px: "<< pi0.px()<<" py: "<< pi0.py() <<" pz: "<< pi0.pz() <<endl;
	HepPoint3D pi0DecPoint;
	HepSymMatrix pi0ErrMatrix;

	setGammaError(p->child(0),IpProfile::position(), IpProfile::position_err_b_life_smeared());
	setGammaError(p->child(1),IpProfile::position(), IpProfile::position_err_b_life_smeared());
	//	cout <<"init mass; " << pi0.mass()<<endl;
	if(!doKmFit(*p,  confLevel,0,m_pi0))
	  {
	    	    continue;
	  }
	//	cout << "after km fit pi0 px: "<< pi0.px()<<" py: "<< pi0.py() <<" pz: "<< pi0.pz() <<endl;
	//	cout <<"mass; " << pi0.mass()<<endl;

	p->userInfo(*(new ParticleInfoMass()));

	ParticleInfoMass& pinf=dynamic_cast<ParticleInfoMass&>(p->userInfo());
	pinf.gammaE1=g1Energy;
	pinf.gammaE2=g2Energy;
	pinf.e9oe25_1=e9oe25_1;
	pinf.e9oe25_2=e9oe25_2;
	pi0Candidates.push_back(p);

	//      v_drAll.push_back();
	if(m_mc)
	  {
	    PIDCorrections::pidRet mPidRet;
	    pair<float, float> ret=pidCorrections.getWeight(-1,111,p->ptot(), p->p3().theta(),expNr,runNr,mPidRet);
	    pidWeight*=ret.first;
	    sPidWeight+=(ret.second*ret.second)/(ret.first*ret.first);
	    treeData.pi0Correction[treeData.numPi0Corr]=ret.first;
	    treeData.sPi0Correction[treeData.numPi0Corr]=ret.second;
	    treeData.pi0MomBin[treeData.numPi0Corr]=mPidRet.momBin;
	    treeData.numPi0Corr++;


	    //	  if(pidCorrections.getWeight(-1,111,p->ptot(), p->p3().theta(),expNr,runNr)>2.0)
		//	    if(pidCorrections.getWeight(-1,111,p->ptot(), p->p3().theta(),expNr,runNr).first<0.0)
	      {
		//		cout <<"pi0 pidweight : " << pidCorrections.getWeight(-1,111,p->ptot(), p->p3().theta(),expNr,runNr).first <<endl;
	      }
	  }
      }

////    cout <<" we have " << pi0Candidates.size() <<" pi0s before: " << endl;
////    for(vector<Particle*>::iterator itp=pi0Candidates.begin();itp!=pi0Candidates.end();itp++)
////      {
////	cout <<"energy: " << (*itp)->p().rho();
////      }
////    cout <<endl;
    cleanPi0s();
////       cout <<" we have " << pi0Candidates.size() <<" pi0s after: " << endl;
////   for(vector<Particle*>::iterator itp=pi0Candidates.begin();itp!=pi0Candidates.end();itp++)
////     {
////	cout <<"energy: " << (*itp)->p().rho();
////     }
////   cout <<endl;
///    ///now we should have all the candidates..
    //    cout <<"we have " << chargedPiCandidates.size() << " pions " << chargedKCandidates.size() <<" kaons " << pi0Candidates.size() <<" pi0 " << KsCandidates.size() <<" Ks " << leptonCandidates.size() <<"leptons " << otherChargedTracks.size() <<" others " <<endl;


    Mdst_gamma_Manager& gamma_mgr=Mdst_gamma_Manager::get_manager();
    int gammaCount=0;
    for(std::vector<Mdst_gamma>::const_iterator i =gamma_mgr.begin();i!=gamma_mgr.end();i++)
      {

	if(passGammaQACuts(i))
	  {
	    const Mdst_gamma& gam=*i;
	    double px=gam.px();
	    double py=gam.py();
	    double pz=gam.pz();
	    
	    double gammaE=sqrt(px*px+py*py+pz*pz);
	    v_gammaE.push_back(gammaE);
	  }
      }
    //        cout <<"num d0cands: "<< D0Candidates.size() <<" chargedDs: "<< chargedDCandidates.size()<< " star: "<< DStarCandidates.size()<<endl;

    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    ///got all candidates, now reconstruct Ds...
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    reconstructD0();
    vector<Particle*> tempParticles;
    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
      {
	double confLevel;
	//		cout <<"mass before d0 mit: "<< (*itD)->p().mag()<<endl;
///!	if(!doKmFit(*(*itD),  confLevel,0))
if(!doKmVtxFit2(*(*itD),  confLevel,0))
	  {
	    delete *itD;
	  }
	else
	  {
	    tempParticles.push_back(*itD);
	    //	    cout <<"mass  d0 mit: "<< (*itD)->p().mag()<<endl;
	  }

      }
    D0Candidates.clear();
    
    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
      {
	D0Candidates.push_back(*itD);
      }
    tempParticles.clear();

    reconstructChargedD();
    for(vector<Particle*>::iterator itD=chargedDCandidates.begin();itD!=chargedDCandidates.end();itD++)
      {
	double confLevel;
       if(!doKmVtxFit2(*(*itD),  confLevel,0))
///!	if(!doKmFit(*(*itD),  confLevel,0))
	  {
	    //	    cout <<" no good charged D " <<endl;
	    delete *itD;
	  }
	else
  {
	    //	    cout <<" good charged D " << endl;
	    tempParticles.push_back(*itD);
	  }
      }
    chargedDCandidates.clear();

    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
      {
	chargedDCandidates.push_back(*itD);
      }
    tempParticles.clear();

    reconstructDStar();
    for(vector<Particle*>::iterator itD=DStarCandidates.begin();itD!=DStarCandidates.end();itD++)
      {
	double confLevel;
	//	cout <<" d star before: "<<	(*itD)->p().mag() <<endl;
	//in principle mass-vertex constrained fit here (not just mass...)
		if(!doKmVtxFit2(*(*itD),  confLevel,0))
///!	if(!doKmFit(*(*itD),  confLevel,0))
  {
    //	    cout <<" no good dstar .. " <<endl;
	    delete *itD;
	  }
	else
	{
	    //	      cout <<"found good dstar .." <<endl;
	    tempParticles.push_back(*itD);
	    //	cout <<" d star after: "<<	(*itD)->p().mag() <<endl;
	  }

      }

    DStarCandidates.clear();
    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
      {
	DStarCandidates.push_back(*itD);
      }
    tempParticles.clear();

    //    cout <<"we have " <<D0Candidates.size()<<" D0 candidates.." <<endl;
    for(int i=0;i<D0Candidates.size();i++)
      {
	histoD0CandidateMass->Fill(D0Candidates[i]->mass());
      }

    if(m_mc)
      {
	//the noDRec checks if the D meson is actually in a decay that we can reconstruct...
	if(numNu==1 && decIds.size()==foundDecIds.size() && foundDPiPi && bestLogProb> -1000000 && !noDRec)
	  {
	    //    cout <<"heureka.." <<endl;
	    treeData.recBToDlNuPiPi=1;
	    int x=0;
	    //	  checkForDPiPi(x,true);	    
	  }
	if(numNu==1 && decIds.size()==foundDecIds.size() && foundSinglePionDecay && bestLogProb> -1000000 && !noDRec)
	  {
	    treeData.recBToDlNuPi=1;
	  }
      }

    //    cout <<"combining Dpil" <<endl;
    ///now put together D and pi's...
    //chargedPiCandidates has all the pions that are not part of the best b candidate
    //we just have to check if it is not part of the D children. Also how much other particles to we allow in the event?
    treeData.size=0;

    DType dtype=dtype_D0;
    treeData.recDDoubleStar=0;
    double bestMNu2=100000;
    int bestDIndex=-1;

    //do the best D separately for D,D* and  n=1,2
    float bestMDiff_D1=100000;
    float bestMDiff_D2=100000;
    float bestMDiff_DStar1=100000;
    float bestMDiff_DStar2=100000;

    int bestDIndex1=-1;
    int bestDStarIndex1=-1;
    int bestDIndex2=-1;
    int bestDStarIndex2=-1;

    bool mcDecaySignature=false;
    bool foundRecDecay=false;
    if(m_mc)
      {
	//if we would want to add the sig_numPion==0 case we should add this case w/o the foundDDoubleStar condition
	//which will obviously not be fulfilled
	mcDecaySignature=(sig_FoundDDoubleStar && sig_numLeptons==1 && sig_numKaons==0 && (sig_numPions==1 || sig_numPions==2) && sig_numPi0==0 && sig_numBaryons==0);
      }

    if(mcDecaySignature)
      {
	//      cout <<"found decay, D** id is: "<< dDoubleStarId <<endl;
      }

    treeData.bestBCharge=bestBcand.charge();

    /////
    treeData.pidCorrection=pidWeight;
    treeData.sPidCorrection=sqrt(sPidWeight)*pidWeight;
    treeData.CrossSectionLumiCorrection=XSectionWeight*lumiWeight;
    treeData.sCrossSectionLumiCorrection=sqrt(sXSectionWeight*sXSectionWeight/(XSectionWeight*XSectionWeight)+sLumiWeight*sLumiWeight/(lumiWeight*lumiWeight))*XSectionWeight*lumiWeight;
    //    cout <<"set sCrossLumi to : " << treeData.sCrossSectionLumiCorrection;
    //    cout <<" sX: "<< sXSectionWeight <<" lumi: "<< sLumiWeight << " weight: "<< treeData.CrossSectionLumiCorrection <<endl;

    treeData.sRelTrackingCorrection=sqrt(sTracking);
    treeData.numTracks=numTracks;



    //loop over found d mesons and check for each candidate if it fulfills the decay signature together with the rest of the event
    while(dtype!=dtype_end)
      {
	vector<Particle*>* dMesons;
	switch(dtype)
	  {
	  case dtype_D0:
	    dMesons=&D0Candidates;
	    break;
	  case dtype_DCharged:
	    dMesons=&chargedDCandidates;
	    break;
	  case dtype_DStar:
	    dMesons=&DStarCandidates;
	    break;

	  }
	dtype=(DType)((int)dtype+1);

	for(vector<Particle*>::iterator itD=dMesons->begin();itD!=dMesons->end();itD++)
	  {
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>((*itD)->userInfo());

	    /////---get the new field DDiff, DStarDiff, hypDMass1, hypDMass2
	    float DDiff=-1;
	    float DStarDiff=-1;
	    float hypDMass1=-1;
	    float hypDMass2=-1;

	    ////----

	    double dCharge=(*itD)->charge();
	    //	    cout <<"dcharge: " << dCharge <<endl;
	    vector<Particle*> localPiCandidates;
	    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	      {
		bool isChild=false;
		if(checkDoubleUse(*(*itD),*(*itP)))
		  isChild=true;

		for(int i=0;i<(*itD)->nChildren();i++)
		  {
		    Particle& child=(*itD)->child(i);
		    if(child.relation().isIdenticalWith((*itP)->relation()))
		      {
			isChild=true;
		      }
		    //need one more level for D*
		    for(int j=0;j<child.nChildren();j++)
		      {
			if(child.child(j).relation().isIdenticalWith((*itP)->relation()))
			  {
			    isChild=true;
			  }
		      }
		  }
		if(!isChild)
		  localPiCandidates.push_back(*itP);

	      }
	    int numOtherTracks=otherChargedTracks.size();
	    int numExtraKaons=0;
	    for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
	      {
		bool isChild=false;
		if(checkDoubleUse(*(*itD),*(*itK)))
		  isChild=true;
		for(int i=0;i<(*itD)->nChildren();i++)
		  {
		    if((*itD)->child(i).relation().isIdenticalWith((*itK)->relation()))
		      {
			isChild=true;
		      }
		    //need one more level for D*?
		    Particle& child=(*itD)->child(i);
		    for(int j=0;j<child.nChildren();j++)
		      {
			if(child.child(j).relation().isIdenticalWith((*itK)->relation()))
			  {
			    isChild=true;
			  }
		      }
		  }
		if(!isChild)
		  numExtraKaons++;
	      }

	    double pionCharge=0;
	    double leptonCharge=0;
	    for(int i=0;i<localPiCandidates.size();i++)
	      {
		pionCharge+=localPiCandidates[i]->charge();
	      }
	    if(leptonCandidates.size()==1)
	      {
		Gen_hepevt hepEvt=get_hepevt(leptonCandidates[0]->mdstCharged());
		//		cout << " we should have : "<<leptonCandidates[0]->pType().name() <<" and have: "<<  hepEvt.idhep() <<endl;
		leptonCharge=leptonCandidates[0]->charge();
	      }
	    int numPions=localPiCandidates.size();
	    //	    cout <<"having " << numPions <<endl;

	    HepLorentzVector pionMom(0.0,0.0,0.0,0.0);
	    for(int i=0;i<numPions;i++)
	      {
		if(i==0)
		  {
		    pionMom=localPiCandidates[i]->p();
		  }
		else
		  {
		    pionMom+=localPiCandidates[i]->p();
		  }
	      }

	    ///in the ==1 case, hypDMass doesn't make sense, since it is equal to mDn anyways...
	    if(numPions==2)
	      {
		Particle& D=*(*itD);
		hypDMass1=(D.p()+localPiCandidates[0]->p()).mag();
		hypDMass2=(D.p()+localPiCandidates[1]->p()).mag();
	      }

	    ///


	    //	    cout <<"pionCharge: " << pionCharge <<" lepton Charge: "<< leptonCharge <<" sum: "<< pionCharge+leptonCharge+dCharge <<endl;

	    if(PRINT)
	      cout <<"numOtherTracks: "<< numOtherTracks <<" numExtraKaons: "<< numExtraKaons <<" numPions: "<< numPions << "lepton Candidates: "<< leptonCandidates.size() << " system charge: "<< pionCharge+leptonCharge+dCharge <<" best B charge: "<<  bestBcand.charge() <<endl;


	    //	    bool recDecaySignature=numOtherTracks==0&&numExtraKaons==0&&numPions>=1&&numPions<=2 &&leptonCandidates.size()==1 && fabs(pionCharge+leptonCharge+dCharge)<=1.0 && bestBcand.charge()==((-1)*(pionCharge+leptonCharge+dCharge));
	    ///	    bool recDecaySignature=numOtherTracks==0&&numExtraKaons==0&&numPions>=1&&numPions<=2 &&leptonCandidates.size()==1 && fabs(pionCharge+leptonCharge+dCharge)<=1.0;
	    //adding zero pions
	    bool recDecaySignature=numOtherTracks==0&&numExtraKaons==0 && numPions>=0 && numPions<=2 &&leptonCandidates.size()==1 && fabs(pionCharge+leptonCharge+dCharge)<=1.0;
	    //	    cout <<"recDecaySignature: " << recDecaySignature<<endl;
	    //	    cout <<"numPions: "<<numPions <<" numExtraKaons: "<< numExtraKaons <<" other: "<< numOtherTracks <<" leptons:" << leptonCandidates.size() <<" charges: "<< pionCharge+leptonCharge+dCharge <<endl;

	    //	    cout <<" using pidWeight: "<< pidWeight <<" B Br weight: " << B_BR_CorrectionFactor << " D BR factor: "<< D_BR_CorrectionFactor<<endl;
	    treeData.systemCharge[treeData.size]=pionCharge+leptonCharge+dCharge;
	    treeData.leptonCharge[treeData.size]=leptonCharge;
	    treeData.dCharge[treeData.size]=dCharge;


	    if(recDecaySignature && (!(excludeSignal && sig_FoundDDoubleStar)))
	      {
		foundRecDecay=true;
		if(PRINT)
		  cout <<"all conditions satisfied ... " << endl;
		Particle& D=*(*itD);
		treeData.recDDoubleStar=1;

		//get this by computing B_sl direction from B_tag. B tag is saved in 'bestBCand' particle, so 
		//I assume that the B_sl momentum is just the mirrored momentum;
		float mNu2=0;		
		float U=0;
		HepLorentzVector bMomentum=bestBcand.p();
		//		cout <<"signal bMomentum vect; "<< signalBMomentum.vect()<<endl;
		//doesn't make sense to just flip the sign of p, since we are not in the CMS
		//	    	    cout <<"best b m: "<< bMomentum.mag() <<endl;

		HepLorentzVector Pxl;
		if(leptonCandidates.size()>0)
		  Pxl=D.p()+leptonCandidates[0]->p()+pionMom;
		else
		  Pxl=D.p()+pionMom;
		float mDnPi=(D.p()+pionMom).mag();
		HepLorentzVector pMDnPi=(D.p()+pionMom);
		if(PRINT)
		  cout <<"mDnPi mass: "<< mDnPi << "mom: "<<pMDnPi.rho()<<"  ("<<pMDnPi.px()<<", " << pMDnPi.py()<< ", " << pMDnPi.pz() <<", "<< pMDnPi.t()<<")" <<endl;
		//	    if(mDnPi<2.0) 
	      
		//	      continue;
		treeData.bestD[treeData.size]=0;
		treeData.numRecPions[treeData.size]=numPions;
		treeData.mDnPi[treeData.size]=(D.p()+pionMom).mag();
		treeData.dDecay[treeData.size]=pinf.decayChannel;
		treeData.dType[treeData.size]=pinf.type;
		//recoil is computed wrt to signal b, so we have to flip the momentum. However, we need the 'real' D from MC (so D** if it is a D** event). Note 
		//that we are only looking at Dlnu events (with D =D(*) or D** 
		//		treeData.w[treeData.size]=


		//FF weight is dependent on q2, and p_l* (momentum of the lepton in the B CMS system) for the D(*) l nu
		// for D** l nu in w and cos(theta_l)

		//		HepLorentzVector p_W = p_l + p_nu;
		//		p_l.boost(-(p_W.boostVector()));
		//		p_D.boost(-(p_W.boostVector()));
		//		float cosTheta = cos( p_l.p3().angle( p_D.p3() ) );

		//		cout <<" b t comp:  " << signalBMomentum.t() <<" D time; "<< D.p().t() <<" D vect: "<< D.p().vect()<<endl;
		//		cout <<" w is: "<< D.p().dot(signalBMomentum)<<" with tag: "<< D.p().dot(bMomentum)<<endl;
		

		//q2 =(p_l+p_\nu)^2
		//		cout <<"or  w is: "<< D.p().t()*signalBMomentum.t()-D.p().vect().dot(signalBMomentum.vect())<<" with tag: "<< D.p().dot(bMomentum)<<endl;
				
		if(leptonCandidates.size()>0)
		  {
		    treeData.leptonMom[treeData.size]=leptonCandidates[0]->p().vect().mag();
		    treeData.leptonTheta[treeData.size]=leptonCandidates[0]->p().vect().theta();
		    treeData.leptonPhi[treeData.size]=leptonCandidates[0]->p().vect().phi();

		  }
		else
		  {
		    treeData.leptonMom[treeData.size]=-999;
		    treeData.leptonTheta[treeData.size]=-999;
		    treeData.leptonPhi[treeData.size]=-999;
		  }
		if(numPions>0)
		  {
		    treeData.pi1Mom[treeData.size]=localPiCandidates[0]->p().vect().mag();
		    treeData.pi1Theta[treeData.size]=localPiCandidates[0]->p().vect().theta();
		    treeData.pi1Phi[treeData.size]=localPiCandidates[0]->p().vect().phi();
		  }
		else
		  {
		    treeData.pi1Mom[treeData.size]=-999;
		    treeData.pi1Theta[treeData.size]=-999;
		    treeData.pi1Phi[treeData.size]=-999;
		  }
		if(numPions>1)
		  {
		    treeData.pi2Mom[treeData.size]=localPiCandidates[1]->p().vect().mag();
		    treeData.pi2Theta[treeData.size]=localPiCandidates[1]->p().vect().theta();
		    treeData.pi2Phi[treeData.size]=localPiCandidates[1]->p().vect().phi();
		  }
		else
		  {
		    treeData.pi2Mom[treeData.size]=-999;
		    treeData.pi2Theta[treeData.size]=-999;
		    treeData.pi2Phi[treeData.size]=-999;
		  }

		//see babar paper
		HepLorentzVector p4Miss=(kinematics::pBeam-bMomentum-Pxl);
		float E_miss=p4Miss.e();
		U=E_miss-p4Miss.vect().mag();

		mNu2=(kinematics::pBeam-bMomentum-Pxl).mag2(); 
		if(mNu2<bestMNu2)
		  {
		    bestMNu2=mNu2;
		    bestDIndex=treeData.size;
		  }

		///best selection based on D Diff;
		//because it was already incremented
		DType tmpDtype=(DType)((int)dtype-1);
		if(numPions==1 && tmpDtype!=dtype_DStar)
		  {
		    if(fabs(pinf.pdgDiff)<bestMDiff_D1)
		      {
			bestMDiff_D1=fabs(pinf.pdgDiff);
			bestDIndex1=treeData.size;
		      }
		  }
		if(numPions==1 && tmpDtype==dtype_DStar)
		  {
		    if(fabs(pinf.pdgDiff)<bestMDiff_DStar1)
		      {
			bestMDiff_DStar1=fabs(pinf.pdgDiff);
			bestDStarIndex1=treeData.size;
		      }

		  }
		if(numPions==2 && tmpDtype!=dtype_DStar)
		  {
		    if(fabs(pinf.pdgDiff)<bestMDiff_D2)
		      {
			bestMDiff_D2=fabs(pinf.pdgDiff);
			bestDIndex2=treeData.size;
		      }

		  }
		if(numPions==2 && tmpDtype==dtype_DStar)
		  {
		    if(fabs(pinf.pdgDiff)<bestMDiff_DStar2)
		      {
			bestMDiff_DStar2=fabs(pinf.pdgDiff);
			bestDStarIndex2=treeData.size;
		      }

		  }

	      
		////

		//		sprintf(bufferAll,"%s && (!foundAnyDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2 || sig_DLNu || sig_DPiLNu|| sig	\
		//_DPiPiLNu || sig_DStarLNu || sig_DStarPiLNu|| sig_DStarPiPiLNu) && !sig_DLNu && !sig_DPiLNu && !sig_DPiPiLNu && !sig_DStarLNu && !sig_DStarPiLNu && !sig_DStarPiPiLNu)",buffer);

		////-------
		////investigate the 'otherBB' shape in the pion==1 case which is a bit too signal like for my taste...
		if(numPions==1 && fabs(mNu2)< 0.2 && leptonCandidates.size()==1 && recDecaySignature)
		  {
		    if(!sig_FoundDDoubleStar|| sig_numKaons!=0 || sig_numPi0!=0 || sig_numBaryons!=0 || sig_numLeptons!=1 || sig_numPions > 2)
		      {
			if(!sigDLNu && !sigDPiLNu && !sigDPiPiLNu && !sigDStarLNu && !sigDStarPiLNu && !sigDStarPiPiLNu)
			  {
			    cout<<" other BB background evt nr "<< evtNr << " run: "<< runNr  << " tag B momentum: " << bMomentum.px() <<", "<< bMomentum.py()<<", " << bMomentum.pz()<<endl;
			    cout <<" numLept: "<< sig_numLeptons<< " kaons: "<< sig_numKaons <<" pions: "<< sig_numPions <<" pi0s: "<< sig_numPi0 <<" num baryos: "<< sig_numBaryons;
			    cout <<" numD: "<< sig_numD <<" DStar " << sig_numDStar <<" num rhos: "<< sig_numRhos <<endl;

			    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
			    int numB=0;
			    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
			      {
				if((fabs(gen_it->idhep())>=500 && fabs(gen_it->idhep()<600)) || (fabs(gen_it->idhep())>=10500 && fabs(gen_it->idhep())<10600))
				  {
				    numB++;
				    recursivePrint(*gen_it,"");
				  }
			      
			      }
			    //			  if(numB<2)
			    //			    {
			    //			      cout <<"no B found ?, event:  " <<  endl;
			    //			      for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
			    //				{
			    //				  Particle p(*gen_it);
			    //				  int lund=fabs(gen_it->idhep());
			    //				  if(lund== 100423 || lund ==100421 ||lund==PY_DStar_2S || lund==100411 || lund==100413 ){
			    //				    cout <<"nd: " << lund << endl;
			    //				  }
			    //				  else{
			    //				    //for some reason some ids make the particle class crash when asked for the name...
			    //				    if(lund<10000)
			    //				      cout <<p.pType().name()<<" (" << gen_it->idhep()<<"), "<<" |p|: " <<p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;
			    //				    else
			    //				      cout  <<" (" << gen_it->idhep()<<"), "<<" |p|: " <<p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;
			    //				  }
			    //
			    //				}
			    //
			    //			    }
			  
			  }
		      }
		    
		  }
	      
		////-------
		////-------
		////-------
	      

		if(PRINT)
		  {

		    if((sigDStarLNu || sigDStarPiLNu || sigDStarPiPiLNu )&& sig_FoundDDoubleStar)
		      {
			//			cout <<"found dlnu and ddouble star again!!!!!" <<endl;
			//			cout <<"mNu2: " << mNu2 <<endl;
		      }
		    if(fabs(mNu2) < 1.0 && numPions==2)
		      {

			if(sig_FoundDDoubleStar)
			  {
			    //			    cout <<"smallMNuTwoPionDDoubleStar" <<endl;
			    foundDDStarFlag=true;
			  }
			if(!sigDStarLNu && !sigDStarPiLNu && !sigDStarPiPiLNu && !sigDLNu && !sigDPiLNu && !sigDPiPiLNu)
			  {
			    if((!sig_FoundDDoubleStar || sig_numPions>2 || sig_numKaons!=0 || sig_numBaryons!=0 || sig_numPi0!=0 || sig_numLeptons!=1))
			      {
				//				cout <<"background event .." <<endl;
				bgFlag=true;
			      }
			  }
			//			cout <<" small mNu2 (" << mNu2 <<" ),  sigDLNu: "<< sigDLNu << " dpilnu: "<< sigDPiLNu <<" dpipilnu: "<< sigDPiPiLNu<<endl;
			//			cout <<"Dstarlnu: " << sigDStarLNu <<" dstarpilnu: "<< sigDStarPiLNu <<" dstarpipilnu: " << sigDStarPiPiLNu <<endl;

			//			cout <<"found any ddouble star? : "<< sig_FoundDDoubleStar <<" num pions: "<< sig_numPions <<" num kaons: ";
			//			cout <<sig_numKaons <<" num pi0:"<< sig_numPi0 <<" num baryon: " << sig_numBaryons <<endl;
			if(!sigDLNu && !sigDPiLNu && !sigDPiPiLNu && !sigDStarLNu && !sigDStarPiLNu && !sigDStarPiPiLNu)
			  {
			    if(!sig_FoundDDoubleStar || (sig_numPions> 2 || sig_numKaons>0 || sig_numBaryons>0 || sig_numPi0>0|| sig_numLeptons!=1)) {
			      cout <<"found strange strange background " <<endl;
			    }
			  }

			Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
			for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
			  {
			    if((fabs(gen_it->idhep())>=500 && fabs(gen_it->idhep()<600)) || (fabs(gen_it->idhep())>=10500 && fabs(gen_it->idhep())<10600))
			      {
				recursivePrint(*gen_it,"");
			      }

			  }
		      }

		    if(fabs(mNu2)<0.02 && !foundDPiPi  && !foundSinglePionDecay)
		      {
			cout <<"small nu mass but nothing found... "<< mNu2 << endl;
		      }
		    if(fabs(mNu2)<0.05 && !sig_FoundDDoubleStar)
		      {
			cout <<"small nu but did not find any" <<endl;
		      }
		  }///end of PRINT
		////////---------
		///////////////////

		if(tmpDtype==dtype_DStar)
		  treeData.DStarDiff[treeData.size]=pinf.pdgDiff;
		else
		  treeData.DDiff[treeData.size]=pinf.pdgDiff;

		treeData.recDType[treeData.size]=dtype;
		//		cout <<"mNu: " << mNu <<endl;
		//		cout <<"lepton mass: "<< leptonCandidates[0]->p().mag()<<endl;
		//				  cout <<"b mass: "<< bMomentum.mag() <<endl;
		//				  cout <<"beam mass: "<< kinematics::pBeam.mag()<<endl;
		//		mNu=(kinematics::pBeam-bMomentum-Pxl).mag();
		//		cout <<"found all we want... with d: " << dtype<<" mNu: "<< mNu <<" m_pxl: "<< Pxl<<endl;
		//	    cout <<" beam; "<< kinematics::pBeam <<endl;
		//	    cout <<" bMomentum: "<< bMomentum <<endl;
		//	    cout <<" Pxl: "<< Pxl <<endl;
		//	    cout <<"mnu: "<< mNu <<endl;
		//	    mNu=(kinematics::pBeam-bMomentum-sigBMom).mag();
		//		float mNu2=(bMomentum-Pxl).mag();
		//	    cout <<"mNu: "<< mNu<< " mNu2: "<< mNu2 <<endl;
		//	    mNu=(bMomentum-Pxl).mag();
		//mNu=(sigBMom-Pxl).m();
		treeData.U[treeData.size]=U;
		treeData.mNu2[treeData.size]=mNu2;
		treeData.mXl[treeData.size]=Pxl.mag();
		treeData.mB[treeData.size]=bMomentum.mag();
		treeData.hypDMass1[treeData.size]=hypDMass1;
		treeData.hypDMass2[treeData.size]=hypDMass2;


		treeData.size++;
		//		cout <<"got mNu: " << mNu <<endl;

	      
		//////////
		//		if(treeData.recBToDlNuPi==1 || treeData.recBToDlNuPiPi==1)
		//		if(treeData.foundDPiPi==1 || foundSinglePionDecay)
		if(fabs(mNu2)<0.1)
		  {
		    if(PRINT)
		      {
			if(numPions==2&& sig_FoundDDoubleStar)
			  {
			    //			    cout <<" found DDstar and two pions" <<endl;
			  }
		      }
		    if(PRINT)
		      cout <<" -------"<<endl<<endl;
		    if(PRINT)
		      {
			if(found_2SD!=sig_d_2S)
			  cout <<" 2SD difference" << found_2SD <<" sig: " << sig_d_2S<<endl;
			if(found_2SD_Star!=sig_dStar_2S)
			  cout <<" 2SD Star difference" << found_2SD_Star <<" sig: "<< sig_dStar_2S<<endl;
			cout <<"decay signature: "<< endl;
			cout <<"found DDouble STAR: "<< sig_FoundDDoubleStar <<" numPions: "<< sig_numPions << " num Kaons: "<< sig_numKaons <<" num leptons:" << sig_numLeptons << " num pi0: "<< sig_numPi0 << " num Baryons: " << sig_numBaryons <<endl;
		      }
		    if(PRINT)
		      cout <<"generated data: "<< endl;

		    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
		    HepLorentzVector b1,b2;
		    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
		      {
			if(gen_it->idhep()!=911)
			  {
			    //			    Particle p(*gen_it);
			    //			    cout <<p.pType().name()<<" mom: "<< p.p().rho()<<endl;
			  }
			//B meson
			if((fabs(gen_it->idhep())>=500 && fabs(gen_it->idhep()<600)) || (fabs(gen_it->idhep())>=10500 && fabs(gen_it->idhep())<10600))
			  {
			    Particle p(*gen_it);
			    if(gen_it->idhep()>0)
			      {
				b1=p.p();
			      }
			    else
			      {
				b2=p.p();
			      }
			    if(PRINT)
			      recursivePrint(*gen_it,"");
			  }
		      }
		    bool tmp;
		    bool tmp2;
		    int iTmp;
	    
	    
		    if(PRINT)
		      {
		
		      }
	    
		    if(PRINT)
		      {
			cout <<"check print: " <<endl;
			checkForDPiPi(iTmp,tmp, tmp2,true);
		      }
		    if(PRINT)
		      cout <<"printing out b decay products: "<< endl;
		    set<int> t1, t2,t3;
		    if(PRINT)
		      getDecayIds(bestBcand,t1,t2,t3,true);
		    if(PRINT)
		      cout <<"D decay products: "<< endl;
		    if(PRINT)
		      getDecayIds(*(*itD),t1,t2,t3,true);

		    HepLorentzVector bSum=b1+b2;
		    if(PRINT)
		      cout <<"done generated " << endl<<endl;
		    if(PRINT)
		      cout << " sum of b vectors: " << bSum.px() <<" " << bSum.py()<<" " << bSum.pz() << " " << bSum.t()<<endl;



		    HepLorentzVector recP;
		    getMassFromDecayParticles(bestBcand,recP);
		    if(PRINT)
		      cout <<"b mass: " << bestBcand.p().mag()<< " or: " << recP.mag()<<endl; 
		    if(PRINT)
		      cout <<"bMomentum: " <<bMomentum.rho() <<" rec: " << recP.rho()<<endl;
		    if(PRINT)
		      cout <<"b 4-vector: "<< bMomentum.px() <<" " << bMomentum.py() << " " << bMomentum.pz() << " " << bMomentum.t() <<endl;
		    if(PRINT)
		      cout <<"rec b 4-vector: "<< recP.px() <<" " << recP.py() << " " << recP.pz() << " " << recP.t() <<endl;
		    if(PRINT)
		      cout <<"deltaE: "<< bestDeltaE <<" bestTagM: "<< bestTagM <<endl;
		    if(PRINT)
		      cout <<"best b cand 4-vector: "<< bestBcand.p().px() <<" " << bestBcand.p().py() << " " << bestBcand.p().pz() << " " << bestBcand.p().t() <<endl;
		    
		    
		    //		mNu=(kinematics::pBeam-bMomentum-Pxl).mag();

		    HepLorentzVector recD;
		    getMassFromDecayParticles(*(*itD),recD);
		    if(PRINT)
		      cout <<"D mass from decay: " << recD.mag() <<" p mass: " <<  (*itD)->pType().mass() <<" or : " << (*itD)->p().mag()<<endl;
		    if(PRINT)
		      cout <<"D mom: " << " ("<<recD.px()<<", " << recD.py()<< ", " << recD.pz() <<", "<< recD.t()<<")" <<endl;
		    if(PRINT)
		      cout <<"or: "<< " ("<<(*itD)->p().px()<<", " << (*itD)->p().py()<< ", " << (*itD)->p().pz() <<", "<< (*itD)->p().t()<<")" <<endl;
		    if(PRINT)		    
		      cout<<"pion ids and momenta: "<< endl;
		    for(int i=0;i<numPions;i++)
		      {
			if(PRINT)
			  cout <<localPiCandidates[i]->mdstCharged().get_ID() <<", " << localPiCandidates[i]->p().rho()<<endl;
		      }
		    if(PRINT)
		      cout <<"lepton candidates: "<< endl;
		    if(PRINT&& leptonCandidates.size()>0)
		      cout <<leptonCandidates[0]->mdstCharged().get_ID() <<", " << leptonCandidates[0]->p().rho()<<endl;
		    if(PRINT)
		      cout <<"pxl 4-vector: "<< Pxl.px() <<" " << Pxl.py() << " " << Pxl.pz() << " " << Pxl.t() <<endl;
		    if(PRINT)
		      cout <<"   D 4-vector: "<< recD.px() <<" " << recD.py() << " " << recD.pz() << " " << recD.t() <<endl;
		    if(PRINT&& leptonCandidates.size()>0)
		      cout <<"   l 4-vector: "<< leptonCandidates[0]->p().px() <<" " << leptonCandidates[0]->p().py() << " " << leptonCandidates[0]->p().pz() << " " << leptonCandidates[0]->p().t() <<endl;
		    if(PRINT)
		      cout <<"pion 4-vector: "<< pionMom.px() <<" " << pionMom.py() << " " << pionMom.pz() << " " << pionMom.t() <<endl;
		    if(PRINT && numPions>0)
		      cout <<" pi1 4-vector: "<< localPiCandidates[0]->p().px() <<" " << localPiCandidates[0]->p().py() << " " << localPiCandidates[0]->p().pz() << " " << localPiCandidates[0]->p().t() <<endl;
		    if(numPions>1){
		      if(PRINT)
			cout <<" pi2 4-vector: "<< localPiCandidates[1]->p().px() <<" " << localPiCandidates[1]->p().py() << " " << localPiCandidates[1]->p().pz() << " " << localPiCandidates[1]->p().t() <<endl;

		    }
		    if(PRINT)
		      cout <<"beam 4-vector: "<< kinematics::pBeam.px() <<" " << kinematics::pBeam.py() << " " << kinematics::pBeam.pz() << " " << kinematics::pBeam.t() <<endl;
		    HepLorentzVector mNu4=kinematics::pBeam-bMomentum-Pxl;
		    if(PRINT)
		      cout <<" mNu 4-vector: "<< mNu4.px() <<" " << mNu4.py() << " " << mNu4.pz() << " " << mNu4.t() <<endl;
		    if(PRINT)
		      cout <<"evt: "<< evtNr <<" runNr: " << runNr <<" mNu2: " << mNu2<<endl;
		  } //if(mNu2<0.1)
		


		//////////////

		
	      } ///if recDecaySignature

	  }/// for all Ds
      } ///for all DTypes
    
    /////--->Don't select best D based on mNu2 (bias), instead using best Ds (separately for n12, D and DStar
    bool singlePionEvent=false;
    bool twoPionEvent=false;
    //priority to DStar
    if(bestDStarIndex1>=0 && bestDStarIndex1<1000)
      {
	treeData.bestD[bestDStarIndex1]=1;
	singlePionEvent=true;
      }
    else
      {
	if(bestDIndex1>=0 && bestDIndex1<1000)
	  {
	    treeData.bestD[bestDIndex1]=1;
	    singlePionEvent=true;
	  }
      }
    //do n=1 and n=2 separately (we can tell by the numPions field)
    if(bestDStarIndex2>=0 && bestDStarIndex2<1000)
      {
	treeData.bestD[bestDStarIndex2]=1;
	twoPionEvent=true;
      }
    else
      {
	if(bestDIndex2>=0 && bestDIndex2<1000)
	  {
	    treeData.bestD[bestDIndex2]=1;
	    twoPionEvent=true;
	  }
      }

    if(singlePionEvent && twoPionEvent)
      treeData.overlapEvent=true;
    else
      treeData.overlapEvent=false;



    //    cout <<"done combining.." <<endl;
    //    if(bestDIndex>=0 && bestDIndex < 1000)
    //      {
    //	treeData.bestD[bestDIndex]=1;
    //      }
////    if(treeData.size>1)
////      {
////	cout<<" we have several D candidates..with " <<endl;
////	for(int i=0;i<treeData.size;i++)
////	  {
////	    	    cout << treeData.numRecPions[i] << ", ";
////	  }
////	cout <<" rec pions " <<endl;
////	printUse();
////      }
    treeData.recDecaySignature=foundRecDecay;
    treeData.mcDecaySignature=mcDecaySignature;
    if((foundRecDecay || sigResDLNu || sigResDPiLNu || sigResDPiPiLNu || sigResDStarLNu || sigResDStarPiLNu || sigResDStarPiPiLNu)  && (!(excludeSignal && sig_FoundDDoubleStar)))
      {
	//	cout <<"saving tree, bgFlag: "<< bgFlag <<endl;
	//	cout <<"saving tree, DD flag: "<<  foundDDStarFlag << endl;
	//	cout <<"saving have pidweight: "<< treeData.pidCorrection<<endl;



    float tempWeight=1.0;
    for(int i1=0;i1<treeData.numKsCorr;i1++)
      {
	tempWeight*=treeData.KsCorrection[i1];
      }
    for(int i1=0;i1<treeData.numChargedCorr;i1++)
      {
	tempWeight*=treeData.ChargedCorrection[i1];
      }
    for(int i1=0;i1<treeData.numPi0Corr;i1++)
      {
	tempWeight*=treeData.pi0Correction[i1];
      }
        cout <<" tempWeight: "<< tempWeight <<" pidWeight: "<< pidWeight <<endl;
    //

	saveTree();
	if(dDoubleStarId==100423)
	  {
	    //	    cout <<"saving found the 2S with mass " << mcMassDDStar<<endl;
	  }

	//	cout <<"indeed foundRec " <<endl;
      }
    //so as not to save twice
    if((mcDecaySignature&& !foundRecDecay) && (!(excludeSignal && sig_FoundDDoubleStar)))
      {
	//	cout <<"saving tree, bgFlag: "<< bgFlag <<endl;
	//	cout <<"saving tree, DD flag: "<<  foundDDStarFlag << endl;
	//		cout <<"mc decay but not foundRec " <<endl;
	
	saveTree();
	if(dDoubleStarId==100423)
	  {
	    //	    cout <<"saving found the 2S with mass " << mcMassDDStar<<endl;
	  }

      }
    //    if(!mcDecaySignature && !foundRecDecay)
    //      cout <<"haven't found anything " <<endl;
    exitEvent();
  
  }

                                                                                                                                                                                                     
  void bToDDoubleStar::getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz, float m_mass)
  {
    //    cout <<"original p: " << (*chr_it).p(0) <<" py: "<< (*chr_it).p(1) << " pz: "<< (*chr_it).p(2) <<endl;
    Mdst_trk& mdsttrk = chr_it->trk();
    Mdst_trk_fit &mdsttrkfit=mdsttrk.mhyp(masshyp);
    HepPoint3D pivot(mdsttrkfit.pivot_x(),mdsttrkfit.pivot_y(),mdsttrkfit.pivot_z());

    HepVector a( 5, 0 );
    a[0] = mdsttrkfit.helix( 0 ); // helix parameters defined at the pivot
    a[1] = mdsttrkfit.helix( 1 );
    a[2] = mdsttrkfit.helix( 2 );
    a[3] = mdsttrkfit.helix( 3 );
    a[4] = mdsttrkfit.helix( 4 );

    HepSymMatrix Ea( 5, 0 );
    Ea[0][0] = mdsttrkfit.error( 0 );
    Ea[1][0] = mdsttrkfit.error( 1 );
    Ea[1][1] = mdsttrkfit.error( 2 );
    Ea[2][0] = mdsttrkfit.error( 3 );
    Ea[2][1] = mdsttrkfit.error( 4 );
    Ea[2][2] = mdsttrkfit.error( 5 );
    Ea[3][0] = mdsttrkfit.error( 6 );
    Ea[3][1] = mdsttrkfit.error( 7 );
    Ea[3][2] = mdsttrkfit.error( 8 );
    Ea[3][3] = mdsttrkfit.error( 9 );
    Ea[4][0] = mdsttrkfit.error( 10 );
    Ea[4][1] = mdsttrkfit.error( 11 );
    Ea[4][2] = mdsttrkfit.error( 12 );
    Ea[4][3] = mdsttrkfit.error( 13 );
    Ea[4][4] = mdsttrkfit.error( 14 );


    Helix helix( pivot, a, Ea );
    helix.pivot( IpProfile::position(1));
    //    cout <<"helix momentum: "<< helix.momentum()<<endl;
    refitPx=helix.momentum().x();
    refitPy=helix.momentum().y();
    refitPz=helix.momentum().z();
    //    cout <<"refit px: "<< refitPx <<" refitPy: "<< refitPy <<" refitPz: "<< refitPz <<endl;
    //    HepLorentzVector boostedVec(helix.momentum(),sqrt(helix.momentum().mag2()+m_mass*m_mass));
    dr  = helix.dr();
    dz  = helix.dz();
    //    HepLorentzVector lzMom(refitPx,refitPy,refitPz,sqrt(helix.momentum().mag2()+m_mass*m_mass));
    //    return Momentum(lzMom);

  }

  void bToDDoubleStar::exitEvent()
  {


#ifdef SAVE_HISTOS
    // saveHistos(allParticlesBoosted, allParticlesNonBoosted);
#endif
    cleanUp();
  }
  // begin_run function

  void bToDDoubleStar::end_run(BelleEvent* evptr, int* status)
  {
    std::cout << "bToDDoubleStar's end_run function" << std::endl;
    //
  }

  void bToDDoubleStar::saveTree()
  {

    //    doBRCorrections();
    //    doFFCorrections();





    pTreeSaver->saveData(&this->treeData);

  }
  void bToDDoubleStar::saveHistos( vector<Hep3Vector>& v_allParticlesBoosted, vector<Hep3Vector>& v_allParticlesNonBoosted)
  {


  }





  void bToDDoubleStar::cleanUp()
  {

    //    cout <<"cleaning up..." <<endl;
    //    cout <<" d0 size"<<D0Candidates.size()<<endl;
    for(int i=0;i<D0Candidates.size();i++){
      delete D0Candidates[i];
    }
    D0Candidates.clear();

    //    cout <<" done D0"<<D0Candidates.size()<<endl;

    //    cout <<"  DStar size " <<DStarCandidates.size()<< endl;
    for(int i=0;i<DStarCandidates.size();i++){
      delete DStarCandidates[i];
    }
    DStarCandidates.clear();
    //    cout <<" done DStar " <<DStarCandidates.size()<< endl;

    //    cout <<" charged D "<<chargedKCandidates.size()<<endl;
    for(int i=0;i<chargedDCandidates.size();i++){

      delete chargedDCandidates[i];
    }
    chargedDCandidates.clear();
    //    cout <<" done charged D "<<chargedDCandidates.size()<<endl;
    //    cout <<" charged K "<<chargedKCandidates.size()<<endl;
    for(int i=0;i<chargedKCandidates.size();i++){
      delete chargedKCandidates[i];
    }
    chargedKCandidates.clear();
    //    cout <<" done charged K"<< chargedKCandidates.size()<<endl;


    //    cout <<"pi0s: "<< pi0Candidates.size()<<endl;
    for(int i=0;i<pi0Candidates.size();i++){
      delete pi0Candidates[i];
    }
    pi0Candidates.clear();
    //    cout <<" done pi0"<<endl;

    for(int i=0;i<leptonCandidates.size();i++){
      delete leptonCandidates[i];
    }
    leptonCandidates.clear();
    //    cout <<" done leptons"<<endl;

    //    cout <<"Ks candidates: " << KsCandidates.size()<<endl;
    for(int i=0;i<KsCandidates.size();i++){
      delete KsCandidates[i];
    }
    KsCandidates.clear();


    //    cout <<" done Ks"<<endl;
    for(int i=0;i<chargedPiCandidates.size();i++){
      delete chargedPiCandidates[i];
    }
    chargedPiCandidates.clear();
    //    cout <<" done charged Pi"<<endl;
    for(int i=0;i<otherChargedTracks.size();i++){
      delete otherChargedTracks[i];
    }
    otherChargedTracks.clear();
    //    cout <<" charged tracks"<<endl;


    v_vertexR.clear();
    v_vertexZ.clear();
    v_pi0GammaE.clear();
    v_gammaE.clear();
    v_asyms.clear();
    //    cout <<"return ..." <<endl;

  }

  // begin_run function
  void bToDDoubleStar::term()
  {
    cout <<"term..." <<endl;
    mesonMassTree->Write();
    histoRecDSpect->Write();
    histoRecD0Spect->Write();
    histoRecDStarSpect->Write();
    histoRecDStarSpectToD0Pi->Write();
    histoRecDStarSpectToD0Pi0->Write();
    histoRecDStarSpectToDPi0->Write();
    histoRecDStarSpectToDPi->Write();

    histoD0CandidateMass->Write();
    histoKs->Write();
    histoChargedBTag_dE->Write();
    histoB0Tag_dE->Write();

    histoChargedBTag_M->Write();
    histoB0Tag_M->Write();

    histoPi0SlowMom->Write();
    histoDStar->Write();
    histoD0Spect->Write();

    //    cout <<"writing file.." <<endl;
    m_file->Write();
    //    cout <<"closing file.." <<endl;
    m_file->Close();
    //    cout <<"on to histo " <<endl;
  }

  void bToDDoubleStar::findDStar(vector<Hep3Vector>& allPB, vector<int>& allPB_Class, vector<int>& allPB_Charge)
  {
    const float D0mass=1.865;
    const float DStarMass=2.010;
    const float pionMass=0.139570;
    const float kaonMass=0.493667;

    kinematics::D0Tag=0;
    kinematics::DStarTag=0;

    if(allPB.size()!=allPB_Class.size() || allPB_Class.size()!=allPB_Charge.size())
      {
	//	cout<<"findstar size problem " << endl;
	exit(1);
      }
    for(int i=0;i<allPB.size();i++)
      {
	//find kaon- (should also work with k+...)
	if(allPB_Class[i]==3 && allPB_Charge[i]<0)
	  { 
	    if(allPB[i].mag()<0.1)  
	      continue;
	    double m_mass=kaonMass;
	    double EK=sqrt(m_mass*m_mass+allPB[i].mag2());

	    //find pion to form D0 (no vertex fit yet...)
	    for(int j=0;j<allPB.size();j++)
	      {
		if(allPB[j].mag()<0.1)
		  continue;

		//found pion, see if we can form D0 mass...
		if(allPB_Class[j]==2 && allPB_Charge[j]>0)
		  {
		    m_mass=pionMass;
		    float EP1=sqrt(m_mass*m_mass+allPB[j].mag2());
		    float d0candidateMass=sqrt((EK+EP1)*(EK+EP1)-(allPB[i]+allPB[j]).mag2());
		    histoD0Spect->Fill(d0candidateMass);
		    if(fabs(d0candidateMass-D0mass)<0.1)
		      {
			kinematics::D0Tag=1;
			Hep3Vector d0candidateMom=(allPB[i]+allPB[j]);
			float d0candidateE=EP1+EK;

			for(int k=0;k<allPB.size();k++)
			  {
			    if(k==j)
			      continue;
			    if(allPB_Class[k]==2 && allPB_Charge[k]>0)
			      {
				if(allPB[k].perp()<0.1)
				  continue;
				float EP2=sqrt(m_mass*m_mass+allPB[k].mag2());
				float dStarcandidateMass=sqrt((d0candidateE+EP2)*(d0candidateE+EP2)-(d0candidateMom+allPB[k]).mag2());
				//martins cut... 
				Hep3Vector dStarCandMom=d0candidateMom+allPB[k];
				if(dStarCandMom.mag()>2.0 && dStarCandMom.mag()<4.9)
				  {
				    histoDStar->Fill(dStarcandidateMass);
				    if(fabs(dStarcandidateMass-DStarMass)<0.05)
				      {
					kinematics::DStarTag=1;
					histoPiSlowMom->Fill(allPB[k].mag());
				      }

				  }
				
			      }

			  }

		      }
		  }


	      }

	  }
      }


  }


  bool bToDDoubleStar::getDDecayProducts(const Gen_hepevt gen_it, int& Kp, int& Km, int& Ks, int& Pip, int& Pim, int& Pi0, int& other)
  {
    int lund=gen_it.idhep();
    if(lund==211)
      {
	Pip++;
	return true;}
    if(lund ==-211)
      {
	Pim++;return true;}
    if(lund==321)
      {
	Kp++;
	return true;
      }
    if(lund==-321)
      {
	Km++;
	return true;
      }
    if(lund==111)
      {
	Pi0++;
	return true;
      }
    if(lund==310)
      {
	Ks++;
	return true;
      }
    genhep_vec* daughters=getDaughters(gen_it);
    //none of the above and not a photon
    if(fabs(lund)!=211 && fabs(lund)!=321 && lund !=111  && lund !=310 && lund !=22 && daughters->size()<=0)    
      {
	other++;
	return true;
      }

    if(daughters->size()<=0)
      {
	return true;
      }


    for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
      {
	//	cout <<"looking at lund; "<< (*it)->idhep()<<endl;
	getDDecayProducts(**it,Kp,Km,Ks,Pip,Pim,Pi0,other);
      }

  }

  void bToDDoubleStar::computeD_BR_CorrectionFactor(double& corrFact, double& corrUncert, int Kp, int Km, int Ks, int Pip,int Pim, int Pi0,int other, int& decType)
  {
    float mc=1.0;
    float data=1.0;
    float error;
    bool found=false;
    //    cout << "D decay: Kp: " << Kp << " Km: "<< Km <<" Ks: "<< Ks <<" Pip: " << Pip << " pim: " << Pim << " pi0: " << Pi0 << " other: "<< other <<endl;
    //Km Pip

    //the error here is on the data, so the uncertainty on the weight is err/mc

    if(Km==1&& Pip==1 && Pim==0 && Ks==0 && Kp==0 && Pi0==0&& other==0)
      {
	decType=0;
	mc=0.0382;
	data=0.0388;
	error=0.0005;
	found=true;

      }
    //km pip pi0
    if(Km==1&& Pip==1 && Pim==0 && Ks==0 && Kp==0 && Pi0==1 && other==0)
      {
	decType=1;
	mc=0.130811;
	data=0.139;
	error=0.005;
	found=true;
      }
    //Km Pip Pip Pim
    if(Km==1&& Pip==1 && Pim==1 && Ks==0 && Kp==0 && Pi0==0 && other==0)
      {
	decType=2;
	mc=0.0708693;
	data=0.0808;
	error=0.002;
	found=true;
      }
    //ks pip pim
    if(Km==0&& Pip==1 && Pim==1 && Ks==1 && Kp==0 && Pi0==0 && other==0)
      {
	decType=3;
	mc=0.0283637;
	data=0.0283;
	error=0.002;
	found=true;
      }
    //ks pip pim pi0
    if(Km==0 && Pip==1 && Pim==1 && Ks==1 && Kp==0 && Pi0==1 && other==0)
      {
	decType=4;
	mc=0.0517367;
	data=0.052;
	error=0.006;
	found=true;
      }
    //KsPi0
    if(Km==0 && Pip==0 && Pim==0 && Ks==1 && Kp==0 && Pi0==1 && other==0)
      {
	decType=5;
	mc=0.0113;
	data=0.0119;
	error=0.0004;
	found=true;

      }

    //kp km
    if(Km==1 && Pip==0 && Pim==0 && Ks==0 && Kp==1 && Pi0==0 && other==0)
      {
	decType=6;
	mc=0.0039;
	data=0.00396;
	error=0.00008;
	found=true;
      }
    //pip pim 
    if(Km==0 && Pip==1 && Pim==1 && Ks==0 && Kp==0 && Pi0==0 && other==0)
      {
	decType=7;
	mc=0.0014;
	data=0.001402;
	error=0.000026;
	found=true;
      }
    //ks ks
    if(Km==0 && Pip==0 && Pim==0 && Ks==2 && Kp==0 && Pi0==0 && other==0)
      {
	decType=8;
	mc=0.0004;
	data=0.00017;
	error=0.00004;
	found=true;
      }
    //pi0 pi0
    if(Km==0 && Pip==0 && Pim==0 && Ks==0 && Kp==0 && Pi0==2 && other==0)
      {
	decType=9;
	mc=0.0008;
	data=0.00082;
	error=0.000035;
	found=true;
      }
    //km pip pip
    if(Km==1 && Pip==2 && Pim==0 && Ks==0 && Kp==0 && Pi0==0 && other==0)
      {
	decType=10;
	mc=0.0950633;
	data=0.0913;
	error=0.0019;
	found=true;
      }
    //km pip pip pi0
    if(Km==1 && Pip==2 && Pim==0 && Ks==0 && Kp==0 && Pi0==1 && other==0)
      {
	decType=11;
	mc=0.0601865;
	data=0.0599;
	error=0.0018;
	found=true;
      }

    //ks pip
    if(Km==0 && Pip==1 && Pim==0 && Ks==1 && Kp==0 && Pi0==0 && other==0)
      {
	decType=12;
	mc=0.0147;
	data=0.0147;
	error=0.0007;
	found=true;
      }

    //ks pip pi0
    if(Km==0 && Pip==1 && Pim==0 && Ks==0 && Kp==0 && Pi0==1 && other==0)
      {
	decType=13;
	mc=0.0650925;
	data=0.0699;
	error=0.0027;
	found=true;
      }
    //kp km pip
    if(Km==1 && Pip==1 && Pim==0 && Ks==0 && Kp==1 && Pi0==0 && other==0)
      {
	decType=14;
	mc=0.00905861;
	data=0.00954;
	error=0.00026;
	found=true;
      }
    //ks kp
    if(Km==0 && Pip==0 && Pim==0 && Ks==1 && Kp==1 && Pi0==0 && other==0)
      {
	decType=15;
	mc=0.00295;
	data=0.00283;
	error=0.00016;
	found=true;
      }
    //ks pip pip pim
    if(Km==0 && Pip==2 && Pim==1 && Ks==0 && Kp==0 && Pi0==0 && other==0)
      {
	decType=16;
	mc=0.0315685;
	data=0.0312;
	error=0.0011;
	found=true;
      }
    //km pip pip pim pi0
    if(Km==1 && Pip==2 && Pim==1 && Ks==0 && Kp==0 && Pi0==1 && other==0)
      {
	decType=17;
	mc=0.0398269;
	data=0.042;
	error=0.004;
	found=true;
      }
    //pip pim pi0

    if(Km==0 && Pip==1 && Pim==1 && Ks==0 && Kp==0 && Pi0==1 && other==0)
      {
	decType=18;
	mc=0.0139744;
	data=0.0143;
	error=0.0006;
	found=true;
      }

    //pip pim pi0 pi0
    if(Km==0 && Pip==1 && Pim==1 && Ks==0 && Kp==0 && Pi0==2 && other==0)
      {
	decType=19;
	mc=0.00528788;
	data=0.01;
	error=0.0009;
	found=true;
      }

    //pip pip pim pi0
    if(Km==0 && Pip==2 && Pim==1 && Ks==0 && Kp==0 && Pi0==1 && other==0)
      {
	decType=20;
	mc=0.0115474;
	data=0.0113;
	error=0.0008;
	found=true;
      }

    //ks pi0 pi0
    if(Km==0 && Pip==0 && Pim==0 && Ks==1 && Kp==0 && Pi0==2 && other==0)
      {
	decType=21;
	mc=0.00927935;
	data=0.0091;
	error=0.0011;
	found=true;
      }

    //pip pi0
    if(Km==0 && Pip==1 && Pim==0 && Ks==0 && Kp==0 && Pi0==1 && other==0)
      {
	decType=22;
	mc=0.0026;
	data=0.00119;
	error=0.00006;
	found=true;
      }
    //pip pip pim
    if(Km==0 && Pip==2 && Pim==1 && Ks==0 && Kp==0 && Pi0==0 && other==0)
      {
	decType=23;
	mc=0.0037091;
	data=0.00318;
	error=0.00018;
	found=true;
      }
    //km pip pip pip pim
    if(Km==1 && Pip==3 && Pim==1 && Ks==0 && Kp==0 && Pi0==0 && other==0)
      {
	decType=24;
	mc=0.00620297;
	data=0.0056;
	error=0.0005;
	found=true;
      }

    //factor to correct mc with
    if(found)
      {
	if(mc==0 || data==0)
	  {
	    cout <<"corr fact data or mc equals zero " << data << " " << mc <<endl;
	    exit(1);
	  }
	corrFact=data/mc;
	corrUncert=error/mc;
      }
    else
      {
	corrFact=1.0;
	corrUncert=0.0;
      }

  }


  //looking at the BR corrections, they only seem to apply for D (so not D* etc)
  bool bToDDoubleStar::isAnyD(int lund)
  {
    if(fabs(lund)==411)
      return true;
    if(fabs(lund)==421)
      return true;

    return false;
  }


  bool bToDDoubleStar::recursivePrint(const Gen_hepevt gen_it, string s)
  {
    //    cout <<"looking at idhep: "<< gen_it.idhep()<<endl;
    int lund=fabs(gen_it.idhep());
    if(lund==911|| lund>9000000 || lund==30343 || lund==100113)
      return false;
    genhep_vec* daughters=getDaughters(gen_it);


    Particle p(gen_it);

    if(lund== 100423 || lund ==100421 ||lund==PY_DStar_2S || lund==100411 || lund==100413 ){
      cout <<s << lund << endl;
    }
    else{
      //for some reason some ids make the particle class crash when asked for the name...
      if(lund<10000)
	cout <<s<<p.pType().name()<<" (" << gen_it.idhep()<<"), "<<" |p|: " <<p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;
      else
	cout  <<s<<" (" << gen_it.idhep()<<"), "<<" |p|: " <<p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;
    }
    if(daughters->size()<=0)
      {
	//	cout <<"has " << p.nChildren()<<" children " <<endl;
	return true;
      }
    //don't print eventual decay products of kaons etc...
    if(lund==211 || lund==321 || lund==13 || lund==111)
      return true;

    for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
      {
	//	cout <<"looking at lund; "<< (*it)->idhep()<<endl;
	recursivePrint(**it,s+("-->"));
      }
    cout <<s<<endl;
  }

  //save all ids from the decay of particle p
  bool bToDDoubleStar::getDecayIds(Particle& p,set<int>& chrId,set<int>& pi0Id, set<int>& gammaId, bool print)
  {

    //pi0s have children, but we still want to have them separately
    if(p.mdstPi0())
      {
	pi0Id.insert(p.mdstPi0().get_ID());
	if(print)
	  {
	    cout <<"getDecayIds::inserting pi0, ID: "<< p.mdstCharged().get_ID() <<" momentum: "<< p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;
	  }
      }


    if(p.nChildren()<=0)
      {
	if(p.mdstCharged())
	  {
	    chrId.insert(p.mdstCharged().get_ID());
	    if(print)
	      {
		cout <<"getDecayIds::inserting charged " <<p.pType().name() << " ,track, ID: "<< p.mdstCharged().get_ID() <<" momentum: "<< p.p().rho()<<endl;
	      }
	  }
	if(p.mdstPi0())
	  {
	    pi0Id.insert(p.mdstPi0().get_ID());
	    if(print)
	      {
		cout <<"getDecayIds::inserting pi0, ID: "<< p.mdstCharged().get_ID() <<" momentum: "<< p.p().rho()<<endl;
	      }
	  }
	if(p.mdstGamma())
	  gammaId.insert(p.mdstGamma().get_ID());
	
	return false;
      }
    for(int i=0;i<p.nChildren();i++)
      {
	getDecayIds(p.child(i),chrId,pi0Id,gammaId,print);
      }
  }
  bool bToDDoubleStar::getMassFromDecayParticles(Particle& p, HepLorentzVector& mom)
  {
    if(p.nChildren()<=0)
      {
	mom+=p.p();
	return true;
      }
    else
      {
	for(int i=0;i<p.nChildren();i++)
	  {
	    getMassFromDecayParticles(p.child(i),mom);
	  }
      }
    return true;
  }

  void bToDDoubleStar::printMCList()
  {

    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
    //      cout <<"iterating" <<endl;
    int motherGenId=0;

    Gen_hepevt ghp;

    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {
	int geantID=abs(gen_it->idhep());//plus is ok, since it is the abs value
	if(gen_it->mother())
	  {
	    motherGenId=gen_it->mother().idhep();
	  }

	Particle* np=new Particle(*gen_it);
	if(np&& np->pType())
	  {
	    //		cout <<"id: "<< np->pType().lund()<<endl;
	    //		cout <<" particle: " << np->pType().name()<<endl;

	    //		cout <<"momentum: " << np->ptot()<<endl;
	    if(np->mc())
	      {
		//		  cout <<"got mc" <<endl;
	      }
	    if(np->mdstTrk())
	      {
		//		   cout <<"got trk" << endl;
	      }
	    if(np->mdstCharged())
	      {
		//		   cout <<"got charged reference " <<endl;
	      }
	  }
	delete np;
      }
    //	cout <<"about to return ..." <<endl;

  }


  //check if this can be reconstructed in the channels that we are looking for
  //at least the D meson.

  //@param gen_it: The b meson id
  //@param ids  the ids of the decay particles
  //@param pids the pids of the decay particles
  //@numNu number of neutrinos in the decay
  void bToDDoubleStar::recCheck(const Gen_hepevt& gen_it, vector<int>& ids, vector<int>& pids, int& numNu, bool onlyCorrectDDecays)
  {
    //special case D's: check if they decay hadronically in a decay thta we can actually check.
    if(gen_it.idhep()==PY_D0)
      {
	//found D0 in the decay we want
	if(!isD0(gen_it) && onlyCorrectDDecays)
	  {
	    //	cout <<"isD0 but not right decay" <<endl;
	    noDRec=true;
	    return;
	  }
	else{
	  //      cout <<"is D0 and right decay" <<endl;
	}
      }


    if((gen_it.idhep()==PY_D || gen_it.idhep()==-PY_D))
      {
	if( !isChargedD(gen_it) && onlyCorrectDDecays)
	  {
	    //	    cout <<"is charged D but not right decay" <<endl;
	    noDRec=true;
	    return;
	  }
	else
	  {
	    //	    cout <<"is charged and right decay" <<endl;
	  }
      }
    if(gen_it.idhep()==PY_DStar || gen_it.idhep()==-PY_DStar || gen_it.idhep()==PY_DStar0)
      {
	if(!isDStar(gen_it) && onlyCorrectDDecays)
	  {
	    //	    cout <<"isD* but not right decay" <<endl;
	    noDRec=true;
	    return;
	  }
	else
	  {
	    //	    cout <<"is D* and correct decay " <<endl;
	  }
      }

    //leptons don't have daughters and checking for them leads to a crash for some reason (mother with ID 0)
    if(abs(gen_it.idhep()) <37)
      {
	if(isNu(gen_it.idhep()))
	  numNu++;
	if(detectable(gen_it.idhep()))
	  {
	    ids.push_back(gen_it.get_ID());
	    pids.push_back(gen_it.idhep());
	  }
	return;
      }
    genhep_vec* daughters=getDaughters(gen_it);
    //    cout <<"looking at " << daughters->size() <<" daughters " <<endl;
    if(daughters->size()==0)
      {
	delete daughters;
	return;
      }
    else
      {
	for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
	  {
	    int daughterId=(*it)->idhep();
	    int daughterHepevtId=(*it)->get_ID();
	    //	    cout << " daughterPID: "<< daughterId <<endl;
	    if(detectable(daughterId))
	      {

		ids.push_back(daughterHepevtId);
		pids.push_back(daughterId);
	      }
	    else
	      {
		//		cout <<"calling rec check with " << daughterId <<endl;
		recCheck(**it,ids,pids,numNu, onlyCorrectDDecays);
	      }
	  }

      }

    delete daughters;

    return;

  }

  //should maybe use the recursive method 
  //this method checks if the Dpipi decay is in the MC and also gets the BR correction factors
  bool bToDDoubleStar::checkForDPiPi(int& bMesonId,bool& foundSinglePionDecay, bool& foundDlNu,bool print)
  {
    foundSinglePionDecay=false;
    overlapFractionCharged=.0;
    overlapFractionPi0=.0;
    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();

    //////////////////////////
    //find best b
    float bestDistance=-1;
    Gen_hepevt_Manager::iterator bestB;
    int numBs=0;
    int oldLund=0;
    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {
	int geantID=abs(gen_it->idhep());//plus is ok, since it is the abs value
	if(geantID==PY_B0 || geantID==PY_B)
	  {
	    if(numBs==1 && oldLund==gen_it->idhep())
	      {
		//		cout <<"found " <<oldLund << " again !! " <<endl;
	      }
	    oldLund=gen_it->idhep();
	    numBs++;
	    //	    cout <<endl<<endl<<"top level B:"<<endl;
	    //	    recursivePrint(*gen_it,"");
	    float distanceToBestB= (gen_it->PX()-bestBPx)*(gen_it->PX()-bestBPx)+(gen_it->PY()-bestBPy)*(gen_it->PY()-bestBPy)+(gen_it->PZ()-bestBPz)*(gen_it->PZ()-bestBPz);
	    //	    cout <<"distance to best b is: "<< distanceToBestB << " (before " << bestDistance <<" ) " <<endl;
	    if(distanceToBestB<bestDistance|| bestDistance<0)
	      {
		//		cout <<"best so far.. " <<endl;
		bestDistance=distanceToBestB;
		bestB=gen_it;
	      }

	  }
      }
    //    cout <<" found " << numBs << " numBs " <<endl;
  
    ////////////////////
    D_BR_CorrectionFactor=1.0;
    sRelD_BR_CorrFactor=0.0;
    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {
	int geantID=abs(gen_it->idhep());//plus is ok, since it is the abs value
	//for D decay branching ratio corrections
	if(isAnyD(geantID))
	  {
	    //	    cout <<"found D in event " <<endl;
	    int Kp=0;
	    int Km=0;
	    int Ks=0;
	    int Pip=0;
	    int Pim=0;
	    int Pi0=0;
	    int other=0;
	    //	    cout <<"found d, lets look at decay " <<endl;
	    //	    recursivePrint(*gen_it,"");
	    getDDecayProducts(*gen_it, Kp,Km,Ks,Pip,Pim,Pi0,other);
	    //	    cout <<" we found " << Kp << " K+ " << ", " << Km <<" K- " << Ks << " Ks, " << Pip <<" Pi+, " << Pim << " Pi- , " << Pi0 << " pi0, " << other<<" others " <<endl;
	    double tmpCorrFact=1.0;
	    double tmpRelErr=0.0;
	    int dDecType=-1;
	    computeD_BR_CorrectionFactor(tmpCorrFact,tmpRelErr,Kp,Km,Ks,Pip,Pim,Pi0,other,dDecType);
	    //do this before we charge conjugate
	    if(dDecType>=0)
	      {
		treeData.D_DecayCorr[treeData.numDDecCorr]=tmpCorrFact;
		treeData.sD_DecayCorr[treeData.numDDecCorr]=tmpRelErr*tmpCorrFact;
		treeData.D_decType[treeData.numDDecCorr]=dDecType;
		treeData.numDDecCorr++;
	      }
	    //	    cout <<"computed corr factor: "<< tmpCorrFact<<endl;
	    //and check for conjugate decay
	    int tmp=Km;
	    Km=Kp;
	    Kp=tmp;
	    tmp=Pim;
	    Pip=Pim;
	    Pim=tmp;
	    dDecType=-1;
	    //	    cout <<"and conjugate.. " <<endl;
	    computeD_BR_CorrectionFactor(tmpCorrFact,tmpRelErr,Kp,Km,Ks,Pip,Pim,Pi0,other,dDecType);
	    if(dDecType>=0)
	      {
		treeData.D_DecayCorr[treeData.numDDecCorr]=tmpCorrFact;
		treeData.sD_DecayCorr[treeData.numDDecCorr]=tmpRelErr*tmpCorrFact;
		treeData.D_decType[treeData.numDDecCorr]=dDecType;
		treeData.numDDecCorr++;
	      }

	    //	    cout <<" and after charge conj:  we found " << Kp << " K+ " << ", " << Km <<" K- " << Ks << " Ks, " << Pip <<" Pi+, " << Pim << " Pi- , " << Pi0 << " pi0, " << other<<" others " <<endl;
	    //	    cout <<" and corr fact now: " << tmpCorrFact<<endl;
	    //might be several D's...
	    //	    cout <<"overall D correction factor " <<D_BR_CorrectionFactor <<endl;

	  }
      }
    //////////////////////
    /////////////done D_BR correction
    /////////////////////////

    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {
	int geantID=abs(gen_it->idhep());//plus is ok, since it is the abs value
	//found B 
	bool foundNu=false;
	bool foundL=false;
	bool foundD=false;
	bool foundDStar=false;
	bool foundPiPlus=false;
	bool foundPiMinus=false;
	bool foundSameChargePions=false;

	if(geantID==PY_B0 || geantID==PY_B)
	  {
	    if(geantID==PY_B0)
	      treeData.mcBCharge=0;
	    else
	      treeData.mcBCharge=1;

	    bool tempFoundDDoubleStar=false;
	    int tempNumLeptons=0;
	    int tempNumPions=0;
	    int tempNumKaons=0;
	    int tempNumRhos=0;
	    int tempNumPi0=0;
	    int tempNumBaryons=0;
	    int tempNumD=0;
	    int tempDCharge=0;
	    int tempNumDStar=0;
	    bool tempNumDStar2S=false;
	    bool tempNumDStarD2S=false;
	    int tempNumNu=0;
	    //general check if we find one of the D(*)npi l nu decays
	    int br_sigs[11];
	    for(int i=0;i<=br_sig_D0Star0;i++)
	      {
		br_sigs[i]=0;
	      }
	    HepLorentzVector p_D,p_l,p_nu;
	    findDecaySignatureForBBRCorrection(*gen_it,tempNumLeptons,tempNumPions,tempNumKaons,tempNumPi0,tempNumBaryons, tempNumNu,br_sigs,p_D,p_l,p_nu);

	    ////////--
	    //	    cout <<" get decay sig:numpions: "<< tempNumPions <<endl;

	    //	        cout <<"found " << tempNumD << " Ds " << tempNumNu <<" Neutrinos " << tempNumPions <<" pions " << tempNumLeptons << " leptons " << tempFoundDDoubleStar <<" doublestar " << tempNumDStar2S << " star 2S " << tempNumKaons <<" kaons " << tempNumPi0 <<" pi0s " << tempNumBaryons <<" baryons " << tempNumDStar << " dstar " <<endl;

	    int foundBR=0;
	    //	    cout <<" heck.." << endl;
	    for(int i=0;i<=br_sig_D0Star0;i++){foundBR+=br_sigs[i];
	      if(br_sigs[i])
		{
		  //		cout <<"found br : "<< i <<", " << foundBR<<endl;
		}
};
	    if(foundBR>0)
	      {
		//		recursivePrint(*gen_it, "");
	      }
	    //	    if(temp>1)
	    //	      cout <<" found more than one D(*(*)) in the decay! " << endl;
	    if(tempNumNu==1 && tempNumLeptons==1 && tempNumPions==0&& tempNumKaons==0 && tempNumPi0==0 && tempNumBaryons==0 && foundBR)
	      {
		    HepLorentzVector signalBMomentum=Particle(*gen_it).p();
		    //copy of the vector to boost into B system
		    HepLorentzVector p_l2=p_l;
		    p_l2.boost(-(signalBMomentum.boostVector()));
		    ff_pL=p_l2.vect().mag();
		    w=(p_D.dot(signalBMomentum))/sqrt(p_D.dot(p_D)*signalBMomentum.dot(signalBMomentum));
		    HepLorentzVector p_W = p_l + p_nu;
		    p_l.boost(-(p_W.boostVector()));
		    p_D.boost(-(p_W.boostVector()));
		    cosTheta = cos( p_l.vect().angle( p_D.vect() ) );
		    q2 =(p_l+p_nu)*(p_l+p_nu);
		    //		    cout <<"p_l: "<< p_l <<" p_nu "<< p_nu <<endl;
		    //		    cout <<"we have pL: "<< ff_pL <<" w: "<< w <<" cosTheta" << cosTheta <<" q2: "<< q2 <<endl;
		    //		    cout <<"dssIdx: " << dssIdx << " ltype: "<< lType <<endl;
		  }
		else{
		  //if we didn't find a decay of interest reset our flags so we can look again...
		  dssIdx=-1;
		  lType=-1;
		}
	  

	    //implement FF corrections here, because we already determined the decay
	    //however, we need the D, l and nu 4-momenta to compute the kinematics needed for the FF
	    if(tempNumNu==1 && tempNumLeptons==1 && tempNumPions==0&& tempNumKaons==0 && tempNumPi0==0 && tempNumBaryons==0)
	      {
		if( br_sigs[br_sig_D0]==1)
		  {
		    br_sig_D0LNu=true;
		  }
		if(br_sigs[br_sig_D]==1)
		  {
		    br_sig_DLNu=true;
		  }

		if(br_sigs[br_sig_DStar]==1)
		  {
		    br_sig_DStarLNu=true;
		  }
		if( br_sigs[br_sig_DStar0]==1)
		  {
		    br_sig_DStar0LNu=true;
		  }
		if( br_sigs[br_sig_D1]==1)
		  {
		    br_sig_D1LNu=true;
		  }
		if( br_sigs[br_sig_D2]==1)
		  {
		    br_sig_D2LNu=true;
		  }
		if( br_sigs[br_sig_D1Prime]==1)
		  {
		    br_sig_D1PrimeLNu=true;
		  }
		if( br_sigs[br_sig_D0Star]==1)
		  {
		    br_sig_D0StarLNu=true;
		  }
		if( br_sigs[br_sig_D10]==1)
		  {
		    br_sig_D10LNu=true;
		  }
		if( br_sigs[br_sig_D1Prime0])
		  {
		    br_sig_D1Prime0LNu=true;
		  }
		if( br_sigs[br_sig_D0Star0]==1)
		  {
		    br_sig_D0Star0LNu=true;
		  }
	      }
	  

	    ///probably have to set the other stuff to zero again..
	    //
	    tempNumLeptons=0;
	    tempNumPions=0;
	    tempNumKaons=0;
	    tempNumRhos=0;
	    tempNumPi0=0;
	    tempNumBaryons=0;
	    tempNumD=0;
	    tempDCharge=0;
	    tempNumDStar=0;
	    tempNumDStar2S=false;
	    tempNumDStarD2S=false;
	    tempNumNu=0;
	    mc_piMom.clear();
	    mc_piTheta.clear();
	    mc_piPhi.clear();
	    mc_piFound.clear();

	    //and the other call for the 'regular' decay signature search where we trace the D decays
	    findDecaySignature(*gen_it,tempFoundDDoubleStar,tempNumLeptons,tempNumPions,tempNumKaons,tempNumRhos,tempNumPi0,tempNumBaryons,tempNumD, tempNumDStar, tempNumNu,tempNumDStar2S, tempNumDStarD2S, tempDCharge);
	if(dDoubleStarId==100423)
	  //	    if(evtNr==74850 &&  runNr== 879)
	      {
		//		recursivePrint(*gen_it,"");
		//		cout <<"tmpNumD: " << tempNumD << " numNu: " << tempNumNu <<" numLep: " << tempNumLeptons<<" Kaons? "<< tempNumKaons;
		//		cout  <<" pions? " << tempNumPions <<" pi0? " << tempNumPi0 << " baryons? " << tempNumBaryons;
		//		cout  <<" D*? " << tempNumDStar  <<" D**? " << tempFoundDDoubleStar <<" D*(2S)? " << tempNumDStar2S <<" D*(2SD) ? " << tempNumDStarD2S << endl;
	      }
	    if((tempNumDStar==1|| tempNumD==1) && tempNumNu==1 && tempNumLeptons==1)
	      {
		if(tempNumKaons==0 && tempNumPi0==0 && tempNumBaryons==0)// && tempNumD==0)
		  {

		    ////////
		    ///////

		    if(tempNumPions==1 || tempNumPions==2)
		      {
			if(tempNumDStar==1)
			  treeData.mcIsDStar=1;
			else
			  treeData.mcIsDStar=0;

			treeData.mcDCharge=tempDCharge;
		      }

		    if(tempNumPions==1)
		      {
			//			if(dDoubleStarId==100423)
			//			  cout <<"found dpppi" <<endl;
			treeData.pi1Mom_mc=mc_piMom[0];
			treeData.pi1Theta_mc=mc_piTheta[0];
			treeData.pi1Phi_mc=mc_piPhi[0];
			treeData.pi1Found=mc_piFound[0];

		      }
		    if(tempNumPions==2)
		      {

			treeData.pi1Mom_mc=mc_piMom[0];
			treeData.pi1Theta_mc=mc_piTheta[0];
			treeData.pi1Phi_mc=mc_piPhi[0];
			treeData.pi1Found=mc_piFound[0];
			treeData.pi2Mom_mc=mc_piMom[1];
			treeData.pi2Theta_mc=mc_piTheta[1];
			treeData.pi2Phi_mc=mc_piPhi[1];
			treeData.pi2Found=mc_piFound[1];
		      }


		  }
	      }

	    if(tempNumD==1 && tempNumNu==1 && tempNumLeptons==1)
	      {
		if(tempNumKaons==0 && tempNumPi0==0 && tempNumBaryons==0 && tempNumDStar==0)
		  {
		    if(tempNumPions<3)
		      {
			//			if(dDoubleStarId==100423)
			//			  cout <<" and found rec " <<endl;
			vector<int> candidateIds;
			vector<int> candidatePids;
			int numNu;
			//to compute overlap with the tag B
			Particle p(*gen_it);
			//get decay id's doesn't work, because the gen_hep B particle class doesn't have the daufthers
			//		    getDecayIds(p,candidateCharged,candidatePi0,candidateGamma,true);
			recCheck(*gen_it, candidateIds,candidatePids, numNu,false);
			set<int> sCandidateIds(candidateIds.begin(),candidateIds.end());
			
			int chargeOverlap=getOverlap(sCandidateIds,chargedIds);
			int pi0Overlap=getOverlap(sCandidateIds,pi0Ids);
			if(chargedIds.size()>0)
			  overlapFractionCharged=(float)(chargeOverlap)/(float)(chargedIds.size());
	
			else
			  overlapFractionCharged=0;
			if(pi0Ids.size()>0)
			  overlapFractionPi0=(float)(pi0Overlap)/(float)(pi0Ids.size());
			else
			  overlapFractionPi0=0;

			//			cout <<"overlapFraction charged: "<< overlapFractionCharged <<" pi0: "<< overlapFractionPi0 <<endl;

			//			cout <<"chargeOverlap: "<< chargeOverlap <<" of " << chargedIds.size()<<endl;
			//			cout <<"pi0Overlap: "<< pi0Overlap <<" of " << pi0Ids.size()<<endl;
		      }

		    //		    if(gen_it==bestB && tempNumPions<3)
		    //		      cout <<"found Dlnu npi decay for tag b " <<endl;
		    //		    else
		    //		      cout <<"not best Dlnu..." <<endl;
		    if(!tempFoundDDoubleStar && !tempNumDStar2S && !tempNumDStarD2S)
		      {

			if(tempNumPions==0)
			  sigDLNu++;
			if(tempNumPions==1)
			  sigDPiLNu++;
			if(tempNumPions==2)
			  sigDPiPiLNu++;
		      }
		    else
		      {

			if(tempNumPions==0)
			  sigResDLNu++;
			if(tempNumPions==1)
			  sigResDPiLNu++;
			if(tempNumPions==2)
			  sigResDPiPiLNu++;
			//			if(dDoubleStarId==100423)
			//			  cout <<"foundres  dd double star, numpions:  " <<tempNumPions<< ", sigres: "<< sigResDPiLNu<<endl;


		      }
		  }

	      }
	  
	    ////add the other relevant hadronic channels

	    ///
	  
	    if(tempNumDStar==1 && tempNumNu==1 && tempNumLeptons==1)
	      {
		if(tempNumKaons==0 && tempNumPi0==0 && tempNumBaryons==0 && tempNumD==0)
		  {
		    if( !tempFoundDDoubleStar && !tempNumDStar2S && !tempNumDStarD2S)
		      {
			if(tempNumPions==0)
			  sigDStarLNu++;
			if(tempNumPions==1)
			  sigDStarPiLNu++;
			if(tempNumPions==2)
			  sigDStarPiPiLNu++;
		      }
		    else
		      {
			//			if(dDoubleStarId==100423)
			//			  cout <<"found dd double star D*2 " << endl;
			if(tempNumPions==0)
			  sigResDStarLNu++;
			if(tempNumPions==1)
			  sigResDStarPiLNu++;
			if(tempNumPions==2)
			  sigResDStarPiPiLNu++;
		      }
		  }
		
	      }
	  

	    //haven't found it yet... (note that this uses different variables. So not the tmp version from above which is only there to check if a certain decay mode is found
	    //the below is only for the case that there is a double star, to check for sig_numpi etc...
	    if(!sig_FoundDDoubleStar)
	      findDecaySignature(*gen_it,sig_FoundDDoubleStar,sig_numLeptons,sig_numPions,sig_numKaons,sig_numRhos,sig_numPi0,sig_numBaryons,sig_numD, sig_numDStar,tempNumNu, sig_dStar_2S, sig_d_2S, sig_DCharge);
	    if(!sig_FoundDDoubleStar)
	      {
		//still none, so reset fields
		sig_numLeptons=0;
		sig_numKaons=0;
		sig_numPions=0;
		sig_numPi0=0;
		sig_numBaryons=0;
		sig_numD=0;
		sig_DCharge=0;
		sig_numDStar=0;
		sig_numRhos=0;
	      }
	    else
	      {
		//		cout<<" found D** for decay of " << gen_it->idhep();
		//		cout <<" numLept: "<< sig_numLeptons<< " kaons: "<< sig_numKaons <<" pions: "<< sig_numPions <<" pi0s: "<< sig_numPi0 <<" num baryos: "<< sig_numBaryons;
		//		cout <<" numD: "<< sig_numD <<" DStar " << sig_numDStar <<" num rhos: "<< sig_numRhos <<endl;
	      }
	  
	    bMesonId= gen_it->get_ID();
	    genhep_vec* daughters=getDaughters(*gen_it);
	    //don't count gammas..
	    //	    int numDaughters=daughters->size();
	    int numDaughters=0;
	    int numGDaughters=0;
	    for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
	      {
		int daughterId=abs((*it)->idhep());
		//ignore photons...
		if(daughterId!=22)
		  numDaughters++;

		//higher order D that will decay... so we have to go down the chain...
		if(daughterId==100421 || daughterId==100423 ||daughterId==100413 || daughterId==100411)
		  {
		    //		    cout <<"found first level ho D" <<endl;

		    if(daughterId==100423 || daughterId==10013)
		      found_2SD_Star=true;
		    if(daughterId==100421 || daughterId==10011)
		      found_2SD=true;
		    genhep_vec* newDaughters=getDaughters(**it);
		    for(genhep_vec::iterator it2=newDaughters->begin();it2!=newDaughters->end();it2++)
		      {
			int tmpDaughterId=abs((*it2)->idhep());
			if(tmpDaughterId==PY_DStar_00|| tmpDaughterId==PY_DStar0Plus|| tmpDaughterId==PY_DStar_0s|| tmpDaughterId==PY_DStar0Plus|| tmpDaughterId==PY_D_1 || tmpDaughterId==PY_D_10 || tmpDaughterId==PY_DStar_2 || tmpDaughterId==PY_DStar_20)
			  {
			    checkDDoubleStarDecay(foundSinglePionDecay,**it2,tmpDaughterId, numGDaughters,foundSameChargePions, foundPiPlus,foundPiMinus,foundD,foundDStar,print);
			  }
			if(tmpDaughterId==100421 || tmpDaughterId==100423 || tmpDaughterId==100411|| tmpDaughterId==100413)
			  {
			    //			    cout <<"found second level ho D" <<endl;
			    genhep_vec* newDaughters2=getDaughters(**it2);
			    for(genhep_vec::iterator it3=newDaughters2->begin();it3!=newDaughters2->end();it3++)
			      {
				int tmpDaughterId2=abs((*it3)->idhep());
				//the check for 2S should never be positive, since we already handle that case...
				if(tmpDaughterId2==PY_DStar_00|| tmpDaughterId2==PY_DStar0Plus|| tmpDaughterId2==PY_DStar_0s|| tmpDaughterId2==PY_DStar0Plus|| tmpDaughterId2==PY_D_1 || tmpDaughterId2==PY_D_10 || tmpDaughterId2==PY_DStar_2 || tmpDaughterId2==PY_DStar_20)
				  {
				    //				    cout <<"found nested D**" <<endl;
				    bool ret=checkDDoubleStarDecay(foundSinglePionDecay,**it3,tmpDaughterId, numGDaughters,foundSameChargePions, foundPiPlus,foundPiMinus,foundD,foundDStar,print);
				    //				    cout <<"found D? " << foundD <<" found D*? " << foundDStar<<endl;
				    

				    if(tmpDaughterId2==100421 || tmpDaughterId2==100411)
				      found_2SD=true;
				  }

			      }

			  }

		      }


		  }
		int decType=-1;
		pair<float,float> ret=getBRCorrection(decType);
		B_BR_CorrectionFactor=ret.first;
		treeData.B_DecayCorr=B_BR_CorrectionFactor;
		treeData.B_decType=decType;
		treeData.sB_DecayCorr=ret.second;
		//		cout <<" sB factor" << ret.second <<" corr: " << ret.first <<endl;
		//		cout <<"setting sB_decaycorr to " << ret.second <<endl;
		//cout <<"B cor factor: "<< B_BR_CorrectionFactor << " D: "<< D_BR_CorrectionFactor <<endl;
		float mom=sqrt((*it)->PX()*(*it)->PX()+(*it)->PY()*(*it)->PY()+(*it)->PZ()*(*it)->PZ());
		if(daughterId==12|| daughterId==14|| daughterId==16)
		  {
		    foundNu=true;
		    continue;
		  }
		if(daughterId==11 || daughterId==13 || daughterId==15)
		  {
		    foundL=true;
		    treeData.leptonP=mom;
		    continue;
		  }

		//D*_00,  D*_0+ D*_0s  (D*_0) D_1+
		if(daughterId==PY_DStar_00|| daughterId==PY_DStar0Plus|| daughterId==PY_DStar_0s|| daughterId==PY_DStar0Plus|| daughterId==PY_D_1 || daughterId==PY_D_10 || daughterId==PY_DStar_2 || daughterId==PY_DStar_20)
		  {
		    treeData.dPID=daughterId;
		    treeData.dMesonP=mom;
		    genhep_vec* grandDaughters=getDaughters(**it);
		    Particle mother(**it);
		    if(daughterId==100423 || daughterId==100421)
		      found_2SD=true;
		    if(print && (daughterId!=100423 && daughterId!=100421))
		      {
			cout << " D is " << mother.pType().name() <<", id: " << (*it)->get_ID()<<" lund: " << mother.lund()<<" mom: " << mother.p().rho()<<" ("<<mother.p().px()<<", " << mother.p().py()<< ", " << mother.p().pz() <<", "<< mother.p().t()<<")" <<endl;
		      }
		    numGDaughters=0;
		    for(genhep_vec::iterator it2=grandDaughters->begin();it2!=grandDaughters->end();it2++)
		      {
			float gdMom=sqrt((*it2)->PX()*(*it2)->PX()+(*it2)->PY()*(*it2)->PY()+(*it2)->PZ()*(*it2)->PZ());
			int gDaughterId=fabs((*it2)->idhep());
			if(gDaughterId!=22)
			  numGDaughters++;
			int signedGDaughterId=(*it2)->idhep();
			//			cout <<"gDaughterid: "<< gDaughterId <<endl;
			
			Particle p(**it2);
			if(print && gDaughterId!=100423 && gDaughterId!=100421)
			  cout <<"daughters: "<< p.pType().name()<<", id: " << (*it2)->get_ID()<<" lund: " << p.lund()<<" mom: " << p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;


			//			  cout <<"grand daughters: "<< p.pType().name()<<" lund: " << p.lund()<<endl;
			//D*, D*0, D+, D0
			if(gDaughterId==PY_DStar || gDaughterId==PY_DStar0 || gDaughterId==PY_D || gDaughterId==PY_D0)
			  {
			    treeData.daughterDPID=abs(p.lund());
			    genhep_vec* ggDaughters=getDaughters(**it2);
			    for(genhep_vec::iterator it3=ggDaughters->begin();it3!=ggDaughters->end();it3++)
			      {
				int ggDaughterId=(*it3)->idhep();
				int signedGGDaughterId=(*it3)->idhep();
				Particle pp(**it3);
				if(abs(pp.lund())==PY_D0 || abs(pp.lund())==PY_D)
				  {
				    genhep_vec* gggDaughters=getDaughters(**it3);
				    for(genhep_vec::iterator it4=gggDaughters->begin();it4!=gggDaughters->end();it4++)
				      {
					int gggDaughterId=(*it4)->idhep();
					int signedGGGDaughterId=(*it4)->idhep();
					Particle ppp(**it4);
					if(print)
					  cout <<"   ->->D daughters: "<< ppp.pType().name()<<", id: " << (*it4)->get_ID()<<" lund: " << ppp.lund()<<" mom: " << ppp.p().rho()<<" ("<<ppp.p().px()<<", " << ppp.p().py()<< ", " << ppp.p().pz() <<", "<< ppp.p().t()<<")" <<endl;
					//the A0 seems to lead to a crash, probably because daughter particles are fishy..
					if(abs(ppp.lund())==213 || abs(ppp.lund())==113|| abs(ppp.lund())==10323 )
					  {
					    //					    cout <<"getting more daughters.." <<endl;
					    genhep_vec* ggggDaughters=getDaughters(**it4);
					    //					    cout <<" gggg done.." <<ggggDaughters->size()<<endl;

					    for(genhep_vec::iterator it5=ggggDaughters->begin();it5!=ggggDaughters->end();it5++)
					      {
						//						cout <<"-->1" <<endl;
						int ggggDaughterId=(*it5)->idhep();
						int signedGGGGDaughterId=(*it5)->idhep();

						Particle pppp(**it5);
						if(print)
						  cout <<"   ->->->D daughters: "<< pppp.pType().name()<<", id: " << (*it5)->get_ID()<<" lund: " << pppp.lund()<<" mom: " << pppp.p().rho()<< " ("<<pppp.p().px()<<", " << pppp.p().py()<< ", " << pppp.p().pz() <<", "<< pppp.p().t()<<")" <<endl;
					      }

					    delete ggggDaughters;
					  }
				      }

				  
				    delete gggDaughters;
				  }
				if(print)
				  cout <<"   ->D daughters: "<< pp.pType().name()<<", id: " << (*it3)->get_ID()<<" lund: " << pp.lund()<<" mom: " << pp.p().rho()<<" ("<<pp.p().px()<<", " << pp.p().py()<< ", " << pp.p().pz() <<", "<< pp.p().t()<<")" <<endl;
			      }
			    delete ggDaughters;
			  }

			//already found piplus...
			if(signedGDaughterId==PY_PI && foundPiPlus)
			  {
			    cout <<"found second piPlus!" <<endl;
			    foundSameChargePions=true;
			  }
			if(signedGDaughterId==PY_PI)
			  {
			    treeData.piPlusP=gdMom;
			    //				cout <<"pion mom: "<< (*it2)->PX() <<", " << (*it2)->PY() <<", " << (*it2)->PZ()<<endl;
			    //				cout <<"id: "<< (*it2)->get_ID()<<endl;
			    //				if(p.mdstVee()!=0)
			    //				  {
			    //  cout <<"found mdst Vee " <<endl;
			    //			

			    foundPiPlus=true;
			    continue;
			  }
			if(signedGDaughterId==-PY_PI && foundPiMinus)
			  {
			    cout <<"found second piMinus!" <<endl;
			    foundSameChargePions=true;
			  }
			if(signedGDaughterId==(-PY_PI))
			  {
			    treeData.piMinusP=gdMom;
			    foundPiMinus=true;
			    continue;
			  }
			if(print)
			  {
			    if( gDaughterId!=100423 && gDaughterId!=100421)
			      cout <<"gDaughter " << p.pType().name()<< ", id: "<< gDaughterId<<", ID: "<<(*it2)->get_ID() <<" mom: " << p.p().rho()<<endl;
			    else
			      cout <<"gDaughter: "<< gDaughterId <<endl;
			  }
			//D
			if(gDaughterId==PY_D || gDaughterId==PY_D0)
			  {
			    foundD=true;
			    continue;
			  }
			//D*
			if(PY_DStar==gDaughterId|| gDaughterId==  PY_DStar0)
			  {
			    foundDStar=true;
			    continue;
			  }
		      }
		    delete grandDaughters;
		  }


	      }

	    //num daughters has to be 3 (D** l,nu), gdaughters has to be 2 (D pi) or 3 (D pi pi)
	    if((numDaughters!=3 || numGDaughters!=3) &&(foundD|| foundDStar) && foundL && foundNu&& foundPiPlus && foundPiMinus && !foundSameChargePions)
	      {
		//		cout <<" found all, but wrong count, daughters are:  ";
		int numDaughters=daughters->size();
		for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
		  {
		    int daughterId=abs((*it)->idhep());
		    Particle dp(**it);
		    //		cout << dp.pType().name()<<endl;
		
		  }
	      }
	    if(numDaughters==3 && numGDaughters==3 &&(foundD|| foundDStar) && foundL && foundNu&& foundPiPlus && foundPiMinus && !foundSameChargePions)
	      {
		if(print)
		  cout <<"=================================="<<endl<<"            found Dpipi!! " <<endl;
		//		cout <<"======================="<<endl;
		delete daughters;
		bMesonId= gen_it->get_ID();
		return true;
	      }
	    else
	      {
		//only one pion 
		if(print)
		  {
		    cout <<"only one pion, num Daughters: " << numDaughters <<" gDaughters: " << numGDaughters;
		    cout <<" foundD: " << foundD <<" DStar? " << foundDStar <<" L?: "<< foundL;
		    cout <<" nu? " << foundNu << " pi-? " << foundPiMinus <<" pi+? " << foundPiPlus<<endl;
		  }

		if(numDaughters==3 && numGDaughters==2&& (foundD|| foundDStar) && foundL && foundNu && (foundPiPlus || foundPiMinus) && !foundSameChargePions && !(foundPiPlus && foundPiMinus))
		  {
		    foundSinglePionDecay=true;
		    if(print)
		      cout <<"foundSinglePionDecay!" <<endl;
		    bMesonId= gen_it->get_ID();
		    //return, since we most likely won't find the two pion decay in the other B
		    return false;
		  }
		if(numDaughters==3 && numGDaughters==2&& (foundD|| foundDStar) && foundL && foundNu && !foundPiPlus && !foundPiMinus && !foundSameChargePions && !(foundPiPlus && foundPiMinus))
		  {
		    foundDlNu=true;
		    if(print)
		      cout <<"foundSinglePionDecay!" <<endl;
		    bMesonId= gen_it->get_ID();
		    //return, since we most likely won't find the two pion decay in the other B
		    return false;
		  }

	      }
	    delete daughters;
	  }
	//reset for the next B
	foundD=false;
	foundDStar=false;
	foundL=false;
	foundNu=false;
	foundPiPlus=false;
	foundPiMinus=false;
	foundSameChargePions=false;
	found_2SD=false;
	found_2SD_Star=false;
      }

    return false;
  }
  ///end of checkfordpipi



  bool bToDDoubleStar::checkDDoubleStarDecay(bool& foundSinglePionDecay,const Gen_hepevt &m_mother,int daughterId,int& numGDaughters, bool& foundSameChargePions, bool& foundPiPlus, bool& foundPiMinus,  bool& foundD, bool& foundDStar, bool print)
  {
    if(print)
      cout <<"check D Double STAR Decay" <<endl;
    float mom=sqrt(m_mother.PX()*m_mother.PX()+m_mother.PY()*m_mother.PY()+m_mother.PZ()*m_mother.PZ());
    treeData.dPID=daughterId;
    treeData.dMesonP=mom;
    genhep_vec* grandDaughters=getDaughters(m_mother);
    Particle mother(m_mother);
    if(print && (daughterId!=100423 && daughterId!=100421))
      cout << " D is " << mother.pType().name() <<", id: " << m_mother.get_ID()<<" lund: " << mother.lund()<<" mom: " << mother.p().rho()<<" ("<<mother.p().px()<<", " << mother.p().py()<< ", " << mother.p().pz() <<", "<< mother.p().t()<<")" <<endl;
    numGDaughters=0;
    for(genhep_vec::iterator it2=grandDaughters->begin();it2!=grandDaughters->end();it2++)
      {
	float gdMom=sqrt((*it2)->PX()*(*it2)->PX()+(*it2)->PY()*(*it2)->PY()+(*it2)->PZ()*(*it2)->PZ());
	int gDaughterId=fabs((*it2)->idhep());
	if(gDaughterId!=22)
	  numGDaughters++;
	int signedGDaughterId=(*it2)->idhep();
	//			cout <<"gDaughterid: "<< gDaughterId <<endl;
			
	Particle p(**it2);
	if(print && gDaughterId!=100423 && gDaughterId!=100421)
	  cout <<"daughters: "<< p.pType().name()<<", id: " << (*it2)->get_ID()<<" lund: " << p.lund()<<" mom: " << p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;


	//			  cout <<"grand daughters: "<< p.pType().name()<<" lund: " << p.lund()<<endl;
	//D*, D*0, D+, D0
	if(gDaughterId==PY_DStar || gDaughterId==PY_DStar0 || gDaughterId==PY_D || gDaughterId==PY_D0)
	  {
	    treeData.daughterDPID=abs(p.lund());
	    genhep_vec* ggDaughters=getDaughters(**it2);
	    for(genhep_vec::iterator it3=ggDaughters->begin();it3!=ggDaughters->end();it3++)
	      {
		int ggDaughterId=(*it3)->idhep();
		int signedGGDaughterId=(*it3)->idhep();
		Particle pp(**it3);
		if(abs(pp.lund())==PY_D0 || abs(pp.lund())==PY_D)
		  {
		    genhep_vec* gggDaughters=getDaughters(**it3);
		    for(genhep_vec::iterator it4=gggDaughters->begin();it4!=gggDaughters->end();it4++)
		      {
			int gggDaughterId=(*it4)->idhep();
			int signedGGGDaughterId=(*it4)->idhep();
			Particle ppp(**it4);
			if(print)
			  cout <<"   ->->D daughters: "<< ppp.pType().name()<<", id: " << (*it4)->get_ID()<<" lund: " << ppp.lund()<<" mom: " << ppp.p().rho()<<" ("<<ppp.p().px()<<", " << ppp.p().py()<< ", " << ppp.p().pz() <<", "<< ppp.p().t()<<")" <<endl;
			//the A0 seems to lead to a crash, probably because daughter particles are fishy..
			if(abs(ppp.lund())==213 || abs(ppp.lund())==113|| abs(ppp.lund())==10323 )
			  {
			    //					    cout <<"getting more daughters.." <<endl;
			    genhep_vec* ggggDaughters=getDaughters(**it4);
			    //					    cout <<" gggg done.." <<ggggDaughters->size()<<endl;

			    for(genhep_vec::iterator it5=ggggDaughters->begin();it5!=ggggDaughters->end();it5++)
			      {
				//						cout <<"-->1" <<endl;
				int ggggDaughterId=(*it5)->idhep();
				int signedGGGGDaughterId=(*it5)->idhep();

				Particle pppp(**it5);
				if(print)
				  cout <<"   ->->->D daughters: "<< pppp.pType().name()<<", id: " << (*it5)->get_ID()<<" lund: " << pppp.lund()<<" mom: " << pppp.p().rho()<< " ("<<pppp.p().px()<<", " << pppp.p().py()<< ", " << pppp.p().pz() <<", "<< pppp.p().t()<<")" <<endl;
			      }

			    delete ggggDaughters;
			  }
		      }

				  
		    delete gggDaughters;
		  }
		if(print)
		  cout <<"   ->D daughters: "<< pp.pType().name()<<", id: " << (*it3)->get_ID()<<" lund: " << pp.lund()<<" mom: " << pp.p().rho()<<" ("<<pp.p().px()<<", " << pp.p().py()<< ", " << pp.p().pz() <<", "<< pp.p().t()<<")" <<endl;
	      }
	    delete ggDaughters;
	  }

	//already found piplus...
	if(signedGDaughterId==PY_PI && foundPiPlus)
	  {
	    cout <<"found second piPlus!" <<endl;
	    foundSameChargePions=true;
	  }
	if(signedGDaughterId==PY_PI)
	  {
	    treeData.piPlusP=gdMom;
	    //				cout <<"pion mom: "<< (*it2)->PX() <<", " << (*it2)->PY() <<", " << (*it2)->PZ()<<endl;
	    //				cout <<"id: "<< (*it2)->get_ID()<<endl;
	    //				if(p.mdstVee()!=0)
	    //				  {
	    //  cout <<"found mdst Vee " <<endl;
	    //			

	    foundPiPlus=true;
	    continue;
	  }
	if(signedGDaughterId==-PY_PI && foundPiMinus)
	  {
	    cout <<"found second piMinus!" <<endl;
	    foundSameChargePions=true;
	  }
	if(signedGDaughterId==(-PY_PI))
	  {
	    treeData.piMinusP=gdMom;
	    foundPiMinus=true;
	    continue;
	  }
	if(print)
	  {
	    if( gDaughterId!=100423 && gDaughterId!=100421)
	      cout <<"gDaughter " << p.pType().name()<< ", id: "<< gDaughterId<<", ID: "<<(*it2)->get_ID() <<" mom: " << p.p().rho()<<endl;
	    else
	      cout <<"gDaughter: "<< gDaughterId <<endl;
	  }
	//D
	if(gDaughterId==PY_D || gDaughterId==PY_D0)
	  {
	    foundD=true;
	    continue;
	  }
	//D*
	if(PY_DStar==gDaughterId|| gDaughterId==  PY_DStar0)
	  {
	    foundDStar=true;
	    continue;
	  }
      }
    delete grandDaughters;
  }


  //difference to the other findDecaySignature is that we don't look at the decay of the D. We are only interested to find B to B Dlnu
  bool bToDDoubleStar::findDecaySignatureForBBRCorrection(const Gen_hepevt &mother,int& numLeptons, int& numPions, int& numKaons, int& numPi0, int& numBaryons, int& numNu,  int* br_Decays, HepLorentzVector& p_D, HepLorentzVector& p_l, HepLorentzVector& p_nu)
  {
    genhep_vec* daughters=getDaughters(mother);
    int lund=fabs(mother.idhep());
    if(lund==911|| lund>9000000)
      return false;
    Particle p(mother);
    //this is a d meson we have been looking for
    if(std::find(D_lundIds.begin(),D_lundIds.end(),lund)!=D_lundIds.end())
      {
	//get momentum
	p_D=p.p();
      }
    if(lund==PY_D0)
      {
	//cout <<"d0" <<endl;
	br_Decays[br_sig_D0]=1;
	lType=2;
	return true;
      }
    if(lund==PY_D)
      {
	//cout <<"d" <<endl;
	br_Decays[br_sig_D]=1;
	lType=2;
	return true;
      }
    if(lund==PY_DStar)
      {
	//cout <<"dstar" <<endl;
	br_Decays[br_sig_DStar]=1;
	lType=1;
	return true;
      }
    if(lund==PY_DStar0)
      {
	br_Decays[br_sig_DStar0]=1;
	lType=1;
	return true;
      }
    //
    if(lund==10413)
      {
	dssIdx=0;
	br_Decays[br_sig_D1]=1;
	return true;
      }
    if(PY_DStar_2==lund)
      {
	dssIdx=1;
	br_Decays[br_sig_D2]=1;
	return true;
      }
    if(lund==20413)
      {
	//I assume that this is the  D_1* from table 9 in bn1335...
	dssIdx=3;
	br_Decays[br_sig_D1Prime]=1;
	return true;
      }
    if(PY_DStar0Plus==lund)
      {
	dssIdx=2;
	br_Decays[br_sig_D0Star]=1;
	return true;
      }
    if(PY_D_10==lund)
      {
	dssIdx=0;
	br_Decays[br_sig_D10]=1;
	return true;
      }
    if(PY_DStar_20==lund)
      {
	dssIdx=1;
	br_Decays[br_sig_D20]=1;
	return true;
      }
    if(20423==lund)
      {
	dssIdx=3;
	br_Decays[br_sig_D1Prime0]=1;
	return true;
      }
    if(PY_DStar_00==lund)
      {
	dssIdx=2;
	br_Decays[br_sig_D0Star0]=1;
	return true;
      }
    ///////
    if(lund==PY_GAMMA)
      {
	//don't follow gamma. It can decay into e+e- which artificially increases the lepton count
	return true;
      }

    if(lund==PY_PI) 
      {
	numPions++;
	return true;
      }
    if(lund==PY_K)
      {
	numKaons++;
	return true;
      }

    if(lund==PY_Pi0)
      {
	numPi0++;
	return true;
      }
    if(lund==PY_E || lund==PY_Mu|| lund==PY_Tau)
      {
	numLeptons++;
	p_l=p.p();
	return true;
      }

    if(lund==PY_NuE || lund==PY_NuMu|| lund==PY_NuTau)
      {
	p_nu=p.p();
	numNu++;

	return true;
      }
    if(lund==2112 || lund==2212)
      {
	numBaryons++;
	return true;
      }
    if(daughters->size()<=0)
      {
	return true;
      }    
    for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
      {
	findDecaySignatureForBBRCorrection(**it,numLeptons,numPions,numKaons,numPi0,numBaryons,numNu,br_Decays,p_D,p_l,p_nu);
      }
  }

  //see if this B has any bDoubleSTar and count number of charged pions, kaons, pi0s (so independent of the actual decay)
  bool bToDDoubleStar::findDecaySignature(const Gen_hepevt &mother,bool& dDoubleStar,int& numLeptons, int& numPions, int& numKaons, int& numRhos,int& numPi0, int& numBaryons,int& numD, int& numDStar, int& numNu, bool& dStar_2S, bool& d_2S, int& dCharge, bool trackRhoDecay)
  {
    genhep_vec* daughters=getDaughters(mother);
    int lund=fabs(mother.idhep());
    if(lund==911|| lund>9000000)
      return false;
    Particle p(mother);
    //the ones with 20* are D'
    if(lund==PY_DStar_00|| lund==PY_DStar0Plus|| lund==PY_DStar_0s|| lund==PY_DStar0Plus|| lund==PY_D_1 || lund==PY_D_10 || lund==PY_DStar_2 || lund==PY_DStar_20 || lund==100413 || lund == 20423 ||lund ==20413|| lund==100421 || lund==100411 || lund==100423)
      {
	dDoubleStar=true;
	dDoubleStarId=lund;
	//// get the reconstructed mass of the D**, use all channels
	double recPx=0;
	double recPy=0;
	double recPz=0;
	double recE=0;
	for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
	  {
	    recPx+=(*it)->PX();
	    recPy+=(*it)->PY();
	    recPz+=(*it)->PZ();
	    recE+=(*it)->E();
	  }
	TLorentzVector tlv(recPx,recPy,recPz,recE);
	mcMassDDStar=tlv.M();
	if(dDoubleStarId==100423)
	  {
	    //	    cout <<"found the 2S with mass " << mcMassDDStar<<endl;
	  }
      }

    if(lund==100421|| lund==100411)
      d_2S=true;
    if(lund==100423|| lund==100413)
      dStar_2S=true;

    if(lund==PY_D || lund==PY_D0)
      {
	if(lund==PY_D0)
	  {
	    dCharge=0;
	  }
	else
	  {
	    dCharge=-1;
	    if(mother.idhep()>0)
	      {
		dCharge=1;
	      }
	  }
	numD++;
	return true;
      }
    if(lund==PY_RHO_0)
      {
	numRhos++;
	//determines if the pions we find from the rho decay are counting into the number of pions or not
	if(!trackRhoDecay)
	  return true;
      }
    if(lund==PY_DStar || lund==PY_DStar0)
      {
	if(lund==PY_DStar0)
	  {
	    dCharge=0;
	  }
	else
	  {
	    dCharge=-1;
	    if(mother.idhep()>0)
	      {
		dCharge=1;
	      }
	  }

	//don't add D decay products to number of pions etc...
	numDStar++;
	return true;
      }

    if(lund==PY_PI) 
      {
	///get mom, theta, phi 
    //	    Gen_hepevt hepEvt=get_hepevt(p->mdstCharged());
	mc_piMom.push_back(p.ptot());
	mc_piPhi.push_back(p.p3().theta());
	mc_piTheta.push_back(p.p3().phi());

	if(foundChargedDecIds.find(p.mdstCharged().get_ID())!=foundChargedDecIds.end())
	  mc_piFound.push_back(true);
	else
	  mc_piFound.push_back(false);


	numPions++;
	return true;
      }
    if(lund==PY_K)
      {
	numKaons++;
	return true;
      }
    if(lund==PY_GAMMA)
      {
	//don't follow gamma. It can decay into e+e- which artificially increases the lepton count
	return true;
      }
    if(lund==PY_Pi0)
      {
	numPi0++;
	return true;
      }
    if(lund==PY_E || lund==PY_Mu|| lund==PY_Tau)
      {
	numLeptons++;
	return true;
      }

    if(lund==PY_NuE || lund==PY_NuMu|| lund==PY_NuTau)
      {

	numNu++;

	return true;
      }
    if(lund==2112 || lund==2212)
      {
	numBaryons++;
	return true;
      }
    if(daughters->size()<=0)
      {
	return true;
      }    
    for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
      {
	findDecaySignature(**it,dDoubleStar,numLeptons,numPions,numKaons,numRhos,numPi0,numBaryons,numD, numDStar,numNu,dStar_2S, d_2S,dCharge,trackRhoDecay);
      }
  }



  genhep_vec* bToDDoubleStar::getDaughters(const Gen_hepevt &mother)
  {
    /* get a vector with all daughters as Gen_hepevt* */

    Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();

    int n_children = mother.daLast() - mother.daFirst() + 1;

    genhep_vec *children = new genhep_vec();

    for(int i=0; i<n_children; i++) {

      Panther_ID ID0(mother.daFirst()+i);
      if(ID0==0)
	{
	  //	  if(PRINT)
	  //	    cout <<"wrong!!!" <<endl;
	  break;
	}
      Gen_hepevt& temp = gen_hepevt_mgr(ID0);

      if (temp) 
	{
	  children->push_back(&temp);
	}
    }

    return children;
  }

  bool bToDDoubleStar::isD0(const Gen_hepevt& gen_it)
  {
    if(foundDPiPi)
      {
	//	    cout <<"checking for D0 decay" <<endl;
      }
    genhep_vec* daughters=getDaughters(gen_it);
    if(daughters->size()==0)
      {
	delete daughters;
	return false;
      }
    else
      {
	bool foundPi0=false;

	int foundKPlus=0;
	int foundPiPlus=0;
	int foundKMinus=0;
	int foundPiMinus=0;
	bool foundKs=false;

	for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
	  {
	    int daughterId=(*it)->idhep();
	    //	    cout <<"decay id: " << daughterId<<endl;
	    switch(daughterId)
	      {
	      case PY_Pi0:
		foundPi0=true;
		break;
	      case PY_PI:
		foundPiPlus++;
		break;
	      case -PY_PI:
		foundPiMinus++;
		break;
	      case PY_K:
		foundKPlus++;
		break;
	      case -PY_K:
		foundKMinus++;
		break;
	      case PY_KS0:
		foundKs=true;
		break;
	      default:
		//	      if(foundDPiPi)
		//	      {
		//		cout <<"found other decay product: "<< daughterId<<endl;
		//}
		//if we found some other decay produce, return
		delete daughters;
		return false;
	      }
	  }

	bool ret=false;
	if(foundKMinus==1 && foundPiPlus==1&& daughters->size()==2)
	  ret=true;
	if(foundKMinus==1&& foundPiPlus==1 && foundPi0 && daughters->size()==3)
	  ret=true;
	if(foundKMinus==1 && foundPiPlus==2 && foundPiMinus==1 && daughters->size()==4)
	  ret=true;
	if(foundKs&& foundPiPlus==1 && foundPiMinus==1 && daughters->size()==3)
	  ret=true;
	if(foundKPlus==1 && foundKMinus==1 && daughters->size()==2)
	  ret=true;
	if(foundKs && foundPi0 && daughters->size()==2)
	  ret=true;

	if(foundDPiPi)
	  {
	    //	  cout << "no good decay ..." <<endl;
	  }
	delete daughters;
	return ret;
      }
  }
  bool bToDDoubleStar::isChargedD(const Gen_hepevt& gen_it)
  {
    int sign=1;
    if(gen_it.idhep()<0)
      sign=-1;
    //for the D- decays we flip the sign of the daughters
    genhep_vec* daughters=getDaughters(gen_it);
    if(daughters->size()==0)
      {
	delete daughters;
	return false;
      }
    else
      {
	bool foundPi0=false;

	int foundKPlus=0;
	int foundPiPlus=0;
	int foundKMinus=0;
	int foundPiMinus=0;
	bool foundKs=false;

	for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
	  {
	    //change the sign according to D+ (I.e. flip for D-)
	    int daughterId=sign*(*it)->idhep();

	    switch(daughterId)
	      {
	      case PY_Pi0:
		foundPi0=true;
		break;
	      case PY_PI:
		foundPiPlus++;
		break;
	      case -PY_PI:
		foundPiMinus++;
		break;
	      case PY_K:
		foundKPlus++;
		break;
	      case -PY_K:
		foundKMinus++;
		break;
	      case PY_KS0:
		foundKs=true;
		break;
	      default:
		//if we found some other decay produce, return
		delete daughters;
		return false;
	      }
	  }
	if(foundKs && foundPiPlus==1 && daughters->size()==2)
	  {
	    delete daughters;
	    return true;
	  }
	if(foundKs&& foundPiPlus==2 && foundPiMinus==1 && daughters->size()==4)
	  {
	    delete daughters;
	    return true;
	  }
	if(foundKMinus==1 && foundPiPlus==2 && daughters->size()==3)
	  {
	    delete daughters;
	    return true;
	  }

	if(foundKPlus==1 && foundKMinus==1  && foundPiPlus==1&& daughters->size()==3)
	  {
	    delete daughters;
	    return true;
	  }
	delete daughters;
	return false;
      }
  }
  bool bToDDoubleStar::isDStar(const Gen_hepevt& gen_it)
  {

    //only look at D* --> pi0 D0
    genhep_vec* daughters=getDaughters(gen_it);
    if(daughters->size()==0)
      {
	delete daughters;
	return false;
      }

    bool ret=false;
    if(daughters->size()==2)
      {
	if((*daughters)[0]->idhep()==PY_Pi0 && (*daughters)[1]->idhep()==PY_D0 &&  gen_it.idhep()==PY_DStar0)
	  ret=true;
	if((*daughters)[1]->idhep()==PY_Pi0 && (*daughters)[0]->idhep()==PY_D0&&  gen_it.idhep()==PY_DStar0)
	  ret=true;

	///D+

	if((*daughters)[0]->idhep()==PY_Pi0 && (*daughters)[1]->idhep()==PY_D&&  gen_it.idhep()==PY_DStar)
	  ret=true;
	if((*daughters)[1]->idhep()==PY_Pi0 && (*daughters)[0]->idhep()==PY_D&&  gen_it.idhep()==PY_DStar)
	  ret=true;

	//charged pion , D0
	if((*daughters)[0]->idhep()==PY_PI && (*daughters)[1]->idhep()==PY_D0&&  gen_it.idhep()==PY_DStar)
	  ret=true;
	if((*daughters)[1]->idhep()==PY_PI && (*daughters)[0]->idhep()==PY_D0&&  gen_it.idhep()==PY_DStar)
	  ret=true;


	//D-
	if((*daughters)[0]->idhep()==PY_Pi0 && (*daughters)[1]->idhep()==-PY_D&&  gen_it.idhep()==-PY_DStar)
	  ret=true;
	if((*daughters)[1]->idhep()==PY_Pi0 && (*daughters)[0]->idhep()==-PY_D&&  gen_it.idhep()==-PY_DStar)
	  ret=true;

	//charged pion , D0
	if((*daughters)[0]->idhep()==-PY_PI && (*daughters)[1]->idhep()==PY_D0&&  gen_it.idhep()==-PY_DStar)
	  ret=true;
	if((*daughters)[1]->idhep()==-PY_PI && (*daughters)[0]->idhep()==PY_D0&&  gen_it.idhep()==-PY_DStar)
	  ret=true;

      }
    delete daughters;
    return ret;

  }

  //reconstruct D0 from pion, kaon, ks..
  void bToDDoubleStar::reconstructD0()
  {
    double m_D0=1.86484;
    double m_d0mass_max=m_D0+0.015;
    double m_d0mass_min=m_D0-0.015;


    //the mass range for which we do a mass/vertex constrained fit
    double m_d0mass_maxLoose=m_D0+3*0.015;
    double m_d0mass_minLoose=m_D0-3*0.015;



    //nominal mass: 1.86484, use +- 0.015, except for K-pi+pi0 (0.025)
    //D-->Kpi
    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
      {
	for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
	  {
	    Particle& pion= *(*itP);
	    Particle& kaon= *(*itK);
	    if(kaon.charge()+pion.charge()!=0) continue;
	    HepLorentzVector p_d0=pion.p()+kaon.p();
	    double m=p_d0.mag();
	    histoRecD0Spect->Fill(m);
	    dMass=m;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    //	    cout <<"1 filling with " << m <<endl;
	    dType=0;
	    dDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m>1.7 && m<2.0)
		  {
		    if(m_mc)
		      {
			//		Gen_hepevt hepEvt=get_hepevt(leptonCandidates[0]->mdstCharged());
			const Mdst_charged &h_kCharged=kaon.mdstCharged();
			const Mdst_charged &h_pCharged=pion.mdstCharged();
			const Gen_hepevt &h_1=get_hepevt(h_kCharged);
			const Gen_hepevt &h_2=get_hepevt(h_pCharged);
			int tempSig=0;
			if(h_1 && h_2 && commonParent(&h_1,&h_2,0,0,0))
			  {
			    tempSig=1;
			  }
			
			dIsSignal=0;
			if(h_1 && h_2 && h_1.mother() && h_2.mother()  && h_1.mother().get_ID()==h_2.mother().get_ID() ){
			  dIsSignal=1;
			}
			if(dIsSignal!=tempSig)
			  {
			    //			    			    cout <<"!!!"<<endl<<endl;
			    //			    			    cout<<"mismatch between parent finding!! : "<< tempSig<<" old metod: "<< dIsSignal<<endl;
			    //			      cout <<"!!!"<<endl<<endl;
			    //			    //			    			      recursivePrint(h_1.mother().mother(),"");
			    //			      cout <<"second " << endl;
			    //			    			      recursivePrint(h_2.mother().mother(),"");
			    //						      cout <<"done " <<endl<<endl;

			  }
		      }
		    mesonMassTree->Fill();
		  }

	      }
	  
	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
	    Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));
	    d0->userInfo(*(new ParticleInfoDecay()));
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(d0->userInfo());
	    pinf.decayChannel=dDecay;
	    pinf.type=dType;
	    pinf.pdgDiff=m-m_D0;
	    pinf.mRec=m;
	    pinf.mDRec=m;

	    //putddecay
	    //	if(!doKmVtxFit2(*(*itD),  confLevel,0))

	    d0->relation().append(kaon);
	    d0->relation().append(pion);
	    if(m_mc)
	      {

		const Mdst_charged &h_kCharged=kaon.mdstCharged();
		const Mdst_charged &h_pCharged=pion.mdstCharged();
		const Gen_hepevt &h_kaon=get_hepevt(h_kCharged);
		const Gen_hepevt &h_pion=get_hepevt(h_pCharged);
		
	
		if(h_kaon && h_pion && h_kaon.mother() && h_pion.mother() && h_kaon.mother().get_ID()==h_pion.mother().get_ID()){
		  d0->relation().genHepevt(h_kaon.mother());
		}
	      } 
	    D0Candidates.push_back(d0);
	  }
      }

    ////--->D to K-pi+pi0
    //uses different mass cut:
    m_d0mass_max=m_D0+0.025;
    m_d0mass_min=m_D0-0.025;
    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
      {
	for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
	  {
	    for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
	      {
		Particle& pi0=*(*itPi0);
		Particle& pion= *(*itP);
		Particle& kaon= *(*itK);
		if(kaon.charge()+pion.charge()!=0) continue;
		HepLorentzVector p_d0=pion.p()+kaon.p()+pi0.p();

		//result of check below: yes, there are mdst set
		//		cout <<"checking mdst for pi0: ";
		//		if(pi0.mdstPi0())
		//		  cout <<"yes " <<endl;
		//		else
		//		  cout <<"no " <<endl;
		
		double m=p_d0.mag();

		histoRecD0Spect->Fill(m);
		dMass=m;
		//	    cout <<"2 filling with " << m <<endl;
		dType=0;
		dDecay=1;
		if(foundDPiPi || foundSinglePionDecay)
		  foundDDoubleStarDecay=1;
		else
		  foundDDoubleStarDecay=0;
		if(SAVE_MESON_MASS_DISTRIBUTIONS)
		  {
		    if(m>1.7 && m<2.0)
		      {
	
			if(m_mc)
			  {

			    const Mdst_charged &h_c1=kaon.mdstCharged();
			    const Mdst_charged &h_c2=pion.mdstCharged();
			    const Mdst_pi0 &h_c3=pi0.mdstPi0();
			    const Gen_hepevt &h_1=get_hepevt(h_c1);
			    const Gen_hepevt &h_2=get_hepevt(h_c2);
			    const Gen_hepevt &h_3=get_hepevt(h_c3);
	
			    dIsSignal=0;
			    if(h_1 && h_2 && h_3)
			      {
				if(commonParent(&h_1,&h_2,&h_3,0,0)){
				  dIsSignal=1;
				}
			      }
			  }
			mesonMassTree->Fill();
		      }
		  }

	      

		if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
		Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));
		d0->userInfo(*(new ParticleInfoDecay()));
		ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(d0->userInfo());
		pinf.decayChannel=dDecay;
		pinf.type=dType;
		pinf.pdgDiff=m-m_D0;
		pinf.mRec=m;
		pinf.mDRec=m;
		d0->relation().append(kaon);
		d0->relation().append(pion);
		d0->relation().append(pi0);
		if(m_mc)
		  {
		    const Mdst_charged &h_c1=kaon.mdstCharged();
		    const Mdst_charged &h_c2=pion.mdstCharged();
		    const Mdst_pi0 &h_c3=pi0.mdstPi0();
		    const Gen_hepevt &h_kaon=get_hepevt(h_c1);
		    const Gen_hepevt &h_pion=get_hepevt(h_c2);
		    const Gen_hepevt &h_pi0=get_hepevt(h_c3);


		    //		    if(h_pi0&&h_kaon && h_pion && h_kaon.mother() && h_pion.mother() && h_kaon.mother().get_ID()==h_pion.mother().get_ID() && h_pi0.mother() && h_pi0.mother().get_ID()==h_kaon.mother().get_ID()){
		    if(h_kaon && h_pion && h_pi0&&commonParent(&h_kaon,&h_pion,&h_pi0,0,0)){
		      d0->relation().genHepevt(h_kaon.mother());
		    }
		  }
		D0Candidates.push_back(d0);
	      }
	  }
      }
    //and set the mass cuts back..
    m_d0mass_max=m_D0+0.015;
    m_d0mass_min=m_D0-0.015;
    ///D->K-pi+pi+pi-
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=(itP1+1);itP2!=chargedPiCandidates.end();itP2++)  
	  {
	    for(vector<Particle*>::iterator itP3=(itP2+1);itP3!=chargedPiCandidates.end();itP3++)  
	      {
		for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
		  {
		    Particle& pion1= *(*itP1);
		    Particle& pion2= *(*itP2);
		    Particle& pion3= *(*itP3);
		    Particle& kaon= *(*itK);
		    if(kaon.charge()+pion1.charge()+pion2.charge()+pion3.charge()!=0) continue;
		    HepLorentzVector p_d0=kaon.p()+pion1.p()+pion2.p()+pion3.p();
		    double m=p_d0.mag();
		    dMass=m;
		    //	    cout <<"3 filling with " << m <<endl;
		    dType=0;
		    dDecay=2;
		    if(foundDPiPi || foundSinglePionDecay)
		      foundDDoubleStarDecay=1;
		    else
		      foundDDoubleStarDecay=0;
		    if(SAVE_MESON_MASS_DISTRIBUTIONS)
		      {
			if(m>1.7 && m<2.0)
			  {
	
			    if(m_mc)
			      {

				const Mdst_charged &h_c1=pion1.mdstCharged();
				const Mdst_charged &h_c2=pion2.mdstCharged();
				const Mdst_charged &h_c3=pion3.mdstCharged();
				const Mdst_charged &h_c4=kaon.mdstCharged();
				const Gen_hepevt &h_1=get_hepevt(h_c1);
				const Gen_hepevt &h_2=get_hepevt(h_c2);
				const Gen_hepevt &h_3=get_hepevt(h_c3);
				const Gen_hepevt &h_4=get_hepevt(h_c4);


				dIsSignal=0;
				if(h_1 && h_2 &&h_3&& h_4&&  h_1.mother() && h_2.mother() && h_3.mother() && h_4.mother()  && h_1.mother().get_ID()==h_2.mother().get_ID() && h_1.mother().get_ID()==h_3.mother().get_ID() && h_1.mother().get_ID()==h_4.mother().get_ID()){
				  dIsSignal=1;
				}
			      } 
			    mesonMassTree->Fill();
			  }
		      }


		      
		    histoRecD0Spect->Fill(m);
		    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
		    Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));
		    d0->userInfo(*(new ParticleInfoDecay()));
		    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(d0->userInfo());
		    pinf.decayChannel=dDecay;
		    pinf.type=dType;
		    pinf.pdgDiff=m-m_D0;
		    pinf.mRec=m;
		    pinf.mDRec=m;
		    d0->relation().append(kaon);
		    d0->relation().append(pion1);
		    d0->relation().append(pion2);
		    d0->relation().append(pion3);
		    if(m_mc)
		      {


			const Mdst_charged &h_c1=pion1.mdstCharged();
			const Mdst_charged &h_c2=pion2.mdstCharged();
			const Mdst_charged &h_c3=pion3.mdstCharged();
			const Mdst_charged &h_c4=kaon.mdstCharged();
			const Gen_hepevt &h_pion1=get_hepevt(h_c1);
			const Gen_hepevt &h_pion2=get_hepevt(h_c2);
			const Gen_hepevt &h_pion3=get_hepevt(h_c3);
			const Gen_hepevt &h_kaon=get_hepevt(h_c4);


		    
			if(h_kaon && h_pion1 && h_pion2&&h_pion3){
			  if(h_kaon.mother() && h_pion1.mother() &&h_pion2.mother() && h_pion3.mother()){
			    if( h_kaon.mother().get_ID()==h_pion1.mother().get_ID() && h_pion1.mother().get_ID()==h_pion2.mother().get_ID() && h_pion2.mother().get_ID()==h_pion2.mother().get_ID()){
			      d0->relation().genHepevt(h_kaon.mother());
			    }
			  }
			}
		      }
		    D0Candidates.push_back(d0);
		  }
	      }
	  }
      }
 
    ////D-->Ks pi+pi-
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=itP1+1;itP2!=chargedPiCandidates.end();itP2++)
	  {
	    for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
	      {
		Particle& pion1=*(*itP1);
		Particle& pion2= *(*itP2);
		Particle& Ks= *(*itKs);
		if((pion1.charge()+pion2.charge())!=0) continue;
		HepLorentzVector p_d0=pion1.p()+pion2.p()+Ks.p();
		double m=p_d0.mag();
		dMass=m;
		//	    cout <<"4 filling with " << m <<endl;
		dType=0;
		dDecay=3;
		if(foundDPiPi || foundSinglePionDecay)
		  foundDDoubleStarDecay=1;
		else
		  foundDDoubleStarDecay=0;
		if(SAVE_MESON_MASS_DISTRIBUTIONS)
		  {
		    if(m>1.85 && m<1.875)
		      {
			if(m_mc)
			  {
			    const Mdst_charged &h_c1=pion1.mdstCharged();
			    const Mdst_charged &h_c2=pion2.mdstCharged();
			    const Mdst_vee2 &h_c3=Ks.mdstVee2();
			    if(h_c3)
			      {
				cout <<"found ks mc particle "<< endl;
			      }
			    const Gen_hepevt &h_1=get_hepevt(h_c1);
			    const Gen_hepevt &h_2=get_hepevt(h_c2);
			    const Gen_hepevt &h_3=get_hepevt(h_c3);


			    dIsSignal=0;
			    //			    cout <<"trying the common parent thing  " << endl;
			    if(h_1 && h_2 && h_3)
			      {
				if(commonParent(&h_1,&h_2,&h_3,0,0))
				  {
				    //				cout <<"signal is !!!" <<endl;
				    dIsSignal=1;
				  }
			      }
			  }
			mesonMassTree->Fill();
		      }
		  }

		  
		histoRecD0Spect->Fill(m);
		if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));
		Particle* d0 =new Particle(p_d0,Ptype("D0"));
		d0->userInfo(*(new ParticleInfoDecay()));
		ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(d0->userInfo());
		pinf.decayChannel=dDecay;
		pinf.type=dType;
		pinf.pdgDiff=m-m_D0;
		pinf.mDRec=m;
		pinf.mRec=m;
		d0->relation().append(pion1);
		d0->relation().append(pion2);
		d0->relation().append(Ks);
		if(m_mc)
		  {

		    const Mdst_charged &h_c1=pion1.mdstCharged();
		    const Mdst_charged &h_c2=pion2.mdstCharged();
		    const Mdst_vee2 &h_c3=Ks.mdstVee2();

		    const Gen_hepevt &h_pion1=get_hepevt(h_c1);
		    const Gen_hepevt &h_pion2=get_hepevt(h_c2);
		    const Gen_hepevt &h_Ks=get_hepevt(h_c3);

		    //		    if(h_pion1&&h_pion2 && h_Ks&& h_pion1.mother() && h_pion2.mother() && h_Ks.mother() && h_pion1.mother().get_ID()==h_pion2.mother().get_ID() && h_Ks.mother().mother()&& h_pion1.mother().get_ID()==h_Ks.mother().mother().get_ID()){
		    if(h_pion1&& h_pion2&& h_Ks&&commonParent(&h_pion1,&h_pion2,&h_Ks,0,0)){
		      d0->relation().genHepevt(h_pion1.mother());
		    }
		  }
		D0Candidates.push_back(d0);
	      }
	  }
      }

    ////D-->K+K-
    for(vector<Particle*>::iterator itK1=chargedKCandidates.begin();itK1!=chargedKCandidates.end();itK1++)
      {
	for(vector<Particle*>::iterator itK2=itK1+1;itK2!=chargedKCandidates.end();itK2++)
	  {
	    Particle& kaon1= *(*itK1);
	    Particle& kaon2= *(*itK2);
	    if(kaon1.charge()+kaon2.charge()!=0) continue;
	    HepLorentzVector p_d0=kaon1.p()+kaon2.p();
	    double m=p_d0.mag();
	    dMass=m;
	    //	    cout <<"5 filling with " << m <<endl;
	    dType=0;
	    dDecay=4;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m>1.7 && m<2.0)
		  {
	
		    if(m_mc)
		      {
			const Mdst_charged &h_c1=kaon1.mdstCharged();
			const Mdst_charged &h_c2=kaon2.mdstCharged();
			const Gen_hepevt &h_1=get_hepevt(h_c1);
			const Gen_hepevt &h_2=get_hepevt(h_c2);


			dIsSignal=0;
			if(h_1 && h_2 && h_1.mother() && h_2.mother()  && h_1.mother().get_ID()==h_2.mother().get_ID()){
			  dIsSignal=1;
			}
		      } 
		    mesonMassTree->Fill();
		  }

	      }
	    
	    histoRecD0Spect->Fill(m);
	    //	    cout <<"found k/k combination, filling with m: "<< m <<endl;
	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
	    Particle* d0 =new Particle(p_d0,Ptype("D0"));
	    d0->userInfo(*(new ParticleInfoDecay()));
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(d0->userInfo());
	    pinf.decayChannel=dDecay;
	    pinf.type=dType;
	    pinf.pdgDiff=m-m_D0;
	    pinf.mDRec=m;
	    pinf.mRec=m;
	    d0->relation().append(kaon1);
	    d0->relation().append(kaon2);
	    if(m_mc)
	      {

		const Mdst_charged &h_c1=kaon1.mdstCharged();
		const Mdst_charged &h_c2=kaon2.mdstCharged();
		const Gen_hepevt &h_kaon1=get_hepevt(h_c1);
		const Gen_hepevt &h_kaon2=get_hepevt(h_c2);


		if(h_kaon1 && h_kaon2 && h_kaon1.mother() && h_kaon2.mother() && h_kaon1.mother().get_ID()==h_kaon2.mother().get_ID()){
		  d0->relation().genHepevt(h_kaon1.mother());
		}
	      } 
	    D0Candidates.push_back(d0);
	  }
      }

    //D-->Kspi0
    for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
      {
	for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
	  {
	    Particle& pi0= *(*itPi0);
	    Particle& Ks= *(*itKs);
	    HepLorentzVector p_d0=pi0.p()+Ks.p();
	    double m=p_d0.mag();
	    dMass=m;
	    //	    cout <<"6 filling with " << m <<endl;
	    dType=0;
	    dDecay=5;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m>1.7 && m<2.0)
		  {
		    if(m_mc)
		      {
			
			const Mdst_pi0 &h_c1=pi0.mdstPi0();
			const Mdst_vee2 &h_c2=Ks.mdstVee2();
			const Gen_hepevt &h_1=get_hepevt(h_c1);
			const Gen_hepevt &h_2=get_hepevt(h_c2);


			dIsSignal=0;
			//			if(h_1 && h_2 && h_1.mother() && h_2.mother()  &&h_2.mother().mother()&& h_1.mother().get_ID()==h_2.mother().mother().get_ID()){
			if(h_1 && h_2)
			  {
			    if(commonParent(&h_1,&h_2,0,0,0)){
			      dIsSignal=1;
			    }
			  }
		      } 
		    mesonMassTree->Fill();
		  }
	      }

	      
	    histoRecD0Spect->Fill(m);
	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
	    Particle* d0 =new Particle(p_d0,Ptype( "D0"));
	    d0->userInfo(*(new ParticleInfoDecay()));
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(d0->userInfo());
	    pinf.decayChannel=dDecay;
	    pinf.type=dType;
	    pinf.pdgDiff=m-m_D0;
	    pinf.mDRec=m;
	    pinf.mRec=m;
	    d0->relation().append(pi0);
	    d0->relation().append(Ks);
	    if(m_mc)
	      {
		const Mdst_pi0 &h_c1=pi0.mdstPi0();
		const Gen_hepevt &h_pi0=get_hepevt(h_c1);
		const Mdst_vee2 &h_c2=Ks.mdstVee2();
		
		const Gen_hepevt &h_Ks=get_hepevt(h_c2);

		//		if(h_pi0 && h_Ks && h_pi0.mother() && h_Ks.mother() && h_Ks.mother().mother()&&h_pi0.mother().get_ID()==h_Ks.mother().mother().get_ID() )
		if(h_pi0&& h_Ks &&commonParent(&h_pi0,&h_Ks,0,0,0))
		  {
		    d0->relation().genHepevt(h_pi0.mother());
		  }
	      } 
	    D0Candidates.push_back(d0);
	  }
      }
  }

  void bToDDoubleStar::reconstructChargedD()
  {
    double m_DPlus=1.86962;
    double m_dPlusmass_max=m_DPlus+0.015;
    double m_dPlusmass_min=m_DPlus-0.015;
    ////D-->KsPi
    for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
      {
	for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	  {
	    Particle& Ks= *(*itKs);
	    Particle& pion= *(*itP);
	    HepLorentzVector p_dPlus=Ks.p()+pion.p();
	    double m=p_dPlus.mag();
	    dMass=m;
	    //	    cout <<"7 filling with " << m <<endl;
	    dType=1;
	    dDecay=0;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m > 1.7 && m<2.0)
		  {
		    if(m_mc)
		      {

			const Mdst_charged &h_c1=pion.mdstCharged();
			const Mdst_vee2 &h_c2=Ks.mdstVee2();
			const Gen_hepevt &h_1=get_hepevt(h_c1);
			const Gen_hepevt &h_2=get_hepevt(h_c2);

			dIsSignal=0;
			//			if(h_1 && h_2 && h_1.mother() && h_2.mother() &&h_2.mother().mother() && h_1.mother().get_ID()==h_2.mother().mother().get_ID()){
			if(h_1 && h_2)
			  {
			    if(commonParent(&h_1,&h_2,0,0,0)){
			      dIsSignal=1;
			    }
			  }
			mesonMassTree->Fill();

		      }
		  }	  
	      }
	      
	    histoRecDSpect->Fill(m);
	      
	    if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
	    Particle* dPlus =new Particle(p_dPlus,Ptype(pion.charge() > 0 ? "D+" : "D-"));
	    dPlus->userInfo(*(new ParticleInfoDecay()));
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(dPlus->userInfo());
	    pinf.decayChannel=dDecay;
	    pinf.type=dType;
	    pinf.pdgDiff=m-m_DPlus;
	    pinf.mDRec=m;
	    pinf.mRec=m;
	    dPlus->relation().append(Ks);
	    dPlus->relation().append(pion);
	    if(m_mc)
	      {


		const Mdst_charged &h_c1=pion.mdstCharged();
		const Mdst_vee2 &h_c2=Ks.mdstVee2();
		const Gen_hepevt &h_pion=get_hepevt(h_c1);
		const Gen_hepevt &h_Ks=get_hepevt(h_c2);

		//		if(h_Ks && h_pion && h_Ks.mother() && h_Ks.mother().mother() && h_pion.mother() && h_Ks.mother().mother().get_ID()==h_pion.mother().get_ID()){
		if(h_pion && h_Ks && commonParent( &h_pion,&h_Ks,0,0,0)){
		  dPlus->relation().genHepevt(h_pion.mother());
		}
	      } 
	    chargedDCandidates.push_back(dPlus);
	  }
      }
    

    ///D->Kspi+pi+pi-
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=(itP1+1);itP2!=chargedPiCandidates.end();itP2++)  
	  {
	    for(vector<Particle*>::iterator itP3=(itP2+1);itP3!=chargedPiCandidates.end();itP3++)  
	      {
		for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
		  {
		    Particle& pion1= *(*itP1);
		    Particle& pion2= *(*itP2);
		    Particle& pion3= *(*itP3);
		    Particle& Ks= *(*itKs);
		    double charge =pion1.charge()+pion2.charge()+pion3.charge();
		    if(fabs(pion1.charge()+pion2.charge()+pion3.charge())>1) continue;
		    HepLorentzVector p_dPlus=Ks.p()+pion1.p()+pion2.p()+pion3.p();
		    double m=p_dPlus.mag();
		    dMass=m;
		    //	    cout <<"8 filling with " << m <<endl;
		    dType=1;
		    dDecay=1;
		    if(foundDPiPi || foundSinglePionDecay)
		      foundDDoubleStarDecay=1;
		    else
		      foundDDoubleStarDecay=0;
		    if(SAVE_MESON_MASS_DISTRIBUTIONS)
		      {
			if(m>1.7 && m < 2.0)
			  {
			    if(m_mc)
			      {

				const Mdst_charged &h_c1=pion1.mdstCharged();
				const Mdst_charged &h_c2=pion2.mdstCharged();
				const Mdst_charged &h_c3=pion3.mdstCharged();
				const Mdst_vee2 &h_c4=Ks.mdstVee2();
				const Gen_hepevt &h_1=get_hepevt(h_c1);
				const Gen_hepevt &h_2=get_hepevt(h_c2);
				const Gen_hepevt &h_3=get_hepevt(h_c3);
				const Gen_hepevt &h_4=get_hepevt(h_c4);



				dIsSignal=0;
				//			if(h_1 && h_2 && h_3  && h_4 && h_1.mother() && h_2.mother()  && h_3.mother() &&h_4.mother() && h_4.mother().mother()&& h_1.mother().get_ID()==h_2.mother().get_ID() && h_1.mother().get_ID()==h_3.mother().get_ID() && h_3.mother().get_ID()==h_4.mother().mother().get_ID()){
				//			    cout <<"trying the common parent thing  again " << endl;
				if(h_1 && h_2 && h_3&& h_4)
				  {
				    if(commonParent(&h_1,&h_2,&h_3,&h_4,0)){
				      dIsSignal=1;
				    }
				  }
			      } 
			    mesonMassTree->Fill();
			  }
		      }
		      

		      
		    histoRecDSpect->Fill(m);
		    if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		    Particle* dPlus =new Particle(p_dPlus,Ptype(charge >0 ? "D+" : "D-"));
		    dPlus->userInfo(*(new ParticleInfoDecay()));
		    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(dPlus->userInfo());
		    pinf.decayChannel=dDecay;
		    pinf.type=dType;
		    pinf.pdgDiff=m-m_DPlus;
		    pinf.mDRec=m;
		    pinf.mRec=m;

		    dPlus->relation().append(Ks);
		    dPlus->relation().append(pion1);
		    dPlus->relation().append(pion2);
		    dPlus->relation().append(pion3);
		    if(m_mc)
		      {

			const Mdst_charged &h_c1=pion1.mdstCharged();
			const Mdst_charged &h_c2=pion2.mdstCharged();
			const Mdst_charged &h_c3=pion3.mdstCharged();
			const Mdst_vee2 &h_c4=Ks.mdstVee2();
			const Gen_hepevt &h_pion1=get_hepevt(h_c1);
			const Gen_hepevt &h_pion2=get_hepevt(h_c2);
			const Gen_hepevt &h_pion3=get_hepevt(h_c3);
			const Gen_hepevt &h_Ks=get_hepevt(h_c4);

		    
			//			if(h_Ks && h_pion1 && h_pion2&&h_pion3){
			//			  if(h_Ks.mother() && h_pion1.mother() &&h_pion2.mother() && h_pion3.mother() && h_Ks.mother().mother()){
			//			    if( h_Ks.mother().mother().get_ID()==h_pion1.mother().get_ID() && h_pion1.mother()==h_pion2.mother().get_ID() && h_pion2.mother().get_ID()==h_pion2.mother().get_ID()){
			if(h_pion1&& h_pion2&& h_pion3 && h_Ks&& commonParent(&h_pion1,&h_pion2,&h_pion3,&h_Ks,0)){
			  dPlus->relation().genHepevt(h_pion1.mother());
			}
		      }
		
		    chargedDCandidates.push_back(dPlus);
		  }
	      }
	  }
      }


    //D-->K-pi+pi+
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=itP1+1;itP2!=chargedPiCandidates.end();itP2++)
	  {
	    for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
	      {
		Particle& pion1=*(*itP1);
		Particle& pion2= *(*itP2);
		Particle& kaon= *(*itK);
		
		if(pion1.charge()!=pion2.charge()) continue;
		if(kaon.charge()==pion1.charge()) continue;
		double charge = pion1.charge()+pion2.charge()+kaon.charge();
		HepLorentzVector p_dPlus=pion1.p()+pion2.p()+kaon.p();
		double m=p_dPlus.mag();
		dMass=m;
		//	    cout <<"9 filling with " << m <<endl;
		dType=1;
		dDecay=2;
		if(foundDPiPi || foundSinglePionDecay)
		  foundDDoubleStarDecay=1;
		else
		  foundDDoubleStarDecay=0;
		if(SAVE_MESON_MASS_DISTRIBUTIONS)
		  {
		    if(m>1.7 && m < 2.0)
		      {
			if(m_mc)
			  {
			    const Mdst_charged &h_c1=pion1.mdstCharged();
			    const Mdst_charged &h_c2=pion2.mdstCharged();
			    const Mdst_charged &h_c3=kaon.mdstCharged();

			    const Gen_hepevt &h_1=get_hepevt(h_c1);
			    const Gen_hepevt &h_2=get_hepevt(h_c2);
			    const Gen_hepevt &h_3=get_hepevt(h_c3);



			    dIsSignal=0;
			    if(h_1 && h_2 && h_3 && h_1.mother() && h_2.mother()  && h_3.mother() && h_1.mother().get_ID()==h_2.mother().get_ID() && h_1.mother().get_ID()==h_3.mother().get_ID()){
			      dIsSignal=1;
			    }
			  }
			mesonMassTree->Fill();
		      }
		  }

		  
		histoRecDSpect->Fill(m);
		if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));
		Particle* dPlus =new Particle(p_dPlus,Ptype(charge > 0 ? "D+": "D-"));
		dPlus->userInfo(*(new ParticleInfoDecay()));
		ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(dPlus->userInfo());
		pinf.decayChannel=dDecay;
		pinf.type=dType;
		pinf.pdgDiff=m-m_DPlus;
		pinf.mDRec=m;
		pinf.mRec=m;

		dPlus->relation().append(pion1);
		dPlus->relation().append(pion2);
		dPlus->relation().append(kaon);
		if(m_mc)
		  {


		    const Mdst_charged &h_c1=pion1.mdstCharged();
		    const Mdst_charged &h_c2=pion2.mdstCharged();
		    const Mdst_charged &h_c3=kaon.mdstCharged();

		    const Gen_hepevt &h_pion1=get_hepevt(h_c1);
		    const Gen_hepevt &h_pion2=get_hepevt(h_c2);
		    const Gen_hepevt &h_kaon=get_hepevt(h_c3);


		    if(h_pion1&&h_pion2 && h_kaon&& h_pion1.mother() && h_pion2.mother() && h_kaon.mother() && h_pion1.mother().get_ID()==h_pion2.mother().get_ID() && h_pion1.mother().get_ID()==h_kaon.mother().get_ID()){
		      dPlus->relation().genHepevt(h_pion1.mother());
		    }
		  }
		chargedDCandidates.push_back(dPlus);
	      }
	  }
      }

    //D-->K+K-pi+
    for(vector<Particle*>::iterator itK1=chargedKCandidates.begin();itK1!=chargedKCandidates.end();itK1++)
      {
	for(vector<Particle*>::iterator itK2=itK1+1;itK2!=chargedKCandidates.end();itK2++)
	  {
	    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	      {
		Particle& kaon1=*(*itK1);
		Particle& kaon2= *(*itK2);
		Particle& pion= *(*itP);
		
		if(kaon1.charge()==kaon2.charge()) continue;
		double charge = kaon1.charge()+kaon2.charge()+pion.charge();
		HepLorentzVector p_dPlus=kaon1.p()+kaon2.p()+pion.p();
		double m=p_dPlus.mag();
		dMass=m;
		//	    cout <<"10 filling with " << m <<endl;
		dType=1;
		dDecay=3;
		if(foundDPiPi || foundSinglePionDecay)
		  foundDDoubleStarDecay=1;
		else
		  foundDDoubleStarDecay=0;
		if(SAVE_MESON_MASS_DISTRIBUTIONS)
		  {
		    if(m>1.7 && m < 2.0)
		      {
			if(m_mc)
			  {

			    const Mdst_charged &h_c1=kaon1.mdstCharged();
			    const Mdst_charged &h_c2=kaon2.mdstCharged();
			    const Mdst_charged &h_c3=pion.mdstCharged();

			    const Gen_hepevt &h_1=get_hepevt(h_c1);
			    const Gen_hepevt &h_2=get_hepevt(h_c2);
			    const Gen_hepevt &h_3=get_hepevt(h_c3);


			    dIsSignal=0;
			    if(h_1 && h_2 && h_3 && h_1.mother() && h_2.mother()  && h_3.mother() && h_1.mother().get_ID()==h_2.mother().get_ID() && h_1.mother().get_ID()==h_3.mother().get_ID()){
			      dIsSignal=1;
			    }
			  }
			mesonMassTree->Fill();
		      }	
		  }
		  	    
		  
		histoRecDSpect->Fill(m);
		if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));
		Particle* dPlus =new Particle(p_dPlus,Ptype(charge > 0 ? "D+": "D-"));
		dPlus->userInfo(*(new ParticleInfoDecay()));
		ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(dPlus->userInfo());
		pinf.decayChannel=dDecay;
		pinf.type=dType;
		pinf.pdgDiff=m-m_DPlus;
		pinf.mDRec=m;
		pinf.mRec=m;

		dPlus->relation().append(kaon1);
		dPlus->relation().append(kaon2);
		dPlus->relation().append(pion);
		if(m_mc)
		  {

		    const Mdst_charged &h_c1=kaon1.mdstCharged();
		    const Mdst_charged &h_c2=kaon2.mdstCharged();
		    const Mdst_charged &h_c3=pion.mdstCharged();

		    const Gen_hepevt &h_kaon1=get_hepevt(h_c1);
		    const Gen_hepevt &h_kaon2=get_hepevt(h_c2);
		    const Gen_hepevt &h_pion=get_hepevt(h_c3);


		    if(h_kaon1&&h_kaon2 && h_pion&& h_kaon1.mother() && h_kaon2.mother() && h_pion.mother() && h_kaon1.mother().get_ID()==h_kaon2.mother().get_ID() && h_kaon1.mother().get_ID()==h_pion.mother().get_ID()){
		      dPlus->relation().genHepevt(h_kaon1.mother());
		    }
		  }
		chargedDCandidates.push_back(dPlus);
	      }
	  }
      }
  }

  void bToDDoubleStar::reconstructDStar()
  {
    double m_DStarPlus=2.01027;
    double m_DStar0=2.00697;
    double m_D0=1.86484;
    double m_DPlus=1.86962;
    //Dmitry does not have max and min on DStar, just the difference (probably due to the uncertainty on the DStar)
    double m_dStarPlusmass_max=m_DStarPlus+10.015;
    double m_dStarPlusmass_min=m_DStarPlus-10.015;

    double m_dStar0mass_max=m_DStar0+10.015;
    double m_dStar0mass_min=m_DStar0-10.015;

    double max_massDifference=0.003;

    //    cout <<" combining " << chargedDCandidates.size() <<" charged Ds with " << pi0Candidates.size() <<" pi0s"<<endl;
    for(vector<Particle*>::iterator itD=chargedDCandidates.begin();itD!=chargedDCandidates.end();itD++)
      {
	//	break;
	for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
	  {
	    Particle& D= *(*itD);
	    Particle& pi0= *(*itPi0);

	    bool doubleUse=false;
	    //make sure that pi0 is not child of D
	    for(int i =0;i<D.nChildren();i++)
	      {
		if(D.child(i).relation().isIdenticalWith(pi0.relation()))
		  {
		    //		    cout <<"found double use 1 .." <<endl;
		    doubleUse=true;
		    break;
		  }
	      }
	    if(doubleUse)
	      continue;

	    HepLorentzVector p_dStar=D.p()+pi0.p();
	    double m=p_dStar.mag();
	    dMass=m;
	    //	    cout <<"11 filling with " << m <<endl;
	    dType=2;
	    dDecay=0;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m>1.9 && m < 2.1)
		  {
		    if(m_mc)
		      {
			const Gen_hepevt &h_D=D.genHepevt();
			const Mdst_pi0 &h_c1=pi0.mdstPi0();
			const Gen_hepevt &h_pi0=get_hepevt(h_c1);

			dIsSignal=0;
			//		if(h_D && h_pi0 && h_D.mother() && h_pi0.mother() && h_D.mother().get_ID()==h_pi0.mother().get_ID()){
			if(h_D && h_pi0 && commonParent(&h_D,&h_pi0,0,0,0))
			  {
			    dIsSignal=1;
			  }
		      } 
		    mesonMassTree->Fill();
		  }

	      }
	    histoRecDStarSpect->Fill(m);
	    histoRecDStarSpectToDPi0->Fill(m);
	    if(m>m_dStarPlusmass_max || m < m_dStarPlusmass_min ||isnan(m)) continue;
	    //	    cout <<"looking at dstar, mass diff: " <<(m_DStarPlus-m_DPlus) <<" vs : " << (m-D.p().mag());
	    //	    cout <<" gives: " << fabs(m-D.p().mag()-(m_DStarPlus-m_DPlus)) <<endl;

	    if(fabs(m-D.p().mag()-(m_DStarPlus-m_DPlus)) > max_massDifference) continue;
	    //	    cout <<"done " <<endl;
	    Particle* dStar =new Particle(p_dStar,Ptype(D.charge()>0 ? "D*+" : "D*-"));
	    dStar->userInfo(*(new ParticleInfoDecay()));
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(dStar->userInfo());
	    ParticleInfoDecay& d_pinf=dynamic_cast<ParticleInfoDecay&>(D.userInfo());
	    pinf.childDecayChannel=d_pinf.decayChannel;
	    pinf.childType=d_pinf.type;
	    pinf.childPdgDiff=d_pinf.pdgDiff;
	    pinf.decayChannel=dDecay;
	    pinf.type=dType;
	    pinf.pdgDiff=m-m_DStarPlus;
	    pinf.mRec=m;
	    pinf.mDRec=d_pinf.m;
	    dStar->relation().append(D);
	    dStar->relation().append(pi0);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Mdst_pi0 &h_c1=pi0.mdstPi0();
		const Gen_hepevt &h_pi0=get_hepevt(h_c1);

		if(h_D && h_pi0 && h_D.mother() && h_pi0.mother() && h_D.mother().get_ID()==h_pi0.mother().get_ID()){
		  dStar->relation().genHepevt(h_D.mother());
		}
	      } 
	    DStarCandidates.push_back(dStar);
	  }
      }


    //    cout <<" combining " << D0Candidates.size() <<"  D0s with " << chargedPiCandidates.size() <<" charged "<<endl;
    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
      {
	for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	  {
	    Particle& D= *(*itD);
	    Particle& pion= *(*itP);

	    bool doubleUse=false;
	    //make sure that pion is not child of D-->what about Ks?
	    if(checkDoubleUse(D,pion))
	      continue;

	    HepLorentzVector p_dStar=D.p()+pion.p();
	    double m=p_dStar.mag();
	    dMass=m;
	    dType=2;
	    dDecay=1;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m>1.9 && m < 2.1)
		  {
		    if(m_mc)
		      {
			const Gen_hepevt &h_D=D.genHepevt();

			const Mdst_charged &h_c1=pion.mdstCharged();
			const Gen_hepevt &h_pion=get_hepevt(h_c1);

			dIsSignal=0;
			//			if(h_D && h_pion && h_D.mother() && h_pion.mother() && h_D.mother().get_ID()==h_pion.mother().get_ID()){
			if(h_D && h_pion && commonParent(&h_D, &h_pion,0,0,0))
			  {
			    dIsSignal=1;
			  }
		      } 
		    mesonMassTree->Fill();
		  }
	      }
	      
	    histoRecDStarSpect->Fill(m);
	    histoRecDStarSpectToD0Pi->Fill(m);
	    if(m>m_dStarPlusmass_max || m < m_dStarPlusmass_min ||isnan(m)) continue;
	    //	    cout <<"D0+pi"<<endl;
	    //	        cout <<"m -D: "<< m-D.p().mag() <<endl;
	    //	    	    cout <<"looking at dstar, mass diff: " <<(m_DStarPlus-m_D0) <<" vs : " << (m-D.p().mag());
	    //	        cout <<" gives: " << fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) <<endl;
	    if(fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) > max_massDifference) continue;
	    //	    cout <<"done" <<endl;
	    Particle* dStar =new Particle(p_dStar,Ptype(pion.charge()>0 ? "D*+" : "D*-"));
	    dStar->userInfo(*(new ParticleInfoDecay()));
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(dStar->userInfo());
	    ParticleInfoDecay& d_pinf=dynamic_cast<ParticleInfoDecay&>(D.userInfo());
	    pinf.childDecayChannel=d_pinf.decayChannel;
	    pinf.childType=d_pinf.type;
	    pinf.childPdgDiff=d_pinf.pdgDiff;
	    pinf.decayChannel=dDecay;
	    pinf.type=dType;
	    pinf.pdgDiff=m-m_DStarPlus;
	    pinf.mRec=m;
	    pinf.mDRec=d_pinf.m;
	    dStar->relation().append(D);
	    dStar->relation().append(pion);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Mdst_charged &h_c1=pion.mdstCharged();
		const Gen_hepevt &h_pion=get_hepevt(h_c1);

		if(h_D && h_pion && h_D.mother() && h_pion.mother() && h_D.mother().get_ID()==h_pion.mother().get_ID()){
		  dStar->relation().genHepevt(h_D.mother());
		}
	      } 
	    DStarCandidates.push_back(dStar);
	  }
      }

    ///D*0 from D0pi0
    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
      {
	//	break;
	for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
	  {
	    Particle& D= *(*itD);
	    Particle& pi0= *(*itPi0);

	    bool doubleUse=false;
	    //make sure that pi0 is not child of D
	    for(int i =0;i<D.nChildren();i++)
	      {
		if(D.child(i).relation().isIdenticalWith(pi0.relation()))
		  {
		    //		    cout <<"found double use 3 .." <<endl;
		    doubleUse=true;
		    break;
		  }
	      }
	    if(doubleUse)
	      continue;


	    HepLorentzVector p_dStar=D.p()+pi0.p();
	    double m=p_dStar.mag();
	    dMass=m;
	    //	    cout <<"12 filling with " << m <<endl;
	    dType=2;
	    dDecay=2;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m>1.9 && m < 2.1)
		  {
		    if(m_mc)
		      {
			const Gen_hepevt &h_D=D.genHepevt();
			const Mdst_pi0 &h_c1=pi0.mdstPi0();
			const Gen_hepevt &h_pi0=get_hepevt(h_c1);

			dIsSignal=0;
			//			if(h_D && h_pion && h_D.mother() && h_pion.mother() && h_D.mother().get_ID()==h_pion.mother().get_ID()){
			if(h_D && h_pi0 && commonParent(&h_D,&h_pi0,0,0,0))
			  {
			    dIsSignal=1;
			  }
		      } 
		    mesonMassTree->Fill();
		  }


	      }
	    histoRecDStarSpect->Fill(m);
	    histoRecDStarSpectToD0Pi0->Fill(m);
	    if(m>m_dStar0mass_max || m < m_dStar0mass_min ||isnan(m)) continue;
	    //	    cout <<"looking at dstar, mass diff: " <<(m_DStar0-m_D0) <<" vs : " << (m-D.p().mag());
	    //	    cout <<" gives: " << fabs(m-D.p().mag()-(m_DStar0-m_D0)) <<endl;
	    if(fabs(m-D.p().mag()-(m_DStar0-m_D0)) > max_massDifference) continue;
	    //	    cout <<" done " <<endl;
	    Particle* dStar =new Particle(p_dStar,Ptype("D*0"));
	    dStar->userInfo(*(new ParticleInfoDecay()));
	    ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>(dStar->userInfo());
	    ParticleInfoDecay& d_pinf=dynamic_cast<ParticleInfoDecay&>(D.userInfo());
	    pinf.childDecayChannel=d_pinf.decayChannel;
	    pinf.childType=d_pinf.type;
	    pinf.childPdgDiff=d_pinf.pdgDiff;
	    pinf.pdgDiff=m-m_DStar0;
	    pinf.mRec=m;
	    pinf.mDRec=d_pinf.m;
	    pinf.decayChannel=dDecay;
	    pinf.type=dType;

	    dStar->relation().append(D);
	    dStar->relation().append(pi0);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();

		const Mdst_pi0 &h_c1=pi0.mdstPi0();
		const Gen_hepevt &h_pi0=get_hepevt(h_c1);

		if(h_D && h_pi0 && h_D.mother() && h_pi0.mother() && h_D.mother().get_ID()==h_pi0.mother().get_ID()){
		  dStar->relation().genHepevt(h_D.mother());
		}
	      } 
	    DStarCandidates.push_back(dStar);
	  }
      }

  }


  //this seems to have some impact...
  unsigned bToDDoubleStar::doKmFit(Particle &p, double& confLevel, int debug, double mass)
  {
    //    return true;
    kmassfitter km;
    if(mass!=0)
      {
	//	cout <<" using mass; "<< mass <<endl;
	km.invariantMass(mass);
      }
    //    km.invariantMass(mass==0 ? p.pType().mass(): mass);
    else
      km.invariantMass(mass==0 ? p.pType().mass(): mass);
    for(unsigned j=0;j<p.relation().nChildren();j++)
      {
	Particle child=p.relation().child(j);
	km.addTrack(child.momentum().p(),child.momentum().x(),child.momentum().dpx(),child.pType().charge(),child.pType().mass());
      }
    km.notDecayPoint();
    unsigned err = km.fit();
    if(err){
      //           cout <<"Err in kmassvertexfitter: "<< err <<endl;
      return 0;
    }
    //    else{cout <<"fit was ok.." <<endl;}
    confLevel=km.cl();
    return makeMother(km,p);
  }

  //dmitries code..
  unsigned bToDDoubleStar::doKmVtxFit2(Particle &p, double& confLevel, int debug, double mass)
  {
    kmassvertexfitter kmvfitter;
    kmvfitter.invariantMass(mass==0 ? p.pType().mass() : mass);
    for(unsigned i=0; i<p.nChildren(); ++i)
      addTrack2fit(kmvfitter, p.child(i));

    if(p.nChildren()>2)
      {
	
      }
    //this is in the example, but probably old interface?
    //    kmvfitter.vertex(IpProfile::position());
    //    kmvfitter.errVertex(IpProfile::position_err());
    kmvfitter.initialVertex(IpProfile::position());
    //no error
    if(!kmvfitter.fit()) {
      makeMother(kmvfitter, p);
      p.momentum().vertex(kmvfitter.vertex(),kmvfitter.errVertex());
      confLevel=kmvfitter.cl();
      return true;
    }
    return false;
  }
  //this seems to have some impact...
  unsigned bToDDoubleStar::doKmVtxFit(Particle &p, double& confLevel, int debug)
  {

    //first get vertex:
    kvertexfitter vtxFit;
    vtxFit.initialVertex(p.x());
    for(unsigned j=0;j<p.relation().nChildren();j++)
      {
	//	Particle child=p.relation().child(j);
	//	vtxFit.addTrackToFit(child.momentum().p(),child.momentum().x(),child.momentum().dpx(),child.pType().charge(),child.pType().mass());
      }

    //    return true;
    kmassvertexfitter km;
    km.initialVertex(p.x());
    km.invariantMass(p.pType().mass());
    for(unsigned j=0;j<p.relation().nChildren();j++)
      {
	Particle child=p.relation().child(j);
	km.addTrack(child.momentum().p(),child.momentum().x(),child.momentum().dpx(),child.pType().charge(),child.pType().mass());
      }
    //    km.notDecayPoint();
    unsigned err = km.fit();
    if(err){
      //           cout <<"Err in kmassvertexfitter: "<< err <<endl;
      return 0;
    }
    //    else{cout <<"fit was ok.." <<endl;}
    confLevel=km.cl();
    return makeMother(km,p);
  }
  int bToDDoubleStar::getOverlap(set<int>& s1, set<int>& s2)
  {
    int ret=0;
    for(set<int>::iterator it=s1.begin();it!=s1.end();it++)
      {
	for(set<int>::iterator it2=s2.begin();it2!=s2.end();it2++)
	  {
	    if((*it)==(*it2))
	      {
		ret++;
	      }
	  }
      }
    return ret;
  }


  //returns correction factor and (absolute)uncertainty
  pair<float,float> bToDDoubleStar::getBRCorrection(int& decType)
  {
    float mc=1.0;
    float data=1.0;
    float error=0.0;

    float corrFactor=1.0;
    float relError=0.0;

    if(br_sig_D0LNu)
      {
	decType=0;
	mc=mcFactors[br_sig_D0];
	data= dataFactors[br_sig_D0];
	error=dataFactorError[br_sig_D0];
	corrFactor*=(data/mc);
	//the division by mc cancels
	relError+=error*error/(data*data);

      }
    if(br_sig_DLNu)
      {
	decType=1;
	//	cout <<"correction... d lnu" <<endl;
	mc= mcFactors[br_sig_D];
	data= dataFactors[br_sig_D];
	error=dataFactorError[br_sig_D];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);

      }
    if(br_sig_DStarLNu)
      {
	decType=2;
	//	cout <<"correction... dstar lnu" <<endl;
	mc= mcFactors[br_sig_DStar];
	data= dataFactors[br_sig_DStar];
	error=dataFactorError[br_sig_DStar];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }

    if(br_sig_DStar0LNu)
      {
	decType=3;
	mc= mcFactors[br_sig_DStar0];
	data= dataFactors[br_sig_DStar0];
	error=dataFactorError[br_sig_DStar0];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(br_sig_D1LNu)
      {
	decType=4;
	mc= mcFactors[br_sig_D1];
	data= dataFactors[br_sig_D1];
	error=dataFactorError[br_sig_D1];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(br_sig_D2LNu)
      {
	decType=5;
	mc= mcFactors[br_sig_D2];
	data= dataFactors[br_sig_D2];
	error=dataFactorError[br_sig_D2];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }

    if(br_sig_D1PrimeLNu)
      {
	decType=6;
	mc= mcFactors[br_sig_D1Prime];
	data= dataFactors[br_sig_D1Prime];
	error=dataFactorError[br_sig_D1Prime];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(br_sig_D0StarLNu)
      {
	decType=7;
	mc= mcFactors[br_sig_D0Star];
	data= dataFactors[br_sig_D0Star];
	error=dataFactorError[br_sig_D0Star];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(br_sig_D10LNu)
      {
	decType=8;
	mc= mcFactors[br_sig_D10];
	data= dataFactors[br_sig_D10];
	error=dataFactorError[br_sig_D10];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(br_sig_D20LNu)
      {
	decType=9;
	mc= mcFactors[br_sig_D20];
	data= dataFactors[br_sig_D20];
	error=dataFactorError[br_sig_D20];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(br_sig_D1Prime0LNu)
      {
	decType=10;
	mc= mcFactors[br_sig_D1Prime0];
	data= dataFactors[br_sig_D1Prime0];
	error=dataFactorError[br_sig_D1Prime0];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(br_sig_D0Star0LNu)
      {
	decType=11;
	mc= mcFactors[br_sig_D0Star0];
	data= dataFactors[br_sig_D0Star0];
	error=dataFactorError[br_sig_D0Star0];
	relError+=error*error/(data*data);
	corrFactor*=(data/mc);
      }
    if(!mc || !data)
      {
	cout <<"B Br no mc or data: "<< mc << data <<endl;
	for(int i=0;i<12;i++)
	  {
	    cout <<" mc corr " << i <<": " << mcFactors[i] <<" data; "<< dataFactors[i] <<endl;
	  }
	exit(1);
      }
    //    return data/mc;

    //this enables to have a multiplicative factor
   return  pair<float,float>(corrFactor,sqrt(relError)*corrFactor);

  }

  float bToDDoubleStar::calcBRCorrection()
  {
    //differentiate in MC and data
    //need to correct for Dlnu and the hadronic D decay
    cout <<"calculating br corrections" << endl;
    //DStarlNu
    mcFactors[br_sig_DStar]=0.0533;
    dataFactors[br_sig_DStar]=0.0493;
    dataFactorError[br_sig_DStar]=0.0011;
    cout <<" dstar: "<< br_sig_DStar <<endl;
    for(int i=0;i<12;i++)
      {
	cout <<"first mc: " << mcFactors[i] << " data: "<< dataFactors[i] <<endl;
      }

    //DStar0lNu
    cout <<" dstar0: "<< br_sig_DStar0 <<endl;
    mcFactors[br_sig_DStar0]=0.0579;
    dataFactors[br_sig_DStar0]=0.0569;
    dataFactorError[br_sig_DStar0]=0.0019;

    //DlNu
    mcFactors[br_sig_D]=0.0219;
    dataFactors[br_sig_D]=0.0219;
    dataFactorError[br_sig_D]=0.0012;

    //D0lNu
    mcFactors[br_sig_D0]=0.0213;
    dataFactors[br_sig_D0]=0.0227;
    dataFactorError[br_sig_D0]=0.001;
    //D1 lNu

    mcFactors[br_sig_D1]=0.0074;
    dataFactors[br_sig_D1]=0.0074;
    dataFactorError[br_sig_D1]=0.0011;
    //D2 lNu
    mcFactors[br_sig_D2]=0.0036;
    dataFactors[br_sig_D2]=0.0047;
    dataFactorError[br_sig_D2]=0.0017;
    //D1Prime lNu
    mcFactors[br_sig_D1Prime]=0.002;
    dataFactors[br_sig_D1Prime]=0.0026;
    dataFactorError[br_sig_D1Prime]=0.0009;
    //D0Star lNu
    mcFactors[br_sig_D0Star]=0.0084;
    dataFactors[br_sig_D0Star]=0.0052;
    dataFactorError[br_sig_D0Star]=0.0022;
    //D10 lNu
    mcFactors[br_sig_D10]=0.0074;
    dataFactors[br_sig_D10]=0.0074;
    dataFactorError[br_sig_D10]=0.0011;
    //D20 lNu
    mcFactors[br_sig_D20]=0.0036;
    dataFactors[br_sig_D20]=0.0047;
    dataFactorError[br_sig_D20]=0.0017;
    //D1Prime0 lnu
    mcFactors[br_sig_D1Prime0]=0.002;
    dataFactors[br_sig_D1Prime0]=0.0026;
    dataFactorError[br_sig_D1Prime0]=0.0009;
    //D0 Star0 lnu
    mcFactors[br_sig_D0Star0]=0.0084;
    dataFactors[br_sig_D0Star0]=0.0052;
    dataFactorError[br_sig_D0Star0]=0.0022;

    /////calculate the D factors    
    cout <<"calced everything, : ";
    for(int i=0;i<12;i++)
      {
	cout <<"mc: " << mcFactors[i] << " data: "<< dataFactors[i] <<endl;
      }
    return 1.0;
  }


  void bToDDoubleStar::addBremsPhoton(Particle* p)
  {
    float minTheta=100;
    float minPx=-1;
    float minPy=-1;
    float minPz=-1;     
    float minE=-1;

    Mdst_gamma_Manager& gamma_mgr=Mdst_gamma_Manager::get_manager();
    for(std::vector<Mdst_gamma>::const_iterator i =gamma_mgr.begin();i!=gamma_mgr.end();i++)
      {

	if(passGammaQACuts(i))
	  {
	    const Mdst_gamma& gam=*i;
	    double px=gam.px();
	    double py=gam.py();
	    double pz=gam.pz();
	    double gammaE=sqrt(px*px+py*py+pz*pz);
	    Hep3Vector pGamma(px,py,pz);
	    float theta=  p->p3().angle(pGamma);
	    float deg=180*(theta/TMath::Pi());
	    if(deg < minTheta)
	      {
		minTheta=deg;
		minPx=px;
		minPy=py;
		minPz=pz;
		minE=gammaE;
	      }
	  }

	
      }
    if(minTheta<5)
      {
	Momentum mom=p->momentum();
	HepLorentzVector lvGamma(minPx,minPy,minPz,minE);
	mom.momentum(p->p()+lvGamma);
      }


  }
  bool bToDDoubleStar::passGammaQACuts(std::vector<Mdst_gamma>::const_iterator i)
  {

    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();
    if(gammaIds.find(i->get_ID())!=gammaIds.end()){
      return false;
    }
    Hep3Vector h(i->px(),i->py(),i->pz());
    const Mdst_gamma& gam=*i;
    int id=(int)gam.get_ID();
    double px=gam.px();
    double py=gam.py();
    double pz=gam.pz();
    //does not make sensee because if we only look in the central region for the
    //computation of the thrust axis, that would change the axis, same for charged
    ///there I take all for the thrust axis and the ones in the central region for
    //the asymmetry
    Mdst_ecl_aux &aux =eclaux_mgr(Panther_ID(gam.ecl().get_ID()));
    if(gam.ecl().quality()!=0)
      {
	//	    cout <<"loosing photon due to  quality " <<endl;
	return false;
      }
    //ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
    double e9oe25 =aux.e9oe25();
    double gammaE=sqrt(px*px+py*py+pz*pz);
    Hep3Vector photVec(px, py,pz);
    //barrel energy cut is the lowest
    if(gammaE<cuts::minGammaEBarrel)
      {
	//	    cout <<" loosing photon in barrel due to energy cut " << gammaE<<endl;
	return false;
      }
    float photTheta= photVec.theta();
    //	cout <<"photTheta: "<< photTheta<<endl;
    photTheta=180*(photTheta/TMath::Pi());
    if(photTheta >cuts::barrelThetaMax) 
      {
	if(gammaE<cuts::minGammaEEndcapBkwd)
	  {
	    //	    cout <<" loosing photon in forwrad endcap due to energy cut " << gammaE<<endl;
	    return false;
	  }
      }
    if(photTheta< cuts::barrelThetaMin)
      {
	if(gammaE<cuts::minGammaEEndcapFwd)
	  {
	    //	    cout <<" loosing photon in backward endcap due to energy cut " << gammaE<<endl;
	    return false;
	  }
      }

    return true;
  }


  bool bToDDoubleStar::commonParent(const Gen_hepevt* h1, const Gen_hepevt* h2, const Gen_hepevt* h3, const Gen_hepevt* h4, int lund)
  {
    int id1,id2,id3,id4;
    //      cout <<"getting first parent.." <<endl;
    if(!h1->mother() || !h2->mother())
      return false;
    //only do this for non charged pion kaon, for those we only go up one level..
    //    if(fabs(h1->idhep())==PY_PI || fabs(h1->idhep())==PY_K)
    {
      //	id1=h1->mother().get_ID();
    }
    //    else
    {
      id1=getParentDId(h1->mother());
    }
    //      cout <<"getting second parent " << endl;
    //    if(fabs(h2->idhep())==PY_PI || fabs(h2->idhep())==PY_K)
    {
      //	id2=h2->mother().get_ID();
    }
    //    else
    {
      id2=getParentDId(h2->mother());
    }
    if(h3 && !h4)
      {
	//  cout <<"checking for h3 mother.. " << endl;
	if(!h3->mother())
	  return false;
	//	  cout <<"trying third.." <<endl;
	//	if(fabs(h3->idhep())==PY_PI || fabs(h3->idhep())==PY_K)
	{
	  //  id3=h3->mother().get_ID();
	}
	//	  else
	{
	  id3=getParentDId(h3->mother());
	}
	if(id1==id2 && id2==id3)
	  return true;
	else
	  return false;
	//	  cout <<"done .. " << endl;
      }
    if(h4)
      {
	if(!h3->mother())
	  return false;

	if(!h4->mother())
	  return false;

	//	  cout <<"getting fourth " << endl;
	//	if(fabs(h3->idhep())==PY_PI || fabs(h3->idhep())==PY_K)
	{
	  //	    id3=h3->mother().get_ID();
	}
	//	else
	{
	  id3=getParentDId(h3->mother());
	}
	//	if(fabs(h4->idhep())==PY_PI || fabs(h4->idhep())==PY_K)
	{
	  //	    id4=h4->mother().get_ID();
	}
	//	else
	{
	  id4=getParentDId(h4->mother());
	}
	if(id1==id2 && id2==id3 && id2==id4)
	  return true;
	else
	  return false;
      }
    if(id1==id2)
      return true;
    else
      return false;
  }


  int bToDDoubleStar::getParentDId(const Gen_hepevt& h)
  {
    //    cout <<" in parent id " << endl;
    int id=fabs(h.idhep());
    //    cout <<" got id " << id <<endl;
    if(id==PY_D|| id==PY_D0 || id==PY_DStar || id==PY_DStar0)
      {
	//	cout <<"returning id " << endl;
	//	cout <<"ID" <<h.get_ID()<<endl;
	return h.get_ID();
      }
    //    cout <<" trying mother .. " << endl;
    if(h.mother())
      {
	//	cout <<"found " << endl;
	return getParentDId(h.mother());
      }
    return -1;
  }


  bool bToDDoubleStar::foundUnwantedDDoubleStar()
  {
    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();

    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {
	if(gen_it->get_ID()==PY_DStar_2S)
	  return true;
	if(gen_it->get_ID()==100421)
	  return true;
	if(gen_it->get_ID()==100413)
	  return true;
	if(gen_it->get_ID()==100411)
	  return true;
	if(gen_it->get_ID()==10431)
	  return true;
	
      }
    return false;
  }

  bool  bToDDoubleStar::checkDoubleUse(Particle& D, Particle& pion)
  {

    for(int i =0;i<D.nChildren();i++)
      {
	int lund=fabs(D.child(i).pType().lund());
	if(lund==PY_PI || lund==PY_Pi0 || lund==PY_K||lund==PY_GAMMA)
	  {
	    if(D.child(i).relation().isIdenticalWith(pion.relation()))
	      return true;
	  }
	else//something with daughters, so either D or KS
	  {
	    if(checkDoubleUse(D.child(i),pion))
	      return true;
	  }

      }
    return false;
  

    ////only for one level
    ///    for(int i =0;i<D.nChildren();i++)
    ///      {
    ///	//gotta check the children of the ks as well
    ///	if(D.child(i).pType().lund()==PY_KS0)
    ///	  {
    ///	    for(int j =0;j<D.child(i).nChildren();j++)
    ///	      {
    ///		if(D.child(i).child(j).relation().isIdenticalWith(pion.relation()))
    ///		  {
    ///		    return true;
    ///		  }
    ///	      }
    ///	    
    ///	  }
    ///	if(D.child(i).relation().isIdenticalWith(pion.relation()))
    ///	  {
    ///	    //		    cout <<"found double use 2 .." <<endl;
    ///	    return true;
    ///
    ///	  }
    ///      }
    ///    return false;
  }


  void  bToDDoubleStar::printUse()
  {

    cout <<"we have " << chargedPiCandidates.size() << "pi cand " << chargedKCandidates.size() <<" k cand ";
    cout <<pi0Candidates.size() <<" pi0 cand " << KsCandidates.size() <<" Ks Can" <<endl;

    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
      {
	ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>((*itD)->userInfo());
	cout <<"D0 decay : " << pinf.decayChannel << " mass diff: "<< pinf.pdgDiff <<endl;
	printChildrenIds(*(*itD));
      }
    for(vector<Particle*>::iterator itD=chargedDCandidates.begin();itD!=chargedDCandidates.end();itD++)
      {
	ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>((*itD)->userInfo());
	cout <<"D+ decay : " << pinf.decayChannel << " mass diff: "<< pinf.pdgDiff <<endl;
	printChildrenIds(*(*itD));
      }
    for(vector<Particle*>::iterator itD=DStarCandidates.begin();itD!=DStarCandidates.end();itD++)
      {
	ParticleInfoDecay& pinf=dynamic_cast<ParticleInfoDecay&>((*itD)->userInfo());
	cout <<"DStar decay : " << pinf.decayChannel << " child: " << pinf.childType<<" child decay: " << pinf.childDecayChannel<<" child pdg diff: " << pinf.childPdgDiff<<endl;
	printChildrenIds(*(*itD));
      }


  }

  void bToDDoubleStar::printChildrenIds(Particle& p)
  {

    for(int i =0;i<p.nChildren();i++)
      {
	int lund=fabs(p.child(i).pType().lund());
	if(lund==PY_PI || lund==PY_Pi0 || lund==PY_K||lund==PY_GAMMA)
	  {
	    cout <<" found " <<p.child(i).pType().name() <<" id: ";
	    if(p.child(i).mdstCharged())
	      {
		cout <<p.child(i).mdstCharged().get_ID() <<" momentum: " << p.child(i).p().rho();

	      }
	    if(p.child(i).mdstPi0())
	      {
		cout <<p.child(i).mdstPi0().get_ID() <<" momentum: " << p.child(i).p().rho();
	      }
	    if(p.child(i).mdstGamma())
	      {
		cout << p.mdstGamma().get_ID() <<" momentum: " << p.child(i).p().rho();
	      }
	    cout <<endl;
	  }
	else//something with daughters, so either D or KS
	  {
	    printChildrenIds(p.child(i));
	  }

      }
    return;
  

  }

  ///same as the gammaQA function
  bool bToDDoubleStar::checkGammaEnergy(float px, float py, float pz, float gammaE)
  {

    Hep3Vector photVec(px, py,pz);
    //barrel energy cut is the lowest
    if(gammaE<cuts::minGammaEBarrel)
      {
	//	    cout <<" loosing photon in barrel due to energy cut " << gammaE<<endl;
	return false;
      }
    float photTheta= photVec.theta();
    //	cout <<"photTheta: "<< photTheta<<endl;
    photTheta=180*(photTheta/TMath::Pi());
    if(photTheta >cuts::barrelThetaMax) 
      {
	if(gammaE<cuts::minGammaEEndcapBkwd)
	  {
	    //	    cout <<" loosing photon in forwrad endcap due to energy cut " << gammaE<<endl;
	    return false;
	  }
      }
    if(photTheta< cuts::barrelThetaMin)
      {
	if(gammaE<cuts::minGammaEEndcapFwd)
	  {
	    //	    cout <<" loosing photon in backward endcap due to energy cut " << gammaE<<endl;
	    return false;
	  }
      }
    return true;

  }

  bool compPi0s(Particle* p1, Particle* p2)
  {
    if(!(p1->mdstPi0() && p2->mdstPi0()))
      {
	cout  <<endl <<" trying to compare non pi0!!!!!!" <<endl<<endl;
	return false;
      }

    ParticleInfoMass& pinf1=dynamic_cast<ParticleInfoMass&>(p1->userInfo());
    ParticleInfoMass& pinf2=dynamic_cast<ParticleInfoMass&>(p2->userInfo());
    float asym1=(pinf1.gammaE1-pinf1.gammaE2)/(pinf1.gammaE1+pinf1.gammaE2);
    float asym2=(pinf2.gammaE1-pinf2.gammaE2)/(pinf2.gammaE1+pinf2.gammaE2);

float gammaE1_1=pinf1.gammaE1;
float gammaE2_1=pinf1.gammaE2;
float gammaE1_2=pinf2.gammaE1;
float gammaE2_2=pinf2.gammaE2;

//    return (asym1<asym2);

float tmp;
//_1 should be the larger one
if(gammaE1_1< gammaE2_1)
  {
    tmp=gammaE2_1;
    gammaE2_1=gammaE1_1;
    gammaE1_1=tmp;
  }
if(gammaE1_2< gammaE2_2)
  {
    tmp=gammaE2_2;
    gammaE2_2=gammaE1_2;
    gammaE1_2=tmp;
  }


//use Robins criterium. We flipped the '>' sign because the sorting is done ascending, but we want descending in energy
if(gammaE1_1!=gammaE1_2)
  return gammaE1_1 > gammaE1_2;
 else
   return gammaE2_1> gammaE2_2;
}


void bToDDoubleStar::cleanPi0s()
{
  //sort pi0s that passed our cuts by energy asymmetry, then go from firt to last and remove stuff that overlaps
  //	if(abs((g1Energy-g2Energy)/(g1Energy+g2Energy))
  sort(pi0Candidates.begin(),pi0Candidates.end(),compPi0s);
  bool deletedOne=false;
  do
  {
    //new round
    deletedOne=false;
    for(vector<Particle*>::iterator itp=pi0Candidates.begin();itp!=pi0Candidates.end();itp++)
      {

	if((*itp)->nChildren()==2 && (*itp)->child(0).mdstGamma() && (*itp)->child(1).mdstGamma())
	  {
	   int id1=(*itp)->child(0).mdstGamma().get_ID();
	   int id2=(*itp)->child(1).mdstGamma().get_ID();
	    for(vector<Particle*>::iterator itp2=itp+1;itp2!=pi0Candidates.end();itp2++)
	      {
		if((*itp2)->nChildren()==2 && (*itp2)->child(0).mdstGamma() && (*itp2)->child(1).mdstGamma())
		  {
		    int id1_2=(*itp2)->child(0).mdstGamma().get_ID();
		    int id2_2=(*itp2)->child(1).mdstGamma().get_ID();
		    
		    if(id1==id1_2 || id1==id2_2 || id2== id1_2 || id2==id2_2)
		      {
			deletedOne=true;
			//	cout <<" found overlap... erase.." <<endl;
			delete *itp2;
			pi0Candidates.erase(itp2);
		    //gotta get out of the loop if we deleted something
			break;
		      }
		    
		  }
		else
		  {
		    cout <<"some problem with getting the gammas of the second pi0s.." <<endl;
		  }

	      }
	  }
	else
	  {
	    cout <<"some problem with getting the gammas of the pi0s.." <<endl;
	  }
	
	if(deletedOne)
	  break;


      }

  }
  while(deletedOne);



}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
