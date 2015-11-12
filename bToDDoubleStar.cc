const bool PRINT=true;
const bool SAVE_MESON_MASS_DISTRIBUTIONS=false;
#include <iomanip>
#include "mdst/findKs.h"
#include "bToDDoubleStar/mc.h"  //one central place to put the define mc
#include "tables/ekpfullrecon_panther.h"
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

#include "fastjet/ClusterSequence.hh"
#include <iostream>

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
    mesonMassTree->Branch("foundDDoubleStarDecay",&foundDDoubleStarDecay,"foundDDoubleStarDecay/I");
    mesonMassTree->Branch("allDTracksFound",&allDTracksFound,"allDTracksFound/I");
    ///end 


    m_mc=false;
#ifdef MC
    m_mc=true;
#endif
    //    cout <<"constructor..." <<endl;
    strcpy(rFileName,"notInitialized.root");
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
    //    cout <<"init.." <<endl;
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
    //    cout <<" begin run  " << endl;
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

  // event function
  void bToDDoubleStar::event(BelleEvent* evptr, int* status)
  {
    bool bgFlag=false;
    bool foundDDStarFlag=false;
    //for the output of the decay signature query
    sig_FoundDDoubleStar=false;
    sig_numLeptons=0;
    sig_numKaons=0;
    sig_numPions=0;
    sig_numD=0;
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


    found_2SD=false;
    found_2SD_Star=false;

    const double m_pi0=0.1349766;
    vector<int> foundDecIds;
    int numNu=0;
    vector<int> decIds;
    vector<int> pids;
    int bMesonId;
    foundSinglePionDecay=false;
    noDRec=false;
    if(m_mc)
      {
	//find b meson id in the simulation
   
	foundDPiPi=checkForDPiPi(bMesonId,foundSinglePionDecay);
	if(foundDPiPi || foundSinglePionDecay)
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
    const double _log_nbout_min=-3;
    //    const double _tag_Mbc=5.25;

    //    cout <<" valid run " << validRun<<endl;
    if(!validRun)
      return;
    int evtNr;
    int runNr;
 
    /////for xcheck
    AnaBrecon brecon;
    evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    runNr=Belle_event_Manager::get_manager().begin()->RunNo();
    //    cout <<" evt: " << evtNr << " runNr " << runNr <<endl;
    kinematics::runNr=runNr;
    kinematics::evtNr=evtNr;

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
      return;
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

	//	cout <<"looking at " << bcand.nChildren() <<" children "<<endl;
	for(int i=0;i<bcand.nChildren();i++)
	  {
	    const Particle& p=bcand.child(i);
	    //	    	    cout <<"child " << i << ": " << p.pType().name() <<" pcode; " << p.pType().lund()<<endl;
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
			//			cout <<"--->--->"<<p3.pType().name()<<endl;
		      }
		  }
	      }
	  }
	//	cout <<"about to get decay ids .." <<endl;

	//	cout <<"got " << chargedIds.size() << " charged Ids, " << pi0Ids.size() <<" pi0s, " << gammaIds.size() <<" gammas " << endl;
      }
    
    //need good b
    if(bestLogProb<_log_nbout_min)
      {
	exitEvent();
	return;
      }
    
    //    cout <<"best log prob: "<< bestLogProb<<endl;
    Particle & bestBcand = const_cast<Particle &>(brecon.getParticle( (int)bestEKP_B.get_ID() ));
    bestBPx=bestBcand.px();
    bestBPy=bestBcand.py();
    bestBPz=bestBcand.pz();
    //    cout <<"looking at b cand" <<endl;

    //    if(bestBcand.mdstVee2())
    //      cout <<"b has vee2 " <<endl;
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
    treeData.found_2SD=found_2SD;
    treeData.found_2SD_Star=found_2SD_Star;
    treeData.foundAnyDDoubleStar=sig_FoundDDoubleStar;
    treeData.sig_numPions=sig_numPions;
    treeData.sig_numD=sig_numD;
    treeData.sig_numDStar=sig_numDStar;
    treeData.sig_numKaons=sig_numKaons;
    treeData.sig_numPi0=sig_numPi0;
    treeData.sig_numLeptons=sig_numLeptons;
    treeData.sig_numBaryons=sig_numBaryons;
    treeData.tagOverlapFractionCharged=overlapFractionCharged;
    treeData.tagOverlapFractionPi0=overlapFractionPi0;
    if((sigDStarLNu || sigDStarPiLNu || sigDStarPiPiLNu )&& sig_FoundDDoubleStar)
      {
	//	cout <<"found dlnu and ddouble star!!!!!" <<endl;

      }

    treeData.sigDLNu=sigDLNu;
    treeData.sigDPiLNu=sigDPiLNu;
    treeData.sigDPiPiLNu=sigDPiPiLNu;

    treeData.sigDStarLNu=sigDStarLNu;
    treeData.sigDStarPiLNu=sigDStarPiLNu;
    treeData.sigDStarPiPiLNu=sigDStarPiPiLNu;

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
	  }
	else
	  {
	    delete p;
	  }

      }

  

    //    cout <<"there are " << mdst_chr_Mgr.size() << " charge tracks in mdst_chr " <<endl;
    for(Mdst_charged_Manager::iterator chr_it=mdst_chr_Mgr.begin();chr_it!=mdst_chr_Mgr.end();chr_it++)
      {
	//part of the tag
	if(chargedIds.find(chr_it->get_ID())!=chargedIds.end()){
	  //	 	  cout <<"charged track is part of track " <<endl;
	  continue;
	}
	//      cout <<"looking at " <<(*chr_it).p(0) <<" " << (*chr_it).p(1) <<" " << (*chr_it).p(2) <<endl;

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
	      strcpy(ptypName,"E+");
	    else
	      strcpy(ptypName,"E-");
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
	      strcpy(ptypName,"MU+");
	    else
	      strcpy(ptypName,"MU-");
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
	cout <<"e_id : "<< e_id <<" mu_id: "<< mu_id <<endl;
	//default pion, good enough
	//	if(!positivelyIdentified)
	//	  continue;

	double dr, dz, refitPx, refitPy, refitPz;
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


	Particle* p=new Particle(*chr_it,string(ptypName));
	  
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
	if(mass>0.15 || mass<0.12)
	  continue;
	float pLab=sqrt(px*px+py*py+pz*pz);
	//      cout <<"pi0mass: "<< mass <<endl;≈

	float g1Energy= sqrt(pi0.gamma(0).px()*pi0.gamma(0).px()+pi0.gamma(0).py()*pi0.gamma(0).py()+pi0.gamma(0).pz()*pi0.gamma(0).pz());
	float g2Energy= sqrt(pi0.gamma(1).px()*pi0.gamma(1).px()+pi0.gamma(1).py()*pi0.gamma(1).py()+pi0.gamma(1).pz()*pi0.gamma(1).pz());
	//	cout <<"pi0 gamma1: "<< g1Energy <<" gamma2: "<< g2Energy <<endl;

	//	if(abs((g1Energy-g2Energy)/(g1Energy+g2Energy))>cuts::maxPi0GAsym)
	//	  continue;
	//let's make this 100 MeV rather than 50 (see Charlotte's study)
	if(g1Energy < 0.1 || g2Energy < 0.1)
	//	if(g1Energy < 0.05 || g2Energy < 0.05)
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

      }
    ///now we should have all the candidates..
    //    cout <<"we have " << chargedPiCandidates.size() << " pions " << chargedKCandidates.size() <<" kaons " << pi0Candidates.size() <<" pi0 " << KsCandidates.size() <<" Ks " << leptonCandidates.size() <<"leptons " << otherChargedTracks.size() <<" others " <<endl;


    Mdst_gamma_Manager& gamma_mgr=Mdst_gamma_Manager::get_manager();
    int gammaCount=0;
    for(std::vector<Mdst_gamma>::const_iterator i =gamma_mgr.begin();i!=gamma_mgr.end();i++)
      {
	if(gammaIds.find(i->get_ID())!=gammaIds.end()){
	  continue;
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
	    continue;
	  }
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25 =aux.e9oe25();
	double gammaE=sqrt(px*px+py*py+pz*pz);
	Hep3Vector photVec(px, py,pz);
	v_gammaE.push_back(gammaE);
	//barrel energy cut is the lowest
	if(gammaE<cuts::minGammaEBarrel)
	  {
	    //	    cout <<" loosing photon in barrel due to energy cut " << gammaE<<endl;
	    continue;
	  }
	float photTheta= photVec.theta();
	if(photTheta >cuts::barrelThetaMax) 
	  {
	    if(gammaE<cuts::minGammaEEndcapBkwd)
	      {
		//	    cout <<" loosing photon in forwrad endcap due to energy cut " << gammaE<<endl;
		continue;
	      }
	  }
	if(photTheta< cuts::barrelThetaMin)
	  {
	    if(gammaE<cuts::minGammaEEndcapFwd)
	      {
		//	    cout <<" loosing photon in backward endcap due to energy cut " << gammaE<<endl;
		continue;
	      }
	  }
      }


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
    bool mcDecaySignature=false;
    bool foundRecDecay=false;
    if(m_mc)
      mcDecaySignature=(sig_FoundDDoubleStar && sig_numLeptons==1 && sig_numKaons==0 && (sig_numPions==1 || sig_numPions==2) && sig_numPi0==0 && sig_numBaryons==0);

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
	    double dCharge=(*itD)->charge();
	    //	    cout <<"dcharge: " << dCharge <<endl;
	    vector<Particle*> localPiCandidates;
	    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	      {
		bool isChild=false;
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
		cout << " we should have : "<<leptonCandidates[0]->pType().name() <<" and have: "<<  hepEvt.idhep() <<endl;

		leptonCharge=leptonCandidates[0]->charge();
	      }
	    int numPions=localPiCandidates.size();


	    HepLorentzVector pionMom;
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
	    //	    cout <<"pionCharge: " << pionCharge <<" lepton Charge: "<< leptonCharge <<" sum: "<< pionCharge+leptonCharge+dCharge <<endl;

	    if(PRINT)
	      cout <<"numOtherTracks: "<< numOtherTracks <<" numExtraKaons: "<< numExtraKaons <<" numPions: "<< numPions << "lepton Candidates: "<< leptonCandidates.size() << " system charge: "<< pionCharge+leptonCharge+dCharge <<" best B charge: "<<  bestBcand.charge() <<endl;


	    bool recDecaySignature=numOtherTracks==0&&numExtraKaons==0&&numPions>=1&&numPions<=2 &&leptonCandidates.size()==1 && fabs(pionCharge+leptonCharge+dCharge)<=1.0 && bestBcand.charge()==((-1)*(pionCharge+leptonCharge+dCharge));

	     if(recDecaySignature)
	      {
		foundRecDecay=true;
		if(PRINT)
		  cout <<"all conditions satisfied ... " << endl;
		Particle& D=*(*itD);
		treeData.recDDoubleStar=1;

		//get this by computing B_sl direction from B_tag. B tag is saved in 'bestBCand' particle, so 
		//I assume that the B_sl momentum is just the mirrored momentum;
		float mNu2=0;		
		HepLorentzVector bMomentum=bestBcand.p();
		//	    cout <<"best b cand px: " << -sigBPx << " py: " << -sigBPy << " pz: "<< -sigBPz <<" e: " << sigE <<endl;
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

		mNu2=(kinematics::pBeam-bMomentum-Pxl).mag2();
		if(mNu2<bestMNu2)
		  {
		    bestMNu2=mNu2;
		    bestDIndex=treeData.size;
		  }
		if(PRINT)
		  {

		    if((sigDStarLNu || sigDStarPiLNu || sigDStarPiPiLNu )&& sig_FoundDDoubleStar)
		      {
			cout <<"found dlnu and ddouble star again!!!!!" <<endl;
			cout <<"mNu2: " << mNu2 <<endl;
		      }
		    if(fabs(mNu2) < 1.0 && numPions==2)
		      {

			if(sig_FoundDDoubleStar)
			  {
			  cout <<"smallMNuTwoPionDDoubleStar" <<endl;
			  foundDDStarFlag=true;
			  }
			if(!sigDStarLNu && !sigDStarPiLNu && !sigDStarPiPiLNu && !sigDLNu && !sigDPiLNu && !sigDPiPiLNu)
			  {
			    if((!sig_FoundDDoubleStar || sig_numPions>2 || sig_numKaons!=0 || sig_numBaryons!=0 || sig_numPi0!=0 || sig_numLeptons!=1))
			      {
				cout <<"background event .." <<endl;
				bgFlag=true;
			      }
			  }
			cout <<" small mNu2 (" << mNu2 <<" ),  sigDLNu: "<< sigDLNu << " dpilnu: "<< sigDPiLNu <<" dpipilnu: "<< sigDPiPiLNu<<endl;
			cout <<"Dstarlnu: " << sigDStarLNu <<" dstarpilnu: "<< sigDStarPiLNu <<" dstarpipilnu: " << sigDStarPiPiLNu <<endl;

			cout <<"found any ddouble star? : "<< sig_FoundDDoubleStar <<" num pions: "<< sig_numPions <<" num kaons: ";
			cout <<sig_numKaons <<" num pi0:"<< sig_numPi0 <<" num baryon: " << sig_numBaryons <<endl;
			if(!sigDLNu && !sigDPiLNu && !sigDPiPiLNu && !sigDStarLNu && !sigDStarPiLNu && !sigDStarPiPiLNu)
			  {
			    if(!sig_FoundDDoubleStar || (sig_numPions> 2 || sig_numKaons>0 || sig_numBaryons>0 || sig_numPi0>0|| sig_numLeptons!=1))
			      cout <<"found strange strange background " <<endl;
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
		  }


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
		treeData.mNu2[treeData.size]=mNu2;
		treeData.mXl[treeData.size]=Pxl.mag();
		treeData.mB[treeData.size]=bMomentum.mag();
		treeData.size++;
		//		cout <<"got mNu: " << mNu <<endl;


		//////////
		//		if(treeData.recBToDlNuPi==1 || treeData.recBToDlNuPiPi==1)
		//		if(treeData.foundDPiPi==1 || foundSinglePionDecay)
		if(fabs(mNu2)<0.05)
		  {
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
		    int iTmp;
	    
	    
		    if(PRINT)
		      {
		
		      }
	    
		    if(PRINT)
		      {
			cout <<"check print: " <<endl;
			checkForDPiPi(iTmp,tmp, true);
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
		  }
		


		//////////////

		
	      }

	  }
      }


    //    cout <<"done combining.." <<endl;
    if(bestDIndex>=0 && bestDIndex < 1000)
      {
	treeData.bestD[bestDIndex]=1;
      }
    treeData.recDecaySignature=foundRecDecay;
    treeData.mcDecaySignature=mcDecaySignature;
    if(foundRecDecay)
      {
	cout <<"saving tree, bgFlag: "<< bgFlag <<endl;
	cout <<"saving tree, DD flag: "<<  foundDDStarFlag << endl;
	saveTree();
	cout <<"indeed foundRec " <<endl;
      }
    if(mcDecaySignature&& !foundRecDecay)
      {
	cout <<"saving tree, bgFlag: "<< bgFlag <<endl;
	cout <<"saving tree, DD flag: "<<  foundDDStarFlag << endl;
	//		cout <<"mc decay but not foundRec " <<endl;
	saveTree();
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



  bool bToDDoubleStar::recursivePrint(const Gen_hepevt gen_it, string s)
  {
    genhep_vec* daughters=getDaughters(gen_it);
    int lund=fabs(gen_it.idhep());

    if(lund==911|| lund>9000000 || lund==30343)
      return false;
    Particle p(gen_it);
    if(lund== 100423 || lund ==100421 ||lund==PY_DStar_2S || lund==100411 || lund==100413 )
      cout <<s << lund << endl;
    else
      cout <<s<<p.pType().name() <<" (" << p.pType().lund()<<"), "<<" |p|: " <<p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;
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
  //this method checks if the Dpipi decay is in the MC
  bool bToDDoubleStar::checkForDPiPi(int& bMesonId,bool& foundSinglePionDecay, bool print)
  {
    foundSinglePionDecay=false;
    overlapFractionCharged=.0;
    overlapFractionPi0=.0;

    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();

    //find best b
    float bestDistance=-1;
    Gen_hepevt_Manager::iterator bestB;
    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {
	int geantID=abs(gen_it->idhep());//plus is ok, since it is the abs value
	if(geantID==PY_B0 || geantID==PY_B)
	  {
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
	    bool tempFoundDDoubleStar=false;
	    int tempNumLeptons=0;
	    int tempNumPions=0;
	    int tempNumKaons=0;
	    int tempNumPi0=0;
	    int tempNumBaryons=0;
	    int tempNumD=0;
	    int tempNumDStar=0;
	    bool tempNumDStar2S=false;
	    bool tempNumDStarD2S=false;
	    int tempNumNu=0;
	    //general check if we find one of the D(*)npi l nu decays

	    findDecaySignature(*gen_it,tempFoundDDoubleStar,tempNumLeptons,tempNumPions,tempNumKaons,tempNumPi0,tempNumBaryons,tempNumD, tempNumDStar, tempNumNu,tempNumDStar2S, tempNumDStarD2S);

	    //	    cout <<" get decay sig:numpions: "<< tempNumPions <<endl;

	    //	    cout <<"found " << tempNumD << " Ds " << tempNumNu <<" Neutrinos " << tempNumPions <<" pions " << tempNumLeptons << " leptons " << tempFoundDDoubleStar <<" doublestar " << tempNumDStar2S << " star 2S " << tempNumKaons <<" kaons " << tempNumPi0 <<" pi0s " << tempNumBaryons <<" baryons " << tempNumDStar << " dstar " <<endl;
	    if(tempNumD==1 && tempNumNu==1 && tempNumLeptons==1 && !tempFoundDDoubleStar && !tempNumDStar2S && !tempNumDStarD2S)
	      {
		if(tempNumKaons==0 && tempNumPi0==0 && tempNumBaryons==0 && tempNumDStar==0)
		  {
		    if(tempNumPions<3)
		      {

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

		    if(tempNumPions==0)
		      sigDLNu++;
		    if(tempNumPions==1)
		      sigDPiLNu++;
		    if(tempNumPions==2)
		      sigDPiPiLNu++;

		  }

	      }

    if(tempNumDStar==1 && tempNumNu==1 && tempNumLeptons==1 && !tempFoundDDoubleStar && !tempNumDStar2S && !tempNumDStarD2S)
	      {
		if(tempNumKaons==0 && tempNumPi0==0 && tempNumBaryons==0 && tempNumD==0)
		  {

		    if(tempNumPions==0)
		      sigDStarLNu++;
		    if(tempNumPions==1)
		      sigDStarPiLNu++;
		    if(tempNumPions==2)
		      sigDStarPiPiLNu++;
		  }
		
	      }
	  
	    //haven't found it yet...
	    if(!sig_FoundDDoubleStar)
	      findDecaySignature(*gen_it,sig_FoundDDoubleStar,sig_numLeptons,sig_numPions,sig_numKaons,sig_numPi0,sig_numBaryons,sig_numD, sig_numDStar,tempNumNu, sig_dStar_2S, sig_d_2S);
	    if(!sig_FoundDDoubleStar)
	      {
		//still none, so reset fields
		sig_numLeptons=0;
		sig_numKaons=0;
		sig_numPions=0;
		sig_numPi0=0;
		sig_numBaryons=0;
		sig_numD=0;
		sig_numDStar=0;
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


  //see if this B has any bDoubleSTar and count number of charged pions, kaons, pi0s (so independent of the actual decay)
  bool bToDDoubleStar::findDecaySignature(const Gen_hepevt &mother,bool& dDoubleStar,int& numLeptons, int& numPions, int& numKaons, int& numPi0, int& numBaryons,int& numD, int& numDStar, int& numNu, bool& dStar_2S, bool& d_2S)
  {
    genhep_vec* daughters=getDaughters(mother);
    int lund=fabs(mother.idhep());
    if(lund==911|| lund>9000000)
      return false;
    Particle p(mother);
    //the ones with 20* are D'
    if(lund==PY_DStar_00|| lund==PY_DStar0Plus|| lund==PY_DStar_0s|| lund==PY_DStar0Plus|| lund==PY_D_1 || lund==PY_D_10 || lund==PY_DStar_2 || lund==PY_DStar_20 || lund==100413 || lund == 20423 ||lund ==20413|| lund==100421 || lund==100411)
      {
	dDoubleStar=true;
      }

    if(lund==100421|| lund==100411)
      d_2S=true;
    if(lund==100423|| lund==100413)
      dStar_2S=true;

    if(lund==PY_D || lund==PY_D0)
      {
	numD++;
	return true;
      }
    if(lund==PY_DStar || lund==PY_DStar0)
      {
	//don't add D decay products to number of pions etc...
	numDStar++;
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
	findDecaySignature(**it,dDoubleStar,numLeptons,numPions,numKaons,numPi0,numBaryons,numD, numDStar,numNu,dStar_2S, d_2S);
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
		if(m>1.0)
		  mesonMassTree->Fill();
	      }


	    Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));
	    //	if(!doKmVtxFit2(*(*itD),  confLevel,0))


	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;

	    d0->relation().append(kaon);
	    d0->relation().append(pion);
	    if(m_mc)
	      {
		const Gen_hepevt &h_kaon=kaon.genHepevt();
		const Gen_hepevt &h_pion=pion.genHepevt();
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
		    if(m>1.0)
		      mesonMassTree->Fill();
		  }

		if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
		Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));
		d0->relation().append(kaon);
		d0->relation().append(pion);
		d0->relation().append(pi0);
		if(m_mc)
		  {
		    const Gen_hepevt &h_kaon=kaon.genHepevt();
		    const Gen_hepevt &h_pion=pion.genHepevt();
		    const Gen_hepevt &h_pi0=pi0.genHepevt();
		    if(h_pi0&&h_kaon && h_pion && h_kaon.mother() && h_pion.mother() && h_kaon.mother().get_ID()==h_pion.mother().get_ID() && h_pi0.mother() && h_pi0.mother()==h_kaon.mother()){
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
			if(m>1.0)
			  mesonMassTree->Fill();
		      }
		    histoRecD0Spect->Fill(m);
		    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
		    Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));
		    d0->relation().append(kaon);
		    d0->relation().append(pion1);
		    d0->relation().append(pion2);
		    d0->relation().append(pion3);
		    if(m_mc)
		      {
			const Gen_hepevt &h_kaon=kaon.genHepevt();
			const Gen_hepevt &h_pion1=pion1.genHepevt();
			const Gen_hepevt &h_pion2=pion2.genHepevt();
			const Gen_hepevt &h_pion3=pion3.genHepevt();
		    
			if(h_kaon && h_pion1 && h_pion2&&h_pion3){
			  if(h_kaon.mother() && h_pion1.mother() &&h_pion2.mother() && h_pion3.mother()){
			    if( h_kaon.mother().get_ID()==h_pion1.mother().get_ID() && h_pion1.mother()==h_pion2.mother() && h_pion2.mother()==h_pion2.mother()){
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
		if(pion1.charge()+pion2.charge()!=0) continue;
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
		    if(m>1.0)
		      mesonMassTree->Fill();
		  }
		histoRecD0Spect->Fill(m);
		if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));
		Particle* d0 =new Particle(p_d0,Ptype("D0"));
		d0->relation().append(pion1);
		d0->relation().append(pion2);
		d0->relation().append(Ks);
		if(m_mc)
		  {
		    const Gen_hepevt &h_pion1=pion1.genHepevt();
		    const Gen_hepevt &h_pion2=pion2.genHepevt();
		    const Gen_hepevt &h_Ks=Ks.genHepevt();
		    if(h_pion1&&h_pion2 && h_Ks&& h_pion1.mother() && h_pion2.mother() && h_Ks.mother() && h_pion1.mother().get_ID()==h_pion2.mother().get_ID() && h_pion1.mother()==h_Ks.mother()){
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
		if(m>1.0)
		  mesonMassTree->Fill();
	      }
	    histoRecD0Spect->Fill(m);
	    //	    cout <<"found k/k combination, filling with m: "<< m <<endl;
	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
	    Particle* d0 =new Particle(p_d0,Ptype("D0"));
	    d0->relation().append(kaon1);
	    d0->relation().append(kaon2);
	    if(m_mc)
	      {
		const Gen_hepevt &h_kaon1=kaon1.genHepevt();
		const Gen_hepevt &h_kaon2=kaon2.genHepevt();
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
		if(m>1.0)
		  mesonMassTree->Fill();
	      }
	    histoRecD0Spect->Fill(m);
	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
	    Particle* d0 =new Particle(p_d0,Ptype( "D0"));
	    d0->relation().append(pi0);
	    d0->relation().append(Ks);
	    if(m_mc)
	      {
		const Gen_hepevt &h_pi0=pi0.genHepevt();
		const Gen_hepevt &h_Ks=Ks.genHepevt();
		if(h_pi0 && h_Ks && h_pi0.mother() && h_Ks.mother() && h_pi0.mother().get_ID()==h_Ks.mother().get_ID()){
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
		if(m>1.0)
		  mesonMassTree->Fill();
	      }
	    histoRecDSpect->Fill(m);
	      
	    if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
	    Particle* dPlus =new Particle(p_dPlus,Ptype(pion.charge() > 0 ? "D+" : "D-"));
	    dPlus->relation().append(Ks);
	    dPlus->relation().append(pion);
	    if(m_mc)
	      {
		const Gen_hepevt &h_Ks=Ks.genHepevt();
		const Gen_hepevt &h_pion=pion.genHepevt();
		if(h_Ks && h_pion && h_Ks.mother() && h_pion.mother() && h_Ks.mother().get_ID()==h_pion.mother().get_ID()){
		  dPlus->relation().genHepevt(h_Ks.mother());
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
			if(m>1.0)
			  mesonMassTree->Fill();
		      }
		    histoRecDSpect->Fill(m);
		    if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		    Particle* dPlus =new Particle(p_dPlus,Ptype(charge >0 ? "D+" : "D-"));
		    dPlus->relation().append(Ks);
		    dPlus->relation().append(pion1);
		    dPlus->relation().append(pion2);
		    dPlus->relation().append(pion3);
		    if(m_mc)
		      {
			const Gen_hepevt &h_Ks=Ks.genHepevt();
			const Gen_hepevt &h_pion1=pion1.genHepevt();
			const Gen_hepevt &h_pion2=pion2.genHepevt();
			const Gen_hepevt &h_pion3=pion3.genHepevt();
		    
			if(h_Ks && h_pion1 && h_pion2&&h_pion3){
			  if(h_Ks.mother() && h_pion1.mother() &&h_pion2.mother() && h_pion3.mother()){
			    if( h_Ks.mother().get_ID()==h_pion1.mother().get_ID() && h_pion1.mother()==h_pion2.mother() && h_pion2.mother()==h_pion2.mother()){
			      dPlus->relation().genHepevt(h_pion1.mother());
			    }
			  }
			}
		      }
		    chargedDCandidates.push_back(dPlus);
		  }
	      }
	  }
      }


    //D-->K-P+P+
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
		    if(m>1.0)
		      mesonMassTree->Fill();
		  }
		histoRecDSpect->Fill(m);
		if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));
		Particle* dPlus =new Particle(p_dPlus,Ptype(charge > 0 ? "D+": "D-"));
		dPlus->relation().append(pion1);
		dPlus->relation().append(pion2);
		dPlus->relation().append(kaon);
		if(m_mc)
		  {
		    const Gen_hepevt &h_pion1=pion1.genHepevt();
		    const Gen_hepevt &h_pion2=pion2.genHepevt();
		    const Gen_hepevt &h_kaon=kaon.genHepevt();
		    if(h_pion1&&h_pion2 && h_kaon&& h_pion1.mother() && h_pion2.mother() && h_kaon.mother() && h_pion1.mother().get_ID()==h_pion2.mother().get_ID() && h_pion1.mother()==h_kaon.mother()){
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
		    if(m>1.0)
		      mesonMassTree->Fill();
		  }
		histoRecDSpect->Fill(m);
		if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));
		Particle* dPlus =new Particle(p_dPlus,Ptype(charge > 0 ? "D+": "D-"));
		dPlus->relation().append(kaon1);
		dPlus->relation().append(kaon2);
		dPlus->relation().append(pion);
		if(m_mc)
		  {
		    const Gen_hepevt &h_kaon1=kaon1.genHepevt();
		    const Gen_hepevt &h_kaon2=kaon2.genHepevt();
		    const Gen_hepevt &h_pion=pion.genHepevt();
		    if(h_kaon1&&h_kaon2 && h_pion&& h_kaon1.mother() && h_kaon2.mother() && h_pion.mother() && h_kaon1.mother().get_ID()==h_kaon2.mother().get_ID() && h_kaon1.mother()==h_pion.mother()){
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
		if(m>1.0)
		  mesonMassTree->Fill();
	      }
	    histoRecDStarSpect->Fill(m);
	    histoRecDStarSpectToDPi0->Fill(m);
	    if(m>m_dStarPlusmass_max || m < m_dStarPlusmass_min ||isnan(m)) continue;
	    //	    cout <<"looking at dstar, mass diff: " <<(m_DStarPlus-m_DPlus) <<" vs : " << (m-D.p().mag());
	    //	    cout <<" gives: " << fabs(m-D.p().mag()-(m_DStarPlus-m_DPlus)) <<endl;

	    if(fabs(m-D.p().mag()-(m_DStarPlus-m_DPlus)) > max_massDifference) continue;
	    //	    cout <<"done " <<endl;
	    Particle* dStar =new Particle(p_dStar,Ptype(D.charge()>0 ? "D*+" : "D*-"));
	    dStar->relation().append(D);
	    dStar->relation().append(pi0);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Gen_hepevt &h_pi0=pi0.genHepevt();
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
	    //make sure that pion is not child of D
	    for(int i =0;i<D.nChildren();i++)
	      {
		if(D.child(i).relation().isIdenticalWith(pion.relation()))
		  {
		    //		    cout <<"found double use 2 .." <<endl;
		    doubleUse=true;
		    break;
		  }
	      }
	    if(doubleUse)
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
		if(m>1.0)
		  mesonMassTree->Fill();
	      }
	    histoRecDStarSpect->Fill(m);
	    histoRecDStarSpectToD0Pi->Fill(m);
	    if(m>m_dStarPlusmass_max || m < m_dStarPlusmass_min ||isnan(m)) continue;
	    //	    cout <<"m -D: "<< m-D.p().mag() <<endl;
	    //	    cout <<"looking at dstar, mass diff: " <<(m_DStarPlus-m_D0) <<" vs : " << (m-D.p().mag());
	    //	    cout <<" gives: " << fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) <<endl;
	    if(fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) > max_massDifference) continue;
	    //	    cout <<"done" <<endl;
	    Particle* dStar =new Particle(p_dStar,Ptype(pion.charge()>0 ? "D*+" : "D*-"));
	    dStar->relation().append(D);
	    dStar->relation().append(pion);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Gen_hepevt &h_pion=pion.genHepevt();
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
	    dDecay=3;
	    if(foundDPiPi || foundSinglePionDecay)
	      foundDDoubleStarDecay=1;
	    else
	      foundDDoubleStarDecay=0;
	    if(SAVE_MESON_MASS_DISTRIBUTIONS)
	      {
		if(m>1.0)
		  mesonMassTree->Fill();
	      }
	    histoRecDStarSpect->Fill(m);
	    histoRecDStarSpectToD0Pi0->Fill(m);
	    if(m>m_dStar0mass_max || m < m_dStar0mass_min ||isnan(m)) continue;
	    //	    cout <<"looking at dstar, mass diff: " <<(m_DStar0-m_D0) <<" vs : " << (m-D.p().mag());
	    //	    cout <<" gives: " << fabs(m-D.p().mag()-(m_DStar0-m_D0)) <<endl;
	    if(fabs(m-D.p().mag()-(m_DStar0-m_D0)) > max_massDifference) continue;
	    //	    cout <<" done " <<endl;
	    Particle* dStar =new Particle(p_dStar,Ptype("D*0"));
	    dStar->relation().append(D);
	    dStar->relation().append(pi0);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Gen_hepevt &h_pi0=pi0.genHepevt();
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
  

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
