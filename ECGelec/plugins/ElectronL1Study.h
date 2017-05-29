#ifndef DEF_ElectronL1Study
#define DEF_ElectronL1Study

#include "EGamma/ECGelec/plugins/ToInclude.h"
#include "EGamma/ECGelec/interface/CommonFunctions.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace IPTools;

class MultiTrajectoryStateMode ;
class EgammaTowerIsolation ;

class ElectronL1Study : public edm::EDAnalyzer {
 public:
  explicit ElectronL1Study(const edm::ParameterSet&);
  ~ElectronL1Study();
	
  typedef math::XYZTLorentzVector LorentzVector ;
	
 private:
  virtual void beginJob(const edm::ParameterSet& conf) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
	
  void Init();
	
  void FillEvent (const edm::Event&, const edm::EventSetup&);
  void FillTrigger (const edm::Event&, const edm::EventSetup&);
  void FillSpikes (const edm::Event&, const edm::EventSetup&);
  void FillEle (const edm::Event&, const edm::EventSetup&);
  void FillSuperClusters(const edm::Event&, const edm::EventSetup&);
  const EcalRecHit getRecHit(DetId id, const EcalRecHitCollection *recHits);

  ///////////
  // SETUP //
  ///////////

  const MultiTrajectoryStateTransform *mtsTransform_;	
  unsigned long long cacheIDTDGeom_;
  unsigned long long cacheIDMagField_;
  edm::ESHandle<MagneticField> theMagField_;
  edm::ESHandle<TrackerGeometry> trackerHandle_;	
		
  const CaloSubdetectorGeometry * theEndcapGeometry_ ;
  const CaloSubdetectorGeometry * theBarrelGeometry_ ;
  const CaloTopology * topology_ ;
  edm::ESHandle<CaloGeometry> theCaloGeom_;
  edm::ESHandle<EcalTrigTowerConstituentsMap> eTTmap_;

  //ADD
  // BeamSpot
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;

  // H/E
  edm::Handle<CaloTowerCollection> * towersH_ ;
  edm::InputTag hcalTowers_ ;
  EgammaTowerIsolation * towerIso1_ ;
  EgammaTowerIsolation * towerIso2_ ;
  double hOverEConeSize_ ;
  double hOverEPtMin_ ;            

  
  //ROKO matching for upgrade
 // for track association / propagation to calo
//   TrackAssociatorParameters parameters_;
//   TrackDetectorAssociator trackAssociator_;;


  //////////////////////
  // OUTPUT VARIABLES //
  //////////////////////

  TTree *mytree_;  
  int nEvent, nRun, nLumi;
  
  // Vertices //
  int _vtx_N;
  double _vtx_x[200], _vtx_y[200], _vtx_z[200];
  double _vtx_normalizedChi2[200], _vtx_ndof[200], _vtx_nTracks[200], _vtx_d0[200];
  GlobalPoint vertexPosition;

  // Trigger Paths //
  int trig_hltInfo[250];
  int _trig_isEleHLTpath;
  int trig_HLT_path[4]; // unbias, EG5, EG8, EG12
  char trig_fired_names[5000];
  
  vector<string> m_HLT_pathsV;
  vector<string> m_HLT_triggered;
  vector<int> m_HLT_pathsV_check;

  // Masking //
  int _trig_nMaskedRCT, _trig_nMaskedCh;
  int _trig_iMaskedRCTeta[100], _trig_iMaskedRCTphi[100], _trig_iMaskedRCTcrate[100], _trig_iMaskedTTeta[100], _trig_iMaskedTTphi[100];

  int _trig_strip_mask_N;
  int _trig_strip_mask_TTieta[1000], _trig_strip_mask_TTiphi[1000], _trig_strip_mask_status[1000],
    _trig_strip_mask_StripID[1000], _trig_strip_mask_PseudoStripID[1000], _trig_strip_mask_TccID[1000], _trig_strip_mask_CCU[1000],
    _trig_strip_mask_xtal_ix[1000][5], _trig_strip_mask_xtal_iy[1000][5], _trig_strip_mask_xtal_iz[1000][5];

  int _trig_xtal_mask_N; // [EB+EE]
  int _trig_xtal_mask_ieta[1000],_trig_xtal_mask_iphi[1000], // for EE : xtal ieta->ix ; iphi -> iy
    _trig_xtal_mask_TTieta[1000],_trig_xtal_mask_TTiphi[1000], // but for EE towers, still ieta, iphi...
    _trig_xtal_mask_Rieta[1000],_trig_xtal_mask_Riphi[1000],
    _trig_xtal_mask_status[1000], _trig_xtal_mask_EBEE[1000]; // EBEE = {0,1} => 0=EB ; 1=EE
  //double _trig_xtal_mask_eT[1000];

  // Trigger primitives //
  //
  // handles to get the TPs
  edm::Handle<EcalTrigPrimDigiCollection> * ecal_tp_;
  edm::Handle<EcalTrigPrimDigiCollection> * ecal_tpM_;
  EcalTrigTowerDetId TPtowid_;
  EcalTrigTowerDetId TPtowidM_;
  //
  // original TP
  int _trig_tower_N, _trig_tower_ieta[4032],_trig_tower_iphi[4032],_trig_tower_adc[4032], _trig_tower_sFGVB[4032]; 
  //
  // cleaned TP
  int _trig_tower_N_M, _trig_tower_ieta_M[4032],_trig_tower_iphi_M[4032],
    _trig_tower_adc_M[4032], _trig_tower_sFGVB_M[4032]; 
  //
  // emulated TP
  int _trig_tower_N_E, _trig_tower_ieta_E[4032],_trig_tower_iphi_E[4032],
    _trig_tower_adc_E[4032][5], _trig_tower_sFGVB_E[4032][5]; 

  // HCAL TP
  int _trig_tower_hcal_N, _trig_tower_hcal_ieta[4032], _trig_tower_hcal_iphi[4032], _trig_tower_hcal_FG[4032],_trig_tower_hcal_et[4032];

  // Level-1 trigger //
  int _trig_L1emIso_N; 
  int _trig_L1emNonIso_N;
  int _trig_L1emIso_ieta[4], _trig_L1emIso_iphi[4], _trig_L1emIso_rank[4]; 
  double _trig_L1emIso_eta[4], _trig_L1emIso_phi[4],_trig_L1emIso_energy[4],_trig_L1emIso_et[4]; 
  int _trig_L1emNonIso_ieta[4], _trig_L1emNonIso_iphi[4],_trig_L1emNonIso_rank[4];
  double _trig_L1emNonIso_eta[4], _trig_L1emNonIso_phi[4], _trig_L1emNonIso_energy[4],_trig_L1emNonIso_et[4];
  //
  // modified collection
  int _trig_L1emIso_N_M; 
  int _trig_L1emNonIso_N_M;
  int _trig_L1emIso_ieta_M[4], _trig_L1emIso_iphi_M[4], _trig_L1emIso_rank_M[4]; 
  double _trig_L1emIso_eta_M[4], _trig_L1emIso_phi_M[4],_trig_L1emIso_energy_M[4],_trig_L1emIso_et_M[4]; 
  int _trig_L1emNonIso_ieta_M[4], _trig_L1emNonIso_iphi_M[4],_trig_L1emNonIso_rank_M[4];
  double _trig_L1emNonIso_eta_M[4], _trig_L1emNonIso_phi_M[4], _trig_L1emNonIso_energy_M[4],_trig_L1emNonIso_et_M[4];
  
  //SLHC collection
  int _trig_L1emIso_N_SLHC; 
  int _trig_L1emNonIso_N_SLHC;
  int _trig_L1emIso_ieta_SLHC[50], _trig_L1emIso_iphi_SLHC[50], _trig_L1emIso_rank_SLHC[50]; 
  double _trig_L1emIso_eta_SLHC[50], _trig_L1emIso_phi_SLHC[50],_trig_L1emIso_energy_SLHC[50],_trig_L1emIso_et_SLHC[50]; 
  int _trig_L1emNonIso_ieta_SLHC[50], _trig_L1emNonIso_iphi_SLHC[50],_trig_L1emNonIso_rank_SLHC[50];
  double _trig_L1emNonIso_eta_SLHC[50], _trig_L1emNonIso_phi_SLHC[50], _trig_L1emNonIso_energy_SLHC[50],_trig_L1emNonIso_et_SLHC[50];

  //modif-alex.
  int _trig_L1tauIso_N_SLHC;
  int _trig_L1tauNonIso_N_SLHC;
  int _trig_L1tauIso_ieta_SLHC[50], _trig_L1tauIso_iphi_SLHC[50], _trig_L1tauIso_rank_SLHC[50];
  double _trig_L1tauIso_eta_SLHC[50], _trig_L1tauIso_phi_SLHC[50],_trig_L1tauIso_energy_SLHC[50],_trig_L1tauIso_et_SLHC[50];
  int _trig_L1tauNonIso_ieta_SLHC[50], _trig_L1tauNonIso_iphi_SLHC[50],_trig_L1tauNonIso_rank_SLHC[50];
  double _trig_L1tauNonIso_eta_SLHC[50], _trig_L1tauNonIso_phi_SLHC[50], _trig_L1tauNonIso_energy_SLHC[50],_trig_L1tauNonIso_et_SLHC[50];


  int  _ele_L1Noniso_SLHC[10],_ele_L1Iso_SLHC[10], _ele_L1Noniso_matchInd_SLHC[10],_ele_L1Iso_matchInd_SLHC[10],
        _ele_L1Iso_SLHC_TTeta[10],_ele_L1Iso_SLHC_TTphi[10],  _ele_L1NonIso_SLHC_TTeta[10],_ele_L1NonIso_SLHC_TTphi[10];

  //modif-alex
  int  _ele_L1TauNoniso_SLHC[10],_ele_L1TauIso_SLHC[10], _ele_L1TauNoniso_matchInd_SLHC[10],_ele_L1TauIso_matchInd_SLHC[10],
        _ele_L1TauIso_SLHC_TTeta[10],_ele_L1TauIso_SLHC_TTphi[10],  _ele_L1TauNonIso_SLHC_TTeta[10],_ele_L1TauNonIso_SLHC_TTphi[10];
  
  
  //
  // L1 prefiring
  int _trig_preL1emIso_N; 
  int _trig_preL1emNonIso_N;
  int _trig_preL1emIso_ieta[4], _trig_preL1emIso_iphi[4], _trig_preL1emIso_rank[4]; 
  int _trig_preL1emNonIso_ieta[4], _trig_preL1emNonIso_iphi[4],_trig_preL1emNonIso_rank[4];
  //
  // L1 postfiring
  int _trig_postL1emIso_N; 
  int _trig_postL1emNonIso_N;
  int _trig_postL1emIso_ieta[4], _trig_postL1emIso_iphi[4], _trig_postL1emIso_rank[4]; 
  int _trig_postL1emNonIso_ieta[4], _trig_postL1emNonIso_iphi[4],_trig_postL1emNonIso_rank[4];
	
  // Spikes //
  int spike_N, spike_SwissCross[5000], spike_TTieta[5000], spike_TTiphi[5000], spike_Rieta[5000], spike_Riphi[5000], spike_severityLevel[5000], spike_outOfTime[5000];
  double  spike_Et[5000], spike_eta[5000], spike_phi[5000], spike_theta[5000], spike_time[5000];

  // Electrons //
  TClonesArray * m_electrons ;
  TLorentzVector myvector ;	
  //
  int ele_N, ele_nSeed;
  int ele_echarge[10], ele_sclNclus[10];
  double ele_he[10] , 
    ele_eseedpout[10] , ele_ep[10] , ele_eseedp[10] , ele_eelepout[10] ,       
    ele_deltaetaseed[10] , ele_deltaetaele[10] , ele_deltaphiseed[10] , ele_deltaphiele[10] , ele_deltaetain[10] , ele_deltaphiin[10] ,
    ele_sigmaietaieta[10] , ele_sigmaetaeta[10] ,
    ele_pin_mode[10] , ele_pout_mode[10] , ele_pin_mean[10] , ele_pout_mean[10] , 
    ele_pTin_mode[10] , ele_pTout_mode[10] , ele_pTin_mean[10] , ele_pTout_mean[10],
    ele_fbrem[10], ele_mva[10], ele_calo_energy[10],
    ele_sclE[10], ele_sclEt[10], ele_sclEta[10], ele_sclPhi[10] ;
  //
  int ele_isbarrel[10] , ele_isendcap[10] , ele_isEBetaGap[10] , 
    ele_isEBphiGap[10] , ele_isEEdeeGap[10] , ele_isEEringGap[10] ,
    ele_eClass[10], ele_isecalDriven[10] , ele_istrackerDriven[10] ;
  double ele_vertex_x[10], ele_vertex_y[10], ele_vertex_z[10];
  int ele_missing_hits[10], ele_lost_hits[10]; 
  double ele_chi2_hits[10];
  double ele_track_x[10], ele_track_y[10], ele_track_z[10];
  double ele_tkSumPt_dr03[10] , ele_ecalRecHitSumEt_dr03[10] , ele_hcalDepth1TowerSumEt_dr03[10] , ele_hcalDepth2TowerSumEt_dr03[10],
    ele_tkSumPt_dr04[10] , ele_ecalRecHitSumEt_dr04[10] , ele_hcalDepth1TowerSumEt_dr04[10] , ele_hcalDepth2TowerSumEt_dr04[10] ;
  int ele_severityLevelSeed[10], ele_severityLevelClusters[10];
  double ele_ECAL_fbrem[10];
  //
  // L1 trigger
  int _ele_TTetaVect[10][50], _ele_TTphiVect[10][50];
  double _ele_TTetVect[10][50];
  int _ele_RCTetaVect[10][10], _ele_RCTphiVect[10][10];
  int _ele_RCTL1isoVect[10][10], _ele_RCTL1nonisoVect[10][10];
  int _ele_RCTL1isoVect_M[10][10], _ele_RCTL1nonisoVect_M[10][10];
  double _ele_RCTetVect[10][10];
  //
  int _ele_RCTeta[10], _ele_RCTphi[10];
  int _ele_RCTL1noniso[10], _ele_RCTL1iso[10];
  int _ele_RCTL1noniso_M[10], _ele_RCTL1iso_M[10];
  //
  // Conversion removal
  int ele_isConversion[10];
  int ele_convFound[10];
  double ele_conv_dcot[10];
  double ele_conv_dist[10];
  int ele_expected_inner_hits[10];
  int ele_expected_inner_hits_aod[10];

  // SuperClusters //
  int _sc_N[2]; 
  int _sc_EB_EE[50];
  double _sc_E[50], _sc_Et[50], _sc_Eta[50], _sc_Phi[50]; 
  int _sc_severityLevelSeed[50];
  double _sc_he[50], _sc_sigmaietaieta[50];
  double _sc_hcalDepth1TowerSumEt_dr03[50], _sc_hcalDepth2TowerSumEt_dr03[50];
  double _sc_ecalRecHitSumEt_dr03[50];
  double _sc_trkiso_dr03[50];
  //
  int _sc_RCTeta[50];
  int _sc_RCTphi[50];
  int _sc_RCTL1iso[50];
  int _sc_RCTL1noniso[50];
  int _sc_RCTL1iso_M[50];
  int _sc_RCTL1noniso_M[50];
  int _sc_TTetaVect[50][50], _sc_TTphiVect[50][50];
  double _sc_TTetVect[50][50];
  int _sc_RCTetaVect[50][10], _sc_RCTphiVect[50][10];
  int _sc_RCTL1isoVect[50][10], _sc_RCTL1nonisoVect[50][10];
  int _sc_RCTL1isoVect_M[50][10], _sc_RCTL1nonisoVect_M[50][10];
  double _sc_RCTetVect[50][10];


  //////////////////////////////////
  //////// INPUT PARAMETERS ////////
  //////////////////////////////////

  // Parameters //
  bool GetL1_ ; // get L1 candidates from data
  bool GetL1M_ ; // get emulated L1 candidates
  bool GetSLHC_L1_;
  bool GetTP_ ; // get the standard trigger primitives
  bool GetTP_Modif_ ; // get the modified collection of trigger primitives (zeroing by hand)
  bool GetTP_Emul_ ; // get the emulated collection of trigger primitives (Jackson-Zabi's sFGVB+zeroing emulator)
  bool GetHcalTP_ ;
  bool GetStripMask_ ;
  bool GetXtalMask_ ;
  bool GetVertices_ ;
  bool PrintDebug_ , PrintDebug_HLT_;
  bool DoFillEle_;
  bool DoFillTrigger_;
  bool DoFillSC_;
  bool DoFillSpikes_;
  bool aod_;	
  std::string type_;

  // Input Tags //
  edm::InputTag EcalRecHitCollectionEB_ ;
  edm::InputTag tpCollectionNormal_ ;
  edm::InputTag tpCollectionModif_ ;
  edm::InputTag tpEmulatorCollection_ ;
  edm::InputTag hcalTPGs_ ;
  edm::InputTag EleTag_;
  edm::InputTag VerticesTag_;
  edm::InputTag dcsTag_;
  edm::InputTag HLTTag_; 
  std::vector<std::string > HLT_ElePaths_;
  std::vector<std::string > HLT_paths_;
  std::string gtRecordCollectionTag_ ;
};

#endif
