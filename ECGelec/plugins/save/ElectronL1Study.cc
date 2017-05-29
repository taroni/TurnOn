#include "EGamma/ECGelec/plugins/ElectronL1Study.h"
#include "EGamma/ECGelec/plugins/CommonFunctions.h"

// ====================================================================================
ElectronL1Study::ElectronL1Study(const edm::ParameterSet& iConfig)
{
  //////////////////////
  // INPUT PARAMETERS //
  //////////////////////

  // RecHits (to get spikes)
  EcalRecHitCollectionEB_ = iConfig.getParameter<edm::InputTag>("EcalRecHitCollectionEB") ;
  // Trigger primitives to get
  GetTP_       = iConfig.getUntrackedParameter<bool>("GetTP");
  GetTP_Modif_ = iConfig.getUntrackedParameter<bool>("GetTPmodif");
  GetTP_Emul_  = iConfig.getUntrackedParameter<bool>("GetTPemul");
  // Trigger primitive tags
  tpCollectionNormal_     = iConfig.getParameter<edm::InputTag> ("TPCollectionNormal") ;
  tpCollectionModif_      = iConfig.getParameter<edm::InputTag> ("TPCollectionModif") ;
  tpEmulatorCollection_   = iConfig.getParameter<edm::InputTag> ("TPEmulatorCollection") ;
  // Compute H/E
  hcalTowers_ = iConfig.getParameter<edm::InputTag>("hcalTowers");
  hOverEPtMin_ = iConfig.getParameter<double>("hOverEPtMin");
  hOverEConeSize_ = iConfig.getParameter<double>("hOverEConeSize");

  // Level-1 candidates
  GetL1_       = iConfig.getUntrackedParameter<bool>("GetL1");
  GetL1M_      = iConfig.getUntrackedParameter<bool>("GetL1M");
  gtRecordCollectionTag_ = iConfig.getParameter<std::string>("GTRecordCollection"); // for pre/post firing
  // Objects tags
  EleTag_ = iConfig.getParameter<edm::InputTag> ("EleTag");
  VerticesTag_ = iConfig.getParameter<edm::InputTag> ("VerticesTag");
  // Printouts
  PrintDebug_  = iConfig.getUntrackedParameter<bool>("PrintDebug");
  PrintDebug_HLT_ = iConfig.getUntrackedParameter<bool>("PrintDebug_HLT");  
  // Processing
  DoFillEle_     = iConfig.getUntrackedParameter<bool>("DoFillEle");
  DoFillTrigger_ = iConfig.getUntrackedParameter<bool>("DoFillTrigger");
  DoFillSC_      = iConfig.getUntrackedParameter<bool>("DoFillSC");
  DoFillSpikes_  = iConfig.getUntrackedParameter<bool>("DoFillSpikes");
  aod_           = iConfig.getUntrackedParameter<bool>("AOD");
  // HLT
  HLTTag_       = iConfig.getParameter<edm::InputTag> ("HLTTag");
  HLT_ElePaths_ = iConfig.getParameter<std::vector<std::string > >("HLTElePaths");
  HLT_paths_ = iConfig.getParameter<std::vector<std::string > >("HLT_paths");
  // Configuration
  dcsTag_       = iConfig.getUntrackedParameter<edm::InputTag>("dcsTag");
  

  /////////////////
  // OUTPUT TREE //
  /////////////////
  edm::Service<TFileService> fs ;
  mytree_  = fs->make <TTree>("eIDSimpleTree","eIDSimpleTree"); 

  // Global
  mytree_->Branch("nEvent",&nEvent,"nEvent/I");
  mytree_->Branch("nRun",&nRun,"nRun/I");
  mytree_->Branch("nLumi",&nLumi,"nLumi/I");
	
  // Vertices
  mytree_->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  mytree_->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[200]/D");
  mytree_->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[200]/D");
  mytree_->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[200]/D");
  mytree_->Branch("vtx_d0",&_vtx_d0,"vtx_d0[200]/D");
  mytree_->Branch("vtx_x",&_vtx_x,"vtx_x[200]/D");
  mytree_->Branch("vtx_y",&_vtx_y,"vtx_y[200]/D");
  mytree_->Branch("vtx_z",&_vtx_z,"vtx_z[200]/D");
	
  // Spikes
  mytree_->Branch("spike_N",&spike_N,"spike_N/I");
  mytree_->Branch("spike_TTieta",&spike_TTieta,"spike_TTieta[5000]/I");
  mytree_->Branch("spike_TTiphi",&spike_TTiphi,"spike_TTiphi[5000]/I");
  mytree_->Branch("spike_Rieta",&spike_Rieta,"spike_Rieta[5000]/I");
  mytree_->Branch("spike_Riphi",&spike_Riphi,"spike_Riphi[5000]/I");
  mytree_->Branch("spike_SwissCross",&spike_SwissCross,"spike_SwissCross[5000]/I");
  mytree_->Branch("spike_severityLevel",&spike_severityLevel,"spike_severityLevel[5000]/I");
  mytree_->Branch("spike_outOfTime",&spike_outOfTime,"spike_outOfTime[5000]/I");
  mytree_->Branch("spike_time", &spike_time,"spike_time[5000]/D");
  mytree_->Branch("spike_Et",&spike_Et,"spike_Et[5000]/D");
  mytree_->Branch("spike_eta",&spike_eta,"spike_eta[5000]/D");
  mytree_->Branch("spike_phi",&spike_phi,"spike_phi[5000]/D");
  mytree_->Branch("spike_theta",&spike_theta,"spike_theta[5000]/D");

  // Towers (original collection)
  mytree_->Branch("trig_tower_N", &_trig_tower_N, "trig_tower_N/I");
  mytree_->Branch("trig_tower_ieta",  &_trig_tower_ieta,  "trig_tower_ieta[trig_tower_N]/I");
  mytree_->Branch("trig_tower_iphi",  &_trig_tower_iphi,  "trig_tower_iphi[trig_tower_N]/I");
  mytree_->Branch("trig_tower_adc",  &_trig_tower_adc,  "trig_tower_adc[trig_tower_N]/I");
  mytree_->Branch("trig_tower_sFGVB",  &_trig_tower_sFGVB,  "trig_tower_sFGVB[trig_tower_N]/I");
	
  // Towers (cleaned collection)
  mytree_->Branch("trig_tower_N_M", &_trig_tower_N_M, "trig_tower_N_M/I");
  mytree_->Branch("trig_tower_ieta_M",  &_trig_tower_ieta_M,  "trig_tower_ieta_M[trig_tower_N_M]/I");
  mytree_->Branch("trig_tower_iphi_M",  &_trig_tower_iphi_M,  "trig_tower_iphi_M[trig_tower_N_M]/I");
  mytree_->Branch("trig_tower_adc_M",  &_trig_tower_adc_M,  "trig_tower_adc_M[trig_tower_N_M]/I");
  mytree_->Branch("trig_tower_sFGVB_M",  &_trig_tower_sFGVB_M,  "trig_tower_sFGVB_M[trig_tower_N_M]/I");
	
  // Towers (emulated)
  mytree_->Branch("trig_tower_N_E", &_trig_tower_N_E, "trig_tower_N_E/I");
  mytree_->Branch("trig_tower_ieta_E",  &_trig_tower_ieta_E,  "trig_tower_ieta_E[trig_tower_N_E]/I");
  mytree_->Branch("trig_tower_iphi_E",  &_trig_tower_iphi_E,  "trig_tower_iphi_E[trig_tower_N_E]/I");
  mytree_->Branch("trig_tower_adc_E",  &_trig_tower_adc_E,  "trig_tower_adc_E[trig_tower_N_E][5]/I");
  mytree_->Branch("trig_tower_sFGVB_E",  &_trig_tower_sFGVB_E,  "trig_tower_sFGVB_E[trig_tower_N_E][5]/I");
		
  // Trigger
  mytree_->Branch("trig_HLT_path",&trig_HLT_path,"trig_HLT_path[4]/I");
  // unbias, EG5, EG8, EG12
  //
  mytree_->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[5000]/C");
  mytree_->Branch("trig_hltInfo",&trig_hltInfo,"trig_hltInfo[250]/I");
  //
  mytree_->Branch("trig_L1emIso_N",     &_trig_L1emIso_N,     "trig_L1emIso_N/I");
  mytree_->Branch("trig_L1emIso_ieta",  &_trig_L1emIso_ieta,  "trig_L1emIso_ieta[4]/I");
  mytree_->Branch("trig_L1emIso_iphi",  &_trig_L1emIso_iphi,  "trig_L1emIso_iphi[4]/I");
  mytree_->Branch("trig_L1emIso_rank",  &_trig_L1emIso_rank,  "trig_L1emIso_rank[4]/I");
  mytree_->Branch("trig_L1emIso_eta",   &_trig_L1emIso_eta,   "trig_L1emIso_eta[4]/D");
  mytree_->Branch("trig_L1emIso_phi",   &_trig_L1emIso_phi,   "trig_L1emIso_phi[4]/D");
  mytree_->Branch("trig_L1emIso_energy",&_trig_L1emIso_energy,"trig_L1emIso_energy[4]/D");
  mytree_->Branch("trig_L1emIso_et",    &_trig_L1emIso_et,    "trig_L1emIso_et[4]/D");
  //
  mytree_->Branch("trig_L1emNonIso_N",     &_trig_L1emNonIso_N,     "trig_L1emNonIso_N/I");
  mytree_->Branch("trig_L1emNonIso_ieta",  &_trig_L1emNonIso_ieta,  "trig_L1emNonIso_ieta[4]/I");
  mytree_->Branch("trig_L1emNonIso_iphi",  &_trig_L1emNonIso_iphi,  "trig_L1emNonIso_iphi[4]/I");
  mytree_->Branch("trig_L1emNonIso_rank",  &_trig_L1emNonIso_rank,  "trig_L1emNonIso_rank[4]/I");
  mytree_->Branch("trig_L1emNonIso_eta",   &_trig_L1emNonIso_eta,   "trig_L1emNonIso_eta[4]/D");
  mytree_->Branch("trig_L1emNonIso_phi",   &_trig_L1emNonIso_phi,   "trig_L1emNonIso_phi[4]/D");
  mytree_->Branch("trig_L1emNonIso_energy",&_trig_L1emNonIso_energy,"trig_L1emNonIso_energy[4]/D");
  mytree_->Branch("trig_L1emNonIso_et",    &_trig_L1emNonIso_et,    "trig_L1emNonIso_et[4]/D");
  
  // L1 candidates : modified collection
  mytree_->Branch("trig_L1emIso_N_M",     &_trig_L1emIso_N_M,     "trig_L1emIso_N_M/I");
  mytree_->Branch("trig_L1emIso_ieta_M",  &_trig_L1emIso_ieta_M,  "trig_L1emIso_ieta_M[4]/I");
  mytree_->Branch("trig_L1emIso_iphi_M",  &_trig_L1emIso_iphi_M,  "trig_L1emIso_iphi_M[4]/I");
  mytree_->Branch("trig_L1emIso_rank_M",  &_trig_L1emIso_rank_M,  "trig_L1emIso_rank_M[4]/I");
  mytree_->Branch("trig_L1emIso_eta_M",   &_trig_L1emIso_eta_M,   "trig_L1emIso_eta_M[4]/D");
  mytree_->Branch("trig_L1emIso_phi_M",   &_trig_L1emIso_phi_M,   "trig_L1emIso_phi_M[4]/D");
  mytree_->Branch("trig_L1emIso_energy_M",&_trig_L1emIso_energy_M,"trig_L1emIso_energy_M[4]/D");
  mytree_->Branch("trig_L1emIso_et_M",    &_trig_L1emIso_et_M,    "trig_L1emIso_et_M[4]/D");
  //
  mytree_->Branch("trig_L1emNonIso_N_M",     &_trig_L1emNonIso_N_M,     "trig_L1emNonIso_N/I");
  mytree_->Branch("trig_L1emNonIso_ieta_M",  &_trig_L1emNonIso_ieta_M,  "trig_L1emNonIso_ieta_M[4]/I");
  mytree_->Branch("trig_L1emNonIso_iphi_M",  &_trig_L1emNonIso_iphi_M,  "trig_L1emNonIso_iphi_M[4]/I");
  mytree_->Branch("trig_L1emNonIso_rank_M",  &_trig_L1emNonIso_rank_M,  "trig_L1emNonIso_rank_M[4]/I");
  mytree_->Branch("trig_L1emNonIso_eta_M",   &_trig_L1emNonIso_eta_M,   "trig_L1emNonIso_eta_M[4]/D");
  mytree_->Branch("trig_L1emNonIso_phi_M",   &_trig_L1emNonIso_phi_M,   "trig_L1emNonIso_phi_M[4]/D");
  mytree_->Branch("trig_L1emNonIso_energy_M",&_trig_L1emNonIso_energy_M,"trig_L1emNonIso_energy_M[4]/D");
  mytree_->Branch("trig_L1emNonIso_et_M",    &_trig_L1emNonIso_et_M,    "trig_L1emNonIso_et_M[4]/D");

  // pre/post - firing
  mytree_->Branch("trig_preL1emIso_N",     &_trig_preL1emIso_N,     "trig_preL1emIso_N/I");
  mytree_->Branch("trig_preL1emIso_ieta",  &_trig_preL1emIso_ieta,  "trig_preL1emIso_ieta[4]/I");
  mytree_->Branch("trig_preL1emIso_iphi",  &_trig_preL1emIso_iphi,  "trig_preL1emIso_iphi[4]/I");
  mytree_->Branch("trig_preL1emIso_rank",  &_trig_preL1emIso_rank,  "trig_preL1emIso_rank[4]/I");
  //
  mytree_->Branch("trig_preL1emNonIso_N",     &_trig_preL1emNonIso_N,     "trig_preL1emNonIso_N/I");
  mytree_->Branch("trig_preL1emNonIso_ieta",  &_trig_preL1emNonIso_ieta,  "trig_preL1emNonIso_ieta[4]/I");
  mytree_->Branch("trig_preL1emNonIso_iphi",  &_trig_preL1emNonIso_iphi,  "trig_preL1emNonIso_iphi[4]/I");
  mytree_->Branch("trig_preL1emNonIso_rank",  &_trig_preL1emNonIso_rank,  "trig_preL1emNonIso_rank[4]/I");
  //
  mytree_->Branch("trig_postL1emIso_N",     &_trig_postL1emIso_N,     "trig_postL1emIso_N/I");
  mytree_->Branch("trig_postL1emIso_ieta",  &_trig_postL1emIso_ieta,  "trig_postL1emIso_ieta[4]/I");
  mytree_->Branch("trig_postL1emIso_iphi",  &_trig_postL1emIso_iphi,  "trig_postL1emIso_iphi[4]/I");
  mytree_->Branch("trig_postL1emIso_rank",  &_trig_postL1emIso_rank,  "trig_postL1emIso_rank[4]/I");
  //
  mytree_->Branch("trig_postL1emNonIso_N",     &_trig_postL1emNonIso_N,     "trig_postL1emNonIso_N/I");
  mytree_->Branch("trig_postL1emNonIso_ieta",  &_trig_postL1emNonIso_ieta,  "trig_postL1emNonIso_ieta[4]/I");
  mytree_->Branch("trig_postL1emNonIso_iphi",  &_trig_postL1emNonIso_iphi,  "trig_postL1emNonIso_iphi[4]/I");
  mytree_->Branch("trig_postL1emNonIso_rank",  &_trig_postL1emNonIso_rank,  "trig_postL1emNonIso_rank[4]/I");
  //
  mytree_->Branch("trig_nMaskedRCT",      &_trig_nMaskedRCT,     "trig_nMaskedRCT/I");      
  mytree_->Branch("trig_iMaskedRCTeta",   &_trig_iMaskedRCTeta,  "trig_iMaskedRCTeta[trig_nMaskedRCT]/I");                                          
  mytree_->Branch("trig_iMaskedRCTcrate", &_trig_iMaskedRCTcrate,"trig_iMaskedRCTcrate[trig_nMaskedRCT]/I");                                                    
  mytree_->Branch("trig_iMaskedRCTphi",   &_trig_iMaskedRCTphi,  "trig_iMaskedRCTphi[trig_nMaskedRCT]/I");
  mytree_->Branch("trig_nMaskedCh",       &_trig_nMaskedCh,      "trig_nMaskedCh/I");    
  mytree_->Branch("trig_iMaskedTTeta",    &_trig_iMaskedTTeta,   "trig_iMaskedTTeta[trig_nMaskedCh]/I");   
  mytree_->Branch("trig_iMaskedTTphi",    &_trig_iMaskedTTphi,   "trig_iMaskedTTphi[trig_nMaskedCh]/I");      	

  // Electrons
  // ele 4V
  m_electrons = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);
  //
  mytree_->Branch("ele_N",&ele_N,"ele_N/I");
  //
  mytree_->Branch("ele_echarge",ele_echarge,"ele_echarge[10]/I");
  mytree_->Branch("ele_he",ele_he,"ele_he[10]/D");
  //
  mytree_->Branch("ele_pin_mode",ele_pin_mode,"ele_pin_mode[10]/D");
  mytree_->Branch("ele_pout_mode",ele_pout_mode,"ele_pout_mode[10]/D");
  mytree_->Branch("ele_pin_mean",ele_pin_mean,"ele_pin_mean[10]/D");
  mytree_->Branch("ele_pout_mean",ele_pout_mean,"ele_pout_mean[10]/D");
  mytree_->Branch("ele_pTin_mode",ele_pTin_mode,"ele_pTin_mode[10]/D");
  mytree_->Branch("ele_pTout_mode",ele_pTout_mode,"ele_pTout_mode[10]/D");
  mytree_->Branch("ele_pTin_mean",ele_pTin_mean,"ele_pTin_mean[10]/D");
  mytree_->Branch("ele_pTout_mean",ele_pTout_mean,"ele_pTout_mean[10]/D");
  //
  mytree_->Branch("ele_eseedpout",ele_eseedpout,"ele_eseedpout[10]/D");
  mytree_->Branch("ele_ep",ele_ep,"ele_ep[10]/D");
  mytree_->Branch("ele_eseedp",ele_eseedp,"ele_eseedp[10]/D");
  mytree_->Branch("ele_eelepout",ele_eelepout,"ele_eelepout[10]/D");
  mytree_->Branch("ele_deltaetaseed",ele_deltaetaseed,"ele_deltaetaseed[10]/D");
  mytree_->Branch("ele_deltaphiseed",ele_deltaphiseed,"ele_deltaphiseed[10]/D");
  mytree_->Branch("ele_deltaetaele",ele_deltaetaele,"ele_deltaetaele[10]/D");
  mytree_->Branch("ele_deltaphiele",ele_deltaphiele,"ele_deltaphiele[10]/D");
  mytree_->Branch("ele_deltaetain",ele_deltaetain,"ele_deltaetain[10]/D");
  mytree_->Branch("ele_deltaphiin",ele_deltaphiin,"ele_deltaphiin[10]/D");
  //
  mytree_->Branch("ele_calo_energy",ele_calo_energy,"ele_calo_energy[10]/D");
  mytree_->Branch("ele_sclE",ele_sclE,"ele_sclE[10]/D");
  mytree_->Branch("ele_sclEt",ele_sclEt,"ele_sclEt[10]/D");
  mytree_->Branch("ele_sclEta",ele_sclEta,"ele_sclEta[10]/D");
  mytree_->Branch("ele_sclPhi",ele_sclPhi,"ele_sclPhi[10]/D");
  //
  mytree_->Branch("ele_sigmaietaieta",ele_sigmaietaieta,"ele_sigmaietaieta[10]/D");
  mytree_->Branch("ele_sigmaetaeta",ele_sigmaetaeta,"ele_sigmaetaeta[10]/D");
  mytree_->Branch("ele_fbrem",ele_fbrem,"ele_fbrem[10]/D");
  mytree_->Branch("ele_ECAL_fbrem",&ele_ECAL_fbrem,"ele_ECAL_fbrem[10]/D");
  mytree_->Branch("ele_mva",ele_mva,"ele_mva[10]/D");
  mytree_->Branch("ele_isbarrel",ele_isbarrel,"ele_isbarrel[10]/I");
  mytree_->Branch("ele_isendcap",ele_isendcap,"ele_isendcap[10]/I");
  mytree_->Branch("ele_isEBetaGap",ele_isEBetaGap,"ele_isEBetaGap[10]/I");
  mytree_->Branch("ele_isEBphiGap",ele_isEBphiGap,"ele_isEBphiGap[10]/I");
  mytree_->Branch("ele_isEEdeeGap",ele_isEEdeeGap,"ele_isEEdeeGap[10]/I");
  mytree_->Branch("ele_isEEringGap",ele_isEEringGap,"ele_isEEringGap[10]/I");
  mytree_->Branch("ele_isecalDriven",ele_isecalDriven,"ele_isecalDriven[10]/I");
  mytree_->Branch("ele_istrackerDriven",ele_istrackerDriven,"ele_istrackerDriven[10]/I");
  mytree_->Branch("ele_eClass",ele_eClass,"ele_eClass[10]/I");
  mytree_->Branch("ele_missing_hits",ele_missing_hits,"ele_missing_hits[10]/I");
  mytree_->Branch("ele_lost_hits",ele_lost_hits,"ele_lost_hits[10]/I");
  mytree_->Branch("ele_chi2_hits",ele_chi2_hits,"ele_chi2_hits[10]/D");	
  //
  mytree_->Branch("ele_severityLevelSeed",ele_severityLevelSeed,"ele_severityLevelSeed[10]/I");
  mytree_->Branch("ele_severityLevelClusters",ele_severityLevelClusters,"ele_severityLevelClusters[10]/I");
	
  //Conversion Removal
  mytree_->Branch("ele_isConversion",ele_isConversion,"ele_isConversion[10]/I");
  mytree_->Branch("ele_convFound",ele_convFound,"ele_convFound[10]/I");
  mytree_->Branch("ele_conv_dist",&ele_conv_dist,"ele_conv_dist[10]/D");
  mytree_->Branch("ele_conv_dcot",&ele_conv_dcot,"ele_conv_dcot[10]/D");
  //
  mytree_->Branch("ele_track_x",ele_track_x,"ele_track_x[10]/D");
  mytree_->Branch("ele_track_y",ele_track_y,"ele_track_y[10]/D");
  mytree_->Branch("ele_track_z",ele_track_z,"ele_track_z[10]/D");	
  //
  mytree_->Branch("ele_vertex_x",ele_vertex_x,"ele_vertex_x[10]/D");
  mytree_->Branch("ele_vertex_y",ele_vertex_y,"ele_vertex_y[10]/D");
  mytree_->Branch("ele_vertex_z",ele_vertex_z,"ele_vertex_z[10]/D"); 
  //
  mytree_->Branch("ele_tkSumPt_dr03",ele_tkSumPt_dr03,"ele_tkSumPt_dr03[10]/D"); 
  mytree_->Branch("ele_ecalRecHitSumEt_dr03",ele_ecalRecHitSumEt_dr03,"ele_ecalRecHitSumEt_dr03[10]/D"); 
  mytree_->Branch("ele_hcalDepth1TowerSumEt_dr03",ele_hcalDepth1TowerSumEt_dr03,"ele_hcalDepth1TowerSumEt_dr03[10]/D"); 
  mytree_->Branch("ele_hcalDepth2TowerSumEt_dr03",ele_hcalDepth2TowerSumEt_dr03,"ele_hcalDepth2TowerSumEt_dr03[10]/D"); 
  //
  mytree_->Branch("ele_tkSumPt_dr04",ele_tkSumPt_dr04,"ele_tkSumPt_dr04[10]/D"); 
  mytree_->Branch("ele_ecalRecHitSumEt_dr04",ele_ecalRecHitSumEt_dr04,"ele_ecalRecHitSumEt_dr04[10]/D"); 
  mytree_->Branch("ele_hcalDepth1TowerSumEt_dr04",ele_hcalDepth1TowerSumEt_dr04,"ele_hcalDepth1TowerSumEt_dr04[10]/D"); 
  mytree_->Branch("ele_hcalDepth2TowerSumEt_dr04",ele_hcalDepth2TowerSumEt_dr04,"ele_hcalDepth2TowerSumEt_dr04[10]/D"); 
  //
  mytree_->Branch("ele_expected_inner_hits",ele_expected_inner_hits,"ele_expected_inner_hits[10]/I");
  mytree_->Branch("ele_sclNclus",ele_sclNclus,"ele_sclNclus[10]/I");
  //
  // L1 matching
  mytree_->Branch("ele_RCTeta",         &_ele_RCTeta,          "ele_RCTeta[10]/I");
  mytree_->Branch("ele_RCTphi",         &_ele_RCTphi,          "ele_RCTphi[10]/I");
  //
  mytree_->Branch("ele_RCTL1iso",       &_ele_RCTL1iso,        "ele_RCTL1iso[10]/I");
  mytree_->Branch("ele_RCTL1noniso",    &_ele_RCTL1noniso,     "ele_RCTL1noniso[10]/I");
  mytree_->Branch("ele_RCTL1iso_M",       &_ele_RCTL1iso_M,        "ele_RCTL1iso_M[10]/I");
  mytree_->Branch("ele_RCTL1noniso_M",    &_ele_RCTL1noniso_M,     "ele_RCTL1noniso_M[10]/I");
  //
  mytree_->Branch("ele_TTetaVect",      &_ele_TTetaVect,       "ele_TTetaVect[10][50]/I");
  mytree_->Branch("ele_TTphiVect",      &_ele_TTphiVect,       "ele_TTphiVect[10][50]/I");
  mytree_->Branch("ele_TTetVect",       &_ele_TTetVect,        "ele_TTetVect[10][50]/D");
  //
  mytree_->Branch("ele_RCTetaVect",     &_ele_RCTetaVect,      "ele_RCTetaVect[10][10]/I");
  mytree_->Branch("ele_RCTphiVect",     &_ele_RCTphiVect,      "ele_RCTphiVect[10][10]/I");
  mytree_->Branch("ele_RCTetVect",      &_ele_RCTetVect,       "ele_RCTetVect[10][10]/D");
  mytree_->Branch("ele_RCTL1isoVect",   &_ele_RCTL1isoVect,    "ele_RCTL1isoVect[10][10]/I");
  mytree_->Branch("ele_RCTL1nonisoVect",&_ele_RCTL1nonisoVect, "ele_RCTL1nonisoVect[10][10]/I");
  mytree_->Branch("ele_RCTL1isoVect_M",   &_ele_RCTL1isoVect_M,    "ele_RCTL1isoVect_M[10][10]/I");
  mytree_->Branch("ele_RCTL1nonisoVect_M",&_ele_RCTL1nonisoVect_M, "ele_RCTL1nonisoVect_M[10][10]/I");
  
  // SuperClusters
  mytree_->Branch("sc_N",   &_sc_N,   "sc_N[2]/I");
  mytree_->Branch("sc_E",   &_sc_E,   "sc_E[50]/D");
  mytree_->Branch("sc_Et",  &_sc_Et,  "sc_Et[50]/D");
  mytree_->Branch("sc_Eta", &_sc_Eta, "sc_Eta[50]/D");
  mytree_->Branch("sc_Phi", &_sc_Phi, "sc_Phi[50]/D");
  mytree_->Branch("sc_severityLevelSeed", &_sc_severityLevelSeed, "sc_severityLevelSeed[50]/I");
  mytree_->Branch("sc_he",&_sc_he, "sc_he[50]/D");
  mytree_->Branch("sc_sigmaietaieta",&_sc_sigmaietaieta, "sc_sigmaietaieta[50]/D");
  mytree_->Branch("sc_hcalDepth1TowerSumEt_dr03", &_sc_hcalDepth1TowerSumEt_dr03, "sc_hcalDepth1TowerSumEt_dr03[50]/D");
  mytree_->Branch("sc_hcalDepth2TowerSumEt_dr03", &_sc_hcalDepth2TowerSumEt_dr03, "sc_hcalDepth2TowerSumEt_dr03[50]/D");
  mytree_->Branch("sc_ecalRecHitSumEt_dr03",      &_sc_ecalRecHitSumEt_dr03,      "sc_ecalRecHitSumEt_dr03[50]/D");
  mytree_->Branch("sc_trkiso_dr03",               &_sc_trkiso_dr03,               "sc_trkiso_dr03[50]/D");
  // SC Level-1 information
  mytree_->Branch("sc_TTetaVect", &_sc_TTetaVect, "_sc_TTetaVect[50][50]/I");
  mytree_->Branch("sc_TTphiVect", &_sc_TTphiVect, "_sc_TTphiVect[50][50]/I");
  mytree_->Branch("sc_TTetVect", &_sc_TTetVect, "_sc_TTetVect[50][50]/D");
  mytree_->Branch("sc_RCTetaVect", &_sc_RCTetaVect, "_sc_RCTetaVect[50][10]/I");
  mytree_->Branch("sc_RCTphiVect", &_sc_RCTphiVect, "_sc_RCTphiVect[50][10]/I");
  mytree_->Branch("sc_RCTetVect", &_sc_RCTetVect, "_sc_RCTetVect[50][10]/D");
  mytree_->Branch("sc_RCTL1isoVect", &_sc_RCTL1isoVect, "_sc_RCTL1isoVect[50][10]/I");
  mytree_->Branch("sc_RCTL1nonisoVect", &_sc_RCTL1nonisoVect, "_sc_RCTL1nonisoVect[50][10]/I");
  mytree_->Branch("sc_RCTL1isoVect_M", &_sc_RCTL1isoVect_M, "_sc_RCTL1isoVect_M[50][10]/I");
  mytree_->Branch("sc_RCTL1nonisoVect_M", &_sc_RCTL1nonisoVect_M, "_sc_RCTL1nonisoVect_M[50][10]/I");
  //
  mytree_->Branch("sc_RCTeta", &_sc_RCTeta, "_sc_RCTeta[50]/I");
  mytree_->Branch("sc_RCTphi", &_sc_RCTphi, "_sc_RCTphi[50]/I");
  mytree_->Branch("sc_RCTL1iso", &_sc_RCTL1iso, "_sc_RCTL1iso[50]/I");
  mytree_->Branch("sc_RCTL1noniso", &_sc_RCTL1noniso, "_sc_RCTL1noniso[50]/I");
  mytree_->Branch("sc_RCTL1iso_M", &_sc_RCTL1iso_M, "_sc_RCTL1iso_M[50]/I");
  mytree_->Branch("sc_RCTL1noniso_M", &_sc_RCTL1noniso_M, "_sc_RCTL1noniso_M[50]/I");
	
}

// ====================================================================================
ElectronL1Study::~ElectronL1Study()
// ====================================================================================
{
  delete m_electrons ;
}

// ====================================================================================
void ElectronL1Study::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // Geometry and magnetic field
  bool updateField(false);
  if (cacheIDMagField_!=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier()){
    updateField = true;
    cacheIDMagField_=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier();
    iSetup.get<IdealMagneticFieldRecord>().get(theMagField_);
  }
	
  bool updateGeometry(false);
  if (cacheIDTDGeom_!=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier()){
    updateGeometry = true;
    cacheIDTDGeom_=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier();
    iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
  }
	
  if(updateField || updateGeometry){
    mtsTransform_ = new MultiTrajectoryStateTransform(trackerHandle_.product(),theMagField_.product());
  }

  unsigned long long cacheIDGeom = 0;
  if(cacheIDGeom!=iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
    cacheIDGeom = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
    iSetup.get<CaloGeometryRecord>().get(theCaloGeom_);
  }

  edm::ESHandle<CaloSubdetectorGeometry> theEndcapGeometry_handle, theBarrelGeometry_handle;
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap",theEndcapGeometry_handle);
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
  iSetup.get<IdealGeometryRecord>().get(eTTmap_);
  theEndcapGeometry_ = &(*theEndcapGeometry_handle);
  theBarrelGeometry_ = &(*theBarrelGeometry_handle);

  // Topology (for sigmaIetaIeta)
  unsigned long long cacheIDTopo_=0;
  edm::ESHandle<CaloTopology> theCaloTopo;
  if (cacheIDTopo_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
    cacheIDTopo_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
    iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  }
  topology_ = theCaloTopo.product() ;
	
  // Tree Maker
  if(PrintDebug_) std::cout << "Init()" << std::endl;
  Init();
  //
  if(PrintDebug_) std::cout << "FillEvent (iEvent, iSetup);" << std::endl;
  FillEvent (iEvent, iSetup);
  //
  if(DoFillTrigger_) FillTrigger (iEvent, iSetup);
  //
  if(PrintDebug_)  std::cout << "m_electrons -> Clear() ;" << std::endl;
  m_electrons -> Clear() ;
  if(DoFillEle_) FillEle (iEvent, iSetup);
  //
  if(PrintDebug_) std::cout << "FillSuperClusters(iEvent, iSetup);" << std::endl;
  if(DoFillSC_) FillSuperClusters(iEvent, iSetup);
  //
  if(DoFillSpikes_) FillSpikes(iEvent, iSetup);
  //
  if(PrintDebug_) std::cout << "mytree_->Fill();" << std::endl;
  mytree_->Fill();
  //
	
} // analyze

// ====================================================================================
void ElectronL1Study::FillEvent (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  nEvent = iEvent.id().event();
  nRun   = iEvent.id().run();
  nLumi  = iEvent.luminosityBlock();

  if(PrintDebug_) cout << "nRun=" << nRun << " | nEvent=" << nEvent << endl;

  // -----------------
  // Vertices
  // -----------------
  Handle<reco::VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(VerticesTag_,recoPrimaryVertexCollection);
	
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByType(recoBeamSpotHandle);
  const reco::BeamSpot bs = *recoBeamSpotHandle;
	
  int vtx_counter=0;
  _vtx_N = recoPrimaryVertexCollection->size();
	
  // select the primary vertex as the one with higest sum of (pt)^2 of tracks                                                                               
  PrimaryVertexSorter PVSorter;
  std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(recoPrimaryVertexCollection.product()) );
	
  if(_vtx_N > 0) {
    GlobalPoint local_vertexPosition(sortedVertices.front().position().x(),
				     sortedVertices.front().position().y(),
				     sortedVertices.front().position().z());
    vertexPosition = local_vertexPosition;
  }
  else {
    GlobalPoint local_vertexPosition(bs.position().x(),
				     bs.position().y(),
				     bs.position().z());
    vertexPosition = local_vertexPosition;
  }
  for( std::vector<reco::Vertex>::const_iterator PV = sortedVertices.begin(); PV != sortedVertices.end(); ++PV){
    if(vtx_counter > 199 ) continue;
		
    _vtx_normalizedChi2[vtx_counter] = PV->normalizedChi2();
    _vtx_ndof[vtx_counter] = PV->ndof();
    _vtx_nTracks[vtx_counter] = PV->tracksSize();
    _vtx_d0[vtx_counter] = PV->position().Rho();
    _vtx_x[vtx_counter] = PV->x();
    _vtx_y[vtx_counter] = PV->y();
    _vtx_z[vtx_counter] = PV->z();
		
    vtx_counter++;
  } // for loop on primary vertices
	
  if(vtx_counter>199) { _vtx_N = 200; cout << "Number of primary vertices>199, vtx_N set to 200" << endl;}
	
}

// ====================================================================================
void ElectronL1Study::FillTrigger (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  /////////
  // HLT //
  /////////

  // Get triggered HLT paths //
  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByLabel (HLTTag_,triggerResultsHandle);
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle);

  const int nTrig = 4;
  string nad_hlt_name[nTrig] = {"HLT_Activity_Ecal_SC","HLT_L1SingleEG5","HLT_L1SingleEG12",
				"HLT_DoubleEle"};

  bool nad_hlt_trig[nTrig];
  for(int iTrig=0 ; iTrig<nTrig ; iTrig++) {
    trig_HLT_path[iTrig] = 0;
    nad_hlt_trig[iTrig] = false;
  }

  if(PrintDebug_HLT_) {
    cout << "triggerResultsHandle->size() = " << triggerResultsHandle->size() << endl
	 << "<----- TRIGGERED HLT PATHS ----->" << endl;
  }

  // loop over HLT results
  strcpy(trig_fired_names,"*");
  for (int iHLT=0; iHLT<static_cast<int>(triggerResultsHandle->size()); iHLT++) {	
    if (triggerResultsHandle->accept (iHLT)) {

      if( PrintDebug_HLT_ )
	cout << " -- " << string(triggerNames.triggerName(iHLT)) << endl;

      // loop over my paths
      for(int iTrig=0 ; iTrig<nTrig ; iTrig++) {
	if( string( triggerNames.triggerName(iHLT) ).find( nad_hlt_name[iTrig] ) != string::npos) {
	  nad_hlt_trig[iTrig] = true;
	  trig_HLT_path[iTrig] = 1;
	}
      }

      // write name of triggered paths
      trig_hltInfo[iHLT] = 1;
      if ( strlen(trig_fired_names) <= 4950) {
	const char* c_str();
	string hlt_string = triggerNames.triggerName(iHLT);
	strcat(trig_fired_names,hlt_string.c_str());
	strcat(trig_fired_names,"*");
      }
      else trig_hltInfo[iHLT] = 0;

    } // endif HLT accept
  } // end loop over results
	  
  // Cross-check
  /*
  for(int iTrig=0 ; iTrig<nTrig ; iTrig++)
    if( nad_hlt_trig[iTrig] )
      trig_HLT_path[iTrig] = 1;
  */
  if(PrintDebug_HLT_) {
    cout << "Cross-check HLT branches" << endl;
    for(int iTrig=0 ; iTrig<nTrig ; iTrig++)
      cout << "HLT name : " << nad_hlt_name[iTrig] << " : " << trig_HLT_path[iTrig] << endl;
  }

  // Paths from list given in .py file //
  UInt_t trigger_size = triggerResultsHandle->size();
  int passEleTrigger  = 0;
	
  // Electron Triggers
  for(int ipath=0;ipath< (int) HLT_ElePaths_.size();ipath++) {
    //cout << " i = " << ipath << " trigger = " << HLT_Paths_[ipath] << endl;
    UInt_t trigger_position = triggerNames.triggerIndex(HLT_ElePaths_[ipath]); //hltpath_);
    if (trigger_position < trigger_size) passEleTrigger = (int)triggerResultsHandle->accept(trigger_position);
    if (passEleTrigger==1) _trig_isEleHLTpath = 1;
  } // for loop on HLT Elepaths


  ///////////////////
  // L1 CANDIDATES //
  ///////////////////
  
  edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl ;
  edm::Handle< l1extra::L1EmParticleCollection > emIsolColl ;
  edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl_M ;
  edm::Handle< l1extra::L1EmParticleCollection > emIsolColl_M ;  

  int counter = 0;

  if( GetL1_ ) {
    if( !GetL1M_ ) {      
      // standard collection ALONE
      iEvent.getByLabel("l1extraParticles","NonIsolated", emNonisolColl ) ;
      iEvent.getByLabel("l1extraParticles","Isolated", emIsolColl ) ;
    } else {
      // standard collection
      iEvent.getByLabel("l1extraParticlesOnline","NonIsolated", emNonisolColl ) ;
      iEvent.getByLabel("l1extraParticlesOnline","Isolated", emIsolColl ) ;
      // modified collection
      iEvent.getByLabel("l1extraParticles","NonIsolated", emNonisolColl_M ) ;
      iEvent.getByLabel("l1extraParticles","Isolated", emIsolColl_M ) ;
    }
  }
  
  // DATA L1 CANDIDATES //
  if( GetL1_ ) {  
    // Isolated candidates
    _trig_L1emIso_N = emIsolColl->size();
    if(PrintDebug_) cout << "N L1 candidate iso : " << _trig_L1emIso_N << endl;
    counter = 0;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl->begin(); emItr != emIsolColl->end() ;++emItr) {
      _trig_L1emIso_eta[counter]    = emItr->eta();
      _trig_L1emIso_phi[counter]    = emItr->phi();
      _trig_L1emIso_energy[counter] = emItr->energy();
      _trig_L1emIso_et[counter]     = emItr->et();
      if(!aod_) {
	_trig_L1emIso_ieta[counter] = emItr->gctEmCand()->regionId().ieta();
	_trig_L1emIso_iphi[counter] = emItr->gctEmCand()->regionId().iphi();
	_trig_L1emIso_rank[counter] = emItr->gctEmCand()->rank(); // 1 rank adc = 1 GeV (in 2011)
      }
      else {
	_trig_L1emIso_ieta[counter] = convertEtaToRCT( emItr->eta() );
	_trig_L1emIso_iphi[counter] = convertPhiToRCT( emItr->phi() );
	_trig_L1emIso_rank[counter] = (int)( emItr->et() );
      }
      counter++;
    }
	  
    // Non Isolated candidates
    _trig_L1emNonIso_N = emNonisolColl->size();
    if(PrintDebug_) cout << "N L1 candidate noniso : " << _trig_L1emNonIso_N << endl;
    counter = 0;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl->begin(); emItr != emNonisolColl->end() ;++emItr) {
      _trig_L1emNonIso_eta[counter]    = emItr->eta();
      _trig_L1emNonIso_phi[counter]    = emItr->phi();
      _trig_L1emNonIso_energy[counter] = emItr->energy();
      _trig_L1emNonIso_et[counter]     = emItr->et();
      if(!aod_) {
	_trig_L1emNonIso_ieta[counter] = emItr->gctEmCand()->regionId().ieta();
	_trig_L1emNonIso_iphi[counter] = emItr->gctEmCand()->regionId().iphi();
	_trig_L1emNonIso_rank[counter] = emItr->gctEmCand()->rank(); // 1 rank adc = 1 GeV (in 2011)
      }
      else {
	_trig_L1emNonIso_ieta[counter] = convertEtaToRCT( emItr->eta() );
	_trig_L1emNonIso_iphi[counter] = convertPhiToRCT( emItr->phi() );
	_trig_L1emNonIso_rank[counter] = (int)( emItr->et() );
      }
      counter++;
    } // for loop on Non Iso cand
  } // endif GetL1_

  // EMULATED L1 CANDIDATES //
  if( GetL1M_ ) {
    // Isolated candidates
    _trig_L1emIso_N_M = emIsolColl_M->size();
    if(PrintDebug_) cout << "N L1 candidate noniso : " << _trig_L1emIso_N_M << endl;
    counter = 0;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl_M->begin(); emItr != emNonisolColl_M->end() ;++emItr) {
      _trig_L1emIso_eta_M[counter]    = emItr->eta();
      _trig_L1emIso_phi_M[counter]    = emItr->phi();
      _trig_L1emIso_energy_M[counter] = emItr->energy();
      _trig_L1emIso_et_M[counter]     = emItr->et();
      if(!aod_) {
	_trig_L1emIso_ieta_M[counter] = emItr->gctEmCand()->regionId().ieta();
	_trig_L1emIso_iphi_M[counter] = emItr->gctEmCand()->regionId().iphi();
	_trig_L1emIso_rank_M[counter] = emItr->gctEmCand()->rank(); // 1 rank adc = 1 GeV (in 2011)
      }
      else {
	_trig_L1emIso_ieta_M[counter] = convertEtaToRCT( emItr->eta() );
	_trig_L1emIso_iphi_M[counter] = convertPhiToRCT( emItr->phi() );
	_trig_L1emIso_rank_M[counter] = (int)( emItr->et() );
      }
      counter++;
    } // loop over isolated candidates

    // Non Isolated candidates
    _trig_L1emNonIso_N_M = emNonisolColl_M->size();
    if(PrintDebug_) cout << "N L1 candidate noniso : " << _trig_L1emNonIso_N_M << endl;
    counter = 0;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl_M->begin(); emItr != emNonisolColl_M->end() ;++emItr) {
      _trig_L1emNonIso_eta_M[counter]    = emItr->eta();
      _trig_L1emNonIso_phi_M[counter]    = emItr->phi();
      _trig_L1emNonIso_energy_M[counter] = emItr->energy();
      _trig_L1emNonIso_et_M[counter]     = emItr->et();
      if(!aod_) {
	_trig_L1emNonIso_ieta_M[counter] = emItr->gctEmCand()->regionId().ieta();
	_trig_L1emNonIso_iphi_M[counter] = emItr->gctEmCand()->regionId().iphi();
	_trig_L1emNonIso_rank_M[counter] = emItr->gctEmCand()->rank(); // 1 rank adc = 1 GeV (in 2011)
      }
      else {
	_trig_L1emNonIso_ieta_M[counter] = convertEtaToRCT( emItr->eta() );
	_trig_L1emNonIso_iphi_M[counter] = convertPhiToRCT( emItr->phi() );
	_trig_L1emNonIso_rank_M[counter] = (int)( emItr->et() );
      }
      counter++;
    } // for loop on Non Iso cand

  } // endif GetL1M_


  ///////////////////////
  // PRE/POST - FIRING //
  ///////////////////////
  if( !aod_ ) {

    if(PrintDebug_) cout << "Pre-Post firing" << endl;
    edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel( edm::InputTag(gtRecordCollectionTag_), gtRecord);
    
    // PRE-FIRING //
    const L1GtPsbWord psb = gtRecord->gtPsbWord(0xbb0d, -1);
    std::vector<int> psbel;
    std::vector<int>::const_iterator ipsbel;

    int rank, iEta, sign, regionEtaRec;
    rank = iEta = sign = regionEtaRec = -999;

    // Non-Iso
    psbel.push_back(psb.aData(4));
    psbel.push_back(psb.aData(5));
    psbel.push_back(psb.bData(4));
    psbel.push_back(psb.bData(5));
    
    counter = 0;
    for(ipsbel=psbel.begin(); ipsbel!=psbel.end(); ipsbel++) {
      rank = iEta = sign = regionEtaRec = -999;
      rank = (*ipsbel)&0x3f; // 1 rank adc = 1 GeV (2011)
      if(rank>0) {
	iEta = int(((*ipsbel)>>6)&7);
	sign = ( ((*ipsbel>>9)&1) ? -1. : 1. ); 
	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (pre, non-iso)"<<std::endl;	
	_trig_preL1emNonIso_ieta[counter] = regionEtaRec;
	_trig_preL1emNonIso_iphi[counter] = int(((*ipsbel)>>10)&0x1f);
	_trig_preL1emNonIso_rank[counter] = rank;
      }
    } // loop non iso
    _trig_preL1emNonIso_N = counter;

    // Iso
    psbel.clear();
    psbel.push_back(psb.aData(6));
    psbel.push_back(psb.aData(7));
    psbel.push_back(psb.bData(6));
    psbel.push_back(psb.bData(7));
    counter = 0;
    for(ipsbel=psbel.begin(); ipsbel!=psbel.end(); ipsbel++) {
     rank = iEta = sign = regionEtaRec = -999;
     rank = (*ipsbel)&0x3f; 
      if(rank>0) {
	iEta = int(((*ipsbel)>>6)&7);
	sign = ( ((*ipsbel>>9)&1) ? -1. : 1. ); 
	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (pre, iso)"<<std::endl;
	_trig_preL1emIso_ieta[counter] = regionEtaRec;
	_trig_preL1emIso_iphi[counter] = int(((*ipsbel)>>10)&0x1f);
	_trig_preL1emIso_rank[counter] = rank;
	counter++;
      }
    }//loop Iso
    _trig_preL1emIso_N = counter;

    // POST-FIRING //
    const L1GtPsbWord psb2 = gtRecord->gtPsbWord(0xbb0d, 1);
    std::vector<int> psbel2;
    std::vector<int>::const_iterator ipsbel2;

    // Non Iso
    psbel2.push_back(psb2.aData(4));
    psbel2.push_back(psb2.aData(5));
    psbel2.push_back(psb2.bData(4));
    psbel2.push_back(psb2.bData(5));
    counter = 0;
    for(ipsbel2=psbel2.begin(); ipsbel2!=psbel2.end(); ipsbel2++) {
      rank = (*ipsbel2)&0x3f; 
      if(rank>0) {
	iEta = int(((*ipsbel2)>>6)&7);
	sign = ( ((*ipsbel2>>9)&1) ? -1. : 1. ); 
	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (post, non-iso)"<<std::endl;
	_trig_postL1emNonIso_ieta[counter] = regionEtaRec;
	_trig_postL1emNonIso_iphi[counter] = int(((*ipsbel2)>>10)&0x1f);
	_trig_postL1emNonIso_rank[counter] = rank;
	counter++;
      }
    }//loop Noniso
    _trig_postL1emNonIso_N = counter;
   
    // Iso
    psbel2.clear();
    psbel2.push_back(psb2.aData(6));
    psbel2.push_back(psb2.aData(7));
    psbel2.push_back(psb2.bData(6));
    psbel2.push_back(psb2.bData(7));
    counter = 0;
    for(ipsbel2=psbel2.begin(); ipsbel2!=psbel2.end(); ipsbel2++) {
      rank = (*ipsbel2)&0x3f; // ET in ADC count... 1 ADC count = 0.5 GeV
      if(rank>0) {
	iEta = int(((*ipsbel2)>>6)&7);
	sign = ( ((*ipsbel2>>9)&1) ? -1. : 1. ); 
       	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (post, iso)"<<std::endl;
	_trig_postL1emIso_ieta[counter] = regionEtaRec;
	_trig_postL1emIso_iphi[counter] = int(((*ipsbel2)>>10)&0x1f);
	_trig_postL1emIso_rank[counter] = rank;
	counter++;
      }
    }//loop Iso
    _trig_postL1emIso_N = counter;
 
  }


  ////////////////////////
  // TRIGGER PRIMITIVES //
  ////////////////////////

  // Data TPs //
  if( GetTP_ ) {
    ecal_tp_ = new edm::Handle<EcalTrigPrimDigiCollection> ;
    iEvent.getByLabel(tpCollectionNormal_,*ecal_tp_);
    if(PrintDebug_) cout << "got by label TPs" << endl;
    _trig_tower_N = ecal_tp_->product()->size();
    if(PrintDebug_) {
      cout << "TP Normal collection size=" << ecal_tp_->product()->size() << endl ;
      cout << "is gonna get the TP data" << endl;
    }
  
    for (int i=0 ; i<_trig_tower_N ; i++) {
      if(PrintDebug_) cout << "loop iteration #" << i << endl;
      EcalTriggerPrimitiveDigi d_ = (*(ecal_tp_->product()))[i]; // EcalTriggerPrimitiveDigi d
      if(PrintDebug_) cout << "got the trigger primitive" << endl;
      TPtowid_ = d_.id(); // const EcalTrigTowerDetId TPtowid
      if(PrintDebug_) cout << "got the tower id" << endl;
      _trig_tower_iphi[i] = TPtowid_.iphi() ;
      _trig_tower_ieta[i] = TPtowid_.ieta() ;
      if(PrintDebug_) cout << "got the ieta and iphi : " << TPtowid_.ieta() << TPtowid_.iphi() << endl;
      _trig_tower_adc[i]  = (d_[0].raw()&0xff) ;  // 0xff  <-> Et(8bits)
      if(PrintDebug_) cout << "got the adc : " << (int)(d_[0].raw()&0xff) << endl;
      _trig_tower_sFGVB[i] = d_[0].sFGVB();       // 0=spike-like / 1=EM-like
      if(PrintDebug_) cout << "got the sFGVB : " << d_[0].sFGVB() << endl;
    }
    if(PrintDebug_) cout << "finished looping" << endl;
  }

  // Zeroed-by-hand TPs //
  if( GetTP_Modif_ ) {
    ecal_tpM_ = new edm::Handle<EcalTrigPrimDigiCollection> ;
    iEvent.getByLabel(tpCollectionModif_,*ecal_tpM_);
    _trig_tower_N_M = ecal_tpM_->product()->size(); 

    for (int i=0 ; i<_trig_tower_N_M ; i++) {
      EcalTriggerPrimitiveDigi dM_ = (*(ecal_tpM_->product()))[i]; // EcalTriggerPrimitiveDigi dM
      TPtowidM_ = dM_.id(); // EcalTrigTowerDetId
      _trig_tower_iphi_M[i] = TPtowidM_.iphi() ;
      _trig_tower_ieta_M[i] = TPtowidM_.ieta() ;
      _trig_tower_adc_M[i]  = (dM_[0].raw()&0xff) ;
      _trig_tower_sFGVB_M[i] = dM_[0].sFGVB(); // 0=spike-like / 1=EM-like
    }
  }

  // Emulated TPs //
  if( GetTP_Emul_ ) {
    ecal_tpM_ = new edm::Handle<EcalTrigPrimDigiCollection> ;
    iEvent.getByLabel(tpEmulatorCollection_, *ecal_tpM_);
    if (PrintDebug_) std::cout<<"TPEmulator collection size="<< ecal_tpM_->product()->size()<<std::endl ;
    _trig_tower_N_E = ecal_tpM_->product()->size();

    for (int i=0 ; i<_trig_tower_N_E ; i++) {
      EcalTriggerPrimitiveDigi dM_ = (*(ecal_tpM_->product()))[i]; //EcalTriggerPrimitiveDigi
      TPtowidM_ = dM_.id();
      _trig_tower_iphi_E[i] = TPtowidM_.iphi() ;
      _trig_tower_ieta_E[i] = TPtowidM_.ieta() ;
    
      if(PrintDebug_)
	cout << "TTieta=" << TPtowidM_.ieta() << " TTiphi=" << TPtowidM_.iphi() << " adcEm=" ;

      for (int j=0 ; j<5 ; j++) {
	_trig_tower_adc_E[i][j] = (dM_[j].raw()&0xff) ;
	_trig_tower_sFGVB_E[i][j] = dM_[j].sFGVB(); 
	if(PrintDebug_) cout << (dM_[j].raw()&0xff) << " " ;
      }
      if(PrintDebug_)
	cout << endl;
    } 
  }


  /////////////
  // MASKING //
  /////////////

  // RCT masking //
  edm::ESHandle<L1RCTChannelMask> channelMask;
  iSetup.get<L1RCTChannelMaskRcd>().get(channelMask);
  const L1RCTChannelMask* cEs = channelMask.product();
  uint n0MaskedRCT = 0;
  for(int i = 0; i< 18; i++)
    for(int j =0; j< 2; j++)
      for(int k =0; k<28; k++)
	if(cEs->ecalMask[i][j][k]){
	  _trig_iMaskedRCTeta[n0MaskedRCT]=k;
	  _trig_iMaskedRCTphi[n0MaskedRCT]=j;
	  _trig_iMaskedRCTcrate[n0MaskedRCT]=i;
	  n0MaskedRCT++;
	}
  _trig_nMaskedRCT = n0MaskedRCT;

  // Tower masking //
  edm::ESHandle<EcalTPGTowerStatus> theEcalTPGTowerStatus_handle;
  iSetup.get<EcalTPGTowerStatusRcd>().get(theEcalTPGTowerStatus_handle);
  const EcalTPGTowerStatus * ecaltpgTowerStatus=theEcalTPGTowerStatus_handle.product();
	
  const EcalTPGTowerStatusMap &towerMap=ecaltpgTowerStatus->getMap();
  EcalTPGTowerStatusMapIterator  ittpg;
	
  uint nMaskedChannels = 0;
  for (ittpg=towerMap.begin();ittpg!=towerMap.end();++ittpg) {
		
    if ((*ittpg).second > 0)
      {
	EcalTrigTowerDetId  ttId((*ittpg).first);
	_trig_iMaskedTTeta[nMaskedChannels] = ttId.ieta();
	_trig_iMaskedTTphi[nMaskedChannels] = ttId.iphi();
	nMaskedChannels++;
      }
  }//loop trigger towers
  _trig_nMaskedCh = nMaskedChannels;  

}

// ====================================================================================
void ElectronL1Study::FillSpikes(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  
  if(PrintDebug_) cout << "starting getting the spike-like rechits " << endl;

  // for 42x	
  unsigned long long cacheSevLevel = 0;
  edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
  if(cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()){
    cacheSevLevel = iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
  }
  const EcalSeverityLevelAlgo* sl=sevLevel.product();

  // Get EB rechits
  edm::Handle<EcalRecHitCollection> rechitsEB; 
  EcalRecHitCollection::const_iterator rechitItr;

  EBDetId id;
  int flag=0;
  uint32_t sev=0;
  double thetaHit, etaHit, phiHit;

  int i=0;
  if (iEvent.getByLabel(EcalRecHitCollectionEB_, rechitsEB) ) {	
    for ( rechitItr = rechitsEB->begin(); rechitItr != rechitsEB->end(); ++rechitItr ) {	
      if(i>=5000) {
	cout << "more than 5000 spikes in run " << iEvent.id().run() << " , event " << iEvent.id().event() << endl;
	break;
      }
      id = rechitItr->id();	      
      sev = 0;
      sev = sl->severityLevel( id, *rechitsEB );
    
      if(sev < 3) continue;

      thetaHit =  (theBarrelGeometry_->getGeometry(id)->getPosition()).theta();
      etaHit =  (theBarrelGeometry_->getGeometry(id)->getPosition()).eta();
      phiHit =  (theBarrelGeometry_->getGeometry(id)->getPosition()).phi();
      const EcalTrigTowerDetId towid = id.tower();
      flag = 0;
      flag = rechitItr->recoFlag();
      if( PrintDebug_ )
	cout << "-- RecHit : "
	     << "flag=" << flag
	     << " | sev=" << sev
	     << " | thetaHit="  << thetaHit 
	     << " | phiHit=" << phiHit
	     << " | etaHit=" << etaHit
	     << " | E="      << rechitItr->energy()
	     << " | Et="     << (rechitItr->energy())*sin(thetaHit)
	     << endl;

      if(flag == EcalRecHit::kOutOfTime) spike_outOfTime[i] = 1;
      else spike_outOfTime[i] = 0;
      spike_severityLevel[i] = sev ;
      spike_time[i] = rechitItr->time();
      spike_Et[i] = (rechitItr->energy())*sin(thetaHit);
      spike_phi[i] = phiHit;
      spike_eta[i] = etaHit;
      spike_theta[i] = thetaHit;
      spike_TTiphi[i] = towid.iphi();
      spike_TTieta[i] = towid.ieta();
      spike_Riphi[i] = getGCTRegionPhi(towid.iphi());
      spike_Rieta[i] = getGCTRegionEta(towid.ieta());
    
      i++ ;
    }
  }
  if(i>=5000) spike_N = 5000;
  else spike_N = i+1;

}

// ====================================================================================
void ElectronL1Study::FillEle(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{

  // SETUP //

  // ECAL Geometry (to find RecHits' theta and tower)

  // RecHits collections
  if(PrintDebug_) cout << "<--- FillEle() --->" << endl
		       << "Get reduced RecHits" << endl;
  const EcalRecHitCollection * reducedRecHits = 0 ;
  edm::Handle< EcalRecHitCollection > reducedEBRecHits;
  edm::Handle< EcalRecHitCollection > reducedEERecHits;
  if(!aod_){
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEE"), reducedEERecHits ) ;
  }
  else{
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEE"), reducedEERecHits ) ;
  }

  // Severity
  if(PrintDebug_) cout << "Get severity handle" << endl;
  unsigned long long cacheSevLevel = 0;
  edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
  if(cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()){
    cacheSevLevel = iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
  }
  const EcalSeverityLevelAlgo* sl=sevLevel.product();

  // Conversion info
  if(PrintDebug_) cout << "Get track/dcs/B handles" << endl;
  edm::Handle<reco::TrackCollection> tracks_h;
  iEvent.getByLabel("generalTracks", tracks_h);
  //
  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel(dcsTag_, dcsHandle);
  //
  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  //
  double evt_bField, currentToBFieldScaleFactor, current;
  if(type_ == "DATA" ) {
    if ((*dcsHandle).size() != 0 ) {	
      currentToBFieldScaleFactor = 2.09237036221512717e-04;
      current = (*dcsHandle)[0].magnetCurrent();
      evt_bField = current*currentToBFieldScaleFactor;
    }
    else {	
      GlobalPoint gPoint(0.,0.,0.);
      evt_bField = magneticField->inTesla(gPoint).z();
    }
  }
  // MC
  else {
    GlobalPoint gPoint(0.,0.,0.);
    evt_bField = magneticField->inTesla(gPoint).z();
  }
  

  // HANDLE ON ELECTRONS //
  if(PrintDebug_) cout << "Get Electrons Handle" << endl;
  edm::Handle<reco::GsfElectronCollection> EleHandle ;
  iEvent.getByLabel (EleTag_.label(),EleHandle) ;


  // LOOP OVER ELECTRONS //
  if(PrintDebug_) cout << "loop over electrons" << endl;
  if(EleHandle->size() < 10 ){ ele_N = EleHandle->size(); }
  else {ele_N = 10;}
  TClonesArray &electrons = *m_electrons;
  int counter = 0;
  int nTow=0;
  int nReg=0;
  
  for(int i=0; i< ele_N; i++) {
    
    // ELECTRON PROPERTIES //
    if(PrintDebug_) cout << "- get 4-vector and ele properties" << endl;
    edm::Ref<reco::GsfElectronCollection> electronEdmRef(EleHandle,i);
    setMomentum (myvector, (*EleHandle)[i].p4());
    new (electrons[counter]) TLorentzVector (myvector);
    
    ele_echarge[counter] = (*EleHandle)[i].charge(); 
    ele_eClass[counter]   = (*EleHandle)[i].classification() ;
    ele_calo_energy[counter] = (*EleHandle)[i].caloEnergy() ;

    // Mesures pT in/out
    ele_pin_mode[counter]    = (*EleHandle)[i].trackMomentumAtVtx().R() ; 
    ele_pout_mode[counter]   = (*EleHandle)[i].trackMomentumOut().R() ; 
    ele_pTin_mode[counter]   = (*EleHandle)[i].trackMomentumAtVtx().Rho() ; 
    ele_pTout_mode[counter]  = (*EleHandle)[i].trackMomentumOut().Rho() ; 
    if(!aod_){
      ele_pin_mean[counter]    = (*EleHandle)[i].gsfTrack()->innerMomentum().R() ; 
      ele_pout_mean[counter]   = (*EleHandle)[i].gsfTrack()->outerMomentum().R(); 
      ele_pTin_mean[counter]   = (*EleHandle)[i].gsfTrack()->innerMomentum().Rho() ; 
      ele_pTout_mean[counter]  = (*EleHandle)[i].gsfTrack()->outerMomentum().Rho() ;
    }

    // Concordances E/p
    ele_eseedpout[counter] = (*EleHandle)[i].eSeedClusterOverPout();
    ele_ep[counter]        = (*EleHandle)[i].eSuperClusterOverP() ;        
    ele_eseedp[counter]    = (*EleHandle)[i].eSeedClusterOverP() ;         
    ele_eelepout[counter]  = (*EleHandle)[i].eEleClusterOverPout() ;       
 
    // Concordance géométrique track/SC
    ele_deltaetaseed[counter]  = (*EleHandle)[i].deltaEtaSeedClusterTrackAtCalo() ; 
    ele_deltaphiseed[counter]  = (*EleHandle)[i].deltaPhiSeedClusterTrackAtCalo() ;  
    ele_deltaetaele[counter]   = (*EleHandle)[i].deltaEtaEleClusterTrackAtCalo() ;  
    ele_deltaphiele[counter]   = (*EleHandle)[i].deltaPhiEleClusterTrackAtCalo() ; 
    ele_deltaetain[counter]    = (*EleHandle)[i].deltaEtaSuperClusterTrackAtVtx();
    ele_deltaphiin[counter]    = (*EleHandle)[i].deltaPhiSuperClusterTrackAtVtx();   

    // Track
    ele_track_x[counter] = (*EleHandle)[i].gsfTrack()->vx();
    ele_track_y[counter] = (*EleHandle)[i].gsfTrack()->vy();
    ele_track_z[counter] = (*EleHandle)[i].gsfTrack()->vz();
    // Vertex
    ele_vertex_x[counter] = (*EleHandle)[i].vertex().x();
    ele_vertex_y[counter] = (*EleHandle)[i].vertex().y();
    ele_vertex_z[counter] = (*EleHandle)[i].vertex().z();
    // Hits
    ele_missing_hits[counter] = (*EleHandle)[i].gsfTrack()->numberOfLostHits();
    ele_lost_hits[counter]    = (*EleHandle)[i].gsfTrack()->numberOfValidHits() ;
    ele_chi2_hits[counter]    = (*EleHandle)[i].gsfTrack()->normalizedChi2() ; 

    // Spread/quality
    ele_he[counter]      = (*EleHandle)[i].hadronicOverEm() ;
    ele_fbrem[counter] = (*EleHandle)[i].fbrem() ;
    ele_mva[counter]   = (*EleHandle)[i].mva() ;
    ele_sigmaietaieta[counter] = (*EleHandle)[i].sigmaIetaIeta() ; 
    if(!aod_)
      ele_expected_inner_hits[counter] = (*EleHandle)[i].gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    // Isolation
    ele_tkSumPt_dr03[counter]              = (*EleHandle)[i].dr03TkSumPt() ;
    ele_ecalRecHitSumEt_dr03[counter]      = (*EleHandle)[i].dr03EcalRecHitSumEt() ;
    ele_hcalDepth1TowerSumEt_dr03[counter] = (*EleHandle)[i].dr03HcalDepth1TowerSumEt() ;
    ele_hcalDepth2TowerSumEt_dr03[counter] = (*EleHandle)[i].dr03HcalDepth2TowerSumEt() ;
    ele_tkSumPt_dr04[counter]              = (*EleHandle)[i].dr04TkSumPt() ;
    ele_ecalRecHitSumEt_dr04[counter]      = (*EleHandle)[i].dr04EcalRecHitSumEt() ;
    ele_hcalDepth1TowerSumEt_dr04[counter] = (*EleHandle)[i].dr04HcalDepth1TowerSumEt() ;
    ele_hcalDepth2TowerSumEt_dr04[counter] = (*EleHandle)[i].dr04HcalDepth2TowerSumEt() ;

    // Conversion Removal
    ele_isConversion[counter] = IsConv( (*EleHandle)[i] );
    ConversionFinder convFinder;
    ConversionInfo convInfo = convFinder.getConversionInfo((*EleHandle)[i], tracks_h, evt_bField);
    ele_conv_dist[counter] = convInfo.dist();
    ele_conv_dcot[counter] = convInfo.dcot();
    ele_convFound[counter] = 0;
    if ( convFinder.isFromConversion(convInfo, 0.02, 0.02) ) ele_convFound[counter] = 1; 

    // SuperCluster
    reco::SuperClusterRef sclRef = (*EleHandle)[i].superCluster();
    math::XYZPoint sclPos = (*EleHandle)[i].superClusterPosition();
    if (!(*EleHandle)[i].ecalDrivenSeed() && (*EleHandle)[i].trackerDrivenSeed()) 
      sclRef = (*EleHandle)[i].pflowSuperCluster();
    double R=TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    double Rt=TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    //
    ele_sclNclus[counter] = sclRef->clustersSize();
    ele_sclE[counter]   = sclRef->energy() ;
    ele_sclEt[counter]  = sclRef->energy()*(Rt/R) ;
    ele_sclEta[counter] = sclRef->eta() ;
    ele_sclPhi[counter] = sclRef->phi() ;
    ele_ECAL_fbrem[counter]   = sclRef->phiWidth()/sclRef->etaWidth();
   
    // ECAL : EB/EE/gaps/driven
    if ((*EleHandle)[i].isEB()) ele_isbarrel[counter] = 1 ; 
    else  ele_isbarrel[counter] = 0 ;
    if ((*EleHandle)[i].isEE()) ele_isendcap[counter] = 1 ; 
    else  ele_isendcap[counter] = 0 ;
    if ((*EleHandle)[i].isEBEtaGap()) ele_isEBetaGap[counter] = 1 ;  
    if ((*EleHandle)[i].isEBPhiGap()) ele_isEBphiGap[counter] = 1 ;  
    if ((*EleHandle)[i].isEEDeeGap()) ele_isEEdeeGap[counter] = 1 ;  
    if ((*EleHandle)[i].isEERingGap()) ele_isEEringGap[counter] = 1 ;
    if ((*EleHandle)[i].ecalDrivenSeed()) ele_isecalDriven[counter] = 1 ;
    if ((*EleHandle)[i].trackerDrivenSeed()) ele_istrackerDriven[counter] = 1 ;

    // L1 CANDIDATE MATCHING //
    if(PrintDebug_) cout << "- match ele to L1 candidates" << endl;
    nTow=0;
    nReg=0;
    // loop over clusters in SC
    for (reco::CaloCluster_iterator clus = sclRef->clustersBegin () ;
	 clus != sclRef->clustersEnd () ; ++clus) {
      
      std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions();

      // loop over RecHits (Xtals) in cluster
      for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ;
	   detitr != clusterDetIds.end () ; ++detitr) {
	
	// check located in ECAL
	if ( (detitr -> first).det () != DetId::Ecal) {
	  std::cout << " det is " << (detitr -> first).det () << std::endl ;
	  continue ;
	}
	
	// Getting this RecHit's informations+tower
	if ((*EleHandle)[i].isEB())  
	  reducedRecHits = reducedEBRecHits.product() ; 
	else 
	  reducedRecHits = reducedEERecHits.product() ;
	//
	EcalRecHitCollection::const_iterator thishit;
	EcalRecHit myhit;
	EcalTrigTowerDetId towid;
	float thetahit;
	//
	// Barrel RecHits
	if ( (detitr -> first).subdetId () == EcalBarrel) {
	  thishit = reducedRecHits->find ( (detitr -> first) ) ;
	  if (thishit == reducedRecHits->end ()) continue;
	  myhit = (*thishit) ;
	  EBDetId detid(thishit->id());
	  towid= detid.tower();
	  thetahit =  theBarrelGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	}
	// Endcaps RecHits
	else {
	  if ( (detitr -> first).subdetId () == EcalEndcap) {
	    thishit = reducedRecHits->find ( (detitr -> first) ) ;
	    if (thishit == reducedRecHits->end ()) continue;
	    myhit = (*thishit) ;
	    EEDetId detid(thishit->id());
	    towid= (*eTTmap_).towerOf(detid);
	    thetahit =  theEndcapGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	  }
	  else continue;
	}
	
	// this RecHit's tower
	int iETA   = towid.ieta();
	int iPHI   = towid.iphi();
	int iReta  = getGCTRegionEta(iETA);
	int iRphi  = getGCTRegionPhi(iPHI);
	double iET = myhit.energy()*sin(thetahit);
	//
	bool newTow = true;
	if(nTow>0) {
	  for (int iTow=0; iTow<nTow; ++iTow) {
	    if(_ele_TTetaVect[counter][iTow] == iETA && _ele_TTphiVect[counter][iTow] == iPHI) {
	      newTow = false;
	      _ele_TTetVect[counter][iTow] +=  iET;
	    }
	  }
	}
	if(newTow) {
	  _ele_TTetaVect[counter][nTow] = iETA;
	  _ele_TTphiVect[counter][nTow] = iPHI;
	  _ele_TTetVect[counter][nTow] =  iET;
	  nTow++;
	}

	// this tower's region
	bool newReg = true;
	if(nReg>0) {
	  for (int iReg=0; iReg<nReg; ++iReg) {
	    if(_ele_RCTetaVect[counter][iReg] == iReta && _ele_RCTphiVect[counter][iReg] == iRphi) {
	      newReg = false;
	      _ele_RCTetVect[counter][iReg] +=  iET;
	    }
	  }
	}

	// Matching with L1 candidates (in this region)
	if(newReg) {
	  _ele_RCTetaVect[counter][nReg] = iReta;
	  _ele_RCTphiVect[counter][nReg] = iRphi;
	  _ele_RCTetVect[counter][nReg]  =  iET;
	  
	  // standard L1 Candidates collection
	  for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
	    if(_trig_L1emIso_iphi[il1] == iRphi && _trig_L1emIso_ieta[il1] == iReta) 
	      _ele_RCTL1isoVect[counter][nReg] = _trig_L1emIso_rank[il1];
	  }
	  for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	    if(_trig_L1emNonIso_iphi[il1] == iRphi && _trig_L1emNonIso_ieta[il1] == iReta) 
	      _ele_RCTL1nonisoVect[counter][nReg] = _trig_L1emNonIso_rank[il1];
	  }
	  
	  // modified L1 Candidates collection
	  for(int il1=0; il1<_trig_L1emIso_N_M; ++il1) {
	    if(_trig_L1emIso_iphi_M[il1] == iRphi && _trig_L1emIso_ieta_M[il1] == iReta) 
	      _ele_RCTL1isoVect_M[counter][nReg] = _trig_L1emIso_rank_M[il1];
	  }
	  for(int il1=0; il1<_trig_L1emNonIso_N_M; ++il1) {
	    if(_trig_L1emNonIso_iphi_M[il1] == iRphi && _trig_L1emNonIso_ieta_M[il1] == iReta) 
	      _ele_RCTL1nonisoVect_M[counter][nReg] = _trig_L1emNonIso_rank_M[il1];
	  }
					
	  nReg++;
	} // if newReg

      } // loop over RecHits in cluster
    } // loop over clusters in SC

    // Determine this electron's main region (containing most energetic tower)
    if(PrintDebug_) cout << "- main ele RCT region => main L1 candidate" << endl;
    double TTetmax2 = 0.;
    int iTTmax2     = -1.;
    int TTetamax2, TTphimax2;

    // reject masked towers
    bool nomaskTT = true; 
    for (int iTow=0; iTow<nTow; ++iTow) {
      nomaskTT = true;
      for(int iMask=0 ; iMask<_trig_nMaskedCh ; iMask++)
	if( _ele_TTetaVect[counter][iTow]==_trig_iMaskedTTeta[iMask] && _ele_TTetaVect[counter][iTow]==_trig_iMaskedTTphi[iMask] )
	  nomaskTT=false;
      if(nomaskTT && _ele_TTetVect[counter][iTow] > TTetmax2) {
	iTTmax2 = iTow;
	TTetmax2 = _ele_TTetVect[counter][iTow];
      } 
    } // loop over towers

    // identify max eT tower's RCT region
    if(iTTmax2>=0) {
      TTetamax2 = getGCTRegionEta(_ele_TTetaVect[counter][iTTmax2]);
      TTphimax2 = getGCTRegionPhi(_ele_TTphiVect[counter][iTTmax2]);
      _ele_RCTeta[counter] = TTetamax2;
      _ele_RCTphi[counter] = TTphimax2;
    }

    // match L1 candidate in this region to this ele
    for(int iReg=0 ; iReg<10 ; iReg++) {
      if( _ele_RCTetaVect[counter][iReg]==_ele_RCTeta[counter] && _ele_RCTphiVect[counter][iReg]==_ele_RCTphi[counter] ) {
	_ele_RCTL1iso[counter] = _ele_RCTL1isoVect[counter][iReg];
	_ele_RCTL1noniso[counter] = _ele_RCTL1nonisoVect[counter][iReg];
      }
      if( _ele_RCTetaVect[counter][iReg]==_ele_RCTeta[counter] && _ele_RCTphiVect[counter][iReg]==_ele_RCTphi[counter] ) {
	_ele_RCTL1iso_M[counter] = _ele_RCTL1isoVect_M[counter][iReg];
	_ele_RCTL1noniso_M[counter] = _ele_RCTL1nonisoVect_M[counter][iReg];
      }
    } // end loop over electron's regions

    // CHECKING SPIKY ELECTRONS //
    if(PrintDebug_) cout << "- check spiky electrons" << endl;
    if ((*EleHandle)[i].isEB()) {
      reducedRecHits = reducedEBRecHits.product() ; 

      // Seed cluster
      const edm::Ptr<reco::CaloCluster> & seedCluster = (*EleHandle)[i].superCluster()->seed();  
      std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(),reducedRecHits);
      ele_severityLevelSeed[counter] = sl->severityLevel(id.first,*reducedRecHits);

      // Other clusters
      for (reco::CaloCluster_iterator bc = (*EleHandle)[i].superCluster()->clustersBegin(); 
	   bc!=(*EleHandle)[i].superCluster()->clustersEnd(); ++bc) {
	if ( seedCluster==(*bc) ) continue;
	std::pair<DetId, float> id = EcalClusterTools::getMaximum((*bc)->hitsAndFractions(),reducedRecHits);
	ele_severityLevelClusters[counter] = sl->severityLevel(id.first,*reducedRecHits);
      }

    } // endif isEB

  } // END LOOP OVER ELECTRONS

  if(PrintDebug_) cout << "--- END FillEle ---" << endl;

}

// ====================================================================================
void ElectronL1Study::FillSuperClusters(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{

  ///////////
  // SETUP //
  ///////////
  
  // Reduced RecHits //
  edm::Handle< EcalRecHitCollection > reducedEBRecHits;
  edm::Handle< EcalRecHitCollection > reducedEERecHits;
  if(!aod_){
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEE"), reducedEERecHits );
  }
  else{
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEE"), reducedEERecHits );
  }
  const EcalRecHitCollection * reducedRecHits;

  // Severity //
  unsigned long long cacheSevLevel = 0;
  edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
  if(cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()){
    cacheSevLevel = iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
  }
  const EcalSeverityLevelAlgo* sl=sevLevel.product();

  // H/E //
  towersH_ = new edm::Handle<CaloTowerCollection>() ;
  if (!iEvent.getByLabel(hcalTowers_,*towersH_))
    { edm::LogError("ElectronHcalHelper::readEvent")<<"failed to get the hcal towers of label "<<hcalTowers_ ; }
  towerIso1_ = new EgammaTowerIsolation(hOverEConeSize_,0.,hOverEPtMin_,1,towersH_->product()) ;
  towerIso2_ = new EgammaTowerIsolation(hOverEConeSize_,0.,hOverEPtMin_,2,towersH_->product()) ;

  // Beam Spot // --> to remove ? 
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  iEvent.getByType(recoBeamSpotHandle) ;
  const reco::BeamSpot bs = *recoBeamSpotHandle ;

  // Isolation
  edm::Handle<TrackCollection> ctfTracksH;  
  iEvent.getByLabel("generalTracks", ctfTracksH); // ctfTracks_

  // Iso Track
  double isolationtrackThresholdB_Barrel = 0.7;     //0.0; 
  double TrackConeOuterRadiusB_Barrel    = 0.3; 
  double TrackConeInnerRadiusB_Barrel    = 0.015;   //0.04;
  double isolationtrackEtaSliceB_Barrel  = 0.015;
  double longImpactParameterB_Barrel     = 0.2;  
  double transImpactParameterB_Barrel    = 999999.; //0.1;  
  
  double isolationtrackThresholdB_Endcap  = 0.7;    // 0.0
  double TrackConeOuterRadiusB_Endcap     = 0.3;
  double TrackConeInnerRadiusB_Endcap     = 0.015;  //0.04;
  double isolationtrackEtaSliceB_Endcap   = 0.015;
  double longImpactParameterB_Endcap      = 0.2;
  double transImpactParameterB_Endcap     = 999999.; //0.1;

  // Iso HCAL
  float egHcalIsoConeSizeOutSmall=0.3;
  int egHcalDepth1=1, egHcalDepth2=2; //float egHcalIsoConeSizeIn=intRadiusHcal_,egHcalIsoPtMin=etMinHcal_;
  double egHcalIsoConeSizeIn = 0.15;  //intRadiusHcal   = 0.15;
  double egHcalIsoPtMin      = 0.0;   //  etMinHcal 
  EgammaTowerIsolation hadDepth1Isolation03(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,towersH_->product()) ;
  EgammaTowerIsolation hadDepth2Isolation03(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,towersH_->product()) ;

  // Iso ECAL
  double egIsoConeSizeInBarrel  = 3.0; // intRadiusEcalBarrel 
  double egIsoConeSizeInEndcap  = 3.0;  // intRadiusEcalEndcaps
  double egIsoJurassicWidth     = 1.5;  // jurassicWidth 
  double egIsoPtMinBarrel       = 0.0;  // etMinBarrel
  double egIsoEMinBarrel        = 0.08; // eMinBarrel
  double egIsoPtMinEndcap       = 0.1;  // etMinEndcaps
  double egIsoEMinEndcap        = 0.0;  // egIsoEMinEndcaps
  bool vetoClustered   = false;  
  bool useNumCrystals  = true;  
  float egIsoConeSizeOutSmall=0.3; //, egIsoConeSizeOutLarge=0.4;


  ////////////
  // BARREL //
  ////////////
  edm::Handle<reco::SuperClusterCollection> sc_coll_EB;
  iEvent.getByLabel(edm::InputTag("correctedHybridSuperClusters"), sc_coll_EB);
  _sc_N[0] = sc_coll_EB->size();
  int index_sc = 0;

  int nTow=0;
  int nReg=0;
  double R, Rt;

  // Loop over SCs //
  for( reco::SuperClusterCollection::const_iterator isc=sc_coll_EB->begin(); isc!=sc_coll_EB->end(); isc++) {
    if(index_sc>24) continue;

    R  = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y() +isc->z()*isc->z());
    Rt = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y());
    
    _sc_E[index_sc]   = isc->energy();
    _sc_Et[index_sc]  = isc->energy()*(Rt/R);
    _sc_Eta[index_sc] = isc->eta();
    _sc_Phi[index_sc] = isc->phi();

    _sc_EB_EE[index_sc] = 0; // this SC is in EB

    reducedRecHits = reducedEBRecHits.product(); 
    
    // Seed cluster's severity
    const edm::Ptr<reco::CaloCluster> & seedCluster = isc->seed(); //(*EleHandle)[i].superCluster()->seed() ;  
    std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(),reducedRecHits);
    _sc_severityLevelSeed[index_sc] =  sl->severityLevel(id.first,*reducedRecHits);

    // Variables of interest
    // H/E
    reco::SuperCluster EmSCCand = *isc;
    double HoE = towerIso1_->getTowerESum(&EmSCCand) + towerIso2_->getTowerESum(&EmSCCand) ;
    HoE /= 	isc->energy() ;     
    _sc_he[index_sc] = HoE ;

    // SigmaIetaIeta
    const reco::CaloCluster & seedCluster1 = *(isc->seed());
    std::vector<float> localCovariances = EcalClusterTools::localCovariances(seedCluster1,reducedRecHits,topology_) ;
    _sc_sigmaietaieta[index_sc]  = sqrt(localCovariances[0]) ;

    // Iso HCAL
    _sc_hcalDepth1TowerSumEt_dr03[index_sc] = hadDepth1Isolation03.getTowerEtSum(&EmSCCand);
    _sc_hcalDepth2TowerSumEt_dr03[index_sc] = hadDepth2Isolation03.getTowerEtSum(&EmSCCand);

    // Iso ECAL
    EcalRecHitMetaCollection ecalBarrelHits(*reducedEBRecHits);
    EgammaRecHitIsolation ecalBarrelIsol03( egIsoConeSizeOutSmall,egIsoConeSizeInBarrel,egIsoJurassicWidth,
					    egIsoPtMinBarrel,egIsoEMinBarrel,theCaloGeom_,&ecalBarrelHits,
					    sevLevel.product(),DetId::Ecal );
    ecalBarrelIsol03.setUseNumCrystals(useNumCrystals);
    ecalBarrelIsol03.setVetoClustered(vetoClustered);

    reco::RecoEcalCandidate * cand = new RecoEcalCandidate();
    math::XYZPoint v(0,0,0); math::XYZVector p = isc->energy() * (isc->position() -v).unit(); double t = sqrt(0. + p.mag2());
    cand->setCharge(0); cand->setVertex(v); cand->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), t));		
    const reco::SuperClusterRef sc_ref(sc_coll_EB, index_sc);
    cand->setSuperCluster(sc_ref);
    
    _sc_ecalRecHitSumEt_dr03[index_sc] = ecalBarrelIsol03.getEtSum(cand);

    // Iso Tracker (Calculate hollow cone track isolation, CONE 0.3)
    reco::Photon * newPhoton = new Photon(); 
    newPhoton->setVertex(v); newPhoton->setCharge(0); newPhoton->setMass(0);
    newPhoton->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), isc->energy()));
    
    double trkiso_hc_03 = 0.;
    int counter = 0;  double ptSum = 0.;
    
    PhotonTkIsolation phoIso(TrackConeOuterRadiusB_Barrel, //RCone, 
			     TrackConeInnerRadiusB_Barrel, //RinnerCone, 
			     isolationtrackEtaSliceB_Barrel, //etaSlice,  
			     isolationtrackThresholdB_Barrel, //pTThresh, 
			     longImpactParameterB_Barrel, //lip , 
			     transImpactParameterB_Barrel, //d0, 
			     ctfTracksH.product(), //trackCollection, ctfTracksH.product(),bs.position()
			     bs.position());       //math::XYZPoint(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0()));
    
    counter  = phoIso.getNumberTracks(newPhoton);
    ptSum    = phoIso.getPtTracks(newPhoton);
    trkiso_hc_03 = ptSum;
    
    _sc_trkiso_dr03[index_sc] = trkiso_hc_03;

    // L1 CANDIDATE MATCHING //
    nTow=0;
    nReg=0;
    // loop over clusters in SC
    for (reco::CaloCluster_iterator clus = isc->clustersBegin () ;
	 clus != isc->clustersEnd () ; ++clus) {
      
      std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions();

      // loop over RecHits (Xtals) in cluster
      for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ;
	   detitr != clusterDetIds.end () ; ++detitr) {
	
	// check located in ECAL
	if ( (detitr -> first).det () != DetId::Ecal) {
	  std::cout << " det is " << (detitr -> first).det () << std::endl ;
	  continue ;
	}
	
	// Get the RecHit and its tower
	EcalRecHitCollection::const_iterator thishit = reducedRecHits->find ( (detitr -> first) ) ;
	if (thishit == reducedRecHits->end ()) continue;
	EcalRecHit myhit = (*thishit) ;
	EBDetId detid(thishit->id());
	EcalTrigTowerDetId towid= detid.tower();
	float thetahit =  theBarrelGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	
	// Tower informations
	int iETA   = towid.ieta();
	int iPHI   = towid.iphi();
	int iReta  = getGCTRegionEta(iETA);
	int iRphi  = getGCTRegionPhi(iPHI);
	double iET = myhit.energy()*sin(thetahit);
	//
	bool newTow = true;
	if(nTow>0) {
	  for (int iTow=0; iTow<nTow; ++iTow) {
	    if(_sc_TTetaVect[index_sc][iTow] == iETA && _sc_TTphiVect[index_sc][iTow] == iPHI) {
	      newTow = false;
	      _sc_TTetVect[index_sc][iTow] +=  iET;
	    }
	  }
	}
	if(newTow) {
	  _sc_TTetaVect[index_sc][nTow] = iETA;
	  _sc_TTphiVect[index_sc][nTow] = iPHI;
	  _sc_TTetVect[index_sc][nTow] =  iET;
	  nTow++;
	}

	// RCT region informations
	bool newReg = true;
	if(nReg>0) {
	  for (int iReg=0; iReg<nReg; ++iReg) {
	    if(_sc_RCTetaVect[index_sc][iReg] == iReta && _sc_RCTphiVect[index_sc][iReg] == iRphi) {
	      newReg = false;
	      _sc_RCTetVect[index_sc][iReg] +=  iET;
	    }
	  }
	}

	// Matching with L1 candidates (in this region)
	if(newReg) {
	  _sc_RCTetaVect[index_sc][nReg] = iReta;
	  _sc_RCTphiVect[index_sc][nReg] = iRphi;
	  _sc_RCTetVect[index_sc][nReg]  =  iET;
	  
	  // standard L1 Candidates collection
	  for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
	    if(_trig_L1emIso_iphi[il1] == iRphi && _trig_L1emIso_ieta[il1] == iReta) 
	      _sc_RCTL1isoVect[index_sc][nReg] = _trig_L1emIso_rank[il1];
	  }
	  for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	    if(_trig_L1emNonIso_iphi[il1] == iRphi && _trig_L1emNonIso_ieta[il1] == iReta) 
	      _sc_RCTL1nonisoVect[index_sc][nReg] = _trig_L1emNonIso_rank[il1];
	  }
	  
	  // modified L1 Candidates collection
	  for(int il1=0; il1<_trig_L1emIso_N_M; ++il1) {
	    if(_trig_L1emIso_iphi_M[il1] == iRphi && _trig_L1emIso_ieta_M[il1] == iReta) 
	      _sc_RCTL1isoVect_M[index_sc][nReg] = _trig_L1emIso_rank_M[il1];
	  }
	  for(int il1=0; il1<_trig_L1emNonIso_N_M; ++il1) {
	    if(_trig_L1emNonIso_iphi_M[il1] == iRphi && _trig_L1emNonIso_ieta_M[il1] == iReta) 
	      _sc_RCTL1nonisoVect_M[index_sc][nReg] = _trig_L1emNonIso_rank_M[il1];
	  }
					
	  nReg++;
	} // if newReg

      } // loop over RecHits in cluster
    } // loop over clusters in SC

    // Determine this SC's main region (containing most energetic tower)
    double TTetmax2 = 0.;
    int iTTmax2     = -1.;
    int TTetamax2, TTphimax2;

    // reject masked towers
    bool nomaskTT = true; 
    for (int iTow=0; iTow<nTow; ++iTow) {
      nomaskTT = true;
      for(int iMask=0 ; iMask<_trig_nMaskedCh ; iMask++)
	if( _sc_TTetaVect[index_sc][iTow]==_trig_iMaskedTTeta[iMask] && _sc_TTetaVect[index_sc][iTow]==_trig_iMaskedTTphi[iMask] )
	  nomaskTT=false;
      if(nomaskTT && _sc_TTetVect[index_sc][iTow] > TTetmax2) {
	iTTmax2 = iTow;
	TTetmax2 = _sc_TTetVect[index_sc][iTow];
      } 
    } // loop over towers

    // identify max eT tower's RCT region
    if(iTTmax2>=0) {
      TTetamax2 = getGCTRegionEta(_sc_TTetaVect[index_sc][iTTmax2]);
      TTphimax2 = getGCTRegionPhi(_sc_TTphiVect[index_sc][iTTmax2]);
      _sc_RCTeta[index_sc] = TTetamax2;
      _sc_RCTphi[index_sc] = TTphimax2;
    }

    // match L1 candidate in this region to this SC
    for(int iReg=0 ; iReg<10 ; iReg++) {
      if( _sc_RCTetaVect[index_sc][iReg]==_sc_RCTeta[index_sc] && _sc_RCTphiVect[index_sc][iReg]==_sc_RCTphi[index_sc] ) {
	_sc_RCTL1iso[index_sc] = _sc_RCTL1isoVect[index_sc][iReg];
	_sc_RCTL1noniso[index_sc] = _sc_RCTL1nonisoVect[index_sc][iReg];
      }
      if( _sc_RCTetaVect[index_sc][iReg]==_sc_RCTeta[index_sc] && _sc_RCTphiVect[index_sc][iReg]==_sc_RCTphi[index_sc] ) {
	_sc_RCTL1iso_M[index_sc] = _sc_RCTL1isoVect_M[index_sc][iReg];
	_sc_RCTL1noniso_M[index_sc] = _sc_RCTL1nonisoVect_M[index_sc][iReg];
      }
    } // end loop over electron's regions

    index_sc++;
  } // end loop over EB SuperClusters
  if(index_sc>24) { _sc_N[0] = 25; cout << "Number of SuperCluster > 25; _sc_N set to 25" << endl;}

  /////////////
  // ENDCAPS //
  /////////////

  edm::Handle<reco::SuperClusterCollection> sc_coll_EE;
  iEvent.getByLabel(edm::InputTag("correctedMulti5x5SuperClustersWithPreshower"), sc_coll_EE);
  _sc_N[1] = sc_coll_EE->size();
  int index_sc_EE = 25;
  
  for( reco::SuperClusterCollection::const_iterator isc=sc_coll_EE->begin(); isc!=sc_coll_EE->end(); isc++) { 
    if(index_sc_EE>49) continue;
    
    reducedRecHits = reducedEERecHits.product() ; 
    R  = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y() +isc->z()*isc->z());
    Rt = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y());
    
    _sc_E[index_sc_EE]   = isc->energy();
    _sc_Et[index_sc_EE]  = isc->energy()*(Rt/R);
    _sc_Eta[index_sc_EE] = isc->eta();
    _sc_Phi[index_sc_EE] = isc->phi();

    _sc_EB_EE[index_sc_EE] = 1; // this SC is in EE

    const EcalRecHitCollection * reducedRecHits = reducedEERecHits.product(); 
    
    // Variables of interest
    // H/E
    reco::SuperCluster EmSCCand = *isc;
    double HoE = towerIso1_->getTowerESum(&EmSCCand) + towerIso2_->getTowerESum(&EmSCCand) ;
    HoE /= 	isc->energy() ;     
    _sc_he[index_sc_EE] = HoE ;

    // SigmaIetaIeta
    const reco::CaloCluster & seedCluster1 = *(isc->seed());
    std::vector<float> localCovariances = EcalClusterTools::localCovariances(seedCluster1,reducedRecHits,topology_) ;
    _sc_sigmaietaieta[index_sc_EE]  = sqrt(localCovariances[0]) ;

    // Iso HCAL
    _sc_hcalDepth1TowerSumEt_dr03[index_sc_EE] = hadDepth1Isolation03.getTowerEtSum(&EmSCCand);
    _sc_hcalDepth2TowerSumEt_dr03[index_sc_EE] = hadDepth2Isolation03.getTowerEtSum(&EmSCCand);

    // Iso ECAL
    EcalRecHitMetaCollection ecalEndcapHits(*reducedEERecHits);
    EgammaRecHitIsolation ecalEndcapIsol03(egIsoConeSizeOutSmall,egIsoConeSizeInEndcap,egIsoJurassicWidth,
					   egIsoPtMinEndcap,egIsoEMinEndcap,theCaloGeom_,&ecalEndcapHits,sevLevel.product(),DetId::Ecal);
    ecalEndcapIsol03.setUseNumCrystals(useNumCrystals);
    ecalEndcapIsol03.setVetoClustered(vetoClustered);
    reco::RecoEcalCandidate * cand = new RecoEcalCandidate();
    math::XYZPoint v(0,0,0); math::XYZVector p = isc->energy() * (isc->position() -v).unit(); double t = sqrt(0. + p.mag2());
    cand->setCharge(0); cand->setVertex(v);
    cand->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), t));		
    const reco::SuperClusterRef sc_ref(sc_coll_EE, index_sc_EE);
    cand->setSuperCluster(sc_ref);
		
    _sc_ecalRecHitSumEt_dr03[index_sc_EE] = ecalEndcapIsol03.getEtSum(cand);
    
    // Track Isolation (calculate hollow cone track isolation, CONE 0.3)
    reco::Photon * newPhoton = new Photon(); 
    newPhoton->setVertex(v); newPhoton->setCharge(0); newPhoton->setMass(0);
    newPhoton->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), isc->energy()));
    
    double trkiso_hc_03 = 0.;
    int counter = 0;  double ptSum = 0.;
    
    PhotonTkIsolation phoIso(TrackConeOuterRadiusB_Endcap, //RCone, 
			     TrackConeInnerRadiusB_Endcap, //RinnerCone, 
			     isolationtrackEtaSliceB_Endcap, //etaSlice,  
			     isolationtrackThresholdB_Endcap, //pTThresh, 
			     longImpactParameterB_Endcap, //lip , 
			     transImpactParameterB_Endcap, //d0, 
			     ctfTracksH.product(), //trackCollection, ctfTracksH.product(),bs.position()
			     bs.position());       //math::XYZPoint(vertexBeamSpot.x0(),vertexBeamSpot.y0(),vertexBeamSpot.z0()));
    
    counter  = phoIso.getNumberTracks(newPhoton);
    ptSum    = phoIso.getPtTracks(newPhoton);
    trkiso_hc_03 = ptSum;
    
    _sc_trkiso_dr03[index_sc_EE] = trkiso_hc_03;

        // L1 CANDIDATE MATCHING //
    nTow=0;
    nReg=0;
    // loop over clusters in SC
    for (reco::CaloCluster_iterator clus = isc->clustersBegin () ;
	 clus != isc->clustersEnd () ; ++clus) {
      
      std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions();

      // loop over RecHits (Xtals) in cluster
      for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ;
	   detitr != clusterDetIds.end () ; ++detitr) {
	
	// check located in ECAL
	if ( (detitr -> first).det () != DetId::Ecal) {
	  std::cout << " det is " << (detitr -> first).det () << std::endl ;
	  continue ;
	}
	
	// Get the RecHit and its tower
	EcalRecHitCollection::const_iterator thishit = reducedRecHits->find ( (detitr -> first) ) ;
	if (thishit == reducedRecHits->end ()) continue;
	EcalRecHit myhit = (*thishit) ;
	EEDetId detid(thishit->id());
	EcalTrigTowerDetId towid = (*eTTmap_).towerOf(detid);
	float thetahit =  theBarrelGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	
	// Tower informations
	int iETA   = towid.ieta();
	int iPHI   = towid.iphi();
	int iReta  = getGCTRegionEta(iETA);
	int iRphi  = getGCTRegionPhi(iPHI);
	double iET = myhit.energy()*sin(thetahit);
	//
	bool newTow = true;
	if(nTow>0) {
	  for (int iTow=0; iTow<nTow; ++iTow) {
	    if(_sc_TTetaVect[index_sc_EE][iTow] == iETA && _sc_TTphiVect[index_sc_EE][iTow] == iPHI) {
	      newTow = false;
	      _sc_TTetVect[index_sc_EE][iTow] +=  iET;
	    }
	  }
	}
	if(newTow) {
	  _sc_TTetaVect[index_sc_EE][nTow] = iETA;
	  _sc_TTphiVect[index_sc_EE][nTow] = iPHI;
	  _sc_TTetVect[index_sc_EE][nTow] =  iET;
	  nTow++;
	}

	// RCT region informations
	bool newReg = true;
	if(nReg>0) {
	  for (int iReg=0; iReg<nReg; ++iReg) {
	    if(_sc_RCTetaVect[index_sc_EE][iReg] == iReta && _sc_RCTphiVect[index_sc_EE][iReg] == iRphi) {
	      newReg = false;
	      _sc_RCTetVect[index_sc_EE][iReg] +=  iET;
	    }
	  }
	}

	// Matching with L1 candidates (in this region)
	if(newReg) {
	  _sc_RCTetaVect[index_sc_EE][nReg] = iReta;
	  _sc_RCTphiVect[index_sc_EE][nReg] = iRphi;
	  _sc_RCTetVect[index_sc_EE][nReg]  =  iET;
	  
	  // standard L1 Candidates collection
	  for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
	    if(_trig_L1emIso_iphi[il1] == iRphi && _trig_L1emIso_ieta[il1] == iReta) 
	      _sc_RCTL1isoVect[index_sc_EE][nReg] = _trig_L1emIso_rank[il1];
	  }
	  for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	    if(_trig_L1emNonIso_iphi[il1] == iRphi && _trig_L1emNonIso_ieta[il1] == iReta) 
	      _sc_RCTL1nonisoVect[index_sc_EE][nReg] = _trig_L1emNonIso_rank[il1];
	  }
	  
	  // modified L1 Candidates collection
	  for(int il1=0; il1<_trig_L1emIso_N_M; ++il1) {
	    if(_trig_L1emIso_iphi_M[il1] == iRphi && _trig_L1emIso_ieta_M[il1] == iReta) 
	      _sc_RCTL1isoVect_M[index_sc_EE][nReg] = _trig_L1emIso_rank_M[il1];
	  }
	  for(int il1=0; il1<_trig_L1emNonIso_N_M; ++il1) {
	    if(_trig_L1emNonIso_iphi_M[il1] == iRphi && _trig_L1emNonIso_ieta_M[il1] == iReta) 
	      _sc_RCTL1nonisoVect_M[index_sc_EE][nReg] = _trig_L1emNonIso_rank_M[il1];
	  }
					
	  nReg++;
	} // if newReg

      } // loop over RecHits in cluster
    } // loop over clusters in SC

    // Determine this SC's main region (containing most energetic tower)
    double TTetmax2 = 0.;
    int iTTmax2     = -1.;
    int TTetamax2, TTphimax2;

    // reject masked towers
    bool nomaskTT = true; 
    for (int iTow=0; iTow<nTow; ++iTow) {
      nomaskTT = true;
      for(int iMask=0 ; iMask<_trig_nMaskedCh ; iMask++)
	if( _sc_TTetaVect[index_sc_EE][iTow]==_trig_iMaskedTTeta[iMask] && _sc_TTetaVect[index_sc_EE][iTow]==_trig_iMaskedTTphi[iMask] )
	  nomaskTT=false;
      if(nomaskTT && _sc_TTetVect[index_sc_EE][iTow] > TTetmax2) {
	iTTmax2 = iTow;
	TTetmax2 = _sc_TTetVect[index_sc_EE][iTow];
      } 
    } // loop over towers

    // identify max eT tower's RCT region
    if(iTTmax2>=0) {
      TTetamax2 = getGCTRegionEta(_sc_TTetaVect[index_sc_EE][iTTmax2]);
      TTphimax2 = getGCTRegionPhi(_sc_TTphiVect[index_sc_EE][iTTmax2]);
      _sc_RCTeta[index_sc_EE] = TTetamax2;
      _sc_RCTphi[index_sc_EE] = TTphimax2;
    }

    // match L1 candidate in this region to this SC
    for(int iReg=0 ; iReg<10 ; iReg++) {
      if( _sc_RCTetaVect[index_sc_EE][iReg]==_sc_RCTeta[index_sc_EE] && _sc_RCTphiVect[index_sc_EE][iReg]==_sc_RCTphi[index_sc_EE] ) {
	_sc_RCTL1iso[index_sc_EE] = _sc_RCTL1isoVect[index_sc_EE][iReg];
	_sc_RCTL1noniso[index_sc_EE] = _sc_RCTL1nonisoVect[index_sc_EE][iReg];
      }
      if( _sc_RCTetaVect[index_sc_EE][iReg]==_sc_RCTeta[index_sc_EE] && _sc_RCTphiVect[index_sc_EE][iReg]==_sc_RCTphi[index_sc_EE] ) {
	_sc_RCTL1iso_M[index_sc_EE] = _sc_RCTL1isoVect_M[index_sc_EE][iReg];
	_sc_RCTL1noniso_M[index_sc_EE] = _sc_RCTL1nonisoVect_M[index_sc_EE][iReg];
      }
    } // end loop over electron's regions

    index_sc_EE++;
  } // end loop over EB SuperClusters
  if(index_sc_EE>24) { _sc_N[1] = 25; cout << "Number of SuperCluster > 25; _sc_N set to 25" << endl;}
 
} 

// ====================================================================================
void ElectronL1Study::Init()
// ====================================================================================
{
  
  // EVENTS //
  nEvent = 0;
  nRun = 0;
  nLumi = 0;


  // VERTICES //
  _vtx_N = 0; 
  for(int iv=0;iv<200;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }


  // TRIGGER //

  // L1 candidates
  _trig_L1emIso_N    = 0; 
  _trig_L1emNonIso_N = 0;
  _trig_L1emIso_N_M  = 0; 
  _trig_L1emNonIso_N_M = 0;
  _trig_preL1emIso_N   = 0; 
  _trig_preL1emNonIso_N  = 0;
  _trig_postL1emIso_N    = 0; 
  _trig_postL1emNonIso_N = 0;

  for(int il1=0 ; il1<4 ; il1++) {
    _trig_L1emIso_ieta[il1] = 0; 
    _trig_L1emIso_iphi[il1] = 0; 
    _trig_L1emIso_rank[il1] = 0; 
    _trig_L1emIso_eta[il1]    = 0.; 
    _trig_L1emIso_phi[il1]   = 0.; 
    _trig_L1emIso_energy[il1] = 0.; 
    _trig_L1emIso_et[il1]     = 0.; 

    _trig_L1emIso_ieta_M[il1] = 0; 
    _trig_L1emIso_iphi_M[il1] = 0; 
    _trig_L1emIso_rank_M[il1] = 0; 
    _trig_L1emIso_eta_M[il1]    = 0.; 
    _trig_L1emIso_phi_M[il1]   = 0.; 
    _trig_L1emIso_energy_M[il1] = 0.; 
    _trig_L1emIso_et_M[il1]     = 0.; 
		
    _trig_L1emNonIso_ieta[il1] = 0; 
    _trig_L1emNonIso_iphi[il1] = 0; 
    _trig_L1emNonIso_rank[il1] = 0; 
    _trig_L1emNonIso_eta[il1]    = 0.; 
    _trig_L1emNonIso_phi[il1]   = 0.; 
    _trig_L1emNonIso_energy[il1] = 0.; 
    _trig_L1emNonIso_et[il1]     = 0.; 

    _trig_L1emNonIso_ieta_M[il1] = 0; 
    _trig_L1emNonIso_iphi_M[il1] = 0; 
    _trig_L1emNonIso_rank_M[il1] = 0; 
    _trig_L1emNonIso_eta_M[il1]    = 0.; 
    _trig_L1emNonIso_phi_M[il1]   = 0.; 
    _trig_L1emNonIso_energy_M[il1] = 0.; 
    _trig_L1emNonIso_et_M[il1]     = 0.; 
		
    _trig_preL1emIso_ieta[il1] = 0; 
    _trig_preL1emIso_iphi[il1] = 0; 
    _trig_preL1emIso_rank[il1] = 0;
    _trig_preL1emNonIso_ieta[il1] = 0; 
    _trig_preL1emNonIso_iphi[il1] = 0; 
    _trig_preL1emNonIso_rank[il1] = 0; 
		
    _trig_postL1emIso_ieta[il1] = 0; 
    _trig_postL1emIso_iphi[il1] = 0; 
    _trig_postL1emIso_rank[il1] = 0;
    _trig_postL1emNonIso_ieta[il1] = 0; 
    _trig_postL1emNonIso_iphi[il1] = 0; 
    _trig_postL1emNonIso_rank[il1] = 0;  
  }

  // Masked Towers
  _trig_nMaskedRCT=0;
  _trig_nMaskedCh=0;
	
  for (int ii=0;ii<100;ii++) {
    _trig_iMaskedRCTeta[ii]   = -999;
    _trig_iMaskedRCTphi[ii]   = -999;
    _trig_iMaskedRCTcrate[ii] = -999;
    _trig_iMaskedTTeta[ii]    = -999;
    _trig_iMaskedTTphi[ii]    = -999;
  }

  // Ecal trigger primitives
  _trig_tower_N = 0;
  _trig_tower_N_M = 0;
  _trig_tower_N_E = 0;

  for(int i=0 ; i<4032 ; i++) {
    _trig_tower_ieta[i]=-999;
    _trig_tower_iphi[i]=-999;
    _trig_tower_adc[i]=-999;
    _trig_tower_sFGVB[i]=-999;

    _trig_tower_ieta_M[i]=-999;
    _trig_tower_iphi_M[i]=-999;
    _trig_tower_adc_M[i]=-999;
    _trig_tower_sFGVB_M[i]=-999;

    _trig_tower_ieta_E[i]=-999;
    _trig_tower_iphi_E[i]=-999;
    for(int j=0 ; j<5 ; j++) {
      _trig_tower_adc_E[i][j]=-999;
      _trig_tower_sFGVB_E[i][j]=-999;
    }
  }

  _trig_isEleHLTpath = 0;

  for(int i=0 ; i<250 ; i++) 
    trig_hltInfo[i] = 0;

  for(int i=0 ; i<4 ; i++)
    trig_HLT_path[i] = 0;


  // ELECTRONS //
  ele_N = 0;
  for(int iEle=0 ; iEle<10 ; iEle++) {
    // L1 Candidate matching
    _ele_RCTeta[iEle]      = -999;
    _ele_RCTphi[iEle]      = -999;
    _ele_RCTL1iso[iEle]    = -999;
    _ele_RCTL1noniso[iEle] = -999;
    _ele_RCTL1iso_M[iEle]    = -999;
    _ele_RCTL1noniso_M[iEle] = -999;

    for(int icc = 0; icc < 50; ++icc) {
      _ele_TTetaVect[iEle][icc] = -999;
      _ele_TTphiVect[iEle][icc] = -999;
      _ele_TTetVect[iEle][icc] = 0.;
    }
    for(int icc = 0; icc < 10; ++icc) {
      _ele_RCTetaVect[iEle][icc] = -999;
      _ele_RCTphiVect[iEle][icc] = -999;
      _ele_RCTetVect[iEle][icc] = 0.;
      _ele_RCTL1isoVect[iEle][icc] = -999;
      _ele_RCTL1nonisoVect[iEle][icc] = -999;
      _ele_RCTL1isoVect_M[iEle][icc] = -999;
      _ele_RCTL1nonisoVect_M[iEle][icc] = -999;
    }

    // Electron properties
    ele_isbarrel[iEle]=0;
    ele_isendcap[iEle]=0;
    ele_isEBetaGap[iEle]=0;
    ele_isEBphiGap[iEle]=0;
    ele_isEEdeeGap[iEle]=0;
    ele_isEEringGap[iEle]=0;
    ele_isecalDriven[iEle]=0;
    ele_istrackerDriven[iEle]=0;
    ele_eClass[iEle]=0;
    //
    ele_echarge[iEle]=0;
    ele_he[iEle]=0; 
    ele_sigmaietaieta[iEle]=0;
    ele_sigmaetaeta[iEle]=0;
    ele_fbrem[iEle]=0 ;
    ele_ECAL_fbrem[iEle]=0; 
    ele_mva[iEle]=0 ;
    ele_missing_hits[iEle]=0;
    //
    ele_isConversion[iEle] = 0;
    ele_convFound[iEle] = 0;
    ele_conv_dist[iEle] = 0.;
    ele_conv_dcot[iEle] = 0.;
    //
    ele_eseedpout[iEle]=0;
    ele_ep[iEle]=0;
    ele_eseedp[iEle]=0;
    ele_eelepout[iEle]=0;  
    ele_pin_mode[iEle]=0;
    ele_pout_mode[iEle]=0;
    ele_pin_mean[iEle]=0;
    ele_pout_mean[iEle]=0; 
    ele_pTin_mode[iEle]=0;
    ele_pTout_mode[iEle]=0; 
    ele_pTin_mean[iEle]=0; ; 
    ele_pTout_mean[iEle]=0;
    //
    ele_deltaetaseed[iEle]=0;
    ele_deltaetaele[iEle]=0;
    ele_deltaphiseed[iEle]=0;
    ele_deltaphiele[iEle]=0;
    ele_deltaetain[iEle]=0;
    ele_deltaphiin[iEle]=0;
    //
    ele_sclE[iEle]=0;
    ele_sclEt[iEle]=0;
    ele_sclEta[iEle]=0;
    ele_sclPhi[iEle]=0;
    ele_calo_energy[iEle]=0;
    //
    ele_track_x[iEle] = 0.;
    ele_track_y[iEle] = 0.;
    ele_track_z[iEle] = 0.;
    ele_vertex_x[iEle]=0; 
    ele_vertex_y[iEle]=0;
    ele_vertex_z[iEle]=0;
    ele_lost_hits[iEle]=0;
    ele_chi2_hits[iEle]=0;
    ele_expected_inner_hits[iEle]=-1;
    //
    ele_tkSumPt_dr03[iEle]=0;
    ele_tkSumPt_dr04[iEle]=0;
    ele_ecalRecHitSumEt_dr03[iEle]=0;
    ele_ecalRecHitSumEt_dr04[iEle]=0;
    ele_hcalDepth1TowerSumEt_dr03[iEle]=0;
    ele_hcalDepth2TowerSumEt_dr03[iEle]=0;
    ele_hcalDepth1TowerSumEt_dr04[iEle]=0;
    ele_hcalDepth2TowerSumEt_dr04[iEle]=0;
    //
    ele_severityLevelSeed[iEle]     = 0.;
    ele_severityLevelClusters[iEle] = 0.;
  }


  // SUPERCLUSTERS //
  for(int iSC=0 ; iSC<2 ; iSC++)
    _sc_N[iSC] = 0;
  for(int isc=0 ; isc<50 ; isc++) {
    _sc_E[isc]   = 0.; 
    _sc_Et[isc]  = 0.; 
    _sc_Eta[isc] = 0.; 
    _sc_Phi[isc] = 0.; 
    _sc_severityLevelSeed[isc] = 0;
    _sc_he[isc]  = -10.;
    _sc_sigmaietaieta[isc] = 0.;
    _sc_hcalDepth1TowerSumEt_dr03[isc] = 0.;
    _sc_hcalDepth2TowerSumEt_dr03[isc] = 0.;
    _sc_ecalRecHitSumEt_dr03[isc]      = 0.;
    _sc_trkiso_dr03[isc]               = 0.;
		
    _sc_RCTeta[isc]=-999;
    _sc_RCTphi[isc]=-999;
    _sc_RCTL1iso[isc]     = -999;
    _sc_RCTL1noniso[isc]  = -999;
		
    for (int li=0;li<50;li++) {
      _sc_TTetaVect[isc][li]=-999;
      _sc_TTphiVect[isc][li]=-999;
      _sc_TTetVect[isc][li]=0.;
    } 
    for (int li=0;li<10;li++) {
      _sc_RCTetaVect[isc][li]=-999;
      _sc_RCTphiVect[isc][li]=-999;
      _sc_RCTetVect[isc][li]=0.;
      _sc_RCTL1isoVect[isc][li]=-999;
      _sc_RCTL1nonisoVect[isc][li]=-999;
    } 
  }

}

// ====================================================================================
void ElectronL1Study::beginJob(const edm::ParameterSet& conf) { }
// ====================================================================================

// ====================================================================================
void ElectronL1Study::endJob() { }
// ====================================================================================
