#include "EGamma/ECGelec/plugins/SimpleNtpleCustom.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace IPTools;
//using namespace math;



// ====================================================================================
SimpleNtpleCustom::SimpleNtpleCustom(const edm::ParameterSet& iConfig) :

  // Nadir study

  /*hcalTowers_ (iConfig.getParameter<edm::InputTag>("hcalTowers")),
  hOverEConeSize_ (iConfig.getParameter<double>("hOverEConeSize")),
  hOverEPtMin_ (iConfig.getParameter<double>("hOverEPtMin")),*/
  
  nadGetL1M_ (iConfig.getUntrackedParameter<bool>("NadL1M")),
  nadGetTP_ (iConfig.getUntrackedParameter<bool>("NadTP")),
  nadGetTP_Modif_ (iConfig.getUntrackedParameter<bool>("NadTPmodif")),
  nadGetTP_Emul_ (iConfig.getUntrackedParameter<bool>("NadTPemul")),
  PrintDebug_ (iConfig.getUntrackedParameter<bool>("PrintDebug")),
  tpCollectionNormal_ (iConfig.getParameter<edm::InputTag> ("TPCollectionNormal") ),
  tpCollectionModif_ (iConfig.getParameter<edm::InputTag> ("TPCollectionModif") ),
  tpEmulatorCollection_ (iConfig.getParameter<edm::InputTag> ("TPEmulatorCollection") ),
  EleIso_TdrHzzTkMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_TdrHzzTkMapTag")),
  EleIso_TdrHzzHcalMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_TdrHzzHcalMapTag")),
  EleIso_Eg4HzzTkMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_Eg4HzzTkMapTag")),
  EleIso_Eg4HzzEcalMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_Eg4HzzEcalMapTag")),
  EleIso_Eg4HzzHcalMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_Eg4HzzHcalMapTag")),
  eleVetoIdMapTag_ (iConfig.getParameter<edm::InputTag>("eleVetoIdMap")),
  eleLooseIdMapTag_ (iConfig.getParameter<edm::InputTag>("eleLooseIdMap")),
  eleMediumIdMapTag_ (iConfig.getParameter<edm::InputTag>("eleMediumIdMap")),
  eleTightIdMapTag_ (iConfig.getParameter<edm::InputTag>("eleTightIdMap")),
  EleTag_ (iConfig.getParameter<edm::InputTag> ("EleTag")),
  PfEleTag_ ("particleFlow"),

  MuonTag_ (iConfig.getParameter<edm::InputTag> ("MuonTag")),
  MuonIso_HzzMapTag_ (iConfig.getParameter<edm::InputTag>("MuonIso_HzzMapTag")),
  MuonIsoTk_HzzMapTag_ (iConfig.getParameter<edm::InputTag>("MuonIsoTk_HzzMapTag")),
  MuonIsoEcal_HzzMapTag_ (iConfig.getParameter<edm::InputTag>("MuonIsoEcal_HzzMapTag")),
  MuonIsoHcal_HzzMapTag_ (iConfig.getParameter<edm::InputTag>("MuonIsoHcal_HzzMapTag")),
  SeedTag_ (iConfig.getParameter<edm::InputTag> ("SeedTag")),
  MCTag_  (iConfig.getParameter<edm::InputTag>("MCTag")),
  TkPTag_ (iConfig.getParameter<edm::InputTag>("TkPTag")),
  CaloJetTag_(iConfig.getParameter<edm::InputTag> ("CaloJetTag")),
  JPTJetTag_(iConfig.getParameter<edm::InputTag> ("JPTJetTag")),
  PFJetTag_(iConfig.getParameter<edm::InputTag> ("PFJetTag")),
  VerticesTag_(iConfig.getParameter<edm::InputTag> ("VerticesTag")),
  dcsTag_ (iConfig.getUntrackedParameter<edm::InputTag>("dcsTag")),
  tracksTag_("generalTracks"),
  // Trigger Stuff
  HLTTag_(iConfig.getParameter<edm::InputTag> ("HLTTag")),
  triggerEventTag_(iConfig.getParameter<edm::InputTag> ("TriggerEventTag")),
  //
  
  PileupSrc_ ("addPileupInfo"),
  RhoCorrection_("kt6PFJets:rho"),
  SigmaRhoCorrection_("kt6PFJets:sigma"),
  type_ (iConfig.getParameter<std::string>("type")),
  aod_ (iConfig.getUntrackedParameter<bool>("AOD")),
  funcname_  (iConfig.getParameter<std::string>("functionName")),
  useBeamSpot_ (iConfig.getParameter<bool>("useBeamSpot")),
  beamSpotTag_ ("offlineBeamSpot"),

  reducedEBRecHitsTag_ ("reducedEcalRecHitsEB"),
  reducedEERecHitsTag_ ("reducedEcalRecHitsEE"),

  trigger_version_emul_ (iConfig.getParameter<std::string>("trigger_version_emul")),

  ///Stage 2 Level 1
  towerTag_ (iConfig.getParameter<edm::InputTag>("towerToken")),
  mpEGTag_ (iConfig.getParameter<edm::InputTag>("mpEGToken")),
  egTag_ (iConfig.getParameter<edm::InputTag>("egToken")),
  
  towerTag_emul_ (iConfig.getParameter<edm::InputTag>("towerToken_emul")),
  clusterTag_emul_ (iConfig.getParameter<edm::InputTag>("clusterToken_emul")),
  mpEGTag_emul_ (iConfig.getParameter<edm::InputTag>("mpEGToken_emul")),
  egTag_emul_ (iConfig.getParameter<edm::InputTag>("egToken_emul")),
    
  ///Stage 1 emulator
  egTag_Stage1_emul_ (iConfig.getParameter<edm::InputTag>("egToken_Stage1_emul"))
		      
		      
{

  //now do what ever initialization is needed
  /*funcbase_ = EcalClusterFunctionFactory::get()->create( funcname_, iConfig ); 
  gtRecordCollectionTag_ = iConfig.getParameter<std::string>("GTRecordCollection") ;
	

  HLT_ElePaths_  = iConfig.getParameter<std::vector<std::string > >("HLTElePaths");
  HLT_MuonPaths_ = iConfig.getParameter<std::vector<std::string > >("HLTMuonPaths");
  HLT_Filters_   = iConfig.getParameter<std::vector<edm::InputTag > >("HLTFilters");
	
  simulation_ = iConfig.getUntrackedParameter<bool>("simulation", false);
  fillsc_     = iConfig.getUntrackedParameter<bool>("FillSC", false);*/
  


  /////////////////////////////////////////////////////////
  ///                 Stage 2 Level 1                   ///
  /////////////////////////////////////////////////////////

  m_towerToken_          = consumes<l1t::CaloTowerBxCollection>(towerTag_);
  m_mpEGToken_           = consumes<l1t::EGammaBxCollection>(mpEGTag_);
  m_egToken_             = consumes<l1t::EGammaBxCollection>(egTag_);

  m_towerToken_emul_     = consumes<l1t::CaloTowerBxCollection>(towerTag_emul_);
  m_clusterToken_emul_   = consumes<l1t::CaloClusterBxCollection>(clusterTag_emul_);
  m_mpEGToken_emul_      = consumes<l1t::EGammaBxCollection>(mpEGTag_emul_);
  m_egToken_emul_        = consumes<l1t::EGammaBxCollection>(egTag_emul_);

  m_egToken_Stage1_emul_ = consumes<l1t::EGammaBxCollection>(egTag_Stage1_emul_);


  PileupSrcToken_ = consumes<vector<PileupSummaryInfo> >(PileupSrc_);
  RhoCorrectionToken_ = consumes<double>(RhoCorrection_);
  SigmaRhoCorrectionToken_ = consumes<double>(SigmaRhoCorrection_);
  BetaCorrectionToken_ = consumes<double>(BetaCorrection_);

  VerticesToken_ = consumes<reco::VertexCollection>(VerticesTag_);
  tracksToken_ = consumes<reco::TrackCollection>(tracksTag_);
  dcsToken_ = consumes<DcsStatusCollection>(dcsTag_);

  HLTToken_ = consumes<edm::TriggerResults>(HLTTag_);
  triggerEventToken_ = consumes<trigger::TriggerEvent>(triggerEventTag_);

  tpCollectionNormalToken_ = consumes<EcalTrigPrimDigiCollection>(tpCollectionNormal_);
  tpCollectionModifToken_ = consumes<EcalTrigPrimDigiCollection>(tpCollectionModif_);
  tpEmulatorCollectionToken_ = consumes<EcalTrigPrimDigiCollection>(tpEmulatorCollection_);

  // standard collection
  emNonisolCollTag_ = iConfig.getParameter<edm::InputTag>("emNonisolCollToken");
  emIsolCollTag_ = iConfig.getParameter<edm::InputTag>("emIsolCollToken");
  // modified collection
  //emNonisolColl_M_Tag_ = iConfig.getParameter<edm::InputTag>("emNonisolColl_M_Token");
  //emIsolColl_M_Tag_ = iConfig.getParameter<edm::InputTag>("emIsolColl_M_Token");


  emNonisolCollToken_ = consumes<l1extra::L1EmParticleCollection>(emNonisolCollTag_);
  emIsolCollToken_ = consumes<l1extra::L1EmParticleCollection>(emIsolCollTag_);
  //emNonisolColl_M_Token_ = consumes<l1extra::L1EmParticleCollection>(emNonisolColl_M_Tag_);
  //emIsolColl_M_Token_ = consumes<l1extra::L1EmParticleCollection>(emIsolColl_M_Tag_);

  EleToken_ = consumes<reco::GsfElectronCollection>(EleTag_);
  PfEleToken_ = consumes<reco::PFCandidateCollection>(PfEleTag_);

  eleVetoIdMapToken_ = consumes<edm::ValueMap<bool> >(eleVetoIdMapTag_);
  eleLooseIdMapToken_ = consumes<edm::ValueMap<bool> >(eleLooseIdMapTag_);
  eleMediumIdMapToken_ = consumes<edm::ValueMap<bool> >(eleMediumIdMapTag_);
  eleTightIdMapToken_ = consumes<edm::ValueMap<bool> >(eleTightIdMapTag_);

  beamSpotToken_ = consumes<reco::BeamSpot>(beamSpotTag_);

  reducedEBRecHitsToken_ = consumes<EcalRecHitCollection>(reducedEBRecHitsTag_);
  reducedEERecHitsToken_ = consumes<EcalRecHitCollection>(reducedEERecHitsTag_);


  if(PrintDebug_) std::cout << "Creating TFileService..." << std::endl;

  edm::Service<TFileService> fs ;
  mytree_  = fs->make <TTree>("eIDSimpleTree","eIDSimpleTree"); 
	
  // Global
  /*mytree_->Branch("nEvent",&nEvent,"nEvent/I");
  mytree_->Branch("nRun",&nRun,"nRun/I");
  mytree_->Branch("nLumi",&nLumi,"nLumi/I");*/
  mytree_->Branch("nEvent",&nEvent,"nEvent/L");
  mytree_->Branch("nRun",&nRun,"nRun/L");
  mytree_->Branch("nLumi",&nLumi,"nLumi/L");
	
  // Pile UP
  mytree_->Branch("PU_N",&_PU_N,"PU_N/I");
  mytree_->Branch("PU_rhoCorr",&_PU_rho,"PU_rhoCorr/D");
  mytree_->Branch("PU_sigmaCorr",&_PU_sigma,"PU_sigmaCorr/D");

  // Vertices
  mytree_->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  mytree_->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[15]/D");
  mytree_->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[15]/D");
  mytree_->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[15]/D");
  mytree_->Branch("vtx_d0",&_vtx_d0,"vtx_d0[15]/D");
  mytree_->Branch("vtx_x",&_vtx_x,"vtx_x[15]/D");
  mytree_->Branch("vtx_y",&_vtx_y,"vtx_y[15]/D");
  mytree_->Branch("vtx_z",&_vtx_z,"vtx_z[15]/D");
	

  // Towers (original collection)
  mytree_->Branch("trig_tower_N", &_trig_tower_N, "trig_tower_N/I");
  mytree_->Branch("trig_tower_ieta",  &_trig_tower_ieta,  "trig_tower_ieta[trig_tower_N]/I");
  mytree_->Branch("trig_tower_iphi",  &_trig_tower_iphi,  "trig_tower_iphi[trig_tower_N]/I");
  mytree_->Branch("trig_tower_adc",  &_trig_tower_adc,  "trig_tower_adc[trig_tower_N]/I");
  mytree_->Branch("trig_tower_sFGVB",  &_trig_tower_sFGVB,  "trig_tower_sFGVB[trig_tower_N]/I");
	
  // Towers (cleaned collection)
  mytree_->Branch("trig_tower_N_modif", &_trig_tower_N_modif, "trig_tower_N_modif/I");
  mytree_->Branch("trig_tower_ieta_modif",  &_trig_tower_ieta_modif,  "trig_tower_ieta_modif[trig_tower_N_modif]/I");
  mytree_->Branch("trig_tower_iphi_modif",  &_trig_tower_iphi_modif,  "trig_tower_iphi_modif[trig_tower_N_modif]/I");
  mytree_->Branch("trig_tower_adc_modif",  &_trig_tower_adc_modif,  "trig_tower_adc_modif[trig_tower_N_modif]/I");
  mytree_->Branch("trig_tower_sFGVB_modif",  &_trig_tower_sFGVB_modif,  "trig_tower_sFGVB_modif[trig_tower_N_modif]/I");
	
  // Towers (emulated)
  mytree_->Branch("trig_tower_N_emul", &_trig_tower_N_emul, "trig_tower_N_emul/I");
  mytree_->Branch("trig_tower_ieta_emul",  &_trig_tower_ieta_emul,  "trig_tower_ieta_emul[trig_tower_N_emul]/I");
  mytree_->Branch("trig_tower_iphi_emul",  &_trig_tower_iphi_emul,  "trig_tower_iphi_emul[trig_tower_N_emul]/I");
  mytree_->Branch("trig_tower_adc_emul",  &_trig_tower_adc_emul,  "trig_tower_adc_emul[trig_tower_N_emul][5]/I");
  mytree_->Branch("trig_tower_sFGVB_emul",  &_trig_tower_sFGVB_emul,  "trig_tower_sFGVB_emul[trig_tower_N_emul][5]/I");
		
  // Trigger
  mytree_->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[5000]/C");
  mytree_->Branch("trig_hltInfo",&trig_hltInfo,"trig_hltInfo[500]/I");
  mytree_->Branch("trig_isUnbiased",&trig_isUnbiased,"trig_isUnbiased/I");
  mytree_->Branch("trig_isPhoton10",&trig_isPhoton10,"trig_isPhoton10/I"); 
  mytree_->Branch("trig_isPhoton15",&trig_isPhoton15,"trig_isPhoton15/I"); 
  mytree_->Branch("trig_isL1SingleEG2",&trig_isL1SingleEG2,"trig_isL1SingleEG2/I"); 
  mytree_->Branch("trig_isL1SingleEG5",&trig_isL1SingleEG5,"trig_isL1SingleEG5/I");
  mytree_->Branch("trig_isL1SingleEG8",&trig_isL1SingleEG8,"trig_isL1SingleEG8/I");
  mytree_->Branch("trig_isEle10_LW",&trig_isEle10_LW,"trig_isEle10_LW/I");
  mytree_->Branch("trig_isEle15_LW",&trig_isEle15_LW,"trig_isEle15_LW/I");
  
  mytree_->Branch("trig_isHLT_Ele30WP60_Ele8_Mass55_v2",&trig_isHLT_Ele30WP60_Ele8_Mass55_v2,"trig_isHLT_Ele30WP60_Ele8_Mass55_v2/I");
  mytree_->Branch("trig_isHLT_Ele30WP60_SC4_Mass55_v3",&trig_isHLT_Ele30WP60_SC4_Mass55_v3,"trig_isHLT_Ele30WP60_SC4_Mass55_v3/I");
  mytree_->Branch("trig_isHLT_Ele27_WPLoose_Gsf_v1",&trig_isHLT_Ele27_WPLoose_Gsf_v1,"trig_isHLT_Ele27_WPLoose_Gsf_v1/I");
  
  //
  mytree_->Branch("trig_isEleHLTpath",  &_trig_isEleHLTpath,  "trig_isEleHLTpath/I");
  mytree_->Branch("trig_isMuonHLTpath", &_trig_isMuonHLTpath, "trig_isMuonHLTpath/I");
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
  //
  mytree_->Branch("trig_HLT_N",      &_trig_HLT_N,     "trig_HLT_N/I");
  mytree_->Branch("trig_HLT_eta",    &_trig_HLT_eta,   "trig_HLT_eta[20]/D");
  mytree_->Branch("trig_HLT_phi",    &_trig_HLT_phi,   "trig_HLT_phi[20]/D");
  mytree_->Branch("trig_HLT_energy", &_trig_HLT_energy,"trig_HLT_energy[20]/D");
  mytree_->Branch("trig_HLT_pt",     &_trig_HLT_pt,    "trig_HLT_pt[20]/D");
  mytree_->Branch("trig_HLT_name",   &_trig_HLT_name,   "trig_HLT_name[20]/I");
	
  // Beam Spot
  mytree_->Branch("BS_x",&BS_x,"BS_x/D");
  mytree_->Branch("BS_y",&BS_y,"BS_y/D");
  mytree_->Branch("BS_z",&BS_z,"BS_z/D");
  mytree_->Branch("BS_dz",&BS_dz,"BS_dz/D");
  mytree_->Branch("BS_dxdz",&BS_dxdz,"BS_dxdz/D");
  mytree_->Branch("BS_dydz",&BS_dydz,"BS_dydz/D");
  mytree_->Branch("BS_bw_x",&BS_bw_x,"BS_bw_x/D");
  mytree_->Branch("BS_bw_y",&BS_bw_y,"BS_bw_y/D");
	
  // MC Properties
  mytree_->Branch("MC_pthat",&_MC_pthat,"MC_pthat/D");
  mytree_->Branch("MC_flavor",&_MC_flavor,"MC_flavor[2]/I");
	
  // MC Truth Matching
  mytree_->Branch("ele_MC_chosenEle_PoP_px",ele_MC_chosenEle_PoP_px,"ele_MC_chosenEle_PoP_px[10]/D");
  mytree_->Branch("ele_MC_chosenEle_PoP_py",ele_MC_chosenEle_PoP_py,"ele_MC_chosenEle_PoP_py[10]/D");
  mytree_->Branch("ele_MC_chosenEle_PoP_pz",ele_MC_chosenEle_PoP_pz,"ele_MC_chosenEle_PoP_pz[10]/D");
  mytree_->Branch("ele_MC_chosenEle_PoP_e",ele_MC_chosenEle_PoP_e,"ele_MC_chosenEle_PoP_e[10]/D");
  mytree_->Branch("ele_MC_chosenPho_PoP_px",ele_MC_chosenPho_PoP_px,"ele_MC_chosenPho_PoP_px[10]/D");
  mytree_->Branch("ele_MC_chosenPho_PoP_py",ele_MC_chosenPho_PoP_py,"ele_MC_chosenPho_PoP_py[10]/D");
  mytree_->Branch("ele_MC_chosenPho_PoP_pz",ele_MC_chosenPho_PoP_pz,"ele_MC_chosenPho_PoP_pz[10]/D");
  mytree_->Branch("ele_MC_chosenPho_PoP_e",ele_MC_chosenPho_PoP_e,"ele_MC_chosenPho_PoP_e[10]/D");
  mytree_->Branch("ele_MC_chosenHad_PoP_px",ele_MC_chosenHad_PoP_px,"ele_MC_chosenHad_PoP_px[10]/D");
  mytree_->Branch("ele_MC_chosenHad_PoP_py",ele_MC_chosenHad_PoP_py,"ele_MC_chosenHad_PoP_py[10]/D");
  mytree_->Branch("ele_MC_chosenHad_PoP_pz",ele_MC_chosenHad_PoP_pz,"ele_MC_chosenHad_PoP_pz[10]/D");
  mytree_->Branch("ele_MC_chosenHad_PoP_e",ele_MC_chosenHad_PoP_e,"ele_MC_chosenHad_PoP_e[10]/D");
  mytree_->Branch("ele_MC_closest_DR_px",ele_MC_closest_DR_px,"ele_MC_closest_DR_px[10]/D");
  mytree_->Branch("ele_MC_closest_DR_py",ele_MC_closest_DR_py,"ele_MC_closest_DR_py[10]/D");
  mytree_->Branch("ele_MC_closest_DR_pz",ele_MC_closest_DR_pz,"ele_MC_closest_DR_pz[10]/D");
  mytree_->Branch("ele_MC_closest_DR_e",ele_MC_closest_DR_e,"ele_MC_closest_DR_e[10]/D");
	
  mytree_->Branch("ele_N",&ele_N,"ele_N/I");
  mytree_->Branch("ele_pT",ele_pT,"ele_pT[10]/D");
  mytree_->Branch("ele_echarge",ele_echarge,"ele_echarge[10]/I");
  mytree_->Branch("ele_he",ele_he,"ele_he[10]/D");
  mytree_->Branch("ele_pin_mode",ele_pin_mode,"ele_pin_mode[10]/D");
  mytree_->Branch("ele_pout_mode",ele_pout_mode,"ele_pout_mode[10]/D");
  mytree_->Branch("ele_pin_mean",ele_pin_mean,"ele_pin_mean[10]/D");
  mytree_->Branch("ele_pout_mean",ele_pout_mean,"ele_pout_mean[10]/D");
  mytree_->Branch("ele_pTin_mode",ele_pTin_mode,"ele_pTin_mode[10]/D");
  mytree_->Branch("ele_pTout_mode",ele_pTout_mode,"ele_pTout_mode[10]/D");
  mytree_->Branch("ele_pTin_mean",ele_pTin_mean,"ele_pTin_mean[10]/D");
  mytree_->Branch("ele_pTout_mean",ele_pTout_mean,"ele_pTout_mean[10]/D");
  mytree_->Branch("ele_calo_energy",ele_calo_energy,"ele_calo_energy[10]/D");
  mytree_->Branch("ele_sclRawE",ele_sclRawE,"els_sclRawE[10]/D");
  mytree_->Branch("ele_sclEpresh",ele_sclEpresh,"els_sclEpresh[10]/D");
  mytree_->Branch("ele_sclE",ele_sclE,"ele_sclE[10]/D");
  mytree_->Branch("ele_sclEt",ele_sclEt,"ele_sclEt[10]/D");
  mytree_->Branch("ele_sclEta",ele_sclEta,"ele_sclEta[10]/D");
  mytree_->Branch("ele_sclPhi",ele_sclPhi,"ele_sclPhi[10]/D");
  mytree_->Branch("ele_sclX",ele_sclX,"ele_sclX[10]/D");
  mytree_->Branch("ele_sclY",ele_sclY,"ele_sclY[10]/D");
  mytree_->Branch("ele_sclZ",ele_sclZ,"ele_sclZ[10]/D");
  mytree_->Branch("ele_sclErr",ele_sclErr,"ele_sclErr[10]/D");
  mytree_->Branch("ele_sclErr_pos",ele_sclErr_pos,"ele_sclErr_pos[10]/D");
  mytree_->Branch("ele_sclErr_neg",ele_sclErr_neg,"ele_sclErr_neg[10]/D");
  mytree_->Branch("ele_trErr",ele_trErr,"ele_trErr[10]/D");
  mytree_->Branch("ele_momErr",ele_momErr,"ele_momErr[10]/D");
  mytree_->Branch("ele_newmom",ele_newmom,"ele_newmom[10]/D");
  mytree_->Branch("ele_newmomErr",ele_newmomErr,"ele_newmomErr[10]/D");
  mytree_->Branch("ele_tr_atcaloX",ele_tr_atcaloX,"ele_tr_atcaloX[10]/D");
  mytree_->Branch("ele_tr_atcaloY",ele_tr_atcaloY,"ele_tr_atcaloY[10]/D");
  mytree_->Branch("ele_tr_atcaloZ",ele_tr_atcaloZ,"ele_tr_atcaloZ[10]/D");
  mytree_->Branch("ele_firsthit_X",ele_firsthit_X,"ele_firsthit_X[10]/D");
  mytree_->Branch("ele_firsthit_Y",ele_firsthit_Y,"ele_firsthit_Y[10]/D");
  mytree_->Branch("ele_firsthit_Z",ele_firsthit_Z,"ele_firsthit_Z[10]/D");
	
  // NEW H/E
  mytree_->Branch("ele_he_00615_0",  _ele_he_00615_0 ,"ele_he_00615_0[10]/D");
  mytree_->Branch("ele_he_005_0",  _ele_he_005_0 ,"ele_he_005_0[10]/D");
  mytree_->Branch("ele_he_005_1",  _ele_he_005_1,"ele_he_005_1[10]/D");
  mytree_->Branch("ele_he_005_15", _ele_he_005_15,"ele_he_005_15[10]/D");
  mytree_->Branch("ele_he_01_0",  _ele_he_01_0,"ele_he_01_0[10]/D");
  mytree_->Branch("ele_he_01_1",  _ele_he_01_1,"ele_he_01_1[10]/D");
  mytree_->Branch("ele_he_01_15", _ele_he_01_15,"ele_he_01_15[10]/D");
  mytree_->Branch("ele_he_015_1", _ele_he_015_1,"ele_he_015_1[10]/D");
  mytree_->Branch("ele_he_015_15",_ele_he_015_15 ,"ele_he_015_15[10]/D");

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
  mytree_->Branch("ele_sigmaietaieta",ele_sigmaietaieta,"ele_sigmaietaieta[10]/D");
  mytree_->Branch("ele_sigmaetaeta",ele_sigmaetaeta,"ele_sigmaetaeta[10]/D");
  mytree_->Branch("ele_e15",ele_e15,"ele_e15[10]/D");
  mytree_->Branch("ele_e25max",ele_e25max,"ele_e25max[10]/D");
  mytree_->Branch("ele_e55",ele_e55,"ele_e55[10]/D");
  mytree_->Branch("ele_e1",ele_e1,"ele_e1[10]/D");
  mytree_->Branch("ele_e33",ele_e33,"ele_e33[10]/D");
  mytree_->Branch("ele_e2overe9",ele_e2overe9,"ele_e2overe9[10]/D");
  mytree_->Branch("ele_fbrem",ele_fbrem,"ele_fbrem[10]/D");
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
  mytree_->Branch("ele_dxy",ele_dxy,"ele_dxy[10]/D");
  mytree_->Branch("ele_dxyB",ele_dxyB,"ele_dxyB[10]/D");
  mytree_->Branch("ele_dz",ele_dz,"ele_dz[10]/D");
  mytree_->Branch("ele_dzB",ele_dzB,"ele_dzB[10]/D");
  mytree_->Branch("ele_dsz",ele_dsz,"ele_dsz[10]/D");
  mytree_->Branch("ele_dszB",ele_dszB,"ele_dszB[10]/D");
	
  mytree_->Branch("ele_dzPV",ele_dzPV,"ele_dzPV[10]/D");
  mytree_->Branch("ele_dzPV_error",ele_dzPV_error,"ele_dzPV_error[10]/D");
  mytree_->Branch("ele_dxyPV",ele_dxyPV,"ele_dxyPV[10]/D");
  mytree_->Branch("ele_dxyPV_error",ele_dxyPV_error,"ele_dxyPV_error[10]/D");
  mytree_->Branch("ele_dszPV",ele_dszPV,"ele_dszPV[10]/D");
  mytree_->Branch("ele_dszPV_error",ele_dszPV_error,"ele_dszPV_error[10]/D");

  mytree_->Branch("ele_VetoIdDecisions",ele_VetoIdDecisions,"ele_VetoIdDecisions[10]/O");
  mytree_->Branch("ele_LooseIdDecisions",ele_LooseIdDecisions,"ele_LooseIdDecisions[10]/O");
  mytree_->Branch("ele_MediumIdDecisions",ele_MediumIdDecisions,"ele_MediumIdDecisions[10]/O");
  mytree_->Branch("ele_TightIdDecisions",ele_TightIdDecisions,"ele_TightIdDecisions[10]/O");
	
  mytree_->Branch("ele_severityLevelSeed",ele_severityLevelSeed,"ele_severityLevelSeed[10]/I");
  mytree_->Branch("ele_severityLevelClusters",ele_severityLevelClusters,"ele_severityLevelClusters[10]/I");
  mytree_->Branch("ele_outOfTimeSeed",ele_outOfTimeSeed,"ele_outOfTimeSeed[10]/I");
  mytree_->Branch("ele_outOfTimeClusters",ele_outOfTimeClusters,"ele_outOfTimeClusters[10]/I");
	
  //Conversion Removal
  mytree_->Branch("ele_isConversion",ele_isConversion,"ele_isConversion[10]/I");
  mytree_->Branch("ele_convFound",ele_convFound,"ele_convFound[10]/I");
  mytree_->Branch("ele_conv_dist",&ele_conv_dist,"ele_conv_dist[10]/D");
  mytree_->Branch("ele_conv_dcot",&ele_conv_dcot,"ele_conv_dcot[10]/D");
	
  mytree_->Branch("ele_track_x",ele_track_x,"ele_track_x[10]/D");
  mytree_->Branch("ele_track_y",ele_track_y,"ele_track_y[10]/D");
  mytree_->Branch("ele_track_z",ele_track_z,"ele_track_z[10]/D");
	
  mytree_->Branch("ele_vertex_x",ele_vertex_x,"ele_vertex_x[10]/D");
  mytree_->Branch("ele_vertex_y",ele_vertex_y,"ele_vertex_y[10]/D");
  mytree_->Branch("ele_vertex_z",ele_vertex_z,"ele_vertex_z[10]/D"); 
  mytree_->Branch("ele_tkSumPt_dr03",ele_tkSumPt_dr03,"ele_tkSumPt_dr03[10]/D"); 
  mytree_->Branch("ele_ecalRecHitSumEt_dr03",ele_ecalRecHitSumEt_dr03,"ele_ecalRecHitSumEt_dr03[10]/D"); 
  mytree_->Branch("ele_hcalDepth1TowerSumEt_dr03",ele_hcalDepth1TowerSumEt_dr03,"ele_hcalDepth1TowerSumEt_dr03[10]/D"); 
  mytree_->Branch("ele_hcalDepth2TowerSumEt_dr03",ele_hcalDepth2TowerSumEt_dr03,"ele_hcalDepth2TowerSumEt_dr03[10]/D"); 
  mytree_->Branch("ele_hcalDepth1plus2TowerSumEt_00615dr03",ele_hcalDepth1plus2TowerSumEt_00615dr03,"ele_hcalDepth1plus2TowerSumEt_00615dr03[10]/D");
  mytree_->Branch("ele_hcalDepth1plus2TowerSumEt_005dr03",ele_hcalDepth1plus2TowerSumEt_005dr03,"ele_hcalDepth1plus2TowerSumEt_005dr03[10]/D");
  mytree_->Branch("ele_hcalDepth1plus2TowerSumEt_0dr03",ele_hcalDepth1plus2TowerSumEt_0dr03,"ele_hcalDepth1plus2TowerSumEt_0dr03[10]/D");

  mytree_->Branch("ele_tkSumPt_dr04",ele_tkSumPt_dr04,"ele_tkSumPt_dr04[10]/D"); 
  mytree_->Branch("ele_ecalRecHitSumEt_dr04",ele_ecalRecHitSumEt_dr04,"ele_ecalRecHitSumEt_dr04[10]/D"); 
  mytree_->Branch("ele_hcalDepth1TowerSumEt_dr04",ele_hcalDepth1TowerSumEt_dr04,"ele_hcalDepth1TowerSumEt_dr04[10]/D"); 
  mytree_->Branch("ele_hcalDepth2TowerSumEt_dr04",ele_hcalDepth2TowerSumEt_dr04,"ele_hcalDepth2TowerSumEt_dr04[10]/D"); 
  mytree_->Branch("ele_tkSumPtTdrHzz_dr025",ele_tkSumPtTdrHzz_dr025,"ele_tkSumPtTdrHzz_dr025[10]/D"); 
  mytree_->Branch("ele_tkSumPtoPtTdrHzz_dr025",ele_tkSumPtoPtTdrHzz_dr025,"ele_tkSumPtoPtTdrHzz_dr025[10]/D"); 	
  mytree_->Branch("ele_hcalSumEtTdrHzz_dr02",ele_hcalSumEtTdrHzz_dr02,"ele_hcalSumEtTdrHzz_dr02[10]/D"); 
  mytree_->Branch("ele_hcalSumEtoPtTdrHzz_dr02",ele_hcalSumEtoPtTdrHzz_dr02,"ele_hcalSumEtoPtTdrHzz_dr02[10]/D"); 
  mytree_->Branch("ele_tkSumPtEg4Hzz_dr03",ele_tkSumPtEg4Hzz_dr03,"ele_tkSumPtEg4Hzz_dr03[10]/D"); 
  mytree_->Branch("ele_tkSumPtoPtEg4Hzz_dr03",ele_tkSumPtoPtEg4Hzz_dr03,"ele_tkSumPtoPtEg4Hzz_dr03[10]/D"); 
  mytree_->Branch("ele_ecalSumEtEg4Hzz_dr03",ele_ecalSumEtEg4Hzz_dr03,"ele_ecalSumEtEg4Hzz_dr03[10]/D"); 
  mytree_->Branch("ele_ecalSumEtoPtEg4Hzz_dr03",ele_ecalSumEtoPtEg4Hzz_dr03,"ele_ecalSumEtoPtEg4Hzz_dr03[10]/D"); 
  mytree_->Branch("ele_hcalSumEtEg4Hzz_dr04",ele_hcalSumEtEg4Hzz_dr04,"ele_hcalSumEtEg4Hzz_dr04[10]/D"); 
  mytree_->Branch("ele_hcalSumEtoPtEg4Hzz_dr04",ele_hcalSumEtoPtEg4Hzz_dr04,"ele_hcalSumEtoPtEg4Hzz_dr04[10]/D");
  mytree_->Branch("ele_ambiguousGsfTracks",    ele_ambiguousGsfTracks,    "ele_ambiguousGsfTracks[10]/I");
  mytree_->Branch("ele_ambiguousGsfTracksdxy", ele_ambiguousGsfTracksdxy, "ele_ambiguousGsfTracksdxy[10][5]/D");
  mytree_->Branch("ele_ambiguousGsfTracksdz",  ele_ambiguousGsfTracksdz,  "ele_ambiguousGsfTracksdz[10][5]/D");
  mytree_->Branch("ele_ambiguousGsfTracksdxyB",ele_ambiguousGsfTracksdxyB,"ele_ambiguousGsfTracksdxyB[10][5]/D");
  mytree_->Branch("ele_ambiguousGsfTracksdzB", ele_ambiguousGsfTracksdzB, "ele_ambiguousGsfTracksdzB[10][5]/D");
  mytree_->Branch("ele_seedSubdet1",ele_seedSubdet1,"ele_seedSubdet1[10]/I");
  mytree_->Branch("ele_seedDphi1Pos",ele_seedDphi1Pos,"ele_seedDphi1Pos[10]/D");
  mytree_->Branch("el_eseedDrz1Pos",ele_seedDrz1Pos,"ele_seedDrz1Pos[10]/D");
  mytree_->Branch("ele_seedDphi1Neg",ele_seedDphi1Neg,"ele_seedDphi1Neg[10]/D");
  mytree_->Branch("el_eseedDrz1Neg",ele_seedDrz1Neg,"ele_seedDrz1Neg[10]/D");
  mytree_->Branch("ele_seedSubdet2",ele_seedSubdet2,"ele_seedSubdet2[10]/I");
  mytree_->Branch("ele_seedDphi2Pos",ele_seedDphi2Pos,"ele_seedDphi2Pos[10]/D");
  mytree_->Branch("ele_seedDrz2Pos",ele_seedDrz2Pos,"ele_seedDrz2Pos[10]/D");
  mytree_->Branch("ele_seedDphi2Neg",ele_seedDphi2Neg,"ele_seedDphi2Neg[10]/D");
  mytree_->Branch("ele_seedDrz2Neg",ele_seedDrz2Neg,"ele_seedDrz2Neg[10]/D");
  mytree_->Branch("ele_isMCEle",ele_isMCEle,"ele_isMCEle[10]/I");
  mytree_->Branch("ele_isMCPhoton",ele_isMCPhoton,"ele_isMCPhoton[10]/I");
  mytree_->Branch("ele_isMCHadron",ele_isMCHadron,"ele_isMCHadron[10]/I");
  mytree_->Branch("ele_isSIM",ele_isSIM,"ele_isSIM[10]/I");
  mytree_->Branch("ele_isSIMEle",ele_isSIMEle,"ele_isSIMEle[10]/I");
  mytree_->Branch("ele_idPDGMatch",ele_idPDGMatch,"ele_idPDGMatch[10]/I");
  mytree_->Branch("ele_idPDGmother_MCEle",ele_idPDGmother_MCEle,"ele_idPDGmother_MCEle[10]/I");
  mytree_->Branch("ele_idPDGMatchSim",ele_idPDGMatchSim,"ele_idPDGMatchSim[10]/I");
	
  mytree_->Branch("ele_nSeed", &ele_nSeed, "ele_nSeed/I");	
  mytree_->Branch("ele_SeedIsEcalDriven",ele_SeedIsEcalDriven,"ele_SeedIsEcalDriven[100]/I");
  mytree_->Branch("ele_SeedIsTrackerDriven",ele_SeedIsTrackerDriven,"ele_SeedIsTrackerDriven[100]/I");
	
  mytree_->Branch("ele_SeedSubdet2",ele_SeedSubdet2,"ele_SeedSubdet2[100]/I");
  mytree_->Branch("ele_SeedDphi2Pos",ele_SeedDphi2Pos,"ele_SeedDphi2Pos[100]/D");
  mytree_->Branch("ele_SeedDrz2Pos",ele_SeedDrz2Pos,"ele_SeedDrz2Pos[100]/D");
  mytree_->Branch("ele_SeedDphi2Neg",ele_SeedDphi2Neg,"ele_SeedDphi2Neg[100]/D");
  mytree_->Branch("ele_SeedDrz2Neg",ele_SeedDrz2Neg,"ele_SeedDrz2Neg[100]/D");
  mytree_->Branch("ele_SeedSubdet1",ele_SeedSubdet1,"ele_SeedSubdet1[100]/I");
  mytree_->Branch("ele_SeedDphi1Pos",ele_SeedDphi1Pos,"ele_SeedDphi1Pos[100]/D");
  mytree_->Branch("ele_SeedDrz1Pos",ele_SeedDrz1Pos,"ele_SeedDrz1Pos[100]/D");
  mytree_->Branch("ele_SeedDphi1Neg",ele_SeedDphi1Neg,"ele_SeedDphi1Neg[100]/D");
  mytree_->Branch("ele_SeedDrz1Neg",ele_SeedDrz1Neg,"ele_SeedDrz1Neg[100]/D");
	
  // For Charge, Clemy's stuff
  mytree_->Branch("ele_expected_inner_hits",ele_expected_inner_hits,"ele_expected_inner_hits[10]/I");
  mytree_->Branch("ele_sclNclus",ele_sclNclus,"ele_sclNclus[10]/I");
	
  mytree_->Branch("ele_chargeGsfSC",ele_chargeGsfSC,"ele_chargeGsfSC[10]/I");
  mytree_->Branch("ele_chargeGsfCtf",ele_chargeGsfCtf,"ele_chargeGsfCtf[10]/I");
  mytree_->Branch("ele_chargeGsfCtfSC",ele_chargeGsfCtfSC,"ele_chargeGsfCtfSC[10]/I");
  mytree_->Branch("ele_chargeDPhiInnEle",ele_chargeDPhiInnEle,"ele_chargeDPhiInnEle[10]/D");
  mytree_->Branch("ele_chargeDPhiInnEleCorr",ele_chargeDPhiInnEleCorr,"ele_chargeDPhiInnEleCorr[10]/D");
  mytree_->Branch("ele_chargeQoverPGsfVtx",ele_chargeQoverPGsfVtx,"ele_chargeQoverPGsfVtx[10]/D");
  mytree_->Branch("ele_chargeQoverPCtf",ele_chargeQoverPCtf,"ele_chargeQoverPCtf[10]/D");
  mytree_->Branch("ele_CtfTrackExists",ele_CtfTrackExists,"ele_CtfTrackExists[10]/I");
	
  // For L1 Trigger, Clemy's stuff
  //modif-alex rct region
  mytree_->Branch("ele_RCTeta",         &_ele_RCTeta,          "ele_RCTeta[10]/I");
  mytree_->Branch("ele_RCTphi",         &_ele_RCTphi,          "ele_RCTphi[10]/I");
  mytree_->Branch("ele_RCTL1iso",       &_ele_RCTL1iso,        "ele_RCTL1iso[10]/I");
  mytree_->Branch("ele_RCTL1noniso",    &_ele_RCTL1noniso,     "ele_RCTL1noniso[10]/I");
  mytree_->Branch("ele_RCTL1iso_M",       &_ele_RCTL1iso_M,        "ele_RCTL1iso_M[10]/I");
  mytree_->Branch("ele_RCTL1noniso_M",    &_ele_RCTL1noniso_M,     "ele_RCTL1noniso_M[10]/I");

  mytree_->Branch("ele_L1Stage2",       &_ele_L1Stage2,        "ele_L1Stage2[10]/I");
  mytree_->Branch("ele_L1Stage2_isoflag",  &_ele_L1Stage2_isoflag,        "ele_L1Stage2_isoflag[10]/I");

  mytree_->Branch("ele_L1Stage2_emul",       &_ele_L1Stage2_emul,        "ele_L1Stage2_emul[10]/I");

  mytree_->Branch("ele_L1Stage2_emul_ieta",         &_ele_L1Stage2_emul_ieta,        "ele_L1Stage2_emul_ieta[10]/I");
  mytree_->Branch("ele_L1Stage2_emul_hwPt",         &_ele_L1Stage2_emul_hwPt,        "ele_L1Stage2_emul_hwPt[10]/I");
  mytree_->Branch("ele_L1Stage2_emul_shape",        &_ele_L1Stage2_emul_shape,        "ele_L1Stage2_emul_shape[10]/I");
  mytree_->Branch("ele_L1Stage2_emul_shapeID",      &_ele_L1Stage2_emul_shapeID,        "ele_L1Stage2_emul_shapeID[10]/I");
  mytree_->Branch("ele_L1Stage2_emul_target",       &_ele_L1Stage2_emul_target,        "ele_L1Stage2_emul_targer[10]/D");
  mytree_->Branch("ele_L1Stage2_emul_hwIsoSum",     &_ele_L1Stage2_emul_hwIsoSum,        "ele_L1Stage2_emul_hwIsoSum[10]/I");
  mytree_->Branch("ele_L1Stage2_emul_nTT",          &_ele_L1Stage2_emul_nTT,        "ele_L1Stage2_emul_nTT[10]/I");
  mytree_->Branch("ele_L1Stage2_emul_hOverERatio",  &_ele_L1Stage2_emul_hOverERatio,        "ele_L1Stage2_emul_hOverERatio[10]/I");
  mytree_->Branch("ele_L1Stage2_emul_isoflag",  &_ele_L1Stage2_emul_isoflag,        "ele_L1Stage2_emul_isoflag[10]/I");

  mytree_->Branch("ele_L1Stage1_emul",       &_ele_L1Stage1_emul,        "ele_L1Stage1_emul[10]/I");
  mytree_->Branch("ele_L1Stage1_emul_isoflag",       &_ele_L1Stage1_emul_isoflag,        "ele_L1Stage1_emul_isoflag[10]/I");
  mytree_->Branch("ele_L1Stage1_emul_eta",       &_ele_L1Stage1_emul_eta,        "ele_L1Stage1_emul_eta[10]/D");
  mytree_->Branch("ele_L1Stage1_emul_phi",       &_ele_L1Stage1_emul_phi,        "ele_L1Stage1_emul_phi[10]/D");

  mytree_->Branch("ele_dR_closest_L1Stage2", &_ele_dR_closest_L1Stage2, "ele_dR_closest_L1Stage2[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_eta", &_ele_closestdR_L1Stage2_eta, "ele_closestdR_L1Stage2_eta[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_phi", &_ele_closestdR_L1Stage2_phi, "ele_closestdR_L1Stage2_phi[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_et", &_ele_closestdR_L1Stage2_et, "ele_closestdR_L1Stage2_et[10]/D");

  mytree_->Branch("ele_dR_closest_L1Stage2_emul", &_ele_dR_closest_L1Stage2_emul, "ele_dR_closest_L1Stage2_emul[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_emul_eta", &_ele_closestdR_L1Stage2_emul_eta, "ele_closestdR_L1Stage2_emul_eta[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_emul_phi", &_ele_closestdR_L1Stage2_emul_phi, "ele_closestdR_L1Stage2_emul_phi[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_mp_emul_eta", &_ele_closestdR_L1Stage2_mp_emul_eta, "ele_closestdR_L1Stage2_mp_emul_eta[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_mp_emul_phi", &_ele_closestdR_L1Stage2_mp_emul_phi, "ele_closestdR_L1Stage2_mp_emul_phi[10]/D");
  mytree_->Branch("ele_closestdR_L1Stage2_emul_et", &_ele_closestdR_L1Stage2_emul_et, "ele_closestdR_L1Stage2_emul_et[10]/D");

  mytree_->Branch("ele_TTetaVect",      &_ele_TTetaVect,       "ele_TTetaVect[10][50]/I");
  mytree_->Branch("ele_TTphiVect",      &_ele_TTphiVect,       "ele_TTphiVect[10][50]/I");
  mytree_->Branch("ele_TTetVect",       &_ele_TTetVect,        "ele_TTetVect[10][50]/D");
  mytree_->Branch("ele_TTetaSeed",      &_ele_TTetaSeed,       "ele_TTetaSeed[10]/I");
  mytree_->Branch("ele_TTphiSeed",      &_ele_TTphiSeed,       "ele_TTphiSeed[10]/I");
  mytree_->Branch("ele_TTetSeed",       &_ele_TTetSeed,        "ele_TTetSeed[10]/D");
  mytree_->Branch("ele_RCTetaVect",     &_ele_RCTetaVect,      "ele_RCTetaVect[10][10]/I");
  mytree_->Branch("ele_RCTphiVect",     &_ele_RCTphiVect,      "ele_RCTphiVect[10][10]/I");
  mytree_->Branch("ele_RCTetVect",      &_ele_RCTetVect,       "ele_RCTetVect[10][10]/D");
  mytree_->Branch("ele_RCTL1isoVect",   &_ele_RCTL1isoVect,    "ele_RCTL1isoVect[10][10]/I");
  mytree_->Branch("ele_RCTL1nonisoVect",&_ele_RCTL1nonisoVect, "ele_RCTL1nonisoVect[10][10]/I");
  mytree_->Branch("ele_RCTL1isoVect_M",   &_ele_RCTL1isoVect_M,    "ele_RCTL1isoVect_M[10][10]/I");
  mytree_->Branch("ele_RCTL1nonisoVect_M",&_ele_RCTL1nonisoVect_M, "ele_RCTL1nonisoVect_M[10][10]/I");
  mytree_->Branch("ele_L1Stage1_emul_isoVect",   &_ele_L1Stage1_emul_isoVect,    "ele_L1Stage1_emul_isoVect[10][10]/I");
  mytree_->Branch("ele_L1Stage1_emul_nonisoVect",&_ele_L1Stage1_emul_nonisoVect, "ele_L1Stage1_emul_nonisoVect[10][10]/I");

	
  //ele TIP/LIP/IP
  mytree_->Branch("ele_Tip",&ele_Tip,"ele_Tip[10]/D");
  mytree_->Branch("ele_Lip",&ele_Lip,"ele_Lip[10]/D");
  mytree_->Branch("ele_STip",&ele_STip,"ele_STip[10]/D");
  mytree_->Branch("ele_SLip",&ele_SLip,"ele_SLip[10]/D");
  mytree_->Branch("ele_TipSignif",&ele_TipSignif,"ele_TipSignif[10]/D");
  mytree_->Branch("ele_LipSignif",&ele_LipSignif,"ele_LipSignif[10]/D");
  mytree_->Branch("ele_Significance3D",&ele_Significance3D,"ele_Significance3D[10]/D");
  mytree_->Branch("ele_Value3D",&ele_Value3D,"ele_Value3D[10]/D");
  mytree_->Branch("ele_Error3D",&ele_Error3D,"ele_Error3D[10]/D");
	

  // fbrem ECAL
  mytree_->Branch("ele_ECAL_fbrem",&ele_ECAL_fbrem,"ele_ECAL_fbrem[10]/D");
  mytree_->Branch("ele_PFcomb",&ele_PFcomb,"ele_PFcomb[10]/D");
  mytree_->Branch("ele_PFcomb_Err",&ele_PFcomb_Err,"ele_PFcomb_Err[10]/D");
  mytree_->Branch("ele_PF_SCenergy",&ele_PF_SCenergy,"ele_PF_SCenergy[10]/D");
  mytree_->Branch("ele_PF_SCenergy_Err",&ele_PF_SCenergy_Err,"ele_PF_SCenergy_Err[10]/D");

  // ele 4V
  m_electrons = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("electrons", "TClonesArray", &m_electrons, 856000,0);
	
  // MET
  mytree_->Branch("met_calo_et",&_met_calo_et,"met_calo_et/D");
  mytree_->Branch("met_calo_px",&_met_calo_px,"met_calo_px/D");
  mytree_->Branch("met_calo_py",&_met_calo_py,"met_calo_py/D");
  mytree_->Branch("met_calo_phi",&_met_calo_phi,"met_calo_phi/D");
  mytree_->Branch("met_calo_set",&_met_calo_set,"met_calo_set/D");
  mytree_->Branch("met_calo_sig",&_met_calo_sig,"met_calo_sig/D");
	
  mytree_->Branch("met_calomu_et",&_met_calomu_et,"met_calomu_et/D");
  mytree_->Branch("met_calomu_px",&_met_calomu_px,"met_calomu_px/D");
  mytree_->Branch("met_calomu_py",&_met_calomu_py,"met_calomu_py/D");
  mytree_->Branch("met_calomu_phi",&_met_calomu_phi,"met_calomu_phi/D");
  mytree_->Branch("met_calomu_set",&_met_calomu_set,"met_calomu_set/D");
  mytree_->Branch("met_calomu_sig",&_met_calomu_sig,"met_calomu_sig/D");
	
  mytree_->Branch("met_tc_et",&_met_tc_et,"met_tc_et/D");
  mytree_->Branch("met_tc_px",&_met_tc_px,"met_tc_px/D");
  mytree_->Branch("met_tc_py",&_met_tc_py,"met_tc_py/D");
  mytree_->Branch("met_tc_phi",&_met_tc_phi,"met_tc_phi/D");
  mytree_->Branch("met_tc_set",&_met_tc_set,"met_tc_set/D");
  mytree_->Branch("met_tc_sig",&_met_tc_sig,"met_tc_sig/D");
	
  mytree_->Branch("met_pf_et",&_met_pf_et,"met_pf_et/D");
  mytree_->Branch("met_pf_px",&_met_pf_px,"met_pf_px/D");
  mytree_->Branch("met_pf_py",&_met_pf_py,"met_pf_py/D");
  mytree_->Branch("met_pf_phi",&_met_pf_phi,"met_pf_phi/D");
  mytree_->Branch("met_pf_set",&_met_pf_set,"met_pf_set/D");
  mytree_->Branch("met_pf_sig",&_met_pf_sig,"met_pf_sig/D");
	
  // Muons
  mytree_->Branch("muons_N",&_muons_N,"muons_N/I");
  m_muons = new TClonesArray ("TLorentzVector");
  mytree_->Branch("muons", "TClonesArray", &m_muons, 256000,0);
  mytree_->Branch("muons_charge",&_muons_charge,"muons_charge[20]/I");
  mytree_->Branch("muons_istracker",&_muons_istracker,"muons_istracker[20]/I");
  mytree_->Branch("muons_isstandalone",&_muons_isstandalone,"muons_isstandalone[20]/I");
  mytree_->Branch("muons_isglobal",&_muons_isglobal,"muons_isglobal[20]/I");
  //
  mytree_->Branch("muons_dxy",&_muons_dxy,"muons_dxy[20]/D");
  mytree_->Branch("muons_dz",&_muons_dz,"muons_dz[20]/D");
  mytree_->Branch("muons_dxyPV",&_muons_dxyPV,"muons_dxyPV[20]/D");
  mytree_->Branch("muons_dzPV",&_muons_dzPV,"muons_dzPV[20]/D");
  mytree_->Branch("muons_normalizedChi2",&_muons_normalizedChi2,"muons_normalizedChi2[20]/D");
  mytree_->Branch("muons_NtrackerHits",&_muons_NtrackerHits,"muons_NtrackerHits[20]/I");
  mytree_->Branch("muons_NpixelHits",&_muons_NpixelHits,"muons_NpixelHits[20]/I");
  mytree_->Branch("muons_NmuonHits",&_muons_NmuonHits,"muons_NmuonHits[20]/I");
  mytree_->Branch("muons_Nmatches",&_muons_Nmatches,"muons_Nmatches[20]/I");
  //
  mytree_->Branch("muons_nTkIsoR03",&_muons_nTkIsoR03,"muons_nTkIsoR03[20]/I");
  mytree_->Branch("muons_nTkIsoR05",&_muons_nTkIsoR05,"muons_nTkIsoR05[20]/I");
  mytree_->Branch("muons_tkIsoR03",&_muons_tkIsoR03,"muons_tkIsoR03[20]/D");
  mytree_->Branch("muons_tkIsoR05",&_muons_tkIsoR05,"muons_tkIsoR05[20]/D");
  mytree_->Branch("muons_emIsoR03",&_muons_emIsoR03,"muons_emIsoR03[20]/D");
  mytree_->Branch("muons_emIsoR05",&_muons_emIsoR05,"muons_emIsoR05[20]/D");
  mytree_->Branch("muons_hadIsoR03",&_muons_hadIsoR03,"muons_hadIsoR03[20]/D");
  mytree_->Branch("muons_hadIsoR05",&_muons_hadIsoR05,"muons_hadIsoR05[20]/D");
	
  //muons TIP/LIP/IP
  mytree_->Branch("muons_Tip",&muons_Tip,"muons_Tip[20]/D");
  mytree_->Branch("muons_Lip",&muons_Lip,"muons_Lip[20]/D");
  mytree_->Branch("muons_STip",&muons_STip,"muons_STip[20]/D");
  mytree_->Branch("muons_SLip",&muons_SLip,"muons_SLip[20]/D");
  mytree_->Branch("muons_TipSignif",&muons_TipSignif,"muons_TipSignif[20]/D");
  mytree_->Branch("muons_LipSignif",&muons_LipSignif,"muons_LipSignif[20]/D");
  mytree_->Branch("muons_Significance3D",&muons_Significance3D,"muons_Significance3D[20]/D");
  mytree_->Branch("muons_Value3D",&muons_Value3D,"muons_Value3D[20]/D");
  mytree_->Branch("muons_Error3D",&muons_Error3D,"muons_Error3D[20]/D");
  //muonID variables for HZZ
  mytree_->Branch("muons_trkDxy",&_muons_trkDxy,"muons_trkDxy[20]/D");
  mytree_->Branch("muons_trkDxyError",&_muons_trkDxyError,"muons_trkDxyError[20]/D");
  mytree_->Branch("muons_trkDxyB",&_muons_trkDxyB,"muons_trkDxyB[20]/D");
  mytree_->Branch("muons_trkDz",&_muons_trkDz,"muons_trkDz[20]/D");
  mytree_->Branch("muons_trkDzError",&_muons_trkDzError,"muons_trkDzError[20]/D");
  mytree_->Branch("muons_trkDzB",&_muons_trkDzB,"muons_trkDzB[20]/D"); 
  mytree_->Branch("muons_trkChi2PerNdof",&_muons_trkChi2PerNdof,"muons_trkChi2PerNdof[20]/D");
  mytree_->Branch("muons_trkCharge",&_muons_trkCharge,"muons_trkCharge[20]/D");
  mytree_->Branch("muons_trkNHits",&_muons_trkNHits,"muons_trkNHits[20]/D");
  mytree_->Branch("muons_trkNPixHits",&_muons_trkNPixHits,"muons_trkNPixHits[20]/D");
  mytree_->Branch("muons_trkmuArbitration",&_muons_trkmuArbitration,"muons_trkmuArbitration[20]/D");
  mytree_->Branch("muons_trkmu2DCompatibilityLoose",&_muons_trkmu2DCompatibilityLoose,"muons_trkmu2DCompatibilityLoose[20]/D");
  mytree_->Branch("muons_trkmu2DCompatibilityTight",&_muons_trkmu2DCompatibilityTight,"muons_trkmu2DCompatibilityTight[20]/D");
  mytree_->Branch("muons_trkmuOneStationLoose",&_muons_trkmuOneStationLoose,"muons_trkmuOneStationLoose[20]/D");
  mytree_->Branch("muons_trkmuOneStationTight",&_muons_trkmuOneStationTight,"muons_trkmuOneStationTight[20]/D");
  mytree_->Branch("muons_trkmuLastStationLoose",&_muons_trkmuLastStationLoose,"muons_trkmuLastStationLoose[20]/D");
  mytree_->Branch("muons_trkmuLastStationTight",&_muons_trkmuLastStationTight,"muons_trkmuLastStationTight[20]/D");
  mytree_->Branch("muons_trkmuOneStationAngLoose",&_muons_trkmuOneStationAngLoose,"muons_trkmuOneStationAngLoose[20]/D");
  mytree_->Branch("muons_trkmuOneStationAngTight",&_muons_trkmuOneStationAngTight,"muons_trkmuOneStationAngTight[20]/D");
  mytree_->Branch("muons_trkmuLastStationAngLoose",&_muons_trkmuLastStationAngLoose,"muons_trkmuLastStationAngLoose[20]/D");
  mytree_->Branch("muons_trkmuLastStationAngTight",&_muons_trkmuLastStationAngTight,"muons_trkmuLastStationAngTight[20]/D");
  mytree_->Branch("muons_trkmuLastStationOptimizedLowPtLoose",&_muons_trkmuLastStationOptimizedLowPtLoose,"muons_trkmuLastStationOptimizedLowPtLoose[20]/D");
  mytree_->Branch("muons_trkmuLastStationOptimizedLowPtTight",&_muons_trkmuLastStationOptimizedLowPtTight,"muons_trkmuLastStationOptimizedLowPtTight[20]/D");
  mytree_->Branch("muons_caloCompatibility",&_muons_caloCompatibility,"muons_caloCompatibility[20]/D");
  mytree_->Branch("muons_segmentCompatibility",&_muons_segmentCompatibility,"muons_segmentCompatibility[20]/D");
  mytree_->Branch("muons_glbmuPromptTight",&_muons_glbmuPromptTight,"muons_glbmuPromptTight[20]/D");	
  mytree_->Branch("muons_hzzIso",&_muons_hzzIso,"muons_hzzIso[20]/D"); 
  mytree_->Branch("muons_hzzIsoTk",&_muons_hzzIsoTk,"muons_hzzIsoTk[20]/D"); 
  mytree_->Branch("muons_hzzIsoEcal",&_muons_hzzIsoEcal,"muons_hzzIsoEcal[20]/D"); 
  mytree_->Branch("muons_hzzIsoHcal",&_muons_hzzIsoHcal,"muons_hzzIsoHcal[20]/D"); 

  // Calo Jets
  _m_jets_calo = new TClonesArray ("TLorentzVector");
  mytree_->Branch("jets_calo_N",&_jets_calo_N,"jets_calo_N/I");
  mytree_->Branch("jets_calo", "TClonesArray", &_m_jets_calo, 256000,0);
	
  // JPT jets
  _m_jets_jpt  = new TClonesArray ("TLorentzVector");
  mytree_->Branch("jets_jpt_N", &_jets_jpt_N, "jets_jpt_N/I");
  mytree_->Branch("jets_jpt",  "TClonesArray", &_m_jets_jpt, 256000,0);
	
  // PF jets
  _m_jets_pf   = new TClonesArray ("TLorentzVector");
  mytree_->Branch("jets_pf_N",  &_jets_pf_N,  "jets_pf_N/I");
  mytree_->Branch ("jets_pf",   "TClonesArray", &_m_jets_pf, 256000,0);
	
  mytree_->Branch ("jets_pf_chargedHadEFrac", &jets_pf_chargedHadEFrac,"jets_pf_chargedHadEFrac[100]/D]");
  mytree_->Branch ("jets_pf_chargedEmEFrac",  &jets_pf_chargedEmEFrac, "jets_pf_chargedEmEFrac[100]/D");
  mytree_->Branch ("jets_pf_chargedMuEFrac",  &jets_pf_chargedMuEFrac, "jets_pf_chargedMuEFrac[100]/D");
	
  mytree_->Branch ("jets_pf_neutralHadEFrac", &jets_pf_neutralHadEFrac, "jets_pf_neutralHadEFrac[100]/D");
  mytree_->Branch ("jets_pf_neutralEmEFrac",  &jets_pf_neutralEmEFrac,  "jets_pf_neutralEmEFrac[100]/D");
  mytree_->Branch ("jets_pf_PhotonEFrac",     &jets_pf_PhotonEFrac,     "jets_pf_PhotonEFrac[100]/D");
	
  mytree_->Branch ("jets_pf_chargedHadMultiplicity", &jets_pf_chargedHadMultiplicity, "jets_pf_chargedHadMultiplicity[100]/I");
  mytree_->Branch ("jets_pf_neutralHadMultiplicity", &jets_pf_neutralHadMultiplicity, "jets_pf_neutralHadMultiplicity[100]/I");
	
  mytree_->Branch ("jets_pf_chargedMultiplicity",    &jets_pf_chargedMultiplicity,    "jets_pf_chargedMultiplicity[100]/I");
  mytree_->Branch ("jets_pf_neutralMultiplicity",    &jets_pf_neutralMultiplicity,    "jets_pf_neutralMultiplicity[100]/I");
	
  mytree_->Branch ("jets_pf_nConstituents",          &jets_pf_nConstituents,          "jets_pf_nConstituents[100]/I");
	

	
  // Generated W,Z's & leptons
  _m_MC_gen_V = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("MC_gen_V", "TClonesArray", &_m_MC_gen_V, 256000,0);
  mytree_->Branch ("MC_gen_V_pdgid",&_MC_gen_V_pdgid, "MC_gen_V_pdgid[10]/D");
  _m_MC_gen_leptons = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("MC_gen_leptons", "TClonesArray", &_m_MC_gen_leptons, 256000,0);
  mytree_->Branch ("MC_gen_leptons_pdgid",&_MC_gen_leptons_pdgid, "MC_gen_leptons_pdgid[10]/D");


  /////////////////////////////////////////////////////////
  ///                 Stage 2 Level 1                   ///
  /////////////////////////////////////////////////////////

  mytree_->Branch ("Stage2_tower_n",&_Stage2_tower_n,"Stage2_tower_n/I");
  mytree_->Branch ("Stage2_tower_hwPt",&_Stage2_tower_hwPt,"Stage2_tower_hwPt[5760]/I");
  mytree_->Branch ("Stage2_tower_hwEta",&_Stage2_tower_hwEta,"Stage2_tower_hwEta[5760]/I");
  mytree_->Branch ("Stage2_tower_hwPhi",&_Stage2_tower_hwPhi,"Stage2_tower_hwPhi[5760]/I");

  mytree_->Branch ("Stage2_mpeg_n",&_Stage2_mpeg_n,"Stage2_mpeg_n/I");
  mytree_->Branch ("Stage2_mpeg_hwPt",&_Stage2_mpeg_hwPt,"Stage2_mpeg_hwPt[12]/I");
  mytree_->Branch ("Stage2_mpeg_hwEta",&_Stage2_mpeg_hwEta,"Stage2_mpeg_hwEta[12]/I");
  mytree_->Branch ("Stage2_mpeg_hwPhi",&_Stage2_mpeg_hwPhi,"Stage2_mpeg_hwPhi[12]/I");

  mytree_->Branch ("Stage2_eg_n",&_Stage2_eg_n,"Stage2_eg_n/I");
  mytree_->Branch ("Stage2_eg_hwPt",&_Stage2_eg_hwPt,"Stage2_eg_hwPt[12]/I");
  mytree_->Branch ("Stage2_eg_hwEta",&_Stage2_eg_hwEta,"Stage2_eg_hwEta[12]/I");
  mytree_->Branch ("Stage2_eg_hwPhi",&_Stage2_eg_hwPhi,"Stage2_eg_hwPhi[12]/I");
  mytree_->Branch ("Stage2_eg_et",&_Stage2_eg_et,"Stage2_eg_et[12]/D");
  mytree_->Branch ("Stage2_eg_eta",&_Stage2_eg_eta,"Stage2_eg_eta[12]/D");
  mytree_->Branch ("Stage2_eg_phi",&_Stage2_eg_phi,"Stage2_eg_phi[12]/D");
  mytree_->Branch ("Stage2_eg_isoflag",&_Stage2_eg_isoflag,"Stage2_eg_isoflag[12]/I");

  mytree_->Branch ("Stage2_tower_emul_n",&_Stage2_tower_emul_n,"Stage2_tower_emul_n/I");
  mytree_->Branch ("Stage2_tower_emul_hwPt",&_Stage2_tower_emul_hwPt,"Stage2_tower_emul_hwPt[5760]/I");
  mytree_->Branch ("Stage2_tower_emul_hwEtEm",&_Stage2_tower_emul_hwEtEm,"Stage2_tower_emul_hwEtEm[5760]/I");
  mytree_->Branch ("Stage2_tower_emul_hwEtHad",&_Stage2_tower_emul_hwEtHad,"Stage2_tower_emul_hwEtHad[5760]/I");
  mytree_->Branch ("Stage2_tower_emul_hwEta",&_Stage2_tower_emul_hwEta,"Stage2_tower_emul_hwEta[5760]/I");
  mytree_->Branch ("Stage2_tower_emul_hwPhi",&_Stage2_tower_emul_hwPhi,"Stage2_tower_emul_hwPhi[5760]/I");

  mytree_->Branch ("Stage2_mpeg_emul_n",&_Stage2_mpeg_emul_n,"Stage2_mpeg_emul_n/I");
  mytree_->Branch ("Stage2_mpeg_emul_hwPt",&_Stage2_mpeg_emul_hwPt,"Stage2_mpeg_emul_hwPt[12]/I");
  mytree_->Branch ("Stage2_mpeg_emul_hwEta",&_Stage2_mpeg_emul_hwEta,"Stage2_mpeg_emul_hwEta[12]/I");
  mytree_->Branch ("Stage2_mpeg_emul_hwPhi",&_Stage2_mpeg_emul_hwPhi,"Stage2_mpeg_emul_hwPhi[12]/I");
  mytree_->Branch ("Stage2_mpeg_emul_et",&_Stage2_mpeg_emul_et,"Stage2_mpeg_emul_et[12]/D");
  mytree_->Branch ("Stage2_mpeg_emul_eta",&_Stage2_mpeg_emul_eta,"Stage2_mpeg_emul_eta[12]/D");
  mytree_->Branch ("Stage2_mpeg_emul_phi",&_Stage2_mpeg_emul_phi,"Stage2_mpeg_emul_phi[12]/D");
  mytree_->Branch ("Stage2_mpeg_emul_shape",&_Stage2_mpeg_emul_shape,"Stage2_mpeg_emul_shape[12]/I");
  mytree_->Branch ("Stage2_mpeg_emul_shapeID",&_Stage2_mpeg_emul_shapeID,"Stage2_mpeg_emul_shapeID[12]/I");

  mytree_->Branch("Stage2_mpeg_emul_hwIsoSum", &_Stage2_mpeg_emul_hwIsoSum,"Stage2_mpeg_emul_hwIsoSum[12]/I");
  mytree_->Branch("Stage2_mpeg_emul_nTT", &_Stage2_mpeg_emul_nTT,"Stage2_mpeg_emul_nTT[12]/I");
  mytree_->Branch("Stage2_mpeg_emul_hOverERatio", &_Stage2_mpeg_emul_hOverERatio,"Stage2_mpeg_emul_hOverERatio[12]/I");
  mytree_->Branch("Stage2_mpeg_emul_isoflag", &_Stage2_mpeg_emul_isoflag,"Stage2_mpeg_emul_isoflag[12]/I");


  mytree_->Branch ("Stage2_eg_emul_n",&_Stage2_eg_emul_n,"Stage2_eg_emul_n/I");
  mytree_->Branch ("Stage2_eg_emul_hwPt",&_Stage2_eg_emul_hwPt,"Stage2_eg_emul_hwPt[12]/I");
  mytree_->Branch ("Stage2_eg_emul_hwEta",&_Stage2_eg_emul_hwEta,"Stage2_eg_emul_hwEta[12]/I");
  mytree_->Branch ("Stage2_eg_emul_hwPhi",&_Stage2_eg_emul_hwPhi,"Stage2_eg_emul_hwPhi[12]/I");
  mytree_->Branch ("Stage2_eg_emul_et",&_Stage2_eg_emul_et,"Stage2_eg_emul_et[12]/D");
  mytree_->Branch ("Stage2_eg_emul_eta",&_Stage2_eg_emul_eta,"Stage2_eg_emul_eta[12]/D");
  mytree_->Branch ("Stage2_eg_emul_phi",&_Stage2_eg_emul_phi,"Stage2_eg_emul_phi[12]/D");

  mytree_->Branch ("Stage1_eg_emul_n",&_Stage1_eg_emul_n,"Stage1_eg_emul_n/I");
  mytree_->Branch ("Stage1_eg_emul_hwPt",&_Stage1_eg_emul_hwPt,"Stage1_eg_emul_hwPt[12]/I");
  mytree_->Branch ("Stage1_eg_emul_hwEta",&_Stage1_eg_emul_hwEta,"Stage1_eg_emul_hwEta[12]/I");
  mytree_->Branch ("Stage1_eg_emul_hwPhi",&_Stage1_eg_emul_hwPhi,"Stage1_eg_emul_hwPhi[12]/I");
  mytree_->Branch("Stage1_eg_emul_isoflag", &_Stage1_eg_emul_isoflag,"Stage1_eg_emul_isoflag[12]/I");
  mytree_->Branch ("Stage1_eg_emul_et",&_Stage1_eg_emul_et,"Stage1_eg_emul_et[12]/D");
  mytree_->Branch ("Stage1_eg_emul_eta",&_Stage1_eg_emul_eta,"Stage1_eg_emul_eta[12]/D");
  mytree_->Branch ("Stage1_eg_emul_phi",&_Stage1_eg_emul_phi,"Stage1_eg_emul_phi[12]/D");


}

// ====================================================================================
SimpleNtpleCustom::~SimpleNtpleCustom()
// ====================================================================================
{
  delete m_electrons ;
  delete m_muons;
  delete _m_jets_calo;
  delete _m_jets_jpt;
  delete _m_jets_pf;
	
  if(type_ == "MC") {
    delete _m_MC_gen_V;
    delete _m_MC_gen_leptons;
  } // if MC
}

// ====================================================================================
void SimpleNtpleCustom::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  
  /*if(PrintDebug_) std::cout << "Update Field" << std::endl;
  // Clemy's Stuff for Charge
  bool updateField(false);
  if (cacheIDMagField_!=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier()){
    updateField = true;
    cacheIDMagField_=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier();
    iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  }
	
  if(PrintDebug_) std::cout << "Update Geo" << std::endl;
  bool updateGeometry(false);
  if (cacheIDTDGeom_!=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier()){
    updateGeometry = true;
    cacheIDTDGeom_=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier();
    iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
  }
  
  if(PrintDebug_) std::cout << "MutliTrajectory" << std::endl;
  if(updateField || updateGeometry){
    mtsTransform_ = new MultiTrajectoryStateTransform(trackerHandle_.product(),theMagField.product());
  }
	
  */
	
  // Tree Maker
  Init();
  if (funcbase_) funcbase_->init(iSetup);
  if(PrintDebug_) std::cout << "FillEvent" << std::endl;
  FillEvent (iEvent, iSetup);
	
  //
  if(PrintDebug_) std::cout << "FillTrigger (iEvent, iSetup);" << std::endl;
  FillTrigger (iEvent, iSetup);

  m_electrons -> Clear() ;
  m_muons -> Clear() ;

  if(type_ == "MC") {
    _m_MC_gen_V->Clear();
    _m_MC_gen_leptons->Clear();
  } // if MC

  if(PrintDebug_) std::cout << "FillEle (iEvent, iSetup);" << std::endl;
  FillEle (iEvent, iSetup);
  //
	
  mytree_->Fill();

	
} // analyze

// ====================================================================================
void SimpleNtpleCustom::FillEvent (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  nEvent = iEvent.id().event();
  nRun   = iEvent.id().run();
  nLumi  = iEvent.luminosityBlock();

  // -----------------
  // Pile-up
  // -----------------
  if(type_ == "MC") {
    Handle<vector<PileupSummaryInfo> > PupInfo;
    //iEvent.getByLabel(PileupSrc_, PupInfo);
    iEvent.getByToken(PileupSrcToken_, PupInfo);
    for (vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand) {
	    
      _PU_N = cand->getPU_NumInteractions();
      
    } // loop on Pile up
  } // if MC

  if(PrintDebug_) std::cout << "Rho Fast Jet" << std::endl;
  // Rho/FastJet Correction
  Handle<double> rhoHandle, sigmaHandle;
  //iEvent.getByLabel(RhoCorrection_, rhoHandle);
  //iEvent.getByLabel(SigmaRhoCorrection_, sigmaHandle);
  iEvent.getByToken(RhoCorrectionToken_, rhoHandle);
  iEvent.getByToken(SigmaRhoCorrectionToken_, sigmaHandle);
  _PU_rho   = *rhoHandle;
  _PU_sigma = *sigmaHandle;
    

  // -----------------
  // Vertices
  // -----------------
  
  if(PrintDebug_) std::cout << "Vertices" << std::endl;
  Handle<reco::VertexCollection> recoPrimaryVertexCollection;
  //iEvent.getByLabel(VerticesTag_,recoPrimaryVertexCollection);
  iEvent.getByToken(VerticesToken_,recoPrimaryVertexCollection);
  
  const reco::VertexCollection & vertices = *recoPrimaryVertexCollection.product();
    
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  
  int vtx_counter=0;
  _vtx_N = recoPrimaryVertexCollection->size();
  
  // select the primary vertex as the one with higest sum of (pt)^2 of tracks                                                                               
  
  
  if(_vtx_N > 0) {
    reco::VertexCollection::const_iterator firstPV = vertices.begin();
    GlobalPoint local_vertexPosition(firstPV->position().x(),
				     firstPV->position().y(),
				     firstPV->position().z());
      vertexPosition = local_vertexPosition;
    }
    else {
      GlobalPoint local_vertexPosition(0,0,0);
    
      vertexPosition = local_vertexPosition;
    }

    
    if(PrintDebug_) std::cout << "Loop on vertices" << std::endl;
    for( std::vector<reco::Vertex>::const_iterator PV = vertices.begin(); PV != vertices.end(); ++PV){
      if(vtx_counter > 14 ) continue;
      
      _vtx_normalizedChi2[vtx_counter] = PV->normalizedChi2();
      _vtx_ndof[vtx_counter] = PV->ndof();
      _vtx_nTracks[vtx_counter] = PV->tracksSize();
      _vtx_d0[vtx_counter] = PV->position().Rho();
      _vtx_x[vtx_counter] = PV->x();
      _vtx_y[vtx_counter] = PV->y();
      _vtx_z[vtx_counter] = PV->z();
      
      vtx_counter++;
    } // for loop on primary vertices
    
    //if(vtx_counter>14) { _vtx_N = 15; cout << "Number of primary vertices>15, vtx_N set to 15" << endl;}
	
}



// ====================================================================================
void SimpleNtpleCustom::FillTrigger (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // ----------------------------------------------
  //  Get HLT info
  // ----------------------------------------------	
	
  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  //iEvent.getByLabel (HLTTag_,triggerResultsHandle);
  iEvent.getByToken (HLTToken_,triggerResultsHandle);
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle);
	
	
  trig_isUnbiased = 0 ;


  
  // LOOP Over Trigger Results
  strcpy(trig_fired_names,"*");
  for (int iHLT=0; iHLT<static_cast<int>(triggerResultsHandle->size()); iHLT++) {	

    if (triggerResultsHandle->accept (iHLT)) {
      trig_hltInfo[iHLT] = 1;
      if ( strlen(trig_fired_names) <= 4950) {
	const char* c_str();
	string hlt_string = triggerNames.triggerName(iHLT);
	cout<<hlt_string<<endl;
	strcat(trig_fired_names,hlt_string.c_str());
	strcat(trig_fired_names,"*");
      }
    }
    else {
      trig_hltInfo[iHLT] = 0;
    }


		
    if(string(triggerNames.triggerName(iHLT)).find("HLT_Activity_Ecal_SC")!= std::string::npos) {
      if (triggerResultsHandle->accept (iHLT)) trig_isUnbiased = 1 ;
    }
		
    if (triggerNames.triggerName(iHLT) == "HLT_L1SingleEG2" ) {
      if (triggerResultsHandle->accept (iHLT)) 
	trig_isL1SingleEG2 = 1 ;
      else 
	trig_isL1SingleEG2 = 0 ;
    }
    if (triggerNames.triggerName(iHLT) == "HLT_L1SingleEG5" ) {
      if (triggerResultsHandle->accept (iHLT)) 
	trig_isL1SingleEG5 = 1 ;
      else 
	trig_isL1SingleEG5 = 0 ;
    }
    if (triggerNames.triggerName(iHLT) == "HLT_L1SingleEG8" ) {
      if (triggerResultsHandle->accept (iHLT)) 
	trig_isL1SingleEG8 = 1 ;
      else 
	trig_isL1SingleEG8 = 0 ;
    }
    if (triggerNames.triggerName(iHLT) == "HLT_Photon10_L1R" ) {
      if (triggerResultsHandle->accept (iHLT)) 
	trig_isPhoton10 = 1 ;
      else 
	trig_isPhoton10 = 0 ;
    }
    if (triggerNames.triggerName(iHLT) == "HLT_Photon15_L1R" ) {
      if (triggerResultsHandle->accept (iHLT)) 
	trig_isPhoton15 = 1 ;
      else 
	trig_isPhoton15 = 0 ;
    }
    if (triggerNames.triggerName(iHLT) == "HLT_Ele10_LW_L1R" ) {
      if (triggerResultsHandle->accept (iHLT)) 
	trig_isEle10_LW = 1 ;
      else 
	trig_isEle10_LW = 0 ;
    }
    if (triggerNames.triggerName(iHLT) == "HLT_Ele15_LW_L1R" ) {
      if (triggerResultsHandle->accept (iHLT)) 
	trig_isEle15_LW = 1 ;
      else 
	trig_isEle15_LW = 0 ;
    } // if HLT Ele

    if (triggerNames.triggerName(iHLT) == "HLT_Ele30WP60_Ele8_Mass55_v2" ) {
      if (triggerResultsHandle->accept (iHLT)){ 
	trig_isHLT_Ele30WP60_Ele8_Mass55_v2 = 1 ;
      }
      else 
	trig_isHLT_Ele30WP60_Ele8_Mass55_v2 = 0 ;
    }
    
    if (triggerNames.triggerName(iHLT) == "HLT_Ele30WP60_SC4_Mass55_v3" ) {
      if (triggerResultsHandle->accept (iHLT)){
	trig_isHLT_Ele30WP60_SC4_Mass55_v3 = 1 ;
      }
      else 
	trig_isHLT_Ele30WP60_SC4_Mass55_v3 = 0 ;
    }

    if (triggerNames.triggerName(iHLT) == "HLT_Ele27_WPLoose_Gsf_v1" ) {
      if (triggerResultsHandle->accept (iHLT)){
	trig_isHLT_Ele27_WPLoose_Gsf_v1 = 1 ;
      }
      else 
	trig_isHLT_Ele27_WPLoose_Gsf_v1 = 0 ;
    }

		
  } // for loop on trigger results 	


  /////////////////////////// 
  // Get TP data  (Nadir)  //
  ///////////////////////////

  
  
  // ORIGINAL TP
  if( nadGetTP_ ) {
    if(PrintDebug_) cout << "create new ecal_tp pointer" << endl;
    edm::Handle<EcalTrigPrimDigiCollection>* ecal_tp_ = new edm::Handle<EcalTrigPrimDigiCollection> ;

    if(PrintDebug_) cout << "..created. get by label the tp collection" << endl;

    //iEvent.getByLabel(tpCollectionNormal_,*ecal_tp_);
    iEvent.getByToken(tpCollectionNormalToken_,*ecal_tp_);

    if(PrintDebug_) cout << "got it" << endl;

    _trig_tower_N = ecal_tp_->product()->size();
    if(PrintDebug_) {
      cout << "TP Normal collection size=" << ecal_tp_->product()->size() << endl ;
      cout << "is gonna get the TP data" << endl;
    }
  
    for (int i=0 ; i<_trig_tower_N ; i++) {
      if(PrintDebug_) cout << "loop iteration #" << i << endl;
      EcalTriggerPrimitiveDigi d_ = (*(ecal_tp_->product()))[i]; // EcalTriggerPrimitiveDigi d
      if(PrintDebug_) cout << "got the trigger primitive" << endl;
      EcalTrigTowerDetId TPtowid_ = d_.id(); // const EcalTrigTowerDetId TPtowid
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


  
  // ZEROING-BY-HAND TP
  if( nadGetTP_Modif_ ) {
	edm::Handle<EcalTrigPrimDigiCollection>* ecal_tpM_ = new edm::Handle<EcalTrigPrimDigiCollection> ;
	//iEvent.getByLabel(tpCollectionModif_,*ecal_tpM_);    
	iEvent.getByToken(tpCollectionModifToken_,*ecal_tpM_);    

	_trig_tower_N_modif = ecal_tpM_->product()->size(); 
	
	for (int i=0 ; i<_trig_tower_N_modif ; i++) {
	  EcalTriggerPrimitiveDigi dM_ = (*(ecal_tpM_->product()))[i]; // EcalTriggerPrimitiveDigi dM
	  EcalTrigTowerDetId TPtowidM_ = dM_.id(); // EcalTrigTowerDetId
	  _trig_tower_iphi_modif[i] = TPtowidM_.iphi() ;
	  _trig_tower_ieta_modif[i] = TPtowidM_.ieta() ;     
	  _trig_tower_adc_modif[i]  = (dM_[0].raw()&0xff) ;      
	  _trig_tower_sFGVB_modif[i] = dM_[0].sFGVB(); // 0=spike-like / 1=EM-like
	  
	}
  }

  // EMULATOR TPs
  if( nadGetTP_Emul_ ) {
    
    edm::Handle<EcalTrigPrimDigiCollection>* ecal_tpM_ = new edm::Handle<EcalTrigPrimDigiCollection> ;
    //iEvent.getByLabel(tpEmulatorCollection_, *ecal_tpM_);    
    iEvent.getByToken(tpEmulatorCollectionToken_, *ecal_tpM_);    
  
    _trig_tower_N_emul = ecal_tpM_->product()->size();

    for (int i=0 ; i<_trig_tower_N_emul ; i++) {
      EcalTriggerPrimitiveDigi dM_ = (*(ecal_tpM_->product()))[i]; //EcalTriggerPrimitiveDigi
      EcalTrigTowerDetId TPtowidM_ = dM_.id();
      _trig_tower_iphi_emul[i] = TPtowidM_.iphi() ;
      _trig_tower_ieta_emul[i] = TPtowidM_.ieta() ;
    
      bool showit = false;
      for(int j=0 ; j<5 ; j++)
	if( (dM_[j].raw()&0xff) > 0 ) showit = true ;
      showit = false;

      if(showit)
	cout << "TTieta=" << TPtowidM_.ieta() << " TTiphi=" << TPtowidM_.iphi() << " adcEm=" ;

      for (int j=0 ; j<5 ; j++) {
	_trig_tower_adc_emul[i][j] = (dM_[j].raw()&0xff) ;	
	_trig_tower_sFGVB_emul[i][j] = dM_[j].sFGVB(); 
	if(showit)
	  cout << (dM_[j].raw()&0xff) << " " ;
      }
      if(showit)
	cout << endl;
    }
    
  }  

  
  // ----------------------------------
  //  Path from list given in .py file
  // ----------------------------------
  UInt_t trigger_size = triggerResultsHandle->size();
  int passEleTrigger  = 0;
  int passMuonTrigger = 0;
	
  // Electron Triggers
  for(int ipath=0;ipath< (int) HLT_ElePaths_.size();ipath++) {    
    UInt_t trigger_position = triggerNames.triggerIndex(HLT_ElePaths_[ipath]);
    if (trigger_position < trigger_size) passEleTrigger = (int)triggerResultsHandle->accept(trigger_position);
    if (passEleTrigger==1) _trig_isEleHLTpath = 1;
  } // for loop on HLT Elepaths
	
  // Muon Triggers
  for(int ipath=0;ipath< (int) HLT_MuonPaths_.size();ipath++) {    
    UInt_t trigger_position = triggerNames.triggerIndex(HLT_MuonPaths_[ipath]);
    if (trigger_position < trigger_size) passMuonTrigger = (int)triggerResultsHandle->accept(trigger_position);
    if (passMuonTrigger==1) _trig_isMuonHLTpath = 1;
  } // for loop on HLT Muonpaths
	


  if(!aod_) {

    // ----------------------
    //  get L1 EM candidate
    // ----------------------
    
    // --- CURRENT BUNCH CROSSING --- //////////////////////////////////////////////////////////////////

    edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl ;
    edm::Handle< l1extra::L1EmParticleCollection > emIsolColl ;
    edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl_M ;
    edm::Handle< l1extra::L1EmParticleCollection > emIsolColl_M ;  

    /*if( !nadGetL1M_ ) {            
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
      }*/

    iEvent.getByToken(emNonisolCollToken_,emNonisolColl);
    iEvent.getByToken(emIsolCollToken_,emIsolColl);

    //if( nadGetL1M_ ) {  
    //iEvent.getByToken(emNonisolColl_M_Token_,emNonisolColl_M);
    // iEvent.getByToken(emIsolColl_M_Token_,emIsolColl_M);
    //}
      


    ///// STANDARD COLLECTION ALONE
    
    // Isolated candidates
    _trig_L1emIso_N = emIsolColl->size();
    if(PrintDebug_) cout << "N L1 candidate iso : " << _trig_L1emIso_N << endl;
    int counter = 0;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl->begin(); emItr != emIsolColl->end() ;++emItr) {
      // Used by Clemy
      _trig_L1emIso_ieta[counter] = emItr->gctEmCand()->regionId().ieta();
      _trig_L1emIso_iphi[counter] = emItr->gctEmCand()->regionId().iphi();
      _trig_L1emIso_rank[counter] = emItr->gctEmCand()->rank(); // ET in ADC count... 1 ADC count = 0.5 GeV
      // From Trigger twiki
      _trig_L1emIso_eta[counter]    = emItr->eta();
      _trig_L1emIso_phi[counter]    = emItr->phi();
      _trig_L1emIso_energy[counter] = emItr->energy();
      _trig_L1emIso_et[counter]     = emItr->et();
      counter++;
    }
	  
    // Non Isolated candidates
    _trig_L1emNonIso_N = emNonisolColl->size();
    if(PrintDebug_) cout << "N L1 candidate noniso : " << _trig_L1emNonIso_N << endl;	  
    counter = 0;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl->begin(); emItr != emNonisolColl->end() ;++emItr){  
      // Used by Clemy
      _trig_L1emNonIso_ieta[counter] = emItr->gctEmCand()->regionId().ieta();
      _trig_L1emNonIso_iphi[counter] = emItr->gctEmCand()->regionId().iphi();
      _trig_L1emNonIso_rank[counter] = emItr->gctEmCand()->rank(); // ET in ADC count... 1 ADC count = 0.5 GeV
      // From Trigger twiki
      _trig_L1emNonIso_eta[counter]    = emItr->eta();
      _trig_L1emNonIso_phi[counter]    = emItr->phi();
      _trig_L1emNonIso_energy[counter] = emItr->energy();
      _trig_L1emNonIso_et[counter]     = emItr->et();
      counter++;
    } // for loop on Non Iso cand
	  
    /*if( nadGetL1M_ ) {
      _trig_L1emIso_N_M = emIsolColl_M->size();
      if(PrintDebug_) cout << "_trig_L1emIso_N_M =" << _trig_L1emIso_N_M << endl;
      counter=0;
      for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl_M->begin();
           emItr != emIsolColl_M->end() ;++emItr) {
           _trig_L1emIso_ieta_M[counter] = emItr->gctEmCand()->regionId().ieta();
           _trig_L1emIso_iphi_M[counter] = emItr->gctEmCand()->regionId().iphi();
           _trig_L1emIso_rank_M[counter] = emItr->gctEmCand()->rank();
           _trig_L1emIso_eta_M[counter]    = emItr->eta();
           _trig_L1emIso_phi_M[counter]    = emItr->phi();
           _trig_L1emIso_energy_M[counter] = emItr->energy();
           _trig_L1emIso_et_M[counter]     = emItr->et();
           counter++;
      }
 
      
      // Non Isolated candidates
      _trig_L1emNonIso_N_M = emNonisolColl_M->size();
      counter = 0;  
      for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonisolColl_M->begin(); 
	   emItr != emNonisolColl_M->end() ;++emItr){  
	// Used by Clemy
	_trig_L1emNonIso_ieta_M[counter] = emItr->gctEmCand()->regionId().ieta();
	_trig_L1emNonIso_iphi_M[counter] = emItr->gctEmCand()->regionId().iphi();
	_trig_L1emNonIso_rank_M[counter] = emItr->gctEmCand()->rank(); 
	// ET in ADC count... 1 ADC count = 0.5 GeV
	// From Trigger twiki
	_trig_L1emNonIso_eta_M[counter]    = emItr->eta();
	_trig_L1emNonIso_phi_M[counter]    = emItr->phi();
	_trig_L1emNonIso_energy_M[counter] = emItr->energy();
	_trig_L1emNonIso_et_M[counter]     = emItr->et();
	counter++;
      } // for loop on Non Iso cand
      }*/
    ///////////////////////////////////////////////////////////////////////////////////////////////

	  
    // --- PRE- AND POST-FIRING ---
	  
    /*edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel( edm::InputTag(gtRecordCollectionTag_), gtRecord);
    //PRE-FIRING
    const L1GtPsbWord psb = gtRecord->gtPsbWord(0xbb0d, -1);
    //psb.print(cout); 
    std::vector<int> psbel;
    psbel.push_back(psb.aData(4));
    psbel.push_back(psb.aData(5));
    psbel.push_back(psb.bData(4));
    psbel.push_back(psb.bData(5));
    counter = 0;
    std::vector<int>::const_iterator ipsbel;
    for(ipsbel=psbel.begin(); ipsbel!=psbel.end(); ipsbel++) {
      int rank = (*ipsbel)&0x3f; // ET in ADC count... 1 ADC count = 0.5 GeV
      if(rank>0) {
	int iEta = int(((*ipsbel)>>6)&7);
	int sign = ( ((*ipsbel>>9)&1) ? -1. : 1. ); 
	int regionEtaRec;
	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (pre, non-iso)"<<std::endl;
	// Used by Clemy
	_trig_preL1emNonIso_ieta[counter] = regionEtaRec;
	_trig_preL1emNonIso_iphi[counter] = int(((*ipsbel)>>10)&0x1f);
	_trig_preL1emNonIso_rank[counter] = rank;
	counter++;
      }
    }//loop Noniso
    _trig_preL1emNonIso_N = counter;
	  
    psbel.clear();
    psbel.push_back(psb.aData(6));
    psbel.push_back(psb.aData(7));
    psbel.push_back(psb.bData(6));
    psbel.push_back(psb.bData(7));
    counter = 0;
    for(ipsbel=psbel.begin(); ipsbel!=psbel.end(); ipsbel++) {
      int rank = (*ipsbel)&0x3f; // ET in ADC count... 1 ADC count = 0.5 GeV
      if(rank>0) {
	int iEta = int(((*ipsbel)>>6)&7);
	int sign = ( ((*ipsbel>>9)&1) ? -1. : 1. ); 
	int regionEtaRec;
	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (pre, iso)"<<std::endl;
	// Used by Clemy
	_trig_preL1emIso_ieta[counter] = regionEtaRec;
	_trig_preL1emIso_iphi[counter] = int(((*ipsbel)>>10)&0x1f);
	_trig_preL1emIso_rank[counter] = rank;
	counter++;
      }
    }//loop Iso
    _trig_preL1emIso_N = counter;
	  
	  
    //POST-FIRING
    const L1GtPsbWord psb2 = gtRecord->gtPsbWord(0xbb0d, 1);
    std::vector<int> psbel2;
    psbel2.push_back(psb2.aData(4));
    psbel2.push_back(psb2.aData(5));
    psbel2.push_back(psb2.bData(4));
    psbel2.push_back(psb2.bData(5));
    counter = 0;
    std::vector<int>::const_iterator ipsbel2;
    for(ipsbel2=psbel2.begin(); ipsbel2!=psbel2.end(); ipsbel2++) {
      int rank = (*ipsbel2)&0x3f; // ET in ADC count... 1 ADC count = 0.5 GeV
      if(rank>0) {
	int iEta = int(((*ipsbel2)>>6)&7);
	int sign = ( ((*ipsbel2>>9)&1) ? -1. : 1. ); 
	int regionEtaRec;
	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (post, non-iso)"<<std::endl;
	// Used by Clemy
	_trig_postL1emNonIso_ieta[counter] = regionEtaRec;
	_trig_postL1emNonIso_iphi[counter] = int(((*ipsbel2)>>10)&0x1f);
	_trig_postL1emNonIso_rank[counter] = rank;
	counter++;
      }
    }//loop Noniso
    _trig_postL1emNonIso_N = counter;
	  
    psbel2.clear();
    psbel2.push_back(psb2.aData(6));
    psbel2.push_back(psb2.aData(7));
    psbel2.push_back(psb2.bData(6));
    psbel2.push_back(psb2.bData(7));
    counter = 0;
    for(ipsbel2=psbel2.begin(); ipsbel2!=psbel2.end(); ipsbel2++) {
      int rank = (*ipsbel2)&0x3f; // ET in ADC count... 1 ADC count = 0.5 GeV
      if(rank>0) {
	int iEta = int(((*ipsbel2)>>6)&7);
	int sign = ( ((*ipsbel2>>9)&1) ? -1. : 1. ); 
	int regionEtaRec;
	if(sign > 0) regionEtaRec = iEta + 11;
	if(sign < 0) regionEtaRec = 10 - iEta;
	if(sign==0) std::cout<<"WEIRD (post, iso)"<<std::endl;
	// Used by Clemy
	_trig_postL1emIso_ieta[counter] = regionEtaRec;
	_trig_postL1emIso_iphi[counter] = int(((*ipsbel2)>>10)&0x1f);
	_trig_postL1emIso_rank[counter] = rank;
	counter++;
      }
    }//loop Iso
    _trig_postL1emIso_N = counter;*/
	  
  } // if AOD
	
	
  // ----------------------
  //  get HLT EM candidate
  // ----------------------
  edm::Handle<trigger::TriggerEvent> trigEvent;
  //iEvent.getByLabel(triggerEventTag_, trigEvent);
  iEvent.getByToken(triggerEventToken_, trigEvent);
	
  const Int_t N_filter(trigEvent->sizeFilters());
  std::vector<Int_t> ID_filter; 
	
	
  int hlt_counter = 0;

  // Loop on user's Filters
  for(int itrig=0;itrig< (int) HLT_Filters_.size();itrig++) {
		
		
    ID_filter.push_back(trigEvent->filterIndex(HLT_Filters_[itrig])); 
		
    const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());
    if( ID_filter[itrig] <  N_filter) {
      const trigger::Keys& keys( trigEvent->filterKeys(ID_filter[itrig])); 
			
      // Loop on HLT objects
      for ( int hlto = 0; hlto < (int) keys.size(); hlto++ ) {
	if(hlt_counter>19) continue;
				
	trigger::size_type hltf = keys[hlto];
	const trigger::TriggerObject& TrigObj(TOC[hltf]);
	_trig_HLT_eta[hlt_counter]    = TrigObj.eta();
	_trig_HLT_phi[hlt_counter]    = TrigObj.phi();
	_trig_HLT_energy[hlt_counter] = TrigObj.energy();
	_trig_HLT_pt[hlt_counter]     = TrigObj.pt();
	_trig_HLT_name[hlt_counter]   = itrig;
	hlt_counter++;
      } // for loop on HLT objects
    } // if idfilter<trigevent size
  } // for loop on filters

  _trig_HLT_N = hlt_counter;
  if(hlt_counter>19) { _trig_HLT_N = 20; cout << "Number of HLT Objects>20, trig_HLT_N set to 20" << endl;}
	


  

  /////////////////////////////////////////////////////////
  ///                 Stage 2 Level 1                   ///
  /////////////////////////////////////////////////////////

  if(type_ != "MC") {

    /*Handle< BXVector<l1t::CaloTower> > towers;
    iEvent.getByToken(m_towerToken_,towers);
    
    int i_tow=0;
    
    for ( int ibx=towers->getFirstBX(); ibx<=towers->getLastBX(); ++ibx) {

      for ( auto itr = towers->begin(ibx); itr !=towers->end(ibx); ++itr ) {
	
        if (itr->hwPt()<=0) continue;
	
        cout << "Stage 2 Tower : " << " BX=" << 0 << " ipt=" << itr->hwPt() << " ieta=" << itr->hwEta() << " iphi=" << itr->hwPhi() << std::endl;
	
	_Stage2_tower_n++;
	_Stage2_tower_hwPt[i_tow] = itr->hwPt();
	_Stage2_tower_hwEta[i_tow] = itr->hwEta();
	_Stage2_tower_hwPhi[i_tow] = itr->hwPhi();
	i_tow++;

      }

      }



  Handle< BXVector<l1t::EGamma> > mpegs;
  iEvent.getByToken(m_mpEGToken_,mpegs);

  int i_mpeg=0;
  
  for ( auto itr = mpegs->begin(0); itr != mpegs->end(0); ++itr ) {

    cout << "Stage 2 MP EG : " << " BX=" << 0 << " ipt=" << itr->hwPt() << " ieta=" << itr->hwEta() << " iphi=" << itr->hwPhi() << std::endl;
    
    if(itr->hwPt() > 0){
      _Stage2_mpeg_hwPt[i_mpeg] = itr->hwPt();
      _Stage2_mpeg_hwEta[i_mpeg] = itr->hwEta();
      _Stage2_mpeg_hwPhi[i_mpeg] = itr->hwPhi();
      i_mpeg++;
      
    }	  
    
  }
   
  _Stage2_mpeg_n = i_mpeg;*/

   


  Handle< BXVector<l1t::EGamma> > egs;
  iEvent.getByToken(m_egToken_,egs);

  int i_eg=0;

  for ( auto itr = egs->begin(0); itr != egs->end(0); ++itr ) {

    cout << "Stage 2 EG : " << " BX=" << 0 << " ipt=" << itr->hwPt() << " ieta=" << itr->hwEta() << " iphi=" << itr->hwPhi() << std::endl;	  

    if(itr->hwPt() > 0){

      _Stage2_eg_hwPt[i_eg] = itr->hwPt();
      _Stage2_eg_hwEta[i_eg] = itr->hwEta();
      _Stage2_eg_hwPhi[i_eg] = itr->hwPhi();

      //Physical units
      _Stage2_eg_et[i_eg] = itr->pt();
      _Stage2_eg_eta[i_eg] = itr->eta();
      _Stage2_eg_phi[i_eg] = itr->phi();
      _Stage2_eg_isoflag[i_eg] = itr->hwIso() & (0x1);

      i_eg++;
    
    }
      
  }

  _Stage2_eg_n = i_eg;

  }
      

  if(trigger_version_emul_=="Stage2"){
  
    Handle< BXVector<l1t::CaloTower> > towers_emul;
    iEvent.getByToken(m_towerToken_emul_,towers_emul);
    vector<l1t::CaloTower> twrs;
    
    int i_tow_emul=0;

    for ( auto itr = towers_emul->begin(0); itr !=towers_emul->end(0); ++itr ) {
    
      twrs.push_back(*itr);
      
      if (itr->hwPt()<=0) continue;
      
      //cout << "Stage 2 Tower Emul: " << " BX=" << 0 << " ipt=" << itr->hwPt() << " ieta=" << itr->hwEta() << " iphi=" << itr->hwPhi() << std::endl;
      //cout << "E=" << tower.hwEtEm() << " H=" << itr->hwEtHad() << std::endl;
      
      _Stage2_tower_emul_n++;
      _Stage2_tower_emul_hwPt[i_tow_emul] = itr->hwPt();
      _Stage2_tower_emul_hwEta[i_tow_emul] = itr->hwEta();
      _Stage2_tower_emul_hwPhi[i_tow_emul] = itr->hwPhi();
      _Stage2_tower_emul_hwEtEm[i_tow_emul] = itr->hwEtEm();
      _Stage2_tower_emul_hwEtHad[i_tow_emul] = itr->hwEtHad();
      i_tow_emul++;
      
    }
    
    int nrTowers = l1t::CaloTools::calNrTowers(-4,4,1,72,twrs,1,999,l1t::CaloTools::CALO);
    
    Handle< BXVector<l1t::CaloCluster> > cls;
    iEvent.getByToken(m_clusterToken_emul_,cls);
    vector<l1t::CaloCluster> clusters;
    for ( auto itr = cls->begin(0); itr != cls->end(0); ++itr ) 
      {
	clusters.push_back(*itr);
	
      }
  
    

    Handle< BXVector<l1t::EGamma> > mpegs_emul;
    iEvent.getByToken(m_mpEGToken_emul_,mpegs_emul);
    
    int i_mpeg_emul=0;
    
    for ( auto itr = mpegs_emul->begin(0); itr != mpegs_emul->end(0); ++itr ) {
      
      cout << "Stage 2 MP EG Emul: " << " BX=" << 0 << " ipt=" << itr->hwPt() << " ieta=" << itr->hwEta() << " iphi=" << itr->hwPhi() << std::endl;
      
      if(itr->hwPt() > 0){
	_Stage2_mpeg_emul_hwPt[i_mpeg_emul] = itr->hwPt();
	_Stage2_mpeg_emul_hwEta[i_mpeg_emul] = itr->hwEta();
	_Stage2_mpeg_emul_hwPhi[i_mpeg_emul] = itr->hwPhi();

	_Stage2_mpeg_emul_et[i_mpeg_emul] = itr->pt();
	_Stage2_mpeg_emul_eta[i_mpeg_emul] = itr->eta();
	_Stage2_mpeg_emul_phi[i_mpeg_emul] = itr->phi();
	
	const l1t::CaloCluster& clus = l1t::CaloTools::getCluster(clusters, itr->hwEta(), itr->hwPhi());
	int shape = 0;
	if( (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_N)) ) shape |= (0x1);
	if( (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_S)) ) shape |= (0x1<<1);
	if( clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT)  && (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_E))  ) shape |= (0x1<<2);
	if( !clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT) && (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_W))  ) shape |= (0x1<<2);
	if( clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT)  && (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_NE)) ) shape |= (0x1<<3);
	if( !clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT) && (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_NW)) ) shape |= (0x1<<3);
	if( clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT)  && (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_SE)) ) shape |= (0x1<<4);
	if( !clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT) && (clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_SW)) ) shape |= (0x1<<4);
	if( clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_NN) ) shape |= (0x1<<5);
	if( clus.checkClusterFlag(l1t::CaloCluster::INCLUDE_SS) ) shape |= (0x1<<6);
	
	_Stage2_mpeg_emul_shape[i_mpeg_emul] = shape;
	
	int qual = itr->hwQual();
	_Stage2_mpeg_emul_shapeID[i_mpeg_emul] = qual>>2 & (0x1);
	
	
	int isoLeftExtension = 2;
	int isoRightExtension = 2;
	
	if(clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT))
	  isoRightExtension++;
	else
	  isoLeftExtension++;
	
	int E9x6 = l1t::CaloTools::calHwEtSum(itr->hwEta(), itr->hwPhi(), twrs,
					      -isoLeftExtension,isoRightExtension, //localEtaMin, localEtaMax
					      -4,4, //localPhiMin, localPhiMax
					      32); //iEtaAbsMax
	
	int etaSide = clus.checkClusterFlag(l1t::CaloCluster::TRIM_LEFT) ? 1 : -1;
	int phiSide = itr->hwEta()>0 ? 1 : -1;
	
	int ecalHwFootPrint = l1t::CaloTools::calHwEtSum(itr->hwEta(),itr->hwPhi(),twrs,
							 0,0,
							 -2,2,
							 32,l1t::CaloTools::ECAL) +
	  l1t::CaloTools::calHwEtSum(itr->hwEta(),itr->hwPhi(),twrs,
				     etaSide,etaSide,
				     -2,2,
				     32,l1t::CaloTools::ECAL);
	int hcalHwFootPrint = l1t::CaloTools::calHwEtSum(itr->hwEta(),itr->hwPhi(),twrs,
							 0,0,
							 0,0,
							 32,l1t::CaloTools::HCAL) +
	  l1t::CaloTools::calHwEtSum(itr->hwEta(),itr->hwPhi(),twrs,
				     0,0,
				     phiSide,phiSide,
				     32,l1t::CaloTools::HCAL);
	int FootPrint = ecalHwFootPrint+hcalHwFootPrint;
	
	_Stage2_mpeg_emul_hwIsoSum[i_mpeg_emul] = E9x6-FootPrint;
	_Stage2_mpeg_emul_nTT[i_mpeg_emul] = nrTowers;
	_Stage2_mpeg_emul_hOverERatio[i_mpeg_emul] = clus.hOverE();
	
	_Stage2_mpeg_emul_isoflag[i_mpeg_emul] = itr->hwIso() & (0x1);
	
	i_mpeg_emul++;
	
      }
      
    }
    
    _Stage2_mpeg_emul_n = i_mpeg_emul;
    
    

    Handle< BXVector<l1t::EGamma> > egs_emul;
    iEvent.getByToken(m_egToken_emul_,egs_emul);
    
    int i_eg_emul=0;
    
    for ( auto itr = egs_emul->begin(0); itr != egs_emul->end(0); ++itr ) {
      
      cout << "Stage 2 EG Emul: " << " BX=" << 0 << " ipt=" << itr->hwPt() << " ieta=" << itr->hwEta() << " iphi=" << itr->hwPhi() << std::endl;	  
      
      if(itr->hwPt() > 0){
      
	_Stage2_eg_emul_hwPt[i_eg_emul] = itr->hwPt();
	_Stage2_eg_emul_hwEta[i_eg_emul] = itr->hwEta();
	_Stage2_eg_emul_hwPhi[i_eg_emul] = itr->hwPhi();
	
	//Physical units
	_Stage2_eg_emul_et[i_eg_emul] = itr->pt();
	_Stage2_eg_emul_eta[i_eg_emul] = itr->eta();
	_Stage2_eg_emul_phi[i_eg_emul] = itr->phi();
	
	i_eg_emul++;
	
      }
      
    }
    
    _Stage2_eg_emul_n = i_eg_emul;
    
  }


  if(trigger_version_emul_=="Stage1"){
    
    Handle< BXVector<l1t::EGamma> > egs_Stage1_emul;
    iEvent.getByToken(m_egToken_Stage1_emul_,egs_Stage1_emul);
    
    int i_eg_Stage1_emul=0;
    
    for ( auto itr = egs_Stage1_emul->begin(0); itr != egs_Stage1_emul->end(0); ++itr ) {
      
      cout << "Stage 1 EG Emul: " << " BX=" << 0 << " ipt=" << itr->hwPt() << " ieta=" << itr->hwEta() << " iphi=" << itr->hwPhi() << std::endl;	  
      
      if(itr->hwPt() > 0){
	
	_Stage1_eg_emul_hwPt[i_eg_Stage1_emul] = itr->hwPt();
	_Stage1_eg_emul_hwEta[i_eg_Stage1_emul] = itr->hwEta();
	_Stage1_eg_emul_hwPhi[i_eg_Stage1_emul] = itr->hwPhi();
	
	//Physical units
	_Stage1_eg_emul_et[i_eg_Stage1_emul] = itr->pt();
	_Stage1_eg_emul_eta[i_eg_Stage1_emul] = itr->eta();
	_Stage1_eg_emul_phi[i_eg_Stage1_emul] = itr->phi();
	
	_Stage1_eg_emul_isoflag[i_eg_Stage1_emul] = itr->hwIso() & (0x1);
	
	i_eg_Stage1_emul++;
	
      }
      
    }
    
    _Stage1_eg_emul_n = i_eg_Stage1_emul;

  }



} // end of FillTrigger



// ====================================================================================
void SimpleNtpleCustom::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // caloMET object (negative vector sum of calorimeter towers)
  edm::Handle< edm::View<reco::CaloMET> > caloMEThandle;
  iEvent.getByLabel("met", caloMEThandle);
	
  // MET object that corrects the basic calorimeter MET for muons
  edm::Handle< edm::View<reco::CaloMET> > muCorrMEThandle;
  iEvent.getByLabel("corMetGlobalMuons", muCorrMEThandle);
	
  // MET object that corrects the basic calorimeter MET for muons and tracks
  edm::Handle< edm::View<reco::MET> > tcMEThandle;
  iEvent.getByLabel("tcMet", tcMEThandle);
	
  // MET object built as the (negative) vector sum of all particles (PFCandidates) reconstructed in the event
  edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
  iEvent.getByLabel("pfMet", pfMEThandle);
	
  // CALO MET
  _met_calo_et  = (caloMEThandle->front() ).et();
  _met_calo_px  = (caloMEThandle->front() ).px();
  _met_calo_py  = (caloMEThandle->front() ).py();
  _met_calo_phi = (caloMEThandle->front() ).phi();
  _met_calo_set = (caloMEThandle->front() ).sumEt();
  _met_calo_sig = (caloMEThandle->front() ).mEtSig();
	
  // CALOMU MET
  _met_calomu_et  = (muCorrMEThandle->front() ).et();
  _met_calomu_px  = (muCorrMEThandle->front() ).px();
  _met_calomu_py  = (muCorrMEThandle->front() ).py();
  _met_calomu_phi = (muCorrMEThandle->front() ).phi();
  _met_calomu_set = (muCorrMEThandle->front() ).sumEt();
  _met_calomu_sig = (muCorrMEThandle->front() ).mEtSig();
	
  // TC MET
  _met_tc_et  = (tcMEThandle->front() ).et();
  _met_tc_px  = (tcMEThandle->front() ).px();
  _met_tc_py  = (tcMEThandle->front() ).py();
  _met_tc_phi = (tcMEThandle->front() ).phi();
  _met_tc_set = (tcMEThandle->front() ).sumEt();
  _met_tc_sig = (tcMEThandle->front() ).mEtSig();
	
  // PFMET
  _met_pf_et  = (pfMEThandle->front() ).et();
  _met_pf_px  = (pfMEThandle->front() ).px();
  _met_pf_py  = (pfMEThandle->front() ).py();
  _met_pf_phi = (pfMEThandle->front() ).phi();
  _met_pf_set = (pfMEThandle->front() ).sumEt();
  _met_pf_sig = (pfMEThandle->front() ).mEtSig();
	
} // end of Fill MET


// ====================================================================================
void SimpleNtpleCustom::FillEle(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  edm::Handle<reco::GsfElectronCollection> EleHandle ;
  //iEvent.getByLabel (EleTag_.label(),EleHandle) ;
  iEvent.getByToken (EleToken_,EleHandle) ;
	
  edm::Handle<reco::PFCandidateCollection> PfEleHandle;
  //iEvent.getByLabel("particleFlow", PfEleHandle);
  iEvent.getByToken(PfEleToken_, PfEleHandle);
  
  if(PrintDebug_) std::cout << "Ele ID..." << std::endl;

  std::vector<edm::Handle<edm::ValueMap<bool> > > eleVIDDecisionHandles(4) ;
  /*iEvent.getByLabel  (eleVetoIdMapTag_ , eleVIDDecisionHandles[0]) ;
  iEvent.getByLabel  (eleLooseIdMapTag_ , eleVIDDecisionHandles[1]) ;
  iEvent.getByLabel  (eleMediumIdMapTag_ , eleVIDDecisionHandles[2]) ;
  iEvent.getByLabel  (eleTightIdMapTag_ , eleVIDDecisionHandles[3]) ;*/
  iEvent.getByToken  (eleVetoIdMapToken_ , eleVIDDecisionHandles[0]) ;
  iEvent.getByToken  (eleLooseIdMapToken_ , eleVIDDecisionHandles[1]) ;
  iEvent.getByToken  (eleMediumIdMapToken_ , eleVIDDecisionHandles[2]) ;
  iEvent.getByToken  (eleTightIdMapToken_ , eleVIDDecisionHandles[3]) ;

  if(PrintDebug_) std::cout << "Beam Spot" << std::endl;
  reco::BeamSpot bs;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  //iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);
  if (beamSpotHandle.isValid() ) bs = *beamSpotHandle;
  else cout << "No Beam spot ! " << endl;
  

  if(PrintDebug_) std::cout << "Calo Topo" << std::endl;
  //calo topology
  const CaloTopology * topology ;
  ///const EcalChannelStatus *chStatus ;
  edm::Handle< EcalRecHitCollection > reducedEBRecHits;
  edm::Handle< EcalRecHitCollection > reducedEERecHits;
	
  unsigned long long cacheIDTopo_=0;
  edm::ESHandle<CaloTopology> theCaloTopo;
  if (cacheIDTopo_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
    cacheIDTopo_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
    iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  }
  topology = theCaloTopo.product() ;
	
  /*edm::ESHandle<EcalChannelStatus> pChannelStatus;
  iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);

  //for42x	
  unsigned long long cacheSevLevel = 0;
  edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
  if(cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()){
    cacheSevLevel = iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
  }
  const EcalSeverityLevelAlgo* sl=sevLevel.product();*/
  

  // geometry (used for L1 trigger)                                                                                                                
  edm::ESHandle<CaloSubdetectorGeometry> theEndcapGeometry_handle, theBarrelGeometry_handle;
	
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap",theEndcapGeometry_handle);
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
	
  iSetup.get<IdealGeometryRecord>().get(eTTmap_);
  theEndcapGeometry_ = &(*theEndcapGeometry_handle);
  theBarrelGeometry_ = &(*theBarrelGeometry_handle);
	
  
  // reduced rechits
  /*if(!aod_){
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEE"), reducedEERecHits ) ;
    }
  else{
  iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEB"), reducedEBRecHits );
  iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEE"), reducedEERecHits ) ;
    }*/

  iEvent.getByToken( reducedEBRecHitsToken_, reducedEBRecHits );
  iEvent.getByToken( reducedEERecHitsToken_, reducedEERecHits ) ;

  edm::Handle<reco::TrackCollection> tracks_h;
  //iEvent.getByLabel("generalTracks", tracks_h);
  iEvent.getByToken(tracksToken_, tracks_h);
	
  edm::Handle<DcsStatusCollection> dcsHandle;
  //iEvent.getByLabel(dcsTag_, dcsHandle);
  iEvent.getByToken(dcsToken_, dcsHandle);
  double evt_bField;
  // need the magnetic field
  //
  // if isData then derive bfield using the
  // magnet current from DcsStatus
  // otherwise take it from the IdealMagneticFieldRecord
  if(PrintDebug_) std::cout << "B Field" << std::endl;
  if(type_ == "DATA" ) 
    {
      // scale factor = 3.801/18166.0 which are
      // average values taken over a stable two
      // week period
      if ((*dcsHandle).size() != 0 ) {	
	float currentToBFieldScaleFactor = 2.09237036221512717e-04;
	float current = (*dcsHandle)[0].magnetCurrent();
	evt_bField = current*currentToBFieldScaleFactor;
      }
      else {	
	edm::ESHandle<MagneticField> magneticField;
	iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
	GlobalPoint gPoint(0.,0.,0.);
	evt_bField = magneticField->inTesla(gPoint).z();
      }
    }
  else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    GlobalPoint gPoint(0.,0.,0.);
    evt_bField = magneticField->inTesla(gPoint).z();
  }
	
  

  // Beam Spot Information
  if(PrintDebug_) std::cout << "Beam Spot again" << std::endl;
  BS_x = bs.position().x();
  BS_y = bs.position().y();
  BS_z = bs.position().z();
	
  BS_dz = bs.sigmaZ();
  BS_dydz = bs.dydz();
  BS_dxdz = bs.dxdz();
	
  BS_bw_x = bs.BeamWidthX();
  BS_bw_y = bs.BeamWidthY();
	

  
	
  // -----------------------------
  // Masked Towers
  // -----------------------------
  /*if(PrintDebug_) std::cout << "Mask Tower" << std::endl;
  // !!! Idealy, should be in FillTrigger... but not possible for now !!!
  //Adding RCT mask  
  // list of RCT channels to mask                                    
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
	_trig_nMaskedRCT = n0MaskedRCT;*/
	
  //Adding TT mask                                                                                                                                                      
  // list of towers masked for trigger                                                                                                                                  
	
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
	
  
  /*	
  // ----------------------------------------------
  //  Get Seeds collection
  // ----------------------------------------------
  if(!aod_){
    edm::Handle<reco::ElectronSeedCollection> elSeeds;
    iEvent.getByLabel(SeedTag_,elSeeds);
	
    if(elSeeds.product()->size() < 100) ele_nSeed = elSeeds.product()->size();
    else ele_nSeed = 100;
	
    if(ele_nSeed > 0){
      reco::ElectronSeedCollection::const_iterator MyS_seed = (*elSeeds).begin();
      for(int counterSeed=0; counterSeed<ele_nSeed; ++counterSeed){
			
	if(MyS_seed->isEcalDriven()) ele_SeedIsEcalDriven[counterSeed] = 1;
	if(MyS_seed->isTrackerDriven()) ele_SeedIsTrackerDriven[counterSeed] = 1;
			
	ele_SeedSubdet2[counterSeed] = int(MyS_seed->subDet2());
	//to avoid some inf values//
	if(fabs(MyS_seed->dPhi2Pos()) < 100.) ele_SeedDphi2Pos[counterSeed] = double(MyS_seed->dPhi2Pos());
	if(fabs(MyS_seed->dRz2Pos()) < 100.)  ele_SeedDrz2Pos[counterSeed]  = double(MyS_seed->dRz2Pos());
	if(fabs(MyS_seed->dPhi2()) < 100.) ele_SeedDphi2Neg[counterSeed] = double(MyS_seed->dPhi2());
	if(fabs(MyS_seed->dRz2()) < 100.)  ele_SeedDrz2Neg[counterSeed]  = double(MyS_seed->dRz2());
			
	ele_SeedSubdet1[counterSeed] = int(MyS_seed->subDet1());
	//to avoid some inf values//
	if(fabs(MyS_seed->dPhi1Pos()) < 100.) ele_SeedDphi1Pos[counterSeed] = double(MyS_seed->dPhi1Pos());
	if(fabs(MyS_seed->dRz1Pos()) < 100.)  ele_SeedDrz1Pos[counterSeed]  = double(MyS_seed->dRz1Pos());
	if(fabs(MyS_seed->dPhi1()) < 100.) ele_SeedDphi1Neg[counterSeed] = double(MyS_seed->dPhi1());
	if(fabs(MyS_seed->dRz1()) < 100.)  ele_SeedDrz1Neg[counterSeed]  = double(MyS_seed->dRz1());
			
	++MyS_seed;
      } // loop on seed
    } // if nSeed>0
  }

 

  // ----------------------------------------------
  //  Get MC information
  // ----------------------------------------------
  edm::Handle<edm::HepMCProduct> HepMCEvt;
  if(type_ == "MC" && aod_ == false) iEvent.getByLabel(MCTag_, HepMCEvt);
  const HepMC::GenEvent* MCEvt = 0; 
  HepMC::GenParticle* genPc = 0;
  HepMC::FourVector pAssSim;
  if(type_ == "MC" && aod_ == false) { MCEvt = HepMCEvt->GetEvent();
    _MC_pthat =  HepMCEvt -> GetEvent() -> event_scale();}
	
  if(type_ == "MC" && aod_ == true) {
    edm::Handle< GenEventInfoProduct > HepMCEvt;
    iEvent.getByLabel(MCTag_, HepMCEvt);
    if(HepMCEvt->hasBinningValues()) _MC_pthat = (HepMCEvt->binningValues())[0];
    else  _MC_pthat = 0.0;
  }

  edm::Handle<reco::GenParticleCollection> genParticlesColl;
  if(aod_ == true && type_ == "MC") iEvent.getByLabel("genParticles", genParticlesColl);

 

	  
  // for H/E

  towersH_ = new edm::Handle<CaloTowerCollection>() ;
  if (!iEvent.getByLabel(hcalTowers_,*towersH_))
    { edm::LogError("ElectronHcalHelper::readEvent")<<"failed to get the hcal towers of label "<<hcalTowers_ ; }
  // H/E, with Rcone = 0.05 & different ET threshold
  EgammaTowerIsolation * towerIso1_00615_0  = new EgammaTowerIsolation(0.0615, 0., 0.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_00615_0  = new EgammaTowerIsolation(0.0615, 0., 0.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_005_0  = new EgammaTowerIsolation(0.05, 0., 0.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_005_0  = new EgammaTowerIsolation(0.05, 0., 0.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_005_1  = new EgammaTowerIsolation(0.05, 0., 1.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_005_1  = new EgammaTowerIsolation(0.05, 0., 1.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_005_15 = new EgammaTowerIsolation(0.05, 0., 1.5, 1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_005_15 = new EgammaTowerIsolation(0.05, 0., 1.5, 2, towersH_->product()) ;
  // H/E, with Rcone = 0.1 & different ET threshold
  EgammaTowerIsolation * towerIso1_01_0   = new EgammaTowerIsolation(0.1, 0., 0.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_01_0   = new EgammaTowerIsolation(0.1, 0., 0.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_01_1   = new EgammaTowerIsolation(0.1, 0., 1.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_01_1   = new EgammaTowerIsolation(0.1, 0., 1.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_01_15  = new EgammaTowerIsolation(0.1, 0., 1.5, 1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_01_15  = new EgammaTowerIsolation(0.1, 0., 1.5, 2, towersH_->product()) ;
  // H/E, with Rcone = 0.15 & different ET threshold
  EgammaTowerIsolation * towerIso1_015_1  = new EgammaTowerIsolation(0.15, 0., 1.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_015_1  = new EgammaTowerIsolation(0.15, 0., 1.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_015_15 = new EgammaTowerIsolation(0.15, 0., 1.5, 1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_015_15 = new EgammaTowerIsolation(0.15, 0., 1.5, 2, towersH_->product()) ;
	  
  //for NEW HCAL Iso 
  EgammaTowerIsolation * towerIso1_00615Ring03_0  = new EgammaTowerIsolation(0.3, 0.0615, 0.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_00615Ring03_0  = new EgammaTowerIsolation(0.3, 0.0615, 0.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_005Ring03_0  = new EgammaTowerIsolation(0.3, 0.05, 0.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_005Ring03_0  = new EgammaTowerIsolation(0.3, 0.05, 0.,  2, towersH_->product()) ;
  EgammaTowerIsolation * towerIso1_0Ring03_0  = new EgammaTowerIsolation(0.3, 0., 0.,  1, towersH_->product()) ;
  EgammaTowerIsolation * towerIso2_0Ring03_0  = new EgammaTowerIsolation(0.3, 0., 0.,  2, towersH_->product()) ;
	  
  

  float deta = -20.;
  float dphi = -20.;*/
	
  if(EleHandle->size() < 10 ){ ele_N = EleHandle->size(); }
  else {ele_N = 10;}
  TClonesArray &electrons = *m_electrons;

  int counter = 0;
	


  // ----------------------------------------------
  //  Loop on Electrons
  // ----------------------------------------------
  int nTow=0;
  int nReg=0;
	
  for(int i=0; i< ele_N; i++){
		
    reco::SuperClusterRef sclRef = (*EleHandle)[i].superCluster();
    math::XYZPoint sclPos = (*EleHandle)[i].superClusterPosition();
    if (!(*EleHandle)[i].ecalDrivenSeed() && (*EleHandle)[i].trackerDrivenSeed()) 
    {
 
      if (!((*EleHandle)[i].parentSuperCluster().isNull()))
      {
   
        sclRef = (*EleHandle)[i].parentSuperCluster();
 
      }
      else
      {
 
             continue;
      }
    }



	//FANBO_setMonmentum	
    edm::Ref<reco::GsfElectronCollection> electronEdmRef(EleHandle,i);
    setMomentum (myvector, (*EleHandle)[i].p4());
    new (electrons[counter]) TLorentzVector (myvector);

    ele_VetoIdDecisions[counter] = (*(eleVIDDecisionHandles[0]))[electronEdmRef]; 
    ele_LooseIdDecisions[counter] = (*(eleVIDDecisionHandles[1]))[electronEdmRef]; 
    ele_MediumIdDecisions[counter] = (*(eleVIDDecisionHandles[2]))[electronEdmRef]; 
    ele_TightIdDecisions[counter] = (*(eleVIDDecisionHandles[3]))[electronEdmRef]; 

    ele_pT[counter] = myvector.Pt();
		
    ele_echarge[counter] = (*EleHandle)[i].charge(); 
    ele_he[counter]      = (*EleHandle)[i].hadronicOverEm() ;
		
    ele_eseedpout[counter] = (*EleHandle)[i].eSeedClusterOverPout();
    ele_ep[counter]        = (*EleHandle)[i].eSuperClusterOverP() ;        
    ele_eseedp[counter]    = (*EleHandle)[i].eSeedClusterOverP() ;         
    ele_eelepout[counter]  = (*EleHandle)[i].eEleClusterOverPout() ;       
		
    ele_pin_mode[counter]    = (*EleHandle)[i].trackMomentumAtVtx().R() ; 
    ele_pout_mode[counter]   = (*EleHandle)[i].trackMomentumOut().R() ; 

    ele_calo_energy[counter] = (*EleHandle)[i].caloEnergy() ;
    ele_pTin_mode[counter]   = (*EleHandle)[i].trackMomentumAtVtx().Rho() ; 
    ele_pTout_mode[counter]  = (*EleHandle)[i].trackMomentumOut().Rho() ; 

    /*if(!aod_){
      ele_pin_mean[counter]    = (*EleHandle)[i].gsfTrack()->innerMomentum().R() ; 
      ele_pout_mean[counter]   = (*EleHandle)[i].gsfTrack()->outerMomentum().R(); 

      ele_pTin_mean[counter]   = (*EleHandle)[i].gsfTrack()->innerMomentum().Rho() ; 
      ele_pTout_mean[counter]  = (*EleHandle)[i].gsfTrack()->outerMomentum().Rho() ;
      }*/

 	
    double R=TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
    double Rt=TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
    ele_sclRawE[counter]   = sclRef->rawEnergy() ;
    ele_sclEpresh[counter] = 0. ;
    if ((*EleHandle)[i].isEE()) ele_sclEpresh[counter]   = sclRef->preshowerEnergy() ;
    ele_sclE[counter]   = sclRef->energy() ;
    ele_sclEt[counter]  = sclRef->energy()*(Rt/R) ;
    ele_sclEta[counter] = sclRef->eta() ;
    ele_sclPhi[counter] = sclRef->phi() ;
    ele_sclX[counter]  = sclPos.X();
    ele_sclY[counter] =  sclPos.Y();
    ele_sclZ[counter] =  sclPos.Z();
    
	

    // NEW H/E
    
    /*reco::SuperCluster EmSCCand = *sclRef;
    double HoE_00615_0  = towerIso1_00615_0->getTowerESum(&EmSCCand)  + towerIso2_00615_0->getTowerESum(&EmSCCand) ;
    double HoE_005_0  = towerIso1_005_0->getTowerESum(&EmSCCand)  + towerIso2_005_0->getTowerESum(&EmSCCand) ;
    double HoE_005_1  = towerIso1_005_1->getTowerESum(&EmSCCand)  + towerIso2_005_1->getTowerESum(&EmSCCand) ;
    double HoE_005_15 = towerIso1_005_15->getTowerESum(&EmSCCand) + towerIso2_005_15->getTowerESum(&EmSCCand) ;

    double HoE_01_0  = towerIso1_01_0->getTowerESum(&EmSCCand)  + towerIso2_01_0->getTowerESum(&EmSCCand) ;
    double HoE_01_1  = towerIso1_01_1->getTowerESum(&EmSCCand)  + towerIso2_01_1->getTowerESum(&EmSCCand) ;
    double HoE_01_15 = towerIso1_01_15->getTowerESum(&EmSCCand) + towerIso2_01_15->getTowerESum(&EmSCCand) ;
    
    double HoE_015_1  = towerIso1_015_1->getTowerESum(&EmSCCand)  + towerIso2_015_1->getTowerESum(&EmSCCand) ;
    double HoE_015_15 = towerIso1_015_15->getTowerESum(&EmSCCand) + towerIso2_015_15->getTowerESum(&EmSCCand) ;
		

    HoE_00615_0  /=         sclRef->energy() ;
    HoE_005_0  /= 	sclRef->energy() ;     
    HoE_005_1  /= 	sclRef->energy() ;     
    HoE_005_15 /= 	sclRef->energy() ;     
    HoE_01_0   /= 	sclRef->energy() ;     
    HoE_01_1   /= 	sclRef->energy() ;     
    HoE_01_15  /= 	sclRef->energy() ;              
    HoE_015_1  /= 	sclRef->energy() ;     
    HoE_015_15 /= 	sclRef->energy() ; 

    _ele_he_00615_0[counter]  = HoE_00615_0 ;
    _ele_he_005_0[counter]  = HoE_005_0 ;
    _ele_he_005_1[counter]  = HoE_005_1 ;	
    _ele_he_005_15[counter] = HoE_005_15 ;
    _ele_he_01_0[counter]   = HoE_01_0 ;
    _ele_he_01_1[counter]   = HoE_01_1 ;	
    _ele_he_01_15[counter]  = HoE_01_15 ;
    _ele_he_015_1[counter]  = HoE_015_1 ;	
    _ele_he_015_15[counter] = HoE_015_15 ;
		
	
    
    // 		total effective uncertainty on energy in GeV
    if (funcbase_) {
      ele_sclErr[counter] = funcbase_->getValue(*sclRef, 0);
      // positive uncertainty
      ele_sclErr_pos[counter] = funcbase_->getValue(*sclRef, 1);      
      // negative uncertainty
      ele_sclErr_neg[counter] = funcbase_->getValue(*sclRef, -1);      
      ele_trErr[counter]=(*EleHandle)[i].trackMomentumError();
      
      //for42X                 
      ele_momErr[counter]=(*EleHandle)[i].p4Error(GsfElectron::P4_COMBINATION);
      

      if (!(*EleHandle)[i].ecalDrivenSeed()) { /// no change if not ecaldriven
	ele_newmom[counter] = (*EleHandle)[i].p4().P();
	ele_newmomErr[counter] = ele_momErr[counter];
      }

      else { //if ecal driven special care for large errors
	if (ele_trErr[counter]/ele_pin_mode[counter] > 0.5 && ele_sclErr[counter]/ele_sclE[counter] <= 0.5) { //take E if sigmaE/E <=0.5 and sigmaP/P >0.5
	  ele_newmom[counter] = ele_sclE[counter];    ele_newmomErr[counter] = ele_sclErr[counter];	 
	}
	else if (ele_trErr[counter]/ele_pin_mode[counter] <= 0.5 && ele_sclErr[counter]/ele_sclE[counter] > 0.5){//take P if sigmaE/E > 0.5 and sigmaP/P <=0.5
	  ele_newmom[counter] = ele_pin_mode[counter];  ele_newmomErr[counter] = ele_trErr[counter];	  
	}
	else if (ele_trErr[counter]/ele_pin_mode[counter] > 0.5 && ele_sclErr[counter]/ele_sclE[counter] > 0.5){//take the lowest sigma/value if sigmaE/E >0.5 and sigmaP/P >0.5
	  if (ele_trErr[counter]/ele_pin_mode[counter] < ele_sclErr[counter]/ele_sclE[counter]) {
	    ele_newmom[counter] = ele_pin_mode[counter]; ele_newmomErr[counter] = ele_trErr[counter];
	  }
	  else{
	    ele_newmom[counter] = ele_sclE[counter]; ele_newmomErr[counter] = ele_sclErr[counter];
	  }
	}
	else { // if sigmaE/E <= 0.5 and sigmaP/P <=0.5 no change
	  ele_newmom[counter] = (*EleHandle)[i].p4().P();
	  ele_newmomErr[counter] = ele_momErr[counter];
	}
      }
    }
    
		
		
    ele_tr_atcaloX[counter] = (*EleHandle)[i].trackPositionAtCalo().x();
    ele_tr_atcaloY[counter] = (*EleHandle)[i].trackPositionAtCalo().y();
    ele_tr_atcaloZ[counter] = (*EleHandle)[i].trackPositionAtCalo().z();    

    if(!aod_){
      ele_firsthit_X[counter] = (*EleHandle)[i].gsfTrack()->innerPosition().x();
      ele_firsthit_Y[counter] = (*EleHandle)[i].gsfTrack()->innerPosition().y();
      ele_firsthit_Z[counter] = (*EleHandle)[i].gsfTrack()->innerPosition().z();
      }  */ 
		

    ele_deltaetaseed[counter] = (*EleHandle)[i].deltaEtaSeedClusterTrackAtCalo() ; 
    ele_deltaphiseed[counter] = (*EleHandle)[i].deltaPhiSeedClusterTrackAtCalo() ;  
    ele_deltaetaele[counter]  = (*EleHandle)[i].deltaEtaEleClusterTrackAtCalo() ;  
    ele_deltaphiele[counter]  = (*EleHandle)[i].deltaPhiEleClusterTrackAtCalo() ; 
    ele_deltaetain[counter]   = (*EleHandle)[i].deltaEtaSuperClusterTrackAtVtx();
    ele_deltaphiin[counter]   = (*EleHandle)[i].deltaPhiSuperClusterTrackAtVtx();   
		
    ele_sigmaietaieta[counter] = (*EleHandle)[i].sigmaIetaIeta() ; 
    ele_sigmaetaeta[counter]   = (*EleHandle)[i].sigmaEtaEta() ;
    ele_e15[counter]           = (*EleHandle)[i].e1x5() ;
    ele_e25max[counter]        = (*EleHandle)[i].e2x5Max() ;
    ele_e55[counter]           = (*EleHandle)[i].e5x5() ;
    const EcalRecHitCollection * reducedRecHits = 0 ;
    if ((*EleHandle)[i].isEB())  
      reducedRecHits = reducedEBRecHits.product() ; 
    else 
      reducedRecHits = reducedEERecHits.product() ;
    const reco::CaloCluster & seedCluster = *(*EleHandle)[i].superCluster()->seed() ;
    ele_e1[counter]            = EcalClusterTools::eMax(seedCluster,reducedRecHits)  ;
    ele_e33[counter]           = EcalClusterTools::e3x3(seedCluster,reducedRecHits,topology)  ;

	
		
    ele_fbrem[counter] = (*EleHandle)[i].fbrem() ;
    ele_mva[counter]   = (*EleHandle)[i].mva_e_pi() ;
		
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
    ele_eClass[counter]   = (*EleHandle)[i].classification() ;
    ele_vertex_x[counter] = (*EleHandle)[i].vertex().x();
    ele_vertex_y[counter] = (*EleHandle)[i].vertex().y();
    ele_vertex_z[counter] = (*EleHandle)[i].vertex().z();
    //if(!aod_){
    ele_missing_hits[counter] = (*EleHandle)[i].gsfTrack()->numberOfLostHits();
    ele_lost_hits[counter]    = (*EleHandle)[i].gsfTrack()->numberOfValidHits() ;
    ele_chi2_hits[counter]    = (*EleHandle)[i].gsfTrack()->normalizedChi2() ;
		
    ele_dxyB[counter] = (*EleHandle)[i].gsfTrack()->dxy(bs.position()) ;
    ele_dxy[counter]  = (*EleHandle)[i].gsfTrack()->dxy() ;
    ele_dzB[counter]  = (*EleHandle)[i].gsfTrack()->dz(bs.position()) ;
    ele_dz[counter]   = (*EleHandle)[i].gsfTrack()->dz() ;
    ele_dszB[counter] = (*EleHandle)[i].gsfTrack()->dsz(bs.position()) ;
    ele_dsz[counter]  = (*EleHandle)[i].gsfTrack()->dsz() ;
		
    ele_dzPV[counter] = (*EleHandle)[i].gsfTrack()->dz(math::XYZPoint(vertexPosition));
    ele_dzPV_error[counter] = (*EleHandle)[i].gsfTrack()->dzError();
    ele_dxyPV[counter] = (*EleHandle)[i].gsfTrack()->dxy(math::XYZPoint(vertexPosition));
    ele_dxyPV_error[counter] = (*EleHandle)[i].gsfTrack()->dxyError();
    ele_dszPV[counter] = (*EleHandle)[i].gsfTrack()->dsz(math::XYZPoint(vertexPosition));
    ele_dszPV_error[counter] = (*EleHandle)[i].gsfTrack()->dszError();
		
    ele_track_x[counter] = (*EleHandle)[i].gsfTrack()->vx();
    ele_track_y[counter] = (*EleHandle)[i].gsfTrack()->vy();
    ele_track_z[counter] = (*EleHandle)[i].gsfTrack()->vz();
   

    // Isolation variables
    ele_tkSumPt_dr03[counter]              = (*EleHandle)[i].dr03TkSumPt() ;
    ele_ecalRecHitSumEt_dr03[counter]      = (*EleHandle)[i].dr03EcalRecHitSumEt() ;
    ele_hcalDepth1TowerSumEt_dr03[counter] = (*EleHandle)[i].dr03HcalDepth1TowerSumEt() ;
    ele_hcalDepth2TowerSumEt_dr03[counter] = (*EleHandle)[i].dr03HcalDepth2TowerSumEt() ;
    ele_tkSumPt_dr04[counter]              = (*EleHandle)[i].dr04TkSumPt() ;
    ele_ecalRecHitSumEt_dr04[counter]      = (*EleHandle)[i].dr04EcalRecHitSumEt() ;
    ele_hcalDepth1TowerSumEt_dr04[counter] = (*EleHandle)[i].dr04HcalDepth1TowerSumEt() ;
    ele_hcalDepth2TowerSumEt_dr04[counter] = (*EleHandle)[i].dr04HcalDepth2TowerSumEt() ;
		
    //NEW HCAL Isolation                                                       
    /*double HcalIso_00615Ring03_0  = (towerIso1_00615Ring03_0->getTowerEtSum(&((*EleHandle)[i]))  +
				     towerIso2_00615Ring03_0->getTowerEtSum(&((*EleHandle)[i])) );
    double HcalIso_005Ring03_0  = (towerIso1_005Ring03_0->getTowerEtSum(&((*EleHandle)[i]))  +
				   towerIso2_005Ring03_0->getTowerEtSum(&((*EleHandle)[i])) );
    double HcalIso_0Ring03_0  = (towerIso1_0Ring03_0->getTowerEtSum(&((*EleHandle)[i]))  +
				 towerIso2_0Ring03_0->getTowerEtSum(&((*EleHandle)[i])) );

    ele_hcalDepth1plus2TowerSumEt_00615dr03[counter] = HcalIso_00615Ring03_0;
    ele_hcalDepth1plus2TowerSumEt_005dr03[counter] = HcalIso_005Ring03_0;
    ele_hcalDepth1plus2TowerSumEt_0dr03[counter] = HcalIso_0Ring03_0;*/


  


    ele_ECAL_fbrem[counter]   = sclRef->phiWidth()/sclRef->etaWidth();
 
		
    // pflow combinaison
    ele_PFcomb[counter] = 0.;
    ele_PFcomb_Err[counter] = 0.;
    std::vector<reco::PFCandidate> candidates = (*PfEleHandle.product());
    for (std::vector<reco::PFCandidate>::iterator it = candidates.begin(); it != candidates.end(); ++it)   {
      reco::PFCandidate::ParticleType type = (*it).particleId();
      // here you can ask for particle type, mu,e,gamma
      if ( type == reco::PFCandidate::e) {

	if (!(*it).gsfTrackRef().isNull() && ((*it).gsfTrackRef()->p() == (*EleHandle)[i].gsfTrack()->p())){	 
	  ele_PFcomb[counter] = (*it).energy();
	}
      }
    }   
 

    ele_PFcomb_Err[counter]   =(*EleHandle)[i].p4Error(GsfElectron::P4_PFLOW_COMBINATION);
    if (!(*EleHandle)[i].parentSuperCluster().isNull()) ele_PF_SCenergy[counter]   = (*EleHandle)[i].parentSuperCluster()->energy();
                
    ele_PF_SCenergy_Err[counter]   = 0 ; //not implemented for the moment


    // Conversion Removal
    ele_isConversion[counter] = IsConv ((*EleHandle)[i]);
    ConversionFinder convFinder;

    ConversionInfo convInfo = convFinder.getConversionInfo((*EleHandle)[i], tracks_h, evt_bField);
    ele_conv_dist[counter] = convInfo.dist();
    ele_conv_dcot[counter] = convInfo.dcot();

    ele_convFound[counter] = 0;
    // Conversion 36X
    if ( convFinder.isFromConversion(convInfo, 0.02, 0.02) ) ele_convFound[counter] = 1;
		


    // ------------------------------
    // For L1 Trigger, Clemy's stuff
    // ------------------------------
    //LOOP MATCHING ON L1 trigger 
		
    //modif-alex l1 matching
    //LOOP MATCHING ON L1 trigger 
		
    nTow=0;
    nReg=0;
    for(int icc = 0; icc < 50; ++icc) {
      _ele_TTetaVect[counter][icc] = -999;
      _ele_TTphiVect[counter][icc] = -999;
      _ele_TTetVect[counter][icc] = 0.;
    }
    _ele_TTetaSeed[counter] = -999;
    _ele_TTphiSeed[counter] = -999;
    _ele_TTetSeed[counter] = 0.;    

    for(int icc = 0; icc < 10; ++icc) {
      _ele_RCTetaVect[counter][icc] = -999;
      _ele_RCTphiVect[counter][icc] = -999;
      _ele_RCTetVect[counter][icc] = 0.;
      _ele_RCTL1isoVect[counter][icc] = -999;
      _ele_RCTL1nonisoVect[counter][icc] = -999;
      _ele_RCTL1isoVect_M[counter][icc] = -999;
      _ele_RCTL1nonisoVect_M[counter][icc] = -999;
      _ele_L1Stage1_emul_isoVect[counter][icc] = -999;
      _ele_L1Stage1_emul_nonisoVect[counter][icc] = -999;
     }
		
    for (reco::CaloCluster_iterator clus = sclRef->clustersBegin () ;
	 clus != sclRef->clustersEnd () ;
	 ++clus){
      std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions() ; //get these from the cluster                                            
      //loop on xtals in cluster                                                                                                                                  
      for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ;
	   detitr != clusterDetIds.end () ;
	   ++detitr)
	{
	  //Here I use the "find" on a digi collection... I have been warned...                                                                                   
	  if ( (detitr -> first).det () != DetId::Ecal)
	    {
	      std::cout << " det is " << (detitr -> first).det () << std::endl ;
	      continue ;
	    }
	  EcalRecHitCollection::const_iterator thishit;
	  EcalRecHit myhit;
	  EcalTrigTowerDetId towid;
	  float thetahit;
	  if ( (detitr -> first).subdetId () == EcalBarrel)
	    {
	      thishit = reducedRecHits->find ( (detitr -> first) ) ;
	      if (thishit == reducedRecHits->end ()) continue;
	      myhit = (*thishit) ;
	      EBDetId detid(thishit->id());
	      towid= detid.tower();
	      thetahit =  theBarrelGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	    }//barrel rechit
	  else {
	    if ( (detitr -> first).subdetId () == EcalEndcap)
	      {
		thishit = reducedRecHits->find ( (detitr -> first) ) ;
		if (thishit == reducedRecHits->end ()) continue;
		myhit = (*thishit) ;
		EEDetId detid(thishit->id());
		towid= (*eTTmap_).towerOf(detid);
		thetahit =  theEndcapGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	      }
	    else continue;
	    }//endcap rechit
				
				
	  int iETA   = towid.ieta();
	  int iPHI   = towid.iphi();
	  int iReta  = getGCTRegionEta(iETA);
	  int iRphi  = getGCTRegionPhi(iPHI);
	  double iET = myhit.energy()*sin(thetahit);
				
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
				
	  bool newReg = true;
	  if(nReg>0) {
	    for (int iReg=0; iReg<nReg; ++iReg) {
	      if(_ele_RCTetaVect[counter][iReg] == iReta && _ele_RCTphiVect[counter][iReg] == iRphi) {
		newReg = false;
		_ele_RCTetVect[counter][iReg] +=  iET;
	      }
	    }
	  }
	  if(newReg) {
	    _ele_RCTetaVect[counter][nReg] = iReta;
	    _ele_RCTphiVect[counter][nReg] = iRphi;
	    _ele_RCTetVect[counter][nReg]  =  iET;
					
	    // standard collection
	    for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
	      if(_trig_L1emIso_iphi[il1] == iRphi && _trig_L1emIso_ieta[il1] == iReta) 
                {
                   _ele_RCTL1isoVect[counter][nReg] = _trig_L1emIso_rank[il1];
                   std::cout<<"matched L1 em Iso energy="<<_trig_L1emIso_rank[il1]<<std::endl;
                   std::cout<<"matched L1 em Iso iphi="<<iRphi<<std::endl;
                   std::cout<<"matched L1 em Iso ieta="<<iReta<<std::endl;
                   std::cout<<"--------******-----end for one----"<<std::endl;



                }
	    }
	    for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	      if(_trig_L1emNonIso_iphi[il1] == iRphi && _trig_L1emNonIso_ieta[il1] == iReta) _ele_RCTL1nonisoVect[counter][nReg] = _trig_L1emNonIso_rank[il1];
	    }
					
	    // modified collection
	    for(int il1=0; il1<_trig_L1emIso_N_M; ++il1) {
	      if(_trig_L1emIso_iphi_M[il1] == iRphi && _trig_L1emIso_ieta_M[il1] == iReta) _ele_RCTL1isoVect_M[counter][nReg] = _trig_L1emIso_rank_M[il1];
	    }
	    for(int il1=0; il1<_trig_L1emNonIso_N_M; ++il1) {
	      if(_trig_L1emNonIso_iphi_M[il1] == iRphi && _trig_L1emNonIso_ieta_M[il1] == iReta) _ele_RCTL1nonisoVect_M[counter][nReg] = _trig_L1emNonIso_rank_M[il1];
	    }


	    // Stage 1 emul
	    for(int il1=0; il1<_Stage1_eg_emul_n; ++il1) {
	      int iRphi_Stage1 = _Stage1_eg_emul_hwPhi[il1];
	      int iReta_Stage1 = (_Stage1_eg_emul_hwEta[il1]<7) ? 11 + _Stage1_eg_emul_hwEta[il1] : 18-_Stage1_eg_emul_hwEta[il1];

	      if(iRphi_Stage1 == iRphi && iReta_Stage1 == iReta){
		if(_Stage1_eg_emul_isoflag[il1]==1) _ele_L1Stage1_emul_isoVect[counter][nReg] = _Stage1_eg_emul_et[il1];
		else _ele_L1Stage1_emul_nonisoVect[counter][nReg] = _Stage1_eg_emul_et[il1];
	      }

	    }

					
	    nReg++;
	  } // if newReg
				
				
	}//loop crystal
    }//loop cluster
		

    double TTetmax2 = 0.;
    int iTTmax2     = -1.;

    for (int iTow=0; iTow<nTow; ++iTow) {
      bool nomaskTT = true;

      for (ittpg=towerMap.begin();ittpg!=towerMap.end();++ittpg) {
	if ((*ittpg).second > 0)
	  {
	    EcalTrigTowerDetId  ttId((*ittpg).first);
	    if(ttId.ieta() == _ele_TTetaVect[counter][iTow] && ttId.iphi() == _ele_TTphiVect[counter][iTow]) {
	      
	      nomaskTT=false;
	    }
	  }
      }//loop trigger towers

      if(nomaskTT && _ele_TTetVect[counter][iTow] > TTetmax2) {
	iTTmax2 = iTow;
	TTetmax2 = _ele_TTetVect[counter][iTow];

	_ele_TTetaSeed[counter] = _ele_TTetaVect[counter][iTow];
	_ele_TTphiSeed[counter] = _ele_TTphiVect[counter][iTow];
	_ele_TTetSeed[counter] = _ele_TTetVect[counter][iTow];

      } // if nomask
    } // for loop on Towers
		

		
		
    if(iTTmax2>=0) {
      int TTetamax2 = getGCTRegionEta(_ele_TTetaVect[counter][iTTmax2]);
      int TTphimax2 = getGCTRegionPhi(_ele_TTphiVect[counter][iTTmax2]);
      _ele_RCTeta[counter] = TTetamax2;
      _ele_RCTphi[counter] = TTphimax2;
			
      // standard collection
      for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
	if(_trig_L1emIso_iphi[il1] == TTphimax2 && _trig_L1emIso_ieta[il1] == TTetamax2)       _ele_RCTL1iso[counter]    = _trig_L1emIso_rank[il1];
      }
      for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	if(_trig_L1emNonIso_iphi[il1] == TTphimax2 && _trig_L1emNonIso_ieta[il1] == TTetamax2) _ele_RCTL1noniso[counter] = _trig_L1emNonIso_rank[il1];
      }

      // modified collection
      for(int il1=0; il1<_trig_L1emIso_N_M; ++il1) {
	if(_trig_L1emIso_iphi_M[il1] == TTphimax2 && _trig_L1emIso_ieta_M[il1] == TTetamax2)       _ele_RCTL1iso_M[counter]    = _trig_L1emIso_rank_M[il1];
      }
      for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	if(_trig_L1emNonIso_iphi_M[il1] == TTphimax2 && _trig_L1emNonIso_ieta_M[il1] == TTetamax2) _ele_RCTL1noniso_M[counter] = _trig_L1emNonIso_rank_M[il1];
      }


      //Stage 1 Level 1

      for(int il1=0; il1<_Stage1_eg_emul_n; ++il1) {
	int iRphi_Stage1 = _Stage1_eg_emul_hwPhi[il1];
	int iReta_Stage1 = (_Stage1_eg_emul_hwEta[il1]<7) ? 11 + _Stage1_eg_emul_hwEta[il1] : 18-_Stage1_eg_emul_hwEta[il1];

	if(iRphi_Stage1 == TTphimax2 && iReta_Stage1 == TTetamax2){
	  _ele_L1Stage1_emul[counter] = _Stage1_eg_emul_et[il1];
	  _ele_L1Stage1_emul_isoflag[counter] = _Stage1_eg_emul_isoflag[il1];
	  _ele_L1Stage1_emul_eta[counter] = _Stage1_eg_emul_eta[il1];
	  _ele_L1Stage1_emul_phi[counter] = _Stage1_eg_emul_phi[il1];
	}
      }

      
      //Stage 2 Level 1

      /*for(int il1=0; il1<_Stage2_mpeg_n; ++il1) {
	
	int delta_phi = abs(_Stage2_mpeg_hwPhi[il1] - _ele_TTphiSeed[counter]);
	if(delta_phi>36)
	   delta_phi = 72 - delta_phi;

	int delta_eta = abs(_Stage2_mpeg_hwEta[il1] - _ele_TTetaSeed[counter]);
	if(_Stage2_mpeg_hwEta[il1]*_ele_TTetaSeed[counter]<0)
	   delta_eta--;

	if(delta_phi<=1 && delta_eta<=1)
	  _ele_L1Stage2[counter] = _Stage2_mpeg_hwPt[il1]/2.;
	  }*/

      for(int il1=0; il1<_Stage2_eg_n; ++il1) {
	
	int CaloPhi = getCaloPhi_fromGTPhi(_Stage2_eg_hwPhi[il1]);
	int delta_phi = abs(CaloPhi - _ele_TTphiSeed[counter]);
	if(delta_phi>36)
	   delta_phi = 72 - delta_phi;

	int CaloEta = getCaloEta_fromGTEta(_Stage2_eg_hwEta[il1]);
	int delta_eta = abs(CaloEta - _ele_TTetaSeed[counter]);
	if(CaloEta*_ele_TTetaSeed[counter]<0)
	   delta_eta--;

	if(delta_phi<=1 && delta_eta<=1){
	  _ele_L1Stage2[counter] = _Stage2_eg_hwPt[il1]/2.;
	  _ele_L1Stage2_isoflag[counter] = _Stage2_eg_isoflag[il1];
	}
      }


      for(int il1=0; il1<_Stage2_mpeg_emul_n; ++il1) {
	
	int delta_phi = abs(_Stage2_mpeg_emul_hwPhi[il1] - _ele_TTphiSeed[counter]);
	if(delta_phi>36)
	   delta_phi = 72 - delta_phi;

	int delta_eta = abs(_Stage2_mpeg_emul_hwEta[il1] - _ele_TTetaSeed[counter]);
	if(_Stage2_mpeg_emul_hwEta[il1]*_ele_TTetaSeed[counter]<0)
	   delta_eta--;

	if(delta_phi<=1 && delta_eta<=1){
	  _ele_L1Stage2_emul[counter] = _Stage2_mpeg_emul_hwPt[il1]/2;

	  _ele_L1Stage2_emul_ieta[counter] = _Stage2_mpeg_emul_hwEta[il1];
	  _ele_L1Stage2_emul_hwPt[counter] = _Stage2_mpeg_emul_hwPt[il1];
	  _ele_L1Stage2_emul_shape[counter] = _Stage2_mpeg_emul_shape[il1];
	  _ele_L1Stage2_emul_shapeID[counter] = _Stage2_mpeg_emul_shapeID[il1];
	  _ele_L1Stage2_emul_target[counter] = ele_sclEt[counter]/(_Stage2_mpeg_emul_hwPt[il1]);

	  _ele_L1Stage2_emul_hwIsoSum[counter] = _Stage2_mpeg_emul_hwIsoSum[il1];
	  _ele_L1Stage2_emul_nTT[counter] = _Stage2_mpeg_emul_nTT[il1];
	  _ele_L1Stage2_emul_hOverERatio[counter] = _Stage2_mpeg_emul_hOverERatio[il1];
	  _ele_L1Stage2_emul_isoflag[counter] = _Stage2_mpeg_emul_isoflag[il1];
	  
	  
	}

      }      

            

    } // if iTTmax2



    /////////////////////////////////////////////////////////
    ///                 Stage 2 Level 1                   ///
    /////////////////////////////////////////////////////////

    // For DeltaR matching
    double DeltaR_min = -1;
    TLorentzVector ele_supercluster;
    ele_supercluster.SetPtEtaPhiM(ele_sclEt[counter], ele_sclEta[counter], ele_sclPhi[counter], 0);


    for(int il1=0; il1<_Stage2_eg_n; ++il1) {
      
      TLorentzVector L1Stage2_eg;
      L1Stage2_eg.SetPtEtaPhiM(_Stage2_eg_et[il1], _Stage2_eg_eta[il1], _Stage2_eg_phi[il1], 0);

      double DeltaR = ele_supercluster.DeltaR(L1Stage2_eg);

      if(DeltaR_min<0 || DeltaR<DeltaR_min){
	DeltaR_min = DeltaR;
	_ele_dR_closest_L1Stage2[counter] = DeltaR;
	_ele_closestdR_L1Stage2_eta[counter] = L1Stage2_eg.Eta();
	_ele_closestdR_L1Stage2_phi[counter] = L1Stage2_eg.Phi();
	_ele_closestdR_L1Stage2_et[counter] = L1Stage2_eg.Et();
      }

    
    }


    double DeltaR_min_emul = -1;

    for(int il1=0; il1<_Stage2_eg_emul_n; ++il1) {
      
      TLorentzVector L1Stage2_eg_emul;
      L1Stage2_eg_emul.SetPtEtaPhiM(_Stage2_eg_emul_et[il1], _Stage2_eg_emul_eta[il1], _Stage2_eg_emul_phi[il1], 0);

      double DeltaR = ele_supercluster.DeltaR(L1Stage2_eg_emul);

      if(DeltaR_min_emul<0 || DeltaR<DeltaR_min_emul){
	DeltaR_min_emul = DeltaR;
	_ele_dR_closest_L1Stage2_emul[counter] = DeltaR;
	_ele_closestdR_L1Stage2_emul_eta[counter] = L1Stage2_eg_emul.Eta();
	_ele_closestdR_L1Stage2_emul_phi[counter] = L1Stage2_eg_emul.Phi();
	_ele_closestdR_L1Stage2_emul_et[counter] = L1Stage2_eg_emul.Et();
      }


      for(int il1=0; il1<_Stage2_mpeg_emul_n; ++il1) {
      
	TLorentzVector L1Stage2_mpeg_emul;
	L1Stage2_mpeg_emul.SetPtEtaPhiM(_Stage2_mpeg_emul_et[il1], _Stage2_mpeg_emul_eta[il1], _Stage2_mpeg_emul_phi[il1], 0);
	
	double DeltaR = ele_supercluster.DeltaR(L1Stage2_mpeg_emul);
	
	if(DeltaR_min_emul<0 || DeltaR<DeltaR_min_emul){
	  DeltaR_min_emul = DeltaR;
	  _ele_closestdR_L1Stage2_mp_emul_eta[counter] = L1Stage2_mpeg_emul.Eta();
	  _ele_closestdR_L1Stage2_mp_emul_phi[counter] = L1Stage2_mpeg_emul.Phi();
	}

      }
    
    }


		
    // ---------------------------
    // For Charge, Clemy's stuff
    // ---------------------------
    if(!aod_)
      ele_expected_inner_hits[counter] = (*EleHandle)[i].gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
		
    ele_sclNclus[counter] = (*EleHandle)[i].superCluster()->clustersSize();
		
    ele_chargeGsfSC[counter]    = 0;
    ele_chargeGsfCtf[counter]   = 0;
    ele_chargeGsfCtfSC[counter] = 0;
    if((*EleHandle)[i].isGsfScPixChargeConsistent())    ele_chargeGsfSC[counter]=1;
    if((*EleHandle)[i].isGsfCtfChargeConsistent())      ele_chargeGsfCtf[counter]=1;
    if((*EleHandle)[i].isGsfCtfScPixChargeConsistent()) ele_chargeGsfCtfSC[counter]=1;
		
		
    if(!aod_)
      ele_chargeQoverPGsfVtx[counter]=float((*EleHandle)[i].gsfTrack()->charge()) / ((*EleHandle)[i].trackMomentumAtVtx().R());
    if(!(((*EleHandle)[i].closestCtfTrackRef()).isNull())) {
      ele_CtfTrackExists[counter] = 1;
      ele_chargeQoverPCtf[counter]=(*EleHandle)[i].closestCtfTrackRef()->qoverp();
    }
    else {
      ele_CtfTrackExists[counter] = 0;
    }
		
    /*if(!aod_){
      TrajectoryStateOnSurface innTSOS_ = mtsTransform_->innerStateOnSurface(*((*EleHandle)[i].gsfTrack()), 
									     *(trackerHandle_.product()), theMagField.product());
      if(innTSOS_.isValid()) {
	GlobalPoint orig(bs.position().x(), bs.position().y(), bs.position().z()) ;
	GlobalPoint scpos(sclRef->position().x(), sclRef->position().y(), sclRef->position().z()) ;
	GlobalPoint scposCorr=scpos;
	if(type_ == "DATA" && (*EleHandle)[i].isEE() && sclRef->eta() > 0) 
	  scposCorr = GlobalPoint(sclRef->position().x()+0.52, sclRef->position().y()-0.81, sclRef->position().z()+0.81) ;
	if(type_ == "DATA" && (*EleHandle)[i].isEE() && sclRef->eta() < 0) 
	  scposCorr = GlobalPoint(sclRef->position().x()-0.02, sclRef->position().y()-0.83, sclRef->position().z()-0.94) ;
	GlobalVector scvect(scpos-orig) ;
	GlobalVector scvectCorr(scposCorr-orig) ;
	GlobalPoint inntkpos = innTSOS_.globalPosition() ;
	GlobalVector inntkvect = GlobalVector(inntkpos-orig) ;
	ele_chargeDPhiInnEle[counter] = inntkvect.phi() - scvect.phi() ;
	ele_chargeDPhiInnEleCorr[counter] = inntkvect.phi() - scvectCorr.phi() ;
      } // is valid
      }*/

 
    /*
    // ------------------
    // Spikes analysis  
    // ------------------
    if ((*EleHandle)[i].isEB()) {
      const EcalRecHitCollection * reducedRecHits = 0 ;
      reducedRecHits = reducedEBRecHits.product() ; 
			
      //seed cluster analysis
      const edm::Ptr<reco::CaloCluster> & seedCluster = (*EleHandle)[i].superCluster()->seed();  
      std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(),reducedRecHits);
      const EcalRecHit & rh = getRecHit(id.first,reducedRecHits);
      int flag = rh.recoFlag();   
      if (flag == EcalRecHit::kOutOfTime) 
	ele_outOfTimeSeed[counter] = 1;   
      else 
	ele_outOfTimeSeed[counter] = 0;   
      
      int sev=sl->severityLevel(id.first,*reducedRecHits);
      ele_severityLevelSeed[counter] = sev ;
			
      int dummyFlag = 0;
      int dummySev = 0;
      for (reco::CaloCluster_iterator bc = (*EleHandle)[i].superCluster()->clustersBegin(); 
	   bc!=(*EleHandle)[i].superCluster()->clustersEnd(); 
	   ++bc) {
			  
	if ( seedCluster==(*bc) ) continue;
	std::pair<DetId, float> id = EcalClusterTools::getMaximum((*bc)->hitsAndFractions(),reducedRecHits);
	const EcalRecHit & rh = getRecHit(id.first,reducedRecHits);
	int flag = rh.recoFlag();   
	if (flag == EcalRecHit::kOutOfTime) 
	  dummyFlag = 1 ;   
	
	int sev=sl->severityLevel(id.first,*reducedRecHits);
	if (sev > dummySev)
	  dummySev = sev ;
      }			  
      ele_severityLevelClusters[counter] = dummySev ;
      ele_outOfTimeClusters[counter] = dummyFlag ;
      ele_e2overe9[counter] = E2overE9( id.first,*reducedRecHits,5,5, true, true);
 
    }
    else {
      ele_e2overe9[counter] = 0;
      ele_severityLevelSeed[counter] = 0 ;
      ele_outOfTimeSeed[counter] = 0 ;	  
      ele_severityLevelClusters[counter] = 0 ;
      ele_outOfTimeClusters[counter] = 0 ;	  
    }
		

    // Ambigous Tracks
    ele_ambiguousGsfTracks[counter] = (*EleHandle)[i].ambiguousGsfTracksSize() ;
		
    reco::GsfTrackRefVector::const_iterator firstTrack = (*EleHandle)[i].ambiguousGsfTracksBegin() ;
    reco::GsfTrackRefVector::const_iterator lastTrack  = (*EleHandle)[i].ambiguousGsfTracksEnd() ;
		
    int jj = 0 ;
    if ( ele_ambiguousGsfTracks[counter] < 5 )
      for (reco::GsfTrackRefVector::const_iterator myTrack = firstTrack ; 
	   myTrack < lastTrack ;
	   ++myTrack)
	{
				
	  ele_ambiguousGsfTracksdxy[counter][jj] = (*myTrack)->dxy() ;	
	  ele_ambiguousGsfTracksdz[counter][jj]  = (*myTrack)->dz() ;    
	  ele_ambiguousGsfTracksdxyB[counter][jj]= (*myTrack)->dxy(bs.position()) ;    
	  ele_ambiguousGsfTracksdzB[counter][jj] = (*myTrack)->dz(bs.position())  ;    
	  ++jj;
	} // for loop on tracks
		
    if(!aod_){
      edm::RefToBase<TrajectorySeed> seed = (*EleHandle)[i].gsfTrack()->extra()->seedRef();
      reco::ElectronSeedRef MyS = seed.castTo<reco::ElectronSeedRef>();
      ele_seedSubdet2[counter] = int(MyS->subDet2());
      if(fabs(MyS->dPhi2Pos()) < 100.) ele_seedDphi2Pos[counter] = double(MyS->dPhi2Pos());
      if(fabs(MyS->dRz2Pos()) < 100.)  ele_seedDrz2Pos[counter]  = double(MyS->dRz2Pos());
      if(fabs(MyS->dPhi2()) < 100.) ele_seedDphi2Neg[counter] = double(MyS->dPhi2());
      if(fabs(MyS->dRz2()) < 100.)  ele_seedDrz2Neg[counter]  = double(MyS->dRz2());
		
      ele_seedSubdet1[counter] = int(MyS->subDet1());
      if(fabs(MyS->dPhi1Pos()) < 100.) ele_seedDphi1Pos[counter] = double(MyS->dPhi1Pos());
      if(fabs(MyS->dRz1Pos()) < 100.)  ele_seedDrz1Pos[counter]  = double(MyS->dRz1Pos());
      if(fabs(MyS->dPhi1()) < 100.) ele_seedDphi1Neg[counter] = double(MyS->dPhi1());
      if(fabs(MyS->dRz1()) < 100.)  ele_seedDrz1Neg[counter]  = double(MyS->dRz1());
    }
		
    if(type_ == "MC"){
			
      double eleOkRatio = 999999.;
      double eleOkRatioE = 999999.;
      double eleOkRatioG = 999999.;
      double eleOkRatioH = 999999.;
			
      int idPDG = 999999 ;
      int idPDG_ele = 999999 ;
      int idPDG_mother_conv = 999999;
      int idPDG_pho = 999999 ;
      int idPDG_had = 999999 ;
			
      HepMC::FourVector MC_chosenEle_PoP_loc;
      HepMC::FourVector MC_chosenPho_PoP_loc;
      HepMC::FourVector MC_chosenHad_PoP_loc;
      HepMC::FourVector MC_closest_DR_loc;
			
      if(aod_ == false){
	HepMC::GenParticle* mother = 0; // for particle gun	
	for (HepMC::GenEvent::particle_const_iterator partIter = MCEvt->particles_begin();
	     partIter != MCEvt->particles_end(); ++partIter) {
	  
	  int idTmpPDG = (*partIter)->pdg_id();  
	  //MC particle
	  genPc = (*partIter);
	  pAssSim = genPc->momentum();
	  //reco electron
	  float eta = (*EleHandle)[i].eta();
	  float phi = (*EleHandle)[i].phi();
	  float p = (*EleHandle)[i].p();
	  //matching conditions
	  dphi = phi-pAssSim.phi();
	  if (fabs(dphi)>CLHEP::pi)
	    dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	  deta = eta - pAssSim.eta();
	  float deltaR = sqrt(pow(deta,2) + pow(dphi,2));
	  //standard
	  if ( deltaR < 0.15 ){                  // in the cone
	    if (  pAssSim.perp() > 1.5 ){
	      double tmpeleRatio = p/pAssSim.t();
						
	      if (idTmpPDG == 11 || idTmpPDG == -11 ){              // looking at Ele
		if ( fabs(tmpeleRatio-1) < fabs(eleOkRatioE-1) ) { //best p/p
		  eleOkRatioE = tmpeleRatio;
		  idPDG_ele = idTmpPDG;
		  
		  mother = *((*partIter)->production_vertex()->particles_begin(HepMC::parents)); 	 
		  if (mother!=0) idPDG_mother_conv = mother->pdg_id();
		  MC_chosenEle_PoP_loc = pAssSim;
		} //best p/p conditions
	      }  // looking at Ele
	      if(idTmpPDG == 22) {                                 // looking at Gamma
		if ( fabs(tmpeleRatio-1) < fabs(eleOkRatioG-1) ) {
		  eleOkRatioG = tmpeleRatio; 
		  idPDG_pho = idTmpPDG;
		  MC_chosenPho_PoP_loc = pAssSim;
		} //best p/p conditions
	      }  // looking at E/Gamma
	      if(abs(idTmpPDG) == 211 || abs(idTmpPDG) == 321){   //looking at had
		if ( fabs(tmpeleRatio-1) < fabs(eleOkRatioH-1) ) {
		  eleOkRatioH = tmpeleRatio; 
		  idPDG_had = idTmpPDG;
		  MC_chosenHad_PoP_loc = pAssSim;
		}  //best p/p
	      }  // looking at had
						
	      if ( fabs(tmpeleRatio-1) < fabs(eleOkRatio-1) ) {   // overall best p/p ratio
		eleOkRatio = tmpeleRatio; 
		idPDG = idTmpPDG;
		MC_closest_DR_loc = pAssSim;
	      }
						
	    }  // p > 1.5 
	  }  // deltaR
				
	  //	}//end loop over vertex
	  //if (okeleFound) continue ;
	}//end loop over MC particles
      } // !AOD
      else{
	const Candidate * mother = 0; //for particlegun
	for (reco::GenParticleCollection::const_iterator partIter = genParticlesColl->begin();
	     partIter != genParticlesColl->end(); ++partIter) {
				

	  int idTmpPDG = partIter->pdgId();  
	  //MC particle	  
	  pAssSim = partIter->p4();
	  //reco electron
	  float eta = (*EleHandle)[i].eta();
	  float phi = (*EleHandle)[i].phi();
	  float p = (*EleHandle)[i].p();
	  //matching conditions
	  dphi = phi-pAssSim.phi();
	  if (fabs(dphi)>CLHEP::pi)
	    dphi = dphi < 0? (CLHEP::twopi) + dphi : dphi - CLHEP::twopi;
	  deta = eta - pAssSim.eta();
	  float deltaR = sqrt(pow(deta,2) + pow(dphi,2));
	  //standard
	  if ( deltaR < 0.15 ){                  // in the cone
	    if (  pAssSim.perp() > 1.5 ){
	      double tmpeleRatio = p/pAssSim.t();
				    
	      if (idTmpPDG == 11 || idTmpPDG == -11 ){              // looking at Ele
		if ( fabs(tmpeleRatio-1) < fabs(eleOkRatioE-1) ) { //best p/p
		  eleOkRatioE = tmpeleRatio;
		  idPDG_ele = idTmpPDG;
		  //for particle gun 		  
		  mother = partIter->mother(); 	 
		  if (mother!=0) idPDG_mother_conv = mother->pdgId();
		  MC_chosenEle_PoP_loc = pAssSim;
		} //best p/p conditions
	      }  // looking at Ele
	      if(idTmpPDG == 22) {                                 // looking at Gamma
		if ( fabs(tmpeleRatio-1) < fabs(eleOkRatioG-1) ) {
		  eleOkRatioG = tmpeleRatio; 
		  idPDG_pho = idTmpPDG;
		  MC_chosenPho_PoP_loc = pAssSim;
		} //best p/p conditions
	      }  // looking at E/Gamma
	      if(abs(idTmpPDG) == 211 || abs(idTmpPDG) == 321){   //looking at had
		if ( fabs(tmpeleRatio-1) < fabs(eleOkRatioH-1) ) {
		  eleOkRatioH = tmpeleRatio; 
		  idPDG_had = idTmpPDG;
		  MC_chosenHad_PoP_loc = pAssSim;
		}  //best p/p
	      }  // looking at had
				    
	      if ( fabs(tmpeleRatio-1) < fabs(eleOkRatio-1) ) {   // overall best p/p ratio
		eleOkRatio = tmpeleRatio; 
		idPDG = idTmpPDG;
		MC_closest_DR_loc = pAssSim;
	      }
				    
	    }  // p > 1.5 
	  }  // deltaR
				

	}//end loop over MC particles

      }

      //real electrons
      if (idPDG_ele == 11 || idPDG_ele == -11) {
	ele_isMCEle[counter] = 1;   
	ele_MC_chosenEle_PoP_px[counter] = MC_chosenEle_PoP_loc.px();
	ele_MC_chosenEle_PoP_py[counter] = MC_chosenEle_PoP_loc.py();
	ele_MC_chosenEle_PoP_pz[counter] = MC_chosenEle_PoP_loc.pz();
	ele_MC_chosenEle_PoP_e[counter] = MC_chosenEle_PoP_loc.e();
	ele_idPDGmother_MCEle[counter] = idPDG_mother_conv;
      }
      //photons (or pi0)
      if(idPDG_pho == 22) {
	ele_isMCPhoton[counter] = 1;
	ele_MC_chosenPho_PoP_px[counter] = MC_chosenPho_PoP_loc.px();
	ele_MC_chosenPho_PoP_py[counter] = MC_chosenPho_PoP_loc.py();
	ele_MC_chosenPho_PoP_pz[counter] = MC_chosenPho_PoP_loc.pz();
	ele_MC_chosenPho_PoP_e[counter] = MC_chosenPho_PoP_loc.e();
      }
      //hadrons
      if(abs(idPDG_had) == 211 || abs(idPDG_had) == 321){
	ele_isMCHadron[counter] = 1; 
	ele_MC_chosenHad_PoP_px[counter] = MC_chosenHad_PoP_loc.px();
	ele_MC_chosenHad_PoP_py[counter] = MC_chosenHad_PoP_loc.py();
	ele_MC_chosenHad_PoP_pz[counter] = MC_chosenHad_PoP_loc.pz();
	ele_MC_chosenHad_PoP_e[counter] = MC_chosenHad_PoP_loc.e(); 
      }
			
      if(idPDG != 999999){
	ele_idPDGMatch[counter] = idPDG;		  
	ele_MC_closest_DR_px[counter] = MC_closest_DR_loc.px();
	ele_MC_closest_DR_py[counter] = MC_closest_DR_loc.py();
	ele_MC_closest_DR_pz[counter] = MC_closest_DR_loc.pz();
	ele_MC_closest_DR_e[counter] = MC_closest_DR_loc.e();
      }
			
      //end standard MC matching
			
      } // end if "MC"*/
		
    ++counter;  
  } // for loop on electrons
  //---------****COUNTER_ELE_FANBO
  ele_N=counter;	
}


/*

// ====================================================================================
void SimpleNtpleCustom::FillMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  // Beam spot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  const reco::BeamSpot bs = *recoBeamSpotHandle ;
  
  // Muon Retrieving
  Handle<View<reco::Muon> > MuonHandle;
  iEvent.getByLabel(MuonTag_, MuonHandle);
  
  TClonesArray &muons = *m_muons;
  int mu_counter = 0;
  
  
  // Get HZZ muon isolation
  edm::Handle<edm::ValueMap<double> > isomuonmap; 
  iEvent.getByLabel(MuonIso_HzzMapTag_, isomuonmap);
  edm::Handle<edm::ValueMap<double> > isoTkmuonmap;
  iEvent.getByLabel(MuonIsoTk_HzzMapTag_, isoTkmuonmap);
  edm::Handle<edm::ValueMap<double> > isoEcalmuonmap;
  iEvent.getByLabel(MuonIsoEcal_HzzMapTag_, isoEcalmuonmap);
  edm::Handle<edm::ValueMap<double> > isoHcalmuonmap; 
  iEvent.getByLabel(MuonIsoHcal_HzzMapTag_, isoHcalmuonmap);
  
  // ----------------------------------------------
  //  Loop over Muons
  // ----------------------------------------------
  _muons_N = MuonHandle->size();
  
  
  for (edm::View<reco::Muon>::const_iterator imuons=MuonHandle->begin(); imuons!=MuonHandle->end(); ++imuons) {  
    if(mu_counter>19) continue;
    
    // 4-vector
    setMomentum (myvector, imuons->p4());
    new (muons[mu_counter]) TLorentzVector (myvector);
    
    _muons_charge[mu_counter] = imuons->charge(); 
    
    // Provenance
    if(imuons->isTrackerMuon())    _muons_istracker[mu_counter]    = 1;
    if(imuons->isStandAloneMuon()) _muons_isstandalone[mu_counter] = 1;
    if(imuons->isGlobalMuon())     _muons_isglobal[mu_counter]     = 1;
    
    // Quality cuts
    reco::TrackRef gm = imuons->globalTrack();
    reco::TrackRef tk = imuons->innerTrack();
    
    if(imuons->isGlobalMuon()==1) {
      _muons_dxy[mu_counter]            = gm->dxy(bs.position()); //beamSpotHandle->position());
      _muons_dz[mu_counter]             = gm->dz(bs.position()); //beamSpotHandle->position());
      _muons_dxyPV[mu_counter]          = gm->dxy(math::XYZPoint(vertexPosition)); //beamSpotHandle->position());
      _muons_dzPV[mu_counter]           = gm->dz(math::XYZPoint(vertexPosition)); //beamSpotHandle->position());
      _muons_normalizedChi2[mu_counter] = gm->normalizedChi2(); 
      _muons_NmuonHits[mu_counter]      = gm->hitPattern().numberOfValidMuonHits(); // muon hit matched to global fit
    } // if Global Track
    
    if(imuons->innerTrack().isAvailable()){
      _muons_trkDxy[mu_counter]=imuons->innerTrack()->dxy();
      _muons_trkDxyError[mu_counter]=imuons->innerTrack()->dxyError();
      _muons_trkDxyB[mu_counter]=imuons->innerTrack()->dxy(bs.position()) ;
      _muons_trkDz[mu_counter]=imuons->innerTrack()->dz();
      _muons_trkDzError[mu_counter]=imuons->innerTrack()->dzError();
      _muons_trkDzB[mu_counter]=imuons->innerTrack()->dz(bs.position());
      _muons_trkChi2PerNdof[mu_counter]=imuons->innerTrack()->normalizedChi2();
      _muons_trkCharge[mu_counter]=imuons->innerTrack()->charge();
      _muons_trkNHits[mu_counter]=imuons->innerTrack()->numberOfValidHits();
      _muons_trkNPixHits[mu_counter]=imuons->innerTrack()->hitPattern().numberOfValidPixelHits();
      // Tracker muon properties
      _muons_trkmuArbitration[mu_counter]=(muon::segmentCompatibility( (*imuons),reco::Muon::SegmentAndTrackArbitration));
      _muons_trkmu2DCompatibilityLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TM2DCompatibilityLoose));
      _muons_trkmu2DCompatibilityTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TM2DCompatibilityTight));
      _muons_trkmuOneStationLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationLoose));
      _muons_trkmuOneStationTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationTight));
      _muons_trkmuLastStationLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationLoose));
      _muons_trkmuLastStationTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationTight));
      _muons_trkmuOneStationAngLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationAngLoose));
      _muons_trkmuOneStationAngTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMOneStationAngTight));
      _muons_trkmuLastStationAngLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationAngLoose));
      _muons_trkmuLastStationAngTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationAngTight));
      _muons_trkmuLastStationOptimizedLowPtLoose[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationOptimizedLowPtLoose));
      _muons_trkmuLastStationOptimizedLowPtTight[mu_counter]=(muon::isGoodMuon( (*imuons),muon::TMLastStationOptimizedLowPtTight));
    }
    
    if(imuons->isGlobalMuon()==1 || imuons->isTrackerMuon()==1) {
      _muons_NtrackerHits[mu_counter]   = tk->hitPattern().numberOfValidTrackerHits();
      _muons_NpixelHits[mu_counter]     = tk->hitPattern().numberOfValidPixelHits();
    } // if Tracker track
    _muons_Nmatches[mu_counter]             = imuons->numberOfMatches(); // number of segments matched to muon stations
    _muons_caloCompatibility[mu_counter]    = imuons->caloCompatibility() ;
    _muons_segmentCompatibility[mu_counter] = ( muon::segmentCompatibility ( (*imuons) , reco::Muon::SegmentAndTrackArbitration) ) ;
    _muons_glbmuPromptTight[mu_counter]     = ( muon::isGoodMuon( (*imuons) , muon::GlobalMuonPromptTight) );
    
    // Isolation
    _muons_nTkIsoR03[mu_counter] = imuons->isolationR03().nTracks; 
    _muons_nTkIsoR05[mu_counter] = imuons->isolationR05().nTracks;
    _muons_tkIsoR03[mu_counter]  = imuons->isolationR03().sumPt;
    _muons_tkIsoR05[mu_counter]  = imuons->isolationR05().sumPt;
    _muons_emIsoR03[mu_counter]  = imuons->isolationR03().emEt;
    _muons_emIsoR05[mu_counter]  = imuons->isolationR05().emEt;
    _muons_hadIsoR03[mu_counter] = imuons->isolationR03().hadEt;
    _muons_hadIsoR05[mu_counter] = imuons->isolationR05().hadEt;
    
    // HZZ Isolation
    edm::Ref<edm::View<reco::Muon> > muref(MuonHandle, mu_counter);
    _muons_hzzIso[mu_counter]     = (*isomuonmap)[muref]; 
    _muons_hzzIsoTk[mu_counter]   = (*isoTkmuonmap)[muref]; 
    _muons_hzzIsoEcal[mu_counter] = (*isoEcalmuonmap)[muref]; 
    _muons_hzzIsoHcal[mu_counter] = (*isoHcalmuonmap)[muref]; 
    
    mu_counter++;
  } // for loop on muons
  if(mu_counter>9) { _muons_N = 10; cout << "Number of muons>100, muons_N set to 10" << endl;}
  
} // end of FillMuons

*/

/*
// ====================================================================================
void SimpleNtpleCustom::FillJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // --------------------------------------------------
  // Calo Jets
  // --------------------------------------------------
  edm::Handle<reco::CaloJetCollection>  calojets;
  iEvent.getByLabel(CaloJetTag_, calojets);
	
  _jets_calo_N = calojets->size();
  int index_calo_jets = 0;
	
  TClonesArray &jets_calo = *_m_jets_calo; 
	
  // Loop on Calo Jets
  for ( reco::CaloJetCollection::const_iterator ijets=calojets->begin(); ijets!=calojets->end(); ijets++) {  
    if (index_calo_jets>99) continue;
		
    setMomentum (myvector, ijets->p4());
    new (jets_calo[index_calo_jets]) TLorentzVector(myvector);
		
		
    index_calo_jets++;
  } // for loop on calo jets
	
  if(index_calo_jets>99) { _jets_calo_N = 100; cout << "Number of calojets>100, RECO_CALOJETS_N set to 100" << endl;}
	
	
  // --------------------------------------------------
  // JPT Jets
  // --------------------------------------------------
  edm::Handle<reco::JPTJetCollection>  jptjets;
  iEvent.getByLabel(JPTJetTag_, jptjets);
	
  _jets_jpt_N = jptjets->size();
  int index_jpt_jets = 0;
	
  TClonesArray &jets_jpt = *_m_jets_jpt; 
	
  // Loop on Jpt Jets
  for ( reco::JPTJetCollection::const_iterator ijets=jptjets->begin(); ijets!=jptjets->end(); ijets++) {  
    if (index_jpt_jets>99) continue;
		
    setMomentum (myvector, ijets->p4());
    new (jets_jpt[index_jpt_jets]) TLorentzVector(myvector);
		
    index_jpt_jets++;
  } // for loop on jpt jets
	
  if(index_jpt_jets>99) { _jets_jpt_N = 100; cout << "Number of jptjets>100, RECO_JPTJETS_N set to 100" << endl;}
	
  // --------------------------------------------------
  // PF Jets
  // --------------------------------------------------
  edm::Handle<reco::PFJetCollection>  pfjets;
  iEvent.getByLabel(PFJetTag_, pfjets);
	
  _jets_pf_N = pfjets->size();
  int index_pf_jets = 0;
	
  TClonesArray &jets_pf = *_m_jets_pf; 
	
  // Loop on Pf Jets
  for ( reco::PFJetCollection::const_iterator ijets=pfjets->begin(); ijets!=pfjets->end(); ijets++) {  
    if (index_pf_jets>99) continue;
		
    setMomentum (myvector, ijets->p4());
    new (jets_pf[index_pf_jets]) TLorentzVector(myvector);
		
    jets_pf_chargedHadEFrac[index_pf_jets] = ijets->chargedHadronEnergyFraction ();
    jets_pf_chargedEmEFrac[index_pf_jets]  = ijets->chargedEmEnergyFraction ();
    jets_pf_chargedMuEFrac[index_pf_jets]  = ijets->chargedMuEnergyFraction ();
		
    jets_pf_neutralHadEFrac[index_pf_jets] = ijets->neutralHadronEnergyFraction ();
    jets_pf_neutralEmEFrac[index_pf_jets]  = ijets->neutralEmEnergyFraction ();
    jets_pf_PhotonEFrac[index_pf_jets]     = ijets->photonEnergyFraction();
		
    jets_pf_chargedHadMultiplicity[index_pf_jets] = ijets->chargedHadronMultiplicity ();
    jets_pf_neutralHadMultiplicity[index_pf_jets] = ijets->neutralHadronMultiplicity ();
		
    jets_pf_chargedMultiplicity[index_pf_jets] = ijets->chargedMultiplicity ();
    jets_pf_neutralMultiplicity[index_pf_jets] = ijets->neutralMultiplicity ();
		
    jets_pf_nConstituents[index_pf_jets]       = ijets->nConstituents();
		
    index_pf_jets++;
  } // for loop on pf jets
	
  if(index_pf_jets>99) { _jets_pf_N = 100; cout << "Number of pfjets>100, RECO_PFJETS_N set to 100" << endl;}
	
} // end of FillJets
*/
 /*
// ====================================================================================
void SimpleNtpleCustom::FillSuperClusters(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
 
  // geometry
  
  unsigned long long cacheIDGeom = 0;
  edm::ESHandle<CaloGeometry> theCaloGeom;
  if(cacheIDGeom!=iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
    cacheIDGeom = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
    iSetup.get<CaloGeometryRecord>().get(theCaloGeom_);
  }
  
  
  // TPGTowerStatus
  //edm::ESHandle<EcalTPGTowerStatus> theEcalTPGTowerStatus_handle;
  //iSetup.get<EcalTPGTowerStatusRcd>().get(theEcalTPGTowerStatus_handle);
  //const EcalTPGTowerStatus * ecaltpgTowerStatus=theEcalTPGTowerStatus_handle.product();
  
  //const EcalTPGTowerStatusMap &towerMap=ecaltpgTowerStatus->getMap();
  //EcalTPGTowerStatusMapIterator  it;
  
  //for42x	
  unsigned long long cacheSevLevel = 0;
  edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
  if(cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()){
    cacheSevLevel = iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
  }
  //const EcalSeverityLevelAlgo* sl=sevLevel.product();

  // for H/E  
  towersH_ = new edm::Handle<CaloTowerCollection>() ;
  if (!iEvent.getByLabel(hcalTowers_,*towersH_))
    { edm::LogError("ElectronHcalHelper::readEvent")<<"failed to get the hcal towers of label "<<hcalTowers_ ; }
  towerIso1_ = new EgammaTowerIsolation(hOverEConeSize_,0.,hOverEPtMin_,1,towersH_->product()) ;
  towerIso2_ = new EgammaTowerIsolation(hOverEConeSize_,0.,hOverEPtMin_,2,towersH_->product()) ;
  

  
  // Beam Spot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  //const reco::BeamSpot bs = *recoBeamSpotHandle ;


  // Isolation
  edm::Handle<TrackCollection> ctfTracksH;  
  iEvent.getByLabel("generalTracks", ctfTracksH); // ctfTracks_

  
 
	
} // FillSuperCluster 

*/
/*
// ====================================================================================
void SimpleNtpleCustom::FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  // get gen particle candidates
  edm::Handle<GenParticleCollection> genCandidatesCollection;
  iEvent.getByLabel("genParticles", genCandidatesCollection);
  
  TClonesArray &MC_gen_V         = *_m_MC_gen_V;
  TClonesArray &MC_gen_leptons   = *_m_MC_gen_leptons;
  
  int counter = 0;
  int counter_daughters = 0;
  
  // ----------------------------
  //      Loop on particles
  // ----------------------------
  for( GenParticleCollection::const_iterator p = genCandidatesCollection->begin();p != genCandidatesCollection->end(); ++ p ) {
    
    if (p->pdgId() == 23 || fabs(p->pdgId())==24) {
      if(p->status()==3) {
	// Fill truth W,Z
	setMomentum (myvector,p->p4());
	new (MC_gen_V[counter]) TLorentzVector(myvector);
	_MC_gen_V_pdgid[counter] = p->pdgId();
	
	//size_t nfirstdau = p->numberOfDaughters();
	
	// Loop on daughters
	for(unsigned int i=0;i<p->numberOfDaughters();i++) {
	  bool islep = false;
	  if(fabs(p->daughter(i)->pdgId())==11) { _MC_flavor[counter] = 0; islep=true;} // electron
	  if(fabs(p->daughter(i)->pdgId())==13) { _MC_flavor[counter] = 1; islep=true;} // muon
	  if(fabs(p->daughter(i)->pdgId())==15) { _MC_flavor[counter] = 2; islep=true;} // taus
	  
	  if(islep) { // p->daughter(i)->status()==1) { ?!
	    setMomentum(myvector, p->daughter(i)->p4());
	    new (MC_gen_leptons[counter_daughters]) TLorentzVector(myvector);
	    _MC_gen_leptons_pdgid[counter_daughters] = p->daughter(i)->pdgId();
	    
	    counter_daughters++;
	  } // if is lepton
	} // for loop on daughters
	counter++;
      } // if status stable
    } // if W or Z
    
  } // for loop on particles
  
} // end of FillTruth
*/
/*
// ====================================================================================
void  SimpleNtpleCustom::FillTipLipIp(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  //Get the B-field
  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  
  //Get Beam Spot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  const reco::BeamSpot bs = *recoBeamSpotHandle ;
  
  GlobalPoint BSVertex(bs.position().x(),bs.position().y(),bs.position().z());
  Basic3DVector<double> BSVertexErr(bs.x0Error(),bs.y0Error(),bs.z0Error());
  
  reco::Vertex::Point BSp(bs.position().x(),bs.position().y(),bs.position().z());
  reco::Vertex::Error BSe;
  
  BSe(0,0) = bs.x0Error()*bs.x0Error();
  BSe(1,1) = bs.y0Error()*bs.y0Error();
  BSe(2,2) = bs.z0Error()*bs.z0Error();
  reco::Vertex BSprimVertex = reco::Vertex(BSp,BSe,1,1,1);
  
  // get the track builder
  ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);
  
  //Get Vertices
  edm::Handle<reco::VertexCollection> recoPrimaryVertexCollection;
  iEvent.getByLabel(VerticesTag_,recoPrimaryVertexCollection);
  reco::Vertex primVertex;
  bool pvfound = (recoPrimaryVertexCollection->size() != 0);
  if (pvfound)
    {
      PrimaryVertexSorter pvs;
      vector<reco::Vertex> sortedList = pvs.sortedList(*(recoPrimaryVertexCollection.product()) );
      primVertex = (sortedList.front());
    } else {
    //creating a dummy PV
    reco::Vertex::Point p(0,0,0);          
    reco::Vertex::Error e;
    e(0,0) = 0.0015*0.0015;
    e(1,1) = 0.0015*0.0015;
    e(2,2) = 15.*15.;
    primVertex = reco::Vertex(p,e,1,1,1);
  }
  //
  GlobalPoint pVertex(primVertex.position().x(),primVertex.position().y(),primVertex.position().z());
  Basic3DVector<double> pVertexErr(primVertex.xError(),primVertex.yError(),primVertex.zError());
  
  //--------------- for propagation
  // electrons:
  const GsfPropagatorAdapter *theGeomPropBw = new GsfPropagatorAdapter(AnalyticalPropagator(magneticField.product(),oppositeToMomentum)); 
  // muons:
  Propagator *thePropBw = new AnalyticalPropagator(magneticField.product(),oppositeToMomentum);                    
  

  float muTip,muLip,muSTip,muSLip,muTipSignif,muLipSignif,muSignificance3D,muValue3D,muError3D ;
  float eleTip,eleLip,eleSTip,eleSLip,eleTipSignif,eleLipSignif,eleSignificance3D,eleValue3D,eleError3D;
  //
  //
  //--track refs
  TrackRef mutrack;
  GsfTrackRef eletrack;
  //--transient tracks
  reco::TransientTrack mutranstrack;
  reco::TransientTrack eletranstrack;

  // =============================================================================
  // Muons
  // =============================================================================
  Handle<View<reco::Muon> > MuonHandle;
  iEvent.getByLabel(MuonTag_, MuonHandle);
  
  unsigned int indexmu=0;
  for (edm::View<reco::Muon>::const_iterator muCand = MuonHandle->begin(); muCand != MuonHandle->end(); ++muCand){
    
    if (indexmu>19) continue;
    mutrack = muCand->get<TrackRef>();
    if (mutrack.isNull()){
      cout <<"tracker trackref is null since muon is STA" << endl;
      mutrack=muCand->get<TrackRef,reco::StandAloneMuonTag>();
    }  

    mutranstrack = trackBuilder->build( mutrack );
    
    TrajectoryStateOnSurface innerMuTSOS;

    if (useBeamSpot_==true){ 
      innerMuTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), BSVertex, mutranstrack.field());
    } 
    else {
      innerMuTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), pVertex, mutranstrack.field());
    } 


    
    if (innerMuTSOS.isValid() && !mutrack.isNull() ){    
      
      //-- get propagated the inner TSOS to the PV:
      TrajectoryStateOnSurface vtxMuTSOS;
      if (useBeamSpot_==true){ 
	vtxMuTSOS = TransverseImpactPointExtrapolator(*thePropBw).extrapolate(innerMuTSOS,BSVertex);
      } 
      else {
	vtxMuTSOS = TransverseImpactPointExtrapolator(*thePropBw).extrapolate(innerMuTSOS,pVertex);
      }
      
      //		     
      if (!vtxMuTSOS.isValid()){		 
	vtxMuTSOS = innerMuTSOS; //-protection for eventual failing extrapolation
      }
      
      //-- get the distances (transverse & longitudinal) between extrapolated position and PV position
      GlobalPoint mimpP = vtxMuTSOS.globalPosition();
      GlobalVector mdistV; 
      GlobalVector direction=vtxMuTSOS.globalDirection();
      
      if (useBeamSpot_==true){ 
	mdistV = mimpP - BSVertex; 
      } 
      else {
	mdistV = mimpP - pVertex; 
      }
      
      GlobalVector transverseIP(mdistV.x(),mdistV.y(),0.); 
      double ps = transverseIP.dot(direction);
      muTip = mdistV.perp()*((ps!=0)?ps/abs(ps):1.); //signed by definition
      muLip = mdistV.z();    // signed by definition
      
      // compute full error calculation:
      // - diagonal terms first:
      AlgebraicSymMatrix33 mvtxerrM; 
      if (useBeamSpot_==true){ 
	mvtxerrM(0,0) = BSVertexErr.x()*BSVertexErr.x(); 
	mvtxerrM(1,1) = BSVertexErr.y()*BSVertexErr.y();
	mvtxerrM(2,2) = BSVertexErr.z()*BSVertexErr.z();
      } 
      else {
	mvtxerrM(0,0) = pVertexErr.x()*pVertexErr.x(); 
	mvtxerrM(1,1) = pVertexErr.y()*pVertexErr.y();
	mvtxerrM(2,2) = pVertexErr.z()*pVertexErr.z();
      }
      
      // - off-diagonal terms:
      AlgebraicSymMatrix33 merrorM = mvtxerrM + vtxMuTSOS.cartesianError().matrix().Sub<AlgebraicSymMatrix33>(0,0);
      AlgebraicVector2 mjacobianTranV;	
      AlgebraicVector1 mjacobianLongV;
      mjacobianTranV[0] = mdistV.x()/mdistV.perp();	
      mjacobianTranV[1] = mdistV.y()/mdistV.perp();
      mjacobianLongV[0] = 1.;	
      //- errors:
      muSTip = sqrt(ROOT::Math::Similarity(merrorM.Sub<AlgebraicSymMatrix22>(0,0),mjacobianTranV));
      muSLip = sqrt(ROOT::Math::Similarity(merrorM.Sub<AlgebraicSymMatrix11>(2,2),mjacobianLongV));
      
      //
      muTipSignif=muTip/muSTip;
      muLipSignif=muLip/muSLip;
      muons_Tip[indexmu] = muTip ;
      muons_Lip[indexmu] = muLip ;
      muons_STip[indexmu] = muSTip ;
      muons_SLip[indexmu] = muSLip ;
      muons_TipSignif[indexmu] = muTipSignif ;
      muons_LipSignif[indexmu] = muLipSignif ;
      
    }
    
 

    // -------------
    //  3DIP & SIP 
    // -------------
    
    TrajectoryStateOnSurface muTSOS;
    if (useBeamSpot_==true){ 
      muTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), BSVertex, mutranstrack.field());
    } 
    else {
      muTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), pVertex, mutranstrack.field());
    }
    
    if (muTSOS.isValid()){
      std::pair<bool,Measurement1D> muIPpair;
      
      if (useBeamSpot_==true){ 		  
	muIPpair = IPTools::absoluteImpactParameter3D(mutranstrack, BSprimVertex);
      } 
      else {	 
	muIPpair = IPTools::absoluteImpactParameter3D(mutranstrack, primVertex);
      }
      
      if (muIPpair.first){
	muSignificance3D = muIPpair.second.significance();
	muValue3D = muIPpair.second.value();
	muError3D = muIPpair.second.error();
	muons_Significance3D[indexmu] = muSignificance3D ;
	muons_Value3D[indexmu] = muValue3D ;
	muons_Error3D[indexmu] = muError3D ;
	
      } 	       
    }
    ++indexmu;
  } //-- muon loop closed
  
  // =============================================================================
  // Electrons
  // =============================================================================
    Handle<edm::View<GsfElectron> > eleCandidates;
    iEvent.getByLabel(EleTag_.label(),eleCandidates);
    unsigned int indexele=0;
    for (edm::View<reco::GsfElectron>::const_iterator eleCand = eleCandidates->begin(); eleCand != eleCandidates->end(); ++eleCand){
      if (indexele>9) continue;
      
      eletrack = eleCand->get<GsfTrackRef>();
      eletranstrack = trackBuilder->build( eletrack ) ;
      //         
      // get initial TSOS:
      TrajectoryStateOnSurface innerTSOS;
      if (useBeamSpot_==true){
	innerTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), BSVertex, eletranstrack.field());
      } 
      else {
	innerTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
      }


      if (innerTSOS.isValid()){
	//-- get propagated the inner TSOS to the PV:
	TrajectoryStateOnSurface vtxTSOS;
	
	if (useBeamSpot_==true){ 
	  vtxTSOS = TransverseImpactPointExtrapolator(*theGeomPropBw).extrapolate(innerTSOS,BSVertex);
	} 
	else {
	  vtxTSOS = TransverseImpactPointExtrapolator(*theGeomPropBw).extrapolate(innerTSOS,pVertex);
	}
	
	//		     
	if (!vtxTSOS.isValid()){		 
	  vtxTSOS = innerTSOS; //-protection for eventual failing extrapolation
	}       
	//
	//-- get the distances (transverse & longitudinal) between extrapolated position and PV position
	GlobalPoint impP = vtxTSOS.globalPosition();
	GlobalVector distV;
	GlobalVector direction=vtxTSOS.globalDirection();
	
	if (useBeamSpot_==true){ 
	  distV = impP - BSVertex; 	
	} 
	else {
	  distV = impP - pVertex; 
	}
	
	GlobalVector transverseIPele(distV.x(),distV.y(),0.); 
	double psele = transverseIPele.dot(direction);
	eleTip = distV.perp()*((psele!=0)?psele/abs(psele):1.); // signed by definition
	eleLip = distV.z();    // signed by definition
	
	// compute full error calculation:
	// - diagonal terms first:
	AlgebraicSymMatrix33 vtxerrM; 
	if (useBeamSpot_==true){ 
	  vtxerrM(0,0) = BSVertexErr.x()*BSVertexErr.x(); 
	  vtxerrM(1,1) = BSVertexErr.y()*BSVertexErr.y();
	  vtxerrM(2,2) = BSVertexErr.z()*BSVertexErr.z();
	} 
	else {
	  vtxerrM(0,0) = pVertexErr.x()*pVertexErr.x(); 
	  vtxerrM(1,1) = pVertexErr.y()*pVertexErr.y();
	  vtxerrM(2,2) = pVertexErr.z()*pVertexErr.z();
	}
	
	// - off-diagonal terms:
	AlgebraicSymMatrix33 errorM = vtxerrM + vtxTSOS.cartesianError().matrix().Sub<AlgebraicSymMatrix33>(0,0);
	AlgebraicVector2 jacobianTranV;	
	AlgebraicVector1 jacobianLongV;
	jacobianTranV[0] = distV.x()/distV.perp();	
	jacobianTranV[1] = distV.y()/distV.perp();
	jacobianLongV[0] = 1.;	
	//- errors:
	eleSTip = sqrt(ROOT::Math::Similarity(errorM.Sub<AlgebraicSymMatrix22>(0,0),jacobianTranV));
	eleSLip = sqrt(ROOT::Math::Similarity(errorM.Sub<AlgebraicSymMatrix11>(2,2),jacobianLongV));
	
	eleTipSignif=eleTip/eleSTip;
	eleLipSignif=eleLip/eleSLip;
	ele_Tip[indexele] = eleTip ;
	ele_Lip[indexele] = eleLip ;
	ele_STip[indexele] = eleSTip ;
	ele_SLip[indexele] = eleSLip ;
	ele_TipSignif[indexele] = eleTipSignif ;
	ele_LipSignif[indexele] = eleLipSignif ;
	
      }

      // -----------------
      //   3DIP & SIP
      // -----------------
      
      TrajectoryStateOnSurface eleTSOS;
      if (useBeamSpot_==true){ 
	eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), BSVertex, eletranstrack.field());
      } 
      else {
	eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
      }
      
      if (eleTSOS.isValid()){
	std::pair<bool,Measurement1D> eleIPpair;
	if (useBeamSpot_==true){ 
	  eleIPpair = IPTools::signedImpactParameter3D(eletranstrack, eleTSOS.globalDirection(), BSprimVertex); 
	}
	else {
	  eleIPpair = IPTools::signedImpactParameter3D(eletranstrack, eleTSOS.globalDirection(), primVertex);
	}
	
	if (eleIPpair.first){
	  eleSignificance3D = eleIPpair.second.significance();
	  eleValue3D = eleIPpair.second.value();
	  eleError3D = eleIPpair.second.error();
	  ele_Significance3D[indexele] = eleSignificance3D ;
	  ele_Value3D[indexele] = eleValue3D ;
	  ele_Error3D[indexele] = eleError3D ;
	  
	} 	
      }            
      //
      ++indexele;
    } //-- ele loop closed
    
	
}*/

// ====================================================================================
void SimpleNtpleCustom::Init()
// ====================================================================================
{
	
  ele_N = 0;
  ele_nSeed = 0;
  nEvent = 0;
  nRun = 0;
  nLumi = 0;

  //Pile-up
  _PU_N     = 0;
  _PU_rho   = 0.;
  _PU_sigma = 0.;
	
	
  // Vertices
  _vtx_N = 0; 
  for(int iv=0;iv<15;iv++) {
    _vtx_normalizedChi2[iv] = 0.;
    _vtx_ndof[iv] = 0.;
    _vtx_nTracks[iv] = 0.;
    _vtx_d0[iv] = 0.;
    _vtx_x[iv] = 0.;
    _vtx_y[iv] = 0.;
    _vtx_z[iv] = 0.;
  }// for loop on vertices
	
  // Beam Spot
  BS_x = 0.;
  BS_y = 0.;
  BS_z = 0.;
	
  BS_dz = 0.;
  BS_dxdz = 0.;
  BS_dydz = 0.;
	
  BS_bw_x = 0.;
  BS_bw_y = 0.;
	
  // MC truth
  _MC_pthat  = 0.;
  _MC_flavor[0] = 10;
  _MC_flavor[1] = 10;
	
  // Trigger towers
  
  _trig_tower_N = 0;
  _trig_tower_N_modif = 0;
  _trig_tower_N_emul = 0;

  for(int i=0 ; i<4032 ; i++) {
    _trig_tower_ieta[i]=-999;
    _trig_tower_iphi[i]=-999;
    _trig_tower_adc[i]=-999;
    _trig_tower_sFGVB[i]=-999;

    _trig_tower_ieta_modif[i]=-999;
    _trig_tower_iphi_modif[i]=-999;
    _trig_tower_adc_modif[i]=-999;
    _trig_tower_sFGVB_modif[i]=-999;

    _trig_tower_ieta_emul[i]=-999;
    _trig_tower_iphi_emul[i]=-999;
    for(int j=0 ; j<5 ; j++) {
      _trig_tower_adc_emul[i][j]=-999;
      _trig_tower_sFGVB_emul[i][j]=-999;
    }
  }

  // Trigger
  for (int i=0;i<600;i++) 
    trig_hltInfo[i] = 0;
	
  trig_isUnbiased = 0 ;
	
  trig_isPhoton10 = 0; 
  trig_isPhoton15 = 0; 
  trig_isL1SingleEG2 = 0; 
  trig_isL1SingleEG5 = 0;
  trig_isL1SingleEG8 = 0;
  trig_isEle10_LW = 0;
  trig_isEle15_LW = 0;
  trig_isHLT_Ele30WP60_Ele8_Mass55_v2 = 0;
  trig_isHLT_Ele30WP60_SC4_Mass55_v3 = 0;
  trig_isHLT_Ele27_WPLoose_Gsf_v1 = 0;
	
  _trig_isEleHLTpath = 0;
  _trig_isMuonHLTpath = 0;
	
  // L1
  _trig_L1emIso_N    = 0; 
  _trig_L1emNonIso_N = 0;
  _trig_L1emIso_N_M    = 0; 
  _trig_L1emNonIso_N_M = 0;
  _trig_preL1emIso_N    = 0; 
  _trig_preL1emNonIso_N = 0;
  _trig_postL1emIso_N    = 0; 
  _trig_postL1emNonIso_N = 0;
  // max set to 4
  for(int il1=0;il1<4;il1++) {
    // Used by Clemy
    _trig_L1emIso_ieta[il1] = 0; 
    _trig_L1emIso_iphi[il1] = 0; 
    _trig_L1emIso_rank[il1] = 0; 
    _trig_L1emIso_ieta_M[il1] = 0; 
    _trig_L1emIso_iphi_M[il1] = 0; 
    _trig_L1emIso_rank_M[il1] = 0; 
     // From Trigger twiki
    _trig_L1emIso_eta[il1]    = 0.; 
    _trig_L1emIso_phi[il1]   = 0.; 
    _trig_L1emIso_energy[il1] = 0.; 
    _trig_L1emIso_et[il1]     = 0.; 
    _trig_L1emIso_eta_M[il1]    = 0.; 
    _trig_L1emIso_phi_M[il1]   = 0.; 
    _trig_L1emIso_energy_M[il1] = 0.; 
    _trig_L1emIso_et_M[il1]     = 0.; 
		
    // Used by Clemy
    _trig_L1emNonIso_ieta[il1] = 0; 
    _trig_L1emNonIso_iphi[il1] = 0; 
    _trig_L1emNonIso_rank[il1] = 0; 
    _trig_L1emNonIso_ieta_M[il1] = 0; 
    _trig_L1emNonIso_iphi_M[il1] = 0; 
    _trig_L1emNonIso_rank_M[il1] = 0; 
    // From Trigger twiki
    _trig_L1emNonIso_eta[il1]    = 0.; 
    _trig_L1emNonIso_phi[il1]   = 0.; 
    _trig_L1emNonIso_energy[il1] = 0.; 
    _trig_L1emNonIso_et[il1]     = 0.; 
    _trig_L1emNonIso_eta_M[il1]    = 0.; 
    _trig_L1emNonIso_phi_M[il1]   = 0.; 
    _trig_L1emNonIso_energy_M[il1] = 0.; 
    _trig_L1emNonIso_et_M[il1]     = 0.; 
		
    // Used by Clemy
    _trig_preL1emIso_ieta[il1] = 0; 
    _trig_preL1emIso_iphi[il1] = 0; 
    _trig_preL1emIso_rank[il1] = 0;
    // Used by Clemy
    _trig_preL1emNonIso_ieta[il1] = 0; 
    _trig_preL1emNonIso_iphi[il1] = 0; 
    _trig_preL1emNonIso_rank[il1] = 0; 
		
    // Used by Clemy
    _trig_postL1emIso_ieta[il1] = 0; 
    _trig_postL1emIso_iphi[il1] = 0; 
    _trig_postL1emIso_rank[il1] = 0;
    // Used by Clemy
    _trig_postL1emNonIso_ieta[il1] = 0; 
    _trig_postL1emNonIso_iphi[il1] = 0; 
    _trig_postL1emNonIso_rank[il1] = 0;  
		
		
  } // for loop on L1 cand
	
  // HLT
  _trig_HLT_N = 0;
  for(int ihlt=0;ihlt<20;ihlt++) {
    _trig_HLT_eta[ihlt]    = 0.; 
    _trig_HLT_phi[ihlt]    = 0.; 
    _trig_HLT_energy[ihlt] = 0.; 
    _trig_HLT_pt[ihlt]     = 0.;
    _trig_HLT_name[ihlt]   = -1;
  } // for loop on hlt
	


  // Stage 2 Level 1

  _Stage2_tower_n = 0;
  for(int i_tow=0; i_tow<5760;i_tow++){
    _Stage2_tower_hwPt[i_tow] = 0;
    _Stage2_tower_hwEta[i_tow] = 0;
    _Stage2_tower_hwPhi[i_tow] = 0;
  }


  _Stage2_mpeg_n = 0;
  for(int i_eg=0; i_eg<12;i_eg++){
    _Stage2_mpeg_hwPt[i_eg] = 0;
    _Stage2_mpeg_hwEta[i_eg] = 0;
    _Stage2_mpeg_hwPhi[i_eg] = 0;
  }

  _Stage2_eg_n = 0;
  for(int i_eg=0; i_eg<12;i_eg++){
    _Stage2_eg_hwPt[i_eg] = 0;
    _Stage2_eg_hwEta[i_eg] = 0;
    _Stage2_eg_hwPhi[i_eg] = 0;
    _Stage2_eg_et[i_eg] = 0;
    _Stage2_eg_eta[i_eg] = 0;
    _Stage2_eg_phi[i_eg] = 0;
    _Stage2_eg_isoflag[i_eg] = 0;
  }

  
  _Stage2_tower_emul_n = 0;
  for(int i_tow=0; i_tow<5760;i_tow++){
    _Stage2_tower_emul_hwPt[i_tow] = 0;
    _Stage2_tower_emul_hwEtEm[i_tow] = 0;
    _Stage2_tower_emul_hwEtHad[i_tow] = 0;
    _Stage2_tower_emul_hwEta[i_tow] = 0;
    _Stage2_tower_emul_hwPhi[i_tow] = 0;
  }

  _Stage2_mpeg_emul_n = 0;
  for(int i_eg=0; i_eg<12;i_eg++){
    _Stage2_mpeg_emul_hwPt[i_eg] = 0;
    _Stage2_mpeg_emul_hwEta[i_eg] = 0;
    _Stage2_mpeg_emul_hwPhi[i_eg] = 0;
    _Stage2_mpeg_emul_et[i_eg] = 0;
    _Stage2_mpeg_emul_eta[i_eg] = 0;
    _Stage2_mpeg_emul_phi[i_eg] = 0;
    _Stage2_mpeg_emul_shape[i_eg] = 0;
    _Stage2_mpeg_emul_shapeID[i_eg] = 0;
    _Stage2_mpeg_emul_hwIsoSum[i_eg] = 0;
    _Stage2_mpeg_emul_nTT[i_eg] = 0;
    _Stage2_mpeg_emul_hOverERatio[i_eg] = 0;
    _Stage2_mpeg_emul_isoflag[i_eg] = 0;
  }

  _Stage2_eg_emul_n = 0;
  for(int i_eg=0; i_eg<12;i_eg++){
    _Stage2_eg_emul_hwPt[i_eg] = 0;
    _Stage2_eg_emul_hwEta[i_eg] = 0;
    _Stage2_eg_emul_hwPhi[i_eg] = 0;
    _Stage2_eg_emul_et[i_eg] = 0;
    _Stage2_eg_emul_eta[i_eg] = 0;
    _Stage2_eg_emul_phi[i_eg] = 0;
  }

  
  _Stage1_eg_emul_n = 0;
  for(int i_eg=0; i_eg<12;i_eg++){
    _Stage1_eg_emul_hwPt[i_eg] = 0;
    _Stage1_eg_emul_hwEta[i_eg] = 0;
    _Stage1_eg_emul_hwPhi[i_eg] = 0;
    _Stage1_eg_emul_et[i_eg] = 0;
    _Stage1_eg_emul_eta[i_eg] = 0;
    _Stage1_eg_emul_phi[i_eg] = 0;
    _Stage1_eg_emul_isoflag[i_eg] = 0;
  }

	
  // Masked Towers
  _trig_nMaskedRCT=0;
  _trig_nMaskedCh=0;
	
  for (int ii=0;ii<100;ii++)
    {
      _trig_iMaskedRCTeta[ii]   = -999;
      _trig_iMaskedRCTphi[ii]   = -999;
      _trig_iMaskedRCTcrate[ii] = -999;
      _trig_iMaskedTTeta[ii]    = -999;
      _trig_iMaskedTTphi[ii]    = -999;
    }//loop  masks
	
  for (int ii=0;ii<10;ii++)
    {
      _ele_RCTeta[ii]      = -999;
      _ele_RCTphi[ii]      = -999;
      _ele_RCTL1iso[ii]    = -999;
      _ele_RCTL1noniso[ii] = -999;
      _ele_RCTL1iso_M[ii]    = -999;
      _ele_RCTL1noniso_M[ii] = -999;
   }//loop electron rct region
	

  for (int ii=0;ii<10;ii++)
    {
      _ele_dR_closest_L1Stage2[ii]      = -999;
      _ele_closestdR_L1Stage2_eta[ii]   = -999;
      _ele_closestdR_L1Stage2_phi[ii]   = -999;
      _ele_closestdR_L1Stage2_et[ii]    = -999;
      _ele_L1Stage2[ii]                 = -999;
      _ele_L1Stage2_isoflag[ii]         = -999;

      _ele_dR_closest_L1Stage2_emul[ii]      = -999;
      _ele_closestdR_L1Stage2_emul_eta[ii]   = -999;
      _ele_closestdR_L1Stage2_emul_phi[ii]   = -999;
      _ele_closestdR_L1Stage2_mp_emul_eta[ii]   = -999;
      _ele_closestdR_L1Stage2_mp_emul_phi[ii]   = -999;
      _ele_closestdR_L1Stage2_emul_et[ii]    = -999;
      _ele_L1Stage2_emul[ii]                 = -999;

      _ele_L1Stage2_emul_ieta[ii]            = -999;
      _ele_L1Stage2_emul_hwPt[ii]            = -999;
      _ele_L1Stage2_emul_shape[ii]           = -999;
      _ele_L1Stage2_emul_shapeID[ii]         = -999;
      _ele_L1Stage2_emul_target[ii]          = -999;
      _ele_L1Stage2_emul_hwIsoSum[ii]        = -999;
      _ele_L1Stage2_emul_nTT[ii]             = -999;
      _ele_L1Stage2_emul_hOverERatio[ii]     = -999;
      _ele_L1Stage2_emul_isoflag[ii]         = -999;

      _ele_L1Stage1_emul[ii]                 = -999;
      _ele_L1Stage1_emul_isoflag[ii]         = -999;
      _ele_L1Stage1_emul_eta[ii]             = -999;
      _ele_L1Stage1_emul_phi[ii]             = -999;



   }//loop Stage2



  // Offline Electrons
  for (int i=0;i<10;i++) 
    {
      ele_MC_chosenEle_PoP_px[i] = 0.;
      ele_MC_chosenEle_PoP_py[i] = 0.;
      ele_MC_chosenEle_PoP_pz[i] = 0.;
      ele_MC_chosenEle_PoP_e[i] = 0.;
		
      ele_MC_chosenPho_PoP_px[i] = 0.;
      ele_MC_chosenPho_PoP_py[i] = 0.;
      ele_MC_chosenPho_PoP_pz[i] = 0.;
      ele_MC_chosenPho_PoP_e[i] = 0.;
		
      ele_MC_chosenHad_PoP_px[i] = 0.;
      ele_MC_chosenHad_PoP_py[i] = 0.;
      ele_MC_chosenHad_PoP_pz[i] = 0.;
      ele_MC_chosenHad_PoP_e[i] = 0.;
		
      ele_MC_closest_DR_px[i] = 0.;
      ele_MC_closest_DR_py[i] = 0.;
      ele_MC_closest_DR_pz[i] = 0.;
      ele_MC_closest_DR_e[i] = 0.;

      _ele_he_00615_0[i]  = 0.; // 
      _ele_he_005_0[i]  = 0.; //
      _ele_he_005_1[i]  = 0.; //HoE_005_1 ;	
      _ele_he_005_15[i] = 0.; //HoE_005_15 ;
      _ele_he_01_0[i]   = 0.; //HoE_01_0 ;
      _ele_he_01_1[i]   = 0.; //HoE_01_1 ;	
      _ele_he_01_15[i]  = 0.; //HoE_01_15 ;
      //_ele_he_015_0[i]  0.; //= HoE_01_0 ;
      _ele_he_015_1[i]  = 0.; //HoE_015_1 ;	
      _ele_he_015_15[i] = 0.; //HoE_015_15 ;
		
      ele_VetoIdDecisions[i] = 0.;
      ele_LooseIdDecisions[i] = 0.;
      ele_MediumIdDecisions[i] = 0.;
      ele_TightIdDecisions[i] = 0.;
		
      ele_pT[i]=0;
      ele_echarge[i]=0;
      ele_he[i]=0; 
      ele_eseedpout[i]=0;
      ele_ep[i]=0;
      ele_eseedp[i]=0;
      ele_eelepout[i]=0;  
      ele_pin_mode[i]=0;
      ele_pout_mode[i]=0;
      ele_pin_mean[i]=0;
      ele_pout_mean[i]=0; 
      ele_calo_energy[i]=0;
      ele_sclRawE[i]=0;
      ele_sclE[i]=0;
      ele_sclEt[i]=0;
      ele_sclEta[i]=0;
      ele_sclPhi[i]=0;
      ele_sclX[i]=0;
      ele_sclY[i]=0;
      ele_sclZ[i]=0;
      ele_sclErr[i]=0;
      ele_sclErr_pos[i]=0;
      ele_sclErr_neg[i]=0;
      ele_trErr[i]=0;
      ele_momErr[i]=0;
      ele_newmomErr[i]=0;
      ele_newmom[i]=0;
      ele_tr_atcaloX[i]=0;
      ele_tr_atcaloY[i]=0;
      ele_tr_atcaloZ[i]=0;
      ele_firsthit_X[i]=0;
      ele_firsthit_Y[i]=0;
      ele_firsthit_Z[i]=0;
		
      ele_pTin_mode[i]=0;
      ele_pTout_mode[i]=0; 
      ele_pTin_mean[i]=0; ; 
      ele_pTout_mean[i]=0;
      ele_deltaetaseed[i]=0;
      ele_deltaetaele[i]=0;
      ele_deltaphiseed[i]=0;
      ele_deltaphiele[i]=0;
      ele_deltaetain[i]=0;
      ele_deltaphiin[i]=0;
      ele_sigmaietaieta[i]=0;
      ele_sigmaetaeta[i]=0;
      ele_e15[i]=0;
      ele_e25max[i]=0;
      ele_e55[i]=0;
      ele_e1[i]=0;
      ele_e33[i]=0;
      ele_e2overe9[i]=0;
      ele_fbrem[i]=0 ;
      ele_mva[i]=0 ;
      ele_isbarrel[i]=0;
      ele_isendcap[i]=0;
      ele_isEBetaGap[i]=0;
      ele_isEBphiGap[i]=0;
      ele_isEEdeeGap[i]=0;
      ele_isEEringGap[i]=0;
      ele_isecalDriven[i]=0;
      ele_istrackerDriven[i]=0;
      ele_eClass[i]=0;
      ele_missing_hits[i]=0;
      ele_dxy[i]=0 ;
      ele_dz[i]=0 ;
      ele_dsz[i]=0.;
      ele_dxyB[i]=0 ;
      ele_dzB[i]=0 ;
      ele_dszB[i]=0 ;
		
      ele_dxyPV[i]=0 ;
      ele_dzPV[i]=0 ;
      ele_dszPV[i]=0 ;
      ele_dxyPV_error[i]=0 ;
      ele_dzPV_error[i]=0 ;
      ele_dszPV_error[i]=0 ;
		
      ele_isConversion[i] = 0;
      ele_convFound[i] = 0;
      ele_conv_dist[i] = 0.;
      ele_conv_dcot[i] = 0.;
		
		
      ele_track_x[i] = 0.;
      ele_track_y[i] = 0.;
      ele_track_z[i] = 0.;
		
      ele_lost_hits[i]=0;
      ele_chi2_hits[i]=0;
      ele_vertex_x[i]=0; 
      ele_vertex_y[i]=0;
      ele_vertex_z[i]=0;
      ele_tkSumPt_dr03[i]=0;
      ele_ecalRecHitSumEt_dr03[i]=0;
      ele_hcalDepth1TowerSumEt_dr03[i]=0;
      ele_hcalDepth2TowerSumEt_dr03[i]=0;
      ele_hcalDepth1plus2TowerSumEt_00615dr03[i]=0;
      ele_hcalDepth1plus2TowerSumEt_005dr03[i]=0;
      ele_hcalDepth1plus2TowerSumEt_0dr03[i]=0;

      ele_tkSumPt_dr04[i]=0;
      ele_ecalRecHitSumEt_dr04[i]=0;
      ele_hcalDepth1TowerSumEt_dr04[i]=0;
      ele_hcalDepth2TowerSumEt_dr04[i]=0;
      ele_tkSumPtTdrHzz_dr025[i]=0;
      ele_tkSumPtoPtTdrHzz_dr025[i]=0;
      ele_hcalSumEtTdrHzz_dr02[i]=0;
      ele_hcalSumEtoPtTdrHzz_dr02[i]=0;
      ele_tkSumPtEg4Hzz_dr03[i]=0;		
      ele_ecalSumEtEg4Hzz_dr03[i]=0;
      ele_hcalSumEtEg4Hzz_dr04[i]=0;
      ele_tkSumPtoPtEg4Hzz_dr03[i]=0;
      ele_ecalSumEtoPtEg4Hzz_dr03[i]=0;		
      ele_hcalSumEtoPtEg4Hzz_dr04[i]=0;
      ele_ambiguousGsfTracks[i]=0;  
		
      ele_ECAL_fbrem[i]=0; 
      ele_PFcomb[i]=0; 
      ele_PFcomb_Err[i]=0; 
      ele_PF_SCenergy[i]=0; 
      ele_PF_SCenergy_Err[i]=0; 


      for (int j=0;j<5;j++) 
	{
	  ele_ambiguousGsfTracksdxy[i][j]=0 ;
	  ele_ambiguousGsfTracksdz[i][j]=0 ;
	  ele_ambiguousGsfTracksdxyB[i][j]=0 ;
	  ele_ambiguousGsfTracksdzB[i][j]=0 ;
	}
      ele_seedSubdet2[i] = -1;
      ele_seedDphi2Pos[i]   = -20.;
      ele_seedDrz2Pos[i]    = -20.;
      ele_seedDphi2Neg[i]   = -20.;
      ele_seedDrz2Neg[i]    = -20.;
		
      ele_seedSubdet1[i] = -1;
      ele_seedDphi1Pos[i]   = -20.;
      ele_seedDrz1Pos[i]    = -20.;
      ele_seedDphi1Neg[i]   = -20.;
      ele_seedDrz1Neg[i]    = -20.;
		
		
      ele_isMCEle[i] = 0;
      ele_isMCPhoton[i] = 0;
      ele_isMCHadron[i] = 0;
      ele_isSIM[i] = 0;
      ele_isSIMEle[i] = 0;
      ele_idPDGMatch[i] = 0;
      ele_idPDGmother_MCEle[i] = 0;
      ele_idPDGMatchSim[i] = 0;
		
      // Flags for Spike, etc...
      ele_severityLevelSeed[i]     = 0.;
      ele_severityLevelClusters[i] = 0.;
      ele_outOfTimeSeed[i]         = 0.;
      ele_outOfTimeClusters[i]     = 0.;
		
      ele_expected_inner_hits[i]=-1;

		
      ele_sclNclus[i]=-1;
		
      ele_chargeGsfSC[i]=-1;
      ele_chargeGsfCtf[i]=-1;
      ele_chargeGsfCtfSC[i]=-1;
      ele_CtfTrackExists[i]=-1;
      ele_chargeDPhiInnEle[i]=-999.;
      ele_chargeDPhiInnEleCorr[i]=-999.;
      ele_chargeQoverPGsfVtx[i]=-999.;
      ele_chargeQoverPCtf[i]=-999.; 
		
		
		
		
    } // for loop on Ele
	
  for(int i=0; i<100; ++i){
    ele_SeedIsEcalDriven[i] = 0;
    ele_SeedIsTrackerDriven[i] = 0;
		
    ele_SeedSubdet1[i] = -1;
    ele_SeedDphi1Pos[i]   = -20.;
    ele_SeedDrz1Pos[i]    = -20.;
    ele_SeedDphi1Neg[i]   = -20.;
    ele_SeedDrz1Neg[i]    = -20.;
		
    ele_SeedSubdet2[i] = -1;
    ele_SeedDphi2Pos[i]   = -20.;
    ele_SeedDrz2Pos[i]    = -20.;
    ele_SeedDphi2Neg[i]   = -20.;
    ele_SeedDrz2Neg[i]    = -20.;
  } // ele Seed
	
	
  // MET
  _met_calo_et  = 0.;
  _met_calo_px  = 0.; 
  _met_calo_py  = 0.; 
  _met_calo_phi = 0.; 
  _met_calo_set = 0.; 
  _met_calo_sig = 0.; 
	
  _met_calomu_et  = 0.;
  _met_calomu_px  = 0.; 
  _met_calomu_py  = 0.;
  _met_calomu_phi = 0.; 
  _met_calomu_set = 0.; 
  _met_calomu_sig = 0; 
	
  _met_tc_et  = 0.;
  _met_tc_px  = 0.; 
  _met_tc_py  = 0.; 
  _met_tc_phi = 0.; 
  _met_tc_set = 0.; 
  _met_tc_sig = 0.; 
	
  _met_pf_et  = 0.;
  _met_pf_px  = 0.; 
  _met_pf_py  = 0.; 
  _met_pf_phi = 0.; 
  _met_pf_set = 0.; 
  _met_pf_sig = 0.; 
	
  // Muons
  _muons_N = 0; 
	
  for(int im=0;im<20;im++) {
    _muons_charge[im] = 0;
    // Provenance
    _muons_istracker[im]    = 0;
    _muons_isstandalone[im] = 0;
    _muons_isglobal[im]     = 0;
    // Quality cuts
    _muons_dxy[im]            = 0.;
    _muons_dz[im]             = 0.;
    _muons_dxyPV[im]            = 0.;
    _muons_dzPV[im]             = 0.;
    _muons_normalizedChi2[im] = 0.;
    _muons_NtrackerHits[im]   = 0; 
    _muons_NpixelHits[im]     = 0; 
    _muons_NmuonHits[im]      = 0; 
    _muons_Nmatches[im]       = 0; 
    // Isolation
    _muons_nTkIsoR03[im] = 0; 
    _muons_nTkIsoR05[im] = 0; 
    _muons_tkIsoR03[im]  = 0.;
    _muons_tkIsoR05[im]  = 0.;
    _muons_emIsoR03[im]  = 0.;
    _muons_emIsoR05[im]  = 0.;
    _muons_hadIsoR03[im] = 0.;
    _muons_hadIsoR05[im] = 0.;
		
    _muons_trkDxy[im] = 0.;
    _muons_trkDxyError[im] = 0.;
    _muons_trkDxyB[im] = 0.;
    _muons_trkDz[im] = 0.;
    _muons_trkDzError[im] = 0.;
    _muons_trkDzB[im] = 0.; 
    _muons_trkChi2PerNdof[im] = 0.;
    _muons_trkCharge[im] = 0.;
    _muons_trkNHits[im] = 0.;
    _muons_trkNPixHits[im] = 0.;
    _muons_trkmuArbitration[im] = 0.;
    _muons_trkmu2DCompatibilityLoose[im] = 0.;
    _muons_trkmu2DCompatibilityTight[im] = 0.;
    _muons_trkmuOneStationLoose[im] = 0.;
    _muons_trkmuOneStationTight[im] = 0.;
    _muons_trkmuLastStationLoose[im] = 0.;
    _muons_trkmuLastStationTight[im] = 0.;
    _muons_trkmuOneStationAngLoose[im] = 0.;
    _muons_trkmuOneStationAngTight[im] = 0.;
    _muons_trkmuLastStationAngLoose[im] = 0.;
    _muons_trkmuLastStationAngTight[im] = 0.;
    _muons_trkmuLastStationOptimizedLowPtLoose[im] = 0.;
    _muons_trkmuLastStationOptimizedLowPtTight[im] = 0.;
    _muons_caloCompatibility[im] = 0.;
    _muons_segmentCompatibility[im] = 0.;
    _muons_glbmuPromptTight[im] = 0.;
		
    _muons_hzzIso[im] = 0.;
    _muons_hzzIsoTk[im] = 0.;
    _muons_hzzIsoEcal[im] = 0.; 
    _muons_hzzIsoHcal[im] = 0.; 
		
    muons_Tip[im] = -999. ;
    muons_Lip[im] = -999. ;
    muons_STip[im] = -999. ;
    muons_SLip[im] = -999. ;
    muons_TipSignif[im] = -999. ;
    muons_LipSignif[im] = -999. ;
    muons_Significance3D[im] = -999. ;
    muons_Value3D[im] = -999. ;
    muons_Error3D[im] = -999. ;
  } // for loop on muons
	
  // Calo Jets
  _jets_calo_N = 0;
	
  // JPT Jets
  _jets_jpt_N = 0;
	
  // PF Jets
  _jets_pf_N = 0;
	
  for(int ipfjet=0;ipfjet<100;ipfjet++) {
    jets_pf_chargedHadEFrac[ipfjet] = 0.;
    jets_pf_chargedEmEFrac[ipfjet]  = 0.;
    jets_pf_chargedMuEFrac[ipfjet]  = 0.;
		
    jets_pf_neutralHadEFrac[ipfjet] = 0.;
    jets_pf_neutralEmEFrac[ipfjet]  = 0.;
    jets_pf_PhotonEFrac[ipfjet]     = 0.;
		
    jets_pf_chargedHadMultiplicity[ipfjet] = 0;
    jets_pf_neutralHadMultiplicity[ipfjet] = 0;
		
    jets_pf_chargedMultiplicity[ipfjet] = 0;
    jets_pf_neutralMultiplicity[ipfjet] = 0;
		
    jets_pf_nConstituents[ipfjet]      = 0;
  } // for loop on PFjets
	
  // SuperClusters
  /*_sc_hybrid_N = 0; 
  for(int isc=0;isc<25;isc++) {
    _sc_hybrid_E[isc]   = 0.; 
    _sc_hybrid_Et[isc]  = 0.; 
    _sc_hybrid_Eta[isc] = 0.; 
    _sc_hybrid_Phi[isc] = 0.; 
    _sc_hybrid_outOfTimeSeed[isc]     = 0;
    _sc_hybrid_severityLevelSeed[isc] = 0;
    _sc_hybrid_e1[isc]  = 0.;
    _sc_hybrid_e33[isc] = 0.;
    _sc_hybrid_he[isc]  = -10.;
    _sc_hybrid_sigmaietaieta[isc] = 0.;
    _sc_hybrid_hcalDepth1TowerSumEt_dr03[isc] = 0.;
    _sc_hybrid_hcalDepth2TowerSumEt_dr03[isc] = 0.;
    _sc_hybrid_ecalRecHitSumEt_dr03[isc]      = 0.;
    _sc_hybrid_trkiso_dr03[isc]               = 0.;
		
    _sc_hybrid_RCTeta[isc]=-999;
    _sc_hybrid_RCTphi[isc]=-999;
    _sc_hybrid_RCTL1iso[isc]     = -999;
    _sc_hybrid_RCTL1noniso[isc]  = -999;
		
    for (int li=0;li<50;li++) {
      _sc_hybrid_TTetaVect[isc][li]=-999;
      _sc_hybrid_TTphiVect[isc][li]=-999;
      _sc_hybrid_TTetVect[isc][li]=0.;
    } // for loop on 50
    for (int li=0;li<10;li++) {
      _sc_hybrid_RCTetaVect[isc][li]=-999;
      _sc_hybrid_RCTphiVect[isc][li]=-999;
      _sc_hybrid_RCTetVect[isc][li]=0.;
      _sc_hybrid_RCTL1isoVect[isc][li]=-999;
      _sc_hybrid_RCTL1nonisoVect[isc][li]=-999;
    } // for loop on 10
		
  } // for loop on EB superclusters
	
  _sc_multi55_N = 0; 
  for(int isc=0;isc<25;isc++) {
    _sc_multi55_E[isc]   = 0.; 
    _sc_multi55_Et[isc]  = 0.; 
    _sc_multi55_Eta[isc] = 0.; 
    _sc_multi55_Phi[isc] = 0.; 
    _sc_multi55_he[isc]  = -10.;
    _sc_multi55_sigmaietaieta[isc] = 0.;
    _sc_multi55_hcalDepth1TowerSumEt_dr03[isc] = 0.;
    _sc_multi55_hcalDepth2TowerSumEt_dr03[isc] = 0.;
    _sc_multi55_ecalRecHitSumEt_dr03[isc]      = 0.;
    _sc_multi55_trkiso_dr03[isc]               = 0.;
		
    _sc_multi55_RCTeta[isc]=-999;
    _sc_multi55_RCTphi[isc]=-999;
    _sc_multi55_RCTL1iso[isc]     = -999;
    _sc_multi55_RCTL1noniso[isc]  = -999;
		
    for (int li=0;li<50;li++) {
      _sc_multi55_TTetaVect[isc][li]=-999;
      _sc_multi55_TTphiVect[isc][li]=-999;
      _sc_multi55_TTetVect[isc][li]=0.;
    } // for loop on 50
    for (int li=0;li<10;li++) {
      _sc_multi55_RCTetaVect[isc][li]=-999;
      _sc_multi55_RCTphiVect[isc][li]=-999;
      _sc_multi55_RCTetVect[isc][li]=0.;
      _sc_multi55_RCTL1isoVect[isc][li]=-999;
      _sc_multi55_RCTL1nonisoVect[isc][li]=-999;
    } // for loop on 10
		
    } // for loop on EE superclusters*/
	
  // Generated W, Z & leptons
  for(int igen=0;igen<10;igen++) {
    _MC_gen_V_pdgid[igen]       = 0.;
    _MC_gen_leptons_pdgid[igen] = 0.;
  } // for loop on igen
	
  for(int i=0;i<10;++i) {
	
    ele_Tip[i] = -999. ;
    ele_Lip[i] = -999. ;
    ele_STip[i] = -999. ;
    ele_SLip[i] = -999. ;
    ele_TipSignif[i] = -999. ;
    ele_LipSignif[i] = -999. ;
    ele_Significance3D[i] = -999. ;
    ele_Value3D[i] = -999. ;
    ele_Error3D[i] = -999. ;
  }
	
}

// ====================================================================================
void SimpleNtpleCustom::beginJob(const edm::ParameterSet& conf)
// ====================================================================================
{

	
}

// ====================================================================================
void SimpleNtpleCustom::endJob() {}
// ====================================================================================

// ====================================================================================
void SimpleNtpleCustom::setMomentum (TLorentzVector &myvector, const LorentzVector & mom)
// ====================================================================================
{
	myvector.SetPx (mom.Px());
	myvector.SetPy (mom.Py());
	myvector.SetPz (mom.Pz());
	myvector.SetE (mom.E());
}

// ====================================================================================
bool SimpleNtpleCustom::IsConv (const reco::GsfElectron & eleRef)
// ====================================================================================
{
	
	bool isAmbiguous = true, isNotFromPixel = true;
	if (eleRef.ambiguousGsfTracksSize() == 0 ) {isAmbiguous = false;}
	
	int  mishits = eleRef.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
	

	if (mishits == 0){isNotFromPixel = false;}
	
	bool is_conversion = false;
	
	if(isAmbiguous || isNotFromPixel) is_conversion = true;
	
	return is_conversion;
	
}

// ====================================================================================
const EcalRecHit SimpleNtpleCustom::getRecHit(DetId id, const EcalRecHitCollection *recHits)
// ====================================================================================
{
	if ( id == DetId(0) ) {
		return EcalRecHit();
	} else {
		EcalRecHitCollection::const_iterator it = recHits->find( id );
		if ( it != recHits->end() ) {
			return (*it);
		} else {
		  return EcalRecHit();
		}
	}
	return EcalRecHit();
}



//modif-alex
//GETTING RCT regions
// ====================================================================================
int SimpleNtpleCustom::getGCTRegionPhi(int ttphi)
// ====================================================================================
{
	int gctphi=0;
	gctphi = (ttphi+1)/4;
	if(ttphi<=2) gctphi=0;
	if(ttphi>=71) gctphi=0;
	
	return gctphi;
}

// ====================================================================================
int SimpleNtpleCustom::getGCTRegionEta(int tteta)
// ====================================================================================
{
	int gcteta = 0;
	
	if(tteta>0) gcteta = (tteta-1)/4 + 11;
	else if(tteta<0) gcteta = (tteta+1)/4 + 10;
	
	return gcteta;
}


// ====================================================================================
int SimpleNtpleCustom::getCaloEta_fromGTEta(int GTEta)
// ====================================================================================
{
	int CaloEta = 0;
	int absGTEta = abs(GTEta);
	
	if(absGTEta<=39) CaloEta = absGTEta/2 + 1;
	else if(absGTEta==40) CaloEta = 20;
	else if(absGTEta>=41 && absGTEta<=42) CaloEta = 21;
	else if(absGTEta>=43 && absGTEta<=44) CaloEta = 22;
	else if(absGTEta>=45 && absGTEta<=47) CaloEta = 23;
	else if(absGTEta>=48 && absGTEta<=50) CaloEta = 24;
	else if(absGTEta>=51 && absGTEta<=53) CaloEta = 25;
	else if(absGTEta>=54 && absGTEta<=57) CaloEta = 26;
	else if(absGTEta>=58 && absGTEta<=61) CaloEta = 27;
	else if(absGTEta>=62) CaloEta = 28;
	
	if(GTEta<0) CaloEta*=-1;

	return CaloEta;
}




// ====================================================================================
int SimpleNtpleCustom::getCaloPhi_fromGTPhi(int GTPhi)
// ====================================================================================
{
	int CaloPhi = GTPhi/2 + 1;
	return CaloPhi;
}



// ====================================================================================
float SimpleNtpleCustom::E2overE9( const DetId id, const EcalRecHitCollection & recHits, 
			     float recHitEtThreshold, float recHitEtThreshold2 , 
			     bool avoidIeta85, bool KillSecondHit)
// ====================================================================================
// taken from CMSSW/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc CMSSW_3_9_0_pre5

{

        // compute e2overe9
        //  
        //   | | | |
        //   +-+-+-+
        //   | |1|2|
        //   +-+-+-+
        //   | | | |
        //
        //   1 - input hit,  2 - highest energy hit in a 3x3 around 1
        // 
        //   rechit 1 must have E_t > recHitEtThreshold
        //   rechit 2 must have E_t > recHitEtThreshold2
        //
        //   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2
        //
        //   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2
        //   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0


        if ( id.subdetId() == EcalBarrel ) {
	  
                EBDetId ebId( id );

                // avoid recHits at |eta|=85 where one side of the neighbours is missing
                if ( abs(ebId.ieta())==85 && avoidIeta85){  return 0;}

                // select recHits with Et above recHitEtThreshold

 
                float e1 = recHitE( id, recHits );
		

                float ete1=recHitApproxEt( id, recHits );


		// check that rechit E_t is above threshold

		if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) { return 0;}
		
		if (ete1 < recHitEtThreshold && !KillSecondHit ) {return 0;}
		

                float e2=-1;
                float ete2=0;
                float s9 = 0;

                // coordinates of 2nd hit relative to central hit
                int e2eta=0;
                int e2phi=0;

		// LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1

                for ( int deta = -1; deta <= +1; ++deta ) {
                   for ( int dphi = -1; dphi <= +1; ++dphi ) {
 
		      // compute 3x3 energy

                      float etmp=recHitE( id, recHits, deta, dphi );
                      s9 += etmp;

                      EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
                      float eapproxet=recHitApproxEt( idtmp, recHits );

                      // remember 2nd highest energy deposit (above threshold) in 3x3 array 
                      if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {

                         e2=etmp;
                         ete2=eapproxet;
                         e2eta=deta;
                         e2phi=dphi;
        
                      }

                   }
                }

                if ( e1 == 0 )  { return 0;}
  
                // return 0 if 2nd hit is below threshold
                if ( e2 == -1 ) {return 0;}

                // compute e2/e9 centered around 1st hit

                float e2nd=e1+e2;
                float e2e9=0;

                if (s9!=0) e2e9=e2nd/s9;
  
                // if central hit has higher energy than 2nd hit
                //  return e2/e9 if 1st hit is above E_t threshold

                if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;

                // if second hit has higher energy than 1st hit

                if ( e2 > e1 ) { 


                  // return 0 if user does not want to flag 2nd hit, or
                  // hits are below E_t thresholds - note here we
		  // now assume the 2nd hit to be the leading hit.

		  if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
		    
                     return 0;
  
                 }


                  else {
 
                    // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2

		    float s92nd=0;
           
                    float e2nd_prime=0;
                    int e2prime_eta=0;
                    int e2prime_phi=0;

                    EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);


                     for ( int deta = -1; deta <= +1; ++deta ) {
                        for ( int dphi = -1; dphi <= +1; ++dphi ) {
 
		           // compute 3x3 energy

                           float etmp=recHitE( secondid, recHits, deta, dphi );
                           s92nd += etmp;

                           if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
			     e2nd_prime=etmp;
                             e2prime_eta=deta;
                             e2prime_phi=dphi;
			   }

			}
		     }

		     // if highest energy hit around E2 is not the same as the input hit, return 0;

		     if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi)) 
		       { 
			 return 0;
		       }


		     // compute E2/E9 around second hit 
		     float e2e9_2=0;
		     if (s92nd!=0) e2e9_2=e2nd/s92nd;
                 
		     //   return the value of E2/E9 calculated around 2nd hit
                   
		     return e2e9_2;


		  }
		  
		}


        } else if ( id.subdetId() == EcalEndcap ) {
	  // only used for EB at the moment
          return 0;
        }
        return 0;
}



// ====================================================================================
float SimpleNtpleCustom::recHitE( const DetId id, const EcalRecHitCollection &recHits )
// ====================================================================================
{
        if ( id == DetId(0) ) {
                return 0;
        } else {
                EcalRecHitCollection::const_iterator it = recHits.find( id );
                if ( it != recHits.end() ) return (*it).energy();
        }
        return 0;
}


// ====================================================================================
float SimpleNtpleCustom::recHitE( const DetId id, const EcalRecHitCollection & recHits,
                                           int di, int dj )
// ====================================================================================
{
        // in the barrel:   di = dEta   dj = dPhi
        // in the endcap:   di = dX     dj = dY
  
        DetId nid;
        if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
        else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );

        return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
}




// ====================================================================================
float SimpleNtpleCustom::recHitApproxEt( const DetId id, const EcalRecHitCollection &recHits )
// ====================================================================================
{
        // for the time being works only for the barrel
        if ( id.subdetId() == EcalBarrel ) {
                return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
        }
        return 0;
}
DEFINE_FWK_MODULE(SimpleNtpleCustom);
