#include "EGamma/ECGelec/plugins/SimpleNtpleAOD.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace IPTools;
//using namespace math;

// ====================================================================================
SimpleNtpleAOD::SimpleNtpleAOD(const edm::ParameterSet& iConfig) :

  // Nadir study
  EcalRecHitCollectionEB_ (iConfig.getParameter<edm::InputTag>("EcalRecHitCollectionEB")),
  tpCollectionNormal_ (iConfig.getParameter<edm::InputTag> ("TPCollectionNormal") ),
  tpCollectionModif_ (iConfig.getParameter<edm::InputTag> ("TPCollectionModif") ),
  tpEmulatorCollection_ (iConfig.getParameter<edm::InputTag> ("TPEmulatorCollection") ),
  
  nadGetL1M_ (iConfig.getUntrackedParameter<bool>("NadL1M")),
  nadGetTP_ (iConfig.getUntrackedParameter<bool>("NadTP")),
  nadGetTP_Modif_ (iConfig.getUntrackedParameter<bool>("NadTPmodif")),
  nadGetTP_Emul_ (iConfig.getUntrackedParameter<bool>("NadTPemul")),
  PrintDebug_ (iConfig.getUntrackedParameter<bool>("PrintDebug")),
  PrintDebug_HLT_ (iConfig.getUntrackedParameter<bool>("PrintDebug_HLT")),  

  hcalTowers_ (iConfig.getParameter<edm::InputTag>("hcalTowers")),
  hOverEConeSize_ (iConfig.getParameter<double>("hOverEConeSize")),
  hOverEPtMin_ (iConfig.getParameter<double>("hOverEPtMin")),
  EleID_VeryLooseTag_ (iConfig.getParameter<edm::InputTag> ("eleID_VeryLooseTag")) ,
  EleID_LooseTag_ (iConfig.getParameter<edm::InputTag> ("eleID_LooseTag")) ,
  EleID_MediumTag_ (iConfig.getParameter<edm::InputTag> ("eleID_MediumTag")) ,
  EleID_TightTag_ (iConfig.getParameter<edm::InputTag> ("eleID_TightTag")) , 
  EleID_SuperTightTag_ (iConfig.getParameter<edm::InputTag> ("eleID_SuperTightTag")) , 
  EleID_HyperTight1Tag_ (iConfig.getParameter<edm::InputTag> ("eleID_HyperTight1Tag")) ,
  EleID_HyperTight2Tag_ (iConfig.getParameter<edm::InputTag> ("eleID_HyperTight2Tag")) ,
  EleID_HyperTight3Tag_ (iConfig.getParameter<edm::InputTag> ("eleID_HyperTight3Tag")) ,
  EleID_HyperTight4Tag_ (iConfig.getParameter<edm::InputTag> ("eleID_HyperTight4Tag")) ,
  EleIso_TdrHzzTkMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_TdrHzzTkMapTag")),
  EleIso_TdrHzzHcalMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_TdrHzzHcalMapTag")),
  EleIso_Eg4HzzTkMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_Eg4HzzTkMapTag")),
  EleIso_Eg4HzzEcalMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_Eg4HzzEcalMapTag")),
  EleIso_Eg4HzzHcalMapTag_ (iConfig.getParameter<edm::InputTag>("eleIso_Eg4HzzHcalMapTag")),
  EleTag_ (iConfig.getParameter<edm::InputTag> ("EleTag")),
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
  // Trigger Stuff
  HLTTag_(iConfig.getParameter<edm::InputTag> ("HLTTag")),
  triggerEventTag_(iConfig.getParameter<edm::InputTag> ("TriggerEventTag")),
  //

  type_ (iConfig.getParameter<std::string>("type")),
  aod_ (iConfig.getUntrackedParameter<bool>("AOD")),
  funcname_  (iConfig.getParameter<std::string>("functionName")),
  useBeamSpot_ (iConfig.getParameter<bool>("useBeamSpot")),
  PileupSrc_ ("addPileupInfo"),
  RhoCorrection_("kt6PFJets:rho"),
  SigmaRhoCorrection_("kt6PFJets:sigma")
				    //PileupSrc_(iConfig.getParameter<edm::InputTag>("Pileupsrc"))
				    //,"addPileupInfo::REDIGI311X"))
				    //
				    // ====================================================================================
{
  //now do what ever initialization is needed
  funcbase_ = EcalClusterFunctionFactory::get()->create( funcname_, iConfig ); 
  gtRecordCollectionTag_ = iConfig.getParameter<std::string>("GTRecordCollection") ;
	
  HLT_ElePaths_  = iConfig.getParameter<std::vector<std::string > >("HLTElePaths");
  HLT_MuonPaths_ = iConfig.getParameter<std::vector<std::string > >("HLTMuonPaths");
  HLT_Filters_   = iConfig.getParameter<std::vector<edm::InputTag > >("HLTFilters");
	
  simulation_ = iConfig.getUntrackedParameter<bool>("simulation", false);
  fillsc_     = iConfig.getUntrackedParameter<bool>("FillSC", false);
	
  edm::Service<TFileService> fs ;
  mytree_  = fs->make <TTree>("eIDSimpleTree","eIDSimpleTree"); 
	
  // Global
  mytree_->Branch("nEvent",&nEvent,"nEvent/I");
  mytree_->Branch("nRun",&nRun,"nRun/I");
  mytree_->Branch("nLumi",&nLumi,"nLumi/I");
	
  // Pile UP
  mytree_->Branch("PU_N",&_PU_N,"PU_N/I");
  mytree_->Branch("PU_rhoCorr",&_PU_rho,"PU_rhoCorr/D");
  mytree_->Branch("PU_sigmaCorr",&_PU_sigma,"PU_sigmaCorr/D");

  // Vertices
  mytree_->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  mytree_->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[200]/D");
  mytree_->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[200]/D");
  mytree_->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[200]/D");
  mytree_->Branch("vtx_d0",&_vtx_d0,"vtx_d0[200]/D");
  mytree_->Branch("vtx_x",&_vtx_x,"vtx_x[200]/D");
  mytree_->Branch("vtx_y",&_vtx_y,"vtx_y[200]/D");
  mytree_->Branch("vtx_z",&_vtx_z,"vtx_z[200]/D");
	
  // Skim
  mytree_->Branch("skim_is1lepton", &_skim_is1lepton, "skim_is1lepton/I");
  mytree_->Branch("skim_is2leptons",&_skim_is2leptons,"skim_is2leptons/I");
  mytree_->Branch("skim_is3leptons",&_skim_is3leptons,"skim_is3leptons/I");

  // Spikes
  mytree_->Branch("spike_N",&spike_N,"spike_N/I");
  mytree_->Branch("spike_TTieta",&spike_TTieta,"spike_TTieta[5000]/I");
  mytree_->Branch("spike_TTiphi",&spike_TTiphi,"spike_TTiphi[5000]/I");
  mytree_->Branch("spike_Rieta",&spike_Rieta,"spike_Rieta[5000]/I");
  mytree_->Branch("spike_Riphi",&spike_Riphi,"spike_Riphi[5000]/I");
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
  mytree_->Branch("trig_HLT_path",&trig_HLT_path,"trig_HLT_path[4]/I");
  // unbias, EG5, EG8, EG12
  //
  mytree_->Branch("trig_fired_names",&trig_fired_names,"trig_fired_names[5000]/C");
  mytree_->Branch("trig_hltInfo",&trig_hltInfo,"trig_hltInfo[250]/I");
  //  mytree_->Branch("trig_isUnbiased",&trig_isUnbiased,"trig_isUnbiased/I");
  mytree_->Branch("trig_isPhoton10",&trig_isPhoton10,"trig_isPhoton10/I"); 
  mytree_->Branch("trig_isPhoton15",&trig_isPhoton15,"trig_isPhoton15/I"); 
  //   mytree_->Branch("trig_isL1SingleEG5",&trig_isL1SingleEG5,"trig_isL1SingleEG5/I");
  //   mytree_->Branch("trig_isL1SingleEG8",&trig_isL1SingleEG8,"trig_isL1SingleEG8/I");
  //   mytree_->Branch("trig_isL1SingleEG12",&trig_isL1SingleEG2,"trig_isL1SingleEG12/I"); 
  mytree_->Branch("trig_isEle10_LW",&trig_isEle10_LW,"trig_isEle10_LW/I");
  mytree_->Branch("trig_isEle15_LW",&trig_isEle15_LW,"trig_isEle15_LW/I");
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
	
  mytree_->Branch("ele_eidVeryLoose",ele_eidVeryLoose,"ele_eidVeryLoose[10]/D");
  mytree_->Branch("ele_eidLoose",ele_eidLoose,"ele_eidLoose[10]/D");
  mytree_->Branch("ele_eidMedium",ele_eidMedium,"ele_eidMedium[10]/D");
  mytree_->Branch("ele_eidTight",ele_eidTight,"ele_eidTight[10]/D"); 
  mytree_->Branch("ele_eidSuperTight",ele_eidSuperTight,"ele_eidSuperTight[10]/D"); 
  mytree_->Branch("ele_eidHyperTight1",ele_eidHyperTight1,"ele_eidHyperTight1[10]/D");
  mytree_->Branch("ele_eidHyperTight2",ele_eidHyperTight2,"ele_eidHyperTight2[10]/D");
  mytree_->Branch("ele_eidHyperTight3",ele_eidHyperTight3,"ele_eidHyperTight3[10]/D");
  mytree_->Branch("ele_eidHyperTight4",ele_eidHyperTight4,"ele_eidHyperTight4[10]/D");
	
	
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
  //mytree_->Branch("tkIso03Rel",tkIso03Rel,"tkIso03Rel[10]/D");
  //mytree_->Branch("ecalIso03Rel",ecalIso03Rel,"ecalIso03Rel[10]/D");
  //mytree_->Branch("hcalIso03Rel",hcalIso03Rel,"hcalIso03Rel[10]/D");
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
  mytree_->Branch("ele_TTetaVect",      &_ele_TTetaVect,       "ele_TTetaVect[10][50]/I");
  mytree_->Branch("ele_TTphiVect",      &_ele_TTphiVect,       "ele_TTphiVect[10][50]/I");
  mytree_->Branch("ele_TTetVect",       &_ele_TTetVect,        "ele_TTetVect[10][50]/D");
  mytree_->Branch("ele_RCTetaVect",     &_ele_RCTetaVect,      "ele_RCTetaVect[10][10]/I");
  mytree_->Branch("ele_RCTphiVect",     &_ele_RCTphiVect,      "ele_RCTphiVect[10][10]/I");
  mytree_->Branch("ele_RCTetVect",      &_ele_RCTetVect,       "ele_RCTetVect[10][10]/D");
  mytree_->Branch("ele_RCTL1isoVect",   &_ele_RCTL1isoVect,    "ele_RCTL1isoVect[10][10]/I");
  mytree_->Branch("ele_RCTL1nonisoVect",&_ele_RCTL1nonisoVect, "ele_RCTL1nonisoVect[10][10]/I");
  mytree_->Branch("ele_RCTL1isoVect_M",   &_ele_RCTL1isoVect_M,    "ele_RCTL1isoVect_M[10][10]/I");
  mytree_->Branch("ele_RCTL1nonisoVect_M",&_ele_RCTL1nonisoVect_M, "ele_RCTL1nonisoVect_M[10][10]/I");
	
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
  mytree_->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);
	
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
	
  // SuperClusters
  // SC EB
  mytree_->Branch("sc_hybrid_N",   &_sc_hybrid_N,   "sc_hybrid_N/I");
  mytree_->Branch("sc_hybrid_E",   &_sc_hybrid_E,   "sc_hybrid_E[25]/D");
  mytree_->Branch("sc_hybrid_Et",  &_sc_hybrid_Et,  "sc_hybrid_Et[25]/D");
  mytree_->Branch("sc_hybrid_Eta", &_sc_hybrid_Eta, "sc_hybrid_Eta[25]/D");
  mytree_->Branch("sc_hybrid_Phi", &_sc_hybrid_Phi, "sc_hybrid_Phi[25]/D");
  mytree_->Branch("sc_hybrid_outOfTimeSeed",     &_sc_hybrid_outOfTimeSeed,     "sc_hybrid_outOfTimeSeed[25]/I");
  mytree_->Branch("sc_hybrid_severityLevelSeed", &_sc_hybrid_severityLevelSeed, "sc_hybrid_severityLevelSeed[25]/I");
  mytree_->Branch("sc_hybrid_e1", &_sc_hybrid_e1, "sc_hybrid_e1[25]/D");
  mytree_->Branch("sc_hybrid_e33", &_sc_hybrid_e33, "sc_hybrid_e33[25]/D");
  mytree_->Branch("sc_hybrid_he",&_sc_hybrid_he, "sc_hybrid_he[25]/D");
  mytree_->Branch("sc_hybrid_sigmaietaieta",&_sc_hybrid_sigmaietaieta, "sc_hybrid_sigmaietaieta[25]/D");
  mytree_->Branch("sc_hybrid_hcalDepth1TowerSumEt_dr03", &_sc_hybrid_hcalDepth1TowerSumEt_dr03, "sc_hybrid_hcalDepth1TowerSumEt_dr03[25]/D");
  mytree_->Branch("sc_hybrid_hcalDepth2TowerSumEt_dr03", &_sc_hybrid_hcalDepth2TowerSumEt_dr03, "sc_hybrid_hcalDepth2TowerSumEt_dr03[25]/D");
  mytree_->Branch("sc_hybrid_ecalRecHitSumEt_dr03",      &_sc_hybrid_ecalRecHitSumEt_dr03,      "sc_hybrid_ecalRecHitSumEt_dr03[25]/D");
  mytree_->Branch("sc_hybrid_trkiso_dr03",               &_sc_hybrid_trkiso_dr03,               "sc_hybrid_trkiso_dr03[25]/D");
  // SC EB : L1 stuff
  mytree_->Branch("sc_hybrid_TTetaVect", &_sc_hybrid_TTetaVect, "_sc_hybrid_TTetaVect[25][50]/I");
  mytree_->Branch("sc_hybrid_TTphiVect", &_sc_hybrid_TTphiVect, "_sc_hybrid_TTphiVect[25][50]/I");
  mytree_->Branch("sc_hybrid_TTetVect", &_sc_hybrid_TTetVect, "_sc_hybrid_TTetVect[25][50]/D");
  mytree_->Branch("sc_hybrid_RCTetaVect", &_sc_hybrid_RCTetaVect, "_sc_hybrid_RCTetaVect[25][10]/I");
  mytree_->Branch("sc_hybrid_RCTphiVect", &_sc_hybrid_RCTphiVect, "_sc_hybrid_RCTphiVect[25][10]/I");
  mytree_->Branch("sc_hybrid_RCTetVect", &_sc_hybrid_RCTetVect, "_sc_hybrid_RCTetVect[25][10]/D");
  mytree_->Branch("sc_hybrid_RCTL1isoVect", &_sc_hybrid_RCTL1isoVect, "_sc_hybrid_RCTL1isoVect[25][10]/I");
  mytree_->Branch("sc_hybrid_RCTL1nonisoVect", &_sc_hybrid_RCTL1nonisoVect, "_sc_hybrid_RCTL1nonisoVect[25][10]/I");
  //
  mytree_->Branch("sc_hybrid_RCTeta", &_sc_hybrid_RCTeta, "_sc_hybrid_RCTeta[25]/I");
  mytree_->Branch("sc_hybrid_RCTphi", &_sc_hybrid_RCTphi, "_sc_hybrid_RCTphi[25]/I");
  mytree_->Branch("sc_hybrid_RCTL1iso", &_sc_hybrid_RCTL1iso, "_sc_hybrid_RCTL1iso[25]/I");
  mytree_->Branch("sc_hybrid_RCTL1noniso", &_sc_hybrid_RCTL1noniso, "_sc_hybrid_RCTL1noniso[25]/I");
  // SC EE 
  mytree_->Branch("sc_multi55_N",   &_sc_multi55_N,   "sc_multi55_N/I");
  mytree_->Branch("sc_multi55_E",   &_sc_multi55_E,   "sc_multi55_E[25]/D");
  mytree_->Branch("sc_multi55_Et",  &_sc_multi55_Et,  "sc_multi55_Et[25]/D");
  mytree_->Branch("sc_multi55_Eta", &_sc_multi55_Eta, "sc_multi55_Eta[25]/D");
  mytree_->Branch("sc_multi55_Phi", &_sc_multi55_Phi, "sc_multi55_Phi[25]/D");
  mytree_->Branch("sc_multi55_he",  &_sc_multi55_he,  "sc_multi55_he[25]/D");
  mytree_->Branch("sc_multi55_sigmaietaieta",&_sc_multi55_sigmaietaieta, "sc_multi55_sigmaietaieta[25]/D");
  mytree_->Branch("sc_multi55_hcalDepth1TowerSumEt_dr03", &_sc_multi55_hcalDepth1TowerSumEt_dr03, "sc_multi55_hcalDepth1TowerSumEt_dr03[25]/D");
  mytree_->Branch("sc_multi55_hcalDepth2TowerSumEt_dr03", &_sc_multi55_hcalDepth2TowerSumEt_dr03, "sc_multi55_hcalDepth2TowerSumEt_dr03[25]/D");
  mytree_->Branch("sc_multi55_ecalRecHitSumEt_dr03",      &_sc_multi55_ecalRecHitSumEt_dr03,      "sc_multi55_ecalRecHitSumEt_dr03[25]/D");
  mytree_->Branch("sc_multi55_trkiso_dr03",               &_sc_multi55_trkiso_dr03,               "sc_multi55_trkiso_dr03[25]/D");
  // SC EE : L1 stuff
  mytree_->Branch("sc_multi55_TTetaVect", &_sc_multi55_TTetaVect, "_sc_multi55_TTetaVect[25][50]/I");
  mytree_->Branch("sc_multi55_TTphiVect", &_sc_multi55_TTphiVect, "_sc_multi55_TTphiVect[25][50]/I");
  mytree_->Branch("sc_multi55_TTetVect", &_sc_multi55_TTetVect, "_sc_multi55_TTetVect[25][50]/D");
  mytree_->Branch("sc_multi55_RCTetaVect", &_sc_multi55_RCTetaVect, "_sc_multi55_RCTetaVect[25][10]/I");
  mytree_->Branch("sc_multi55_RCTphiVect", &_sc_multi55_RCTphiVect, "_sc_multi55_RCTphiVect[25][10]/I");
  mytree_->Branch("sc_multi55_RCTetVect", &_sc_multi55_RCTetVect, "_sc_multi55_RCTetVect[25][10]/D");
  mytree_->Branch("sc_multi55_RCTL1isoVect", &_sc_multi55_RCTL1isoVect, "_sc_multi55_RCTL1isoVect[25][10]/I");
  mytree_->Branch("sc_multi55_RCTL1nonisoVect", &_sc_multi55_RCTL1nonisoVect, "_sc_multi55_RCTL1nonisoVect[25][10]/I");
  //
  mytree_->Branch("sc_multi55_RCTeta", &_sc_multi55_RCTeta, "_sc_multi55_RCTeta[25]/I");
  mytree_->Branch("sc_multi55_RCTphi", &_sc_multi55_RCTphi, "_sc_multi55_RCTphi[25]/I");
  mytree_->Branch("sc_multi55_RCTL1iso", &_sc_multi55_RCTL1iso, "_sc_multi55_RCTL1iso[25]/I");
  mytree_->Branch("sc_multi55_RCTL1noniso", &_sc_multi55_RCTL1noniso, "_sc_multi55_RCTL1noniso[25]/I");
	
  // Generated W,Z's & leptons
  _m_MC_gen_V = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("MC_gen_V", "TClonesArray", &_m_MC_gen_V, 256000,0);
  mytree_->Branch ("MC_gen_V_pdgid",&_MC_gen_V_pdgid, "MC_gen_V_pdgid[10]/D");
	
  _m_MC_gen_leptons = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("MC_gen_leptons", "TClonesArray", &_m_MC_gen_leptons, 256000,0);
  mytree_->Branch ("MC_gen_leptons_pdgid",&_MC_gen_leptons_pdgid, "MC_gen_leptons_pdgid[10]/D");
}

// ====================================================================================
SimpleNtpleAOD::~SimpleNtpleAOD()
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
void SimpleNtpleAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // Clemy's Stuff for Charge
  bool updateField(false);
  if (cacheIDMagField_!=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier()){
    updateField = true;
    cacheIDMagField_=iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier();
    iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
  }
	
  bool updateGeometry(false);
  if (cacheIDTDGeom_!=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier()){
    updateGeometry = true;
    cacheIDTDGeom_=iSetup.get<TrackerDigiGeometryRecord>().cacheIdentifier();
    iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
  }
	
  if(updateField || updateGeometry){
    mtsTransform_ = new MultiTrajectoryStateTransform(trackerHandle_.product(),theMagField.product());
  }
	
  //  if (cacheIDTopo_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
  //     cacheIDTopo_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
  //     iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  //   }
	
  // Tree Maker
  //if(PrintDebug_) std::cout << "Init()" << std::endl;
  Init();
  if (funcbase_) funcbase_->init(iSetup);
  //std::cout << "FillEvent (iEvent, iSetup);" << std::endl;
  FillEvent (iEvent, iSetup);
	
  // for Skimming
  AnalysisUtils * skim = new AnalysisUtils();
  _skim_is1lepton  = skim->doSkim(iEvent, iSetup, true, true,   20., 20., 1, 1);
  _skim_is2leptons = skim->doSkim(iEvent, iSetup, false, false, 10., 15., 2, 1);
  _skim_is3leptons = skim->doSkim(iEvent, iSetup, false, false,  5., 10., 3, 2);
  //bool isEleID_, bool isMuonID_, double lep_ptLow_, double lep_ptHigh_, int nLep_ptLow_, int nLep_ptHigh_);
	
  //
  FillTrigger (iEvent, iSetup);
  //std::cout << "m_electrons -> Clear() ;" << std::endl;
  m_electrons -> Clear() ;
  //std::cout << "muons" << std::endl;
  m_muons -> Clear() ;
  //std::cout << "gen V" << std::endl;
  if(type_ == "MC") {
    _m_MC_gen_V->Clear();
    //cout << "gen leptons" << endl;
    _m_MC_gen_leptons->Clear();
  } // if MC
  //std::cout << "FillEle (iEvent, iSetup);" << std::endl;
  FillEle (iEvent, iSetup);
  //
	
  //std::cout << "FillMuon (iEvent, iSetup);" << std::endl;
  //FillMuons (iEvent, iSetup);
  //std::cout << "FillMET (iEvent, iSetup);" << std::endl;
  //FillMET (iEvent, iSetup);
  //std::cout << "FillJets(iEvent, iSetup);" << std::endl;
  //FillJets(iEvent, iSetup);
	
  //std::cout << "if(fillsc_) FillSuperClusters(iEvent, iSetup);" << std::endl;
  //if(fillsc_) FillSuperClusters(iEvent, iSetup);
  //std::cout << "FillTruth(iEvent, iSetup);" << std::endl;
  //if(type_ == "MC") FillTruth(iEvent, iSetup);
	
  //std::cout << "FillTipLipIp(iEvent, iSetup);" << std::endl;
  //if(!aod_) 
  //FillTipLipIp(iEvent, iSetup);
	
  //std::cout << "mytree_->Fill();" << std::endl;
  mytree_->Fill();
	
} // analyze

// ====================================================================================
void SimpleNtpleAOD::FillEvent (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  nEvent = iEvent.id().event();
  nRun   = iEvent.id().run();
  nLumi  = iEvent.luminosityBlock();

  if(PrintDebug_) cout << "nRun=" << nRun << " | nEvent=" << nEvent << endl;

  // -----------------
  // Pile-up
  // -----------------
  if(type_ == "MC") {
    Handle<vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByLabel(PileupSrc_, PupInfo);
    for (vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand) {
	    
      _PU_N = cand->getPU_NumInteractions();
      //cout << " PU = "<< _PU_N << endl;
    } // loop on Pile up
  } // if MC

  // Rho/FastJet Correction
  Handle<double> rhoHandle, sigmaHandle;
  iEvent.getByLabel(RhoCorrection_, rhoHandle);
  iEvent.getByLabel(SigmaRhoCorrection_, sigmaHandle);
  _PU_rho   = *rhoHandle;
  _PU_sigma = *sigmaHandle;
	
  // 	cout << "Rho, Sigma: " << _rho << "   " << _sigma << endl;

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
void SimpleNtpleAOD::FillTrigger (const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
	
  // ----------------------------------------------
  //  Get HLT info
  // ----------------------------------------------
	
  // Get HLTTag when running
  //Handle<trigger::TriggerEvent> triggerEventHLT;
  //iEvent.getByLabel("hltTriggerSummaryAOD", triggerEventHLT);
  //cout << " HLT = " << triggerEventHLT.provenance()->processName()  << endl;
	
  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByLabel (HLTTag_,triggerResultsHandle);
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResultsHandle);
	  
  const int nTrig = 4;

  string nad_hlt_name[nTrig] = {"HLT_Activity_Ecal_SC","HLT_L1SingleEG5","HLT_L1SingleEG8","HLT_L1SingleEG12"};
  bool nad_hlt_trig[nTrig];

  for(int iTrig=0 ; iTrig<nTrig ; iTrig++) {
    trig_HLT_path[iTrig] = 0;
    nad_hlt_trig[iTrig] = false;
  }

  if(PrintDebug_HLT_) {
    cout << "triggerResultsHandle->size() = " << triggerResultsHandle->size() << endl
	 << "<----- HLT MENU ----->" << endl;
  }

  // LOOP Over Trigger Results
  strcpy(trig_fired_names,"*");

  for (int iHLT=0; iHLT<static_cast<int>(triggerResultsHandle->size()); iHLT++) {	

    if( string(triggerNames.triggerName(iHLT)).find("HLT_L1SingleEG")!= string::npos || 
	string(triggerNames.triggerName(iHLT)).find("HLT_Activity_Ecal")!= string::npos ) {
	
      if( PrintDebug_HLT_ ) {
	cout << " -- " << string(triggerNames.triggerName(iHLT)) ;
      }
      if (triggerResultsHandle->accept (iHLT)) {
	if( PrintDebug_HLT_ ) { cout << "   ==========> TRIGGERED !!!!!!!!!!! --" << endl; }
	for(int iTrig=0 ; iTrig<nTrig ; iTrig++) {
	  if( string( triggerNames.triggerName(iHLT) ).find( nad_hlt_name[iTrig] ) != string::npos) {
	    nad_hlt_trig[iTrig] = true;
	    trig_HLT_path[iTrig] = 1;
	  }
	}
      }
      else if(PrintDebug_HLT_) cout << endl;
    }

    if (triggerResultsHandle->accept (iHLT)) {
      trig_hltInfo[iHLT] = 1;
      if ( strlen(trig_fired_names) <= 4950) {
	const char* c_str();
	string hlt_string = triggerNames.triggerName(iHLT);
	strcat(trig_fired_names,hlt_string.c_str());
	strcat(trig_fired_names,"*");
      }
    }
    else {
      trig_hltInfo[iHLT] = 0;
    }
	
  } // loop on trigger results

  for(int iTrig=0 ; iTrig<nTrig ; iTrig++)
    if( nad_hlt_trig[iTrig] )
      trig_HLT_path[iTrig] = 1;

  /////////////////////////// 
  // Get TP data  (Nadir)  //
  ///////////////////////////

  // commented to treat non modif reco SingleElectron dataset
  // bool PrintDebug_ = true;
  //   if(PrintDebug_) std::cout << "" << endl;
  
  // ORIGINAL TP
  if( nadGetTP_ ) {
    if(PrintDebug_) cout << "create new ecal_tp pointer" << endl;
    //edm::Handle<EcalTrigPrimDigiCollection> tp;
    ecal_tp_ = new edm::Handle<EcalTrigPrimDigiCollection> ;

    if(PrintDebug_) cout << "..created. get by label the tp collection" << endl;

    iEvent.getByLabel(tpCollectionNormal_,*ecal_tp_);
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
      TPtowid_ = d_.id(); // const EcalTrigTowerDetId TPtowid
      if(PrintDebug_) cout << "got the tower id" << endl;
      _trig_tower_iphi[i] = TPtowid_.iphi() ;
      _trig_tower_ieta[i] = TPtowid_.ieta() ;
      if(PrintDebug_) cout << "got the ieta and iphi : " << TPtowid_.ieta() << TPtowid_.iphi() << endl;
      //_trig_tower_adc[i]  = (d[0].raw()&0xfff) ;  0xfff <-> TTF(3bits)+FG(1bit)+Et(8bits)
      _trig_tower_adc[i]  = (d_[0].raw()&0xff) ;  // 0xff  <-> Et(8bits)
      if(PrintDebug_) cout << "got the adc : " << (int)(d_[0].raw()&0xff) << endl;
      //if(_trig_tower_adc[i]>0)
      //cout << _trig_tower_adc[i] << "   " ;
      _trig_tower_sFGVB[i] = d_[0].sFGVB();       // 0=spike-like / 1=EM-like
      if(PrintDebug_) cout << "got the sFGVB : " << d_[0].sFGVB() << endl;
      //_trig_tower_sFGVB[i] = d[0].l1aSpike();
      //if(d[0].l1aSpike()!=0) cout << "sFGVB=" << d[0].l1aSpike() << endl;
    }
    if(PrintDebug_) cout << "finished looping" << endl;
  }

  // ZEROING-BY-HAND TP
  if( nadGetTP_Modif_ ) {
    ecal_tpM_ = new edm::Handle<EcalTrigPrimDigiCollection> ;
    iEvent.getByLabel(tpCollectionModif_,*ecal_tpM_);
    //std::cout << "TP Modif collection size=" << tpM.product()->size() << std::endl ;

    _trig_tower_N_modif = ecal_tpM_->product()->size(); 

    for (int i=0 ; i<_trig_tower_N_modif ; i++) {
      EcalTriggerPrimitiveDigi dM_ = (*(ecal_tpM_->product()))[i]; // EcalTriggerPrimitiveDigi dM
      TPtowidM_ = dM_.id(); // EcalTrigTowerDetId
      _trig_tower_iphi_modif[i] = TPtowidM_.iphi() ;
      _trig_tower_ieta_modif[i] = TPtowidM_.ieta() ;
      //_trig_tower_adc_modif[i]  = (dM[0].raw()&0xfff) ;
      _trig_tower_adc_modif[i]  = (dM_[0].raw()&0xff) ;
      //if(_trig_tower_adc_modif[i]>0)
      //cout << _trig_tower_adc_modif[i] << "   " ;
      _trig_tower_sFGVB_modif[i] = dM_[0].sFGVB(); // 0=spike-like / 1=EM-like
      //_trig_tower_sFGVB_modif[i] = dM[0].l1aSpike();
    }
  }

  // EMULATOR TPs
  if( nadGetTP_Emul_ ) {
    
    ecal_tpM_ = new edm::Handle<EcalTrigPrimDigiCollection> ;
    iEvent.getByLabel(tpEmulatorCollection_, *ecal_tpM_);
    //if (print_) std::cout<<"TPEmulator collection size="<<tpEmul.product()->size()<<std::endl ;
  
    _trig_tower_N_emul = ecal_tpM_->product()->size();

    for (int i=0 ; i<_trig_tower_N_emul ; i++) {
      EcalTriggerPrimitiveDigi dM_ = (*(ecal_tpM_->product()))[i]; //EcalTriggerPrimitiveDigi
      TPtowidM_ = dM_.id();
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
	//_trig_tower_sFGVB_emul[i][j] = d[j].l1aSpike(); 
	_trig_tower_sFGVB_emul[i][j] = dM_[j].sFGVB(); 
	if(showit)
	  cout << (dM_[j].raw()&0xff) << " " ;
      }
      if(showit)
	cout << endl;
    }
    
  }  

  //------------------------//
  // GET THE SPIKES (Nadir) //
  //------------------------//

  // geometry (used for L1 trigger)    
  //cout << "get the geometry" << endl;

  edm::ESHandle<CaloSubdetectorGeometry> theBarrelGeometry_handle;     
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
  //iSetup.get<IdealGeometryRecord>().get(eTTmap_);
  theBarrelGeometry_ = &(*theBarrelGeometry_handle);

  //cout << "starting getting the spike-like rechits " << endl;

  // channel status : old way to do it
  /*
  edm::ESHandle<EcalChannelStatus> pChannelStatus;
  iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
  const EcalChannelStatus *chaStatus = pChannelStatus.product();
  */
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
  //cout << "looking" << endl;
  if (iEvent.getByLabel(EcalRecHitCollectionEB_, rechitsEB) ) {
    //cout << "looping" << endl;
    for ( rechitItr = rechitsEB->begin(); rechitItr != rechitsEB->end(); ++rechitItr ) {	
      //cout << "i=" << i << endl;
      if(i>=5000) {
	cout << "more than 5000 spikes in run " << iEvent.id().run() << " , event " << iEvent.id().event() << endl;
	break;
      }
      id = rechitItr->id();	      
      //sev = EcalSeverityLevelAlgo::severityLevel( id, *rechitsEB, *chaStatus );
      //sev = EcalSeverityLevelAlgo::severityLevel(id,*rechitsEB,*chaStatus, 5., 
      //EcalSeverityLevelAlgo::kSwissCross,0.95) ;
      sev = 0;
      sev = sl->severityLevel( id, *rechitsEB );
      //if(sev>0) cout << "got sev=" << sev << endl;

      if (sev >= 1) {
	//cout << " severity youpi " << endl;

	thetaHit =  (theBarrelGeometry_->getGeometry(id)->getPosition()).theta();
	etaHit =  (theBarrelGeometry_->getGeometry(id)->getPosition()).eta();
	phiHit =  (theBarrelGeometry_->getGeometry(id)->getPosition()).phi();
	const EcalTrigTowerDetId towid = id.tower();
	flag = 0;
	flag = rechitItr->recoFlag();
	if( PrintDebug_ && false )
	  cout << "-- RecHit : "
	       << "flag=" << flag
	       << " | sev=" << sev
	       << " | thetaHit="  << thetaHit 
	       << " | phiHit=" << phiHit
	       << " | etaHit=" << etaHit
	       << " | E="      << rechitItr->energy()
	       << " | Et="     << (rechitItr->energy())*sin(thetaHit)
	       << endl;

	//cout << "is gonna fill the variables" << endl;

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
  }
  if(i>=5000) spike_N = 5000;
  else spike_N = i+1;


  // ----------------------------------
  //  Path from list given in .py file
  // ----------------------------------
  UInt_t trigger_size = triggerResultsHandle->size();
  int passEleTrigger  = 0;
  int passMuonTrigger = 0;
	
  // Electron Triggers
  for(int ipath=0;ipath< (int) HLT_ElePaths_.size();ipath++) {
    //cout << " i = " << ipath << " trigger = " << HLT_Paths_[ipath] << endl;
    UInt_t trigger_position = triggerNames.triggerIndex(HLT_ElePaths_[ipath]); //hltpath_);
    if (trigger_position < trigger_size) passEleTrigger = (int)triggerResultsHandle->accept(trigger_position);
    if (passEleTrigger==1) _trig_isEleHLTpath = 1;
  } // for loop on HLT Elepaths
	
  // Muon Triggers
  for(int ipath=0;ipath< (int) HLT_MuonPaths_.size();ipath++) {
    //cout << " i = " << ipath << " trigger = " << HLT_Paths_[ipath] << endl;
    UInt_t trigger_position = triggerNames.triggerIndex(HLT_MuonPaths_[ipath]); //hltpath_);
    if (trigger_position < trigger_size) passMuonTrigger = (int)triggerResultsHandle->accept(trigger_position);
    if (passMuonTrigger==1) _trig_isMuonHLTpath = 1;
  } // for loop on HLT Muonpaths
  
  
  // ----------------------
  //  get L1 EM candidate
  // ----------------------
  
  // --- CURRENT BUNCH CROSSING --- //////////////////////////////////////////////////////////////////

  edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl ;
  edm::Handle< l1extra::L1EmParticleCollection > emIsolColl ;
  edm::Handle< l1extra::L1EmParticleCollection > emNonisolColl_M ;
  edm::Handle< l1extra::L1EmParticleCollection > emIsolColl_M ;  
  
  if( !nadGetL1M_ ) {      
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
  
  ///// STANDARD COLLECTION ALONE
  
  // Isolated candidates
  _trig_L1emIso_N = emIsolColl->size();
  if(PrintDebug_) cout << "N L1 candidate iso : " << _trig_L1emIso_N << endl;
  int counter = 0;
  for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl->begin(); emItr != emIsolColl->end() ;++emItr) {
    if( !aod_ ) {
      _trig_L1emIso_ieta[counter] = emItr->gctEmCand()->regionId().ieta();
      _trig_L1emIso_iphi[counter] = emItr->gctEmCand()->regionId().iphi();
      _trig_L1emIso_rank[counter] = emItr->gctEmCand()->rank(); // 1 adc = 1 GeV (2011 setting)
    }
    else {
      //_trig_L1emIso_ieta[counter] = getRCT_ieta( emItr->eta() );
      //_trig_L1emIso_iphi[counter] = getRCT_iphi( emItr->phi() );
      _trig_L1emIso_rank[counter] = (int)(emItr->et());
    }
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
      if( !aod_ ) {
	_trig_L1emNonIso_ieta[counter] = emItr->gctEmCand()->regionId().ieta();
	_trig_L1emNonIso_iphi[counter] = emItr->gctEmCand()->regionId().iphi();
	_trig_L1emNonIso_rank[counter] = emItr->gctEmCand()->rank(); // ET in ADC count... 1 ADC count = 0.5 GeV
      }
      else {
	//_trig_L1emNonIso_ieta[counter] = getRCT_ieta( emItr->eta() );
	//_trig_L1emNonIso_iphi[counter] = getRCT_iphi( emItr->phi() );
	_trig_L1emNonIso_rank[counter] = (int)(emItr->et());
      }
      _trig_L1emNonIso_eta[counter]    = emItr->eta();
      _trig_L1emNonIso_phi[counter]    = emItr->phi();
      _trig_L1emNonIso_energy[counter] = emItr->energy();
      _trig_L1emNonIso_et[counter]     = emItr->et();
      counter++;
    } // for loop on Non Iso cand
	  
    
    ///// MODIFIED COLLECTION IF ASKED
    if( nadGetL1M_ ) {

      // Isolated candidates
      _trig_L1emIso_N_M = emIsolColl_M->size();
      if(PrintDebug_) cout << "_trig_L1emIso_N_M =" << _trig_L1emIso_N_M << endl;
      counter=0;
      for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl_M->begin(); 
	   emItr != emIsolColl_M->end() ;++emItr) {
	// Used by Clemy
	_trig_L1emIso_ieta_M[counter] = emItr->gctEmCand()->regionId().ieta();
	_trig_L1emIso_iphi_M[counter] = emItr->gctEmCand()->regionId().iphi();
	_trig_L1emIso_rank_M[counter] = emItr->gctEmCand()->rank(); 
	// ET in ADC count... 1 ADC count = 0.5 GeV
	// From Trigger twiki
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
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////

    
    // --- PRE- AND POST-FIRING ---
    
    if(PrintDebug_) cout << "Pre-Post firing" << endl;

    edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
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
    _trig_postL1emIso_N = counter;
    
    // } //if AOD
	
	
  // ----------------------
  //  get HLT EM candidate
  // ----------------------
  edm::Handle<trigger::TriggerEvent> trigEvent;
  iEvent.getByLabel(triggerEventTag_, trigEvent);
	
  const Int_t N_filter(trigEvent->sizeFilters());
  std::vector<Int_t> ID_filter; 
	
  // Print Official Filters
  //for(int ifi=0;ifi<N_filter;ifi++) {
  //cout << "filter tag " << ifi << " = " << trigEvent->filterTag(ifi) << endl;
  //} // for loop on filters
	
  int hlt_counter = 0;

  // Loop on user's Filters
  for(int itrig=0;itrig< (int) HLT_Filters_.size();itrig++) {
		
		
    ID_filter.push_back(trigEvent->filterIndex(HLT_Filters_[itrig])); 
		
    const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());
    if( ID_filter[itrig] <  N_filter) { // !!! To be checked !!! trigEvent->size() ) {
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
	

} // end of FillTrigger

// ====================================================================================
void SimpleNtpleAOD::FillMET (const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
void SimpleNtpleAOD::FillEle(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  edm::Handle<reco::GsfElectronCollection> EleHandle ;
  iEvent.getByLabel (EleTag_.label(),EleHandle) ;
	
  edm::Handle<reco::PFCandidateCollection> PfEleHandle;
  iEvent.getByLabel("particleFlow", PfEleHandle);
	
  std::vector<edm::Handle<edm::ValueMap<float> > > eleIdCutHandles(9) ;
  iEvent.getByLabel  (EleID_VeryLooseTag_ , eleIdCutHandles[0]) ;
  iEvent.getByLabel  (EleID_LooseTag_ , eleIdCutHandles[1]) ;
  iEvent.getByLabel  (EleID_MediumTag_ , eleIdCutHandles[2]) ;
  iEvent.getByLabel  (EleID_TightTag_ , eleIdCutHandles[3]) ;
  iEvent.getByLabel  (EleID_SuperTightTag_ , eleIdCutHandles[4]) ; 
  iEvent.getByLabel  (EleID_HyperTight1Tag_ , eleIdCutHandles[5]) ;
  iEvent.getByLabel  (EleID_HyperTight2Tag_ , eleIdCutHandles[6]) ;
  iEvent.getByLabel  (EleID_HyperTight3Tag_ , eleIdCutHandles[7]) ;
  iEvent.getByLabel  (EleID_HyperTight4Tag_ , eleIdCutHandles[8]) ;
	

  //std::cout << " FillEle recoBeamSpotHandle " << std::endl;

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  iEvent.getByType(recoBeamSpotHandle) ;
  const reco::BeamSpot bs = *recoBeamSpotHandle ;
	
  //std::cout << " FillEle calo topology " << std::endl;

  //calo topology
  const CaloTopology * topology ;
  const EcalChannelStatus *chStatus ;
  edm::Handle< EcalRecHitCollection > reducedEBRecHits;
  edm::Handle< EcalRecHitCollection > reducedEERecHits;
	
  unsigned long long cacheIDTopo_=0;
  edm::ESHandle<CaloTopology> theCaloTopo;
  if (cacheIDTopo_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
    cacheIDTopo_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
    iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  }
  topology = theCaloTopo.product() ;
	
  edm::ESHandle<EcalChannelStatus> pChannelStatus;
  iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
  chStatus = pChannelStatus.product();

  //for42x	
  unsigned long long cacheSevLevel = 0;
  edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
  if(cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()){
    cacheSevLevel = iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
  }
  const EcalSeverityLevelAlgo* sl=sevLevel.product();
  //std::cout << " FillEle geometry (used for L1 trigger) " << std::endl;

  // geometry (used for L1 trigger)                                                                                                                
  edm::ESHandle<CaloSubdetectorGeometry> theEndcapGeometry_handle, theBarrelGeometry_handle;
	
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap",theEndcapGeometry_handle);
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
	
  iSetup.get<IdealGeometryRecord>().get(eTTmap_);
  theEndcapGeometry_ = &(*theEndcapGeometry_handle);
  theBarrelGeometry_ = &(*theBarrelGeometry_handle);
	
  //std::cout << " FillEle reduced rechits " << std::endl;
  // reduced rechits
  if(!aod_){
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEE"), reducedEERecHits ) ;
  }
  else{
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEE"), reducedEERecHits ) ;
  }

  edm::Handle<reco::TrackCollection> tracks_h;
  iEvent.getByLabel("generalTracks", tracks_h);
	
  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel(dcsTag_, dcsHandle);
  double evt_bField;
  // need the magnetic field
  //
  // if isData then derive bfield using the
  // magnet current from DcsStatus
  // otherwise take it from the IdealMagneticFieldRecord
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
	
  //std::cout << " FillEle Beam Spot Information " << std::endl;

  // Beam Spot Information
  BS_x = bs.position().x();
  BS_y = bs.position().y();
  BS_z = bs.position().z();
	
  BS_dz = bs.sigmaZ();
  BS_dydz = bs.dydz();
  BS_dxdz = bs.dxdz();
	
  BS_bw_x = bs.BeamWidthX();
  BS_bw_y = bs.BeamWidthY();
	

  //std::cout << " FillEle Masked Towers " << std::endl;
	
  // -----------------------------
  // Masked Towers
  // -----------------------------
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
	  //cout << "ECAL masked channel: RCT crate " << i << " iphi " << j <<" ieta " <<k <<endl;                                                                      
	  _trig_iMaskedRCTeta[n0MaskedRCT]=k;
	  _trig_iMaskedRCTphi[n0MaskedRCT]=j;
	  _trig_iMaskedRCTcrate[n0MaskedRCT]=i;
	  n0MaskedRCT++;
	}
  _trig_nMaskedRCT = n0MaskedRCT;
	
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
	
  //std::cout << " FillEle Seeds collection " << std::endl;
	
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

  //std::cout << " FillEle Get MC information " << std::endl;

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

  //std::cout << " FillEle Get TrackingParticles info " << std::endl;

  // ----------------------------------------------
  //  Get TrackingParticles info   
  // ----------------------------------------------
	
  /*
  // Get simulated       
  edm::Handle<TrackingParticleCollection> simCollection;
  if(type_ == "MC" && simulation_ == true) iEvent.getByLabel(TkPTag_, simCollection);
	
  // Get reconstructed
  edm::Handle<edm::View<reco::Track> >  recCollection;
  iEvent.getByLabel("electronGsfTracks", recCollection);
	
  edm::RefToBaseVector<reco::Track> tc(recCollection);
  for (unsigned int j=0; j<recCollection->size();j++)
  tc.push_back(edm::RefToBase<reco::Track>(recCollection,j));
	
	
	
  // Get associator 
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  if(type_ == "MC" && simulation_ == true) iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", theHitsAssociator);
	
  TrackAssociatorBase* theAssociatorByHits;
  reco::RecoToSimCollection recoToSim;
  if(type_ == "MC" && simulation_ == true) {
  theAssociatorByHits = (TrackAssociatorBase*)theHitsAssociator.product();
  recoToSim = theAssociatorByHits->associateRecoToSim(recCollection, simCollection,&iEvent);
  //     recoToSim = theAssociatorByHits->associateRecoToSim(tc, simCollection,&iEvent);
		
  if(recoToSim.size() > 0) std::cout << "size = " << recoToSim.size() << std::endl;
  }
  */	

	  
  // Get OLD HZZ isolation (commented... NOT USED ANYMORE)
  //edm::Handle<edm::ValueMap<float> > isoTkelemap;
  //if(!aod_) iEvent.getByLabel(EleIso_TdrHzzTkMapTag_, isoTkelemap);
  //edm::Handle<edm::ValueMap<float> > isoHadelemap;
  //if(!aod_) iEvent.getByLabel(EleIso_TdrHzzHcalMapTag_, isoHadelemap);
	  
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
  //EgammaTowerIsolation * towerIso1_015_0  = new EgammaTowerIsolation(0.15, 0., 0.,  1, towersH_->product()) ;
  //EgammaTowerIsolation * towerIso2_015_0  = new EgammaTowerIsolation(0.15, 0., 0.,  2, towersH_->product()) ;
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
	  
  // Get e/g for HZZ Isolation (e/g isolation optimized for HZZ)
  edm::Handle<edm::ValueMap<double> > egmisoTkelemap;
  iEvent.getByLabel(EleIso_Eg4HzzTkMapTag_ , egmisoTkelemap);
  edm::Handle<edm::ValueMap<double> > egmisoEcalelemap;
  iEvent.getByLabel(EleIso_Eg4HzzEcalMapTag_ , egmisoEcalelemap);
  edm::Handle<edm::ValueMap<double> > egmisoHcalelemap;
  iEvent.getByLabel(EleIso_Eg4HzzHcalMapTag_ , egmisoHcalelemap);
  float deta = -20.;
  float dphi = -20.;
	
  if(EleHandle->size() < 10 ){ ele_N = EleHandle->size(); }
  else {ele_N = 10;}
  TClonesArray &electrons = *m_electrons;
  int counter = 0;
	
  //std::cout << " FillEle Loop on Electrons " << std::endl;

  // ----------------------------------------------
  //  Loop on Electrons
  // ----------------------------------------------
  int nTow=0;
  int nReg=0;
	
  for(int i=0; i< ele_N; i++){
		
    edm::Ref<reco::GsfElectronCollection> electronEdmRef(EleHandle,i);
    setMomentum (myvector, (*EleHandle)[i].p4());
    new (electrons[counter]) TLorentzVector (myvector);
		
    ele_eidVeryLoose[counter] = (*(eleIdCutHandles[0]))[electronEdmRef]; 
    ele_eidLoose[counter] = (*(eleIdCutHandles[1]))[electronEdmRef]; 
    ele_eidMedium[counter] = (*(eleIdCutHandles[2]))[electronEdmRef]; 
    ele_eidTight[counter] = (*(eleIdCutHandles[3]))[electronEdmRef]; 
    ele_eidSuperTight[counter] = (*(eleIdCutHandles[4]))[electronEdmRef]; 
    ele_eidHyperTight1[counter] = (*(eleIdCutHandles[5]))[electronEdmRef]; 
    ele_eidHyperTight2[counter] = (*(eleIdCutHandles[6]))[electronEdmRef]; 
    ele_eidHyperTight3[counter] = (*(eleIdCutHandles[7]))[electronEdmRef]; 
    ele_eidHyperTight4[counter] = (*(eleIdCutHandles[8]))[electronEdmRef]; 
		
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

    if(!aod_){
      ele_pin_mean[counter]    = (*EleHandle)[i].gsfTrack()->innerMomentum().R() ; 
      ele_pout_mean[counter]   = (*EleHandle)[i].gsfTrack()->outerMomentum().R(); 

      ele_pTin_mean[counter]   = (*EleHandle)[i].gsfTrack()->innerMomentum().Rho() ; 
      ele_pTout_mean[counter]  = (*EleHandle)[i].gsfTrack()->outerMomentum().Rho() ;
    }

    //std::cout << " FillEle Get SuperCluster Informations " << std::endl;

    // Get SuperCluster Informations
    reco::SuperClusterRef sclRef = (*EleHandle)[i].superCluster();
    math::XYZPoint sclPos = (*EleHandle)[i].superClusterPosition();
    if (!(*EleHandle)[i].ecalDrivenSeed() && (*EleHandle)[i].trackerDrivenSeed()) 
      sclRef = (*EleHandle)[i].pflowSuperCluster();
		
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
    //cout << "scl X Y Z : "<< sclPos.X() <<", "<<sclPos.Y()<<", "<<sclPos.Z()<<endl;
    //cout << "scl X Y Z : "<< sclX[counter] <<", "<<sclY[counter]<<", "<<sclZ[counter]<<endl;
	

    // NEW H/E
    //reco::SuperCluster EmSCCand = *isc;
    reco::SuperCluster EmSCCand = *sclRef;
    double HoE_00615_0  = towerIso1_00615_0->getTowerESum(&EmSCCand)  + towerIso2_00615_0->getTowerESum(&EmSCCand) ;
    double HoE_005_0  = towerIso1_005_0->getTowerESum(&EmSCCand)  + towerIso2_005_0->getTowerESum(&EmSCCand) ;
    double HoE_005_1  = towerIso1_005_1->getTowerESum(&EmSCCand)  + towerIso2_005_1->getTowerESum(&EmSCCand) ;
    double HoE_005_15 = towerIso1_005_15->getTowerESum(&EmSCCand) + towerIso2_005_15->getTowerESum(&EmSCCand) ;

    double HoE_01_0  = towerIso1_01_0->getTowerESum(&EmSCCand)  + towerIso2_01_0->getTowerESum(&EmSCCand) ;
    double HoE_01_1  = towerIso1_01_1->getTowerESum(&EmSCCand)  + towerIso2_01_1->getTowerESum(&EmSCCand) ;
    double HoE_01_15 = towerIso1_01_15->getTowerESum(&EmSCCand) + towerIso2_01_15->getTowerESum(&EmSCCand) ;

    //double HoE_015_0  = towerIso1_01_0->getTowerESum(&EmSCCand)  + towerIso2_01_0->getTowerESum(&EmSCCand) ;
    double HoE_015_1  = towerIso1_015_1->getTowerESum(&EmSCCand)  + towerIso2_015_1->getTowerESum(&EmSCCand) ;
    double HoE_015_15 = towerIso1_015_15->getTowerESum(&EmSCCand) + towerIso2_015_15->getTowerESum(&EmSCCand) ;
		

    HoE_00615_0  /=         sclRef->energy() ;
    HoE_005_0  /= 	sclRef->energy() ;     
    HoE_005_1  /= 	sclRef->energy() ;     
    HoE_005_15 /= 	sclRef->energy() ;     
    HoE_01_0   /= 	sclRef->energy() ;     
    HoE_01_1   /= 	sclRef->energy() ;     
    HoE_01_15  /= 	sclRef->energy() ;     
    //HoE_015_0  /= 	sclRef->energy() ;     
    HoE_015_1  /= 	sclRef->energy() ;     
    HoE_015_15 /= 	sclRef->energy() ; 

    _ele_he_00615_0[counter]  = HoE_00615_0 ;
    _ele_he_005_0[counter]  = HoE_005_0 ;
    _ele_he_005_1[counter]  = HoE_005_1 ;	
    _ele_he_005_15[counter] = HoE_005_15 ;
    _ele_he_01_0[counter]   = HoE_01_0 ;
    _ele_he_01_1[counter]   = HoE_01_1 ;	
    _ele_he_01_15[counter]  = HoE_01_15 ;
    //_ele_he_015_0[counter]  = HoE_01_0 ;
    _ele_he_015_1[counter]  = HoE_015_1 ;	
    _ele_he_015_15[counter] = HoE_015_15 ;
		
	
    //	cout<<"debug 0 "<<counter<<endl;
    // 		total effective uncertainty on energy in GeV
    if (funcbase_) {
      ele_sclErr[counter] = funcbase_->getValue(*sclRef, 0); //
      //	cout<<"debug 1 "<<counter<<endl;
      // positive uncertainty
      ele_sclErr_pos[counter] = funcbase_->getValue(*sclRef, 1);
      //cout<<"debug 2 "<<counter<<endl;
      // negative uncertainty
      ele_sclErr_neg[counter] = funcbase_->getValue(*sclRef, -1);
      //	cout<<"debug 3 "<<counter<<endl;
      ele_trErr[counter]=(*EleHandle)[i].trackMomentumError();
      //cout<<"debug 4 "<<counter<<endl;
      //for42X           
      //			ele_momErr[counter]=(*EleHandle)[i].electronMomentumError();
      ele_momErr[counter]=(*EleHandle)[i].p4Error(GsfElectron::P4_COMBINATION);
      //cout<<"debug 5 "<<counter<<endl;

      //  cout<<" ele_trErr[counter]/ele_pin_mode[counter] = "<<ele_trErr[counter]/ele_pin_mode[counter]<<" ele_sclErr[counter]/ele_sclE[counter] = "<<ele_sclErr[counter]/ele_sclE[counter]<<endl;
      if (!(*EleHandle)[i].ecalDrivenSeed()) { /// no change if not ecaldriven
	// cout<<" no change not ecal driven"<<endl;
	ele_newmom[counter] = (*EleHandle)[i].p4().P();
	ele_newmomErr[counter] = ele_momErr[counter];
      }
      else { //if ecal driven special care for large errors
	if (ele_trErr[counter]/ele_pin_mode[counter] > 0.5 && ele_sclErr[counter]/ele_sclE[counter] <= 0.5) { //take E if sigmaE/E <=0.5 and sigmaP/P >0.5
	  ele_newmom[counter] = ele_sclE[counter];    ele_newmomErr[counter] = ele_sclErr[counter];
	  // cout<<" E choice new mom= "<<ele_newmom[counter]<<" new err= "<<ele_newmomErr[counter]<<endl;
	}
	else if (ele_trErr[counter]/ele_pin_mode[counter] <= 0.5 && ele_sclErr[counter]/ele_sclE[counter] > 0.5){//take P if sigmaE/E > 0.5 and sigmaP/P <=0.5
	  ele_newmom[counter] = ele_pin_mode[counter];  ele_newmomErr[counter] = ele_trErr[counter];
	  // cout<<" P choice new mom= "<<ele_newmom[counter]<<" new err= "<<ele_newmomErr[counter]<<endl;
	}
	else if (ele_trErr[counter]/ele_pin_mode[counter] > 0.5 && ele_sclErr[counter]/ele_sclE[counter] > 0.5){//take the lowest sigma/value if sigmaE/E >0.5 and sigmaP/P >0.5
	  if (ele_trErr[counter]/ele_pin_mode[counter] < ele_sclErr[counter]/ele_sclE[counter]) {
	    ele_newmom[counter] = ele_pin_mode[counter]; ele_newmomErr[counter] = ele_trErr[counter];
	    //	cout<<" P choice new mom= "<<ele_newmom[counter]<<" new err= "<<ele_newmomErr[counter]<<endl;
	  }
	  else{
	    ele_newmom[counter] = ele_sclE[counter]; ele_newmomErr[counter] = ele_sclErr[counter];
	    //	cout<<" E choice new mom= "<<ele_newmom[counter]<<" new err= "<<ele_newmomErr[counter]<<endl;
	  }
	}
	else { // if sigmaE/E <= 0.5 and sigmaP/P <=0.5 no change
	  // cout<<" no change"<<endl;
	  ele_newmom[counter] = (*EleHandle)[i].p4().P();
	  ele_newmomErr[counter] = ele_momErr[counter];
	}
      }
    }
    //	else  cout<<"no function base for ecal errors"<<endl;
    //	cout<<" new mom= "<< ele_newmom[counter]<<" new error= "<<ele_newmomErr[counter] <<endl;
		
		
    ele_tr_atcaloX[counter] = (*EleHandle)[i].trackPositionAtCalo().x();
    ele_tr_atcaloY[counter] = (*EleHandle)[i].trackPositionAtCalo().y();
    ele_tr_atcaloZ[counter] = (*EleHandle)[i].trackPositionAtCalo().z();
    //cout<<"track @calo: "<<(*EleHandle)[i].trackPositionAtCalo().x()<<" "<<(*EleHandle)[i].trackPositionAtCalo().y()<<" "<<(*EleHandle)[i].trackPositionAtCalo().z()<<endl;

    if(!aod_){
      ele_firsthit_X[counter] = (*EleHandle)[i].gsfTrack()->innerPosition().x();
      ele_firsthit_Y[counter] = (*EleHandle)[i].gsfTrack()->innerPosition().y();
      ele_firsthit_Z[counter] = (*EleHandle)[i].gsfTrack()->innerPosition().z();
    }
    // cout<<"first hit: "<<firsthit_X[counter]<<" "<<firsthit_Y[counter]<<" "<<firsthit_Z[counter]<<endl;
		
		

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
    ele_mva[counter]   = (*EleHandle)[i].mva() ;
		
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
    //} // if AOD

    //std::cout << " FillEle Get  Isolation variables " << std::endl;

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
    double HcalIso_00615Ring03_0  = (towerIso1_00615Ring03_0->getTowerEtSum(&((*EleHandle)[i]))  +
				     towerIso2_00615Ring03_0->getTowerEtSum(&((*EleHandle)[i])) );
    double HcalIso_005Ring03_0  = (towerIso1_005Ring03_0->getTowerEtSum(&((*EleHandle)[i]))  +
				   towerIso2_005Ring03_0->getTowerEtSum(&((*EleHandle)[i])) );
    double HcalIso_0Ring03_0  = (towerIso1_0Ring03_0->getTowerEtSum(&((*EleHandle)[i]))  +
				 towerIso2_0Ring03_0->getTowerEtSum(&((*EleHandle)[i])) );

    ele_hcalDepth1plus2TowerSumEt_00615dr03[counter] = HcalIso_00615Ring03_0;
    ele_hcalDepth1plus2TowerSumEt_005dr03[counter] = HcalIso_005Ring03_0;
    ele_hcalDepth1plus2TowerSumEt_0dr03[counter] = HcalIso_0Ring03_0;


    //std::cout << " FillEle Get from HZZ isolation  " << std::endl;
    //from HZZ isolation 
		
    //SumPt
    //if(!aod_){
    //ele_tkSumPtTdrHzz_dr025[counter]       = (*isoTkelemap)[electronEdmRef] ;
    //ele_hcalSumEtTdrHzz_dr02[counter]      = (*isoHadelemap)[electronEdmRef];
    ele_tkSumPtEg4Hzz_dr03[counter]        = (*egmisoTkelemap)[electronEdmRef];
    ele_ecalSumEtEg4Hzz_dr03[counter]      = (*egmisoEcalelemap)[electronEdmRef];
    ele_hcalSumEtEg4Hzz_dr04[counter]      = (*egmisoHcalelemap)[electronEdmRef];
    //}
    //SumPt/Pt normalization
		
    double pT_new = (*EleHandle)[i].p4().Pt(); //* ele_newmom[counter] / (*EleHandle)[i].p4().P();

    //if(!aod_){
    //ele_tkSumPtoPtTdrHzz_dr025[counter]    = (*isoTkelemap)[electronEdmRef]/pT_new;
    //ele_hcalSumEtoPtTdrHzz_dr02[counter]   = (*isoHadelemap)[electronEdmRef]/pT_new;
    ele_tkSumPtoPtEg4Hzz_dr03[counter]     = (*egmisoTkelemap)[electronEdmRef]/pT_new;
    ele_ecalSumEtoPtEg4Hzz_dr03[counter]   = (*egmisoEcalelemap)[electronEdmRef]/pT_new;
    ele_hcalSumEtoPtEg4Hzz_dr04[counter]   = (*egmisoHcalelemap)[electronEdmRef]/pT_new;
    //}
    //		std::cout << " FillEle Conversion Removal " << std::endl;


    ele_ECAL_fbrem[counter]   = sclRef->phiWidth()/sclRef->etaWidth();
    //cout<< "ecalfbrem = "<<ele_ECAL_fbrem[counter]<<" fbrem= "<<ele_fbrem[counter]<<endl;
		
    // pflow combinaison
    ele_PFcomb[counter] = 0.;
    ele_PFcomb_Err[counter] = 0.;
    std::vector<reco::PFCandidate> candidates = (*PfEleHandle.product());
    for (std::vector<reco::PFCandidate>::iterator it = candidates.begin(); it != candidates.end(); ++it)   {
      reco::PFCandidate::ParticleType type = (*it).particleId();
      // here you can ask for particle type, mu,e,gamma
      if ( type == reco::PFCandidate::e) {
	// gsfElectronRef pas encore implémentée, utilisation de gsfTrackRef instead
	//if (!(*it).gsfElectronRef().isNull() && ((*it).gsfElectronRef()->p4().e() == (*EleHandle)[i].p4().e())){
	//FIXME !aod_ added for DATA 
	if (!(*it).gsfTrackRef().isNull() && ((*it).gsfTrackRef()->p() == (*EleHandle)[i].gsfTrack()->p())){
	  //std::cout << "[E-p combi] found corresponding Gsfelectron " << std::endl;
	  //std::cout << "[E-p combi] EG comb " << (*EleHandle)[i].p4().P() << " PF comb "<< (*it).energy()<< std::endl;
	  ele_PFcomb[counter] = (*it).energy();
	}
      }
    }   
    //ele_PFcomb[counter]   = (*EleHandle)[i].p4(GsfElectron::P4_PFLOW_COMBINATION).P();
    //cout<<" PF comb = "<< ele_PFcomb[counter] << " EG comb= "<<(*EleHandle)[i].p4().P()<<endl;
    ele_PFcomb_Err[counter]   =(*EleHandle)[i].p4Error(GsfElectron::P4_PFLOW_COMBINATION);
    //cout<<" PF comberr = "<< ele_PFcomb_Err[counter]<<endl;
    if (!(*EleHandle)[i].pflowSuperCluster().isNull()) ele_PF_SCenergy[counter]   = (*EleHandle)[i].pflowSuperCluster()->energy();
    //cout<<" PF SC energy = "<< ele_PF_SCenergy[counter]<<endl;
                
    ele_PF_SCenergy_Err[counter]   = 0 ; //not implemented for the moment


    // Conversion Removal
    ele_isConversion[counter] = IsConv ((*EleHandle)[i]);
    ConversionFinder convFinder;

    //		std::cout << " FillEle Conversion >=36X  " << std::endl;

    // Conversion >=36X
    //if ( convFinder.isElFromConversion((*EleHandle)[i], tracks_h, evt_bField, 0.02, 0.02, 0.45)) ele_convFound[counter] = 1;
    ConversionInfo convInfo = convFinder.getConversionInfo((*EleHandle)[i], tracks_h, evt_bField);
    ele_conv_dist[counter] = convInfo.dist();
    ele_conv_dcot[counter] = convInfo.dcot();

    ele_convFound[counter] = 0;
    // Conversion 36X
    if ( convFinder.isFromConversion(convInfo, 0.02, 0.02) ) ele_convFound[counter] = 1; 
		

    //std::cout << " FillEle For L1 Trigger, Clemy's stuff   " << std::endl;

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
    for(int icc = 0; icc < 10; ++icc) {
      _ele_RCTetaVect[counter][icc] = -999;
      _ele_RCTphiVect[counter][icc] = -999;
      _ele_RCTetVect[counter][icc] = 0.;
      _ele_RCTL1isoVect[counter][icc] = -999;
      _ele_RCTL1nonisoVect[counter][icc] = -999;
      _ele_RCTL1isoVect_M[counter][icc] = -999;
      _ele_RCTL1nonisoVect_M[counter][icc] = -999;
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
				
	  //XTAL max                                                                                 
	  //if(myhit.energy() > rechitmax){
	  //rechitmax = myhit.energy();
	  //towidmax=towid;
	  //}
				
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
	      if(_trig_L1emIso_iphi[il1] == iRphi && _trig_L1emIso_ieta[il1] == iReta) _ele_RCTL1isoVect[counter][nReg] = _trig_L1emIso_rank[il1];
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

					
	    nReg++;
	  } // if newReg
				
				
	}//loop crystal
    }//loop cluster
		
    //double TTetmax  = 0.;
    //int iTTmax      = -1.;
    double TTetmax2 = 0.;
    int iTTmax2     = -1.;
    //     std::cout<<"ele: "<<counter<<std::endl;
    //     std::cout<<"SC ET, eta, phi: "<<sclEt[counter]<<" "<<sclEta[counter]<<" "<<sclPhi[counter]<<std::endl;
    //     std::cout<<" towers"<<std::endl;
    for (int iTow=0; iTow<nTow; ++iTow) {
      bool nomaskTT = true;
      //if(eTTetVect[counter][iTow] > 2) 	std::cout<<"TT ET, eta, phi: "<<eTTetVect[counter][iTow]<<" "<<eTTetaVect[counter][iTow]<<" "<<eTTphiVect[counter][iTow]<<std::endl;
      for (ittpg=towerMap.begin();ittpg!=towerMap.end();++ittpg) {
	if ((*ittpg).second > 0)
	  {
	    EcalTrigTowerDetId  ttId((*ittpg).first);
	    if(ttId.ieta() == _ele_TTetaVect[counter][iTow] && ttId.iphi() == _ele_TTphiVect[counter][iTow]) {
	      //if(eTTetVect[counter][iTow] > 2) std::cout<<"*** masked ***"<<std::endl;
	      nomaskTT=false;
	    }
	  }
      }//loop trigger towers
      //if(_ele_eTTetVect[counter][iTow] > TTetmax) {
      //iTTmax = iTow;
      //TTetmax = _ele_eTTetVect[counter][iTow];
      //}
      if(nomaskTT && _ele_TTetVect[counter][iTow] > TTetmax2) {
	iTTmax2 = iTow;
	TTetmax2 = _ele_TTetVect[counter][iTow];
      } // if nomask
    } // for loop on Towers
		
    //     std::cout<<" regions"<<std::endl;
    //     for (int iReg=0; iReg<nReg; ++iReg) {
    //       if(_ele_eRCTetVect[counter][iReg] > 2)     std::cout<<"RCT ET, eta, phi: "<<_ele_eRCTetVect[counter][iReg]<<" "<<_ele_eRCTetaVect[counter][iReg]<<" "<<_ele_eRCTphiVect[counter][iReg]<<std::endl;
    //     }
		
    //int TTetamax = getGCTRegionEta(_ele_eTTetaVect[counter][iTTmax]);
    //int TTphimax = getGCTRegionPhi(_ele_eTTphiVect[counter][iTTmax]);
    //_ele_eleRCTeta[counter]=TTetamax;
    //_ele_eleRCTphi[counter]=TTphimax;
		
    //for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
    //if(_trig_L1emIso_iphi[il1] == TTphimax && _trig_L1emIso_ieta[il1] == TTetamax) _ele_eleRCTL1iso[counter] = _trig_L1emIso_rank[il1];
    //}
    //for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
    //if(_trig_L1emNonIso_iphi[il1] == TTphimax && _trig_L1emNonIso_ieta[il1] == TTetamax) _ele_eleRCTL1noniso[counter] = _trig_L1emNonIso_rank[il1];
    //}
		
		
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
            

    } // if iTTmax2
		
    //     std::cout<<"main TT eta, phi: "<<eTTetaVect[counter][iTTmax]<<" "<<eTTphiVect[counter][iTTmax]<<std::endl;
		
    // ---------------------------
    // For Charge, Clemy's stuff
    // ---------------------------
    if(!aod_)
      ele_expected_inner_hits[counter] = (*EleHandle)[i].gsfTrack()->trackerExpectedHitsInner().numberOfHits();
		
    //tkIso03Rel[counter] = (*EleHandle)[i].dr03TkSumPt()/(*EleHandle)[i].p4().Pt();
    //ecalIso03Rel[counter] = (*EleHandle)[i].dr03EcalRecHitSumEt()/(*EleHandle)[i].p4().Pt();
    //hcalIso03Rel[counter] = (*EleHandle)[i].dr03HcalTowerSumEt()/(*EleHandle)[i].p4().Pt();
		
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
		
    if(!aod_){
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
    }

    //std::cout << " Spikes analysis " << std::endl;

    // ------------------
    // Spikes analysis  
    // ------------------
    if ((*EleHandle)[i].isEB()) { /*spikes are in the barrel*/
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
      //for42X
      //			int sev = EcalSeverityLevelAlgo::severityLevel(id.first,*reducedRecHits,*chStatus, 5., EcalSeverityLevelAlgo::kSwissCross,0.95) ;
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
	//for42X
	//			  int sev = EcalSeverityLevelAlgo::severityLevel(id.first,*reducedRecHits,*chStatus, 5., EcalSeverityLevelAlgo::kSwissCross,0.95);
	int sev=sl->severityLevel(id.first,*reducedRecHits);
	if (sev > dummySev)
	  dummySev = sev ;
      }			  
      ele_severityLevelClusters[counter] = dummySev ;
      ele_outOfTimeClusters[counter] = dummyFlag ;
      ele_e2overe9[counter] = E2overE9( id.first,*reducedRecHits,5,5, true, true);
      //	cout<<"e1/e9= "<<ele_e1[counter]/ele_e33[counter]<<" e2/e9= "<<ele_e2overe9[counter]<<endl;
    }
    else {
      ele_e2overe9[counter] = 0;
      ele_severityLevelSeed[counter] = 0 ;
      ele_outOfTimeSeed[counter] = 0 ;	  
      ele_severityLevelClusters[counter] = 0 ;
      ele_outOfTimeClusters[counter] = 0 ;	  
    }
		
    //std::cout << " Ambigous Tracks " << std::endl;

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
      //		bool okeleFound_conv = false; 
      //		bool okeleFound = false; 
			
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
	  // 	for (HepMC::GenEvent::vertex_const_iterator vertIter = MCEvt->vertices_begin();
	  // 	     vertIter != MCEvt->vertices_end(); ++vertIter) {
				
	  //				HepMC::FourVector momentum = (*partIter)->momentum();
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
		  //for particle gun 
		  //idPDG_mother_conv = (*((*partIter)->production_vertex()->particles_begin(HepMC::parents)))->pdg_id();
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
				
	  //			    HepMC::FourVector momentum = (*partIter)->momentum();
	  int idTmpPDG = partIter->pdgId();  
	  //MC particle
	  //genPc = (*partIter);
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
		  //idPDG_mother_conv = partIter->mother()->pdgId();
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
				
	  //	}//end loop over vertex
	  //if (okeleFound) continue ;
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
			
    } // end if "MC"
		
    ++counter;  
  } // for loop on electrons
	
}

// ====================================================================================
void SimpleNtpleAOD::FillMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  // Beam spot
  //Handle<reco::BeamSpot> beamSpotHandle;
  //iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle);
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  iEvent.getByType(recoBeamSpotHandle) ;
  const reco::BeamSpot bs = *recoBeamSpotHandle ;
  
  // Muon Retrieving
  Handle<View<reco::Muon> > MuonHandle;
  iEvent.getByLabel(MuonTag_, MuonHandle);
  
  TClonesArray &muons = *m_muons;
  int mu_counter = 0;
  
  
  // Get HZZ muon isolation
  edm::Handle<edm::ValueMap<double> > isomuonmap;
  //  if(!aod_) 
  iEvent.getByLabel(MuonIso_HzzMapTag_, isomuonmap);
  edm::Handle<edm::ValueMap<double> > isoTkmuonmap;
  //if(!aod_) 
  iEvent.getByLabel(MuonIsoTk_HzzMapTag_, isoTkmuonmap);
  edm::Handle<edm::ValueMap<double> > isoEcalmuonmap;
  //if(!aod_) 
  iEvent.getByLabel(MuonIsoEcal_HzzMapTag_, isoEcalmuonmap);
  edm::Handle<edm::ValueMap<double> > isoHcalmuonmap;
  //if(!aod_) 
  iEvent.getByLabel(MuonIsoHcal_HzzMapTag_, isoHcalmuonmap);
  
  // ----------------------------------------------
  //  Loop over Muons
  // ----------------------------------------------
  _muons_N = MuonHandle->size();
  
  
  for (edm::View<reco::Muon>::const_iterator imuons=MuonHandle->begin(); imuons!=MuonHandle->end(); ++imuons) {  
    if(mu_counter>19) continue;
    
    // 4-vector
    //edm::Ref<reco::MuonCollection> muonEdmRef(MuonHandle,i);
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
    //if(!aod_) {
    edm::Ref<edm::View<reco::Muon> > muref(MuonHandle, mu_counter);
    _muons_hzzIso[mu_counter]     = (*isomuonmap)[muref]; 
    _muons_hzzIsoTk[mu_counter]   = (*isoTkmuonmap)[muref]; 
    _muons_hzzIsoEcal[mu_counter] = (*isoEcalmuonmap)[muref]; 
    _muons_hzzIsoHcal[mu_counter] = (*isoHcalmuonmap)[muref]; 
    //} // if AOD
    
    mu_counter++;
  } // for loop on muons
  if(mu_counter>9) { _muons_N = 10; cout << "Number of muons>100, muons_N set to 10" << endl;}
  
} // end of FillMuons

// ====================================================================================
void SimpleNtpleAOD::FillJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
		
    // 		_jets_calo_E[index_calo_jets]   = ijets->energy();
    // 		_jets_calo_pT[index_calo_jets]  = ijets->pt();
    // 		_jets_calo_px[index_calo_jets]  = ijets->px();
    // 		_jets_calo_py[index_calo_jets]  = ijets->py();
    // 		_jets_calo_pz[index_calo_jets]  = ijets->pz();
    // 		_jets_calo_eta[index_calo_jets] = ijets->eta();
    // 		_jets_calo_phi[index_calo_jets] = ijets->phi();
    // 		// will add more variables in the future...
		
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

// ====================================================================================
void SimpleNtpleAOD::FillSuperClusters(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  //std::cout << "FillSuperClusters  geometry " << std::endl;
  // geometry
  const CaloGeometry * geometry;
  unsigned long long cacheIDGeom = 0;
  edm::ESHandle<CaloGeometry> theCaloGeom;
  if(cacheIDGeom!=iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
    cacheIDGeom = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
    iSetup.get<CaloGeometryRecord>().get(theCaloGeom);
  }
  geometry = theCaloGeom.product() ;
  
  //std::cout << "FillSuperClusters  TPGTowerStatus " << std::endl;
  
  // TPGTowerStatus
  edm::ESHandle<EcalTPGTowerStatus> theEcalTPGTowerStatus_handle;
  iSetup.get<EcalTPGTowerStatusRcd>().get(theEcalTPGTowerStatus_handle);
  const EcalTPGTowerStatus * ecaltpgTowerStatus=theEcalTPGTowerStatus_handle.product();
  
  const EcalTPGTowerStatusMap &towerMap=ecaltpgTowerStatus->getMap();
  EcalTPGTowerStatusMapIterator  it;
  
  //for42x	
  unsigned long long cacheSevLevel = 0;
  edm::ESHandle<EcalSeverityLevelAlgo> sevLevel;
  if(cacheSevLevel != iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier()){
    cacheSevLevel = iSetup.get<EcalSeverityLevelAlgoRcd>().cacheIdentifier();
    iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevel);
  }
  const EcalSeverityLevelAlgo* sl=sevLevel.product();

  // for H/E
  //if (towerIso1_) delete towerIso1_ ; towerIso1_ = 0 ;
  //if (towerIso2_) delete towerIso2_ ; towerIso2_ = 0 ;
  //if (towersH_) delete towersH_ ; towersH_ = 0 ;
  towersH_ = new edm::Handle<CaloTowerCollection>() ;
  if (!iEvent.getByLabel(hcalTowers_,*towersH_))
    { edm::LogError("ElectronHcalHelper::readEvent")<<"failed to get the hcal towers of label "<<hcalTowers_ ; }
  towerIso1_ = new EgammaTowerIsolation(hOverEConeSize_,0.,hOverEPtMin_,1,towersH_->product()) ;
  towerIso2_ = new EgammaTowerIsolation(hOverEConeSize_,0.,hOverEPtMin_,2,towersH_->product()) ;
  
  //std::cout << "FillSuperClusters  Beam Spot " << std::endl;
  
  // Beam Spot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  iEvent.getByType(recoBeamSpotHandle) ;
  const reco::BeamSpot bs = *recoBeamSpotHandle ;
  
  //std::cout << "FillSuperClusters  Isolation " << std::endl;	

  // Isolation
  edm::Handle<TrackCollection> ctfTracksH;  
  iEvent.getByLabel("generalTracks", ctfTracksH); // ctfTracks_
  //   //get the tracks
  //   edm::Handle<reco::TrackCollection> tracks;
  //   e.getByLabel(trackInputTag_,tracks);
  //   if(!tracks.isValid()) {
  //     return;
  //   }
  //   const reco::TrackCollection* trackCollection = tracks.product();
  
  //std::cout << "FillSuperClusters  Iso Track " << std::endl;
  
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
  
  //  double intRadiusBarrel = 0.015; 
  //   double intRadiusEndcap = 0.015; 
  //   double stripBarrel     = 0.015; 
  //   double stripEndcap     = 0.015; 
  //   double ptMin           = 0.7; 
  //   double maxVtxDist      = 0.2; 
  //   double drb             = 999999999.;  //  maxDrbTk 
  
  
  //std::cout << "FillSuperClusters  Iso HCAL " << std::endl;
  // Iso HCAL
  float egHcalIsoConeSizeOutSmall=0.3;
  //float egHcalIsoConeSizeOutLarge=0.4;
  int egHcalDepth1=1, egHcalDepth2=2; //float egHcalIsoConeSizeIn=intRadiusHcal_,egHcalIsoPtMin=etMinHcal_;
  double egHcalIsoConeSizeIn = 0.15;  //intRadiusHcal   = 0.15;
  double egHcalIsoPtMin      = 0.0;   //  etMinHcal 
  EgammaTowerIsolation hadDepth1Isolation03(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth1,towersH_->product()) ;
  EgammaTowerIsolation hadDepth2Isolation03(egHcalIsoConeSizeOutSmall,egHcalIsoConeSizeIn,egHcalIsoPtMin,egHcalDepth2,towersH_->product()) ;
  
  //std::cout << "FillSuperClusters  Iso ECAL " << std::endl;
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
  // for SpikeRemoval -- not in 361p4 -- 
  //int severityLevelCut           = 4;
  //double severityRecHitThreshold = 5.0;
  //double spikeIdThreshold        = 0.95;
  //string spId                    = "kSwissCrossBordersIncluded";  // ikeIdString 
  //for42x
  //EcalSeverityLevelAlgo::SpikeId spId = EcalSeverityLevelAlgo::kSwissCrossBordersIncluded;
  //
  //float extRadiusSmall=0.3, extRadiusLarge=0.4 ;
  float egIsoConeSizeOutSmall=0.3; //, egIsoConeSizeOutLarge=0.4;
  
  
  //std::cout << "FillSuperClusters  EB SuperCluster " << std::endl;
  
  // -----------------
  //  EB SuperCluster
  // -----------------
  // Retrieve SuperCluster
  edm::Handle<reco::SuperClusterCollection> sc_coll_EB;
  iEvent.getByLabel(edm::InputTag("correctedHybridSuperClusters"), sc_coll_EB);
  
  //cout << " size EB = " << sc_coll_EB->size() << endl;
  
  _sc_hybrid_N = sc_coll_EB->size();
  int index_sc = 0;
  
  //std::cout << "FillSuperClusters  SpikeRemoval " << std::endl;
  
  // Define stuff for SpikeRemoval
  const CaloTopology * topology ;
  const EcalChannelStatus *chStatus ;
  edm::Handle< EcalRecHitCollection > reducedEBRecHits;
  edm::Handle< EcalRecHitCollection > reducedEERecHits;
  
  unsigned long long cacheIDTopo_=0;
  edm::ESHandle<CaloTopology> theCaloTopo;
  if (cacheIDTopo_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
    cacheIDTopo_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
    iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  }
  topology = theCaloTopo.product() ;
  
  edm::ESHandle<EcalChannelStatus> pChannelStatus;
  iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
  chStatus = pChannelStatus.product();

  // reduced rechits
  if(!aod_){
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEE"), reducedEERecHits );
  }
  else{
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEE"), reducedEERecHits );
  }
  
  // For L1
  int nTow=0;
  int nReg=0;
  
  //std::cout << "FillSuperClusters  Loop EB SuperCluster " << std::endl;
  
  // --------------------------
  // Loop on SuperClusters EB
  // --------------------------
  for( reco::SuperClusterCollection::const_iterator isc=sc_coll_EB->begin(); isc!=sc_coll_EB->end(); isc++) {
    if(index_sc>24) continue;
    double R  = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y() +isc->z()*isc->z());
    double Rt = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y());
    
    _sc_hybrid_E[index_sc]   = isc->energy();
    _sc_hybrid_Et[index_sc]  = isc->energy()*(Rt/R);
    _sc_hybrid_Eta[index_sc] = isc->eta();
    _sc_hybrid_Phi[index_sc] = isc->phi();
    
    const EcalRecHitCollection * reducedRecHits = 0 ;
    reducedRecHits = reducedEBRecHits.product() ; 
    
    //seed cluster analysis
    const edm::Ptr<reco::CaloCluster> & seedCluster = isc->seed(); //(*EleHandle)[i].superCluster()->seed() ;  
    std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(),reducedRecHits);
    const EcalRecHit & rh = getRecHit(id.first,reducedRecHits);
    int flag = rh.recoFlag();   
    
    // Out of time
    if (flag == EcalRecHit::kOutOfTime) 
      _sc_hybrid_outOfTimeSeed[index_sc] = 1;   
    else 
      _sc_hybrid_outOfTimeSeed[index_sc] = 0;   
    
    // Severity Level
    //for42X
    //    int sev = EcalSeverityLevelAlgo::severityLevel(id.first,*reducedRecHits,*chStatus, 5., EcalSeverityLevelAlgo::kSwissCross,0.95) ;
    int sev=sl->severityLevel(id.first,*reducedRecHits);
    _sc_hybrid_severityLevelSeed[index_sc] = sev ;
    
    // Old SpikeRemoval e1/e9
    const reco::CaloCluster & seedCluster1 = *(isc->seed());
    _sc_hybrid_e1[index_sc]   = EcalClusterTools::eMax(seedCluster1,reducedRecHits)  ;
    _sc_hybrid_e33[index_sc]  = EcalClusterTools::e3x3(seedCluster1,reducedRecHits,topology)  ;
    
    // H/E
    reco::SuperCluster EmSCCand = *isc;
    double HoE = towerIso1_->getTowerESum(&EmSCCand) + towerIso2_->getTowerESum(&EmSCCand) ;
    HoE /= 	isc->energy() ;     
    _sc_hybrid_he[index_sc] = HoE ;
    
    // SigmaIetaIeta
    std::vector<float> localCovariances = EcalClusterTools::localCovariances(seedCluster1,reducedRecHits,topology) ;
    _sc_hybrid_sigmaietaieta[index_sc]  = sqrt(localCovariances[0]) ;
    
    //	std::cout << "FillSuperClusters  Iso Track " << std::endl;
    
    // Iso Track
    //ElectronTkIsolation tkIsolation03(extRadiusSmall,intRadiusBarrel,intRadiusEndcap,stripBarrel,stripEndcap,ptMin,maxVtxDist,drb,ctfTracksH.product(),bs.position()) ;
    //_sc_hybrid_tkSumPt_dr03[index_sc] = tkIsolation03.getPtTracks(isc);
    
    //	std::cout << "FillSuperClusters  Iso HCAL " << std::endl;
    
    // Iso HCAL
    _sc_hybrid_hcalDepth1TowerSumEt_dr03[index_sc] = hadDepth1Isolation03.getTowerEtSum(&EmSCCand);
    _sc_hybrid_hcalDepth2TowerSumEt_dr03[index_sc] = hadDepth2Isolation03.getTowerEtSum(&EmSCCand);
    
    //std::cout << "FillSuperClusters  Iso ECAL " << std::endl;
    
    // Iso ECAL
    EcalRecHitMetaCollection ecalBarrelHits(*reducedEBRecHits);
    //for42x
    //    EgammaRecHitIsolation ecalBarrelIsol03(egIsoConeSizeOutSmall,egIsoConeSizeInBarrel,egIsoJurassicWidth,egIsoPtMinBarrel,egIsoEMinBarrel,theCaloGeom,&ecalBarrelHits,DetId::Ecal);
    EgammaRecHitIsolation ecalBarrelIsol03(egIsoConeSizeOutSmall,egIsoConeSizeInBarrel,egIsoJurassicWidth,egIsoPtMinBarrel,egIsoEMinBarrel,theCaloGeom,&ecalBarrelHits,sevLevel.product(),DetId::Ecal);
    ecalBarrelIsol03.setUseNumCrystals(useNumCrystals);
    ecalBarrelIsol03.setVetoClustered(vetoClustered);
    //for42x
    //   // !!! Spike Removal... not in 361p4! Have to add it after !!!
    //    ecalBarrelIsol03.doSpikeRemoval(reducedEBRecHits.product(),pChannelStatus.product(),severityLevelCut,severityRecHitThreshold,spId,spikeIdThreshold);
    
    //std::cout << "FillSuperClusters  ugly " << std::endl;
    
    // ugly...
    reco::RecoEcalCandidate * cand = new RecoEcalCandidate();
    math::XYZPoint v(0,0,0); math::XYZVector p = isc->energy() * (isc->position() -v).unit(); double t = sqrt(0. + p.mag2());
    cand->setCharge(0); cand->setVertex(v); cand->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), t));		
    const reco::SuperClusterRef sc_ref(sc_coll_EB, index_sc);
    cand->setSuperCluster(sc_ref);
    //reco::SuperClusterRef sc = cand->get<reco::SuperClusterRef>();
    
    _sc_hybrid_ecalRecHitSumEt_dr03[index_sc] = ecalBarrelIsol03.getEtSum(cand);
    
    //	std::cout << "FillSuperClusters  Track Isolation " << std::endl;
    
    // Track Isolation
    // Calculate hollow cone track isolation, CONE 0.3
    reco::Photon * newPhoton = new Photon(); 
    newPhoton->setVertex(v); newPhoton->setCharge(0); newPhoton->setMass(0);
    newPhoton->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), isc->energy()));
    
    //int ntrk_03 = 0.; 
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
    
    _sc_hybrid_trkiso_dr03[index_sc] = trkiso_hc_03;
    
    //	std::cout << "FillSuperClusters SC EB -- modif-alex l1 matching " << std::endl;

    if(!aod_){
      // ____________________________
      // SC EB -- modif-alex l1 matching
      // LOOP MATCHING ON L1 trigger 
      // ____________________________
      nTow=0;
      nReg=0;
      for(int icc = 0; icc < 50; ++icc) {
	_sc_hybrid_TTetaVect[index_sc][icc] = -999;
	_sc_hybrid_TTphiVect[index_sc][icc] = -999;
	_sc_hybrid_TTetVect[index_sc][icc]  = 0.;
      }
      for(int icc = 0; icc < 10; ++icc) {
	_sc_hybrid_RCTetaVect[index_sc][icc]      = -999;
	_sc_hybrid_RCTphiVect[index_sc][icc]      = -999;
	_sc_hybrid_RCTetVect[index_sc][icc]       = 0.;
	_sc_hybrid_RCTL1isoVect[index_sc][icc]    = -999;
	_sc_hybrid_RCTL1nonisoVect[index_sc][icc] = -999;
      }
      
      for (reco::CaloCluster_iterator clus = isc->clustersBegin () ;
	   clus != isc->clustersEnd () ;
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
	    //if ( (detitr -> first).subdetId () == EcalBarrel)
	    //{
	    thishit = reducedRecHits->find ( (detitr -> first) ) ;
	    if (thishit == reducedRecHits->end ()) continue;
	    myhit = (*thishit) ;
	    EBDetId detid(thishit->id());
	    towid= detid.tower();
	    thetahit =  theBarrelGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	    //}//barrel rechit
	    //  else {
	    // 	    if ( (detitr -> first).subdetId () == EcalEndcap)
	    // 	      {
	    // 		thishit = reducedRecHits->find ( (detitr -> first) ) ;
	    // 		if (thishit == reducedRecHits->end ()) continue;
	    // 		myhit = (*thishit) ;
	    // 		EEDetId detid(thishit->id());
	    // 		towid= (*eTTmap_).towerOf(detid);
	    // 		thetahit =  theEndcapGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	    // 	      }
	    // 	    else continue;
	    // 	  }//endcap rechit
	    
	    int iETA=towid.ieta();
	    int iPHI=towid.iphi();
	    int iReta=getGCTRegionEta(iETA);
	    int iRphi=getGCTRegionPhi(iPHI);
	    double iET=myhit.energy()*sin(thetahit);
	    
	    bool newTow = true;
	    if(nTow>0) {
	      for (int iTow=0; iTow<nTow; ++iTow) {
		if(_sc_hybrid_TTetaVect[index_sc][iTow] == iETA && _sc_hybrid_TTphiVect[index_sc][iTow] == iPHI) {
		  newTow = false;
		  _sc_hybrid_TTetVect[index_sc][iTow] +=  iET;
		}
	      }
	    } // if nTow>0
	    if(newTow) {
	      _sc_hybrid_TTetaVect[index_sc][nTow] = iETA;
	      _sc_hybrid_TTphiVect[index_sc][nTow] = iPHI;
	      _sc_hybrid_TTetVect[index_sc][nTow] =  iET;
	      nTow++;
	    } // if newTow
	    
	    bool newReg = true;
	    if(nReg>0) {
	      for (int iReg=0; iReg<nReg; ++iReg) {
		if(_sc_hybrid_RCTetaVect[index_sc][iReg] == iReta && _sc_hybrid_RCTphiVect[index_sc][iReg] == iRphi) {
		  newReg = false;
		  _sc_hybrid_RCTetVect[index_sc][iReg] +=  iET;
		}
	      }
	    } // if newreg>0
	    
	    if(newReg) {
	      _sc_hybrid_RCTetaVect[index_sc][nReg] = iReta;
	      _sc_hybrid_RCTphiVect[index_sc][nReg] = iRphi;
	      _sc_hybrid_RCTetVect[index_sc][nReg] =  iET;
	      
	      for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
		if(_trig_L1emIso_iphi[il1] == iRphi && _trig_L1emIso_ieta[il1] == iReta) _sc_hybrid_RCTL1isoVect[index_sc][nReg] = _trig_L1emIso_rank[il1];
	      }
	      for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
		if(_trig_L1emNonIso_iphi[il1] == iRphi && _trig_L1emNonIso_ieta[il1] == iReta) _sc_hybrid_RCTL1nonisoVect[index_sc][nReg] = _trig_L1emNonIso_rank[il1];
	      }
	      nReg++;
	    } // if newReg
	  }//loop crystal
      }//loop cluster
      
      //double TTetmax=0.;
      //int iTTmax=-1;
      double TTetmax2 = 0.;
      int iTTmax2     = -1;
      
      for (int iTow=0; iTow<nTow; ++iTow) {
	bool nomaskTT = true;
	for (it=towerMap.begin();it!=towerMap.end();++it) {
	  if ((*it).second > 0) {
	    EcalTrigTowerDetId  ttId((*it).first);
	    if(ttId.ieta() == _sc_hybrid_TTetaVect[index_sc][iTow] && ttId.iphi() == _sc_hybrid_TTphiVect[index_sc][iTow]) {
	      nomaskTT=false;
	    } // if ttId ieta
	  } // if ut.second>0
	}//loop trigger towers
	
	if(nomaskTT && _sc_hybrid_TTetVect[index_sc][iTow] > TTetmax2) {
	  iTTmax2 = iTow;
	  TTetmax2 = _sc_hybrid_TTetVect[index_sc][iTow];
	} // if nomaskTT
      } // for loop on towers
      
      //int TTetamax = getGCTRegionEta(_sc_hybrid_TTetaVect[index_sc][iTTmax]);
      //int TTphimax = getGCTRegionPhi(_sc_hybrid_TTphiVect[index_sc][iTTmax]);
      //_sc_hybrid_RCTeta[index_sc]=TTetamax;
      //_sc_hybrid_RCTphi[index_sc]=TTphimax;
      
      //for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
      //if(_trig_L1emIso_iphi[il1] == TTphimax && _trig_L1emIso_ieta[il1] == TTetamax) _sc_hybrid_RCTL1iso[index_sc] = _trig_L1emIso_rank[il1];
      //}
      //for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
      //if(_trig_L1emNonIso_iphi[il1] == TTphimax && _trig_L1emNonIso_ieta[il1] == TTetamax) _sc_hybrid_RCTL1noniso[index_sc] = _trig_L1emNonIso_rank[il1];
      //}
      
      if(iTTmax2>=0) {
	int TTetamax2 = getGCTRegionEta(_sc_hybrid_TTetaVect[index_sc][iTTmax2]);
	int TTphimax2 = getGCTRegionPhi(_sc_hybrid_TTphiVect[index_sc][iTTmax2]);
	_sc_hybrid_RCTeta[index_sc] = TTetamax2;
	_sc_hybrid_RCTphi[index_sc] = TTphimax2;
	
	for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
	  if(_trig_L1emIso_iphi[il1] == TTphimax2 && _trig_L1emIso_ieta[il1] == TTetamax2) _sc_hybrid_RCTL1iso[index_sc] = _trig_L1emIso_rank[il1];
	}
	for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	  if(_trig_L1emNonIso_iphi[il1] == TTphimax2 && _trig_L1emNonIso_ieta[il1] == TTetamax2) _sc_hybrid_RCTL1noniso[index_sc] = _trig_L1emNonIso_rank[il1];
	}
      } // if iTTmax2
    }//!AOD	
    index_sc++;
  } // for loop on super clusters
  
  if(index_sc>24) { _sc_hybrid_N = 25; cout << "Number of SuperCluster > 25; _sc_hybrid_N set to 25" << endl;}
  
  //	std::cout << "FillSuperClusters  EE SuperCluster  " << std::endl;
  
  // -----------------
  //  EE SuperCluster
  // -----------------
  edm::Handle<reco::SuperClusterCollection> sc_coll_EE;
  iEvent.getByLabel(edm::InputTag("correctedMulti5x5SuperClustersWithPreshower"), sc_coll_EE);
  
  _sc_multi55_N = sc_coll_EE->size();
  
  //cout << " size EE = " << sc_coll_EE->size() << endl;
  
  int index_sc_EE = 0;
  
  for( reco::SuperClusterCollection::const_iterator isc=sc_coll_EE->begin(); isc!=sc_coll_EE->end(); isc++) { 
    if(index_sc_EE>24) continue;
    
    const EcalRecHitCollection * reducedRecHits = 0 ;
    reducedRecHits = reducedEERecHits.product() ; 
    //seed cluster analysis
    const reco::CaloCluster & seedCluster1 = *(isc->seed());
    
    // 4-vector
    double R  = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y() +isc->z()*isc->z());
    double Rt = TMath::Sqrt(isc->x()*isc->x() + isc->y()*isc->y());
    
    _sc_multi55_E[index_sc_EE]   = isc->energy();
    _sc_multi55_Et[index_sc_EE]  = isc->energy()*(Rt/R);
    _sc_multi55_Eta[index_sc_EE] = isc->eta();
    _sc_multi55_Phi[index_sc_EE] = isc->phi();
    
    // H/E
    reco::SuperCluster EmSCCand = *isc;
    double HoE = towerIso1_->getTowerESum(&EmSCCand) + towerIso2_->getTowerESum(&EmSCCand) ;
    HoE /= 	isc->energy() ;     
    _sc_multi55_he[index_sc_EE] = HoE;
    
    // SigmaIetaIeta
    std::vector<float> localCovariances = EcalClusterTools::localCovariances(seedCluster1,reducedRecHits,topology) ;
    _sc_multi55_sigmaietaieta[index_sc_EE]  = sqrt(localCovariances[0]) ;
		
    // Iso HCAL
    _sc_multi55_hcalDepth1TowerSumEt_dr03[index_sc_EE] = hadDepth1Isolation03.getTowerEtSum(&EmSCCand);
    _sc_multi55_hcalDepth2TowerSumEt_dr03[index_sc_EE] = hadDepth2Isolation03.getTowerEtSum(&EmSCCand);
    
    // Iso ECAL
    EcalRecHitMetaCollection ecalEndcapHits(*reducedEERecHits);
    //for42x
    //    EgammaRecHitIsolation ecalEndcapIsol03(egIsoConeSizeOutSmall,egIsoConeSizeInEndcap,egIsoJurassicWidth,egIsoPtMinEndcap,egIsoEMinEndcap,theCaloGeom,&ecalEndcapHits,DetId::Ecal);
    EgammaRecHitIsolation ecalEndcapIsol03(egIsoConeSizeOutSmall,egIsoConeSizeInEndcap,egIsoJurassicWidth,egIsoPtMinEndcap,egIsoEMinEndcap,theCaloGeom,&ecalEndcapHits,sevLevel.product(),DetId::Ecal);
    ecalEndcapIsol03.setUseNumCrystals(useNumCrystals);
    ecalEndcapIsol03.setVetoClustered(vetoClustered);
    // ugly...
    reco::RecoEcalCandidate * cand = new RecoEcalCandidate();
    math::XYZPoint v(0,0,0); math::XYZVector p = isc->energy() * (isc->position() -v).unit(); double t = sqrt(0. + p.mag2());
    cand->setCharge(0); cand->setVertex(v);
    cand->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), t));		
    const reco::SuperClusterRef sc_ref(sc_coll_EE, index_sc_EE);
    cand->setSuperCluster(sc_ref);
		
    _sc_multi55_ecalRecHitSumEt_dr03[index_sc_EE] = ecalEndcapIsol03.getEtSum(cand);
    
    // Track Isolation
    // Calculate hollow cone track isolation, CONE 0.3
    reco::Photon * newPhoton = new Photon(); 
    newPhoton->setVertex(v); newPhoton->setCharge(0); newPhoton->setMass(0);
    newPhoton->setP4(reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), isc->energy()));
    
    //int ntrk_03 = 0.; 
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
    
    _sc_multi55_trkiso_dr03[index_sc_EE] = trkiso_hc_03;
    
    //std::cout << "FillSuperClusters  SC EE -- modif-alex l1 matching  " << std::endl;
    
    
    // ____________________________
    // SC EE -- modif-alex l1 matching
    // LOOP MATCHING ON L1 trigger 
    // ____________________________
    if(!aod_){
      nTow = 0;
      nReg = 0;
      for(int icc = 0; icc < 50; ++icc) {
	_sc_multi55_TTetaVect[index_sc_EE][icc] = -999;
	_sc_multi55_TTphiVect[index_sc_EE][icc] = -999;
	_sc_multi55_TTetVect[index_sc_EE][icc]  = 0.;
      }
      for(int icc = 0; icc < 10; ++icc) {
	_sc_multi55_RCTetaVect[index_sc_EE][icc]      = -999;
	_sc_multi55_RCTphiVect[index_sc_EE][icc]      = -999;
	_sc_multi55_RCTetVect[index_sc_EE][icc]       = 0.;
	_sc_multi55_RCTL1isoVect[index_sc_EE][icc]    = -999;
	_sc_multi55_RCTL1nonisoVect[index_sc_EE][icc] = -999;
      }
      
      for (reco::CaloCluster_iterator clus = isc->clustersBegin () ;
	   clus != isc->clustersEnd () ;
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
	    // 	    if ( (detitr -> first).subdetId () == EcalEndcap)
	    // 	      {
	    thishit = reducedRecHits->find ( (detitr -> first) ) ;
	    if (thishit == reducedRecHits->end ()) continue;
	    myhit = (*thishit) ;
	    EEDetId detid(thishit->id());
	    towid= (*eTTmap_).towerOf(detid);
	    thetahit =  theEndcapGeometry_->getGeometry((detitr -> first))->getPosition().theta();
	    // 	      }
	    // 	    else continue;
	    // 	  }//endcap rechit
	    
	    int iETA=towid.ieta();
	    int iPHI=towid.iphi();
	    int iReta=getGCTRegionEta(iETA);
	    int iRphi=getGCTRegionPhi(iPHI);
	    double iET=myhit.energy()*sin(thetahit);
	    
	    bool newTow = true;
	    if(nTow>0) {
	      for (int iTow=0; iTow<nTow; ++iTow) {
		if(_sc_multi55_TTetaVect[index_sc_EE][iTow] == iETA && _sc_multi55_TTphiVect[index_sc_EE][iTow] == iPHI) {
		  newTow = false;
		  _sc_multi55_TTetVect[index_sc_EE][iTow] +=  iET;
		}
	      }
	    } // if nTow>0
	    if(newTow) {
	      _sc_multi55_TTetaVect[index_sc_EE][nTow] = iETA;
	      _sc_multi55_TTphiVect[index_sc_EE][nTow] = iPHI;
	      _sc_multi55_TTetVect[index_sc_EE][nTow] =  iET;
	      nTow++;
	    } // if newTow
	    
	    bool newReg = true;
	    if(nReg>0) {
	      for (int iReg=0; iReg<nReg; ++iReg) {
		if(_sc_multi55_RCTetaVect[index_sc_EE][iReg] == iReta && _sc_multi55_RCTphiVect[index_sc_EE][iReg] == iRphi) {
		  newReg = false;
		  _sc_multi55_RCTetVect[index_sc_EE][iReg] +=  iET;
		}
	      }
	    } // if newreg>0
	    
	    if(newReg) {
	      _sc_multi55_RCTetaVect[index_sc_EE][nReg] = iReta;
	      _sc_multi55_RCTphiVect[index_sc_EE][nReg] = iRphi;
	      _sc_multi55_RCTetVect[index_sc_EE][nReg] =  iET;
	      
	      for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
		if(_trig_L1emIso_iphi[il1] == iRphi && _trig_L1emIso_ieta[il1] == iReta) _sc_multi55_RCTL1isoVect[index_sc_EE][nReg] = _trig_L1emIso_rank[il1];
	      }
	      for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
		if(_trig_L1emNonIso_iphi[il1] == iRphi && _trig_L1emNonIso_ieta[il1] == iReta) _sc_multi55_RCTL1nonisoVect[index_sc_EE][nReg] = _trig_L1emNonIso_rank[il1];
	      }
	      nReg++;
	    } // if newReg
	  }//loop crystal
      }//loop cluster
      
      //double TTetmax=0.;
      //int iTTmax=-1;
      double TTetmax2 = 0.;
      int iTTmax2     = -1;
      
      for (int iTow=0; iTow<nTow; ++iTow) {
	bool nomaskTT = true;
	for (it=towerMap.begin();it!=towerMap.end();++it) {
	  if ((*it).second > 0) {
	    EcalTrigTowerDetId  ttId((*it).first);
	    if(ttId.ieta() == _sc_multi55_TTetaVect[index_sc_EE][iTow] && ttId.iphi() == _sc_multi55_TTphiVect[index_sc_EE][iTow]) {
	      nomaskTT=false;
	    } // if ttId ieta
	  } // if ut.second>0
	}//loop trigger towers
	
	if(nomaskTT && _sc_multi55_TTetVect[index_sc_EE][iTow] > TTetmax2) {
	  iTTmax2 = iTow;
	  TTetmax2 = _sc_multi55_TTetVect[index_sc_EE][iTow];
	} // if nomaskTT
      } // for loop on towers
      
      if(iTTmax2>=0) {
	int TTetamax2 = getGCTRegionEta(_sc_multi55_TTetaVect[index_sc_EE][iTTmax2]);
	int TTphimax2 = getGCTRegionPhi(_sc_multi55_TTphiVect[index_sc_EE][iTTmax2]);
	_sc_multi55_RCTeta[index_sc_EE] = TTetamax2;
	_sc_multi55_RCTphi[index_sc_EE] = TTphimax2;
	
	for(int il1=0; il1<_trig_L1emIso_N; ++il1) {
	  if(_trig_L1emIso_iphi[il1] == TTphimax2 && _trig_L1emIso_ieta[il1] == TTetamax2) _sc_multi55_RCTL1iso[index_sc_EE] = _trig_L1emIso_rank[il1];
	}
	for(int il1=0; il1<_trig_L1emNonIso_N; ++il1) {
	  if(_trig_L1emNonIso_iphi[il1] == TTphimax2 && _trig_L1emNonIso_ieta[il1] == TTetamax2) _sc_multi55_RCTL1noniso[index_sc_EE] = _trig_L1emNonIso_rank[il1];
	}
      } // if iTTmax2
    }
    
    index_sc_EE++;
  } // for loop on super clusters
  
  if(index_sc_EE>24) { _sc_multi55_N = 25; cout << "Number of SuperCluster > 25; _sc_multi55_N set to 25" << endl;}
	
} // FillSuperCluster 

// ====================================================================================
void SimpleNtpleAOD::FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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


// ====================================================================================
void  SimpleNtpleAOD::FillTipLipIp(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// ====================================================================================
{
  //Get the B-field
  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  
  //Get Beam Spot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
  iEvent.getByType(recoBeamSpotHandle) ;
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
  //const MultiTrajectoryStateTransform *theMtsTransform = new MultiTrajectoryStateTransform;
  const GsfPropagatorAdapter *theGeomPropBw = new GsfPropagatorAdapter(AnalyticalPropagator(magneticField.product(),oppositeToMomentum)); 
  // muons:
  //const TrajectoryStateTransform *theTransform = new TrajectoryStateTransform;
  Propagator *thePropBw = new AnalyticalPropagator(magneticField.product(),oppositeToMomentum);                    
  // 
  //edm::ESHandle<TrackerGeometry> trackerHandle;
  //iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle); 
  //
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
      //innerMuTSOS = mutranstrack.stateOnSurface(BSVertex);
      innerMuTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), BSVertex, mutranstrack.field());
    } 
    else {
      //innerMuTSOS = mutranstrack.stateOnSurface(pVertex);
      innerMuTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), pVertex, mutranstrack.field());
    } 

    // get initial TSOS (now protected against STA muons):
    //if (!mutrack.isNull())
    //{
    //innerMuTSOS = theTransform->innerStateOnSurface(*mutrack, *trackerHandle.product(), magneticField.product());
    //}	
    
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
    
    //if (mutrack.isNull()){
    //mutrack=muCand->get<TrackRef,reco::StandAloneMuonTag>();
    //}  
    //mutranstrack = trackBuilder->build( mutrack ) ;

    // -------------
    //  3DIP & SIP 
    // -------------
    
    TrajectoryStateOnSurface muTSOS;
    if (useBeamSpot_==true){ 
      //muTSOS = mutranstrack.stateOnSurface(BSVertex);
      muTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), BSVertex, mutranstrack.field());
    } 
    else {
      //muTSOS = mutranstrack.stateOnSurface(pVertex);
      muTSOS = IPTools::transverseExtrapolate(mutranstrack.impactPointState(), pVertex, mutranstrack.field());
    }
    
    if (muTSOS.isValid()){
      std::pair<bool,Measurement1D> muIPpair;
      
      if (useBeamSpot_==true){ 		  
	//muIPpair = IPTools::signedImpactParameter3D(mutranstrack, muTSOS.globalDirection(), BSprimVertex);
	muIPpair = IPTools::absoluteImpactParameter3D(mutranstrack, BSprimVertex);
      } 
      else {	 
	//muIPpair = IPTools::signedImpactParameter3D(mutranstrack, muTSOS.globalDirection(), primVertex);
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
	//innerTSOS = eletranstrack.stateOnSurface(BSVertex);
	innerTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), BSVertex, eletranstrack.field());
      } 
      else {
	//innerTSOS = eletranstrack.stateOnSurface(pVertex);
	innerTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
      }

      //= theMtsTransform->innerStateOnSurface(*eletrack, *trackerHandle.product(), magneticField.product());
      //
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
      //eletranstrack = trackBuilder->build( eletrack ) ;
      
      TrajectoryStateOnSurface eleTSOS;
      if (useBeamSpot_==true){ 
	//eleTSOS = eletranstrack.stateOnSurface(BSVertex);
	eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), BSVertex, eletranstrack.field());
      } 
      else {
	//eleTSOS = eletranstrack.stateOnSurface(pVertex);
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
    
	
}

// ====================================================================================
void SimpleNtpleAOD::Init()
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
	
  // Skim
  _skim_is1lepton  = 0;
  _skim_is2leptons = 0;
  _skim_is3leptons = 0;
	
  // Vertices
  _vtx_N = 0; 
  for(int iv=0;iv<200;iv++) {
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
	
  // Spikes
  for(int isp=0 ; isp<5000 ; isp++) {
    spike_outOfTime[isp] = -999;
    spike_severityLevel[isp] = -999 ;
    spike_Et[isp] = -999 ;
    spike_phi[isp] = -999 ;
    spike_eta[isp] = -999 ;
    spike_theta[isp] = -999 ;
    spike_TTiphi[isp] = -999 ;
    spike_TTieta[isp] = -999 ;
    spike_Riphi[isp] = -999 ;
    spike_Rieta[isp] = -999 ;
  }


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
  for (int i=0;i<250;i++) 
    trig_hltInfo[i] = 0;

  for(int iTrig=0 ; iTrig<4 ; iTrig++)
    trig_HLT_path[iTrig] = 0 ;
	
  trig_isPhoton10 = 0; 
  trig_isPhoton15 = 0; 
  trig_isEle10_LW = 0;
  trig_isEle15_LW = 0;
	
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
		
		
      ele_eidVeryLoose[i] = 0.; 
      ele_eidLoose[i] = 0.;
      ele_eidMedium[i] = 0.;
      ele_eidTight[i] = 0.; 
      ele_eidSuperTight[i] = 0.; 
      ele_eidHyperTight1[i] = 0.;
      ele_eidHyperTight2[i] = 0.;
      ele_eidHyperTight3[i] = 0.;
      ele_eidHyperTight4[i] = 0.;
		
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
		
      //tkIso03Rel[i]=-999.;
      //		ecalIso03Rel[i]=-999.;
      //hcalIso03Rel[i]=-999.;
		
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
  _sc_hybrid_N = 0; 
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
		
  } // for loop on EE superclusters
	
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
void SimpleNtpleAOD::beginJob(const edm::ParameterSet& conf)
// ====================================================================================
{
  //hcalhelper_ = new ElectronHcalHelper(conf);
  //edm::Ref<reco::GsfElectronCollection> electronEdmRef(EleHandle,i);
  
  //	reco::SuperCluster EmSCCand1; // = *isc;
  //reco::RecoEcalCandidate EcalCand;
  
  //sc_struct = new converter::SuperClusterToCandidate(conf);
	
}

// ====================================================================================
void SimpleNtpleAOD::endJob() {}
// ====================================================================================

// ====================================================================================
void SimpleNtpleAOD::setMomentum (TLorentzVector &myvector, const LorentzVector & mom)
// ====================================================================================
{
	myvector.SetPx (mom.Px());
	myvector.SetPy (mom.Py());
	myvector.SetPz (mom.Pz());
	myvector.SetE (mom.E());
}

// ====================================================================================
bool SimpleNtpleAOD::IsConv (const reco::GsfElectron & eleRef) //edm::Ref<reco::GsfElectronCollection> eleRef)
// ====================================================================================
{
	
	bool isAmbiguous = true, isNotFromPixel = true;
	if (eleRef.ambiguousGsfTracksSize() == 0 ) {isAmbiguous = false;}
	
	/*
	 TrackingRecHitRef rhit =  eleRef->gsfTrack()->extra()->recHit(0);
	 int subdetId = rhit->geographicalId().subdetId();
	 int layerId  = 0;
	 DetId id = rhit->geographicalId();
	 if (id.subdetId()==3) layerId = ((TIBDetId)(id)).layer();
	 if (id.subdetId()==5) layerId = ((TOBDetId)(id)).layer();
	 if (id.subdetId()==1) layerId = ((PXBDetId)(id)).layer();
	 if (id.subdetId()==4) layerId = ((TIDDetId)(id)).wheel();
	 if (id.subdetId()==6) layerId = ((TECDetId)(id)).wheel();
	 if (id.subdetId()==2) layerId = ((PXFDetId)(id)).disk();
	 //std::cout << " subdetIdele layerIdele = " << id.subdetId() << "   " << layerId << std::endl;
	 
	 if ((id.subdetId()==1 && layerId == 1) || (id.subdetId()==2 && layerId == 1)) {isNotFromPixel = false;}
	 */
	
	int  mishits = eleRef.gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 
	
	//std::cout << "mishits = " << mishits << std::endl;
	if (mishits == 0){isNotFromPixel = false;}
	
	// 
	bool is_conversion = false;
	
	if(isAmbiguous || isNotFromPixel) is_conversion = true;
	
	return is_conversion;
	
}

// ====================================================================================
const EcalRecHit SimpleNtpleAOD::getRecHit(DetId id, const EcalRecHitCollection *recHits)
// ====================================================================================
{
	if ( id == DetId(0) ) {
		return EcalRecHit();
	} else {
		EcalRecHitCollection::const_iterator it = recHits->find( id );
		if ( it != recHits->end() ) {
			return (*it);
		} else {
			//throw cms::Exception("EcalRecHitNotFound") << "The recHit corresponding to the DetId" << id.rawId() << " not found in the EcalRecHitCollection";
			// the recHit is not in the collection (hopefully zero suppressed)
			return EcalRecHit();
		}
	}
	return EcalRecHit();
}



//modif-alex
//GETTING RCT regions
// ====================================================================================
int SimpleNtpleAOD::getGCTRegionPhi(int ttphi)
// ====================================================================================
{
	int gctphi=0;
	gctphi = (ttphi+1)/4;
	if(ttphi<=2) gctphi=0;
	if(ttphi>=71) gctphi=0;
	
	return gctphi;
}

// ====================================================================================
int SimpleNtpleAOD::getGCTRegionEta(int tteta)
// ====================================================================================
{
	int gcteta = 0;
	
	if(tteta>0) gcteta = (tteta-1)/4 + 11;
	else if(tteta<0) gcteta = (tteta+1)/4 + 10;
	
	return gcteta;
}
/*
// ===============================================================================================
// unified acces to isolations
std::pair<int,double> ElectronTkIsolation::getIso(const reco::GsfElectron* electron) const  
// ===============================================================================================
{
	int counter  =0 ;
	double ptSum =0.;
	//Take the electron track
	reco::GsfTrackRef tmpTrack = electron->gsfTrack() ;
	math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).momentum () ; 
	double tmpElectronEtaAtVertex = (*tmpTrack).eta();
	
	
	for ( reco::TrackCollection::const_iterator itrTr  = (*trackCollection_).begin() ; 
		 itrTr != (*trackCollection_).end()   ; 
		 ++itrTr ) {
		
		math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).momentum () ; 
		
		double this_pt  = (*itrTr).pt();
		if ( this_pt < ptLow_ ) continue;
		
		double dzCut = 0;
		switch( dzOption_ ) {
			case egammaisolation::EgammaTrackSelector::dz : dzCut = fabs( (*itrTr).dz() - (*tmpTrack).dz() ); break;
			case egammaisolation::EgammaTrackSelector::vz : dzCut = fabs( (*itrTr).vz() - (*tmpTrack).vz() ); break;
			case egammaisolation::EgammaTrackSelector::bs : dzCut = fabs( (*itrTr).dz(beamPoint_) - (*tmpTrack).dz(beamPoint_) ); break;
			case egammaisolation::EgammaTrackSelector::vtx: dzCut = fabs( (*itrTr).dz(tmpTrack->vertex()) ); break;
			default : dzCut = fabs( (*itrTr).vz() - (*tmpTrack).vz() ); break;
		}
		if (dzCut > lip_ ) continue;
		if (fabs( (*itrTr).dxy(beamPoint_) ) > drb_   ) continue;
		double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(),tmpElectronMomentumAtVtx) ;
		double deta = (*itrTr).eta() - tmpElectronEtaAtVertex;
		if (fabs(tmpElectronEtaAtVertex) < 1.479) { 
			if ( fabs(dr) < extRadius_ && fabs(dr) >= intRadiusBarrel_ && fabs(deta) >= stripBarrel_)
			{
				++counter ;
				ptSum += this_pt;
			}
		}
		else {
			if ( fabs(dr) < extRadius_ && fabs(dr) >= intRadiusEndcap_ && fabs(deta) >= stripEndcap_)
			{
				++counter ;
				ptSum += this_pt;
			}
		}
		
	}//end loop over tracks                 
	
	std::pair<int,double> retval;
	retval.first  = counter;
	retval.second = ptSum;
	
	return retval;
} // end of get TrkIso
*/



// ====================================================================================
float SimpleNtpleAOD::E2overE9( const DetId id, const EcalRecHitCollection & recHits, 
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
float SimpleNtpleAOD::recHitE( const DetId id, const EcalRecHitCollection &recHits )
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
float SimpleNtpleAOD::recHitE( const DetId id, const EcalRecHitCollection & recHits,
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
float SimpleNtpleAOD::recHitApproxEt( const DetId id, const EcalRecHitCollection &recHits )
// ====================================================================================
{
        // for the time being works only for the barrel
        if ( id.subdetId() == EcalBarrel ) {
                return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
        }
        return 0;
}
