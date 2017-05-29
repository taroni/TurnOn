// C++
#include <memory>
#include <iostream>

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TVector3.h"

// CMSSW
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

// Pile UP
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
// Vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
// L1 Trigger
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "L1Trigger/L1ExtraFromDigis/interface/L1ExtraParticleMapProd.h"
// RCT
#include "CondFormats/L1TObjects/interface/L1RCTChannelMask.h"
#include "CondFormats/DataRecord/interface/L1RCTChannelMaskRcd.h"
// TPG
#include "CondFormats/DataRecord/interface/EcalTPGTowerStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTPGTowerStatus.h"
#include "CondFormats/DataRecord/interface/EcalTPGCrystalStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTPGCrystalStatus.h"
// TPG (Nadir study)
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
// Electron/SuperCluster
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
// PF electron
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
// For Photon Iso
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "RecoEgamma/PhotonIdentification/interface/PhotonIsolationCalculator.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/PhotonTkIsolation.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
// Calo Jets
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
// JPT Jets
#include "DataFormats/JetReco/interface/JPTJet.h"
// PF Jets
#include "DataFormats/JetReco/interface/PFJet.h"
// MET
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

// TrackingParticles
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
//#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include <vector>

//Clusters 
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h" 
// For H/E - Iso on SC
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
// For Skim
#include "EGamma/ECGelec/interface/AnalysisUtils.h"
//
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "TLorentzVector.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
// Transient tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/GsfTools/interface/GSUtilities.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
//
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

// Other specific
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


/// among the includes ///
class MultiTrajectoryStateMode ;
// for H/E
class EgammaTowerIsolation ;

//
// class declaration
//

class SimpleNtpleAOD : public edm::EDAnalyzer {
 public:
  explicit SimpleNtpleAOD(const edm::ParameterSet&);
  ~SimpleNtpleAOD();
	
  typedef math::XYZTLorentzVector LorentzVector ;
  typedef edm::View<reco::Track> trackCollection ;
	
 private:
  virtual void beginJob(const edm::ParameterSet& conf) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
	
  void Init();
	
  void FillEvent (const edm::Event&, const edm::EventSetup&);
  void FillTrigger (const edm::Event&, const edm::EventSetup&);
  void FillEle (const edm::Event&, const edm::EventSetup&);
  void FillMuons (const edm::Event&, const edm::EventSetup&);
  void FillMET (const edm::Event&, const edm::EventSetup&);
  void FillJets(const edm::Event&, const edm::EventSetup&);
  void FillSuperClusters(const edm::Event&, const edm::EventSetup&);
  void FillTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void FillTipLipIp(const edm::Event& iEvent, const edm::EventSetup& iSetup);
		
  void setMomentum (TLorentzVector &myvector, const LorentzVector & mom) ;
  bool IsConv (const reco::GsfElectron & eleRef);
  // For Trigger
  int getGCTRegionPhi(int ttphi) ; //modif-alex
  int getGCTRegionEta(int tteta) ; //modif-alex

  // e2overe9 new anti-spike variable steph
  float E2overE9( const DetId , const EcalRecHitCollection & , float , float ,  bool , bool);
  float recHitE( const DetId , const EcalRecHitCollection & );
  float recHitE( const DetId , const EcalRecHitCollection & ,  int , int );
  float recHitApproxEt( const DetId , const EcalRecHitCollection &);
  // end e2overe9

  // ----------member data ---------------------------
  ElectronHcalHelper * hcalhelper_;
	
  TTree *mytree_;
	
  int nEvent, nRun, nLumi;
	
  // Vertices
  int _vtx_N;
  double _vtx_x[200], _vtx_y[200], _vtx_z[200];
  double _vtx_normalizedChi2[200], _vtx_ndof[200], _vtx_nTracks[200], _vtx_d0[200];

  //Pile-up
  int _PU_N;
  double _PU_rho, _PU_sigma;  //corrections from FastJets

  // Skim
  int _skim_is1lepton;
  int _skim_is2leptons;
  int _skim_is3leptons;

  // original TP
  int _trig_tower_N, _trig_tower_ieta[4032],_trig_tower_iphi[4032],_trig_tower_adc[4032], _trig_tower_sFGVB[4032]; 

  // cleaned TP
  int _trig_tower_N_modif, _trig_tower_ieta_modif[4032],_trig_tower_iphi_modif[4032],
    _trig_tower_adc_modif[4032], _trig_tower_sFGVB_modif[4032]; 

  // emulated TP
  int _trig_tower_N_emul, _trig_tower_ieta_emul[4032],_trig_tower_iphi_emul[4032],
    _trig_tower_adc_emul[4032][5], _trig_tower_sFGVB_emul[4032][5]; 

  // Trigger Paths
  int trig_hltInfo[250];
  int trig_isPhoton10, trig_isPhoton15, trig_isEle10_LW, trig_isEle15_LW;
  int _trig_isEleHLTpath, _trig_isMuonHLTpath;
  int trig_HLT_path[4]; // unbias, EG5, EG8, EG12
  char trig_fired_names[5000];
  // L1
  int _trig_L1emIso_N; 
  int _trig_L1emNonIso_N;
  int _trig_L1emIso_ieta[4], _trig_L1emIso_iphi[4], _trig_L1emIso_rank[4]; 
  double _trig_L1emIso_eta[4], _trig_L1emIso_phi[4],_trig_L1emIso_energy[4],_trig_L1emIso_et[4]; 
  int _trig_L1emNonIso_ieta[4], _trig_L1emNonIso_iphi[4],_trig_L1emNonIso_rank[4];
  double _trig_L1emNonIso_eta[4], _trig_L1emNonIso_phi[4], _trig_L1emNonIso_energy[4],_trig_L1emNonIso_et[4];

  // L1 modif
  int _trig_L1emIso_N_M; 
  int _trig_L1emNonIso_N_M;
  int _trig_L1emIso_ieta_M[4], _trig_L1emIso_iphi_M[4], _trig_L1emIso_rank_M[4]; 
  double _trig_L1emIso_eta_M[4], _trig_L1emIso_phi_M[4],_trig_L1emIso_energy_M[4],_trig_L1emIso_et_M[4]; 
  int _trig_L1emNonIso_ieta_M[4], _trig_L1emNonIso_iphi_M[4],_trig_L1emNonIso_rank_M[4];
  double _trig_L1emNonIso_eta_M[4], _trig_L1emNonIso_phi_M[4], _trig_L1emNonIso_energy_M[4],_trig_L1emNonIso_et_M[4];

  // L1 prefiring
  int _trig_preL1emIso_N; 
  int _trig_preL1emNonIso_N;
  int _trig_preL1emIso_ieta[4], _trig_preL1emIso_iphi[4], _trig_preL1emIso_rank[4]; 
  int _trig_preL1emNonIso_ieta[4], _trig_preL1emNonIso_iphi[4],_trig_preL1emNonIso_rank[4];
  // L1 postfiring
  int _trig_postL1emIso_N; 
  int _trig_postL1emNonIso_N;
  int _trig_postL1emIso_ieta[4], _trig_postL1emIso_iphi[4], _trig_postL1emIso_rank[4]; 
  int _trig_postL1emNonIso_ieta[4], _trig_postL1emNonIso_iphi[4],_trig_postL1emNonIso_rank[4];
	
  int _trig_nMaskedRCT, _trig_nMaskedCh;
  int _trig_iMaskedRCTeta[100], _trig_iMaskedRCTphi[100], _trig_iMaskedRCTcrate[100], _trig_iMaskedTTeta[100], _trig_iMaskedTTphi[100];

  // HLT
  int _trig_HLT_N;
  double _trig_HLT_eta[20], _trig_HLT_phi[20], _trig_HLT_energy[20], _trig_HLT_pt[20];
  int _trig_HLT_name[20];
  //std::vector <std::string> _trig_HLT_name;
	
  // Beam Spot
  double BS_x, BS_y, BS_z, BS_dydz, BS_dxdz, BS_dz, BS_bw_x, BS_bw_y;
  GlobalPoint vertexPosition;

  // Gen Particles
  TClonesArray * _m_MC_gen_V;
  TClonesArray * _m_MC_gen_leptons;
  double _MC_gen_V_pdgid[10];
  double _MC_gen_leptons_pdgid[10];

  // MC Properties
  double _MC_pthat;
  int _MC_flavor[2];

  // Spikes info
  int spike_N, spike_SwissCross[5000], spike_TTieta[5000], spike_TTiphi[5000], spike_Rieta[5000], spike_Riphi[5000], spike_severityLevel[5000], spike_outOfTime[5000];
  double  spike_Et[5000], spike_eta[5000], spike_phi[5000], spike_theta[5000], spike_time[5000];

  // Electrons
  //electrons;
  int ele_N, ele_nSeed;
  int ele_echarge[10];
  double ele_he[10] , 
    ele_eseedpout[10] , ele_ep[10] , ele_eseedp[10] , ele_eelepout[10] ,       
    ele_deltaetaseed[10] , ele_deltaetaele[10] , ele_deltaphiseed[10] , ele_deltaphiele[10] , ele_deltaetain[10] , ele_deltaphiin[10] ,
    ele_sigmaietaieta[10] , ele_sigmaetaeta[10] , ele_e15[10] , ele_e25max[10] , ele_e55[10] , ele_e1[10] , ele_e33[10] , ele_e2overe9[10],
    ele_pin_mode[10] , ele_pout_mode[10] , ele_pin_mean[10] , ele_pout_mean[10] , 
    ele_calo_energy[10] ,
    ele_pTin_mode[10] , ele_pTout_mode[10] , ele_pTin_mean[10] , ele_pTout_mean[10] ;
  double ele_fbrem[10], ele_mva[10] ;
  double ele_sclE[10], ele_sclEt[10], ele_sclEta[10], ele_sclPhi[10] ;
  double ele_sclRawE[10], ele_sclErr[10], ele_sclErr_pos[10], ele_sclErr_neg[10], ele_trErr[10], ele_momErr[10], ele_newmom[10], ele_newmomErr[10] ;
  double ele_sclEpresh[10] ;
  double ele_sclX[10], ele_sclY[10], ele_sclZ[10];
  double ele_tr_atcaloX[10],ele_tr_atcaloY[10],ele_tr_atcaloZ[10];
  double ele_firsthit_X[10],ele_firsthit_Y[10],ele_firsthit_Z[10];
  int ele_isbarrel[10] , ele_isendcap[10] , ele_isEBetaGap[10] , 
    ele_isEBphiGap[10] , ele_isEEdeeGap[10] , ele_isEEringGap[10] ,
    ele_eClass[10], ele_isecalDriven[10] , ele_istrackerDriven[10] ;
  double ele_vertex_x[10], ele_vertex_y[10], ele_vertex_z[10];
  int ele_missing_hits[10], ele_lost_hits[10]; 
  double ele_chi2_hits[10], ele_dxyB[10], ele_dxy[10], ele_dzB[10], ele_dz[10], ele_dszB[10], ele_dsz[10];              
  double ele_dzPV[10], ele_dzPV_error[10], ele_dxyPV[10], ele_dxyPV_error[10], ele_dszPV[10], ele_dszPV_error[10];
  double ele_track_x[10], ele_track_y[10], ele_track_z[10];
  double ele_tkSumPt_dr03[10] , ele_ecalRecHitSumEt_dr03[10] , ele_hcalDepth1TowerSumEt_dr03[10] , ele_hcalDepth2TowerSumEt_dr03[10] ,
    ele_tkSumPt_dr04[10] , ele_ecalRecHitSumEt_dr04[10] , ele_hcalDepth1TowerSumEt_dr04[10] , ele_hcalDepth2TowerSumEt_dr04[10] ;
  double ele_hcalDepth1plus2TowerSumEt_00615dr03[10], ele_hcalDepth1plus2TowerSumEt_005dr03[10], ele_hcalDepth1plus2TowerSumEt_0dr03[10];
  int ele_ambiguousGsfTracks[10] ;
  double ele_ambiguousGsfTracksdxy[10][5], ele_ambiguousGsfTracksdz[10][5], 
    ele_ambiguousGsfTracksdxyB[10][5], ele_ambiguousGsfTracksdzB[10][5] ;
  int ele_seedSubdet1[10];
  double ele_seedDphi1Pos[10], ele_seedDrz1Pos[10], ele_seedDphi1Neg[10], ele_seedDrz1Neg[10];
  int ele_seedSubdet2[10];
  double ele_seedDphi2Pos[10], ele_seedDrz2Pos[10], ele_seedDphi2Neg[10], ele_seedDrz2Neg[10];
  //eID CiC
  double ele_eidVeryLoose[10], ele_eidLoose[10], ele_eidMedium[10], 
    ele_eidTight[10], ele_eidSuperTight[10], 
    ele_eidHyperTight1[10], ele_eidHyperTight2[10], ele_eidHyperTight3[10], ele_eidHyperTight4[10] ;
  double ele_tkSumPtTdrHzz_dr025[10], ele_tkSumPtoPtTdrHzz_dr025[10], ele_hcalSumEtTdrHzz_dr02[10], ele_hcalSumEtoPtTdrHzz_dr02[10]  ;
  double ele_tkSumPtEg4Hzz_dr03[10], ele_ecalSumEtEg4Hzz_dr03[10], ele_hcalSumEtEg4Hzz_dr04[10];
  double ele_tkSumPtoPtEg4Hzz_dr03[10], ele_ecalSumEtoPtEg4Hzz_dr03[10], ele_hcalSumEtoPtEg4Hzz_dr04[10];
	
  int ele_severityLevelSeed[10], ele_severityLevelClusters[10], ele_outOfTimeSeed[10], ele_outOfTimeClusters[10];

  int ele_isMCEle[10], ele_isMCPhoton[10], ele_isMCHadron[10], ele_isSIM[10], ele_isSIMEle[10];
  int ele_idPDGMatch[10], ele_idPDGmother_MCEle[10], ele_idPDGMatchSim[10];
	
  double ele_ECAL_fbrem[10], ele_PFcomb[10], ele_PFcomb_Err[10], ele_PF_SCenergy[10], ele_PF_SCenergy_Err[10];
	

  //Seeds general collection
  int ele_SeedSubdet1[100];
  double ele_SeedDphi1Pos[100], ele_SeedDrz1Pos[100], ele_SeedDphi1Neg[100], ele_SeedDrz1Neg[100];

  int ele_SeedSubdet2[100];
  double ele_SeedDphi2Pos[100], ele_SeedDrz2Pos[100], ele_SeedDphi2Neg[100], ele_SeedDrz2Neg[100];
    
  int ele_SeedIsEcalDriven[100], ele_SeedIsTrackerDriven[100];	


  // MC matching
  double ele_MC_chosenEle_PoP_px[10], ele_MC_chosenEle_PoP_py[10], ele_MC_chosenEle_PoP_pz[10], ele_MC_chosenEle_PoP_e[10];
  double ele_MC_chosenPho_PoP_px[10], ele_MC_chosenPho_PoP_py[10], ele_MC_chosenPho_PoP_pz[10], ele_MC_chosenPho_PoP_e[10];
  double ele_MC_chosenHad_PoP_px[10], ele_MC_chosenHad_PoP_py[10], ele_MC_chosenHad_PoP_pz[10], ele_MC_chosenHad_PoP_e[10];
  double ele_MC_closest_DR_px[10], ele_MC_closest_DR_py[10], ele_MC_closest_DR_pz[10], ele_MC_closest_DR_e[10];

  // Conversion removal
  int ele_isConversion[10];
  int ele_convFound[10];
  double ele_conv_dcot[10];
  double ele_conv_dist[10];

  // For Charge (Clemy's stuff)
  int ele_chargeGsfSC[10], ele_chargeGsfCtf[10], ele_chargeGsfCtfSC[10], ele_CtfTrackExists[10];
  double ele_chargeDPhiInnEle[10], ele_chargeDPhiInnEleCorr[10], ele_chargeQoverPGsfVtx[10], ele_chargeQoverPCtf[10];
  //double tkIso03Rel[10], ecalIso03Rel[10], hcalIso03Rel[10] ;
  int  ele_expected_inner_hits[10], ele_sclNclus[10];

  const MultiTrajectoryStateTransform *mtsTransform_;
	
  edm::ESHandle<MagneticField> theMagField;
  edm::ESHandle<TrackerGeometry> trackerHandle_;
	
  unsigned long long cacheIDTDGeom_;
  unsigned long long cacheIDMagField_;

  // For L1 Trigger (Clemy's stuff)
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
		
  const CaloSubdetectorGeometry * theEndcapGeometry_ ;
  const CaloSubdetectorGeometry * theBarrelGeometry_ ;
  edm::ESHandle<EcalTrigTowerConstituentsMap> eTTmap_;

  // Vector for electrons
  TClonesArray * m_electrons ;
	
  TLorentzVector myvector ;    

  // MET
  double _met_calo_et,_met_calo_px, _met_calo_py, _met_calo_phi, _met_calo_set, _met_calo_sig; 
  double _met_calomu_et,_met_calomu_px, _met_calomu_py, _met_calomu_phi, _met_calomu_set, _met_calomu_sig; 
  double _met_tc_et,_met_tc_px, _met_tc_py, _met_tc_phi, _met_tc_set, _met_tc_sig; 
  double _met_pf_et,_met_pf_px, _met_pf_py, _met_pf_phi, _met_pf_set, _met_pf_sig; 

  // Muons
  int _muons_N;
  TClonesArray * m_muons;
  int _muons_charge[20];
  int _muons_istracker[20], _muons_isstandalone[20], _muons_isglobal[20];
  double _muons_dxy[20], _muons_dz[20], _muons_dxyPV[20], _muons_dzPV[20], _muons_normalizedChi2[20];
  int  _muons_NtrackerHits[20], _muons_NpixelHits[20], _muons_NmuonHits[20], _muons_Nmatches[20];
  int _muons_nTkIsoR03[20], _muons_nTkIsoR05[20];
  double _muons_tkIsoR03[20],_muons_tkIsoR05[20],_muons_emIsoR03[20],_muons_emIsoR05[20],_muons_hadIsoR03[20],_muons_hadIsoR05[20];
	
  double _muons_trkDxy[20], _muons_trkDxyError[20], _muons_trkDxyB[20],
    _muons_trkDz[20], _muons_trkDzError[20], _muons_trkDzB[20], _muons_trkChi2PerNdof[20], 
    _muons_trkCharge[20],_muons_trkNHits[20],_muons_trkNPixHits[20];
  // Tracker muon properties
  double _muons_trkmuArbitration[20],
    _muons_trkmu2DCompatibilityLoose[20],
    _muons_trkmu2DCompatibilityTight[20],
    _muons_trkmuOneStationLoose[20],
    _muons_trkmuOneStationTight[20],
    _muons_trkmuLastStationLoose[20],
    _muons_trkmuLastStationTight[20],
    _muons_trkmuOneStationAngLoose[20],
    _muons_trkmuOneStationAngTight[20],
    _muons_trkmuLastStationAngLoose[20],
    _muons_trkmuLastStationAngTight[20],
    _muons_trkmuLastStationOptimizedLowPtLoose[20],
    _muons_trkmuLastStationOptimizedLowPtTight[20];
	
  double _muons_caloCompatibility[20], _muons_segmentCompatibility[20], _muons_glbmuPromptTight[20] ; 
	
  double _muons_hzzIso[20], _muons_hzzIsoTk[20], _muons_hzzIsoEcal[20], _muons_hzzIsoHcal[20];
	
  // Calo Jets 
  TClonesArray * _m_jets_calo ;
  int _jets_calo_N;
  //double _jets_calo_E[100], _jets_calo_pT[100], _jets_calo_px[100], _jets_calo_py[100], _jets_calo_pz[100], _jets_calo_eta[100], _jets_calo_phi[100];
	
  // JPT Jets
  TClonesArray * _m_jets_jpt ;
  int _jets_jpt_N;
	
  // PF Jets
  TClonesArray * _m_jets_pf;
  int _jets_pf_N;
	
  double jets_pf_chargedHadEFrac[100], jets_pf_chargedEmEFrac[100], jets_pf_chargedMuEFrac[100]; 
  double jets_pf_neutralHadEFrac[100], jets_pf_neutralEmEFrac[100], jets_pf_PhotonEFrac[100];
	  
  int jets_pf_chargedHadMultiplicity[100], jets_pf_neutralHadMultiplicity[100];
  int jets_pf_chargedMultiplicity[100], jets_pf_neutralMultiplicity[100];
  int jets_pf_nConstituents[100];

  // SuperClusters
  //converter::SuperClusterToCandidate * sc_struct;
  // SC EB
  int _sc_hybrid_N; 
  double _sc_hybrid_E[25], _sc_hybrid_Et[25], _sc_hybrid_Eta[25], _sc_hybrid_Phi[25]; 
  int _sc_hybrid_outOfTimeSeed[25],_sc_hybrid_severityLevelSeed[25];
  double _sc_hybrid_e1[25], _sc_hybrid_e33[25];
  double _sc_hybrid_he[25], _sc_hybrid_sigmaietaieta[25];
  double _sc_hybrid_hcalDepth1TowerSumEt_dr03[25], _sc_hybrid_hcalDepth2TowerSumEt_dr03[25];
  double _sc_hybrid_ecalRecHitSumEt_dr03[25];
  double _sc_hybrid_trkiso_dr03[25];

  int _sc_hybrid_RCTeta[25];
  int _sc_hybrid_RCTphi[25];
  int _sc_hybrid_RCTL1iso[25];
  int _sc_hybrid_RCTL1noniso[25];
  int _sc_hybrid_TTetaVect[25][50], _sc_hybrid_TTphiVect[25][50];
  double _sc_hybrid_TTetVect[25][50];
  int _sc_hybrid_RCTetaVect[25][10], _sc_hybrid_RCTphiVect[25][10], _sc_hybrid_RCTL1isoVect[25][10], _sc_hybrid_RCTL1nonisoVect[25][10];
  double _sc_hybrid_RCTetVect[25][10];

  // SC EE
  int _sc_multi55_N; 
  double _sc_multi55_E[25], _sc_multi55_Et[25], _sc_multi55_Eta[25], _sc_multi55_Phi[25];
  double _sc_multi55_he[25], _sc_multi55_sigmaietaieta[25]; 
  double _sc_multi55_hcalDepth1TowerSumEt_dr03[25], _sc_multi55_hcalDepth2TowerSumEt_dr03[25];
  double _sc_multi55_ecalRecHitSumEt_dr03[25];
  double _sc_multi55_trkiso_dr03[25];

  int _sc_multi55_RCTeta[25];
  int _sc_multi55_RCTphi[25];
  int _sc_multi55_RCTL1iso[25];
  int _sc_multi55_RCTL1noniso[25];
  int _sc_multi55_TTetaVect[25][50], _sc_multi55_TTphiVect[25][50];
  double _sc_multi55_TTetVect[25][50];
  int _sc_multi55_RCTetaVect[25][10], _sc_multi55_RCTphiVect[25][10], _sc_multi55_RCTL1isoVect[25][10], _sc_multi55_RCTL1nonisoVect[25][10];
  double _sc_multi55_RCTetVect[25][10];

  //add TIP/LIP/IP variables
  double muons_Tip[20],muons_Lip[20],muons_STip[20],muons_SLip[20],muons_TipSignif[20],muons_LipSignif[20],muons_Significance3D[20],muons_Value3D[20],muons_Error3D[20] ;
  double ele_Tip[10],ele_Lip[10],ele_STip[10],ele_SLip[10],ele_TipSignif[10],ele_LipSignif[10],ele_Significance3D[10],ele_Value3D[10],ele_Error3D[10];

  // NEW H/E
  double _ele_he_00615_0[10], _ele_he_005_0[10],_ele_he_005_1[10],_ele_he_005_15[10], _ele_he_01_0[10], _ele_he_01_1[10], _ele_he_01_15[10], _ele_he_015_1[10],_ele_he_015_15[10];
  //_ele_he_015_0[10] 
	

  // for H/E
  edm::Handle<CaloTowerCollection> * towersH_ ;
  edm::InputTag hcalTowers_ ;
  EgammaTowerIsolation * towerIso1_ ;
  EgammaTowerIsolation * towerIso2_ ;
  double hOverEConeSize_ ;
  double hOverEPtMin_ ;            

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////    NADIR STUFF ////////
  ////////////////////////////////

  // Booleans
  bool nadGetL1M_ ; // true => get two collections of L1 candidates (standard / "cleaned") ; false => get only the standard one
  bool nadGetTP_ ; // true => get the standard trigger primitives
  bool nadGetTP_Modif_ ; // true => get the modified collection of trigger primitives (zeroing by hand)
  bool nadGetTP_Emul_ ; // true => get the emulated collection of trigger primitives (Jackson-Zabi's sFGVB+zeroing emulator)
  bool PrintDebug_ , PrintDebug_HLT_;

  // handles to get the TPs
  edm::Handle<EcalTrigPrimDigiCollection> * ecal_tp_;
  edm::Handle<EcalTrigPrimDigiCollection> * ecal_tpM_;
  EcalTrigTowerDetId TPtowid_;
  EcalTrigTowerDetId TPtowidM_;

  // tags
  edm::InputTag EcalRecHitCollectionEB_ ;
  edm::InputTag tpCollectionNormal_ ;
  edm::InputTag tpCollectionModif_ ;
  edm::InputTag tpEmulatorCollection_ ;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  edm::InputTag EleID_VeryLooseTag_ ;
  edm::InputTag EleID_LooseTag_ ;
  edm::InputTag EleID_MediumTag_ ;
  edm::InputTag EleID_TightTag_ ; 
  edm::InputTag EleID_SuperTightTag_ ; 
  edm::InputTag EleID_HyperTight1Tag_ ;
  edm::InputTag EleID_HyperTight2Tag_ ;
  edm::InputTag EleID_HyperTight3Tag_ ;
  edm::InputTag EleID_HyperTight4Tag_ ;
  edm::InputTag EleIso_TdrHzzTkMapTag_, EleIso_TdrHzzHcalMapTag_ ;
  edm::InputTag EleIso_Eg4HzzTkMapTag_, EleIso_Eg4HzzEcalMapTag_,EleIso_Eg4HzzHcalMapTag_ ;
	
  edm::InputTag EleTag_;
  edm::InputTag MuonTag_;
  edm::InputTag MuonIso_HzzMapTag_ ;
  edm::InputTag MuonIsoTk_HzzMapTag_ ;
  edm::InputTag MuonIsoEcal_HzzMapTag_ ;
  edm::InputTag MuonIsoHcal_HzzMapTag_ ;

  edm::InputTag SeedTag_;
  edm::InputTag MCTag_;
  edm::InputTag TkPTag_;
  edm::InputTag CaloJetTag_;
  edm::InputTag JPTJetTag_;
  edm::InputTag PFJetTag_;
  edm::InputTag VerticesTag_;
  edm::InputTag dcsTag_;
	
  // Trigger Stuff
  edm::InputTag HLTTag_; 
  edm::InputTag triggerEventTag_;
  std::vector<std::string > HLT_ElePaths_;
  std::vector<std::string > HLT_MuonPaths_;
  std::vector<edm::InputTag > HLT_Filters_;

  //Pile-up
  edm::InputTag PileupSrc_;
  edm::InputTag RhoCorrection_, SigmaRhoCorrection_, BetaCorrection_;

  std::string type_;	
  bool aod_;	
  bool simulation_;
  bool fillsc_;

  const EcalRecHit getRecHit(DetId id, const EcalRecHitCollection *recHits);
	
  std::string gtRecordCollectionTag_ ;
  EcalClusterFunctionBaseClass* funcbase_;
  std::string funcname_;
  bool useBeamSpot_ ;
};
