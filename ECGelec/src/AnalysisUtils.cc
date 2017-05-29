#include "EGamma/ECGelec/interface/AnalysisUtils.h"
// Electron
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// ===================================================================================================
//  Constructor
// ===================================================================================================
AnalysisUtils::AnalysisUtils() {

  _is_1lepton  = false;
  _is_2leptons = false;
  _is_3leptons = false;

} // end of constructor

// ===================================================================================================
//  Destructor
// ===================================================================================================
AnalysisUtils::~AnalysisUtils() {}



// ===================================================================================================
bool AnalysisUtils::doSkim(const edm::Event& iEvent, const edm::EventSetup& iSetup,
			   bool isEleID_, bool isMuonID_, double lep_ptLow_, double lep_ptHigh_, int nLep_ptLow_, int nLep_ptHigh_) 
// ===================================================================================================
{

  // ------------------------------------------------
  // GSF Electron	
  // ------------------------------------------------
  edm::Handle<reco::GsfElectronCollection> EleHandle;
  iEvent.getByLabel("gsfElectrons", EleHandle);
  
  // loop on electrons
  int nEle_ptLow = 0;
  int nEle_ptHigh = 0;
  for(unsigned int eleIt = 0; eleIt < EleHandle->size(); ++eleIt)
    {
      
      reco::GsfElectronRef eleRef(EleHandle, eleIt);
      
      //electron ID
      bool eleIDBool = true ;
      if (isEleID_) {
	//barrel
	if (eleRef->isEB()) {
	  if ( fabs(eleRef->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) eleIDBool = false ;
	  if ( fabs(eleRef->deltaPhiSuperClusterTrackAtVtx()) > 0.08 )  eleIDBool = false ;
	  if ( eleRef->hadronicOverEm()  > 0.15 ) eleIDBool = false ;
	  if ( eleRef->sigmaIetaIeta() > 0.01 ) eleIDBool = false ;
	  
	}
	//endcap
	if (eleRef->isEE()) {
	  if ( fabs(eleRef->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) eleIDBool = false ;
	  if ( fabs(eleRef->deltaPhiSuperClusterTrackAtVtx()) > 0.07)  eleIDBool = false ;
	  if ( eleRef->hadronicOverEm()  > 0.07 ) eleIDBool = false ;
	  if ( eleRef->sigmaIetaIeta() > 0.03 ) eleIDBool = false ;
	  
	}
		}
      // count electrons
      if( eleIDBool && eleRef->pt() > lep_ptLow_)  ++nEle_ptLow;
      if( eleIDBool && eleRef->pt() > lep_ptHigh_) ++nEle_ptHigh;
      
    } // end loop on electrons
  
  // ------------------------------------------------
  // Muons
  // ------------------------------------------------
  edm::Handle<reco::MuonCollection> MuonHandle;
  iEvent.getByLabel("muons", MuonHandle);
  
  // loop on muons
  int nMuon_ptLow = 0;
  int nMuon_ptHigh = 0;
  for(unsigned int muonIt = 0; muonIt < MuonHandle->size(); ++muonIt)
    { 	
      edm::Ref<reco::MuonCollection> muonRef(MuonHandle,muonIt);
      
      //muon ID
      bool muonIDBool = true ;
      if (isMuonID_) {
	if ( muonRef->isGlobalMuon() == 0 ) muonIDBool = false ;
		}
      // count muons
      if( muonIDBool && muonRef->pt() > lep_ptLow_) ++nMuon_ptLow;
      if( muonIDBool && muonRef->pt() > lep_ptHigh_) ++nMuon_ptHigh;
      
    } // end loop on muons
  
  
  // print out results
  if ( ( nEle_ptLow + nMuon_ptLow >= nLep_ptLow_ ) &&  
       ( nEle_ptHigh + nMuon_ptHigh >= nLep_ptHigh_ )  ){
    //std::cout << ">>>SkimLeptonStudies::run=" << runId << "::lumi=" 
    //<< lumiId << "::eventId=" 
    //<< eventId << "::Found at least "
    //<< nLep_ptLow_ << " electrons with pt > " << lep_ptLow_ << " AND " 
    //<< nLep_ptHigh_ << " electrons with pt > " << lep_ptHigh_ << std::endl; 
    return true;
  }
  
  return false;
  
} // end of doSkim
