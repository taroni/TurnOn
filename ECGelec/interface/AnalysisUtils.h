#ifndef AnalysisUtils_h
#define AnalysisUtils_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//// Electron
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//// Muons
//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"


class AnalysisUtils
{
 public:

  AnalysisUtils();

  ~AnalysisUtils();

  bool doSkim(const edm::Event& iEvent, const edm::EventSetup& iSetup,
	      bool isEleID_, bool isMuonID_, double lep_ptLow_, double lep_ptHigh_, int nLep_ptLow_, int nLep_ptHigh_);

  bool Is1Lepton();
  bool Is2Leptons();
  bool Is3Leptons();

 private:

  bool _is_1lepton;
  bool _is_2leptons;
  bool _is_3leptons;

};

#endif

