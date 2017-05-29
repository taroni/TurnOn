#ifndef MuonHzzIsolationProducer_h
#define MuonHzzIsolationProducer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

class MuonHzzIsolationProducer : public edm::EDProducer {

 public:
  explicit MuonHzzIsolationProducer(const edm::ParameterSet&);
  ~MuonHzzIsolationProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  bool testTrackerTrack(edm::View<reco::Track>::const_iterator, reco::Muon, math::XYZPoint*);

  edm::InputTag PVtag_;
  edm::InputTag muonsTag_;
  edm::InputTag theECALIsoDepositLabel;    // EM calorimeter Isolation deposit label
  edm::InputTag theHCALIsoDepositLabel;    // Hadron calorimeter Isolation deposit label
  edm::InputTag theHOCALIsoDepositLabel;   // Outer calorimeter Isolation deposit label
  edm::InputTag theTrackerIsoDepositLabel; // Tracker Isolation deposit label 
  edm::InputTag electronTag_;
  edm::InputTag tracksTag_;
  double mainConeSize;
  double isocut;
  //  double ptMin;
  double vetoConeSize;
  // double maxDz;
  double trkIsoW;   
  double ecalW;     
  double hcalW; 

};

#endif
