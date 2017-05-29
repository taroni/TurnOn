#ifndef SPIKE_REMOVAL
#define SPIKE_REMOVAL

#include <memory>
#include <algorithm>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"


#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"




struct SpikeRemoval{
  public:
  explicit SpikeRemoval(const edm::ParameterSet& iConfig) {}

  virtual ~SpikeRemoval() {}

  // Collections to be selected
  typedef reco::GsfElectronCollection collection;
  typedef std::vector<reco::GsfElectronRef> container ; 
  //typedef std::vector<reco::GsfElectron> container ; 
  typedef container::const_iterator const_iterator;

  //define iterators with above typedef
  const_iterator begin () const { return selected_.begin () ; }
  const_iterator end () const { return  selected_.end () ; }

  void select(edm::Handle<reco::GsfElectronCollection> electrons, 
	      const edm::Event& iEvent , 
	      const edm::EventSetup& iEs)
  {
    selected_.clear();

    //calo topology
    unsigned long long cacheIDTopo_ = 0;
    edm::ESHandle<CaloTopology> theCaloTopo;
    if (cacheIDTopo_!=iEs.get<CaloTopologyRecord>().cacheIdentifier()){
      cacheIDTopo_=iEs.get<CaloTopologyRecord>().cacheIdentifier();
      iEs.get<CaloTopologyRecord>().get(theCaloTopo);
    }
    topology = theCaloTopo.product() ;
    
    // reduced rechits
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEB"), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag("ecalRecHit:EcalRecHitsEE"), reducedEERecHits ) ;


    // Loop over electrons
    unsigned int i = 0 ;
    for ( reco::GsfElectronCollection::const_iterator eleIt = electrons->begin () ;
	  eleIt != electrons->end () ;
	  ++eleIt )
      {
	edm::Ref<reco::GsfElectronCollection> electronRef(electrons,i);
	if (selection(electronRef))
	  selected_.push_back (electronRef) ;
	//selected_.push_back (*eleIt) ;
	++i;
      }
  }
  
  private:
  container selected_ ;
  const CaloTopology * topology ;
  edm::Handle< EcalRecHitCollection > reducedEBRecHits;
  edm::Handle< EcalRecHitCollection > reducedEERecHits;

   
  bool selection(edm::Ref<reco::GsfElectronCollection> eleRef);
  
};

#endif

