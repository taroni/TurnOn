#include "EGamma/ECGelec/plugins/MuonHzzIsolationProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

// Tracker tracks
#include "DataFormats/TrackReco/interface/Track.h"

// Electrons
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

// Primary Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//
#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

// constructor
MuonHzzIsolationProducer::MuonHzzIsolationProducer(const edm::ParameterSet& pset) {

  //primary vertex
  PVtag_                    = pset.getParameter<edm::InputTag>("PVLabel");     
  //deposits
  theECALIsoDepositLabel    = pset.getParameter<edm::InputTag>("ECALIsoDepositLabel");        
  theHCALIsoDepositLabel    = pset.getParameter<edm::InputTag>("HCALIsoDepositLabel");        
  theHOCALIsoDepositLabel   = pset.getParameter<edm::InputTag>("HOCALIsoDepositLabel");      
  theTrackerIsoDepositLabel = pset.getParameter<edm::InputTag>("TrackerIsoDepositLabel");   
  // input objects
  muonsTag_                 = pset.getParameter<edm::InputTag>("MuonsLabel");
  tracksTag_                = pset.getParameter<edm::InputTag>("TracksLabel");
  electronTag_              = pset.getParameter<edm::InputTag>("ElectronsLabel");
  // algo parameters
  mainConeSize              = pset.getParameter<double>("isolationCone");
  vetoConeSize              = pset.getParameter<double>("isolationConeVeto");
  trkIsoW                   = pset.getParameter<double>("trkIsoWeight");           
  ecalW                     = pset.getParameter<double>("ecalWeight");             
  hcalW                     = pset.getParameter<double>("hcalWeight");
  // isolation cut
  isocut                    = pset.getParameter<double>("isolationcut");

  string alias;
  string iName = "MuonIsolation";

  produces<vector<double> >( alias = iName + "X" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "CalIso" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "ECalIso" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "HCalIso" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "SumpT" ).setBranchAlias( alias );
  produces<reco::MuonCollection>();
  produces<edm::ValueMap<double> >();

  produces<edm::ValueMap<double> >(alias = "Tk" ).setBranchAlias( alias);
  produces<edm::ValueMap<double> >(alias = "Ecal" ).setBranchAlias( alias);
  produces<edm::ValueMap<double> >(alias = "Hcal" ).setBranchAlias( alias);

 //  produces<edm::ValueMap<double> >( alias = "TrkIsoWeight" ).setBranchAlias( alias );           
//   produces<edm::ValueMap<double> >( alias = "EcalWeight" ).setBranchAlias( alias ); 
//   produces<edm::ValueMap<double> >( alias = "HcalWeight" ).setBranchAlias( alias ); 
//   produces<edm::ValueMap<double> >( alias = "MainConeSize" ).setBranchAlias( alias );
//   produces<edm::ValueMap<double> >( alias = "VetoConeSize" ).setBranchAlias( alias );
  

}


// destructor
MuonHzzIsolationProducer::~MuonHzzIsolationProducer() {

}


void MuonHzzIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  /*auto_ptr<vector<double> > XAfterVetos( new vector<double> );
  auto_ptr<vector<double> > CalIsoAfterVetos( new vector<double> );
  auto_ptr<vector<double> > EcalAfterVetos( new vector<double> );
  auto_ptr<vector<double> > HcalAfterVetos( new vector<double> );
  auto_ptr<vector<double> > sumPtAfterVetos( new vector<double> );
  auto_ptr<reco::MuonCollection> Isolmuon( new reco::MuonCollection );
  auto_ptr<edm::ValueMap<double> > IsoMuMap(new edm::ValueMap<double> ());
  auto_ptr<edm::ValueMap<double> > IsoMuTkMap(new edm::ValueMap<double> ());
  auto_ptr<edm::ValueMap<double> > IsoMuEcalMap(new edm::ValueMap<double> ());
  auto_ptr<edm::ValueMap<double> > IsoMuHcalMap(new edm::ValueMap<double> ());*/
  
  std::unique_ptr<vector<double> > XAfterVetos( new vector<double> );
  std::unique_ptr<vector<double> > CalIsoAfterVetos( new vector<double> );
  std::unique_ptr<vector<double> > EcalAfterVetos( new vector<double> );
  std::unique_ptr<vector<double> > HcalAfterVetos( new vector<double> );
  std::unique_ptr<vector<double> > sumPtAfterVetos( new vector<double> );
  std::unique_ptr<reco::MuonCollection> Isolmuon( new reco::MuonCollection );
  std::unique_ptr<edm::ValueMap<double> > IsoMuMap(new edm::ValueMap<double> ());
  std::unique_ptr<edm::ValueMap<double> > IsoMuTkMap(new edm::ValueMap<double> ());
  std::unique_ptr<edm::ValueMap<double> > IsoMuEcalMap(new edm::ValueMap<double> ());
  std::unique_ptr<edm::ValueMap<double> > IsoMuHcalMap(new edm::ValueMap<double> ());

  edm::ValueMap<double>::Filler filler(*IsoMuMap);
  edm::ValueMap<double>::Filler fillerTk(*IsoMuTkMap);
  edm::ValueMap<double>::Filler fillerEcal(*IsoMuEcalMap);
  edm::ValueMap<double>::Filler fillerHcal(*IsoMuHcalMap);

 //  auto_ptr<edm::ValueMap<double> > TrkIsoWeightMap(new edm::ValueMap<double> ());       
//   edm::ValueMap<double>::Filler TrkIsoWeightFiller(*TrkIsoWeightMap);                         
//   auto_ptr<edm::ValueMap<double> > EcalWeightMap(new edm::ValueMap<double> ());       
//   edm::ValueMap<double>::Filler EcalWeightFiller(*EcalWeightMap);  
//   auto_ptr<edm::ValueMap<double> > HcalWeightMap(new edm::ValueMap<double> ());       
//   edm::ValueMap<double>::Filler HcalWeightFiller(*HcalWeightMap); 
//   auto_ptr<edm::ValueMap<double> > MainConeSizeMap(new edm::ValueMap<double> ());       
//   edm::ValueMap<double>::Filler MainConeSizeFiller(*MainConeSizeMap); 
//   auto_ptr<edm::ValueMap<double> > VetoConeSizeMap(new edm::ValueMap<double> ());       
//   edm::ValueMap<double>::Filler VetoConeSizeFiller(*VetoConeSizeMap); 

  // Get tracks
  Handle<edm::View<reco::Track> > tracks;
  iEvent.getByLabel(tracksTag_.label(), tracks);

  //get electrons
  Handle<edm::View<reco::GsfElectron> > EleCandidates;
  iEvent.getByLabel(electronTag_.label(), EleCandidates);

  // Get muons used to build a map
  Handle<edm::View<reco::Muon> > AllMuons0;
  iEvent.getByLabel(muonsTag_.label(), AllMuons0);
  vector<const reco::Muon*> Muons;   
  edm::View<reco::Muon>::const_iterator Muon_it;  

  for(Muon_it = AllMuons0 -> begin(); Muon_it != AllMuons0 -> end(); ++Muon_it){
    Muons.push_back(&(*Muon_it));
  }

  //**************** MAPS ******************
  //to get the deposits inside an isolation cone
  
  Handle<reco::IsoDepositMap> ecalIso;    
  iEvent.getByLabel(theECALIsoDepositLabel, ecalIso);

  Handle<reco::IsoDepositMap> hcalIso;
  iEvent.getByLabel(theHCALIsoDepositLabel, hcalIso);
  
  //Handle<reco::IsoDepositMap> hocalIso;
  //iEvent.getByLabel(theHOCALIsoDepositLabel, hocalIso);
  
  Handle<reco::IsoDepositMap> trackerIso;
  iEvent.getByLabel(theTrackerIsoDepositLabel, trackerIso);

  //************** weights ******************

  double delta_R = 0.;
  double depVeto03Trk = 0., depVeto03Ecal = 0., depVeto03Hcal = 0., depVeto03HOcal = 0., depVeto03CalIso = 0., depVeto03X = 0.;
  double TrkIso = 0., CalIso = 0., X1 = 0., X2 = 0., X3 = 0., X4 = 0., X5 = 0.;
  
  vector<double> sumPtOverMuPtAfterVetos, HOcalAfterVetos; 
  vector<double> X1Vec, X2Vec, X3Vec, X4Vec, X5Vec;

  unsigned int index=0;

  size_t n = AllMuons0->size();
  std::vector<float> iso(n), isoTk(n), isoEcal(n), isoHcal(n);
  //std::vector<float> trkIsoWeight(n), ecalWeight(n), hcalWeight(n), mainCone(n), vetoCone(n);

  for(unsigned int init=0; init<n; init++){
    
    iso[init]          = 0.;
    //  trkIsoWeight[init] = 0.;
    //     ecalWeight[init]   = 0.;
    //     hcalWeight[init]   = 0.;
    //     mainCone[init]     = 0.;
    //     vetoCone[init]     = 0.;
  }

  // primary vertex
  edm::Handle<reco::VertexCollection> recoPVCollection;
  iEvent.getByLabel(PVtag_, recoPVCollection);
  reco::Vertex primVertex;
  math::XYZPoint pVertex(0., 0., 0.);
  bool pvfound = (recoPVCollection->size() != 0);
  if(pvfound) {
    primVertex = recoPVCollection->front();
    pVertex = math::XYZPoint(primVertex.position().x(), primVertex.position().y(), primVertex.position().z());
  }


  for (unsigned h = 0 ; h < AllMuons0->size() ; h++){ 

    //cout << "Candidate muon found for tight isolation" << endl;
    
    //********* veto against the other muons in an isolation cone ************
    edm::Ref<edm::View<reco::Muon> > muonRef(AllMuons0,h);
        
    reco::IsoDeposit depTracker((*trackerIso)[muonRef]); //get sumPt around h-th muon
    reco::IsoDeposit depEcal((*ecalIso)[muonRef]);       //get ECal dep around h-th muon
    reco::IsoDeposit depHcal((*hcalIso)[muonRef]);       //get HCal dep around h-th muon
    //reco::IsoDeposit depHOcal =((*hocalIso)[muonRef]);     //get HoCal dep around h-th muon	
    
    reco::IsoDeposit::Vetos myVetoTrkVec;                              //vector of vetos    
    reco::IsoDeposit::Direction dirMu1 = depTracker.direction();      //direction of the muon you're isolating
    
    for(unsigned int k = 0; k < AllMuons0->size(); k++) {               //inner loop over muons
      edm::Ref<edm::View<reco::Muon> > muonRefbis(AllMuons0,k);

      if(k == h)
	continue;      //skip same muons
      reco::IsoDeposit depTracker2((*trackerIso)[muonRefbis]); //get sumPt around h-th muon
      reco::IsoDeposit::Direction dirMu2 = depTracker2.direction();                     
      
      delta_R = dirMu1.deltaR(dirMu2);         //calculate delta_R
      
      if(delta_R < mainConeSize) {             //i.e. if another muon falls into the cone
	
	IsoDeposit::Veto myVetoTrk(dirMu2, vetoConeSize);
	myVetoTrkVec.push_back(myVetoTrk);   //fill veto vector for tracker
      }
    }

    for (edm::View<reco::GsfElectron>::const_iterator cand = EleCandidates->begin(); cand != EleCandidates->end(); ++cand) {

      reco::Particle::LorentzVector trackele(cand ->px(), cand ->py(), cand -> pz(), cand -> p());      
      reco::IsoDeposit::Direction dirTrackele(trackele.eta(), trackele.phi());

      float dRele = dirMu1.deltaR(dirTrackele);
      
      if(dRele < mainConeSize ) {
	
	IsoDeposit::Veto myVetoTrk3(dirTrackele, vetoConeSize);
	myVetoTrkVec.push_back(myVetoTrk3);   //fill veto vector for tracker	
      }
    }
    
    for (edm::View<reco::Track>::const_iterator itTrack = tracks -> begin(); itTrack != tracks -> end(); itTrack++) {    
      
      bool goodTrack = testTrackerTrack(itTrack, AllMuons0->at(h), &pVertex);
      
      reco::Particle::LorentzVector track(itTrack -> px(), itTrack -> py(), itTrack -> pz(), itTrack -> p());      
      reco::IsoDeposit::Direction dirTrack(track.eta(), track.phi());
      
      float dR = dirMu1.deltaR(dirTrack);
      
      if(dR < mainConeSize && goodTrack == false) {
	
	IsoDeposit::Veto myVetoTrk2(dirTrack, vetoConeSize);
	myVetoTrkVec.push_back(myVetoTrk2);   //fill veto vector for tracker	
      }
    }


    //****** redefining all the isolation variables after vetos
    depVeto03Trk = depTracker.depositWithin(mainConeSize, myVetoTrkVec, false);    //sumPt
    depVeto03Ecal = depEcal.depositWithin(mainConeSize, myVetoTrkVec, false);      //ECal
    depVeto03Hcal = depHcal.depositWithin(mainConeSize, myVetoTrkVec, false);      //HCal
    //depVeto03HOcal = depHOcal.depositWithin(mainConeSize, myVetoTrkVec, false);    //HOCal
    depVeto03CalIso = ecalW * depVeto03Ecal + hcalW * depVeto03Hcal;               //CalIso   
    depVeto03X = trkIsoW * depVeto03Trk + depVeto03CalIso;                         //X        
    
    TrkIso = depVeto03Trk;
    CalIso = depVeto03CalIso;
    
    X1 = TrkIso + CalIso;                                //linear1
    X2 = 1.5 * TrkIso + 0.5 * CalIso;                    //linear2
    X3 = TrkIso * CalIso;                                //hyperbolic
    X4 = sqrt( (TrkIso*TrkIso) + (CalIso*CalIso));       //circular
    X5 = sqrt( 4*(TrkIso*TrkIso) + (CalIso*CalIso) );    //elliptical
    
    double ptmu = AllMuons0->at(h).pt();  
    if(ptmu != 0.) {                  //do not divide by zero!
      double sumptovermupt = depVeto03Trk/ptmu;                   //sumPt/pT_mu
      sumPtOverMuPtAfterVetos.push_back(sumptovermupt);
    }
    
    //cout << "X variable:" << depVeto03X << endl;
	//if ( depVeto03X < isocut ){
      sumPtAfterVetos->push_back(depVeto03Trk);        //sumPt
      EcalAfterVetos->push_back(depVeto03Ecal);        //ECal
      HcalAfterVetos->push_back(depVeto03Hcal);        //HCal
      HOcalAfterVetos.push_back(depVeto03HOcal);       //HOCal
      CalIsoAfterVetos->push_back(depVeto03CalIso);    //CalIso
      XAfterVetos->push_back(depVeto03X);              //X
      X1Vec.push_back(X1);                             //X1
      X2Vec.push_back(X2);                             //X2
      X3Vec.push_back(X3);                             //X3
      X4Vec.push_back(X4);                             //X4
      X5Vec.push_back(X5);                             //X5
      
      Isolmuon->push_back(AllMuons0->at(h));

      iso[index] = float(depVeto03X);
      isoTk[index] = float(depVeto03Trk);
      isoEcal[index] = float(depVeto03Ecal);
      isoHcal[index] = float(depVeto03Hcal);

      //trkIsoWeight[index] = float(trkIsoW); 
      //ecalWeight[index] = float(ecalW);     
      //hcalWeight[index] = float(hcalW); 
      //mainCone[index] = float(mainConeSize);     
      //vetoCone[index] = float(vetoConeSize);
      
      //cout << "trkIsoW, ecalW, hcalW, mainCone, VetoCone in MuonIsolationProducer = " 
	  // << trkIsoW << ", " << ecalW << ", " << hcalW << ", "
	  // << mainConeSize << ", " << vetoConeSize << endl;
      //cout << "iso[" << index << "] = " << iso[index] << endl;
    //}
    index++;
  }

  // filling map
  //cout << "Size of iso=" << iso.size() << endl;
  filler.insert(AllMuons0, iso.begin(), iso.end());
  filler.fill();  

  fillerTk.insert(AllMuons0, isoTk.begin(), isoTk.end());
  fillerTk.fill();  

  fillerEcal.insert(AllMuons0, isoEcal.begin(), isoEcal.end());
  fillerEcal.fill();  

  fillerHcal.insert(AllMuons0, isoHcal.begin(), isoHcal.end());
  fillerHcal.fill();  

  // TrkIsoWeightFiller.insert(AllMuons0, trkIsoWeight.begin(), trkIsoWeight.end());           
//   TrkIsoWeightFiller.fill();                                                                    
  
//   EcalWeightFiller.insert(AllMuons0, ecalWeight.begin(), ecalWeight.end());           
//   EcalWeightFiller.fill(); 

//   HcalWeightFiller.insert(AllMuons0, hcalWeight.begin(), hcalWeight.end());           
//   HcalWeightFiller.fill(); 

//   MainConeSizeFiller.insert(AllMuons0, mainCone.begin(), mainCone.end());           
//   MainConeSizeFiller.fill();

//   VetoConeSizeFiller.insert(AllMuons0, vetoCone.begin(), vetoCone.end());           
//   VetoConeSizeFiller.fill(); 
  

  const string & isoName = "MuonIsolation";
  /*iEvent.put( XAfterVetos,      isoName + "X" );
  iEvent.put( CalIsoAfterVetos, isoName + "CalIso" );
  iEvent.put( EcalAfterVetos,   isoName + "ECalIso" );  
  iEvent.put( HcalAfterVetos,   isoName + "HCalIso" );
  iEvent.put( sumPtAfterVetos,  isoName + "SumpT" );*/
  iEvent.put( std::move(XAfterVetos),       isoName + "X" );
  iEvent.put( std::move(CalIsoAfterVetos),  isoName + "CalIso" );
  iEvent.put( std::move(EcalAfterVetos),    isoName + "ECalIso" );  
  iEvent.put( std::move(HcalAfterVetos),    isoName + "HCalIso" );
  iEvent.put( std::move(sumPtAfterVetos),   isoName + "SumpT" );

  const string iName = "";
  iEvent.put( std::move(Isolmuon), iName );
  iEvent.put( std::move(IsoMuMap), iName );
  iEvent.put( std::move(IsoMuTkMap), "Tk" );
  iEvent.put( std::move(IsoMuEcalMap), "Ecal" );
  iEvent.put( std::move(IsoMuHcalMap), "Hcal" );
//   iEvent.put( TrkIsoWeightMap, "TrkIsoWeight");                       
//   iEvent.put( EcalWeightMap, "EcalWeight");     
//   iEvent.put( HcalWeightMap, "HcalWeight");  
//   iEvent.put( MainConeSizeMap, "MainConeSize"); 
//   iEvent.put( VetoConeSizeMap, "VetoConeSize"); 

}


bool MuonHzzIsolationProducer::testTrackerTrack(edm::View<reco::Track>::const_iterator iterTrack,  reco::Muon muon, math::XYZPoint* pv) {

  /************* quality track requirements ************/  

  bool keep = true;
  // Extract properties at Vertex
  float vz  = iterTrack -> vz();
  float edz = iterTrack -> dzError();
  float d0  = iterTrack -> dxy(*pv);   //d0();
  float ed0 = iterTrack -> dxyError(); //d0Error();
  // Difference with lepton/Z vertex
  float dz  = muon.vertex().z() - vz;
  
  dz  = fabs(dz);                //impact parameter in the (r, z) plane
  edz = fabs(edz);               //error on dz 
  d0  = fabs(d0);                //impact parameter in the (r, phi) plane
  ed0 = fabs(ed0);               //error on d0   

  reco::Particle::LorentzVector track(iterTrack -> px(), iterTrack -> py(), iterTrack -> pz(), iterTrack -> p());
  float track_pt =  track.pt();
  
  const reco::HitPattern& p = iterTrack->hitPattern();
//  int nhits = p.numberOfHits(reco::HitPattern::MISSING_INNER_HITS); 
  int nhits = p.numberOfHits(reco::HitPattern::TRACK_HITS); 
  //int nhits = iterTrack -> recHitsSize(); 

  if ( nhits < 8 ) {
    if ( track_pt < 1.) return false;
    if ( d0 > 0.04) return false;
    if ( dz > 0.50) return false;
    if ( d0 / ed0 > 7.) return false;
    if ( dz / edz > 7.) return false;
  }
  else if ( nhits < 10 ) {
    if ( track_pt < 1.) return false;
    if ( d0 > 0.20) return false;
    if ( dz > 2.00) return false;
    if ( d0 / ed0 > 10.) return false;
    if ( dz / edz > 10.) return false;
  }
  else {
    if ( track_pt < 1.) return false;
    if ( d0 > 1.00) return false;
    if ( dz > 5.00) return false;
  }
  return keep;
}
