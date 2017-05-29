#include "EGamma/ECGelec/plugins/SpikeRemoval.h"


#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

using namespace std;
bool SpikeRemoval::selection(edm::Ref<reco::GsfElectronCollection> eleRef)
{
  const reco::CaloCluster & seedCluster = *(eleRef->superCluster()->seed()) ;  
  const EcalRecHitCollection * reducedRecHits = 0 ;
  if (eleRef->isEB())
    { reducedRecHits = reducedEBRecHits.product() ; }
  else
    { reducedRecHits = reducedEERecHits.product() ; }
  
  double e1 = EcalClusterTools::eMax(seedCluster,reducedRecHits)  ;
  //double e2x2 = EcalClusterTools::e2x2(seedCluster,reducedRecHits,topology)  ;
  double e3x3 = EcalClusterTools::e3x3(seedCluster,reducedRecHits,topology)  ;
 

  //  cout << "e1 = " << e1 << " e1/e3*3 = " << e1/e3x3 << endl;
  if (fabs(e1)>5. && e1/e3x3>0.95) { /*cout<< "false" << endl;*/ return false; }
  
  //  cout<< "true" <<endl;
  return true;
}
