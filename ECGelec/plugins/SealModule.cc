#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/Framework/src/MakeModuleHelper.h" //ADD
#include "FWCore/Framework/src/WorkerMaker.h" //ADD
#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h" //ADD
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h" //ADD

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "EGamma/ECGelec/plugins/RunSelect.h"
//#include "EGamma/ECGelec/plugins/SpikeRemoval.h"
#include "EGamma/ECGelec/plugins/ElectronL1Study.h"
#include "EGamma/ECGelec/plugins/ElectronJetsL1Study.h"
#include "EGamma/ECGelec/plugins/ElectronL1Study_TrigOnly.h"
#include "EGamma/ECGelec/plugins/ElectronL1Study_dumpHLT.h"


// typedef ObjectSelector<
//   SpikeRemoval
//   > SpikeRemovalSelector;

DEFINE_FWK_MODULE(RunSelect);
//DEFINE_FWK_MODULE(SpikeRemovalSelector);
DEFINE_FWK_MODULE(ElectronL1Study);
DEFINE_FWK_MODULE(ElectronJetsL1Study);
DEFINE_FWK_MODULE(ElectronL1Study_TrigOnly);
DEFINE_FWK_MODULE(ElectronL1Study_dumpHLT);
