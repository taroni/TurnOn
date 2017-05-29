#ifndef DEF_CommonFunctions
#define DEF_CommonFunctions

#include "EGamma/ECGelec/plugins/ToInclude.h"

typedef math::XYZTLorentzVector LorentzVector ;
#define arrondi(a) (floor((a) + 0.5))

bool VBTFcuts(TString asked,
	      double pt, double eT, double eta, double etaEle, double ele_tkSumPt_dr03, double ele_ecalRecHitSumEt_dr03,
	      double ele_hcalDepth1TowerSumEt_dr03, double ele_hcalDepth2TowerSumEt_dr03, double ele_expected_inner_hits,
	      double ele_deltaphiin, double ele_deltaetain, double ele_he, double ele_sigmaietaieta,
	      double ele_conv_dist, double ele_conv_dcot, double ele_fbrem, int ele_isConversion);

int convertEtaToTT(double eta);
int convertPhiToTT(double phi);

int convertEtaToRCT(double eta);
int convertPhiToRCT(double phi);
int getGCTRegionPhi(int ttphi);
int getGCTRegionEta(int tteta);
void setMomentum(TLorentzVector & myvector , const LorentzVector & mom);
bool IsConv (const reco::GsfElectron & eleRef);

#endif
