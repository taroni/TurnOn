#include "EGamma/ECGelec/interface/CommonFunctions.h"

// ====================================================================================
// VBTFcuts : determine wether a given VBTF set of cuts is passed by the electron
// ====================================================================================
bool VBTFcuts(TString asked,
			       double pt, double eT, double eta, double etaEle, double ele_tkSumPt_dr03, double ele_ecalRecHitSumEt_dr03, 
			       double ele_hcalDepth1TowerSumEt_dr03, double ele_hcalDepth2TowerSumEt_dr03, double ele_expected_inner_hits,
			       double ele_deltaphiin, double ele_deltaetain, double ele_he, double ele_sigmaietaieta,
			       double ele_conv_dist, double ele_conv_dcot, double ele_fbrem, int ele_isConversion)
{
  //isolation variables
  double trackIsoRel03 = ele_tkSumPt_dr03 / pt;
  double ecalIsoRel03 = ele_ecalRecHitSumEt_dr03 / pt;
  double hcalIsoRel03 = (ele_hcalDepth1TowerSumEt_dr03+ele_hcalDepth2TowerSumEt_dr03) / pt;

  // define cuts
  // expectedInnerHits ; EB:|DPhiIn|,|DEtaIn|,H/E,|tkIso|,|hIso|,|eIso| ; EE:idem
  double cuts[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};

  double cuts95[13] = {1 , 0.8 , 0.007 , 0.15 , 0.15 , 0.12 , 2.00 , 0.7 , 0.01 , 0.07 , 0.08 , 0.05 , 0.06} ;
  double cuts80[13] = {0 , 0.06 , 0.004 , 0.04 , 0.09 , 0.10 , 0.07 , 0.03 , 0.007 , 0.025 , 0.04 , 0.025 , 0.05} ;
  double cuts60[13] = {0 , 0.025 , 0.004 , 0.025 , 0.04 , 0.03 , 0.04 , 0.02 , 0.005 , 0.025 , 0.025 , 0.02 , 0.02} ;

  cuts95[9] = cuts80[9] = cuts60[9] = 0.15 ; // 2011  

  // define HLT_Ele27 cuts
  //double trackIsoRel03_H = ele_tkSumPt_dr03 / eT;
  //double ecalIsoRel03_H = ele_ecalRecHitSumEt_dr03 / eT;
  //double hcalIsoRel03_H = (ele_hcalDepth1TowerSumEt_dr03+ele_hcalDepth2TowerSumEt_dr03) / eT;

  // HoE , sigmaietaieta , ecaliso/et , hcaliso/et , deta, dphi, tkiso/et
  double cutsHLT_Ele27_IdIso_EB[7] = {0.05 , 0.011 , 0.15 , 0.15 , 0.008 , 0.07 , 0.15};
  double cutsHLT_Ele27_IdIso_EE[7] = {0.05 , 0.031 , 0.075 , 0.075 , 0.008 , 0.05 , 0.075};

  double values[7] = { ele_he , ele_sigmaietaieta , ecalIsoRel03 , hcalIsoRel03 , 
		       fabs(ele_deltaetain) , fabs(ele_deltaphiin) , trackIsoRel03 } ;

  // "conversion" selection
  if(asked=="conv") {
    // ID
    if ( fabs(ele_deltaphiin) >= 0.1 ) return false;
    if ( ele_he >= 0.1 ) return false;
    if ( ele_fbrem <= -0.1 ) return false;
    if ( fabs(etaEle)<1.479 ) {
      if ( ele_sigmaietaieta <= 0.008 ) return false;
    }
    else {
      return false;
      if ( ele_sigmaietaieta <= 0.02 ) return false;
    }
    if( fabs(ele_deltaetain) >= 0.01 ) return false;

    // ISO
    if ( fabs(ele_tkSumPt_dr03) >= 0.5 ) return false;
    if ( fabs(ele_hcalDepth1TowerSumEt_dr03) >= 0.5 ) return false;

    // "conversion ID"
    if( ele_isConversion != 1 ) return false;

    // if passed all cuts...
    return true;
  }

  // HLT-like selection
  else if(asked=="HLT_Ele27") {
    for(int i=0 ; i<7 ; i++) {
      if( etaEle<1.479 ) {
	if( values[i] >= cutsHLT_Ele27_IdIso_EB[i] ) return false;
      }
      else {
	if( values[i] >= cutsHLT_Ele27_IdIso_EE[i] ) return false;
      }
    }
    return true;
  }
  
  // VBTF selection
  else if(asked=="WP95" || asked=="WP80" || asked=="WP60") {

    if(asked=="WP95") {
      for(int i=0 ; i<13 ; i++)
	cuts[i] = cuts95[i] ;
    }
    else if(asked=="WP80") {
      for(int i=0 ; i<13 ; i++)
	cuts[i] = cuts80[i] ;
    }
    else if(asked=="WP60") {
      for(int i=0 ; i<13 ; i++)
	cuts[i] = cuts60[i] ;
    }

    // missing hits
    if ( ele_expected_inner_hits > cuts[0] ) return false ;
	  
    if( fabs(etaEle)<1.479 ){ // barrel
      // idiEletification
      if ( fabs(ele_deltaphiin) >= cuts[1] )   return false ;
      if ( fabs(ele_deltaetain) >= cuts[2] ) return false ;
      if ( ele_he >= cuts[3] )                 return false ;
      if ( ele_sigmaietaieta >= 0.01 ) return false ;
	    
      // Relative isolation
      if ( fabs(trackIsoRel03)>= cuts[4] ) return false ;
      if ( fabs(hcalIsoRel03) >= cuts[5] ) return false ;
      if ( fabs(ecalIsoRel03) >= cuts[6] ) return false;

      // Combined Isolation
      
    } else { // endcap
      //return false; // remove endcap electrons
      // identification
      if ( fabs(ele_deltaphiin) >= cuts[7] )   return false ;
      if ( fabs(ele_deltaetain) >= cuts[8] ) return false ;
      if ( ele_he >= cuts[9] )                 return false ;
      if ( ele_sigmaietaieta >= 0.03 ) return false ;

      // Relative isolation
      if ( fabs(trackIsoRel03)>= cuts[10] ) return false ;
      if ( fabs(hcalIsoRel03) >= cuts[11] ) return false ;
      if ( fabs(ecalIsoRel03) >= cuts[12] ) return false;
    }

    //fiducial isolation
    //if ( fabs(eta) >= 2.5 || (fabs(eta)< 1.566 && fabs(eta)>1.4442)) return false ;
 
    // conversion rejection
    if(asked=="WP80")
      if ( fabs(ele_conv_dist) <= 0.02 && fabs(ele_conv_dcot) <=0.02 ) return false;
    if(asked=="WP60")
      if ( fabs(ele_conv_dist) <= 0.02 && fabs(ele_conv_dcot) <=0.02 ) return false;

    // if passed all cuts...
    return true;
  }

  else {
    return false;
  }

}

// int convertEtaToTT(double eta) {
// 
//   const int nTT=28;
//   float deta=0.087;
//   double etaTT[29]={0*deta, 1*deta, 2*deta, 3*deta, 4*deta, 5*deta, 
// 		   6*deta, 7*deta, 8*deta, 9*deta, 10*deta, 11*deta, 
// 		   12*deta, 13*deta, 14*deta, 15*deta, 16*deta, 17*deta, 
// 		   18*deta, 19*deta, 20*deta,1.83,1.93,2.043,2.172,2.322,2.50,2.65,3};
// 
//   if( fabs(eta)>3 ) return -999;
//   
// // // //   TMP
//   if(fabs(eta)> 1.479 ) return -999; not implemented for endcap yet.
// 
//   int sign = (eta < 0)? -1 : 1; 
// 
//   for (int i=0 ; i<nTT ; i++) {
//       if (fabs(eta)>=etaTT[i] && fabs(eta)<etaTT[i+1]) return sign*(i+1)
//   }
// }
int convertEtaToTT(double eta) {

  int sign = (eta < 0)? -1 : 1; //modif-alex

  if( fabs(eta)>3 ) return -999; //no ecal anymore
  //return (eta > 0)? (int)ceil(eta/0.087) : (int)floor(eta/0.087); 

  //Alex fixing geometry/ size of towers from ieta 21 are changing
  if( fabs(eta) < 1.740 )
    return (eta > 0)? (int)ceil(eta/0.087) : (int)floor(eta/0.087); 
  else{
    if( fabs(eta) >= 1.740 && fabs(eta) < 1.830) return sign*21; //size 0.09
    if( fabs(eta) >= 1.830 && fabs(eta) < 1.930) return sign*22; //size 0.1
    if( fabs(eta) >= 1.930 && fabs(eta) < 2.043) return sign*23; //size 0.113
    if( fabs(eta) >= 2.043 && fabs(eta) < 2.172) return sign*24; //size 0.129
    if( fabs(eta) >= 2.172 && fabs(eta) < 2.322) return sign*25; //size 0.15
    if( fabs(eta) >= 2.322 && fabs(eta) < 2.500) return sign*26; //size 0.178
    if( fabs(eta) >= 2.500 && fabs(eta) < 2.650) return sign*27; //size 0.15
    if( fabs(eta) >= 2.650 && fabs(eta) <= 3.000) return sign*28; //size 0.35
  }

  return -999;
}

int convertPhiToTT(double phi) {
  
  if( fabs(phi)<3.15) {
    int result = (int) ceil( phi/0.087 );
    if( phi >= 0 ) return result;
    else return (result+72);
  }

  return -999;
}



int convertEtaToRCT(double eta) {

  const int nReg=14;
  double etaRCT[nReg]={-2.586,-1.956,-1.566,-1.218,-0.8695,-0.5215,-0.174,0.174,0.5215,0.8695,1.218,1.566,1.956,2.586};

  if( eta<-3 || eta>3 ) return -999;

  double dist=10;
  int iFound=-999;

  for(int i=0 ; i<nReg ; i++) {
    if( fabs( eta - etaRCT[i] ) < dist ) {
      iFound = i;
      dist = fabs( eta - etaRCT[i] );
    }
  }

  return (4+iFound);
}

int convertPhiToRCT(double phi) {
  
  if( phi>-3.15 && phi<3.15 ) {
    int result = (int) arrondi( phi/0.349066 );
    if( phi >= 0 ) return result;
    else return (result+18);
  }

  return -999;
}

// ====================================================================================
int getGCTRegionPhi(int ttphi)
// ====================================================================================
{
  int gctphi=0;
  gctphi = (ttphi+1)/4;
  if(ttphi<=2) gctphi=0;
  if(ttphi>=71) gctphi=0;
  
  return gctphi;
}

// ====================================================================================
int getGCTRegionEta(int tteta)
// ====================================================================================
{
  int gcteta = 0;
  
  if(tteta>0) gcteta = (tteta-1)/4 + 11;
  else if(tteta<0) gcteta = (tteta+1)/4 + 10;
  
  return gcteta;
}

// ====================================================================================
void setMomentum(TLorentzVector & myvector , const LorentzVector & mom)
// ====================================================================================
{
  myvector.SetPx (mom.Px());
  myvector.SetPy (mom.Py());
  myvector.SetPz (mom.Pz());
  myvector.SetE (mom.E());
}

// ====================================================================================
bool IsConv (const reco::GsfElectron & eleRef)
// ====================================================================================
{
  bool isAmbiguous = true, isNotFromPixel = true;
  if (eleRef.ambiguousGsfTracksSize() == 0 ) {isAmbiguous = false;}

  //  int  mishits = eleRef.gsfTrack()->trackerExpectedHitsInner().numberOfHits();   
  int  mishits = eleRef.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  if (mishits == 0){isNotFromPixel = false;}
  
  bool is_conversion = false;
  if(isAmbiguous || isNotFromPixel) is_conversion = true;
  
  return is_conversion;
}
