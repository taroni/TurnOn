// -*- C++ -*-
//
// Package:   RunSelect
// Class:     RunSelect
//
//

#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <set>

// user include files
#include "EGamma/ECGelec/plugins/RunSelect.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;
using namespace std;

RunSelect::RunSelect(const edm::ParameterSet& iConfig)
{
  requireNoTimeScan = iConfig.getUntrackedParameter<bool>("requireNoTimeScan",true);
  requireCollidingBX = iConfig.getUntrackedParameter<bool>("requireCollidingBX",true);
  requireNoLumiScan = iConfig.getUntrackedParameter<bool>("requireNoLumiScan",true);
  requireLumiScan = iConfig.getUntrackedParameter<bool>("requireLumiScan",false);
  debug = iConfig.getUntrackedParameter<bool>("debug",false);
}

RunSelect::~RunSelect()
{
}

bool RunSelect::filter( edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  // list of runs in reprocessed minbias sample from dec. data
  // only interesting runs were reprocessed
  
/* 
 - 22/04/2010: lumi scans, time scans et runs douteux d'apres Clementine et la run/registery 
 - 10/05/2010: mise a jour d'apres run registry au 10/05, en demandant egalement L1 good 
 (HLT, Track, et Egam ne changent rien)
 - 10/05/2010: migration de la suppression des runs/LS correspondants depuis le cfg vers ce code .cc pour
 %G�%@viter les duplicated events  
*/
  
 int nruns = 71;
 int Runs[] = {
   132440, 132473, 132476, 132477, 132596, 132598, 132599, 132601, 132602, 132605, 
   132606, 132656, 132658, 132659, 132661, 132662, 132716, 132959, 132960, 132961, 
   132965, 132968, 133029, 133030, 133031, 133034, 133035, 133036, 133046, 133082, 
   133158, 133321, 133446, 133448, 133450, 133474, 133483, 133509, 133874, 133875,
   133876, 133877, 133881, 133885, 133927, 133928, 135059, 135149, 135175, 135445,
   135521, 135523, 135525, 135528, 135535, 135537, 135573, 135575, 135735,//, 135059
   // 09/06/2010
   136066, 136080, 136082, 136087, 136088, 136097, 136098, 137027, 137028,
   // 02/07/2010: Update of Run Selection
   138564, 138571, 138572, 138737, 138738, 138739, 138742, 138744, 138745, 138746,
   138747, 138750, 138751, 138919, 138920, 138921, 138923, 138924, 138937, 138939,
   139020, 139096, 139098, 139100, 139102, 139103

 };

  // Colliding BX
  vector<int> noBX ; 
  noBX.push_back(-1); // if info not present set to -1
  //vector<int> BX1 ; 
  //BX1.push_back(51); BX1.push_back(2724);
 
  vector<int> CollidingBX[] = {
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,     
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,//, noBX
    // 09/06/2010: Update of Run Selection
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,
    // 02/07/2010: Update of Run Selection
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,
    noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX, noBX,
    noBX, noBX, noBX, noBX, noBX, noBX
  };

  // lumi scans
  int LumiScan_min[] = {
    // -1 when no lumi scan
  -1,  -1,  -1,  -1,  -1,   1,  -1,1090,   1, 740, 
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
  -1,   1, 223,  -1,  -1,  -1,  -1,  -1,  -1, 165, 
 450,   1,  -1,  74,  -1,  -1, 305,  -1,  -1,  -1,
  -1,  -1,  -1,  -1,  -1,   1,  -1,  -1,  -1,  -1, 
  -1,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,//,  -1 
  // 09/06/2010: Update of Run Selection
  // WARNING: Lumi Scan should be checked...
  -1,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
  // 02/07/2010: Update of Run Selection
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
  -1,  -1,  -1,  -1,  -1,  -1
  };
  
  int LumiScan_max[] = {
    // -1 when no lumi scan
  -1,  -1,  -1,  -1,  -1,9999,  -1,9999,9999,9999, 
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
  -1,9999,9999,  -1,  -1,  -1,  -1,  -1,  -1, 230, 
 600,9999,  -1,9999,  -1,  -1, 370,  -1,  -1,  -1,
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, 
  -1,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,//,  -1 
  // 09/06/2010: Update of Run Selection
  // WARNING: Lumi Scan should be checked...
  -1,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
  // 02/07/2010: Update of Run Selection
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
  -1,  -1,  -1,  -1,  -1,  -1
  };
  

  bool TimeScan[] = {  /* not only ECAL and also including pixel bias scan */
    0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,//, 0
    // 09/06/2010: Update of Run Selection
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    // 02/07/2010: Update of Run Selection
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
  };
  
  //info on the event
  //int ievt = iEvent.id().event();
  int irun = iEvent.id().run();
  int ils = iEvent.luminosityBlock();
  int ibx = iEvent.bunchCrossing();
  
  // find run in the list
  int index=-1;
  for (int i=0; i<nruns; i++) {
    if (Runs[i]==irun) index = i;
  }
  if (index==-1) {
    //std::cout << "run " << irun << " not found in list of runs " << std::endl;
    return false;
  }

  //some debug cout
  if (debug) {
    // for (int index2=0;index2<nruns;index2++){
    int index2= index;
    cout<< "irun = "<<irun<<" index= "<<index2<< " run[index]= "<<Runs[index2]<<endl;
    cout<<" LS= "<<ils<<endl;
    cout<<" BX= "<<ibx<<" Colliding BX size="<< CollidingBX[index2].size()<<" value = ";
    for (unsigned int i=0;i<CollidingBX[index2].size();i++) cout<<" "<<CollidingBX[index2][i];
    cout<<endl;
    cout<<" lumiscan from "<<LumiScan_min[index2]<<" to "<<LumiScan_max[index2] <<endl;
    if (TimeScan[index2]) cout <<" TimeScan"<<endl; else cout <<" no TimeScan"<<endl;

    // } //debug all list
  }

  //  cout<<"run= "<<Runs[index]<<" LS= "<<ils<<endl;

  // CC 10/05 good run/LS by hand, version run registry du 10/05 en demandant
  // L1:GOOD, Pix:GOOD, Strip:GOOD, Ecal:GOOD, Hcal:GOOD, Es:GOOD, Track:GOOD, Egam:GOOD
  
//   if (irun==132440 && (ils<85 || (ils>138 && ils<141) || (ils>401))) return false;
//   if (irun==132473 && (ils>29)) return false;
//   if (irun==132476 && (ils<23 || (ils>28 && ils<54) || (ils>57))) return false;
//   if (irun==132477 && ((ils>5 && ils<34) || (ils>35 && ils<63)|| (ils>64 && ils<90)|| (ils>93 && ils<118)|| (ils>121 &&
//    ils<148)|| (ils>149 && ils<176)|| (ils>179 && ils<225)|| (ils>236 && ils<368)|| (ils>384 && ils<517) || ils>520)) return false;
//   if (irun==132596 && (ils<382 || (ils>383 && ils<447) || (ils>453))) return false;
//   if (irun==132598 && (ils<80 || (ils>82 && ils<174) || (ils>188))) return false;
//   if (irun==132599 && ((ils>379 && ils<381) || (ils>538))) return false;
//   if (irun==132601 && ((ils>207 && ils<209) || (ils>259 && ils<261)|| (ils>1131))) return false;
//   if (irun==132602 && ((ils>83))) return false;
//   if (irun==132605 && ((ils>444 && ils<446) || (ils>622 && ils<624)|| (ils>829 && ils<831)|| (ils>968 ))) return false;
//   if (irun==132606 && ((ils>37))) return false;
//   if (irun==132656 && ((ils>140))) return false;
//   if (irun==132658 && ((ils>177))) return false;
//   if (irun==132661 && ((ils>130))) return false;
//   if (irun==132662 && ((ils>130 && ils<132) || (ils>217))) return false;
//   if (irun==132716 && (ils<220 || (ils>591 && ils<593) || (ils>640))) return false;
//   if (irun==132959 && ((ils>276 && ils<278) || (ils>417))) return false;
//   if (irun==132659 && ((ils>84))) return false;
//   if (irun==132960 && ((ils>190))) return false;
//   if (irun==132961 && ((ils>427))) return false;
//   if (irun==132965 && ((ils>107))) return false;
//   if (irun==132968 && ((ils>173))) return false;
//   if (irun==133029 && (ils<101 || (ils>115 && ils<129) || (ils>350))) return false;
//   if (irun==133030 && ((ils>29))) return false;
//   if (irun==133031 && ((ils>18))) return false;
//   if (irun==133034 && ((ils<131) || (ils>325))) return false;
//   if (irun==133035 && ((ils>306))) return false;
//   if (irun==133036 && ((ils>225))) return false;
//   if (irun==133046 && ((ils>43 && ils<45) || (ils>323))) return false;
//   if (irun==133082 && ((ils>336 && ils<523) || (ils>592 && ils<595) ||(ils>608))) return false;
//   if (irun==133158 && ((ils<65) || (ils>786))) return false;
//   if (irun==133321 && ((ils>383))) return false;
//   if (irun==133446 && ((ils<105) || (ils>273))) return false;
//   if (irun==133448 && ((ils>516))) return false;
//   if (irun==133450 && ((ils>329 && ils<332) || (ils>658))) return false;
//   if (irun==133474 && ((ils>95 && ils<157) || (ils>189))) return false;
//   if (irun==133483 && ((ils<94) || (ils>159 && ils<161) ||(ils>591 && ils<652) ||(ils>658))) return false;

//   if (irun==133509 && ((ils<60) || (ils>75))) return false;
//   if (irun==133874 && ((ils<166) || (ils>814 && ils<817) ||(ils>875))) return false;
//   if (irun==133875 && ((ils>20 && ils<22) || (ils>49))) return false;
//   if (irun==133876 && ((ils>330))) return false;
//   if (irun==133877 && ((ils>1640 && ils<1643) || (ils>1997))) return false;
//   if (irun==133881 && ((ils>71 && ils<74) || (ils>223 && ils<225) ||(ils>562))) return false;
//   if (irun==133885 && ((ils>132 && ils<134) || (ils>728))) return false;
//   if (irun==133927 && ((ils>57))) return false;
//   if (irun==133928 && ((ils>645))) return false;

//   if (irun==135059 && ((ils<59) || (ils>67))) return false;


  // SB 19/05/10 good run/LS by hand, version run registry du 19/05 en demandant
  // QFLAGS=L1t:GOOD,Pix:GOOD,Strip:GOOD,Ecal:GOOD,Hcal:GOOD,Es:GOOD,Track:GOOD,Egam:GOOD
  // DCS=Bpix,Fpix,Tibtid,TecM,TecP,Tob,Ebminus,Ebplus,EeMinus,EePlus,EsMinus,EsPlus,HbheA,HbheB,HbheC,H0,Hf
  // energy = 3500 GeV

//   if (irun==132440 && (ils<85 || (ils>138 && ils<141) || (ils>401))) return false;
//   if (irun==132473 && (ils>29)) return false;
//   if (irun==132476 && (ils<23 || (ils>28 && ils<54) || (ils>57))) return false;
//   if (irun==132477 && ((ils>5 && ils<34) || (ils>35 && ils<63)|| (ils>64 && ils<90)|| (ils>93 && ils<118)|| (ils>121 && ils<148)|| (ils>149 && ils<176)|| (ils>179 && ils<225)|| (ils>236 && ils<368)|| (ils>384 && ils<517) || ils>520)) return false;
//   if (irun==132596 && (ils<382 || (ils>383 && ils<447) || (ils>453))) return false;
//   if (irun==132598 && (ils<80 || (ils>82 && ils<174) || (ils>188))) return false;
//   if (irun==132599 && ((ils>379 && ils<381) || (ils>538))) return false;
//   if (irun==132601 && ((ils>207 && ils<209) || (ils>259 && ils<261)|| (ils>1131))) return false;
//   if (irun==132602 && ((ils>83))) return false;
//   if (irun==132605 && ((ils>444 && ils<446) || (ils>622 && ils<624)|| (ils>829 && ils<831)|| (ils>968 ))) return false;
//   if (irun==132606 && ((ils>37))) return false;
//   if (irun==132656 && ((ils>140))) return false;
//   if (irun==132658 && ((ils>177))) return false;
//   if (irun==132661 && ((ils>130))) return false;
//   if (irun==132662 && ((ils>130 && ils<132) || (ils>217))) return false;
//   if (irun==132716 && (ils<220 || (ils>591 && ils<593) || (ils>640))) return false;
//   if (irun==132959 && ((ils>276 && ils<278) || (ils>417))) return false;
//   if (irun==132960 && ((ils>190))) return false;
//   if (irun==132961 && ((ils>427))) return false;
//   if (irun==132965 && ((ils>107))) return false;
//   if (irun==132968 && ((ils>173))) return false;
//   if (irun==133029 && (ils<101 || (ils>115 && ils<129) || (ils>350))) return false;
//   if (irun==133030 && ((ils>29))) return false;
//   if (irun==133031 && ((ils>18))) return false;
//   if (irun==133034 && ((ils<131) || (ils>325))) return false;
//   if (irun==133035 && ((ils>306))) return false;
//   if (irun==133036 && ((ils>225))) return false;
//   if (irun==133046 && ((ils>43 && ils<45) || (ils>323))) return false;
//   if (irun==133082 && ((ils>336 && ils<523) || (ils>592 && ils<595) ||(ils>608))) return false;
//   if (irun==133158 && ((ils<65) || (ils>786))) return false;
//   if (irun==133321 && ((ils>383))) return false;
//   if (irun==133446 && ((ils<105) || (ils>273))) return false;
//   if (irun==133448 && ((ils>516))) return false;
//   if (irun==133450 && ((ils>329 && ils<332) || (ils>658))) return false;
//   if (irun==133474 && ((ils>95 && ils<157) || (ils>189))) return false;
//   if (irun==133483 && ((ils<94) || (ils>159 && ils<161) ||(ils>591 && ils<652) ||(ils>658))) return false;
// 
//   if (irun==133509 && ((ils<60) || (ils>75))) return false;
//   if (irun==133874 && ((ils<166) || (ils>814 && ils<817) ||(ils>875))) return false;
//   if (irun==133875 && ((ils>20 && ils<22) || (ils>49))) return false;
//   if (irun==133876 && ((ils>330))) return false;
//   if (irun==133877 && ((ils>1640 && ils<1643) || (ils>1997))) return false;
//   if (irun==133881 && ((ils>71 && ils<74) || (ils>223 && ils<225) ||(ils>562))) return false;
//   if (irun==133885 && ((ils>132 && ils<134) || (ils>728))) return false;
//   if (irun==133927 && ((ils>57))) return false;
//   if (irun==133928 && ((ils>645))) return false;

  //  if (irun==135059 && ((ils<59) || (ils>67))) return false;
  
  // CC 28/05/10 good run/LS by hand, version run registry du 28/05 jusqu'a 135735 en demandant
  // QFLAGS=L1t:GOOD,Pix:GOOD,Strip:GOOD,Ecal:GOOD,Hcal:GOOD,Es:GOOD,Track:GOOD,Egam:GOOD
  // DCS=Bpix,Fpix,Tibtid,TecM,TecP,Tob,Ebminus,Ebplus,EeMinus,EePlus,EsMinus,EsPlus,HbheA,HbheB,HbheC,H0,Hf
  // energy = 3500 GeV

  if (irun==132440 && (ils<85 || (ils>138 && ils<141) || (ils>401))) return false;
  if (irun==132473 && (ils>29)) return false;
  if (irun==132476 && (ils<23 || (ils>28 && ils<54) || (ils>57))) return false;
  if (irun==132477 && ((ils>5 && ils<34) || (ils>35 && ils<63)|| (ils>64 && ils<90)|| (ils>93 && ils<118)|| (ils>121 && ils<148)|| (ils>149 && ils<176)|| (ils>179 && ils<225)|| (ils>236 && ils<368)|| (ils>384 && ils<517) || ils>520)) return false;
  if (irun==132596 && (ils<382 || (ils>383 && ils<447) || (ils>453))) return false;
  if (irun==132598 && (ils<80 || (ils>82 && ils<174) || (ils>188))) return false;
  if (irun==132599 && ((ils>379 && ils<381) || (ils>538))) return false;
  if (irun==132601 && ((ils>207 && ils<209) || (ils>259 && ils<261)|| (ils>1131))) return false;
  if (irun==132602 && ((ils>83))) return false;
  if (irun==132605 && ((ils>444 && ils<446) || (ils>622 && ils<624)|| (ils>829 && ils<831)|| (ils>968 ))) return false;
  if (irun==132606 && ((ils>37))) return false;
  if (irun==132656 && ((ils>140))) return false;
  if (irun==132658 && ((ils>177))) return false;
  if (irun==132661 && ((ils>130))) return false;
  if (irun==132662 && ((ils>130 && ils<132) || (ils>217))) return false;
  if (irun==132716 && (ils<220 || (ils>591 && ils<593) || (ils>640))) return false;
  if (irun==132959 && ((ils>276 && ils<278) || (ils>417))) return false;
  if (irun==132960 && ((ils>190))) return false;
  if (irun==132961 && ((ils>427))) return false;
  if (irun==132965 && ((ils>107))) return false;
  if (irun==132968 && ((ils>173))) return false;

  if (irun==133029 && (ils<101 || (ils>115 && ils<129) || (ils>350))) return false;
  if (irun==133030 && ((ils>29))) return false;
  if (irun==133031 && ((ils>18))) return false;
  if (irun==133034 && ((ils<131) || (ils>325))) return false;
  if (irun==133035 && ((ils>306))) return false;
  if (irun==133036 && ((ils>225))) return false;
  if (irun==133046 && ((ils>43 && ils<45) || (ils>323))) return false;
  if (irun==133082 && ((ils>336 && ils<523) || (ils>592 && ils<595) ||(ils>608))) return false;
  if (irun==133158 && ((ils<65) || (ils>786))) return false;
  if (irun==133321 && ((ils>383))) return false;
  if (irun==133446 && ((ils<105) || (ils>273))) return false;
  if (irun==133448 && ((ils>516))) return false;
  if (irun==133450 && ((ils>329 && ils<332) || (ils>658))) return false;
  if (irun==133474 && ((ils>95 && ils<157) || (ils>189))) return false;
  if (irun==133483 && ((ils<94) || (ils>159 && ils<161) ||(ils>591 && ils<652) ||(ils>658))) return false;
  if (irun==133509 && ((ils<60) || (ils>75))) return false;
  if (irun==133874 && ((ils<166) || (ils>814 && ils<817) ||(ils>875))) return false;
  if (irun==133875 && ((ils>20 && ils<22) || (ils>49))) return false;
  if (irun==133876 && ((ils>330))) return false;
  if (irun==133877 && ((ils>1640 && ils<1643) || (ils>1997))) return false;
  if (irun==133881 && ((ils>71 && ils<74) || (ils>223 && ils<225) ||(ils>562))) return false;
  if (irun==133885 && ((ils>132 && ils<134) || (ils>728))) return false;
  if (irun==133927 && ((ils>57))) return false;
  if (irun==133928 && ((ils>645))) return false;
  
  if (irun==135059 && ((ils<59) || (ils>67))) return false;
  if (irun==135149 && ((ils<294) || (ils>337 && ils<339) ||(ils>754 && ils<756) ||(ils>932 && ils<934)||(ils>1808 && ils<1811) || ils>3382)) return false;
  if (irun==135175 && ((ils<55) || (ils>561 && ils<563) || (ils>790 && ils<792) ||(ils>1082))) return false;
  if (irun==135445 && ((ils<997) || (ils>1329 && ils<1332) ||(ils>1827))) return false;
  if (irun==135521 && ((ils<60) || (ils>440 && ils<442) ||(ils>524))) return false;
  if (irun==135523 && (((ils>124 && ils<126) ||(ils>225)))) return false;
  if (irun==135525 && (((ils>3 && ils<6) ||(ils>457)))) return false;
  if (irun==135528 && ((ils<98) || (ils>147 && ils<149) || (ils>813 && ils<816)||(ils>924 && ils<926) || ils>1436)) return false;
  if (irun==135535 && (((ils<75) ||(ils>167 && ils<169) || (ils>246)))) return false;
  if (irun==135537 && ((ils<39) || (ils>69))) return false;
  if (irun==135573 && ((ils<102) || (ils>163))) return false;
  if (irun==135575 && ((ils<2) || (ils>210 && ils<213) ||(ils>381 && ils<384) ||(ils>638 && ils<645)||(ils>11561 && ils<1163) || ils>1253)) return false;
  if (irun==135735 && ((ils<31) || (ils>333))) return false;

  // 09/06/2010: Update of Run Selection
  if (irun==136066 && ((ils<181) || (ils>297 && ils<299) || (ils>348 && ils<350) || (ils>529 &&  ils<532) || (ils>595 && ils<597) || (ils>1184))) return false;
  if (irun==136080 && ((ils<249) || (ils>262))) return false;
  if (irun==136082 && (( (ils>422 && ils<477) || (ils>506)))) return false;
  if (irun==136087 && ((ils<250) || (ils>315 && ils<335) ||  ils>354)) return false;
  if (irun==136088 && ((ils>262))) return false;
  if (irun==136097 && ((ils>91)))  return false;
  if (irun==136098 && ((ils>25))) return false;
  if (irun==137027 && ((ils<98) && (ils>200))) return false;
  if (irun==137028 && ((ils>484))) return false;

  // 02/07/2010: Update of Run Selection
  if (irun==138564 && (ils>19)) return false;
  if (irun==138571 && (ils>23)) return false;
  if (irun==138572 && (ils>238)) return false;
  if (irun==138737 && (ils>88)) return false;
  if (irun==138738 && (ils>94)) return false;
  if (irun==138739 && (ils>13)) return false;
  if (irun==138742 && (((ils>20)&&(ils<22))||(ils>54))) return false;
  if (irun==138744 && (ils>39)) return false;
  if (irun==138745 && (ils>17)) return false;
  if (irun==138746 && (ils>145)) return false;
  if (irun==138747 && ((ils>71 &&ils<73)||(ils>131))) return false;
  if (irun==138750 && ((ils>46 && ils<49) || (ils>208 && ils<210) || (ils>623 && ils<626) || (ils>715 && ils<717) || (ils>726))) return false;
  if (irun==138751 && ((ils>110 && ils<112)||(ils>147))) return false;
  if (irun==138919 && ((ils<62) || (ils>168))) return false;
  if (irun==138920 && (ils>74)) return false;
  if (irun==138921 && ((ils>181 && ils<183)||(ils>203))) return false;
  if (irun==138923 && (ils>17)) return false;
  if (irun==138924 && (ils>86)) return false;
  if (irun==138937 && ((ils>29 && ils<31)||(ils>47))) return false;
  if (irun==138939 && (ils>36)) return false;
  if (irun==139020 && ((ils<227) || (ils>316 && ils<319) || (ils>617))) return false;
  if (irun==139096 && ((ils<193) || (ils>280))) return false;
  if (irun==139098 && (ils>201)) return false;
  if (irun==139100 && (ils>315)) return false;
  if (irun==139102 && ((ils>12 && ils<14)||(ils>55))) return false;
  if (irun==139103 && ((ils<7) || (ils>449))) return false;

  if (requireNoTimeScan  && TimeScan[index]){
    if (debug) {cout <<" no time scan requierement failed"<<endl; cout<<endl; }
    return false;
  }
  
  if (requireCollidingBX) { 
    bool collBX = false;
    for (unsigned int i=0;i<CollidingBX[index].size();i++) 
      if (ibx==CollidingBX[index][i]) collBX=true; 
    if (!collBX) {
      if (debug) {cout <<"Colliding BX requierement failed"<<endl; cout<<endl;}
      return false;
    }
  }
  if (requireNoLumiScan  && (ils>=LumiScan_min[index] && ils<=LumiScan_max[index])  ) {
    if (debug) {cout <<" no lumi scan requierement failed"<<endl; cout<<endl; }
    return false;
  }
  if (requireLumiScan && !(ils>=LumiScan_min[index] && ils<=LumiScan_max[index])){
    if (debug) {cout <<" LumiScan positive requirement failed"<<endl; cout<<endl;}
    return false;
  }
  
  if (debug) {
    cout<<"all requirements OK"<<endl;
    cout<<endl;
  }

  return true;

}
