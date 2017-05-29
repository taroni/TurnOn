///g++ -o    selectPairs selectPairs.cc   `root-config --cflags --libs` -L $ROOFITSYS/lib -lRooFit -lRooFitCore
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFrame.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include "baseFuncNad.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <glob.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>

vector<string> globVector(const string& pattern){
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> files;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    files.push_back(string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return files;
}


void process(TChain* myChain, string dirOut, int iEff, int iEff_M, TString tagCuts, double tagEt, TString probeCuts, double probeEt, double massCut1, double massCut2,int nEntries, ofstream& outlog, TString tag_hlt , bool checkCharge, bool data2011,  TString eg_menu, bool debug, string fileName, bool isMC, bool isStage2);

void selectPairs(vector<string> files, int nEntries, int iEff, int iEff_M, TString eg_menu, string dirOut, string dirIn, TString tagCuts, double tagEt, TString probeCuts, double probeEt, double massCut1, double massCut2, TString tag_hlt, TString nameChain, bool checkCharge, bool data2011, bool debug, bool isMC, bool isStage2);



using namespace std;
string myjobid;
//-------
int main(int argc, char**argv){
  myjobid=argv[1];
  stringstream inputDir, outputDir;
  inputDir.str("");
  outputDir.str(""); 
  inputDir <<"/eos/cms/store/user/taroni/DoubleEG/turnOnNtuple_EEFG_20ADC_simTP/makePairs";
  outputDir<<"selectPairsDir";

  vector<string> vinput;
  string line;
  ifstream myfile (myjobid);
  if(myfile.is_open())
    {
      while (getline(myfile, line))
	{
	  stringstream myline; 
	  myline.str("");
	  myline << "root://eoscms/"<< line; 
	  //myline <<  line; 
	    
	  vinput.push_back(myline.str().c_str()); 
	  
	}
      myfile.close();
    }
  for (int i =0; i< vinput.size(); i++){
    cout<< "my input file " << vinput[i]<< endl; 
  }
  
  cout << outputDir.str()<< endl;
  selectPairs(vinput, -1,1,1 , "L1_EG_MENU_2011",
	      outputDir.str(),
	      inputDir.str(),
	      "Medium", 30, "Loose", 5, //RunII ID
	      //"WP80", 20, "WP95", 5, 
	      //	      60, 120, "HLT_Ele23_Ele12", "ElePairs",
	      75, 105, "HLT_Ele23_Ele12", "ElePairs",
	      true, true, false,false,true);
  
  return 0;
}
//------





void selectPairs(vector<string>filenames, int nEntries=-1, int iEff=1, int iEff_M=1, TString eg_menu="L1_EG_MENU_2011",
		 string dirOut="selectPairsDir",
		 string dirIn="root://eoscms//eos/cms/store/user/taroni/DoubleEG/turnOnNtuple_EEFG_20ADC_simTP/makePairs",
		 //TString tagCuts="WP80", double tagEt=20, TString probeCuts="WP95", double probeEt=5,
		 TString tagCuts="Medium", double tagEt=20, TString probeCuts="Veto", double probeEt=5,  
		 double massCut1=60, double massCut2=120, TString tag_hlt="noHLT", TString nameChain="ElePairs",
		 bool checkCharge=true, bool data2011=true, bool debug=false, bool isMC=false, bool isStage2=true)
{
  // Output Log
  stringstream logfile; 
  logfile.str("");
  logfile << dirOut << "/log_selectPairs.txt";
  ofstream outlog(logfile.str().c_str(),ios::out);
  // Input trees
  cout << "number of files: " << filenames.size()<< endl; 
  for (int i =0 ; i< filenames.size() ; i++){
    cout << filenames[i] << endl;
  }
  for (int ifile=0; ifile<filenames.size() ; ifile++){
    TChain * myChain = new TChain(nameChain);
    stringstream myfile; 
    myfile.str("");
    
    string file = filenames[ifile];
    cout << " file name " << file << endl ; 

    myChain->Add(file.c_str()); 
    
    if(debug) cout << myChain->GetEntries() << " entries" << endl;
    cout << "myChain name " << myChain->GetName() << endl;
    // Process the tree
    if(debug) cout << "process the tree" << endl;
    
    process(myChain, dirOut, iEff, iEff_M, tagCuts, tagEt, probeCuts, probeEt, massCut1, massCut2,
	    nEntries, outlog, tag_hlt, checkCharge, data2011, eg_menu, debug, filenames[ifile], isMC, isStage2);
    
    cout << "DONE, file " << ifile << endl;
  }
}

void process(TChain* myChain, string dirOut, int iEff, int iEff_M,
	     TString tagCuts, double tagEt, TString probeCuts, double probeEt, double massCut1, double massCut2,
	     int nEntries, ofstream& outlog, TString tag_hlt , bool checkCharge, bool data2011, 
	     TString eg_menu, bool debug, string fileName, bool isMC, bool isStage2)  {

  stringstream outlogstr; 
  outlogstr.str(""); 
  outlogstr << dirOut << "/elepairs.json";
  ofstream out_elepairs_json(outlogstr.str().c_str(),ios::out);
  outlogstr.str(""); 
  outlogstr << dirOut << "/notToZero.txt";  
  ofstream outNotToZero(outlogstr.str().c_str(), ios::out);
  outlogstr.str(""); 
  outlogstr << dirOut << "/probes_trigEG15.txt";    
  ofstream out_probes_trigEG15(outlogstr.str().c_str(), ios::out);
  outlogstr.str(""); 
  outlogstr << dirOut << "/probes_notrigEG15.txt";  
  ofstream out_probes_notrigEG15(outlogstr.str().c_str(),ios::out);
  outlogstr.str(""); 
  outlogstr << dirOut << "/probes_notrigEG15_RLE.txt";  
  ofstream out_notrigEG15_RLE(outlogstr.str().c_str(),ios::out);
  outlogstr.str(""); 
  outlogstr << dirOut << "/probes_notrigEG15_R.txt";  
  ofstream out_notrigEG15_R(outlogstr.str().c_str(),ios::out);
  outlogstr.str(""); 
  outlogstr << dirOut << "/inefficiency.json";
  ofstream out_custom_json(outlogstr.str().c_str(),ios::out);

  out_custom_json << "{" << endl;
  out_elepairs_json << "{" << endl;
  // OUTPUT TREES //
  if(debug) cout << "declare output trees" << endl;
  // L1 EG menu
  vector<int> menu;
  if(eg_menu=="L1_EG_MENU_2010") {
    int EG_2010[10] = {0,2,4,10,16,20,24,30,40,60};
    for(int i=0 ; i<10 ; i++)
      menu.push_back(EG_2010[i]);
  }
  else if(eg_menu=="L1_EG_MENU_2011") {
    for(int i=0; i<41 ; i++)
      menu.push_back(5*i);
  }
  const int nEG=menu.size();

  TString EGval[nEG];
  ostringstream ossi;
  for(int iEG=0;iEG<nEG;iEG++) {
    ossi.str("");
    ossi << menu[iEG] ;
    EGval[iEG] = ossi.str();
  }
  ////

  const int nECAL = 3; // EB / EE / EBEE
  const int nColl = 10; // _N / _M / _S2 / _S2E /_S2EI / _S1I / _S1E / _S1EI / _S2I / _S2test

  TString name_ecal[nECAL] = {"_EB","_EE", "_EBEE"};
  TString name_coll[nColl] = {"_N","_M", "_S2", "_S2E", "_S2EI", "_S1I", "_S1E", "_S1EI", "_S2I", "_S2test"};  

  TTree* treenew[nECAL][nColl];
  float sc_et[nECAL][nColl];
  float ele_pt[nECAL][nColl];
  float sc_eta[nECAL][nColl];
  float sc_dr[nECAL][nColl];
  int n_vtx;
  double PU_rhoCorr;
  int PU_N;
  int sc_L1_nTT[nECAL][nColl];


  Int_t l1_bin[nECAL][nColl][nEG];

  int nHLT;
  vector<TString> HLT_names;
  TString HLT_names_2011[10] = {"unbias","EG2","EG5","EG8","EG12","EG12old","Ele15","Ele17","Ele8calIso","Ele8tkIso"} ;
  TString HLT_names_2010[4] = {"unbias","EG2","EG5","EG8"} ;

  HLT_names.clear();
  if(data2011) {
    nHLT = 10;
    for(int i=0 ; i<nHLT ; i++) 
      HLT_names.push_back(HLT_names_2011[i]);
  }
  else {
    nHLT = 4;
    for(int i=0 ; i<nHLT ; i++) 
      HLT_names.push_back(HLT_names_2010[i]);
  }



  TString name, name_scet, name_pt, name_sceta, name_scdr, name_vtx_N, name_PU, name_PU_MC, name_sc_L1_nTT, name_l1bin;
  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {

      name = "treenew" + name_ecal[iECAL] + name_coll[iColl];
      treenew[iECAL][iColl] = new TTree(name,name);

      name_scet = "sc_et"+name_ecal[iECAL]+name_coll[iColl];
      name_pt = "ele_pt"+name_ecal[iECAL]+name_coll[iColl];
      name_sceta = "sc_eta"+name_ecal[iECAL]+name_coll[iColl];
      name_scdr = "sc_dr"+name_ecal[iECAL]+name_coll[iColl]; 
      name_vtx_N = "n_vtx";
      name_PU = "PU_rhoCorr";
      name_PU_MC = "PU_N";
      name_sc_L1_nTT = "sc_L1_nTT"+name_ecal[iECAL]+name_coll[iColl];


      treenew[iECAL][iColl]->Branch(name_scet, &sc_et[iECAL][iColl], name_scet+"/F");
      treenew[iECAL][iColl]->Branch(name_pt, &ele_pt[iECAL][iColl], name_pt+"/F");
      treenew[iECAL][iColl]->Branch(name_scdr, &sc_dr[iECAL][iColl], name_scdr+"/F");	
      treenew[iECAL][iColl]->Branch(name_sceta, &sc_eta[iECAL][iColl], name_sceta+"/F");
      treenew[iECAL][iColl]->Branch(name_vtx_N, &n_vtx, name_vtx_N+"/I");
      treenew[iECAL][iColl]->Branch(name_PU, &PU_rhoCorr, name_PU+"/D");
      treenew[iECAL][iColl]->Branch(name_PU_MC, &PU_N, name_PU_MC+"/I");      
      treenew[iECAL][iColl]->Branch(name_sc_L1_nTT, &sc_L1_nTT[iECAL][iColl], name_sc_L1_nTT+"/I");



      for(int iEG=0 ; iEG<nEG ; iEG++) {
	name_l1bin = "l1_"+EGval[iEG]+name_ecal[iECAL]+name_coll[iColl];
	treenew[iECAL][iColl]->Branch(name_l1bin , &l1_bin[iECAL][iColl][iEG] , name_l1bin+"/I");	
      }
    }
  }

  if(debug) cout << "get the input tree" << endl;

  // INPUT TREE //
  int nEvent, nRun, nLumi ;
  double _PU_rhoCorr;
  int _PU_N;
  int vtx_N;
  int trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8, 
    trig_isL1SingleEG12, trig_isL1SingleEG12_old, trig_is_HLT_Ele15, trig_is_HLT_Ele17, trig_is_HLT_Ele8_calIso, trig_is_HLT_Ele8_tkIso ;   
  
  trig_isUnbiased =  trig_isL1SingleEG2 =  trig_isL1SingleEG5 =  trig_isL1SingleEG8 =  trig_isL1SingleEG12 =  trig_isL1SingleEG12_old =  trig_is_HLT_Ele15 =  trig_is_HLT_Ele17 =  trig_is_HLT_Ele8_calIso =  trig_is_HLT_Ele8_tkIso = 0;
  
  // TP info
  int trig_tower_N,trig_tower_ieta[4032],trig_tower_iphi[4032],trig_tower_adc[4032],trig_tower_sFGVB[4032];
  int trig_tower_N_modif,trig_tower_ieta_modif[4032],trig_tower_iphi_modif[4032],trig_tower_adc_modif[4032],trig_tower_sFGVB_modif[4032]; 
  int trig_tower_N_emul,trig_tower_ieta_emul[4032],trig_tower_iphi_emul[4032],trig_tower_adc_emul[4032][5],trig_tower_sFGVB_emul[4032][5]; 

  trig_tower_N = trig_tower_N_modif = trig_tower_N_emul = 0;
  for(int iT=0;iT<4032;iT++) {
    trig_tower_ieta[iT]=trig_tower_iphi[iT]=trig_tower_adc[iT]=trig_tower_sFGVB[iT]=-999;
    trig_tower_ieta_modif[iT]=trig_tower_iphi_modif[iT]=trig_tower_adc_modif[iT]=trig_tower_sFGVB_modif[iT]=-999;
    trig_tower_ieta_emul[iT]=trig_tower_iphi_emul[iT]=-999;
    for(int iEm=0;iEm<5;iEm++)
      trig_tower_adc_emul[iT][iEm]=trig_tower_sFGVB_emul[iT][iEm]=-999;
  }

  // L1 candidates info
  int trig_L1emIso_N, trig_L1emNonIso_N, trig_L1emIso_N_M, trig_L1emNonIso_N_M;
  int trig_L1emIso_ieta[4],trig_L1emIso_energy[4], trig_L1emIso_iphi[4], trig_L1emIso_rank[4];
  int trig_L1emNonIso_ieta[4],trig_L1emNonIso_energy[4], trig_L1emNonIso_iphi[4], trig_L1emNonIso_rank[4];
  int trig_L1emIso_ieta_M[4], trig_L1emIso_iphi_M[4], trig_L1emIso_rank_M[4];
  int trig_L1emNonIso_ieta_M[4], trig_L1emNonIso_iphi_M[4], trig_L1emNonIso_rank_M[4];

  // L1 prefiring
  int trig_preL1emIso_N; 
  int trig_preL1emNonIso_N;
  int trig_preL1emIso_ieta[4], trig_preL1emIso_iphi[4], trig_preL1emIso_rank[4]; 
  int trig_preL1emNonIso_ieta[4], trig_preL1emNonIso_iphi[4],trig_preL1emNonIso_rank[4];
  // L1 postfiring
  int trig_postL1emIso_N; 
  int trig_postL1emNonIso_N;
  int trig_postL1emIso_ieta[4], trig_postL1emIso_iphi[4], trig_postL1emIso_rank[4]; 
  int trig_postL1emNonIso_ieta[4], trig_postL1emNonIso_iphi[4],trig_postL1emNonIso_rank[4];
  
  // Masking
  int trig_nMaskedRCT, trig_nMaskedCh;
  int trig_iMaskedRCTeta[100], trig_iMaskedRCTphi[100], trig_iMaskedRCTcrate[100], trig_iMaskedTTeta[100], trig_iMaskedTTphi[100];


  // Pairs info
  double pair_M;
  double pair_eta[2],pair_sclPhi[2], pair_sclEta[2], pair_sclEt[2], pair_phi[2], pair_pT[2], pair_eT[2], pair_E[2];
  int pair_Run2IdLevel[2], pair_HLT_Ele27_cut[2], pair_fidu[2], pair_charge[2], pair_RCTeta[2], pair_RCTphi[2], 
    //int pair_cuts[2], pair_HLT_Ele27_cut[2], pair_fidu[2], pair_charge[2], pair_RCTeta[2], pair_RCTphi[2], 
    pair_L1iso[2], pair_L1noniso[2], pair_L1iso_M[2], pair_L1noniso_M[2];
  int pair_L1Stage2[2], pair_L1Stage2_emul[2];
  int pair_L1Stage2_isoflag[2], pair_L1Stage2_emul_isoflag[2];
  int pair_L1Stage2_emul_nTT[2];

  int pair_L1Stage1_emul[2];
  int pair_L1Stage1_emul_isoflag[2];

  int pair_RCTetaVect[2][10], pair_RCTphiVect[2][10], 
    pair_L1isoVect[2][10], pair_L1nonisoVect[2][10], 
    pair_L1isoVect_M[2][10], pair_L1nonisoVect_M[2][10],
    pair_L1Stage1_emul_isoVect[2][10], pair_L1Stage1_emul_nonisoVect[2][10];
  int pair_TTetaVect[2][50], pair_TTphiVect[2][50];
  double pair_TTetVect[2][50];
  double pair_RCTetVect[2][10];

  // Global
  myChain->SetBranchAddress("nEvent",&nEvent);
  myChain->SetBranchAddress("nRun",&nRun);
  myChain->SetBranchAddress("nLumi",&nLumi);

  myChain->SetBranchAddress("PU_rhoCorr",&_PU_rhoCorr);
  myChain->SetBranchAddress("PU_N",&_PU_N);

  myChain->SetBranchAddress("vtx_N",&vtx_N);

  // L1 candidates
  myChain->SetBranchAddress("trig_L1emIso_N", &trig_L1emIso_N);
  myChain->SetBranchAddress("trig_L1emIso_ieta", &trig_L1emIso_ieta);
  myChain->SetBranchAddress("trig_L1emIso_energy", &trig_L1emIso_energy);
  myChain->SetBranchAddress("trig_L1emIso_iphi", &trig_L1emIso_iphi);
  myChain->SetBranchAddress("trig_L1emIso_rank", &trig_L1emIso_rank);

  myChain->SetBranchAddress("trig_L1emNonIso_energy", &trig_L1emNonIso_energy);
  myChain->SetBranchAddress("trig_L1emNonIso_N", &trig_L1emNonIso_N);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta", &trig_L1emNonIso_ieta);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi", &trig_L1emNonIso_iphi); 
  myChain->SetBranchAddress("trig_L1emNonIso_rank", &trig_L1emNonIso_rank);

  myChain->SetBranchAddress("trig_L1emIso_N_M", &trig_L1emIso_N_M);
  myChain->SetBranchAddress("trig_L1emIso_ieta_M", &trig_L1emIso_ieta_M);
  myChain->SetBranchAddress("trig_L1emIso_iphi_M", &trig_L1emIso_iphi_M);
  myChain->SetBranchAddress("trig_L1emIso_rank_M", &trig_L1emIso_rank_M);

  myChain->SetBranchAddress("trig_L1emNonIso_N_M", &trig_L1emNonIso_N_M);
  myChain->SetBranchAddress("trig_L1emNonIso_ieta_M", &trig_L1emNonIso_ieta_M);
  myChain->SetBranchAddress("trig_L1emNonIso_iphi_M", &trig_L1emNonIso_iphi_M);
  myChain->SetBranchAddress("trig_L1emNonIso_rank_M", &trig_L1emNonIso_rank_M);

  // pre/post-firing L1 candidates
  myChain->SetBranchAddress("trig_preL1emIso_N",     &trig_preL1emIso_N);
  myChain->SetBranchAddress("trig_preL1emIso_ieta",  &trig_preL1emIso_ieta);
  myChain->SetBranchAddress("trig_preL1emIso_iphi",  &trig_preL1emIso_iphi);
  myChain->SetBranchAddress("trig_preL1emIso_rank",  &trig_preL1emIso_rank);
  //
  myChain->SetBranchAddress("trig_preL1emNonIso_N",     &trig_preL1emNonIso_N);
  myChain->SetBranchAddress("trig_preL1emNonIso_ieta",  &trig_preL1emNonIso_ieta);
  myChain->SetBranchAddress("trig_preL1emNonIso_iphi",  &trig_preL1emNonIso_iphi);
  myChain->SetBranchAddress("trig_preL1emNonIso_rank",  &trig_preL1emNonIso_rank);
  //
  myChain->SetBranchAddress("trig_postL1emIso_N",     &trig_postL1emIso_N);
  myChain->SetBranchAddress("trig_postL1emIso_ieta",  &trig_postL1emIso_ieta);
  myChain->SetBranchAddress("trig_postL1emIso_iphi",  &trig_postL1emIso_iphi);
  myChain->SetBranchAddress("trig_postL1emIso_rank",  &trig_postL1emIso_rank);
  //
  myChain->SetBranchAddress("trig_postL1emNonIso_N",     &trig_postL1emNonIso_N);
  myChain->SetBranchAddress("trig_postL1emNonIso_ieta",  &trig_postL1emNonIso_ieta);
  myChain->SetBranchAddress("trig_postL1emNonIso_iphi",  &trig_postL1emNonIso_iphi);
  myChain->SetBranchAddress("trig_postL1emNonIso_rank",  &trig_postL1emNonIso_rank);


  // Trigger Towers
 
  // Masking
  myChain->SetBranchAddress("trig_nMaskedRCT",      &trig_nMaskedRCT);      
  myChain->SetBranchAddress("trig_iMaskedRCTeta",   &trig_iMaskedRCTeta);                                          
  myChain->SetBranchAddress("trig_iMaskedRCTcrate", &trig_iMaskedRCTcrate);
  myChain->SetBranchAddress("trig_iMaskedRCTphi",   &trig_iMaskedRCTphi);
  myChain->SetBranchAddress("trig_nMaskedCh",       &trig_nMaskedCh);    
  myChain->SetBranchAddress("trig_iMaskedTTeta",    &trig_iMaskedTTeta);   
  myChain->SetBranchAddress("trig_iMaskedTTphi",    &trig_iMaskedTTphi);      	


  // Pair Info
  myChain->SetBranchAddress("pair_M",&pair_M);
  //
  myChain->SetBranchAddress("pair_Run2IdLevel",&pair_Run2IdLevel); 
  // 0 : notID | 1 : Veto | 2 : Loose | 3 : Medium | 4 : Tight
  //myChain->SetBranchAddress("pair_cuts",&pair_cuts); 
  // 0 : noCut | 1 : VBTF 95 | 2 : VBTF 80 | 3 : VBTF 60
  myChain->SetBranchAddress("pair_HLT_Ele27_cut",&pair_HLT_Ele27_cut); 
  myChain->SetBranchAddress("pair_fidu",&pair_fidu);
  //
  myChain->SetBranchAddress("pair_eta",&pair_eta);
  myChain->SetBranchAddress("pair_sclEta",&pair_sclEta);
  myChain->SetBranchAddress("pair_sclPhi",&pair_sclPhi);
  myChain->SetBranchAddress("pair_phi",&pair_phi);
  myChain->SetBranchAddress("pair_RCTeta",&pair_RCTeta);
  myChain->SetBranchAddress("pair_RCTphi",&pair_RCTphi);
  //
  myChain->SetBranchAddress("pair_charge",&pair_charge);
  myChain->SetBranchAddress("pair_pT",&pair_pT);
  myChain->SetBranchAddress("pair_eT",&pair_eT);
  myChain->SetBranchAddress("pair_sclEt",&pair_sclEt);
  myChain->SetBranchAddress("pair_E",&pair_E);

  myChain->SetBranchAddress("pair_TTetaVect",&pair_TTetaVect);
  myChain->SetBranchAddress("pair_TTphiVect",&pair_TTphiVect);
  myChain->SetBranchAddress("pair_TTetVect",&pair_TTetVect);  

  myChain->SetBranchAddress("pair_L1iso",&pair_L1iso);
  myChain->SetBranchAddress("pair_L1noniso",&pair_L1noniso);
  myChain->SetBranchAddress("pair_L1iso_M",&pair_L1iso_M);
  myChain->SetBranchAddress("pair_L1noniso_M",&pair_L1noniso_M);
  myChain->SetBranchAddress("pair_L1Stage2",&pair_L1Stage2);
  myChain->SetBranchAddress("pair_L1Stage2_isoflag",&pair_L1Stage2_isoflag);
  myChain->SetBranchAddress("pair_L1Stage2_emul",&pair_L1Stage2_emul);
  myChain->SetBranchAddress("pair_L1Stage2_emul_isoflag",&pair_L1Stage2_emul_isoflag);
  myChain->SetBranchAddress("pair_L1Stage2_emul_nTT",&pair_L1Stage2_emul_nTT);
  myChain->SetBranchAddress("pair_L1Stage1_emul",&pair_L1Stage1_emul);
  myChain->SetBranchAddress("pair_L1Stage1_emul_isoflag",&pair_L1Stage1_emul_isoflag);

  //
  myChain->SetBranchAddress("pair_RCTetVect",&pair_RCTetVect);
  myChain->SetBranchAddress("pair_RCTetaVect",&pair_RCTetaVect);
  myChain->SetBranchAddress("pair_RCTphiVect",&pair_RCTphiVect);
  myChain->SetBranchAddress("pair_L1isoVect",&pair_L1isoVect);
  myChain->SetBranchAddress("pair_L1nonisoVect",&pair_L1nonisoVect);
  myChain->SetBranchAddress("pair_L1isoVect_M",&pair_L1isoVect_M);
  myChain->SetBranchAddress("pair_L1nonisoVect_M",&pair_L1nonisoVect_M);
  myChain->SetBranchAddress("pair_L1Stage1_emul_isoVect",&pair_L1Stage1_emul_isoVect);
  myChain->SetBranchAddress("pair_L1Stage1_emul_nonisoVect",&pair_L1Stage1_emul_nonisoVect);

  if (debug) cout << "DEFINITION OF THE TAG AND THE BE" << endl;
  // DEFINITION OF THE TAG AND THE PROBE //
  int iTag, iProbe;
  vector<int> goodTags, goodProbes;
  int askedIdLevelIndex[2]; // 0:tag ; 1:probe //previously 'cutsDef'
  //int cutsDef[2]; // 0:tag ; 1:probe
  bool passingCuts[2][2]; // [iEle][iTPDefCuts] iTPDefCuts=0(1) <-> Tag(Probe)
  //TString cutsName[5] = {"noCuts","WP95","WP80","WP60","HLT_Ele27"};
  TString idLevels[5] = {"notID","Veto","Loose","Medium","Tight"}; //previously 'cutsName'
  double etDef[2] = {tagEt,probeEt};
  TString askedIdLevel[2] = {tagCuts,probeCuts}; //previously 'asked'
  for(int iEle=0 ; iEle<2 ; iEle++) {
    for(int iCut=0 ; iCut<5 ; iCut++) {
      if(askedIdLevel[iEle]==idLevels[iCut]) askedIdLevelIndex[iEle] = iCut ;
    }
  }

  /*TString asked[2] = {tagCuts,probeCuts};
  for(int iEle=0 ; iEle<2 ; iEle++) {
    for(int iCut=0 ; iCut<5 ; iCut++) {
      if(asked[iEle]==cutsName[iCut]) cutsDef[iEle] = iCut ;
    }
    }*/

  if (debug) cout << "USEFUL VARIABLES"<< endl;
  // USEFUL VARIABLES //
  int b_hlt, nCurrentRun,nCurrentEvt;
  nCurrentRun = nCurrentEvt = 0;
  TString filename = "";
  TString flag = "";
  vector<int> IdxPt;
  TString nameTagProbe[2] = {"TAG","PROBE"} ;

  MAPTT adcTT;
  MAPTT::iterator iterTT;
  pair<int,int> coords;
  vector<int> ietaTT, iphiTT;


  // TRIGGERING PROPERTIES //
  vector<int> trigHLT;
  vector<int> firedEG[nColl][2]; // [N/M/S2/S2E/S2EI/S1I/S1E/S1EI/S2I][tag/probe]
  int iTrig[nColl][2]; // [N/M/S2/S2E/S2EI/S1I/S1E/S1EI/S2I][tag/probe]
  int tagEG, tagIsoEG, tagIsoEGer;

  // COUNTERS //
  int nCases[3] = {0,0,0}; // all ; same ; different
  int nZeroStatus[3] = {0,0,0}; // ToBeZeroed ; Zeroed ; NotZeroed
  int nEB_NotToZero_Zeroed = 0;
  int nEE_NotToZero_Zeroed = 0;
  int n_fails_fidu = 0;

  if (debug) cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  int nEle[nECAL][nColl], nMatched[nECAL][nColl][nEG], nLost[nECAL][nColl][nEG], nAppeared[nECAL][nColl][nEG];
  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {
      nEle[iECAL][iColl] = 0;
      for(int iEG=0 ; iEG<nEG ; iEG++) {
	nMatched[iECAL][iColl][iEG] = nLost[iECAL][iColl][iEG] = nAppeared[iECAL][iColl][iEG] = 0;
      }
    }
  }

  int nProbesTrig15, nProbesNoTrig15, nCurrentRun_ineff, nCurrentLumi_ineff,nCurrentRun_elepairs, nCurrentLumi_elepairs;
  nProbesTrig15 = nProbesNoTrig15 = 0;
  nCurrentRun_ineff = nCurrentLumi_ineff = -1;
  nCurrentRun_elepairs = nCurrentLumi_elepairs = -1;

  // Z MASS PEAK PLOT //
  TH1F* h_Mee[3];
  TString nameZ[3] = {"h_EBEB","h_EEEE","h_EBEE"};
  for(int i=0;i<3;i++) {
    h_Mee[i] = new TH1F(nameZ[i],nameZ[i],80,50,130);
  }

  if (debug) cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;

  // ETA/PHI EFFICIENCY MAP
  TH2F* h_map_trig[nEG];
  TH2F* h_map_all;

  for(int iEG=0 ; iEG<nEG ; iEG++) {
    h_map_trig[iEG] = new TH2F("h_map_trig_EG"+EGval[iEG],"eta/phi map ele trig EG"+EGval[iEG],100,-2.5,2.5,100,-3.14,3.14);
  }
  h_map_all = new TH2F("h_map_all","eta/phi map all ele",100,-2.5,2.5,100,-3.14,3.14);
  if (debug) cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;

  // ETA/ET PHI/ET ETA/PHI plots
  TH2F* h_Map[3];
  TString name_Map[3] = {"h_EtaEt","h_PhiEt","h_EtaPhi"};
  h_Map[0] = new TH2F(name_Map[0],name_Map[0],50,-2.5,2.5,50,0,100);
  h_Map[1] = new TH2F(name_Map[1],name_Map[1],63,-3.14,3.14,50,0,100);
  h_Map[2] = new TH2F(name_Map[2],name_Map[2],50,-2.5,2.5,63,-3.14,3.14);

  TH1F* h_Map1D[3];
  TString name_Map1D[3] = {"h_Eta","h_Phi","h_Et"};
  h_Map1D[0] = new TH1F(name_Map1D[0],name_Map1D[0],50,-2.5,2.5); // eta
  h_Map1D[1] = new TH1F(name_Map1D[1],name_Map1D[1],63,-3.14,3.14); // phi
  h_Map1D[2] = new TH1F(name_Map1D[2],name_Map1D[2],50,0,100); // eT
 
  if (debug) cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  
  // JSON FILE READER //
  
  string jsonDir = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/" ;

  vector < string > jsonFile;
  
  // jsonFile.push_back(jsonDir + "Collisions11/7TeV/Reprocessing/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt");
  // jsonFile.push_back(jsonDir + "Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt");
  // jsonFile.push_back(jsonDir + "Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt");
  // jsonFile.push_back(jsonDir + "Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt");
  //jsonFile.push_back(jsonDir + "Collisions12/8TeV/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt");
  jsonFile.push_back(jsonDir + "Collisions16/13TeV/Cert_271036-282092_13TeV_PromptReco_Collisions16_JSON.txt");

  map<int, vector<pair<int, int> > > jsonMap[jsonFile.size()] ;
 
  // for(int i=0 ; i<jsonFile.size() ; i++)
  //   jsonMap[i] = readJSONFile(jsonFile[i]);

  if (debug) cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  





  // LOOP OVER PAIRS //
  
  int numEntries = myChain->GetEntries () ;
  int nProcess = numEntries;
  if(nEntries>=0 && nEntries<numEntries)
    nProcess = nEntries;
  cout << "will process " << nProcess << "/" << numEntries << " pairs" << endl;

  for (int iPair = 0 ; iPair < nProcess ; iPair++ ) { 
    
    if (iPair > 0 && iPair/1000. ==  int(iPair/1000.) ) cout << "processing " << iPair << "th entry" << endl;
    
    if(debug) cout << "is getting the entry #" << iPair << endl;
    myChain->GetEntry (iPair) ;
    if(debug) cout << "got it" << endl;

    // Get HLT information //
    int trigHLT_2011[10] = {trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8,
			    trig_isL1SingleEG12, trig_isL1SingleEG12_old, trig_is_HLT_Ele15, trig_is_HLT_Ele17,
			    trig_is_HLT_Ele8_calIso, trig_is_HLT_Ele8_tkIso} ;
    int trigHLT_2010[4] = {trig_isUnbiased, trig_isL1SingleEG2, trig_isL1SingleEG5, trig_isL1SingleEG8} ;

    trigHLT.clear();

  
    for(int iHLT=0 ; iHLT<nHLT ; iHLT++) {
      if(data2011) 
	trigHLT.push_back( trigHLT_2011[iHLT] ) ;
      else 
	trigHLT.push_back( trigHLT_2010[iHLT] ) ;
    }


    // LOG : FILE / RUN / EVENT / HLT / L1 //
    if( iPair==0 || filename != myChain->GetFile()->GetName() ) {
      filename = myChain->GetFile()->GetName() ;
      outlog << "<------ File : " << filename << " ------>" << endl;
    }

    if( iPair==0 || nRun!=nCurrentRun ) {
      nCurrentRun = nRun ;
      outlog << endl << "---- nRun=" << nRun << " ----" << endl ;
    }
    

    // HLT SELECTION //
    /*b_hlt=1;
    if(tag_hlt=="HLT_EG2")
      b_hlt = trig_isL1SingleEG2 ;
    if(tag_hlt=="HLT_EG5")
      b_hlt = trig_isL1SingleEG5 ;
    if(tag_hlt=="HLT_EG8")
      b_hlt = trig_isL1SingleEG8 ;
    if(tag_hlt=="HLT_unbias")
      b_hlt = trig_isUnbiased ;
    if(tag_hlt=="noHLT")
      b_hlt = 1;
    if(b_hlt != 1) {      
      continue;
      }*/

    if(debug) cout << "passed HLT selection" << endl;

    // MASS CUT //
    if( massCut1 <= massCut2 ) {
      if( pair_M < massCut1 || pair_M > massCut2 ) continue;
    }
    else {
      if( pair_M < massCut2 || pair_M > massCut1 ) continue;
    }
    

    if(debug) cout << "passed mass cut" << endl;

    // CHARGE OPPOSITION //
    if( checkCharge )
      if( pair_charge[0]*pair_charge[1] != -1 )
	{
	  continue;
	}
    
    if(debug) cout << "passed opposite-charge selection" << endl;
    
    // CHECK TRIGGERING PROPERTIES //
    for(int iEle=0 ; iEle<2 ; iEle++) {
      for(int iColl=0 ; iColl<nColl ; iColl++) {
	firedEG[iColl][iEle].clear();
	firedEG[iColl][iEle].resize(nEG,0);
      }
    }

  
    for(int iEle=0 ; iEle<2 ; iEle++) {
      for(int iR=0 ; iR < 10 ; iR++) {
	globalFireL1_Normal( pair_RCTetVect[iEle][iR] , pair_L1nonisoVect[iEle][iR],
			     pair_L1isoVect[iEle][iR], firedEG[0][iEle], menu);  //menu); 
	
	globalFireL1_Normal( pair_RCTetVect[iEle][iR] , pair_L1nonisoVect_M[iEle][iR],
			     pair_L1isoVect_M[iEle][iR], firedEG[1][iEle], menu);//menu); 

	globalFireL1_Normal( pair_RCTetVect[iEle][iR] , pair_L1Stage1_emul_nonisoVect[iEle][iR],
			     pair_L1Stage1_emul_isoVect[iEle][iR], firedEG[6][iEle], menu);//menu); 
	

	//if(debug) cout << "the pairRCTvalue for emul*****"<<endl<<pair_RCTetVect[iEle][iR] << " " << pair_L1nonisoVect_M[iEle][iR] << " " << pair_L1isoVect_M[iEle][iR] << endl;

	//T&P with isolated candidates, iso only for probe
	/*if(iEle==0)
	  globalFireL1_Normal( pair_RCTetVect[iEle][iR] , pair_L1nonisoVect[iEle][iR],
			       pair_L1isoVect[iEle][iR], firedEG[5][iEle], menu);
			       else*/
	  globalFireL1_Normal( pair_RCTetVect[iEle][iR] , 0,
			       pair_L1isoVect[iEle][iR], firedEG[5][iEle], menu);

	  /*if(iEle==0)
	  globalFireL1_Normal( pair_RCTetVect[iEle][iR] , pair_L1Stage1_emul_nonisoVect[iEle][iR],
			       pair_L1Stage1_emul_isoVect[iEle][iR], firedEG[7][iEle], menu);
			       else*/
	  globalFireL1_Normal( pair_RCTetVect[iEle][iR] , 0,
			       pair_L1Stage1_emul_isoVect[iEle][iR], firedEG[7][iEle], menu);


      }

      globalFireL1_Normal( 20 , pair_L1Stage2[iEle],
			 pair_L1Stage2[iEle], firedEG[2][iEle], menu); //No reason to use RCT energy for Stage 2, 20 dummy instead? No Iso, non iso for now

      globalFireL1_Normal( 20 , pair_L1Stage2_emul[iEle],
			   pair_L1Stage2_emul[iEle], firedEG[3][iEle], menu);



      //T&P with isolated candidates, iso only for probe
      /*if(iEle==0){
	globalFireL1_Normal( 20 , 0,
			     pair_L1Stage2_emul[iEle], firedEG[4][iEle], menu);
	globalFireL1_Normal( 20 , 0,
			     pair_L1Stage2[iEle], firedEG[8][iEle], menu);
      }
      else{*/
	globalFireL1_Normal( 20 , 0,
			     (pair_L1Stage2_emul_isoflag[iEle]==1)*pair_L1Stage2_emul[iEle], firedEG[4][iEle], menu);
	globalFireL1_Normal( 20 , 0,
			     (pair_L1Stage2_isoflag[iEle]==1)*pair_L1Stage2[iEle], firedEG[8][iEle], menu);
	//}

	firedEG[9][iEle][0] = (pair_L1Stage2[iEle]>=40) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=24) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=22 && fabs(pair_sclEta[iEle])<2.1) ;
	firedEG[9][iEle][1] = (pair_L1Stage2[iEle]>=34) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=26) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=24 && fabs(pair_sclEta[iEle])<2.1) ;
	firedEG[9][iEle][2] = (pair_L1Stage2[iEle]>=30) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=28) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=26 && fabs(pair_sclEta[iEle])<2.1) ;
	firedEG[9][iEle][3] = (pair_L1Stage2[iEle]>=32) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=25) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=25 && fabs(pair_sclEta[iEle])<2.1) ;
	firedEG[9][iEle][4] = (pair_L1Stage2[iEle]>=40);
	firedEG[9][iEle][5] = (pair_L1Stage2[iEle]>=24 && pair_L1Stage2_isoflag[iEle]==1);
	firedEG[9][iEle][6] = (pair_L1Stage2[iEle]>=32) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=30) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=28 && fabs(pair_sclEta[iEle])<2.1) ;
	firedEG[9][iEle][7] = (pair_L1Stage2[iEle]>=35) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=30) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=26 && fabs(pair_sclEta[iEle])<2.1) ;
	firedEG[9][iEle][8] = (pair_L1Stage2[iEle]>=34) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=32) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=30 && fabs(pair_sclEta[iEle])<2.1) ;
	firedEG[9][iEle][9] = (pair_L1Stage2[iEle]>=36) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=30) || (pair_L1Stage2_isoflag[iEle]==1 && pair_L1Stage2[iEle]>=28 && fabs(pair_sclEta[iEle])<2.1) ;

    }

    


    if(debug) {
      //for(int iColl=1 ; iColl<nColl ; iColl++) {
      int iColl=6;
	for(int iEle=0 ; iEle<2 ; iEle++) {
	  cout << "Ele" << iEle << " :: " ;
	  for(int i=0 ; i<(int)firedEG[iColl][iEle].size() ; i++) {
	    cout << "EG" << i << " : " << firedEG[iColl][iEle][i] << " | ";
	  }
	  cout << endl;
	}      
	//}
      }

    for(int iColl=0 ; iColl<nColl ; iColl++)
      for(int iEle=0 ; iEle<2 ; iEle++)
	iTrig[iColl][iEle] = 0;
    
    for(int iColl=0 ; iColl<nColl ; iColl++)    
      for(int iEle=0 ; iEle<2 ; iEle++)
	for(int iEG=0 ; iEG<nEG ; iEG++)
	  if(firedEG[iColl][iEle][iEG]==1) iTrig[iColl][iEle]= EGval[iEG].Atoi();//iEG;
	
    if(debug) cout << "TRIG emul fanbo------********* : ele0=" << iTrig[1][0] << " | ele1=" << iTrig[1][1] << endl;

    if(tag_hlt=="HLT_Ele30_Mass55")
      tagEG=40, tagIsoEG=24, tagIsoEGer=22; // EG40, IsoEG24, IsoEG22er
    else if(tag_hlt=="HLT_Ele23_Ele12")
      tagEG=30, tagIsoEG=28, tagIsoEGer=26; // EG40, IsoEG24, IsoEG22er
    

    // DECIDE WHO IS THE TAG/PROBE //
    iTag = iProbe = -1 ;
    goodTags.clear();
    goodProbes.clear();
    
    for(int i=0 ; i<2 ; i++)
      for(int j=0 ; j<2 ; j++)
	passingCuts[i][j] = false;

    for(int iEle=0 ; iEle<2 ; iEle++) {
      for(int iTP=0 ; iTP<2 ; iTP++) {
	if( pair_Run2IdLevel[iEle] >= askedIdLevelIndex[iTP] ) { 
	  if( pair_fidu[iEle]==1 ) {
	    passingCuts[iEle][iTP]=true; 
	    if(debug) cout<<"pass cuts"<<endl; 
	  }
	}
      }
    }

    /*for(int iEle=0 ; iEle<2 ; iEle++) {
      for(int iTP=0 ; iTP<2 ; iTP++) {
	if( cutsDef[iTP]==4 ) {
	  if( pair_HLT_Ele27_cut[iEle]==1 ) { 
	    if( pair_fidu[iEle]==1 ) {
	      passingCuts[iEle][iTP]=true; 
	      if(debug) cout<<"pass cuts" << endl; 
	    }
	  }
	}
	else {
	  if( pair_cuts[iEle] >= cutsDef[iTP] ) { 
	    if( pair_fidu[iEle]==1 ) {
	      passingCuts[iEle][iTP]=true; 
	      if(debug) cout<<"pass cuts"<<endl; 
	    }
	  }
	}
      }
      }*/



    for( int iEle=0 ; iEle<2 ; iEle++ ) {
      if ( passingCuts[iEle][0] ) {
	if(debug) cout << "ele pass cuts tag" << endl;
	if(pair_sclEt[iEle] >= etDef[0]) {
	  if(debug) cout << "et is ok : " << pair_sclEt[iEle] << endl; 
	  if (debug) cout << __LINE__ << " "<< iTrig[0][iEle] << " " << tagEG << endl; 
	  if(isMC && isStage2 && (iTrig[3][iEle] >= tagEG || iTrig[4][iEle] >= tagIsoEG || (iTrig[4][iEle] >= tagIsoEGer && abs(pair_sclEta[iEle])<2.1) ) )
	    //T&P with Stage-2 in MC
	    goodTags.push_back(iEle);
	  else if(isMC && !isStage2 && iTrig[6][iEle] >= tagEG)
	    //T&P with Stage-1 in MC
	    goodTags.push_back(iEle);
	  else if(!isMC && !isStage2 && iTrig[0][iEle] >= tagEG ) {
	    //This requirement has to be applied to the trigger used to take data
	    if(debug) cout << "triggering " << iTrig[0][iEle] << endl;
	    goodTags.push_back(iEle);
	  }
	  else if(!isMC && isStage2 && (iTrig[2][iEle] >= tagEG || iTrig[8][iEle] >= tagIsoEG || (iTrig[8][iEle] >= tagIsoEGer && abs(pair_sclEta[iEle])<2.1) ) ) {
	    //This requirement has to be applied to the trigger used to take data
	    if(debug) cout << "triggering " << iTrig[2][iEle] << endl;
	    goodTags.push_back(iEle);
	  }

	}
      }
      
      if ( passingCuts[iEle][1] ) {
	if(debug) cout << "ele pass cuts probe" << endl;
	if(pair_sclEt[iEle] >= etDef[1] ) {
	  if(debug) cout << "et is ok : " << pair_sclEt[iEle] << endl;
	  goodProbes.push_back(iEle);
	}
      }
    }
   
    if(debug) cout << "goodTags.size()=" << goodTags.size() << " | goodProbes.size()=" << goodProbes.size() << endl;

    if( goodTags.size()==0 || goodProbes.size()==0 ) continue; // no tag or no probe
 
    
    // Look every cases
    for(int iCaseT=0 ; iCaseT<(int)goodTags.size() ; iCaseT++) {
      
      iTag = goodTags[iCaseT];
      
	for(int iCaseP=0 ; iCaseP<(int)goodProbes.size() ; iCaseP++) {
	  iProbe = goodProbes[iCaseP];
	  if(iProbe==iTag) continue;
	  
	  // CUSTOM ELEPAIRS JSON
	  if( nRun!=nCurrentRun_elepairs ) {
	    if(nCurrentRun_elepairs!=-1) {
	      out_elepairs_json << "]," << endl;
	    }
	    nCurrentLumi_elepairs=-1;
	    out_elepairs_json << '"' << nRun << '"' << ": [" ;
	    nCurrentRun_elepairs = nRun;
	  }
	  if( nLumi!=nCurrentLumi_elepairs ) {
	    if( nCurrentLumi_elepairs!=-1 ) {
	      out_elepairs_json << "," ;
	    }
	    out_elepairs_json << "[" << nLumi << "," << nLumi << "]" ;
	    nCurrentLumi_elepairs=nLumi;
	  }
	  
	  
	  // FILL TRIGGER TREE FOR EFFICIENCY ///////////////////////
	  int idxECAL=-1;
	  if( fabs(pair_eta[iProbe])<1.479 ) idxECAL=0;
	  else idxECAL = 1;
	  
	  for(int iColl=0 ; iColl<nColl ; iColl++) {
	    sc_et[idxECAL][iColl] = pair_sclEt[iProbe];
	    ele_pt[idxECAL][iColl] = pair_pT[iProbe];	    
	    sc_eta[idxECAL][iColl] = pair_sclEta[iProbe];
	    sc_dr[idxECAL][iColl] = 0;
	    sc_L1_nTT[idxECAL][iColl] = pair_L1Stage2_emul_nTT[iProbe];
	    
	    PU_rhoCorr = _PU_rhoCorr;
	    PU_N = _PU_N;
	    n_vtx = vtx_N;
	    
	    if(iColl<2){
	      if(pair_L1iso[iProbe]>0 || pair_L1noniso[iProbe]>0)
		sc_dr[idxECAL][iColl] = 1;
	    }
	    else if(iColl==2){
	    if(pair_L1Stage2[iProbe]>0)
	      sc_dr[idxECAL][iColl] = 1;
	    }
	    else if(iColl==3){
	      if(pair_L1Stage2_emul[iProbe]>0)
		sc_dr[idxECAL][iColl] = 1;
	    }
	    else if(iColl==4){
	      if(pair_L1Stage2_emul[iProbe]>0)
		sc_dr[idxECAL][iColl] = 1;
	    }
	    else if(iColl==5){
	      if(pair_L1iso[iProbe]>0 || pair_L1noniso[iProbe]>0)
		sc_dr[idxECAL][iColl] = 1;
	    }
	    else if(iColl>5 && iColl<8){
	      if(pair_L1Stage1_emul[iProbe]>0)
		sc_dr[idxECAL][iColl] = 1;
	    }
	    else if(iColl==8){
	      if(pair_L1Stage2[iProbe]>0)
		sc_dr[idxECAL][iColl] = 1;
	    }
	    
	    for(int iEG=0 ; iEG<nEG ; iEG++)
	      l1_bin[idxECAL][iColl][iEG]  = firedEG[iColl][iProbe][iEG];
	    
	    
	  }
	  
	  for(int iColl=0 ; iColl<nColl ; iColl++)
	    treenew[idxECAL][iColl]->Fill();
	  ////////////////////////////////////////////////////////////
	  //EBEE tree
	  for(int iColl=0 ; iColl<nColl ; iColl++) {
	    sc_et[2][iColl] = pair_sclEt[iProbe];
	    sc_eta[2][iColl] = pair_sclEta[iProbe];
	    ele_pt[2][iColl] = pair_pT[iProbe];
	    sc_dr[2][iColl] = 0;
	    sc_L1_nTT[2][iColl] = pair_L1Stage2_emul_nTT[iProbe];
	    
	    PU_rhoCorr = _PU_rhoCorr;
	    PU_N = _PU_N;
	    n_vtx = vtx_N;
	    
	    if(iColl<2){
	      if(pair_L1iso[iProbe]>0 || pair_L1noniso[iProbe]>0)
		sc_dr[2][iColl] = 1;
	    }
	    else if(iColl==2){
	    if(pair_L1Stage2[iProbe]>0)
	      sc_dr[2][iColl] = 1;
	    }
	    else if(iColl==3){
	      if(pair_L1Stage2_emul[iProbe]>0)
		sc_dr[2][iColl] = 1;
	    }
	    else if(iColl==4){
	      if(pair_L1Stage2_emul[iProbe]>0)
		sc_dr[2][iColl] = 1;
	    }
	    else if(iColl==5){
	      if(pair_L1iso[iProbe]>0 || pair_L1noniso[iProbe]>0)
		sc_dr[2][iColl] = 1;
	    }
	    else if(iColl>5 && iColl<8){
	      if(pair_L1Stage1_emul[iProbe]>0)
		sc_dr[2][iColl] = 1;
	    }
	    else if(iColl==8){
	      if(pair_L1Stage2[iProbe]>0)
		sc_dr[2][iColl] = 1;
	    }
	    
	    for(int iEG=0 ; iEG<nEG ; iEG++)
	      l1_bin[2][iColl][iEG]  = firedEG[iColl][iProbe][iEG];
	    
	    
	  }
	  
	  for(int iColl=0 ; iColl<nColl ; iColl++)
	    treenew[2][iColl]->Fill();



	  ////////////////////////////////////////////////////////////

	  
	  // FILL EFFICIENCY MAP TOOLS
	  h_map_all->Fill( pair_eta[iProbe] , pair_phi[iProbe] );
	  for(int iEG=0 ; iEG<nEG ; iEG++) {
	    if( firedEG[0][iProbe][iEG]==1 ) {
	      h_map_trig[iEG]->Fill( pair_eta[iProbe] , pair_phi[iProbe] );
	    }
	  }
	  
	  // Z MASS PEAK PLOT // EBEB EEEE EBEE
	  if( abs(pair_eta[iProbe])<1.479 && abs(pair_eta[iTag])<1.479 ) h_Mee[0]->Fill(pair_M);
	  if( abs(pair_eta[iProbe])>1.479 && abs(pair_eta[iTag])>1.479 ) h_Mee[1]->Fill(pair_M);
	  else h_Mee[2]->Fill(pair_M);
	  
	  // OUTPUT JSON //
	  if( nRun!=nCurrentRun_ineff ) {
	    if(nCurrentRun_ineff!=-1) {
	      out_custom_json << "]," << endl;
	    }
	    nCurrentLumi_ineff=-1;
	    out_custom_json << '"' << nRun << '"' << ": [" ;
	    nCurrentRun_ineff = nRun;
	  }
	  if( nLumi!=nCurrentLumi_ineff ) {
	    if( nCurrentLumi_ineff!=-1 ) {
	      out_custom_json << "," ;
	    }
	    out_custom_json << "[" << nLumi << "," << nLumi << "]" ;
	    nCurrentLumi_ineff=nLumi;
	  }
	  
	  
	  // LOG PROBES RUN:EVENT //
	  if( pair_sclEt[iProbe]>10. && firedEG[0][iProbe][20]==iEff ) {
	    out_probes_trigEG15 << "'" << nRun << ":" << nEvent << "'" ;
	    //out_probes_trigEG15 << "'" << nRun ;
	    if(iPair < nProcess-1 ) out_probes_trigEG15 << "," ;
	    nProbesTrig15++ ;
	    
	    
	    
	    out_probes_notrigEG15 << "'" << nRun << ":" << nEvent << "'" ;
	    out_notrigEG15_RLE << "'" << nRun << ":" << nLumi << ":" << nEvent << "'" ;
	    out_notrigEG15_R  << nRun ;
	    if(iPair < nProcess-1 ) {
	      out_probes_notrigEG15 << "," ;
	      out_notrigEG15_RLE << "," ;
	      out_notrigEG15_R << "," ;
	    }
	    nProbesNoTrig15++ ;
	  }

	  // COUNTING //
	  flag = "" ;
	  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
	    for(int iColl=0 ; iColl<nColl ; iColl++) {
	      nEle[iECAL][iColl]++ ;
	      for(int iEG=0 ; iEG<nEG ; iEG++) {
		if(firedEG[0][iProbe][iEG]==1)  nMatched[iECAL][0][iEG]++ ;
		
	      }
	    }
	  }
	  
	  // MAP THE TP //
	  adcTT.clear();
	  for(int t=0 ; t<trig_tower_N ; t++) {
	    coords = make_pair( trig_tower_ieta[t] , trig_tower_iphi[t] );
	    adcTT[coords].first = t ;
	  }
	  for(int t=0 ; t<trig_tower_N_emul ; t++) {
	    coords = make_pair( trig_tower_ieta_emul[t] , trig_tower_iphi_emul[t] );
	    iterTT = adcTT.find( coords );
	    if( iterTT != adcTT.end() )
	      adcTT[coords].second = t;
	    else {
	      outlog << "mapping problem" << endl;
	      adcTT[coords].first = -777;
	      adcTT[coords].second = t;
	    }
	  }
	  

	  if( (pair_sclEt[iProbe]>9.&&pair_sclEt[iProbe]<18.) && firedEG[0][iProbe][20]==iEff && firedEG[1][iProbe][20]==iEff_M&&fabs(pair_sclEta[iProbe]) > 1.566&&fabs(pair_sclEta[iProbe]) <2.5) {
	    
	    
            outlog<<"******------"<<endl;
            outlog<<"The Pair Number   "<<iPair<<endl;
            if(firedEG[0][iProbe][20]==1)
	      //std::cout<< "fired EG20 trigger"<<std::endl;
	      outlog<<"TAG   SCL(ET,ETA,PHI)    "<<"L1ISO(ET,ETA,PHi)     "<<"L1NONISO(ET,ETA,PHi)    "<<"ELECAN(ET,ETA,PHi)"<<endl
		    <<"("<<pair_sclEt[iTag]<<", "<<pair_sclEta[iTag]<<","<<pair_sclPhi[iTag]<<")"<<"      "
		    <<"("<<pair_L1iso[iTag]<<", "<<pair_RCTeta[iTag]<<", "<<pair_RCTphi[iTag]<<")"<<"      " 
		    <<"("<<pair_L1noniso[iTag]<<", "<<pair_RCTeta[iTag]<<", "<<pair_RCTphi[iTag]<<")"<<"      " 
		    <<"("<<pair_eT[iTag]<<", "<<", "<<pair_eta[iTag]<<", "<<pair_phi[iTag]<<")"<<"  "
		    <<"---------"<<endl;
            outlog<<"Prob   SCL(ET,ETA,PHI)    "<<"L1ISO(ET,ETA,PHi)     "<<"L1NONISO(ET,ETA,PHi)    "<<"ELECAN(ET,ETA,PHi)"<<endl
                  <<"("<<pair_sclEt[iProbe]<<", "<<pair_sclEta[iProbe]<<","<<pair_sclPhi[iProbe]<<")"<<"      "
                  <<"("<<pair_L1iso[iProbe]<<", "<<pair_RCTeta[iProbe]<<", "<<pair_RCTphi[iProbe]<<")"<<"      " 
                  <<"("<<pair_L1noniso[iProbe]<<", "<<pair_RCTeta[iProbe]<<", "<<pair_RCTphi[iProbe]<<")"<<"      " 
                  <<"("<<pair_eT[iProbe]<<", "<<", "<<pair_eta[iProbe]<<", "<<pair_phi[iProbe]<<")"<<endl;
	    
	    outlog << "Pair M=" << pair_M << endl;
	    
	    
	    
	    int IdxTagProbe[2] = {iTag , iProbe} ;
	    
	    // FILLING THE MAPS
	    h_Map[0]->Fill( pair_eta[IdxTagProbe[1]] , pair_eT[IdxTagProbe[1]] ); // eta/et
	    h_Map[1]->Fill( pair_phi[IdxTagProbe[1]] , pair_eT[IdxTagProbe[1]] ); // phi/et
	    h_Map[2]->Fill( pair_eta[IdxTagProbe[1]] , pair_phi[IdxTagProbe[1]] ); // eta/phi

	    h_Map1D[0]->Fill( pair_eta[IdxTagProbe[1]] ); // eta
	    h_Map1D[1]->Fill( pair_phi[IdxTagProbe[1]] ); // phi
	    h_Map1D[2]->Fill( pair_eT[IdxTagProbe[1]] ); // et
	    
	    
	  } //endif probe inefficient
	}
	// loop over goodProbes
    } // loop over goodTags
    
  } // loop over pairs
    


  // OUTPUT FILE //
  stringstream treeout; 
  treeout.str(""); 
  if(debug) cout << "test" << endl;
  cout<<"fileName="<<fileName<<endl;
  cout << "tree position " << fileName.find("tree"); 

  treeout <<dirOut << "/effi_TagProbe_" <<  fileName.substr(fileName.find("tree")); 
  if(debug) cout << "test" << endl;


  cout << "output file: " << treeout.str() << endl;
  TFile *outfile = new TFile(treeout.str().c_str(),"RECREATE");
  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {
      treenew[iECAL][iColl]->Write();
    }
  }
  outfile->Close();
  stringstream copyoutfile;
  copyoutfile.str("");
  int lenght = fileName.find("makePairs")-fileName.find("/store");
  copyoutfile << "cmsStage -f "<< treeout.str().c_str() <<  " "<< fileName.substr(fileName.find("/store"), lenght) << treeout.str().c_str();
  cout << copyoutfile.str()<< endl;
  system(copyoutfile.str().c_str());
  cout << "...filled !" << endl;

  // PLOTTING Z MASS

  stringstream outfilestr; 


  // PLOTTING MAPS

  cout << "Plotting maps" << endl;
  TCanvas* c_Maps[3];
  TString name_c_Maps[3] = {"_EtaEt","_PhiEt","_EtaPhi"};
  gStyle->SetPalette(1);
  for(int iC=0 ; iC<3 ; iC++) {
    outfilestr.str("");
    outfilestr << dirOut << "/maps"<< name_c_Maps[iC]<< ".gif";
    c_Maps[iC] = new TCanvas("c"+name_c_Maps[iC],"c"+name_c_Maps[iC],0,0,800,600);
    h_Map[iC]->Draw("colz");
    c_Maps[iC]->Print(outfilestr.str().c_str(),"gif");
  }

  TCanvas* c_Maps1D[3];
  TString name_c_Maps1D[3] = {"_Eta","_Phi","_eT"};
  gStyle->SetPalette(1);
  for(int iC=0 ; iC<3 ; iC++) {
    outfilestr.str("");
    outfilestr << dirOut << "/maps"<< name_c_Maps1D[iC]<< ".gif";
    c_Maps1D[iC] = new TCanvas("c"+name_c_Maps1D[iC],"c"+name_c_Maps1D[iC],0,0,800,600);
    h_Map1D[iC]->Draw("");
    c_Maps1D[iC]->Print(outfilestr.str().c_str(),"gif");
  }
  cout << "plotted" << endl;


  // Efficiency eta/phi map
  double n_ele_all[10000];
  double n_ele_trig[10000][nEG];
  double bin_eta[10000], bin_phi[10000], effi_val[10000];

  for(int iBinEta=0 ; iBinEta<100 ; iBinEta++) {
    for(int iBinPhi=0 ; iBinPhi<100 ; iBinPhi++) {

      bin_eta[ iBinEta*100 + iBinPhi ] = -2.5  + iBinEta*(5./100.);
      bin_phi[ iBinEta*100 + iBinPhi ] = -3.14 + iBinPhi*(6.28/100.);

      n_ele_all[ iBinEta*100 + iBinPhi ] = (double)( h_map_all->GetBinContent(iBinEta+1,iBinPhi+1) );

      for(int iEG=0 ; iEG<nEG ; iEG++) {
	n_ele_trig[ iBinEta*100 + iBinPhi ][iEG] = (double)( h_map_trig[iEG]->GetBinContent(iBinEta+1,iBinPhi+1) );

	if( n_ele_all[ iBinEta*100 + iBinPhi ] != 0 ) {
	  effi_val[ iBinEta*100 + iBinPhi ] = n_ele_trig[ iBinEta*100 + iBinPhi ][iEG] / n_ele_all[ iBinEta*100 + iBinPhi ];
	}
	else {
	  effi_val[ iBinEta*100 + iBinPhi ] = -0.1;
	}

      }

    }
  }


  TGraph2D* g_effi_map[nEG];
  TCanvas*  c_effi_map[nEG];

  //for(int iEG=0 ; iEG<nEG ; iEG++) {
  int jEG=15;
  g_effi_map[jEG] = new TGraph2D("g_effi_map_EG"+EGval[jEG],"Efficiency Eta/Phi Map : EG"+EGval[jEG],10000,bin_eta,bin_phi,effi_val);
  c_effi_map[jEG] = new TCanvas("c_effi_map_EG"+EGval[jEG],"Efficiency Eta/Phi Map : EG"+EGval[jEG],0,0,800,600);
  gStyle->SetPalette(1);
  g_effi_map[jEG]->Draw("colz");
  outfilestr.str("");
  outfilestr << dirOut << "/g_effi_map_EG"<<EGval[jEG]<<".gif"; 
  c_effi_map[jEG]->Print(outfilestr.str().c_str());
  outfilestr.str("");
  outfilestr << dirOut << "/g_effi_map_EG"<<EGval[jEG]<<".C"; 
  c_effi_map[jEG]->Print(outfilestr.str().c_str());
    //}

  // OUTPUT LOGS //
  if(debug) cout << "Filling output logs..." << endl;
  outfilestr.str("");
  outfilestr << dirOut << "/counters_selectPairs.csv";
  ofstream outCounters(outfilestr.str().c_str(), ios::out);

  for(int iECAL=0 ; iECAL<nECAL ; iECAL++ ) {
    for(int iColl=0 ; iColl<nColl ; iColl++) {
      outCounters << nEle[iECAL][iColl] << " | " ;
      for(int iEG=0 ; iEG<nEG ; iEG++) {
	outCounters << nMatched[iECAL][iColl][iEG] << "," ;
      }
    }

    outCounters << endl;
  }    

  outCounters << endl << "nProbesTrig15=" << nProbesTrig15 
	      << " | nProbesNoTrig15=" << nProbesNoTrig15 << endl;

  // MASS PLOTS output
  stringstream myoutputfile;
  myoutputfile.str("");
  myoutputfile <<dirOut<<  "/mass_Z_plot_"<< fileName.substr(fileName.find("tree"));
  TFile *outplot = new TFile(myoutputfile.str().c_str(),"RECREATE");
  for(int i=0 ; i<3 ; i++)
    h_Mee[i]->Write();
  outplot->Close();
  

  out_custom_json << ']' << endl << '}' ;
  out_elepairs_json << ']' << endl << '}' ;

  outlog << endl << endl;
  for(int i=0 ; i<3 ; i++)
    outlog << "-----" ;
  outlog << " RESULTS " ;
  outlog << endl
	 << "nCases = " << nCases[0] 
	 << endl;
  for(int i=0 ; i<9 ; i++)
    outlog << "-----" ;

  
  for (int i=0; i<3 ; i++){
    h_Mee[i]->~TH1F();
  }
  for(int iEG=0 ; iEG<nEG ; iEG++) {
    h_map_trig[iEG] -> ~TH2F(); 
  }
  h_map_all->~TH2F();

  h_Map[0]->~TH2F();
  h_Map[1]->~TH2F();
  h_Map[2]->~TH2F();

  h_Map1D[0] ->~TH1F();
  h_Map1D[1] ->~TH1F();
  h_Map1D[2] ->~TH1F();
    

  

}
