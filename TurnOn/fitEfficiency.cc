//g++  -o fitEfficiency FuncCB.cc fitEfficiency.cc  `root-config --cflags --libs` -L  $ROOFITSYS/lib -lRooFit -lRooFitCore -I$ROOFITSYS/include
//g++ -o  fitEfficiency fitEfficiency.cc  `root-config --cflags --libs` -L $ROOFITSYS/lib -lRooFit -lRooFitCore -I$ROOFITSYS/include
/////////////////////////////////////////////////////
// FITTING WITH A ROOFIT-USER DEFINED CRYSTAL BALL //
/////////////////////////////////////////////////////

#ifndef DEF_FASTEFF
#define DEF_FASTEFF

// General C++
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <glob.h>
#include <cstdlib>
#include <stdio.h> 
#include <stdlib.h>
#include <vector>
#include <glob.h>

// RooFit headers
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooEfficiency.h"
#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooHist.h"
#include "RooWorkspace.h"

// Root headers
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
#include "TTree.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"
#include "RooVoigtian.h"


// Personal headers
#include "FuncCB.h"
//#include "RooCruijff.h"
#include "tdrstyle.h"

#define DEBUG 0

using namespace RooFit ;

TString fileIn;

using namespace std;
vector<string> globVector(const string& pattern){
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> files;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    //cout << string(glob_result.gl_pathv[i]) << endl;
    files.push_back(string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return files;
}

bool dirExists(const char *path)
{
    struct stat info;

    if(stat( path, &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
      return true;
    else
        return false;
}


void  fastEfficiencySilvia(unsigned int iEG1, unsigned int iEG2, int iECAL1, int iColl1, int iECAL2, int iColl2, 
			     TString fileIn0, TString fileIn1, TString label0, TString label1,
			     string dirIn, int run_number, int lumi, int nCPU, int color1, int style1, int color2, int style2,
			   TString probe, TString tag, bool logx);

void  fastEfficiencySilvia_3plots(unsigned int iEG1, unsigned int iEG2, unsigned int iEG3, int iECAL1, int iColl1, int iECAL2, int iColl2, int iECAL3, int iColl3,
				  TString fileIn0, TString fileIn1, TString fileIn2, TString label0, TString label1, TString label2,
			     string dirIn, int run_number, int lumi, int nCPU, int color1, int style1, int color2, int style2, int color3, int style3,
				  TString probe, TString tag, bool logx);

void  fastEfficiencySilvia_4plots(unsigned int iEG1, unsigned int iEG2, unsigned int iEG3, unsigned int iEG4, int iECAL1, int iColl1, int iECAL2, int iColl2, int iECAL3, int iColl3, int iECAL4, int iColl4,
				  TString fileIn0, TString fileIn1, TString fileIn2, TString fileIn3, TString label0, TString label1, TString label2, TString label3,
			     string dirIn, int run_number, int lumi, int nCPU, int color1, int style1, int color2, int style2, int color3, int style3, int color4, int style4,
				  TString probe, TString tag, bool logx);


void  fitManyFiles(unsigned int iEG, int iECAL1, int iColl1, int iECAL2, int iColl2, 
			     TString fileIn0, TString label0, TString label1,
			     string dirIn, TString lumi, int nCPU, int color1, int style1, int color2, int style2,
			     TString probe, TString tag);


void loadPresentationStyle();
string myjobid; 
stringstream mydir ; 

//-------
int main(int argc, char**argv){

  if (argc<2) {
    cout << "argument needed"<< endl;
    return 1; 
  }
  myjobid = argv[1];
  mydir.str("");
  mydir << myjobid <<"/"<< "selectPairsDir"; 
  if (DEBUG) cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  TApplication app("App",&argc, argv);

  fastEfficiencySilvia_4plots(30,30,30,30,1,2,1,3,1,4,1,8, "effi_TagProbe_tree.root","effi_TagProbe_tree.root","effi_TagProbe_tree.root","effi_TagProbe_tree.root",
  			   "SingleEG30 prod", "SingleEG30 emul", "SingleIsoEG30 emul","SingleIsoEG30 prod",mydir.str(), 0,585, 4, 
  			      kBlack, kFullCircle, kRed, kOpenSquare, kBlue, kOpenTriangleUp, kGreen, kFullTriangleDown, "Loose", "Medium",false);
  app.Run(); 
  
  return 0;
}

void loadPresentationStyle(){
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(0.01);
  gStyle->SetLabelOffset(0.005, "XYZ");
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetHistLineWidth(2);
  setTDRStyle();
}

void fastEfficiencySilvia_4plots(unsigned int iEG1, unsigned int iEG2, unsigned int iEG3, unsigned int iEG4, int iECAL1, int iColl1, int iECAL2, int iColl2, int iECAL3, int iColl3, int iECAL4, int iColl4,
			  TString fileIn0="effi_TagProbe_tree_changed.root",TString fileIn1="effi_TagProbe_tree_changed.root", TString fileIn2="effi_TagProbe_tree_changed.root", TString fileIn3="effi_TagProbe_tree_changed.root", 
			  TString label0 = "FGVB ratio theshold 0.90", TString label1 = "FGVB ratio threshold 0.95", TString label2 = "FGVB ratio threshold 0.95", TString label3 = "FGVB ratio threshold 0.95",
			  string dirIn="result/"//=mydir.str()
				 , int run_number=0, int lumi=0, int nCPU=4, 
			 int color1=kBlack, int style1=kFullCircle, int color2=kRed, int style2=kOpenSquare, int color3=kBlue, int style3=kOpenTriangleUp, int color4=kGreen, int style4=kOpenTriangleDown,
				 TString probe="WP80", TString tag="WP80",
				 bool logx=true){
  // STYLE //
  if (DEBUG) cout << __LINE__ << endl;
  gROOT->Reset();
  loadPresentationStyle();  
  gROOT->ForceStyle();
  
  iEG1= iEG1/5;
  iEG2= iEG2/5;
  iEG3= iEG3/5;
  iEG4= iEG4/5;

  const int nEG = 71;
  double thres[nEG];
  for(int i=0 ; i<nEG ; i++) thres[i]=5*i;
  if (DEBUG) cout << __LINE__ << endl;

  TString names[nEG];
  ostringstream ossi;
  for(int i=0;i<(int)nEG;i++) {
    ossi.str("");
    ossi << thres[i] ;
    names[i] = ossi.str();
  }
  if (DEBUG) cout << __LINE__ << endl;

  // NAMES //
  const int nECAL=3;
  const int nColl=10;

  TString name_leg_ecal[nECAL] = {"Barrel","Endcaps","Inclusive"};
  TString name_leg_coll[nColl] = {"Online","Emulation"};  
  if (DEBUG) cout << __LINE__ << endl;

  TString name_ecal[nECAL] = {"_EB","_EE","_EBEE"};
  TString name_coll[nColl] = {"_N","_M", "_S2", "_S2E", "_S2EI", "_S1I", "_S1E","_S1EI", "_S2I", "_S2test"};
  if (DEBUG) cout << __LINE__ << endl;

  stringstream dirResults;
  dirResults.str("");
  dirResults <<  dirIn <<  "/turnons/EG"<<names[iEG1]<<"/" ;
  stringstream createDir; 
  createDir.str(""); 
  if(DEBUG) createDir<< dirIn << "/turnons"; 
  if (dirExists(createDir.str().c_str())==false) {
    createDir.str("");
    createDir << "mkdir "<<  dirIn << "/turnons";
    system(createDir.str().c_str());
  }
  if(DEBUG) cout << __LINE__ << " " << dirExists(dirResults.str().c_str()) << endl;
  createDir.str(""); 
  createDir<< "mkdir " << dirResults.str().c_str()<< ""; 
  if (dirExists(dirResults.str().c_str())==false) system(createDir.str().c_str());
  cout << __LINE__ << " " << dirExists(dirResults.str().c_str()) << endl;
  cout << "DIR CREATED" << endl; 
  stringstream name_image ;
  name_image.str("");
  name_image << dirResults.str() << "eff_EG"<<names[iEG1]<<"_EG"<<names[iEG2]<<"_EG"<<names[iEG3]<<"_tag"<<tag<<"_probe"<<probe<<name_ecal[iECAL1]<<name_coll[iColl1]<<"_vs"<<name_ecal[iECAL2]<<name_coll[iColl2] << "_vs"<<name_ecal[iECAL3]<<name_coll[iColl3]<< "_vs"<<name_ecal[iECAL4]<<name_coll[iColl4];
  if (DEBUG) cout << __LINE__ << endl;

  // Output log //
  if(logx)
    name_image << "_logx.txt"; 
  else
    name_image << ".txt"; 
  ofstream fichier(name_image.str().c_str(), ios::out);


  // BINNING //
  const int nbins[nEG] = {20,20,20,19,20,20,20,20,21,21,
			  21,22,29,29};

  if (DEBUG) cout << __LINE__ << endl;
  Double_t bins_20[20] = {10, 14,  18,  20,  22, 26, 30, 35,  40, 45, 50, 60, 70, 100, 150, 200, 400, 600, 800, 1000}; // EG20
  Double_t bins_21[20] = {9, 11, 13, 15, 17, 19, 20, 21, 22, 23, 25, 27, 30, 35, 40, 45, 50, 60, 70, 150}; // EG21
  Double_t bins_22[20] = {10, 12, 14, 16, 18, 20, 21, 22, 23, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 150}; // EG22

  Double_t bins_23[19] = {11, 13, 15, 17, 19, 21, 22, 23, 24, 25, 27, 30, 35, 40, 45, 50, 60, 70, 150}; // EG23

  Double_t bins_24[20] = {10, 12, 14, 16, 18, 20, 22, 23, 24, 25, 26, 28, 30, 35, 40, 45, 50, 60, 70, 100}; // EG24
  Double_t bins_25[20] = {11, 13, 15, 17, 19, 21, 23, 24, 25, 26, 27, 29, 30, 35, 40, 45, 50, 60, 70, 100}; // EG25
  Double_t bins_26[20] = {10, 12, 14, 16, 18, 20, 22, 24, 25, 26, 27, 28, 30, 35, 40, 45, 50, 60, 70, 150}; // EG26
  Double_t bins_27[20] = {11, 13, 15, 17, 19, 21, 23, 25, 26, 27, 28, 29, 33, 35, 40, 45, 50, 60, 70, 150}; // EG27
  if (DEBUG) cout << __LINE__ << endl;

  Double_t bins_28[21] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 27, 28, 29, 30, 32, 35, 40, 45, 50, 60, 70, 150}; // EG28
  Double_t bins_29[21] = {11, 13, 15, 17, 19, 21, 23, 25, 27, 28, 29, 30, 31, 33, 35, 40, 45, 50, 60, 70, 150}; // EG29
  Double_t bins_30[21] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 29, 30, 31, 32, 35, 40, 45, 50, 60, 70, 100}; // EG30

  Double_t bins_40[22] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 38, 39, 40, 42, 45, 50, 60, 70, 150}; // EG40
  Double_t bins_50[29] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 48, 50, 55, 60, 70, 90, 110, 130, 150, 170, 190}; // EG50
  Double_t bins_60[27] = {10,  15,  20,  25,  30,  34,  36,  38,  40, 42,
			  44,  48,  52,  56,  60, 64, 68, 90, 110, 130, 
			  150, 170, 190, 210, 230, 250, 270 }; // EG60

  if (DEBUG) cout << __LINE__ << endl;

  vector< Double_t* > bins;
  bins.push_back( bins_20 ); bins.push_back( bins_21 ); bins.push_back( bins_22 ); bins.push_back( bins_23 ); bins.push_back( bins_24 ); 
  bins.push_back( bins_25 ); bins.push_back( bins_26 ); bins.push_back( bins_27 ); bins.push_back( bins_28 ); bins.push_back( bins_29 ); 
  if (DEBUG) cout << __LINE__ << endl;

  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_30 );
  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_40 );
  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_50 );
  for(int iV=0 ; iV<10 ; iV++) bins.push_back( bins_60 );
   bins.push_back( bins_30 );
   bins.push_back( bins_40 );
   bins.push_back( bins_50 );
   bins.push_back( bins_60 );
   // cout << __LINE__  << " " << nbins[13] << " " ; 
   // for (int ibin =0; ibin < nbins[13] ; ibin++){
   //   cout << bins[13][ibin] << " " ;
   // }
   cout << endl;
   // RooBinning binning = RooBinning(26,  bins_60, "binning");
   RooBinning binning = RooBinning(21,  bins_40, "binning");
  //RooBinning binning2 = RooBinning(nbins[iEG2]-1, bins[iEG2], "binning");
  //RooBinning binning3 = RooBinning(nbins[iEG2]-1, bins[iEG2], "binning");

  if (DEBUG) cout << __LINE__ << endl;


  //INPUT DATA //
  stringstream filename ;
  filename.str(""); 
  filename << dirIn << "/" <<fileIn0; 
  TFile* f0 = TFile::Open(filename.str().c_str());
  filename.str(""); 
  filename << dirIn << "/" <<fileIn1;   
  TFile* f1 = TFile::Open(filename.str().c_str());
  filename.str(""); 
  filename << dirIn << "/" <<fileIn2;   
  TFile* f2 = TFile::Open(filename.str().c_str());
  filename.str(""); 
  filename << dirIn << "/" <<fileIn3;   
  TFile* f3 = TFile::Open(filename.str().c_str());

  TTree* treenew;
  TTree* treenew_2;
  TTree* treenew_3;
  TTree* treenew_4;

  treenew = (TTree*) f0->Get( "treenew"+name_ecal[iECAL1]+name_coll[iColl1] ) ;
  treenew_2 = (TTree*) f1->Get( "treenew"+name_ecal[iECAL2]+name_coll[iColl2] ) ;
  treenew_3 = (TTree*) f2->Get( "treenew"+name_ecal[iECAL3]+name_coll[iColl3] ) ;
  treenew_4 = (TTree*) f3->Get( "treenew"+name_ecal[iECAL4]+name_coll[iColl4] ) ;

  if (DEBUG) cout << __LINE__ << endl;

  TString name_scet[4], name_scdr[4], name_l1bin[4];
  name_scet[0] = "sc_et"+name_ecal[iECAL1]+name_coll[iColl1];
  name_scet[1] = "sc_et"+name_ecal[iECAL2]+name_coll[iColl2];
  name_scet[2] = "sc_et"+name_ecal[iECAL3]+name_coll[iColl3];
  name_scet[3] = "sc_et"+name_ecal[iECAL4]+name_coll[iColl4];

  name_scdr[0] = "sc_dr"+name_ecal[iECAL1]+name_coll[iColl1];
  name_scdr[1] = "sc_dr"+name_ecal[iECAL2]+name_coll[iColl2];
  name_scdr[2] = "sc_dr"+name_ecal[iECAL3]+name_coll[iColl3];
  name_scdr[3] = "sc_dr"+name_ecal[iECAL4]+name_coll[iColl4];
 
  name_l1bin[0] = "l1_"+names[iEG1]+name_ecal[iECAL1]+name_coll[iColl1];
  name_l1bin[1] = "l1_"+names[iEG2]+name_ecal[iECAL2]+name_coll[iColl2];
  name_l1bin[2] = "l1_"+names[iEG3]+name_ecal[iECAL3]+name_coll[iColl3];
  name_l1bin[3] = "l1_"+names[iEG4]+name_ecal[iECAL4]+name_coll[iColl4];
  cout << __LINE__ << " " << name_l1bin[0] << " " << name_l1bin[1] <<  " " << name_l1bin[2] <<  " " << name_l1bin[3] << endl;
  
  RooRealVar et_plot(name_scet[0],name_scet[0],0,1000) ;
  RooRealVar dr(name_scdr[0],name_scdr[0],0.5,1.5) ; 
  RooRealVar et_plot2(name_scet[1],name_scet[1],0,1000) ;
  RooRealVar dr2(name_scdr[1],name_scdr[1],0.5,1.5) ;
  RooRealVar et_plot3(name_scet[2],name_scet[2],0,1000) ;
  RooRealVar dr3(name_scdr[2],name_scdr[2],0.5,1.5) ;
  RooRealVar et_plot4(name_scet[3],name_scet[3],0,1000) ;
  RooRealVar dr4(name_scdr[3],name_scdr[3],0.5,1.5) ;
  if (DEBUG) cout << __LINE__ << endl;

  // Acceptance state cut (1 or 0)
  RooCategory cut(name_l1bin[0],name_l1bin[0]) ;
  cut.defineType("accept",1) ;
  cut.defineType("reject",0) ;
  RooCategory cut2(name_l1bin[1],name_l1bin[1]) ;
  cut2.defineType("accept",1) ;
  cut2.defineType("reject",0) ;
  RooCategory cut3(name_l1bin[2],name_l1bin[2]) ;
  cut3.defineType("accept",1) ;
  cut3.defineType("reject",0) ;
  RooCategory cut4(name_l1bin[3],name_l1bin[3]) ;
  cut4.defineType("accept",1) ;
  cut4.defineType("reject",0) ;
  
  // PARAMETRES ROOFIT CRYSTAL BALL
  // PARAMETRES ROOFIT CRYSTAL BALL
  RooRealVar norm("norm","N",1,0.6,1.5);
  RooRealVar alpha("alpha","#alpha",7.7,5.5,8.8);
  RooRealVar n("n","n",1.5,1.01,2.1);
  RooRealVar mean("mean","mean",5*iEG1,20*iEG1/10.,5*iEG1+20*iEG1/10.);
  RooRealVar sigma("sigma","#sigma",4.0,3.01,5);
  RooRealVar sigmaR("sigmaR","#sigma",4.0,0.01,5);
  RooRealVar alphaR("alphaR","#alpha",5.5,3.5,10.);

  RooRealVar norm2("norm2","N",1,0.6,1.5);
  RooRealVar alpha2("alpha2","#alpha",6.,5.01,10.);
  RooRealVar n2("n2","n",1.5,.51,36);
  RooRealVar mean2("mean2","mean",5*iEG2,20*iEG2/10.,5*iEG2+20*iEG2/10.);
  //mean.setVal(thres[iEG]);
  RooRealVar sigma2("sigma2","#sigma",4.0,3.01,5);
  RooRealVar sigmaR2("sigmaR2","#sigma",4.0,0.01,5);
  RooRealVar alphaR2("alphaR2","#alpha",5.5,3.5,10.);

  RooRealVar norm3("norm3","N",1.,0.62,1.6);
  RooRealVar alpha3("alpha3","#alpha",6.,0.05,10.);
  RooRealVar n3("n3","n",1.5,.5,35);
  RooRealVar mean3("mean3","mean",5*iEG3,20*iEG3/10.,10*iEG3);
  //mean.setVal(thres[iEG]);
  RooRealVar sigma3("sigma3","#sigma",4.,0.01,5);
  RooRealVar sigmaR3("sigmaR3","#sigma",4.0,0.01,5);
  RooRealVar alphaR3("alphaR3","#alpha",5.5,3.5,10.);

  RooRealVar norm4("norm4","N",1.,0.62,1.5);
  RooRealVar alpha4("alpha4","#alpha",6.,0.05,10.);
  RooRealVar n4("n4","n",1.5,.5,35.);
  RooRealVar mean4("mean4","mean",5*iEG4,20*iEG4/10.,10*iEG4);
  //mean.setVal(thres[iEG]);
  RooRealVar sigma4("sigma4","#sigma",4.,0.01,5);
  RooRealVar sigmaR4("sigmaR4","#sigma",4.0,0.01,5);
  RooRealVar alphaR4("alphaR4","#alpha",4.,1.,10.);

  FuncCB cb("cb","Crystal Ball Integree",et_plot,mean,sigma,alpha,n,norm) ;
  FuncCB cb2("cb2","Crystal Ball Integree",et_plot2,mean2,sigma2,alpha2,n2,norm2) ;
  FuncCB cb3("cb3","Crystal Ball Integree",et_plot3,mean3,sigma3,alpha3,n3,norm3) ;
  FuncCB cb4("cb4","Crystal Ball Integree",et_plot4,mean4,sigma4,alpha4,n4,norm4) ;

  RooRealVar offset("offset","offset",.5,-1.01,1.01);
  RooRealVar slope("slope","slope",.5,0.05,10.);
  RooRealVar scale("scale", "scale", 1., 0.2, 2.); 

  
  RooFormulaVar effFormula("effFormula", "0.5*(TMath::Erf((et_plot-1)/0.5)+1)",et_plot) ;


 if (DEBUG) cout << __LINE__ << endl;
  
  // EFFICIENCY //
 RooEfficiency eff("eff","efficiency",cb,cut,"accept");
 //RooEfficiency eff("eff","efficiency",effFormula, cut,"accept" ) ;
  RooEfficiency eff2("eff2","efficiency",cb2,cut2,"accept");
  RooEfficiency eff3("eff3","efficiency",cb3,cut3,"accept");
  RooEfficiency eff4("eff4","efficiency",cb4,cut4,"accept");

  RooDataSet dataSet("data","data", RooArgSet(et_plot, cut), Import(*treenew));//,Import(*treenew)); 
  if (DEBUG) cout << __LINE__ << endl;
  RooDataSet dataSet2("data2","data2", RooArgSet(et_plot2, cut2), Import(*treenew_2));//,Import(*treenew_2));
  if (DEBUG) cout << __LINE__ << endl;
  RooDataSet dataSet3("data3","data3", RooArgSet(et_plot3, cut3), Import(*treenew_3));//,Import(*treenew_2));
  if (DEBUG) cout << __LINE__ << endl;
  RooDataSet dataSet4("data4","data4", RooArgSet(et_plot4, cut4), Import(*treenew_4));
  cout << __LINE__ << endl;
  dataSet.Print();
  dataSet2.Print();
  dataSet3.Print();
  dataSet4.Print();
  if (DEBUG) cout << __PRETTY_FUNCTION__ << " " << __LINE__ << endl;
  
  // PLOT //
  RooPlot* frame  = et_plot.frame(Bins(180000),Title("Fitted efficiency")) ;
  RooPlot* frame2 = et_plot2.frame(Bins(180000),Title("Fitted efficiency")) ;
  RooPlot* frame3 = et_plot3.frame(Bins(180000),Title("Fitted efficiency")) ;
  RooPlot* frame4 = et_plot4.frame(Bins(180000),Title("Fitted efficiency")) ;

  dataSet.plotOn(frame, Binning(binning), Efficiency(cut), DrawOption("LP"), MarkerColor(color1), LineColor(color1), MarkerStyle(style1), LineWidth(2) );
  dataSet2.plotOn(frame2, Binning(binning), Efficiency(cut2), DrawOption("LP"),MarkerColor(color2), LineColor(color2), MarkerStyle(style2), LineWidth(2) );
  dataSet3.plotOn(frame3, Binning(binning), Efficiency(cut3), DrawOption("LP"),MarkerColor(color3), LineColor(color3), MarkerStyle(style3), LineWidth(2) );
  dataSet4.plotOn(frame4, Binning(binning), Efficiency(cut4), DrawOption("LP"),MarkerColor(color4), LineColor(color4), MarkerStyle(style4), LineWidth(2) );
  if (DEBUG) cout << __LINE__ << endl;


  /////////////////////// FITTING /////////////////////////////

  double fit_cuts_min = 25;
  double fit_cuts_max = 1000;
  if (DEBUG) cout << __LINE__ << endl;

  et_plot.setRange("interesting",fit_cuts_min,fit_cuts_max);
  et_plot2.setRange("interesting",fit_cuts_min,fit_cuts_max);
  et_plot3.setRange("interesting",fit_cuts_min,fit_cuts_max);
  et_plot4.setRange("interesting",fit_cuts_min,fit_cuts_max);

  RooFitResult* roofitres1 = new RooFitResult("roofitres1","roofitres1");
  RooFitResult* roofitres2 = new RooFitResult("roofitres2","roofitres2");
  RooFitResult* roofitres3 = new RooFitResult("roofitres3","roofitres3");
  RooFitResult* roofitres4 = new RooFitResult("roofitres4","roofitres4");
  if (DEBUG) cout << __LINE__ << endl;

  fichier << "Fit characteristics :"   << endl ;
  fichier << "EG "     << names[iEG1] << endl ;
  fichier << "Fit Range , EB Coll : [" << fit_cuts_min << "," << fit_cuts_max << "]" << endl ;
  fichier << "Fit Range , EE Coll : [" << fit_cuts_min << "," << fit_cuts_max << "]" << endl ;
  fichier << "----------------------"  << endl ;
  if (DEBUG) cout << __LINE__ << endl;

  // // Fit #1 //
  // roofitres1 = eff.fitTo(dataSet,ConditionalObservables(et_plot),Range("interesting"), Minos(kFALSE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE));
  // if (DEBUG) cout << __LINE__ << endl;

  // cb.plotOn(frame,LineColor(color1),LineWidth(2));

  // double res_norm1  = norm.getVal();
  // double err_norm1  = norm.getErrorLo();
  // double res_mean1  = mean.getVal();
  // double err_mean1  = mean.getError();
  // double res_sigma1 = sigma.getVal();
  // double err_sigma1 = sigma.getError();
  // double res_n1     = n.getVal();
  // double err_n1     = n.getError();
  // double res_alpha1 = alpha.getVal();
  // double err_alpha1 = alpha.getError();

  // fichier << "<----------------- EB ----------------->" << endl
  // 	  << "double res_mean="  << res_mean1  << "; "
  // 	  << "double res_sigma=" << res_sigma1 << "; "
  // 	  << "double res_alpha=" << res_alpha1 << "; "
  // 	  << "double res_n="     << res_n1     << "; "
  // 	  << "double res_norm="  << res_norm1  << "; "
  // 	  << endl
  // 	  << "double err_mean="  << err_mean1  << "; "
  // 	  << "double err_sigma=" << err_sigma1 << "; "
  // 	  << "double err_alpha=" << err_alpha1 << "; "
  // 	  << "double err_n="     << err_n1     << "; "
  // 	  << "double err_norm="  << err_norm1  << "; "
  // 	  << endl;

  // // Fit #2 //
  // roofitres2 = eff2.fitTo(dataSet2,ConditionalObservables(et_plot2),Range("interesting"),Minos(kFALSE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE));
 
  // cb2.plotOn(frame2,LineColor(color2),LineWidth(2));

  // double res_norm2  = norm2.getVal();
  // double err_norm2  = norm2.getErrorLo();
  // double res_mean2  = mean2.getVal();
  // double err_mean2  = mean2.getError();
  // double res_sigma2 = sigma2.getVal();
  // double err_sigma2 = sigma2.getError();
  // double res_n2     = n2.getVal();
  // double err_n2     = n2.getError();
  // double res_alpha2 = alpha2.getVal();
  // double err_alpha2 = alpha2.getError();

  // fichier << "<----------------- EE ----------------->" << endl
  // 	  << "double res_mean="  << res_mean2  << "; "
  // 	  << "double res_sigma=" << res_sigma2 << "; "
  // 	  << "double res_alpha=" << res_alpha2 << "; "
  // 	  << "double res_n="     << res_n2     << "; "
  // 	  << "double res_norm="  << res_norm2  << "; "
  // 	  << endl
  // 	  << "double err_mean="  << err_mean2  << "; "
  // 	  << "double err_sigma=" << err_sigma2 << "; "
  // 	  << "double err_alpha=" << err_alpha2 << "; "
  // 	  << "double err_n="     << err_n2     << "; "
  // 	  << "double err_norm="  << err_norm2  << "; "
  // 	  << endl;

  // // Fit #3 //
  // roofitres3 = eff3.fitTo(dataSet3,ConditionalObservables(et_plot3),Range("interesting"),Minos(kFALSE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE));
 
  // cb3.plotOn(frame3,LineColor(color3),LineWidth(2));

  // double res_norm3  = norm3.getVal();
  // double err_norm3  = norm3.getErrorLo();
  // double res_mean3  = mean3.getVal();
  // double err_mean3  = mean3.getError();
  // double res_sigma3 = sigma3.getVal();
  // double err_sigma3 = sigma3.getError();
  // double res_n3     = n3.getVal();
  // double err_n3     = n3.getError();
  // double res_alpha3 = alpha3.getVal();
  // double err_alpha3 = alpha3.getError();

  // fichier << "<----------------- EE ----------------->" << endl
  // 	  << "double res_mean="  << res_mean3  << "; "
  // 	  << "double res_sigma=" << res_sigma3 << "; "
  // 	  << "double res_alpha=" << res_alpha3 << "; "
  // 	  << "double res_n="     << res_n3     << "; "
  // 	  << "double res_norm="  << res_norm3  << "; "
  // 	  << endl
  // 	  << "double err_mean="  << err_mean3  << "; "
  // 	  << "double err_sigma=" << err_sigma3 << "; "
  // 	  << "double err_alpha=" << err_alpha3 << "; "
  // 	  << "double err_n="     << err_n3     << "; "
  // 	  << "double err_norm="  << err_norm3  << "; "
  // 	  << endl;

  // // Fit #4 //
  // roofitres4 = eff4.fitTo(dataSet4,ConditionalObservables(et_plot4),Range("interesting"),Minos(kFALSE),Warnings(kFALSE),NumCPU(nCPU),Save(kTRUE));
 
  // cb4.plotOn(frame4,LineColor(color4),LineWidth(2));

  // double res_norm4  = norm4.getVal();
  // double err_norm4  = norm4.getErrorLo();
  // double res_mean4  = mean4.getVal();
  // double err_mean4  = mean4.getError();
  // double res_sigma4 = sigma4.getVal();
  // double err_sigma4 = sigma4.getError();
  // double res_n4     = n4.getVal();
  // double err_n4     = n4.getError();
  // double res_alpha4 = alpha4.getVal();
  // double err_alpha4 = alpha4.getError();

  // fichier << "<----------------- EE ----------------->" << endl
  // 	  << "double res_mean="  << res_mean4  << "; "
  // 	  << "double res_sigma=" << res_sigma4 << "; "
  // 	  << "double res_alpha=" << res_alpha4 << "; "
  // 	  << "double res_n="     << res_n4     << "; "
  // 	  << "double res_norm="  << res_norm4  << "; "
  // 	  << endl
  // 	  << "double err_mean="  << err_mean4  << "; "
  // 	  << "double err_sigma=" << err_sigma4 << "; "
  // 	  << "double err_alpha=" << err_alpha4 << "; "
  // 	  << "double err_n="     << err_n4     << "; "
  // 	  << "double err_norm="  << err_norm4  << "; "
  // 	  << endl;
    

  ////////////////////////////  DRAWING PLOTS AND LEGENDS /////////////////////////////////
  TCanvas* ca = new TCanvas("ca","Trigger Efficiency") ;

  ca->SetGridx();
  ca->SetGridy();
  ca->cd();


  if(logx)
    gPad->SetLogx();
  gPad->SetObjectStat(1);

  frame->GetYaxis()->SetRangeUser(0,1.05);
  frame->GetXaxis()->SetRangeUser(5,150.);// set to 1000 to see the plateau
  frame->GetYaxis()->SetTitle("Efficiency");
  frame->GetXaxis()->SetTitle("E_{T} [GeV]");
  frame->Draw() ;

  frame2->GetYaxis()->SetRangeUser(0,1.05);
  frame2->GetXaxis()->SetRangeUser(1,150.);// set to 1000 to see the plateau
  frame2->GetYaxis()->SetTitle("Efficiency");
  frame2->GetXaxis()->SetTitle("E_{T} [GeV]");
  frame2->Draw("same") ;

  frame3->GetYaxis()->SetRangeUser(0,1.05);
  frame3->GetXaxis()->SetRangeUser(1,150.);// set to 1000 to see the plateau
  frame3->GetYaxis()->SetTitle("Efficiency");
  frame3->GetXaxis()->SetTitle("E_{T} [GeV]");
  frame3->Draw("same") ;

  frame4->GetYaxis()->SetRangeUser(0,1.05);
  frame4->GetXaxis()->SetRangeUser(1,150.);// set to 1000 to see the plateau
  frame4->GetYaxis()->SetTitle("Efficiency");
  frame4->GetXaxis()->SetTitle("E_{T} [GeV]");
  frame4->Draw("same") ;

  TH1F *SCeta1 = new TH1F("SCeta1","SCeta1",50,-2.5,2.5);
  TH1F *SCeta2 = new TH1F("SCeta2","SCeta2",50,-2.5,2.5);
  TH1F *SCeta3 = new TH1F("SCeta3","SCeta3",50,-2.5,2.5);
  TH1F *SCeta4 = new TH1F("SCeta4","SCeta4",50,-2.5,2.5);

  SCeta1->SetLineColor(color1) ;
  SCeta1->SetMarkerColor(color1);
  SCeta1->SetMarkerStyle(style1);

  SCeta2->SetLineColor(color2) ;
  SCeta2->SetMarkerColor(color2);
  SCeta2->SetMarkerStyle(style2);

  SCeta3->SetLineColor(color3) ;
  SCeta3->SetMarkerColor(color3);
  SCeta3->SetMarkerStyle(style3);
  
  SCeta4->SetLineColor(color4) ;
  SCeta4->SetMarkerColor(color4);
  SCeta4->SetMarkerStyle(style4);

  //TLegend *leg = new TLegend(0.2,0.435,0.6,0.560,NULL,"brNDC"); // mid : x=353.5
  TLegend * leg = new TLegend(0.5,0.25,0.85,0.55,NULL,"brNDC");
  if(logx && iEG1>15)
    leg = new TLegend(0.19,0.435,0.44,0.560,NULL,"brNDC");

  leg->SetLineColor(1);
  leg->SetTextColor(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.025);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetMargin(0.2);
  leg->SetHeader("#bf{Endcaps}");
  leg->AddEntry(SCeta1,label0,"p");
  leg->AddEntry(SCeta2,label1,"p");
  leg->AddEntry(SCeta3,label2,"p");
  leg->AddEntry(SCeta4,label3,"p");
  leg->Draw();

  TLegend * leg2 = new TLegend(0.16,0.725,0.58,0.905,NULL,"brNDC");
  leg2->SetBorderSize(0);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.03);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  
  leg2->AddEntry("NULL","L1 Trigger EG"+names[iEG1],"h");

  TLatex *texl = new TLatex(10,1.06,"CMS Preliminary, #sqrt{s}=13 TeV");
  // if(run_number>0 && lumi>0)
  //   texl = new TLatex(5,1.06,Form("CMS Preliminary, #sqrt{s}=13 TeV, Run %i %i pb^{-1}",run_number,lumi));
  // else if(lumi>0)
  //   texl = new TLatex(5,1.06,Form("CMS Preliminary, #sqrt{s}=13 TeV, %i pb^{-1}",lumi));
  texl->SetTextSize(0.03);
  //texl->Draw("same");


//  TPaveText *pt2 = new TPaveText(0.220,0.605,0.487,0.685,"brNDC"); // mid : x=353.5                                         
  TPaveText *pt2 = new TPaveText(0.170,0.605,0.437,0.685,"brNDC"); 
  pt2->SetLineColor(1);
  pt2->SetTextColor(1);
  pt2->SetTextFont(42);
  pt2->SetTextSize(0.03);
  pt2->SetFillColor(kWhite);
  pt2->SetShadowColor(kWhite);
  pt2->AddText("L1 E/Gamma Trigger");
  pt2->AddText("Electrons from Z");
  //pt2->Draw();
  
  //TString name_image="eff_EG20_2012_12fb";

  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, ".cxx").c_str());
  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, ".png").c_str());//+".png","png");
  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, ".gif").c_str());//+".gif","gif");
  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, ".pdf").c_str());//+".pdf","pdf");
  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, ".ps", 0, 3).c_str()); //,+".ps","ps");
  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, ".eps").c_str());//,+".eps","eps");

  frame->GetXaxis()->SetRangeUser(30,100.);
  //texl->DrawLatex(30,1.06, "CMS Preliminary, #sqrt{s}=13 TeV");
  ca->Update();
  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, "_zoomed.pdf").c_str());
  ca->Print(name_image.str().replace(name_image.str().find(".txt"), 4, "_zoomed.png").c_str());
  

  /////////////////////////////
  // SAVE THE ROO FIT RESULT //
  /////////////////////////////
  TFile * lastoutputFile = new TFile(name_image.str().replace(name_image.str().find(".txt"), 4, "_fitres.root", 0,12).c_str(), "RECREATE"); //name_image+"_fitres.root", "RECREATE"); 
  lastoutputFile->cd();
  // cout << "old configuration, max eff  " << cb.value(100.) << endl;
  // cout << "new configuration, max eff  " << cb2.value(100.) << endl;
  // cout << "old configuration, 0.95 eff " << 0.95*(cb.value(100.)) << " at " << cb.find095Et(cb.value(100.)) << endl;
  // cout << "new configuration, 0.95 eff " << 0.95*(cb2.value(100.))<< " at " << cb2.find095Et(cb2.value(100.)) << endl;

  // double oldMaxEff, old95Eff, old95Et, newMaxEff, new95Eff, new95Et;
  // oldMaxEff= cb.value(100.);
  // newMaxEff= cb2.value(100.);
  // old95Et  = cb.find095Et(cb.value(100.));
  // new95Et  = cb2.find095Et(cb2.value(100.));
  // old95Eff = cb.value(old95Et);
  // new95Eff = cb2.value(new95Et);
  
  TTree * lastTree = new TTree("Efficiency","Efficiency");
  // lastTree->Branch("oldMaxEff",&oldMaxEff,"oldMaxEff/D");
  // lastTree->Branch("old95Eff" ,&old95Eff ,"old95Eff/D" );
  // lastTree->Branch("old95Et"  ,&old95Et  ,"old95Et/D"  );
  // lastTree->Branch("newMaxEff",&newMaxEff,"newMaxEff/D");
  // lastTree->Branch("new95Eff" ,&new95Eff ,"new95Eff/D" );
  // lastTree->Branch("new95Et"  ,&new95Et  ,"new95Et/D"  );

  lastTree->Fill();
  lastTree->Write();
  
  // RooWorkspace *w = new RooWorkspace("workspace","workspace") ;

  // w->import(dataSet);
  // w->import(dataSet2);
  // w->import(dataSet3);
  // w->import(dataSet4);  

  // w->import(*roofitres1,"roofitres1");
  // w->import(*roofitres2,"roofitres2");
  // w->import(*roofitres3,"roofitres3");
  // w->import(*roofitres4,"roofitres4");

  // cout << "CREATES WORKSPACE : " << endl;
  // w->Print();
  // w->Write(); 

  ca->Write(); 

  lastoutputFile->Close();
  f1->Close();

  gSystem->Exit(0);
}




#endif
