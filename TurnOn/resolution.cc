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
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include "baseFuncNad.h"
#include <sys/types.h>
#include <sys/stat.h>
//#include <boost/filesystem.hpp>
#include <iostream>
//
struct stat info;

//void resolution(string nameChain, string file,int iFile,string dirOut, TString dirIn);
void resolution(string nameChain, string file,int iFile,TString RunPhase);
string myjobid ;
using namespace std;

int main()
{


   resolution("produceNtuple/eIDSimpleTree","/afs/cern.ch/work/f/fmeng/FANBOWORKINGAREA/CMSSW_7_4_12/src/EGamma/ECGelec/test/ratioScan/tree_testFastSim_TPG_checknow_2015Cno_833.root",1,"2015C");




}

//void resolution(string nameChain, string file,int iFile,string dirOut, TString dirIn)
void resolution(string nameChain, string file,int iFile,TString RunPhase)
{
//   std::size_t pos = file.find("tree");
//   std::string str2 = file.substr(pos);
   TChain * myChain = new TChain(nameChain.c_str());
   myChain->Add(file.c_str());
   int numEntries = myChain->GetEntries () ;
   int ele_N;
   TClonesArray * electrons = new TClonesArray ("TLorentzVector");
   int ele_severityLevelSeed[10];
   TLorentzVector* cand[2];
   //TLorentzVector* cand_sure;
   bool cut_HLT_Ele27[2];
   double ele_sclEta[10], ele_sclEt[10],ele_sclPhi[10];
   double ele_ecalRecHitSumEt_dr03[10], ele_tkSumPt_dr03[10];
   double ele_hcalDepth1TowerSumEt_dr03[10], ele_hcalDepth2TowerSumEt_dr03[10];
   int ele_expected_inner_hits[10];
   double ele_deltaphiin[10], ele_deltaetain[10];
   double ele_he[10], ele_sigmaietaieta[10];
   double ele_conv_dist[10], ele_conv_dcot[10];
   double ele_fbrem[10];
   int ele_isConversion[10];
   int cutEle[2],fidu[2];
   int goodelectron[10];
   double ele_RCTetaVect[10][10], ele_RCTphiVect[10][10];
   double ele_RCTetVect[10][10];
   TH1F *resolu=new TH1F("resolution","barrel resolution ",20,-1,1);
   TLorentzVector total;


   myChain->SetBranchAddress("electrons",&electrons);
   myChain->SetBranchAddress("ele_N",   &ele_N);
   myChain->SetBranchAddress("ele_severityLevelSeed", &ele_severityLevelSeed);
   myChain->SetBranchAddress("ele_sclEta",&ele_sclEta);
   myChain->SetBranchAddress("ele_sclPhi",&ele_sclPhi);
   myChain->SetBranchAddress("ele_tkSumPt_dr03",&ele_tkSumPt_dr03);
   myChain->SetBranchAddress("ele_ecalRecHitSumEt_dr03", &ele_ecalRecHitSumEt_dr03);
   myChain->SetBranchAddress("ele_hcalDepth1TowerSumEt_dr03", &ele_hcalDepth1TowerSumEt_dr03);
   myChain->SetBranchAddress("ele_hcalDepth2TowerSumEt_dr03", &ele_hcalDepth2TowerSumEt_dr03);
   myChain->SetBranchAddress("ele_expected_inner_hits",&ele_expected_inner_hits);
   myChain->SetBranchAddress("ele_deltaphiin",&ele_deltaphiin);
   myChain->SetBranchAddress("ele_deltaetain",&ele_deltaetain);
   myChain->SetBranchAddress("ele_he",&ele_he); 
   myChain->SetBranchAddress("ele_sigmaietaieta",&ele_sigmaietaieta);
   myChain->SetBranchAddress("ele_conv_dist",&ele_conv_dist);
   myChain->SetBranchAddress("ele_conv_dcot",&ele_conv_dcot);
   myChain->SetBranchAddress("ele_fbrem",&ele_fbrem);
   myChain->SetBranchAddress("ele_isConversion",&ele_isConversion);
   myChain->SetBranchAddress("ele_RCTetaVect", &ele_RCTetaVect);
   myChain->SetBranchAddress("ele_RCTphiVect", &ele_RCTphiVect);
   myChain->SetBranchAddress("ele_RCTetVect", &ele_RCTetVect);




   for (int iEvent = 0 ; iEvent <numEntries; iEvent++ )
   {
      for (int number=0;number<10;number++)
           goodelectron[number]=0;
      myChain->GetEntry (iEvent);
//      cout << "<-- ele_N=" << ele_N << endl<< "--- electrons.size=" << electrons->GetSize() << endl;
      for( int iEle1=0 ; iEle1<ele_N ; iEle1++ )
      {
         if( ele_severityLevelSeed[iEle1] >= 3 ) continue;

         cand[0] = (TLorentzVector*) (electrons->At (iEle1)) ;
         cut_HLT_Ele27[0] = VBTFcuts( "HLT_Ele27", RunPhase,
                                      cand[0]->Pt(), cand[0]->Et(), ele_sclEta[iEle1], cand[0]->Eta(), ele_tkSumPt_dr03[iEle1], ele_ecalRecHitSumEt_dr03[iEle1],                                    ele_hcalDepth1TowerSumEt_dr03[iEle1], ele_hcalDepth2TowerSumEt_dr03[iEle1], ele_expected_inner_hits[iEle1],
                                      ele_deltaphiin[iEle1], ele_deltaetain[iEle1], ele_he[iEle1], ele_sigmaietaieta[iEle1],
                                      ele_conv_dist[iEle1], ele_conv_dcot[iEle1], ele_fbrem[iEle1], ele_isConversion[iEle1] ) ;
      //  cout<<"cut HLT_Ele27 value="<<cut_HLT_Ele27<<endl;
      // check if ele is a good tag candidate : pass VBTF 95 and has pT>5 GeV
         cutEle[0]= 0;
         cutEle[0] = whichCuts( RunPhase, cand[0]->Pt(), cand[0]->Et(), ele_sclEta[iEle1], cand[0]->Eta(), ele_tkSumPt_dr03[iEle1], ele_ecalRecHitSumEt_dr03[iEle1],
                               ele_hcalDepth1TowerSumEt_dr03[iEle1], ele_hcalDepth2TowerSumEt_dr03[iEle1], ele_expected_inner_hits[iEle1],
                               ele_deltaphiin[iEle1], ele_deltaetain[iEle1], ele_he[iEle1], ele_sigmaietaieta[iEle1],
                               ele_conv_dist[iEle1], ele_conv_dcot[iEle1], ele_fbrem[iEle1], ele_isConversion[iEle1] ) ;
        
         fidu[0]= 0;
         if ( fabs(ele_sclEta[iEle1]) < 2.5 && ( fabs(ele_sclEta[iEle1]) > 1.566 || fabs(ele_sclEta[iEle1])<1.4442 ) )
         fidu[0]= 1 ;
         if( cutEle>0 && cand[0]->Et()>=5. ) 
         {
         //   if(debug) cout << "--- ele #" << iEle1 << " is a good tag candidate" << endl;


            for( int iEle2=0 ; iEle2<ele_N ; iEle2++ ) {
             //  if(debug) cout << "----- looks Ele #" << iEle2 << endl;

               if(iEle2<=iEle1 ) continue;
               cand[1] = (TLorentzVector*) (electrons->At (iEle2)) ;
               if( ele_severityLevelSeed[iEle2] >= 3 ) continue;
               cut_HLT_Ele27[1] = VBTFcuts( "HLT_Ele27", RunPhase,
                                            cand[1]->Pt(), cand[1]->Et(), ele_sclEta[iEle1], cand[1]->Eta(), ele_tkSumPt_dr03[iEle1], ele_ecalRecHitSumEt_dr03[iEle1],
                                            ele_hcalDepth1TowerSumEt_dr03[iEle1], ele_hcalDepth2TowerSumEt_dr03[iEle1], ele_expected_inner_hits[iEle1],
                                            ele_deltaphiin[iEle1], ele_deltaetain[iEle1], ele_he[iEle1], ele_sigmaietaieta[iEle1],
                                            ele_conv_dist[iEle1], ele_conv_dcot[iEle1], ele_fbrem[iEle1], ele_isConversion[iEle1] ) ;

               cutEle[1] = whichCuts( RunPhase, cand[1]->Pt(), cand[0]->Et(), ele_sclEta[iEle2], cand[1]->Eta(), ele_tkSumPt_dr03[iEle2], ele_ecalRecHitSumEt_dr03[iEle2],
                                      ele_hcalDepth1TowerSumEt_dr03[iEle2], ele_hcalDepth2TowerSumEt_dr03[iEle2], ele_expected_inner_hits[iEle2],
                                      ele_deltaphiin[iEle2], ele_deltaetain[iEle2], ele_he[iEle2], ele_sigmaietaieta[iEle2],
                                      ele_conv_dist[iEle2], ele_conv_dcot[iEle2], ele_fbrem[iEle2], ele_isConversion[iEle2] ) ;
               fidu[1] = 0;
               if ( fabs(ele_sclEta[iEle2]) < 2.5 && ( fabs(ele_sclEta[iEle2]) > 1.566 || fabs(ele_sclEta[iEle2])<1.4442 ) )
                  fidu[1] = 1 ;

               if( cutEle[1]>0 && iEle2<=iEle1 ) continue;

             //  if(debug) cout << "---> OK to form a pre-selected pair <--" << endl;
               total = (*cand[0]) + (*cand[1]) ;
               if( total.M() < 30. ) continue;
               goodelectron[iEle1]=1;        
               goodelectron[iEle2]=1;        

            }
    
         }  

      } 
    // electron and L1 EG eta matching
       for(int number1=0;number1<10;number1++)
       {
         if(goodelectron[number1]>0)
         { 
      //      cand_sure=(TLorentzVector*) (electrons->At (number1)) ; 
//      myChain->SetBranchAddress("ele_RCTetaVect", &ele_RCTetaVect);
//   myChain->SetBranchAddress("ele_RCTphiVect", &ele_RCTphiVect);
            for(int icc1=0;icc1<10;icc1++)
            { 
               // if((cand_sure->Eta()==ele_RCTetaVect[number1][icc])&&(cand_sure->Phi()==ele_RCTphiVect[number1][icc])&&(ele_RCTphiVect[number1][icc]>-800))           
              //  cout<<"there is something*******"<<endl;
              //  the supercluse match the L1 in which the highest TT is within the supercluster
               for(int icc2=0;icc2<10;icc2++){
                  if(ele_RCTetaVect[number1][icc2]==ele_sclEta[icc1]&&ele_RCTphiVect[number1][icc2]==ele_sclPhi[icc1])   
                       cout<<"there is something*******"<<endl;
                  }

            }
         }


       } 
     cout<<"event number"<<iEvent<<endl; 
   }

}
