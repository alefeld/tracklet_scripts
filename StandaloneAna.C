#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <TCanvas.h>
#include <TProfile.h>

#include "FPGAEvent.h"

TTree *fpgaTree;
FPGAEvent *fpgaEvent; 
int numEvts;


bool checkType(int type, bool filterOnGood) {

  bool pass=false;
  
  if(filterOnGood) {

    // filter by the ones that are good
    pass=false;
    if(fabs(type)==11) return true; //electron
    if(fabs(type)==13) return true; //muon
    if(fabs(type)==211) return true; //charged pions
    if(fabs(type)==2212) return true; //protons
    if(fabs(type)==321) return true; //charged kaons

  
  } else {
    // filter by vetoing the ones that are bad
    pass=true;
    if(fabs(type)==12 || fabs(type)==14 || fabs(type)==16) return false;
    if(fabs(type)==17) return false; //tau
    if(fabs(type)<=10) return false; //quarks itrkCurv
    if(fabs(type)==130 || fabs(type)==310 || fabs(type)==311) return false; //neutral kaons
    if(fabs(type)==111) return false; //pi zero 
    if(fabs(type)==2112) return false; //neutrons
    if(fabs(type)==3212 || fabs(type)==3122 || fabs(type)==3322 ) return false; //sigma zero and lambda zero
    if(fabs(type)==3222 || fabs(type)==3112 || fabs(type)==3312) return false; //charged sigma
    if(fabs(type)>20 && fabs(type)<100) return false; //bosons plus exotics
    if(fabs(type)==521 || fabs(type)==511) return false; //b mesons
    if(fabs(type)==421 || fabs(type)==411) return false; //c mesons

  }

  return pass;

}

void loadFPGATree(TString fileName = "myTest.root") {
 
  TFile *tFile = new TFile(fileName.Data());
  fpgaEvent = 0;
  fpgaTree = (TTree*)tFile->Get("FPGAEvent");

  fpgaTree->SetBranchAddress("Event", &fpgaEvent);
  numEvts = (int)fpgaTree->GetEntries();
 
  printf("Number of Entries %i\n",numEvts);

}

void analyzeEvent(TString histFileName="test.root", TString rootFileName="myTest.root") {

  const bool stubsTxtFile = false;
  const double Maxpt = 200.0;
  const double Minpt = 0.0;

  const int minStubs = 4;
  const int minIndStubs = 3;

  ofstream stubid ("StubIDMulti.txt");

  if(fpgaEvent==0) loadFPGATree(rootFileName);

  //Setup histogram locations
  TFile *histFile = new TFile(histFileName,"RECREATE");
  histFile->cd();
  TDirectory *barrelOnly = histFile->mkdir("barrelOnly");
  TDirectory *diskOnly = histFile->mkdir("diskOnly");
  TDirectory *overlap = histFile->mkdir("overlap");
  histFile->cd();

  //Initialize histograms NOT separated by region

  TH1F* h_mcTrkPt_tot  = new TH1F("h_mcTrkPt_tot", "mc trk Pt",200,0.,200.);
  TH1F* h_mcTrkPt_low_tot  = new TH1F("h_mcTrkPt_low_tot", "mc trk Pt",50,0.,10.);
  TH1F* h_mcTrkPhi_tot = new TH1F("h_mcTrkPhi_tot","mc trk Phi",200,-TMath::Pi(),TMath::Pi());
  TH1F* h_mcTrkZ0_tot  = new TH1F("h_mcTrkZ0_tot", "mc trk Z0",200,-50.,50.);
  TH1F* h_mcTrkVr_tot  = new TH1F("h_mcTrkVr_tot", "mc sq(vx^2+vy^2)",2000,0.,2.);
  TH1F* h_mcTrkEta_tot = new TH1F("h_mcTrkEta_tot","mc trk Eta",100,-5.,5.);
  TH1F* h_mcType_tot   = new TH1F("h_mcType_tot",  "mc type",801,-400.5,400.5);

  TH1F* h_trkPtWODup_tot   = new TH1F("h_trkPtWODup_tot","trk pt",100,0.,25.);
  TH1F* h_trkPtWODup_low_tot   = new TH1F("h_trkPtWODup_low_tot","trk pt",50,0.,10.);
  TH1F* h_trkPhiWODup_tot  = new TH1F("h_trkPhiWODup_tot","trk (global) phi",100,0.,TMath::TwoPi());
  TH1F* h_trkEtaWODup_tot  = new TH1F("h_trkEtaWODup_tot","trk eta",100,-4.,4.);
  TH1F* h_trkZ0WODup_tot   = new TH1F("h_trkZ0WODup_tot","trk z0",200,-50.,50.);

  TH1F* h_trkPtMulti_tot     = new TH1F("h_trkPtMulti_tot","multiple trk pt",200,-25,25);
  TH1F* h_trkPhiMulti_tot    = new TH1F("h_trkPhiMulti_tot","multiple trk phi",100,0,TMath::TwoPi());
  TH1F* h_trkEtaMulti_tot    = new TH1F("h_trkEta_Multi_tot","multiple trk Eta",100,-5.,5.);
  TH1F* h_trkZ0Multi_tot     = new TH1F("h_trkZ0_Multi_tot","multiple trk Z0",200,-50.,50.);
  TH1F* h_itrkPhiMulti_tot   = new TH1F("h_itrkPhiMulti_tot","multiple trk iphi",1000,-50000.,120000.);
  TH1F* h_trkPhiMultiLayered_tot = new TH1F("h_trkPhiMultiLayered_tot","multiple trk phi layered",1000,0,0.2243995);

  TH1F* h_nBarrelTrkOverlap       = new TH1F("h_nBarrelTrkOverlap","Num of barrel trks",12,0.5,12.5);
  TH1F* h_nBarrelTrkOverlapWODup  = new TH1F("h_nBarrelTrkOverlapWODup","Num of barrel trks WO dup",12,0.5,12.5);
  TH1F* h_nDiskTrkOverlap         = new TH1F("h_nDiskTrkOverlap","Num of disk trks",12,0.5,12.5);
  TH1F* h_nDiskTrkOverlapWODup    = new TH1F("h_nDiskTrkOverlapWODup","Num of disk trks WO dup",12,0.5,12.5);
  TH1F* h_nOverlapTrkOverlap      = new TH1F("h_nOverlapTrkOverlap","Num of overlap trks",12,0.5,12.5);
  TH1F* h_nOverlapTrkOverlapWODup = new TH1F("h_nOverlapTrkOverlapWODup","Num of overlap trks WO dup",12,0.5,12.5);

  TH1F* h_nTrkSec_eta1      = new TH1F("h_nTrkSec_eta1","Num Tracks per Sector",26,-0.5,25.5);
  TH1F* h_nTrkSec_eta2      = new TH1F("h_nTrkSec_eta2","Num Tracks per Sector",26,-0.5,25.5);
  TH1F* h_nTrkSec_eta3      = new TH1F("h_nTrkSec_eta3","Num Tracks per Sector",26,-0.5,25.5);
  TH1F* h_nTrkSec_eta4      = new TH1F("h_nTrkSec_eta4","Num Tracks per Sector",26,-0.5,25.5);

  TH1F* h_chisq_afterPD     = new TH1F("h_chisq_afterPD","chisq of tracks",100,0,7);
  TH1F* h_chisq_beforePD    = new TH1F("h_chisq_beforePD","chisq of tracks",100,0,7);
    
  //For Single Muons
//  TH1F* h_nMCTrkEvt_tot    = new TH1F("h_nMCTrkEvt_tot","Num MC Tracks per Event",15,0.5,15.5);
//  TH1F* h_nTrkEvt_tot      = new TH1F("h_nTrkEvt_tot","Num Tracks per Event",15,0.5,15.5);
//  TH1F* h_nTrkEvtWODup_tot = new TH1F("h_nTrkEvtWODup_tot","Num Tracks per Event after duplicate removal",10,0.5,10.5);
//  TH1F* h_nDupTrk_tot      = new TH1F("h_nDupTrk_tot","Num of duplicates removed",10,0.5,10.5);
//  TH1F* h_nMCTrkSec_tot    = new TH1F("h_nMCTrkSec_tot","Num MC Tracks per Sector",15,0.5,15.5);
//  TH1F* h_nTrkSec_tot      = new TH1F("h_nTrkSec_tot","Num Tracks per Sector",15,0.5,15.5);
//  TH1F* h_nTrkSecTail_tot  = new TH1F("h_nTrkSecTail_tot","Num Tracks per Sector",15,15.5,30.5);
//  TH1F* h_nTrkSecWODup_tot = new TH1F("h_nTrkSecWODup_tot","Num Tracks per Sector (after dup)",15,0.5,15.5);
//  TH1F* h_nTrkSecWODupTail_tot = new TH1F("h_nTrkSecWODupTail_tot","Num Tracks per Sector (after dup)",15,15.5,30.5);
//  TH2F* h_NDupVSNTrkWO_tot = new TH2F("h_NDupVSNTrkWO_tot","Num Dup vs Num Trk (after)",10,0.5,10.5,10,0.5,10.5);
    

  //For ttbar
  TH1F* h_nMCTrkEvt_tot    = new TH1F("h_nMCTrkEvt_tot","Num MC Tracks per Event",51,-0.5,200.5);
  TH1F* h_nTrkEvt_tot      = new TH1F("h_nTrkEvt_tot","Num Tracks per Event",51,-0.5,400.5);
  TH1F* h_nTrkEvtWODup_tot = new TH1F("h_nTrkEvtWODup_tot","Num Tracks per Event after duplicate removal",51,-0.5,300.5);
  TH1F* h_nDupTrk_tot      = new TH1F("h_nDupTrk_tot","Num of duplicates removed",51,-0.5,200.5);
  TH1F* h_nMCTrkSec_tot    = new TH1F("h_nMCTrkSec_tot","Num MC Tracks per Sector",51,-0.5,50.5);
  TH1F* h_nTrkSec_tot      = new TH1F("h_nTrkSec_tot","Num Tracks per Sector",51,-0.5,50.5);
  TH1F* h_nTrkSecTail_tot  = new TH1F("h_nTrkSecTail_tot","Num Tracks per Sector",51,49.5,100.5);
  TH1F* h_nTrkSecWODup_tot = new TH1F("h_nTrkSecWODup_tot","Num Tracks per Sector (after dup)",51,-0.5,50.5);
  TH1F* h_nTrkSecWODupTail_tot = new TH1F("h_nTrkSecWODupTail_tot","Num Tracks per Sector (after dup)",51,49.5,100.5);
  TH2F* h_NDupVSNTrkWO_tot = new TH2F("h_NDupVSNTrkWO_tot","Num Dup vs Num Trk (after)",10,0.5,10.5,10,0.5,10.5);


  TH1F* h_nStub = new TH1F("h_nStub","nStub per Event",31,-0.5,30.5);
  TH1F* h_nStubTrk     = new TH1F("h_nStubTrk","Num Stubs per Track",9,-0.5,8.5);
  TH1F* h_nStub_ghost = new TH1F("h_nStub_ghost","nStub per Event with ghosts",31,-0.5,30.5);
  TH1F* h_iLTkt_tot = new TH1F("h_iLTkt_tot","Track seed",45,-22.5,22.5);
  TH1F* h_iLTktWODup_tot = new TH1F("h_iLTktWODup_tot","Track seed W/O duplicates",45,-22.5,22.5);
  TH1F* h_iLTkt_ghost = new TH1F("h_iLTkt_ghost","Track seed in ghost events",45,-22.5,22.5);

  TH1F* h_matchDeltaPhi_AllM_tot   = new TH1F("h_matchDeltaPhi_AllM_tot","fpga-mc delta phi (all)",100,-0.003,0.003);
  TH1F* h_matchDeltaEta_AllM_tot   = new TH1F("h_matchDeltaEta_AllM_tot","fpga-mc delta eta (all)",100,-0.025,0.025);
  TH1F* h_matchDeltaR_AllM_tot     = new TH1F("h_matchDeltaR_AllM_tot"  ,"fpga-mc delta R (all)"  ,100,0.,0.025);
  TH1F* h_matchDeltaZ0_AllM_tot    = new TH1F("h_matchDeltaZ0_AllM_tot"  ,"fpga-mc delta Z0 (all)"  ,100,-1.5,1.5);
  TH1F* h_matchDeltaPtOPt_AllM_tot = new TH1F("h_matchDeltaPtOPt_AllM_tot" ,"fpga-mc (delta Pt)/Pt (all)" ,100,-0.1,0.1);

  TH2F* h_matchDeltaPhi_chi2_AllM   = new TH2F("h_matchDeltaPhi_chi2_AllM","fpga-mc delta phi v. chi2",40,-0.003,0.003,40,0,20);
  TH2F* h_matchDeltaEta_chi2_AllM   = new TH2F("h_matchDeltaEta_chi2_AllM","fpga-mc delta eta v. chi2",40,-0.025,0.025,40,0,20);
  TH2F* h_matchDeltaR_chi2_AllM     = new TH2F("h_matchDeltaR_chi2_AllM"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025,40,0,20);
  TH2F* h_matchDeltaZ0_chi2_AllM    = new TH2F("h_matchDeltaZ0_chi2_AllM"  ,"fpga-mc delta Z0 v. chi2"  ,40,-1.5,1.5,40,0,20);
  TH2F* h_matchDeltaPtOPt_chi2_AllM = new TH2F("h_matchDeltaPtOPt_chi2_AllM" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1,40,0,20);

  TH1F* h_matchDeltaPhi_chi2_0to1   = new TH1F("h_matchDeltaPhi_chi2_0to1","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_0to1   = new TH1F("h_matchDeltaEta_chi2_0to1","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_0to1     = new TH1F("h_matchDeltaR_chi2_0to1"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_0to1    = new TH1F("h_matchDeltaZ0_chi2_0to1"  ,"fpga-mc delta Z0 v. chi2"  ,40,-1.5,1.5);
  TH1F* h_matchDeltaPtOPt_chi2_0to1 = new TH1F("h_matchDeltaPtOPt_chi2_0to1" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_matchDeltaPhi_chi2_1to2   = new TH1F("h_matchDeltaPhi_chi2_1to2","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_1to2   = new TH1F("h_matchDeltaEta_chi2_1to2","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_1to2     = new TH1F("h_matchDeltaR_chi2_1to2"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_1to2    = new TH1F("h_matchDeltaZ0_chi2_1to2"  ,"fpga-mc delta Z0 v. chi2"  ,40,-1.5,1.5);
  TH1F* h_matchDeltaPtOPt_chi2_1to2 = new TH1F("h_matchDeltaPtOPt_chi2_1to2" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_matchDeltaPhi_chi2_2to3   = new TH1F("h_matchDeltaPhi_chi2_2to3","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_2to3   = new TH1F("h_matchDeltaEta_chi2_2to3","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_2to3     = new TH1F("h_matchDeltaR_chi2_2to3"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_2to3    = new TH1F("h_matchDeltaZ0_chi2_2to3"  ,"fpga-mc delta Z0 v. chi2"  ,40,-1.5,1.5);
  TH1F* h_matchDeltaPtOPt_chi2_2to3 = new TH1F("h_matchDeltaPtOPt_chi2_2to3" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_matchDeltaPhi_chi2_3to10   = new TH1F("h_matchDeltaPhi_chi2_3to10","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_3to10   = new TH1F("h_matchDeltaEta_chi2_3to10","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_3to10     = new TH1F("h_matchDeltaR_chi2_3to10"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_3to10    = new TH1F("h_matchDeltaZ0_chi2_3to10"  ,"fpga-mc delta Z0 v. chi2"  ,40,-1.5,1.5);
  TH1F* h_matchDeltaPtOPt_chi2_3to10 = new TH1F("h_matchDeltaPtOPt_chi2_3to10" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH2F* h_nTrkVsmcTrkPt_tot      = new TH2F("h_nTrkVsmcTrkPt_tot", "nTrk vs mc trk Pt",600,0.,150.,10,-0.5,9.5);
  TH2F* h_nTrkWODupVsmcTrkPt_tot = new TH2F("h_nTrkWODupVsmcTrkPt_tot", "NTrk vs mc trk Pt (after)",600,0.,150.,10,-0.5,9.5);

  TH1F* h_L1Multi = new TH1F("h_L1Multi", "nL1Stubs",10,0.5,10.5);
  TH1F* h_L2Multi = new TH1F("h_L2Multi", "nL2Stubs",10,0.5,10.5);
  TH1F* h_L3Multi = new TH1F("h_L3Multi", "nL3Stubs",10,0.5,10.5);
  TH1F* h_L4Multi = new TH1F("h_L4Multi", "nL4Stubs",10,0.5,10.5);
  TH1F* h_L5Multi = new TH1F("h_L5Multi", "nL5Stubs",10,0.5,10.5);
  TH1F* h_L6Multi = new TH1F("h_L6Multi", "nL6Stubs",10,0.5,10.5);
  TH1F* h_F1Multi = new TH1F("h_F1Multi", "nF1Stubs",10,0.5,10.5);
  TH1F* h_F2Multi = new TH1F("h_F2Multi", "nF2Stubs",10,0.5,10.5);
  TH1F* h_F3Multi = new TH1F("h_F3Multi", "nF3Stubs",10,0.5,10.5);
  TH1F* h_F4Multi = new TH1F("h_F4Multi", "nF4Stubs",10,0.5,10.5);
  TH1F* h_F5Multi = new TH1F("h_F5Multi", "nF5Stubs",10,0.5,10.5);
  TH1F* h_B1Multi = new TH1F("h_B1Multi", "nB1Stubs",10,0.5,10.5);
  TH1F* h_B2Multi = new TH1F("h_B2Multi", "nB2Stubs",10,0.5,10.5);
  TH1F* h_B3Multi = new TH1F("h_B3Multi", "nB3Stubs",10,0.5,10.5);
  TH1F* h_B4Multi = new TH1F("h_B4Multi", "nB4Stubs",10,0.5,10.5);
  TH1F* h_B5Multi = new TH1F("h_B5Multi", "nB5Stubs",10,0.5,10.5);

  TH1F* h_L1Multi_ghost = new TH1F("h_L1Multi_ghost", "nL1Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_L2Multi_ghost = new TH1F("h_L2Multi_ghost", "nL2Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_L3Multi_ghost = new TH1F("h_L3Multi_ghost", "nL3Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_L4Multi_ghost = new TH1F("h_L4Multi_ghost", "nL4Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_L5Multi_ghost = new TH1F("h_L5Multi_ghost", "nL5Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_L6Multi_ghost = new TH1F("h_L6Multi_ghost", "nL6Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_F1Multi_ghost = new TH1F("h_F1Multi_ghost", "nF1Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_F2Multi_ghost = new TH1F("h_F2Multi_ghost", "nF2Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_F3Multi_ghost = new TH1F("h_F3Multi_ghost", "nF3Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_F4Multi_ghost = new TH1F("h_F4Multi_ghost", "nF4Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_F5Multi_ghost = new TH1F("h_F5Multi_ghost", "nF5Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_B1Multi_ghost = new TH1F("h_B1Multi_ghost", "nB1Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_B2Multi_ghost = new TH1F("h_B2Multi_ghost", "nB2Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_B3Multi_ghost = new TH1F("h_B3Multi_ghost", "nB3Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_B4Multi_ghost = new TH1F("h_B4Multi_ghost", "nB4Stubs w/ghosts",10,0.5,10.5);
  TH1F* h_B5Multi_ghost = new TH1F("h_B5Multi_ghost", "nB5Stubs w/ghosts",10,0.5,10.5);

  TH1F* h_L1L2Multi = new TH1F("h_L1L2Multi", "n L1L2 tracks per event",10,0.5,10.5);
  TH1F* h_L3L4Multi = new TH1F("h_L3L4Multi", "n L3L4 tracks per event",10,0.5,10.5);
  TH1F* h_L5L6Multi = new TH1F("h_L5L6Multi", "n L5L6 tracks per event",10,0.5,10.5);
  TH1F* h_B1B2Multi = new TH1F("h_B1B2Multi", "n B1B2 tracks per event",10,0.5,10.5);
  TH1F* h_B3B4Multi = new TH1F("h_B3B4Multi", "n B3B4 tracks per event",10,0.5,10.5);
  TH1F* h_F1F2Multi = new TH1F("h_F1F2Multi", "n F1F2 tracks per event",10,0.5,10.5);
  TH1F* h_F3F4Multi = new TH1F("h_F3F4Multi", "n F3F4 tracks per event",10,0.5,10.5);
  TH1F* h_B1LMulti = new TH1F("h_B1LMulti", "n L1B1 tracks per event",10,0.5,10.5);
  TH1F* h_F1LMulti = new TH1F("h_F1LMulti", "n L1F1 tracks per event",10,0.5,10.5);

  TH1F* h_L1L2Multi_ghost = new TH1F("h_L1L2Multi_ghost", "n L1L2 tracks per event",10,0.5,10.5);
  TH1F* h_L3L4Multi_ghost = new TH1F("h_L3L4Multi_ghost", "n L3L4 tracks per event",10,0.5,10.5);
  TH1F* h_L5L6Multi_ghost = new TH1F("h_L5L6Multi_ghost", "n L5L6 tracks per event",10,0.5,10.5);
  TH1F* h_B1B2Multi_ghost = new TH1F("h_B1B2Multi_ghost", "n B1B2 tracks per event",10,0.5,10.5);
  TH1F* h_B3B4Multi_ghost = new TH1F("h_B3B4Multi_ghost", "n B3B4 tracks per event",10,0.5,10.5);
  TH1F* h_F1F2Multi_ghost = new TH1F("h_F1F2Multi_ghost", "n F1F2 tracks per event",10,0.5,10.5);
  TH1F* h_F3F4Multi_ghost = new TH1F("h_F3F4Multi_ghost", "n F3F4 tracks per event",10,0.5,10.5);
  TH1F* h_B1LMulti_ghost = new TH1F("h_B1LMulti_ghost", "n L1B1 tracks per event",10,0.5,10.5);
  TH1F* h_F1LMulti_ghost = new TH1F("h_F1LMulti_ghost", "n L1F1 tracks per event",10,0.5,10.5);

  TH2F* h_IDvsiZ = new TH2F("h_IDvsiZ", "Stub ID vs iZ",200,-200,200,64,-0.5,512);
  TH2F* h_IDvsZ = new TH2F("h_IDvsZ", "Stub ID vs Z",110,-110,200,64,-0.5,512);
  TH2F* h_IDvsLayer = new TH2F("h_IDvsLayer", "Stub ID vs Layer",41,-20.5,20.5,64,-0.5,512);
  
  // For Resolution vs. nStub by region

  TH1F* h_matchDeltaPhi_trklet_nS4[3];
  TH1F* h_matchDeltaEta_trklet_nS4[3];
  TH1F* h_matchDeltaR_trklet_nS4[3];
  TH1F* h_matchDeltaZ0_trklet_nS4[3];
  TH1F* h_matchDeltaPtOPt_trklet_nS4[3];

  TH1F* h_matchDeltaPhi_trklet_nS5[3];
  TH1F* h_matchDeltaEta_trklet_nS5[3];
  TH1F* h_matchDeltaR_trklet_nS5[3];
  TH1F* h_matchDeltaZ0_trklet_nS5[3];
  TH1F* h_matchDeltaPtOPt_trklet_nS5[3];

  TH1F* h_matchDeltaPhi_trklet_nS6[3];
  TH1F* h_matchDeltaEta_trklet_nS6[3];
  TH1F* h_matchDeltaR_trklet_nS6[3];
  TH1F* h_matchDeltaZ0_trklet_nS6[3];
  TH1F* h_matchDeltaPtOPt_trklet_nS6[3];

  for(int r=0; r<3; r++) {

    if(r==0) barrelOnly->cd();
    if(r==1) diskOnly->cd();
    if(r==2) overlap->cd();

    h_matchDeltaPhi_trklet_nS4[r]   = new TH1F("h_matchDeltaPhi_nS4","fpga-mc delta phi (nstub=4)",100,-0.003,0.003);
    h_matchDeltaEta_trklet_nS4[r]   = new TH1F("h_matchDeltaEta_nS4","fpga-mc delta eta (nstub=4)",100,-0.025,0.025);
    h_matchDeltaR_trklet_nS4[r]     = new TH1F("h_matchDeltaR_nS4"  ,"fpga-mc delta R (nstub=4)"  ,100,0.,0.025);
    h_matchDeltaZ0_trklet_nS4[r]    = new TH1F("h_matchDeltaZ0_nS4"  ,"fpga-mc delta Z0 (nstub=4)"  ,100,-1.5,1.5);
    h_matchDeltaPtOPt_trklet_nS4[r] = new TH1F("h_matchDeltaPtOPt_nS4" ,"fpga-mc (delta Pt)/Pt (nstub=4)" ,100,-0.1,0.1);

    h_matchDeltaPhi_trklet_nS5[r]   = new TH1F("h_matchDeltaPhi_nS5","fpga-mc delta phi (nstub=5)",100,-0.003,0.003);
    h_matchDeltaEta_trklet_nS5[r]   = new TH1F("h_matchDeltaEta_nS5","fpga-mc delta eta (nstub=5)",100,-0.025,0.025);
    h_matchDeltaR_trklet_nS5[r]     = new TH1F("h_matchDeltaR_nS5"  ,"fpga-mc delta R (nstub=5)"  ,100,0.,0.025);
    h_matchDeltaZ0_trklet_nS5[r]    = new TH1F("h_matchDeltaZ0_nS5"  ,"fpga-mc delta Z0 (nstub=5)"  ,100,-1.5,1.5);
    h_matchDeltaPtOPt_trklet_nS5[r] = new TH1F("h_matchDeltaPtOPt_nS5" ,"fpga-mc (delta Pt)/Pt (nstub=5)" ,100,-0.1,0.1);

    h_matchDeltaPhi_trklet_nS6[r]   = new TH1F("h_matchDeltaPhi_nS6","fpga-mc delta phi (nstub=6)",100,-0.003,0.003);
    h_matchDeltaEta_trklet_nS6[r]   = new TH1F("h_matchDeltaEta_nS6","fpga-mc delta eta (nstub=6)",100,-0.025,0.025);
    h_matchDeltaR_trklet_nS6[r]     = new TH1F("h_matchDeltaR_nS6"  ,"fpga-mc delta R (nstub=6)"  ,100,0.,0.025);
    h_matchDeltaZ0_trklet_nS6[r]    = new TH1F("h_matchDeltaZ0_nS6"  ,"fpga-mc delta Z0 (nstub=6)"  ,100,-1.5,1.5);
    h_matchDeltaPtOPt_trklet_nS6[r] = new TH1F("h_matchDeltaPtOPt_nS6" ,"fpga-mc (delta Pt)/Pt (nstub=6)" ,100,-0.1,0.1);
  }

  // Loop over events
  for(int i=0; i<numEvts; i++) {
//    cout << "Analyzing event  " << i << "..." << endl;
    fpgaTree->GetEntry(i);

    // Initialize map to count stubs
    std::map<int,int> nStubEvt;
    int nStubTot=0;
    for(int iLayers=1; iLayers<=6; iLayers++) nStubEvt[iLayers]=0;
    for(int iDisks=1; iDisks<=5; iDisks++) {
      nStubEvt[10+iDisks]=0;
      nStubEvt[-10-iDisks]=0;
    }

    // Plot stubs
    for(int iStub=0; iStub<(int)fpgaEvent->stubs.size(); iStub++) {
      FPGAEventStub stub = fpgaEvent->stubs.at(iStub);
      h_IDvsiZ->Fill(stub.iz_,stub.stubID_);
      h_IDvsZ->Fill(stub.z_,stub.stubID_);
      h_IDvsLayer->Fill(min(stub.layer_,stub.disk_),stub.stubID_);
      nStubEvt[min(stub.layer_,stub.disk_)]++;
      nStubTot++;
    }
    h_nStub->Fill(nStubTot);
    h_L1Multi->Fill(nStubEvt[1]);
    h_L2Multi->Fill(nStubEvt[2]);
    h_L3Multi->Fill(nStubEvt[3]);
    h_L4Multi->Fill(nStubEvt[4]);
    h_L5Multi->Fill(nStubEvt[5]);
    h_L6Multi->Fill(nStubEvt[6]);
    h_F1Multi->Fill(nStubEvt[10+1]);
    h_F2Multi->Fill(nStubEvt[10+2]);
    h_F3Multi->Fill(nStubEvt[10+3]);
    h_F4Multi->Fill(nStubEvt[10+4]);
    h_F5Multi->Fill(nStubEvt[10+5]);
    h_B1Multi->Fill(nStubEvt[-10-1]);
    h_B2Multi->Fill(nStubEvt[-10-2]);
    h_B3Multi->Fill(nStubEvt[-10-3]);
    h_B4Multi->Fill(nStubEvt[-10-4]);
    h_B5Multi->Fill(nStubEvt[-10-5]);
    
    //Loop through events
    int nDupTrk_tot=0;

    int nTrkSec_tot[28]={0};
    int nTrkSec_eta[28][4]={{0}};
    int nTrkSecWODup_tot[28]={0};
    int nTrkWODup_tot=0;

    std::map<int,int> nILTkt;
    nILTkt[1] = 0; nILTkt[3] = 0; nILTkt[5] = 0;
    nILTkt[-11] = 0; nILTkt[-13] = 0; nILTkt[11] = 0; nILTkt[13] = 0;
    nILTkt[-21] = 0; nILTkt[21] = 0;

    for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
      FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
      if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt) {

        nILTkt[track.seed_]++;
        h_iLTkt_tot->Fill(track.seed_);

        h_chisq_beforePD->Fill(track.chisq_);
        h_nStubTrk->Fill(track.stubID_.size());

        nTrkSec_tot[track.sector_]++;

        if( track.eta_<=-1.0 ) nTrkSec_eta[track.sector_][0]++;
        if( track.eta_>-1.0 && track.eta_<=0) nTrkSec_eta[track.sector_][1]++;
        if( track.eta_<1.0 && track.eta_>0) nTrkSec_eta[track.sector_][2]++;
        if( track.eta_>=1.0 ) nTrkSec_eta[track.sector_][3]++;

        // Catch 7-stub tracks
        if(track.stubID_.size()>6) {
          cout << "Track with 7 stubs! Seed =" << track.seed_ << endl << "layers = ";
          for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
            cout << st->first << ", ";
          }
          cout << endl;
        }

        // Track dump -- matches fpga.cc
//          printf("Track Parameters: \n");
//          printf("irinv = %i \n", track.irinv_ );
//          printf("iphi0 = %i \n", track.iphi0_ );
//          printf("iz0   = %i \n", track.iz0_ );
//          printf("it    = %i \n", track.it_ );
//          printf("stubID=");
//          std::map<int, int> stubs = track.stubID_;
//          for(std::map<int, int>::iterator sb=stubs.begin(); sb!=stubs.end(); sb++) printf(" %i -- %i ",sb->first,sb->second);
//          printf("\n");
//          printf("dup   = %i\n \n", track.duplicate_);
//		      printf("chisq   = %f\n \n", track.chisq_);

        if(track.duplicate_) {
          nDupTrk_tot++;
        } else {
          //if track isn't a duplicate
          nTrkWODup_tot++;
          nTrkSecWODup_tot[track.sector_]++;

          h_iLTktWODup_tot->Fill(track.seed_);

          h_trkPtWODup_tot->Fill(abs(track.pt_));
          h_trkPtWODup_low_tot->Fill(abs(track.pt_));
          h_trkPhiWODup_tot->Fill(track.phi0_);
          h_trkEtaWODup_tot->Fill(track.eta_);
          h_trkZ0WODup_tot->Fill(track.z0_);

          h_chisq_afterPD->Fill(track.chisq_);

        }//end duplicates loop
      }//end if nStub >= nStubs
    }//end loop over tracks

    //make event level plots, but suppress events with no tracks passing requirements
    if(nTrkWODup_tot>0) {

      //total plots
      h_nTrkEvt_tot->Fill(nTrkWODup_tot+nDupTrk_tot);
      h_nTrkEvtWODup_tot->Fill(nTrkWODup_tot);
      h_nDupTrk_tot->Fill(nDupTrk_tot);
      h_NDupVSNTrkWO_tot->Fill(nTrkWODup_tot,nDupTrk_tot);

      h_L1L2Multi->Fill(nILTkt[1]);
      h_L3L4Multi->Fill(nILTkt[3]);
      h_L5L6Multi->Fill(nILTkt[5]);
      h_B1B2Multi->Fill(nILTkt[-11]);
      h_B3B4Multi->Fill(nILTkt[-13]);
      h_F1F2Multi->Fill(nILTkt[11]);
      h_F3F4Multi->Fill(nILTkt[13]);
      h_B1LMulti->Fill(nILTkt[-21]);
      h_F1LMulti->Fill(nILTkt[21]);
            
      //make plots for events with multiple unique tracks, output track information to .txt file
      if(nTrkWODup_tot>=2) {
        //Stub layer plots
        h_nStub_ghost->Fill(nStubTot);
        h_L1Multi_ghost->Fill(nStubEvt[1]);
        h_L2Multi_ghost->Fill(nStubEvt[2]);
        h_L3Multi_ghost->Fill(nStubEvt[3]);
        h_L4Multi_ghost->Fill(nStubEvt[4]);
        h_L5Multi_ghost->Fill(nStubEvt[5]);
        h_L6Multi_ghost->Fill(nStubEvt[6]);
        h_F1Multi_ghost->Fill(nStubEvt[10+1]);
        h_F2Multi_ghost->Fill(nStubEvt[10+2]);
        h_F3Multi_ghost->Fill(nStubEvt[10+3]);
        h_F4Multi_ghost->Fill(nStubEvt[10+4]);
        h_F5Multi_ghost->Fill(nStubEvt[10+5]);
        h_B1Multi_ghost->Fill(nStubEvt[-10-1]);
        h_B2Multi_ghost->Fill(nStubEvt[-10-2]);
        h_B3Multi_ghost->Fill(nStubEvt[-10-3]);
        h_B4Multi_ghost->Fill(nStubEvt[-10-4]);
        h_B5Multi_ghost->Fill(nStubEvt[-10-5]);

        h_L1L2Multi_ghost->Fill(nILTkt[1]);
        h_L3L4Multi_ghost->Fill(nILTkt[3]);
        h_L5L6Multi_ghost->Fill(nILTkt[5]);
        h_F1F2Multi_ghost->Fill(nILTkt[11]);
        h_F3F4Multi_ghost->Fill(nILTkt[13]);
        h_B1B2Multi_ghost->Fill(nILTkt[-11]);
        h_B3B4Multi_ghost->Fill(nILTkt[-13]);
        h_F1LMulti_ghost->Fill(nILTkt[21]);
        h_B1LMulti_ghost->Fill(nILTkt[-21]);

        stubid << Form("Event %i:",i);
        for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++) {
          FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
          if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt && !track.duplicate_) {

            h_iLTkt_ghost->Fill(track.seed_);

            //Plot simple variables
            h_trkPhiMulti_tot->Fill(track.phi0_);
            h_trkEtaMulti_tot->Fill(track.eta_);
            h_trkZ0Multi_tot->Fill(track.z0_);
            h_itrkPhiMulti_tot->Fill(track.iphi0_);
            h_trkPtMulti_tot->Fill(abs(track.pt_));

            //Make a superposition of the sectors
            if(track.phi0_<0) {
              h_trkPhiMultiLayered_tot->Fill(track.phi0_+TMath::TwoPi()/28);
            } else {
              int nPhiSec=(track.phi0_/(TMath::TwoPi()/28));
              h_trkPhiMultiLayered_tot->Fill(track.phi0_-nPhiSec*TMath::TwoPi()/28);
            }
          }
        }

        // Write out track information to .txt file
        if(stubsTxtFile) {
          for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
            FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
            if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt && !track.duplicate_) {
              //Output stubids to .txt file
              int nShareGhost=0;
              stubid << Form("\nTrack %i: ",ntrk);
              stubid << Form(" Phi=%f, Eta=%f, Z0=%f, Pt=%f, iLTkt=%i, Sector=%i", track.phi0_,track.eta_,track.z0_,abs(track.pt_),track.seed_,track.sector_);
              stubid << Form("\n    Stub ids: ");
              for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
                stubid << Form("  %i,%i",st->first,st->second);
                for(int ntrk2=0; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
                  FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);
                  if(!track2.duplicate_ && track2.stubID_.size() >= minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && ntrk!=ntrk2) {
                    if(track2.stubID_.find(st->first) != track2.stubID_.end()) {
                      if(st->second == track2.stubID_[st->first] && st->second != 63) {
                        stubid << "*";
                        nShareGhost++;
                      }
                    }
                  }//end if duplicate
                }//end loop over secondary tracks
              }//end loop over stubs
              if (nTrkWODup_tot==2)stubid << Form("\n    %i Shared, %i Independent", nShareGhost, (int)track.stubID_.size()-nShareGhost);
            }//end if duplicate
          }//end loop over tracks
          stubid << "\n\n";
        }//end if stubsTxtFile

      }//end if nTrkWODup>=2
    }//end if nTrkWODup>0


    // Now look at the monte carlo distributions
    // only look at the events where we have a good FPGA track

    int nMCTrkEvt_tot=0;
    int nMCTrkSec_tot[28]={0};
    for(int ntrk=0; ntrk<(int)fpgaEvent->mcTracks.size(); ntrk++){

      FPGAEventMCTrack mcTrack = fpgaEvent->mcTracks.at(ntrk);

      //apply some track requirement
      if(abs(mcTrack.pt_)>Minpt && abs(mcTrack.pt_)<Maxpt
         && sqrt(mcTrack.vy_*mcTrack.vy_ + mcTrack.vx_*mcTrack.vx_)<0.1  
         && checkType(mcTrack.type_,true)
         && fabs(mcTrack.eta_)<2.5
         && fabs(mcTrack.vz_)<20
         ) {

        //if(!checkType(mcTrack.type_,false)) printf("Type neither good nor bad %i\n",mcTrack.type_);
        h_mcTrkPt_tot->Fill(abs(mcTrack.pt_));
        h_mcTrkPt_low_tot->Fill(abs(mcTrack.pt_));
        h_mcTrkPhi_tot->Fill(mcTrack.phi_); 
        h_mcTrkZ0_tot->Fill(mcTrack.vz_);
        h_mcTrkVr_tot->Fill(sqrt(mcTrack.vy_*mcTrack.vy_ + mcTrack.vx_*mcTrack.vx_));   
        h_mcTrkEta_tot->Fill(mcTrack.eta_);
        h_mcType_tot->Fill(mcTrack.type_);

        //Tally number of tracks per event and sector (artifically)
        nMCTrkEvt_tot++;
        nMCTrkSec_tot[(int)((mcTrack.phi_+TMath::Pi())/(TMath::TwoPi()/28))]++;

        // Try to see if this particle can be matched to a found track
        FPGAEventTrackMatch *bestMatch = 0;
        for(int fntrk=0; fntrk<(int)fpgaEvent->tracks.size(); fntrk++){
          FPGAEventTrack fpgaTrack = fpgaEvent->tracks.at(fntrk);

          if(fpgaTrack.stubID_.size() >= minStubs && !fpgaTrack.duplicate_) {
            FPGAEventTrackMatch *match = new FPGAEventTrackMatch(mcTrack,fpgaTrack);

            h_matchDeltaPhi_AllM_tot->Fill(match->deltaPhi_);
            h_matchDeltaEta_AllM_tot->Fill(match->deltaEta_);
            h_matchDeltaR_AllM_tot  ->Fill(match->deltaR_);
            h_matchDeltaZ0_AllM_tot ->Fill(match->deltaZ0_);
            h_matchDeltaPtOPt_AllM_tot ->Fill(match->deltaPt_/abs(mcTrack.pt_));

            h_matchDeltaPhi_chi2_AllM->Fill(match->deltaPhi_,fpgaTrack.chisq_);
            h_matchDeltaEta_chi2_AllM->Fill(match->deltaEta_,fpgaTrack.chisq_);
            h_matchDeltaR_chi2_AllM->Fill(match->deltaR_,fpgaTrack.chisq_);
            h_matchDeltaZ0_chi2_AllM->Fill(match->deltaZ0_,fpgaTrack.chisq_);
            h_matchDeltaPtOPt_chi2_AllM->Fill(match->deltaPt_/abs(mcTrack.pt_),fpgaTrack.chisq_);

            if(fpgaTrack.chisq_ >= 0.0 && fpgaTrack.chisq_ < 1.0) {
              h_matchDeltaPhi_chi2_0to1->Fill(match->deltaPhi_);
              h_matchDeltaEta_chi2_0to1->Fill(match->deltaEta_);
              h_matchDeltaR_chi2_0to1->Fill(match->deltaR_);
              h_matchDeltaZ0_chi2_0to1->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_chi2_0to1->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(fpgaTrack.chisq_ >= 1.0 && fpgaTrack.chisq_ < 2.0) {
              h_matchDeltaPhi_chi2_1to2->Fill(match->deltaPhi_);
              h_matchDeltaEta_chi2_1to2->Fill(match->deltaEta_);
              h_matchDeltaR_chi2_1to2->Fill(match->deltaR_);
              h_matchDeltaZ0_chi2_1to2->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_chi2_1to2->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(fpgaTrack.chisq_ >= 2.0 && fpgaTrack.chisq_ < 3.0) {
              h_matchDeltaPhi_chi2_2to3->Fill(match->deltaPhi_);
              h_matchDeltaEta_chi2_2to3->Fill(match->deltaEta_);
              h_matchDeltaR_chi2_2to3->Fill(match->deltaR_);
              h_matchDeltaZ0_chi2_2to3->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_chi2_2to3->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(fpgaTrack.chisq_ >= 3.0 && fpgaTrack.chisq_ < 10.0) {
              h_matchDeltaPhi_chi2_3to10->Fill(match->deltaPhi_);
              h_matchDeltaEta_chi2_3to10->Fill(match->deltaEta_);
              h_matchDeltaR_chi2_3to10->Fill(match->deltaR_);
              h_matchDeltaZ0_chi2_3to10->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_chi2_3to10->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }

            int iregion = abs(fpgaTrack.seed_/10);
            if(fpgaTrack.stubID_.size() == 4) {
              h_matchDeltaPhi_trklet_nS4[iregion]->Fill(match->deltaPhi_);
              h_matchDeltaEta_trklet_nS4[iregion]->Fill(match->deltaEta_);
              h_matchDeltaR_trklet_nS4[iregion]  ->Fill(match->deltaR_);
              h_matchDeltaZ0_trklet_nS4[iregion] ->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_trklet_nS4[iregion] ->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(fpgaTrack.stubID_.size() == 5) {
              h_matchDeltaPhi_trklet_nS5[iregion]->Fill(match->deltaPhi_);
              h_matchDeltaEta_trklet_nS5[iregion]->Fill(match->deltaEta_);
              h_matchDeltaR_trklet_nS5[iregion]  ->Fill(match->deltaR_);
              h_matchDeltaZ0_trklet_nS5[iregion] ->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_trklet_nS5[iregion] ->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(fpgaTrack.stubID_.size() >= 6) {
              h_matchDeltaPhi_trklet_nS6[iregion]->Fill(match->deltaPhi_);
              h_matchDeltaEta_trklet_nS6[iregion]->Fill(match->deltaEta_);
              h_matchDeltaR_trklet_nS6[iregion]  ->Fill(match->deltaR_);
              h_matchDeltaZ0_trklet_nS6[iregion] ->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_trklet_nS6[iregion] ->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(bestMatch == 0 ) {
              //minimal requirements for the first match
              if(match->deltaR_ < 0.2) bestMatch = match;
            } else {
              if( match->deltaR_ < bestMatch->deltaR_ )
                bestMatch = match;
            } //picked best match
          } //end minStubs requirement
        } //end loop over fpga tracks
      }//end if particle passes selection
    }//end loop over MC tracks

    //Make total sector plots
    h_nMCTrkEvt_tot->Fill(nMCTrkEvt_tot);
    for(int nSec=0; nSec<28; nSec++) {
      h_nTrkSec_eta1->Fill(min(25,nTrkSec_eta[nSec][0]));
      h_nTrkSec_eta2->Fill(min(25,nTrkSec_eta[nSec][1]));
      h_nTrkSec_eta3->Fill(min(25,nTrkSec_eta[nSec][2]));
      h_nTrkSec_eta4->Fill(min(25,nTrkSec_eta[nSec][3]));

      h_nTrkSec_tot->Fill(min(50,nTrkSec_tot[nSec]));
      h_nTrkSecTail_tot->Fill(nTrkSec_tot[nSec]);
      h_nTrkSecWODup_tot->Fill(min(50,nTrkSecWODup_tot[nSec]));
      h_nTrkSecWODupTail_tot->Fill(nTrkSecWODup_tot[nSec]);
      h_nMCTrkSec_tot->Fill(min(50,nMCTrkSec_tot[nSec]));
    }

  }//end loop over events
  
  // Write and close file
  histFile->Write();
  histFile->Close();

  //drop pointer to force reloading of input file
  fpgaEvent = 0;

  stubid.close();

}

