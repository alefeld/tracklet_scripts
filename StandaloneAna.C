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

  const bool prebin = false;
  const bool adjSector = true;
  const bool stubsTxtFile = false;
  const double Maxpt = 200.0;
  const double Minpt = 2.0;

  const int minStubs = 4;

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
  TH1F* h_mcTrkPhi_tot = new TH1F("h_mcTrkPhi_tot","mc trk Phi",200,-TMath::Pi(),TMath::Pi());
  TH1F* h_mcTrkZ0_tot  = new TH1F("h_mcTrkZ0_tot", "mc trk Z0",200,-50.,50.);
  TH1F* h_mcTrkVr_tot  = new TH1F("h_mcTrkVr_tot", "mc sq(vx^2+vy^2)",2000,0.,2.);
  TH1F* h_mcTrkEta_tot = new TH1F("h_mcTrkEta_tot","mc trk Eta",100,-5.,5.);
  TH1F* h_mcType_tot   = new TH1F("h_mcType_tot",  "mc type",801,-400.5,400.5);

  TH1F* h_trkPtWODup_tot   = new TH1F("h_trkPtWODup_tot","trk pt",200,-25.,25.);
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

    
  //For Single Muons
  TH1F* h_nMCTrkEvt_tot    = new TH1F("h_nMCTrkEvt_tot","Num MC Tracks per Event",15,0.5,15.5);
  TH1F* h_nTrkEvt_tot      = new TH1F("h_nTrkEvt_tot","Num Tracks per Event",15,0.5,15.5);
  TH1F* h_nTrkEvtWODup_tot = new TH1F("h_nTrkEvtWODup_tot","Num Tracks per Event after duplicate removal",10,0.5,10.5);
  TH1F* h_nDupTrk_tot      = new TH1F("h_nDupTrk_tot","Num of duplicates removed",10,0.5,10.5);
  TH1F* h_nMCTrkSec_tot    = new TH1F("h_nMCTrkSec_tot","Num MC Tracks per Sector",15,0.5,15.5);
  TH1F* h_nTrkSec_tot      = new TH1F("h_nTrkSec_tot","Num Tracks per Sector",15,0.5,15.5);
  TH1F* h_nTrkSecTail_tot  = new TH1F("h_nTrkSecTail_tot","Num Tracks per Sector",15,15.5,30.5);
  TH1F* h_nTrkSecWODup_tot = new TH1F("h_nTrkSecWODup_tot","Num Tracks per Sector (after dup)",15,0.5,15.5);
  TH1F* h_nTrkSecWODupTail_tot = new TH1F("h_nTrkSecWODupTail_tot","Num Tracks per Sector (after dup)",15,15.5,30.5);
  TH2F* h_NDupVSNTrkWO_tot = new TH2F("h_NDupVSNTrkWO_tot","Num Dup vs Num Trk (after)",10,0.5,10.5,10,0.5,10.5);
    
  /*
  //For ttbar
  TH1F* h_nMCTrkEvt_tot    = new TH1F("h_nMCTrkEvt_tot","Num MC Tracks per Event",51,-0.5,50.5);
  TH1F* h_nTrkEvt_tot      = new TH1F("h_nTrkEvt_tot","Num Tracks per Event",51,-0.5,50.5);
  TH1F* h_nTrkEvtWODup_tot = new TH1F("h_nTrkEvtWODup_tot","Num Tracks per Event after duplicate removal",51,-0.5,50.5);
  TH1F* h_nDupTrk_tot      = new TH1F("h_nDupTrk_tot","Num of duplicates removed",51,-0.5,50.5);
  TH1F* h_nMCTrkSec_tot    = new TH1F("h_nMCTrkSec_tot","Num MC Tracks per Sector",51,-0.5,50.5);
  TH1F* h_nTrkSec_tot      = new TH1F("h_nTrkSec_tot","Num Tracks per Sector",51,-0.5,50.5);
  TH1F* h_nTrkSecTail_tot  = new TH1F("h_nTrkSecTail_tot","Num Tracks per Sector",51,49.5,100.5);
  TH1F* h_nTrkSecWODup_tot = new TH1F("h_nTrkSecWODup_tot","Num Tracks per Sector (after dup)",51,-0.5,50.5);
  TH1F* h_nTrkSecWODupTail_tot = new TH1F("h_nTrkSecWODupTail_tot","Num Tracks per Sector (after dup)",51,49.5,100.5);
  TH2F* h_NDupVSNTrkWO_tot = new TH2F("h_NDupVSNTrkWO_tot","Num Dup vs Num Trk (after)",10,0.5,10.5,10,0.5,10.5);
  */

  TH1F* h_nStub = new TH1F("h_nStub","nStub per Event",31,-0.5,30.5);
  TH1F* h_nStub_ghost = new TH1F("h_nStub_ghost","nStub per Event with ghosts",31,-0.5,30.5);
  TH1F* h_iLTkt_tot = new TH1F("h_iLTkt_tot","Track seed",45,-22.5,22.5);
  TH1F* h_iLTktWODup_tot = new TH1F("h_iLTktWODup_tot","Track seed W/O duplicates",45,-22.5,22.5);
  TH1F* h_iLTkt_ghost = new TH1F("h_iLTkt_ghost","Track seed in ghost events",45,-22.5,22.5);

  TH1F* h_matchDeltaPhi_AllM_tot   = new TH1F("h_matchDeltaPhi_AllM_tot","fpga-mc delta phi (all)",100,-0.003,0.003);
  TH1F* h_matchDeltaEta_AllM_tot   = new TH1F("h_matchDeltaEta_AllM_tot","fpga-mc delta eta (all)",100,-0.025,0.025);
  TH1F* h_matchDeltaR_AllM_tot     = new TH1F("h_matchDeltaR_AllM_tot"  ,"fpga-mc delta R (all)"  ,100,0.,0.025);
  TH1F* h_matchDeltaZ0_AllM_tot    = new TH1F("h_matchDeltaZ0_AllM_tot"  ,"fpga-mc delta Z0 (all)"  ,100,-2.5,2.5);
  TH1F* h_matchDeltaPtOPt_AllM_tot = new TH1F("h_matchDeltaPtOPt_AllM_tot" ,"fpga-mc (delta Pt)/Pt (all)" ,100,-0.1,0.1);

  TH2F* h_matchDeltaPhi_chi2_AllM   = new TH2F("h_matchDeltaPhi_chi2_AllM","fpga-mc delta phi v. chi2",40,-0.003,0.003,40,0,20);
  TH2F* h_matchDeltaEta_chi2_AllM   = new TH2F("h_matchDeltaEta_chi2_AllM","fpga-mc delta eta v. chi2",40,-0.025,0.025,40,0,20);
  TH2F* h_matchDeltaR_chi2_AllM     = new TH2F("h_matchDeltaR_chi2_AllM"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025,40,0,20);
  TH2F* h_matchDeltaZ0_chi2_AllM    = new TH2F("h_matchDeltaZ0_chi2_AllM"  ,"fpga-mc delta Z0 v. chi2"  ,40,-2.5,2.5,40,0,20);
  TH2F* h_matchDeltaPtOPt_chi2_AllM = new TH2F("h_matchDeltaPtOPt_chi2_AllM" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1,40,0,20);

  TH1F* h_matchDeltaPhi_chi2_1to2   = new TH1F("h_matchDeltaPhi_chi2_1to2","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_1to2   = new TH1F("h_matchDeltaEta_chi2_1to2","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_1to2     = new TH1F("h_matchDeltaR_chi2_1to2"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_1to2    = new TH1F("h_matchDeltaZ0_chi2_1to2"  ,"fpga-mc delta Z0 v. chi2"  ,40,-2.5,2.5);
  TH1F* h_matchDeltaPtOPt_chi2_1to2 = new TH1F("h_matchDeltaPtOPt_chi2_1to2" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_matchDeltaPhi_chi2_2to3   = new TH1F("h_matchDeltaPhi_chi2_2to3","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_2to3   = new TH1F("h_matchDeltaEta_chi2_2to3","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_2to3     = new TH1F("h_matchDeltaR_chi2_2to3"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_2to3    = new TH1F("h_matchDeltaZ0_chi2_2to3"  ,"fpga-mc delta Z0 v. chi2"  ,40,-2.5,2.5);
  TH1F* h_matchDeltaPtOPt_chi2_2to3 = new TH1F("h_matchDeltaPtOPt_chi2_2to3" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_matchDeltaPhi_chi2_3to4   = new TH1F("h_matchDeltaPhi_chi2_3to4","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_3to4   = new TH1F("h_matchDeltaEta_chi2_3to4","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_3to4     = new TH1F("h_matchDeltaR_chi2_3to4"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_3to4    = new TH1F("h_matchDeltaZ0_chi2_3to4"  ,"fpga-mc delta Z0 v. chi2"  ,40,-2.5,2.5);
  TH1F* h_matchDeltaPtOPt_chi2_3to4 = new TH1F("h_matchDeltaPtOPt_chi2_3to4" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_matchDeltaPhi_chi2_4to10   = new TH1F("h_matchDeltaPhi_chi2_4to10","fpga-mc delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_matchDeltaEta_chi2_4to10   = new TH1F("h_matchDeltaEta_chi2_4to10","fpga-mc delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_matchDeltaR_chi2_4to10     = new TH1F("h_matchDeltaR_chi2_4to10"  ,"fpga-mc delta R v. chi2"  ,40,0.,0.025);
  TH1F* h_matchDeltaZ0_chi2_4to10    = new TH1F("h_matchDeltaZ0_chi2_4to10"  ,"fpga-mc delta Z0 v. chi2"  ,40,-2.5,2.5);
  TH1F* h_matchDeltaPtOPt_chi2_4to10 = new TH1F("h_matchDeltaPtOPt_chi2_4to10" ,"fpga-mc (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

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
    h_matchDeltaZ0_trklet_nS4[r]    = new TH1F("h_matchDeltaZ0_nS4"  ,"fpga-mc delta Z0 (nstub=4)"  ,100,-2.5,2.5);
    h_matchDeltaPtOPt_trklet_nS4[r] = new TH1F("h_matchDeltaPtOPt_nS4" ,"fpga-mc (delta Pt)/Pt (nstub=4)" ,100,-0.1,0.1);

    h_matchDeltaPhi_trklet_nS5[r]   = new TH1F("h_matchDeltaPhi_nS5","fpga-mc delta phi (nstub=5)",100,-0.003,0.003);
    h_matchDeltaEta_trklet_nS5[r]   = new TH1F("h_matchDeltaEta_nS5","fpga-mc delta eta (nstub=5)",100,-0.025,0.025);
    h_matchDeltaR_trklet_nS5[r]     = new TH1F("h_matchDeltaR_nS5"  ,"fpga-mc delta R (nstub=5)"  ,100,0.,0.025);
    h_matchDeltaZ0_trklet_nS5[r]    = new TH1F("h_matchDeltaZ0_nS5"  ,"fpga-mc delta Z0 (nstub=5)"  ,100,-2.5,2.5);
    h_matchDeltaPtOPt_trklet_nS5[r] = new TH1F("h_matchDeltaPtOPt_nS5" ,"fpga-mc (delta Pt)/Pt (nstub=5)" ,100,-0.1,0.1);

    h_matchDeltaPhi_trklet_nS6[r]   = new TH1F("h_matchDeltaPhi_nS6","fpga-mc delta phi (nstub=6)",100,-0.003,0.003);
    h_matchDeltaEta_trklet_nS6[r]   = new TH1F("h_matchDeltaEta_nS6","fpga-mc delta eta (nstub=6)",100,-0.025,0.025);
    h_matchDeltaR_trklet_nS6[r]     = new TH1F("h_matchDeltaR_nS6"  ,"fpga-mc delta R (nstub=6)"  ,100,0.,0.025);
    h_matchDeltaZ0_trklet_nS6[r]    = new TH1F("h_matchDeltaZ0_nS6"  ,"fpga-mc delta Z0 (nstub=6)"  ,100,-2.5,2.5);
    h_matchDeltaPtOPt_trklet_nS6[r] = new TH1F("h_matchDeltaPtOPt_nS6" ,"fpga-mc (delta Pt)/Pt (nstub=6)" ,100,-0.1,0.1);
  }

  // Loop over events
  for(int i=0; i<numEvts; i++) {
  //for(int i=0; i<21; i++) {

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
    
    // Duplicate Removal Section
    // Map of track ID (number) to duplicates
    std::map<int,bool> duplicateMap;
    for(int j=0; j<(int)fpgaEvent->tracks.size(); j++) duplicateMap.insert(std::pair<int,bool>(j,false) );
    
    int nShare[1000][1000]={{0}};
    int nStub[1000]={0};

    // Loop over tracks in this event to count stubs
    for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++) {

      // Get this track, count stubs
      FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
      nStub[ntrk] = track.stubID_.size();
      /*
      cout << track.nLayerMatch_ << "+" << track.nDiskMatch_;
      cout << ", " << track.iLTkt_ << ", " <<track.stubID_.size() << endl;
      for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
        cout << st->first << "," << st->second << " ";
      }
      cout << endl;
      */
    }
    /*
    cout << "Event: " << i << endl;
      */

    // Loop over tracks in this event
    for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++) {

      // Get this track
      FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);

      // Pass nStub cut, pt cut, duplicate cut
      if(nStub[ntrk]>=minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt && !duplicateMap.at(ntrk)) {

        //Loop over other tracks
        for(int ntrk2=ntrk+1; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
          FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);

          // Pass nStub cut, pt cut, duplicate cut
          if(nStub[ntrk2]>=minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && !duplicateMap.at(ntrk2) ) {

            //Determine how many stubs these two tracks share
            for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
              if(track2.stubID_.find(st->first) != track2.stubID_.end()) {
                if(st->second == track2.stubID_[st->first] && st->second != 63) {
                  nShare[ntrk][ntrk2]++;
                  nShare[ntrk2][ntrk]++;
                }
              }
            }
            //cout << "nStub[" << ntrk2 << "] = " << nStub[ntrk2] << ", nShare[" << ntrk2 << "] = " << nShare[ntrk2] << endl;
          }
        } //end loop of secondary tracks

        // We've counted stubs and nshare, so remove duplicates in the same sector
        for(int ntrk2=ntrk+1; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
          FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);
          if(nStub[ntrk2]>=minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && (track2.nSector_ == track.nSector_) && !duplicateMap.at(ntrk2) ) {
            //Decide whether we want to flag track as duplicate

            //chi2 method
            if((nStub[ntrk]-nShare[ntrk][ntrk2])<3 || (nStub[ntrk2]-nShare[ntrk][ntrk2])<3) {
              if(track.chiS_ > track2.chiS_) duplicateMap.at(ntrk)=true;
              if(track.chiS_ <= track2.chiS_) duplicateMap.at(ntrk2)=true;
            }

            //nstub method
            //Remove the track with the least number of stubs, then by seeding layer
//            if((nStub[ntrk]-nShare[ntrk][ntrk2])<3 && nStub[ntrk2]>nStub[ntrk]){
//              duplicateMap.at(ntrk)=true;
//            }
//            if((nStub[ntrk2]-nShare[ntrk][ntrk2])<3 && nStub[ntrk]>=nStub[ntrk2]){
//              duplicateMap.at(ntrk2)=true;
//            }

          }//end track2 selection
        }//end loop over secondary tracks
      }//end track1 selection
    }//end loop over primary tracks
    
    // Adjacent sector removal
    if(adjSector) {
      for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++) {
        // Get this track, pass nStub cut, pt cut, duplicate cut
        FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
        if(nStub[ntrk]>=minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt && !duplicateMap.at(ntrk)) {
          //Loop over other tracks to look for duplicates
          for(int ntrk2=0; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
            // Pass nStub cut, pt cut, duplicate cut
            // Select adjacent-sector tracks for removal (N+1 sector only)
            FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);
            if(nStub[ntrk2]>=minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && !duplicateMap.at(ntrk2)) {
              if((track2.nSector_ == track.nSector_+1) || ((track.nSector_==27) && (track2.nSector_==0))) {
//              if(abs(track2.nSector_-track.nSector_)==1 || abs(track2.nSector_-track.nSector_)==27) {
                // Only try to flag secondary track as duplicate
                if((nStub[ntrk2]-nShare[ntrk][ntrk2])<3 || (nStub[ntrk]-nShare[ntrk][ntrk2])<3) duplicateMap.at(ntrk2)=true;
//                if((nStub[ntrk2]-nShare[ntrk][ntrk2])<3 && nStub[ntrk] >= nStub[ntrk2]) duplicateMap.at(ntrk2)=true;
//                if((nStub[ntrk]-nShare[ntrk][ntrk2])<3 && nStub[ntrk] < nStub[ntrk2]) duplicateMap.at(ntrk)=true;
              }// end track2 sector selection
            }// end track2 selection
          }//end loop over secondary tracks
        }//end track1 selection
      }//end loop over primary tracks
    }//end adjacent sector removal
    
    //Loop through tracks again to find out how many were flagged as duplicates
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
      if(nStub[ntrk] >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt) {

        nILTkt[track.iLTkt_]++;

        nTrkSec_tot[track.nSector_]++;
        if( track.eta_<=-1.0 ) nTrkSec_eta[track.nSector_][0]++;
        if( track.eta_>-1.0 && track.eta_<=0) nTrkSec_eta[track.nSector_][1]++;
        if( track.eta_<1.0 && track.eta_>0) nTrkSec_eta[track.nSector_][2]++;
        if( track.eta_>=1.0 ) nTrkSec_eta[track.nSector_][3]++;

        h_iLTkt_tot->Fill(track.iLTkt_);

        if(duplicateMap.at(ntrk)) {
          nDupTrk_tot++;
        } else {
          //if track isn't a duplicate
          nTrkWODup_tot++;
          nTrkSecWODup_tot[track.nSector_]++;

          h_iLTktWODup_tot->Fill(track.iLTkt_);

          h_trkPtWODup_tot->Fill(abs(track.pt_));
          h_trkPhiWODup_tot->Fill(track.phi_);
          h_trkEtaWODup_tot->Fill(track.eta_);
          h_trkZ0WODup_tot->Fill(track.z0_);

        }//end duplicates loop
      }//end if nStub >= nStubs
    }//end loop over tracks

    //make event level plots, but suppress events with no tracks passing requirements
    if(nTrkWODup_tot>0) {

      //total plots
      h_nTrkEvt_tot->Fill(TMath::Min(nTrkWODup_tot+nDupTrk_tot,99));
      h_nTrkEvtWODup_tot->Fill(TMath::Min(nTrkWODup_tot,99));
      h_nDupTrk_tot->Fill(TMath::Min(nDupTrk_tot,99));
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
          if(nStub[ntrk] >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt && !duplicateMap.at(ntrk)) {

            h_iLTkt_ghost->Fill(track.iLTkt_);

            //Plot simple variables
            h_trkPhiMulti_tot->Fill(track.phi_);
            h_trkEtaMulti_tot->Fill(track.eta_);
            h_trkZ0Multi_tot->Fill(track.z0_);
            h_itrkPhiMulti_tot->Fill(track.iphi_);
            h_trkPtMulti_tot->Fill(abs(track.pt_));

            //Make a superposition of the sectors
            if(track.phi_<0) {
              h_trkPhiMultiLayered_tot->Fill(track.phi_+TMath::TwoPi()/28);
            } else {
              int nPhiSec=(track.phi_/(TMath::TwoPi()/28));
              h_trkPhiMultiLayered_tot->Fill(track.phi_-nPhiSec*TMath::TwoPi()/28);
            }
          }
        }

        // Write out track information to .txt file
        if(stubsTxtFile) {
          for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
            FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
            if(nStub[ntrk] >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt && !duplicateMap.at(ntrk)) {
              //Output stubids to .txt file
              int nShareGhost=0;
              stubid << Form("\nTrack %i: ",ntrk);
              stubid << Form(" Phi=%f, Eta=%f, Z0=%f, Pt=%f, iLTkt=%i, Sector=%i", track.phi_,track.eta_,track.z0_,abs(track.pt_),track.iLTkt_,track.nSector_);
              stubid << Form("\n    Stub ids: ");
              for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
                stubid << Form("  %i,%i",st->first,st->second);
                for(int ntrk2=0; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
                  FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);
                  if(!duplicateMap.at(ntrk2) && nStub[ntrk2] >= minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && ntrk!=ntrk2) {
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
    if(nTrkWODup_tot>0) {
      for(int ntrk=0; ntrk<(int)fpgaEvent->mcTracks.size(); ntrk++){

        FPGAEventMCTrack track = fpgaEvent->mcTracks.at(ntrk);

        //apply some track requirement
        if(abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt
           && sqrt(track.vy_*track.vy_ + track.vx_*track.vx_)<0.1  
           && checkType(track.type_,true)
           && fabs(track.eta_)<2.5
           && fabs(track.vz_)<20
           ) {

          //if(!checkType(track.type_,false)) printf("Type neither good nor bad %i\n",track.type_);
          h_mcTrkPt_tot->Fill(abs(track.pt_));  
          h_mcTrkPhi_tot->Fill(track.phi_); 
          h_mcTrkZ0_tot->Fill(track.vz_);
          h_mcTrkVr_tot->Fill(sqrt(track.vy_*track.vy_ + track.vx_*track.vx_));   
          h_mcTrkEta_tot->Fill(track.eta_);
          h_mcType_tot->Fill(track.type_);

          //Tally number of tracks per event and sector (artifically)
          nMCTrkEvt_tot++;
          nMCTrkSec_tot[(int)((track.phi_+TMath::Pi())/(TMath::TwoPi()/28))]++;

          // Try to see if this particle can be matched to a found track
          FPGAEventTrackMatch *bestMatch = 0;
          for(int fntrk=0; fntrk<(int)fpgaEvent->tracks.size(); fntrk++){
            if(!duplicateMap.at(fntrk)) {
              FPGAEventTrack fpgaTrack = fpgaEvent->tracks.at(fntrk);

              int nStubMC = fpgaTrack.stubID_.size();

              //discriminate between tracklet regions
              int r = -1;
              if(fpgaTrack.iLTkt_ == 1 || fpgaTrack.iLTkt_ == 3 || fpgaTrack.iLTkt_ == 5)
                r = 0; // barrel only tracklet region
              if(abs(fpgaTrack.iLTkt_) == 11 || abs(fpgaTrack.iLTkt_) == 13)
                r = 1; // disk only tracklet regions
              if(abs(fpgaTrack.iLTkt_) == 21 || abs(fpgaTrack.iLTkt_) == 22) {
                r = 2; // overlap tracklet region
//                cout << "Stubs: ";
//                for(std::map<int, int>::iterator  st=fpgaTrack.stubID_.begin(); st!=fpgaTrack.stubID_.end(); st++) {
//                  cout << st->first << " ";
//                }
//                cout << endl;
              }

              if(nStubMC >= minStubs && !duplicateMap.at(fntrk)) {
                FPGAEventTrackMatch *match = new FPGAEventTrackMatch(track,fpgaTrack);

                h_matchDeltaPhi_AllM_tot->Fill(match->deltaPhi_);
                h_matchDeltaEta_AllM_tot->Fill(match->deltaEta_);
                h_matchDeltaR_AllM_tot  ->Fill(match->deltaR_);
                h_matchDeltaZ0_AllM_tot ->Fill(match->deltaZ0_);
                h_matchDeltaPtOPt_AllM_tot ->Fill(match->deltaPt_/abs(track.pt_));

                h_matchDeltaPhi_chi2_AllM->Fill(match->deltaPhi_,fpgaTrack.chiS_);
                h_matchDeltaEta_chi2_AllM->Fill(match->deltaEta_,fpgaTrack.chiS_);
                h_matchDeltaR_chi2_AllM->Fill(match->deltaR_,fpgaTrack.chiS_);
                h_matchDeltaZ0_chi2_AllM->Fill(match->deltaZ0_,fpgaTrack.chiS_);
                h_matchDeltaPtOPt_chi2_AllM->Fill(match->deltaPt_/abs(track.pt_),fpgaTrack.chiS_);

                if(fpgaTrack.chiS_ >= 1.0 && fpgaTrack.chiS_ < 2.0) {
                  h_matchDeltaPhi_chi2_1to2->Fill(match->deltaPhi_);
                  h_matchDeltaEta_chi2_1to2->Fill(match->deltaEta_);
                  h_matchDeltaR_chi2_1to2->Fill(match->deltaR_);
                  h_matchDeltaZ0_chi2_1to2->Fill(match->deltaZ0_);
                  h_matchDeltaPtOPt_chi2_1to2->Fill(match->deltaPt_/abs(track.pt_));
                }

                if(fpgaTrack.chiS_ >= 2.0 && fpgaTrack.chiS_ < 3.0) {
                  h_matchDeltaPhi_chi2_2to3->Fill(match->deltaPhi_);
                  h_matchDeltaEta_chi2_2to3->Fill(match->deltaEta_);
                  h_matchDeltaR_chi2_2to3->Fill(match->deltaR_);
                  h_matchDeltaZ0_chi2_2to3->Fill(match->deltaZ0_);
                  h_matchDeltaPtOPt_chi2_2to3->Fill(match->deltaPt_/abs(track.pt_));
                }

                if(fpgaTrack.chiS_ >= 3.0 && fpgaTrack.chiS_ < 4.0) {
                  h_matchDeltaPhi_chi2_3to4->Fill(match->deltaPhi_);
                  h_matchDeltaEta_chi2_3to4->Fill(match->deltaEta_);
                  h_matchDeltaR_chi2_3to4->Fill(match->deltaR_);
                  h_matchDeltaZ0_chi2_3to4->Fill(match->deltaZ0_);
                  h_matchDeltaPtOPt_chi2_3to4->Fill(match->deltaPt_/abs(track.pt_));
                }

                if(fpgaTrack.chiS_ >= 4.0 && fpgaTrack.chiS_ < 10.0) {
                  h_matchDeltaPhi_chi2_4to10->Fill(match->deltaPhi_);
                  h_matchDeltaEta_chi2_4to10->Fill(match->deltaEta_);
                  h_matchDeltaR_chi2_4to10->Fill(match->deltaR_);
                  h_matchDeltaZ0_chi2_4to10->Fill(match->deltaZ0_);
                  h_matchDeltaPtOPt_chi2_4to10->Fill(match->deltaPt_/abs(track.pt_));
                }

                if(nStubMC == 4) {
                  h_matchDeltaPhi_trklet_nS4[r]->Fill(match->deltaPhi_);
                  h_matchDeltaEta_trklet_nS4[r]->Fill(match->deltaEta_);
                  h_matchDeltaR_trklet_nS4[r]  ->Fill(match->deltaR_);
                  h_matchDeltaZ0_trklet_nS4[r] ->Fill(match->deltaZ0_);
                  h_matchDeltaPtOPt_trklet_nS4[r] ->Fill(match->deltaPt_/abs(track.pt_));
                }
                if(nStubMC == 5) {
                  h_matchDeltaPhi_trklet_nS5[r]->Fill(match->deltaPhi_);
                  h_matchDeltaEta_trklet_nS5[r]->Fill(match->deltaEta_);
                  h_matchDeltaR_trklet_nS5[r]  ->Fill(match->deltaR_);
                  h_matchDeltaZ0_trklet_nS5[r] ->Fill(match->deltaZ0_);
                  h_matchDeltaPtOPt_trklet_nS5[r] ->Fill(match->deltaPt_/abs(track.pt_));
                }
                if(nStubMC >= 6) {
                  h_matchDeltaPhi_trklet_nS6[r]->Fill(match->deltaPhi_);
                  h_matchDeltaEta_trklet_nS6[r]->Fill(match->deltaEta_);
                  h_matchDeltaR_trklet_nS6[r]  ->Fill(match->deltaR_);
                  h_matchDeltaZ0_trklet_nS6[r] ->Fill(match->deltaZ0_);
                  h_matchDeltaPtOPt_trklet_nS6[r] ->Fill(match->deltaPt_/abs(track.pt_));
                }

                if(bestMatch == 0 ) {
                  //minimal requirements for the first match
                  if(match->deltaR_ < 0.2 && !duplicateMap.at(fntrk) ) bestMatch = match;
                } else {
                  if( match->deltaR_ < bestMatch->deltaR_ && !duplicateMap.at(fntrk) )
                    bestMatch = match;
                } //picked best match
              } //end minStubs requirement
            } //end not a duplicate requirement
          } //end loop over fpga tracks
        }//end if particle passes selection
      }//end loop over MC tracks
    }//end if event has FPGA track

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

