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
FPGAEvent *fpgaEvent=0;
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

void analyzeEvent(TString rootFileName="myTest.root",bool Hourglass = true) {

  TString histFileName = "ana_"+rootFileName;
  cout << "Output file will be named ana_" << rootFileName << endl;

  const bool stubsTxtFile = false;
  const double Maxpt = 200.0;
  const double Minpt = 2.0;

  const unsigned int NSector=Hourglass?9:27;  // 9 hourglass, 27 standard
  //const unsigned int NSector = 27;
  //if(Hourglass) NSector = 9;
  const int nEtaBins = 1;
  const double etaBinSize = 2.4/nEtaBins*2;
  const int nPhiBins = 1;
  const double phiBinSize = 2*M_PI/(NSector*nPhiBins);

  const int minStubs = 4;
  const int minIndStubs = 3;

  ofstream stubid ("StubIDMulti.txt");
  ofstream asymstubidplus ("AsymStubsPlus.txt");
  ofstream asymstubidminus ("AsymStubsMinus.txt");

//  loadFPGATree(rootFileName);
  if(fpgaEvent==0) loadFPGATree(rootFileName);

  //Setup histogram locations
  TFile *histFile = new TFile(histFileName,"RECREATE");
  histFile->cd();
  TDirectory *barrelOnly = histFile->mkdir("barrelOnly");
  TDirectory *diskOnly = histFile->mkdir("diskOnly");
  TDirectory *overlap = histFile->mkdir("overlap");
  histFile->cd();

  //Initialize histograms NOT separated by region
    
  //For Single Muons
//  TH1F* h_nMCTrkEvt_tot    = new TH1F("h_nMCTrkEvt_tot","Num MC Tracks per Event",15,0.5,15.5);
//  TH1F* h_nTrkEvt_tot      = new TH1F("h_nTrkEvt_tot","Num Tracks per Event",15,0.5,15.5);
//  TH1F* h_nTrkEvtWODup_tot = new TH1F("h_nTrkEvtWODup_tot","Num Tracks per Event after duplicate removal",10,0.5,10.5);
//  TH1F* h_nTrkEvtWODupHybrid = new TH1F("h_nTrkEvtWODupHybrid","Num Tracks per Event after hybrid-like duplicate removal",10,0.5,10.5);
//  TH1F* h_nDupTrk_tot      = new TH1F("h_nDupTrk_tot","Num of duplicates removed",10,0.5,10.5);
//  TH1F* h_nMCTrkSec_tot    = new TH1F("h_nMCTrkSec_tot","Num MC Tracks per Sector",15,0.5,15.5);
//  TH1F* h_nTrkSec_tot      = new TH1F("h_nTrkSec_tot","Num Tracks per Sector",15,0.5,15.5);
//  TH1F* h_nTrkSecTail_tot  = new TH1F("h_nTrkSecTail_tot","Num Tracks per Sector",15,15.5,30.5);
//  TH1F* h_nTrkSecWODup_tot = new TH1F("h_nTrkSecWODup_tot","Num Tracks per Sector (after dup)",15,0.5,15.5);
//  TH1F* h_nTrkSecWODupHybrid = new TH1F("h_nTrkSecWODupHybrid","Num Tracks per Sector hybrid (after PD)",15,0.5,15.5);
//  TH1F* h_nTrkSecWODupTail_tot = new TH1F("h_nTrkSecWODupTail_tot","Num Tracks per Sector (after dup)",15,15.5,30.5); 
//  TH2F* h_NDupVSNTrkWO_tot = new TH2F("h_NDupVSNTrkWO_tot","Num Dup vs Num Trk (after)",10,0.5,10.5,10,0.5,10.5);

  //For ttbar
  TH1F* h_nMCTrkEvt_tot    = new TH1F("h_nMCTrkEvt_tot","Num MC Tracks per Event",100,-0.5,300.5);
  TH1F* h_nTrkEvt_tot      = new TH1F("h_nTrkEvt_tot","Num Tracks per Event",160,-0.5,800.5);
  TH1F* h_nTrkEvtWODup_tot = new TH1F("h_nTrkEvtWODup_tot","Num Tracks per Event after duplicate removal",80,-0.5,400.5);
  TH1F* h_nTrkEvtWODupHybrid = new TH1F("h_nTrkEvtWODupHybrid","Num Tracks per Event after hybrid-like duplicate removal",80,-0.5,400.5);
  TH1F* h_nDupTrk_tot      = new TH1F("h_nDupTrk_tot","Num of duplicates removed",80,-0.5,550.5);
  TH1F* h_nMCTrkSec_tot    = new TH1F("h_nMCTrkSec_tot","Num MC Tracks per Sector",51,-0.5,50.5);
  TH1F* h_nTrkSec_tot      = new TH1F("h_nTrkSec_tot","Num Tracks per Sector",50,-0.5,149.5);
  TH1F* h_nTrkSecTail_tot  = new TH1F("h_nTrkSecTail_tot","Num Tracks per Sector",50,149.5,249.5);
  TH1F* h_nTrkSecWODup_tot = new TH1F("h_nTrkSecWODup_tot","Num Tracks per Sector (after dup)",51,-0.5,50.5);
  TH1F* h_nTrkSecWODupHybrid = new TH1F("h_nTrkSecWODupHybrid","Num Tracks per Sector hybrid (after PD)",51,-0.5,50.5);
  TH1F* h_nTrkSecWODupTail_tot = new TH1F("h_nTrkSecWODupTail_tot","Num Tracks per Sector (after PD)",51,49.5,100.5);
  TH2F* h_NDupVSNTrkWO_tot = new TH2F("h_NDupVSNTrkWO_tot","Num Dup vs Num Trk (after)",10,0.5,10.5,10,0.5,10.5);

  TH1F* h_nStubTrkHybridResim = new TH1F("h_nStubTrkHybridResim","Num stubs per track (hybrid)",31,-0.5,30.5);
  TH1F* h_nStubTrkHybridResim_after = new TH1F("h_nStubTrkHybridResim_after","Num stubs per track (hybrid) after PD",31,-0.5,30.5);

  TH1F* h_StubID_presel = new TH1F("h_StubID_presel","StubID distribution eta",150,-100,250);
  TH1F* h_StubID_tot = new TH1F("h_StubID_tot","StubID distribution eta",150,-100,250);
  TH1F* h_StubID_plus = new TH1F("h_StubID_plus","StubID distribution eta>1.4",150,-100,250);
  TH1F* h_StubID_minus = new TH1F("h_StubID_minus","StubID distribution eta<-1.4",150,-100,250);

  TH1F* h_nStub = new TH1F("h_nStub","nStub per Event",31,-0.5,30.5);
  TH1F* h_nStubLayer = new TH1F("h_nStubLayer","nStub per Layer",49,-24.5,24.5);
  TH1F* h_nStubTrk_premerge = new TH1F("h_nStubTrk_premerge","Num Stubs per Track before merging",21,-0.5,20.5);
  TH1F* h_nStubTrk_prefit = new TH1F("h_nStubTrk_prefit","Num Stubs per Track before fitting",21,-0.5,20.5);
  TH1F* h_nStubTrk_postfit = new TH1F("h_nStubTrk_postfit","Num Stubs per Track after fitting",21,-0.5,20.5);
  TH1F* h_nStubTrk = new TH1F("h_nStubTrk","Num Stubs per Track",21,-0.5,20.5);
  TH1F* h_nStubTrkWODup = new TH1F("h_nStubTrkWODup","Num Stubs per Track (no duplicates)",21,-0.5,20.5);
  TH1F* h_nStub_ghost = new TH1F("h_nStub_ghost","nStub per Event with ghosts",31,-0.5,30.5);
  TH1F* h_iLTkt_tot = new TH1F("h_iLTkt_tot","Track seed",45,-22.5,22.5);
  TH1F* h_iLTktWODup_tot = new TH1F("h_iLTktWODup_tot","Track seed W/O duplicates",45,-22.5,22.5);
  TH1F* h_iLTkt_ghost = new TH1F("h_iLTkt_ghost","Track seed in ghost events",45,-22.5,22.5);


  TH1F* h_nTrkSec2_tot;
  TH1F* h_secLatency;
  TH1F* h_secComparisons;
  TH1F* h_secComparisons_evt1;

  // Standard Sector
  if(!Hourglass) {
    h_nTrkSec2_tot = new TH1F("h_nTrkSec2_tot","Num Tracks per Sector",51,-0.5,50.5);
    h_secLatency = new TH1F("h_secLatency","Latency per sector",50,0.5,49.5);
    h_secComparisons = new TH1F("h_secComparisons","Number of comparisons per event (per sector)",100,0.5,499.5);
    h_secComparisons_evt1 = new TH1F("h_secComparisons_evt1","Number of comparisons per clock (evt 1)",50,0.5,49.5);
  }

  // Hourglass Sector
  if(Hourglass) {
    h_nTrkSec2_tot = new TH1F("h_nTrkSec2_tot","Num Tracks per Sector",51,-0.5,150.5);
    h_secLatency = new TH1F("h_secLatency","Latency per sector",75,0.5,74.5);
    h_secComparisons = new TH1F("h_secComparisons","Number of comparisons per event (per sector)",100,0.5,4000.5);
    h_secComparisons_evt1 = new TH1F("h_secComparisons_evt1","Number of comparisons per clock (on sector)",50,0.5,149.5);
  }

  TH1F* h_nTrkBin;
  TH1F* h_binLatency;
  TH1F* h_binComparisons;
  TH1F* h_binComparisons_evt1;

  // Binned Standard
  if(!Hourglass) {
    h_nTrkBin = new TH1F("h_nTrkBin","Num Tracks per Bin",51,-0.5,50.5);
    h_binLatency = new TH1F("h_binLatency","Latency per bin",50,0.5,49.5);
    h_binComparisons = new TH1F("h_binComparisons","Number of comparisons per bin",100,0.5,2000.5);
    h_binComparisons_evt1 = new TH1F("h_binComparisons_evt1","Number of comparisons per clock (one bin)",50,0.5,149.5);
  }

  // Binned Hourglass
  if(Hourglass) {
    h_nTrkBin = new TH1F("h_nTrkBin","Num Tracks per Bin",51,-0.5,150.5);
    h_binLatency = new TH1F("h_binLatency","Latency per bin",75,0.5,74.5);
    h_binComparisons = new TH1F("h_binComparisons","Number of comparisons per bin",100,0.5,3000.5);
    h_binComparisons_evt1 = new TH1F("h_binComparisons_evt1","Number of comparisons per clock (one bin)",50,0.5,149.5);
  }


  TH1F* h_mcTrkPt_tot  = new TH1F("h_mcTrkPt_tot", "mc trk Pt",200,0.,200.);
  TH1F* h_mcTrkPt_low_tot  = new TH1F("h_mcTrkPt_low_tot", "mc trk Pt",50,0.,10.);
  TH1F* h_mcTrkPhi_tot = new TH1F("h_mcTrkPhi_tot","mc trk Phi",200,-TMath::Pi(),TMath::Pi());
  TH1F* h_mcTrkZ0_tot  = new TH1F("h_mcTrkZ0_tot", "mc trk Z0",200,-50.,50.);
  TH1F* h_mcTrkVr_tot  = new TH1F("h_mcTrkVr_tot", "mc sq(vx^2+vy^2)",2000,0.,2.);
  TH1F* h_mcTrkEta_tot = new TH1F("h_mcTrkEta_tot","mc trk Eta",100,-5.,5.);
  TH1F* h_mcType_tot   = new TH1F("h_mcType_tot",  "mc type",801,-400.5,400.5);

  TH1F* h_trkPt_tot   = new TH1F("h_trkPt_tot","trk pt",100,0.,25.);
  TH1F* h_trkPt_low_tot   = new TH1F("h_trkPt_low_tot","trk pt",50,0.,10.);
  TH1F* h_trkPhi_tot  = new TH1F("h_trkPhi_tot","trk (global) phi",100,0.,TMath::TwoPi());
  TH1F* h_trkEta_tot  = new TH1F("h_trkEta_tot","trk eta",100,-4.,4.);
  TH1F* h_trkZ0_tot   = new TH1F("h_trkZ0_tot","trk z0",200,-50.,50.);

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

  TH2F* h_chisq_ichisq    = new TH2F("h_chisq_ichisq","ichisq of tracks",100,0,50,256,0,1000);
  //TH2F* h_chisq_ichisq    = new TH2F("h_chisq_ichisq","ichisq of tracks",256,0,256,256,0,256);
  TH2F* h_chisqDiff16     = new TH2F("h_chisqDiff16","ichisq/16 vs. diff",100,-50,50,100,0,50);
  TH2F* h_chisqDiff16Zoom = new TH2F("h_chisqDiff16Zoom","ichisq/16 vs. diff (zoom)",50,-5,5,30,0,5);
  TH1F* h_chisq_beforePD    = new TH1F("h_chisq_beforePD","chisq of tracks",100,0,128);
  TH1F* h_chisq_afterPD     = new TH1F("h_chisq_afterPD","chisq of tracks",100,0,128);
  TH1F* h_ichisq_beforePD    = new TH1F("h_ichisq_beforePD","ichisq of tracks",256,0,1023);
  TH1F* h_ichisq_afterPD     = new TH1F("h_ichisq_afterPD","ichisq of tracks",256,0,1023);
  //TH1F* h_chisqDiff16_proj[10];
  /*for(int i=0; i<10; i++) {
    h_chisqDiff16_proj[i] = new TH1F("h_chisqDiff16_proj"+str(i),"Projection",10,-5,5);
  }*/

  TH2F* h_gridMapPlus = new TH2F("h_gridMapPlus","Track location (+ eta)",40,-0.5,39,35,-0.5,34.5);
  TH2F* h_gridMapPlusWODup = new TH2F("h_gridMapPlusWODup","Track location (+ eta) noDup",40,-0.5,39,35,-0.5,34.5);
  TH2F* h_gridMapMinus = new TH2F("h_gridMapMinus","Track location (- eta)",40,-0.5,39,35,-0.5,34.5);
  TH2F* h_gridMapMinusWODup = new TH2F("h_gridMapMinusWODup","Track location (- eta) noDup",40,-0.5,39,35,-0.5,34.5);

  TH2F* h_gridMap = new TH2F("h_gridMap","Track location in parameter space",40,-0.5,39,35,-0.5,34.5);
  TH2F* h_gridMapWODup = new TH2F("h_gridMapWODup","Track location in parameter space (no duplicates)",40,-0.5,39,35,-0.5,34.5);

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

  TH2F* h_IDvsiZ = new TH2F("h_IDvsiZ", "Stub ID vs iZ",200,-200,200,64,-0.5,512);
  TH2F* h_IDvsZ = new TH2F("h_IDvsZ", "Stub ID vs Z",110,-110,200,64,-0.5,512);
  TH2F* h_IDvsLayer = new TH2F("h_IDvsLayer", "Stub ID vs Layer",41,-20.5,20.5,64,-0.5,512);
  
  // For Resolution vs. nStub by region

  TH1F* h_matchDeltaPhi_4stub[3];
  TH1F* h_matchDeltaEta_4stub[3];
  TH1F* h_matchDeltaR_4stub[3];
  TH1F* h_matchDeltaZ0_4stub[3];
  TH1F* h_matchDeltaPtOPt_4stub[3];

  TH1F* h_matchDeltaPhi_5stub[3];
  TH1F* h_matchDeltaEta_5stub[3];
  TH1F* h_matchDeltaR_5stub[3];
  TH1F* h_matchDeltaZ0_5stub[3];
  TH1F* h_matchDeltaPtOPt_5stub[3];

  TH1F* h_matchDeltaPhi_6stub[3];
  TH1F* h_matchDeltaEta_6stub[3];
  TH1F* h_matchDeltaR_6stub[3];
  TH1F* h_matchDeltaZ0_6stub[3];
  TH1F* h_matchDeltaPtOPt_6stub[3];

  for(int r=0; r<3; r++) {

    if(r==0) barrelOnly->cd();
    if(r==1) diskOnly->cd();
    if(r==2) overlap->cd();

    h_matchDeltaPhi_4stub[r]   = new TH1F("h_matchDeltaPhi_4stub","fpga-mc delta phi (nstub=4)",100,-0.003,0.003);
    h_matchDeltaEta_4stub[r]   = new TH1F("h_matchDeltaEta_4stub","fpga-mc delta eta (nstub=4)",100,-0.025,0.025);
    h_matchDeltaR_4stub[r]     = new TH1F("h_matchDeltaR_4stub","fpga-mc delta R (nstub=4)"  ,100,0.,0.025);
    h_matchDeltaZ0_4stub[r]    = new TH1F("h_matchDeltaZ0_4stub","fpga-mc delta Z0 (nstub=4)"  ,100,-1.5,1.5);
    h_matchDeltaPtOPt_4stub[r] = new TH1F("h_matchDeltaPtOPt_4stub","fpga-mc (delta Pt)/Pt (nstub=4)" ,100,-0.1,0.1);

    h_matchDeltaPhi_5stub[r]   = new TH1F("h_matchDeltaPhi_5stub","fpga-mc delta phi (nstub=5)",100,-0.003,0.003);
    h_matchDeltaEta_5stub[r]   = new TH1F("h_matchDeltaEta_5stub","fpga-mc delta eta (nstub=5)",100,-0.025,0.025);
    h_matchDeltaR_5stub[r]     = new TH1F("h_matchDeltaR_5stub"  ,"fpga-mc delta R (nstub=5)"  ,100,0.,0.025);
    h_matchDeltaZ0_5stub[r]    = new TH1F("h_matchDeltaZ0_5stub"  ,"fpga-mc delta Z0 (nstub=5)"  ,100,-1.5,1.5);
    h_matchDeltaPtOPt_5stub[r] = new TH1F("h_matchDeltaPtOPt_5stub" ,"fpga-mc (delta Pt)/Pt (nstub=5)" ,100,-0.1,0.1);

    h_matchDeltaPhi_6stub[r]   = new TH1F("h_matchDeltaPhi_6stub","fpga-mc delta phi (nstub=6)",100,-0.003,0.003);
    h_matchDeltaEta_6stub[r]   = new TH1F("h_matchDeltaEta_6stub","fpga-mc delta eta (nstub=6)",100,-0.025,0.025);
    h_matchDeltaR_6stub[r]     = new TH1F("h_matchDeltaR_6stub"  ,"fpga-mc delta R (nstub=6)"  ,100,0.,0.025);
    h_matchDeltaZ0_6stub[r]    = new TH1F("h_matchDeltaZ0_6stub"  ,"fpga-mc delta Z0 (nstub=6)"  ,100,-1.5,1.5);
    h_matchDeltaPtOPt_6stub[r] = new TH1F("h_matchDeltaPtOPt_6stub" ,"fpga-mc (delta Pt)/Pt (nstub=6)" ,100,-0.1,0.1);
  }

  histFile->cd();

  int nGridEvts=0;
  TH2F* h_gridMap_manyTracks[10];
  TH2F* h_gridMapBig_manyTracks[10];
  TH2F* h_gridMapGlobal_manyTracks[10];
  TH2F* h_gridMapGlobalWODup_manyTracks[10];
  for(int g=0; g<10; g++) {
    TString name = Form("h_gridMap_manyTracks_%i",g);
    TString nameBig = Form("h_gridMapBig_manyTracks_%i",g);
    TString nameGlobal = Form("h_gridMapGlobal_manyTracks_%i",g);
    TString nameGlobalWODup = Form("h_gridMapGlobalWODup_manyTracks_%i",g);
    h_gridMap_manyTracks[g] = new TH2F(name, "Grid map (many tracks)",40,-0.5,39,19,-0.5,18.5);
    h_gridMapBig_manyTracks[g] = new TH2F(nameBig, "Grid map (many tracks)",40,-0.5,39,35,-0.5,34.5);
    h_gridMapGlobal_manyTracks[g] = new TH2F(nameGlobal, "Grid map (many tracks)",450,-0.5,449.5,40,-0.5,39);
    h_gridMapGlobalWODup_manyTracks[g] = new TH2F(nameGlobalWODup, "Grid map (many tracks)",450,-0.5,449.5,40,-0.5,39);
  }

  // Loop over events
  for(int i=0; i<numEvts; i++) {
//  for(int i=0; i<1; i++) {
    if(!(i%100)) cout << "Analyzing events  " << i << "-" << i+99 << "..." << endl;
    fpgaTree->GetEntry(i);

//    if((fpgaEvent->mcTracks[0].pt_ < 2) || (fpgaEvent->mcTracks[0].pt_ > 3)) continue;
//    if((fpgaEvent->mcTracks[0].pt_ < 9) || (fpgaEvent->mcTracks[0].pt_ > 10)) continue;

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

    ///////////////////////////////////
    // Duplicate Removal Reenactment //
    ///////////////////////////////////

    vector<FPGAEventTrack> secTracks[NSector];
    int nTrkWODupHybrid=0;

    unsigned int numTrk = fpgaEvent->tracks.size();
    vector<unsigned int> nComparisons[NSector];

    // Initialize duplicate arrays to false and fill sector tracks vectors
    for(unsigned int itrk=0; itrk<numTrk; itrk++) {
      FPGAEventTrack track = fpgaEvent->tracks.at(itrk);
      if(track.stubID_.size() < minStubs || abs(track.pt_)<Minpt || abs(track.pt_)>Maxpt) continue;
      secTracks[track.sector_].push_back(track);

    }

    //
    // Full sector removal
    //

    // Loop through sectors
    for(unsigned int nSec=0; nSec<NSector; nSec++) {
      uint numSecTrk = secTracks[nSec].size();
      h_nTrkSec2_tot->Fill((int)secTracks[nSec].size());
      //cout << "numSecTrk = " << numSecTrk << endl;

      // Initialize some objects
      unsigned int latency = 0;
      bool dupArray[numSecTrk]; // Array of duplicate flags for primary track
      for(unsigned int itrk=0; itrk<numSecTrk; itrk++) {
        dupArray[itrk] = false;
      }

      bool dupMap[numTrk][numTrk]; // 2D array of tracks being duplicates to other tracks
      for(unsigned int itrk=0; itrk<numTrk; itrk++) {
        for(unsigned int jtrk=0; jtrk<numTrk; jtrk++) {
          dupMap[itrk][jtrk] = false;
        }
      }

      // Fill some plots if there's only 1 track
      if(numSecTrk<=1) {
        h_secLatency->Fill(numSecTrk);
        if(i==0 && nSec==0) h_secComparisons_evt1->Fill(0); // For first event only
        h_secComparisons->Fill(0);
        continue;
      }

      // Loop through tracks to find number of stubs in common
      for(unsigned int itrk=0; itrk<numSecTrk-1; itrk++) { // numSecTrk-1 since last track has no other to compare to
        FPGAEventTrack track1 = secTracks[nSec].at(itrk);
          
//        // If primary track is a duplicate, it cannot veto any...move on
//        if(dupArray[itrk]==1) continue;
        bool primDup = dupArray[itrk];

        // Firmware analogies
        if(primDup==false) {
          latency++;
          nComparisons[nSec].push_back(secTracks[nSec].size()-itrk);
        }

        // Initialize structures for counting shared stubs
        int nStubP = 0;
        vector<int> nStubS(numSecTrk);
        vector<int> nShare(numSecTrk);

        // Get and count primary stubs
        std::map<int, int> stubsTrk1 = track1.stubID_;
        nStubP = stubsTrk1.size();

        // Loop over secondary tracks
        for(unsigned int jtrk=itrk+1; jtrk<numSecTrk; jtrk++) {
          FPGAEventTrack track2 = secTracks[nSec].at(jtrk);

  //        // Skip duplicate tracks
  //        if(dupArray[jtrk]==1) continue;

          // Get and count secondary stubs
          std::map<int, int> stubsTrk2 = track2.stubID_;
          nStubS[jtrk] = stubsTrk2.size();

          // Count shared stubs
          for(std::map<int, int>::iterator  st=stubsTrk1.begin(); st!=stubsTrk1.end(); st++) {
            if(stubsTrk2.find(st->first) != stubsTrk2.end()) {
              if(st->second == stubsTrk2[st->first]) nShare[jtrk]++;
            }
          }
        }

        // Tag duplicates
        for(unsigned int jtrk=itrk+1; jtrk<numSecTrk; jtrk++) {
          FPGAEventTrack track2 = secTracks[nSec].at(jtrk);

          // Skip duplicate tracks
          //if(dupArray[jtrk]==true) continue;
      
          // If removal condition is satisfied
          if((nStubP-nShare[jtrk] < minIndStubs) || (nStubS[jtrk]-nShare[jtrk] < minIndStubs)) {
            // Save map of which tracks are duplicates to which
            dupMap[itrk][jtrk] = true;
            dupMap[jtrk][itrk] = true;

            // Don't go through with removal if one of the tracks was already a duplicate
            if(primDup==true || dupArray[jtrk]==true) continue;
            if(track1.ichisq_ > track2.ichisq_) {
              dupArray[itrk] = true;
            }
            else if(track1.ichisq_ <= track2.ichisq_) {
              dupArray[jtrk] = true;
            }
            else cout << "Error: Didn't tag either track in duplicate pair." << endl;
          }
        }
      }

//      cout << endl << "Event " << i << " Sector " << nSec << endl;
//      for(uint j=0; j<numSecTrk; j++) {
//        for(uint k=0; k<numSecTrk; k++) {
//          cout << dupMap[j][k] << " ";
//        }
//        cout << endl;
//      }

      h_secLatency->Fill((int)latency);
      unsigned int nComparisons_tot = 0;
      for(unsigned int clock=0; clock<nComparisons[nSec].size(); clock++) {
        if(i==0 && nSec==0) h_secComparisons_evt1->Fill((int)nComparisons[nSec].at(clock)); // For first event only
        nComparisons_tot += nComparisons[nSec].at(clock);
      }
      h_secComparisons->Fill((int)nComparisons_tot);

      // For hybrid event plots
      int dupArrayHybrid[numSecTrk];
      int nTrkSecWODupHybrid=0;
      for(uint j=0; j<numSecTrk; j++) {
        dupArrayHybrid[j]=false;
        //cout << endl << "Event " << i << endl;
        for(uint k=j+1; k<numSecTrk; k++) {
          //cout << dupMap[j][k] << " ";
        }
        for(uint k=j+1; k<numSecTrk; k++) {
          if(dupMap[j][k]) {
            dupArrayHybrid[j]=true;
            for(uint m=0; m<numSecTrk; m++) {
              if(dupMap[j][m] && j!=m) dupMap[k][m]=true;
            }
            break;
          }
        }
        if(!dupArrayHybrid[j]) {
          nTrkWODupHybrid++;
          nTrkSecWODupHybrid++;
        }
      }
      for(uint j=0; j<numSecTrk; j++) {
        FPGAEventTrack track1 = secTracks[nSec].at(j);
        h_nStubTrkHybridResim->Fill(track1.stubID_.size());
        if(dupArrayHybrid[j]) continue;

        vector<pair<int,int>> trackStubs;

        for(std::map<int,int>::iterator iStub=track1.stubID_.begin(); iStub!=track1.stubID_.end(); iStub++) {
          trackStubs.push_back(pair<int,int>(iStub->first,iStub->second));
        }
        for(uint k=0; k<numSecTrk; k++) {
          if(dupMap[j][k] && j!=k) {
            FPGAEventTrack track2 = secTracks[nSec].at(k);
            for(std::map<int,int>::iterator iStub=track2.stubID_.begin(); iStub!=track2.stubID_.end(); iStub++) {
              pair<int,int> stubpair = pair<int,int>(iStub->first,iStub->second);
              if(std::find(trackStubs.begin(), trackStubs.end(), stubpair) == trackStubs.end()) trackStubs.push_back(stubpair);
            }
          }
        }
        h_nStubTrkHybridResim_after->Fill(trackStubs.size());
      }

      h_nTrkSecWODupHybrid->Fill(nTrkSecWODupHybrid);

    } // end sector loop

    // Hybrid event plot
    h_nTrkEvtWODupHybrid->Fill(nTrkWODupHybrid);



    //
    // Binned sector removal & firmware analogies
    //

    // Sector loop
    for(unsigned int nSec=0; nSec<NSector; nSec++) {
      uint numSecTrk = secTracks[nSec].size();

      // Initialize variables
      vector<FPGAEventTrack> binTracks[nEtaBins][nPhiBins];
      unsigned int binLatency[nEtaBins][nPhiBins] = {0};
      vector<unsigned int> nBinComparisons[nEtaBins][nPhiBins];

      // Sort tracks into bins
      for(unsigned int itrk=0; itrk<numSecTrk; itrk++) {
        FPGAEventTrack track = secTracks[nSec].at(itrk);
        int etaBin = min(nEtaBins-1.,max(0.,(track.eta_+2.4)/etaBinSize));
        double phi_temp = track.phi0_;
        if(Hourglass) phi_temp = phi_temp+2/M_PI;
        double phiReduced = 0;
        if(!Hourglass) {
          if(track.sector_<=12) phiReduced = phi_temp-track.sector_*2*M_PI/NSector;
          else phiReduced = phi_temp+(27-track.sector_)*2*M_PI/NSector;
        } else {
          if(track.sector_<=3) phiReduced = phi_temp-M_PI/NSector-track.sector_*2*M_PI/NSector;
          else phiReduced = phi_temp+M_PI/NSector+(8-track.sector_)*2*M_PI/NSector;
        }
        int phiBin = min(nPhiBins-1.,max(0.,(phiReduced)/phiBinSize));
        //cout << track.sector_ << " " << phi_temp << " " << phiReduced << " " << phiBin << endl;
        binTracks[etaBin][phiBin].push_back(track);
      }

      // Loop over bins for removal
      for(unsigned int nEtaBin=0; nEtaBin<nEtaBins; nEtaBin++) {
        for(unsigned int nPhiBin=0; nPhiBin<nPhiBins; nPhiBin++) {
          unsigned int numTrkBin = binTracks[nEtaBin][nPhiBin].size();
          //cout << numTrkBin << endl;
          //h_nTrkBin->Fill((int)binTracks[nEtaBin][nPhiBin].size());
          h_nTrkBin->Fill(numTrkBin);
          //cout << "numTrkBin = " << numTrkBin << endl;
          bool dupArray[numTrkBin]; // Array of duplicate flags for primary track
          for(unsigned int itrk=0; itrk<numTrkBin; itrk++) {
            dupArray[itrk] = false;
          }
          // If 1 or 0 tracks in bin, can take a shortcut here
          if(numTrkBin<=1) {
            h_binLatency->Fill(numTrkBin);
            if(i==0 && nSec==0 && nEtaBin==0 && nPhiBin==0) h_binComparisons_evt1->Fill(0); // For first event only
            h_binComparisons->Fill(0);
            continue;
          }

          // Loop over tracks and compare pairwise
          for(unsigned int itrk=0; itrk<numTrkBin-1; itrk++) { // numTrkBin-1 since last track has no other to compare to
            FPGAEventTrack track1 = binTracks[nEtaBin][nPhiBin].at(itrk);

            // If primary track is a duplicate, it cannot veto any...move on
            if(dupArray[itrk]==1) continue; // Need to continue for plots to be made
            //bool primDup = dupArray[itrk];

            // For firmware analogies
            //if(primDup==false) {
              binLatency[nEtaBin][nPhiBin]++;
              nBinComparisons[nEtaBin][nPhiBin].push_back(binTracks[nEtaBin][nPhiBin].size()-itrk);
            //}

            // Initialize comparison variables
            int nStubP = 0;
            vector<int> nStubS(numTrkBin);
            vector<int> nShare(numTrkBin);

            // Get and count primary stubs
            std::map<int, int> stubsTrk1 = track1.stubID_;
            nStubP = stubsTrk1.size();

            // Loop over secondary track
            for(unsigned int jtrk=itrk+1; jtrk<numTrkBin; jtrk++) {
              FPGAEventTrack track2 = binTracks[nEtaBin][nPhiBin].at(jtrk);

      //        // Skip duplicate tracks
      //        if(dupArray[jtrk]==1) continue;

              // Get and count secondary stubs
              std::map<int, int> stubsTrk2 = track2.stubID_;
              nStubS[jtrk] = stubsTrk2.size();

              // Count shared stubs
              for(std::map<int, int>::iterator  st=stubsTrk1.begin(); st!=stubsTrk1.end(); st++) {
                if(stubsTrk2.find(st->first) != stubsTrk2.end()) {
                  if(st->second == stubsTrk2[st->first]) nShare[jtrk]++;
                }
              }
            }

            // Tag duplicates
            for(unsigned int jtrk=itrk+1; jtrk<numTrkBin; jtrk++) {
              FPGAEventTrack track2 = binTracks[nEtaBin][nPhiBin].at(jtrk);
              // Skip duplicate tracks
              if(dupArray[jtrk]==true) continue;

              // Duplicate criteria
              if((nStubP-nShare[jtrk] < minIndStubs) || (nStubS[jtrk]-nShare[jtrk] < minIndStubs)) {

                // Don't go through with removal if one of the tracks was already a duplicate
                //if(primDup==true || dupArray[jtrk]==true) continue;
                if(track1.ichisq_ > track2.ichisq_) {
                  dupArray[itrk] = true;
                }
                else if(track1.ichisq_ <= track2.ichisq_) {
                  dupArray[jtrk] = true;
                }
                else cout << "Error: Didn't tag either track in duplicate pair." << endl;
              }
            }
          }

          h_binLatency->Fill((int)binLatency[nEtaBin][nPhiBin]);
          unsigned int nBinComparisons_tot = 0;
          for(unsigned int clock=0; clock<nBinComparisons[nEtaBin][nPhiBin].size(); clock++) {
            if(i==0 && nSec==0 && nEtaBin==0 && nPhiBin==0) h_binComparisons_evt1->Fill((int)nBinComparisons[nEtaBin][nPhiBin].at(clock)); // For first event only
            nBinComparisons_tot += nBinComparisons[nEtaBin][nPhiBin].at(clock);
          }
          h_binComparisons->Fill((int)nBinComparisons_tot);
        } // end phi bin loop
      } // end eta bin loop
    } // end sector loop

    ///////////////////////////
    // Grid Removal - Events //
    ///////////////////////////

    // Plots for 10 events with large numbers of tracks
    if(nGridEvts<10 && fpgaEvent->tracks.size()>=8) {
      for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
        FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
  //      if(track.eta_ > 1.0) continue;
        if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt) {

          double phiBin = (track.phi0_-2*M_PI/27*track.sector_)/(2*M_PI/9/50) + 1;
          phiBin = std::max(phiBin,0.);
          phiBin = std::min(phiBin,18.);

          double phiBinBig = (track.phi0_-2*M_PI/27*track.sector_)/(2*M_PI/9/50) + 9;
          phiBin = std::max(phiBin,0.);
          phiBin = std::min(phiBin,34.);

          // Correct phi domain to 0 < phi < 2*PI for global map  
          double phiGlobal;
          if(track.phi0_ < 0) {
            phiGlobal = 2*M_PI+track.phi0_;
          } else if(track.phi0_ > 2*M_PI) {
            phiGlobal = track.phi0_-2*M_PI;
          } else phiGlobal = track.phi0_;
          double phiBinGlobal = phiGlobal/(2*M_PI/9/50);

          double ptBin = 1/track.pt_*40+20;
          ptBin = std::max(ptBin,0.);
          ptBin = std::min(ptBin,39.);

          h_gridMap_manyTracks[nGridEvts]->Fill((int)ptBin,(int)phiBin);
          h_gridMapBig_manyTracks[nGridEvts]->Fill((int)ptBin,(int)phiBinBig);
          h_gridMapGlobal_manyTracks[nGridEvts]->Fill((int)phiBinGlobal,(int)ptBin);
          if(!track.duplicate_) h_gridMapGlobalWODup_manyTracks[nGridEvts]->Fill((int)phiBinGlobal,(int)ptBin);

        }
      }
      nGridEvts++;
    }

    ////////////////////////
    // Grid Removal Plots //
    ////////////////////////

    int gridSector=28;

    // Loop over tracks
    for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
      FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
//        if(track.seed_ != 1 && track.seed_ != 3 && track.seed_ != 5) continue;
//        if(abs(track.seed_) != 11 && abs(track.seed_) != 13) continue;
//        if(abs(track.seed_) != 21 && abs(track.seed_) != 22) continue;
      if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt) {

        nILTkt[track.seed_]++;
        h_iLTkt_tot->Fill(track.seed_);

        h_chisq_ichisq->Fill(track.chisq_, track.ichisq_);
        h_chisqDiff16->Fill(track.ichisq_/16.0-track.chisq_, track.ichisq_/16.0);
        h_chisqDiff16Zoom->Fill(track.ichisq_/16.0-track.chisq_, track.ichisq_/16.0);
//        h_chisqDiff16->Fill(track.ichisq_/2.0-track.chisq_, track.ichisq_/2.0);
//        h_chisqDiff16Zoom->Fill(track.ichisq_/2.0-track.chisq_, track.ichisq_/2.0);
        h_chisq_beforePD->Fill(track.chisq_);
        h_ichisq_beforePD->Fill(track.ichisq_);
        h_nStubTrk_premerge->Fill(track.stubIDpremerge_.size());
        if(!track.duplicate_) h_nStubTrk_prefit->Fill(track.stubIDprefit_.size());
        if(!track.duplicate_) h_nStubTrk_postfit->Fill(track.nStubpostfit_);
        h_nStubTrk->Fill(track.stubID_.size()); // After fit, important for hybrid project
        if(!track.duplicate_) h_nStubTrkWODup->Fill(track.stubID_.size());

        nTrkSec_tot[track.sector_]++;

        if(i==0) { //Only 1 event
          if(gridSector==28) gridSector=track.sector_;
//            if(track.sector_==gridSector) {

          double phiBin = (track.phi0_-2*M_PI/27*track.sector_)/(2*M_PI/9/50) + 9;
          phiBin = std::max(phiBin,0.);
          phiBin = std::min(phiBin,34.);

          double ptBin = 1/track.pt_*40+20;
          ptBin = std::max(ptBin,0.);
          ptBin = std::min(ptBin,39.);

          h_gridMap->Fill((int)ptBin,(int)phiBin);
          if(!track.duplicate_) h_gridMapWODup->Fill((int)ptBin,(int)phiBin);

          if(track.eta_ > 0) h_gridMapPlus->Fill((int)ptBin,(int)phiBin);
          if(track.eta_ > 0 && !track.duplicate_) h_gridMapPlusWODup->Fill((int)ptBin,(int)phiBin);
          if(track.eta_ < 0) h_gridMapMinus->Fill((int)ptBin,(int)phiBin);
          if(track.eta_ < 0 && !track.duplicate_) h_gridMapMinusWODup->Fill((int)ptBin,(int)phiBin);

//          }
        }

        if( track.eta_<=-1.0 ) nTrkSec_eta[track.sector_][0]++;
        if( track.eta_>-1.0 && track.eta_<=0) nTrkSec_eta[track.sector_][1]++;
        if( track.eta_<1.0 && track.eta_>0) nTrkSec_eta[track.sector_][2]++;
        if( track.eta_>=1.0 ) nTrkSec_eta[track.sector_][3]++;

        // Plotting stubids
//        std::map<int,std::vector<int>> Ids;
//        for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
//          bool exists=false;
//          for(std::vector<int>::iterator eachid=Ids->find(st->second).begin(); eachid!=Ids->find(st->second).end(); eachid++) {
//            if(eachid == st->second()) exists=true;
//          }
//          if(!exists) h_StubID_tot->Fill(st->second);
//          Ids.insert(
//        }
        if( track.eta_ > 1.4 ) {
          for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
            h_StubID_plus->Fill(st->second);
          }
        }
        if( track.eta_ < -1.4 ) {
          for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
            h_StubID_minus->Fill(st->second);
          }
        }
      }
    }

    /////////////////
    // Track Plots //
    /////////////////

    // Loop over tracks
    for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
      FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
//        if(track.seed_ != 1 && track.seed_ != 3 && track.seed_ != 5) continue;
//        if(abs(track.seed_) != 11 && abs(track.seed_) != 13) continue;
//        if(abs(track.seed_) != 21 && abs(track.seed_) != 22) continue;
      if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt) {

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
//          printf("chisq   = %f\n \n", track.chisq_);

        h_trkPt_tot->Fill(abs(track.pt_));
        h_trkPt_low_tot->Fill(abs(track.pt_));
        h_trkPhi_tot->Fill(track.phi0_);
        h_trkEta_tot->Fill(track.eta_);
        h_trkZ0_tot->Fill(track.z0_);

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
          h_ichisq_afterPD->Fill(track.ichisq_);

        }//end duplicates loop
      }//end if nStub >= nStubs
    }//end loop over tracks

    /////////////////////////////
    // Event Plots & Stub File //
    /////////////////////////////

    //make event level plots, but suppress events with no tracks passing requirements
    if(nTrkWODup_tot>0) {

      //total plots
      h_nTrkEvt_tot->Fill(nTrkWODup_tot+nDupTrk_tot);
      h_nTrkEvtWODup_tot->Fill(nTrkWODup_tot);
      h_nDupTrk_tot->Fill(nDupTrk_tot);
      h_NDupVSNTrkWO_tot->Fill(nTrkWODup_tot,nDupTrk_tot);
            
      //make plots for events with multiple unique tracks, output track information to .txt file
      if(nTrkWODup_tot>=2) {
        //Stub layer plots
        h_nStub_ghost->Fill(nStubTot);

        if(stubsTxtFile) stubid << Form("Event %i:",i);
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
            if(track.eta_ > -1.4) continue;
            if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt && !track.duplicate_) {
              //Output stubids to .txt file
              int nShareGhost=0;
              stubid << Form("\nTrack %i: ",ntrk);
              stubid << Form(" Phi=%f, Eta=%f, Z0=%f, Pt=%f, iChisq=%d, iLTkt=%i, Sector=%i", track.phi0_,track.eta_,track.z0_,abs(track.pt_),track.ichisq_,track.seed_,track.sector_);
              stubid << Form("\n    Stub ids: ");
              for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
                stubid << Form("  %i,%i",st->first,st->second);
                for(int ntrk2=0; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
                  FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);
                  if(!track2.duplicate_ && track2.stubID_.size() >= minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && ntrk!=ntrk2) {
                    if(track2.stubID_.find(st->first) != track2.stubID_.end()) {
                      if(st->second == track2.stubID_[st->first]) {
//                      if(st->second == track2.stubID_[st->first] && st->second != 63) {
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


    ///////////////////////////////////
    // Asymmetry Issues Plots (Long) //
    ///////////////////////////////////

    int asymGridPlus[40][450]={{0}};
    int asymGridMinus[40][450]={{0}};

    // Asymmetry investigations
    std::vector<std::pair<int,int>> stubs; //<sector,stubid (bottom 9 bits)
    // Bin first
    for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
      FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
      //if(track.duplicate_) continue;

      // Correct phi domain to 0 < phi < 2*PI for global map  
      double phiGlobal;
      if(track.phi0_ < 0) {
        phiGlobal = 2*M_PI+track.phi0_;
      } else if(track.phi0_ > 2*M_PI) {
        phiGlobal = track.phi0_-2*M_PI;
      } else phiGlobal = track.phi0_;
      double phiBinGlobal = phiGlobal/(2*M_PI/9/50);
      
      double ptBin = 1/track.pt_*40+20;
      ptBin = std::max(ptBin,0.);
      ptBin = std::min(ptBin,39.);
      
      if(track.eta_ > 1.4) asymGridPlus[(int)ptBin][(int)phiBinGlobal]++;
      if(track.eta_ < -1.4) asymGridMinus[(int)ptBin][(int)phiBinGlobal]++;

      // Fill stub histogram
      for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
        int stubsector = track.sector_+(st->second>>9)-1;
        if(stubsector==28) stubsector=0;
        if(stubsector==-1) stubsector=26;
        int stubid = st->second&((1<<9)-1);
        if(std::find(stubs.begin(),stubs.end(),std::pair<int,int>(stubsector,stubid))==stubs.end()) {
          stubs.push_back(std::pair<int,int>(stubsector,stubid));
          h_nStubLayer->Fill(st->first);
        }
      }
    }
    // Now chose a bin with many tracks
    int chosenphibinplus, chosenphibinminus = -1;
    int chosenptbinplus, chosenptbinminus = -1;
    for(int ptbin=0; ptbin<40; ptbin++) {
      for(int phibin=0; phibin<450; phibin++) {
        if(asymGridPlus[ptbin][phibin] >= 4 && chosenphibinplus == -1) {
          chosenphibinplus = phibin;
          chosenptbinplus = ptbin;
        }
        if(asymGridMinus[ptbin][phibin] >= 4 && chosenphibinminus == -1) {
          chosenphibinminus = phibin;
          chosenptbinminus = ptbin;
        }
      }
    }
    // Write out track information to .txt file.
    // Plus loop:
    if(stubsTxtFile) {
      asymstubidplus << Form("Event %i:",i);
      for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
        FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
        
        double phiGlobal;
        if(track.phi0_ < 0) {
          phiGlobal = 2*M_PI+track.phi0_;
        } else if(track.phi0_ > 2*M_PI) {
          phiGlobal = track.phi0_-2*M_PI;
        } else phiGlobal = track.phi0_;
        double phiBinGlobal = phiGlobal/(2*M_PI/9/50);
        
        double ptBin = 1/track.pt_*40+20;
        ptBin = std::max(ptBin,0.);
        ptBin = std::min(ptBin,39.);
        
        if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt) {
          //if(track.duplicate_) continue;
          if((int)phiBinGlobal!=chosenphibinplus || (int)ptBin!=chosenptbinplus) continue;
          
          //Output stubids to .txt file
          //First track information
          asymstubidplus << Form("\nTrack %i: ",ntrk);
          asymstubidplus << Form(" Phi=%f, Eta=%f, Z0=%f, Pt=%f, iChisq=%d, iLTkt=%i, Sector=%i", track.phi0_,track.eta_,track.z0_,abs(track.pt_),track.ichisq_,track.seed_,track.sector_);
          //Then stubids
          asymstubidplus << Form("\n    Stub ids: ");
          for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
            asymstubidplus << Form("  %i,%i",st->first,st->second);
            
            // Find other tracks with the same stub
            for(int ntrk2=0; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
              FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);
              if(track2.stubID_.size() >= minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && ntrk!=ntrk2) {
                
                double phiGlobal2;
                if(track2.phi0_ < 0) {
                  phiGlobal2 = 2*M_PI+track2.phi0_;
                } else if(track2.phi0_ > 2*M_PI) {
                  phiGlobal2 = track2.phi0_-2*M_PI;
                } else phiGlobal2 = track2.phi0_;
                double phiBinGlobal2 = phiGlobal2/(2*M_PI/9/50);
                
                double ptBin2 = 1/track2.pt_*40+20;
                ptBin2 = std::max(ptBin2,0.);
                ptBin2 = std::min(ptBin2,39.);
                
                if((int)phiBinGlobal2!=chosenphibinplus || (int)ptBin2!=chosenptbinplus) continue;
                //if(track2.duplicate_) continue;
                if(track2.stubID_.find(st->first) != track2.stubID_.end()) {
                  if(st->second == track2.stubID_[st->first]) {
                    asymstubidplus << "*";
                  }
                }
              }//end if duplicate
            }//end loop over secondary tracks
          }//end loop over stubs
        }//end if duplicate
      }//end loop over tracks
      asymstubidplus << "\n\n";
    }//end if stubsTxtFile
    
    // Minus loop:
    if(stubsTxtFile) {
      asymstubidplus << Form("Event %i:",i);
      for(int ntrk=0; ntrk<(int)fpgaEvent->tracks.size(); ntrk++){
        FPGAEventTrack track = fpgaEvent->tracks.at(ntrk);
        
        double phiGlobal;
        if(track.phi0_ < 0) {
          phiGlobal = 2*M_PI+track.phi0_;
        } else if(track.phi0_ > 2*M_PI) {
          phiGlobal = track.phi0_-2*M_PI;
        } else phiGlobal = track.phi0_;
        double phiBinGlobal = phiGlobal/(2*M_PI/9/50);
        
        double ptBin = 1/track.pt_*40+20;
        ptBin = std::max(ptBin,0.);
        ptBin = std::min(ptBin,39.);
        
        if(track.stubID_.size() >= minStubs && abs(track.pt_)>Minpt && abs(track.pt_)<Maxpt) {
          //if(track.duplicate_) continue;
          if((int)phiBinGlobal!=chosenphibinminus || (int)ptBin!=chosenptbinminus) continue;

          //Output stubids to .txt file
          //First track information
          asymstubidminus << Form("\nTrack %i: ",ntrk);
          asymstubidminus << Form(" Phi=%f, Eta=%f, Z0=%f, Pt=%f, iChisq=%d, iLTkt=%i, Sector=%i", track.phi0_,track.eta_,track.z0_,abs(track.pt_),track.ichisq_,track.seed_,track.sector_);
          //Then stubids
          asymstubidminus << Form("\n    Stub ids: ");
          for(std::map<int, int>::iterator  st=track.stubID_.begin(); st!=track.stubID_.end(); st++) {
            asymstubidminus << Form("  %i,%i",st->first,st->second);
            
            // Find other tracks with the same stub
            for(int ntrk2=0; ntrk2<(int)fpgaEvent->tracks.size(); ntrk2++) {
              FPGAEventTrack track2 = fpgaEvent->tracks.at(ntrk2);
              if(track2.stubID_.size() >= minStubs && abs(track2.pt_)>Minpt && abs(track2.pt_)<Maxpt && ntrk!=ntrk2) {

                double phiGlobal2;
                if(track2.phi0_ < 0) {
                  phiGlobal2 = 2*M_PI+track2.phi0_;
                } else if(track2.phi0_ > 2*M_PI) {
                  phiGlobal2 = track2.phi0_-2*M_PI;
                } else phiGlobal2 = track2.phi0_;
                double phiBinGlobal2 = phiGlobal2/(2*M_PI/9/50);
                
                double ptBin2 = 1/track2.pt_*40+20;
                ptBin2 = std::max(ptBin2,0.);
                ptBin2 = std::min(ptBin2,39.);

                if((int)phiBinGlobal2!=chosenphibinminus || (int)ptBin2!=chosenptbinminus) continue;
                //if(track2.duplicate_) continue;
                if(track2.stubID_.find(st->first) != track2.stubID_.end()) {
                  if(st->second == track2.stubID_[st->first]) {
                    asymstubidminus << "*";
                  }
                }
              }//end if duplicate
            }//end loop over secondary tracks
          }//end loop over stubs
        }//end if duplicate
      }//end loop over tracks
      asymstubidminus << "\n\n";
    }//end if stubsTxtFile

    ///////////////////
    // MC Track Loop //
    ///////////////////

    // Now look at the monte carlo distributions
    // only look at the events where we have a good FPGA track

    int nMCTrkEvt_tot=0;
    int nMCTrkSec_tot[27]={0};
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
        //nMCTrkSec_tot[(int)((mcTrack.phi_+TMath::Pi())/(TMath::TwoPi()/28))]++;

        // Try to see if this particle can be matched to a found track
        FPGAEventTrackMatch *bestMatch = 0;
        for(int fntrk=0; fntrk<(int)fpgaEvent->tracks.size(); fntrk++){
          FPGAEventTrack fpgaTrack = fpgaEvent->tracks.at(fntrk);
          if(fpgaTrack.duplicate_) continue;

          if(fpgaTrack.stubID_.size() >= minStubs) {
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
              h_matchDeltaPhi_4stub[iregion]->Fill(match->deltaPhi_);
              h_matchDeltaEta_4stub[iregion]->Fill(match->deltaEta_);
              h_matchDeltaR_4stub[iregion]  ->Fill(match->deltaR_);
              h_matchDeltaZ0_4stub[iregion] ->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_4stub[iregion] ->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(fpgaTrack.stubID_.size() == 5) {
              h_matchDeltaPhi_5stub[iregion]->Fill(match->deltaPhi_);
              h_matchDeltaEta_5stub[iregion]->Fill(match->deltaEta_);
              h_matchDeltaR_5stub[iregion]  ->Fill(match->deltaR_);
              h_matchDeltaZ0_5stub[iregion] ->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_5stub[iregion] ->Fill(match->deltaPt_/abs(mcTrack.pt_));
            }
            if(fpgaTrack.stubID_.size() >= 6) {
              h_matchDeltaPhi_6stub[iregion]->Fill(match->deltaPhi_);
              h_matchDeltaEta_6stub[iregion]->Fill(match->deltaEta_);
              h_matchDeltaR_6stub[iregion]  ->Fill(match->deltaR_);
              h_matchDeltaZ0_6stub[iregion] ->Fill(match->deltaZ0_);
              h_matchDeltaPtOPt_6stub[iregion] ->Fill(match->deltaPt_/abs(mcTrack.pt_));
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

    h_nMCTrkEvt_tot->Fill(nMCTrkEvt_tot);

    //////////////////
    // Sector Plots //
    //////////////////

    //Make total sector plots
    for(uint nSec=0; nSec<NSector; nSec++) {
      h_nTrkSec_eta1->Fill(min(25,nTrkSec_eta[nSec][0]));
      h_nTrkSec_eta2->Fill(min(25,nTrkSec_eta[nSec][1]));
      h_nTrkSec_eta3->Fill(min(25,nTrkSec_eta[nSec][2]));
      h_nTrkSec_eta4->Fill(min(25,nTrkSec_eta[nSec][3]));

      h_nTrkSec_tot->Fill(min(150,nTrkSec_tot[nSec]));
      h_nTrkSecTail_tot->Fill(nTrkSec_tot[nSec]);
      h_nTrkSecWODup_tot->Fill(min(50,nTrkSecWODup_tot[nSec]));
      h_nTrkSecWODupTail_tot->Fill(nTrkSecWODup_tot[nSec]);
      h_nMCTrkSec_tot->Fill(min(50,nMCTrkSec_tot[nSec]));
    }

  }//end loop over events


  stubid.close();
  asymstubidplus.close();
  asymstubidminus.close();
  
  // Write and close file
  histFile->Write();
  histFile->Close();

  //drop pointer to force reloading of input file
  fpgaEvent = 0;

}

