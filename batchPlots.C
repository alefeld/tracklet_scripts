#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <THStack.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <iostream>
#include <fstream>
#include <iomanip>


void allPlots(TString rootFileName="totalMuFPGAAnalysis.root", TString directory="") {
  TFile *batchRootFile = new TFile(rootFileName);
  TDirectory *batchD = (TDirectory*)batchRootFile->GetDirectory(directory);

  for(auto&& keyAsObj : *batchD->GetListOfKeys()) {
    auto key = (TKey*) keyAsObj;
//    cout << key->GetName() << " " << key->GetClassName() << endl;
    if(strncmp(key->GetClassName(),"TH1F",4)==0) makeSingleHist(rootFileName, key->GetName(), directory);
    if(strncmp(key->GetClassName(),"TH2F",4)==0) make2DHist(rootFileName, key->GetName(), directory);
  }
}

void allPlotsCSV(TString rootFileName="totalMuFPGAAnalysis.root", TString directory="") {
  TFile *batchRootFile = new TFile(rootFileName);
  TDirectory *batchD = (TDirectory*)batchRootFile->GetDirectory(directory);

  for(auto&& keyAsObj : *batchD->GetListOfKeys()) {
    auto key = (TKey*) keyAsObj;
//    cout << key->GetName() << " " << key->GetClassName() << endl;
    if(strncmp(key->GetClassName(),"TH1F",4)==0) makeSingleHist(rootFileName, key->GetName(), directory);
    if(strncmp(key->GetClassName(),"TH2F",4)==0) make2DHistCSV(rootFileName, key->GetName(), directory);
  }
}

void makeProfilePlots(TString rootFile="totalMu10FPGAAnalysis.root") {

    makeProfile(rootFile, "h_NTrkVsmcTrkPt", "h_NTrkWODupVsmcTrkPt", "mcEtaLess1", "Profile_NTrkvsmcTrkPt_EtaLess1");
    makeProfile(rootFile, "h_NTrkVsmcTrkPt", "h_NTrkWODupVsmcTrkPt", "mcEta1to1.7", "Profile_NTrkvsmcTrkPt_Eta1to1.7");
    makeProfile(rootFile, "h_NTrkVsmcTrkPt", "h_NTrkWODupVsmcTrkPt", "mcEtaMore1.7", "Profile_NTrkvsmcTrkPt_EtaMore1.7");

}

void makeTriplePlots() {

    makeTripleHistSpecial("totalMu10_samesecFPGAAnalysis.root", "totalMu10_samesecFPGAAnalysis.root", "totalMu10_adjsecFPGAAnalysis.root", "h_nTrkEvt_tot", "h_nTrkWODup_tot", "h_nTrkEvtWODup_tot", "total", "totalMu10nTrkEvt");
    makeTripleHistSpecial("totalMu150_samesecFPGAAnalysis.root", "totalMu150_samesecFPGAAnalysis.root", "totalMu150_adjsecFPGAAnalysis.root", "h_nTrkEvt_tot", "h_nTrkWODup_tot", "h_nTrkEvtWODup_tot", "total", "totalMu150nTrkEvt");
    makeTripleHistSpecial("ttbar_samesecFPGAAnalysis.root", "ttbar_samesecFPGAAnalysis.root", "ttbar_adjsecFPGAAnalysis.root", "h_nTrkSec_tot", "h_nTrkSecWODup_tot", "h_nTrkSecWODup_tot", "total", "ttbar_nTrkSec_tot");
    makeTripleHistSpecial("ttbar_samesecFPGAAnalysis.root", "ttbar_samesecFPGAAnalysis.root", "ttbar_adjsecFPGAAnalysis.root", "h_nTrkSecTail_tot", "h_nTrkSecWODupTail_tot", "h_nTrkSecWODupTail_tot", "total", "ttbar_nTrkSecTail_tot");

}

void makeChisqPlots(TString rootFile1="mu_standalone_Ana.root",TString rootFile2="mu_CMSSW_Ana.root") {

  makeSingleHist(rootFile1, "h_chisq_afterPD", "h_chisq_standalone", false);
  makeSingleHist(rootFile2, "h_chisq", "h_chisq_CMSSW", false);
  makeDoubleHist(rootFile1, rootFile2, "h_chisq_afterPD", "h_chisq", "h_chisq_compare", false);
  makeSingleHist(rootFile1, "h_chisq_afterPD", "h_chisq_standalone_log", true);
  makeSingleHist(rootFile2, "h_chisq", "h_chisq_CMSSW_log", true);
  makeDoubleHist(rootFile1, rootFile2, "h_chisq_afterPD", "h_chisq", "h_chisq_compare_log", true);

}

//void makeResPlots(TString rootFile="totalMuFPGAAnalysis.root") {

//    resolutionPlots(rootFile,"h_matchDeltaPhi_AllM","mcEtaLess1","DeltaPhi_allMatch_etaLess1");
//    resolutionPlots(rootFile,"h_matchDeltaEta_AllM","mcEtaLess1","DeltaEta_allMatch_etaLess1");
//    resolutionPlots(rootFile,"h_matchDeltaR_AllM","mcEtaLess1","DeltaR_allMatch_etaLess1");
//    resolutionPlots(rootFile,"h_matchDeltaPtOPt_AllM","mcEtaLess1","DeltaPtOPt_allMatch_etaLess1");
//    resolutionPlots(rootFile,"h_matchDeltaPhi_AllM","mcEta1to1.7","DeltaPhi_allMatch_eta1to1.7");
//    resolutionPlots(rootFile,"h_matchDeltaEta_AllM","mcEta1to1.7","DeltaEta_allMatch_eta1to1.7");
//    resolutionPlots(rootFile,"h_matchDeltaR_AllM","mcEta1to1.7","DeltaR_allMatch_eta1to1.7");
//    resolutionPlots(rootFile,"h_matchDeltaPtOPt_AllM","mcEta1to1.7","DeltaPtOPt_allMatch_eta1to1.7");
//    resolutionPlots(rootFile,"h_matchDeltaPhi_AllM","mcEtaMore1.7","DeltaPhi_allMatch_etaMore1.7");
//    resolutionPlots(rootFile,"h_matchDeltaEta_AllM","mcEtaMore1.7","DeltaEta_allMatch_etaMore1.7");
//    resolutionPlots(rootFile,"h_matchDeltaR_AllM","mcEtaMore1.7","DeltaR_allMatch_etaMore1.7");
//    resolutionPlots(rootFile,"h_matchDeltaPtOPt_AllM","mcEtaMore1.7","DeltaPtOPt_allMatch_etaMore1.7");

  //resolutionPlots(rootFile,"h_matchDeltaPhi_AllM","barrelOnly","DeltaPhi_allMatch_barrelOnly");
  //resolutionPlots(rootFile,"h_matchDeltaEta_AllM","barrelOnly","DeltaEta_allMatch_barrelOnly");
  //resolutionPlots(rootFile,"h_matchDeltaR_AllM","barrelOnly","DeltaR_allMatch_barrelOnly");
  //resolutionPlots(rootFile,"h_matchDeltaPtOPt_AllM","barrelOnly","DeltaPtOPt_allMatch_barrelOnly");
  //resolutionPlots(rootFile,"h_matchDeltaPhi_AllM","diskOnly","DeltaPhi_allMatch_diskonly");
  //resolutionPlots(rootFile,"h_matchDeltaEta_AllM","diskOnly","DeltaEta_allMatch_diskonly");
  //resolutionPlots(rootFile,"h_matchDeltaR_AllM","diskOnly","DeltaR_allMatch_diskonly");
  //resolutionPlots(rootFile,"h_matchDeltaPtOPt_AllM","diskOnly","DeltaPtOPt_allMatch_diskonly");
  //resolutionPlots(rootFile,"h_matchDeltaPhi_AllM","overlap","DeltaPhi_allMatch_overlap");
  //resolutionPlots(rootFile,"h_matchDeltaEta_AllM","overlap","DeltaEta_allMatch_overlap");
  //resolutionPlots(rootFile,"h_matchDeltaR_AllM","overlap","DeltaR_allMatch_overlap");
  //resolutionPlots(rootFile,"h_matchDeltaPtOPt_AllM","overlap","DeltaPtOPt_allMatch_overlap"); 
  /*
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS4","mcEtaLess1","DeltaZ0_nS4_mcEtaLess1");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS5","mcEtaLess1","DeltaZ0_nS5_mcEtaLess1");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS6","mcEtaLess1","DeltaZ0_nS6_mcEtaLess1");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS4","mcEta1to1.7","DeltaZ0_nS4_mcEta1to1.7");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS5","mcEta1to1.7","DeltaZ0_nS5_mcEta1to1.7");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS6","mcEta1to1.7","DeltaZ0_nS6_mcEta1to1.7");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS4","mcEtaMore1.7","DeltaZ0_nS4_mcEtaMore1.7");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS5","mcEtaMore1.7","DeltaZ0_nS5_mcEtaMore1.7");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS6","mcEtaMore1.7","DeltaZ0_nS6_mcEtaMore1.7");
*/
  /*
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS4","barrelOnly","DeltaZ0_nS4_diskOnly");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS5","barrelOnly","DeltaZ0_nS5_diskOnly");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS6","barrelOnly","DeltaZ0_nS6_diskOnly");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS4","diskOnly","DeltaZ0_nS4_diskOnly");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS5","diskOnly","DeltaZ0_nS5_diskOnly");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS6","diskOnly","DeltaZ0_nS6_diskOnly");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS4","overlap","DeltaZ0_nS4_overlap");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS5","overlap","DeltaZ0_nS5_overlap");
    resolutionPlots(rootFile,"h_matchDeltaZ0_nS6","overlap","DeltaZ0_nS6_overlap");
  */
//    combinationResolution(rootFile,"Phi","tracklet","PhiResVsNStub_trklet");
//    combinationResolution(rootFile,"Eta","tracklet","EtaResVsNStub_trklet");
//    combinationResolution(rootFile,"Z0","tracklet","Z0ResVsNStub_trklet");
//    combinationResolution(rootFile,"PtOPt","tracklet","PtOPtResVsNStub_trklet");
//  
//    //combinationResolution(rootFile,"Phi","eta","PhiResVsNStub_eta");
//    //combinationResolution(rootFile,"Eta","eta","EtaResVsNStub_eta");
//    //combinationResolution(rootFile,"Z0","eta","Z0ResVsNStub_eta");
//    //combinationResolution(rootFile,"PtOPt","eta","PtOPtResVsNStub_eta");



//}

void makeDoubleResPlots(TString rootFile1="plusMu10FPGAAnalysis.root", TString rootFile2="minusMu10FPGAAnalysis.root") {

    doubleResPlots(rootFile1,rootFile2,"PtOPt","total","double_mcTrkPt");
    doubleResPlots(rootFile1,rootFile2,"Phi","total","double_mcTrkPhi");
    doubleResPlots(rootFile1,rootFile2,"Z0","total","double_mcTrkZ0");
    doubleResPlots(rootFile1,rootFile2,"Eta","total","double_mcTrkEta");

}

void makeEffPlots(TString rootFile="totalMuFPGAAnalysis.root") {

    efficiencyPlot(rootFile,"Phi","mcEtaLess1","PhiEff_mcEtaLess1");
    efficiencyPlot(rootFile,"Eta","mcEtaLess1","EtaEff_mcEtaLess1");
    efficiencyPlot(rootFile,"Z0","mcEtaLess1","Z0Eff_mcEtaLess1");
    efficiencyPlot(rootFile,"Pt","mcEtaLess1","PtEff_mcEtaLess1");

    efficiencyPlot(rootFile,"Phi","mcEta1to1.7","PhiEff_mcEta1to1.7");
    efficiencyPlot(rootFile,"Eta","mcEta1to1.7","EtaEff_mcEta1to1.7");
    efficiencyPlot(rootFile,"Z0","mcEta1to1.7","Z0Eff_mcEta1to1.7");
    efficiencyPlot(rootFile,"Pt","mcEta1to1.7","PtEff_mcEta1to1.7");

    efficiencyPlot(rootFile,"Phi","mcEtaMore1.7","PhiEff_mcEtaMore1.7");
    efficiencyPlot(rootFile,"Eta","mcEtaMore1.7","EtaEff_mcEtaMore1.7");
    efficiencyPlot(rootFile,"Z0","mcEtaMore1.7","Z0Eff_mcEtaMore1.7");
    efficiencyPlot(rootFile,"Pt","mcEtaMore1.7","PtEff_mcEtaMore1.7");

    efficiencyPlot(rootFile,"Phi","total","PhiEff_total");
    efficiencyPlot(rootFile,"Eta","total","EtaEff_total");
    efficiencyPlot(rootFile,"Z0","total","Z0Eff_total");
    efficiencyPlot(rootFile,"Pt","total","PtEff_total");

}

void makeMultiplicityPlots(TString rootFile="mu10Analysis.root") {

  makeBarrelStubHist(rootFile,"");
  makeBarrelStubHist(rootFile,"_ghost");
  makeDiskStubHist(rootFile,"F","");
  makeDiskStubHist(rootFile,"F","_ghost");
  makeDiskStubHist(rootFile,"B","");
  makeDiskStubHist(rootFile,"B","_ghost");
  makeBarrelTrackHist(rootFile,"");
  makeBarrelTrackHist(rootFile,"_ghost");
  makeDiskTrackHist(rootFile,"F","");
  makeDiskTrackHist(rootFile,"F","_ghost");
  makeDiskTrackHist(rootFile,"B","");
  makeDiskTrackHist(rootFile,"B","_ghost");
  makeSingleHist(rootFile,"h_iLTkt_tot");
  makeSingleHist(rootFile,"h_iLTktWODup_tot");
  makeSingleHist(rootFile,"h_iLTkt_ghost");

}


void DecemberFixedPlots() {
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_nTrkEvt_tot","barrel_fixed/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_nDupTrk_tot","barrel_fixed/h_nDupTrk_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_nTrkEvtWODup_tot","barrel_fixed/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_nMCTrkEvt_tot","barrel_fixed/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_trkPt_tot","barrel_fixed/h_trkPt_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_trkPhi_tot","barrel_fixed/h_trkPhi_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_trkEta_tot","barrel_fixed/h_trkEta_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_trkPtWODup_tot","barrel_fixed/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_trkPhiWODup_tot","barrel_fixed/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_barrel_fixed_ana.root","h_trkEtaWODup_tot","barrel_fixed/h_trkEtaWODup_tot");

  makeSingleHist("ttbar140_disk_fixed_ana.root","h_nTrkEvt_tot","disk_fixed/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_nDupTrk_tot","disk_fixed/h_nDupTrk_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_nTrkEvtWODup_tot","disk_fixed/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_nMCTrkEvt_tot","disk_fixed/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_trkPt_tot","disk_fixed/h_trkPt_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_trkPhi_tot","disk_fixed/h_trkPhi_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_trkEta_tot","disk_fixed/h_trkEta_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_trkPtWODup_tot","disk_fixed/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_trkPhiWODup_tot","disk_fixed/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_disk_fixed_ana.root","h_trkEtaWODup_tot","disk_fixed/h_trkEtaWODup_tot");

  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_nTrkEvt_tot","hybrid_fixed/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_nDupTrk_tot","hybrid_fixed/h_nDupTrk_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_nTrkEvtWODup_tot","hybrid_fixed/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_nMCTrkEvt_tot","hybrid_fixed/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_trkPt_tot","hybrid_fixed/h_trkPt_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_trkPhi_tot","hybrid_fixed/h_trkPhi_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_trkEta_tot","hybrid_fixed/h_trkEta_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_trkPtWODup_tot","hybrid_fixed/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_trkPhiWODup_tot","hybrid_fixed/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_hybrid_fixed_ana.root","h_trkEtaWODup_tot","hybrid_fixed/h_trkEtaWODup_tot");

  makeSingleHist("ttbar140_all_fixed_ana.root","h_nTrkEvt_tot","all_fixed/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_nDupTrk_tot","all_fixed/h_nDupTrk_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_nTrkEvtWODup_tot","all_fixed/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_nMCTrkEvt_tot","all_fixed/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_trkPt_tot","all_fixed/h_trkPt_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_trkPhi_tot","all_fixed/h_trkPhi_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_trkEta_tot","all_fixed/h_trkEta_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_trkPtWODup_tot","all_fixed/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_trkPhiWODup_tot","all_fixed/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_all_fixed_ana.root","h_trkEtaWODup_tot","all_fixed/h_trkEtaWODup_tot");

  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_nTrkEvt_tot","hybridplus_fixed/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_nDupTrk_tot","hybridplus_fixed/h_nDupTrk_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_nTrkEvtWODup_tot","hybridplus_fixed/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_nMCTrkEvt_tot","hybridplus_fixed/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_trkPt_tot","hybridplus_fixed/h_trkPt_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_trkPhi_tot","hybridplus_fixed/h_trkPhi_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_trkEta_tot","hybridplus_fixed/h_trkEta_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_trkPtWODup_tot","hybridplus_fixed/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_trkPhiWODup_tot","hybridplus_fixed/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_hybridplus_fixed_ana.root","h_trkEtaWODup_tot","hybridplus_fixed/h_trkEtaWODup_tot");
}

void DecemberOldPlots() {
  makeSingleHist("ttbar140_barrel_old_ana.root","h_nTrkEvt_tot","barrel_old/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_nDupTrk_tot","barrel_old/h_nDupTrk_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_nTrkEvtWODup_tot","barrel_old/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_nMCTrkEvt_tot","barrel_old/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_trkPt_tot","barrel_old/h_trkPt_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_trkPhi_tot","barrel_old/h_trkPhi_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_trkEta_tot","barrel_old/h_trkEta_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_trkPtWODup_tot","barrel_old/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_trkPhiWODup_tot","barrel_old/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_barrel_old_ana.root","h_trkEtaWODup_tot","barrel_old/h_trkEtaWODup_tot");

  makeSingleHist("ttbar140_disk_old_ana.root","h_nTrkEvt_tot","disk_old/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_nDupTrk_tot","disk_old/h_nDupTrk_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_nTrkEvtWODup_tot","disk_old/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_nMCTrkEvt_tot","disk_old/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_trkPt_tot","disk_old/h_trkPt_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_trkPhi_tot","disk_old/h_trkPhi_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_trkEta_tot","disk_old/h_trkEta_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_trkPtWODup_tot","disk_old/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_trkPhiWODup_tot","disk_old/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_disk_old_ana.root","h_trkEtaWODup_tot","disk_old/h_trkEtaWODup_tot");

  makeSingleHist("ttbar140_hybrid_old_ana.root","h_nTrkEvt_tot","hybrid_old/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_nDupTrk_tot","hybrid_old/h_nDupTrk_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_nTrkEvtWODup_tot","hybrid_old/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_nMCTrkEvt_tot","hybrid_old/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_trkPt_tot","hybrid_old/h_trkPt_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_trkPhi_tot","hybrid_old/h_trkPhi_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_trkEta_tot","hybrid_old/h_trkEta_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_trkPtWODup_tot","hybrid_old/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_trkPhiWODup_tot","hybrid_old/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_hybrid_old_ana.root","h_trkEtaWODup_tot","hybrid_old/h_trkEtaWODup_tot");

  makeSingleHist("ttbar140_all_old_ana.root","h_nTrkEvt_tot","all_old/h_nTrkEvt_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_nDupTrk_tot","all_old/h_nDupTrk_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_nTrkEvtWODup_tot","all_old/h_nTrkEvtWODup_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_nMCTrkEvt_tot","all_old/h_nMCTrkEvt_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_trkPt_tot","all_old/h_trkPt_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_trkPhi_tot","all_old/h_trkPhi_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_trkEta_tot","all_old/h_trkEta_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_trkPtWODup_tot","all_old/h_trkPtWODup_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_trkPhiWODup_tot","all_old/h_trkPhiWODup_tot");
  makeSingleHist("ttbar140_all_old_ana.root","h_trkEtaWODup_tot","all_old/h_trkEtaWODup_tot");

}

void DoubleDecemberPlots() {
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_nTrkEvt_tot","h_nTrkEvt_tot","barrel_double/h_nTrkEvt");
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_nTrkEvtWODup_tot","h_nTrkEvtWODup_tot","barrel_double/h_nTrkEvtWODup");
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_nDupTrk_tot","h_nDupTrk_tot","barrel_double/h_nDupTrk");
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_nMCTrkEvt_tot","h_nMCTrkEvt_tot","barrel_double/h_nMCTrkEvt");
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_trkPt_tot","h_trkPt_tot","barrel_double/h_trkPt",false,true);
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_trkPhi_tot","h_trkPhi_tot","barrel_double/h_trkPhi",false,true);
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_trkEta_tot","h_trkEta_tot","barrel_double/h_trkEta",false,true);
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_trkPtWODup_tot","h_trkPtWODup_tot","barrel_double/h_trkPtWODup",false,true);
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_trkPhiWODup_tot","h_trkPhiWODup_tot","barrel_double/h_trkPhiWODup",false,true);
  makeDoubleHist("ttbar140_barrel_fixed_ana.root","ttbar140_barrel_old_ana.root","h_trkEtaWODup_tot","h_trkEtaWODup_tot","barrel_double/h_trkEtaWODup",false,true);

  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_nTrkEvt_tot","h_nTrkEvt_tot","disk_double/h_nTrkEvt");
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_nTrkEvtWODup_tot","h_nTrkEvtWODup_tot","disk_double/h_nTrkEvtWODup");
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_nDupTrk_tot","h_nDupTrk_tot","disk_double/h_nDupTrk");
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_nMCTrkEvt_tot","h_nMCTrkEvt_tot","disk_double/h_nMCTrkEvt");
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_trkPt_tot","h_trkPt_tot","disk_double/h_trkPt",false,true);
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_trkPhi_tot","h_trkPhi_tot","disk_double/h_trkPhi",false,true);
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_trkEta_tot","h_trkEta_tot","disk_double/h_trkEta",false,true);
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_trkPtWODup_tot","h_trkPtWODup_tot","disk_double/h_trkPtWODup",false,true);
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_trkPhiWODup_tot","h_trkPhiWODup_tot","disk_double/h_trkPhiWODup",false,true);
  makeDoubleHist("ttbar140_disk_fixed_ana.root","ttbar140_disk_old_ana.root","h_trkEtaWODup_tot","h_trkEtaWODup_tot","disk_double/h_trkEtaWODup",false,true);

  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_nTrkEvt_tot","h_nTrkEvt_tot","hybrid_double/h_nTrkEvt");
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_nTrkEvtWODup_tot","h_nTrkEvtWODup_tot","hybrid_double/h_nTrkEvtWODup");
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_nDupTrk_tot","h_nDupTrk_tot","hybrid_double/h_nDupTrk");
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_nMCTrkEvt_tot","h_nMCTrkEvt_tot","hybrid_double/h_nMCTrkEvt");
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_trkPt_tot","h_trkPt_tot","hybrid_double/h_trkPt",false,true);
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_trkPhi_tot","h_trkPhi_tot","hybrid_double/h_trkPhi",false,true);
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_trkEta_tot","h_trkEta_tot","hybrid_double/h_trkEta",false,true);
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_trkPtWODup_tot","h_trkPtWODup_tot","hybrid_double/h_trkPtWODup",false,true);
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_trkPhiWODup_tot","h_trkPhiWODup_tot","hybrid_double/h_trkPhiWODup",false,true);
  makeDoubleHist("ttbar140_hybrid_fixed_ana.root","ttbar140_hybrid_old_ana.root","h_trkEtaWODup_tot","h_trkEtaWODup_tot","hybrid_double/h_trkEtaWODup",false,true);

  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_nTrkEvt_tot","h_nTrkEvt_tot","all_double/h_nTrkEvt");
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_nTrkEvtWODup_tot","h_nTrkEvtWODup_tot","all_double/h_nTrkEvtWODup");
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_nDupTrk_tot","h_nDupTrk_tot","all_double/h_nDupTrk");
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_nMCTrkEvt_tot","h_nMCTrkEvt_tot","all_double/h_nMCTrkEvt");
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_trkPt_tot","h_trkPt_tot","all_double/h_trkPt",false,true);
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_trkPhi_tot","h_trkPhi_tot","all_double/h_trkPhi",false,true);
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_trkEta_tot","h_trkEta_tot","all_double/h_trkEta",false,true);
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_trkPtWODup_tot","h_trkPtWODup_tot","all_double/h_trkPtWODup",false,true);
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_trkPhiWODup_tot","h_trkPhiWODup_tot","all_double/h_trkPhiWODup",false,true);
  makeDoubleHist("ttbar140_all_fixed_ana.root","ttbar140_all_old_ana.root","h_trkEtaWODup_tot","h_trkEtaWODup_tot","all_double/h_trkEtaWODup",false,true);

}