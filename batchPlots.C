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

void feb20Plots(TString rootFile1="ttbar_newtune_ana.root") {
  make2DHist(rootFile1, "h_chisq_ichisq", "", "chisq_correlation");
  makeSingleHist(rootFile1, "h_ichisq_beforePD", "", "ichisq_beforePD");
  makeSingleHist(rootFile1, "h_ichisq_afterPD", "", "ichisq_afterPD");
  makeSingleHist(rootFile1, "h_chisq_beforePD", "", "chisq_beforePD");
  makeSingleHist(rootFile1, "h_chisq_afterPD", "", "chisq_afterPD");

}

void makeChisqPlots(TString rootFile1="ana_ttbar_11.root") {

  makeSingleHist(rootFile1, "h_chisq_beforePD", "", "h_chisq_before", false);
  makeSingleHist(rootFile1, "h_ichisq_beforePD", "", "h_ichisq_before", false);
  makeSingleHist(rootFile1, "h_chisq_afterPD", "", "h_chisq_after", false);
  makeSingleHist(rootFile1, "h_ichisq_afterPD", "", "h_ichisq_after", false);
  makeDoubleHist(rootFile1, rootFile1, "h_ichisq_beforePD", "h_ichisq_afterPD", "h_ichisq_double");
  makeDoubleHist(rootFile1, rootFile1, "h_chisq_beforePD", "h_chisq_afterPD", "h_chisq_double");
  make2DHist(rootFile1, "h_chisq_ichisq", "", "h_ichiVSchi");
  make2DHist(rootFile1, "h_chisqDiff16", "", "h_ichi16VSDiff");
  make2DHist(rootFile1, "h_chisqDiff16Zoom", "", "h_ichi16VSDiffZoom");
  for(int slice=1; slice<13; slice++) {
    makeXSliceHist(rootFile1, "h_chisqDiff16Zoom", "", slice, slice, "h_ichi16VSDiffZoom_Slice"+to_string(slice));
  }

}

void makeEtaAsym() {

  makeSingleHist("ana_ttbar_ichi.root","h_trkEta_tot","","trkEta_before");
  makeSingleHist("ana_ttbar_ichi.root","h_trkEtaWODup_tot","","trkEta_ichi_after");
  makeSingleHist("ana_ttbar_grid.root","h_trkEtaWODup_tot","","trkEta_grid_after");

}

void makeResPlots(TString rootFile="ana_myTest.root") {

//  resolutionPlot(rootFile, "h_matchDeltaPhi_5stub", "", "PhiRes_5stub");

  combinationResolution(rootFile,"","Phi", "PhiResVsNStub");
  combinationResolution(rootFile,"","Eta","EtaResVsNStub");
  combinationResolution(rootFile,"","Z0","Z0ResVsNStub");
  combinationResolution(rootFile,"","PtOPt","PtOPtResVsNStub");
//  
//    //combinationResolution(rootFile,"Phi","eta","PhiResVsNStub_eta");
//    //combinationResolution(rootFile,"Eta","eta","EtaResVsNStub_eta");
//    //combinationResolution(rootFile,"Z0","eta","Z0ResVsNStub_eta");
//    //combinationResolution(rootFile,"PtOPt","eta","PtOPtResVsNStub_eta");

}

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

void StandardMuPlots() {
  makeDoubleHist("ana_Mu_10000.root","ana_Mu_10000.root","h_nTrkEvt_tot","h_nTrkEvtWODup_tot","Mu50 PD Performance",true,true);
  makeSingleHist("ana_MuMinus_ichi_adj.root","h_nDupTrk_tot","","nDupTrack");
  makeSingleHist("ana_MuMinus_ichi_adj.root","h_nMCTrkEvt_tot","","nMCTrk");
  makeDoubleHist("ana_MuMinus_ichi_adj.root","ana_MuMinus_ichi_adj.root","h_ichisq_beforePD","h_ichisq_afterPD","iChisquare Before&After",false,true);
}

void StandardMuPlotsMCpt() {
  makeDoubleHist("myTest_2to3_ana.root","myTest_2to3_ana.root","h_nTrkEvt_tot","h_nTrkEvtWODup_tot","Mu10_PD_Performance_2to3",true,true);
  makeDoubleHist("myTest_9to10_ana.root","myTest_9to10_ana.root","h_nTrkEvt_tot","h_nTrkEvtWODup_tot","Mu10_PD_Performance_9to10",true,true);
}

void StandardTTbarPlots() {
  makeTripleHist("output_ttbar_D17_PU200_noPD_300.root", "output_ttbar_D17_PU200_nstubPD_300.root", "output_ttbar_D17_PU200_condensed_300.root", "eff_eta", "eff_eta", "eff_eta", "", "ttbar_efficiencies");
  makeTripleHist("output_ttbar_D17_PU200_condensed.root", "output_ttbar_D17_PU200_truncate.root", "output_ttbar_D17_PU200_extreme.root", "eff_eta", "eff_eta", "eff_eta", "", "ttbar_ichi_compare");
}

void GridPlots() {
  make2DHist("ana_MuMinus_grid_adj.root","h_gridMapBig_manyTracks_0","","Mu10_gridMapBig_evt0");
  make2DHist("ana_MuMinus_grid_adj.root","h_gridMapBig_manyTracks_1","","Mu10_gridMapBig_evt1");
  make2DHist("ana_MuMinus_grid_adj.root","h_gridMapBig_manyTracks_2","","Mu10_gridMapBig_evt2");
  make2DHist("ana_MuMinus_grid_adj.root","h_gridMapBig_manyTracks_3","","Mu10_gridMapBig_evt3");
  make2DHist("ana_MuMinus_grid_adj.root","h_gridMapBig_manyTracks_4","","Mu10_gridMapBig_evt4");
  make2DHist("ana_MuMinus_grid_adj.root","h_gridMapBig_manyTracks_5","","Mu10_gridMapBig_evt5");

  make2DHist("ana_ttbar_grid_adj.root","h_gridMap","","ttbar_gridMapSmall");
  make2DHist("ana_ttbar_grid_adj.root","h_gridMapBig","","ttbar_gridMapBig");
  make2DHist("ana_ttbar_grid_adj.root","h_gridMapBigWODup","","ttbar_gridMapBigWODup");
  make2DHist("ana_ttbar_grid_adj.root","h_gridMapGlobal_manyTracks_0","","ttbar_gridMapGlobal");
  make2DHist("ana_ttbar_grid_adj.root","h_gridMapGlobalWODup_manyTracks_0","","ttbar_gridMapGlobalWODup");
}
