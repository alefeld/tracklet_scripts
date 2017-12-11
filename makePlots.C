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

void makeSingleHist(TString rootFileName="mu10Analysis.root", TString hname="h_iLTkt", TString outputname="", bool log=false) {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    TString outputdir = "Plots/";
    TDirectory *d = (TDirectory*)rootFile->GetDirectory("");

    TCanvas* c = new TCanvas(hname);

    // Get histograms
    TH1F* h1 = (TH1F*)d->Get(hname);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);

    // Setup the pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->Draw();
        p1->cd();

    // Draw histogram
    THStack* hs = new THStack("hs",hname);
        hs->SetTitle(hname);
        hs->Add(h1,"s");
        hs->Draw("nostack");

    // Make stat box
    p1->Update();
    TPaveStats* st1 = (TPaveStats*)h1->FindObject("stats");
        st1->SetOptStat(110111);
        st1->SetX1NDC(0.8);
        st1->SetX2NDC(0.99);
        st1->SetY1NDC(0.8);
        st1->SetY2NDC(1.0);
        st1->SetTextColor(h1->GetLineColor());

    if(log) {
      p1->SetLogy(); //LOG option
      hs->SetMinimum(0.001); //LOG option
    }
    p1->Update();

    // Save and close
    if(outputname=="") {
      c->Print(outputdir+hname+".png","png");
    } else {
      c->Print(outputdir+outputname+".png","png");
    }
    c->Close();

}

void makeDoubleHist(TString rootFileName1="totalMuFPGAAnalysis.root", TString rootFileName2="totalMuFPGAAnalysis.root", TString hname1="h_itrkCurv", TString hname2="h_itrkCurvWODup", TString outputname="", bool log=false, bool norm=false, bool plotRatio=false) {

    // Load File and prep for loading hists/saving
    TFile *rootFile1 = new TFile(rootFileName1);
    TFile *rootFile2 = new TFile(rootFileName2);
    rootFileName1.ReplaceAll("","");
    rootFileName2.ReplaceAll("","");
    TString outputdir = "Plots/";
    TDirectory *d1;
    TDirectory *d2;
    d1 = (TDirectory*)rootFile1->GetDirectory("");
    d2 = (TDirectory*)rootFile2->GetDirectory("");

    TCanvas* c = new TCanvas("c_comb_"+hname1);

    // Get histograms
    TH1F* h1 = (TH1F*)d1->Get(hname1);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);
        if(norm) h1->Scale(1/h1->Integral());

    TH1F* h2 = (TH1F*)d2->Get(hname2);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);
        if(norm) h2->Scale(1/h2->Integral());

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs",hname1+" & "+hname2);
        hs->SetTitle(hname1+" & "+hname2);

        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Draw("nostack");

    // Make stat boxes
    p1->Update();
    TPaveStats* st1 = (TPaveStats*)h1->FindObject("stats");
        st1->SetOptStat(110111);
        st1->SetX1NDC(0.8);
        st1->SetX2NDC(0.99);
        st1->SetY1NDC(0.8);
        st1->SetY2NDC(1.0);
        st1->SetTextColor(h1->GetLineColor());
    TPaveStats* st2 = (TPaveStats*)h2->FindObject("stats");
        st2->SetOptStat(110111);
        st2->SetX1NDC(0.8);
        st2->SetX2NDC(0.99);
        st2->SetY1NDC(0.6);
        st2->SetY2NDC(0.8);
        st2->SetTextColor(h2->GetLineColor());

    //Log options
    if(log) {
      p1->SetLogy();
      hs->SetMinimum(0.001);
      hs->SetMinimum(1);
    }

    p1->Update();

    if(plotRatio) {
        // Make room for ratio plot
        p1->SetPad(0,0.225,1,1.0);
        p1->Modified();
        p1->Update();

        // Setup the ratio plot
        c->cd();
        TPad *p2 = new TPad("p2","p2",0,0.05,1,0.3);
            p2->SetTopMargin(0);
            p2->SetBottomMargin(0.25);
            p2->SetGrid();
            p2->Draw();
            p2->cd();

        // Plot ratio
        TH1F *h3 = (TH1F*)h1->Clone("h3");
            h3->Reset("ICE");
            h3->Divide(h1,h2);
            h3->Draw();

            // Plot settings
            h3->SetTitle("");
            h3->SetStats(0);
            h3->SetLineColor(kBlue);
            h3->SetMarkerStyle(6);
            h3->SetMarkerColor(kBlue);
            h3->SetMarkerSize(1.0);  

            // X-axis parameters
            h3->SetLabelSize(0.12,"X");
            h3->SetLabelOffset(0.04,"X");
            h3->SetTitleOffset(1.0,"X");
            h3->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
            h3->GetXaxis()->SetTitleColor(kBlack);
            h3->SetTitleFont(62,"X");
            h3->SetTitleSize(0.12,"X");
            // Y-axis parameters
            h3->GetYaxis()->SetNdivisions(50000+404);
            h3->SetMinimum(0);
            h3->GetYaxis()->SetLabelSize(0.1);
            h3->GetYaxis()->SetTitle("Ratio");
            h3->GetYaxis()->SetTitleSize(0.15);
            h3->GetYaxis()->SetTitleOffset(.15);
	  
        p2->Modified();
        p2->Update();        
        
    }

    // Save and close
    if(outputname=="") {
        c->Print(outputdir+"comb_"+hname1+".png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }
//    c->Close();

}

void makeDoubleDoubleHist(TString rootFileName1="totalMuFPGAAnalysis.root", TString rootFileName2="totalMuFPGAAnalysis.root", TString hname1="h_itrkCurv", TString hname2="h_itrkCurvWODup", TString outputname="", bool log=false, bool plotRatio=false) {

    // Load File and prep for loading hists/saving
    TFile *rootFile1 = new TFile(rootFileName1);
    TFile *rootFile2 = new TFile(rootFileName2);
    rootFileName1.ReplaceAll("","");
    rootFileName2.ReplaceAll("","");
    TString outputdir = "Plots/";
    TDirectory *d1;
    TDirectory *d2;
    d1 = (TDirectory*)rootFile1->GetDirectory("");
    d2 = (TDirectory*)rootFile2->GetDirectory("");

    TCanvas* c = new TCanvas("c_comb_"+hname1);

    // Get histograms
    TH1F* h1 = (TH1F*)d1->Get(hname1);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);

    TH1F* h2 = (TH1F*)d2->Get(hname1);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);

    TH1F* h3 = (TH1F*)d1->Get(hname2);
        h3->SetLineColor(kBlack);
        h3->SetLineWidth(3);
        h3->SetLineStyle(2);

    TH1F* h4 = (TH1F*)d2->Get(hname2);
        h4->SetLineColor(kRed);
        h4->SetLineWidth(3);
        h4->SetLineStyle(2);

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs",hname1+" & "+hname2);
        hs->SetTitle(hname1+" & "+hname2);

        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Add(h3,"s");
        hs->Add(h4,"s");
        hs->Draw("nostack");

    // Make stat boxes
    p1->Update();
    TPaveStats* st1 = (TPaveStats*)h1->FindObject("stats");
        st1->SetOptStat(110111);
        st1->SetX1NDC(0.8);
        st1->SetX2NDC(0.99);
        st1->SetY1NDC(0.8);
        st1->SetY2NDC(1.0);
        st1->SetTextColor(h1->GetLineColor());
    TPaveStats* st2 = (TPaveStats*)h2->FindObject("stats");
        st2->SetOptStat(110111);
        st2->SetX1NDC(0.8);
        st2->SetX2NDC(0.99);
        st2->SetY1NDC(0.6);
        st2->SetY2NDC(0.8);
        st2->SetTextColor(h2->GetLineColor());
    h3->SetStats(0);
    h4->SetStats(0);

    //Log options
    if(log) {
      p1->SetLogy();
      hs->SetMinimum(0.001);
      hs->SetMinimum(1);
    }

    p1->Update();

    // Save and close
    if(outputname=="") {
        c->Print(outputdir+"comb_"+hname1+"_"+hname2+".png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }
//    c->Close();

}

void makeProfile(TString rootFileName="totalMuFPGAAnalysis.root", TString hname1="h_NTrkVsmcTrkPt_tot", TString hname2="h_NTrkWODupVsmcTrkPt_tot", TString region="total", TString outputname="0") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("FPGAAnalysis.root","");
    TString outputdir = "Plots/";
    TDirectory *d;
    if(region=="total") {
        d = (TDirectory*)rootFile->GetDirectory("");
    } else {
        d = (TDirectory*)rootFile->GetDirectory(region);
    }

    // Get histograms
    TH2F* h1 = (TH2F*)d->Get(hname1);
    TH2F* h2 = (TH2F*)d->Get(hname2);

    // Hardcode specific plot ranges
    if(rootFileName=="totalMu10_adjsec" && (hname1=="h_NTrkVsmcTrkPt" || hname1=="h_NTrkVsmcTrkPt_tot")) {
        const float xmax=10;
        const float xmin=0;
        h1->SetAxisRange(xmin,xmax,"X");
        h2->SetAxisRange(xmin,xmax,"X");
    }



    // Make profiles
    TProfile* p1 = h1->ProfileX("p_"+hname1);
        p1->SetMarkerSize(0.8);
    TProfile* p2 = h2->ProfileX("p_"+hname2);
        p2->SetMarkerColor(2);
        p2->SetMarkerSize(0.8);
        p2->SetLineColor(2);
    if(rootFileName=="totalMu150_adjsec" && (hname1=="h_NTrkVsmcTrkPt" || hname1=="h_NTrkVsmcTrkPt_tot")) {
        p1->SetBins(150,0,150);
        p2->SetBins(150,0,150);
    }

    // Draw profiles
	TCanvas* c = new TCanvas("c_comb_"+hname1);
        p1->Draw();
        p1->SetTitle("");
        p1->SetStats(0);
        p2->Draw("same");
        c->SetGridx();
        c->SetGridy();

    // Save and close
    if(outputname=="0") {
        c->Print(outputdir+"comb_"+hname1+".png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }
    c->Close();

}

void make2DHist(TString rootFileName="totalMuFPGAAnalysis.root", TString hname="h_PhiTanThDiff", TString region="barrelOnly", TString outputname="0") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("FPGAAnalysis.root","");
    TString outputdir = "Plots/" + rootFileName + "/";
    TDirectory *d;
    if(region=="total") {
        d = (TDirectory*)rootFile->GetDirectory("");
    } else {
        d = (TDirectory*)rootFile->GetDirectory(region);
    }

	TCanvas* c = new TCanvas(hname);
        c->SetGrid();

    // Get histogram
    TH2F* h = (TH2F*)d->Get(hname);

        // hardcode specific plot ranges
        if(rootFileName=="totalMu10" && (hname=="h_NTrkVsmcTrkPt" || hname=="h_NTrkVsmcTrkPt_tot" || hname=="h_NTrkWODupVsmcTrkPt" || hname=="h_NTrkWODupVsmcTrkPt_tot")) {
            const float xmin=0.;
            const float xmax=10.;
            h->SetAxisRange(xmin,xmax,"X");
        }

        // Draw histogram
        h->Draw("colz");
        c->Update();

        // Draw stats box
        TPaveStats* st = (TPaveStats*)h->FindObject("stats");
            st->SetX1NDC(0.71);
            st->SetX2NDC(0.89);
            st->SetY1NDC(0.8);
            st->SetY2NDC(.99);

    // Save and close
    if(outputname=="0") {
        c->Print(outputdir+hname+"_"+region+".png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }

    c->Close();

}

void makeTripleHist(TString rootFileName="totalMuFPGAAnalysis.root", TString hname1="h_nTrkSec_tot", TString hname2="h_nTrkSecWODup_tot", TString hname3="h_nTrkSecMC_tot", TString region="total", TString outputname="0") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("FPGAAnalysis.root","");
    TString outputdir = "Plots/" + rootFileName + "/";
    TDirectory *d;
    if(region=="total") {
        d = (TDirectory*)rootFile->GetDirectory("");
    } else {
        d = (TDirectory*)rootFile->GetDirectory(region);
    }

    TCanvas* c = new TCanvas("c_comb_"+hname1);

    // Get histograms
    TH1F* h1 = (TH1F*)d->Get(hname1);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);

    TH1F* h2 = (TH1F*)d->Get(hname2);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);

    TH1F* h3 = (TH1F*)d->Get(hname3);
        h3->SetLineColor(kYellow+1);
        h3->SetLineWidth(3);

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs",hname1+" & "+hname2+" & "+hname3);
        hs->SetTitle(hname1+" & "+hname2+" & "+hname3);
        hs->SetMinimum(0);
        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Add(h3,"s");
        hs->Draw("nostack");

    // Make stat boxes
    p1->Update();
    TPaveStats* st1 = (TPaveStats*)h1->FindObject("stats");
        st1->SetOptStat(110111);
        st1->SetX1NDC(0.8);
        st1->SetX2NDC(0.99);
        st1->SetY1NDC(0.8);
        st1->SetY2NDC(1.0);
        st1->SetTextColor(h1->GetLineColor());
    TPaveStats* st2 = (TPaveStats*)h2->FindObject("stats");
        st2->SetOptStat(110111);
        st2->SetX1NDC(0.8);
        st2->SetX2NDC(0.99);
        st2->SetY1NDC(0.6);
        st2->SetY2NDC(0.8);
        st2->SetTextColor(h2->GetLineColor());
    TPaveStats* st3 = (TPaveStats*)h3->FindObject("stats");
        st3->SetOptStat(110111);
        st3->SetX1NDC(0.8);
        st3->SetX2NDC(0.99);
        st3->SetY1NDC(0.4);
        st3->SetY2NDC(0.6);
        st3->SetTextColor(h3->GetLineColor());

    p1->Update();

    // Save and close
    if(outputname=="0") {
        c->Print(outputdir+"comb_"+hname1+"_"+region+".png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }
    c->Close();

}

void makeBarrelStubHist(TString rootFileName="mu10Analysis.root", TString ghost="") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("Analysis.root","");
    TString outputdir = "Plots/";
    TDirectory *d;
    d = (TDirectory*)rootFile->GetDirectory("");

    TCanvas* c = new TCanvas("BarrelStubs");

    // Get histograms
    TH1F* h1 = (TH1F*)d->Get("h_L1Multi"+ghost);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);
        h1->SetLineStyle(2);
        h1->Scale(1/h1->GetEntries());
        h1->SetStats(0);

    TH1F* h2 = (TH1F*)d->Get("h_L2Multi"+ghost);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);
        h2->SetLineStyle(2);
        h2->Scale(1/h2->GetEntries());
        h2->SetStats(0);

    TH1F* h3 = (TH1F*)d->Get("h_L3Multi"+ghost);
        h3->SetLineColor(kYellow+1);
        h3->SetLineWidth(3);
        h3->SetLineStyle(2);
        h3->Scale(1/h3->GetEntries());
        h3->SetStats(0);

    TH1F* h4 = (TH1F*)d->Get("h_L4Multi"+ghost);
        h4->SetLineColor(kBlue+1);
        h4->SetLineWidth(3);
        h4->SetLineStyle(2);
        h4->Scale(1/h4->GetEntries());
        h4->SetStats(0);

    TH1F* h5 = (TH1F*)d->Get("h_L5Multi"+ghost);
        h5->SetLineColor(kMagenta+2);
        h5->SetLineWidth(3);
        h5->SetLineStyle(2);
        h5->Scale(1/h5->GetEntries());
        h5->SetStats(0);

    TH1F* h6 = (TH1F*)d->Get("h_L6Multi"+ghost);
        h6->SetLineColor(kOrange+7);
        h6->SetLineWidth(3);
        h6->SetLineStyle(2);
        h6->Scale(1/h6->GetEntries());
        h6->SetStats(0);

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->SetLogy();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs","Barrel Stub Multiplicity");
        hs->SetTitle("Barrel Stub Multiplicity");
        hs->SetMinimum(0.000125);
        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Add(h3,"s");
        hs->Add(h4,"s");
        hs->Add(h5,"s");
        hs->Add(h6,"s");
        hs->Draw("nostack");

    // Make legend
    p1->Update();

    TLegend* leg = new TLegend(0.80,0.55,0.9,0.9);
    //leg->SetMargin(0.035);
    leg->SetTextSize(0.035);
    leg->AddEntry(h1, "L1", "lp");
    leg->AddEntry(h2, "L2", "lp");
    leg->AddEntry(h3, "L3", "lp");
    leg->AddEntry(h4, "L4", "lp");
    leg->AddEntry(h5, "L5", "lp");
    leg->AddEntry(h6, "L6", "lp");
    leg->Draw();

    p1->Update();

    // Save and close
    c->Print(outputdir+"BarrelStubMultiplicity"+ghost+".png","png");
    c->Close();

}

void makeDiskStubHist(TString rootFileName="mu10Analysis.root", TString half="F", TString ghost="") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("Analysis.root","");
    TString outputdir = "Plots/";
    TDirectory *d;
    d = (TDirectory*)rootFile->GetDirectory("");

    TCanvas* c = new TCanvas("BackwardDiskStubs");

    // Get histograms
    TH1F* h1 = (TH1F*)d->Get("h_"+half+"1Multi"+ghost);
    h1->SetLineColorAlpha(kBlack,0.2);
        h1->SetLineWidth(3);
        h1->SetLineStyle(2);
        h1->Scale(1/h1->GetEntries());
        h1->SetStats(0);

    TH1F* h2 = (TH1F*)d->Get("h_"+half+"2Multi"+ghost);
    h2->SetLineColorAlpha(kRed,0.2);
        h2->SetLineWidth(3);
        h2->SetLineStyle(2);
        h2->Scale(1/h2->GetEntries());
        h2->SetStats(0);

    TH1F* h3 = (TH1F*)d->Get("h_"+half+"3Multi"+ghost);
    h3->SetLineColorAlpha(kYellow+1,0.2);
        h3->SetLineWidth(3);
        h3->SetLineStyle(2);
        h3->Scale(1/h3->GetEntries());
        h3->SetStats(0);

    TH1F* h4 = (TH1F*)d->Get("h_"+half+"4Multi"+ghost);
    h4->SetLineColorAlpha(kBlue+1,0.2);
        h4->SetLineWidth(3);
        h4->SetLineStyle(2);
        h4->Scale(1/h4->GetEntries());
        h4->SetStats(0);

    TH1F* h5 = (TH1F*)d->Get("h_"+half+"5Multi"+ghost);
    h5->SetLineColorAlpha(kMagenta+2,0.2);
        h5->SetLineWidth(3);
        h5->SetLineStyle(2);
        h5->Scale(1/h5->GetEntries());
        h5->SetStats(0);

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->SetLogy();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs","Disk Stub Multiplicity");
        hs->SetTitle(half+" Disk Stub Multiplicity");
        hs->SetMinimum(0.000125);
        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Add(h3,"s");
        hs->Add(h4,"s");
        hs->Add(h5,"s");
        hs->Draw("nostack");

    // Make legend
    p1->Update();

    TLegend* leg = new TLegend(0.80,0.55,0.9,0.9);
    //leg->SetMargin(0.035);
    leg->SetTextSize(0.035);
    leg->AddEntry(h1, half+"1", "lp");
    leg->AddEntry(h2, half+"2", "lp");
    leg->AddEntry(h3, half+"3", "lp");
    leg->AddEntry(h4, half+"4", "lp");
    leg->AddEntry(h5, half+"5", "lp");
    leg->Draw();

    p1->Update();

    // Save and close
    c->Print(outputdir+half+"DiskStubMultiplicity"+ghost+".png","png");
    c->Close();

}

void makeBarrelTrackHist(TString rootFileName="mu10Analysis.root", TString ghost="") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("Analysis.root","");
    TString outputdir = "Plots/";
    TDirectory *d;
    d = (TDirectory*)rootFile->GetDirectory("");

    TCanvas* c = new TCanvas("BarrelTracks");

    // Get histograms
    TH1F* h1 = (TH1F*)d->Get("h_L1L2Multi"+ghost);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);
        h1->SetLineStyle(2);
        h1->Scale(1/h1->GetEntries());
        h1->SetStats(0);

    TH1F* h2 = (TH1F*)d->Get("h_L3L4Multi"+ghost);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);
        h2->SetLineStyle(2);
        h2->Scale(1/h2->GetEntries());
        h2->SetStats(0);

    TH1F* h3 = (TH1F*)d->Get("h_L5L6Multi"+ghost);
        h3->SetLineColor(kYellow+1);
        h3->SetLineWidth(3);
        h3->SetLineStyle(2);
        h3->Scale(1/h3->GetEntries());
        h3->SetStats(0);

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->SetLogy();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs","Barrel Track Multiplicity");
        hs->SetTitle("Barrel Track Multiplicity");
        hs->SetMinimum(0.000125);
        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Add(h3,"s");
        hs->Draw("nostack");

    // Make legend
    p1->Update();

    TLegend* leg = new TLegend(0.80,0.55,0.9,0.9);
    //leg->SetMargin(0.035);
    leg->SetTextSize(0.035);
    leg->AddEntry(h1, "L1L2", "lp");
    leg->AddEntry(h2, "L3L4", "lp");
    leg->AddEntry(h3, "L5L6", "lp");
    leg->Draw();

    p1->Update();

    // Save and close
    c->Print(outputdir+"BarrelTrackMultiplicity"+ghost+".png","png");
    c->Close();

}

void makeDiskTrackHist(TString rootFileName="mu10Analysis.root", TString half="F", TString ghost="") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("Analysis.root","");
    TString outputdir = "Plots/";
    TDirectory *d;
    d = (TDirectory*)rootFile->GetDirectory("");

    TCanvas* c = new TCanvas("DiskTracks");

    // Get histograms
    TH1F* h1 = (TH1F*)d->Get("h_"+half+"1"+half+"2"+"Multi"+ghost);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);
        h1->SetLineStyle(2);
        h1->Scale(1/h1->GetEntries());
        h1->SetStats(0);

    TH1F* h2 = (TH1F*)d->Get("h_"+half+"3"+half+"4"+"Multi"+ghost);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);
        h2->SetLineStyle(2);
        h2->Scale(1/h2->GetEntries());
        h2->SetStats(0);

    TH1F* h3 = (TH1F*)d->Get("h_"+half+"1L"+"Multi"+ghost);
        h3->SetLineColor(kYellow+1);
        h3->SetLineWidth(3);
        h3->SetLineStyle(2);
        h3->Scale(1/h3->GetEntries());
        h3->SetStats(0);

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->SetLogy();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs","Disk Track Multiplicity");
        hs->SetTitle("Disk Track Multiplicity");
        hs->SetMinimum(0.000125);
        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Add(h3,"s");
        hs->Draw("nostack");

    // Make legend
    p1->Update();

    TLegend* leg = new TLegend(0.80,0.55,0.9,0.9);
    //leg->SetMargin(0.035);
    leg->SetTextSize(0.035);
    leg->AddEntry(h1, half+"1"+half+"2", "lp");
    leg->AddEntry(h2, half+"3"+half+"4", "lp");
    leg->AddEntry(h3, half+"1L", "lp");
    leg->Draw();

    p1->Update();

    // Save and close
    c->Print(outputdir+half+"DiskTrackMultiplicity"+ghost+".png","png");
    c->Close();

}


void resolutionPlot(TString rootFileName="totalMuFPGAAnalysis.root", TString h1name="h_matchDeltaPhi_AllM", TString outputname="") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    TString outputdir = "Plots/Resolution/";
    TDirectory *d;
    d = (TDirectory*)rootFile->GetDirectory("");

	TCanvas* c = new TCanvas(h1name);
        c->SetGrid();

    // Get histogram and draw the fit curve on top
    TH1F* h1 = (TH1F*)d->Get(h1name);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);
        h1->Draw();
        h1->Fit("gaus","Q");

    TF1* g1 = h1->GetFunction("gaus");
        g1->SetLineColor(kRed);

    // Save and close
    if(outputname=="") {
        c->Print(outputdir+h1name+"_res.png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }
    c->Close();
}

void doubleResPlots(TString rootFileName1="plusMu10FPGAAnalysis.root", TString rootFileName2="minusMu10FPGAAnalysis.root", TString variable="Pt", TString region="total", TString outputname="0") {

    // Load File and prep for loading hists/saving
    TFile *rootFile1 = new TFile(rootFileName1);
    TFile *rootFile2 = new TFile(rootFileName2);
    TString outputdir = "Plots/";
    TDirectory *d1;
    TDirectory *d2;
    TString hname;
    if(region=="total") {
        d1 = (TDirectory*)rootFile1->GetDirectory("");
        d2 = (TDirectory*)rootFile2->GetDirectory("");
        hname = "h_matchDelta"+variable+"_AllM_tot";
    } else {
        d1 = (TDirectory*)rootFile1->GetDirectory(region);
        d2 = (TDirectory*)rootFile2->GetDirectory(region);
        hname = "h_matchDelta"+variable+"_AllM";
    }

	TCanvas* c = new TCanvas(hname);
        c->SetGrid();

    // Get histograms and draw the fit curves on top
    TH1F* h1 = (TH1F*)d1->Get(hname);
        h1->SetLineColor(kRed);
        h1->SetLineWidth(3);
        h1->Fit("gaus","Q");

    TH1F* h2 = (TH1F*)d2->Get(hname);
        h2->SetLineColor(kBlue);
        h2->SetLineWidth(3);
        h2->Fit("gaus","Q");

    TH1F* h3 = new TH1F(*h1);
        h3->Add(h2);
        h3->SetLineColor(kBlack);
        h3->SetLineWidth(3);
        h3->Fit("gaus","Q");

    TF1* g1 = h1->GetFunction("gaus");
        g1->SetLineColor(46);
    TF1* g2 = h2->GetFunction("gaus");
        g2->SetLineColor(38);
    TF1* g3 = h3->GetFunction("gaus");
        g3->SetLineColor(13);

    c->Clear();

    h3->Draw();
    h1->Draw("sames");
    h2->Draw("sames");

    c->Update();

    TPaveStats* st1 = (TPaveStats*)h1->FindObject("stats");
        st1->SetX1NDC(0.8);
        st1->SetX2NDC(0.995);
        st1->SetY1NDC(0.8);
        st1->SetY2NDC(1.0);
        st1->SetTextColor(h1->GetLineColor());
    TPaveStats* st2 = (TPaveStats*)h2->FindObject("stats");
        st2->SetX1NDC(0.8);
        st2->SetX2NDC(0.995);
        st2->SetY1NDC(0.6);
        st2->SetY2NDC(0.8);
        st2->SetTextColor(h2->GetLineColor());
    TPaveStats* st3 = (TPaveStats*)h3->FindObject("stats");
        st3->SetX1NDC(0.8);
        st3->SetX2NDC(0.995);
        st3->SetY1NDC(0.4);
        st3->SetY2NDC(0.6);
        st3->SetTextColor(h3->GetLineColor());


    // Save and close
    c->Modified();
    c->Update();
    if(outputname=="0") {
        c->Print(outputdir+hname+"_"+region+".png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }
    c->Close();
}

void combinationResolution(TString rootFileName="totalMuFPGAAnalysis.root", TString variable="Phi", TString regionType="tracklet", TString outputname="plotName") {

    // Load File and prep for loading hists/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("Analysis.root","");
    TString outputdir = "Plots/" + rootFileName + "/Resolution/";
    TDirectory *d;

    // Create objects to fill
    TCanvas* c = new TCanvas();
    c->SetGrid();

    TMultiGraph* mg = new TMultiGraph;
    mg->Draw("apl");

    TGraphErrors* g[4];

    TLegend* leg = new TLegend(0.525,0.75,0.9,0.9);
    leg->SetMargin(0.035);
    leg->SetTextSize(0.035);
    leg->Draw();

    TH1F* htot[3];

    // Get the resolution for each region and plot, separated by nStub
    for(int r=0; r<3; r++) {
        TString region;
        if(regionType=="tracklet") {
            if(r==0) region = "barrelOnly";
            if(r==1) region = "diskOnly";
            if(r==2) region = "overlap";
        }
        if(regionType=="eta") {
            if(r==0) region = "mcEtaLess1";
            if(r==1) region = "mcEta1to1.7";
            if(r==2) region = "mcEtaMore1.7";
        }
        d = (TDirectory*)rootFile->GetDirectory(region);
        g[r] = new TGraphErrors();
        g[r]->SetMarkerStyle(8);
        g[r]->SetMarkerColor(r+2);
        g[r]->SetLineColor(r+2);
        g[r]->SetMarkerSize(1.25);

        for(int nStub=4; nStub<=6; nStub++) {
            TString hname = Form("h_matchDelta"+variable+"_nS%i",nStub);
            TH1F* h = (TH1F*)d->Get(hname);
            if(h->GetEntries()==0) continue;

            // Sum over regions to make total hists
            if(r==0) htot[nStub-4]=h;
            else htot[nStub-4]->Add(h,1);

            // Get resolution points except for 4 stub overlap (bad statistics)
            if(regionType!="tracklet" || nStub!=4 || r!=2) {
              cout << "test1" << endl;
                h->Fit("gaus","Q0");
                cout << "test2" << endl;
                TF1* gaus = h->GetFunction("gaus");
                g[r]->SetPoint(g[r]->GetN(),nStub,gaus->GetParameter(2));
                g[r]->SetPointError(g[r]->GetN()-1,0,gaus->GetParError(2));
            }
        }
        // Plot all regional resolution graphs
        mg->Add(g[r]);
        leg->AddEntry(g[r],variable+" resolution, "+region+" tracks","p");
    }

    // Plot total hist resolutions
    g[3] = new TGraphErrors();
    for(int nStub=4; nStub<=6; nStub++) {
        htot[nStub-4]->Fit("gaus","Q0");
        TF1* gaus = htot[nStub-4]->GetFunction("gaus");
        g[3]->SetPoint(g[3]->GetN(),nStub,gaus->GetParameter(2));
        g[3]->SetPointError(g[3]->GetN()-1,0,gaus->GetParError(2));
    }
    g[3]->SetMarkerStyle(8);
    g[3]->SetMarkerColor(1);
    g[3]->SetLineColor(1);
    g[3]->SetMarkerSize(1.25);
    mg->Add(g[3]);
    leg->AddEntry(g[3],variable+" resolution, all tracks","p");        

    mg->GetXaxis()->SetTitle("Number of Stubs");

    // Save and close
    c->Update();
    c->Print(outputdir+outputname+".png","png");
    c->Close();

}

void efficiencyPlot(TString rootFileName="totalMuFPGAAnalysis.root", TString variable="Phi", TString region="etaLess1", TString outputname="plotName") {

    // Load File and prep for loading hist/saving
    TFile *rootFile = new TFile(rootFileName);
    rootFileName.ReplaceAll("FPGAAnalysis.root","");
    TString outputdir = "Plots/" + rootFileName + "/";
	TCanvas* c = new TCanvas(outputname);
        c->SetGrid();

    // If total region, get total efficiency plot by summing 3 regions
    TDirectory* d;
    TEfficiency* e;
    if(region=="total") {
        TH1F* h1tot = new TH1F;
        TH1F* h2tot = new TH1F;

        for(int r=0; r<3; r++) {
            if(r==0) region = "mcEtaLess1";
            if(r==1) region = "mcEta1to1.7";
            if(r==2) region = "mcEtaMore1.7";
    
            d = (TDirectory*)rootFile->GetDirectory(region);
            TH1F* h1 = (TH1F*)d->Get("h_mcTrk"+variable+"_Fnd");
            TH1F* h2 = (TH1F*)d->Get("h_mcTrk"+variable);

            // Initialize h1tot, ht2tot the first time, then add histograms on top
            if(r==0) {
                h1tot=h1;
                h2tot=h2;
            } else {
                h1tot->Add(h1);
                h2tot->Add(h2);
            }
        }

        e = new TEfficiency(*h1tot,*h2tot);

    // Otherwise, make regional efficiency plot
    } else {

        d = (TDirectory*)rootFile->GetDirectory(region);
        TH1F* h1 = (TH1F*)d->Get("h_mcTrk"+variable+"_Fnd");
        TH1F* h2 = (TH1F*)d->Get("h_mcTrk"+variable);

        e = new TEfficiency(*h1,*h2);

    }
    e->SetTitle(";"+variable+";"); // x-axis title cut off, so this is currently useless
    e->Draw();

    // Save and close
    c->Update();
    c->Print(outputdir+outputname+".png","png");
    c->Close();

}

void makeTripleHistSpecial(TString rootFileName1="totalMuFPGAAnalysis.root", TString rootFileName2="totalMuFPGAAnalysis.root", TString rootFileName3="totalMuFPGAAnalysis.root", TString hname1="h_nTrkSec_tot", TString hname2="h_nTrkSecWODup_tot", TString hname3="h_nTrkSecWODup_tot", TString region="total", TString outputname="0") {

    // Load File and prep for loading hists/saving
    TFile *rootFile1 = new TFile(rootFileName1);
    TFile *rootFile2 = new TFile(rootFileName2);
    TFile *rootFile3 = new TFile(rootFileName3);
    rootFileName1.ReplaceAll("FPGAAnalysis.root","");
    rootFileName2.ReplaceAll("FPGAAnalysis.root","");
    rootFileName3.ReplaceAll("FPGAAnalysis.root","");
    TString outputdir = "Plots/";
    TDirectory *d1;
    TDirectory *d2;
    TDirectory *d3;
    if(region=="total") {
        d1 = (TDirectory*)rootFile1->GetDirectory("");
        d2 = (TDirectory*)rootFile2->GetDirectory("");
        d3 = (TDirectory*)rootFile3->GetDirectory("");
    } else {
        d1 = (TDirectory*)rootFile1->GetDirectory(region);
        d2 = (TDirectory*)rootFile2->GetDirectory(region);
        d3 = (TDirectory*)rootFile3->GetDirectory(region);
    }

    TCanvas* c = new TCanvas("c_"+hname1);

    // Get histograms
    TH1F* h1 = (TH1F*)d1->Get(hname1);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);
        h1->SetName("Before Removal");
//        h1->GetXaxis()->SetRangeUser(1,50);
//        h1->SetName("No Ordering");

    TH1F* h2 = (TH1F*)d2->Get(hname2);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);
//        h2->GetXaxis()->SetRangeUser(1,50);
        h2->SetName("Same Sector Removal");
//        h2->SetName("Descending Order");

    TH1F* h3 = (TH1F*)d3->Get(hname3);
        h3->SetLineColor(kBlue);
        h3->SetLineWidth(3);
//        h3->GetXaxis()->SetRangeUser(1,50);
        h3->SetName("Adjacent Sector Removal");
//        h3->SetName("Ascending Order");

    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->Draw();
        p1->cd();

    // Draw histograms
    THStack* hs = new THStack("hs",hname1+" & "+hname2+" & "+hname3);
//        hs->SetTitle(hname1+" & "+hname2+" & "+hname3);
        hs->SetTitle("Number of Tracks per Event");
        hs->SetMinimum(0.00001);
        p1->SetLogy();
        h1->Scale(1/h1->GetEntries());
        h2->Scale(1/h2->GetEntries());
        h3->Scale(1/h3->GetEntries());
        hs->Add(h1,"s");
        hs->Add(h2,"s");
        hs->Add(h3,"s");
        hs->Draw("nostack");
//        hs->GetXaxis()->SetRangeUser(1,20);
//        hs->SetMaximum(40000);


    // Make stat boxes
    p1->Update();
    TPaveStats* st1 = (TPaveStats*)h1->FindObject("stats");
        st1->SetX1NDC(0.8);
        st1->SetX2NDC(0.99);
        st1->SetY1NDC(0.8);
        st1->SetY2NDC(1.0);
        st1->SetTextColor(h1->GetLineColor());
    TPaveStats* st2 = (TPaveStats*)h2->FindObject("stats");
        st2->SetX1NDC(0.8);
        st2->SetX2NDC(0.99);
        st2->SetY1NDC(0.6);
        st2->SetY2NDC(0.8);
        st2->SetTextColor(h2->GetLineColor());
        st2->AddText("SameSec");
    TPaveStats* st3 = (TPaveStats*)h3->FindObject("stats");
        st3->SetX1NDC(0.8);
        st3->SetX2NDC(0.99);
        st3->SetY1NDC(0.4);
        st3->SetY2NDC(0.6);
        st3->SetTextColor(h3->GetLineColor());
        st3->AddText("AdjSec");

    p1->Update();

    // Save and close
    if(outputname=="0") {
        c->Print(outputdir+"comb_"+hname1+"_"+region+".png","png");
    } else {
        c->Print(outputdir+outputname+".png","png");
    }
    c->Close();

}

void makeTDRPlot(TString rootFileName1="mu_40000_2to10_adj_analyze.root", TString rootFileName2="mu_40000_2to10_ss_analyze.root", TString rootFileName3="mu_40000_2to10_adj_analyze.root", TString hname1="h_nTrkEvt_tot", TString hname2="h_nTrkEvtWODup_tot", TString hname3="h_nTrkEvtWODup_tot", TString outputname="DuplicateRemovalComparison2-10GeVMuons") {

    // Load File and prep for loading hists/saving
    TFile *rootFile1 = new TFile(rootFileName1);
    TFile *rootFile2 = new TFile(rootFileName2);
    TFile *rootFile3 = new TFile(rootFileName3);
    TString outputdir = "Plots/";
    TDirectory *d1;
    TDirectory *d2;
    TDirectory *d3;
    d1 = (TDirectory*)rootFile1->GetDirectory("");
    d2 = (TDirectory*)rootFile2->GetDirectory("");
    d3 = (TDirectory*)rootFile3->GetDirectory("");

    TCanvas* c = new TCanvas("c_"+hname1);

    // Get histograms
    TH1F* h1 = (TH1F*)d1->Get(hname1);
        h1->SetLineColor(kBlack);
        h1->SetLineWidth(3);
        h1->Scale(1/(h1->Integral()));

    TH1F* h2 = (TH1F*)d2->Get(hname2);
        h2->SetLineColor(kRed);
        h2->SetLineWidth(3);
        h2->Scale(1/(h2->Integral()));

    TH1F* h3 = (TH1F*)d3->Get(hname3);
        h3->SetLineColor(kBlue);
        h3->SetLineWidth(3);
        h3->Scale(1/(h3->Integral()));


    // Setup the superimposed histograms pad
    TPad *p1 = new TPad("p1","p1",0,0.0,1,1.0);
        p1->SetGrid();
        p1->Draw();
        p1->cd();

    // Draw histograms  
    THStack* hs = new THStack("hs","");
        hs->SetMinimum(0);
        hs->Add(h1,"");
        hs->Add(h2,"");
        hs->Add(h3,"");
        hs->Draw("nostack");
        p1->SetPad(0,0.025,1,1);
        hs->GetXaxis()->Set(10,0.5,10.5);
        c->Update();


    // Change Style
    TLegend* leg = new TLegend(0.475,0.7,0.9,0.9,"Single muons 2<P_{T}<10GeV");
    leg->SetMargin(0.13);
    leg->SetTextSize(0.035);
    leg->SetFillColor(0);
    leg->Draw();
    leg->AddEntry(h1,"Before duplicate removal","l");
    leg->AddEntry(h2,"Duplicate removal within sectors","l");
    leg->AddEntry(h3,"+ Duplicate removal between sectors","l");
    hs->SetMinimum(0.001); //LOG SETTINGS
    p1->SetLogy(); // LOG SETTINGS
    p1->Update();

    // Add Labels
//    TLatex* title = new TLatex(2.5,1.05,"CMS Preliminary Simulation, Phase-2"); //LINEAR SETTINGS
    TLatex* title = new TLatex(3.0,1.9,"CMS Preliminary Simulation, Phase-2"); //LOG SETTINGS
    title->SetTextSize(0.035);
    title->Draw();
//    TLatex* xaxis = new TLatex(2.5,-0.1,"Number of tracks (#geq1) per event"); //LINEAR SETTINGS
    TLatex* xaxis = new TLatex(2.5,0.0005,"Number of tracks (#geq1) per event"); //LOG SETTINGS
    xaxis->Draw();
    TLatex* yaxis = new TLatex(-0.12,0.03,"Fraction of events"); //LOG SETTINGS
    yaxis->SetTextAngle(90);
    yaxis->Draw();
    p1->Update();

    // Save and close
    c->Print(outputdir+outputname+".pdf","pdf");
//    c->Close();

}

void allPlots(TString rootFile="totalMuFPGAAnalysis.root") {

//    makeDoubleHist(rootFile, rootFile, "h_NstubShare", "h_NstubShare_Near", "barrelOnly", "stubShare_barrel", true);
//    makeDoubleHist(rootFile, rootFile, "h_NstubShare", "h_NstubShare_Near", "diskOnly", "stubShare_disk", true);
//    makeDoubleHist(rootFile, rootFile, "h_NstubShare", "h_NstubShare_Near", "overlap", "stubShare_overlap", true);
//    makeDoubleHist(rootFile, rootFile, "h_nTrkSec", "h_nTrkSecWODup", "barrelOnly", "nTrkSec_barrel", false);
//    makeDoubleHist(rootFile, rootFile, "h_nTrkSec", "h_nTrkSecWODup", "diskOnly", "nTrkSec_disk", false);
//    makeDoubleHist(rootFile, rootFile, "h_nTrkSec", "h_nTrkSecWODup", "overlap", "nTrkSec_overlap", false);
//    makeDoubleHist(rootFile, rootFile, "h_nTrkSec_tot", "h_nTrkSecWODup_tot", "total", "nTrkSec_total", false);
//    makeDoubleHist(rootFile, rootFile, "h_itrkCurv", "h_itrkCurvWODup", "barrelOnly", "iCurv_barrel", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkCurv", "h_itrkCurvWODup", "diskOnly", "iCurv_disk", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkCurv", "h_itrkCurvWODup", "overlap", "iCurv_overlap", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkPt", "h_trkPtWODup", "barrelOnly", "Pt_barrel", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkPt", "h_trkPtWODup", "diskOnly", "Pt_disk", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkPt", "h_trkPtWODup", "overlap", "Pt_overlap", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkPhi", "h_itrkPhiWODup", "barrelOnly", "iPhi_barrel", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkPhi", "h_itrkPhiWODup", "diskOnly", "iPhi_disk", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkPhi", "h_itrkPhiWODup", "overlap", "iPhi_overlap", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkTanTh", "h_itrkTanThWODup", "barrelOnly", "TanTh_barrel", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkTanTh", "h_itrkTanThWODup", "diskOnly", "TanTh_disk", true);
//    makeDoubleHist(rootFile, rootFile, "h_itrkTanTh", "h_itrkTanThWODup", "overlap", "TanTh_overlap", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkEta", "h_trkEtaWODup", "barrelOnly", "Eta_barrel", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkEta", "h_trkEtaWODup", "diskOnly", "Eta_disk", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkEta", "h_trkEtaWODup", "overlap", "Eta_overlap", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkPhi", "h_trkPhiWODup", "barrelOnly", "Phi_barrel", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkPhi", "h_trkPhiWODup", "diskOnly", "Phi_disk", true);
//    makeDoubleHist(rootFile, rootFile, "h_trkPhi", "h_trkPhiWODup", "overlap", "Phi_overlap", true);

    makeTripleHist(rootFile, "h_nTrkSec_tot", "h_nTrkSecWODup_tot", "h_nMCTrkSec_tot", "total", "nTrkSecComb_tot");
    makeTripleHist(rootFile, "h_nTrkSec", "h_nTrkSecWODup", "h_nMCTrkSec", "mcEtaLess1", "nTrkSecComb_etaLess1");
    makeTripleHist(rootFile, "h_nTrkSec", "h_nTrkSecWODup", "h_nMCTrkSec", "mcEta1to1.7", "nTrkSecComb_eta1to1.7");
    makeTripleHist(rootFile, "h_nTrkSec", "h_nTrkSecWODup", "h_nMCTrkSec", "mcEtaMore1.7", "nTrkSecComb_etaMore1.7");

    make2DHist(rootFile, "h_PhiTanThDiff", "barrelOnly", "PhiTanThDiff_barrel");
    make2DHist(rootFile, "h_PhiTanThDiff", "diskOnly", "PhiTanThDiff_disk");
    make2DHist(rootFile, "h_PhiTanThDiff", "overlap", "PhiTanThDiff_overlap");
    make2DHist(rootFile, "h_NTrkVsmcTrkPt", "mcEtaLess1", "NTrkvsmcPt_EtaLess1");
    make2DHist(rootFile, "h_NTrkVsmcTrkPt", "mcEta1to1.7", "NTrkvsmcPt_Eta1to1.7");
    make2DHist(rootFile, "h_NTrkVsmcTrkPt", "mcEtaMore1.7", "NTrkvsmcPt_EtaMore1.7");
    make2DHist(rootFile, "h_NTrkVsmcTrkPt_tot", "total", "NTrkvsmcPt_total");
    make2DHist(rootFile, "h_NTrkWODupVsmcTrkPt", "mcEtaLess1", "NTrkWODupvsmcPt_EtaLess1");
    make2DHist(rootFile, "h_NTrkWODupVsmcTrkPt", "mcEta1to1.7", "NTrkWODupvsmcPt_Eta1to1.7");
    make2DHist(rootFile, "h_NTrkWODupVsmcTrkPt", "mcEtaMore1.7", "NTrkWODupvsmcPt_EtaMore1.7");
    make2DHist(rootFile, "h_NTrkWODupVsmcTrkPt_tot", "total", "NTrkWODupvsmcPt_total");
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


void DecemberPlots(TString rootFile="ttbar140_barrel_ana.root") {
  makeSingleHist(rootFile,"h_nTrkEvt_tot");
  makeSingleHist(rootFile,"h_nDupTrk_tot");
  makeSingleHist(rootFile,"h_nTrkEvtWODup_tot");
  makeSingleHist(rootFile,"h_nMCTrkEvt_tot");
  makeSingleHist(rootFile,"h_trkPtWODup_tot");
  makeSingleHist(rootFile,"h_trkPhiWODup_tot");
  makeSingleHist(rootFile,"h_trkEtaWODup_tot");

}

void DoubleDecemberPlots(TString rootFile1="ttbar140_barrel_fixedIDs_ana.root", TString rootFile2="ttbar140_barrel_ana.root") {
  makeDoubleDoubleHist(rootFile1,rootFile2,"h_nTrkEvt_tot","h_nTrkEvtWODup_tot");
  makeDoubleHist(rootFile1,rootFile2,"h_nTrkEvt_tot","h_nTrkEvt_tot");
  makeDoubleHist(rootFile1,rootFile2,"h_nTrkEvtWODup_tot","h_nTrkEvtWODup_tot");
  makeDoubleHist(rootFile1,rootFile2,"h_nDupTrk_tot","h_nDupTrk_tot");
  makeDoubleHist(rootFile1,rootFile2,"h_nMCTrkEvt_tot","h_nMCTrkEvt_tot");
  makeDoubleHist(rootFile1,rootFile2,"h_trkPtWODup_tot","h_trkPtWODup_tot","",false,true);
  makeDoubleHist(rootFile1,rootFile2,"h_trkPhiWODup_tot","h_trkPhiWODup_tot","",false,true);
  makeDoubleHist(rootFile1,rootFile2,"h_trkEtaWODup_tot","h_trkEtaWODup_tot","",false,true);

}
