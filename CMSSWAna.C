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

TTree *cmsswTree;
int numCMSSWEvts;


void loadCMSSWTree(TString fileName = "myTest.root") {
 
  cout << "Loading root tree " << fileName << endl;
  TFile *tFile = new TFile(fileName.Data());
//  tFile->cd("L1TrackNtuple");
  cmsswTree = (TTree*)tFile->Get("L1TrackNtuple/eventTree");

//  cmssw->SetBranchAddress("eventTree", &Event);
  numCMSSWEvts = (int)cmsswTree->GetEntries();
 
  printf("Number of Entries %i\n",numCMSSWEvts);

}

void analyzeCMSSWEvent(TString histFileName="test.root", TString rootFileName="myCMSSWTest.root") {

  const double Maxpt = 200.0;
  const double Minpt = 0.0;
  const int minStubs = 4;

  loadCMSSWTree(rootFileName);

  //Setup histogram locations
  TFile *histFile = new TFile(histFileName,"RECREATE");
  histFile->cd();

  // define leafs & branches

  // tracking particles
  vector<float>* tp_pt;
  vector<float>* tp_eta;
  vector<float>* tp_phi;
  vector<float>* tp_dxy;
  vector<float>* tp_z0;
  vector<float>* tp_d0;
  vector<int>*   tp_pdgid;
  vector<int>*   tp_nmatch;
  vector<int>*   tp_nstub;
  vector<int>*   tp_eventid;
  
  // *L1 track* properties, for tracking particles matched to a L1 track
  vector<float>* matchtrk_pt;
  vector<float>* matchtrk_eta;
  vector<float>* matchtrk_phi;
  vector<float>* matchtrk_d0;
  vector<float>* matchtrk_z0;
  vector<float>* matchtrk_chi2;
  vector<int>*   matchtrk_nstub;

  // *L1 track* properties, for *all* tracking particles matched to a L1 track
  vector<vector<float>>* allmatchtrk_pt;
  vector<vector<float>>* allmatchtrk_eta;
  vector<vector<float>>* allmatchtrk_phi;
  vector<vector<float>>* allmatchtrk_z0;
  vector<vector<float>>* allmatchtrk_chi2;
  vector<vector<int>>*   allmatchtrk_nstub;

  // all L1 tracks
  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<float>* trk_z0;
  vector<float>* trk_chi2;
  vector<int>*   trk_nstub;

  TBranch* b_tp_pt;
  TBranch* b_tp_eta;
  TBranch* b_tp_phi;
  TBranch* b_tp_dxy;
  TBranch* b_tp_z0;
  TBranch* b_tp_d0;
  TBranch* b_tp_pdgid;
  TBranch* b_tp_nmatch;
  TBranch* b_tp_nstub;
  TBranch* b_tp_eventid;

  TBranch* b_matchtrk_pt;
  TBranch* b_matchtrk_eta;
  TBranch* b_matchtrk_phi;
  TBranch* b_matchtrk_d0;
  TBranch* b_matchtrk_z0;
  TBranch* b_matchtrk_chi2; 
  TBranch* b_matchtrk_nstub;

  TBranch* b_allmatchtrk_pt;
  TBranch* b_allmatchtrk_eta;
  TBranch* b_allmatchtrk_phi;
  TBranch* b_allmatchtrk_z0;
  TBranch* b_allmatchtrk_chi2; 
  TBranch* b_allmatchtrk_nstub;

  TBranch* b_trk_pt; 
  TBranch* b_trk_eta; 
  TBranch* b_trk_phi;
  TBranch* b_trk_z0; 
  TBranch* b_trk_chi2; 
  TBranch* b_trk_nstub; 

  tp_pt  = 0;
  tp_eta = 0;
  tp_phi = 0;
  tp_dxy = 0;
  tp_z0  = 0;
  tp_d0  = 0;
  tp_pdgid = 0;
  tp_nmatch = 0;
  tp_nstub = 0;
  tp_eventid = 0;

  matchtrk_pt  = 0;
  matchtrk_eta = 0;
  matchtrk_phi = 0;
  matchtrk_d0  = 0;
  matchtrk_z0  = 0;
  matchtrk_chi2  = 0; 
  matchtrk_nstub = 0;

  allmatchtrk_pt  = 0;
  allmatchtrk_eta = 0;
  allmatchtrk_phi = 0;
  allmatchtrk_z0  = 0;
  allmatchtrk_chi2  = 0; 
  allmatchtrk_nstub = 0;

  trk_pt = 0; 
  trk_eta = 0; 
  trk_phi = 0; 
  trk_z0 = 0; 
  trk_chi2 = 0; 
  trk_nstub = 0; 

  cmsswTree->SetBranchAddress("tp_pt",     &tp_pt,     &b_tp_pt);
  cmsswTree->SetBranchAddress("tp_eta",    &tp_eta,    &b_tp_eta);
  cmsswTree->SetBranchAddress("tp_phi",    &tp_phi,    &b_tp_phi);
  cmsswTree->SetBranchAddress("tp_dxy",    &tp_dxy,    &b_tp_dxy);
  cmsswTree->SetBranchAddress("tp_z0",     &tp_z0,     &b_tp_z0);
  cmsswTree->SetBranchAddress("tp_d0",     &tp_d0,     &b_tp_d0);
  cmsswTree->SetBranchAddress("tp_pdgid",  &tp_pdgid,  &b_tp_pdgid);
  cmsswTree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
  cmsswTree->SetBranchAddress("tp_nstub",      &tp_nstub,      &b_tp_nstub);
  cmsswTree->SetBranchAddress("tp_eventid",    &tp_eventid,    &b_tp_eventid);

  cmsswTree->SetBranchAddress("matchtrk_pt",    &matchtrk_pt,    &b_matchtrk_pt);
  cmsswTree->SetBranchAddress("matchtrk_eta",   &matchtrk_eta,   &b_matchtrk_eta);
  cmsswTree->SetBranchAddress("matchtrk_phi",   &matchtrk_phi,   &b_matchtrk_phi);
  cmsswTree->SetBranchAddress("matchtrk_d0",    &matchtrk_d0,    &b_matchtrk_d0);
  cmsswTree->SetBranchAddress("matchtrk_z0",    &matchtrk_z0,    &b_matchtrk_z0);
  cmsswTree->SetBranchAddress("matchtrk_chi2",  &matchtrk_chi2,  &b_matchtrk_chi2);
  cmsswTree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);

  cmsswTree->SetBranchAddress("allmatchtrk_pt",    &allmatchtrk_pt,    &b_allmatchtrk_pt);
  cmsswTree->SetBranchAddress("allmatchtrk_eta",   &allmatchtrk_eta,   &b_allmatchtrk_eta);
  cmsswTree->SetBranchAddress("allmatchtrk_phi",   &allmatchtrk_phi,   &b_allmatchtrk_phi);
  cmsswTree->SetBranchAddress("allmatchtrk_z0",    &allmatchtrk_z0,    &b_allmatchtrk_z0);
  cmsswTree->SetBranchAddress("allmatchtrk_chi2",  &allmatchtrk_chi2,  &b_allmatchtrk_chi2);
  cmsswTree->SetBranchAddress("allmatchtrk_nstub", &allmatchtrk_nstub, &b_allmatchtrk_nstub);

  cmsswTree->SetBranchAddress("trk_pt",   &trk_pt,   &b_trk_pt);
  cmsswTree->SetBranchAddress("trk_eta",  &trk_eta,  &b_trk_eta);
  cmsswTree->SetBranchAddress("trk_phi",  &trk_phi,  &b_trk_phi);
  cmsswTree->SetBranchAddress("trk_z0",   &trk_z0,   &b_trk_z0);
  cmsswTree->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
  cmsswTree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);


  //Initialize histograms

  TH1F* h_tpPt  = new TH1F("h_tpPt", "tp Pt",200,0.,200.);
  TH1F* h_tpPt_low  = new TH1F("h_tpPt_low", "tp Pt",50,0.,10.);
  TH1F* h_tpPhi = new TH1F("h_tpPhi","tp Phi",100,-TMath::Pi(),TMath::Pi());
  TH1F* h_tpEta = new TH1F("h_tpEta","tp Eta",100,-4.,4.);
  TH1F* h_tpZ0  = new TH1F("h_tpZ0", "tp Z0",200,-50.,50.);
  TH1F* h_tpVr  = new TH1F("h_tpVr", "tp sq(vx^2+vy^2)",2000,0.,2.);
  TH1F* h_tpType   = new TH1F("h_tpType",  "tp type",801,-400.5,400.5);

  TH1F* h_trkPt_presel = new TH1F("h_trkPt_presel","L1 track pt",200,0.,200);
  TH1F* h_trkPt_presel_low = new TH1F("h_trkPt_presel_low","L1 track pt",50,0.,10);
  TH1F* h_trkPhi_presel = new TH1F("h_trkPhi_presel","L1 track phi",100,-TMath::Pi(),TMath::Pi());
  TH1F* h_trkEta_presel = new TH1F("h_trkEta_presel","L1 track Eta",100,-4.,4.);
  TH1F* h_trkZ0_presel = new TH1F("h_trkZ0_presel","L1 track Z0",200,-50.,50.);

  TH1F* h_trkPt        = new TH1F("h_trkPt","L1 track pt",200,0,200);
  TH1F* h_trkPt_low    = new TH1F("h_trkPt_low","L1 track pt",50,0,10);
  TH1F* h_trkPhi       = new TH1F("h_trkPhi","L1 track (global) phi",100,-TMath::Pi(),TMath::Pi());
  TH1F* h_trkEta       = new TH1F("h_trkEta","L1 track eta",100,-4.,4.);
  TH1F* h_trkZ0        = new TH1F("h_trkZ0","L1 track z0",200,-50.,50.);

  TH1F* h_chisq             = new TH1F("h_chisq","chisq of L1 tracks",100,0,7);
  TH1F* h_ichisq            = new TH1F("h_ichisq","ichisq of L1 tracks",100,0,7);
    
  TH1F* h_nStubTrk     = new TH1F("h_nStubTrk","Num Stubs per Track",8,-0.5,8.5);

  //For Single Muons
  TH1F* h_nTrkEvt      = new TH1F("h_nTrkEvt","Num Tracks per Event",16,-0.5,15.5);
    
  /*
  //For ttbar
  TH1F* h_nTrkEvt      = new TH1F("h_nTrkEvt","Num Tracks per Event",51,-0.5,50.5);
  */

  TH1F* h_allmatchDeltaPhi   = new TH1F("h_allmatchDeltaPhi_AllM","L1-tp delta phi",100,-0.003,0.003);
  TH1F* h_allmatchDeltaEta   = new TH1F("h_allmatchDeltaEta_AllM","L1-tp delta eta",100,-0.025,0.025);
  TH1F* h_allmatchDeltaZ0    = new TH1F("h_allmatchDeltaZ0_AllM"  ,"L1-tp delta Z0"  ,100,-2.5,2.5);
  TH1F* h_allmatchDeltaPtOPt = new TH1F("h_allmatchDeltaPtOPt_AllM" ,"L1-tp (delta Pt)/Pt" ,100,-0.1,0.1);

  TH1F* h_bestmatchDeltaPhi   = new TH1F("h_bestmatchDeltaPhi","L1-tp delta phi",100,-0.003,0.003);
  TH1F* h_bestmatchDeltaEta   = new TH1F("h_bestmatchDeltaEta","L1-tp delta eta",100,-0.025,0.025);
  TH1F* h_bestmatchDeltaZ0    = new TH1F("h_bestmatchDeltaZ0"  ,"L1-tp delta Z0"  ,100,-2.5,2.5);
  TH1F* h_bestmatchDeltaPtOPt = new TH1F("h_bestmatchDeltaPtOPt" ,"L1-tp (delta Pt)/Pt" ,100,-0.1,0.1);

  TH2F* h_allmatchDeltaPhi_chi2   = new TH2F("h_allmatchDeltaPhi_chi2_AllM","L1-tp delta phi v. chi2",40,-0.003,0.003,40,0,20);
  TH2F* h_allmatchDeltaEta_chi2   = new TH2F("h_allmatchDeltaEta_chi2_AllM","L1-tp delta eta v. chi2",40,-0.025,0.025,40,0,20);
  TH2F* h_allmatchDeltaZ0_chi2    = new TH2F("h_allmatchDeltaZ0_chi2_AllM"  ,"L1-tp delta Z0 v. chi2"  ,40,-4,4,40,0,20);
  TH2F* h_allmatchDeltaPtOPt_chi2 = new TH2F("h_allmatchDeltaPtOPt_chi2_AllM" ,"L1-tp (delta Pt)/Pt v. chi2" ,40,-0.1,0.1,40,0,20);

  TH1F* h_allmatchDeltaPhi_chi2_0to2   = new TH1F("h_allmatchDeltaPhi_chi2_0to2","L1-tp delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_allmatchDeltaEta_chi2_0to2   = new TH1F("h_allmatchDeltaEta_chi2_0to2","L1-tp delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_allmatchDeltaZ0_chi2_0to2    = new TH1F("h_allmatchDeltaZ0_chi2_0to2"  ,"L1-tp delta Z0 v. chi2"  ,40,-4,4);
  TH1F* h_allmatchDeltaPtOPt_chi2_0to2 = new TH1F("h_allmatchDeltaPtOPt_chi2_0to2" ,"L1-tp (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_allmatchDeltaPhi_chi2_2to3   = new TH1F("h_allmatchDeltaPhi_chi2_2to3","L1-tp delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_allmatchDeltaEta_chi2_2to3   = new TH1F("h_allmatchDeltaEta_chi2_2to3","L1-tp delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_allmatchDeltaZ0_chi2_2to3    = new TH1F("h_allmatchDeltaZ0_chi2_2to3"  ,"L1-tp delta Z0 v. chi2"  ,40,-4,4);
  TH1F* h_allmatchDeltaPtOPt_chi2_2to3 = new TH1F("h_allmatchDeltaPtOPt_chi2_2to3" ,"L1-tp (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_allmatchDeltaPhi_chi2_3to4   = new TH1F("h_allmatchDeltaPhi_chi2_3to4","L1-tp delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_allmatchDeltaEta_chi2_3to4   = new TH1F("h_allmatchDeltaEta_chi2_3to4","L1-tp delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_allmatchDeltaZ0_chi2_3to4    = new TH1F("h_allmatchDeltaZ0_chi2_3to4"  ,"L1-tp delta Z0 v. chi2"  ,40,-4,4);
  TH1F* h_allmatchDeltaPtOPt_chi2_3to4 = new TH1F("h_allmatchDeltaPtOPt_chi2_3to4" ,"L1-tp (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH1F* h_allmatchDeltaPhi_chi2_gr4   = new TH1F("h_allmatchDeltaPhi_chi2_gr4","L1-tp delta phi v. chi2",40,-0.003,0.003);
  TH1F* h_allmatchDeltaEta_chi2_gr4   = new TH1F("h_allmatchDeltaEta_chi2_gr4","L1-tp delta eta v. chi2",40,-0.025,0.025);
  TH1F* h_allmatchDeltaZ0_chi2_gr4    = new TH1F("h_allmatchDeltaZ0_chi2_gr4"  ,"L1-tp delta Z0 v. chi2"  ,40,-4,4);
  TH1F* h_allmatchDeltaPtOPt_chi2_gr4 = new TH1F("h_allmatchDeltaPtOPt_chi2_gr4" ,"L1-tp (delta Pt)/Pt v. chi2" ,40,-0.1,0.1);

  TH2F* h_nMatchTrkVsTpPt      = new TH2F("h_nMatchTrkVsTpPt", "nMatchTrk vs tp Pt",600,0.,150.,10,-0.5,9.5);
  TH2F* h_nMatchTrkVsTpPt_low  = new TH2F("h_nMatchTrkVsTpPt_low", "nMatchTrk vs tp Pt",40,0.,10.,10,-0.5,9.5);

  // Loop over events
  for(int i=0; i<numCMSSWEvts; i++) {
    cmsswTree->GetEntry(i);

    unsigned int nTrkEvt=trk_pt->size();
    h_nTrkEvt->Fill(nTrkEvt);

    // Loop over L1 tracks
    for(unsigned int j=0; j<nTrkEvt; j++) {

      h_trkPt_presel->Fill(trk_pt->at(j));
      h_trkPt_presel_low->Fill(trk_pt->at(j));
      h_trkEta_presel->Fill(trk_eta->at(j));
      h_trkPhi_presel->Fill(trk_phi->at(j));
      h_trkZ0_presel->Fill(trk_z0->at(j));

      if(trk_nstub->at(j) < minStubs) continue;
      if(trk_pt->at(j) < Minpt || trk_pt->at(j) > Maxpt) continue;

      h_trkPt->Fill(trk_pt->at(j));
      h_trkPt_low->Fill(trk_pt->at(j));
      h_trkEta->Fill(trk_eta->at(j));
      h_trkPhi->Fill(trk_phi->at(j));
      h_trkZ0->Fill(trk_z0->at(j));

      h_nStubTrk->Fill(trk_nstub->at(j));
      h_chisq->Fill(trk_chi2->at(j));

    }

    // Loop over tracking particles
    unsigned int nTpEvt = tp_pt->size();
    for(unsigned int j=0; j<nTpEvt; j++) {

      h_tpPt->Fill(tp_pt->at(j));
      h_tpPt_low->Fill(tp_pt->at(j));
      h_tpPhi->Fill(tp_phi->at(j));
      h_tpEta->Fill(tp_eta->at(j));
      h_tpZ0->Fill(tp_z0->at(j));
      h_tpVr->Fill(tp_dxy->at(j));
      h_tpType->Fill(tp_pdgid->at(j));

      h_nMatchTrkVsTpPt->Fill(tp_pt->at(j), allmatchtrk_pt->at(j).size());
      h_nMatchTrkVsTpPt_low->Fill(tp_pt->at(j), allmatchtrk_pt->at(j).size());

      h_bestmatchDeltaPtOPt->Fill((matchtrk_pt->at(j) - tp_pt->at(j))/tp_pt->at(j));
      h_bestmatchDeltaPhi->Fill(matchtrk_phi->at(j) - tp_phi->at(j));
      h_bestmatchDeltaEta->Fill(matchtrk_eta->at(j) - tp_eta->at(j));
      h_bestmatchDeltaZ0->Fill(matchtrk_z0->at(j) - tp_z0->at(j));

      // Loop over all trks matched to tp
      unsigned int nMatchTrk = allmatchtrk_pt->at(j).size();
      for(unsigned int k=0; k<nMatchTrk; k++) {

        h_allmatchDeltaPtOPt->Fill((allmatchtrk_pt->at(j).at(k) - tp_pt->at(j))/tp_pt->at(j));
        h_allmatchDeltaPhi->Fill(allmatchtrk_phi->at(j).at(k) - tp_phi->at(j));
        h_allmatchDeltaEta->Fill(allmatchtrk_eta->at(j).at(k) - tp_eta->at(j));
        h_allmatchDeltaZ0->Fill(allmatchtrk_z0->at(j).at(k) - tp_z0->at(j));

        h_allmatchDeltaPtOPt_chi2->Fill(allmatchtrk_chi2->at(j).at(k),(allmatchtrk_pt->at(j).at(k) - tp_pt->at(j))/tp_pt->at(j));
        h_allmatchDeltaPhi_chi2->Fill(allmatchtrk_chi2->at(j).at(k),allmatchtrk_phi->at(j).at(k) - tp_phi->at(j));
        h_allmatchDeltaEta_chi2->Fill(allmatchtrk_chi2->at(j).at(k),allmatchtrk_eta->at(j).at(k) - tp_eta->at(j));
        h_allmatchDeltaZ0_chi2->Fill(allmatchtrk_chi2->at(j).at(k),allmatchtrk_z0->at(j).at(k) - tp_z0->at(j));

        if(allmatchtrk_chi2->at(j).at(k) > 0 && allmatchtrk_chi2->at(j).at(k) < 2) {
          h_allmatchDeltaPtOPt_chi2_0to2->Fill((allmatchtrk_pt->at(j).at(k) - tp_pt->at(j))/tp_pt->at(j));
          h_allmatchDeltaPhi_chi2_0to2->Fill(allmatchtrk_phi->at(j).at(k) - tp_phi->at(j));
          h_allmatchDeltaEta_chi2_0to2->Fill(allmatchtrk_eta->at(j).at(k) - tp_eta->at(j));
          h_allmatchDeltaZ0_chi2_0to2->Fill(allmatchtrk_z0->at(j).at(k) - tp_z0->at(j));
        }
        if(allmatchtrk_chi2->at(j).at(k) > 2 && allmatchtrk_chi2->at(j).at(k) < 3) {
          h_allmatchDeltaPtOPt_chi2_2to3->Fill((allmatchtrk_pt->at(j).at(k) - tp_pt->at(j))/tp_pt->at(j));
          h_allmatchDeltaPhi_chi2_2to3->Fill(allmatchtrk_phi->at(j).at(k) - tp_phi->at(j));
          h_allmatchDeltaEta_chi2_2to3->Fill(allmatchtrk_eta->at(j).at(k) - tp_eta->at(j));
          h_allmatchDeltaZ0_chi2_2to3->Fill(allmatchtrk_z0->at(j).at(k) - tp_z0->at(j));
        }
        if(allmatchtrk_chi2->at(j).at(k) > 3 && allmatchtrk_chi2->at(j).at(k) < 4) {
          h_allmatchDeltaPtOPt_chi2_3to4->Fill((allmatchtrk_pt->at(j).at(k) - tp_pt->at(j))/tp_pt->at(j));
          h_allmatchDeltaPhi_chi2_3to4->Fill(allmatchtrk_phi->at(j).at(k) - tp_phi->at(j));
          h_allmatchDeltaEta_chi2_3to4->Fill(allmatchtrk_eta->at(j).at(k) - tp_eta->at(j));
          h_allmatchDeltaZ0_chi2_3to4->Fill(allmatchtrk_z0->at(j).at(k) - tp_z0->at(j));
        }
        if(allmatchtrk_chi2->at(j).at(k) > 4) {
          h_allmatchDeltaPtOPt_chi2_gr4->Fill((allmatchtrk_pt->at(j).at(k) - tp_pt->at(j))/tp_pt->at(j));
          h_allmatchDeltaPhi_chi2_gr4->Fill(allmatchtrk_phi->at(j).at(k) - tp_phi->at(j));
          h_allmatchDeltaEta_chi2_gr4->Fill(allmatchtrk_eta->at(j).at(k) - tp_eta->at(j));
          h_allmatchDeltaZ0_chi2_gr4->Fill(allmatchtrk_z0->at(j).at(k) - tp_z0->at(j));
        }

      }
    }

  }//end loop over events
  
  // Write and close file
  histFile->Write();
  histFile->Close();

}

