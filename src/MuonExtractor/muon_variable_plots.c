#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void muon_variable_plots() {
    TFile *file = TFile::Open("Muons.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "TTree 'Events' not found in file!" << std::endl;
        return;
    }

    // Branch variables
    UInt_t nMuon;
    Float_t Muon_pt[10];
    Float_t Muon_eta[10];
    Float_t Muon_phi[10];
    Float_t Muon_mass[10];
    Float_t Muon_dxy[10];
    Float_t Muon_dz[10];

    // Set branch addresses
    tree->SetBranchAddress("nMuon", &nMuon);
    tree->SetBranchAddress("Muon_pt", Muon_pt);
    tree->SetBranchAddress("Muon_eta", Muon_eta);
    tree->SetBranchAddress("Muon_phi", Muon_phi);
    tree->SetBranchAddress("Muon_mass", Muon_mass);
    tree->SetBranchAddress("Muon_dxy", Muon_dxy);
    tree->SetBranchAddress("Muon_dz", Muon_dz);

    // Histograms
    TH1F *h_muon_pt     = new TH1F("h_muon_pt", "Muon pT; pT [GeV]; Entries", 50, 0, 100);
    TH1F *h_muon_eta    = new TH1F("h_muon_eta", "Muon eta; #eta; Entries", 50, -3, 3);
    TH1F *h_muon_phi    = new TH1F("h_muon_phi", "Muon phi; #phi; Entries", 50, -3.5, 3.5);
    TH1F *h_muon_mass   = new TH1F("h_muon_mass", "Muon mass; mass [GeV]; Entries", 50, 0, 0.2);
    TH1F *h_dxy         = new TH1F("h_dxy", "Muon dxy; dxy [cm]; Entries", 50, -0.5, 0.5);
    TH1F *h_dz          = new TH1F("h_dz", "Muon dz; dz [cm]; Entries", 50, -1.0, 1.0);

    TH2F *h_eta_phi     = new TH2F("h_eta_phi", "Muon eta vs phi; #eta; #phi", 50, -3, 3, 50, -3.5, 3.5);
    TH2F *h_dxy_dz      = new TH2F("h_dxy_dz", "Muon dxy vs dz; dxy [cm]; dz [cm]", 50, -0.5, 0.5, 50, -1.0, 1.0);

    // Loop over events and muons
    Long64_t nentries = tree->GetEntries();
    int passed_events = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        if (nMuon > 0) passed_events++;  //  at least one muon

        for (UInt_t j = 0; j < nMuon; ++j) {
            h_muon_pt->Fill(Muon_pt[j]);
            h_muon_eta->Fill(Muon_eta[j]);
            h_muon_phi->Fill(Muon_phi[j]);
            h_muon_mass->Fill(Muon_mass[j]);
            h_dxy->Fill(Muon_dxy[j]);
            h_dz->Fill(Muon_dz[j]);
            h_eta_phi->Fill(Muon_eta[j], Muon_phi[j]);
            h_dxy_dz->Fill(Muon_dxy[j], Muon_dz[j]);
        }
    }

    std::cout << "Total events: " << nentries << std::endl;
    std::cout << "Events passed: " << passed_events << std::endl;
    std::cout<<"\n\n\n\n"<<std::endl;

    // odf output
    TCanvas *c1 = new TCanvas("c1", "1D Histograms", 1000, 800);
    c1->Divide(3, 2);
    c1->cd(1); h_muon_pt->Draw();
    c1->cd(2); h_muon_eta->Draw();
    c1->cd(3); h_muon_phi->Draw();
    c1->cd(4); h_muon_mass->Draw();
    c1->cd(5); h_dxy->Draw();
    c1->cd(6); h_dz->Draw();

    TCanvas *c2 = new TCanvas("c2", "2D Histograms", 1000, 400);
    c2->Divide(2, 1);
    c2->cd(1); h_eta_phi->Draw("COLZ");
    c2->cd(2); h_dxy_dz->Draw("COLZ");

    // Saving pdf file
    c1->Print("muon_histograms.pdf[");  
    c1->Print("muon_histograms.pdf");   
    c2->Print("muon_histograms.pdf");   
    c2->Print("muon_histograms.pdf]");  

    // Save histograms to ROOT file
    TFile *outfile = new TFile("histograms.root", "RECREATE");
    h_muon_pt->Write();
    h_muon_eta->Write();
    h_muon_phi->Write();
    h_muon_mass->Write();
    h_dxy->Write();
    h_dz->Write();
    h_eta_phi->Write();
    h_dxy_dz->Write();
    outfile->Close();

    std::cout << "\n\n\nHistograms saved!!!" << std::endl;
}

