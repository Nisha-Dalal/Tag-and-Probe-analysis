#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h> 
#include <TMath.h>
#include <cmath>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>

void tag_and_probe_analysis() {
    // uses the Muons.root file
    TFile *file = TFile::Open("Muons.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening Muons.root!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "TTree 'Events' not found!" << std::endl;
        return;
    }

    // Connect branches 
    UInt_t nMuon;
    Float_t Muon_pt[10];
    Float_t Muon_eta[10];
    Float_t Muon_phi[10];
    Float_t Muon_mass[10];
    Float_t Muon_dxy[10];
    Float_t Muon_dz[10];
    Float_t Muon_pfRelIso03_all[10];
    Int_t Muon_charge[10];                  
    Bool_t Muon_tightId[10];                
    Bool_t Muon_mediumId[10];               // Adding medium ID
    Bool_t Muon_looseId[10];                // Adding loose ID
    Int_t Muon_genPartIdx[10];            
    Int_t GenPart_pdgId[200];               
    Int_t GenPart_genPartIdxMother[200];
    Int_t PV_npvs;

    // Set branch addresses
    tree->SetBranchAddress("nMuon", &nMuon);
    tree->SetBranchAddress("Muon_pt", Muon_pt);
    tree->SetBranchAddress("Muon_eta", Muon_eta);
    tree->SetBranchAddress("Muon_phi", Muon_phi);
    tree->SetBranchAddress("Muon_mass", Muon_mass);
    tree->SetBranchAddress("Muon_dxy", Muon_dxy);
    tree->SetBranchAddress("Muon_dz", Muon_dz);
    tree->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all);
    tree->SetBranchAddress("Muon_charge", Muon_charge);
    tree->SetBranchAddress("Muon_tightId", Muon_tightId);
    tree->SetBranchAddress("Muon_mediumId", Muon_mediumId);    // Medium ID
    tree->SetBranchAddress("Muon_looseId", Muon_looseId);      // Loose ID
    tree->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx);
    tree->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    tree->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    tree->SetBranchAddress("PV_npvs", &PV_npvs);


    // Define histograms with 10 bins from 10 to 110 GeV
    TH1F *h_totalprobe_pT = new TH1F("h_totalprobe_pT", "Total Probe Muon pT;pT [GeV];Entries", 20, 0, 100);
    TH1F *h_totalprobe_eta = new TH1F("h_totalprobe_eta", "Total Probe Muon eta;#eta;Entries", 40, -2.5, 2.5);
    TH1F *h_totalprobe_phi = new TH1F("h_totalprobe_phi", "Total Probe Muon phi;#phi;Entries", 30, -3.2, 3.2);
    TH1F *h_totalprobe_nvtx = new TH1F("h_totalprobe_nvtx", "Total Probe Muon PV;# of Vertices;Entries", 50, 0, 100);
    TH1F *h_totalprobe_Mll = new TH1F("h_totalprobe_Mll", "Total Probe Mll;M_{ll} [GeV];Entries", 60, 60, 120);

    // Tight ID histograms
    TH1F *h_passprobe_tight_pT = new TH1F("h_passprobe_tight_pT", "Passing Probe Muon (Tight ID) pT;pT [GeV];Entries", 20, 0, 100);
    TH1F *h_passprobe_tight_eta = new TH1F("h_passprobe_tight_eta", "Passing Probe Muon (Tight ID) eta;#eta;Entries", 40, -2.5, 2.5);
    TH1F *h_passprobe_tight_phi = new TH1F("h_passprobe_tight_phi", "Passing Probe Muon (Tight ID) phi;#phi;Entries", 30, -3.2, 3.2);
    TH1F *h_passprobe_tight_nvtx = new TH1F("h_passprobe_tight_nvtx", "Passing Probe Muon (Tight ID) PV;# of Vertices;Entries", 50, 0, 100);
    TH1F *h_passprobe_tight_Mll = new TH1F("h_passprobe_tight_Mll", "Passing Probe (Tight ID) Mll;M_{ll} [GeV];Entries", 60, 60, 120);

    // Medium ID histograms
    TH1F *h_passprobe_medium_pT = new TH1F("h_passprobe_medium_pT", "Passing Probe Muon (Medium ID) pT;pT [GeV];Entries", 20, 0, 100);
    TH1F *h_passprobe_medium_eta = new TH1F("h_passprobe_medium_eta", "Passing Probe Muon (Medium ID) eta;#eta;Entries", 40, -2.5, 2.5);
    TH1F *h_passprobe_medium_phi = new TH1F("h_passprobe_medium_phi", "Passing Probe Muon (Medium ID) phi;#phi;Entries", 30, -3.2, 3.2);
    TH1F *h_passprobe_medium_nvtx = new TH1F("h_passprobe_medium_nvtx", "Passing Probe Muon (Medium ID) PV;# of Vertices;Entries", 50, 0, 100);
    TH1F *h_passprobe_medium_Mll = new TH1F("h_passprobe_medium_Mll", "Passing Probe (Medium ID) Mll;M_{ll} [GeV];Entries", 60, 60, 120);

    // Loose ID histograms
    TH1F *h_passprobe_loose_pT = new TH1F("h_passprobe_loose_pT", "Passing Probe Muon (Loose ID) pT;pT [GeV];Entries",20, 0, 100);
    TH1F *h_passprobe_loose_eta = new TH1F("h_passprobe_loose_eta", "Passing Probe Muon (Loose ID) eta;#eta;Entries", 40, -2.5, 2.5);
    TH1F *h_passprobe_loose_phi = new TH1F("h_passprobe_loose_phi", "Passing Probe Muon (Loose ID) phi;#phi;Entries", 30, -3.2, 3.2);
    TH1F *h_passprobe_loose_nvtx = new TH1F("h_passprobe_loose_nvtx", "Passing Probe Muon (Loose ID) PV;# of Vertices;Entries", 50, 0, 100);
    TH1F *h_passprobe_loose_Mll = new TH1F("h_passprobe_loose_Mll", "Passing Probe (Loose ID) Mll;M_{ll} [GeV];Entries", 60, 60, 120);

    // 2D histograms
    TH2F *h_totalprobe_pt_eta = new TH2F("h_totalprobe_pt_eta", "Total Probe pT vs |#eta|;|#eta|;pT [GeV]", 20, 0, 2.5, 20, 0, 50);
    TH2F *h_passprobe_tight_pt_eta  = new TH2F("h_passprobe_tight_pt_eta",  "Passing Probe (Tight ID) pT vs |#eta|;|#eta|;pT [GeV]", 20, 0, 2.5, 20, 0, 50);
    TH2F *h_passprobe_medium_pt_eta = new TH2F("h_passprobe_medium_pt_eta", "Passing Probe (Medium ID) pT vs |#eta|;|#eta|;pT [GeV]", 20, 0, 2.5, 20, 0, 50);
    TH2F *h_passprobe_loose_pt_eta  = new TH2F("h_passprobe_loose_pt_eta",  "Passing Probe (Loose ID) pT vs |#eta|;|#eta|;pT [GeV]",  20, 0, 2.5, 20, 0, 50);

    // Isolation histograms (with ID requirements)
    TH1F *h_passprobe_tightiso_pT   = new TH1F("h_passprobe_tightiso_pT", "Passing Probe (Tight ID+Iso) pT;pT [GeV];Entries", 20, 0, 100);
    TH1F *h_passprobe_tightiso_eta  = new TH1F("h_passprobe_tightiso_eta", "Passing Probe (Tight ID+Iso) eta;#eta;Entries", 40, -2.5, 2.5);
    TH1F *h_passprobe_tightiso_phi  = new TH1F("h_passprobe_tightiso_phi", "Passing Probe (Tight ID+Iso) phi;#phi;Entries", 30, -3.2, 3.2);
    TH1F *h_passprobe_tightiso_nvtx = new TH1F("h_passprobe_tightiso_nvtx", "Passing Probe (Tight ID+Iso) PV;# of Vertices;Entries", 50, 0, 100);
    
    TH1F *h_passprobe_mediumiso_pT   = new TH1F("h_passprobe_mediumiso_pT", "Passing Probe (Medium ID+Iso) pT;pT [GeV];Entries",20, 0, 100);
    TH1F *h_passprobe_mediumiso_eta  = new TH1F("h_passprobe_mediumiso_eta", "Passing Probe (Medium ID+Iso) eta;#eta;Entries", 40, -2.5, 2.5);
    TH1F *h_passprobe_mediumiso_phi  = new TH1F("h_passprobe_mediumiso_phi", "Passing Probe (Medium ID+Iso) phi;#phi;Entries", 30, -3.2, 3.2);
    TH1F *h_passprobe_mediumiso_nvtx = new TH1F("h_passprobe_mediumiso_nvtx", "Passing Probe (Medium ID+Iso) PV;# of Vertices;Entries", 50, 0, 100);

    TH1F *h_passprobe_looseiso_pT   = new TH1F("h_passprobe_looseiso_pT", "Passing Probe (Loose ID+Iso) pT;pT [GeV];Entries", 20, 0, 100);
    TH1F *h_passprobe_looseiso_eta  = new TH1F("h_passprobe_looseiso_eta", "Passing Probe (Loose ID+Iso) eta;#eta;Entries", 40, -2.5, 2.5);
    TH1F *h_passprobe_looseiso_phi  = new TH1F("h_passprobe_looseiso_phi", "Passing Probe (Loose ID+Iso) phi;#phi;Entries", 30, -3.2, 3.2);
    TH1F *h_passprobe_looseiso_nvtx = new TH1F("h_passprobe_looseiso_nvtx", "Passing Probe (Loose ID+Iso) PV;# of Vertices;Entries", 50, 0, 100);  // Failed Probe
    TH1F *h_failedprobe_Mll = new TH1F("h_failedprobe_Mll", "Failed Probe Mll;M_{ll} [GeV];Entries", 50, 50, 150);

    // Event counters  
    Long64_t total_events = 0;
    Long64_t tight_pass_events = 0;
    Long64_t medium_pass_events = 0;
    Long64_t loose_pass_events = 0;
    Long64_t tightiso_pass_events = 0;
    Long64_t mediumiso_pass_events = 0;
    Long64_t looseiso_pass_events = 0;

    // Event loop
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (nMuon < 2) continue;

        int tagIndex = -1;
        for (UInt_t j = 0; j < nMuon; ++j) {
            int genIdx = Muon_genPartIdx[j];
            if (Muon_pt[j] > 25 &&
                Muon_tightId[j] == 1 &&
                std::abs(Muon_eta[j]) < 2.4 &&
                Muon_pfRelIso03_all[j] < 0.15 &&
                genIdx >= 0 &&
                GenPart_genPartIdxMother[genIdx] >= 0 &&
                std::abs(GenPart_pdgId[GenPart_genPartIdxMother[genIdx]]) == 23) {
                tagIndex = j;
                break;
            }
        }

        if (tagIndex == -1) continue;

        for (UInt_t k = 0; k < nMuon; ++k) {
            if (k == tagIndex) continue;
            if (Muon_charge[tagIndex] * Muon_charge[k] >= 0) continue;

            int genIdx = Muon_genPartIdx[k];
            if (Muon_pt[k] > 1 &&
                std::abs(Muon_eta[k]) < 2.4 &&
                genIdx >= 0 &&
                GenPart_genPartIdxMother[genIdx] >= 0 &&
                std::abs(GenPart_pdgId[GenPart_genPartIdxMother[genIdx]]) == 23) {

                const float m_mu = 0.105;
                float px1 = Muon_pt[tagIndex] * cos(Muon_phi[tagIndex]);
                float py1 = Muon_pt[tagIndex] * sin(Muon_phi[tagIndex]);
                float pz1 = Muon_pt[tagIndex] * sinh(Muon_eta[tagIndex]);
                float e1  = sqrt(px1*px1 + py1*py1 + pz1*pz1 + m_mu*m_mu);

                float px2 = Muon_pt[k] * cos(Muon_phi[k]);
                float py2 = Muon_pt[k] * sin(Muon_phi[k]);
                float pz2 = Muon_pt[k] * sinh(Muon_eta[k]);
                float e2  = sqrt(px2*px2 + py2*py2 + pz2*pz2 + m_mu*m_mu);

                float Mll = sqrt((e1 + e2)*(e1 + e2) -
                                 (px1 + px2)*(px1 + px2) -
                                 (py1 + py2)*(py1 + py2) -
                                 (pz1 + pz2)*(pz1 + pz2));

                if (Mll > 60 && Mll < 120) {
                    total_events++;
                    
                    // Fill total probe histograms
                    h_totalprobe_pT->Fill(Muon_pt[k]);
                    h_totalprobe_eta->Fill(Muon_eta[k]);
                    h_totalprobe_phi->Fill(Muon_phi[k]);
                    h_totalprobe_nvtx->Fill(PV_npvs);
                    h_totalprobe_Mll->Fill(Mll);
                    h_totalprobe_pt_eta->Fill(Muon_eta[k], Muon_pt[k]);

                    bool passedAny = false;

                    // Tight ID selection
                    if (Muon_tightId[k] == 1) {
                        passedAny = true;
                        tight_pass_events++;
                        h_passprobe_tight_pT->Fill(Muon_pt[k]);
                        h_passprobe_tight_eta->Fill(Muon_eta[k]);
                        h_passprobe_tight_phi->Fill(Muon_phi[k]);
                        h_passprobe_tight_nvtx->Fill(PV_npvs);
                        h_passprobe_tight_Mll->Fill(Mll);
                        h_passprobe_tight_pt_eta->Fill(Muon_eta[k], Muon_pt[k]);
                        
                        if (Muon_pfRelIso03_all[k] < 0.15) {
                            tightiso_pass_events++;
                            h_passprobe_tightiso_pT->Fill(Muon_pt[k]);
                            h_passprobe_tightiso_eta->Fill(Muon_eta[k]);
                            h_passprobe_tightiso_phi->Fill(Muon_phi[k]);
                            h_passprobe_tightiso_nvtx->Fill(PV_npvs);
                        }
                    }

                    // Medium ID selection
                    if (Muon_mediumId[k] == 1) {
                        passedAny = true;
                        medium_pass_events++;
                        h_passprobe_medium_pT->Fill(Muon_pt[k]);
                        h_passprobe_medium_eta->Fill(Muon_eta[k]);
                        h_passprobe_medium_phi->Fill(Muon_phi[k]);
                        h_passprobe_medium_nvtx->Fill(PV_npvs);
                        h_passprobe_medium_Mll->Fill(Mll);
                        h_passprobe_medium_pt_eta->Fill(Muon_eta[k], Muon_pt[k]);
                        
                        if (Muon_pfRelIso03_all[k] < 0.15) {
                            mediumiso_pass_events++;
                            h_passprobe_mediumiso_pT->Fill(Muon_pt[k]);
                            h_passprobe_mediumiso_eta->Fill(Muon_eta[k]);
                            h_passprobe_mediumiso_phi->Fill(Muon_phi[k]);
                            h_passprobe_mediumiso_nvtx->Fill(PV_npvs);
                        }
                    }

                    // Loose ID selection
                    if (Muon_looseId[k] == 1) {
                        passedAny = true;
                        loose_pass_events++;
                        h_passprobe_loose_pT->Fill(Muon_pt[k]);
                        h_passprobe_loose_eta->Fill(Muon_eta[k]);
                        h_passprobe_loose_phi->Fill(Muon_phi[k]);
                        h_passprobe_loose_nvtx->Fill(PV_npvs);
                        h_passprobe_loose_Mll->Fill(Mll);
                        h_passprobe_loose_pt_eta->Fill(Muon_eta[k], Muon_pt[k]);
                        
                        if (Muon_pfRelIso03_all[k] < 0.15) {
                            looseiso_pass_events++;
                            h_passprobe_looseiso_pT->Fill(Muon_pt[k]);
                            h_passprobe_looseiso_eta->Fill(Muon_eta[k]);
                            h_passprobe_looseiso_phi->Fill(Muon_phi[k]);
                            h_passprobe_looseiso_nvtx->Fill(PV_npvs);
                        }
                    }
                    if (!passedAny) {
                        h_failedprobe_Mll->Fill(Mll);
                    }
                }
            }
        }
    }

    // Print event counts
    std::cout << "\n===== Event Counts =====\n";
    std::cout << "Total probe events: " << total_events << "\n";
    std::cout << "Tight ID passing events: " << tight_pass_events << " (" << 100.0*tight_pass_events/total_events << "%)\n";
    std::cout << "Medium ID passing events: " << medium_pass_events << " (" << 100.0*medium_pass_events/total_events << "%)\n";
    std::cout << "Loose ID passing events: " << loose_pass_events << " (" << 100.0*loose_pass_events/total_events << "%)\n";
    std::cout << "Tight ID+ISO passing events: " << tightiso_pass_events << " (" << 100.0*tightiso_pass_events/total_events << "%)\n";
    std::cout << "Medium ID+ISO passing events: " << mediumiso_pass_events << " (" << 100.0*mediumiso_pass_events/total_events << "%)\n";
    std::cout << "Loose ID+ISO passing events: " << looseiso_pass_events << " (" << 100.0*looseiso_pass_events/total_events << "%)\n";

      // Create efficiency graphs using TGraphAsymmErrors
    TGraphAsymmErrors *g_eff_tight_pT = new TGraphAsymmErrors(h_passprobe_tight_pT, h_totalprobe_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tight_eta = new TGraphAsymmErrors(h_passprobe_tight_eta, h_totalprobe_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tight_phi = new TGraphAsymmErrors(h_passprobe_tight_phi, h_totalprobe_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tight_nvtx = new TGraphAsymmErrors(h_passprobe_tight_nvtx, h_totalprobe_nvtx, "cl=0.683 b(1,1) mode");

    TGraphAsymmErrors *g_eff_medium_pT = new TGraphAsymmErrors(h_passprobe_medium_pT, h_totalprobe_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_medium_eta = new TGraphAsymmErrors(h_passprobe_medium_eta, h_totalprobe_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_medium_phi = new TGraphAsymmErrors(h_passprobe_medium_phi, h_totalprobe_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_medium_nvtx = new TGraphAsymmErrors(h_passprobe_medium_nvtx, h_totalprobe_nvtx, "cl=0.683 b(1,1) mode");

    TGraphAsymmErrors *g_eff_loose_pT = new TGraphAsymmErrors(h_passprobe_loose_pT, h_totalprobe_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_loose_eta = new TGraphAsymmErrors(h_passprobe_loose_eta, h_totalprobe_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_loose_phi = new TGraphAsymmErrors(h_passprobe_loose_phi, h_totalprobe_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_loose_nvtx = new TGraphAsymmErrors(h_passprobe_loose_nvtx, h_totalprobe_nvtx, "cl=0.683 b(1,1) mode");

    // Isolation efficiency graphs (relative to ID passing)
    TGraphAsymmErrors *g_eff_tightiso_pT = new TGraphAsymmErrors(h_passprobe_tightiso_pT, h_passprobe_tight_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tightiso_eta = new TGraphAsymmErrors(h_passprobe_tightiso_eta, h_passprobe_tight_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tightiso_phi = new TGraphAsymmErrors(h_passprobe_tightiso_phi, h_passprobe_tight_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tightiso_nvtx = new TGraphAsymmErrors(h_passprobe_tightiso_nvtx, h_passprobe_tight_nvtx, "cl=0.683 b(1,1) mode");

    TGraphAsymmErrors *g_eff_mediumiso_pT = new TGraphAsymmErrors(h_passprobe_mediumiso_pT, h_passprobe_medium_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_mediumiso_eta = new TGraphAsymmErrors(h_passprobe_mediumiso_eta, h_passprobe_medium_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_mediumiso_phi = new TGraphAsymmErrors(h_passprobe_mediumiso_phi, h_passprobe_medium_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_mediumiso_nvtx = new TGraphAsymmErrors(h_passprobe_mediumiso_nvtx, h_passprobe_medium_nvtx, "cl=0.683 b(1,1) mode");

    TGraphAsymmErrors *g_eff_looseiso_pT = new TGraphAsymmErrors(h_passprobe_looseiso_pT, h_passprobe_loose_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_looseiso_eta = new TGraphAsymmErrors(h_passprobe_looseiso_eta, h_passprobe_loose_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_looseiso_phi = new TGraphAsymmErrors(h_passprobe_looseiso_phi, h_passprobe_loose_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_looseiso_nvtx = new TGraphAsymmErrors(h_passprobe_looseiso_nvtx, h_passprobe_loose_nvtx, "cl=0.683 b(1,1) mode");

    // Set graph properties
    auto styleGraph = [](TGraphAsymmErrors* graph, Color_t color, Style_t marker, const char* title) {
        graph->SetLineColor(color);
        graph->SetMarkerColor(color);
        graph->SetMarkerStyle(marker);
        graph->SetMarkerSize(1.2);
        graph->SetLineWidth(2);
        graph->SetTitle(title);
        graph->GetYaxis()->SetTitle("Efficiency");
        graph->GetYaxis()->SetRangeUser(0, 1.1);
    };

    styleGraph(g_eff_tight_pT, kBlue+1, 20, "Tight ID Efficiency");
    styleGraph(g_eff_tight_eta, kBlue+1, 20, "Tight ID Efficiency");
    styleGraph(g_eff_tight_phi, kBlue+1, 20, "Tight ID Efficiency");
    styleGraph(g_eff_tight_nvtx, kBlue+1, 20, "Tight ID Efficiency");

    styleGraph(g_eff_medium_pT, kRed+1, 21, "Medium ID Efficiency");
    styleGraph(g_eff_medium_eta, kRed+1, 21, "Medium ID Efficiency");
    styleGraph(g_eff_medium_phi, kRed+1, 21, "Medium ID Efficiency");
    styleGraph(g_eff_medium_nvtx, kRed+1, 21, "Medium ID Efficiency");

    styleGraph(g_eff_loose_pT, kGreen+2, 22, "Loose ID Efficiency");
    styleGraph(g_eff_loose_eta, kGreen+2, 22, "Loose ID Efficiency");
    styleGraph(g_eff_loose_phi, kGreen+2, 22, "Loose ID Efficiency");
    styleGraph(g_eff_loose_nvtx, kGreen+2, 22, "Loose ID Efficiency");

    styleGraph(g_eff_tightiso_pT, kMagenta+1, 23, "Tight ID+ISO Efficiency");
    styleGraph(g_eff_tightiso_eta, kMagenta+1, 23, "Tight ID+ISO Efficiency");
    styleGraph(g_eff_tightiso_phi, kMagenta+1, 23, "Tight ID+ISO Efficiency");
    styleGraph(g_eff_tightiso_nvtx, kMagenta+1, 23, "Tight ID+ISO Efficiency");

    styleGraph(g_eff_mediumiso_pT, kOrange+7, 24, "Medium ID+ISO Efficiency");
    styleGraph(g_eff_mediumiso_eta, kOrange+7, 24, "Medium ID+ISO Efficiency");
    styleGraph(g_eff_mediumiso_phi, kOrange+7, 24, "Medium ID+ISO Efficiency");
    styleGraph(g_eff_mediumiso_nvtx, kOrange+7, 24, "Medium ID+ISO Efficiency");

    styleGraph(g_eff_looseiso_pT, kCyan+1, 25, "Loose ID+ISO Efficiency");
    styleGraph(g_eff_looseiso_eta, kCyan+1, 25, "Loose ID+ISO Efficiency");
    styleGraph(g_eff_looseiso_phi, kCyan+1, 25, "Loose ID+ISO Efficiency");
    styleGraph(g_eff_looseiso_nvtx, kCyan+1, 25, "Loose ID+ISO Efficiency");

    // Save histograms and graphs to file
    TFile *outfile = new TFile("D:/ROOTFiles/DataFiles/TagAndProbePlots.root", "RECREATE");
    
    // Write total probe histograms
    h_totalprobe_pT->Write();      
    h_totalprobe_eta->Write();     
    h_totalprobe_phi->Write();     
    h_totalprobe_nvtx->Write();    
    h_totalprobe_Mll->Write();     
    h_totalprobe_pt_eta->Write();  

    // Write tight ID histograms
    h_passprobe_tight_pT->Write();      
    h_passprobe_tight_eta->Write();     
    h_passprobe_tight_phi->Write();     
    h_passprobe_tight_nvtx->Write();    
    h_passprobe_tight_Mll->Write();     
    h_passprobe_tight_pt_eta->Write();  

    // Write medium ID histograms
    h_passprobe_medium_pT->Write();      
    h_passprobe_medium_eta->Write();     
    h_passprobe_medium_phi->Write();     
    h_passprobe_medium_nvtx->Write();    
    h_passprobe_medium_Mll->Write();     
    h_passprobe_medium_pt_eta->Write();  

    // Write loose ID histograms
    h_passprobe_loose_pT->Write();      
    h_passprobe_loose_eta->Write();     
    h_passprobe_loose_phi->Write();     
    h_passprobe_loose_nvtx->Write();    
    h_passprobe_loose_Mll->Write();     
    h_passprobe_loose_pt_eta->Write();  

    // Write failed probe histogram
    h_failedprobe_Mll->Write();

    // Write isolation histograms
    h_passprobe_tightiso_pT->Write();   
    h_passprobe_tightiso_eta->Write();  
    h_passprobe_tightiso_phi->Write();  
    h_passprobe_tightiso_nvtx->Write();

    h_passprobe_mediumiso_pT->Write();   
    h_passprobe_mediumiso_eta->Write();  
    h_passprobe_mediumiso_phi->Write();  
    h_passprobe_mediumiso_nvtx->Write();

    h_passprobe_looseiso_pT->Write();   
    h_passprobe_looseiso_eta->Write();  
    h_passprobe_looseiso_phi->Write();  
    h_passprobe_looseiso_nvtx->Write();

    // Write efficiency graphs
    g_eff_tight_pT->Write("g_eff_tight_pT");
    g_eff_tight_eta->Write("g_eff_tight_eta");
    g_eff_tight_phi->Write("g_eff_tight_phi");
    g_eff_tight_nvtx->Write("g_eff_tight_nvtx");

    g_eff_medium_pT->Write("g_eff_medium_pT");
    g_eff_medium_eta->Write("g_eff_medium_eta");
    g_eff_medium_phi->Write("g_eff_medium_phi");
    g_eff_medium_nvtx->Write("g_eff_medium_nvtx");

    g_eff_loose_pT->Write("g_eff_loose_pT");
    g_eff_loose_eta->Write("g_eff_loose_eta");
    g_eff_loose_phi->Write("g_eff_loose_phi");
    g_eff_loose_nvtx->Write("g_eff_loose_nvtx");

    g_eff_tightiso_pT->Write("g_eff_tightiso_pT");
    g_eff_tightiso_eta->Write("g_eff_tightiso_eta");
    g_eff_tightiso_phi->Write("g_eff_tightiso_phi");
    g_eff_tightiso_nvtx->Write("g_eff_tightiso_nvtx");

    g_eff_mediumiso_pT->Write("g_eff_mediumiso_pT");
    g_eff_mediumiso_eta->Write("g_eff_mediumiso_eta");
    g_eff_mediumiso_phi->Write("g_eff_mediumiso_phi");
    g_eff_mediumiso_nvtx->Write("g_eff_mediumiso_nvtx");

    g_eff_looseiso_pT->Write("g_eff_looseiso_pT");
    g_eff_looseiso_eta->Write("g_eff_looseiso_eta");
    g_eff_looseiso_phi->Write("g_eff_looseiso_phi");
    g_eff_looseiso_nvtx->Write("g_eff_looseiso_nvtx");

    outfile->Close();
    std::cout << "Histograms and efficiency graphs saved to TagAndProbePlots.root\n";

    // Function to plot individual efficiency graphs
    auto plotSingleEfficiency = [](TGraphAsymmErrors* graph, const std::string& filename, 
        const std::string& title, bool setPtRange = false) {
        TCanvas *c_single = new TCanvas("c_single", "Efficiency Canvas", 800, 600);
        c_single->SetGrid();

        graph->Draw("AP");
        graph->SetTitle(title.c_str());
        graph->GetYaxis()->SetRangeUser(0, 1.1);

        if(setPtRange) {
            graph->GetXaxis()->SetRangeUser(10, 110); // Set explicit range for pT plots
        }

        // Customize appearance
        graph->GetXaxis()->SetTitleSize(0.04);
        graph->GetYaxis()->SetTitleSize(0.04);
        graph->GetXaxis()->SetLabelSize(0.035);
        graph->GetYaxis()->SetLabelSize(0.035);

        c_single->SaveAs(filename.c_str());
        delete c_single;
    };

    // Plot individual efficiency graphs for tight ID
    plotSingleEfficiency(g_eff_tight_pT, "D:/ROOTFiles/DataFiles/Eff_Tight_pT.png", 
       "Tight ID Efficiency vs pT;pT [GeV];Efficiency", true);
    plotSingleEfficiency(g_eff_tight_eta, "D:/ROOTFiles/DataFiles/Eff_Tight_eta.png", 
       "Tight ID Efficiency vs #eta;#eta;Efficiency");
    plotSingleEfficiency(g_eff_tight_phi, "D:/ROOTFiles/DataFiles/Eff_Tight_phi.png", 
       "Tight ID Efficiency vs #phi;#phi;Efficiency");
    plotSingleEfficiency(g_eff_tight_nvtx, "D:/ROOTFiles/DataFiles/Eff_Tight_nvtx.png", 
       "Tight ID Efficiency vs N_{vertices};N_{vertices};Efficiency");

    // Plot individual efficiency graphs for medium ID
    plotSingleEfficiency(g_eff_medium_pT, "D:/ROOTFiles/DataFiles/Eff_Medium_pT.png", 
       "Medium ID Efficiency vs pT;pT [GeV];Efficiency", true);
    plotSingleEfficiency(g_eff_medium_eta, "D:/ROOTFiles/DataFiles/Eff_Medium_eta.png", 
       "Medium ID Efficiency vs #eta;#eta;Efficiency");
    plotSingleEfficiency(g_eff_medium_phi, "D:/ROOTFiles/DataFiles/Eff_Medium_phi.png", 
       "Medium ID Efficiency vs #phi;#phi;Efficiency");
    plotSingleEfficiency(g_eff_medium_nvtx, "D:/ROOTFiles/DataFiles/Eff_Medium_nvtx.png", 
       "Medium ID Efficiency vs N_{vertices};N_{vertices};Efficiency");

    // Plot individual efficiency graphs for loose ID
    plotSingleEfficiency(g_eff_loose_pT, "D:/ROOTFiles/DataFiles/Eff_Loose_pT.png", 
       "Loose ID Efficiency vs pT;pT [GeV];Efficiency", true);
    plotSingleEfficiency(g_eff_loose_eta, "D:/ROOTFiles/DataFiles/Eff_Loose_eta.png", 
       "Loose ID Efficiency vs #eta;#eta;Efficiency");
    plotSingleEfficiency(g_eff_loose_phi, "D:/ROOTFiles/DataFiles/Eff_Loose_phi.png", 
       "Loose ID Efficiency vs #phi;#phi;Efficiency");
    plotSingleEfficiency(g_eff_loose_nvtx, "D:/ROOTFiles/DataFiles/Eff_Loose_nvtx.png", 
       "Loose ID Efficiency vs N_{vertices};N_{vertices};Efficiency");

    // Plot individual efficiency graphs for tight ID+ISO
    plotSingleEfficiency(g_eff_tightiso_pT, "D:/ROOTFiles/DataFiles/Eff_TightIso_pT.png", 
       "Tight ID+ISO Efficiency vs pT;pT [GeV];Efficiency", true);
    plotSingleEfficiency(g_eff_tightiso_eta, "D:/ROOTFiles/DataFiles/Eff_TightIso_eta.png", 
       "Tight ID+ISO Efficiency vs #eta;#eta;Efficiency");
    plotSingleEfficiency(g_eff_tightiso_phi, "D:/ROOTFiles/DataFiles/Eff_TightIso_phi.png", 
       "Tight ID+ISO Efficiency vs #phi;#phi;Efficiency");
    plotSingleEfficiency(g_eff_tightiso_nvtx, "D:/ROOTFiles/DataFiles/Eff_TightIso_nvtx.png", 
       "Tight ID+ISO Efficiency vs N_{vertices};N_{vertices};Efficiency");

    // Plot individual efficiency graphs for medium ID+ISO
    plotSingleEfficiency(g_eff_mediumiso_pT, "D:/ROOTFiles/DataFiles/Eff_MediumIso_pT.png", 
       "Medium ID+ISO Efficiency vs pT;pT [GeV];Efficiency", true);
    plotSingleEfficiency(g_eff_mediumiso_eta, "D:/ROOTFiles/DataFiles/Eff_MediumIso_eta.png", 
       "Medium ID+ISO Efficiency vs #eta;#eta;Efficiency");
    plotSingleEfficiency(g_eff_mediumiso_phi, "D:/ROOTFiles/DataFiles/Eff_MediumIso_phi.png", 
       "Medium ID+ISO Efficiency vs #phi;#phi;Efficiency");
    plotSingleEfficiency(g_eff_mediumiso_nvtx, "D:/ROOTFiles/DataFiles/Eff_MediumIso_nvtx.png", 
       "Medium ID+ISO Efficiency vs N_{vertices};N_{vertices};Efficiency");

    // Plot individual efficiency graphs for loose ID+ISO
    plotSingleEfficiency(g_eff_looseiso_pT, "D:/ROOTFiles/DataFiles/Eff_LooseIso_pT.png", 
       "Loose ID+ISO Efficiency vs pT;pT [GeV];Efficiency", true);
    plotSingleEfficiency(g_eff_looseiso_eta, "D:/ROOTFiles/DataFiles/Eff_LooseIso_eta.png", 
       "Loose ID+ISO Efficiency vs #eta;#eta;Efficiency");
    plotSingleEfficiency(g_eff_looseiso_phi, "D:/ROOTFiles/DataFiles/Eff_LooseIso_phi.png", 
       "Loose ID+ISO Efficiency vs #phi;#phi;Efficiency");
    plotSingleEfficiency(g_eff_looseiso_nvtx, "D:/ROOTFiles/DataFiles/Eff_LooseIso_nvtx.png", 
       "Loose ID+ISO Efficiency vs N_{vertices};N_{vertices};Efficiency");
    
    // Plotting
    TCanvas *c = new TCanvas("c", "Canvas", 1200, 800);
    gStyle->SetOptStat(0);  // Turn off stat box for all histograms
    gStyle->SetErrorX(0.5);

    // Updated lambda with 3 parameters (3rd has default value)
    auto styleHistogram = [](TH1F* hist, Color_t color, bool isTotal = false) {
        hist->SetLineColor(color);
        hist->SetFillColorAlpha(color, isTotal ? 0.2 : 0.3);
        hist->SetLineWidth(isTotal ? 1 : 2);
        hist->SetMarkerStyle(isTotal ? 1 : 20);
        hist->SetMarkerColor(color);
        hist->SetMarkerSize(0.8);
        hist->SetStats(0);  // Turn off stats for this histogram
    };

    // Plot total probe histograms (now with isTotal=true)
    std::vector<std::pair<TH1F*, std::string>> total_hists = {
        {h_totalprobe_pT,   "D:/ROOTFiles/DataFiles/TotalProbe_pT.png"},
        {h_totalprobe_eta,  "D:/ROOTFiles/DataFiles/TotalProbe_eta.png"},
        {h_totalprobe_phi,  "D:/ROOTFiles/DataFiles/TotalProbe_phi.png"},
        {h_totalprobe_nvtx, "D:/ROOTFiles/DataFiles/TotalProbe_nvtx.png"},
        {h_totalprobe_Mll,  "D:/ROOTFiles/DataFiles/TotalProbe_Mll.png"}
    };

    for (auto &[hist, filename] : total_hists) {
        styleHistogram(hist, kRed, true); 
        hist->Draw("HIST");
        c->SaveAs(filename.c_str());
    }

    // plotting failed probe
    {
        TCanvas *c_failed = new TCanvas("c_failed", "Failed Probes", 800, 600);
        
        // Style the failed probe histogram differently
        h_failedprobe_Mll->SetLineColor(kRed);
        h_failedprobe_Mll->SetLineWidth(2);
        h_failedprobe_Mll->SetFillColorAlpha(kRed, 0.3);
        h_failedprobe_Mll->SetStats(0);
        
        h_failedprobe_Mll->Draw("HIST");
        
        c_failed->SaveAs("D:/ROOTFiles/DataFiles/FailedProbes_Mll.png");
        delete c_failed;
    }

    // Plot passing probe histograms (separate for each ID type)
    std::vector<std::tuple<TH1F*, std::string, Color_t>> pass_hists = {
        {h_passprobe_tight_pT,   "D:/ROOTFiles/DataFiles/PassProbe_Tight_pT.png", kBlue+3},
        {h_passprobe_tight_eta,  "D:/ROOTFiles/DataFiles/PassProbe_Tight_eta.png", kBlue+3},
        {h_passprobe_tight_phi,  "D:/ROOTFiles/DataFiles/PassProbe_Tight_phi.png", kBlue+3},
        {h_passprobe_tight_nvtx, "D:/ROOTFiles/DataFiles/PassProbe_Tight_nvtx.png", kBlue+3},
        {h_passprobe_tight_Mll,  "D:/ROOTFiles/DataFiles/PassProbe_Tight_Mll.png", kBlue+3},
        
        {h_passprobe_medium_pT,   "D:/ROOTFiles/DataFiles/PassProbe_Medium_pT.png", kRed+2},
        {h_passprobe_medium_eta,  "D:/ROOTFiles/DataFiles/PassProbe_Medium_eta.png", kRed+2},
        {h_passprobe_medium_phi,  "D:/ROOTFiles/DataFiles/PassProbe_Medium_phi.png", kRed+2},
        {h_passprobe_medium_nvtx, "D:/ROOTFiles/DataFiles/PassProbe_Medium_nvtx.png", kRed+2},
        {h_passprobe_medium_Mll,  "D:/ROOTFiles/DataFiles/PassProbe_Medium_Mll.png", kRed+2},
        
        {h_passprobe_loose_pT,   "D:/ROOTFiles/DataFiles/PassProbe_Loose_pT.png", kGreen+2},
        {h_passprobe_loose_eta,  "D:/ROOTFiles/DataFiles/PassProbe_Loose_eta.png", kGreen+2},
        {h_passprobe_loose_phi,  "D:/ROOTFiles/DataFiles/PassProbe_Loose_phi.png", kGreen+2},
        {h_passprobe_loose_nvtx, "D:/ROOTFiles/DataFiles/PassProbe_Loose_nvtx.png", kGreen+2},
        {h_passprobe_loose_Mll,  "D:/ROOTFiles/DataFiles/PassProbe_Loose_Mll.png", kGreen+2}
    };

    for (auto &[hist, filename, color] : pass_hists) {
        styleHistogram(hist, color);  // Using default isTotal=false
        hist->Draw("HIST");
        c->SaveAs(filename.c_str());
    }

    // Plot efficiency graphs (comparison of different ID types)
    auto plotEfficiencyComparison = [&](const std::vector<TGraphAsymmErrors*>& graphs, 
        const std::string& filename,
        const std::string& title) {
        TCanvas *c_eff = new TCanvas("c_eff", "Efficiency Canvas", 1200, 800);
        c_eff->SetGrid();

        // Draw the first graph to set up the axes
        graphs[0]->Draw("AP");
        graphs[0]->SetTitle(title.c_str());
        graphs[0]->GetYaxis()->SetRangeUser(0, 1.10);

        // Draw the other graphs
        for (size_t i = 1; i < graphs.size(); ++i) {
            graphs[i]->Draw("P SAME");
        }

        // Create legend in bottom right with box
        TLegend *leg = new TLegend(0.65, 0.15, 0.85, 0.35);  // Bottom right coordinates
        leg->SetBorderSize(1);  // Add border
        leg->SetFillColor(kWhite);  // White background
        leg->SetTextSize(0.035);

        // Customize legend entries - Only Medium and Loose now
        leg->AddEntry(g_eff_medium_pT, "Medium ID", "lp");
        leg->AddEntry(g_eff_loose_pT, "Loose ID", "lp");
        leg->Draw();

        c_eff->SaveAs(filename.c_str());
        delete c_eff;
    };

    // Plot ID efficiency comparisons (Medium vs Loose only)
    std::vector<std::tuple<std::vector<TGraphAsymmErrors*>, std::string, std::string>> eff_comparisons = {
        {{g_eff_medium_pT, g_eff_loose_pT}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompare_pT.png", "Muon ID Efficiency vs pT"},
         
        {{g_eff_medium_eta, g_eff_loose_eta}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompare_eta.png", "Muon ID Efficiency vs #eta"},
         
        {{g_eff_medium_phi, g_eff_loose_phi}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompare_phi.png", "Muon ID Efficiency vs #phi"},
         
        {{g_eff_medium_nvtx, g_eff_loose_nvtx}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompare_nvtx.png", "Muon ID Efficiency vs N_{vertices}"}
    };

    for (auto &[graphs, filename, title] : eff_comparisons) {
        plotEfficiencyComparison(graphs, filename, title);
    }

    // Plot isolation efficiency comparisons (Medium vs Loose only)
    std::vector<std::tuple<std::vector<TGraphAsymmErrors*>, std::string, std::string>> iso_eff_comparisons = {
        {{g_eff_mediumiso_pT, g_eff_looseiso_pT}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIso_pT.png", "Muon Isolation Efficiency vs pT (relative to ID)"},
         
        {{g_eff_mediumiso_eta, g_eff_looseiso_eta}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIso_eta.png", "Muon Isolation Efficiency vs #eta (relative to ID)"},
         
        {{g_eff_mediumiso_phi, g_eff_looseiso_phi}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIso_phi.png", "Muon Isolation Efficiency vs #phi (relative to ID)"},
         
        {{g_eff_mediumiso_nvtx, g_eff_looseiso_nvtx}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIso_nvtx.png", "Muon Isolation Efficiency vs N_{vertices} (relative to ID)"}
    };

    for (auto &[graphs, filename, title] : iso_eff_comparisons) {
        plotEfficiencyComparison(graphs, filename, title);
    }

    // Create ID+ISO efficiency graphs (relative to total probes)
    TGraphAsymmErrors *g_eff_tightIDISO_pT = new TGraphAsymmErrors(h_passprobe_tightiso_pT, h_totalprobe_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tightIDISO_eta = new TGraphAsymmErrors(h_passprobe_tightiso_eta, h_totalprobe_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tightIDISO_phi = new TGraphAsymmErrors(h_passprobe_tightiso_phi, h_totalprobe_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_tightIDISO_nvtx = new TGraphAsymmErrors(h_passprobe_tightiso_nvtx, h_totalprobe_nvtx, "cl=0.683 b(1,1) mode");

    TGraphAsymmErrors *g_eff_mediumIDISO_pT = new TGraphAsymmErrors(h_passprobe_mediumiso_pT, h_totalprobe_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_mediumIDISO_eta = new TGraphAsymmErrors(h_passprobe_mediumiso_eta, h_totalprobe_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_mediumIDISO_phi = new TGraphAsymmErrors(h_passprobe_mediumiso_phi, h_totalprobe_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_mediumIDISO_nvtx = new TGraphAsymmErrors(h_passprobe_mediumiso_nvtx, h_totalprobe_nvtx, "cl=0.683 b(1,1) mode");

    TGraphAsymmErrors *g_eff_looseIDISO_pT = new TGraphAsymmErrors(h_passprobe_looseiso_pT, h_totalprobe_pT, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_looseIDISO_eta = new TGraphAsymmErrors(h_passprobe_looseiso_eta, h_totalprobe_eta, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_looseIDISO_phi = new TGraphAsymmErrors(h_passprobe_looseiso_phi, h_totalprobe_phi, "cl=0.683 b(1,1) mode");
    TGraphAsymmErrors *g_eff_looseIDISO_nvtx = new TGraphAsymmErrors(h_passprobe_looseiso_nvtx, h_totalprobe_nvtx, "cl=0.683 b(1,1) mode");

    // Set graph properties for ID+ISO
    styleGraph(g_eff_tightIDISO_pT, kMagenta+1, 23, "Tight ID+ISO Efficiency");
    styleGraph(g_eff_tightIDISO_eta, kMagenta+1, 23, "Tight ID+ISO Efficiency");
    styleGraph(g_eff_tightIDISO_phi, kMagenta+1, 23, "Tight ID+ISO Efficiency");
    styleGraph(g_eff_tightIDISO_nvtx, kMagenta+1, 23, "Tight ID+ISO Efficiency");

    styleGraph(g_eff_mediumIDISO_pT, kOrange+7, 24, "Medium ID+ISO Efficiency");
    styleGraph(g_eff_mediumIDISO_eta, kOrange+7, 24, "Medium ID+ISO Efficiency");
    styleGraph(g_eff_mediumIDISO_phi, kOrange+7, 24, "Medium ID+ISO Efficiency");
    styleGraph(g_eff_mediumIDISO_nvtx, kOrange+7, 24, "Medium ID+ISO Efficiency");

    styleGraph(g_eff_looseIDISO_pT, kCyan+1, 25, "Loose ID+ISO Efficiency");
    styleGraph(g_eff_looseIDISO_eta, kCyan+1, 25, "Loose ID+ISO Efficiency");
    styleGraph(g_eff_looseIDISO_phi, kCyan+1, 25, "Loose ID+ISO Efficiency");
    styleGraph(g_eff_looseIDISO_nvtx, kCyan+1, 25, "Loose ID+ISO Efficiency");

    // Write ID+ISO efficiency graphs to file
    g_eff_tightIDISO_pT->Write("g_eff_tightIDISO_pT");
    g_eff_tightIDISO_eta->Write("g_eff_tightIDISO_eta");
    g_eff_tightIDISO_phi->Write("g_eff_tightIDISO_phi");
    g_eff_tightIDISO_nvtx->Write("g_eff_tightIDISO_nvtx");

    g_eff_mediumIDISO_pT->Write("g_eff_mediumIDISO_pT");
    g_eff_mediumIDISO_eta->Write("g_eff_mediumIDISO_eta");
    g_eff_mediumIDISO_phi->Write("g_eff_mediumIDISO_phi");
    g_eff_mediumIDISO_nvtx->Write("g_eff_mediumIDISO_nvtx");

    g_eff_looseIDISO_pT->Write("g_eff_looseIDISO_pT");
    g_eff_looseIDISO_eta->Write("g_eff_looseIDISO_eta");
    g_eff_looseIDISO_phi->Write("g_eff_looseIDISO_phi");
    g_eff_looseIDISO_nvtx->Write("g_eff_looseIDISO_nvtx");

    // Plot ID+ISO efficiency comparisons (relative to total probes)
    std::vector<std::tuple<std::vector<TGraphAsymmErrors*>, std::string, std::string>> idiso_eff_comparisons = {
        {{g_eff_mediumIDISO_pT, g_eff_looseIDISO_pT}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIDISO_pT.png", "Muon ID+ISO Efficiency vs pT"},
         
        {{g_eff_mediumIDISO_eta, g_eff_looseIDISO_eta}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIDISO_eta.png", "Muon ID+ISO Efficiency vs #eta"},
         
        {{g_eff_mediumIDISO_phi, g_eff_looseIDISO_phi}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIDISO_phi.png", "Muon ID+ISO Efficiency vs #phi"},
         
        {{g_eff_mediumIDISO_nvtx, g_eff_looseIDISO_nvtx}, 
         "D:/ROOTFiles/DataFiles/EfficiencyCompareIDISO_nvtx.png", "Muon ID+ISO Efficiency vs N_{vertices}"}
    };

    for (auto &[graphs, filename, title] : idiso_eff_comparisons) {
        plotEfficiencyComparison(graphs, filename, title);
    }
    // Draw and save 2D efficiency plots for each ID type
    auto plot2DEfficiency = [&](TH2F* pass_hist, TH2F* total_hist, const std::string& filename, const std::string& title) {
        // Clone and calculate efficiency
        TH2F* h_eff = (TH2F*)pass_hist->Clone();
        h_eff->Divide(total_hist);
        h_eff->SetTitle(title.c_str());
        h_eff->GetZaxis()->SetTitle("Efficiency");
        h_eff->SetMinimum(0);
        h_eff->SetMaximum(1);
        
        // Create canvas with larger size
        TCanvas *c2d = new TCanvas("c2d", "2D Efficiency", 1200, 900);
        gStyle->SetPaintTextFormat(".2f"); // Set 2 decimal precision for text
        gStyle->SetOptStat(0); // Remove statistics box
        
        // Draw with improved styling
        h_eff->Draw("COLZ"); // Color plot first
        h_eff->Draw("TEXT SAME"); // Overlay text
        
        // Improve axis labels
        h_eff->GetXaxis()->SetTitleSize(0.04);
        h_eff->GetYaxis()->SetTitleSize(0.04);
        h_eff->GetZaxis()->SetTitleSize(0.04);
        
        // Adjust text attributes
        h_eff->SetMarkerSize(1.2); // Larger text
        h_eff->SetMarkerColor(kBlack); // High contrast text
        
        // Add padding around the plot
        c2d->SetLeftMargin(0.12);
        c2d->SetRightMargin(0.14);
        c2d->SetBottomMargin(0.12);
        
        c2d->Update();
        c2d->SaveAs(filename.c_str());
        
        // Cleanup
        delete c2d;
        delete h_eff;
    };
    
    // Generate plots 
    plot2DEfficiency(h_passprobe_tight_pt_eta, h_totalprobe_pt_eta, 
        "D:/ROOTFiles/DataFiles/Efficiency2D_Tight.png", 
        "Tight ID Efficiency pT vs |#eta|;|#eta|;pT [GeV]");
    
    plot2DEfficiency(h_passprobe_medium_pt_eta, h_totalprobe_pt_eta, 
        "D:/ROOTFiles/DataFiles/Efficiency2D_Medium.png", 
        "Medium ID Efficiency pT vs |#eta|;|#eta|;pT [GeV]");
    
    plot2DEfficiency(h_passprobe_loose_pt_eta, h_totalprobe_pt_eta, 
        "D:/ROOTFiles/DataFiles/Efficiency2D_Loose.png", 
        "Loose ID Efficiency pT vs |#eta|;|#eta|;pT [GeV]");
    
    // Cleanup
    delete c;
    delete h_totalprobe_pt_eta;
    delete h_passprobe_tight_pt_eta;
    delete h_passprobe_medium_pt_eta;
    delete h_passprobe_loose_pt_eta;
    delete file;
}