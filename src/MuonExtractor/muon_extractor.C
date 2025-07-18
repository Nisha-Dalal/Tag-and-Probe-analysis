#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <vector>
#include <string>

/*
    This code extracts muon-related branches along with useful event metadata,
    MET information, primary vertex count, and gen-level information
    from multiple input ROOT files.
    The selected branches are saved into "Muons.root".
*/

void muon_extractor() {
    std::vector<std::string> branches_to_keep = {
        // Event-level metadata
        "run", "event", "luminosityBlock",

        // Muon kinematics and ID
        "nMuon",
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
        "Muon_looseId", "Muon_mediumId", "Muon_tightId",
        "Muon_pfRelIso03_all", "Muon_pfRelIso04_all",
        "Muon_dxy", "Muon_dxyErr", "Muon_dz", "Muon_dzErr",
        "Muon_ip3d", "Muon_sip3d", "Muon_genPartIdx",
        "Muon_highPtId", "Muon_inTimeMuon", "Muon_isGlobal",
        "Muon_isPFcand", "Muon_isTracker", "Muon_mediumPromptId",
        "Muon_miniIsoId", "Muon_multiIsoId", "Muon_mvaId",
        "Muon_mvaLowPt", "Muon_mvaTTH",
        "Muon_softId", "Muon_softMva", "Muon_tkRelIso",

        // Event-level variables
        "PV_npvs",
        "MET_phi", "PuppiMET_phi",
        "MET_pt", "PuppiMET_pt",

        // Gen-particle info
        "GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_status", "GenPart_statusFlags"
    };
    std::vector<std::string> input_files = {
        "D:/ROOTFiles/DataFiles/1B11D47A-B6B6-2444-9A83-1649F6B54F3A.root",
        "D:/ROOTFiles/DataFiles/14B6A8AE-C9FE-D744-80A4-DDE5D008C1CD.root",
        "D:/ROOTFiles/DataFiles/202FEE1A-6266-ED4D-A0B2-F3D0C5B6EE1B.root",
        "D:/ROOTFiles/DataFiles/216C2A60-D4EC-B647-8421-BAB0C6A4247C.root",
        "D:/ROOTFiles/DataFiles/1542B5DE-398C-0A4E-A970-9CD38709B98F.root",
        "D:/ROOTFiles/DataFiles/1544A357-D3B8-B84C-A854-C38E93BCC67A.root",
        "D:/ROOTFiles/DataFiles/6C40BF66-B8A7-534D-BA55-03A9EE560117.root",
        "D:/ROOTFiles/DataFiles/53A2A4DA-60DD-1E48-86E7-466B4B2A78FC.root",
        "D:/ROOTFiles/DataFiles/4B9EE242-BFEC-D14F-9D4A-B6682A4E1B7C.root",
        "D:/ROOTFiles/DataFiles/3411EB3A-03C1-6A4F-A0E7-570107EEFD0C.root",
        "D:/ROOTFiles/DataFiles/488285BF-DA04-914A-8D07-72E5342923DD.root",
        "D:/ROOTFiles/DataFiles/22766655-2217-B848-8979-AD0EA0128682.root",
        "D:/ROOTFiles/DataFiles/59183214-61E6-A64C-AEDF-8EB180E8CFC6.root",
      
    };
    // Load the ROOT files into a TChain
    TChain *chain = new TChain("Events");
    for (const auto& file : input_files) {
        chain->Add(file.c_str());
    }

    if (chain->GetNtrees() == 0 || chain->GetEntries() == 0) {
        std::cerr << "Error: No valid input files or entries found!" << std::endl;
        delete chain;
        return;
    }

    Long64_t total_events = chain->GetEntries();
    std::cout << "Total number of events: " << total_events << std::endl;

    // Disable all branches, then enable only selected ones
    chain->SetBranchStatus("*", 0);
    for (const auto& branch : branches_to_keep) {
        chain->SetBranchStatus(branch.c_str(), 1);
    }

    // Optional: set cache size for speedup
    chain->SetCacheSize(512 * 1024);

    // Create output file
    TFile *output_file = new TFile("D:/ROOTFiles/DataFiles/Muons.root", "RECREATE");
    if (!output_file || output_file->IsZombie()) {
    std::cerr << "Error: Could not create Muons.root at D:/ROOTFiles/DataFiles/" << std::endl;
    delete chain;
    return;
}


    output_file->cd();

    // Create a new tree to save passed events
    TTree *output_tree = chain->CloneTree(0);  // Empty tree structure with same branches

    // Define variables for basic selection
    UInt_t nMuon = 0;
    chain->SetBranchAddress("nMuon", &nMuon);

    // Count passing events
    Long64_t passed_events = 0;

    for (Long64_t i = 0; i < total_events; ++i) {
        chain->GetEntry(i);

        // Apply basic cut: at least one muon
        if (nMuon > 0) {
            output_tree->Fill();
            ++passed_events;
        }

        // Optional: progress bar
        if (i % 500000 == 0) {
            std::cout << "Processed " << i << " / " << total_events << " events...\r" << std::flush;
        }
    }

    std::cout << "\n Number of events passed basic cuts: " << passed_events << std::endl;

    output_tree->SetBasketSize("*", 8000);
    output_tree->Write("", TObject::kOverwrite);
    output_file->Close();

    delete chain;
    delete output_file;

    std::cout << "\n Muon branches successfully extracted to Muons.root\n";
}
