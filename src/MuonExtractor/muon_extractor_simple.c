#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>

void muon_extractor_simple() {
    TChain *chain = new TChain("Events");
    chain->Add("D:/ROOTFiles/DataFiles/14B6A8AE-C9FE-D744-80A4-DDE5D008C1CD.root");

    if (chain->GetNtrees() == 0 || chain->GetEntries() == 0) {
        std::cerr << "Error: No valid input files or entries found!" << std::endl;
        delete chain;
        return;
    }

    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("run", 1);
    chain->SetBranchStatus("event", 1);
    chain->SetBranchStatus("nMuon", 1);
    chain->SetBranchStatus("Muon_pt", 1);

    TFile *output_file = new TFile("Muons_simple.root", "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        std::cerr << "Error: Could not create Muons_simple.root" << std::endl;
        delete chain;
        return;
    }

    output_file->cd();
    TTree *output_tree = chain->CloneTree(-1, "fast");
    if (!output_tree) {
        std::cerr << "Error: Failed to clone tree" << std::endl;
        delete chain;
        delete output_file;
        return;
    }

    output_tree->Write("", TObject::kOverwrite);
    output_file->Close();

    delete chain;
    delete output_file;

    std::cout << "\nMuon branches successfully extracted to Muons_simple.root\n";
}
