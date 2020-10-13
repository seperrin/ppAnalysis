#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "AliAnalysisTaskMyMuonTree_AOD.h"

//int GetNumberOfFiles(TString path, TString name) {
//    /*
//     returns the number of fils that exist in a folder (the path) containing the Tstring name
//     */
//    Int_t num = 0;
//    struct dirent **namelist;
//    Int_t n = scandir(path, &namelist, 0, alphasort);
//    if (n < 1) {cout << "empty folder" << endl; return 0;}
//    else {
//        while (n--) { if (strstr(namelist[n]->d_name, name) != NULL) num++;}
//        return num;
//    }
//}

int MergeRoots() {
    
    gROOT->ProcessLine("hadd merged.root GoodLHC16o/muonGrid.root GoodLHC16p/muonGrid.root");
    
//    TString path = "~/../../Volumes/Transcend2/ppAnalysis/Scripts/";
//    TString periods[] = {"16o", "16p"};
//
//    TString filename = "GoodLHC";
//    TString outfilename = "MuonMergedGoodLHC";
//
//    int numberOfFiles = sizeof(periods)/sizeof(TString); //GetNumberOfFiles(path, filename);
//    std::cout << "number of files = " << numberOfFiles << std::endl;
//    if (numberOfFiles == 0) {std::cout << "no files to add, terminating" << std::endl; return 0;}
//
//    TString outputName = path + outfilename+".root";
//    TFile* fOut = new TFile(outputName, "RECREATE");
//
//    //Initialisation du Tree Mergé
//
//    // Initialisation of the new TTree tAvalanche
//    MyEventLight* event = 0;
//    TClonesArray *fCorrelations = 0;
//    TClonesArray *fTracklets = 0;
//    TClonesArray *fDimuons = 0;
//    MyEventLight* eventInput = 0;
//    TTree *tree = new TTree("MyMuonTree","Muon tree");
//    tree->Branch("event", event);
//
//
//    // On commence par remplir le TTree tAvalanche
//    for (int i = 0; i<numberOfFiles; i++) {
//        //for (int i = 2; i<3; i++) {
//        TString inputName = path + filename + Form("%s/muonGrid.root", periods[i].Data());
//        std::cout << "including " << inputName << std::endl;
//        TFile* fIn = TFile::Open(inputName, "READ");
//        TTree* theTree = (TTree*)fIn->Get("MyMuonTree");
//
//
//
//        theTree->SetBranchAddress("event", &eventInput);
//
//        int nIn = theTree->GetEntries();
//
//        for (int l = 0; l<nIn; l++) {
//            theTree->GetEntry(l);
//
//            event = eventInput;
//            fCorrelations = event->fCorrelations;
//            fTracklets = event->fTracklets;
//            fDimuons = event->fDimuons;
//
//            tree->Fill();
//        }
//
//    }
//    fOut->cd();
//    tree->Write();
//
//    fOut->Close();
    
    return 0;
    
}
