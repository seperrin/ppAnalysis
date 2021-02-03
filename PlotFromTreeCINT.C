#include <TBrowser.h>
#include <TBufferFile.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TMath.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TFitResultPtr.h>
#include <TPaveText.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <string>

//#include "TTreeReader.h"
//#include "Event.h"

# include "TF1.h"
# include "TF2.h"
# include "TProfile.h"
# include "TVirtualFitter.h"
# include "TFitter.h"
# include "TProfile.h"
# include "TMinuit.h"
# include "TRandom.h"
# include <iostream>
# include <fstream>
# include <string>
# include "TH1F.h"
# include "TH1I.h"
# include "TH1D.h"
# include "TH2F.h"
# include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"


// FUNCTIONS

void PlotFromTreeCINT();


void PlotFromTreeCINT(){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    
// ************************************
// Définitions de paramètres          *
// ************************************

    double ZvtxCut = 10;
    double DPhiCut = 0.01;
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 99999;
    double MinMultCentral = 37;
    double MaxMultPeriph = 23;
    double CentSPDTrackletsCentral = 1;
    double CentSPDTrackletsPeriph = 7;

    double LowJPsiMass = 3.0;
    double HighJPsiMass = 3.3;

    double MinInvMass = 2.0;
    double MaxInvMass = 4.5;
    const int NbinsInvMass = 10;
    double SizeBinInvMass = (MaxInvMass-MinInvMass)/NbinsInvMass;

    double MinDeltaPhi = -TMath::Pi()/2;
    double MaxDeltaPhi = 1.5*TMath::Pi();
    const int NbinsDeltaPhi = 12;
    double SizeBinDeltaPhi = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhi;
    
    double MinDeltaEta = 0;
    double MaxDeltaEta = 6;
    const int NbinsDeltaEta = 20;
    double SizeBinDeltaEta = (MaxDeltaEta-MinDeltaEta)/NbinsDeltaEta;
    int barWidth = 50;
    
    Char_t Group_Period[50] = "Group12_LHC18m_CINT7CENT";
   // Char_t *arrayOfPeriods[] = {"Group5_LHC17l_CINT","Group5_LHC17m_CINT","Group5_LHC17o_CINT","Group5_LHC17r_CINT","Group5_LHC18c_CINT","Group5_LHC18d_CINT","Group5_LHC18e_CINT","Group5_LHC18f_CINT"};
    Char_t *arrayOfPeriods[] = {"Group12_LHC18m_CINT7CENT"};
    int numberOfPeriods = sizeof(arrayOfPeriods) / sizeof(arrayOfPeriods[0]);
    
    const double binsCent[6] = {0,1,10,20,40,100};
    Char_t filePMLim[200];
    Char_t filePMLimSaveName[200];
    Char_t filePM[200];
    Char_t fileNmean[200];
    Char_t fileNmeanROOT[200];
    Char_t fileInLoc[200];
// *************************
// Initialiser les graphes *
// *************************
    
    TH1F* hnseg(NULL);
    TH1I* hnseg_num(NULL);
    TH1* PMSliced[20]={NULL};
    TH2F* PercentileMethod(NULL);
    TH2F* PercentileMethodLimits(NULL);
    TH1F* hnseg_raw(NULL);
    TH1I* hnseg_raw_num(NULL);
    TH1I* hnseg_den(NULL);
    TH1F* hnseg_spdtkl(NULL);
    TH1I* hnseg_spdtkl_num(NULL);
    TH1F* hnseg_spdtklval(NULL);
    TH1I* hnseg_spdtklval_num(NULL);
    TH1I* hnseg_spdtklval_dist(NULL);
    TH1I* hnseg_Ntkl_dist(NULL);
    TProfile* hn(NULL);
    
    
    hnseg = new TH1F("hnseg",
                     "Mean Ntkl wrt Zvtx",
                     81,-10,10);
    hnseg->SetXTitle("Zvtx");
    hnseg->SetYTitle("Mean Ntkl");
    
    hnseg_num = new TH1I("hnseg_num",
                     "Mean Ntkl wrt Zvtx - Numerator",
                     81,-10,10);
    hnseg_num->SetXTitle("Zvtx");
    hnseg_num->SetYTitle("Sum Ntkl");
    
    PercentileMethod = new TH2F("PercentileMethod",
                     "Ntkl wrt Zvtx",
                     20,-10,10,100,0,100);
    PercentileMethod->SetXTitle("Zvtx");
    PercentileMethod->SetYTitle("Ntkl");
    
    PercentileMethodLimits = new TH2F("PercentileMethodLimits",
                     "Ntkl Limit wrt Zvtx and percentile",
                                      20,-10,10,5,binsCent);
    PercentileMethodLimits->SetXTitle("Zvtx");
    PercentileMethodLimits->SetYTitle("Percentile");
    
    hnseg_raw = new TH1F("hnseg_raw",
                     "Mean raw Ntkl wrt Zvtx",
                     81,-10,10);
    hnseg_raw->SetXTitle("Zvtx");
    hnseg_raw->SetYTitle("Mean raw Ntkl");
    
    hnseg_raw_num = new TH1I("hnseg_raw_num",
                     "Mean raw Ntkl wrt Zvtx - Numerator",
                     81,-10,10);
    hnseg_raw_num->SetXTitle("Zvtx");
    hnseg_raw_num->SetYTitle("Mean raw Ntkl");
    
    hnseg_den = new TH1I("hnseg_den",
                     "Zvtx distribution",
                     81,-10,10);
    hnseg_den->SetXTitle("Zvtx");
    hnseg_den->SetYTitle("Count");
    
    hnseg_spdtkl = new TH1F("hnseg_spdtkl",
                     "Mean SPDTkl wrt Zvtx",
                     81,-10,10);
    hnseg_spdtkl->SetXTitle("Zvtx");
    hnseg_spdtkl->SetYTitle("Mean SPDTkl");
    
    hnseg_spdtkl_num = new TH1I("hnseg_spdtkl_num",
                     "Mean SPDTkl wrt Zvtx - Numerator",
                     81,-10,10);
    hnseg_spdtkl_num->SetXTitle("Zvtx");
    hnseg_spdtkl_num->SetYTitle("Mean SPDTkl");
    
    hnseg_spdtklval = new TH1F("hnseg_spdtklval",
                     "Mean SPDTklVal wrt Zvtx",
                     81,-10,10);
    hnseg_spdtklval->SetXTitle("Zvtx");
    hnseg_spdtklval->SetYTitle("Mean SPDTklVal");
    
    hnseg_spdtklval_num = new TH1I("hnseg_spdtklval_num",
                        "Mean SPDTklVal wrt Zvtx - Numerator",
                        81,-10,10);
       hnseg_spdtklval_num->SetXTitle("Zvtx");
       hnseg_spdtklval_num->SetYTitle("Mean SPDTklVal");
    
    hnseg_spdtklval_dist = new TH1I("hnseg_spdtklval_dist",
                     "Distribution of SPDTkl value",
                     160,0,160);
    hnseg_spdtklval_dist->SetXTitle("SPDTkl value");
    hnseg_spdtklval_dist->SetYTitle("Count");
    
    hnseg_Ntkl_dist = new TH1I("hnseg_Ntkl_dist",
                     "Distribution of Ntkl value",
                     160,0,160);
    hnseg_Ntkl_dist->SetXTitle("Ntkl value");
    hnseg_Ntkl_dist->SetYTitle("Count");
    
    hn = new TProfile("hn",
                     "Mean SPDTklVal wrt Zvtx TPro",
                     81,-10,10);
    hn->SetXTitle("Zvtx");
    hn->SetYTitle("Mean SPDTklVal");

// *************************
// Analyse                 *
// *************************
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        sprintf(filePMLim,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/%s/%s_PMLim.pdf",Group_Period,Group_Period);
        sprintf(filePMLimSaveName,"PercentileMethodStudy/AccountNtklNonZero/%s/%s_PMLim.txt",Group_Period,Group_Period);
        sprintf(filePM,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/%s/%s_PM.pdf",Group_Period,Group_Period);
        sprintf(fileNmean,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/%s/%s_Nmean.pdf",Group_Period,Group_Period);
        sprintf(fileNmeanROOT,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/%s/%s_Nmean.root",Group_Period,Group_Period);
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        
    // TFile fileIn("~/../../Volumes/Transcend2/ppAnalysis/Scripts/muonGrid.root");
 //   TFile fileIn("~/../../Volumes/Transcend2/ppAnalysis/Scripts/Group12_LHC18m_CINT7CENT/muonGrid.root");
    TFile fileIn(fileInLoc);
   // TFile fileIn("~/../../Volumes/Transcend2/ppAnalysis/Scripts/CINT_16o_PS_CutsEvent_Val_CENT/muonGrid.root");
    std::cout << "1" <<std::endl;
        TTree* theTree = NULL;
        fileIn.GetObject("MyMuonTree",theTree);
    std::cout << "2" <<std::endl;
        
    MyEventLight* fEvent = 0;
    TClonesArray *fCorrelations = 0;
    TClonesArray *fTracklets = 0;
    TClonesArray *fDimuons = 0;
    CorrelationLight *correl = 0; //= new CorrelationLight(); //object must be created before
    TrackletLight *trac = 0;
    DimuonLight *dimu = 0;
        //setting the branch address
    std::cout << "3" <<std::endl;
  //  theTree->SetBranchAddress("event", &fEvent);
    TBranch *branch = theTree->GetBranch("event");
//        auto bntrack = theTree->GetBranch("tracks");
        //auto branch  = theTree->GetBranch("mcparticles.fPx");
    std::cout << "3.5" <<std::endl;
        branch->SetAddress(&fEvent);
    std::cout << "4" <<std::endl;
    auto nevent = theTree->GetEntries();
    
    std::cout << "5" <<std::endl;
    fCorrelations = fEvent->fCorrelations;
    fTracklets = fEvent->fTracklets;
    fDimuons = fEvent->fDimuons;
//            auto nevent = bntrack->GetEntries();
    //cout << " theTree->GetEntries() " << theTree->GetEntries() << endl;
    
    
// ************************************
// Boucle sur les evenements, notés i *
// ************************************
    
  //  Double_t px[1000], py[1000];
        for (Int_t i=0;i<nevent;i++) {
            if(i%100000 == 0){
                    std::cout << "[";
                    double portion = double(i)/nevent;
                    long pos = barWidth * portion;
                    for (int k = 0; k < barWidth; ++k) {
                        if (k < pos) std::cout << "=";
                        else if (k == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                    std::cout << "] " << long(100 * portion) << "%     " << i << "/" << nevent << " Tree " << tree_idx << "/" << numberOfPeriods;
                    std::cout.flush();
                std::cout << std::endl;
            }
            theTree->GetEvent(i);
     //       cout << "Looking at Event " << i << endl;
//            cout << "fPassPhysicsSelection : " << fEvent->fPassPhysicsSelection << endl;
//            cout << "fPassTriggerSelection : " << fEvent->fPassTriggerSelection << endl;
//            cout << "fTriggerCMUL7 : " << fEvent->fTriggerCMUL7 << endl;
//            cout << "fTriggerCINT7 : " << fEvent->fTriggerCINT7 << endl;
            //cout << "tracks->GetEntries() " << tracks->GetEntries() << endl;
          //  cout << "Computing the number of tracklets in Event " << i << endl;
            
          //  cout << "Reading Event " << i <<endl;
            // Apply a cut on the TrackletEta at +-1
            int NumberCloseEtaTracklets = 0;
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                if((TMath::Abs(trac->fEta) < 1)){
                        NumberCloseEtaTracklets++;
                    
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
                hnseg_spdtkl_num->Fill(fEvent->fVertexZ, fEvent->fCentralitySPDTracklets);
               hnseg_spdtklval_num->Fill(fEvent->fVertexZ, fEvent->fSPDTrackletsValue);
                hn->Fill(fEvent->fVertexZ, fEvent->fSPDTrackletsValue);
                hnseg_den->Fill(fEvent->fVertexZ);
                hnseg_raw_num->Fill(fEvent->fVertexZ, fEvent->fNTracklets);
                hnseg_num->Fill(fEvent->fVertexZ, NumberCloseEtaTracklets);
                PercentileMethod->Fill(fEvent->fVertexZ,NumberCloseEtaTracklets);
                hnseg_Ntkl_dist->Fill(NumberCloseEtaTracklets);
               }
            hnseg_spdtklval_dist->Fill(fEvent->fSPDTrackletsValue);
            
            // hnseg->Fill(fEvent->fVertexZ, NumberCloseEtaTracklets);
           //  hnseg_raw->Fill(fEvent->fVertexZ, fEvent->fNTracklets);
            
            
        }
    }
    
    TCanvas*cPMproj=new TCanvas();
    cPMproj->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMproj->cd(i);
        PMSliced[i] = PercentileMethod->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSliced[i]->Draw("e");
    }
    
    cout << "PMSliced[0] info: " << PMSliced[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSliced[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSliced[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimits->Fill((i-10)+0.5, 0.5);
                if(ratio<0.90){
                    PercentileMethodLimits->Fill((i-10)+0.5, 5.5);
                    if(ratio<0.80){
                        PercentileMethodLimits->Fill((i-10)+0.5, 15);
                        if(ratio<0.60){
                            PercentileMethodLimits->Fill((i-10)+0.5, 30);
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimSave(filePMLimSaveName, std::ofstream::out);
    
    filePMLimSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimSave << "{";
        for(int cent_idx=1; cent_idx<6; cent_idx++){
            filePMLimSave << PercentileMethodLimits->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=5){
                filePMLimSave << " ";
            }
        }
        if(z_idx!=20){
            filePMLimSave << "},\n";
        }
        else{
            filePMLimSave << "}   }\n";
        }
    }
    cout <<filePMLimSave.is_open() <<endl;
    
    filePMLimSave.close();

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    
    
    
// ***********************************
// Tracer les graphes sur les canevas *
// ***********************************


    TCanvas*c1=new TCanvas();
    c1->Divide(2,1);
    c1->cd(1);
    hnseg->Divide(hnseg_num, hnseg_den);
    hnseg->SetFillColor(18);
    hnseg->Draw();
//    c1->cd(2);
//    hnseg_raw->Divide(hnseg_den);
//    hnseg_raw->SetFillColor(18);
//    hnseg_raw->Draw();
    c1->cd(2);
    hnseg_den->SetFillColor(18);
    hnseg_den->Draw();
    
    TCanvas*c1solo=new TCanvas();
        hnseg->Draw();
    c1solo->SaveAs(fileNmean);
    
    TFile *outputFileNmean;
        outputFileNmean = new TFile(fileNmeanROOT,"RECREATE");
        hnseg->Write();
    
    TCanvas*c2=new TCanvas();
    hnseg_raw->Divide(hnseg_raw_num, hnseg_den);
    hnseg_raw->SetFillColor(18);
    hnseg_raw->Draw();
    
    TCanvas*c3=new TCanvas();
    hnseg_spdtkl->Divide(hnseg_spdtkl_num, hnseg_den);
    hnseg_spdtkl->SetFillColor(18);
    hnseg_spdtkl->Draw();
    
    TCanvas*c4=new TCanvas();
    hnseg_spdtklval->Divide(hnseg_spdtklval_num, hnseg_den);
    hnseg_spdtklval->SetFillColor(18);
    hnseg_spdtklval->Draw();
    
    TCanvas*c4Pro=new TCanvas();
    hn->SetFillColor(18);
    hn->Draw();
    
    TCanvas*c45=new TCanvas();
    hnseg_spdtklval_dist->SetFillColor(18);
    hnseg_spdtklval_dist->Draw("e");
    
    TCanvas*c456=new TCanvas();
    hnseg_Ntkl_dist->SetFillColor(18);
    hnseg_Ntkl_dist->Draw("e");
    
    TCanvas*cPM=new TCanvas();
    PercentileMethod->SetFillColor(18);
    PercentileMethod->Draw("colz");
    cPM->SaveAs(filePM);
    
    TCanvas*cPMLim=new TCanvas();
    PercentileMethodLimits->SetFillColor(18);
    PercentileMethodLimits->Draw("colz");
    cPMLim->SaveAs(filePMLim);
    
    
}
    
