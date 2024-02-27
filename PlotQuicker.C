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


#include "TAxis.h"
#include "TLegend.h"
#include "TH1.h"
#include "TArrayD.h"
# include "TF1.h"
# include "TF2.h"
# include "TProfile.h"
# include "TVirtualFitter.h"
# include "TGraphPolar.h"
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

//Code pour faire des plots rapides d'observables
// FUNCTIONS

void PlotQuicker();
void PlotQuickest();
void PlotQuickie();
void PlotQuickie2();
void PlotComparison(); //V0M SPDtracklets
void PlotSaucisson();
void PlotMCdistrib();
void PlotCorrelEstimators();
void PlotPolar();
void CorrelEstMult();
void PlotPHD();
void CheckPYTHIA();
Double_t myFunc(double x);
Double_t myFunc2(double x);
Double_t myEtaLow(double x);
Double_t myEtaHigh(double x);
TH2I* hPlotDimuEtaPhi(NULL);
TH2I* hPlotTklEtaPhi(NULL);

void PlotQuicker(){
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
    
    Char_t Group_Period[50] = "Group1_LHC16h";
   // Char_t *arrayOfPeriods[] = {"Group5_LHC17l_CINT","Group5_LHC17m_CINT","Group5_LHC17o_CINT","Group5_LHC17r_CINT","Group5_LHC18c_CINT","Group5_LHC18d_CINT","Group5_LHC18e_CINT","Group5_LHC18f_CINT"};
    //Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
    Char_t *arrayOfPeriods[] = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l","Group2_LHC17h","Group3_LHC17h","Group4_LHC17k","Group4_LHC18l","Group4_LHC18m","Group4_LHC18o","Group4_LHC18p","Group5_LHC17l","Group5_LHC17m","Group5_LHC17o","Group5_LHC17r","Group5_LHC18c","Group5_LHC18d","Group5_LHC18e","Group5_LHC18f","Group6_LHC18m","Group7_LHC18m","Group8_LHC18m","Group9_LHC18m","Group10_LHC18m","Group11_LHC18m","Group12_LHC18m"};
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
    TH2F* hPlotXY(NULL);
    TH2F* hEtaWrtZ(NULL);
    TProfile* hn(NULL);
    TH2I* hPlotTklTklDeltaEtaDeltaPhi(NULL);
    TH2I* hPlotDimuTklDeltaEtaDeltaPhi(NULL);
    
    
    hPlotXY = new TH2F("hPlotXY",
                        "X and Y of vertex SPD",
                        500,-0.5,0.5,500,-0.5,0.5);
       hPlotXY->SetXTitle("X");
       hPlotXY->SetYTitle("Y");
    
    hPlotDimuEtaPhi = new TH2I("hPlotDimuEtaPhi",
                     "EtaPhi Distribution for dimuons",
                     720,0,TMath::Pi()*2,100,-6,0);
    hPlotDimuEtaPhi->SetXTitle("Phi");
    hPlotDimuEtaPhi->SetYTitle("Eta");
    
    hPlotTklEtaPhi = new TH2I("hPlotTklEtaPhi",
                        "Eta Phi Distribution for tracklets",
                        720,0,TMath::Pi()*2,100,-1,1);
       hPlotTklEtaPhi->SetXTitle("Phi");
       hPlotTklEtaPhi->SetYTitle("Eta");
    
    hPlotTklTklDeltaEtaDeltaPhi = new TH2I("hPlotTklTklDeltaEtaDeltaPhi",
                     "DeltaEta DeltaPhi Distribution for TklTkl",
                     72,-TMath::Pi()/2,TMath::Pi()*1.5,100,0,2);
    hPlotTklTklDeltaEtaDeltaPhi->SetXTitle("DeltaPhi");
    hPlotTklTklDeltaEtaDeltaPhi->SetYTitle("DeltaEta");
    
    hPlotDimuTklDeltaEtaDeltaPhi = new TH2I("hPlotDimuTklDeltaEtaDeltaPhi",
                     "DeltaEta DeltaPhi Distribution for DimuTkl",
                     72,-TMath::Pi()/2,TMath::Pi()*1.5,100,0,6);
    hPlotDimuTklDeltaEtaDeltaPhi->SetXTitle("DeltaPhi");
    hPlotDimuTklDeltaEtaDeltaPhi->SetYTitle("DeltaEta");
    
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
    
    hEtaWrtZ = new TH2F("hEtaWrtZ",
                      "Tracklet eta wrt z_vertex",
                      100,-15,15,100,-2.5,2.5);
    hEtaWrtZ->SetXTitle("Z_{vertex} (cm)");
    hEtaWrtZ->SetYTitle("Tracklet Eta");

// *************************
// Analyse                 *
// *************************
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
        
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
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                hEtaWrtZ->Fill(fEvent->fVertexZ, trac->fEta);
            }
            
            continue;
            
          //  cout << "Reading Event " << i <<endl;
            // Apply a cut on the TrackletEta at +-1
            int NumberCloseEtaTracklets = 0;
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                    if(fEvent->fVertexZ>0 && fEvent->fVertexZ<1){
                        hPlotTklEtaPhi->Fill(trac->fPhi,trac->fEta);
                    }
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
                hnseg_spdtkl_num->Fill(fEvent->fVertexZ, fEvent->fCentralitySPDTracklets);
               hnseg_spdtklval_num->Fill(fEvent->fVertexZ, fEvent->fSPDTrackletsValue);
                hn->Fill(fEvent->fVertexZ, fEvent->fSPDTrackletsValue);
                hnseg_den->Fill(fEvent->fVertexZ);
                hnseg_raw_num->Fill(fEvent->fVertexZ, fEvent->fNTracklets);
                hnseg_num->Fill(fEvent->fVertexZ, NumberCloseEtaTracklets);
                PercentileMethod->Fill(fEvent->fVertexZ,NumberCloseEtaTracklets);
                hnseg_Ntkl_dist->Fill(NumberCloseEtaTracklets);
                
                
                if(fEvent->fVertexZ>0 && fEvent->fVertexZ<1){
                    for (Int_t j=0; j<fEvent->fNDimuons; j++) {
                    dimu = (DimuonLight*)fDimuons->At(j);
                        if((dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut) && (dimu->fInvMass > MinInvMass) && (dimu->fInvMass < MaxInvMass)){
                            hPlotDimuEtaPhi->Fill(dimu->fPhi,dimu->fEta);
                        }
                    }
                    
                  //  hPlotXY->Fill(fEvent->fVertexXSPD,fEvent->fVertexYSPD);
                }
                
                
                
                
               }
            hnseg_spdtklval_dist->Fill(fEvent->fSPDTrackletsValue);
            
            // hnseg->Fill(fEvent->fVertexZ, NumberCloseEtaTracklets);
           //  hnseg_raw->Fill(fEvent->fVertexZ, fEvent->fNTracklets);
            
            
        }
        
        
    }
    
    TCanvas* cDimuEtaPhi=new TCanvas();
    hPlotDimuEtaPhi->Draw("colz");
    
    TCanvas* cTklEtaPhi=new TCanvas();
    hPlotTklEtaPhi->Draw("colz");
    
    TCanvas* cXY = new TCanvas();
    hPlotXY->Draw("colz");
    
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
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    
    double tkl1phi = 0;
    double tkl1eta = 0;
    double tkl2phi = 0;
    double tkl2eta = 0;
    double tkl3phi = 0;
    double tkl3eta = 0;
    double dimuonphi = 0;
    double dimuoneta = 0;
    
    TF2 *fPhiEtaTkl = new TF2("fPhiEtaTkl","myFunc(x,y)",0,TMath::Pi()*2,-1,1);
    TF2 *fPhiEtaDimu = new TF2("fPhiEtaDimu","myFunc2(x,y)",0,TMath::Pi()*2,-6,0);
    
    for(int essai=0; essai<1000000;essai++){
        if(essai%1000000 == 0){
            cout << essai << "/1000000000"<<endl;
        }
        fPhiEtaTkl->GetRandom2(tkl1phi, tkl1eta);
        fPhiEtaTkl->GetRandom2(tkl2phi, tkl2eta);
        fPhiEtaTkl->GetRandom2(tkl3phi, tkl3eta);
        fPhiEtaDimu->GetRandom2(dimuonphi, dimuoneta);
        
        
        double DeltaPhiTklTkl = tkl1phi-tkl2phi;
        double DeltaPhiDimuTkl = tkl3phi-dimuonphi;
        double DeltaEtaTklTkl = TMath::Abs(tkl1eta-tkl2eta);
        double DeltaEtaDimuTkl = tkl3eta-dimuoneta;
        
        if(DeltaPhiTklTkl < -TMath::Pi()/2){
            DeltaPhiTklTkl += 2* TMath::Pi();
        }
        if(DeltaPhiTklTkl > 1.5*TMath::Pi()){
            DeltaPhiTklTkl -= 2* TMath::Pi();
        }
        if(DeltaPhiDimuTkl < -TMath::Pi()/2){
            DeltaPhiDimuTkl += 2* TMath::Pi();
        }
        if(DeltaPhiDimuTkl > 1.5*TMath::Pi()){
            DeltaPhiDimuTkl -= 2* TMath::Pi();
        }
        
        if(DeltaEtaTklTkl>1.2){
            hPlotTklTklDeltaEtaDeltaPhi->Fill(DeltaPhiTklTkl,DeltaEtaTklTkl);
        }
        hPlotDimuTklDeltaEtaDeltaPhi->Fill(DeltaPhiDimuTkl,DeltaEtaDimuTkl);
    }
    
    TCanvas*cDeltaEtaDeltaPhiTklTkl=new TCanvas();
    hPlotTklTklDeltaEtaDeltaPhi->Draw("colz");
    
    TCanvas*cDeltaPhiDeltaEtaDimuTkl=new TCanvas();
    hPlotDimuTklDeltaEtaDeltaPhi->Draw("colz");
    
    
    
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
    
    TCanvas*cPMLim=new TCanvas();
    PercentileMethodLimits->SetFillColor(18);
    PercentileMethodLimits->Draw("colz");
    
    
    TF1 *fEtaLow = new TF1("fEtaLow","myEtaLow(x)",-10,10);
    TF1 *fEtaHigh = new TF1("fEtaHigh","myEtaHigh(x)",-10,10);
    
    TLine *l=new TLine(-10,-2.5,-10,2.5);
    l->SetLineColor(kBlack);
    l->SetLineWidth(3);
    l->SetLineStyle(9);
    
    TLine *l2=new TLine(10,-2.5,10,2.5);
    l2->SetLineColor(kBlack);
    l2->SetLineWidth(3);
    l2->SetLineStyle(9);
    
    TCanvas*cAA=new TCanvas();
    hEtaWrtZ->Draw("colz");
    l->Draw("same");
    l2->Draw("same");
    fEtaLow->SetLineWidth(3);
    fEtaLow->Draw("same");
    fEtaHigh->SetLineWidth(3);
    fEtaHigh->Draw("same");
    
    
}

void PlotQuickest(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
    Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
   // Char_t *arrayOfPeriods[] = {"Group5_LHC17l_CINT","Group5_LHC17m_CINT","Group5_LHC17o_CINT","Group5_LHC17r_CINT","Group5_LHC18c_CINT","Group5_LHC18d_CINT","Group5_LHC18e_CINT","Group5_LHC18f_CINT"};
    Char_t *arrayOfPeriods[] = {"Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose"};//, "Group1_LHC17k", "Group4_LHC17k", "Group4_LHC18o"};
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
    
    TH2F* V0M_vs_SPD(NULL);
    TH1F* NewEst_norm(NULL);
    TH1F* NewEst_norm_pond(NULL);
    TH1F* V0M_Dist(NULL);
    TH1F* V0M_Dist_Groups(NULL);
    TH1F* V0M_Dist_Groups_norm(NULL);
    TH1F* SPDTracklets_Dist(NULL);
    TH1F* SPDTracklets_Dist_Groups(NULL);
    TH1F* SPDTracklets_Dist_Groups_norm(NULL);
    TH2F* V0M_vs_SPD_norm(NULL);
    TH1F* V0M_Dist_norm(NULL);
    TH1F* SPDTracklets_Dist_norm(NULL);
    
    
    
    V0M_Dist_Groups = new TH1F("V0M_Dist_Groups",
                     "V0M_Dist_Groups",
                     1000,-0.5, 999.5);
    V0M_Dist_Groups->SetXTitle("V0M");
    V0M_Dist_Groups->SetYTitle("Count");
    
    V0M_Dist_Groups_norm = new TH1F("V0M_Dist_Groups_norm",
                     "V0M_Dist_Groups_norm",
                     1001,-0.05, 100.05);
    V0M_Dist_Groups_norm->SetXTitle("V0M_norm");
    V0M_Dist_Groups_norm->SetYTitle("Count");
    
    SPDTracklets_Dist_Groups = new TH1F("SPDTracklets_Dist_Groups",
                     "SPDTracklets_Dist_Groups",
                     1001,-0.05, 100.05);
    SPDTracklets_Dist_Groups->SetXTitle("SPDTracklets");
    SPDTracklets_Dist_Groups->SetYTitle("Count");
    
    SPDTracklets_Dist_Groups_norm = new TH1F("SPDTracklets_Dist_Groups_norm",
                     "SPDTracklets_Dist_Groups_norm",
                     1001,-0.05, 100.05);
    SPDTracklets_Dist_Groups_norm->SetXTitle("SPDTracklets_norm");
    SPDTracklets_Dist_Groups_norm->SetYTitle("Count");
    
   
// *************************
// Analyse                 *
// *************************
    
    TCanvas* cV0M_Dist_Groups=new TCanvas("V0M Distribution Groups");
    TCanvas* cSPDTracklets_Dist_Groups=new TCanvas("SPDTracklets Distribution Groups");
    TCanvas* cV0M_Dist_Groups_norm=new TCanvas("V0M norm Distribution Groups");
    TCanvas* cSPDTracklets_Dist_Groups_norm=new TCanvas("SPDTracklets norm Distribution Groups");
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        
        V0M_vs_SPD = new TH2F("V0M_vs_SPD",
                            "V0M_vs_SPD",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD->SetXTitle("V0M");
           V0M_vs_SPD->SetYTitle("SPDTracklets");
        
        V0M_Dist = new TH1F("V0M_Dist",
                         "V0M_Dist",
                         1000,-0.5, 999.5);
        V0M_Dist->SetXTitle("V0M");
        V0M_Dist->SetYTitle("Count");
        
        SPDTracklets_Dist = new TH1F("SPDTracklets_Dist",
                         "SPDTracklets_Dist",
                         1000,-0.5, 999.5);
        SPDTracklets_Dist->SetXTitle("SPDTracklets");
        SPDTracklets_Dist->SetYTitle("Count");
        
        V0M_vs_SPD_norm = new TH2F("V0M_vs_SPD_norm",
                            "V0M_vs_SPD_norm",
                            1001,-0.05, 100.05,1001,-0.05, 100.05);
           V0M_vs_SPD_norm->SetXTitle("V0M_norm");
           V0M_vs_SPD_norm->SetYTitle("SPDTracklets_norm");
        
        NewEst_norm = new TH1F("NewEst_norm",
                         "NewEst_norm",
                         1000,0, 100);
        NewEst_norm->SetXTitle("NewEst_norm");
        NewEst_norm->SetYTitle("Count");
        
        NewEst_norm_pond = new TH1F("NewEst_norm_pond",
                         "NewEst_norm_pond",
                         1000,0, 100);
        NewEst_norm_pond->SetXTitle("NewEst_norm_pond");
        NewEst_norm_pond->SetYTitle("Count");
        
        V0M_Dist_norm = new TH1F("V0M_Dist_norm",
                         "V0M_Dist_norm",
                         1000,0, 100);
        V0M_Dist_norm->SetXTitle("V0M_norm");
        V0M_Dist_norm->SetYTitle("Count");
        
        SPDTracklets_Dist_norm = new TH1F("SPDTracklets_Dist_norm",
                         "SPDTracklets_Dist_norm",
                         100,-0.5, 99.5);
        SPDTracklets_Dist_norm->SetXTitle("SPDTracklets_norm");
        SPDTracklets_Dist_norm->SetYTitle("Count");
        
        
        
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
            
            if(i%10 != 0){
                continue;
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
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
                V0M_Dist->Fill(fEvent->fV0MValue);
                SPDTracklets_Dist->Fill(fEvent->fSPDTrackletsValue);
                V0M_vs_SPD->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                V0M_Dist_norm->Fill((fEvent->fV0MValue/medianV0M16h)*10);
                SPDTracklets_Dist_norm->Fill((fEvent->fSPDTrackletsValue/medianSPD16h)*10,0.1);
                V0M_vs_SPD_norm->Fill((fEvent->fV0MValue/medianV0M16h)*10,((fEvent->fSPDTrackletsValue/medianSPD16h)*10));
                
                //V0M_Dist_norm->Scale(1./entries16h);
                //SPDTracklets_Dist_norm->Scale(1./10);
                NewEst_norm->Fill(((4.3*((fEvent->fV0MValue/medianV0M16h)*10)+(2.*((fEvent->fSPDTrackletsValue/medianSPD16h)*10)))/6.3));
                NewEst_norm_pond->Fill(((4.3*(5.39/5.62)*((fEvent->fV0MValue/medianV0M16h)*10)+(2.*((fEvent->fSPDTrackletsValue/medianSPD16h)*10)))/6.3));

               }
           
        }
        
        //Getting quantiles
        int cumul = 0;
        double q_90 = 0;
        double q_80 = 0;
        double q_70 = 0;
        double q_60 = 0;
        double q_50 = 0;
        double q_40 = 0;
        double q_30 = 0;
        double q_20 = 0;
        double q_15 = 0;
        double q_10 = 0;
        double q_05 = 0;
        double q_03 = 0;
        double q_01 = 0;
 
        entries16h = NewEst_norm->GetEntries();
                
        for(int bin_idx=NewEst_norm->GetNbinsX(); bin_idx>=1; bin_idx--){
            if(cumul <= entries16h*0.9){
                q_90 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                if(cumul <= entries16h*0.8){
                    q_80 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                    if(cumul <= entries16h*0.7){
                        q_70 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                        if(cumul <= entries16h*0.6){
                            q_60 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                            if(cumul <= entries16h*0.5){
                                q_50 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                if(cumul <= entries16h*0.4){
                                    q_40 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                    if(cumul <= entries16h*0.3){
                                        q_30 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                        if(cumul <= entries16h*0.2){
                                            q_20 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                            if(cumul <= entries16h*0.15){
                                                q_15 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                                if(cumul <= entries16h*0.1){
                                                    q_10 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                                    if(cumul <= entries16h*0.05){
                                                        q_05 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                                        if(cumul <= entries16h*0.03){
                                                            q_03 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                                            if(cumul <= entries16h*0.01){
                                                                q_01 = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            cumul+=NewEst_norm->GetBinContent(bin_idx);
        }
        
        double maxNewEst = 0;
                    for(int bin_idx=1; bin_idx<=NewEst_norm->GetNbinsX(); bin_idx++){
                        if(NewEst_norm->GetBinContent(bin_idx)>0){
                            maxNewEst = NewEst_norm->GetXaxis()->GetBinCenter(bin_idx);
                        }
                    }
        
        cout << "0-1% is above "<< q_01<<endl;
        cout << "1-3% is above "<< q_03<<endl;
        cout << "3-5% is above "<< q_05<<endl;
        cout << "5-10% is above "<< q_10<<endl;
        cout << "10-15% is above "<< q_15<<endl;
        cout << "15-20% is above "<< q_20<<endl;
        cout << "20-30% is above "<< q_30<<endl;
        cout << "30-40% is above "<< q_40<<endl;
        cout << "40-50% is above "<< q_50<<endl;
        cout << "50-60% is above "<< q_60<<endl;
        cout << "60-70% is above "<< q_70<<endl;
        cout << "70-80% is above "<< q_80<<endl;
        cout << "80-90% is above "<< q_90<<endl;
        cout << "90-100% is the rest from 0 to "<<q_90<<endl;
        
        
        
        //Getting quantiles for Average = Density
        int cumul_pond = 0;
        double q_90_pond = 0;
        double q_80_pond = 0;
        double q_70_pond = 0;
        double q_60_pond = 0;
        double q_50_pond = 0;
        double q_40_pond = 0;
        double q_30_pond = 0;
        double q_20_pond = 0;
        double q_15_pond = 0;
        double q_10_pond = 0;
        double q_05_pond = 0;
        double q_03_pond = 0;
        double q_01_pond = 0;
        
        for(int bin_idx=NewEst_norm_pond->GetNbinsX(); bin_idx>=1; bin_idx--){
            if(cumul_pond <= entries16h*0.9){
                q_90_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                if(cumul_pond <= entries16h*0.8){
                    q_80_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                    if(cumul_pond <= entries16h*0.7){
                        q_70_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                        if(cumul_pond <= entries16h*0.6){
                            q_60_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                            if(cumul_pond <= entries16h*0.5){
                                q_50_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                if(cumul_pond <= entries16h*0.4){
                                    q_40_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                    if(cumul_pond <= entries16h*0.3){
                                        q_30_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                        if(cumul_pond <= entries16h*0.2){
                                            q_20_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                            if(cumul_pond <= entries16h*0.15){
                                                q_15_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                                if(cumul_pond <= entries16h*0.1){
                                                    q_10_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                                    if(cumul_pond <= entries16h*0.05){
                                                        q_05_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                                        if(cumul_pond <= entries16h*0.03){
                                                            q_03_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                                            if(cumul_pond <= entries16h*0.01){
                                                                q_01_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            cumul_pond+=NewEst_norm_pond->GetBinContent(bin_idx);
        }
        
        double maxNewEst_pond = 0;
                    for(int bin_idx=1; bin_idx<=NewEst_norm_pond->GetNbinsX(); bin_idx++){
                        if(NewEst_norm_pond->GetBinContent(bin_idx)>0){
                            maxNewEst_pond = NewEst_norm_pond->GetXaxis()->GetBinCenter(bin_idx);
                        }
                    }
        
        cout << "Takind density into account"<<endl;
        cout << "0-1% is above "<< q_01_pond<<endl;
        cout << "1-3% is above "<< q_03_pond<<endl;
        cout << "3-5% is above "<< q_05_pond<<endl;
        cout << "5-10% is above "<< q_10_pond<<endl;
        cout << "10-15% is above "<< q_15_pond<<endl;
        cout << "15-20% is above "<< q_20_pond<<endl;
        cout << "20-30% is above "<< q_30_pond<<endl;
        cout << "30-40% is above "<< q_40_pond<<endl;
        cout << "40-50% is above "<< q_50_pond<<endl;
        cout << "50-60% is above "<< q_60_pond<<endl;
        cout << "60-70% is above "<< q_70_pond<<endl;
        cout << "70-80% is above "<< q_80_pond<<endl;
        cout << "80-90% is above "<< q_90_pond<<endl;
        cout << "90-100% is the rest from 0 to "<<q_90_pond<<endl;
        
        

            //Getting distributions infos

            //Finding the medians
            int entriesV0M = V0M_Dist->GetEntries();
            int cumulativeV0M = 0;
            double medianV0M = 0;
            for(int bin_idx=1; bin_idx<=V0M_Dist->GetNbinsX(); bin_idx++){
                if(cumulativeV0M < entriesV0M/2){
                    cumulativeV0M+= V0M_Dist->GetBinContent(bin_idx);
                    if(cumulativeV0M >= entriesV0M/2){
                        medianV0M = V0M_Dist->GetXaxis()->GetBinCenter(bin_idx);
                    }
                }
            }

            int entriesSPD = SPDTracklets_Dist->GetEntries();
            int cumulativeSPD = 0;
            double medianSPD = 0;
            for(int bin_idx=1; bin_idx<=SPDTracklets_Dist->GetNbinsX(); bin_idx++){
                if(cumulativeSPD < entriesSPD/2){
                    cumulativeSPD+= SPDTracklets_Dist->GetBinContent(bin_idx);
                    if(cumulativeSPD >= entriesSPD/2){
                        medianSPD = SPDTracklets_Dist->GetXaxis()->GetBinCenter(bin_idx);
                    }
                }
            }

            //Finding minimum and maximum

            double maxV0M = 0;
            for(int bin_idx=1; bin_idx<=V0M_Dist->GetNbinsX(); bin_idx++){
                if(V0M_Dist->GetBinContent(bin_idx)>0){
                    maxV0M = V0M_Dist->GetXaxis()->GetBinCenter(bin_idx);
                }
            }

            double maxSPD = 0;
            for(int bin_idx=1; bin_idx<=SPDTracklets_Dist->GetNbinsX(); bin_idx++){
                if(SPDTracklets_Dist->GetBinContent(bin_idx)>0){
                    maxSPD = SPDTracklets_Dist->GetXaxis()->GetBinCenter(bin_idx);
                }
            }



        cout << "Entries V0M: " << entriesV0M<<endl;
        cout << "Maximum V0M: " << maxV0M<<endl;
        cout << "Median V0M: " << medianV0M<<endl;
        cout << "Entries SPD: " << entriesSPD<<endl;
        cout << "Maximum SPD: " << maxSPD<<endl;
        cout << "Median SPD: " << medianSPD<<endl;
        
        char hname[100];
        sprintf(hname, "V0M_Dist - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M=new TCanvas(hname);
        V0M_Dist->SetTitle(hname);
        V0M_Dist->Draw();
        
        sprintf(hname, "SPDTracklets_Dist - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cSPDTracklets=new TCanvas(hname);
        SPDTracklets_Dist->SetTitle(hname);
        SPDTracklets_Dist->Draw();
        
        sprintf(hname, "V0M_vs_SPD - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD = new TCanvas(hname);
        V0M_vs_SPD->SetTitle(hname);
        V0M_vs_SPD->Draw("colz");
        
        sprintf(hname, "V0M_Dist_norm - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_norm=new TCanvas(hname);
        V0M_Dist_norm->SetTitle(hname);
        V0M_Dist_norm->Draw();
        
        sprintf(hname, "SPDTracklets_Dist_norm - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cSPDTracklets_norm=new TCanvas(hname);
        SPDTracklets_Dist_norm->SetTitle(hname);
        SPDTracklets_Dist_norm->Draw();
        
        sprintf(hname, "NewEst_norm - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cNewEst_norm=new TCanvas(hname);
        cNewEst_norm->cd();
        NewEst_norm->SetTitle(hname);
        NewEst_norm->SetLineColor(kGreen);
        NewEst_norm->SetLineWidth(3);
        NewEst_norm->Draw();
        V0M_Dist_norm->SetLineColor(kBlue);
        V0M_Dist_norm->SetLineWidth(3);
        V0M_Dist_norm->Draw("same");
        SPDTracklets_Dist_norm->SetLineColor(kRed);
        SPDTracklets_Dist_norm->SetLineWidth(3);
        SPDTracklets_Dist_norm->Draw("same");
        NewEst_norm_pond->SetLineColor(kMagenta);
        NewEst_norm_pond->SetLineWidth(3);
        NewEst_norm_pond->Draw("same");
        
        sprintf(hname, "V0M_vs_SPD_norm - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_norm = new TCanvas(hname);
        V0M_vs_SPD_norm->SetTitle(hname);
        V0M_vs_SPD_norm->Draw("colz");
        
        char message[80];
        sprintf(message, "%s", arrayOfPeriods[tree_idx]);
        legend->AddEntry(V0M_Dist,message);
        
        cV0M_Dist_Groups->cd();
      //  V0M_Dist->SetLineColor(tree_idx+1);
        V0M_Dist->Draw("same");
        
        cV0M_Dist_Groups_norm->cd();
       // V0M_Dist_norm->SetLineColor(tree_idx+1);
        V0M_Dist_norm->Draw("same");
        
        cSPDTracklets_Dist_Groups->cd();
       // SPDTracklets_Dist->SetLineColor(tree_idx+1);
        SPDTracklets_Dist->Draw("same");
        
        cSPDTracklets_Dist_Groups_norm->cd();
        //SPDTracklets_Dist_norm->SetLineColor(tree_idx+1);
        SPDTracklets_Dist_norm->Draw("same");
        
        
    }
    
    cV0M_Dist_Groups->cd();
    legend->Draw();
    cV0M_Dist_Groups_norm->cd();
    legend->Draw();
    cSPDTracklets_Dist_Groups->cd();
    legend->Draw();
    cSPDTracklets_Dist_Groups_norm->cd();
    legend->Draw();

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}

void PlotComparison(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;

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
    
    Char_t Group_Period[50] = "Group1_LHC16h";
   // Char_t *arrayOfPeriods[] = {"Group5_LHC17l_CINT","Group5_LHC17m_CINT","Group5_LHC17o_CINT","Group5_LHC17r_CINT","Group5_LHC18c_CINT","Group5_LHC18d_CINT","Group5_LHC18e_CINT","Group5_LHC18f_CINT"};
    Char_t *arrayOfPeriods[] = {"Group1_LHC16h", "Group1_LHC16o"};//, "Group1_LHC17k", "Group4_LHC17k", "Group4_LHC18o"};
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
    
    TH2F* V0M_vs_SPD_All(NULL);
    TH2F* Tracklets_All(NULL);
    TH1F* Tracklets_Dist_All(NULL);
    TH1F* Z_All(NULL);
    TH2F* V0M_vs_SPD_BothCentral(NULL);
    TH2F* Tracklets_BothCentral(NULL);
    TH1F* Tracklets_Dist_BothCentral(NULL);
    TH1F* Z_BothCentral(NULL);
    TH2F* V0M_vs_SPD_BothPeripheral(NULL);
    TH2F* Tracklets_BothPeripheral(NULL);
    TH1F* Tracklets_Dist_BothPeripheral(NULL);
    TH1F* Z_BothPeripheral(NULL);
    TH2F* V0M_vs_SPD_OnlyV0MCentral(NULL);
    TH2F* Tracklets_OnlyV0MCentral(NULL);
    TH1F* Tracklets_Dist_OnlyV0MCentral(NULL);
    TH1F* Z_OnlyV0MCentral(NULL);
    TH2F* V0M_vs_SPD_OnlyV0MPeripheral(NULL);
    TH2F* Tracklets_OnlyV0MPeripheral(NULL);
    TH1F* Tracklets_Dist_OnlyV0MPeripheral(NULL);
    TH1F* Z_OnlyV0MPeripheral(NULL);
    TH2F* V0M_vs_SPD_OnlySPDCentral(NULL);
    TH2F* Tracklets_OnlySPDCentral(NULL);
    TH1F* Tracklets_Dist_OnlySPDCentral(NULL);
    TH1F* Z_OnlySPDCentral(NULL);
    TH2F* V0M_vs_SPD_OnlySPDPeripheral(NULL);
    TH2F* Tracklets_OnlySPDPeripheral(NULL);
    TH1F* Tracklets_Dist_OnlySPDPeripheral(NULL);
    TH1F* Z_OnlySPDPeripheral(NULL);
    
    TH2F* V0M_vs_SPD_V0MCentral(NULL);
    TH2F* Tracklets_V0MCentral(NULL);
    TH1F* Tracklets_Dist_V0MCentral(NULL);
    TH1F* Tracklets_Dist_V0MCentral_Dimu(NULL);
    TH1F* Z_V0MCentral(NULL);
    TH2F* V0M_vs_SPD_V0MPeripheral(NULL);
    TH2F* Tracklets_V0MPeripheral(NULL);
    TH1F* Tracklets_Dist_V0MPeripheral(NULL);
    TH1F* Tracklets_Dist_V0MPeripheral_Dimu(NULL);
    TH1F* Z_V0MPeripheral(NULL);
    TH2F* V0M_vs_SPD_SPDCentral(NULL);
    TH2F* Tracklets_SPDCentral(NULL);
    TH1F* Tracklets_Dist_SPDCentral(NULL);
    TH1F* Tracklets_Dist_SPDCentral_Dimu(NULL);
    TH1F* Z_SPDCentral(NULL);
    TH2F* V0M_vs_SPD_SPDPeripheral(NULL);
    TH2F* Tracklets_SPDPeripheral(NULL);
    TH1F* Tracklets_Dist_SPDPeripheral(NULL);
    TH1F* Tracklets_Dist_SPDPeripheral_Dimu(NULL);
    TH1F* Z_SPDPeripheral(NULL);
    
   
// *************************
// Analyse                 *
// *************************
    
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
       // sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
        
        V0M_vs_SPD_All = new TH2F("V0M_vs_SPD_All",
                            "V0M_vs_SPD_All",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD_All->SetXTitle("V0M");
           V0M_vs_SPD_All->SetYTitle("SPDTracklets");
        
        Z_All = new TH1F("Z_All",
                         "Z_All",
                         200,-10,10);
        Z_All->SetXTitle("V0M");
        Z_All->SetYTitle("SPDTracklets");
        
        Tracklets_All = new TH2F("Tracklets_All",
                         "Tracklets_All",
                         120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
        Tracklets_All->SetXTitle("Phi");
        Tracklets_All->SetYTitle("Eta");
        
        Tracklets_Dist_All = new TH1F("Tracklets_Dist_All",
                         "Tracklets_Dist_All",
                         100,-0.5,99.5);
        Tracklets_Dist_All->SetXTitle("Tkl");
        
        V0M_vs_SPD_BothCentral = new TH2F("V0M_vs_SPD_BothCentral",
                            "V0M_vs_SPD_BothCentral",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD_BothCentral->SetXTitle("V0M");
           V0M_vs_SPD_BothCentral->SetYTitle("SPDTracklets");
        
        Z_BothCentral = new TH1F("Z_BothCentral",
                         "Z_BothCentral",
                         200,-10,10);
        Z_BothCentral->SetXTitle("V0M");
        Z_BothCentral->SetYTitle("SPDTracklets");
        
        Tracklets_BothCentral = new TH2F("Tracklets_BothCentral",
                         "Tracklets_BothCentral",
                         120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
        Tracklets_BothCentral->SetXTitle("Phi");
        Tracklets_BothCentral->SetYTitle("Eta");
        
        Tracklets_Dist_BothCentral = new TH1F("Tracklets_Dist_BothCentral",
                         "Tracklets_Dist_BothCentral",
                         100,-0.5,99.5);
        Tracklets_Dist_BothCentral->SetXTitle("Tkl");
        
        V0M_vs_SPD_BothPeripheral = new TH2F("V0M_vs_SPD_BothPeripheral",
                            "V0M_vs_SPD_BothPeripheral",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD_BothPeripheral->SetXTitle("V0M");
           V0M_vs_SPD_BothPeripheral->SetYTitle("SPDTracklets");
        
        Z_BothPeripheral = new TH1F("Z_BothPeripheral",
                         "Z_BothPeripheral",
                         200,-10,10);
        Z_BothPeripheral->SetXTitle("V0M");
        Z_BothPeripheral->SetYTitle("SPDTracklets");
        
        Tracklets_BothPeripheral = new TH2F("Tracklets_BothPeripheral",
                         "Tracklets_BothPeripheral",
                         120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
        Tracklets_BothPeripheral->SetXTitle("Phi");
        Tracklets_BothPeripheral->SetYTitle("Eta");
        
        Tracklets_Dist_BothPeripheral = new TH1F("Tracklets_Dist_BothPeripheral",
                         "Tracklets_Dist_BothPeripheral",
                         100,-0.5,99.5);
        Tracklets_Dist_BothPeripheral->SetXTitle("Tkl");
        
        V0M_vs_SPD_OnlyV0MCentral = new TH2F("V0M_vs_SPD_OnlyV0MCentral",
                            "V0M_vs_SPD_OnlyV0MCentral",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD_OnlyV0MCentral->SetXTitle("V0M");
           V0M_vs_SPD_OnlyV0MCentral->SetYTitle("SPDTracklets");
        
        Z_OnlyV0MCentral = new TH1F("Z_OnlyV0MCentral",
                         "Z_OnlyV0MCentral",
                         200,-10,10);
        Z_OnlyV0MCentral->SetXTitle("V0M");
        Z_OnlyV0MCentral->SetYTitle("SPDTracklets");
        
        Tracklets_OnlyV0MCentral = new TH2F("Tracklets_OnlyV0MCentral",
                         "Tracklets_OnlyV0MCentral",
                         120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
        Tracklets_OnlyV0MCentral->SetXTitle("Phi");
        Tracklets_OnlyV0MCentral->SetYTitle("Eta");
        
        Tracklets_Dist_OnlyV0MCentral = new TH1F("Tracklets_Dist_OnlyV0MCentral",
                         "Tracklets_Dist_OnlyV0MCentral",
                         100,-0.5,99.5);
        Tracklets_Dist_OnlyV0MCentral->SetXTitle("Tkl");
        
        V0M_vs_SPD_OnlyV0MPeripheral = new TH2F("V0M_vs_SPD_OnlyV0MPeripheral",
                            "V0M_vs_SPD_OnlyV0MPeripheral",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD_OnlyV0MPeripheral->SetXTitle("V0M");
           V0M_vs_SPD_OnlyV0MPeripheral->SetYTitle("SPDTracklets");
        
        Z_OnlyV0MPeripheral = new TH1F("Z_OnlyV0MPeripheral",
                         "Z_OnlyV0MPeripheral",
                         200,-10,10);
        Z_OnlyV0MPeripheral->SetXTitle("V0M");
        Z_OnlyV0MPeripheral->SetYTitle("SPDTracklets");
        
        Tracklets_OnlyV0MPeripheral = new TH2F("Tracklets_OnlyV0MPeripheral",
                         "Tracklets_OnlyV0MPeripheral",
                         120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
        Tracklets_OnlyV0MPeripheral->SetXTitle("Phi");
        Tracklets_OnlyV0MPeripheral->SetYTitle("Eta");
        
        Tracklets_Dist_OnlyV0MPeripheral = new TH1F("Tracklets_Dist_OnlyV0MPeripheral",
                         "Tracklets_Dist_OnlyV0MPeripheral",
                         100,-0.5,99.5);
        Tracklets_Dist_OnlyV0MPeripheral->SetXTitle("Tkl");
        
        V0M_vs_SPD_OnlySPDCentral = new TH2F("V0M_vs_SPD_OnlySPDCentral",
                            "V0M_vs_SPD_OnlySPDCentral",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD_OnlySPDCentral->SetXTitle("V0M");
           V0M_vs_SPD_OnlySPDCentral->SetYTitle("SPDTracklets");
        
        Z_OnlySPDCentral = new TH1F("Z_OnlySPDCentral",
                         "Z_OnlySPDCentral",
                         200,-10,10);
        Z_OnlySPDCentral->SetXTitle("V0M");
        Z_OnlySPDCentral->SetYTitle("SPDTracklets");
        
        Tracklets_OnlySPDCentral = new TH2F("Tracklets_OnlySPDCentral",
                         "Tracklets_OnlySPDCentral",
                         120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
        Tracklets_OnlySPDCentral->SetXTitle("Phi");
        Tracklets_OnlySPDCentral->SetYTitle("Eta");
        
        Tracklets_Dist_OnlySPDCentral = new TH1F("Tracklets_Dist_OnlySPDCentral",
                         "Tracklets_Dist_OnlySPDCentral",
                         100,-0.5,99.5);
        Tracklets_Dist_OnlySPDCentral->SetXTitle("Tkl");
        
        V0M_vs_SPD_OnlySPDPeripheral = new TH2F("V0M_vs_SPD_OnlySPDPeripheral",
                            "V0M_vs_SPD_OnlySPDPeripheral",
                            1000,-0.5,999.5,150,-0.5,149.5);
           V0M_vs_SPD_OnlySPDPeripheral->SetXTitle("V0M");
           V0M_vs_SPD_OnlySPDPeripheral->SetYTitle("SPDTracklets");
        
        Z_OnlySPDPeripheral = new TH1F("Z_OnlySPDPeripheral",
                         "Z_OnlySPDPeripheral",
                         200,-10,10);
        Z_OnlySPDPeripheral->SetXTitle("V0M");
        Z_OnlySPDPeripheral->SetYTitle("SPDTracklets");
        
        Tracklets_OnlySPDPeripheral = new TH2F("Tracklets_OnlySPDPeripheral",
                         "Tracklets_OnlySPDPeripheral",
                         120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
        Tracklets_OnlySPDPeripheral->SetXTitle("Phi");
        Tracklets_OnlySPDPeripheral->SetYTitle("Eta");
        
        Tracklets_Dist_OnlySPDPeripheral = new TH1F("Tracklets_Dist_OnlySPDPeripheral",
                         "Tracklets_Dist_OnlySPDPeripheral",
                         100,-0.5,99.5);
        Tracklets_Dist_OnlySPDPeripheral->SetXTitle("Tkl");
        
        
        V0M_vs_SPD_V0MCentral = new TH2F("V0M_vs_SPD_V0MCentral",
                                   "V0M_vs_SPD_V0MCentral",
                                   1000,-0.5,999.5,150,-0.5,149.5);
                  V0M_vs_SPD_V0MCentral->SetXTitle("V0M");
                  V0M_vs_SPD_V0MCentral->SetYTitle("SPDTracklets");
               
               Z_V0MCentral = new TH1F("Z_V0MCentral",
                                "Z_V0MCentral",
                                200,-10,10);
               Z_V0MCentral->SetXTitle("V0M");
               Z_V0MCentral->SetYTitle("SPDTracklets");
               
               Tracklets_V0MCentral = new TH2F("Tracklets_V0MCentral",
                                "Tracklets_V0MCentral",
                                120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
               Tracklets_V0MCentral->SetXTitle("Phi");
               Tracklets_V0MCentral->SetYTitle("Eta");
               
               Tracklets_Dist_V0MCentral = new TH1F("Tracklets_Dist_V0MCentral",
                                "Tracklets_Dist_V0MCentral",
                                100,-0.5,99.5);
               Tracklets_Dist_V0MCentral->SetXTitle("Tkl");
        
        Tracklets_Dist_V0MCentral_Dimu = new TH1F("Tracklets_Dist_V0MCentral_Dimu",
                         "Tracklets_Dist_V0MCentral_Dimu",
                         100,-0.5,99.5);
        Tracklets_Dist_V0MCentral_Dimu->SetXTitle("Tkl");
               
               V0M_vs_SPD_V0MPeripheral = new TH2F("V0M_vs_SPD_V0MPeripheral",
                                   "V0M_vs_SPD_V0MPeripheral",
                                   1000,-0.5,999.5,150,-0.5,149.5);
                  V0M_vs_SPD_V0MPeripheral->SetXTitle("V0M");
                  V0M_vs_SPD_V0MPeripheral->SetYTitle("SPDTracklets");
               
               Z_V0MPeripheral = new TH1F("Z_V0MPeripheral",
                                "Z_V0MPeripheral",
                                200,-10,10);
               Z_V0MPeripheral->SetXTitle("V0M");
               Z_V0MPeripheral->SetYTitle("SPDTracklets");
               
               Tracklets_V0MPeripheral = new TH2F("Tracklets_V0MPeripheral",
                                "Tracklets_V0MPeripheral",
                                120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
               Tracklets_V0MPeripheral->SetXTitle("Phi");
               Tracklets_V0MPeripheral->SetYTitle("Eta");
               
               Tracklets_Dist_V0MPeripheral = new TH1F("Tracklets_Dist_V0MPeripheral",
                                "Tracklets_Dist_V0MPeripheral",
                                100,-0.5,99.5);
               Tracklets_Dist_V0MPeripheral->SetXTitle("Tkl");
        
        Tracklets_Dist_V0MPeripheral_Dimu = new TH1F("Tracklets_Dist_V0MPeripheral_Dimu",
                         "Tracklets_Dist_V0MPeripheral_Dimu",
                         100,-0.5,99.5);
        Tracklets_Dist_V0MPeripheral_Dimu->SetXTitle("Tkl");
               
               V0M_vs_SPD_SPDCentral = new TH2F("V0M_vs_SPD_SPDCentral",
                                   "V0M_vs_SPD_SPDCentral",
                                   1000,-0.5,999.5,150,-0.5,149.5);
                  V0M_vs_SPD_SPDCentral->SetXTitle("V0M");
                  V0M_vs_SPD_SPDCentral->SetYTitle("SPDTracklets");
               
               Z_SPDCentral = new TH1F("Z_SPDCentral",
                                "Z_SPDCentral",
                                200,-10,10);
               Z_SPDCentral->SetXTitle("V0M");
               Z_SPDCentral->SetYTitle("SPDTracklets");
               
               Tracklets_SPDCentral = new TH2F("Tracklets_SPDCentral",
                                "Tracklets_SPDCentral",
                                120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
               Tracklets_SPDCentral->SetXTitle("Phi");
               Tracklets_SPDCentral->SetYTitle("Eta");
               
               Tracklets_Dist_SPDCentral = new TH1F("Tracklets_Dist_SPDCentral",
                                "Tracklets_Dist_SPDCentral",
                                100,-0.5,99.5);
               Tracklets_Dist_SPDCentral->SetXTitle("Tkl");
        
        Tracklets_Dist_SPDCentral_Dimu = new TH1F("Tracklets_Dist_SPDCentral_Dimu",
                         "Tracklets_Dist_SPDCentral_Dimu",
                         100,-0.5,99.5);
        Tracklets_Dist_SPDCentral_Dimu->SetXTitle("Tkl");
               
               V0M_vs_SPD_SPDPeripheral = new TH2F("V0M_vs_SPD_SPDPeripheral",
                                   "V0M_vs_SPD_SPDPeripheral",
                                   1000,-0.5,999.5,150,-0.5,149.5);
                  V0M_vs_SPD_SPDPeripheral->SetXTitle("V0M");
                  V0M_vs_SPD_SPDPeripheral->SetYTitle("SPDTracklets");
               
               Z_SPDPeripheral = new TH1F("Z_SPDPeripheral",
                                "Z_SPDPeripheral",
                                200,-10,10);
               Z_SPDPeripheral->SetXTitle("V0M");
               Z_SPDPeripheral->SetYTitle("SPDTracklets");
               
               Tracklets_SPDPeripheral = new TH2F("Tracklets_SPDPeripheral",
                                "Tracklets_SPDPeripheral",
                                120,MinDeltaPhi,MaxDeltaPhi,100,-6,6);
               Tracklets_SPDPeripheral->SetXTitle("Phi");
               Tracklets_SPDPeripheral->SetYTitle("Eta");
               
               Tracklets_Dist_SPDPeripheral = new TH1F("Tracklets_Dist_SPDPeripheral",
                                "Tracklets_Dist_SPDPeripheral",
                                100,-0.5,99.5);
               Tracklets_Dist_SPDPeripheral->SetXTitle("Tkl");
        
        Tracklets_Dist_SPDPeripheral_Dimu = new TH1F("Tracklets_Dist_SPDPeripheral_Dimu",
                         "Tracklets_Dist_SPDPeripheral_Dimu",
                         100,-0.5,99.5);
        Tracklets_Dist_SPDPeripheral_Dimu->SetXTitle("Tkl");
        
        
        
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
            
//            if(i%100 != 0){
//                continue;
//            }
            theTree->GetEvent(i);
            
            bool isV0MPeripheral = kFALSE;
            bool isSPDPeripheral = kFALSE;
            bool isV0MCentral = kFALSE;
            bool isSPDCentral = kFALSE;
            
            if(fEvent->fCentralityV0M > 40){
                isV0MPeripheral = kTRUE;
            }
            if(fEvent->fCentralitySPDTracklets > 40){
                isSPDPeripheral = kTRUE;
            }
            if(fEvent->fCentralityV0M < 5){
                isV0MCentral = kTRUE;
            }
            if(fEvent->fCentralitySPDTracklets < 5){
                isSPDCentral = kTRUE;
            }
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
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
                
                
                // === ALL EVENTS ===
                
                V0M_vs_SPD_All->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                    trac = (TrackletLight*)fTracklets->At(j);
                    if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                            Tracklets_All->Fill(trac->fPhi, trac->fEta);
                    }
                    
                }
                Tracklets_Dist_All->Fill(NumberCloseEtaTracklets);
                Z_All->Fill(fEvent->fVertexZ);
                
                // === EVENTS BOTH CENTRAL ===
                
                if(isV0MCentral == kTRUE && isSPDCentral == kTRUE){
                    V0M_vs_SPD_BothCentral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                    for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                        trac = (TrackletLight*)fTracklets->At(j);
                        if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                Tracklets_BothCentral->Fill(trac->fPhi, trac->fEta);
                        }
                        
                    }
                    Tracklets_Dist_BothCentral->Fill(NumberCloseEtaTracklets);
                    Z_BothCentral->Fill(fEvent->fVertexZ);
                }
                
                // == EVENTS BOTH PERIPHERAL ===
                
                if(isV0MPeripheral == kTRUE && isSPDPeripheral == kTRUE){
                    V0M_vs_SPD_BothPeripheral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                    for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                        trac = (TrackletLight*)fTracklets->At(j);
                        if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                Tracklets_BothPeripheral->Fill(trac->fPhi, trac->fEta);
                        }
                        
                    }
                    Tracklets_Dist_BothPeripheral->Fill(NumberCloseEtaTracklets);
                    Z_BothPeripheral->Fill(fEvent->fVertexZ);
                }
                
                // === V0M CENTRAL not SPD ===
                
                if(isV0MCentral == kTRUE && isSPDCentral == kFALSE){
                    V0M_vs_SPD_OnlyV0MCentral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                    for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                        trac = (TrackletLight*)fTracklets->At(j);
                        if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                Tracklets_OnlyV0MCentral->Fill(trac->fPhi, trac->fEta);
                        }
                        
                    }
                    Tracklets_Dist_OnlyV0MCentral->Fill(NumberCloseEtaTracklets);
                    Z_OnlyV0MCentral->Fill(fEvent->fVertexZ);
                }
                
                if(isV0MCentral == kTRUE){
                                   V0M_vs_SPD_V0MCentral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                                   for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                                       trac = (TrackletLight*)fTracklets->At(j);
                                       if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                               Tracklets_V0MCentral->Fill(trac->fPhi, trac->fEta);
                                       }
                                       
                                   }
                                   Tracklets_Dist_V0MCentral->Fill(NumberCloseEtaTracklets);
                                   Z_V0MCentral->Fill(fEvent->fVertexZ);
                                    if(fEvent->fNDimuons>0){
                                        Tracklets_Dist_V0MCentral_Dimu->Fill(NumberCloseEtaTracklets);
                                    }
                               }
                
                // === V0M PERIPH not SPD ===
                
                if(isV0MPeripheral == kTRUE && isSPDPeripheral == kFALSE){
                    V0M_vs_SPD_OnlyV0MPeripheral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                    for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                        trac = (TrackletLight*)fTracklets->At(j);
                        if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                Tracklets_OnlyV0MPeripheral->Fill(trac->fPhi, trac->fEta);
                        }
                        
                    }
                    Tracklets_Dist_OnlyV0MPeripheral->Fill(NumberCloseEtaTracklets);
                    Z_OnlyV0MPeripheral->Fill(fEvent->fVertexZ);
                }
                
                if(isV0MPeripheral == kTRUE){
                    V0M_vs_SPD_V0MPeripheral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                    for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                        trac = (TrackletLight*)fTracklets->At(j);
                        if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                Tracklets_V0MPeripheral->Fill(trac->fPhi, trac->fEta);
                        }
                        
                    }
                    Tracklets_Dist_V0MPeripheral->Fill(NumberCloseEtaTracklets);
                    Z_V0MPeripheral->Fill(fEvent->fVertexZ);
                    if(fEvent->fNDimuons>0){
                        Tracklets_Dist_V0MPeripheral_Dimu->Fill(NumberCloseEtaTracklets);
                    }
                }
                
                // === SPD CENTRAL not V0M ===
                
                if(isV0MCentral == kFALSE && isSPDCentral == kTRUE){
                    V0M_vs_SPD_OnlySPDCentral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                    for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                        trac = (TrackletLight*)fTracklets->At(j);
                        if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                Tracklets_OnlySPDCentral->Fill(trac->fPhi, trac->fEta);
                        }
                        
                    }
                    Tracklets_Dist_OnlySPDCentral->Fill(NumberCloseEtaTracklets);
                    Z_OnlySPDCentral->Fill(fEvent->fVertexZ);
                }
                
                if(isSPDCentral == kTRUE){
                                   V0M_vs_SPD_SPDCentral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                                   for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                                       trac = (TrackletLight*)fTracklets->At(j);
                                       if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                               Tracklets_SPDCentral->Fill(trac->fPhi, trac->fEta);
                                       }
                                       
                                   }
                                   Tracklets_Dist_SPDCentral->Fill(NumberCloseEtaTracklets);
                                   Z_SPDCentral->Fill(fEvent->fVertexZ);
                                    if(fEvent->fNDimuons>0){
                                        Tracklets_Dist_SPDCentral_Dimu->Fill(NumberCloseEtaTracklets);
                                    }
                               }
                
                // === SPD PERIPH not V0M ===
                
                if(isV0MPeripheral == kFALSE && isSPDPeripheral == kTRUE){
                    V0M_vs_SPD_OnlySPDPeripheral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                    for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                        trac = (TrackletLight*)fTracklets->At(j);
                        if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                Tracklets_OnlySPDPeripheral->Fill(trac->fPhi, trac->fEta);
                        }
                        
                    }
                    Tracklets_Dist_OnlySPDPeripheral->Fill(NumberCloseEtaTracklets);
                    Z_OnlySPDPeripheral->Fill(fEvent->fVertexZ);
                }
                
                if(isSPDPeripheral == kTRUE){
                                   V0M_vs_SPD_SPDPeripheral->Fill(fEvent->fV0MValue, fEvent->fSPDTrackletsValue);
                                   for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                                       trac = (TrackletLight*)fTracklets->At(j);
                                       if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                                               Tracklets_SPDPeripheral->Fill(trac->fPhi, trac->fEta);
                                       }
                                       
                                   }
                                   Tracklets_Dist_SPDPeripheral->Fill(NumberCloseEtaTracklets);
                                   Z_SPDPeripheral->Fill(fEvent->fVertexZ);
                                    if(fEvent->fNDimuons>0){
                                        Tracklets_Dist_SPDPeripheral_Dimu->Fill(NumberCloseEtaTracklets);
                                    }
                               }
                
                
                
                

               }
           
        }
        
        char hname[100];
       
        sprintf(hname, "V0M_vs_SPD_All - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_All = new TCanvas(hname);
        V0M_vs_SPD_All->SetTitle(hname);
        V0M_vs_SPD_All->Draw("colz");
        
        sprintf(hname, "Tracklets_All - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_All = new TCanvas(hname);
        Tracklets_All->SetTitle(hname);
        Tracklets_All->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_All - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_All = new TCanvas(hname);
        Tracklets_Dist_All->SetTitle(hname);
        Tracklets_Dist_All->Draw("colz");
        
        sprintf(hname, "Z_All - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_All = new TCanvas(hname);
        Z_All->SetTitle(hname);
        Z_All->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_BothCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_BothCentral = new TCanvas(hname);
        V0M_vs_SPD_BothCentral->SetTitle(hname);
        V0M_vs_SPD_BothCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_BothCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_BothCentral = new TCanvas(hname);
        Tracklets_BothCentral->SetTitle(hname);
        Tracklets_BothCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_BothCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_BothCentral = new TCanvas(hname);
        Tracklets_Dist_BothCentral->SetTitle(hname);
        Tracklets_Dist_BothCentral->Draw("colz");
        
        sprintf(hname, "Z_BothCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_BothCentral = new TCanvas(hname);
        Z_BothCentral->SetTitle(hname);
        Z_BothCentral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_BothPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_BothPeripheral = new TCanvas(hname);
        V0M_vs_SPD_BothPeripheral->SetTitle(hname);
        V0M_vs_SPD_BothPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_BothPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_BothPeripheral = new TCanvas(hname);
        Tracklets_BothPeripheral->SetTitle(hname);
        Tracklets_BothPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_BothPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_BothPeripheral = new TCanvas(hname);
        Tracklets_Dist_BothPeripheral->SetTitle(hname);
        Tracklets_Dist_BothPeripheral->Draw("colz");
        
        sprintf(hname, "Z_BothPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_BothPeripheral = new TCanvas(hname);
        Z_BothPeripheral->SetTitle(hname);
        Z_BothPeripheral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_OnlyV0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_OnlyV0MCentral = new TCanvas(hname);
        V0M_vs_SPD_OnlyV0MCentral->SetTitle(hname);
        V0M_vs_SPD_OnlyV0MCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_OnlyV0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_OnlyV0MCentral = new TCanvas(hname);
        Tracklets_OnlyV0MCentral->SetTitle(hname);
        Tracklets_OnlyV0MCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_OnlyV0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_OnlyV0MCentral = new TCanvas(hname);
        Tracklets_Dist_OnlyV0MCentral->SetTitle(hname);
        Tracklets_Dist_OnlyV0MCentral->Draw("colz");
        
        sprintf(hname, "Z_OnlyV0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_OnlyV0MCentral = new TCanvas(hname);
        Z_OnlyV0MCentral->SetTitle(hname);
        Z_OnlyV0MCentral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_OnlyV0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_OnlyV0MPeripheral = new TCanvas(hname);
        V0M_vs_SPD_OnlyV0MPeripheral->SetTitle(hname);
        V0M_vs_SPD_OnlyV0MPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_OnlyV0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_OnlyV0MPeripheral = new TCanvas(hname);
        Tracklets_OnlyV0MPeripheral->SetTitle(hname);
        Tracklets_OnlyV0MPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_OnlyV0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_OnlyV0MPeripheral = new TCanvas(hname);
        Tracklets_Dist_OnlyV0MPeripheral->SetTitle(hname);
        Tracklets_Dist_OnlyV0MPeripheral->Draw("colz");
        
        sprintf(hname, "Z_OnlyV0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_OnlyV0MPeripheral = new TCanvas(hname);
        Z_OnlyV0MPeripheral->SetTitle(hname);
        Z_OnlyV0MPeripheral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_OnlySPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_OnlySPDCentral = new TCanvas(hname);
        V0M_vs_SPD_OnlySPDCentral->SetTitle(hname);
        V0M_vs_SPD_OnlySPDCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_OnlySPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_OnlySPDCentral = new TCanvas(hname);
        Tracklets_OnlySPDCentral->SetTitle(hname);
        Tracklets_OnlySPDCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_OnlySPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_OnlySPDCentral = new TCanvas(hname);
        Tracklets_Dist_OnlySPDCentral->SetTitle(hname);
        Tracklets_Dist_OnlySPDCentral->Draw("colz");
        
        sprintf(hname, "Z_OnlySPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_OnlySPDCentral = new TCanvas(hname);
        Z_OnlySPDCentral->SetTitle(hname);
        Z_OnlySPDCentral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_OnlySPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_OnlySPDPeripheral = new TCanvas(hname);
        V0M_vs_SPD_OnlySPDPeripheral->SetTitle(hname);
        V0M_vs_SPD_OnlySPDPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_OnlySPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_OnlySPDPeripheral = new TCanvas(hname);
        Tracklets_OnlySPDPeripheral->SetTitle(hname);
        Tracklets_OnlySPDPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_OnlySPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_OnlySPDPeripheral = new TCanvas(hname);
        Tracklets_Dist_OnlySPDPeripheral->SetTitle(hname);
        Tracklets_Dist_OnlySPDPeripheral->Draw("colz");
        
        sprintf(hname, "Z_OnlySPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_OnlySPDPeripheral = new TCanvas(hname);
        Z_OnlySPDPeripheral->SetTitle(hname);
        Z_OnlySPDPeripheral->Draw("colz");
        
        
        sprintf(hname, "V0M_vs_SPD_V0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_V0MCentral = new TCanvas(hname);
        V0M_vs_SPD_V0MCentral->SetTitle(hname);
        V0M_vs_SPD_V0MCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_V0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_V0MCentral = new TCanvas(hname);
        Tracklets_V0MCentral->SetTitle(hname);
        Tracklets_V0MCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_V0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_V0MCentral = new TCanvas(hname);
        Tracklets_Dist_V0MCentral->SetTitle(hname);
        Tracklets_Dist_V0MCentral->Scale(1./Tracklets_Dist_V0MCentral->GetEntries());
        Tracklets_Dist_V0MCentral->Draw("colz");
        Tracklets_Dist_V0MCentral_Dimu->SetLineColor(kRed);
        Tracklets_Dist_V0MCentral_Dimu->Scale(1./Tracklets_Dist_V0MCentral_Dimu->GetEntries());
        Tracklets_Dist_V0MCentral_Dimu->Draw("SAME");
        
        sprintf(hname, "Z_V0MCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_V0MCentral = new TCanvas(hname);
        Z_V0MCentral->SetTitle(hname);
        Z_V0MCentral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_V0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_V0MPeripheral = new TCanvas(hname);
        V0M_vs_SPD_V0MPeripheral->SetTitle(hname);
        V0M_vs_SPD_V0MPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_V0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_V0MPeripheral = new TCanvas(hname);
        Tracklets_V0MPeripheral->SetTitle(hname);
        Tracklets_V0MPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_V0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_V0MPeripheral = new TCanvas(hname);
        Tracklets_Dist_V0MPeripheral->SetTitle(hname);
        Tracklets_Dist_V0MPeripheral->Scale(1./Tracklets_Dist_V0MPeripheral->GetEntries());
        Tracklets_Dist_V0MPeripheral->Draw("colz");
        Tracklets_Dist_V0MPeripheral_Dimu->SetLineColor(kRed);
        Tracklets_Dist_V0MPeripheral_Dimu->Scale(1./Tracklets_Dist_V0MPeripheral_Dimu->GetEntries());
        Tracklets_Dist_V0MPeripheral_Dimu->Draw("SAME");
        
        sprintf(hname, "Z_V0MPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_V0MPeripheral = new TCanvas(hname);
        Z_V0MPeripheral->SetTitle(hname);
        Z_V0MPeripheral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_SPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_SPDCentral = new TCanvas(hname);
        V0M_vs_SPD_SPDCentral->SetTitle(hname);
        V0M_vs_SPD_SPDCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_SPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_SPDCentral = new TCanvas(hname);
        Tracklets_SPDCentral->SetTitle(hname);
        Tracklets_SPDCentral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_SPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_SPDCentral = new TCanvas(hname);
        Tracklets_Dist_SPDCentral->SetTitle(hname);
        Tracklets_Dist_SPDCentral->Scale(1./Tracklets_Dist_SPDCentral->GetEntries());
        Tracklets_Dist_SPDCentral->Draw("colz");
        Tracklets_Dist_SPDCentral_Dimu->SetLineColor(kRed);
        Tracklets_Dist_SPDCentral_Dimu->Scale(1./Tracklets_Dist_SPDCentral_Dimu->GetEntries());
        Tracklets_Dist_SPDCentral_Dimu->Draw("SAME");
        
        sprintf(hname, "Z_SPDCentral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_SPDCentral = new TCanvas(hname);
        Z_SPDCentral->SetTitle(hname);
        Z_SPDCentral->Draw("colz");
        
        sprintf(hname, "V0M_vs_SPD_SPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M_vs_SPD_SPDPeripheral = new TCanvas(hname);
        V0M_vs_SPD_SPDPeripheral->SetTitle(hname);
        V0M_vs_SPD_SPDPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_SPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_SPDPeripheral = new TCanvas(hname);
        Tracklets_SPDPeripheral->SetTitle(hname);
        Tracklets_SPDPeripheral->Draw("colz");
        
        sprintf(hname, "Tracklets_Dist_SPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cTracklets_Dist_SPDPeripheral = new TCanvas(hname);
        Tracklets_Dist_SPDPeripheral->SetTitle(hname);
        Tracklets_Dist_SPDPeripheral->Scale(1./Tracklets_Dist_SPDPeripheral->GetEntries());
        Tracklets_Dist_SPDPeripheral->Draw("colz");
        Tracklets_Dist_SPDPeripheral_Dimu->SetLineColor(kRed);
        Tracklets_Dist_SPDPeripheral_Dimu->Scale(1./Tracklets_Dist_SPDPeripheral_Dimu->GetEntries());
        Tracklets_Dist_SPDPeripheral_Dimu->Draw("SAME");
        
        sprintf(hname, "Z_SPDPeripheral - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cZ_SPDPeripheral = new TCanvas(hname);
        Z_SPDPeripheral->SetTitle(hname);
        Z_SPDPeripheral->Draw("colz");
        
        
    }

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}
    
Double_t myFunc(double x, double y){
    return hPlotTklEtaPhi->GetBinContent(floor(720*x/(2*TMath::Pi())),floor(100*(y+1)/2)); }

Double_t myFunc2(double x, double y){
    return hPlotDimuEtaPhi->GetBinContent(floor(720*x/(2*TMath::Pi())),floor(100*(y+6)/6)); }

Double_t myEtaLow(double x){
    return TMath::Log(TMath::Tan(0.5*TMath::ATan(7.6/(14.1+(floor(x)+0.5))))) + 0.1; }

Double_t myEtaHigh(double x){
return -(TMath::Log(TMath::Tan(0.5*TMath::ATan(7.6/(14.1-(floor(x)+0.5)))))) - 0.05; }


void PlotMass(){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    
// ************************************
// Définitions de paramètres          *
// ************************************

    double ZvtxCut = 10;
    double SigmaZvtxCut = 0.25;
    double DPhiCut = 0.01;
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 12;
    double MinMultCentral = 37;
    double MaxMultPeriph = 23;
    double CentSPDTrackletsCentral = 1;
    double CentSPDTrackletsPeriph = 7;

    double LowJPsiMass = 3.0;
    double HighJPsiMass = 3.3;

    double MinInvMass = 1.0;
    double MaxInvMass = 5;
    const int NbinsDimuInvMass = 500;
    double SizeBinInvMass = (MaxInvMass-MinInvMass)/NbinsDimuInvMass;

    double MinDeltaPhi = -TMath::Pi()/2;
    double MaxDeltaPhi = 1.5*TMath::Pi();
    const int NbinsDeltaPhi = 12;
    double SizeBinDeltaPhi = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhi;
    
    double MinDeltaEta = 0;
    double MaxDeltaEta = 6;
    const int NbinsDeltaEta = 20;
    double SizeBinDeltaEta = (MaxDeltaEta-MinDeltaEta)/NbinsDeltaEta;
    int barWidth = 50;
    
    Char_t FitFileName[500];
    //   Char_t AssociateFileName[200];

       sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_Run2_Manuel.root");
    
    //Char_t Group_Period[50] = "Group1_LHC16h";
   // Char_t *arrayOfPeriods[] = {"Group5_LHC17l_CINT","Group5_LHC17m_CINT","Group5_LHC17o_CINT","Group5_LHC17r_CINT","Group5_LHC18c_CINT","Group5_LHC18d_CINT","Group5_LHC18e_CINT","Group5_LHC18f_CINT"};
    Char_t *arrayOfPeriods[] = {"Cvetan_LHC16r_PU2"};
   // Char_t *arrayOfPeriods[] = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l","Group2_LHC17h","Group3_LHC17h","Group4_LHC17k","Group4_LHC18l","Group4_LHC18m","Group4_LHC18o","Group4_LHC18p","Group5_LHC17l","Group5_LHC17m","Group5_LHC17o","Group5_LHC17r","Group5_LHC18c","Group5_LHC18d","Group5_LHC18e","Group5_LHC18f","Group6_LHC18m","Group7_LHC18m","Group8_LHC18m","Group9_LHC18m","Group10_LHC18m","Group11_LHC18m","Group12_LHC18m"};
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
    
    TH1I* hnseg(NULL);
    TH2I* hPU(NULL);
    TF1 *fPUCut = new TF1("fPUCut","100+4*x",0,200);
    
    
    hnseg = new TH1I("hnseg",
                     "Invariant mass of dimuon",
                     NbinsDimuInvMass,MinInvMass,MaxInvMass);
    hnseg->SetXTitle("Mass of dimuon (GeV/c^{2})");
    hPU = new TH2I("hPU",
                     "Clusters vs tracklets",
                     200,0,200,1000,0,1000);
    hPU->SetXTitle("Tracklets");
    hPU->SetYTitle("Clusters");
    

// *************************
// Analyse                 *
// *************************
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        
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
            
            if(fEvent->fIsPileupFromSPDMultBins || fEvent->fVertexNC < 1 || fEvent->fNDimuons <1 || fEvent->fPassPhysicsSelection == 0 || fEvent->fIsPileupClustVsTkl == 1){//}|| fEvent->fIsPileupClustVsTkl == 0 || fEvent->fPassPhysicsSelection == 0){ //|| fEvent->fNPileupVtx != 0
                continue;
            }
            
            hPU->Fill(fEvent->fNTracklets, fEvent->fSPDClustersValue);
     //       cout << "Looking at Event " << i << endl;
//            cout << "fPassPhysicsSelection : " << fEvent->fPassPhysicsSelection << endl;
//            cout << "fPassTriggerSelection : " << fEvent->fPassTriggerSelection << endl;
//            cout << "fTriggerCMUL7 : " << fEvent->fTriggerCMUL7 << endl;
//            cout << "fTriggerCINT7 : " << fEvent->fTriggerCINT7 << endl;
            //cout << "tracks->GetEntries() " << tracks->GetEntries() << endl;
          //  cout << "Computing the number of tracklets in Event " << i << endl;
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
//            for (Int_t j=0; j<fEvent->fNDimuons; j++) {
//                dimu = (DimuonLight*)fDimuons->At(j);
//                if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut) && (dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)){
//
//                    hnseg->Fill(dimu->fInvMass);
//
//
//                }
//            }

            
            
        }
        
        
    }
    
    cout << "Reading of events finished" <<endl;
    
//    TFile *f = new TFile(FitFileName,"UPDATE");
//
//    hPU->Write();
//    f->Close();
//
    hPU->Draw("colz");
    fPUCut->Draw("same");

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    
    
    
    
}


void PlotQuickie(){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    
// ************************************
// Définitions de paramètres          *
// ************************************

    double ZvtxCut = 10;
    double DPhiCut = 0.01;
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 12;
    double MinMultCentral = 37;
    double MaxMultPeriph = 23;
    double CentSPDTrackletsCentral = 1;
    double CentSPDTrackletsPeriph = 7;
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
   // Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
  // Char_t *arrayOfPeriods[] = {"Group1_LHC16h_CINT_AllEst","Group1_LHC16j_CINT_AllEst","Group1_LHC16p_CINT_AllEst","Group1_LHC17i_CINT_AllEst","Group2_LHC17h_CINT_AllEst","Group3_LHC17h_CINT_AllEst","Group4_LHC18o_CINT_AllEst","Group5_LHC17l_CINT_AllEst"};
    Char_t *arrayOfPeriods[] = {"Group1_LHC16h_CINT_AllEst"};
   // Char_t *arrayOfPeriods[] = {"Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose"};//, "Group1_LHC17k", "Group4_LHC17k", "Group4_LHC18o"};
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
    
    TH1F* SPDTracklets_Dist(NULL);
    TH1F* Z_Dist(NULL);
    TProfile* SPDTracklets_vs_Z(NULL);

    
    
   
   
// *************************
// Analyse                 *
// *************************
    
    TCanvas* multiplecan=new TCanvas("V0M Distribution Groups");
    TCanvas* cSPDTracklets_Dist_Groups=new TCanvas("SPDTracklets Distribution Groups");
    TCanvas* cV0M_Dist_Groups_norm=new TCanvas("V0M norm Distribution Groups");
    TCanvas* cSPDTracklets_Dist_Groups_norm=new TCanvas("SPDTracklets norm Distribution Groups");
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Sauvegarde\ /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        //sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        
        
        SPDTracklets_Dist = new TH1F("SPDTracklets_Dist",
                         "SPDTracklets_Dist",
                         1000,-0.5, 999.5);
        SPDTracklets_Dist->SetXTitle("SPDTracklets");
        SPDTracklets_Dist->SetYTitle("Count");
        
        Z_Dist = new TH1F("Z_Dist",
                         "Z_Dist",
                         20,-10, 10);
        Z_Dist->SetXTitle("Z");
        Z_Dist->SetYTitle("Count");
        
        SPDTracklets_vs_Z = new TProfile("SPDTracklets_vs_Z",
                            "SPDTracklets_vs_Z",
                            81,-10, 10);
           SPDTracklets_vs_Z->SetXTitle("Z");
           SPDTracklets_vs_Z->SetYTitle("SPDTracklets");
        
        
        
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
            
            if(i%100 != 0){
                continue;
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
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
                SPDTracklets_Dist->Fill(fEvent->fSPDTrackletsValue);
                Z_Dist->Fill(fEvent->fVertexZ);
                SPDTracklets_vs_Z->Fill(fEvent->fVertexZ, fEvent->fSPDTrackletsValue);
               }
           
        }
      
        
        //SPDTracklets_vs_Z->Divide(Z_Dist);
        
        multiplecan->cd();

        SPDTracklets_vs_Z->SetLineColor(tree_idx+1);
        SPDTracklets_vs_Z->Draw("same");
        
        char message[80];
        sprintf(message, "%s", arrayOfPeriods[tree_idx]);
        legend->AddEntry(SPDTracklets_vs_Z,message);
        
        
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    
    
    multiplecan->cd();
    legend->Draw();

}

void PlotQuickie2(){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    
// ************************************
// Définitions de paramètres          *
// ************************************

    double ZvtxCut = 10;
    double DPhiCut = 0.01;
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 12;
    double MinMultCentral = 37;
    double MaxMultPeriph = 23;
    double CentSPDTrackletsCentral = 1;
    double CentSPDTrackletsPeriph = 7;
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
   // Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
   Char_t *arrayOfPeriods[] = {"Group1_LHC16h_CINT_AllEst","Group1_LHC16j_CINT_AllEst"};
   // Char_t *arrayOfPeriods[] = {"Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose"};//, "Group1_LHC17k", "Group4_LHC17k", "Group4_LHC18o"};
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
    
    TH2F* SPDTracklets_vs_Cent(NULL);

    
    
   
   
// *************************
// Analyse                 *
// *************************
    
    TCanvas* multiplecan=new TCanvas("SPDTracklets 2D");
    TCanvas* cSPDTracklets_Dist_Groups=new TCanvas("SPDTracklets Distribution Groups");
    TCanvas* cV0M_Dist_Groups_norm=new TCanvas("V0M norm Distribution Groups");
    TCanvas* cSPDTracklets_Dist_Groups_norm=new TCanvas("SPDTracklets norm Distribution Groups");
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    SPDTracklets_vs_Cent = new TH2F("SPDTracklets_vs_Cent",
                     "SPDTracklets_vs_Cent",
                     200,-0.5, 199.5,100,0,100);
    SPDTracklets_vs_Cent->SetXTitle("SPDTracklets");
    SPDTracklets_vs_Cent->SetYTitle("Centrality");
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Sauvegarde\ /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        //sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        
    

        
        
        
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
            
            if(i%100 != 0){
                continue;
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
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
                SPDTracklets_vs_Cent->Fill(fEvent->fSPDTrackletsValue, fEvent->fCentralitySPDTracklets);
               }
           
        }
      
        
        //SPDTracklets_vs_Z->Divide(Z_Dist);
        
        
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    
    
    multiplecan->cd();

    SPDTracklets_vs_Cent->Draw("colz");
    multiplecan->SaveAs("~/Mult.root");

}


void PlotSaucisson(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
    //Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
    Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
  //  Char_t *arrayOfPeriods[] = {"Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose"};//, "Group1_LHC17k", "Group4_LHC17k", "Group4_LHC18o"};
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
    
    TH1F* V0M_0_100(NULL);
    TH1F* V0M_0_1(NULL);
    TH1F* V0M_1_3(NULL);
    TH1F* V0M_3_5(NULL);
    TH1F* V0M_5_10(NULL);
    TH1F* V0M_10_15(NULL);
    TH1F* V0M_15_20(NULL);
    TH1F* V0M_20_30(NULL);
    TH1F* V0M_30_40(NULL);
    TH1F* V0M_40_50(NULL);
    TH1F* V0M_50_60(NULL);
    TH1F* V0M_60_70(NULL);
    TH1F* V0M_70_80(NULL);
    TH1F* V0M_80_90(NULL);
    TH1F* V0M_90_100(NULL);
    TH1F* SPDT_0_100(NULL);
    TH1F* SPDT_0_1(NULL);
    TH1F* SPDT_1_3(NULL);
    TH1F* SPDT_3_5(NULL);
    TH1F* SPDT_5_10(NULL);
    TH1F* SPDT_10_15(NULL);
    TH1F* SPDT_15_20(NULL);
    TH1F* SPDT_20_30(NULL);
    TH1F* SPDT_30_40(NULL);
    TH1F* SPDT_40_50(NULL);
    TH1F* SPDT_50_60(NULL);
    TH1F* SPDT_60_70(NULL);
    TH1F* SPDT_70_80(NULL);
    TH1F* SPDT_80_90(NULL);
    TH1F* SPDT_90_100(NULL);
    
    
    
    
    
    V0M_0_100 = new TH1F("V0M_0_100",
                     "V0M_0_100",
                     200,-0.5, 199.5);
    V0M_0_100->SetXTitle("Ntkl");
    V0M_0_100->SetYTitle("Count");
    V0M_0_1 = new TH1F("V0M_0_1",
                     "V0M_0_1",
                     200,-0.5, 199.5);
    V0M_0_1->SetXTitle("Ntkl");
    V0M_0_1->SetYTitle("Count");
    V0M_1_3 = new TH1F("V0M_1_3",
                     "V0M_1_3",
                     200,-0.5, 199.5);
    V0M_1_3->SetXTitle("Ntkl");
    V0M_1_3->SetYTitle("Count");
    V0M_3_5 = new TH1F("V0M_3_5",
                     "V0M_3_5",
                     200,-0.5, 199.5);
    V0M_3_5->SetXTitle("Ntkl");
    V0M_3_5->SetYTitle("Count");
    V0M_5_10 = new TH1F("V0M_5_10",
                     "V0M_5_10",
                     200,-0.5, 199.5);
    V0M_5_10->SetXTitle("Ntkl");
    V0M_5_10->SetYTitle("Count");
    V0M_10_15 = new TH1F("V0M_10_15",
                     "V0M_10_15",
                     200,-0.5, 199.5);
    V0M_10_15->SetXTitle("Ntkl");
    V0M_10_15->SetYTitle("Count");
    V0M_15_20 = new TH1F("V0M_15_20",
                     "V0M_15_20",
                     200,-0.5, 199.5);
    V0M_15_20->SetXTitle("Ntkl");
    V0M_15_20->SetYTitle("Count");
    V0M_20_30 = new TH1F("V0M_20_30",
                     "V0M_20_30",
                     200,-0.5, 199.5);
    V0M_20_30->SetXTitle("Ntkl");
    V0M_20_30->SetYTitle("Count");
    V0M_30_40 = new TH1F("V0M_30_40",
                     "V0M_30_40",
                     200,-0.5, 199.5);
    V0M_30_40->SetXTitle("Ntkl");
    V0M_30_40->SetYTitle("Count");
    V0M_40_50 = new TH1F("V0M_40_50",
                     "V0M_40_50",
                     200,-0.5, 199.5);
    V0M_40_50->SetXTitle("Ntkl");
    V0M_40_50->SetYTitle("Count");
    V0M_50_60 = new TH1F("V0M_50_60",
                     "V0M_50_60",
                     200,-0.5, 199.5);
    V0M_50_60->SetXTitle("Ntkl");
    V0M_50_60->SetYTitle("Count");
    V0M_60_70 = new TH1F("V0M_60_70",
                     "V0M_60_70",
                     200,-0.5, 199.5);
    V0M_60_70->SetXTitle("Ntkl");
    V0M_60_70->SetYTitle("Count");
    V0M_70_80 = new TH1F("V0M_70_80",
                     "V0M_70_80",
                     200,-0.5, 199.5);
    V0M_70_80->SetXTitle("Ntkl");
    V0M_70_80->SetYTitle("Count");
    V0M_80_90 = new TH1F("V0M_80_90",
                     "V0M_80_90",
                     200,-0.5, 199.5);
    V0M_80_90->SetXTitle("Ntkl");
    V0M_80_90->SetYTitle("Count");
    V0M_90_100 = new TH1F("V0M_90_100",
                     "V0M_90_100",
                     200,-0.5, 199.5);
    V0M_90_100->SetXTitle("Ntkl");
    V0M_90_100->SetYTitle("Count");
    
    SPDT_0_100 = new TH1F("SPDT_0_100",
                     "SPDT_0_100",
                     200,-0.5, 199.5);
    SPDT_0_100->SetXTitle("Ntkl");
    SPDT_0_100->SetYTitle("Count");
    SPDT_0_1 = new TH1F("SPDT_0_1",
                     "SPDT_0_1",
                     200,-0.5, 199.5);
    SPDT_0_1->SetXTitle("Ntkl");
    SPDT_0_1->SetYTitle("Count");
    SPDT_1_3 = new TH1F("SPDT_1_3",
                     "SPDT_1_3",
                     200,-0.5, 199.5);
    SPDT_1_3->SetXTitle("Ntkl");
    SPDT_1_3->SetYTitle("Count");
    SPDT_3_5 = new TH1F("SPDT_3_5",
                     "SPDT_3_5",
                     200,-0.5, 199.5);
    SPDT_3_5->SetXTitle("Ntkl");
    SPDT_3_5->SetYTitle("Count");
    SPDT_5_10 = new TH1F("SPDT_5_10",
                     "SPDT_5_10",
                     200,-0.5, 199.5);
    SPDT_5_10->SetXTitle("Ntkl");
    SPDT_5_10->SetYTitle("Count");
    SPDT_10_15 = new TH1F("SPDT_10_15",
                     "SPDT_10_15",
                     200,-0.5, 199.5);
    SPDT_10_15->SetXTitle("Ntkl");
    SPDT_10_15->SetYTitle("Count");
    SPDT_15_20 = new TH1F("SPDT_15_20",
                     "SPDT_15_20",
                     200,-0.5, 199.5);
    SPDT_15_20->SetXTitle("Ntkl");
    SPDT_15_20->SetYTitle("Count");
    SPDT_20_30 = new TH1F("SPDT_20_30",
                     "SPDT_20_30",
                     200,-0.5, 199.5);
    SPDT_20_30->SetXTitle("Ntkl");
    SPDT_20_30->SetYTitle("Count");
    SPDT_30_40 = new TH1F("SPDT_30_40",
                     "SPDT_30_40",
                     200,-0.5, 199.5);
    SPDT_30_40->SetXTitle("Ntkl");
    SPDT_30_40->SetYTitle("Count");
    SPDT_40_50 = new TH1F("SPDT_40_50",
                     "SPDT_40_50",
                     200,-0.5, 199.5);
    SPDT_40_50->SetXTitle("Ntkl");
    SPDT_40_50->SetYTitle("Count");
    SPDT_50_60 = new TH1F("SPDT_50_60",
                     "SPDT_50_60",
                     200,-0.5, 199.5);
    SPDT_50_60->SetXTitle("Ntkl");
    SPDT_50_60->SetYTitle("Count");
    SPDT_60_70 = new TH1F("SPDT_60_70",
                     "SPDT_60_70",
                     200,-0.5, 199.5);
    SPDT_60_70->SetXTitle("Ntkl");
    SPDT_60_70->SetYTitle("Count");
    SPDT_70_80 = new TH1F("SPDT_70_80",
                     "SPDT_70_80",
                     200,-0.5, 199.5);
    SPDT_70_80->SetXTitle("Ntkl");
    SPDT_70_80->SetYTitle("Count");
    SPDT_80_90 = new TH1F("SPDT_80_90",
                     "SPDT_80_90",
                     200,-0.5, 199.5);
    SPDT_80_90->SetXTitle("Ntkl");
    SPDT_80_90->SetYTitle("Count");
    SPDT_90_100 = new TH1F("SPDT_90_100",
                     "SPDT_90_100",
                     200,-0.5, 199.5);
    SPDT_90_100->SetXTitle("Ntkl");
    SPDT_90_100->SetYTitle("Count");
    
    
   
// *************************
// Analyse                 *
// *************************
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
       // sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
        
// /Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT

sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);        
        
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
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
            
                V0M_0_100->Fill(NumberCloseEtaTracklets);
                SPDT_0_100->Fill(NumberCloseEtaTracklets);
                
                if(fEvent->fCentralityV0M<1){
                    V0M_0_1->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=1 && fEvent->fCentralityV0M<3){
                    V0M_1_3->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=3 && fEvent->fCentralityV0M<5){
                    V0M_3_5->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=5 && fEvent->fCentralityV0M<10){
                    V0M_5_10->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=10 && fEvent->fCentralityV0M<15){
                    V0M_10_15->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=15 && fEvent->fCentralityV0M<20){
                    V0M_15_20->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=20 && fEvent->fCentralityV0M<30){
                    V0M_20_30->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=30 && fEvent->fCentralityV0M<40){
                    V0M_30_40->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=40 && fEvent->fCentralityV0M<50){
                    V0M_40_50->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=50 && fEvent->fCentralityV0M<60){
                    V0M_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=50 && fEvent->fCentralityV0M<60){
                    V0M_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=60 && fEvent->fCentralityV0M<70){
                    V0M_60_70->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=70 && fEvent->fCentralityV0M<80){
                    V0M_70_80->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=80 && fEvent->fCentralityV0M<90){
                    V0M_80_90->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=90){
                    V0M_90_100->Fill(NumberCloseEtaTracklets);
                }
                
                if(fEvent->fCentralitySPDTracklets<1){
                    SPDT_0_1->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=1 && fEvent->fCentralitySPDTracklets<3){
                    SPDT_1_3->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=3 && fEvent->fCentralitySPDTracklets<5){
                    SPDT_3_5->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=5 && fEvent->fCentralitySPDTracklets<10){
                    SPDT_5_10->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=10 && fEvent->fCentralitySPDTracklets<15){
                    SPDT_10_15->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=15 && fEvent->fCentralitySPDTracklets<20){
                    SPDT_15_20->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=20 && fEvent->fCentralitySPDTracklets<30){
                    SPDT_20_30->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=30 && fEvent->fCentralitySPDTracklets<40){
                    SPDT_30_40->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=40 && fEvent->fCentralitySPDTracklets<50){
                    SPDT_40_50->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=50 && fEvent->fCentralitySPDTracklets<60){
                    SPDT_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=50 && fEvent->fCentralitySPDTracklets<60){
                    SPDT_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=60 && fEvent->fCentralitySPDTracklets<70){
                    SPDT_60_70->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=70 && fEvent->fCentralitySPDTracklets<80){
                    SPDT_70_80->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=80 && fEvent->fCentralitySPDTracklets<90){
                    SPDT_80_90->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=90){
                    SPDT_90_100->Fill(NumberCloseEtaTracklets);
                }
                

               }
           
        }
        
        

        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
            V0M_0_1->SetLineColor(gStyle->GetColorPalette(int(nColors*0)));
            V0M_0_1->SetFillColor(gStyle->GetColorPalette(int(nColors*0)));
            SPDT_0_1->SetLineColor(gStyle->GetColorPalette(int(nColors*0)));
            SPDT_0_1->SetFillColor(gStyle->GetColorPalette(int(nColors*0)));
        V0M_1_3->SetLineColor(gStyle->GetColorPalette(int(nColors*0.05)));
        V0M_1_3->SetFillColor(gStyle->GetColorPalette(int(nColors*0.05)));
        SPDT_1_3->SetLineColor(gStyle->GetColorPalette(int(nColors*0.05)));
        SPDT_1_3->SetFillColor(gStyle->GetColorPalette(int(nColors*0.05)));
        V0M_3_5->SetLineColor(gStyle->GetColorPalette(int(nColors*0.1)));
        V0M_3_5->SetFillColor(gStyle->GetColorPalette(int(nColors*0.1)));
        SPDT_3_5->SetLineColor(gStyle->GetColorPalette(int(nColors*0.1)));
        SPDT_3_5->SetFillColor(gStyle->GetColorPalette(int(nColors*0.1)));
        V0M_5_10->SetLineColor(gStyle->GetColorPalette(int(nColors*0.15)));
        V0M_5_10->SetFillColor(gStyle->GetColorPalette(int(nColors*0.15)));
        SPDT_5_10->SetLineColor(gStyle->GetColorPalette(int(nColors*0.15)));
        SPDT_5_10->SetFillColor(gStyle->GetColorPalette(int(nColors*0.15)));
        V0M_10_15->SetLineColor(gStyle->GetColorPalette(int(nColors*0.2)));
        V0M_10_15->SetFillColor(gStyle->GetColorPalette(int(nColors*0.2)));
        SPDT_10_15->SetLineColor(gStyle->GetColorPalette(int(nColors*0.2)));
        SPDT_10_15->SetFillColor(gStyle->GetColorPalette(int(nColors*0.2)));
        V0M_15_20->SetLineColor(gStyle->GetColorPalette(int(nColors*0.25)));
        V0M_15_20->SetFillColor(gStyle->GetColorPalette(int(nColors*0.25)));
        SPDT_15_20->SetLineColor(gStyle->GetColorPalette(int(nColors*0.25)));
        SPDT_15_20->SetFillColor(gStyle->GetColorPalette(int(nColors*0.25)));
        V0M_20_30->SetLineColor(gStyle->GetColorPalette(int(nColors*0.3)));
        V0M_20_30->SetFillColor(gStyle->GetColorPalette(int(nColors*0.3)));
        SPDT_20_30->SetLineColor(gStyle->GetColorPalette(int(nColors*0.3)));
        SPDT_20_30->SetFillColor(gStyle->GetColorPalette(int(nColors*0.3)));
        V0M_30_40->SetLineColor(gStyle->GetColorPalette(int(nColors*0.35)));
        V0M_30_40->SetFillColor(gStyle->GetColorPalette(int(nColors*0.35)));
        SPDT_30_40->SetLineColor(gStyle->GetColorPalette(int(nColors*0.35)));
        SPDT_30_40->SetFillColor(gStyle->GetColorPalette(int(nColors*0.35)));
        V0M_40_50->SetLineColor(gStyle->GetColorPalette(int(nColors*0.4)));
        V0M_40_50->SetFillColor(gStyle->GetColorPalette(int(nColors*0.4)));
        SPDT_40_50->SetLineColor(gStyle->GetColorPalette(int(nColors*0.4)));
        SPDT_40_50->SetFillColor(gStyle->GetColorPalette(int(nColors*0.4)));
        V0M_50_60->SetLineColor(gStyle->GetColorPalette(int(nColors*0.45)));
        V0M_50_60->SetFillColor(gStyle->GetColorPalette(int(nColors*0.45)));
        SPDT_50_60->SetLineColor(gStyle->GetColorPalette(int(nColors*0.45)));
        SPDT_50_60->SetFillColor(gStyle->GetColorPalette(int(nColors*0.45)));
        V0M_60_70->SetLineColor(gStyle->GetColorPalette(int(nColors*0.5)));
        V0M_60_70->SetFillColor(gStyle->GetColorPalette(int(nColors*0.5)));
        SPDT_60_70->SetLineColor(gStyle->GetColorPalette(int(nColors*0.5)));
        SPDT_60_70->SetFillColor(gStyle->GetColorPalette(int(nColors*0.5)));
        V0M_70_80->SetLineColor(gStyle->GetColorPalette(int(nColors*0.55)));
        V0M_70_80->SetFillColor(gStyle->GetColorPalette(int(nColors*0.55)));
        SPDT_70_80->SetLineColor(gStyle->GetColorPalette(int(nColors*0.55)));
        SPDT_70_80->SetFillColor(gStyle->GetColorPalette(int(nColors*0.55)));
        V0M_80_90->SetLineColor(gStyle->GetColorPalette(int(nColors*0.6)));
        V0M_80_90->SetFillColor(gStyle->GetColorPalette(int(nColors*0.6)));
        SPDT_80_90->SetLineColor(gStyle->GetColorPalette(int(nColors*0.6)));
        SPDT_80_90->SetFillColor(gStyle->GetColorPalette(int(nColors*0.6)));
        V0M_90_100->SetLineColor(gStyle->GetColorPalette(int(nColors*0.65)));
        V0M_90_100->SetFillColor(gStyle->GetColorPalette(int(nColors*0.65)));
        SPDT_90_100->SetLineColor(gStyle->GetColorPalette(int(nColors*0.65)));
        SPDT_90_100->SetFillColor(gStyle->GetColorPalette(int(nColors*0.65)));
        
               
        
        char hname[100];
        sprintf(hname, "V0M Distribution - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M=new TCanvas(hname);
        V0M_0_100->SetTitle(hname);
        V0M_0_100->Draw();
        V0M_0_1->Draw("same");
        V0M_1_3->Draw("same");
        V0M_3_5->Draw("same");
        V0M_5_10->Draw("same");
        V0M_10_15->Draw("same");
        V0M_15_20->Draw("same");
        V0M_20_30->Draw("same");
        V0M_30_40->Draw("same");
        V0M_40_50->Draw("same");
        V0M_50_60->Draw("same");
        V0M_60_70->Draw("same");
        V0M_70_80->Draw("same");
        V0M_80_90->Draw("same");
        V0M_90_100->Draw("same");
        
        sprintf(hname, "SPDT Distribution - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cSPDT=new TCanvas(hname);
        SPDT_0_100->SetTitle(hname);
        SPDT_0_100->Draw();
        SPDT_0_1->Draw("same");
        SPDT_1_3->Draw("same");
        SPDT_3_5->Draw("same");
        SPDT_5_10->Draw("same");
        SPDT_10_15->Draw("same");
        SPDT_15_20->Draw("same");
        SPDT_20_30->Draw("same");
        SPDT_30_40->Draw("same");
        SPDT_40_50->Draw("same");
        SPDT_50_60->Draw("same");
        SPDT_60_70->Draw("same");
        SPDT_70_80->Draw("same");
        SPDT_80_90->Draw("same");
        SPDT_90_100->Draw("same");
        
        
        
    }
    
    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}

void PlotMCdistrib(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
    //Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
    Char_t *arrayOfPeriods[] = {"Group1_LHC16hMCno"};
  //  Char_t *arrayOfPeriods[] = {"Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose"};//, "Group1_LHC17k", "Group4_LHC17k", "Group4_LHC18o"};
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
    
    TH1F* V0M_0_100(NULL);
    TH1F* V0M_0_1(NULL);
    TH1F* V0M_1_3(NULL);
    TH1F* V0M_3_5(NULL);
    TH1F* V0M_5_10(NULL);
    TH1F* V0M_10_15(NULL);
    TH1F* V0M_15_20(NULL);
    TH1F* V0M_20_30(NULL);
    TH1F* V0M_30_40(NULL);
    TH1F* V0M_40_50(NULL);
    TH1F* V0M_50_60(NULL);
    TH1F* V0M_60_70(NULL);
    TH1F* V0M_70_80(NULL);
    TH1F* V0M_80_90(NULL);
    TH1F* V0M_90_100(NULL);
    TH1F* SPDT_0_100(NULL);
    TH1F* SPDT_0_1(NULL);
    TH1F* SPDT_1_3(NULL);
    TH1F* SPDT_3_5(NULL);
    TH1F* SPDT_5_10(NULL);
    TH1F* SPDT_10_15(NULL);
    TH1F* SPDT_15_20(NULL);
    TH1F* SPDT_20_30(NULL);
    TH1F* SPDT_30_40(NULL);
    TH1F* SPDT_40_50(NULL);
    TH1F* SPDT_50_60(NULL);
    TH1F* SPDT_60_70(NULL);
    TH1F* SPDT_70_80(NULL);
    TH1F* SPDT_80_90(NULL);
    TH1F* SPDT_90_100(NULL);
    
    
    
    
    
    V0M_0_100 = new TH1F("V0M_0_100",
                     "V0M_0_100",
                     100,0, 100);
    V0M_0_100->SetXTitle("Percentile");
    V0M_0_100->SetYTitle("Count");
    V0M_0_1 = new TH1F("V0M_0_1",
                     "V0M_0_1",
                     100,0, 100);
    V0M_0_1->SetXTitle("Percentile");
    V0M_0_1->SetYTitle("Count");
    V0M_1_3 = new TH1F("V0M_1_3",
                     "V0M_1_3",
                     200,-0.5, 199.5);
    V0M_1_3->SetXTitle("Ntkl");
    V0M_1_3->SetYTitle("Count");
    V0M_3_5 = new TH1F("V0M_3_5",
                     "V0M_3_5",
                     200,-0.5, 199.5);
    V0M_3_5->SetXTitle("Ntkl");
    V0M_3_5->SetYTitle("Count");
    V0M_5_10 = new TH1F("V0M_5_10",
                     "V0M_5_10",
                     200,-0.5, 199.5);
    V0M_5_10->SetXTitle("Ntkl");
    V0M_5_10->SetYTitle("Count");
    V0M_10_15 = new TH1F("V0M_10_15",
                     "V0M_10_15",
                     200,-0.5, 199.5);
    V0M_10_15->SetXTitle("Ntkl");
    V0M_10_15->SetYTitle("Count");
    V0M_15_20 = new TH1F("V0M_15_20",
                     "V0M_15_20",
                     200,-0.5, 199.5);
    V0M_15_20->SetXTitle("Ntkl");
    V0M_15_20->SetYTitle("Count");
    V0M_20_30 = new TH1F("V0M_20_30",
                     "V0M_20_30",
                     200,-0.5, 199.5);
    V0M_20_30->SetXTitle("Ntkl");
    V0M_20_30->SetYTitle("Count");
    V0M_30_40 = new TH1F("V0M_30_40",
                     "V0M_30_40",
                     200,-0.5, 199.5);
    V0M_30_40->SetXTitle("Ntkl");
    V0M_30_40->SetYTitle("Count");
    V0M_40_50 = new TH1F("V0M_40_50",
                     "V0M_40_50",
                     200,-0.5, 199.5);
    V0M_40_50->SetXTitle("Ntkl");
    V0M_40_50->SetYTitle("Count");
    V0M_50_60 = new TH1F("V0M_50_60",
                     "V0M_50_60",
                     200,-0.5, 199.5);
    V0M_50_60->SetXTitle("Ntkl");
    V0M_50_60->SetYTitle("Count");
    V0M_60_70 = new TH1F("V0M_60_70",
                     "V0M_60_70",
                     200,-0.5, 199.5);
    V0M_60_70->SetXTitle("Ntkl");
    V0M_60_70->SetYTitle("Count");
    V0M_70_80 = new TH1F("V0M_70_80",
                     "V0M_70_80",
                     200,-0.5, 199.5);
    V0M_70_80->SetXTitle("Ntkl");
    V0M_70_80->SetYTitle("Count");
    V0M_80_90 = new TH1F("V0M_80_90",
                     "V0M_80_90",
                     200,-0.5, 199.5);
    V0M_80_90->SetXTitle("Ntkl");
    V0M_80_90->SetYTitle("Count");
    V0M_90_100 = new TH1F("V0M_90_100",
                     "V0M_90_100",
                     200,-0.5, 199.5);
    V0M_90_100->SetXTitle("Ntkl");
    V0M_90_100->SetYTitle("Count");
    
    SPDT_0_100 = new TH1F("SPDT_0_100",
                     "SPDT_0_100",
                     200,-0.5, 199.5);
    SPDT_0_100->SetXTitle("Ntkl");
    SPDT_0_100->SetYTitle("Count");
    SPDT_0_1 = new TH1F("SPDT_0_1",
                     "SPDT_0_1",
                     200,-0.5, 199.5);
    SPDT_0_1->SetXTitle("Ntkl");
    SPDT_0_1->SetYTitle("Count");
    SPDT_1_3 = new TH1F("SPDT_1_3",
                     "SPDT_1_3",
                     200,-0.5, 199.5);
    SPDT_1_3->SetXTitle("Ntkl");
    SPDT_1_3->SetYTitle("Count");
    SPDT_3_5 = new TH1F("SPDT_3_5",
                     "SPDT_3_5",
                     200,-0.5, 199.5);
    SPDT_3_5->SetXTitle("Ntkl");
    SPDT_3_5->SetYTitle("Count");
    SPDT_5_10 = new TH1F("SPDT_5_10",
                     "SPDT_5_10",
                     200,-0.5, 199.5);
    SPDT_5_10->SetXTitle("Ntkl");
    SPDT_5_10->SetYTitle("Count");
    SPDT_10_15 = new TH1F("SPDT_10_15",
                     "SPDT_10_15",
                     200,-0.5, 199.5);
    SPDT_10_15->SetXTitle("Ntkl");
    SPDT_10_15->SetYTitle("Count");
    SPDT_15_20 = new TH1F("SPDT_15_20",
                     "SPDT_15_20",
                     200,-0.5, 199.5);
    SPDT_15_20->SetXTitle("Ntkl");
    SPDT_15_20->SetYTitle("Count");
    SPDT_20_30 = new TH1F("SPDT_20_30",
                     "SPDT_20_30",
                     200,-0.5, 199.5);
    SPDT_20_30->SetXTitle("Ntkl");
    SPDT_20_30->SetYTitle("Count");
    SPDT_30_40 = new TH1F("SPDT_30_40",
                     "SPDT_30_40",
                     200,-0.5, 199.5);
    SPDT_30_40->SetXTitle("Ntkl");
    SPDT_30_40->SetYTitle("Count");
    SPDT_40_50 = new TH1F("SPDT_40_50",
                     "SPDT_40_50",
                     200,-0.5, 199.5);
    SPDT_40_50->SetXTitle("Ntkl");
    SPDT_40_50->SetYTitle("Count");
    SPDT_50_60 = new TH1F("SPDT_50_60",
                     "SPDT_50_60",
                     200,-0.5, 199.5);
    SPDT_50_60->SetXTitle("Ntkl");
    SPDT_50_60->SetYTitle("Count");
    SPDT_60_70 = new TH1F("SPDT_60_70",
                     "SPDT_60_70",
                     200,-0.5, 199.5);
    SPDT_60_70->SetXTitle("Ntkl");
    SPDT_60_70->SetYTitle("Count");
    SPDT_70_80 = new TH1F("SPDT_70_80",
                     "SPDT_70_80",
                     200,-0.5, 199.5);
    SPDT_70_80->SetXTitle("Ntkl");
    SPDT_70_80->SetYTitle("Count");
    SPDT_80_90 = new TH1F("SPDT_80_90",
                     "SPDT_80_90",
                     200,-0.5, 199.5);
    SPDT_80_90->SetXTitle("Ntkl");
    SPDT_80_90->SetYTitle("Count");
    SPDT_90_100 = new TH1F("SPDT_90_100",
                     "SPDT_90_100",
                     200,-0.5, 199.5);
    SPDT_90_100->SetXTitle("Ntkl");
    SPDT_90_100->SetYTitle("Count");
    
    
   
// *************************
// Analyse                 *
// *************************
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
        
// /Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT

//sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
        
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
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
          //  if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
            
                V0M_0_100->Fill(fEvent->fCentralityV0M);
                SPDT_0_100->Fill(NumberCloseEtaTracklets);
                
                if(fEvent->fNTracklets==0){
                    V0M_0_1->Fill(fEvent->fCentralityV0M);
                }
                if(fEvent->fCentralityV0M>=1 && fEvent->fCentralityV0M<3){
                    V0M_1_3->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=3 && fEvent->fCentralityV0M<5){
                    V0M_3_5->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=5 && fEvent->fCentralityV0M<10){
                    V0M_5_10->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=10 && fEvent->fCentralityV0M<15){
                    V0M_10_15->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=15 && fEvent->fCentralityV0M<20){
                    V0M_15_20->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=20 && fEvent->fCentralityV0M<30){
                    V0M_20_30->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=30 && fEvent->fCentralityV0M<40){
                    V0M_30_40->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=40 && fEvent->fCentralityV0M<50){
                    V0M_40_50->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=50 && fEvent->fCentralityV0M<60){
                    V0M_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=50 && fEvent->fCentralityV0M<60){
                    V0M_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=60 && fEvent->fCentralityV0M<70){
                    V0M_60_70->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=70 && fEvent->fCentralityV0M<80){
                    V0M_70_80->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=80 && fEvent->fCentralityV0M<90){
                    V0M_80_90->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralityV0M>=90){
                    V0M_90_100->Fill(NumberCloseEtaTracklets);
                }
                
                if(fEvent->fCentralitySPDTracklets<1){
                    SPDT_0_1->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=1 && fEvent->fCentralitySPDTracklets<3){
                    SPDT_1_3->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=3 && fEvent->fCentralitySPDTracklets<5){
                    SPDT_3_5->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=5 && fEvent->fCentralitySPDTracklets<10){
                    SPDT_5_10->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=10 && fEvent->fCentralitySPDTracklets<15){
                    SPDT_10_15->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=15 && fEvent->fCentralitySPDTracklets<20){
                    SPDT_15_20->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=20 && fEvent->fCentralitySPDTracklets<30){
                    SPDT_20_30->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=30 && fEvent->fCentralitySPDTracklets<40){
                    SPDT_30_40->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=40 && fEvent->fCentralitySPDTracklets<50){
                    SPDT_40_50->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=50 && fEvent->fCentralitySPDTracklets<60){
                    SPDT_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=50 && fEvent->fCentralitySPDTracklets<60){
                    SPDT_50_60->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=60 && fEvent->fCentralitySPDTracklets<70){
                    SPDT_60_70->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=70 && fEvent->fCentralitySPDTracklets<80){
                    SPDT_70_80->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=80 && fEvent->fCentralitySPDTracklets<90){
                    SPDT_80_90->Fill(NumberCloseEtaTracklets);
                }
                if(fEvent->fCentralitySPDTracklets>=90){
                    SPDT_90_100->Fill(NumberCloseEtaTracklets);
                }
                

              // }
           
        }
        
        

        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
            V0M_0_1->SetLineColor(gStyle->GetColorPalette(int(nColors*0)));
            V0M_0_1->SetFillColor(gStyle->GetColorPalette(int(nColors*0)));
            SPDT_0_1->SetLineColor(gStyle->GetColorPalette(int(nColors*0)));
            SPDT_0_1->SetFillColor(gStyle->GetColorPalette(int(nColors*0)));
        V0M_1_3->SetLineColor(gStyle->GetColorPalette(int(nColors*0.05)));
        V0M_1_3->SetFillColor(gStyle->GetColorPalette(int(nColors*0.05)));
        SPDT_1_3->SetLineColor(gStyle->GetColorPalette(int(nColors*0.05)));
        SPDT_1_3->SetFillColor(gStyle->GetColorPalette(int(nColors*0.05)));
        V0M_3_5->SetLineColor(gStyle->GetColorPalette(int(nColors*0.1)));
        V0M_3_5->SetFillColor(gStyle->GetColorPalette(int(nColors*0.1)));
        SPDT_3_5->SetLineColor(gStyle->GetColorPalette(int(nColors*0.1)));
        SPDT_3_5->SetFillColor(gStyle->GetColorPalette(int(nColors*0.1)));
        V0M_5_10->SetLineColor(gStyle->GetColorPalette(int(nColors*0.15)));
        V0M_5_10->SetFillColor(gStyle->GetColorPalette(int(nColors*0.15)));
        SPDT_5_10->SetLineColor(gStyle->GetColorPalette(int(nColors*0.15)));
        SPDT_5_10->SetFillColor(gStyle->GetColorPalette(int(nColors*0.15)));
        V0M_10_15->SetLineColor(gStyle->GetColorPalette(int(nColors*0.2)));
        V0M_10_15->SetFillColor(gStyle->GetColorPalette(int(nColors*0.2)));
        SPDT_10_15->SetLineColor(gStyle->GetColorPalette(int(nColors*0.2)));
        SPDT_10_15->SetFillColor(gStyle->GetColorPalette(int(nColors*0.2)));
        V0M_15_20->SetLineColor(gStyle->GetColorPalette(int(nColors*0.25)));
        V0M_15_20->SetFillColor(gStyle->GetColorPalette(int(nColors*0.25)));
        SPDT_15_20->SetLineColor(gStyle->GetColorPalette(int(nColors*0.25)));
        SPDT_15_20->SetFillColor(gStyle->GetColorPalette(int(nColors*0.25)));
        V0M_20_30->SetLineColor(gStyle->GetColorPalette(int(nColors*0.3)));
        V0M_20_30->SetFillColor(gStyle->GetColorPalette(int(nColors*0.3)));
        SPDT_20_30->SetLineColor(gStyle->GetColorPalette(int(nColors*0.3)));
        SPDT_20_30->SetFillColor(gStyle->GetColorPalette(int(nColors*0.3)));
        V0M_30_40->SetLineColor(gStyle->GetColorPalette(int(nColors*0.35)));
        V0M_30_40->SetFillColor(gStyle->GetColorPalette(int(nColors*0.35)));
        SPDT_30_40->SetLineColor(gStyle->GetColorPalette(int(nColors*0.35)));
        SPDT_30_40->SetFillColor(gStyle->GetColorPalette(int(nColors*0.35)));
        V0M_40_50->SetLineColor(gStyle->GetColorPalette(int(nColors*0.4)));
        V0M_40_50->SetFillColor(gStyle->GetColorPalette(int(nColors*0.4)));
        SPDT_40_50->SetLineColor(gStyle->GetColorPalette(int(nColors*0.4)));
        SPDT_40_50->SetFillColor(gStyle->GetColorPalette(int(nColors*0.4)));
        V0M_50_60->SetLineColor(gStyle->GetColorPalette(int(nColors*0.45)));
        V0M_50_60->SetFillColor(gStyle->GetColorPalette(int(nColors*0.45)));
        SPDT_50_60->SetLineColor(gStyle->GetColorPalette(int(nColors*0.45)));
        SPDT_50_60->SetFillColor(gStyle->GetColorPalette(int(nColors*0.45)));
        V0M_60_70->SetLineColor(gStyle->GetColorPalette(int(nColors*0.5)));
        V0M_60_70->SetFillColor(gStyle->GetColorPalette(int(nColors*0.5)));
        SPDT_60_70->SetLineColor(gStyle->GetColorPalette(int(nColors*0.5)));
        SPDT_60_70->SetFillColor(gStyle->GetColorPalette(int(nColors*0.5)));
        V0M_70_80->SetLineColor(gStyle->GetColorPalette(int(nColors*0.55)));
        V0M_70_80->SetFillColor(gStyle->GetColorPalette(int(nColors*0.55)));
        SPDT_70_80->SetLineColor(gStyle->GetColorPalette(int(nColors*0.55)));
        SPDT_70_80->SetFillColor(gStyle->GetColorPalette(int(nColors*0.55)));
        V0M_80_90->SetLineColor(gStyle->GetColorPalette(int(nColors*0.6)));
        V0M_80_90->SetFillColor(gStyle->GetColorPalette(int(nColors*0.6)));
        SPDT_80_90->SetLineColor(gStyle->GetColorPalette(int(nColors*0.6)));
        SPDT_80_90->SetFillColor(gStyle->GetColorPalette(int(nColors*0.6)));
        V0M_90_100->SetLineColor(gStyle->GetColorPalette(int(nColors*0.65)));
        V0M_90_100->SetFillColor(gStyle->GetColorPalette(int(nColors*0.65)));
        SPDT_90_100->SetLineColor(gStyle->GetColorPalette(int(nColors*0.65)));
        SPDT_90_100->SetFillColor(gStyle->GetColorPalette(int(nColors*0.65)));
        
               
        
        char hname[100];
        sprintf(hname, "V0M Distribution - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cV0M=new TCanvas(hname);
        V0M_0_100->SetTitle(hname);
        V0M_0_100->Draw();
        V0M_0_1->Draw("same");
//        V0M_1_3->Draw("same");
//        V0M_3_5->Draw("same");
//        V0M_5_10->Draw("same");
//        V0M_10_15->Draw("same");
//        V0M_15_20->Draw("same");
//        V0M_20_30->Draw("same");
//        V0M_30_40->Draw("same");
//        V0M_40_50->Draw("same");
//        V0M_50_60->Draw("same");
//        V0M_60_70->Draw("same");
//        V0M_70_80->Draw("same");
//        V0M_80_90->Draw("same");
//        V0M_90_100->Draw("same");
        
        sprintf(hname, "SPDT Distribution - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* cSPDT=new TCanvas(hname);
        SPDT_0_100->SetTitle(hname);
        SPDT_0_100->Draw();
        SPDT_0_1->Draw("same");
        SPDT_1_3->Draw("same");
        SPDT_3_5->Draw("same");
        SPDT_5_10->Draw("same");
        SPDT_10_15->Draw("same");
        SPDT_15_20->Draw("same");
        SPDT_20_30->Draw("same");
        SPDT_30_40->Draw("same");
        SPDT_40_50->Draw("same");
        SPDT_50_60->Draw("same");
        SPDT_60_70->Draw("same");
        SPDT_70_80->Draw("same");
        SPDT_80_90->Draw("same");
        SPDT_90_100->Draw("same");
        
        
        
    }
    
    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}

void PlotCorrelEstimators(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
   // Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
  //  Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
    Char_t *arrayOfPeriods[] = {"Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose"};//, "Group1_LHC17k", "Group4_LHC17k", "Group4_LHC18o"};
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
    
    TH2F* V0M_SPDT(NULL);
    TH2F* V0M_SPDC(NULL);
    TH2F* V0M_Mixed(NULL);
    TH2F* SPDT_SPDC(NULL);
    TH2F* SPDT_Mixed(NULL);
    TH2F* SPDC_Mixed(NULL);
    TH1F* Mixed_Dist(NULL);
    
    
    
    
    
    V0M_SPDT = new TH2F("V0M_SPDT",
                     "V0M_SPDT",
                     1000,0, 1000, 1000,0, 1000);
    V0M_SPDT->SetXTitle("V0M");
    V0M_SPDT->SetYTitle("SPDT");
    V0M_SPDC = new TH2F("V0M_SPDC",
                     "V0M_SPDC",
                     1000,0, 1000, 1000,0, 1000);
    V0M_SPDC->SetXTitle("V0M");
    V0M_SPDC->SetYTitle("SPDC");
    V0M_Mixed = new TH2F("V0M_Mixed",
                     "V0M_Mixed",
                     1000,0, 1000, 1000,0, 1000);
    V0M_Mixed->SetXTitle("V0M");
    V0M_Mixed->SetYTitle("Mixed");
    SPDT_SPDC = new TH2F("SPDT_SPDC",
                     "SPDT_SPDC",
                     1000,0, 1000, 1000,0, 1000);
    SPDT_SPDC->SetXTitle("SPDT");
    SPDT_SPDC->SetYTitle("SPDC");
    SPDT_Mixed = new TH2F("SPDT_Mixed",
                     "SPDT_Mixed",
                     1000,0, 1000, 1000,0, 1000);
    SPDT_Mixed->SetXTitle("SPDT");
    SPDT_Mixed->SetYTitle("Mixed");
    SPDC_Mixed = new TH2F("SPDC_Mixed",
                     "SPDC_Mixed",
                     1000,0, 1000, 1000,0, 1000);
    SPDC_Mixed->SetXTitle("SPDC");
    SPDC_Mixed->SetYTitle("Mixed");
    Mixed_Dist = new TH1F("Mixed_Dist",
                     "Mixed_Dist",
                     15000,0, 150);
    Mixed_Dist->SetXTitle("Mixed");
    Mixed_Dist->SetYTitle("Count");

    int neventos=0;
    double limitesmixed[100]{0};
    
    
   
// *************************
// Analyse                 *
// *************************
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
            
            sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
            
    // /Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT

    //sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
            
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
                    if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                            NumberCloseEtaTracklets++;
                    }
                    
                }
                
            //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
                if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
                    
                    double mixedvalue = ((4.3*(5.39/5.62)*((fEvent->fV0MValue/medianV0M16h)*10)+(2.*((fEvent->fSPDTrackletsValue/medianSPD16h)*10)))/6.3);
                    
                    Mixed_Dist->Fill(mixedvalue);
                    neventos++;
                   }
               
            }
                   TCanvas* c37=new TCanvas();
                   Mixed_Dist->Draw("");
    }
    
    cout << neventos<< " events seen"<<endl;
    
    int quantile = 1;
    int runningsum = 0;
    for(int binx = Mixed_Dist->GetNbinsX(); binx>=0; binx--){
        
        runningsum+=Mixed_Dist->GetBinContent(binx+1);
        cout << "binx: "<<binx<<" runningsum: "<<runningsum<<endl;
        if (runningsum > neventos*(quantile/100.)){
            cout << "above " << neventos*(quantile/100.)<<endl;
            limitesmixed[quantile-1] = (binx+1)*150/15000.; // limitmixed[x] is the last bin in which the amount of estimator was strictly less than x% of the distrib
            cout << "limites "<< quantile-1<< " fixed at "<<(binx+1)*150/15000.<<endl;
            quantile++;
        }
    }
    
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        
// /Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT

//sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
        
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
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ //&& fEvent->fNPileupVtx == 0
                
                double mixedvalue = ((4.3*(5.39/5.62)*((fEvent->fV0MValue/medianV0M16h)*10)+(2.*((fEvent->fSPDTrackletsValue/medianSPD16h)*10)))/6.3);
                int mixedperc = 0;
                for(int quantos=0; quantos<100; quantos++){
                    if(mixedvalue > limitesmixed[quantos]){
                        mixedperc = quantos+1;
                        break;
                    }
                }
                
                V0M_SPDT->Fill(fEvent->fCentralityV0M, fEvent->fCentralitySPDTracklets);
                V0M_SPDC->Fill(fEvent->fCentralityV0M, fEvent->fCentralitySPDClusters);
                V0M_Mixed->Fill(fEvent->fCentralityV0M, mixedperc);
                SPDT_SPDC->Fill(fEvent->fCentralitySPDTracklets, fEvent->fCentralitySPDClusters);
                SPDT_Mixed->Fill(fEvent->fCentralitySPDTracklets, mixedperc);
                SPDC_Mixed->Fill(fEvent->fCentralitySPDClusters, mixedperc);

               }
           
        }
        
               
        
        char hname[100];
        sprintf(hname, "V0M vs SPDT - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* c1=new TCanvas(hname);
        V0M_SPDT->SetTitle(hname);
        V0M_SPDT->Draw("colz");
       
        sprintf(hname, "V0M vs SPDC - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* c2=new TCanvas(hname);
        V0M_SPDC->SetTitle(hname);
        V0M_SPDC->Draw("colz");
        
        sprintf(hname, "V0M vs Mixed - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* c3=new TCanvas(hname);
        V0M_Mixed->SetTitle(hname);
        V0M_Mixed->Draw("colz");
        
        sprintf(hname, "SPDT vs SPDC - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* c4=new TCanvas(hname);
        SPDT_SPDC->SetTitle(hname);
        SPDT_SPDC->Draw("colz");
        
        sprintf(hname, "SPDT vs Mixed - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* c5=new TCanvas(hname);
        SPDT_Mixed->SetTitle(hname);
        SPDT_Mixed->Draw("colz");
        
        sprintf(hname, "SPDC vs Mixed - %s" ,arrayOfPeriods[tree_idx]);
        TCanvas* c6=new TCanvas(hname);
        SPDC_Mixed->SetTitle(hname);
        SPDC_Mixed->Draw("colz");
        
        
        
    }
    
    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}

void PlotPolar(){

    auto c = new TCanvas("myCanvas","myCanvas",1600,400);
    c->Divide(6,1);
    Double_t rmin=0.;
    Double_t rmax=TMath::Pi()*2.;
    const Int_t npoints=1000;
    Double_t r[npoints];
    Double_t theta[npoints];
    
    gStyle->SetPalette(kRainBow);
    int nColors = gStyle->GetNumberOfColors();
    
    
    for(int order=0;order<=6;order++){
    for (Int_t ipt = 0; ipt < npoints; ipt++) {
        r[ipt] = ipt*(rmax-rmin)/npoints + rmin;
        theta[ipt] = TMath::Cos(order*r[ipt]);
    }
    c->cd(order+1);
        char hname[20];
        sprintf(hname, "n=%i", order);
    TGraphPolar grP1 (npoints,r,theta);
    grP1.SetTitle(hname);
    grP1.SetLineWidth(6);
    grP1.SetLineColor(gStyle->GetColorPalette(int(nColors-(nColors*(order)/6.))));
   // grP1.SetLineColor(2);
        if(order>0){
    grP1.DrawClone("L");
        }
        c->Update();
    }

}

void CorrelEstMult(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
   // Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
  //  Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
    Char_t *arrayOfPeriods[] = {"Group1_LHC16hMCno2Chg"};//
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
    
    TH2F* V0M_Nch(NULL);
    TH2F* SPDT_Nch(NULL);
    TH2F* SPDC_Nch(NULL);
    
    
    
    
    
    V0M_Nch = new TH2F("V0M_Nch",
                     "V0M_Nch",
                     1000,0, 1000, 1000,0, 1000);
    V0M_Nch->SetXTitle("V0M");
    V0M_Nch->SetYTitle("Nch");
    SPDT_Nch = new TH2F("SPDT_Nch",
                     "SPDT_Nch",
                     1000,0, 1000, 1000,0, 1000);
    SPDT_Nch->SetXTitle("SPDT");
    SPDT_Nch->SetYTitle("Nch");
    SPDC_Nch = new TH2F("SPDC_Nch",
                     "SPDC_Nch",
                     1000,0, 1000, 1000,0, 1000);
    SPDC_Nch->SetXTitle("SPDC");
    SPDC_Nch->SetYTitle("Nch");
    
    int neventos=0;
    double limitesmixed[100]{0};
    
    
   
// *************************
// Analyse                 *
// *************************
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
            
            sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
            
    // /Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT

    //sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
            
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
        
        TClonesArray *fMCDimuons = 0;
        MCDimuonLight *mcdimuon = 0;
        
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
        fMCDimuons = fEvent->fMCDimuons;
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
                
                if(i%10!=0){
                    continue;
                }
                
                theTree->GetEvent(i);
         //       cout << "Looking at Event " << i << endl;
    //            cout << "fPassPhysicsSelection : " << fEvent->fPassPhysicsSelection << endl;
    //            cout << "fPassTriggerSelection : " << fEvent->fPassTriggerSelection << endl;
    //            cout << "fTriggerCMUL7 : " << fEvent->fTriggerCMUL7 << endl;
    //            cout << "fTriggerCINT7 : " << fEvent->fTriggerCINT7 << endl;
                //cout << "tracks->GetEntries() " << tracks->GetEntries() << endl;
              //  cout << "Computing the number of tracklets in Event " << i << endl;
                
             int NumberAcceptedEtaTracklets = 0;
              int NumberOfTrackletsForCentrality = 0;
              for (Int_t j=0; j<fEvent->fNMCDimuons; j++) {
                  mcdimuon = (MCDimuonLight*)fMCDimuons->At(j);
                  if((TMath::Abs(mcdimuon->fEta) < 2.5) && (mcdimuon->fPt > 0.4) && (mcdimuon->fCharge != 0) && (TMath::Abs(mcdimuon->fPDG)==211 || TMath::Abs(mcdimuon->fPDG)==321 || TMath::Abs(mcdimuon->fPDG)==2212)){
                      NumberOfTrackletsForCentrality++;
                  }
                  if((mcdimuon->fEta < myEtaHigh(fEvent->fVertexZ)) && (mcdimuon->fEta > myEtaLow(fEvent->fVertexZ)) && (mcdimuon->fPt > 0.5) && (mcdimuon->fPt < 5.0) && (mcdimuon->fCharge != 0) && (TMath::Abs(mcdimuon->fPDG)==211 || TMath::Abs(mcdimuon->fPDG)==321 || TMath::Abs(mcdimuon->fPDG)==2212)){
                      NumberAcceptedEtaTracklets++;
                  }
                  
              }
                
                V0M_Nch->Fill(fEvent->fV0MValue,NumberOfTrackletsForCentrality);
                SPDT_Nch->Fill(fEvent->fSPDTrackletsValue,NumberOfTrackletsForCentrality);
                SPDC_Nch->Fill(fEvent->fSPDClustersValue,NumberOfTrackletsForCentrality);
                
    }
        
        
    }

    
    TCanvas* cV0M_Nch=new TCanvas();
    cV0M_Nch->cd();
    V0M_Nch->Draw("colz");

    
    TCanvas* cSPDT_Nch=new TCanvas();
    cSPDT_Nch->cd();
    SPDT_Nch->Draw("colz");

    
    TCanvas* cSPDC_Nch=new TCanvas();
    cSPDC_Nch->cd();
    SPDC_Nch->Draw("colz");

    
 
    
    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}

void PlotPHD(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
   // Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
    Char_t *arrayOfPeriods[] = {"Group4_LHC18o"};
   // Char_t *arrayOfPeriods[] = {"Group1_LHC16hMCno2Chg"};//
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
    
    TH2F* hnsegSigma(NULL);
    
    
    
    
    
    hnsegSigma = new TH2F("hnsegSigma",
                      "z_{vtx} resolution wrt number of contributors",
                      60,0,60,1000,0,10);
    hnsegSigma->SetXTitle("Contributors");
    hnsegSigma->SetYTitle("z_{vtx} Resolution (SPD)");
    
    int neventos=0;
    double limitesmixed[100]{0};
    
    
   
// *************************
// Analyse                 *
// *************************
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
            
            sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
            
    // /Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT

    //sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
            
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
        
        TClonesArray *fMCDimuons = 0;
        MCDimuonLight *mcdimuon = 0;
        
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
        fMCDimuons = fEvent->fMCDimuons;
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
                
//                if(i%10!=0){
//                    continue;
//                }
                
                theTree->GetEvent(i);
         //       cout << "Looking at Event " << i << endl;
    //            cout << "fPassPhysicsSelection : " << fEvent->fPassPhysicsSelection << endl;
    //            cout << "fPassTriggerSelection : " << fEvent->fPassTriggerSelection << endl;
    //            cout << "fTriggerCMUL7 : " << fEvent->fTriggerCMUL7 << endl;
    //            cout << "fTriggerCINT7 : " << fEvent->fTriggerCINT7 << endl;
                //cout << "tracks->GetEntries() " << tracks->GetEntries() << endl;
              //  cout << "Computing the number of tracklets in Event " << i << endl;
                
             
        
    hnsegSigma->Fill(fEvent->fVertexNC, fEvent->fSPDVertexSigmaZ); //SIGMA

    
            }

    }
 
    TCanvas* chnsegsigma = new TCanvas;
       chnsegsigma->cd();
       hnsegSigma->Draw("colz");
    
    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}

void CheckPYTHIA(){
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
    
    double meanV0M16h = 107.5;
    double meanSPDTkl16h = 16.5;
    double entries16h = 680844;
    double medianV0M16h = 82;
    double maxV0M16h = 868;
    double medianSPD16h = 12;
    double maxSPD16h = 139;

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
    
   // Char_t Group_Period[1000] = "Manu16hCINT_PhySelTRUE_SelNtkl1_Vertexer_CorrelVeryLoose";
  //  Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
    Char_t *arrayOfPeriods[] = {"Group1_LHC16hMCnospe", "Group1_LHC16hMCno2Chg"};//
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
    
    TH1F* Nch(NULL);
    TH1F* Nchexplicit(NULL);
    TH1F* PDG(NULL);
    TH1F* Charge(NULL);
    TH1F* Pt(NULL);
    TH1F* Eta(NULL);
    
    
    
    
    Nch = new TH1F("Nch",
                     "Nch",
                     1000,0, 1000);
    Nch->SetXTitle("Nch");
    Nch->SetYTitle("Count");
    
    Nchexplicit = new TH1F("Nchexplicit",
                     "Nchexplicit",
                     1000,0, 1000);
    Nchexplicit->SetXTitle("Nch");
    Nchexplicit->SetYTitle("Count");

PDG = new TH1F("PDG",
                     "PDG",
                     10000,-5000, 5000);
    PDG->SetXTitle("PDG");
    PDG->SetYTitle("Count");
    
    Charge = new TH1F("Charge",
                     "Charge",
                     10,-5, 5);
    Charge->SetXTitle("Charge");
    Charge->SetYTitle("Count");
    
    Eta = new TH1F("Eta",
                     "Eta",
                     1000,-5, 5);
    Eta->SetXTitle("Eta");
    Eta->SetYTitle("Count");
    
    Pt = new TH1F("Pt",
                     "Pt",
                     1000,0, 20);
    Pt->SetXTitle("Pt");
    Pt->SetYTitle("Count");
    
    int neventos=0;
    double limitesmixed[100]{0};
    int countery = 0;
    
    
   
// *************************
// Analyse                 *
// *************************
    
    TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    TCanvas* cNch=new TCanvas();
TCanvas* cPDG=new TCanvas();
    TCanvas* cCharge=new TCanvas();
    TCanvas* cPt=new TCanvas();
    TCanvas* cEta=new TCanvas();

    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
            
            sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
            
    // /Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT

    //sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
            
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
        
        TClonesArray *fMCDimuons = 0;
        MCDimuonLight *mcdimuon = 0;
        
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
        fMCDimuons = fEvent->fMCDimuons;
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
                
                if(i%10!=0){
                    continue;
                }
                
                theTree->GetEvent(i);
         //       cout << "Looking at Event " << i << endl;
    //            cout << "fPassPhysicsSelection : " << fEvent->fPassPhysicsSelection << endl;
    //            cout << "fPassTriggerSelection : " << fEvent->fPassTriggerSelection << endl;
    //            cout << "fTriggerCMUL7 : " << fEvent->fTriggerCMUL7 << endl;
    //            cout << "fTriggerCINT7 : " << fEvent->fTriggerCINT7 << endl;
                //cout << "tracks->GetEntries() " << tracks->GetEntries() << endl;
              //  cout << "Computing the number of tracklets in Event " << i << endl;
                
             int NumberAcceptedEtaTracklets = 0;
              int NumberOfTrackletsForCentrality = 0;
                int NumberOfTrackletsForCentralityExplicit = 0;
              for (Int_t j=0; j<fEvent->fNMCDimuons; j++) {
                  mcdimuon = (MCDimuonLight*)fMCDimuons->At(j);
//                  if((TMath::Abs(mcdimuon->fEta) < 2.5) && (mcdimuon->fPt > 0.4) && (mcdimuon->fCharge != 0) && (TMath::Abs(mcdimuon->fPDG)==211 || TMath::Abs(mcdimuon->fPDG)==321 || TMath::Abs(mcdimuon->fPDG)==2212)){
//                      NumberOfTrackletsForCentrality++;
//                  }
//                  if((mcdimuon->fEta < myEtaHigh(fEvent->fVertexZ)) && (mcdimuon->fEta > myEtaLow(fEvent->fVertexZ)) && (mcdimuon->fPt > 0.5) && (mcdimuon->fPt < 5.0) && (mcdimuon->fCharge != 0) && (TMath::Abs(mcdimuon->fPDG)==211 || TMath::Abs(mcdimuon->fPDG)==321 || TMath::Abs(mcdimuon->fPDG)==2212)){
//                      NumberAcceptedEtaTracklets++;
//                  }
                  if((TMath::Abs(mcdimuon->fEta) < 2.5) && (mcdimuon->fPt > 0.4) && (mcdimuon->fCharge != 0)){
                                       NumberOfTrackletsForCentrality++;
                    PDG->Fill(mcdimuon->fPDG);
                      Charge->Fill(mcdimuon->fCharge);
                      Pt->Fill(mcdimuon->fPt);
                      Eta->Fill(mcdimuon->fEta);
  if((TMath::Abs(mcdimuon->fPDG)==211 || TMath::Abs(mcdimuon->fPDG)==321 || TMath::Abs(mcdimuon->fPDG)==2212));
                          NumberOfTrackletsForCentralityExplicit++;
                      }
                  
                  if(NumberOfTrackletsForCentrality != NumberOfTrackletsForCentralityExplicit){
                      cout << "DIFFERENCE"<<endl;
                  }
                  
                  Nch->Fill(NumberOfTrackletsForCentrality);
                 // Nchexplicit->Fill(NumberOfTrackletsForCentralityExplicit);
                  
              }
    }
        
Nch->Scale(1./nevent);
        PDG->Scale(1./nevent);
        Charge->Scale(1./nevent);
        Pt->Scale(1./nevent);
        Eta->Scale(1./nevent);
        
        
        
        if(countery==0){
            cNch->cd();
            Nch->DrawCopy();
            cPDG->cd();
            PDG->DrawCopy();
            cCharge->cd();
            Charge->DrawCopy();
            cPt->cd();
            Pt->DrawCopy();
            cEta->cd();
            Eta->DrawCopy();
            countery++;
        }
        
        else if(countery>0){
        cNch->cd();
        Nch->SetLineColor(kRed);
        Nch->Draw("same");
            cPDG->cd();
            PDG->SetLineColor(kRed);
            PDG->Draw("same");
            cCharge->cd();
            Charge->SetLineColor(kRed);
            Charge->Draw("same");
            cPt->cd();
            Pt->SetLineColor(kRed);
            Pt->Draw("same");
            cEta->cd();
            Eta->SetLineColor(kRed);
            Eta->Draw("same");
        }
       Nch->Clear();
        PDG->Clear();
        Charge->Clear();
        Pt->Clear();
        Eta->Clear();
        
    }

 
    
    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    

}
