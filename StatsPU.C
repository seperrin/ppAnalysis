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

//Code pour regarder les statistiques et les effets des différents cuts (notamment pour comprendre l'effetd es PU cuts)


// FUNCTIONS

void StatsPU(Char_t* GroupeChoisi);
void MultipleStatsPU();

void MultipleStatsPU(){
     Char_t *arrayPeriods[] = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l","Group2_LHC17h","Group3_LHC17h","Group4_LHC17k","Group4_LHC18l","Group4_LHC18m","Group4_LHC18o","Group4_LHC18p","Group5_LHC17l","Group5_LHC17m","Group5_LHC17o","Group5_LHC17r","Group5_LHC18c","Group5_LHC18d","Group5_LHC18e","Group5_LHC18f","Group6_LHC18m","Group7_LHC18m","Group8_LHC18m","Group9_LHC18m","Group10_LHC18m","Group11_LHC18m","Group12_LHC18m"};
    int nOfPeriods = sizeof(arrayPeriods) / sizeof(arrayPeriods[0]);
    
    for(int idx=0; idx<nOfPeriods; idx++){
        StatsPU(arrayPeriods[idx]);
    }
}

void StatsPU(Char_t* GroupeChoisi){
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
    Char_t *arrayOfPeriods[] = {GroupeChoisi};
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
    
    
    TH1I* hDistDimuons(NULL);
    TH1I* hDistCloseTkl(NULL);
    TH1I* hDistVertexNC(NULL);
    TH1I* hDistNPileupVtx(NULL);
    TH1I* hDistIsPileupFromSPDMultBins(NULL);
    TH1I* hDistSigmaZ(NULL);
    TH1I* hDistZ(NULL);
    TH2I* hDistNclustNspd(NULL);
    TH2I* hDistNclustNspdCut(NULL);
    TH1I* hPhysicsSelection(NULL);
    
    
    
    hDistDimuons = new TH1I("hDistDimuons",
                        "Number of dimuons per event",
                        5,-0.5,4.5);
       hDistDimuons->SetXTitle("Nunber of dimuons");
       hDistDimuons->SetYTitle("Count");
    
    hDistCloseTkl = new TH1I("hDistCloseTkl",
                     "Number of close tracklets per event",
                     120,-0.5,119.5);
    hDistCloseTkl->SetXTitle("Number of close tracklets");
    hDistCloseTkl->SetYTitle("Count");
    
    hDistVertexNC = new TH1I("hDistVertexNC",
                     "Number of VertexNC per event",
                     120,-20.5,99.5);
    hDistVertexNC->SetXTitle("Number of VertexNC");
    hDistVertexNC->SetYTitle("Count");
    
    hDistNPileupVtx = new TH1I("hDistNPileupVtx",
                     "Number of NPileupVtx per event",
                     140,-20.5,119.5);
    hDistNPileupVtx->SetXTitle("Number of NPileupVtx");
    hDistNPileupVtx->SetYTitle("Count");
    
    hDistIsPileupFromSPDMultBins = new TH1I("hDistIsPileupFromSPDMultBins",
                     "IsPileupFromSPDMultBins per event",
                     2,-0.5,1.5);
    hDistIsPileupFromSPDMultBins->SetXTitle("IsPileupFromSPDMultBins");
    hDistIsPileupFromSPDMultBins->SetYTitle("Count");
    
    hDistSigmaZ = new TH1I("hDistSigmaZ",
                     "SigmaZ per event",
                     1000,0,10);
    hDistSigmaZ->SetXTitle("SigmaZ");
    hDistSigmaZ->SetYTitle("Count");
    
    hDistZ = new TH1I("hDistZ",
                        "Z per event",
                        1000,-20,20);
       hDistZ->SetXTitle("Z");
       hDistZ->SetYTitle("Count");
    
    hDistNclustNspd= new TH2I("hDistNclustNspd",
                     "Distribution of Nclusters VS SpdTracklets",
                     100,0,100, 500, 0, 500);
    hDistNclustNspd->SetXTitle("SpdTracklets");
    hDistNclustNspd->SetYTitle("Nclusters");
    
    hDistNclustNspdCut= new TH2I("hDistNclustNspdCut",
                     "Distribution of Nclusters VS SpdTracklets - Cut",
                     100,0,100, 500, 0, 500);
    hDistNclustNspdCut->SetXTitle("SpdTracklets");
    hDistNclustNspdCut->SetYTitle("Nclusters");
    
    hPhysicsSelection = new TH1I("hPhysicsSelection",
                     "PhysicsSelection per event",
                     2,-0.5,1.5);
    hPhysicsSelection->SetXTitle("PhysicsSelection");
    hPhysicsSelection->SetYTitle("Count");

// *************************
// Analyse                 *
// *************************
    
    int RejectedBySigmaZ = 0;
    int RejectedByZ = 0;
    int SafeEventsWOPU = 0;
    int SafeEventsRemovedByPU = 0;
    int ManuelCount = 0;
    int RejectedByTklNoMrad = 0;
    int RejectedByNC = 0;
    int RejectedByTklMrad = 0;
    int EventCountMrad = 0;
    int HarshCount = 0;
    int EventsInTTree = 0;
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        //sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
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
        EventsInTTree = nevent;
    
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
            
            // Event counter à la Manuel
            
            int NumberCloseEtaTrackletsManuel = 0;
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                if((TMath::Abs(trac->fEta) < 1)){
                        NumberCloseEtaTrackletsManuel++;
                }
                
            }
            
            if(NumberCloseEtaTrackletsManuel ==0){
                RejectedByTklNoMrad++;
            }
            
            if( fEvent->fVertexNC <= 0){
                RejectedByNC++;
            }
            
            if(!fEvent->fIsPileupFromSPDMultBins && fEvent->fSPDVertexSigmaZ<=0.25 && TMath::Abs(fEvent->fVertexZ)<10 && NumberCloseEtaTrackletsManuel != 0 && fEvent->fVertexNC > 0){
                ManuelCount++;
            }
            
            if(!fEvent->fIsPileupFromSPDMultBins && fEvent->fSPDVertexSigmaZ<=0.25 && TMath::Abs(fEvent->fVertexZ)<10 && NumberCloseEtaTrackletsManuel != 0 && fEvent->fVertexNC > 0 && fEvent->fNPileupVtx == 0){
                HarshCount++;
            }
            
            
            
            hDistDimuons->Fill(fEvent->fNDimuons);
            
            int NumberCloseEtaTracklets = 0;
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                }
                
            }
            
            hPhysicsSelection->Fill(fEvent->fPassPhysicsSelection);
            hDistCloseTkl->Fill(NumberCloseEtaTracklets);
            hDistVertexNC->Fill(fEvent->fVertexNC);
             hDistNPileupVtx->Fill(fEvent->fNPileupVtx);
            hDistIsPileupFromSPDMultBins->Fill(fEvent->fIsPileupFromSPDMultBins);
            hDistSigmaZ->Fill(fEvent->fSPDVertexSigmaZ);
            hDistZ->Fill(fEvent->fVertexZ);
            
            if(TMath::Abs(fEvent->fVertexZ)>=10){
                RejectedByZ++;
            }
            
            if(fEvent->fSPDVertexSigmaZ>0.25){
                           RejectedBySigmaZ++;
                }
            
            if(NumberCloseEtaTracklets == 0){
                RejectedByTklMrad++;
            }
            
            if(!fEvent->fIsPileupFromSPDMultBins && fEvent->fSPDVertexSigmaZ<=0.25 && TMath::Abs(fEvent->fVertexZ)<10 && NumberCloseEtaTracklets != 0 && fEvent->fVertexNC > 0){
                EventCountMrad++;
            }
            
            if(NumberCloseEtaTracklets>0 && TMath::Abs(fEvent->fVertexZ)<10 && fEvent->fSPDVertexSigmaZ<=0.25){
                SafeEventsWOPU++;
                hDistNclustNspd->Fill(fEvent->fSPDTrackletsValue, fEvent->fSPDClustersValue);
                if(fEvent->fSPDClustersValue <= (100 + 4*(fEvent->fSPDTrackletsValue)) && fEvent->fVertexNC > 0 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0){
                    hDistNclustNspdCut->Fill(fEvent->fSPDTrackletsValue, fEvent->fSPDClustersValue);
                }
                if(fEvent->fVertexNC < 1 || fEvent->fNPileupVtx != 0 || fEvent->fIsPileupFromSPDMultBins != 0 ){
                    SafeEventsRemovedByPU++;
                }
            }
            
            
            
            
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC > 0 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){

                
                
                    for (Int_t j=0; j<fEvent->fNDimuons; j++) {
                    dimu = (DimuonLight*)fDimuons->At(j);
                        if((dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut) && (dimu->fInvMass > MinInvMass) && (dimu->fInvMass < MaxInvMass)){
                        }
                    }
                    
                
                
                
                
               }
            
            
        }
        
        cout << "Total number of events in the TTree: " << nevent<<endl;
    }

    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    
    
    
// ***********************************
// Tracer les graphes sur les canevas *
// ***********************************
    
    TCanvas*c1=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c1->SetTitle("Distribution of dimuons per event");
    hDistDimuons->Draw();
    
    TCanvas*c2=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c2->SetTitle("Distribution of close tracklets per event");
    hDistCloseTkl->Draw();
    
    cout << hDistCloseTkl->GetBinContent(1) << " events have 0 close tracklets, will be discarded"<<endl;
    
    TCanvas*c6=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c6->SetTitle("Distribution of SigmaZ per event");
    hDistSigmaZ->Draw();
    
    cout << RejectedBySigmaZ << " events SigmaZ > 0.25, will be discarded"<<endl;
    
    TCanvas*c7=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c7->SetTitle("Distribution of Z per event");
    hDistZ->Draw();
    
    cout << RejectedByZ << " events Abs(Z) >= 10, will be discarded"<<endl;
    
    
    
    cout << "There are " << SafeEventsWOPU << " events passing cuts on Z, SigmaZ, and Close Tracklets. Events safe without PU cuts for now."<<endl;
    
    cout <<"Out of all events..."<<endl;
    TCanvas*c3=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c3->SetTitle("Distribution of VertexNC per event");
    hDistVertexNC->Draw();
    
    int rejectVNC = 0;
    for(int binno=1; binno<22; binno++){
        rejectVNC += hDistVertexNC->GetBinContent(binno);
    }
    
    cout << rejectVNC << " events have 0 vertex contributors, will be discarded bu PU cuts"<<endl;
    
    TCanvas*c4=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c4->SetTitle("Distribution of NPileupVtx per event");
    hDistNPileupVtx->Draw();
    
    cout << hDistNPileupVtx->GetEntries()-hDistNPileupVtx->GetBinContent(21) << " events have PU vertices, will be discarded bu PU cuts"<<endl;
    
    TCanvas*c5=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c5->SetTitle("Distribution of IsPileupFromSPDMultBins per event");
    hDistIsPileupFromSPDMultBins->Draw();
    
    cout << hDistIsPileupFromSPDMultBins->GetBinContent(2) << " events have PU from multiplicity bins, will be discarded bu PU cuts"<<endl;
    
    cout <<"Out of the " << SafeEventsWOPU << " events passing non-PU cuts, " << SafeEventsRemovedByPU << " have been discarded by PU cuts" <<endl;
    cout <<"That represents " << double((SafeEventsRemovedByPU*100/SafeEventsWOPU)) << "% of events remobed solely by PU cuts"<<endl;
    
    
    TCanvas*c8=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c8->SetTitle("Distribution of NClust vs NSPDTkl per event");
    hDistNclustNspd->Draw("colz");
    
    TCanvas*c9=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c9->SetTitle("Distribution of NClust vs NSPDTkl per event - Cut Nclusters <= 100 + 4Ntracklets");
    hDistNclustNspdCut->Draw("colz");
    
    TCanvas*c10=new TCanvas();
      // ROOT THAT Extraction 2 plot
      c10->SetTitle("PhysicsSelection Distribution");
    hPhysicsSelection->Draw();
    
    cout << "Applying Manuel's Analysis Cuts, we arrive at: " << ManuelCount << "#CMUL7 events"<<endl;
    
    cout << "+++++++++++++++++++++++++++++++++" <<endl;
    
    cout << "Group: " << arrayOfPeriods[0] <<endl;
    cout << "Number of events in TTree: " << EventsInTTree <<endl;
    cout << "Number of events not passing PU mult cut: " << hDistIsPileupFromSPDMultBins->GetBinContent(2)<<endl;
    cout << "Number of events not having at least a vertex contributor: " << RejectedByNC<<endl;
    cout << "Number of events not passing SigmaZ cut: "<< RejectedBySigmaZ<<endl;
    cout << "Number of events not passing Z cut: " << RejectedByZ<<endl;
    cout << "Number of events not having valid tracklets (Eta acceptance): " << RejectedByTklNoMrad<<endl;
    cout << "Valid CMUL7 events [Manuel-like selection]: " << ManuelCount<<endl;
    cout << "Number of events not having valid tracklets (Eta acceptance and 10 mrad cut): " <<RejectedByTklMrad<<endl;
    cout << "Valid CMUL7 events [10 mrad additional Tkl cut]: " << EventCountMrad<<endl;
    cout << "Number of events removed from harsh PU cut :"<< long(hDistNPileupVtx->GetEntries()-hDistNPileupVtx->GetBinContent(21))<<endl;
    cout << "Valid CMUL7 events [Harsh PU, no mrad]: "<<HarshCount<<endl;
    
    
}


