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
 #include <TGraphErrors.h>
 #include <TMatrixD.h>
 #include <TLorentzVector.h>
 #include <string>
 #include <vector>

 //#include "TTreeReader.h"
 //#include "Event.h"

 # include "TF1.h"
 # include "TF2.h"
 # include "TProfile.h"
 # include "TVirtualFitter.h"
 # include "TBackCompFitter.h"
 # include "TGraphErrors.h"
 # include "TFitter.h"
 # include "TMinuit.h"
 # include "TRandom.h"
 # include <iostream>
 # include <string>
 # include "TH1F.h"
 # include "TH1D.h"
 # include "TH2F.h"
 # include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"


 // FUNCTIONS

 void FitExt2UncA();
 Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par);
 Double_t BackFcnV2(Double_t *x,Double_t *par);
 Double_t SignalFcnJPsiV2(Double_t *x,Double_t *par);
 Double_t FourierV2(Double_t *x,Double_t *par);
 Double_t FourierV5(Double_t *x,Double_t *par);
 Double_t TwoCBE2E(Double_t *x,Double_t *par);
 Double_t ExpBkg(Double_t *x,Double_t *par);
 Double_t JPsiCrystalBallExtended(Double_t *x,Double_t *par);
 Double_t Psi2SCrystalBallExtended(Double_t *x,Double_t *par);
 TFitResultPtr FittingAllInvMass(const char *histoname, TCanvas *canvas);
 TFitResultPtr FittingAllInvMassBin(const char *histoname, TCanvas *canvas, int i);
 void FcnCombinedAllMass(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  );
 int GetCent(double cent);

Double_t mJpsi =  3.096916;
 Double_t mPsip =  3.686108;
// Double_t ratMass = 1.01;
 Double_t ratSigma = 1.05;
Int_t npfits;

TH1F* hnseg(NULL);
TH1F* V2JPsiTkl(NULL);

Char_t FitFileName[200];

void FitExt2UncA(){
    
    sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFilesRun2.root");
    
    TH1::SetDefaultSumw2();
        bool doTracklets = kFALSE;
        bool doMixedEvents = kFALSE;
    bool Extraction2Combined = kTRUE;
        bool CombineFits = kTRUE;
        
    // ************************************
    // Définitions de paramètres          *
    // ************************************

        double ZvtxCut = 10;
        double SigmaZvtxCut = 0.25;
        double DPhiCut = 0.01;
        double TklEtaCut  = 1;
        double LowDimuYCut = -4;
        double HighDimuYCut = -2.5;
        double LowDimuPtCut = 0;
        double HighDimuPtCut = 99999;
    //    double MinMultCentral = 37;
    //    double MaxMultPeriph = 23;
        double CentSPDTrackletsCentral = 1;
        double CentSPDTrackletsPeriph = 40;

        double LowJPsiMass = 3.0;
        double HighJPsiMass = 3.3;

        double MinInvMass = 2.1;
        double MaxInvMass = 5.1;
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
        
        double MinDeltaEtaTKL = 1.2;
        double MaxDeltaEtaTKL = 2.4;
        const int NbinsDeltaEtaTKL = 24;
        double SizeBinDeltaEtaTKL = (MaxDeltaEtaTKL-MinDeltaEtaTKL)/NbinsDeltaEtaTKL;
        
        const int NbinsDeltaPhiTKL = 24;
        double SizeBinDeltaPhiTKL = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhiTKL;
        
        const int NbBinsCent = 5;
        const int NbBinsZvtx = 20;
    
    // *************************
    // Initialiser les graphes *
    // *************************
    
    TH1F* YieldValue(NULL);
    TH1F* YieldError(NULL);
    
    YieldValue = new TH1F("YieldValue",
    "YieldValue",
    10000,0,100);
    YieldError = new TH1F("YieldError",
    "YieldError",
    10000,0,100);
        
      //  TH1F* hnseg(NULL);
        TH1F* hnseg2(NULL);
        TH1F* hnseg3(NULL);
        TH1F* hnseg4(NULL);
        TH1F* hnseg5(NULL);
        TH1F* hnseg6(NULL);
        TH1F* hnseg7(NULL);
        TH2F* hnseg8(NULL);
        TH2F* hnsegSigma(NULL);
        TH1F* hnseg8_proj_tampon(NULL);
        TH1F* hnseg8Tkl_proj_tampon(NULL);
        TH1F* hnseg8_proj(NULL);
        TH1F* hnseg8Tkl_proj(NULL);
        TH2F* hnseg8ME(NULL);
        TH1F* hnseg8ME_proj_tampon(NULL);
        TH1F* hnseg8TklME_proj_tampon(NULL);
        TH1D* hnseg8ME_proj(NULL);
        TH1D* hnseg8TklME_proj(NULL);
        TH1F* ME_proj_tampon(NULL);
        TH1F* ProjCopy(NULL);
        TH1F* ProjCopy2(NULL);
        TH1F* ProjCopyTkl(NULL);
        TH1F* ProjCopy2Tkl(NULL);
        TH1F* ME_proj_Tkl_tampon(NULL);
        TH1F* ME_proj(NULL);
        TH1F* SE_proj_tampon(NULL);
        TH1F* SE_proj_Tkl_tampon(NULL);
        TH1F* SE_proj(NULL);
        
        TH1F* Sijk(NULL);
        TH1F* Mijk(NULL);
        TH1F* SoverMijk(NULL);
        TH1F* SoverMik(NULL);
        
        TH1F* Sij(NULL);
        TH1F* Mij(NULL);
        TH1F* SoverMij(NULL);
        TH1F* SoverMi(NULL);
        
        TH1F* Yields[NbBinsCent][NbinsInvMass]{ NULL };
        TH1F* YieldsTkl[NbBinsCent]{ NULL };
        TH2F* Correlations[NbBinsCent][NbBinsZvtx][NbinsInvMass]{ NULL };
        TH2F* CorrelationsTkl[NbBinsCent][NbBinsZvtx]{ NULL };
        TH2F* CorrelationsME[NbBinsCent][NbBinsZvtx][NbinsInvMass]{ NULL };
        TH2F* CorrelationsTklME[NbBinsCent][NbBinsZvtx]{ NULL };
        TH1F* Yield_tampon(NULL);
        TH1F* YieldTkl_tampon(NULL);
        TH1F* YieldWrtMass_tampon(NULL);
        TH1F* Yield_allC(NULL);
        TH1F* Yield_Central(NULL);
        TH1F* Yield_Periph(NULL);
        TH1F* Yield_Difference(NULL);
        TH1F* YieldTkl_allC(NULL);
        TH1F* YieldTkl_Central(NULL);
        TH1F* YieldTkl_Periph(NULL);
        TH1F* YieldTkl_Difference(NULL);
        TH1F* Yield_Central_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yield_Periph_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yield_Difference_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yields_PhiBin[NbBinsCent][NbinsDeltaPhi]{ NULL };
        TH1F* YieldWrtMass_allC[NbinsDeltaPhi]{ NULL };
        TH1F* YieldWrtMass_Central[NbinsDeltaPhi]{ NULL };
        TH1F* YieldWrtMass_Periph[NbinsDeltaPhi]{ NULL };
        
        // Try plots 2D DetaDphi
        
        TH2F* YCentral(NULL);
        TH2F* YPeriph(NULL);
        TH2F* YDifference(NULL);
        TH2F* YTklCentral(NULL);
        TH2F* YTklCentralME(NULL);
        TH2F* YTklPeriph(NULL);
        TH2F* YTklPeriphME(NULL);
        TH2F* YTklDifference(NULL);
        TH1F* YTklCentral_proj_tampon(NULL);
        TH1F* YTklPeriph_proj_tampon(NULL);
        TH1F* YTklDifference_proj_tampon(NULL);
        
        YCentral = new TH2F("YCentral",
                          "Yield delta eta wrt delta phi - Central",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        YCentral->SetXTitle("Correlation DeltaPhi (rad)");
        YCentral->SetYTitle("Correlation DeltaEta");
        YPeriph = new TH2F("YPeriph",
                          "Yield delta eta wrt delta phi - Periph",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        YPeriph->SetXTitle("Correlation DeltaPhi (rad)");
        YPeriph->SetYTitle("Correlation DeltaEta");
        YDifference = new TH2F("YDifference",
                          "Yield delta eta wrt delta phi - Difference",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        YDifference->SetXTitle("Correlation DeltaPhi (rad)");
        YDifference->SetYTitle("Correlation DeltaEta");
        YTklCentral = new TH2F("YTklCentral",
                          "Yield delta eta wrt delta phi - Tkl Central",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        YTklCentral->SetXTitle("Correlation DeltaPhi (rad)");
        YTklCentral->SetYTitle("Correlation DeltaEta");
        YTklPeriph = new TH2F("YTklPeriph",
                          "Yield delta eta wrt delta phi - Tkl Periph",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        YTklPeriph->SetXTitle("Correlation DeltaPhi (rad)");
        YTklPeriph->SetYTitle("Correlation DeltaEta");
        YTklDifference = new TH2F("YTklDifference",
                          "Yield delta eta wrt delta phi - Tkl Difference",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        YTklDifference->SetXTitle("Correlation DeltaPhi (rad)");
        YTklDifference->SetYTitle("Correlation DeltaEta");
        YTklCentralME = new TH2F("YTklCentralME",
                          "Yield delta eta wrt delta phi - Tkl Central",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        YTklCentralME->SetXTitle("Correlation DeltaPhi (rad)");
        YTklCentralME->SetYTitle("Correlation DeltaEta");
        YTklPeriphME = new TH2F("YTklPeriphME",
                          "Yield delta eta wrt delta phi - Tkl Periph",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        YTklPeriphME->SetXTitle("Correlation DeltaPhi (rad)");
        YTklPeriphME->SetYTitle("Correlation DeltaEta");
        // AE
        
        TH2F* hnseg9(NULL);
        TH2F* hnseg10(NULL);
        TH1F* hDPhi(NULL);
        TH2F* hPtWrtMassInv[3]{NULL};
        TH1 *hPtWrtMassInvSliced[3][6]= {NULL};
        
        TH1F* CentV0M(NULL);
        TH1F* CentTKL(NULL);
        TH1F* CentCL0(NULL);
        TH1F* CentCL1(NULL);
        TH1F* CentSPDTracklets(NULL);
        TH1F* CentSPDClusters(NULL);
        TH2F* CentV0Mwrttkl(NULL);
        TH2F* CentTKLwrttkl(NULL);
        TH2F* CentCL0wrttkl(NULL);
        TH2F* CentCL1wrttkl(NULL);
        TH2F* CentSPDTrackletswrttkl(NULL);
        TH2F* CentSPDClusterswrttkl(NULL);
        
        TH1F* Yields_Central_1(NULL);
        TH1F* Yields_Periph_1(NULL);
        TH1F* Yields_Difference_1(NULL);
        
        TH1F* coefficients0(NULL);
        TH1F* coefficients1(NULL);
        TH1F* coefficients2(NULL);
        TH1F* baselines0(NULL);
        TH1F* c2b0(NULL);
      //  TH1F* V2JPsiTkl(NULL);
        
        ProjCopy = new TH1F("ProjCopy",
                          "ProjCopy",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        ProjCopy->SetXTitle("Correlation DeltaPhi (rad)");
        ProjCopy->SetYTitle("Correlation DeltaEta");
        ProjCopy2 = new TH1F("ProjCopy2",
                          "ProjCopy2",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        ProjCopy2->SetXTitle("Correlation DeltaPhi (rad)");
        ProjCopy2->SetYTitle("Correlation DeltaEta");
        
        ProjCopyTkl = new TH1F("ProjCopyTkl",
                          "ProjCopyTkl",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        ProjCopyTkl->SetXTitle("Correlation DeltaPhi (rad)");
        ProjCopyTkl->SetYTitle("Correlation DeltaEta");
        ProjCopy2Tkl = new TH1F("ProjCopy2Tkl",
                          "ProjCopy2Tkl",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        ProjCopy2Tkl->SetXTitle("Correlation DeltaPhi (rad)");
        ProjCopy2Tkl->SetYTitle("Correlation DeltaEta");
        
        
        Yields_Central_1 = new TH1F("Yields_Central_1",
                         "Yield of JPsi-tkl in Central collisions wrt Phi",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Central_1->SetXTitle("Phi of the correlation (rad)");
        Yields_Central_1->SetYTitle("Yield");
        
        Yields_Periph_1 = new TH1F("Yields_Periph_1",
                         "Yield of JPsi-tkl in Periph collisions wrt Phi",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Periph_1->SetXTitle("Phi of the correlation (rad)");
        Yields_Periph_1->SetYTitle("Yield");
        
        Yields_Difference_1 = new TH1F("Yields_Difference_1",
                         "Yield of JPsi-tkl in (Central-Periph) collisions wrt Phi",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Difference_1->SetXTitle("Phi of the correlation (rad)");
        Yields_Difference_1->SetYTitle("Yield");
        
        coefficients0 = new TH1F("coefficients0",
                            "Coefficient Fourrier 0",
                            NbinsInvMass,MinInvMass,MaxInvMass);
           coefficients0->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients0->SetYTitle("Coefficient 0 (Fourier)");
        coefficients1 = new TH1F("coefficients1",
                               "Coefficient Fourrier 1",
                               NbinsInvMass,MinInvMass,MaxInvMass);
              coefficients1->SetXTitle("Mass of dimuon (GeV/c^{2})");
              coefficients1->SetYTitle("Coefficient 1 (Fourier)");
        coefficients2 = new TH1F("coefficients2",
                               "Coefficient Fourrier 2",
                               NbinsInvMass,MinInvMass,MaxInvMass);
              coefficients2->SetXTitle("Mass of dimuon (GeV/c^{2})");
              coefficients2->SetYTitle("Coefficient 2 (Fourier)");
        baselines0 = new TH1F("baselines0",
                         "Baseline 0",
                         NbinsInvMass,MinInvMass,MaxInvMass);
        baselines0->SetXTitle("Mass of dimuon (GeV/c^{2})");
        baselines0->SetYTitle("Baseline 0 (YieldPeriph deltaPhi=0)");
        c2b0 = new TH1F("c2b0",
                         "Coefficient 2 + Baseline 0",
                         NbinsInvMass,MinInvMass,MaxInvMass);
        c2b0->SetXTitle("Mass of dimuon (GeV/c^{2})");
        c2b0->SetYTitle("Coefficient 2 + Baseline 0");
        V2JPsiTkl = new TH1F("V2JPsiTkl",
                               "V2JPsiTkl wrt Mass",
                               NbinsInvMass,MinInvMass,MaxInvMass);
              V2JPsiTkl->SetXTitle("Mass of dimuon (GeV/c^{2})");
              V2JPsiTkl->SetYTitle("V2JPsiTkl");
        
        char hname[100];
        char hname2[100];
        
        for (int j=0; j <NbinsInvMass; j++){
           sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
           Yield_Central_MassBin[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           }
           
           for (int j=0; j <NbinsInvMass; j++){
           sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
           Yield_Periph_MassBin[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           }
           
           for (int j=0; j <NbinsInvMass; j++){
           sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
           Yield_Difference_MassBin[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           }
        
        TH1F* InvMass_Central(NULL);
        TH1F* InvMass_Periph(NULL);
        
    //    TH1F* YieldWrtMass_PhiBin_Central[NbinsDeltaPhi] = { NULL };
    //    TH1F* YieldWrtMass_PhiBin_Periph[NbinsDeltaPhi] = { NULL };
    //
    //    for (int j=0; j<NbinsDeltaPhi; j++){
    //        sprintf(hname,"Yield wrt Mass, Phi: %d pi/6 to %d pi/6 - Central",j-3,j-2);
    //        sprintf(hname2,"YieldWrtMassPhiBin_Central_%d",j);
    //        YieldWrtMass_PhiBin_Central[j] = new TH1F(hname2, hname,NbinsInvMass,MinInvMass,MaxInvMass);
    //    }
    //
    //    for (int j=0; j<NbinsDeltaPhi; j++){
    //        sprintf(hname,"Yield wrt Mass, Phi: %d pi/6 to %d pi/6 - Periph",j-3,j-2);
    //        sprintf(hname2,"YieldWrtMassPhiBin_Periph_%d",j);
    //        YieldWrtMass_PhiBin_Periph[j] = new TH1F(hname2, hname,NbinsInvMass,MinInvMass,MaxInvMass);
    //    }
        
        
        for(int i=0; i<NbBinsCent; i++){
            for(int k=0; k<NbinsInvMass; k++){
                char hname[100];
                sprintf(hname,"Yields %d %d ",i,k);
                Yields[i][k] = new TH1F(hname,
                                  "Yields Correlation wrt delta phi",
                                  NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
                Yields[i][k]->SetXTitle("Correlation DeltaPhi (rad)");
                for(int j=0; j<NbBinsZvtx; j++){
                    sprintf(hname,"Correlations %d %d %d ",i,j,k);
                    Correlations[i][j][k] = new TH2F(hname,
                                      "Correlation delta eta wrt delta phi",
                                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                    Correlations[i][j][k]->SetXTitle("Correlation DeltaPhi (rad)");
                    Correlations[i][j][k]->SetYTitle("Correlation DeltaEta");
                    sprintf(hname,"CorrelationsME %d %d %d ",i,j,k);
                    CorrelationsME[i][j][k] = new TH2F(hname,
                                      "MixedEvent Correlation delta eta wrt delta phi",
                                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                    CorrelationsME[i][j][k]->SetXTitle("Correlation DeltaPhi (rad)");
                    CorrelationsME[i][j][k]->SetYTitle("Correlation DeltaEta");
                }
            }
            
        }
        
        for(int i=0; i<NbBinsCent; i++){
            char hname[100];
            sprintf(hname,"Yields Tkl %d ",i);
            YieldsTkl[i] = new TH1F(hname,
                              "Yields Correlation Tkl wrt delta phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
            YieldsTkl[i]->SetXTitle("Correlation DeltaPhi (rad)");
            for(int j=0; j<NbBinsZvtx; j++){
                sprintf(hname,"Correlations Tkl %d %d ",i,j);
                CorrelationsTkl[i][j] = new TH2F(hname,
                                  "Correlation Tkl delta eta wrt delta phi",
                                  NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
                CorrelationsTkl[i][j]->SetXTitle("Correlation DeltaPhi (rad)");
                CorrelationsTkl[i][j]->SetYTitle("Correlation DeltaEta");
                sprintf(hname,"CorrelationsME Tkl %d %d ",i,j);
                CorrelationsTklME[i][j] = new TH2F(hname,
                                  "MixedEvent Correlation Tkl delta eta wrt delta phi",
                                  NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
                CorrelationsTklME[i][j]->SetXTitle("Correlation DeltaPhi (rad)");
                CorrelationsTklME[i][j]->SetYTitle("Correlation DeltaEta");
            }
        }
        
        for(int i=0; i<NbBinsCent; i++){
               for(int p=0; p<NbinsDeltaPhi; p++){
                   char hname[100];
                   sprintf(hname,"Yields PhiBin %d %d ",i,p);
                   Yields_PhiBin[i][p] = new TH1F(hname,
                                     "Yields Correlation wrt mass",
                                     NbinsInvMass,MinInvMass,MaxInvMass);
                   Yields_PhiBin[i][p]->SetXTitle("Correlation Inv Mass (GeV)");
           }
        }
            
        for(int p=0; p<NbinsDeltaPhi; p++){
                char hname[100];
                char hname2[100];
                sprintf(hname,"Yields PhiBin %d All C",p);
                YieldWrtMass_allC[p] = new TH1F(hname,
                                  "Yields Correlation wrt mass, all C",
                                  NbinsInvMass,MinInvMass,MaxInvMass);
                YieldWrtMass_allC[p]->SetXTitle("Correlation Inv Mass (GeV)");
            
            sprintf(hname,"Yields PhiBin %d Periph ",p);
            sprintf(hname2,"Yields Correlation wrt mass, Periph, Phi: %d pi/6 to %d pi/6",p-3, p-2);
            YieldWrtMass_Periph[p] = new TH1F(hname,
                              hname2,
                              NbinsInvMass,MinInvMass,MaxInvMass);
            YieldWrtMass_Periph[p]->SetXTitle("Correlation Inv Mass (GeV)");
            
            sprintf(hname,"Yields PhiBin %d Central",p);
            sprintf(hname2,"Yields Correlation wrt mass, Central, Phi: %d pi/6 to %d pi/6",p-3, p-2);
            YieldWrtMass_Central[p] = new TH1F(hname,
                              hname2,
                              NbinsInvMass,MinInvMass,MaxInvMass);
            YieldWrtMass_Central[p]->SetXTitle("Correlation Inv Mass (GeV)");
        }
            
        
        Yield_allC = new TH1F("Yield_allC",
                             "Yields Correlation wrt delta phi, all C, all mass",
                             NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Yield_allC->SetXTitle("Correlation DeltaPhi (rad)");
        Yield_Central = new TH1F("Yield_Central",
                          "Yields Correlation wrt delta phi, Central, 3.0-3.25 GeV",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yield_Central->SetXTitle("Correlation DeltaPhi (rad)");
        Yield_Periph = new TH1F("Yield_Periph",
                          "Yields Correlation wrt delta phi, Periph, 3.0-3.25 GeV",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yield_Periph->SetXTitle("Correlation DeltaPhi (rad)");
        Yield_Difference = new TH1F("Yield_Difference",
                          "Yields Correlation wrt delta phi, Central-Periph, all mass",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yield_Difference->SetXTitle("Correlation DeltaPhi (rad)");
           Yield_tampon = new TH1F("Yield_tampon",
                             "Yields Correlation wrt delta phi, all C, all mass",
                             NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Yield_tampon->SetXTitle("Correlation DeltaPhi (rad)");
        YieldWrtMass_tampon = new TH1F("YieldWrtMass_tampon",
                          "Yields Correlation wrt inv mass",
                          NbinsInvMass,MinInvMass,MaxInvMass);
        YieldWrtMass_tampon->SetXTitle("Correlation Inv Mass (GeV)");
        
        YieldTkl_allC = new TH1F("YieldTkl_allC",
                             "Yields Correlation Tkl wrt delta phi, all C, all mass",
                             NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
           YieldTkl_allC->SetXTitle("Correlation DeltaPhi (rad)");
        YieldTkl_Central = new TH1F("YieldTkl_Central",
                          "Yields Correlation Tkl wrt delta phi, Central",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldTkl_Central->SetXTitle("Correlation DeltaPhi (rad)");
        YieldTkl_Periph = new TH1F("YieldTkl_Periph",
                          "Yields Correlation Tkl wrt delta phi, Periph",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldTkl_Periph->SetXTitle("Correlation DeltaPhi (rad)");
        YieldTkl_Difference = new TH1F("YieldTkl_Difference",
                          "Yields Correlation Tkl wrt delta phi, Central-Periph, all mass",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldTkl_Difference->SetXTitle("Correlation DeltaPhi (rad)");
        YieldTkl_tampon = new TH1F("YieldTkl_tampon",
                          "Yields Correlation Tkl wrt delta phi, all C, all mass",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldTkl_tampon->SetXTitle("Correlation DeltaPhi (rad)");
        
        
        hnseg = new TH1F("hnseg",
                         "Invariant mass of dimuon",
                         250,MinInvMass,MaxInvMass);
        hnseg->SetXTitle("Mass of dimuon (GeV/c^{2})");
        hnseg->SetYTitle("Count");
        hDPhi = new TH1F("hDPhi",
                         "DeltaPhi of Tracklet",
                         10000,-TMath::Pi(),TMath::Pi());
        hDPhi->SetXTitle("DPhi tracklet");
        hDPhi->SetYTitle("Count");
        InvMass_Central = new TH1F("InvMass_Central",
                         "Invariant mass of dimuon - Central",
                         250,MinInvMass,MaxInvMass);
        InvMass_Central->SetXTitle("Mass of dimuon (GeV/c^{2})");
        InvMass_Central->SetYTitle("Count");
        InvMass_Periph = new TH1F("InvMass_Periph",
                         "Invariant mass of dimuon - Periph",
                         250,MinInvMass,MaxInvMass);
        InvMass_Periph->SetXTitle("Mass of dimuon (GeV/c^{2})");
        InvMass_Periph->SetYTitle("Count");
        hnseg2 = new TH1F("hnseg2",
                         "Pt of dimuon",
                         100,0,10);
        hnseg2->SetXTitle("Dimuon p_{t} (GeV/c)");
        hnseg2->SetYTitle("Count");
        hnseg3 = new TH1F("hnseg3",
                         "Eta of dimuon",
                         1000,-10,0);
        hnseg3->SetXTitle("Dimuon Eta");
        hnseg3->SetYTitle("Count");
        hnseg4 = new TH1F("hnseg4",
                          "Y of dimuon",
                          1000,-2.4,-4.2);
        hnseg4->SetXTitle("Dimuon Y");
        hnseg4->SetYTitle("Count");
        hnseg5 = new TH1F("hnseg5",
                          "Phi of dimuon",
                          1000,0, 2*TMath::Pi());
        hnseg5->SetXTitle("Dimuon Phi (rad)");
        hnseg5->SetYTitle("Count");
        hnseg6 = new TH1F("hnseg6",
                          "Correlation delta phi",
                          1000,MinDeltaPhi,MaxDeltaPhi);
        hnseg6->SetXTitle("Correlation DeltaPhi (rad)");
        hnseg6->SetYTitle("Count");
        hnseg7 = new TH1F("hnseg7",
                          "Correlation delta eta",
                          1000,MinDeltaEta,MaxDeltaEta);
        hnseg7->SetXTitle("Correlation DeltaEta");
        hnseg7->SetYTitle("Count");
        hnsegSigma = new TH2F("hnsegSigma",
                          "Zvtx resolution wrt NContributors",
                          60,0,60,1000,0,10);
        hnsegSigma->SetXTitle("NContributors");
        hnsegSigma->SetYTitle("Zvtx Resolution (SPD)");
        hnseg8 = new TH2F("hnseg8",
                          "Correlation delta eta wrt delta phi",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        hnseg8->SetXTitle("Correlation DeltaPhi (rad)");
        hnseg8->SetYTitle("Correlation DeltaEta");
        hnseg8ME = new TH2F("hnseg8ME",
                          "MixedEvent Correlation delta eta wrt delta phi",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        hnseg8ME->SetXTitle("Correlation DeltaPhi (rad)");
        hnseg8ME->SetYTitle("Correlation DeltaEta");
        hnseg8_proj = new TH1F("hnseg8_proj", "Correlation wrt delta eta, projection",NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        hnseg8_proj->SetXTitle("Correlation DeltaEta");
        hnseg8_proj->SetYTitle("Correlation mean");
        hnseg8Tkl_proj = new TH1F("hnseg8Tkl_proj", "Correlation Tkl wrt delta eta, projection",NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        hnseg8Tkl_proj->SetXTitle("Correlation DeltaEta");
        hnseg8Tkl_proj->SetYTitle("Correlation mean");
        hnseg8ME_proj = new TH1D("hnseg8ME_proj", "MixedEvent Correlation wrt delta eta, projection",NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        hnseg8ME_proj->SetXTitle("Correlation DeltaEta");
        hnseg8ME_proj->SetYTitle("Correlation mean");
        hnseg8TklME_proj = new TH1D("hnseg8TklME_proj", "MixedEvent Correlation Tkl wrt delta eta, projection",NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        hnseg8TklME_proj->SetXTitle("Correlation DeltaEta");
        hnseg8TklME_proj->SetYTitle("Correlation mean");
     
        Sijk = new TH1F("Sijk",
                                       "Sijk",
                                       NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Sijk->SetXTitle("Correlation DeltaPhi (rad)");
        Mijk = new TH1F("Mijk",
                                       "Mijk",
                                       NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Mijk->SetXTitle("Correlation DeltaPhi (rad)");
        SoverMijk = new TH1F("SoverMijk",
                                       "SoverMijk",
                                       NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        SoverMijk->SetXTitle("Correlation DeltaPhi (rad)");
        SoverMik = new TH1F("SoverMik",
                                       "SoverMik",
                                       NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        SoverMik->SetXTitle("Correlation DeltaPhi (rad)");
        
        Sij = new TH1F("Sij",
                                       "Sij",
                                       NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        Sij->SetXTitle("Correlation DeltaPhi (rad)");
        Mij = new TH1F("Mij",
                                       "Mij",
                                       NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        Mij->SetXTitle("Correlation DeltaPhi (rad)");
        SoverMij = new TH1F("SoverMij",
                                       "SoverMij",
                                       NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        SoverMij->SetXTitle("Correlation DeltaPhi (rad)");
        SoverMi = new TH1F("SoverMi",
                                       "SoverMi",
                                       NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        SoverMi->SetXTitle("Correlation DeltaPhi (rad)");

        hnseg9 = new TH2F("hnseg9",
                          "Tracklet eta wrt z_vertex",
                          100,-15,15,100,-2.5,2.5);
        hnseg9->SetXTitle("Z_{vertex} (cm)");
        hnseg9->SetYTitle("Tracklet Eta");
        hnseg10 = new TH2F("hnseg10",
                          "Number of tracklets wrt V0M percentile",
                          100,0,100,100,0,400);
        hnseg10->SetXTitle("V0M Percentile");
        hnseg10->SetYTitle("Tracklet count");
        Double_t binsPt[7] = {0,2,3,4,6,8,12};
        hPtWrtMassInv[0] = new TH2F("hPtWrtMassInv_0",
                          "Pt wrt Mass Inv - All centralities",
                          250,MinInvMass,MaxInvMass,6,binsPt);
        hPtWrtMassInv[0]->SetXTitle("Mass Inv");
        hPtWrtMassInv[0]->SetYTitle("Pt");
        hPtWrtMassInv[1] = new TH2F("hPtWrtMassInv_1",
                          "Pt wrt Mass Inv - Central",
                          250,MinInvMass,MaxInvMass,6,binsPt);
        hPtWrtMassInv[1]->SetXTitle("Mass Inv");
        hPtWrtMassInv[1]->SetYTitle("Pt");
        hPtWrtMassInv[2] = new TH2F("hPtWrtMassInv_2",
                          "Pt wrt Mass Inv - Periph",
                          250,MinInvMass,MaxInvMass,6,binsPt);
        hPtWrtMassInv[2]->SetXTitle("Mass Inv");
        hPtWrtMassInv[2]->SetYTitle("Pt");
        
        
        CentV0M = new TH1F("CentV0M",
                         "Distribution of CentV0M",
                         100,0,100);
        CentV0M->SetXTitle("CentV0M");
        CentV0M->SetYTitle("Count");
        CentTKL = new TH1F("CentTKL",
                         "Distribution of CentTKL",
                         100,0,100);
        CentTKL->SetXTitle("CentTKL");
        CentTKL->SetYTitle("Count");
        CentCL0 = new TH1F("CentCL0",
                         "Distribution of CentCL0",
                         100,0,100);
        CentCL0->SetXTitle("CentCL0");
        CentCL0->SetYTitle("Count");
        CentCL1 = new TH1F("CentCL1",
                         "Distribution of CentCL1",
                         100,0,100);
        CentCL1->SetXTitle("CentCL1");
        CentCL1->SetYTitle("Count");
        CentSPDTracklets = new TH1F("CentSPDTracklets",
                         "Distribution of CentSPDTracklets",
                         100,0,100);
        CentSPDTracklets->SetXTitle("CentSPDTracklets");
        CentSPDTracklets->SetYTitle("Count");
        CentSPDClusters = new TH1F("CentSPDClusters",
                         "Distribution of CentSPDClusters",
                         100,0,100);
        CentSPDClusters->SetXTitle("CentSPDClusters");
        CentSPDClusters->SetYTitle("Count");
        
        CentV0Mwrttkl = new TH2F("CentV0Mwrttkl",
                         "Distribution of CentV0Mwrttkl",
                         100,0,100, 100, 0, 100);
        CentV0Mwrttkl->SetXTitle("CentV0Mwrttkl");
        CentV0Mwrttkl->SetYTitle("NTracklets");
        CentTKLwrttkl = new TH2F("CentTKLwrttkl",
                         "Distribution of CentTKLwrttkl",
                         100,0,100, 100, 0, 100);
        CentTKLwrttkl->SetXTitle("CentTKLwrttkl");
        CentTKLwrttkl->SetYTitle("NTracklets");
        CentCL0wrttkl = new TH2F("CentCL0wrttkl",
                         "Distribution of CentCL0wrttkl",
                         100,0,100, 100, 0, 100);
        CentCL0wrttkl->SetXTitle("CentCL0wrttkl");
        CentCL0wrttkl->SetYTitle("NTracklets");
        CentCL1wrttkl = new TH2F("CentCL1wrttkl",
                         "Distribution of CentCL1wrttkl",
                         100,0,100, 100, 0, 100);
        CentCL1wrttkl->SetXTitle("CentCL1wrttkl");
        CentCL1wrttkl->SetYTitle("NTracklets");
        CentSPDTrackletswrttkl = new TH2F("CentSPDTrackletswrttkl",
                         "Distribution of CentSPDTrackletswrttkl",
                         100,0,100, 100, 0, 100);
        CentSPDTrackletswrttkl->SetXTitle("CentSPDTrackletswrttkl");
        CentSPDTrackletswrttkl->SetYTitle("NTracklets");
        CentSPDClusterswrttkl = new TH2F("CentSPDClusterswrttkl",
                         "Distribution of CentSPDClusterswrttkl",
                         100,0,100, 100, 0, 100);
        CentSPDClusterswrttkl->SetXTitle("CentSPDClusterswrttkl");
        CentSPDClusterswrttkl->SetYTitle("NTracklets");
    // *************************
    // Analyse                 *
    // *************************
        
        
        int DimuCentralSeen = 0;
        int DimuPeriphSeen = 0;
        int DimuSeenMassCut = 0;
        int DimuSeenNoMassCut = 0;
        int EventRejected = 0;
        int EventNC = 0;
        int EventPileUpMult = 0;
        int EventPileUpVtx = 0;
        int RefTracklets = 0;
        int RefTrackletsCentral = 0;
        int RefTrackletsPeriph = 0;
        int cmul = 0;
        int barWidth = 50;
        std::vector <double> Pools[NbBinsCent][NbBinsZvtx]; //[12]
        int PoolsSize[NbBinsCent][NbBinsZvtx] = {0}; //[12]
        std::vector <double> PoolsTkl[NbBinsCent][NbBinsZvtx];
        int PoolsSizeTkl[NbBinsCent][NbBinsZvtx] = {0};
        int DimuonCounter[NbBinsCent][NbBinsZvtx][NbinsInvMass] = {0};
        int DimuonCounterZint[NbBinsCent][NbinsInvMass] = {0};
        int RefTklCounter[NbBinsCent][NbBinsZvtx] = {0};
        int RefTklCounterZint[NbBinsCent] = {0};
        double NormME[NbBinsCent][NbBinsZvtx][NbinsInvMass] = {0};
        double NormMETkl[NbBinsCent][NbBinsZvtx] = {0};
        
        int DimuC = 0;
        int DimuP = 0;
        int TklC = 0;
        int TklP = 0;
        int NormTklCentral = 0;
        int NormTklPeriph = 0;
        int countsigma = 0;
    
    
    TCanvas *cinvmass = new TCanvas("cinvmass","Fitting All Phi - Different Centralities",10,10,700,500);
    cinvmass->SetTitle("Inv Mass Fits");
    // ROOT THAT Inv mass fits
    cinvmass->Divide(1,3);
    
    // Récup histos
    
    TFile *filerec = new TFile(FitFileName);
    hnseg = (TH1F*)filerec->Get("hnseg");
    InvMass_Central = (TH1F*)filerec->Get("InvMass_Central");
    InvMass_Periph = (TH1F*)filerec->Get("InvMass_Periph");
    V2JPsiTkl = (TH1F*)filerec->Get("V2JPsiTkl");
    coefficients0 = (TH1F*)filerec->Get("coefficients0");
    coefficients1 = (TH1F*)filerec->Get("coefficients1");
    coefficients2 = (TH1F*)filerec->Get("coefficients2");
    baselines0 = (TH1F*)filerec->Get("baselines0");
    
    
    for(int p=0; p<NbinsDeltaPhi; p++){
        sprintf(hname,"Yields PhiBin %d Central",p);
        YieldWrtMass_Central[p] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"Yields PhiBin %d Periph",p);
        YieldWrtMass_Periph[p] = (TH1F*)filerec->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    Yield_Central_MassBin[j] = (TH1F*)filerec->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    Yield_Periph_MassBin[j] = (TH1F*)filerec->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    Yield_Difference_MassBin[j] = (TH1F*)filerec->Get(hname);
    }
    
    
    cout << "START ALL"<<endl;
    
    
    
    
    
    
    
    
    
    
    
    TFitResultPtr res;
    TFitResultPtr rescent;
    TFitResultPtr resperiph;
    TFitResultPtr resu;
    
   ROOT::Math::Minimizer* minim{NULL};
    
    TCanvas*c16=new TCanvas();
    // ROOT THAT Extraction 2 plot
    c16->SetTitle("Extraction method 2");
    
    char histoname[50];
    int niterations = 100;
    TRandom rand = 0;
    rand.SetSeed(1234);
    
    sprintf(histoname,"hnseg");
        res = FittingAllInvMassBin(histoname, cinvmass, 0);
    double par[16];
    double parerr[16];
    
    if(!Extraction2Combined){
        niterations = 1;
    }
    double par12[niterations];
    double parerr12[niterations];

    for(int trial = 0; trial<niterations; trial++){
                // Parameter setting (either from mass fit either start values)
                cout << "Parameter Setting"<<endl;
        if(!Extraction2Combined){
                for(int i=0; i <16; i++){
                        par[i] = res->Parameter(i);
                        if(i>11){
                            par[i] = 0.1;
                        }
            
                }
        }
        if(Extraction2Combined){
            for(int i=0; i <12; i++){
                        parerr[i] = res->ParError(i);
                        par[i] = res->Parameter(i);
                cout << "par " << i << " = " << par[i] << " before Gauss" << endl;
                        par[i] += rand.Gaus(0,parerr[i]);
                cout << "par " << i << " = " << par[i] << " after Gauss" << endl;
                }
            for(int i=12; i<16; i++){
                par[i] = 0.1;
            }
        }

                //Definition and design of function on V2 plot
                c16->cd();
                TF1 *fitV2_2 = new TF1("fitV2_2",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
                fitV2_2->SetNpx(500);
                fitV2_2->SetLineWidth(4);
                fitV2_2->SetLineColor(kMagenta);
                fitV2_2->SetParameters(par);
                TF1 *fitJustV2_2 = new TF1("fitJustV2_2",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,4);
                fitJustV2_2->SetNpx(500);
                fitJustV2_2->SetLineWidth(0);
                fitJustV2_2->SetLineColor(kWhite);
                TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2,2.1,5.1,16);
                backFcnV2_2->SetLineColor(kRed);
                TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,2.1,5.1,16);
                signalFcnJPsiV2_2->SetLineColor(kBlue);
                signalFcnJPsiV2_2->SetNpx(500);
        
        if(trial==niterations-1){
            V2JPsiTkl->GetListOfFunctions()->Clear();
            V2JPsiTkl->GetListOfFunctions()->Add(fitV2_2);
        }


                   TVirtualFitter::Fitter(V2JPsiTkl)->SetMaxIterations(10000);
                   TVirtualFitter::Fitter(V2JPsiTkl)->SetPrecision();
                    for(int i=0; i<=11; i++){
                        fitV2_2->FixParameter(i,par[i]);
                    }
                

                   fitV2_2->SetParName(0,"Norm_{JPsi}");
                   fitV2_2->SetParName(1,"M_{JPsi}");
                   fitV2_2->SetParName(2,"Sigma_{JPsi}");
                   fitV2_2->SetParName(3,"a_{1}");
                   fitV2_2->SetParName(4,"n_{1}");
                   fitV2_2->SetParName(5,"a_{2}");
                   fitV2_2->SetParName(6,"n_{2}");
                   fitV2_2->SetParName(7,"Norm_{Psi2S}");
                   fitV2_2->SetParName(8,"Norm_{TailLowM}");
                   fitV2_2->SetParName(9,"Exp_{TailLowM}");
                   fitV2_2->SetParName(10,"Norm_{TailHighM}");
                   fitV2_2->SetParName(11,"Exp_{TailHighM}");
                    fitV2_2->SetParName(12,"V2_2 JPsi");
                    fitV2_2->SetParName(13,"V2_2 Bkg M2");
                fitV2_2->SetParName(14,"V2_2 Bkg M1");
                fitV2_2->SetParName(15,"V2_2 Nkg M0");
//
                fitJustV2_2->SetParName(0,"V2_2 JPsi");
                    fitJustV2_2->SetParName(1,"V2_2 Bkg M2");
                fitJustV2_2->SetParName(2,"V2_2 Bkg M1");
                fitJustV2_2->SetParName(3,"V2_2 Nkg M0");


                //Fit of V2
                    resu = V2JPsiTkl->Fit("fitV2_2","SBMERI+","ep");
                    Double_t para[16];
                    fitV2_2->GetParameters(para);
                    backFcnV2_2->SetParameters(para);
                    signalFcnJPsiV2_2->SetParameters(para);
        
        if(resu->CovMatrixStatus() == 3){
        par12[trial] = para[12];
        parerr12[trial] = resu->ParError(12);
        }
        else{
            par12[trial]=parerr12[trial]=0;
        }

        if(trial==niterations-1){
                fitV2_2->Draw("same");
                signalFcnJPsiV2_2->Draw("same");
                backFcnV2_2->Draw("same");
                  // draw the legend
                  TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                  legend->SetTextFont(72);
                  legend->SetTextSize(0.04);
                Char_t message[80];
                sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
                 legend->AddEntry(fitV2_2,message);
                    if(resu->CovMatrixStatus() == 3){
                           sprintf(message,"The fit is a success");
                    }
                    else{
                           sprintf(message,"The fit is a failure");
                    }
                
                       legend->AddEntry(fitV2_2,message);
                legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
                legend->AddEntry(backFcnV2_2,"Background");
                  legend->AddEntry(V2JPsiTkl,"Data","lpe");
                  legend->Draw();
        }
        
    }
    
    if(Extraction2Combined){
        double moypar12 = 0;
        double moyparerr12 = 0;
        double typeApar12 = 0;
        double sigparerr12 = 0;
        int failures = 0;
        int NewOutliers = 1;
        double SigmaCutOutliers = 10;
        
        while(NewOutliers !=0){
            moypar12 = 0;
            moyparerr12 = 0;
            typeApar12 = 0;
            sigparerr12 = 0;
            int failures = 0;
            failures = 0;
            NewOutliers = 0;
            for(int z=0; z<niterations; z++){
                moypar12 += par12[z];
                moyparerr12 += parerr12[z];
                if(par12[z]==0 && parerr12[z]==0){
                    failures++;
                }
            }
            moypar12/=(niterations-failures);   // Moyenne de par12
            moyparerr12/=(niterations-failures); //Incertitude fit 2
            for(int z=0; z<niterations; z++){
                if(par12[z]!=0 || parerr12[z]!=0){
                    typeApar12 += pow(moypar12-par12[z],2);
                }
                if(par12[z]!=0 || parerr12[z]!=0){
                    sigparerr12 += pow(moyparerr12-parerr12[z],2);
                }
            }
            
            typeApar12/=((niterations-failures)*(niterations-1-failures));
            typeApar12 = sqrt(typeApar12); // Incertitude de type A (répétabitilé) issue du fit 1 propagé
            
            sigparerr12/=((niterations-failures)*(niterations-1-failures));
            sigparerr12 = sqrt(sigparerr12); // Ecart-type valeurs des erreurs.
            
            std::cout << "moypar12: " << moypar12 << ", moyparerr12: " << moyparerr12 << ", typeApar12: "<<typeApar12 << ", sigparerr12: "<<sigparerr12<<endl;

            for(int z=0; z<niterations; z++){
                if(par12[z]!=0 || parerr12[z]!=0){
                    if(TMath::Abs(par12[z]-moypar12)>(SigmaCutOutliers*typeApar12) && TMath::Abs(parerr12[z]-moyparerr12)>(SigmaCutOutliers*sigparerr12)){
                        NewOutliers++;
                        std::cout << "Iteration with Yield Value: " << par12[z] << " and Error: " << parerr12[z] << " is an outlier. Suppressed"<<endl;
                    
                        par12[z] =0;
                        parerr12[z] =0;
                    }
                }
            }

            std::cout << "There were " << failures << " failures"<<endl;
            std::cout << "moypar12: " << moypar12 << ", moyparerr12: " << moyparerr12 << ", typeApar12: "<<typeApar12<<endl;
            
        }
        
        cout << "Value of V2 is " << moypar12 <<endl;
        cout << "Value of V2 error is " << sqrt(pow(moyparerr12,2)+pow(typeApar12,2)) <<endl;
    }
    
}

    
    
    
    
    
    
//FITTING METHODS

Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par)
// Par 0->7: signal, 8->11 Bkg, 12: v2 JPsi, 13->15: V2 bkg
{ return (JPsiCrystalBallExtended(x,par)*par[12] + (ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*x[0]*x[0] + par[14]*x[0] + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t BackFcnV2(Double_t *x,Double_t *par)
{return ((ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*x[0]*x[0] + par[14]*x[0] + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t SignalFcnJPsiV2(Double_t *x,Double_t *par)
{return (JPsiCrystalBallExtended(x,par)*par[12])/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t FourierV2(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);}

Double_t FourierV5(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);} //+ 2*par[3]*cos(3*x[0]) + 2*par[4]*cos(4*x[0]);} + 2*par[5]*cos(5*x[0]) + 2*par[6]*cos(6*x[0]) + 2*par[7]*cos(7*x[0]) + 2*par[8]*cos(8*x[0]) + 2*par[9]*cos(9*x[0]) + 2*par[10]*cos(10*x[0]) + 2*par[11]*cos(11*x[0]);}

Double_t TwoCBE2E(Double_t *x,Double_t *par)

{ return JPsiCrystalBallExtended(x,par) + Psi2SCrystalBallExtended(x,par) + ExpBkg(x,&par[8]);}

Double_t ExpBkg(Double_t *x,Double_t *par)

{ return par[0]*(exp(x[0]*par[1]*(-1))) + par[2]*(exp(x[0]*par[3]*(-1))); }

Double_t JPsiCrystalBallExtended(Double_t *x,Double_t *par)
{
    
    Double_t sum = 0;
    
    Double_t t = (x[0]-par[1])/par[2];
    if (par[3] < 0) t = -t;
    
    Double_t absAlpha = fabs((Double_t)par[3]);
    Double_t absAlpha2 = fabs((Double_t)par[5]);
    
    
    if (t < -absAlpha) //left tail
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        
        sum += (par[0]/(par[2]*sqrt(TMath::Pi()*2)))*(a/TMath::Power(b - t, par[4]));
    } else if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += (par[0]/(par[2]*sqrt(TMath::Pi()*2)))*(exp(-0.5*t*t));
    } else if (t >= absAlpha2) //right
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        
        sum += (par[0]/(par[2]*sqrt(TMath::Pi()*2)))*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
    return sum ;
}

Double_t Psi2SCrystalBallExtended(Double_t *x,Double_t *par)
{
      Double_t sum = 0;
    
    Double_t t = (x[0]-(par[1]*mPsip/mJpsi))/(par[2]*ratSigma);
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t < -absAlpha) //left tail
  {
      Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
      Double_t b = par[4]/absAlpha - absAlpha;
      
      sum += (par[7]/(par[2]*ratSigma*sqrt(2*TMath::Pi())))*(a/TMath::Power(b - t, par[4]));
  }
  else if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
      sum += (par[7]/(par[2]*ratSigma*sqrt(2*TMath::Pi())))*(exp(-0.5*t*t));
  }
  else if (t >= absAlpha2) //right tail
  {
      Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
      Double_t d = par[6]/absAlpha2 - absAlpha2;
      
      sum += (par[7]/(par[2]*ratSigma*sqrt(2*TMath::Pi())))*(c/TMath::Power(d + t, par[6]));
  } else
      sum += 0;
    
    return sum ;

}
    
    TFitResultPtr FittingAllInvMass(const char *histoname, TCanvas *cinvmass){
        TFitResultPtr res = FittingAllInvMassBin(histoname, cinvmass, 0);
        return res;
    }

TFitResultPtr FittingAllInvMassBin(const char *histoname, TCanvas *cinvmass, int i){
 //Bevington Exercise by Peter Malzacher, modified by Rene Brun, TO BE MODIFIED FOR OUR USE
    // PARAMETERS
    /*
     ------ Signal CB2
     0: Normalisation Gauss
     1: Mean x
     2: Sig Gauss
     3: Alpha 1
     4: Power 1
     5: Alpha 2
     6: Power 2
     7: Normalisation Psi2s, JPsi
     ------ Background 2Exp
     8: Normalisation Before
     9: Exponent Before
     10: Norlmalisation After
     11: Exponent After
     */
    TFile *file0 = new TFile(FitFileName);
    
    cinvmass->cd(i+1);
   cinvmass->SetFillColor(33);
   cinvmass->SetFrameFillColor(41);
   cinvmass->SetGrid();
   TH1F *histo = (TH1F*)file0->Get(histoname);
   // create a TF1 with the range from 0 to 3 and 6 parameters
   TF1 *fitFcn = new TF1("fitFcn",TwoCBE2E,2.1,5.1,12);
   fitFcn->SetNpx(500);
   fitFcn->SetLineWidth(4);
   fitFcn->SetLineColor(kMagenta);
   // first try without starting values for the parameters
   // This defaults to 1 for each param.
   // this results in an ok fit for the polynomial function
   // however the non-linear part (lorenzian) does not
   // respond well.
    Double_t params[12] = {10,1,1,1,1,1,1,1,1,1,1,1};
   fitFcn->SetParameters(params);
    TVirtualFitter::Fitter(histo)->SetMaxIterations(100000);
    TVirtualFitter::Fitter(histo)->SetPrecision();
 //  histo->Fit("fitFcn","0");
   // First fit, fix MJPsi and Sigma
    
    fitFcn->SetParLimits(0,0.1,100000000);
    fitFcn->SetParLimits(1,3.05,3.15);
    fitFcn->SetParLimits(2,0.03,0.1);
    fitFcn->SetParLimits(3,0.7,1.1);
    fitFcn->SetParLimits(4,1,30);
    fitFcn->SetParLimits(5,1,5);
    fitFcn->SetParLimits(6,5,25);
    fitFcn->SetParLimits(7,0.001,1000000);
    fitFcn->SetParLimits(8,0.1,1000000);
    fitFcn->SetParLimits(9,0.01,20);
    fitFcn->SetParLimits(10,0.01,1000000);
    fitFcn->SetParLimits(11,0.01,50);
    
   fitFcn->FixParameter(1,3.096916); // Mean x core
    fitFcn->FixParameter(2,0.07);
    fitFcn->FixParameter(3,0.9);
    fitFcn->FixParameter(4,10);
    fitFcn->FixParameter(5,2);
    fitFcn->FixParameter(6,15);
    fitFcn->SetParameter(7,100);
    fitFcn->SetParameter(8,1000);
    fitFcn->SetParameter(9,0.5);
    fitFcn->SetParameter(11,10);

    
    fitFcn->SetParName(0,"Norm_{JPsi}");
    fitFcn->SetParName(1,"M_{JPsi}");
    fitFcn->SetParName(2,"Sigma_{JPsi}");
    fitFcn->SetParName(3,"a_{1}");
    fitFcn->SetParName(4,"n_{1}");
    fitFcn->SetParName(5,"a_{2}");
    fitFcn->SetParName(6,"n_{2}");
    fitFcn->SetParName(7,"Norm_{Psi2S}");
    fitFcn->SetParName(8,"Norm_{TailLowM}");
    fitFcn->SetParName(9,"Exp_{TailLowM}");
    fitFcn->SetParName(10,"Norm_{TailHighM}");
    fitFcn->SetParName(11,"Exp_{TailHighM}");
    
   TFitResultPtr res = histo->Fit("fitFcn","SBMER","ep");
   res = histo->Fit("fitFcn","SBMERQ","ep");
      fitFcn->ReleaseParameter(1);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    fitFcn->ReleaseParameter(2);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    fitFcn->ReleaseParameter(7);
   res = histo->Fit("fitFcn","SBMERQ","ep");
   res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(3);
//    fitFcn->ReleaseParameter(4);
//   res = histo->Fit("fitFcn","SBMER","ep");
//   res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(5);
//    fitFcn->ReleaseParameter(6);
//   res = histo->Fit("fitFcn","SBMER","ep");
//    res = histo->Fit("fitFcn","SBMER","ep");

    
   // improve the picture:
   TF1 *backFcn = new TF1("backFcn",ExpBkg,2.1,5.1,4);
   backFcn->SetLineColor(kRed);
   TF1 *signalFcnJPsi = new TF1("signalFcnJPsi",JPsiCrystalBallExtended,2.1,5.1,8);
   TF1 *signalFcnPsi2S = new TF1("signalFcnPsi2S",Psi2SCrystalBallExtended,2.1,5.1,8);
   TPaveText *pave = new TPaveText(0.15,0.5,0.3,0.65,"brNDC");
   signalFcnJPsi->SetLineColor(kBlue);
   signalFcnJPsi->SetNpx(500);
    signalFcnPsi2S->SetLineColor(kGreen);
    signalFcnPsi2S->SetNpx(500);
   Double_t par[12];
   // writes the fit results into the par array
    gStyle->SetOptFit(1011);
   fitFcn->GetParameters(par);
   signalFcnJPsi->SetParameters(par);
    signalFcnPsi2S->SetParameters(par);
//   auto covMatrix = res->GetCovarianceMatrix();
//   std::cout << "Covariance matrix from the fit ";
//   covMatrix.Print();
   Double_t integral = (signalFcnJPsi->Integral(2.1,3.45))/0.01;
    auto covtot = res->GetCovarianceMatrix();
    auto covsgn = covtot.GetSub(0,8,0,8);
    std::cout << "Matrice totale" <<endl;
    covtot.Print();
    std::cout << "Matrice réduite" <<endl;
    covsgn.Print();
    std::cout << "STATUS COV " << res->CovMatrixStatus() <<endl;
   Double_t integralerror = (signalFcnJPsi->IntegralError(2.1,3.45,signalFcnJPsi->GetParameters(), res->GetCovarianceMatrix().GetSub(0,8,0,8).GetMatrixArray() ))/0.01;
    std::cout << "Erreur integrale " << integralerror <<endl;
    
    std::cout << "Fitted " << histoname << std::endl;
    std::cout << "Nb JPsi total: " << integral << std::endl;
 //   std::cout << "integral error: " << integralerror << std::endl;
   signalFcnJPsi->Draw("same");
    signalFcnPsi2S->Draw("same");
   backFcn->SetParameters(&par[8]);
   backFcn->Draw("same");
   // draw the legend
    char str[50];
    sprintf(str, "N_{JPsi} %f +/- %f", integral, integralerror);
   pave->AddText(str);
    sprintf(str, "M_{Psi2S} = %f, Sig_{Psi2S} = %f", mPsip, ratSigma*par[2]);
    pave->AddText(str);
   pave->Draw();
   TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
   legend->SetTextFont(72);
   legend->SetTextSize(0.03);
    Char_t message[80];
    sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcn->GetChisquare(),fitFcn->GetNDF());
     legend->AddEntry(fitFcn,message);
    if(res->CovMatrixStatus() == 3){
               sprintf(message,"The fit is a success");
           }
           else{
               sprintf(message,"The fit is a failure");
           }
           legend->AddEntry(fitFcn,message);
   legend->AddEntry(histo,"Data","lpe");
   legend->AddEntry(backFcn,"Background fit","l");
   legend->AddEntry(signalFcnJPsi,"JPsi Signal fit","l");
    legend->AddEntry(signalFcnPsi2S,"Psi2S Signal fit","l");
   legend->Draw();
    
//    cout << "Histogramme: " << histoname << endl;
//    for(int j=0; j<10; j++){
//        cout << "Bin de masse numero " << j << " - Signal : " << (signalFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 << " et Background : " << (backFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 <<endl;
//    }
    
    return res;
}

void FcnCombinedAllMass(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  TAxis *xaxis1  = hnseg->GetXaxis();
  TAxis *xaxis2  = V2JPsiTkl->GetXaxis();

  int nbinX1 = hnseg->GetNbinsX();
  int nbinX2 = V2JPsiTkl->GetNbinsX();

  double chi2 = 0;
  double x[1];
  double tmp;
  npfits = 0;
  for (int ix = 1; ix <= nbinX1; ++ix) {
    x[0] = xaxis1->GetBinCenter(ix);
      if ( hnseg->GetBinError(ix) > 0 ) {
        tmp = (hnseg->GetBinContent(ix) - TwoCBE2E(x,p))/hnseg->GetBinError(ix);
        chi2 += tmp*tmp;
        npfits++;
      }
  }
  for (int ix = 1; ix <= nbinX2; ++ix) {
     x[0] = xaxis2->GetBinCenter(ix);
      if ( V2JPsiTkl->GetBinError(ix) > 0 ) {
        tmp = (V2JPsiTkl->GetBinContent(ix) - FourierV2_WrtInvMass(x,p))/V2JPsiTkl->GetBinError(ix);
        chi2 += tmp*tmp;
        npfits++;
      }
  }
  fval = chi2;
}

int GetCent(double cent){
    if(cent <= 1){
        return 0;
    }
    else if(cent <= 10){
        return 1;
    }
    else if(cent <= 20){
        return 2;
    }
    else if(cent <= 40){
        return 3;
    }
    else if(cent <= 100){
        return 4;
    }
//    else if(cent <= 30){
//        return 6;
//    }
//    else if(cent <= 40){
//        return 7;
//    }
//    else if(cent <= 50){
//        return 8;
//    }
//    else if(cent <= 60){
//        return 9;
//    }
//    else if(cent <= 70){
//        return 10;
//    }
//    else if(cent <= 80){
//        return 11;
//    }
//    else if(cent <= 90){
//        return 12;
//    }
//    else{
//        return 13;
//    }
}


