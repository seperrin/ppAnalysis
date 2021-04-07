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

 void FitTrainingPtBinned();
 Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par);
 Double_t CvetanF(Double_t *x,Double_t *par);
void ChisquareCvetanF(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t CvetanFPtBinned(Double_t *x,Double_t *par);
void ChisquareCvetanFPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t CvetanFTKL(Double_t *x,Double_t *par);
Double_t ZYAM(Double_t *x,Double_t *par);
void ChisquareZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t ZYAMPtBinned(Double_t *x,Double_t *par);
void ChisquareZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate(Double_t *x,Double_t *par);
void ChisquarePRLTemplate(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplatePtBinned(Double_t *x,Double_t *par);
void ChisquarePRLTemplatePtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate_RidgeAndZero(Double_t *x,Double_t *par);
Double_t PRLTemplate_RidgeAndZeroPtBinned(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphAndG(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphAndGPtBinned(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphZYAM(Double_t *x,Double_t *par);
void ChisquarePRLTemplate_PeriphZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate_PeriphZYAMPtBinned(Double_t *x,Double_t *par);
void ChisquarePRLTemplate_PeriphZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
 Double_t BackFcnV2(Double_t *x,Double_t *par);
Double_t BackFcnV2Poly(Double_t *x,Double_t *par);
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

TH1F* Yields_Central_1(NULL);
TH1F* Yields_Periph_1(NULL);
TH1F* Baseline_Central_1(NULL);
TH1F* Baseline_Periph_1(NULL);
TH1F* Yields_Periph_1_MinusBaseline(NULL);
TH1F* Yields_Central_1_MinusBaseline(NULL);

TH1F* YieldTkl_Central(NULL);
TH1F* YieldTkl_Periph(NULL);
TH1F* BaselineTkl_Central(NULL);
TH1F* BaselineTkl_Periph(NULL);
TH1F* YieldTkl_Periph_MinusBaseline(NULL);
TH1F* YieldTkl_Central_MinusBaseline(NULL);

double baseline_periph = 1;
double errbaseline_periph = 1;
double baseline_central = 99999;
double errbaseline_central = 99999;

double baselineTKL_periph = 1;
double errbaselineTKL_periph = 1;
double baselineTKL_central = 9999;
double errbaselineTKL_central = 1;

double PtBins[] = {0,2,4,6,8,12};
const int NbPtBins = 5;
double LowDimuPtCut = 0;
double HighDimuPtCut = 12;

double V2_Ext1[NbPtBins] = {0};
double V2_Ext2[NbPtBins] = {0};
double V2_CvetanQuentin[NbPtBins] = {0};
double V2_CvetanQuentinMe[NbPtBins] = {0};
double V2_ZYAM[NbPtBins] = {0};
double V2_PRL[NbPtBins] = {0};
double V2_PRLPeriphZYAM[NbPtBins] = {0};
double errV2_Ext1[NbPtBins] = {0};
double errV2_Ext2[NbPtBins] = {0};
double errV2_CvetanQuentin[NbPtBins] = {0};
double errV2_CvetanQuentinMe[NbPtBins] = {0};
double errV2_ZYAM[NbPtBins] = {0};
double errV2_PRL[NbPtBins] = {0};
double errV2_PRLPeriphZYAM[NbPtBins] = {0};


double F_CvetanQuentin[NbPtBins] = {0};
double F_PRL[NbPtBins] = {0};
double F_PRLPeriphZYAM[NbPtBins] = {0};
double errF_CvetanQuentin[NbPtBins] = {0};
double errF_PRL[NbPtBins] = {0};
double errF_PRLPeriphZYAM[NbPtBins] = {0};


TH1F* Yields_Central_1PtBinned[NbPtBins] = {NULL};
TH1F* Yields_Periph_1PtBinned[NbPtBins] = {NULL};
TH1F* Baseline_Central_1PtBinned[NbPtBins] = {NULL};
TH1F* Baseline_Periph_1PtBinned[NbPtBins] = {NULL};
TH1F* Yields_Periph_1_MinusBaselinePtBinned[NbPtBins] = {NULL};
TH1F* Yields_Central_1_MinusBaselinePtBinned[NbPtBins] = {NULL};

double baseline_periphPtBinned[NbPtBins] = {NULL};
double errbaseline_periphPtBinned[NbPtBins] = {NULL};
double baseline_centralPtBinned[NbPtBins] = {NULL};
double errbaseline_centralPtBinned[NbPtBins] = {NULL};

Char_t FitFileName[200];

void FitTrainingPtBinned(){
    
    sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile1_Run2_0-5_40-90_pt0-2-4-6-8-12.root");
    
    TH1::SetDefaultSumw2();
        bool doTracklets = kFALSE;
        bool doMixedEvents = kFALSE;
    bool Extraction2Combined = kFALSE;
        bool CombineFits = kFALSE;
    bool PtBinned = kTRUE;
        
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
        int BinZeroLeftTKL = floor((0-MinDeltaPhi)*NbinsDeltaPhiTKL/(2*TMath::Pi()));
        
        const int NbBinsCent = 14;
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
    TH2F* CorrelationsMEScaled[NbBinsCent][NbBinsZvtx][NbinsInvMass]{ NULL };
    TH2F* CorrelationsTklMEScaled[NbBinsCent][NbBinsZvtx]{ NULL };
        TH1F* Yield_tampon(NULL);
        TH1F* YieldTkl_tampon(NULL);
        TH1F* YieldWrtMass_tampon(NULL);
        TH1F* Yield_allC(NULL);
        TH1F* Yield_Central(NULL);
        TH1F* Yield_Periph(NULL);
        TH1F* Yield_Difference(NULL);
        TH1F* YieldTkl_allC(NULL);
//        TH1F* YieldTkl_Central(NULL);
//        TH1F* YieldTkl_Periph(NULL);
        TH1F* YieldTkl_Difference(NULL);
        TH1F* Yield_Central_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yield_Periph_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yield_Difference_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yields_PhiBin[NbBinsCent][NbinsDeltaPhi]{ NULL };
        TH1F* YieldWrtMass_allC[NbinsDeltaPhi]{ NULL };
        TH1F* YieldWrtMass_Central[NbinsDeltaPhi]{ NULL };
        TH1F* YieldWrtMass_Periph[NbinsDeltaPhi]{ NULL };
        
    
    // Plots added for Pt Binning
    TH1F* YieldsPtBinned[NbBinsCent][NbinsInvMass][NbPtBins]{ NULL };//
    TH2I* CorrelationsMEPtBinned[NbBinsCent][NbBinsZvtx][NbinsInvMass][NbPtBins]{ NULL };//
    TH2F* CorrelationsPtBinned[NbBinsCent][NbBinsZvtx][NbinsInvMass][NbPtBins]{ NULL };//
    TH2F* CorrelationsMEScaledPtBinned[NbBinsCent][NbBinsZvtx][NbinsInvMass][NbPtBins]{ NULL };//
    TH1F* Yield_Central_MassBinPtBinned[NbinsInvMass][NbPtBins] = { NULL };//
    TH1F* Yield_Periph_MassBinPtBinned[NbinsInvMass][NbPtBins] = { NULL };//
    TH1F* Yield_Difference_MassBinPtBinned[NbinsInvMass][NbPtBins] = { NULL };//
    TH1F* Yields_PhiBinPtBinned[NbBinsCent][NbinsDeltaPhi][NbPtBins]{ NULL };//
    TH1F* YieldWrtMass_allCPtBinned[NbinsDeltaPhi][NbPtBins]{ NULL };//
    TH1F* YieldWrtMass_CentralPtBinned[NbinsDeltaPhi][NbPtBins]{ NULL };//
    TH1F* YieldWrtMass_PeriphPtBinned[NbinsDeltaPhi][NbPtBins]{ NULL };//
    TH1F* baselines0PtBinned[NbPtBins]{ NULL };
    TH1F* coefficients0PtBinned[NbPtBins]{ NULL };
    TH1F* coefficients1PtBinned[NbPtBins]{ NULL };
    TH1F* coefficients2PtBinned[NbPtBins]{ NULL };
    TH1F* c2b0PtBinned[NbPtBins]{ NULL };
    TH1F* V2JPsiTklPtBinned[NbPtBins]{ NULL };
    
//    TH1F* Yields_Central_1PtBinned[NbPtBins]{ NULL };
//       TH1F* Yields_Periph_1PtBinned[NbPtBins]{ NULL };
       TH1F* Yields_Difference_1PtBinned[NbPtBins]{ NULL };
    
    TH1F *Yields_Central_1_CvetanPtBinned[NbPtBins]{ NULL };
    TH1F *Yields_Central_1_CvetanMePtBinned[NbPtBins]{ NULL };
    TH1F *Yields_Central_1_ZYAMPtBinned[NbPtBins]{ NULL };
    TH1F *Yields_Central_1_PRLTemplatePtBinned[NbPtBins]{ NULL };
    TH1F *Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[NbPtBins]{ NULL };
    
    TH1F *YieldTkl_Central_CvetanPtBinned(NULL);
    TH1F *YieldTkl_Central_CvetanMePtBinned(NULL);
    TH1F *YieldTkl_Central_ZYAMPtBinned(NULL);
    TH1F *YieldTkl_Central_PRLTemplatePtBinned(NULL);
    TH1F *YieldTkl_Central_PRLTemplate_PeriphZYAMPtBinned(NULL);
    
    
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
        TH1 *hPtWrtMassInvSliced[3][NbPtBins]= {NULL};
        TH1 *hPtWrtMassInvSlicedRebinned[3][NbPtBins]= {NULL};
        
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
    
    Baseline_Central_1 = new TH1F("Baseline_Central_1",
                     "Baseline of JPsi-tkl in Central collisions wrt Phi",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Baseline_Central_1->SetXTitle("Phi of the correlation (rad)");
    Baseline_Central_1->SetYTitle("Baseline");
    
    Baseline_Periph_1 = new TH1F("Baseline_Periph_1",
                     "Baseline of JPsi-tkl in Periph collisions wrt Phi",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Baseline_Periph_1->SetXTitle("Phi of the correlation (rad)");
    Baseline_Periph_1->SetYTitle("Baseline");
    
    Yields_Central_1_MinusBaseline = new TH1F("Yields_Central_1_MinusBaseline",
                     "Yield of JPsi-tkl in Central collisions wrt Phi _MinusBaseline",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yields_Central_1_MinusBaseline->SetXTitle("Phi of the correlation (rad)");
    Yields_Central_1_MinusBaseline->SetYTitle("Yield_MinusBaseline");
    
    Yields_Periph_1_MinusBaseline = new TH1F("Yields_Periph_1_MinusBaseline",
                     "Yield of JPsi-tkl in Periph collisions wrt Phi _MinusBaseline",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yields_Periph_1_MinusBaseline->SetXTitle("Phi of the correlation (rad)");
    Yields_Periph_1_MinusBaseline->SetYTitle("Yield_MinusBaseline");
        
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
                    CorrelationsMEScaled[i][j][k] = new TH2F(hname,
                                                     "MixedEvent Correlation delta eta wrt delta phi - Scaled",
                                                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                                   CorrelationsMEScaled[i][j][k]->SetXTitle("Correlation DeltaPhi (rad)");
                                   CorrelationsMEScaled[i][j][k]->SetYTitle("Correlation DeltaEta");
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
                CorrelationsTklMEScaled[i][j] = new TH2F(hname,
                                  "MixedEvent Correlation Tkl delta eta wrt delta phi - Scaled",
                                  NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
                CorrelationsTklMEScaled[i][j]->SetXTitle("Correlation DeltaPhi (rad)");
                CorrelationsTklMEScaled[i][j]->SetYTitle("Correlation DeltaEta");
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
    
    BaselineTkl_Central = new TH1F("BaselineTkl_Central",
                     "Baseline of tkl-tkl in Central collisions wrt Phi",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    BaselineTkl_Central->SetXTitle("Phi of the correlation (rad)");
    BaselineTkl_Central->SetYTitle("Baseline");
    
    BaselineTkl_Periph = new TH1F("BaselineTkl_Periph",
                     "Baseline of tkl-tkl in Periph collisions wrt Phi",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    BaselineTkl_Periph->SetXTitle("Phi of the correlation (rad)");
    BaselineTkl_Periph->SetYTitle("Baseline");
    
    YieldTkl_Central_MinusBaseline = new TH1F("YieldTkl_Central_MinusBaseline",
                     "Yield of tkl-tkl in Central collisions wrt Phi _MinusBaseline",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_Central_MinusBaseline->SetXTitle("Phi of the correlation (rad)");
    YieldTkl_Central_MinusBaseline->SetYTitle("Yield_MinusBaseline");
    
    YieldTkl_Periph_MinusBaseline = new TH1F("YieldTkl_Periph_MinusBaseline",
                     "Yield of tkl-tkl in Periph collisions wrt Phi _MinusBaseline",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_Periph_MinusBaseline->SetXTitle("Phi of the correlation (rad)");
    YieldTkl_Periph_MinusBaseline->SetYTitle("Yield_MinusBaseline");
        
        
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
    
    
    
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for(int i=0; i<NbBinsCent; i++){
               for(int k=0; k<NbinsInvMass; k++){
                   char hname[100];
                   sprintf(hname,"Yields %d %d %d ",i,k,ptbin);
                   YieldsPtBinned[i][k][ptbin] = new TH1F(hname,
                                     "YieldsPtBinned Correlation wrt delta phi",
                                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
                   YieldsPtBinned[i][k][ptbin]->SetXTitle("Correlation DeltaPhi (rad)");
                   for(int j=0; j<NbBinsZvtx; j++){
                       sprintf(hname,"CorrelationsPtBinned %d %d %d %d ",i,j,k,ptbin);
                       CorrelationsPtBinned[i][j][k][ptbin] = new TH2F(hname,
                                         "CorrelationPtBinned delta eta wrt delta phi",
                                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                       CorrelationsPtBinned[i][j][k][ptbin]->SetXTitle("Correlation DeltaPhi (rad)");
                       CorrelationsPtBinned[i][j][k][ptbin]->SetYTitle("Correlation DeltaEta");
                       sprintf(hname,"CorrelationsMEPtBinned %d %d %d %d ",i,j,k,ptbin);
                       CorrelationsMEPtBinned[i][j][k][ptbin] = new TH2I(hname,
                                         "MixedEvent CorrelationPtBinned delta eta wrt delta phi",
                                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                       CorrelationsMEPtBinned[i][j][k][ptbin]->SetXTitle("Correlation DeltaPhi (rad)");
                       CorrelationsMEPtBinned[i][j][k][ptbin]->SetYTitle("Correlation DeltaEta");
                       
                       sprintf(hname,"CorrelationsMEScaledPtBinned %d %d %d %d ",i,j,k,ptbin);
                       CorrelationsMEScaledPtBinned[i][j][k][ptbin] = new TH2F(hname,
                                         "MixedEvent CorrelationPtBinned delta eta wrt delta phi - Scaled",
                                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                       CorrelationsMEScaledPtBinned[i][j][k][ptbin]->SetXTitle("Correlation DeltaPhi (rad)");
                       CorrelationsMEScaledPtBinned[i][j][k][ptbin]->SetYTitle("Correlation DeltaEta");
                   }
               }
               
           }
        
        for (int j=0; j <NbinsInvMass; j++){
        char hname[100];
        sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central - PtBinned %d",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1),ptbin);
        Yield_Central_MassBinPtBinned[j][ptbin] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        }
        
        for (int j=0; j <NbinsInvMass; j++){
        char hname[100];
        sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph - PtBinned %d",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1),ptbin);
        Yield_Periph_MassBinPtBinned[j][ptbin] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        }
        
        for (int j=0; j <NbinsInvMass; j++){
        char hname[100];
        sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference - PtBinned %d",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1),ptbin);
        Yield_Difference_MassBinPtBinned[j][ptbin] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        }
        
        for(int i=0; i<NbBinsCent; i++){
                  for(int p=0; p<NbinsDeltaPhi; p++){
                      char hname[100];
                      sprintf(hname,"Yields PhiBinPtBinned %d %d %d",i,p,ptbin);
                      Yields_PhiBinPtBinned[i][p][ptbin] = new TH1F(hname,
                                        "Yields CorrelationPtBinned wrt mass",
                                        NbinsInvMass,MinInvMass,MaxInvMass);
                      Yields_PhiBinPtBinned[i][p][ptbin]->SetXTitle("Correlation Inv Mass (GeV)");
              }
           }
        
        for(int p=0; p<NbinsDeltaPhi; p++){
                   char hname[100];
                   char hname2[100];
                   sprintf(hname,"Yields PhiBinPtBinned %d All C %d",p,ptbin);
                   YieldWrtMass_allCPtBinned[p][ptbin] = new TH1F(hname,
                                     "Yields CorrelationPtBinned wrt mass, all C",
                                     NbinsInvMass,MinInvMass,MaxInvMass);
                   YieldWrtMass_allCPtBinned[p][ptbin]->SetXTitle("Correlation Inv Mass (GeV)");
               
               sprintf(hname,"Yields PhiBinPtBinned %d Periph %d",p,ptbin);
               sprintf(hname2,"Yields CorrelationPtBinned wrt mass, Periph, Phi: %d pi/6 to %d pi/6",p-3, p-2);
               YieldWrtMass_PeriphPtBinned[p][ptbin] = new TH1F(hname,
                                 hname2,
                                 NbinsInvMass,MinInvMass,MaxInvMass);
               YieldWrtMass_PeriphPtBinned[p][ptbin]->SetXTitle("Correlation Inv Mass (GeV)");
               
               sprintf(hname,"Yields PhiBinPtBinned %d Central %d",p,ptbin);
               sprintf(hname2,"Yields CorrelationPtBinned wrt mass, Central, Phi: %d pi/6 to %d pi/6",p-3, p-2);
               YieldWrtMass_CentralPtBinned[p][ptbin] = new TH1F(hname,
                                 hname2,
                                 NbinsInvMass,MinInvMass,MaxInvMass);
               YieldWrtMass_CentralPtBinned[p][ptbin]->SetXTitle("Correlation Inv Mass (GeV)");
           }
        
        char hname[100];
        char hname2[100];
        sprintf(hname,"Yields_Central_1PtBinned %d",ptbin);
        sprintf(hname2,"Yield of JPsi-tkl in Central collisions wrt Phi - PtBinned %d",ptbin);
        Yields_Central_1PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Central_1PtBinned[ptbin]->SetXTitle("Phi of the correlation (rad)");
        Yields_Central_1PtBinned[ptbin]->SetYTitle("Yield");
        
        sprintf(hname,"Yields_Periph_1PtBinned %d",ptbin);
        sprintf(hname2,"Yield of JPsi-tkl in Periph collisions wrt Phi - PtBinned %d",ptbin);
        Yields_Periph_1PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Periph_1PtBinned[ptbin]->SetXTitle("Phi of the correlation (rad)");
        Yields_Periph_1PtBinned[ptbin]->SetYTitle("Yield");
        
        sprintf(hname,"Yields_Difference_1PtBinned %d",ptbin);
        sprintf(hname2,"Yield of JPsi-tkl in (Central-Periph) collisions wrt Phi - PtBinned %d",ptbin);
        Yields_Difference_1PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Difference_1PtBinned[ptbin]->SetXTitle("Phi of the correlation (rad)");
        Yields_Difference_1PtBinned[ptbin]->SetYTitle("Yield");
        
        sprintf(hname,"Baseline_Central_1PtBinned %d",ptbin);
               sprintf(hname2,"Baseline of JPsi-tkl in Central collisions wrt Phi - PtBinned %d",ptbin);
        Baseline_Central_1PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Baseline_Central_1PtBinned[ptbin]->SetXTitle("Phi of the correlation (rad)");
           Baseline_Central_1PtBinned[ptbin]->SetYTitle("Baseline");
           
           sprintf(hname,"Baseline_Periph_1PtBinned %d",ptbin);
                  sprintf(hname2,"Baseline of JPsi-tkl in Periph collisions wrt Phi - PtBinned %d",ptbin);
           Baseline_Periph_1PtBinned[ptbin] = new TH1F(hname,
                               hname2,
                               NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
              Baseline_Periph_1PtBinned[ptbin]->SetXTitle("Phi of the correlation (rad)");
              Baseline_Periph_1PtBinned[ptbin]->SetYTitle("Baseline");
        
        sprintf(hname,"Yields_Central_1_MinusBaselinePtBinned %d",ptbin);
        sprintf(hname2,"Yield of JPsi-tkl in Central collisions wrt Phi _MinusBaseline - PtBinned %d",ptbin);
        Yields_Central_1_MinusBaselinePtBinned[ptbin] = new TH1F(hname,hname2,
                            NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Yields_Central_1_MinusBaselinePtBinned[ptbin]->SetXTitle("Phi of the correlation (rad)");
           Yields_Central_1_MinusBaselinePtBinned[ptbin]->SetYTitle("Yield_MinusBaseline");
        
        sprintf(hname,"Yields_Periph_1_MinusBaselinePtBinned %d",ptbin);
        sprintf(hname2,"Yield of JPsi-tkl in Periph collisions wrt Phi _MinusBaseline - PtBinned %d",ptbin);
        Yields_Periph_1_MinusBaselinePtBinned[ptbin] = new TH1F(hname,hname2,
                            NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Yields_Periph_1_MinusBaselinePtBinned[ptbin]->SetXTitle("Phi of the correlation (rad)");
           Yields_Periph_1_MinusBaselinePtBinned[ptbin]->SetYTitle("Yield_MinusBaseline");
        
        sprintf(hname,"coefficients0PtBinned %d",ptbin);
        sprintf(hname2,"Coefficient Fourrier 0 - PtBinned %d",ptbin);
        coefficients0PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsInvMass,MinInvMass,MaxInvMass);
           coefficients0PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients0PtBinned[ptbin]->SetYTitle("Coefficient 0 (Fourier)");
        
        sprintf(hname,"coefficients1PtBinned %d",ptbin);
        sprintf(hname2,"Coefficient Fourrier 1 - PtBinned %d",ptbin);
        coefficients1PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsInvMass,MinInvMass,MaxInvMass);
           coefficients1PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients1PtBinned[ptbin]->SetYTitle("Coefficient 1 (Fourier)");
        
        sprintf(hname,"coefficients2PtBinned %d",ptbin);
        sprintf(hname2,"Coefficient Fourrier 2 - PtBinned %d",ptbin);
        coefficients2PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsInvMass,MinInvMass,MaxInvMass);
           coefficients2PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients2PtBinned[ptbin]->SetYTitle("Coefficient 2 (Fourier)");
        
        sprintf(hname,"baselines0PtBinned %d",ptbin);
        sprintf(hname2,"Baseline 0 - PtBinned %d",ptbin);
        baselines0PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsInvMass,MinInvMass,MaxInvMass);
        baselines0PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
        baselines0PtBinned[ptbin]->SetYTitle("Baseline 0 (YieldPeriph deltaPhi=0)");
        
        sprintf(hname,"c2b0PtBinned %d",ptbin);
        sprintf(hname2,"Coefficient 2 + Baseline 0 - PtBinned %d",ptbin);
        c2b0PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsInvMass,MinInvMass,MaxInvMass);
        c2b0PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
        c2b0PtBinned[ptbin]->SetYTitle("Coefficient 2 + Baseline 0");
        
        sprintf(hname,"V2JPsiTklPtBinned %d",ptbin);
        sprintf(hname2,"V2JPsiTkl wrt Mass - PtBinned %d",ptbin);
        V2JPsiTklPtBinned[ptbin] = new TH1F(hname,
                               hname2,
                               NbinsInvMass,MinInvMass,MaxInvMass);
              V2JPsiTklPtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
              V2JPsiTklPtBinned[ptbin]->SetYTitle("V2JPsiTkl");
    }
    
    
    
    
    
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
        int DimuonCounter[NbBinsCent][NbBinsZvtx][NbinsInvMass] = {0};
        int DimuonCounterZint[NbBinsCent][NbinsInvMass] = {0};
        int RefTklCounter[NbBinsCent][NbBinsZvtx] = {0};
        int RefTklCounterZint[NbBinsCent] = {0};
        double NormME[NbBinsCent][NbBinsZvtx][NbinsInvMass] = {0};
        double NormMETkl[NbBinsCent][NbBinsZvtx] = {0};
        
    //For PtBinning
       int DimuSeenMassCutPtBinned[NbPtBins] = {0};
       int DimuonCounterPtBinned[NbBinsCent][NbBinsZvtx][NbinsInvMass][NbPtBins] = {0};
       int DimuonCounterZintPtBinned[NbBinsCent][NbinsInvMass][NbPtBins] = {0};
       double NormMEPtBinned[NbBinsCent][NbBinsZvtx][NbinsInvMass][NbPtBins] = {0};
    
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
    
    YieldTkl_Central = (TH1F*)filerec->Get("YieldTkl_Central");
    YieldTkl_Periph = (TH1F*)filerec->Get("YieldTkl_Periph");
    YieldTkl_Difference = (TH1F*)filerec->Get("YieldTkl_Difference");
    
    
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
    
    
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for(int j=0; j<NbinsInvMass; j++){
           sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central - PtBinned %d",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1),ptbin);
            Yield_Central_MassBinPtBinned[j][ptbin] = (TH1F*)filerec->Get(hname);
            
             sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph - PtBinned %d",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1),ptbin);
                   Yield_Periph_MassBinPtBinned[j][ptbin] = (TH1F*)filerec->Get(hname);
           
            sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference - PtBinned %d",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1),ptbin);
            Yield_Difference_MassBinPtBinned[j][ptbin] = (TH1F*)filerec->Get(hname);

        }
    }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for(int p=0; p<NbinsDeltaPhi; p++){
             sprintf(hname,"Yields PhiBinPtBinned %d Periph %d",p,ptbin);
                          sprintf(hname2,"Yields CorrelationPtBinned wrt mass, Periph, Phi: %d pi/6 to %d pi/6",p-3, p-2);
                          YieldWrtMass_PeriphPtBinned[p][ptbin] = (TH1F*)filerec->Get(hname);
                          
                          sprintf(hname,"Yields PhiBinPtBinned %d Central %d",p,ptbin);
                          sprintf(hname2,"Yields CorrelationPtBinned wrt mass, Central, Phi: %d pi/6 to %d pi/6",p-3, p-2);
                          YieldWrtMass_CentralPtBinned[p][ptbin] = (TH1F*)filerec->Get(hname);

        }
    }
    
    for(int j=0;j<3;j++){
        for (int i=0;i<NbPtBins;i++){
            sprintf(hname,"bin%d_%d",i+1, j);
            hPtWrtMassInvSliced[j][i] = (TH1F*)filerec->Get(hname);
        }
    }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        sprintf(hname,"coefficients0PtBinned %d",ptbin);
        coefficients0PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"coefficients1PtBinned %d",ptbin);
        coefficients1PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"coefficients2PtBinned %d",ptbin);
        coefficients2PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"baselines0PtBinned %d",ptbin);
        baselines0PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"V2JPsiTklPtBinned %d",ptbin);
        V2JPsiTklPtBinned[ptbin] = (TH1F*)filerec->Get(hname);
    }
    
    TCanvas*c16PtBinned=new TCanvas();
    // ROOT THAT Extraction 2 plot
    c16PtBinned->SetTitle("Extraction method 2 PtBinned");
    c16PtBinned->Divide(2,3);
    //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        c16PtBinned->cd(ptbin+1);
        c2b0PtBinned[ptbin]->Add(coefficients0PtBinned[ptbin], baselines0PtBinned[ptbin]);
        V2JPsiTklPtBinned[ptbin]->Divide(coefficients2PtBinned[ptbin], c2b0PtBinned[ptbin]);
        
        V2JPsiTklPtBinned[ptbin]->Draw();
    }
    
    
    
    
    cout << "START ALL"<<endl;
    
    
    
    
    
    
    
    
    
    
    
    TFitResultPtr res;
    TFitResultPtr rescent;
    TFitResultPtr resperiph;
    TFitResultPtr resu;
    
   ROOT::Math::Minimizer* minim{NULL};
    
        char histoname[50];
        //Values OK for 20M evts
    //    double startvalues[16] = {800, 3.096916, 0.07, 0.9, 10, 2, 15, 15,   10,     0.1,          50000,   0.01,        1,         1,   1,   1};
    //    double lowvalues[16] =   {1,    3.05,    0.03, 0.7, 1,  1, 5,  0.01, 0.0001,  0.000000001,  10,    0.0000000001, 0.000001, -10, -10, -10};
    //    double highvalues[16] = {100000, 3.15,   0.1,  1.1, 30, 5, 25, 50,   100000, 10,           100000, 10,     10,        10,  10,  10};
        
        double startvalues[16] = {1100, 3.096916, 0.07,  0.9, 10, 2, 15, 20,    90,      0.01,     80000,      1.4,      1,         1,   1,   1};
        double lowvalues[16] =   {10,    3.05,     0.03, 0.7, 1,  1, 5,  0.001,  1,    0.00001,   1000,        0.0001, 0.00000001, -10, -10, -10};
        double highvalues[16] = {100000, 3.15,     1,  1.1, 30, 5, 25, 10000, 1000000, 100,      10000000,  100,    10,        10,  10,  10};
        
    //    for(int index=0; index<16; index++){
    //        if(index==0 || index==7 || index==8 || index==10){
    //            startvalues[index]*=nevent1/20000000;
    //            lowvalues[index]*=nevent1/20000000;
    //            highvalues[index]*=nevent1/20000000;
    //        }
    //    }
        //Values adapted to the number of events considered
        sprintf(histoname,"hnseg");
        res = FittingAllInvMassBin(histoname, cinvmass, 0);
        double par[16];
        
        // Parameter setting (either from mass fit either start values)
        for(int i=0; i <16; i++){
                par[i] = res->Parameter(i);
                if(i>11){
                    par[i] = 1;
                }
        }
    
            TCanvas*c16=new TCanvas();
           // ROOT THAT Extraction 2 plot
           c16->SetTitle("Extraction method 2");
           //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
           c2b0->Add(coefficients0, baselines0);
           V2JPsiTkl->Divide(coefficients2, c2b0);
           
           V2JPsiTkl->Draw();
        
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
        TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2Poly,2.1,5.1,16);
        backFcnV2_2->SetLineColor(kRed);
        TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,2.1,5.1,16);
        signalFcnJPsiV2_2->SetLineColor(kBlue);
        signalFcnJPsiV2_2->SetNpx(500);
        
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
        
        fitJustV2_2->SetParName(0,"V2_2 JPsi");
            fitJustV2_2->SetParName(1,"V2_2 Bkg M2");
        fitJustV2_2->SetParName(2,"V2_2 Bkg M1");
        fitJustV2_2->SetParName(3,"V2_2 Nkg M0");
        
           double minParams[16];
           double parErrors[16];
        
        //Fit of V2
            res = V2JPsiTkl->Fit("fitV2_2","SBMERI+","ep");
            Double_t para[16];
            fitV2_2->GetParameters(para);
            backFcnV2_2->SetParameters(para);
            signalFcnJPsiV2_2->SetParameters(para);
        
        
        fitV2_2->Draw("same");
      //  signalFcnJPsiV2_2->Draw("same");
        backFcnV2_2->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
          legend->SetTextFont(61);
          legend->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
         legend->AddEntry(fitV2_2,message);

            if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
            }
            else{
                   sprintf(message,"The fit is a failure");
            }
        
               legend->AddEntry(fitV2_2,message);
       // legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
        legend->AddEntry(backFcnV2_2,"Background");
          legend->AddEntry(V2JPsiTkl,"Data","lpe");
          legend->Draw();
    
    
    //Ext 2 Pt Binned
    
     TCanvas *cinvmassPtBinned = new TCanvas("cinvmassPtBinned","Fitting All Phi - Different Centralities - Pt Binned",10,10,700,500);
        cinvmassPtBinned->SetTitle("Inv Mass Fits - Pt Binned");
        // ROOT THAT Inv mass fits
        cinvmassPtBinned->Divide(NbPtBins,3);
        
    if(PtBinned){
        for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    
        sprintf(histoname,Form("bin%d_0",ptbin+1)); // ZZZZZZZ
        res = FittingAllInvMassBin(histoname, cinvmassPtBinned, ptbin);
        double par[16];
        
        // Parameter setting (either from mass fit either start values)
        for(int i=0; i <16; i++){
                par[i] = res->Parameter(i);
                if(i>11){
                    par[i] = 1;
                }
        }
        
        //Definition and design of function on V2 plot
        c16PtBinned->cd(ptbin+1);
        TF1 *fitV2_2PtBinned = new TF1("fitV2_2PtBinned",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
        fitV2_2PtBinned->SetNpx(500);
        fitV2_2PtBinned->SetLineWidth(4);
        fitV2_2PtBinned->SetLineColor(kMagenta);
        fitV2_2PtBinned->SetParameters(par);
        TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2Poly,2.1,5.1,16);
        backFcnV2_2->SetLineColor(kRed);
        TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,2.1,5.1,16);
        signalFcnJPsiV2_2->SetLineColor(kBlue);
        signalFcnJPsiV2_2->SetNpx(500);
        
           TVirtualFitter::Fitter(V2JPsiTklPtBinned[ptbin])->SetMaxIterations(10000);
           TVirtualFitter::Fitter(V2JPsiTklPtBinned[ptbin])->SetPrecision();
            for(int i=0; i<=11; i++){
                fitV2_2PtBinned->FixParameter(i,par[i]);
            }
           
           fitV2_2PtBinned->SetParName(0,"Norm_{JPsi}");
           fitV2_2PtBinned->SetParName(1,"M_{JPsi}");
           fitV2_2PtBinned->SetParName(2,"Sigma_{JPsi}");
           fitV2_2PtBinned->SetParName(3,"a_{1}");
           fitV2_2PtBinned->SetParName(4,"n_{1}");
           fitV2_2PtBinned->SetParName(5,"a_{2}");
           fitV2_2PtBinned->SetParName(6,"n_{2}");
           fitV2_2PtBinned->SetParName(7,"Norm_{Psi2S}");
           fitV2_2PtBinned->SetParName(8,"Norm_{TailLowM}");
           fitV2_2PtBinned->SetParName(9,"Exp_{TailLowM}");
           fitV2_2PtBinned->SetParName(10,"Norm_{TailHighM}");
           fitV2_2PtBinned->SetParName(11,"Exp_{TailHighM}");
            fitV2_2PtBinned->SetParName(12,"V2_2 JPsi");
            fitV2_2PtBinned->SetParName(13,"V2_2 Bkg M2");
        fitV2_2PtBinned->SetParName(14,"V2_2 Bkg M1");
        fitV2_2PtBinned->SetParName(15,"V2_2 Nkg M0");

        
           double minParams[16];
           double parErrors[16];
        
        //Fit of V2
            res = V2JPsiTklPtBinned[ptbin]->Fit("fitV2_2PtBinned","SBMERI+","ep");
            Double_t para[16];
            fitV2_2PtBinned->GetParameters(para);
            backFcnV2_2->SetParameters(para);
            signalFcnJPsiV2_2->SetParameters(para);
        
        
        fitV2_2PtBinned->Draw("same");
      //  signalFcnJPsiV2_2->Draw("same");
        backFcnV2_2->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
          legend->SetTextFont(61);
          legend->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2PtBinned->GetChisquare(),fitV2_2PtBinned->GetNDF());
         legend->AddEntry(fitV2_2,message);

            if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
            }
            else{
                   sprintf(message,"The fit is a failure");
            }
        
               legend->AddEntry(fitV2_2PtBinned,message);
       // legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
        legend->AddEntry(backFcnV2_2,"Background");
          legend->AddEntry(V2JPsiTklPtBinned[ptbin],"Data","lpe");
          legend->Draw();
            
            V2_Ext2[ptbin] = para[12];
            errV2_Ext2[ptbin] = fitV2_2PtBinned->GetParError(12);
        }
    }
    
    
    
    
    
    
    
    //Combining or not fit for Central-Periph method Ext1
    
    
    int niterations = 1;
    cinvmass->cd();
    sprintf(histoname,"InvMass_Central");
   res = FittingAllInvMassBin(histoname, cinvmass, 1); //ZZZZZZZ

    TCanvas* chnsegsigma = new TCanvas;
    chnsegsigma->cd();
    hnsegSigma->Draw("colz");

    double parerr[12];
    TCanvas* c10_fit = new TCanvas;
       c10_fit->Divide(6,4);
    for(int j=0; j<12; j++){
        parerr[j] = res->ParError(j);
    }
    for(int j=0; j <16; j++){
        par[j] = res->Parameter(j);
        if(j>11){
            par[j] = 3;
        }
    }

    double par12[niterations];
    double parerr12[niterations];
    
    Char_t fileLogName[200] = "Log_FitTraining.txt";
    std::ofstream fileLog(fileLogName, std::ofstream::out);
    
    for(int i=0; i<NbinsDeltaPhi; i++){ //Change
        fileLog << "Bin " << i << " in DeltaPhi Central" <<endl;
        
        for(int t=0; t<niterations; t++){
            fileLog << "Interation " << t <<endl;
        TVirtualPad* subpad = c10_fit->cd(i+1);
            c10_fit->cd(i+1);
            if(t==niterations-1){
                subpad->Clear();
            }
        TF1 *fitY_1Central = new TF1("fitY_1Central",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
          fitY_1Central->SetNpx(500);
          fitY_1Central->SetLineWidth(4);
          fitY_1Central->SetLineColor(kMagenta);
          fitY_1Central->SetParameters(par);
           TVirtualFitter::Fitter(YieldWrtMass_Central[i])->SetMaxIterations(10000);
           TVirtualFitter::Fitter(YieldWrtMass_Central[i])->SetPrecision();
        //  histo->Fit("fitFcn","0");
          // second try: set start values for some parameters
        for(int k=0; k<=11; k++){
                fitY_1Central->FixParameter(k,par[k]);
        }
            fitY_1Central->SetParLimits(12,0,10000);

           fitY_1Central->SetParName(0,"Norm_{JPsi}");
           fitY_1Central->SetParName(1,"M_{JPsi}");
           fitY_1Central->SetParName(2,"Sigma_{JPsi}");
           fitY_1Central->SetParName(3,"a_{1}");
           fitY_1Central->SetParName(4,"n_{1}");
           fitY_1Central->SetParName(5,"a_{2}");
           fitY_1Central->SetParName(6,"n_{2}");
           fitY_1Central->SetParName(7,"Norm_{Psi2S}");
           fitY_1Central->SetParName(8,"Norm_{TailLowM}");
           fitY_1Central->SetParName(9,"Exp_{TailLowM}");
           fitY_1Central->SetParName(10,"Norm_{TailHighM}");
           fitY_1Central->SetParName(11,"Exp_{TailHighM}");
            fitY_1Central->SetParName(12,"Y_1 JPsi");
            fitY_1Central->SetParName(13,"Y_1 Bkg M2");
        fitY_1Central->SetParName(14,"Y_1 Bkg M1");
        fitY_1Central->SetParName(15,"Y_1 Nkg M0");

          rescent = YieldWrtMass_Central[i]->Fit("fitY_1Central","SMERIQ+","ep");
            rescent->Print("V");
            fileLog << "Chi2: " << rescent->Chi2() <<endl;
            fileLog << "Status: " << int(rescent->Status()) <<endl;
            fileLog << "Prob: " << rescent->Prob() <<endl;
            fileLog << "IsValid: " << rescent->IsValid() <<endl;
            fileLog << "HasMinosError(12): " << rescent->HasMinosError(12) <<endl;
            
          Double_t param[16];
            Double_t paramerrs[16];
          fitY_1Central->GetParameters(param);
            for(int i=0; i<16; i++){
                paramerrs[i] = fitY_1Central->GetParError(i);
            }

            if(rescent->CovMatrixStatus() == 3){
            par12[t] = param[12];
            parerr12[t] = paramerrs[12];
                YieldValue->Fill(param[12]);
                YieldError->Fill(paramerrs[12]);
            }
            else{
                par12[t]=parerr12[t]=0;
            }

            std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
            std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;
            
            fileLog << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
            fileLog << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

        TF1 *backFcnY_1Central = new TF1("backFcnY_1Central",BackFcnV2Poly,2.1,5.1,16);
        backFcnY_1Central->SetLineColor(kRed);
        TF1 *signalFcnJPsiY_1Central = new TF1("signalFcnJPsiY_1Central",SignalFcnJPsiV2,2.1,5.1,16);
          // writes the fit results into the par array

            if(t==niterations-1){
                YieldWrtMass_Central[i]->GetListOfFunctions()->Clear();
                YieldWrtMass_Central[i]->GetListOfFunctions()->Add(fitY_1Central);
            }
        signalFcnJPsiY_1Central->SetLineColor(kBlue);
        signalFcnJPsiY_1Central->SetNpx(500);
        backFcnY_1Central->SetParameters(param);
        signalFcnJPsiY_1Central->SetParameters(param);
       // signalFcnJPsiY_1Central->Draw("same");
        backFcnY_1Central->Draw("same");

//                if(i==2){
//                    baseline_central = param[12];
//                    errbaseline_central = fitY_1Central->GetParError(12);
//                }
//                if(i==3){
//                    baseline_central += param[12];
//                    baseline_central /= 2;
//                    errbaseline_central = sqrt(pow(errbaseline_central,2)+pow(fitY_1Central->GetParError(12),2));
//                }
                
                if(param[12]<baseline_central){ //En central baseline = minimum
                    baseline_central = param[12];
                    errbaseline_central = fitY_1Central->GetParError(12);
                }
                
                Yields_Central_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[12]);
                Yields_Central_1->SetBinError(i+1,fitY_1Central->GetParError(12));
            

          // writes the fit results into the par array
           gStyle->SetOptFit(1011);
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetTextFont(72);
          legend->SetTextSize(0.04);
        Char_t message[80];
        sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
         legend->AddEntry(fitY_1Central,message);
        if(rescent->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitY_1Central,message);
      //  legend->AddEntry(signalFcnJPsiY_1Central,"JPsi signal");
        legend->AddEntry(backFcnY_1Central,"Background");
          legend->AddEntry(YieldWrtMass_Central[i],"Data","lpe");
          legend->Draw();
        }
        
//        TCanvas* cyo = new TCanvas;
//        cyo->Divide(2,1);
//        cyo->cd(1);
//        YieldValue->Draw();
//        cyo->cd(2);
//        YieldError->Draw();
        
    }
    
    for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
        Baseline_Central_1->SetBinContent(phi_idx+1, baseline_central);
        Baseline_Central_1->SetBinError(phi_idx+1, errbaseline_central);
    }
    
    Yields_Central_1_MinusBaseline->Add(Yields_Central_1,Baseline_Central_1,1,-1);
    

        
    
    cinvmass->cd();
    sprintf(histoname,"InvMass_Periph");
    res = FittingAllInvMassBin(histoname, cinvmass, 2); //ZZZZZ

       for(int j=0; j<12; j++){
           parerr[j] = res->ParError(j);
       }
       for(int j=0; j <16; j++){
           par[j] = res->Parameter(j);
           if(j>11){
               par[j] = 3;
           }
       }

//    double baseline = 1;
//    double errbaseline = 1;

       for(int i=0; i<NbinsDeltaPhi; i++){ //Change
           fileLog << "Bin " << i << " in DeltaPhi Periph" <<endl;
           for(int t=0; t<niterations; t++){
               fileLog << "Interation " << t <<endl;
               TVirtualPad* subpad = c10_fit->cd(NbinsDeltaPhi+i+1); //Change
               c10_fit->cd(NbinsDeltaPhi+i+1); //Change
               if(t==niterations-1){
                   subpad->Clear();
               }
           TF1 *fitY_1Periph = new TF1("fitY_1Periph",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
             fitY_1Periph->SetNpx(500);
             fitY_1Periph->SetLineWidth(4);
             fitY_1Periph->SetLineColor(kMagenta);
             fitY_1Periph->SetParameters(par);
              TVirtualFitter::Fitter(YieldWrtMass_Periph[i])->SetMaxIterations(10000);
              TVirtualFitter::Fitter(YieldWrtMass_Periph[i])->SetPrecision();
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters
            for(int k=0; k<=11; k++){
                    fitY_1Periph->FixParameter(k,par[k]);
            }
              fitY_1Periph->SetParLimits(12,0,10000);

              fitY_1Periph->SetParName(0,"Norm_{JPsi}");
              fitY_1Periph->SetParName(1,"M_{JPsi}");
              fitY_1Periph->SetParName(2,"Sigma_{JPsi}");
              fitY_1Periph->SetParName(3,"a_{1}");
              fitY_1Periph->SetParName(4,"n_{1}");
              fitY_1Periph->SetParName(5,"a_{2}");
              fitY_1Periph->SetParName(6,"n_{2}");
              fitY_1Periph->SetParName(7,"Norm_{Psi2S}");
              fitY_1Periph->SetParName(8,"Norm_{TailLowM}");
              fitY_1Periph->SetParName(9,"Exp_{TailLowM}");
              fitY_1Periph->SetParName(10,"Norm_{TailHighM}");
              fitY_1Periph->SetParName(11,"Exp_{TailHighM}");
               fitY_1Periph->SetParName(12,"Y_1 JPsi");
               fitY_1Periph->SetParName(13,"Y_1 Bkg M2");
           fitY_1Periph->SetParName(14,"Y_1 Bkg M1");
           fitY_1Periph->SetParName(15,"Y_1 Nkg M0");

             resperiph = YieldWrtMass_Periph[i]->Fit("fitY_1Periph","SMERIQ+","ep");
               fileLog << "Chi2: " << resperiph->Chi2() <<endl;
               fileLog << "Status: " << int(resperiph->Status()) <<endl;
               fileLog << "Prob: " << resperiph->Prob() <<endl;
               fileLog << "IsValid: " << resperiph->IsValid() <<endl;
               fileLog << "HasMinosError(12): " << resperiph->HasMinosError(12) <<endl;
               Double_t param[16];
               Double_t paramerrs[16];
             fitY_1Periph->GetParameters(param);

               for(int i=0; i<16; i++){
                               paramerrs[i] = fitY_1Periph->GetParError(i);
                           }

                           if(resperiph->CovMatrixStatus() == 3){
                           par12[t] = param[12];
                           parerr12[t] = paramerrs[12];
                           }
                           else{
                               par12[t]=parerr12[t]=0;
                           }
               std::cout << "COV status Periph t=" << t << " : " << resperiph->CovMatrixStatus()<<endl;
               std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;
               
               fileLog << "COV status Periph t=" << t << " : " << resperiph->CovMatrixStatus()<<endl;
               fileLog << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

           TF1 *backFcnY_1Periph = new TF1("backFcnY_1Periph",BackFcnV2Poly,2.1,5.1,16);
           backFcnY_1Periph->SetLineColor(kRed);
           TF1 *signalFcnJPsiY_1Periph = new TF1("signalFcnJPsiY_1Periph",SignalFcnJPsiV2,2.1,5.1,16);
             // writes the fit results into the par array

               if(t==niterations-1){
                   YieldWrtMass_Periph[i]->GetListOfFunctions()->Clear();
                   YieldWrtMass_Periph[i]->GetListOfFunctions()->Add(fitY_1Periph);
               }

           signalFcnJPsiY_1Periph->SetLineColor(kBlue);
           signalFcnJPsiY_1Periph->SetNpx(500);
           backFcnY_1Periph->SetParameters(param);
           signalFcnJPsiY_1Periph->SetParameters(param);
         //  signalFcnJPsiY_1Periph->Draw("same");
           backFcnY_1Periph->Draw("same");


                   if(i==2){
                       baseline_periph = param[12];
                       errbaseline_periph = fitY_1Periph->GetParError(12);
                   }
                   if(i==3){
                       baseline_periph += param[12];
                       baseline_periph /= 2;
                       errbaseline_periph = sqrt(pow(errbaseline_periph,2)+pow(fitY_1Periph->GetParError(12),2));
                   }
                 Yields_Periph_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[12]);
                 Yields_Periph_1->SetBinError(i+1,fitY_1Periph->GetParError(12));
               

             // writes the fit results into the par array
              gStyle->SetOptFit(1011);
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
           Char_t message[80];
           sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
            legend->AddEntry(fitY_1Periph,message);
           if(resperiph->CovMatrixStatus() == 3){
                      sprintf(message,"The fit is a success");
                  }
                  else{
                      sprintf(message,"The fit is a failure");
                  }
                  legend->AddEntry(fitY_1Periph,message);
       //    legend->AddEntry(signalFcnJPsiY_1Periph,"JPsi signal");
           legend->AddEntry(backFcnY_1Periph,"Background");
             legend->AddEntry(YieldWrtMass_Periph[i],"Data","lpe");
             legend->Draw();
           }

       }
    
    for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
        Baseline_Periph_1->SetBinContent(phi_idx+1, baseline_periph);
        Baseline_Periph_1->SetBinError(phi_idx+1, errbaseline_periph);
    }
    
    Yields_Periph_1_MinusBaseline->Add(Yields_Periph_1,Baseline_Periph_1,1,-1);
    
    //Fill Baseline Periph and Y-B
    TCanvas* cYBcentral = new TCanvas;
    cYBcentral->SetTitle("Central Yield-Baseline");
    cYBcentral->Divide(1,3);
    cYBcentral->cd(1);
    Yields_Central_1->DrawCopy();
    cYBcentral->cd(2);
    Baseline_Central_1->DrawCopy();
    cYBcentral->cd(3);
    Yields_Central_1_MinusBaseline->DrawCopy();
    
    TCanvas* cYBperiph = new TCanvas;
    cYBperiph->SetTitle("Periph Yield-Baseline");
    cYBperiph->Divide(1,3);
    cYBperiph->cd(1);
    Yields_Periph_1->DrawCopy();
    cYBperiph->cd(2);
    Baseline_Periph_1->DrawCopy();
    cYBperiph->cd(3);
    Yields_Periph_1_MinusBaseline->DrawCopy();
    
    TH1F *Yields_Central_1_Cvetan = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_Cvetan");
    TH1F *Yields_Central_1_CvetanMe = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_CvetanMe");
    TH1F *Yields_Central_1_ZYAM = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_ZYAM");
    TH1F *Yields_Central_1_PRLTemplate = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_PRLTemplate");
    TH1F *Yields_Central_1_PRLTemplate_PeriphZYAM = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_PRLTemplate_PeriphZYAM");

    {
        TCanvas* c17 = new TCanvas;
        c17->Divide(1,3);
        //Creer canvas pour imprimer les plots Periph yield et Central yield wrt phi et leur difference
        c17->cd(1);
        Yields_Central_1->SetStats(kTRUE);
        Yields_Central_1->DrawCopy();
        c17->cd(2);
        Yields_Periph_1->SetStats(kTRUE);
        Yields_Periph_1->DrawCopy();
        c17->cd(3);
        Yields_Difference_1->Add(Yields_Central_1,Yields_Periph_1,1,-1);
        Yields_Difference_1->SetStats(kTRUE);
        Yields_Difference_1->DrawCopy();

             TF1 *fitFcnV2_2 = new TF1("fitFcnV2_2",FourierV2,-TMath::Pi()/2,1.5*TMath::Pi(),3);
             fitFcnV2_2->SetNpx(500);
             fitFcnV2_2->SetLineWidth(4);
             fitFcnV2_2->SetLineColor(kMagenta);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[3] = {1,1,1};
             fitFcnV2_2->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Difference_1)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Difference_1)->SetPrecision();
                gStyle->SetOptFit(1011);
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_2->SetParName(0,"a0");
              fitFcnV2_2->SetParName(1,"a1");
              fitFcnV2_2->SetParName(2,"a2");

             TFitResultPtr res = Yields_Difference_1->Fit("fitFcnV2_2","SBMERI+","ep");
            Double_t par[3];
            fitFcnV2_2->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_2->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
               Char_t message[80];
               sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
               legend->AddEntry(fitFcnV2_2,message);
        sprintf(message,"V2_1 JPsi-Tkl: %.4f / (%.4f + %.4f) = %.4f +- %.4f",par[2],par[0], baseline_periph,par[2]/(par[0] + baseline_periph),(par[2]/(par[0] + baseline_periph)*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(fitFcnV2_2->GetParError(0)/par[0],2)+pow(errbaseline_periph/baseline_periph,2))));
        legend->AddEntry(fitFcnV2_2,message);
        if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitFcnV2_2,message);
             legend->AddEntry(Yields_Difference_1,"Data","lpe");
             legend->Draw();

    }
    
    
    
    // Fit Cvetan
    
    TCanvas* cCvetan = new TCanvas;
       cCvetan->SetTitle("Cvetan fit");
       cCvetan->Divide(1,1);
       cCvetan->cd(1);
       Yields_Central_1->DrawCopy();
//       cCvetan->cd(2);
//       Baseline_Central_1->DrawCopy();
//       cCvetan->cd(3);
//       Yields_Central_1_MinusBaseline->DrawCopy();
    
    
    {
        cCvetan->cd(1);
             TF1 *fitFcnV2_Cvetan = new TF1("fitFcnV2_Cvetan",CvetanF,-TMath::Pi()/2,1.5*TMath::Pi(),4);
             fitFcnV2_Cvetan->SetNpx(500);
             fitFcnV2_Cvetan->SetLineWidth(4);
             fitFcnV2_Cvetan->SetLineColor(kBlue);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[4] = {1,0,0.01,1};
             fitFcnV2_Cvetan->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->SetPrecision();
                TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->SetFCN(ChisquareCvetanF);
            gStyle->SetOptFit(1011);
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_Cvetan->SetParName(0,"V0");
                fitFcnV2_Cvetan->SetParName(1,"V1");
                fitFcnV2_Cvetan->SetParName(2,"V2");
              fitFcnV2_Cvetan->SetParName(3,"F");
        
        fitFcnV2_Cvetan->SetParLimits(3,0.1,50);
        
        fitFcnV2_Cvetan->FixParameter(0,1);
        fitFcnV2_Cvetan->FixParameter(1,0);
      //  fitFcnV2_Cvetan->FixParameter(3,1);
              

             TFitResultPtr res = Yields_Central_1_Cvetan->Fit("fitFcnV2_Cvetan","USBMERI+","ep");
            Double_t par[4];
        
                double chi2, edm, errdef;
               int nvpar, nparx;
                TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->GetStats(chi2,edm,errdef,nvpar,nparx);
                fitFcnV2_Cvetan->SetChisquare(chi2);
                int ndf = npfits-nvpar;
                fitFcnV2_Cvetan->SetNDF(ndf);
        
            fitFcnV2_Cvetan->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_Cvetan->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
               Char_t message[80];
               sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_Cvetan->GetChisquare(),fitFcnV2_Cvetan->GetNDF());
               legend->AddEntry(fitFcnV2_Cvetan,message);
        if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitFcnV2_Cvetan,message);
             legend->AddEntry(Yields_Central_1_Cvetan,"Data","lpe");
             legend->Draw();

    }
    
    // Fit Cvetan putting my constraints (F=1, V1 free, V0 not 1)
    
    TCanvas* cCvetanMe = new TCanvas;
       cCvetanMe->SetTitle("Cvetan fit - My constraints F=1, v1 free, v0 free");
       cCvetanMe->Divide(1,1);
       cCvetanMe->cd(1);
       Yields_Central_1->DrawCopy();
//       cCvetanMe->cd(2);
//       Baseline_Central_1->DrawCopy();
//       cCvetanMe->cd(3);
//       Yields_Central_1_MinusBaseline->DrawCopy();
    
    
    {
        cCvetanMe->cd(1);
             TF1 *fitFcnV2_CvetanMe = new TF1("fitFcnV2_CvetanMe",CvetanF,-TMath::Pi()/2,1.5*TMath::Pi(),4);
             fitFcnV2_CvetanMe->SetNpx(500);
             fitFcnV2_CvetanMe->SetLineWidth(4);
             fitFcnV2_CvetanMe->SetLineColor(kBlue);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[4] = {1,0,0.01,1};
             fitFcnV2_CvetanMe->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->SetPrecision();
            TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->SetFCN(ChisquareCvetanF);
            gStyle->SetOptFit(1011);
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_CvetanMe->SetParName(0,"V0");
                fitFcnV2_CvetanMe->SetParName(1,"V1");
                fitFcnV2_CvetanMe->SetParName(2,"V2");
              fitFcnV2_CvetanMe->SetParName(3,"F");
        
        fitFcnV2_CvetanMe->SetParLimits(3,0.1,50);
        
     //   fitFcnV2_CvetanMe->FixParameter(0,1);
       // fitFcnV2_CvetanMe->FixParameter(1,0);
        fitFcnV2_CvetanMe->FixParameter(3,1);
              

             TFitResultPtr res = Yields_Central_1_CvetanMe->Fit("fitFcnV2_CvetanMe","USBMERI+","ep");
        
        double chi2, edm, errdef;
        int nvpar, nparx;
         TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->GetStats(chi2,edm,errdef,nvpar,nparx);
         fitFcnV2_CvetanMe->SetChisquare(chi2);
         int ndf = npfits-nvpar;
         fitFcnV2_CvetanMe->SetNDF(ndf);
        
            Double_t par[4];
            fitFcnV2_CvetanMe->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_CvetanMe->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
               Char_t message[80];
               sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_CvetanMe->GetChisquare(),fitFcnV2_CvetanMe->GetNDF());
               legend->AddEntry(fitFcnV2_CvetanMe,message);
        if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitFcnV2_CvetanMe,message);
             legend->AddEntry(Yields_Central_1_CvetanMe,"Data","lpe");
             legend->Draw();

    }
    
    
    // Fit ZYAM Yc = Bc + (Yp-Bp) + a0  + 2a1cos + 2a2cos
    
    
     TCanvas* cZYAM = new TCanvas;
        cZYAM->SetTitle("ZYAM fit");
        cZYAM->Divide(1,1);
        cZYAM->cd(1);
        Yields_Central_1->DrawCopy();
    //    cZYAM->cd(2);
    //    Baseline_Central_1->DrawCopy();
    //    cZYAM->cd(3);
    //    Yields_Central_1_MinusBaseline->DrawCopy();
    
    
    {
        cZYAM->cd(1);
             TF1 *fitFcnV2_ZYAM = new TF1("fitFcnV2_ZYAM",ZYAM,-TMath::Pi()/2,1.5*TMath::Pi(),3);
             fitFcnV2_ZYAM->SetNpx(500);
             fitFcnV2_ZYAM->SetLineWidth(4);
             fitFcnV2_ZYAM->SetLineColor(kGreen);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[3] = {1,1,1};
             fitFcnV2_ZYAM->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->SetPrecision();
            TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->SetFCN(ChisquareZYAM);
            gStyle->SetOptFit(1011);
            //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_ZYAM->SetParName(0,"a0");
                fitFcnV2_ZYAM->SetParName(1,"a1");
                fitFcnV2_ZYAM->SetParName(2,"a2");
            
              

             TFitResultPtr res = Yields_Central_1_ZYAM->Fit("fitFcnV2_ZYAM","USBMERI+","ep");
        
            double chi2, edm, errdef;
            int nvpar, nparx;
             TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->GetStats(chi2,edm,errdef,nvpar,nparx);
             fitFcnV2_ZYAM->SetChisquare(chi2);
             int ndf = npfits-nvpar;
             fitFcnV2_ZYAM->SetNDF(ndf);
        
            Double_t par[3];
            fitFcnV2_ZYAM->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_ZYAM->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
               Char_t message[80];
               sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_ZYAM->GetChisquare(),fitFcnV2_ZYAM->GetNDF());
               legend->AddEntry(fitFcnV2_ZYAM,message);
        if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitFcnV2_ZYAM,message);
             legend->AddEntry(Yields_Central_1_ZYAM,"Data","lpe");
             legend->Draw();

    }
    
    // Fit PRL Template
    
    TCanvas* cPRLTemplate = new TCanvas;
       cPRLTemplate->SetTitle("PRL Template fit");
       cPRLTemplate->Divide(1,1);
       cPRLTemplate->cd(1);
       Yields_Central_1->DrawCopy();
//       cPRLTemplate->cd(2);
//       Baseline_Central_1->DrawCopy();
//       cPRLTemplate->cd(3);
//       Yields_Central_1_MinusBaseline->DrawCopy();
       
       
       {
           cPRLTemplate->cd(1);
                TF1 *fitFcnV2_PRLTemplate = new TF1("fitFcnV2_PRLTemplate",PRLTemplate,-TMath::Pi()/2,1.5*TMath::Pi(),2);
                fitFcnV2_PRLTemplate->SetNpx(500);
                fitFcnV2_PRLTemplate->SetLineWidth(4);
                fitFcnV2_PRLTemplate->SetLineColor(kRed);
           TF1 *fitFcnV2_PRLTemplate_RidgeAndZero = new TF1("fitFcnV2_PRLTemplate_RidgeAndZero",PRLTemplate_RidgeAndZero,-TMath::Pi()/2,1.5*TMath::Pi(),2);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetNpx(500);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetLineWidth(2);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetLineStyle(9);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetLineColor(kRed);
           TF1 *fitFcnV2_PRLTemplate_PeriphAndG = new TF1("fitFcnV2_PRLTemplate_PeriphAndG",PRLTemplate_PeriphAndG,-TMath::Pi()/2,1.5*TMath::Pi(),2);
           fitFcnV2_PRLTemplate_PeriphAndG->SetNpx(500);
           fitFcnV2_PRLTemplate_PeriphAndG->SetLineWidth(2);
           fitFcnV2_PRLTemplate_PeriphAndG->SetLineStyle(9);
           fitFcnV2_PRLTemplate_PeriphAndG->SetLineColor(kRed);
                // first try without starting values for the parameters
                // This defaults to 1 for each param.
                // this results in an ok fit for the polynomial function
                // however the non-linear part (lorenzian) does not
                // respond well.
                 Double_t params[2] = {1,1};
                fitFcnV2_PRLTemplate->SetParameters(params);
                 TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->SetMaxIterations(10000);
                 TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->SetPrecision();
                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->SetFCN(ChisquarePRLTemplate);
               gStyle->SetOptFit(1011);
               //  histo->Fit("fitFcn","0");
                // second try: set start values for some parameters

                 fitFcnV2_PRLTemplate->SetParName(0,"V2");
                   fitFcnV2_PRLTemplate->SetParName(1,"F");
           
         //  fitFcnV2_PRLTemplate->SetParLimits(0,0.0001,5000000);
        //   fitFcnV2_PRLTemplate->FixParameter(0,0.002);
           fitFcnV2_PRLTemplate->SetParLimits(1,0.1,50);
                 //  fitFcnV2_PRLTemplate->FixParameter(1,2.5);
                 

                TFitResultPtr res = Yields_Central_1_PRLTemplate->Fit("fitFcnV2_PRLTemplate","USBMERI+","ep");
           
               double chi2, edm, errdef;
               int nvpar, nparx;
                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->GetStats(chi2,edm,errdef,nvpar,nparx);
                fitFcnV2_PRLTemplate->SetChisquare(chi2);
                int ndf = npfits-nvpar;
                fitFcnV2_PRLTemplate->SetNDF(ndf);
           
               Double_t par[2];
               fitFcnV2_PRLTemplate->GetParameters(par);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetParameters(par);
           fitFcnV2_PRLTemplate_PeriphAndG->SetParameters(par);
                // improve the pictu
              //   std::cout << "integral error: " << integralerror << std::endl;
                fitFcnV2_PRLTemplate->Draw("same");
           fitFcnV2_PRLTemplate_RidgeAndZero->Draw("same");
           fitFcnV2_PRLTemplate_PeriphAndG->Draw("same");
                // draw the legend
                TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                legend->SetTextFont(72);
                legend->SetTextSize(0.04);
                  Char_t message[80];
                  sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate->GetChisquare(),fitFcnV2_PRLTemplate->GetNDF());
                  legend->AddEntry(fitFcnV2_PRLTemplate,message);
           if(res->CovMatrixStatus() == 3){
                      sprintf(message,"The fit is a success");
                  }
                  else{
                      sprintf(message,"The fit is a failure");
                  }
                  legend->AddEntry(fitFcnV2_PRLTemplate,message);
                legend->AddEntry(Yields_Central_1_PRLTemplate,"Data","lpe");
                legend->Draw();

       }
    
    // Fit PRL Template + ZYAM Periph
       
       TCanvas* cPRLTemplate_PeriphZYAM = new TCanvas;
          cPRLTemplate_PeriphZYAM->SetTitle("PRL Template fit Periph ZYAM");
          cPRLTemplate_PeriphZYAM->Divide(1,1);
          cPRLTemplate_PeriphZYAM->cd(1);
          Yields_Central_1->DrawCopy();
//          cPRLTemplate_PeriphZYAM->cd(2);
//          Baseline_Central_1->DrawCopy();
//          cPRLTemplate_PeriphZYAM->cd(3);
//          Yields_Central_1_MinusBaseline->DrawCopy();
          
          
          {
              cPRLTemplate_PeriphZYAM->cd(1);
                   TF1 *fitFcnV2_PRLTemplate_PeriphZYAM = new TF1("fitFcnV2_PRLTemplate_PeriphZYAM",PRLTemplate_PeriphZYAM,-TMath::Pi()/2,1.5*TMath::Pi(),2);
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetNpx(500);
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetLineWidth(4);
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetLineColor(kBlack);
                   // first try without starting values for the parameters
                   // This defaults to 1 for each param.
                   // this results in an ok fit for the polynomial function
                   // however the non-linear part (lorenzian) does not
                   // respond well.
                    Double_t params[2] = {1,1};
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetParameters(params);
                    TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->SetMaxIterations(10000);
                    TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->SetPrecision();
                    TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->SetFCN(ChisquarePRLTemplate_PeriphZYAM);
                  gStyle->SetOptFit(1011);
                  //  histo->Fit("fitFcn","0");
                   // second try: set start values for some parameters

                    fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(0,"V2");
                      fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(1,"F");
              
            //  fitFcnV2_PRLTemplate->SetParLimits(0,0.0001,5000000);
           //   fitFcnV2_PRLTemplate->FixParameter(0,0.002);
              fitFcnV2_PRLTemplate_PeriphZYAM->SetParLimits(1,0.1,5000);
                    //  fitFcnV2_PRLTemplate->FixParameter(1,2.5);
                    

                   TFitResultPtr res = Yields_Central_1_PRLTemplate_PeriphZYAM->Fit("fitFcnV2_PRLTemplate_PeriphZYAM","USBMERI+","ep");
              
              double chi2, edm, errdef;
              int nvpar, nparx;
               TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->GetStats(chi2,edm,errdef,nvpar,nparx);
               fitFcnV2_PRLTemplate_PeriphZYAM->SetChisquare(chi2);
               int ndf = npfits-nvpar;
               fitFcnV2_PRLTemplate_PeriphZYAM->SetNDF(ndf);
              
                  Double_t par[2];
                  fitFcnV2_PRLTemplate_PeriphZYAM->GetParameters(par);
                   // improve the pictu
                 //   std::cout << "integral error: " << integralerror << std::endl;
                   fitFcnV2_PRLTemplate_PeriphZYAM->Draw("same");
                   // draw the legend
                   TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                   legend->SetTextFont(72);
                   legend->SetTextSize(0.04);
                     Char_t message[80];
                     sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate_PeriphZYAM->GetChisquare(),fitFcnV2_PRLTemplate_PeriphZYAM->GetNDF());
                     legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
              if(res->CovMatrixStatus() == 3){
                         sprintf(message,"The fit is a success");
                     }
                     else{
                         sprintf(message,"The fit is a failure");
                     }
                     legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                   legend->AddEntry(Yields_Central_1_PRLTemplate_PeriphZYAM,"Data","lpe");
                   legend->Draw();
              

          }
    
    
    
    
    
    
    
    
    
    
       //DimuTkl Ext1 Pt Binned
    
    
    if(PtBinned){
        
        for(int ptbin=0;ptbin<NbPtBins;ptbin++){
            
            baseline_centralPtBinned[ptbin] = 9999.;
    
        int niterations = 1;
        cinvmassPtBinned->cd();
        sprintf(histoname,Form("bin%d_1",ptbin+1));
        res = FittingAllInvMassBin(histoname, cinvmassPtBinned, NbPtBins + ptbin);

       double parerr[12];
        TCanvas* c10_fitPtBinned = new TCanvas;
           c10_fitPtBinned->Divide(NbinsDeltaPhi/2,4);
        TRandom rand = 0;
        for(int j=0; j<12; j++){
            parerr[j] = res->ParError(j);
        }
        for(int j=0; j <16; j++){
            par[j] = res->Parameter(j);
            if(j>11){
                par[j] = 3;
            }
        }
            //PLOT ICI

        double par12[niterations];
        double parerr12[niterations];

        
        for(int i=0; i<NbinsDeltaPhi; i++){ //NbPhiBins
            for(int t=0; t<niterations; t++){
                TVirtualPad* subpad = c10_fitPtBinned->cd(i+1);
                c10_fitPtBinned->cd(i+1);
                if(t==niterations-1){
                    subpad->Clear();
                }
            TF1 *fitY_1Central = new TF1("fitY_1Central",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
              fitY_1Central->SetNpx(500);
              fitY_1Central->SetLineWidth(4);
              fitY_1Central->SetLineColor(kMagenta);
              fitY_1Central->SetParameters(par);
               TVirtualFitter::Fitter(YieldWrtMass_CentralPtBinned[i][ptbin])->SetMaxIterations(10000);
               TVirtualFitter::Fitter(YieldWrtMass_CentralPtBinned[i][ptbin])->SetPrecision();
            //  histo->Fit("fitFcn","0");
              // second try: set start values for some parameters
            for(int k=0; k<=11; k++){

                    fitY_1Central->FixParameter(k,par[k]);
                
            }
                fitY_1Central->SetParLimits(12,0,10000);

               fitY_1Central->SetParName(0,"Norm_{JPsi}");
               fitY_1Central->SetParName(1,"M_{JPsi}");
               fitY_1Central->SetParName(2,"Sigma_{JPsi}");
               fitY_1Central->SetParName(3,"a_{1}");
               fitY_1Central->SetParName(4,"n_{1}");
               fitY_1Central->SetParName(5,"a_{2}");
               fitY_1Central->SetParName(6,"n_{2}");
               fitY_1Central->SetParName(7,"Norm_{Psi2S}");
               fitY_1Central->SetParName(8,"Norm_{TailLowM}");
               fitY_1Central->SetParName(9,"Exp_{TailLowM}");
               fitY_1Central->SetParName(10,"Norm_{TailHighM}");
               fitY_1Central->SetParName(11,"Exp_{TailHighM}");
                fitY_1Central->SetParName(12,"Y_1 JPsi");
                fitY_1Central->SetParName(13,"Y_1 Bkg M2");
            fitY_1Central->SetParName(14,"Y_1 Bkg M1");
            fitY_1Central->SetParName(15,"Y_1 Nkg M0");

              rescent = YieldWrtMass_CentralPtBinned[i][ptbin]->Fit("fitY_1Central","SBMERIQ+","ep");
              Double_t param[16];
                Double_t paramerrs[16];
              fitY_1Central->GetParameters(param);
                for(int i=0; i<16; i++){
                    paramerrs[i] = fitY_1Central->GetParError(i);
                }
                
                if(rescent->CovMatrixStatus() == 3){
                par12[t] = param[12];
                parerr12[t] = paramerrs[12];
                }
                else{
                    par12[t]=parerr12[t]=0;
                }
                
                std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
                std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

            TF1 *backFcnY_1Central = new TF1("backFcnY_1Central",BackFcnV2Poly,2.1,5.1,16);
            backFcnY_1Central->SetLineColor(kRed);
            TF1 *signalFcnJPsiY_1Central = new TF1("signalFcnJPsiY_1Central",SignalFcnJPsiV2,2.1,5.1,16);
              // writes the fit results into the par array
            
                if(t==niterations-1){
                    YieldWrtMass_CentralPtBinned[i][ptbin]->GetListOfFunctions()->Clear();
                    YieldWrtMass_CentralPtBinned[i][ptbin]->GetListOfFunctions()->Add(fitY_1Central);
                }
            signalFcnJPsiY_1Central->SetLineColor(kBlue);
            signalFcnJPsiY_1Central->SetNpx(500);
            backFcnY_1Central->SetParameters(param);
            signalFcnJPsiY_1Central->SetParameters(param);
          //  signalFcnJPsiY_1Central->Draw("same");
            backFcnY_1Central->Draw("same");
                
                if(param[12]<baseline_centralPtBinned[ptbin]){ //En central baseline = minimum
                    baseline_centralPtBinned[ptbin] = param[12];
                    errbaseline_centralPtBinned[ptbin] = fitY_1Central->GetParError(12);
                }

                    Yields_Central_1PtBinned[ptbin]->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[12]); //AEUGG
                    Yields_Central_1PtBinned[ptbin]->SetBinError(i+1,fitY_1Central->GetParError(12));
                
                
              // writes the fit results into the par array
               gStyle->SetOptFit(1011);
              // draw the legend
              TLegend *legend=new TLegend(0.12,0.75,0.60,0.90);
              legend->SetTextFont(61);
              legend->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitY_1Central->GetChisquare(),fitY_1Central->GetNDF());
             legend->AddEntry(fitY_1Central,message);
            if(rescent->CovMatrixStatus() == 3){
                       sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                   }
                   legend->AddEntry(fitY_1Central,message);
         //   legend->AddEntry(signalFcnJPsiY_1Central,"JPsi signal");
            legend->AddEntry(backFcnY_1Central,"Background");
              legend->AddEntry(YieldWrtMass_CentralPtBinned[i][ptbin],"Data","lpe");
              legend->Draw();
            }
            
        }
            
            for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
                Baseline_Central_1PtBinned[ptbin]->SetBinContent(phi_idx+1, baseline_centralPtBinned[ptbin]);
                Baseline_Central_1PtBinned[ptbin]->SetBinError(phi_idx+1, errbaseline_centralPtBinned[ptbin]);
            }
            
            Yields_Central_1_MinusBaselinePtBinned[ptbin]->Add(Yields_Central_1PtBinned[ptbin],Baseline_Central_1PtBinned[ptbin],1,-1);
        
            
        cinvmassPtBinned->cd(); //aeugg
        sprintf(histoname,Form("bin%d_2",ptbin+1));
        res = FittingAllInvMassBin(histoname, cinvmassPtBinned, NbPtBins+NbPtBins+ptbin);
        
           for(int j=0; j<12; j++){
               parerr[j] = res->ParError(j);
           }
           for(int j=0; j <16; j++){
               par[j] = res->Parameter(j);
               if(j>11){
                   par[j] = 3;
               }
           }
        
           
           for(int i=0; i<NbinsDeltaPhi; i++){ //NbPhiBins
               for(int t=0; t<niterations; t++){
                   TVirtualPad* subpad = c10_fitPtBinned->cd(NbinsDeltaPhi+i+1);
                   c10_fitPtBinned->cd(NbinsDeltaPhi+i+1);
                   if(t==niterations-1){
                       subpad->Clear();
                   }
               TF1 *fitY_1Periph = new TF1("fitY_1Periph",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
                 fitY_1Periph->SetNpx(500);
                 fitY_1Periph->SetLineWidth(4);
                 fitY_1Periph->SetLineColor(kMagenta);
                 fitY_1Periph->SetParameters(par);
                  TVirtualFitter::Fitter(YieldWrtMass_PeriphPtBinned[i][ptbin])->SetMaxIterations(10000);
                  TVirtualFitter::Fitter(YieldWrtMass_PeriphPtBinned[i][ptbin])->SetPrecision();
               //  histo->Fit("fitFcn","0");
                 // second try: set start values for some parameters
                for(int k=0; k<=11; k++){
                        fitY_1Periph->FixParameter(k,par[k]);
                }
                  fitY_1Periph->SetParLimits(12,0,10000);
                  
                  fitY_1Periph->SetParName(0,"Norm_{JPsi}");
                  fitY_1Periph->SetParName(1,"M_{JPsi}");
                  fitY_1Periph->SetParName(2,"Sigma_{JPsi}");
                  fitY_1Periph->SetParName(3,"a_{1}");
                  fitY_1Periph->SetParName(4,"n_{1}");
                  fitY_1Periph->SetParName(5,"a_{2}");
                  fitY_1Periph->SetParName(6,"n_{2}");
                  fitY_1Periph->SetParName(7,"Norm_{Psi2S}");
                  fitY_1Periph->SetParName(8,"Norm_{TailLowM}");
                  fitY_1Periph->SetParName(9,"Exp_{TailLowM}");
                  fitY_1Periph->SetParName(10,"Norm_{TailHighM}");
                  fitY_1Periph->SetParName(11,"Exp_{TailHighM}");
                   fitY_1Periph->SetParName(12,"Y_1 JPsi");
                   fitY_1Periph->SetParName(13,"Y_1 Bkg M2");
               fitY_1Periph->SetParName(14,"Y_1 Bkg M1");
               fitY_1Periph->SetParName(15,"Y_1 Nkg M0");
                  
                 resperiph = YieldWrtMass_PeriphPtBinned[i][ptbin]->Fit("fitY_1Periph","SBMERIQ+","ep");
                   Double_t param[16];
                   Double_t paramerrs[16];
                 fitY_1Periph->GetParameters(param);
                   
                   for(int i=0; i<16; i++){
                                   paramerrs[i] = fitY_1Periph->GetParError(i);
                               }
                               
                               if(resperiph->CovMatrixStatus() == 3){
                               par12[t] = param[12];
                               parerr12[t] = paramerrs[12];
                               }
                               else{
                                   par12[t]=parerr12[t]=0;
                               }
                   std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
                    std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

               TF1 *backFcnY_1Periph = new TF1("backFcnY_1Periph",BackFcnV2Poly,2.1,5.1,16);
               backFcnY_1Periph->SetLineColor(kRed);
               TF1 *signalFcnJPsiY_1Periph = new TF1("signalFcnJPsiY_1Periph",SignalFcnJPsiV2,2.1,5.1,16);
                 // writes the fit results into the par array
                   
                   if(t==niterations-1){
                       YieldWrtMass_PeriphPtBinned[i][ptbin]->GetListOfFunctions()->Clear();
                       YieldWrtMass_PeriphPtBinned[i][ptbin]->GetListOfFunctions()->Add(fitY_1Periph);
                   }
               
               signalFcnJPsiY_1Periph->SetLineColor(kBlue);
               signalFcnJPsiY_1Periph->SetNpx(500);
               backFcnY_1Periph->SetParameters(param);
               signalFcnJPsiY_1Periph->SetParameters(param);
            //   signalFcnJPsiY_1Periph->Draw("same");
               backFcnY_1Periph->Draw("same");
                   
                      if(i==2){
                           baseline_periphPtBinned[ptbin] = param[12];
                           errbaseline_periphPtBinned[ptbin] = fitY_1Periph->GetParError(12);
                       }
                       if(i==3){
                           baseline_periphPtBinned[ptbin] += param[12];
                           baseline_periphPtBinned[ptbin] /= 2;
                           errbaseline_periphPtBinned[ptbin] = sqrt(pow(errbaseline_periphPtBinned[ptbin],2)+pow(fitY_1Periph->GetParError(12),2));
                       }
                   
                     Yields_Periph_1PtBinned[ptbin]->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[12]);
                     Yields_Periph_1PtBinned[ptbin]->SetBinError(i+1,fitY_1Periph->GetParError(12));
                   
                   
                 // writes the fit results into the par array
                  gStyle->SetOptFit(1011);
                 // draw the legend
                 TLegend *legend=new TLegend(0.12,0.75,0.60,0.90);
                 legend->SetTextFont(61);
                 legend->SetTextSize(0.03);
               Char_t message[80];
               sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitY_1Periph->GetChisquare(),fitY_1Periph->GetNDF());
                legend->AddEntry(fitY_1Periph,message);
               if(resperiph->CovMatrixStatus() == 3){
                          sprintf(message,"The fit is a success");
                      }
                      else{
                          sprintf(message,"The fit is a failure");
                      }
                      legend->AddEntry(fitY_1Periph,message);
          //     legend->AddEntry(signalFcnJPsiY_1Periph,"JPsi signal");
               legend->AddEntry(backFcnY_1Periph,"Background");
                 legend->AddEntry(YieldWrtMass_Periph[i],"Data","lpe");
                 legend->Draw();
               }
               
           }
            
            for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
                Baseline_Periph_1PtBinned[ptbin]->SetBinContent(phi_idx+1, baseline_periphPtBinned[ptbin]);
                Baseline_Periph_1PtBinned[ptbin]->SetBinError(phi_idx+1, errbaseline_periphPtBinned[ptbin]);
            }
            
            Yields_Periph_1_MinusBaselinePtBinned[ptbin]->Add(Yields_Periph_1PtBinned[ptbin],Baseline_Periph_1PtBinned[ptbin],1,-1);
            
            //Fill Baseline Periph and Y-B
            char hname[100];
            TCanvas* cYBcentralPtBinned = new TCanvas;
            sprintf(hname, "Central Yield-Baseline - PtBinned %d", ptbin);
            cYBcentralPtBinned->SetTitle(hname);
            cYBcentralPtBinned->Divide(1,3);
            cYBcentralPtBinned->cd(1);
            Yields_Central_1PtBinned[ptbin]->DrawCopy();
            cYBcentralPtBinned->cd(2);
            Baseline_Central_1PtBinned[ptbin]->DrawCopy();
            cYBcentralPtBinned->cd(3);
            Yields_Central_1_MinusBaselinePtBinned[ptbin]->DrawCopy();
            
            TCanvas* cYBperiphPtBinned = new TCanvas;
            sprintf(hname, "Periph Yield-Baseline - PtBinned %d", ptbin);
            cYBperiphPtBinned->SetTitle(hname);
            cYBperiphPtBinned->Divide(1,3);
            cYBperiphPtBinned->cd(1);
            Yields_Periph_1PtBinned[ptbin]->DrawCopy();
            cYBperiphPtBinned->cd(2);
            Baseline_Periph_1PtBinned[ptbin]->DrawCopy();
            cYBperiphPtBinned->cd(3);
            Yields_Periph_1_MinusBaselinePtBinned[ptbin]->DrawCopy();
            
            sprintf(hname, "Yields_Central_1_CvetanPtBinned - PtBinned %d", ptbin);
            Yields_Central_1_CvetanPtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_CvetanMePtBinned - PtBinned %d", ptbin);
            Yields_Central_1_CvetanMePtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_ZYAMPtBinned - PtBinned %d", ptbin);
            Yields_Central_1_ZYAMPtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_PRLTemplatePtBinned - PtBinned %d", ptbin);
            Yields_Central_1_PRLTemplatePtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned - PtBinned %d", ptbin);
            Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            
            {
                TCanvas* c17PtBinned = new TCanvas;
                c17PtBinned->Divide(1,3);
                //Creer canvas pour imprimer les plots Periph yield et Central yield wrt phi et leur difference
                c17PtBinned->cd(1);
                Yields_Central_1PtBinned[ptbin]->Draw();
                c17PtBinned->cd(2);
                Yields_Periph_1PtBinned[ptbin]->Draw();
                c17PtBinned->cd(3);
                Yields_Difference_1PtBinned[ptbin]->Add(Yields_Central_1PtBinned[ptbin],Yields_Periph_1PtBinned[ptbin],1,-1);
                Yields_Difference_1PtBinned[ptbin]->Draw();
            
                     TF1 *fitFcnV2_2 = new TF1("fitFcnV2_2",FourierV2,-TMath::Pi()/2,1.5*TMath::Pi(),3);
                     fitFcnV2_2->SetNpx(500);
                     fitFcnV2_2->SetLineWidth(4);
                     fitFcnV2_2->SetLineColor(kMagenta);
                     // first try without starting values for the parameters
                     // This defaults to 1 for each param.
                     // this results in an ok fit for the polynomial function
                     // however the non-linear part (lorenzian) does not
                     // respond well.
                      Double_t params[3] = {1,1,1};
                     fitFcnV2_2->SetParameters(params);
                      TVirtualFitter::Fitter(Yields_Difference_1PtBinned[ptbin])->SetMaxIterations(10000);
                      TVirtualFitter::Fitter(Yields_Difference_1PtBinned[ptbin])->SetPrecision();
                   //  histo->Fit("fitFcn","0");
                     // second try: set start values for some parameters
                      
                      fitFcnV2_2->SetParName(0,"a0");
                      fitFcnV2_2->SetParName(1,"a1");
                      fitFcnV2_2->SetParName(2,"a2");
                      
                     TFitResultPtr res = Yields_Difference_1PtBinned[ptbin]->Fit("fitFcnV2_2","SBMERI+","ep");
                    Double_t par[3];
                    fitFcnV2_2->GetParameters(par);
                     // improve the pictu
                   //   std::cout << "integral error: " << integralerror << std::endl;
                     fitFcnV2_2->Draw("same");
                     // draw the legend
                     TLegend *legend=new TLegend(0.15,0.70,0.4,0.90);
                     legend->SetTextFont(61);
                     legend->SetTextSize(0.03);
                       Char_t message[80];
                       sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
                       legend->AddEntry(fitFcnV2_2,message);
                sprintf(message,"V2_1 JPsi-Tkl: %.4f / (%.4f + %.4f) = %.4f +- %.4f",par[2],par[0], baseline_periphPtBinned[ptbin],par[2]/(par[0] + baseline_periphPtBinned[ptbin]),(par[2]/(par[0] + baseline_periphPtBinned[ptbin])*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(fitFcnV2_2->GetParError(0)/par[0],2)+pow(errbaseline_periphPtBinned[ptbin]/baseline_periphPtBinned[ptbin],2))));
                legend->AddEntry(fitFcnV2_2,message);
                if(res->CovMatrixStatus() == 3){
                           sprintf(message,"The fit is a success");
                       }
                       else{
                           sprintf(message,"The fit is a failure");
                       }
                       legend->AddEntry(fitFcnV2_2,message);
                     legend->AddEntry(Yields_Difference_1,"Data","lpe");
                     legend->Draw();
                
                V2_Ext1[ptbin] = par[2]/(par[0] + baseline_periphPtBinned[ptbin]);
                errV2_Ext1[ptbin] = (par[2]/(par[0] + baseline_periphPtBinned[ptbin])*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(fitFcnV2_2->GetParError(0)/par[0],2)+pow(errbaseline_periphPtBinned[ptbin]/baseline_periphPtBinned[ptbin],2)));
                   
                
            }
            
            
             // Fit Cvetan PtBinned
                
                TCanvas* cCvetanPtBinned = new TCanvas;
            sprintf(hname, "Cvetan fit - PtBinned %d", ptbin);
                   cCvetanPtBinned->SetTitle(hname);
                   cCvetanPtBinned->Divide(1,1);
                   cCvetanPtBinned->cd(1);
                   Yields_Central_1->DrawCopy();
            //       cCvetan->cd(2);
            //       Baseline_Central_1->DrawCopy();
            //       cCvetan->cd(3);
            //       Yields_Central_1_MinusBaseline->DrawCopy();
            
            
                
               {
                    cCvetanPtBinned->cd(1);
                         TF1 *fitFcnV2_Cvetan = new TF1("fitFcnV2_Cvetan",CvetanFPtBinned,-TMath::Pi()/2,1.5*TMath::Pi(),5);
                         fitFcnV2_Cvetan->SetNpx(500);
                         fitFcnV2_Cvetan->SetLineWidth(4);
                         fitFcnV2_Cvetan->SetLineColor(kBlue);
                         // first try without starting values for the parameters
                         // This defaults to 1 for each param.
                         // this results in an ok fit for the polynomial function
                         // however the non-linear part (lorenzian) does not
                         // respond well.
                          Double_t params[5] = {static_cast<Double_t>(ptbin),1,0,0.01,1};
                         fitFcnV2_Cvetan->SetParameters(params);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanPtBinned[ptbin])->SetMaxIterations(10000);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanPtBinned[ptbin])->SetPrecision();
                        TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetFCN(ChisquareCvetanFPtBinned);
                        gStyle->SetOptFit(1011);
                       //  histo->Fit("fitFcn","0");
                         // second try: set start values for some parameters

                        fitFcnV2_Cvetan->SetParName(0,"ptbin");
                          fitFcnV2_Cvetan->SetParName(1,"V0");
                            fitFcnV2_Cvetan->SetParName(2,"V1");
                            fitFcnV2_Cvetan->SetParName(3,"V2");
                          fitFcnV2_Cvetan->SetParName(4,"F");
                    
                    fitFcnV2_Cvetan->SetParLimits(4,0.1,50);
                    
                    fitFcnV2_Cvetan->FixParameter(0,ptbin);
                    fitFcnV2_Cvetan->FixParameter(1,1);
                    fitFcnV2_Cvetan->FixParameter(2,0);
                  //  fitFcnV2_Cvetan->FixParameter(4,1);
                          

                         TFitResultPtr res = Yields_Central_1_CvetanPtBinned[ptbin]->Fit("fitFcnV2_Cvetan","USBMERI+","ep");
                        Double_t par[5];
                    double chi2, edm, errdef;
                    int nvpar, nparx;
                        TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                    
                         fitFcnV2_Cvetan->SetChisquare(chi2);
                         int ndf = npfits-nvpar;
                         fitFcnV2_Cvetan->SetNDF(ndf);
                    
                    
                        fitFcnV2_Cvetan->GetParameters(par);
                         // improve the pictu
                       //   std::cout << "integral error: " << integralerror << std::endl;
                         fitFcnV2_Cvetan->Draw("same");
                         // draw the legend
                         TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                         legend->SetTextFont(72);
                         legend->SetTextSize(0.04);
                           Char_t message[80];
                           sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_Cvetan->GetChisquare(),fitFcnV2_Cvetan->GetNDF());
                           legend->AddEntry(fitFcnV2_Cvetan,message);
                    if(res->CovMatrixStatus() == 3){
                               sprintf(message,"The fit is a success");
                           }
                           else{
                               sprintf(message,"The fit is a failure");
                           }
                           legend->AddEntry(fitFcnV2_Cvetan,message);
                         legend->AddEntry(Yields_Central_1_CvetanPtBinned[ptbin],"Data","lpe");
                         legend->Draw();

                    V2_CvetanQuentin[ptbin] = par[3];
                    errV2_CvetanQuentin[ptbin] = fitFcnV2_Cvetan->GetParError(3);
                    F_CvetanQuentin[ptbin] = par[4];
                    errF_CvetanQuentin[ptbin] = fitFcnV2_Cvetan->GetParError(4);
                    
                }
            
            
            // Fit Cvetan putting my constraints (F=1, V1 free, V0 not 1) - PtBinned
                
                TCanvas* cCvetanMePtBinned = new TCanvas;
            sprintf(hname, "Cvetan fit - My constraints F=1, v1 free, v0 free - PtBinned %d", ptbin);
                   cCvetanMePtBinned->SetTitle(hname);
                   cCvetanMePtBinned->Divide(1,1);
                   cCvetanMePtBinned->cd(1);
                   Yields_Central_1->DrawCopy();
            //       cCvetanMe->cd(2);
            //       Baseline_Central_1->DrawCopy();
            //       cCvetanMe->cd(3);
            //       Yields_Central_1_MinusBaseline->DrawCopy();
                
                
                {
                    cCvetanMePtBinned->cd(1);
                         TF1 *fitFcnV2_CvetanMe = new TF1("fitFcnV2_CvetanMe",CvetanFPtBinned,-TMath::Pi()/2,1.5*TMath::Pi(),5);
                         fitFcnV2_CvetanMe->SetNpx(500);
                         fitFcnV2_CvetanMe->SetLineWidth(4);
                         fitFcnV2_CvetanMe->SetLineColor(kBlue);
                         // first try without starting values for the parameters
                         // This defaults to 1 for each param.
                         // this results in an ok fit for the polynomial function
                         // however the non-linear part (lorenzian) does not
                         // respond well.
                          Double_t params[5] = {static_cast<Double_t>(ptbin),1,0,0.01,1};
                         fitFcnV2_CvetanMe->SetParameters(params);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetMaxIterations(10000);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetPrecision();
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetFCN(ChisquareCvetanFPtBinned);
                        gStyle->SetOptFit(1011);
                       //  histo->Fit("fitFcn","0");
                         // second try: set start values for some parameters

                        fitFcnV2_CvetanMe->SetParName(0,"ptbin");
                          fitFcnV2_CvetanMe->SetParName(1,"V0");
                            fitFcnV2_CvetanMe->SetParName(2,"V1");
                            fitFcnV2_CvetanMe->SetParName(3,"V2");
                          fitFcnV2_CvetanMe->SetParName(4,"F");
                    
                    fitFcnV2_CvetanMe->SetParLimits(4,0.1,50);
                    
                    fitFcnV2_CvetanMe->FixParameter(0,ptbin);
                 //   fitFcnV2_CvetanMe->FixParameter(1,1);
                   // fitFcnV2_CvetanMe->FixParameter(2,0);
                    fitFcnV2_CvetanMe->FixParameter(4,1);
                          

                         TFitResultPtr res = Yields_Central_1_CvetanMePtBinned[ptbin]->Fit("fitFcnV2_CvetanMe","USBMERI+","ep");
                    
                    double chi2, edm, errdef;
                    int nvpar, nparx;
                     TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                     fitFcnV2_CvetanMe->SetChisquare(chi2);
                     int ndf = npfits-nvpar;
                     fitFcnV2_CvetanMe->SetNDF(ndf);
                    
                        Double_t par[5];
                        fitFcnV2_CvetanMe->GetParameters(par);
                         // improve the pictu
                       //   std::cout << "integral error: " << integralerror << std::endl;
                         fitFcnV2_CvetanMe->Draw("same");
                         // draw the legend
                         TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                         legend->SetTextFont(72);
                         legend->SetTextSize(0.04);
                           Char_t message[80];
                           sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_CvetanMe->GetChisquare(),fitFcnV2_CvetanMe->GetNDF());
                           legend->AddEntry(fitFcnV2_CvetanMe,message);
                    if(res->CovMatrixStatus() == 3){
                               sprintf(message,"The fit is a success");
                           }
                           else{
                               sprintf(message,"The fit is a failure");
                           }
                           legend->AddEntry(fitFcnV2_CvetanMe,message);
                         legend->AddEntry(Yields_Central_1_CvetanMePtBinned[ptbin],"Data","lpe");
                         legend->Draw();
                    
                    V2_CvetanQuentinMe[ptbin] = par[3]/par[1];
                    errV2_CvetanQuentinMe[ptbin] = (par[3]/(par[1])*sqrt(pow(fitFcnV2_CvetanMe->GetParError(3)/par[3],2)+pow(fitFcnV2_CvetanMe->GetParError(1)/par[1],2)));

                }
            
            
            // Fit ZYAM Yc = Bc + (Yp-Bp) + a0  + 2a1cos + 2a2cos  -  PtBinned
            
            
             TCanvas* cZYAMPtBinned = new TCanvas;
            sprintf(hname, "ZYAM fit - PtBinned %d", ptbin);
                cZYAMPtBinned->SetTitle(hname);
                cZYAMPtBinned->Divide(1,1);
                cZYAMPtBinned->cd(1);
                Yields_Central_1PtBinned[ptbin]->DrawCopy();
            //    cZYAM->cd(2);
            //    Baseline_Central_1->DrawCopy();
            //    cZYAM->cd(3);
            //    Yields_Central_1_MinusBaseline->DrawCopy();
            
            
            {
                cZYAMPtBinned->cd(1);
                     TF1 *fitFcnV2_ZYAM = new TF1("fitFcnV2_ZYAM",ZYAMPtBinned,-TMath::Pi()/2,1.5*TMath::Pi(),4);
                     fitFcnV2_ZYAM->SetNpx(500);
                     fitFcnV2_ZYAM->SetLineWidth(4);
                     fitFcnV2_ZYAM->SetLineColor(kGreen);
                     // first try without starting values for the parameters
                     // This defaults to 1 for each param.
                     // this results in an ok fit for the polynomial function
                     // however the non-linear part (lorenzian) does not
                     // respond well.
                      Double_t params[4] = {static_cast<Double_t>(ptbin),1,1,1};
                     fitFcnV2_ZYAM->SetParameters(params);
                      TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->SetMaxIterations(10000);
                      TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->SetPrecision();
                    TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->SetFCN(ChisquareZYAMPtBinned);
                    gStyle->SetOptFit(1011);
                    //  histo->Fit("fitFcn","0");
                     // second try: set start values for some parameters

                    fitFcnV2_ZYAM->SetParName(0,"ptbin");
                      fitFcnV2_ZYAM->SetParName(1,"a0");
                        fitFcnV2_ZYAM->SetParName(2,"a1");
                        fitFcnV2_ZYAM->SetParName(3,"a2");
                    
                       fitFcnV2_ZYAM->FixParameter(0,ptbin);

                     TFitResultPtr res = Yields_Central_1_ZYAMPtBinned[ptbin]->Fit("fitFcnV2_ZYAM","USBMERI+","ep");
                
                double chi2, edm, errdef;
                int nvpar, nparx;
                 TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                 fitFcnV2_ZYAM->SetChisquare(chi2);
                 int ndf = npfits-nvpar;
                 fitFcnV2_ZYAM->SetNDF(ndf);
                
                    Double_t par[4];
                    fitFcnV2_ZYAM->GetParameters(par);
                     // improve the pictu
                   //   std::cout << "integral error: " << integralerror << std::endl;
                     fitFcnV2_ZYAM->Draw("same");
                     // draw the legend
                     TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                     legend->SetTextFont(72);
                     legend->SetTextSize(0.04);
                       Char_t message[80];
                       sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_ZYAM->GetChisquare(),fitFcnV2_ZYAM->GetNDF());
                       legend->AddEntry(fitFcnV2_ZYAM,message);
                if(res->CovMatrixStatus() == 3){
                           sprintf(message,"The fit is a success");
                       }
                       else{
                           sprintf(message,"The fit is a failure");
                       }
                       legend->AddEntry(fitFcnV2_ZYAM,message);
                     legend->AddEntry(Yields_Central_1_ZYAMPtBinned[ptbin],"Data","lpe");
                     legend->Draw();
                
                
                V2_ZYAM[ptbin] = par[3]/(par[1]+baseline_centralPtBinned[ptbin]);
                errV2_ZYAM[ptbin] = (par[3]/(par[1] + baseline_centralPtBinned[ptbin])*sqrt(pow(fitFcnV2_ZYAM->GetParError(3)/par[3],2)+pow(fitFcnV2_ZYAM->GetParError(1)/par[1],2)+pow(errbaseline_centralPtBinned[ptbin]/baseline_centralPtBinned[ptbin],2)));

            }
            
            
             // Fit PRL Template - PtBinned
                
                TCanvas* cPRLTemplatePtBinned = new TCanvas;
                sprintf(hname, "PRL Template fit - PtBinned %d", ptbin);
                   cPRLTemplatePtBinned->SetTitle(hname);
                   cPRLTemplatePtBinned->Divide(1,1);
                   cPRLTemplatePtBinned->cd(1);
                   Yields_Central_1PtBinned[ptbin]->DrawCopy();
            //       cPRLTemplate->cd(2);
            //       Baseline_Central_1->DrawCopy();
            //       cPRLTemplate->cd(3);
            //       Yields_Central_1_MinusBaseline->DrawCopy();
                   
                   
                   {
                       cPRLTemplatePtBinned->cd(1);
                            TF1 *fitFcnV2_PRLTemplate = new TF1("fitFcnV2_PRLTemplate",PRLTemplatePtBinned,-TMath::Pi()/2,1.5*TMath::Pi(),3);
                            fitFcnV2_PRLTemplate->SetNpx(500);
                            fitFcnV2_PRLTemplate->SetLineWidth(4);
                            fitFcnV2_PRLTemplate->SetLineColor(kRed);
                       TF1 *fitFcnV2_PRLTemplate_RidgeAndZero = new TF1("fitFcnV2_PRLTemplate_RidgeAndZero",PRLTemplate_RidgeAndZeroPtBinned,-TMath::Pi()/2,1.5*TMath::Pi(),3);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetNpx(500);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetLineWidth(2);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetLineStyle(9);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetLineColor(kRed);
                       TF1 *fitFcnV2_PRLTemplate_PeriphAndG = new TF1("fitFcnV2_PRLTemplate_PeriphAndG",PRLTemplate_PeriphAndGPtBinned,-TMath::Pi()/2,1.5*TMath::Pi(),3);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetNpx(500);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetLineWidth(2);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetLineStyle(9);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetLineColor(kRed);
                            // first try without starting values for the parameters
                            // This defaults to 1 for each param.
                            // this results in an ok fit for the polynomial function
                            // however the non-linear part (lorenzian) does not
                            // respond well.
                             Double_t params[3] = {static_cast<Double_t>(ptbin),1,1};
                            fitFcnV2_PRLTemplate->SetParameters(params);
                             TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->SetMaxIterations(10000);
                             TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->SetPrecision();
                            TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->SetFCN(ChisquarePRLTemplatePtBinned);
                           gStyle->SetOptFit(1011);
                           //  histo->Fit("fitFcn","0");
                            // second try: set start values for some parameters

                            fitFcnV2_PRLTemplate->SetParName(0,"ptbin");
                             fitFcnV2_PRLTemplate->SetParName(1,"V2");
                               fitFcnV2_PRLTemplate->SetParName(2,"F");
                       

                       fitFcnV2_PRLTemplate->FixParameter(0,ptbin);
                       fitFcnV2_PRLTemplate->SetParLimits(2,0.1,50);
                             

                            TFitResultPtr res = Yields_Central_1_PRLTemplatePtBinned[ptbin]->Fit("fitFcnV2_PRLTemplate","USBMERI+","ep");
                       
                       double chi2, edm, errdef;
                       int nvpar, nparx;
                        TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                        fitFcnV2_PRLTemplate->SetChisquare(chi2);
                        int ndf = npfits-nvpar;
                        fitFcnV2_PRLTemplate->SetNDF(ndf);
                       
                           Double_t par[3];
                           fitFcnV2_PRLTemplate->GetParameters(par);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetParameters(par);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetParameters(par);
                            // improve the pictu
                          //   std::cout << "integral error: " << integralerror << std::endl;
                            fitFcnV2_PRLTemplate->Draw("same");
                       fitFcnV2_PRLTemplate_RidgeAndZero->Draw("same");
                       fitFcnV2_PRLTemplate_PeriphAndG->Draw("same");
                            // draw the legend
                            TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                            legend->SetTextFont(72);
                            legend->SetTextSize(0.04);
                              Char_t message[80];
                              sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate->GetChisquare(),fitFcnV2_PRLTemplate->GetNDF());
                              legend->AddEntry(fitFcnV2_PRLTemplate,message);
                       if(res->CovMatrixStatus() == 3){
                                  sprintf(message,"The fit is a success");
                              }
                              else{
                                  sprintf(message,"The fit is a failure");
                              }
                              legend->AddEntry(fitFcnV2_PRLTemplate,message);
                            legend->AddEntry(Yields_Central_1_PRLTemplatePtBinned[ptbin],"Data","lpe");
                            legend->Draw();
                       
                       V2_PRL[ptbin] = par[1];
                       errV2_PRL[ptbin] = fitFcnV2_PRLTemplate->GetParError(1);
                       F_PRL[ptbin] = par[2];
                       errF_PRL[ptbin] = fitFcnV2_PRLTemplate->GetParError(2);

                   }
            
            
             // Fit PRL Template + ZYAM Periph - PtBinned
                   
                   TCanvas* cPRLTemplate_PeriphZYAMPtBinned = new TCanvas;
            sprintf(hname, "PRL Template fit Periph ZYAM - PtBinned %d", ptbin);
                      cPRLTemplate_PeriphZYAMPtBinned->SetTitle(hname);
                      cPRLTemplate_PeriphZYAMPtBinned->Divide(1,1);
                      cPRLTemplate_PeriphZYAMPtBinned->cd(1);
                      Yields_Central_1PtBinned[ptbin]->DrawCopy();
            //          cPRLTemplate_PeriphZYAM->cd(2);
            //          Baseline_Central_1->DrawCopy();
            //          cPRLTemplate_PeriphZYAM->cd(3);
            //          Yields_Central_1_MinusBaseline->DrawCopy();
                      
                      
                      {
                          cPRLTemplate_PeriphZYAMPtBinned->cd(1);
                               TF1 *fitFcnV2_PRLTemplate_PeriphZYAM = new TF1("fitFcnV2_PRLTemplate_PeriphZYAM",PRLTemplate_PeriphZYAMPtBinned,-TMath::Pi()/2,1.5*TMath::Pi(),3);
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetNpx(500);
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetLineWidth(4);
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetLineColor(kBlack);
                               // first try without starting values for the parameters
                               // This defaults to 1 for each param.
                               // this results in an ok fit for the polynomial function
                               // however the non-linear part (lorenzian) does not
                               // respond well.
                                Double_t params[3] = {static_cast<Double_t>(ptbin),1,1};
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetParameters(params);
                                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->SetMaxIterations(10000);
                                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->SetPrecision();
                                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->SetFCN(ChisquarePRLTemplate_PeriphZYAMPtBinned);
                              gStyle->SetOptFit(1011);
                              //  histo->Fit("fitFcn","0");
                               // second try: set start values for some parameters

                                fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(0,"ptbin");
                                fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(1,"V2");
                                  fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(2,"F");
                          
                            fitFcnV2_PRLTemplate_PeriphZYAM->FixParameter(0,ptbin);
                          fitFcnV2_PRLTemplate_PeriphZYAM->SetParLimits(2,0.1,5000);
                                

                               TFitResultPtr res = Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin]->Fit("fitFcnV2_PRLTemplate_PeriphZYAM","USBMERI+","ep");
                          
                          double chi2, edm, errdef;
                          int nvpar, nparx;
                           TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                           fitFcnV2_PRLTemplate_PeriphZYAM->SetChisquare(chi2);
                           int ndf = npfits-nvpar;
                           fitFcnV2_PRLTemplate_PeriphZYAM->SetNDF(ndf);
                          
                              Double_t par[3];
                              fitFcnV2_PRLTemplate_PeriphZYAM->GetParameters(par);
                               // improve the pictu
                             //   std::cout << "integral error: " << integralerror << std::endl;
                               fitFcnV2_PRLTemplate_PeriphZYAM->Draw("same");
                               // draw the legend
                               TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                               legend->SetTextFont(72);
                               legend->SetTextSize(0.04);
                                 Char_t message[80];
                                 sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate_PeriphZYAM->GetChisquare(),fitFcnV2_PRLTemplate_PeriphZYAM->GetNDF());
                                 legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                          if(res->CovMatrixStatus() == 3){
                                     sprintf(message,"The fit is a success");
                                 }
                                 else{
                                     sprintf(message,"The fit is a failure");
                                 }
                                 legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                               legend->AddEntry(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin],"Data","lpe");
                               legend->Draw();
                          
                          V2_PRLPeriphZYAM[ptbin] = par[1];
                          errV2_PRLPeriphZYAM[ptbin] = fitFcnV2_PRLTemplate_PeriphZYAM->GetParError(1);
                          F_PRLPeriphZYAM[ptbin] = par[2];
                          errF_PRLPeriphZYAM[ptbin] = fitFcnV2_PRLTemplate_PeriphZYAM->GetParError(2);
                          

                      }
            
            
            
            
            
            
            
    
        }
        
    }

    
    
    double PtMiddle[NbPtBins] = {0};
    double PtErrorSize[NbPtBins] = {0};
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        PtMiddle[ptbin] = (PtBins[ptbin]+PtBins[ptbin+1])/2 - 0.02*3;
        PtErrorSize[ptbin] = (PtBins[ptbin+1]-PtBins[ptbin])/2;
    }
    
    TCanvas* cV2Pt = new TCanvas;
    cV2Pt->cd();
    
    TGraphErrors *grV2_wrt_Pt1 = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext1,PtErrorSize,errV2_Ext1);
             // TGraph *gr3 = new TGraph (n, K3, chi);
              grV2_wrt_Pt1->SetTitle("V2 JPsi-Tracklet wrt Pt for different extraction methods");
              grV2_wrt_Pt1->GetXaxis()->SetTitle("Pt (GeV/c)");
              grV2_wrt_Pt1->GetYaxis()->SetTitle("V2 (JPsi-tracklet)");
            grV2_wrt_Pt1->SetMarkerColor(4);
            grV2_wrt_Pt1->SetLineColor(4);
            grV2_wrt_Pt1->SetMarkerStyle(5);
              grV2_wrt_Pt1->Draw("AP");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_Pt2 = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext2,PtErrorSize,errV2_Ext2);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_Pt2->SetMarkerColor(3);
    grV2_wrt_Pt2->SetLineColor(3);
    grV2_wrt_Pt2->SetMarkerStyle(4);
      grV2_wrt_Pt2->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_PtCvetanQuentin = new TGraphErrors(NbPtBins,PtMiddle,V2_CvetanQuentin,PtErrorSize,errV2_CvetanQuentin);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtCvetanQuentin->SetMarkerColor(28);
    grV2_wrt_PtCvetanQuentin->SetLineColor(28);
    grV2_wrt_PtCvetanQuentin->SetMarkerStyle(3);
      grV2_wrt_PtCvetanQuentin->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }

//    TGraphErrors *grV2_wrt_PtCvetanQuentinMe = new TGraphErrors(NbPtBins,PtMiddle,V2_CvetanQuentinMe,PtErrorSize,errV2_CvetanQuentinMe);
//     // TGraph *gr3 = new TGraph (n, K3, chi);
//    grV2_wrt_PtCvetanQuentinMe->SetMarkerColor(6);
//    grV2_wrt_PtCvetanQuentinMe->SetLineColor(6);
//    grV2_wrt_PtCvetanQuentinMe->SetMarkerStyle(4);
//      grV2_wrt_PtCvetanQuentinMe->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_PtZYAM = new TGraphErrors(NbPtBins,PtMiddle,V2_ZYAM,PtErrorSize,errV2_ZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtZYAM->SetMarkerColor(7);
    grV2_wrt_PtZYAM->SetLineColor(7);
    grV2_wrt_PtZYAM->SetMarkerStyle(42);
      grV2_wrt_PtZYAM->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_PtPRL = new TGraphErrors(NbPtBins,PtMiddle,V2_PRL,PtErrorSize,errV2_PRL);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtPRL->SetMarkerColor(2);
    grV2_wrt_PtPRL->SetLineColor(2);
    grV2_wrt_PtPRL->SetMarkerStyle(26);
      grV2_wrt_PtPRL->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_PtPRLPeriphZYAM = new TGraphErrors(NbPtBins,PtMiddle,V2_PRLPeriphZYAM,PtErrorSize,errV2_PRLPeriphZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtPRLPeriphZYAM->SetMarkerColor(1);
    grV2_wrt_PtPRLPeriphZYAM->SetLineColor(1);
    grV2_wrt_PtPRLPeriphZYAM->SetMarkerStyle(25);
      grV2_wrt_PtPRLPeriphZYAM->Draw("P");
    
    TLegend *legendo=new TLegend(0.12,0.80,0.40,0.90);
             legendo->SetTextFont(41);
             legendo->SetTextSize(0.02);
           legendo->AddEntry(grV2_wrt_Pt1,"V2_Ext1 (pPb-like)");
            legendo->AddEntry(grV2_wrt_Pt2,"V2_Ext2 (pPb-like)");
            legendo->AddEntry(grV2_wrt_PtCvetanQuentin,"V2_CvetanQuentin (Template + ZYAM Periph)");
          //  legendo->AddEntry(grV2_wrt_PtCvetanQuentinMe,"V2_CvetanQuentinMe");
            legendo->AddEntry(grV2_wrt_PtZYAM,"V2_ZYAM");
            legendo->AddEntry(grV2_wrt_PtPRL,"V2_ATLAS Template");
            legendo->AddEntry(grV2_wrt_PtPRLPeriphZYAM,"V2_ATLAS Template and ZYAM)");
             legendo->Draw();
    
    cV2Pt->Update();
   // TLine *l=new TLine(cV2Pt->GetUxmin(),0.0,cV2Pt->GetUxmax(),0.0);
    TLine *l=new TLine(0.0,0.0,12.0,0.0);
    l->SetLineColor(kBlack);
    l->SetLineWidth(1);
    l->SetLineStyle(9);
    l->Draw();
    
    
    
    
    
    double PtMiddleF[NbPtBins] = {0};
    double PtErrorSizeF[NbPtBins] = {0};
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        PtMiddleF[ptbin] = (PtBins[ptbin]+PtBins[ptbin+1])/2 - 0.02*1;
        PtErrorSizeF[ptbin] = (PtBins[ptbin+1]-PtBins[ptbin])/2;
    }
    
    TCanvas* cFPt = new TCanvas;
    cFPt->cd();
    
    TGraphErrors *grF_wrt_PtCvetanQuentin = new TGraphErrors(NbPtBins,PtMiddleF,F_CvetanQuentin,PtErrorSize,errF_CvetanQuentin);
             // TGraph *gr3 = new TGraph (n, K3, chi);
              grF_wrt_PtCvetanQuentin->SetTitle("F of V2 JPsi-Tracklet extraction wrt Pt for different extraction methods");
              grF_wrt_PtCvetanQuentin->GetXaxis()->SetTitle("Pt (GeV/c)");
              grF_wrt_PtCvetanQuentin->GetYaxis()->SetTitle("F");
            grF_wrt_PtCvetanQuentin->SetMarkerColor(28);
            grF_wrt_PtCvetanQuentin->SetLineColor(28);
            grF_wrt_PtCvetanQuentin->SetMarkerStyle(3);
              grF_wrt_PtCvetanQuentin->Draw("AP");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddleF[ptbin] += 0.02;
    }
    
    TGraphErrors *grF_wrt_PtPRL = new TGraphErrors(NbPtBins,PtMiddleF,F_PRL,PtErrorSize,errF_PRL);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grF_wrt_PtPRL->SetMarkerColor(2);
    grF_wrt_PtPRL->SetLineColor(2);
    grF_wrt_PtPRL->SetMarkerStyle(26);
      grF_wrt_PtPRL->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddleF[ptbin] += 0.02;
    }
    
    TGraphErrors *grF_wrt_PtPRLPeriphZYAM = new TGraphErrors(NbPtBins,PtMiddleF,F_PRLPeriphZYAM,PtErrorSize,errF_PRLPeriphZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grF_wrt_PtPRLPeriphZYAM->SetMarkerColor(1);
    grF_wrt_PtPRLPeriphZYAM->SetLineColor(1);
    grF_wrt_PtPRLPeriphZYAM->SetMarkerStyle(25);
      grF_wrt_PtPRLPeriphZYAM->Draw("P");
    
    TLegend *legendf=new TLegend(0.12,0.80,0.40,0.90);
             legendf->SetTextFont(41);
             legendf->SetTextSize(0.02);
            legendf->AddEntry(grF_wrt_PtCvetanQuentin,"F_CvetanQuentin (Template + ZYAM Periph)");
            legendf->AddEntry(grF_wrt_PtPRL,"F_ATLAS Template");
            legendf->AddEntry(grF_wrt_PtPRLPeriphZYAM,"F_ATLAS Template and ZYAM");
             legendf->Draw();
    
    cFPt->Update();
    //TLine *lF=new TLine(cFPt->GetUxmin(),1.0,cFPt->GetUxmax(),1.0);
    TLine *lF=new TLine(0.,1.0,12.0,1.0);
    lF->SetLineColor(kBlack);
    lF->SetLineWidth(1);
    lF->SetLineStyle(9);
    lF->Draw();
    
    
    
    // Prendre les résultats sur V2JPsi-Tkl et utiliser ceux sur v2tkl pour trouver le v2JPsi pour chaque méthode
    
    //Résultats TKL 0-5% 40-90% - Baseline central 0
    
    double v2ClassiqueTKL = 0.065334;
    double errv2ClassiqueTKL = 0.000541066;
    double v2CvetanQuentinTKL = 0.0662997;
    double errv2CvetanQuentinTKL = 0.000409345;
    double v2CvetanQuentinMeTKL = 0.0654186;
    double errv2CvetanQuentinMeTKL = 0.000414581;
    double v2ZYAMTKL = 0.0653334;
    double errv2ZYAMTKL = 0.00104418;
    double v2PRLTKL = 0.131845;
    double errv2PRLTKL = 0.00186927;
    double v2PRLPeriphZYAMTKL = 0.0673054;
    double errv2PRLPeriphZYAMTKL = 0.000405759;

    //
    
    double v2_Ext1[NbPtBins] = {0};
    double v2_Ext2[NbPtBins] = {0};
    double v2_CvetanQuentin[NbPtBins] = {0};
    double v2_CvetanQuentinMe[NbPtBins] = {0};
    double v2_ZYAM[NbPtBins] = {0};
    double v2_PRL[NbPtBins] = {0};
    double v2_PRLPeriphZYAM[NbPtBins] = {0};
    double errv2_Ext1[NbPtBins] = {0};
    double errv2_Ext2[NbPtBins] = {0};
    double errv2_CvetanQuentin[NbPtBins] = {0};
    double errv2_CvetanQuentinMe[NbPtBins] = {0};
    double errv2_ZYAM[NbPtBins] = {0};
    double errv2_PRL[NbPtBins] = {0};
    double errv2_PRLPeriphZYAM[NbPtBins] = {0};
    
    for(int pt_idx=0; pt_idx < NbPtBins; pt_idx++){
        v2_Ext1[pt_idx] = V2_Ext1[pt_idx]/v2ClassiqueTKL;
        v2_Ext2[pt_idx] = V2_Ext2[pt_idx]/v2ClassiqueTKL;
        v2_CvetanQuentin[pt_idx] = V2_CvetanQuentin[pt_idx]/v2CvetanQuentinTKL;
        v2_CvetanQuentinMe[pt_idx] = V2_CvetanQuentinMe[pt_idx]/v2CvetanQuentinMeTKL;
        v2_ZYAM[pt_idx] = V2_ZYAM[pt_idx]/v2ZYAMTKL;
        v2_PRL[pt_idx] = V2_PRL[pt_idx]/v2PRLTKL;
        v2_PRLPeriphZYAM[pt_idx] = V2_PRLPeriphZYAM[pt_idx]/v2PRLPeriphZYAMTKL;
    }
    
    for(int pt_idx=0; pt_idx < NbPtBins; pt_idx++){
        errv2_Ext1[pt_idx] = abs(sqrt(pow(errV2_Ext1[pt_idx]/V2_Ext1[pt_idx],2)+pow(errv2ClassiqueTKL/v2ClassiqueTKL,2))*v2_Ext1[pt_idx]);
        errv2_Ext2[pt_idx] = abs(sqrt(pow(errV2_Ext2[pt_idx]/V2_Ext2[pt_idx],2)+pow(errv2ClassiqueTKL/v2ClassiqueTKL,2))*v2_Ext2[pt_idx]);
        errv2_CvetanQuentin[pt_idx] = abs(sqrt(pow(errV2_CvetanQuentin[pt_idx]/V2_CvetanQuentin[pt_idx],2)+pow(errv2CvetanQuentinTKL/v2CvetanQuentinTKL,2))*v2_CvetanQuentin[pt_idx]);
        errv2_CvetanQuentinMe[pt_idx] = abs(sqrt(pow(errV2_CvetanQuentinMe[pt_idx]/V2_CvetanQuentinMe[pt_idx],2)+pow(errv2CvetanQuentinMeTKL/v2CvetanQuentinMeTKL,2))*v2_CvetanQuentinMe[pt_idx]);
        errv2_ZYAM[pt_idx] = abs(sqrt(pow(errV2_ZYAM[pt_idx]/V2_ZYAM[pt_idx],2)+pow(errv2ZYAMTKL/v2ZYAMTKL,2))*v2_ZYAM[pt_idx]);
        errv2_PRL[pt_idx] = abs(sqrt(pow(errV2_PRL[pt_idx]/V2_PRL[pt_idx],2)+pow(errv2PRLTKL/v2PRLTKL,2))*v2_PRL[pt_idx]);
        errv2_PRLPeriphZYAM[pt_idx] = abs(sqrt(pow(errV2_PRLPeriphZYAM[pt_idx]/V2_PRLPeriphZYAM[pt_idx],2)+pow(errv2PRLPeriphZYAMTKL/v2PRLPeriphZYAMTKL,2))*v2_PRLPeriphZYAM[pt_idx]);
    }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
         PtMiddle[ptbin] = (PtBins[ptbin]+PtBins[ptbin+1])/2 - 0.02*3;
         PtErrorSize[ptbin] = (PtBins[ptbin+1]-PtBins[ptbin])/2;
     }
     
     TCanvas* cv2Pt = new TCanvas;
     cv2Pt->cd();
     
     TGraphErrors *grv2_wrt_Pt1 = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext1,PtErrorSize,errv2_Ext1);
              // TGraph *gr3 = new TGraph (n, K3, chi);
               grv2_wrt_Pt1->SetTitle("v2 JPsi wrt Pt for different extraction methods");
               grv2_wrt_Pt1->GetXaxis()->SetTitle("Pt (GeV/c)");
               grv2_wrt_Pt1->GetYaxis()->SetTitle("v2 (JPsi)");
             grv2_wrt_Pt1->SetMarkerColor(4);
             grv2_wrt_Pt1->SetLineColor(4);
             grv2_wrt_Pt1->SetMarkerStyle(5);
               grv2_wrt_Pt1->Draw("AP");
     
     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.02;
     }
     
     TGraphErrors *grv2_wrt_Pt2 = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext2,PtErrorSize,errv2_Ext2);
      // TGraph *gr3 = new TGraph (n, K3, chi);
     grv2_wrt_Pt2->SetMarkerColor(3);
     grv2_wrt_Pt2->SetLineColor(3);
     grv2_wrt_Pt2->SetMarkerStyle(4);
       grv2_wrt_Pt2->Draw("P");
     
     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.02;
     }
     
     TGraphErrors *grv2_wrt_PtCvetanQuentin = new TGraphErrors(NbPtBins,PtMiddle,v2_CvetanQuentin,PtErrorSize,errv2_CvetanQuentin);
      // TGraph *gr3 = new TGraph (n, K3, chi);
     grv2_wrt_PtCvetanQuentin->SetMarkerColor(28);
     grv2_wrt_PtCvetanQuentin->SetLineColor(28);
     grv2_wrt_PtCvetanQuentin->SetMarkerStyle(3);
       grv2_wrt_PtCvetanQuentin->Draw("P");
     
     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.02;
     }

//     TGraphErrors *grv2_wrt_PtCvetanQuentinMe = new TGraphErrors(NbPtBins,PtMiddle,v2_CvetanQuentinMe,PtErrorSize,errv2_CvetanQuentinMe);
//      // TGraph *gr3 = new TGraph (n, K3, chi);
//     grv2_wrt_PtCvetanQuentinMe->SetMarkerColor(6);
//     grv2_wrt_PtCvetanQuentinMe->SetLineColor(6);
//     grv2_wrt_PtCvetanQuentinMe->SetMarkerStyle(4);
//       grv2_wrt_PtCvetanQuentinMe->Draw("P");
     
     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.02;
     }
     
     TGraphErrors *grv2_wrt_PtZYAM = new TGraphErrors(NbPtBins,PtMiddle,v2_ZYAM,PtErrorSize,errv2_ZYAM);
      // TGraph *gr3 = new TGraph (n, K3, chi);
     grv2_wrt_PtZYAM->SetMarkerColor(7);
     grv2_wrt_PtZYAM->SetLineColor(7);
     grv2_wrt_PtZYAM->SetMarkerStyle(42);
       grv2_wrt_PtZYAM->Draw("P");
     
     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.02;
     }
     
     TGraphErrors *grv2_wrt_PtPRL = new TGraphErrors(NbPtBins,PtMiddle,v2_PRL,PtErrorSize,errv2_PRL);
      // TGraph *gr3 = new TGraph (n, K3, chi);
     grv2_wrt_PtPRL->SetMarkerColor(2);
     grv2_wrt_PtPRL->SetLineColor(2);
     grv2_wrt_PtPRL->SetMarkerStyle(26);
       grv2_wrt_PtPRL->Draw("P");
     
     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.02;
     }
     
     TGraphErrors *grv2_wrt_PtPRLPeriphZYAM = new TGraphErrors(NbPtBins,PtMiddle,v2_PRLPeriphZYAM,PtErrorSize,errv2_PRLPeriphZYAM);
      // TGraph *gr3 = new TGraph (n, K3, chi);
     grv2_wrt_PtPRLPeriphZYAM->SetMarkerColor(1);
     grv2_wrt_PtPRLPeriphZYAM->SetLineColor(1);
     grv2_wrt_PtPRLPeriphZYAM->SetMarkerStyle(25);
       grv2_wrt_PtPRLPeriphZYAM->Draw("P");
     
     TLegend *legendov=new TLegend(0.12,0.80,0.40,0.90);
              legendov->SetTextFont(41);
              legendov->SetTextSize(0.02);
            legendov->AddEntry(grv2_wrt_Pt1,"v2_Ext1 (pPb-like)");
             legendov->AddEntry(grv2_wrt_Pt2,"v2_Ext2 (pPb-like)");
             legendov->AddEntry(grv2_wrt_PtCvetanQuentin,"v2_CvetanQuentin (Template + ZYAM Periph)");
           //  legendov->AddEntry(grv2_wrt_PtCvetanQuentinMe,"v2_CvetanQuentinMe");
             legendov->AddEntry(grv2_wrt_PtZYAM,"v2_ZYAM");
             legendov->AddEntry(grv2_wrt_PtPRL,"v2_ATLAS Template");
             legendov->AddEntry(grv2_wrt_PtPRLPeriphZYAM,"v2_ATLAS Template and ZYAM");
              legendov->Draw();
     
     cv2Pt->Update();
    // TLine *l=new TLine(cV2Pt->GetUxmin(),0.0,cV2Pt->GetUxmax(),0.0);
     TLine *lv=new TLine(0.0,0.0,12.0,0.0);
     lv->SetLineColor(kBlack);
     lv->SetLineWidth(1);
     lv->SetLineStyle(9);
     lv->Draw();
    

    if(doTracklets){

        TCanvas* c6TKL = new TCanvas;
            //Tracklets yield DeltaEta wrt DeltaPhi TH2 -> Projected for Periph and Central

            c6TKL->Divide(2,2);
            c6TKL->cd(1);
            YieldTkl_allC->Draw("E");
            c6TKL->cd(2);
            YieldTkl_Difference->Draw("E");
            c6TKL->cd(3);
            YieldTkl_Central->Draw("E");
            c6TKL->cd(4);
            YieldTkl_Periph->Draw("E");
            c6TKL->Draw();
            c6TKL->Modified();
            c6TKL->ForceUpdate();

        baselineTKL_periph = (YieldTkl_Periph->GetBinContent(BinZeroLeftTKL) + YieldTkl_Periph->GetBinContent(BinZeroLeftTKL+1))/2;
        errbaselineTKL_periph = sqrt(pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL),2) + pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL+1),2));
        
        for(int bin_idx = 1; bin_idx<=NbinsDeltaPhiTKL; bin_idx++){
            if(YieldTkl_Central->GetBinContent(bin_idx)<baselineTKL_central){
                baselineTKL_central = YieldTkl_Central->GetBinContent(bin_idx);
                errbaselineTKL_central = YieldTkl_Central->GetBinError(bin_idx);
            }
        }
        
        
        for(int phi_idx = 0; phi_idx<NbinsDeltaPhiTKL; phi_idx++){
            BaselineTkl_Periph->SetBinContent(phi_idx+1, baselineTKL_periph);
            BaselineTkl_Periph->SetBinError(phi_idx+1, errbaselineTKL_periph);
        }
        
        YieldTkl_Periph_MinusBaseline->Add(YieldTkl_Periph,BaselineTkl_Periph,1,-1);
        
        for(int phi_idx = 0; phi_idx<NbinsDeltaPhiTKL; phi_idx++){
            BaselineTkl_Central->SetBinContent(phi_idx+1, baselineTKL_central);
            BaselineTkl_Central->SetBinError(phi_idx+1, errbaselineTKL_central);
        }
        
        YieldTkl_Central_MinusBaseline->Add(YieldTkl_Central,BaselineTkl_Central,1,-1);
        
        
        TCanvas*cTKLCentralMinusBaseline=new TCanvas();
        cTKLCentralMinusBaseline->Divide(1,3);
        cTKLCentralMinusBaseline->cd(1);
        YieldTkl_Central->DrawCopy();
        cTKLCentralMinusBaseline->cd(2);
        BaselineTkl_Central->DrawCopy();
        cTKLCentralMinusBaseline->cd(3);
        YieldTkl_Central_MinusBaseline->DrawCopy();
        
        TCanvas*cTKLPeriphMinusBaseline=new TCanvas();
        cTKLPeriphMinusBaseline->Divide(1,3);
        cTKLPeriphMinusBaseline->cd(1);
        YieldTkl_Periph->DrawCopy();
        cTKLPeriphMinusBaseline->cd(2);
        BaselineTkl_Periph->DrawCopy();
        cTKLPeriphMinusBaseline->cd(3);
        YieldTkl_Periph_MinusBaseline->DrawCopy();

        TCanvas*c14TKL=new TCanvas();
        //Tracklets Yield difference wrt Phi fit
        
        
        TH1F *YieldTkl_Central_Cvetan = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_Cvetan");
        TH1F *YieldTkl_Central_CvetanMe = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_CvetanMe");
        TH1F *YieldTkl_Central_ZYAM = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_ZYAM");
        TH1F *YieldTkl_Central_PRLTemplate = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_PRLTemplate");
        TH1F *YieldTkl_Central_PRLTemplate_PeriphZYAM = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_PRLTemplate_PeriphZYAM");
        
        //Methode classique C-P
        {
        c14TKL->cd();
        // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
            TH1F *histo = YieldTkl_Difference;
          // create a TF1 with the range from 0 to 3 and 6 parameters
          TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV5,-TMath::Pi()/2,1.5*TMath::Pi(),3);
          fitFcnV2TKL->SetNpx(500);
          fitFcnV2TKL->SetLineWidth(4);
          fitFcnV2TKL->SetLineColor(kMagenta);
          // first try without starting values for the parameters
          // This defaults to 1 for each param.
          // this results in an ok fit for the polynomial function
          // however the non-linear part (lorenzian) does not
          // respond well.
           Double_t params[3] = {1,0.01,0.01};
          fitFcnV2TKL->SetParameters(params);
           TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
           TVirtualFitter::Fitter(histo)->SetPrecision();
        //  histo->Fit("fitFcn","0");
          // second try: set start values for some parameters

           fitFcnV2TKL->SetParName(0,"a0");
           fitFcnV2TKL->SetParName(1,"a1");
           fitFcnV2TKL->SetParName(2,"a2");
//            fitFcnV2TKL->SetParName(3,"a3");
//            fitFcnV2TKL->SetParName(4,"a4");
//        fitFcnV2TKL->SetParName(5,"a5");
//        fitFcnV2TKL->SetParName(6,"a6");
//        fitFcnV2TKL->SetParName(7,"a7");
//        fitFcnV2TKL->SetParName(8,"a8");
//        fitFcnV2TKL->SetParName(9,"a9");
//        fitFcnV2TKL->SetParName(10,"a10");
//        fitFcnV2TKL->SetParName(11,"a11");
//        fitFcnV2TKL->SetParLimits(6,-0.1,0.1);

          TFitResultPtr res = histo->Fit("fitFcnV2TKL","SBMERI+","ep");
          // improve the pictu
        //   std::cout << "integral error: " << integralerror << std::endl;
            Double_t par[3];
            fitFcnV2TKL->GetParameters(par);
          fitFcnV2TKL->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetTextFont(72);
          legend->SetTextSize(0.04);
            Char_t message[80];
            sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
            legend->AddEntry(fitFcnV2TKL,message);
        sprintf(message,"V2 Tkl-Tkl: %.4f / (%.4f + %.4f) = %.4f +- %.4f",par[2],par[0], baselineTKL_periph,par[2]/(par[0] + baselineTKL_periph),(par[2]/(par[0] + baselineTKL_periph)*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)+pow(errbaselineTKL_periph/baselineTKL_periph,2))));
        legend->AddEntry(fitFcnV2TKL,message);
        if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitFcnV2TKL,message);
          legend->AddEntry(histo,"Data","lpe");
          legend->Draw();
        }
        
        
        
        
        //Cvetan-Quentin fit
        
        TCanvas* cTKLCvetan = new TCanvas;
        cTKLCvetan->SetTitle("TKL Cvetan-Quentin Fit");
        cTKLCvetan->Divide(1,1);
        cTKLCvetan->cd(1);
        YieldTkl_Central->DrawCopy();
        
        {
            cTKLCvetan->cd(1);
                 TF1 *fitFcnV2_Cvetan = new TF1("fitFcnV2_Cvetan",CvetanFTKL,-TMath::Pi()/2,1.5*TMath::Pi(),4);
                 fitFcnV2_Cvetan->SetNpx(500);
                 fitFcnV2_Cvetan->SetLineWidth(4);
                 fitFcnV2_Cvetan->SetLineColor(kBlue);
                 // first try without starting values for the parameters
                 // This defaults to 1 for each param.
                 // this results in an ok fit for the polynomial function
                 // however the non-linear part (lorenzian) does not
                 // respond well.
                  Double_t params[4] = {1,0,0.01,1};
                 fitFcnV2_Cvetan->SetParameters(params);
                  TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->SetMaxIterations(10000);
                  TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->SetPrecision();
                gStyle->SetOptFit(1011);
               //  histo->Fit("fitFcn","0");
                 // second try: set start values for some parameters

                  fitFcnV2_Cvetan->SetParName(0,"V0");
                    fitFcnV2_Cvetan->SetParName(1,"V1");
                    fitFcnV2_Cvetan->SetParName(2,"V2");
                  fitFcnV2_Cvetan->SetParName(3,"F");
            
            fitFcnV2_Cvetan->SetParLimits(3,0.1,50);
            
            fitFcnV2_Cvetan->FixParameter(0,1);
            fitFcnV2_Cvetan->FixParameter(1,0);
          //  fitFcnV2_Cvetan->FixParameter(3,1);
                  

                 TFitResultPtr res = YieldTkl_Central_Cvetan->Fit("fitFcnV2_Cvetan","SBMERI+","ep");
                Double_t par[4];
                fitFcnV2_Cvetan->GetParameters(par);
                 // improve the pictu
               //   std::cout << "integral error: " << integralerror << std::endl;
                 fitFcnV2_Cvetan->Draw("same");
                 // draw the legend
                 TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                 legend->SetTextFont(72);
                 legend->SetTextSize(0.04);
                   Char_t message[80];
                   sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_Cvetan->GetChisquare(),fitFcnV2_Cvetan->GetNDF());
                   legend->AddEntry(fitFcnV2_Cvetan,message);
            if(res->CovMatrixStatus() == 3){
                       sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                   }
                   legend->AddEntry(fitFcnV2_Cvetan,message);
                 legend->AddEntry(YieldTkl_Central_Cvetan,"Data","lpe");
                 legend->Draw();

        }
        
        

    }
    
    
    
//        for (int i=0;i<6;i++){
//               sprintf(histoname,Form("bin%d_0",i+1)); // ZZZZZZZ
//               res = FittingAllInvMassBin(histoname, c11b_0, i);
//           }
//    for (int i=0;i<6;i++){
//        sprintf(histoname,Form("bin%d_1",i+1));
//        res = FittingAllInvMassBin(histoname, c11b_1, i);
//    }
//    for (int i=0;i<6;i++){
//        sprintf(histoname,Form("bin%d_2",i+1));
//        res = FittingAllInvMassBin(histoname, c11b_2, i);
//    }
    
    
    
    
    
//    for(int i=0; i<12; i++){
//  //  histoname = "Inv Mass, Phi: " + std::to_string(0) + " pi/6 to " + std::to_string(1) + " pi/6";
//        sprintf(histoname,"InvMassPhiBin_%d",i);
//        FittingAllInvMassPhiBin(histoname, cinvmassphibin, i);
//    }
    std::cout << "===== Dimuons counted =====" <<std::endl;
    std::cout<< "DimuCentralSeen " << DimuCentralSeen <<std::endl;
    std::cout<< "DimuPeriphSeen " << DimuPeriphSeen <<std::endl;
    std::cout<< "DimuSeenMassCut " << DimuSeenMassCut <<std::endl;
    std::cout<< "DimuSeenNoMassCut " << DimuSeenNoMassCut <<std::endl;
    std::cout << "EventRejected: " << EventRejected <<std::endl;
    std::cout << "EventPileUpVtx: " << EventPileUpVtx <<std::endl;
    std::cout << "EventPileUpMult: " << EventPileUpMult <<std::endl;
    std::cout << "EventNC: " << EventNC <<std::endl;
    std::cout << "cmul  " << cmul <<std::endl;
    std::cout <<"countsigma :" << countsigma <<std::endl;
}
    









// FITTING INVARIANT MASS METHODS

Double_t CvetanF(Double_t *x,Double_t *par)

{   int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double YMinusBp = Yields_Periph_1_MinusBaseline->GetBinContent(bintolook+1);
    return baseline_central*(par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + par[3]*YMinusBp; }

void ChisquareCvetanF(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-CvetanF(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t CvetanFPtBinned(Double_t *x,Double_t *par)

{   int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int ptbin = par[0];
    double YMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinContent(bintolook+1);
    return baseline_centralPtBinned[ptbin]*(par[1] + 2*par[2]*cos(x[0]) + 2*par[3]*cos(2*x[0])) + par[4]*YMinusBp; }

void ChisquareCvetanFPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-CvetanFPtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t CvetanFTKL(Double_t *x,Double_t *par)

{   int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double YMinusBp = YieldTkl_Periph_MinusBaseline->GetBinContent(bintolook+1);
    return baselineTKL_central*(par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + par[3]*YMinusBp; }

Double_t ZYAM(Double_t *x,Double_t *par)

{   int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double YMinusBp = Yields_Periph_1_MinusBaseline->GetBinContent(bintolook+1);
    return baseline_central + (par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + 1*YMinusBp; }

void ChisquareZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-ZYAM(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t ZYAMPtBinned(Double_t *x,Double_t *par)

{   int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int ptbin = par[0];
    double YMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinContent(bintolook+1);
    return baseline_centralPtBinned[ptbin] + (par[1] + 2*par[2]*cos(x[0]) + 2*par[3]*cos(2*x[0])) + 1*YMinusBp; }

void ChisquareZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-ZYAMPtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplate(Double_t *x,Double_t *par)

{   double integral_Yperiph = Yields_Periph_1->Integral(1,Yields_Periph_1->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double Yperiphbin = Yields_Periph_1->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[1]*Yperiphbin + (1+2*par[0]*cos(2*x[0]))*G; }

void ChisquarePRLTemplate(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-PRLTemplate(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplatePtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_Yperiph = Yields_Periph_1PtBinned[ptbin]->Integral(1,Yields_Periph_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double Yperiphbin = Yields_Periph_1PtBinned[ptbin]->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[2]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[2]*Yperiphbin + (1+2*par[1]*cos(2*x[0]))*G; }

void ChisquarePRLTemplatePtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-PRLTemplatePtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplate_RidgeAndZero(Double_t *x,Double_t *par)

{   double integral_Yperiph = Yields_Periph_1->Integral(1,Yields_Periph_1->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double Yperiphbin = baseline_periph;
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[1]*Yperiphbin + (1+2*par[0]*cos(2*x[0]))*G; }

Double_t PRLTemplate_RidgeAndZeroPtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_Yperiph = Yields_Periph_1PtBinned[ptbin]->Integral(1,Yields_Periph_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double Yperiphbin = baseline_periphPtBinned[ptbin];
    double G = (integral_Yreal-(par[2]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[2]*Yperiphbin + (1+2*par[1]*cos(2*x[0]))*G; }

Double_t PRLTemplate_PeriphAndG(Double_t *x,Double_t *par)

{   double integral_Yperiph = Yields_Periph_1->Integral(1,Yields_Periph_1->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double Yperiphbin = Yields_Periph_1->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[1]*Yperiphbin + G; }

Double_t PRLTemplate_PeriphAndGPtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_Yperiph = Yields_Periph_1PtBinned[ptbin]->Integral(1,Yields_Periph_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double Yperiphbin = Yields_Periph_1PtBinned[ptbin]->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[2]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[2]*Yperiphbin + G; }

Double_t PRLTemplate_PeriphZYAM(Double_t *x,Double_t *par)

{   double integral_YperiphMinusBp = Yields_Periph_1_MinusBaseline->Integral(1,Yields_Periph_1_MinusBaseline->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double YMinusBp = Yields_Periph_1_MinusBaseline->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_YperiphMinusBp))/(2*TMath::Pi());
    
    return par[1]*YMinusBp + (1+2*par[0]*cos(2*x[0]))*G; }

void ChisquarePRLTemplate_PeriphZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-PRLTemplate_PeriphZYAM(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}


Double_t PRLTemplate_PeriphZYAMPtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_YperiphMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->Integral(1,Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    double YMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[2]*integral_YperiphMinusBp))/(2*TMath::Pi());
    
    return par[2]*YMinusBp + (1+2*par[1]*cos(2*x[0]))*G; }

void ChisquarePRLTemplate_PeriphZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-PRLTemplate_PeriphZYAMPtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}


Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par)
// Par 0->7: signal, 8->11 Bkg, 12: v2 JPsi, 13->15: V2 bkg
{ return (JPsiCrystalBallExtended(x,par)*par[12] + (ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*x[0]*x[0] + par[14]*x[0] + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t BackFcnV2(Double_t *x,Double_t *par)
{return ((ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*x[0]*x[0] + par[14]*x[0] + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t BackFcnV2Poly(Double_t *x,Double_t *par)
{return (par[13]*x[0]*x[0] + par[14]*x[0] + par[15]) ;}

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
    fitFcn->SetParLimits(8,0.1,10000000);
    fitFcn->SetParLimits(9,0.01,20);
    fitFcn->SetParLimits(10,0.01,10000000);
    fitFcn->SetParLimits(11,0.01,50);
    
//   fitFcn->FixParameter(1,3.096916); // Mean x core
//    fitFcn->FixParameter(2,0.07);
//    fitFcn->FixParameter(3,0.9);
//    fitFcn->FixParameter(4,10);
//    fitFcn->FixParameter(5,2);
//    fitFcn->FixParameter(6,15);
//    fitFcn->SetParameter(7,100);
//    fitFcn->SetParameter(8,1000);
//    fitFcn->SetParameter(9,0.5);
//    fitFcn->SetParameter(11,10);
    
    fitFcn->FixParameter(1,3.096916); // Mean x core
    fitFcn->FixParameter(2,0.07);
    fitFcn->FixParameter(3,0.883);
    fitFcn->FixParameter(4,9.940);
    fitFcn->FixParameter(5,1.832);
    fitFcn->FixParameter(6,15.323);
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
    Double_t par[12];
    
   TFitResultPtr res = histo->Fit("fitFcn","SBMER","ep");
 //  res = histo->Fit("fitFcn","SBMERQ","ep");
      fitFcn->ReleaseParameter(1);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    fitFcn->ReleaseParameter(2);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    fitFcn->ReleaseParameter(7);
   res = histo->Fit("fitFcn","SBMERQ","ep");
   res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(3);
//       res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(5);
//   res = histo->Fit("fitFcn","SBMER","ep");
//   fitFcn->ReleaseParameter(4);
//    res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(6);
//   res = histo->Fit("fitFcn","SBMER","ep");
//    res = histo->Fit("fitFcn","SBMER","ep");

    
   // improve the picture:
    histo->SetStats(kTRUE);
   TF1 *backFcn = new TF1("backFcn",ExpBkg,2.1,5.1,4);
   backFcn->SetLineColor(kRed);
   TF1 *signalFcnJPsi = new TF1("signalFcnJPsi",JPsiCrystalBallExtended,2.1,5.1,8);
   TF1 *signalFcnPsi2S = new TF1("signalFcnPsi2S",Psi2SCrystalBallExtended,2.1,5.1,8);
   TPaveText *pave = new TPaveText(0.15,0.5,0.3,0.65,"brNDC");
   signalFcnJPsi->SetLineColor(kBlue);
   signalFcnJPsi->SetNpx(500);
    signalFcnPsi2S->SetLineColor(kGreen);
    signalFcnPsi2S->SetNpx(500);
   //Double_t par[12];
   // writes the fit results into the par array
   // gStyle->SetOptFit(1011);
    gStyle->SetOptFit(0);
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
    sprintf(str, "N_{JPsi} %i +/- %i", int(integral), int(integralerror));
   pave->AddText(str);
    sprintf(str, "M_{Psi2S} = %f, Sig_{Psi2S} = %f", int(mPsip*1000)/1000., int(ratSigma*par[2]*1000)/1000.);
  //  pave->AddText(str);
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
