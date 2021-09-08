 #include <TBrowser.h>
 #include <TBufferFile.h>
 #include <TCanvas.h>
 #include <TClass.h>
 #include <TMath.h>
#include <TLatex.h>
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
#include <sys/types.h>
#include <sys/stat.h>

 //#include "TTreeReader.h"
 //#include "Event.h"

 # include "TF1.h"
 # include "TF2.h"
 # include "TProfile.h"
# include "TGraph2D.h"
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

//Re-tourner le post-processing et les fits à partir de données TKL déjà enregistrées

 // FUNCTIONS

 void FitTrainingTKL();
 
Double_t CvetanFTKL(Double_t *x,Double_t *par);
void ChisquareCvetanFTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
 Double_t FourierV2(Double_t *x,Double_t *par);
 Double_t FourierV5(Double_t *x,Double_t *par);
Double_t ZYAMTKL(Double_t *x,Double_t *par);
void ChisquareZYAMTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplateTKL(Double_t *x,Double_t *par);
void ChisquarePRLTemplateTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate_RidgeAndZeroTKL(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphAndGTKL(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphZYAMTKL(Double_t *x,Double_t *par);
void ChisquarePRLTemplate_PeriphZYAMTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate_PeriphZYAM_RidgeTKL(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphZYAM_PeriphAndGTKL(Double_t *x,Double_t *par);
Double_t DiffPRL(Double_t *x,Double_t *par);

Double_t mJpsi =  3.096916;
 Double_t mPsip =  3.686108;
// Double_t ratMass = 1.01;
 Double_t ratSigma = 1.05;
Int_t npfits;


TH1F* YieldTkl_Central(NULL);
TH1F* YieldTkl_Periph(NULL);
TH1F* BaselineTkl_Central(NULL);
TH1F* BaselineTkl_Periph(NULL);
TH1F* YieldTkl_Periph_MinusBaseline(NULL);
TH1F* YieldTkl_Central_MinusBaseline(NULL);

double baselineTKL_periph = 1;
double errbaselineTKL_periph = 1;
double baselineTKL_central = 9999;
double errbaselineTKL_central = 1;

double V2ClassiqueTKL = 0;
double errV2ClassiqueTKL = 0;
double V2ClassiqueTKL_noZYAM = 0;
double errV2ClassiqueTKL_noZYAM = 0;
double V2CvetanTKL = 0;
double errV2CvetanTKL = 0;
double V2CvetanMeTKL = 0;
double errV2CvetanMeTKL = 0;
double V2ZYAMTKL = 0;
double errV2ZYAMTKL = 0;
double V2PRLTKL = 0;
double errV2PRLTKL = 0;
double V2PRL_PeriphZYAMTKL = 0;
double errV2PRL_PeriphZYAMTKL = 0;

double FCvetanTKL = 0;
double errFCvetanTKL = 0;
double FPRLTKL = 0;
double errFPRLTKL = 0;
double FPRL_PeriphZYAMTKL = 0;
double errFPRL_PeriphZYAMTKL = 0;


double v2ClassiqueTKL = 0;
double errv2ClassiqueTKL = 0;
double v2ClassiqueTKL_noZYAM = 0;
double errv2ClassiqueTKL_noZYAM = 0;
double v2CvetanTKL = 0;
double errv2CvetanTKL = 0;
double v2CvetanMeTKL = 0;
double errv2CvetanMeTKL = 0;
double v2ZYAMTKL = 0;
double errv2ZYAMTKL = 0;
double v2PRLTKL = 0;
double errv2PRLTKL = 0;
double v2PRL_PeriphZYAMTKL = 0;
double errv2PRL_PeriphZYAMTKL = 0;

double fourierCentral[3] = {0,0,0};
double fourierPeriph[3] = {0,0,0};

Char_t FitFileName[200];

void FitTrainingTKL(){
    
    Char_t FolderName[200];
    Char_t CanvasName[200];
    sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_TKL_16h_PercentileMethodSPDTracklets_0-5_40-100_pt0-12_Eta1.6.root");
        sprintf(FolderName,"~/Desktop/Images JavierAnalysis/2021 mai/NewAnalysis_16h10_TKL_QGPFrance_0-5_40-100_pt0-12/FitTrainingTKLa");
    
    TH1::SetDefaultSumw2();
        bool doTracklets = kTRUE;
        bool doMixedEvents = kTRUE;
    bool Extraction2Combined = kFALSE;
        bool CombineFits = kFALSE;
    bool PtBinned = kFALSE;
        
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
        
        const int NbinsDeltaPhiTKL = 48;
        double SizeBinDeltaPhiTKL = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhiTKL;
        int BinZeroLeftTKL = floor((0-MinDeltaPhi)*NbinsDeltaPhiTKL/(2*TMath::Pi()));
        
        const int NbBinsCent = 14;
        const int NbBinsZvtx = 20;
    
    // *************************
    // Initialiser les graphes *
    // *************************
    
    
        TH1F* YieldTkl_allC(NULL);
//        TH1F* YieldTkl_Central(NULL);
//        TH1F* YieldTkl_Periph(NULL);
        TH1F* YieldTkl_Difference(NULL);
    TH1F* YieldTkl_FDifference(NULL);
    
    TH1F *YieldTkl_Central_CvetanPtBinned(NULL);
    TH1F *YieldTkl_Central_CvetanMePtBinned(NULL);
    TH1F *YieldTkl_Central_ZYAMPtBinned(NULL);
    TH1F *YieldTkl_Central_PRLTemplatePtBinned(NULL);
    TH1F *YieldTkl_Central_PRLTemplate_PeriphZYAMPtBinned(NULL);

            

        YieldTkl_allC = new TH1F("YieldTkl_allC",
                             "Yields Correlation tkl-tkl wrt #Delta#phi, all C, all mass",
                             NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
           YieldTkl_allC->SetXTitle("Correlation #Delta#phi (rad)");
            YieldTkl_allC->SetYTitle("Yields");
        YieldTkl_Central = new TH1F("YieldTkl_Central",
                          "Yields Correlation tkl-tkl wrt #Delta#phi, Central",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldTkl_Central->SetXTitle("Correlation #Delta#phi (rad)");
        YieldTkl_Central->SetYTitle("Yields_{Central}");
        YieldTkl_Periph = new TH1F("YieldTkl_Periph",
                          "Yields Correlation tkl-tkl wrt #Delta#phi, Periph",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldTkl_Periph->SetXTitle("Correlation #Delta#phi (rad)");
        YieldTkl_Periph->SetYTitle("Yields_{Periph}");
        YieldTkl_Difference = new TH1F("YieldTkl_Difference",
                          "Subtracted Yields Correlation tkl-tkl wrt #Delta#phi, Central-Periph, all mass",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldTkl_Difference->SetXTitle("Correlation #Delta#phi (rad)");
        YieldTkl_Difference->SetYTitle("Yields_{Subtracted}");
    
    YieldTkl_FDifference = new TH1F("YieldTkl_FDifference",
                      "Subtracted Yields Correlation tkl-tkl wrt #Delta#phi, Central-F*Periph, all mass",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_FDifference->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_FDifference->SetYTitle("Yields_{Central} - F*Yields_{Periph}");
    
    BaselineTkl_Central = new TH1F("BaselineTkl_Central",
                     "Baseline of tkl-tkl in Central collisions wrt #Delta#phi",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    BaselineTkl_Central->SetXTitle("#Delta#phi (rad)");
    BaselineTkl_Central->SetYTitle("Baseline_{Central}");
    
    BaselineTkl_Periph = new TH1F("BaselineTkl_Periph",
                     "Baseline of tkl-tkl in Periph collisions wrt #Delta#phi",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    BaselineTkl_Periph->SetXTitle("#Delta#phi (rad)");
    BaselineTkl_Periph->SetYTitle("Baseline_{Periph}");
    
    YieldTkl_Central_MinusBaseline = new TH1F("YieldTkl_Central_MinusBaseline",
                     "Yield of tkl-tkl in Central collisions wrt #Delta#phi - Subtracted Baseline",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_Central_MinusBaseline->SetXTitle("#Delta#phi (rad)");
    YieldTkl_Central_MinusBaseline->SetYTitle("Yield_{Central} - Baseline_{Central}");
    
    YieldTkl_Periph_MinusBaseline = new TH1F("YieldTkl_Periph_MinusBaseline",
                     "Yield of tkl-tkl in Periph collisions wrt #Delta#phi - Subtracted Baseline",
                     NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_Periph_MinusBaseline->SetXTitle("#Delta#phi (rad)");
    YieldTkl_Periph_MinusBaseline->SetYTitle("Yield_{Periph} - Baseline_{Periph}");
    
    
    
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
        
    
        int DimuC = 0;
        int DimuP = 0;
        int TklC = 0;
        int TklP = 0;
        int NormTklCentral = 0;
        int NormTklPeriph = 0;
        int countsigma = 0;
    
    // Récup histos
    
    TFile *filerec = new TFile(FitFileName);
    
    YieldTkl_Central = (TH1F*)filerec->Get("YieldTkl_Central");
    YieldTkl_Periph = (TH1F*)filerec->Get("YieldTkl_Periph");
    YieldTkl_Difference = (TH1F*)filerec->Get("YieldTkl_Difference");
    
    YieldTkl_FDifference->Add(YieldTkl_Central, YieldTkl_Periph, 1, -2.678);


    
    cout << "START ALL"<<endl;


        TCanvas* c6TKL = new TCanvas;
            //Tracklets yield DeltaEta wrt DeltaPhi TH2 -> Projected for Periph and Central
    gStyle->SetOptStat("n");

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

    //Calculs baselinesTKL
    
    
        baselineTKL_periph = (YieldTkl_Periph->GetBinContent(BinZeroLeftTKL) + YieldTkl_Periph->GetBinContent(BinZeroLeftTKL+1))/2;
        errbaselineTKL_periph = sqrt(pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL),2) + pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL+1),2));
//
//            baselineTKL_central = (YieldTkl_Central->GetBinContent(BinZeroLeftTKL) + YieldTkl_Central->GetBinContent(BinZeroLeftTKL+1))/2;
//           errbaselineTKL_central = sqrt(pow(YieldTkl_Central->GetBinError(BinZeroLeftTKL),2) + pow(YieldTkl_Central->GetBinError(BinZeroLeftTKL+1),2));
        
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
        TCanvas*c14TKLCentral=new TCanvas();
        TCanvas*c14TKLPeriph=new TCanvas();
        //Tracklets Yield difference wrt Phi fit
        
        
        TH1F *YieldTkl_Central_Cvetan = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_Cvetan");
        TH1F *YieldTkl_Central_CvetanMe = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_CvetanMe");
        TH1F *YieldTkl_Central_ZYAM = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_ZYAM");
        TH1F *YieldTkl_Central_PRLTemplate = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_PRLTemplate");
        TH1F *YieldTkl_Central_PRLTemplate_PeriphZYAM = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_PRLTemplate_PeriphZYAM");
    
    {
            c14TKLCentral->cd();
            // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
                TH1F *histo = YieldTkl_Central;
              // create a TF1 with the range from 0 to 3 and 6 parameters
              TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV2,-TMath::Pi()/2,1.5*TMath::Pi(),3);
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
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
                st->SetX1NDC(0.8); //new x start position
                st->SetY1NDC(0.8); //new x end position
              // improve the pictu
            //   std::cout << "integral error: " << integralerror << std::endl;
                fitFcnV2TKL->GetParameters(fourierCentral);
              fitFcnV2TKL->Draw("same");
              // draw the legend
              TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
              legend->SetFillColorAlpha(kWhite, 0.);
              legend->SetBorderSize(0);
                legend->SetTextFont(42);
              legend->SetTextSize(0.03);
                Char_t message[80];
                sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
                legend->AddEntry(fitFcnV2TKL,message);
            legend->AddEntry(fitFcnV2TKL,message);
            if(res->CovMatrixStatus() == 3){
                      // sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                       legend->AddEntry(fitFcnV2TKL,message);
                   }
              legend->AddEntry(histo,"Data","lpe");
              legend->Draw();
                
            }
    
    {
            c14TKLPeriph->cd();
            // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
                TH1F *histo = YieldTkl_Periph;
              // create a TF1 with the range from 0 to 3 and 6 parameters
              TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV2,-TMath::Pi()/2,1.5*TMath::Pi(),3);
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
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
                st->SetX1NDC(0.8); //new x start position
                st->SetY1NDC(0.8); //new x end position
              // improve the pictu
            //   std::cout << "integral error: " << integralerror << std::endl;
                fitFcnV2TKL->GetParameters(fourierPeriph);
              fitFcnV2TKL->Draw("same");
              // draw the legend
              TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
              legend->SetFillColorAlpha(kWhite, 0.);
              legend->SetBorderSize(0);
                legend->SetTextFont(42);
              legend->SetTextSize(0.03);
                Char_t message[80];
                sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
                legend->AddEntry(fitFcnV2TKL,message);
            legend->AddEntry(fitFcnV2TKL,message);
            if(res->CovMatrixStatus() == 3){
                     //  sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                       legend->AddEntry(fitFcnV2TKL,message);
                   }
              legend->AddEntry(histo,"Data","lpe");
              legend->Draw();
                
            }
        
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
        //    fitFcnV2TKL->FixParameter(1,0);

          TFitResultPtr res = histo->Fit("fitFcnV2TKL","SBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
            st->SetX1NDC(0.8); //new x start position
            st->SetY1NDC(0.8); //new x end position
          // improve the pictu
        //   std::cout << "integral error: " << integralerror << std::endl;
            Double_t par[3];
            fitFcnV2TKL->GetParameters(par);
          fitFcnV2TKL->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetFillColorAlpha(kWhite, 0.);
          legend->SetBorderSize(0);
            legend->SetTextFont(42);
          legend->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
            legend->AddEntry(fitFcnV2TKL,message);
        sprintf(message,"V_{2,tkl-tkl}: #frac{%.4f}{%.4f + %.4f} = %.5f +- %.5f",par[2],par[0], baselineTKL_periph,par[2]/(par[0] + baselineTKL_periph),(par[2]/(par[0] + baselineTKL_periph)*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)+pow(errbaselineTKL_periph/baselineTKL_periph,2))));
        legend->AddEntry(fitFcnV2TKL,message);
//        sprintf(message,"V2 Tkl-Tkl: %.4f / (%.4f) = %.5f +- %.5f",par[2],par[0],par[2]/(par[0]),(par[2]/(par[0])*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2))));
//        legend->AddEntry(fitFcnV2TKL,message);
        if(res->CovMatrixStatus() == 3){
                  // sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
                   legend->AddEntry(fitFcnV2TKL,message);
               }
          legend->AddEntry(histo,"Data","lpe");
          legend->Draw();
            
            V2ClassiqueTKL = par[2]/(par[0] + baselineTKL_periph);
            errV2ClassiqueTKL = (par[2]/(par[0] + baselineTKL_periph)*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)+pow(errbaselineTKL_periph/baselineTKL_periph,2)));
            
            v2ClassiqueTKL = sqrt(V2ClassiqueTKL);
            errv2ClassiqueTKL = 0.5*(errV2ClassiqueTKL/V2ClassiqueTKL)*v2ClassiqueTKL;
            
            
            V2ClassiqueTKL_noZYAM = par[2]/(par[0]);
            errV2ClassiqueTKL_noZYAM = (par[2]/(par[0])*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)));
            
            v2ClassiqueTKL_noZYAM = sqrt(V2ClassiqueTKL_noZYAM);
            errv2ClassiqueTKL_noZYAM = 0.5*(errV2ClassiqueTKL_noZYAM/V2ClassiqueTKL_noZYAM)*v2ClassiqueTKL_noZYAM;
            
            cout << "===== Méthode Classique =====" <<endl;
            cout << "V2 TKL +/- err V2 TKL: " << V2ClassiqueTKL << " +/- " << errV2ClassiqueTKL <<endl;
            cout << "v2 TKL +/- err v2 TKL: " << v2ClassiqueTKL << " +/- " << errv2ClassiqueTKL <<endl;
            cout << "=============================" <<endl;
            
            cout << "===== Méthode Classique no ZYAM Periph =====" <<endl;
            cout << "V2 TKL +/- err V2 TKL: " << V2ClassiqueTKL_noZYAM << " +/- " << errV2ClassiqueTKL_noZYAM <<endl;
            cout << "v2 TKL +/- err v2 TKL: " << v2ClassiqueTKL_noZYAM << " +/- " << errv2ClassiqueTKL_noZYAM <<endl;
            cout << "=============================" <<endl;
        }
        
    sprintf(CanvasName,"%s/YieldsTKL_Difference.pdf",FolderName);
    c14TKL->SaveAs(CanvasName);
    sprintf(CanvasName,"%s/YieldsTKL_Central.pdf",FolderName);
    c14TKLCentral->SaveAs(CanvasName);
    sprintf(CanvasName,"%s/YieldsTKL_Periph.pdf",FolderName);
    c14TKLPeriph->SaveAs(CanvasName);
    
        
        
        
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
                 fitFcnV2_Cvetan->SetLineColor(kOrange+3);
                 // first try without starting values for the parameters
                 // This defaults to 1 for each param.
                 // this results in an ok fit for the polynomial function
                 // however the non-linear part (lorenzian) does not
                 // respond well.
                  Double_t params[4] = {1,0,0.01,1};
                 fitFcnV2_Cvetan->SetParameters(params);
                  TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->SetMaxIterations(10000);
                  TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->SetPrecision();
                TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->SetFCN(ChisquareCvetanFTKL);
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
        //    fitFcnV2_Cvetan->FixParameter(3,1);
                  

                 TFitResultPtr res = YieldTkl_Central_Cvetan->Fit("fitFcnV2_Cvetan","USBMERI+","ep");
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)YieldTkl_Central_Cvetan->FindObject("stats");
                st->SetX1NDC(0.8); //new x start position
                st->SetY1NDC(0.8); //new x end position
            
            double chi2, edm, errdef;
            int nvpar, nparx;
             TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->GetStats(chi2,edm,errdef,nvpar,nparx);
             fitFcnV2_Cvetan->SetChisquare(chi2);
             int ndf = npfits-nvpar;
             fitFcnV2_Cvetan->SetNDF(ndf);
            
                Double_t par[4];
                fitFcnV2_Cvetan->GetParameters(par);
                 // improve the pictu
               //   std::cout << "integral error: " << integralerror << std::endl;
                 fitFcnV2_Cvetan->Draw("same");
                 // draw the legend
                 TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                 legend->SetFillColorAlpha(kWhite, 0.);
                 legend->SetBorderSize(0);
                   legend->SetTextFont(42);
                 legend->SetTextSize(0.03);
                   Char_t message[80];
                   sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_Cvetan->GetChisquare(),fitFcnV2_Cvetan->GetNDF());
                   legend->AddEntry(fitFcnV2_Cvetan,message);
            if(res->CovMatrixStatus() == 3){
                     //  sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                       legend->AddEntry(fitFcnV2_Cvetan,message);
                   }
                 legend->AddEntry(YieldTkl_Central_Cvetan,"Data","lpe");
                 legend->Draw();
            
            V2CvetanTKL = par[2];
            errV2CvetanTKL = fitFcnV2_Cvetan->GetParError(2);
            
            FCvetanTKL = par[3];
            errFCvetanTKL = fitFcnV2_Cvetan->GetParError(3);
            
            v2CvetanTKL = sqrt(V2CvetanTKL);
            errv2CvetanTKL = 0.5*(errV2CvetanTKL/V2CvetanTKL)*v2CvetanTKL;
            
            cout << "===== Cvetan Fit =====" <<endl;
            cout << "V2 TKL +/- err V2 TKL: " << V2CvetanTKL << " +/- " << errV2CvetanTKL <<endl;
            cout << "F +/- err F: " << FCvetanTKL << " +/- " << errFCvetanTKL <<endl;
            cout << "v2 TKL +/- err v2 TKL: " << v2CvetanTKL << " +/- " << errv2CvetanTKL <<endl;
            cout << "======================" <<endl;

        }
    
    sprintf(CanvasName,"%s/Fit_Cvetan.pdf",FolderName);
    cTKLCvetan->SaveAs(CanvasName);
    
    
    // Fit Cvetan putting my constraints (F=1, V1 free, V0 not 1)
        
        TCanvas* cTKLCvetanMe = new TCanvas;
           cTKLCvetanMe->SetTitle("TKL Cvetan fit - My constraints F=1, v1 free, v0 free");
           cTKLCvetanMe->Divide(1,1);
           cTKLCvetanMe->cd(1);
           YieldTkl_Central->DrawCopy();
    //       cCvetanMe->cd(2);
    //       Baseline_Central_1->DrawCopy();
    //       cCvetanMe->cd(3);
    //       Yields_Central_1_MinusBaseline->DrawCopy();
        
        
        {
            cTKLCvetanMe->cd(1);
                 TF1 *fitFcnV2_CvetanMe = new TF1("fitFcnV2_CvetanMe",CvetanFTKL,-TMath::Pi()/2,1.5*TMath::Pi(),4);
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
                  TVirtualFitter::Fitter(YieldTkl_Central_CvetanMe)->SetMaxIterations(10000);
                  TVirtualFitter::Fitter(YieldTkl_Central_CvetanMe)->SetPrecision();
                TVirtualFitter::Fitter(YieldTkl_Central_CvetanMe)->SetFCN(ChisquareCvetanFTKL);
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
                  

                 TFitResultPtr res = YieldTkl_Central_CvetanMe->Fit("fitFcnV2_CvetanMe","USBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)YieldTkl_Central_CvetanMe->FindObject("stats");
            st->SetX1NDC(0.8); //new x start position
            st->SetY1NDC(0.8); //new x end position
            
            double chi2, edm, errdef;
            int nvpar, nparx;
             TVirtualFitter::Fitter(YieldTkl_Central_CvetanMe)->GetStats(chi2,edm,errdef,nvpar,nparx);
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
                 legend->SetFillColorAlpha(kWhite, 0.);
                 legend->SetBorderSize(0);
                   legend->SetTextFont(42);
                 legend->SetTextSize(0.03);
                   Char_t message[80];
                   sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_CvetanMe->GetChisquare(),fitFcnV2_CvetanMe->GetNDF());
                   legend->AddEntry(fitFcnV2_CvetanMe,message);
            if(res->CovMatrixStatus() == 3){
                     //  sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                       legend->AddEntry(fitFcnV2_CvetanMe,message);
                   }
                 legend->AddEntry(YieldTkl_Central_CvetanMe,"Data","lpe");
                 legend->Draw();

            V2CvetanMeTKL = par[2];
            errV2CvetanMeTKL = fitFcnV2_CvetanMe->GetParError(2);
            
            V2CvetanMeTKL = par[2]/par[0];
            errV2CvetanMeTKL = V2CvetanMeTKL*sqrt(pow(fitFcnV2_CvetanMe->GetParError(2)/par[2],2)+pow(fitFcnV2_CvetanMe->GetParError(0)/par[0],2));
            
            v2CvetanMeTKL = sqrt(V2CvetanMeTKL);
            errv2CvetanMeTKL = 0.5*(errV2CvetanMeTKL/V2CvetanMeTKL)*v2CvetanMeTKL;
            
            
            cout << "===== CvetanMe Fit =====" <<endl;
            cout << "V2 TKL +/- err V2 TKL: " << V2CvetanMeTKL << " +/- " << errV2CvetanMeTKL <<endl;
            cout << "v2 TKL +/- err v2 TKL: " << v2CvetanMeTKL << " +/- " << errv2CvetanMeTKL <<endl;
            cout << "======================" <<endl;

        }
    
    sprintf(CanvasName,"%s/Fit_CvetanMe.pdf",FolderName);
    cTKLCvetanMe->SaveAs(CanvasName);
    
    // Fit ZYAM Yc = Bc + (Yp-Bp) + a0  + 2a1cos + 2a2cos
       
       
        TCanvas* cZYAMTKL = new TCanvas;
           cZYAMTKL->SetTitle("TKL ZYAM fit");
           cZYAMTKL->Divide(1,1);
           cZYAMTKL->cd(1);
           YieldTkl_Central->DrawCopy();
       //    cZYAM->cd(2);
       //    Baseline_Central_1->DrawCopy();
       //    cZYAM->cd(3);
       //    Yields_Central_1_MinusBaseline->DrawCopy();
       
       
       {
           cZYAMTKL->cd(1);
                TF1 *fitFcnV2_ZYAM = new TF1("fitFcnV2_ZYAM",ZYAMTKL,-TMath::Pi()/2,1.5*TMath::Pi(),3);
                fitFcnV2_ZYAM->SetNpx(500);
                fitFcnV2_ZYAM->SetLineWidth(4);
                fitFcnV2_ZYAM->SetLineColor(kCyan-7);
                // first try without starting values for the parameters
                // This defaults to 1 for each param.
                // this results in an ok fit for the polynomial function
                // however the non-linear part (lorenzian) does not
                // respond well.
                 Double_t params[3] = {1,1,1};
                fitFcnV2_ZYAM->SetParameters(params);
                TVirtualFitter::Fitter(YieldTkl_Central_ZYAM)->SetMaxIterations(10000);
                 TVirtualFitter::Fitter(YieldTkl_Central_ZYAM)->SetPrecision();
                TVirtualFitter::Fitter(YieldTkl_Central_ZYAM)->SetFCN(ChisquareZYAMTKL);
               //  histo->Fit("fitFcn","0");
                // second try: set start values for some parameters

                 fitFcnV2_ZYAM->SetParName(0,"a0");
                   fitFcnV2_ZYAM->SetParName(1,"a1");
                   fitFcnV2_ZYAM->SetParName(2,"a2");
               
                 

                TFitResultPtr res = YieldTkl_Central_ZYAM->Fit("fitFcnV2_ZYAM","USBMERI+","ep");
                gStyle->SetOptStat("n");
               gStyle->SetOptFit(1011);
               TPaveStats *st = (TPaveStats*)YieldTkl_Central_ZYAM->FindObject("stats");
               st->SetX1NDC(0.8); //new x start position
               st->SetY1NDC(0.8); //new x end position
           
                    double chi2, edm, errdef;
                      int nvpar, nparx;
                       TVirtualFitter::Fitter(YieldTkl_Central_ZYAM)->GetStats(chi2,edm,errdef,nvpar,nparx);
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
                legend->SetFillColorAlpha(kWhite, 0.);
                legend->SetBorderSize(0);
                  legend->SetTextFont(42);
                legend->SetTextSize(0.03);
                  Char_t message[80];
                  sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_ZYAM->GetChisquare(),fitFcnV2_ZYAM->GetNDF());
                  legend->AddEntry(fitFcnV2_ZYAM,message);
           if(res->CovMatrixStatus() == 3){
                     // sprintf(message,"The fit is a success");
                  }
                  else{
                      sprintf(message,"The fit is a failure");
                      legend->AddEntry(fitFcnV2_ZYAM,message);
                  }
                legend->AddEntry(YieldTkl_Central_ZYAM,"Data","lpe");
                legend->Draw();

           
           V2ZYAMTKL = par[2]/(par[0]+baselineTKL_central);
           errV2ZYAMTKL = (par[2]/(par[0] + baselineTKL_central)*sqrt(pow(fitFcnV2_ZYAM->GetParError(2)/par[2],2)+pow(fitFcnV2_ZYAM->GetParError(0)/par[0],2)+pow(errbaselineTKL_central/baselineTKL_central,2)));
           
           v2ZYAMTKL = sqrt(V2ZYAMTKL);
           errv2ZYAMTKL = 0.5*(errV2ZYAMTKL/V2ZYAMTKL)*v2ZYAMTKL;
           
           
           cout << "===== ZYAM Fit =====" <<endl;
           cout << "V2 TKL +/- err V2 TKL: " << V2ZYAMTKL << " +/- " << errV2ZYAMTKL <<endl;
           cout << "v2 TKL +/- err v2 TKL: " << v2ZYAMTKL << " +/- " << errv2ZYAMTKL <<endl;
           cout << "====================" <<endl;
           
       }
    
    sprintf(CanvasName,"%s/Fit_ZYAM.pdf",FolderName);
    cZYAMTKL->SaveAs(CanvasName);
    
    // Fit PRL Template
        
        TCanvas* cPRLTemplateTKL = new TCanvas;
           cPRLTemplateTKL->SetTitle("TKL PRL Template fit");
           cPRLTemplateTKL->Divide(1,1);
           cPRLTemplateTKL->cd(1);
           YieldTkl_Central->DrawCopy();
    //       cPRLTemplate->cd(2);
    //       Baseline_Central_1->DrawCopy();
    //       cPRLTemplate->cd(3);
    //       Yields_Central_1_MinusBaseline->DrawCopy();
    
    
//    Double_t x, y, z;
//    int counterid=0;
//      Int_t np = 200*200;
//      TGraph2D *Chi2Evol = new TGraph2D();
//      Chi2Evol->SetTitle("Chi2/ndf as a function of fixed V2 and F - ATLAS fit; V2; F factor; Value of Chi2/ndf");
           
    {
        
        cPRLTemplateTKL->cd(1);
                 TF1 *fitFcnV2_PRLTemplate = new TF1("fitFcnV2_PRLTemplate",PRLTemplateTKL,-TMath::Pi()/2,1.5*TMath::Pi(),2);
                 fitFcnV2_PRLTemplate->SetNpx(500);
                 fitFcnV2_PRLTemplate->SetLineWidth(4);
                 fitFcnV2_PRLTemplate->SetLineColor(kRed+1);
            TF1 *fitFcnV2_PRLTemplate_RidgeAndZero = new TF1("fitFcnV2_PRLTemplate_RidgeAndZero",PRLTemplate_RidgeAndZeroTKL,-TMath::Pi()/2,1.5*TMath::Pi(),2);
            fitFcnV2_PRLTemplate_RidgeAndZero->SetNpx(500);
            fitFcnV2_PRLTemplate_RidgeAndZero->SetLineWidth(2);
            fitFcnV2_PRLTemplate_RidgeAndZero->SetLineStyle(9);
            fitFcnV2_PRLTemplate_RidgeAndZero->SetLineColor(kRed);
            TF1 *fitFcnV2_PRLTemplate_PeriphAndG = new TF1("fitFcnV2_PRLTemplate_PeriphAndG",PRLTemplate_PeriphAndGTKL,-TMath::Pi()/2,1.5*TMath::Pi(),2);
            fitFcnV2_PRLTemplate_PeriphAndG->SetNpx(500);
            fitFcnV2_PRLTemplate_PeriphAndG->SetLineWidth(2);
            fitFcnV2_PRLTemplate_PeriphAndG->SetLineStyle(9);
            fitFcnV2_PRLTemplate_PeriphAndG->SetLineColor(kBlack);
                 // first try without starting values for the parameters
                 // This defaults to 1 for each param.
                 // this results in an ok fit for the polynomial function
                 // however the non-linear part (lorenzian) does not
                 // respond well.
                  Double_t params[2] = {0.005,6};
                 fitFcnV2_PRLTemplate->SetParameters(params);
                  TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate)->SetMaxIterations(10000);
                  TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate)->SetPrecision();
                 TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate)->SetFCN(ChisquarePRLTemplateTKL);
                //  histo->Fit("fitFcn","0");
                 // second try: set start values for some parameters

                  fitFcnV2_PRLTemplate->SetParName(0,"V2");
                    fitFcnV2_PRLTemplate->SetParName(1,"F");
            
          //  fitFcnV2_PRLTemplate->SetParLimits(0,0.0001,5000000);
           // fitFcnV2_PRLTemplate->FixParameter(0,0.005);
          //  fitFcnV2_PRLTemplate->FixParameter(1,1);
         //   fitFcnV2_PRLTemplate->SetParLimits(1,0.1,50);
        //    fitFcnV2_PRLTemplate->SetParLimits(0,0.001,0.008);
                  //  fitFcnV2_PRLTemplate->FixParameter(1,2.5);
           
//    for(double V2_idx=0.0; V2_idx<0.03; V2_idx+=0.00015){ //200
//        for(double F_idx=0.5; F_idx<4.0; F_idx+=0.0175){ //200
//
//            params[0] = V2_idx;
//            params[1] = F_idx;
//            fitFcnV2_PRLTemplate->SetParameters(params);
               
            
//            fitFcnV2_PRLTemplate->FixParameter(0,V2_idx);
//            fitFcnV2_PRLTemplate->FixParameter(1,F_idx);
                     

                    TFitResultPtr res = YieldTkl_Central_PRLTemplate->Fit("fitFcnV2_PRLTemplate","USBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)YieldTkl_Central_PRLTemplate->FindObject("stats");
                    st->SetX1NDC(0.8); //new x start position
                    st->SetY1NDC(0.8); //new x end position
               
               double chi2, edm, errdef;
                int nvpar, nparx;
                 TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate)->GetStats(chi2,edm,errdef,nvpar,nparx);
                 fitFcnV2_PRLTemplate->SetChisquare(chi2);
                 int ndf = npfits-nvpar;
                 fitFcnV2_PRLTemplate->SetNDF(ndf);
               
                   Double_t par[2];
                   fitFcnV2_PRLTemplate->GetParameters(par);
               fitFcnV2_PRLTemplate_RidgeAndZero->SetParameters(par);
               fitFcnV2_PRLTemplate_PeriphAndG->SetParameters(par);
               
               V2PRLTKL = par[0];
               errV2PRLTKL = fitFcnV2_PRLTemplate->GetParError(0);
               
               FPRLTKL = par[1];
               errFPRLTKL = fitFcnV2_PRLTemplate->GetParError(1);
               
               v2PRLTKL = sqrt(V2PRLTKL);
               errv2PRLTKL = 0.5*(errV2PRLTKL/V2PRLTKL)*v2PRLTKL;
            
//
//            x = V2PRLTKL;
//               y = FPRLTKL;
//               z = chi2/ndf;
//               Chi2Evol->SetPoint(counterid,x,y,z);
            
               
               
               cout << "===== PRL Fit =====" <<endl;
               cout << "V2 TKL +/- err V2 TKL: " << V2PRLTKL << " +/- " << errV2PRLTKL <<endl;
               cout << "F +/- err F: " << FPRLTKL << " +/- " << errFPRLTKL <<endl;
               cout << "v2 TKL +/- err v2 TKL: " << v2PRLTKL << " +/- " << errv2PRLTKL <<endl;
               cout << "===================" <<endl;
//            counterid++;
//           }
//    }
        
        // improve the pictu
                        //   std::cout << "integral error: " << integralerror << std::endl;
                          fitFcnV2_PRLTemplate->Draw("same");
                     fitFcnV2_PRLTemplate_RidgeAndZero->Draw("same");
                     fitFcnV2_PRLTemplate_PeriphAndG->Draw("same");
                          // draw the legend
                          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                          legend->SetFillColorAlpha(kWhite, 0.);
                          legend->SetBorderSize(0);
                            legend->SetTextFont(42);
                          legend->SetTextSize(0.03);
                            Char_t message[80];
                            sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate->GetChisquare(),fitFcnV2_PRLTemplate->GetNDF());
                            legend->AddEntry(fitFcnV2_PRLTemplate,message);
                     if(res->CovMatrixStatus() == 3){
//                                sprintf(message,"The fit is a success");
                            }
                            else{
                                sprintf(message,"The fit is a failure");
                                legend->AddEntry(fitFcnV2_PRLTemplate,message);
                            }
                          legend->AddEntry(YieldTkl_Central_PRLTemplate,"Data","lpe");
                            legend->AddEntry(fitFcnV2_PRLTemplate_RidgeAndZero,"F*Baseline_{Periph} + Ridge*G", "L");
                            legend->AddEntry(fitFcnV2_PRLTemplate_PeriphAndG,"F*Yield_{Periph} + G", "L");
                          legend->Draw();
        
    sprintf(CanvasName,"%s/Fit_PRLTemplate.pdf",FolderName);
    cPRLTemplateTKL->SaveAs(CanvasName);
    
//    TCanvas* cChi2Evol = new TCanvas;
//    cChi2Evol->SetTitle("Chi2/ndf evolution");
//    Chi2Evol->Draw("surf1");
    }
    
    // Fit PRL Template + ZYAM Periph
    
    TCanvas* cdiffPRL = new TCanvas;
    cdiffPRL->SetTitle("Difference PRL Methods");
    cdiffPRL->Divide(1,1);
           
           TCanvas* cPRLTemplate_PeriphZYAMTKL = new TCanvas;
              cPRLTemplate_PeriphZYAMTKL->SetTitle("TKL PRL Template fit Periph ZYAM");
              cPRLTemplate_PeriphZYAMTKL->Divide(1,1);
              cPRLTemplate_PeriphZYAMTKL->cd(1);
              YieldTkl_Central->DrawCopy();
    //          cPRLTemplate_PeriphZYAM->cd(2);
    //          Baseline_Central_1->DrawCopy();
    //          cPRLTemplate_PeriphZYAM->cd(3);
    //          Yields_Central_1_MinusBaseline->DrawCopy();
              
              
              {
                  cPRLTemplate_PeriphZYAMTKL->cd(1);
                       TF1 *fitFcnV2_PRLTemplate_PeriphZYAM = new TF1("fitFcnV2_PRLTemplate_PeriphZYAM",PRLTemplate_PeriphZYAMTKL,-TMath::Pi()/2,1.5*TMath::Pi(),2);
                       fitFcnV2_PRLTemplate_PeriphZYAM->SetNpx(500);
                       fitFcnV2_PRLTemplate_PeriphZYAM->SetLineWidth(4);
                       fitFcnV2_PRLTemplate_PeriphZYAM->SetLineColor(kBlack);
                      TF1 *fitFcnV2_PRLTemplate_RidgeAndZero = new TF1("fitFcnV2_PRLTemplate_RidgeAndZero",PRLTemplate_PeriphZYAM_RidgeTKL,-TMath::Pi()/2,1.5*TMath::Pi(),2);
                      fitFcnV2_PRLTemplate_RidgeAndZero->SetNpx(500);
                      fitFcnV2_PRLTemplate_RidgeAndZero->SetLineWidth(2);
                      fitFcnV2_PRLTemplate_RidgeAndZero->SetLineStyle(9);
                      fitFcnV2_PRLTemplate_RidgeAndZero->SetLineColor(kRed);
                      TF1 *fitFcnV2_PRLTemplate_PeriphAndG = new TF1("fitFcnV2_PRLTemplate_PeriphAndG",PRLTemplate_PeriphZYAM_PeriphAndGTKL,-TMath::Pi()/2,1.5*TMath::Pi(),2);
                      fitFcnV2_PRLTemplate_PeriphAndG->SetNpx(500);
                      fitFcnV2_PRLTemplate_PeriphAndG->SetLineWidth(2);
                      fitFcnV2_PRLTemplate_PeriphAndG->SetLineStyle(9);
                      fitFcnV2_PRLTemplate_PeriphAndG->SetLineColor(kBlack);
                       // first try without starting values for the parameters
                       // This defaults to 1 for each param.
                       // this results in an ok fit for the polynomial function
                       // however the non-linear part (lorenzian) does not
                       // respond well.
                        Double_t params[2] = {1,1};
                       fitFcnV2_PRLTemplate_PeriphZYAM->SetParameters(params);
                        TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate_PeriphZYAM)->SetMaxIterations(10000);
                        TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate_PeriphZYAM)->SetPrecision();
                        TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate_PeriphZYAM)->SetFCN(ChisquarePRLTemplate_PeriphZYAMTKL);
                      //  histo->Fit("fitFcn","0");
                       // second try: set start values for some parameters

                        fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(0,"V2");
                          fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(1,"F");
                  
                //  fitFcnV2_PRLTemplate->SetParLimits(0,0.0001,5000000);
               //   fitFcnV2_PRLTemplate->FixParameter(0,0.002);
                  fitFcnV2_PRLTemplate_PeriphZYAM->SetParLimits(1,0.1,5000);
                        //  fitFcnV2_PRLTemplate->FixParameter(1,2.5);
                        

                       TFitResultPtr res = YieldTkl_Central_PRLTemplate_PeriphZYAM->Fit("fitFcnV2_PRLTemplate_PeriphZYAM","USBMERI+","ep");
                  gStyle->SetOptStat("n");
                  gStyle->SetOptFit(1011);
                  TPaveStats *st = (TPaveStats*)YieldTkl_Central_PRLTemplate_PeriphZYAM->FindObject("stats");
                  st->SetX1NDC(0.8); //new x start position
                  st->SetY1NDC(0.8); //new x end position
                  
                  double chi2, edm, errdef;
                  int nvpar, nparx;
                   TVirtualFitter::Fitter(YieldTkl_Central_PRLTemplate_PeriphZYAM)->GetStats(chi2,edm,errdef,nvpar,nparx);
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetChisquare(chi2);
                   int ndf = npfits-nvpar;
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetNDF(ndf);
                  
                      Double_t par[2];
                      fitFcnV2_PRLTemplate_PeriphZYAM->GetParameters(par);
                      fitFcnV2_PRLTemplate_RidgeAndZero->SetParameters(par);
                      fitFcnV2_PRLTemplate_PeriphAndG->SetParameters(par);
                       // improve the pictu
                     //   std::cout << "integral error: " << integralerror << std::endl;
                       fitFcnV2_PRLTemplate_PeriphZYAM->Draw("same");
                      fitFcnV2_PRLTemplate_RidgeAndZero->Draw("same");
                      fitFcnV2_PRLTemplate_PeriphAndG->Draw("same");
                       // draw the legend
                       TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                       legend->SetFillColorAlpha(kWhite, 0.);
                       legend->SetBorderSize(0);
                         legend->SetTextFont(42);
                       legend->SetTextSize(0.03);
                         Char_t message[80];
                         sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate_PeriphZYAM->GetChisquare(),fitFcnV2_PRLTemplate_PeriphZYAM->GetNDF());
                         legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                  if(res->CovMatrixStatus() == 3){
                  //          sprintf(message,"The fit is a success");
                         }
                         else{
                             sprintf(message,"The fit is a failure");
                             legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                         }
                       legend->AddEntry(YieldTkl_Central_PRLTemplate_PeriphZYAM,"Data","lpe");
                        legend->AddEntry(fitFcnV2_PRLTemplate_RidgeAndZero,"Ridge*G", "L");
                        legend->AddEntry(fitFcnV2_PRLTemplate_PeriphAndG,"F*(Yield_{Periph}-Baseline_{Periph}) + G", "L");
                       legend->Draw();
                  
                  V2PRL_PeriphZYAMTKL = par[0];
                  errV2PRL_PeriphZYAMTKL = fitFcnV2_PRLTemplate_PeriphZYAM->GetParError(0);
                  
                  FPRL_PeriphZYAMTKL = par[1];
                  errFPRL_PeriphZYAMTKL = fitFcnV2_PRLTemplate_PeriphZYAM->GetParError(1);
                  
                  v2PRL_PeriphZYAMTKL = sqrt(V2PRL_PeriphZYAMTKL);
                  errv2PRL_PeriphZYAMTKL = 0.5*(errV2PRL_PeriphZYAMTKL/V2PRL_PeriphZYAMTKL)*v2PRL_PeriphZYAMTKL;
                  
                  
                  cout << "===== PRL ZYAM Fit =====" <<endl;
                  cout << "V2 TKL +/- err V2 TKL: " << V2PRL_PeriphZYAMTKL << " +/- " << errV2PRL_PeriphZYAMTKL <<endl;
                  cout << "F +/- err F: " << FPRL_PeriphZYAMTKL << " +/- " << errFPRL_PeriphZYAMTKL <<endl;
                  cout << "v2 TKL +/- err v2 TKL: " << v2PRL_PeriphZYAMTKL << " +/- " << errv2PRL_PeriphZYAMTKL <<endl;
                  cout << "===================" <<endl;
                  

                  cdiffPRL->cd(1);
                  
                  TF1 *diffPRL = new TF1("diffPRL",DiffPRL,-TMath::Pi()/2,1.5*TMath::Pi(),4);
                  diffPRL->SetNpx(500);
                  diffPRL->SetLineWidth(4);
                  diffPRL->SetLineColor(kBlack);
                 // Double_t paramoss[4] = {V2PRLTKL,FPRLTKL,V2PRL_PeriphZYAMTKL,FPRL_PeriphZYAMTKL};
                  Double_t paramoss[4] = {V2PRLTKL,FPRLTKL,V2PRLTKL,FPRLTKL};
                  diffPRL->SetParameters(paramoss);
                  diffPRL->Draw();
              }
    
    sprintf(CanvasName,"%s/Fit_PRLTamplate_ZYAM.pdf",FolderName);
    cPRLTemplate_PeriphZYAMTKL->SaveAs(CanvasName);
    
    TCanvas* cV2_TKL_Methods = new TCanvas;
    const char *methods[6] = {"p-Pb like", "p-Pb like - no ZYAM" , "Central and Periph ZYAM", "Cvetan-Quentin (Template + ZYAM)", "ATLAS Template", "ATLAS Template and ZYAM"};
    const char *methodsF[3] = {"Cvetan-Quentin", "ATLAS Template", "ATLAS Template and ZYAM"};
    double V2Valeurs[6] = {V2ClassiqueTKL, V2ClassiqueTKL_noZYAM, V2ZYAMTKL, V2CvetanTKL, V2PRLTKL, V2PRL_PeriphZYAMTKL};
    double V2Erreurs[6] = {errV2ClassiqueTKL, errV2ClassiqueTKL_noZYAM, errV2ZYAMTKL, errV2CvetanTKL, errV2PRLTKL, errV2PRL_PeriphZYAMTKL};
    double v2Valeurs[6] = {v2ClassiqueTKL, v2ClassiqueTKL_noZYAM, v2ZYAMTKL, v2CvetanTKL, v2PRLTKL, v2PRL_PeriphZYAMTKL};
    double v2Erreurs[6] = {errv2ClassiqueTKL, errv2ClassiqueTKL_noZYAM, errv2ZYAMTKL, errv2CvetanTKL, errv2PRLTKL, errv2PRL_PeriphZYAMTKL};
    double FValeurs[3] = {FCvetanTKL, FPRLTKL, FPRL_PeriphZYAMTKL};
    double FErreurs[3] = {errFCvetanTKL, errFPRLTKL, errFPRL_PeriphZYAMTKL};
    double accroche[6] = {0.5,1.5,2.5,3.5,4.5,5.5};
    double erraccroche[6] = {0.5,0.5,0.5,0.5,0.5,0.5};
    double accrocheF[3] = {0.5,1.5,2.5};
    double erraccrocheF[3] = {0.5,0.5,0.5};
        
    TGraphErrors *grV2_TKL_Methods = new TGraphErrors(6,accroche,V2Valeurs,erraccroche,V2Erreurs);
     // TGraph *gr3 = new TGraph (n, K3, chi);
      grV2_TKL_Methods->SetTitle("V_{2,tkl-tkl} for different extraction methods");
      grV2_TKL_Methods->GetXaxis()->SetTitle("Method");
      grV2_TKL_Methods->GetYaxis()->SetTitle("V_{2,tkl-tkl}");
    for (int i=1;i<=6;i++) grV2_TKL_Methods->GetXaxis()->SetBinLabel(15*i-4,methods[i-1]);
    grV2_TKL_Methods->GetXaxis()->LabelsOption("u");
    grV2_TKL_Methods->GetXaxis()->SetLabelSize(0.025);
   grV2_TKL_Methods->GetYaxis()->SetRangeUser( -0.005, 0.015);
    grV2_TKL_Methods->SetMarkerColor(1);
    grV2_TKL_Methods->SetLineColor(1);
    grV2_TKL_Methods->SetMarkerStyle(4);
      grV2_TKL_Methods->Draw("AP");
    TLine *lv=new TLine(0.0,0.0,6.,0.0);
    lv->SetLineColor(kBlack);
    lv->SetLineWidth(1);
    lv->SetLineStyle(9);
    lv->Draw();
    
    TPaveText *legend1=new TPaveText(0.1,0.1,0.2215,0.9, "brNDC");
    legend1->SetFillColorAlpha(kAzure-3, 0.4);
    legend1->SetFillStyle(3354);
    legend1->SetBorderSize(0);
    legend1->Draw();
    TPaveText *legend2=new TPaveText(0.2215,0.1,0.343,0.9, "brNDC");
    legend2->SetFillColorAlpha(kAzure-3, 0.4);
    legend2->SetBorderSize(0);
    legend2->Draw();
    TPaveText *legend3=new TPaveText(0.343,0.1,0.4645,0.9, "brNDC");
    legend3->SetFillColorAlpha(kCyan-7, 0.4);
    legend3->SetFillStyle(3354);
    legend3->SetBorderSize(0);
    legend3->Draw();
    TPaveText *legend4=new TPaveText(0.4645,0.1,0.586,0.9, "brNDC");
    legend4->SetFillColorAlpha(kOrange+3, 0.4);
    legend4->SetFillStyle(3354);
    legend4->SetBorderSize(0);
    legend4->Draw();
    TPaveText *legend5=new TPaveText(0.586,0.1,0.7075,0.9, "brNDC");
    legend5->SetFillColorAlpha(kRed+1, 0.4);
    legend5->SetBorderSize(0);
    legend5->Draw();
    TPaveText *legend6=new TPaveText(0.7075,0.1,0.829,0.9, "brNDC");
    legend6->SetFillColorAlpha(kBlack, 0.4);
    legend6->SetFillStyle(3354);
    legend6->SetBorderSize(0);
    legend6->Draw();
    
    
    sprintf(CanvasName,"%s/V2_TKL_Methods.pdf",FolderName);
    cV2_TKL_Methods->SaveAs(CanvasName);
    
    TCanvas* cF_TKL_Methods = new TCanvas;
    TGraphErrors *grF_TKL_Methods = new TGraphErrors(3,accrocheF,FValeurs,erraccrocheF,FErreurs);
      // TGraph *gr3 = new TGraph (n, K3, chi);
       grF_TKL_Methods->SetTitle("F for different extraction methods");
       grF_TKL_Methods->GetXaxis()->SetTitle("Method");
       grF_TKL_Methods->GetYaxis()->SetTitle("F_{tkl-tkl}");
     for (int i=1;i<=3;i++) grF_TKL_Methods->GetXaxis()->SetBinLabel(30*i-10,methodsF[i-1]);
     grF_TKL_Methods->GetXaxis()->LabelsOption("u");
    grF_TKL_Methods->GetYaxis()->SetRangeUser( 0.8, 3.);
     grF_TKL_Methods->SetMarkerColor(1);
     grF_TKL_Methods->SetLineColor(1);
     grF_TKL_Methods->SetMarkerStyle(4);
       grF_TKL_Methods->Draw("AP");
    TLine *lF=new TLine(0.0,1.,3.,1.);
       lF->SetLineColor(kBlack);
       lF->SetLineWidth(1);
       lF->SetLineStyle(9);
       lF->Draw();
    
    TPaveText *legendF4=new TPaveText(0.1,0.1,0.343,0.9, "brNDC");
    legendF4->SetFillColorAlpha(kOrange+3, 0.4);
    legendF4->SetFillStyle(3354);
    legendF4->SetBorderSize(0);
    legendF4->Draw();
    TPaveText *legendF5=new TPaveText(0.343,0.1,0.586,0.9, "brNDC");
    legendF5->SetFillColorAlpha(kRed+1, 0.4);
    legendF5->SetBorderSize(0);
    legendF5->Draw();
    TPaveText *legendF6=new TPaveText(0.586,0.1,0.829,0.9, "brNDC");
    legendF6->SetFillColorAlpha(kBlack, 0.4);
    legendF6->SetFillStyle(3354);
    legendF6->SetBorderSize(0);
    legendF6->Draw();
    
    sprintf(CanvasName,"%s/F_TKL_Methods.pdf",FolderName);
    cF_TKL_Methods->SaveAs(CanvasName);
    
    
    TCanvas* cv2_TKL_Methods = new TCanvas;
    TGraphErrors *grv2_TKL_Methods = new TGraphErrors(6,accroche,v2Valeurs,erraccroche,v2Erreurs);
      // TGraph *gr3 = new TGraph (n, K3, chi);
       grv2_TKL_Methods->SetTitle("v_{2,tkl-tkl} for different extraction methods");
       grv2_TKL_Methods->GetXaxis()->SetTitle("Method");
       grv2_TKL_Methods->GetYaxis()->SetTitle("v_{2,tkl-tkl}");
     for (int i=1;i<=6;i++) grv2_TKL_Methods->GetXaxis()->SetBinLabel(15*i-4,methods[i-1]);
     grv2_TKL_Methods->GetXaxis()->LabelsOption("u");
    grv2_TKL_Methods->GetXaxis()->SetLabelSize(0.025);
    grv2_TKL_Methods->GetYaxis()->SetRangeUser( -0.05, 0.15);
     grv2_TKL_Methods->SetMarkerColor(1);
     grv2_TKL_Methods->SetLineColor(1);
     grv2_TKL_Methods->SetMarkerStyle(4);
       grv2_TKL_Methods->Draw("AP");
    lv->Draw();
    legend1->Draw();
    legend2->Draw();
    legend3->Draw();
    legend4->Draw();
    legend5->Draw();
    legend6->Draw();
    
    sprintf(CanvasName,"%s/Smallv2_TKL_Methods.pdf",FolderName);
    cv2_TKL_Methods->SaveAs(CanvasName);
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


Double_t CvetanFTKL(Double_t *x,Double_t *par)

{   int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
  //  double YMinusaprim0 = YieldTkl_Periph->GetBinContent(bintolook+1) - fourierPeriph[0];
    double YMinusBp = YieldTkl_Periph_MinusBaseline->GetBinContent(bintolook+1);
    return baselineTKL_central*(par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + par[3]*YMinusBp; }
   // return fourierCentral[0]*(par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + par[3]*YMinusaprim0; }

void ChisquareCvetanFTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = YieldTkl_Central->GetXaxis();
    
    for(int bin_idx=0; bin_idx<YieldTkl_Central->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (YieldTkl_Central->GetBinContent(bin_idx+1)-CvetanFTKL(x,par))/(sqrt(pow(YieldTkl_Central->GetBinError(bin_idx+1),2)+pow(par[3]*YieldTkl_Periph_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t ZYAMTKL(Double_t *x,Double_t *par)

{   int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
    double YMinusBp = YieldTkl_Periph_MinusBaseline->GetBinContent(bintolook+1);
    return baselineTKL_central + (par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + 1*YMinusBp; }

void ChisquareZYAMTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = YieldTkl_Central->GetXaxis();
    
    for(int bin_idx=0; bin_idx<YieldTkl_Central->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (YieldTkl_Central->GetBinContent(bin_idx+1)-ZYAMTKL(x,par))/(sqrt(pow(YieldTkl_Central->GetBinError(bin_idx+1),2)+pow(errbaselineTKL_central,2)+pow(YieldTkl_Periph_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplateTKL(Double_t *x,Double_t *par)

{   double integral_Yperiph = YieldTkl_Periph->Integral(1,YieldTkl_Periph->GetNbinsX()+1,"width");
    double integral_Yreal = YieldTkl_Central->Integral(1,YieldTkl_Central->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
    double Yperiphbin = YieldTkl_Periph->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(2*TMath::Pi());
    
  //  cout << "G PRL = " << G<<endl;
    
    return par[1]*Yperiphbin + (1+2*par[0]*cos(2*x[0]))*G; }

void ChisquarePRLTemplateTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = YieldTkl_Central->GetXaxis();
    
    for(int bin_idx=0; bin_idx<YieldTkl_Central->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (YieldTkl_Central->GetBinContent(bin_idx+1)-PRLTemplateTKL(x,par))/(sqrt(pow(YieldTkl_Central->GetBinError(bin_idx+1),2)+pow(par[1]*YieldTkl_Periph->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplate_RidgeAndZeroTKL(Double_t *x,Double_t *par)

{   double integral_Yperiph = YieldTkl_Periph->Integral(1,YieldTkl_Periph->GetNbinsX()+1,"width");
    double integral_Yreal = YieldTkl_Central->Integral(1,YieldTkl_Central->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
    double Yperiphbin = baselineTKL_periph;
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[1]*Yperiphbin + (1+2*par[0]*cos(2*x[0]))*G; }

Double_t PRLTemplate_PeriphAndGTKL(Double_t *x,Double_t *par)

{   double integral_Yperiph = YieldTkl_Periph->Integral(1,YieldTkl_Periph->GetNbinsX()+1,"width");
    double integral_Yreal = YieldTkl_Central->Integral(1,YieldTkl_Central->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
    double Yperiphbin = YieldTkl_Periph->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(2*TMath::Pi());
    
    return par[1]*Yperiphbin + G; }

Double_t PRLTemplate_PeriphZYAMTKL(Double_t *x,Double_t *par)

{   double integral_YperiphMinusBp = YieldTkl_Periph_MinusBaseline->Integral(1,YieldTkl_Periph_MinusBaseline->GetNbinsX()+1,"width");
    double integral_Yreal = YieldTkl_Central->Integral(1,YieldTkl_Central->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
    double YMinusBp = YieldTkl_Periph_MinusBaseline->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_YperiphMinusBp))/(2*TMath::Pi());
    
    //cout << "G PRL ZYAM = " << G<<endl;
    return par[1]*YMinusBp + (1+2*par[0]*cos(2*x[0]))*G; }

void ChisquarePRLTemplate_PeriphZYAMTKL(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = YieldTkl_Central->GetXaxis();
    
    for(int bin_idx=0; bin_idx<YieldTkl_Central->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (YieldTkl_Central->GetBinContent(bin_idx+1)-PRLTemplate_PeriphZYAMTKL(x,par))/(sqrt(pow(YieldTkl_Central->GetBinError(bin_idx+1),2)+pow(par[1]*YieldTkl_Periph_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplate_PeriphZYAM_RidgeTKL(Double_t *x,Double_t *par)

{   double integral_YperiphMinusBp = YieldTkl_Periph_MinusBaseline->Integral(1,YieldTkl_Periph_MinusBaseline->GetNbinsX()+1,"width");
    double integral_Yreal = YieldTkl_Central->Integral(1,YieldTkl_Central->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
    double YMinusBp = YieldTkl_Periph_MinusBaseline->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_YperiphMinusBp))/(2*TMath::Pi());
    return (1+2*par[0]*cos(2*x[0]))*G; }

Double_t PRLTemplate_PeriphZYAM_PeriphAndGTKL(Double_t *x,Double_t *par)

{   double integral_YperiphMinusBp = YieldTkl_Periph_MinusBaseline->Integral(1,YieldTkl_Periph_MinusBaseline->GetNbinsX()+1,"width");
    double integral_Yreal = YieldTkl_Central->Integral(1,YieldTkl_Central->GetNbinsX()+1,"width");
    int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *48));
    double YMinusBp = YieldTkl_Periph_MinusBaseline->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_YperiphMinusBp))/(2*TMath::Pi());
    
    return par[1]*YMinusBp + G; }

Double_t DiffPRL(Double_t *x,Double_t *par)

{
   return PRLTemplateTKL(x,par)-PRLTemplate_PeriphZYAMTKL(x,&par[2]);
}


Double_t FourierV2(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);}

Double_t FourierV5(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);} //+ 2*par[3]*cos(3*x[0]) + 2*par[4]*cos(4*x[0]);} + 2*par[5]*cos(5*x[0]) + 2*par[6]*cos(6*x[0]) + 2*par[7]*cos(7*x[0]) + 2*par[8]*cos(8*x[0]) + 2*par[9]*cos(9*x[0]) + 2*par[10]*cos(10*x[0]) + 2*par[11]*cos(11*x[0]);}

