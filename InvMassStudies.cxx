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
#include <algorithm>

//#include "TTreeReader.h"
//#include "Event.h"


#include "TAxis.h"
#include "TH1.h"
#include "TArrayD.h"
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

//Code pour essayer différents fits de masse invariante et pulls
// FUNCTIONS

void InvMassStudies();
void MassFitsAndPulls(string SignalF, string BackgroundF, double min, double max, double ratsigma);
Double_t TwoCBE2E(Double_t *x,Double_t *par);
Double_t TwoNA602E(Double_t *x,Double_t *par);
Double_t TwoCBEVWG(Double_t *x,Double_t *par);
Double_t TwoNA60VWG(Double_t *x,Double_t *par);
Double_t TwoCBEPOL1POL2(Double_t *x,Double_t *par);
Double_t TwoNA60POL1POL2(Double_t *x,Double_t *par);
Double_t ExpBkg(Double_t *x,Double_t *par);
Double_t Pol1Pol2(Double_t *x,Double_t *par);
Double_t VWGaussian(Double_t *x,Double_t *par);
Double_t JPsiCrystalBallExtended(Double_t *x,Double_t *par);
Double_t Psi2SCrystalBallExtended(Double_t *x,Double_t *par);
Double_t JPsiNA60(Double_t *x,Double_t *par);
Double_t Psi2SNA60(Double_t *x,Double_t *par);
TFitResultPtr FittingAllInvMassBinStudy(const char *histoname, TCanvas *canvas, int i, string SignalF, string BackgroundF, double min, double max, double ratsigma);
Double_t NumberOfParameters(string SignalF, string BackgroundF);




TH1F* hnseg(NULL);
TH1F* hpool(NULL);
char histoname[50];

Double_t mJpsi =  3.096916;
 Double_t mPsip =  3.686108;
// Double_t ratMass = 1.01;
 //Double_t ratSigma = 1.05;
Int_t npfits;

double MinInvMass = 2.1;
double MaxInvMass = 5.1;

const int NbinsDimuInvMass = 250;

bool PtBinned = kFALSE;
double PtBins[] = {0,2,4,6,8,12};
const int NbPtBins = 5;
double LowDimuPtCut = 0;
double HighDimuPtCut = 12;

Char_t FitFileName[200];
Char_t FolderName[200];



string VarSignals[3] = {"CB2-Run2", "CB2-MC", "NA60-MC"}; //FIXME NA60 breaks VWG ET POL1POL2 IGNOBLESS
string VarBackgrounds[3] = {"VWG", "POL1POL2", "DoubleExpo"};
double VarMinInvMass[2] = {2.3, 2.4};
double VarMaxInvMass[2] = {4.9, 4.7};
double VarRatSigma[2] = {1.0, 1.05};


void InvMassStudies(){
//    for(int signal=0; signal<3; signal++){
//        for(int background=0; background<3; background++){
//            for(int massrange=0; massrange<2; massrange++){
//                for(int ratidx=0;ratidx<2;ratidx++){
//                    MassFitsAndPulls(VarSignals[signal], VarBackgrounds[background], VarMinInvMass[massrange], VarMaxInvMass[massrange], VarRatSigma[ratidx]);
//                }
//            }
//        }
//    }
    
    for(int signal=0; signal<1; signal++){
           for(int background=0; background<3; background++){
                   for(int massrange=0; massrange<2; massrange++){
                                   for(int ratidx=0;ratidx<2;ratidx++){
                                       MassFitsAndPulls(VarSignals[signal], VarBackgrounds[background], VarMinInvMass[massrange], VarMaxInvMass[massrange], VarRatSigma[ratidx]);
                                   }
                               }
           }
       }
    
  //  MassFitsAndPulls(VarSignals[0], VarBackgrounds[2], 2.1, 5.1, 1.05);
    
}

void MassFitsAndPulls(string SignalF, string BackgroundF, double min, double max, double ratsigma){
    
    
    TFitResultPtr res;
    
    TCanvas *cinvmass = new TCanvas("cinvmass","Fitting All Phi - Different Centralities",10,10,700,500);
    cinvmass->SetTitle("Inv Mass Fits");
    // ROOT THAT Inv mass fits
    cinvmass->Divide(1,3);
    
    sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysis_Run2_QGPFrance_0-5_40-100_pt0-2-4-6-8-12.root");
    sprintf(FolderName,"~/Desktop/Images JavierAnalysis/2021 mai/NewAnalysis_16h10_TKL_QGPFrance_0-5_40-100_pt0-12/FitTrainingTKLa");
    
    TFile *filerec = new TFile(FitFileName);
    hnseg = (TH1F*)filerec->Get("hnseg");
//    InvMass_Central = (TH1F*)filerec->Get("InvMass_Central");
//    InvMass_Periph = (TH1F*)filerec->Get("InvMass_Periph");
    
    
//    TCanvas *ctest0 = new TCanvas("ctest0","Test NA60",10,10,700,500);
//    ctest0->SetTitle("Test NA60");
//
//    TF1 *fitFcntest = NULL;
//
//                 fitFcntest = new TF1("fitFcntest",JPsiNA60,min,max,11);
//    double paramstest[11] = {10,mJpsi,0.07,-0.4061,0.2302,1.2048,0.0039,2.0627,0.1836,1.2989,0.00643};
//    fitFcntest->SetNpx(500);
//    fitFcntest->SetLineWidth(4);
//    fitFcntest->SetLineColor(kMagenta);
//    fitFcntest->SetParameters(paramstest);
//
//    hnseg->Draw();
//    fitFcntest->Draw("same");
//    return;
    
    hpool = new TH1F("hpull",
                            "Invariant mass of dimuon - Pulls",
                            NbinsDimuInvMass,MinInvMass,MaxInvMass);
           hpool->SetXTitle("Mass of dimuon (GeV/c^{2})");
           hpool->SetYTitle("Pull (Sigma units)");
    hpool->GetYaxis()->SetRangeUser(-10., 10.);
    hpool->SetStats(kFALSE);
    
    TLine *l0=new TLine(2.1,0,5.1,0);
    TLine *l1=new TLine(2.1,1,5.1,1);
    TLine *lm1=new TLine(2.1,-1,5.1,-1);
    l1->SetLineColor(8);
    l1->SetLineStyle(10);
    lm1->SetLineColor(8);
    lm1->SetLineStyle(10);
    TLine *l3=new TLine(2.1,3,5.1,3);
    TLine *lm3=new TLine(2.1,-3,5.1,-3);
    l3->SetLineColor(46);
    l3->SetLineStyle(10);
    lm3->SetLineColor(46);
    lm3->SetLineStyle(10);

    
    hpool->GetListOfFunctions()->Add(l0);
    hpool->GetListOfFunctions()->Add(l1);
    hpool->GetListOfFunctions()->Add(lm1);
    hpool->GetListOfFunctions()->Add(l3);
    hpool->GetListOfFunctions()->Add(lm3);
    
    int numParameters = NumberOfParameters(SignalF, BackgroundF);
    
    sprintf(histoname,"hnseg");
        res = FittingAllInvMassBinStudy(histoname, cinvmass, 0, SignalF, BackgroundF, min, max, ratsigma);
        double par[numParameters];
        
        // Parameter setting (either from mass fit either start values)
        for(int i=0; i <numParameters; i++){
                par[i] = res->Parameter(i);
        }
    
    //FIXED: Add the pulls in mass inv fit
    
    {
        ofstream myfiletxt;
        myfiletxt.open("/tmp/pullstot.txt");
        double params[numParameters];
        for(int i=0; i<numParameters; i++){
            params[i] = par[i];
        }
        
        TF1 *fitFcn = NULL;
        
        if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
          if(BackgroundF =="DoubleExpo"){
              fitFcn = new TF1("fitFcn",TwoCBE2E,min,max,numParameters);
          }
          else if(BackgroundF =="POL1POL2"){
              fitFcn = new TF1("fitFcn",TwoCBEPOL1POL2,min,max,numParameters);
          }
          else if(BackgroundF =="VWG"){
              fitFcn = new TF1("fitFcn",TwoCBEVWG,min,max,numParameters);
          }
        }
        else if(SignalF=="NA60-MC"){
            if(BackgroundF =="DoubleExpo"){
                fitFcn = new TF1("fitFcn",TwoNA602E,min,max,numParameters);
            }
            else if(BackgroundF =="POL1POL2"){
                fitFcn = new TF1("fitFcn",TwoNA60POL1POL2,min,max,numParameters);
            }
            else if(BackgroundF =="VWG"){
                fitFcn = new TF1("fitFcn",TwoNA60VWG,min,max,numParameters);
            }
        }
        
          fitFcn->SetNpx(500);
          fitFcn->SetLineWidth(4);
          fitFcn->SetLineColor(kMagenta);
          fitFcn->SetParameters(params);

        //Pour chaue bin
            //Trouver la valeur du centre du bin pour le fit et pour les données
            //Calculer le pool
            // Ajouter le point à un histo
        
        myfiletxt << "POUET" <<endl;
         myfiletxt << "POUET " << hnseg->GetNbinsX() <<endl;
        for(int bin=0; bin<hnseg->GetNbinsX(); bin++){
            if(MinInvMass + 3*(0.5+bin)/NbinsDimuInvMass > min && MinInvMass + 3*(0.5+bin)/NbinsDimuInvMass < max){
                double fit_prediction = fitFcn->Eval(MinInvMass + 3*(0.5+bin)/NbinsDimuInvMass);
                double data = hnseg->GetBinContent(bin+1);
                double error = hnseg->GetBinError(bin+1);
                if(error>0){
                    double pool = (data-fit_prediction)/error;

                    hpool->SetBinContent(bin+1,pool);
                    
                   }
            }
        }
        myfiletxt.close();

        TCanvas*cpoolall = new TCanvas();
        cpoolall->Divide(1,1);
        cpoolall->SetTitle("Deviations inv mass fit, all p_{T}");
        hpool->SetTitle("Pulls from invariant mass fit, all p_{T}");
        hpool->DrawCopy();
        
        TCanvas*cmassandpoolall = new TCanvas();
        TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
        pad1->SetBottomMargin(0.00001);
        pad1->SetBorderMode(0);
        pad2->SetTopMargin(0.00001);
        pad2->SetBottomMargin(0.1);
        pad2->SetBorderMode(0);
        pad1->Draw();
        pad2->Draw();
        cmassandpoolall->SetTitle("Invariant mass fit and pulls, all p_{T}");
        pad1->cd();
        TVirtualPad *padvirt = cinvmass->GetPad(1);
        double x1 = padvirt->GetAbsXlowNDC();
        double x2 = x1 + padvirt->GetAbsWNDC();
        double y1 = padvirt->GetAbsYlowNDC();
        double y2 = y1 + padvirt->GetAbsHNDC();
        padvirt->SetPad(0,0.33,1,1);
        padvirt->DrawClone();
        padvirt->SetPad(x1,y1,x2,y2);
        pad2->cd();
        hpool->DrawCopy();
        
    }
    
    TCanvas *cinvmassPtBinned = new TCanvas("cinvmassPtBinned","Fitting All Phi - Different Centralities - Pt Binned",10,10,700,500);
        cinvmassPtBinned->SetTitle("Inv Mass Fits - Pt Binned");
        // ROOT THAT Inv mass fits
        cinvmassPtBinned->Divide(NbPtBins,3);
        
    if(PtBinned){
        TCanvas*cpool = new TCanvas();
        cpool->SetTitle("Deviation inv mass fit");
        cpool->Divide(2,4);
        
        for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    
        sprintf(histoname,Form("bin%d_0",ptbin+1)); // ZZZZZZZ
            res = FittingAllInvMassBinStudy(histoname, cinvmassPtBinned, ptbin, SignalF, BackgroundF, min, max, ratsigma);
        double par[numParameters];
        
        // Parameter setting (either from mass fit either start values)
        for(int i=0; i <numParameters; i++){
                par[i] = res->Parameter(i);
        }
            
        {
            ofstream myfiletxt;
            myfiletxt.open("/tmp/pulls.txt");
            TFile *file0 = new TFile(FitFileName);
            TH1F *histo = (TH1F*)file0->Get(histoname);
            double params[numParameters];
            for(int i=0; i<numParameters; i++){
                params[i] = par[i];
            }
           TF1 *fitFcn = NULL;
            
            if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
              if(BackgroundF =="DoubleExpo"){
                  fitFcn = new TF1("fitFcn",TwoCBE2E,min,max,numParameters);
              }
              else if(BackgroundF =="POL1POL2"){
                  fitFcn = new TF1("fitFcn",TwoCBEPOL1POL2,min,max,numParameters);
              }
              else if(BackgroundF =="VWG"){
                  fitFcn = new TF1("fitFcn",TwoCBEVWG,min,max,numParameters);
              }
            }
            else if(SignalF=="NA60-MC"){
                if(BackgroundF =="DoubleExpo"){
                    fitFcn = new TF1("fitFcn",TwoNA602E,min,max,numParameters);
                }
                else if(BackgroundF =="POL1POL2"){
                    fitFcn = new TF1("fitFcn",TwoNA60POL1POL2,min,max,numParameters);
                }
                else if(BackgroundF =="VWG"){
                    fitFcn = new TF1("fitFcn",TwoNA60VWG,min,max,numParameters);
                }
            }
            
              fitFcn->SetNpx(500);
              fitFcn->SetLineWidth(4);
              fitFcn->SetLineColor(kMagenta);
              fitFcn->SetParameters(params);
            
            //Pour chaue bin
                //Trouver la valeur du centre du bin pour le fit et pour les données
                //Calculer le pool
                // Ajouter le point à un histo

            for(int bin=0; bin<histo->GetNbinsX(); bin++){
                if(MinInvMass + 3*(0.5+bin)/NbinsDimuInvMass > min && MinInvMass + 3*(0.5+bin)/NbinsDimuInvMass < max){
                    double fit_prediction = fitFcn->Eval(MinInvMass + 3*(0.5+bin)/NbinsDimuInvMass);
                    double data = histo->GetBinContent(bin+1);
                    double error = histo->GetBinError(bin+1);
                    double pool = (data-fit_prediction)/error;
                    
                    
                    hpool->SetBinContent(bin+1,pool);
                }
            }
            myfiletxt.close();
            
            cpool->cd(ptbin+1);
            char title[50];
            sprintf(title,Form("Pulls_PtBin_[%.2f,%.2f] GeV",PtBins[ptbin],PtBins[ptbin+1]));
            hpool->SetTitle(title);
            hpool->DrawCopy();
            
            TCanvas*cmassandpoolbinned = new TCanvas();
            TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
            TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
            pad1->SetBottomMargin(0.00001);
            pad1->SetBorderMode(0);
            pad2->SetTopMargin(0.00001);
            pad2->SetBottomMargin(0.1);
            pad2->SetBorderMode(0);
            pad1->Draw();
            pad2->Draw();
            sprintf(title,Form("Mass fit and Pulls PtBin_[%.2f,%.2f] GeV",PtBins[ptbin],PtBins[ptbin+1]));
            cmassandpoolbinned->SetTitle(title);
            pad1->cd();
            TVirtualPad *padvirt = cinvmassPtBinned->GetPad(ptbin+1);
            double x1 = padvirt->GetAbsXlowNDC();
            double x2 = x1 + padvirt->GetAbsWNDC();
            double y1 = padvirt->GetAbsYlowNDC();
            double y2 = y1 + padvirt->GetAbsHNDC();
            padvirt->SetPad(0,0.33,1,1);
            padvirt->DrawClone();
            padvirt->SetPad(x1,y1,x2,y2);
            pad2->cd();
            hpool->DrawCopy();
            
        }
        }
       
    }
    
    
    //Get the inv mass plots (pt, centrality...)
    //Fit with various forms
    //Get statistics and pulls
}


TFitResultPtr FittingAllInvMassBinStudy(const char *histoname, TCanvas *cinvmass, int i, string SignalF, string BackgroundF, double min, double max, double ratsigma){
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
    
    ofstream myfiletxt;
    myfiletxt.open("/tmp/debugging.txt");
    
    int numParameters = NumberOfParameters(SignalF, BackgroundF);
    int nextPar = 0;
    
    cinvmass->cd(i+1);
   cinvmass->SetFillColor(33);
   cinvmass->SetFrameFillColor(41);
   cinvmass->SetGrid();
   TH1F *histo = (TH1F*)file0->Get(histoname);
   // create a TF1 with the range from 0 to 3 and 6 parameters
    
    TF1 *fitFcn = NULL;
    
    if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
      if(BackgroundF =="DoubleExpo"){
          fitFcn = new TF1("fitFcn",TwoCBE2E,min,max,numParameters);
      }
      else if(BackgroundF =="POL1POL2"){
          fitFcn = new TF1("fitFcn",TwoCBEPOL1POL2,min,max,numParameters);
      }
      else if(BackgroundF =="VWG"){
          fitFcn = new TF1("fitFcn",TwoCBEVWG,min,max,numParameters);
      }
    }
    else if(SignalF=="NA60-MC"){
        if(BackgroundF =="DoubleExpo"){
            fitFcn = new TF1("fitFcn",TwoNA602E,min,max,numParameters);
        }
        else if(BackgroundF =="POL1POL2"){
            fitFcn = new TF1("fitFcn",TwoNA60POL1POL2,min,max,numParameters);
        }
        else if(BackgroundF =="VWG"){
            fitFcn = new TF1("fitFcn",TwoNA60VWG,min,max,numParameters);
        }
    }
    
   fitFcn->SetNpx(500);
   fitFcn->SetLineWidth(4);
   fitFcn->SetLineColor(kMagenta);
   // first try without starting values for the parameters
   // This defaults to 1 for each param.
   // this results in an ok fit for the polynomial function
   // however the non-linear part (lorenzian) does not
   // respond well.
    Double_t params[numParameters];
    cout << "numParameters "<<numParameters<<endl;
    std::fill_n(params, numParameters, 1);
    params[0] = 10;
    cout << "params[3] "<<params[3]<<endl;
    cout << "SignalF "<<SignalF<<endl;
    cout << "BackgroundF "<<BackgroundF<<endl;
   fitFcn->SetParameters(params);
    TVirtualFitter::Fitter(histo)->SetMaxIterations(100000);
    TVirtualFitter::Fitter(histo)->SetPrecision();
 //  histo->Fit("fitFcn","0");
   // First fit, fix MJPsi and Sigma
    
     if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
         nextPar = 9;
         
         fitFcn->SetParName(0,"Norm_{J/#psi}");
         fitFcn->SetParName(1,"M_{J/#psi}");
         fitFcn->SetParName(2,"Sigma_{J/#psi}");
         fitFcn->SetParName(3,"a_{1}");
         fitFcn->SetParName(4,"n_{1}");
         fitFcn->SetParName(5,"a_{2}");
         fitFcn->SetParName(6,"n_{2}");
         fitFcn->SetParName(7,"Norm_{#Psi(2S)}");
         fitFcn->SetParName(8,"Sigma_{ratio}");
         
         
           fitFcn->SetParLimits(0,0.1,100000000);
           fitFcn->SetParLimits(1,3.05,3.15);
           fitFcn->SetParLimits(2,0.03,0.1);
           fitFcn->SetParLimits(3,0.7,1.1);
           fitFcn->SetParLimits(4,1,30);
           fitFcn->SetParLimits(5,1,5);
           fitFcn->SetParLimits(6,5,25);
           fitFcn->SetParLimits(7,0.001,1000000);
         fitFcn->SetParLimits(8,0.1,2);
         
         fitFcn->SetParameter(7,100);
         
         fitFcn->FixParameter(1,mJpsi); // Mean x core
         fitFcn->FixParameter(2,0.07);
         fitFcn->FixParameter(8,ratsigma);
         
         if(SignalF=="CB2-Run2"){
              fitFcn->FixParameter(3,0.883);
              fitFcn->FixParameter(4,9.940);
              fitFcn->FixParameter(5,1.832);
              fitFcn->FixParameter(6,15.323);
         }
         else if(SignalF=="CB2-MC"){
              fitFcn->FixParameter(3,0.993);
              fitFcn->FixParameter(4,2.9075);
              fitFcn->FixParameter(5,2.182);
              fitFcn->FixParameter(6,3.122);
         }
         
       }
       if(SignalF=="NA60-MC"){
           
           nextPar = 13;
           
           fitFcn->SetParName(0,"Norm_{J/#psi}");
           fitFcn->SetParName(1,"M_{J/#psi}");
           fitFcn->SetParName(2,"Sigma_{J/#psi}");
           fitFcn->SetParName(3,"a_{L}");
           fitFcn->SetParName(4,"p1_{L}");
           fitFcn->SetParName(5,"p2_{L}");
           fitFcn->SetParName(6,"p3_{L}");
           fitFcn->SetParName(7,"a_{R}");
           fitFcn->SetParName(8,"p1_{R}");
           fitFcn->SetParName(9,"p2_{R}");
           fitFcn->SetParName(10,"p3_{r}");
           fitFcn->SetParName(11,"Norm_{#Psi(2S)}");
           fitFcn->SetParName(12,"Sigma_{ratio}");
           
           fitFcn->SetParLimits(0,0.1,10000000);
           fitFcn->SetParLimits(1,3.05,3.15);
           fitFcn->SetParLimits(2,0.03,0.1);
           fitFcn->SetParLimits(3,-1.,1.);
           fitFcn->SetParLimits(4,0,10);
           fitFcn->SetParLimits(5,0,10);
           fitFcn->SetParLimits(6,0,10);
           fitFcn->SetParLimits(7,0,10);
           fitFcn->SetParLimits(8,0,10);
           fitFcn->SetParLimits(9,0,10);
           fitFcn->SetParLimits(10,0,10);
           fitFcn->SetParLimits(11,0.001,1000000);
           fitFcn->SetParLimits(12,0.1,2);
           
           fitFcn->SetParameter(0,10000);
           fitFcn->SetParameter(11,100);
           
           fitFcn->FixParameter(1,mJpsi); // Mean x core
           fitFcn->FixParameter(2,0.07);
           fitFcn->FixParameter(12,ratsigma);
           
           fitFcn->FixParameter(3,-0.4061);
           fitFcn->FixParameter(4,0.2302);
           fitFcn->FixParameter(5,1.2048);
           fitFcn->FixParameter(6,0.0390);
           fitFcn->FixParameter(7,2.0627);
           fitFcn->FixParameter(8,0.1836);
           fitFcn->FixParameter(9,1.2989);
           fitFcn->FixParameter(10,0.0643);
           
       }
    
    if(BackgroundF =="DoubleExpo"){
        fitFcn->SetParName(nextPar+0,"Norm_{TailLowM}");
        fitFcn->SetParName(nextPar+1,"Exp_{TailLowM}");
        fitFcn->SetParName(nextPar+2,"Norm_{TailHighM}");
        fitFcn->SetParName(nextPar+3,"Exp_{TailHighM}");
        
        fitFcn->SetParLimits(nextPar+0,0.1,10000000);
        fitFcn->SetParLimits(nextPar+1,0.01,20);
        fitFcn->SetParLimits(nextPar+2,0.01,100000000);
        fitFcn->SetParLimits(nextPar+3,0.01,50);
        
        fitFcn->SetParameter(nextPar+0,1000);
        fitFcn->SetParameter(nextPar+1,0.5);
        fitFcn->SetParameter(nextPar+3,10);
    }
       if(BackgroundF =="POL1POL2"){
           
           fitFcn->SetParName(nextPar+0,"Norm_{Bkg}");
           fitFcn->SetParName(nextPar+1,"a_{1}");
           fitFcn->SetParName(nextPar+2,"b_{1}");
           fitFcn->SetParName(nextPar+3,"b_{2}");
           
           fitFcn->SetParLimits(nextPar+0,0.001,10000000);
           fitFcn->SetParLimits(nextPar+1,-100,100);
           fitFcn->SetParLimits(nextPar+2,-100,100);
           fitFcn->SetParLimits(nextPar+3,-100,100);
           
           fitFcn->SetParameter(nextPar+0,1000);
           fitFcn->SetParameter(nextPar+1,0.1);
           fitFcn->SetParameter(nextPar+2,0.1);
           fitFcn->SetParameter(nextPar+3,0.1);
       }
    
    if(BackgroundF =="VWG"){
        
        fitFcn->SetParName(nextPar+0,"Norm_{Bkg}");
        fitFcn->SetParName(nextPar+1,"Mean_{Bkg}");
        fitFcn->SetParName(nextPar+2,"Alpha");
        fitFcn->SetParName(nextPar+3,"Beta");
        
        fitFcn->SetParLimits(nextPar+0,0.001,100000000);
        fitFcn->SetParLimits(nextPar+1,0.001,max);
        fitFcn->SetParLimits(nextPar+2,0.001,100);
        fitFcn->SetParLimits(nextPar+3,0.001,100);
        
        fitFcn->SetParameter(nextPar+0,100);
        fitFcn->SetParameter(nextPar+1,1);
        fitFcn->SetParameter(nextPar+2,1);
        fitFcn->SetParameter(nextPar+3,1);
    }
    
    Double_t par[numParameters];
    
   TFitResultPtr res = histo->Fit("fitFcn","SBMER","ep");
 //  res = histo->Fit("fitFcn","SBMERQ","ep");
      fitFcn->ReleaseParameter(1);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    fitFcn->ReleaseParameter(2);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
        fitFcn->ReleaseParameter(7);
    }
    if(SignalF=="NA60-MC"){
        fitFcn->ReleaseParameter(11);
    }
   res = histo->Fit("fitFcn","SBMERQ","ep");
   res = histo->Fit("fitFcn","SBMER","ep");
    gStyle->SetOptStat("ne");
    gStyle->SetOptFit(1012);
    TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
    st->SetX1NDC(0.8); //new x start position
    st->SetY1NDC(0.3); //new x end position
//    fitFcn->ReleaseParameter(3);
//       res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(5);
//   res = histo->Fit("fitFcn","SBMER","ep");
//   fitFcn->ReleaseParameter(4);
//    res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(6);
//   res = histo->Fit("fitFcn","SBMER","ep");
//    res = histo->Fit("fitFcn","SBMER","ep");

    
    TF1 *backFcn = NULL;
           
             if(BackgroundF =="DoubleExpo"){
                 backFcn = new TF1("backFcn",ExpBkg,min,max,4);
             }
             else if(BackgroundF =="POL1POL2"){
                 backFcn = new TF1("backFcn",Pol1Pol2,min,max,4);
             }
             else if(BackgroundF =="VWG"){
                 backFcn = new TF1("backFcn",VWGaussian,min,max,4);
             }
    
    TF1 *signalFcnJPsi = NULL;
    TF1 *signalFcnPsi2S = NULL;
    
        if(SignalF =="CB2-Run2" || SignalF =="CB2-MC"){
            signalFcnJPsi = new TF1("signalFcnJPsi",JPsiCrystalBallExtended,min,max,numParameters-4);
            signalFcnPsi2S = new TF1("signalFcnPsi2S",Psi2SCrystalBallExtended,min,max,numParameters-4);
        }
        else if(SignalF =="NA60-MC"){
            signalFcnJPsi = new TF1("signalFcnJPsi",JPsiNA60,min,max,numParameters-4);
            signalFcnPsi2S = new TF1("signalFcnPsi2S",Psi2SNA60,min,max,numParameters-4);
        }
        
    
   backFcn->SetLineColor(kRed);
   TPaveText *pave = new TPaveText(0.5,0.5,0.7,0.6,"brNDC");
   signalFcnJPsi->SetLineColor(kBlue);
   signalFcnJPsi->SetNpx(500);
    signalFcnPsi2S->SetLineColor(kGreen);
    signalFcnPsi2S->SetNpx(500);

   fitFcn->GetParameters(par);
   signalFcnJPsi->SetParameters(par);
    signalFcnPsi2S->SetParameters(par);

   Double_t integral = (signalFcnJPsi->Integral(min,3.45))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
    auto covtot = res->GetCovarianceMatrix();
    auto covsgn = covtot.GetSub(0,numParameters-4,0,numParameters-4);
    std::cout << "Matrice totale" <<endl;
    covtot.Print();
    std::cout << "Matrice réduite" <<endl;
    covsgn.Print();
    std::cout << "STATUS COV " << res->CovMatrixStatus() <<endl;
   Double_t integralerror = (signalFcnJPsi->IntegralError(min,3.45,signalFcnJPsi->GetParameters(), res->GetCovarianceMatrix().GetSub(0,numParameters-4,0,numParameters-4).GetMatrixArray() ))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
    std::cout << "Erreur integrale " << integralerror <<endl;
    
    std::cout << "Fitted " << histoname << std::endl;
    std::cout << "Nb JPsi total: " << integral << std::endl;
 //   std::cout << "integral error: " << integralerror << std::endl;
   signalFcnJPsi->Draw("same");
    signalFcnPsi2S->Draw("same");
   backFcn->SetParameters(&par[numParameters-4]);
   backFcn->Draw("same");
   // draw the legend
    char str[50];
    string fitPerformed = SignalF+BackgroundF;
    int n = fitPerformed.length();
    char char_fitPerformed[n + 1];
       strcpy(char_fitPerformed, fitPerformed.c_str());
    
    sprintf(str, "N_{J/#psi} = %i +/- %i  ", int(integral), int(integralerror));
   pave->AddText(str);
    pave->AddText(char_fitPerformed);
    pave->SetTextFont(42);
    pave->SetTextSize(0.04);
    pave->SetBorderSize(0);
    pave->SetFillStyle(0);
    sprintf(str, "M_{#Psi(2S)} = %f, Sig_{#Psi(2S)} = %f", int(mPsip*1000)/1000., int(par[numParameters-4]*par[2]*1000)/1000.);
  //  pave->AddText(str);
   pave->Draw();
   TLegend *legend=new TLegend(0.4,0.6,0.6,0.8);
   legend->SetTextFont(42);
   legend->SetTextSize(0.03);
    legend->SetFillColorAlpha(kWhite, 0.);
    legend->SetBorderSize(0);
    Char_t message[80];
    sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcn->GetChisquare(),fitFcn->GetNDF());
     legend->AddEntry(fitFcn,message);
    if(res->CovMatrixStatus() == 3){
              // sprintf(message,"The fit is a success");
           }
           else{
               sprintf(message,"The fit is a failure");
               legend->AddEntry(fitFcn,message);
           }
   legend->AddEntry(histo,"Data","lpe");
   legend->AddEntry(backFcn,"Background fit","l");
   legend->AddEntry(signalFcnJPsi,"J/#psi Signal fit","l");
    legend->AddEntry(signalFcnPsi2S,"#Psi(2S) Signal fit","l");
   legend->Draw();
    
//    cout << "Histogramme: " << histoname << endl;
//    for(int j=0; j<10; j++){
//        cout << "Bin de masse numero " << j << " - Signal : " << (signalFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 << " et Background : " << (backFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 <<endl;
//    }
     myfiletxt.close();
    return res;
}


//Functions templates

Double_t NumberOfParameters(string SignalF, string BackgroundF){
    
    int numberParameters = 1; //Start at 1 to account for fixed sigmaratio parameter
    
    if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
        numberParameters += 8;
    }
    if(SignalF=="NA60-MC"){
        numberParameters += 12;
    }
    if(BackgroundF =="DoubleExpo" || BackgroundF =="POL1POL2" || BackgroundF =="VWG"){
        numberParameters += 4;
    }
    
    return numberParameters;
    
}

Double_t TwoCBE2E(Double_t *x,Double_t *par)

{ return JPsiCrystalBallExtended(x,par) + Psi2SCrystalBallExtended(x,par) + ExpBkg(x,&par[9]);}

Double_t TwoNA602E(Double_t *x,Double_t *par)

{ return JPsiNA60(x,par) + Psi2SNA60(x,par) + ExpBkg(x,&par[13]);}

Double_t TwoCBEVWG(Double_t *x,Double_t *par)

{ return JPsiCrystalBallExtended(x,par) + Psi2SCrystalBallExtended(x,par) + VWGaussian(x,&par[9]);}

Double_t TwoNA60VWG(Double_t *x,Double_t *par)

{ return JPsiNA60(x,par) + Psi2SNA60(x,par) + VWGaussian(x,&par[13]);}

Double_t TwoCBEPOL1POL2(Double_t *x,Double_t *par)

{ return JPsiCrystalBallExtended(x,par) + Psi2SCrystalBallExtended(x,par) + Pol1Pol2(x,&par[9]);}

Double_t TwoNA60POL1POL2(Double_t *x,Double_t *par)

{ return JPsiNA60(x,par) + Psi2SNA60(x,par) + Pol1Pol2(x,&par[13]);}


Double_t ExpBkg(Double_t *x,Double_t *par)

{ return par[0]*(exp(x[0]*par[1]*(-1))) + par[2]*(exp(x[0]*par[3]*(-1))); }

Double_t Pol1Pol2(Double_t *x,Double_t *par)

{ return par[0]*((1+par[1]*x[0])/(par[2]*x[0] + par[3]*(x[0]*x[0]))); }

Double_t VWGaussian(Double_t *x,Double_t *par)

{ double sigma = par[2] + par[3]*((x[0]-par[1])/par[1]);
    return par[0]*exp(-1.0*pow(x[0]-par[1],2)/(2*pow(sigma,2))); }

//FIXME : Coder NA60

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
    
    Double_t t = (x[0]-(par[1]*mPsip/mJpsi))/(par[2]*par[8]);
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t < -absAlpha) //left tail
  {
      Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
      Double_t b = par[4]/absAlpha - absAlpha;
      
      sum += (par[7]/(par[2]*par[8]*sqrt(2*TMath::Pi())))*(a/TMath::Power(b - t, par[4]));
  }
  else if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
      sum += (par[7]/(par[2]*par[8]*sqrt(2*TMath::Pi())))*(exp(-0.5*t*t));
  }
  else if (t >= absAlpha2) //right tail
  {
      Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
      Double_t d = par[6]/absAlpha2 - absAlpha2;
      
      sum += (par[7]/(par[2]*par[8]*sqrt(2*TMath::Pi())))*(c/TMath::Power(d + t, par[6]));
  } else
      sum += 0;
    
    return sum ;

}

Double_t JPsiNA60(Double_t *x,Double_t *par)
{
    
    Double_t sum = 0;
    
    Double_t t0 = 1;
    
    Double_t t = (x[0]-par[1])/par[2];
    
    if (t < par[3]) //left tail
    {
        t0 = 1+pow(par[4]*(par[3]-t),par[5]-(par[6]*sqrt(par[3]-t)));

        sum += par[0]*exp(-0.5*pow(t/t0,2));
    }
    else if (t > par[7]) //right
    {
        t0 = 1+pow(par[8]*(t-par[7]),par[9]-(par[10]*sqrt(t-par[7])));

        sum += par[0]*exp(-0.5*pow(t/t0,2));
    }
    else{
        t0 = 1;
        
        sum += par[0]*exp(-0.5*pow(t/t0,2));
    }
    
    return sum ;
}

Double_t Psi2SNA60(Double_t *x,Double_t *par)
{
    
    Double_t sum = 0;
    
    Double_t t0 = 1;
    
    Double_t t = (x[0]-(par[1]*mPsip/mJpsi))/(par[2]*par[12]);
    
    if (t < par[3]) //left tail
    {
        t0 = 1+pow(par[4]*(par[3]-t),par[5]-(par[6]*sqrt(par[3]-t)));
        
        sum += par[11]*exp(-0.5*pow(t/t0,2));
    }
    else if (t > par[7]) //right
    {
        t0 = 1+pow(par[8]*(t-par[7]),par[9]-(par[10]*sqrt(t-par[7])));
        
        sum += par[11]*exp(-0.5*pow(t/t0,2));
    }
    else{
        t0 = 1;
        
        sum += par[11]*exp(-0.5*pow(t/t0,2));
    }
    
    return sum ;
}
