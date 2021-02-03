// Macro qui récupère les histogrammes de CINT et les superpose


#include <fstream>
#include <cstdio>


#include "TMath.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <cmath>


void PlotTogetherCINT(){
    
//    TH2F* YTklCentral(NULL);
//    TH1F* YLOL(NULL);
//    TH2F* YTklNum(NULL);
//    TH2F* YTklDen(NULL);
//
//    YTklNum = new TH2F("YTklNum",
//                      "Yield delta eta wrt delta phi - Tkl YTklNum",
//                      10,0,10,20,0,20);
//    YTklNum->SetXTitle("Correlation DeltaPhi (rad)");
//    YTklNum->SetYTitle("Correlation DeltaEta");
//
//    YTklDen = new TH2F("YTklDen",
//                      "Yield delta eta wrt delta phi - Tkl YTklDen",
//                      10,0,10,20,0,20);
//    YTklDen->SetXTitle("Correlation DeltaPhi (rad)");
//    YTklDen->SetYTitle("Correlation DeltaEta");
//
//    YTklCentral = new TH2F("YTklCentral",
//                      "Yield delta eta wrt delta phi - Tkl Central",
//                      10,0,10,20,0,20);
//    YTklCentral->SetXTitle("Correlation DeltaPhi (rad)");
//    YTklCentral->SetYTitle("Correlation DeltaEta");
//
//    for(int i=0; i<10; i++){
//        for(int j=0; j<20; j++){
//            YTklNum->SetBinContent(i+1,j+1,5*i+j);
//            YTklDen->SetBinContent(i+1,j+1,40-i+j);
//        }
//    }
//
//    YTklCentral->Divide(YTklNum,YTklDen);
//
//    TCanvas *cavt = new TCanvas("cavt","AEavt",0,0,600,600);
//    YTklCentral->DrawCopy("colz");
//
//    YLOL = (TH1F*)YTklCentral->ProjectionX("_px",0,-1,"e");
//     TCanvas *cdiv = new TCanvas("cdiv","AEdiv",0,0,600,600);
//    cdiv->Divide(1,2);
//    cdiv->cd(1);
//    YTklNum->Draw("colz");
//    cdiv->cd(2);
//    YTklDen->Draw("colz");
//
//    TCanvas *clel = new TCanvas("clel","AE2",0,0,600,600);
//    YLOL->DrawCopy("e");
//
//    TCanvas *cTogether = new TCanvas("cTogether","AE",0,0,600,600);
//    YTklCentral->DrawCopy("colz");

    auto legend = new TLegend(0.55,0.4,0.85,0.7);
    legend->SetHeader("DataSets"); // option "C" allows to center the header

    TH1F* hTogether(NULL);
    TH1F* hTampon(NULL);
    Char_t fileInput[200];

    Char_t *Group_Periods[] = {"Group1_CINT","Group2_LHC17h_CINT","Group3_LHC17h_CINT","Group4_CINT","Group5_CINT","Group6_LHC18m_CINT","Group7_LHC18m_CINT","Group8_LHC18m_CINT","Group9_LHC18m_CINT","Group10_LHC18m_CINT","Group11_LHC18m_CINT","Group12_LHC18m_CINT7CENT"};
    int numberOfPeriods = sizeof(Group_Periods) / sizeof(Group_Periods[0]);

    hTogether = new TH1F("hTogether",
                     "Mean Ntkl wrt Zvtx",
                     81,-10,10);
    hTogether->SetXTitle("Zvtx");
    hTogether->SetYTitle("Mean Ntkl");

    hTampon = new TH1F("hTampon",
                     "Mean Ntkl wrt Zvtx",
                     81,-10,10);
    hTampon->SetXTitle("Zvtx");
    hTampon->SetYTitle("Mean Ntkl");

    TCanvas *cTogether = new TCanvas("cTogether","Mean Ntkl wrt z_vtx (CINT7CENT)",0,0,600,600);

   for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
       sprintf(fileInput,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/%s/%s_Nmean.root",Group_Periods[tree_idx],Group_Periods[tree_idx]);
       TFile *fROOT = new TFile(fileInput);


       hTampon = (TH1F*)fROOT->Get("hnseg;1");
       if(tree_idx!=9){
           hTampon->SetLineColor(tree_idx+1);
       }
       else{
           hTampon->SetLineColor(46);
       }
       hTampon->Draw("same");

       legend->AddEntry(hTampon, Group_Periods[tree_idx] ,"lep");
   }

    legend->Draw();

    
}

