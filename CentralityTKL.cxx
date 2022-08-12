// Macro qui plot la centralité tkl


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
 # include <iostream>
# include <sstream>
 # include <string>
# include <iterator>

void CentraliteTKLV0M5();
void CentraliteTKLV0M10();
void CentraliteTKLSPDT5();
void CentraliteTKLSPDT10();


void CentraliteTKLV0M10(){
    // Figure V0M 10mrad
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

        
       // Associate cut 0-100 ATLAS
        
       const Int_t numsituations=4;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.5;
       errsit[i]=0;
       };
        
       Double_t v2pPb[numsituations] = {
       0.0589,//
       0.06,//
       0.062,//
       0.0614//
       };

       Double_t errv2pPb[numsituations] = {
       0.0007,//
       0.0006,//
       0.0011,//
       0.0023//
       };

       std::string names[numsituations] = {
       "0-1",
       "1-3",
       "3-5",
        "5-10"
       };

       TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        v2syst->SetTitle("v_{2,tkl} for different central classes");
       TAxis *ax = v2syst->GetHistogram()->GetXaxis();
       Double_t x1 = ax->GetBinLowEdge(1);
       Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
       v2syst->GetHistogram()->GetXaxis()->Set(4,x1-0.5,x2+0.5);
        v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

       for(Int_t k=0;k<numsituations;k++){
       v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
       }

       v2syst->SetMarkerStyle(21);
        v2syst->SetMarkerColor(kRed);
        v2syst->SetLineColor(kRed);
       v2syst->Draw("AP");
        
        
        
        // Associate cut 20-100 ATLAS
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.45;
        errsit[i]=0;
        };

        Double_t v2pPb2[numsituations] = {
        0.0604,//
        0.0613,//
        0.0625,//
        0.0625//
        };

        Double_t errv2pPb2[numsituations] = {
        0.0005,//
        0.0004,//
        0.0005,//
        0.0007//
        };

        TGraphErrors *v2syst2 = new TGraphErrors(numsituations,sit,v2pPb2,errsit,errv2pPb2);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst2->SetMarkerStyle(21);
         v2syst2->SetMarkerColor(kOrange);
         v2syst2->SetLineColor(kOrange);
        v2syst2->Draw("P");
        
        
        // Associate cut 40-100 ATLAS
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.4;
            errsit[i]=0;
            };

            Double_t v2pPb3[numsituations] = {
            0.0622,//
            0.0633,//
            0.0647,//
            0.0652//
            };

            Double_t errv2pPb3[numsituations] = {
            0.0005,//
            0.0004,//
            0.0005,//
            0.0006//
            };

            TGraphErrors *v2syst3 = new TGraphErrors(numsituations,sit,v2pPb3,errsit,errv2pPb3);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst3->SetMarkerStyle(21);
             v2syst3->SetMarkerColor(kGreen);
             v2syst3->SetLineColor(kGreen);
            v2syst3->Draw("P");
    
    // Associate cut 60-100 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.35;
        errsit[i]=0;
        };

        Double_t v2pPb4[numsituations] = {
        0.0651,//
        0.0664,//
        0.0682,//
        0.0692//
        };

        Double_t errv2pPb4[numsituations] = {
        0.0005,//
        0.0005,//
        0.0007,//
        0.0008//
        };

        TGraphErrors *v2syst4 = new TGraphErrors(numsituations,sit,v2pPb4,errsit,errv2pPb4);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst4->SetMarkerStyle(21);
         v2syst4->SetMarkerColor(kBlue);
         v2syst4->SetLineColor(kBlue);
        v2syst4->Draw("P");
    
    // Associate cut 20-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.3;
        errsit[i]=0;
        };

        Double_t v2pPb5[numsituations] = {
        0.0601,//
        0.0610,//
        0.0622,//
        0.0620//
        };

        Double_t errv2pPb5[numsituations] = {
        0.0005,//
        0.0004,//
        0.0006,//
        0.0007//
        };

        TGraphErrors *v2syst5 = new TGraphErrors(numsituations,sit,v2pPb5,errsit,errv2pPb5);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst5->SetMarkerStyle(25);
         v2syst5->SetMarkerColor(kOrange);
         v2syst5->SetLineColor(kOrange);
        v2syst5->Draw("P");
        
    // Associate cut 40-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.25;
        errsit[i]=0;
        };

        Double_t v2pPb6[numsituations] = {
        0.0616,//
        0.0626,//
        0.0638,//
        0.0641//
        };

        Double_t errv2pPb6[numsituations] = {
        0.0005,//
        0.0005,//
        0.0006,//
        0.0007//
        };

        TGraphErrors *v2syst6 = new TGraphErrors(numsituations,sit,v2pPb6,errsit,errv2pPb6);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst6->SetMarkerStyle(25);
         v2syst6->SetMarkerColor(kGreen);
         v2syst6->SetLineColor(kGreen);
        v2syst6->Draw("P");
    
    // Associate cut 60-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.2;
        errsit[i]=0;
        };

        Double_t v2pPb7[numsituations] = {
        0.0637,//
        0.0648,//
        0.0663,//
        0.0671//
        };

        Double_t errv2pPb7[numsituations] = {
        0.0006,//
        0.0007,//
        0.0008,//
        0.0010//
        };

        TGraphErrors *v2syst7 = new TGraphErrors(numsituations,sit,v2pPb7,errsit,errv2pPb7);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst7->SetMarkerStyle(25);
         v2syst7->SetMarkerColor(kBlue);
         v2syst7->SetLineColor(kBlue);
        v2syst7->Draw("P");
    
    
    // Associate cut 0-100 ALICE

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.15;
       errsit[i]=0;
       };
        
       Double_t v2pPb8[numsituations] = {
       0.0387,//
       0.0363,//
       0.0314,//
       0.0227//
       };

       Double_t errv2pPb8[numsituations] = {
       0.0004,//
       0.0003,//
       0.0004,//
       0.0006//
       };

       TGraphErrors *v2syst8 = new TGraphErrors(numsituations,sit,v2pPb8,errsit,errv2pPb8);
    
       v2syst8->SetMarkerStyle(20);
        v2syst8->SetMarkerColor(kRed);
        v2syst8->SetLineColor(kRed);
       v2syst8->Draw("P");
        
        
        
        // Associate cut 20-100 ALICE
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.1;
        errsit[i]=0;
        };

        Double_t v2pPb9[numsituations] = {
        0.0485,//
        0.0476,//
        0.0458,//
        0.0426//
        };

        Double_t errv2pPb9[numsituations] = {
        0.0003,//
        0.0002,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst9 = new TGraphErrors(numsituations,sit,v2pPb9,errsit,errv2pPb9);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst9->SetMarkerStyle(20);
         v2syst9->SetMarkerColor(kOrange);
         v2syst9->SetLineColor(kOrange);
        v2syst9->Draw("P");
        
        
        // Associate cut 40-100 ALICE
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.05;
            errsit[i]=0;
            };

            Double_t v2pPb10[numsituations] = {
            0.0531,//
            0.0528,//
            0.0520,//
            0.0502//
            };

            Double_t errv2pPb10[numsituations] = {
            0.0003,//
            0.0003,//
            0.0003,//
            0.0004//
            };

            TGraphErrors *v2syst10 = new TGraphErrors(numsituations,sit,v2pPb10,errsit,errv2pPb10);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst10->SetMarkerStyle(20);
             v2syst10->SetMarkerColor(kGreen);
             v2syst10->SetLineColor(kGreen);
            v2syst10->Draw("P");
    
    // Associate cut 60-100 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i;
        errsit[i]=0;
        };

        Double_t v2pPb11[numsituations] = {
        0.0570,//
        0.0571,//
        0.0571,//
        0.0563//
        };

        Double_t errv2pPb11[numsituations] = {
        0.0003,//
        0.0003,//
        0.0004,//
        0.0004//
        };

        TGraphErrors *v2syst11 = new TGraphErrors(numsituations,sit,v2pPb11,errsit,errv2pPb11);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst11->SetMarkerStyle(20);
         v2syst11->SetMarkerColor(kBlue);
         v2syst11->SetLineColor(kBlue);
        v2syst11->Draw("P");
    
    // Associate cut 20-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.05;
        errsit[i]=0;
        };

        Double_t v2pPb12[numsituations] = {
        0.0482,//
        0.0473,//
        0.0453,//
        0.0420//
        };

        Double_t errv2pPb12[numsituations] = {
        0.0003,//
        0.0002,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst12 = new TGraphErrors(numsituations,sit,v2pPb12,errsit,errv2pPb12);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst12->SetMarkerStyle(24);
         v2syst12->SetMarkerColor(kOrange);
         v2syst12->SetLineColor(kOrange);
        v2syst12->Draw("P");
        
    // Associate cut 40-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb13[numsituations] = {
        0.0524,//
        0.0520,//
        0.0511,//
        0.0491//
        };

        Double_t errv2pPb13[numsituations] = {
        0.0003,//
        0.0003,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst13 = new TGraphErrors(numsituations,sit,v2pPb13,errsit,errv2pPb13);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst13->SetMarkerStyle(24);
         v2syst13->SetMarkerColor(kGreen);
         v2syst13->SetLineColor(kGreen);
        v2syst13->Draw("P");
    
    // Associate cut 60-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.15;
        errsit[i]=0;
        };

        Double_t v2pPb14[numsituations] = {
        0.0556,//
        0.0556,//
        0.0553,//
        0.0542//
        };

        Double_t errv2pPb14[numsituations] = {
        0.0004,//
        0.0004,//
        0.0005,//
        0.0015//
        };

        TGraphErrors *v2syst14 = new TGraphErrors(numsituations,sit,v2pPb14,errsit,errv2pPb14);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst14->SetMarkerStyle(24);
         v2syst14->SetMarkerColor(kBlue);
         v2syst14->SetLineColor(kBlue);
        v2syst14->Draw("P");
        
    
        
    
    
        
        TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
        legend5->SetHeader("V0M - 10mrad");
          legend5->SetFillColorAlpha(kWhite, 0.);
          legend5->SetBorderSize(0);
            legend5->SetTextFont(42);
          legend5->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Template fit 0-100");
            legend5->AddEntry(v2syst,message);
       sprintf(message,"Template fit 20-100");
        legend5->AddEntry(v2syst2,message);
        sprintf(message,"Template fit 40-100");
        legend5->AddEntry(v2syst3,message);
    sprintf(message,"Template fit 60-100");
    legend5->AddEntry(v2syst4,message);
    sprintf(message,"Template fit 20-80");
    legend5->AddEntry(v2syst5,message);
    sprintf(message,"Template fit 40-80");
    legend5->AddEntry(v2syst6,message);
    sprintf(message,"Template fit 60-80");
    legend5->AddEntry(v2syst7,message);
    sprintf(message,"Yield sub 0-100");
            legend5->AddEntry(v2syst8,message);
       sprintf(message,"Yield sub 20-100");
        legend5->AddEntry(v2syst9,message);
        sprintf(message,"Yield sub 40-100");
        legend5->AddEntry(v2syst10,message);
    sprintf(message,"Yield sub 60-100");
    legend5->AddEntry(v2syst11,message);
    sprintf(message,"Yield sub 20-80");
    legend5->AddEntry(v2syst12,message);
    sprintf(message,"Yield sub 40-80");
    legend5->AddEntry(v2syst13,message);
    sprintf(message,"Yield sub 60-80");
    legend5->AddEntry(v2syst14,message);
          legend5->Draw();
}


void CentraliteTKLSPDT10(){
    // Figure SPDT 10mrad
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

        
       // Associate cut 0-100 ATLAS
        
       const Int_t numsituations=4;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.5;
       errsit[i]=0;
       };
        
       Double_t v2pPb[numsituations] = {
       0.0752,//
       0.0768,//
       0.0805,//
       0.0901//
       };

       Double_t errv2pPb[numsituations] = {
       0.0009,//
       0.0009,//
       0.0017,//
       0.0045//
       };

       std::string names[numsituations] = {
       "0-1",
       "1-3",
       "3-5",
        "5-10"
       };

       TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        v2syst->SetTitle("v_{2,tkl} for different central classes");
       TAxis *ax = v2syst->GetHistogram()->GetXaxis();
       Double_t x1 = ax->GetBinLowEdge(1);
       Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
       v2syst->GetHistogram()->GetXaxis()->Set(4,x1-0.5,x2+0.5);
        v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

       for(Int_t k=0;k<numsituations;k++){
       v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
       }

       v2syst->SetMarkerStyle(21);
        v2syst->SetMarkerColor(kRed);
        v2syst->SetLineColor(kRed);
       v2syst->Draw("AP");
        
        
        
        // Associate cut 20-100 ATLAS
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.45;
        errsit[i]=0;
        };

        Double_t v2pPb2[numsituations] = {
        0.0820,//
        0.0839,//
        0.0873,//
        0.0909//
        };

        Double_t errv2pPb2[numsituations] = {
        0.0006,//
        0.0006,//
        0.0008,//
        0.0011//
        };

        TGraphErrors *v2syst2 = new TGraphErrors(numsituations,sit,v2pPb2,errsit,errv2pPb2);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst2->SetMarkerStyle(21);
         v2syst2->SetMarkerColor(kOrange);
         v2syst2->SetLineColor(kOrange);
        v2syst2->Draw("P");
        
        
        // Associate cut 40-100 ATLAS
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.4;
            errsit[i]=0;
            };

            Double_t v2pPb3[numsituations] = {
            0.0968,//
            0.0999,//
            0.1057,//
            0.1115//
            };

            Double_t errv2pPb3[numsituations] = {
            0.0009,//
            0.0009,//
            0.0012,//
            0.0014//
            };

            TGraphErrors *v2syst3 = new TGraphErrors(numsituations,sit,v2pPb3,errsit,errv2pPb3);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst3->SetMarkerStyle(21);
             v2syst3->SetMarkerColor(kGreen);
             v2syst3->SetLineColor(kGreen);
            v2syst3->Draw("P");
    
    // Associate cut 60-100 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.35;
        errsit[i]=0;
        };

        Double_t v2pPb4[numsituations] = {
        0.1105,//
        0.1140,//
        0.1203,//
        0.1261//
        };

        Double_t errv2pPb4[numsituations] = {
        0.0014,//
        0.0015,//
        0.0017,//
        0.0020//
        };

        TGraphErrors *v2syst4 = new TGraphErrors(numsituations,sit,v2pPb4,errsit,errv2pPb4);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst4->SetMarkerStyle(21);
         v2syst4->SetMarkerColor(kBlue);
         v2syst4->SetLineColor(kBlue);
        v2syst4->Draw("P");
    
    // Associate cut 20-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.3;
        errsit[i]=0;
        };

        Double_t v2pPb5[numsituations] = {
        0.0815,//
        0.0833,//
        0.0866,//
        0.0901//
        };

        Double_t errv2pPb5[numsituations] = {
        0.0006,//
        0.0006,//
        0.0009,//
        0.0011//
        };

        TGraphErrors *v2syst5 = new TGraphErrors(numsituations,sit,v2pPb5,errsit,errv2pPb5);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst5->SetMarkerStyle(25);
         v2syst5->SetMarkerColor(kOrange);
         v2syst5->SetLineColor(kOrange);
        v2syst5->Draw("P");
        
    // Associate cut 40-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.25;
        errsit[i]=0;
        };

        Double_t v2pPb6[numsituations] = {
        0.0959,//
        0.0990,//
        0.1048,//
        0.1108//
        };

        Double_t errv2pPb6[numsituations] = {
        0.0009,//
        0.0010,//
        0.0012,//
        0.0015//
        };

        TGraphErrors *v2syst6 = new TGraphErrors(numsituations,sit,v2pPb6,errsit,errv2pPb6);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst6->SetMarkerStyle(25);
         v2syst6->SetMarkerColor(kGreen);
         v2syst6->SetLineColor(kGreen);
        v2syst6->Draw("P");
    
    // Associate cut 60-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.2;
        errsit[i]=0;
        };

        Double_t v2pPb7[numsituations] = {
        0.1114,//
        0.1153,//
        0.1224,//
        0.1291//
        };

        Double_t errv2pPb7[numsituations] = {
        0.0017,//
        0.0019,//
        0.0023,//
        0.0027//
        };

        TGraphErrors *v2syst7 = new TGraphErrors(numsituations,sit,v2pPb7,errsit,errv2pPb7);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst7->SetMarkerStyle(25);
         v2syst7->SetMarkerColor(kBlue);
         v2syst7->SetLineColor(kBlue);
        v2syst7->Draw("P");
    
    
    // Associate cut 0-100 ALICE

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.15;
       errsit[i]=0;
       };
        
       Double_t v2pPb8[numsituations] = {
       0.0475,//
       0.0442,//
       0.0370,//
       0.0268//
       };

       Double_t errv2pPb8[numsituations] = {
       0.0003,//
       0.0002,//
       0.0003,//
       0.0005//
       };

       TGraphErrors *v2syst8 = new TGraphErrors(numsituations,sit,v2pPb8,errsit,errv2pPb8);
    
       v2syst8->SetMarkerStyle(20);
        v2syst8->SetMarkerColor(kRed);
        v2syst8->SetLineColor(kRed);
       v2syst8->Draw("P");
        
        
        
        // Associate cut 20-100 ALICE
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.1;
        errsit[i]=0;
        };

        Double_t v2pPb9[numsituations] = {
        0.0575,//
        0.0563,//
        0.0538,//
        0.0511//
        };

        Double_t errv2pPb9[numsituations] = {
        0.0002,//
        0.0002,//
        0.0002,//
        0.0003//
        };

        TGraphErrors *v2syst9 = new TGraphErrors(numsituations,sit,v2pPb9,errsit,errv2pPb9);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst9->SetMarkerStyle(20);
         v2syst9->SetMarkerColor(kOrange);
         v2syst9->SetLineColor(kOrange);
        v2syst9->Draw("P");
        
        
        // Associate cut 40-100 ALICE
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.05;
            errsit[i]=0;
            };

            Double_t v2pPb10[numsituations] = {
            0.0624,//
            0.0620,//
            0.0613,//
            0.0606//
            };

            Double_t errv2pPb10[numsituations] = {
            0.0002,//
            0.0002,//
            0.0002,//
            0.0003//
            };

            TGraphErrors *v2syst10 = new TGraphErrors(numsituations,sit,v2pPb10,errsit,errv2pPb10);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst10->SetMarkerStyle(20);
             v2syst10->SetMarkerColor(kGreen);
             v2syst10->SetLineColor(kGreen);
            v2syst10->Draw("P");
    
    // Associate cut 60-100 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i;
        errsit[i]=0;
        };

        Double_t v2pPb11[numsituations] = {
        0.0664,//
        0.0667,//
        0.0671,//
        0.0678//
        };

        Double_t errv2pPb11[numsituations] = {
        0.0002,//
        0.0002,//
        0.0003,//
        0.0003//
        };

        TGraphErrors *v2syst11 = new TGraphErrors(numsituations,sit,v2pPb11,errsit,errv2pPb11);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst11->SetMarkerStyle(20);
         v2syst11->SetMarkerColor(kBlue);
         v2syst11->SetLineColor(kBlue);
        v2syst11->Draw("P");
    
    // Associate cut 20-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.05;
        errsit[i]=0;
        };

        Double_t v2pPb12[numsituations] = {
        0.0572,//
        0.0560,//
        0.0535,//
        0.0506//
        };

        Double_t errv2pPb12[numsituations] = {
        0.0002,//
        0.0002,//
        0.0002,//
        0.0003//
        };

        TGraphErrors *v2syst12 = new TGraphErrors(numsituations,sit,v2pPb12,errsit,errv2pPb12);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst12->SetMarkerStyle(24);
         v2syst12->SetMarkerColor(kOrange);
         v2syst12->SetLineColor(kOrange);
        v2syst12->Draw("P");
        
    // Associate cut 40-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb13[numsituations] = {
        0.0619,//
        0.0616,//
        0.0607,//
        0.0598//
        };

        Double_t errv2pPb13[numsituations] = {
        0.0002,//
        0.0002,//
        0.0002,//
        0.0003//
        };

        TGraphErrors *v2syst13 = new TGraphErrors(numsituations,sit,v2pPb13,errsit,errv2pPb13);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst13->SetMarkerStyle(24);
         v2syst13->SetMarkerColor(kGreen);
         v2syst13->SetLineColor(kGreen);
        v2syst13->Draw("P");
    
    // Associate cut 60-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.15;
        errsit[i]=0;
        };

        Double_t v2pPb14[numsituations] = {
        0.0656,//
        0.0658,//
        0.0661,//
        0.0665//
        };

        Double_t errv2pPb14[numsituations] = {
        0.0002,//
        0.0002,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst14 = new TGraphErrors(numsituations,sit,v2pPb14,errsit,errv2pPb14);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst14->SetMarkerStyle(24);
         v2syst14->SetMarkerColor(kBlue);
         v2syst14->SetLineColor(kBlue);
        v2syst14->Draw("P");
        
    
        
    
    
        
        TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
        legend5->SetHeader("SPDTracklets - 10mrad");
          legend5->SetFillColorAlpha(kWhite, 0.);
          legend5->SetBorderSize(0);
            legend5->SetTextFont(42);
          legend5->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Template fit 0-100");
            legend5->AddEntry(v2syst,message);
       sprintf(message,"Template fit 20-100");
        legend5->AddEntry(v2syst2,message);
        sprintf(message,"Template fit 40-100");
        legend5->AddEntry(v2syst3,message);
    sprintf(message,"Template fit 60-100");
    legend5->AddEntry(v2syst4,message);
    sprintf(message,"Template fit 20-80");
    legend5->AddEntry(v2syst5,message);
    sprintf(message,"Template fit 40-80");
    legend5->AddEntry(v2syst6,message);
    sprintf(message,"Template fit 60-80");
    legend5->AddEntry(v2syst7,message);
    sprintf(message,"Yield sub 0-100");
            legend5->AddEntry(v2syst8,message);
       sprintf(message,"Yield sub 20-100");
        legend5->AddEntry(v2syst9,message);
        sprintf(message,"Yield sub 40-100");
        legend5->AddEntry(v2syst10,message);
    sprintf(message,"Yield sub 60-100");
    legend5->AddEntry(v2syst11,message);
    sprintf(message,"Yield sub 20-80");
    legend5->AddEntry(v2syst12,message);
    sprintf(message,"Yield sub 40-80");
    legend5->AddEntry(v2syst13,message);
    sprintf(message,"Yield sub 60-80");
    legend5->AddEntry(v2syst14,message);
          legend5->Draw();
}

void CentraliteTKLV0M5(){
    // Figure V0M 5mrad
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

        
       // Associate cut 0-100 ATLAS
        
       const Int_t numsituations=4;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.5;
       errsit[i]=0;
       };
        
       Double_t v2pPb[numsituations] = {
       0.0684,//
       0.0709,//
       0.0775,//
       0.0791//
       };

       Double_t errv2pPb[numsituations] = {
       0.0009,//
       0.0008,//
       0.0013,//
       0.0027//
       };

       std::string names[numsituations] = {
       "0-1",
       "1-3",
       "3-5",
        "5-10"
       };

       TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        v2syst->SetTitle("v_{2,tkl} for different central classes");
       TAxis *ax = v2syst->GetHistogram()->GetXaxis();
       Double_t x1 = ax->GetBinLowEdge(1);
       Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
       v2syst->GetHistogram()->GetXaxis()->Set(4,x1-0.5,x2+0.5);
        v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

       for(Int_t k=0;k<numsituations;k++){
       v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
       }

       v2syst->SetMarkerStyle(21);
        v2syst->SetMarkerColor(kRed);
        v2syst->SetLineColor(kRed);
       v2syst->Draw("AP");
        
        
        
        // Associate cut 20-100 ATLAS
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.45;
        errsit[i]=0;
        };

        Double_t v2pPb2[numsituations] = {
        0.0720,//
        0.0741,//
        0.0778,//
        0.0783//
        };

        Double_t errv2pPb2[numsituations] = {
        0.0006,//
        0.0005,//
        0.0006,//
        0.0008//
        };

        TGraphErrors *v2syst2 = new TGraphErrors(numsituations,sit,v2pPb2,errsit,errv2pPb2);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst2->SetMarkerStyle(21);
         v2syst2->SetMarkerColor(kOrange);
         v2syst2->SetLineColor(kOrange);
        v2syst2->Draw("P");
        
        
        // Associate cut 40-100 ATLAS
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.4;
            errsit[i]=0;
            };

            Double_t v2pPb3[numsituations] = {
            0.0768,//
            0.0792,//
            0.0831,//
            0.0847//
            };

            Double_t errv2pPb3[numsituations] = {
            0.0006,//
            0.0005,//
            0.0006,//
            0.0008//
            };

            TGraphErrors *v2syst3 = new TGraphErrors(numsituations,sit,v2pPb3,errsit,errv2pPb3);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst3->SetMarkerStyle(21);
             v2syst3->SetMarkerColor(kGreen);
             v2syst3->SetLineColor(kGreen);
            v2syst3->Draw("P");
    
    // Associate cut 60-100 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.35;
        errsit[i]=0;
        };

        Double_t v2pPb4[numsituations] = {
        0.0815,//
        0.0842,//
        0.0884,//
        0.0906//
        };

        Double_t errv2pPb4[numsituations] = {
        0.0006,//
        0.0006,//
        0.0008,//
        0.0009//
        };

        TGraphErrors *v2syst4 = new TGraphErrors(numsituations,sit,v2pPb4,errsit,errv2pPb4);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst4->SetMarkerStyle(21);
         v2syst4->SetMarkerColor(kBlue);
         v2syst4->SetLineColor(kBlue);
        v2syst4->Draw("P");
    
    // Associate cut 20-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.3;
        errsit[i]=0;
        };

        Double_t v2pPb5[numsituations] = {
        0.0716,//
        0.0736,//
        0.0773,//
        0.0776//
        };

        Double_t errv2pPb5[numsituations] = {
        0.0006,//
        0.0005,//
        0.0007,//
        0.0008//
        };

        TGraphErrors *v2syst5 = new TGraphErrors(numsituations,sit,v2pPb5,errsit,errv2pPb5);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst5->SetMarkerStyle(25);
         v2syst5->SetMarkerColor(kOrange);
         v2syst5->SetLineColor(kOrange);
        v2syst5->Draw("P");
        
    // Associate cut 40-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.25;
        errsit[i]=0;
        };

        Double_t v2pPb6[numsituations] = {
        0.0759,//
        0.0783,//
        0.0821,//
        0.0835//
        };

        Double_t errv2pPb6[numsituations] = {
        0.0006,//
        0.0006,//
        0.0007,//
        0.0008//
        };

        TGraphErrors *v2syst6 = new TGraphErrors(numsituations,sit,v2pPb6,errsit,errv2pPb6);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst6->SetMarkerStyle(25);
         v2syst6->SetMarkerColor(kGreen);
         v2syst6->SetLineColor(kGreen);
        v2syst6->Draw("P");
    
    // Associate cut 60-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.2;
        errsit[i]=0;
        };

        Double_t v2pPb7[numsituations] = {
        0.0801,//
        0.0827,//
        0.0868,//
        0.0887//
        };

        Double_t errv2pPb7[numsituations] = {
        0.0007,//
        0.0008,//
        0.0009,//
        0.0011//
        };

        TGraphErrors *v2syst7 = new TGraphErrors(numsituations,sit,v2pPb7,errsit,errv2pPb7);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst7->SetMarkerStyle(25);
         v2syst7->SetMarkerColor(kBlue);
         v2syst7->SetLineColor(kBlue);
        v2syst7->Draw("P");
    
    
    // Associate cut 0-100 ALICE

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.15;
       errsit[i]=0;
       };
        
       Double_t v2pPb8[numsituations] = {
       0.0456,//
       0.0435,//
       0.0391,//
       0.0285//
       };

       Double_t errv2pPb8[numsituations] = {
       0.0005,//
       0.0004,//
       0.0005,//
       0.0008//
       };

       TGraphErrors *v2syst8 = new TGraphErrors(numsituations,sit,v2pPb8,errsit,errv2pPb8);
    
       v2syst8->SetMarkerStyle(20);
        v2syst8->SetMarkerColor(kRed);
        v2syst8->SetLineColor(kRed);
       v2syst8->Draw("P");
        
        
        
        // Associate cut 20-100 ALICE
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.1;
        errsit[i]=0;
        };

        Double_t v2pPb9[numsituations] = {
        0.0585,//
        0.0582,//
        0.0574,//
        0.0538//
        };

        Double_t errv2pPb9[numsituations] = {
        0.0004,//
        0.0003,//
        0.0004,//
        0.0005//
        };

        TGraphErrors *v2syst9 = new TGraphErrors(numsituations,sit,v2pPb9,errsit,errv2pPb9);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst9->SetMarkerStyle(20);
         v2syst9->SetMarkerColor(kOrange);
         v2syst9->SetLineColor(kOrange);
        v2syst9->Draw("P");
        
        
        // Associate cut 40-100 ALICE
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.05;
            errsit[i]=0;
            };

            Double_t v2pPb10[numsituations] = {
            0.0654,//
            0.0660,//
            0.0666,//
            0.0650//
            };

            Double_t errv2pPb10[numsituations] = {
            0.0004,//
            0.0003,//
            0.0004,//
            0.0005//
            };

            TGraphErrors *v2syst10 = new TGraphErrors(numsituations,sit,v2pPb10,errsit,errv2pPb10);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst10->SetMarkerStyle(20);
             v2syst10->SetMarkerColor(kGreen);
             v2syst10->SetLineColor(kGreen);
            v2syst10->Draw("P");
    
    // Associate cut 60-100 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i;
        errsit[i]=0;
        };

        Double_t v2pPb11[numsituations] = {
        0.0713,//
        0.0725,//
        0.0741,//
        0.0740//
        };

        Double_t errv2pPb11[numsituations] = {
        0.0004,//
        0.0004,//
        0.0005,//
        0.0006//
        };

        TGraphErrors *v2syst11 = new TGraphErrors(numsituations,sit,v2pPb11,errsit,errv2pPb11);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst11->SetMarkerStyle(20);
         v2syst11->SetMarkerColor(kBlue);
         v2syst11->SetLineColor(kBlue);
        v2syst11->Draw("P");
    
    // Associate cut 20-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.05;
        errsit[i]=0;
        };

        Double_t v2pPb12[numsituations] = {
        0.0581,//
        0.0577,//
        0.0568,//
        0.0530//
        };

        Double_t errv2pPb12[numsituations] = {
        0.0004,//
        0.0003,//
        0.0004,//
        0.0005//
        };

        TGraphErrors *v2syst12 = new TGraphErrors(numsituations,sit,v2pPb12,errsit,errv2pPb12);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst12->SetMarkerStyle(24);
         v2syst12->SetMarkerColor(kOrange);
         v2syst12->SetLineColor(kOrange);
        v2syst12->Draw("P");
        
    // Associate cut 40-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb13[numsituations] = {
        0.0646,//
        0.0650,//
        0.0655,//
        0.0636//
        };

        Double_t errv2pPb13[numsituations] = {
        0.0004,//
        0.0004,//
        0.0004,//
        0.0005//
        };

        TGraphErrors *v2syst13 = new TGraphErrors(numsituations,sit,v2pPb13,errsit,errv2pPb13);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst13->SetMarkerStyle(24);
         v2syst13->SetMarkerColor(kGreen);
         v2syst13->SetLineColor(kGreen);
        v2syst13->Draw("P");
    
    // Associate cut 60-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.15;
        errsit[i]=0;
        };

        Double_t v2pPb14[numsituations] = {
        0.0697,//
        0.0707,//
        0.0721,//
        0.0716//
        };

        Double_t errv2pPb14[numsituations] = {
        0.0005,//
        0.0005,//
        0.0006,//
        0.0007//
        };

        TGraphErrors *v2syst14 = new TGraphErrors(numsituations,sit,v2pPb14,errsit,errv2pPb14);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst14->SetMarkerStyle(24);
         v2syst14->SetMarkerColor(kBlue);
         v2syst14->SetLineColor(kBlue);
        v2syst14->Draw("P");
        
    
        
    
    
        
        TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
        legend5->SetHeader("V0M - 5mrad");
          legend5->SetFillColorAlpha(kWhite, 0.);
          legend5->SetBorderSize(0);
            legend5->SetTextFont(42);
          legend5->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Template fit 0-100");
            legend5->AddEntry(v2syst,message);
       sprintf(message,"Template fit 20-100");
        legend5->AddEntry(v2syst2,message);
        sprintf(message,"Template fit 40-100");
        legend5->AddEntry(v2syst3,message);
    sprintf(message,"Template fit 60-100");
    legend5->AddEntry(v2syst4,message);
    sprintf(message,"Template fit 20-80");
    legend5->AddEntry(v2syst5,message);
    sprintf(message,"Template fit 40-80");
    legend5->AddEntry(v2syst6,message);
    sprintf(message,"Template fit 60-80");
    legend5->AddEntry(v2syst7,message);
    sprintf(message,"Yield sub 0-100");
            legend5->AddEntry(v2syst8,message);
       sprintf(message,"Yield sub 20-100");
        legend5->AddEntry(v2syst9,message);
        sprintf(message,"Yield sub 40-100");
        legend5->AddEntry(v2syst10,message);
    sprintf(message,"Yield sub 60-100");
    legend5->AddEntry(v2syst11,message);
    sprintf(message,"Yield sub 20-80");
    legend5->AddEntry(v2syst12,message);
    sprintf(message,"Yield sub 40-80");
    legend5->AddEntry(v2syst13,message);
    sprintf(message,"Yield sub 60-80");
    legend5->AddEntry(v2syst14,message);
          legend5->Draw();
}

void CentraliteTKLSPDT5(){
    // Figure SPDT 5mrad
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

        
       // Associate cut 0-100 ATLAS
        
       const Int_t numsituations=4;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.5;
       errsit[i]=0;
       };
        
       Double_t v2pPb[numsituations] = {
       0.0881,//
       0.0914,//
       0.1012,//
       0.1161//
       };

       Double_t errv2pPb[numsituations] = {
       0.0011,//
       0.0010,//
       0.0018,//
       0.0048//
       };

       std::string names[numsituations] = {
       "0-1",
       "1-3",
       "3-5",
        "5-10"
       };

       TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        v2syst->SetTitle("v_{2,tkl} for different central classes");
       TAxis *ax = v2syst->GetHistogram()->GetXaxis();
       Double_t x1 = ax->GetBinLowEdge(1);
       Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
       v2syst->GetHistogram()->GetXaxis()->Set(4,x1-0.5,x2+0.5);
        v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

       for(Int_t k=0;k<numsituations;k++){
       v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
       }

       v2syst->SetMarkerStyle(21);
        v2syst->SetMarkerColor(kRed);
        v2syst->SetLineColor(kRed);
       v2syst->Draw("AP");
        
        
        
        // Associate cut 20-100 ATLAS
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.45;
        errsit[i]=0;
        };

        Double_t v2pPb2[numsituations] = {
        0.0997,//
        0.1028,//
        0.1090,//
        0.1138//
        };

        Double_t errv2pPb2[numsituations] = {
        0.0007,//
        0.0007,//
        0.0009,//
        0.0012//
        };

        TGraphErrors *v2syst2 = new TGraphErrors(numsituations,sit,v2pPb2,errsit,errv2pPb2);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst2->SetMarkerStyle(21);
         v2syst2->SetMarkerColor(kOrange);
         v2syst2->SetLineColor(kOrange);
        v2syst2->Draw("P");
        
        
        // Associate cut 40-100 ATLAS
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.4;
            errsit[i]=0;
            };

            Double_t v2pPb3[numsituations] = {
            0.1180,//
            0.1220,//
            0.1299,//
            0.1370//
            };

            Double_t errv2pPb3[numsituations] = {
            0.0009,//
            0.0010,//
            0.0012,//
            0.0014//
            };

            TGraphErrors *v2syst3 = new TGraphErrors(numsituations,sit,v2pPb3,errsit,errv2pPb3);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst3->SetMarkerStyle(21);
             v2syst3->SetMarkerColor(kGreen);
             v2syst3->SetLineColor(kGreen);
            v2syst3->Draw("P");
    
    // Associate cut 60-100 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.35;
        errsit[i]=0;
        };

        Double_t v2pPb4[numsituations] = {
        0.1328,//
        0.1369,//
        0.1445,//
        0.1515//
        };

        Double_t errv2pPb4[numsituations] = {
        0.0013,//
        0.0014,//
        0.0016,//
        0.0019//
        };

        TGraphErrors *v2syst4 = new TGraphErrors(numsituations,sit,v2pPb4,errsit,errv2pPb4);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst4->SetMarkerStyle(21);
         v2syst4->SetMarkerColor(kBlue);
         v2syst4->SetLineColor(kBlue);
        v2syst4->Draw("P");
    
    // Associate cut 20-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.3;
        errsit[i]=0;
        };

        Double_t v2pPb5[numsituations] = {
        0.0987,//
        0.1017,//
        0.1078,//
        0.1124//
        };

        Double_t errv2pPb5[numsituations] = {
        0.0007,//
        0.0007,//
        0.0009,//
        0.0012//
        };

        TGraphErrors *v2syst5 = new TGraphErrors(numsituations,sit,v2pPb5,errsit,errv2pPb5);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst5->SetMarkerStyle(25);
         v2syst5->SetMarkerColor(kOrange);
         v2syst5->SetLineColor(kOrange);
        v2syst5->Draw("P");
        
    // Associate cut 40-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.25;
        errsit[i]=0;
        };

        Double_t v2pPb6[numsituations] = {
        0.1163,//
        0.1203,//
        0.1282,//
        0.1354//
        };

        Double_t errv2pPb6[numsituations] = {
        0.0010,//
        0.0010,//
        0.0013,//
        0.0016//
        };

        TGraphErrors *v2syst6 = new TGraphErrors(numsituations,sit,v2pPb6,errsit,errv2pPb6);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst6->SetMarkerStyle(25);
         v2syst6->SetMarkerColor(kGreen);
         v2syst6->SetLineColor(kGreen);
        v2syst6->Draw("P");
    
    // Associate cut 60-80 ATLAS
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.2;
        errsit[i]=0;
        };

        Double_t v2pPb7[numsituations] = {
        0.1319,//
        0.1363,//
        0.1447,//
        0.1524//
        };

        Double_t errv2pPb7[numsituations] = {
        0.0017,//
        0.0019,//
        0.0022,//
        0.0026//
        };

        TGraphErrors *v2syst7 = new TGraphErrors(numsituations,sit,v2pPb7,errsit,errv2pPb7);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst7->SetMarkerStyle(25);
         v2syst7->SetMarkerColor(kBlue);
         v2syst7->SetLineColor(kBlue);
        v2syst7->Draw("P");
    
    
    // Associate cut 0-100 ALICE

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.15;
       errsit[i]=0;
       };
        
       Double_t v2pPb8[numsituations] = {
       0.0572,//
       0.0540,//
       0.0469,//
       0.0340//
       };

       Double_t errv2pPb8[numsituations] = {
       0.0003,//
       0.0003,//
       0.0004,//
       0.0007//
       };

       TGraphErrors *v2syst8 = new TGraphErrors(numsituations,sit,v2pPb8,errsit,errv2pPb8);
    
       v2syst8->SetMarkerStyle(20);
        v2syst8->SetMarkerColor(kRed);
        v2syst8->SetLineColor(kRed);
       v2syst8->Draw("P");
        
        
        
        // Associate cut 20-100 ALICE
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i-0.1;
        errsit[i]=0;
        };

        Double_t v2pPb9[numsituations] = {
        0.0706,//
        0.0701,//
        0.0688,//
        0.0658//
        };

        Double_t errv2pPb9[numsituations] = {
        0.0003,//
        0.0002,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst9 = new TGraphErrors(numsituations,sit,v2pPb9,errsit,errv2pPb9);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst9->SetMarkerStyle(20);
         v2syst9->SetMarkerColor(kOrange);
         v2syst9->SetLineColor(kOrange);
        v2syst9->Draw("P");
        
        
        // Associate cut 40-100 ALICE
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i-0.05;
            errsit[i]=0;
            };

            Double_t v2pPb10[numsituations] = {
            0.0775,//
            0.0781,//
            0.0791,//
            0.0789//
            };

            Double_t errv2pPb10[numsituations] = {
            0.0003,//
            0.0002,//
            0.0003,//
            0.0003//
            };

            TGraphErrors *v2syst10 = new TGraphErrors(numsituations,sit,v2pPb10,errsit,errv2pPb10);
        //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
        //    Double_t x12 = ax2->GetBinLowEdge(1);
        //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
        //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
        //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
        //
        //    for(Int_t k=0;k<numsituations;k++){
        //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        //    }

            v2syst10->SetMarkerStyle(20);
             v2syst10->SetMarkerColor(kGreen);
             v2syst10->SetLineColor(kGreen);
            v2syst10->Draw("P");
    
    // Associate cut 60-100 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i;
        errsit[i]=0;
        };

        Double_t v2pPb11[numsituations] = {
        0.0843,//
        0.0859,//
        0.0888,//
        0.0909//
        };

        Double_t errv2pPb11[numsituations] = {
        0.0003,//
        0.0003,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst11 = new TGraphErrors(numsituations,sit,v2pPb11,errsit,errv2pPb11);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst11->SetMarkerStyle(20);
         v2syst11->SetMarkerColor(kBlue);
         v2syst11->SetLineColor(kBlue);
        v2syst11->Draw("P");
    
    // Associate cut 20-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.05;
        errsit[i]=0;
        };

        Double_t v2pPb12[numsituations] = {
        0.0702,//
        0.0696,//
        0.0682,//
        0.0650//
        };

        Double_t errv2pPb12[numsituations] = {
        0.0003,//
        0.0002,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst12 = new TGraphErrors(numsituations,sit,v2pPb12,errsit,errv2pPb12);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst12->SetMarkerStyle(24);
         v2syst12->SetMarkerColor(kOrange);
         v2syst12->SetLineColor(kOrange);
        v2syst12->Draw("P");
        
    // Associate cut 40-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb13[numsituations] = {
        0.0766,//
        0.0770,//
        0.0777,//
        0.0772//
        };

        Double_t errv2pPb13[numsituations] = {
        0.0003,//
        0.0002,//
        0.0003,//
        0.0004//
        };

        TGraphErrors *v2syst13 = new TGraphErrors(numsituations,sit,v2pPb13,errsit,errv2pPb13);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst13->SetMarkerStyle(24);
         v2syst13->SetMarkerColor(kGreen);
         v2syst13->SetLineColor(kGreen);
        v2syst13->Draw("P");
    
    // Associate cut 60-80 ALICE
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.15;
        errsit[i]=0;
        };

        Double_t v2pPb14[numsituations] = {
        0.0821,//
        0.0834,//
        0.0857,//
        0.0870//
        };

        Double_t errv2pPb14[numsituations] = {
        0.0003,//
        0.0003,//
        0.0004,//
        0.0005//
        };

        TGraphErrors *v2syst14 = new TGraphErrors(numsituations,sit,v2pPb14,errsit,errv2pPb14);
    //    TAxis *ax2 = v2syst2->GetHistogram()->GetXaxis();
    //    Double_t x12 = ax2->GetBinLowEdge(1);
    //    Double_t x22 = ax2->GetBinUpEdge(ax2->GetNbins());
    //    v2syst2->GetHistogram()->GetXaxis()->Set(3,x12,x22);
    //     v2syst2->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    //
    //    for(Int_t k=0;k<numsituations;k++){
    //    v2syst2->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    //    }

        v2syst14->SetMarkerStyle(24);
         v2syst14->SetMarkerColor(kBlue);
         v2syst14->SetLineColor(kBlue);
        v2syst14->Draw("P");
        
    
        
    
    
        
        TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
        legend5->SetHeader("SPDTracklets - 5mrad");
          legend5->SetFillColorAlpha(kWhite, 0.);
          legend5->SetBorderSize(0);
            legend5->SetTextFont(42);
          legend5->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Template fit 0-100");
            legend5->AddEntry(v2syst,message);
       sprintf(message,"Template fit 20-100");
        legend5->AddEntry(v2syst2,message);
        sprintf(message,"Template fit 40-100");
        legend5->AddEntry(v2syst3,message);
    sprintf(message,"Template fit 60-100");
    legend5->AddEntry(v2syst4,message);
    sprintf(message,"Template fit 20-80");
    legend5->AddEntry(v2syst5,message);
    sprintf(message,"Template fit 40-80");
    legend5->AddEntry(v2syst6,message);
    sprintf(message,"Template fit 60-80");
    legend5->AddEntry(v2syst7,message);
    sprintf(message,"Yield sub 0-100");
            legend5->AddEntry(v2syst8,message);
       sprintf(message,"Yield sub 20-100");
        legend5->AddEntry(v2syst9,message);
        sprintf(message,"Yield sub 40-100");
        legend5->AddEntry(v2syst10,message);
    sprintf(message,"Yield sub 60-100");
    legend5->AddEntry(v2syst11,message);
    sprintf(message,"Yield sub 20-80");
    legend5->AddEntry(v2syst12,message);
    sprintf(message,"Yield sub 40-80");
    legend5->AddEntry(v2syst13,message);
    sprintf(message,"Yield sub 60-80");
    legend5->AddEntry(v2syst14,message);
          legend5->Draw();
}
