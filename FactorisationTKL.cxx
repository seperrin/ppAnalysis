// Macro qui plot la factorisation tkl


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

void FactorisationTKLV0M();
void FactorisationTKL();
void FactorisationTKLSPDT();
void FactorisationTKLSPDC();


void FactorisationTKLV0M(){
    // Figure V0M
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

        
       // Associate cut 1mrad
        
       const Int_t numsituations=4;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.1;
       errsit[i]=0;
       };
        
       Double_t v2pPb[numsituations] = {
       0.0887559,//
       0.0677103,//
       0.0674184,//
       0.0694063//
       };

       Double_t errv2pPb[numsituations] = {
       0.00376338,//
       0.00290572,//
       0.0016006,//
       0.00119471//
       };

       std::string names[numsituations] = {
       "1mrad",
       "2mrad",
       "5mrad",
        "10mrad"
       };

       TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        v2syst->SetTitle("v_{2,tkl} for different reference tracklet cuts");
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
        
        
        
        // Associate cut 2mrad
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i;
        errsit[i]=0;
        };

        Double_t v2pPb2[numsituations] = {
        0.0566688,//
        0.0631085,//
        0.068743,//
        0.0650263//
        };

        Double_t errv2pPb2[numsituations] = {
        0.00387948,//
        0.00204661,//
        0.00102938,//
        0.000838474//
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
        
        
        // Associate cut 5mrad
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i+0.1;
            errsit[i]=0;
            };

            Double_t v2pPb3[numsituations] = {
            0.05657,//
            0.0653024,//
            0.066157,//
            0.0601991//
            };

            Double_t errv2pPb3[numsituations] = {
            0.0023643,//
            0.00121174,//
            0.00065597,//
            0.000556342//
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
             v2syst3->SetMarkerColor(kYellow);
             v2syst3->SetLineColor(kYellow);
            v2syst3->Draw("P");
    
    // Associate cut 10mrad
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.2;
        errsit[i]=0;
        };

        Double_t v2pPb4[numsituations] = {
        0.0619827,//
        0.0630695,//
        0.0596104,//
        0.052475//
        };

        Double_t errv2pPb4[numsituations] = {
        0.00176867,//
        0.00103186,//
        0.000597737,//
        0.000523736//
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
         v2syst4->SetMarkerColor(kGreen);
         v2syst4->SetLineColor(kGreen);
        v2syst4->Draw("P");
    
    // Associate cut No
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.3;
        errsit[i]=0;
        };

        Double_t v2pPb5[numsituations] = {
        0.0573849,//
        0.055888,//
        0.0521362,//
        0.0457097//
        };

        Double_t errv2pPb5[numsituations] = {
        0.00165586,//
        0.00101373,//
        0.00059436,//
        0.000522758//
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

        v2syst5->SetMarkerStyle(21);
         v2syst5->SetMarkerColor(kBlue);
         v2syst5->SetLineColor(kBlue);
        v2syst5->Draw("P");
        
    
        
        TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
        legend5->SetHeader("V0M");
          legend5->SetFillColorAlpha(kWhite, 0.);
          legend5->SetBorderSize(0);
            legend5->SetTextFont(42);
          legend5->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Associate 1mrad");
            legend5->AddEntry(v2syst,message);
       sprintf(message,"Associate 2mrad");
        legend5->AddEntry(v2syst2,message);
        sprintf(message,"Associate 5mrad");
        legend5->AddEntry(v2syst3,message);
    sprintf(message,"Associate 10mrad");
    legend5->AddEntry(v2syst4,message);
    sprintf(message,"Associate None");
    legend5->AddEntry(v2syst5,message);
          legend5->Draw();
}

void FactorisationTKLSPDT(){
    // Figure SPDT
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

        
       // Associate cut 1mrad
        
       const Int_t numsituations=4;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.1;
       errsit[i]=0;
       };
        
       Double_t v2pPb[numsituations] = {
       0.117968,
       0.0915954,
       0.082218,
       0.0808574
       };

       Double_t errv2pPb[numsituations] = {
       0.0027357,
       0.0019517,
       0.00114161,
       0.000888346
       };

       std::string names[numsituations] = {
       "1mrad",
       "2mrad",
       "5mrad",
        "10mrad"
       };

       TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        v2syst->SetTitle("v_{2,tkl} for different reference tracklet cuts");
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
        
        
        
        // Associate cut 2mrad
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i;
        errsit[i]=0;
        };

        Double_t v2pPb2[numsituations] = {
        0.0874985,
        0.0857152,
        0.082124,
        0.0776924
        };

        Double_t errv2pPb2[numsituations] = {
        0.00236422,
        0.00136307,
        0.000748283,
        0.0006078
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
        
        
        // Associate cut 5mrad
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i+0.1;
            errsit[i]=0;
            };

            Double_t v2pPb3[numsituations] = {
            0.080152,
            0.0826714,
            0.0784327,
            0.0710069
            };

            Double_t errv2pPb3[numsituations] = {
            0.00151523,
            0.000850408,
            0.000476329,
            0.000405232
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
             v2syst3->SetMarkerColor(kYellow);
             v2syst3->SetLineColor(kYellow);
            v2syst3->Draw("P");
    
    // Associate cut 10mrad
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.2;
        errsit[i]=0;
        };

        Double_t v2pPb4[numsituations] = {
        0.079434,
        0.078431,
        0.0706797,
        0.0617831
        };

        Double_t errv2pPb4[numsituations] = {
        0.00123112,
        0.000729209,
        0.000430546,
        0.000379367
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
         v2syst4->SetMarkerColor(kGreen);
         v2syst4->SetLineColor(kGreen);
        v2syst4->Draw("P");
    
    // Associate cut No
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.3;
        errsit[i]=0;
        };

        Double_t v2pPb5[numsituations] = {
        0.0724977,
        0.069429,
        0.0631626,
        0.0542276
        };

        Double_t errv2pPb5[numsituations] = {
        0.00116038,
        0.000711612,
        0.000416423,
        0.000373213
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

        v2syst5->SetMarkerStyle(21);
         v2syst5->SetMarkerColor(kBlue);
         v2syst5->SetLineColor(kBlue);
        v2syst5->Draw("P");
        
    
        
        TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
        legend5->SetHeader("SPDTracklets");
          legend5->SetFillColorAlpha(kWhite, 0.);
          legend5->SetBorderSize(0);
            legend5->SetTextFont(42);
          legend5->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Associate 1mrad");
            legend5->AddEntry(v2syst,message);
       sprintf(message,"Associate 2mrad");
        legend5->AddEntry(v2syst2,message);
        sprintf(message,"Associate 5mrad");
        legend5->AddEntry(v2syst3,message);
    sprintf(message,"Associate 10mrad");
    legend5->AddEntry(v2syst4,message);
    sprintf(message,"Associate None");
    legend5->AddEntry(v2syst5,message);
          legend5->Draw();
}

void FactorisationTKLSPDC(){
    // Figure SPDC
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

        
       // Associate cut 1mrad
        
       const Int_t numsituations=4;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i-0.1;
       errsit[i]=0;
       };
        
       Double_t v2pPb[numsituations] = {
       0.111693,//
       0.0872431,//
       0.0808892,//
       0.0788889//
       };

       Double_t errv2pPb[numsituations] = {
       0.00292278,//
       0.0020944,//
       0.00118886,//
       0.000930723//
       };

       std::string names[numsituations] = {
       "1mrad",
       "2mrad",
       "5mrad",
        "10mrad"
       };

       TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        v2syst->SetTitle("v_{2,tkl} for different reference tracklet cuts");
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
        
        
        
        // Associate cut 2mrad
        
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i;
        errsit[i]=0;
        };

        Double_t v2pPb2[numsituations] = {
        0.0837124,//
        0.0825727,//
        0.0801325,//
        0.0769878//
        };

        Double_t errv2pPb2[numsituations] = {
        0.00252027,//
        0.00144142,//
        0.000782649,//
        0.000626415//
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
        
        
        // Associate cut 5mrad
            
            
            for (Int_t i=0;i<numsituations;i++){
            sit[i]=i+0.1;
            errsit[i]=0;
            };

            Double_t v2pPb3[numsituations] = {
            0.0774592,//
            0.0799717,//
            0.0773282,//
            0.0695026//
            };

            Double_t errv2pPb3[numsituations] = {
            0.00160177,//
            0.000898507,//
            0.000494277,//
            0.000423023//
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
             v2syst3->SetMarkerColor(kYellow);
             v2syst3->SetLineColor(kYellow);
            v2syst3->Draw("P");
    
    // Associate cut 10mrad
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.2;
        errsit[i]=0;
        };

        Double_t v2pPb4[numsituations] = {
        0.0766938,//
        0.0781517,//
        0.0697926,//
        0.0606523//
        };

        Double_t errv2pPb4[numsituations] = {
        0.00130882,//
        0.000750821,//
        0.000446891,//
        0.000395664//
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
         v2syst4->SetMarkerColor(kGreen);
         v2syst4->SetLineColor(kGreen);
        v2syst4->Draw("P");
    
    // Associate cut No
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.3;
        errsit[i]=0;
        };

        Double_t v2pPb5[numsituations] = {
        0.0706388,//
        0.0684815,//
        0.061056,//
        0.0529027//
        };

        Double_t errv2pPb5[numsituations] = {
        0.00122506,//
        0.000740518,//
        0.000441521,//
        0.000391809//
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

        v2syst5->SetMarkerStyle(21);
         v2syst5->SetMarkerColor(kBlue);
         v2syst5->SetLineColor(kBlue);
        v2syst5->Draw("P");
        
    
        
        TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
        legend5->SetHeader("SPDClusters");
          legend5->SetFillColorAlpha(kWhite, 0.);
          legend5->SetBorderSize(0);
            legend5->SetTextFont(42);
          legend5->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Associate 1mrad");
            legend5->AddEntry(v2syst,message);
       sprintf(message,"Associate 2mrad");
        legend5->AddEntry(v2syst2,message);
        sprintf(message,"Associate 5mrad");
        legend5->AddEntry(v2syst3,message);
    sprintf(message,"Associate 10mrad");
    legend5->AddEntry(v2syst4,message);
    sprintf(message,"Associate None");
    legend5->AddEntry(v2syst5,message);
          legend5->Draw();
}
