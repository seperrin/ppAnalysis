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
 # include <fstream>
# include <sstream>
 # include <string>
# include <iterator>
 # include "TH1F.h"
 # include "TH1D.h"
 # include "TH2F.h"
 # include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"


void SystematicsPlotTKLV0MpPb();
void SystematicsPlotTKLSPDTpPb();
void SystematicsPlotTKLSPDCpPb();
void SystematicsPlotTKLV0MATLAS();
void SystematicsPlotTKLSPDTATLAS();
void SystematicsPlotTKLSPDCATLAS();
void SystematicsPlotTKLV0MpPbDPHI();
void SystematicsPlotTKLSPDTpPbDPHI();
void SystematicsPlotTKLSPDCpPbDPHI();

void SystematicsPlotTKLV0MpPb()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };
    
   Double_t v2pPb[numsituations] = {
   0.052475,
   0.066157,
   0.0570693
   };

   Double_t errv2pPb[numsituations] = {
   0.000523736,
   0.00065597,
   0.000654495
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.05,0.07);
    v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.0617764,
    0.0776326,
    0.0637896
    };

    Double_t errv2pPb2[numsituations] = {
    0.000616509,
    0.000769559,
    0.000731431
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.0636416,
        0.0803551,
        0.0668897
        };

        Double_t errv2pPb3[numsituations] = {
        0.000839039,
        0.000994038,
        0.0010611
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
     //   v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("V0M - Subtraction or Template fit");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - Subtraction");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - Template fit");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - Away-side Gaussian");
    //legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}

void SystematicsPlotTKLSPDTpPb()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };

   Double_t v2pPb[numsituations] = {
   0.0617831,
   0.0784327,
   0.0668407
   };

   Double_t errv2pPb[numsituations] = {
   0.000379367,
   0.000476329,
   0.000448392
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.101491,
    0.124272,
    0.115903
    };

    Double_t errv2pPb2[numsituations] = {
    0.00175828,
    0.00182721,
    0.0027521
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.0658283,
        0.086481,
        0.0934327
        };

        Double_t errv2pPb3[numsituations] = {
        0.000364606,
        0.000456387,
        0.000550942
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
      //  v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("SPDTracklets - Subtraction or Template fit");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - Subtraction");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - Template fit");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - Away-side Gaussian");
   // legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}

void SystematicsPlotTKLSPDCpPb()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };

   Double_t v2pPb[numsituations] = {
   0.0606523,
   0.0773282,
   0.0660357
   };

   Double_t errv2pPb[numsituations] = {
   0.000395664,
   0.000494277,
   0.000341794
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.0921614,
    0.115051,
    0.106755
    };

    Double_t errv2pPb2[numsituations] = {
    0.001505,
    0.00157482,
    0.00140208
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.0623351,
        0.0829773,
        0
        };

        Double_t errv2pPb3[numsituations] = {
        0.000390879,
        0.00047783,
        99
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
       // v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("SPDClusters - Subtraction or Template fit");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - Subtraction");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - Template fit");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - Away-side Gaussian");
   // legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}

void SystematicsPlotTKLV0MATLAS()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };

   Double_t v2pPb[numsituations] = {
   0.032528,
   0.0248328,
   0.03505
   };

   Double_t errv2pPb[numsituations] = {
   0.00305687,
   0.0010619,
   0.00189249
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.0232528,
    0.05248328,
    0.033505
    };

    Double_t errv2pPb2[numsituations] = {
    0.003605687,
    0.00710619,
    0.00689249
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.05232528,
        0.045248328,
        0.0633505
        };

        Double_t errv2pPb3[numsituations] = {
        0.0013605687,
        0.006710619,
        0.004689249
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
        v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("V0M - Template fit - z cut");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - 10cm");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - 8cm");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - 12cm");
    legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}

void SystematicsPlotTKLSPDTATLAS()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };

   Double_t v2pPb[numsituations] = {
   0.032528,
   0.0248328,
   0.03505
   };

   Double_t errv2pPb[numsituations] = {
   0.00305687,
   0.0010619,
   0.00189249
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.0232528,
    0.05248328,
    0.033505
    };

    Double_t errv2pPb2[numsituations] = {
    0.003605687,
    0.00710619,
    0.00689249
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.05232528,
        0.045248328,
        0.0633505
        };

        Double_t errv2pPb3[numsituations] = {
        0.0013605687,
        0.006710619,
        0.004689249
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
        v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("SPDTracklets - Template fit - z cut");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - 10cm");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - 8cm");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - 12cm");
    legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}

void SystematicsPlotTKLSPDCATLAS()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };

   Double_t v2pPb[numsituations] = {
   0.032528,
   0.0248328,
   0.03505
   };

   Double_t errv2pPb[numsituations] = {
   0.00305687,
   0.0010619,
   0.00189249
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.0232528,
    0.05248328,
    0.033505
    };

    Double_t errv2pPb2[numsituations] = {
    0.003605687,
    0.00710619,
    0.00689249
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.05232528,
        0.045248328,
        0.0633505
        };

        Double_t errv2pPb3[numsituations] = {
        0.0013605687,
        0.006710619,
        0.004689249
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
        v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("SPDClusters - Template fit - z cut");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - 10cm");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - 8cm");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - 12cm");
    legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}


void SystematicsPlotTKLV0MpPbDPHI()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // V0M
    
    
    
   const Int_t numsituations=5;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };
    
   Double_t v2pPb[numsituations] = {
    0.0887559,
    0.0631085,
   0.066157,
   0.052475,
   0.0391958
   };

   Double_t errv2pPb[numsituations] = {
   0.00376338,
   0.00204661,
   0.00065597,
   0.000523736,
   0.000515732
   };

   std::string names[numsituations] = {
   "1mrad",
   "2mrad",
   "5mrad",
    "10mrad",
       "No cut",
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("DPhiCut impact on v_{2,tkl}");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(5,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //SPD T
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.117968,
    0.0857152,
    0.0784327,
        0.0617831,
        0.0468189
    };

    Double_t errv2pPb2[numsituations] = {
    0.0027357,
    0.00136307,
    0.000476329,
        0.000379367,
        0.000366046
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //SPDC
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.111693,
        0.0825727,
        0.0773282,
            0.0606523,
            0.045706
        };

        Double_t errv2pPb3[numsituations] = {
        0.00292278,
        0.00144142,
        0.000494277,
            0.000395664,
            0.000383984
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
        v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("Centrality Estimator");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"V0M");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"SPDTracklets");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"SPDClusters");
    legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}

void SystematicsPlotTKLSPDTpPbDPHI()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };

   Double_t v2pPb[numsituations] = {
   0.0617831,
   0.0784327,
   0.0668407
   };

   Double_t errv2pPb[numsituations] = {
   0.000379367,
   0.000476329,
   0.000448392
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.101491,
    0.124272,
    0.115903
    };

    Double_t errv2pPb2[numsituations] = {
    0.00175828,
    0.00182721,
    0.0027521
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.0658283,
        0.086481,
        0.0934327
        };

        Double_t errv2pPb3[numsituations] = {
        0.000364606,
        0.000456387,
        0.000550942
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
      //  v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("SPDTracklets - Subtraction or Template fit");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - Subtraction");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - Template fit");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - Away-side Gaussian");
   // legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}

void SystematicsPlotTKLSPDCpPbDPHI()
{
   TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    
   // Default case
    
    
    
   const Int_t numsituations=3;

   Double_t sit[numsituations];
   Double_t errsit[numsituations];

   for (Int_t i=0;i<numsituations;i++){
   sit[i]=i;
   errsit[i]=0;
   };

   Double_t v2pPb[numsituations] = {
   0.0606523,
   0.0773282,
   0.0660357
   };

   Double_t errv2pPb[numsituations] = {
   0.000395664,
   0.000494277,
   0.000341794
   };

   std::string names[numsituations] = {
   "0-5 40-100 10mrad",
   "0-5 40-100 5mrad",
   "0-5 60-100 10mrad"
   };

   TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
    v2syst->SetTitle("Systematics study");
   TAxis *ax = v2syst->GetHistogram()->GetXaxis();
   Double_t x1 = ax->GetBinLowEdge(1);
   Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
   v2syst->GetHistogram()->GetXaxis()->Set(3,x1-0.5,x2+0.5);
    v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
    v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");

   for(Int_t k=0;k<numsituations;k++){
   v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
   }

   v2syst->SetMarkerStyle(21);
    v2syst->SetMarkerColor(kBlue);
    v2syst->SetLineColor(kBlue);
   v2syst->Draw("AP");
    
    
    
    //Alternative 1
    
    
    
    for (Int_t i=0;i<numsituations;i++){
    sit[i]=i-0.1;
    errsit[i]=0;
    };

    Double_t v2pPb2[numsituations] = {
    0.0921614,
    0.115051,
    0.106755
    };

    Double_t errv2pPb2[numsituations] = {
    0.001505,
    0.00157482,
    0.00140208
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
     v2syst2->SetMarkerColor(kRed);
     v2syst2->SetLineColor(kRed);
    v2syst2->Draw("P");
    
    
    //Alternative 2
        
        
        for (Int_t i=0;i<numsituations;i++){
        sit[i]=i+0.1;
        errsit[i]=0;
        };

        Double_t v2pPb3[numsituations] = {
        0.0623351,
        0.0829773,
        0
        };

        Double_t errv2pPb3[numsituations] = {
        0.000390879,
        0.00047783,
        99
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
         v2syst3->SetMarkerColor(kBlack);
         v2syst3->SetLineColor(kBlack);
       // v2syst3->Draw("P");
    
    
    TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
    legend5->SetHeader("SPDClusters - Subtraction or Template fit");
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Default - Subtraction");
        legend5->AddEntry(v2syst,message);
   sprintf(message,"Alternative - Template fit");
    legend5->AddEntry(v2syst2,message);
    sprintf(message,"Alternative - Away-side Gaussian");
   // legend5->AddEntry(v2syst3,message);
      legend5->Draw();
}
