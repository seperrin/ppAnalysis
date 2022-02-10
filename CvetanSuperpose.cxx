// Macro qui superpose mon cross-check et les plots du papier de Cvetan


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


void CvetanSuperpose(){
    
//    TH2F* YTklCentral(NULL);
   // TGraphMultiErrors* Fig5(NULL);
//    TH2F* YTklNum(NULL);
//    TH2F* YTklDen(NULL);

    double PtBin[] = {0,2,3,4,6,8};
    int NbPtBins = 5;
    
// Figure 5

       auto c5 = new TCanvas("c5","c5",200,10,600,400);
       double ax[]      = {1, 2.5, 3.5, 5, 7};
       double ay[]      = {-0.014, 0.004, 0.038, 0.092, 0.033};
       double aexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
       double aexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
       double* aeylstat = new double[5]  {0.023, 0.026, 0.029, 0.026, 0.039};
       double* aeyhstat = new double[5]  {0.023, 0.026, 0.029, 0.026, 0.039};
       double* aeylsys  = new double[5]  {0.015, 0.016, 0.016, 0.013, 0.019};
       double* aeyhsys  = new double[5]   {0.015, 0.016, 0.016, 0.013, 0.019};
    
    double by[]      = {-0.016, 0.015, 0.037, 0.094, 0.022};
    double* beylstat = new double[5]  {0.024, 0.026, 0.031, 0.026, 0.040};
    double* beyhstat = new double[5]  {0.024, 0.026, 0.031, 0.026, 0.040};
     
       TGraphMultiErrors* Fig5 = new TGraphMultiErrors("Fig5", "v_{2,Jpsi} wrt p_{T}", 5, ax, ay, aexl, aexh, aeylstat, aeyhstat);
        Fig5->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Fig5->GetYaxis()->SetTitle("v_{2,Jpsi}");
       Fig5->AddYError(5, aeylsys, aeyhsys);
       Fig5->SetMarkerStyle(20);
    Fig5->SetMarkerColor(kBlue);
       Fig5->SetLineColorAlpha(46, 0.);
       Fig5->GetAttLine(0)->SetLineColor(kBlue);
       Fig5->GetAttLine(1)->SetLineColor(kBlue);
    Fig5->GetAttFill(1)->SetFillStyle(0);
    
    TGraphMultiErrors* Fig5me = new TGraphMultiErrors("Fig5me", "v_{2,Jpsi} wrt p_{T}", 5, ax, by, aexl, aexh, beylstat, beyhstat);
        Fig5me->GetXaxis()->CenterTitle("p_{T}");
    Fig5me->GetYaxis()->CenterTitle("v_{2,Jpsi}");
       Fig5me->SetMarkerStyle(20);
    Fig5me->SetMarkerColor(kRed);
       Fig5me->SetLineColorAlpha(46, 0.);
       Fig5me->GetAttLine(0)->SetLineColor(kRed);
       
       Fig5->Draw("a p s ; ; 5 s=0.5");
    Fig5me->Draw("p s ; ; 5 s=0.5");
    
    TLine *l=new TLine(0.,0.0,8.,0.0);
       l->SetLineColor(kBlack);
       l->SetLineWidth(1);
       l->SetLineStyle(9);
       l->Draw("same");
    
    TLegend *legend5=new TLegend(0.12,0.80,0.60,0.90);
      legend5->SetFillColorAlpha(kWhite, 0.);
      legend5->SetBorderSize(0);
        legend5->SetTextFont(42);
      legend5->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"p-Pb analysis data points");
        legend5->AddEntry(Fig5,message);
    sprintf(message,"This analysis cross-check data points");
    legend5->AddEntry(Fig5me,message);
      legend5->Draw();
    
    
// Figure 4
    
    auto c4 = new TCanvas("c4","c4",200,10,600,400);
       double ax4[]      = {1.7, 2.1, 2.5, 2.85, 3.15, 3.5, 3.9, 4.3};
       double ay4[]      = {0.00278, 0.00147, 0.00183, 0.00312, 0.00375, 0.000235, 0.00157, -0.000984};
    double aexl4[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double aexh4[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
       double* aeylstat4 = new double[8]  {0.000994,0.00111,0.00128,0.00138,0.000965,0.00197,0.00246,0.00292};
       double* aeyhstat4 = new double[8]  {0.000994,0.00111,0.00128,0.00138,0.000965,0.00197,0.00246,0.00292};
       

    double by4[]      = {0.002558, 0.000947, 0.00188, 0.002813, 0.003463, -0.0002968, 0.001428, -0.0003251};
    double* beylstat4 = new double[8]  {0.001074,0.001159,0.001329,0.001413,0.001102,0.002064,0.002544,0.003025};
    double* beyhstat4 = new double[8]  {0.001074,0.001159,0.001329,0.001413,0.001102,0.002064,0.002544,0.003025};
     
       TGraphMultiErrors* Fig4 = new TGraphMultiErrors("Fig4", "V_{2,Jpsi-Tkl} wrt m_{dimuon}", 8, ax4, ay4, aexl4, aexh4, aeylstat4, aeyhstat4);
        Fig4->GetXaxis()->SetTitle("m_{dimuon} (GeV/c^{2})");
    Fig4->GetYaxis()->SetTitle("V_{2,Jpsi-Tkl}");
       Fig4->SetMarkerStyle(20);
    Fig4->SetMarkerColor(kBlue);
       Fig4->SetLineColorAlpha(46, 0.);
       Fig4->GetAttLine(0)->SetLineColor(kBlue);
    
    TGraphMultiErrors* Fig4me = new TGraphMultiErrors("Fig4me", "V_{2,Jpsi-Tkl} wrt m_{dimuon}", 8, ax4, by4, aexl4, aexh4, beylstat4, beyhstat4);
        Fig4me->GetXaxis()->SetTitle("m_{dimuon} (GeV/c^{2})");
        Fig4me->GetYaxis()->SetTitle("V_{2,Jpsi-Tkl}");
       Fig4me->SetMarkerStyle(20);
    Fig4me->SetMarkerColor(kRed);
       Fig4me->SetLineColorAlpha(46, 0.);
       Fig4me->GetAttLine(0)->SetLineColor(kRed);
       
       Fig4->Draw("a p s ; ; 5 s=0.5");
    Fig4me->Draw("p s ; ; 5 s=0.5");
    
    TLine *l4=new TLine(1.,0.0,5.,0.0);
       l4->SetLineColor(kBlack);
       l4->SetLineWidth(1);
       l4->SetLineStyle(9);
       l4->Draw("same");
    
    TLegend *legend4=new TLegend(0.12,0.80,0.60,0.90);
      legend4->SetFillColorAlpha(kWhite, 0.);
      legend4->SetBorderSize(0);
        legend4->SetTextFont(42);
      legend4->SetTextSize(0.03);
        Char_t message4[80];
        sprintf(message4,"p-Pb analysis data points");
        legend4->AddEntry(Fig4,message4);
    sprintf(message4,"This analysis cross-check data points");
    legend4->AddEntry(Fig4me,message4);
      legend4->Draw();
    
    
// Figure 3 central
    
    auto c3central = new TCanvas("c3central","c3central",200,10,600,400);
       double ax3central[]      = {0.5*3.14/6, 1.5*3.14/6, 2.5*3.14/6, 3.5*3.14/6, 4.5*3.14/6, 5.5*3.14/6};
       double ay3central[]      = {3.1434, 3.1301, 3.126, 3.1768, 3.2081, 3.2857};
    double aexl3central[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double aexh3central[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double* aeylstat3central = new double[6]  {0.011519, 0.011992, 0.011677, 0.011834, 0.011677, 0.012308};
       double* aeyhstat3central = new double[6]  {0.011519, 0.011992, 0.011677, 0.011834, 0.011677, 0.012308};
       

    double by3central[]      = {30.58, 30.48, 30.45, 30.95, 31.28, 32.05};
    double* beylstat3central = new double[6]  {0.12,0.12,0.12,0.12,0.12,0.12};
    double* beyhstat3central = new double[6]  {0.12,0.12,0.12,0.12,0.12,0.12};
    
    for(int idx=0; idx<6; idx++){
        by3central[idx]/=9.721;
        beylstat3central[idx]/=9.721;
        beyhstat3central[idx]/=9.721;
    }
     
       TGraphMultiErrors* Fig3central = new TGraphMultiErrors("Fig3central", "Central Yields wrt Delta phi", 6, ax3central, ay3central, aexl3central, aexh3central, aeylstat3central, aeyhstat3central);
        Fig3central->GetXaxis()->SetTitle("Delta phi (rad)");
    Fig3central->GetYaxis()->SetTitle("Yield_{Central} (a.u.)");
       Fig3central->SetMarkerStyle(20);
    Fig3central->SetMarkerColor(kBlue);
       Fig3central->SetLineColorAlpha(46, 0.);
       Fig3central->GetAttLine(0)->SetLineColor(kBlue);
    
    TGraphMultiErrors* Fig3centralme = new TGraphMultiErrors("Fig3centralme", "Central Yields wrt Delta phi", 6, ax3central, by3central, aexl3central, aexh3central, beylstat3central, beyhstat3central);
        Fig3centralme->GetXaxis()->SetTitle("Delta phi (rad)");
        Fig3centralme->GetYaxis()->SetTitle("Yield_{Central} (a.u.)");
       Fig3centralme->SetMarkerStyle(20);
    Fig3centralme->SetMarkerColor(kRed);
       Fig3centralme->SetLineColorAlpha(46, 0.);
       Fig3centralme->GetAttLine(0)->SetLineColor(kRed);
       
       Fig3central->Draw("a p s ; ; 5 s=0.5");
    Fig3centralme->Draw("p s ; ; 5 s=0.5");
    
    TLegend *legend3central=new TLegend(0.12,0.80,0.60,0.90);
      legend3central->SetFillColorAlpha(kWhite, 0.);
      legend3central->SetBorderSize(0);
        legend3central->SetTextFont(42);
      legend3central->SetTextSize(0.03);
        Char_t message3central[80];
        sprintf(message3central,"p-Pb analysis data points");
        legend3central->AddEntry(Fig3central,message3central);
    sprintf(message3central,"This analysis cross-check data points");
    legend3central->AddEntry(Fig3centralme,message3central);
      legend3central->Draw();
    
    // Figure 3 periph
    
    auto c3periph = new TCanvas("c3periph","c3periph",200,10,600,400);
       double ax3periph[]      = {0.5*3.14/6, 1.5*3.14/6, 2.5*3.14/6, 3.5*3.14/6, 4.5*3.14/6, 5.5*3.14/6};
       double ay3periph[]      = {1., 1.013, 1.033, 1.053, 1.091, 1.116};
    double aexl3periph[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double aexh3periph[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double* aeylstat3periph = new double[6]  {0.0075, 0.0075, 0.0079, 0.0083, 0.0083, 0.008};
       double* aeyhstat3periph = new double[6]  {0.0075, 0.0075, 0.0079, 0.0083, 0.0083, 0.008};
       

    double by3periph[]      = {9.721, 9.86, 9.994, 10.21, 10.56, 10.77};
    double* beylstat3periph = new double[6]  {0.073,0.07,0.074,0.07,0.08,0.08};
    double* beyhstat3periph = new double[6]  {0.073,0.07,0.074,0.07,0.08,0.08};
    
    for(int idx=0; idx<6; idx++){
        by3periph[idx]/=9.721;
        beylstat3periph[idx]/=9.721;
        beyhstat3periph[idx]/=9.721;
    }
     
       TGraphMultiErrors* Fig3periph = new TGraphMultiErrors("Fig3periph", "Peripheral Yields wrt Delta phi", 6, ax3periph, ay3periph, aexl3periph, aexh3periph, aeylstat3periph, aeyhstat3periph);
        Fig3periph->GetXaxis()->SetTitle("Delta phi (rad)");
    Fig3periph->GetYaxis()->SetTitle("Yield_{periph} (a.u.)");
       Fig3periph->SetMarkerStyle(20);
    Fig3periph->SetMarkerColor(kBlue);
       Fig3periph->SetLineColorAlpha(46, 0.);
       Fig3periph->GetAttLine(0)->SetLineColor(kBlue);
    
    TGraphMultiErrors* Fig3periphme = new TGraphMultiErrors("Fig3periphme", "Peripheral Yields wrt Delta phi", 6, ax3periph, by3periph, aexl3periph, aexh3periph, beylstat3periph, beyhstat3periph);
        Fig3periphme->GetXaxis()->SetTitle("Delta phi (rad)");
        Fig3periphme->GetYaxis()->SetTitle("Yield_{periph} (a.u.)");
       Fig3periphme->SetMarkerStyle(20);
    Fig3periphme->SetMarkerColor(kRed);
       Fig3periphme->SetLineColorAlpha(46, 0.);
       Fig3periphme->GetAttLine(0)->SetLineColor(kRed);
       
       Fig3periph->Draw("a p s ; ; 5 s=0.5");
    Fig3periphme->Draw("p s ; ; 5 s=0.5");
    
    TLegend *legend3periph=new TLegend(0.12,0.80,0.60,0.90);
      legend3periph->SetFillColorAlpha(kWhite, 0.);
      legend3periph->SetBorderSize(0);
        legend3periph->SetTextFont(42);
      legend3periph->SetTextSize(0.03);
        Char_t message3periph[80];
        sprintf(message3periph,"p-Pb analysis data points");
        legend3periph->AddEntry(Fig3periph,message3periph);
    sprintf(message3periph,"This analysis cross-check data points");
    legend3periph->AddEntry(Fig3periphme,message3periph);
      legend3periph->Draw();
    
// Figure 3 sub
    
    auto c3sub = new TCanvas("c3sub","c3sub",200,10,600,400);
       double ax3sub[]      = {0.5*3.14/6, 1.5*3.14/6, 2.5*3.14/6, 3.5*3.14/6, 4.5*3.14/6, 5.5*3.14/6};
       double ay3sub[]      = {1., 1.013, 1.033, 1.053, 1.091, 1.116};
    double aexl3sub[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double aexh3sub[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double* aeylstat3sub = new double[6]  {0.0075, 0.0075, 0.0079, 0.0083, 0.0083, 0.008};
       double* aeyhstat3sub = new double[6]  {0.0075, 0.0075, 0.0079, 0.0083, 0.0083, 0.008};
       

    double by3sub[]      = {9.721, 9.86, 9.994, 10.21, 10.56, 10.77};
    double* beylstat3sub = new double[6]  {0.073,0.07,0.074,0.07,0.08,0.08};
    double* beyhstat3sub = new double[6]  {0.073,0.07,0.074,0.07,0.08,0.08};
    
    for(int idx=0; idx<6; idx++){
        ay3sub[idx] = ay3central[idx] - ay3periph[idx];
        aeylstat3sub[idx] = sqrt(pow(aeylstat3central[idx],2)+pow(aeylstat3periph[idx],2));
        aeyhstat3sub[idx] = sqrt(pow(aeyhstat3central[idx],2)+pow(aeyhstat3periph[idx],2));
        by3sub[idx] = by3central[idx] - by3periph[idx];
        beylstat3sub[idx] = sqrt(pow(beylstat3central[idx],2)+pow(beylstat3periph[idx],2));
        beyhstat3sub[idx] = sqrt(pow(beyhstat3central[idx],2)+pow(beyhstat3periph[idx],2));
        
    }
     
       TGraphMultiErrors* Fig3sub = new TGraphMultiErrors("Fig3sub", "Subtracted Yields wrt Delta phi", 6, ax3sub, ay3sub, aexl3sub, aexh3sub, aeylstat3sub, aeyhstat3sub);
        Fig3sub->GetXaxis()->SetTitle("Delta phi (rad)");
    Fig3sub->GetYaxis()->SetTitle("Yield_{sub} (a.u.)");
       Fig3sub->SetMarkerStyle(20);
    Fig3sub->SetMarkerColor(kBlue);
       Fig3sub->SetLineColorAlpha(46, 0.);
       Fig3sub->GetAttLine(0)->SetLineColor(kBlue);
    
    TGraphMultiErrors* Fig3subme = new TGraphMultiErrors("Fig3subme", "Subtracted Yields wrt Delta phi", 6, ax3sub, by3sub, aexl3sub, aexh3sub, beylstat3sub, beyhstat3sub);
        Fig3subme->GetXaxis()->SetTitle("Delta phi (rad)");
        Fig3subme->GetYaxis()->SetTitle("Yield_{sub} (a.u.)");
       Fig3subme->SetMarkerStyle(20);
    Fig3subme->SetMarkerColor(kRed);
       Fig3subme->SetLineColorAlpha(46, 0.);
       Fig3subme->GetAttLine(0)->SetLineColor(kRed);
       
       Fig3sub->Draw("a p s ; ; 5 s=0.5");
    Fig3subme->Draw("p s ; ; 5 s=0.5");
    
    TLegend *legend3sub=new TLegend(0.12,0.80,0.60,0.90);
      legend3sub->SetFillColorAlpha(kWhite, 0.);
      legend3sub->SetBorderSize(0);
        legend3sub->SetTextFont(42);
      legend3sub->SetTextSize(0.03);
        Char_t message3sub[80];
        sprintf(message3sub,"p-Pb analysis data points");
        legend3sub->AddEntry(Fig3sub,message3sub);
    sprintf(message3sub,"This analysis cross-check data points");
    legend3sub->AddEntry(Fig3subme,message3sub);
      legend3sub->Draw();
    
    
// Figure 2 central
    
    auto c2central = new TCanvas("c2central","c2central",200,10,600,400);
       double ax2central[]      = {1.7, 2.1, 2.5, 2.85, 3.15, 3.5, 3.9, 4.3};
       double ay2central[]      = {2.855, 2.849, 2.867, 2.841, 2.830, 2.849, 2.86, 2.83};
    double aexl2central[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double aexh2central[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double* aeylstat2central = new double[8]  {0.00688, 0.007568, 0.008772, 0.009631, 0.007138, 0.0135, 0.0164, 0.0196};
       double* aeyhstat2central = new double[8]  {0.00688, 0.007568, 0.008772, 0.009631, 0.007138, 0.0135, 0.0164, 0.0196};
       

    double by2central[]      = {31, 30.89, 31.08, 30.84, 30.63, 30.88,0.,30.82};
    double* beylstat2central = new double[8]  {0.075,0.077,0.066,0.107,0.079,0.147,0.,0.219};
    double* beyhstat2central = new double[8]  {0.075,0.077,0.066,0.107,0.079,0.147,0.,0.219};
    
    for(int idx = 0; idx<8; idx++){
        ay2central[idx]/=2.855;
        aeylstat2central[idx]/=2.855;
        aeyhstat2central[idx]/=2.855;
        by2central[idx]/=31;
        beylstat2central[idx]/=31;
        beyhstat2central[idx]/=31;
    }
    

     
       TGraphMultiErrors* Fig2central = new TGraphMultiErrors("Fig2central", "Central Yields wrt dimuon mass", 8, ax2central, ay2central, aexl2central, aexh2central, aeylstat2central, aeyhstat2central);
        Fig2central->GetXaxis()->SetTitle("m_{dimuon} (GeV/c^{2})");
    Fig2central->GetYaxis()->SetTitle("Yield_{Central} (a.u.)");
       Fig2central->SetMarkerStyle(20);
    Fig2central->SetMarkerColor(kBlue);
       Fig2central->SetLineColorAlpha(46, 0.);
       Fig2central->GetAttLine(0)->SetLineColor(kBlue);
    
    TGraphMultiErrors* Fig2centralme = new TGraphMultiErrors("Fig2centralme", "Central Yields wrt dimuon mass", 8, ax2central, by2central, aexl2central, aexh2central, beylstat2central, beyhstat2central);
        Fig2centralme->GetXaxis()->SetTitle("m_{dimuon} (GeV/c^{2})");
        Fig2centralme->GetYaxis()->SetTitle("Yield_{Central} (a.u.)");
       Fig2centralme->SetMarkerStyle(20);
    Fig2centralme->SetMarkerColor(kRed);
       Fig2centralme->SetLineColorAlpha(46, 0.);
       Fig2centralme->GetAttLine(0)->SetLineColor(kRed);
       
       Fig2central->Draw("a p s ; ; 5 s=0.5");
    Fig2centralme->Draw("p s ; ; 5 s=0.5");
    
    TLegend *legend2central=new TLegend(0.12,0.80,0.60,0.90);
      legend2central->SetFillColorAlpha(kWhite, 0.);
      legend2central->SetBorderSize(0);
        legend2central->SetTextFont(42);
      legend2central->SetTextSize(0.03);
        Char_t message2central[80];
        sprintf(message2central,"p-Pb analysis data points");
        legend2central->AddEntry(Fig2central,message2central);
    sprintf(message2central,"This analysis cross-check data points");
    legend2central->AddEntry(Fig2centralme,message2central);
      legend2central->Draw();
    
    // Figure 2 periph
    
    auto c2periph = new TCanvas("c2periph","c2periph",200,10,600,400);
       double ax2periph[]      = {1.7, 2.1, 2.5, 2.85, 3.15, 3.5, 3.9, 4.3};
       double ay2periph[]      = {0.9288, 0.9283, 0.9349, 0.9206, 0.9178, 0.9478, 0.9466, 0.9045};
    double aexl2periph[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double aexh2periph[]      = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double* aeylstat2periph = new double[8]  {0.006192, 0.007118, 0.008044, 0.008218, 0.00544, 0.01267, 0.01563, 0.01823};
       double* aeyhstat2periph = new double[8]  {0.006192, 0.007118, 0.008044, 0.008218, 0.00544, 0.01267, 0.01563, 0.01823};
       

    double by2periph[]      = {10.01, 9.964, 10.13, 9.946, 9.907, 10.13, 10.06, 9.759};
    double* beylstat2periph = new double[8]  {0.0719,0.0801,0.0923,0.0935,0.0620,0.145,0.177,0.209};
    double* beyhstat2periph = new double[8]  {0.0719,0.0801,0.0923,0.0935,0.0620,0.145,0.177,0.209};
    
    for(int idx = 0; idx<8; idx++){
        ay2periph[idx]/=0.9288;
        aeylstat2periph[idx]/=0.9288;
        aeyhstat2periph[idx]/=0.9288;
        by2periph[idx]/=10.01;
        beylstat2periph[idx]/=10.01;
        beyhstat2periph[idx]/=10.01;
    }
     
       TGraphMultiErrors* Fig2periph = new TGraphMultiErrors("Fig2periph", "Peripheral Yields wrt dimuon mass", 8, ax2periph, ay2periph, aexl2periph, aexh2periph, aeylstat2periph, aeyhstat2periph);
        Fig2periph->GetXaxis()->SetTitle("m_{dimuon} (GeV/c^{2})");
    Fig2periph->GetYaxis()->SetTitle("Yield_{periph} (a.u.)");
       Fig2periph->SetMarkerStyle(20);
    Fig2periph->SetMarkerColor(kBlue);
       Fig2periph->SetLineColorAlpha(46, 0.);
       Fig2periph->GetAttLine(0)->SetLineColor(kBlue);
    
    TGraphMultiErrors* Fig2periphme = new TGraphMultiErrors("Fig2periphme", "Peripheral Yields wrt dimuon mass", 8, ax2periph, by2periph, aexl2periph, aexh2periph, beylstat2periph, beyhstat2periph);
        Fig2periphme->GetXaxis()->SetTitle("m_{dimuon} (GeV/c^{2})");
        Fig2periphme->GetYaxis()->SetTitle("Yield_{periph} (a.u.)");
       Fig2periphme->SetMarkerStyle(20);
    Fig2periphme->SetMarkerColor(kRed);
       Fig2periphme->SetLineColorAlpha(46, 0.);
       Fig2periphme->GetAttLine(0)->SetLineColor(kRed);
       
       Fig2periph->Draw("a p s ; ; 5 s=0.5");
    Fig2periphme->Draw("p s ; ; 5 s=0.5");
    
    TLegend *legend2periph=new TLegend(0.12,0.80,0.60,0.90);
      legend2periph->SetFillColorAlpha(kWhite, 0.);
      legend2periph->SetBorderSize(0);
        legend2periph->SetTextFont(42);
      legend2periph->SetTextSize(0.03);
        Char_t message2periph[80];
        sprintf(message2periph,"p-Pb analysis data points");
        legend2periph->AddEntry(Fig2periph,message2periph);
    sprintf(message2periph,"This analysis cross-check data points");
    legend2periph->AddEntry(Fig2periphme,message2periph);
      legend2periph->Draw();
    
    
}
