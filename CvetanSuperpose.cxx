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
    
    
    
    {
    
    // Figure Combine systems

       auto csys = new TCanvas("csys","csys",200,10,800,600);
       double qax[]      = {1, 2.5, 3.5, 5, 7};
       double qay[]      = {0,-0.015,0.031,0.082,0.033};
       double qaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
       double qaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
       double* qaeylstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
       double* qaeyhstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
       double* qaeylsys  = new double[5]  {0.008,0.011,0.013,0.011,0.017};
       double* qaeyhsys  = new double[5]   {0.008,0.011,0.013,0.011,0.017};
    
    double qbx[]      = {0.64, 1.49, 2.47, 3.46, 4.45, 5.45, 6.819};
    double qby[]      = {0.011, 0.043, 0.074, 0.088,0.085,0.103, 0.083};
        double qbexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        double qbexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double* qbeylstat = new double[7]  {0.0085, 0.0069, 0.0069, 0.0077, 0.009, 0.011, 0.011};
    double* qbeyhstat = new double[7]  {0.0085, 0.0069, 0.0069, 0.0077, 0.009, 0.011, 0.011};
    double* qbeylsys  = new double[7]  {0.0038391, 0.0036633, 0.004898, 0.0035068, 0.0037855,0.0029726,  0.0036802};
    double* qbeyhsys  = new double[7]   {0.0038391, 0.0036633, 0.004898, 0.0035068, 0.0037855,0.0029726,  0.0036802};
    
    double qcx[]      = {1, 2.5, 3.5, 5, 7};
    double qcy[]      = {-0.00584, 0.01802, -0.01266, 0.03064, -0.05892};
    double* qceylstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503};
    double* qceyhstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503};
   // double* qceylsys  = new double[5]  {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
  //  double* qceyhsys  = new double[5]   {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
double* qceylsys  = new double[5]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258};
double* qceyhsys  = new double[5]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258};    
    
        for(int idx=0; idx<5; idx++){
            qceylsys[idx] = abs(qcy[idx])*sqrt(pow(qceylsys[idx]/qcy[idx],2)+pow(5.9/100.,2));
            qceyhsys[idx] = qceylsys[idx];
        }
    
    //0.0550667 0.000573161
     
       TGraphMultiErrors* Figsys = new TGraphMultiErrors("Figsys", "#it{v}_{2,J/#psi} wrt #it{p}_{T} - Systems", 5, qax, qay, qaexl, qaexh, qaeylstat, qaeyhstat);
        Figsys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    Figsys->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
       Figsys->AddYError(5, qaeylsys, qaeyhsys);
       Figsys->SetMarkerStyle(24);
    Figsys->SetMarkerColor(kGray+2);
       Figsys->SetLineColorAlpha(46, 0.);
       Figsys->GetAttLine(0)->SetLineColor(kGray+2);
       Figsys->GetAttLine(1)->SetLineColor(kGray+2);
        Figsys->GetAttLine(0)->SetLineWidth(2);
        Figsys->GetAttLine(1)->SetLineWidth(2);
    Figsys->GetAttFill(1)->SetFillStyle(0);
    
    TGraphMultiErrors* Fig5PbPb = new TGraphMultiErrors("Fig5PbPb", "#it{#v}_{2,J/#psi} wrt #it{p}_{T}", 7, qbx, qby, qbexl, qbexh, qbeylstat, qbeyhstat);
        Fig5PbPb->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
    Fig5PbPb->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
       Fig5PbPb->AddYError(7, qbeylsys, qbeyhsys);
          Fig5PbPb->SetMarkerStyle(21);
       Fig5PbPb->SetMarkerColor(kBlack);
          Fig5PbPb->SetLineColorAlpha(46, 0.);
          Fig5PbPb->GetAttLine(0)->SetLineColor(kBlack);
          Fig5PbPb->GetAttLine(1)->SetLineColor(kBlack);
        Fig5PbPb->GetAttLine(0)->SetLineWidth(2);
        Fig5PbPb->GetAttLine(1)->SetLineWidth(2);
       Fig5PbPb->GetAttFill(1)->SetFillStyle(0);
        
    TGraphMultiErrors* Fig5pp = new TGraphMultiErrors("Fig5pp", "#it{#v}_{2,J/#psi} wrt #it{p}_{T}", 5, qcx, qcy, qaexl, qaexh, qceylstat, qceyhstat);
        Fig5pp->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
    Fig5pp->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
    Fig5pp->AddYError(5, qceylsys, qceyhsys);
       Fig5pp->SetMarkerStyle(20);
    Fig5pp->SetMarkerColor(kAzure-3);
       Fig5pp->SetLineColorAlpha(46, 0.);
       Fig5pp->GetAttLine(0)->SetLineColor(kAzure-3);
    Fig5pp->GetAttLine(1)->SetLineColor(kAzure-3);
        Fig5pp->GetAttLine(0)->SetLineWidth(2);
        Fig5pp->GetAttLine(1)->SetLineWidth(2);
     Fig5pp->GetAttFill(1)->SetFillStyle(0);
       
       Figsys->Draw("a p s ; z ; 5 s=0.5");
    Fig5PbPb->Draw("p s ; z ; 5 s=0.5");
    Fig5pp->Draw("p s ; z ; 5 s=0.5");
    
    TLine *ls=new TLine(0.,0.0,8.,0.0);
       ls->SetLineColor(kBlack);
       ls->SetLineWidth(1);
       ls->SetLineStyle(9);
       ls->Draw("same");
        
        TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
                  legendov2->SetFillColorAlpha(kWhite, 0.);
                  legendov2->SetBorderSize(0);
                   legendov2->SetTextFont(42);
                   legendov2->SetTextSize(0.05);
        Char_t messago[80];
        sprintf(messago,"ALICE Preliminary");
         legendov2->AddEntry(Fig5pp,messago,"");
                  legendov2->Draw();
        
        TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
                     legendov4->SetFillColorAlpha(kWhite, 0.);
                     legendov4->SetBorderSize(0);
                      legendov4->SetTextFont(62);
                      legendov4->SetTextSize(0.04);
        legendov4->SetTextColor(kAzure-3);
           sprintf(messago,"5.9%% global syst. uncertainty");
            legendov4->AddEntry(Fig5pp,messago,"");
                     legendov4->Draw();
    
    TLegend *legendsys=new TLegend(0.12,0.80,0.60,0.90);
      legendsys->SetFillColorAlpha(kWhite, 0.);
      legendsys->SetBorderSize(0);
        legendsys->SetTextFont(42);
      legendsys->SetTextSize(0.03);
    sprintf(message,"Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, (30-50%%) (JHEP 10 (2020) 141)");
    legendsys->AddEntry(Fig5PbPb,message);
        sprintf(message,"2.5 < #it{y}_{cms} < 4.0");
        legendsys->AddEntry(Fig5PbPb,message,"");
        sprintf(message,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02, 8.16 TeV, (0-20%%)-(40-100%%) (PLB 780 (2018) 7-20)");
        legendsys->AddEntry(Figsys,message);
        sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.03 < #it{y}_{cms} < 3.53");
        legendsys->AddEntry(Figsys,message,"");
    sprintf(message,"pp, #sqrt{#it{s}_{NN}} = 13 TeV, (0-5%%)-(40-100%%)");
    legendsys->AddEntry(Fig5pp,message);
        sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.5 < #it{y}_{cms} < 4.0");
        legendsys->AddEntry(Fig5pp,message,"");
      legendsys->Draw();
    
    
}
    
    
    {
        
        // Figure Combine systems

           auto csys = new TCanvas("csys2","csys2",200,10,800,600);
           double qax[]      = {1, 2.5, 3.5, 5, 7};
           double qay[]      = {0,-0.015,0.031,0.082,0.033};
           double qaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
           double qaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
           double* qaeylstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
           double* qaeyhstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
           double* qaeylsys  = new double[5]  {0.008,0.011,0.013,0.011,0.017};
           double* qaeyhsys  = new double[5]   {0.008,0.011,0.013,0.011,0.017};
        
        for(int idx=0; idx<5; idx++){
            qay[idx]/=0.068;
            qaeylstat[idx]/=0.068;
            qaeyhstat[idx]/=0.068;
            qaeylsys[idx]/=0.068;
            qaeyhsys[idx]/=0.068;
        }
        
        double qbx[]      = {1, 2.5, 3.5, 5, 7};
        double qby[]      = {0.023,0.016,0.074,0.069,0.079};
        double* qbeylstat = new double[5]  {0.015, 0.018, 0.022, 0.020, 0.034};
        double* qbeyhstat = new double[5]  {0.015, 0.018, 0.022, 0.020, 0.034};
        double* qbeylsys  = new double[5]  {0.014, 0.014, 0.015, 0.015, 0.022};
        double* qbeyhsys  = new double[5]   {0.014, 0.014, 0.015, 0.015, 0.022};
        
        for(int idx=0; idx<5; idx++){
            qby[idx]/=0.068;
            qbeylstat[idx]/=0.068;
            qbeyhstat[idx]/=0.068;
            qbeylsys[idx]/=0.068;
            qbeyhsys[idx]/=0.068;
        }
        
        double qcx[]      = {1, 2.5, 3.5, 5, 7};
        double oldqcy[]      = {-0.00584, 0.01802, -0.01266, 0.03064, -0.05892}; //10mrad
        double qcy[]  {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684};
        double* qceylstat = new double[5]  {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613};
                double* qceyhstat = new double[5]  {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613};
      //  double* qceylsys  = new double[5]  {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
    //    double* qceyhsys  = new double[5]   {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
double* qceylsys  = new double[5]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258};
        double* qceyhsys  = new double[5]   {0.0099, 0.0077, 0.0140, 0.0085, 0.0258};
//       double* qceylstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503};
//        double* qceyhstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503}; 10mrad
//        double* qceylsys  = new double[5]  {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
//        double* qceyhsys  = new double[5]   {0.0099, 0.0137, 0.0168, 0.0189, 0.0380}; 10 mrad
//        {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684}; 5mr
//        {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613}; 5mr

        for(int idx=0; idx<5; idx++){
            qcy[idx]/=0.066;
            qceylstat[idx]/=0.066;
            qceyhstat[idx]/=0.066;
            qceylsys[idx]/=0.066;
            qceylsys[idx]/=abs(qcy[idx]*0.066/oldqcy[idx]);
            qceyhsys[idx]/=0.066;
            qceyhsys[idx]/=abs(qcy[idx]*0.066/oldqcy[idx]);
        }
        
        //0.0550667 0.000573161
         
           TGraphMultiErrors* Figsys = new TGraphMultiErrors("Figsys", "#it{#epsilon}_{2,J/#psi} wrt #it{p}_{T} - Systems", 5, qax, qay, qaexl, qaexh, qaeylstat, qaeyhstat);
            Figsys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        Figsys->GetYaxis()->SetTitle("#it{#epsilon}_{2,J/#psi}");
           Figsys->AddYError(5, qaeylsys, qaeyhsys);
           Figsys->SetMarkerStyle(24);
        Figsys->SetMarkerColor(kGray+2);
           Figsys->SetLineColorAlpha(46, 0.);
           Figsys->GetAttLine(0)->SetLineColor(kGray+2);
           Figsys->GetAttLine(1)->SetLineColor(kGray+2);
            Figsys->GetAttLine(0)->SetLineWidth(2);
            Figsys->GetAttLine(1)->SetLineWidth(2);
        Figsys->GetAttFill(1)->SetFillStyle(0);
        
        TGraphMultiErrors* Fig5Pbp = new TGraphMultiErrors("Fig5PbPb", "#it{#epsilon}_{2,J/#psi} wrt #it{p}_{T}", 5, qbx, qby, qaexl, qaexh, qbeylstat, qbeyhstat);
            Fig5Pbp->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
        Fig5Pbp->GetYaxis()->CenterTitle("#it{#epsilon}_{2,J/#psi}");
           Fig5Pbp->AddYError(5, qbeylsys, qbeyhsys);
              Fig5Pbp->SetMarkerStyle(24);
           Fig5Pbp->SetMarkerColor(kBlack);
              Fig5Pbp->SetLineColorAlpha(46, 0.);
              Fig5Pbp->GetAttLine(0)->SetLineColor(kBlack);
              Fig5Pbp->GetAttLine(1)->SetLineColor(kBlack);
            Fig5Pbp->GetAttLine(0)->SetLineWidth(2);
            Fig5Pbp->GetAttLine(1)->SetLineWidth(2);
           Fig5Pbp->GetAttFill(1)->SetFillStyle(0);
            
        TGraphMultiErrors* Fig5pp = new TGraphMultiErrors("Fig5pp", "#it{#epsilon}_{2,J/#psi} wrt #it{p}_{T}", 5, qcx, qcy, qaexl, qaexh, qceylstat, qceyhstat);
            Fig5pp->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
        Fig5pp->GetYaxis()->CenterTitle("#it{#epsilon}_{2,J/#psi}");
        Fig5pp->AddYError(5, qceylsys, qceyhsys);
           Fig5pp->SetMarkerStyle(20);
        Fig5pp->SetMarkerColor(kAzure-3);
           Fig5pp->SetLineColorAlpha(46, 0.);
           Fig5pp->GetAttLine(0)->SetLineColor(kAzure-3);
        Fig5pp->GetAttLine(1)->SetLineColor(kAzure-3);
            Fig5pp->GetAttLine(0)->SetLineWidth(2);
            Fig5pp->GetAttLine(1)->SetLineWidth(2);
         Fig5pp->GetAttFill(1)->SetFillStyle(0);
           
           Figsys->Draw("a p s ; z ; 5 s=0.5");
        Fig5Pbp->Draw("p s ; z ; 5 s=0.5");
        Fig5pp->Draw("p s ; z ; 5 s=0.5");
        
        TLine *ls=new TLine(0.,0.0,8.,0.0);
           ls->SetLineColor(kBlack);
           ls->SetLineWidth(1);
           ls->SetLineStyle(9);
           ls->Draw("same");
            
            TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
                      legendov2->SetFillColorAlpha(kWhite, 0.);
                      legendov2->SetBorderSize(0);
                       legendov2->SetTextFont(42);
                       legendov2->SetTextSize(0.05);
            Char_t messago[80];
            sprintf(messago,"ALICE Preliminary");
             legendov2->AddEntry(Fig5pp,messago,"");
                      legendov2->Draw();
            
//            TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
//                         legendov4->SetFillColorAlpha(kWhite, 0.);
//                         legendov4->SetBorderSize(0);
//                          legendov4->SetTextFont(62);
//                          legendov4->SetTextSize(0.04);
//            legendov4->SetTextColor(kAzure-3);
//               sprintf(messago,"5.9%% global syst. uncertainty");
//                legendov4->AddEntry(Fig5pp,messago,"");
//                         legendov4->Draw();
        
        TLegend *legendsys=new TLegend(0.12,0.80,0.60,0.90);
          legendsys->SetFillColorAlpha(kWhite, 0.);
          legendsys->SetBorderSize(0);
            legendsys->SetTextFont(42);
          legendsys->SetTextSize(0.03);
         sprintf(message,"Pb-p, #sqrt{#it{s}_{NN}} = 5.02, 8.16 TeV, (0-20%%)-(40-100%%) (PLB 780 (2018) 7-20)");
                   legendsys->AddEntry(Fig5Pbp,message);
                   sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, -4.46 < #it{y}_{cms} < 2.96");
               legendsys->AddEntry(Fig5Pbp,message,"");
            sprintf(message,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02, 8.16 TeV, (0-20%%)-(40-100%%) (PLB 780 (2018) 7-20)");
            legendsys->AddEntry(Figsys,message);
            sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.03 < #it{y}_{cms} < 3.53");
        legendsys->AddEntry(Figsys,message,"");
        sprintf(message,"pp, #sqrt{#it{s}_{NN}} = 13 TeV, (0-5%%)-(40-100%%)");
        legendsys->AddEntry(Fig5pp,message);
            sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.5 < #it{y}_{cms} < 4.0");
            legendsys->AddEntry(Fig5pp,message,"");
          legendsys->Draw();
        
        
    }
       
       // Final v2 V0M
    {
    
    double yax[]      = {1, 2.5, 3.5, 5, 7, 10, 12};
    double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
    double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
    double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
    double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

          auto cv0m = new TCanvas("cv0m","cv0m",200,10,800,600);
          double vax[]      = {1, 2.5, 3.5, 5, 7, 10};
          double vay[]      = {-0.00584,0.01802,-0.01266,0.03064,-0.05892,0.03671};
          double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
          double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
          double* vaeylstat = new double[6]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503,0.02098};
          double* vaeyhstat = new double[6]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503,0.02098};
         // double* vaeylsys  = new double[6]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258,0.0234};
        double* vaeylsys  = new double[6] {0.0093, 0.0077, 0.0149, 0.0094, 0.0621, 0.1097};
          double* vaeyhsys  = new double[6] {0.0093, 0.0077, 0.0149, 0.0094, 0.0621, 0.1097};
       
       double vbx[]      = {1.2, 2.7, 3.7, 5.2, 7.2, 10.2};
       double vby[]      = {-0.00441, 0.01912, -0.01151, 0.03532,-0.04245,0.04257};
       double* vbeylstat = new double[6]  {0.00820, 0.00972, 0.01101, 0.00962,0.01472,0.02161};
       double* vbeyhstat = new double[6]  {0.00820, 0.00972, 0.01101, 0.00962,0.01472,0.02161};
       //double* vbeylsys  = new double[6]  {0.0093, 0.0074, 0.0136, 0.0104, 0.0245,0.0277};
        double* vbeylsys  = new double[6]   {0.0089, 0.0071, 0.0146, 0.0107, 0.0575,0.1666};
       double* vbeyhsys  = new double[6]   {0.0089, 0.0071, 0.0146, 0.0107, 0.0575,0.1666};
    
    TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
           Figempty->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
       Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
        Figempty->GetYaxis()->SetRangeUser(-0.15, 0.20);
        Figempty->GetXaxis()->SetRangeUser(0, 12.);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
          Figempty->AddYError(7, yaeylsys, yaeyhsys);
    Figempty->GetXaxis()->SetRangeUser(0,15);
        
          TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "", 6, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);//#it{v}_{2,J/#psi} wrt #it{p}_{T}
           Figv0malice->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
       Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
        Figv0malice->GetYaxis()->SetRangeUser(-0.15, 0.20);
        Figv0malice->GetXaxis()->SetRangeUser(0, 12.);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetPadBottomMargin(0.15);
        gStyle->SetPadLeftMargin(0.15);
          Figv0malice->AddYError(6, vaeylsys, vaeyhsys);
          Figv0malice->SetMarkerStyle(8);
       Figv0malice->SetMarkerColor(kAzure-3);
          Figv0malice->SetLineColorAlpha(46, 0.);
          Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
          Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
    Figv0malice->GetAttLine(0)->SetLineWidth(2);
    Figv0malice->GetAttLine(1)->SetLineWidth(2);
       Figv0malice->GetAttFill(1)->SetFillStyle(0);
    Figv0malice->GetXaxis()->SetRangeUser(0,15);
       
       TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "", 6, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
           Figv0matlas->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
       Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
        Figv0matlas->GetYaxis()->SetRangeUser(-0.15, 0.20);
        Figv0matlas->GetXaxis()->SetRangeUser(0, 12.);
        gStyle->SetPadTickX(0);
        gStyle->SetPadTickY(0);
          Figv0matlas->AddYError(6, vbeylsys, vbeyhsys);
             Figv0matlas->SetMarkerStyle(21);
        Figv0matlas->SetMarkerColor(kRed+1);
             Figv0matlas->SetLineColorAlpha(46, 0.);
             Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
             Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
    Figv0matlas->GetAttLine(0)->SetLineWidth(2);
    Figv0matlas->GetAttLine(1)->SetLineWidth(2);
          Figv0matlas->GetAttFill(1)->SetFillStyle(0);
       
          Figempty->Draw("a p s ; z ; 5 s=0.5");
          Figv0malice->Draw("p s ; z ; 5 s=0.5");
       Figv0matlas->Draw("p s ; z ; 5 s=0.5");
       
       TLine *ls2=new TLine(0.,0.0,12.,0.0);
          ls2->SetLineColor(kBlack);
          ls2->SetLineWidth(1);
          ls2->SetLineStyle(9);
          ls2->Draw("same");
       
       TLegend *legendov=new TLegend(0.12,0.65,0.40,0.80);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.035);
           Char_t messagi[80];
          // sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
         //  legendov->AddEntry(Figv0malice,messagi,"");
                   legendov->AddEntry(Figv0malice,"Subtracted yields method");

                    legendov->AddEntry(Figv0matlas,"Template fit method");

                     legendov->Draw();
    
           TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
                         legendov2->SetFillColorAlpha(kWhite, 0.);
                         legendov2->SetBorderSize(0);
                          legendov2->SetTextFont(42);
                          legendov2->SetTextSize(0.05);
               Char_t messago[80];
               sprintf(messago,"ALICE Preliminary");
                legendov2->AddEntry(Figv0matlas,messago,"");
                         legendov2->Draw();
           
           TLegend *legendov3=new TLegend(0.5,0.67,0.80,0.87);
                     legendov3->SetFillColorAlpha(kWhite, 0.);
                     legendov3->SetBorderSize(0);
                      legendov3->SetTextFont(42);
                      legendov3->SetTextSize(0.04);
           sprintf(messago,"V0M (0-5%%)-(40-100%%)");
            legendov3->AddEntry(Figv0matlas,messago,"");
           sprintf(messago,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
           legendov3->AddEntry(Figv0matlas,messago,"");
            sprintf(messago,"2.5 < #it{y}_{cms} < 4.0");
            legendov3->AddEntry(Figv0matlas,messago,"");
           sprintf(messago,"1.5 < |#it{#Delta#eta}| < 5.0");
           legendov3->AddEntry(Figv0matlas,messago,"");
                     legendov3->Draw();
           
           TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
                     legendov4->SetFillColorAlpha(kWhite, 0.);
                     legendov4->SetBorderSize(0);
                      legendov4->SetTextFont(62);
                      legendov4->SetTextSize(0.04);
        legendov4->SetTextColor(kAzure-3);
           sprintf(messago,"5.9%% global syst. uncertainty");
            legendov4->AddEntry(Figv0matlas,messago,"");
                     legendov4->Draw();
        
        TLegend *legendov42=new TLegend(0.08,0.12,0.40,0.17);
                  legendov42->SetFillColorAlpha(kWhite, 0.);
                  legendov42->SetBorderSize(0);
                   legendov42->SetTextFont(62);
                   legendov42->SetTextSize(0.04);
        legendov42->SetTextColor(kRed+1);
        sprintf(messago,"2.7%% global syst. uncertainty");
         legendov42->AddEntry(Figv0matlas,messago,"");
                  legendov42->Draw();
    
    
    }
    
    // Final v2 SPDT
    {
    
    double yax[]      = {1, 2.5, 3.5, 5, 7, 10, 12};
    double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
    double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
    double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
    double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

          auto cspdt = new TCanvas("cspdt","cspdt",200,10,600,400);
          double vax[]      = {1, 2.5, 3.5, 5, 7, 10};
          double vay[]      = {-0.01032,0.01880,-0.00154,0.04516,0.02240,0.13548};
          double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
          double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
          double* vaeylstat = new double[6]  {0.00644, 0.00744, 0.00825, 0.00728, 0.01110,0.01487};
          double* vaeyhstat = new double[6]  {0.00644, 0.00744, 0.00825, 0.00728, 0.01110,0.01487};
          //double* vaeylsys  = new double[6]  {0.0087,0.0042,0.0120,0.0058,0.0191,0.0324};
        double* vaeylsys  = new double[6]  {0.0081,0.0037,0.0108,0.0058,0.0153,0.0195};
          double* vaeyhsys  = new double[6] {0.0081,0.0037,0.0108,0.0058,0.0153,0.0195};
       
       double vbx[]      = {1.2, 2.7, 3.7, 5.2, 7.2, 10.2};
       double vby[]      = {-0.00876, 0.00805, -0.01061, 0.01805,-0.01715,0.05059};
       double* vbeylstat = new double[6]  {0.00488, 0.00596, 0.00698, 0.00647,0.01044,0.01645};
       double* vbeyhstat = new double[6]  {0.00488, 0.00596, 0.00698, 0.00647,0.01044,0.01645};
      // double* vbeylsys  = new double[6]  {0.0069, 0.0038, 0.0091, 0.0057,0.0145,0.0412};
        double* vbeylsys  = new double[6]  {0.0064, 0.0037, 0.0086, 0.0059,0.0129,0.0275};
       double* vbeyhsys  = new double[6]  {0.0064, 0.0037, 0.0086, 0.0059,0.0129,0.0275};
    
    TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
           Figempty->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
       Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
          Figempty->AddYError(7, yaeylsys, yaeyhsys);
    Figempty->GetXaxis()->SetRangeUser(0,15);
        
          TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 6, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);
           Figv0malice->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
       Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
          Figv0malice->AddYError(6, vaeylsys, vaeyhsys);
          Figv0malice->SetMarkerStyle(8);
       Figv0malice->SetMarkerColor(kAzure-3);
          Figv0malice->SetLineColorAlpha(46, 0.);
          Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
          Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
    Figv0malice->GetAttLine(0)->SetLineWidth(2);
    Figv0malice->GetAttLine(1)->SetLineWidth(2);
       Figv0malice->GetAttFill(1)->SetFillStyle(0);
    Figv0malice->GetXaxis()->SetRangeUser(0,15);
       
       TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 6, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
           Figv0matlas->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/c)");
       Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
          Figv0matlas->AddYError(6, vbeylsys, vbeyhsys);
             Figv0matlas->SetMarkerStyle(21);
        Figv0matlas->SetMarkerColor(kRed+1);
             Figv0matlas->SetLineColorAlpha(46, 0.);
             Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
             Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
    Figv0matlas->GetAttLine(0)->SetLineWidth(2);
    Figv0matlas->GetAttLine(1)->SetLineWidth(2);
          Figv0matlas->GetAttFill(1)->SetFillStyle(0);
       
          Figempty->Draw("a p s ; ; 5 s=0.5");
          Figv0malice->Draw("p s ; ; 5 s=0.5");
       Figv0matlas->Draw("p s ; ; 5 s=0.5");
       
       TLine *ls2=new TLine(0.,0.0,12.,0.0);
          ls2->SetLineColor(kBlack);
          ls2->SetLineWidth(1);
          ls2->SetLineStyle(9);
          ls2->Draw("same");
       
       TLegend *legendov=new TLegend(0.12,0.70,0.40,0.90);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.022);
           Char_t messagi[80];
           sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
           legendov->AddEntry(Figv0malice,messagi,"");
                   legendov->AddEntry(Figv0malice,"Subtracted yields default: Y_{C}-Y_{P}=a_{0}+2a_{1}cos(#Delta#phi)+2a_{2}cos(2#Delta#phi)");

                    legendov->AddEntry(Figv0matlas,"ATLAS Template fit: Y_{C}-F(Y_{P}-Y_{P}(0))=G(1+#it{v}_{2,2}cos(2#Delta#phi))");

                     legendov->Draw();
    
           TLegend *legendov2=new TLegend(0.12,0.60,0.40,0.70);
                         legendov2->SetFillColorAlpha(kWhite, 0.);
                         legendov2->SetBorderSize(0);
                          legendov2->SetTextFont(72);
                          legendov2->SetTextSize(0.03);
               Char_t messago[80];
               sprintf(messago,"ALICE Requested for preliminary");
                legendov2->AddEntry(Figv0matlas,messago,"");
                         legendov2->Draw();
           
           TLegend *legendov3=new TLegend(0.6,0.70,0.90,0.90);
                     legendov3->SetFillColorAlpha(kWhite, 0.);
                     legendov3->SetBorderSize(0);
                      legendov3->SetTextFont(62);
                      legendov3->SetTextSize(0.03);
           sprintf(messago,"SPDTracklets (0-5%)-(40-100%)");
            legendov3->AddEntry(Figv0matlas,messago,"");
           sprintf(messago,"ALICE pp, #sqrt{s_{NN}}=13 TeV");
           legendov3->AddEntry(Figv0matlas,messago,"");
           sprintf(messago,"1.5<|#Delta#eta|<5.0");
           legendov3->AddEntry(Figv0matlas,messago,"");
                     legendov3->Draw();
           
           TLegend *legendov4=new TLegend(0.45,0.10,0.85,0.20);
                     legendov4->SetFillColorAlpha(kWhite, 0.);
                     legendov4->SetBorderSize(0);
                      legendov4->SetTextFont(62);
                      legendov4->SetTextSize(0.03);
                        legendov4->SetTextColor(kAzure-3);
           sprintf(messago,"7.9 percent global syst. uncertainty");
            legendov4->AddEntry(Figv0matlas,messago,"");
                     legendov4->Draw();
        
        TLegend *legendov42=new TLegend(0.45,0.00,0.85,0.10);
                  legendov42->SetFillColorAlpha(kWhite, 0.);
                  legendov42->SetBorderSize(0);
                   legendov42->SetTextFont(62);
                   legendov42->SetTextSize(0.03);
                     legendov42->SetTextColor(kRed+1);
        sprintf(messago,"10.9 percent global syst. uncertainty");
         legendov42->AddEntry(Figv0matlas,messago,"");
                  legendov42->Draw();
    
    
    }
    
    
    double pp[] = {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684};
    double ppstat[] = {0.00842,0.00965,0.01082, 0.00967,0.01503};
   // double ppsyst[] = {0.0099, 0.0077, 0.0140, 0.0085, 0.0258};
    double ppsyst[] = {0.0093, 0.0077, 0.0149, 0.0094, 0.0621};
    
    for(int idx=0; idx<5; idx++){
        pp[idx]/=0.066;
        ppstat[idx]/=0.066;
        ppsyst[idx]/=0.066;
    }
    
    double pPb[] = {-0.000504,-0.014926,0.031617, 0.082274, 0.033703};
    double pPbstat[] = {0.018020,0.020871,0.023864, 0.021273,0.032535};
    double pPbsyst[] = {0.007717, 0.011314, 0.013313, 0.011523, 0.017635};
    
    for(int idx=0; idx<5; idx++){
        pPb[idx]/=0.068;
        pPbstat[idx]/=0.068;
        pPbsyst[idx]/=0.068;
    }
    
    double Pbp[] = {0.022974,0.016185,0.073875, 0.069017, 0.078325};
    double Pbpstat[] = {0.015113,0.018057,0.021716, 0.020426,0.034275};
    double Pbpsyst[] = {0.014058, 0.013872, 0.014531, 0.014579, 0.021298};
    
    for(int idx=0; idx<5; idx++){
        Pbp[idx]/=0.068;
        Pbpstat[idx]/=0.068;
        Pbpsyst[idx]/=0.068;
    }
    
    double avgmedium[] = {0.022974,0.016185,0.073875, 0.069017, 0.078325};
    double avgmediumstat[] = {0.015113,0.018057,0.021716, 0.020426,0.034275};
    double avgmediumsyst[] = {0.014058, 0.013872, 0.014531, 0.014579, 0.021298};
    
    for(int idx=0; idx<5; idx++){
        avgmedium[idx]=(Pbp[idx]+pPb[idx])/2.;
        avgmediumstat[idx] = sqrt(pow(Pbpstat[idx],2)+pow(pPbstat[idx],2))/2.;
        avgmediumsyst[idx] = sqrt(pow(Pbpsyst[idx],2)+pow(pPbsyst[idx],2))/2.;
    }
    
    double chi2=0;
    
    for(int idx=0; idx<5; idx++){
        chi2+=pow((pPb[idx]-pp[idx]),2)/(pow(ppstat[idx],2)+pow(pPbstat[idx],2)+pow(ppsyst[idx],2)+pow(pPbsyst[idx],2));
    }
    cout << "chi2 pp vs p-Pb: " << chi2 <<endl;
  //  cout << "Sigma deviation: " << sqrt(chi2/5) <<endl;
    cout << "Probability compatibility: " << TMath::Prob(chi2,5)/2.<<endl;
    cout << "Sigma deviation: " << 1.42*TMath::ErfInverse(1-(TMath::Prob(chi2,5))) <<endl;
    
    chi2=0;
    
    for(int idx=0; idx<5; idx++){
        chi2+=pow((Pbp[idx]-pp[idx]),2)/(pow(ppstat[idx],2)+pow(Pbpstat[idx],2)+pow(ppsyst[idx],2)+pow(Pbpsyst[idx],2));
    }
    cout << "chi2 pp vs Pb-p: " << chi2 <<endl;
   // cout << "Sigma deviation: " << sqrt(chi2/5) <<endl;
    cout << "Probability compatibility: " << TMath::Prob(chi2,5)/2.<<endl;
    cout << "Sigma deviation: " << 1.42*TMath::ErfInverse(1-(TMath::Prob(chi2,5))) <<endl;
    
    chi2=0;
    
    for(int idx=0; idx<5; idx++){
        chi2+=pow((avgmedium[idx]-pp[idx]),2)/(pow(ppstat[idx],2)+pow(avgmediumstat[idx],2)+pow(ppsyst[idx],2)+pow(avgmediumsyst[idx],2));
    }
    cout << "chi2 pp vs Pb-p+p-Pb avg.: " << chi2 <<endl;
   // cout << "Sigma deviation: " << sqrt(chi2/5) <<endl;
    cout << "Probability compatibility: " << TMath::Prob(chi2,5)/2.<<endl;
    cout << "Sigma deviation: " << 1.42*TMath::ErfInverse(1-(TMath::Prob(chi2,5))) <<endl;
    
    
    // v2 TKL wrt Centrality V0M
    {
    
    double yax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
    double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
    double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
    double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
    double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
    double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

          auto cv0m = new TCanvas("cv0m","cv0m",200,10,800,600);
          double vax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
          double vay[]      = {0.05314, 0.05250, 0.05122, 0.04958, 0.04801, 0.04739, 0.04210};
          double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
          double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
          double* vaeylstat = new double[7]  {0.00061, 0.0007, 0.00087, 0.00084, 0.00109, 0.00134, 0.00156};
          double* vaeyhstat = new double[7]  {0.00061, 0.0007, 0.00087, 0.00084, 0.00109, 0.00134, 0.00156};
         // double* vaeylsys  = new double[6]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258,0.0234};
//        double* vaeylsys  = new double[7] {0.0119, 0.0148, 0.0175, 0.0195, 0.0230, 0.0270, 0.0531}; PYTHIA Full
//          double* vaeyhsys  = new double[7] {0.0119, 0.0148, 0.0175, 0.0195, 0.0230, 0.0270, 0.0531};
        
                double* vaeylsys  = new double[7] {0.0059, 0.0071, 0.0082, 0.0090, 0.0101, 0.0112, 0.0137}; //PYTHIA Half top bottom
                  double* vaeyhsys  = new double[7] {0.0054, 0.0064, 0.0072, 0.0078, 0.0086, 0.0091, 0.0105};
       
       double vbx[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
       double vby[]      = {0.06229, 0.06460, 0.06496, 0.06543, 0.06719, 0.07017, 0.06855};
       double* vbeylstat = new double[7]  {0.00087, 0.00107, 0.00135, 0.00140, 0.00188, 0.00236, 0.00313};
       double* vbeyhstat = new double[7]  {0.00087, 0.00107, 0.00135, 0.00140, 0.00188, 0.00236, 0.00313};
       //double* vbeylsys  = new double[6]  {0.0093, 0.0074, 0.0136, 0.0104, 0.0245,0.0277};
//        double* vbeylsys  = new double[7]   {0.0090, 0.0122, 0.0143, 0.0145, 0.0180, 0.0220, 0.0295}; PYTHIA Full
//       double* vbeyhsys  = new double[7]   {0.0090, 0.0122, 0.0143, 0.0145, 0.0180, 0.0220, 0.0295};
        
                double* vbeylsys  = new double[7]   {0.0044, 0.0058, 0.0068, 0.0070, 0.0085, 0.0101, 0.0128}; //PYTHIA Half top bottom
               double* vbeyhsys  = new double[7]   {0.0041, 0.0054, 0.0062, 0.0064, 0.0077, 0.0088, 0.0109};
    
    TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,tkl} wrt Ceentrality", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
           Figempty->GetXaxis()->SetTitle("Centrality - V0M (Percent)");
       Figempty->GetYaxis()->SetTitle("#it{v}_{2,tkl}");
        Figempty->GetYaxis()->SetRangeUser(-0.15, 0.20);
        Figempty->GetXaxis()->SetRangeUser(0, 30.);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
          Figempty->AddYError(7, yaeylsys, yaeyhsys);
    Figempty->GetXaxis()->SetRangeUser(0,30);
        
          TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "", 7, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);//#it{v}_{2,J/#psi} wrt #it{p}_{T}
           Figv0malice->GetXaxis()->SetTitle("Centrality - V0M (Percent)");
       Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,tkl}");
        Figv0malice->GetYaxis()->SetRangeUser(-0.15, 0.20);
        Figv0malice->GetXaxis()->SetRangeUser(0, 30.);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetPadBottomMargin(0.15);
        gStyle->SetPadLeftMargin(0.15);
          Figv0malice->AddYError(7, vaeylsys, vaeyhsys);
          Figv0malice->SetMarkerStyle(8);
       Figv0malice->SetMarkerColor(kAzure-3);
          Figv0malice->SetLineColorAlpha(46, 0.);
          Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
          Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
    Figv0malice->GetAttLine(0)->SetLineWidth(2);
    Figv0malice->GetAttLine(1)->SetLineWidth(2);
       Figv0malice->GetAttFill(1)->SetFillStyle(0);
    Figv0malice->GetXaxis()->SetRangeUser(0,30);
       
       TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "", 7, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
           Figv0matlas->GetXaxis()->CenterTitle("Centrality - V0M (Percent)");
       Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,tkl}");
        Figv0matlas->GetYaxis()->SetRangeUser(-0.15, 0.20);
        Figv0matlas->GetXaxis()->SetRangeUser(0, 30.);
        gStyle->SetPadTickX(0);
        gStyle->SetPadTickY(0);
          Figv0matlas->AddYError(7, vbeylsys, vbeyhsys);
             Figv0matlas->SetMarkerStyle(21);
        Figv0matlas->SetMarkerColor(kRed+1);
             Figv0matlas->SetLineColorAlpha(46, 0.);
             Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
             Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
    Figv0matlas->GetAttLine(0)->SetLineWidth(2);
    Figv0matlas->GetAttLine(1)->SetLineWidth(2);
          Figv0matlas->GetAttFill(1)->SetFillStyle(0);
       
          Figempty->Draw("a p s ; z ; 5 s=0.5");
          Figv0malice->Draw("p s ; z ; 5 s=0.5");
       Figv0matlas->Draw("p s ; z ; 5 s=0.5");
       
       TLine *ls2=new TLine(0.,0.0,30.,0.0);
          ls2->SetLineColor(kBlack);
          ls2->SetLineWidth(1);
          ls2->SetLineStyle(9);
          ls2->Draw("same");
       
       TLegend *legendov=new TLegend(0.12,0.65,0.40,0.80);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.035);
           Char_t messagi[80];
          // sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
         //  legendov->AddEntry(Figv0malice,messagi,"");
                   legendov->AddEntry(Figv0malice,"Subtracted yields method");

                    legendov->AddEntry(Figv0matlas,"Template fit method");

                     legendov->Draw();
    
//           TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
//                         legendov2->SetFillColorAlpha(kWhite, 0.);
//                         legendov2->SetBorderSize(0);
//                          legendov2->SetTextFont(42);
//                          legendov2->SetTextSize(0.05);
//               Char_t messago[80];
//               sprintf(messago,"ALICE Preliminary");
//                legendov2->AddEntry(Figv0matlas,messago,"");
//                         legendov2->Draw();
//
//           TLegend *legendov3=new TLegend(0.5,0.67,0.80,0.87);
//                     legendov3->SetFillColorAlpha(kWhite, 0.);
//                     legendov3->SetBorderSize(0);
//                      legendov3->SetTextFont(42);
//                      legendov3->SetTextSize(0.04);
//           sprintf(messago,"V0M (0-5%%)-(40-100%%)");
//            legendov3->AddEntry(Figv0matlas,messago,"");
//           sprintf(messago,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
//           legendov3->AddEntry(Figv0matlas,messago,"");
//            sprintf(messago,"2.5 < #it{y}_{cms} < 4.0");
//            legendov3->AddEntry(Figv0matlas,messago,"");
//           sprintf(messago,"1.5 < |#it{#Delta#eta}| < 5.0");
//           legendov3->AddEntry(Figv0matlas,messago,"");
//                     legendov3->Draw();
//
//           TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
//                     legendov4->SetFillColorAlpha(kWhite, 0.);
//                     legendov4->SetBorderSize(0);
//                      legendov4->SetTextFont(62);
//                      legendov4->SetTextSize(0.04);
//        legendov4->SetTextColor(kAzure-3);
//           sprintf(messago,"5.9%% global syst. uncertainty");
//            legendov4->AddEntry(Figv0matlas,messago,"");
//                     legendov4->Draw();
//
//        TLegend *legendov42=new TLegend(0.08,0.12,0.40,0.17);
//                  legendov42->SetFillColorAlpha(kWhite, 0.);
//                  legendov42->SetBorderSize(0);
//                   legendov42->SetTextFont(62);
//                   legendov42->SetTextSize(0.04);
//        legendov42->SetTextColor(kRed+1);
//        sprintf(messago,"2.7%% global syst. uncertainty");
//         legendov42->AddEntry(Figv0matlas,messago,"");
//                  legendov42->Draw();
//
//
    }
    
    // v2 TKL wrt Centrality SPDT
    {
        
        double yax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
        double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
        double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
        double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
        double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
        double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
        double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
        double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

              auto cv0m2 = new TCanvas("cv0m2","cv0m2",200,10,800,600);
              double vax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
              double vay[]      = {0.06241, 0.06172, 0.06069, 0.06055, 0.05994, 0.05813, 0.05676};
              double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
              double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
              double* vaeylstat = new double[7]  {0.00045, 0.00054, 0.00069, 0.00066, 0.00087, 0.00114, 0.00124};
              double* vaeyhstat = new double[7]  {0.00045, 0.00054, 0.00069, 0.00066, 0.00087, 0.00114, 0.00124};
             // double* vaeylsys  = new double[6]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258,0.0234};
//            double* vaeylsys  = new double[7] {0.0313, 0.0336, 0.0341, 0.0367, 0.0553, 0.0788, 0.0897}; PYTHIA full
//              double* vaeyhsys  = new double[7] {0.0313, 0.0336, 0.0341, 0.0367, 0.0553, 0.0788, 0.0897};
        
                    double* vaeylsys  = new double[7] {0.0137, 0.0144, 0.0145, 0.0152, 0.0180, 0.0203, 0.0244}; //PYTHIA half top bottom
                      double* vaeyhsys  = new double[7] {0.0116, 0.0120, 0.0121, 0.0125, 0.0142, 0.0154, 0.0172};
           
           double vbx[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
           double vby[]      = {0.09687, 0.10436, 0.10816, 0.11416, 0.12519, 0.13274, 0.14855};
           double* vbeylstat = new double[7]  {0.00157, 0.00205, 0.00251, 0.00283, 0.00401, 0.00546, 0.00780};
           double* vbeyhstat = new double[7]  {0.00157, 0.00205, 0.00251, 0.00283, 0.00401, 0.00546, 0.00780};
           //double* vbeylsys  = new double[6]  {0.0093, 0.0074, 0.0136, 0.0104, 0.0245,0.0277};
         //   double* vbeylsys  = new double[7]   {0.0246, 0.0254, 0.0261, 0.0290, 0.0363, 0.0431, 0.0558}; PYTHIA full
           //double* vbeyhsys  = new double[7]   {0.0246, 0.0254, 0.0261, 0.0290, 0.0363, 0.0431, 0.0558};
        
        double* vbeylsys  = new double[7]   {0.0137, 0.0150, 0.0160, 0.0184, 0.0228, 0.0269, 0.0327}; //PYTHIA half top bottom
        double* vbeyhsys  = new double[7]   {0.0129, 0.0142, 0.0153, 0.0177, 0.0218, 0.0256, 0.0307};
        
        TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,tkl} wrt Ceentrality", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
               Figempty->GetXaxis()->SetTitle("Centrality - SPDTracklets (Percent)");
           Figempty->GetYaxis()->SetTitle("#it{v}_{2,tkl}");
            Figempty->GetYaxis()->SetRangeUser(-0.15, 0.20);
            Figempty->GetXaxis()->SetRangeUser(0, 30.);
            gStyle->SetPadTickX(1);
            gStyle->SetPadTickY(1);
              Figempty->AddYError(7, yaeylsys, yaeyhsys);
        Figempty->GetXaxis()->SetRangeUser(0,30);
            
              TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "", 7, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);//#it{v}_{2,J/#psi} wrt #it{p}_{T}
               Figv0malice->GetXaxis()->SetTitle("Centrality - SPDTracklets (Percent)");
           Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,tkl}");
            Figv0malice->GetYaxis()->SetRangeUser(-0.15, 0.20);
            Figv0malice->GetXaxis()->SetRangeUser(0, 30.);
            gStyle->SetPadTickX(1);
            gStyle->SetPadTickY(1);
            gStyle->SetPadBottomMargin(0.15);
            gStyle->SetPadLeftMargin(0.15);
              Figv0malice->AddYError(7, vaeylsys, vaeyhsys);
              Figv0malice->SetMarkerStyle(8);
           Figv0malice->SetMarkerColor(kAzure-3);
              Figv0malice->SetLineColorAlpha(46, 0.);
              Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
              Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
        Figv0malice->GetAttLine(0)->SetLineWidth(2);
        Figv0malice->GetAttLine(1)->SetLineWidth(2);
           Figv0malice->GetAttFill(1)->SetFillStyle(0);
        Figv0malice->GetXaxis()->SetRangeUser(0,30);
           
           TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "", 7, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
               Figv0matlas->GetXaxis()->CenterTitle("Centrality - V0M (Percent)");
           Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,tkl}");
            Figv0matlas->GetYaxis()->SetRangeUser(-0.15, 0.20);
            Figv0matlas->GetXaxis()->SetRangeUser(0, 30.);
            gStyle->SetPadTickX(0);
            gStyle->SetPadTickY(0);
              Figv0matlas->AddYError(7, vbeylsys, vbeyhsys);
                 Figv0matlas->SetMarkerStyle(21);
            Figv0matlas->SetMarkerColor(kRed+1);
                 Figv0matlas->SetLineColorAlpha(46, 0.);
                 Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
                 Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
        Figv0matlas->GetAttLine(0)->SetLineWidth(2);
        Figv0matlas->GetAttLine(1)->SetLineWidth(2);
              Figv0matlas->GetAttFill(1)->SetFillStyle(0);
           
              Figempty->Draw("a p s ; z ; 5 s=0.5");
              Figv0malice->Draw("p s ; z ; 5 s=0.5");
           Figv0matlas->Draw("p s ; z ; 5 s=0.5");
           
           TLine *ls2=new TLine(0.,0.0,30.,0.0);
              ls2->SetLineColor(kBlack);
              ls2->SetLineWidth(1);
              ls2->SetLineStyle(9);
              ls2->Draw("same");
           
           TLegend *legendov=new TLegend(0.12,0.65,0.40,0.80);
                         legendov->SetFillColorAlpha(kWhite, 0.);
                         legendov->SetBorderSize(0);
                          legendov->SetTextFont(42);
                          legendov->SetTextSize(0.035);
               Char_t messagi[80];
              // sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
             //  legendov->AddEntry(Figv0malice,messagi,"");
                       legendov->AddEntry(Figv0malice,"Subtracted yields method");

                        legendov->AddEntry(Figv0matlas,"Template fit method");

                         legendov->Draw();
        
    //           TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
    //                         legendov2->SetFillColorAlpha(kWhite, 0.);
    //                         legendov2->SetBorderSize(0);
    //                          legendov2->SetTextFont(42);
    //                          legendov2->SetTextSize(0.05);
    //               Char_t messago[80];
    //               sprintf(messago,"ALICE Preliminary");
    //                legendov2->AddEntry(Figv0matlas,messago,"");
    //                         legendov2->Draw();
    //
    //           TLegend *legendov3=new TLegend(0.5,0.67,0.80,0.87);
    //                     legendov3->SetFillColorAlpha(kWhite, 0.);
    //                     legendov3->SetBorderSize(0);
    //                      legendov3->SetTextFont(42);
    //                      legendov3->SetTextSize(0.04);
    //           sprintf(messago,"V0M (0-5%%)-(40-100%%)");
    //            legendov3->AddEntry(Figv0matlas,messago,"");
    //           sprintf(messago,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
    //           legendov3->AddEntry(Figv0matlas,messago,"");
    //            sprintf(messago,"2.5 < #it{y}_{cms} < 4.0");
    //            legendov3->AddEntry(Figv0matlas,messago,"");
    //           sprintf(messago,"1.5 < |#it{#Delta#eta}| < 5.0");
    //           legendov3->AddEntry(Figv0matlas,messago,"");
    //                     legendov3->Draw();
    //
    //           TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
    //                     legendov4->SetFillColorAlpha(kWhite, 0.);
    //                     legendov4->SetBorderSize(0);
    //                      legendov4->SetTextFont(62);
    //                      legendov4->SetTextSize(0.04);
    //        legendov4->SetTextColor(kAzure-3);
    //           sprintf(messago,"5.9%% global syst. uncertainty");
    //            legendov4->AddEntry(Figv0matlas,messago,"");
    //                     legendov4->Draw();
    //
    //        TLegend *legendov42=new TLegend(0.08,0.12,0.40,0.17);
    //                  legendov42->SetFillColorAlpha(kWhite, 0.);
    //                  legendov42->SetBorderSize(0);
    //                   legendov42->SetTextFont(62);
    //                   legendov42->SetTextSize(0.04);
    //        legendov42->SetTextColor(kRed+1);
    //        sprintf(messago,"2.7%% global syst. uncertainty");
    //         legendov42->AddEntry(Figv0matlas,messago,"");
    //                  legendov42->Draw();
    //
    //
        }

    
    // Final v2 V0M PYTHIA syst
       {
       
       double yax[]      = {1, 2.5, 3.5, 5, 7, 10, 12};
       double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
       double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
       double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
       double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
       double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
       double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
       double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

             auto cv0mprop = new TCanvas("cv0mprop","cv0mprop",200,10,800,600);
             double vax[]      = {1, 2.5, 3.5, 5, 7, 10};
             double vay[]      = {-0.00584,0.01802,-0.01266,0.03064,-0.05892,0.03671};
             double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
             double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
             double* vaeylstat = new double[6]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503,0.02098};
             double* vaeyhstat = new double[6]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503,0.02098};
            // double* vaeylsys  = new double[6]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258,0.0234};
         //  double* vaeylsys  = new double[6] {0.0092, 0.0088, 0.0141, 0.0123, 0.0248, 0.0253}; //PYTHIA Full
           //  double* vaeyhsys  = new double[6] {0.0092, 0.0088, 0.0141, 0.0123, 0.0248, 0.0253};
           
           double* vaeylsys  = new double[6] {0.0090, 0.0072, 0.0136, 0.0088, 0.0183, 0.0230}; //PYTHIA Full non prop tkl
           double* vaeyhsys  = new double[6] {0.0090, 0.0072, 0.0136, 0.0088, 0.0183, 0.0230};
           
//           double* vaeylsys  = new double[6] {0.0090, 0.0076, 0.0137, 0.0096, 0.0201, 0.0235}; //PYTHIA Half top bottom
//           double* vaeyhsys  = new double[6] {0.0090, 0.0076, 0.0137, 0.0098, 0.0198, 0.0236};
          
          double vbx[]      = {1.2, 2.7, 3.7, 5.2, 7.2, 10.2};
          double vby[]      = {-0.00441, 0.01912, -0.01151, 0.03532,-0.04245,0.04257};
          double* vbeylstat = new double[6]  {0.00820, 0.00972, 0.01101, 0.00962,0.01472,0.02161};
          double* vbeyhstat = new double[6]  {0.00820, 0.00972, 0.01101, 0.00962,0.01472,0.02161};
          //double* vbeylsys  = new double[6]  {0.0093, 0.0074, 0.0136, 0.0104, 0.0245,0.0277};
        //   double* vbeylsys  = new double[6]   {0.0087, 0.0075, 0.0133, 0.0124, 0.0186,0.0271};// PYTHIA Full
         // double* vbeyhsys  = new double[6]   {0.0087, 0.0075, 0.0133, 0.0124, 0.0186,0.0271};
           
           double* vbeylsys  = new double[6]   {0.0087, 0.0066, 0.0131, 0.0104, 0.0167,0.0258};// PYTHIA Full no tkl prop
           double* vbeyhsys  = new double[6]   {0.0087, 0.0066, 0.0131, 0.0104, 0.0167,0.0258};
           
//           double* vbeylsys  = new double[6]   {0.0087, 0.0068, 0.0132, 0.0109, 0.0172,0.0261}; //PYTHIA Half top bottom
//           double* vbeyhsys  = new double[6]   {0.0087, 0.0068, 0.0132, 0.0109, 0.0171,0.0262};
       
       TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
              Figempty->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
          Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
           Figempty->GetYaxis()->SetRangeUser(-0.15, 0.20);
           Figempty->GetXaxis()->SetRangeUser(0, 12.);
           gStyle->SetPadTickX(1);
           gStyle->SetPadTickY(1);
             Figempty->AddYError(7, yaeylsys, yaeyhsys);
       Figempty->GetXaxis()->SetRangeUser(0,15);
           
             TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "", 6, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);//#it{v}_{2,J/#psi} wrt #it{p}_{T}
              Figv0malice->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
          Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
           Figv0malice->GetYaxis()->SetRangeUser(-0.15, 0.20);
           Figv0malice->GetXaxis()->SetRangeUser(0, 12.);
           gStyle->SetPadTickX(1);
           gStyle->SetPadTickY(1);
           gStyle->SetPadBottomMargin(0.15);
           gStyle->SetPadLeftMargin(0.15);
             Figv0malice->AddYError(6, vaeylsys, vaeyhsys);
             Figv0malice->SetMarkerStyle(8);
          Figv0malice->SetMarkerColor(kAzure-3);
             Figv0malice->SetLineColorAlpha(46, 0.);
             Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
             Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
       Figv0malice->GetAttLine(0)->SetLineWidth(2);
       Figv0malice->GetAttLine(1)->SetLineWidth(2);
          Figv0malice->GetAttFill(1)->SetFillStyle(0);
       Figv0malice->GetXaxis()->SetRangeUser(0,15);
          
          TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "", 6, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
              Figv0matlas->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
          Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
           Figv0matlas->GetYaxis()->SetRangeUser(-0.15, 0.20);
           Figv0matlas->GetXaxis()->SetRangeUser(0, 12.);
           gStyle->SetPadTickX(0);
           gStyle->SetPadTickY(0);
             Figv0matlas->AddYError(6, vbeylsys, vbeyhsys);
                Figv0matlas->SetMarkerStyle(21);
           Figv0matlas->SetMarkerColor(kRed+1);
                Figv0matlas->SetLineColorAlpha(46, 0.);
                Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
                Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
       Figv0matlas->GetAttLine(0)->SetLineWidth(2);
       Figv0matlas->GetAttLine(1)->SetLineWidth(2);
             Figv0matlas->GetAttFill(1)->SetFillStyle(0);
          
             Figempty->Draw("a p s ; z ; 5 s=0.5");
             Figv0malice->Draw("p s ; z ; 5 s=0.5");
          Figv0matlas->Draw("p s ; z ; 5 s=0.5");
          
          TLine *ls2=new TLine(0.,0.0,12.,0.0);
             ls2->SetLineColor(kBlack);
             ls2->SetLineWidth(1);
             ls2->SetLineStyle(9);
             ls2->Draw("same");
          
          TLegend *legendov=new TLegend(0.12,0.65,0.40,0.80);
                        legendov->SetFillColorAlpha(kWhite, 0.);
                        legendov->SetBorderSize(0);
                         legendov->SetTextFont(42);
                         legendov->SetTextSize(0.035);
              Char_t messagi[80];
             // sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
            //  legendov->AddEntry(Figv0malice,messagi,"");
                      legendov->AddEntry(Figv0malice,"Subtracted yields method");

                       legendov->AddEntry(Figv0matlas,"Template fit method");

                        legendov->Draw();
       
              TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
                            legendov2->SetFillColorAlpha(kWhite, 0.);
                            legendov2->SetBorderSize(0);
                             legendov2->SetTextFont(42);
                             legendov2->SetTextSize(0.05);
                  Char_t messago[80];
                  sprintf(messago,"ALICE Preliminary");
                   legendov2->AddEntry(Figv0matlas,messago,"");
                            legendov2->Draw();
              
              TLegend *legendov3=new TLegend(0.5,0.67,0.80,0.87);
                        legendov3->SetFillColorAlpha(kWhite, 0.);
                        legendov3->SetBorderSize(0);
                         legendov3->SetTextFont(42);
                         legendov3->SetTextSize(0.04);
              sprintf(messago,"V0M (0-5%%)-(40-100%%)");
               legendov3->AddEntry(Figv0matlas,messago,"");
              sprintf(messago,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
              legendov3->AddEntry(Figv0matlas,messago,"");
               sprintf(messago,"2.5 < #it{y}_{cms} < 4.0");
               legendov3->AddEntry(Figv0matlas,messago,"");
              sprintf(messago,"1.5 < |#it{#Delta#eta}| < 5.0");
              legendov3->AddEntry(Figv0matlas,messago,"");
                        legendov3->Draw();
              
              TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
                        legendov4->SetFillColorAlpha(kWhite, 0.);
                        legendov4->SetBorderSize(0);
                         legendov4->SetTextFont(62);
                         legendov4->SetTextSize(0.04);
           legendov4->SetTextColor(kAzure-3);
              sprintf(messago,"28%% global syst. uncertainty");
               legendov4->AddEntry(Figv0matlas,messago,"");
                        legendov4->Draw();
           
           TLegend *legendov42=new TLegend(0.08,0.12,0.40,0.17);
                     legendov42->SetFillColorAlpha(kWhite, 0.);
                     legendov42->SetBorderSize(0);
                      legendov42->SetTextFont(62);
                      legendov42->SetTextSize(0.04);
           legendov42->SetTextColor(kRed+1);
           sprintf(messago,"19%% global syst. uncertainty");
            legendov42->AddEntry(Figv0matlas,messago,"");
                     legendov42->Draw();
       
       
       }
       
       // Final v2 SPDT PYTHIA syst
       {
       
       double yax[]      = {1, 2.5, 3.5, 5, 7, 10, 12};
       double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
       double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
       double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
       double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
       double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
       double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
       double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

             auto cspdtprop = new TCanvas("cspdtprop","cspdtprop",200,10,600,400);
             double vax[]      = {1, 2.5, 3.5, 5, 7, 10};
             double vay[]      = {-0.01032,0.01880,-0.00154,0.04516,0.02240,0.13548};
             double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
             double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
             double* vaeylstat = new double[6]  {0.00644, 0.00744, 0.00825, 0.00728, 0.01110,0.01487};
             double* vaeyhstat = new double[6]  {0.00644, 0.00744, 0.00825, 0.00728, 0.01110,0.01487};
             //double* vaeylsys  = new double[6]  {0.0087,0.0042,0.0120,0.0058,0.0191,0.0324};
      //     double* vaeylsys  = new double[6]  {0.0098,0.0105,0.0106,0.0247,0.0195,0.0742}; //PYTHIA Full
        //     double* vaeyhsys  = new double[6] {0.0098,0.0105,0.0106,0.0247,0.0195,0.0742};
           
           double* vaeylsys  = new double[6]  {0.0084,0.0049,0.0106,0.0105,0.0161,0.0314}; //PYTHIA half top bottom
           double* vaeyhsys  = new double[6] {0.0083,0.0055,0.0106,0.0119,0.0163,0.0357};
          
          double vbx[]      = {1.2, 2.7, 3.7, 5.2, 7.2, 10.2};
          double vby[]      = {-0.00876, 0.00805, -0.01061, 0.01805,-0.01715,0.05059};
          double* vbeylstat = new double[6]  {0.00488, 0.00596, 0.00698, 0.00647,0.01044,0.01645};
          double* vbeyhstat = new double[6]  {0.00488, 0.00596, 0.00698, 0.00647,0.01044,0.01645};
         // double* vbeylsys  = new double[6]  {0.0069, 0.0038, 0.0091, 0.0057,0.0145,0.0412};
      //     double* vbeylsys  = new double[6]  {0.0068, 0.0039, 0.0087, 0.0074,0.0144,0.0261}; //PYTHIA Full
        //  double* vbeyhsys  = new double[6]  {0.0068, 0.0039, 0.0087, 0.0074,0.0144,0.0261};
           
           double* vbeylsys  = new double[6]  {0.0066, 0.0035, 0.0084, 0.0063,0.0139,0.0238}; //PYTHIA Half top bottom
           double* vbeyhsys  = new double[6]  {0.0066, 0.0035, 0.0084, 0.0064,0.0140,0.0239};
       
       TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
              Figempty->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
          Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
             Figempty->AddYError(7, yaeylsys, yaeyhsys);
       Figempty->GetXaxis()->SetRangeUser(0,15);
           
             TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 6, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);
              Figv0malice->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
          Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
             Figv0malice->AddYError(6, vaeylsys, vaeyhsys);
             Figv0malice->SetMarkerStyle(8);
          Figv0malice->SetMarkerColor(kAzure-3);
             Figv0malice->SetLineColorAlpha(46, 0.);
             Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
             Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
       Figv0malice->GetAttLine(0)->SetLineWidth(2);
       Figv0malice->GetAttLine(1)->SetLineWidth(2);
          Figv0malice->GetAttFill(1)->SetFillStyle(0);
       Figv0malice->GetXaxis()->SetRangeUser(0,15);
          
          TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 6, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
              Figv0matlas->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/c)");
          Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
             Figv0matlas->AddYError(6, vbeylsys, vbeyhsys);
                Figv0matlas->SetMarkerStyle(21);
           Figv0matlas->SetMarkerColor(kRed+1);
                Figv0matlas->SetLineColorAlpha(46, 0.);
                Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
                Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
       Figv0matlas->GetAttLine(0)->SetLineWidth(2);
       Figv0matlas->GetAttLine(1)->SetLineWidth(2);
             Figv0matlas->GetAttFill(1)->SetFillStyle(0);
          
             Figempty->Draw("a p s ; ; 5 s=0.5");
             Figv0malice->Draw("p s ; ; 5 s=0.5");
          Figv0matlas->Draw("p s ; ; 5 s=0.5");
          
          TLine *ls2=new TLine(0.,0.0,12.,0.0);
             ls2->SetLineColor(kBlack);
             ls2->SetLineWidth(1);
             ls2->SetLineStyle(9);
             ls2->Draw("same");
          
          TLegend *legendov=new TLegend(0.12,0.70,0.40,0.90);
                        legendov->SetFillColorAlpha(kWhite, 0.);
                        legendov->SetBorderSize(0);
                         legendov->SetTextFont(42);
                         legendov->SetTextSize(0.022);
              Char_t messagi[80];
              sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
              legendov->AddEntry(Figv0malice,messagi,"");
                      legendov->AddEntry(Figv0malice,"Subtracted yields default: Y_{C}-Y_{P}=a_{0}+2a_{1}cos(#Delta#phi)+2a_{2}cos(2#Delta#phi)");

                       legendov->AddEntry(Figv0matlas,"ATLAS Template fit: Y_{C}-F(Y_{P}-Y_{P}(0))=G(1+#it{v}_{2,2}cos(2#Delta#phi))");

                        legendov->Draw();
       
              TLegend *legendov2=new TLegend(0.12,0.60,0.40,0.70);
                            legendov2->SetFillColorAlpha(kWhite, 0.);
                            legendov2->SetBorderSize(0);
                             legendov2->SetTextFont(72);
                             legendov2->SetTextSize(0.03);
                  Char_t messago[80];
                  sprintf(messago,"ALICE Requested for preliminary");
                   legendov2->AddEntry(Figv0matlas,messago,"");
                            legendov2->Draw();
              
              TLegend *legendov3=new TLegend(0.6,0.70,0.90,0.90);
                        legendov3->SetFillColorAlpha(kWhite, 0.);
                        legendov3->SetBorderSize(0);
                         legendov3->SetTextFont(62);
                         legendov3->SetTextSize(0.03);
              sprintf(messago,"SPDTracklets (0-5%)-(40-100%)");
               legendov3->AddEntry(Figv0matlas,messago,"");
              sprintf(messago,"ALICE pp, #sqrt{s_{NN}}=13 TeV");
              legendov3->AddEntry(Figv0matlas,messago,"");
              sprintf(messago,"1.5<|#Delta#eta|<5.0");
              legendov3->AddEntry(Figv0matlas,messago,"");
                        legendov3->Draw();
              
              TLegend *legendov4=new TLegend(0.45,0.10,0.85,0.20);
                        legendov4->SetFillColorAlpha(kWhite, 0.);
                        legendov4->SetBorderSize(0);
                         legendov4->SetTextFont(62);
                         legendov4->SetTextSize(0.03);
                           legendov4->SetTextColor(kAzure-3);
              sprintf(messago,"7.9 percent global syst. uncertainty");
               legendov4->AddEntry(Figv0matlas,messago,"");
                        legendov4->Draw();
           
           TLegend *legendov42=new TLegend(0.45,0.00,0.85,0.10);
                     legendov42->SetFillColorAlpha(kWhite, 0.);
                     legendov42->SetBorderSize(0);
                      legendov42->SetTextFont(62);
                      legendov42->SetTextSize(0.03);
                        legendov42->SetTextColor(kRed+1);
           sprintf(messago,"10.9 percent global syst. uncertainty");
            legendov42->AddEntry(Figv0matlas,messago,"");
                     legendov42->Draw();
       
       
       }
    
    
    
    // Final v2 SPDC PYTHIA syst
     {
     
     double yax[]      = {1, 2.5, 3.5, 5, 7, 10, 12};
     double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
     double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
     double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
     double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
     double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
     double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
     double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

           auto cspdcprop = new TCanvas("cspdcprop","cspdcprop",200,10,600,400);
           double vax[]      = {1, 2.5, 3.5, 5, 7, 10};
           double vay[]      = {-0.0072621,0.0158289,0.011473,0.0347869,0.0348071,0.136301};
           double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
           double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
           double* vaeylstat = new double[6]  {0.00544949, 0.00650557, 0.00741776, 0.00679115, 0.0107516,0.0149132};
           double* vaeyhstat = new double[6]  {0.00544949, 0.00650557, 0.00741776, 0.00679115, 0.0107516,0.0149132};
           //double* vaeylsys  = new double[6]  {0.0087,0.0042,0.0120,0.0058,0.0191,0.0324};

                  double* vaeylsys  = new double[6]  {0.0090,0.0145,0.0222,0.0222,0.0369,0.0854}; //PYTHIA full
                  double* vaeyhsys  = new double[6] {0.0090,0.0145,0.0222,0.0222,0.0369,0.0854};
         
//         double* vaeylsys  = new double[6]  {0.0086,0.0130,0.0217,0.0171,0.0341,0.0649}; //PYTHIA half top bottom
//         double* vaeyhsys  = new double[6] {0.0085,0.0131,0.0217,0.0175,0.0343,0.0664};
        
        double vbx[]      = {1.2, 2.7, 3.7, 5.2, 7.2, 10.2};
        double vby[]      = {-0.006617, 0.00673282, -0.0006353, 0.00748165,-0.0077531,0.0530079};
        double* vbeylstat = new double[6]  {0.00435074, 0.00562463, 0.00697922, 0.00686317,0.011278,0.0184729};
        double* vbeyhstat = new double[6]  {0.00435074, 0.00562463, 0.00697922, 0.00686317,0.011278,0.0184729};
       // double* vbeylsys  = new double[6]  {0.0069, 0.0038, 0.0091, 0.0057,0.0145,0.0412};
         
         
                  double* vbeylsys  = new double[6]  {0.0067, 0.0128, 0.0157, 0.0240,0.0183,0.0578}; //PYTHIA full
                  double* vbeyhsys  = new double[6]  {0.0067, 0.0128, 0.0157, 0.0240,0.0183,0.0578};
         
         
//         double* vbeylsys  = new double[6]  {0.0066, 0.0127, 0.0157, 0.0239,0.0183,0.0571}; //PYTHIA Half top bottom
//         double* vbeyhsys  = new double[6]  {0.0066, 0.0127, 0.0157, 0.0239,0.0183,0.0571};
     
     TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
            Figempty->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
           Figempty->AddYError(7, yaeylsys, yaeyhsys);
     Figempty->GetXaxis()->SetRangeUser(0,15);
         
           TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 6, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);
            Figv0malice->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
           Figv0malice->AddYError(6, vaeylsys, vaeyhsys);
           Figv0malice->SetMarkerStyle(8);
        Figv0malice->SetMarkerColor(kAzure-3);
           Figv0malice->SetLineColorAlpha(46, 0.);
           Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
           Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
     Figv0malice->GetAttLine(0)->SetLineWidth(2);
     Figv0malice->GetAttLine(1)->SetLineWidth(2);
        Figv0malice->GetAttFill(1)->SetFillStyle(0);
     Figv0malice->GetXaxis()->SetRangeUser(0,15);
        
        TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 6, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
            Figv0matlas->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/c)");
        Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
           Figv0matlas->AddYError(6, vbeylsys, vbeyhsys);
              Figv0matlas->SetMarkerStyle(21);
         Figv0matlas->SetMarkerColor(kRed+1);
              Figv0matlas->SetLineColorAlpha(46, 0.);
              Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
              Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
     Figv0matlas->GetAttLine(0)->SetLineWidth(2);
     Figv0matlas->GetAttLine(1)->SetLineWidth(2);
           Figv0matlas->GetAttFill(1)->SetFillStyle(0);
        
           Figempty->Draw("a p s ; ; 5 s=0.5");
           Figv0malice->Draw("p s ; ; 5 s=0.5");
        Figv0matlas->Draw("p s ; ; 5 s=0.5");
        
        TLine *ls2=new TLine(0.,0.0,12.,0.0);
           ls2->SetLineColor(kBlack);
           ls2->SetLineWidth(1);
           ls2->SetLineStyle(9);
           ls2->Draw("same");
        
        TLegend *legendov=new TLegend(0.12,0.70,0.40,0.90);
                      legendov->SetFillColorAlpha(kWhite, 0.);
                      legendov->SetBorderSize(0);
                       legendov->SetTextFont(42);
                       legendov->SetTextSize(0.022);
            Char_t messagi[80];
            sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
            legendov->AddEntry(Figv0malice,messagi,"");
                    legendov->AddEntry(Figv0malice,"Subtracted yields default: Y_{C}-Y_{P}=a_{0}+2a_{1}cos(#Delta#phi)+2a_{2}cos(2#Delta#phi)");

                     legendov->AddEntry(Figv0matlas,"ATLAS Template fit: Y_{C}-F(Y_{P}-Y_{P}(0))=G(1+#it{v}_{2,2}cos(2#Delta#phi))");

                      legendov->Draw();
     
            TLegend *legendov2=new TLegend(0.12,0.60,0.40,0.70);
                          legendov2->SetFillColorAlpha(kWhite, 0.);
                          legendov2->SetBorderSize(0);
                           legendov2->SetTextFont(72);
                           legendov2->SetTextSize(0.03);
                Char_t messago[80];
                sprintf(messago,"ALICE Requested for preliminary");
                 legendov2->AddEntry(Figv0matlas,messago,"");
                          legendov2->Draw();
            
            TLegend *legendov3=new TLegend(0.6,0.70,0.90,0.90);
                      legendov3->SetFillColorAlpha(kWhite, 0.);
                      legendov3->SetBorderSize(0);
                       legendov3->SetTextFont(62);
                       legendov3->SetTextSize(0.03);
            sprintf(messago,"SPDClusters (0-5%)-(40-100%)");
             legendov3->AddEntry(Figv0matlas,messago,"");
            sprintf(messago,"ALICE pp, #sqrt{s_{NN}}=13 TeV");
            legendov3->AddEntry(Figv0matlas,messago,"");
            sprintf(messago,"1.5<|#Delta#eta|<5.0");
            legendov3->AddEntry(Figv0matlas,messago,"");
                      legendov3->Draw();
            
            TLegend *legendov4=new TLegend(0.45,0.10,0.85,0.20);
                      legendov4->SetFillColorAlpha(kWhite, 0.);
                      legendov4->SetBorderSize(0);
                       legendov4->SetTextFont(62);
                       legendov4->SetTextSize(0.03);
                         legendov4->SetTextColor(kAzure-3);
            sprintf(messago,"7.9 percent global syst. uncertainty");
             legendov4->AddEntry(Figv0matlas,messago,"");
                      legendov4->Draw();
         
         TLegend *legendov42=new TLegend(0.45,0.00,0.85,0.10);
                   legendov42->SetFillColorAlpha(kWhite, 0.);
                   legendov42->SetBorderSize(0);
                    legendov42->SetTextFont(62);
                    legendov42->SetTextSize(0.03);
                      legendov42->SetTextColor(kRed+1);
         sprintf(messago,"10.9 percent global syst. uncertainty");
          legendov42->AddEntry(Figv0matlas,messago,"");
                   legendov42->Draw();
     
     
     }



    
    {
            
            // Figure Combine systems

               auto csysprop = new TCanvas("csysprop","csysprop",200,10,800,600);
               double qax[]      = {1, 2.5, 3.5, 5, 7};
               double qay[]      = {0,-0.015,0.031,0.082,0.033};
               double qaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
               double qaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
               double* qaeylstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
               double* qaeyhstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
               double* qaeylsys  = new double[5]  {0.008,0.011,0.013,0.011,0.017};
               double* qaeyhsys  = new double[5]   {0.008,0.011,0.013,0.011,0.017};
            
            for(int idx=0; idx<5; idx++){
                qay[idx]/=0.068;
                qaeylstat[idx]/=0.068;
                qaeyhstat[idx]/=0.068;
                qaeylsys[idx]/=0.068;
                qaeyhsys[idx]/=0.068;
            }
            
            double qbx[]      = {1, 2.5, 3.5, 5, 7};
            double qby[]      = {0.023,0.016,0.074,0.069,0.079};
            double* qbeylstat = new double[5]  {0.015, 0.018, 0.022, 0.020, 0.034};
            double* qbeyhstat = new double[5]  {0.015, 0.018, 0.022, 0.020, 0.034};
            double* qbeylsys  = new double[5]  {0.014, 0.014, 0.015, 0.015, 0.022};
            double* qbeyhsys  = new double[5]   {0.014, 0.014, 0.015, 0.015, 0.022};
            
            for(int idx=0; idx<5; idx++){
                qby[idx]/=0.068;
                qbeylstat[idx]/=0.068;
                qbeyhstat[idx]/=0.068;
                qbeylsys[idx]/=0.068;
                qbeyhsys[idx]/=0.068;
            }
            
            double qcx[]      = {1, 2.5, 3.5, 5, 7};
            double oldqcy[]      = {-0.00584, 0.01802, -0.01266, 0.03064, -0.05892}; //10mrad
            double qcy[]  {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684};
            double* qceylstat = new double[5]  {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613};
                    double* qceyhstat = new double[5]  {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613};
          //  double* qceylsys  = new double[5]  {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
        //    double* qceyhsys  = new double[5]   {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
        
  //  double* qceylsys  = new double[5]  {0.0092, 0.0088, 0.0141, 0.0123, 0.0248};
   //         double* qceyhsys  = new double[5]  {0.0092, 0.0088, 0.0141, 0.0123, 0.0248};
        
        double* qceylsys  = new double[5] {0.0090, 0.0076, 0.0137, 0.0096, 0.0201}; //PYTHIA Half top bottom
        double* qceyhsys  = new double[5] {0.0090, 0.0076, 0.0137, 0.0098, 0.0198};
        
        
        
        
    //       double* qceylstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503};
    //        double* qceyhstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503}; 10mrad
    //        double* qceylsys  = new double[5]  {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
    //        double* qceyhsys  = new double[5]   {0.0099, 0.0137, 0.0168, 0.0189, 0.0380}; 10 mrad
    //        {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684}; 5mr
    //        {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613}; 5mr

            for(int idx=0; idx<5; idx++){
                qcy[idx]/=0.066;
                qceylstat[idx]/=0.066;
                qceyhstat[idx]/=0.066;
                qceylsys[idx]/=0.066;
                qceylsys[idx]/=abs(qcy[idx]*0.066/oldqcy[idx]);
                qceyhsys[idx]/=0.066;
                qceyhsys[idx]/=abs(qcy[idx]*0.066/oldqcy[idx]);
            }
            
            //0.0550667 0.000573161
             
               TGraphMultiErrors* Figsys = new TGraphMultiErrors("Figsys", "#it{#epsilon}_{2,J/#psi} wrt #it{p}_{T} - Systems", 5, qax, qay, qaexl, qaexh, qaeylstat, qaeyhstat);
                Figsys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            Figsys->GetYaxis()->SetTitle("#it{#epsilon}_{2,J/#psi}");
               Figsys->AddYError(5, qaeylsys, qaeyhsys);
               Figsys->SetMarkerStyle(24);
            Figsys->SetMarkerColor(kGray+2);
               Figsys->SetLineColorAlpha(46, 0.);
               Figsys->GetAttLine(0)->SetLineColor(kGray+2);
               Figsys->GetAttLine(1)->SetLineColor(kGray+2);
                Figsys->GetAttLine(0)->SetLineWidth(2);
                Figsys->GetAttLine(1)->SetLineWidth(2);
            Figsys->GetAttFill(1)->SetFillStyle(0);
            
            TGraphMultiErrors* Fig5Pbp = new TGraphMultiErrors("Fig5PbPb", "#it{#epsilon}_{2,J/#psi} wrt #it{p}_{T}", 5, qbx, qby, qaexl, qaexh, qbeylstat, qbeyhstat);
                Fig5Pbp->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
            Fig5Pbp->GetYaxis()->CenterTitle("#it{#epsilon}_{2,J/#psi}");
               Fig5Pbp->AddYError(5, qbeylsys, qbeyhsys);
                  Fig5Pbp->SetMarkerStyle(24);
               Fig5Pbp->SetMarkerColor(kBlack);
                  Fig5Pbp->SetLineColorAlpha(46, 0.);
                  Fig5Pbp->GetAttLine(0)->SetLineColor(kBlack);
                  Fig5Pbp->GetAttLine(1)->SetLineColor(kBlack);
                Fig5Pbp->GetAttLine(0)->SetLineWidth(2);
                Fig5Pbp->GetAttLine(1)->SetLineWidth(2);
               Fig5Pbp->GetAttFill(1)->SetFillStyle(0);
                
            TGraphMultiErrors* Fig5pp = new TGraphMultiErrors("Fig5pp", "#it{#epsilon}_{2,J/#psi} wrt #it{p}_{T}", 5, qcx, qcy, qaexl, qaexh, qceylstat, qceyhstat);
                Fig5pp->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
            Fig5pp->GetYaxis()->CenterTitle("#it{#epsilon}_{2,J/#psi}");
            Fig5pp->AddYError(5, qceylsys, qceyhsys);
               Fig5pp->SetMarkerStyle(20);
            Fig5pp->SetMarkerColor(kAzure-3);
               Fig5pp->SetLineColorAlpha(46, 0.);
               Fig5pp->GetAttLine(0)->SetLineColor(kAzure-3);
            Fig5pp->GetAttLine(1)->SetLineColor(kAzure-3);
                Fig5pp->GetAttLine(0)->SetLineWidth(2);
                Fig5pp->GetAttLine(1)->SetLineWidth(2);
             Fig5pp->GetAttFill(1)->SetFillStyle(0);
               
               Figsys->Draw("a p s ; z ; 5 s=0.5");
            Fig5Pbp->Draw("p s ; z ; 5 s=0.5");
            Fig5pp->Draw("p s ; z ; 5 s=0.5");
            
            TLine *ls=new TLine(0.,0.0,8.,0.0);
               ls->SetLineColor(kBlack);
               ls->SetLineWidth(1);
               ls->SetLineStyle(9);
               ls->Draw("same");
                
                TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
                          legendov2->SetFillColorAlpha(kWhite, 0.);
                          legendov2->SetBorderSize(0);
                           legendov2->SetTextFont(42);
                           legendov2->SetTextSize(0.05);
                Char_t messago[80];
                sprintf(messago,"ALICE Preliminary");
                 legendov2->AddEntry(Fig5pp,messago,"");
                          legendov2->Draw();
                
    //            TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
    //                         legendov4->SetFillColorAlpha(kWhite, 0.);
    //                         legendov4->SetBorderSize(0);
    //                          legendov4->SetTextFont(62);
    //                          legendov4->SetTextSize(0.04);
    //            legendov4->SetTextColor(kAzure-3);
    //               sprintf(messago,"5.9%% global syst. uncertainty");
    //                legendov4->AddEntry(Fig5pp,messago,"");
    //                         legendov4->Draw();
            
            TLegend *legendsys=new TLegend(0.12,0.80,0.60,0.90);
              legendsys->SetFillColorAlpha(kWhite, 0.);
              legendsys->SetBorderSize(0);
                legendsys->SetTextFont(42);
              legendsys->SetTextSize(0.03);
             sprintf(message,"Pb-p, #sqrt{#it{s}_{NN}} = 5.02, 8.16 TeV, (0-20%%)-(40-100%%) (PLB 780 (2018) 7-20)");
                       legendsys->AddEntry(Fig5Pbp,message);
                       sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, -4.46 < #it{y}_{cms} < 2.96");
                   legendsys->AddEntry(Fig5Pbp,message,"");
                sprintf(message,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02, 8.16 TeV, (0-20%%)-(40-100%%) (PLB 780 (2018) 7-20)");
                legendsys->AddEntry(Figsys,message);
                sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.03 < #it{y}_{cms} < 3.53");
            legendsys->AddEntry(Figsys,message,"");
            sprintf(message,"pp, #sqrt{#it{s}_{NN}} = 13 TeV, (0-5%%)-(40-100%%)");
            legendsys->AddEntry(Fig5pp,message);
                sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.5 < #it{y}_{cms} < 4.0");
                legendsys->AddEntry(Fig5pp,message,"");
              legendsys->Draw();
            
            
        }
    
    
    {
        
        // Figure Combine systems

           auto csys = new TCanvas("csys","csys",200,10,800,600);
           double qax[]      = {1, 2.5, 3.5, 5, 7};
           double qay[]      = {0,-0.015,0.031,0.082,0.033};
           double qaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
           double qaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5};
           double* qaeylstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
           double* qaeyhstat = new double[5]  {0.018, 0.021, 0.024, 0.021, 0.033};
           double* qaeylsys  = new double[5]  {0.008,0.011,0.013,0.011,0.017};
           double* qaeyhsys  = new double[5]   {0.008,0.011,0.013,0.011,0.017};
        
        double qbx[]      = {0.64, 1.49, 2.47, 3.46, 4.45, 5.45, 6.819};
        double qby[]      = {0.011, 0.043, 0.074, 0.088,0.085,0.103, 0.083};
            double qbexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
            double qbexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        double* qbeylstat = new double[7]  {0.0085, 0.0069, 0.0069, 0.0077, 0.009, 0.011, 0.011};
        double* qbeyhstat = new double[7]  {0.0085, 0.0069, 0.0069, 0.0077, 0.009, 0.011, 0.011};
        double* qbeylsys  = new double[7]  {0.0038391, 0.0036633, 0.004898, 0.0035068, 0.0037855,0.0029726,  0.0036802};
        double* qbeyhsys  = new double[7]   {0.0038391, 0.0036633, 0.004898, 0.0035068, 0.0037855,0.0029726,  0.0036802};
        
       double qcx[]      = {1, 2.5, 3.5, 5, 7};
                double oldqcy[]      = {-0.00584, 0.01802, -0.01266, 0.03064, -0.05892}; //10mrad
                double qcy[]  {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684};
                double* qceylstat = new double[5]  {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613};
                        double* qceyhstat = new double[5]  {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613};
              //  double* qceylsys  = new double[5]  {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
            //    double* qceyhsys  = new double[5]   {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
        
        //  double* qceylsys  = new double[5]  {0.0092, 0.0088, 0.0141, 0.0123, 0.0248};
        //         double* qceyhsys  = new double[5]  {0.0092, 0.0088, 0.0141, 0.0123, 0.0248};
             
             double* qceylsys  = new double[5] {0.0090, 0.0076, 0.0137, 0.0096, 0.0201}; //PYTHIA Half top bottom
             double* qceyhsys  = new double[5] {0.0090, 0.0076, 0.0137, 0.0098, 0.0198};
        
        //       double* qceylstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503};
        //        double* qceyhstat = new double[5]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503}; 10mrad
        //        double* qceylsys  = new double[5]  {0.0099, 0.0137, 0.0168, 0.0189, 0.0380};
        //        double* qceyhsys  = new double[5]   {0.0099, 0.0137, 0.0168, 0.0189, 0.0380}; 10 mrad
        //        {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684}; 5mr
        //        {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613}; 5mr

                for(int idx=0; idx<5; idx++){
                    qceylsys[idx]/=abs(qcy[idx]/oldqcy[idx]);
                    qceyhsys[idx]/=abs(qcy[idx]/oldqcy[idx]);
                }
        
        //0.0550667 0.000573161
         
           TGraphMultiErrors* Figsys = new TGraphMultiErrors("Figsys", "#it{v}_{2,J/#psi} wrt #it{p}_{T} - Systems (Updated)", 5, qax, qay, qaexl, qaexh, qaeylstat, qaeyhstat);
            Figsys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        Figsys->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
           Figsys->AddYError(5, qaeylsys, qaeyhsys);
           Figsys->SetMarkerStyle(24);
        Figsys->SetMarkerColor(kGray+2);
           Figsys->SetLineColorAlpha(46, 0.);
           Figsys->GetAttLine(0)->SetLineColor(kGray+2);
           Figsys->GetAttLine(1)->SetLineColor(kGray+2);
            Figsys->GetAttLine(0)->SetLineWidth(2);
            Figsys->GetAttLine(1)->SetLineWidth(2);
        Figsys->GetAttFill(1)->SetFillStyle(0);
        
        TGraphMultiErrors* Fig5PbPb = new TGraphMultiErrors("Fig5PbPb", "#it{#v}_{2,J/#psi} wrt #it{p}_{T}", 7, qbx, qby, qbexl, qbexh, qbeylstat, qbeyhstat);
            Fig5PbPb->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
        Fig5PbPb->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
           Fig5PbPb->AddYError(7, qbeylsys, qbeyhsys);
              Fig5PbPb->SetMarkerStyle(21);
           Fig5PbPb->SetMarkerColor(kBlack);
              Fig5PbPb->SetLineColorAlpha(46, 0.);
              Fig5PbPb->GetAttLine(0)->SetLineColor(kBlack);
              Fig5PbPb->GetAttLine(1)->SetLineColor(kBlack);
            Fig5PbPb->GetAttLine(0)->SetLineWidth(2);
            Fig5PbPb->GetAttLine(1)->SetLineWidth(2);
           Fig5PbPb->GetAttFill(1)->SetFillStyle(0);
            
        TGraphMultiErrors* Fig5pp = new TGraphMultiErrors("Fig5pp", "#it{#v}_{2,J/#psi} wrt #it{p}_{T}", 5, qcx, qcy, qaexl, qaexh, qceylstat, qceyhstat);
            Fig5pp->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
        Fig5pp->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
        Fig5pp->AddYError(5, qceylsys, qceyhsys);
           Fig5pp->SetMarkerStyle(20);
        Fig5pp->SetMarkerColor(kAzure-3);
           Fig5pp->SetLineColorAlpha(46, 0.);
           Fig5pp->GetAttLine(0)->SetLineColor(kAzure-3);
        Fig5pp->GetAttLine(1)->SetLineColor(kAzure-3);
            Fig5pp->GetAttLine(0)->SetLineWidth(2);
            Fig5pp->GetAttLine(1)->SetLineWidth(2);
         Fig5pp->GetAttFill(1)->SetFillStyle(0);
           
           Figsys->Draw("a p s ; z ; 5 s=0.5");
        Fig5PbPb->Draw("p s ; z ; 5 s=0.5");
        Fig5pp->Draw("p s ; z ; 5 s=0.5");
        
        TLine *ls=new TLine(0.,0.0,8.,0.0);
           ls->SetLineColor(kBlack);
           ls->SetLineWidth(1);
           ls->SetLineStyle(9);
           ls->Draw("same");
            
            TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
                      legendov2->SetFillColorAlpha(kWhite, 0.);
                      legendov2->SetBorderSize(0);
                       legendov2->SetTextFont(42);
                       legendov2->SetTextSize(0.05);
            Char_t messago[80];
            sprintf(messago,"ALICE Preliminary");
             legendov2->AddEntry(Fig5pp,messago,"");
                      legendov2->Draw();
            
            TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
                         legendov4->SetFillColorAlpha(kWhite, 0.);
                         legendov4->SetBorderSize(0);
                          legendov4->SetTextFont(62);
                          legendov4->SetTextSize(0.04);
            legendov4->SetTextColor(kAzure-3);
               sprintf(messago,"5.9%% global syst. uncertainty");
                legendov4->AddEntry(Fig5pp,messago,"");
                         legendov4->Draw();
        
        TLegend *legendsys=new TLegend(0.12,0.80,0.60,0.90);
          legendsys->SetFillColorAlpha(kWhite, 0.);
          legendsys->SetBorderSize(0);
            legendsys->SetTextFont(42);
          legendsys->SetTextSize(0.03);
        sprintf(message,"Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, (30-50%%) (JHEP 10 (2020) 141)");
        legendsys->AddEntry(Fig5PbPb,message);
            sprintf(message,"2.5 < #it{y}_{cms} < 4.0");
            legendsys->AddEntry(Fig5PbPb,message,"");
            sprintf(message,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02, 8.16 TeV, (0-20%%)-(40-100%%) (PLB 780 (2018) 7-20)");
            legendsys->AddEntry(Figsys,message);
            sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.03 < #it{y}_{cms} < 3.53");
            legendsys->AddEntry(Figsys,message,"");
        sprintf(message,"pp, #sqrt{#it{s}_{NN}} = 13 TeV, (0-5%%)-(40-100%%)");
        legendsys->AddEntry(Fig5pp,message);
            sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0, 2.5 < #it{y}_{cms} < 4.0");
            legendsys->AddEntry(Fig5pp,message,"");
          legendsys->Draw();
        
        
    }
    
    
    double aaold[]      = {-0.00584, 0.01802, -0.01266, 0.03064, -0.05892}; //10mrad

      double pp2[] = {-0.00705551,0.0222758,-0.0156563, 0.028918, -0.0591684};//5mrad
      double ppstat2[] = {0.00831353,0.0095233,0.0106528, 0.00952946,0.0147613};//5mrad
     // double ppsyst[] = {0.0099, 0.0077, 0.0140, 0.0085, 0.0258};
     // double ppsyst2[] = {0.0092, 0.0088, 0.0141, 0.0123, 0.0248};
    double ppsyst2[] = {0.0090, 0.0076, 0.0137, 0.0098, 0.0198}; //PYTHIA Half top
      
      for(int idx=0; idx<5; idx++){
          
          ppsyst2[idx]/=abs(pp2[idx]/aaold[idx]);
          pp2[idx]/=0.066;
          ppstat2[idx]/=0.066;
          ppsyst2[idx]/=0.066;
      }
      
      double pPb2[] = {-0.000504,-0.014926,0.031617, 0.082274, 0.033703};
      double pPbstat2[] = {0.018020,0.020871,0.023864, 0.021273,0.032535};
      double pPbsyst2[] = {0.007717, 0.011314, 0.013313, 0.011523, 0.017635};
      
      for(int idx=0; idx<5; idx++){
          pPb2[idx]/=0.068;
          pPbstat2[idx]/=0.068;
          pPbsyst2[idx]/=0.068;
      }
      
      double Pbp2[] = {0.022974,0.016185,0.073875, 0.069017, 0.078325};
      double Pbpstat2[] = {0.015113,0.018057,0.021716, 0.020426,0.034275};
      double Pbpsyst2[] = {0.014058, 0.013872, 0.014531, 0.014579, 0.021298};
      
      for(int idx=0; idx<5; idx++){
          Pbp2[idx]/=0.068;
          Pbpstat2[idx]/=0.068;
          Pbpsyst2[idx]/=0.068;
      }
      
      double avgmedium2[] = {0.022974,0.016185,0.073875, 0.069017, 0.078325};
      double avgmediumstat2[] = {0.015113,0.018057,0.021716, 0.020426,0.034275};
      double avgmediumsyst2[] = {0.014058, 0.013872, 0.014531, 0.014579, 0.021298};
      
      for(int idx=0; idx<5; idx++){
          avgmedium2[idx]=(Pbp2[idx]+pPb2[idx])/2.;
          avgmediumstat2[idx] = sqrt(pow(Pbpstat2[idx],2)+pow(pPbstat2[idx],2))/2.;
          avgmediumsyst2[idx] = sqrt(pow(Pbpsyst2[idx],2)+pow(pPbsyst2[idx],2))/2.;
      }
      
      double chi22=0;
      
      for(int idx=0; idx<5; idx++){
          chi22+=pow((pPb2[idx]-pp2[idx]),2)/(pow(ppstat2[idx],2)+pow(pPbstat2[idx],2)+pow(ppsyst2[idx],2)+pow(pPbsyst2[idx],2));
      }
      cout << "chi2 pp vs p-Pb: " << chi22 <<endl;
    //  cout << "Sigma deviation: " << sqrt(chi2/2) <<endl;
      cout << "Probability compatibility: " << TMath::Prob(chi22,5)/2.<<endl;
      cout << "Sigma deviation: " << 1.42*TMath::ErfInverse(1-(TMath::Prob(chi22,5))) <<endl;
      
      chi22=0;
      
      for(int idx=0; idx<5; idx++){
          chi22+=pow((Pbp2[idx]-pp2[idx]),2)/(pow(ppstat2[idx],2)+pow(Pbpstat2[idx],2)+pow(ppsyst2[idx],2)+pow(Pbpsyst2[idx],2));
      }
      cout << "chi2 pp vs Pb-p: " << chi22 <<endl;
     // cout << "Sigma deviation: " << sqrt(chi2/2) <<endl;
      cout << "Probability compatibility: " << TMath::Prob(chi22,5)/2.<<endl;
      cout << "Sigma deviation: " << 1.42*TMath::ErfInverse(1-(TMath::Prob(chi22,5))) <<endl;
      
      chi22=0;
      
      for(int idx=0; idx<5; idx++){
          chi22+=pow((avgmedium2[idx]-pp2[idx]),2)/(pow(ppstat2[idx],2)+pow(avgmediumstat2[idx],2)+pow(ppsyst2[idx],2)+pow(avgmediumsyst2[idx],2));
      }
      cout << "chi2 pp vs Pb-p+p-Pb avg.: " << chi22 <<endl;
     // cout << "Sigma deviation: " << sqrt(chi22/2) <<endl;
      cout << "Probability compatibility: " << TMath::Prob(chi22,5)/2.<<endl;
      cout << "Sigma deviation: " << 1.42*TMath::ErfInverse(1-(TMath::Prob(chi22,5))) <<endl;
      
      

    
       
        // v2 DIMU wrt Centrality V0M
        {
        
        double yax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
        double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
        double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
        double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
        double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
        double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
        double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
        double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

              auto cv0 = new TCanvas("cv0","cv0",200,10,800,600);
              double vax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
              double vay[]      = {0.01596, 0.01253, 0.04085, 0.02387, 0.01582, 0.03770, 0.00523};
              double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
              double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
              double* vaeylstat = new double[7]  {0.0078, 0.0070, 0.00824, 0.00655, 0.00801, 0.00937, 0.00961};
              double* vaeyhstat = new double[7]  {0.0078, 0.0070, 0.00824, 0.00655, 0.00801, 0.00937, 0.00961};
             // double* vaeylsys  = new double[6]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258,0.0234};
            double* vaeylsys  = new double[7] {0.0131, 0.0126, 0.0266, 0.0134, 0.0141, 0.0285, 0.0227}; //PYTHIA Full
              double* vaeyhsys  = new double[7] {0.0131, 0.0126, 0.0266, 0.0134, 0.0141, 0.0285, 0.0227};
            
//            double* vaeylsys  = new double[7] {0.0127, 0.0122, 0.0234, 0.0103, 0.0122, 0.0200, 0.0217}; //PYTHIA Half top bottom
//            double* vaeyhsys  = new double[7] {0.0128, 0.0122, 0.0236, 0.0105, 0.0123, 0.0207, 0.0218};
           
           double vbx[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
           double vby[]      = {0.01143, 0.00467, 0.04030, 0.01954, 0.01001, 0.03448, -0.00803};
           double* vbeylstat = new double[7]  {0.00829, 0.00783, 0.00875, 0.00746, 0.00946, 0.01174, 0.01320};
           double* vbeyhstat = new double[7]  {0.00829, 0.00783, 0.00875, 0.00746, 0.00946, 0.01174, 0.01320};
           //double* vbeylsys  = new double[6]  {0.0093, 0.0074, 0.0136, 0.0104, 0.0245,0.0277};
            double* vbeylsys  = new double[7]   {0.0109, 0.0138, 0.0263, 0.0073, 0.0099, 0.0232, 0.0252}; //PYTHIA Full
           double* vbeyhsys  = new double[7]   {0.0109, 0.0138, 0.0263, 0.0073, 0.0099, 0.0232, 0.0252};
            
//            double* vbeylsys  = new double[7]   {0.0108, 0.0138, 0.0250, 0.0062, 0.0096, 0.0210, 0.0250}; //PYTHIA Half top bottom
//            double* vbeyhsys  = new double[7]   {0.0108, 0.0138, 0.0251, 0.0063, 0.0096, 0.0211, 0.0250};
        
        TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt Centrality", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
               Figempty->GetXaxis()->SetTitle("Centrality - V0M (Percent)");
           Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
            Figempty->GetYaxis()->SetRangeUser(-0.15, 0.20);
            Figempty->GetXaxis()->SetRangeUser(0, 30.);
            gStyle->SetPadTickX(1);
            gStyle->SetPadTickY(1);
              Figempty->AddYError(7, yaeylsys, yaeyhsys);
        Figempty->GetXaxis()->SetRangeUser(0,30);
            
              TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "", 7, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);//#it{v}_{2,J/#psi} wrt #it{p}_{T}
               Figv0malice->GetXaxis()->SetTitle("Centrality - V0M (Percent)");
           Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
            Figv0malice->GetYaxis()->SetRangeUser(-0.15, 0.20);
            Figv0malice->GetXaxis()->SetRangeUser(0, 30.);
            gStyle->SetPadTickX(1);
            gStyle->SetPadTickY(1);
            gStyle->SetPadBottomMargin(0.15);
            gStyle->SetPadLeftMargin(0.15);
              Figv0malice->AddYError(7, vaeylsys, vaeyhsys);
              Figv0malice->SetMarkerStyle(8);
           Figv0malice->SetMarkerColor(kAzure-3);
              Figv0malice->SetLineColorAlpha(46, 0.);
              Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
              Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
        Figv0malice->GetAttLine(0)->SetLineWidth(2);
        Figv0malice->GetAttLine(1)->SetLineWidth(2);
           Figv0malice->GetAttFill(1)->SetFillStyle(0);
        Figv0malice->GetXaxis()->SetRangeUser(0,30);
           
           TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "", 7, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
               Figv0matlas->GetXaxis()->CenterTitle("Centrality - V0M (Percent)");
           Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
            Figv0matlas->GetYaxis()->SetRangeUser(-0.15, 0.20);
            Figv0matlas->GetXaxis()->SetRangeUser(0, 30.);
            gStyle->SetPadTickX(0);
            gStyle->SetPadTickY(0);
              Figv0matlas->AddYError(7, vbeylsys, vbeyhsys);
                 Figv0matlas->SetMarkerStyle(21);
            Figv0matlas->SetMarkerColor(kRed+1);
                 Figv0matlas->SetLineColorAlpha(46, 0.);
                 Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
                 Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
        Figv0matlas->GetAttLine(0)->SetLineWidth(2);
        Figv0matlas->GetAttLine(1)->SetLineWidth(2);
              Figv0matlas->GetAttFill(1)->SetFillStyle(0);
           
              Figempty->Draw("a p s ; z ; 5 s=0.5");
              Figv0malice->Draw("p s ; z ; 5 s=0.5");
           Figv0matlas->Draw("p s ; z ; 5 s=0.5");
           
           TLine *ls2=new TLine(0.,0.0,30.,0.0);
              ls2->SetLineColor(kBlack);
              ls2->SetLineWidth(1);
              ls2->SetLineStyle(9);
              ls2->Draw("same");
           
           TLegend *legendov=new TLegend(0.12,0.65,0.40,0.80);
                         legendov->SetFillColorAlpha(kWhite, 0.);
                         legendov->SetBorderSize(0);
                          legendov->SetTextFont(42);
                          legendov->SetTextSize(0.035);
               Char_t messagi[80];
              // sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
             //  legendov->AddEntry(Figv0malice,messagi,"");
                       legendov->AddEntry(Figv0malice,"Subtracted yields method");

                        legendov->AddEntry(Figv0matlas,"Template fit method");

                         legendov->Draw();
        
    //           TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
    //                         legendov2->SetFillColorAlpha(kWhite, 0.);
    //                         legendov2->SetBorderSize(0);
    //                          legendov2->SetTextFont(42);
    //                          legendov2->SetTextSize(0.05);
    //               Char_t messago[80];
    //               sprintf(messago,"ALICE Preliminary");
    //                legendov2->AddEntry(Figv0matlas,messago,"");
    //                         legendov2->Draw();
    //
    //           TLegend *legendov3=new TLegend(0.5,0.67,0.80,0.87);
    //                     legendov3->SetFillColorAlpha(kWhite, 0.);
    //                     legendov3->SetBorderSize(0);
    //                      legendov3->SetTextFont(42);
    //                      legendov3->SetTextSize(0.04);
    //           sprintf(messago,"V0M (0-5%%)-(40-100%%)");
    //            legendov3->AddEntry(Figv0matlas,messago,"");
    //           sprintf(messago,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
    //           legendov3->AddEntry(Figv0matlas,messago,"");
    //            sprintf(messago,"2.5 < #it{y}_{cms} < 4.0");
    //            legendov3->AddEntry(Figv0matlas,messago,"");
    //           sprintf(messago,"1.5 < |#it{#Delta#eta}| < 5.0");
    //           legendov3->AddEntry(Figv0matlas,messago,"");
    //                     legendov3->Draw();
    //
    //           TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
    //                     legendov4->SetFillColorAlpha(kWhite, 0.);
    //                     legendov4->SetBorderSize(0);
    //                      legendov4->SetTextFont(62);
    //                      legendov4->SetTextSize(0.04);
    //        legendov4->SetTextColor(kAzure-3);
    //           sprintf(messago,"5.9%% global syst. uncertainty");
    //            legendov4->AddEntry(Figv0matlas,messago,"");
    //                     legendov4->Draw();
    //
    //        TLegend *legendov42=new TLegend(0.08,0.12,0.40,0.17);
    //                  legendov42->SetFillColorAlpha(kWhite, 0.);
    //                  legendov42->SetBorderSize(0);
    //                   legendov42->SetTextFont(62);
    //                   legendov42->SetTextSize(0.04);
    //        legendov42->SetTextColor(kRed+1);
    //        sprintf(messago,"2.7%% global syst. uncertainty");
    //         legendov42->AddEntry(Figv0matlas,messago,"");
    //                  legendov42->Draw();
    //
    //
        }
        
        // v2 DIMU wrt Centrality SPDT
        {
            
            double yax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
            double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
            double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
            double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
            double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
            double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
            double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
            double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

                  auto cv = new TCanvas("cv","cv",200,10,800,600);
                  double vax[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
                  double vay[]      = {0.02236, 0.03450, 0.02851, 0.04882, 0.02843, 0.03669, 0.03467};
                  double vaexl[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
                  double vaexh[]    = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
                  double* vaeylstat = new double[7]  {0.00594, 0.00554, 0.00658, 0.00513, 0.00628, 0.00762, 0.00716};
                  double* vaeyhstat = new double[7]  {0.00594, 0.00554, 0.00658, 0.00513, 0.00628, 0.00762, 0.00716};
                 // double* vaeylsys  = new double[6]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258,0.0234};
         //      double* vaeylsys  = new double[7] {0.0162, 0.0192, 0.0177, 0.0328, 0.0275, 0.0499, 0.0550}; //PYTHIA Full
           //       double* vaeyhsys  = new double[7] {0.0162, 0.0192, 0.0177, 0.0328, 0.0275, 0.0499, 0.0550};
            
            double* vaeylsys  = new double[7] {0.0124, 0.0077, 0.0094, 0.0173, 0.0106, 0.0106, 0.0117}; //PYTHIA Half top bottom
            double* vaeyhsys  = new double[7] {0.0127, 0.0089, 0.0101, 0.0187, 0.0119, 0.0135, 0.0158};
               
               double vbx[]      = {0.5, 2.0, 4.0, 7.5, 12.5, 17.5, 25.};
               double vby[]      = {0.00425, 0.01300, 0.00767, 0.02202, 0.00406, 0.00806, 0.00398};
               double* vbeylstat = new double[7]  {0.00515, 0.00479, 0.00560, 0.00487, 0.00580, 0.00746, 0.00793};
               double* vbeyhstat = new double[7]  {0.00515, 0.00479, 0.00560, 0.00487, 0.00580, 0.00746, 0.00793};
               //double* vbeylsys  = new double[6]  {0.0093, 0.0074, 0.0136, 0.0104, 0.0245,0.0277};
       //         double* vbeylsys  = new double[7]   {0.0308, 0.0047, 0.0050, 0.0124, 0.0071, 0.0049, 0.0059}; //PYTHIA Full
         //      double* vbeyhsys  = new double[7]   {0.0308, 0.0047, 0.0050, 0.0124, 0.0071, 0.0049, 0.0059};
            
            double* vbeylsys  = new double[7]   {0.0308, 0.0040, 0.0048, 0.0116, 0.0071, 0.0044, 0.0058}; //PYTHIA Half top bottom
            double* vbeyhsys  = new double[7]   {0.0308, 0.0040, 0.0048, 0.0116, 0.0071, 0.0044, 0.0058};
            
            TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt Centrality", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
                   Figempty->GetXaxis()->SetTitle("Centrality - SPDTracklets (Percent)");
               Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
                Figempty->GetYaxis()->SetRangeUser(-0.15, 0.20);
                Figempty->GetXaxis()->SetRangeUser(0, 30.);
                gStyle->SetPadTickX(1);
                gStyle->SetPadTickY(1);
                  Figempty->AddYError(7, yaeylsys, yaeyhsys);
            Figempty->GetXaxis()->SetRangeUser(0,30);
                
                  TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "", 7, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);//#it{v}_{2,J/#psi} wrt #it{p}_{T}
                   Figv0malice->GetXaxis()->SetTitle("Centrality - SPDTracklets (Percent)");
               Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
                Figv0malice->GetYaxis()->SetRangeUser(-0.15, 0.20);
                Figv0malice->GetXaxis()->SetRangeUser(0, 30.);
                gStyle->SetPadTickX(1);
                gStyle->SetPadTickY(1);
                gStyle->SetPadBottomMargin(0.15);
                gStyle->SetPadLeftMargin(0.15);
                  Figv0malice->AddYError(7, vaeylsys, vaeyhsys);
                  Figv0malice->SetMarkerStyle(8);
               Figv0malice->SetMarkerColor(kAzure-3);
                  Figv0malice->SetLineColorAlpha(46, 0.);
                  Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
                  Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
            Figv0malice->GetAttLine(0)->SetLineWidth(2);
            Figv0malice->GetAttLine(1)->SetLineWidth(2);
               Figv0malice->GetAttFill(1)->SetFillStyle(0);
            Figv0malice->GetXaxis()->SetRangeUser(0,30);
               
               TGraphMultiErrors* Figv0matlas = new TGraphMultiErrors("Figv0matlas", "", 7, vbx, vby, vaexl, vaexh, vbeylstat, vbeyhstat);
                   Figv0matlas->GetXaxis()->CenterTitle("Centrality - V0M (Percent)");
               Figv0matlas->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
                Figv0matlas->GetYaxis()->SetRangeUser(-0.15, 0.20);
                Figv0matlas->GetXaxis()->SetRangeUser(0, 30.);
                gStyle->SetPadTickX(0);
                gStyle->SetPadTickY(0);
                  Figv0matlas->AddYError(7, vbeylsys, vbeyhsys);
                     Figv0matlas->SetMarkerStyle(21);
                Figv0matlas->SetMarkerColor(kRed+1);
                     Figv0matlas->SetLineColorAlpha(46, 0.);
                     Figv0matlas->GetAttLine(0)->SetLineColor(kRed+1);
                     Figv0matlas->GetAttLine(1)->SetLineColor(kRed+1);
            Figv0matlas->GetAttLine(0)->SetLineWidth(2);
            Figv0matlas->GetAttLine(1)->SetLineWidth(2);
                  Figv0matlas->GetAttFill(1)->SetFillStyle(0);
               
                  Figempty->Draw("a p s ; z ; 5 s=0.5");
                  Figv0malice->Draw("p s ; z ; 5 s=0.5");
               Figv0matlas->Draw("p s ; z ; 5 s=0.5");
               
               TLine *ls2=new TLine(0.,0.0,30.,0.0);
                  ls2->SetLineColor(kBlack);
                  ls2->SetLineWidth(1);
                  ls2->SetLineStyle(9);
                  ls2->Draw("same");
               
               TLegend *legendov=new TLegend(0.12,0.65,0.40,0.80);
                             legendov->SetFillColorAlpha(kWhite, 0.);
                             legendov->SetBorderSize(0);
                              legendov->SetTextFont(42);
                              legendov->SetTextSize(0.035);
                   Char_t messagi[80];
                  // sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
                 //  legendov->AddEntry(Figv0malice,messagi,"");
                           legendov->AddEntry(Figv0malice,"Subtracted yields method");

                            legendov->AddEntry(Figv0matlas,"Template fit method");

                             legendov->Draw();
            
        //           TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
        //                         legendov2->SetFillColorAlpha(kWhite, 0.);
        //                         legendov2->SetBorderSize(0);
        //                          legendov2->SetTextFont(42);
        //                          legendov2->SetTextSize(0.05);
        //               Char_t messago[80];
        //               sprintf(messago,"ALICE Preliminary");
        //                legendov2->AddEntry(Figv0matlas,messago,"");
        //                         legendov2->Draw();
        //
        //           TLegend *legendov3=new TLegend(0.5,0.67,0.80,0.87);
        //                     legendov3->SetFillColorAlpha(kWhite, 0.);
        //                     legendov3->SetBorderSize(0);
        //                      legendov3->SetTextFont(42);
        //                      legendov3->SetTextSize(0.04);
        //           sprintf(messago,"V0M (0-5%%)-(40-100%%)");
        //            legendov3->AddEntry(Figv0matlas,messago,"");
        //           sprintf(messago,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
        //           legendov3->AddEntry(Figv0matlas,messago,"");
        //            sprintf(messago,"2.5 < #it{y}_{cms} < 4.0");
        //            legendov3->AddEntry(Figv0matlas,messago,"");
        //           sprintf(messago,"1.5 < |#it{#Delta#eta}| < 5.0");
        //           legendov3->AddEntry(Figv0matlas,messago,"");
        //                     legendov3->Draw();
        //
        //           TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
        //                     legendov4->SetFillColorAlpha(kWhite, 0.);
        //                     legendov4->SetBorderSize(0);
        //                      legendov4->SetTextFont(62);
        //                      legendov4->SetTextSize(0.04);
        //        legendov4->SetTextColor(kAzure-3);
        //           sprintf(messago,"5.9%% global syst. uncertainty");
        //            legendov4->AddEntry(Figv0matlas,messago,"");
        //                     legendov4->Draw();
        //
        //        TLegend *legendov42=new TLegend(0.08,0.12,0.40,0.17);
        //                  legendov42->SetFillColorAlpha(kWhite, 0.);
        //                  legendov42->SetBorderSize(0);
        //                   legendov42->SetTextFont(62);
        //                   legendov42->SetTextSize(0.04);
        //        legendov42->SetTextColor(kRed+1);
        //        sprintf(messago,"2.7%% global syst. uncertainty");
        //         legendov42->AddEntry(Figv0matlas,messago,"");
        //                  legendov42->Draw();
        //
        //
            }
    
    
    
    {
        double ppc[] = {0.01596,0.01253,0.04085, 0.02387, 0.01582, 0.03770, 0.00523};
          double ppstatc[] = {0.0078,0.007,0.00824, 0.00655,0.00801, 0.00937, 0.00961};

        double ppsystc[] = {0.0126, 0.0121, 0.0226, 0.0096, 0.0119, 0.0187, 0.0217};

          
          double chi22c=0;
          
          for(int idx=0; idx<7; idx++){
              chi22c+=pow(ppc[idx],2)/(pow(ppstatc[idx],2)+pow(ppsystc[idx],2));
          }
          cout << "chi2 pp cent : " << chi22c <<endl;
        //  cout << "Sigma deviation: " << sqrt(chi2/2) <<endl;
          cout << "Probability compatibility: " << TMath::Prob(chi22c,7)/2.<<endl;
          cout << "Sigma deviation: " << 1.42*TMath::ErfInverse(1-(TMath::Prob(chi22c,7))) <<endl;
          
          
    }
    
    //Final v2 V0M + CMS
    {
         
         double yax[]      = {1, 2.5, 3.5, 5, 7, 10, 12};
         double yay[]      = {0.,0.,0.,0.,0.,0.,0.};
         double yaexl[]    = {0.,0.,0.,0.,0.,0.,0.};
         double yaexh[]    = {0.,0.,0.,0.,0.,0.,0.};
         double* yaeylstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
         double* yaeyhstat = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
         double* yaeylsys  = new double[7]  {0.,0.,0.,0.,0.,0.,0.};
         double* yaeyhsys  = new double[7]   {0.,0.,0.,0.,0.,0.,0.};

               auto cv0mpropcms = new TCanvas("cv0mpropcms","cv0mpropcms",200,10,800,600);
               double vax[]      = {1, 2.5, 3.5, 5, 7, 10};
               double vay[]      = {-0.00584,0.01802,-0.01266,0.03064,-0.05892,0.03671};
               double vaexl[]    = {1., 0.5, 0.5, 1., 1., 2.};
               double vaexh[]    = {1., 0.5, 0.5, 1., 1., 2.};
               double* vaeylstat = new double[6]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503,0.02098};
               double* vaeyhstat = new double[6]  {0.00842, 0.00965, 0.01082, 0.00967, 0.01503,0.02098};
              // double* vaeylsys  = new double[6]  {0.0099, 0.0077, 0.0140, 0.0085, 0.0258,0.0234};
             double* vaeylsys  = new double[6] {0.0092, 0.0088, 0.0141, 0.0123, 0.0248, 0.0253}; //PYTHIA Full
               double* vaeyhsys  = new double[6] {0.0092, 0.0088, 0.0141, 0.0123, 0.0248, 0.0253};
             
//             double* vaeylsys  = new double[6] {0.0090, 0.0076, 0.0137, 0.0096, 0.0201, 0.0235}; //PYTHIA Half top bottom
//             double* vaeyhsys  = new double[6] {0.0090, 0.0076, 0.0137, 0.0098, 0.0198, 0.0236};
//
            double vbx[]      = {3., 5., 7.};
            double vby[]      = {0.061, 0.04, -0.043};
        double vaexlD0[]    = {1.,1.,1.};
        double vaexhD0[]    = {1.,1.,1.};
            double* vbeylstat = new double[3]  {0.018, 0.02, 0.031};
            double* vbeyhstat = new double[3]  {0.018, 0.02, 0.031};
            //double* vbeylsys  = new double[6]  {0.0093, 0.0074, 0.0136, 0.0104, 0.0245,0.0277};
    //         double* vbeylsys  = new double[6]   {0.0087, 0.0075, 0.0133, 0.0124, 0.0186,0.0271}; PYTHIA Full
      //      double* vbeyhsys  = new double[6]   {0.0087, 0.0075, 0.0133, 0.0124, 0.0186,0.0271};
             
             double* vbeylsys  = new double[3]  {0.013, 0.015, 0.014};
             double* vbeyhsys  = new double[3]  {0.013, 0.015, 0.014};
        
        
        double vbxh[]      = {0.296, 0.494, 0.694, 0.895, 1.18, 1.58, 1.98, 2.46, 3.14, 4.02, 5.17};
                double vbyh[]     = {0.0207, 0.0303, 0.0385, 0.0468, 0.0577, 0.0695, 0.0789, 0.0845, 0.0929, 0.0867, 0.0703};
            double vaexlh[]    = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            double vaexhh[]    = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
                double* vbeylstath = new double[11]  {0.000279175, 0.000307433, 0.00036263, 0.000433623, 0.000400655, 0.000559635, 0.000771996, 0.000917843, 0.00124438, 0.00179962, 0.00263416};
                double* vbeyhstath = new double[11]  {0.000279175, 0.000307433, 0.00036263, 0.000433623, 0.000400655, 0.000559635, 0.000771996, 0.000917843, 0.00124438, 0.00179962, 0.00263416};
                 
                 double* vbeylsysh  = new double[11]  {0.00205845, 0.00301128, 0.00383291, 0.00465451, 0.00574147, 0.00691361, 0.00784575, 0.0115267, 0.0126731, 0.0170246, 0.0260422};
                 double* vbeyhsysh  = new double[11]  {0.00205845, 0.00301128, 0.00383291, 0.00465451, 0.00574147, 0.00691361, 0.00784575, 0.0115267, 0.0126731, 0.0170246, 0.0260422};
        
        
         
         TGraphMultiErrors* Figempty = new TGraphMultiErrors("Figempty", "#it{v}_{2,J/#psi} wrt #it{p}_{T}", 7, yax, yay, yaexl, yaexh, yaeylstat, yaeyhstat);
                Figempty->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            Figempty->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
             Figempty->GetYaxis()->SetRangeUser(-0.15, 0.20);
             Figempty->GetXaxis()->SetRangeUser(0, 12.);
             gStyle->SetPadTickX(1);
             gStyle->SetPadTickY(1);
               Figempty->AddYError(7, yaeylsys, yaeyhsys);
         Figempty->GetXaxis()->SetRangeUser(0,15);
             
               TGraphMultiErrors* Figv0malice = new TGraphMultiErrors("Figv0malice", "", 6, vax, vay, vaexl, vaexh, vaeylstat, vaeyhstat);//#it{v}_{2,J/#psi} wrt #it{p}_{T}
                Figv0malice->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
            Figv0malice->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
             Figv0malice->GetYaxis()->SetRangeUser(-0.15, 0.20);
             Figv0malice->GetXaxis()->SetRangeUser(0, 12.);
             gStyle->SetPadTickX(1);
             gStyle->SetPadTickY(1);
             gStyle->SetPadBottomMargin(0.15);
             gStyle->SetPadLeftMargin(0.15);
               Figv0malice->AddYError(6, vaeylsys, vaeyhsys);
               Figv0malice->SetMarkerStyle(8);
            Figv0malice->SetMarkerColor(kAzure-3);
               Figv0malice->SetLineColor(kAzure-3);
        Figv0malice->SetLineWidth(2);
               Figv0malice->GetAttLine(0)->SetLineColor(kAzure-3);
               Figv0malice->GetAttLine(1)->SetLineColor(kAzure-3);
         Figv0malice->GetAttLine(0)->SetLineWidth(2);
         Figv0malice->GetAttLine(1)->SetLineWidth(2);
            Figv0malice->GetAttFill(1)->SetFillStyle(0);
         Figv0malice->GetXaxis()->SetRangeUser(0,15);
            
            TGraphMultiErrors* Figv0D0 = new TGraphMultiErrors("Figv0D0", "", 3, vbx, vby, vaexlD0, vaexhD0, vbeylstat, vbeyhstat);
                Figv0D0->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
            Figv0D0->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
             Figv0D0->GetYaxis()->SetRangeUser(-0.15, 0.20);
             Figv0D0->GetXaxis()->SetRangeUser(0, 12.);
             gStyle->SetPadTickX(0);
             gStyle->SetPadTickY(0);
               Figv0D0->AddYError(3, vbeylsys, vbeyhsys);
                  Figv0D0->SetMarkerStyle(21);
             Figv0D0->SetMarkerColor(kRed+1);
                  Figv0D0->SetLineColor(kRed+1);
                  Figv0D0->GetAttLine(0)->SetLineColor(kRed+1);
                  Figv0D0->GetAttLine(1)->SetLineColor(kRed+1);
        Figv0D0->SetLineWidth(2);
         Figv0D0->GetAttLine(0)->SetLineWidth(2);
         Figv0D0->GetAttLine(1)->SetLineWidth(2);
               Figv0D0->GetAttFill(1)->SetFillStyle(0);
        
        TGraphMultiErrors* Figv0CMSH = new TGraphMultiErrors("Figv0CMSH", "", 11, vbxh, vbyh, vaexlh, vaexhh, vbeylstath, vbeyhstath);
                Figv0CMSH->GetXaxis()->CenterTitle("#it{p}_{T} (GeV/#it{c})");
            Figv0CMSH->GetYaxis()->CenterTitle("#it{v}_{2,J/#psi}");
             Figv0CMSH->GetYaxis()->SetRangeUser(-0.15, 0.20);
             Figv0CMSH->GetXaxis()->SetRangeUser(0, 12.);
             gStyle->SetPadTickX(0);
             gStyle->SetPadTickY(0);
               Figv0CMSH->AddYError(11, vbeylsysh, vbeyhsysh);
                  Figv0CMSH->SetMarkerStyle(21);
             Figv0CMSH->SetMarkerColor(kGreen+1);
                  Figv0CMSH->SetLineColorAlpha(46,0);
                  Figv0CMSH->GetAttLine(0)->SetLineColor(kGreen+1);
                  Figv0CMSH->GetAttLine(1)->SetLineColor(kGreen+1);
        Figv0CMSH->SetLineWidth(2);
         Figv0CMSH->GetAttLine(0)->SetLineWidth(2);
         Figv0CMSH->GetAttLine(1)->SetLineWidth(2);
               Figv0CMSH->GetAttFill(1)->SetFillStyle(0);
            
               Figempty->Draw("a p s ; z ; 5 z s=0.2");
               Figv0malice->Draw("p s ; z ; 5 z s=0.2");
            Figv0D0->Draw("p s ; z ; 5 z s=0.2");
        Figv0CMSH->Draw("p s ; z ; 5 z s=0.2");
            
            TLine *ls2=new TLine(0.,0.0,12.,0.0);
               ls2->SetLineColor(kBlack);
               ls2->SetLineWidth(1);
               ls2->SetLineStyle(9);
               ls2->Draw("same");
            
            TLegend *legendov=new TLegend(0.12,0.65,0.40,0.80);
                          legendov->SetFillColorAlpha(kWhite, 0.);
                          legendov->SetBorderSize(0);
                           legendov->SetTextFont(42);
                           legendov->SetTextSize(0.035);
                Char_t messagi[80];
               // sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
              //  legendov->AddEntry(Figv0malice,messagi,"");
                        legendov->AddEntry(Figv0malice,"J/#psi - This analysis, 2.5 < #it{y}_{cms} < 4.0, V0M (0-5%)-(40-100%)");

                         legendov->AddEntry(Figv0D0,"Prompt D^{0} - CMS, |#it{y}_{cms}| < 1.0, N_{trk} > 100");
        
        legendov->AddEntry(Figv0CMSH,"Charged hadrons - CMS, |#it{#Delta#eta}| > 2.0, (150 > N_{trk} > 105) - (20 > N_{trk} > 10)");

                          legendov->Draw();
         
//                TLegend *legendov2=new TLegend(0.12,0.75,0.40,0.90);
//                              legendov2->SetFillColorAlpha(kWhite, 0.);
//                              legendov2->SetBorderSize(0);
//                               legendov2->SetTextFont(42);
//                               legendov2->SetTextSize(0.05);
//                    Char_t messago[80];
//                    sprintf(messago,"ALICE Preliminary");
//                     legendov2->AddEntry(Figv0matlas,messago,"");
//                              legendov2->Draw();
//
//                TLegend *legendov3=new TLegend(0.5,0.67,0.80,0.87);
//                          legendov3->SetFillColorAlpha(kWhite, 0.);
//                          legendov3->SetBorderSize(0);
//                           legendov3->SetTextFont(42);
//                           legendov3->SetTextSize(0.04);
//                sprintf(messago,"V0M (0-5%%)-(40-100%%)");
//                 legendov3->AddEntry(Figv0matlas,messago,"");
//                sprintf(messago,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
//                legendov3->AddEntry(Figv0matlas,messago,"");
//                 sprintf(messago,"2.5 < #it{y}_{cms} < 4.0");
//                 legendov3->AddEntry(Figv0matlas,messago,"");
//                sprintf(messago,"1.5 < |#it{#Delta#eta}| < 5.0");
//                legendov3->AddEntry(Figv0matlas,messago,"");
//                          legendov3->Draw();
                
//                TLegend *legendov4=new TLegend(0.08,0.17,0.40,0.22);
//                          legendov4->SetFillColorAlpha(kWhite, 0.);
//                          legendov4->SetBorderSize(0);
//                           legendov4->SetTextFont(62);
//                           legendov4->SetTextSize(0.04);
//             legendov4->SetTextColor(kAzure-3);
//                sprintf(messago,"5.9%% global syst. uncertainty");
//                 legendov4->AddEntry(Figv0matlas,messago,"");
//                          legendov4->Draw();
//
//             TLegend *legendov42=new TLegend(0.08,0.12,0.40,0.17);
//                       legendov42->SetFillColorAlpha(kWhite, 0.);
//                       legendov42->SetBorderSize(0);
//                        legendov42->SetTextFont(62);
//                        legendov42->SetTextSize(0.04);
//             legendov42->SetTextColor(kRed+1);
//             sprintf(messago,"2.7%% global syst. uncertainty");
//              legendov42->AddEntry(Figv0matlas,messago,"");
//                       legendov42->Draw();
//
         
         }

}
