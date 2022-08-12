// Macro qui plot la factorisation dimu


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

void FactorisationDimu();

class CSVRow
{
    public:
        std::string_view operator[](std::size_t index) const
        {
            return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(',', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

void FactorisationDimu(){
    
    
    
    
    CSVRow row;
     
     Char_t FolderCSVName[500];
     Char_t SystematicsFileName[500];
    const int numsituations = 6;
         sprintf(FolderCSVName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrierDimu/NewAnalysisAllEst_Run2_V0MPercentile_0-5_40-100_pt0-2-3-4-6-8-12_DPhiCut1mrad");
         cout<<"FolderCSVName"<<FolderCSVName<<endl;
     
     
     sprintf(SystematicsFileName,"%s/SystematicsFile.csv", FolderCSVName);
     cout<<"SystematicsFileName"<<SystematicsFileName<<endl;
     std::ifstream filecsv(SystematicsFileName);
     
    Double_t v2pPb[numsituations];
    Double_t errv2pPb[numsituations];
    
     while(filecsv >> row)
     {
         
         for(int ptbin=0; ptbin<numsituations; ptbin++){
             v2pPb[ptbin] = stod(std::string(row[28+(18*ptbin)]));
             errv2pPb[ptbin] = stod(std::string(row[29+(18*ptbin)]));
         }
         
     }
    
    
    sprintf(FolderCSVName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrierDimu/NewAnalysisAllEst_Run2_V0MPercentile_0-5_40-100_pt0-2-3-4-6-8-12_DPhiCut2mrad");
         cout<<"FolderCSVName"<<FolderCSVName<<endl;
     
     
     sprintf(SystematicsFileName,"%s/SystematicsFile.csv", FolderCSVName);
     cout<<"SystematicsFileName"<<SystematicsFileName<<endl;
     std::ifstream filecsv2(SystematicsFileName);
     
    Double_t v2pPb2[numsituations];
    Double_t errv2pPb2[numsituations];
    
     while(filecsv2 >> row)
     {
         
         for(int ptbin=0; ptbin<numsituations; ptbin++){
             v2pPb2[ptbin] = stod(std::string(row[28+(18*ptbin)]));
             errv2pPb2[ptbin] = stod(std::string(row[29+(18*ptbin)]));
         }
         
     }
    
    sprintf(FolderCSVName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrierDimu/NewAnalysisAllEst_Run2_V0MPercentile_0-5_40-100_pt0-2-3-4-6-8-12_DPhiCut5mrad");
         cout<<"FolderCSVName"<<FolderCSVName<<endl;
     
     
     sprintf(SystematicsFileName,"%s/SystematicsFile.csv", FolderCSVName);
     cout<<"SystematicsFileName"<<SystematicsFileName<<endl;
     std::ifstream filecsv3(SystematicsFileName);
     
    Double_t v2pPb3[numsituations];
    Double_t errv2pPb3[numsituations];
    
     while(filecsv3 >> row)
     {
         
         for(int ptbin=0; ptbin<numsituations; ptbin++){
             v2pPb3[ptbin] = stod(std::string(row[28+(18*ptbin)]));
             errv2pPb3[ptbin] = stod(std::string(row[29+(18*ptbin)]));
         }
         
     }
    
    sprintf(FolderCSVName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrierDimu/NewAnalysisAllEst_Run2_V0MPercentile_0-5_40-100_pt0-2-3-4-6-8-12");
         cout<<"FolderCSVName"<<FolderCSVName<<endl;
     
     
     sprintf(SystematicsFileName,"%s/SystematicsFile.csv", FolderCSVName);
     cout<<"SystematicsFileName"<<SystematicsFileName<<endl;
     std::ifstream filecsv4(SystematicsFileName);
     
    Double_t v2pPb4[numsituations];
    Double_t errv2pPb4[numsituations];
    
     while(filecsv4 >> row)
     {
         
         for(int ptbin=0; ptbin<numsituations; ptbin++){
             v2pPb4[ptbin] = stod(std::string(row[28+(18*ptbin)]));
             errv2pPb4[ptbin] = stod(std::string(row[29+(18*ptbin)]));
         }
         
     }
    
        // Figure V0M
        
        TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);

    Double_t pT[numsituations] = {1,2.5,3.5,5,7,10};
           Double_t errsit[numsituations] = {0,0,0,0,0,0};

    //1mrad

           TGraphErrors *v2syst = new TGraphErrors(numsituations,pT,v2pPb,errsit,errv2pPb);
            v2syst->SetTitle("v_{2,J/#psi} for different reference tracklet cuts");
            v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.,0.15);
            v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");

           v2syst->SetMarkerStyle(21);
            v2syst->SetMarkerColor(kRed);
            v2syst->SetLineColor(kRed);
           v2syst->Draw("AP");
            
            
            
            // Associate cut 2mrad
            
            
            for (Int_t i=0;i<numsituations;i++){
            pT[i]+=0.05;
            };

            TGraphErrors *v2syst2 = new TGraphErrors(numsituations,pT,v2pPb2,errsit,errv2pPb2);
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
                pT[i]+=0.05;
                };

                TGraphErrors *v2syst3 = new TGraphErrors(numsituations,pT,v2pPb3,errsit,errv2pPb3);
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
        
        // Associate cut 10mrad
            
            
            for (Int_t i=0;i<numsituations;i++){
            pT[i]+=0.05;
            };


            TGraphErrors *v2syst4 = new TGraphErrors(numsituations,pT,v2pPb4,errsit,errv2pPb4);
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
        
            
            TLegend *legend5=new TLegend(0.70,0.70,0.90,0.90);
            legend5->SetHeader("V0M");
              legend5->SetFillColorAlpha(kWhite, 0.);
              legend5->SetBorderSize(0);
                legend5->SetTextFont(42);
              legend5->SetTextSize(0.03);
                Char_t message[80];
                sprintf(message,"1mrad");
                legend5->AddEntry(v2syst,message);
           sprintf(message,"2mrad");
            legend5->AddEntry(v2syst2,message);
            sprintf(message,"5mrad");
            legend5->AddEntry(v2syst3,message);
        sprintf(message,"10mrad");
        legend5->AddEntry(v2syst4,message);
              legend5->Draw();
    
    
}
