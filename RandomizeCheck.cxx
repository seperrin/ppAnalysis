#include <TBrowser.h>
#include <TBufferFile.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TMath.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TLorentzVector.h>

# include "TF1.h"
# include "TF2.h"
# include "TProfile.h"
# include "TVirtualFitter.h"
# include "TFitter.h"
# include "TMinuit.h"
# include "TRandom.h"
# include "TObjArray.h"
# include <iostream>
# include "TH1F.h"
# include "TH1D.h"
# include "TH2F.h"

void RandomizeCheck(){
    TH1F* hsource(NULL);
    TrackletLight* trac1(NULL);
    TrackletLight* trac2(NULL);
    TrackletLight* tracklet(NULL);
    TClonesArray* fTracklets(NULL);
    
    double finalmean[6] = {NULL};
    double finalmeanerror[6] = {NULL};
    double rms[6] = {NULL};
    double nrand[6] = {1,2,3,5,10,20};
    
    for(int values=0; values<6; values++){
        
        double means[10] = {NULL};
        
        for(int trial = 0; trial<10;trial++){
            
            hsource = new TH1F("hsource","Histo source",201,-100.5,100.5);
            hsource->SetXTitle("source axe x");
            
        
            for(int loop=0; loop<10000;loop++){

                fTracklets = new TClonesArray("TrackletLight", 0);

                for (Int_t itrklt = 0; itrklt < 100; itrklt++) {
                    TClonesArray &tracklets = *fTracklets;
                    TrackletLight *tracklet = new (tracklets[itrklt]) TrackletLight();
                   tracklet->fEta = itrklt+1;
                }

                fTracklets->Randomize(nrand[values]);

                for(int idx_1=0; idx_1<100; idx_1++){
                    for(int idx_2=idx_1+1; idx_2<100; idx_2++){
                        trac1 = (TrackletLight*)fTracklets->At(idx_1);
                        trac2 = (TrackletLight*)fTracklets->At(idx_2);
                        hsource->Fill(trac2->fEta-trac1->fEta);
                    }
                }
                
            }
            TCanvas* c10 = new TCanvas;
            hsource->Draw();
            finalmean[values] += hsource->GetMean(1);
            finalmeanerror[values] += hsource->GetMeanError(1);
            means[trial] = hsource->GetMean(1);
            hsource->Clear();
        }
        
        finalmean[values] /= 10;
        finalmeanerror[values] /= 10;
        for(int idx=0; idx<10; idx++){
            rms[values] = pow(means[idx]-finalmean[values],2);
        }
        rms[values] = sqrt(rms[values]);
        
        cout << "Final mean: "<<finalmean[values] << " +/- " <<finalmeanerror[values]<< " +/- " <<rms[values]<<endl;
        
    }
    
    
    
    
    
    auto c5 = new TCanvas("c5","c5",200,10,600,400);
       double aexl[]    = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
     
       TGraphMultiErrors* Fig5 = new TGraphMultiErrors("Fig_Randomize", "Mean wrt number of randomizings", 6, nrand, finalmean, aexl, aexl, finalmeanerror, finalmeanerror);
        Fig5->GetXaxis()->SetTitle("Number of Randomize()");
    Fig5->GetYaxis()->SetTitle("Mean of correlations");
       Fig5->AddYError(6, rms, rms);
       Fig5->SetMarkerStyle(20);
    Fig5->SetMarkerColor(kBlue);
       Fig5->SetLineColorAlpha(46, 0.);
       Fig5->GetAttLine(0)->SetLineColor(kBlue);
       Fig5->GetAttLine(1)->SetLineColor(kBlue);
    Fig5->GetAttFill(1)->SetFillStyle(0);
    
    Fig5->Draw("a p s ; ; 5 s=0.5");
    
    
    
}

