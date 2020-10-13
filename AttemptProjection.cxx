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
# include <iostream>
# include "TH1F.h"
# include "TH1D.h"
# include "TH2F.h"

void AttemptProjection(){
    TH2F* hsource(NULL);
    hsource = new TH2F("hsource","Histo source",5,0,5,10,0,10);
    hsource->SetXTitle("source axe x");
    for(int i=0; i<6; i++){
        for(int j=0; j<11; j++){
            hsource->SetBinContent(i,j,i+j);
        }
    }
    TH1D* hprojection = hsource->ProjectionX("hprojection",0,-1,"e");
    for(int i=0; i<6; i++){
        hprojection->SetBinContent(i, hprojection->GetBinContent(i)/10);
    }
    
    TCanvas* c10 = new TCanvas;
    c10->Divide(1,2);
    c10->cd(1);
    hsource->Draw("colz");
    c10->cd(2);
    hprojection->Draw();
}
