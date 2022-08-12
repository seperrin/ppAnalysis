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
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

//#include "TTreeReader.h"
//#include "Event.h"

# include "TF1.h"
# include "TF2.h"
# include "TProfile.h"
# include "TVirtualFitter.h"
# include "TBackCompFitter.h"
# include "TGraphErrors.h"
# include "TFitter.h"
# include "TMinuit.h"
# include "TRandom.h"
# include <iostream>
# include <fstream>
# include <string>
# include "TH1F.h"
# include "TH1D.h"
# include "TH2F.h"
# include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"

void PetitTest();
void AutreTest();

int pouet = 1;

TH2F* deuxD(NULL);
TH1D* heinD(NULL);
TH1F* unD(NULL);
TH1F* unD2(NULL);

void PetitTest(){
    TH1::SetDefaultSumw2();
    deuxD = new TH2F("deuxD",
                     "deuxD",
                     2,0,2, 5, 0, 5);
    deuxD->Fill(0.5,0.5);
    deuxD->Fill(0.5,0.5);
    deuxD->Fill(1.5,0.5);
    deuxD->Fill(1.5,0.5);
    deuxD->Fill(1.5,0.5);
    deuxD->SetBinError(1,1,5);
    
    heinD = (TH1D*)(deuxD->ProjectionY("aaa",1,2,"e"));
    
    TCanvas* ccan=new TCanvas();
    ccan->Divide(1,3);
    ccan->cd(1);
    deuxD->DrawCopy("colz");
    ccan->cd(2);
    heinD->DrawCopy();
    ccan->cd(3);
    heinD = (TH1D*)(deuxD->ProjectionY("_py",1,2));
    heinD->DrawCopy();
}

void PetitTestAdd(){
    unD = new TH1F("unD",
                     "unD",
                     2,0,2);
    for(int i=0; i<5;i++){
        unD->Fill(0.5);
    }
    for(int i=0; i<2;i++){
        unD->Fill(1.5);
    }
    
    cout << "Histo 1:"<<endl;
    cout << "Bin 1: Contenu " << unD->GetBinContent(1) <<" Erreur " << unD->GetBinError(1) <<" Erreur relative " << unD->GetBinError(1)/unD->GetBinContent(1) <<endl;
    cout << "Bin 2: Contenu " << unD->GetBinContent(2) <<" Erreur " << unD->GetBinError(2) <<" Erreur relative " << unD->GetBinError(2)/unD->GetBinContent(2) <<endl;
    
    unD2 = new TH1F("unD2",
                     "unD2",
                     2,0,2);
    unD2->Fill(0.5);
    unD2->Fill(1.5);
    unD2->Fill(1.5);
    cout << "Histo 2:"<<endl;
    cout << "Bin 1: Contenu " << unD2->GetBinContent(1) <<" Erreur " << unD2->GetBinError(1) <<" Erreur relative " << unD2->GetBinError(1)/unD2->GetBinContent(1) <<endl;
    cout << "Bin 2: Contenu " << unD2->GetBinContent(2) <<" Erreur " << unD2->GetBinError(2) <<" Erreur relative " << unD2->GetBinError(2)/unD2->GetBinContent(2) <<endl;
    
    unD->Add(unD2);
    
    cout << "Histo sommé Add:"<<endl;
    cout << "Bin 1: Contenu " << unD->GetBinContent(1) <<" Erreur " << unD->GetBinError(1) <<" Erreur relative " << unD->GetBinError(1)/unD->GetBinContent(1) <<endl;
    cout << "Bin 2: Contenu " << unD->GetBinContent(2) <<" Erreur " << unD->GetBinError(2) <<" Erreur relative " << unD->GetBinError(2)/unD->GetBinContent(2) <<endl;
    
    unD->SetBinContent(1,unD->GetBinContent(1)*0.25);
    unD->SetBinContent(2,unD->GetBinContent(2)*0.25);
  //  unD->SetBinError(1,unD->GetBinError(1)*0.25);
    //unD->SetBinError(2,unD->GetBinError(2)*0.25);
    
    cout << "Histo scaled down:"<<endl;
    cout << "Bin 1: Contenu " << unD->GetBinContent(1) <<" Erreur " << unD->GetBinError(1) <<" Erreur relative " << unD->GetBinError(1)/unD->GetBinContent(1) <<endl;
    cout << "Bin 2: Contenu " << unD->GetBinContent(2) <<" Erreur " << unD->GetBinError(2) <<" Erreur relative " << unD->GetBinError(2)/unD->GetBinContent(2) <<endl;

}

void AutreTest(){
    cout << "On entre dans la fonction AutreTest, pouet a pour valeur " << pouet <<endl;
}

