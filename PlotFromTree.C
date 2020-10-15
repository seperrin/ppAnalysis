#include <TBrowser.h>
#include <TBufferFile.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TMath.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TFitResultPtr.h>
#include <TPaveText.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <string>

//#include "TTreeReader.h"
//#include "Event.h"

# include "TF1.h"
# include "TF2.h"
# include "TProfile.h"
# include "TVirtualFitter.h"
# include "TGraphErrors.h"
# include "TFitter.h"
# include "TMinuit.h"
# include "TRandom.h"
# include <iostream>
# include <string>
# include "TH1F.h"
# include "TH1D.h"
# include "TH2F.h"
# include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"


// FUNCTIONS

void PlotFromTree();
Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par);
Double_t FourierV2(Double_t *x,Double_t *par);
Double_t TwoCBE2E(Double_t *x,Double_t *par);
Double_t ExpBkg(Double_t *x,Double_t *par);
Double_t TwoCrystalBallExtended(Double_t *x,Double_t *par);
TFitResultPtr FittingAllInvMass(const char *histoname, TCanvas *canvas);
TFitResultPtr FittingAllInvMassBin(const char *histoname, TCanvas *canvas, int i);

void PlotFromTree(){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    
    bool doTracklets = kTRUE;
    
// ************************************
// Définitions de paramètres          *
// ************************************

    double ZvtxCut = 10;
    double DPhiCut = 0.01;
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 99999;
//    double MinMultCentral = 37;
//    double MaxMultPeriph = 23;
    double CentSPDTrackletsCentral = 0.75;
    double CentSPDTrackletsPeriph = 7;

    double LowJPsiMass = 3.0;
    double HighJPsiMass = 3.3;

    double MinInvMass = 2.0;
    double MaxInvMass = 4.5;
    const int NbinsInvMass = 10;
    double SizeBinInvMass = (MaxInvMass-MinInvMass)/NbinsInvMass;

    double MinDeltaPhi = -TMath::Pi()/2;
    double MaxDeltaPhi = 1.5*TMath::Pi();
    const int NbinsDeltaPhi = 12;
    double SizeBinDeltaPhi = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhi;
    
    double MinDeltaEta = 0;
    double MaxDeltaEta = 6;
    const int NbinsDeltaEta = 20;
    double SizeBinDeltaEta = (MaxDeltaEta-MinDeltaEta)/NbinsDeltaEta;
    
    double MinDeltaEtaTKL = 1.2;
    double MaxDeltaEtaTKL = 2.4;
    const int NbinsDeltaEtaTKL = 24;
    double SizeBinDeltaEtaTKL = (MaxDeltaEtaTKL-MinDeltaEtaTKL)/NbinsDeltaEtaTKL;
    
// *************************
// Initialiser les graphes *
// *************************
    
    TH1F* hnseg(NULL);
    TH1F* hnseg2(NULL);
    TH1F* hnseg3(NULL);
    TH1F* hnseg4(NULL);
    TH1F* hnseg5(NULL);
    TH1F* hnseg6(NULL);
    TH1F* hnseg7(NULL);
    TH2F* hnseg8(NULL);
    TH2F* hnseg8_masscut_allC(NULL);
    TH2F* hnseg8_masscut_central(NULL);
    TH2F* hnseg8_masscut_periph(NULL);
    TH2F* hnseg8TKL(NULL);
    TH2F* hnseg8TKL_central(NULL);
    TH2F* hnseg8TKL_periph(NULL);
    TH2F* hnseg8yTKL(NULL);
    TH2F* hnseg8yTKL_central(NULL);
    TH2F* hnseg8yTKL_periph(NULL);
    TH2F* yield_differenceTKL(NULL);
    TH1F* yield_centralTKL_proj_tampon(NULL);
    TH1F* yield_centralTKL_proj(NULL);
    TH1F* yield_periphTKL_proj_tampon(NULL);
    TH1F* yield_periphTKL_proj(NULL);
    TH1F* yield_differenceTKL_proj_tampon(NULL);
    TH1F* yield_differenceTKL_proj(NULL);
    
    TH2F* hnseg8y(NULL);
    TH2F* hnseg8y_masscut_allC(NULL);
    TH2F* hnseg8y_masscut_central(NULL);
    TH1F* yield_central_proj_tampon(NULL);
    TH1F* yield_central_proj(NULL);
    TH2F* hnseg8y_masscut_periph(NULL);
    TH1F* yield_periph_proj_tampon(NULL);
    TH1F* yield_periph_proj(NULL);
    TH2F* yield_difference(NULL);
    TH1F* yield_difference_proj_tampon(NULL);
    TH1F* yield_difference_proj(NULL);
    TH2F* hnseg9(NULL);
    TH2F* hnseg10(NULL);
    TH1F* hDPhi(NULL);
    TH2F* hPtWrtMassInv[3]{NULL};
    TH1 *hPtWrtMassInvSliced[3][6]= {NULL};
    
    TH1F* CentV0M(NULL);
    TH1F* CentTKL(NULL);
    TH1F* CentCL0(NULL);
    TH1F* CentCL1(NULL);
    TH1F* CentSPDTracklets(NULL);
    TH1F* CentSPDClusters(NULL);
    TH2F* CentV0Mwrttkl(NULL);
    TH2F* CentTKLwrttkl(NULL);
    TH2F* CentCL0wrttkl(NULL);
    TH2F* CentCL1wrttkl(NULL);
    TH2F* CentSPDTrackletswrttkl(NULL);
    TH2F* CentSPDClusterswrttkl(NULL);
    
//    Float_t fCentralityV0M;
//    Float_t fCentralityTKL;
//      Float_t fCentralityCL0;
//      Float_t fCentralityCL1;
//      Float_t fCentralitySPDTracklets;
//      Float_t fCentralitySPDClusters;
    
    TH2F* YieldsWrtDeltaPhiMassBin_Central[NbinsInvMass] = { NULL };
    TH2F* YieldsWrtDeltaPhiMassBin_Periph[NbinsInvMass] = { NULL };
    TH2F* YieldsWrtDeltaPhiMassBin_Difference[NbinsInvMass] = { NULL };
    TH1F* YieldsWrtDeltaPhiMassBin_CentralProj[NbinsInvMass] = { NULL };
    TH1F* YieldsWrtDeltaPhiMassBin_PeriphProj[NbinsInvMass] = { NULL };
    TH1F* YieldsWrtDeltaPhiMassBin_DifferenceProj[NbinsInvMass] = { NULL };
    TH1F* YieldsWrtDeltaPhiMassBin_CentralProj_Tampon(NULL);
    TH1F* YieldsWrtDeltaPhiMassBin_PeriphProj_Tampon(NULL);
    TH1F* YieldsWrtDeltaPhiMassBin_DifferenceProj_Tampon(NULL);
    
    TH1F* Yields_Central_1(NULL);
    TH1F* Yields_Periph_1(NULL);
    TH1F* Yields_Difference_1(NULL);
    
    TH1F* coefficients0(NULL);
    TH1F* coefficients1(NULL);
    TH1F* coefficients2(NULL);
    TH1F* baselines0(NULL);
    TH1F* c2b0(NULL);
    TH1F* V2JPsiTkl(NULL);
    
    Yields_Central_1 = new TH1F("Yields_Central_1",
                     "Yield of JPsi-tkl in Central collisions wrt Phi",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yields_Central_1->SetXTitle("Phi of the correlation (rad)");
    Yields_Central_1->SetYTitle("Yield");
    
    Yields_Periph_1 = new TH1F("Yields_Periph_1",
                     "Yield of JPsi-tkl in Periph collisions wrt Phi",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yields_Periph_1->SetXTitle("Phi of the correlation (rad)");
    Yields_Periph_1->SetYTitle("Yield");
    
    Yields_Difference_1 = new TH1F("Yields_Difference_1",
                     "Yield of JPsi-tkl in (Central-Periph) collisions wrt Phi",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yields_Difference_1->SetXTitle("Phi of the correlation (rad)");
    Yields_Difference_1->SetYTitle("Yield");
    
    coefficients0 = new TH1F("coefficients0",
                        "Coefficient Fourrier 0",
                        NbinsInvMass,MinInvMass,MaxInvMass);
       coefficients0->SetXTitle("Mass of dimuon (GeV/c^{2})");
       coefficients0->SetYTitle("Coefficient 0 (Fourier)");
    coefficients1 = new TH1F("coefficients1",
                           "Coefficient Fourrier 1",
                           NbinsInvMass,MinInvMass,MaxInvMass);
          coefficients1->SetXTitle("Mass of dimuon (GeV/c^{2})");
          coefficients1->SetYTitle("Coefficient 1 (Fourier)");
    coefficients2 = new TH1F("coefficients2",
                           "Coefficient Fourrier 2",
                           NbinsInvMass,MinInvMass,MaxInvMass);
          coefficients2->SetXTitle("Mass of dimuon (GeV/c^{2})");
          coefficients2->SetYTitle("Coefficient 2 (Fourier)");
    baselines0 = new TH1F("baselines0",
                     "Baseline 0",
                     NbinsInvMass,MinInvMass,MaxInvMass);
    baselines0->SetXTitle("Mass of dimuon (GeV/c^{2})");
    baselines0->SetYTitle("Baseline 0 (YieldPeriph deltaPhi=0)");
    c2b0 = new TH1F("c2b0",
                     "Coefficient 2 + Baseline 0",
                     NbinsInvMass,MinInvMass,MaxInvMass);
    c2b0->SetXTitle("Mass of dimuon (GeV/c^{2})");
    c2b0->SetYTitle("Coefficient 2 + Baseline 0");
    V2JPsiTkl = new TH1F("V2JPsiTkl",
                           "V2JPsiTkl wrt Mass",
                           NbinsInvMass,MinInvMass,MaxInvMass);
          V2JPsiTkl->SetXTitle("Mass of dimuon (GeV/c^{2})");
          V2JPsiTkl->SetYTitle("V2JPsiTkl");
    
    char hname[100];
    char hname2[100];
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Yield in Mass Bin %f GeV to %f GeV - Central",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    YieldsWrtDeltaPhiMassBin_Central[j] = new TH2F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Yield in Mass Bin %f GeV to %f GeV - Periph",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    YieldsWrtDeltaPhiMassBin_Periph[j] = new TH2F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Yield in Mass Bin %f GeV to %f GeV - Difference",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    YieldsWrtDeltaPhiMassBin_Difference[j] = new TH2F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    YieldsWrtDeltaPhiMassBin_CentralProj[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    YieldsWrtDeltaPhiMassBin_PeriphProj[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
    YieldsWrtDeltaPhiMassBin_DifferenceProj[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    }
    
    
    
    TH1F* InvMass_Central(NULL);
    TH1F* InvMass_Periph(NULL);
    
    
    
    
    TH1F* YieldWrtMass_PhiBin_Central[NbinsDeltaPhi] = { NULL };
    TH1F* YieldWrtMass_PhiBin_Periph[NbinsDeltaPhi] = { NULL };
    
    for (int j=0; j<NbinsDeltaPhi; j++){
        sprintf(hname,"Yield wrt Mass, Phi: %d pi/6 to %d pi/6 - Central",j-3,j-2);
        sprintf(hname2,"YieldWrtMassPhiBin_Central_%d",j);
        YieldWrtMass_PhiBin_Central[j] = new TH1F(hname2, hname,NbinsInvMass,MinInvMass,MaxInvMass);
    }
    
    for (int j=0; j<NbinsDeltaPhi; j++){
        sprintf(hname,"Yield wrt Mass, Phi: %d pi/6 to %d pi/6 - Periph",j-3,j-2);
        sprintf(hname2,"YieldWrtMassPhiBin_Periph_%d",j);
        YieldWrtMass_PhiBin_Periph[j] = new TH1F(hname2, hname,NbinsInvMass,MinInvMass,MaxInvMass);
    }
    
    
    
    
    hnseg = new TH1F("hnseg",
                     "Invariant mass of dimuon",
                     250,2,4.5);
    hnseg->SetXTitle("Mass of dimuon (GeV/c^{2})");
    hnseg->SetYTitle("Count");
    hDPhi = new TH1F("hDPhi",
                     "DeltaPhi of Tracklet",
                     10000,-TMath::Pi(),TMath::Pi());
    hDPhi->SetXTitle("DPhi tracklet");
    hDPhi->SetYTitle("Count");
    InvMass_Central = new TH1F("InvMass_Central",
                     "Invariant mass of dimuon - Central",
                     250,2,4.5);
    InvMass_Central->SetXTitle("Mass of dimuon (GeV/c^{2})");
    InvMass_Central->SetYTitle("Count");
    InvMass_Periph = new TH1F("InvMass_Periph",
                     "Invariant mass of dimuon - Periph",
                     250,2,4.5);
    InvMass_Periph->SetXTitle("Mass of dimuon (GeV/c^{2})");
    InvMass_Periph->SetYTitle("Count");
    hnseg2 = new TH1F("hnseg2",
                     "Pt of dimuon",
                     100,0,10);
    hnseg2->SetXTitle("Dimuon p_{t} (GeV/c)");
    hnseg2->SetYTitle("Count");
    hnseg3 = new TH1F("hnseg3",
                     "Eta of dimuon",
                     1000,-10,0);
    hnseg3->SetXTitle("Dimuon Eta");
    hnseg3->SetYTitle("Count");
    hnseg4 = new TH1F("hnseg4",
                      "Y of dimuon",
                      1000,-2.4,-4.2);
    hnseg4->SetXTitle("Dimuon Y");
    hnseg4->SetYTitle("Count");
    hnseg5 = new TH1F("hnseg5",
                      "Phi of dimuon",
                      1000,0, 2*TMath::Pi());
    hnseg5->SetXTitle("Dimuon Phi (rad)");
    hnseg5->SetYTitle("Count");
    hnseg6 = new TH1F("hnseg6",
                      "Correlation delta phi",
                      1000,MinDeltaPhi,MaxDeltaPhi);
    hnseg6->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg6->SetYTitle("Count");
    hnseg7 = new TH1F("hnseg7",
                      "Correlation delta eta",
                      1000,MinDeltaEta,MaxDeltaEta);
    hnseg7->SetXTitle("Correlation DeltaEta");
    hnseg7->SetYTitle("Count");
    hnseg8 = new TH2F("hnseg8",
                      "Correlation delta eta wrt delta phi",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8->SetYTitle("Correlation DeltaEta");
    
    hnseg8TKL = new TH2F("hnseg8TKL",
                      "Tracklets delta eta wrt delta phi",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8TKL->SetXTitle("Tracklets DeltaPhi (rad)");
    hnseg8TKL->SetYTitle("Tracklets DeltaEta");
    hnseg8TKL_central = new TH2F("hnseg8TKL_central",
                      "Tracklets delta eta wrt delta phi - Central",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8TKL_central->SetXTitle("Tracklets DeltaPhi (rad)");
    hnseg8TKL_central->SetYTitle("Tracklets DeltaEta");
    hnseg8TKL_periph = new TH2F("hnseg8TKL_periph",
                      "Tracklets delta eta wrt delta phi - Periph",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8TKL_periph->SetXTitle("Tracklets DeltaPhi (rad)");
    hnseg8TKL_periph->SetYTitle("Tracklets DeltaEta");
    hnseg8yTKL = new TH2F("hnseg8yTKL",
                      "Tracklets delta eta wrt delta phi, tracklet yields",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8yTKL->SetXTitle("Tracklets DeltaPhi (rad)");
    hnseg8yTKL->SetYTitle("Tracklets DeltaEta");
    hnseg8yTKL_central = new TH2F("hnseg8yTKL_central",
                                   "Tracklets delta eta wrt delta phi, mass cut, central, tracklet yields",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8yTKL_central->SetXTitle("Tracklets DeltaPhi (rad)");
    hnseg8yTKL_central->SetYTitle("Tracklets DeltaEta");
    hnseg8yTKL_periph = new TH2F("hnseg8yTKL_periph",
                                   "Tracklets delta eta wrt delta phi, mass cut, peripheral, tracklet yields",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8yTKL_periph->SetXTitle("Tracklets DeltaPhi (rad)");
    hnseg8yTKL_periph->SetYTitle("Tracklets DeltaEta");
    yield_differenceTKL = new TH2F("yield_differenceTKL", "Tracklets delta eta wrt delta phi, mass cut, central-peripheral, tracklet yields",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    yield_differenceTKL->SetXTitle("Tracklets DeltaPhi (rad)");
    yield_differenceTKL->SetYTitle("Tracklets DeltaEta");
    yield_centralTKL_proj = new TH1F("yield_centralTKL_proj", "Tracklets delta eta wrt delta phi, mass cut, central, tracklet yields, projection",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    yield_centralTKL_proj->SetXTitle("Tracklets DeltaPhi (rad)");
    yield_centralTKL_proj->SetYTitle("Tracklets Yield (Central)");
    yield_periphTKL_proj = new TH1F("yield_periphTKL_proj", "Tracklets delta eta wrt delta phi, mass cut, periph, tracklet yields, projection",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    yield_periphTKL_proj->SetXTitle("Tracklets DeltaPhi (rad)");
    yield_periphTKL_proj->SetYTitle("Tracklets Yield (Periph)");
    yield_differenceTKL_proj = new TH1F("yield_differenceTKL_proj", "Tracklets delta eta wrt delta phi, mass cut, central-peripheral, tracklet yields, projection",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
   yield_differenceTKL_proj->SetXTitle("Tracklets DeltaPhi (rad)");
   yield_differenceTKL_proj->SetYTitle("Tracklets Yield (Central - Peripheral)");
    
    hnseg8_masscut_allC = new TH2F("hnseg8_masscut_allC",
                                   "Correlation delta eta wrt delta phi, mass cut, all centralities",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8_masscut_allC->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8_masscut_allC->SetYTitle("Correlation DeltaEta");
    hnseg8_masscut_central = new TH2F("hnseg8_masscut_central",
                                   "Correlation delta eta wrt delta phi, mass cut, central",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8_masscut_central->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8_masscut_central->SetYTitle("Correlation DeltaEta");
    hnseg8_masscut_periph = new TH2F("hnseg8_masscut_periph",
                                   "Correlation delta eta wrt delta phi, mass cut, peripheral",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8_masscut_periph->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8_masscut_periph->SetYTitle("Correlation DeltaEta");
    hnseg8y = new TH2F("hnseg8y",
                      "Correlation delta eta wrt delta phi, tracklet yields",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8y->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8y->SetYTitle("Correlation DeltaEta");
    hnseg8y_masscut_allC = new TH2F("hnseg8y_masscut_allC",
                                   "Correlation delta eta wrt delta phi, mass cut, all centralities, tracklet yields",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8y_masscut_allC->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8y_masscut_allC->SetYTitle("Correlation DeltaEta");
    hnseg8y_masscut_central = new TH2F("hnseg8y_masscut_central",
                                   "Correlation delta eta wrt delta phi, mass cut, central, tracklet yields",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8y_masscut_central->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8y_masscut_central->SetYTitle("Correlation DeltaEta");
    yield_central_proj = new TH1F("yield_central_proj", "Correlation delta eta wrt delta phi, mass cut, central, tracklet yields, projection",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    yield_central_proj->SetXTitle("Correlation DeltaPhi (rad)");
    yield_central_proj->SetYTitle("Tracklet Yield (Central)");
    hnseg8y_masscut_periph = new TH2F("hnseg8y_masscut_periph",
                                   "Correlation delta eta wrt delta phi, mass cut, peripheral, tracklet yields",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8y_masscut_periph->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8y_masscut_periph->SetYTitle("Correlation DeltaEta");
    yield_periph_proj = new TH1F("yield_periph_proj", "Correlation delta eta wrt delta phi, mass cut, periph, tracklet yields, projection",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    yield_periph_proj->SetXTitle("Correlation DeltaPhi (rad)");
    yield_periph_proj->SetYTitle("Tracklet Yield (Periph)");
    yield_difference = new TH2F("yield_difference", "Correlation delta eta wrt delta phi, mass cut, central-peripheral, tracklet yields",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    yield_difference->SetXTitle("Correlation DeltaPhi (rad)");
    yield_difference->SetYTitle("Correlation DeltaEta");
    yield_difference_proj = new TH1F("yield_difference_proj", "Correlation delta eta wrt delta phi, mass cut, central-peripheral, tracklet yields, projection",NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    yield_difference_proj->SetXTitle("Correlation DeltaPhi (rad)");
    yield_difference_proj->SetYTitle("Tracklet Yield (Central - Peripheral)");
    
    hnseg9 = new TH2F("hnseg9",
                      "Tracklet eta wrt z_vertex",
                      100,-15,15,100,-2.5,2.5);
    hnseg9->SetXTitle("Z_{vertex} (cm)");
    hnseg9->SetYTitle("Tracklet Eta");
    hnseg10 = new TH2F("hnseg10",
                      "Number of tracklets wrt V0M percentile",
                      100,0,100,100,0,400);
    hnseg10->SetXTitle("V0M Percentile");
    hnseg10->SetYTitle("Tracklet count");
    Double_t binsPt[7] = {0,2,3,4,6,8,12};
    hPtWrtMassInv[0] = new TH2F("hPtWrtMassInv_0",
                      "Pt wrt Mass Inv - All centralities",
                      250,2,4.5,6,binsPt);
    hPtWrtMassInv[0]->SetXTitle("Mass Inv");
    hPtWrtMassInv[0]->SetYTitle("Pt");
    hPtWrtMassInv[1] = new TH2F("hPtWrtMassInv_1",
                      "Pt wrt Mass Inv - Central",
                      250,2,4.5,6,binsPt);
    hPtWrtMassInv[1]->SetXTitle("Mass Inv");
    hPtWrtMassInv[1]->SetYTitle("Pt");
    hPtWrtMassInv[2] = new TH2F("hPtWrtMassInv_2",
                      "Pt wrt Mass Inv - Periph",
                      250,2,4.5,6,binsPt);
    hPtWrtMassInv[2]->SetXTitle("Mass Inv");
    hPtWrtMassInv[2]->SetYTitle("Pt");
    
    
    CentV0M = new TH1F("CentV0M",
                     "Distribution of CentV0M",
                     100,0,100);
    CentV0M->SetXTitle("CentV0M");
    CentV0M->SetYTitle("Count");
    CentTKL = new TH1F("CentTKL",
                     "Distribution of CentTKL",
                     100,0,100);
    CentTKL->SetXTitle("CentTKL");
    CentTKL->SetYTitle("Count");
    CentCL0 = new TH1F("CentCL0",
                     "Distribution of CentCL0",
                     100,0,100);
    CentCL0->SetXTitle("CentCL0");
    CentCL0->SetYTitle("Count");
    CentCL1 = new TH1F("CentCL1",
                     "Distribution of CentCL1",
                     100,0,100);
    CentCL1->SetXTitle("CentCL1");
    CentCL1->SetYTitle("Count");
    CentSPDTracklets = new TH1F("CentSPDTracklets",
                     "Distribution of CentSPDTracklets",
                     100,0,100);
    CentSPDTracklets->SetXTitle("CentSPDTracklets");
    CentSPDTracklets->SetYTitle("Count");
    CentSPDClusters = new TH1F("CentSPDClusters",
                     "Distribution of CentSPDClusters",
                     100,0,100);
    CentSPDClusters->SetXTitle("CentSPDClusters");
    CentSPDClusters->SetYTitle("Count");
    
    CentV0Mwrttkl = new TH2F("CentV0Mwrttkl",
                     "Distribution of CentV0Mwrttkl",
                     100,0,100, 100, 0, 100);
    CentV0Mwrttkl->SetXTitle("CentV0Mwrttkl");
    CentV0Mwrttkl->SetYTitle("NTracklets");
    CentTKLwrttkl = new TH2F("CentTKLwrttkl",
                     "Distribution of CentTKLwrttkl",
                     100,0,100, 100, 0, 100);
    CentTKLwrttkl->SetXTitle("CentTKLwrttkl");
    CentTKLwrttkl->SetYTitle("NTracklets");
    CentCL0wrttkl = new TH2F("CentCL0wrttkl",
                     "Distribution of CentCL0wrttkl",
                     100,0,100, 100, 0, 100);
    CentCL0wrttkl->SetXTitle("CentCL0wrttkl");
    CentCL0wrttkl->SetYTitle("NTracklets");
    CentCL1wrttkl = new TH2F("CentCL1wrttkl",
                     "Distribution of CentCL1wrttkl",
                     100,0,100, 100, 0, 100);
    CentCL1wrttkl->SetXTitle("CentCL1wrttkl");
    CentCL1wrttkl->SetYTitle("NTracklets");
    CentSPDTrackletswrttkl = new TH2F("CentSPDTrackletswrttkl",
                     "Distribution of CentSPDTrackletswrttkl",
                     100,0,100, 100, 0, 100);
    CentSPDTrackletswrttkl->SetXTitle("CentSPDTrackletswrttkl");
    CentSPDTrackletswrttkl->SetYTitle("NTracklets");
    CentSPDClusterswrttkl = new TH2F("CentSPDClusterswrttkl",
                     "Distribution of CentSPDClusterswrttkl",
                     100,0,100, 100, 0, 100);
    CentSPDClusterswrttkl->SetXTitle("CentSPDClusterswrttkl");
    CentSPDClusterswrttkl->SetYTitle("NTracklets");
// *************************
// Analyse                 *
// *************************
    
    
        
    TFile fileIn("~/../../Volumes/Transcend2/ppAnalysis/Scripts/merge16hj_PS_CutsEvent.root");
    std::cout << "1" <<std::endl;
        TTree* theTree = NULL;
        fileIn.GetObject("MyMuonTree",theTree);
    std::cout << "2" <<std::endl;
        
    MyEventLight* fEvent = 0;
    TClonesArray *fCorrelations = 0;
    TClonesArray *fTracklets = 0;
    TClonesArray *fDimuons = 0;
    CorrelationLight *correl = 0; //= new CorrelationLight(); //object must be created before
    TrackletLight *trac = 0;
    TrackletLight *tracklet1 = 0;
    TrackletLight *tracklet2 = 0;
    DimuonLight *dimu = 0;
        //setting the branch address
    std::cout << "3" <<std::endl;
  //  theTree->SetBranchAddress("event", &fEvent);
    TBranch *branch = theTree->GetBranch("event");
//        auto bntrack = theTree->GetBranch("tracks");
        //auto branch  = theTree->GetBranch("mcparticles.fPx");
    std::cout << "3.5" <<std::endl;
        branch->SetAddress(&fEvent);
    std::cout << "4" <<std::endl;
    auto nevent = theTree->GetEntries();
    
    std::cout << "5" <<std::endl;
    fCorrelations = fEvent->fCorrelations;
    fTracklets = fEvent->fTracklets;
    fDimuons = fEvent->fDimuons;
//            auto nevent = bntrack->GetEntries();
    //cout << " theTree->GetEntries() " << theTree->GetEntries() << endl;
    
    int DimuCentralSeen = 0;
    int DimuPeriphSeen = 0;
    int DimuSeenMassCut = 0;
    int DimuSeenNoMassCut = 0;
    int EventWithPileupMult = 0;
    int RefTracklets = 0;
    int RefTrackletsCentral = 0;
    int RefTrackletsPeriph = 0;
    int cmul = 0;
    int barWidth = 50;
    int DimuCentralSeenMassBin[10];
    int DimuPeriphSeenMassBin[10];
    
    for(int i=0; i<NbinsInvMass; i++){
        DimuCentralSeenMassBin[i] = 0;
        DimuPeriphSeenMassBin[i] = 0;
    }
    
// ************************************
// Boucle sur les evenements, notés i *
// ************************************
    
  //  Double_t px[1000], py[1000];
        for (Int_t i=0;i<nevent;i++) {
            if(i%100000 == 0){
                    std::cout << "[";
                    int pos = barWidth * i / nevent;
                    for (int k = 0; k < barWidth; ++k) {
                        if (k < pos) std::cout << "=";
                        else if (k == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                    std::cout << "] " << int(i * 100 / nevent) << "%     " << i << "/" << nevent;
                    std::cout.flush();
                std::cout << std::endl;
            }
            if(doTracklets && (i%10!=0)){
                continue;
            }
            theTree->GetEvent(i);
     //       cout << "Looking at Event " << i << endl;
//            cout << "fPassPhysicsSelection : " << fEvent->fPassPhysicsSelection << endl;
//            cout << "fPassTriggerSelection : " << fEvent->fPassTriggerSelection << endl;
//            cout << "fTriggerCMUL7 : " << fEvent->fTriggerCMUL7 << endl;
//            cout << "fTriggerCINT7 : " << fEvent->fTriggerCINT7 << endl;
            //cout << "tracks->GetEntries() " << tracks->GetEntries() << endl;
          //  cout << "Computing the number of tracklets in Event " << i << endl;
            
            if(fEvent->fIsPileupFromSPDMultBins){
               // cout << "Event " << i << " has mult pileup" << endl;
                EventWithPileupMult++;
                continue;
            }
          //  cout << "Reading Event " << i <<endl;
            // Apply a cut on the TrackletEta at +-1
            int NumberCloseEtaTracklets = 0;
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                hnseg9->Fill(fEvent->fVertexZ, trac->fEta);
                if((TMath::Abs(trac->fEta) < 1)){
                    hDPhi->Fill(trac->fDPhi);
                    if((TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                    }
                }
                
            }
            
            CentV0M->Fill(fEvent->fCentralityV0M);
            CentTKL->Fill(fEvent->fCentralityTKL);
            CentCL0->Fill(fEvent->fCentralityCL0);
            CentCL1->Fill(fEvent->fCentralityCL1);
            CentSPDTracklets->Fill(fEvent->fCentralitySPDTracklets);
            CentSPDClusters->Fill(fEvent->fCentralitySPDClusters);
            
            CentV0Mwrttkl->Fill(fEvent->fCentralityV0M, NumberCloseEtaTracklets);
            CentTKLwrttkl->Fill(fEvent->fCentralityTKL, NumberCloseEtaTracklets);
            CentCL0wrttkl->Fill(fEvent->fCentralityCL0, NumberCloseEtaTracklets);
            CentCL1wrttkl->Fill(fEvent->fCentralityCL1, NumberCloseEtaTracklets);
            CentSPDTrackletswrttkl->Fill(fEvent->fCentralitySPDTracklets, NumberCloseEtaTracklets);
            CentSPDClusterswrttkl->Fill(fEvent->fCentralitySPDClusters, NumberCloseEtaTracklets);
            
         //   cout << "There are " << NumberCloseEtaTracklets << " selected tracklets (on Eta, DPhi) in Event " << i << endl;
//
//            if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut)){ //(fEvent->fNPileupVtx == 0) &&
//                Eventpassingcut++;
//            }
            
        //    cout << "Counting the dimuons in Event " << i << endl;
         //   cout << "In total, raw: " << fEvent->fNDimuons << " dimuons in Event " << i << endl;
            //Finding the number of dimuons of each type
            for (Int_t j=0; j<fEvent->fNDimuons; j++) {
             //   cout << "Dimuon " << j << " out of " << fEvent->fNDimuons << " in Event " << i << endl;
                dimu = (DimuonLight*)fDimuons->At(j);
               // hnseg->Fill(dimu->fInvMass);
                if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)){
                    //(fEvent->fNPileupVtx == 0) &&
                    //  (TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (dimu->fEta < HighDimuEtaCut ) && (dimu->fEta > LowDimuEtaCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)
                    // fEvent->fNPileupVtx == 0) && (TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)
               //     cout << "Dimuon " << j << " passed the cuts ( Zvtx, Eta, US, Pt) " << endl;
                 //   cout << "Mass " << dimu->fInvMass << endl;
                    DimuSeenNoMassCut++;
                   // cout << "From the beginning we have now seen in total " << DimuSeenNoMassCut << " dimuons" <<  endl;
                    hnseg->Fill(dimu->fInvMass);
                    hPtWrtMassInv[0]->Fill(dimu->fInvMass, dimu->fPt);
               //     cout << "Event " << i << " has " << NumberCloseEtaTracklets << " selected tracklets" << endl;
                    if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
          //          if(NumberCloseEtaTracklets >= MinMultCentral){
                //        cout << "The event is central" <<  endl;
                        InvMass_Central->Fill(dimu->fInvMass);
                        hPtWrtMassInv[1]->Fill(dimu->fInvMass, dimu->fPt);
                    }
                    if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
           //         if(NumberCloseEtaTracklets <= MaxMultPeriph){
                 //       cout << "The event is peripheral" <<  endl;
                        InvMass_Periph->Fill(dimu->fInvMass);
                        hPtWrtMassInv[2]->Fill(dimu->fInvMass, dimu->fPt);
                    }
                    if((dimu->fInvMass > LowJPsiMass) && (dimu->fInvMass < HighJPsiMass)){
             //           cout << "The dimuon " << j << " is in the range for the JPsi mass" <<  endl;
                        DimuSeenMassCut++;
                //        cout << "From the beginning we have now seen in total " << DimuSeenMassCut << " dimuons in the mass range of JPsi" <<  endl;
                        if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
                //       if(NumberCloseEtaTracklets >= MinMultCentral){
                     //      cout << "And, the dimuon " << j << " is in a central event" <<  endl;
                           DimuCentralSeen++;
                     //      cout << "From the beginning we have now seen in total " << DimuCentralSeen << " dimuons in the mass range of JPsi in central events" <<  endl;
                       }
                        if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
                //       if(NumberCloseEtaTracklets <= MaxMultPeriph){
                      //     cout << "And, the dimuon " << j << " is in a peripheral event" <<  endl;
                           DimuPeriphSeen++;
                      //     cout << "From the beginning we have now seen in total " << DimuPeriphSeen << " dimuons in the mass range of JPsi in peripheral events" <<  endl;
                       }
                    }
                    for(int bin = 0; bin<NbinsInvMass ; bin++){
                        if((dimu->fInvMass > MinInvMass + SizeBinInvMass*bin) && (dimu->fInvMass < MinInvMass + SizeBinInvMass*(bin+1))){
                            if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
                 //          if(NumberCloseEtaTracklets >= MinMultCentral){
                               DimuCentralSeenMassBin[bin]++;
                           }
                            if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
                   //        if(NumberCloseEtaTracklets <= MaxMultPeriph){
                               DimuPeriphSeenMassBin[bin]++;
                           }
                        }
                    }
                }
            }
            
// ***********************************
// Boucle sur les correlations, notés j    *
// ***********************************
      //      cout << "Counting the correlations in Event " << i << endl;
       //     cout << "In total, raw: " << fEvent->fNCorrelations << " correlations in Event " << i << endl;
            for (Int_t j=0; j<fEvent->fNCorrelations; j++) {
                correl = (CorrelationLight*)fCorrelations->At(j);
             //   cout << "Correlation " << j << " out of " << fEvent->fNCorrelations << " in Event " << i << endl;
                            //cout << truc->Pt() << endl;
                if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (correl->fDimuonY < HighDimuYCut ) && (correl->fDimuonY > LowDimuYCut) && (correl->fDimuonCharge == 0) && (correl->fDimuonPt > LowDimuPtCut) && (correl->fDimuonPt < HighDimuPtCut) && (TMath::Abs(correl->fTrackletEta) < 1) && (TMath::Abs(correl->fTrackletDPhi) < DPhiCut)){   //Cuts    (fEvent->fNPileupVtx == 0) &&
            //        cout << "Correlation " << j << " passed the cuts on dimuon and tracklets" << endl;
                    Float_t DeltaPhi = correl->fDeltaPhi;
             //       cout << "Correlation " << j << " has DeltaPhi of " << correl->fDeltaPhi << endl;
                    if(DeltaPhi < -TMath::Pi()/2){
                        DeltaPhi += 2* TMath::Pi();
                    }
                    if(DeltaPhi > 1.5*TMath::Pi()){
                        DeltaPhi -= 2* TMath::Pi();
                    }
            //        cout << "That is " << DeltaPhi << " modulo 2pi" <<endl;
                    hnseg2->Fill(correl->fDimuonPt);
                    hnseg3->Fill(correl->fDimuonEta);
                    hnseg4->Fill(correl->fDimuonY);
                    hnseg5->Fill(correl->fDimuonPhi);
                    hnseg6->Fill(DeltaPhi);
                    hnseg7->Fill(correl->fDeltaEta);
                    hnseg8->Fill(DeltaPhi, correl->fDeltaEta);
             //       cout << "The dimuon in correlation " << j << " has a mass of " << correl->fDimuonInvMass << endl;
                    if((correl->fDimuonInvMass > LowJPsiMass) && (correl->fDimuonInvMass < HighJPsiMass)){
            //            cout << "The dimuon mass is in the range of the JPsi mass" << endl;
                        hnseg8_masscut_allC->Fill(DeltaPhi, correl->fDeltaEta);
             //           cout << "Event " << i << " still has " << NumberCloseEtaTracklets << " selected tracklets" << endl;
                        if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
            //            if(NumberCloseEtaTracklets >= MinMultCentral){
             //               cout << "It is a central event" << endl;
                            hnseg8_masscut_central->Fill(DeltaPhi, correl->fDeltaEta);
                        }
                        if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
              //          if(NumberCloseEtaTracklets <= MaxMultPeriph){
             //               cout << "It is a peripheral event" << endl;
                            hnseg8_masscut_periph->Fill(DeltaPhi, correl->fDeltaEta);
                        }
                    }
                    
                    for(int bin = 0; bin<NbinsInvMass ; bin++){
                        if((correl->fDimuonInvMass > MinInvMass + SizeBinInvMass*bin) && (correl->fDimuonInvMass < MinInvMass + SizeBinInvMass*(bin+1))){
                            if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
              //             if(NumberCloseEtaTracklets >= MinMultCentral){
                               YieldsWrtDeltaPhiMassBin_Central[bin]->Fill(DeltaPhi, correl->fDeltaEta);
                           }
                            if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
               //            if(NumberCloseEtaTracklets <= MaxMultPeriph){
                               YieldsWrtDeltaPhiMassBin_Periph[bin]->Fill(DeltaPhi, correl->fDeltaEta);
                           }
                        }
                    }
                    
                    for(int bin = 0; bin<NbinsDeltaPhi ; bin++){
                        if((DeltaPhi >= MinDeltaPhi + bin*SizeBinDeltaPhi) && (DeltaPhi < MinDeltaPhi + (bin+1)*SizeBinDeltaPhi)){
                            if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
               //            if(NumberCloseEtaTracklets >= MinMultCentral){
                               YieldWrtMass_PhiBin_Central[bin]->Fill(correl->fDimuonInvMass);
                           }
                            if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
                 //          if(NumberCloseEtaTracklets <= MaxMultPeriph){
                               YieldWrtMass_PhiBin_Periph[bin]->Fill(correl->fDimuonInvMass);
                           }
                        }
                    }
                }
            }
            
            if(TMath::Abs(fEvent->fVertexZ) < ZvtxCut){
                hnseg10->Fill(fEvent->fCentralityV0M, fEvent->fNTracklets);
            }
            
            
            
            
            
            // TREATEMENT OF TRACKLET V2
            if(doTracklets){
         //   cout << "There are " << fEvent->fNTracklets << " tracklets" << endl;
                RefTracklets += fEvent->fNTracklets;
                if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
       //            if(NumberCloseEtaTracklets >= MinMultCentral){
        //               cout << "It is a central event" << endl;
                       RefTrackletsCentral += fEvent->fNTracklets;
                   }
               if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
     //          if(NumberCloseEtaTracklets <= MaxMultPeriph){
    //               cout << "It is a peripheral event" << endl;
                   RefTrackletsPeriph += fEvent->fNTracklets;
               }
                
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
        //        cout << "Tracklet 1 " << j << endl;
                tracklet1 = (TrackletLight*)fTracklets->At(j);
                for (Int_t k=j+1; k<fEvent->fNTracklets; k++){
             //       cout << "Tracklet 2 " << k << endl;
                    tracklet2 = (TrackletLight*)fTracklets->At(k);
                    if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(tracklet1->fDPhi) < DPhiCut) && (TMath::Abs(tracklet2->fDPhi) < DPhiCut)){   //Cuts
                            Float_t DeltaPhi = tracklet1->fPhi - tracklet2->fPhi;
                            if(DeltaPhi < -TMath::Pi()/2){
                                DeltaPhi += 2* TMath::Pi();
                            }
                            if(DeltaPhi > 1.5*TMath::Pi()){
                                DeltaPhi -= 2* TMath::Pi();
                            }
                            Float_t DeltaEta = TMath::Abs(tracklet1->fEta - tracklet2->fEta);
                            
                            hnseg8TKL->Fill(DeltaPhi, DeltaEta);
                        
                            if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
                                       //            if(NumberCloseEtaTracklets >= MinMultCentral){
                                        //               cout << "It is a central event" << endl;
                                                       hnseg8TKL_central->Fill(DeltaPhi, DeltaEta);
                                                   }
                           if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
                 //          if(NumberCloseEtaTracklets <= MaxMultPeriph){
                //               cout << "It is a peripheral event" << endl;
                               hnseg8TKL_periph->Fill(DeltaPhi, DeltaEta);
                           }
                        
                        }
                }
            }
        }
            
            
            
            
            
        }
    cout << "Reading of events finished" <<endl;
    
    cout << "Converting count histograms into yields, by dividing by the number of dimuons seen in each case" <<endl;
    for(int i=0; i<=hnseg8->GetNbinsX()+1; i++){
        for(int j=0; j<=hnseg8->GetNbinsY()+1; j++){
         //   cout << "Histogram hnseg8y looks at all dimuons regardless of mass, so division by " << DimuSeenNoMassCut << endl;
//            cout << "Histogram hnseg8y_masscut_allC looks at dimuons in mass range of JPsi, so division by " << DimuSeenMassCut << endl;
//            cout << "Histogram hnseg8y_masscut_central looks at central dimuons in mass range of JPsi, so division by " << DimuCentralSeen << endl;
//            cout << "Histogram hnseg8_masscut_periph looks at peripheral dimuons in mass range of JPsi, so division by " << DimuPeriphSeen << endl;
            hnseg8y->SetBinContent(i,j, hnseg8->GetBinContent(i,j)/DimuSeenNoMassCut);
            hnseg8y_masscut_allC->SetBinContent(i,j, hnseg8_masscut_allC->GetBinContent(i,j)/DimuSeenMassCut);
            hnseg8y_masscut_central->SetBinContent(i,j, hnseg8_masscut_central->GetBinContent(i,j)/DimuCentralSeen);
            hnseg8y_masscut_periph->SetBinContent(i,j, hnseg8_masscut_periph->GetBinContent(i,j)/DimuPeriphSeen);
            yield_difference->SetBinContent(i,j, hnseg8y_masscut_central->GetBinContent(i,j) - hnseg8y_masscut_periph->GetBinContent(i,j));
            
            for(int massbin=0; massbin<NbinsInvMass; massbin++){
                 YieldsWrtDeltaPhiMassBin_Central[massbin]->SetBinContent(i,j, YieldsWrtDeltaPhiMassBin_Central[massbin]->GetBinContent(i,j)/DimuCentralSeenMassBin[massbin]);
                YieldsWrtDeltaPhiMassBin_Periph[massbin]->SetBinContent(i,j, YieldsWrtDeltaPhiMassBin_Periph[massbin]->GetBinContent(i,j)/DimuPeriphSeenMassBin[massbin]);
                YieldsWrtDeltaPhiMassBin_Difference[massbin]->SetBinContent(i,j, YieldsWrtDeltaPhiMassBin_Central[massbin]->GetBinContent(i,j) - YieldsWrtDeltaPhiMassBin_Periph[massbin]->GetBinContent(i,j));
                
            }
            
        }
    }
    if(doTracklets){
    // For TklTkl Analysis
    for(int i=0; i<=hnseg8TKL->GetNbinsX()+1; i++){
            for(int j=0; j<=hnseg8TKL->GetNbinsY()+1; j++){
                hnseg8yTKL->SetBinContent(i,j, hnseg8TKL->GetBinContent(i,j)/RefTracklets);
                hnseg8yTKL_central->SetBinContent(i,j, hnseg8TKL_central->GetBinContent(i,j)/RefTrackletsCentral);
                hnseg8yTKL_periph->SetBinContent(i,j, hnseg8TKL_periph->GetBinContent(i,j)/RefTrackletsPeriph);
                yield_differenceTKL->SetBinContent(i,j, hnseg8yTKL_central->GetBinContent(i,j) - hnseg8yTKL_periph->GetBinContent(i,j));
            }
        }
    }
    
   // TH1F* InvMass_CentralRebinned(NULL);
    TH1F *InvMass_CentralRebinned = dynamic_cast<TH1F*>(InvMass_Central->Rebin(25,"InvMass_CentralRebinned"));
    TH1F *InvMass_PeriphRebinned = dynamic_cast<TH1F*>(InvMass_Periph->Rebin(25,"InvMass_PeriphRebinned"));
    double errYieldswrtMass_Central[12][11];
    double errYieldswrtMass_Periph[12][11];
    for(int phibin=0; phibin<NbinsDeltaPhi; phibin++){
        for(int i=0; i<=NbinsInvMass; i++){
            errYieldswrtMass_Central[phibin][i] = YieldWrtMass_PhiBin_Central[phibin]->GetBinError(i);
            errYieldswrtMass_Periph[phibin][i] = YieldWrtMass_PhiBin_Periph[phibin]->GetBinError(i);
        }
    }
    
   // for(int i=0; i<=NbinsInvMass; i++){
        for(int phibin=0; phibin<NbinsDeltaPhi; phibin++){
                        YieldWrtMass_PhiBin_Central[phibin]->Divide(InvMass_CentralRebinned);
                        YieldWrtMass_PhiBin_Periph[phibin]->Divide(InvMass_PeriphRebinned);
                     //  YieldWrtMass_PhiBin_Periph[phibin]->SetBinContent(i, YieldWrtMass_PhiBin_Periph[phibin]->GetBinContent(i)/DimuPeriphSeenMassBin[i]);
                   }
   // }
    for(int phibin=0; phibin<NbinsDeltaPhi; phibin++){
        for(int i=0; i<=NbinsInvMass; i++){
            YieldWrtMass_PhiBin_Central[phibin]->SetBinError(i+1, errYieldswrtMass_Central[phibin][i]/InvMass_CentralRebinned->GetBinContent(i));
            YieldWrtMass_PhiBin_Periph[phibin]->SetBinError(i+1, errYieldswrtMass_Periph[phibin][i]/InvMass_PeriphRebinned->GetBinContent(i));
        }
    }
    
    std::cout << "===== Dimuons counted in total =====" <<std::endl;
    std::cout<< "DimuCentralSeen (dimuons central in JPsi mass range) " << DimuCentralSeen <<std::endl;
    std::cout<< "DimuPeriphSeen (dimuons peripheral in JPsi mass range) " << DimuPeriphSeen <<std::endl;
    std::cout<< "DimuSeenMassCut (all dimuons in JPsi mass range) " << DimuSeenMassCut <<std::endl;
    std::cout<< "DimuSeenNoMassCut (all dimuons regardless of mass) " << DimuSeenNoMassCut <<std::endl;
    for(int i=0; i<NbinsInvMass; i++){
        std::cout<< "DimuCentralSeenMassBin[" << i << "] : " << DimuCentralSeenMassBin[i] << std::endl;
        std::cout<< "DimuPeriphSeenMassBin[" << i << "] : " << DimuPeriphSeenMassBin[i] << std::endl;
    }
    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    //}
    
    
    
// ***********************************
// Tracer les graphes sur les canevas *
// ***********************************


    TCanvas*c1=new TCanvas();
    //Canvas avec plein de petits plots généraux
    c1->Divide(4,2);
    c1->cd(1);
    hnseg->SetFillColor(18);
    hnseg->Draw();
    c1->cd(2);
    hnseg2->SetFillColor(18);
    hnseg2->Draw();
    c1->cd(3);
    hnseg3->SetFillColor(18);
    hnseg3->Draw();
    c1->cd(4);
    hnseg4->SetFillColor(18);
    hnseg4->Draw();
    c1->cd(5);
    hnseg5->SetFillColor(18);
    hnseg5->Draw();
    c1->cd(6);
    hnseg6->SetFillColor(18);
    hnseg6->Draw();
    c1->cd(7);
    hnseg7->SetFillColor(18);
    hnseg7->Draw();
    c1->cd(8);
    hnseg8->Draw("colz");
    
    TCanvas*cDPhi=new TCanvas();
    //Plot DPhi distribution
    hDPhi->SetFillColor(18);
    hDPhi->Draw();

    TCanvas*c2=new TCanvas();
    //Plot TrkEtaWrtZvtx
    hnseg9->Draw("colz");
    
    TCanvas*c3=new TCanvas();
    //Plut Ntrk wrt V0M
    hnseg10->Draw("colz");
    
    TCanvas*c4=new TCanvas();
    //Correlation count DeltaEta wrt DeltaPhi TH2
    c4->Divide(2,2);
    c4->cd(1);
    hnseg8->Draw("colz");
    c4->cd(2);
    hnseg8_masscut_allC->SetFillColor(18);
    hnseg8_masscut_allC->Draw("colz");
    c4->cd(3);
    hnseg8_masscut_central->SetFillColor(18);
    hnseg8_masscut_central->Draw("colz");
    c4->cd(4);
    hnseg8_masscut_periph->SetFillColor(18);
    hnseg8_masscut_periph->Draw("colz");
    
    TCanvas*c5=new TCanvas();
     //Correlation yield DeltaEta wrt DeltaPhi TH2
       c5->Divide(2,2);
       c5->cd(1);
       hnseg8y->SetFillColor(18);
       hnseg8y->Draw("colz");
       c5->cd(2);
       hnseg8y_masscut_allC->SetFillColor(18);
       hnseg8y_masscut_allC->Draw("colz");
       c5->cd(3);
       hnseg8y_masscut_central->SetFillColor(18);
       hnseg8y_masscut_central->Draw("colz");
       c5->cd(4);
       hnseg8y_masscut_periph->SetFillColor(18);
       hnseg8y_masscut_periph->Draw("colz");
    
     TCanvas* c6 = new TCanvas;
    //Correlation yield DeltaEta wrt DeltaPhi TH2 -> Projected for Periph and Central
        yield_central_proj_tampon = (TH1F*)(hnseg8y_masscut_central->ProjectionX("yield_central_proj",0,-1,"e"));
        for(int i=1; i<yield_central_proj->GetNbinsX()+1; i++){
            std::cout<< "Bin " << i << " de la projection Ctrl: " << yield_central_proj_tampon->GetBinContent(i) <<std::endl;
            yield_central_proj->SetBinContent(i, yield_central_proj_tampon->GetBinContent(i));
    //        double error=0;
    //       for(int j=1; j<yield_difference->GetNbinsY()+1; j++){
    //            error += pow((yield_difference_proj->GetBinContent(i) - yield_difference->GetBinContent(i,j)),2);
    //        }
    //        error /= yield_difference->GetNbinsY()-1;
    //        error = TMath::Sqrt(error);
            yield_central_proj->SetBinError(i, sqrt(yield_central_proj->GetBinContent(i))/sqrt(DimuCentralSeen));
        }
        yield_periph_proj_tampon = (TH1F*)(hnseg8y_masscut_periph->ProjectionX("yield_periph_proj",0,-1,"e"));
            for(int i=1; i<yield_periph_proj->GetNbinsX()+1; i++){
                std::cout<< "Bin " << i << " de la projection Periph: " << yield_periph_proj_tampon->GetBinContent(i) <<std::endl;
                yield_periph_proj->SetBinContent(i, yield_periph_proj_tampon->GetBinContent(i));
        //        double error=0;
        //       for(int j=1; j<yield_difference->GetNbinsY()+1; j++){
        //            error += pow((yield_difference_proj->GetBinContent(i) - yield_difference->GetBinContent(i,j)),2);
        //        }
        //        error /= yield_difference->GetNbinsY()-1;
        //        error = TMath::Sqrt(error);
                yield_periph_proj->SetBinError(i, sqrt(yield_periph_proj->GetBinContent(i))/sqrt(DimuPeriphSeen));
            }
    c6->Divide(2,2);
    c6->cd(1);
    hnseg8y_masscut_central->Draw("colz");
    c6->cd(2);
    hnseg8y_masscut_periph->Draw("colz");
    c6->cd(3);
    yield_central_proj->Draw("E");
    c6->cd(4);
    yield_periph_proj->Draw("E");
    c6->Draw();
    c6->Modified();
    c6->ForceUpdate();
    
    TCanvas* c7 = new TCanvas;
    //Correlation yield DeltaEta wrt DeltaPhi TH2 -> Projected Difference
    yield_difference_proj_tampon = (TH1F*)(yield_difference->ProjectionX("yield_difference_proj",0,-1,"e"));
    for(int i=1; i<yield_difference_proj->GetNbinsX()+1; i++){
        std::cout<< "Bin " << i << " de la projection Diffce: " << yield_difference_proj_tampon->GetBinContent(i) <<std::endl;
        yield_difference_proj->SetBinContent(i, yield_difference_proj_tampon->GetBinContent(i));
//        double error=0;
//       for(int j=1; j<yield_difference->GetNbinsY()+1; j++){
//            error += pow((yield_difference_proj->GetBinContent(i) - yield_difference->GetBinContent(i,j)),2);
//        }
//        error /= yield_difference->GetNbinsY()-1;
//        error = TMath::Sqrt(error);
        yield_difference_proj->SetBinError(i, TMath::Sqrt( pow(yield_central_proj->GetBinError(i),2)+pow(yield_periph_proj->GetBinError(i),2)) );
    }
    c7->Divide(1,2);
    c7->cd(1);
    yield_difference->Draw("colz");
    c7->cd(2);
    yield_difference_proj->Draw("E");
    c7->Draw();
    c7->Modified();
    c7->ForceUpdate();
    
    TCanvas* c8 = new TCanvas;
    //Correlation yield DeltaEta wrt DeltaPhi TH2 -> Central, Periph, Difference [MassBins]
    c8->Divide(NbinsInvMass,3);
    for(int i=1; i<=NbinsInvMass; i++){
        c8->cd(i);
        YieldsWrtDeltaPhiMassBin_Central[i-1]->Draw("colz");
        c8->cd(NbinsInvMass+i);
        YieldsWrtDeltaPhiMassBin_Periph[i-1]->Draw("colz");
        c8->cd(2*NbinsInvMass+i);
        YieldsWrtDeltaPhiMassBin_Difference[i-1]->Draw("colz");
    }
    
    TCanvas* c9 = new TCanvas;
    //Correlation yield DeltaEta wrt DeltaPhi TH2 -> PROJECTION Central, Periph, Difference [MassBins]
    c9->Divide(NbinsInvMass,3);
    for(int i=1; i<NbinsInvMass+1; i++){
        YieldsWrtDeltaPhiMassBin_CentralProj_Tampon = (TH1F*)(YieldsWrtDeltaPhiMassBin_Central[i-1]->ProjectionX("yieldbin_central_proj",0,-1,"e"));
            for(int j=1; j<NbinsDeltaPhi+1; j++){
                YieldsWrtDeltaPhiMassBin_CentralProj[i-1]->SetBinContent(j, YieldsWrtDeltaPhiMassBin_CentralProj_Tampon->GetBinContent(j));
                YieldsWrtDeltaPhiMassBin_CentralProj[i-1]->SetBinError(j, sqrt(YieldsWrtDeltaPhiMassBin_CentralProj[i-1]->GetBinContent(j))/sqrt(DimuCentralSeenMassBin[i-1]));
            }
        c9->cd(i);
        YieldsWrtDeltaPhiMassBin_CentralProj[i-1]->Draw("E");
        
        YieldsWrtDeltaPhiMassBin_PeriphProj_Tampon = (TH1F*)(YieldsWrtDeltaPhiMassBin_Periph[i-1]->ProjectionX("yieldbin_periph_proj",0,-1,"e"));
            for(int j=1; j<NbinsDeltaPhi+1; j++){
                YieldsWrtDeltaPhiMassBin_PeriphProj[i-1]->SetBinContent(j, YieldsWrtDeltaPhiMassBin_PeriphProj_Tampon->GetBinContent(j));
                YieldsWrtDeltaPhiMassBin_PeriphProj[i-1]->SetBinError(j, sqrt(YieldsWrtDeltaPhiMassBin_PeriphProj[i-1]->GetBinContent(j))/sqrt(DimuPeriphSeenMassBin[i-1]));
                if(j==4){
                    baselines0->Fill(MinInvMass + (-0.5+i)*SizeBinInvMass, YieldsWrtDeltaPhiMassBin_PeriphProj[i-1]->GetBinContent(j));
                    baselines0->SetBinError(i, YieldsWrtDeltaPhiMassBin_PeriphProj[i-1]->GetBinError(j));
                }
            }
        c9->cd(NbinsInvMass+i);
        YieldsWrtDeltaPhiMassBin_PeriphProj[i-1]->Draw("E");
        
            YieldsWrtDeltaPhiMassBin_DifferenceProj_Tampon = (TH1F*)(YieldsWrtDeltaPhiMassBin_Difference[i-1]->ProjectionX("yieldbin_difference_proj",0,-1,"e"));
                for(int j=1; j<NbinsDeltaPhi+1; j++){
                    if(YieldsWrtDeltaPhiMassBin_DifferenceProj[i-1]){
                 //       std::cout << "OK histo found, i = " << i <<std::endl;
                    }
                    double tampon = 0;
                    if(YieldsWrtDeltaPhiMassBin_DifferenceProj_Tampon){
                  //                         std::cout << "OK histo tampon found, i = " << i <<std::endl;
                    }
                    tampon = YieldsWrtDeltaPhiMassBin_DifferenceProj_Tampon->GetBinContent(j);
                    (YieldsWrtDeltaPhiMassBin_DifferenceProj[i-1])->SetBinContent(j, tampon);
                    (YieldsWrtDeltaPhiMassBin_DifferenceProj[i-1])->SetBinError(j, TMath::Sqrt( pow(YieldsWrtDeltaPhiMassBin_CentralProj[i-1]->GetBinError(j),2)+pow(YieldsWrtDeltaPhiMassBin_PeriphProj[i-1]->GetBinError(j),2)) );
                }
            c9->cd(2*NbinsInvMass+i);
            YieldsWrtDeltaPhiMassBin_DifferenceProj[i-1]->Draw("E");
    }
    
    TCanvas* c10 = new TCanvas;
    // Yields wrt mass -> Periph, Central [Phi bins]
    c10->Divide(6,4);
    for(int i=1; i<NbinsDeltaPhi+1; i++){
        c10->cd(i);
        YieldWrtMass_PhiBin_Central[i-1]->Draw("E");
        c10->cd(NbinsDeltaPhi+i);
        YieldWrtMass_PhiBin_Periph[i-1]->Draw("E");
    }
    
    TCanvas* c11a = new TCanvas;
    //Pt wrt Mass -> Central, Periph, All
    c11a->Divide(3,1);
    for(int i=0; i<3; i++){
        c11a->cd(i+1);
        hPtWrtMassInv[i]->Draw("colz");
    }
    TCanvas* c11b_0 = new TCanvas;
    //Mass -> All C [Pt bins]
    c11b_0->Divide(2,3);
    for (int i=0;i<6;i++){
        hPtWrtMassInvSliced[0][i] = hPtWrtMassInv[0]->ProjectionX(Form("bin%d_0",i+1),i+1,i+1);
    }
    TCanvas* c11b_1 = new TCanvas;
    //Mass -> Central [Pt bins]
    c11b_1->Divide(2,3);
    for (int i=0;i<6;i++){
        hPtWrtMassInvSliced[1][i] = hPtWrtMassInv[1]->ProjectionX(Form("bin%d_1",i+1),i+1,i+1);
    }
    TCanvas* c11b_2 = new TCanvas;
    //Mass -> Periph [Pt bins]
    c11b_2->Divide(2,3);
    for (int i=0;i<6;i++){
        hPtWrtMassInvSliced[2][i] = hPtWrtMassInv[2]->ProjectionX(Form("bin%d_2",i+1),i+1,i+1);
    }
    
    TCanvas*c12=new TCanvas();
    // Estimators
    c12->Divide(3,2);
    c12->cd(1);
    CentV0M->SetFillColor(18);
    CentV0M->Draw();
    c12->cd(2);
    CentTKL->SetFillColor(18);
    CentTKL->Draw();
    c12->cd(3);
    CentCL0->SetFillColor(18);
    CentCL0->Draw();
    c12->cd(4);
    CentCL1->SetFillColor(18);
    CentCL1->Draw();
    c12->cd(5);
    CentSPDTracklets->SetFillColor(18);
    CentSPDTracklets->Draw();
    c12->cd(6);
    CentSPDClusters->SetFillColor(18);
    CentSPDClusters->Draw();
    
    TCanvas*c13=new TCanvas();
    //Estimators
    c13->Divide(3,2);
    c13->cd(1);
    CentV0Mwrttkl->SetFillColor(18);
    CentV0Mwrttkl->Draw();
    c13->cd(2);
    CentTKLwrttkl->SetFillColor(18);
    CentTKLwrttkl->Draw();
    c13->cd(3);
    CentCL0wrttkl->SetFillColor(18);
    CentCL0wrttkl->Draw();
    c13->cd(4);
    CentCL1wrttkl->SetFillColor(18);
    CentCL1wrttkl->Draw();
    c13->cd(5);
    CentSPDTrackletswrttkl->SetFillColor(18);
    CentSPDTrackletswrttkl->Draw();
    c13->cd(6);
    CentSPDClusterswrttkl->SetFillColor(18);
    CentSPDClusterswrttkl->Draw();
    
    
    
    TCanvas*c14=new TCanvas();
    //Yield difference wrt Phi [Mass bins] fits
    c14->Divide(5,2);
    for(int i=1; i<=10; i++){
        c14->cd(i);
    // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
        TH1F *histo = YieldsWrtDeltaPhiMassBin_DifferenceProj[i-1];
      // create a TF1 with the range from 0 to 3 and 6 parameters
      TF1 *fitFcnV2_2 = new TF1("fitFcnV2_2",FourierV2,-TMath::Pi()/2,1.5*TMath::Pi(),3);
      fitFcnV2_2->SetNpx(500);
      fitFcnV2_2->SetLineWidth(4);
      fitFcnV2_2->SetLineColor(kMagenta);
      // first try without starting values for the parameters
      // This defaults to 1 for each param.
      // this results in an ok fit for the polynomial function
      // however the non-linear part (lorenzian) does not
      // respond well.
       Double_t params[3] = {1,1,1};
      fitFcnV2_2->SetParameters(params);
       TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
       TVirtualFitter::Fitter(histo)->SetPrecision();
    //  histo->Fit("fitFcn","0");
      // second try: set start values for some parameters
       
       fitFcnV2_2->SetParName(0,"a0");
       fitFcnV2_2->SetParName(1,"a1");
       fitFcnV2_2->SetParName(2,"a2");
       
      TFitResultPtr res = histo->Fit("fitFcnV2_2","SBMERI+","ep");
      // improve the pictu
    //   std::cout << "integral error: " << integralerror << std::endl;
        Double_t par[3];
        fitFcnV2_2->GetParameters(par);
        coefficients0->Fill(MinInvMass + (-0.5+i)*SizeBinInvMass, par[0]);
        coefficients1->Fill(MinInvMass + (-0.5+i)*SizeBinInvMass, par[1]);
        coefficients2->Fill(MinInvMass + (-0.5+i)*SizeBinInvMass, par[2]);
        coefficients0->SetBinError(i,fitFcnV2_2->GetParError(0));
        coefficients1->SetBinError(i,fitFcnV2_2->GetParError(1));
        coefficients2->SetBinError(i,fitFcnV2_2->GetParError(2));
      fitFcnV2_2->Draw("same");
      // draw the legend
      TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
      legend->SetTextFont(72);
      legend->SetTextSize(0.04);
        Char_t message[80];
        sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
        legend->AddEntry(fitFcnV2_2,message);
      legend->AddEntry(histo,"Data","lpe");
      legend->Draw();
    
    }
    
    TCanvas*c15=new TCanvas();
    // Fourier coefficients of Yields wrt Phi [Mass bins] -> Extraction method 2
    c15->Divide(2,2);
    c15->cd(1);
    coefficients0->Draw();
    c15->cd(2);
    baselines0->Draw();
    c15->cd(3);
    coefficients1->Draw();
    c15->cd(4);
    coefficients2->Draw();
            

    // Plot du V2 JPsi Tracklet
    TCanvas*c16=new TCanvas();
    //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
    c2b0->Add(coefficients0, baselines0);
    V2JPsiTkl->Divide(coefficients2, c2b0);
        
//
//        double V2_2[10];
//        double eV2_2[10];
//        double mass[10];
//        double emass[10];
//        for(int i=0; i<10; i++){
//            V2_2[i] = coefficients[i][2]/coefficients[i][0];
//            eV2_2[i] = abs(V2_2[i]*(sqrt(pow((ecoefficients[i][0]/coefficients[i][0]),2)+pow((ecoefficients[i][2]/coefficients[i][2]),2))));
//            mass[i] = MinInvMass + (0.5+i)*SizeBinInvMass;
//            emass[i] = 0;
//        }
//
//        V2JPsiTkl = new TGraphErrors(10,mass,V2_2,emass,eV2_2);
//        V2JPsiTkl->SetTitle("V2 Jpsi/tkl wrt Mass");
//        V2JPsiTkl->SetMarkerColor(4);
//        V2JPsiTkl->SetMarkerStyle(21);
    V2JPsiTkl->Draw();
        
        

    
    
    
    
    
    
    
    
    TFile *outputFile;
     outputFile = new TFile("/Volumes/SEBUSB/InvMass.root","RECREATE");
     hnseg->Write();
    InvMass_Central->Write();
    InvMass_Periph->Write();
    for(int j=0;j<3;j++){
        for (int i=0;i<6;i++){
            hPtWrtMassInvSliced[j][i]->Write();
        }
    }
    
//    for(int i=0; i<12; i++){
//        InvMass_PhiBin[i]->Write();
//    }
     outputFile->Close();
    
    
    
//    TCanvas *cinvmassphibin = new TCanvas("cinvmassphibin","Fitting PhiBins",10,10,700,500);
//    cinvmassphibin->Divide(4,3);
    TCanvas *cinvmass = new TCanvas("cinvmass","Fitting All Phi - Different Centralities",10,10,700,500);
    cinvmass->Divide(1,3);
//
    
    TFitResultPtr res;
    char histoname[50];
    sprintf(histoname,"hnseg");
    res = FittingAllInvMassBin(histoname, cinvmass, 0);
    
    double par[16];
    for(int i=0; i <16; i++){
        par[i] = res->Parameter(i);
        if(i>11){
            par[i] = 1;
        }
    }
    
    c16->cd();
    TF1 *fitV2_2 = new TF1("fitV2_2",FourierV2_WrtInvMass,2,4.5,16);
      fitV2_2->SetNpx(500);
      fitV2_2->SetLineWidth(4);
      fitV2_2->SetLineColor(kMagenta);
      fitV2_2->SetParameters(par);
       TVirtualFitter::Fitter(V2JPsiTkl)->SetMaxIterations(10000);
       TVirtualFitter::Fitter(V2JPsiTkl)->SetPrecision();
    //  histo->Fit("fitFcn","0");
      // second try: set start values for some parameters
    for(int i=0; i<=11; i++){
        fitV2_2->FixParameter(i,par[i]);
    }
       
       fitV2_2->SetParName(0,"N_{JPsi} Gaussian approx.");
       fitV2_2->SetParName(1,"M_{JPsi}");
       fitV2_2->SetParName(2,"Sigma_{JPsi}");
       fitV2_2->SetParName(3,"a_{1}");
       fitV2_2->SetParName(4,"n_{1}");
       fitV2_2->SetParName(5,"a_{2}");
       fitV2_2->SetParName(6,"n_{2}");
       fitV2_2->SetParName(7,"Norm_{Tails/Core}");
       fitV2_2->SetParName(8,"Norm_{TailLowM}");
       fitV2_2->SetParName(9,"Exp_{TailLowM}");
       fitV2_2->SetParName(10,"Norm_{TailHighM}");
       fitV2_2->SetParName(11,"Exp_{TailHighM}");
        fitV2_2->SetParName(12,"V2_2 JPsi");
        fitV2_2->SetParName(13,"V2_2 Bkg M2");
    fitV2_2->SetParName(14,"V2_2 Bkg M1");
    fitV2_2->SetParName(15,"V2_2 Nkg M0");
       
      res = V2JPsiTkl->Fit("fitV2_2","SBMERI+","ep");
      Double_t para[12];
      // writes the fit results into the par array
       gStyle->SetOptFit(1011);
      // draw the legend
      TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
      legend->SetTextFont(72);
      legend->SetTextSize(0.04);
    Char_t message[80];
    sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
     legend->AddEntry(fitV2_2,message);
      legend->AddEntry(V2JPsiTkl,"Data","lpe");
      legend->Draw();
    
    
    
    cinvmass->cd();
    sprintf(histoname,"InvMass_Central");
    res = FittingAllInvMassBin(histoname, cinvmass, 1);
    
    TCanvas* c10_fit = new TCanvas;
       c10_fit->Divide(6,4);
    for(int i=0; i <16; i++){
        par[i] = res->Parameter(i);
        if(i>11){
            par[i] = 1;
        }
    }
    
    for(int i=0; i<NbinsDeltaPhi; i++){
        c10_fit->cd(i+1);
        TF1 *fitY_1Central = new TF1("fitY_1Central",FourierV2_WrtInvMass,2,4.5,16);
          fitY_1Central->SetNpx(500);
          fitY_1Central->SetLineWidth(4);
          fitY_1Central->SetLineColor(kMagenta);
          fitY_1Central->SetParameters(par);
           TVirtualFitter::Fitter(YieldWrtMass_PhiBin_Central[i])->SetMaxIterations(10000);
           TVirtualFitter::Fitter(YieldWrtMass_PhiBin_Central[i])->SetPrecision();
        //  histo->Fit("fitFcn","0");
          // second try: set start values for some parameters
        for(int i=0; i<=11; i++){
            fitY_1Central->FixParameter(i,par[i]);
        }
           
           fitY_1Central->SetParName(0,"N_{JPsi} Gaussian approx.");
           fitY_1Central->SetParName(1,"M_{JPsi}");
           fitY_1Central->SetParName(2,"Sigma_{JPsi}");
           fitY_1Central->SetParName(3,"a_{1}");
           fitY_1Central->SetParName(4,"n_{1}");
           fitY_1Central->SetParName(5,"a_{2}");
           fitY_1Central->SetParName(6,"n_{2}");
           fitY_1Central->SetParName(7,"Norm_{Tails/Core}");
           fitY_1Central->SetParName(8,"Norm_{TailLowM}");
           fitY_1Central->SetParName(9,"Exp_{TailLowM}");
           fitY_1Central->SetParName(10,"Norm_{TailHighM}");
           fitY_1Central->SetParName(11,"Exp_{TailHighM}");
            fitY_1Central->SetParName(12,"Y_1 JPsi");
            fitY_1Central->SetParName(13,"Y_1 Bkg M2");
        fitY_1Central->SetParName(14,"Y_1 Bkg M1");
        fitY_1Central->SetParName(15,"Y_1 Nkg M0");
           
          res = YieldWrtMass_PhiBin_Central[i]->Fit("fitY_1Central","SBMERQ+","ep");
          Double_t par[16];
          fitY_1Central->GetParameters(par);
          Yields_Central_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, par[12]);
          Yields_Central_1->SetBinError(i+1,fitY_1Central->GetParError(12));
          // writes the fit results into the par array
           gStyle->SetOptFit(1011);
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetTextFont(72);
          legend->SetTextSize(0.04);
        Char_t message[80];
        sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
         legend->AddEntry(fitY_1Central,message);
          legend->AddEntry(YieldWrtMass_PhiBin_Central[i],"Data","lpe");
          legend->Draw();
    }
    
    cinvmass->cd();
    sprintf(histoname,"InvMass_Periph");
    res = FittingAllInvMassBin(histoname, cinvmass, 2);
    
       for(int i=0; i <16; i++){
           par[i] = res->Parameter(i);
           if(i>11){
               par[i] = 1;
           }
       }
    
    double baseline = 1;
    double errbaseline = 1;
       
       for(int i=0; i<NbinsDeltaPhi; i++){
           c10_fit->cd(NbinsDeltaPhi+i+1);
           TF1 *fitY_1Periph = new TF1("fitY_1Periph",FourierV2_WrtInvMass,2,4.5,16);
             fitY_1Periph->SetNpx(500);
             fitY_1Periph->SetLineWidth(4);
             fitY_1Periph->SetLineColor(kMagenta);
             fitY_1Periph->SetParameters(par);
              TVirtualFitter::Fitter(YieldWrtMass_PhiBin_Periph[i])->SetMaxIterations(10000);
              TVirtualFitter::Fitter(YieldWrtMass_PhiBin_Periph[i])->SetPrecision();
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters
           for(int j=0; j<=11; j++){
               fitY_1Periph->FixParameter(j,par[j]);
           }
              
              fitY_1Periph->SetParName(0,"N_{JPsi} Gaussian approx.");
              fitY_1Periph->SetParName(1,"M_{JPsi}");
              fitY_1Periph->SetParName(2,"Sigma_{JPsi}");
              fitY_1Periph->SetParName(3,"a_{1}");
              fitY_1Periph->SetParName(4,"n_{1}");
              fitY_1Periph->SetParName(5,"a_{2}");
              fitY_1Periph->SetParName(6,"n_{2}");
              fitY_1Periph->SetParName(7,"Norm_{Tails/Core}");
              fitY_1Periph->SetParName(8,"Norm_{TailLowM}");
              fitY_1Periph->SetParName(9,"Exp_{TailLowM}");
              fitY_1Periph->SetParName(10,"Norm_{TailHighM}");
              fitY_1Periph->SetParName(11,"Exp_{TailHighM}");
               fitY_1Periph->SetParName(12,"Y_1 JPsi");
               fitY_1Periph->SetParName(13,"Y_1 Bkg M2");
           fitY_1Periph->SetParName(14,"Y_1 Bkg M1");
           fitY_1Periph->SetParName(15,"Y_1 Nkg M0");
              
             res = YieldWrtMass_PhiBin_Periph[i]->Fit("fitY_1Periph","SBMERQ+","ep");
             Double_t par[16];
           if(i==3){
               baseline = par[12];
               errbaseline = fitY_1Periph->GetParError(12);
           }
             fitY_1Periph->GetParameters(par);
             Yields_Periph_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, par[12]);
             Yields_Periph_1->SetBinError(i+1,fitY_1Periph->GetParError(12));
             // writes the fit results into the par array
              gStyle->SetOptFit(1011);
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
           Char_t message[80];
           sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
            legend->AddEntry(fitY_1Periph,message);
             legend->AddEntry(YieldWrtMass_PhiBin_Periph[i],"Data","lpe");
             legend->Draw();
       }
    
    {
        TCanvas* c17 = new TCanvas;
        c17->Divide(1,3);
        //Creer canvas pour imprimer les plots Periph yield et Central yield wrt phi et leur difference
        c17->cd(1);
        Yields_Central_1->Draw();
        c17->cd(2);
        Yields_Periph_1->Draw();
        c17->cd(3);
        Yields_Difference_1->Add(Yields_Central_1,Yields_Periph_1,1,-1);
        Yields_Difference_1->Draw();
    
             TF1 *fitFcnV2_2 = new TF1("fitFcnV2_2",FourierV2,-TMath::Pi()/2,1.5*TMath::Pi(),3);
             fitFcnV2_2->SetNpx(500);
             fitFcnV2_2->SetLineWidth(4);
             fitFcnV2_2->SetLineColor(kMagenta);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[3] = {1,1,1};
             fitFcnV2_2->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Difference_1)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Difference_1)->SetPrecision();
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters
              
              fitFcnV2_2->SetParName(0,"a0");
              fitFcnV2_2->SetParName(1,"a1");
              fitFcnV2_2->SetParName(2,"a2");
              
             TFitResultPtr res = Yields_Difference_1->Fit("fitFcnV2_2","SBMERI+","ep");
            Double_t par[3];
            fitFcnV2_2->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_2->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
               Char_t message[80];
               sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
               legend->AddEntry(fitFcnV2_2,message);
        sprintf(message,"V2_1 JPsi-Tkl: %.4f / (%.4f + %.4f) = %.4f +- %.4f",par[2],par[0], baseline,par[2]/(par[0] + baseline),(par[2]/(par[0] + baseline)*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(fitFcnV2_2->GetParError(0)/par[0],2)+pow(errbaseline/baseline,2))));
        legend->AddEntry(fitFcnV2_2,message);
             legend->AddEntry(Yields_Difference_1,"Data","lpe");
             legend->Draw();
           
    }
    
    
    double baselineTKL = 1;
    double errbaselineTKL = 1;
    
    if(doTracklets){
        TCanvas* c18 = new TCanvas;
        c18->cd();
        hnseg8TKL->Draw("colz");
        
        TCanvas*c4TKL=new TCanvas();
        //Tracklets count DeltaEta wrt DeltaPhi TH2
        c4TKL->Divide(2,2);
        c4TKL->cd(1);
        hnseg8TKL->Draw("colz");
        c4TKL->cd(3);
        hnseg8TKL_central->SetFillColor(18);
        hnseg8TKL_central->Draw("colz");
        c4TKL->cd(4);
        hnseg8TKL_periph->SetFillColor(18);
        hnseg8TKL_periph->Draw("colz");
        
        TCanvas*c5TKL=new TCanvas();
         //Correlation yield DeltaEta wrt DeltaPhi TH2
           c5TKL->Divide(2,2);
           c5TKL->cd(1);
           hnseg8yTKL->SetFillColor(18);
           hnseg8yTKL->Draw("colz");
           c5TKL->cd(3);
           hnseg8yTKL_central->SetFillColor(18);
           hnseg8yTKL_central->Draw("colz");
           c5TKL->cd(4);
           hnseg8yTKL_periph->SetFillColor(18);
           hnseg8yTKL_periph->Draw("colz");
        
        TCanvas* c6TKL = new TCanvas;
            //Tracklets yield DeltaEta wrt DeltaPhi TH2 -> Projected for Periph and Central
                yield_centralTKL_proj_tampon = (TH1F*)(hnseg8yTKL_central->ProjectionX("yield_centralTKL_proj",0,-1,"e"));
                for(int i=1; i<yield_centralTKL_proj->GetNbinsX()+1; i++){
                    yield_centralTKL_proj->SetBinContent(i, yield_centralTKL_proj_tampon->GetBinContent(i));
                    yield_centralTKL_proj->SetBinError(i, sqrt(yield_centralTKL_proj->GetBinContent(i))/sqrt(RefTrackletsCentral));
                }
                yield_periphTKL_proj_tampon = (TH1F*)(hnseg8yTKL_periph->ProjectionX("yield_periphTKL_proj",0,-1,"e"));
                    for(int i=1; i<yield_periphTKL_proj->GetNbinsX()+1; i++){
                        yield_periphTKL_proj->SetBinContent(i, yield_periphTKL_proj_tampon->GetBinContent(i));
                        yield_periphTKL_proj->SetBinError(i, sqrt(yield_periphTKL_proj->GetBinContent(i))/sqrt(RefTrackletsPeriph));
                        if(i==4){
                            baselineTKL = yield_periphTKL_proj->GetBinContent(i);
                            errbaselineTKL = yield_periphTKL_proj->GetBinError(i);
                        }
                    }
            c6TKL->Divide(2,2);
            c6TKL->cd(1);
            hnseg8yTKL_central->Draw("colz");
            c6TKL->cd(2);
            hnseg8yTKL_periph->Draw("colz");
            c6TKL->cd(3);
            yield_centralTKL_proj->Draw("E");
            c6TKL->cd(4);
            yield_periphTKL_proj->Draw("E");
            c6TKL->Draw();
            c6TKL->Modified();
            c6TKL->ForceUpdate();
            
            TCanvas* c7TKL = new TCanvas;
            //Tracklets yield DeltaEta wrt DeltaPhi TH2 -> Projected Difference
            yield_differenceTKL_proj_tampon = (TH1F*)(yield_differenceTKL->ProjectionX("yield_differenceTKL_proj",0,-1,"e"));
            for(int i=1; i<yield_differenceTKL_proj->GetNbinsX()+1; i++){
                yield_differenceTKL_proj->SetBinContent(i, yield_differenceTKL_proj_tampon->GetBinContent(i));
                yield_differenceTKL_proj->SetBinError(i, TMath::Sqrt( pow(yield_centralTKL_proj->GetBinError(i),2)+pow(yield_periphTKL_proj->GetBinError(i),2)) );
            }
            c7TKL->Divide(1,2);
            c7TKL->cd(1);
            yield_differenceTKL->Draw("colz");
            c7TKL->cd(2);
            yield_differenceTKL_proj->Draw("E");
            c7TKL->Draw();
            c7TKL->Modified();
            c7TKL->ForceUpdate();
        
        TCanvas*c14TKL=new TCanvas();
        //Tracklets Yield difference wrt Phi fit
        c14TKL->cd();
        // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
            TH1F *histo = yield_differenceTKL_proj;
          // create a TF1 with the range from 0 to 3 and 6 parameters
          TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV2,-TMath::Pi()/2,1.5*TMath::Pi(),3);
          fitFcnV2TKL->SetNpx(500);
          fitFcnV2TKL->SetLineWidth(4);
          fitFcnV2TKL->SetLineColor(kMagenta);
          // first try without starting values for the parameters
          // This defaults to 1 for each param.
          // this results in an ok fit for the polynomial function
          // however the non-linear part (lorenzian) does not
          // respond well.
           Double_t params[3] = {1,1,1};
          fitFcnV2TKL->SetParameters(params);
           TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
           TVirtualFitter::Fitter(histo)->SetPrecision();
        //  histo->Fit("fitFcn","0");
          // second try: set start values for some parameters
           
           fitFcnV2TKL->SetParName(0,"a0");
           fitFcnV2TKL->SetParName(1,"a1");
           fitFcnV2TKL->SetParName(2,"a2");
           
          TFitResultPtr res = histo->Fit("fitFcnV2TKL","SBMERI+","ep");
          // improve the pictu
        //   std::cout << "integral error: " << integralerror << std::endl;
            Double_t par[3];
            fitFcnV2TKL->GetParameters(par);
          fitFcnV2TKL->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetTextFont(72);
          legend->SetTextSize(0.04);
            Char_t message[80];
            sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
            legend->AddEntry(fitFcnV2TKL,message);
        sprintf(message,"V2 Tkl-Tkl: %.4f / (%.4f + %.4f) = %.4f +- %.4f",par[2],par[0], baselineTKL,par[2]/(par[0] + baselineTKL),(par[2]/(par[0] + baselineTKL)*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)+pow(errbaselineTKL/baselineTKL,2))));
        legend->AddEntry(fitFcnV2TKL,message);
          legend->AddEntry(histo,"Data","lpe");
          legend->Draw();
        
    }
    
    
    
        for (int i=0;i<6;i++){
               sprintf(histoname,Form("bin%d_0",i+1));
               res = FittingAllInvMassBin(histoname, c11b_0, i);
           }
    for (int i=0;i<6;i++){
        sprintf(histoname,Form("bin%d_1",i+1));
        res = FittingAllInvMassBin(histoname, c11b_1, i);
    }
    for (int i=0;i<6;i++){
        sprintf(histoname,Form("bin%d_2",i+1));
        res = FittingAllInvMassBin(histoname, c11b_2, i);
    }
//    for(int i=0; i<12; i++){
//  //  histoname = "Inv Mass, Phi: " + std::to_string(0) + " pi/6 to " + std::to_string(1) + " pi/6";
//        sprintf(histoname,"InvMassPhiBin_%d",i);
//        FittingAllInvMassPhiBin(histoname, cinvmassphibin, i);
//    }
    std::cout << "===== Dimuons counted =====" <<std::endl;
    std::cout<< "DimuCentralSeen " << DimuCentralSeen <<std::endl;
    std::cout<< "DimuPeriphSeen " << DimuPeriphSeen <<std::endl;
    std::cout<< "DimuSeenMassCut " << DimuSeenMassCut <<std::endl;
    std::cout<< "DimuSeenNoMassCut " << DimuSeenNoMassCut <<std::endl;
    for(int i=0; i<NbinsInvMass; i++){
        std::cout<< "DimuCentralSeenMassBin[" << i << "] : " << DimuCentralSeenMassBin[i] << std::endl;
        std::cout<< "DimuPeriphSeenMassBin[" << i << "] : " << DimuPeriphSeenMassBin[i] << std::endl;
    }
    std::cout << "EventWithPileupMult  " << EventWithPileupMult <<std::endl;
    std::cout << "cmul  " << cmul <<std::endl;
}
    









// FITTING INVARIANT MASS METHODS

Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par)
// Par 0->7: signal, 8->11 Bkg, 12: v2 JPsi, 13->15: V2 bkg
{ return (TwoCrystalBallExtended(x,par)*par[12] + ExpBkg(x,&par[8])*(par[13]*x[0]*x[0] + par[14]*x[0] + par[15]))/(TwoCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t FourierV2(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);}

Double_t TwoCBE2E(Double_t *x,Double_t *par)

{ return TwoCrystalBallExtended(x,par) + ExpBkg(x,&par[8]);}

Double_t ExpBkg(Double_t *x,Double_t *par)

{ return par[0]*(exp(x[0]*par[1]*(-1))) + par[2]*(exp(x[0]*par[3]*(-1))); }

Double_t TwoCrystalBallExtended(Double_t *x,Double_t *par)
{
    
    Double_t mJpsi =  3.096916;
    Double_t mPsip =  3.686108;
  Double_t ratMass = 1.01;
  Double_t ratSigma = 1.02;
    Double_t sum = 0;
    
    Double_t t = (x[0]-par[1])/par[2];
    if (par[3] < 0) t = -t;
    
    Double_t absAlpha = fabs((Double_t)par[3]);
    Double_t absAlpha2 = fabs((Double_t)par[5]);
    
    
    if (t < -absAlpha) //left tail
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        
        sum += (par[0]/(par[2]*sqrt(TMath::Pi()*2)))*(a/TMath::Power(b - t, par[4]));
    } else if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += (par[0]/(par[2]*sqrt(TMath::Pi()*2)))*(exp(-0.5*t*t));
    } else if (t >= absAlpha2) //right
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        
        sum += (par[0]/(par[2]*sqrt(TMath::Pi()*2)))*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
//  t = (x[0]-(par[1]*mPsip/mJpsi))/(par[2]*mPsip/mJpsi);
  t = (x[0]-(par[1]*mPsip/mJpsi))/(par[2]*ratSigma);
    if (par[3] < 0) t = -t;
    
    absAlpha = fabs((Double_t)par[3]);
    absAlpha2 = fabs((Double_t)par[5]);
    
    if (t < -absAlpha) //left tail
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        
        sum += (par[7]/(par[2]*ratSigma*sqrt(2*TMath::Pi())))*(a/TMath::Power(b - t, par[4]));
    }
    else if (t >= -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += (par[7]/(par[2]*ratSigma*sqrt(2*TMath::Pi())))*(exp(-0.5*t*t));
    }
    else if (t >= absAlpha2) //right tail
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        
        sum += (par[7]/(par[2]*ratSigma*sqrt(2*TMath::Pi())))*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
    return sum ;
}
    
    TFitResultPtr FittingAllInvMass(const char *histoname, TCanvas *cinvmass){
        TFitResultPtr res = FittingAllInvMassBin(histoname, cinvmass, 0);
        return res;
    }

TFitResultPtr FittingAllInvMassBin(const char *histoname, TCanvas *cinvmass, int i){
 //Bevington Exercise by Peter Malzacher, modified by Rene Brun, TO BE MODIFIED FOR OUR USE
    // PARAMETERS
    /*
     ------ Signal CB2
     0: Normalisation Gauss
     1: Mean x
     2: Sig Gauss
     3: Alpha 1
     4: Power 1
     5: Alpha 2
     6: Power 2
     7: Normalisation Tails-Core
     ------ Background 2Exp
     8: Normalisation Before
     9: Exponent Before
     10: Norlmalisation After
     11: Exponent After
     */
    TFile *file0 = new TFile("/Volumes/SEBUSB/InvMass.root");
    
    cinvmass->cd(i+1);
   cinvmass->SetFillColor(33);
   cinvmass->SetFrameFillColor(41);
   cinvmass->SetGrid();
   TH1F *histo = (TH1F*)file0->Get(histoname);
   // create a TF1 with the range from 0 to 3 and 6 parameters
   TF1 *fitFcn = new TF1("fitFcn",TwoCBE2E,2,4.5,12);
   fitFcn->SetNpx(500);
   fitFcn->SetLineWidth(4);
   fitFcn->SetLineColor(kMagenta);
   // first try without starting values for the parameters
   // This defaults to 1 for each param.
   // this results in an ok fit for the polynomial function
   // however the non-linear part (lorenzian) does not
   // respond well.
    Double_t params[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
   fitFcn->SetParameters(params);
    TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
    TVirtualFitter::Fitter(histo)->SetPrecision();
 //  histo->Fit("fitFcn","0");
   // second try: set start values for some parameters
   fitFcn->FixParameter(1,3.096916); // Mean x core
    fitFcn->SetParameter(2,0.05);
    fitFcn->SetParameter(3,0.9);
    fitFcn->SetParameter(4,10);
    fitFcn->SetParameter(5,2);
    fitFcn->SetParameter(6,15);
    
    fitFcn->SetParLimits(0,0.1,1000000);
  //  fitFcn->SetParLimits(1,15);
    fitFcn->SetParLimits(2,0.01,0.2);
    fitFcn->SetParLimits(3,0.1,5);
    fitFcn->SetParLimits(4,1,40);
    fitFcn->SetParLimits(5,0.1,20);
    fitFcn->SetParLimits(6,0.01,50);
    fitFcn->SetParLimits(7,0.1,1000000);
    fitFcn->SetParLimits(8,0.1,1000000);
    fitFcn->SetParLimits(9,0.1,20);
    fitFcn->SetParLimits(10,0.1,100000000);
    fitFcn->SetParLimits(11,0.1,200);
    
    fitFcn->SetParName(0,"N_{JPsi} Gaussian approx.");
    fitFcn->SetParName(1,"M_{JPsi}");
    fitFcn->SetParName(2,"Sigma_{JPsi}");
    fitFcn->SetParName(3,"a_{1}");
    fitFcn->SetParName(4,"n_{1}");
    fitFcn->SetParName(5,"a_{2}");
    fitFcn->SetParName(6,"n_{2}");
    fitFcn->SetParName(7,"Norm_{Tails/Core}");
    fitFcn->SetParName(8,"Norm_{TailLowM}");
    fitFcn->SetParName(9,"Exp_{TailLowM}");
    fitFcn->SetParName(10,"Norm_{TailHighM}");
    fitFcn->SetParName(11,"Exp_{TailHighM}");
    
   TFitResultPtr res = histo->Fit("fitFcn","SBMERQ+","ep");
   // improve the picture:
   TF1 *backFcn = new TF1("backFcn",ExpBkg,2,4.5,4);
   backFcn->SetLineColor(kRed);
   TF1 *signalFcn = new TF1("signalFcn",TwoCrystalBallExtended,2,4.5,8);
   TPaveText *pave = new TPaveText(0.15,0.5,0.3,0.65,"brNDC");
   signalFcn->SetLineColor(kBlue);
   signalFcn->SetNpx(500);
   Double_t par[12];
   // writes the fit results into the par array
    gStyle->SetOptFit(1011);
   fitFcn->GetParameters(par);
   signalFcn->SetParameters(par);
   Double_t integral = (signalFcn->Integral(2,3.45))/0.01;
   Double_t integralerror = (fitFcn->IntegralError(2,3.45,res->GetParams(), res->GetCovarianceMatrix().GetMatrixArray() ))/0.01;
    std::cout << "Fitted " << histoname << std::endl;
    std::cout << "Nb JPsi total: " << integral << std::endl;
 //   std::cout << "integral error: " << integralerror << std::endl;
   signalFcn->Draw("same");
   backFcn->SetParameters(&par[8]);
   backFcn->Draw("same");
   // draw the legend
    char str[50];
    sprintf(str, "N_{JPsi} %f +/- %f", integral, integralerror);
   pave->AddText(str);
   pave->Draw();
   TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
   legend->SetTextFont(72);
   legend->SetTextSize(0.04);
    Char_t message[80];
    sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcn->GetChisquare(),fitFcn->GetNDF());
     legend->AddEntry(fitFcn,message);
   legend->AddEntry(histo,"Data","lpe");
   legend->AddEntry(backFcn,"Background fit","l");
   legend->AddEntry(signalFcn,"Signal fit","l");
   legend->Draw();
    
//    cout << "Histogramme: " << histoname << endl;
//    for(int j=0; j<10; j++){
//        cout << "Bin de masse numero " << j << " - Signal : " << (signalFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 << " et Background : " << (backFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 <<endl;
//    }
    
    return res;
}

