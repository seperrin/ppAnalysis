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
# include "TFitter.h"
# include "TProfile.h"
# include "TMinuit.h"
# include "TRandom.h"
# include <iostream>
# include <fstream>
# include <string>
# include "TH1F.h"
# include "TH1I.h"
# include "TH1D.h"
# include "TH2F.h"
# include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"

//Plot de différentes observables CINT et Obtention des limites PercentileMethod
// FUNCTIONS

void PlotFromTreeCINT();


void PlotFromTreeCINT(){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    
// ************************************
// Définitions de paramètres          *
// ************************************

    double ZvtxCut = 10;
    double DPhiCut = 0.01;
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 99999;
    double MinMultCentral = 37;
    double MaxMultPeriph = 23;
    double CentSPDTrackletsCentral = 1;
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
    int barWidth = 50;
    
    Char_t Group_Period[50] = "Group1_CINT_AllEst";//"Group12_LHC18m_CINT7CENT"
    Char_t *arrayOfPeriods[] = {"Group1_LHC16h_CINT_AllEst","Group1_LHC16j_CINT_AllEst","Group1_LHC16k_CINT_AllEst_P1","Group1_LHC16k_CINT_AllEst_Part2","Group1_LHC16o_CINT_AllEst","Group1_LHC16p_CINT_AllEst","Group1_LHC17i_CINT_AllEst","Group1_LHC17k_CINT_AllEst","Group1_LHC17l_CINT_AllEst"};
    //Char_t *arrayOfPeriods[] = {"Group6_LHC18m_CINT_AllEst"};
    int numberOfPeriods = sizeof(arrayOfPeriods) / sizeof(arrayOfPeriods[0]);
    
    const double binsCent[15] = {0,1,3,5,10,15,20,30,40,50,60,70,80,90,100};
    Char_t filePMLim[300];
    Char_t filePMLimSaveName[300];
    Char_t filePM[300];
    Char_t filePMLimSPDTracklets[300];
    Char_t filePMLimSPDTrackletsSaveName[300];
    Char_t filePMSPDTracklets[300];
    Char_t filePMLimSPDClusters[300];
    Char_t filePMLimSPDClustersSaveName[300];
    Char_t filePMSPDClusters[300];
    Char_t filePMLimZNC[300];
    Char_t filePMLimZNCSaveName[300];
    Char_t filePMZNC[300];
    Char_t filePMLimZNA[300];
    Char_t filePMLimZNASaveName[300];
    Char_t filePMZNA[300];
    Char_t filePMLimZNAC[300];
    Char_t filePMLimZNACSaveName[300];
    Char_t filePMZNAC[300];
    Char_t filePMLimV0C[300];
    Char_t filePMLimV0CSaveName[300];
    Char_t filePMV0C[300];
    Char_t filePMLimV0A[300];
    Char_t filePMLimV0ASaveName[300];
    Char_t filePMV0A[300];
    Char_t filePMLimV0M[300];
    Char_t filePMLimV0MSaveName[300];
    Char_t filePMV0M[300];
    Char_t filePMLimADC[300];
    Char_t filePMLimADCSaveName[300];
    Char_t filePMADC[300];
    Char_t filePMLimADA[300];
    Char_t filePMLimADASaveName[300];
    Char_t filePMADA[300];
    Char_t filePMLimADM[300];
    Char_t filePMLimADMSaveName[300];
    Char_t filePMADM[300];
    Char_t fileNmean[300];
    Char_t fileNmeanROOT[300];
    Char_t filespdtkl[300];
    Char_t fileMeanSPDClusters[300];
    Char_t fileMeanZNC[300];
    Char_t fileMeanZNA[300];
    Char_t fileMeanZNAC[300];
    Char_t fileMeanV0C[300];
    Char_t fileMeanV0A[300];
    Char_t fileMeanV0M[300];
    Char_t fileMeanADC[300];
    Char_t fileMeanADA[300];
    Char_t fileMeanADM[300];
    Char_t fileInLoc[300];
// *************************
// Initialiser les graphes *
// *************************
    
    TH1* PMSliced[20]={NULL};
    TH1* PMSlicedSPDTracklets[20]={NULL};
    TH1* PMSlicedSPDClusters[20]={NULL};
    TH1* PMSlicedZNC[20]={NULL};
    TH1* PMSlicedZNA[20]={NULL};
    TH1* PMSlicedZNAC[20]={NULL};
    TH1* PMSlicedV0C[20]={NULL};
    TH1* PMSlicedV0A[20]={NULL};
    TH1* PMSlicedV0M[20]={NULL};
    TH1* PMSlicedADC[20]={NULL};
    TH1* PMSlicedADA[20]={NULL};
    TH1* PMSlicedADM[20]={NULL};
    TH2F* PercentileMethodSPDTracklets(NULL);
    TH2F* PercentileMethodLimitsSPDTracklets(NULL);
    TH2F* PercentileMethod(NULL);
    TH2F* PercentileMethodLimits(NULL);
    TH2F* PercentileMethodSPDClusters(NULL);
    TH2F* PercentileMethodLimitsSPDClusters(NULL);
    TH2F* PercentileMethodZNC(NULL);
    TH2F* PercentileMethodLimitsZNC(NULL);
    TH2F* PercentileMethodZNA(NULL);
    TH2F* PercentileMethodLimitsZNA(NULL);
    TH2F* PercentileMethodZNAC(NULL);
    TH2F* PercentileMethodLimitsZNAC(NULL);
    TH2F* PercentileMethodV0C(NULL);
    TH2F* PercentileMethodLimitsV0C(NULL);
    TH2F* PercentileMethodV0A(NULL);
    TH2F* PercentileMethodLimitsV0A(NULL);
    TH2F* PercentileMethodV0M(NULL);
    TH2F* PercentileMethodLimitsV0M(NULL);
    TH2F* PercentileMethodADC(NULL);
    TH2F* PercentileMethodLimitsADC(NULL);
    TH2F* PercentileMethodADA(NULL);
    TH2F* PercentileMethodLimitsADA(NULL);
    TH2F* PercentileMethodADM(NULL);
    TH2F* PercentileMethodLimitsADM(NULL);
   // TH1F* hnseg_raw(NULL);
   // TH1I* hnseg_raw_num(NULL);
    TH1F* hnseg(NULL);
    TH1I* hnseg_num(NULL);
    TH1I* hnseg_den(NULL);
   // TH1F* hnseg_spdtkl(NULL);
    //TH1I* hnseg_spdtkl_num(NULL);
    //TH1F* hnseg_spdtklval(NULL);
    //TH1I* hnseg_spdtklval_num(NULL);
    //TH1I* hnseg_spdtklval_dist(NULL);
    //TH1I* hnseg_Ntkl_dist(NULL);
    TProfile* hn_MeanSPDTracklets(NULL);
    TProfile* hn_MeanSPDClusters(NULL);
    TProfile* hn_MeanZNC(NULL);
    TProfile* hn_MeanZNA(NULL);
    TProfile* hn_MeanZNAC(NULL);
    TProfile* hn_MeanV0C(NULL);
    TProfile* hn_MeanV0A(NULL);
    TProfile* hn_MeanV0M(NULL);
    TProfile* hn_MeanADC(NULL);
    TProfile* hn_MeanADA(NULL);
    TProfile* hn_MeanADM(NULL);
    
    
    hnseg = new TH1F("hnseg",
                     "Mean Ntkl wrt Zvtx",
                     81,-10,10);
    hnseg->SetXTitle("Zvtx");
    hnseg->SetYTitle("Mean Ntkl");
    
    hnseg_num = new TH1I("hnseg_num",
                     "Mean Ntkl wrt Zvtx - Numerator",
                     81,-10,10);
    hnseg_num->SetXTitle("Zvtx");
    hnseg_num->SetYTitle("Sum Ntkl");
    
    hnseg_den = new TH1I("hnseg_den",
                     "Zvtx distribution",
                     81,-10,10);
    hnseg_den->SetXTitle("Zvtx");
    hnseg_den->SetYTitle("Count");
    
    
    
    PercentileMethod = new TH2F("PercentileMethod",
                     "Ntkl wrt Zvtx",
                     20,-10,10,150,0,150);
    PercentileMethod->SetXTitle("Zvtx");
    PercentileMethod->SetYTitle("Ntkl");
    
    PercentileMethodLimits = new TH2F("PercentileMethodLimits",
                     "Ntkl Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimits->SetXTitle("Zvtx");
    PercentileMethodLimits->SetYTitle("Percentile");
    
    PercentileMethodSPDTracklets = new TH2F("PercentileMethodSPDTracklets",
                     "SPDTracklets wrt Zvtx",
                     20,-10,10,150,0,150);
    PercentileMethodSPDTracklets->SetXTitle("Zvtx");
    PercentileMethodSPDTracklets->SetYTitle("SPDTracklets");
    
    PercentileMethodLimitsSPDTracklets = new TH2F("PercentileMethodLimitsSPDTracklets",
                     "SPDTracklets Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsSPDTracklets->SetXTitle("Zvtx");
    PercentileMethodLimitsSPDTracklets->SetYTitle("Percentile");
    
    PercentileMethodSPDClusters = new TH2F("PercentileMethodSPDClusters",
                     "SPDClusters wrt Zvtx",
                     20,-10,10,1000,0,1000);
    PercentileMethodSPDClusters->SetXTitle("Zvtx");
    PercentileMethodSPDClusters->SetYTitle("SPDClusters");
    
    PercentileMethodLimitsSPDClusters = new TH2F("PercentileMethodLimitsSPDClusters",
                     "SPDClusters Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsSPDClusters->SetXTitle("Zvtx");
    PercentileMethodLimitsSPDClusters->SetYTitle("Percentile");
    
    PercentileMethodZNC = new TH2F("PercentileMethodZNC",
                     "ZNC wrt Zvtx",
                     20,-10,10,2000,0,2000);
    PercentileMethodZNC->SetXTitle("Zvtx");
    PercentileMethodZNC->SetYTitle("ZNC");
    
    PercentileMethodLimitsZNC = new TH2F("PercentileMethodLimitsZNC",
                     "ZNC Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsZNC->SetXTitle("Zvtx");
    PercentileMethodLimitsZNC->SetYTitle("Percentile");
    
    PercentileMethodZNA = new TH2F("PercentileMethodZNA",
                     "ZNA wrt Zvtx",
                     20,-10,10,2000,0,2000);
    PercentileMethodZNA->SetXTitle("Zvtx");
    PercentileMethodZNA->SetYTitle("ZNA");
    
    PercentileMethodLimitsZNA = new TH2F("PercentileMethodLimitsZNA",
                     "ZNA Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsZNA->SetXTitle("Zvtx");
    PercentileMethodLimitsZNA->SetYTitle("Percentile");
    
    PercentileMethodZNAC = new TH2F("PercentileMethodZNAC",
                     "ZNAC wrt Zvtx",
                     20,-10,10,2000,0,2000);
    PercentileMethodZNAC->SetXTitle("Zvtx");
    PercentileMethodZNAC->SetYTitle("ZNAC");
    
    PercentileMethodLimitsZNAC = new TH2F("PercentileMethodLimitsZNAC",
                     "ZNAC Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsZNAC->SetXTitle("Zvtx");
    PercentileMethodLimitsZNAC->SetYTitle("Percentile");
    
    PercentileMethodV0C = new TH2F("PercentileMethodV0C",
                     "V0C wrt Zvtx",
                     20,-10,10,2000,0,2000);
    PercentileMethodV0C->SetXTitle("Zvtx");
    PercentileMethodV0C->SetYTitle("V0C");
    
    PercentileMethodLimitsV0C = new TH2F("PercentileMethodLimitsV0C",
                     "V0C Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsV0C->SetXTitle("Zvtx");
    PercentileMethodLimitsV0C->SetYTitle("Percentile");
    
    PercentileMethodV0A = new TH2F("PercentileMethodV0A",
                     "V0A wrt Zvtx",
                     20,-10,10,2000,0,2000);
    PercentileMethodV0A->SetXTitle("Zvtx");
    PercentileMethodV0A->SetYTitle("V0A");
    
    PercentileMethodLimitsV0A = new TH2F("PercentileMethodLimitsV0A",
                     "V0A Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsV0A->SetXTitle("Zvtx");
    PercentileMethodLimitsV0A->SetYTitle("Percentile");
    
    PercentileMethodV0M = new TH2F("PercentileMethodV0M",
                     "V0M wrt Zvtx",
                     20,-10,10,2000,0,2000);
    PercentileMethodV0M->SetXTitle("Zvtx");
    PercentileMethodV0M->SetYTitle("V0M");
    
    PercentileMethodLimitsV0M = new TH2F("PercentileMethodLimitsV0M",
                     "V0M Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsV0M->SetXTitle("Zvtx");
    PercentileMethodLimitsV0M->SetYTitle("Percentile");
    
    PercentileMethodADC = new TH2F("PercentileMethodADC",
                     "ADC wrt Zvtx",
                     20,-10,10,10000,0,10000);
    PercentileMethodADC->SetXTitle("Zvtx");
    PercentileMethodADC->SetYTitle("ADC");
    
    PercentileMethodLimitsADC = new TH2F("PercentileMethodLimitsADC",
                     "ADC Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsADC->SetXTitle("Zvtx");
    PercentileMethodLimitsADC->SetYTitle("Percentile");
    
    PercentileMethodADA = new TH2F("PercentileMethodADA",
                     "ADA wrt Zvtx",
                     20,-10,10,10000,0,10000);
    PercentileMethodADA->SetXTitle("Zvtx");
    PercentileMethodADA->SetYTitle("ADA");
    
    PercentileMethodLimitsADA = new TH2F("PercentileMethodLimitsADA",
                     "ADA Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsADA->SetXTitle("Zvtx");
    PercentileMethodLimitsADA->SetYTitle("Percentile");
    
    PercentileMethodADM = new TH2F("PercentileMethodADM",
                     "ADM wrt Zvtx",
                     20,-10,10,10000,0,10000);
    PercentileMethodADM->SetXTitle("Zvtx");
    PercentileMethodADM->SetYTitle("ADM");
    
    PercentileMethodLimitsADM = new TH2F("PercentileMethodLimitsADM",
                     "ADM Limit wrt Zvtx and percentile",
                                      20,-10,10,14,binsCent);
    PercentileMethodLimitsADM->SetXTitle("Zvtx");
    PercentileMethodLimitsADM->SetYTitle("Percentile");
    
    
    
//    hnseg_raw = new TH1F("hnseg_raw",
//                     "Mean raw Ntkl wrt Zvtx",
//                     81,-10,10);
//    hnseg_raw->SetXTitle("Zvtx");
//    hnseg_raw->SetYTitle("Mean raw Ntkl");
//
//    hnseg_raw_num = new TH1I("hnseg_raw_num",
//                     "Mean raw Ntkl wrt Zvtx - Numerator",
//                     81,-10,10);
//    hnseg_raw_num->SetXTitle("Zvtx");
//    hnseg_raw_num->SetYTitle("Mean raw Ntkl");
    
//    hnseg_spdtkl = new TH1F("hnseg_spdtkl",
//                     "Mean SPDTracklets wrt Zvtx",
//                     81,-10,10);
//    hnseg_spdtkl->SetXTitle("Zvtx");
//    hnseg_spdtkl->SetYTitle("Mean SPDTracklets");
//
//    hnseg_spdtkl_num = new TH1I("hnseg_spdtkl_num",
//                     "Mean SPDTracklets wrt Zvtx - Numerator",
//                     81,-10,10);
//    hnseg_spdtkl_num->SetXTitle("Zvtx");
//    hnseg_spdtkl_num->SetYTitle("Mean SPDTracklets");
    
//    hnseg_spdtklval = new TH1F("hnseg_spdtklval",
//                     "Mean SPDTrackletsValue wrt Zvtx",
//                     81,-10,10);
//    hnseg_spdtklval->SetXTitle("Zvtx");
//    hnseg_spdtklval->SetYTitle("Mean SPDTrackletsValue");
//
//    hnseg_spdtklval_num = new TH1I("hnseg_spdtklval_num",
//                        "Mean SPDTrackletsValue wrt Zvtx - Numerator",
//                        81,-10,10);
//       hnseg_spdtklval_num->SetXTitle("Zvtx");
//       hnseg_spdtklval_num->SetYTitle("Mean SPDTrackletsValue");
    
//    hnseg_spdtklval_dist = new TH1I("hnseg_spdtklval_dist",
//                     "Distribution of SPDTrackletsValue",
//                     160,0,160);
//    hnseg_spdtklval_dist->SetXTitle("SPDTrackletsValue");
//    hnseg_spdtklval_dist->SetYTitle("Count");
//
//    hnseg_Ntkl_dist = new TH1I("hnseg_Ntkl_dist",
//                     "Distribution of Ntkl value",
//                     160,0,160);
//    hnseg_Ntkl_dist->SetXTitle("Ntkl value");
//    hnseg_Ntkl_dist->SetYTitle("Count");
    
    hn_MeanSPDTracklets = new TProfile("hn_MeanSPDTracklets",
                     "Mean SPDTrackletsValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanSPDTracklets->SetXTitle("Zvtx");
    hn_MeanSPDTracklets->SetYTitle("Mean SPDTrackletsValue");
    
    hn_MeanSPDClusters = new TProfile("hn_MeanSPDClusters",
                     "Mean SPDClustersValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanSPDClusters->SetXTitle("Zvtx");
    hn_MeanSPDClusters->SetYTitle("Mean SPDClustersValue");
    
    hn_MeanZNC = new TProfile("hn_MeanZNC",
                     "Mean ZNCValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanZNC->SetXTitle("Zvtx");
    hn_MeanZNC->SetYTitle("Mean ZNCValue");
    
    hn_MeanZNA = new TProfile("hn_MeanZNA",
                     "Mean ZNAValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanZNA->SetXTitle("Zvtx");
    hn_MeanZNA->SetYTitle("Mean ZNAValue");
    
    hn_MeanZNAC = new TProfile("hn_MeanZNAC",
                     "Mean ZNACValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanZNAC->SetXTitle("Zvtx");
    hn_MeanZNAC->SetYTitle("Mean ZNACValue");
    
    hn_MeanV0C = new TProfile("hn_MeanV0C",
                     "Mean V0CValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanV0C->SetXTitle("Zvtx");
    hn_MeanV0C->SetYTitle("Mean V0CValue");
    
    hn_MeanV0A = new TProfile("hn_MeanV0A",
                     "Mean V0AValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanV0A->SetXTitle("Zvtx");
    hn_MeanV0A->SetYTitle("Mean V0AValue");
    
    hn_MeanV0M = new TProfile("hn_MeanV0M",
                     "Mean V0MValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanV0M->SetXTitle("Zvtx");
    hn_MeanV0M->SetYTitle("Mean V0MValue");
    
    hn_MeanADC = new TProfile("hn_MeanADC",
                     "Mean ADCValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanADC->SetXTitle("Zvtx");
    hn_MeanADC->SetYTitle("Mean ADCValue");
    
    hn_MeanADA = new TProfile("hn_MeanADA",
                     "Mean ADAValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanADA->SetXTitle("Zvtx");
    hn_MeanADA->SetYTitle("Mean ADAValue");
    
    hn_MeanADM = new TProfile("hn_MeanADM",
                     "Mean ADMValue wrt Zvtx TPro",
                     81,-10,10);
    hn_MeanADM->SetXTitle("Zvtx");
    hn_MeanADM->SetYTitle("Mean ADMValue");

// *************************
// Analyse                 *
// *************************
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        sprintf(filePMLim,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLim.pdf",Group_Period,Group_Period);
        sprintf(filePMLimSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLim.txt",Group_Period,Group_Period);
        sprintf(filePM,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PM.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimSPDTracklets,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimSPDTracklets.pdf",Group_Period,Group_Period);
        sprintf(filePMLimSPDTrackletsSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimSPDTracklets.txt",Group_Period,Group_Period);
        sprintf(filePMSPDTracklets,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMSPDTracklets.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimSPDClusters,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimSPDClusters.pdf",Group_Period,Group_Period);
        sprintf(filePMLimSPDClustersSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimSPDClusters.txt",Group_Period,Group_Period);
        sprintf(filePMSPDClusters,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMSPDClusters.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimZNC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimZNC.pdf",Group_Period,Group_Period);
        sprintf(filePMLimZNCSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimZNC.txt",Group_Period,Group_Period);
        sprintf(filePMZNC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMZNC.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimZNA,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimZNA.pdf",Group_Period,Group_Period);
        sprintf(filePMLimZNASaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimZNA.txt",Group_Period,Group_Period);
        sprintf(filePMZNA,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMZNA.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimZNAC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimZNAC.pdf",Group_Period,Group_Period);
        sprintf(filePMLimZNACSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimZNAC.txt",Group_Period,Group_Period);
        sprintf(filePMZNAC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMZNAC.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimV0C,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimV0C.pdf",Group_Period,Group_Period);
        sprintf(filePMLimV0CSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimV0C.txt",Group_Period,Group_Period);
        sprintf(filePMV0C,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMV0C.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimV0A,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimV0A.pdf",Group_Period,Group_Period);
        sprintf(filePMLimV0ASaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimV0A.txt",Group_Period,Group_Period);
        sprintf(filePMV0A,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMV0A.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimV0M,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimV0M.pdf",Group_Period,Group_Period);
        sprintf(filePMLimV0MSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimV0M.txt",Group_Period,Group_Period);
        sprintf(filePMV0M,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMV0M.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimADC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimADC.pdf",Group_Period,Group_Period);
        sprintf(filePMLimADCSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimADC.txt",Group_Period,Group_Period);
        sprintf(filePMADC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMADC.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimADA,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimADA.pdf",Group_Period,Group_Period);
        sprintf(filePMLimADASaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimADA.txt",Group_Period,Group_Period);
        sprintf(filePMADA,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMADA.pdf",Group_Period,Group_Period);
        
        sprintf(filePMLimADM,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimADM.pdf",Group_Period,Group_Period);
        sprintf(filePMLimADMSaveName,"PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMLimADM.txt",Group_Period,Group_Period);
        sprintf(filePMADM,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_PMADM.pdf",Group_Period,Group_Period);
        
        
        
        sprintf(fileNmean,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_Nmean.pdf",Group_Period,Group_Period);
        sprintf(fileNmeanROOT,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_Nmean.root",Group_Period,Group_Period);
        
        sprintf(filespdtkl,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanSPDTrackletsValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanSPDClusters,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanSPDClustersValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanZNC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanZNCValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanZNA,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanZNAValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanZNAC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanZNACValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanV0C,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanV0CValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanV0A,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanV0AValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanV0M,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanV0MValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanADC,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanADCValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanADA,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanADAValue.pdf",Group_Period,Group_Period);
        
        sprintf(fileMeanADM,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/PercentileMethodStudy/AccountNtklNonZero/NewAnalysis/%s/%s_MeanADMValue.pdf",Group_Period,Group_Period);
        
       
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CINT/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        
    // TFile fileIn("~/../../Volumes/Transcend2/ppAnalysis/Scripts/muonGrid.root");
 //   TFile fileIn("~/../../Volumes/Transcend2/ppAnalysis/Scripts/Group12_LHC18m_CINT7CENT/muonGrid.root");
    TFile fileIn(fileInLoc);
   // TFile fileIn("~/../../Volumes/Transcend2/ppAnalysis/Scripts/CINT_16o_PS_CutsEvent_Val_CENT/muonGrid.root");
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
    
    
// ************************************
// Boucle sur les evenements, notés i *
// ************************************
    
  //  Double_t px[1000], py[1000];
        for (Int_t i=0;i<nevent;i++) {
            if(i%100000 == 0){
                    std::cout << "[";
                    double portion = double(i)/nevent;
                    long pos = barWidth * portion;
                    for (int k = 0; k < barWidth; ++k) {
                        if (k < pos) std::cout << "=";
                        else if (k == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                    std::cout << "] " << long(100 * portion) << "%     " << i << "/" << nevent << " Tree " << tree_idx << "/" << numberOfPeriods;
                    std::cout.flush();
                std::cout << std::endl;
            }
            theTree->GetEvent(i);
     //       cout << "Looking at Event " << i << endl;
//            cout << "fPassPhysicsSelection : " << fEvent->fPassPhysicsSelection << endl;
//            cout << "fPassTriggerSelection : " << fEvent->fPassTriggerSelection << endl;
//            cout << "fTriggerCMUL7 : " << fEvent->fTriggerCMUL7 << endl;
//            cout << "fTriggerCINT7 : " << fEvent->fTriggerCINT7 << endl;
            //cout << "tracks->GetEntries() " << tracks->GetEntries() << endl;
          //  cout << "Computing the number of tracklets in Event " << i << endl;
            
          //  cout << "Reading Event " << i <<endl;
            // Apply a cut on the TrackletEta at +-1
            int NumberCloseEtaTracklets = 0;
            
        //    if(NumberCloseEtaTracklets>0 && fEvent->fVertexNC >= 1 && fEvent->fNPileupVtx == 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<0.25 && (TMath::Abs(fEvent->fVertexZ))<10){
            if(fEvent->fVertexNC > 0 && fEvent->fIsPileupFromSPDMultBins == 0 && fEvent->fSPDVertexSigmaZ<=0.25 && (TMath::Abs(fEvent->fVertexZ))<10){ // && fEvent->fNPileupVtx == 0
                
                for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                    trac = (TrackletLight*)fTracklets->At(j);
                    if((TMath::Abs(trac->fEta) < 1) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                            NumberCloseEtaTracklets++;
                        
                    }
                    
                }
                
                if(NumberCloseEtaTracklets>0){
                
                
              //  hnseg_spdtkl_num->Fill(fEvent->fVertexZ, fEvent->fCentralitySPDTracklets);
               //hnseg_spdtklval_num->Fill(fEvent->fVertexZ, fEvent->fSPDTrackletsValue);
                hn_MeanSPDTracklets->Fill(fEvent->fVertexZ, fEvent->fSPDTrackletsValue);
                    hn_MeanSPDClusters->Fill(fEvent->fVertexZ, fEvent->fSPDClustersValue);
                    hn_MeanZNC->Fill(fEvent->fVertexZ, fEvent->fZNCValue);
                    hn_MeanZNA->Fill(fEvent->fVertexZ, fEvent->fZNAValue);
                    hn_MeanZNAC->Fill(fEvent->fVertexZ, fEvent->fZNACValue);
                    hn_MeanV0C->Fill(fEvent->fVertexZ, fEvent->fV0CValue);
                    hn_MeanV0A->Fill(fEvent->fVertexZ, fEvent->fV0AValue);
                    hn_MeanV0M->Fill(fEvent->fVertexZ, fEvent->fV0MValue);
                    hn_MeanADC->Fill(fEvent->fVertexZ, fEvent->fADCValue);
                    hn_MeanADA->Fill(fEvent->fVertexZ, fEvent->fADAValue);
                    hn_MeanADM->Fill(fEvent->fVertexZ, fEvent->fADMValue);
                hnseg_den->Fill(fEvent->fVertexZ);
                //hnseg_raw_num->Fill(fEvent->fVertexZ, fEvent->fNTracklets);
                hnseg_num->Fill(fEvent->fVertexZ, NumberCloseEtaTracklets);
                PercentileMethod->Fill(fEvent->fVertexZ,NumberCloseEtaTracklets);
                PercentileMethodSPDTracklets->Fill(fEvent->fVertexZ,fEvent->fSPDTrackletsValue);
                    PercentileMethodSPDClusters->Fill(fEvent->fVertexZ,fEvent->fSPDClustersValue);
                    PercentileMethodZNC->Fill(fEvent->fVertexZ,fEvent->fZNCValue);
                    PercentileMethodZNA->Fill(fEvent->fVertexZ,fEvent->fZNAValue);
                    PercentileMethodZNAC->Fill(fEvent->fVertexZ,fEvent->fZNACValue);
                    PercentileMethodV0C->Fill(fEvent->fVertexZ,fEvent->fV0CValue);
                    PercentileMethodV0A->Fill(fEvent->fVertexZ,fEvent->fV0AValue);
                    PercentileMethodV0M->Fill(fEvent->fVertexZ,fEvent->fV0MValue);
                    PercentileMethodADC->Fill(fEvent->fVertexZ,fEvent->fADCValue);
                    PercentileMethodADA->Fill(fEvent->fVertexZ,fEvent->fADAValue);
                    PercentileMethodADM->Fill(fEvent->fVertexZ,fEvent->fADMValue);
                    
                //hnseg_Ntkl_dist->Fill(NumberCloseEtaTracklets);
               }
           // hnseg_spdtklval_dist->Fill(fEvent->fSPDTrackletsValue);
            
            // hnseg->Fill(fEvent->fVertexZ, NumberCloseEtaTracklets);
           //  hnseg_raw->Fill(fEvent->fVertexZ, fEvent->fNTracklets);
            
            }
        }
    }
    
    TCanvas*cPMproj=new TCanvas();
    cPMproj->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMproj->cd(i);
        PMSliced[i] = PercentileMethod->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSliced[i]->Draw("e");
    }
    
    cout << "PMSliced[0] info: " << PMSliced[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSliced[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSliced[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimits->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimits->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimits->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimits->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimits->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimits->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimits->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimits->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimits->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimits->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimits->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimits->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimits->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimSave(filePMLimSaveName, std::ofstream::out);
    
    filePMLimSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimSave << PercentileMethodLimits->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimSave << "},\n";
        }
        else{
            filePMLimSave << "}   }\n";
        }
    }
    cout <<filePMLimSave.is_open() <<endl;
    
    filePMLimSave.close();
    
    
    
    
    
    // Percentile Method SPDTracklets
    
    TCanvas*cPMprojSPDTracklets=new TCanvas();
    cPMprojSPDTracklets->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojSPDTracklets->cd(i);
        PMSlicedSPDTracklets[i] = PercentileMethodSPDTracklets->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedSPDTracklets[i]->Draw("e");
    }
    
    cout << "PMSlicedSPDTracklets[0] info: " << PMSlicedSPDTracklets[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedSPDTracklets[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedSPDTracklets[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsSPDTracklets->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimSPDTrackletsSave(filePMLimSPDTrackletsSaveName, std::ofstream::out);
    
    filePMLimSPDTrackletsSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimSPDTrackletsSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimSPDTrackletsSave << PercentileMethodLimitsSPDTracklets->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimSPDTrackletsSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimSPDTrackletsSave << "},\n";
        }
        else{
            filePMLimSPDTrackletsSave << "}   }\n";
        }
    }
    cout <<filePMLimSPDTrackletsSave.is_open() <<endl;
    
    filePMLimSPDTrackletsSave.close();
    
    
    
    // Percentile Method SPDClusters
    
    TCanvas*cPMprojSPDClusters=new TCanvas();
    cPMprojSPDClusters->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojSPDClusters->cd(i);
        PMSlicedSPDClusters[i] = PercentileMethodSPDClusters->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedSPDClusters[i]->Draw("e");
    }
    
    cout << "PMSlicedSPDClusters[0] info: " << PMSlicedSPDClusters[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedSPDClusters[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedSPDClusters[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsSPDClusters->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimSPDClustersSave(filePMLimSPDClustersSaveName, std::ofstream::out);
    
    filePMLimSPDClustersSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimSPDClustersSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimSPDClustersSave << PercentileMethodLimitsSPDClusters->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimSPDClustersSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimSPDClustersSave << "},\n";
        }
        else{
            filePMLimSPDClustersSave << "}   }\n";
        }
    }
    cout <<filePMLimSPDClustersSave.is_open() <<endl;
    
    filePMLimSPDClustersSave.close();
    
    
    
    // Percentile Method ZNC
    
    TCanvas*cPMprojZNC=new TCanvas();
    cPMprojZNC->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojZNC->cd(i);
        PMSlicedZNC[i] = PercentileMethodZNC->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedZNC[i]->Draw("e");
    }
    
    cout << "PMSlicedZNC[0] info: " << PMSlicedZNC[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedZNC[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedZNC[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsZNC->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsZNC->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsZNC->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsZNC->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsZNC->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsZNC->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsZNC->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsZNC->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsZNC->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsZNC->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsZNC->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsZNC->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsZNC->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimZNCSave(filePMLimZNCSaveName, std::ofstream::out);
    
    filePMLimZNCSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimZNCSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimZNCSave << PercentileMethodLimitsZNC->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimZNCSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimZNCSave << "},\n";
        }
        else{
            filePMLimZNCSave << "}   }\n";
        }
    }
    cout <<filePMLimZNCSave.is_open() <<endl;
    
    filePMLimZNCSave.close();
    
    
    
    // Percentile Method ZNA
    
    TCanvas*cPMprojZNA=new TCanvas();
    cPMprojZNA->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojZNA->cd(i);
        PMSlicedZNA[i] = PercentileMethodZNA->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedZNA[i]->Draw("e");
    }
    
    cout << "PMSlicedZNA[0] info: " << PMSlicedZNA[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedZNA[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedZNA[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsZNA->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsZNA->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsZNA->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsZNA->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsZNA->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsZNA->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsZNA->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsZNA->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsZNA->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsZNA->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsZNA->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsZNA->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsZNA->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimZNASave(filePMLimZNASaveName, std::ofstream::out);
    
    filePMLimZNASave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimZNASave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimZNASave << PercentileMethodLimitsZNA->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimZNASave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimZNASave << "},\n";
        }
        else{
            filePMLimZNASave << "}   }\n";
        }
    }
    cout <<filePMLimZNASave.is_open() <<endl;
    
    filePMLimZNASave.close();
    
    
    
    // Percentile Method ZNAC
    
    TCanvas*cPMprojZNAC=new TCanvas();
    cPMprojZNAC->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojZNAC->cd(i);
        PMSlicedZNAC[i] = PercentileMethodZNAC->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedZNAC[i]->Draw("e");
    }
    
    cout << "PMSlicedZNAC[0] info: " << PMSlicedZNAC[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedZNAC[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedZNAC[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsZNAC->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimZNACSave(filePMLimZNACSaveName, std::ofstream::out);
    
    filePMLimZNACSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimZNACSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimZNACSave << PercentileMethodLimitsZNAC->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimZNACSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimZNACSave << "},\n";
        }
        else{
            filePMLimZNACSave << "}   }\n";
        }
    }
    cout <<filePMLimZNACSave.is_open() <<endl;
    
    filePMLimZNACSave.close();
    
    
    // Percentile Method V0C
    
    TCanvas*cPMprojV0C=new TCanvas();
    cPMprojV0C->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojV0C->cd(i);
        PMSlicedV0C[i] = PercentileMethodV0C->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedV0C[i]->Draw("e");
    }
    
    cout << "PMSlicedV0C[0] info: " << PMSlicedV0C[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedV0C[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedV0C[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsV0C->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsV0C->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsV0C->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsV0C->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsV0C->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsV0C->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsV0C->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsV0C->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsV0C->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsV0C->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsV0C->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsV0C->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsV0C->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimV0CSave(filePMLimV0CSaveName, std::ofstream::out);
    
    filePMLimV0CSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimV0CSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimV0CSave << PercentileMethodLimitsV0C->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimV0CSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimV0CSave << "},\n";
        }
        else{
            filePMLimV0CSave << "}   }\n";
        }
    }
    cout <<filePMLimV0CSave.is_open() <<endl;
    
    filePMLimV0CSave.close();
    
    
    
    // Percentile Method V0A
    
    TCanvas*cPMprojV0A=new TCanvas();
    cPMprojV0A->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojV0A->cd(i);
        PMSlicedV0A[i] = PercentileMethodV0A->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedV0A[i]->Draw("e");
    }
    
    cout << "PMSlicedV0A[0] info: " << PMSlicedV0A[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedV0A[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedV0A[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsV0A->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsV0A->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsV0A->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsV0A->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsV0A->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsV0A->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsV0A->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsV0A->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsV0A->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsV0A->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsV0A->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsV0A->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsV0A->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimV0ASave(filePMLimV0ASaveName, std::ofstream::out);
    
    filePMLimV0ASave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimV0ASave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimV0ASave << PercentileMethodLimitsV0A->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimV0ASave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimV0ASave << "},\n";
        }
        else{
            filePMLimV0ASave << "}   }\n";
        }
    }
    cout <<filePMLimV0ASave.is_open() <<endl;
    
    filePMLimV0ASave.close();
    
    
    
    
    // Percentile Method V0M
    
    TCanvas*cPMprojV0M=new TCanvas();
    cPMprojV0M->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojV0M->cd(i);
        PMSlicedV0M[i] = PercentileMethodV0M->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedV0M[i]->Draw("e");
    }
    
    cout << "PMSlicedV0M[0] info: " << PMSlicedV0M[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedV0M[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedV0M[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsV0M->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsV0M->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsV0M->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsV0M->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsV0M->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsV0M->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsV0M->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsV0M->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsV0M->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsV0M->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsV0M->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsV0M->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsV0M->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimV0MSave(filePMLimV0MSaveName, std::ofstream::out);
    
    filePMLimV0MSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimV0MSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimV0MSave << PercentileMethodLimitsV0M->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimV0MSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimV0MSave << "},\n";
        }
        else{
            filePMLimV0MSave << "}   }\n";
        }
    }
    cout <<filePMLimV0MSave.is_open() <<endl;
    
    filePMLimV0MSave.close();
    
    
    
    // Percentile Method ADC
    
    TCanvas*cPMprojADC=new TCanvas();
    cPMprojADC->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojADC->cd(i);
        PMSlicedADC[i] = PercentileMethodADC->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedADC[i]->Draw("e");
    }
    
    cout << "PMSlicedADC[0] info: " << PMSlicedADC[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedADC[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedADC[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsADC->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsADC->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsADC->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsADC->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsADC->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsADC->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsADC->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsADC->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsADC->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsADC->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsADC->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsADC->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsADC->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimADCSave(filePMLimADCSaveName, std::ofstream::out);
    
    filePMLimADCSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimADCSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimADCSave << PercentileMethodLimitsADC->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimADCSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimADCSave << "},\n";
        }
        else{
            filePMLimADCSave << "}   }\n";
        }
    }
    cout <<filePMLimADCSave.is_open() <<endl;
    
    filePMLimADCSave.close();
    
    
    
    
    // Percentile Method ADA
    
    TCanvas*cPMprojADA=new TCanvas();
    cPMprojADA->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojADA->cd(i);
        PMSlicedADA[i] = PercentileMethodADA->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedADA[i]->Draw("e");
    }
    
    cout << "PMSlicedADA[0] info: " << PMSlicedADA[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedADA[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedADA[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsADA->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsADA->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsADA->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsADA->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsADA->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsADA->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsADA->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsADA->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsADA->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsADA->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsADA->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsADA->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsADA->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimADASave(filePMLimADASaveName, std::ofstream::out);
    
    filePMLimADASave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimADASave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimADASave << PercentileMethodLimitsADA->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimADASave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimADASave << "},\n";
        }
        else{
            filePMLimADASave << "}   }\n";
        }
    }
    cout <<filePMLimADASave.is_open() <<endl;
    
    filePMLimADASave.close();
    
    
    
    // Percentile Method ADM
    
    TCanvas*cPMprojADM=new TCanvas();
    cPMprojADM->Divide(4,5);
    for(int i=0;i<20;i++){
        cPMprojADM->cd(i);
        PMSlicedADM[i] = PercentileMethodADM->ProjectionY(Form("bin%d_0",i+1),i+1,i+1);
        PMSlicedADM[i]->Draw("e");
    }
    
    cout << "PMSlicedADM[0] info: " << PMSlicedADM[0]->GetEntries() <<endl;
    
    for(int i=0; i<20; i++){
        int nentries = PMSlicedADM[i]->GetEntries();
        int sum = 0;
        int j = 1;
        while(sum<nentries){
            sum += PMSlicedADM[i]->GetBinContent(j);
            double ratio = double(sum)/nentries;
            if(i==0){
                cout << "i: "<< i<<" j: " << j << " sum: " << sum << " ratio: " << ratio <<endl;
            }
            if(ratio<0.99){
                PercentileMethodLimitsADM->Fill((i-10)+0.5, 0.5);
                if(ratio<0.97){
                    PercentileMethodLimitsADM->Fill((i-10)+0.5, 2);
                    if(ratio<0.95){
                        PercentileMethodLimitsADM->Fill((i-10)+0.5, 4);
                        if(ratio<0.90){
                            PercentileMethodLimitsADM->Fill((i-10)+0.5, 7.5);
                            if(ratio<0.85){
                                PercentileMethodLimitsADM->Fill((i-10)+0.5, 12.5);
                                if(ratio<0.80){
                                    PercentileMethodLimitsADM->Fill((i-10)+0.5, 17.5);
                                    if(ratio<0.70){
                                        PercentileMethodLimitsADM->Fill((i-10)+0.5, 25);
                                        if(ratio<0.60){
                                            PercentileMethodLimitsADM->Fill((i-10)+0.5, 35);
                                            if(ratio<0.50){
                                                PercentileMethodLimitsADM->Fill((i-10)+0.5, 45);
                                                if(ratio<0.40){
                                                    PercentileMethodLimitsADM->Fill((i-10)+0.5, 55);
                                                    if(ratio<0.30){
                                                        PercentileMethodLimitsADM->Fill((i-10)+0.5, 65);
                                                        if(ratio<0.20){
                                                            PercentileMethodLimitsADM->Fill((i-10)+0.5, 75);
                                                            if(ratio<0.10){
                                                                PercentileMethodLimitsADM->Fill((i-10)+0.5, 85);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            j++;
        }
    }
    

    
    //Ajouter une boucle pour lire les résultats finaux de Percentile Method et les écrire dans un fichier
    
    cout << "Reading of events finished" <<endl;
    
    std::ofstream filePMLimADMSave(filePMLimADMSaveName, std::ofstream::out);
    
    filePMLimADMSave << "{   ";
    for(int z_idx=1; z_idx<21; z_idx++){
        filePMLimADMSave << "{";
        for(int cent_idx=1; cent_idx<15; cent_idx++){
            filePMLimADMSave << PercentileMethodLimitsADM->GetBinContent(z_idx, cent_idx);
            if(cent_idx!=14){
                filePMLimADMSave << ", ";
            }
        }
        if(z_idx!=20){
            filePMLimADMSave << "},\n";
        }
        else{
            filePMLimADMSave << "}   }\n";
        }
    }
    cout <<filePMLimADMSave.is_open() <<endl;
    
    filePMLimADMSave.close();

    std::cout << "==================== Analysis Terminated ====================" << std::endl;
    
    
    
// ***********************************
// Tracer les graphes sur les canevas *
// ***********************************


    TCanvas*c1=new TCanvas();
    c1->Divide(2,1);
    c1->cd(1);
    hnseg->Divide(hnseg_num, hnseg_den);
    hnseg->SetFillColor(18);
    hnseg->Draw();
//    c1->cd(2);
//    hnseg_raw->Divide(hnseg_den);
//    hnseg_raw->SetFillColor(18);
//    hnseg_raw->Draw();
    c1->cd(2);
    hnseg_den->SetFillColor(18);
    hnseg_den->Draw();
    
    TCanvas*c1solo=new TCanvas();
        hnseg->Draw();
    c1solo->SaveAs(fileNmean);
    
    TFile *outputFileNmean;
        outputFileNmean = new TFile(fileNmeanROOT,"RECREATE");
        hnseg->Write();
    
//    TCanvas*c2=new TCanvas();
//    hnseg_raw->Divide(hnseg_raw_num, hnseg_den);
//    hnseg_raw->SetFillColor(18);
//    hnseg_raw->Draw();
    
//    TCanvas*c3=new TCanvas();
//    hnseg_spdtkl->Divide(hnseg_spdtkl_num, hnseg_den);
//    hnseg_spdtkl->SetFillColor(18);
//    hnseg_spdtkl->Draw();
    
//    TCanvas*c4=new TCanvas();
//    hnseg_spdtklval->Divide(hnseg_spdtklval_num, hnseg_den);
//    hnseg_spdtklval->SetFillColor(18);
//    hnseg_spdtklval->Draw();
    
    TCanvas*c4Pro=new TCanvas();
    hn_MeanSPDTracklets->SetFillColor(18);
    hn_MeanSPDTracklets->Draw();
    c4Pro->SaveAs(filespdtkl);
    
    TCanvas*c4SPDClusters=new TCanvas();
    hn_MeanSPDClusters->SetFillColor(18);
    hn_MeanSPDClusters->Draw();
    c4SPDClusters->SaveAs(fileMeanSPDClusters);
    
    TCanvas*c4ZNC=new TCanvas();
    hn_MeanZNC->SetFillColor(18);
    hn_MeanZNC->Draw();
    c4ZNC->SaveAs(fileMeanZNC);
    
    TCanvas*c4ZNA=new TCanvas();
    hn_MeanZNA->SetFillColor(18);
    hn_MeanZNA->Draw();
    c4ZNA->SaveAs(fileMeanZNA);
    
    TCanvas*c4ZNAC=new TCanvas();
    hn_MeanZNAC->SetFillColor(18);
    hn_MeanZNAC->Draw();
    c4ZNAC->SaveAs(fileMeanZNAC);
    
    TCanvas*c4V0C=new TCanvas();
    hn_MeanV0C->SetFillColor(18);
    hn_MeanV0C->Draw();
    c4V0C->SaveAs(fileMeanV0C);
    
    TCanvas*c4V0A=new TCanvas();
    hn_MeanV0A->SetFillColor(18);
    hn_MeanV0A->Draw();
    c4V0A->SaveAs(fileMeanV0A);
    
    TCanvas*c4V0M=new TCanvas();
    hn_MeanV0M->SetFillColor(18);
    hn_MeanV0M->Draw();
    c4V0M->SaveAs(fileMeanV0M);
    
    TCanvas*c4ADC=new TCanvas();
    hn_MeanADC->SetFillColor(18);
    hn_MeanADC->Draw();
    c4ADC->SaveAs(fileMeanADC);
    
    TCanvas*c4ADA=new TCanvas();
    hn_MeanADA->SetFillColor(18);
    hn_MeanADA->Draw();
    c4ADA->SaveAs(fileMeanADA);
    
    TCanvas*c4ADM=new TCanvas();
    hn_MeanADM->SetFillColor(18);
    hn_MeanADM->Draw();
    c4ADM->SaveAs(fileMeanADM);
    
//    TCanvas*c45=new TCanvas();
//    hnseg_spdtklval_dist->SetFillColor(18);
//    hnseg_spdtklval_dist->Draw("e");
    
//    TCanvas*c456=new TCanvas();
//    hnseg_Ntkl_dist->SetFillColor(18);
//    hnseg_Ntkl_dist->Draw("e");
    
    TCanvas*cPM=new TCanvas();
    PercentileMethod->SetFillColor(18);
    PercentileMethod->Draw("colz");
    cPM->SaveAs(filePM);
    
    TCanvas*cPMLim=new TCanvas();
    PercentileMethodLimits->SetFillColor(18);
    PercentileMethodLimits->Draw("colz");
    cPMLim->SaveAs(filePMLim);
    
    TCanvas*cPMSPDTracklets=new TCanvas();
    PercentileMethodSPDTracklets->SetFillColor(18);
    PercentileMethodSPDTracklets->Draw("colz");
    cPMSPDTracklets->SaveAs(filePMSPDTracklets);
    
    TCanvas*cPMLimSPDTracklets=new TCanvas();
    PercentileMethodLimitsSPDTracklets->SetFillColor(18);
    PercentileMethodLimitsSPDTracklets->Draw("colz");
    cPMLimSPDTracklets->SaveAs(filePMLimSPDTracklets);
    
    TCanvas*cPMSPDClusters=new TCanvas();
    PercentileMethodSPDClusters->SetFillColor(18);
    PercentileMethodSPDClusters->Draw("colz");
    cPMSPDClusters->SaveAs(filePMSPDClusters);
    
    TCanvas*cPMLimSPDClusters=new TCanvas();
    PercentileMethodLimitsSPDClusters->SetFillColor(18);
    PercentileMethodLimitsSPDClusters->Draw("colz");
    cPMLimSPDClusters->SaveAs(filePMLimSPDClusters);
    
    TCanvas*cPMZNC=new TCanvas();
    PercentileMethodZNC->SetFillColor(18);
    PercentileMethodZNC->Draw("colz");
    cPMZNC->SaveAs(filePMZNC);
    
    TCanvas*cPMLimZNC=new TCanvas();
    PercentileMethodLimitsZNC->SetFillColor(18);
    PercentileMethodLimitsZNC->Draw("colz");
    cPMLimZNC->SaveAs(filePMLimZNC);
    
    TCanvas*cPMZNA=new TCanvas();
    PercentileMethodZNA->SetFillColor(18);
    PercentileMethodZNA->Draw("colz");
    cPMZNA->SaveAs(filePMZNA);
    
    TCanvas*cPMLimZNA=new TCanvas();
    PercentileMethodLimitsZNA->SetFillColor(18);
    PercentileMethodLimitsZNA->Draw("colz");
    cPMLimZNA->SaveAs(filePMLimZNA);
    
    TCanvas*cPMZNAC=new TCanvas();
    PercentileMethodZNAC->SetFillColor(18);
    PercentileMethodZNAC->Draw("colz");
    cPMZNAC->SaveAs(filePMZNAC);
    
    TCanvas*cPMLimZNAC=new TCanvas();
    PercentileMethodLimitsZNAC->SetFillColor(18);
    PercentileMethodLimitsZNAC->Draw("colz");
    cPMLimZNAC->SaveAs(filePMLimZNAC);
    
    TCanvas*cPMV0C=new TCanvas();
    PercentileMethodV0C->SetFillColor(18);
    PercentileMethodV0C->Draw("colz");
    cPMV0C->SaveAs(filePMV0C);
    
    TCanvas*cPMLimV0C=new TCanvas();
    PercentileMethodLimitsV0C->SetFillColor(18);
    PercentileMethodLimitsV0C->Draw("colz");
    cPMLimV0C->SaveAs(filePMLimV0C);
    
    TCanvas*cPMV0A=new TCanvas();
    PercentileMethodV0A->SetFillColor(18);
    PercentileMethodV0A->Draw("colz");
    cPMV0A->SaveAs(filePMV0A);
    
    TCanvas*cPMLimV0A=new TCanvas();
    PercentileMethodLimitsV0A->SetFillColor(18);
    PercentileMethodLimitsV0A->Draw("colz");
    cPMLimV0A->SaveAs(filePMLimV0A);
    
    TCanvas*cPMV0M=new TCanvas();
    PercentileMethodV0M->SetFillColor(18);
    PercentileMethodV0M->Draw("colz");
    cPMV0M->SaveAs(filePMV0M);
    
    TCanvas*cPMLimV0M=new TCanvas();
    PercentileMethodLimitsV0M->SetFillColor(18);
    PercentileMethodLimitsV0M->Draw("colz");
    cPMLimV0M->SaveAs(filePMLimV0M);
    
    TCanvas*cPMADC=new TCanvas();
    PercentileMethodADC->SetFillColor(18);
    PercentileMethodADC->Draw("colz");
    cPMADC->SaveAs(filePMADC);
    
    TCanvas*cPMLimADC=new TCanvas();
    PercentileMethodLimitsADC->SetFillColor(18);
    PercentileMethodLimitsADC->Draw("colz");
    cPMLimADC->SaveAs(filePMLimADC);
    
    TCanvas*cPMADA=new TCanvas();
    PercentileMethodADA->SetFillColor(18);
    PercentileMethodADA->Draw("colz");
    cPMADA->SaveAs(filePMADA);
    
    TCanvas*cPMLimADA=new TCanvas();
    PercentileMethodLimitsADA->SetFillColor(18);
    PercentileMethodLimitsADA->Draw("colz");
    cPMLimADA->SaveAs(filePMLimADA);
    
    TCanvas*cPMADM=new TCanvas();
    PercentileMethodADM->SetFillColor(18);
    PercentileMethodADM->Draw("colz");
    cPMADM->SaveAs(filePMADM);
    
    TCanvas*cPMLimADM=new TCanvas();
    PercentileMethodLimitsADM->SetFillColor(18);
    PercentileMethodLimitsADM->Draw("colz");
    cPMLimADM->SaveAs(filePMLimADM);
    
    
}
    
