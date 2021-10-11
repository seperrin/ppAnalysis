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

// Code d'analyse complet

// FUNCTIONS

void PlotFromTreeTKL(Char_t radical[200]);
Double_t FourierV2(Double_t *x,Double_t *par);
Double_t FourierV5(Double_t *x,Double_t *par);
int GetCent(float cent);
int GetCentPM(int Ntkl, float SPDTrackletsPer, float SPDClustersPer, float V0MPer, int zvtx_idx, int groupnumber);
int GetCentCvetan(float V0MPer);
bool isCentral(int centint);
bool isPeripheral(int centint);
Double_t myEtaLow(double x);
Double_t myEtaHigh(double x);

// Centrality bins:
//0: 0-1%
//1: 1-3%
//2: 3-5%
//3: 5-10%
//4: 10-15%
//5: 15-20%
//6: 20-30%
//7: 30-40%
//8: 40-50%
//9: 50-60%
//10: 60-70%
//11: 70-80%
//12: 80-90%
//13: 90-100%

int CentralLowBound = 0;
int CentralHighBound = 5;
int PeripheralLowBound = 8;
int PeripheralHighBound = 13;

string centralityMethod = "V0MPercentile";

bool isMultiplicityStudy = kFALSE;

int NtklCentralLowBound = 6;
int NtklCentralHighBound = 12;
int NtklPeripheralLowBound = 0;
int NtklPeripheralHighBound = 6;


bool isCentral(int centint){
    if((centint>=CentralLowBound)&&(centint<=CentralHighBound))
        return kTRUE;
    else{
        return kFALSE;
    }
}

bool isPeripheral(int centint){
    if((centint>=PeripheralLowBound)&&(centint<=PeripheralHighBound))
        return kTRUE;
    else{
        return kFALSE;
    }
}

int LimitsPM_SPDTracklets_Uncal_DataDriven[12][20][14] =
{       {   {27, 21, 18, 15, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}, //Group1
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {30, 23, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 23, 18, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {33, 27, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {34, 27, 23, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {35, 28, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {35, 28, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {36, 28, 25, 20, 16, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {36, 29, 25, 20, 17, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {36, 29, 25, 20, 17, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {37, 29, 25, 20, 17, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {36, 28, 25, 20, 16, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0}   },
        {   {26, 21, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}, //Group2
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {32, 25, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {34, 27, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {34, 27, 23, 19, 15, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {34, 27, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {35, 27, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {32, 26, 23, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0}   },
        {   {27, 21, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}, //Group3
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {32, 25, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {34, 27, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {34, 27, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {35, 27, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {35, 28, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {35, 28, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0}   },
        {   {25, 20, 17, 14, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0}, //Group4
            {27, 21, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 23, 18, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {34, 27, 23, 18, 15, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {32, 25, 22, 17, 15, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 22, 20, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {27, 21, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0}, //Group5
            {29, 23, 20, 16, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {33, 26, 23, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 1, 0},
            {33, 26, 23, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {33, 27, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {35, 28, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {36, 29, 25, 20, 16, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {37, 29, 25, 20, 17, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {38, 30, 26, 20, 17, 15, 11, 9, 7, 5, 4, 3, 2, 0},
            {38, 30, 26, 21, 17, 15, 11, 9, 7, 5, 4, 3, 2, 0},
            {39, 31, 27, 21, 18, 15, 11, 9, 7, 5, 4, 3, 2, 0},
            {39, 31, 27, 21, 18, 15, 11, 9, 7, 5, 4, 3, 2, 0},
            {38, 30, 26, 21, 17, 15, 11, 9, 7, 5, 4, 3, 2, 0},
            {36, 29, 25, 20, 16, 14, 11, 8, 6, 5, 4, 3, 2, 0},
            {34, 27, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {32, 25, 22, 17, 15, 12, 9, 7, 6, 4, 3, 2, 1, 0}   },
        {   {23, 19, 16, 13, 11, 9, 7, 5, 4, 3, 2, 2, 1, 0}, //Group6
            {25, 20, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {26, 20, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 22, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 20, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 23, 20, 16, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 15, 13, 9, 7, 6, 4, 3, 2, 1, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {34, 27, 23, 19, 15, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {34, 27, 24, 19, 16, 13, 10, 8, 6, 5, 4, 3, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {32, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 22, 20, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0}   },
       {   {24, 18, 16, 13, 11, 9, 7, 5, 4, 3, 3, 2, 1, 0}, //Group7
            {25, 20, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {26, 20, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 22, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 20, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 15, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {32, 25, 22, 17, 15, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {27, 21, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {24, 19, 16, 13, 11, 9, 7, 5, 4, 3, 3, 2, 1, 0}, //Group8
            {25, 20, 17, 14, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {26, 20, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 21, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 20, 16, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 22, 20, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 15, 13, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 21, 18, 15, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {23, 18, 16, 13, 11, 9, 7, 5, 4, 3, 2, 2, 1, 0}, //Group9
            {25, 20, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {26, 20, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 22, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 4, 3, 2, 2, 0},
            {33, 26, 23, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {32, 26, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 22, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {24, 19, 16, 13, 11, 9, 7, 5, 4, 3, 3, 2, 1, 0}, //Group10
            {25, 20, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {26, 21, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 22, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 20, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 23, 20, 16, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 22, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 20, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {25, 20, 17, 14, 12, 10, 8, 6, 4, 3, 3, 2, 1, 0}, //Group11
            {27, 21, 18, 15, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 22, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 20, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {25, 20, 17, 14, 12, 10, 8, 6, 5, 3, 3, 2, 1, 0}, //Group12
            {27, 21, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {29, 23, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 24, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 15, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 25, 22, 17, 15, 13, 10, 7, 6, 4, 3, 2, 1, 0},
            {31, 25, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 24, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 22, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 21, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   }

};


Int_t npfits;


void PlotFromTreeTKL(Char_t radical[200]){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    TH1::SetDefaultSumw2();
    bool doTracklets = kTRUE;
    bool isCMSMethod = kFALSE;
    bool doMixedEvents = kTRUE;
    bool KeepOnlyOne = kFALSE;
    int valueOnlyOne = 3;
    bool AdditionalCutNtkl = kFALSE;
    int valueAdditionalCutNtkl = 3; //<=
    int preciseZbinFocus = 10;
    int preciseCentFocus = 10;
    
    
// ************************************
// Définitions de paramètres          *
// ************************************

    double SigmaZvtxCut = 0.25;
    double DPhiCut = 0.01; //0.01
    double TklEtaCut  = 1;
    
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 12;

    double CentSPDTrackletsCentral = 1;
    double CentSPDTrackletsPeriph = 40;
    
   double MinDeltaEtaTKL = -2.4; //1.2  //2.4
    double MaxDeltaEtaTKL = 2.4;
    const int NbinsDeltaEtaTKL = 48; //24
    double SizeBinDeltaEtaTKL = (MaxDeltaEtaTKL-MinDeltaEtaTKL)/NbinsDeltaEtaTKL;
    
    double MinDeltaPhi = -TMath::Pi()/2;
    double MaxDeltaPhi = 1.5*TMath::Pi();
    const int NbinsDeltaPhiTKL = 48;
    double SizeBinDeltaPhiTKL = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhiTKL;
    int BinZeroLeftTKL = floor((0-MinDeltaPhi)*NbinsDeltaPhiTKL/(2*TMath::Pi()));
    
    const int NbBinsCent = 14;
    
    
    int DeltaEtaLongCMS = 2.0;
    int DeltaEtaShortCMS = 1.0;
    
    double R_SPD1 = 3.9;
    double R_SPD2 = 7.6;
    double corrFactorDeltaPhi = R_SPD1/(R_SPD2-R_SPD1);
    
    
    // Getting info from radical name
    
    
    Char_t DPhiCutText[20] = "10mrad";
    Char_t CentEstText[50] = "PercentileMethodSPDTracklets";
    Char_t CentralClassText[20] = "0-5";
    Char_t PeriphClassText[20] = "40-100";
    Char_t DeltaEtaGapText[20] = "1.2";
    Char_t ZvtxCutText[20] = "10";
    
    
    Char_t EMNormText[20] = "Method1";
    
    
    Char_t EMPoolMaxText[20] = "100";
    Char_t EMPoolThresholdText[20] = "10";
    Char_t ScalingPeriphText[20] = "1";
    
    
    Char_t SummationText[20] = "Method2";
    
    bool hasChangedDPhiCut = kFALSE;
    bool hasChangedEtaGap = kFALSE;
    bool hasChangedZvtxCut = kFALSE;
    bool hasChangedEMNorm = kFALSE;
    bool hasChangedEMMax = kFALSE;
    bool hasChangedEMThreshold = kFALSE;
    bool hasChangedPeriphScaling = kFALSE;
    bool hasChangedSummationZvtx = kFALSE;
    
    
    Char_t dataUsed[50];
    Char_t estimator[50];
    Char_t rest[50];
    Char_t speTampon[50];
    Char_t specificity[100];
    Char_t lessCentral[5];
    Char_t mostCentral[5];
    Char_t lessPeriph[5];
    Char_t mostPeriph[5];
    int lessCentralNum;
    int mostCentralNum;
    int lessPeriphNum;
    int mostPeriphNum;
    Char_t minpt[5];
    Char_t maxpt[5];
    
    double DPhiCutNum;
    double DeltaEtaGapNum;
    double ZvtxCutNum;
    
    sscanf(radical, "NewAnalysisAllEst_TKL_%[^_]_%[^_]_%[^-]-%[^_]_%[^-]-%[^_]_pt%[^-]-%s", dataUsed, estimator, mostCentral,lessCentral, lessPeriph, mostPeriph, minpt, rest);
    
    string srest;
    stringstream ss;
    ss << rest;
    ss >> srest;

    if (srest.find("_") != std::string::npos) {
        std::cout << "Underscore found: There is a specificity" << endl;
        sscanf(rest,"%[^_]_%s", maxpt, specificity);
    }
    else{
        cout << "No underscore found: No specificity"<<endl;
        sprintf(maxpt, "%s", rest);
    }
    
    string sspe;
    stringstream ss2;
    ss2 << specificity;
    ss2 >> sspe;
    
    while (sspe.find("_") != std::string::npos){
        std::cout << "Underscore found: There is another specificity" << endl;
        sscanf(specificity,"%[^_]_%s", speTampon, specificity);
        
        stringstream ss3;
        ss3 << specificity;
        ss3 >> sspe;
        
        string sspecT;
        stringstream ss3T;
        ss3T << speTampon;
        ss3T >> sspecT;
        
        if (sspecT.find("DPhiCut") != std::string::npos) {
            std::cout << "Change on DPhiCut" << endl;
            if (sspecT.find("DPhiCutNo") != std::string::npos) {
                cout << "There is no DPhi cut"<<endl;
                sprintf(DPhiCutText,"None");
                hasChangedDPhiCut=kTRUE;
                DPhiCutNum = 9999.;
            }
            else{
                sscanf(speTampon, "DPhiCut%[^mrad]", DPhiCutText);
                cout << "DPhiCutText set to " << DPhiCutText<< " mrad"<<endl;
                DPhiCutNum = stod(DPhiCutText);
                sprintf(DPhiCutText, "%smrad",DPhiCutText);
                hasChangedDPhiCut = kTRUE;
            }
        }
        
        if (sspecT.find("Eta") != std::string::npos) {
            std::cout << "Change on Eta" << endl;
           sscanf(speTampon, "Eta%s", DeltaEtaGapText);
           cout << "DeltaEtaGapText set to " << DeltaEtaGapText<<endl;
            hasChangedEtaGap = kTRUE;
        }
        
        if (sspecT.find("Zvtx") != std::string::npos) {
            std::cout << "Precisions on Zvtx" << endl;
           sscanf(speTampon, "Zvtx%s", ZvtxCutText);
           cout << "ZvtxCutText set to " << ZvtxCutText<<endl;
            hasChangedZvtxCut = kTRUE;
        }
        
        if (sspecT.find("EMNorm") != std::string::npos) {
            std::cout << "Precisions on EMNorm" << endl;
           sscanf(speTampon, "EMNorm%s", EMNormText);
           cout << "EMNormText set to " << EMNormText<<endl;
            hasChangedEMNorm = kTRUE;
        }
        
        if (sspecT.find("EMMax") != std::string::npos) {
            std::cout << "Precisions on EMMax" << endl;
           sscanf(speTampon, "EMMax%s", EMPoolMaxText);
           cout << "EMPoolMaxText set to " << EMPoolMaxText<<endl;
            hasChangedEMMax = kTRUE;
        }
        
        if (sspecT.find("EMThreshold") != std::string::npos) {
            std::cout << "Precisions on EMThreshold" << endl;
           sscanf(speTampon, "EMThreshold%s", EMPoolThresholdText);
           cout << "EMPoolThresholdText set to " << EMPoolThresholdText<<endl;
            hasChangedEMThreshold = kTRUE;
        }
        
        if (sspecT.find("PeriphScaling") != std::string::npos) {
            std::cout << "Precisions on PeriphScaling" << endl;
           sscanf(speTampon, "PeriphScaling%s", ScalingPeriphText);
           cout << "ScalingPeriphText set to " << ScalingPeriphText<<endl;
            hasChangedPeriphScaling = kTRUE;
        }
        
        if (sspecT.find("SummationZvtx") != std::string::npos) {
            std::cout << "Precisions on SummationZvtx" << endl;
           sscanf(speTampon, "SummationZvtx%s", SummationText);
           cout << "SummationText set to " << SummationText<<endl;
            hasChangedSummationZvtx= kTRUE;
        }
        
    }
    
    
    if (sspe.find("DPhiCut") != std::string::npos) {
        std::cout << "Change on DPhiCut" << endl;
        if (sspe.find("DPhiCutNo") != std::string::npos) {
            cout << "There is no DPhi cut"<<endl;
            sprintf(DPhiCutText,"None");
            hasChangedDPhiCut=kTRUE;
            DPhiCutNum = 9999.;
        }
        else{
            sscanf(specificity, "DPhiCut%[^mrad]", DPhiCutText);
            cout << "DPhiCutText set to " << DPhiCutText<< " mrad"<<endl;
            DPhiCutNum = stod(DPhiCutText);
            sprintf(DPhiCutText, "%smrad",DPhiCutText);
            hasChangedDPhiCut = kTRUE;
        }
    }
    
    if (sspe.find("Eta") != std::string::npos) {
        std::cout << "Change on Eta" << endl;
       sscanf(specificity, "Eta%s", DeltaEtaGapText);
       cout << "DeltaEtaGapText set to " << DeltaEtaGapText<<endl;
        hasChangedEtaGap = kTRUE;
    }
    
    if (sspe.find("Zvtx") != std::string::npos) {
        std::cout << "Precisions on Zvtx" << endl;
       sscanf(specificity, "Zvtx%s", ZvtxCutText);
       cout << "ZvtxCutText set to " << ZvtxCutText<<endl;
        hasChangedZvtxCut = kTRUE;
    }
    
    if (sspe.find("EMNorm") != std::string::npos) {
        std::cout << "Precisions on EMNorm" << endl;
       sscanf(specificity, "EMNorm%s", EMNormText);
       cout << "EMNormText set to " << EMNormText<<endl;
        hasChangedEMNorm = kTRUE;
    }
    
    if (sspe.find("EMMax") != std::string::npos) {
        std::cout << "Precisions on EMMax" << endl;
       sscanf(specificity, "EMMax%s", EMPoolMaxText);
       cout << "EMPoolMaxText set to " << EMPoolMaxText<<endl;
        hasChangedEMMax = kTRUE;
    }
    
    if (sspe.find("EMThreshold") != std::string::npos) {
        std::cout << "Precisions on EMThreshold" << endl;
       sscanf(specificity, "EMThreshold%s", EMPoolThresholdText);
       cout << "EMPoolThresholdText set to " << EMPoolThresholdText<<endl;
        hasChangedEMThreshold = kTRUE;
    }
    
    if (sspe.find("PeriphScaling") != std::string::npos) {
        std::cout << "Precisions on PeriphScaling" << endl;
       sscanf(specificity, "PeriphScaling%s", ScalingPeriphText);
       cout << "ScalingPeriphText set to " << ScalingPeriphText<<endl;
        hasChangedPeriphScaling = kTRUE;
    }
    
    if (sspe.find("SummationZvtx") != std::string::npos) {
        std::cout << "Precisions on SummationZvtx" << endl;
       sscanf(specificity, "SummationZvtx%s", SummationText);
       cout << "SummationText set to " << SummationText<<endl;
        hasChangedSummationZvtx= kTRUE;
    }
    
    stringstream intValue1(lessCentral);
    intValue1 >> lessCentralNum;
    stringstream intValue2(mostCentral);
    intValue2 >> mostCentralNum;
    stringstream intValue3(lessPeriph);
    intValue3 >> lessPeriphNum;
    stringstream intValue4(mostPeriph);
    intValue4 >> mostPeriphNum;
    
    sprintf(CentEstText,"%s", estimator);
    cout << "CentEstText " << CentEstText<<endl;
    sprintf(CentralClassText, "%i-%i", mostCentralNum, lessCentralNum);
    cout << "CentralClassText " << CentralClassText<<endl;
    sprintf(PeriphClassText, "%i-%i", lessPeriphNum, mostPeriphNum);
    cout << "PeriphClassText " << PeriphClassText<<endl;
    cout << "min pt " << minpt << " max pt " << maxpt<<endl;
    
    // Là on a toutes les linformations de radical et on sait qui a été changé
    
    //MAJ des informations sur la centralité (estimateur et classes)
    
    CentralLowBound = GetCent(mostCentralNum);
    CentralHighBound = GetCent(lessCentralNum);
    PeripheralLowBound = GetCent(lessCentralNum)+1;
    PeripheralHighBound = GetCent(mostPeriphNum);
    
    stringstream scent;
    scent << CentEstText;
    scent >> centralityMethod;
    
    
    double ZvtxCut = stod(ZvtxCutText);
    if(hasChangedDPhiCut){
        DPhiCut = DPhiCutNum;
    }

    double DeltaEtaTKLCut = stod(DeltaEtaGapText); //1.2  //1.2
     
    int NbBinsZvtx = int(2*ZvtxCut);
    
    int EMPoolMax = stoi(EMPoolMaxText);
    int EMPoolThreshold = stoi(EMPoolThresholdText);
    double ScalingPeriph = stod(ScalingPeriphText);
    
    //FIXME__: Adapter les déclarations avec NBBinsZvtx qui n'est plus const
    //FIXME__: Ajouter les valeurs Num des 5 dernières variables systematiques
    //FIXME__: Adapter arrayofPeriods
    //FIXME__: Adapter % =!
    
    
    
    
    
    
    
    Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
    string sdataUsed;
    stringstream strim;
    strim << dataUsed;
    strim >> sdataUsed;
    
    bool is16h25 = kFALSE;
    int dataPercentage = 100;
    
    if(strncmp (dataUsed,"Run2",10) == 0){
//        arrayOfPeriods = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l","Group2_LHC17h","Group3_LHC17h","Group4_LHC17k","Group4_LHC18l","Group4_LHC18m","Group4_LHC18o","Group4_LHC18p","Group5_LHC17l","Group5_LHC17m","Group5_LHC17o","Group5_LHC17r","Group5_LHC18c","Group5_LHC18d","Group5_LHC18e","Group5_LHC18f","Group6_LHC18m","Group7_LHC18m","Group8_LHC18m","Group9_LHC18m","Group10_LHC18m","Group11_LHC18m","Group12_LHC18m"};
    }
    else if(sdataUsed.find("16hj") != std::string::npos){
//        arrayOfPeriods = {"Group1_LHC16h","Group1_LHC16j"};
        if(strncmp (dataUsed,"16hj",10) != 0){
            sscanf(dataUsed, "16hj%i", &dataPercentage);
        }
    }
    else if(sdataUsed.find("16h") != std::string::npos){
//        arrayOfPeriods = {"Group1_LHC16h"};
        if(strncmp (dataUsed,"16h25",10) == 0){
            is16h25 = kTRUE;
        }
        else if(strncmp (dataUsed,"16h",10) != 0){
            sscanf(dataUsed, "16h%i", &dataPercentage);
        }
    }
    else{
        cout << "The data set you are asking for is not a possibility yet. Put it by hand"<<endl;
        return;
    }
  //  Char_t Group_Period[50] = "Group1";
  // Char_t *arrayOfPeriods[] = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l","Group2_LHC17h","Group3_LHC17h","Group4_LHC17k","Group4_LHC18l","Group4_LHC18m","Group4_LHC18o","Group4_LHC18p","Group5_LHC17l","Group5_LHC17m","Group5_LHC17o","Group5_LHC17r","Group5_LHC18c","Group5_LHC18d","Group5_LHC18e","Group5_LHC18f","Group6_LHC18m","Group7_LHC18m","Group8_LHC18m","Group9_LHC18m","Group10_LHC18m","Group11_LHC18m","Group12_LHC18m"};
  //  Char_t *arrayOfPeriods[] = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l"};
   // Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
   // Char_t *arrayOfPeriods[] = {"Cvetan_LHC16r"};
  //  Char_t *arrayOfPeriods[] = {"Group1_LHC16h","Group1_LHC16j"};
  // Char_t *arrayOfPeriods[] = {"Group8_LHC18m_CvetanPU_OnlyMuonTrackCutsApplied"};
    int numberOfPeriods = sizeof(arrayOfPeriods) / sizeof(arrayOfPeriods[0]);
    
    Char_t fileInLoc[200];
    Char_t FolderName[200];
    Char_t CanvasName[200];
    Char_t FitFileName[200];
 //   Char_t AssociateFileName[200];
    
    Char_t RadicalName[200];
    sprintf(RadicalName,"%s",radical);

    sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_%s.root",RadicalName);
//    sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/Systematics/FitFile_TKL_DPhi10mrad_PercentileMethodSPDTracklets_EtaGap1.2_Zcut10.root");
//    sprintf(AssociateFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/Systematics/AssociateFile_TKL_DPhi10mrad_PercentileMethodSPDTracklets_EtaGap1.2_Zcut10.txt");
    sprintf(FolderName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2021octobre/%s",RadicalName);
    
    int rc = mkdir(FolderName, 0777);
    cout << "rc " <<rc<<endl;
    
    

// *************************
// Initialiser les graphes *
// *************************
    
    TH1F* hnseg;
    TH1F* hnseg2(NULL);
    TH1F* hnseg3(NULL);
    TH1F* hnseg4(NULL);
    TH1F* hnseg5(NULL);
    TH1F* hnseg6(NULL);
    TH1F* hnseg7(NULL);
    TH2F* hnseg8(NULL);
    TH2F* hnsegSigma(NULL);
    TH1F* hnseg8Tkl_proj_tampon(NULL);
    TH1F* hnseg8Tkl_proj(NULL);
    TH1F* hnseg8TklME_proj_tampon(NULL);
    TH1D* hnseg8TklME_proj(NULL);
    TH1F* ProjCopyTkl(NULL);
    TH1F* ProjCopy2Tkl(NULL);
    TH1F* ME_proj_Tkl_tampon(NULL);
    TH1I* SE_proj_Tkl_tampon(NULL);
    
    TH1F* Sijk(NULL);
    TH1F* Mijk(NULL);
    TH1F* SoverMijk(NULL);
    TH1F* SoverMik(NULL);
    
    TH1F* Sij(NULL);
    TH1F* Mij(NULL);
    TH1F* SoverMij(NULL);
    TH1F* SoverMi(NULL);
    
    TH1F* YieldsTkl[NbBinsCent]{ NULL };
    TH1F* YieldsTklShortCMS[NbBinsCent]{ NULL };
    TH1F* YieldsTklLongCMS[NbBinsCent]{ NULL };
    TH2I* CorrelationsTkl[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTkl, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2I*));
    TH2I* CorrelationsTklShortCMS[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklShortCMS, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2I*));
    TH2I* CorrelationsTklLongCMS[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklLongCMS, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2I*));
    TH2I* CorrelationsTklME[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklME, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2I*));
    TH2I* CorrelationsTklMEShortCMS[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklMEShortCMS, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2I*));
    TH2I* CorrelationsTklMELongCMS[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklMELongCMS, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2I*));
    TH2F* CorrelationsTklMEScaled[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklMEScaled, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2F*));
    TH2F* CorrelationsTklMEScaledShortCMS[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklMEScaledShortCMS, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2F*));
    TH2F* CorrelationsTklMEScaledLongCMS[NbBinsCent][NbBinsZvtx];
    memset( CorrelationsTklMEScaledLongCMS, 0, NbBinsCent*NbBinsZvtx*sizeof(TH2F*));
    TH1F* Yield_tampon(NULL);
    TH1F* YieldTkl_tampon(NULL);
    TH1F* Yield_allC(NULL);
    TH1F* YieldTkl_allC(NULL);
    TH1F* YieldTkl_Central(NULL);
    TH1F* YieldTkl_Periph(NULL);
    TH1F* YieldTkl_Difference(NULL);
    TH1F* YieldTkl_allCShortCMS(NULL);
    TH1F* YieldTkl_CentralShortCMS(NULL);
    TH1F* YieldTkl_PeriphShortCMS(NULL);
    TH1F* YieldTkl_DifferenceShortCMS(NULL);
    TH1F* YieldTkl_allCLongCMS(NULL);
    TH1F* YieldTkl_CentralLongCMS(NULL);
    TH1F* YieldTkl_PeriphLongCMS(NULL);
    TH1F* YieldTkl_DifferenceLongCMS(NULL);

    
    // Try plots 2D DetaDphi
    
    TH2F* YCentral(NULL);
    TH2F* YPeriph(NULL);
    TH2F* YDifference(NULL);
    TH2F* YTklCentral(NULL);
    TH2I* YTklCentralME(NULL);
    TH2F* YTklPeriph(NULL);
    TH2I* YTklPeriphME(NULL);
    TH2F* YTklDifference(NULL);
    TH1D* YTklCentral_proj_tampon(NULL);
    TH1D* Correl_tampon(NULL);
    TH1D* YTklPeriph_proj_tampon(NULL);
    TH1D* YTklDifference_proj_tampon(NULL);
    TH1D* YTklCentralME_proj_tampon(NULL);
    TH1D* YTklPeriphME_proj_tampon(NULL);
    
   
    YTklCentral = new TH2F("YTklCentral",
                      "Yield #Delta#eta wrt #Delta#phi - Tkl Central",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklCentral->SetXTitle("Correlation #Delta#phi (rad)");
    YTklCentral->SetYTitle("Correlation #Delta#eta");
    YTklPeriph = new TH2F("YTklPeriph",
                      "Yield #Delta#eta wrt #Delta#phi - Tkl Periph",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklPeriph->SetXTitle("Correlation #Delta#phi (rad)");
    YTklPeriph->SetYTitle("Correlation #Delta#eta");
    YTklDifference = new TH2F("YTklDifference",
                      "Yield #Delta#eta wrt #Delta#phi - Tkl Difference",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklDifference->SetXTitle("Correlation #Delta#phi (rad)");
    YTklDifference->SetYTitle("Correlation #Delta#eta");
    YTklCentralME = new TH2I("YTklCentralME",
                      "Yield #Delta#eta wrt #Delta#phi - Tkl Central ME",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklCentralME->SetXTitle("Correlation #Delta#phi (rad)");
    YTklCentralME->SetYTitle("Correlation #Delta#eta");
    YTklPeriphME = new TH2I("YTklPeriphME",
                      "Yield #Delta#eta wrt #Delta#phi - Tkl Periph ME",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklPeriphME->SetXTitle("Correlation #Delta#phi (rad)");
    YTklPeriphME->SetYTitle("Correlation #Delta#eta");
    // AE
    
    TH2F* hnseg9(NULL);
    TH2F* hnseg10(NULL);
    TH1F* hDPhi(NULL);
   
    
    
    ProjCopyTkl = new TH1F("ProjCopyTkl",
                      "ProjCopyTkl",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    ProjCopyTkl->SetXTitle("Correlation #Delta#phi (rad)");
    ProjCopyTkl->SetYTitle("Correlation #Delta#eta");
    ProjCopy2Tkl = new TH1F("ProjCopy2Tkl",
                      "ProjCopy2Tkl",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    ProjCopy2Tkl->SetXTitle("Correlation #Delta#phi (rad)");
    ProjCopy2Tkl->SetYTitle("Correlation #Delta#eta");
    
    char hname[200];
    char hname1[200];
    char hname2[200];
    
    
    for(int i=0; i<NbBinsCent; i++){
        char hname[200];
        sprintf(hname,"Yields Tkl %d ",i);
        YieldsTkl[i] = new TH1F(hname,
                          "Yields Correlation Tkl wrt #Delta#phi",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldsTkl[i]->SetXTitle("Correlation #Delta#phi (rad)");
        sprintf(hname,"Yields Tkl ShortCMS %d ",i);
        YieldsTklShortCMS[i] = new TH1F(hname,
                          "Yields Correlation Tkl ShortCMS wrt #Delta#phi",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldsTklShortCMS[i]->SetXTitle("Correlation #Delta#phi (rad)");
        sprintf(hname,"Yields Tkl LongCMS %d ",i);
        YieldsTklLongCMS[i] = new TH1F(hname,
                          "Yields Correlation Tkl LongCMS wrt #Delta#phi",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldsTklLongCMS[i]->SetXTitle("Correlation #Delta#phi (rad)");
        for(int j=0; j<NbBinsZvtx; j++){
            sprintf(hname,"Correlations Tkl %d %d ",i,j);
            CorrelationsTkl[i][j] = new TH2I(hname,
                              "Correlation Tkl #Delta#eta wrt #Delta#phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTkl[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTkl[i][j]->SetYTitle("Correlation #Delta#eta");
            sprintf(hname,"CorrelationsME Tkl %d %d ",i,j);
            CorrelationsTklME[i][j] = new TH2I(hname,
                              "MixedEvent Correlation Tkl #Delta#eta wrt #Delta#phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklME[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklME[i][j]->SetYTitle("Correlation #Delta#eta");
            sprintf(hname,"CorrelationsMEScaled Tkl %d %d ",i,j);
            CorrelationsTklMEScaled[i][j] = new TH2F(hname,
                              "MixedEvent Correlation Tkl #Delta#eta wrt #Delta#phi - Scaled",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklMEScaled[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklMEScaled[i][j]->SetYTitle("Correlation #Delta#eta");
            
            sprintf(hname,"Correlations Tkl ShortCMS %d %d ",i,j);
            CorrelationsTklShortCMS[i][j] = new TH2I(hname,
                              "Correlation Tkl ShortCMS #Delta#eta wrt #Delta#phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklShortCMS[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklShortCMS[i][j]->SetYTitle("Correlation #Delta#eta");
            sprintf(hname,"CorrelationsME Tkl ShortCMS %d %d ",i,j);
            CorrelationsTklMEShortCMS[i][j] = new TH2I(hname,
                              "MixedEvent Correlation Tkl ShortCMS #Delta#eta wrt #Delta#phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklMEShortCMS[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklMEShortCMS[i][j]->SetYTitle("Correlation #Delta#eta");
            sprintf(hname,"CorrelationsME Tkl ShortCMS Scaled %d %d ",i,j);
            CorrelationsTklMEScaledShortCMS[i][j] = new TH2F(hname,
                              "MixedEvent Correlation Tkl ShortCMS #Delta#eta wrt #Delta#phi - Scaled",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklMEScaledShortCMS[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklMEScaledShortCMS[i][j]->SetYTitle("Correlation #Delta#eta");
            
            sprintf(hname,"Correlations Tkl LongCMS %d %d ",i,j);
            CorrelationsTklLongCMS[i][j] = new TH2I(hname,
                              "Correlation Tkl LongCMS #Delta#eta wrt #Delta#phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklLongCMS[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklLongCMS[i][j]->SetYTitle("Correlation #Delta#eta");
            sprintf(hname,"CorrelationsME Tkl LongCMS %d %d ",i,j);
            CorrelationsTklMELongCMS[i][j] = new TH2I(hname,
                              "MixedEvent Correlation Tkl LongCMS #Delta#eta wrt #Delta#phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklMELongCMS[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklMELongCMS[i][j]->SetYTitle("Correlation #Delta#eta");
            sprintf(hname,"CorrelationsME Tkl LongCMS Scaled %d %d ",i,j);
            CorrelationsTklMEScaledLongCMS[i][j] = new TH2F(hname,
                              "MixedEvent Correlation Tkl LongCMS #Delta#eta wrt #Delta#phi - Scaled",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklMEScaledLongCMS[i][j]->SetXTitle("Correlation #Delta#phi (rad)");
            CorrelationsTklMEScaledLongCMS[i][j]->SetYTitle("Correlation #Delta#eta");
        }
    }
    
     YieldTkl_allC = new TH1F("YieldTkl_allC",
                            "Yields tkl-tkl wrt #Delta#phi, all C, all mass",
                            NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
          YieldTkl_allC->SetXTitle("#Delta#phi (rad)");
       YieldTkl_allC->SetXTitle("Yields");
       YieldTkl_Central = new TH1F("YieldTkl_Central",
                         "Yields tkl-tkl wrt #Delta#phi, Central",
                         NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
       YieldTkl_Central->SetXTitle("#Delta#phi (rad)");
       YieldTkl_Central->SetYTitle("Yields_{Central}");
       YieldTkl_Periph = new TH1F("YieldTkl_Periph",
                         "Yields tkl-tkl wrt #Delta#phi, Periph",
                         NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
       YieldTkl_Periph->SetXTitle("#Delta#phi (rad)");
       YieldTkl_Periph->SetYTitle("Yields_{Periph}");
       YieldTkl_Difference = new TH1F("YieldTkl_Difference",
                         "Yields tkl-tkl wrt #Delta#phi, Central-Periph, all mass",
                         NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
       YieldTkl_Difference->SetXTitle("#Delta#phi (rad)");
       YieldTkl_Difference->SetYTitle("Yields_{Subtracted}");
       YieldTkl_tampon = new TH1F("YieldTkl_tampon",
                         "Yields tkl-tkl wrt #Delta#phi, all C, all mass",
                         NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
       YieldTkl_tampon->SetXTitle("#Delta#phi (rad)");
       YieldTkl_tampon->SetYTitle("Yields");
    
    YieldTkl_allCShortCMS = new TH1F("YieldTkl_allCShortCMS",
                         "Yields Correlation Tkl ShortCMS wrt #Delta#phi, all C, all mass",
                         NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
       YieldTkl_allCShortCMS->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_CentralShortCMS = new TH1F("YieldTkl_CentralShortCMS",
                      "Yields Correlation Tkl ShortCMS wrt #Delta#phi, Central",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_CentralShortCMS->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_PeriphShortCMS = new TH1F("YieldTkl_PeriphShortCMS",
                      "Yields Correlation Tkl ShortCMS wrt #Delta#phi, Periph",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_PeriphShortCMS->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_DifferenceShortCMS = new TH1F("YieldTkl_DifferenceShortCMS",
                      "Yields Correlation Tkl ShortCMS wrt #Delta#phi, Central-Periph, all mass",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_DifferenceShortCMS->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_allCLongCMS = new TH1F("YieldTkl_allCLongCMS",
                         "Yields Correlation Tkl LongCMS wrt #Delta#phi, all C, all mass",
                         NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
       YieldTkl_allCLongCMS->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_CentralLongCMS = new TH1F("YieldTkl_CentralLongCMS",
                      "Yields Correlation Tkl LongCMS wrt #Delta#phi, Central",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_CentralLongCMS->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_PeriphLongCMS = new TH1F("YieldTkl_PeriphLongCMS",
                      "Yields Correlation Tkl LongCMS wrt #Delta#phi, Periph",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_PeriphLongCMS->SetXTitle("Correlation #Delta#phi (rad)");
    YieldTkl_DifferenceLongCMS = new TH1F("YieldTkl_DifferenceLongCMS",
                      "Yields Correlation Tkl LongCMS wrt #Delta#phi, Central-Periph, all mass",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_DifferenceLongCMS->SetXTitle("Correlation #Delta#phi (rad)");
    
    
    hDPhi = new TH1F("hDPhi",
                     "#Delta#Phi of Tracklet",
                     10000,-TMath::Pi(),TMath::Pi());
    hDPhi->SetXTitle("Tracklet #Delta#Phi (rad)");
    hDPhi->SetYTitle("Count");
    
    hnseg8Tkl_proj = new TH1F("hnseg8Tkl_proj", "Correlation tkl-tkl wrt #Delta#eta, projection",NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8Tkl_proj->SetXTitle("Correlation #Delta#eta");
    hnseg8Tkl_proj->SetYTitle("Correlation mean");
    hnseg8TklME_proj = new TH1D("hnseg8TklME_proj", "MixedEvent Correlation tkl-tkl wrt #Delta#eta, projection",NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8TklME_proj->SetXTitle("Correlation #Delta#eta");
    hnseg8TklME_proj->SetYTitle("Correlation mean");
    
    Sij = new TH1F("Sij",
                                   "Sij",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    Sij->SetXTitle("Correlation #Delta#phi (rad)");
    Mij = new TH1F("Mij",
                                   "Mij",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    Mij->SetXTitle("Correlation #Delta#phi (rad)");
    SoverMij = new TH1F("SoverMij",
                                   "SoverMij",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    SoverMij->SetXTitle("Correlation #Delta#phi (rad)");
    SoverMi = new TH1F("SoverMi",
                                   "SoverMi",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    SoverMi->SetXTitle("Correlation #Delta#phi (rad)");

    hnseg9 = new TH2F("hnseg9",
                      "Tracklet #eta wrt z_{vtx}",
                      100,-15,15,100,-2.5,2.5);
    hnseg9->SetXTitle("Z_{vtx} (cm)");
    hnseg9->SetYTitle("Tracklet #eta");
    hnseg10 = new TH2F("hnseg10",
                      "Number of tracklets wrt V0M percentile",
                      100,0,100,100,0,400);
    hnseg10->SetXTitle("V0M Percentile");
    hnseg10->SetYTitle("Tracklet count");
    
    
    
    
    
// *************************
// Analyse                 *
// *************************
    
    
   
    int EventRejected = 0;
    int EventNC = 0;
    int EventPileUpMult = 0;
    int EventPileUpVtx = 0;
    int RefTracklets = 0;
    int RefTrackletsCentral = 0;
    int RefTrackletsPeriph = 0;
    int cmul = 0;
    int barWidth = 50;
    int RefTklCounter[NbBinsCent][NbBinsZvtx];
    memset( RefTklCounter, 0, NbBinsCent*NbBinsZvtx*sizeof(int));
    int RefTklCounterZint[NbBinsCent] = {0};
    double NormMETkl[NbBinsCent][NbBinsZvtx];
    memset( NormMETkl, 0, NbBinsCent*NbBinsZvtx*sizeof(double));
    double NormMETklZint[NbBinsCent] = {0};
    
    //Nassoc
    int Nassoc_Central_TKL = 0;
    int Nassoc_Periph_TKL = 0;
    
    
    
    int TklC = 0;
    int TklP = 0;
    int NormTklCentral = 0;
    int NormTklPeriph = 0;
    int countsigma = 0;
    
    
    MyEventLight* fEvent = 0;
    TClonesArray *fCorrelations = 0;
    TClonesArray *fTracklets = 0;
    TClonesArray *fDimuons = 0;
    CorrelationLight *correl = 0; //= new CorrelationLight(); //object must be created before
    TrackletLight *trac = 0;
    TrackletLight *tracklet1 = 0;
    TrackletLight *tracklet2 = 0;
    TrackletLight *trackletME = 0;
    DimuonLight *dimu = 0;
    TTree* theTree = NULL;
    
    
// ************************************
// Boucle sur les evenements, notés i *
// ************************************
    
    cout << "=== NOW EVENT STUDY ===" <<endl;
    cout << "Switching back to 1st TTree" <<endl;
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
   //     sprintf(fileInLoc,"~/../../Volumes/Sauvegarde /LegacySebAnalysepp/NewAnalysis/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CMUL/%s_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
     //   sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/NewAnalysis_AllEst/CINT/%s_CINT_AllEst/muonGrid.root",arrayOfPeriods[tree_idx]);
    TFile fileIn(fileInLoc);
        
        char str [50];
        int GroupNum;
        
        std::vector <double> PoolsTkl[NbBinsCent][NbBinsZvtx];
        std::vector <int> PoolsTklEventTracker[NbBinsCent][NbBinsZvtx];
        int PoolsSizeTkl[NbBinsCent][NbBinsZvtx];
        memset( PoolsSizeTkl, 0, NbBinsCent*NbBinsZvtx*sizeof(int));

        sscanf (arrayOfPeriods[tree_idx],"Group%d_%s",&GroupNum,str);
    
    std::cout << "1" <<std::endl;
            fileIn.GetObject("MyMuonTree",theTree);
        std::cout << "2" <<std::endl;
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
        fTracklets = fEvent->fTracklets;
      //  fDimuons = fEvent->fDimuons;
    
    branch = theTree->GetBranch("event");
    branch->SetAddress(&fEvent);
    fTracklets = fEvent->fTracklets;

        int TklMEcounter =0;
        
        ofstream myfiletxt;
        myfiletxt.open("/tmp/memorycheckbin.txt");
    
        for (int i=0;i<nevent;i++) {
            if(i%10000 == 0){
                    std::cout << "[";
                    double portion = double(i)/nevent;
                    long pos = barWidth * portion;
                    for (int k = 0; k < barWidth; ++k) {
                        if (k < pos) std::cout << "=";
                        else if (k == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                    std::cout << "] " << long(100 * portion) << "%     " << i << "/" << nevent << " Tree " << tree_idx+1 << "/" << numberOfPeriods << " " <<radical;
                    std::cout.flush();
                std::cout << std::endl;
            }
            
            if(is16h25){
                if(doTracklets && (i%4!=0)){
                    continue;
                }
            }
            
            if(doTracklets && (i%100 >= dataPercentage)){
                continue;
            }
                theTree->GetEvent(i);
            if(fEvent->fIsPileupFromSPDMultBins){
                EventPileUpMult++;
            }
//            if(fEvent->fNPileupVtx != 0){
//                EventPileUpVtx++;
//            }
            if(fEvent->fVertexNC < 1){
                EventNC++;
            }
            if(fEvent->fIsPileupFromSPDMultBins || fEvent->fVertexNC < 1){ //|| fEvent->fNPileupVtx != 0
                EventRejected++;
                continue;
            }
            int NumberCloseEtaTracklets = 0;
            int NumberOfTrackletsPassingEtaCut = 0;
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                if((TMath::Abs(trac->fEta) < TklEtaCut) && (TMath::Abs(trac->fDPhi) < DPhiCut)){
                    NumberOfTrackletsPassingEtaCut++;
                    hDPhi->Fill(trac->fDPhi);
                    if((TMath::Abs(trac->fDPhi) < DPhiCut)){
                        NumberCloseEtaTracklets++;
                    }
                }
                
            }
            
            if(NumberOfTrackletsPassingEtaCut==0){
                continue;
            }
            
            if(KeepOnlyOne){
                if(valueOnlyOne!=NumberOfTrackletsPassingEtaCut){
                    continue;
                }
            }
            
            if(AdditionalCutNtkl){
                if(NumberOfTrackletsPassingEtaCut <= valueAdditionalCutNtkl){
                    continue;
                }
            }
            
            if((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut)){
                hnseg10->Fill(fEvent->fCentralityV0M, NumberCloseEtaTracklets); //fEvent->fNTracklets
            }
            
//            bool isEventDimu = kFALSE;
//
//
//            for (Int_t j=0; j<fEvent->fNDimuons; j++) {
//            dimu = (DimuonLight*)fDimuons->At(j);
//            if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut) && (dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)){
//                isEventDimu = kTRUE;
//            }
//                if(!isEventDimu){
//                    continue;
//                }
//            }
            
// TREATEMENT OF TRACKLET V2
            if(doTracklets){
                
                if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut)){
                 // double cent = fEvent->fCentralitySPDTracklets;
                //  int centint = GetCent(cent);
                  double zv = fEvent->fVertexZ;
                  int zvint = floor(zv) + ZvtxCut;
              //  cout << NumberOfTrackletsPassingEtaCut << " " << zvint << " " <<GroupNum<<endl;
                  int centint = GetCentPM(NumberOfTrackletsPassingEtaCut, fEvent->fCentralitySPDTracklets, fEvent->fCentralitySPDClusters, fEvent->fCentralityV0M, zvint, GroupNum);
                  //  int centint = GetCentCvetan(fEvent->fCentralityV0M);
                    if(NumberCloseEtaTracklets>=2){
                   RefTklCounter[centint][zvint] += NumberCloseEtaTracklets-1; //fEvent->fNTracklets //FIXME ok
                    fTracklets->Randomize(); //Moved here to avoid randomising everytime
                    }
              //  if(cent <= CentSPDTrackletsCentral){
                if(isCentral(centint)){
                    TklC++;
                }
               // else if(cent > CentSPDTrackletsPeriph){
                if(isPeripheral(centint)){
                    TklP++;
                }
                }
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                tracklet1 = (TrackletLight*)fTracklets->At(j);
                for (Int_t k=j+1; k<fEvent->fNTracklets; k++){
                    tracklet2 = (TrackletLight*)fTracklets->At(k);
                    if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut) && (TMath::Abs(tracklet1->fDPhi) < DPhiCut) && (TMath::Abs(tracklet2->fDPhi) < DPhiCut) && (tracklet1->fEta < myEtaHigh(fEvent->fVertexZ)) && (tracklet1->fEta > myEtaLow(fEvent->fVertexZ)) && (tracklet2->fEta < myEtaHigh(fEvent->fVertexZ)) && (tracklet2->fEta > myEtaLow(fEvent->fVertexZ))){   //Cuts FIXME ok
                        Float_t DeltaPhi = (tracklet1->fPhi + corrFactorDeltaPhi*(tracklet1->fDPhi)) - (tracklet2->fPhi + corrFactorDeltaPhi*(tracklet2->fDPhi));
                            if(DeltaPhi < -TMath::Pi()/2){
                                DeltaPhi += 2* TMath::Pi();
                            }
                            if(DeltaPhi > 1.5*TMath::Pi()){
                                DeltaPhi -= 2* TMath::Pi();
                            }
                            Float_t DeltaEta = tracklet1->fEta - tracklet2->fEta; //DeltaEtaAbs TMath::Abs(
                        double zv = fEvent->fVertexZ;
                        int zvint = floor(zv) + ZvtxCut;
                         int centint = GetCentPM(NumberOfTrackletsPassingEtaCut, fEvent->fCentralitySPDTracklets, fEvent->fCentralitySPDClusters, fEvent->fCentralityV0M, zvint, GroupNum);
                      //  int centint = GetCentCvetan(fEvent->fCentralityV0M);
                        
                        if(isCMSMethod){
                            if(TMath::Abs(DeltaEta)<DeltaEtaShortCMS){
                                 CorrelationsTklShortCMS[centint][zvint]->Fill(DeltaPhi,DeltaEta);
                             }
                            if(TMath::Abs(DeltaEta)>DeltaEtaLongCMS && TMath::Abs(DeltaEta)<MaxDeltaEtaTKL){
                                CorrelationsTklLongCMS[centint][zvint]->Fill(DeltaPhi,DeltaEta);
                            }
                        }
                        if(isCentral(centint)){
                              Nassoc_Central_TKL++;
                          }
                          if(isPeripheral(centint)){
                              Nassoc_Periph_TKL++;
                          }
                        
                        if((TMath::Abs(DeltaEta)<DeltaEtaTKLCut) || (TMath::Abs(DeltaEta)>MaxDeltaEtaTKL)){
                            continue;
                        }
                            //   double cent = fEvent->fCentralitySPDTracklets;
                          //     int centint = GetCent(cent);
                               CorrelationsTkl[centint][zvint]->Fill(DeltaPhi,DeltaEta);
                              //  if(cent <= CentSPDTrackletsCentral){
                                if(isCentral(centint)){
                                    YTklCentral->Fill(DeltaPhi,DeltaEta);
                                }
                              //  else if(cent > CentSPDTrackletsPeriph){
                                if(isPeripheral(centint)){
                                    YTklPeriph->Fill(DeltaPhi,DeltaEta);
                                }
                   }
                }
            }
                
    // Event mixing
                
                if(doMixedEvents){
                    if((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut)){
                        //   double centME = fEvent->fCentralitySPDTracklets;
                        //  int centintME = GetCent(centME);
                          double zvME = fEvent->fVertexZ;
                          int zvintME = floor(zvME) + ZvtxCut;
                            int centintME = GetCentPM(NumberOfTrackletsPassingEtaCut, fEvent->fCentralitySPDTracklets, fEvent->fCentralitySPDClusters, fEvent->fCentralityV0M, zvintME, GroupNum);
                     //   int centintME = GetCentCvetan(fEvent->fCentralityV0M);
                        TklMEcounter++;
                        
                        if(PoolsSizeTkl[centintME][zvintME]>=10){
                            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                              tracklet1 = (TrackletLight*)fTracklets->At(j);
                                if((tracklet1->fEta < myEtaHigh(zvME)) && (tracklet1->fEta > myEtaLow(zvME)) && (TMath::Abs(tracklet1->fDPhi) < DPhiCut)){ //FIXME ok

                                  for(int k=0; k< PoolsTkl[centintME][zvintME].size(); k+=2){
                                              double correlTklMEPhi = (tracklet1->fPhi + corrFactorDeltaPhi*(tracklet1->fDPhi)) - PoolsTkl[centintME][zvintME].at(k);
                                              if(correlTklMEPhi < -TMath::Pi()/2){
                                                  correlTklMEPhi += 2* TMath::Pi();
                                              }
                                              if(correlTklMEPhi > 1.5*TMath::Pi()){
                                                  correlTklMEPhi -= 2* TMath::Pi();
                                              }
                                      double correlTklMEEta = tracklet1->fEta - PoolsTkl[centintME][zvintME].at(k+1); //DeltaEtaAbs
                                          if(TMath::Abs(correlTklMEEta)>DeltaEtaTKLCut && TMath::Abs(correlTklMEEta)<MaxDeltaEtaTKL){
                                               CorrelationsTklME[centintME][zvintME]->Fill(correlTklMEPhi,correlTklMEEta);
//                                              myfiletxt << "CorrelationsTklME[" <<centintME<<"]["<<zvintME<<"] has entries " << CorrelationsTklME[centintME][zvintME]->GetEntries()<<endl;
//                                              myfiletxt << "Updated bin has entries" << CorrelationsTklME[centintME][zvintME]->GetBinContent(CorrelationsTklME[centintME][zvintME]->GetXaxis()->FindBin(correlTklMEPhi),CorrelationsTklME[centintME][zvintME]->GetYaxis()->FindBin(correlTklMEEta))<<endl;
                                           }
                                      
                                      if(isCMSMethod){
                                          if(TMath::Abs(correlTklMEEta)<DeltaEtaShortCMS){
                                                   CorrelationsTklMEShortCMS[centintME][zvintME]->Fill(correlTklMEPhi,correlTklMEEta);
                                               }
                                          if(TMath::Abs(correlTklMEEta)>DeltaEtaLongCMS){
                                              CorrelationsTklMELongCMS[centintME][zvintME]->Fill(correlTklMEPhi,correlTklMEEta);
                                          }
                                      }
                                      
                                            if(correlTklMEPhi > -SizeBinDeltaPhiTKL/2 && correlTklMEPhi < SizeBinDeltaPhiTKL/2 && TMath::Abs(correlTklMEEta) < SizeBinDeltaEtaTKL/2){
                                                NormMETkl[centintME][zvintME] += 1.;
                                                NormMETklZint[centintME] += 1.;
                                            }
                                      
                                      
                                  //        if(cent <= CentSPDTrackletsCentral){
                                      if(isCentral(centintME)){
                                              if(TMath::Abs(correlTklMEEta)>DeltaEtaTKLCut && TMath::Abs(correlTklMEEta)<MaxDeltaEtaTKL){
                                                  YTklCentralME->Fill(correlTklMEPhi,correlTklMEEta);
                                              }
                                              if(correlTklMEPhi > 0 && correlTklMEPhi < SizeBinDeltaPhiTKL && TMath::Abs(correlTklMEEta) < SizeBinDeltaEtaTKL){
                                                  NormTklCentral += 1;
                                              }
                                          }
                                 //         else if(cent > CentSPDTrackletsPeriph){
                                      if(isPeripheral(centintME)){
                                             if(TMath::Abs(correlTklMEEta)>DeltaEtaTKLCut && TMath::Abs(correlTklMEEta)<MaxDeltaEtaTKL){
                                                  YTklPeriphME->Fill(correlTklMEPhi,correlTklMEEta);
                                              }
                                              if(correlTklMEPhi > 0 && correlTklMEPhi < SizeBinDeltaPhiTKL && TMath::Abs(correlTklMEEta) < SizeBinDeltaEtaTKL){
                                                  NormTklPeriph += 1;
                                              }
                                          }
                                          
                                  }
                                }
                            }
                        }
                        
                        if(PoolsSizeTkl[centintME][zvintME] == 100){
                         //    cout << "The Tkl pool is full" <<endl;
                            int valueDiscarded = PoolsTklEventTracker[centintME][zvintME].front();
                           //  cout << "Discarding events with index " << valueDiscarded<<endl;
                            while(PoolsTklEventTracker[centintME][zvintME].front()==valueDiscarded){
                          //      cout << "Event number at front: " << PoolsTklEventTracker[centintME].front()<<endl;
                                PoolsTklEventTracker[centintME][zvintME].erase(PoolsTklEventTracker[centintME][zvintME].begin(),PoolsTklEventTracker[centintME][zvintME].begin()+2);
                              //  cout << "Discarded Event tracker front - New front number : "<< PoolsTklEventTracker[centintME].front() << endl;
                             //   cout << "Four first elements in PoolsTkl: " << PoolsTkl[centintME].at(0) << " " << PoolsTkl[centintME].at(1) << " " << PoolsTkl[centintME].at(2) << " " << PoolsTkl[centintME].at(3) << endl;
                              //  cout << "Size " << PoolsTkl[centintME].size()<<endl;
                                PoolsTkl[centintME][zvintME].erase(PoolsTkl[centintME][zvintME].begin(),PoolsTkl[centintME][zvintME].begin()+2);
                              //  cout << "Two first discarded - New Four first elements in PoolsTkl: " << PoolsTkl[centintME].at(0) << " " << PoolsTkl[centintME].at(1) << " " << PoolsTkl[centintME].at(2) << " " << PoolsTkl[centintME].at(3) << endl;
                            }
                            
                            PoolsSizeTkl[centintME][zvintME] -= 1;
                            
                            bool hasBeenFilled = kFALSE;
                           //  cout << "Will now add new event"<<endl;
                            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                               trackletME = (TrackletLight*)fTracklets->At(j);
                                //FIXME ok
                               if((trackletME->fEta < myEtaHigh(zvME)) && (trackletME->fEta > myEtaLow(zvME)) && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
                                   double trackletMEPhi = (Float_t)trackletME->fPhi + corrFactorDeltaPhi*(trackletME->fDPhi);
                                   double trackletMEEta = trackletME->fEta;
                                   PoolsTkl[centintME][zvintME].push_back(trackletMEPhi);
                                   PoolsTkl[centintME][zvintME].push_back(trackletMEEta);
                                   PoolsTklEventTracker[centintME][zvintME].push_back(TklMEcounter);
                                   PoolsTklEventTracker[centintME][zvintME].push_back(TklMEcounter);
                                   hasBeenFilled = kTRUE;
                               }

                           }
                            
                            if(hasBeenFilled){
                                PoolsSizeTkl[centintME][zvintME] += 1;
                            }
                            
                          //   cout << "New event added"<<endl;
                             
                        }
                        
                        else if(PoolsSizeTkl[centintME][zvintME] < 100){
                           //   cout << "Pool is not full - Will add event"<<endl;
                            bool hasBeenFilled = kFALSE;
                            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                                trackletME = (TrackletLight*)fTracklets->At(j);
                                //FIXME ok
                                if((trackletME->fEta < myEtaHigh(zvME)) && (trackletME->fEta > myEtaLow(zvME)) && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
                                    double trackletMEPhi = (Float_t)trackletME->fPhi + corrFactorDeltaPhi*(trackletME->fDPhi);
                                    double trackletMEEta = trackletME->fEta;
                                    PoolsTkl[centintME][zvintME].push_back(trackletMEPhi);
                                    PoolsTkl[centintME][zvintME].push_back(trackletMEEta);
                                    PoolsTklEventTracker[centintME][zvintME].push_back(TklMEcounter);
                                    PoolsTklEventTracker[centintME][zvintME].push_back(TklMEcounter);
                                    hasBeenFilled = kTRUE;
                                }

                            }
                            if(hasBeenFilled){
                                PoolsSizeTkl[centintME][zvintME] += 1;
                            }
                            //  cout << "Event added"<<endl;
                        }
                        
                        
                        
                        
                    }
                    
                }
            
            }
            
            
            
    }
        myfiletxt.close();
    }
    cout << "Reading of events finished" <<endl;

    TCanvas* cCorrMETKL=new TCanvas();
    cCorrMETKL->Divide(2,2);
    cCorrMETKL->SetTitle("Correlations Tkl-Tkl [0:Central][10: z_vtx = 0]");
    cCorrMETKL->cd(1);
    CorrelationsTkl[preciseCentFocus][preciseZbinFocus]->DrawCopy("colz");
    cCorrMETKL->cd(2);
    CorrelationsTklME[preciseCentFocus][preciseZbinFocus]->DrawCopy("colz");
    cCorrMETKL->cd(3);
    Correl_tampon = (TH1D*)(CorrelationsTkl[preciseCentFocus][preciseZbinFocus]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    cCorrMETKL->cd(4);
    Correl_tampon = (TH1D*)(CorrelationsTklME[preciseCentFocus][preciseZbinFocus]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    sprintf(CanvasName,"%s/CorrelationsTKL SE ME.pdf",FolderName);
    cCorrMETKL->SaveAs(CanvasName);

    
    if(doTracklets){
        TCanvas*ctestemichTklSE = new TCanvas();
        ctestemichTklSE->SetTitle("Naïve Tkl-Tkl SE Yield definition TH2 Correlations / Nb ref tracklets");
        ctestemichTklSE->Divide(1,3);
        
        YTklCentral->Scale(1./TklC);
        YTklPeriph->Scale(1./TklP);
        YTklDifference->Add(YTklCentral,YTklPeriph,1,-1);
        
        ctestemichTklSE->cd(1);
        YTklCentral->DrawCopy("colz");
        ctestemichTklSE->cd(2);
        YTklPeriph->DrawCopy("colz");
        ctestemichTklSE->cd(3);
        YTklDifference->DrawCopy("colz");
        sprintf(CanvasName,"%s/Naive TKL SE.pdf",FolderName);
        ctestemichTklSE->SaveAs(CanvasName);
        
        TCanvas*ctestemichTklSEProj = new TCanvas();
        ctestemichTklSEProj->SetTitle("Naïve Tkl-Tkl SE Yield definition TH2 Correlations / Nb ref tracklets - Projected");
        ctestemichTklSEProj->Divide(1,3);
        YTklCentral_proj_tampon = (TH1D*)(YTklCentral->ProjectionX("_px",1,-1,"e"));
        YTklPeriph_proj_tampon = (TH1D*)(YTklPeriph->ProjectionX("_px",1,-1,"e"));
        YTklDifference_proj_tampon = (TH1D*)(YTklDifference->ProjectionX("_px",1,-1,"e"));
        ctestemichTklSEProj->cd(1);
        
        YTklCentral_proj_tampon->DrawCopy("e");
        ctestemichTklSEProj->cd(2);
        YTklPeriph_proj_tampon->DrawCopy("e");
        ctestemichTklSEProj->cd(3);
        YTklDifference_proj_tampon->DrawCopy("e");
        sprintf(CanvasName,"%s/Naive TKL SE Proj.pdf",FolderName);
        ctestemichTklSEProj->SaveAs(CanvasName);
        
        if(doMixedEvents){
            YTklCentral->Scale(TklC);
            YTklPeriph->Scale(TklP);
            YTklCentralME->Scale(1./NormTklCentral);
            YTklPeriphME->Scale(1./NormTklPeriph);
            YTklCentral->Divide(YTklCentral,YTklCentralME);
            YTklPeriph->Divide(YTklPeriph,YTklPeriphME);
            YTklCentral->Scale(1./TklC);
            YTklPeriph->Scale(1./TklP);
            
            TCanvas*ctestemichTklME = new TCanvas();
            ctestemichTklME->SetTitle("Naïve Tkl-Tkl ME Yield definition TH2 Correlations ME / Nb ref tracklets");
            ctestemichTklME->Divide(1,3);
            ctestemichTklME->cd(1);
            YTklCentralME->DrawCopy("colz");
            ctestemichTklME->cd(2);
            YTklPeriphME->DrawCopy("colz");
            sprintf(CanvasName,"%s/Naive TKL ME.pdf",FolderName);
            ctestemichTklME->SaveAs(CanvasName);
            
            TCanvas*ctestemichTklMEProj = new TCanvas();
            ctestemichTklMEProj->SetTitle("Naïve Tkl-Tkl ME Yield definition TH2 Correlations ME / Nb ref tracklets - Projected");
            ctestemichTklMEProj->Divide(1,3);
            YTklCentralME_proj_tampon = (TH1D*)(YTklCentralME->ProjectionX("_px",1,-1,"e"));
            YTklPeriphME_proj_tampon = (TH1D*)(YTklPeriphME->ProjectionX("_px",1,-1,"e"));
            ctestemichTklMEProj->cd(1);
            
            YTklCentralME_proj_tampon->DrawCopy("e");
            ctestemichTklMEProj->cd(2);
            YTklPeriphME_proj_tampon->DrawCopy("e");
            sprintf(CanvasName,"%s/Naive TKL ME Proj.pdf",FolderName);
            ctestemichTklMEProj->SaveAs(CanvasName);
            
            TCanvas*ctestemichTklSEMEDiv = new TCanvas();
            ctestemichTklSEMEDiv->SetTitle("Naïve Tkl-Tkl SE/ME Yield definition TH2");
            ctestemichTklSEMEDiv->Divide(1,3);
            ctestemichTklSEMEDiv->cd(1);
            YTklCentral->DrawCopy("colz");
            ctestemichTklSEMEDiv->cd(2);
            YTklPeriph->DrawCopy("colz");
            sprintf(CanvasName,"%s/Naive TKL Yields.pdf",FolderName);
            ctestemichTklSEMEDiv->SaveAs(CanvasName);
            
        }
        
        YTklCentral_proj_tampon = (TH1D*)(YTklCentral->ProjectionX("_px",1,-1,"e"));
        YTklPeriph_proj_tampon = (TH1D*)(YTklPeriph->ProjectionX("_px",1,-1,"e"));
        YTklDifference->Add(YTklCentral,YTklPeriph,1,-1);
        YTklDifference_proj_tampon = (TH1D*)(YTklDifference->ProjectionX("YTklDifference_proj",1,-1,"e"));
        
        TCanvas*ctestemichTkl1D = new TCanvas();
        ctestemichTkl1D->SetTitle("Final Yields (SE/ME) Tkl-Tkl");
        ctestemichTkl1D->Divide(3,2);
        ctestemichTkl1D->cd(1);
        YTklCentral_proj_tampon->DrawCopy("E");
        ctestemichTkl1D->cd(2);
        YTklPeriph_proj_tampon->DrawCopy("E");
        ctestemichTkl1D->cd(3);
        YTklDifference_proj_tampon->DrawCopy("E");
        ctestemichTkl1D->cd(4);
        YTklCentral->DrawCopy("colz");
        ctestemichTkl1D->cd(5);
        YTklPeriph->DrawCopy("colz");
        ctestemichTkl1D->cd(6);
        YTklDifference->DrawCopy("colz");
        sprintf(CanvasName,"%s/Naive TKL Yields and Proj.pdf",FolderName);
        ctestemichTkl1D->SaveAs(CanvasName);
    }
    
    
    
    if(doMixedEvents){
    
        cout << "Normalisation of CorrelationsTklME based on max of Eta projected  applied" <<endl;
    
        //Scaling all CorrelationsMETkl
        cout << "Scaling CorrelationsMETkl" <<endl;
        
        if(doTracklets){
        for(int i=0; i<NbBinsCent; i++){
               for(int j=0; j<NbBinsZvtx; j++){
                 //  cout << i << " " << j <<endl;
                //   cout << "NormMETkl[i][j]:" << NormMETkl[i][j] <<endl;
                   if(NormMETkl[i][j]>0){
                        // CorrelationsTklME[i][j]->Scale(1./NormMETkl[i][j]);
                       for(int binx=1; binx<(1+CorrelationsTklME[i][j]->GetNbinsX()); binx++){
                           for(int biny=1; biny<(1+CorrelationsTklME[i][j]->GetNbinsY()); biny++){
                               CorrelationsTklMEScaled[i][j]->SetBinContent(binx,biny, (CorrelationsTklME[i][j]->GetBinContent(binx,biny))/NormMETkl[i][j]);
                                CorrelationsTklMEScaled[i][j]->SetBinError(binx,biny, (CorrelationsTklME[i][j]->GetBinError(binx,biny))/NormMETkl[i][j]);
                           }
                       }
                   }
                   
               }
        }
        }
        
        if(doTracklets && isCMSMethod){
        for(int i=0; i<NbBinsCent; i++){
               for(int j=0; j<NbBinsZvtx; j++){
                 //  cout << i << " " << j <<endl;
                //   cout << "NormMETkl[i][j]:" << NormMETkl[i][j] <<endl;
                   if(NormMETkl[i][j]>0){
                        // CorrelationsTklME[i][j]->Scale(1./NormMETkl[i][j]);
                       for(int binx=1; binx<(1+CorrelationsTklMEShortCMS[i][j]->GetNbinsX()); binx++){
                           for(int biny=1; biny<(1+CorrelationsTklMEShortCMS[i][j]->GetNbinsY()); biny++){
                               CorrelationsTklMEScaledShortCMS[i][j]->SetBinContent(binx,biny, (CorrelationsTklMEShortCMS[i][j]->GetBinContent(binx,biny))/NormMETkl[i][j]);
                                CorrelationsTklMEScaledShortCMS[i][j]->SetBinError(binx,biny, (CorrelationsTklMEShortCMS[i][j]->GetBinError(binx,biny))/NormMETkl[i][j]);
                           }
                       }
                       for(int binx=1; binx<(1+CorrelationsTklMELongCMS[i][j]->GetNbinsX()); binx++){
                           for(int biny=1; biny<(1+CorrelationsTklMELongCMS[i][j]->GetNbinsY()); biny++){
                               CorrelationsTklMEScaledLongCMS[i][j]->SetBinContent(binx,biny, (CorrelationsTklMELongCMS[i][j]->GetBinContent(binx,biny))/NormMETkl[i][j]);
                                CorrelationsTklMEScaledLongCMS[i][j]->SetBinError(binx,biny, (CorrelationsTklMELongCMS[i][j]->GetBinError(binx,biny))/NormMETkl[i][j]);
                           }
                       }
                   }

               }
        }
        }
}
    
    if(doTracklets){
        
        cout << "Calculation of number of tracklets seen integrated on z [Cent]" <<endl;
        for(int i=0; i<NbBinsCent; i++){
            for(int j=0; j<NbBinsZvtx; j++){
                RefTklCounterZint[i] += RefTklCounter[i][j];
            }
        }
        
        TCanvas*ctesteTkl = new TCanvas();
        ctesteTkl->SetTitle("SE Projected, ME projected and division Tkl-Tkl");
        ctesteTkl->Divide(1,3);
            
            TCanvas*ctestesumTkl = new TCanvas();
            ctestesumTkl->SetTitle("Sommation des yields selon z Tkl-Tkl");
            ctestesumTkl->Divide(5,2);
            
            TCanvas*ctestefinTkl = new TCanvas();
            ctestefinTkl->SetTitle("Yield[0] final sommé selon z Tkl-Tkl");
            
            TCanvas*ctesteprojTkl = new TCanvas();
            ctesteprojTkl->SetTitle("Projction SE et ME Tkl-Tkl");
            ctesteprojTkl->Divide(2,2);
        
        ctesteprojTkl->cd(1);
        CorrelationsTkl[preciseCentFocus][preciseZbinFocus]->DrawCopy("colz");
        ctesteprojTkl->cd(2);
        CorrelationsTklME[preciseCentFocus][preciseZbinFocus]->DrawCopy("colz");
        
        cout << "Calculation of YieldsTkl[Cent]" <<endl;
        for(int i=0; i<NbBinsCent; i++){
                      SoverMi->Reset();
                      for(int j=0; j<NbBinsZvtx; j++){
                          SoverMij->Reset();
                          ProjCopyTkl->Reset();
                          ProjCopy2Tkl->Reset();
                          if(doMixedEvents){
                              ME_proj_Tkl_tampon = (TH1F*)(CorrelationsTklMEScaled[i][j]->ProjectionX("ME_proj_Tkl",1,-1,"e")); //Change underflow
                          }
                          SE_proj_Tkl_tampon = (TH1I*)(CorrelationsTkl[i][j]->ProjectionX("SE_proj_Tkl",1,-1,"e")); //Change underflow
                          
                          ProjCopyTkl->Add(SE_proj_Tkl_tampon);
                          //   ProjCopy->Sumw2();
                          if(doMixedEvents){
                              ProjCopy2Tkl->Add(ME_proj_Tkl_tampon);
                          }
                          //   ProjCopy2->Sumw2();
                          if(doMixedEvents){
                              SoverMij->Divide(ProjCopyTkl,ProjCopy2Tkl);
                          }
                          else if(!doMixedEvents){
                              SoverMij->Add(ProjCopyTkl);
                          }
                             if(i ==preciseCentFocus && j==preciseZbinFocus){
                                 ctesteprojTkl->cd(3);
                                 ProjCopyTkl->DrawCopy();
                                 ctesteprojTkl->cd(4);
                                 ProjCopy2Tkl->DrawCopy();
                                 ctesteTkl->cd(1);
                                 ProjCopyTkl->DrawCopy();
                                 cout << "Relative error on SE Tkl bin 4 is: " << ProjCopyTkl->GetBinError(4)/ProjCopyTkl->GetBinContent(4) << endl;
                                  cout << "ProjCopyTkl->GetBinError(4): " << ProjCopyTkl->GetBinError(4) << endl;
                                 cout << "ProjCopyTkl->GetBinContent(4): " << ProjCopyTkl->GetBinContent(4) << endl;
                                 ctesteTkl->cd(2);
                                 if(doMixedEvents){
                                     ProjCopy2Tkl->DrawCopy();
                                     cout << "Relative error on ME Tkl bin 4 is: " << ProjCopy2Tkl->GetBinError(4)/ProjCopy2Tkl->GetBinContent(4) << endl;
                                     cout << "ProjCopy2Tkl->GetBinError(4): " << ProjCopy2Tkl->GetBinError(4) << endl;
                                     cout << "ProjCopy2Tkl->GetBinContent(4): " << ProjCopy2Tkl->GetBinContent(4) << endl;
                                 }
                                 ctesteTkl->cd(3);
                                 SoverMij->DrawCopy();
                                 cout << "Relative error on SE/ME Tkl bin 4 is: " << SoverMij->GetBinError(4)/SoverMij->GetBinContent(4) << endl;
                             }
                             if(i ==preciseCentFocus && j<5){
                                 ctestesumTkl->cd(j+1);
                                 SoverMij->DrawCopy();
                                 cout << "Relative error on Mij Tkl " << j << ", bin 4 is: " << SoverMij->GetBinError(4)/SoverMij->GetBinContent(4) << endl;
                                 cout << "Error on Mij Tkl " << j << ", bin 4 is: " << SoverMij->GetBinError(4) << endl;
                             }
                            SoverMi->Add(SoverMij);
                             if(i ==preciseCentFocus && j<5){
                                 ctestesumTkl->cd(5+j+1);
                                 SoverMi->DrawCopy();
                                 cout << "Relative error on Mi Tkl " << j << ", bin 4 is: " << SoverMi->GetBinError(4)/SoverMi->GetBinContent(4) << endl;
                                 cout << "Error on Mi Tkl " << j << ", bin 4 is: " << SoverMi->GetBinError(4) << endl;
                             }
                        
                          RefTracklets+=RefTklCounter[i][j];
                      }
                      YieldsTkl[i]->Add(SoverMi);
                      if(RefTklCounterZint[i] >0){
                        //  YieldsTkl[i]->Scale(1./RefTklCounterZint[i]); // ATTACHER ERREUR
                          for(int binx=1; binx<(1+YieldsTkl[i]->GetNbinsX()); binx++){

                                    double oldContent = YieldsTkl[i]->GetBinContent(binx);
                                    double oldError = YieldsTkl[i]->GetBinError(binx);
                                    YieldsTkl[i]->SetBinContent(binx, (YieldsTkl[i]->GetBinContent(binx))/RefTklCounterZint[i]);
                                    double newContent = YieldsTkl[i]->GetBinContent(binx);
                                    YieldsTkl[i]->SetBinError(binx, newContent*sqrt(pow(oldError/oldContent,2)+(1./RefTklCounterZint[i])));

                            }
                      }
              
            if(i==preciseCentFocus){
                ctestefinTkl->cd();
                YieldsTkl[i]->DrawCopy("colz");
            }
            
            TFile *f = new TFile(FitFileName,"UPDATE");
            YieldsTkl[i]->Write();
            f->Close();
            
//            ofstream myassociatetxt;
//            myassociatetxt.open(AssociateFileName);
//            myassociatetxt<< "RefTklCounterZint[" << i << "] = " << RefTklCounterZint[i] <<endl;
//            myassociatetxt.close();
            
          }
        
        sprintf(CanvasName,"%s/TKL SE ME Division.pdf",FolderName);
        ctesteTkl->SaveAs(CanvasName);
        sprintf(CanvasName,"%s/Z_SummationCheckTKL.pdf",FolderName);
        ctestesumTkl->SaveAs(CanvasName);
        sprintf(CanvasName,"%s/YieldTkl[0].pdf",FolderName);
        ctestefinTkl->SaveAs(CanvasName);
        printf(CanvasName,"%s/TKL SE ME.pdf",FolderName);
        ctesteprojTkl->SaveAs(CanvasName);
        
        cout << "Calculation of YieldsTkl for allC, central, periph" <<endl;
        for(int i=0; i<NbBinsCent; i++){
                YieldTkl_tampon->Reset();
                YieldTkl_tampon->Add(YieldsTkl[i]);
                YieldTkl_tampon->Scale(RefTklCounterZint[i]);
                YieldTkl_allC->Add(YieldTkl_tampon);
        }
        YieldTkl_allC->Scale(1./RefTracklets);
        
        int RefTklCnt=0;
        for(int i=CentralLowBound; i<CentralHighBound+1; i++){ //CentPeriph
            YieldTkl_tampon->Reset();
            YieldTkl_tampon->Add(YieldsTkl[i]);
            YieldTkl_tampon->Scale(RefTklCounterZint[i]);
            YieldTkl_Central->Add(YieldTkl_tampon);
            RefTklCnt += RefTklCounterZint[i];
        }
         YieldTkl_Central->Scale(1./RefTklCnt);
        RefTklCnt=0;
        for(int i=PeripheralLowBound; i<PeripheralHighBound+1; i++){
            cout << "i= " << i <<endl;
            YieldTkl_tampon->Reset();
            YieldTkl_tampon->Add(YieldsTkl[i]);
            cout << "Taken YieldsTkl[i]" <<endl;
            YieldTkl_tampon->Scale(RefTklCounterZint[i]);
            cout << "Scaled by RefTklCounterZint[i]: " << RefTklCounterZint[i] <<endl;
            YieldTkl_Periph->Add(YieldTkl_tampon);
            RefTklCnt += RefTklCounterZint[i];
            cout << "Now RefTklCnt is : " << RefTklCnt <<endl;
        }
        YieldTkl_Periph->Scale(1./RefTklCnt);
        YieldTkl_Difference->Add(YieldTkl_Central,YieldTkl_Periph,1,-1);
        
    }
    
    if(doTracklets && isCMSMethod){


           cout << "Calculation of YieldsTklShortCMS[Cent]" <<endl;
           for(int i=0; i<NbBinsCent; i++){
                         SoverMi->Reset();
                         for(int j=0; j<NbBinsZvtx; j++){
                             SoverMij->Reset();
                             ProjCopyTkl->Reset();
                             ProjCopy2Tkl->Reset();
                             if(doMixedEvents){
                                 ME_proj_Tkl_tampon = (TH1F*)(CorrelationsTklMEScaledShortCMS[i][j]->ProjectionX("ME_proj_Tkl",1,-1,"e")); //Change underflow
                             }
                             SE_proj_Tkl_tampon = (TH1I*)(CorrelationsTklShortCMS[i][j]->ProjectionX("SE_proj_Tkl",1,-1,"e")); //Change underflow

                             ProjCopyTkl->Add(SE_proj_Tkl_tampon);
                             //   ProjCopy->Sumw2();
                             if(doMixedEvents){
                                 ProjCopy2Tkl->Add(ME_proj_Tkl_tampon);
                             }
                             //   ProjCopy2->Sumw2();
                             if(doMixedEvents){
                                 SoverMij->Divide(ProjCopyTkl,ProjCopy2Tkl);
                             }
                             else if(!doMixedEvents){
                                 SoverMij->Add(ProjCopyTkl);
                             }


                               SoverMi->Add(SoverMij);

                         }
                         YieldsTklShortCMS[i]->Add(SoverMi);
                         if(RefTklCounterZint[i] >0){
                           //  YieldsTkl[i]->Scale(1./RefTklCounterZint[i]); // ATTACHER ERREUR
                             for(int binx=1; binx<(1+YieldsTklShortCMS[i]->GetNbinsX()); binx++){

                                       double oldContent = YieldsTklShortCMS[i]->GetBinContent(binx);
                                       double oldError = YieldsTklShortCMS[i]->GetBinError(binx);
                                       YieldsTklShortCMS[i]->SetBinContent(binx, (YieldsTklShortCMS[i]->GetBinContent(binx))/RefTklCounterZint[i]);
                                       double newContent = YieldsTklShortCMS[i]->GetBinContent(binx);
                                       YieldsTklShortCMS[i]->SetBinError(binx, newContent*sqrt(pow(oldError/oldContent,2)+(1./RefTklCounterZint[i])));

                               }
                         }


             }


           cout << "Calculation of YieldsTklShortCMS for allC, central, periph" <<endl;
           for(int i=0; i<NbBinsCent; i++){
                   YieldTkl_tampon->Reset();
                   YieldTkl_tampon->Add(YieldsTklShortCMS[i]);
                   YieldTkl_tampon->Scale(RefTklCounterZint[i]);
                   YieldTkl_allCShortCMS->Add(YieldTkl_tampon);
           }
           YieldTkl_allCShortCMS->Scale(1./RefTracklets);

           int RefTklCnt=0;
           for(int i=CentralLowBound; i<CentralHighBound+1; i++){ //CentPeriph
               YieldTkl_tampon->Reset();
               YieldTkl_tampon->Add(YieldsTklShortCMS[i]);
               YieldTkl_tampon->Scale(RefTklCounterZint[i]);
               YieldTkl_CentralShortCMS->Add(YieldTkl_tampon);
               RefTklCnt += RefTklCounterZint[i];
           }
            YieldTkl_CentralShortCMS->Scale(1./RefTklCnt);
           RefTklCnt=0;
           for(int i=PeripheralLowBound; i<PeripheralHighBound+1; i++){
               YieldTkl_tampon->Reset();
               YieldTkl_tampon->Add(YieldsTklShortCMS[i]);
               YieldTkl_tampon->Scale(RefTklCounterZint[i]);
               YieldTkl_PeriphShortCMS->Add(YieldTkl_tampon);
               RefTklCnt += RefTklCounterZint[i];
           }
           YieldTkl_PeriphShortCMS->Scale(1./RefTklCnt);
           YieldTkl_DifferenceShortCMS->Add(YieldTkl_CentralShortCMS,YieldTkl_PeriphShortCMS,1,-1);

       }
    
    if(doTracklets && isCMSMethod){


        cout << "Calculation of YieldsTklLongCMS[Cent]" <<endl;
        for(int i=0; i<NbBinsCent; i++){
                      SoverMi->Reset();
                      for(int j=0; j<NbBinsZvtx; j++){
                          SoverMij->Reset();
                          ProjCopyTkl->Reset();
                          ProjCopy2Tkl->Reset();
                          if(doMixedEvents){
                              ME_proj_Tkl_tampon = (TH1F*)(CorrelationsTklMEScaledLongCMS[i][j]->ProjectionX("ME_proj_Tkl",1,-1,"e")); //Change underflow
                          }
                          SE_proj_Tkl_tampon = (TH1I*)(CorrelationsTklLongCMS[i][j]->ProjectionX("SE_proj_Tkl",1,-1,"e")); //Change underflow

                          ProjCopyTkl->Add(SE_proj_Tkl_tampon);
                          //   ProjCopy->Sumw2();
                          if(doMixedEvents){
                              ProjCopy2Tkl->Add(ME_proj_Tkl_tampon);
                          }
                          //   ProjCopy2->Sumw2();
                          if(doMixedEvents){
                              SoverMij->Divide(ProjCopyTkl,ProjCopy2Tkl);
                          }
                          else if(!doMixedEvents){
                              SoverMij->Add(ProjCopyTkl);
                          }


                            SoverMi->Add(SoverMij);

                      }
                      YieldsTklLongCMS[i]->Add(SoverMi);
                      if(RefTklCounterZint[i] >0){
                        //  YieldsTkl[i]->Scale(1./RefTklCounterZint[i]); // ATTACHER ERREUR
                          for(int binx=1; binx<(1+YieldsTklLongCMS[i]->GetNbinsX()); binx++){

                                    double oldContent = YieldsTklLongCMS[i]->GetBinContent(binx);
                                    double oldError = YieldsTklLongCMS[i]->GetBinError(binx);
                                    YieldsTklLongCMS[i]->SetBinContent(binx, (YieldsTklLongCMS[i]->GetBinContent(binx))/RefTklCounterZint[i]);
                                    double newContent = YieldsTklLongCMS[i]->GetBinContent(binx);
                                    YieldsTklLongCMS[i]->SetBinError(binx, newContent*sqrt(pow(oldError/oldContent,2)+(1./RefTklCounterZint[i])));

                            }
                      }


          }


        cout << "Calculation of YieldsTklLongCMS for allC, central, periph" <<endl;
        for(int i=0; i<NbBinsCent; i++){
                YieldTkl_tampon->Reset();
                YieldTkl_tampon->Add(YieldsTklLongCMS[i]);
                YieldTkl_tampon->Scale(RefTklCounterZint[i]);
                YieldTkl_allCLongCMS->Add(YieldTkl_tampon);
        }
        YieldTkl_allCLongCMS->Scale(1./RefTracklets);

        int RefTklCnt=0;
        for(int i=CentralLowBound; i<CentralHighBound+1; i++){ //CentPeriph
            YieldTkl_tampon->Reset();
            YieldTkl_tampon->Add(YieldsTklLongCMS[i]);
            YieldTkl_tampon->Scale(RefTklCounterZint[i]);
            YieldTkl_CentralLongCMS->Add(YieldTkl_tampon);
            RefTklCnt += RefTklCounterZint[i];
        }
         YieldTkl_CentralLongCMS->Scale(1./RefTklCnt);
        RefTklCnt=0;
        for(int i=PeripheralLowBound; i<PeripheralHighBound+1; i++){
            YieldTkl_tampon->Reset();
            YieldTkl_tampon->Add(YieldsTklLongCMS[i]);
            YieldTkl_tampon->Scale(RefTklCounterZint[i]);
            YieldTkl_PeriphLongCMS->Add(YieldTkl_tampon);
            RefTklCnt += RefTklCounterZint[i];
        }
        YieldTkl_PeriphLongCMS->Scale(1./RefTklCnt);
        YieldTkl_DifferenceLongCMS->Add(YieldTkl_CentralLongCMS,YieldTkl_PeriphLongCMS,1,-1);

    }
    
    std::cout << "==================== Analysis Terminated - PLOTTING TIME ====================" << std::endl;
    //}
    
    
// ***********************************
// Tracer les graphes sur les canevas *
// ***********************************
    
    double baselineTKL = 1;
    double errbaselineTKL = 1;
    
    if(doTracklets){
        
        TCanvas* c6TKL = new TCanvas;
            //Tracklets yield DeltaEta wrt DeltaPhi TH2 -> Projected for Periph and Central

            c6TKL->Divide(2,2);
            c6TKL->cd(1);
            YieldTkl_allC->Draw("E");
            c6TKL->cd(2);
            YieldTkl_Difference->Draw("E");
            c6TKL->cd(3);
            YieldTkl_Central->Draw("E");
            c6TKL->cd(4);
            YieldTkl_Periph->Draw("E");
            c6TKL->Draw();
            c6TKL->Modified();
            c6TKL->ForceUpdate();
        
        sprintf(CanvasName,"%s/V2TKL_Sub.pdf",FolderName);
        c6TKL->SaveAs(CanvasName);
        
        TFile *f = new TFile(FitFileName,"UPDATE");
        YieldTkl_allC->Write();
        YieldTkl_Periph->Write();
        YieldTkl_Central->Write();
        YieldTkl_Difference->Write();
        f->Close();
        
        baselineTKL = (YieldTkl_Periph->GetBinContent(BinZeroLeftTKL) + YieldTkl_Periph->GetBinContent(BinZeroLeftTKL+1))/2;
        errbaselineTKL = sqrt(pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL),2) + pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL+1),2));
            
        TCanvas*c14TKL=new TCanvas();
        //Tracklets Yield difference wrt Phi fit
        c14TKL->cd();
        // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
            TH1F *histo = YieldTkl_Difference;
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
           Double_t params[3] = {1,0.01,0.01};
          fitFcnV2TKL->SetParameters(params);
           TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
           TVirtualFitter::Fitter(histo)->SetPrecision();
        //  histo->Fit("fitFcn","0");
          // second try: set start values for some parameters
           
           fitFcnV2TKL->SetParName(0,"a0");
           fitFcnV2TKL->SetParName(1,"a1");
           fitFcnV2TKL->SetParName(2,"a2");
            fitFcnV2TKL->SetParName(3,"a3");
//            fitFcnV2TKL->SetParName(4,"a4");
//        fitFcnV2TKL->SetParName(5,"a5");
//        fitFcnV2TKL->SetParName(6,"a6");
//        fitFcnV2TKL->SetParName(7,"a7");
//        fitFcnV2TKL->SetParName(8,"a8");
//        fitFcnV2TKL->SetParName(9,"a9");
//        fitFcnV2TKL->SetParName(10,"a10");
//        fitFcnV2TKL->SetParName(11,"a11");
//        fitFcnV2TKL->SetParLimits(6,-0.1,0.1);
           
          TFitResultPtr res = histo->Fit("fitFcnV2TKL","SBMERI+","ep");
        gStyle->SetOptStat("n");
        gStyle->SetOptFit(1011);
        if(!doTracklets){
        TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
        st->SetX1NDC(0.8); //new x start position
        st->SetY1NDC(0.8); //new x end position
        }
          // improve the pictu
        //   std::cout << "integral error: " << integralerror << std::endl;
            Double_t par[3];
            fitFcnV2TKL->GetParameters(par);
          fitFcnV2TKL->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
          legend->SetFillColorAlpha(kWhite, 0.);
          legend->SetBorderSize(0);
            legend->SetTextFont(42);
          legend->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
            legend->AddEntry(fitFcnV2TKL,message);
        sprintf(message,"V_{2,tkl-tkl}: #frac{%.4f}{%.4f + %.4f} = %.4f +- %.4f",par[2],par[0], baselineTKL,par[2]/(par[0] + baselineTKL),(par[2]/(par[0] + baselineTKL)*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)+pow(errbaselineTKL/baselineTKL,2))));
        legend->AddEntry(fitFcnV2TKL,message);
        if(res->CovMatrixStatus() == 3){
                 //  sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
                   legend->AddEntry(fitFcnV2TKL,message);
               }
          legend->AddEntry(histo,"Data","lpe");
          legend->Draw();
        
        sprintf(CanvasName,"%s/V2TKL.pdf",FolderName);
        c14TKL->SaveAs(CanvasName);
        
    }
    
    Double_t Acoef_Periph[4];
    Double_t Acoef_Central[4];
    Double_t errAcoef_Periph[4];
    Double_t errAcoef_Central[4];
    
    if (doTracklets && isCMSMethod){
            
            TCanvas*c14TKLperiph=new TCanvas();
            //Tracklets Yield difference wrt Phi fit
            c14TKLperiph->cd();
            // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
                TH1F *histo = YieldTkl_Periph;
              // create a TF1 with the range from 0 to 3 and 6 parameters
              TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV5,-TMath::Pi()/2,1.5*TMath::Pi(),4);
              fitFcnV2TKL->SetNpx(500);
              fitFcnV2TKL->SetLineWidth(4);
              fitFcnV2TKL->SetLineColor(kMagenta);
              // first try without starting values for the parameters
              // This defaults to 1 for each param.
              // this results in an ok fit for the polynomial function
              // however the non-linear part (lorenzian) does not
              // respond well.
               Double_t params[4] = {1,0.01,0.01, 0.01};
              fitFcnV2TKL->SetParameters(params);
               TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
               TVirtualFitter::Fitter(histo)->SetPrecision();
            //  histo->Fit("fitFcn","0");
              // second try: set start values for some parameters
               
               fitFcnV2TKL->SetParName(0,"a0");
               fitFcnV2TKL->SetParName(1,"a1");
               fitFcnV2TKL->SetParName(2,"a2");
                fitFcnV2TKL->SetParName(3,"a3");
    //            fitFcnV2TKL->SetParName(4,"a4");
    //        fitFcnV2TKL->SetParName(5,"a5");
    //        fitFcnV2TKL->SetParName(6,"a6");
    //        fitFcnV2TKL->SetParName(7,"a7");
    //        fitFcnV2TKL->SetParName(8,"a8");
    //        fitFcnV2TKL->SetParName(9,"a9");
    //        fitFcnV2TKL->SetParName(10,"a10");
    //        fitFcnV2TKL->SetParName(11,"a11");
    //        fitFcnV2TKL->SetParLimits(6,-0.1,0.1);
               
              TFitResultPtr res = histo->Fit("fitFcnV2TKL","SBMERI+","ep");
        gStyle->SetOptStat("n");
        gStyle->SetOptFit(1011);
        if(!doTracklets){
        TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
        st->SetX1NDC(0.8); //new x start position
        st->SetY1NDC(0.8); //new x end position
        }
              // improve the pictu
            //   std::cout << "integral error: " << integralerror << std::endl;
                fitFcnV2TKL->GetParameters(Acoef_Periph);
        for(int index=0; index<4;index++){
            errAcoef_Periph[index] = fitFcnV2TKL->GetParError(index);
        }
              fitFcnV2TKL->Draw("same");
              // draw the legend
              TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
              legend->SetFillColorAlpha(kWhite, 0.);
              legend->SetBorderSize(0);
                legend->SetTextFont(42);
              legend->SetTextSize(0.03);
                Char_t message[80];
                sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
                legend->AddEntry(fitFcnV2TKL,message);
            sprintf(message,"V_{2,tkl-tkl} Periph: #frac {%.4f}{%.4f} = %.4f +- %.4f",Acoef_Periph[2],Acoef_Periph[0],Acoef_Periph[2]/(Acoef_Periph[0]),(Acoef_Periph[2]/(Acoef_Periph[0])*sqrt(pow(fitFcnV2TKL->GetParError(2)/Acoef_Periph[2],2)+pow(fitFcnV2TKL->GetParError(0)/Acoef_Periph[0],2))));
            legend->AddEntry(fitFcnV2TKL,message);
            if(res->CovMatrixStatus() == 3){
                       //sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                       legend->AddEntry(fitFcnV2TKL,message);
                   }
              legend->AddEntry(histo,"Data","lpe");
              legend->Draw();
        
        }
    
    if (doTracklets && isCMSMethod){
            
            TCanvas*c14TKLcent=new TCanvas();
            //Tracklets Yield difference wrt Phi fit
            c14TKLcent->cd();
            // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
                TH1F *histo = YieldTkl_Central;
              // create a TF1 with the range from 0 to 3 and 6 parameters
              TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV5,-TMath::Pi()/2,1.5*TMath::Pi(),4);
              fitFcnV2TKL->SetNpx(500);
              fitFcnV2TKL->SetLineWidth(4);
              fitFcnV2TKL->SetLineColor(kMagenta);
              // first try without starting values for the parameters
              // This defaults to 1 for each param.
              // this results in an ok fit for the polynomial function
              // however the non-linear part (lorenzian) does not
              // respond well.
               Double_t params[4] = {1,0.01,0.01, 0.01};
              fitFcnV2TKL->SetParameters(params);
               TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
               TVirtualFitter::Fitter(histo)->SetPrecision();
            //  histo->Fit("fitFcn","0");
              // second try: set start values for some parameters
               
               fitFcnV2TKL->SetParName(0,"a0");
               fitFcnV2TKL->SetParName(1,"a1");
               fitFcnV2TKL->SetParName(2,"a2");
                fitFcnV2TKL->SetParName(3,"a3");
    //            fitFcnV2TKL->SetParName(4,"a4");
    //        fitFcnV2TKL->SetParName(5,"a5");
    //        fitFcnV2TKL->SetParName(6,"a6");
    //        fitFcnV2TKL->SetParName(7,"a7");
    //        fitFcnV2TKL->SetParName(8,"a8");
    //        fitFcnV2TKL->SetParName(9,"a9");
    //        fitFcnV2TKL->SetParName(10,"a10");
    //        fitFcnV2TKL->SetParName(11,"a11");
    //        fitFcnV2TKL->SetParLimits(6,-0.1,0.1);
               
              TFitResultPtr res = histo->Fit("fitFcnV2TKL","SBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
        if(!doTracklets){
            TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
            st->SetX1NDC(0.8); //new x start position
            st->SetY1NDC(0.8); //new x end position
        }
              // improve the pictu
            //   std::cout << "integral error: " << integralerror << std::endl;
                fitFcnV2TKL->GetParameters(Acoef_Central);
                for(int index=0; index<4;index++){
                    errAcoef_Central[index] = fitFcnV2TKL->GetParError(index);
                }
              fitFcnV2TKL->Draw("same");
              // draw the legend
              TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
              legend->SetFillColorAlpha(kWhite, 0.);
              legend->SetBorderSize(0);
                legend->SetTextFont(42);
              legend->SetTextSize(0.03);
                Char_t message[80];
                sprintf(message,"Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
                legend->AddEntry(fitFcnV2TKL,message);
            sprintf(message,"V_{2,tkl-tkl} Central: #frac{%.4f}{%.4f} = %.4f +- %.4f",Acoef_Central[2],Acoef_Central[0],Acoef_Central[2]/(Acoef_Central[0]),(Acoef_Central[2]/(Acoef_Central[0])*sqrt(pow(fitFcnV2TKL->GetParError(2)/Acoef_Central[2],2)+pow(fitFcnV2TKL->GetParError(0)/Acoef_Central[0],2))));
            legend->AddEntry(fitFcnV2TKL,message);
            if(res->CovMatrixStatus() == 3){
                     //  sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                       legend->AddEntry(fitFcnV2TKL,message);
                   }
              legend->AddEntry(histo,"Data","lpe");
              legend->Draw();
        
        }
    
    if(doTracklets && isCMSMethod){
        Double_t V2_CMS_Periph[4];
        Double_t V2_CMS_Central[4];
        Double_t V2_CMS_Difference[4];
        Double_t errV2_CMS_Periph[4];
        Double_t errV2_CMS_Central[4];
        Double_t errV2_CMS_Difference[4];
        Double_t errCompound[4];
        
        for(int index=0; index<4; index++){
            cout << "Acoef_Periph_" << index << " = " << Acoef_Periph[index] <<" +/- "<<errAcoef_Periph[index]<<endl;
        }
        for(int index=0; index<4; index++){
                   cout << "Acoef_Central_" << index << " = " << Acoef_Central[index] <<" +/- "<< errAcoef_Central[index] <<endl;
               }
        
        //Calculation V2 Periph and Central
        for(int index=0; index<4; index++){
            V2_CMS_Periph[index]=Acoef_Periph[index]/Acoef_Periph[0];
            V2_CMS_Central[index]=Acoef_Central[index]/Acoef_Central[0];
            
            errV2_CMS_Periph[index] = V2_CMS_Periph[index]*sqrt(pow(errAcoef_Periph[index]/Acoef_Periph[index],2)+pow(errAcoef_Periph[0]/Acoef_Periph[0],2));
            errV2_CMS_Central[index] = V2_CMS_Central[index]*sqrt(pow(errAcoef_Central[index]/Acoef_Central[index],2)+pow(errAcoef_Central[0]/Acoef_Central[0],2));
        }
        
        for(int index=0; index<4; index++){
                   cout << "V2_Periph_" << index << " = " << V2_CMS_Periph[index] <<" +/- "<<errV2_CMS_Periph[index]<<endl;
               }
               for(int index=0; index<4; index++){
                          cout << "V2_Central_" << index << " = " << V2_CMS_Central[index] <<" +/- "<< errV2_CMS_Central[index] <<endl;
                      }
        
        TCanvas* ctry=new TCanvas();
        ctry->Divide(1,2);
        ctry->cd(1);
        YieldTkl_CentralLongCMS->DrawCopy();
        ctry->cd(2);
        YieldTkl_CentralShortCMS->DrawCopy();
        
        TCanvas* ctry2=new TCanvas();
        ctry2->Divide(1,2);
        ctry2->cd(1);
        YieldTkl_PeriphLongCMS->DrawCopy();
        ctry2->cd(2);
        YieldTkl_PeriphShortCMS->DrawCopy();
        
        //Calculation Y_jet_Central and Y_jet_Periph
        Double_t e1;
        Double_t e2;
        Double_t e3;
        Double_t e4;
        
        Double_t Y_jet_Central = -(YieldTkl_CentralLongCMS->IntegralAndError(1,YieldTkl_CentralLongCMS->GetNbinsX(),e1) - YieldTkl_CentralShortCMS->IntegralAndError(1,YieldTkl_CentralShortCMS->GetNbinsX(),e2));
        Double_t Y_jet_Periph = -(YieldTkl_PeriphLongCMS->IntegralAndError(1,YieldTkl_PeriphLongCMS->GetNbinsX(),e3) - YieldTkl_PeriphShortCMS->IntegralAndError(1,YieldTkl_PeriphShortCMS->GetNbinsX(),e4)) ;
        Double_t errY_jet_Central = sqrt(pow(e1,2)+pow(e2,2));
        Double_t errY_jet_Periph = sqrt(pow(e3,2)+pow(e4,2));
        
        cout << "Y_jet_Central = " << Y_jet_Central << " +/- " << errY_jet_Central <<endl;
        cout << "Y_jet_Periph = " << Y_jet_Periph << " +/- " << errY_jet_Periph <<endl;
        
        cout << "Nassoc_Central_TKL = " << Nassoc_Central_TKL<<endl;
        cout << "Nassoc_Periph_TKL = " << Nassoc_Periph_TKL<<endl;
        
        Double_t factorCMS = static_cast<Double_t>(double(Nassoc_Periph_TKL)/Nassoc_Central_TKL)*(Y_jet_Central/Y_jet_Periph);
        Double_t errfactorCMS = factorCMS*sqrt(pow(1/sqrt(Nassoc_Periph_TKL),2) + pow(1/sqrt(Nassoc_Central_TKL),2) + pow(errY_jet_Central/Y_jet_Central,2) + pow(errY_jet_Periph/Y_jet_Periph,2));
        
        cout << "factorCMS = " << factorCMS<<endl;
        
        for(int index=0; index<4; index++){
            V2_CMS_Difference[index] = V2_CMS_Central[index] - V2_CMS_Periph[index]*factorCMS;
            errCompound[index] = V2_CMS_Periph[index]*factorCMS*sqrt(pow(errfactorCMS/factorCMS,2) + pow(errV2_CMS_Periph[index]/V2_CMS_Periph[index],2));
            errV2_CMS_Difference[index] = sqrt(pow(errV2_CMS_Central[index],2)+pow(errCompound[index],2));
            
            cout << "V_CMSDiff_" << index << " = " << V2_CMS_Difference[index] << " +/- " <<  errV2_CMS_Difference[index]<<endl;
        }

    }
    
    
}
    









// FITTING INVARIANT MASS METHODS



Double_t FourierV2(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);}

Double_t FourierV5(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]) + 2*par[3]*cos(3*x[0]);}// + 2*par[4]*cos(4*x[0]);}// + 2*par[5]*cos(5*x[0]) + 2*par[6]*cos(6*x[0]) + 2*par[7]*cos(7*x[0]) + 2*par[8]*cos(8*x[0]) + 2*par[9]*cos(9*x[0]) + 2*par[10]*cos(10*x[0]) + 2*par[11]*cos(11*x[0]);}


int GetCent(float cent){
    if(cent <= 1){
        return 0;
    }
    else if(cent <= 3){
        return 1;
    }
    else if(cent <= 5){
        return 2;
    }
    else if(cent <= 10){
        return 3;
    }
    else if(cent <= 15){
        return 4;
    }
    else if(cent <= 20){
        return 5;
    }
    else if(cent <= 30){
        return 6;
    }
    else if(cent <= 40){
        return 7;
    }
    else if(cent <= 50){
        return 8;
    }
    else if(cent <= 60){
        return 9;
    }
    else if(cent <= 70){
        return 10;
    }
    else if(cent <= 80){
        return 11;
    }
    else if(cent <= 90){
        return 12;
    }
    else{
        return 13;
    }
}

int GetCentPM(int Ntkl, float SPDTrackletsPer, float SPDClustersPer, float V0MPer, int zvtx_idx, int groupnumber){
    
    float estimator = 0;

    if(centralityMethod == "PercentileMethodSPDTracklets"){
        if(!isMultiplicityStudy){
            if(Ntkl==0){
                return 13;
            }

            for(int cent_index=0;cent_index<14;cent_index++){
                if(LimitsPM_SPDTracklets_Uncal_DataDriven[groupnumber-1][zvtx_idx][cent_index] < Ntkl){
                    return cent_index;
                }
            }
        }
        else if(isMultiplicityStudy){
            if(Ntkl<=NtklPeripheralHighBound){
                           return 13;
                       }
            else if(Ntkl>NtklCentralLowBound){
                return 0;
            }
            else{
                return int((PeripheralHighBound+CentralLowBound)*0.5);
            }
            
        }
    }
    else if(centralityMethod == "SPDTrackletsPercentile"){
        estimator = SPDTrackletsPer;
    }
    else if(centralityMethod == "SPDClustersPercentile"){
        estimator = SPDClustersPer;
    }
    else if(centralityMethod == "V0MPercentile"){
        estimator = V0MPer;
    }
    else{
        return -1;
    }
    return GetCent(estimator);
    
}

int GetCentCvetan(float V0MPer){
    
    float estimator = 0;

    if(centralityMethod == "V0MPercentile"){
        estimator = V0MPer;
    }
    else{
        return -1;
    }
    return GetCent(estimator);
    
}

Double_t myEtaLow(double x){
    return TMath::Log(TMath::Tan(0.5*TMath::ATan(7.6/(14.1+(floor(x)+0.5))))) + 0.1; }

Double_t myEtaHigh(double x){
return -(TMath::Log(TMath::Tan(0.5*TMath::ATan(7.6/(14.1-(floor(x)+0.5)))))) - 0.05; }

