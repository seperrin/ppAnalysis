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
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <string>
#include <vector>

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
# include <string>
# include "TH1F.h"
# include "TH1D.h"
# include "TH2F.h"
# include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"


// FUNCTIONS

void PlotSingleBin();
Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par);
Double_t BackFcnV2(Double_t *x,Double_t *par);
Double_t BackFcnV2Poly(Double_t *x,Double_t *par);
Double_t SignalFcnJPsiV2(Double_t *x,Double_t *par);
Double_t FourierV2(Double_t *x,Double_t *par);
Double_t FourierV5(Double_t *x,Double_t *par);
Double_t TwoCBE2E(Double_t *x,Double_t *par);
Double_t ExpBkg(Double_t *x,Double_t *par);
Double_t JPsiCrystalBallExtended(Double_t *x,Double_t *par);
Double_t Psi2SCrystalBallExtended(Double_t *x,Double_t *par);
TFitResultPtr FittingAllInvMass(const char *histoname, TCanvas *canvas);
TFitResultPtr FittingAllInvMassBin(const char *histoname, TCanvas *canvas, int i);
void FcnCombinedAllMass(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  );
int GetCent(double cent);
int GetCentPM(int Ntkl, int Zvtx, int groupnumber);
bool isCentral(int centint);
bool isPeripheral(int centint);

// Centrality bins:
//0: 0-1%
//1: 1-5%
//2: 5-10%
//3: 10-15%
//4: 15-20%
//5: 20-30%
//6: 30-40%
//7: 40-50%
//8: 50-60%
//9: 60-70%
//10: 70-80%
//11: 80-90%
//12: 90-100&

int CentralLowBound = 0;
int CentralHighBound = 0;
int PeripheralLowBound = 7; //7
int PeripheralHighBound = 11; //11

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


int LimitsPM[12][20][13] =
{       {   {24, 17, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0}, //Group1
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {31, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {31, 22, 17, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {31, 22, 17, 15, 13, 10, 7, 6, 5, 3, 3, 2, 0},
            {32, 22, 18, 15, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {32, 22, 18, 15, 13, 10, 8, 6, 5, 3, 3, 2, 0},
            {31, 22, 17, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {24, 17, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0}, //Group2
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 16, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {31, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {24, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0}, //Group3
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 16, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {31, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {31, 22, 17, 14, 12, 9, 7, 6, 4, 3, 2, 2, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {23, 16, 12, 10, 9, 7, 5, 4, 3, 3, 2, 1, 0}, //Group4
            {24, 17, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 16, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {25, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}, //Group5
            {27, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {31, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {32, 22, 17, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {32, 22, 18, 15, 13, 10, 8, 6, 5, 3, 2, 2, 0},
            {33, 23, 18, 15, 13, 10, 8, 6, 5, 4, 3, 2, 0},
            {34, 23, 18, 15, 13, 10, 8, 6, 5, 4, 3, 2, 0},
            {34, 24, 19, 16, 13, 10, 8, 6, 5, 4, 3, 2, 0},
            {34, 24, 19, 16, 14, 10, 8, 6, 5, 4, 3, 2, 0},
            {34, 23, 18, 15, 13, 10, 8, 6, 5, 4, 3, 2, 0},
            {32, 22, 18, 15, 13, 10, 7, 6, 4, 3, 2, 2, 0},
            {30, 21, 17, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0}   },
        {   {21, 15, 12, 10, 8, 6, 5, 4, 3, 2, 2, 1, 0}, //Group6
            {22, 15, 12, 10, 9, 7, 5, 4, 3, 2, 2, 1, 0},
            {23, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {30, 21, 17, 14, 12, 9, 7, 6, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0}   },
        {   {21, 15, 12, 10, 8, 6, 5, 4, 3, 2, 2, 1, 0}, //Group7
            {22, 16, 12, 10, 9, 7, 5, 4, 3, 3, 2, 1, 0},
            {23, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0}   },
        {   {21, 15, 12, 10, 8, 6, 5, 4, 3, 2, 2, 1, 0}, //Group8
            {22, 15, 12, 10, 9, 7, 5, 4, 3, 3, 2, 1, 0},
            {23, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {24, 17, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0}   },
        {   {21, 15, 12, 10, 8, 6, 5, 4, 3, 2, 2, 1, 0}, //Group9
            {22, 15, 12, 10, 9, 7, 5, 4, 3, 3, 2, 1, 0},
            {23, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 14, 12, 9, 7, 5, 4, 3, 2, 1, 0},
            {29, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0}   },
        {   {22, 15, 12, 10, 8, 7, 5, 4, 3, 2, 2, 1, 0}, //Group10
            {22, 16, 12, 10, 9, 7, 5, 4, 3, 3, 2, 1, 0},
            {23, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0},
            {23, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0}   },
        {   {23, 16, 13, 11, 9, 7, 5, 4, 3, 3, 2, 1, 0}, //Group11
            {24, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {24, 17, 13, 11, 10, 7, 6, 5, 4, 3, 2, 1, 0},
            {23, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0}   },
        {   {23, 16, 13, 11, 9, 7, 5, 4, 3, 3, 2, 1, 0}, //Group12
            {24, 17, 13, 11, 10, 7, 6, 4, 3, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {26, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {27, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 20, 16, 13, 11, 9, 7, 5, 4, 3, 2, 1, 0},
            {28, 19, 15, 13, 11, 8, 7, 5, 4, 3, 2, 1, 0},
            {27, 18, 15, 12, 11, 8, 6, 5, 4, 3, 2, 1, 0},
            {25, 17, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0},
            {24, 16, 13, 11, 9, 7, 6, 4, 3, 3, 2, 1, 0}   }

};

Double_t mJpsi =  3.096916;
 Double_t mPsip =  3.686108;
// Double_t ratMass = 1.01;
 Double_t ratSigma = 1.05;
Int_t npfits;

TH1F* hnseg(NULL);
TH1F* V2JPsiTkl(NULL);

void PlotSingleBin(){
 //   freopen( "logPlotFromTreeJavier16h_NoDPhiCut_NoConstraintPileUp.txt", "w", stdout );
    TH1::SetDefaultSumw2();
    bool doTracklets = kTRUE;
    bool doMixedEvents = kTRUE;
    bool CombineFits = kFALSE;
    bool KeepOnlyOne = kFALSE;
    int valueOnlyOne = 3;
    bool AdditionalCutNtkl = kFALSE;
    int valueAdditionalCutNtkl = 3; //<=
    int preciseCentFocus = 0;
    
// ************************************
// Définitions de paramètres          *
// ************************************

    
    double ZvtxBin = 10;
    double ZvtxCut = 10;
    double SigmaZvtxCut = 0.25;
    double DPhiCut = 0.01;
    double TklEtaCut  = 1;
    double LowDimuYCut = -4;
    double HighDimuYCut = -2.5;
    double LowDimuPtCut = 0;
    double HighDimuPtCut = 9999;
//    double MinMultCentral = 37;
//    double MaxMultPeriph = 23;
    double CentSPDTrackletsCentral = 1;
    double CentSPDTrackletsPeriph = 40;

    double LowJPsiMass = 3.0;
    double HighJPsiMass = 3.3;

    double MinInvMass = 2.1;
    double MaxInvMass = 5.1;
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
    
    double MinDeltaEtaTKL = -2.4; //1.2
    double MaxDeltaEtaTKL = 2.4;
    const int NbinsDeltaEtaTKL = 48; //24
    double DeltaEtaTKLCut = 1.2; //1.2
    double SizeBinDeltaEtaTKL = (MaxDeltaEtaTKL-MinDeltaEtaTKL)/NbinsDeltaEtaTKL;
    
    const int NbinsDeltaPhiTKL = 48;
    double SizeBinDeltaPhiTKL = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhiTKL;
    
    const int NbBinsCent = 13;
    const int NbBinsZvtx = 20;
    
    
  //  Char_t Group_Period[50] = "Group1";
//    Char_t *arrayOfPeriods[] = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l","Group2_LHC17h","Group3_LHC17h","Group4_LHC17k","Group4_LHC18l","Group4_LHC18m","Group4_LHC18o","Group4_LHC18p","Group5_LHC17l","Group5_LHC17m","Group5_LHC17o","Group5_LHC17r","Group5_LHC18c","Group5_LHC18d","Group5_LHC18e","Group5_LHC18f","Group6_LHC18m","Group7_LHC18m","Group8_LHC18m","Group9_LHC18m","Group10_LHC18m","Group11_LHC18m","Group12_LHC18m"};
    Char_t *arrayOfPeriods[] = {"Group1_LHC16h","Group1_LHC16j","Group1_LHC16k","Group1_LHC16o","Group1_LHC16p","Group1_LHC17i","Group1_LHC17k","Group1_LHC17l"};
  //  Char_t *arrayOfPeriods[] = {"Group1_LHC16h"};
    int numberOfPeriods = sizeof(arrayOfPeriods) / sizeof(arrayOfPeriods[0]);
    
  //  const double binsCent[6] = {0,1,10,20,40,100};
    Char_t fileInLoc[200];
    Char_t FitFileName[200];
    
    sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFilesTest.root");

// *************************
// Initialiser les graphes *
// *************************
    
    TH2F* hsingletrac(NULL);
    TH1F* hnseg(NULL);
    TH1F* hnseg2(NULL);
    TH1F* hnseg3(NULL);
    TH1F* hnseg4(NULL);
    TH1F* hnseg5(NULL);
    TH1F* hnseg6(NULL);
    TH1F* hnseg7(NULL);
    TH2F* hnseg8(NULL);
    TH2F* hnsegSigma(NULL);
    TH1F* hnseg8_proj_tampon(NULL);
    TH1F* hnseg8Tkl_proj_tampon(NULL);
    TH1F* hnseg8_proj(NULL);
    TH1F* hnseg8Tkl_proj(NULL);
    TH2F* hnseg8ME(NULL);
    TH1F* hnseg8ME_proj_tampon(NULL);
    TH1F* hnseg8TklME_proj_tampon(NULL);
    TH1D* hnseg8ME_proj(NULL);
    TH1D* hnseg8TklME_proj(NULL);
    TH1F* ME_proj_tampon(NULL);
    TH1F* ProjCopy(NULL);
    TH1F* ProjCopy2(NULL);
    TH1F* ProjCopyTkl(NULL);
    TH1F* ProjCopy2Tkl(NULL);
    TH1F* ME_proj_Tkl_tampon(NULL);
    TH1F* ME_proj(NULL);
    TH1F* SE_proj_tampon(NULL);
    TH1F* SE_proj_Tkl_tampon(NULL);
    TH1F* SE_proj(NULL);
    
    TH1F* Sijk(NULL);
    TH1F* Mijk(NULL);
    TH1F* SoverMijk(NULL);
    TH1F* SoverMik(NULL);
    
    TH1F* Sij(NULL);
    TH1F* Mij(NULL);
    TH1F* SoverMij(NULL);
    TH1F* SoverMi(NULL);
    
    TH1F* Yields[NbBinsCent][NbinsInvMass]{ NULL };
    TH1F* YieldsTkl[NbBinsCent]{ NULL };
    TH2F* Correlations[NbBinsCent][NbinsInvMass]{ NULL };
    TH2F* CorrelationsTkl[NbBinsCent]{ NULL };
    TH2I* CorrelationsME[NbBinsCent][NbinsInvMass]{ NULL };
    TH2I* CorrelationsMEMassSummed[NbBinsCent]{ NULL };
    TH2I* CorrelationsTklME[NbBinsCent]{ NULL };
    TH2F* CorrelationsMEScaled[NbBinsCent][NbinsInvMass]{ NULL };
    TH2F* CorrelationsMEMassSummedScaled[NbBinsCent]{ NULL };
    TH2F* CorrelationsTklMEScaled[NbBinsCent]{ NULL };
    TH1F* Yield_tampon(NULL);
    TH1F* YieldTkl_tampon(NULL);
    TH1F* YieldWrtMass_tampon(NULL);
    TH1F* Yield_allC(NULL);
    TH1F* Yield_Central(NULL);
    TH1F* Yield_Periph(NULL);
    TH1F* Yield_Difference(NULL);
    TH1F* YieldTkl_allC(NULL);
    TH1F* YieldTkl_Central(NULL);
    TH1F* YieldTkl_Periph(NULL);
    TH1F* YieldTkl_Difference(NULL);
    TH1F* Yield_Central_MassBin[NbinsInvMass] = { NULL };
    TH1F* Yield_Periph_MassBin[NbinsInvMass] = { NULL };
    TH1F* Yield_Difference_MassBin[NbinsInvMass] = { NULL };
    TH1F* Yields_PhiBin[NbBinsCent][NbinsDeltaPhi]{ NULL };
    TH1F* YieldWrtMass_allC[NbinsDeltaPhi]{ NULL };
    TH1F* YieldWrtMass_Central[NbinsDeltaPhi]{ NULL };
    TH1F* YieldWrtMass_Periph[NbinsDeltaPhi]{ NULL };
    
    // Try plots 2D DetaDphi
    
    TH2F* YCentral(NULL);
    TH2F* YPeriph(NULL);
    TH2F* YDifference(NULL);
    TH2F* YTklCentral(NULL);
    TH2F* YTklCentralME(NULL);
    TH2F* YTklPeriph(NULL);
    TH2F* YTklPeriphME(NULL);
    TH2F* YTklDifference(NULL);
    TH1D* YTklCentral_proj_tampon(NULL);
    TH1D* Correl_tampon(NULL);
    TH1D* YTklPeriph_proj_tampon(NULL);
    TH1D* YTklDifference_proj_tampon(NULL);
    TH1D* YTklCentralME_proj_tampon(NULL);
    TH1D* YTklPeriphME_proj_tampon(NULL);
    
    
    hsingletrac = new TH2F("hsingletrac",
                      "Single tracklet EtaPhi distribution",
                           NbinsDeltaPhi,0,MaxDeltaPhi-MinDeltaPhi,40,-2,2);
    hsingletrac->SetXTitle("Tkl Phi (rad)");
    hsingletrac->SetYTitle("Tkl Eta");
    
    YCentral = new TH2F("YCentral",
                      "Yield delta eta wrt delta phi - Central",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    YCentral->SetXTitle("Correlation DeltaPhi (rad)");
    YCentral->SetYTitle("Correlation DeltaEta");
    YPeriph = new TH2F("YPeriph",
                      "Yield delta eta wrt delta phi - Periph",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    YPeriph->SetXTitle("Correlation DeltaPhi (rad)");
    YPeriph->SetYTitle("Correlation DeltaEta");
    YDifference = new TH2F("YDifference",
                      "Yield delta eta wrt delta phi - Difference",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    YDifference->SetXTitle("Correlation DeltaPhi (rad)");
    YDifference->SetYTitle("Correlation DeltaEta");
    YTklCentral = new TH2F("YTklCentral",
                      "Yield delta eta wrt delta phi - Tkl Central",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklCentral->SetXTitle("Correlation DeltaPhi (rad)");
    YTklCentral->SetYTitle("Correlation DeltaEta");
    YTklPeriph = new TH2F("YTklPeriph",
                      "Yield delta eta wrt delta phi - Tkl Periph",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklPeriph->SetXTitle("Correlation DeltaPhi (rad)");
    YTklPeriph->SetYTitle("Correlation DeltaEta");
    YTklDifference = new TH2F("YTklDifference",
                      "Yield delta eta wrt delta phi - Tkl Difference",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklDifference->SetXTitle("Correlation DeltaPhi (rad)");
    YTklDifference->SetYTitle("Correlation DeltaEta");
    YTklCentralME = new TH2F("YTklCentralME",
                      "Yield delta eta wrt delta phi - Tkl Central",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklCentralME->SetXTitle("Correlation DeltaPhi (rad)");
    YTklCentralME->SetYTitle("Correlation DeltaEta");
    YTklPeriphME = new TH2F("YTklPeriphME",
                      "Yield delta eta wrt delta phi - Tkl Periph",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    YTklPeriphME->SetXTitle("Correlation DeltaPhi (rad)");
    YTklPeriphME->SetYTitle("Correlation DeltaEta");
    // AE
    
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
    
    TH1F* Yields_Central_1(NULL);
    TH1F* Yields_Periph_1(NULL);
    TH1F* Yields_Difference_1(NULL);
    
    TH1F* coefficients0(NULL);
    TH1F* coefficients1(NULL);
    TH1F* coefficients2(NULL);
    TH1F* baselines0(NULL);
    TH1F* c2b0(NULL);
    TH1F* V2JPsiTkl(NULL);
    
    ProjCopy = new TH1F("ProjCopy",
                      "ProjCopy",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    ProjCopy->SetXTitle("Correlation DeltaPhi (rad)");
    ProjCopy->SetYTitle("Correlation DeltaEta");
    ProjCopy2 = new TH1F("ProjCopy2",
                      "ProjCopy2",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    ProjCopy2->SetXTitle("Correlation DeltaPhi (rad)");
    ProjCopy2->SetYTitle("Correlation DeltaEta");
    
    ProjCopyTkl = new TH1F("ProjCopyTkl",
                      "ProjCopyTkl",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    ProjCopyTkl->SetXTitle("Correlation DeltaPhi (rad)");
    ProjCopyTkl->SetYTitle("Correlation DeltaEta");
    ProjCopy2Tkl = new TH1F("ProjCopy2Tkl",
                      "ProjCopy2Tkl",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    ProjCopy2Tkl->SetXTitle("Correlation DeltaPhi (rad)");
    ProjCopy2Tkl->SetYTitle("Correlation DeltaEta");
    
    
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
       sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
       Yield_Central_MassBin[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
       }
       
       for (int j=0; j <NbinsInvMass; j++){
       sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
       Yield_Periph_MassBin[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
       }
       
       for (int j=0; j <NbinsInvMass; j++){
       sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference",MinInvMass+SizeBinInvMass*j,MinInvMass+SizeBinInvMass*(j+1));
       Yield_Difference_MassBin[j] = new TH1F(hname, hname,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
       }
    
    TH1F* InvMass_Central(NULL);
    TH1F* InvMass_Periph(NULL);
    
//    TH1F* YieldWrtMass_PhiBin_Central[NbinsDeltaPhi] = { NULL };
//    TH1F* YieldWrtMass_PhiBin_Periph[NbinsDeltaPhi] = { NULL };
//
//    for (int j=0; j<NbinsDeltaPhi; j++){
//        sprintf(hname,"Yield wrt Mass, Phi: %d pi/6 to %d pi/6 - Central",j-3,j-2);
//        sprintf(hname2,"YieldWrtMassPhiBin_Central_%d",j);
//        YieldWrtMass_PhiBin_Central[j] = new TH1F(hname2, hname,NbinsInvMass,MinInvMass,MaxInvMass);
//    }
//
//    for (int j=0; j<NbinsDeltaPhi; j++){
//        sprintf(hname,"Yield wrt Mass, Phi: %d pi/6 to %d pi/6 - Periph",j-3,j-2);
//        sprintf(hname2,"YieldWrtMassPhiBin_Periph_%d",j);
//        YieldWrtMass_PhiBin_Periph[j] = new TH1F(hname2, hname,NbinsInvMass,MinInvMass,MaxInvMass);
//    }
    
    
    for(int i=0; i<NbBinsCent; i++){
        for(int k=0; k<NbinsInvMass; k++){
            char hname[100];
            sprintf(hname,"Yields %d %d ",i,k);
            Yields[i][k] = new TH1F(hname,
                              "Yields Correlation wrt delta phi",
                              NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
            Yields[i][k]->SetXTitle("Correlation DeltaPhi (rad)");
                sprintf(hname,"Correlations %d %d ",i,k);
                Correlations[i][k] = new TH2F(hname,
                                  "Correlation delta eta wrt delta phi",
                                  NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                Correlations[i][k]->SetXTitle("Correlation DeltaPhi (rad)");
                Correlations[i][k]->SetYTitle("Correlation DeltaEta");
                sprintf(hname,"CorrelationsME %d %d ",i,k);
                CorrelationsME[i][k] = new TH2I(hname,
                                  "MixedEvent Correlation delta eta wrt delta phi",
                                  NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                CorrelationsME[i][k]->SetXTitle("Correlation DeltaPhi (rad)");
                CorrelationsME[i][k]->SetYTitle("Correlation DeltaEta");
                sprintf(hname,"CorrelationsME Scaled %d %d ",i,k);
                CorrelationsMEScaled[i][k] = new TH2F(hname,
                                  "MixedEvent Correlation delta eta wrt delta phi - Scaled",
                                  NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
                CorrelationsMEScaled[i][k]->SetXTitle("Correlation DeltaPhi (rad)");
                CorrelationsMEScaled[i][k]->SetYTitle("Correlation DeltaEta");
            
        }
        char hname2[100];
        sprintf(hname2,"CorrelationsMEMassSummed %d ",i);
        CorrelationsMEMassSummed[i] = new TH2I(hname2,
                          "MixedEvent Correlation delta eta wrt delta phi - All Masses",
                          NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
        CorrelationsMEMassSummed[i]->SetXTitle("Correlation DeltaPhi (rad)");
        CorrelationsMEMassSummed[i]->SetYTitle("Correlation DeltaEta");
    }
    
    for(int i=0; i<NbBinsCent; i++){
        char hname[100];
        sprintf(hname,"Yields Tkl %d ",i);
        YieldsTkl[i] = new TH1F(hname,
                          "Yields Correlation Tkl wrt delta phi",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
        YieldsTkl[i]->SetXTitle("Correlation DeltaPhi (rad)");
            sprintf(hname,"Correlations Tkl %d ",i);
            CorrelationsTkl[i] = new TH2F(hname,
                              "Correlation Tkl delta eta wrt delta phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTkl[i]->SetXTitle("Correlation DeltaPhi (rad)");
            CorrelationsTkl[i]->SetYTitle("Correlation DeltaEta");
            sprintf(hname,"CorrelationsME Tkl %d ",i);
            CorrelationsTklME[i] = new TH2I(hname,
                              "MixedEvent Correlation Tkl delta eta wrt delta phi",
                              NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
            CorrelationsTklME[i]->SetXTitle("Correlation DeltaPhi (rad)");
            CorrelationsTklME[i]->SetYTitle("Correlation DeltaEta");
        sprintf(hname,"CorrelationsME Tkl Scaled %d ",i);
        CorrelationsTklMEScaled[i] = new TH2F(hname,
                          "MixedEvent Correlation Tkl delta eta wrt delta phi",
                          NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
        CorrelationsTklMEScaled[i]->SetXTitle("Correlation DeltaPhi (rad)");
        CorrelationsTklMEScaled[i]->SetYTitle("Correlation DeltaEta");
    
    }
    
    for(int i=0; i<NbBinsCent; i++){
           for(int p=0; p<NbinsDeltaPhi; p++){
               char hname[100];
               sprintf(hname,"Yields PhiBin %d %d ",i,p);
               Yields_PhiBin[i][p] = new TH1F(hname,
                                 "Yields Correlation wrt mass",
                                 NbinsInvMass,MinInvMass,MaxInvMass);
               Yields_PhiBin[i][p]->SetXTitle("Correlation Inv Mass (GeV)");
       }
    }
        
    for(int p=0; p<NbinsDeltaPhi; p++){
            char hname[100];
            char hname2[100];
            sprintf(hname,"Yields PhiBin %d All C",p);
            YieldWrtMass_allC[p] = new TH1F(hname,
                              "Yields Correlation wrt mass, all C",
                              NbinsInvMass,MinInvMass,MaxInvMass);
            YieldWrtMass_allC[p]->SetXTitle("Correlation Inv Mass (GeV)");
        
        sprintf(hname,"Yields PhiBin %d Periph ",p);
        sprintf(hname2,"Yields Correlation wrt mass, Periph, Phi: %d pi/6 to %d pi/6",p-3, p-2);
        YieldWrtMass_Periph[p] = new TH1F(hname,
                          hname2,
                          NbinsInvMass,MinInvMass,MaxInvMass);
        YieldWrtMass_Periph[p]->SetXTitle("Correlation Inv Mass (GeV)");
        
        sprintf(hname,"Yields PhiBin %d Central",p);
        sprintf(hname2,"Yields Correlation wrt mass, Central, Phi: %d pi/6 to %d pi/6",p-3, p-2);
        YieldWrtMass_Central[p] = new TH1F(hname,
                          hname2,
                          NbinsInvMass,MinInvMass,MaxInvMass);
        YieldWrtMass_Central[p]->SetXTitle("Correlation Inv Mass (GeV)");
    }
        
    
    Yield_allC = new TH1F("Yield_allC",
                         "Yields Correlation wrt delta phi, all C, all mass",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
       Yield_allC->SetXTitle("Correlation DeltaPhi (rad)");
    Yield_Central = new TH1F("Yield_Central",
                      "Yields Correlation wrt delta phi, Central, 3.0-3.25 GeV",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yield_Central->SetXTitle("Correlation DeltaPhi (rad)");
    Yield_Periph = new TH1F("Yield_Periph",
                      "Yields Correlation wrt delta phi, Periph, 3.0-3.25 GeV",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yield_Periph->SetXTitle("Correlation DeltaPhi (rad)");
    Yield_Difference = new TH1F("Yield_Difference",
                      "Yields Correlation wrt delta phi, Central-Periph, all mass",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yield_Difference->SetXTitle("Correlation DeltaPhi (rad)");
       Yield_tampon = new TH1F("Yield_tampon",
                         "Yields Correlation wrt delta phi, all C, all mass",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
       Yield_tampon->SetXTitle("Correlation DeltaPhi (rad)");
    YieldWrtMass_tampon = new TH1F("YieldWrtMass_tampon",
                      "Yields Correlation wrt inv mass",
                      NbinsInvMass,MinInvMass,MaxInvMass);
    YieldWrtMass_tampon->SetXTitle("Correlation Inv Mass (GeV)");
    
    YieldTkl_allC = new TH1F("YieldTkl_allC",
                         "Yields Correlation Tkl wrt delta phi, all C, all mass",
                         NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
       YieldTkl_allC->SetXTitle("Correlation DeltaPhi (rad)");
    YieldTkl_Central = new TH1F("YieldTkl_Central",
                      "Yields Correlation Tkl wrt delta phi, Central",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_Central->SetXTitle("Correlation DeltaPhi (rad)");
    YieldTkl_Periph = new TH1F("YieldTkl_Periph",
                      "Yields Correlation Tkl wrt delta phi, Periph",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_Periph->SetXTitle("Correlation DeltaPhi (rad)");
    YieldTkl_Difference = new TH1F("YieldTkl_Difference",
                      "Yields Correlation Tkl wrt delta phi, Central-Periph, all mass",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_Difference->SetXTitle("Correlation DeltaPhi (rad)");
    YieldTkl_tampon = new TH1F("YieldTkl_tampon",
                      "Yields Correlation Tkl wrt delta phi, all C, all mass",
                      NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    YieldTkl_tampon->SetXTitle("Correlation DeltaPhi (rad)");
    
    
    hnseg = new TH1F("hnseg",
                     "Invariant mass of dimuon",
                     250,MinInvMass,MaxInvMass);
    hnseg->SetXTitle("Mass of dimuon (GeV/c^{2})");
    hnseg->SetYTitle("Count");
    hDPhi = new TH1F("hDPhi",
                     "DeltaPhi of Tracklet",
                     10000,-TMath::Pi(),TMath::Pi());
    hDPhi->SetXTitle("DPhi tracklet");
    hDPhi->SetYTitle("Count");
    InvMass_Central = new TH1F("InvMass_Central",
                     "Invariant mass of dimuon - Central",
                     250,MinInvMass,MaxInvMass);
    InvMass_Central->SetXTitle("Mass of dimuon (GeV/c^{2})");
    InvMass_Central->SetYTitle("Count");
    InvMass_Periph = new TH1F("InvMass_Periph",
                     "Invariant mass of dimuon - Periph",
                     250,MinInvMass,MaxInvMass);
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
    hnsegSigma = new TH2F("hnsegSigma",
                      "Zvtx resolution wrt NContributors",
                      60,0,60,1000,0,10);
    hnsegSigma->SetXTitle("NContributors");
    hnsegSigma->SetYTitle("Zvtx Resolution (SPD)");
    hnseg8 = new TH2F("hnseg8",
                      "Correlation delta eta wrt delta phi",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8->SetYTitle("Correlation DeltaEta");
    hnseg8ME = new TH2F("hnseg8ME",
                      "MixedEvent Correlation delta eta wrt delta phi",
                      NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi,NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8ME->SetXTitle("Correlation DeltaPhi (rad)");
    hnseg8ME->SetYTitle("Correlation DeltaEta");
    hnseg8_proj = new TH1F("hnseg8_proj", "Correlation wrt delta eta, projection",NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8_proj->SetXTitle("Correlation DeltaEta");
    hnseg8_proj->SetYTitle("Correlation mean");
    hnseg8Tkl_proj = new TH1F("hnseg8Tkl_proj", "Correlation Tkl wrt delta eta, projection",NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8Tkl_proj->SetXTitle("Correlation DeltaEta");
    hnseg8Tkl_proj->SetYTitle("Correlation mean");
    hnseg8ME_proj = new TH1D("hnseg8ME_proj", "MixedEvent Correlation wrt delta eta, projection",NbinsDeltaEta,MinDeltaEta,MaxDeltaEta);
    hnseg8ME_proj->SetXTitle("Correlation DeltaEta");
    hnseg8ME_proj->SetYTitle("Correlation mean");
    hnseg8TklME_proj = new TH1D("hnseg8TklME_proj", "MixedEvent Correlation Tkl wrt delta eta, projection",NbinsDeltaEtaTKL,MinDeltaEtaTKL,MaxDeltaEtaTKL);
    hnseg8TklME_proj->SetXTitle("Correlation DeltaEta");
    hnseg8TklME_proj->SetYTitle("Correlation mean");
 
    Sijk = new TH1F("Sijk",
                                   "Sijk",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Sijk->SetXTitle("Correlation DeltaPhi (rad)");
    Mijk = new TH1F("Mijk",
                                   "Mijk",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Mijk->SetXTitle("Correlation DeltaPhi (rad)");
    SoverMijk = new TH1F("SoverMijk",
                                   "SoverMijk",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    SoverMijk->SetXTitle("Correlation DeltaPhi (rad)");
    SoverMik = new TH1F("SoverMik",
                                   "SoverMik",
                                   NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    SoverMik->SetXTitle("Correlation DeltaPhi (rad)");
    
    Sij = new TH1F("Sij",
                                   "Sij",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    Sij->SetXTitle("Correlation DeltaPhi (rad)");
    Mij = new TH1F("Mij",
                                   "Mij",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    Mij->SetXTitle("Correlation DeltaPhi (rad)");
    SoverMij = new TH1F("SoverMij",
                                   "SoverMij",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    SoverMij->SetXTitle("Correlation DeltaPhi (rad)");
    SoverMi = new TH1F("SoverMi",
                                   "SoverMi",
                                   NbinsDeltaPhiTKL,MinDeltaPhi,MaxDeltaPhi);
    SoverMi->SetXTitle("Correlation DeltaPhi (rad)");

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
                      250,MinInvMass,MaxInvMass,6,binsPt);
    hPtWrtMassInv[0]->SetXTitle("Mass Inv");
    hPtWrtMassInv[0]->SetYTitle("Pt");
    hPtWrtMassInv[1] = new TH2F("hPtWrtMassInv_1",
                      "Pt wrt Mass Inv - Central",
                      250,MinInvMass,MaxInvMass,6,binsPt);
    hPtWrtMassInv[1]->SetXTitle("Mass Inv");
    hPtWrtMassInv[1]->SetYTitle("Pt");
    hPtWrtMassInv[2] = new TH2F("hPtWrtMassInv_2",
                      "Pt wrt Mass Inv - Periph",
                      250,MinInvMass,MaxInvMass,6,binsPt);
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
    
    
    int DimuCentralSeen = 0;
    int DimuPeriphSeen = 0;
    int DimuSeenMassCut = 0;
    int DimuSeenNoMassCut = 0;
    int EventRejected = 0;
    int EventNC = 0;
    int EventPileUpMult = 0;
    int EventPileUpVtx = 0;
    int RefTracklets = 0;
    int RefTrackletsCentral = 0;
    int RefTrackletsPeriph = 0;
    int cmul = 0;
    int barWidth = 50;
    int DimuonCounter[NbBinsCent][NbinsInvMass] = {0};
    int DimuonCounterZint[NbBinsCent][NbinsInvMass] = {0};
    int RefTklCounter[NbBinsCent] = {0};
    int RefTklCounterZint[NbBinsCent] = {0};
    double NormME[NbBinsCent][NbinsInvMass] = {0};
    double NormMETkl[NbBinsCent] = {0};
    
    int DimuC = 0;
    int DimuP = 0;
    int TklC = 0;
    int TklP = 0;
    int NormTklCentral = 0;
    int NormTklPeriph = 0;
    int countsigma = 0;
    
        
//    TFile fileIn1("~/../../Volumes/Transcend2/ppAnalysis/Scripts/merge16All17hikl_PS_CutsEvent.root");
//    TFile fileIn2("~/../../Volumes/Transcend2/ppAnalysis/Scripts/merge17mor_PS_CutsEvent.root");
//    TFile fileIn3("~/../../Volumes/Transcend2/ppAnalysis/Scripts/merge18All_PS_CutsEvent.root");
//    TFile fileIn1("~/../../Volumes/Transcend2/ppAnalysis/Scripts/GoodLHC16h_pass1_PS_CutsEvent/muonGrid.root");
//    TFile fileIn2("~/../../Volumes/Transcend2/ppAnalysis/Scripts/GoodLHC16h_pass1_PS_CutsEvent/muonGrid.root");
//    TFile fileIn3("~/../../Volumes/Transcend2/ppAnalysis/Scripts/GoodLHC16h_pass1_PS_CutsEvent/muonGrid.root");
//               TFile fileIn1("~/../../Volumes/Transcend2/ppAnalysis/Scripts/Group1_LHC16h/muonGrid.root");
//               TFile fileIn2("~/../../Volumes/Transcend2/ppAnalysis/Scripts/Group1_LHC16k/muonGrid.root");
//               TFile fileIn3("~/../../Volumes/Transcend2/ppAnalysis/Scripts/Group1_LHC16h/muonGrid.root");
  //  TFile fileIn1("~/../../Volumes/Transcend2/ppAnalysis/Scripts/GoodLHC17r_muoncalopass1_PS_CutsEvent/muonGrid.root");
//    TFile fileIn2("~/../../Volumes/Transcend2/ppAnalysis/Scripts/GoodLHC17r_muoncalopass1_PS_CutsEvent/muonGrid.root");
    
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
    
//    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
//
//           sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
//       TFile fileIn(fileInLoc);
//
//        char str [20];
//        int GroupNum;
//
//        sscanf (arrayOfPeriods[tree_idx],"Group%d_%s",&GroupNum,str);
//
//
//
//    std::cout << "1" <<std::endl;
//        fileIn.GetObject("MyMuonTree",theTree);
//    std::cout << "2" <<std::endl;
//        //setting the branch address
//    std::cout << "3" <<std::endl;
//  //  theTree->SetBranchAddress("event", &fEvent);
//    TBranch *branch = theTree->GetBranch("event");
////        auto bntrack = theTree->GetBranch("tracks");
//        //auto branch  = theTree->GetBranch("mcparticles.fPx");
//    std::cout << "3.5" <<std::endl;
//        branch->SetAddress(&fEvent);
//    std::cout << "4" <<std::endl;
//    auto nevent = theTree->GetEntries();
//
//    std::cout << "5" <<std::endl;
//    fCorrelations = fEvent->fCorrelations;
//    fTracklets = fEvent->fTracklets;
//    fDimuons = fEvent->fDimuons;
////            auto nevent = bntrack->GetEntries();
//    //cout << " theTree->GetEntries() " << theTree->GetEntries() << endl;
//
//
//// ************************************
//// Boucle sur les evenements pour pools mixed-event, notés i *
//// ************************************
//    if(doMixedEvents){
//        cout << "MIXED EVENTS IS WANTED - WILL NOW PROCESS EVENTS AND ORGANISE POOLS" << endl;
//      //  Double_t px[1000], py[1000];
//            for (int i=0;i<nevent;i++) {
//                if(i%100000 == 0){
//                        std::cout << "[";
//                        double portion = double(i)/nevent;
//                        long pos = barWidth * portion;
//                        for (int k = 0; k < barWidth; ++k) {
//                            if (k < pos) std::cout << "=";
//                            else if (k == pos) std::cout << ">";
//                            else std::cout << " ";
//                        }
//                        std::cout << "] " << long(100 * portion) << "%     " << i << "/" << nevent << " Tree " << tree_idx+1 << "/" << numberOfPeriods << " Pooling";
//                        std::cout.flush();
//                    std::cout << std::endl;
//                }
////                if(doTracklets && (i%10!=0)){
////                    continue;
////                }
//                    theTree->GetEvent(i);
//                cout << "=== EVENT " << i << " ==="<<endl;
//                if(fEvent->fIsPileupFromSPDMultBins || fEvent->fNPileupVtx != 0 || fEvent->fVertexNC < 1){
//                    cout << "Rejected Pileup"<<endl;
//                    continue;
//                }
//                if(floor(fEvent->fVertexZ) + ZvtxCut != ZvtxBin){
//                    cout << "Rejected Zcut"<<endl;
//                    continue;
//                }
//
//                int NumberOfTrackletsPassingEtaCut = 0;
//                for (Int_t j=0; j<fEvent->fNTracklets; j++) {
//                    trackletME = (TrackletLight*)fTracklets->At(j);
//                    if(TMath::Abs(trackletME->fEta) < TklEtaCut && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
//                        NumberOfTrackletsPassingEtaCut++;
//                    }
//
//                }
//
//                if(NumberOfTrackletsPassingEtaCut==0){
//                    cout << "Rejected Tkl Cuts"<<endl;
//                    continue;
//                }
//
//                if(KeepOnlyOne){
//                    if(NumberOfTrackletsPassingEtaCut != valueOnlyOne){
//                        continue;
//                    }
//                }
//
//                if(AdditionalCutNtkl){
//                   if(NumberOfTrackletsPassingEtaCut <= valueAdditionalCutNtkl){
//                       cout << "Rejected Additional tkl cut"<<endl;
//                       continue;
//                   }
//               }
//
//                for (Int_t j=0; j<fEvent->fNDimuons; j++) {
//                   dimu = (DimuonLight*)fDimuons->At(j);
//                   if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut) && (dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut) && (dimu->fInvMass > MinInvMass) && (dimu->fInvMass < MaxInvMass)){
//
//                       double cent = fEvent->fCentralitySPDTracklets;
//                       int zv = floor(fEvent->fVertexZ) + ZvtxCut;
//              //         int centint = GetCent(cent);
//                       int centint = GetCentPM(NumberOfTrackletsPassingEtaCut, zv, GroupNum);
//                       int phint = 0;
//                       double phi = dimu->fPhi;
//                       if(phi < -TMath::Pi()/2){
//                        phi += 2* TMath::Pi();
//                          }
//                          if(phi > 1.5*TMath::Pi()){
//                              phi -= 2* TMath::Pi();
//                          }
//                       phi = phi*6/TMath::Pi();
//                       phint = floor(phi) + 3;
//                       if(PoolsSize[centint]<100){ //[phint]
//                           for (Int_t j=0; j<fEvent->fNTracklets; j++) {
//                               trackletME = (TrackletLight*)fTracklets->At(j);
//                               if(TMath::Abs(trackletME->fEta) < TklEtaCut && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
//                                   double trackletMEPhi = (Float_t)trackletME->fPhi;
//                                   double trackletMEEta = trackletME->fEta;
//                                   Pools[centint].push_back(trackletMEPhi); //[phint]
//                                   Pools[centint].push_back(trackletMEEta); //[phint]
//                               }
//
//                           }
//                           PoolsSize[centint] += 1; //[phint]
//                       }
//                   }
//                }
//
//                if(doTracklets){
//
//                   if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut)){
//                       double cent = fEvent->fCentralitySPDTracklets;
//                       int zv = floor(fEvent->fVertexZ) + ZvtxCut;
//                   //    int centint = GetCent(cent);
//                       int centint = GetCentPM(NumberOfTrackletsPassingEtaCut, zv, GroupNum);
//                       if(PoolsSizeTkl[centint]<100){
//                           for (Int_t j=0; j<fEvent->fNTracklets; j++) {
//                               trackletME = (TrackletLight*)fTracklets->At(j);
//                               if(TMath::Abs(trackletME->fEta) < TklEtaCut && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
//                                   double trackletMEPhi = (Float_t)trackletME->fPhi;
//                                   double trackletMEEta = trackletME->fEta;
//                                   PoolsTkl[centint].push_back(trackletMEPhi);
//                                   PoolsTkl[centint].push_back(trackletMEEta);
//                                   PoolsTklEventTracker[centint].push_back(i);
//                                   PoolsTklEventTracker[centint].push_back(i);
//                                       cout << "i " << i << " trackletMEPhi " << trackletMEPhi << " trackletMEEta " <<trackletMEEta<<endl;
//
//                               }
//
//                           }
//                           PoolsSizeTkl[centint] += 1;
//                       }
//                   }
//
//                }
//
//
//            }
//    }
//
//    }
    
// ************************************
// Boucle sur les evenements, notés i *
// ************************************
    
    cout << "=== NOW EVENT STUDY ===" <<endl;
    cout << "Switching back to 1st TTree" <<endl;
    
    for(int tree_idx=0; tree_idx<numberOfPeriods; tree_idx++){
        
        sprintf(fileInLoc,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/%s/muonGrid.root",arrayOfPeriods[tree_idx]);
    TFile fileIn(fileInLoc);
        
        char str [20];
        int GroupNum;
        
        std::vector <double> Pools[NbBinsCent]; //[12]
        std::vector <int> PoolsEventTracker[NbBinsCent];
        int PoolsSize[NbBinsCent] = {0}; //[12]
        std::vector <double> PoolsTkl[NbBinsCent];
        std::vector <int> PoolsTklEventTracker[NbBinsCent];
        int PoolsSizeTkl[NbBinsCent] = {0};

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
        fCorrelations = fEvent->fCorrelations;
        fTracklets = fEvent->fTracklets;
        fDimuons = fEvent->fDimuons;
    
    branch = theTree->GetBranch("event");
    branch->SetAddress(&fEvent);
    fCorrelations = fEvent->fCorrelations;
    fTracklets = fEvent->fTracklets;
    fDimuons = fEvent->fDimuons;
        int DimuMEcounter =0;
        int TklMEcounter =0;
    
        for (int i=0;i<nevent;i++) {
            if(i%100000 == 0){
                    std::cout << "[";
                    double portion = double(i)/nevent;
                    long pos = barWidth * portion;
                    for (int k = 0; k < barWidth; ++k) {
                        if (k < pos) std::cout << "=";
                        else if (k == pos) std::cout << ">";
                        else std::cout << " ";
                    }
                    std::cout << "] " << long(100 * portion) << "%     " << i << "/" << nevent << " Tree " << tree_idx+1 << "/" << numberOfPeriods;
                    std::cout.flush();
                std::cout << std::endl;
            }
//            if(doTracklets && (i%10!=0)){
//                continue;
//            }
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
            if(floor(fEvent->fVertexZ) + ZvtxCut != ZvtxBin){
                continue;
            }
            
            int NumberCloseEtaTracklets = 0;
            int NumberOfTrackletsPassingEtaCut = 0;
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                trac = (TrackletLight*)fTracklets->At(j);
                hsingletrac->Fill(trac->fPhi, trac->fEta);
                if((TMath::Abs(trac->fEta) < TklEtaCut && (TMath::Abs(trac->fDPhi) < DPhiCut))){
                    NumberOfTrackletsPassingEtaCut++;
                    hnseg9->Fill(fEvent->fVertexZ, trac->fEta);
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
            
            CentV0M->Fill(fEvent->fCentralityV0M);
//            CentTKL->Fill(fEvent->fCentralityTKL);
//            CentCL0->Fill(fEvent->fCentralityCL0);
//            CentCL1->Fill(fEvent->fCentralityCL1);
            CentSPDTracklets->Fill(fEvent->fCentralitySPDTracklets);
//            CentSPDClusters->Fill(fEvent->fCentralitySPDClusters);
            
            CentV0Mwrttkl->Fill(fEvent->fCentralityV0M, NumberCloseEtaTracklets);
//            CentTKLwrttkl->Fill(fEvent->fCentralityTKL, NumberCloseEtaTracklets);
//            CentCL0wrttkl->Fill(fEvent->fCentralityCL0, NumberCloseEtaTracklets);
//            CentCL1wrttkl->Fill(fEvent->fCentralityCL1, NumberCloseEtaTracklets);
            CentSPDTrackletswrttkl->Fill(fEvent->fCentralitySPDTracklets, NumberCloseEtaTracklets);
//            CentSPDClusterswrttkl->Fill(fEvent->fCentralitySPDClusters, NumberCloseEtaTracklets);

//            hnsegSigma->Fill(fEvent->fVertexNC, fEvent->fSPDVertexSigmaZ); //SIGMA
//
//            if(sqrt(fEvent->fSPDVertexSigmaZ)>0.25){
//                countsigma++;
//            }
            
            for (Int_t j=0; j<fEvent->fNDimuons; j++) {
                dimu = (DimuonLight*)fDimuons->At(j);
                if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut) && (dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)){
                    //(fEvent->fNPileupVtx == 0) &&
                    //  (TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (dimu->fEta < HighDimuEtaCut ) && (dimu->fEta > LowDimuEtaCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)
                    // fEvent->fNPileupVtx == 0) && (TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (dimu->fY < HighDimuYCut ) && (dimu->fY > LowDimuYCut) && (dimu->fCharge == 0) && (dimu->fPt > LowDimuPtCut) && (dimu->fPt < HighDimuPtCut)
                    DimuSeenNoMassCut++;
                    
                    int centint = 99;
                    
                    if((dimu->fInvMass > MinInvMass) && (dimu->fInvMass < MaxInvMass)){
                        double mass = dimu->fInvMass;
                        int massint = int((mass-MinInvMass)/SizeBinInvMass);
                       double cent = fEvent->fCentralitySPDTracklets;
                       double zv = fEvent->fVertexZ;
                       int zvint = floor(zv) + ZvtxCut;
                  //      int centint = GetCent(cent);
                        centint = GetCentPM(NumberOfTrackletsPassingEtaCut, zvint, GroupNum);
                        DimuonCounter[centint][massint]++;
                     //   if(cent <= CentSPDTrackletsCentral){
                        if(isCentral(centint)){
                            DimuC++;
                        }
                      //  else if(cent > CentSPDTrackletsPeriph){
                        else if(isPeripheral(centint)){ //CentPeriph
                           DimuP++;
                        }
                    }
                    
                    hnseg->Fill(dimu->fInvMass);
                    hPtWrtMassInv[0]->Fill(dimu->fInvMass, dimu->fPt);
                  //  if(fEvent->fCentralitySPDTracklets<=CentSPDTrackletsCentral){
                    if(isCentral(centint)){
                        InvMass_Central->Fill(dimu->fInvMass);
                        hPtWrtMassInv[1]->Fill(dimu->fInvMass, dimu->fPt);
                    }
                  //  if(fEvent->fCentralitySPDTracklets>=CentSPDTrackletsPeriph){
                    if(isPeripheral(centint)){ //CentPeriph
                        InvMass_Periph->Fill(dimu->fInvMass);
                        hPtWrtMassInv[2]->Fill(dimu->fInvMass, dimu->fPt);
                    }
                                
// Event mixing
                    if(doMixedEvents){
                        if((dimu->fInvMass > MinInvMass) && (dimu->fInvMass < MaxInvMass)){
                        double massME = dimu->fInvMass;
                        int massintME = int((massME-MinInvMass)/SizeBinInvMass);
                     //  double centME = fEvent->fCentralitySPDTracklets;
                     //   int centintME = GetCent(centME);
                       double zvME = fEvent->fVertexZ;
                       int zvintME = floor(zvME) + ZvtxCut;
                        int centintME = GetCentPM(NumberOfTrackletsPassingEtaCut, zvintME, GroupNum);
                       int phintME = 0;
                       double phiME = dimu->fPhi;
                       if(phiME < -TMath::Pi()/2){
                        phiME += 2* TMath::Pi();
                          }
                          if(phiME > 1.5*TMath::Pi()){
                              phiME -= 2* TMath::Pi();
                          }
                       phintME = floor(phiME*6/TMath::Pi()) + 3;
                            DimuMEcounter++;
                            
                            if(PoolsSize[centintME]>=10){

                               for(int k=0; k< Pools[centintME].size(); k+=2){ //[phintME]
                                           double correlMEPhi = Pools[centintME].at(k) - dimu->fPhi; //[phintME]
                             //      cout << ls[centintME][zvintME].at(k) " << k << " : " << Pools[centintME][zvintME].at(k) <<endl; //[phintME]
                                           if(correlMEPhi < -TMath::Pi()/2){
                                               correlMEPhi += 2* TMath::Pi();
                                           }
                                           if(correlMEPhi > 1.5*TMath::Pi()){
                                               correlMEPhi -= 2* TMath::Pi();
                                           }
                                           double correlMEEta = TMath::Abs(dimu->fEta - Pools[centintME].at(k+1)); //[phintME]
                                                CorrelationsME[centintME][massintME]->Fill(correlMEPhi,correlMEEta);
                                                CorrelationsMEMassSummed[centintME]->Fill(correlMEPhi,correlMEEta);
                                    }
                            }
                            if(PoolsSize[centintME] == 100){
                                int valueDiscarded = PoolsEventTracker[centintME].front();
                                while(PoolsEventTracker[centintME].front()==valueDiscarded){
                                    PoolsEventTracker[centintME].erase(PoolsEventTracker[centintME].begin(),PoolsEventTracker[centintME].begin()+2);
                                    Pools[centintME].erase(Pools[centintME].begin(),Pools[centintME].begin()+2);
                                }
                                for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                                   trackletME = (TrackletLight*)fTracklets->At(j);
                                   if(TMath::Abs(trackletME->fEta) < TklEtaCut && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
                                       double trackletMEPhi = (Float_t)trackletME->fPhi;
                                       double trackletMEEta = trackletME->fEta;
                                       Pools[centintME].push_back(trackletMEPhi);
                                       Pools[centintME].push_back(trackletMEEta);
                                       PoolsEventTracker[centintME].push_back(DimuMEcounter);
                                       PoolsEventTracker[centintME].push_back(DimuMEcounter);
                                   }
    
                               }
                            }
                            else if(PoolsSize[centintME] < 100){
                               for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                                   trackletME = (TrackletLight*)fTracklets->At(j);
                                   if(TMath::Abs(trackletME->fEta) < TklEtaCut && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
                                       double trackletMEPhi = (Float_t)trackletME->fPhi;
                                       double trackletMEEta = trackletME->fEta;
                                       Pools[centintME].push_back(trackletMEPhi);
                                       Pools[centintME].push_back(trackletMEEta);
                                       PoolsEventTracker[centintME].push_back(DimuMEcounter);
                                       PoolsEventTracker[centintME].push_back(DimuMEcounter);
                                   }
    
                               }
                               PoolsSize[centintME] += 1;
                           }
                        }
                        
                    }
                                       
                                       
                }
            }
           
            
// ***********************************
// Boucle sur les correlations, notés j    *
// ***********************************
            int counter_one = 0;
            for (Int_t j=0; j<fEvent->fNCorrelations; j++) {
                correl = (CorrelationLight*)fCorrelations->At(j);
                if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut) && (correl->fDimuonY < HighDimuYCut ) && (correl->fDimuonY > LowDimuYCut) && (correl->fDimuonCharge == 0) && (correl->fDimuonPt > LowDimuPtCut) && (correl->fDimuonPt < HighDimuPtCut) && (TMath::Abs(correl->fTrackletEta) < 1) && (TMath::Abs(correl->fTrackletDPhi) < DPhiCut) && (TMath::Abs(correl->fTrackletEta) < TklEtaCut)){   //Cuts    (fEvent->fNPileupVtx == 0) &&
                    counter_one++;
                    Float_t DeltaPhi = correl->fDeltaPhi;
                    if(DeltaPhi < -TMath::Pi()/2){
                        DeltaPhi += 2* TMath::Pi();
                    }
                    if(DeltaPhi > 1.5*TMath::Pi()){
                        DeltaPhi -= 2* TMath::Pi();
                    }
                    if (counter_one==1){
                        hnseg2->Fill(correl->fDimuonPt);
                        hnseg3->Fill(correl->fDimuonEta);
                        hnseg4->Fill(correl->fDimuonY);
                        hnseg5->Fill(correl->fDimuonPhi);
                    }
                    hnseg6->Fill(DeltaPhi);
                    hnseg7->Fill(correl->fDeltaEta);
                    hnseg8->Fill(DeltaPhi, correl->fDeltaEta);
                    
                    if((correl->fDimuonInvMass > MinInvMass) && (correl->fDimuonInvMass < MaxInvMass)){
                    double mass = correl->fDimuonInvMass;
                     int massint = int((mass-MinInvMass)/SizeBinInvMass);
                   // double cent = fEvent->fCentralitySPDTracklets;
                  //  int centint = GetCent(cent);
                    double zv = fEvent->fVertexZ;
                    int zvint = floor(zv) + ZvtxCut;
                    int centint = GetCentPM(NumberOfTrackletsPassingEtaCut, zvint, GroupNum);
                    //    cout << centint << " " << zvint << " " << massint <<endl;
                    Correlations[centint][massint]->Fill(DeltaPhi,correl->fDeltaEta);
                      //  if(cent <= CentSPDTrackletsCentral){
                        if(isCentral(centint)){
                            YCentral->Fill(DeltaPhi,correl->fDeltaEta);
                        }
                      //  else if(cent > CentSPDTrackletsPeriph){
                        else if(isPeripheral(centint)){ //CentPeriph
                            YPeriph->Fill(DeltaPhi,correl->fDeltaEta);
                        }
                    }
                }
            }
            
            if((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut)){
                hnseg10->Fill(fEvent->fCentralityV0M, NumberCloseEtaTracklets); //fEvent->fNTracklets
            }
            
            
            
            
            
// TREATEMENT OF TRACKLET V2
            if(doTracklets){
                
                if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut)){
                 // double cent = fEvent->fCentralitySPDTracklets;
                //  int centint = GetCent(cent);
                  double zv = fEvent->fVertexZ;
                  int zvint = floor(zv) + ZvtxCut;
              //  cout << NumberOfTrackletsPassingEtaCut << " " << zvint << " " <<GroupNum<<endl;
                  int centint = GetCentPM(NumberOfTrackletsPassingEtaCut, zvint, GroupNum);
                   RefTklCounter[centint] += NumberCloseEtaTracklets; //fEvent->fNTracklets
              //  if(cent <= CentSPDTrackletsCentral){
                if(isCentral(centint)){
                    TklC++;
                }
               // else if(cent > CentSPDTrackletsPeriph){
                else if(isPeripheral(centint)){ //CentPeriph
                    TklP++;
                }
                }
            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                fTracklets->Randomize();
                tracklet1 = (TrackletLight*)fTracklets->At(j);
                for (Int_t k=j+1; k<fEvent->fNTracklets; k++){
                    tracklet2 = (TrackletLight*)fTracklets->At(k);
                    if ((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut) && (TMath::Abs(tracklet1->fDPhi) < DPhiCut) && (TMath::Abs(tracklet2->fDPhi) < DPhiCut) && (TMath::Abs(tracklet1->fEta) < TklEtaCut) && (TMath::Abs(tracklet2->fEta) < TklEtaCut)){   //Cuts
                            Float_t DeltaPhi = tracklet1->fPhi - tracklet2->fPhi;
                            if(DeltaPhi < -TMath::Pi()/2){
                                DeltaPhi += 2* TMath::Pi();
                            }
                            if(DeltaPhi > 1.5*TMath::Pi()){
                                DeltaPhi -= 2* TMath::Pi();
                            }
                            Float_t DeltaEta = tracklet1->fEta - tracklet2->fEta; //DeltaEtaAbs TMath::Abs(
                        if(TMath::Abs(DeltaEta)<DeltaEtaTKLCut){
                            continue;
                        }
                            //   double cent = fEvent->fCentralitySPDTracklets;
                          //     int centint = GetCent(cent);
                               double zv = fEvent->fVertexZ;
                               int zvint = floor(zv) + ZvtxCut;
                                int centint = GetCentPM(NumberOfTrackletsPassingEtaCut, zvint, GroupNum);
                               CorrelationsTkl[centint]->Fill(DeltaPhi,DeltaEta);
                              //  if(cent <= CentSPDTrackletsCentral){
                                if(isCentral(centint)){
                                    YTklCentral->Fill(DeltaPhi,DeltaEta);
                                }
                              //  else if(cent > CentSPDTrackletsPeriph){
                                else if(isPeripheral(centint)){ //CentPeriph
                                    YTklPeriph->Fill(DeltaPhi,DeltaEta);
                                }
                   }
                }
            }
                
    // Event mixing
                
                if(doMixedEvents){
                  //  cout << "ME asked" <<endl;
                    if((TMath::Abs(fEvent->fVertexZ) < ZvtxCut) && (TMath::Abs(fEvent->fSPDVertexSigmaZ) < SigmaZvtxCut)){
                        //   double centME = fEvent->fCentralitySPDTracklets;
                         //  int centintME = GetCent(centME);
                           double zvME = fEvent->fVertexZ;
                           int zvintME = floor(zvME) + ZvtxCut;
                             int centintME = GetCentPM(NumberOfTrackletsPassingEtaCut, zvintME, GroupNum);
                        TklMEcounter++;
                        //cout << "centintME " <<centintME<<endl;
                        //cout << " PoolsSizeTkl[centintME] " << PoolsSizeTkl[centintME]<<endl;
                        if(PoolsSizeTkl[centintME]>=10){
                        //    cout << ">=10" <<endl;
                            for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                              tracklet1 = (TrackletLight*)fTracklets->At(j);
                                if(TMath::Abs(tracklet1->fEta) < TklEtaCut && (TMath::Abs(tracklet1->fDPhi) < DPhiCut)){
                                
                                  for(int k=0; k< PoolsTkl[centintME].size(); k+=2){
                                              double correlTklMEPhi = tracklet1->fPhi - PoolsTkl[centintME].at(k);
                                  //    cout << " PoolsTkl[centintME].at(k) " << PoolsTkl[centintME].at(k)<<endl;
                                              if(correlTklMEPhi < -TMath::Pi()/2){
                                                  correlTklMEPhi += 2* TMath::Pi();
                                              }
                                              if(correlTklMEPhi > 1.5*TMath::Pi()){
                                                  correlTklMEPhi -= 2* TMath::Pi();
                                              }
                                      double correlTklMEEta = tracklet1->fEta - PoolsTkl[centintME].at(k+1); //DeltaEtaAbs
                                      if(TMath::Abs(correlTklMEEta)>DeltaEtaTKLCut){
                                            CorrelationsTklME[centintME]->Fill(correlTklMEPhi,correlTklMEEta);
                                        }
                                      if(correlTklMEPhi > -SizeBinDeltaPhiTKL/2 && correlTklMEPhi < SizeBinDeltaPhiTKL/2 && TMath::Abs(correlTklMEEta) < SizeBinDeltaEtaTKL/2){
                                                NormMETkl[centintME] += 1.;
                                            }
                                      
                                      
                                  //        if(cent <= CentSPDTrackletsCentral){
                                      if(isCentral(centintME)){
                                              if(TMath::Abs(correlTklMEEta)>DeltaEtaTKLCut){
                                                  YTklCentralME->Fill(correlTklMEPhi,correlTklMEEta);
                                              }
                                              if(correlTklMEPhi > 0 && correlTklMEPhi < SizeBinDeltaPhiTKL && TMath::Abs(correlTklMEEta) < SizeBinDeltaEtaTKL){
                                                  NormTklCentral += 1;
                                              }
                                          }
                                 //         else if(cent > CentSPDTrackletsPeriph){
                                      else if(isPeripheral(centintME)){ //CentPeriph
                                            if(TMath::Abs(correlTklMEEta)>DeltaEtaTKLCut){
                                                YTklPeriphME->Fill(correlTklMEPhi,correlTklMEEta);
                                            }
                                              if(correlTklMEPhi > 0 && correlTklMEPhi < SizeBinDeltaPhiTKL && TMath::Abs(correlTklMEEta) < SizeBinDeltaEtaTKL){
                                                  NormTklPeriph += 1;
                                              }
                                          }
                                          
                                  }
                                }
                            }
                          //  cout << "Tkl ME has been calculated" << endl;
                        }
                        
                        
                        if(PoolsSizeTkl[centintME] == 100){
                         //   cout << "The Tkl pool is full" <<endl;
                           int valueDiscarded = PoolsTklEventTracker[centintME].front();
                          //  cout << "Discarding events with index " << valueDiscarded<<endl;
                           while(PoolsTklEventTracker[centintME].front()==valueDiscarded){
                         //      cout << "Event number at front: " << PoolsTklEventTracker[centintME].front()<<endl;
                               PoolsTklEventTracker[centintME].erase(PoolsTklEventTracker[centintME].begin(),PoolsTklEventTracker[centintME].begin()+2);
                             //  cout << "Discarded Event tracker front - New front number : "<< PoolsTklEventTracker[centintME].front() << endl;
                            //   cout << "Four first elements in PoolsTkl: " << PoolsTkl[centintME].at(0) << " " << PoolsTkl[centintME].at(1) << " " << PoolsTkl[centintME].at(2) << " " << PoolsTkl[centintME].at(3) << endl;
                             //  cout << "Size " << PoolsTkl[centintME].size()<<endl;
                               PoolsTkl[centintME].erase(PoolsTkl[centintME].begin(),PoolsTkl[centintME].begin()+2);
                             //  cout << "Two first discarded - New Four first elements in PoolsTkl: " << PoolsTkl[centintME].at(0) << " " << PoolsTkl[centintME].at(1) << " " << PoolsTkl[centintME].at(2) << " " << PoolsTkl[centintME].at(3) << endl;
                           }
                          //  cout << "Will now add new event"<<endl;
                           for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                              trackletME = (TrackletLight*)fTracklets->At(j);
                              if(TMath::Abs(trackletME->fEta) < TklEtaCut && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
                                  double trackletMEPhi = (Float_t)trackletME->fPhi;
                                  double trackletMEEta = trackletME->fEta;
                                  PoolsTkl[centintME].push_back(trackletMEPhi);
                                  PoolsTkl[centintME].push_back(trackletMEEta);
                                  PoolsTklEventTracker[centintME].push_back(TklMEcounter);
                                  PoolsTklEventTracker[centintME].push_back(TklMEcounter);
                              }

                          }
                         //   cout << "New event added"<<endl;
                            
                       }
                        
                        
                        else if(PoolsSizeTkl[centintME] < 100){
                         //   cout << "Pool is not full - Will add event"<<endl;
                          for (Int_t j=0; j<fEvent->fNTracklets; j++) {
                              trackletME = (TrackletLight*)fTracklets->At(j);
                              if(TMath::Abs(trackletME->fEta) < TklEtaCut && (TMath::Abs(trackletME->fDPhi) < DPhiCut)){
                                  double trackletMEPhi = (Float_t)trackletME->fPhi;
                                  double trackletMEEta = trackletME->fEta;
                                  PoolsTkl[centintME].push_back(trackletMEPhi);
                                  PoolsTkl[centintME].push_back(trackletMEEta);
                                  PoolsTklEventTracker[centintME].push_back(TklMEcounter);
                                  PoolsTklEventTracker[centintME].push_back(TklMEcounter);
                              }

                          }
                          PoolsSizeTkl[centintME] += 1;
                          //  cout << "Event added"<<endl;
                      }
                        
                    }
                    
                }
            
            }
            
            
            
    }
    }
    
    cout << "Reading of events finished" <<endl;
  //  cout << "CorrelationsME, CorrelationsTklME are being set to 1" <<endl;
    
//
//    for(int i=0; i<NbBinsCent; i++){
//        for(int j=0; j<NbBinsZvtx; j++){
//            for(int k=0; k<10; k++){
//                for(int l=0; l<240; l++){
//                    CorrelationsME[i][j][k]->SetBinContent(l+1,1);
//                }
//            }
//        }
//    }
//
//    for(int i=0; i<NbBinsCent; i++){
//        for(int j=0; j<NbBinsZvtx; j++){
//                for(int l=0; l<288; l++){
//                    CorrelationsTklME[i][j]->SetBinContent(l+1,1);
//                }
//        }
//    }
    
    
    TCanvas* cCorrMEDimuTkl=new TCanvas();
    cCorrMEDimuTkl->Divide(2,2);
    cCorrMEDimuTkl->SetTitle("Correlations Dimu-Tkl [0:Central][3: m de 3.0 à 3.3 GeV]");
    cCorrMEDimuTkl->cd(1);
    Correlations[preciseCentFocus][3]->DrawCopy("colz");
    cCorrMEDimuTkl->cd(2);
    CorrelationsME[preciseCentFocus][3]->DrawCopy("colz");
    cCorrMEDimuTkl->cd(3);
    Correl_tampon = (TH1D*)(Correlations[preciseCentFocus][3]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    cCorrMEDimuTkl->cd(4);
    Correl_tampon = (TH1D*)(CorrelationsME[preciseCentFocus][3]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    
    TCanvas* cCorrMEDimuTkl2=new TCanvas();
    cCorrMEDimuTkl2->Divide(2,2);
    cCorrMEDimuTkl2->SetTitle("Correlations ME Dimu-Tkl [Ctrl/Periph][allmass]");
    cCorrMEDimuTkl2->cd(2);
    CorrelationsMEMassSummed[0]->DrawCopy("colz");
    cCorrMEDimuTkl2->cd(4);
    Correl_tampon = (TH1D*)(CorrelationsMEMassSummed[0]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    cCorrMEDimuTkl2->cd(1);
    CorrelationsMEMassSummed[4]->DrawCopy("colz");
    cCorrMEDimuTkl2->cd(3);
    Correl_tampon = (TH1D*)(CorrelationsMEMassSummed[4]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    
    TCanvas* cCorrMETKL=new TCanvas();
    cCorrMETKL->Divide(2,2);
    cCorrMETKL->SetTitle("Correlations Tkl-Tkl [0:Central]");
    cCorrMETKL->cd(1);
    CorrelationsTkl[preciseCentFocus]->DrawCopy("colz");
    cCorrMETKL->cd(2);
    CorrelationsTklME[preciseCentFocus]->DrawCopy("colz");
    cCorrMETKL->cd(3);
    Correl_tampon = (TH1D*)(CorrelationsTkl[preciseCentFocus]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    cCorrMETKL->cd(4);
    Correl_tampon = (TH1D*)(CorrelationsTklME[preciseCentFocus]->ProjectionX("_px",1,-1,"e"));
    Correl_tampon->DrawCopy("e");
    
    //TEST Michael
    YCentral->Scale(1./DimuC);
    YPeriph->Scale(1./DimuP);
    YDifference->Add(YCentral,YPeriph,1,-1);
    
    TCanvas*ctestemich = new TCanvas();
    ctestemich->SetTitle("Naïve Dimu-Tkl Yield definition TH2 Correlations / Nb Dimuons");
    ctestemich->Divide(1,3);
    ctestemich->cd(1);
    YCentral->DrawCopy("colz");
    ctestemich->cd(2);
    YPeriph->DrawCopy("colz");
    ctestemich->cd(3);
    YDifference->DrawCopy("colz");
    
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
            
            TCanvas*ctestemichTklMEProj = new TCanvas();
            ctestemichTklMEProj->SetTitle("Naïve Tkl-Tkl ME Yield definition TH2 Correlations ME / Nb ref tracklets - Projected");
            ctestemichTklMEProj->Divide(1,3);
            YTklCentralME_proj_tampon = (TH1D*)(YTklCentralME->ProjectionX("_px",1,-1,"e"));
            YTklPeriphME_proj_tampon = (TH1D*)(YTklPeriphME->ProjectionX("_px",1,-1,"e"));
            ctestemichTklMEProj->cd(1);
            
            YTklCentralME_proj_tampon->DrawCopy("e");
            ctestemichTklMEProj->cd(2);
            YTklPeriphME_proj_tampon->DrawCopy("e");
            
            TCanvas*ctestemichTklSEMEDiv = new TCanvas();
            ctestemichTklSEMEDiv->SetTitle("Naïve Tkl-Tkl SE/ME Yield definition TH2");
            ctestemichTklSEMEDiv->Divide(1,3);
            ctestemichTklSEMEDiv->cd(1);
            YTklCentral->DrawCopy("colz");
            ctestemichTklSEMEDiv->cd(2);
            YTklPeriph->DrawCopy("colz");
            
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
    }
    
    
    
    if(doMixedEvents){
    
        TCanvas*cteste0before = new TCanvas();
           cteste0before->SetTitle("Test on norm before, CorrelationsME[0][10][3]");
        cteste0before->cd();
        CorrelationsME[preciseCentFocus][3]->DrawCopy("colz");
        
        
        cout << "Normalisation of CorrelationsME based on max of Eta projected applied" <<endl;
    // Mixed events normalisation
    //        float max = 0;
    //        int maxindex = 0;
            for(int i=0; i<NbBinsCent; i++){
                
                    for(int k=0; k<10; k++){
                        float maxME = 0;
                        int maxindexME = 0;
                        hnseg8ME_proj_tampon = (TH1F*)(CorrelationsME[i][k]->ProjectionY("hnseg8ME_proj",1,-1,"e")); //Change underflow
                        for(int l=1; l<hnseg8ME_proj->GetNbinsX()+1; l++){
                            float valueME = hnseg8ME_proj_tampon->GetBinContent(l);
                           // cout << "valueME " << valueME <<endl;
                            if (valueME > maxME){
                                maxME = valueME;
                                maxindexME = l;
                            }
                        }
                        if(maxME>0){
                            NormME[i][k] = 1/maxME;
                        }
    //                    hnseg8_proj_tampon = (TH1F*)(Correlations[i][j][k]->ProjectionY("hnseg8_proj",0,-1,"e"));
    //                    for(int l=1; l<hnseg8_proj->GetNbinsX()+1; l++){
    //                        float value = hnseg8_proj_tampon->GetBinContent(l);
    //                       // cout << "value " << value <<endl;
    //                        if (value > max){
    //                            max = value;
    //                            maxindex = l;
    //                        }
    //                    }
                    }
                
            }
    
        
            //Scaling all CorrelationsME
            
        cout << "CorrelationsME[0][3] bin 6,10 : " << CorrelationsME[0][3]->GetBinContent(6,10) << endl;
        cout << "CorrelationsME[0][3] bin error 6,10 : " << CorrelationsME[0][3]->GetBinError(6,10) << endl;
        cout << "NormME[0][3] : " << NormME[0][3] << endl;
        
            for(int i=0; i<NbBinsCent; i++){
                
                       for(int k=0; k<10; k++){
                       //    CorrelationsME[i][k]->Scale(NormME[i][k]);
                           for(int binx=1; binx<(1+CorrelationsME[i][k]->GetNbinsX()); binx++){
                               for(int biny=1; biny<(1+CorrelationsME[i][k]->GetNbinsY()); biny++){
                                   if(NormME[i][k]>0){
                                   CorrelationsMEScaled[i][k]->SetBinContent(binx,biny, (CorrelationsME[i][k]->GetBinContent(binx,biny))*NormME[i][k]);
                                    CorrelationsMEScaled[i][k]->SetBinError(binx,biny, (CorrelationsME[i][k]->GetBinError(binx,biny))*NormME[i][k]);
                                   }
                               }
                           }
                       }
                   
            }
        cout << "Scaling done" <<endl;
        cout << "CorrelationsMEScaled[0][3] bin 6,10 : " << CorrelationsMEScaled[0][3]->GetBinContent(6,10) << endl;
        cout << "CorrelationsMEScaled[0][3] bin error 6,10 : " << CorrelationsMEScaled[0][3]->GetBinError(6,10) << endl;
        
        
    
//   TCanvas*cteste0after = new TCanvas();
//       cteste0after->SetTitle("Test on norm after");
//    cteste0after->cd();
//    CorrelationsME[0][0][3]->DrawCopy("colz");
    
        cout << "Normalisation of CorrelationsTklME based on max of Eta projected  applied" <<endl;
    // Mixed events normalisation Tkl
//        float maxTkl = 0;
//        int maxindexTkl = 0;
//        for(int i=0; i<NbBinsCent; i++){
//            for(int j=0; j<NbBinsZvtx; j++){
//                float maxMETkl = 0;
//                int maxindexMETkl = 0;
//                    hnseg8TklME_proj_tampon = (TH1F*)(CorrelationsTklME[i][j]->ProjectionY("hnseg8TklME_proj",0,-1,"e"));
//                    for(int l=1; l<hnseg8TklME_proj->GetNbinsX()+1; l++){
//                        float valueMETkl = hnseg8TklME_proj_tampon->GetBinContent(l);
//                       // cout << "valueME " << valueME <<endl;
//                        if (valueMETkl > maxMETkl){
//                            maxMETkl = valueMETkl;
//                            maxindexMETkl = l;
//                        }
//                    }
//                    if(maxMETkl>0){
//                        NormMETkl[i][j] = 1/maxMETkl;
//                    }
    //                hnseg8Tkl_proj_tampon = (TH1F*)(CorrelationsTkl[i][j]->ProjectionY("hnseg8Tkl_proj",0,-1,"e"));
    //                for(int l=1; l<hnseg8Tkl_proj->GetNbinsX()+1; l++){
    //                    float valueTkl = hnseg8Tkl_proj_tampon->GetBinContent(l);
    //                   // cout << "value " << value <<endl;
    //                    if (valueTkl > maxTkl){
    //                        maxTkl = valueTkl;
    //                        maxindexTkl = l;
    //                    }
    //                }
//            }
//        }
        
        
    //    if(maxMETkl > 0){
    //        NormMETkl = maxTkl/maxMETkl;
    //        cout << "NormMETkl " << NormMETkl <<endl;
    //    }
        
        //Scaling all CorrelationsMETkl
        
        for(int i=0; i<NbBinsCent; i++){
               
                     //  CorrelationsTklME[i]->Scale(1./NormMETkl[i]);
            for(int binx=1; binx<(1+CorrelationsTklME[i]->GetNbinsX()); binx++){
                for(int biny=1; biny<(1+CorrelationsTklME[i]->GetNbinsY()); biny++){
                    if(NormMETkl[i]>0){
                    CorrelationsTklMEScaled[i]->SetBinContent(binx,biny, (CorrelationsTklME[i]->GetBinContent(binx,biny))/NormMETkl[i]);
                     CorrelationsTklMEScaled[i]->SetBinError(binx,biny, (CorrelationsTklME[i]->GetBinError(binx,biny))/NormMETkl[i]);
                    }
                }
            }
               
        }
}
    
    

    cout << "Calculation of DimuonCOunterZint[Cent][mass]" <<endl;
    for(int i=0; i<NbBinsCent; i++){
              for(int k=0; k<NbinsInvMass; k++){
                  
                      DimuonCounterZint[i][k] += DimuonCounter[i][k];
                  
              }
    }
    
    {
    
    TCanvas*cteste = new TCanvas();
    cteste->SetTitle("SE Projected, ME projected and division Dimu-Tkl");
    cteste->Divide(1,3);
        
        TCanvas*ctestesum = new TCanvas();
        ctestesum->SetTitle("Sommation des yields selon z");
        ctestesum->Divide(5,2);
        
        TCanvas*ctestefin = new TCanvas();
        ctestefin->SetTitle("Yield[0][3] final sommé selon z");
        
        TCanvas*ctesteproj = new TCanvas();
        ctesteproj->SetTitle("Projction SE et ME Dimu-Tkl");
        ctesteproj->Divide(2,2);
    
        cout << "Converting count histograms into yields, by dividing by the number of dimuons seen in each case and Mixed Event and summing on Z" <<endl;
        cout << "Calculation of Yields[Cent][mass]" <<endl;
    for(int i=0; i<NbBinsCent; i++){
           for(int k=0; k<NbinsInvMass; k++){
               SoverMik->Reset();
             //  SoverMik->Sumw2();
               
               
               
                   SoverMijk->Reset();
                //   SoverMijk->Sumw2();
                   ProjCopy->Reset();
                   ProjCopy2->Reset();
                   if(doMixedEvents){
                       ME_proj_tampon = (TH1F*)(CorrelationsMEScaled[i][k]->ProjectionX("ME_proj",1,-1,"e")); //Change underflow
                   }
                   if(i ==preciseCentFocus && k==3){
                       ctesteproj->cd(1);
                       Correlations[i][k]->DrawCopy("colz");
                       ctesteproj->cd(2);
                       if(doMixedEvents){
                           CorrelationsMEScaled[i][k]->DrawCopy("colz");
                       }
                       
                   }
                   SE_proj_tampon = (TH1F*)(Correlations[i][k]->ProjectionX("SE_proj",1,-1,"e")); //Change underflow
                   ProjCopy->Add(SE_proj_tampon);
                //   ProjCopy->Sumw2();
                   if(doMixedEvents){
                       ProjCopy2->Add(ME_proj_tampon);
                   }
                //   ProjCopy2->Sumw2();
                   if(doMixedEvents){
                       SoverMijk->Divide(ProjCopy,ProjCopy2);
                   }
                   else if(!doMixedEvents){
                       SoverMijk->Add(ProjCopy);
                   }
                   if(i ==preciseCentFocus && k==3){
                       ctesteproj->cd(3);
                       ProjCopy->DrawCopy();
                       ctesteproj->cd(4);
                        if(doMixedEvents){
                            ProjCopy2->DrawCopy();
                        }
                       cteste->cd(1);
                       ProjCopy->DrawCopy();
                       cout << "Relative error on SE bin 4 is: " << ProjCopy->GetBinError(4)/ProjCopy->GetBinContent(4) << endl;
                       cteste->cd(2);
                       if(doMixedEvents){
                           ProjCopy2->DrawCopy();
                           cout << "Relative error on ME bin 4 is: " << ProjCopy2->GetBinError(4)/ProjCopy2->GetBinContent(4) << endl;
                       }
                       cteste->cd(3);
                       SoverMijk->DrawCopy();
                       cout << "Relative error on SE/ME bin 4 is: " << SoverMijk->GetBinError(4)/SoverMijk->GetBinContent(4) << endl;
                   }
                   if(i ==preciseCentFocus && k==3){
                       ctestesum->cd(1);
                       SoverMijk->DrawCopy();
                       cout << "Relative error on Mijk " << 0 << ", bin 4 is: " << SoverMijk->GetBinError(4)/SoverMijk->GetBinContent(4) << endl;
                       cout << "Error on Mijk " << 0 << ", bin 4 is: " << SoverMijk->GetBinError(4) << endl;
                   }
                   SoverMik->Add(SoverMijk);
                   if(i ==preciseCentFocus && k==3){
                       ctestesum->cd(5+1);
                       SoverMik->DrawCopy();
                       cout << "Relative error on Mik " << 0 << ", bin 4 is: " << SoverMik->GetBinError(4)/SoverMik->GetBinContent(4) << endl;
                       cout << "Error on Mik " << 0 << ", bin 4 is: " << SoverMik->GetBinError(4) << endl;
                   }
                   DimuSeenMassCut+=DimuonCounter[i][k];
               
                   
                   
                   
              // Yields[i][k]->Sumw2();
               Yields[i][k]->Add(SoverMik);
               if(DimuonCounterZint[i][k] >0){
                   Yields[i][k]->Scale(1./DimuonCounterZint[i][k]);
               }
               if(i==preciseCentFocus && k==3){
               ctestefin->cd();
               Yields[i][k]->DrawCopy("colz");
               }
    }
   }
    }
    
//    for(int i=0; i<NbBinsCent; i++){
//        for(int k=0; k<10; k++){
//            for(int p=0; p<12; p++){
//                Yields[i][k]->SetBinError(p+1, sqrt(Yields[i][k]->GetBinContent(p+1))/DimuonCounterZint[i][k]);
//            }
//        }
//    }
   
  
    cout << "DimuSeenMassCut " << DimuSeenMassCut <<endl;
    cout << "Calculation of Yields in all C, central, periph and difference" <<endl;
    
    for(int i=0; i<NbBinsCent; i++){
        for(int k=0; k<NbinsInvMass; k++){
            Yield_tampon->Reset();
            Yield_tampon->Add(Yields[i][k]);
            Yield_tampon->Scale(DimuonCounterZint[i][k]);
            Yield_allC->Add(Yield_tampon);
        }
    }
    Yield_allC->Scale(1./DimuSeenMassCut);
    
    int DimuCnt=0;
    for(int i=CentralLowBound; i<CentralHighBound+1; i++){ //CentPeriph
        Yield_tampon->Reset();
        Yield_tampon->Add(Yields[i][3]);
        Yield_tampon->Scale(DimuonCounterZint[i][3]);
        Yield_Central->Add(Yield_tampon);
        DimuCnt += DimuonCounterZint[i][3];
    }
     Yield_Central->Scale(1./DimuCnt);
    DimuCnt=0;
    for(int i=PeripheralLowBound; i<PeripheralHighBound+1; i++){
        Yield_tampon->Reset();
        Yield_tampon->Add(Yields[i][3]);
        Yield_tampon->Scale(DimuonCounterZint[i][3]);
        Yield_Periph->Add(Yield_tampon);
        DimuCnt += DimuonCounterZint[i][3];
    }
    Yield_Periph->Scale(1./DimuCnt);
    Yield_Difference->Add(Yield_Central,Yield_Periph,1,-1);
    
    
    TCanvas*cYi = new TCanvas();
    cYi->SetTitle("Yields 3.0-3.3 GeV");
    cYi->Divide(2,2);
    cYi->cd(1);
    Yield_allC->Draw("same");
    cYi->cd(2);
    Yield_Difference->Draw("same");
    cYi->cd(3);
    Yield_Central->Draw("same");
    cYi->cd(4);
    Yield_Periph->Draw("same");
    
    cout << "Calculation of Yields_MassBin[mass]" <<endl;
    for(int k=0; k<NbinsInvMass; k++){
        DimuCnt=0;
        for(int i=CentralLowBound; i<CentralHighBound+1; i++){ //CentPeriph
            Yield_tampon->Reset();
            Yield_tampon->Add(Yields[i][k]);
            Yield_tampon->Scale(DimuonCounterZint[i][k]);
            Yield_Central_MassBin[k]->Add(Yield_tampon);
            DimuCnt += DimuonCounterZint[i][k];
        }
         Yield_Central_MassBin[k]->Scale(1./DimuCnt);
        DimuCnt=0;
        for(int i=PeripheralLowBound; i<PeripheralHighBound+1; i++){
            Yield_tampon->Reset();
            Yield_tampon->Add(Yields[i][k]);
            Yield_tampon->Scale(DimuonCounterZint[i][k]);
            Yield_Periph_MassBin[k]->Add(Yield_tampon);
            DimuCnt += DimuonCounterZint[i][k];
        }
        Yield_Periph_MassBin[k]->Scale(1./DimuCnt);
        Yield_Difference_MassBin[k]->Add(Yield_Central_MassBin[k],Yield_Periph_MassBin[k],1,-1);
              
    }
        cout << "Rebinning of InvMass plots to prepare for division" <<endl;
    TH1F *InvMass_TotRebinned = dynamic_cast<TH1F*>(hnseg->Rebin(25,"InvMass_TotRebinned"));
    TH1F *InvMass_CentralRebinned = dynamic_cast<TH1F*>(InvMass_Central->Rebin(25,"InvMass_CentralRebinned"));
    TH1F *InvMass_PeriphRebinned = dynamic_cast<TH1F*>(InvMass_Periph->Rebin(25,"InvMass_PeriphRebinned"));
    
    
    cout << "Reorganisaton of Yields[Cent][mass] to create Yields_PhiBin wrt inv mass and calculation them bor all C, periph, central" <<endl;
        for(int i=0; i<NbBinsCent; i++){
            for(int k=0; k<NbinsInvMass; k++){
                Yield_tampon->Reset();
                Yield_tampon->Add(Yields[i][k]);
                Yield_tampon->Scale(DimuonCounterZint[i][k]);
                for(int deltaphibin=0; deltaphibin<NbinsDeltaPhi; deltaphibin++){
                    Yields_PhiBin[i][deltaphibin]->SetBinContent(k+1,Yield_tampon->GetBinContent(deltaphibin+1));
                    Yields_PhiBin[i][deltaphibin]->SetBinError(k+1,Yield_tampon->GetBinError(deltaphibin+1));
                }
            }
        }

    
        
        for(int p=0; p<NbinsDeltaPhi; p++){
            for(int i=0; i<NbBinsCent; i++){
                YieldWrtMass_tampon->Reset();
                YieldWrtMass_tampon->Add(Yields_PhiBin[i][p]);
                YieldWrtMass_allC[p]->Add(YieldWrtMass_tampon);
            }
            YieldWrtMass_allC[p]->Divide(InvMass_TotRebinned);
        }
        
         for(int p=0; p<NbinsDeltaPhi; p++){
            for(int i=CentralLowBound; i<CentralHighBound+1; i++){ //CentPeriph
                YieldWrtMass_tampon->Reset();
                YieldWrtMass_tampon->Add(Yields_PhiBin[i][p]);
                YieldWrtMass_Central[p]->Add(YieldWrtMass_tampon);
            }
            YieldWrtMass_Central[p]->Divide(InvMass_CentralRebinned);
        }
        
        for(int p=0; p<NbinsDeltaPhi; p++){
            for(int i=PeripheralLowBound; i<PeripheralHighBound+1; i++){
                YieldWrtMass_tampon->Reset();
                YieldWrtMass_tampon->Add(Yields_PhiBin[i][p]);
                YieldWrtMass_Periph[p]->Add(YieldWrtMass_tampon);
            }
            YieldWrtMass_Periph[p]->Divide(InvMass_PeriphRebinned);
        }
        
    
    
    
    
    if(doTracklets){
        
        cout << "Calculation of number of tracklets seen integrated on z [Cent]" <<endl;
        for(int i=0; i<NbBinsCent; i++){
            
                RefTklCounterZint[i] += RefTklCounter[i];
            
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
        CorrelationsTkl[preciseCentFocus]->DrawCopy("colz");
        ctesteprojTkl->cd(2);
        CorrelationsTklMEScaled[preciseCentFocus]->DrawCopy("colz");
        
        cout << "Calculation of YieldsTkl[Cent]" <<endl;
        for(int i=0; i<NbBinsCent; i++){
                      SoverMi->Reset();
                      
            
                          SoverMij->Reset();
                          ProjCopyTkl->Reset();
                          ProjCopy2Tkl->Reset();
                          if(doMixedEvents){
                              ME_proj_Tkl_tampon = (TH1F*)(CorrelationsTklMEScaled[i]->ProjectionX("ME_proj_Tkl",1,-1,"e")); //Change underflow
                          }
                          SE_proj_Tkl_tampon = (TH1F*)(CorrelationsTkl[i]->ProjectionX("SE_proj_Tkl",1,-1,"e")); //Change underflow
                          
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
                             if(i ==preciseCentFocus){
                                 ctesteprojTkl->cd(3);
                                 ProjCopyTkl->DrawCopy();
                                 ctesteprojTkl->cd(4);
                                 ProjCopy2Tkl->DrawCopy();
                                 ctesteTkl->cd(1);
                                 ProjCopyTkl->DrawCopy();
                                 cout << "Relative error on SE Tkl bin 4 is: " << ProjCopyTkl->GetBinError(4)/ProjCopyTkl->GetBinContent(4) << endl;
                                 ctesteTkl->cd(2);
                                 if(doMixedEvents){
                                     ProjCopy2Tkl->DrawCopy();
                                     cout << "Relative error on ME Tkl bin 4 is: " << ProjCopy2Tkl->GetBinError(4)/ProjCopy2Tkl->GetBinContent(4) << endl;
                                 }
                                 ctesteTkl->cd(3);
                                 SoverMij->DrawCopy();
                                 cout << "Relative error on SE/ME Tkl bin 4 is: " << SoverMij->GetBinError(4)/SoverMij->GetBinContent(4) << endl;
                             }
                             if(i ==preciseCentFocus){
                                 ctestesumTkl->cd(1);
                                 SoverMij->DrawCopy();
                                 cout << "Relative error on Mij Tkl " << 0 << ", bin 4 is: " << SoverMij->GetBinError(4)/SoverMij->GetBinContent(4) << endl;
                                 cout << "Error on Mij Tkl " << 0 << ", bin 4 is: " << SoverMij->GetBinError(4) << endl;
                             }
                            SoverMi->Add(SoverMij);
                             if(i ==preciseCentFocus){
                                 ctestesumTkl->cd(5+0+1);
                                 SoverMi->DrawCopy();
                                 cout << "Relative error on Mi Tkl " << 0 << ", bin 4 is: " << SoverMi->GetBinError(4)/SoverMik->GetBinContent(4) << endl;
                                 cout << "Error on Mi Tkl " << 0 << ", bin 4 is: " << SoverMi->GetBinError(4) << endl;
                             }
                        
                          RefTracklets+=RefTklCounter[i];
                      
                          
                      YieldsTkl[i]->Add(SoverMi);
                      if(RefTklCounterZint[i] >0){
                          YieldsTkl[i]->Scale(1./RefTklCounterZint[i]);
                      }
              
            if(i==preciseCentFocus){
                ctestefinTkl->cd();
                YieldsTkl[i]->DrawCopy("colz");
            }
            
          }
        
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
        for(int i=PeripheralLowBound; i<PeripheralHighBound; i++){
            YieldTkl_tampon->Reset();
            YieldTkl_tampon->Add(YieldsTkl[i]);
            YieldTkl_tampon->Scale(RefTklCounterZint[i]);
            YieldTkl_Periph->Add(YieldTkl_tampon);
            RefTklCnt += RefTklCounterZint[i];
        }
        YieldTkl_Periph->Scale(1./RefTklCnt);
        YieldTkl_Difference->Add(YieldTkl_Central,YieldTkl_Periph,1,-1);
        
    }
    
    std::cout << "==================== Analysis Terminated - PLOTTING TIME ====================" << std::endl;
    //}
    
    
    
// ***********************************
// Tracer les graphes sur les canevas *
// ***********************************


    TCanvas*c1=new TCanvas();
    //Canvas avec plein de petits plots généraux
    c1->SetTitle("General Plots");
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
    
    TCanvas*cJav=new TCanvas();
    hsingletrac->Draw("colz");
    
    TCanvas*cDPhi=new TCanvas();
    cDPhi->SetTitle("DPhi distribution");
    //Plot DPhi distribution
    hDPhi->SetFillColor(18);
    hDPhi->Draw();

    TCanvas*c2=new TCanvas();
    c2->SetTitle("TklEtaWrtZvtx");
    //Plot TrkEtaWrtZvtx
    hnseg9->Draw("colz");
    
    TCanvas*c3=new TCanvas();
    c3->SetTitle("Ntrk wrt V0M");
    //Plut Ntrk wrt V0M
    hnseg10->Draw("colz");
    
    TCanvas* c8zoom = new TCanvas;
    c8zoom->SetTitle("Yields in Mass Bins CtrlPeriphDiff zoom");
    //Correlation yield DeltaEta wrt DeltaPhi TH2 -> Central, Periph, Difference [MassBins]
    c8zoom->Divide(3,3);
    for(int i=1; i<=3; i++){
        c8zoom->cd(i);
        Yield_Central_MassBin[i+1]->Draw("E");
        c8zoom->cd(3+i);
        Yield_Periph_MassBin[i+1]->Draw("E");
        c8zoom->cd(2*3+i);
        Yield_Difference_MassBin[i+1]->Draw("E");
    }
    
    TFile *FitFile;
    FitFile = new TFile(FitFileName,"RECREATE");
    
    TCanvas* c8 = new TCanvas;
    // ROOT THAT Yield in bins of mass
    c8->SetTitle("Yields in Mass Bins CtrlPeriphDiff");
    //Correlation yield DeltaEta wrt DeltaPhi TH2 -> Central, Periph, Difference [MassBins]
    c8->Divide(NbinsInvMass,3);
    for(int i=1; i<=NbinsInvMass; i++){
        c8->cd(i);
        Yield_Central_MassBin[i-1]->Draw("E");
        Yield_Central_MassBin[i-1]->Write();
        c8->cd(NbinsInvMass+i);
        Yield_Periph_MassBin[i-1]->Draw("E");
        Yield_Periph_MassBin[i-1]->Write();
        c8->cd(2*NbinsInvMass+i);
        Yield_Difference_MassBin[i-1]->Draw("E");
        Yield_Difference_MassBin[i-1]->Write();
    }
    
    for(int i=1; i<NbinsInvMass+1; i++){
        baselines0->Fill(MinInvMass + (-0.5+i)*SizeBinInvMass, (Yield_Periph_MassBin[i-1]->GetBinContent(3)+Yield_Periph_MassBin[i-1]->GetBinContent(4))/2);
        baselines0->SetBinError(i, sqrt(pow(Yield_Periph_MassBin[i-1]->GetBinError(3),2)+pow(Yield_Periph_MassBin[i-1]->GetBinError(4),2))/2);
    }
    
    TCanvas* c10 = new TCanvas;
    // ROOT THAT yields phi bins
    c10->SetTitle("Yields wrt mass Central Periph");
    // Yields wrt mass -> Periph, Central [Phi bins]
    c10->Divide(6,4);
    for(int i=1; i<NbinsDeltaPhi+1; i++){
        c10->cd(i);
        YieldWrtMass_Central[i-1]->Draw("E");
        YieldWrtMass_Central[i-1]->Write();
        c10->cd(NbinsDeltaPhi+i);
        YieldWrtMass_Periph[i-1]->Draw("E");
        YieldWrtMass_Periph[i-1]->Write();
    }
    
    FitFile->Close();
    
    TCanvas* c11a = new TCanvas;
    c11a->SetTitle("Pt wrt Mass CentralPeriphAll");
    //Pt wrt Mass -> Central, Periph, All
    c11a->Divide(3,1);
    for(int i=0; i<3; i++){
        c11a->cd(i+1);
        hPtWrtMassInv[i]->Draw("colz");
    }
    TCanvas* c11b_0 = new TCanvas;
    c11b_0->SetTitle("Mass all C, pT bins");
    //Mass -> All C [Pt bins]
    c11b_0->Divide(2,3);
    for (int i=0;i<6;i++){
        hPtWrtMassInvSliced[0][i] = hPtWrtMassInv[0]->ProjectionX(Form("bin%d_0",i+1),i+1,i+1);
    }
    TCanvas* c11b_1 = new TCanvas;
    c11b_1->SetTitle("Mass Central, pT bins");
    //Mass -> Central [Pt bins]
    c11b_1->Divide(2,3);
    for (int i=0;i<6;i++){
        hPtWrtMassInvSliced[1][i] = hPtWrtMassInv[1]->ProjectionX(Form("bin%d_1",i+1),i+1,i+1);
    }
    TCanvas* c11b_2 = new TCanvas;
    c11b_2->SetTitle("Mass Periph, pT bins");
    //Mass -> Periph [Pt bins]
    c11b_2->Divide(2,3);
    for (int i=0;i<6;i++){
        hPtWrtMassInvSliced[2][i] = hPtWrtMassInv[2]->ProjectionX(Form("bin%d_2",i+1),i+1,i+1);
    }
    
    TCanvas*c12=new TCanvas();
    c12->SetTitle("Centrality Estimators");
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
    c13->SetTitle("Estimators wrt tkl");
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
    c14->SetTitle("Yield wrt Phi difference fits");
    //Yield difference wrt Phi [Mass bins] fits
    c14->Divide(5,2);
    for(int i=1; i<=10; i++){
        c14->cd(i);
    // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
        TH1F *histo = Yield_Difference_MassBin[i-1];
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
        cout << "STATUS COV : " << res->CovMatrixStatus() <<endl;
      // draw the legend
      TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
      legend->SetTextFont(72);
      legend->SetTextSize(0.04);
        Char_t message[80];
        sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
        legend->AddEntry(fitFcnV2_2,message);
        if(res->CovMatrixStatus() == 3){
            sprintf(message,"The fit is a success");
        }
        else{
            sprintf(message,"The fit is a failure");
        }
        legend->AddEntry(fitFcnV2_2,message);
      legend->AddEntry(histo,"Data","lpe");
      legend->Draw();
    
    }
    
    TFile *f = new TFile(FitFileName,"UPDATE");
    
    TCanvas*c15=new TCanvas();
    // ROOT THAT Fouriers
    c15->SetTitle("Fourier coefs of yields wrt Phi");
    // Fourier coefficients of Yields wrt Phi [Mass bins] -> Extraction method 2
    c15->Divide(2,2);
    c15->cd(1);
    coefficients0->Draw();
    coefficients0->Write();
    c15->cd(2);
    baselines0->Draw();
    baselines0->Write();
    c15->cd(3);
    coefficients1->Draw();
    coefficients1->Write();
    c15->cd(4);
    coefficients2->Draw();
    coefficients2->Write();
            

    // Plot du V2 JPsi Tracklet
    TCanvas*c16=new TCanvas();
    // ROOT THAT Extraction 2 plot
    c16->SetTitle("Extraction method 2");
    //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
    c2b0->Add(coefficients0, baselines0);
    V2JPsiTkl->Divide(coefficients2, c2b0);
    
    V2JPsiTkl->Draw();
    V2JPsiTkl->Write();
    
    hnseg->Write();
    InvMass_Central->Write();
    InvMass_Periph->Write();

    f->Close();
    
    
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
    cinvmass->SetTitle("Inv Mass Fits");
    // ROOT THAT Inv mass fits
    cinvmass->Divide(1,3);
    
    TFitResultPtr res;
    TFitResultPtr rescent;
    TFitResultPtr resperiph;
    TFitResult resu;
    char histoname[50];
    //Values OK for 20M evts
//    double startvalues[16] = {800, 3.096916, 0.07, 0.9, 10, 2, 15, 15,   10,     0.1,          50000,   0.01,        1,         1,   1,   1};
//    double lowvalues[16] =   {1,    3.05,    0.03, 0.7, 1,  1, 5,  0.01, 0.0001,  0.000000001,  10,    0.0000000001, 0.000001, -10, -10, -10};
//    double highvalues[16] = {100000, 3.15,   0.1,  1.1, 30, 5, 25, 50,   100000, 10,           100000, 10,     10,        10,  10,  10};
    
    double startvalues[16] = {1100, 3.096916, 0.07,  0.9, 10, 2, 15, 20,    90,      0.01,     80000,      1.4,      1,         1,   1,   1};
    double lowvalues[16] =   {10,    3.05,     0.03, 0.7, 1,  1, 5,  0.001,  1,    0.00001,   1000,        0.0001, 0.00000001, -10, -10, -10};
    double highvalues[16] = {100000, 3.15,     1,  1.1, 30, 5, 25, 10000, 1000000, 100,      10000000,  100,    10,        10,  10,  10};
    
//    for(int index=0; index<16; index++){
//        if(index==0 || index==7 || index==8 || index==10){
//            startvalues[index]*=nevent1/20000000;
//            lowvalues[index]*=nevent1/20000000;
//            highvalues[index]*=nevent1/20000000;
//        }
//    }
    //Values adapted to the number of events considered
    
    sprintf(histoname,"hnseg");
    if(!CombineFits){
        res = FittingAllInvMassBin(histoname, cinvmass, 0);
    }
    double par[16];
    
    // Parameter setting (either from mass fit either start values)
    for(int i=0; i <16; i++){
        if(!CombineFits){
            par[i] = res->Parameter(i);
            if(i>11){
                par[i] = 1;
            }
        }
        if(CombineFits){
            par[i] = startvalues[i];
        }
    }
    
    //Definition and design of function on V2 plot
    c16->cd();
    TF1 *fitV2_2 = new TF1("fitV2_2",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
    fitV2_2->SetNpx(500);
    fitV2_2->SetLineWidth(4);
    fitV2_2->SetLineColor(kMagenta);
    fitV2_2->SetParameters(par);
    TF1 *fitJustV2_2 = new TF1("fitJustV2_2",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,4);
    fitJustV2_2->SetNpx(500);
    fitJustV2_2->SetLineWidth(0);
    fitJustV2_2->SetLineColor(kWhite);
    TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2Poly,2.1,5.1,16);
    backFcnV2_2->SetLineColor(kRed);
    TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,2.1,5.1,16);
    signalFcnJPsiV2_2->SetLineColor(kBlue);
    signalFcnJPsiV2_2->SetNpx(500);
    
    if(!CombineFits){
       TVirtualFitter::Fitter(V2JPsiTkl)->SetMaxIterations(10000);
       TVirtualFitter::Fitter(V2JPsiTkl)->SetPrecision();
        for(int i=0; i<=11; i++){
            fitV2_2->FixParameter(i,par[i]);
        }
    }
       
       fitV2_2->SetParName(0,"Norm_{JPsi}");
       fitV2_2->SetParName(1,"M_{JPsi}");
       fitV2_2->SetParName(2,"Sigma_{JPsi}");
       fitV2_2->SetParName(3,"a_{1}");
       fitV2_2->SetParName(4,"n_{1}");
       fitV2_2->SetParName(5,"a_{2}");
       fitV2_2->SetParName(6,"n_{2}");
       fitV2_2->SetParName(7,"Norm_{Psi2S}");
       fitV2_2->SetParName(8,"Norm_{TailLowM}");
       fitV2_2->SetParName(9,"Exp_{TailLowM}");
       fitV2_2->SetParName(10,"Norm_{TailHighM}");
       fitV2_2->SetParName(11,"Exp_{TailHighM}");
        fitV2_2->SetParName(12,"V2_2 JPsi");
        fitV2_2->SetParName(13,"V2_2 Bkg M2");
    fitV2_2->SetParName(14,"V2_2 Bkg M1");
    fitV2_2->SetParName(15,"V2_2 Nkg M0");
    
    fitJustV2_2->SetParName(0,"V2_2 JPsi");
        fitJustV2_2->SetParName(1,"V2_2 Bkg M2");
    fitJustV2_2->SetParName(2,"V2_2 Bkg M1");
    fitJustV2_2->SetParName(3,"V2_2 Nkg M0");
    
       double minParams[16];
       double parErrors[16];
    
    // Combined fit of inv mass and V2_2
            if(CombineFits){
                TBackCompFitter * virminuit = (TBackCompFitter *) TVirtualFitter::Fitter(0,16);
                   for (int i = 0; i < 16; ++i) {
                     virminuit->SetParameter(i, fitV2_2->GetParName(i), fitV2_2->GetParameter(i), 0.01, lowvalues[i],highvalues[i]);
                   }
                virminuit->FixParameter(1);
                virminuit->FixParameter(2);
                virminuit->FixParameter(3);
                virminuit->FixParameter(4);
                virminuit->FixParameter(5);
                virminuit->FixParameter(6);
                   virminuit->SetFCN(FcnCombinedAllMass);

                   double arglist[100];
                   arglist[0] = 1;
                   // set print level
                   virminuit->ExecuteCommand("SET PRINT",arglist,2);

                   // minimize
                   arglist[0] = 5000; // number of function calls
                   arglist[1] = 0.01; // tolerance
                   virminuit->ExecuteCommand("MIGRAD",arglist,2);
                
                gStyle->SetOptFit(1011);
                
                   virminuit->ReleaseParameter(1);
                virminuit->ExecuteCommand("MIGRAD",arglist,2);
                 virminuit->ReleaseParameter(2);
                virminuit->ExecuteCommand("MIGRAD",arglist,2);
                 virminuit->ReleaseParameter(7);
                virminuit->ExecuteCommand("MIGRAD",arglist,2);
                virminuit->ExecuteCommand("MIGRAD",arglist,2);
                resu = virminuit->GetFitResult();
              //  resu = minuit->GetFitResult();
                std::cout << "pront";
                resu.GetCovarianceMatrix().Print();

                   //get result
                   for (int i = 0; i < 16; ++i) {
                     minParams[i] = virminuit->GetParameter(i);
                     parErrors[i] = virminuit->GetParError(i);
                   }
                   double chi2, edm, errdef;
                   int nvpar, nparx;
                   virminuit->GetStats(chi2,edm,errdef,nvpar,nparx);

                   fitV2_2->SetParameters(minParams);
                   fitV2_2->SetParErrors(parErrors);
                fitJustV2_2->SetParameters(&minParams[12]);
                fitJustV2_2->SetParErrors(&parErrors[12]);
                   fitV2_2->SetChisquare(chi2);
                   int ndf = npfits-nvpar;
                   fitV2_2->SetNDF(ndf);
            }
    //Fit of V2
    if(!CombineFits){
        res = V2JPsiTkl->Fit("fitV2_2","SBMERI+","ep");
        Double_t para[16];
        fitV2_2->GetParameters(para);
        backFcnV2_2->SetParameters(para);
        signalFcnJPsiV2_2->SetParameters(para);
    }
    
    if(CombineFits){
        backFcnV2_2->SetParameters(minParams);
        signalFcnJPsiV2_2->SetParameters(minParams);
    }
    
    fitV2_2->Draw("same");
    if(CombineFits){
        V2JPsiTkl->GetListOfFunctions()->Add(fitJustV2_2);
    }
  //  signalFcnJPsiV2_2->Draw("same");
    backFcnV2_2->Draw("same");
      // draw the legend
      TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
      legend->SetTextFont(72);
      legend->SetTextSize(0.04);
    Char_t message[80];
    sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
     legend->AddEntry(fitV2_2,message);
    if(CombineFits){
        if(resu.CovMatrixStatus() == 3){
               sprintf(message,"The fit is a success");
        }
        else{
               sprintf(message,"The fit is a failure");
        }
    }
    if(!CombineFits){
        if(res->CovMatrixStatus() == 3){
               sprintf(message,"The fit is a success");
        }
        else{
               sprintf(message,"The fit is a failure");
        }
    }
           legend->AddEntry(fitV2_2,message);
  //  legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
    legend->AddEntry(backFcnV2_2,"Background");
      legend->AddEntry(V2JPsiTkl,"Data","lpe");
      legend->Draw();
    
    //Plot of inv mass
    if(CombineFits){
        TFile *file0 = new TFile("/Volumes/SEBUSB/InvMass.root");
            
            cinvmass->cd(1);
           cinvmass->SetFillColor(33);
           cinvmass->SetFrameFillColor(41);
           cinvmass->SetGrid();
           TH1F *histo = (TH1F*)file0->Get(histoname);
           // create a TF1 with the range from 0 to 3 and 6 parameters
           TF1 *fitFcn = new TF1("fitFcn",TwoCBE2E,2.1,5.1,12);
           fitFcn->SetNpx(500);
           fitFcn->SetLineWidth(4);
           fitFcn->SetLineColor(kMagenta);

            fitFcn->SetParName(0,"Norm_{JPsi}");
            fitFcn->SetParName(1,"M_{JPsi}");
            fitFcn->SetParName(2,"Sigma_{JPsi}");
            fitFcn->SetParName(3,"a_{1}");
            fitFcn->SetParName(4,"n_{1}");
            fitFcn->SetParName(5,"a_{2}");
            fitFcn->SetParName(6,"n_{2}");
            fitFcn->SetParName(7,"Norm_{Psi2S}");
            fitFcn->SetParName(8,"Norm_{TailLowM}");
            fitFcn->SetParName(9,"Exp_{TailLowM}");
            fitFcn->SetParName(10,"Norm_{TailHighM}");
            fitFcn->SetParName(11,"Exp_{TailHighM}");

           // improve the picture:
           TF1 *backFcn = new TF1("backFcn",ExpBkg,2.1,5.1,4);
           backFcn->SetLineColor(kRed);
           TF1 *signalFcnJPsi = new TF1("signalFcnJPsi",JPsiCrystalBallExtended,2.1,5.1,8);
           TF1 *signalFcnPsi2S = new TF1("signalFcnPsi2S",Psi2SCrystalBallExtended,2.1,5.1,8);
           TPaveText *pave = new TPaveText(0.15,0.5,0.3,0.65,"brNDC");
           signalFcnJPsi->SetLineColor(kBlue);
           signalFcnJPsi->SetNpx(500);
            signalFcnPsi2S->SetLineColor(kGreen);
            signalFcnPsi2S->SetNpx(500);
            histo->GetListOfFunctions()->Add(fitFcn);
           // writes the fit results into the par array
            gStyle->SetOptFit(1011);
            fitFcn->SetParameters(minParams);
            fitFcn->SetParErrors(parErrors);
           signalFcnJPsi->SetParameters(minParams);
            signalFcnPsi2S->SetParameters(minParams);
         //  auto covMatrix = res->GetCovarianceMatrix();
           std::cout << "Covariance matrix from the fit ";
           resu.GetCovarianceMatrix().Print();
           Double_t integral = (signalFcnJPsi->Integral(2.1,3.45))/0.01;
            std::cout << "STATUS COV " << resu.CovMatrixStatus() <<endl;
           Double_t integralerror = (signalFcnJPsi->IntegralError(2.1,3.45,signalFcnJPsi->GetParameters(), resu.GetCovarianceMatrix().GetSub(0,8,0,8).GetMatrixArray() ))/0.01;
            std::cout << "Erreur integrale " << integralerror <<endl;
            
            std::cout << "Fitted " << histoname << std::endl;
            std::cout << "Nb JPsi total: " << integral << std::endl;
         //   std::cout << "integral error: " << integralerror << std::endl;
        histo->Draw();
        fitFcn->Draw("same");
           signalFcnJPsi->Draw("same");
            signalFcnPsi2S->Draw("same");
           backFcn->SetParameters(&minParams[8]);
           backFcn->Draw("same");
           // draw the legend
            char str[50];
            sprintf(str, "N_{JPsi} %f +/- %f", integral, integralerror);
           pave->AddText(str);
            sprintf(str, "M_{Psi2S} = %f, Sig_{Psi2S} = %f", mPsip, ratSigma*minParams[2]);
            pave->AddText(str);
           pave->Draw();
           TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
           legend->SetTextFont(72);
           legend->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
             legend->AddEntry(fitV2_2,message);
            if(resu.CovMatrixStatus() == 3){
                       sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                   }
                   legend->AddEntry(fitV2_2,message);
           legend->AddEntry(histo,"Data","lpe");
           legend->AddEntry(backFcn,"Background fit","l");
           legend->AddEntry(signalFcnJPsi,"JPsi Signal fit","l");
            legend->AddEntry(signalFcnPsi2S,"Psi2S Signal fit","l");
           legend->Draw();
    }
    
    
    //Combining or not fit for Central-Periph method Ext1
    
    
    int niterations = 1;
    if(CombineFits){
        niterations = 1;
    }
    cinvmass->cd();
    sprintf(histoname,"InvMass_Central");
   res = FittingAllInvMassBin(histoname, cinvmass, 1); //ZZZZZZZ
    
    TCanvas* chnsegsigma = new TCanvas;
    chnsegsigma->cd();
    hnsegSigma->Draw("colz");

    double parerr[12];
    TCanvas* c10_fit = new TCanvas;
       c10_fit->Divide(6,4);
    TRandom rand = 0;
    for(int j=0; j<12; j++){
        parerr[j] = res->ParError(j);
    }
    for(int j=0; j <16; j++){
        par[j] = res->Parameter(j);
        if(j>11){
            par[j] = 3;
        }
    }

    double par12[niterations];
    double parerr12[niterations];
    
    for(int i=0; i<NbinsDeltaPhi; i++){
        for(int t=0; t<niterations; t++){
        TVirtualPad* subpad = c10_fit->cd(i+1);
            c10_fit->cd(i+1);
            if(t==niterations-1){
                subpad->Clear();
            }
        TF1 *fitY_1Central = new TF1("fitY_1Central",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
          fitY_1Central->SetNpx(500);
          fitY_1Central->SetLineWidth(4);
          fitY_1Central->SetLineColor(kMagenta);
          fitY_1Central->SetParameters(par);
           TVirtualFitter::Fitter(YieldWrtMass_Central[i])->SetMaxIterations(10000);
           TVirtualFitter::Fitter(YieldWrtMass_Central[i])->SetPrecision();
        //  histo->Fit("fitFcn","0");
          // second try: set start values for some parameters
        for(int k=0; k<=11; k++){
            if(CombineFits){
                double valuepar = par[k] + rand.Gaus(0,parerr[k]);
                fitY_1Central->FixParameter(k,valuepar);
                std::cout<<"valuepar: " << valuepar << " i " << k <<endl;
            }
            if(!CombineFits){
                fitY_1Central->FixParameter(k,par[k]);
            }
        }
            fitY_1Central->SetParLimits(12,0,10000);

           fitY_1Central->SetParName(0,"Norm_{JPsi}");
           fitY_1Central->SetParName(1,"M_{JPsi}");
           fitY_1Central->SetParName(2,"Sigma_{JPsi}");
           fitY_1Central->SetParName(3,"a_{1}");
           fitY_1Central->SetParName(4,"n_{1}");
           fitY_1Central->SetParName(5,"a_{2}");
           fitY_1Central->SetParName(6,"n_{2}");
           fitY_1Central->SetParName(7,"Norm_{Psi2S}");
           fitY_1Central->SetParName(8,"Norm_{TailLowM}");
           fitY_1Central->SetParName(9,"Exp_{TailLowM}");
           fitY_1Central->SetParName(10,"Norm_{TailHighM}");
           fitY_1Central->SetParName(11,"Exp_{TailHighM}");
            fitY_1Central->SetParName(12,"Y_1 JPsi");
            fitY_1Central->SetParName(13,"Y_1 Bkg M2");
        fitY_1Central->SetParName(14,"Y_1 Bkg M1");
        fitY_1Central->SetParName(15,"Y_1 Nkg M0");

          rescent = YieldWrtMass_Central[i]->Fit("fitY_1Central","SBMERIQ+","ep");
          Double_t param[16];
            Double_t paramerrs[16];
          fitY_1Central->GetParameters(param);
            for(int i=0; i<16; i++){
                paramerrs[i] = fitY_1Central->GetParError(i);
            }
            
            if(rescent->CovMatrixStatus() == 3){
            par12[t] = param[12];
            parerr12[t] = paramerrs[12];
            }
            else{
                par12[t]=parerr12[t]=0;
            }
            
            std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
            std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

        TF1 *backFcnY_1Central = new TF1("backFcnY_1Central",BackFcnV2Poly,2.1,5.1,16);
        backFcnY_1Central->SetLineColor(kRed);
        TF1 *signalFcnJPsiY_1Central = new TF1("signalFcnJPsiY_1Central",SignalFcnJPsiV2,2.1,5.1,16);
          // writes the fit results into the par array
        
            if(t==niterations-1){
                YieldWrtMass_Central[i]->GetListOfFunctions()->Clear();
                YieldWrtMass_Central[i]->GetListOfFunctions()->Add(fitY_1Central);
            }
        signalFcnJPsiY_1Central->SetLineColor(kBlue);
        signalFcnJPsiY_1Central->SetNpx(500);
        backFcnY_1Central->SetParameters(param);
        signalFcnJPsiY_1Central->SetParameters(param);
      //  signalFcnJPsiY_1Central->Draw("same");
        backFcnY_1Central->Draw("same");
            
            if(!CombineFits){
                Yields_Central_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[12]);
                Yields_Central_1->SetBinError(i+1,fitY_1Central->GetParError(12));
            }
            
          // writes the fit results into the par array
           gStyle->SetOptFit(1011);
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetTextFont(72);
          legend->SetTextSize(0.04);
        Char_t message[80];
        sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
         legend->AddEntry(fitY_1Central,message);
        if(rescent->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitY_1Central,message);
      //  legend->AddEntry(signalFcnJPsiY_1Central,"JPsi signal");
        legend->AddEntry(backFcnY_1Central,"Background");
          legend->AddEntry(YieldWrtMass_Central[i],"Data","lpe");
          legend->Draw();
        }
        
        if(CombineFits){
            double moypar12 = 0;
            double moyparerr12 = 0;
            double typeApar12 = 0;
            int failures = 0;
            for(int z=0; z<niterations; z++){
                moypar12 += par12[z];
                moyparerr12 += parerr12[z];
                if(par12[z]==0 && parerr12[z]==0){
                    failures++;
                }
            }
            moypar12/=(niterations-failures);   // Moyenne de par12
            moyparerr12/=(niterations-failures); //Incertitude fit 2
            for(int z=0; z<niterations; z++){
                if(par12[z]!=0 || parerr12[z]!=0){
                    typeApar12 += pow(moypar12-par12[z],2);
                }
            }
            typeApar12/=((niterations-failures)*(niterations-1-failures));
            typeApar12 = sqrt(typeApar12); // Incertitude de type A (répétabitilé) issue du fit 1 propagé
            
            std::cout << "There were " << failures << " failures"<<endl;
            std::cout << "moypar12: " << moypar12 << ", moyparerr12: " << moyparerr12 << ", typeApar12: "<<typeApar12<<endl;
            Yields_Central_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, moypar12);
            Yields_Central_1->SetBinError(i+1,sqrt(pow(moyparerr12,2)+pow(typeApar12,2)));
        }
    }
    
    cinvmass->cd();
    sprintf(histoname,"InvMass_Periph");
    res = FittingAllInvMassBin(histoname, cinvmass, 2); //ZZZZZ
    
       for(int j=0; j<12; j++){
           parerr[j] = res->ParError(j);
       }
       for(int j=0; j <16; j++){
           par[j] = res->Parameter(j);
           if(j>11){
               par[j] = 3;
           }
       }
    
    double baseline = 1;
    double errbaseline = 1;
       
       for(int i=0; i<NbinsDeltaPhi; i++){
           for(int t=0; t<niterations; t++){
               TVirtualPad* subpad = c10_fit->cd(NbinsDeltaPhi+i+1);
               c10_fit->cd(NbinsDeltaPhi+i+1);
               if(t==niterations-1){
                   subpad->Clear();
               }
           TF1 *fitY_1Periph = new TF1("fitY_1Periph",FourierV2_WrtInvMass,MinInvMass,MaxInvMass,16);
             fitY_1Periph->SetNpx(500);
             fitY_1Periph->SetLineWidth(4);
             fitY_1Periph->SetLineColor(kMagenta);
             fitY_1Periph->SetParameters(par);
              TVirtualFitter::Fitter(YieldWrtMass_Periph[i])->SetMaxIterations(10000);
              TVirtualFitter::Fitter(YieldWrtMass_Periph[i])->SetPrecision();
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters
            for(int k=0; k<=11; k++){
                if(CombineFits){
                      double valuepar = par[k] + rand.Gaus(0,parerr[k]);
                      fitY_1Periph->FixParameter(k,valuepar);
//                      std::cout<<"valuepar: " << valuepar << " k " << k <<endl;
                }
                if(!CombineFits){
                    fitY_1Periph->FixParameter(k,par[k]);
                }
            }
              fitY_1Periph->SetParLimits(12,0,10000);
              
              fitY_1Periph->SetParName(0,"Norm_{JPsi}");
              fitY_1Periph->SetParName(1,"M_{JPsi}");
              fitY_1Periph->SetParName(2,"Sigma_{JPsi}");
              fitY_1Periph->SetParName(3,"a_{1}");
              fitY_1Periph->SetParName(4,"n_{1}");
              fitY_1Periph->SetParName(5,"a_{2}");
              fitY_1Periph->SetParName(6,"n_{2}");
              fitY_1Periph->SetParName(7,"Norm_{Psi2S}");
              fitY_1Periph->SetParName(8,"Norm_{TailLowM}");
              fitY_1Periph->SetParName(9,"Exp_{TailLowM}");
              fitY_1Periph->SetParName(10,"Norm_{TailHighM}");
              fitY_1Periph->SetParName(11,"Exp_{TailHighM}");
               fitY_1Periph->SetParName(12,"Y_1 JPsi");
               fitY_1Periph->SetParName(13,"Y_1 Bkg M2");
           fitY_1Periph->SetParName(14,"Y_1 Bkg M1");
           fitY_1Periph->SetParName(15,"Y_1 Nkg M0");
              
             resperiph = YieldWrtMass_Periph[i]->Fit("fitY_1Periph","SBMERIQ+","ep");
               Double_t param[16];
               Double_t paramerrs[16];
             fitY_1Periph->GetParameters(param);
               
               for(int i=0; i<16; i++){
                               paramerrs[i] = fitY_1Periph->GetParError(i);
                           }
                           
                           if(resperiph->CovMatrixStatus() == 3){
                           par12[t] = param[12];
                           parerr12[t] = paramerrs[12];
                           }
                           else{
                               par12[t]=parerr12[t]=0;
                           }
               std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
                std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

           TF1 *backFcnY_1Periph = new TF1("backFcnY_1Periph",BackFcnV2Poly,2.1,5.1,16);
           backFcnY_1Periph->SetLineColor(kRed);
           TF1 *signalFcnJPsiY_1Periph = new TF1("signalFcnJPsiY_1Periph",SignalFcnJPsiV2,2.1,5.1,16);
             // writes the fit results into the par array
               
               if(t==niterations-1){
                   YieldWrtMass_Periph[i]->GetListOfFunctions()->Clear();
                   YieldWrtMass_Periph[i]->GetListOfFunctions()->Add(fitY_1Periph);
               }
           
           signalFcnJPsiY_1Periph->SetLineColor(kBlue);
           signalFcnJPsiY_1Periph->SetNpx(500);
           backFcnY_1Periph->SetParameters(param);
           signalFcnJPsiY_1Periph->SetParameters(param);
         //  signalFcnJPsiY_1Periph->Draw("same");
           backFcnY_1Periph->Draw("same");
               
               if(!CombineFits){
                   if(i==2){
                       baseline = param[12];
                       errbaseline = fitY_1Periph->GetParError(12);
                   }
                   if(i==3){
                       baseline += param[12];
                       baseline /= 2;
                       errbaseline = sqrt(pow(errbaseline,2)+pow(fitY_1Periph->GetParError(12),2));
                   }
                 Yields_Periph_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[12]);
                 Yields_Periph_1->SetBinError(i+1,fitY_1Periph->GetParError(12));
               }
               
             // writes the fit results into the par array
              gStyle->SetOptFit(1011);
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetTextFont(72);
             legend->SetTextSize(0.04);
           Char_t message[80];
           sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
            legend->AddEntry(fitY_1Periph,message);
           if(resperiph->CovMatrixStatus() == 3){
                      sprintf(message,"The fit is a success");
                  }
                  else{
                      sprintf(message,"The fit is a failure");
                  }
                  legend->AddEntry(fitY_1Periph,message);
        //   legend->AddEntry(signalFcnJPsiY_1Periph,"JPsi signal");
           legend->AddEntry(backFcnY_1Periph,"Background");
             legend->AddEntry(YieldWrtMass_Periph[i],"Data","lpe");
             legend->Draw();
           }
           if(CombineFits){
                  double moypar12 = 0;
                  double moyparerr12 = 0;
                  double typeApar12 = 0;
                  int failures = 0;
                  for(int z=0; z<niterations; z++){
                      moypar12 += par12[z];
                      moyparerr12 += parerr12[z];
                      if(par12[z]==0 && parerr12[z]==0){
                          failures++;
                      }
                  }
                  moypar12/=(niterations-failures);   // Moyenne de par12
                  moyparerr12/=(niterations-failures); //Incertitude fit 2
                  for(int z=0; z<niterations; z++){
                      if(par12[z]!=0 || parerr12[z]!=0){
                          typeApar12 += pow(moypar12-par12[z],2);
                      }
                  }
                  typeApar12/=((niterations-failures)*(niterations-1-failures));
                  typeApar12 = sqrt(typeApar12); // Incertitude de type A (répétabitilé) issue du fit 1 propagé
                  
                  std::cout << "There were " << failures << " failures"<<endl;
                  std::cout << "moypar12: " << moypar12 << ", moyparerr12: " << moyparerr12 << ", typeApar12: "<<typeApar12<<endl;
                  Yields_Periph_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, moypar12);
                  Yields_Periph_1->SetBinError(i+1,sqrt(pow(moyparerr12,2)+pow(typeApar12,2)));
           
               if(i==2){
                   baseline = moypar12;
                   errbaseline = sqrt(pow(moyparerr12,2)+pow(typeApar12,2));
               }
               if(i==3){
                   baseline += moypar12;
                   baseline /= 2;
                   errbaseline = sqrt(pow(errbaseline,2)+pow(sqrt(pow(moyparerr12,2)+pow(typeApar12,2)),2));
               }
           }
           
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
        if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitFcnV2_2,message);
             legend->AddEntry(Yields_Difference_1,"Data","lpe");
             legend->Draw();
           
    }
    
    
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
        
        baselineTKL = (YieldTkl_Periph->GetBinContent(3) + YieldTkl_Periph->GetBinContent(4))/2;
        errbaselineTKL = sqrt(pow(YieldTkl_Periph->GetBinError(3),2) + pow(YieldTkl_Periph->GetBinError(4),2));
            
        TCanvas*c14TKL=new TCanvas();
        //Tracklets Yield difference wrt Phi fit
        c14TKL->cd();
        // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
            TH1F *histo = YieldTkl_Difference;
          // create a TF1 with the range from 0 to 3 and 6 parameters
          TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV5,-TMath::Pi()/2,1.5*TMath::Pi(),3);
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
          //  fitFcnV2TKL->SetParName(3,"a3");
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
          // improve the pictu
        //   std::cout << "integral error: " << integralerror << std::endl;
        histo->SetStats(kFALSE);
            Double_t par[3];
            fitFcnV2TKL->GetParameters(par);
          fitFcnV2TKL->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetTextFont(61);
          legend->SetTextSize(0.02);
            Char_t message[80];
            sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
            legend->AddEntry(fitFcnV2TKL,message);
        sprintf(message,"V2 Tkl-Tkl: %.4f / (%.4f + %.4f) = %.4f +- %.4f",par[2],par[0], baselineTKL,par[2]/(par[0] + baselineTKL),(par[2]/(par[0] + baselineTKL)*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)+pow(errbaselineTKL/baselineTKL,2))));
        legend->AddEntry(fitFcnV2TKL,message);
        if(res->CovMatrixStatus() == 3){
                   sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
               }
               legend->AddEntry(fitFcnV2TKL,message);
          legend->AddEntry(histo,"Data","lpe");
          legend->Draw();
        
    }
    
    
    
//        for (int i=0;i<6;i++){
//               sprintf(histoname,Form("bin%d_0",i+1)); // ZZZZZZZ
//               res = FittingAllInvMassBin(histoname, c11b_0, i);
//           }
//    for (int i=0;i<6;i++){
//        sprintf(histoname,Form("bin%d_1",i+1));
//        res = FittingAllInvMassBin(histoname, c11b_1, i);
//    }
//    for (int i=0;i<6;i++){
//        sprintf(histoname,Form("bin%d_2",i+1));
//        res = FittingAllInvMassBin(histoname, c11b_2, i);
//    }
    
    
    
    
    
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
    std::cout << "EventRejected: " << EventRejected <<std::endl;
    std::cout << "EventPileUpVtx: " << EventPileUpVtx <<std::endl;
    std::cout << "EventPileUpMult: " << EventPileUpMult <<std::endl;
    std::cout << "EventNC: " << EventNC <<std::endl;
    std::cout << "cmul  " << cmul <<std::endl;
    std::cout <<"countsigma :" << countsigma <<std::endl;
}
    









// FITTING INVARIANT MASS METHODS

Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par)
// Par 0->7: signal, 8->11 Bkg, 12: v2 JPsi, 13->15: V2 bkg
{ return (JPsiCrystalBallExtended(x,par)*par[12] + (ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*x[0]*x[0] + par[14]*x[0] + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t BackFcnV2(Double_t *x,Double_t *par)
{return ((ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*x[0]*x[0] + par[14]*x[0] + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t BackFcnV2Poly(Double_t *x,Double_t *par)
{return (par[13]*x[0]*x[0] + par[14]*x[0] + par[15]) ;}

Double_t SignalFcnJPsiV2(Double_t *x,Double_t *par)
{return (JPsiCrystalBallExtended(x,par)*par[12])/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t FourierV2(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);}

Double_t FourierV5(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);}// + 2*par[3]*cos(3*x[0]);}// + 2*par[4]*cos(4*x[0]);}// + 2*par[5]*cos(5*x[0]) + 2*par[6]*cos(6*x[0]) + 2*par[7]*cos(7*x[0]) + 2*par[8]*cos(8*x[0]) + 2*par[9]*cos(9*x[0]) + 2*par[10]*cos(10*x[0]) + 2*par[11]*cos(11*x[0]);}

Double_t TwoCBE2E(Double_t *x,Double_t *par)

{ return JPsiCrystalBallExtended(x,par) + Psi2SCrystalBallExtended(x,par) + ExpBkg(x,&par[8]);}

Double_t ExpBkg(Double_t *x,Double_t *par)

{ return par[0]*(exp(x[0]*par[1]*(-1))) + par[2]*(exp(x[0]*par[3]*(-1))); }

Double_t JPsiCrystalBallExtended(Double_t *x,Double_t *par)
{
    
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
    
    return sum ;
}

Double_t Psi2SCrystalBallExtended(Double_t *x,Double_t *par)
{
      Double_t sum = 0;
    
    Double_t t = (x[0]-(par[1]*mPsip/mJpsi))/(par[2]*ratSigma);
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
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
     7: Normalisation Psi2s, JPsi
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
   TF1 *fitFcn = new TF1("fitFcn",TwoCBE2E,2.1,5.1,12);
   fitFcn->SetNpx(500);
   fitFcn->SetLineWidth(4);
   fitFcn->SetLineColor(kMagenta);
   // first try without starting values for the parameters
   // This defaults to 1 for each param.
   // this results in an ok fit for the polynomial function
   // however the non-linear part (lorenzian) does not
   // respond well.
    Double_t params[12] = {10,1,1,1,1,1,1,1,1,1,1,1};
   fitFcn->SetParameters(params);
    TVirtualFitter::Fitter(histo)->SetMaxIterations(100000);
    TVirtualFitter::Fitter(histo)->SetPrecision();
 //  histo->Fit("fitFcn","0");
   // First fit, fix MJPsi and Sigma
    
    fitFcn->SetParLimits(0,0.1,100000000);
    fitFcn->SetParLimits(1,3.05,3.15);
    fitFcn->SetParLimits(2,0.03,0.1);
    fitFcn->SetParLimits(3,0.7,1.1);
    fitFcn->SetParLimits(4,1,30);
    fitFcn->SetParLimits(5,1,5);
    fitFcn->SetParLimits(6,5,25);
    fitFcn->SetParLimits(7,0.001,1000000);
    fitFcn->SetParLimits(8,0.1,1000000);
    fitFcn->SetParLimits(9,0.01,20);
    fitFcn->SetParLimits(10,0.01,1000000);
    fitFcn->SetParLimits(11,0.01,50);
    
   fitFcn->FixParameter(1,3.096916); // Mean x core
    fitFcn->FixParameter(2,0.07);
    fitFcn->FixParameter(3,0.9);
    fitFcn->FixParameter(4,10);
    fitFcn->FixParameter(5,2);
    fitFcn->FixParameter(6,15);
    fitFcn->SetParameter(7,100);
    fitFcn->SetParameter(8,1000);
    fitFcn->SetParameter(9,0.5);
    fitFcn->SetParameter(11,10);

    
    fitFcn->SetParName(0,"Norm_{JPsi}");
    fitFcn->SetParName(1,"M_{JPsi}");
    fitFcn->SetParName(2,"Sigma_{JPsi}");
    fitFcn->SetParName(3,"a_{1}");
    fitFcn->SetParName(4,"n_{1}");
    fitFcn->SetParName(5,"a_{2}");
    fitFcn->SetParName(6,"n_{2}");
    fitFcn->SetParName(7,"Norm_{Psi2S}");
    fitFcn->SetParName(8,"Norm_{TailLowM}");
    fitFcn->SetParName(9,"Exp_{TailLowM}");
    fitFcn->SetParName(10,"Norm_{TailHighM}");
    fitFcn->SetParName(11,"Exp_{TailHighM}");
    
   TFitResultPtr res = histo->Fit("fitFcn","SBMER","ep");
   res = histo->Fit("fitFcn","SBMERQ","ep");
      fitFcn->ReleaseParameter(1);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    fitFcn->ReleaseParameter(2);
   res = histo->Fit("fitFcn","SBMERQ","ep");
    fitFcn->ReleaseParameter(7);
   res = histo->Fit("fitFcn","SBMERQ","ep");
   res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(3);
//    fitFcn->ReleaseParameter(4);
//   res = histo->Fit("fitFcn","SBMER","ep");
//   res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(5);
//    fitFcn->ReleaseParameter(6);
//   res = histo->Fit("fitFcn","SBMER","ep");
//    res = histo->Fit("fitFcn","SBMER","ep");

    
   // improve the picture:
   TF1 *backFcn = new TF1("backFcn",ExpBkg,2.1,5.1,4);
   backFcn->SetLineColor(kRed);
   TF1 *signalFcnJPsi = new TF1("signalFcnJPsi",JPsiCrystalBallExtended,2.1,5.1,8);
   TF1 *signalFcnPsi2S = new TF1("signalFcnPsi2S",Psi2SCrystalBallExtended,2.1,5.1,8);
   TPaveText *pave = new TPaveText(0.15,0.5,0.3,0.65,"brNDC");
   signalFcnJPsi->SetLineColor(kBlue);
   signalFcnJPsi->SetNpx(500);
    signalFcnPsi2S->SetLineColor(kGreen);
    signalFcnPsi2S->SetNpx(500);
   Double_t par[12];
   // writes the fit results into the par array
    gStyle->SetOptFit(1011);
   fitFcn->GetParameters(par);
   signalFcnJPsi->SetParameters(par);
    signalFcnPsi2S->SetParameters(par);
//   auto covMatrix = res->GetCovarianceMatrix();
//   std::cout << "Covariance matrix from the fit ";
//   covMatrix.Print();
   Double_t integral = (signalFcnJPsi->Integral(2.1,3.45))/0.01;
    auto covtot = res->GetCovarianceMatrix();
    auto covsgn = covtot.GetSub(0,8,0,8);
    std::cout << "Matrice totale" <<endl;
    covtot.Print();
    std::cout << "Matrice réduite" <<endl;
    covsgn.Print();
    std::cout << "STATUS COV " << res->CovMatrixStatus() <<endl;
   Double_t integralerror = (signalFcnJPsi->IntegralError(2.1,3.45,signalFcnJPsi->GetParameters(), res->GetCovarianceMatrix().GetSub(0,8,0,8).GetMatrixArray() ))/0.01;
    std::cout << "Erreur integrale " << integralerror <<endl;
    
    std::cout << "Fitted " << histoname << std::endl;
    std::cout << "Nb JPsi total: " << integral << std::endl;
 //   std::cout << "integral error: " << integralerror << std::endl;
   signalFcnJPsi->Draw("same");
    signalFcnPsi2S->Draw("same");
   backFcn->SetParameters(&par[8]);
   backFcn->Draw("same");
   // draw the legend
    char str[50];
    sprintf(str, "N_{JPsi} %f +/- %f", integral, integralerror);
   pave->AddText(str);
    sprintf(str, "M_{Psi2S} = %f, Sig_{Psi2S} = %f", mPsip, ratSigma*par[2]);
    pave->AddText(str);
   pave->Draw();
   TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
   legend->SetTextFont(72);
   legend->SetTextSize(0.03);
    Char_t message[80];
    sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcn->GetChisquare(),fitFcn->GetNDF());
     legend->AddEntry(fitFcn,message);
    if(res->CovMatrixStatus() == 3){
               sprintf(message,"The fit is a success");
           }
           else{
               sprintf(message,"The fit is a failure");
           }
           legend->AddEntry(fitFcn,message);
   legend->AddEntry(histo,"Data","lpe");
   legend->AddEntry(backFcn,"Background fit","l");
   legend->AddEntry(signalFcnJPsi,"JPsi Signal fit","l");
    legend->AddEntry(signalFcnPsi2S,"Psi2S Signal fit","l");
   legend->Draw();
    
//    cout << "Histogramme: " << histoname << endl;
//    for(int j=0; j<10; j++){
//        cout << "Bin de masse numero " << j << " - Signal : " << (signalFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 << " et Background : " << (backFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 <<endl;
//    }
    
    return res;
}

void FcnCombinedAllMass(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  TAxis *xaxis1  = hnseg->GetXaxis();
  TAxis *xaxis2  = V2JPsiTkl->GetXaxis();

  int nbinX1 = hnseg->GetNbinsX();
  int nbinX2 = V2JPsiTkl->GetNbinsX();

  double chi2 = 0;
  double x[1];
  double tmp;
  npfits = 0;
  for (int ix = 1; ix <= nbinX1; ++ix) {
    x[0] = xaxis1->GetBinCenter(ix);
      if ( hnseg->GetBinError(ix) > 0 ) {
        tmp = (hnseg->GetBinContent(ix) - TwoCBE2E(x,p))/hnseg->GetBinError(ix);
        chi2 += tmp*tmp;
        npfits++;
      }
  }
  for (int ix = 1; ix <= nbinX2; ++ix) {
     x[0] = xaxis2->GetBinCenter(ix);
      if ( V2JPsiTkl->GetBinError(ix) > 0 ) {
        tmp = (V2JPsiTkl->GetBinContent(ix) - FourierV2_WrtInvMass(x,p))/V2JPsiTkl->GetBinError(ix);
        chi2 += tmp*tmp;
        npfits++;
      }
  }
  fval = chi2;
}

int GetCent(double cent){
    if(cent <= 1){
        return 0;
    }
    else if(cent <= 10){
        return 1;
    }
    else if(cent <= 20){
        return 2;
    }
    else if(cent <= 40){
        return 3;
    }
    else if(cent <= 100){
        return 4;
    }
//    else if(cent <= 30){
//        return 6;
//    }
//    else if(cent <= 40){
//        return 7;
//    }
//    else if(cent <= 50){
//        return 8;
//    }
//    else if(cent <= 60){
//        return 9;
//    }
//    else if(cent <= 70){
//        return 10;
//    }
//    else if(cent <= 80){
//        return 11;
//    }
//    else if(cent <= 90){
//        return 12;
//    }
//    else{
//        return 13;
//    }
}

int GetCentPM(int Ntkl, int zvtx_idx, int groupnumber){

    if(Ntkl==0){
        return 12;
    }
    
//    if(Ntkl>40){
//        return 0;
//    }

    for(int cent_index=0;cent_index<13;cent_index++){
        if(LimitsPM[groupnumber-1][zvtx_idx][cent_index] < Ntkl){
            return cent_index;
        }
    }
}

