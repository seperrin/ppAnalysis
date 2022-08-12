 #include <TBrowser.h>
 #include <TBufferFile.h>
 #include <TCanvas.h>
 #include <TClass.h>
 #include <TPad.h>
 #include <TMath.h>
 #include <TLatex.h>
 #include <TObject.h>
 #include <TClonesArray.h>
 #include <TFile.h>
 #include <TTree.h>
 #include <TFitResultPtr.h>
 #include <TPaveText.h>
 #include <TGraph.h>
 #include <TGraphErrors.h>
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
# include "TPaveStats.h"
 # include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"

// Re-tourner le post-processing et les fits à partir de données déjà enregistrées

 // FUNCTIONS

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

void Runner(Char_t radical[500], bool isFull);

 void FitTrainingPtBinned(Char_t radical[500], string SignalF, string BackgroundF, double minmass, double maxmass, double ratsigma, double minv2, double maxv2, string BackgroundV2F);
 Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par);
 Double_t CvetanF(Double_t *x,Double_t *par);
void ChisquareCvetanF(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t CvetanFPtBinned(Double_t *x,Double_t *par);
void ChisquareCvetanFPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t CvetanFTKL(Double_t *x,Double_t *par);
Double_t ZYAM(Double_t *x,Double_t *par);
void ChisquareZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t ZYAMPtBinned(Double_t *x,Double_t *par);
void ChisquareZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate(Double_t *x,Double_t *par);
void ChisquarePRLTemplate(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplatePtBinned(Double_t *x,Double_t *par);
void ChisquarePRLTemplatePtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplatePtBinnedAlt(Double_t *x,Double_t *par);
void ChisquarePRLTemplatePtBinnedAlt(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate_RidgeAndZero(Double_t *x,Double_t *par);
Double_t PRLTemplate_RidgeAndZeroPtBinned(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphAndG(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphAndGPtBinned(Double_t *x,Double_t *par);
Double_t PRLTemplate_PeriphZYAM(Double_t *x,Double_t *par);
void ChisquarePRLTemplate_PeriphZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
Double_t PRLTemplate_PeriphZYAMPtBinned(Double_t *x,Double_t *par);
void ChisquarePRLTemplate_PeriphZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  );
 Double_t BackFcnV2(Double_t *x,Double_t *par);
Double_t BackFcnV2Poly(Double_t *x,Double_t *par);
 Double_t SignalFcnJPsiV2(Double_t *x,Double_t *par);
 Double_t FourierV2(Double_t *x,Double_t *par);
 Double_t FourierV3(Double_t *x,Double_t *par);
 Double_t FourierV5(Double_t *x,Double_t *par);
 Double_t MassFunction(Double_t *x,Double_t *par);
 Double_t ExpBkg(Double_t *x,Double_t *par);
 Double_t Pol1Pol2(Double_t *x,Double_t *par);
 Double_t VWGaussian(Double_t *x,Double_t *par);
 Double_t Tchebychev(Double_t *x,Double_t *par);
 Double_t SmolTcheby(Double_t x,Double_t *par, int n);
Double_t ExpBkgHole(Double_t *x,Double_t *par);
Double_t Pol1Pol2Hole(Double_t *x,Double_t *par);
Double_t VWGaussianHole(Double_t *x,Double_t *par);
Double_t TchebychevHole(Double_t *x,Double_t *par);
 Double_t JPsiCrystalBallExtended(Double_t *x,Double_t *par);
 Double_t Psi2SCrystalBallExtended(Double_t *x,Double_t *par);
Double_t JPsiNA60(Double_t *x,Double_t *par);
Double_t Psi2SNA60(Double_t *x,Double_t *par);
Double_t UsedJPsiSignal(Double_t *x,Double_t *par);
Double_t UsedPsi2SSignal(Double_t *x,Double_t *par);
Double_t UsedBackground(Double_t *x,Double_t *par);
Double_t UsedBackgroundHole(Double_t *x,Double_t *par);
// TFitResultPtr FittingAllInvMass(const char *histoname, TCanvas *canvas);
 TFitResultPtr FittingAllInvMassBin(const char *histoname, TCanvas *cinvmass, int i, string SignalF, string BackgroundF, double min, double max, double ratsigma);
 void FcnCombinedAllMass(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  );
 int GetCent(double cent);
int NumberOfParameters(string SignalF, string BackgroundF);
int NumberOfParametersBkg(string BackgroundF);
int NumberOfParametersBkgV2(string BackgroundV2F);
string GetParameterInfo(int i, string SignalF, string BackgroundF, string BackgroundV2F, string prefix);

Double_t mJpsi =  3.096916;
 Double_t mPsip =  3.686108;
// Double_t ratMass = 1.01;
 Double_t ratSigma = 1.05;
Double_t mDiff = 0.589188;
Int_t npfits;

Int_t actuelptbin = 0;

double MassBins[] = {1.0,1.5,1.9,2.3,2.7,3.0,3.3,3.7,4.1,4.5,5.0};
const int NbinsInvMass = 10;
double MinInvMass = 1.;
double MaxInvMass = 5.;
double MinInvMassFit = 2; //2
double MaxInvMassFit = 5.; //4.5
double MinV2Fit = 1.;
double MaxV2Fit = 5.;

const int NbinsDimuInvMass = 400; //300

int numLoops = 0;

int numParameters;
int numParametersBkg;
int numParametersBkgV2;

TH1I* hnseg(NULL);
TH1F* hpool(NULL);
TH1F* V2JPsiTkl(NULL);

TH1F* Yields_Central_1(NULL);
TH1F* Yields_Periph_1(NULL);
TH1F* Baseline_Central_1(NULL);
TH1F* Baseline_Periph_1(NULL);
TH1F* Yields_Periph_1_MinusBaseline(NULL);
TH1F* Yields_Central_1_MinusBaseline(NULL);

TH1F* YieldTkl_Central(NULL);
TH1F* YieldTkl_Periph(NULL);
TH1F* BaselineTkl_Central(NULL);
TH1F* BaselineTkl_Periph(NULL);
TH1F* YieldTkl_Periph_MinusBaseline(NULL);
TH1F* YieldTkl_Central_MinusBaseline(NULL);

double baseline_periph = 1;
double errbaseline_periph = 1;
double baseline_central = 99999;
double errbaseline_central = 99999;

double baselineTKL_periph = 1;
double errbaselineTKL_periph = 1;
double baselineTKL_central = 9999;
double errbaselineTKL_central = 1;

double PtBins[] = {0,2,3,4,6,8,12};
int NbPtBins = 6;
double LowDimuPtCut = 0;
double HighDimuPtCut = 12;

TH1F* V2JPsiTklPtBinned[10]{ NULL };
TH1F* V22ATLASPtBinned[10]{ NULL };
TH1F* FATLASPtBinned[10]{ NULL };
TH1 *hPtWrtMassInvSliced[3][10]= {NULL};

TH1F* Yield_Central_MassBinPtBinned[10][10]= {NULL};
TH1F* Yield_Periph_MassBinPtBinned[10][10]= {NULL};

double V2_Ext1[10] = {0};
double V2_Ext1_V3[10] = {0};
double V2_Ext2[10] = {0};
double V2_Ext3[10] = {0};
double V2_Ext1_noZYAM[10] = {0};
double V2_Ext2_noZYAM[10] = {0};
double V2_Ext3_noZYAM[10] = {0};
double V2_CvetanQuentin[10] = {0};
double V2_CvetanQuentinMe[10] = {0};
double V2_ZYAM[10] = {0};
double V2_PRL[10] = {0};
double V2_PRLPeriphZYAM[10] = {0};
double errV2_Ext1[10] = {0};
double errV2_Ext1_V3[10] = {0};
double errV2_Ext2[10] = {0};
double errV2_Ext3[10] = {0};
double errV2_Ext1_noZYAM[10] = {0};
double errV2_Ext2_noZYAM[10] = {0};
double errV2_Ext3_noZYAM[10] = {0};
double errV2_CvetanQuentin[10] = {0};
double errV2_CvetanQuentinMe[10] = {0};
double errV2_ZYAM[10] = {0};
double errV2_PRL[10] = {0};
double errV2_PRLPeriphZYAM[10] = {0};
double V2_ATLASExt2[10] = {0};
double errV2_ATLASExt2[10] = {0};

double F_CvetanQuentin[10] = {0};
double F_PRL[10] = {0};
double F_PRLPeriphZYAM[10] = {0};
double errF_CvetanQuentin[10] = {0};
double errF_PRL[10] = {0};
double errF_PRLPeriphZYAM[10] = {0};


TH1F* Yields_Central_1PtBinned[10] = {NULL};
TH1F* Yields_Periph_1PtBinned[10] = {NULL};
TH1F* Baseline_Central_1PtBinned[10] = {NULL};
TH1F* Baseline_Periph_1PtBinned[10] = {NULL};
TH1F* Yields_Periph_1_MinusBaselinePtBinned[10] = {NULL};
TH1F* Yields_Central_1_MinusBaselinePtBinned[10] = {NULL};

double baseline_periphPtBinned[10] = {NULL};
double errbaseline_periphPtBinned[10] = {NULL};
double baseline_centralPtBinned[10] = {NULL};
double errbaseline_centralPtBinned[10] = {NULL};

Char_t FitFileName[500];
string VarSignals[3] = {"CB2-Run2", "CB2-MC","NA60-MC"};
string VarBackgrounds[4] = {"DoubleExpo", "POL1POL2", "VWG","Tchebychev"};
double VarMinInvMass[3] = {2., 2.3, 2.4};
double VarMaxInvMass[3] = {5., 4.9, 4.7};
double VarRatSigma[2] = {1.05, 1.0};

double VarMinV2fit[2] = {1.,1.5};
double VarMaxV2fit[2] = {5.,4.5};
string VarV2Backgrounds[3] = {"Pol2","Pol1","Pol3"};

double MinDeltaPhi;
double MaxDeltaPhi;
const int NbinsDeltaPhi = 6;

string GlobalSignal = "";
string GlobalBackground = "";
double GlobalMinMass;
double GlobalMaxMass;
double GlobalRatSigma;
double GlobalMinV2;
double GlobalMaxV2;
string GlobalBackgroundV2 = "";

void Runner(Char_t radical[500], bool isFullStudy = kFALSE){
    
    if(!isFullStudy){
    cout << "=== Will run FitTrainingPtBinned default ===" <<endl;
        FitTrainingPtBinned(radical, "CB2-Run2", "DoubleExpo", 2., 5., 1.05, 1., 5., "Pol2");
    }
    if(isFullStudy){
        cout << "=== Will run FitTrainingPtBinned default and fit variants ===" <<endl;
        for(int signal=2; signal<3; signal++){//3
            for(int background=2; background<4; background++){//4
                for(int massrange=0; massrange<3; massrange++){//3
                    for(int ratsigma=0; ratsigma<2; ratsigma++){//2
                        for(int v2range=0; v2range<1; v2range++){//2
                            for(int v2bkg=0; v2bkg<1; v2bkg++){//3
                                FitTrainingPtBinned(radical, VarSignals[signal], VarBackgrounds[background], VarMinInvMass[massrange], VarMaxInvMass[massrange], VarRatSigma[ratsigma], VarMinV2fit[v2range], VarMaxV2fit[v2range], VarV2Backgrounds[v2bkg]);
                                numLoops+=1;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return;
}

void FitTrainingPtBinned(Char_t radical[500], string SignalF, string BackgroundF, double minmass, double maxmass, double ratsigma, double minv2, double maxv2, string BackgroundV2F){
    
     numParameters = NumberOfParameters(SignalF, BackgroundF);
     numParametersBkg = NumberOfParametersBkg(BackgroundF);
     numParametersBkgV2 = NumberOfParametersBkgV2(BackgroundV2F);
    
     GlobalSignal = SignalF;
     GlobalBackground = BackgroundF;
     GlobalMinMass = minmass;
     GlobalMaxMass = maxmass;
     GlobalRatSigma = ratsigma;
     GlobalMinV2 = minv2;
     GlobalMaxV2 = maxv2;
     GlobalBackgroundV2 = BackgroundV2F;

    
    bool isCentralityStudy = kFALSE;
    cout << "Preparing the analysis of " << radical<<endl;
    cout << "isCentralityStudy ? : " << isCentralityStudy <<endl;
    
    if(isCentralityStudy){
        cout << "ATTENTION: LES RESULTATS UTILISANT LES FOURIERS (Ext2 et Ext3) SONT INUTILISABLES"<<endl;
    }
    
    MinInvMassFit = minmass; //2
    MaxInvMassFit = maxmass; //4.5
    MinV2Fit = minv2;
    MaxV2Fit = maxv2;
    
  //  sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile1_Run2_0-5_40-90_pt0-2-4-6-8-12.root");
    //sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFiles/Septembre2021-Run2ReferencesTKLSystematicsStart/FitFile_CrossCheckQGPFrance_Run2_PercentileMethodSPDTracklets_0-5_40-100_pt0-2-4-6-8-12.root");
   // sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysis_Run2_QGPFrance_0-5_40-100_pt0-2-4-6-8-12.root"); //BEST
  //  sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_Run2_V0MPercentile_0-5_40-100_pt0-2-3-4-6-8-12.root");
    //sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_Run2_V0MPercentile_0-5_40-100_pt0-2-3-4-6-8-12.root");
   // sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysis_16h_QGPFrance_0-5_40-100_pt0-4-12_Test.root");
    //sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFiles/FinAvril2021-PUfixedMUTCutsfixedButBadEtaCutsNoPhiSPDCorr/FitFile_GoodPU_Run2_0-5_40-100_pt0-2-4-6-8-12.root");
  //  sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_GoodPU_Gr1_0-5_40-100_pt8-12.root");
    
    Char_t RadicalName[500];
    Char_t FolderName[500];
    Char_t CanvasFolderName[500];
    Char_t FitCentralFileName[500];
    Char_t FitPeriphFileName[500];
    Char_t SystematicsFileName[500];
    Char_t CanvasName[500];
    
    sprintf(RadicalName,"%s",radical);
   
    if(!isCentralityStudy){
        sprintf(FitFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_%s.root",RadicalName);
    }
    
    sprintf(FolderName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrierDimu/%s", radical);
    
    sprintf(SystematicsFileName,"%s/SystematicsFile.csv", FolderName);
    
    if(SignalF == VarSignals[0] && BackgroundF == VarBackgrounds[0] && minmass == VarMinInvMass[0] && maxmass == VarMaxInvMass[0] && ratsigma == VarRatSigma[0] && minv2 == VarMinV2fit[0] && maxv2 == VarMaxV2fit[0] && BackgroundV2F == VarV2Backgrounds[0]){
    sprintf(CanvasFolderName,"%s",FolderName);
    }
    else{
            sprintf(SystematicsFileName,"%s/SystematicsFile_%s_%s_%.2f-%.2f_%.2f_%.2f-%.2f_%s.csv", FolderName,SignalF.c_str(), BackgroundF.c_str(), minmass, maxmass, ratsigma, minv2, maxv2, BackgroundV2F.c_str());
        

        sprintf(CanvasFolderName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrierDimu/%s/%s_%s_%.2f-%.2f_%.2f_%.2f-%.2f_%s", radical, SignalF.c_str(), BackgroundF.c_str(), minmass, maxmass, ratsigma, minv2, maxv2, BackgroundV2F.c_str());
        int rc = mkdir(CanvasFolderName, 0777);
        cout << "rc " <<rc<<endl;
    }
    
    
    // Getting info from radical name
    
    
    Char_t DPhiCutText[20] = "10mrad";
    Char_t CentEstText[50] = "PercentileMethodSPDTracklets";
    Char_t CentralClassText[20] = "0-5";
    Char_t PeriphClassText[20] = "40-100";
    Char_t DeltaEtaGapText[20] = "1.5";
    Char_t DeltaEtaGapMaxText[20] = "5.0";
    Char_t ZvtxCutText[20] = "10";
    
    
    Char_t EMNormText[20] = "Method1";
    
    
    Char_t EMPoolMaxText[20] = "100";
    Char_t EMPoolThresholdText[20] = "10";
    Char_t ScalingPeriphText[20] = "1";
    Char_t EMChangeText[20] = "None";
    
    
    Char_t SummationText[20] = "Method1c";
    
    bool hasChangedDPhiCut = kFALSE;
    bool hasChangedEtaGap = kFALSE;
    bool hasChangedEtaGapMax = kFALSE;
    bool hasChangedZvtxCut = kFALSE;
    bool hasChangedEMNorm = kFALSE;
    bool hasChangedEMMax = kFALSE;
    bool hasChangedEM = kFALSE;
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
    Char_t binspt[30];
    Char_t binpt[30];
    
    double DPhiCutNum;
    double DeltaEtaGapNum;
    double DeltaEtaGapMaxNum;
    double ZvtxCutNum;
    
    sscanf(radical, "NewAnalysisAllEst_%[^_]_%[^_]_%[^-]-%[^_]_%[^-]-%[^_]_pt%[^-]-%s", dataUsed, estimator, mostCentral,lessCentral, lessPeriph, mostPeriph, minpt, rest);
    
    cout << "Radical is: "<< radical<<endl;
    cout << "dataUsed is: "<< dataUsed<<endl;
    cout << "estimator is: "<< estimator<<endl;
    cout << "mostCentral is: "<< mostCentral<<endl;
    cout << "lessCentral is: "<< lessCentral<<endl;
    cout << "lessPeriph is: "<< lessPeriph<<endl;
    cout << "mostPeriph is: "<< mostPeriph<<endl;
    cout << "minpt is: "<< minpt<<endl;
    cout << "rest is: "<< rest<<endl;
    
    
    string srest;
    stringstream ss;
    ss << rest;
    ss >> srest;
    
    PtBins[0] = stoi(minpt);
    
    cout << "Set PtBins[0] to minpt: PtBins[0] = " << PtBins[0]<<endl;
    cout << "Looking into rest:"<<endl;

    if (srest.find("_") != std::string::npos) {
        std::cout << "Underscore found: There is a specificity" << endl;
        sscanf(rest,"%[^_]_%s", binspt, specificity);
    }
    else{
        cout << "No underscore found: No specificity"<<endl;
        sprintf(binspt, "%s", rest);
    }
    
    string sbinspt;
    stringstream ssbins;
    ssbins << binspt;
    ssbins >> sbinspt;
    
    cout << "binspt is: " <<binspt<<endl;
    
    int nbptbins = 0;

    if (sbinspt.find("-") != std::string::npos) {
        while(sbinspt.find("-") != std::string::npos){
            nbptbins++;
            std::cout << "Tiret found: There is another ptbin" << endl;
            sscanf(binspt,"%[^-]-%s", binpt, binspt);
            PtBins[nbptbins] = stoi(binpt);
            cout << "For bin " <<nbptbins << " pt limit is "<<PtBins[nbptbins]<<endl;
            stringstream ssbinsT;
            ssbinsT << binspt;
            ssbinsT >> sbinspt;
            
            cout << "sbinspt is: " << sbinspt<<endl;
            //Convert binpt to int and addit to an array
        }
    }
    if (sbinspt.find("-") == std::string::npos) {
        nbptbins++;
        cout << "No tiret found: No new ptbin"<<endl;
        if (sbinspt.find("_") != std::string::npos) {
            std::cout << "Underscore found: There is a specificity" << endl;
            sscanf(binspt,"%[^_]_%s", binspt, specificity);
        }
        else{
            cout << "No underscore found: No specificity"<<endl;
        }
        PtBins[nbptbins] = stoi(binspt);
        cout << "For bin " <<nbptbins << " pt limit is "<<PtBins[nbptbins]<<endl;
    }
    
    NbPtBins = nbptbins;
    double LowDimuPtCut = PtBins[0];
    double HighDimuPtCut = PtBins[nbptbins];
    
    cout << "NbPtBins is: " <<NbPtBins<<endl;
    cout << "LowDimuPtCut is: " <<LowDimuPtCut<<endl;
    cout << "HighDimuPtCut is: " <<HighDimuPtCut<<endl;
    
    
    string sspe;
    stringstream ss2;
    ss2 << specificity;
    ss2 >> sspe;
    
    cout << "Looking at specificity: "<<specificity<<endl;
    
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
        
        if (sspecT.find("EtaMin") != std::string::npos) {
            std::cout << "Change on EtaMin" << endl;
           sscanf(speTampon, "EtaMin%s", DeltaEtaGapText);
           cout << "DeltaEtaGapText set to " << DeltaEtaGapText<<endl;
            hasChangedEtaGap = kTRUE;
        }
        
        if (sspecT.find("EtaMax") != std::string::npos) {
            std::cout << "Change on EtaMax" << endl;
           sscanf(speTampon, "EtaMax%s", DeltaEtaGapMaxText);
           cout << "DeltaEtaGapMaxText set to " << DeltaEtaGapMaxText<<endl;
            hasChangedEtaGapMax = kTRUE;
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
        
        if (sspecT.find("EMChange") != std::string::npos) {
            std::cout << "Precisions on EMChange" << endl;
           sscanf(speTampon, "EMChange%s", EMChangeText);
           cout << "EMChangeText set to " << EMChangeText<<endl;
            hasChangedEM = kTRUE;
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
    
    if (sspe.find("EtaMin") != std::string::npos) {
        std::cout << "Change on EtaMin" << endl;
       sscanf(specificity, "EtaMin%s", DeltaEtaGapText);
       cout << "DeltaEtaGapText set to " << DeltaEtaGapText<<endl;
        hasChangedEtaGap = kTRUE;
    }
    
    if (sspe.find("EtaMax") != std::string::npos) {
        std::cout << "Change on EtaMax" << endl;
       sscanf(specificity, "EtaMax%s", DeltaEtaGapMaxText);
       cout << "DeltaEtaGapMaxText set to " << DeltaEtaGapMaxText<<endl;
        hasChangedEtaGapMax = kTRUE;
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
    
    if (sspe.find("EMChange") != std::string::npos) {
        std::cout << "Precisions on EMChange" << endl;
       sscanf(specificity, "EMChange%s", EMChangeText);
       cout << "EMChangeText set to " << EMChangeText<<endl;
        hasChangedEM = kTRUE;
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
    cout << "min pt " << minpt<<endl;
    
    // Là on a toutes les linformations de radical et on sait qui a été changé
    
    //MAJ des informations sur la centralité (estimateur et classes)
    
    int CentralLowBound = GetCent(mostCentralNum);
    int CentralHighBound = GetCent(lessCentralNum);
    int PeripheralLowBound = GetCent(lessPeriphNum)+1;
    int PeripheralHighBound = GetCent(mostPeriphNum);
    
    
    double ZvtxCut = stod(ZvtxCutText);

    double DeltaEtaDimuCut = stod(DeltaEtaGapText); //1.2  //1.2
    double MaxDeltaEta = stod(DeltaEtaGapMaxText);
     
    int NbBinsZvtx = int(2*ZvtxCut);
    
    int EMPoolMax = stoi(EMPoolMaxText);
    int EMPoolThreshold = stoi(EMPoolThresholdText);
    double ScalingPeriph = stod(ScalingPeriphText);
    
    
    if(isCentralityStudy){
        sscanf(radical, "NewAnalysisAllEst_%[^_]_%[^_]_%[^-]-%[^_]_%[^-]-%[^_]_pt%[^-]-%s", dataUsed, estimator, mostCentral,lessCentral, lessPeriph, mostPeriph, minpt, rest);
        cout << "For centrality study, classes are: " << mostCentral << "-" << lessCentral <<" and " << lessPeriph << "-" << mostPeriph<<endl;
        if(lessCentralNum == 1){
            sprintf(FitCentralFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_%s-%s_20-100_pt%s-%s.root", dataUsed, estimator, mostCentral,lessCentral, minpt, rest);
        }
        if(lessCentralNum == 3 && mostCentralNum == 0){
            sprintf(FitCentralFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts//FitFile_NewAnalysisAllEst_%s_%s_%s-%s_60-100_pt%s-%s.root", dataUsed, estimator, mostCentral,lessCentral, minpt, rest);
        }
        if(lessCentralNum == 3 && mostCentralNum == 1){
            sprintf(FitCentralFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts//FitFile_NewAnalysisAllEst_%s_%s_%s-%s_40-80_pt%s-%s.root", dataUsed, estimator, mostCentral,lessCentral, minpt, rest);
        }
        if(lessCentralNum == 5 && mostCentralNum == 0){
            sprintf(FitCentralFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts//FitFile_NewAnalysisAllEst_%s_%s_%s-%s_40-100_pt%s-%s.root", dataUsed, estimator, mostCentral,lessCentral, minpt, rest);
        }
        if(lessCentralNum == 5 && mostCentralNum == 3){
            sprintf(FitCentralFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts//FitFile_NewAnalysisAllEst_%s_%s_%s-%s_60-80_pt%s-%s.root", dataUsed, estimator, mostCentral,lessCentral, minpt, rest);
        }
        if(lessCentralNum == 10 && mostCentralNum == 0){
            sprintf(FitCentralFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts//FitFile_NewAnalysisAllEst_%s_%s_%s-%s_20-80_pt%s-%s.root", dataUsed, estimator, mostCentral,lessCentral, minpt, rest);
        }
        if(lessCentralNum == 10 && mostCentralNum == 5){
            sprintf(FitCentralFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts//FitFile_NewAnalysisAllEst_%s_%s_%s-%s_0-100_pt%s-%s.root", dataUsed, estimator, mostCentral,lessCentral, minpt, rest);
        }
        if(lessPeriphNum == 20 && mostPeriphNum == 100){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_0-1_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 30 && mostPeriphNum == 100){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_0-5_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 50 && mostPeriphNum == 100){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_0-5_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 60 && mostPeriphNum == 100){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_0-3_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 20 && mostPeriphNum == 80){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_0-10_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 40 && mostPeriphNum == 80){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_1-3_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 60 && mostPeriphNum == 80){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_3-5_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 0 && mostPeriphNum == 100){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_5-10_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        if(lessPeriphNum == 40 && mostPeriphNum == 100){
            sprintf(FitPeriphFileName,"~/../../Volumes/Transcend2/ppAnalysis/Scripts/FitFile_NewAnalysisAllEst_%s_%s_0-5_%s-%s_pt%s-%s.root", dataUsed, estimator, lessPeriph, mostPeriph, minpt, rest);
        }
        
        sprintf(FitFileName,"%s",FitCentralFileName);
    }
    
    
    TH1::SetDefaultSumw2();
        bool doTracklets = kFALSE;
        bool doMixedEvents = kTRUE;
    bool Extraction2Combined = kFALSE;
        bool CombineFits = kFALSE;
    
    bool PtBinned = kTRUE;
        
    // ************************************
    // Définitions de paramètres          *
    // ************************************

        //double ZvtxCut = 10;
        double SigmaZvtxCut = 0.25;
        double DPhiCut = 0.01;
        if(hasChangedDPhiCut){
            DPhiCut = DPhiCutNum/1000.;
        }
    
        double TklEtaCut  = 1;
        double LowDimuYCut = -4;
        double HighDimuYCut = -2.5;
       // double LowDimuPtCut = 0;
        //double HighDimuPtCut = 99999;
    //    double MinMultCentral = 37;
    //    double MaxMultPeriph = 23;
        double CentSPDTrackletsCentral = 1;
        double CentSPDTrackletsPeriph = 40;

        double LowJPsiMass = 3.0;
        double HighJPsiMass = 3.3;

//        const int NbinsInvMass = 10;
//        double SizeBinInvMass = (MaxInvMass-MinInvMass)/NbinsInvMass;

    MinDeltaPhi = 0;//-TMath::Pi()/2;
    MaxDeltaPhi = TMath::Pi();;//1.5*TMath::Pi();
    //12;
        double SizeBinDeltaPhi = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhi;
        
        double MinDeltaEta = 0;
       // double MaxDeltaEta = 6;
        int NbinsDeltaEta = int(2*MaxDeltaEta);
        double SizeBinDeltaEta = (MaxDeltaEta-MinDeltaEta)/NbinsDeltaEta;
        
        double MinDeltaEtaTKL = 1.2;
        double MaxDeltaEtaTKL = 2.4;
        const int NbinsDeltaEtaTKL = 24;
        double SizeBinDeltaEtaTKL = (MaxDeltaEtaTKL-MinDeltaEtaTKL)/NbinsDeltaEtaTKL;
        
        const int NbinsDeltaPhiTKL = 24;
        double SizeBinDeltaPhiTKL = (MaxDeltaPhi-MinDeltaPhi)/NbinsDeltaPhiTKL;
    int BinZeroLeftTKL = 1;//floor((0-MinDeltaPhi)*NbinsDeltaPhiTKL/(2*TMath::Pi()));
        
        const int NbBinsCent = 14;
       // const int NbBinsZvtx = 20;
    
    // *************************
    // Initialiser les graphes *
    // *************************
    
    TH1F* YieldValue(NULL);
    TH1F* YieldError(NULL);
    
    YieldValue = new TH1F("YieldValue",
    "YieldValue",
    10000,0,100);
    YieldError = new TH1F("YieldError",
    "YieldError",
    10000,0,100);

        TH1F* Yield_Central_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yield_Periph_MassBin[NbinsInvMass] = { NULL };
        TH1F* Yield_Difference_MassBin[NbinsInvMass] = { NULL };
        TH1F* YieldWrtMass_Central[NbinsDeltaPhi]{ NULL };
        TH1F* YieldWrtMass_Periph[NbinsDeltaPhi]{ NULL };
        
    
    // Plots added for Pt Binning
//    TH1F* Yield_Central_MassBinPtBinned[NbinsInvMass][NbPtBins];//
//    memset( Yield_Central_MassBinPtBinned, 0, NbinsInvMass*NbPtBins*sizeof(TH1F*));
//    TH1F* Yield_Periph_MassBinPtBinned[NbinsInvMass][NbPtBins];//
//    memset( Yield_Periph_MassBinPtBinned, 0, NbinsInvMass*NbPtBins*sizeof(TH1F*));
    TH1F* Yield_Difference_MassBinPtBinned[NbinsInvMass][NbPtBins];//
    memset( Yield_Difference_MassBinPtBinned, 0, NbinsInvMass*NbPtBins*sizeof(TH1F*));
    TH1F* YieldWrtMass_CentralPtBinned[NbinsDeltaPhi][NbPtBins];//
    memset( YieldWrtMass_CentralPtBinned, 0, NbinsDeltaPhi*NbPtBins*sizeof(TH1F*));
    TH1F* YieldWrtMass_PeriphPtBinned[NbinsDeltaPhi][NbPtBins];//
    memset( YieldWrtMass_PeriphPtBinned, 0, NbinsDeltaPhi*NbPtBins*sizeof(TH1F*));
    TH1F* baselines0PtBinned[NbPtBins];
    memset( baselines0PtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F* coefficients0PtBinned[NbPtBins];
    memset( coefficients0PtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F* coefficients1PtBinned[NbPtBins];
    memset( coefficients1PtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F* coefficients2PtBinned[NbPtBins];
    memset( coefficients2PtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F* c2b0PtBinned[NbPtBins];
    memset( c2b0PtBinned, 0, NbPtBins*sizeof(TH1F*));
    //TH1F* V2JPsiTklPtBinned[NbPtBins]{ NULL }; CombineFits
    TH1F* V2JPsiTklPtBinned_noZYAM[NbPtBins];
    memset( V2JPsiTklPtBinned_noZYAM, 0, NbPtBins*sizeof(TH1F*));
    
//    TH1F* Yields_Central_1PtBinned[NbPtBins]{ NULL };
//       TH1F* Yields_Periph_1PtBinned[NbPtBins]{ NULL };
       TH1F* Yields_Difference_1PtBinned[NbPtBins];
    memset( Yields_Difference_1PtBinned, 0, NbPtBins*sizeof(TH1F*));
    
    TH1F *Yields_Central_1_CvetanPtBinned[NbPtBins];
    memset( Yields_Central_1_CvetanPtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F *Yields_Central_1_CvetanMePtBinned[NbPtBins];
    memset( Yields_Central_1_CvetanMePtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F *Yields_Central_1_ZYAMPtBinned[NbPtBins];
    memset( Yields_Central_1_ZYAMPtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F *Yields_Central_1_PRLTemplatePtBinned[NbPtBins];
    memset( Yields_Central_1_PRLTemplatePtBinned, 0, NbPtBins*sizeof(TH1F*));
    TH1F *Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[NbPtBins];
    memset( Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned, 0, NbPtBins*sizeof(TH1F*));

    
        // AE
        TH2F* hPtWrtMassInv[3]{NULL};

        
        TH1F* Yields_Difference_1(NULL);
        
        TH1F* coefficients0(NULL);
        TH1F* coefficients1(NULL);
        TH1F* coefficients2(NULL);
        TH1F* baselines0(NULL);
        TH1F* c2b0(NULL);
      //  TH1F* V2JPsiTkl(NULL);
    TH1F* V2JPsiTkl_noZYAM(NULL);
        
        Yields_Central_1 = new TH1F("Yields_Central_1",
                         "Yield of J/#psi-tkl in Central collisions wrt #Delta#phi",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Central_1->SetXTitle("#Delta#phi of the correlation (rad)");
        Yields_Central_1->SetYTitle("Yield_{Central}");
        
        Yields_Periph_1 = new TH1F("Yields_Periph_1",
                         "Yield of J/#psi-tkl in Periph collisions wrt #Delta#phi",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Periph_1->SetXTitle("#Delta#phi of the correlation (rad)");
        Yields_Periph_1->SetYTitle("Yield_{Periph}");
    
    Baseline_Central_1 = new TH1F("Baseline_Central_1",
                     "Baseline of J/#psi-tkl in Central collisions wrt #Delta#phi",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Baseline_Central_1->SetXTitle("#Delta#phi of the correlation (rad)");
    Baseline_Central_1->SetYTitle("Baseline_{Central}");
    
    Baseline_Periph_1 = new TH1F("Baseline_Periph_1",
                     "Baseline of J/#psi-tkl in Periph collisions wrt #Delta#phi",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Baseline_Periph_1->SetXTitle("#Delta#phi of the correlation (rad)");
    Baseline_Periph_1->SetYTitle("Baseline_{Periph}");
    
    Yields_Central_1_MinusBaseline = new TH1F("Yields_Central_1_MinusBaseline",
                     "Yield of J/#psi-tkl in Central collisions wrt #Delta#phi - Subtracted Baseline",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yields_Central_1_MinusBaseline->SetXTitle("#Delta#phi of the correlation (rad)");
    Yields_Central_1_MinusBaseline->SetYTitle("Yield_{Central} - Baseline_{Central}");
    
    Yields_Periph_1_MinusBaseline = new TH1F("Yields_Periph_1_MinusBaseline",
                     "Yield of J/#psi-tkl in Periph collisions wrt #Delta#phi - Subtracted Baseline",
                     NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
    Yields_Periph_1_MinusBaseline->SetXTitle("#Delta#phi of the correlation (rad)");
    Yields_Periph_1_MinusBaseline->SetYTitle("Yield_{Periph} - Baseline_{Periph}");
        
        Yields_Difference_1 = new TH1F("Yields_Difference_1",
                         "Subtracted yield of J/#psi-tkl (Central-Periph) wrt #Delta#phi",
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Difference_1->SetXTitle("#Delta#phi of the correlation (rad)");
        Yields_Difference_1->SetYTitle("Yield_{Subtracted}");
        
        coefficients0 = new TH1F("coefficients0",
                            "0^{th} Fourier coefficient wrt dimuon mass",
                            NbinsInvMass,MassBins);
           coefficients0->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients0->SetYTitle("0^{th} Fourier coefficient");
        coefficients1 = new TH1F("coefficients1",
                               "1^{st} Fourier coefficient wrt dimuon mass",
                               NbinsInvMass,MassBins);
              coefficients1->SetXTitle("Mass of dimuon (GeV/c^{2})");
              coefficients1->SetYTitle("1^{st} Fourier coefficient");
        coefficients2 = new TH1F("coefficients2",
                               "2^{nd} Fourier coefficient wrt dimuon mass",
                               NbinsInvMass,MassBins);
              coefficients2->SetXTitle("Mass of dimuon (GeV/c^{2})");
              coefficients2->SetYTitle("2^{nd} Fourier coefficient");
        baselines0 = new TH1F("baselines0",
                         "Baseline_{Periph} wrt dimuon mass",
                         NbinsInvMass,MassBins);
        baselines0->SetXTitle("Mass of dimuon (GeV/c^{2})");
        baselines0->SetYTitle("Baseline_{Periph}");
        c2b0 = new TH1F("c2b0",
                         "2^{nd} Fourier coefficient + Baseline_{Periph} wrt dimuon mass",
                         NbinsInvMass,MassBins);
        c2b0->SetXTitle("Mass of dimuon (GeV/c^{2})");
        c2b0->SetYTitle("2^{nd} Fourier coefficient + Baseline_{Periph}");
        V2JPsiTkl = new TH1F("V2JPsiTkl",
                               "V_{2,J/#psi-tkl} wrt dimuon mass",
                               NbinsInvMass,MassBins);
              V2JPsiTkl->SetXTitle("Mass of dimuon (GeV/c^{2})");
              V2JPsiTkl->SetYTitle("V_{2,J/#psi-tkl}");
    V2JPsiTkl_noZYAM = new TH1F("V2JPsiTkl_noZYAM",
                     "V_{2,J/#psi-tkl} wrt dimuon mass - no ZYAM",
                     NbinsInvMass,MassBins);
    V2JPsiTkl_noZYAM->SetXTitle("Mass of dimuon (GeV/c^{2})");
    V2JPsiTkl_noZYAM->SetYTitle("V_{2,J/#psi-tkl} - no ZYAM");
        
        char hname[200];
        char hname1[200];
        char hname2[200];
        
        for (int j=0; j <NbinsInvMass; j++){
           sprintf(hname1,"Projected yield in Mass Bin %f GeV to %f GeV - Central",MassBins[j],MassBins[j+1]);
           sprintf(hname2,"Projected yield in Mass Bin, #m_{#mu#mu} #in [%.2f,%.2f] GeV/c^{2} - Central",MassBins[j],MassBins[j+1]);
            Yield_Central_MassBin[j] = new TH1F(hname1, hname2,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           }
           
           for (int j=0; j <NbinsInvMass; j++){
        sprintf(hname1,"Projected yield in Mass Bin %f GeV to %f GeV - Periph",MassBins[j],MassBins[j+1]);
        sprintf(hname2,"Projected yield in Mass Bin, #m_{#mu#mu} #in [%.2f,%.2f] GeV/c^{2} - Periph",MassBins[j],MassBins[j+1]);
           Yield_Periph_MassBin[j] = new TH1F(hname1, hname2,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           }
           
           for (int j=0; j <NbinsInvMass; j++){
               sprintf(hname1,"Projected yield in Mass Bin %f GeV to %f GeV - Difference",MassBins[j],MassBins[j+1]);
               sprintf(hname2,"Projected yield in Mass Bin, #m_{#mu#mu} #in [%.2f,%.2f] GeV/c^{2} - Difference",MassBins[j],MassBins[j+1]);
           Yield_Difference_MassBin[j] = new TH1F(hname1, hname2,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           }
        
        TH1I* InvMass_Central(NULL);
        TH1I* InvMass_Periph(NULL);
        
            
        for(int p=0; p<NbinsDeltaPhi; p++){
                char hname[200];
                char hname2[200];
            
            sprintf(hname,"Yields PhiBin %d Periph",p);
            sprintf(hname2,"Yields dimuon-tkl wrt wrt dimuon mass, Periph, #Delta#phi #in [#frac{%d#pi}{6},#frac{%d#pi}{6}]",p-3, p-2);
            YieldWrtMass_Periph[p] = new TH1F(hname,
                              hname2,
                              NbinsInvMass,MassBins);
            YieldWrtMass_Periph[p]->SetXTitle("Correlation dimuon Inv Mass (GeV/c^{2})");
            
            sprintf(hname,"Yields PhiBin %d Central",p);
            sprintf(hname2,"Yields dimuon-tkl wrt wrt dimuon mass, Central, #Delta#phi #in [#frac{%d#pi}{6},#frac{%d#pi}{6}]",p-3, p-2);
            YieldWrtMass_Central[p] = new TH1F(hname,
                              hname2,
                              NbinsInvMass,MassBins);
            YieldWrtMass_Central[p]->SetXTitle("Correlation dimuon Inv Mass (GeV/c^{2})");
        }
    

        
        
    char haxis[200];
    sprintf(haxis, "Count within bin of %d MeV/c^{2}", 1000*(MaxInvMass-MinInvMass)/NbinsDimuInvMass);
        hnseg = new TH1I("hnseg",
                         "Invariant mass of dimuon",
                         NbinsDimuInvMass,MinInvMass,MaxInvMass);
        hnseg->SetXTitle("Mass of dimuon (GeV/c^{2})");
        hnseg->SetYTitle(haxis);
    
        hpool = new TH1F("hpull",
                            "Invariant mass of dimuon - Pulls",
                            NbinsDimuInvMass,MinInvMass,MaxInvMass);
           hpool->SetXTitle("Mass of dimuon (GeV/c^{2})");
           hpool->SetYTitle("Pull (Sigma units)");
    hpool->GetYaxis()->SetRangeUser(-10., 10.);
    hpool->SetStats(kFALSE);
    
    TLine *l0=new TLine(minmass,0,maxmass,0);
    TLine *l1=new TLine(minmass,1,maxmass,1);
    TLine *lm1=new TLine(minmass,-1,maxmass,-1);
    l1->SetLineColor(8);
    l1->SetLineStyle(10);
    lm1->SetLineColor(8);
    lm1->SetLineStyle(10);
    TLine *l3=new TLine(minmass,3,maxmass,3);
    TLine *lm3=new TLine(minmass,-3,maxmass,-3);
    l3->SetLineColor(46);
    l3->SetLineStyle(10);
    lm3->SetLineColor(46);
    lm3->SetLineStyle(10);

    
    hpool->GetListOfFunctions()->Add(l0);
    hpool->GetListOfFunctions()->Add(l1);
    hpool->GetListOfFunctions()->Add(lm1);
    hpool->GetListOfFunctions()->Add(l3);
    hpool->GetListOfFunctions()->Add(lm3);

        InvMass_Central = new TH1I("InvMass_Central",
                         "Invariant mass of dimuon - Central",
                         NbinsDimuInvMass,MinInvMass,MaxInvMass);
        InvMass_Central->SetXTitle("Mass of dimuon (GeV/c^{2})");
        InvMass_Central->SetYTitle(haxis);
        InvMass_Periph = new TH1I("InvMass_Periph",
                         "Invariant mass of dimuon - Periph",
                         NbinsDimuInvMass,MinInvMass,MaxInvMass);
        InvMass_Periph->SetXTitle("Mass of dimuon (GeV/c^{2})");
    InvMass_Periph->SetYTitle(haxis);

        hPtWrtMassInv[0] = new TH2F("hPtWrtMassInv_0",
                          "p_{T} wrt dimuon inv mass - All centralities",
                          NbinsDimuInvMass,MinInvMass,MaxInvMass,NbPtBins,PtBins);
        hPtWrtMassInv[0]->SetXTitle("Dimuon inv mass (GeV/c^{2})");
        hPtWrtMassInv[0]->SetYTitle("p_{T} (GeV/c)");
        hPtWrtMassInv[1] = new TH2F("hPtWrtMassInv_1",
                          "p_{T} wrt dimuon inv mass - Central",
                          NbinsDimuInvMass,MinInvMass,MaxInvMass,NbPtBins,PtBins);
        hPtWrtMassInv[1]->SetXTitle("Dimuon inv mass (GeV/c^{2})");
        hPtWrtMassInv[1]->SetYTitle("p_{T} (GeV/c)");
        hPtWrtMassInv[2] = new TH2F("hPtWrtMassInv_2",
                          "p_{T} wrt dimuon inv mass - Periph",
                          NbinsDimuInvMass,MinInvMass,MaxInvMass,NbPtBins,PtBins);
        hPtWrtMassInv[2]->SetXTitle("Dimuon inv mass (GeV/c^{2})");
        hPtWrtMassInv[2]->SetYTitle("p_{T} (GeV/c)");

    
    
    
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for (int j=0; j <NbinsInvMass; j++){
        char hname1[200];
        char hname2[200];
        sprintf(hname1,"Projected yield in Mass Bin %f GeV to %f GeV - Central - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
        sprintf(hname2,"Projected yield in Mass Bin, #m_{#mu#mu} #in [%.2f,%.2f] GeV/c^{2} - Central - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",MassBins[j],MassBins[j+1],PtBins[ptbin], PtBins[ptbin+1]);
        Yield_Central_MassBinPtBinned[j][ptbin] = new TH1F(hname1, hname2,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        }
        
        for (int j=0; j <NbinsInvMass; j++){
        char hname1[200];
        char hname2[200];
        sprintf(hname1,"Projected yield in Mass Bin %f GeV to %f GeV - Periph - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
        sprintf(hname2,"Projected yield in Mass Bin, #m_{#mu#mu} #in [%.2f,%.2f] GeV/c^{2} - Periph - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",MassBins[j],MassBins[j+1],PtBins[ptbin], PtBins[ptbin+1]);
        Yield_Periph_MassBinPtBinned[j][ptbin] = new TH1F(hname1, hname2,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        }
        
        for (int j=0; j <NbinsInvMass; j++){
        char hname1[200];
        char hname2[200];
        sprintf(hname1,"Projected yield in Mass Bin %f GeV to %f GeV - Difference - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
        sprintf(hname2,"Projected yield in Mass Bin, #m_{#mu#mu} #in [%.2f,%.2f] GeV/c^{2} - Difference - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",MassBins[j],MassBins[j+1],PtBins[ptbin], PtBins[ptbin+1]);
        Yield_Difference_MassBinPtBinned[j][ptbin] = new TH1F(hname1, hname2,NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        }
        
        for(int p=0; p<NbinsDeltaPhi; p++){
                   char hname[200];
                   char hname2[200];
               
               sprintf(hname,"Yields PhiBinPtBinned %d Periph %d",p,ptbin);
               sprintf(hname2,"Yields CorrelationPtBinned wrt dimuon mass, Periph, #Delta#phi #in [#frac{%d#pi}{6},#frac{%d#pi}{6}]",p-3, p-2);
               YieldWrtMass_PeriphPtBinned[p][ptbin] = new TH1F(hname,
                                 hname2,
                                 NbinsInvMass,MassBins);
               YieldWrtMass_PeriphPtBinned[p][ptbin]->SetXTitle("Correlation dimuon Inv Mass (GeV/c^{2})");
               
               sprintf(hname,"Yields PhiBinPtBinned %d Central %d",p,ptbin);
               sprintf(hname2,"Yields CorrelationPtBinned wrt dimuon mass, Central, #Delta#phi #in [#frac{%d#pi}{6},#frac{%d#pi}{6}]",p-3, p-2);
               YieldWrtMass_CentralPtBinned[p][ptbin] = new TH1F(hname,
                                 hname2,
                                 NbinsInvMass,MassBins);
               YieldWrtMass_CentralPtBinned[p][ptbin]->SetXTitle("Correlation dimuon Inv Mass (GeV/c^{2})");
           }
        
        char hname[200];
        char hname2[200];
        sprintf(hname,"Yields_Central_1PtBinned %d",ptbin);
        sprintf(hname2,"Yield of J/#psi-tkl in Central collisions wrt #Delta#phi - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        Yields_Central_1PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Central_1PtBinned[ptbin]->SetXTitle("#Delta#phi (rad)");
        Yields_Central_1PtBinned[ptbin]->SetYTitle("Yield_{Central}");
        
        sprintf(hname,"Yields_Periph_1PtBinned %d",ptbin);
        sprintf(hname2,"Yield of J/#psi-tkl in Periph collisions wrt #Delta#phi - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        Yields_Periph_1PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Periph_1PtBinned[ptbin]->SetXTitle("#Delta#phi (rad)");
        Yields_Periph_1PtBinned[ptbin]->SetYTitle("Yield_{Periph}");
        
        sprintf(hname,"Yields_Difference_1PtBinned %d",ptbin);
        sprintf(hname2,"Subtracted Yield of J/#psi-tkl wrt #Delta#phi - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        Yields_Difference_1PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
        Yields_Difference_1PtBinned[ptbin]->SetXTitle("#Delta#phi (rad)");
        Yields_Difference_1PtBinned[ptbin]->SetYTitle("Yield_{Subtracted}");
        
        sprintf(hname,"Baseline_Central_1PtBinned %d",ptbin);
               sprintf(hname2,"Baseline of J/#psi-tkl in Central collisions wrt #Delta#phi - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        Baseline_Central_1PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Baseline_Central_1PtBinned[ptbin]->SetXTitle("#Delta#phi (rad)");
           Baseline_Central_1PtBinned[ptbin]->SetYTitle("Baseline_{Central}");
           
           sprintf(hname,"Baseline_Periph_1PtBinned %d",ptbin);
                  sprintf(hname2,"Baseline of J/#psi-tkl in Periph collisions wrt #Delta#phi - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
           Baseline_Periph_1PtBinned[ptbin] = new TH1F(hname,
                               hname2,
                               NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
              Baseline_Periph_1PtBinned[ptbin]->SetXTitle("#Delta#phi (rad)");
              Baseline_Periph_1PtBinned[ptbin]->SetYTitle("Baseline_{Periph}");
        
        sprintf(hname,"Yields_Central_1_MinusBaselinePtBinned %d",ptbin);
        sprintf(hname2,"Yield of J/#psi-tkl in Central collisions wrt #Delta#phi - Subtracted Baseline - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        Yields_Central_1_MinusBaselinePtBinned[ptbin] = new TH1F(hname,hname2,
                            NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Yields_Central_1_MinusBaselinePtBinned[ptbin]->SetXTitle("#Delta#phi (rad)");
           Yields_Central_1_MinusBaselinePtBinned[ptbin]->SetYTitle("Yield_{Central} - Baseline_{Central}");
        
        sprintf(hname,"Yields_Periph_1_MinusBaselinePtBinned %d",ptbin);
        sprintf(hname2,"Yield of J/#psi-tkl in Periph collisions wrt #Delta#phi - Subtracted Baseline - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        Yields_Periph_1_MinusBaselinePtBinned[ptbin] = new TH1F(hname,hname2,
                            NbinsDeltaPhi,MinDeltaPhi,MaxDeltaPhi);
           Yields_Periph_1_MinusBaselinePtBinned[ptbin]->SetXTitle("#Delta#phi (rad)");
           Yields_Periph_1_MinusBaselinePtBinned[ptbin]->SetYTitle("Yield_{Periph} - Baseline_{Periph}");
        
        sprintf(hname,"coefficients0PtBinned %d",ptbin);
        sprintf(hname2,"0^{th} Fourier coefficient - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        coefficients0PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsInvMass,MassBins);
           coefficients0PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients0PtBinned[ptbin]->SetYTitle("0^{th} Fourier coefficient");
        
        sprintf(hname,"coefficients1PtBinned %d",ptbin);
        sprintf(hname2,"1^{st} Fourier coefficient - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        coefficients1PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsInvMass,MassBins);
           coefficients1PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients1PtBinned[ptbin]->SetYTitle("1^{st} Fourier coefficient");
        
        sprintf(hname,"coefficients2PtBinned %d",ptbin);
        sprintf(hname2,"2^{nd} Fourier coefficient - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        coefficients2PtBinned[ptbin] = new TH1F(hname,
                            hname2,
                            NbinsInvMass,MassBins);
           coefficients2PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
           coefficients2PtBinned[ptbin]->SetYTitle("2^{nd} Fourier coefficient");
        
        sprintf(hname,"baselines0PtBinned %d",ptbin);
        sprintf(hname2,"Baseline_{Periph} - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        baselines0PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsInvMass,MassBins);
        baselines0PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
        baselines0PtBinned[ptbin]->SetYTitle("Baseline_{Periph}");
        
        sprintf(hname,"c2b0PtBinned %d",ptbin);
        sprintf(hname2,"2^{nd} Fourier coefficient + Baseline_{Periph} - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        c2b0PtBinned[ptbin] = new TH1F(hname,
                         hname2,
                         NbinsInvMass,MassBins);
        c2b0PtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
        c2b0PtBinned[ptbin]->SetYTitle("2^{nd} Fourier coefficient + Baseline_{Periph}");
        
        sprintf(hname,"V2JPsiTklPtBinned %d",ptbin);
        sprintf(hname2,"V_{2,J/#psi-tkl} wrt dimuon mass - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        V2JPsiTklPtBinned[ptbin] = new TH1F(hname,
                               hname2,
                               NbinsInvMass,MassBins);
              V2JPsiTklPtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
              V2JPsiTklPtBinned[ptbin]->SetYTitle("V_{2,J/#psi-tkl}");
        
        sprintf(hname,"V2JPsiTklPtBinned no ZYAM %d",ptbin);
        sprintf(hname2,"V_{2,J/#psi-tkl} wrt dimuon mass - PtBinned - no ZYAM, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        V2JPsiTklPtBinned_noZYAM[ptbin] = new TH1F(hname,
                               hname2,
                               NbinsInvMass,MassBins);
              V2JPsiTklPtBinned_noZYAM[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
              V2JPsiTklPtBinned_noZYAM[ptbin]->SetYTitle("V_{2,J/#psi-tkl} - no ZYAM");
        
        sprintf(hname,"V22ATLASPtBinned %d",ptbin);
        sprintf(hname2,"V_{2,J/#psi-tkl} wrt dimuon mass ATLAS 2 - PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        V22ATLASPtBinned[ptbin] = new TH1F(hname,
                               hname2,
                               NbinsInvMass,MassBins);
              V22ATLASPtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
              V22ATLASPtBinned[ptbin]->SetYTitle("V_{2,J/#psi-tkl}");
        
        sprintf(hname,"FATLASPtBinned %d",ptbin);
        sprintf(hname2,"F wrt dimuon mass ATLAS 2- PtBinned, p_{T} #in [%.2f,%.2f] GeV/c",PtBins[ptbin], PtBins[ptbin+1]);
        FATLASPtBinned[ptbin] = new TH1F(hname,
                               hname2,
                               NbinsInvMass,MassBins);
              FATLASPtBinned[ptbin]->SetXTitle("Mass of dimuon (GeV/c^{2})");
              FATLASPtBinned[ptbin]->SetYTitle("V_{2,J/#psi-tkl}");
    }
    
    
    
    
    
    // *************************
    // Analyse                 *
    // *************************
    
    
    TCanvas *cinvmass = new TCanvas("cinvmass","Fitting All Phi - Different Centralities",10,10,700,500);
    cinvmass->SetTitle("Inv Mass Fits");
    // ROOT THAT Inv mass fits
    cinvmass->Divide(1,3);
    
    // Récup histos
    if(!isCentralityStudy){
        
    TFile *filerec = new TFile(FitFileName);
    hnseg = (TH1I*)filerec->Get("hnseg");
    InvMass_Central = (TH1I*)filerec->Get("InvMass_Central");
    InvMass_Periph = (TH1I*)filerec->Get("InvMass_Periph");
    V2JPsiTkl = (TH1F*)filerec->Get("V2JPsiTkl");
    coefficients0 = (TH1F*)filerec->Get("coefficients0");
    coefficients1 = (TH1F*)filerec->Get("coefficients1");
    coefficients2 = (TH1F*)filerec->Get("coefficients2");
    baselines0 = (TH1F*)filerec->Get("baselines0");
    
    for(int index=0; index<NbinsInvMass; index++){
        cout << "Content "<< baselines0->GetBinContent(index+1)<<endl;
        cout << "Error "<< baselines0->GetBinError(index+1)<<endl;
    }
    
    
    for(int p=0; p<NbinsDeltaPhi; p++){
        sprintf(hname,"Yields PhiBin %d Central",p);
        YieldWrtMass_Central[p] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"Yields PhiBin %d Periph",p);
        YieldWrtMass_Periph[p] = (TH1F*)filerec->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central",MassBins[j],MassBins[j+1]);
    Yield_Central_MassBin[j] = (TH1F*)filerec->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph",MassBins[j],MassBins[j+1]);
    Yield_Periph_MassBin[j] = (TH1F*)filerec->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference",MassBins[j],MassBins[j+1]);
    Yield_Difference_MassBin[j] = (TH1F*)filerec->Get(hname);
    }
    
    
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for(int j=0; j<NbinsInvMass; j++){
           sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
            Yield_Central_MassBinPtBinned[j][ptbin] = (TH1F*)filerec->Get(hname);
            
             sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
                   Yield_Periph_MassBinPtBinned[j][ptbin] = (TH1F*)filerec->Get(hname);
           
            sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
            Yield_Difference_MassBinPtBinned[j][ptbin] = (TH1F*)filerec->Get(hname);

        }
    }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for(int p=0; p<NbinsDeltaPhi; p++){
             sprintf(hname,"Yields PhiBinPtBinned %d Periph %d",p,ptbin);
                          YieldWrtMass_PeriphPtBinned[p][ptbin] = (TH1F*)filerec->Get(hname);
                          
                          sprintf(hname,"Yields PhiBinPtBinned %d Central %d",p,ptbin);
                          YieldWrtMass_CentralPtBinned[p][ptbin] = (TH1F*)filerec->Get(hname);

        }
    }
    
    for(int j=0;j<3;j++){
        for (int i=0;i<NbPtBins;i++){
            sprintf(hname,"bin%d_%d",i+1, j);
            hPtWrtMassInvSliced[j][i] = (TH1F*)filerec->Get(hname);
        }
    }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        sprintf(hname,"coefficients0PtBinned %d",ptbin);
        coefficients0PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"coefficients1PtBinned %d",ptbin);
        coefficients1PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"coefficients2PtBinned %d",ptbin);
        coefficients2PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"baselines0PtBinned %d",ptbin);
        baselines0PtBinned[ptbin] = (TH1F*)filerec->Get(hname);
        sprintf(hname,"V2JPsiTklPtBinned %d",ptbin);
        V2JPsiTklPtBinned[ptbin] = (TH1F*)filerec->Get(hname);
    }
        
    }
    
    if(isCentralityStudy){
        
    TFile *filereccent = new TFile(FitCentralFileName);
    TFile *filerecperiph = new TFile(FitPeriphFileName);
    hnseg = (TH1I*)filereccent->Get("hnseg");
    InvMass_Central = (TH1I*)filereccent->Get("InvMass_Central");
    InvMass_Periph = (TH1I*)filerecperiph->Get("InvMass_Periph");
    V2JPsiTkl = (TH1F*)filereccent->Get("V2JPsiTkl");
    coefficients0 = (TH1F*)filereccent->Get("coefficients0");
    coefficients1 = (TH1F*)filereccent->Get("coefficients1");
    coefficients2 = (TH1F*)filereccent->Get("coefficients2");
    baselines0 = (TH1F*)filerecperiph->Get("baselines0");
    
    for(int index=0; index<NbinsInvMass; index++){
        cout << "Content "<< baselines0->GetBinContent(index+1)<<endl;
        cout << "Error "<< baselines0->GetBinError(index+1)<<endl;
    }
    
    
    for(int p=0; p<NbinsDeltaPhi; p++){
        sprintf(hname,"Yields PhiBin %d Central",p);
        YieldWrtMass_Central[p] = (TH1F*)filereccent->Get(hname);
        sprintf(hname,"Yields PhiBin %d Periph",p);
        YieldWrtMass_Periph[p] = (TH1F*)filerecperiph->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central",MassBins[j],MassBins[j+1]);
    Yield_Central_MassBin[j] = (TH1F*)filereccent->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph",MassBins[j],MassBins[j+1]);
    Yield_Periph_MassBin[j] = (TH1F*)filerecperiph->Get(hname);
    }
    
    for (int j=0; j <NbinsInvMass; j++){
    sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference",MassBins[j],MassBins[j+1]);
    Yield_Difference_MassBin[j] = (TH1F*)filereccent->Get(hname);
    }
    
    
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for(int j=0; j<NbinsInvMass; j++){
           sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Central - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
            Yield_Central_MassBinPtBinned[j][ptbin] = (TH1F*)filereccent->Get(hname);
            
             sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Periph - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
                   Yield_Periph_MassBinPtBinned[j][ptbin] = (TH1F*)filerecperiph->Get(hname);
           
            sprintf(hname,"Projected yield in Mass Bin %f GeV to %f GeV - Difference - PtBinned %d",MassBins[j],MassBins[j+1],ptbin);
            Yield_Difference_MassBinPtBinned[j][ptbin] = (TH1F*)filereccent->Get(hname);

        }
    }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        for(int p=0; p<NbinsDeltaPhi; p++){
             sprintf(hname,"Yields PhiBinPtBinned %d Periph %d",p,ptbin);
                          YieldWrtMass_PeriphPtBinned[p][ptbin] = (TH1F*)filerecperiph->Get(hname);
                          
                          sprintf(hname,"Yields PhiBinPtBinned %d Central %d",p,ptbin);
                          YieldWrtMass_CentralPtBinned[p][ptbin] = (TH1F*)filereccent->Get(hname);

        }
    }
    
    for(int j=0;j<3;j+=2){
        for (int i=0;i<NbPtBins;i++){
            sprintf(hname,"bin%d_%d",i+1, j);
            hPtWrtMassInvSliced[j][i] = (TH1F*)filereccent->Get(hname);//Bigclasses
        }
    }
        for (int i=0;i<NbPtBins;i++){
            sprintf(hname,"bin%d_%d",i+1, 1);
            hPtWrtMassInvSliced[1][i] = (TH1F*)filerecperiph->Get(hname);//Bigclasses
        }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        sprintf(hname,"coefficients0PtBinned %d",ptbin);
        coefficients0PtBinned[ptbin] = (TH1F*)filereccent->Get(hname);
        sprintf(hname,"coefficients1PtBinned %d",ptbin);
        coefficients1PtBinned[ptbin] = (TH1F*)filereccent->Get(hname);
        sprintf(hname,"coefficients2PtBinned %d",ptbin);
        coefficients2PtBinned[ptbin] = (TH1F*)filereccent->Get(hname);
        sprintf(hname,"baselines0PtBinned %d",ptbin);
        baselines0PtBinned[ptbin] = (TH1F*)filerecperiph->Get(hname);
        sprintf(hname,"V2JPsiTklPtBinned %d",ptbin);
        V2JPsiTklPtBinned[ptbin] = (TH1F*)filereccent->Get(hname);
    }
        
    }
    
    TCanvas*c16PtBinned=new TCanvas();
    // ROOT THAT Extraction 2 plot
    c16PtBinned->SetTitle("Extraction method 2 PtBinned");
    c16PtBinned->Divide(2,3);
    //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        c16PtBinned->cd(ptbin+1);
        c2b0PtBinned[ptbin]->Add(coefficients0PtBinned[ptbin], baselines0PtBinned[ptbin]);
        V2JPsiTklPtBinned[ptbin]->Divide(coefficients2PtBinned[ptbin], c2b0PtBinned[ptbin]);
        
        V2JPsiTklPtBinned[ptbin]->Draw();
    }
    
    TCanvas*c16PtBinned_noZYAM=new TCanvas();
    // ROOT THAT Extraction 2 plot
    c16PtBinned_noZYAM->SetTitle("Extraction method 2 PtBinned - no ZYAM");
    c16PtBinned_noZYAM->Divide(2,3);
    //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        c16PtBinned_noZYAM->cd(ptbin+1);
        V2JPsiTklPtBinned_noZYAM[ptbin]->Divide(coefficients2PtBinned[ptbin], coefficients0PtBinned[ptbin]);
        
        V2JPsiTklPtBinned_noZYAM[ptbin]->Draw();
    }
    
    
    
    cout << "START ALL"<<endl;
    
    
    
    
    
    
    
    
    
    
    
    TFitResultPtr res;
    TFitResultPtr rescent;
    TFitResultPtr resperiph;
    ROOT::Fit::FitResult resu;
    
   ROOT::Math::Minimizer* minim{NULL};
    
        char histoname[50];
        //Values OK for 20M evts
    //    double startvalues[16] = {800, 3.096916, 0.07, 0.9, 10, 2, 15, 15,   10,     0.1,          50000,   0.01,        1,         1,   1,   1};
    //    double lowvalues[16] =   {1,    3.05,    0.03, 0.7, 1,  1, 5,  0.01, 0.0001,  0.000000001,  10,    0.0000000001, 0.000001, -10, -10, -10};
    //    double highvalues[16] = {100000, 3.15,   0.1,  1.1, 30, 5, 25, 50,   100000, 10,           100000, 10,     10,        10,  10,  10};
    
    
        //Fixme handle startvalues not needed now, just no combine fits
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
        res = FittingAllInvMassBin(histoname, cinvmass, 0, SignalF, BackgroundF, minmass, maxmass, ratsigma);
        double par[numParameters];
        
        // Parameter setting (either from mass fit either start values)
        for(int i=0; i <numParameters+numParametersBkgV2+1; i++){
                par[i] = res->Parameter(i);
                if(i>=numParameters){
                    par[i] = 1;
                }
        }
    
    //FIXME: Add the pulls in mass inv fit FIXED
    
    {
        ofstream myfiletxt;
        myfiletxt.open("/tmp/pullstot.txt");
        double params[numParameters];
        for(int i=0; i<numParameters; i++){
            params[i] = par[i];
        }
        TF1 *fitFcn = new TF1("fitFcn",MassFunction,MinInvMassFit,MaxInvMassFit,numParameters);
          fitFcn->SetNpx(500);
          fitFcn->SetLineWidth(4);
          fitFcn->SetLineColor(kMagenta);
          fitFcn->SetParameters(params);

        //Pour chaue bin
            //Trouver la valeur du centre du bin pour le fit et pour les données
            //Calculer le pool
            // Ajouter le point à un histo
        
        myfiletxt << "POUET" <<endl;
         myfiletxt << "POUET " << hnseg->GetNbinsX() <<endl;
        for(int bin=0; bin<hnseg->GetNbinsX(); bin++){
            if(MinInvMass + (MaxInvMass-MinInvMass)*(0.5+bin)/NbinsDimuInvMass < MinInvMassFit || MinInvMass + (MaxInvMass-MinInvMass)*(0.5+bin)/NbinsDimuInvMass > MaxInvMassFit){
                continue;
            }
            double fit_prediction = fitFcn->Eval(MinInvMass + (MaxInvMass-MinInvMass)*(0.5+bin)/NbinsDimuInvMass);
            double data = hnseg->GetBinContent(bin+1);
            double error = hnseg->GetBinError(bin+1);
            if(error>0){
                double pool = (data-fit_prediction)/error;

                hpool->SetBinContent(bin+1,pool);
                
//                if(bin < 10){
//                myfiletxt << "Bin " << bin << " : Fit - " << fit_prediction << " Data - " << data << " Error - " << error << " Pool value - " << pool <<endl;
//                }
               }
        }
        myfiletxt.close();

        TCanvas*cpoolall = new TCanvas();
        cpoolall->Divide(1,1);
        cpoolall->SetTitle("Deviations inv mass fit, all p_{T}");
        hpool->SetTitle("Pulls from invariant mass fit, all p_{T}");
        hpool->DrawCopy();
        
        TCanvas*cmassandpoolall = new TCanvas();
        TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
        pad1->SetBottomMargin(0.00001);
        pad1->SetBorderMode(0);
        pad2->SetTopMargin(0.00001);
        pad2->SetBottomMargin(0.1);
        pad2->SetBorderMode(0);
        pad1->Draw();
        pad2->Draw();
        cmassandpoolall->SetTitle("Invariant mass fit and pulls, all p_{T}");
        pad1->cd();
        TVirtualPad *padvirt = cinvmass->GetPad(1);
        double x1 = padvirt->GetAbsXlowNDC();
        double x2 = x1 + padvirt->GetAbsWNDC();
        double y1 = padvirt->GetAbsYlowNDC();
        double y2 = y1 + padvirt->GetAbsHNDC();
        padvirt->SetPad(0,0.33,1,1);
        padvirt->DrawClone();
        padvirt->SetPad(x1,y1,x2,y2);
        pad2->cd();
        hpool->DrawCopy();
        
        sprintf(CanvasName,"%s/Pulls_InvMassFit.pdf",CanvasFolderName);
        cpoolall->SaveAs(CanvasName);
        
    }
    
{
    
            TCanvas*c16=new TCanvas();
           // ROOT THAT Extraction 2 plot
           c16->SetTitle("Extraction method 2");
           //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
           c2b0->Add(coefficients0, baselines0);
           V2JPsiTkl->Divide(coefficients2, c2b0);
           
           V2JPsiTkl->Draw();
        
        //Definition and design of function on V2 plot
        c16->cd();
        TF1 *fitV2_2 = new TF1("fitV2_2",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        fitV2_2->SetNpx(500);
        fitV2_2->SetLineWidth(4);
        fitV2_2->SetLineColor(kMagenta);
        fitV2_2->SetParameters(par);
        TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        backFcnV2_2->SetLineColor(kRed);
        TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        signalFcnJPsiV2_2->SetLineColor(kBlue);
        signalFcnJPsiV2_2->SetNpx(500);
        
           TVirtualFitter::Fitter(V2JPsiTkl)->SetMaxIterations(10000);
           TVirtualFitter::Fitter(V2JPsiTkl)->SetPrecision();
            for(int i=0; i<numParameters; i++){
                fitV2_2->FixParameter(i,par[i]);
            }
    
    for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
        string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "V2_2");
        const char *titular = title.c_str();
        fitV2_2->SetParName(index,titular);
    }
           
//           fitV2_2->SetParName(0,"Norm_{J/#psi}");
//           fitV2_2->SetParName(1,"M_{J/#psi}");
//           fitV2_2->SetParName(2,"Sigma_{J/#psi}");
//           fitV2_2->SetParName(3,"a_{1}");
//           fitV2_2->SetParName(4,"n_{1}");
//           fitV2_2->SetParName(5,"a_{2}");
//           fitV2_2->SetParName(6,"n_{2}");
//           fitV2_2->SetParName(7,"Norm_{#Psi(2S)}");
//           fitV2_2->SetParName(8,"Norm_{TailLowM}");
//           fitV2_2->SetParName(9,"Exp_{TailLowM}");
//           fitV2_2->SetParName(10,"Norm_{TailHighM}");
//           fitV2_2->SetParName(11,"Exp_{TailHighM}");
//            fitV2_2->SetParName(12,"V2_2 J/#psi-tkl");
//            fitV2_2->SetParName(13,"V2_2 Bkg M2");
//        fitV2_2->SetParName(14,"V2_2 Bkg M1");
//        fitV2_2->SetParName(15,"V2_2 Bkg M0");
        
           double minParams[numParameters+numParametersBkgV2+1];
           double parErrors[numParameters+numParametersBkgV2+1];
        
        //Fit of V2
    gStyle->SetOptFit(1011);
        gStyle->SetOptStat("n");
            res = V2JPsiTkl->Fit("fitV2_2","SBMERI+","ep");
            TPaveStats *st = (TPaveStats*)V2JPsiTkl->FindObject("stats");
            st->SetX1NDC(0.75); //new x start position
            st->SetY1NDC(0.75); //new x end position
            Double_t para[numParameters+numParametersBkgV2+1];
            fitV2_2->GetParameters(para);
            backFcnV2_2->SetParameters(para);
            signalFcnJPsiV2_2->SetParameters(para);
        
        
        fitV2_2->Draw("same");
      //  signalFcnJPsiV2_2->Draw("same");
        backFcnV2_2->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
        legend->SetFillColorAlpha(kWhite, 0.);
        legend->SetBorderSize(0);
          legend->SetTextFont(42);
          legend->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
         legend->AddEntry(fitV2_2,message);

            if(res->CovMatrixStatus() == 3){
                  // sprintf(message,"The fit is a success");
            }
            else{
                   sprintf(message,"The fit is a failure");
                legend->AddEntry(fitV2_2,message);
            }
       // legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
        legend->AddEntry(backFcnV2_2,"Background");
          legend->AddEntry(V2JPsiTkl,"Data","lpe");
          legend->Draw();
    
    sprintf(CanvasName,"%s/V2_Ext2.pdf",CanvasFolderName);
    c16->SaveAs(CanvasName);
    
}
    
    
    { //Ext 2 no ZYAM
        
                TCanvas*c16_noZYAM=new TCanvas();
               // ROOT THAT Extraction 2 plot
               c16_noZYAM->SetTitle("Extraction method 2 - no ZYAM");
               //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 2)
               V2JPsiTkl_noZYAM->Divide(coefficients2, coefficients0);
               
               V2JPsiTkl_noZYAM->Draw();
            
            //Definition and design of function on V2 plot
            c16_noZYAM->cd();
            TF1 *fitV2_2 = new TF1("fitV2_2",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
            fitV2_2->SetNpx(500);
            fitV2_2->SetLineWidth(4);
            fitV2_2->SetLineColor(kMagenta);
            fitV2_2->SetParameters(par);
            TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
            backFcnV2_2->SetLineColor(kRed);
            TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
            signalFcnJPsiV2_2->SetLineColor(kBlue);
            signalFcnJPsiV2_2->SetNpx(500);
            
               TVirtualFitter::Fitter(V2JPsiTkl_noZYAM)->SetMaxIterations(10000);
               TVirtualFitter::Fitter(V2JPsiTkl_noZYAM)->SetPrecision();
                for(int i=0; i<numParameters; i++){
                    fitV2_2->FixParameter(i,par[i]);
                }
        
        for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
            string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "V2_2");
            const char *titular = title.c_str();
            fitV2_2->SetParName(index,titular);
        }
        
             
//               fitV2_2->SetParName(0,"Norm_{J/#psi}");
//               fitV2_2->SetParName(1,"M_{J/#psi}");
//               fitV2_2->SetParName(2,"Sigma_{J/#psi}");
//               fitV2_2->SetParName(3,"a_{1}");
//               fitV2_2->SetParName(4,"n_{1}");
//               fitV2_2->SetParName(5,"a_{2}");
//               fitV2_2->SetParName(6,"n_{2}");
//               fitV2_2->SetParName(7,"Norm_{#Psi(2S)}");
//               fitV2_2->SetParName(8,"Norm_{TailLowM}");
//               fitV2_2->SetParName(9,"Exp_{TailLowM}");
//               fitV2_2->SetParName(10,"Norm_{TailHighM}");
//               fitV2_2->SetParName(11,"Exp_{TailHighM}");
//                fitV2_2->SetParName(12,"V2_2 J/#psi");
//                fitV2_2->SetParName(13,"V2_2 Bkg M2");
//            fitV2_2->SetParName(14,"V2_2 Bkg M1");
//            fitV2_2->SetParName(15,"V2_2 Nkg M0");
            
               double minParams[numParameters+numParametersBkgV2+1];
               double parErrors[numParameters+numParametersBkgV2+1];
            
            //Fit of V2
                res = V2JPsiTkl_noZYAM->Fit("fitV2_2","SBMERI+","ep");
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)V2JPsiTkl_noZYAM->FindObject("stats");
                st->SetX1NDC(0.75); //new x start position
                st->SetY1NDC(0.75); //new x end position
                Double_t para[numParameters+numParametersBkgV2+1];
                fitV2_2->GetParameters(para);
                backFcnV2_2->SetParameters(para);
                signalFcnJPsiV2_2->SetParameters(para);
            
            
            fitV2_2->Draw("same");
          //  signalFcnJPsiV2_2->Draw("same");
            backFcnV2_2->Draw("same");
              // draw the legend
              TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
                legend->SetFillColorAlpha(kWhite, 0.);
                legend->SetBorderSize(0);
                  legend->SetTextFont(42);
              legend->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitV2_2->GetChisquare(),fitV2_2->GetNDF());
             legend->AddEntry(fitV2_2,message);

                if(res->CovMatrixStatus() == 3){
                     //  sprintf(message,"The fit is a success");
                }
                else{
                       sprintf(message,"The fit is a failure");
                    legend->AddEntry(fitV2_2,message);
                }
           // legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
            legend->AddEntry(backFcnV2_2,"Background");
              legend->AddEntry(V2JPsiTkl,"Data","lpe");
              legend->Draw();
        
        sprintf(CanvasName,"%s/V2_Ext2_noZYAM.pdf",CanvasFolderName);
        c16_noZYAM->SaveAs(CanvasName);
    }
    
    // Ext 3
    
    double b0_JPsi = 0;
    double c0_JPsi = 0;
    double c2_JPsi = 0;
    double eb0_JPsi = 0;
    double ec0_JPsi = 0;
    double ec2_JPsi = 0;
    
    TCanvas*c316=new TCanvas();
               c316->SetTitle("Extraction method 3");
               //V2_2 Jpsi-tkl wrt Mass fit (Extraction method 3)
            
            //Definition and design of function on V2 plot
           c316->Divide(2,2);
    
{
         c316->cd(1);
         TF1 *fitb0_3 = new TF1("fitb0_3",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
         fitb0_3->SetNpx(500);
         fitb0_3->SetLineWidth(4);
         fitb0_3->SetLineColor(kMagenta);
         fitb0_3->SetParameters(par);
         TF1 *backFcnb0_3 = new TF1("backFcnb0_3",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
         backFcnb0_3->SetLineColor(kRed);
         TF1 *signalFcnJPsib0_3 = new TF1("signalFcnJPsib0_3",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
         signalFcnJPsib0_3->SetLineColor(kBlue);
         signalFcnJPsib0_3->SetNpx(500);
         
            TVirtualFitter::Fitter(baselines0)->SetMaxIterations(10000);
            TVirtualFitter::Fitter(baselines0)->SetPrecision();
             for(int i=0; i<numParameters; i++){
                 fitb0_3->FixParameter(i,par[i]);
             }
    
    for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
        string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "B0");
        const char *titular = title.c_str();
        fitb0_3->SetParName(index,titular);
    }
            
//            fitb0_3->SetParName(0,"Norm_{J/#psi}");
//            fitb0_3->SetParName(1,"M_{J/#psi}");
//            fitb0_3->SetParName(2,"Sigma_{J/#psi}");
//            fitb0_3->SetParName(3,"a_{1}");
//            fitb0_3->SetParName(4,"n_{1}");
//            fitb0_3->SetParName(5,"a_{2}");
//            fitb0_3->SetParName(6,"n_{2}");
//            fitb0_3->SetParName(7,"Norm_{#Psi(2S)}");
//            fitb0_3->SetParName(8,"Norm_{TailLowM}");
//            fitb0_3->SetParName(9,"Exp_{TailLowM}");
//            fitb0_3->SetParName(10,"Norm_{TailHighM}");
//            fitb0_3->SetParName(11,"Exp_{TailHighM}");
//             fitb0_3->SetParName(12,"B0 J/#psi");
//             fitb0_3->SetParName(13,"B0 Bkg M2");
//         fitb0_3->SetParName(14,"B0 Bkg M1");
//         fitb0_3->SetParName(15,"B0 Nkg M0");
    
         
         //Fit of b0
             res = baselines0->Fit("fitb0_3","SBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)baselines0->FindObject("stats");
            st->SetX1NDC(0.75); //new x start position
            st->SetY1NDC(0.75); //new x end position
             Double_t para3[numParameters+numParametersBkgV2+1];
             fitb0_3->GetParameters(para3);
             backFcnb0_3->SetParameters(para3);
             signalFcnJPsib0_3->SetParameters(para3);
         
        b0_JPsi = para3[numParameters];
        eb0_JPsi = fitb0_3->GetParError(numParameters);
         
         fitb0_3->Draw("same");
       //  signalFcnJPsib0_2->Draw("same");
         backFcnb0_3->Draw("same");
           // draw the legend
           TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
           legend->SetFillColorAlpha(kWhite, 0.);
           legend->SetBorderSize(0);
             legend->SetTextFont(42);
           legend->SetTextSize(0.03);
         Char_t message[80];
         sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitb0_3->GetChisquare(),fitb0_3->GetNDF());
          legend->AddEntry(fitb0_3,message);

             if(res->CovMatrixStatus() == 3){
                   // sprintf(message,"The fit is a success");
             }
             else{
                    sprintf(message,"The fit is a failure");
                 legend->AddEntry(fitb0_3,message);
             }
        // legend->AddEntry(signalFcnJPsib0_2,"JPsi signal");
         legend->AddEntry(backFcnb0_3,"Background");
           legend->AddEntry(baselines0,"Data","lpe");
           legend->Draw();
    
    
    
    
}
    
    {
             c316->cd(2);
             TF1 *fitc0_3 = new TF1("fitc0_3",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             fitc0_3->SetNpx(500);
             fitc0_3->SetLineWidth(4);
             fitc0_3->SetLineColor(kMagenta);
             fitc0_3->SetParameters(par);
             TF1 *backFcnc0_3 = new TF1("backFcnc0_3",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             backFcnc0_3->SetLineColor(kRed);
             TF1 *signalFcnJPsic0_3 = new TF1("signalFcnJPsic0_3",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             signalFcnJPsic0_3->SetLineColor(kBlue);
             signalFcnJPsic0_3->SetNpx(500);
             
                TVirtualFitter::Fitter(coefficients0)->SetMaxIterations(10000);
                TVirtualFitter::Fitter(coefficients0)->SetPrecision();
                 for(int i=0; i<numParameters; i++){
                     fitc0_3->FixParameter(i,par[i]);
                 }
        
        for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
            string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "C0");
            const char *titular = title.c_str();
            fitc0_3->SetParName(index,titular);
        }
        
//                fitc0_3->SetParName(0,"Norm_{J/#psi}");
//                fitc0_3->SetParName(1,"M_{J/#psi}");
//                fitc0_3->SetParName(2,"Sigma_{J/#psi}");
//                fitc0_3->SetParName(3,"a_{1}");
//                fitc0_3->SetParName(4,"n_{1}");
//                fitc0_3->SetParName(5,"a_{2}");
//                fitc0_3->SetParName(6,"n_{2}");
//                fitc0_3->SetParName(7,"Norm_{#Psi(2S)}");
//                fitc0_3->SetParName(8,"Norm_{TailLowM}");
//                fitc0_3->SetParName(9,"Exp_{TailLowM}");
//                fitc0_3->SetParName(10,"Norm_{TailHighM}");
//                fitc0_3->SetParName(11,"Exp_{TailHighM}");
//                 fitc0_3->SetParName(12,"c0 J/#psi");
//                 fitc0_3->SetParName(13,"c0 Bkg M2");
//             fitc0_3->SetParName(14,"c0 Bkg M1");
//             fitc0_3->SetParName(15,"c0 Nkg M0");
             
             
             //Fit of c0
                 res = coefficients0->Fit("fitc0_3","SBMERI+","ep");
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)coefficients0->FindObject("stats");
                st->SetX1NDC(0.75); //new x start position
                st->SetY1NDC(0.75); //new x end position
                 Double_t para3[numParameters+numParametersBkgV2+1];
                 fitc0_3->GetParameters(para3);
                 backFcnc0_3->SetParameters(para3);
                 signalFcnJPsic0_3->SetParameters(para3);
             
            c0_JPsi = para3[numParameters];
            ec0_JPsi = fitc0_3->GetParError(numParameters);
             
             fitc0_3->Draw("same");
           //  signalFcnJPsic0_2->Draw("same");
             backFcnc0_3->Draw("same");
               // draw the legend
               TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
               legend->SetFillColorAlpha(kWhite, 0.);
               legend->SetBorderSize(0);
                 legend->SetTextFont(42);
               legend->SetTextSize(0.03);
             Char_t message[80];
             sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitc0_3->GetChisquare(),fitc0_3->GetNDF());
              legend->AddEntry(fitc0_3,message);

                 if(res->CovMatrixStatus() == 3){
                       // sprintf(message,"The fit is a success");
                 }
                 else{
                        sprintf(message,"The fit is a failure");
                     legend->AddEntry(fitc0_3,message);
                 }
            // legend->AddEntry(signalFcnJPsic0_2,"JPsi signal");
             legend->AddEntry(backFcnc0_3,"Background");
               legend->AddEntry(coefficients0,"Data","lpe");
               legend->Draw();
        
        
        
        
    }
    
    {
             c316->cd(3);
             TF1 *fitc2_3 = new TF1("fitc2_3",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             fitc2_3->SetNpx(500);
             fitc2_3->SetLineWidth(4);
             fitc2_3->SetLineColor(kMagenta);
             fitc2_3->SetParameters(par);
             TF1 *backFcnc2_3 = new TF1("backFcnc2_3",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             backFcnc2_3->SetLineColor(kRed);
             TF1 *signalFcnJPsic2_3 = new TF1("signalFcnJPsic2_3",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             signalFcnJPsic2_3->SetLineColor(kBlue);
             signalFcnJPsic2_3->SetNpx(500);
             
                TVirtualFitter::Fitter(coefficients2)->SetMaxIterations(10000);
                TVirtualFitter::Fitter(coefficients2)->SetPrecision();
                 for(int i=0; i<numParameters; i++){
                     fitc2_3->FixParameter(i,par[i]);
                 }

        for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
            string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "C2");
            const char *titular = title.c_str();
            fitc2_3->SetParName(index,titular);
        }
                
//                fitc2_3->SetParName(0,"Norm_{J/#psi}");
//                fitc2_3->SetParName(1,"M_{J/#psi}");
//                fitc2_3->SetParName(2,"Sigma_{J/#psi}");
//                fitc2_3->SetParName(3,"a_{1}");
//                fitc2_3->SetParName(4,"n_{1}");
//                fitc2_3->SetParName(5,"a_{2}");
//                fitc2_3->SetParName(6,"n_{2}");
//                fitc2_3->SetParName(7,"Norm_{#Psi(2S)}");
//                fitc2_3->SetParName(8,"Norm_{TailLowM}");
//                fitc2_3->SetParName(9,"Exp_{TailLowM}");
//                fitc2_3->SetParName(10,"Norm_{TailHighM}");
//                fitc2_3->SetParName(11,"Exp_{TailHighM}");
//                 fitc2_3->SetParName(12,"c2 J/#psi");
//                 fitc2_3->SetParName(13,"c2 Bkg M2");
//             fitc2_3->SetParName(14,"c2 Bkg M1");
//             fitc2_3->SetParName(15,"c2 Nkg M0");
             
             //Fit of c2
                 res = coefficients2->Fit("fitc2_3","SBMERI+","ep");
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)coefficients2->FindObject("stats");
                st->SetX1NDC(0.75); //new x start position
                st->SetY1NDC(0.75); //new x end position
                 Double_t para3[numParameters+numParametersBkgV2+1];
                 fitc2_3->GetParameters(para3);
                 backFcnc2_3->SetParameters(para3);
                 signalFcnJPsic2_3->SetParameters(para3);
             
            c2_JPsi = para3[numParameters];
            ec2_JPsi = fitc2_3->GetParError(numParameters);
             
             fitc2_3->Draw("same");
           //  signalFcnJPsic2_2->Draw("same");
             backFcnc2_3->Draw("same");
               // draw the legend
               TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
               legend->SetFillColorAlpha(kWhite, 0.);
               legend->SetBorderSize(0);
                 legend->SetTextFont(42);
               legend->SetTextSize(0.03);
             Char_t message[80];
             sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitc2_3->GetChisquare(),fitc2_3->GetNDF());
              legend->AddEntry(fitc2_3,message);

                 if(res->CovMatrixStatus() == 3){
                     //   sprintf(message,"The fit is a success");
                 }
                 else{
                        sprintf(message,"The fit is a failure");
                     legend->AddEntry(fitc2_3,message);
                 }
            // legend->AddEntry(signalFcnJPsic2_2,"JPsi signal");
             legend->AddEntry(backFcnc2_3,"Background");
               legend->AddEntry(coefficients2,"Data","lpe");
               legend->Draw();
        
        
        
        
    }
    
    {
        c316->cd(4);
        
         Char_t message[80];
          sprintf(message,"V2_3 J/#psi-tkl = #frac{%.4f}{%.4f + %.4f} = %.5f +- %.5f",c2_JPsi,c0_JPsi, b0_JPsi,c2_JPsi/(c0_JPsi + b0_JPsi),(c2_JPsi/(c0_JPsi + b0_JPsi)*sqrt(pow(ec2_JPsi/c2_JPsi,2)+pow(sqrt(pow(ec0_JPsi,2)+pow(eb0_JPsi,2))/(c0_JPsi + b0_JPsi),2))));
        TPaveText *pave = new TPaveText(0.2,0.2,0.8,0.8,"brNDC");
        pave->AddText(message);
        pave->SetTextFont(42);
        pave->SetTextSize(0.05);
        pave->SetBorderSize(0);
        pave->SetFillStyle(0);
           pave->Draw();
        
    }
    
    sprintf(CanvasName,"%s/V2_Ext3.pdf",CanvasFolderName);
    c316->SaveAs(CanvasName);

    //Ext 2 Pt Binned
    
     TCanvas *cinvmassPtBinned = new TCanvas("cinvmassPtBinned","Fitting All Phi - Different Centralities - Pt Binned",10,10,700,500);
        cinvmassPtBinned->SetTitle("Inv Mass Fits - Pt Binned");
        // ROOT THAT Inv mass fits
        cinvmassPtBinned->Divide(NbPtBins,3);
     
    TCanvas *cinvmassPtBinnedCombine = new TCanvas("cinvmassPtBinnedCombine","Fitting All Phi - Different Centralities - Pt Binned Combine",10,10,700,500);
    cinvmassPtBinnedCombine->SetTitle("Inv Mass Fits - Pt Binned Combine");
    // ROOT THAT Inv mass fits
    cinvmassPtBinnedCombine->Divide(NbPtBins,3);
        
    if(PtBinned){
        TCanvas*cpool = new TCanvas();
        cpool->SetTitle("Deviation inv mass fit");
        cpool->Divide(2,4);
        
        for(int ptbin=0;ptbin<NbPtBins;ptbin++){
            cout << "ptbincombine "<< ptbin;
        sprintf(histoname,Form("bin%d_0",ptbin+1)); // ZZZZZZZ
            if(!CombineFits){
        res = FittingAllInvMassBin(histoname, cinvmassPtBinned, ptbin, SignalF, BackgroundF, minmass, maxmass, ratsigma);
            }
        double par[16];
        
        // Parameter setting (either from mass fit either start values)
        for(int i=0; i <numParameters+numParametersBkgV2+1; i++){
            if(!CombineFits){
                par[i] = res->Parameter(i);
                if(i>=numParameters){
                    par[i] = 1;
                }
            }
            if(CombineFits){
                par[i] = startvalues[i];
            }
        }
            
        {
            ofstream myfiletxt;
            myfiletxt.open("/tmp/pulls.txt");
            TFile *file0 = new TFile(FitFileName);
            TH1F *histo = (TH1F*)file0->Get(histoname);
            double params[numParameters];
            for(int i=0; i<numParameters; i++){
                params[i] = par[i];
            }
            TF1 *fitFcn = new TF1("fitFcn",MassFunction,MinInvMassFit,MaxInvMassFit,numParameters);
              fitFcn->SetNpx(500);
              fitFcn->SetLineWidth(4);
              fitFcn->SetLineColor(kMagenta);
              fitFcn->SetParameters(params);
            
            //Pour chaue bin
                //Trouver la valeur du centre du bin pour le fit et pour les données
                //Calculer le pool
                // Ajouter le point à un histo

            for(int bin=0; bin<histo->GetNbinsX(); bin++){
                double fit_prediction = fitFcn->Eval(MinInvMass + (MaxInvMass-MinInvMass)*(0.5+bin)/NbinsDimuInvMass);
                double data = histo->GetBinContent(bin+1);
                double error = histo->GetBinError(bin+1);
                double pool = (data-fit_prediction)/error;
                
                if(bin < 10){
                myfiletxt << "Bin " << bin << " : Fit - " << fit_prediction << " Data - " << data << " Error - " << error << " Pool value - " << pool <<endl;
                }
                
                hpool->SetBinContent(bin+1,pool);
            }
            myfiletxt.close();
            
            cpool->cd(ptbin+1);
            char title[50];
            sprintf(title,Form("Pulls_PtBin_[%.2f,%.2f] GeV",PtBins[ptbin],PtBins[ptbin+1]));
            hpool->SetTitle(title);
            hpool->DrawCopy();
            
            TCanvas*cmassandpoolbinned = new TCanvas();
            TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
            TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
            pad1->SetBottomMargin(0.00001);
            pad1->SetBorderMode(0);
            pad2->SetTopMargin(0.00001);
            pad2->SetBottomMargin(0.1);
            pad2->SetBorderMode(0);
            pad1->Draw();
            pad2->Draw();
            sprintf(title,Form("Mass fit and Pulls PtBin_[%.2f,%.2f] GeV",PtBins[ptbin],PtBins[ptbin+1]));
            cmassandpoolbinned->SetTitle(title);
            pad1->cd();
            TVirtualPad *padvirt = cinvmassPtBinned->GetPad(ptbin+1);
            double x1 = padvirt->GetAbsXlowNDC();
            double x2 = x1 + padvirt->GetAbsWNDC();
            double y1 = padvirt->GetAbsYlowNDC();
            double y2 = y1 + padvirt->GetAbsHNDC();
            padvirt->SetPad(0,0.33,1,1);
            padvirt->DrawClone();
            padvirt->SetPad(x1,y1,x2,y2);
            pad2->cd();
            hpool->DrawCopy();
            
        }
        
        //Definition and design of function on V2 plot
        c16PtBinned->cd(ptbin+1);
        TF1 *fitV2_2PtBinned = new TF1("fitV2_2PtBinned",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        fitV2_2PtBinned->SetNpx(500);
        fitV2_2PtBinned->SetLineWidth(2);
        fitV2_2PtBinned->SetLineColor(kAzure+3);
        fitV2_2PtBinned->SetParameters(par);
        TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        backFcnV2_2->SetLineColor(kRed);
            backFcnV2_2->SetLineStyle(9);
            backFcnV2_2->SetLineWidth(2);
        TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        signalFcnJPsiV2_2->SetLineColor(kBlue);
        signalFcnJPsiV2_2->SetNpx(500);
            
            gStyle->SetErrorX(0.);

            V2JPsiTklPtBinned[ptbin]->SetMarkerColor(kBlack);
            V2JPsiTklPtBinned[ptbin]->SetMarkerStyle(20);
            V2JPsiTklPtBinned[ptbin]->SetLineColor(kBlack);
            V2JPsiTklPtBinned[ptbin]->GetYaxis()->SetRangeUser(-0.004, 0.009);
            V2JPsiTklPtBinned[ptbin]->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/c^{2})");
            V2JPsiTklPtBinned[ptbin]->GetYaxis()->SetTitle("#it{V}_{2}{dimuon-tracklet,sub}");
            
        
          if(!CombineFits){
           TVirtualFitter::Fitter(V2JPsiTklPtBinned[ptbin])->SetMaxIterations(10000);
           TVirtualFitter::Fitter(V2JPsiTklPtBinned[ptbin])->SetPrecision();
            for(int i=0; i<numParameters; i++){
                fitV2_2PtBinned->FixParameter(i,par[i]);
            }
            }
            
            for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "V2_2");
                const char *titular = title.c_str();
                fitV2_2PtBinned->SetParName(index,titular);
            }
            
//           fitV2_2PtBinned->SetParName(0,"Norm_{J/#psi}");
//           fitV2_2PtBinned->SetParName(1,"M_{J/#psi}");
//           fitV2_2PtBinned->SetParName(2,"Sigma_{J/#psi}");
//           fitV2_2PtBinned->SetParName(3,"a_{1}");
//           fitV2_2PtBinned->SetParName(4,"n_{1}");
//           fitV2_2PtBinned->SetParName(5,"a_{2}");
//           fitV2_2PtBinned->SetParName(6,"n_{2}");
//           fitV2_2PtBinned->SetParName(7,"Norm_{#Psi(2S)}");
//           fitV2_2PtBinned->SetParName(8,"Norm_{TailLowM}");
//           fitV2_2PtBinned->SetParName(9,"Exp_{TailLowM}");
//           fitV2_2PtBinned->SetParName(10,"Norm_{TailHighM}");
//           fitV2_2PtBinned->SetParName(11,"Exp_{TailHighM}");
//            fitV2_2PtBinned->SetParName(12,"V2_2 J/#psi");
//            fitV2_2PtBinned->SetParName(13,"V2_2 Bkg M2");
//        fitV2_2PtBinned->SetParName(14,"V2_2 Bkg M1");
//        fitV2_2PtBinned->SetParName(15,"V2_2 Bkg M0");

        
           double minParams[numParameters+numParametersBkgV2+1];
           double parErrors[numParameters+numParametersBkgV2+1];
            Double_t para[numParameters+numParametersBkgV2+1];
            
            // Combined fit of inv mass and V2_2
            if(CombineFits){
                //Fixme startvalues et tout ça // Just not use combine fits for now
                double startvalues2[16] = {100, 3.096916, 0.07,  0.9, 10, 2, 15, 1,    10,      0.01,     100,      1.4,      1,         1,   1,   1};
                double lowvalues2[16] =   {0.1,    3.05,     0.03, 0.7, 1,  1, 5,  0.001,  0.1,    -20,   0.001,        0.01, -10, -10, -10, -10};
                double highvalues2[16] = {1000000, 3.15,     0.1,  1.1, 30, 5, 25, 10000, 1000000, 100,      10000000,  100,    10,        10,  10,  10};
                
                TBackCompFitter * virminuit = (TBackCompFitter *) TVirtualFitter::Fitter(0,numParameters+numParametersBkgV2+1);
                   for (int i = 0; i < numParameters+numParametersBkgV2+1; ++i) {
                     virminuit->SetParameter(i, fitV2_2PtBinned->GetParName(i), fitV2_2PtBinned->GetParameter(i), startvalues2[i], lowvalues2[i],highvalues2[i]);
                   }
                virminuit->FixParameter(1);
                virminuit->FixParameter(2);
                virminuit->FixParameter(3);
                virminuit->FixParameter(4);
                virminuit->FixParameter(5);
                virminuit->FixParameter(6);
                actuelptbin = ptbin;
                   virminuit->SetFCN(FcnCombinedAllMass);

                   double arglist[100];
                   arglist[0] = 1;
                   // set print level
                   virminuit->ExecuteCommand("SET PRINT",arglist,2);

                   // minimize
                   arglist[0] = 5000; // number of function calls
                   arglist[1] = 0.01; // tolerance
                   virminuit->ExecuteCommand("MIGRAD",arglist,2);
                
                gStyle->SetOptStat("n");
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
              //  resu.GetCovarianceMatrix().Print();

                   //get result
                   for (int i = 0; i < numParameters+numParametersBkgV2+1; ++i) {
                     para[i] = virminuit->GetParameter(i);
                     parErrors[i] = virminuit->GetParError(i);
                   }
                   double chi2, edm, errdef;
                   int nvpar, nparx;
                   virminuit->GetStats(chi2,edm,errdef,nvpar,nparx);

                   fitV2_2PtBinned->SetParameters(para);
                   fitV2_2PtBinned->SetParErrors(parErrors);
//                fitJustV2_2->SetParameters(&minParams[12]);
//                fitJustV2_2->SetParErrors(&parErrors[12]);
                   fitV2_2PtBinned->SetChisquare(chi2);
                   int ndf = npfits-nvpar;
                   fitV2_2PtBinned->SetNDF(ndf);
            }
        
        //Fit of V2
            if(!CombineFits){
            res = V2JPsiTklPtBinned[ptbin]->Fit("fitV2_2PtBinned","SBMERI+","ep");
            
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)V2JPsiTklPtBinned[ptbin]->FindObject("stats");
            st->SetX1NDC(0.75); //new x start position
            st->SetY1NDC(0.75); //new x end position
            fitV2_2PtBinned->GetParameters(para);
            backFcnV2_2->SetParameters(para);
            signalFcnJPsiV2_2->SetParameters(para);
            }
        
            if(CombineFits){
                backFcnV2_2->SetParameters(para);
                signalFcnJPsiV2_2->SetParameters(para);
            }
        
        fitV2_2PtBinned->Draw("same");
            if(CombineFits){
              //  V2JPsiTkl->GetListOfFunctions()->Add(fitJustV2_2);
            }
      //  signalFcnJPsiV2_2->Draw("same");
        backFcnV2_2->Draw("same");
            gStyle->SetOptFit(1011);
          // draw the legend
          TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
          legend->SetFillColorAlpha(kWhite, 0.);
          legend->SetBorderSize(0);
            legend->SetTextFont(42);
          legend->SetTextSize(0.035);
            legend->AddEntry(V2JPsiTklPtBinned[ptbin],"Data","lpe");
        Char_t message[80];
            sprintf(message,"Total fit");
        //sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitV2_2PtBinned->GetChisquare(),fitV2_2PtBinned->GetNDF());
         legend->AddEntry(fitV2_2PtBinned,message);
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
                  // sprintf(message,"The fit is a success");
            }
            else{
                   sprintf(message,"The fit is a failure");
                legend->AddEntry(fitV2_2PtBinned,message);
            }
            }
       // legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
        legend->AddEntry(backFcnV2_2,"Background");
            sprintf(message,"#it{V}_{2,2} = %.3e #pm %.3e",para[numParameters],fitV2_2PtBinned->GetParError(numParameters));
            legend->AddEntry(fitV2_2PtBinned,message,"");
          legend->Draw();
            
            TLegend *legendov3=new TLegend(0.6,0.70,0.90,0.90);
                      legendov3->SetFillColorAlpha(kWhite, 0.);
                      legendov3->SetBorderSize(0);
                       legendov3->SetTextFont(42);
                       legendov3->SetTextSize(0.04);
            sprintf(message,"V0M (0-5%%)-(40-100%%)");
             legendov3->AddEntry(fitV2_2PtBinned,message,"");
            sprintf(message,"pp, #sqrt{#it{s}_{NN}} = 13 TeV");
            legendov3->AddEntry(fitV2_2PtBinned,message,"");
            sprintf(message,"2.5 < #it{y}_cms < 4.0");
            legendov3->AddEntry(fitV2_2PtBinned,message,"");
            sprintf(message,"1.5 < |#it{#Delta#eta}| < 5.0");
            legendov3->AddEntry(fitV2_2PtBinned,message,"");
            sprintf(message,"3 < #it{p}_{T} < 4 GeV/#it{c}");
            legendov3->AddEntry(fitV2_2PtBinned,message,"");
                      legendov3->Draw();
            
            TLegend *legendov4=new TLegend(0.6,0.70,0.90,0.90);
                      legendov4->SetFillColorAlpha(kWhite, 0.);
                      legendov4->SetBorderSize(0);
                       legendov4->SetTextFont(42);
                       legendov4->SetTextSize(0.04);
            sprintf(message,"ALICE Preliminary");
             legendov4->AddEntry(fitV2_2PtBinned,message,"");
                      legendov4->Draw();

            V2_Ext2[ptbin] = para[numParameters];
            errV2_Ext2[ptbin] = fitV2_2PtBinned->GetParError(numParameters);
            
            sprintf(CanvasName,"%s/V2_Ext2_PtBinned.pdf",CanvasFolderName);
            c16PtBinned->SaveAs(CanvasName);
            
            //Plot of inv mass
            if(CombineFits){
                //AEUGGG
                    cinvmassPtBinnedCombine->cd(ptbin+1);
                   cinvmassPtBinnedCombine->SetFillColor(33);
                   cinvmassPtBinnedCombine->SetFrameFillColor(41);
                   cinvmassPtBinnedCombine->SetGrid();
                   TH1F *histo = (TH1F*)hPtWrtMassInvSliced[0][ptbin];
                   // create a TF1 with the range from 0 to 3 and 6 parameters
                   TF1 *fitFcn = new TF1("fitFcn",MassFunction,MinInvMassFit,MinInvMassFit,numParameters);
                   fitFcn->SetNpx(500);
                   fitFcn->SetLineWidth(4);
                   fitFcn->SetLineColor(kMagenta);

                
                for(int index=0; index<numParameters; index++){
                    string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "V2_2");
                    const char *titular = title.c_str();
                    fitFcn->SetParName(index,titular);
                }
                
//                    fitFcn->SetParName(0,"Norm_{JPsi}");
//                    fitFcn->SetParName(1,"M_{JPsi}");
//                    fitFcn->SetParName(2,"Sigma_{JPsi}");
//                    fitFcn->SetParName(3,"a_{1}");
//                    fitFcn->SetParName(4,"n_{1}");
//                    fitFcn->SetParName(5,"a_{2}");
//                    fitFcn->SetParName(6,"n_{2}");
//                    fitFcn->SetParName(7,"Norm_{Psi2S}");
//                    fitFcn->SetParName(8,"Norm_{TailLowM}");
//                    fitFcn->SetParName(9,"Exp_{TailLowM}");
//                    fitFcn->SetParName(10,"Norm_{TailHighM}");
//                    fitFcn->SetParName(11,"Exp_{TailHighM}");

                   // improve the picture:
                   TF1 *backFcn = new TF1("backFcn",UsedBackground,MinInvMassFit,MinInvMassFit,numParameters+numParametersBkgV2+1);
                   backFcn->SetLineColor(kRed);
                   TF1 *signalFcnJPsi = new TF1("signalFcnJPsi",UsedJPsiSignal,MinInvMassFit,MinInvMassFit,numParameters-numParametersBkg);
                   TF1 *signalFcnPsi2S = new TF1("signalFcnPsi2S",UsedPsi2SSignal,MinInvMassFit,MinInvMassFit,numParameters-numParametersBkg);
                   TPaveText *pave = new TPaveText(0.15,0.5,0.3,0.65,"brNDC");
                   signalFcnJPsi->SetLineColor(kBlue);
                   signalFcnJPsi->SetNpx(500);
                    signalFcnPsi2S->SetLineColor(kGreen);
                    signalFcnPsi2S->SetNpx(500);
                    histo->GetListOfFunctions()->Add(fitFcn);
                   // writes the fit results into the par array
                    gStyle->SetOptFit(1011);
                    fitFcn->SetParameters(para);
                    fitFcn->SetParErrors(parErrors);
                   signalFcnJPsi->SetParameters(para);
                    signalFcnPsi2S->SetParameters(para);
                 //  auto covMatrix = res->GetCovarianceMatrix();
                   std::cout << "Covariance matrix from the fit ";
                  // resu.GetCovarianceMatrix().Print();
                   Double_t integral = (signalFcnJPsi->Integral(2.1,3.45))*(MaxInvMass-MinInvMass)/NbinsDimuInvMass;
                    std::cout << "STATUS COV " << resu.CovMatrixStatus() <<endl;
                Double_t integralerror = 0;//(signalFcnJPsi->IntegralError(2.1,3.45,signalFcnJPsi->GetParameters(), resu.GetCovarianceMatrix().GetSub(0,8,0,8).GetMatrixArray() ))/0.01;
                  //  std::cout << "Erreur integrale " << integralerror <<endl;
                    
                    std::cout << "Fitted " << histoname << std::endl;
                    std::cout << "Nb JPsi total: " << integral << std::endl;
                 //   std::cout << "integral error: " << integralerror << std::endl;
                histo->Draw();
                fitFcn->Draw("same");
                   signalFcnJPsi->Draw("same");
                    signalFcnPsi2S->Draw("same");
                   backFcn->SetParameters(para);
                   backFcn->Draw("same");
                   // draw the legend
                    char str[50];
                    sprintf(str, "N_{JPsi} %f +/- %f", integral, integralerror);
                   pave->AddText(str);
                    sprintf(str, "M_{Psi2S} = %f, Sig_{Psi2S} = %f", mPsip, ratSigma*para[2]);
                    pave->AddText(str);
                   pave->Draw();
                   TLegend *legend=new TLegend(0.12,0.75,0.60,0.90);
                   legend->SetTextFont(61);
                   legend->SetTextSize(0.03);
                    Char_t message[80];
                    sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcn->GetChisquare(),fitFcn->GetNDF());
                     legend->AddEntry(fitFcn,message);
                    if(resu.CovMatrixStatus() == 3){
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
            }
            
            
            
            
        }
    }
    
    //FIXME EXT 2 ATLAS
    {
        {
            
            TCanvas*c16PtBinnedFin=new TCanvas();
            // ROOT THAT Extraction 2 plot
            c16PtBinnedFin->SetTitle("Extraction method 2 PtBinned ATLAS End");
            c16PtBinnedFin->Divide(2,3);
            
            
                for(int ptbin=0;ptbin<NbPtBins;ptbin++){
                    
                    // Plot du V2 JPsi Tracklet PtBinned
                    TCanvas*c16PtBinnedATLAS=new TCanvas();
                    // ROOT THAT Extraction 2 plot
                    c16PtBinnedATLAS->SetTitle("Extraction method 2 ATLAS PtBinned");
                    c16PtBinnedATLAS->Divide(5,2);
                    
                    for(int i=1; i<=NbinsInvMass; i++){
                        
                        TF1 *fitFcnV2_PRLTemplateAlt = new TF1("fitFcnV2_PRLTemplateAlt",PRLTemplatePtBinnedAlt,0,TMath::Pi(),4);//-TMath::Pi()/2,1.5*TMath::Pi(),2);
                        fitFcnV2_PRLTemplateAlt->SetNpx(500);
                        fitFcnV2_PRLTemplateAlt->SetLineWidth(4);
                        fitFcnV2_PRLTemplateAlt->SetLineColor(kRed+1);
                        
                        c16PtBinnedATLAS->cd(i);
                        cout << "Fit of the difference to extract Fourier wrt mass. Mass bin: " << i <<" Pt bin: "<<ptbin <<endl;
                    // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
                        TH1F *histo = Yield_Central_MassBinPtBinned[i-1][ptbin];
                        
                        Double_t params[4] = { static_cast<Double_t>( ptbin ),1,1,static_cast<Double_t>( i-1 )};
                                fitFcnV2_PRLTemplateAlt->SetParameters(params);
                                 TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
                                 TVirtualFitter::Fitter(histo)->SetPrecision();
                                TVirtualFitter::Fitter(histo)->SetFCN(ChisquarePRLTemplatePtBinnedAlt);
                              // gStyle->SetOptFit(1011);
                               //  histo->Fit("fitFcn","0");
                                // second try: set start values for some parameters

                        fitFcnV2_PRLTemplateAlt->SetParName(0,"ptbin");
                                 fitFcnV2_PRLTemplateAlt->SetParName(1,"V2");
                                   fitFcnV2_PRLTemplateAlt->SetParName(2,"F");
                        fitFcnV2_PRLTemplateAlt->SetParName(3,"MassBin");
                           
                         //  fitFcnV2_PRLTemplate->SetParLimits(0,0.0001,5000000);
                           fitFcnV2_PRLTemplateAlt->FixParameter(0,ptbin);
                      //  fitFcnV2_PRLTemplateAlt->FixParameter(1,0);
                      //  fitFcnV2_PRLTemplateAlt->SetParLimits(1,-1,1);
                           fitFcnV2_PRLTemplateAlt->SetParLimits(2,0.1,2);
                                   fitFcnV2_PRLTemplateAlt->FixParameter(3,i-1);
                        
                        TFitResultPtr res = histo->Fit("fitFcnV2_PRLTemplateAlt","USBMERI+","ep");
                             gStyle->SetOptStat("n");
                             gStyle->SetOptFit(1011);
                             TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
                             st->SetX1NDC(0.8); //new x start position
                             st->SetY1NDC(0.8); //new x end position
                        
                            double chi2, edm, errdef;
                            int nvpar, nparx;
                             TVirtualFitter::Fitter(histo)->GetStats(chi2,edm,errdef,nvpar,nparx);
                             fitFcnV2_PRLTemplateAlt->SetChisquare(chi2);
                             int ndf = npfits-nvpar;
                             fitFcnV2_PRLTemplateAlt->SetNDF(ndf);
                        
                            Double_t par[4];
                            fitFcnV2_PRLTemplateAlt->GetParameters(par);
                             // improve the pictu
                           //   std::cout << "integral error: " << integralerror << std::endl;
                             fitFcnV2_PRLTemplateAlt->Draw("same");
                             // draw the legend
                             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                             legend->SetFillColorAlpha(kWhite, 0.);
                             legend->SetBorderSize(0);
                               legend->SetTextFont(42);
                             legend->SetTextSize(0.03);
                               Char_t message[80];
                               sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplateAlt->GetChisquare(),fitFcnV2_PRLTemplateAlt->GetNDF());
                               legend->AddEntry(fitFcnV2_PRLTemplateAlt,message);
                        if(res->CovMatrixStatus() == 3){
                                //   sprintf(message,"The fit is a success");
                               }
                               else{
                                   sprintf(message,"The fit is a failure");
                                   legend->AddEntry(fitFcnV2_PRLTemplateAlt,message);
                               }
                             legend->AddEntry(histo,"Data","lpe");
                             legend->Draw();
                        
                        V22ATLASPtBinned[ptbin]->SetBinContent(i, par[1]);
                        V22ATLASPtBinned[ptbin]->SetBinError(i,fitFcnV2_PRLTemplateAlt->GetParError(1));
                        FATLASPtBinned[ptbin]->SetBinContent(i, par[2]);
                        FATLASPtBinned[ptbin]->SetBinError(i,fitFcnV2_PRLTemplateAlt->GetParError(2));
                    //  fitFcnV2_2->Draw("same");
                        cout << "STATUS COV : " << res->CovMatrixStatus() <<endl;
                      // draw the legend
        //              TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
        //              legend->SetTextFont(61);
        //              legend->SetTextSize(0.03);
        //                Char_t message[80];
        //                sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
        //                legend->AddEntry(fitFcnV2_2,message);
        //                if(res->CovMatrixStatus() == 3){
        //                    sprintf(message,"The fit is a success");
        //                }
        //                else{
        //                    sprintf(message,"The fit is a failure");
        //                }
        //                legend->AddEntry(fitFcnV2_2,message);
        //              legend->AddEntry(histo,"Data","lpe");
        //              legend->Draw();
                    //FIXME Ajouter d'autres méthodes ok
                    }
                    sprintf(CanvasName,"%s/ATLASExt2_PtBin%i.pdf",CanvasFolderName,ptbin);
                    c16PtBinnedATLAS->SaveAs(CanvasName);
                }
            
            
            TCanvas *cinvmassPtBinnedATLAS = new TCanvas("cinvmassPtBinnedATLAS","Fitting All Phi - Different Centralities - Pt Binned",10,10,700,500);
            cinvmassPtBinnedATLAS->SetTitle("Inv Mass Fits - Pt Binned ATLAS");
            // ROOT THAT Inv mass fits
            cinvmassPtBinnedATLAS->Divide(NbPtBins,3);
            
            for(int ptbin=0;ptbin<NbPtBins;ptbin++){
                sprintf(histoname,Form("bin%d_0",ptbin+1)); // ZZZZZZZ
                TFitResultPtr res = FittingAllInvMassBin(histoname, cinvmassPtBinned, ptbin, SignalF, BackgroundF, minmass, maxmass, ratsigma);
                           
                       double par[16];
                       
                       // Parameter setting (either from mass fit either start values)
                       for(int i=0; i <numParameters+numParametersBkgV2+1; i++){
                               par[i] = res->Parameter(i);
                               if(i>=numParameters){
                                   par[i] = 1;
                               }
                       }
                
                
                TF1 *fitV2_2PtBinnedAlt = new TF1("fitV2_2PtBinnedAlt",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                fitV2_2PtBinnedAlt->SetNpx(500);
                fitV2_2PtBinnedAlt->SetLineWidth(4);
                fitV2_2PtBinnedAlt->SetLineColor(kMagenta);
                fitV2_2PtBinnedAlt->SetParameters(par);
                TF1 *backFcnV2_2Alt = new TF1("backFcnV2_2Alt",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                backFcnV2_2Alt->SetLineColor(kRed);
                TF1 *signalFcnJPsiV2_2Alt = new TF1("signalFcnJPsiV2_2Alt",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                signalFcnJPsiV2_2Alt->SetLineColor(kBlue);
                signalFcnJPsiV2_2Alt->SetNpx(500);
                
                 for(int i=0; i<numParameters; i++){
                     fitV2_2PtBinnedAlt->FixParameter(i,par[i]);
                 }
                
                for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                    string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "V2_2");
                    const char *titular = title.c_str();
                    fitV2_2PtBinnedAlt->SetParName(index,titular);
                }
                
                
                
                    //Fit of V2
                        
                           double minParams[numParameters+numParametersBkgV2+1];
                           double parErrors[numParameters+numParametersBkgV2+1];
                            Double_t para[numParameters+numParametersBkgV2+1];
                
                c16PtBinnedFin->cd(ptbin+1);
                
                TVirtualFitter::Fitter(V22ATLASPtBinned[ptbin])->SetMaxIterations(10000);
                TVirtualFitter::Fitter(V22ATLASPtBinned[ptbin])->SetPrecision();
                   
                          res = V22ATLASPtBinned[ptbin]->Fit("fitV2_2PtBinnedAlt","SBMERI+","ep");
                          
                          gStyle->SetOptFit(1011);
                          TPaveStats *st = (TPaveStats*)V22ATLASPtBinned[ptbin]->FindObject("stats");
                          st->SetX1NDC(0.75); //new x start position
                          st->SetY1NDC(0.75); //new x end position
                          fitV2_2PtBinnedAlt->GetParameters(para);
                          backFcnV2_2Alt->SetParameters(para);
                          signalFcnJPsiV2_2Alt->SetParameters(para);
        
                              backFcnV2_2Alt->SetParameters(para);
                              signalFcnJPsiV2_2Alt->SetParameters(para);
                      
                      fitV2_2PtBinnedAlt->Draw("same");
                           
                    //  signalFcnJPsiV2_2Alt->Draw("same");
                      backFcnV2_2Alt->Draw("same");
                          gStyle->SetOptFit(1011);
                        // draw the legend
                        TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
                        legend->SetFillColorAlpha(kWhite, 0.);
                        legend->SetBorderSize(0);
                          legend->SetTextFont(42);
                        legend->SetTextSize(0.03);
                      Char_t message[80];
                      sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitV2_2PtBinnedAlt->GetChisquare(),fitV2_2PtBinnedAlt->GetNDF());
                       legend->AddEntry(fitV2_2PtBinnedAlt,message);
                          sprintf(message,"V2,2 = %.4e / %.4e",para[numParameters],fitV2_2PtBinnedAlt->GetParError(numParameters));
                          legend->AddEntry(fitV2_2PtBinnedAlt,message);
                          if(res->CovMatrixStatus() == 3){
                                // sprintf(message,"The fit is a success");
                          }
                          else{
                                 sprintf(message,"The fit is a failure");
                              legend->AddEntry(fitV2_2PtBinnedAlt,message);
                          }
                          
                     // legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
                      legend->AddEntry(backFcnV2_2Alt,"Background");
                        legend->AddEntry(V22ATLASPtBinned[ptbin],"Data","lpe");
                        legend->Draw();
                          
                       V2_ATLASExt2[ptbin] = para[numParameters];
                          errV2_ATLASExt2[ptbin] = fitV2_2PtBinnedAlt->GetParError(numParameters);
                        
                }
            
            sprintf(CanvasName,"%s/V2_ATLASExt2_PtBinned.pdf",CanvasFolderName);
            c16PtBinnedFin->SaveAs(CanvasName);
            
            }
    }
    
    
    //Ext 2 Pt Binned no ZYAM
    
     TCanvas *cinvmassPtBinned_noZYAM = new TCanvas("cinvmassPtBinned_noZYAM","Fitting All Phi - Different Centralities - Pt Binned",10,10,700,500);
        cinvmassPtBinned_noZYAM->SetTitle("Inv Mass Fits - Pt Binned - no ZYAM (useless plot)");
        // ROOT THAT Inv mass fits
        cinvmassPtBinned_noZYAM->Divide(NbPtBins,3);
        
    if(PtBinned){
        for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    
        sprintf(histoname,Form("bin%d_0",ptbin+1)); // ZZZZZZZ
        res = FittingAllInvMassBin(histoname, cinvmassPtBinned_noZYAM, ptbin, SignalF, BackgroundF, minmass, maxmass, ratsigma);
        double par[numParameters+numParametersBkgV2+1];
        
        // Parameter setting (either from mass fit either start values)
        for(int i=0; i <numParameters+numParametersBkgV2+1; i++){
                par[i] = res->Parameter(i);
                if(i>=numParameters){
                    par[i] = 1;
                }
        }
        
        //Definition and design of function on V2 plot
        c16PtBinned_noZYAM->cd(ptbin+1);
        TF1 *fitV2_2PtBinned = new TF1("fitV2_2PtBinned",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        fitV2_2PtBinned->SetNpx(500);
        fitV2_2PtBinned->SetLineWidth(4);
        fitV2_2PtBinned->SetLineColor(kMagenta);
        fitV2_2PtBinned->SetParameters(par);
        TF1 *backFcnV2_2 = new TF1("backFcnV2_2",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        backFcnV2_2->SetLineColor(kRed);
        TF1 *signalFcnJPsiV2_2 = new TF1("signalFcnJPsiV2_2",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        signalFcnJPsiV2_2->SetLineColor(kBlue);
        signalFcnJPsiV2_2->SetNpx(500);
        
           TVirtualFitter::Fitter(V2JPsiTklPtBinned[ptbin])->SetMaxIterations(10000);
           TVirtualFitter::Fitter(V2JPsiTklPtBinned[ptbin])->SetPrecision();
            for(int i=0; i<numParameters; i++){
                fitV2_2PtBinned->FixParameter(i,par[i]);
            }
            
            for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "V2_2");
                const char *titular = title.c_str();
                fitV2_2PtBinned->SetParName(index,titular);
            }
            
//           fitV2_2PtBinned->SetParName(0,"Norm_{J/#psi}");
//           fitV2_2PtBinned->SetParName(1,"M_{J/#psi}");
//           fitV2_2PtBinned->SetParName(2,"Sigma_{J/#psi}");
//           fitV2_2PtBinned->SetParName(3,"a_{1}");
//           fitV2_2PtBinned->SetParName(4,"n_{1}");
//           fitV2_2PtBinned->SetParName(5,"a_{2}");
//           fitV2_2PtBinned->SetParName(6,"n_{2}");
//           fitV2_2PtBinned->SetParName(7,"Norm_{#Psi(2S)}");
//           fitV2_2PtBinned->SetParName(8,"Norm_{TailLowM}");
//           fitV2_2PtBinned->SetParName(9,"Exp_{TailLowM}");
//           fitV2_2PtBinned->SetParName(10,"Norm_{TailHighM}");
//           fitV2_2PtBinned->SetParName(11,"Exp_{TailHighM}");
//            fitV2_2PtBinned->SetParName(12,"V2_2 J/#psi");
//            fitV2_2PtBinned->SetParName(13,"V2_2 Bkg M2");
//        fitV2_2PtBinned->SetParName(14,"V2_2 Bkg M1");
//        fitV2_2PtBinned->SetParName(15,"V2_2 Bkg M0");

        
           double minParams[numParameters+numParametersBkgV2+1];
           double parErrors[numParameters+numParametersBkgV2+1];
        
        //Fit of V2
            res = V2JPsiTklPtBinned_noZYAM[ptbin]->Fit("fitV2_2PtBinned","SBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)V2JPsiTklPtBinned_noZYAM[ptbin]->FindObject("stats");
            st->SetX1NDC(0.75); //new x start position
            st->SetY1NDC(0.75); //new x end position
            Double_t para[numParameters+numParametersBkgV2+1];
            fitV2_2PtBinned->GetParameters(para);
            backFcnV2_2->SetParameters(para);
            signalFcnJPsiV2_2->SetParameters(para);
        
        
        fitV2_2PtBinned->Draw("same");
      //  signalFcnJPsiV2_2->Draw("same");
        backFcnV2_2->Draw("same");
          // draw the legend
          TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
          legend->SetFillColorAlpha(kWhite, 0.);
          legend->SetBorderSize(0);
            legend->SetTextFont(42);
          legend->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitV2_2PtBinned->GetChisquare(),fitV2_2PtBinned->GetNDF());
         legend->AddEntry(fitV2_2PtBinned,message);

            if(res->CovMatrixStatus() == 3){
                 //  sprintf(message,"The fit is a success");
            }
            else{
                   sprintf(message,"The fit is a failure");
                legend->AddEntry(fitV2_2PtBinned,message);
            }
       // legend->AddEntry(signalFcnJPsiV2_2,"JPsi signal");
        legend->AddEntry(backFcnV2_2,"Background");
          legend->AddEntry(V2JPsiTklPtBinned[ptbin],"Data","lpe");
          legend->Draw();
            
            V2_Ext2_noZYAM[ptbin] = para[numParameters];
            errV2_Ext2_noZYAM[ptbin] = fitV2_2PtBinned->GetParError(numParameters);
        }
        
        sprintf(CanvasName,"%s/V2_Ext2_noZYAM_PtBinned.pdf",CanvasFolderName);
        c16PtBinned_noZYAM->SaveAs(CanvasName);
    }
    
    
    // Ext 3 PtBinned
    
    
    double b0_JPsiPtBinned[NbPtBins];
    memset( b0_JPsiPtBinned, 0, NbPtBins*sizeof(double));
        double c0_JPsiPtBinned[NbPtBins];
    memset( c0_JPsiPtBinned, 0, NbPtBins*sizeof(double));
        double c2_JPsiPtBinned[NbPtBins];
    memset( c2_JPsiPtBinned, 0, NbPtBins*sizeof(double));
        double eb0_JPsiPtBinned[NbPtBins];
    memset( eb0_JPsiPtBinned, 0, NbPtBins*sizeof(double));
        double ec0_JPsiPtBinned[NbPtBins];
    memset( ec0_JPsiPtBinned, 0, NbPtBins*sizeof(double));
        double ec2_JPsiPtBinned[NbPtBins];
    memset( ec2_JPsiPtBinned, 0, NbPtBins*sizeof(double));
    
    TCanvas*c316PtBinned=new TCanvas();
    // ROOT THAT Extraction 2 plot
    c316PtBinned->SetTitle("Extraction method 3 PtBinned");
    c316PtBinned->Divide(NbPtBins,4);
    
    if(PtBinned){
           for(int ptbin=0;ptbin<NbPtBins;ptbin++){
       
           sprintf(histoname,Form("bin%d_0",ptbin+1)); // ZZZZZZZ
           res = FittingAllInvMassBin(histoname, cinvmassPtBinned, ptbin, SignalF, BackgroundF, minmass, maxmass, ratsigma);
           double par[numParameters+numParametersBkgV2+1];
           
           // Parameter setting (either from mass fit either start values)
           for(int i=0; i <numParameters+numParametersBkgV2+1; i++){
                   par[i] = res->Parameter(i);
                   if(i>=numParameters){
                       par[i] = 1;
                   }
           }
        
    {
             c316PtBinned->cd(1+ptbin);
             TF1 *fitb0_3 = new TF1("fitb0_3",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             fitb0_3->SetNpx(500);
             fitb0_3->SetLineWidth(4);
             fitb0_3->SetLineColor(kMagenta);
             fitb0_3->SetParameters(par);
             TF1 *backFcnb0_3 = new TF1("backFcnb0_3",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             backFcnb0_3->SetLineColor(kRed);
             TF1 *signalFcnJPsib0_3 = new TF1("signalFcnJPsib0_3",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             signalFcnJPsib0_3->SetLineColor(kBlue);
             signalFcnJPsib0_3->SetNpx(500);
             
                TVirtualFitter::Fitter(baselines0PtBinned[ptbin])->SetMaxIterations(10000);
                TVirtualFitter::Fitter(baselines0PtBinned[ptbin])->SetPrecision();
                 for(int i=0; i<numParameters; i++){
                     fitb0_3->FixParameter(i,par[i]);
                 }
        
        for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
            string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "B0");
            const char *titular = title.c_str();
            fitb0_3->SetParName(index,titular);
        }
        
//                fitb0_3->SetParName(0,"Norm_{J/#psi}");
//                fitb0_3->SetParName(1,"M_{J/#psi}");
//                fitb0_3->SetParName(2,"Sigma_{J/#psi}");
//                fitb0_3->SetParName(3,"a_{1}");
//                fitb0_3->SetParName(4,"n_{1}");
//                fitb0_3->SetParName(5,"a_{2}");
//                fitb0_3->SetParName(6,"n_{2}");
//                fitb0_3->SetParName(7,"Norm_{#Psi(2S)}");
//                fitb0_3->SetParName(8,"Norm_{TailLowM}");
//                fitb0_3->SetParName(9,"Exp_{TailLowM}");
//                fitb0_3->SetParName(10,"Norm_{TailHighM}");
//                fitb0_3->SetParName(11,"Exp_{TailHighM}");
//                 fitb0_3->SetParName(12,"B0 J/#psi");
//                 fitb0_3->SetParName(13,"B0 Bkg M2");
//             fitb0_3->SetParName(14,"B0 Bkg M1");
//             fitb0_3->SetParName(15,"B0 Nkg M0");
             
             //Fit of b0
                 res = baselines0PtBinned[ptbin]->Fit("fitb0_3","SBMERI+","ep");
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)baselines0PtBinned[ptbin]->FindObject("stats");
                st->SetX1NDC(0.75); //new x start position
                st->SetY1NDC(0.75); //new x end position
                 Double_t para3[numParameters+numParametersBkgV2+1];
                 fitb0_3->GetParameters(para3);
                 backFcnb0_3->SetParameters(para3);
                 signalFcnJPsib0_3->SetParameters(para3);
             
            b0_JPsiPtBinned[ptbin] = para3[numParameters];
            eb0_JPsiPtBinned[ptbin] = fitb0_3->GetParError(numParameters);
             
             fitb0_3->Draw("same");
           //  signalFcnJPsib0_2->Draw("same");
             backFcnb0_3->Draw("same");
               // draw the legend
               TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
               legend->SetFillColorAlpha(kWhite, 0.);
               legend->SetBorderSize(0);
                 legend->SetTextFont(42);
               legend->SetTextSize(0.03);
             Char_t message[80];
             sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitb0_3->GetChisquare(),fitb0_3->GetNDF());
              legend->AddEntry(fitb0_3,message);

                 if(res->CovMatrixStatus() == 3){
                      //  sprintf(message,"The fit is a success");
                 }
                 else{
                        sprintf(message,"The fit is a failure");
                    legend->AddEntry(fitb0_3,message);
                 }
            // legend->AddEntry(signalFcnJPsib0_2,"JPsi signal");
             legend->AddEntry(backFcnb0_3,"Background");
               legend->AddEntry(baselines0PtBinned[ptbin],"Data","lpe");
               legend->Draw();
        
        
        
        
    }
        
        {
                 c316PtBinned->cd(1+ptbin + 1*NbPtBins);
                 TF1 *fitc0_3 = new TF1("fitc0_3",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 fitc0_3->SetNpx(500);
                 fitc0_3->SetLineWidth(4);
                 fitc0_3->SetLineColor(kMagenta);
                 fitc0_3->SetParameters(par);
                 TF1 *backFcnc0_3 = new TF1("backFcnc0_3",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 backFcnc0_3->SetLineColor(kRed);
                 TF1 *signalFcnJPsic0_3 = new TF1("signalFcnJPsic0_3",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 signalFcnJPsic0_3->SetLineColor(kBlue);
                 signalFcnJPsic0_3->SetNpx(500);
                 
                    TVirtualFitter::Fitter(coefficients0PtBinned[ptbin])->SetMaxIterations(10000);
                    TVirtualFitter::Fitter(coefficients0PtBinned[ptbin])->SetPrecision();
                     for(int i=0; i<numParameters; i++){
                         fitc0_3->FixParameter(i,par[i]);
                     }
            
            for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "C0");
                const char *titular = title.c_str();
                fitc0_3->SetParName(index,titular);
            }
//                    fitc0_3->SetParName(0,"Norm_{J/#psi}");
//                    fitc0_3->SetParName(1,"M_{J/#psi}");
//                    fitc0_3->SetParName(2,"Sigma_{J/#psi}");
//                    fitc0_3->SetParName(3,"a_{1}");
//                    fitc0_3->SetParName(4,"n_{1}");
//                    fitc0_3->SetParName(5,"a_{2}");
//                    fitc0_3->SetParName(6,"n_{2}");
//                    fitc0_3->SetParName(7,"Norm_{#Psi(2S)}");
//                    fitc0_3->SetParName(8,"Norm_{TailLowM}");
//                    fitc0_3->SetParName(9,"Exp_{TailLowM}");
//                    fitc0_3->SetParName(10,"Norm_{TailHighM}");
//                    fitc0_3->SetParName(11,"Exp_{TailHighM}");
//                     fitc0_3->SetParName(12,"c0 J/#psi");
//                     fitc0_3->SetParName(13,"c0 Bkg M2");
//                 fitc0_3->SetParName(14,"c0 Bkg M1");
//                 fitc0_3->SetParName(15,"c0 Nkg M0");
        
                 
                 //Fit of c0
                     res = coefficients0PtBinned[ptbin]->Fit("fitc0_3","SBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)coefficients0PtBinned[ptbin]->FindObject("stats");
                    st->SetX1NDC(0.75); //new x start position
                    st->SetY1NDC(0.75); //new x end position
                     Double_t para3[numParameters+numParametersBkgV2+1];
                     fitc0_3->GetParameters(para3);
                     backFcnc0_3->SetParameters(para3);
                     signalFcnJPsic0_3->SetParameters(para3);
                 
                c0_JPsiPtBinned[ptbin] = para3[numParameters];
                ec0_JPsiPtBinned[ptbin] = fitc0_3->GetParError(numParameters);
                 
                 fitc0_3->Draw("same");
               //  signalFcnJPsic0_2->Draw("same");
                 backFcnc0_3->Draw("same");
                   // draw the legend
                   TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
                   legend->SetFillColorAlpha(kWhite, 0.);
                   legend->SetBorderSize(0);
                     legend->SetTextFont(42);
                   legend->SetTextSize(0.03);
                 Char_t message[80];
                 sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitc0_3->GetChisquare(),fitc0_3->GetNDF());
                  legend->AddEntry(fitc0_3,message);

                     if(res->CovMatrixStatus() == 3){
                        //    sprintf(message,"The fit is a success");
                     }
                     else{
                            sprintf(message,"The fit is a failure");
                         legend->AddEntry(fitc0_3,message);
                     }
                // legend->AddEntry(signalFcnJPsic0_2,"JPsi signal");
                 legend->AddEntry(backFcnc0_3,"Background");
                   legend->AddEntry(coefficients0PtBinned[ptbin],"Data","lpe");
                   legend->Draw();
            
            
            
            
        }
        
        {
                 c316PtBinned->cd(1+ptbin + 2*NbPtBins);
                 TF1 *fitc2_3 = new TF1("fitc2_3",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 fitc2_3->SetNpx(500);
                 fitc2_3->SetLineWidth(4);
                 fitc2_3->SetLineColor(kMagenta);
                 fitc2_3->SetParameters(par);
                 TF1 *backFcnc2_3 = new TF1("backFcnc2_3",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 backFcnc2_3->SetLineColor(kRed);
                 TF1 *signalFcnJPsic2_3 = new TF1("signalFcnJPsic2_3",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 signalFcnJPsic2_3->SetLineColor(kBlue);
                 signalFcnJPsic2_3->SetNpx(500);
                 
                    TVirtualFitter::Fitter(coefficients2PtBinned[ptbin])->SetMaxIterations(10000);
                    TVirtualFitter::Fitter(coefficients2PtBinned[ptbin])->SetPrecision();
                     for(int i=0; i<numParameters; i++){
                         fitc2_3->FixParameter(i,par[i]);
                     }
            
            for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "C2");
                const char *titular = title.c_str();
                fitc2_3->SetParName(index,titular);
            }
//                    fitc2_3->SetParName(0,"Norm_{J/#psi}");
//                    fitc2_3->SetParName(1,"M_{J/#psi}");
//                    fitc2_3->SetParName(2,"Sigma_{J/#psi}");
//                    fitc2_3->SetParName(3,"a_{1}");
//                    fitc2_3->SetParName(4,"n_{1}");
//                    fitc2_3->SetParName(5,"a_{2}");
//                    fitc2_3->SetParName(6,"n_{2}");
//                    fitc2_3->SetParName(7,"Norm_{#Psi(2S)}");
//                    fitc2_3->SetParName(8,"Norm_{TailLowM}");
//                    fitc2_3->SetParName(9,"Exp_{TailLowM}");
//                    fitc2_3->SetParName(10,"Norm_{TailHighM}");
//                    fitc2_3->SetParName(11,"Exp_{TailHighM}");
//                     fitc2_3->SetParName(12,"c2 J/#psi");
//                     fitc2_3->SetParName(13,"c2 Bkg M2");
//                 fitc2_3->SetParName(14,"c2 Bkg M1");
//                 fitc2_3->SetParName(15,"c2 Bkg M0");
                 
                 //Fit of c2
                     res = coefficients2PtBinned[ptbin]->Fit("fitc2_3","SBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)coefficients2PtBinned[ptbin]->FindObject("stats");
                    st->SetX1NDC(0.75); //new x start position
                    st->SetY1NDC(0.75); //new x end position
                     Double_t para3[numParameters+numParametersBkgV2+1];
                     fitc2_3->GetParameters(para3);
                     backFcnc2_3->SetParameters(para3);
                     signalFcnJPsic2_3->SetParameters(para3);
                 
                c2_JPsiPtBinned[ptbin] = para3[numParameters];
                ec2_JPsiPtBinned[ptbin] = fitc2_3->GetParError(numParameters);
                 
                 fitc2_3->Draw("same");
               //  signalFcnJPsic2_2->Draw("same");
                 backFcnc2_3->Draw("same");
                   // draw the legend
                   TLegend *legend=new TLegend(0.12,0.80,0.60,0.90);
                   legend->SetFillColorAlpha(kWhite, 0.);
                   legend->SetBorderSize(0);
                     legend->SetTextFont(42);
                   legend->SetTextSize(0.03);
                 Char_t message[80];
                 sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitc2_3->GetChisquare(),fitc2_3->GetNDF());
                  legend->AddEntry(fitc2_3,message);

                     if(res->CovMatrixStatus() == 3){
                           // sprintf(message,"The fit is a success");
                     }
                     else{
                            sprintf(message,"The fit is a failure");
                            legend->AddEntry(fitc2_3,message);
                     }
                // legend->AddEntry(signalFcnJPsic2_2,"JPsi signal");
                 legend->AddEntry(backFcnc2_3,"Background");
                   legend->AddEntry(coefficients2PtBinned[ptbin],"Data","lpe");
                   legend->Draw();
            
            
            
            
        }
        
        {
            c316PtBinned->cd(1+ptbin + 3*NbPtBins);
            
            V2_Ext3[ptbin] = c2_JPsiPtBinned[ptbin]/(c0_JPsiPtBinned[ptbin] + b0_JPsiPtBinned[ptbin]);
           // errV2_Ext3[ptbin] = (c2_JPsiPtBinned[ptbin]/(c0_JPsiPtBinned[ptbin] + b0_JPsiPtBinned[ptbin])*sqrt(pow(ec2_JPsiPtBinned[ptbin]/c2_JPsiPtBinned[ptbin],2)+pow(ec0_JPsiPtBinned[ptbin]/c0_JPsiPtBinned[ptbin],2)+pow(eb0_JPsiPtBinned[ptbin]/b0_JPsiPtBinned[ptbin],2)));
            
            errV2_Ext3[ptbin] = (c2_JPsiPtBinned[ptbin]/(c0_JPsiPtBinned[ptbin] + b0_JPsiPtBinned[ptbin])*sqrt(pow(ec2_JPsiPtBinned[ptbin]/c2_JPsiPtBinned[ptbin],2)+pow(sqrt(pow(ec0_JPsiPtBinned[ptbin],2)+pow(eb0_JPsiPtBinned[ptbin],2))/(c0_JPsiPtBinned[ptbin] + b0_JPsiPtBinned[ptbin]),2)));
            
            V2_Ext3_noZYAM[ptbin] = c2_JPsiPtBinned[ptbin]/(c0_JPsiPtBinned[ptbin]);
            errV2_Ext3_noZYAM[ptbin] = (c2_JPsiPtBinned[ptbin]/(c0_JPsiPtBinned[ptbin])*sqrt(pow(ec2_JPsiPtBinned[ptbin]/c2_JPsiPtBinned[ptbin],2)+pow(ec0_JPsiPtBinned[ptbin]/c0_JPsiPtBinned[ptbin],2)));
            
             Char_t message[80];
              sprintf(message,"V2_3 J/#psi-tkl: #frac{%.4f}{%.4f + %.4f} = %.5f +- %.5f",c2_JPsiPtBinned[ptbin],c0_JPsiPtBinned[ptbin], b0_JPsiPtBinned[ptbin],V2_Ext3[ptbin],errV2_Ext3[ptbin]);
            TPaveText *pave = new TPaveText(0.2,0.2,0.8,0.8,"brNDC");
            pave->AddText(message);
            pave->SetTextFont(42);
            pave->SetTextSize(0.05);
            pave->SetBorderSize(0);
            pave->SetFillStyle(0);
               pave->Draw();
            
        }
           }
        
        sprintf(CanvasName,"%s/V2_Ext3_PtBinned.pdf",CanvasFolderName);
        c316PtBinned->SaveAs(CanvasName);
    }
    
    
    
    //Combining or not fit for Central-Periph method Ext1
    
    
    int niterations = 1;
    cinvmass->cd();
    sprintf(histoname,"InvMass_Central");
    res = FittingAllInvMassBin(histoname, cinvmass, 1, SignalF, BackgroundF, minmass, maxmass, ratsigma);

    double parerr[numParameters];
    TCanvas* c10_fit = new TCanvas;
       c10_fit->Divide(6,4);
    for(int j=0; j<numParameters; j++){
        parerr[j] = res->ParError(j);
    }
    for(int j=0; j <numParameters+numParametersBkgV2+1; j++){
        par[j] = res->Parameter(j);
        if(j>=numParameters){
            par[j] = 3;
        }
    }

    double par12[niterations];
    double parerr12[niterations];
    
    Char_t fileLogName[200] = "Log_FitTraining.txt";
    std::ofstream fileLog(fileLogName, std::ofstream::out);
    
    for(int i=0; i<NbinsDeltaPhi; i++){ //Change
        fileLog << "Bin " << i << " in DeltaPhi Central" <<endl;
        
        for(int t=0; t<niterations; t++){
            fileLog << "Interation " << t <<endl;
        TVirtualPad* subpad = c10_fit->cd(i+1);
            c10_fit->cd(i+1);
            if(t==niterations-1){
                subpad->Clear();
            }
        TF1 *fitY_1Central = new TF1("fitY_1Central",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
          fitY_1Central->SetNpx(500);
          fitY_1Central->SetLineWidth(4);
          fitY_1Central->SetLineColor(kMagenta);
          fitY_1Central->SetParameters(par);
           TVirtualFitter::Fitter(YieldWrtMass_Central[i])->SetMaxIterations(10000);
           TVirtualFitter::Fitter(YieldWrtMass_Central[i])->SetPrecision();
        //  histo->Fit("fitFcn","0");
          // second try: set start values for some parameters
        for(int k=0; k<numParameters; k++){
                fitY_1Central->FixParameter(k,par[k]);
        }
            
            fitY_1Central->SetParLimits(numParameters,0,10000);
            
            for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "Y_1");
                const char *titular = title.c_str();
                fitY_1Central->SetParName(index,titular);
            }

//           fitY_1Central->SetParName(0,"Norm_{J/#psi}");
//           fitY_1Central->SetParName(1,"M_{J/#psi}");
//           fitY_1Central->SetParName(2,"Sigma_{J/#psi}");
//           fitY_1Central->SetParName(3,"a_{1}");
//           fitY_1Central->SetParName(4,"n_{1}");
//           fitY_1Central->SetParName(5,"a_{2}");
//           fitY_1Central->SetParName(6,"n_{2}");
//           fitY_1Central->SetParName(7,"Norm_{#Psi(2S)}");
//           fitY_1Central->SetParName(8,"Norm_{TailLowM}");
//           fitY_1Central->SetParName(9,"Exp_{TailLowM}");
//           fitY_1Central->SetParName(10,"Norm_{TailHighM}");
//           fitY_1Central->SetParName(11,"Exp_{TailHighM}");
//            fitY_1Central->SetParName(12,"Y_1 J/#psi");
//            fitY_1Central->SetParName(13,"Y_1 Bkg M2");
//        fitY_1Central->SetParName(14,"Y_1 Bkg M1");
//        fitY_1Central->SetParName(15,"Y_1 Bkg M0");

          rescent = YieldWrtMass_Central[i]->Fit("fitY_1Central","SMERIQ+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)YieldWrtMass_Central[i]->FindObject("stats");
            st->SetX1NDC(0.75); //new x start position
            st->SetY1NDC(0.75); //new x end position
            rescent->Print("V");
            fileLog << "Chi2: " << rescent->Chi2() <<endl;
            fileLog << "Status: " << int(rescent->Status()) <<endl;
            fileLog << "Prob: " << rescent->Prob() <<endl;
            fileLog << "IsValid: " << rescent->IsValid() <<endl;
            fileLog << "HasMinosError(numParameters): " << rescent->HasMinosError(numParameters) <<endl;
            
          Double_t param[numParameters+numParametersBkgV2+1];
            Double_t paramerrs[numParameters+numParametersBkgV2+1];
          fitY_1Central->GetParameters(param);
            for(int i=0; i<numParameters+numParametersBkgV2+1; i++){
                paramerrs[i] = fitY_1Central->GetParError(i);
            }

            if(rescent->CovMatrixStatus() == 3){
            par12[t] = param[numParameters];
            parerr12[t] = paramerrs[numParameters];
                YieldValue->Fill(param[numParameters]);
                YieldError->Fill(paramerrs[numParameters]);
            }
            else{
                par12[t]=parerr12[t]=0;
            }

            std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
            std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;
            
            fileLog << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
            fileLog << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

        TF1 *backFcnY_1Central = new TF1("backFcnY_1Central",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
        backFcnY_1Central->SetLineColor(kRed);
        TF1 *signalFcnJPsiY_1Central = new TF1("signalFcnJPsiY_1Central",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
          // writes the fit results into the par array

            if(t==niterations-1){
                YieldWrtMass_Central[i]->GetListOfFunctions()->Clear();
                YieldWrtMass_Central[i]->GetListOfFunctions()->Add(fitY_1Central);
            }
        signalFcnJPsiY_1Central->SetLineColor(kBlue);
        signalFcnJPsiY_1Central->SetNpx(500);
        backFcnY_1Central->SetParameters(param);
        signalFcnJPsiY_1Central->SetParameters(param);
       // signalFcnJPsiY_1Central->Draw("same");
        backFcnY_1Central->Draw("same");

//                if(i==2){
//                    baseline_central = param[12];
//                    errbaseline_central = fitY_1Central->GetParError(12);
//                }
//                if(i==3){
//                    baseline_central += param[12];
//                    baseline_central /= 2;
//                    errbaseline_central = sqrt(pow(errbaseline_central,2)+pow(fitY_1Central->GetParError(12),2));
//                }
                
                if(param[numParameters]<baseline_central){ //En central baseline = minimum
                    baseline_central = param[numParameters];
                    errbaseline_central = fitY_1Central->GetParError(numParameters);
                }
                
                Yields_Central_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[numParameters]);
                Yields_Central_1->SetBinError(i+1,fitY_1Central->GetParError(numParameters));
            

//          // writes the fit results into the par array
//           gStyle->SetOptFit(1011);
          // draw the legend
          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
          legend->SetFillColorAlpha(kWhite, 0.);
          legend->SetBorderSize(0);
            legend->SetTextFont(42);
          legend->SetTextSize(0.03);
        Char_t message[80];
        sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitY_1Central->GetChisquare(),fitY_1Central->GetNDF());
         legend->AddEntry(fitY_1Central,message);
        if(rescent->CovMatrixStatus() == 3){
                  // sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
                   legend->AddEntry(fitY_1Central,message);
               }
      //  legend->AddEntry(signalFcnJPsiY_1Central,"JPsi signal");
        legend->AddEntry(backFcnY_1Central,"Background");
          legend->AddEntry(YieldWrtMass_Central[i],"Data","lpe");
          legend->Draw();
        }
        
//        TCanvas* cyo = new TCanvas;
//        cyo->Divide(2,1);
//        cyo->cd(1);
//        YieldValue->Draw();
//        cyo->cd(2);
//        YieldError->Draw();
        
    }
    
    for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
        Baseline_Central_1->SetBinContent(phi_idx+1, baseline_central);
        Baseline_Central_1->SetBinError(phi_idx+1, errbaseline_central);
    }
    
    Yields_Central_1_MinusBaseline->Add(Yields_Central_1,Baseline_Central_1,1,-1);
    

        
    
    cinvmass->cd();
    sprintf(histoname,"InvMass_Periph");
    res = FittingAllInvMassBin(histoname, cinvmass, 2, SignalF, BackgroundF, minmass, maxmass, ratsigma);

       for(int j=0; j<numParameters; j++){
           parerr[j] = res->ParError(j);
       }
       for(int j=0; j <numParameters+numParametersBkgV2+1; j++){
           par[j] = res->Parameter(j);
           if(j>=numParameters){
               par[j] = 3;
           }
       }

//    double baseline = 1;
//    double errbaseline = 1;

       for(int i=0; i<NbinsDeltaPhi; i++){ //Change
           fileLog << "Bin " << i << " in DeltaPhi Periph" <<endl;
           for(int t=0; t<niterations; t++){
               fileLog << "Interation " << t <<endl;
               TVirtualPad* subpad = c10_fit->cd(NbinsDeltaPhi+i+1); //Change
               c10_fit->cd(NbinsDeltaPhi+i+1); //Change
               if(t==niterations-1){
                   subpad->Clear();
               }
           TF1 *fitY_1Periph = new TF1("fitY_1Periph",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
             fitY_1Periph->SetNpx(500);
             fitY_1Periph->SetLineWidth(4);
             fitY_1Periph->SetLineColor(kMagenta);
             fitY_1Periph->SetParameters(par);
              TVirtualFitter::Fitter(YieldWrtMass_Periph[i])->SetMaxIterations(10000);
              TVirtualFitter::Fitter(YieldWrtMass_Periph[i])->SetPrecision();
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters
            for(int k=0; k<numParameters; k++){
                    fitY_1Periph->FixParameter(k,par[k]);
            }
              fitY_1Periph->SetParLimits(numParameters,0,10000);

              fitY_1Periph->SetParName(0,"Norm_{J/#psi}");
              fitY_1Periph->SetParName(1,"M_{J/#psi}");
              fitY_1Periph->SetParName(2,"Sigma_{J/#psi}");
              fitY_1Periph->SetParName(3,"a_{1}");
              fitY_1Periph->SetParName(4,"n_{1}");
              fitY_1Periph->SetParName(5,"a_{2}");
              fitY_1Periph->SetParName(6,"n_{2}");
              fitY_1Periph->SetParName(7,"Norm_{#Psi(2S)}");
              fitY_1Periph->SetParName(8,"Norm_{TailLowM}");
              fitY_1Periph->SetParName(9,"Exp_{TailLowM}");
              fitY_1Periph->SetParName(10,"Norm_{TailHighM}");
              fitY_1Periph->SetParName(11,"Exp_{TailHighM}");
               fitY_1Periph->SetParName(12,"Y_1 J/#psi");
               fitY_1Periph->SetParName(13,"Y_1 Bkg M2");
           fitY_1Periph->SetParName(14,"Y_1 Bkg M1");
           fitY_1Periph->SetParName(15,"Y_1 Bkg M0");

             resperiph = YieldWrtMass_Periph[i]->Fit("fitY_1Periph","SMERIQ+","ep");
               gStyle->SetOptStat("n");
               gStyle->SetOptFit(1011);
               TPaveStats *st = (TPaveStats*)YieldWrtMass_Periph[i]->FindObject("stats");
               st->SetX1NDC(0.75); //new x start position
               st->SetY1NDC(0.75); //new x end position
               fileLog << "Chi2: " << resperiph->Chi2() <<endl;
               fileLog << "Status: " << int(resperiph->Status()) <<endl;
               fileLog << "Prob: " << resperiph->Prob() <<endl;
               fileLog << "IsValid: " << resperiph->IsValid() <<endl;
               fileLog << "HasMinosError(numParameters): " << resperiph->HasMinosError(numParameters) <<endl;
               Double_t param[numParameters+numParametersBkgV2+1];
               Double_t paramerrs[numParameters+numParametersBkgV2+1];
             fitY_1Periph->GetParameters(param);

               for(int i=0; i<numParameters+numParametersBkgV2+1; i++){
                               paramerrs[i] = fitY_1Periph->GetParError(i);
                           }

                           if(resperiph->CovMatrixStatus() == 3){
                           par12[t] = param[numParameters];
                           parerr12[t] = paramerrs[numParameters];
                           }
                           else{
                               par12[t]=parerr12[t]=0;
                           }
               std::cout << "COV status Periph t=" << t << " : " << resperiph->CovMatrixStatus()<<endl;
               std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;
               
               fileLog << "COV status Periph t=" << t << " : " << resperiph->CovMatrixStatus()<<endl;
               fileLog << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

           TF1 *backFcnY_1Periph = new TF1("backFcnY_1Periph",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
           backFcnY_1Periph->SetLineColor(kRed);
           TF1 *signalFcnJPsiY_1Periph = new TF1("signalFcnJPsiY_1Periph",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
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


                   if(i==0){
                       baseline_periph = param[numParameters];
                       errbaseline_periph = fitY_1Periph->GetParError(numParameters);
                   }
//                   if(i==3){
//                       baseline_periph += param[numParameters];
//                       baseline_periph /= 2;
//                       errbaseline_periph = sqrt(pow(errbaseline_periph,2)+pow(fitY_1Periph->GetParError(numParameters),2));
//                       errbaseline_periph /= 2;
//                   }
                 Yields_Periph_1->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[numParameters]);
                 Yields_Periph_1->SetBinError(i+1,fitY_1Periph->GetParError(numParameters));
               

             // writes the fit results into the par array
             // gStyle->SetOptFit(1011);
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetFillColorAlpha(kWhite, 0.);
             legend->SetBorderSize(0);
               legend->SetTextFont(42);
             legend->SetTextSize(0.04);
           Char_t message[80];
           sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitY_1Periph->GetChisquare(),fitY_1Periph->GetNDF());
            legend->AddEntry(fitY_1Periph,message);
           if(resperiph->CovMatrixStatus() == 3){
                     // sprintf(message,"The fit is a success");
                  }
                  else{
                      sprintf(message,"The fit is a failure");
                      legend->AddEntry(fitY_1Periph,message);
                  }
       //    legend->AddEntry(signalFcnJPsiY_1Periph,"JPsi signal");
           legend->AddEntry(backFcnY_1Periph,"Background");
             legend->AddEntry(YieldWrtMass_Periph[i],"Data","lpe");
             legend->Draw();
           }

       }
    
     sprintf(CanvasName,"%s/Idk.pdf",CanvasFolderName);
     c10_fit->SaveAs(CanvasName);
    
    for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
        Baseline_Periph_1->SetBinContent(phi_idx+1, baseline_periph);
        Baseline_Periph_1->SetBinError(phi_idx+1, errbaseline_periph);
    }
    
    Yields_Periph_1_MinusBaseline->Add(Yields_Periph_1,Baseline_Periph_1,1,-1);
    
    //Fill Baseline Periph and Y-B
    TCanvas* cYBcentral = new TCanvas;
    cYBcentral->SetTitle("Central Yield-Baseline");
    cYBcentral->Divide(1,3);
    cYBcentral->cd(1);
    Yields_Central_1->DrawCopy();
    cYBcentral->cd(2);
    Baseline_Central_1->DrawCopy();
    cYBcentral->cd(3);
    Yields_Central_1_MinusBaseline->DrawCopy();
    
    sprintf(CanvasName,"%s/CentralYield-Baseline.pdf",CanvasFolderName);
    cYBcentral->SaveAs(CanvasName);
    
    TCanvas* cYBperiph = new TCanvas;
    cYBperiph->SetTitle("Periph Yield-Baseline");
    cYBperiph->Divide(1,3);
    cYBperiph->cd(1);
    Yields_Periph_1->DrawCopy();
    cYBperiph->cd(2);
    Baseline_Periph_1->DrawCopy();
    cYBperiph->cd(3);
    Yields_Periph_1_MinusBaseline->DrawCopy();
    
    sprintf(CanvasName,"%s/PeriphYield-Baseline.pdf",CanvasFolderName);
    cYBperiph->SaveAs(CanvasName);
    
    TH1F *Yields_Central_1_Cvetan = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_Cvetan");
    TH1F *Yields_Central_1_CvetanMe = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_CvetanMe");
    TH1F *Yields_Central_1_ZYAM = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_ZYAM");
    TH1F *Yields_Central_1_PRLTemplate = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_PRLTemplate");
    TH1F *Yields_Central_1_PRLTemplate_PeriphZYAM = (TH1F*)Yields_Central_1->Clone("Yields_Central_1_PRLTemplate_PeriphZYAM");

    {
        TCanvas* c17 = new TCanvas;
        c17->Divide(1,3);
        //Creer canvas pour imprimer les plots Periph yield et Central yield wrt phi et leur difference
        c17->cd(1);
        Yields_Central_1->SetStats(kTRUE);
        Yields_Central_1->DrawCopy();
        c17->cd(2);
        Yields_Periph_1->SetStats(kTRUE);
        Yields_Periph_1->DrawCopy();
        c17->cd(3);
        Yields_Difference_1->Add(Yields_Central_1,Yields_Periph_1,1,-1);
        Yields_Difference_1->SetStats(kTRUE);
        Yields_Difference_1->DrawCopy();

        TF1 *fitFcnV2_2 = new TF1("fitFcnV2_2",FourierV2,0,TMath::Pi(),3);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
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
              //  gStyle->SetOptFit(1011);
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_2->SetParName(0,"a0");
              fitFcnV2_2->SetParName(1,"a1");
              fitFcnV2_2->SetParName(2,"a2");

             TFitResultPtr res = Yields_Difference_1->Fit("fitFcnV2_2","SBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)Yields_Difference_1->FindObject("stats");
            st->SetX1NDC(0.8); //new x start position
            st->SetY1NDC(0.8); //new x end position
            Double_t par[3];
            fitFcnV2_2->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_2->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetFillColorAlpha(kWhite, 0.);
             legend->SetBorderSize(0);
               legend->SetTextFont(42);
             legend->SetTextSize(0.03);
               Char_t message[80];
               sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
               legend->AddEntry(fitFcnV2_2,message);
        sprintf(message,"V2_1 J/#psi-tkl: #frac{%.4f}{%.4f + %.4f} = %.4f +- %.4f",par[2],par[0], baseline_periph,par[2]/(par[0] + baseline_periph),abs(par[2]/(par[0] + baseline_periph)*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(sqrt(pow(fitFcnV2_2->GetParError(0),2)+pow(errbaseline_periph,2))/(par[0] + baseline_periph),2))));
        legend->AddEntry(fitFcnV2_2,message);
        if(res->CovMatrixStatus() == 3){
                 //  sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
                   legend->AddEntry(fitFcnV2_2,message);
               }
             legend->AddEntry(Yields_Difference_1,"Data","lpe");
             legend->Draw();

        sprintf(CanvasName,"%s/V2_Ext1.pdf",CanvasFolderName);
        c17->SaveAs(CanvasName);
    }
    
    
    
    // Fit Cvetan
    
    TCanvas* cCvetan = new TCanvas;
       cCvetan->SetTitle("Cvetan fit");
       cCvetan->Divide(1,1);
       cCvetan->cd(1);
       Yields_Central_1->DrawCopy();
//       cCvetan->cd(2);
//       Baseline_Central_1->DrawCopy();
//       cCvetan->cd(3);
//       Yields_Central_1_MinusBaseline->DrawCopy();
    
    
    {
        cCvetan->cd(1);
             TF1 *fitFcnV2_Cvetan = new TF1("fitFcnV2_Cvetan",CvetanF,0,TMath::Pi(),4);//-TMath::Pi()/2,1.5*TMath::Pi(),4);
             fitFcnV2_Cvetan->SetNpx(500);
             fitFcnV2_Cvetan->SetLineWidth(4);
             fitFcnV2_Cvetan->SetLineColor(kOrange+3);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[4] = {1,0,0.01,1};
             fitFcnV2_Cvetan->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->SetPrecision();
                TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->SetFCN(ChisquareCvetanF);
           // gStyle->SetOptFit(1011);
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_Cvetan->SetParName(0,"V0");
                fitFcnV2_Cvetan->SetParName(1,"V1");
                fitFcnV2_Cvetan->SetParName(2,"V2");
              fitFcnV2_Cvetan->SetParName(3,"F");
        
        fitFcnV2_Cvetan->SetParLimits(3,0.1,50);
        
        fitFcnV2_Cvetan->FixParameter(0,1);
        fitFcnV2_Cvetan->FixParameter(1,0);
      //  fitFcnV2_Cvetan->FixParameter(3,1);
              

             TFitResultPtr res = Yields_Central_1_Cvetan->Fit("fitFcnV2_Cvetan","USBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)Yields_Central_1_Cvetan->FindObject("stats");
            st->SetX1NDC(0.8); //new x start position
            st->SetY1NDC(0.8); //new x end position
            Double_t par[4];
        
                double chi2, edm, errdef;
               int nvpar, nparx;
                TVirtualFitter::Fitter(Yields_Central_1_Cvetan)->GetStats(chi2,edm,errdef,nvpar,nparx);
                fitFcnV2_Cvetan->SetChisquare(chi2);
                int ndf = npfits-nvpar;
                fitFcnV2_Cvetan->SetNDF(ndf);
        
            fitFcnV2_Cvetan->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_Cvetan->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetFillColorAlpha(kWhite, 0.);
             legend->SetBorderSize(0);
               legend->SetTextFont(42);
             legend->SetTextSize(0.03);
               Char_t message[80];
               sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_Cvetan->GetChisquare(),fitFcnV2_Cvetan->GetNDF());
               legend->AddEntry(fitFcnV2_Cvetan,message);
        if(res->CovMatrixStatus() == 3){
                 //  sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
                   legend->AddEntry(fitFcnV2_Cvetan,message);
               }
             legend->AddEntry(Yields_Central_1_Cvetan,"Data","lpe");
             legend->Draw();

        sprintf(CanvasName,"%s/V2_Cvetan.pdf",CanvasFolderName);
        cCvetan->SaveAs(CanvasName);
    }
    
    // Fit Cvetan putting my constraints (F=1, V1 free, V0 not 1)
    
    TCanvas* cCvetanMe = new TCanvas;
       cCvetanMe->SetTitle("Cvetan fit - My constraints F=1, v1 free, v0 free");
       cCvetanMe->Divide(1,1);
       cCvetanMe->cd(1);
       Yields_Central_1->DrawCopy();
//       cCvetanMe->cd(2);
//       Baseline_Central_1->DrawCopy();
//       cCvetanMe->cd(3);
//       Yields_Central_1_MinusBaseline->DrawCopy();
    
    
    {
        cCvetanMe->cd(1);
             TF1 *fitFcnV2_CvetanMe = new TF1("fitFcnV2_CvetanMe",CvetanF,0,TMath::Pi(),4);//-TMath::Pi()/2,1.5*TMath::Pi(),4);
             fitFcnV2_CvetanMe->SetNpx(500);
             fitFcnV2_CvetanMe->SetLineWidth(4);
             fitFcnV2_CvetanMe->SetLineColor(kBlue);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[4] = {1,0,0.01,1};
             fitFcnV2_CvetanMe->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->SetPrecision();
            TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->SetFCN(ChisquareCvetanF);
           // gStyle->SetOptFit(1011);
           //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_CvetanMe->SetParName(0,"V0");
                fitFcnV2_CvetanMe->SetParName(1,"V1");
                fitFcnV2_CvetanMe->SetParName(2,"V2");
              fitFcnV2_CvetanMe->SetParName(3,"F");
        
        fitFcnV2_CvetanMe->SetParLimits(3,0.1,50);
        
     //   fitFcnV2_CvetanMe->FixParameter(0,1);
       // fitFcnV2_CvetanMe->FixParameter(1,0);
        fitFcnV2_CvetanMe->FixParameter(3,1);
              

             TFitResultPtr res = Yields_Central_1_CvetanMe->Fit("fitFcnV2_CvetanMe","USBMERI+","ep");
            gStyle->SetOptStat("n");
           gStyle->SetOptFit(1011);
           TPaveStats *st = (TPaveStats*)Yields_Central_1_CvetanMe->FindObject("stats");
           st->SetX1NDC(0.8); //new x start position
           st->SetY1NDC(0.8); //new x end position
        
        double chi2, edm, errdef;
        int nvpar, nparx;
         TVirtualFitter::Fitter(Yields_Central_1_CvetanMe)->GetStats(chi2,edm,errdef,nvpar,nparx);
         fitFcnV2_CvetanMe->SetChisquare(chi2);
         int ndf = npfits-nvpar;
         fitFcnV2_CvetanMe->SetNDF(ndf);
        
            Double_t par[4];
            fitFcnV2_CvetanMe->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_CvetanMe->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetFillColorAlpha(kWhite, 0.);
             legend->SetBorderSize(0);
               legend->SetTextFont(42);
             legend->SetTextSize(0.03);
               Char_t message[80];
               sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_CvetanMe->GetChisquare(),fitFcnV2_CvetanMe->GetNDF());
               legend->AddEntry(fitFcnV2_CvetanMe,message);
        if(res->CovMatrixStatus() == 3){
                 //  sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
                   legend->AddEntry(fitFcnV2_CvetanMe,message);
               }
             legend->AddEntry(Yields_Central_1_CvetanMe,"Data","lpe");
             legend->Draw();

    }
    
    
    // Fit ZYAM Yc = Bc + (Yp-Bp) + a0  + 2a1cos + 2a2cos
    
    
     TCanvas* cZYAM = new TCanvas;
        cZYAM->SetTitle("ZYAM fit");
        cZYAM->Divide(1,1);
        cZYAM->cd(1);
        Yields_Central_1->DrawCopy();
    //    cZYAM->cd(2);
    //    Baseline_Central_1->DrawCopy();
    //    cZYAM->cd(3);
    //    Yields_Central_1_MinusBaseline->DrawCopy();
    
    
    {
        cZYAM->cd(1);
             TF1 *fitFcnV2_ZYAM = new TF1("fitFcnV2_ZYAM",ZYAM,0,TMath::Pi(),3);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
             fitFcnV2_ZYAM->SetNpx(500);
             fitFcnV2_ZYAM->SetLineWidth(4);
             fitFcnV2_ZYAM->SetLineColor(kCyan-7);
             // first try without starting values for the parameters
             // This defaults to 1 for each param.
             // this results in an ok fit for the polynomial function
             // however the non-linear part (lorenzian) does not
             // respond well.
              Double_t params[3] = {1,1,1};
             fitFcnV2_ZYAM->SetParameters(params);
              TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->SetMaxIterations(10000);
              TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->SetPrecision();
            TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->SetFCN(ChisquareZYAM);
           // gStyle->SetOptFit(1011);
            //  histo->Fit("fitFcn","0");
             // second try: set start values for some parameters

              fitFcnV2_ZYAM->SetParName(0,"a0");
                fitFcnV2_ZYAM->SetParName(1,"a1");
                fitFcnV2_ZYAM->SetParName(2,"a2");
            
              

             TFitResultPtr res = Yields_Central_1_ZYAM->Fit("fitFcnV2_ZYAM","USBMERI+","ep");
            gStyle->SetOptStat("n");
            gStyle->SetOptFit(1011);
            TPaveStats *st = (TPaveStats*)Yields_Central_1_ZYAM->FindObject("stats");
            st->SetX1NDC(0.8); //new x start position
            st->SetY1NDC(0.8); //new x end position
        
            double chi2, edm, errdef;
            int nvpar, nparx;
             TVirtualFitter::Fitter(Yields_Central_1_ZYAM)->GetStats(chi2,edm,errdef,nvpar,nparx);
             fitFcnV2_ZYAM->SetChisquare(chi2);
             int ndf = npfits-nvpar;
             fitFcnV2_ZYAM->SetNDF(ndf);
        
            Double_t par[3];
            fitFcnV2_ZYAM->GetParameters(par);
             // improve the pictu
           //   std::cout << "integral error: " << integralerror << std::endl;
             fitFcnV2_ZYAM->Draw("same");
             // draw the legend
             TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
             legend->SetFillColorAlpha(kWhite, 0.);
             legend->SetBorderSize(0);
               legend->SetTextFont(42);
             legend->SetTextSize(0.03);
               Char_t message[80];
               sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_ZYAM->GetChisquare(),fitFcnV2_ZYAM->GetNDF());
               legend->AddEntry(fitFcnV2_ZYAM,message);
        if(res->CovMatrixStatus() == 3){
                 //  sprintf(message,"The fit is a success");
               }
               else{
                   sprintf(message,"The fit is a failure");
                   legend->AddEntry(fitFcnV2_ZYAM,message);
               }
             legend->AddEntry(Yields_Central_1_ZYAM,"Data","lpe");
             legend->Draw();

        sprintf(CanvasName,"%s/V2_ZYAM.pdf",CanvasFolderName);
        cZYAM->SaveAs(CanvasName);
    }
    
    // Fit PRL Template
    
    TCanvas* cPRLTemplate = new TCanvas;
       cPRLTemplate->SetTitle("PRL Template fit");
       cPRLTemplate->Divide(1,1);
       cPRLTemplate->cd(1);
       Yields_Central_1->DrawCopy();
//       cPRLTemplate->cd(2);
//       Baseline_Central_1->DrawCopy();
//       cPRLTemplate->cd(3);
//       Yields_Central_1_MinusBaseline->DrawCopy();
       
       
       {
           cPRLTemplate->cd(1);
                TF1 *fitFcnV2_PRLTemplate = new TF1("fitFcnV2_PRLTemplate",PRLTemplate,0,TMath::Pi(),2);//-TMath::Pi()/2,1.5*TMath::Pi(),2);
                fitFcnV2_PRLTemplate->SetNpx(500);
                fitFcnV2_PRLTemplate->SetLineWidth(4);
                fitFcnV2_PRLTemplate->SetLineColor(kRed+1);
           TF1 *fitFcnV2_PRLTemplate_RidgeAndZero = new TF1("fitFcnV2_PRLTemplate_RidgeAndZero",PRLTemplate_RidgeAndZero,0,TMath::Pi(),2);//-TMath::Pi()/2,1.5*TMath::Pi(),2);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetNpx(500);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetLineWidth(2);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetLineStyle(9);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetLineColor(kRed);
           TF1 *fitFcnV2_PRLTemplate_PeriphAndG = new TF1("fitFcnV2_PRLTemplate_PeriphAndG",PRLTemplate_PeriphAndG,0,TMath::Pi(),2);//-TMath::Pi()/2,1.5*TMath::Pi(),2);
           fitFcnV2_PRLTemplate_PeriphAndG->SetNpx(500);
           fitFcnV2_PRLTemplate_PeriphAndG->SetLineWidth(2);
           fitFcnV2_PRLTemplate_PeriphAndG->SetLineStyle(9);
           fitFcnV2_PRLTemplate_PeriphAndG->SetLineColor(kBlack);
                // first try without starting values for the parameters
                // This defaults to 1 for each param.
                // this results in an ok fit for the polynomial function
                // however the non-linear part (lorenzian) does not
                // respond well.
                 Double_t params[2] = {1,1};
                fitFcnV2_PRLTemplate->SetParameters(params);
                 TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->SetMaxIterations(10000);
                 TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->SetPrecision();
                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->SetFCN(ChisquarePRLTemplate);
              // gStyle->SetOptFit(1011);
               //  histo->Fit("fitFcn","0");
                // second try: set start values for some parameters

                 fitFcnV2_PRLTemplate->SetParName(0,"V2");
                   fitFcnV2_PRLTemplate->SetParName(1,"F");
           
         //  fitFcnV2_PRLTemplate->SetParLimits(0,0.0001,5000000);
        //   fitFcnV2_PRLTemplate->FixParameter(0,0.002);
           fitFcnV2_PRLTemplate->SetParLimits(1,0.1,50);
                 //  fitFcnV2_PRLTemplate->FixParameter(1,2.5);
                 

                TFitResultPtr res = Yields_Central_1_PRLTemplate->Fit("fitFcnV2_PRLTemplate","USBMERI+","ep");
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)Yields_Central_1_PRLTemplate->FindObject("stats");
                st->SetX1NDC(0.8); //new x start position
                st->SetY1NDC(0.8); //new x end position
           
               double chi2, edm, errdef;
               int nvpar, nparx;
                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate)->GetStats(chi2,edm,errdef,nvpar,nparx);
                fitFcnV2_PRLTemplate->SetChisquare(chi2);
                int ndf = npfits-nvpar;
                fitFcnV2_PRLTemplate->SetNDF(ndf);
           
               Double_t par[2];
               fitFcnV2_PRLTemplate->GetParameters(par);
           fitFcnV2_PRLTemplate_RidgeAndZero->SetParameters(par);
           fitFcnV2_PRLTemplate_PeriphAndG->SetParameters(par);
                // improve the pictu
              //   std::cout << "integral error: " << integralerror << std::endl;
                fitFcnV2_PRLTemplate->Draw("same");
           fitFcnV2_PRLTemplate_RidgeAndZero->Draw("same");
           fitFcnV2_PRLTemplate_PeriphAndG->Draw("same");
                // draw the legend
                TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                legend->SetFillColorAlpha(kWhite, 0.);
                legend->SetBorderSize(0);
                  legend->SetTextFont(42);
                legend->SetTextSize(0.03);
                  Char_t message[80];
                  sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate->GetChisquare(),fitFcnV2_PRLTemplate->GetNDF());
                  legend->AddEntry(fitFcnV2_PRLTemplate,message);
           if(res->CovMatrixStatus() == 3){
                   //   sprintf(message,"The fit is a success");
                  }
                  else{
                      sprintf(message,"The fit is a failure");
                      legend->AddEntry(fitFcnV2_PRLTemplate,message);
                  }
                legend->AddEntry(Yields_Central_1_PRLTemplate,"Data","lpe");
                legend->AddEntry(fitFcnV2_PRLTemplate_RidgeAndZero,"F*Baseline_{Periph} + Ridge*G", "L");
                legend->AddEntry(fitFcnV2_PRLTemplate_PeriphAndG,"F*Yield_{Periph} + G", "L");
                legend->Draw();

           sprintf(CanvasName,"%s/V2_ATLAS.pdf",CanvasFolderName);
           cPRLTemplate->SaveAs(CanvasName);
       }
    
    // Fit PRL Template + ZYAM Periph
       
       TCanvas* cPRLTemplate_PeriphZYAM = new TCanvas;
          cPRLTemplate_PeriphZYAM->SetTitle("PRL Template fit Periph ZYAM");
          cPRLTemplate_PeriphZYAM->Divide(1,1);
          cPRLTemplate_PeriphZYAM->cd(1);
          Yields_Central_1->DrawCopy();
//          cPRLTemplate_PeriphZYAM->cd(2);
//          Baseline_Central_1->DrawCopy();
//          cPRLTemplate_PeriphZYAM->cd(3);
//          Yields_Central_1_MinusBaseline->DrawCopy();
          
          
          {
              cPRLTemplate_PeriphZYAM->cd(1);
                   TF1 *fitFcnV2_PRLTemplate_PeriphZYAM = new TF1("fitFcnV2_PRLTemplate_PeriphZYAM",PRLTemplate_PeriphZYAM,0,TMath::Pi(),2);//-TMath::Pi()/2,1.5*TMath::Pi(),2);
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetNpx(500);
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetLineWidth(4);
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetLineColor(kBlack);
                   // first try without starting values for the parameters
                   // This defaults to 1 for each param.
                   // this results in an ok fit for the polynomial function
                   // however the non-linear part (lorenzian) does not
                   // respond well.
                    Double_t params[2] = {1,1};
                   fitFcnV2_PRLTemplate_PeriphZYAM->SetParameters(params);
                    TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->SetMaxIterations(10000);
                    TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->SetPrecision();
                    TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->SetFCN(ChisquarePRLTemplate_PeriphZYAM);
                // gStyle->SetOptFit(1011);
                  //  histo->Fit("fitFcn","0");
                   // second try: set start values for some parameters

                    fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(0,"V2");
                      fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(1,"F");
              
            //  fitFcnV2_PRLTemplate->SetParLimits(0,0.0001,5000000);
           //   fitFcnV2_PRLTemplate->FixParameter(0,0.002);
              fitFcnV2_PRLTemplate_PeriphZYAM->SetParLimits(1,0.1,5000);
                    //  fitFcnV2_PRLTemplate->FixParameter(1,2.5);
                    

                   TFitResultPtr res = Yields_Central_1_PRLTemplate_PeriphZYAM->Fit("fitFcnV2_PRLTemplate_PeriphZYAM","USBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)Yields_Central_1_PRLTemplate_PeriphZYAM->FindObject("stats");
                    st->SetX1NDC(0.8); //new x start position
                    st->SetY1NDC(0.8); //new x end position
              
              double chi2, edm, errdef;
              int nvpar, nparx;
               TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAM)->GetStats(chi2,edm,errdef,nvpar,nparx);
               fitFcnV2_PRLTemplate_PeriphZYAM->SetChisquare(chi2);
               int ndf = npfits-nvpar;
               fitFcnV2_PRLTemplate_PeriphZYAM->SetNDF(ndf);
              
                  Double_t par[2];
                  fitFcnV2_PRLTemplate_PeriphZYAM->GetParameters(par);
                   // improve the pictu
                 //   std::cout << "integral error: " << integralerror << std::endl;
                   fitFcnV2_PRLTemplate_PeriphZYAM->Draw("same");
                   // draw the legend
                   TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                   legend->SetFillColorAlpha(kWhite, 0.);
                   legend->SetBorderSize(0);
                     legend->SetTextFont(42);
                   legend->SetTextSize(0.03);
                     Char_t message[80];
                     sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate_PeriphZYAM->GetChisquare(),fitFcnV2_PRLTemplate_PeriphZYAM->GetNDF());
                     legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
              if(res->CovMatrixStatus() == 3){
                       //  sprintf(message,"The fit is a success");
                     }
                     else{
                         sprintf(message,"The fit is a failure");
                         legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                     }
                   legend->AddEntry(Yields_Central_1_PRLTemplate_PeriphZYAM,"Data","lpe");
                   legend->Draw();
              
              sprintf(CanvasName,"%s/V2_ATLASZYAM.pdf",CanvasFolderName);
              cPRLTemplate_PeriphZYAM->SaveAs(CanvasName);
          }
    
    
    
    
    
    
    
    
    
    
       //DimuTkl Ext1 Pt Binned
    
    
    if(PtBinned){
        
        for(int ptbin=0;ptbin<NbPtBins;ptbin++){
            
            baseline_centralPtBinned[ptbin] = 9999.;
    
        int niterations = 1;
        cinvmassPtBinned->cd();
        sprintf(histoname,Form("bin%d_1",ptbin+1));
        res = FittingAllInvMassBin(histoname, cinvmassPtBinned, NbPtBins + ptbin, SignalF, BackgroundF, minmass, maxmass, ratsigma);
            cout << "Out of this"<<endl;
       double parerr[numParameters];
        TCanvas* c10_fitPtBinned = new TCanvas;
           c10_fitPtBinned->Divide(NbinsDeltaPhi/2,4);
        TRandom rand = 0;
        for(int j=0; j<numParameters; j++){
            cout << "Out of this1"<<endl;
            parerr[j] = res->ParError(j);
        }
            cout << "Out of this2"<<endl;
        for(int j=0; j <numParameters+numParametersBkgV2+1; j++){
            par[j] = res->Parameter(j);
            if(j>=numParameters){
                par[j] = 3;
            }
        }
            cout << "Out of this3"<<endl;
            //PLOT ICI

        double par12[niterations];
        double parerr12[niterations];

        
        for(int i=0; i<NbinsDeltaPhi; i++){ //NbPhiBins
            for(int t=0; t<niterations; t++){
                TVirtualPad* subpad = c10_fitPtBinned->cd(i+1);
                c10_fitPtBinned->cd(i+1);
                if(t==niterations-1){
                    subpad->Clear();
                }
            TF1 *fitY_1Central = new TF1("fitY_1Central",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
              fitY_1Central->SetNpx(500);
              fitY_1Central->SetLineWidth(4);
              fitY_1Central->SetLineColor(kMagenta);
              fitY_1Central->SetParameters(par);
               TVirtualFitter::Fitter(YieldWrtMass_CentralPtBinned[i][ptbin])->SetMaxIterations(10000);
               TVirtualFitter::Fitter(YieldWrtMass_CentralPtBinned[i][ptbin])->SetPrecision();
            //  histo->Fit("fitFcn","0");
              // second try: set start values for some parameters
            for(int k=0; k<numParameters; k++){

                    fitY_1Central->FixParameter(k,par[k]);
                
            }
                fitY_1Central->SetParLimits(numParameters,0,10000);
                
                for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                    string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "Y_1");
                    const char *titular = title.c_str();
                    fitY_1Central->SetParName(index,titular);
                }
                
//               fitY_1Central->SetParName(0,"Norm_{J/#psi}");
//               fitY_1Central->SetParName(1,"M_{J/#psi}");
//               fitY_1Central->SetParName(2,"Sigma_{J/#psi}");
//               fitY_1Central->SetParName(3,"a_{1}");
//               fitY_1Central->SetParName(4,"n_{1}");
//               fitY_1Central->SetParName(5,"a_{2}");
//               fitY_1Central->SetParName(6,"n_{2}");
//               fitY_1Central->SetParName(7,"Norm_{#Psi(2S)}");
//               fitY_1Central->SetParName(8,"Norm_{TailLowM}");
//               fitY_1Central->SetParName(9,"Exp_{TailLowM}");
//               fitY_1Central->SetParName(10,"Norm_{TailHighM}");
//               fitY_1Central->SetParName(11,"Exp_{TailHighM}");
//                fitY_1Central->SetParName(12,"Y_1 J/#psi");
//                fitY_1Central->SetParName(13,"Y_1 Bkg M2");
//            fitY_1Central->SetParName(14,"Y_1 Bkg M1");
//            fitY_1Central->SetParName(15,"Y_1 Bkg M0");

              rescent = YieldWrtMass_CentralPtBinned[i][ptbin]->Fit("fitY_1Central","SBMERIQ+","ep");
                gStyle->SetOptStat("n");
                gStyle->SetOptFit(1011);
                TPaveStats *st = (TPaveStats*)YieldWrtMass_CentralPtBinned[i][ptbin]->FindObject("stats");
                st->SetX1NDC(0.8); //new x start position
                st->SetY1NDC(0.8); //new x end position
              Double_t param[numParameters+numParametersBkgV2+1];
                Double_t paramerrs[numParameters+numParametersBkgV2+1];
              fitY_1Central->GetParameters(param);
                for(int i=0; i<numParameters+numParametersBkgV2+1; i++){
                    paramerrs[i] = fitY_1Central->GetParError(i);
                }
                
                if(rescent->CovMatrixStatus() == 3){
                par12[t] = param[numParameters];
                parerr12[t] = paramerrs[numParameters];
                }
                else{
                    par12[t]=parerr12[t]=0;
                }
                
                std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
                std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

            TF1 *backFcnY_1Central = new TF1("backFcnY_1Central",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
            backFcnY_1Central->SetLineColor(kRed);
            TF1 *signalFcnJPsiY_1Central = new TF1("signalFcnJPsiY_1Central",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
              // writes the fit results into the par array
            
                if(t==niterations-1){
                    YieldWrtMass_CentralPtBinned[i][ptbin]->GetListOfFunctions()->Clear();
                    YieldWrtMass_CentralPtBinned[i][ptbin]->GetListOfFunctions()->Add(fitY_1Central);
                }
            signalFcnJPsiY_1Central->SetLineColor(kBlue);
            signalFcnJPsiY_1Central->SetNpx(500);
            backFcnY_1Central->SetParameters(param);
            signalFcnJPsiY_1Central->SetParameters(param);
          //  signalFcnJPsiY_1Central->Draw("same");
            backFcnY_1Central->Draw("same");
// Calcul baseline central
                
                if(param[numParameters]<baseline_centralPtBinned[ptbin]){ //En central baseline = minimum
                    baseline_centralPtBinned[ptbin] = param[numParameters];
                    errbaseline_centralPtBinned[ptbin] = fitY_1Central->GetParError(numParameters);
                }
                
//                if(i==2){ // Si baseline = ... DeltaPhi = 0
//                    baseline_centralPtBinned[ptbin] = param[12];
//                    errbaseline_centralPtBinned[ptbin] = fitY_1Central->GetParError(12);
//                }
//                if(i==3){
//                    baseline_centralPtBinned[ptbin] += param[12];
//                    baseline_centralPtBinned[ptbin] /= 2;
//                    errbaseline_centralPtBinned[ptbin] = sqrt(pow(errbaseline_centralPtBinned[ptbin],2)+pow(fitY_1Central->GetParError(12),2));
//                }
                

                    Yields_Central_1PtBinned[ptbin]->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[numParameters]); //AEUGG
                    Yields_Central_1PtBinned[ptbin]->SetBinError(i+1,fitY_1Central->GetParError(numParameters));
                
                
              // writes the fit results into the par array
              // gStyle->SetOptFit(1011);
              // draw the legend
              TLegend *legend=new TLegend(0.12,0.75,0.60,0.90);
              legend->SetFillColorAlpha(kWhite, 0.);
              legend->SetBorderSize(0);
                legend->SetTextFont(42);
              legend->SetTextSize(0.03);
            Char_t message[80];
            sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitY_1Central->GetChisquare(),fitY_1Central->GetNDF());
             legend->AddEntry(fitY_1Central,message);
            if(rescent->CovMatrixStatus() == 3){
                     //  sprintf(message,"The fit is a success");
                   }
                   else{
                       sprintf(message,"The fit is a failure");
                       legend->AddEntry(fitY_1Central,message);
                   }
         //   legend->AddEntry(signalFcnJPsiY_1Central,"JPsi signal");
            legend->AddEntry(backFcnY_1Central,"Background");
              legend->AddEntry(YieldWrtMass_CentralPtBinned[i][ptbin],"Data","lpe");
              legend->Draw();
            }
            
        }
            
            for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
                Baseline_Central_1PtBinned[ptbin]->SetBinContent(phi_idx+1, baseline_centralPtBinned[ptbin]);
                Baseline_Central_1PtBinned[ptbin]->SetBinError(phi_idx+1, errbaseline_centralPtBinned[ptbin]);
            }
            
            Yields_Central_1_MinusBaselinePtBinned[ptbin]->Add(Yields_Central_1PtBinned[ptbin],Baseline_Central_1PtBinned[ptbin],1,-1);
            
        cinvmassPtBinned->cd(); //aeugg
        sprintf(histoname,Form("bin%d_2",ptbin+1));
        res = FittingAllInvMassBin(histoname, cinvmassPtBinned, NbPtBins + NbPtBins + ptbin, SignalF, BackgroundF, minmass, maxmass, ratsigma);
        
           for(int j=0; j<numParameters; j++){
               parerr[j] = res->ParError(j);
           }
           for(int j=0; j <numParameters+numParametersBkgV2+1; j++){
               par[j] = res->Parameter(j);
               if(j>=numParameters){
                   par[j] = 3;
               }
           }
        
           
           for(int i=0; i<NbinsDeltaPhi; i++){ //NbPhiBins
               for(int t=0; t<niterations; t++){
                   TVirtualPad* subpad = c10_fitPtBinned->cd(NbinsDeltaPhi+i+1);
                   c10_fitPtBinned->cd(NbinsDeltaPhi+i+1);
                   if(t==niterations-1){
                       subpad->Clear();
                   }
               TF1 *fitY_1Periph = new TF1("fitY_1Periph",FourierV2_WrtInvMass,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 fitY_1Periph->SetNpx(500);
                 fitY_1Periph->SetLineWidth(4);
                 fitY_1Periph->SetLineColor(kMagenta);
                 fitY_1Periph->SetParameters(par);
                  TVirtualFitter::Fitter(YieldWrtMass_PeriphPtBinned[i][ptbin])->SetMaxIterations(10000);
                  TVirtualFitter::Fitter(YieldWrtMass_PeriphPtBinned[i][ptbin])->SetPrecision();
               //  histo->Fit("fitFcn","0");
                 // second try: set start values for some parameters
                for(int k=0; k<numParameters; k++){
                        fitY_1Periph->FixParameter(k,par[k]);
                }
                  fitY_1Periph->SetParLimits(numParameters,0,10000);

                   for(int index=0; index<numParameters+numParametersBkgV2+1; index++){
                       string title = GetParameterInfo(index, GlobalSignal, GlobalBackground, GlobalBackgroundV2, "Y_1");
                       const char *titular = title.c_str();
                       fitY_1Periph->SetParName(index,titular);
                   }
                   
//                  fitY_1Periph->SetParName(0,"Norm_{J/#psi}");
//                  fitY_1Periph->SetParName(1,"M_{J/#psi}");
//                  fitY_1Periph->SetParName(2,"Sigma_{J/#psi}");
//                  fitY_1Periph->SetParName(3,"a_{1}");
//                  fitY_1Periph->SetParName(4,"n_{1}");
//                  fitY_1Periph->SetParName(5,"a_{2}");
//                  fitY_1Periph->SetParName(6,"n_{2}");
//                  fitY_1Periph->SetParName(7,"Norm_{#Psi(2S)}");
//                  fitY_1Periph->SetParName(8,"Norm_{TailLowM}");
//                  fitY_1Periph->SetParName(9,"Exp_{TailLowM}");
//                  fitY_1Periph->SetParName(10,"Norm_{TailHighM}");
//                  fitY_1Periph->SetParName(11,"Exp_{TailHighM}");
//                   fitY_1Periph->SetParName(12,"Y_1 J/#psi");
//                   fitY_1Periph->SetParName(13,"Y_1 Bkg M2");
//               fitY_1Periph->SetParName(14,"Y_1 Bkg M1");
//               fitY_1Periph->SetParName(15,"Y_1 Bkg M0");
                  
                 resperiph = YieldWrtMass_PeriphPtBinned[i][ptbin]->Fit("fitY_1Periph","SBMERIQ+","ep");
                   gStyle->SetOptStat("n");
                   gStyle->SetOptFit(1011);
                   TPaveStats *st = (TPaveStats*)YieldWrtMass_PeriphPtBinned[i][ptbin]->FindObject("stats");
                   st->SetX1NDC(0.8); //new x start position
                   st->SetY1NDC(0.8); //new x end position
                   Double_t param[numParameters+numParametersBkgV2+1];
                   Double_t paramerrs[numParameters+numParametersBkgV2+1];
                 fitY_1Periph->GetParameters(param);
                   
                   for(int i=0; i<numParameters+numParametersBkgV2+1; i++){
                                   paramerrs[i] = fitY_1Periph->GetParError(i);
                               }
                               
                               if(resperiph->CovMatrixStatus() == 3){
                               par12[t] = param[numParameters];
                               parerr12[t] = paramerrs[numParameters];
                               }
                               else{
                                   par12[t]=parerr12[t]=0;
                               }
                   std::cout << "COV status Central t=" << t << " : " << rescent->CovMatrixStatus()<<endl;
                    std::cout << "par12: " << par12[t] << ", parerr12: " << parerr12[t] <<endl;

               TF1 *backFcnY_1Periph = new TF1("backFcnY_1Periph",BackFcnV2Poly,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
               backFcnY_1Periph->SetLineColor(kRed);
               TF1 *signalFcnJPsiY_1Periph = new TF1("signalFcnJPsiY_1Periph",SignalFcnJPsiV2,MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
                 // writes the fit results into the par array
                   
                   if(t==niterations-1){
                       YieldWrtMass_PeriphPtBinned[i][ptbin]->GetListOfFunctions()->Clear();
                       YieldWrtMass_PeriphPtBinned[i][ptbin]->GetListOfFunctions()->Add(fitY_1Periph);
                   }
               
               signalFcnJPsiY_1Periph->SetLineColor(kBlue);
               signalFcnJPsiY_1Periph->SetNpx(500);
               backFcnY_1Periph->SetParameters(param);
               signalFcnJPsiY_1Periph->SetParameters(param);
            //   signalFcnJPsiY_1Periph->Draw("same");
               backFcnY_1Periph->Draw("same");
                   
                      if(i==0){
                           baseline_periphPtBinned[ptbin] = param[numParameters];
                           errbaseline_periphPtBinned[ptbin] = fitY_1Periph->GetParError(numParameters);
                       }
//                       if(i==3){
//                           baseline_periphPtBinned[ptbin] += param[numParameters];
//                           baseline_periphPtBinned[ptbin] /= 2;
//                           errbaseline_periphPtBinned[ptbin] = sqrt(pow(errbaseline_periphPtBinned[ptbin],2)+pow(fitY_1Periph->GetParError(numParameters),2));
//                           errbaseline_periphPtBinned[ptbin] /= 2;
//                       }
                   
                     Yields_Periph_1PtBinned[ptbin]->Fill(MinDeltaPhi + (i+0.5)*SizeBinDeltaPhi, param[numParameters]);
                     Yields_Periph_1PtBinned[ptbin]->SetBinError(i+1,fitY_1Periph->GetParError(numParameters));
                   
                   
                 // writes the fit results into the par array
                 // gStyle->SetOptFit(1011);
                 // draw the legend
                 TLegend *legend=new TLegend(0.12,0.75,0.60,0.90);
                 legend->SetFillColorAlpha(kWhite, 0.);
                 legend->SetBorderSize(0);
                   legend->SetTextFont(42);
                 legend->SetTextSize(0.03);
               Char_t message[80];
               sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitY_1Periph->GetChisquare(),fitY_1Periph->GetNDF());
                legend->AddEntry(fitY_1Periph,message);
               if(resperiph->CovMatrixStatus() == 3){
                         // sprintf(message,"The fit is a success");
                      }
                      else{
                          sprintf(message,"The fit is a failure");
                          legend->AddEntry(fitY_1Periph,message);
                      }
                      
          //     legend->AddEntry(signalFcnJPsiY_1Periph,"JPsi signal");
               legend->AddEntry(backFcnY_1Periph,"Background");
                 legend->AddEntry(YieldWrtMass_Periph[i],"Data","lpe");
                 legend->Draw();
               }
               
           }
            
            for(int phi_idx = 0; phi_idx<NbinsDeltaPhi; phi_idx++){
                Baseline_Periph_1PtBinned[ptbin]->SetBinContent(phi_idx+1, baseline_periphPtBinned[ptbin]);
                Baseline_Periph_1PtBinned[ptbin]->SetBinError(phi_idx+1, errbaseline_periphPtBinned[ptbin]);
            }
            
            Yields_Periph_1_MinusBaselinePtBinned[ptbin]->Add(Yields_Periph_1PtBinned[ptbin],Baseline_Periph_1PtBinned[ptbin],1,-1);
            
            sprintf(CanvasName,"%s/Idk_PtBin%i.pdf",CanvasFolderName,ptbin);
            c10_fitPtBinned->SaveAs(CanvasName);
            
            //Fill Baseline Periph and Y-B
            char hname[200];
            TCanvas* cYBcentralPtBinned = new TCanvas;
            sprintf(hname, "Central Yield-Baseline - PtBinned %d", ptbin);
            cYBcentralPtBinned->SetTitle(hname);
            cYBcentralPtBinned->Divide(1,3);
            cYBcentralPtBinned->cd(1);
            Yields_Central_1PtBinned[ptbin]->DrawCopy();
            cYBcentralPtBinned->cd(2);
            Baseline_Central_1PtBinned[ptbin]->DrawCopy();
            cYBcentralPtBinned->cd(3);
            Yields_Central_1_MinusBaselinePtBinned[ptbin]->DrawCopy();
            
            sprintf(CanvasName,"%s/CentralYield-BaselinePtBin%i.pdf",CanvasFolderName,ptbin);
            cYBcentralPtBinned->SaveAs(CanvasName);
            
            TCanvas* cYBperiphPtBinned = new TCanvas;
            sprintf(hname, "Periph Yield-Baseline - PtBinned %d", ptbin);
            cYBperiphPtBinned->SetTitle(hname);
            cYBperiphPtBinned->Divide(1,3);
            cYBperiphPtBinned->cd(1);
            Yields_Periph_1PtBinned[ptbin]->DrawCopy();
            cYBperiphPtBinned->cd(2);
            Baseline_Periph_1PtBinned[ptbin]->DrawCopy();
            cYBperiphPtBinned->cd(3);
            Yields_Periph_1_MinusBaselinePtBinned[ptbin]->DrawCopy();
            
            sprintf(CanvasName,"%s/PeriphYield-BaselinePtBin%i.pdf",CanvasFolderName,ptbin);
            cYBperiphPtBinned->SaveAs(CanvasName);
            
            sprintf(hname, "Yields_Central_1_CvetanPtBinned - PtBinned %d", ptbin);
            Yields_Central_1_CvetanPtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_CvetanMePtBinned - PtBinned %d", ptbin);
            Yields_Central_1_CvetanMePtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_ZYAMPtBinned - PtBinned %d", ptbin);
            Yields_Central_1_ZYAMPtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_PRLTemplatePtBinned - PtBinned %d", ptbin);
            Yields_Central_1_PRLTemplatePtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            sprintf(hname, "Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned - PtBinned %d", ptbin);
            Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin] = (TH1F*)Yields_Central_1PtBinned[ptbin]->Clone(hname);
            
            {
                TCanvas* c17PtBinned = new TCanvas;
                c17PtBinned->Divide(1,3);
                //Creer canvas pour imprimer les plots Periph yield et Central yield wrt phi et leur difference
                c17PtBinned->cd(1);
                Yields_Central_1PtBinned[ptbin]->Draw();
                c17PtBinned->cd(2);
                Yields_Periph_1PtBinned[ptbin]->Draw();
                c17PtBinned->cd(3);
                Yields_Difference_1PtBinned[ptbin]->Add(Yields_Central_1PtBinned[ptbin],Yields_Periph_1PtBinned[ptbin],1,-1);
                Yields_Difference_1PtBinned[ptbin]->Draw();
            
                     TF1 *fitFcnV2_2 = new TF1("fitFcnV2_2",FourierV2,0,TMath::Pi(),3);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
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
                      TVirtualFitter::Fitter(Yields_Difference_1PtBinned[ptbin])->SetMaxIterations(10000);
                      TVirtualFitter::Fitter(Yields_Difference_1PtBinned[ptbin])->SetPrecision();
                   //  histo->Fit("fitFcn","0");
                     // second try: set start values for some parameters
                      
                      fitFcnV2_2->SetParName(0,"a0");
                      fitFcnV2_2->SetParName(1,"a1");
                      fitFcnV2_2->SetParName(2,"a2");
                      
                     TFitResultPtr res = Yields_Difference_1PtBinned[ptbin]->Fit("fitFcnV2_2","SBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)Yields_Difference_1PtBinned[ptbin]->FindObject("stats");
                    st->SetX1NDC(0.8); //new x start position
                    st->SetY1NDC(0.8); //new x end position
                    Double_t par[3];
                    fitFcnV2_2->GetParameters(par);
                     // improve the pictu
                   //   std::cout << "integral error: " << integralerror << std::endl;
                     fitFcnV2_2->Draw("same");
                     // draw the legend
                     TLegend *legend=new TLegend(0.15,0.70,0.4,0.90);
                     legend->SetFillColorAlpha(kWhite, 0.);
                     legend->SetBorderSize(0);
                       legend->SetTextFont(42);
                     legend->SetTextSize(0.03);
                       Char_t message[80];
                       sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
                       legend->AddEntry(fitFcnV2_2,message);
                sprintf(message,"V2_1 J/#psi-tkl: #frac{%.4f}{%.4f + %.4f} = %.4f +- %.4f",par[2],par[0], baseline_periphPtBinned[ptbin],par[2]/(par[0] + baseline_periphPtBinned[ptbin]),abs(par[2]/(par[0] + baseline_periphPtBinned[ptbin])*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(sqrt(pow(fitFcnV2_2->GetParError(0),2)+pow(errbaseline_periphPtBinned[ptbin],2))/(par[0] + baseline_periphPtBinned[ptbin]),2))));
                legend->AddEntry(fitFcnV2_2,message);
                if(res->CovMatrixStatus() == 3){
                          // sprintf(message,"The fit is a success");
                       }
                       else{
                           sprintf(message,"The fit is a failure");
                           legend->AddEntry(fitFcnV2_2,message);
                       }
                       
                     legend->AddEntry(Yields_Difference_1,"Data","lpe");
                     legend->Draw();
                
                V2_Ext1[ptbin] = par[2]/(par[0] + baseline_periphPtBinned[ptbin]);
                errV2_Ext1[ptbin] = abs(par[2]/(par[0] + baseline_periphPtBinned[ptbin])*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(sqrt(pow(fitFcnV2_2->GetParError(0),2)+pow(errbaseline_periphPtBinned[ptbin],2))/(par[0] + baseline_periphPtBinned[ptbin]),2)));
                
                V2_Ext1_noZYAM[ptbin] = par[2]/(par[0]);
                errV2_Ext1_noZYAM[ptbin] = abs(par[2]/(par[0])*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(fitFcnV2_2->GetParError(0)/par[0],2)));
                   
                sprintf(CanvasName,"%s/V2_Ext1_PtBin%i.pdf",CanvasFolderName,ptbin);
                c17PtBinned->SaveAs(CanvasName);
            }
            // Fit Ext1 V3
            {
                TCanvas* c17PtBinned3 = new TCanvas;
                c17PtBinned3->Divide(1,3);
                //Creer canvas pour imprimer les plots Periph yield et Central yield wrt phi et leur difference
                c17PtBinned3->cd(1);
                Yields_Central_1PtBinned[ptbin]->Draw();
                c17PtBinned3->cd(2);
                Yields_Periph_1PtBinned[ptbin]->Draw();
                c17PtBinned3->cd(3);
                Yields_Difference_1PtBinned[ptbin]->Add(Yields_Central_1PtBinned[ptbin],Yields_Periph_1PtBinned[ptbin],1,-1);
                Yields_Difference_1PtBinned[ptbin]->Draw();
            
                     TF1 *fitFcnV2_2 = new TF1("fitFcnV2_2",FourierV3,0,TMath::Pi(),4);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
                     fitFcnV2_2->SetNpx(500);
                     fitFcnV2_2->SetLineWidth(4);
                     fitFcnV2_2->SetLineColor(kMagenta);
                     // first try without starting values for the parameters
                     // This defaults to 1 for each param.
                     // this results in an ok fit for the polynomial function
                     // however the non-linear part (lorenzian) does not
                     // respond well.
                      Double_t params[4] = {1,1,1,1};
                     fitFcnV2_2->SetParameters(params);
                      TVirtualFitter::Fitter(Yields_Difference_1PtBinned[ptbin])->SetMaxIterations(10000);
                      TVirtualFitter::Fitter(Yields_Difference_1PtBinned[ptbin])->SetPrecision();
                   //  histo->Fit("fitFcn","0");
                     // second try: set start values for some parameters
                      
                      fitFcnV2_2->SetParName(0,"a0");
                      fitFcnV2_2->SetParName(1,"a1");
                      fitFcnV2_2->SetParName(2,"a2");
                    fitFcnV2_2->SetParName(3,"a3");
                      
                     TFitResultPtr res = Yields_Difference_1PtBinned[ptbin]->Fit("fitFcnV2_2","SBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)Yields_Difference_1PtBinned[ptbin]->FindObject("stats");
                    st->SetX1NDC(0.8); //new x start position
                    st->SetY1NDC(0.8); //new x end position
                    Double_t par[3];
                    fitFcnV2_2->GetParameters(par);
                     // improve the pictu
                   //   std::cout << "integral error: " << integralerror << std::endl;
                     fitFcnV2_2->Draw("same");
                     // draw the legend
                     TLegend *legend=new TLegend(0.15,0.70,0.4,0.90);
                     legend->SetFillColorAlpha(kWhite, 0.);
                     legend->SetBorderSize(0);
                       legend->SetTextFont(42);
                     legend->SetTextSize(0.03);
                       Char_t message[80];
                       sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_2->GetChisquare(),fitFcnV2_2->GetNDF());
                       legend->AddEntry(fitFcnV2_2,message);
                sprintf(message,"V2_1 J/#psi-tkl: #frac{%.4f}{%.4f + %.4f} = %.4f +- %.4f",par[2],par[0], baseline_periphPtBinned[ptbin],par[2]/(par[0] + baseline_periphPtBinned[ptbin]),abs(par[2]/(par[0] + baseline_periphPtBinned[ptbin])*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(sqrt(pow(fitFcnV2_2->GetParError(0),2)+pow(errbaseline_periphPtBinned[ptbin],2))/(par[0] + baseline_periphPtBinned[ptbin]),2))));
                legend->AddEntry(fitFcnV2_2,message);
                if(res->CovMatrixStatus() == 3){
                          // sprintf(message,"The fit is a success");
                       }
                       else{
                           sprintf(message,"The fit is a failure");
                           legend->AddEntry(fitFcnV2_2,message);
                       }
                       
                     legend->AddEntry(Yields_Difference_1,"Data","lpe");
                     legend->Draw();
                
                V2_Ext1_V3[ptbin] = par[2]/(par[0] + baseline_periphPtBinned[ptbin]);
                errV2_Ext1_V3[ptbin] = abs(par[2]/(par[0] + baseline_periphPtBinned[ptbin])*sqrt(pow(fitFcnV2_2->GetParError(2)/par[2],2)+pow(sqrt(pow(fitFcnV2_2->GetParError(0),2)+pow(errbaseline_periphPtBinned[ptbin],2))/(par[0] + baseline_periphPtBinned[ptbin]),2)));
                   
                sprintf(CanvasName,"%s/V2_Ext1_V3_PtBin%i.pdf",CanvasFolderName,ptbin);
                c17PtBinned3->SaveAs(CanvasName);
            }
            
            
             // Fit Cvetan PtBinned
                
                TCanvas* cCvetanPtBinned = new TCanvas;
            sprintf(hname, "Cvetan fit - PtBinned %d", ptbin);
                   cCvetanPtBinned->SetTitle(hname);
                   cCvetanPtBinned->Divide(1,1);
                   cCvetanPtBinned->cd(1);
                   Yields_Central_1->DrawCopy();
            //       cCvetan->cd(2);
            //       Baseline_Central_1->DrawCopy();
            //       cCvetan->cd(3);
            //       Yields_Central_1_MinusBaseline->DrawCopy();
            
            
                
               {
                    cCvetanPtBinned->cd(1);
                         TF1 *fitFcnV2_Cvetan = new TF1("fitFcnV2_Cvetan",CvetanFPtBinned,0,TMath::Pi(),5);//-TMath::Pi()/2,1.5*TMath::Pi(),5);
                         fitFcnV2_Cvetan->SetNpx(500);
                         fitFcnV2_Cvetan->SetLineWidth(4);
                         fitFcnV2_Cvetan->SetLineColor(kOrange+3);
                         // first try without starting values for the parameters
                         // This defaults to 1 for each param.
                         // this results in an ok fit for the polynomial function
                         // however the non-linear part (lorenzian) does not
                         // respond well.
                          Double_t params[5] = {static_cast<Double_t>(ptbin),1,0,0.01,1};
                         fitFcnV2_Cvetan->SetParameters(params);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanPtBinned[ptbin])->SetMaxIterations(10000);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanPtBinned[ptbin])->SetPrecision();
                        TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetFCN(ChisquareCvetanFPtBinned);
                       
                       //  histo->Fit("fitFcn","0");
                         // second try: set start values for some parameters

                        fitFcnV2_Cvetan->SetParName(0,"ptbin");
                          fitFcnV2_Cvetan->SetParName(1,"V0");
                            fitFcnV2_Cvetan->SetParName(2,"V1");
                            fitFcnV2_Cvetan->SetParName(3,"V2");
                          fitFcnV2_Cvetan->SetParName(4,"F");
                    
                    fitFcnV2_Cvetan->SetParLimits(4,0.1,50);
                    
                    fitFcnV2_Cvetan->FixParameter(0,ptbin);
                    fitFcnV2_Cvetan->FixParameter(1,1);
                    fitFcnV2_Cvetan->FixParameter(2,0);
                  //  fitFcnV2_Cvetan->FixParameter(4,1);
                          

                         TFitResultPtr res = Yields_Central_1_CvetanPtBinned[ptbin]->Fit("fitFcnV2_Cvetan","USBMERI+","ep");
                   gStyle->SetOptStat("n");
                      gStyle->SetOptFit(1011);
                      TPaveStats *st = (TPaveStats*)Yields_Central_1_CvetanPtBinned[ptbin]->FindObject("stats");
                      st->SetX1NDC(0.8); //new x start position
                      st->SetY1NDC(0.8); //new x end position
                        Double_t par[5];
                    double chi2, edm, errdef;
                    int nvpar, nparx;
                        TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                    
                         fitFcnV2_Cvetan->SetChisquare(chi2);
                         int ndf = npfits-nvpar;
                         fitFcnV2_Cvetan->SetNDF(ndf);
                    
                    
                        fitFcnV2_Cvetan->GetParameters(par);
                         // improve the pictu
                       //   std::cout << "integral error: " << integralerror << std::endl;
                         fitFcnV2_Cvetan->Draw("same");
                         // draw the legend
                         TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                         legend->SetFillColorAlpha(kWhite, 0.);
                         legend->SetBorderSize(0);
                           legend->SetTextFont(42);
                         legend->SetTextSize(0.03);
                           Char_t message[80];
                           sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_Cvetan->GetChisquare(),fitFcnV2_Cvetan->GetNDF());
                           legend->AddEntry(fitFcnV2_Cvetan,message);
                    if(res->CovMatrixStatus() == 3){
                             //  sprintf(message,"The fit is a success");
                           }
                           else{
                               sprintf(message,"The fit is a failure");
                               legend->AddEntry(fitFcnV2_Cvetan,message);
                           }
                         legend->AddEntry(Yields_Central_1_CvetanPtBinned[ptbin],"Data","lpe");
                         legend->Draw();

                    V2_CvetanQuentin[ptbin] = par[3];
                    errV2_CvetanQuentin[ptbin] = fitFcnV2_Cvetan->GetParError(3);
                    F_CvetanQuentin[ptbin] = par[4];
                    errF_CvetanQuentin[ptbin] = fitFcnV2_Cvetan->GetParError(4);
                    
                   
                   sprintf(CanvasName,"%s/V2_Cvetan_PtBin%i.pdf",CanvasFolderName,ptbin);
                   cCvetanPtBinned->SaveAs(CanvasName);
                }
            
            
            // Fit Cvetan putting my constraints (F=1, V1 free, V0 not 1) - PtBinned
                
                TCanvas* cCvetanMePtBinned = new TCanvas;
            sprintf(hname, "Cvetan fit - My constraints F=1, v1 free, v0 free - PtBinned %d", ptbin);
                   cCvetanMePtBinned->SetTitle(hname);
                   cCvetanMePtBinned->Divide(1,1);
                   cCvetanMePtBinned->cd(1);
                   Yields_Central_1->DrawCopy();
            //       cCvetanMe->cd(2);
            //       Baseline_Central_1->DrawCopy();
            //       cCvetanMe->cd(3);
            //       Yields_Central_1_MinusBaseline->DrawCopy();
                
                
                {
                    cCvetanMePtBinned->cd(1);
                         TF1 *fitFcnV2_CvetanMe = new TF1("fitFcnV2_CvetanMe",CvetanFPtBinned,0,TMath::Pi(),5);//-TMath::Pi()/2,1.5*TMath::Pi(),5);
                         fitFcnV2_CvetanMe->SetNpx(500);
                         fitFcnV2_CvetanMe->SetLineWidth(4);
                         fitFcnV2_CvetanMe->SetLineColor(kBlue);
                         // first try without starting values for the parameters
                         // This defaults to 1 for each param.
                         // this results in an ok fit for the polynomial function
                         // however the non-linear part (lorenzian) does not
                         // respond well.
                          Double_t params[5] = {static_cast<Double_t>(ptbin),1,0,0.01,1};
                         fitFcnV2_CvetanMe->SetParameters(params);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetMaxIterations(10000);
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetPrecision();
                          TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->SetFCN(ChisquareCvetanFPtBinned);
                      //  gStyle->SetOptFit(1011);
                       //  histo->Fit("fitFcn","0");
                         // second try: set start values for some parameters

                        fitFcnV2_CvetanMe->SetParName(0,"ptbin");
                          fitFcnV2_CvetanMe->SetParName(1,"V0");
                            fitFcnV2_CvetanMe->SetParName(2,"V1");
                            fitFcnV2_CvetanMe->SetParName(3,"V2");
                          fitFcnV2_CvetanMe->SetParName(4,"F");
                    
                    fitFcnV2_CvetanMe->SetParLimits(4,0.1,50);
                    
                    fitFcnV2_CvetanMe->FixParameter(0,ptbin);
                 //   fitFcnV2_CvetanMe->FixParameter(1,1);
                   // fitFcnV2_CvetanMe->FixParameter(2,0);
                    fitFcnV2_CvetanMe->FixParameter(4,1);
                          

                         TFitResultPtr res = Yields_Central_1_CvetanMePtBinned[ptbin]->Fit("fitFcnV2_CvetanMe","USBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)Yields_Central_1_CvetanMePtBinned[ptbin]->FindObject("stats");
                    st->SetX1NDC(0.8); //new x start position
                    st->SetY1NDC(0.8); //new x end position
                    
                    double chi2, edm, errdef;
                    int nvpar, nparx;
                     TVirtualFitter::Fitter(Yields_Central_1_CvetanMePtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                     fitFcnV2_CvetanMe->SetChisquare(chi2);
                     int ndf = npfits-nvpar;
                     fitFcnV2_CvetanMe->SetNDF(ndf);
                    
                        Double_t par[5];
                        fitFcnV2_CvetanMe->GetParameters(par);
                         // improve the pictu
                       //   std::cout << "integral error: " << integralerror << std::endl;
                         fitFcnV2_CvetanMe->Draw("same");
                         // draw the legend
                         TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                         legend->SetFillColorAlpha(kWhite, 0.);
                         legend->SetBorderSize(0);
                           legend->SetTextFont(42);
                         legend->SetTextSize(0.03);
                           Char_t message[80];
                           sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_CvetanMe->GetChisquare(),fitFcnV2_CvetanMe->GetNDF());
                           legend->AddEntry(fitFcnV2_CvetanMe,message);
                    if(res->CovMatrixStatus() == 3){
                            //   sprintf(message,"The fit is a success");
                           }
                           else{
                               sprintf(message,"The fit is a failure");
                               legend->AddEntry(fitFcnV2_CvetanMe,message);
                           }
                         legend->AddEntry(Yields_Central_1_CvetanMePtBinned[ptbin],"Data","lpe");
                         legend->Draw();
                    
                    V2_CvetanQuentinMe[ptbin] = par[3]/par[1];
                    errV2_CvetanQuentinMe[ptbin] = (par[3]/(par[1])*sqrt(pow(fitFcnV2_CvetanMe->GetParError(3)/par[3],2)+pow(fitFcnV2_CvetanMe->GetParError(1)/par[1],2)));

                }
            
            
            // Fit ZYAM Yc = Bc + (Yp-Bp) + a0  + 2a1cos + 2a2cos  -  PtBinned
            
            
             TCanvas* cZYAMPtBinned = new TCanvas;
            sprintf(hname, "ZYAM fit - PtBinned %d", ptbin);
                cZYAMPtBinned->SetTitle(hname);
                cZYAMPtBinned->Divide(1,1);
                cZYAMPtBinned->cd(1);
                Yields_Central_1PtBinned[ptbin]->DrawCopy();
            //    cZYAM->cd(2);
            //    Baseline_Central_1->DrawCopy();
            //    cZYAM->cd(3);
            //    Yields_Central_1_MinusBaseline->DrawCopy();
            
            
            {
                cZYAMPtBinned->cd(1);
                     TF1 *fitFcnV2_ZYAM = new TF1("fitFcnV2_ZYAM",ZYAMPtBinned,0,TMath::Pi(),4);//-TMath::Pi()/2,1.5*TMath::Pi(),4);
                     fitFcnV2_ZYAM->SetNpx(500);
                     fitFcnV2_ZYAM->SetLineWidth(4);
                     fitFcnV2_ZYAM->SetLineColor(kCyan-7);
                     // first try without starting values for the parameters
                     // This defaults to 1 for each param.
                     // this results in an ok fit for the polynomial function
                     // however the non-linear part (lorenzian) does not
                     // respond well.
                      Double_t params[4] = {static_cast<Double_t>(ptbin),1,1,1};
                     fitFcnV2_ZYAM->SetParameters(params);
                      TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->SetMaxIterations(10000);
                      TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->SetPrecision();
                    TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->SetFCN(ChisquareZYAMPtBinned);
                    //gStyle->SetOptFit(1011);
                    //  histo->Fit("fitFcn","0");
                     // second try: set start values for some parameters

                    fitFcnV2_ZYAM->SetParName(0,"ptbin");
                      fitFcnV2_ZYAM->SetParName(1,"a0");
                        fitFcnV2_ZYAM->SetParName(2,"a1");
                        fitFcnV2_ZYAM->SetParName(3,"a2");
                    
                       fitFcnV2_ZYAM->FixParameter(0,ptbin);

                     TFitResultPtr res = Yields_Central_1_ZYAMPtBinned[ptbin]->Fit("fitFcnV2_ZYAM","USBMERI+","ep");
                    gStyle->SetOptStat("n");
                    gStyle->SetOptFit(1011);
                    TPaveStats *st = (TPaveStats*)Yields_Central_1_ZYAMPtBinned[ptbin]->FindObject("stats");
                    st->SetX1NDC(0.8); //new x start position
                    st->SetY1NDC(0.8); //new x end position
                
                double chi2, edm, errdef;
                int nvpar, nparx;
                 TVirtualFitter::Fitter(Yields_Central_1_ZYAMPtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                 fitFcnV2_ZYAM->SetChisquare(chi2);
                 int ndf = npfits-nvpar;
                 fitFcnV2_ZYAM->SetNDF(ndf);
                
                    Double_t par[4];
                    fitFcnV2_ZYAM->GetParameters(par);
                     // improve the pictu
                   //   std::cout << "integral error: " << integralerror << std::endl;
                     fitFcnV2_ZYAM->Draw("same");
                     // draw the legend
                     TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                     legend->SetFillColorAlpha(kWhite, 0.);
                     legend->SetBorderSize(0);
                       legend->SetTextFont(42);
                     legend->SetTextSize(0.03);
                       Char_t message[80];
                       sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_ZYAM->GetChisquare(),fitFcnV2_ZYAM->GetNDF());
                       legend->AddEntry(fitFcnV2_ZYAM,message);
                if(res->CovMatrixStatus() == 3){
                        //   sprintf(message,"The fit is a success");
                       }
                       else{
                           sprintf(message,"The fit is a failure");
                           legend->AddEntry(fitFcnV2_ZYAM,message);
                       }
                     legend->AddEntry(Yields_Central_1_ZYAMPtBinned[ptbin],"Data","lpe");
                     legend->Draw();
                
                
                V2_ZYAM[ptbin] = par[3]/(par[1]+baseline_centralPtBinned[ptbin]);
              //  errV2_ZYAM[ptbin] = (par[3]/(par[1] + baseline_centralPtBinned[ptbin])*sqrt(pow(fitFcnV2_ZYAM->GetParError(3)/par[3],2)+pow(fitFcnV2_ZYAM->GetParError(1)/par[1],2)+pow(errbaseline_centralPtBinned[ptbin]/baseline_centralPtBinned[ptbin],2)));
                
                errV2_ZYAM[ptbin] = (par[3]/(par[1] + baseline_centralPtBinned[ptbin])*sqrt(pow(fitFcnV2_ZYAM->GetParError(3)/par[3],2)+pow(sqrt(pow(fitFcnV2_ZYAM->GetParError(1),2)+pow(errbaseline_centralPtBinned[ptbin],2))/(par[1] + baseline_centralPtBinned[ptbin]),2)));

                sprintf(CanvasName,"%s/V2_ZYAM_PtBin%i.pdf",CanvasFolderName,ptbin);
                cZYAMPtBinned->SaveAs(CanvasName);
            }
            
            
             // Fit PRL Template - PtBinned
                
                TCanvas* cPRLTemplatePtBinned = new TCanvas;
                sprintf(hname, "PRL Template fit - PtBinned %d", ptbin);
                   cPRLTemplatePtBinned->SetTitle(hname);
                   cPRLTemplatePtBinned->Divide(1,1);
                   cPRLTemplatePtBinned->cd(1);
                   Yields_Central_1PtBinned[ptbin]->DrawCopy();
            //       cPRLTemplate->cd(2);
            //       Baseline_Central_1->DrawCopy();
            //       cPRLTemplate->cd(3);
            //       Yields_Central_1_MinusBaseline->DrawCopy();
                   
                   
                   {
                       cPRLTemplatePtBinned->cd(1);
                            TF1 *fitFcnV2_PRLTemplate = new TF1("fitFcnV2_PRLTemplate",PRLTemplatePtBinned,0,TMath::Pi(),3);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
                            fitFcnV2_PRLTemplate->SetNpx(500);
                            fitFcnV2_PRLTemplate->SetLineWidth(4);
                            fitFcnV2_PRLTemplate->SetLineColor(kRed+1);
                       TF1 *fitFcnV2_PRLTemplate_RidgeAndZero = new TF1("fitFcnV2_PRLTemplate_RidgeAndZero",PRLTemplate_RidgeAndZeroPtBinned,0,TMath::Pi(),3);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetNpx(500);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetLineWidth(2);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetLineStyle(9);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetLineColor(kRed);
                       TF1 *fitFcnV2_PRLTemplate_PeriphAndG = new TF1("fitFcnV2_PRLTemplate_PeriphAndG",PRLTemplate_PeriphAndGPtBinned,0,TMath::Pi(),3);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetNpx(500);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetLineWidth(2);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetLineStyle(9);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetLineColor(kBlack);
                            // first try without starting values for the parameters
                            // This defaults to 1 for each param.
                            // this results in an ok fit for the polynomial function
                            // however the non-linear part (lorenzian) does not
                            // respond well.
                             Double_t params[3] = {static_cast<Double_t>(ptbin),0.01,1};
                            fitFcnV2_PRLTemplate->SetParameters(params);
                             TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->SetMaxIterations(10000);
                             TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->SetPrecision();
                            TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->SetFCN(ChisquarePRLTemplatePtBinned);
                         //  gStyle->SetOptFit(1011);
                           //  histo->Fit("fitFcn","0");
                            // second try: set start values for some parameters

                            fitFcnV2_PRLTemplate->SetParName(0,"ptbin");
                             fitFcnV2_PRLTemplate->SetParName(1,"V2");
                               fitFcnV2_PRLTemplate->SetParName(2,"F");
                       

                       fitFcnV2_PRLTemplate->FixParameter(0,ptbin);
                       fitFcnV2_PRLTemplate->SetParLimits(2,0.1,50);
                             

                            TFitResultPtr res = Yields_Central_1_PRLTemplatePtBinned[ptbin]->Fit("fitFcnV2_PRLTemplate","USBMERI+","ep");
                       gStyle->SetOptStat("n");
                       gStyle->SetOptFit(1011);
                       TPaveStats *st = (TPaveStats*)Yields_Central_1_PRLTemplatePtBinned[ptbin]->FindObject("stats");
                       st->SetX1NDC(0.8); //new x start position
                       st->SetY1NDC(0.8); //new x end position
                       
                       double chi2, edm, errdef;
                       int nvpar, nparx;
                        TVirtualFitter::Fitter(Yields_Central_1_PRLTemplatePtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                        fitFcnV2_PRLTemplate->SetChisquare(chi2);
                        int ndf = npfits-nvpar;
                        fitFcnV2_PRLTemplate->SetNDF(ndf);
                       
                           Double_t par[3];
                           fitFcnV2_PRLTemplate->GetParameters(par);
                       fitFcnV2_PRLTemplate_RidgeAndZero->SetParameters(par);
                       fitFcnV2_PRLTemplate_PeriphAndG->SetParameters(par);
                            // improve the pictu
                          //   std::cout << "integral error: " << integralerror << std::endl;
                            fitFcnV2_PRLTemplate->Draw("same");
                       fitFcnV2_PRLTemplate_RidgeAndZero->Draw("same");
                       fitFcnV2_PRLTemplate_PeriphAndG->Draw("same");
                            // draw the legend
                            TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                            legend->SetFillColorAlpha(kWhite, 0.);
                            legend->SetBorderSize(0);
                              legend->SetTextFont(42);
                            legend->SetTextSize(0.03);
                              Char_t message[80];
                              sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate->GetChisquare(),fitFcnV2_PRLTemplate->GetNDF());
                              legend->AddEntry(fitFcnV2_PRLTemplate,message);
                       if(res->CovMatrixStatus() == 3){
                                 // sprintf(message,"The fit is a success");
                              }
                              else{
                                  sprintf(message,"The fit is a failure");
                                  legend->AddEntry(fitFcnV2_PRLTemplate,message);
                              }
                            legend->AddEntry(Yields_Central_1_PRLTemplatePtBinned[ptbin],"Data","lpe");
                           legend->AddEntry(fitFcnV2_PRLTemplate_RidgeAndZero,"F*Baseline_{Periph} + Ridge*G", "L");
                           legend->AddEntry(fitFcnV2_PRLTemplate_PeriphAndG,"F*Yield_{Periph} + G", "L");
                            legend->Draw();
                       
                       V2_PRL[ptbin] = par[1];
                       errV2_PRL[ptbin] = fitFcnV2_PRLTemplate->GetParError(1);
                       F_PRL[ptbin] = par[2];
                       errF_PRL[ptbin] = fitFcnV2_PRLTemplate->GetParError(2);

                       sprintf(CanvasName,"%s/V2_ATLAS_PtBin%i.pdf",CanvasFolderName,ptbin);
                       cPRLTemplatePtBinned->SaveAs(CanvasName);
                   }
            
            
             // Fit PRL Template + ZYAM Periph - PtBinned
                   
                   TCanvas* cPRLTemplate_PeriphZYAMPtBinned = new TCanvas;
            sprintf(hname, "PRL Template fit Periph ZYAM - PtBinned %d", ptbin);
                      cPRLTemplate_PeriphZYAMPtBinned->SetTitle(hname);
                      cPRLTemplate_PeriphZYAMPtBinned->Divide(1,1);
                      cPRLTemplate_PeriphZYAMPtBinned->cd(1);
                      Yields_Central_1PtBinned[ptbin]->DrawCopy();
            //          cPRLTemplate_PeriphZYAM->cd(2);
            //          Baseline_Central_1->DrawCopy();
            //          cPRLTemplate_PeriphZYAM->cd(3);
            //          Yields_Central_1_MinusBaseline->DrawCopy();
                      
                      
                      {
                          cPRLTemplate_PeriphZYAMPtBinned->cd(1);
                               TF1 *fitFcnV2_PRLTemplate_PeriphZYAM = new TF1("fitFcnV2_PRLTemplate_PeriphZYAM",PRLTemplate_PeriphZYAMPtBinned,0,TMath::Pi(),3);//-TMath::Pi()/2,1.5*TMath::Pi(),3);
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetNpx(500);
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetLineWidth(4);
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetLineColor(kBlack);
                               // first try without starting values for the parameters
                               // This defaults to 1 for each param.
                               // this results in an ok fit for the polynomial function
                               // however the non-linear part (lorenzian) does not
                               // respond well.
                                Double_t params[3] = {static_cast<Double_t>(ptbin),1,1};
                               fitFcnV2_PRLTemplate_PeriphZYAM->SetParameters(params);
                                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->SetMaxIterations(10000);
                                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->SetPrecision();
                                TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->SetFCN(ChisquarePRLTemplate_PeriphZYAMPtBinned);
                            //  gStyle->SetOptFit(1011);
                              //  histo->Fit("fitFcn","0");
                               // second try: set start values for some parameters

                                fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(0,"ptbin");
                                fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(1,"V2");
                                  fitFcnV2_PRLTemplate_PeriphZYAM->SetParName(2,"F");
                          
                            fitFcnV2_PRLTemplate_PeriphZYAM->FixParameter(0,ptbin);
                          fitFcnV2_PRLTemplate_PeriphZYAM->SetParLimits(2,0.1,5000);
                                

                               TFitResultPtr res = Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin]->Fit("fitFcnV2_PRLTemplate_PeriphZYAM","USBMERI+","ep");
                          gStyle->SetOptStat("n");
                          gStyle->SetOptFit(1011);
                          TPaveStats *st = (TPaveStats*)Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin]->FindObject("stats");
                          st->SetX1NDC(0.8); //new x start position
                          st->SetY1NDC(0.8); //new x end position
                          
                          double chi2, edm, errdef;
                          int nvpar, nparx;
                           TVirtualFitter::Fitter(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin])->GetStats(chi2,edm,errdef,nvpar,nparx);
                           fitFcnV2_PRLTemplate_PeriphZYAM->SetChisquare(chi2);
                           int ndf = npfits-nvpar;
                           fitFcnV2_PRLTemplate_PeriphZYAM->SetNDF(ndf);
                          
                              Double_t par[3];
                              fitFcnV2_PRLTemplate_PeriphZYAM->GetParameters(par);
                               // improve the pictu
                             //   std::cout << "integral error: " << integralerror << std::endl;
                               fitFcnV2_PRLTemplate_PeriphZYAM->Draw("same");
                               // draw the legend
                               TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
                               legend->SetFillColorAlpha(kWhite, 0.);
                               legend->SetBorderSize(0);
                                 legend->SetTextFont(42);
                               legend->SetTextSize(0.03);
                                 Char_t message[80];
                                 sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d",fitFcnV2_PRLTemplate_PeriphZYAM->GetChisquare(),fitFcnV2_PRLTemplate_PeriphZYAM->GetNDF());
                                 legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                          if(res->CovMatrixStatus() == 3){
                                   //  sprintf(message,"The fit is a success");
                                 }
                                 else{
                                     sprintf(message,"The fit is a failure");
                                     legend->AddEntry(fitFcnV2_PRLTemplate_PeriphZYAM,message);
                                 }
                               legend->AddEntry(Yields_Central_1_PRLTemplate_PeriphZYAMPtBinned[ptbin],"Data","lpe");
                               legend->Draw();
                          
                          V2_PRLPeriphZYAM[ptbin] = par[1];
                          errV2_PRLPeriphZYAM[ptbin] = fitFcnV2_PRLTemplate_PeriphZYAM->GetParError(1);
                          F_PRLPeriphZYAM[ptbin] = par[2];
                          errF_PRLPeriphZYAM[ptbin] = fitFcnV2_PRLTemplate_PeriphZYAM->GetParError(2);
                          
                          sprintf(CanvasName,"%s/V2_ATLASZYAM_PtBin%i.pdf",CanvasFolderName,ptbin);
                          cPRLTemplate_PeriphZYAMPtBinned->SaveAs(CanvasName);
                      }
            
            
            
            
            
            
            
    
        }
        
    }

    
    
    double PtMiddle[NbPtBins];
    memset( PtMiddle, 0, NbPtBins*sizeof(double));
    double PtErrorSize[NbPtBins];
    memset( PtErrorSize, 0, NbPtBins*sizeof(double));
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        PtMiddle[ptbin] = (PtBins[ptbin]+PtBins[ptbin+1])/2 - 0.02*4;
        PtErrorSize[ptbin] = (PtBins[ptbin+1]-PtBins[ptbin])/2;
    }
    
    TCanvas* cV2Pt = new TCanvas;
    cV2Pt->cd();
    
    TGraphErrors *grV2_wrt_Pt1 = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext1,PtErrorSize,errV2_Ext1);
             // TGraph *gr3 = new TGraph (n, K3, chi);
              grV2_wrt_Pt1->SetTitle("V_{2,J/#psi-tkl} wrt p_{T} for different extraction methods");
              grV2_wrt_Pt1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
              grV2_wrt_Pt1->GetYaxis()->SetTitle("V_{2,J/#psi-tkl}");
            grV2_wrt_Pt1->SetMarkerColor(kAzure-3);
            grV2_wrt_Pt1->SetLineColor(kAzure-3);
            grV2_wrt_Pt1->SetLineStyle(9);
            grV2_wrt_Pt1->SetMarkerStyle(4);
            grV2_wrt_Pt1->GetYaxis()->SetRangeUser(-0.004,0.014);
            grV2_wrt_Pt1->GetXaxis()->SetRangeUser(PtBins[0],PtBins[NbPtBins]);
    
              grV2_wrt_Pt1->Draw("AP");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_Pt1_noZYAM = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext1_noZYAM,PtErrorSize,errV2_Ext1_noZYAM);
             // TGraph *gr3 = new TGraph (n, K3, chi);
            grV2_wrt_Pt1_noZYAM->SetMarkerColor(kAzure-3);
            grV2_wrt_Pt1_noZYAM->SetLineColor(kAzure-3);
            grV2_wrt_Pt1_noZYAM->SetMarkerStyle(8);
              grV2_wrt_Pt1_noZYAM->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_Pt2 = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext2,PtErrorSize,errV2_Ext2);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_Pt2->SetMarkerColor(kGreen-7);
    grV2_wrt_Pt2->SetLineColor(kGreen-7);
    grV2_wrt_Pt2->SetLineStyle(9);
    grV2_wrt_Pt2->SetMarkerStyle(32);
      grV2_wrt_Pt2->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_Pt2_noZYAM = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext2_noZYAM,PtErrorSize,errV2_Ext2_noZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_Pt2_noZYAM->SetMarkerColor(kGreen-7);
    grV2_wrt_Pt2_noZYAM->SetLineColor(kGreen-7);
    grV2_wrt_Pt2_noZYAM->SetMarkerStyle(23);
      grV2_wrt_Pt2_noZYAM->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_Pt3 = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext3,PtErrorSize,errV2_Ext3);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_Pt3->SetMarkerColor(kMagenta-7);
    grV2_wrt_Pt3->SetLineColor(kMagenta-7);
    grV2_wrt_Pt3->SetLineStyle(9);
    grV2_wrt_Pt3->SetMarkerStyle(26);
      grV2_wrt_Pt3->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_Pt3_noZYAM = new TGraphErrors(NbPtBins,PtMiddle,V2_Ext3_noZYAM,PtErrorSize,errV2_Ext3_noZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_Pt3_noZYAM->SetMarkerColor(kMagenta-7);
    grV2_wrt_Pt3_noZYAM->SetLineColor(kMagenta-7);
    grV2_wrt_Pt3_noZYAM->SetMarkerStyle(22);
      grV2_wrt_Pt3_noZYAM->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_PtCvetanQuentin = new TGraphErrors(NbPtBins,PtMiddle,V2_CvetanQuentin,PtErrorSize,errV2_CvetanQuentin);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtCvetanQuentin->SetMarkerColor(kOrange+3);
    grV2_wrt_PtCvetanQuentin->SetLineColor(kOrange+3);
    grV2_wrt_PtCvetanQuentin->SetLineStyle(9);
    grV2_wrt_PtCvetanQuentin->SetMarkerStyle(28);
      grV2_wrt_PtCvetanQuentin->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }

//    TGraphErrors *grV2_wrt_PtCvetanQuentinMe = new TGraphErrors(NbPtBins,PtMiddle,V2_CvetanQuentinMe,PtErrorSize,errV2_CvetanQuentinMe);
//     // TGraph *gr3 = new TGraph (n, K3, chi);
//    grV2_wrt_PtCvetanQuentinMe->SetMarkerColor(6);
//    grV2_wrt_PtCvetanQuentinMe->SetLineColor(6);
//    grV2_wrt_PtCvetanQuentinMe->SetMarkerStyle(4);
//      grV2_wrt_PtCvetanQuentinMe->Draw("P");
//
//    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//    PtMiddle[ptbin] += 0.02;
//    }
    
    TGraphErrors *grV2_wrt_PtZYAM = new TGraphErrors(NbPtBins,PtMiddle,V2_ZYAM,PtErrorSize,errV2_ZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtZYAM->SetMarkerColor(kCyan-7);
    grV2_wrt_PtZYAM->SetLineColor(kCyan-7);
    grV2_wrt_PtZYAM->SetLineStyle(9);
    grV2_wrt_PtZYAM->SetMarkerStyle(27);
      grV2_wrt_PtZYAM->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_PtPRL = new TGraphErrors(NbPtBins,PtMiddle,V2_PRL,PtErrorSize,errV2_PRL);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtPRL->SetMarkerColor(kRed+1);
    grV2_wrt_PtPRL->SetLineColor(kRed+1);
    grV2_wrt_PtPRL->SetMarkerStyle(21);
      grV2_wrt_PtPRL->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_Pt2ATLAS = new TGraphErrors(NbPtBins,PtMiddle,V2_ATLASExt2,PtErrorSize,errV2_ATLASExt2);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_Pt2ATLAS->SetMarkerColor(kYellow);
    grV2_wrt_Pt2ATLAS->SetLineColor(kYellow);
    grV2_wrt_Pt2ATLAS->SetLineStyle(9);
    grV2_wrt_Pt2ATLAS->SetMarkerStyle(32);
      grV2_wrt_Pt2ATLAS->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.02;
    }
    
    TGraphErrors *grV2_wrt_PtPRLPeriphZYAM = new TGraphErrors(NbPtBins,PtMiddle,V2_PRLPeriphZYAM,PtErrorSize,errV2_PRLPeriphZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grV2_wrt_PtPRLPeriphZYAM->SetMarkerColor(kBlack);
    grV2_wrt_PtPRLPeriphZYAM->SetLineColor(kBlack);
    grV2_wrt_PtPRLPeriphZYAM->SetLineStyle(9);
    grV2_wrt_PtPRLPeriphZYAM->SetMarkerStyle(25);
      grV2_wrt_PtPRLPeriphZYAM->Draw("P");
    
    TLegend *legendo=new TLegend(0.12,0.60,0.40,0.90);
            legendo->SetFillColorAlpha(kWhite, 0.);
            legendo->SetBorderSize(0);
             legendo->SetTextFont(42);
             legendo->SetTextSize(0.02);
           legendo->AddEntry(grV2_wrt_Pt1,"V2_Ext1 (pPb-like)");
    legendo->AddEntry(grV2_wrt_Pt1_noZYAM,"V2_Ext1_noZYAM (pPb-like no ZYAM)");
            legendo->AddEntry(grV2_wrt_Pt2,"V2_Ext2 (pPb-like)");
    legendo->AddEntry(grV2_wrt_Pt2_noZYAM,"V2_Ext2_noZYAM (pPb-like no ZYAM)");
            legendo->AddEntry(grV2_wrt_Pt3,"V2_Ext3 (pPb-like)");
    legendo->AddEntry(grV2_wrt_Pt3_noZYAM,"V2_Ext3_noZYAM (pPb-like no ZYAM)");
            legendo->AddEntry(grV2_wrt_PtZYAM,"V2_ZYAM");
            legendo->AddEntry(grV2_wrt_PtCvetanQuentin,"V2_CvetanQuentin (Template + ZYAM Periph)");
          //  legendo->AddEntry(grV2_wrt_PtCvetanQuentinMe,"V2_CvetanQuentinMe");
            legendo->AddEntry(grV2_wrt_PtPRL,"V2_ATLAS Template");
    legendo->AddEntry(grV2_wrt_Pt2ATLAS,"V2_ATLAS Template 2");
            legendo->AddEntry(grV2_wrt_PtPRLPeriphZYAM,"V2_ATLAS Template and ZYAM)");
             legendo->Draw();
    
    cV2Pt->Update();
   // TLine *l=new TLine(cV2Pt->GetUxmin(),0.0,cV2Pt->GetUxmax(),0.0);
    TLine *l=new TLine(0.0,0.0,12.0,0.0);
    l->SetLineColor(kBlack);
    l->SetLineWidth(1);
    l->SetLineStyle(9);
    l->Draw();
    
    
    
    sprintf(CanvasName,"%s/V2wrtPt.pdf",CanvasFolderName);
    cV2Pt->SaveAs(CanvasName);
    
    double PtMiddleF[NbPtBins];
    memset( PtMiddleF, 0, NbPtBins*sizeof(double));
    double PtErrorSizeF[NbPtBins];
    memset( PtErrorSizeF, 0, NbPtBins*sizeof(double));
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        PtMiddleF[ptbin] = (PtBins[ptbin]+PtBins[ptbin+1])/2 - 0.02*1;
        PtErrorSizeF[ptbin] = (PtBins[ptbin+1]-PtBins[ptbin])/2;
    }
    
    TCanvas* cFPt = new TCanvas;
    cFPt->cd();
    
    TGraphErrors *grF_wrt_PtCvetanQuentin = new TGraphErrors(NbPtBins,PtMiddleF,F_CvetanQuentin,PtErrorSize,errF_CvetanQuentin);
             // TGraph *gr3 = new TGraph (n, K3, chi);
              grF_wrt_PtCvetanQuentin->SetTitle("F of V_{2,J/#psi-tkl} extraction wrt p_{T} for different extraction methods");
              grF_wrt_PtCvetanQuentin->GetXaxis()->SetTitle("p_{T} (GeV/c)");
              grF_wrt_PtCvetanQuentin->GetYaxis()->SetTitle("F");
            grF_wrt_PtCvetanQuentin->SetMarkerColor(kOrange+3);
            grF_wrt_PtCvetanQuentin->SetLineColor(kOrange+3);
            grF_wrt_PtCvetanQuentin->SetLineStyle(9);
            grF_wrt_PtCvetanQuentin->SetMarkerStyle(28);
            grF_wrt_PtCvetanQuentin->GetYaxis()->SetRangeUser(0.8,2.5);
            grF_wrt_PtCvetanQuentin->GetXaxis()->SetRangeUser(PtBins[0],PtBins[NbPtBins]);
              grF_wrt_PtCvetanQuentin->Draw("AP");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddleF[ptbin] += 0.02;
    }
    
    TGraphErrors *grF_wrt_PtPRL = new TGraphErrors(NbPtBins,PtMiddleF,F_PRL,PtErrorSize,errF_PRL);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grF_wrt_PtPRL->SetMarkerColor(kRed+1);
    grF_wrt_PtPRL->SetLineColor(kRed+1);
    grF_wrt_PtPRL->SetMarkerStyle(21);
      grF_wrt_PtPRL->Draw("P");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddleF[ptbin] += 0.02;
    }
    
    TGraphErrors *grF_wrt_PtPRLPeriphZYAM = new TGraphErrors(NbPtBins,PtMiddleF,F_PRLPeriphZYAM,PtErrorSize,errF_PRLPeriphZYAM);
     // TGraph *gr3 = new TGraph (n, K3, chi);
    grF_wrt_PtPRLPeriphZYAM->SetMarkerColor(kBlack);
    grF_wrt_PtPRLPeriphZYAM->SetLineColor(kBlack);
    grF_wrt_PtPRLPeriphZYAM->SetLineStyle(9);
    grF_wrt_PtPRLPeriphZYAM->SetMarkerStyle(25);
      grF_wrt_PtPRLPeriphZYAM->Draw("P");
    
    TLegend *legendf=new TLegend(0.12,0.75,0.40,0.90);
    legendf->SetFillColorAlpha(kWhite, 0.);
    legendf->SetBorderSize(0);
     legendf->SetTextFont(42);
     legendf->SetTextSize(0.02);
            legendf->AddEntry(grF_wrt_PtCvetanQuentin,"F_CvetanQuentin (Template + ZYAM Periph)");
            legendf->AddEntry(grF_wrt_PtPRL,"F_ATLAS Template");
            legendf->AddEntry(grF_wrt_PtPRLPeriphZYAM,"F_ATLAS Template and ZYAM");
             legendf->Draw();
    
    cFPt->Update();
    //TLine *lF=new TLine(cFPt->GetUxmin(),1.0,cFPt->GetUxmax(),1.0);
    TLine *lF=new TLine(0.,1.0,12.0,1.0);
    lF->SetLineColor(kBlack);
    lF->SetLineWidth(1);
    lF->SetLineStyle(9);
    lF->Draw();
    
    sprintf(CanvasName,"%s/FwrtPt.pdf",CanvasFolderName);
    cFPt->SaveAs(CanvasName);
    
    CSVRow row;
    
    Char_t RadicalTKLName[500];
    Char_t FolderCSVTKLName[500];
    Char_t SystematicsFileTKLName[500];
    
    cout << "Will get tkl csv"<<endl;
    cout << "radical " <<radical <<endl;
   
   sscanf(radical, "NewAnalysisAllEstCentBvr_%[^_]_%[^_]_%[^-]-%[^_]_%[^-]-%[^_]_pt%[^-]-%s", dataUsed, estimator, mostCentral,lessCentral, lessPeriph, mostPeriph, minpt, rest);
    cout << "dataUsed " <<dataUsed <<endl;
    cout << "estimator " <<estimator <<endl;
    cout << "mostCentral " <<mostCentral <<endl;
    cout << "lessCentral " <<lessCentral <<endl;
    cout << "lessPeriph " <<lessPeriph <<endl;
    cout << "mostPeriph " <<mostPeriph <<endl;
    cout << "minpt " <<minpt <<endl;
    cout << "rest " <<rest <<endl;
    
    sprintf(RadicalTKLName, "NewAnalysisAllEst_TKL_16h_%s_0-5_40-100_pt0-12", estimator);
    cout<<"RadicalTKLName"<<RadicalTKLName<<endl;
    sprintf(FolderCSVTKLName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrier/%s", RadicalTKLName);
    cout<<"FolderCSVTKLName"<<FolderCSVTKLName<<endl;
    bool useTklAdapt = kFALSE;
    if(isCentralityStudy && useTklAdapt){
  //  if(kTRUE){
        sprintf(FolderCSVTKLName,"/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrier/NewAnalysisAllEst_TKL_16h_%s_%i-%i_%i-%i_pt0-12", estimator, stoi(mostCentral),stoi(lessCentral), stoi(lessPeriph), stoi(mostPeriph));
        cout<<"FolderCSVTKLName"<<FolderCSVTKLName<<endl;
    }
    
    
    sprintf(SystematicsFileTKLName,"%s/SystematicsFile.csv", FolderCSVTKLName);
    cout<<"SystematicsFileTKLName"<<SystematicsFileTKLName<<endl;
    std::ifstream filecsv(SystematicsFileTKLName);
    
    double v2ClassiqueTKL;
    double errv2ClassiqueTKL;
    double v2ClassiqueTKL_noZYAM;
    double errv2ClassiqueTKL_noZYAM;
    double v2CvetanQuentinTKL;
    double errv2CvetanQuentinTKL;
    double v2ZYAMTKL;
    double errv2ZYAMTKL;
    double v2PRLTKL;
    double errv2PRLTKL;
    double v2PRLPeriphZYAMTKL;
    double errv2PRLPeriphZYAMTKL;
    
    while(filecsv >> row)
    {
         v2ClassiqueTKL = stod(std::string(row[11]));
         errv2ClassiqueTKL = stod(std::string(row[12]));
         v2ClassiqueTKL_noZYAM = stod(std::string(row[13]));
         errv2ClassiqueTKL_noZYAM = stod(std::string(row[14]));
         v2CvetanQuentinTKL = stod(std::string(row[15]));
         errv2CvetanQuentinTKL = stod(std::string(row[16]));
         v2ZYAMTKL = stod(std::string(row[17]));
         errv2ZYAMTKL = stod(std::string(row[18]));
         v2PRLTKL = stod(std::string(row[19]));
         errv2PRLTKL = stod(std::string(row[20]));
         v2PRLPeriphZYAMTKL = stod(std::string(row[21]));
         errv2PRLPeriphZYAMTKL = stod(std::string(row[22]));
    }
    
    
    // Prendre les résultats sur V2JPsi-Tkl et utiliser ceux sur v2tkl pour trouver le v2JPsi pour chaque méthode
    
    //Résultats TKL 0-5% 40-90% - Baseline central 0 Yp unc. NOT PROPAGATED
    
//            double v2ClassiqueTKL = 0.065334;
//            double errv2ClassiqueTKL = 0.000541066;
//            double v2CvetanQuentinTKL = 0.0662997;
//            double errv2CvetanQuentinTKL = 0.000409345;
//            double v2CvetanQuentinMeTKL = 0.0654186;
//            double errv2CvetanQuentinMeTKL = 0.000414581;
//            double v2ZYAMTKL = 0.0653334;
//            double errv2ZYAMTKL = 0.00104418;
//            double v2PRLTKL = 0.131845;
//            double errv2PRLTKL = 0.00186927;
//            double v2PRLPeriphZYAMTKL = 0.0673054;
//            double errv2PRLPeriphZYAMTKL = 0.000405759;
    
    //Résultats TKL 0-5% 40-90% - Baseline central 0 Yp unc.  PROPAGATED
       
//               double v2ClassiqueTKL = 0.065334;
//               double errv2ClassiqueTKL = 0.000541066;
//               double v2CvetanQuentinTKL = 0.0663438;
//               double errv2CvetanQuentinTKL = 0.000677986;
//               double v2CvetanQuentinMeTKL = 0.065326;
//               double errv2CvetanQuentinMeTKL = 0.000688108;
//               double v2ZYAMTKL = 0.0652409;
//               double errv2ZYAMTKL = 0.0017316;
//               double v2PRLTKL = 0.131677;
//               double errv2PRLTKL = 0.00237045;
//               double v2PRLPeriphZYAMTKL = 0.0672543;
//               double errv2PRLPeriphZYAMTKL = 0.000672808;
    
    //Résultats TKL 0-5% 40-90% - Baseline central min Yp unc. PROPAGATED
    
//            double v2ClassiqueTKL = 0.065334;
//            double errv2ClassiqueTKL = 0.000541066;
//            double v2CvetanQuentinTKL = 0.067723;
//            double errv2CvetanQuentinTKL = 0.000668852;
//            double v2CvetanQuentinMeTKL = 0.0655555;
//            double errv2CvetanQuentinMeTKL = 0.000690527;
//            double v2ZYAMTKL = 0.0652409;
//            double errv2ZYAMTKL = 0.000812056;
//            double v2PRLTKL = 0.131677;
//            double errv2PRLTKL = 0.00237045;
//            double v2PRLPeriphZYAMTKL = 0.0672543;
//            double errv2PRLPeriphZYAMTKL = 0.000672808;
    
    
    //Résultats TKL 0-5% 40-100% - Baseline central min Yp unc. PROPAGATED FitFile_NewAnalysis_16hjo10_TKL_0-5_40-100_pt0-12.root
        
//            double v2ClassiqueTKL = 0.0612282;
//            double errv2ClassiqueTKL = 0.000638029;
//            double v2ClassiqueTKL_noZYAM = 0.0702399;
//            double errv2ClassiqueTKL_noZYAM = 0.00144388;
//            double v2CvetanQuentinTKL = 0.0622471;
//            double errv2CvetanQuentinTKL = 0.00187047;
//            double v2CvetanQuentinMeTKL = 0.0611412;
//            double errv2CvetanQuentinMeTKL = 0.000849029;
//            double v2ZYAMTKL = 0.0611414;
//            double errv2ZYAMTKL = 0.00129501;
//            double v2PRLTKL = 0.108204;
//            double errv2PRLTKL = 0.00371707;
//            double v2PRLPeriphZYAMTKL = 0.0628502;
//            double errv2PRLPeriphZYAMTKL = 0.00230936;
    
    //Résultats TKL 0-5% 40-100% - Baseline central min Yp unc. NewAnalysisAllEst_TKL_16h_V0MPercentile_0-5_40-100_pt0-12_SummationZvtxMethod1c
    
//        double v2ClassiqueTKL = 0.0524698;
//        double errv2ClassiqueTKL = 0.000302415;
//        double v2ClassiqueTKL_noZYAM = 0.0617702;
//        double errv2ClassiqueTKL_noZYAM = 0.000355983;
//        double v2CvetanQuentinTKL = 0.0509369;
//        double errv2CvetanQuentinTKL = 0.000464702;
//        double v2CvetanQuentinMeTKL = 0.0611412; //no
//        double errv2CvetanQuentinMeTKL = 0.000849029; //no
//        double v2ZYAMTKL = 0.0523923;
//        double errv2ZYAMTKL = 0.000424847;
//        double v2PRLTKL = 0.0636312;
//        double errv2PRLTKL = 0.000484392;
//        double v2PRLPeriphZYAMTKL = 0.0504003;
//        double errv2PRLPeriphZYAMTKL = 0.000501764;
    

    //
    
    double v2_Ext1[NbPtBins];
    memset( v2_Ext1, 0, NbPtBins*sizeof(double));
    double v2_Ext2[NbPtBins];
    memset( v2_Ext2, 0, NbPtBins*sizeof(double));
    double v2_Ext3[NbPtBins];
    memset( v2_Ext3, 0, NbPtBins*sizeof(double));
    double v2_Ext1_noZYAM[NbPtBins];
    memset( v2_Ext1_noZYAM, 0, NbPtBins*sizeof(double));
    double v2_Ext2_noZYAM[NbPtBins];
    memset( v2_Ext2_noZYAM, 0, NbPtBins*sizeof(double));
    double v2_Ext3_noZYAM[NbPtBins];
    memset( v2_Ext3_noZYAM, 0, NbPtBins*sizeof(double));
    double v2_CvetanQuentin[NbPtBins];
    memset( v2_CvetanQuentin, 0, NbPtBins*sizeof(double));
    double v2_CvetanQuentinMe[NbPtBins];
    memset( v2_CvetanQuentinMe, 0, NbPtBins*sizeof(double));
    double v2_ZYAM[NbPtBins];
    memset( v2_ZYAM, 0, NbPtBins*sizeof(double));
    double v2_PRL[NbPtBins];
    memset( v2_PRL, 0, NbPtBins*sizeof(double));
    double v2_ATLASExt2[NbPtBins];
    memset( v2_ATLASExt2, 0, NbPtBins*sizeof(double));
    double v2_PRLPeriphZYAM[NbPtBins];
    memset( v2_PRLPeriphZYAM, 0, NbPtBins*sizeof(double));
    double v2_Ext1_V3[NbPtBins];
    memset( v2_Ext1_V3, 0, NbPtBins*sizeof(double));
    double errv2_Ext1[NbPtBins];
    memset( errv2_Ext1, 0, NbPtBins*sizeof(double));
    double errv2_Ext2[NbPtBins];
    memset( errv2_Ext2, 0, NbPtBins*sizeof(double));
    double errv2_Ext3[NbPtBins];
    memset( errv2_Ext3, 0, NbPtBins*sizeof(double));
    double errv2_Ext1_noZYAM[NbPtBins];
    memset( errv2_Ext1_noZYAM, 0, NbPtBins*sizeof(double));
    double errv2_Ext2_noZYAM[NbPtBins];
    memset( errv2_Ext2_noZYAM, 0, NbPtBins*sizeof(double));
    double errv2_Ext3_noZYAM[NbPtBins];
    memset( errv2_Ext3_noZYAM, 0, NbPtBins*sizeof(double));
    double errv2_CvetanQuentin[NbPtBins];
    memset( errv2_CvetanQuentin, 0, NbPtBins*sizeof(double));
    double errv2_CvetanQuentinMe[NbPtBins];
    memset( errv2_CvetanQuentinMe, 0, NbPtBins*sizeof(double));
    double errv2_ZYAM[NbPtBins];
    memset( errv2_ZYAM, 0, NbPtBins*sizeof(double));
    double errv2_PRL[NbPtBins];
    memset( errv2_PRL, 0, NbPtBins*sizeof(double));
    double errv2_ATLASExt2[NbPtBins];
    memset( errv2_ATLASExt2, 0, NbPtBins*sizeof(double));
    double errv2_PRLPeriphZYAM[NbPtBins];
    memset( errv2_PRLPeriphZYAM, 0, NbPtBins*sizeof(double));
    double errv2_Ext1_V3[NbPtBins];
    memset( errv2_Ext1_V3, 0, NbPtBins*sizeof(double));
    
//    double v2_PbPb[NbPtBins-1] = {0.01,0.06,0.075,0.05};
//    double v2_pPb[NbPtBins] = {0,-0.015,0.035,0.08,0.04};
//    double errv2_PbPb[NbPtBins-1] = {0.015,0.015,0.02,0.03};
//    double errv2_pPb[NbPtBins] = {0.015,0.02,0.02,0.015,0.04};
//    double PbPbMiddle[NbPtBins-1] = {1,3,5,7};
//    double PbPbPtErrorSize[NbPtBins-1] = {1,1,1,1};
//    double pPbMiddle[NbPtBins] = {1,2.5,3.5,5,7};
//    double pPbPtErrorSize[NbPtBins] = {1,0.5,0.5,1,1};
    
    for(int pt_idx=0; pt_idx < NbPtBins; pt_idx++){
        v2_Ext1[pt_idx] = V2_Ext1[pt_idx]/v2ClassiqueTKL;
        v2_Ext2[pt_idx] = V2_Ext2[pt_idx]/v2ClassiqueTKL;
        v2_Ext3[pt_idx] = V2_Ext3[pt_idx]/v2ClassiqueTKL;
        v2_Ext1_noZYAM[pt_idx] = V2_Ext1_noZYAM[pt_idx]/v2ClassiqueTKL_noZYAM;
        v2_Ext2_noZYAM[pt_idx] = V2_Ext2_noZYAM[pt_idx]/v2ClassiqueTKL_noZYAM;
        v2_Ext3_noZYAM[pt_idx] = V2_Ext3_noZYAM[pt_idx]/v2ClassiqueTKL_noZYAM;
        v2_CvetanQuentin[pt_idx] = V2_CvetanQuentin[pt_idx]/v2CvetanQuentinTKL;
       // v2_CvetanQuentinMe[pt_idx] = V2_CvetanQuentinMe[pt_idx]/v2CvetanQuentinMeTKL;
        v2_ZYAM[pt_idx] = V2_ZYAM[pt_idx]/v2ZYAMTKL;
        v2_PRL[pt_idx] = V2_PRL[pt_idx]/v2PRLTKL;
        v2_ATLASExt2[pt_idx] = V2_ATLASExt2[pt_idx]/v2PRLTKL;
        v2_PRLPeriphZYAM[pt_idx] = V2_PRLPeriphZYAM[pt_idx]/v2PRLPeriphZYAMTKL;
        v2_Ext1_V3[pt_idx] = V2_Ext1_V3[pt_idx]/v2ClassiqueTKL;
    }
    
    for(int pt_idx=0; pt_idx < NbPtBins; pt_idx++){
        errv2_Ext1[pt_idx] = abs(sqrt(pow(errV2_Ext1[pt_idx]/V2_Ext1[pt_idx],2)+pow(errv2ClassiqueTKL/v2ClassiqueTKL,2))*v2_Ext1[pt_idx]);
        errv2_Ext2[pt_idx] = abs(sqrt(pow(errV2_Ext2[pt_idx]/V2_Ext2[pt_idx],2)+pow(errv2ClassiqueTKL/v2ClassiqueTKL,2))*v2_Ext2[pt_idx]);
        errv2_Ext3[pt_idx] = abs(sqrt(pow(errV2_Ext3[pt_idx]/V2_Ext3[pt_idx],2)+pow(errv2ClassiqueTKL/v2ClassiqueTKL,2))*v2_Ext3[pt_idx]);
        errv2_Ext1_noZYAM[pt_idx] = abs(sqrt(pow(errV2_Ext1_noZYAM[pt_idx]/V2_Ext1_noZYAM[pt_idx],2)+pow(errv2ClassiqueTKL_noZYAM/v2ClassiqueTKL_noZYAM,2))*v2_Ext1_noZYAM[pt_idx]);
        errv2_Ext2_noZYAM[pt_idx] = abs(sqrt(pow(errV2_Ext2_noZYAM[pt_idx]/V2_Ext2_noZYAM[pt_idx],2)+pow(errv2ClassiqueTKL_noZYAM/v2ClassiqueTKL_noZYAM,2))*v2_Ext2_noZYAM[pt_idx]);
        errv2_Ext3_noZYAM[pt_idx] = abs(sqrt(pow(errV2_Ext3_noZYAM[pt_idx]/V2_Ext3_noZYAM[pt_idx],2)+pow(errv2ClassiqueTKL_noZYAM/v2ClassiqueTKL_noZYAM,2))*v2_Ext3_noZYAM[pt_idx]);
        errv2_CvetanQuentin[pt_idx] = abs(sqrt(pow(errV2_CvetanQuentin[pt_idx]/V2_CvetanQuentin[pt_idx],2)+pow(errv2CvetanQuentinTKL/v2CvetanQuentinTKL,2))*v2_CvetanQuentin[pt_idx]);
      //  errv2_CvetanQuentinMe[pt_idx] = abs(sqrt(pow(errV2_CvetanQuentinMe[pt_idx]/V2_CvetanQuentinMe[pt_idx],2)+pow(errv2CvetanQuentinMeTKL/v2CvetanQuentinMeTKL,2))*v2_CvetanQuentinMe[pt_idx]);
        errv2_ZYAM[pt_idx] = abs(sqrt(pow(errV2_ZYAM[pt_idx]/V2_ZYAM[pt_idx],2)+pow(errv2ZYAMTKL/v2ZYAMTKL,2))*v2_ZYAM[pt_idx]);
        errv2_PRL[pt_idx] = abs(sqrt(pow(errV2_PRL[pt_idx]/V2_PRL[pt_idx],2)+pow(errv2PRLTKL/v2PRLTKL,2))*v2_PRL[pt_idx]);
        errv2_ATLASExt2[pt_idx] = abs(sqrt(pow(errV2_ATLASExt2[pt_idx]/V2_ATLASExt2[pt_idx],2)+pow(errv2PRLTKL/v2PRLTKL,2))*v2_PRL[pt_idx]);
        errv2_PRLPeriphZYAM[pt_idx] = abs(sqrt(pow(errV2_PRLPeriphZYAM[pt_idx]/V2_PRLPeriphZYAM[pt_idx],2)+pow(errv2PRLPeriphZYAMTKL/v2PRLPeriphZYAMTKL,2))*v2_PRLPeriphZYAM[pt_idx]);
        errv2_Ext1_V3[pt_idx] = abs(sqrt(pow(errV2_Ext1_V3[pt_idx]/V2_Ext1_V3[pt_idx],2)+pow(errv2ClassiqueTKL/v2ClassiqueTKL,2))*v2_Ext1_V3[pt_idx]);
    }
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
        PtMiddle[ptbin] = (PtBins[ptbin]+PtBins[ptbin+1])/2 - 0.1*1;
        PtErrorSize[ptbin] = 0;//(PtBins[ptbin+1]-PtBins[ptbin])/2;
     }
     
     TCanvas* cv2Pt = new TCanvas;
     cv2Pt->cd();
    
    TGraphErrors *grv2_wrt_PtPRL = new TGraphErrors(NbPtBins,PtMiddle,v2_PRL,PtErrorSize,errv2_PRL);
        // TGraph *gr3 = new TGraph (n, K3, chi);
    grv2_wrt_PtPRL->SetTitle("#it{v}_{2,J/#psi} wrt #it{p}_{T,J/#psi} for different extraction methods");
    grv2_wrt_PtPRL->GetXaxis()->SetTitle("#it{p}_{T,J/#psi} (GeV/c)");
    grv2_wrt_PtPRL->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
        grv2_wrt_PtPRL->SetMarkerColor(kRed+1);
          grv2_wrt_PtPRL->SetLineColor(kRed+1);
      grv2_wrt_PtPRL->SetLineWidth(2);
          grv2_wrt_PtPRL->SetMarkerStyle(21);
    grv2_wrt_PtPRL->GetYaxis()->SetRangeUser(-0.10,0.15);
    grv2_wrt_PtPRL->GetXaxis()->SetRangeUser(PtBins[0],15);
            grv2_wrt_PtPRL->Draw("AP");
    
    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
    PtMiddle[ptbin] += 0.1;
    }
     
     TGraphErrors *grv2_wrt_Pt1 = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext1,PtErrorSize,errv2_Ext1);
              // TGraph *gr3 = new TGraph (n, K3, chi);
//               grv2_wrt_Pt1->SetTitle("#it{v}_{2,J/#psi} wrt #it{p}_{T,J/#psi} for different extraction methods");
//               grv2_wrt_Pt1->GetXaxis()->SetTitle("#it{p}_{T,J/#psi} (GeV/c)");
//               grv2_wrt_Pt1->GetYaxis()->SetTitle("#it{v}_{2,J/#psi}");
             grv2_wrt_Pt1->SetMarkerColor(kAzure-3);
             grv2_wrt_Pt1->SetLineColor(kAzure-3);
            // grv2_wrt_Pt1->SetLineStyle(9);
                grv2_wrt_Pt1->SetLineWidth(2);
             grv2_wrt_Pt1->SetMarkerStyle(8);
//             grv2_wrt_Pt1->GetYaxis()->SetRangeUser(-0.10,0.15);
//             grv2_wrt_Pt1->GetXaxis()->SetRangeUser(PtBins[0],15);
             grv2_wrt_Pt1->Draw("P");

     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.1;
     }
    
//    TGraphErrors *grv2_wrt_Pt1_noZYAM = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext1_noZYAM,PtErrorSize,errv2_Ext1_noZYAM);
//     // TGraph *gr3 = new TGraph (n, K3, chi);
//    grv2_wrt_Pt1_noZYAM->SetMarkerColor(kAzure-3);
//       grv2_wrt_Pt1_noZYAM->SetLineColor(kAzure-3);
//       grv2_wrt_Pt1_noZYAM->SetMarkerStyle(8);
//         grv2_wrt_Pt1_noZYAM->Draw("P");
//
//    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//    PtMiddle[ptbin] += 0.02;
//    }
     
     TGraphErrors *grv2_wrt_Pt2 = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext2,PtErrorSize,errv2_Ext2);
      // TGraph *gr3 = new TGraph (n, K3, chi);
      grv2_wrt_Pt2->SetMarkerColor(kGreen+2);
        grv2_wrt_Pt2->SetLineColor(kGreen+2);
    grv2_wrt_Pt2->SetLineWidth(2);
       // grv2_wrt_Pt2->SetLineStyle(9);
        grv2_wrt_Pt2->SetMarkerStyle(33);
          grv2_wrt_Pt2->Draw("P");

     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
     PtMiddle[ptbin] += 0.1;
     }
    
//    TGraphErrors *grv2_wrt_Pt2_noZYAM = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext2_noZYAM,PtErrorSize,errv2_Ext2_noZYAM);
//     // TGraph *gr3 = new TGraph (n, K3, chi);
//    grv2_wrt_Pt2_noZYAM->SetMarkerColor(kGreen-7);
//    grv2_wrt_Pt2_noZYAM->SetLineColor(kGreen-7);
//    grv2_wrt_Pt2_noZYAM->SetMarkerStyle(23);
//      grv2_wrt_Pt2_noZYAM->Draw("P");
//
//    for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//    PtMiddle[ptbin] += 0.02;
//    }
    
//    TGraphErrors *grv2_wrt_Pt3 = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext3,PtErrorSize,errv2_Ext3);
//         // TGraph *gr3 = new TGraph (n, K3, chi);
//        grv2_wrt_Pt3->SetMarkerColor(kMagenta-7);
//        grv2_wrt_Pt3->SetLineColor(kMagenta-7);
//        grv2_wrt_Pt3->SetLineStyle(9);
//        grv2_wrt_Pt3->SetMarkerStyle(26);
//          grv2_wrt_Pt3->Draw("P");
//
//        for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//        PtMiddle[ptbin] += 0.02;
//        }
    
//    TGraphErrors *grv2_wrt_Pt3_noZYAM = new TGraphErrors(NbPtBins,PtMiddle,v2_Ext3_noZYAM,PtErrorSize,errv2_Ext3_noZYAM);
//            // TGraph *gr3 = new TGraph (n, K3, chi);
//           grv2_wrt_Pt3_noZYAM->SetMarkerColor(kMagenta-7);
//           grv2_wrt_Pt3_noZYAM->SetLineColor(kMagenta-7);
//           grv2_wrt_Pt3_noZYAM->SetMarkerStyle(22);
//             grv2_wrt_Pt3_noZYAM->Draw("P");
//
//           for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//           PtMiddle[ptbin] += 0.02;
//           }
//     
//     TGraphErrors *grv2_wrt_PtCvetanQuentin = new TGraphErrors(NbPtBins,PtMiddle,v2_CvetanQuentin,PtErrorSize,errv2_CvetanQuentin);
//      // TGraph *gr3 = new TGraph (n, K3, chi);
//    grv2_wrt_PtCvetanQuentin->SetMarkerColor(kOrange+3);
//     grv2_wrt_PtCvetanQuentin->SetLineColor(kOrange+3);
//     grv2_wrt_PtCvetanQuentin->SetLineStyle(9);
//     grv2_wrt_PtCvetanQuentin->SetMarkerStyle(28);
//       grv2_wrt_PtCvetanQuentin->Draw("P");
//
//     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//     PtMiddle[ptbin] += 0.02;
//     }

//     TGraphErrors *grv2_wrt_PtCvetanQuentinMe = new TGraphErrors(NbPtBins,PtMiddle,v2_CvetanQuentinMe,PtErrorSize,errv2_CvetanQuentinMe);
//      // TGraph *gr3 = new TGraph (n, K3, chi);
//     grv2_wrt_PtCvetanQuentinMe->SetMarkerColor(6);
//     grv2_wrt_PtCvetanQuentinMe->SetLineColor(6);
//     grv2_wrt_PtCvetanQuentinMe->SetMarkerStyle(4);
//       grv2_wrt_PtCvetanQuentinMe->Draw("P");
//
//     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//     PtMiddle[ptbin] += 0.02;
//     }

//     TGraphErrors *grv2_wrt_PtZYAM = new TGraphErrors(NbPtBins,PtMiddle,v2_ZYAM,PtErrorSize,errv2_ZYAM);
//      // TGraph *gr3 = new TGraph (n, K3, chi);
//     grv2_wrt_PtZYAM->SetMarkerColor(kCyan-7);
//     grv2_wrt_PtZYAM->SetLineColor(kCyan-7);
//     grv2_wrt_PtZYAM->SetLineStyle(9);
//     grv2_wrt_PtZYAM->SetMarkerStyle(27);
//       grv2_wrt_PtZYAM->Draw("P");
//
//     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//     PtMiddle[ptbin] += 0.02;
//     }
     
//     TGraphErrors *grv2_wrt_PtPRL = new TGraphErrors(NbPtBins,PtMiddle,v2_PRL,PtErrorSize,errv2_PRL);
//      // TGraph *gr3 = new TGraph (n, K3, chi);
//  //  grv2_wrt_PtPRL->SetTitle("v2 JPsi wrt Pt for different systems");
//                  // grv2_wrt_PtPRL->GetXaxis()->SetTitle("Pt (GeV/c)");
//                  // grv2_wrt_PtPRL->GetYaxis()->SetTitle("v2 (JPsi)");
//      grv2_wrt_PtPRL->SetMarkerColor(kRed+1);
//        grv2_wrt_PtPRL->SetLineColor(kRed+1);
//    grv2_wrt_PtPRL->SetLineWidth(2);
//        grv2_wrt_PtPRL->SetMarkerStyle(21);
//          grv2_wrt_PtPRL->Draw("P");
    
//    TGraphErrors *grv2_wrt_PtPbPb = new TGraphErrors(NbPtBins-1,PbPbMiddle,v2_PbPb,PbPbPtErrorSize,errv2_PbPb);
//        // TGraph *gr3 = new TGraph (n, K3, chi);
//       grv2_wrt_PtPbPb->SetMarkerColor(kBlack);
//       grv2_wrt_PtPbPb->SetLineColor(kBlack);
//       grv2_wrt_PtPbPb->SetMarkerStyle(26);
//         grv2_wrt_PtPbPb->Draw("P");
//
//    TGraphErrors *grv2_wrt_PtpPb = new TGraphErrors(NbPtBins,pPbMiddle,v2_pPb,pPbPtErrorSize,errv2_pPb);
//     // TGraph *gr3 = new TGraph (n, K3, chi);
//    grv2_wrt_PtpPb->SetMarkerColor(kBlue);
//    grv2_wrt_PtpPb->SetLineColor(kBlue);
//    grv2_wrt_PtpPb->SetMarkerStyle(26);
//      grv2_wrt_PtpPb->Draw("P");
     
//     for(int ptbin=0;ptbin<NbPtBins;ptbin++){
//     PtMiddle[ptbin] += 0.02;
//     }
//
//     TGraphErrors *grv2_wrt_PtPRLPeriphZYAM = new TGraphErrors(NbPtBins,PtMiddle,v2_PRLPeriphZYAM,PtErrorSize,errv2_PRLPeriphZYAM);
//      // TGraph *gr3 = new TGraph (n, K3, chi);
//    grv2_wrt_PtPRLPeriphZYAM->SetMarkerColor(kBlack);
//     grv2_wrt_PtPRLPeriphZYAM->SetLineColor(kBlack);
//     grv2_wrt_PtPRLPeriphZYAM->SetLineStyle(9);
//     grv2_wrt_PtPRLPeriphZYAM->SetMarkerStyle(25);
//       grv2_wrt_PtPRLPeriphZYAM->Draw("P");
//
     TLegend *legendov=new TLegend(0.12,0.70,0.40,0.90);
              legendov->SetFillColorAlpha(kWhite, 0.);
              legendov->SetBorderSize(0);
               legendov->SetTextFont(42);
               legendov->SetTextSize(0.022);
    Char_t messagi[80];
    sprintf(messagi,"#it{v}_{2,J/#psi} Extraction methods:");
    legendov->AddEntry(grv2_wrt_Pt1,messagi,"");
            legendov->AddEntry(grv2_wrt_Pt1,"Subtracted yields default: Y_{C}-Y_{P}=a_{0}+2a_{1}cos(#Delta#phi)+2a_{2}cos(2#Delta#phi)");
  //  legendov->AddEntry(grv2_wrt_Pt1_noZYAM,"v2_Ext1_noZYAM (pPb-like no ZYAM)");
             legendov->AddEntry(grv2_wrt_Pt2,"Subtracted yields alternative");
  //  legendov->AddEntry(grv2_wrt_Pt2_noZYAM,"v2_Ext2_noZYAM (pPb-like no ZYAM)");
      //      legendov->AddEntry(grv2_wrt_Pt3,"v_{2,Ext3} (pPb-like)");
  //  legendov->AddEntry(grv2_wrt_Pt3_noZYAM,"v2_Ext3_noZYAM (pPb-like no ZYAM)");
        //    legendov->AddEntry(grv2_wrt_PtZYAM,"v_{2,ZYAM} : (Y_{C}-Y_{C}(0))-(Y_{P}-Y_{P}(0))=a_{0}+2a_{1}cos(#Delta#phi)+2a_{2}cos(2#Delta#phi)");
          //   legendov->AddEntry(grv2_wrt_PtCvetanQuentin,"v_{2,Template+PeriphZYAM}: Y_{C}-F(Y_{P}-Y_{P}(0))=Y_{C,0}(1+v_{2,2}cos(2#Delta#phi))");
            // legendov->AddEntry(grv2_wrt_PtCvetanQuentinMe,"v2_CvetanQuentinMe");
             legendov->AddEntry(grv2_wrt_PtPRL,"ATLAS Template fit: Y_{C}-F(Y_{P}-Y_{P}(0))=G(1+#it{v}_{2,2}cos(2#Delta#phi))");
//    legendov->AddEntry(grv2_wrt_PtPbPb,"Pb-Pb, 5.02 TeV, 5-20%");
//    legendov->AddEntry(grv2_wrt_PtpPb,"p-Pb, 5.02-8.16 TeV, y>0 direction");
             //legendov->AddEntry(grv2_wrt_PtPRLPeriphZYAM,"v_{2,TemplateATLAS+ZYAM}: (Y_{C}-Y_{C}(0))-F(Y_{P}-Y_{P}(0))=G(1+v_{2,2}cos(2#Delta#phi))");
              legendov->Draw();
    TLegend *legendov2=new TLegend(0.12,0.60,0.40,0.70);
                  legendov2->SetFillColorAlpha(kWhite, 0.);
                  legendov2->SetBorderSize(0);
                   legendov2->SetTextFont(72);
                   legendov2->SetTextSize(0.03);
        Char_t messago[80];
        sprintf(messago,"ALICE Requested for preliminary");
       //  legendov2->AddEntry(grv2_wrt_PtPRL,messago,"");
                  legendov2->Draw();
    
    TLegend *legendov3=new TLegend(0.6,0.70,0.90,0.90);
              legendov3->SetFillColorAlpha(kWhite, 0.);
              legendov3->SetBorderSize(0);
               legendov3->SetTextFont(62);
               legendov3->SetTextSize(0.03);
    sprintf(messago,"V0M (0-5%)-(40-100%)");
    // legendov3->AddEntry(grv2_wrt_PtPRL,messago,"");
    sprintf(messago,"ALICE pp, #sqrt{s_{NN}}=13 TeV");
    legendov3->AddEntry(grv2_wrt_PtPRL,messago,"");
    sprintf(messago,"1.5<|#Delta#eta|<5.0");
    //legendov3->AddEntry(grv2_wrt_PtPRL,messago,"");
              legendov3->Draw();
    
    TLegend *legendov4=new TLegend(0.45,0.10,0.85,0.20);
              legendov4->SetFillColorAlpha(kWhite, 0.);
              legendov4->SetBorderSize(0);
               legendov4->SetTextFont(62);
               legendov4->SetTextSize(0.03);
    sprintf(messago,"3.5 percent global syst. uncertainty");
   //  legendov4->AddEntry(grv2_wrt_PtPRL,messago,"");
              legendov4->Draw();
     
     cv2Pt->Update();
    // TLine *l=new TLine(cV2Pt->GetUxmin(),0.0,cV2Pt->GetUxmax(),0.0);
     TLine *lv=new TLine(0.0,0.0,12.0,0.0);
     lv->SetLineColor(kBlack);
     lv->SetLineWidth(1);
     lv->SetLineStyle(9);
     lv->Draw();
    
    sprintf(CanvasName,"%s/Smallv2wrtPt.pdf",CanvasFolderName);
    cv2Pt->SaveAs(CanvasName);
    
    TCanvas*call = new TCanvas();
    TPad *pad1 = new TPad("pad1","pad1",0,0.5,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.5);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    pad1->SetBottomMargin(0.00001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    call->SetTitle("Preli2");
    pad1->cd();
    TVirtualPad *padvirt = c16PtBinned->GetPad(3);
    double x1 = padvirt->GetAbsXlowNDC();
    double x2 = x1 + padvirt->GetAbsWNDC();
    double y1 = padvirt->GetAbsYlowNDC();
    double y2 = y1 + padvirt->GetAbsHNDC();
    padvirt->SetPad(0,0.5,1,1);
    padvirt->DrawClone();
    padvirt->SetPad(x1,y1,x2,y2);
    pad2->cd();
    TVirtualPad *padvirt2 = cinvmassPtBinned->GetPad(3);
    double x12 = padvirt2->GetAbsXlowNDC();
    double x22 = x12 + padvirt2->GetAbsWNDC();
    double y12 = padvirt2->GetAbsYlowNDC();
    double y22 = y12 + padvirt2->GetAbsHNDC();
    padvirt2->SetPad(0,0,1,0.5);
    padvirt2->DrawClone();
    padvirt2->SetPad(x12,y12,x22,y22);
    
    sprintf(CanvasName,"%s/Preli2.pdf",CanvasFolderName);
    call->SaveAs(CanvasName);
    

//    if(doTracklets){
//
//        TCanvas* c6TKL = new TCanvas;
//            //Tracklets yield DeltaEta wrt DeltaPhi TH2 -> Projected for Periph and Central
//
//            c6TKL->Divide(2,2);
//            c6TKL->cd(1);
//            YieldTkl_allC->Draw("E");
//            c6TKL->cd(2);
//            YieldTkl_Difference->Draw("E");
//            c6TKL->cd(3);
//            YieldTkl_Central->Draw("E");
//            c6TKL->cd(4);
//            YieldTkl_Periph->Draw("E");
//            c6TKL->Draw();
//            c6TKL->Modified();
//            c6TKL->ForceUpdate();
//
//        baselineTKL_periph = (YieldTkl_Periph->GetBinContent(BinZeroLeftTKL) + YieldTkl_Periph->GetBinContent(BinZeroLeftTKL+1))/2;
//        errbaselineTKL_periph = sqrt(pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL),2) + pow(YieldTkl_Periph->GetBinError(BinZeroLeftTKL+1),2));
//
//        for(int bin_idx = 1; bin_idx<=NbinsDeltaPhiTKL; bin_idx++){
//            if(YieldTkl_Central->GetBinContent(bin_idx)<baselineTKL_central){
//                baselineTKL_central = YieldTkl_Central->GetBinContent(bin_idx);
//                errbaselineTKL_central = YieldTkl_Central->GetBinError(bin_idx);
//            }
//        }
//
//
//        for(int phi_idx = 0; phi_idx<NbinsDeltaPhiTKL; phi_idx++){
//            BaselineTkl_Periph->SetBinContent(phi_idx+1, baselineTKL_periph);
//            BaselineTkl_Periph->SetBinError(phi_idx+1, errbaselineTKL_periph);
//        }
//
//        YieldTkl_Periph_MinusBaseline->Add(YieldTkl_Periph,BaselineTkl_Periph,1,-1);
//
//        for(int phi_idx = 0; phi_idx<NbinsDeltaPhiTKL; phi_idx++){
//            BaselineTkl_Central->SetBinContent(phi_idx+1, baselineTKL_central);
//            BaselineTkl_Central->SetBinError(phi_idx+1, errbaselineTKL_central);
//        }
//
//        YieldTkl_Central_MinusBaseline->Add(YieldTkl_Central,BaselineTkl_Central,1,-1);
//
//
//        TCanvas*cTKLCentralMinusBaseline=new TCanvas();
//        cTKLCentralMinusBaseline->Divide(1,3);
//        cTKLCentralMinusBaseline->cd(1);
//        YieldTkl_Central->DrawCopy();
//        cTKLCentralMinusBaseline->cd(2);
//        BaselineTkl_Central->DrawCopy();
//        cTKLCentralMinusBaseline->cd(3);
//        YieldTkl_Central_MinusBaseline->DrawCopy();
//
//        TCanvas*cTKLPeriphMinusBaseline=new TCanvas();
//        cTKLPeriphMinusBaseline->Divide(1,3);
//        cTKLPeriphMinusBaseline->cd(1);
//        YieldTkl_Periph->DrawCopy();
//        cTKLPeriphMinusBaseline->cd(2);
//        BaselineTkl_Periph->DrawCopy();
//        cTKLPeriphMinusBaseline->cd(3);
//        YieldTkl_Periph_MinusBaseline->DrawCopy();
//
//        TCanvas*c14TKL=new TCanvas();
//        //Tracklets Yield difference wrt Phi fit
//
//
//        TH1F *YieldTkl_Central_Cvetan = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_Cvetan");
//        TH1F *YieldTkl_Central_CvetanMe = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_CvetanMe");
//        TH1F *YieldTkl_Central_ZYAM = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_ZYAM");
//        TH1F *YieldTkl_Central_PRLTemplate = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_PRLTemplate");
//        TH1F *YieldTkl_Central_PRLTemplate_PeriphZYAM = (TH1F*)YieldTkl_Central->Clone("YieldTkl_Central_PRLTemplate_PeriphZYAM");
//
//        //Methode classique C-P
//        {
//        c14TKL->cd();
//        // Ici on fit YieldsWrtDeltaPhiMassBin_DifferenceProj
//            TH1F *histo = YieldTkl_Difference;
//          // create a TF1 with the range from 0 to 3 and 6 parameters
//          TF1 *fitFcnV2TKL = new TF1("fitFcnV2TKL",FourierV5,-TMath::Pi()/2,1.5*TMath::Pi(),3);
//          fitFcnV2TKL->SetNpx(500);
//          fitFcnV2TKL->SetLineWidth(4);
//          fitFcnV2TKL->SetLineColor(kMagenta);
//          // first try without starting values for the parameters
//          // This defaults to 1 for each param.
//          // this results in an ok fit for the polynomial function
//          // however the non-linear part (lorenzian) does not
//          // respond well.
//           Double_t params[3] = {1,0.01,0.01};
//          fitFcnV2TKL->SetParameters(params);
//           TVirtualFitter::Fitter(histo)->SetMaxIterations(10000);
//           TVirtualFitter::Fitter(histo)->SetPrecision();
//        //  histo->Fit("fitFcn","0");
//          // second try: set start values for some parameters
//
//           fitFcnV2TKL->SetParName(0,"a0");
//           fitFcnV2TKL->SetParName(1,"a1");
//           fitFcnV2TKL->SetParName(2,"a2");
////            fitFcnV2TKL->SetParName(3,"a3");
////            fitFcnV2TKL->SetParName(4,"a4");
////        fitFcnV2TKL->SetParName(5,"a5");
////        fitFcnV2TKL->SetParName(6,"a6");
////        fitFcnV2TKL->SetParName(7,"a7");
////        fitFcnV2TKL->SetParName(8,"a8");
////        fitFcnV2TKL->SetParName(9,"a9");
////        fitFcnV2TKL->SetParName(10,"a10");
////        fitFcnV2TKL->SetParName(11,"a11");
////        fitFcnV2TKL->SetParLimits(6,-0.1,0.1);
//
//          TFitResultPtr res = histo->Fit("fitFcnV2TKL","SBMERI+","ep");
//          // improve the pictu
//        //   std::cout << "integral error: " << integralerror << std::endl;
//            Double_t par[3];
//            fitFcnV2TKL->GetParameters(par);
//          fitFcnV2TKL->Draw("same");
//          // draw the legend
//          TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
//          legend->SetTextFont(72);
//          legend->SetTextSize(0.04);
//            Char_t message[80];
//            sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2TKL->GetChisquare(),fitFcnV2TKL->GetNDF());
//            legend->AddEntry(fitFcnV2TKL,message);
//        sprintf(message,"V2 Tkl-Tkl: %.4f / (%.4f + %.4f) = %.4f +- %.4f",par[2],par[0], baselineTKL_periph,par[2]/(par[0] + baselineTKL_periph),(par[2]/(par[0] + baselineTKL_periph)*sqrt(pow(fitFcnV2TKL->GetParError(2)/par[2],2)+pow(fitFcnV2TKL->GetParError(0)/par[0],2)+pow(errbaselineTKL_periph/baselineTKL_periph,2))));
//        legend->AddEntry(fitFcnV2TKL,message);
//        if(res->CovMatrixStatus() == 3){
//                   sprintf(message,"The fit is a success");
//               }
//               else{
//                   sprintf(message,"The fit is a failure");
//               }
//               legend->AddEntry(fitFcnV2TKL,message);
//          legend->AddEntry(histo,"Data","lpe");
//          legend->Draw();
//        }
//
//
//
//
//        //Cvetan-Quentin fit
//
//        TCanvas* cTKLCvetan = new TCanvas;
//        cTKLCvetan->SetTitle("TKL Cvetan-Quentin Fit");
//        cTKLCvetan->Divide(1,1);
//        cTKLCvetan->cd(1);
//        YieldTkl_Central->DrawCopy();
//
//        {
//            cTKLCvetan->cd(1);
//                 TF1 *fitFcnV2_Cvetan = new TF1("fitFcnV2_Cvetan",CvetanFTKL,-TMath::Pi()/2,1.5*TMath::Pi(),4);
//                 fitFcnV2_Cvetan->SetNpx(500);
//                 fitFcnV2_Cvetan->SetLineWidth(4);
//                 fitFcnV2_Cvetan->SetLineColor(kBlue);
//                 // first try without starting values for the parameters
//                 // This defaults to 1 for each param.
//                 // this results in an ok fit for the polynomial function
//                 // however the non-linear part (lorenzian) does not
//                 // respond well.
//                  Double_t params[4] = {1,0,0.01,1};
//                 fitFcnV2_Cvetan->SetParameters(params);
//                  TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->SetMaxIterations(10000);
//                  TVirtualFitter::Fitter(YieldTkl_Central_Cvetan)->SetPrecision();
//               // gStyle->SetOptFit(1011);
//               //  histo->Fit("fitFcn","0");
//                 // second try: set start values for some parameters
//
//                  fitFcnV2_Cvetan->SetParName(0,"V0");
//                    fitFcnV2_Cvetan->SetParName(1,"V1");
//                    fitFcnV2_Cvetan->SetParName(2,"V2");
//                  fitFcnV2_Cvetan->SetParName(3,"F");
//
//            fitFcnV2_Cvetan->SetParLimits(3,0.1,50);
//
//            fitFcnV2_Cvetan->FixParameter(0,1);
//            fitFcnV2_Cvetan->FixParameter(1,0);
//          //  fitFcnV2_Cvetan->FixParameter(3,1);
//
//
//                 TFitResultPtr res = YieldTkl_Central_Cvetan->Fit("fitFcnV2_Cvetan","SBMERI+","ep");
//                Double_t par[4];
//                fitFcnV2_Cvetan->GetParameters(par);
//                 // improve the pictu
//               //   std::cout << "integral error: " << integralerror << std::endl;
//                 fitFcnV2_Cvetan->Draw("same");
//                 // draw the legend
//                 TLegend *legend=new TLegend(0.15,0.65,0.3,0.85);
//                 legend->SetTextFont(72);
//                 legend->SetTextSize(0.04);
//                   Char_t message[80];
//                   sprintf(message,"Global Fit : #chi^{2}/NDF = %.2f / %d",fitFcnV2_Cvetan->GetChisquare(),fitFcnV2_Cvetan->GetNDF());
//                   legend->AddEntry(fitFcnV2_Cvetan,message);
//            if(res->CovMatrixStatus() == 3){
//                       sprintf(message,"The fit is a success");
//                   }
//                   else{
//                       sprintf(message,"The fit is a failure");
//                   }
//                   legend->AddEntry(fitFcnV2_Cvetan,message);
//                 legend->AddEntry(YieldTkl_Central_Cvetan,"Data","lpe");
//                 legend->Draw();
//
//        }
//
//
//
//    }
    
    
    
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
    std::cout << "===== Finished: Generate .csv =====" <<std::endl;
    
    
    
    std::ofstream myfiletx(SystematicsFileName, std::ofstream::out);
    
    if(!myfiletx)
       {
           cout << "Couldn't open syst file" << endl;
       }
//    Char_t ptbinschar[500] = "";
//    for(int bin=0; bin<NbPtBins; bin++){
//        cout <<"ptbinschar avant tour "<< bin << ": "<<ptbinschar<<endl;
//        cout << "PtBins[bin] = "<<PtBins[bin]<<endl;
//        sprintf(ptbinschar,"%s%f,", ptbinschar,PtBins[bin]);
//        cout <<"ptbinschar apres tour "<< bin << ": "<<ptbinschar<<endl;
//    }
//
//    cout <<"ptbinschar"<<ptbinschar<<endl;
//
//    Char_t resultsv2[500] = "";
//    for(int bin=0; bin<NbPtBins; bin++){
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_Ext1[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_Ext1[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_Ext2[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_Ext2[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_Ext3[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_Ext3[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_Ext1_noZYAM[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_Ext1_noZYAM[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_CvetanQuentin[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_CvetanQuentin[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_ZYAM[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_ZYAM[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_PRL[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_PRL[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_PRLPeriphZYAM[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_PRLPeriphZYAM[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,v2_Ext1_V3[bin]);
//        sprintf(resultsv2,"%s%0.7f,", resultsv2,errv2_Ext1_V3[bin]);
//    }
//
//    cout <<"resultsv2"<<resultsv2<<endl;
    
    myfiletx << DPhiCutText << "," << CentEstText << "," << CentralClassText << "," << PeriphClassText << "," << int(NbPtBins) << ",";
    
    for(int bin=0; bin<NbPtBins+1; bin++){
        myfiletx << PtBins[bin] << ",";
    }
    
    myfiletx << DeltaEtaDimuCut << "," << MaxDeltaEta << "," << ZvtxCutText << "," << EMNormText << "," << EMPoolMaxText << "," << EMPoolThresholdText << "," << EMChangeText << "," << SummationText << "," << BackgroundV2F << "," << minv2 << "," << maxv2 << "," << minmass << "," << maxmass << "," << SignalF << "," << BackgroundF << "," << ratsigma << ",";
    
    for(int bin=0; bin<NbPtBins; bin++){
        myfiletx << v2_Ext1[bin] << "," << errv2_Ext1[bin] << "," << v2_Ext2[bin] << "," << errv2_Ext2[bin] << "," << v2_Ext3[bin] << "," << errv2_Ext3[bin] << "," << v2_Ext1_noZYAM[bin] << "," << errv2_Ext1_noZYAM[bin] << "," << v2_CvetanQuentin[bin] << "," << errv2_CvetanQuentin[bin] << "," << v2_ZYAM[bin] << "," << errv2_ZYAM[bin] << "," << v2_PRL[bin] << "," << errv2_PRL[bin] << "," << v2_PRLPeriphZYAM[bin] << "," << errv2_PRLPeriphZYAM[bin] << "," << v2_Ext1_V3[bin] << "," << errv2_Ext1_V3[bin];
//        if(bin<NbPtBins-1){
            myfiletx << ",";
//        }
    }
    
    for(int bin=0; bin<NbPtBins; bin++){
        myfiletx << v2_ATLASExt2[bin] << "," << errv2_ATLASExt2[bin];
        if(bin<NbPtBins-1){
            myfiletx << ",";
        }
    }
    
    myfiletx.close();
    
    sprintf(CanvasName,"%s/InvMass.pdf",CanvasFolderName);
    cinvmass->SaveAs(CanvasName);
    sprintf(CanvasName,"%s/InvMassBinnedPt.pdf",CanvasFolderName);
    cinvmassPtBinned->SaveAs(CanvasName);
    sprintf(CanvasName,"%s/InvMassBinnedPt_noZYAM.pdf",CanvasFolderName);
    cinvmassPtBinned_noZYAM->SaveAs(CanvasName);
}
    












// FITTING INVARIANT MASS METHODS

Double_t CvetanF(Double_t *x,Double_t *par)

{   //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double YMinusBp = Yields_Periph_1_MinusBaseline->GetBinContent(bintolook+1);
    return baseline_central*(par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + par[3]*YMinusBp; }

void ChisquareCvetanF(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-CvetanF(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(par[3]*Yields_Periph_1_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t CvetanFPtBinned(Double_t *x,Double_t *par)

{   //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    int ptbin = par[0];
    double YMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinContent(bintolook+1);
    return baseline_centralPtBinned[ptbin]*(par[1] + 2*par[2]*cos(x[0]) + 2*par[3]*cos(2*x[0])) + par[4]*YMinusBp; }

void ChisquareCvetanFPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-CvetanFPtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(par[4]*Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t CvetanFTKL(Double_t *x,Double_t *par)

{   //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double YMinusBp = YieldTkl_Periph_MinusBaseline->GetBinContent(bintolook+1);
    return baselineTKL_central*(par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + par[3]*YMinusBp; }

Double_t ZYAM(Double_t *x,Double_t *par)

{   //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double YMinusBp = Yields_Periph_1_MinusBaseline->GetBinContent(bintolook+1);
    return baseline_central + (par[0] + 2*par[1]*cos(x[0]) + 2*par[2]*cos(2*x[0])) + 1*YMinusBp; }

void ChisquareZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-ZYAM(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(errbaseline_central,2)+pow(Yields_Periph_1_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t ZYAMPtBinned(Double_t *x,Double_t *par)

{   //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    int ptbin = par[0];
    double YMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinContent(bintolook+1);
    return baseline_centralPtBinned[ptbin] + (par[1] + 2*par[2]*cos(x[0]) + 2*par[3]*cos(2*x[0])) + 1*YMinusBp; }

void ChisquareZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-ZYAMPtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(errbaseline_centralPtBinned[ptbin],2)+pow(Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplate(Double_t *x,Double_t *par)

{   double integral_Yperiph = Yields_Periph_1->Integral(1,Yields_Periph_1->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double Yperiphbin = Yields_Periph_1->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(MaxDeltaPhi-MinDeltaPhi);//factor 2 deno
    
    return par[1]*Yperiphbin + (1+2*par[0]*cos(2*x[0]))*G; }

void ChisquarePRLTemplate(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-PRLTemplate(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(par[1]*Yields_Periph_1->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplatePtBinnedAlt(Double_t *x,Double_t *par)

{   int ptbin = int(par[0]);
    int massbin = int(par[3]);
    double integral_Yperiph = Yield_Periph_MassBinPtBinned[massbin][ptbin]->Integral(1,Yield_Periph_MassBinPtBinned[massbin][ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yield_Central_MassBinPtBinned[massbin][ptbin]->Integral(1,Yield_Central_MassBinPtBinned[massbin][ptbin]->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double Yperiphbin = Yield_Periph_MassBinPtBinned[massbin][ptbin]->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[2]*integral_Yperiph))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[2]*Yperiphbin + (1+2*par[1]*cos(2*x[0]))*G; }

void ChisquarePRLTemplatePtBinnedAlt(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = int(par[0]);
    int massbin = int(par[3]);
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yield_Central_MassBinPtBinned[massbin][ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yield_Central_MassBinPtBinned[massbin][ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yield_Central_MassBinPtBinned[massbin][ptbin]->GetBinContent(bin_idx+1)-PRLTemplatePtBinnedAlt(x,par))/(sqrt(pow(Yield_Central_MassBinPtBinned[massbin][ptbin]->GetBinError(bin_idx+1),2)+pow(par[2]*Yield_Periph_MassBinPtBinned[massbin][ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplatePtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_Yperiph = Yields_Periph_1PtBinned[ptbin]->Integral(1,Yields_Periph_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double Yperiphbin = Yields_Periph_1PtBinned[ptbin]->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[2]*integral_Yperiph))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[2]*Yperiphbin + (1+2*par[1]*cos(2*x[0]))*G; }

void ChisquarePRLTemplatePtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-PRLTemplatePtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(par[2]*Yields_Periph_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}

Double_t PRLTemplate_RidgeAndZero(Double_t *x,Double_t *par)

{   double integral_Yperiph = Yields_Periph_1->Integral(1,Yields_Periph_1->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double Yperiphbin = baseline_periph;
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[1]*Yperiphbin + (1+2*par[0]*cos(2*x[0]))*G; }

Double_t PRLTemplate_RidgeAndZeroPtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_Yperiph = Yields_Periph_1PtBinned[ptbin]->Integral(1,Yields_Periph_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double Yperiphbin = baseline_periphPtBinned[ptbin];
    double G = (integral_Yreal-(par[2]*integral_Yperiph))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[2]*Yperiphbin + (1+2*par[1]*cos(2*x[0]))*G; }

Double_t PRLTemplate_PeriphAndG(Double_t *x,Double_t *par)

{   double integral_Yperiph = Yields_Periph_1->Integral(1,Yields_Periph_1->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double Yperiphbin = Yields_Periph_1->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_Yperiph))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[1]*Yperiphbin + G; }

Double_t PRLTemplate_PeriphAndGPtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_Yperiph = Yields_Periph_1PtBinned[ptbin]->Integral(1,Yields_Periph_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double Yperiphbin = Yields_Periph_1PtBinned[ptbin]->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[2]*integral_Yperiph))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[2]*Yperiphbin + G; }

Double_t PRLTemplate_PeriphZYAM(Double_t *x,Double_t *par)

{   double integral_YperiphMinusBp = Yields_Periph_1_MinusBaseline->Integral(1,Yields_Periph_1_MinusBaseline->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1->Integral(1,Yields_Central_1->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double YMinusBp = Yields_Periph_1_MinusBaseline->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[1]*integral_YperiphMinusBp))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[1]*YMinusBp + (1+2*par[0]*cos(2*x[0]))*G; }

void ChisquarePRLTemplate_PeriphZYAM(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1->GetBinContent(bin_idx+1)-PRLTemplate_PeriphZYAM(x,par))/(sqrt(pow(Yields_Central_1->GetBinError(bin_idx+1),2)+pow(par[1]*Yields_Periph_1_MinusBaseline->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}


Double_t PRLTemplate_PeriphZYAMPtBinned(Double_t *x,Double_t *par)

{   int ptbin = par[0];
    double integral_YperiphMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->Integral(1,Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetNbinsX()+1,"width");
    double integral_Yreal = Yields_Central_1PtBinned[ptbin]->Integral(1,Yields_Central_1PtBinned[ptbin]->GetNbinsX()+1,"width");
    //int bintolook = floor((   (  (x[0]+(TMath::Pi()/2))  /   (2*TMath::Pi())  )    *12));
    int bintolook = floor((   (  (x[0]-MinDeltaPhi)  /   (MaxDeltaPhi-MinDeltaPhi)  )    *NbinsDeltaPhi));
    double YMinusBp = Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinContent(bintolook+1);
    double G = (integral_Yreal-(par[2]*integral_YperiphMinusBp))/(MaxDeltaPhi-MinDeltaPhi);//deno x2
    
    return par[2]*YMinusBp + (1+2*par[1]*cos(2*x[0]))*G; }

void ChisquarePRLTemplate_PeriphZYAMPtBinned(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *par, Int_t /*iflag */  )

{   int ptbin = par[0];
    npfits = 0;
    double tampon = 0;
    double chi = 0;
    double x[1];
    
    TAxis *xaxis1  = Yields_Central_1PtBinned[ptbin]->GetXaxis();
    
    for(int bin_idx=0; bin_idx<Yields_Central_1PtBinned[ptbin]->GetNbinsX(); bin_idx++){
        x[0] = xaxis1->GetBinCenter(bin_idx+1);
        tampon = (Yields_Central_1PtBinned[ptbin]->GetBinContent(bin_idx+1)-PRLTemplate_PeriphZYAMPtBinned(x,par))/(sqrt(pow(Yields_Central_1PtBinned[ptbin]->GetBinError(bin_idx+1),2)+pow(par[2]*Yields_Periph_1_MinusBaselinePtBinned[ptbin]->GetBinError(bin_idx+1),2)));
        chi += tampon*tampon;
        npfits++;
    }
    fval = chi;}


Double_t FourierV2_WrtInvMass(Double_t *x,Double_t *par)
// Par 0->7: signal, 8->11 Bkg, 12: v2 JPsi, 13->15: V2 bkg
//{ return (JPsiCrystalBallExtended(x,par)*par[12] + (ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[14]*(x[0]-mJpsi) + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}
{   Double_t xx = x[0];
    Double_t iBin = V2JPsiTkl->FindBin(xx);
    Double_t binLow = V2JPsiTkl->GetXaxis()->GetBinLowEdge(iBin);
    Double_t binHigh = V2JPsiTkl->GetXaxis()->GetBinUpEdge(iBin);
    
    TF1 *jpsi = new TF1("jpsi", UsedJPsiSignal, MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
    jpsi->SetParameters(par);
    TF1 *psip = new TF1("psip", UsedPsi2SSignal, MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
    psip->SetParameters(par);
    TF1 *bkg = new TF1("bkg", UsedBackground, MinV2Fit,MaxV2Fit,numParameters+numParametersBkgV2+1);
    bkg->SetParameters(par);
    
    double jpsinum = jpsi->Integral(binLow,binHigh);
    double psipnum = psip->Integral(binLow,binHigh);
    double bkgnum = bkg->Integral(binLow,binHigh);
    
    delete jpsi;
    delete psip;
    delete bkg;
    
    if(GlobalBackgroundV2=="Pol1"){
      return (jpsinum*par[numParameters+numParametersBkgV2-2] + (bkgnum+psipnum)*(par[numParameters+numParametersBkgV2-1]*(x[0]-mJpsi) + par[numParameters+numParametersBkgV2]))/(jpsinum+psipnum+bkgnum);
      }
    if(GlobalBackgroundV2=="Pol2"){
    return (jpsinum*par[numParameters+numParametersBkgV2-3] + (bkgnum+psipnum)*(par[numParameters+numParametersBkgV2-2]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[numParameters+numParametersBkgV2-1]*(x[0]-mJpsi) + par[numParameters+numParametersBkgV2]))/(jpsinum+psipnum+bkgnum);
    }
      if(GlobalBackgroundV2=="Pol3"){
      return (jpsinum*par[numParameters+numParametersBkgV2-4] + (bkgnum+psipnum)*(par[numParameters+numParametersBkgV2-3]*(x[0]-mJpsi)*(x[0]-mJpsi)*(x[0]-mJpsi) + par[numParameters+numParametersBkgV2-2]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[numParameters+numParametersBkgV2-1]*(x[0]-mJpsi) + par[numParameters+numParametersBkgV2]))/(jpsinum+psipnum+bkgnum);
      }
      else{
          return 0;
      }
}
    //return (jpsi->Integral(binLow,binHigh)*par[12] + (bkg->Integral(binLow,binHigh))*(par[13]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[14]*(x[0]-mJpsi) + par[15]))/(jpsi->Integral(binLow,binHigh)+bkg->Integral(binLow,binHigh)) ;}

//Double_t BackFcnV2(Double_t *x,Double_t *par)
//{return ((ExpBkg(x,&par[8])+Psi2SCrystalBallExtended(x,par))*(par[13]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[14]*(x[0]-mJpsi) + par[15]))/(JPsiCrystalBallExtended(x,par)+Psi2SCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}
////{return ((ExpBkg(x,&par[8]))*(par[13]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[14]*(x[0]-mJpsi) + par[15]))/(JPsiCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t BackFcnV2Poly(Double_t *x,Double_t *par)
{ if(GlobalBackgroundV2=="Pol1"){
    return (par[numParameters+1]*(x[0]-mJpsi) + par[numParameters+2]);
    }
  if(GlobalBackgroundV2=="Pol2"){
  return (par[numParameters+1]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[numParameters+2]*(x[0]-mJpsi) + par[numParameters+3]);
  }
    if(GlobalBackgroundV2=="Pol3"){
    return (par[numParameters+1]*(x[0]-mJpsi)*(x[0]-mJpsi)*(x[0]-mJpsi) + par[numParameters+2]*(x[0]-mJpsi)*(x[0]-mJpsi) + par[numParameters+3]*(x[0]-mJpsi) + par[numParameters+4]);
    }
    else{
        return 0;
    }
}

Double_t SignalFcnJPsiV2(Double_t *x,Double_t *par)
{return (UsedJPsiSignal(x,par)*par[numParameters])/(UsedJPsiSignal(x,par)+UsedPsi2SSignal(x,par)+UsedBackground(x,par)) ;}
//{return (JPsiCrystalBallExtended(x,par)*par[12])/(JPsiCrystalBallExtended(x,par)+ExpBkg(x,&par[8])) ;}

Double_t UsedJPsiSignal(Double_t *x,Double_t *par)
{
    if(GlobalSignal=="CB2-Run2" || GlobalSignal=="CB2-MC" || GlobalSignal=="CB2-FREE"){
        return JPsiCrystalBallExtended(x,par);
    }
    else if(GlobalSignal=="NA60-MC" || GlobalSignal=="NA60-FREE"){
        return JPsiNA60(x,par);
    }
    else{
        return 0;
    }
}

Double_t UsedPsi2SSignal(Double_t *x,Double_t *par)
{
    if(GlobalSignal=="CB2-Run2" || GlobalSignal=="CB2-MC" || GlobalSignal=="CB2-FREE"){
        return Psi2SCrystalBallExtended(x,par);
    }
    else if(GlobalSignal=="NA60-MC" || GlobalSignal=="NA60-FREE"){
        return Psi2SNA60(x,par);
    }
    else{
        return 0;
    }
}

Double_t UsedBackground(Double_t *x,Double_t *par)
{
    if(GlobalBackground =="DoubleExpo"){
        return ExpBkg(x,&par[numParameters-numParametersBkg]);
    }
    else if(GlobalBackground =="POL1POL2"){
        return Pol1Pol2(x,&par[numParameters-numParametersBkg]);
    }
    else if(GlobalBackground =="VWG"){
        return VWGaussian(x,&par[numParameters-numParametersBkg]);
    }
      else if(GlobalBackground =="Tchebychev"){
          return Tchebychev(x,&par[numParameters-numParametersBkg]);
      }
      else{
          return 0;
      }
}

Double_t UsedBackgroundHole(Double_t *x,Double_t *par)
{
    if(GlobalBackground =="DoubleExpo"){
        return ExpBkgHole(x,&par[numParameters-numParametersBkg]);
    }
    else if(GlobalBackground =="POL1POL2"){
        return Pol1Pol2Hole(x,&par[numParameters-numParametersBkg]);
    }
    else if(GlobalBackground =="VWG"){
        return VWGaussianHole(x,&par[numParameters-numParametersBkg]);
    }
      else if(GlobalBackground =="Tchebychev"){
          return TchebychevHole(x,&par[numParameters-numParametersBkg]);
      }
      else{
          return 0;
      }
}


Double_t FourierV2(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);}

Double_t FourierV3(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]) + 2*par[3]*cos(3*x[0]);}

Double_t FourierV5(Double_t *x,Double_t *par)

{ return par[0] + 2*par[2]*cos(2*x[0]) + 2*par[1]*cos(x[0]);} //+ 2*par[3]*cos(3*x[0]) + 2*par[4]*cos(4*x[0]);} + 2*par[5]*cos(5*x[0]) + 2*par[6]*cos(6*x[0]) + 2*par[7]*cos(7*x[0]) + 2*par[8]*cos(8*x[0]) + 2*par[9]*cos(9*x[0]) + 2*par[10]*cos(10*x[0]) + 2*par[11]*cos(11*x[0]);}

Double_t MassFunction(Double_t *x,Double_t *par)

{ return UsedJPsiSignal(x,par) + UsedPsi2SSignal(x,par) + UsedBackground(x,par);}
//{ return JPsiCrystalBallExtended(x,par) + ExpBkg(x,&par[8]);}

Double_t ExpBkg(Double_t *x,Double_t *par)

{ return exp(x[0]*par[1]*(-1)) + par[2]*(exp(x[0]*par[3]*(-1))); }

Double_t Pol1Pol2(Double_t *x,Double_t *par)

{ return par[0]*(((1+par[1]*x[0]))/(1+par[2]*x[0] + par[3]*(x[0]*x[0]))); }

Double_t VWGaussian(Double_t *x,Double_t *par)

{ double sigma = par[2] + par[3]*((x[0]-par[1])/par[1]);
    return par[0]*exp(-1.0*pow(x[0]-par[1],2)/(2*pow(sigma,2))); }

Double_t Tchebychev(Double_t *x,Double_t *par)

{   Double_t xx = (x[0]-par[5])/(par[6]-par[5]);
    
    Double_t poly = SmolTcheby(xx,par, 0);;
    for(int index=1; index<5;index++){
        poly += par[index]*SmolTcheby(xx,par, index);
    }
    poly *= par[0];
    return poly;
}

Double_t SmolTcheby(Double_t x,Double_t *par, int n)

{ if(n==0){
    return 1;
    }
    if(n==1){
        return x;
    }
    if(n>=2){
        return 2*x*SmolTcheby(x, par, n-1) - SmolTcheby(x, par, n-2);
    }
    else{
        return 0;
    }
}


Double_t JPsiCrystalBallExtended(Double_t *x,Double_t *par)
{
    
    Double_t sum = 0;
    
    Double_t t = (x[0]-par[1])/par[2];
    if (par[3] < 0) t = -t;
    
    Double_t absAlpha = fabs((Double_t)par[3]);
    Double_t absAlpha2 = fabs((Double_t)par[5]);
    
    
    if (t <= -absAlpha) //left tail
    {
        Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
        Double_t b = par[4]/absAlpha - absAlpha;
        
        sum += (par[0])*(a/TMath::Power(b - t, par[4]));
    } else if (t > -absAlpha && t < absAlpha2) // gaussian core
    {
        sum += (par[0])*(exp(-0.5*t*t));
    } else if (t >= absAlpha2) //right
    {
        Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
        Double_t d = par[6]/absAlpha2 - absAlpha2;
        
        sum += (par[0])*(c/TMath::Power(d + t, par[6]));
    } else
        sum += 0;
    
    return sum ;
}

Double_t Psi2SCrystalBallExtended(Double_t *x,Double_t *par)
{
      Double_t sum = 0;
    
    Double_t t = (x[0]-(par[1]+mDiff))/(par[2]*par[8]);
  if (par[3] < 0) t = -t;
  
  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  
  if (t < -absAlpha) //left tail
  {
      Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
      Double_t b = par[4]/absAlpha - absAlpha;
      
      sum += (par[7])*(a/TMath::Power(b - t, par[4]));
  }
  else if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
      sum += (par[7])*(exp(-0.5*t*t));
  }
  else if (t >= absAlpha2) //right tail
  {
      Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
      Double_t d = par[6]/absAlpha2 - absAlpha2;
      
      sum += (par[7])*(c/TMath::Power(d + t, par[6]));
  } else
      sum += 0;
    
    return sum ;

}


Double_t JPsiNA60(Double_t *x,Double_t *par)
{
    
    Double_t sum = 0;
    
    Double_t t0 = 1;
    
    Double_t t = (x[0]-par[1])/par[2];
    
    if (t < par[3]) //left tail
    {
        t0 = 1+pow(par[4]*(par[3]-t),par[5]-(par[6]*sqrt(par[3]-t)));

        sum += par[0]*exp(-0.5*pow(t/t0,2));
    }
    else if (t > par[7]) //right
    {
        t0 = 1+pow(par[8]*(t-par[7]),par[9]-(par[10]*sqrt(t-par[7])));

        sum += par[0]*exp(-0.5*pow(t/t0,2));
    }
    else{
        t0 = 1;
        
        sum += par[0]*exp(-0.5*pow(t/t0,2));
    }
    
    return sum ;
}

Double_t Psi2SNA60(Double_t *x,Double_t *par)
{
    
    Double_t sum = 0;
    
    Double_t t0 = 1;
    
    Double_t t = (x[0]-(par[1]+mDiff))/(par[2]*par[12]);
    
    if (t < par[3]) //left tail
    {
        t0 = 1+pow(par[4]*(par[3]-t),par[5]-(par[6]*sqrt(par[3]-t)));
        
        sum += par[11]*exp(-0.5*pow(t/t0,2));
    }
    else if (t > par[7]) //right
    {
        t0 = 1+pow(par[8]*(t-par[7]),par[9]-(par[10]*sqrt(t-par[7])));
        
        sum += par[11]*exp(-0.5*pow(t/t0,2));
    }
    else{
        t0 = 1;
        
        sum += par[11]*exp(-0.5*pow(t/t0,2));
    }
    
    return sum ;
}

    
//    TFitResultPtr FittingAllInvMass(const char *histoname, TCanvas *cinvmass){
//        TFitResultPtr res = FittingAllInvMassBin(histoname, cinvmass, 0);
//        return res;
//    }

TFitResultPtr FittingAllInvMassBin(const char *histoname, TCanvas *cinvmass, int i, string SignalF, string BackgroundF, double min, double max, double ratsigma){
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
    TFile *file0 = new TFile(FitFileName); //POUET
    
    //ofstream myfiletxt;
    //myfiletxt.open("/tmp/debugging.txt");
    
    int nextPar = 0;
    
    cinvmass->cd(i+1);
   cinvmass->SetFillColor(33);
   cinvmass->SetFrameFillColor(41);
   cinvmass->SetGrid();
   TH1F *histo = (TH1F*)file0->Get(histoname);
  //  THnSparse *histoa = (THnSparse*)file0->Get("CMUL_fhDimuon");
    //  TH1F* histo = (TH1F*)histoa->Projection(0);
  //  create a TF1 with the range from 0 to 3 and 6 parameters
    
    TF1 *fitFcn = NULL;
    TF1 *backFcn1 = NULL;
     TFitResultPtr res2;
    bool isFitRetried = kFALSE;
    
    fitFcn = new TF1("fitFcn",MassFunction,min,max,numParameters);
    backFcn1 = new TF1("backFcn1",UsedBackgroundHole,1.,5.,numParameters);
    
   fitFcn->SetNpx(500);
   fitFcn->SetLineWidth(2);
   fitFcn->SetLineColor(kAzure+3);
   // first try without starting values for the parameters
   // This defaults to 1 for each param.
   // this results in an ok fit for the polynomial function
   // however the non-linear part (lorenzian) does not
   // respond well.
    Double_t params[numParameters];
    cout << "numParameters "<<numParameters<<endl;
    std::fill_n(params, numParameters, 1);
    params[0] = 10;
    cout << "params[3] "<<params[3]<<endl;
    cout << "SignalF "<<SignalF<<endl;
    cout << "BackgroundF "<<BackgroundF<<endl;
   fitFcn->SetParameters(params);
    TVirtualFitter::Fitter(histo)->SetMaxIterations(10000000);
    TVirtualFitter::Fitter(histo)->SetPrecision();
 //  histo->Fit("fitFcn","0");
   // First fit, fix MJPsi and Sigma
    
     if(SignalF=="CB2-Run2" || SignalF=="CB2-MC" || SignalF=="CB2-FREE"){
         nextPar = 9;
         
         fitFcn->SetParName(0,"Norm_{J/#psi}");
         fitFcn->SetParName(1,"M_{J/#psi}");
         fitFcn->SetParName(2,"Sigma_{J/#psi}");
         fitFcn->SetParName(3,"a_{1}");
         fitFcn->SetParName(4,"n_{1}");
         fitFcn->SetParName(5,"a_{2}");
         fitFcn->SetParName(6,"n_{2}");
         fitFcn->SetParName(7,"Norm_{#Psi(2S)}");
         fitFcn->SetParName(8,"Sigma_{ratio}");
         
         
           fitFcn->SetParLimits(0,0.0001,100000000);
           fitFcn->SetParLimits(1,3.05,3.15);
           fitFcn->SetParLimits(2,0.06,0.1);
           fitFcn->SetParLimits(3,0.7,1.1);
           fitFcn->SetParLimits(4,1,50);
           fitFcn->SetParLimits(5,1,10);
           fitFcn->SetParLimits(6,5,50);
           fitFcn->SetParLimits(7,0.001,1000000);
         fitFcn->SetParLimits(8,0.1,2);
         
         fitFcn->SetParameter(7,100);
         fitFcn->SetParameter(0,3300);
         
         fitFcn->FixParameter(1,mJpsi); // Mean x core //mJpsi
         fitFcn->FixParameter(2,0.0708);
         fitFcn->FixParameter(8,ratsigma);
         
         if(SignalF=="CB2-Run2"){
              fitFcn->FixParameter(3,0.883);
              fitFcn->FixParameter(4,9.940);
              fitFcn->FixParameter(5,1.832);
              fitFcn->FixParameter(6,15.323);
         }
         else if(SignalF=="CB2-MC"){
              fitFcn->FixParameter(3,0.993);
              fitFcn->FixParameter(4,2.9075);
              fitFcn->FixParameter(5,2.182);
              fitFcn->FixParameter(6,3.122);
         }
         
         else if(SignalF=="CB2-FREE"){
              fitFcn->FixParameter(3,0.8842);
              fitFcn->FixParameter(4,14.54);
              fitFcn->FixParameter(5,1.855);
              fitFcn->FixParameter(6,21.46);
         }
         
       }
       if(SignalF=="NA60-MC" || SignalF=="NA60-FREE"){
           
           nextPar = 13;
           
           fitFcn->SetParName(0,"Norm_{J/#psi}");
           fitFcn->SetParName(1,"M_{J/#psi}");
           fitFcn->SetParName(2,"Sigma_{J/#psi}");
           fitFcn->SetParName(3,"a_{L}");
           fitFcn->SetParName(4,"p1_{L}");
           fitFcn->SetParName(5,"p2_{L}");
           fitFcn->SetParName(6,"p3_{L}");
           fitFcn->SetParName(7,"a_{R}");
           fitFcn->SetParName(8,"p1_{R}");
           fitFcn->SetParName(9,"p2_{R}");
           fitFcn->SetParName(10,"p3_{r}");
           fitFcn->SetParName(11,"Norm_{#Psi(2S)}");
           fitFcn->SetParName(12,"Sigma_{ratio}");
           
           fitFcn->SetParLimits(0,0.0001,10000000);
           fitFcn->SetParLimits(1,3.05,3.15);
           fitFcn->SetParLimits(2,0.06,0.1);
           fitFcn->SetParLimits(3,-1.,1.);
           fitFcn->SetParLimits(4,0,10);
           fitFcn->SetParLimits(5,0,10);
           fitFcn->SetParLimits(6,0,10);
           fitFcn->SetParLimits(7,0,10);
           fitFcn->SetParLimits(8,0,10);
           fitFcn->SetParLimits(9,0,10);
           fitFcn->SetParLimits(10,0,10);
           fitFcn->SetParLimits(11,0.001,1000000);
           fitFcn->SetParLimits(12,0.1,2);
           
           fitFcn->SetParameter(0,10000);
           fitFcn->SetParameter(11,100);
           
           fitFcn->FixParameter(1,mJpsi); // Mean x core
           fitFcn->FixParameter(2,0.07);
           fitFcn->FixParameter(12,ratsigma);
           
           if(SignalF=="NA60-MC"){
               fitFcn->FixParameter(3,-0.4061);
               fitFcn->FixParameter(4,0.2302);
               fitFcn->FixParameter(5,1.2048);
               fitFcn->FixParameter(6,0.0390);
               fitFcn->FixParameter(7,2.0627);
               fitFcn->FixParameter(8,0.1836);
               fitFcn->FixParameter(9,1.2989);
               fitFcn->FixParameter(10,0.0643);
           }
           
       }
    
    if(BackgroundF =="DoubleExpo"){
        fitFcn->SetParName(nextPar+0,"Norm_{TailLowM}");
        fitFcn->SetParName(nextPar+1,"Exp_{TailLowM}");
        fitFcn->SetParName(nextPar+2,"Norm_{TailHighM}");
        fitFcn->SetParName(nextPar+3,"Exp_{TailHighM}");
        
        fitFcn->SetParLimits(nextPar+0,0.1,10000000);
        fitFcn->SetParLimits(nextPar+1,-20,20);
        fitFcn->SetParLimits(nextPar+2,0.01,100000000);
        fitFcn->SetParLimits(nextPar+3,0.01,50);
        
        fitFcn->FixParameter(nextPar+0,1);
        fitFcn->SetParameter(nextPar+1,0.01);
        fitFcn->SetParameter(nextPar+2,100);
        fitFcn->SetParameter(nextPar+3,1.4);
    }
       if(BackgroundF =="POL1POL2"){
           
           fitFcn->SetParName(nextPar+0,"Norm_{Bkg}");
           fitFcn->SetParName(nextPar+1,"a_{1}");
           fitFcn->SetParName(nextPar+2,"b_{1}");
           fitFcn->SetParName(nextPar+3,"b_{2}");
           
           fitFcn->SetParLimits(nextPar+0,10,40000);
           fitFcn->SetParLimits(nextPar+1,-0.5,0.5);
           fitFcn->SetParLimits(nextPar+2,-2,0);
           fitFcn->SetParLimits(nextPar+3,0,1);

           
           fitFcn->SetParameter(nextPar+0,30000);
           fitFcn->SetParameter(nextPar+1,-0.16);
           fitFcn->SetParameter(nextPar+2,-1.);
           fitFcn->SetParameter(nextPar+3,0.5);
       }
    
    if(BackgroundF =="VWG"){
        fitFcn->SetParName(nextPar+0,"Norm_{Bkg}");
        fitFcn->SetParName(nextPar+1,"Mean_{Bkg}");
        fitFcn->SetParName(nextPar+2,"Alpha");
        fitFcn->SetParName(nextPar+3,"Beta");
        
        fitFcn->SetParLimits(nextPar+0,0.001,10000000);
        fitFcn->SetParLimits(nextPar+1,0.001,3);
        fitFcn->SetParLimits(nextPar+2,0.001,10);
        fitFcn->SetParLimits(nextPar+3,0.001,10);
        
        fitFcn->SetParameter(nextPar+0,30000);
        fitFcn->SetParameter(nextPar+1,1.5);
        fitFcn->SetParameter(nextPar+2,0.5);
        fitFcn->SetParameter(nextPar+3,0.5);
    }
    
    if(BackgroundF =="Tchebychev"){
        fitFcn->SetParName(nextPar+0,"N");
        fitFcn->SetParName(nextPar+1,"c_{1}");
        fitFcn->SetParName(nextPar+2,"c_{2}");
        fitFcn->SetParName(nextPar+3,"c_{3}");
        fitFcn->SetParName(nextPar+4,"c_{4}");
        fitFcn->SetParName(nextPar+5,"min");
        fitFcn->SetParName(nextPar+6,"max");
        
        fitFcn->SetParLimits(nextPar+0,0,1000000);
        fitFcn->SetParLimits(nextPar+1,-10,10);
        fitFcn->SetParLimits(nextPar+2,-1,1);
        fitFcn->SetParLimits(nextPar+3,-1,1);
        fitFcn->SetParLimits(nextPar+4,-1,1);
        
        fitFcn->SetParameter(nextPar+0,195000);
        fitFcn->SetParameter(nextPar+1,-0.8);
        fitFcn->SetParameter(nextPar+2,0.15);
        fitFcn->SetParameter(nextPar+3,-0.009);
        fitFcn->SetParameter(nextPar+4,0.0002);
        fitFcn->FixParameter(nextPar+5,min);
        fitFcn->FixParameter(nextPar+6,max);
    }
    
    Double_t par[numParameters];
    
   TFitResultPtr res = histo->Fit("fitFcn","SBMER","ep");
   //res = histo->Fit("fitFcn","SBMER","ep");
     fitFcn->ReleaseParameter(1);
   res = histo->Fit("fitFcn","SBMER","ep");
     fitFcn->ReleaseParameter(2);
   res = histo->Fit("fitFcn","SBMER","ep");
    if(SignalF=="CB2-Run2" || SignalF=="CB2-MC" || SignalF=="CB2-FREE"){
        fitFcn->ReleaseParameter(7);
        res = histo->Fit("fitFcn","SBMER","ep");
    }
    if(SignalF=="NA60-MC" || SignalF=="NA60-FREE"){
        fitFcn->ReleaseParameter(11);
        res = histo->Fit("fitFcn","SBMER","ep");
    }
  //  if(res->CovMatrixStatus() !=3){
       res = histo->Fit("fitFcn","SBMER","ep");
       res = histo->Fit("fitFcn","SBMER","ep");
  //  }
    
//    for(int idx = 0; idx<NumberOfParameters(SignalF,BackgroundF);idx++){
//        if(SignalF=="NA60-MC"){
//            if((idx>=3 && idx<=10)|| idx==12){ // if((idx>=3 && idx<=10)|| idx==12){
//                continue;
//            }
//        }
//        if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
//            if((idx>=3 && idx<=6)||idx==8){ //if((idx>=3 && idx<=6)||idx==8){
//                continue;
//            }
//        }
//        fitFcn->ReleaseParameter(idx);
//    }
 //   res = histo->Fit("fitFcn","SBMER","ep");
    
    //If fit failed first time, try to set bkg first
    
   // if(kFALSE){
    if(res->CovMatrixStatus() !=3){
        
        cout << "\n\n The fit is retried \n\n"<<endl;
        
        isFitRetried = kTRUE;
        
        if(SignalF=="CB2-Run2" || SignalF=="CB2-MC" || SignalF=="CB2-FREE"){
             nextPar = 9;
             
             fitFcn->SetParName(0,"Norm_{J/#psi}");
             fitFcn->SetParName(1,"M_{J/#psi}");
             fitFcn->SetParName(2,"Sigma_{J/#psi}");
             fitFcn->SetParName(3,"a_{1}");
             fitFcn->SetParName(4,"n_{1}");
             fitFcn->SetParName(5,"a_{2}");
             fitFcn->SetParName(6,"n_{2}");
             fitFcn->SetParName(7,"Norm_{#Psi(2S)}");
             fitFcn->SetParName(8,"Sigma_{ratio}");
             
             
               fitFcn->SetParLimits(0,0.00001,100000000);
               fitFcn->SetParLimits(1,3.05,3.15);
               fitFcn->SetParLimits(2,0.06,0.1);
               fitFcn->SetParLimits(3,0.7,1.1);
               fitFcn->SetParLimits(4,1,50);
               fitFcn->SetParLimits(5,1,10);
               fitFcn->SetParLimits(6,5,50);
               fitFcn->SetParLimits(7,0.0001,1000000);
             fitFcn->SetParLimits(8,0.1,2);
             
             fitFcn->SetParameter(7,100);
             fitFcn->SetParameter(0,3300);
             
             fitFcn->FixParameter(1,mJpsi); // Mean x core
             fitFcn->FixParameter(2,0.07);
             fitFcn->FixParameter(8,ratsigma);
             
             if(SignalF=="CB2-Run2"){
                  fitFcn->FixParameter(3,0.883);
                  fitFcn->FixParameter(4,9.940);
                  fitFcn->FixParameter(5,1.832);
                  fitFcn->FixParameter(6,15.323);
             }
             else if(SignalF=="CB2-MC"){
                  fitFcn->FixParameter(3,0.993);
                  fitFcn->FixParameter(4,2.9075);
                  fitFcn->FixParameter(5,2.182);
                  fitFcn->FixParameter(6,3.122);
             }
            else if(SignalF=="CB2-FREE"){
                 fitFcn->FixParameter(3,0.8842);
                 fitFcn->FixParameter(4,14.54);
                 fitFcn->FixParameter(5,1.855);
                 fitFcn->FixParameter(6,21.46);
            }
             
           }
           if(SignalF=="NA60-MC" || SignalF=="NA60-FREE"){
               
               nextPar = 13;
               
               fitFcn->SetParName(0,"Norm_{J/#psi}");
               fitFcn->SetParName(1,"M_{J/#psi}");
               fitFcn->SetParName(2,"Sigma_{J/#psi}");
               fitFcn->SetParName(3,"a_{L}");
               fitFcn->SetParName(4,"p1_{L}");
               fitFcn->SetParName(5,"p2_{L}");
               fitFcn->SetParName(6,"p3_{L}");
               fitFcn->SetParName(7,"a_{R}");
               fitFcn->SetParName(8,"p1_{R}");
               fitFcn->SetParName(9,"p2_{R}");
               fitFcn->SetParName(10,"p3_{r}");
               fitFcn->SetParName(11,"Norm_{#Psi(2S)}");
               fitFcn->SetParName(12,"Sigma_{ratio}");
               
               fitFcn->SetParLimits(0,0.0001,10000000);
               fitFcn->SetParLimits(1,3.05,3.15);
               fitFcn->SetParLimits(2,0.06,0.1);
               fitFcn->SetParLimits(3,-1.,1.);
               fitFcn->SetParLimits(4,0,10);
               fitFcn->SetParLimits(5,0,10);
               fitFcn->SetParLimits(6,0,10);
               fitFcn->SetParLimits(7,0,10);
               fitFcn->SetParLimits(8,0,10);
               fitFcn->SetParLimits(9,0,10);
               fitFcn->SetParLimits(10,0,10);
               fitFcn->SetParLimits(11,0.001,1000000);
               fitFcn->SetParLimits(12,0.1,2);
               
               fitFcn->SetParameter(0,10000);
               fitFcn->SetParameter(11,100);
               
               fitFcn->FixParameter(1,mJpsi); // Mean x core
               fitFcn->FixParameter(2,0.07);
               fitFcn->FixParameter(12,ratsigma);
               
               if(SignalF=="NA60-MC"){
                   fitFcn->FixParameter(3,-0.4061);
                   fitFcn->FixParameter(4,0.2302);
                   fitFcn->FixParameter(5,1.2048);
                   fitFcn->FixParameter(6,0.0390);
                   fitFcn->FixParameter(7,2.0627);
                   fitFcn->FixParameter(8,0.1836);
                   fitFcn->FixParameter(9,1.2989);
                   fitFcn->FixParameter(10,0.0643);
               }
               
           }
        
        if(BackgroundF =="DoubleExpo"){
            backFcn1->SetParName(nextPar+0,"Norm_{TailLowM}");
            backFcn1->SetParName(nextPar+1,"Exp_{TailLowM}");
            backFcn1->SetParName(nextPar+2,"Norm_{TailHighM}");
            backFcn1->SetParName(nextPar+3,"Exp_{TailHighM}");
            
            backFcn1->SetParLimits(nextPar+0,0.01,20000);
            backFcn1->SetParLimits(nextPar+1,-20,20);
            backFcn1->SetParLimits(nextPar+2,0.001,100000000);
            backFcn1->SetParLimits(nextPar+3,0.001,50);
            
            backFcn1->FixParameter(nextPar+0,1);
            backFcn1->SetParameter(nextPar+1,0.01);
            backFcn1->SetParameter(nextPar+2,100);
            backFcn1->SetParameter(nextPar+3,10);
            
            TFitResultPtr resb = histo->Fit("backFcn1","SBMER","ep");
            double parb[20];
            backFcn1->GetParameters(parb);
            fitFcn->SetParName(nextPar+0,"Norm_{TailLowM}");
            fitFcn->SetParName(nextPar+1,"Exp_{TailLowM}");
            fitFcn->SetParName(nextPar+2,"Norm_{TailHighM}");
            fitFcn->SetParName(nextPar+3,"Exp_{TailHighM}");
            
            fitFcn->FixParameter(nextPar+0,parb[nextPar+0]);
            fitFcn->FixParameter(nextPar+1,parb[nextPar+1]);
            fitFcn->FixParameter(nextPar+2,parb[nextPar+2]);
            fitFcn->FixParameter(nextPar+3,parb[nextPar+3]);
            
        }
           if(BackgroundF =="POL1POL2"){
               
               
               backFcn1->SetParName(nextPar+0,"Norm_{Bkg}");
               backFcn1->SetParName(nextPar+1,"a_{1}");
               backFcn1->SetParName(nextPar+2,"b_{1}");
               backFcn1->SetParName(nextPar+3,"b_{2}");
               
               backFcn1->SetParLimits(nextPar+0,1000,40000);
               backFcn1->SetParLimits(nextPar+1,-0.5,0.5);
               backFcn1->SetParLimits(nextPar+2,-2,0);
               backFcn1->SetParLimits(nextPar+3,0,1);
               
               backFcn1->SetParameter(nextPar+0,3000);
               backFcn1->SetParameter(nextPar+1,-0.16);
               backFcn1->SetParameter(nextPar+2,-1.);
               backFcn1->SetParameter(nextPar+3,0.5);
            
               
               TFitResultPtr resb = histo->Fit("backFcn1","SBMER","ep");
               double parb[20];
               backFcn1->GetParameters(parb);
               fitFcn->SetParName(nextPar+0,"Norm_{Bkg}");
               fitFcn->SetParName(nextPar+1,"a_{1}");
               fitFcn->SetParName(nextPar+2,"b_{1}");
               fitFcn->SetParName(nextPar+3,"b_{2}");
               
               fitFcn->FixParameter(nextPar+0,parb[nextPar+0]);
               fitFcn->FixParameter(nextPar+1,parb[nextPar+1]);
               fitFcn->FixParameter(nextPar+2,parb[nextPar+2]);
               fitFcn->FixParameter(nextPar+3,parb[nextPar+3]);
               
           }
        
        if(BackgroundF =="VWG"){
            backFcn1->SetParName(nextPar+0,"Norm_{Bkg}");
            backFcn1->SetParName(nextPar+1,"Mean_{Bkg}");
            backFcn1->SetParName(nextPar+2,"Alpha");
            backFcn1->SetParName(nextPar+3,"Beta");
            
            backFcn1->SetParLimits(nextPar+0,0.0001,1000000);
            backFcn1->SetParLimits(nextPar+1,0.,2.5);
            backFcn1->SetParLimits(nextPar+2,0.00001,10);
            backFcn1->SetParLimits(nextPar+3,0.00001,10);
            
            backFcn1->SetParameter(nextPar+0,3000);
            backFcn1->SetParameter(nextPar+1,1.5);
            backFcn1->SetParameter(nextPar+2,0.5);
            backFcn1->SetParameter(nextPar+3,0.5);
            
            TFitResultPtr resb = histo->Fit("backFcn1","SBMER","ep");
            double parb[20];
            backFcn1->GetParameters(parb);
            fitFcn->SetParName(nextPar+0,"Norm_{Bkg}");
            fitFcn->SetParName(nextPar+1,"Mean_{Bkg}");
            fitFcn->SetParName(nextPar+2,"Alpha");
            fitFcn->SetParName(nextPar+3,"Beta");
            
            fitFcn->FixParameter(nextPar+0,parb[nextPar+0]);
            fitFcn->FixParameter(nextPar+1,parb[nextPar+1]);
            fitFcn->FixParameter(nextPar+2,parb[nextPar+2]);
            fitFcn->FixParameter(nextPar+3,parb[nextPar+3]);
        }
        
        if(BackgroundF =="Tchebychev"){
            backFcn1->SetParName(nextPar+0,"N");
            backFcn1->SetParName(nextPar+1,"c_{1}");
            backFcn1->SetParName(nextPar+2,"c_{2}");
            backFcn1->SetParName(nextPar+3,"c_{3}");
            backFcn1->SetParName(nextPar+4,"c_{4}");
            backFcn1->SetParName(nextPar+5,"min");
            backFcn1->SetParName(nextPar+6,"max");

            
            backFcn1->SetParLimits(nextPar+0,0,1000000);
            backFcn1->SetParLimits(nextPar+1,-1,1);
            backFcn1->SetParLimits(nextPar+2,-1,1);
            backFcn1->SetParLimits(nextPar+3,-1,1);
            backFcn1->SetParLimits(nextPar+4,-1,1);

            
            backFcn1->SetParameter(nextPar+0,195000);
            backFcn1->SetParameter(nextPar+1,-0.1);
            backFcn1->SetParameter(nextPar+2,0.15);
            backFcn1->SetParameter(nextPar+3,-0.01);
            backFcn1->SetParameter(nextPar+4,0);
            backFcn1->FixParameter(nextPar+5,min);
            backFcn1->FixParameter(nextPar+6,max);
            
            TFitResultPtr resb = histo->Fit("backFcn1","SBMER","ep");
            double parb[20];
            backFcn1->GetParameters(parb);
            fitFcn->SetParName(nextPar+0,"N");
            fitFcn->SetParName(nextPar+1,"c_{1}");
            fitFcn->SetParName(nextPar+2,"c_{2}");
            fitFcn->SetParName(nextPar+3,"c_{3}");
            fitFcn->SetParName(nextPar+4,"c_{4}");
           // fitFcn->SetParName(nextPar+5,"c_{5}");
            
            fitFcn->FixParameter(nextPar+0,parb[nextPar+0]);
            fitFcn->FixParameter(nextPar+1,parb[nextPar+1]);
            fitFcn->FixParameter(nextPar+2,parb[nextPar+2]);
            fitFcn->FixParameter(nextPar+3,parb[nextPar+3]);
            fitFcn->FixParameter(nextPar+4,parb[nextPar+4]);
           // fitFcn->FixParameter(nextPar+5,parb[nextPar+5]);
        }
        
        res2 = histo->Fit("fitFcn","SBMER","ep");
     //   res2 = histo->Fit("fitFcn","SBMER","ep");
           fitFcn->ReleaseParameter(1);
        res2 = histo->Fit("fitFcn","SBMER","ep");
         fitFcn->ReleaseParameter(2);
        res2 = histo->Fit("fitFcn","SBMER","ep");
         if(SignalF=="CB2-Run2" || SignalF=="CB2-MC" || SignalF=="CB2-FREE"){
             fitFcn->ReleaseParameter(7);
             res2 = histo->Fit("fitFcn","SBMER","ep");
         }
         if(SignalF=="NA60-MC" || SignalF=="NA60-FREE"){
             fitFcn->ReleaseParameter(11);
             res2 = histo->Fit("fitFcn","SBMER","ep");
         }
         if(res2->CovMatrixStatus() !=3){
            res2 = histo->Fit("fitFcn","SBMER","ep");
            res2 = histo->Fit("fitFcn","SBMER","ep");
         }
        
//        for(int idx = 0; idx<NumberOfParameters(SignalF,BackgroundF);idx++){
//            if(SignalF=="CB2-Run2" || SignalF=="CB2-MC"){
//               if((idx>=0 && idx<=6)|| idx==8){ //if((idx>=3 && idx<=6)||idx==8){
//                    continue;
//                }
//            }
//            fitFcn->ReleaseParameter(idx);
//        }
//        fitFcn->SetParLimits(9,0.001,10000000);
//        fitFcn->SetParLimits(10,0.,2.5);
//        fitFcn->SetParLimits(11,0.001,100);
//        fitFcn->SetParLimits(12,0.001,100);
//        res2 = histo->Fit("fitFcn","SBMER","ep");
        
    }
    
    
    
    gStyle->SetOptStat("ne");
    gStyle->SetOptFit(1012);
    TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
    st->SetX1NDC(0.8); //new x start position
    st->SetY1NDC(0.3); //new x end position
//    fitFcn->ReleaseParameter(3);
//       res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(5);
//   res = histo->Fit("fitFcn","SBMER","ep");
//   fitFcn->ReleaseParameter(4);
//    res = histo->Fit("fitFcn","SBMER","ep");
//    fitFcn->ReleaseParameter(6);
//   res = histo->Fit("fitFcn","SBMER","ep");
//    res = histo->Fit("fitFcn","SBMER","ep");

    
    TF1 *backFcn = NULL;
    
    backFcn = new TF1("backFcn",UsedBackground,min,max,numParameters);
    
    
    gStyle->SetErrorX(0.);

    histo->SetMarkerColor(kBlack);
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(0.25);
    histo->SetLineColor(kBlack);
    histo->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/c^{2})");
    histo->GetYaxis()->SetTitle("Entries per 10 MeV/c^{2}");
    

    
    TF1 *signalFcnJPsi = NULL;
    TF1 *signalFcnPsi2S = NULL;
    
    signalFcnJPsi = new TF1("signalFcnJPsi",UsedJPsiSignal,0,max,numParameters-NumberOfParametersBkg(BackgroundF));
    signalFcnPsi2S = new TF1("signalFcnPsi2S",UsedPsi2SSignal,0,max,numParameters-NumberOfParametersBkg(BackgroundF));
    
    backFcn->SetLineStyle(9);
    backFcn->SetLineWidth(2);
   backFcn->SetLineColor(kRed+1);
   TPaveText *pave = new TPaveText(0.5,0.5,0.7,0.6,"brNDC");
    signalFcnJPsi->SetLineStyle(9);
    signalFcnJPsi->SetLineWidth(2);
   signalFcnJPsi->SetLineColor(kBlue);
   signalFcnJPsi->SetNpx(500);
    signalFcnPsi2S->SetLineStyle(9);
    signalFcnPsi2S->SetLineWidth(2);
    signalFcnPsi2S->SetLineColor(kGreen-3);
    signalFcnPsi2S->SetNpx(500);

   fitFcn->GetParameters(par);
   signalFcnJPsi->SetParameters(par);
    signalFcnPsi2S->SetParameters(par);
    
    Double_t integral;
    Double_t integralerror;
    Double_t integralPsi2s;
    Double_t integralerrorPsi2s;
    
    cout << "===== Fit results =====" <<endl;
    if(!isFitRetried){
        cout << "First try" <<endl;
   integral = (signalFcnJPsi->Integral(0,5))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
        integralPsi2s = (signalFcnPsi2S->Integral(0,10))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
    auto covtot = res->GetCovarianceMatrix();
    auto covsgn = covtot.GetSub(0,numParameters-NumberOfParametersBkg(BackgroundF)-1,0,numParameters-NumberOfParametersBkg(BackgroundF)-1);
    std::cout << "Matrice totale" <<endl;
    covtot.Print();
    std::cout << "Matrice réduite" <<endl;
    covsgn.Print();
    std::cout << "STATUS COV " << res->CovMatrixStatus() <<endl;
    integralerror = (signalFcnJPsi->IntegralError(0,5,signalFcnJPsi->GetParameters(), res->GetCovarianceMatrix().GetSub(0,numParameters-NumberOfParametersBkg(BackgroundF)-1,0,numParameters-NumberOfParametersBkg(BackgroundF)-1).GetMatrixArray() ))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
    integralerrorPsi2s = (signalFcnPsi2S->IntegralError(0,10,signalFcnPsi2S->GetParameters(), res->GetCovarianceMatrix().GetSub(0,numParameters-NumberOfParametersBkg(BackgroundF)-1,0,numParameters-NumberOfParametersBkg(BackgroundF)-1).GetMatrixArray() ))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
    std::cout << "Erreur integrale Jpsi " << integralerror <<endl;
        std::cout << "Erreur integrale Psi2s " << integralerrorPsi2s <<endl;
    
    std::cout << "Fitted " << histoname << std::endl;
    std::cout << "Nb JPsi total: " << integral << std::endl;
 //   std::cout << "integral error: " << integralerror << std::endl;
    }
    if(isFitRetried){
           cout << "Second try" <<endl;
      integral = (signalFcnJPsi->Integral(0,5))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
        integralPsi2s = (signalFcnPsi2S->Integral(0,10))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
       auto covtot = res2->GetCovarianceMatrix();
       auto covsgn = covtot.GetSub(0,numParameters-NumberOfParametersBkg(BackgroundF)-1,0,numParameters-NumberOfParametersBkg(BackgroundF)-1);
       std::cout << "Matrice totale" <<endl;
       covtot.Print();
       std::cout << "Matrice réduite" <<endl;
       covsgn.Print();
       std::cout << "STATUS COV " << res2->CovMatrixStatus() <<endl;
      integralerror = (signalFcnJPsi->IntegralError(0,5,signalFcnJPsi->GetParameters(), res->GetCovarianceMatrix().GetSub(0,numParameters-NumberOfParametersBkg(BackgroundF)-1,0,numParameters-NumberOfParametersBkg(BackgroundF)-1).GetMatrixArray() ))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
      integralerrorPsi2s = (signalFcnPsi2S->IntegralError(0,10,signalFcnPsi2S->GetParameters(), res->GetCovarianceMatrix().GetSub(0,numParameters-NumberOfParametersBkg(BackgroundF)-1,0,numParameters-NumberOfParametersBkg(BackgroundF)-1).GetMatrixArray() ))*NbinsDimuInvMass/(MaxInvMass-MinInvMass);
       std::cout << "Erreur integrale " << integralerror <<endl;
       
       std::cout << "Fitted " << histoname << std::endl;
       std::cout << "Nb JPsi total: " << integral << std::endl;
    //   std::cout << "integral error: " << integralerror << std::endl;
       }
    
    
   signalFcnJPsi->Draw("same");
    signalFcnPsi2S->Draw("same");
   backFcn->SetParameters(par);
   backFcn->Draw("same");
    
    
   // draw the legend
    char str[50];
    string fitPerformed = SignalF+BackgroundF;
    int n = fitPerformed.length();
    char char_fitPerformed[n + 1];
       strcpy(char_fitPerformed, fitPerformed.c_str());
    
    sprintf(str, "N_{J/#psi} = %i +/- %i  ", int(integral), int(integralerror));
   pave->AddText(str);
    sprintf(str, "N_{#Psi(2S)} = %i +/- %i  ", int(integralPsi2s), int(integralerrorPsi2s));
    pave->AddText(str);
    pave->AddText(char_fitPerformed);
    pave->SetTextFont(42);
    pave->SetTextSize(0.04);
    pave->SetBorderSize(0);
    pave->SetFillStyle(0);
    sprintf(str, "M_{#Psi(2S)} = %f, Sig_{#Psi(2S)} = %f", int((mJpsi+mDiff)*1000)/1000., int(par[numParameters-NumberOfParametersBkg(BackgroundF)-2]*par[2]*1000)/1000.);
  //  pave->AddText(str);
  // pave->Draw();
   TLegend *legend=new TLegend(0.4,0.6,0.6,0.8);
   legend->SetTextFont(42);
   legend->SetTextSize(0.035);
    legend->SetFillColorAlpha(kWhite, 0.);
    legend->SetBorderSize(0);
    legend->AddEntry(histo,"Data","lpe");
    Char_t message[80];
    sprintf(message,"Total fit");
    //sprintf(message,"Fit, #chi^{2}/NDF = %.2f / %d = %.2f",fitFcn->GetChisquare(),fitFcn->GetNDF(),fitFcn->GetChisquare()/fitFcn->GetNDF());
     legend->AddEntry(fitFcn,message);//<<<
    
    cout << "First try: " << res->CovMatrixStatus() <<endl;
    if(isFitRetried){
        cout<<"Second try: " << res2->CovMatrixStatus()<<endl;
    }
    if(res->CovMatrixStatus() == 3){
        cout << "Status " << res->IsValid() <<endl;
           }
    else if(isFitRetried && res2->CovMatrixStatus() == 3){
        cout <<"The fit is ok on second try"<<endl;
               sprintf(message,"The fit is ok on second try");
           }
       else{
           cout <<"The fit is a failure"<<endl;
           sprintf(message,"The fit is a failure");
           legend->AddEntry(fitFcn,message);
       }
   legend->AddEntry(backFcn,"Background","l");
   legend->AddEntry(signalFcnJPsi,"J/#psi","l");
    legend->AddEntry(signalFcnPsi2S,"#Psi(2S)","l");
   legend->Draw();
    
    cout << "Test1"<<endl;
//    cout << "Histogramme: " << histoname << endl;
//    for(int j=0; j<10; j++){
//        cout << "Bin de masse numero " << j << " - Signal : " << (signalFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 << " et Background : " << (backFcn->Integral(2 + j*0.25,2 + (j+1)*0.25))/0.01 <<endl;
//    }
    // myfiletxt.close();
    cout << "Test2"<<endl;
    
    if(!isFitRetried){
        cout << "Test3"<<endl;
        return res;
    }
    if(isFitRetried){
        cout << "Test4"<<endl;
        return res2;
    }
    else{
        cout << "Test5"<<endl;
        return 0;
    }
}

void FcnCombinedAllMass(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
   // cout << "acutelptbin " << actuelptbin<<endl;
  TAxis *xaxis1  = hPtWrtMassInvSliced[0][actuelptbin]->GetXaxis();
  TAxis *xaxis2  = V2JPsiTklPtBinned[actuelptbin]->GetXaxis();

  int nbinX1 = hPtWrtMassInvSliced[0][actuelptbin]->GetNbinsX();
  int nbinX2 = V2JPsiTklPtBinned[actuelptbin]->GetNbinsX();

  double chi2 = 0;
  double x[1];
  double tmp;
  npfits = 0;
  for (int ix = 1; ix <= nbinX1; ++ix) {
    x[0] = xaxis1->GetBinCenter(ix);
//      if(x[0]<1.5 || x[0]>4.5){
//          continue;
//      }
      if ( hPtWrtMassInvSliced[0][actuelptbin]->GetBinError(ix) > 0 ) {
        tmp = (hPtWrtMassInvSliced[0][actuelptbin]->GetBinContent(ix) - MassFunction(x,p))/hPtWrtMassInvSliced[0][actuelptbin]->GetBinError(ix);
        chi2 += tmp*tmp;
        npfits++;
      }
  }
  for (int ix = 1; ix <= nbinX2; ++ix) {
     x[0] = xaxis2->GetBinCenter(ix);
//      if(x[0]<1.5 || x[0]>4.5){
//          continue;
//      }
      if ( V2JPsiTklPtBinned[actuelptbin]->GetBinError(ix) > 0 ) {
        tmp = (V2JPsiTklPtBinned[actuelptbin]->GetBinContent(ix) - FourierV2_WrtInvMass(x,p))/V2JPsiTklPtBinned[actuelptbin]->GetBinError(ix);
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
    
//    if(cent<=20){
//        return 0;
//    }
//    else if(cent<=40){
//        return 1;
//    }
//    else{
//        return 2;
//    }
}


int NumberOfParameters(string SignalF, string BackgroundF){
    
    int numberParameters = 1; //Start at 1 to account for fixed sigmaratio parameter
    
    if(SignalF=="CB2-Run2" || SignalF=="CB2-MC" || SignalF=="CB2-FREE"){
        numberParameters += 8;
    }
    if(SignalF=="NA60-MC" || SignalF=="NA60-FREE"){
        numberParameters += 12;
    }
    if(BackgroundF =="VWG"){
        numberParameters += 4;
    }
    if(BackgroundF =="DoubleExpo" || BackgroundF =="POL1POL2"){
        numberParameters += 4;
    }
    if(BackgroundF =="Tchebychev"){
        numberParameters += 7;
    }
    
    return numberParameters;
    
}

int NumberOfParametersBkg(string BackgroundF){
    
    int numberParameters = 0;
    
    if(BackgroundF =="VWG"){
        numberParameters += 4;
    }
    if(BackgroundF =="DoubleExpo" || BackgroundF =="POL1POL2"){
        numberParameters += 4;
    }
    if(BackgroundF =="Tchebychev"){
        numberParameters += 7;
    }
    
    
    return numberParameters;
    
}

int NumberOfParametersBkgV2(string BackgroundV2F){
    
    int numberParametersV2 = 0;
    
    if(BackgroundV2F =="Pol2"){
        numberParametersV2 += 3;
    }
    if(BackgroundV2F =="Pol1"){
        numberParametersV2 += 2;
    }
    if(BackgroundV2F =="Pol3"){
        numberParametersV2 += 4;
    }
    
    
    return numberParametersV2;
    
}

Double_t ExpBkgHole(Double_t *x,Double_t *par)

{
    if(x[0]>2.2 && x[0]<3.9){
        TF1::RejectPoint();
        return 0;
    }
    return exp(x[0]*par[1]*(-1)) + par[2]*(exp(x[0]*par[3]*(-1))); }

Double_t Pol1Pol2Hole(Double_t *x,Double_t *par)

{if(x[0]>2.3 && x[0]<3.9){
        TF1::RejectPoint();
        return 0;
    }
    
    return par[0]*(((1+par[1]*x[0]))/(1+par[2]*x[0] + par[3]*(x[0]*x[0]))); }

Double_t VWGaussianHole(Double_t *x,Double_t *par)

{ double sigma = par[2] + par[3]*((x[0]-par[1])/par[1]);
    if(x[0]>2.4 && x[0]<3.9){
        TF1::RejectPoint();
        return 0;
    }
    return par[0]*exp(-1.0*pow(x[0]-par[1],2)/(2*pow(sigma,2))); }

Double_t TchebychevHole(Double_t *x,Double_t *par)

{Double_t xx = (x[0]-par[5])/(par[6]-par[5]);
    
    Double_t poly = SmolTcheby(xx,par, 0);;
    if(x[0]>2.2 && x[0]<3.9){
        TF1::RejectPoint();
        return 0;
    }
for(int index=1; index<5;index++){
    poly += par[index]*SmolTcheby(xx,par, index);
}
poly *= par[0];
return poly;
}

string GetParameterInfo(int idx, string SignalF, string BackgroundF, string BackgroundV2F, string prefix){
    int nextParameter = 0;
    if(SignalF=="CB2-Run2" || SignalF=="CB2-MC" || SignalF=="CB2-FREE"){
    nextParameter = 9;
        if(idx==0){
            return "Norm_{J/#psi}";
        }
        if(idx==1){
            return "M_{J/#psi}";
        }
        if(idx==2){
            return "Sigma_{J/#psi}";
        }
        if(idx==3){
            return "a_{1}";
        }
        if(idx==4){
            return "n_{1}";
        }
        if(idx==5){
            return "a_{2}";
        }
        if(idx==6){
            return "n_{2}";
        }
        if(idx==7){
            return "Norm_{#Psi(2S)}";
        }
        if(idx==8){
            return "Sigma_{ratio}";
        }
    }
    if(SignalF=="NA60-MC" || SignalF=="NA60-FREE"){
        nextParameter = 13;
        if(idx==0){
            return "Norm_{J/#psi}";
        }
        if(idx==1){
            return "M_{J/#psi}";
        }
        if(idx==2){
            return "Sigma_{J/#psi}";
        }
        if(idx==3){
            return "a_{L}";
        }
        if(idx==4){
            return "p1_{L}";
        }
        if(idx==5){
            return "p2_{L}";
        }
        if(idx==6){
            return "p3_{L}";
        }
        if(idx==7){
            return "a_{R}";
        }
        if(idx==8){
            return "p1_{R}";
        }
        if(idx==9){
            return "p2_{R}";
        }
        if(idx==10){
            return "p3_{R}";
        }
        if(idx==11){
            return "Norm_{#Psi(2S)}";
        }
        if(idx==12){
            return "Sigma_{ratio}";
        }
    }
    if(BackgroundF =="DoubleExpo"){
        if(idx==nextParameter+0){
            return "Norm_{TailLowM}";
        }
        if(idx==nextParameter+1){
            return "Exp_{TailLowM}";
        }
        if(idx==nextParameter+2){
            return "Norm_{TailHighM}";
        }
        if(idx==nextParameter+3){
            return "Exp_{TailHighM}";
        }
        nextParameter += 4;
    }
    if(BackgroundF =="POL1POL2"){
        if(idx==nextParameter+0){
            return "Norm_{Bkg}";
        }
        if(idx==nextParameter+1){
            return "a_{1}";
        }
        if(idx==nextParameter+2){
            return "b_{1}";
        }
        if(idx==nextParameter+3){
            return "b_{2}";
        }
        nextParameter += 4;
    }
    if(BackgroundF =="VWG"){
        if(idx==nextParameter+0){
            return "Norm_{Bkg}";
        }
        if(idx==nextParameter+1){
            return "Mean_{Bkg}";
        }
        if(idx==nextParameter+2){
            return "Alpha";
        }
        if(idx==nextParameter+3){
            return "Beta";
        }
        nextParameter += 4;
    }
    if(BackgroundF =="Tchebychev"){
        if(idx==nextParameter+0){
            return "N";
        }
        if(idx==nextParameter+1){
            return "c_{1}";
        }
        if(idx==nextParameter+2){
            return "c_{2}";
        }
        if(idx==nextParameter+3){
            return "c_{3}";
        }
        if(idx==nextParameter+4){
            return "c_{4}";
        }
        if(idx==nextParameter+5){
            return "min";
        }
        if(idx==nextParameter+6){
            return "max";
        }
        nextParameter += 7;
    }
    
    if(idx==nextParameter){
        return prefix + " J/#psi";
    }
    nextParameter++;
    
    if(BackgroundV2F =="Pol1"){
        if(idx==nextParameter+0){
            return prefix + " Bkg_1";
        }
        if(idx==nextParameter+1){
            return prefix + " Bkg_0";
        }
    }
    if(BackgroundV2F =="Pol2"){
        if(idx==nextParameter+0){
            return prefix + " Bkg_2";
        }
        if(idx==nextParameter+1){
            return prefix + " Bkg_1";
        }
        if(idx==nextParameter+2){
            return prefix + " Bkg_0";
        }
    }
    if(BackgroundV2F =="Pol3"){
        if(idx==nextParameter+0){
            return prefix + " Bkg_3";
        }
        if(idx==nextParameter+1){
            return prefix + " Bkg_2";
        }
        if(idx==nextParameter+2){
            return prefix + " Bkg_1";
        }
        if(idx==nextParameter+3){
            return prefix + " Bkg_0";
        }
    }
    return 0;

}
