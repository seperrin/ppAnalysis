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
 #include <TGraphErrors.h>
 #include <TMatrixD.h>
 #include <TLorentzVector.h>
 #include <vector>
#include <sys/types.h>
#include <sys/stat.h>

 //#include "TTreeReader.h"
 //#include "Event.h"

 # include "TF1.h"
 # include "TF2.h"
 # include "TProfile.h"
 # include "TLegend.h"
# include "TGraph2D.h"
 # include "TVirtualFitter.h"
 # include "TBackCompFitter.h"
 # include "TGraphErrors.h"
 # include "TFitter.h"
 # include "TMinuit.h"
 # include "TRandom.h"
 # include <iostream>
 # include <fstream>
# include <sstream>
 # include <string>
# include <iterator>
 # include "TH1F.h"
 # include "TH1D.h"
 # include "TH2F.h"
 # include "Scripts/AliAnalysisTaskMyMuonTree_AOD.h"

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

void SystematicsDimu(Char_t ListeCSVFiles[100], Char_t SystematicsFocus[50]);
bool isRowSuitableForSystematics(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]);
void CentralityDimu(Char_t ListeCSVFiles[100], Char_t SystematicsFocus[50]);
void CentralityPlotDimu(Char_t ListeCSVFiles[100], Char_t SystematicsFocus[50]);
bool isRowSuitableForCentrality(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]);

int NbPtBins = 6;
double PtBins[] = {0,2,3,4,6,8,12};
//int NbPtBins = 2;
//double PtBins[] = {0,1,12};


void SystematicsDimu(Char_t listeCSVFiles[100], Char_t SystematicsFocus[50])
{
    
    std::ifstream file(listeCSVFiles);
    CSVRow row;
    std::string str;
    
    cout << "LOLILOL"<<endl;
    
    double v2pPbSPDT[10][100]={0};
    double errv2pPbSPDT[10][100]={0};
    double v2pPbSPDC[10][100]={0};
    double errv2pPbSPDC[10][100]={0};
    double v2pPbV0M[10][100]={0};
    double errv2pPbV0M[10][100]={0};
    
    double v2ATLASSPDT[10][100]={0};
    double errv2ATLASSPDT[10][100]={0};
    double v2ATLASSPDC[10][100]={0};
    double errv2ATLASSPDC[10][100]={0};
    double v2ATLASV0M[10][100]={0};
    double errv2ATLASV0M[10][100]={0};
    
    cout << "LOLILOL"<<endl;
    
    double v2ppSystV0M[10]={0};
    double errv2ppSystV0M[10]={0};
    double v2ppSystV0Malt[10][100]={0};
    double errv2ppSystV0Malt[10][100]={0};
    double v2ppSystV0Malto[10]={0};
    double errv2ppSystV0Malto[10]={0};
    
    double v2ppSystSPDT[10]={0};
    double errv2ppSystSPDT[10]={0};
    double v2ppSystSPDTalt[10][100]={0};
    double errv2ppSystSPDTalt[10][100]={0};
    double v2ppSystSPDTalto[10]={0};
    double errv2ppSystSPDTalto[10]={0};
    
    double v2ppSystSPDC[10]={0};
    double errv2ppSystSPDC[10]={0};
    double v2ppSystSPDCalt[10][100]={0};
    double errv2ppSystSPDCalt[10][100]={0};
    double v2ppSystSPDCalto[10]={0};
    double errv2ppSystSPDCalto[10]={0};
    
    double v2ppATLASSystV0M[10]={0};
    double errv2ppATLASSystV0M[10]={0};
    double v2ppATLASSystV0Malt[10][100]={0};
    double errv2ppATLASSystV0Malt[10][100]={0};
    double v2ppATLASSystV0Malto[10]={0};
    double errv2ppATLASSystV0Malto[10]={0};
    
    double v2ppATLASSystSPDT[10]={0};
    double errv2ppATLASSystSPDT[10]={0};
    double v2ppATLASSystSPDTalt[10][100]={0};
    double errv2ppATLASSystSPDTalt[10][100]={0};
    double v2ppATLASSystSPDTalto[10]={0};
    double errv2ppATLASSystSPDTalto[10]={0};
    
    double v2ppATLASSystSPDC[10]={0};
    double errv2ppATLASSystSPDC[10]={0};
    double v2ppATLASSystSPDCalt[10][100]={0};
    double errv2ppATLASSystSPDCalt[10][100]={0};
    double v2ppATLASSystSPDCalto[10]={0};
    double errv2ppATLASSystSPDCalto[10]={0};
    
    cout << "LOLILOL"<<endl;
    
    
    int IndexOfVariable;
    int IndexOfVariableSignal;
    int IndexOfVariableBackground;
    int IndexOfVariableRatSigma;
    int IndexOfVariableMinMass;
    int IndexOfVariableMaxMass;
    int IndexOfVariableMinV2;
    int IndexOfVariableMaxV2;
    int IndexOfStorage;
    
    string variableDPhiCut[5] = {"1mrad", "2mrad", "5mrad", "10mrad", "None"};
    string variableZvtxCut[3] = {"8","10","12"};
    string variableEtaMin[3] = {"1.3","1.5","1.7"};
    string variableEtaMax[3] = {"4.8","5","5.2"};
    string variableEMNorm[3] = {"Method1","Method2","Method3"};
    string variableSummationZvtx[3] = {"Method1c","Method1a","Method2"};
    string variablePooling[3] = {"None","PhUniOpti"};
    string variableSignal[3] = {"CB2-Run2", "CB2-MC","NA60-MC"};
    string variableBackground[4] = {"DoubleExpo", "POL1POL2", "VWG","Tchebychev"};
    string variableRatSigma[2] = {"1.05","1"};
    string variableMinMass[3] = {"2","2.3","2.4"};
    string variableMaxMass[3] = {"5","4.9","4.7"};
    string variableExtractionMethod[3] = {"Method1","Method2","Method3"};
    string variableBackgroundV2[3] = {"Pol1","Pol2","Pol3"};
    string variableMinV2[2] = {"1","1.5"};
    string variableMaxV2[2] = {"5","4.5"};
    string variableZYAM[2] = {"ZYAM","noZYAM"};
    string variableATLAS[2] = {"Cvetan","ATLAS"};
    cout << "LOLILOL"<<endl;
    
    string naming[100];
    
    cout << "LOLILOL"<<endl;
    
    int NombreEssais= 0;
    int default_idx = 0;
    
    bool isClassicVariation = kFALSE;
    
    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
        IndexOfVariable = 0;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
        IndexOfVariable = 4 + NbPtBins + 1 + 3;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "EtaMin") == 0){
        IndexOfVariable = 4 + NbPtBins + 1 + 1;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "EtaMax") == 0){
        IndexOfVariable = 4 + NbPtBins + 1 + 2;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "EMNorm") == 0){
        IndexOfVariable = 4 + NbPtBins + 1 + 4;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
        IndexOfVariable = 4 + NbPtBins + 1 + 8;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "Pooling") == 0){
        IndexOfVariable = 4 + NbPtBins + 1 + 7;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "InvMass") == 0){
      IndexOfVariableSignal = 4 + NbPtBins + 1 + 14;
      IndexOfVariableBackground = 4 + NbPtBins + 1 + 15;
      IndexOfVariableRatSigma = 4 + NbPtBins + 1 + 16;
      IndexOfVariableMinMass = 4 + NbPtBins + 1 + 12;
      IndexOfVariableMaxMass = 4 + NbPtBins + 1 + 13;
        
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
        IndexOfVariable = 4 + NbPtBins + 1 + 9;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "RangeV2") == 0){
        IndexOfVariableMinV2 = 4 + NbPtBins + 1 + 10;
        IndexOfVariableMaxV2 = 4 + NbPtBins + 1 + 11;
        isClassicVariation = kTRUE;
    }
    
    cout << "We will run a simple read of csv files" <<endl;
    
    while (std::getline(file, str))
    {
        std::ifstream filecsv(str);
        while(filecsv >> row)
        {
           // cout << "File looked at " << str<<endl;
            bool isVariationFound = kFALSE;
                if(isRowSuitableForSystematics(row, SystematicsFocus, "V0MPercentile")){
                    cout << "Keeping file " << str <<endl;
                    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
                        for(int idx=0; idx<5; idx++){
                            cout << "idx "<< idx << " variableDPhiCut[idx] "<< variableDPhiCut[idx]<<endl;
                            cout << "idx "<< idx << " row[IndexOfVariable] "<< row[IndexOfVariable]<<endl;
                            if (variableDPhiCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                cout << "idx "<< idx << "  variableDPhiCut[idx] "<<variableDPhiCut[idx]<< " row[IndexOfVariable] "<< row[IndexOfVariable]<<endl;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableZvtxCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EtaMin") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaMin[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EtaMax") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaMax[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EMNorm") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEMNorm[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableSummationZvtx[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "Pooling") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variablePooling[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "InvMass") == 0){
                        int indexSignal = 0;
                        int indexBackground = 0;
                        int indexRatSigma = 0;
                        int indexMass = 0;
                      for(int idx=0; idx<3; idx++){
                          if (variableSignal[idx] == row[IndexOfVariableSignal]){
                              indexSignal = idx;
                          }
                      }
                        for(int idx=0; idx<4; idx++){
                            if (variableBackground[idx] == row[IndexOfVariableBackground]){
                                indexBackground = idx;
                            }
                        }
                        for(int idx=0; idx<2; idx++){
                            if (variableRatSigma[idx] == row[IndexOfVariableRatSigma]){
                                indexRatSigma = idx;
                            }
                        }
                        for(int idx=0; idx<3; idx++){
                            if (variableMinMass[idx] == row[IndexOfVariableMinMass] && variableMaxMass[idx] == row[IndexOfVariableMaxMass]){
                                indexMass = idx;
                            }
                        }
                        IndexOfStorage = indexMass + 3*indexRatSigma + 6*indexBackground + 24*indexSignal;
                        isVariationFound = kTRUE;
                    }
                    if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableBackgroundV2[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "RangeV2") == 0){
                      for(int idx=0; idx<2; idx++){
                          if (variableMinV2[idx] == row[IndexOfVariableMinV2] && variableMaxV2[idx] == row[IndexOfVariableMaxV2]){
                              IndexOfStorage = idx;
                              isVariationFound = kTRUE;
                              continue;
                          }
                      }
                    }
                    
                    for(int ptbin=0; ptbin<NbPtBins; ptbin++){
                        if(isVariationFound){//FIx indices
                            v2pPbV0M[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbV0M[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2ATLASV0M[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASV0M[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ZYAM") == 0){
                            v2pPbV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2pPbV0M[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 23 + 18*ptbin]));
                            errv2pPbV0M[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 24 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ATLAS") == 0){
                            v2pPbV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2ATLASV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ExtractionMethod") == 0){//Fixme
                            v2pPbV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2pPbV0M[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 19 + 18*ptbin]));
                            errv2pPbV0M[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 20 + 18*ptbin]));
                            v2pPbV0M[ptbin][2] = stod(std::string(row[4 + NbPtBins + 1 + 21 + 18*ptbin]));
                            errv2pPbV0M[ptbin][2] = stod(std::string(row[4 + NbPtBins + 1 + 22 + 18*ptbin]));
                            v2ATLASV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASV0M[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                            v2ATLASV0M[ptbin][1] = stod(std::string(row[5 + NbPtBins + 1 + 16 + 18*NbPtBins + 2*ptbin]));
                            errv2ATLASV0M[ptbin][1] = stod(std::string(row[6 + NbPtBins + 1 + 16 + 18*NbPtBins + 2*ptbin]));
                        }
                    }
                }
                else if(isRowSuitableForSystematics(row, SystematicsFocus, "SPDTrackletsPercentile")){
                    cout << "Keeping file " << str <<endl;
                    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
                        for(int idx=0; idx<5; idx++){
                            cout << "idx "<< idx << " variableDPhiCut[idx] "<< variableDPhiCut[idx]<<endl;
                            cout << "idx "<< idx << " row[IndexOfVariable] "<< row[IndexOfVariable]<<endl;
                            if (variableDPhiCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                cout << "idx "<< idx << "  variableDPhiCut[idx] "<<variableDPhiCut[idx]<< " row[IndexOfVariable] "<< row[IndexOfVariable]<<endl;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableZvtxCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EtaMin") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaMin[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EtaMax") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaMax[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EMNorm") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEMNorm[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableSummationZvtx[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "Pooling") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variablePooling[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "InvMass") == 0){
                        int indexSignal = 0;
                        int indexBackground = 0;
                        int indexRatSigma = 0;
                        int indexMass = 0;
                      for(int idx=0; idx<3; idx++){
                          if (variableSignal[idx] == row[IndexOfVariableSignal]){
                              indexSignal = idx;
                          }
                      }
                        for(int idx=0; idx<4; idx++){
                            if (variableBackground[idx] == row[IndexOfVariableBackground]){
                                indexBackground = idx;
                            }
                        }
                        for(int idx=0; idx<2; idx++){
                            if (variableRatSigma[idx] == row[IndexOfVariableRatSigma]){
                                indexRatSigma = idx;
                            }
                        }
                        for(int idx=0; idx<3; idx++){
                            if (variableMinMass[idx] == row[IndexOfVariableMinMass] && variableMaxMass[idx] == row[IndexOfVariableMaxMass]){
                                indexMass = idx;
                            }
                        }
                        IndexOfStorage = indexMass + 3*indexRatSigma + 6*indexBackground + 24*indexSignal;
                        isVariationFound = kTRUE;
                    }
                    if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableBackgroundV2[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "RangeV2") == 0){
                      for(int idx=0; idx<2; idx++){
                          if (variableMinV2[idx] == row[IndexOfVariableMinV2] && variableMaxV2[idx] == row[IndexOfVariableMaxV2]){
                              IndexOfStorage = idx;
                              isVariationFound = kTRUE;
                              continue;
                          }
                      }
                    }
                    
                    for(int ptbin=0; ptbin<NbPtBins; ptbin++){
                        if(isVariationFound){//FIx indices
                            v2pPbSPDT[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDT[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2ATLASSPDT[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASSPDT[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ZYAM") == 0){
                            v2pPbSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2pPbSPDT[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 23 + 18*ptbin]));
                            errv2pPbSPDT[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 24 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ATLAS") == 0){
                            v2pPbSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2ATLASSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ExtractionMethod") == 0){//Fixme
                            v2pPbSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2pPbSPDT[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 19 + 18*ptbin]));
                            errv2pPbSPDT[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 20 + 18*ptbin]));
                            v2pPbSPDT[ptbin][2] = stod(std::string(row[4 + NbPtBins + 1 + 21 + 18*ptbin]));
                            errv2pPbSPDT[ptbin][2] = stod(std::string(row[4 + NbPtBins + 1 + 22 + 18*ptbin]));
                            v2ATLASSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASSPDT[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                            v2ATLASSPDT[ptbin][1] = stod(std::string(row[5 + NbPtBins + 1 + 16 + 18*NbPtBins + 2*ptbin]));
                            errv2ATLASSPDT[ptbin][1] = stod(std::string(row[6 + NbPtBins + 1 + 16 + 18*NbPtBins + 2*ptbin]));
                        }
                    }
                }
                else if(isRowSuitableForSystematics(row, SystematicsFocus, "SPDClustersPercentile")){
                    cout << "Keeping file " << str <<endl;
                    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
                        for(int idx=0; idx<5; idx++){
                            cout << "idx "<< idx << " variableDPhiCut[idx] "<< variableDPhiCut[idx]<<endl;
                            cout << "idx "<< idx << " row[IndexOfVariable] "<< row[IndexOfVariable]<<endl;
                            if (variableDPhiCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                cout << "idx "<< idx << "  variableDPhiCut[idx] "<<variableDPhiCut[idx]<< " row[IndexOfVariable] "<< row[IndexOfVariable]<<endl;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableZvtxCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EtaMin") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaMin[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EtaMax") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaMax[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "EMNorm") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEMNorm[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableSummationZvtx[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "Pooling") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variablePooling[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "InvMass") == 0){
                        int indexSignal = 0;
                        int indexBackground = 0;
                        int indexRatSigma = 0;
                        int indexMass = 0;
                      for(int idx=0; idx<3; idx++){
                          if (variableSignal[idx] == row[IndexOfVariableSignal]){
                              indexSignal = idx;
                          }
                      }
                        for(int idx=0; idx<4; idx++){
                            if (variableBackground[idx] == row[IndexOfVariableBackground]){
                                indexBackground = idx;
                            }
                        }
                        for(int idx=0; idx<2; idx++){
                            if (variableRatSigma[idx] == row[IndexOfVariableRatSigma]){
                                indexRatSigma = idx;
                            }
                        }
                        for(int idx=0; idx<3; idx++){
                            if (variableMinMass[idx] == row[IndexOfVariableMinMass] && variableMaxMass[idx] == row[IndexOfVariableMaxMass]){
                                indexMass = idx;
                            }
                        }
                        IndexOfStorage = indexMass + 3*indexRatSigma + 6*indexBackground + 24*indexSignal;
                        isVariationFound = kTRUE;
                    }
                    if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableBackgroundV2[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    if(strcmp(SystematicsFocus, "RangeV2") == 0){
                      for(int idx=0; idx<2; idx++){
                          if (variableMinV2[idx] == row[IndexOfVariableMinV2] && variableMaxV2[idx] == row[IndexOfVariableMaxV2]){
                              IndexOfStorage = idx;
                              isVariationFound = kTRUE;
                              continue;
                          }
                      }
                    }
                    
                    for(int ptbin=0; ptbin<NbPtBins; ptbin++){
                        if(isVariationFound){//FIx indices
                            v2pPbSPDC[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDC[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2ATLASSPDC[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASSPDC[ptbin][IndexOfStorage] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ZYAM") == 0){
                            v2pPbSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2pPbSPDC[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 23 + 18*ptbin]));
                            errv2pPbSPDC[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 24 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ATLAS") == 0){
                            v2pPbSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2ATLASSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                        }
                        else if(strcmp(SystematicsFocus, "ExtractionMethod") == 0){//Fixme
                            v2pPbSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 17 + 18*ptbin]));
                            errv2pPbSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 18 + 18*ptbin]));
                            v2pPbSPDC[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 19 + 18*ptbin]));
                            errv2pPbSPDC[ptbin][1] = stod(std::string(row[4 + NbPtBins + 1 + 20 + 18*ptbin]));
                            v2pPbSPDC[ptbin][2] = stod(std::string(row[4 + NbPtBins + 1 + 21 + 18*ptbin]));
                            errv2pPbSPDC[ptbin][2] = stod(std::string(row[4 + NbPtBins + 1 + 22 + 18*ptbin]));
                            v2ATLASSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 29 + 18*ptbin]));
                            errv2ATLASSPDC[ptbin][0] = stod(std::string(row[4 + NbPtBins + 1 + 30 + 18*ptbin]));
                            v2ATLASSPDC[ptbin][1] = stod(std::string(row[5 + NbPtBins + 1 + 16 + 18*NbPtBins + 2*ptbin]));
                            errv2ATLASSPDC[ptbin][1] = stod(std::string(row[6 + NbPtBins + 1 + 16 + 18*NbPtBins + 2*ptbin]));
                        }
                    }
                }
                else{
                    continue;
                }
            }
        
    }
    
    cout << "\n\n\n ==================== \n Systematics results \n ==================== \n\n\n"<<endl;
    
    for(int ptbin=0; ptbin<NbPtBins; ptbin++){
        
        cout << "\n\n === PtBin #"<<ptbin<<" from "<< PtBins[ptbin] << " to " << PtBins[ptbin+1] << " GeV/c === \n\n"<<endl;
        
    if(isClassicVariation){
        if(strcmp(SystematicsFocus, "DPhiCut") == 0){
            default_idx = 3;
        }
        if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
            default_idx = 1;
            NombreEssais = 3;
        }
        if(strcmp(SystematicsFocus, "EtaMin") == 0){
            default_idx = 1;
            NombreEssais = 3;
        }
        if(strcmp(SystematicsFocus, "EtaMax") == 0){
            default_idx = 1;
            NombreEssais = 3;
        }
        if(strcmp(SystematicsFocus, "EMNorm") == 0){
            default_idx = 0;
            NombreEssais = 3;
        }
        if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
            default_idx = 0;
            NombreEssais = 3;
        }
        if(strcmp(SystematicsFocus, "Pooling") == 0){
            default_idx = 0;
            NombreEssais = 2;
        }
        if(strcmp(SystematicsFocus, "InvMass") == 0){
          default_idx = 0;
            NombreEssais = 72;
        }
        if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
            default_idx = 1;
            NombreEssais = 3;
        }
        if(strcmp(SystematicsFocus, "RangeV2") == 0){
            default_idx = 0;
            NombreEssais = 2;
        }
        cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
        //Fixme ptbin boucles
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
            int nbRejected = 0;
            
            
            v2ppSystSPDT[ptbin]=v2pPbSPDT[ptbin][default_idx];//AEUGH
            errv2ppSystSPDT[ptbin]=errv2pPbSPDT[ptbin][default_idx];
            
            for(int idxo = 0; idxo<100; idxo++){
                if(idxo<default_idx){
                    v2ppSystSPDTalt[ptbin][idxo]=v2pPbSPDT[ptbin][idxo];
                    errv2ppSystSPDTalt[ptbin][idxo]=errv2pPbSPDT[ptbin][idxo];
                    
                }
                else if(idxo>default_idx){
                    v2ppSystSPDTalt[ptbin][idxo-1]=v2pPbSPDT[ptbin][idxo];
                    errv2ppSystSPDTalt[ptbin][idxo-1]=errv2pPbSPDT[ptbin][idxo];
                }
            }
                
        
        for(int idx=0; idx<100; idx++){
            if(v2pPbSPDT[ptbin][idx]!=0 || errv2pPbSPDT[ptbin][idx]!=0){
            cout << "v2: "<<v2pPbSPDT[ptbin][idx]<< " +/- "<<errv2pPbSPDT[ptbin][idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbSPDT[ptbin][idx];
                if(v2pPbSPDT[ptbin][idx]>v2_max){
                    v2_max = v2pPbSPDT[ptbin][idx];
                }
                if(v2pPbSPDT[ptbin][idx]<v2_min){
                    v2_min = v2pPbSPDT[ptbin][idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<100; idx++){
            if(v2pPbSPDT[ptbin][idx]!=0 || errv2pPbSPDT[ptbin][idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbSPDT[ptbin][idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbSPDT[ptbin][default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            
        if(strcmp(SystematicsFocus, "InvMass") == 0){
            
            for(int idx=0; idx<100; idx++){
                if(v2pPbSPDT[ptbin][idx]!=0 || errv2pPbSPDT[ptbin][idx]!=0){
                    
                    if(abs(v2pPbSPDT[ptbin][idx]-v2_moy)/v2_syst > 3.){
                        cout << "v2: "<<v2pPbSPDT[ptbin][idx]<< " +/- "<<errv2pPbSPDT[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2pPbSPDT[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                        v2pPbSPDT[ptbin][idx]=0.;
                        errv2pPbSPDT[ptbin][idx]=0.;
                        nbRejected++;
                    }

                }
            }
            
            cout << nbRejected << " cases were rejected because outliers"<<endl;
            
            while(nbRejected > 0){
                nbRejected = 0;
                
                v2_moy = 0;
                v2_syst = 0;
                nbMesures = 0;
                        
                for(int idx=0; idx<100; idx++){
                    if(v2pPbSPDT[ptbin][idx]!=0 || errv2pPbSPDT[ptbin][idx]!=0){
                        nbMesures++;
                        v2_moy+=v2pPbSPDT[ptbin][idx];
                    }
                }
                v2_moy/=nbMesures;
                for(int idx=0; idx<100; idx++){
                    if(v2pPbSPDT[ptbin][idx]!=0 || errv2pPbSPDT[ptbin][idx]!=0){
                        v2_syst+=pow(v2_moy-v2pPbSPDT[ptbin][idx],2);
                    }
                }
                    v2_syst/=nbMesures;
                v2_syst = sqrt(v2_syst);
                    v2_default = v2pPbSPDT[ptbin][default_idx];
                cout << "After removal of outliers - New values are:"<<endl;
                cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                
                for(int idx=0; idx<100; idx++){
                    if(v2pPbSPDT[ptbin][idx]!=0 || errv2pPbSPDT[ptbin][idx]!=0){
                        
                        if(abs(v2pPbSPDT[ptbin][idx]-v2_moy)/v2_syst > 3.){
                            cout << "v2: "<<v2pPbSPDT[ptbin][idx]<< " +/- "<<errv2pPbSPDT[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2pPbSPDT[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                            v2pPbSPDT[ptbin][idx]=0.;
                            errv2pPbSPDT[ptbin][idx]=0.;
                            nbRejected++;
                        }

                    }
                }
                
                cout << nbRejected << " cases were rejected because outliers"<<endl;
            
            }
            
            cout << "Everything has been tidied. The final values after outliers removal are:"<<endl;
            cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            
        }
            
            
            
        //Plot et calculer la systematique
        }
        
        cout << "ATLAS Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
            int nbRejected = 0;
            
            
            v2ppATLASSystSPDT[ptbin]=v2ATLASSPDT[ptbin][default_idx];//AEUGH
            errv2ppATLASSystSPDT[ptbin]=errv2ATLASSPDT[ptbin][default_idx];
            
            for(int idxo = 0; idxo<100; idxo++){
                if(idxo<default_idx){
                    v2ppATLASSystSPDTalt[ptbin][idxo]=v2ATLASSPDT[ptbin][idxo];
                    errv2ppATLASSystSPDTalt[ptbin][idxo]=errv2ATLASSPDT[ptbin][idxo];
                    
                }
                else if(idxo>default_idx){
                    v2ppATLASSystSPDTalt[ptbin][idxo-1]=v2ATLASSPDT[ptbin][idxo];
                    errv2ppATLASSystSPDTalt[ptbin][idxo-1]=errv2ATLASSPDT[ptbin][idxo];
                }
            }
        
        for(int idx=0; idx<100; idx++){
            if(v2ATLASSPDT[ptbin][idx]!=0 || errv2ATLASSPDT[ptbin][idx]!=0){
            cout << "v2: "<<v2ATLASSPDT[ptbin][idx]<< " +/- "<<errv2ATLASSPDT[ptbin][idx] <<endl;
                nbMesures++;
                v2_moy+=v2ATLASSPDT[ptbin][idx];
                if(v2ATLASSPDT[ptbin][idx]>v2_max){
                    v2_max = v2ATLASSPDT[ptbin][idx];
                }
                if(v2ATLASSPDT[ptbin][idx]<v2_min){
                    v2_min = v2ATLASSPDT[ptbin][idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<100; idx++){
            if(v2ATLASSPDT[ptbin][idx]!=0 || errv2ATLASSPDT[ptbin][idx]!=0){
                v2_syst+=pow(v2_moy-v2ATLASSPDT[ptbin][idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2ATLASSPDT[ptbin][default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            
            
            if(strcmp(SystematicsFocus, "InvMass") == 0){
                
                for(int idx=0; idx<100; idx++){
                    if(v2ATLASSPDT[ptbin][idx]!=0 || errv2ATLASSPDT[ptbin][idx]!=0){
                        
                        if(abs(v2ATLASSPDT[ptbin][idx]-v2_moy)/v2_syst > 3.){
                            cout << "v2: "<<v2ATLASSPDT[ptbin][idx]<< " +/- "<<errv2ATLASSPDT[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2ATLASSPDT[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                            v2ATLASSPDT[ptbin][idx]=0.;
                            errv2ATLASSPDT[ptbin][idx]=0.;
                            nbRejected++;
                        }

                    }
                }
                
                cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                while(nbRejected > 0){
                    nbRejected = 0;
                    
                    v2_moy = 0;
                    v2_syst = 0;
                    nbMesures = 0;
                            
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASSPDT[ptbin][idx]!=0 || errv2ATLASSPDT[ptbin][idx]!=0){
                            nbMesures++;
                            v2_moy+=v2ATLASSPDT[ptbin][idx];
                        }
                    }
                    v2_moy/=nbMesures;
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASSPDT[ptbin][idx]!=0 || errv2ATLASSPDT[ptbin][idx]!=0){
                            v2_syst+=pow(v2_moy-v2ATLASSPDT[ptbin][idx],2);
                        }
                    }
                        v2_syst/=nbMesures;
                    v2_syst = sqrt(v2_syst);
                        v2_default = v2ATLASSPDT[ptbin][default_idx];
                    cout << "After removal of outliers - New values are:"<<endl;
                    cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                    
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASSPDT[ptbin][idx]!=0 || errv2ATLASSPDT[ptbin][idx]!=0){
                            
                            if(abs(v2ATLASSPDT[ptbin][idx]-v2_moy)/v2_syst > 3.){
                                cout << "v2: "<<v2ATLASSPDT[ptbin][idx]<< " +/- "<<errv2ATLASSPDT[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2ATLASSPDT[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                                v2ATLASSPDT[ptbin][idx]=0.;
                                errv2ATLASSPDT[ptbin][idx]=0.;
                                nbRejected++;
                            }

                        }
                    }
                    
                    cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                }
                
                cout << "Everything has been tidied. The final values after outliers removal are:"<<endl;
                cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                
            }
            
            
        }
        
        
        cout << "Centrality estimator : SPDClustersPercentile"<<endl;
        
        //Fixme ptbin boucles
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
            int nbRejected = 0;
                
        v2ppSystSPDC[ptbin]=v2pPbSPDC[ptbin][default_idx];//AEUGH
        errv2ppSystSPDC[ptbin]=errv2pPbSPDC[ptbin][default_idx];
        
        for(int idxo = 0; idxo<100; idxo++){
            if(idxo<default_idx){
                v2ppSystSPDCalt[ptbin][idxo]=v2pPbSPDC[ptbin][idxo];
                errv2ppSystSPDCalt[ptbin][idxo]=errv2pPbSPDC[ptbin][idxo];
                
            }
            else if(idxo>default_idx){
                v2ppSystSPDCalt[ptbin][idxo-1]=v2pPbSPDC[ptbin][idxo];
                errv2ppSystSPDCalt[ptbin][idxo-1]=errv2pPbSPDC[ptbin][idxo];
            }
        }
            
        for(int idx=0; idx<100; idx++){
            if(v2pPbSPDC[ptbin][idx]!=0 || errv2pPbSPDC[ptbin][idx]!=0){
            cout << "v2: "<<v2pPbSPDC[ptbin][idx]<< " +/- "<<errv2pPbSPDC[ptbin][idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbSPDC[ptbin][idx];
                if(v2pPbSPDC[ptbin][idx]>v2_max){
                    v2_max = v2pPbSPDC[ptbin][idx];
                }
                if(v2pPbSPDC[ptbin][idx]<v2_min){
                    v2_min = v2pPbSPDC[ptbin][idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<100; idx++){
            if(v2pPbSPDC[ptbin][idx]!=0 || errv2pPbSPDC[ptbin][idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbSPDC[ptbin][idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbSPDC[ptbin][default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            
            if(strcmp(SystematicsFocus, "InvMass") == 0){
                
                for(int idx=0; idx<100; idx++){
                    if(v2pPbSPDC[ptbin][idx]!=0 || errv2pPbSPDC[ptbin][idx]!=0){
                        
                        if(abs(v2pPbSPDC[ptbin][idx]-v2_moy)/v2_syst > 3.){
                            cout << "v2: "<<v2pPbSPDC[ptbin][idx]<< " +/- "<<errv2pPbSPDC[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2pPbSPDC[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                            v2pPbSPDC[ptbin][idx]=0.;
                            errv2pPbSPDC[ptbin][idx]=0.;
                            nbRejected++;
                        }

                    }
                }
                
                cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                while(nbRejected > 0){
                    nbRejected = 0;
                    
                    v2_moy = 0;
                    v2_syst = 0;
                    nbMesures = 0;
                            
                    for(int idx=0; idx<100; idx++){
                        if(v2pPbSPDC[ptbin][idx]!=0 || errv2pPbSPDC[ptbin][idx]!=0){
                            nbMesures++;
                            v2_moy+=v2pPbSPDC[ptbin][idx];
                        }
                    }
                    v2_moy/=nbMesures;
                    for(int idx=0; idx<100; idx++){
                        if(v2pPbSPDC[ptbin][idx]!=0 || errv2pPbSPDC[ptbin][idx]!=0){
                            v2_syst+=pow(v2_moy-v2pPbSPDC[ptbin][idx],2);
                        }
                    }
                        v2_syst/=nbMesures;
                    v2_syst = sqrt(v2_syst);
                        v2_default = v2pPbSPDC[ptbin][default_idx];
                    cout << "After removal of outliers - New values are:"<<endl;
                    cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                    
                    for(int idx=0; idx<100; idx++){
                        if(v2pPbSPDC[ptbin][idx]!=0 || errv2pPbSPDC[ptbin][idx]!=0){
                            
                            if(abs(v2pPbSPDC[ptbin][idx]-v2_moy)/v2_syst > 3.){
                                cout << "v2: "<<v2pPbSPDC[ptbin][idx]<< " +/- "<<errv2pPbSPDC[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2pPbSPDC[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                                v2pPbSPDC[ptbin][idx]=0.;
                                errv2pPbSPDC[ptbin][idx]=0.;
                                nbRejected++;
                            }

                        }
                    }
                    
                    cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                }
                
                cout << "Everything has been tidied. The final values after outliers removal are:"<<endl;
                cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                
            }
            
            
            
        //Plot et calculer la systematique
        }
        
        cout << "ATLAS Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
            int nbRejected = 0;
            
            v2ppATLASSystSPDC[ptbin]=v2ATLASSPDC[ptbin][default_idx];//AEUGH
            errv2ppATLASSystSPDC[ptbin]=errv2ATLASSPDC[ptbin][default_idx];
            
            for(int idxo = 0; idxo<100; idxo++){
                if(idxo<default_idx){
                    v2ppATLASSystSPDCalt[ptbin][idxo]=v2ATLASSPDC[ptbin][idxo];
                    errv2ppATLASSystSPDCalt[ptbin][idxo]=errv2ATLASSPDC[ptbin][idxo];
                    
                }
                else if(idxo>default_idx){
                    v2ppATLASSystSPDCalt[ptbin][idxo-1]=v2ATLASSPDC[ptbin][idxo];
                    errv2ppATLASSystSPDCalt[ptbin][idxo-1]=errv2ATLASSPDC[ptbin][idxo];
                }
            }
        
        for(int idx=0; idx<100; idx++){
            if(v2ATLASSPDC[ptbin][idx]!=0 || errv2ATLASSPDC[ptbin][idx]!=0){
            cout << "v2: "<<v2ATLASSPDC[ptbin][idx]<< " +/- "<<errv2ATLASSPDC[ptbin][idx] <<endl;
                nbMesures++;
                v2_moy+=v2ATLASSPDC[ptbin][idx];
                if(v2ATLASSPDC[ptbin][idx]>v2_max){
                    v2_max = v2ATLASSPDC[ptbin][idx];
                }
                if(v2ATLASSPDC[ptbin][idx]<v2_min){
                    v2_min = v2ATLASSPDC[ptbin][idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<100; idx++){
            if(v2ATLASSPDC[ptbin][idx]!=0 || errv2ATLASSPDC[ptbin][idx]!=0){
                v2_syst+=pow(v2_moy-v2ATLASSPDC[ptbin][idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2ATLASSPDC[ptbin][default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            
            if(strcmp(SystematicsFocus, "InvMass") == 0){
                
                for(int idx=0; idx<100; idx++){
                    if(v2ATLASSPDC[ptbin][idx]!=0 || errv2ATLASSPDC[ptbin][idx]!=0){
                        
                        if(abs(v2ATLASSPDC[ptbin][idx]-v2_moy)/v2_syst > 3.){
                            cout << "v2: "<<v2ATLASSPDC[ptbin][idx]<< " +/- "<<errv2ATLASSPDC[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2ATLASSPDC[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                            v2ATLASSPDC[ptbin][idx]=0.;
                            errv2ATLASSPDC[ptbin][idx]=0.;
                            nbRejected++;
                        }

                    }
                }
                
                cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                while(nbRejected > 0){
                    nbRejected = 0;
                    
                    v2_moy = 0;
                    v2_syst = 0;
                    nbMesures = 0;
                            
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASSPDC[ptbin][idx]!=0 || errv2ATLASSPDC[ptbin][idx]!=0){
                            nbMesures++;
                            v2_moy+=v2ATLASSPDC[ptbin][idx];
                        }
                    }
                    v2_moy/=nbMesures;
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASSPDC[ptbin][idx]!=0 || errv2ATLASSPDC[ptbin][idx]!=0){
                            v2_syst+=pow(v2_moy-v2ATLASSPDC[ptbin][idx],2);
                        }
                    }
                        v2_syst/=nbMesures;
                    v2_syst = sqrt(v2_syst);
                        v2_default = v2ATLASSPDC[ptbin][default_idx];
                    cout << "After removal of outliers - New values are:"<<endl;
                    cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                    
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASSPDC[ptbin][idx]!=0 || errv2ATLASSPDC[ptbin][idx]!=0){
                            
                            if(abs(v2ATLASSPDC[ptbin][idx]-v2_moy)/v2_syst > 3.){
                                cout << "v2: "<<v2ATLASSPDC[ptbin][idx]<< " +/- "<<errv2ATLASSPDC[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2ATLASSPDC[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                                v2ATLASSPDC[ptbin][idx]=0.;
                                errv2ATLASSPDC[ptbin][idx]=0.;
                                nbRejected++;
                            }

                        }
                    }
                    
                    cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                }
                
                cout << "Everything has been tidied. The final values after outliers removal are:"<<endl;
                cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                
            }
            
            
            
        }
        
        
        cout << "Centrality estimator : V0MPercentile"<<endl;
        
        //Fixme ptbin boucles
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
            int nbRejected = 0;
            
            
            v2ppSystV0M[ptbin]=v2pPbV0M[ptbin][default_idx];//AEUGH
            errv2ppSystV0M[ptbin]=errv2pPbV0M[ptbin][default_idx];
            
            for(int idxo = 0; idxo<100; idxo++){
                if(idxo<default_idx){
                    v2ppSystV0Malt[ptbin][idxo]=v2pPbV0M[ptbin][idxo];
                    errv2ppSystV0Malt[ptbin][idxo]=errv2pPbV0M[ptbin][idxo];
                    
                }
                else if(idxo>default_idx){
                    v2ppSystV0Malt[ptbin][idxo-1]=v2pPbV0M[ptbin][idxo];
                    errv2ppSystV0Malt[ptbin][idxo-1]=errv2pPbV0M[ptbin][idxo];
                }
            }
            
            for(int indice=0; indice<NombreEssais; indice++){
                if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
                    naming[indice] = variableZvtxCut[indice];
                }
                if(strcmp(SystematicsFocus, "EtaMin") == 0){
                    naming[indice] = variableEtaMin[indice];
                }
                if(strcmp(SystematicsFocus, "EtaMax") == 0){
                    naming[indice] = variableEtaMax[indice];
                }
                if(strcmp(SystematicsFocus, "EMNorm") == 0){
                    naming[indice] = variableEMNorm[indice];
                }
                if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
                    naming[indice] = variableSummationZvtx[indice];
                }
                if(strcmp(SystematicsFocus, "Pooling") == 0){
                    naming[indice] = variablePooling[indice];
                }
                if(strcmp(SystematicsFocus, "InvMass") == 0){
                    naming[indice] = variableSignal[int((indice-(indice%24))/24.)]+"_"+variableBackground[int(((indice-(indice%6))%24)/6.)]+"_"+variableRatSigma[int(((indice-(indice%3))%6)/3.)]+"_"+variableMinMass[int(indice%3)]+"_"+variableMaxMass[int(indice%3)];
                }
                if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
                    naming[indice] = variableBackgroundV2[indice];
                }
                if(strcmp(SystematicsFocus, "RangeV2") == 0){
                    naming[indice] = variableMinV2[indice]+"_"+variableMaxV2[indice];
                }
            }
                
        
        for(int idx=0; idx<100; idx++){
            if(v2pPbV0M[ptbin][idx]!=0 || errv2pPbV0M[ptbin][idx]!=0){
            cout << "v2: "<<v2pPbV0M[ptbin][idx]<< " +/- "<<errv2pPbV0M[ptbin][idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbV0M[ptbin][idx];
                if(v2pPbV0M[ptbin][idx]>v2_max){
                    v2_max = v2pPbV0M[ptbin][idx];
                }
                if(v2pPbV0M[ptbin][idx]<v2_min){
                    v2_min = v2pPbV0M[ptbin][idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<100; idx++){
            if(v2pPbV0M[ptbin][idx]!=0 || errv2pPbV0M[ptbin][idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbV0M[ptbin][idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbV0M[ptbin][default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            
            if(strcmp(SystematicsFocus, "InvMass") == 0){
                
                for(int idx=0; idx<100; idx++){
                    if(v2pPbV0M[ptbin][idx]!=0 || errv2pPbV0M[ptbin][idx]!=0){
                        
                        if(abs(v2pPbV0M[ptbin][idx]-v2_moy)/v2_syst > 3.){
                            cout << "v2: "<<v2pPbV0M[ptbin][idx]<< " +/- "<<errv2pPbV0M[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2pPbV0M[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                            v2pPbV0M[ptbin][idx]=0.;
                            errv2pPbV0M[ptbin][idx]=0.;
                            nbRejected++;
                        }

                    }
                }
                
                cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                while(nbRejected > 0){
                    nbRejected = 0;
                    
                    v2_moy = 0;
                    v2_syst = 0;
                    nbMesures = 0;
                            
                    for(int idx=0; idx<100; idx++){
                        if(v2pPbV0M[ptbin][idx]!=0 || errv2pPbV0M[ptbin][idx]!=0){
                            nbMesures++;
                            v2_moy+=v2pPbV0M[ptbin][idx];
                        }
                    }
                    v2_moy/=nbMesures;
                    for(int idx=0; idx<100; idx++){
                        if(v2pPbV0M[ptbin][idx]!=0 || errv2pPbV0M[ptbin][idx]!=0){
                            v2_syst+=pow(v2_moy-v2pPbV0M[ptbin][idx],2);
                        }
                    }
                        v2_syst/=nbMesures;
                    v2_syst = sqrt(v2_syst);
                        v2_default = v2pPbV0M[ptbin][default_idx];
                    cout << "After removal of outliers - New values are:"<<endl;
                    cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                    
                    for(int idx=0; idx<100; idx++){
                        if(v2pPbV0M[ptbin][idx]!=0 || errv2pPbV0M[ptbin][idx]!=0){
                            
                            if(abs(v2pPbV0M[ptbin][idx]-v2_moy)/v2_syst > 3.){
                                cout << "v2: "<<v2pPbV0M[ptbin][idx]<< " +/- "<<errv2pPbV0M[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2pPbV0M[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                                v2pPbV0M[ptbin][idx]=0.;
                                errv2pPbV0M[ptbin][idx]=0.;
                                nbRejected++;
                            }

                        }
                    }
                    
                    cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                }
                
                cout << "Everything has been tidied. The final values after outliers removal are:"<<endl;
                cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                
            }
            
            
            
        //Plot et calculer la systematique
        }
        
        cout << "ATLAS Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
            int nbRejected = 0;
            
            v2ppATLASSystV0M[ptbin]=v2ATLASV0M[ptbin][default_idx];//AEUGH
                errv2ppATLASSystV0M[ptbin]=errv2ATLASV0M[ptbin][default_idx];
                
                for(int idxo = 0; idxo<100; idxo++){
                    if(idxo<default_idx){
                        v2ppATLASSystV0Malt[ptbin][idxo]=v2ATLASV0M[ptbin][idxo];
                        errv2ppATLASSystV0Malt[ptbin][idxo]=errv2ATLASV0M[ptbin][idxo];
                        
                    }
                    else if(idxo>default_idx){
                        v2ppATLASSystV0Malt[ptbin][idxo-1]=v2ATLASV0M[ptbin][idxo];
                        errv2ppATLASSystV0Malt[ptbin][idxo-1]=errv2ATLASV0M[ptbin][idxo];
                    }
                }
            
        
        for(int idx=0; idx<100; idx++){
            if(v2ATLASV0M[ptbin][idx]!=0 || errv2ATLASV0M[ptbin][idx]!=0){
            cout << "v2: "<<v2ATLASV0M[ptbin][idx]<< " +/- "<<errv2ATLASV0M[ptbin][idx] <<endl;
                nbMesures++;
                v2_moy+=v2ATLASV0M[ptbin][idx];
                if(v2ATLASV0M[ptbin][idx]>v2_max){
                    v2_max = v2ATLASV0M[ptbin][idx];
                }
                if(v2ATLASV0M[ptbin][idx]<v2_min){
                    v2_min = v2ATLASV0M[ptbin][idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<100; idx++){
            if(v2ATLASV0M[ptbin][idx]!=0 || errv2ATLASV0M[ptbin][idx]!=0){
                v2_syst+=pow(v2_moy-v2ATLASV0M[ptbin][idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2ATLASV0M[ptbin][default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            
            if(strcmp(SystematicsFocus, "InvMass") == 0){
                
                for(int idx=0; idx<100; idx++){
                    if(v2ATLASV0M[ptbin][idx]!=0 || errv2ATLASV0M[ptbin][idx]!=0){
                        
                        if(abs(v2ATLASV0M[ptbin][idx]-v2_moy)/v2_syst > 3.){
                            cout << "v2: "<<v2ATLASV0M[ptbin][idx]<< " +/- "<<errv2ATLASV0M[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2ATLASV0M[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                            v2ATLASV0M[ptbin][idx]=0.;
                            errv2ATLASV0M[ptbin][idx]=0.;
                            nbRejected++;
                        }

                    }
                }
                
                cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                while(nbRejected > 0){
                    nbRejected = 0;
                    
                    v2_moy = 0;
                    v2_syst = 0;
                    nbMesures = 0;
                            
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASV0M[ptbin][idx]!=0 || errv2ATLASV0M[ptbin][idx]!=0){
                            nbMesures++;
                            v2_moy+=v2ATLASV0M[ptbin][idx];
                        }
                    }
                    v2_moy/=nbMesures;
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASV0M[ptbin][idx]!=0 || errv2ATLASV0M[ptbin][idx]!=0){
                            v2_syst+=pow(v2_moy-v2ATLASV0M[ptbin][idx],2);
                        }
                    }
                        v2_syst/=nbMesures;
                    v2_syst = sqrt(v2_syst);
                        v2_default = v2ATLASV0M[ptbin][default_idx];
                    cout << "After removal of outliers - New values are:"<<endl;
                    cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                    
                    for(int idx=0; idx<100; idx++){
                        if(v2ATLASV0M[ptbin][idx]!=0 || errv2ATLASV0M[ptbin][idx]!=0){
                            
                            if(abs(v2ATLASV0M[ptbin][idx]-v2_moy)/v2_syst > 3.){
                                cout << "v2: "<<v2ATLASV0M[ptbin][idx]<< " +/- "<<errv2ATLASV0M[ptbin][idx] << " too far from "<< v2_moy << " by " << abs(v2ATLASV0M[ptbin][idx]-v2_moy)/v2_syst << " sigmas - REJECTED"<< endl;
                                v2ATLASV0M[ptbin][idx]=0.;
                                errv2ATLASV0M[ptbin][idx]=0.;
                                nbRejected++;
                            }

                        }
                    }
                    
                    cout << nbRejected << " cases were rejected because outliers"<<endl;
                
                }
                
                cout << "Everything has been tidied. The final values after outliers removal are:"<<endl;
                cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                
            }
            
            
        }
    }
    else if(strcmp(SystematicsFocus, "ZYAM") == 0){
        cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbSPDT[ptbin][0];
            v2_syst = abs(v2pPbSPDT[ptbin][1]-v2pPbSPDT[ptbin][0]);
            cout << "v2_noZYAM: "<<v2pPbSPDT[ptbin][1]<< " +/- "<<errv2pPbSPDT[ptbin][1] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : SPDClustersPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbSPDC[ptbin][0];
            v2_syst = abs(v2pPbSPDC[ptbin][1]-v2pPbSPDC[ptbin][0]);
            cout << "v2_noZYAM: "<<v2pPbSPDC[ptbin][1]<< " +/- "<<errv2pPbSPDC[ptbin][1] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : V0MPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbV0M[ptbin][0];
            v2_syst = abs(v2pPbV0M[ptbin][1]-v2pPbV0M[ptbin][0]);
            cout << "v2_noZYAM: "<<v2pPbV0M[ptbin][1]<< " +/- "<<errv2pPbV0M[ptbin][1] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
    }
    
    else if(strcmp(SystematicsFocus, "ATLAS") == 0){
        cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbSPDT[ptbin][0];
            v2_syst = abs(v2ATLASSPDT[ptbin][0]-v2pPbSPDT[ptbin][0]);
            cout << "v2_ATLAS: "<<v2ATLASSPDT[ptbin][0]<< " +/- "<<errv2ATLASSPDT[ptbin][0] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : SPDClustersPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbSPDC[ptbin][0];
            v2_syst = abs(v2ATLASSPDC[ptbin][0]-v2pPbSPDC[ptbin][0]);
            cout << "v2_ATLAS: "<<v2ATLASSPDC[ptbin][0]<< " +/- "<<errv2ATLASSPDC[ptbin][0] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : V0MPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbV0M[ptbin][0];
            v2_syst = abs(v2ATLASV0M[ptbin][0]-v2pPbV0M[ptbin][0]);
            cout << "v2_ATLAS: "<<v2ATLASV0M[ptbin][0]<< " +/- "<<errv2ATLASV0M[ptbin][0] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
    }
        
        else if(strcmp(SystematicsFocus, "ExtractionMethod") == 0){//Fixme
            cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
            
            cout << "p-Pb Method"<<endl;
            {
            double v2_default = 0;
            double v2_moy = 0;
            double v2_syst = 0;
                double v2_min = 999.;
                double v2_max = -999.;
            int nbMesures = 0;
                
                for(int idx=0; idx<3; idx++){
                        if(v2pPbSPDT[ptbin][idx]>v2_max){
                            v2_max = v2pPbSPDT[ptbin][idx];
                        }
                        if(v2pPbSPDT[ptbin][idx]<v2_min){
                            v2_min = v2pPbSPDT[ptbin][idx];
                        }
                }
                
                v2_default = v2pPbSPDT[ptbin][0];
                v2_moy = (v2pPbSPDT[ptbin][0]+v2pPbSPDT[ptbin][1]+v2pPbSPDT[ptbin][2])/3.;
                v2_syst = sqrt((pow(v2pPbSPDT[ptbin][0]-v2_moy,2)+pow(v2pPbSPDT[ptbin][1]-v2_moy,2)+pow(v2pPbSPDT[ptbin][2]-v2_moy,2))/3);
                cout << "v2_Ext1: "<<v2pPbSPDT[ptbin][0]<< " +/- "<<errv2pPbSPDT[ptbin][0] <<endl;
                cout << "v2_Ext2: "<<v2pPbSPDT[ptbin][1]<< " +/- "<<errv2pPbSPDT[ptbin][1] <<endl;
                cout << "v2_Ext3: "<<v2pPbSPDT[ptbin][2]<< " +/- "<<errv2pPbSPDT[ptbin][2] <<endl;
            cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            //Plot et calculer la systematique
            }
            
            cout << "ATLAS Method"<<endl;
            {
            double v2_default = 0;
            double v2_moy = 0;
            double v2_syst = 0;
                double v2_min = 999.;
                double v2_max = -999.;
            int nbMesures = 0;
                
                for(int idx=0; idx<2; idx++){
                        if(v2ATLASSPDT[ptbin][idx]>v2_max){
                            v2_max = v2ATLASSPDT[ptbin][idx];
                        }
                        if(v2ATLASSPDT[ptbin][idx]<v2_min){
                            v2_min = v2ATLASSPDT[ptbin][idx];
                        }
                }
                
                v2_default = v2ATLASSPDT[ptbin][0];
                v2_moy = (v2ATLASSPDT[ptbin][0]+v2ATLASSPDT[ptbin][1])/2.;
                v2_syst = sqrt((pow(v2ATLASSPDT[ptbin][0]-v2_moy,2)+pow(v2ATLASSPDT[ptbin][1]-v2_moy,2))/2);
                cout << "v2_Ext1: "<<v2ATLASSPDT[ptbin][0]<< " +/- "<<errv2ATLASSPDT[ptbin][0] <<endl;
                cout << "v2_Ext2: "<<v2ATLASSPDT[ptbin][1]<< " +/- "<<errv2ATLASSPDT[ptbin][1] <<endl;
            cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            //Plot et calculer la systematique
            }
            
            cout << "Centrality estimator : SPDClustersPercentile"<<endl;
            
            cout << "p-Pb Method"<<endl;
            {
            double v2_default = 0;
                double v2_moy = 0;
            double v2_syst = 0;
                double v2_min = 999.;
                double v2_max = -999.;
            int nbMesures = 0;
                
                
                for(int idx=0; idx<3; idx++){
                        if(v2pPbSPDC[ptbin][idx]>v2_max){
                            v2_max = v2pPbSPDC[ptbin][idx];
                        }
                        if(v2pPbSPDC[ptbin][idx]<v2_min){
                            v2_min = v2pPbSPDC[ptbin][idx];
                        }
                }
                
                
                v2_default = v2pPbSPDC[ptbin][0];
                v2_moy = (v2pPbSPDC[ptbin][0]+v2pPbSPDC[ptbin][1]+v2pPbSPDC[ptbin][2])/3.;
                v2_syst = sqrt((pow(v2pPbSPDC[ptbin][0]-v2_moy,2)+pow(v2pPbSPDC[ptbin][1]-v2_moy,2)+pow(v2pPbSPDC[ptbin][2]-v2_moy,2))/3);
                cout << "v2_Ext1: "<<v2pPbSPDC[ptbin][0]<< " +/- "<<errv2pPbSPDC[ptbin][0] <<endl;
                cout << "v2_Ext2: "<<v2pPbSPDC[ptbin][1]<< " +/- "<<errv2pPbSPDC[ptbin][1] <<endl;
                cout << "v2_Ext3: "<<v2pPbSPDC[ptbin][2]<< " +/- "<<errv2pPbSPDC[ptbin][2] <<endl;
            cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            //Plot et calculer la systematique
            }
            
            cout << "ATLAS Method"<<endl;
            {
            double v2_default = 0;
                double v2_moy = 0;
            double v2_syst = 0;
                double v2_min = 999.;
                double v2_max = -999.;
            int nbMesures = 0;
                
                
                for(int idx=0; idx<2; idx++){
                        if(v2ATLASSPDC[ptbin][idx]>v2_max){
                            v2_max = v2ATLASSPDC[ptbin][idx];
                        }
                        if(v2ATLASSPDC[ptbin][idx]<v2_min){
                            v2_min = v2ATLASSPDC[ptbin][idx];
                        }
                }
                
                
                v2_default = v2ATLASSPDC[ptbin][0];
                v2_moy = (v2ATLASSPDC[ptbin][0]+v2ATLASSPDC[ptbin][1])/2.;
                v2_syst = sqrt((pow(v2ATLASSPDC[ptbin][0]-v2_moy,2)+pow(v2ATLASSPDC[ptbin][1]-v2_moy,2))/2);
                cout << "v2_Ext1: "<<v2ATLASSPDC[ptbin][0]<< " +/- "<<errv2ATLASSPDC[ptbin][0] <<endl;
                cout << "v2_Ext2: "<<v2ATLASSPDC[ptbin][1]<< " +/- "<<errv2ATLASSPDC[ptbin][1] <<endl;
            cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            //Plot et calculer la systematique
            }
            
            cout << "Centrality estimator : V0MPercentile"<<endl;
            
            cout << "p-Pb Method"<<endl;
            {
            double v2_default = 0;
                double v2_moy = 0;
            double v2_syst = 0;
                double v2_min = 999.;
                double v2_max = -999.;
            int nbMesures = 0;
                
                for(int idx=0; idx<3; idx++){
                        if(v2pPbV0M[ptbin][idx]>v2_max){
                            v2_max = v2pPbV0M[ptbin][idx];
                        }
                        if(v2pPbV0M[ptbin][idx]<v2_min){
                            v2_min = v2pPbV0M[ptbin][idx];
                        }
                }
                
                
                v2_default = v2pPbV0M[ptbin][0];
                v2_moy = (v2pPbV0M[ptbin][0]+v2pPbV0M[ptbin][1]+v2pPbV0M[ptbin][2])/3.;
                v2_syst = sqrt((pow(v2pPbV0M[ptbin][0]-v2_moy,2)+pow(v2pPbV0M[ptbin][1]-v2_moy,2)+pow(v2pPbV0M[ptbin][2]-v2_moy,2))/3);
                cout << "v2_Ext1: "<<v2pPbV0M[ptbin][0]<< " +/- "<<errv2pPbV0M[ptbin][0] <<endl;
                cout << "v2_Ext2: "<<v2pPbV0M[ptbin][1]<< " +/- "<<errv2pPbV0M[ptbin][1] <<endl;
                cout << "v2_Ext3: "<<v2pPbV0M[ptbin][2]<< " +/- "<<errv2pPbV0M[ptbin][2] <<endl;
            cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            //Plot et calculer la systematique
            }
            
            cout << "ATLAS Method"<<endl;
            {
            double v2_default = 0;
                double v2_moy = 0;
            double v2_syst = 0;
                double v2_min = 999.;
                double v2_max = -999.;
            int nbMesures = 0;
                
                for(int idx=0; idx<2; idx++){
                        if(v2ATLASV0M[ptbin][idx]>v2_max){
                            v2_max = v2ATLASV0M[ptbin][idx];
                        }
                        if(v2ATLASV0M[ptbin][idx]<v2_min){
                            v2_min = v2ATLASV0M[ptbin][idx];
                        }
                }
                
                
                v2_default = v2ATLASV0M[ptbin][0];
                v2_moy = (v2ATLASV0M[ptbin][0]+v2ATLASV0M[ptbin][1])/2.;
                v2_syst = sqrt((pow(v2ATLASV0M[ptbin][0]-v2_moy,2)+pow(v2ATLASV0M[ptbin][1]-v2_moy,2))/2);
                cout << "v2_Ext1: "<<v2ATLASV0M[ptbin][0]<< " +/- "<<errv2ATLASV0M[ptbin][0] <<endl;
                cout << "v2_Ext2: "<<v2ATLASV0M[ptbin][1]<< " +/- "<<errv2ATLASV0M[ptbin][1] <<endl;
            cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
                cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
            //Plot et calculer la systematique
            }
        }
   
}
    
    {
    
    Double_t sit[NbPtBins];//AEUGH
    Double_t errsit[NbPtBins];
    
    int counter = 0;

    for (Int_t i=0;i<NbPtBins;i++){
        sit[i]=(PtBins[i+1]+PtBins[i])/2.;
    errsit[i]=0;
    };
    
    TCanvas *clel = new TCanvas("c1", "c1",0,0,2000,1000);
    clel->cd();
    
    TGraphErrors *v2systdimu = new TGraphErrors(NbPtBins,sit,v2ppSystV0M,errsit,errv2ppSystV0M);
    Char_t TitlePlot[500];
    sprintf(TitlePlot,"Sytematics study - YieldSub - V0M");
     v2systdimu->SetTitle(TitlePlot);
     v2systdimu->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
     v2systdimu->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimu->GetHistogram()->GetXaxis()->SetTitle("p_T");
    v2systdimu->SetMarkerStyle(20);
    v2systdimu->SetMarkerColor(kAzure-3);
    v2systdimu->SetLineColor(kAzure-3);
    
    v2systdimu->Draw("AP");
    
    TLegend *legendov=new TLegend(0.12,0.60,0.40,0.85);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.035);
    
    Char_t message[80];
    sprintf(message,"%s",naming[default_idx].c_str());
    
    legendov->AddEntry(v2systdimu,message);

                     
    
    
    for(int essai = 0; essai<NombreEssais-1;essai++){
        
        bool isValid = kFALSE;
        for(int ptbin=0; ptbin<NbPtBins; ptbin++){
        v2ppSystV0Malto[ptbin] = v2ppSystV0Malt[ptbin][essai];
            if(v2ppSystV0Malto[ptbin]!=0){
                isValid = kTRUE;
            }
        errv2ppSystV0Malto[ptbin] = errv2ppSystV0Malt[ptbin][essai];
            double variation;
            
            if(essai<NombreEssais/2){
                variation = -0.2 + 0.2*essai/NombreEssais;
            }
            else{
                variation = 0.2 - 0.2*(NombreEssais-essai-1)/NombreEssais;
            }
            sit[ptbin] = (PtBins[ptbin+1]+PtBins[ptbin])/2. + variation;
        }
        if(!isValid){
            continue;
        }
        
        
        
        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
        TGraphErrors *v2systdimuessai = new TGraphErrors(NbPtBins,sit,v2ppSystV0Malto,errsit,errv2ppSystV0Malto);
        Char_t TitlePlot[500];
        sprintf(TitlePlot,"Sstematics study - YieldSub - V0M");
         v2systdimuessai->SetTitle(TitlePlot);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimuessai->GetHistogram()->GetXaxis()->SetTitle("p_T");
        v2systdimuessai->SetMarkerStyle(33);
        v2systdimuessai->SetLineColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
        v2systdimuessai->SetMarkerColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
       // v2systdimuessai->SetMarkerColor(kBlack);
        //Markerstyle et marker color
        v2systdimuessai->Draw("P same");
        
        Char_t messago[80];
        
        if(essai<default_idx){
            sprintf(messago,"%s",naming[essai].c_str());
        legendov->AddEntry(v2systdimuessai,messago);
        }
        else if (essai>=default_idx){
            sprintf(messago,"%s",naming[essai+1].c_str());
            legendov->AddEntry(v2systdimuessai,messago);
        }
    }
    //legendov->Draw();
        
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop/SystematicsDimu");
        sprintf(CanvasName,"%s/SystDimu_YieldSub_V0M.pdf",FolderName);
        clel->SaveAs(CanvasName);
        
    }
    
    {
    
    Double_t sit[NbPtBins];//AEUGH
    Double_t errsit[NbPtBins];
    
    int counter = 0;

    for (Int_t i=0;i<NbPtBins;i++){
        sit[i]=(PtBins[i+1]+PtBins[i])/2.;
    errsit[i]=0;
    };
    
    TCanvas *clel = new TCanvas("c1", "c1",0,0,2000,1000);
    clel->cd();
    
    TGraphErrors *v2systdimu = new TGraphErrors(NbPtBins,sit,v2ppSystSPDT,errsit,errv2ppSystSPDT);
    Char_t TitlePlot[500];
    sprintf(TitlePlot,"Sytematics study - YieldSub - SPDT");
     v2systdimu->SetTitle(TitlePlot);
     v2systdimu->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
     v2systdimu->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimu->GetHistogram()->GetXaxis()->SetTitle("p_T");
    v2systdimu->SetMarkerStyle(20);
    v2systdimu->SetMarkerColor(kAzure-3);
    v2systdimu->SetLineColor(kAzure-3);
    
    v2systdimu->Draw("AP");
    
   TLegend *legendov=new TLegend(0.12,0.60,0.40,0.85);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.035);
    
    Char_t message[80];
    sprintf(message,"%s",naming[default_idx].c_str());
    
    legendov->AddEntry(v2systdimu,message);

                     
    
    
    for(int essai = 0; essai<NombreEssais-1;essai++){
        bool isValid = kFALSE;
        
        for(int ptbin=0; ptbin<NbPtBins; ptbin++){
        v2ppSystSPDTalto[ptbin] = v2ppSystSPDTalt[ptbin][essai];
            if(v2ppSystSPDTalto[ptbin]!=0){
                isValid=kTRUE;
            }
        errv2ppSystSPDTalto[ptbin] = errv2ppSystSPDTalt[ptbin][essai];
            double variation;
            
            if(essai<NombreEssais/2){
                variation = -0.2 + 0.2*essai/NombreEssais;
            }
            else{
                variation = 0.2 - 0.2*(NombreEssais-essai-1)/NombreEssais;
            }
            sit[ptbin] = (PtBins[ptbin+1]+PtBins[ptbin])/2. + variation;
        }
        
        if(!isValid){
            continue;
        }
        
        
        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
        TGraphErrors *v2systdimuessai = new TGraphErrors(NbPtBins,sit,v2ppSystSPDTalto,errsit,errv2ppSystSPDTalto);
        Char_t TitlePlot[500];
        sprintf(TitlePlot,"Sstematics study - YieldSub - SPDT");
         v2systdimuessai->SetTitle(TitlePlot);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimuessai->GetHistogram()->GetXaxis()->SetTitle("p_T");
        v2systdimuessai->SetMarkerStyle(33);
        v2systdimuessai->SetLineColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
        v2systdimuessai->SetMarkerColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
       // v2systdimuessai->SetMarkerColor(kBlack);
        //Markerstyle et marker color
        v2systdimuessai->Draw("P same");
        
        Char_t messago[80];
        
        if(essai<default_idx){
            sprintf(messago,"%s",naming[essai].c_str());
        legendov->AddEntry(v2systdimuessai,messago);
        }
        else if (essai>=default_idx){
            sprintf(messago,"%s",naming[essai+1].c_str());
            legendov->AddEntry(v2systdimuessai,messago);
        }
    }
    //legendov->Draw();
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop/SystematicsDimu");
        sprintf(CanvasName,"%s/SystDimu_YieldSub_SPDT.pdf",FolderName);
        clel->SaveAs(CanvasName);
    }
    
    
    {
    
    Double_t sit[NbPtBins];//AEUGH
    Double_t errsit[NbPtBins];
    
    int counter = 0;

    for (Int_t i=0;i<NbPtBins;i++){
        sit[i]=(PtBins[i+1]+PtBins[i])/2.;
    errsit[i]=0;
    };
    
    TCanvas *clel = new TCanvas("c1", "c1",0,0,2000,1000);
    clel->cd();
    
    TGraphErrors *v2systdimu = new TGraphErrors(NbPtBins,sit,v2ppSystSPDC,errsit,errv2ppSystSPDC);
    Char_t TitlePlot[500];
    sprintf(TitlePlot,"Sytematics study - YieldSub - SPDC");
     v2systdimu->SetTitle(TitlePlot);
     v2systdimu->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
     v2systdimu->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimu->GetHistogram()->GetXaxis()->SetTitle("p_T");
    v2systdimu->SetMarkerStyle(20);
    v2systdimu->SetMarkerColor(kAzure-3);
    v2systdimu->SetLineColor(kAzure-3);
    
    v2systdimu->Draw("AP");
    
    TLegend *legendov=new TLegend(0.12,0.60,0.40,0.85);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.035);
    
    Char_t message[80];
    sprintf(message,"%s",naming[default_idx].c_str());
    
    legendov->AddEntry(v2systdimu,message);

                     
    
    
    for(int essai = 0; essai<NombreEssais-1;essai++){
        
        bool isValid = kFALSE;
        for(int ptbin=0; ptbin<NbPtBins; ptbin++){
        v2ppSystSPDCalto[ptbin] = v2ppSystSPDCalt[ptbin][essai];
            if(v2ppSystSPDCalto[ptbin]!=0){
                isValid = kTRUE;
            }
        errv2ppSystSPDCalto[ptbin] = errv2ppSystSPDCalt[ptbin][essai];
            double variation;
            
            if(essai<NombreEssais/2){
                variation = -0.2 + 0.2*essai/NombreEssais;
            }
            else{
                variation = 0.2 - 0.2*(NombreEssais-essai-1)/NombreEssais;
            }
            sit[ptbin] = (PtBins[ptbin+1]+PtBins[ptbin])/2. + variation;
        }
        if(!isValid){
            continue;
        }
        
        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
        TGraphErrors *v2systdimuessai = new TGraphErrors(NbPtBins,sit,v2ppSystSPDCalto,errsit,errv2ppSystSPDCalto);
        Char_t TitlePlot[500];
        sprintf(TitlePlot,"Sstematics study - YieldSub - SPDC");
         v2systdimuessai->SetTitle(TitlePlot);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimuessai->GetHistogram()->GetXaxis()->SetTitle("p_T");
        v2systdimuessai->SetMarkerStyle(33);
        v2systdimuessai->SetLineColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
        v2systdimuessai->SetMarkerColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
       // v2systdimuessai->SetMarkerColor(kBlack);
        //Markerstyle et marker color
        v2systdimuessai->Draw("P same");
        
        Char_t messago[80];
        
        if(essai<default_idx){
            sprintf(messago,"%s",naming[essai].c_str());
        legendov->AddEntry(v2systdimuessai,messago);
        }
        else if (essai>=default_idx){
            sprintf(messago,"%s",naming[essai+1].c_str());
            legendov->AddEntry(v2systdimuessai,messago);
        }
    }
    //legendov->Draw();
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop/SystematicsDimu");
        sprintf(CanvasName,"%s/SystDimu_YieldSub_SPDC.pdf",FolderName);
        clel->SaveAs(CanvasName);
        
    }
    
    
    //ATLAS AEUGH
    
    {
    
    Double_t sit[NbPtBins];//AEUGH
    Double_t errsit[NbPtBins];
    
    int counter = 0;

    for (Int_t i=0;i<NbPtBins;i++){
        sit[i]=(PtBins[i+1]+PtBins[i])/2.;
    errsit[i]=0;
    };
    
    TCanvas *clel = new TCanvas("c1", "c1",0,0,2000,1000);
    clel->cd();
    
    TGraphErrors *v2systdimu = new TGraphErrors(NbPtBins,sit,v2ppATLASSystV0M,errsit,errv2ppATLASSystV0M);
    Char_t TitlePlot[500];
    sprintf(TitlePlot,"Sytematics study - Template fit - V0M");
     v2systdimu->SetTitle(TitlePlot);
     v2systdimu->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
     v2systdimu->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimu->GetHistogram()->GetXaxis()->SetTitle("p_T");
    v2systdimu->SetMarkerStyle(21);
    v2systdimu->SetMarkerColor(kRed+1);
    v2systdimu->SetLineColor(kRed+1);
    
    v2systdimu->Draw("AP");
    
    TLegend *legendov=new TLegend(0.12,0.60,0.40,0.85);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.035);
    
    Char_t message[80];
    sprintf(message,"%s",naming[default_idx].c_str());
    
    legendov->AddEntry(v2systdimu,message);

                     
    
    
    for(int essai = 0; essai<NombreEssais-1;essai++){
        
        bool isValid=kFALSE;
        for(int ptbin=0; ptbin<NbPtBins; ptbin++){
        v2ppATLASSystV0Malto[ptbin] = v2ppATLASSystV0Malt[ptbin][essai];
            if(v2ppATLASSystV0Malto[ptbin]!=0){
                isValid = kTRUE;
            }
        errv2ppATLASSystV0Malto[ptbin] = errv2ppATLASSystV0Malt[ptbin][essai];
            double variation;
            
            if(essai<NombreEssais/2){
                variation = -0.2 + 0.2*essai/NombreEssais;
            }
            else{
                variation = 0.2 - 0.2*(NombreEssais-essai-1)/NombreEssais;
            }
            sit[ptbin] = (PtBins[ptbin+1]+PtBins[ptbin])/2. + variation;
        }
        
        if(!isValid){
            continue;
        }
        
        
        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
        TGraphErrors *v2systdimuessai = new TGraphErrors(NbPtBins,sit,v2ppATLASSystV0Malto,errsit,errv2ppATLASSystV0Malto);
        Char_t TitlePlot[500];
        sprintf(TitlePlot,"Sytematics study - Template Fit - V0M");
         v2systdimuessai->SetTitle(TitlePlot);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimuessai->GetHistogram()->GetXaxis()->SetTitle("p_T");
        v2systdimuessai->SetMarkerStyle(33);
        v2systdimuessai->SetLineColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
        v2systdimuessai->SetMarkerColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
       // v2systdimuessai->SetMarkerColor(kBlack);
        //Markerstyle et marker color
        v2systdimuessai->Draw("P same");
        
        Char_t messago[80];
        
        if(essai<default_idx){
            sprintf(messago,"%s",naming[essai].c_str());
        legendov->AddEntry(v2systdimuessai,messago);
        }
        else if (essai>=default_idx){
            sprintf(messago,"%s",naming[essai+1].c_str());
            legendov->AddEntry(v2systdimuessai,messago);
        }
    }
    //legendov->Draw();
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop/SystematicsDimu");
        sprintf(CanvasName,"%s/SystDimu_ATLAS_V0M.pdf",FolderName);
        clel->SaveAs(CanvasName);
    }
    
    {
    
    Double_t sit[NbPtBins];//AEUGH
    Double_t errsit[NbPtBins];
    
    int counter = 0;

    for (Int_t i=0;i<NbPtBins;i++){
        sit[i]=(PtBins[i+1]+PtBins[i])/2.;
    errsit[i]=0;
    };
    
    TCanvas *clel = new TCanvas("c1", "c1",0,0,2000,1000);
    clel->cd();
    
    TGraphErrors *v2systdimu = new TGraphErrors(NbPtBins,sit,v2ppATLASSystSPDT,errsit,errv2ppATLASSystSPDT);
    Char_t TitlePlot[500];
    sprintf(TitlePlot,"Sytematics study - Template fit - SPDT");
     v2systdimu->SetTitle(TitlePlot);
     v2systdimu->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
     v2systdimu->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimu->GetHistogram()->GetXaxis()->SetTitle("p_T");
    v2systdimu->SetMarkerStyle(21);
    v2systdimu->SetMarkerColor(kRed+1);
    v2systdimu->SetLineColor(kRed+1);
    
    v2systdimu->Draw("AP");
    
    TLegend *legendov=new TLegend(0.12,0.55,0.40,0.88);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.005);
    
    Char_t message[80];
    sprintf(message,"%s",naming[default_idx].c_str());
    
    legendov->AddEntry(v2systdimu,message);

                     
    
    
    for(int essai = 0; essai<NombreEssais-1;essai++){
        
        bool isValid=kFALSE;
        for(int ptbin=0; ptbin<NbPtBins; ptbin++){
        v2ppATLASSystSPDTalto[ptbin] = v2ppATLASSystSPDTalt[ptbin][essai];
            if(v2ppATLASSystSPDTalto[ptbin]!=0){
                isValid=kTRUE;
            }
        errv2ppATLASSystSPDTalto[ptbin] = errv2ppATLASSystSPDTalt[ptbin][essai];
            double variation;
            
            if(essai<NombreEssais/2){
                variation = -0.2 + 0.2*essai/NombreEssais;
            }
            else{
                variation = 0.2 - 0.2*(NombreEssais-essai-1)/NombreEssais;
            }
            sit[ptbin] = (PtBins[ptbin+1]+PtBins[ptbin])/2. + variation;
        }
        
        if(!isValid){
            continue;
        }
        
        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
        TGraphErrors *v2systdimuessai = new TGraphErrors(NbPtBins,sit,v2ppATLASSystSPDTalto,errsit,errv2ppATLASSystSPDTalto);
        Char_t TitlePlot[500];
        sprintf(TitlePlot,"Sstematics study - Template fit - SPDT");
         v2systdimuessai->SetTitle(TitlePlot);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimuessai->GetHistogram()->GetXaxis()->SetTitle("p_T");
        v2systdimuessai->SetMarkerStyle(33);
        v2systdimuessai->SetLineColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
        v2systdimuessai->SetMarkerColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
       // v2systdimuessai->SetMarkerColor(kBlack);
        //Markerstyle et marker color
        v2systdimuessai->Draw("P same");
        
        Char_t messago[80];
        
        if(essai<default_idx){
            sprintf(messago,"%s",naming[essai].c_str());
        legendov->AddEntry(v2systdimuessai,messago);
        }
        else if (essai>=default_idx){
            sprintf(messago,"%s",naming[essai+1].c_str());
            legendov->AddEntry(v2systdimuessai,messago);
        }
    }
    //legendov->Draw();
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop/SystematicsDimu");
        sprintf(CanvasName,"%s/SystDimu_ATLAS_SPDT.pdf",FolderName);
        clel->SaveAs(CanvasName);
    }
    
    
    {
    
    Double_t sit[NbPtBins];//AEUGH
    Double_t errsit[NbPtBins];
    
    int counter = 0;

    for (Int_t i=0;i<NbPtBins;i++){
        sit[i]=(PtBins[i+1]+PtBins[i])/2.;
    errsit[i]=0;
    };
    
    TCanvas *clel = new TCanvas("c1", "c1",0,0,2000,1000);
    clel->cd();
    
    TGraphErrors *v2systdimu = new TGraphErrors(NbPtBins,sit,v2ppATLASSystSPDC,errsit,errv2ppATLASSystSPDC);
    Char_t TitlePlot[500];
    sprintf(TitlePlot,"Sytematics study - Template fit - SPDC");
     v2systdimu->SetTitle(TitlePlot);
     v2systdimu->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
     v2systdimu->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimu->GetHistogram()->GetXaxis()->SetTitle("p_T");
    v2systdimu->SetMarkerStyle(21);
    v2systdimu->SetMarkerColor(kRed+1);
    v2systdimu->SetLineColor(kRed+1);
    
    v2systdimu->Draw("AP");
    
    TLegend *legendov=new TLegend(0.12,0.60,0.40,0.85);
                     legendov->SetFillColorAlpha(kWhite, 0.);
                     legendov->SetBorderSize(0);
                      legendov->SetTextFont(42);
                      legendov->SetTextSize(0.035);
    
    Char_t message[80];
    sprintf(message,"%s",naming[default_idx].c_str());
    
    legendov->AddEntry(v2systdimu,message);

                     
    
    
    for(int essai = 0; essai<NombreEssais-1;essai++){
        
        bool isValid=kFALSE;
        for(int ptbin=0; ptbin<NbPtBins; ptbin++){
        v2ppATLASSystSPDCalto[ptbin] = v2ppATLASSystSPDCalt[ptbin][essai];
            if( v2ppATLASSystSPDCalto[ptbin]!=0){
                isValid=kTRUE;
            }
        errv2ppATLASSystSPDCalto[ptbin] = errv2ppATLASSystSPDCalt[ptbin][essai];
            double variation;
            
            if(essai<NombreEssais/2){
                variation = -0.2 + 0.2*essai/NombreEssais;
            }
            else{
                variation = 0.2 - 0.2*(NombreEssais-essai-1)/NombreEssais;
            }
            sit[ptbin] = (PtBins[ptbin+1]+PtBins[ptbin])/2. + variation;
        }
        
        if(!isValid){
            continue;
        }
        
        
        gStyle->SetPalette(kRainBow);
           int nColors = gStyle->GetNumberOfColors();
        
        TGraphErrors *v2systdimuessai = new TGraphErrors(NbPtBins,sit,v2ppATLASSystSPDCalto,errsit,errv2ppATLASSystSPDCalto);
        Char_t TitlePlot[500];
        sprintf(TitlePlot,"Sstematics study - Template fit - SPDC");
         v2systdimuessai->SetTitle(TitlePlot);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
         v2systdimuessai->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2systdimuessai->GetHistogram()->GetXaxis()->SetTitle("p_T");
        v2systdimuessai->SetMarkerStyle(33);
        v2systdimuessai->SetLineColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
        v2systdimuessai->SetMarkerColor(gStyle->GetColorPalette(int(nColors*essai/NombreEssais)));
       // v2systdimuessai->SetMarkerColor(kBlack);
        //Markerstyle et marker color
        v2systdimuessai->Draw("P same");
        
        Char_t messago[80];
        
        if(essai<default_idx){
            sprintf(messago,"%s",naming[essai].c_str());
        legendov->AddEntry(v2systdimuessai,messago);
        }
        else if (essai>=default_idx){
            sprintf(messago,"%s",naming[essai+1].c_str());
            legendov->AddEntry(v2systdimuessai,messago);
        }
    }
    //legendov->Draw();
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop/SystematicsDimu");
        sprintf(CanvasName,"%s/SystDimu_ATLAS_SPDC.pdf",FolderName);
        clel->SaveAs(CanvasName);
    }
    
}



bool isRowSuitableForSystematics(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]){
    //Dire si la row est utile pour la systematique que l'on veut
    
    bool isSuitable = kTRUE;
    //Default values
    
    Char_t DPhiCut[10] = "10mrad"; //0
    Char_t CentralClass[10] = "0-5";  //2
    Char_t PeriphClass[10] = "40-100"; //3
    //NbPtBins 4
    //PtLimits [5,NbPtBins+5]
    Char_t EtaMin[10] = "1.5"; //NbPtBins+5 +1
    Char_t EtaMax[10] = "5"; //NbPtBins+5 +2
    Char_t ZvtxCut[10] = "10"; //NbPtBins+5 +3
    Char_t EMNorm[10] = "Method1"; //NbPtBins+5 +4
    Char_t EMMax[10] = "100"; //NbPtBins+5 +5
    Char_t EMThreshold[10] = "10"; //NbPtBins+5 +6
    Char_t EMChange[10] = "None"; //NbPtBins+5 +7
    Char_t SummationMethod[10] = "Method1c"; //NbPtBins+5 +8
    Char_t BackgroundV2[10] = "Pol2"; //NbPtBins+5 +9
    Char_t V2Min[10] = "1"; //NbPtBins+5 +10
    Char_t V2Max[10] = "5"; //NbPtBins+5 +11
    Char_t MinMass[10] = "2"; //NbPtBins+5 +12
    Char_t MaxMass[10] = "5"; //NbPtBins+5 +13
    Char_t SignalF[10] = "CB2-Run2"; //NbPtBins+5 +14
    Char_t BackgroundF[20] = "DoubleExpo"; //NbPtBins+5 +15
    Char_t RatSigma[10] = "1.05"; //NbPtBins+5 +16
    
    Char_t *DefaultConditions[21];
    memset(DefaultConditions,0,21*sizeof(Char_t*));
    
    Char_t charray[10];
    
    DefaultConditions[0] = DPhiCut;
    DefaultConditions[1] = CentralityEstimator;
    DefaultConditions[2] = CentralClass;
    DefaultConditions[3] = PeriphClass;
  //  cout << NbPtBins << " is NbPtBins"<<endl;
    sprintf(charray, "%i", int(NbPtBins));
    DefaultConditions[4] = charray;
    DefaultConditions[5] = EtaMin;
    DefaultConditions[6] = EtaMax;
    DefaultConditions[7] = ZvtxCut;
    DefaultConditions[8] = EMNorm;
    DefaultConditions[9] = EMMax;
    DefaultConditions[10] = EMThreshold;
    DefaultConditions[11] = EMChange;
    DefaultConditions[12] = SummationMethod;
    DefaultConditions[13] = BackgroundV2;
    DefaultConditions[14] = V2Min;
    DefaultConditions[15] = V2Max;
    DefaultConditions[16] = MinMass;
    DefaultConditions[17] = MaxMass;
    DefaultConditions[18] = SignalF;
    DefaultConditions[19] = BackgroundF;
    DefaultConditions[20] = RatSigma;
    
    //Char_t *DefaultConditions[11] = {DPhiCut,CentralityEstimator,CentralClass,PeriphClass,NbPtbins,        EtaMin,EtaMax,ZvtxCut,EMNorm,EMMax,EMThreshold,EMChange,SummationMethod,BackgroundV2,V2Min,V2Max,MinMass,MaxMass,SignalF,BackgroundF,RatSigma};
    
    for(int index=0;index<NbPtBins+1;index++){
        if(stoi(std::string(row[5+index])) != PtBins[index]){
        //    cout << "row[" << 5+index << "] = " << row[5+index] << " is not the same as " << PtBins[index]<< " ABORT" <<endl;
            isSuitable = kFALSE;
            return isSuitable;
        }
    }
    
    for(int index=0; index<21; index++){
        if (strcmp(SystematicsFocus, "DPhiCut") == 0){
            if(index == 0){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
            if(index == 7){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EtaMin") == 0){
            if(index == 5){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EtaMax") == 0){
            if(index == 6){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EMNorm") == 0){
            if(index == 8){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
            if(index == 12){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "Pooling") == 0){
            if(index == 11){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "InvMass") == 0){
            if(index >= 16 && index <= 20){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
            if(index == 13){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "RangeV2") == 0){
            if(index >= 14 && index <= 15){
                continue;
            }
        }
        if(index<=4){
            if(row[index] != DefaultConditions[index]){
               // cout << index<< "row[" << index << "] = " << row[index] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
                isSuitable = kFALSE;
                return isSuitable;
            }
        }
        else{
            if(row[index+NbPtBins+1] != DefaultConditions[index]){
               // cout << index<< "row[" << index+NbPtBins+1 << "] = " << row[index+NbPtBins+1] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
                isSuitable = kFALSE;
                return isSuitable;
            }
        }
    }
    
    return isSuitable;
    
}


bool isRowSuitableForCentrality(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]){
    //Dire si la row est utile pour la systematique que l'on veut
    
    bool isSuitable = kTRUE;
    //Default values
    
    Char_t DPhiCut[10] = "10mrad"; //0
    Char_t CentralClass[10] = "0-5";  //2
    Char_t PeriphClass[10] = "40-100"; //3
    //NbPtBins 4
    //PtLimits [5,NbPtBins+5]
    Char_t EtaMin[10] = "1.5"; //NbPtBins+5 +1
    Char_t EtaMax[10] = "5"; //NbPtBins+5 +2
    Char_t ZvtxCut[10] = "10"; //NbPtBins+5 +3
    Char_t EMNorm[10] = "Method1"; //NbPtBins+5 +4
    Char_t EMMax[10] = "100"; //NbPtBins+5 +5
    Char_t EMThreshold[10] = "10"; //NbPtBins+5 +6
    Char_t EMChange[10] = "None"; //NbPtBins+5 +7
    Char_t SummationMethod[10] = "Method1c"; //NbPtBins+5 +8
    Char_t BackgroundV2[10] = "Pol2"; //NbPtBins+5 +9
    Char_t V2Min[10] = "1"; //NbPtBins+5 +10
    Char_t V2Max[10] = "5"; //NbPtBins+5 +11
    Char_t MinMass[10] = "2"; //NbPtBins+5 +12
    Char_t MaxMass[10] = "5"; //NbPtBins+5 +13
    Char_t SignalF[10] = "CB2-Run2"; //NbPtBins+5 +14
    Char_t BackgroundF[20] = "DoubleExpo"; //NbPtBins+5 +15
    Char_t RatSigma[10] = "1.05"; //NbPtBins+5 +16
    
    Char_t *DefaultConditions[21];
    memset(DefaultConditions,0,21*sizeof(Char_t*));
    
    Char_t charray[10];
    
    DefaultConditions[0] = DPhiCut;
    DefaultConditions[1] = CentralityEstimator;
    DefaultConditions[2] = CentralClass;
    DefaultConditions[3] = PeriphClass;
  //  cout << NbPtBins << " is NbPtBins"<<endl;
    sprintf(charray, "%i", int(NbPtBins));
    DefaultConditions[4] = charray;
    DefaultConditions[5] = EtaMin;
    DefaultConditions[6] = EtaMax;
    DefaultConditions[7] = ZvtxCut;
    DefaultConditions[8] = EMNorm;
    DefaultConditions[9] = EMMax;
    DefaultConditions[10] = EMThreshold;
    DefaultConditions[11] = EMChange;
    DefaultConditions[12] = SummationMethod;
    DefaultConditions[13] = BackgroundV2;
    DefaultConditions[14] = V2Min;
    DefaultConditions[15] = V2Max;
    DefaultConditions[16] = MinMass;
    DefaultConditions[17] = MaxMass;
    DefaultConditions[18] = SignalF;
    DefaultConditions[19] = BackgroundF;
    DefaultConditions[20] = RatSigma;
    
    //Char_t *DefaultConditions[11] = {DPhiCut,CentralityEstimator,CentralClass,PeriphClass,NbPtbins,        EtaMin,EtaMax,ZvtxCut,EMNorm,EMMax,EMThreshold,EMChange,SummationMethod,BackgroundV2,V2Min,V2Max,MinMass,MaxMass,SignalF,BackgroundF,RatSigma};
    
    for(int index=0;index<NbPtBins+1;index++){
//0 1
//2 3        
        
        if(index==0){
                    if(std::string(row[2]) != "0-3" && std::string(row[2]) != "0-5" && std::string(row[2]) != "0-10"){
                        isSuitable = kFALSE;
                        return isSuitable;
                    }
                }
                if(index==1){
                    if(std::string(row[3]) != "40-100" && std::string(row[3]) != "40-80" && std::string(row[3]) != "30-100" && std::string(row[3]) != "50-100"){
                        isSuitable = kFALSE;
                        return isSuitable;
                    }
                }
            
        
//        if(index==2){
  //          if(std::string(row[2]) != "0-1" && std::string(row[2]) != "1-3" && std::string(row[2]) != "3-5" && std::string(row[2]) != "5-10"){
  //              isSuitable = kFALSE;
   //             return isSuitable;
    //        }
    //    }
    //    if(index==3){
      //      if(std::string(row[3]) != "0-100" && std::string(row[3]) != "20-100" && std::string(row[3]) != "40-100" && std::string(row[3]) != "60-100" && std::string(row[3]) != "20-80" && std::string(row[3]) != "40-80" && std::string(row[3]) != "60-80"){
        //        isSuitable = kFALSE;
          //      return isSuitable;
      //      }
      //  }
        if(stoi(std::string(row[5+index])) != PtBins[index]){
            cout << "row[" << 5+index << "] = " << row[5+index] << " is not the same as " << PtBins[index]<< " ABORT" <<endl;
            isSuitable = kFALSE;
            return isSuitable;
        }
    }
    
    for(int index=0; index<21; index++){
        if (strcmp(SystematicsFocus, "DPhiCut") == 0){
            if(index == 0){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
            if(index == 7){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EtaMin") == 0){
            if(index == 5){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EtaMax") == 0){
            if(index == 6){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EMNorm") == 0){
            if(index == 8){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
            if(index == 12){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "Pooling") == 0){
            if(index == 11){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "InvMass") == 0){
            if(index >= 16 && index <= 20){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "BackgroundV2") == 0){
            if(index == 13){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "RangeV2") == 0){
            if(index >= 14 && index <= 15){
                continue;
            }
        }
        if(index<=1 || index==4){
            if(row[index] != DefaultConditions[index]){
                cout << index<< "row[" << index << "] = " << row[index] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
                isSuitable = kFALSE;
                return isSuitable;
            }
        }
        else if(index>4){
            if(row[index+NbPtBins+1] != DefaultConditions[index]){
                cout << index<< "row[" << index+NbPtBins+1 << "] = " << row[index+NbPtBins+1] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
                isSuitable = kFALSE;
                return isSuitable;
            }
        }
    }
    
    return isSuitable;
    
}


void CentralityDimu(Char_t listeCSVFiles[200], Char_t SystematicsFocus[50])
{
    
    cout << "Lele"<<endl;
       
       TCanvas *c1 = new TCanvas("c1", "c1",15,49,500,800);
       const Int_t numsituations=12;
       const Int_t ptbins=6;

       Double_t sit[ptbins*numsituations];
       Double_t errsit[ptbins*numsituations];
       int counter = 0;

       for (Int_t i=0;i<numsituations*ptbins;i++){
       sit[i]=i;
       errsit[i]=0;
       };
       
       cout << "Lele1"<<endl;
        
       Double_t v2pPb[ptbins*numsituations] = {NULL};
       Double_t errv2pPb[ptbins*numsituations] = {NULL};

       cout << "Lele2"<<endl;

       std::string names[ptbins*numsituations] = {""};
       
       
       cout << "We will run a simple read of csv files -- SPDTrackletsPercentile" <<endl;
    
    
    for(int ptbin=0; ptbin<ptbins;ptbin++){
int ptcounter =0;
        cout << "\n\nPTBIN " << ptbin<<"\n\n"<<endl;
    std::ifstream file(listeCSVFiles);
    CSVRow row;
    std::string str;
    
    while (std::getline(file, str))
    {
        std::ifstream filecsv(str);
        
        while(filecsv >> row)//fixme adapt ptbin
        {
            cout << "Looking at "<<str<<endl;
            bool isVariationFound = kFALSE;
                if(isRowSuitableForCentrality(row, SystematicsFocus, "SPDClustersPercentile")){
                    cout << "\nFile: " << str <<endl;
                    cout << "ptbin: " << ptbin <<endl;
                    cout << "v2_ATLAS = " << int(10000*stod(std::string(row[4 + ptbins + 1 + 29 + 18*ptbin])))/10000. << " +- " << int(10000*stod(std::string(row[4 + ptbins + 1 + 30 + 18*ptbin])))/10000.<<endl;;
                
            v2pPb[counter] = int(10000*stod(std::string(row[4 + ptbins + 1 + 29 + 18*ptbin])))/10000.;
            errv2pPb[counter] = int(10000*stod(std::string(row[4 + ptbins + 1 + 30 + 18*ptbin])))/10000.;
            names[counter] = std::string(row[2])+"_"+std::string(row[3]);
                    cout << v2pPb[counter]<<endl;
                    cout << errv2pPb[counter]<<endl;
            cout << names[counter]<<endl;
            counter++;
		ptcounter++;
                }
            cout << "Finished looking at "<<str<<endl;
            }
    }
    double v_moy = 0;
    double v_syst = 0;
    for(int idx=counter-ptcounter; idx<counter;idx++){
        v_moy+=v2pPb[idx];
    }
    v_moy/=ptcounter;
    
    for(int idx=counter-ptcounter; idx<counter;idx++){
        v_syst+=pow(v2pPb[idx]-v_moy,2);
    }
    v_syst/=ptcounter;
    v_syst=sqrt(v_syst);
    
    cout << "Ptbin " << ptbin << " Moyenne " << v_moy<<endl;
    cout << "Ptbin " << ptbin << " Syst " << v_syst<<endl;
    }
    
    TGraphErrors *v2syst = new TGraphErrors(ptbins*numsituations,sit,v2pPb,errsit,errv2pPb);
         v2syst->SetTitle("Centrality study - ATLAS - SPDC");
        TAxis *ax = v2syst->GetHistogram()->GetXaxis();
        Double_t x1 = ax->GetBinLowEdge(1);
        Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
        v2syst->GetHistogram()->GetXaxis()->Set(ptbins*numsituations,x1-0.5,x2+0.5);
         v2syst->GetHistogram()->GetYaxis()->SetRangeUser(-0.1,0.1);
         v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");
        v2syst->SetMarkerStyle(20);

        for(Int_t k=0;k<ptbins*numsituations;k++){
        v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
        }
        
        v2syst->Draw("AP");
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop");
        sprintf(CanvasName,"%s/CentralitiesDimu_ATLAS_SPDC_ChecknopropAAA.pdf",FolderName);
        c1->SaveAs(CanvasName);
    
}

void CentralityPlotDimu(Char_t listeCSVFiles[200], Char_t SystematicsFocus[50])
{
    
    cout << "Lele"<<endl;
       
       const Int_t numsituations=7;
       const Int_t ptbins=2;

       Double_t sit[numsituations];
       Double_t errsit[numsituations];
       int counter = 0;

       for (Int_t i=0;i<numsituations;i++){
       sit[i]=i;
       errsit[i]=0;
       };
       
       cout << "Lele1"<<endl;
        
       Double_t v2pPb[numsituations] = {NULL};
       Double_t errv2pPb[numsituations] = {NULL};

       cout << "Lele2"<<endl;

       std::string names[numsituations] = {""};
       
       
       cout << "We will run a simple read of csv files -- SPDTrackletsPercentile" <<endl;
    
    
    for(int ptbin=0; ptbin<ptbins;ptbin++){
        counter = 0;
        TCanvas *c1 = new TCanvas("c1", "c1",0,0,2000,1000);
        c1 -> SetBottomMargin(0.2);
        cout << "\n\nPTBIN " << ptbin<<"\n\n"<<endl;
    std::ifstream file(listeCSVFiles);
    CSVRow row;
    std::string str;
    
    while (std::getline(file, str))
    {
        std::ifstream filecsv(str);
        
        while(filecsv >> row)//fixme adapt ptbin
        {
            cout << "Looking at "<<str<<endl;
            bool isVariationFound = kFALSE;
                if(isRowSuitableForCentrality(row, SystematicsFocus, "SPDTrackletsPercentile")){
                    cout << "\nFile: " << str <<endl;//17 18
                    cout << "ptbin: " << ptbin <<endl;//29 30
                    cout << "v2_ATLAS = " << int(10000*stod(std::string(row[4 + ptbins + 1 + 17 + 18*ptbin])))/10000. << " +- " << int(10000*stod(std::string(row[4 + ptbins + 1 + 18 + 18*ptbin])))/10000.<<endl;;
                
            v2pPb[counter] = int(10000*stod(std::string(row[4 + ptbins + 1 + 17 + 18*ptbin])))/10000.;
            errv2pPb[counter] = int(10000*stod(std::string(row[4 + ptbins + 1 + 18 + 18*ptbin])))/10000.;
            names[counter] = std::string(row[2])+"_"+std::string(row[3]);
                    cout << v2pPb[counter]<<endl;
                    cout << errv2pPb[counter]<<endl;
            cout << names[counter]<<endl;
            counter++;
                }
           // cout << "Finished looking at "<<str<<endl;
            }
    }
            swap(names[2], names[5]);
            swap(names[3], names[6]);
            swap(names[4], names[5]);
            swap(names[5], names[6]);

        swap(v2pPb[2], v2pPb[5]);
        swap(v2pPb[3], v2pPb[6]);
        swap(v2pPb[4], v2pPb[5]);
        swap(v2pPb[5], v2pPb[6]);
        
        swap(errv2pPb[2], errv2pPb[5]);
        swap(errv2pPb[3], errv2pPb[6]);
        swap(errv2pPb[4], errv2pPb[5]);
        swap(errv2pPb[5], errv2pPb[6]);
        
        
        TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
        Char_t TitlePlot[500];
        sprintf(TitlePlot,"Centrality study - YieldSub - SPDTracklets - PtBin %i", ptbin);
         v2syst->SetTitle(TitlePlot);
        TAxis *ax = v2syst->GetHistogram()->GetXaxis();
        Double_t x1 = ax->GetBinLowEdge(1);
        Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
        v2syst->GetHistogram()->GetXaxis()->Set(numsituations*5,x1-0.5,x2+0.5);
         v2syst->GetHistogram()->GetYaxis()->SetRangeUser(-0.2,0.2);
         v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,J/#psi}");
        v2syst->SetMarkerStyle(20);

        for(Int_t k=0;k<numsituations;k++){
            int bin_index = v2syst->GetHistogram()->GetXaxis()->FindBin(k);
        v2syst->GetHistogram()->GetXaxis()->SetBinLabel(bin_index,names[k].c_str());
        }
        
        v2syst->GetHistogram()->GetXaxis()->LabelsOption("v");
        v2syst->GetHistogram()->GetXaxis()->SetLabelSize(0.04);
        
        v2syst->Draw("AP");
        
        TLine *l=new TLine(0.,0.0,numsituations,0.0);
        l->SetLineColor(kBlack);
        l->SetLineWidth(1);
        l->SetLineStyle(9);
        l->Draw("same");
        
        for(int idx=1; idx<=3; idx++){
        
        TLine *lv=new TLine(idx*7-0.5,-0.2,idx*7-0.5,0.2);
        lv->SetLineColor(kGray);
        lv->SetLineWidth(1);
        lv->SetLineStyle(1);
        lv->Draw("same");
            
        }
        
        
        Char_t CanvasName[500];
        Char_t FolderName[500];
        sprintf(FolderName,"~/Desktop");
        sprintf(CanvasName,"%s/CentralitiesStudyFullestDimu_YieldSub_intSPDT_PtBin%i.pdf",FolderName, ptbin);
        c1->SaveAs(CanvasName);
        
    }
    
}





