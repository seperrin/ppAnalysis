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
# include "TGraph2D.h"
 # include "TVirtualFitter.h"
 # include "TBackCompFitter.h"
 # include "TMultiGraph.h"
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

void SystematicsTKL(Char_t ListeCSVFiles[100], Char_t SystematicsFocus[50]);
bool isRowSuitableForSystematics(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]);
void CentralityStudyTKL(Char_t ListeCSVFiles[100], Char_t SystematicsFocus[50]);
bool isRowSuitableForCentrality(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]);
void CentralityStudyTKLPYTHIA(Char_t ListeCSVFiles[100], Char_t SystematicsFocus[50]);
bool isRowSuitableForCentralityPYTHIA(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]);

void SystematicsTKL(Char_t listeCSVFiles[100], Char_t SystematicsFocus[50])
{
    std::ifstream file(listeCSVFiles);
    CSVRow row;
    std::string str;
    
    double v2pPbSPDT[10]={0};
    double errv2pPbSPDT[10]={0};
    double v2pPbSPDC[10]={0};
    double errv2pPbSPDC[10]={0};
    double v2pPbV0M[10]={0};
    double errv2pPbV0M[10]={0};
    double v2pPbMixed[10]={0};
    double errv2pPbMixed[10]={0};
    
    double v2ATLASSPDT[10]={0};
    double errv2ATLASSPDT[10]={0};
    double v2ATLASSPDC[10]={0};
    double errv2ATLASSPDC[10]={0};
    double v2ATLASV0M[10]={0};
    double errv2ATLASV0M[10]={0};
    double v2ATLASMixed[10]={0};
    double errv2ATLASMixed[10]={0};
    
    int IndexOfVariable;
    int IndexOfStorage;
    
    string variableDPhiCut[5] = {"1mrad", "2mrad", "5mrad", "10mrad", "None"};
    string variableZvtxCut[3] = {"8","10","12"};
    string variableEtaGap[3] = {"1.0","1.2","1.4"};
    string variableEMNorm[3] = {"Method1","Method2"};
    string variableSummationZvtx[3] = {"Method1c","Method1a","Method2"};
   // string variableExtractionMethod[3] = {"Method1","Method2","Method3"};
    string variableZYAM[2] = {"ZYAM","noZYAM"};
    string variableFitMethod[3] = {"V2","V3","Gauss"};
    string variableATLAS[2] = {"Cvetan","ATLAS"};
    
    bool isClassicVariation = kFALSE;
    
    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
        IndexOfVariable = 0;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
        IndexOfVariable = 5;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "EtaGap") == 0){
        IndexOfVariable = 4;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "EMNorm") == 0){
        IndexOfVariable = 6;
        isClassicVariation = kTRUE;
    }
    if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
        IndexOfVariable = 10;
        isClassicVariation = kTRUE;
    }
    
    cout << "We will run a simple read of csv files" <<endl;
    
    while (std::getline(file, str))
    {
        std::ifstream filecsv(str);
        while(filecsv >> row)
        {
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
                    if(strcmp(SystematicsFocus, "EtaGap") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaGap[idx] == row[IndexOfVariable]){
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
                    
                    if(isVariationFound){
                        v2pPbV0M[IndexOfStorage] = stod(std::string(row[11]));
                        errv2pPbV0M[IndexOfStorage] = stod(std::string(row[12]));
                        v2ATLASV0M[IndexOfStorage] = stod(std::string(row[19]));
                        errv2ATLASV0M[IndexOfStorage] = stod(std::string(row[20]));
                    }
                    else if(strcmp(SystematicsFocus, "ZYAM") == 0){
                        v2pPbV0M[0] = stod(std::string(row[11]));
                        errv2pPbV0M[0] = stod(std::string(row[12]));
                        v2pPbV0M[1] = stod(std::string(row[13]));
                        errv2pPbV0M[1] = stod(std::string(row[14]));
                    }
                    else if(strcmp(SystematicsFocus, "ATLAS") == 0){
                        v2pPbV0M[0] = stod(std::string(row[11]));
                        errv2pPbV0M[0] = stod(std::string(row[12]));
                        v2ATLASV0M[0] = stod(std::string(row[19]));
                        errv2ATLASV0M[0] = stod(std::string(row[20]));
                    }
                    else if(strcmp(SystematicsFocus, "FitMethod") == 0){
                        v2pPbV0M[0] = stod(std::string(row[11]));
                        errv2pPbV0M[0] = stod(std::string(row[12]));
                        v2pPbV0M[1] = stod(std::string(row[23]));
                        errv2pPbV0M[1] = stod(std::string(row[24]));
                        v2pPbV0M[2] = stod(std::string(row[25]));
                        errv2pPbV0M[2] = stod(std::string(row[26]));
                    }
                }
                else if(isRowSuitableForSystematics(row, SystematicsFocus, "SPDTrackletsPercentile")){
                    cout << "Keeping file " << str <<endl;
                    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
                        for(int idx=0; idx<5; idx++){
                            if (variableDPhiCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
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
                    if(strcmp(SystematicsFocus, "EtaGap") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaGap[idx] == row[IndexOfVariable]){
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
                    
                    if(isVariationFound){
                        v2pPbSPDT[IndexOfStorage] = stod(std::string(row[11]));
                        errv2pPbSPDT[IndexOfStorage] = stod(std::string(row[12]));
                        v2ATLASSPDT[IndexOfStorage] = stod(std::string(row[19]));
                        errv2ATLASSPDT[IndexOfStorage] = stod(std::string(row[20]));
                    }
                    else if(strcmp(SystematicsFocus, "ZYAM") == 0){
                        v2pPbSPDT[0] = stod(std::string(row[11]));
                        errv2pPbSPDT[0] = stod(std::string(row[12]));
                        v2pPbSPDT[1] = stod(std::string(row[13]));
                        errv2pPbSPDT[1] = stod(std::string(row[14]));
                    }
                    else if(strcmp(SystematicsFocus, "ATLAS") == 0){
                       v2pPbSPDT[0] = stod(std::string(row[11]));
                       errv2pPbSPDT[0] = stod(std::string(row[12]));
                       v2ATLASSPDT[0] = stod(std::string(row[19]));
                       errv2ATLASSPDT[0] = stod(std::string(row[20]));
                   }
                    else if(strcmp(SystematicsFocus, "FitMethod") == 0){
                        v2pPbSPDT[0] = stod(std::string(row[11]));
                        errv2pPbSPDT[0] = stod(std::string(row[12]));
                        v2pPbSPDT[1] = stod(std::string(row[23]));
                        errv2pPbSPDT[1] = stod(std::string(row[24]));
                        v2pPbSPDT[2] = stod(std::string(row[25]));
                        errv2pPbSPDT[2] = stod(std::string(row[26]));
                    }
                }
                else if(isRowSuitableForSystematics(row, SystematicsFocus, "SPDClustersPercentile")){
                    cout << "Keeping file " << str <<endl;
                    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
                        for(int idx=0; idx<5; idx++){
                            if (variableDPhiCut[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
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
                    if(strcmp(SystematicsFocus, "EtaGap") == 0){
                        for(int idx=0; idx<3; idx++){
                            if (variableEtaGap[idx] == row[IndexOfVariable]){
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
                    
                    if(isVariationFound){
                        v2pPbSPDC[IndexOfStorage] = stod(std::string(row[11]));
                        errv2pPbSPDC[IndexOfStorage] = stod(std::string(row[12]));
                        v2ATLASSPDC[IndexOfStorage] = stod(std::string(row[19]));
                        errv2ATLASSPDC[IndexOfStorage] = stod(std::string(row[20]));
                    }
                    else if(strcmp(SystematicsFocus, "ZYAM") == 0){
                        v2pPbSPDC[0] = stod(std::string(row[11]));
                        errv2pPbSPDC[0] = stod(std::string(row[12]));
                        v2pPbSPDC[1] = stod(std::string(row[13]));
                        errv2pPbSPDC[1] = stod(std::string(row[14]));
                    }
                    else if(strcmp(SystematicsFocus, "ATLAS") == 0){
                        v2pPbSPDC[0] = stod(std::string(row[11]));
                        errv2pPbSPDC[0] = stod(std::string(row[12]));
                        v2ATLASSPDC[0] = stod(std::string(row[19]));
                        errv2ATLASSPDC[0] = stod(std::string(row[20]));
                    }
                    else if(strcmp(SystematicsFocus, "FitMethod") == 0){
                        v2pPbSPDC[0] = stod(std::string(row[11]));
                        errv2pPbSPDC[0] = stod(std::string(row[12]));
                        v2pPbSPDC[1] = stod(std::string(row[23]));
                        errv2pPbSPDC[1] = stod(std::string(row[24]));
                        v2pPbSPDC[2] = stod(std::string(row[25]));
                        errv2pPbSPDC[2] = stod(std::string(row[26]));
                    }
                }
            else if(isRowSuitableForSystematics(row, SystematicsFocus, "V0MPlusSPD16h")){
                cout << "Keeping file " << str <<endl;
                if (strcmp(SystematicsFocus, "DPhiCut") == 0){
                    for(int idx=0; idx<5; idx++){
                        if (variableDPhiCut[idx] == row[IndexOfVariable]){
                            IndexOfStorage = idx;
                            isVariationFound = kTRUE;
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
                if(strcmp(SystematicsFocus, "EtaGap") == 0){
                    for(int idx=0; idx<3; idx++){
                        if (variableEtaGap[idx] == row[IndexOfVariable]){
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
                
                if(isVariationFound){
                    v2pPbMixed[IndexOfStorage] = stod(std::string(row[11]));
                    errv2pPbMixed[IndexOfStorage] = stod(std::string(row[12]));
                    v2ATLASMixed[IndexOfStorage] = stod(std::string(row[19]));
                    errv2ATLASMixed[IndexOfStorage] = stod(std::string(row[20]));
                }
                else if(strcmp(SystematicsFocus, "ZYAM") == 0){
                    v2pPbMixed[0] = stod(std::string(row[11]));
                    errv2pPbMixed[0] = stod(std::string(row[12]));
                    v2pPbMixed[1] = stod(std::string(row[13]));
                    errv2pPbMixed[1] = stod(std::string(row[14]));
                }
                else if(strcmp(SystematicsFocus, "ATLAS") == 0){
                    v2pPbMixed[0] = stod(std::string(row[11]));
                    errv2pPbMixed[0] = stod(std::string(row[12]));
                    v2ATLASMixed[0] = stod(std::string(row[19]));
                    errv2ATLASMixed[0] = stod(std::string(row[20]));
                }
                else if(strcmp(SystematicsFocus, "FitMethod") == 0){
                    v2pPbMixed[0] = stod(std::string(row[11]));
                    errv2pPbMixed[0] = stod(std::string(row[12]));
                    v2pPbMixed[1] = stod(std::string(row[23]));
                    errv2pPbMixed[1] = stod(std::string(row[24]));
                    v2pPbMixed[2] = stod(std::string(row[25]));
                    errv2pPbMixed[2] = stod(std::string(row[26]));
                }
            }
                else{
                    continue;
                }
            }
        
    }
    
    cout << "\n\n\n ==================== \n Systematics results \n ==================== \n\n\n"<<endl;
    if(isClassicVariation){
        int default_idx = 0;
        if(strcmp(SystematicsFocus, "DPhiCut") == 0){
            default_idx = 3;
        }
        if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
            default_idx = 1;
        }
        if(strcmp(SystematicsFocus, "EtaGap") == 0){
            default_idx = 1;
        }
        if(strcmp(SystematicsFocus, "EMNorm") == 0){
            default_idx = 0;
        }
        if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
            default_idx = 0;
        }
        cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
        double v2_min = 999.;
        double v2_max = -999.;
        int nbMesures = 0;
                
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDT[idx]!=0 || errv2pPbSPDT[idx]!=0){
            cout << "v2: "<<v2pPbSPDT[idx]<< " +/- "<<errv2pPbSPDT[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbSPDT[idx];
                if(v2pPbSPDT[idx]>v2_max){
                    v2_max = v2pPbSPDT[idx];
                }
                if(v2pPbSPDT[idx]<v2_min){
                    v2_min = v2pPbSPDT[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDT[idx]!=0 || errv2pPbSPDT[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbSPDT[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbSPDT[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "ATLAS Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2ATLASSPDT[idx]!=0 || errv2ATLASSPDT[idx]!=0){
            cout << "v2: "<<v2ATLASSPDT[idx]<< " +/- "<<errv2ATLASSPDT[idx] <<endl;
                nbMesures++;
                v2_moy+=v2ATLASSPDT[idx];
                if(v2ATLASSPDT[idx]>v2_max){
                    v2_max = v2ATLASSPDT[idx];
                }
                if(v2ATLASSPDT[idx]<v2_min){
                    v2_min = v2ATLASSPDT[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2ATLASSPDT[idx]!=0 || errv2ATLASSPDT[idx]!=0){
                v2_syst+=pow(v2_moy-v2ATLASSPDT[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2ATLASSPDT[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        }
        
        
        
        cout << "Centrality estimator : MixedEstimator"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
        double v2_min = 999.;
        double v2_max = -999.;
        int nbMesures = 0;
                
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbMixed[idx]!=0 || errv2pPbMixed[idx]!=0){
            cout << "v2: "<<v2pPbMixed[idx]<< " +/- "<<errv2pPbMixed[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbMixed[idx];
                if(v2pPbMixed[idx]>v2_max){
                    v2_max = v2pPbMixed[idx];
                }
                if(v2pPbMixed[idx]<v2_min){
                    v2_min = v2pPbMixed[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbMixed[idx]!=0 || errv2pPbMixed[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbMixed[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbMixed[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "ATLAS Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2ATLASMixed[idx]!=0 || errv2ATLASMixed[idx]!=0){
            cout << "v2: "<<v2ATLASMixed[idx]<< " +/- "<<errv2ATLASMixed[idx] <<endl;
                nbMesures++;
                v2_moy+=v2ATLASMixed[idx];
                if(v2ATLASMixed[idx]>v2_max){
                    v2_max = v2ATLASMixed[idx];
                }
                if(v2ATLASMixed[idx]<v2_min){
                    v2_min = v2ATLASMixed[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2ATLASMixed[idx]!=0 || errv2ATLASMixed[idx]!=0){
                v2_syst+=pow(v2_moy-v2ATLASMixed[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2ATLASMixed[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        }
        
        
        cout << "Centrality estimator : SPDClustersPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDC[idx]!=0 || errv2pPbSPDC[idx]!=0){
            cout << "v2: "<<v2pPbSPDC[idx]<< " +/- "<<errv2pPbSPDC[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbSPDC[idx];
                if(v2pPbSPDC[idx]>v2_max){
                    v2_max = v2pPbSPDC[idx];
                }
                if(v2pPbSPDC[idx]<v2_min){
                    v2_min = v2pPbSPDC[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDC[idx]!=0 || errv2pPbSPDC[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbSPDC[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbSPDC[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "ATLAS Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2ATLASSPDC[idx]!=0 || errv2ATLASSPDC[idx]!=0){
            cout << "v2: "<<v2ATLASSPDC[idx]<< " +/- "<<errv2ATLASSPDC[idx] <<endl;
                nbMesures++;
                v2_moy+=v2ATLASSPDC[idx];
                if(v2ATLASSPDC[idx]>v2_max){
                    v2_max = v2ATLASSPDC[idx];
                }
                if(v2ATLASSPDC[idx]<v2_min){
                    v2_min = v2ATLASSPDC[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2ATLASSPDC[idx]!=0 || errv2ATLASSPDC[idx]!=0){
                v2_syst+=pow(v2_moy-v2ATLASSPDC[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2ATLASSPDC[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        }
        
        
        cout << "Centrality estimator : V0MPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbV0M[idx]!=0 || errv2pPbV0M[idx]!=0){
            cout << "v2: "<<v2pPbV0M[idx]<< " +/- "<<errv2pPbV0M[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbV0M[idx];
                if(v2pPbV0M[idx]>v2_max){
                    v2_max = v2pPbV0M[idx];
                }
                if(v2pPbV0M[idx]<v2_min){
                    v2_min = v2pPbV0M[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbV0M[idx]!=0 || errv2pPbV0M[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbV0M[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbV0M[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "ATLAS Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2ATLASV0M[idx]!=0 || errv2ATLASV0M[idx]!=0){
            cout << "v2: "<<v2ATLASV0M[idx]<< " +/- "<<errv2ATLASV0M[idx] <<endl;
                nbMesures++;
                v2_moy+=v2ATLASV0M[idx];
                if(v2ATLASV0M[idx]>v2_max){
                    v2_max = v2ATLASV0M[idx];
                }
                if(v2ATLASV0M[idx]<v2_min){
                    v2_min = v2ATLASV0M[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2ATLASV0M[idx]!=0 || errv2ATLASV0M[idx]!=0){
                v2_syst+=pow(v2_moy-v2ATLASV0M[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2ATLASV0M[default_idx];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        }
    }
    else if(strcmp(SystematicsFocus, "ZYAM") == 0){
        cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbSPDT[0];
            v2_syst = abs(v2pPbSPDT[1]-v2pPbSPDT[0]);
            cout << "v2_noZYAM: "<<v2pPbSPDT[1]<< " +/- "<<errv2pPbSPDT[1] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : SPDClustersPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbSPDC[0];
            v2_syst = abs(v2pPbSPDC[1]-v2pPbSPDC[0]);
            cout << "v2_noZYAM: "<<v2pPbSPDC[1]<< " +/- "<<errv2pPbSPDC[1] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : V0MPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbV0M[0];
            v2_syst = abs(v2pPbV0M[1]-v2pPbV0M[0]);
            cout << "v2_noZYAM: "<<v2pPbV0M[1]<< " +/- "<<errv2pPbV0M[1] <<endl;
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
            v2_default = v2pPbSPDT[0];
            v2_syst = abs(v2ATLASSPDT[0]-v2pPbSPDT[0]);
            cout << "v2_ATLAS: "<<v2ATLASSPDT[0]<< " +/- "<<errv2ATLASSPDT[0] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : SPDClustersPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbSPDC[0];
            v2_syst = abs(v2ATLASSPDC[0]-v2pPbSPDC[0]);
            cout << "v2_ATLAS: "<<v2ATLASSPDC[0]<< " +/- "<<errv2ATLASSPDC[0] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : V0MPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_default = 0;
        double v2_syst = 0;
        int nbMesures = 0;
            v2_default = v2pPbV0M[0];
            v2_syst = abs(v2ATLASV0M[0]-v2pPbV0M[0]);
            cout << "v2_ATLAS: "<<v2ATLASV0M[0]<< " +/- "<<errv2ATLASV0M[0] <<endl;
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
    }
    
    else if(strcmp(SystematicsFocus, "FitMethod") == 0){
        cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDT[idx]!=0 || errv2pPbSPDT[idx]!=0){
            cout << "v2: "<<v2pPbSPDT[idx]<< " +/- "<<errv2pPbSPDT[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbSPDT[idx];
                if(v2pPbSPDT[idx]>v2_max){
                    v2_max = v2pPbSPDT[idx];
                }
                if(v2pPbSPDT[idx]<v2_min){
                    v2_min = v2pPbSPDT[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDT[idx]!=0 || errv2pPbSPDT[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbSPDT[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbSPDT[0];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }

        cout << "Centrality estimator : SPDClustersPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDC[idx]!=0 || errv2pPbSPDC[idx]!=0){
            cout << "v2: "<<v2pPbSPDC[idx]<< " +/- "<<errv2pPbSPDC[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbSPDC[idx];
                if(v2pPbSPDC[idx]>v2_max){
                    v2_max = v2pPbSPDC[idx];
                }
                if(v2pPbSPDC[idx]<v2_min){
                    v2_min = v2pPbSPDC[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbSPDC[idx]!=0 || errv2pPbSPDC[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbSPDC[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbSPDC[0];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
        
        cout << "Centrality estimator : MixedEstimator"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbMixed[idx]!=0 || errv2pPbMixed[idx]!=0){
            cout << "v2: "<<v2pPbMixed[idx]<< " +/- "<<errv2pPbMixed[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbMixed[idx];
                if(v2pPbMixed[idx]>v2_max){
                    v2_max = v2pPbMixed[idx];
                }
                if(v2pPbMixed[idx]<v2_min){
                    v2_min = v2pPbMixed[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbMixed[idx]!=0 || errv2pPbMixed[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbMixed[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbMixed[0];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }

        cout << "Centrality estimator : V0MPercentile"<<endl;
        
        cout << "p-Pb Method"<<endl;
        {
        double v2_moy = 0;
        double v2_syst = 0;
            double v2_min = 999.;
            double v2_max = -999.;
        int nbMesures = 0;
        
        for(int idx=0; idx<10; idx++){
            if(v2pPbV0M[idx]!=0 || errv2pPbV0M[idx]!=0){
            cout << "v2: "<<v2pPbV0M[idx]<< " +/- "<<errv2pPbV0M[idx] <<endl;
                nbMesures++;
                v2_moy+=v2pPbV0M[idx];
                if(v2pPbV0M[idx]>v2_max){
                    v2_max = v2pPbV0M[idx];
                }
                if(v2pPbV0M[idx]<v2_min){
                    v2_min = v2pPbV0M[idx];
                }
            }
        }
        v2_moy/=nbMesures;
        for(int idx=0; idx<10; idx++){
            if(v2pPbV0M[idx]!=0 || errv2pPbV0M[idx]!=0){
                v2_syst+=pow(v2_moy-v2pPbV0M[idx],2);
            }
        }
            v2_syst/=nbMesures;
        v2_syst = sqrt(v2_syst);
            double v2_default = v2pPbV0M[0];
        cout << "Systematics: v2_default = " <<v2_default<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_default)*100<<"%)"<<endl;
            cout << "Systematics Distance: v2_default = " <<v2_default<<" systematic = "<<v2_max-v2_min<<" ("<< ((v2_max-v2_min)/v2_default)*100<<"%)"<<endl;
        //Plot et calculer la systematique
        }
    }
}

void CentralityStudyTKL(Char_t listeCSVFiles[200], Char_t SystematicsFocus[50])
{
    std::ifstream file(listeCSVFiles);
    CSVRow row;
    std::string str;
    
    
    cout << "Lele"<<endl;
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);
    const Int_t numsituations=12;

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
    
    
    cout << "We will run a simple read of csv files -- Mixed" <<endl;
    
    while (std::getline(file, str))
    {
        std::ifstream filecsv(str);
        
        while(filecsv >> row)
        {
           // cout << "Looking at "<<str<<endl;
            bool isVariationFound = kFALSE;
                if(isRowSuitableForCentrality(row, SystematicsFocus, "V0MPlusSPD16h")){
                    cout << "\n\nFile: " << str <<endl;
                    cout << "v2_ATLAS = " << int(10000*stod(std::string(row[11])))/10000. << " +- " << int(10000*stod(std::string(row[12])))/10000.<<endl;;
                
            v2pPb[counter] = int(10000*stod(std::string(row[11])))/10000.;
            errv2pPb[counter] = int(10000*stod(std::string(row[12])))/10000.;
            names[counter] = std::string(row[2])+"_"+std::string(row[3]);
                    cout << v2pPb[counter]<<endl;
                    cout << errv2pPb[counter]<<endl;
            cout << names[counter]<<endl;
            counter++;
                }
           // cout << "Finished looking at "<<str<<endl;
            }
    }
    double v_moy = 0;
    double v_syst = 0;
    for(int idx=0; idx<counter;idx++){
        v_moy+=v2pPb[idx];
    }
    v_moy/=counter;
    
    for(int idx=0; idx<counter;idx++){
        v_syst+=pow(v2pPb[idx]-v_moy,2);
    }
    v_syst/=counter;
    v_syst=sqrt(v_syst);
    
    cout << "Moyenne " << v_moy<<endl;
    cout << "Syst " << v_syst<<endl;
    
    
    TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,v2pPb,errsit,errv2pPb);
     v2syst->SetTitle("Centrality study - YieldSub - Mixed");
    TAxis *ax = v2syst->GetHistogram()->GetXaxis();
    Double_t x1 = ax->GetBinLowEdge(1);
    Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
    v2syst->GetHistogram()->GetXaxis()->Set(numsituations,x1-0.5,x2+0.5);
     v2syst->GetHistogram()->GetYaxis()->SetRangeUser(0.04,0.15);
     v2syst->GetHistogram()->GetYaxis()->SetTitle("v_{2,tkl}");
    v2syst->SetMarkerStyle(20);

    for(Int_t k=0;k<numsituations;k++){
    v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
    }
    
    v2syst->Draw("AP");
    
    Char_t CanvasName[500];
    Char_t FolderName[500];
    sprintf(FolderName,"~/Desktop");
    sprintf(CanvasName,"%s/CentralitiesTKL_YieldSub_Mixed.pdf",FolderName);
    c1->SaveAs(CanvasName);
    
}







bool isRowSuitableForSystematics(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]){
    //Dire si la row est utile pour la systematique que l'on veut
    
    bool isSuitable = kTRUE;
    //Default values
    
    Char_t DPhiCut[10] = "10mrad"; //0
    Char_t CentralClass[10] = "0-5";  //2
    Char_t PeriphClass[10] = "40-100"; //3
    Char_t EtaGap[10] = "1.2"; //4
    Char_t ZvtxCut[10] = "10"; //5
    Char_t EMNorm[10] = "Method1"; //6
    Char_t EMMax[10] = "30"; //7
    Char_t EMThreshold[10] = "10"; //8
    Char_t PeriphScale[10] = "1"; //9
    Char_t SummationMethod[10] = "Method1c"; //10
    
    Char_t *DefaultConditions[11] = {DPhiCut,CentralityEstimator,CentralClass,PeriphClass,EtaGap,ZvtxCut,EMNorm,EMMax,EMThreshold,PeriphScale,SummationMethod};
    
    for(int index=0; index<11; index++){
        if (strcmp(SystematicsFocus, "DPhiCut") == 0){
            if(index == 0){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
            if(index == 5){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EtaGap") == 0){
            if(index == 4){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "EMNorm") == 0){
            if(index == 6){
                continue;
            }
        }
        if(strcmp(SystematicsFocus, "SummationZvtx") == 0){
            if(index == 10){
                continue;
            }
        }
        if(row[index] != DefaultConditions[index]){
          //  cout << "row[" << index << "] = " << row[index] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
            isSuitable = kFALSE;
            return isSuitable;
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
    Char_t EtaGap[10] = "1.2"; //4
    Char_t ZvtxCut[10] = "10"; //5
    Char_t EMNorm[10] = "Method1"; //6
    Char_t EMMax[10] = "30"; //7
    Char_t EMThreshold[10] = "10"; //8
    Char_t PeriphScale[10] = "1"; //9
    Char_t SummationMethod[10] = "Method1c"; //10
    
    Char_t *DefaultConditions[11] = {DPhiCut,CentralityEstimator,CentralClass,PeriphClass,EtaGap,ZvtxCut,EMNorm,EMMax,EMThreshold,PeriphScale,SummationMethod};
    
    for(int index=0; index<11; index++){
       // cout << "row[" << index << "] = " << row[index] << " comapred with " << DefaultConditions[index]<< " " <<endl;
        if(index == 2){
            if(row[index] != "0-3" && row[index] != "0-5" && row[index] != "0-10"){
            isSuitable = kFALSE;
            return isSuitable;
            }
        }
        else if(index == 3){
            if(row[index] != "30-100" && row[index] != "40-100" && row[index] != "50-100" && row[index] != "40-80"){
            isSuitable = kFALSE;
            return isSuitable;
            }
        }
        else if(row[index] != DefaultConditions[index]){
           // cout << "row[" << index << "] = " << row[index] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
            isSuitable = kFALSE;
            return isSuitable;
        }
    }
    
    return isSuitable;
    
}

bool isRowSuitableForCentralityPYTHIA(CSVRow row, Char_t SystematicsFocus[50], Char_t CentralityEstimator[50]){
    //Dire si la row est utile pour la systematique que l'on veut
    
    bool isSuitable = kTRUE;
    //Default values
    
    Char_t DPhiCut[10] = "10mrad"; //0
    Char_t CentralClass[10] = "140-130";  //2
    Char_t PeriphClass[10] = "10-0"; //3
    Char_t EtaGap[10] = "1.2"; //4
    Char_t ZvtxCut[10] = "10"; //5
    Char_t EMNorm[10] = "Method1"; //6
    Char_t EMMax[10] = "30"; //7
    Char_t EMThreshold[10] = "10"; //8
    Char_t PeriphScale[10] = "1"; //9
    Char_t SummationMethod[10] = "Method1c"; //10
    
    Char_t *DefaultConditions[11] = {DPhiCut,CentralityEstimator,CentralClass,PeriphClass,EtaGap,ZvtxCut,EMNorm,EMMax,EMThreshold,PeriphScale,SummationMethod};
    
    //Char_t *DefaultConditions[11] = {DPhiCut,CentralityEstimator,CentralClass,PeriphClass,NbPtbins,        EtaMin,EtaMax,ZvtxCut,EMNorm,EMMax,EMThreshold,EMChange,SummationMethod,BackgroundV2,V2Min,V2Max,MinMass,MaxMass,SignalF,BackgroundF,RatSigma};
    
    for(int index=0; index<11; index++){
       // cout << "row[" << index << "] = " << row[index] << " comapred with " << DefaultConditions[index]<< " " <<endl;
        if(index == 2){
            if(row[index] != "140-130" && row[index] != "130-120" && row[index] != "120-110" && row[index] != "110-100" && row[index] != "100-90" && row[index] != "90-80" && row[index] != "80-70" && row[index] != "70-60" && row[index] != "60-50" && row[index] != "50-40" && row[index] != "40-30" && row[index] != "30-20"){
            isSuitable = kFALSE;
            return isSuitable;
            }
        }
        else if(index == 3){
            if(row[index] != "5-0"){
            isSuitable = kFALSE;
            return isSuitable;
            }
        }
        else if(row[index] != DefaultConditions[index]){
           // cout << "row[" << index << "] = " << row[index] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
            isSuitable = kFALSE;
            return isSuitable;
        }
    }
    
    return isSuitable;
    
}


void CentralityStudyTKLPYTHIA(Char_t listeCSVFiles[200], Char_t SystematicsFocus[50])
{
    std::ifstream file(listeCSVFiles);
    CSVRow row;
    std::string str;
    
    
    cout << "Lele"<<endl;
    
    TCanvas *c1 = new TCanvas("c1", "c1",15,49,1051,500);
    const Int_t numsituations=12;

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
    Double_t V22pPb[numsituations] = {NULL};
    Double_t errV22pPb[numsituations] = {NULL};


    cout << "Lele2"<<endl;

    std::string names[numsituations] = {""};
    
    
    cout << "We will run a simple read of csv files -- SPDTrackletsPercentile" <<endl;
    
    while (std::getline(file, str))
    {
        std::ifstream filecsv(str);
        
        while(filecsv >> row)
        {
           // cout << "Looking at "<<str<<endl;
            bool isVariationFound = kFALSE;
                if(isRowSuitableForCentralityPYTHIA(row, SystematicsFocus, "V0MPercentile")){
                    cout << "\n\nFile: " << str <<endl;
                    cout << "v2_ATLAS = " << int(10000*stod(std::string(row[19])))/10000. << " +- " << int(10000*stod(std::string(row[20])))/10000.<<endl; //ATLAS 19 20
                
            v2pPb[counter] = int(10000*stod(std::string(row[19])))/10000.;
            errv2pPb[counter] = int(10000*stod(std::string(row[20])))/10000.;
                    V22pPb[counter] = int(10000*stod(std::string(row[19])))/10000.;
                    errV22pPb[counter] = int(10000*stod(std::string(row[20])))/10000.;
            names[counter] = std::string(row[2])+"_"+std::string(row[3]);
                    cout << v2pPb[counter]<<endl;
                    cout << errv2pPb[counter]<<endl;
            cout << names[counter]<<endl;
                    
                    Char_t lessCentral[5];
                            Char_t mostCentral[5];
                           // sprintf(tempy, "%s", std::string(row[2]))
                    sscanf(std::string(row[2]).c_str(), "%[^-]-%s", mostCentral,lessCentral);
                            int mostCentralNum;
                            stringstream intValue2(mostCentral);
                            intValue2 >> mostCentralNum;
                            int lessCentralNum;
                            stringstream intValue3(lessCentral);
                            intValue3 >> lessCentralNum;
                            
                    sit[counter] = (mostCentralNum+lessCentralNum)/2.;
                    errsit[counter] = TMath::Abs((mostCentralNum-lessCentralNum)/2.);
                    if(V22pPb[counter]==0){
//                        if(TMath::Abs(sit[counter]-25) <0.1){
//                            V22pPb[counter] = -0.002787;
//                            errV22pPb[counter] = 0.004059;
//                        }
//                        else if(TMath::Abs(sit[counter]-35) <0.1){
//                            V22pPb[counter] = -0.001009;
//                            errV22pPb[counter] = 0.001304;
//                        }
//                        else if(TMath::Abs(sit[counter]-45) <0.1){
//                            V22pPb[counter] = -0.0001187;
//                            errV22pPb[counter] = 0.0008000;
//                        }
//                        else if(TMath::Abs(sit[counter]-55) <0.1){
//                            V22pPb[counter] = -0.00002289;
//                            errV22pPb[counter] = 0.0005914;
//                        }
//                        else if(TMath::Abs(sit[counter]-75) <0.1){
//                            V22pPb[counter] = -0.00001888;
//                            errV22pPb[counter] = 0.0004244;
//                        }
                        cout<< "ALERT"<<endl;
                    }
            counter++;
                }
           // cout << "Finished looking at "<<str<<endl;
            }
    }
    double v_moy = 0;
    double v_syst = 0;
    for(int idx=0; idx<counter;idx++){
        v_moy+=v2pPb[idx];
    }
    v_moy/=counter;
    
    for(int idx=0; idx<counter;idx++){
        v_syst+=pow(v2pPb[idx]-v_moy,2);
    }
    v_syst/=counter;
    v_syst=sqrt(v_syst);
    
    cout << "Moyenne " << v_moy<<endl;
    cout << "Syst " << v_syst<<endl;
    
    
    TGraphErrors *v2syst = new TGraphErrors(numsituations,sit,V22pPb,errsit,errV22pPb);
     v2syst->SetTitle("Centrality study PYTHIA - ATLAS - Periph: 5-0");
    TAxis *ax = v2syst->GetHistogram()->GetXaxis();
    Double_t x1 = ax->GetBinLowEdge(1);
    Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
    v2syst->GetHistogram()->GetXaxis()->Set(numsituations,x1-0.5,x2+0.5);
     v2syst->GetHistogram()->GetYaxis()->SetRangeUser(-0.005,0.005);
     v2syst->GetHistogram()->GetYaxis()->SetTitle("V_{22,tkl}");
    v2syst->SetMarkerStyle(20);
    v2syst->SetMarkerColor(kMagenta);
    v2syst->SetLineColor(kMagenta);
    
    TLine *l=new TLine(0.,0.0,140,0.0);
    l->SetLineColor(kBlack);
    l->SetLineWidth(1);
    l->SetLineStyle(9);

//    for(Int_t k=0;k<numsituations;k++){
//    v2syst->GetHistogram()->GetXaxis()->SetBinLabel(k+1,names[k].c_str());
//    }
    
    //v2syst->Draw("AP");
    //l->Draw("same");

    
    
    
    TMultiGraph* mg = new TMultiGraph();

    //loop on your files reading each TGraph and add it to the TMultiGraph
    TFile *filerec5 = new TFile("~/Desktop/5-0.root");
    TFile *filerec10 = new TFile("~/Desktop/10-0.root");
    TFile *filerec20 = new TFile("~/Desktop/20-10.root");
    TFile *filerec30 = new TFile("~/Desktop/30-20.root");
    TGraph* g = (TGraph*)filerec30->Get("AP");
    mg->Add(g);
    g = (TGraph*)filerec20->Get("AP");
    mg->Add(g);
    g = (TGraph*)filerec10->Get("AP");
    mg->Add(g);
    g = (TGraph*)filerec5->Get("AP");
    mg->Add(g);
    
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("N^{ch}");
    mg->GetYaxis()->SetTitle("V_{22,tkl}");
    mg->GetYaxis()->SetRangeUser(-0.006,0.01);
    l->Draw("same");

Char_t CanvasName[500];
Char_t FolderName[500];
sprintf(FolderName,"~/Desktop");
sprintf(CanvasName,"%s/PYTHIACheckFullerFinal.pdf",FolderName);
c1->SaveAs(CanvasName);
    
    
    TFile *filerec = new TFile("~/Desktop/Cumulative.root","RECREATE");
    mg->Write("Cumul");
    
    //hnseg = (TH1F*)filerec->Get("hnseg");
    
}





