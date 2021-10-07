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

void SystematicsTKL(Char_t listeCSVFiles[100], Char_t SystematicsFocus[50])
{
    std::ifstream file(listeCSVFiles);
    CSVRow row;
    std::string str;
    
    double v2pPbPMSPDT[10]={0};
    double errv2pPbPMSPDT[10]={0};
    double v2pPbSPDT[10]={0};
    double errv2pPbSPDT[10]={0};
    double v2pPbSPDC[10]={0};
    double errv2pPbSPDC[10]={0};
    double v2pPbV0M[10]={0};
    double errv2pPbV0M[10]={0};
    
    double v2ATLASPMSPDT[10]={0};
    double errv2ATLASPMSPDT[10]={0};
    double v2ATLASSPDT[10]={0};
    double errv2ATLASSPDT[10]={0};
    double v2ATLASSPDC[10]={0};
    double errv2ATLASSPDC[10]={0};
    double v2ATLASV0M[10]={0};
    double errv2ATLASV0M[10]={0};
    
    int IndexOfVariable;
    int IndexOfStorage;
    
    string variableDPhiCut[5] = {"1mrad", "2mrad", "5mrad", "10mrad", "None"};
    string variableZvtxCut[3] = {"8","10","12"};
    string variableEtaGap[4] = {"1.0","1.2","1.6","2.0"};
    
    
    if (strcmp(SystematicsFocus, "DPhiCut") == 0){
        IndexOfVariable = 0;
    }
    if(strcmp(SystematicsFocus, "ZvtxCut") == 0){
        IndexOfVariable = 5;
    }
    if(strcmp(SystematicsFocus, "EtaGap") == 0){
        IndexOfVariable = 4;
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
                        for(int idx=0; idx<4; idx++){
                            if (variableEtaGap[idx] == row[IndexOfVariable]){
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
                }
                else if(isRowSuitableForSystematics(row, SystematicsFocus, "PercentileMethodSPDTracklets")){
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
                        for(int idx=0; idx<4; idx++){
                            if (variableEtaGap[idx] == row[IndexOfVariable]){
                                IndexOfStorage = idx;
                                isVariationFound = kTRUE;
                                continue;
                            }
                        }
                    }
                    
                    if(isVariationFound){
                        v2pPbPMSPDT[IndexOfStorage] = stod(std::string(row[11]));
                        errv2pPbPMSPDT[IndexOfStorage] = stod(std::string(row[12]));
                        v2ATLASPMSPDT[IndexOfStorage] = stod(std::string(row[19]));
                        errv2ATLASPMSPDT[IndexOfStorage] = stod(std::string(row[20]));
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
                        for(int idx=0; idx<4; idx++){
                            if (variableEtaGap[idx] == row[IndexOfVariable]){
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
                        for(int idx=0; idx<4; idx++){
                            if (variableEtaGap[idx] == row[IndexOfVariable]){
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
                }
                else{
                    continue;
                }
            }
        
    }
    
    cout << "\n\n\n ==================== \n Systematics results \n ==================== \n\n\n"<<endl;
    
    cout << "Centrality estimator : PercentileMethodSPDTracklets"<<endl;
    
    cout << "p-Pb Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2pPbPMSPDT[idx]!=0 || errv2pPbPMSPDT[idx]!=0){
        cout << "v2: "<<v2pPbPMSPDT[idx]<< " +/- "<<errv2pPbPMSPDT[idx] <<endl;
            nbMesures++;
            v2_moy+=v2pPbPMSPDT[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2pPbPMSPDT[idx]!=0 || errv2pPbPMSPDT[idx]!=0){
            v2_syst+=pow(v2_moy-v2pPbPMSPDT[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    //Plot et calculer la systematique
    }
    
    cout << "ATLAS Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2ATLASPMSPDT[idx]!=0 || errv2ATLASPMSPDT[idx]!=0){
        cout << "v2: "<<v2ATLASPMSPDT[idx]<< " +/- "<<errv2ATLASPMSPDT[idx] <<endl;
            nbMesures++;
            v2_moy+=v2ATLASPMSPDT[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2ATLASPMSPDT[idx]!=0 || errv2ATLASPMSPDT[idx]!=0){
            v2_syst+=pow(v2_moy-v2ATLASPMSPDT[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    }
    
    
    
    
    cout << "Centrality estimator : SPDTrackletsPercentile"<<endl;
    
    cout << "p-Pb Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2pPbSPDT[idx]!=0 || errv2pPbSPDT[idx]!=0){
        cout << "v2: "<<v2pPbSPDT[idx]<< " +/- "<<errv2pPbSPDT[idx] <<endl;
            nbMesures++;
            v2_moy+=v2pPbSPDT[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2pPbSPDT[idx]!=0 || errv2pPbSPDT[idx]!=0){
            v2_syst+=pow(v2_moy-v2pPbSPDT[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    //Plot et calculer la systematique
    }
    
    cout << "ATLAS Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2ATLASSPDT[idx]!=0 || errv2ATLASSPDT[idx]!=0){
        cout << "v2: "<<v2ATLASSPDT[idx]<< " +/- "<<errv2ATLASSPDT[idx] <<endl;
            nbMesures++;
            v2_moy+=v2ATLASSPDT[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2ATLASSPDT[idx]!=0 || errv2ATLASSPDT[idx]!=0){
            v2_syst+=pow(v2_moy-v2ATLASSPDT[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    }
    
    
    cout << "Centrality estimator : SPDClustersPercentile"<<endl;
    
    cout << "p-Pb Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2pPbSPDC[idx]!=0 || errv2pPbSPDC[idx]!=0){
        cout << "v2: "<<v2pPbSPDC[idx]<< " +/- "<<errv2pPbSPDC[idx] <<endl;
            nbMesures++;
            v2_moy+=v2pPbSPDC[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2pPbSPDC[idx]!=0 || errv2pPbSPDC[idx]!=0){
            v2_syst+=pow(v2_moy-v2pPbSPDC[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    //Plot et calculer la systematique
    }
    
    cout << "ATLAS Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2ATLASSPDC[idx]!=0 || errv2ATLASSPDC[idx]!=0){
        cout << "v2: "<<v2ATLASSPDC[idx]<< " +/- "<<errv2ATLASSPDC[idx] <<endl;
            nbMesures++;
            v2_moy+=v2ATLASSPDC[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2ATLASSPDC[idx]!=0 || errv2ATLASSPDC[idx]!=0){
            v2_syst+=pow(v2_moy-v2ATLASSPDC[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    }
    
    
    cout << "Centrality estimator : V0MPercentile"<<endl;
    
    cout << "p-Pb Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2pPbV0M[idx]!=0 || errv2pPbV0M[idx]!=0){
        cout << "v2: "<<v2pPbV0M[idx]<< " +/- "<<errv2pPbV0M[idx] <<endl;
            nbMesures++;
            v2_moy+=v2pPbV0M[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2pPbV0M[idx]!=0 || errv2pPbV0M[idx]!=0){
            v2_syst+=pow(v2_moy-v2pPbV0M[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    //Plot et calculer la systematique
    }
    
    cout << "ATLAS Method"<<endl;
    {
    double v2_moy = 0;
    double v2_syst = 0;
    int nbMesures = 0;
    
    for(int idx=0; idx<10; idx++){
        if(v2ATLASV0M[idx]!=0 || errv2ATLASV0M[idx]!=0){
        cout << "v2: "<<v2ATLASV0M[idx]<< " +/- "<<errv2ATLASV0M[idx] <<endl;
            nbMesures++;
            v2_moy+=v2ATLASV0M[idx];
        }
    }
    v2_moy/=nbMesures;
    for(int idx=0; idx<10; idx++){
        if(v2ATLASV0M[idx]!=0 || errv2ATLASV0M[idx]!=0){
            v2_syst+=pow(v2_moy-v2ATLASV0M[idx],2);
        }
    }
    v2_syst = sqrt(v2_syst);
    cout << "Systematics: v2_moy = " <<v2_moy<<" systematic = "<<v2_syst<<" ("<< (v2_syst/v2_moy)*100<<"%)"<<endl;
    }
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
    Char_t EMMax[10] = "100"; //7
    Char_t EMThreshold[10] = "10"; //8
    Char_t PeriphScale[10] = "1"; //9
    Char_t SummationMethod[10] = "Method2"; //10
    
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
        if(row[index] != DefaultConditions[index]){
          //  cout << "row[" << index << "] = " << row[index] << " is not the same as " << DefaultConditions[index]<< " ABORT" <<endl;
            isSuitable = kFALSE;
            return isSuitable;
        }
    }
    
    return isSuitable;
    
}
