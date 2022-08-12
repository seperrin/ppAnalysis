#!/bin/bash


SCRIPTDIR=${pwd}

if [ -f "${pwd}ListCSV.txt" ]
        then mv "${pwd}ListCSV.txt" "${pwd}OldListCSV.txt"
fi

for fo in /Users/sperrin/Desktop/ImagesJavierAnalysis/2022aoutPYTHIA/NewAnalysisAllEstMCATLASCentBvrNospe_TKL_16h_*/
        do
        focsv="${fo}SystematicsFile.csv"
        if [ -f "${focsv}" ]
                then echo "${focsv}" >> ListCSV.txt
        fi
        done
root -l <<-EOF
.L ${SCRIPTDIR}../AliAnalysisTaskMyMuonTree_AOD.cxx++
.L ${SCRIPTDIR}../../SystematicsTKL.C
CentralityStudyTKLPYTHIA("${pwd}ListCSV.txt","")
sleep(1000)
EOF
read -p "Want to close ?  (Y/n) " -n 1 -r
