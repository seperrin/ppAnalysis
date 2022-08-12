#!/bin/bash


SCRIPTDIR=${pwd}

if [ -f "${pwd}ListCSV.txt" ]
	then mv "${pwd}ListCSV.txt" "${pwd}OldListCSV.txt"
fi

for fo in /Users/sperrin/Desktop/ImagesJavierAnalysis/2022juinMixed/NewAnalysisAllEstMCNO_TKL_16h_*/
	do
	focsv="${fo}SystematicsFile.csv"
	if [ -f "${focsv}" ]
		then echo "${focsv}" >> ListCSV.txt
	fi
	done
root -l <<-EOF
.L ${SCRIPTDIR}../AliAnalysisTaskMyMuonTree_AOD.cxx++
.L ${SCRIPTDIR}../../SystematicsTKL.C
SystematicsTKL("${pwd}ListCSV.txt","ZvtxCut")
SystematicsTKL("${pwd}ListCSV.txt","EtaGap")
SystematicsTKL("${pwd}ListCSV.txt","SummationZvtx")
SystematicsTKL("${pwd}ListCSV.txt","FitMethod")
EOF
read -p "Want to close ?  (Y/n) " -n 1 -r 
