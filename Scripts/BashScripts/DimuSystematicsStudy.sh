#!/bin/bash


SCRIPTDIR=${pwd}

if [ -f "${pwd}ListCSV.txt" ]
        then mv "${pwd}ListCSV.txt" "${pwd}OldListCSV.txt"
fi

for fo in /Users/sperrin/Desktop/ImagesJavierAnalysis/2022juinDimu/NewAnalysisAllEstCentBvr_Run2_V0M*/*.csv
        do
        focsv="${fo}"
        if [ -f "${focsv}" ]
                then echo "${focsv}" >> ListCSV.txt
        fi
        done
root -l <<-EOF
.L ${SCRIPTDIR}../AliAnalysisTaskMyMuonTree_AOD.cxx++
.L ${SCRIPTDIR}../../SystematicsDimu.C

SystematicsDimu("${pwd}ListCSV.txt","ZvtxCut")
SystematicsDimu("${pwd}ListCSV.txt","EtaMin")
SystematicsDimu("${pwd}ListCSV.txt","EtaMax")
SystematicsDimu("${pwd}ListCSV.txt","SummationZvtx")
SystematicsDimu("${pwd}ListCSV.txt","Pooling")
SystematicsDimu("${pwd}ListCSV.txt","BackgroundV2")
SystematicsDimu("${pwd}ListCSV.txt","RangeV2")
SystematicsDimu("${pwd}ListCSV.txt","ExtractionMethod")

sleep(10000)
EOF
read -p "Want to close ?  (Y/n) " -n 1 -r
