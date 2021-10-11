#!/bin/bash

FOLDER="/Users/sperrin/Desktop/ImagesJavierAnalysis/2021octobre/*/"
SCRIPTDIR="${pwd}"
prefix="/Users/sperrin/Desktop/ImagesJavierAnalysis/2021octobre/"
suffix="/"
for fo in /Users/sperrin/Desktop/ImagesJavierAnalysis/2021octobre/NewAnalysisAllEst_TKL_16h25_*/ 
	do
	echo "Looking into folder ${fo}"
	focsv="${fo}SystematicsFile.csv"
	echo ${focsv}
	if [ -f "${focsv}" ] 
		then echo "There is already a .csv file in here"
	else
		echo "Will run GenerateCSV on ${fo}" 
		read -p "Are you sure? (Y/n) " -n 1 -r
		echo    # (optional) move to a new line
		if [[ ! $REPLY =~ ^[Yy]$ ]]
		then
    			continue
		fi
		foo=${fo#"$prefix"}
		foo=${foo%"$suffix"}
		echo ${foo}
		root -b -l <<-EOF
		.L ${SCRIPTDIR}../AliAnalysisTaskMyMuonTree_AOD.cxx++
		.L ${SCRIPTDIR}../../FitTrainingTKL.C
		FitTrainingTKL("${foo}")
		EOF
	fi
	done
