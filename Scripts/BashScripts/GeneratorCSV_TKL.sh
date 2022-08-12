#!/bin/bash

FOLDER="/Users/sperrin/Desktop/ImagesJavierAnalysis/2022aoutPYTHIA/*/"
SCRIPTDIR="${pwd}"
prefix="/Users/sperrin/Desktop/ImagesJavierAnalysis/2022aoutPYTHIA/"
suffix="/"
for fo in /Users/sperrin/Desktop/ImagesJavierAnalysis/2022aoutPYTHIA/NewAnalysisAllEstMCATLASCentBvrNospe_TKL_16h_*/ 
	do
	echo "Looking into folder ${fo}"
	focsv="${fo}SystematicsFile.csv"
	echo ${focsv}
	if [ -f "${focsv}" ] 
		then echo "There is already a .csv file in here"
	else
		echo "Will run GenerateCSV on ${fo}" 
		sleep 5
		echo    # (optional) move to a new line
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
