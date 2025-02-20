#!/bin/bash

SCRIPTDIR=$(pwd)
DOSSIER="/Users/sperrin/Desktop/ImagesJavierAnalysis/2022fevrier/"

sh $(pwd)/GetListOfRadicals.sh

runAnalysisTKL() {
        echo "Will now run the analysis of ${1}"
        sh $(pwd)/RunRootTKL.sh $1 >& "${DOSSIER}${1}/log.txt"
	sleep 2
}

runAnalysisDimu() {
        echo "Will now run the analysis of ${1}"
        sh $(pwd)/RunRootDimu.sh $1 >& "${DOSSIER}${1}/log.txt"
        sleep 2
}

echo "HMM"

if [ -f $(pwd)/ListRadicals.txt ]
	then
	for file in "$(pwd)/ListRadicals.txt"
                do
                while IFS="" read -r line || [ -n "${line}" ]
                        do
			echo "Identified this analysis to run from list of radicals: ${line}"
			if [ -f "${DOSSIER}${line}/SistematicsFile.csv" ]
                        	then echo "Le dossier existe deja, l'analyse a sans doute déjà tourne"
			else
				mkdir "${DOSSIER}${line}"
				runAnalysisTKL "${line}"
			fi
                        done < "${pwd}ListRadicals.txt"
                done
fi

