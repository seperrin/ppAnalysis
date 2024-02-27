#!/bin/bash

SCRIPTDIR="${pwd}"

if [ -f "${pwd}ListRadicals.txt" ]
        then mv "${pwd}ListRadicals.txt" "${pwd}OldListRadicals.txt"
fi

touch "${pwd}ListRadicals.txt"

declare -a ZvtxCutArray=("8" "10" "12")
declare -a DPhiCutArray=("5mrad")
declare -a EtaGapArray=("1.0" "1.2" "1.4" "1.6" "1.8" "2.0")
declare -a CentralityEstimatorArray=("V0MPercentile")
declare -a CentralClassArray=("0-3" "0-5" "0-10")
declare -a PeriphClassArray=("30-100" "40-100" "40-80" "50-100")
declare -a EMNormArray=("Method1" "Method2" "Method3" "Method4")
declare -a EMMaxArray=("50" "100" "500")
declare -a EMThresholdArray=("5" "10" "50")
declare -a PeriphScalingArray=("0.9" "1" "1.2")
declare -a SummationZvtxArray=("Method1" "Method2")

DefaultZvtxCut="10"
DefaultDPhiCut="10mrad"
DefaultEtaGap="1.2"
DefaultCentralClass="0-5"
DefaultPeriphClass="40-100"
DefaultEMNorm="Method1"
DefaultEMMax="30"
DefaultEMThreshold="10"
DefaultPeriphScaling="1"
DefaultSummationZvtx="Method1c"






DefaultDataUsed="Run2"
declare -a SystematicsFocusArray=("")


Radical="NewAnalysisAllEstCentBvr_"
Radical="${Radical}${DefaultDataUsed}"



echo ${Radical} >> "${pwd}ListRadicals.txt"

contains() {
    echo $1 | grep -w -q $2
}

addParamMaybe(){
	if [[ "${SystematicsFocusArray[*]}" =~ "$1" ]]
		then echo "Trouve dans la liste des focus à faire"
		addParam "$@"
	fi
}
 
addParam() {
	if [[ "${SystematicsFocusArray[*]}" =~ "$1" ]]
		then addVariable "$2" "$3" "${@:4}"
	else
		addDefault "$2" "$3"
	fi
}

addDefault() {		
	mv "${pwd}ListRadicals.txt" "${pwd}OldListRadicals.txt"
	touch "${pwd}ListRadicals.txt"
	for file in "${pwd}OldListRadicals.txt"
		do
		while IFS="" read -r line || [ -n "${line}" ]
			do
  			echo "${line}_$1$2" >> "${pwd}ListRadicals.txt"
			done < "${pwd}OldListRadicals.txt"
		done
}

addVariable() {
	mv "${pwd}ListRadicals.txt" "${pwd}OldListRadicals.txt"
	touch "${pwd}ListRadicals.txt"
	for variable in "${@:3}"
		do 
		if [ "${variable}" != "$2" ]
			then
			for file in "${pwd}OldListRadicals.txt"
                		do
                		while IFS="" read -r line || [ -n "${line}" ]
                        		do
                        		echo "${line}_$1${variable}" >> "${pwd}ListRadicals.txt"
                        		done < "${pwd}OldListRadicals.txt"
                		done
		fi
		done
}

addVariable "" "" "${CentralityEstimatorArray[@]}"
addParam "CentralClass" "" ${DefaultCentralClass} "${CentralClassArray[@]}"
addParam "PeriphClass" "" ${DefaultPeriphClass} "${PeriphClassArray[@]}"
addDefault "pt" "0-3-6"

addParamMaybe "ZvtxCut" "Zvtx" ${DefaultZvtxCut} "${ZvtxCutArray[@]}"
addParamMaybe "DPhiCut" "DPhiCut" ${DefaultDPhiCut} "${DPhiCutArray[@]}"
addParamMaybe "EtaGap" "Eta" ${DefaultEtaGap} "${EtaGapArray[@]}"
addParamMaybe "EMNorm" "EMNorm" ${DefaultEMNorm} "${EMNormArray[@]}"
addParamMaybe "EMMax" "EMMax" ${DefaultEMMax} "${EMMaxArray[@]}"
addParamMaybe "EMThreshold" "EMThreshold" ${DefaultEMThreshold} "${EMThresholdArray[@]}"
addParamMaybe "PeriphScaling" "PeriphScaling" ${DefaultPeriphScaling} "${PeriphScalingArray[@]}"
addParamMaybe "SummationZvtx" "SummationZvtx" ${DefaultSummationZvtx} "${SummationZvtxArray[@]}"

