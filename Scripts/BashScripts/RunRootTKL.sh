#!/bin/bash

root -l <<-EOF
.L ${SCRIPTDIR}../AliAnalysisTaskMyMuonTree_AOD.cxx++
.L ${SCRIPTDIR}../../FitTrainingTKL.C
FitTrainingTKL("${1}")
EOF
