#!/bin/bash

root -l <<-EOF
.L ${SCRIPTDIR}../AliAnalysisTaskMyMuonTree_AOD.cxx++
.L ${SCRIPTDIR}../../FitTrainingPtBinned.C
Runner("${1}")
EOF
