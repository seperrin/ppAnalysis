#!/bin/bash

root -l <<-EOF
.L ${SCRIPTDIR}../AliAnalysisTaskMyMuonTree_AOD.cxx++
.L ${SCRIPTDIR}../../PlotFromTreeTKL.C
PlotFromTreeTKL("${1}")

EOF
