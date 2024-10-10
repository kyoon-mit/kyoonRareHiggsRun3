#!/bin/bash

export HRARE_DIR=$(dirname $(realpath -s $BASH_SOURCE[0]))
if [[ -z "${CMSSW_BASE}" ]]; then
    echo "Please run cmsenv in the top directory first."
else
    cd ${CMSSW_BASE}/src/kyoonRareHiggsRun3/CustomCharmJets
    scram b -j 4
fi