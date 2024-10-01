#!/bin/bash

export HRARE_DIR=$(dirname $(realpath -s $BASH_SOURCE[0]))
if [[ -z "${CMSSW_BASE}" ]]; then
    echo "Please run cmsenv in the top directory first."
else
    cp -ru ${HRARE_DIR}/Custom ${CMSSW_BASE}/src
    cd ${CMSSW_BASE}/src/Custom
    scram b -j 4
fi