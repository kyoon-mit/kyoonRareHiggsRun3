<!-- README -->
# Rare Higgs Decays at the CMS Detector using Run 2+3 Data
> Author: Kyungseop Yoon

### CMSSW version: CMSSW_14_1_0_pre4

This version is intended to be used with [Combine v10.0.2](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#combine-v10-recommended-version).

### Setup
1. Obtain CMSSW_14_1_0 and build Combine following the [instruction](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#combine-v10-recommended-version).
2. Run `. env.sh`.

## Analyses
> 1. $H \rightarrow J/\Psi + c\bar{c}$

    branch: main
    .
    |-- JPsiCC        # H -> ccbar + J/Psi
        |-- analysis
        |-- housekeeping
        |-- interface
        |-- json
        |-- src
        |-- studies
            |-- genstudies
    |-- python-package
