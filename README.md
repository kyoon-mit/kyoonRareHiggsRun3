<!-- README -->
# Rare Higgs Decays at the CMS Detector using Run 2+3 Data
> Author: Kyungseop Yoon

### CMSSW version: CMSSW_13_3_0

### Setup
1. Open `hrare.yml`. Scroll to the bottom and modify `prefix` to reflect where you want the environment to be installed. Create the conda environment `hrare` from the `.yml` file.
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
