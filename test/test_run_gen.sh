# From https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=GluGluHToBB_M-125_TuneCP5_13p6TeV_powheg-pythia8&page=0&shown=127

# Adapted from TSG-Run3Winter24wmLHEGS-00080
# Detailed command -> https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_setup/TSG-Run3Winter24wmLHEGS-00080
cmsDriver.py fragments/fragment_ggh-hphigamma.py \
             --python_filename test_out_fragment.py \
             --eventcontent RAWSIM,LHE \
             --datatier GEN-SIM,LHC \
             --fileout test.root \
             --conditions 133X_mcRun3_2024_realistic_v7 \
             --beamspot Realistic25ns13p6TeVEarly2023Collision \
             --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed="int(125125)" \
             --step LHE,GEN,SIM \
             --geometry DB:Extended \
             --era Run3_2023 \
             --mc \
             -n 1000 \
> log.txt