my_test_input_file=step6_1.root
cmsDriver.py mc_2018UL \
    --mc \
    --eventcontent NANOAODSIM \
    --datatier NANOAODSIM \
    --fileout file:nanoaod_test.root \
    --conditions 106X_upgrade2018_realistic_v16_L1v1 \
    --step NANO \
    --filein file:${my_test_input_file} \
    --era Run2_2018,run2_nanoAOD_106Xv2 \
    --customise Custom/nanoaod_test/nanoaod_ky_cff.nanoAOD_customizeJets \
    -n -1

root -l -b -q 'output/getbranchlist.C("nanoaod_test.root", "Events")' > output/branchlist.txt