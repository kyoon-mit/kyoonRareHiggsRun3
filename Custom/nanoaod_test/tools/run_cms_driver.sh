nanoaod_dir="${CMSSW_BASE}/src/Custom/nanoaod_test"
input_file="${nanoaod_dir}/input/step6_1.root"
output_file="${nanoaod_dir}/output/nanoaod_test.root"
output_branches="${nanoaod_dir}/output/branchlist_nanoaod_test.txt"
cmsDriver.py mc_2018UL \
    --mc \
    --eventcontent NANOAODSIM \
    --datatier NANOAODSIM \
    --fileout file:${output_file} \
    --conditions 106X_upgrade2018_realistic_v16_L1v1 \
    --step NANO \
    --filein file:${input_file} \
    --era Run2_2018,run2_nanoAOD_106Xv2 \
    --customise Custom/nanoaod_test/nanoaod_ky_cff.nanoAOD_customizeJets \
    -n -1

run_root="root -l -b -q '${nanoaod_dir}/output/getbranchlist.C(\"${output_file}\", \"Events\")' > ${output_branches}"
eval $run_root