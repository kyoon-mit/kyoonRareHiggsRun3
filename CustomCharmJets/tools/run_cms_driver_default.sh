nanoaod_dir="${CMSSW_BASE}/src/kyoonRareHiggsRun3/CustomCharmJets"
input_file="${nanoaod_dir}/input/step6_1.root"
output_file="${nanoaod_dir}/output/nanoaod_default_1.root"
output_branches="${nanoaod_dir}/output/branchlist_nanoaod_default_1.txt"

cmsDriver.py mc_default_2018UL \
    --python_filename ${nanoaod_dir}/python/step7_cfg.py \
    --mc \
    --eventcontent NANOAODSIM \
    --datatier NANOAODSIM \
    --filein file:${input_file} \
    --fileout file:${output_file} \
    --conditions 106X_upgrade2018_realistic_v16_L1v1 \
    --step NANO \
    --era Run2_2018,run2_nanoAOD_106Xv2 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    -n -1

run_root="root -l -b -q '${nanoaod_dir}/tools/getbranchlist.C(\"${output_file}\", \"Events\")' > ${output_branches}"
eval $run_root