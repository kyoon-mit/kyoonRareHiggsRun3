# From https://dalfonso.web.cern.ch/step7_RUN0_cfg.py
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --python_filename step7_cfg.py --eventcontent NANOAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --fileout file:step7.root --conditions 106X_upgrade2018_realistic_v16_L1v1 --step NANO --filein file:step6.root --era Run2_2018,run2_nanoAOD_106Xv2 --no_exec --mc -n -1
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Modifier_run2_nanoAOD_106Xv2_cff import run2_nanoAOD_106Xv2

process = cms.Process('NANO',Run2_2018,run2_nanoAOD_106Xv2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        '/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hphigamma-powheg/UL2018-MINIAODSIMv9/220102_170607/0000/step6_0.root'
#        '/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hrhogamma-powheg/UL2018-MINIAODSIMv9/220102_170758/0000/step6_0.root'
#        '/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hphiKLKSgamma-powheg/UL2018-MINIAODSIMv9/220207_215847/0000/step6_0.root' 
#        '/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/ggh-homegagamma-powheg/UL2018-MINIAODSIMv9/230123_193725/0001/step6_0.root' 
#        '/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/ggh-homegagamma-powheg-fixgen/UL2018-MINIAODSIMv9/231026_212158/0000/step6_0.root' 
#        'file:/afs/cern.ch/work/s/selvaggi/public/4Maria/hjpsicc/step6.root'
#        '/store/cmst3/group/vhcc/jp/samples/GluGluH_HJPsiCC_M125_13TeV_powheg_pythia8_MINIAOD/UL2018_test/240312_155903/0000/step6_0.root'
#        '/store/cmst3/group/vhcc/jp/samples/GluGluH_HJPsiCC_MUMUCC_KC1p0_M125_13TeV_powheg_pythia8_MINIAOD_/UL2018_V1/240720_142323/0003/step6_0.root'
        '/store/cmst3/group/vhcc/jp/samples/ZJets_ZJPsiCC_13TeV_amcatnlo_pythia8_MINIAOD/UL2018_test/240312_155651/0000/step6_0.root'
        
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('--python_filename nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
#    fileName = cms.untracked.string('root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/vbf-hphigamma-powheg/NANOAOD_00/step7_VBS_Phigamma_0.root'),
#    fileName = cms.untracked.string('root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hphigamma-powheg/NANOAOD_01/step7_VBS_Phigamma_0.root'),   
#    fileName = cms.untracked.string('root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hrhogamma-powheg/NANOAOD_01/step7_VBS_Rhogamma_0.root'),   
#    fileName = cms.untracked.string('root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hphiKLKSgamma-powheg/NANOAOD_01/step7_VBS_PhiKLKSgamma_0.root'
#    fileName = cms.untracked.string('root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/ggh-homegagamma-powheg-fixgen/NANOAOD_03_test6/step7_ggH_OmegaGamma_0.root'
#    fileName = cms.untracked.string('/tmp/dalfonso/step7.root'
#    fileName = cms.untracked.string('root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/JpsiTest/GluGluH_HJPsiCC/NANOAOD_test3/step7_GluGluH_HJPsiCC_0.root'
    fileName = cms.untracked.string('root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/JpsiTest/ZJets_ZJPsiCC/NANOAOD_D04/nano_ZJets_ZJPsiCC_0.root'
    ),   
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# Automatic addition of the customisation function from Hrare.NanoAOD.nano_cff

from Hrare.NanoAOD.nano_cff import nanoAOD_customizeMesons

#call to customisation function nanoAOD_customizeMesons imported from Hrare.NanoAOD.nano_cff
process = nanoAOD_customizeMesons(process)

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion