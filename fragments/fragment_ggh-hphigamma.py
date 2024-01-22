import FWCore.ParameterSet.Config as cms
# Link to cards
externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/RunIII/13p6TeV/slc7_amd64_gcc10/Powheg/V2/gg_H_quark-mass-effects_slc7_amd64_gcc10_CMSSW_12_4_6_my_ggH_M125.tgz'),
    nEvents = cms.untracked.uint32(500),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh'),
    generateConcurrently = cms.untracked.bool(False)
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import * ## Found here: https://github.com/cms-sw/cmssw/blob/f5cbbf1ab5806ffc54ccaab6a77201f35d63fca7/Configuration/Generator/python/MCTunesRun3ECM13p6TeV/PythiaCP5TuneUpSettings_cfi.py
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

generator = cms.EDFilter("Pythia8ConcurrentHadronizerFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13600.),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PSweightsSettingsBlock,
        pythia8PowhegEmissionVetoSettingsBlock,
        processParameters = cms.vstring(
            'JetMatching:setMad = off',
            'JetMatching:scheme = 1',
            'JetMatching:merge = on',
            'JetMatching:jetAlgorithm = 2',
            'JetMatching:etaJetMax = 999.',
            'JetMatching:coneRadius = 1.',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = 30.', # this is the actual merging scale
            'JetMatching:doFxFx = on',
            'JetMatching:qCutME = 10.',# this must match the ptj cut in the lhe generation step
            'JetMatching:nQmatch = 5', # 4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
            'JetMatching:nJetMax = 2', # number of partons in born matrix element for highest multiplicity
            'SLHA:useDecayTable = off',
            'POWHEG:nFinal = 1',   ## Number of final state particles
                                   ## (BEFORE THE DECAYS) in the LHE
                                   ## other than emitted extra parton
            '25:m0 = 125.0', ## higgs mass
            '25:onMode = off', ## switch off all decay channels
            '25:addChannel = 1 1.00 103 22 333' ## 22=photon, 333=phi(1020) \\ not sure what the others are
            '25:onIfMatch = 321 -321' ## K+ K-
        ),
    parameterSets = cms.vstring('pythia8CommonSettings',
                                'pythia8CP5Settings',
                                'pythia8PSweightsSettings',
                                'pythia8aMCatNLOSettings',
                                'pythia8PowhegEmissionVetoSettings',
                                'processParameters',
                                )
    )
)

ProductionFilterSequence = cms.Sequence(generator)