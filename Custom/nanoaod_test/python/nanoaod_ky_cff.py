import FWCore.ParameterSet.Config as cms
from Custom.nanoaod_test.custom_jetCollectionTools import RecoJetAdder
from Custom.nanoaod_test.custom_jme_cff import ReclusterAK4CHSJets

def nanoAOD_customizeJets(process):
    process.load('Custom.nanoaod_test.MesonsReco_cff')
    process.load('Custom.nanoaod_test.DiMuonReco_cff')
    # process.finalJetsAK4Constituents = cms.EDProducer('PatJetConstituentPtrSelector',
    #                                                   src = cms.InputTag('slimmedJetsPuppi'),
    #                                                   cut = cms.string('')
    #                                                  )
    
    # process.customAK4ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
    #                                                 candidates = candInput,
    #                                                 jets = cms.InputTag("finalJetsPuppi"), # was finalJets before
    #                                                 jet_radius = cms.double(0.4),
    #                                                 name = cms.string("JetPFCands"),
    #                                                 idx_name = cms.string("pFCandsIdx"),
    #                                                 nameSV = cms.string("JetSVs"),
    #                                                 idx_nameSV = cms.string("sVIdx"),
    #                                                 )
    # process.customizedPFCandsTask.add(process.customAK4ConstituentsTable)
    # muonCollection = cms.InputTag('linkedObjects', 'muons')
    process.nanoSequence = cms.Sequence(process.nanoSequence + process.V0Sequence + process.V0Tables + process.DiMuProdSequence + process.DiMuTables)

    # Jet customization
    runOnMC = False
    recoJA = RecoJetAdder(runOnMC=runOnMC)
    process = ReclusterAK4CHSJets(process, recoJA, runOnMC)
    return process