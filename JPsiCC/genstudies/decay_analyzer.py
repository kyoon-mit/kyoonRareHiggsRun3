'''
This module examines the decay of a particle using the gen-level mother and daughter information.
'''

from kytools import jsonreader
import os
import ROOT

# Load macro
ROOT.gSystem.CompileMacro(os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'src', 'GenAnalyzer.cc'), 'k')

# Disable multithreading for now
ROOT.DisableImplicitMT()

# Get TChain of generator-level events in the NANOAOD format
genevents = jsonreader.get_chain_from_json(anpath='JPsiCC',
                                           jsonname='NANOAOD.json',
                                           keys=['mariadlf', '20240321', 'GEN-signal'],
                                           treename='Events')

# Create RDF from TChain
rdf = ROOT.RDataFrame(genevents)

# Find Higgs Daughter
rdf = rdf.Define('GenPart_HiggsDaughters_pdgId',
                 'HiggsDaughtersPDG(GenPart_pdgId, GenPart_genPartIdxMother)')\
         .Define('GenPart_HiggsDaughters_idx',
                 'HiggsDaughtersIdx(GenPart_pdgId, GenPart_genPartIdxMother)')\
         .Define('GenPart_HiggsGrandDaughters_pdgId',
                 'GenericDaughtersPDG(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_HiggsDaughters_idx)')\
         .Define('GenPart_HiggsGrandDaughters_idx',
                 'GenericDaughtersIdx(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_HiggsDaughters_idx)')

rdf.Snapshot('GenEvents', 'genevents.root', ['GenPart_HiggsDaughters_pdgId',
                                             'GenPart_HiggsDaughters_idx',
                                             'GenPart_HiggsGrandDaughters_pdgId',
                                             'GenPart_HiggsGrandDaughters_idx'])

# Get some histograms
# vars_of_interest = ['GenPart_HiggsDaughters_pdgId']

# my_hists = []
# for var in vars_of_interest:
#     print(var)
#     my_hists.append(rdf.Histo1D(var))

# Draw
# save_dir = os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'plots', 'genstudies')
# for hist in my_hists:
#     c = ROOT.TCanvas()
#     hist.Draw()
#     c.SaveAs(os.path.join(save_dir, f'{hist.GetTitle()}.png'))
#     c.Close()