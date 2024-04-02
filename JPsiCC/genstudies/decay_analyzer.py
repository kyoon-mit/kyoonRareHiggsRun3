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

# Find Higgs
rdf = (rdf.Define('GenPart_Higgs_idx',
                'HiggsIdx(GenPart_pdgId, GenPart_genPartIdxMother)')
        .Define('Higgs_energy',
                'GenPart_energy[GenPart_Higgs_idx]')
        .Define('Higgs_eta',
                'GenPart_eta[GenPart_Higgs_idx]')
        .Define('Higgs_mass',
                'GenPart_mass[GenPart_Higgs_idx]')
        .Define('Higgs_phi',
                'GenPart_phi[GenPart_Higgs_idx]')
        .Define('Higgs_pt',
                'GenPart_pt[GenPart_Higgs_idx]')
        .Define('Higgs_px',
                'GenPart_px[GenPart_Higgs_idx]')
        .Define('Higgs_py',
                'GenPart_py[GenPart_Higgs_idx]')
        .Define('Higgs_pz',
                'GenPart_pz[GenPart_Higgs_idx]')
      )

# Find Higgs Daughter
rdf = (rdf.Define('GenPart_HiggsDaughters_pdgId',
                'HiggsDaughtersPDG(GenPart_pdgId, GenPart_genPartIdxMother)')
        .Define('GenPart_HiggsDaughters_idx',
                'HiggsDaughtersIdx(GenPart_pdgId, GenPart_genPartIdxMother)')
        .Define('GenPart_HiggsGrandDaughters_pdgId',
                'GenericDaughtersPDG(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_HiggsDaughters_idx)')
        .Define('GenPart_HiggsGrandDaughters_idx',
                'GenericDaughtersIdx(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_energy',
                'SelectByIdx(GenPart_energy, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_eta',
                'SelectByIdx(GenPart_eta, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_mass',
                'SelectByIdx(GenPart_mass, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_phi',
                'SelectByIdx(GenPart_phi, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_pt',
                'SelectByIdx(GenPart_pt, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_px',
                'SelectByIdx(GenPart_energy, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_py',
                'SelectByIdx(GenPart_py, GenPart_HiggsDaughters_idx)')
        .Define('HiggsDaughters_pz',
                'SelectByIdx(GenPart_pz, GenPart_HiggsDaughters_idx)')
      )

# Snapshot
branches_of_interest =\
    ['Higgs_energy',
     'Higgs_eta',
     'Higgs_mass',
     'Higgs_phi',
     'Higgs_pt',
     'Higgs_px',
     'Higgs_py',
     'Higgs_pz',
     'GenPart_HiggsDaughters_pdgId',
     'GenPart_HiggsGrandDaughters_pdgId',
     'HiggsDaughters_energy',
     'HiggsDaughters_eta',
     'HiggsDaughters_mass',
     'HiggsDaughters_phi',
     'HiggsDaughters_pt',
     'HiggsDaughters_px',
     'HiggsDaughters_py',
     'HiggsDaughters_pz']

branches_of_second_interest =\
    ['GenPart_Higgs_idx',
     'GenPart_HiggsDaughters_idx',
     'GenPart_HiggsGrandDaughters_idx']

rdf.Snapshot('GenEvents', 'genevents.root', (branches_of_interest + branches_of_second_interest))

# Get some histograms

my_hists = []
for var in branches_of_interest:
    print(var)
    my_hists.append(rdf.Histo1D(var))

# Draw
save_dir = os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'plots', 'genstudies')
for hist in my_hists:
    c = ROOT.TCanvas()
    hist.Draw()
    c.SaveAs(os.path.join(save_dir, f'{hist.GetTitle()}.png'))
    c.Close()