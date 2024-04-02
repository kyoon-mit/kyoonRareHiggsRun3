'''
This module analyzes various physics properties of the Higgs decay products from
a generator-level perspective.

The input file is produced by decay_analyzer.py. It is assumed that the file
contains the following information at the generator level
    - PDG ID of the Higgs daughters.
    - Energy, mass, momentum, and angle of the Higgs.
    - Energy, mass, momentum, and angle of the Higgs daughters.
'''
import sys, os
import ROOT

# Get input file
if len(sys.argv) == 1: input_file_name = 'genevents.root'
else: input_file_name = sys.argv[1]
input_file = ROOT.TFile(input_file_name, 'READ')

# Load macro
ROOT.gSystem.CompileMacro(os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'src', 'GenAnalyzer.cc'), 'k')

# Disable multithreading for now
ROOT.DisableImplicitMT()

# Create RDF
rdf = ROOT.RDataFrame('GenEvents', input_file)

# Get muons
rdf = (rdf.Define('muon1_index', 'GenPart_HiggsDaughters_pdgId==13')
          .Define('muon2_index', 'GenPart_HiggsDaughters_pdgId==-13')
          .Define('muon1_mass', 'HiggsDaughters_mass[muon1_index][0]')
          .Define('muon1_energy', 'HiggsDaughters_energy[muon1_index][0]')
          .Define('muon1_px', 'HiggsDaughters_px[muon1_index][0]')
          .Define('muon1_py', 'HiggsDaughters_py[muon1_index][0]')
          .Define('muon1_pz', 'HiggsDaughters_pz[muon1_index][0]')
          .Define('muon2_mass', 'HiggsDaughters_mass[muon2_index][0]')
          .Define('muon2_energy', 'HiggsDaughters_energy[muon2_index][0]')
          .Define('muon2_px', 'HiggsDaughters_px[muon2_index][0]')
          .Define('muon2_py', 'HiggsDaughters_py[muon2_index][0]')
          .Define('muon2_pz', 'HiggsDaughters_pz[muon2_index][0]')
      )

# Get J/Psi
rdf = (rdf.Define('JPsi_cand1',
                  'SumPxPyPzE(muon1_px, muon1_py, muon1_pz, muon1_energy,\
                              muon2_px, muon2_py, muon2_pz, muon2_energy)')
          .Define('JPsi_cand2',
                  'SumPxPyPzM(muon1_px, muon1_py, muon1_pz, muon1_mass,\
                              muon2_px, muon2_py, muon2_pz, muon2_mass)')
          .Define('JPsi_cand1_mass', 'JPsi_cand1.M()')
          .Define('JPsi_cand2_mass', 'JPsi_cand2.M()')
      )

#
# rdf = (rdf.Define('DiffVector',
#                   'DiffVector()'))

# Snapshot
branches_of_interest =\
    ['muon1_mass',
     'JPsi_cand1_mass',
     'JPsi_cand2_mass',
    ]

rdf.Snapshot('GenEvents', 'daughters.root', (branches_of_interest))