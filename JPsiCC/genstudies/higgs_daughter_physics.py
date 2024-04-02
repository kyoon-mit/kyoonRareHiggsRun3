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
          .Define('muon1_phi', 'HiggsDaughters_phi[muon1_index][0]')
          .Define('muon1_eta', 'HiggsDaughters_eta[muon1_index][0]')
          .Define('muon1_pt', 'HiggsDaughters_pt[muon1_index][0]')
          .Define('muon1_px', 'HiggsDaughters_px[muon1_index][0]')
          .Define('muon1_py', 'HiggsDaughters_py[muon1_index][0]')
          .Define('muon1_pz', 'HiggsDaughters_pz[muon1_index][0]')
          .Define('muon2_mass', 'HiggsDaughters_mass[muon2_index][0]')
          .Define('muon2_energy', 'HiggsDaughters_energy[muon2_index][0]')
          .Define('muon2_phi', 'HiggsDaughters_phi[muon2_index][0]')
          .Define('muon2_eta', 'HiggsDaughters_eta[muon2_index][0]')
          .Define('muon2_pt', 'HiggsDaughters_pt[muon2_index][0]')
          .Define('muon2_px', 'HiggsDaughters_px[muon2_index][0]')
          .Define('muon2_py', 'HiggsDaughters_py[muon2_index][0]')
          .Define('muon2_pz', 'HiggsDaughters_pz[muon2_index][0]')
      )

# Get muon-muon separation
rdf = (rdf.Define('dR_muplus_muminus', 'DeltaR(muon1_eta, muon1_phi, muon2_eta, muon2_phi)'))

# Get J/Psi
rdf = (rdf.Define('JPsi_cand1',
                  'SumPxPyPzE(muon1_px, muon1_py, muon1_pz, muon1_energy,\
                              muon2_px, muon2_py, muon2_pz, muon2_energy)')
          .Define('JPsi_cand2',
                  'SumPxPyPzM(muon1_px, muon1_py, muon1_pz, muon1_mass,\
                              muon2_px, muon2_py, muon2_pz, muon2_mass)')
          .Define('JPsi_cand3',
                  'SumPtEtaPhiE(muon1_pt, muon1_eta, muon1_phi, muon1_energy,\
                                muon2_pt, muon2_eta, muon2_phi, muon2_energy)')
          .Define('JPsi_cand1_mass', 'JPsi_cand1.M()')
          .Define('JPsi_cand1_eta', 'JPsi_cand1.Eta()')
          .Define('JPsi_cand1_phi', 'JPsi_cand1.Phi()')
          .Define('JPsi_cand1_p', 'JPsi_cand1.P()')
          .Define('JPsi_cand1_pt', 'JPsi_cand1.Pt()')
          .Define('JPsi_cand1_mt', 'JPsi_cand1.Mt()')
          .Define('JPsi_cand1_energy', 'JPsi_cand1.E()')
          .Define('JPsi_cand2_mass', 'JPsi_cand2.M()')
          .Define('JPsi_cand2_eta', 'JPsi_cand2.Eta()')
          .Define('JPsi_cand2_phi', 'JPsi_cand2.Phi()')
          .Define('JPsi_cand2_p', 'JPsi_cand2.P()')
          .Define('JPsi_cand2_pt', 'JPsi_cand2.Pt()')
          .Define('JPsi_cand2_mt', 'JPsi_cand2.Mt()')
          .Define('JPsi_cand2_energy', 'JPsi_cand2.E()')
          .Define('JPsi_cand3_mass', 'JPsi_cand3.M()')
          .Define('JPsi_cand3_eta', 'JPsi_cand3.Eta()')
          .Define('JPsi_cand3_phi', 'JPsi_cand3.Phi()')
          .Define('JPsi_cand3_p', 'JPsi_cand3.P()')
          .Define('JPsi_cand3_pt', 'JPsi_cand3.Pt()')
          .Define('JPsi_cand3_mt', 'JPsi_cand3.Mt()')
          .Define('JPsi_cand3_energy', 'JPsi_cand3.E()')
      )

# Get charms
rdf = (rdf.Define('charm1_index', 'GenPart_HiggsDaughters_pdgId==4')
          .Define('charm2_index', 'GenPart_HiggsDaughters_pdgId==-4')
          .Define('charm1_mass', 'HiggsDaughters_mass[charm1_index][0]')
          .Define('charm1_energy', 'HiggsDaughters_energy[charm1_index][0]')
          .Define('charm1_phi', 'HiggsDaughters_phi[charm1_index][0]')
          .Define('charm1_eta', 'HiggsDaughters_eta[charm1_index][0]')
          .Define('charm1_pt', 'HiggsDaughters_pt[charm1_index][0]')
          .Define('charm1_px', 'HiggsDaughters_px[charm1_index][0]')
          .Define('charm1_py', 'HiggsDaughters_py[charm1_index][0]')
          .Define('charm1_pz', 'HiggsDaughters_pz[charm1_index][0]')
          .Define('charm2_mass', 'HiggsDaughters_mass[charm2_index][0]')
          .Define('charm2_energy', 'HiggsDaughters_energy[charm2_index][0]')
          .Define('charm2_phi', 'HiggsDaughters_phi[charm2_index][0]')
          .Define('charm2_eta', 'HiggsDaughters_eta[charm2_index][0]')
          .Define('charm2_pt', 'HiggsDaughters_pt[charm2_index][0]')
          .Define('charm2_px', 'HiggsDaughters_px[charm2_index][0]')
          .Define('charm2_py', 'HiggsDaughters_py[charm2_index][0]')
          .Define('charm2_pz', 'HiggsDaughters_pz[charm2_index][0]')
      )

# Get charm-charm separation
rdf = (rdf.Define('dR_c_cbar', 'DeltaR(charm1_eta, charm1_phi, charm2_eta, charm2_phi)'))

# Get di-charm
rdf = (rdf.Define('dicharm_cand1',
                  'SumPxPyPzE(charm1_px, charm1_py, charm1_pz, charm1_energy,\
                              charm2_px, charm2_py, charm2_pz, charm2_energy)')
          .Define('dicharm_cand2',
                  'SumPxPyPzM(charm1_px, charm1_py, charm1_pz, charm1_mass,\
                              charm2_px, charm2_py, charm2_pz, charm2_mass)')
          .Define('dicharm_cand3',
                  'SumPtEtaPhiE(charm1_pt, charm1_eta, charm1_phi, charm1_energy,\
                                charm2_pt, charm2_eta, charm2_phi, charm2_energy)')
          .Define('dicharm_cand1_mass', 'dicharm_cand1.M()')
          .Define('dicharm_cand1_eta', 'dicharm_cand1.Eta()')
          .Define('dicharm_cand1_phi', 'dicharm_cand1.Phi()')
          .Define('dicharm_cand1_p', 'dicharm_cand1.P()')
          .Define('dicharm_cand1_pt', 'dicharm_cand1.Pt()')
          .Define('dicharm_cand1_mt', 'dicharm_cand1.Mt()')
          .Define('dicharm_cand1_energy', 'dicharm_cand1.E()')
          .Define('dicharm_cand2_mass', 'dicharm_cand2.M()')
          .Define('dicharm_cand2_eta', 'dicharm_cand2.Eta()')
          .Define('dicharm_cand2_phi', 'dicharm_cand2.Phi()')
          .Define('dicharm_cand2_p', 'dicharm_cand2.P()')
          .Define('dicharm_cand2_pt', 'dicharm_cand2.Pt()')
          .Define('dicharm_cand2_mt', 'dicharm_cand2.Mt()')
          .Define('dicharm_cand2_energy', 'dicharm_cand2.E()')
          .Define('dicharm_cand3_mass', 'dicharm_cand3.M()')
          .Define('dicharm_cand3_eta', 'dicharm_cand3.Eta()')
          .Define('dicharm_cand3_phi', 'dicharm_cand3.Phi()')
          .Define('dicharm_cand3_p', 'dicharm_cand3.P()')
          .Define('dicharm_cand3_pt', 'dicharm_cand3.Pt()')
          .Define('dicharm_cand3_mt', 'dicharm_cand3.Mt()')
          .Define('dicharm_cand3_energy', 'dicharm_cand3.E()')
      )

# Get charm-JPsi separation
rdf = (rdf.Define('dR_c_JPsi', 'DeltaR(charm1_eta, charm1_phi, JPsi_cand1_eta, JPsi_cand1_phi)'))
rdf = (rdf.Define('dR_cbar_JPsi', 'DeltaR(charm2_eta, charm2_phi, JPsi_cand1_eta, JPsi_cand1_phi)'))

# Get dicharm-JPsi separation
rdf = (rdf.Define('dR_dicharm_JPsi', 'DeltaR(dicharm_cand1_eta, dicharm_cand1_phi, JPsi_cand1_eta, JPsi_cand1_phi)'))

# Get Higgs
rdf = (rdf.Define('Higgs_cand1',
                  'SumPtEtaPhiE(dicharm_cand1_pt, dicharm_cand1_eta, dicharm_cand1_phi, dicharm_cand1_energy,\
                                JPsi_cand1_pt, JPsi_cand1_eta, JPsi_cand1_phi, JPsi_cand1_energy)')
          .Define('Higgs_cand2',
                  'SumPtEtaPhiE(dicharm_cand2_pt, dicharm_cand2_eta, dicharm_cand2_phi, dicharm_cand2_energy,\
                                JPsi_cand2_pt, JPsi_cand2_eta, JPsi_cand2_phi, JPsi_cand2_energy)')
          .Define('Higgs_cand3',
                  'SumPtEtaPhiE(dicharm_cand3_pt, dicharm_cand3_eta, dicharm_cand3_phi, dicharm_cand3_energy,\
                                JPsi_cand3_pt, JPsi_cand3_eta, JPsi_cand3_phi, JPsi_cand3_energy)')
          .Define('Higgs_cand1_mass', 'Higgs_cand1.M()')
          .Define('Higgs_cand1_eta', 'Higgs_cand1.Eta()')
          .Define('Higgs_cand1_phi', 'Higgs_cand1.Phi()')
          .Define('Higgs_cand1_p', 'Higgs_cand1.P()')
          .Define('Higgs_cand1_pt', 'Higgs_cand1.Pt()')
          .Define('Higgs_cand1_mt', 'Higgs_cand1.Mt()')
          .Define('Higgs_cand1_energy', 'Higgs_cand1.E()')
          .Define('Higgs_cand2_mass', 'Higgs_cand2.M()')
          .Define('Higgs_cand2_eta', 'Higgs_cand2.Eta()')
          .Define('Higgs_cand2_phi', 'Higgs_cand2.Phi()')
          .Define('Higgs_cand2_p', 'Higgs_cand2.P()')
          .Define('Higgs_cand2_pt', 'Higgs_cand2.Pt()')
          .Define('Higgs_cand2_mt', 'Higgs_cand2.Mt()')
          .Define('Higgs_cand2_energy', 'Higgs_cand2.E()')
          .Define('Higgs_cand3_mass', 'Higgs_cand3.M()')
          .Define('Higgs_cand3_eta', 'Higgs_cand3.Eta()')
          .Define('Higgs_cand3_phi', 'Higgs_cand3.Phi()')
          .Define('Higgs_cand3_p', 'Higgs_cand3.P()')
          .Define('Higgs_cand3_pt', 'Higgs_cand3.Pt()')
          .Define('Higgs_cand3_mt', 'Higgs_cand3.Mt()')
          .Define('Higgs_cand3_energy', 'Higgs_cand3.E()')
      )

# Snapshot
branches_of_interest =\
    ['JPsi_cand1_mass',
     'JPsi_cand1_eta',
     'JPsi_cand1_phi',
     'JPsi_cand1_p',
     'JPsi_cand1_pt',
     'JPsi_cand1_mt',
     'JPsi_cand1_energy',
     'JPsi_cand2_mass',
     'JPsi_cand2_eta',
     'JPsi_cand2_phi',
     'JPsi_cand2_p',
     'JPsi_cand2_pt',
     'JPsi_cand2_mt',
     'JPsi_cand2_energy',
     'JPsi_cand3_mass',
     'JPsi_cand3_eta',
     'JPsi_cand3_phi',
     'JPsi_cand3_p',
     'JPsi_cand3_pt',
     'JPsi_cand3_mt',
     'JPsi_cand3_energy',
     'dicharm_cand1_mass',
     'dicharm_cand1_eta',
     'dicharm_cand1_phi',
     'dicharm_cand1_p',
     'dicharm_cand1_pt',
     'dicharm_cand1_mt',
     'dicharm_cand1_energy',
     'dicharm_cand2_mass',
     'dicharm_cand2_eta',
     'dicharm_cand2_phi',
     'dicharm_cand2_p',
     'dicharm_cand2_pt',
     'dicharm_cand2_mt',
     'dicharm_cand2_energy',
     'dicharm_cand3_mass',
     'dicharm_cand3_eta',
     'dicharm_cand3_phi',
     'dicharm_cand3_p',
     'dicharm_cand3_pt',
     'dicharm_cand3_mt',
     'dicharm_cand3_energy',
     'Higgs_cand1_mass',
     'Higgs_cand1_eta',
     'Higgs_cand1_phi',
     'Higgs_cand1_p',
     'Higgs_cand1_pt',
     'Higgs_cand1_mt',
     'Higgs_cand1_energy',
     'Higgs_cand2_mass',
     'Higgs_cand2_eta',
     'Higgs_cand2_phi',
     'Higgs_cand2_p',
     'Higgs_cand2_pt',
     'Higgs_cand2_mt',
     'Higgs_cand2_energy',
     'Higgs_cand3_mass',
     'Higgs_cand3_eta',
     'Higgs_cand3_phi',
     'Higgs_cand3_p',
     'Higgs_cand3_pt',
     'Higgs_cand3_mt',
     'Higgs_cand3_energy',
     'dR_muplus_muminus',
     'dR_c_JPsi',
     'dR_cbar_JPsi',
     'dR_dicharm_JPsi']

branches_of_second_interest =\
    ['muon1_mass',
     'muon1_energy',
     'muon1_phi',
     'muon1_eta',
     'muon1_px',
     'muon1_py',
     'muon1_pz',
     'muon2_mass',
     'muon2_energy',
     'muon2_phi',
     'muon2_eta',
     'muon2_px',
     'muon2_py',
     'muon2_pz',
     'charm1_mass',
     'charm1_energy',
     'charm1_phi',
     'charm1_eta',
     'charm1_px',
     'charm1_py',
     'charm1_pz',
     'charm2_mass',
     'charm2_energy',
     'charm2_phi',
     'charm2_eta',
     'charm2_px',
     'charm2_py',
     'charm2_pz']

rdf.Snapshot('GenEvents', 'daughters.root', (branches_of_interest + branches_of_second_interest))

# Get histograms
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