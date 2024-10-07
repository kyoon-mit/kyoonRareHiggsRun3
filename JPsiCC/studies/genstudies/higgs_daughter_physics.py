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
rdf = (rdf.Define('muminus_index', 'IndexFindPDG(GenPart_HiggsDaughters_pdgId, 13)') # muon
          .Define('muplus_index', 'IndexFindPDG(GenPart_HiggsDaughters_pdgId, -13)') # anti-muon
          .Define('muminus_mass', 'HiggsDaughters_mass[muminus_index]')
          .Define('muminus_energy', 'HiggsDaughters_energy[muminus_index]')
          .Define('muminus_phi', 'HiggsDaughters_phi[muminus_index]')
          .Define('muminus_eta', 'HiggsDaughters_eta[muminus_index]')
          .Define('muminus_pt', 'HiggsDaughters_pt[muminus_index]')
          .Define('muminus_px', 'HiggsDaughters_px[muminus_index]')
          .Define('muminus_py', 'HiggsDaughters_py[muminus_index]')
          .Define('muminus_pz', 'HiggsDaughters_pz[muminus_index]')
          .Define('muplus_mass', 'HiggsDaughters_mass[muplus_index]')
          .Define('muplus_energy', 'HiggsDaughters_energy[muplus_index]')
          .Define('muplus_phi', 'HiggsDaughters_phi[muplus_index]')
          .Define('muplus_eta', 'HiggsDaughters_eta[muplus_index]')
          .Define('muplus_pt', 'HiggsDaughters_pt[muplus_index]')
          .Define('muplus_px', 'HiggsDaughters_px[muplus_index]')
          .Define('muplus_py', 'HiggsDaughters_py[muplus_index]')
          .Define('muplus_pz', 'HiggsDaughters_pz[muplus_index]')
      )

# Get muon-muon separation
rdf = (rdf.Define('dR_muminus_muplus', 'DeltaR(muminus_eta, muminus_phi, muplus_eta, muplus_phi)')
          .Define('deta_muminus_muplus', 'muminus_eta - muplus_eta')
          .Define('dphi_muminus_muplus', 'muminus_phi - muplus_phi')
          .Define('dpt_muminus_muplus', 'muminus_pt - muplus_pt')
          .Define('dE_muminus_muplus', 'muminus_energy - muplus_energy')
      )

# Get J/Psi
rdf = (rdf.Define('JPsi_cand',
                  'SumPxPyPzE(muminus_px, muminus_py, muminus_pz, muminus_energy,\
                              muplus_px, muplus_py, muplus_pz, muplus_energy)')
          .Define('JPsi_cand_mass', 'JPsi_cand.M()')
          .Define('JPsi_cand_eta', 'JPsi_cand.Eta()')
          .Define('JPsi_cand_phi', 'JPsi_cand.Phi()')
          .Define('JPsi_cand_p', 'JPsi_cand.P()')
          .Define('JPsi_cand_pt', 'JPsi_cand.Pt()')
          .Define('JPsi_cand_mt', 'JPsi_cand.Mt()')
          .Define('JPsi_cand_energy', 'JPsi_cand.E()')
      )

# Get charms
rdf = (rdf.Define('charm_index', 'IndexFindPDG(GenPart_HiggsDaughters_pdgId, 4)')  # charm
          .Define('anticharm_index', 'IndexFindPDG(GenPart_HiggsDaughters_pdgId, -4)') # anti-charm
          .Define('charm_mass', 'HiggsDaughters_mass[charm_index]')
          .Define('charm_energy', 'HiggsDaughters_energy[charm_index]')
          .Define('charm_phi', 'HiggsDaughters_phi[charm_index]')
          .Define('charm_eta', 'HiggsDaughters_eta[charm_index]')
          .Define('charm_pt', 'HiggsDaughters_pt[charm_index]')
          .Define('charm_px', 'HiggsDaughters_px[charm_index]')
          .Define('charm_py', 'HiggsDaughters_py[charm_index]')
          .Define('charm_pz', 'HiggsDaughters_pz[charm_index]')
          .Define('anticharm_mass', 'HiggsDaughters_mass[anticharm_index]')
          .Define('anticharm_energy', 'HiggsDaughters_energy[anticharm_index]')
          .Define('anticharm_phi', 'HiggsDaughters_phi[anticharm_index]')
          .Define('anticharm_eta', 'HiggsDaughters_eta[anticharm_index]')
          .Define('anticharm_pt', 'HiggsDaughters_pt[anticharm_index]')
          .Define('anticharm_px', 'HiggsDaughters_px[anticharm_index]')
          .Define('anticharm_py', 'HiggsDaughters_py[anticharm_index]')
          .Define('anticharm_pz', 'HiggsDaughters_pz[anticharm_index]')
          .Define('sorted_charm_indices', 'SortByPt(charm_index, anticharm_index, HiggsDaughters_pt)')
          .Define('leadcharm_index', 'sorted_charm_indices[0]') # leading charm
          .Define('subcharm_index', 'sorted_charm_indices[1]')  # sub-leading charm
          .Define('leadcharm_pdgId', 'GenPart_HiggsDaughters_pdgId[leadcharm_index]')
          .Define('leadcharm_mass', 'HiggsDaughters_mass[leadcharm_index]')
          .Define('leadcharm_energy', 'HiggsDaughters_energy[leadcharm_index]')
          .Define('leadcharm_phi', 'HiggsDaughters_phi[leadcharm_index]')
          .Define('leadcharm_eta', 'HiggsDaughters_eta[leadcharm_index]')
          .Define('leadcharm_pt', 'HiggsDaughters_pt[leadcharm_index]')
          .Define('leadcharm_px', 'HiggsDaughters_px[leadcharm_index]')
          .Define('leadcharm_py', 'HiggsDaughters_py[leadcharm_index]')
          .Define('leadcharm_pz', 'HiggsDaughters_pz[leadcharm_index]')
          .Define('subcharm_pdgId', 'GenPart_HiggsDaughters_pdgId[subcharm_index]')
          .Define('subcharm_mass', 'HiggsDaughters_mass[subcharm_index]')
          .Define('subcharm_energy', 'HiggsDaughters_energy[subcharm_index]')
          .Define('subcharm_phi', 'HiggsDaughters_phi[subcharm_index]')
          .Define('subcharm_eta', 'HiggsDaughters_eta[subcharm_index]')
          .Define('subcharm_pt', 'HiggsDaughters_pt[subcharm_index]')
          .Define('subcharm_px', 'HiggsDaughters_px[subcharm_index]')
          .Define('subcharm_py', 'HiggsDaughters_py[subcharm_index]')
          .Define('subcharm_pz', 'HiggsDaughters_pz[subcharm_index]')
      )

# Get charm-charm separation
rdf = (rdf.Define('dR_charm_anticharm', 'DeltaR(charm_eta, charm_phi, anticharm_eta, anticharm_phi)')
          .Define('deta_charm_anticharm', 'charm_eta - anticharm_eta')
          .Define('dphi_charm_anticharm', 'charm_phi - anticharm_phi')
          .Define('dpt_charm_anticharm', 'charm_pt - anticharm_pt')
          .Define('dpt_leadcharm_subcharm', 'leadcharm_pt - subcharm_pt')
      )

# Get di-charm
rdf = (rdf.Define('dicharm_cand',
                  'SumPxPyPzE(charm_px, charm_py, charm_pz, charm_energy,\
                              anticharm_px, anticharm_py, anticharm_pz, anticharm_energy)')
          .Define('dicharm_cand_mass', 'dicharm_cand.M()')
          .Define('dicharm_cand_eta', 'dicharm_cand.Eta()')
          .Define('dicharm_cand_phi', 'dicharm_cand.Phi()')
          .Define('dicharm_cand_p', 'dicharm_cand.P()')
          .Define('dicharm_cand_pt', 'dicharm_cand.Pt()')
          .Define('dicharm_cand_mt', 'dicharm_cand.Mt()')
          .Define('dicharm_cand_energy', 'dicharm_cand.E()')
      )

# Get charm-JPsi separation
rdf = (rdf.Define('dR_charm_JPsi', 'DeltaR(charm_eta, charm_phi, JPsi_cand_eta, JPsi_cand_phi)')
          .Define('deta_charm_JPsi', 'charm_eta - JPsi_cand_eta')
          .Define('dphi_charm_JPsi', 'charm_phi - JPsi_cand_phi')
          .Define('dpt_charm_JPsi', 'charm_pt - JPsi_cand_pt')
          .Define('dR_anticharm_JPsi', 'DeltaR(anticharm_eta, anticharm_phi, JPsi_cand_eta, JPsi_cand_phi)')
          .Define('deta_anticharm_JPsi', 'anticharm_eta - JPsi_cand_eta')
          .Define('dphi_anticharm_JPsi', 'anticharm_phi - JPsi_cand_phi')
          .Define('dpt_anticharm_JPsi', 'anticharm_pt - JPsi_cand_pt')
          .Define('dR_leadcharm_JPsi', 'DeltaR(leadcharm_eta, leadcharm_phi, JPsi_cand_eta, JPsi_cand_phi)')
          .Define('deta_leadcharm_JPsi', 'leadcharm_eta - JPsi_cand_eta')
          .Define('dphi_leadcharm_JPsi', 'leadcharm_phi - JPsi_cand_phi')
          .Define('dpt_leadcharm_JPsi', 'leadcharm_pt - JPsi_cand_pt')
          .Define('dR_subcharm_JPsi', 'DeltaR(subcharm_eta, subcharm_phi, JPsi_cand_eta, JPsi_cand_phi)')
          .Define('deta_subcharm_JPsi', 'subcharm_eta - JPsi_cand_eta')
          .Define('dphi_subcharm_JPsi', 'subcharm_phi - JPsi_cand_phi')
          .Define('dpt_subcharm_JPsi', 'subcharm_pt - JPsi_cand_pt')
      )

# Get dicharm-JPsi separation
rdf = (rdf.Define('dR_dicharm_JPsi', 'DeltaR(dicharm_cand_eta, dicharm_cand_phi, JPsi_cand_eta, JPsi_cand_phi)')
          .Define('deta_dicharm_JPsi', 'dicharm_cand_eta - JPsi_cand_eta')
          .Define('dphi_dicharm_JPsi', 'dicharm_cand_phi - JPsi_cand_phi')
          .Define('dpt_dicharm_JPsi', 'dicharm_cand_pt - JPsi_cand_pt')
      )

# Get Higgs
rdf = (rdf.Define('Higgs_cand',
                  'SumPtEtaPhiE(dicharm_cand_pt, dicharm_cand_eta, dicharm_cand_phi, dicharm_cand_energy,\
                                JPsi_cand_pt, JPsi_cand_eta, JPsi_cand_phi, JPsi_cand_energy)')
          .Define('Higgs_cand_mass', 'Higgs_cand.M()')
          .Define('Higgs_cand_eta', 'Higgs_cand.Eta()')
          .Define('Higgs_cand_phi', 'Higgs_cand.Phi()')
          .Define('Higgs_cand_p', 'Higgs_cand.P()')
          .Define('Higgs_cand_pt', 'Higgs_cand.Pt()')
          .Define('Higgs_cand_mt', 'Higgs_cand.Mt()')
          .Define('Higgs_cand_energy', 'Higgs_cand.E()')
      )

# Snapshot
branches_of_interest =\
    ['JPsi_cand_mass',
     'JPsi_cand_eta',
     'JPsi_cand_phi',
     'JPsi_cand_p',
     'JPsi_cand_pt',
     'JPsi_cand_mt',
     'JPsi_cand_energy',
     'dicharm_cand_mass',
     'dicharm_cand_eta',
     'dicharm_cand_phi',
     'dicharm_cand_p',
     'dicharm_cand_pt',
     'dicharm_cand_mt',
     'dicharm_cand_energy',
     'Higgs_cand_mass',
     'Higgs_cand_eta',
     'Higgs_cand_phi',
     'Higgs_cand_p',
     'Higgs_cand_pt',
     'Higgs_cand_mt',
     'Higgs_cand_energy',
     'genWeight'
    ]

branches_separation =\
    ['dR_muminus_muplus',
     'deta_muminus_muplus',
     'dphi_muminus_muplus',
     'dpt_muminus_muplus',
     'dE_muminus_muplus',
     'dR_charm_anticharm',
     'deta_charm_anticharm',
     'dphi_charm_anticharm',
     'dpt_charm_anticharm',
     'dpt_leadcharm_subcharm',
     'dR_charm_JPsi',
     'deta_charm_JPsi',
     'dphi_charm_JPsi',
     'dpt_charm_JPsi',
     'dR_anticharm_JPsi',
     'deta_anticharm_JPsi',
     'dphi_anticharm_JPsi',
     'dpt_anticharm_JPsi',
     'dR_leadcharm_JPsi',
     'deta_leadcharm_JPsi',
     'dphi_leadcharm_JPsi',
     'dpt_leadcharm_JPsi',
     'dR_subcharm_JPsi',
     'deta_subcharm_JPsi',
     'dphi_subcharm_JPsi',
     'dpt_subcharm_JPsi',
     'dR_dicharm_JPsi',
     'deta_dicharm_JPsi',
     'dphi_dicharm_JPsi',
     'dpt_dicharm_JPsi'
    ]

branches_individual =\
    ['muminus_mass',
     'muminus_energy',
     'muminus_phi',
     'muminus_eta',
     'muminus_px',
     'muminus_py',
     'muminus_pz',
     'muplus_mass',
     'muplus_energy',
     'muplus_phi',
     'muplus_eta',
     'muplus_px',
     'muplus_py',
     'muplus_pz',
     'charm_mass',
     'charm_energy',
     'charm_phi',
     'charm_eta',
     'charm_px',
     'charm_py',
     'charm_pz',
     'anticharm_mass',
     'anticharm_energy',
     'anticharm_phi',
     'anticharm_eta',
     'anticharm_px',
     'anticharm_py',
     'anticharm_pz',
     'leadcharm_pdgId',
     'leadcharm_mass',
     'leadcharm_energy',
     'leadcharm_phi',
     'leadcharm_eta',
     'leadcharm_px',
     'leadcharm_py',
     'leadcharm_pz',
     'subcharm_pdgId',
     'subcharm_mass',
     'subcharm_energy',
     'subcharm_phi',
     'subcharm_eta',
     'subcharm_px',
     'subcharm_py',
     'subcharm_pz',
    ]

branches_snapshot = branches_of_interest + branches_separation + branches_individual
branches_plot = branches_snapshot

rdf.Snapshot('GenEvents', 'daughters.root', (branches_snapshot))

# Get histograms
my_hists = {}
for var in branches_plot:
    print(var)
    hist = None
    if var.endswith('pt') or var.endswith('mt'):
        hist = rdf.Histo1D(ROOT.RDF.TH1DModel(var, var, 60, 0, 300), var, 'genWeight')
    elif var.endswith('p') or var.endswith('energy'):
        hist = rdf.Histo1D(ROOT.RDF.TH1DModel(var, var, 400, 0, 2000), var, 'genWeight')
    elif var.endswith('eta'):
        hist = rdf.Histo1D(ROOT.RDF.TH1DModel(var, var, 100, -10., 10.), var, 'genWeight')
    elif var.endswith('phi'):
        hist = rdf.Histo1D(ROOT.RDF.TH1DModel(var, var, 40, -4., 4.), var, 'genWeight')
    elif var.startswith('dR'):
        hist = rdf.Histo1D(ROOT.RDF.TH1DModel(var, var, 50, 0., 10.), var, 'genWeight')
    elif var.startswith('deta'):
        hist = rdf.Histo1D(ROOT.RDF.TH1DModel(var, var, 100, -10., 10.), var, 'genWeight')
    elif var.startswith('dphi'):
        hist = rdf.Histo1D(ROOT.RDF.TH1DModel(var, var, 100, -10., 10.), var, 'genWeight')
    else:
        hist = rdf.Histo1D(var, 'genWeight')
    my_hists[var] = hist

# Draw
save_dir = os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'plots', 'genstudies', 'genvar')
for key, hist in my_hists.items():
    c = ROOT.TCanvas()
    hist.Draw()
    c.SaveAs(os.path.join(save_dir, f'{key}.png'))
    c.Close()