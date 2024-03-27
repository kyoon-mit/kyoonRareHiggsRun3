from kytools import jsonreader
import os
import ROOT

# Disable multithreading for now
ROOT.DisableImplicitMT()

# Get TChain of generator-level events in the NANOAOD format
genevents = jsonreader.get_chain_from_json(anpath='JPsiCC',
                                           jsonname='NANOAOD.json',
                                           keys=['mariadlf', '20240321', 'GEN-signal'],
                                           treename='Events')

# Create RDF from TChain
rdf = ROOT.RDataFrame(genevents)

# Get some histograms
vars_of_interest = ['Jpsi_doca',
                    'Jpsi_gen_mass',
                    'Jpsi_gen_pt',
                    'Jpsi_gen_trk1_pt',
                    'Jpsi_gen_trk2_pt',
                    'Jpsi_iso',
                    'Jpsi_kin_cosAlphaXY',
                    'Jpsi_kin_eta',
                    'Jpsi_kin_lxy',
                    'Jpsi_kin_mass',
                    'Jpsi_kin_massErr',
                    'Jpsi_kin_phi',
                    'Jpsi_kin_pt',
                    'Jpsi_kin_sipBS',
                    'Jpsi_kin_sipPV',
                    'Jpsi_kin_slxy',
                    'Jpsi_kin_vtx_chi2dof',
                    'Jpsi_kin_vtx_prob',
                    'Jpsi_mass',
                    'Jpsi_muon1_eta',
                    'Jpsi_muon1_phi',
                    'Jpsi_muon1_pt',
                    'Jpsi_muon2_eta',
                    'Jpsi_muon2_phi',
                    'Jpsi_muon2_pt',
                    'Jpsi_gen_pdgId',
                    'Jpsi_gen_trk1_mpdgId',
                    'Jpsi_gen_trk1_pdgId',
                    'Jpsi_gen_trk2_mpdgId',
                    'Jpsi_gen_trk2_pdgId',
                    'Jpsi_kin_valid',
                    'Jpsi_muon1_isMediumMuon',
                    'Jpsi_muon1_isTighMuon',
                    'Jpsi_muon2_isMediumMuon',
                    'Jpsi_muon2_isTighMuon',
                    'nMuon']
my_hists = []
for var in vars_of_interest:
    print(var)
    my_hists.append(rdf.Histo1D(var))

# Draw
save_dir = os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'plots', 'genstudies')
for hist in my_hists:
    c = ROOT.TCanvas()
    hist.Draw()
    c.SaveAs(os.path.join(save_dir, f'{hist.GetTitle()}.png'))
    c.Close()