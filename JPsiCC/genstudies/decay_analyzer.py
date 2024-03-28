'''
This module examines the decay of a particle using the gen-level mother and daughter information.
'''

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

# Find Higgs Daughter
ROOT.gInterpreter.Declare(
"""
using namespace ROOT;
RVec<int> HiggsDaughtersPDG(
    const RVecI& GenPart_pdgId,
    const RVecI& GenPart_genPartIdxMother
) {
    RVec<int> daughters_pdg;
    int idx_low = 0;
    int idx_high = GenPart_pdgId.size();
    for (int i=idx_high; i>=idx_low; i--) {
        if (GenPart_genPartIdxMother[i] < 0 || GenPart_genPartIdxMother[i] > idx_high) continue; // TODO: why these nonsensical indices?
        if (GenPart_pdgId[GenPart_genPartIdxMother[i]]==25 && // if mother is Higgs
            GenPart_pdgId[i]!=25 // but itself is not the Higgs
            ) {
            daughters_pdg.push_back(GenPart_pdgId[i]);
        }
    }
    return daughters_pdg;
}
""")

rdf = rdf.Define('GenPart_HiggsDaughters_pdgId',
                 'HiggsDaughtersPDG(GenPart_pdgId, GenPart_genPartIdxMother)')
rdf.Snapshot('GenEvents', 'genevents.root', ['GenPart_HiggsDaughters_pdgId'])

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