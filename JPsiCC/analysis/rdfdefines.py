'''
rdfdefines
+++++++++++++++++++++++++++++++++++++++++++++

Library of definitions that contain analysis category-specific RDF defines and filters.
'''
from kytools import jsonreader
import ROOT

def load_functions():
    ROOT.gSystem.AddDynamicPath("./.")
    ROOT.gROOT.ProcessLine(".include ./.")
    ROOT.gInterpreter.AddIncludePath("./.")
    ROOT.gInterpreter.Declare('#include "./config/functions.h"')

def get_rdf_branches(key):
    '''Get RDF branches to snapshot for a specific key in rdf_hists.json.

    Args:
        key (str): Key to the histogram definitions in the json file.

    Returns:
        branchlist (str): List of RDF branches in the histogram definitions.
    '''
    hist_defs = jsonreader.get_object_from_json(anpath='JPsiCC',
                                                    jsonname='rdf_hists.json',
                                                    keys=['NANOAOD_to_RDF'])
    branchlist = []
    for _, hdef in hist_defs[key].items():
        branchlist.append(hdef['name'])
    return branchlist

def compute_sum_weights(rdf_runs):
    '''Compute the weights.

    Only works if the RDF has columns named 'genEventSumw' and 'genEventCount'.

    Args:
        rdf_runs (ROOT.RDataFrame): ROOT.RDataFrame object created from the 'Runs' tree.

    Returns:
        gen_event_sumw (float): Sum of 'genEventSumw' column values.
        gen_event_count (float): Sum of 'genEventCount' column values.
    '''
    gen_event_sumw = rdf_runs.Sum('genEventSumw').GetValue()
    gen_event_count = rdf_runs.Sum('genEventCount').GetValue()
    if not gen_event_count == gen_event_sumw:
        print('WARNING in compute_sum_weights:\ngen_event_sumw not equal to gen_events_count')
    return gen_event_sumw, gen_event_count

def rdf_def_sample_meta(rdf_events):
    '''RDF define columns pertinent to the meta information.

    Args:
        rdf_events (ROOT.RDataFrame): The RDataFrame object created from the 'Events' tree
            from a JSON specification file.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): List of TBranches pertinent to the definitions.
    '''
    new_rdf = (rdf_events.DefinePerSample('sample', 'rdfsampleinfo_.GetSampleName()')
                         .DefinePerSample('sample_category', 'rdfsampleinfo_.GetS("sample_category")')
                         .DefinePerSample('xsec', 'rdfsampleinfo_.GetD("xsec")')
                         .DefinePerSample('lumi', 'rdfsampleinfo_.GetD("lumi")')
    )
    branches = ['sample', 'sample_category', 'xsec', 'lumi']
    return new_rdf, branches

def rdf_def_weights(rdf_events, sum_weights, data=False):
    '''RDF define a column with weights.

    Args:
        rdf_events (ROOT.RDataFrame): The ROOT.RDataFrame object created from the 'Events' tree.
        dict_sum_weights (float): Sum of weights (computed from compute_sum_weights).
        data (bool, optional): Whether the provided RDF is data.
            Defaults to False.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): List of TBranches pertinent to the definitions.
    '''
    # def_sum_weights = '''double def_sum_weights(unsigned int slot, const ROOT::RDF::RSampleInfo &id) {\n'''
    # counter = 0
    # for key, value in dict_sum_weights.items():
    #     if counter == 0:
    #         def_sum_weights += f'    if (id.GetS("sample_category") == "{key}") {{return {value};}}\n'
    #     else:
    #         def_sum_weights += f'    else if (id.GetS("sample_category") == "{key}") {{return {value};}}\n'
    #     counter += 1
    # def_sum_weights += '    else {return 0.0;}\n}'
    # ROOT.gInterpreter.Declare(def_sum_weights)
    if not data: weight_exp = 'xsec*genWeight*lumi/sum_weights'
    else: weight_exp = '1.0' # 'xsec*lumi/sum_weights'
    new_rdf = (rdf_events.Define('sum_weights', f'{sum_weights}')
                         .Define('w', weight_exp)
    )
    branches = ['sum_weights', 'w']
    return new_rdf, branches

def rdf_filter_triggers(rdf, CAT, YEAR):
    '''RDF define and filter with triggers.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        CAT (str): Analysis category.
        YEAR (int): Year of data-taking.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): List of TBranches pertinent to the definitions.

    Raises:
        ValueError: If the combination of 'CAT' and 'YEAR' does not correspond
            to an existing trigger.
    '''
    match f'{CAT}_{YEAR}':
        case 'GF_2018':
            trigger = '(HLT_Dimuon20_Jpsi_Barrel_Seagulls || HLT_Dimuon25_Jpsi)'
        case _:
            raise ValueError(f'Trigger does not exist for the combination of {CAT} and {YEAR}.')
    new_rdf = (rdf.Define('trigger', trigger)
                  .Filter('(trigger==1)'))
    branches = ['trigger']
    return new_rdf, branches

def rdf_def_goodjets(rdf, CAT, YEAR):
    '''RDF define with good jets definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        CAT (str): Analysis category.
        YEAR (int): Year of data-taking.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): List of TBranches pertinent to the definitions.

    Raises:
        ValueError: If the combination of 'CAT' and 'YEAR' does not correspond
            to an existing trigger.
    '''
    match f'{CAT}_{YEAR}':
        case 'GF_2018':
            goodjets = '(Jet_pt > 20 && abs(Jet_eta) < 2.5 && Jet_btagDeepCvL > -1)'
        case _:
            raise ValueError(f'Good-jets definition does not exist for the combination of {CAT} and {YEAR}.')
    new_rdf = (rdf.Define('goodjets', goodjets))
                #   .Filter('goodjets==true'))
    branches = ['goodjets']
    return new_rdf, branches

def rdf_def_jpsi(rdf):
    '''RDF definition for J/Psi.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
        branches (list(str)): List of TBranches pertinent to the definitions.
    '''
    new_rdf = (rdf.Filter('nJpsi>0', 'Event must contain at least one J/psi with the given purity criteria.')
                  )
    branches = get_rdf_branches('Jpsi')
    return new_rdf, branches

def rdf_def_jets(rdf):
    '''Jets RDF definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
    '''
    load_functions()
    TRIGGER = 'HLT_Dimuon20_Jpsi_Barrel_Seagulls or HLT_Dimuon25_Jpsi' # TODO: both data and MC (only 2018)
    GOODJETS = '(Jet_pt>20 and abs(Jet_eta)<2.5 and Jet_btagDeepCvL>-1)'
    new_rdf= (rdf.Filter('nJpsi>0', 'at least one Jpsi')
                # .Define('triggerAna', f'{TRIGGER}') # TODO: bool -> filter
                .Define('goodJets', f'{GOODJETS}')
                .Define('nGoodJets', 'Sum(goodJets)*1.0f').Filter('Sum(goodJets)>1', 'two jets for cc')
                .Define('mJJ', 'Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets])')
                # .Define('massHiggs', 'Minv3massive(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)')
                # .Define('ptHiggs', 'Ptinv3massive(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)')
                .Define('jet1Pt', 'Jet_pt[goodJets][0]')
                .Define('jet2Pt', 'Jet_pt[goodJets][1]')
                .Define('jet1Eta', 'Jet_eta[goodJets][0]')
                .Define('jet2Eta', 'Jet_eta[goodJets][1]')
                .Define('jet1CvL', 'Jet_btagDeepCvL[goodJets][0]')
                .Define('jet2CvL', 'Jet_btagDeepCvL[goodJets][1]')
                .Define('jet1partonFlavour', 'abs(Jet_partonFlavour[goodJets][0])')
                .Define('jet2partonFlavour', 'abs(Jet_partonFlavour[goodJets][1])')
                .Define('index_CloseFar', 'jetCloseFar(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets],Jpsi_kin_pt,Jpsi_kin_eta, Jpsi_kin_phi,Jpsi_kin_mass)')
                .Filter('index_CloseFar[0]!= -1', 'at least one close jet')
                .Define('jetCloseCvL', 'Jet_btagDeepCvL[goodJets][index_CloseFar[0]]')
                .Define('jetFarCvL', 'Jet_btagDeepCvL[goodJets][index_CloseFar[1]]')
                .Define('jetClosePt', 'Jet_pt[goodJets][index_CloseFar[0]]')
                .Define('jetFarPt', 'Jet_pt[goodJets][index_CloseFar[1]]')
                .Define('jetFarPtCharm', '(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==4) ? Jet_cRegCorr[goodJets][index_CloseFar[1]] * Jet_pt[goodJets][index_CloseFar[1]]:-1.')
                .Define('jetFarPtGluon', '(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==21) ? Jet_cRegCorr[goodJets][index_CloseFar[1]] * Jet_pt[goodJets][index_CloseFar[1]]:-1.')
                .Define('jetFarCvLCharm', '(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==4) ? Jet_btagDeepCvL[goodJets][index_CloseFar[1]]:-1.')
                .Define('jetFarCvLGluon', '(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==21) ? Jet_btagDeepCvL[goodJets][index_CloseFar[1]]:-1.')
                .Define('jetCloseScale', 'Jet_pt[goodJets][index_CloseFar[0]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]')
                .Define('jetFarScale', 'Jet_pt[goodJets][index_CloseFar[1]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]')
                .Define('jetClosecRegCorrScale', 'Jet_pt[goodJets][index_CloseFar[0]]*Jet_cRegCorr[goodJets][index_CloseFar[0]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]')
                .Define('jetFarcRegCorrScale', 'Jet_pt[goodJets][index_CloseFar[1]]*Jet_cRegCorr[goodJets][index_CloseFar[1]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]')     
                .Define('jetClosecRegCorr', 'Jet_cRegCorr[goodJets][index_CloseFar[0]]')
                .Define('jetFarcRegCorr', 'Jet_cRegCorr[goodJets][index_CloseFar[1]]')
                .Define('jetClosepartonFlavour', 'abs(Jet_partonFlavour[goodJets][index_CloseFar[0]])')
                .Define('jetFarpartonFlavour', 'abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])')
                .Define('jetClosenConst', 'Jet_nConstituents[goodJets][index_CloseFar[0]]')
                .Define('jetFarnConst', 'Jet_nConstituents[goodJets][index_CloseFar[1]]')
                .Define('jetCloseJPsiRatio', 'Jpsi_kin_pt[0]/Jet_pt[goodJets][index_CloseFar[0]]')
                .Define('jetClosenMuons', 'Jet_nMuons[goodJets][index_CloseFar[0]]')
                .Define('jetFarnMuons', 'Jet_nMuons[goodJets][index_CloseFar[1]]')
                .Define('jetClosenElectrons', 'Jet_nElectrons[goodJets][index_CloseFar[0]]')
                .Define('jetFarnElectrons', 'Jet_nElectrons[goodJets][index_CloseFar[1]]')
                .Define('massHiggsCorr', 'Minv3massiveCorr(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jet_cRegCorr[goodJets], Jet_muonIdx1[goodJets][index_CloseFar[0]], Jet_muonIdx2[goodJets][index_CloseFar[0]], Muon_pt, Muon_eta, Muon_phi, Jpsi_muon1_pt, Jpsi_muon1_eta, Jpsi_muon1_phi, Jpsi_muon2_pt, Jpsi_muon2_eta, Jpsi_muon2_phi, Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0])')
                .Define('minDRjpsi', 'std::min(deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]),deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]))')
                .Filter('minDRjpsi<0.3', 'jPsi very close to 1 charm-jet')
            )
    return new_rdf