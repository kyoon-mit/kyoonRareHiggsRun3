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
        print(f'gen_event_sumw: {gen_event_sumw}, gen_event_count: {gen_event_count}')
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
    if not data: weight_exp = 'xsec*genWeight*lumi/sum_weights'
    else: weight_exp = '1.0'
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
    try: new_rdf = rdf.Define('trigger', trigger)
    except: new_rdf = rdf
    new_rdf = new_rdf.Filter('(trigger==1)')
    branches = ['trigger']
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

def rdf_def_muons(rdf):
    '''Muons RDF definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
    '''
    new_rdf = rdf
    branches = get_rdf_branches('muon')
    return new_rdf, branches

def rdf_def_jets(rdf, CAT, YEAR):
    '''Jets RDF definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        CAT (str): Analysis category.
        YEAR (int): Year of data-taking.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
        branches (list(str)): List of TBranches pertinent to the definitions.

    Raises:
        ValueError: If the combination of 'CAT' and 'YEAR' does not correspond
            to an existing trigger.
    '''
    load_functions()
    match f'{CAT}_{YEAR}':
        case 'GF_2018':
            goodjets = '(Jet_pt > 20 && abs(Jet_eta) < 2.5 && Jet_btagDeepCvL > -1)'
        case _:
            raise ValueError(f'Good-jets definition does not exist for the combination of {CAT} and {YEAR}.')
    new_rdf = (rdf.Define('goodjets', goodjets))
    branches = ['goodjets']
    branches += get_rdf_branches('jet')
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