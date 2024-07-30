'''
rdfdefines
+++++++++++++++++++++++++++++++++++++++++++++

Library of definitions that contain analysis category-specific RDF defines and filters.
'''
from kytools import jsonreader
import ROOT, os
from datetime import datetime

ROOT.EnableImplicitMT

def load_functions(mode='reco'):
    '''Load user-defined functions.

    Args:
        mode (str, optional): Mode of loading.
            Options are 'reco' and 'gen'.

    Returns:
        (None)
    '''
    match mode:
        case 'reco':
            ROOT.gSystem.CompileMacro(os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'src', 'functions.cc'), 'k')
        case 'gen':
            ROOT.gSystem.CompileMacro(os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'src', 'GenAnalyzer.cc'), 'k')

def add_cut_flow(cut_flow_dict, cut_str, cut_nentries):
    '''Add the cut flow str and value to the dictionary containing the cut flow.

    Each cut_str is modified so that the Nth entry will have a prefix 'N.'.

    Args:
        cut_flow_dict (dict(str: float)): The dictionary containing the cut flow.
        cut_str (str): The string describing the cut.
        cut_nentries (int or float): The number of entries after the cut.

    Returns:
        cut_flow_dict (dict(str: float)): Updated dictionary.

    Raises:
        TypeError: If cut_str is not str.
        TypeError: If cut_netries is not numerical.
    '''
    if not isinstance(cut_str, str): raise TypeError('cut_str must be str.')
    if not isinstance(cut_nentries, (int, float)): raise TypeError('cut_nentries must be numerical.')
    nth_entry = len(cut_flow_dict)
    cut_flow_dict[f'#{nth_entry} {cut_str}'] = float(cut_nentries)
    return cut_flow_dict

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

def compute_sum_weights(rdf_runs, sample_name):
    '''Compute the weights.

    Only works if the RDF has columns named 'genEventSumw' and 'genEventCount'.

    Args:
        rdf_runs (ROOT.RDataFrame): ROOT.RDataFrame object created from the 'Runs' tree.
        sample_name (str): Name of the sample.

    Returns:
        gen_event_sumw (float): Sum of 'genEventSumw' column values.
        gen_event_count (float): Sum of 'genEventCount' column values.
    '''
    gen_event_sumw = rdf_runs.Filter(f'sample=="{sample_name}"').Sum('genEventSumw').GetValue()
    gen_event_count = rdf_runs.Filter(f'sample=="{sample_name}"').Sum('genEventCount').GetValue()
    if not gen_event_count == gen_event_sumw:
        print('WARNING in compute_sum_weights:\ngen_event_sumw not equal to gen_events_count')
        print(f'gen_event_sumw: {gen_event_sumw}, gen_event_count: {gen_event_count}')
    return gen_event_sumw, gen_event_count

def rdf_def_sample_meta(rdf_events, branches=[]):
    '''RDF define columns pertinent to the meta information.

    Args:
        rdf_events (ROOT.RDataFrame): The RDataFrame object created from the 'Events' tree
            from a JSON specification file.
        branches (list(str), optional): List of TBranches.
            Defaults to [].

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): Updated list of TBranches.
    '''
    new_rdf = (rdf_events.DefinePerSample('sample', 'rdfsampleinfo_.GetSampleName()')
                        .DefinePerSample('sample_category', 'rdfsampleinfo_.GetS("sample_category")')
                        .DefinePerSample('xsec', 'rdfsampleinfo_.GetD("xsec")')
                        .DefinePerSample('lumi', 'rdfsampleinfo_.GetD("lumi")')
    )
    branches += ['sample', 'sample_category', 'xsec', 'lumi']
    return new_rdf, branches

def rdf_def_weights(rdf_runs, rdf_events, branches=[], cut_flow_dict={}, data=False):
    '''RDF define a column with weights.

    Args:
        rdf_runs (ROOT.RDataFrame): ROOT.RDataFrame object created from the 'Runs' tree.
            Only works if the RDF has columns named 'genEventSumw' and 'genEventCount'.
            Can be (None) if data=True.
        rdf_events (ROOT.RDataFrame): The ROOT.RDataFrame object created from the 'Events' tree.
        branches (list(str), optional): List of TBranches.
            Defaults to [].
        cut_flow_dict (dict(str: float), optional): The dictionary containing the cut flow.
            Defaults to {}.
        data (bool, optional): Whether the provided RDF is data.
            Defaults to False.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): Updated list of TBranches.
        cut_flow_dict (dict(str: float)): Updated dictionary containing the cut flow.
    '''
    hash = datetime.now().strftime("%H%M%S%f")
    if not data:
        sample_col =rdf_events.Take['string']('sample')
        sample_names = list(set(sample_col.GetValue()))
        list_sum_weights = []
        for sname in sample_names:
            list_sum_weights.append((sname, compute_sum_weights(rdf_runs, sname)))
        func_sum_weights = f'float func_sum_weights_{hash}(unsigned int slot, const ROOT::RDF::RSampleInfo &id) {{\n'
        max_idx = len(list_sum_weights)
        for idx in range(max_idx):
            item = list_sum_weights[idx]
            if idx == 0:
                func_sum_weights += f'  if (id.Contains("{item[0]}")) {{return {item[1][0]};}}\n'
            else:
                func_sum_weights += f'  else if (id.Contains("{item[0]}")) {{return {item[1][0]};}}\n'
        func_sum_weights += f'  else {{return 1.;}}\n}}'
        weight_exp = 'xsec*genWeight*lumi/sum_weights'
    else:
        weight_exp = '1.0'
        func_sum_weights = f'float func_sum_weights_{hash}(unsigned int slot, const ROOT::RDF::RSampleInfo &id) {{return 1.;}}'
    ROOT.gInterpreter.Declare(func_sum_weights)
    new_rdf = (rdf_events.DefinePerSample('sum_weights', f'func_sum_weights_{hash}(rdfslot_, rdfsampleinfo_)')
                        .Define('w', weight_exp)
    )
    cut_flow_dict = add_cut_flow(cut_flow_dict, 'w', new_rdf.Sum('w').GetValue())
    branches += ['sum_weights', 'w']
    if not data: branches.append('genWeight')
    return new_rdf, branches, cut_flow_dict

def rdf_filter_triggers(rdf, CAT, YEAR, branches=[], cut_flow_dict={}):
    '''RDF define and filter with triggers.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        CAT (str): Analysis category.
        YEAR (int): Year of data-taking.
        branches (list(str), optional): List of TBranches.
            Defaults to [].
        cut_flow_dict (dict(str: float), optional): The dictionary containing the cut flow.
            Defaults to {}.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): Updated list of TBranches.
        cut_flow_dict (dict(str: float)): Updated dictionary containing the cut flow.

    Raises:
        ValueError: If the combination of 'CAT' and 'YEAR' does not correspond
            to an existing trigger.
    '''
    match f'{CAT}_{YEAR}':
        case 'GF_2018':
            # trigger = '(HLT_Dimuon20_Jpsi_Barrel_Seagulls || HLT_Dimuon25_Jpsi)'
            trigger = 'HLT_Dimuon25_Jpsi'
            print(f'The following triggers are applied: {trigger}')
        case _:
            raise ValueError(f'Trigger does not exist for the combination of {CAT} and {YEAR}.')
    try: new_rdf = rdf.Define('trigger', trigger)
    except: new_rdf = rdf
    new_rdf = new_rdf.Filter('trigger>0')
    cut_flow_dict = add_cut_flow(cut_flow_dict, 'trigger', new_rdf.Sum('w').GetValue())
    branches += ['trigger']
    return new_rdf, branches, cut_flow_dict

def rdf_def_jpsi(rdf, branches=[], cut_flow_dict={}):
    '''RDF definition for J/Psi.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        branches (list(str), optional): List of TBranches.
            Defaults to [].
        cut_flow_dict (dict(str: float), optional): The dictionary containing the cut flow.
            Defaults to {}.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
        branches (list(str)): Updated list of TBranches.
        cut_flow_dict (dict(str: float)): Updated dictionary containing the cut flow.
    '''
    new_rdf = (rdf.Filter('nJpsi>0', 'Event must contain at least one J/psi with the given purity criteria.')
    )
    cut_flow_dict = add_cut_flow(cut_flow_dict, 'nJpsi>0', new_rdf.Sum('w').GetValue())
    branches += get_rdf_branches('Jpsi')
    return new_rdf, branches, cut_flow_dict

def rdf_def_muons(rdf, branches=[], cut_flow_dict={}):
    '''Muons RDF definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        branches (list(str), optional): List of TBranches.
            Defaults to [].
        cut_flow_dict (dict(str: float), optional): The dictionary containing the cut flow.
            Defaults to {}.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
        branches (list(str)): Updated list of TBranches.
        cut_flow_dict (dict(str: float)): Updated dictionary containing the cut flow.
    '''
    new_rdf = rdf
    branches += get_rdf_branches('muon')
    return new_rdf, branches, cut_flow_dict

def rdf_def_vertex(rdf, branches=[], cut_flow_dict={}):
    '''Vertex RDF definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        branches (list(str), optional): List of TBranches.
            Defaults to [].
        cut_flow_dict (dict(str: float), optional): The dictionary containing the cut flow.
            Defaults to {}.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
        branches (list(str)): Updated list of TBranches.
        cut_flow_dict (dict(str: float)): Updated dictionary containing the cut flow.
    '''
    new_rdf = (rdf.Filter('PV_npvsGood>0', 'Number of PVs must be at least 1.')
    )
    cut_flow_dict = add_cut_flow(cut_flow_dict, 'PV_npvsGood>0', new_rdf.Sum('w').GetValue())
    branches += get_rdf_branches('vertex')
    return new_rdf, branches, cut_flow_dict

def rdf_def_jets(rdf, CAT, YEAR, branches=[], cut_flow_dict={}, data=False):
    '''Jets RDF definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        CAT (str): Analysis category.
        YEAR (int): Year of data-taking.
        branches (list(str), optional): List of TBranches.
            Defaults to [].
        cut_flow_dict (dict(str: float), optional): The dictionary containing the cut flow.
            Defaults to {}.
        data (bool, optional): Whether the provided RDF is data.
            Defaults to False.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
        branches (list(str)): Updated list of TBranches.
        cut_flow_dict (dict(str: float)): Updated dictionary containing the cut flow.

    Raises:
        ValueError: If the combination of 'CAT' and 'YEAR' does not correspond
            to an existing trigger.
    '''
    load_functions('reco')
    match f'{CAT}_{YEAR}':
        case 'GF_2018':
            goodjets = '(Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_btagDeepCvL > -1)'
        case _:
            raise ValueError(f'Good-jets definition does not exist for the combination of {CAT} and {YEAR}.')
    new_rdf = (rdf.Define('goodJets', goodjets)
                .Define('nGoodJets', 'Sum(goodJets)')
                .Filter('nGoodJets>=2', 'Events must contain two jets from the charm quarks.')
    )
    cut_flow_dict = add_cut_flow(cut_flow_dict, 'nGoodJets>=2', new_rdf.Sum('w').GetValue())
    new_rdf = (new_rdf.Define('goodJetPt', 'Jet_pt[goodJets]')
                .Define('goodJetEta', 'Jet_eta[goodJets]')
                .Define('goodJetPhi', 'Jet_phi[goodJets]')
                .Define('goodJetMass', 'Jet_mass[goodJets]')
                .Define('goodJetCvL', 'Jet_btagDeepCvL[goodJets]')
    )
    new_rdf = (new_rdf.Define('jet1_pt', 'goodJetPt[0]')
                .Define('jet2_pt', 'goodJetPt[1]')
                .Define('jet1_eta', 'goodJetEta[0]')
                .Define('jet2_eta', 'goodJetEta[1]')
                .Define('jet1_phi', 'goodJetPhi[0]')
                .Define('jet2_phi', 'goodJetPhi[1]')
                .Define('jet1_mass', 'goodJetMass[0]')
                .Define('jet2_mass', 'goodJetMass[1]')
                .Define('jet1_CvL', 'goodJetCvL[0]')
                .Define('jet2_CvL', 'goodJetCvL[1]')
                .Define('index_CloseFar', 'jetCloseFar(goodJetPt, goodJetEta, goodJetPhi, goodJetMass,'
                                        'Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)')
                .Filter('index_CloseFar[0]!= -1', 'at least one close jet')
                .Define('jetClose_nMuons','Jet_nMuons[goodJets][index_CloseFar[0]]')
                .Define('jetFar_nMuons','Jet_nMuons[goodJets][index_CloseFar[1]]')
                .Define('jetClose_nElectrons','Jet_nElectrons[goodJets][index_CloseFar[0]]')
                .Define('jetFar_nElectrons','Jet_nElectrons[goodJets][index_CloseFar[1]]')
                .Define('jetClose_pt','goodJetPt[index_CloseFar[0]]')
                .Define('jetFar_pt','goodJetPt[index_CloseFar[1]]')
                .Define('jetClose_cRegCorr','Jet_cRegCorr[goodJets][index_CloseFar[0]]')
                .Define('jetFar_cRegCorr','Jet_cRegCorr[goodJets][index_CloseFar[1]]')
                .Define('jetClose_CvL','goodJetCvL[index_CloseFar[0]]')
                .Define('jetFar_CvL','goodJetCvL[index_CloseFar[1]]')
                .Define('jetClose_nConst','Jet_nConstituents[goodJets][index_CloseFar[0]]')
                .Define('jetFar_nConst','Jet_nConstituents[goodJets][index_CloseFar[1]]')
                .Define('jetClose_JpsiPtRatio','Jpsi_kin_pt[0]/goodJetPt[index_CloseFar[0]]')
    )
    cut_flow_dict = add_cut_flow(cut_flow_dict, 'nCloseJets>0', new_rdf.Sum('w').GetValue())
    branches = ['goodJets', 'nGoodJets', 'goodJetPt', 'goodJetEta', 'goodJetPhi', 'goodJetMass', 'goodJetCvL']
    branches += get_rdf_branches('jet')
    if not data:
        new_rdf = (new_rdf.Define('goodPartonFlavour', 'Jet_partonFlavour[goodJets]')
                    .Define('jet1_partonFlavour', 'abs(goodPartonFlavour[0])')
                    .Define('jet2_partonFlavour', 'abs(goodPartonFlavour[1])')
                    .Define('jetFar_ptCharm','(abs(goodPartonFlavour[index_CloseFar[1]])==4) ? Jet_cRegCorr[goodJets][index_CloseFar[1]] * goodJetPt[index_CloseFar[1]]:-1.')
                    .Define('jetFar_ptGluon','(abs(goodPartonFlavour[index_CloseFar[1]])==21) ? Jet_cRegCorr[goodJets][index_CloseFar[1]] * goodJetPt[index_CloseFar[1]]:-1.')
                    .Define('jetFar_CvLCharm','(abs(goodPartonFlavour[index_CloseFar[1]])==4) ? goodJetCvL[index_CloseFar[1]]:-1.')
                    .Define('jetFar_CvLGluon','(abs(goodPartonFlavour[index_CloseFar[1]])==21) ? goodJetCvL[index_CloseFar[1]]:-1.')
                    .Define('jetClose_partonFlavour','abs(goodPartonFlavour[index_CloseFar[0]])')
                    .Define('jetFar_partonFlavour','abs(goodPartonFlavour[index_CloseFar[1]])')
                    .Define('jetClose_Scale','goodJetPt[index_CloseFar[0]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]')
                    .Define('jetFar_Scale','goodJetPt[index_CloseFar[1]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]')
                    .Define('jetClose_cRegCorrScale','goodJetPt[index_CloseFar[0]]*Jet_cRegCorr[goodJets][index_CloseFar[0]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]')
                    .Define('jetFar_cRegCorrScale','goodJetPt[index_CloseFar[1]]*Jet_cRegCorr[goodJets][index_CloseFar[1]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]')     
        )
        branches += ['goodPartonFlavour']
        branches += get_rdf_branches('jet_mconly')
    return new_rdf, branches, cut_flow_dict

def rdf_def_higgs(rdf, branches=[], cut_flow_dict={}):
    '''Higgs RDF definition.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        branches (list(str), optional): List of TBranches.
            Defaults to [].
        cut_flow_dict (dict(str: float), optional): The dictionary containing the cut flow.
            Defaults to {}.

    Returns:
        new_rdf (ROOT.RDataFrame): Modified ROOT.RDataFrame.
        branches (list(str)): Updated list of TBranches.
        cut_flow_dict (dict(str: float)): Updated dictionary containing the cut flow.
    '''
    load_functions('reco')
    new_rdf = rdf.Define('higgs_mass_corr',
                        ('Minv3massiveCorr(goodJetPt, goodJetEta, goodJetPhi, goodJetMass,' +\
                        'Jet_cRegCorr[goodJets], Jet_muonIdx1[goodJets][index_CloseFar[0]],' +\
                        'Jet_muonIdx2[goodJets][index_CloseFar[0]],' +\
                        'Muon_pt, Muon_eta, Muon_phi,' +\
                        'Jpsi_muon1_pt, Jpsi_muon1_eta, Jpsi_muon1_phi,' +\
                        'Jpsi_muon2_pt, Jpsi_muon2_eta, Jpsi_muon2_phi,' +\
                        'Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0])'))
    branches += get_rdf_branches('higgs')
    return new_rdf, branches, cut_flow_dict

def rdf_def_genpart(rdf, branches=[]): # TODO: rename branches so that they start with 'gen'
    '''RDF define gen particles.

    Args:
        rdf (ROOT.RDataFrame): The ROOT.RDataFrame object.
        branches (list(str), optional): List of TBranches.
            Defaults to [].

    Returns:
        new_rdf (ROOT.RDataFrame): Modified RDataFrame.
        branches (list(str)): Updated list of TBranches.
    '''
    load_functions(mode='gen')

    # Find Higgs
    new_rdf = (rdf.Define('GenPart_Higgs_idx',
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
                # Find Higgs Daughter
                .Define('GenPart_HiggsDaughters_pdgId',
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
                        'SelectByIdx(GenPart_px, GenPart_HiggsDaughters_idx)')
                .Define('HiggsDaughters_py',
                        'SelectByIdx(GenPart_py, GenPart_HiggsDaughters_idx)')
                .Define('HiggsDaughters_pz',
                        'SelectByIdx(GenPart_pz, GenPart_HiggsDaughters_idx)')   
    )

    # Get muons
    new_rdf = (new_rdf.Define('muminus_index', 'IndexFindPDG(GenPart_HiggsDaughters_pdgId, 13)') # muon
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
    new_rdf = (new_rdf.Define('dR_muminus_muplus', 'DeltaR(muminus_eta, muminus_phi, muplus_eta, muplus_phi)')
                    .Define('deta_muminus_muplus', 'muminus_eta - muplus_eta')
                    .Define('dphi_muminus_muplus', 'muminus_phi - muplus_phi')
                    .Define('dpt_muminus_muplus', 'muminus_pt - muplus_pt')
                    .Define('dE_muminus_muplus', 'muminus_energy - muplus_energy')
    )

    # Get J/Psi
    new_rdf = (new_rdf.Define('JPsi_cand',
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
    new_rdf = (new_rdf.Define('charm_index', 'IndexFindPDG(GenPart_HiggsDaughters_pdgId, 4)')  # charm
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
    new_rdf = (new_rdf.Define('dR_charm_anticharm', 'DeltaR(charm_eta, charm_phi, anticharm_eta, anticharm_phi)')
                    .Define('deta_charm_anticharm', 'charm_eta - anticharm_eta')
                    .Define('dphi_charm_anticharm', 'charm_phi - anticharm_phi')
                    .Define('dpt_charm_anticharm', 'charm_pt - anticharm_pt')
                    .Define('dpt_leadcharm_subcharm', 'leadcharm_pt - subcharm_pt')
    )

    # Get di-charm
    new_rdf = (new_rdf.Define('dicharm_cand',
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
    new_rdf = (new_rdf.Define('dR_charm_JPsi', 'DeltaR(charm_eta, charm_phi, JPsi_cand_eta, JPsi_cand_phi)')
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
    new_rdf = (new_rdf.Define('dR_dicharm_JPsi', 'DeltaR(dicharm_cand_eta, dicharm_cand_phi, JPsi_cand_eta, JPsi_cand_phi)')
                    .Define('deta_dicharm_JPsi', 'dicharm_cand_eta - JPsi_cand_eta')
                    .Define('dphi_dicharm_JPsi', 'dicharm_cand_phi - JPsi_cand_phi')
                    .Define('dpt_dicharm_JPsi', 'dicharm_cand_pt - JPsi_cand_pt')
    )

    # Get Higgs
    new_rdf = (new_rdf.Define('Higgs_cand',
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

    # Branches
    branches += ['GenPart_Higgs_idx',
                'Higgs_energy',
                'Higgs_eta',
                'Higgs_mass',
                'Higgs_phi',
                'Higgs_pt',
                'Higgs_px',
                'Higgs_py',
                'Higgs_pz']
    
    branches += ['GenPart_HiggsDaughters_pdgId',
                'GenPart_HiggsDaughters_idx',
                'GenPart_HiggsGrandDaughters_pdgId',
                'GenPart_HiggsGrandDaughters_idx',
                'HiggsDaughters_energy',
                'HiggsDaughters_eta',
                'HiggsDaughters_mass',
                'HiggsDaughters_phi',
                'HiggsDaughters_pt',
                'HiggsDaughters_px',
                'HiggsDaughters_py',
                'HiggsDaughters_pz']
    
    branches += ['JPsi_cand_mass',
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
                'Higgs_cand_energy']
    
    branches += ['dR_muminus_muplus',
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
                'dpt_dicharm_JPsi']
    
    branches += ['muminus_mass',
                'muminus_energy',
                'muminus_phi',
                'muminus_eta',
                'muminus_pt',
                'muminus_px',
                'muminus_py',
                'muminus_pz',
                'muplus_mass',
                'muplus_energy',
                'muplus_phi',
                'muplus_eta',
                'muplus_pt',
                'muplus_px',
                'muplus_py',
                'muplus_pz',
                'charm_mass',
                'charm_energy',
                'charm_phi',
                'charm_eta',
                'charm_pt',
                'charm_px',
                'charm_py',
                'charm_pz',
                'anticharm_mass',
                'anticharm_energy',
                'anticharm_phi',
                'anticharm_eta',
                'anticharm_pt',
                'anticharm_px',
                'anticharm_py',
                'anticharm_pz',
                'leadcharm_pdgId',
                'leadcharm_mass',
                'leadcharm_energy',
                'leadcharm_phi',
                'leadcharm_eta',
                'leadcharm_pt',
                'leadcharm_px',
                'leadcharm_py',
                'leadcharm_pz',
                'subcharm_pdgId',
                'subcharm_mass',
                'subcharm_energy',
                'subcharm_phi',
                'subcharm_eta',
                'subcharm_pt',
                'subcharm_px',
                'subcharm_py',
                'subcharm_pz']
    
    return new_rdf, branches