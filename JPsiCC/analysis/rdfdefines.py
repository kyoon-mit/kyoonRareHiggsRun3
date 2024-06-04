'''
rdfdefines
+++++++++++++++++++++++++++++++++++++++++++++

Library of definitions that contain analysis category-specific RDF defines and filters.
'''
import ROOT

def load_functions():
    ROOT.gSystem.AddDynamicPath("./.")
    ROOT.gROOT.ProcessLine(".include ./.")
    ROOT.gInterpreter.AddIncludePath("./.")
    ROOT.gInterpreter.Declare('#include "./config/functions.h"')

def compute_weights(rdf, lumi, data):
    '''Compute the weights.

    Args:
        rdf (RDataFrame): The RDataFrame object.
        lumi (float): Luminosity.
        data (bool): Whether to use DATA.

    Returns:
        weight (str): Weight expression.
    '''
    # gen_event_sumw = rdf.Sum('genEventSumw').GetValue()
    # gen_event_count = rdf.Sum('genEventCount').GenValue()
    # weight = 1./gen_event_sumw
    # print(f'{gen_event_sumw}, {gen_event_count}, {weight}')
    if not data:
        weight = f'{lumi}*genWeight*100000'
    else:
        weight = 1
    return weight

def rdf_def_generic(rdf, lumi, data):
    '''Generic rdf definition.

    Args:
        rdf (RDataFrame): The RDataFrame object.
        lumi (float): Luminosity.
        data (bool): Whether to use DATA.

    Returns:
        new_rdf (RDataFrame): Modified RDataFrame.
    '''
    load_functions()
    TRIGGER = 'HLT_Dimuon20_Jpsi_Barrel_Seagulls or HLT_Dimuon25_Jpsi'
    GOODJETS = '(Jet_pt>20 and abs(Jet_eta)<2.5 and Jet_btagDeepCvL>-1)'
    weight = compute_weights(rdf, lumi, data)
    new_rdf= (rdf.Filter('nJpsi>0', 'at least one Jpsi')
                .Define('w', f'{weight}')
                # .Define('triggerAna', f'{TRIGGER}')
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