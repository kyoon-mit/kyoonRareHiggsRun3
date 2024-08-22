import os
import ROOT
from jpsicc import analyzer as an

def applyTriggers(SAMP, YEAR, VERS, CAT, file_name, trigger_list, hist_names):
    '''Apply triggers and make histograms.

    Args:
        SAMP (str): Either one of the following options.
            'DATA_BKG', 'MC_BKG', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        file_name (str): Name of the ROOT file which contains the events along
            with the columns of triggers.
        trigger_list (list(str)): List of triggers to apply.
        hist_names (list(str)): List of TH1 names to create.

    Raises:
        TypeError: If file_names and hist_names are not (list(str)).

    Returns:
        (None)
    '''
    if not isinstance(file_name, str)\
        and not all(isinstance(h, str) for h in hist_names):
        raise TypeError('The value provided for file_names or hist_names is not list(str).')
    loader = an.JPsiCCLoader(SAMP, YEAR, VERS, CAT)
    for i in range(len(trigger_list)):
        col_name = f'trigger_{i+1:02d}'
        loader.readSnapshot(file_name)
        loader.cut(f'{col_name}>0', define_col=col_name, define_exp=trigger_list[i])
        loader.makeHistos(genplots=True, cdf=True, normSIG=False, save=True, draw=True,
                          keys=hist_names, user_sfx=trigger_list[i].replace(' || ', '_OR_'))

def makeEfficiencyHisto(file_names, hist_names):
    '''Make efficiency histograms.

    The first index in the arguments are taken as the 'total'.

    Args:
        file_names (list(str)): List of ROOT file names.
        hist_names (list(str)): List of TH1 names in each of the file names.
    
    Raises:
        TypeError: If file_names and hist_names are not (list(str)).
        IndexError: If the lengths of file_names and hist_names are different.
        ValueError: If the lengths are less than 2.

    Returns:
        eff_list (list(ROOT.TEfficiency)): List of the TEfficiency objects.
    '''
    if not all(isinstance(f, str) for f in file_names)\
        and not all(isinstance(h, str) for h in hist_names):
        raise TypeError('The value provided for file_names or hist_names is not list(str).')
    elif len(file_names) != len(hist_names):
        raise IndexError('The lengths of file_names and hist_names must be the same.')
    elif len(file_names) < 2:
        raise ValueError('There must be at least two items in the file_names and hist_names each.')
    f_total = ROOT.TFile.Open(file_names[0], 'READ')
    h_total = f_total.Get(hist_names[0])
    h_total.SetDirectory(0)
    eff_list = []
    for i in range(1, len(file_names)):
        f_passed = ROOT.TFile.Open(file_names[i], 'READ')
        h_passed = f_passed.Get(hist_names[i])
        h_passed.SetDirectory(0)
        p_eff = ROOT.TEfficiency(h_passed, h_total)
        p_eff.SetDirectory(0)
        eff_list.append(p_eff)
    return eff_list

if __name__=='__main__':
    print('Choose preset (eff, hist): ', end='')
    preset = input()
    match preset:
        case 'eff':
            filedirectory = '/work/submit/kyoon/CMSSW_13_3_0/src/RareHiggsRun3/JPsiCC/analysis'
            file_names = ['histos_MC_SIG_2018_GF_v202407_20240822_WEIGHT_no_filter.root',
                        'histos_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon25_Jpsi.root',
                        'histos_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon20_Jpsi_Barrel_Seagulls.root',
                        'histos_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon25_Jpsi_OR_HLT_Dimuon20_Jpsi_Barrel_Seagulls.root',
                        'histos_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon0_Jpsi.root',
                        'histos_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Mu8.root']
            file_names = [os.path.join(filedirectory, f) for f in file_names]
            hist_names = ['Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240822_WEIGHT_no_filter_cumulative',
                        'Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon25_Jpsi_cumulative',
                        'Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon20_Jpsi_Barrel_Seagulls_cumulative',
                        'Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon25_Jpsi_OR_HLT_Dimuon20_Jpsi_Barrel_Seagulls_cumulative',
                        'Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Dimuon0_Jpsi_cumulative',
                        'Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240822_WEIGHT_HLT_Mu8_cumulative']
            eff_list = makeEfficiencyHisto(file_names, hist_names)
            for i in range(len(hist_names)-1):
                c = ROOT.TCanvas()
                fname = f'trigger_efficiency_{hist_names[i]}'
                eff_list[i].SetMarkerColor(ROOT.kBlue)
                eff_list[i].SetLineColor(ROOT.kBlue)
                eff_list[i].SetMarkerStyle(8)
                eff_list[i].SetMarkerSize(0.5)
                eff_list[i].SetLineWidth(0)
                eff_list[i].Draw('P')
                c.SetTitle(fname) # TODO: How do I set the title?
                c.SaveAs(os.path.join(filedirectory, f'{fname}.png'))
                c.Close()
        
        case 'hist':
            applyTriggers(SAMP='MC_SIG', YEAR=2018, VERS='202407', CAT='GF',
                        file_name='snapshot_MC_SIG_2018_GF_v202407_20240822_WEIGHT_no_filter.root',
                        trigger_list=['HLT_Dimuon25_Jpsi',
                                        # 'HLT_Dimuon25_Jpsi_noCorrL1',
                                        'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
                                        'HLT_Dimuon25_Jpsi || HLT_Dimuon20_Jpsi_Barrel_Seagulls',
                                        # 'HLT_Dimuon25_Jpsi || HLT_Mu8',
                                        # 'HLT_Dimuon25_Jpsi || HLT_Dimuon0_Jpsi',
                                        'HLT_Mu8',
                                        # 'HLT_Mu12',
                                        # 'HLT_Mu15',
                                        # 'HLT_Mu17',
                                        # 'HLT_Mu19',
                                        # 'HLT_Mu8 || HLT_Dimuon0_Jpsi',
                                        'HLT_Dimuon0_Jpsi'],
                                        # 'HLT_Dimuon0_Jpsi_NoVertexing',
                                        # 'HLT_Dimuon0_Jpsi_L1_NoOS',
                                        # 'HLT_Dimuon0_Jpsi_NoVertexing_NoOS',
                                        # 'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05',
                                        # 'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05 || HLT_DoubleMu4_Jpsi_NoVertexing',
                                        # 'HLT_DoubleMu4_Jpsi_NoVertexing'],
                        hist_names=['nMuon', 'Jpsi_kin_pt', 'Jpsi_kin_mass', 'Jpsi_kin_eta'
                                    'gen_muplus_pt', 'gen_muminus_pt', 'gen_muplus_eta', 'gen_muminus_eta'])