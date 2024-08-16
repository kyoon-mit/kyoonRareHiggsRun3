import os
import ROOT

def makeEfficiencyHisto(filenames, histonames):
    '''Make efficiency histograms.

    The first index in the arguments are taken as the 'total'.

    Args:
        filenames (list(str)): List of ROOT file names.
        histonames (list(str)): List of TH1 names corresponding to the filenames.
    
    Raises:
        TypeError: If filenames and histonames are not (list(str)).
        IndexError: If the lengths of filenames and histonames are different.
        ValueError: If the lengths are less than 2.

    Returns:
        eff_list (list(ROOT.TEfficiency)): List of the TEfficiency objects.
    '''
    if not all(isinstance(f, str) for f in filenames)\
        and not all(isinstance(h, str) for h in histonames):
        raise TypeError('The value provided for filenames or histonames is not list(str).')
    elif len(filenames) != len(histonames):
        raise IndexError('The lengths of filenames and histonames must be the same.')
    elif len(filenames) < 2:
        raise ValueError('There must be at least two items in the filenames and histonames each.')
    f_total = ROOT.TFile.Open(filenames[0], 'READ')
    h_total = f_total.Get(histonames[0])
    h_total.SetDirectory(0)
    eff_list = []
    for i in range(1, len(filenames)):
        f_passed = ROOT.TFile.Open(filenames[i], 'READ')
        h_passed = f_passed.Get(histonames[i])
        h_passed.SetDirectory(0)
        p_eff = ROOT.TEfficiency(h_passed, h_total)
        p_eff.SetDirectory(0)
        eff_list.append(p_eff)
    return eff_list

if __name__=='__main__':
    filedirectory = '/work/submit/kyoon/CMSSW_13_3_0/src/RareHiggsRun3/JPsiCC/analysis'
    filenames = ['histos_MC_SIG_2018_GF_v202407_20240816_WEIGHT_no_filter.root',
                 'histos_MC_SIG_2018_GF_v202407_20240816_WEIGHT_HLT_Dimuon25_Jpsi.root']
    filenames = [os.path.join(filedirectory, f) for f in filenames]
    histonames = ['Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240816_WEIGHT_no_filter_cumulative',
                  'Jpsi_kin_pt_MC_SIG_2018_GF_v202407_20240816_WEIGHT_HLT_Dimuon25_Jpsi_cumulative']
    eff_list = makeEfficiencyHisto(filenames, histonames)
    c = ROOT.TCanvas()
    eff_list[0].SetMarkerColor(ROOT.kBlue)
    eff_list[0].SetLineColor(ROOT.kBlue)
    eff_list[0].SetMarkerStyle(8)
    eff_list[0].SetMarkerSize(0.5)
    eff_list[0].SetLineWidth(0)
    eff_list[0].Draw('P')
    c.SaveAs(os.path.join(filedirectory, 'test.png'))