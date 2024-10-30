import os
from jpsicc import limits

def fitSigBkgIndiv(YEAR, VERS, CAT, CMSSW, weights,
                   samp_dict, var_min, var_max, SR_low, SR_high):
    '''Fit signal and background individually and save to a workspace. Plots are created.

    Args:
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        weights (bool): Whether weights were used.
        suffix (str, optional): Suffix of the workspace file name.
        samp_dict (dict): Dictionary containing the configuration for each sample.
            It must follow the following format.
            samp_dict[samp_key] = {'SAMP': SAMP,
                                   'filename': filename,
                                   'treename': treename,
                                   'col_name': col_name,
                                   'bkg_or_sig': 'bkg' or 'sig',
                                   'decay_mother': 'Z', 'H', or 'both'.
                                   'pdf_type': pdf_type}
        var_min (int or float): Lower bound of the variable, also of the CR.
        var_max (int or float): Upper bound of the variable, also of the CR.
        SR_low (int or float): Lower bound of the SR.
        SR_high (int or float: Upper bound of the SR.

    Returns:
        (None)

    Raises:
        KeyError: If samp_dict is not in the proper format.
    '''
    nbins = int(var_max - var_min)
    rwc = limits.RooWorkspaceCreator(YEAR, VERS, CAT, CMSSW, weights, f'{var_min:.0f}_{var_max:.0f}')
    rwc.addVar(var_name='m_mumucc', var_title='m_{#mu^{+}#mu^{-}c#bar{c}}', value=125, var_min=var_min, var_max=var_max, unit='GeV/c^{2}')
    rwc.defineRegions(SR_low=SR_low, SR_high=SR_high, CR_low=var_min, CR_high=var_max)
    for val in samp_dict.values():
        rwc.addDataHist(SAMP=val['SAMP'], filename=val['filename'], treename=val['treename'], var_name='m_mumucc', col_name=val['col_name'], nbins=nbins, weight_name='w')
        rwc.addPDF(SAMP=val['SAMP'], pdf_type=val['pdf_type'], var_name='m_mumucc')
    for decay_mother in ('H', 'Z'):
        var_savename = f'm{decay_mother}'
        for val in samp_dict.values():
            if val['decay_mother'] in (decay_mother, 'both'):
                rwc.configPlot(SAMP=val['SAMP'], preset=f'{val["bkg_or_sig"]}_pdf', pdf_type=val['pdf_type'])
                rwc.configPlot(SAMP=val['SAMP'], preset='data')
            rwc.makePlot(var_name='m_mumucc', plot_name=f'{val["SAMP"]}__{val["pdf_type"]}__{var_savename}_{var_min:.0f}_{var_max:.0f}', plot_title='H#rightarrow J/#psi + c#bar{c}')
    return

if __name__=='__main__':
    snapshot_dir = '/work/submit/mariadlf/Hrare_JPsiCC/OCT25/'
    samp_dict = {}
    samp_dict['bkg_JpsiToMuMu'] =\
        {'SAMP': 'MC_BKG_JpsiToMuMu',
        'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_10_2018.root'),
        'treename': 'events',
        'col_name': 'massHiggsCorr',
        'bkg_or_sig': 'bkg',
        'decay_mother': 'both',
        'pdf_type': 'gaussian_X_exponential'}
    samp_dict['bkg_BToJpsi'] =\
        {'SAMP': 'MC_BKG_BToJpsi',
        'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_11_2018.root'),
        'treename': 'events',
        'col_name': 'massHiggsCorr',
        'bkg_or_sig': 'bkg',
        'decay_mother': 'both',
        'pdf_type': 'gaussian_X_exponential'}
    samp_dict['sig_H'] =\
        {'SAMP': 'MC_SIG_H',
        'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_1000_2018.root'),
        'treename': 'events',
        'col_name': 'massHiggsCorr',
        'bkg_or_sig': 'sig',
        'decay_mother': 'H',
        'pdf_type': 'gaussian'}
    samp_dict['sig_Z'] =\
        {'SAMP': 'MC_SIG_Z',
        'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_1001_2018.root'),
        'treename': 'events',
        'col_name': 'massHiggsCorr',
        'bkg_or_sig': 'sig',
        'decay_mother': 'Z',
        'pdf_type': 'gaussian'}
    
    fitSigBkgIndiv(2018, '202410', 'GF', 'ROOT_6_33_01', weights=True,
                   samp_dict=samp_dict, var_min=0, var_max=200, SR_low=60, SR_high=160)

    # var_min, var_max = 90, 200
    # nbins = int(var_max-var_min)
    # rwc = RooWorkspaceCreator(2018, '202407', 'GF', 'ROOT_6_33_01', suffix='90-200')
    # rwc.addVar(var_name='m_mumucc', var_title='m_{#mu^{+}#mu^{-}c#bar{c}}', value=125, var_min=var_min, var_max=var_max, unit='GeV/c^{2}')
    # rwc.defineRegions(SR_low=110., SR_high=130., CR_low=var_min, CR_high=var_max)
    # rwc.addDataHist(SAMP='MC_SIG_H', filename=fname_sig_H, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addDataHist(SAMP='MC_SIG_Z', filename=fname_sig_Z, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addDataHist(SAMP='MC_BKG_JpsiToMuMu', filename=fname_bkg_JpsiToMuMu, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addDataHist(SAMP='MC_BKG_BToJpsi', filename=fname_bkg_JpsiToMuMu, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addPDF(SAMP='MC_SIG_H', pdf_type='gaussian', var_name='m_mumucc')
    # rwc.addPDF(SAMP='MC_SIG_Z', pdf_type='gaussian', var_name='m_mumucc')
    # rwc.addPDF(SAMP='MC_BKG_JpsiToMuMu', pdf_type='exponential', var_name='m_mumucc')
    # rwc.addPDF(SAMP='MC_BKG_BToJpsi', pdf_type='exponential', var_name='m_mumucc')

    # rwc.configPlot(SAMP='MC_BKG_JpsiToMuMu', preset='bkg_pdf', pdf_type='exponential')
    # rwc.configPlot(SAMP='MC_BKG_JpsiToMuMu', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_BKG_JpsiToMuMu__exponential__mH_{var_min}_{var_max}', plot_title='H#rightarrow J/#psi + c#bar{c}')

    # rwc.clearConfigPlot()
    # rwc.configPlot(SAMP='MC_BKG_BToJpsi', preset='bkg_pdf', pdf_type='exponential')
    # rwc.configPlot(SAMP='MC_BKG_BToJpsi', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_BKG_BToJpsi__exponential__mH_{var_min}_{var_max}', plot_title='H#rightarrow J/#psi + c#bar{c}')

    # rwc.clearConfigPlot()
    # rwc.configPlot(SAMP='MC_SIG_H', preset='sig_pdf', pdf_type='gaussian')
    # rwc.configPlot(SAMP='MC_SIG_H', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_SIG_H__gaussian__mH_{var_min}_{var_max}', plot_title='H#rightarrow J/#psi + c#bar{c}')

    # rwc.clearConfigPlot()
    # rwc.configPlot(SAMP='MC_SIG_Z', preset='sig_pdf', pdf_type='gaussian')
    # rwc.configPlot(SAMP='MC_SIG_Z', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_SIG_Z__gaussian__mZ_{var_min}_{var_max}', plot_title='Z#rightarrow J/#psi + c#bar{c}')


    # var_min, var_max = 0, 200
    # nbins = int(var_max-var_min)
    # rwc = RooWorkspaceCreator(2018, '202407', 'GF', 'ROOT_6_33_01', suffix="0-200")
    # rwc.addVar(var_name='m_mumucc', var_title='m_{#mu^{+}#mu^{-}c#bar{c}}', value=125, var_min=var_min, var_max=var_max, unit='GeV/c^{2}')
    # rwc.defineRegions(SR_low=70., SR_high=140., CR_low=var_min, CR_high=var_max)
    # rwc.addDataHist(SAMP='MC_SIG_H', filename=fname_sig_H, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addDataHist(SAMP='MC_SIG_Z', filename=fname_sig_Z, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addDataHist(SAMP='MC_BKG_JpsiToMuMu', filename=fname_bkg_JpsiToMuMu, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addDataHist(SAMP='MC_BKG_BToJpsi', filename=fname_bkg_JpsiToMuMu, treename='events', var_name='m_mumucc', col_name='massHiggsCorr', nbins=nbins, weight_name='w')
    # rwc.addPDF(SAMP='MC_SIG_H', pdf_type='gaussian', var_name='m_mumucc')
    # rwc.addPDF(SAMP='MC_SIG_Z', pdf_type='gaussian', var_name='m_mumucc')
    # rwc.addPDF(SAMP='MC_BKG_JpsiToMuMu', pdf_type='gaussian_X_exponential', var_name='m_mumucc')
    # rwc.addPDF(SAMP='MC_BKG_BToJpsi', pdf_type='gaussian_X_exponential', var_name='m_mumucc')

    # rwc.configPlot(SAMP='MC_BKG_JpsiToMuMu', preset='bkg_pdf', pdf_type='gaussian_X_exponential')
    # rwc.configPlot(SAMP='MC_BKG_JpsiToMuMu', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_BKG_JpsiToMuMu__gaussian_X_exponential__mH_{var_min: }_{var_max}', plot_title='H#rightarrow J/#psi + c#bar{c}')

    # rwc.clearConfigPlot()
    # rwc.configPlot(SAMP='MC_BKG_BToJpsi', preset='bkg_pdf', pdf_type='gaussian_X_exponential')
    # rwc.configPlot(SAMP='MC_BKG_BToJpsi', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_BKG_BToJpsi__gaussian_X_exponential__mH_{var_min}_{var_max}', plot_title='H#rightarrow J/#psi + c#bar{c}')

    # rwc.clearConfigPlot()
    # rwc.configPlot(SAMP='MC_SIG_H', preset='sig_pdf', pdf_type='gaussian')
    # rwc.configPlot(SAMP='MC_SIG_H', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_SIG_H__gaussian__mH_{var_min}_{var_max}', plot_title='H#rightarrow J/#psi + c#bar{c}')

    # rwc.clearConfigPlot()
    # rwc.configPlot(SAMP='MC_SIG_Z', preset='sig_pdf', pdf_type='gaussian')
    # rwc.configPlot(SAMP='MC_SIG_Z', preset='data')
    # rwc.makePlot(var_name='m_mumucc', plot_name=f'MC_SIG_Z__gaussian__mZ_{var_min}_{var_max}', plot_title='Z#rightarrow J/#psi + c#bar{c}')