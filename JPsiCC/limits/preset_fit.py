import os
from jpsicc import limits

def fitSigBkgIndiv(YEAR, VERS, CAT, CMSSW, weights,
                   samp_dict, var_min, var_max, SR_low, SR_high, suffix='',
                   sig_pdf_type='', bkg_pdf_type='', mean_low=None, mean_high=None,
                   sigma_low=None, sigma_high=None, step_low=None, step_high=None):
    '''Fit signal and background individually and save to a workspace. Plots are created.

    Args:
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        weights (bool): Whether weights were used.
        samp_dict (dict): Dictionary containing the configuration for each sample.
            It must follow the following format.
            samp_dict[samp_key] = {'SAMP': SAMP,
                                   'filename': filename,
                                   'treename': treename,
                                   'col_name': col_name,
                                   'bkg_or_sig': 'bkg' or 'sig',
                                   'decay_mother': 'Z', 'H', or 'both',
                                   'pdf_type': pdf_type,
                                   'strategy': fit strategy}
        var_min (int or float): Lower bound of the variable, also of the CR.
        var_max (int or float): Upper bound of the variable, also of the CR.
        SR_low (int or float): Lower bound of the SR.
        SR_high (int or float: Upper bound of the SR.
        suffix (str, optional): Suffix of the workspace file name.
            Defaults to ''.
        sig_pdf_type (str, optional): If specified, override pdf_type of sig in the samp_dict.
            Defaults to ''.
        bkg_pdf_type (str, optional): If specified, override pdf_type of bkg in the samp_dict.
            Defaults to ''.
        mean_low (float, optional): Lower range of the mean.
                Defaults to SR_low.
            mean_high (float, optional): Upper range of the mean.
                Defaults to SR_high.
            sigma_low (float, optional): Lower range of the Gaussian sigma.
                Defaults to SR_low.
            sigma_high (float, optional): Upper range of the Gaussian sigma.
                Defaults to SR_high.
            step_low (float, optional): Lower range of the turnon.
                Defaults to SR_low.
            step_high (float, optional): Upper range of the turnon.
                Defaults to SR_high.

    Returns:
        (None)

    Raises:
        KeyError: If samp_dict is not in the proper format.
    '''
    if mean_low is None: mean_low = SR_low
    if mean_high is None: mean_high = SR_high
    if sigma_low is None: sigma_low = SR_low
    if sigma_high is None: sigma_high = SR_high
    if step_low is None: step_low = SR_low
    if step_high is None: step_high = SR_high
    nbins = int(var_max - var_min)
    rwc = limits.RooWorkspaceCreator(YEAR, VERS, CAT, CMSSW, weights, suffix=f'{suffix}_m_mumucc_{var_min:.0f}_{var_max:.0f}')
    rwc.addVar(var_name='m_mumucc', var_title='m_{#mu^{+}#mu^{-}c#bar{c}}', value=125, var_min=var_min, var_max=var_max, unit='GeV/c^{2}')
    rwc.defineRegions(SR_low=SR_low, SR_high=SR_high, CR_low=var_min, CR_high=var_max)
    for val in samp_dict.values():
        if sig_pdf_type and val['bkg_or_sig']=='sig':
            val['pdf_type'] = sig_pdf_type
        elif bkg_pdf_type and val['bkg_or_sig']=='bkg':
            val['pdf_type'] = bkg_pdf_type
        if val['bkg_or_sig']=='bkg': mean_low, mean_high = 105, 115
        rwc.addDataHist(SAMP=val['SAMP'], filename=val['filename'], treename=val['treename'], var_name='m_mumucc', col_name=val['col_name'], nbins=nbins, weight_name='w')
        rwc.addPDF(SAMP=val['SAMP'], pdf_type=val['pdf_type'], var_name='m_mumucc', strategy=val['strategy'],
                   mean_low=mean_low, mean_high=mean_high, sigma_low=sigma_low, sigma_high=sigma_high,
                   step_low=step_low, step_high=step_high)
    for decay_mother in ('H', 'Z'):
        var_savename = f'm{decay_mother}'
        for val in samp_dict.values():
            if val['decay_mother'] in (decay_mother, 'both'):
                rwc.configPlot(SAMP=val['SAMP'], preset=f'{val["bkg_or_sig"]}_pdf', pdf_type=val['pdf_type'])
                rwc.configPlot(SAMP=val['SAMP'], preset='data')
                rwc.makePlot(var_name='m_mumucc', plot_name=f'{val["SAMP"]}__{val["pdf_type"]}__{var_savename}_{var_min:.0f}_{var_max:.0f}',
                             plot_title=f'{"Z" if CAT=="Z" else "H"} #rightarrow J/#psi + c#bar{{c}}')
    return

if __name__=='__main__':
    # snapshot_dir = '/work/submit/mariadlf/Hrare_JPsiCC/OCT25/'
    snapshot_dir = '/ceph/submit/data/user/k/kyoon/snapshots/mariadlf_Hrare_JPsiCC_20241025'
    bkg_list = ['gaussian_X_step_bernstein_3rd_order',
                'gaussian_X_step_bernstein_4th_order']
    # bkg_list = ['bernstein_4th_order']
    sig_list = ['crystal_ball']
    
    for bkg in bkg_list:
        samp_dict = {}
        samp_dict['bkg_JpsiToMuMu+BToJpsi'] =\
            {'SAMP': 'MC_BKG_JpsiToMuMu+BToJpsi',
            'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_bkg_comb_2018.root'),
            'treename': 'events',
            'col_name': 'massHiggsCorr',
            'bkg_or_sig': 'bkg',
            'decay_mother': 'both',
            'pdf_type': bkg,
            'strategy': 2}
        # samp_dict['bkg_JpsiToMuMu'] =\
        #     {'SAMP': 'MC_BKG_JpsiToMuMu',
        #     'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_10_2018.root'),
        #     'treename': 'events',
        #     'col_name': 'massHiggsCorr',
        #     'bkg_or_sig': 'bkg',
        #     'decay_mother': 'both',
        #     'pdf_type': bkg,
        #     'strategy': 2}
        # samp_dict['bkg_BToJpsi'] =\
        #     {'SAMP': 'MC_BKG_BToJpsi',
        #     'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_11_2018.root'),
        #     'treename': 'events',
        #     'col_name': 'massHiggsCorr',
        #     'bkg_or_sig': 'bkg',
        #     'decay_mother': 'both',
        #     'pdf_type': bkg,
        #     'strategy': 2}
        samp_dict_Z = samp_dict
        fitSigBkgIndiv(2018, '202410', 'GF', 'ROOT_6_33_01', weights=True, suffix=bkg,
                        samp_dict=samp_dict, var_min=55, var_max=160, SR_low=55, SR_high=160,
                        mean_low=85, mean_high=115, sigma_low=10, sigma_high=60,
                        step_low=95, step_high=105)
        # fitSigBkgIndiv(2018, '202410', 'Z', 'ROOT_6_33_01', weights=True,
        #             samp_dict=samp_dict_Z, var_min=55, var_max=160, SR_low=55, SR_high=120)

    # for sig in sig_list:
    #     samp_dict, samp_dict_Z = {}, {}
    #     samp_dict['sig_H'] =\
    #         {'SAMP': 'MC_SIG_H',
    #         'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_1000_2018.root'),
    #         'treename': 'events',
    #         'col_name': 'massHiggsCorr',
    #         'bkg_or_sig': 'sig',
    #         'decay_mother': 'H',
    #         'pdf_type': sig,
    #         'strategy': 2}
    #     samp_dict_Z['sig_Z'] =\
    #         {'SAMP': 'MC_SIG_Z',
    #         'filename': os.path.join(snapshot_dir, 'snapshotJpsiCC_1001_2018.root'),
    #         'treename': 'events',
    #         'col_name': 'massHiggsCorr',
    #         'bkg_or_sig': 'sig',
    #         'decay_mother': 'Z',
    #         'pdf_type': sig,
    #         'strategy': 2}
    #     fitSigBkgIndiv(2018, '202410', 'GF', 'ROOT_6_33_01', weights=True,
    #                 samp_dict=samp_dict, var_min=60, var_max=160, SR_low=60, SR_high=150)
    #     fitSigBkgIndiv(2018, '202410', 'Z', 'ROOT_6_33_01', weights=True,
    #                 samp_dict=samp_dict_Z, var_min=60, var_max=160, SR_low=60, SR_high=120)
    