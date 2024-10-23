from jpsicc import analyzer as an

if __name__=='__main__':
    print('Choose preset (mcbkg, databkg, mcsig, allsnap, gen, allplots, noweight, cut): ', end='')
    preset = input()
    match preset:
        case 'mcbkg':
            trig = 'HLT_Dimuon0_Jpsi'
            mcbkg = an.JPsiCCLoader('MC_BKG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                                   filter_jpsi=False, filter_jets=False, filter_higgs=False)
            mcbkg.cut('trigger_05>0',
                      define_col='trigger_05',
                      define_exp=trig)
            mcbkg.snapshotRDF(user_sfx=trig)
            # mcbkg.readSnapshot(f'snapshot_MC_BKG_2018_GF_v202407_20240823_WEIGHT_{trig}.root')
            mcbkg.makeHistos(genplots=True, cdf=True, save=True, draw=True, normSIG=False,
                             keys=['nMuon', 'Jpsi_kin_pt', 'Jpsi_kin_mass', 'Jpsi_kin_eta', 'Muon_pt'],
                                #    'gen_muplus_pt', 'gen_muminus_pt', 'gen_muplus_eta', 'gen_muminus_eta'],
                             user_sfx=trig)

        case 'databkg':
            databkg = an.JPsiCCLoader('DATA_BKG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
            databkg.createDataRDF('event_spec_data_bkg_skim.json')
            databkg.defineColumnsRDF()
            databkg.snapshotRDF()
            databkg.makeCutFlowPlot()

        case 'mcsig':
            # verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                                   filter_jpsi=False, filter_jets=False, filter_higgs=False)
            mcsig.snapshotRDF(user_sfx='no_filter')
            mcsig.readSnapshot('snapshot_MC_SIG_2018_GF_v202407_20241023_WEIGHT_no_filter.root')
            # # mcsig.cut('nJpsi > 0')
            mcsig.makeHistos(genplots=False, cdf=True, save=True, draw=True, normSIG=False,
                             user_sfx='no_filter')
            # mcsig.snapshotRDF('snapshot_MC_SIG_2018_GF_v202407_20240821_WEIGHT_HLT_Dimuon20_Jpsi_Barrel_Seagulls.root')
        
        case 'allsnap':
            mcbkg = an.JPsiCCLoader('MC_BKG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF()
            mcbkg.snapshotRDF()
            mcbkg.makeCutFlowPlot()
            
            databkg = an.JPsiCCLoader('DATA_BKG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
            databkg.createDataRDF('event_spec_data_bkg_skim.json')
            databkg.defineColumnsRDF()
            databkg.snapshotRDF()
            databkg.makeCutFlowPlot()
            
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF()
            mcsig.snapshotRDF()
            mcsig.makeCutFlowPlot()

        case 'gen':
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
            mcsig.readSnapshot('snapshot_MC_SIG_2018_GF_v202407_20240816_WEIGHT_no_filter.root')
            mcsig.cut('nMuon<2')
            mcsig.cut('gen_muplus_pt>3.5')
            mcsig.cut('gen_muminus_pt>3.5')
            # mcsig.cut('gen_muplus_eta>-2.4 && gen_muplus_eta<2.4')
            # mcsig.cut('gen_muminus_eta>-2.4 && gen_muminus_eta<2.4')
            mcsig.makeHistos(genplots=True, plot=True, draw=True, normSIG=False, cdf=False, save=True, user_sfx='nmuon<2+gen_muplus_pt>3.5+gen_muminus_pt>3.5')

        case 'allplots' | 'noweight' | 'cut':
            if preset=='noweight': use_weight, draw_indiv = False, False
            else: use_weight, draw_indiv = True, True
            
            mcsig = an.JPsiCCAnalyzer('MC_SIG', 2018, '202407', 'GF', 'CMSSW_13_3_0', weights=use_weight)
            # an.readSnapshotSAMP('MC_BKG1', 'snapshot_MC_BKG_2018_GF_v202406_20240802_WEIGHT.root',
            #                     sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            # an.readSnapshotSAMP('MC_BKG2', 'snapshot_MC_BKG_2018_GF_v202406_20240802_WEIGHT.root',
            #                     sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            # an.readSnapshotSAMP('DATA_BKG', 'snapshot_DATA_BKG_2018_GF_v202406_20240802_WEIGHT.root')
            mcsig.readSnapshotSAMP('MC_SIG', 'snapshot_MC_SIG_2018_GF_v202407_20240802_WEIGHT.root')

            if preset != 'cut': mcsig.stackMultiHistos(draw_indiv=draw_indiv)
            else:
                mcsig.addCutOptions('Jpsi_kin_pt[0]>5') # TODO: make function to detect & flatten RVec
                mcsig.addCutOptions('Jpsi_kin_pt[0]>10')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>15')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>20')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>25')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>30')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>30')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>35')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>40')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>45')
                mcsig.addCutOptions('Jpsi_kin_pt[0]>50')
                mcsig.addCutOptions('Jpsi_iso[0]>0.5')
                mcsig.addCutOptions('Jpsi_iso[0]>0.6')
                mcsig.addCutOptions('Jpsi_iso[0]>0.7')
                mcsig.addCutOptions('Jpsi_iso[0]>0.8')
                mcsig.addCutOptions('Jpsi_kin_sipPV[0]<2')
                mcsig.addCutOptions('Jpsi_kin_sipPV[0]<1.75')
                mcsig.addCutOptions('Jpsi_kin_sipPV[0]<1.5')
                mcsig.addCutOptions('Jpsi_kin_sipPV[0]<1.25')
                mcsig.makeCutTable(modify=False, plot=True)
                mcsig.makeCutFlowPlotAll()
                # mcsig.stackMultiHistos(draw_indiv=draw_indiv, user_sfx='cut')