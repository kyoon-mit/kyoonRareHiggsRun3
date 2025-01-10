from jpsicc import analyzer as an
from datetime import date

if __name__=='__main__':
    from sys import argv
    print(argv)
    print('Choose preset (full, stack, mcsig, mcbkg, data, allsnap, gen, allplots, noweight, cut): ', end='')
    # preset = input()
    preset = argv[1]
    print(preset)
    today = date.today()
    today_format = f'{today.year}{today.month:02}{today.day:02}'
    sfx = 'synctest'
    match preset:
        case 'full':
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF(filter_trigger=True, filter_vertex=True, filter_muons=True,
                                   filter_jpsi=True, filter_jets=True, filter_higgs=True)
            mcsig.snapshotRDF(user_sfx=sfx)
            
            mcbkg = an.JPsiCCLoader('MC_BKG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF(filter_trigger=True, filter_vertex=True, filter_muons=True,
                                   filter_jpsi=True, filter_jets=True, filter_higgs=True)
            mcbkg.snapshotRDF(user_sfx=sfx)

            data = an.JPsiCCLoader('DATA', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            data.createDataRDF('event_spec_data_skim.json')
            data.defineColumnsRDF(filter_trigger=True, filter_vertex=True, filter_muons=True,
                                    filter_jpsi=True, filter_jets=True, filter_higgs=True)
            data.snapshotRDF(user_sfx=sfx)

            analyzer = an.JPsiCCAnalyzer(2018, '202410', 'GF', 'CMSSW_13_3_0', weights=True)
            analyzer.readSnapshotSAMP('MC_SIG', f'snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root')
            analyzer.readSnapshotSAMP('MC_BKG1', f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root',
                                        sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            analyzer.readSnapshotSAMP('MC_BKG2', f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root',
                                        sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            analyzer.readSnapshotSAMP('DATA', f'snapshot_DATA_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root')
            # analyzer.addCutOptions('')
            analyzer.stackMultiHistos(draw_indiv=False,
                                      user_sfx=sfx)
            
        case 'stack':
            analyzer = an.JPsiCCAnalyzer(2018, '202410', 'GF', 'CMSSW_13_3_0', weights=True)
            analyzer.readSnapshotSAMP('MC_SIG', f'snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root')
            analyzer.readSnapshotSAMP('MC_BKG1', f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root',
                                        sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            analyzer.readSnapshotSAMP('MC_BKG2', f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root',
                                        sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            analyzer.readSnapshotSAMP('DATA', f'snapshot_DATA_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root')
            analyzer.addCutOptions('')
            analyzer.stackMultiHistos(keys=['jetPt',
                                            'goodJetPt',
                                            'jetClose_pt',
                                            'jetFar_pt'],
                                      user_sfx=sfx)

        case 'mcsig':
            # verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF(filter_trigger=True, filter_vertex=True, filter_muons=True,
                                   filter_jpsi=True, filter_jets=True, filter_higgs=True)
            mcsig.snapshotRDF(user_sfx=sfx)
            mcsig.readSnapshot(f'snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root')
            mcsig.makeHistos(genplots=False, cdf=False, save=True, draw=True, normSIG=False,
                             user_sfx=sfx)
            
        case 'mcbkg':
            mcbkg = an.JPsiCCLoader('MC_BKG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                                   filter_jpsi=False, filter_jets=True, filter_higgs=False)
            mcbkg.snapshotRDF(user_sfx=sfx)
            mcbkg.readSnapshot(f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root')
            mcbkg.makeHistos(genplots=False, cdf=False, save=True, draw=True, normSIG=False,
                             user_sfx=sfx)
            
        case 'data':
            databkg = an.JPsiCCLoader('DATA', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            databkg.createDataRDF('event_spec_data_skim.json')
            databkg.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                                   filter_jpsi=False, filter_jets=True, filter_higgs=False)
            databkg.snapshotRDF(user_sfx=sfx)
            databkg.readSnapshot(f'snapshot_DATA_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root')
            databkg.makeHistos(genplots=False, cdf=False, save=True, draw=True, normSIG=False,
                               user_sfx=sfx)

        case 'allsnap':
            mcbkg = an.JPsiCCLoader('MC_BKG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF()
            mcbkg.snapshotRDF()
            mcbkg.makeCutFlowPlot()
            
            databkg = an.JPsiCCLoader('DATA', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            databkg.createDataRDF('event_spec_DATA_skim.json')
            databkg.defineColumnsRDF()
            databkg.snapshotRDF()
            databkg.makeCutFlowPlot()
            
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF()
            mcsig.snapshotRDF()
            mcsig.makeCutFlowPlot()

        case 'gen':
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcsig.readSnapshot(f'snapshot_MC_SIG_2018_GF_v202410_{today_format}_WEIGHT_{sfx}.root')
            mcsig.cut('nMuon<2')
            mcsig.cut('gen_muplus_pt>3.5')
            mcsig.cut('gen_muminus_pt>3.5')
            # mcsig.cut('gen_muplus_eta>-2.4 && gen_muplus_eta<2.4')
            # mcsig.cut('gen_muminus_eta>-2.4 && gen_muminus_eta<2.4')
            mcsig.makeHistos(genplots=True, plot=True, draw=True, normSIG=False, cdf=False, save=True, user_sfx='nmuon<2+gen_muplus_pt>3.5+gen_muminus_pt>3.5')

        case 'allplots' | 'noweight' | 'cut':
            if preset=='noweight': use_weight, draw_indiv = False, False
            else: use_weight, draw_indiv = True, True
            
            mcsig = an.JPsiCCAnalyzer('MC_SIG', 2018, '202410', 'GF', 'CMSSW_13_3_0', weights=use_weight)
            mcbkg1 = an.JPsiCCAnalyzer('MC_BKG1', 2018, '202410', 'GF', 'CMSSW_13_3_0', weights=use_weight)
            mcbkg2 = an.JPsiCCAnalyzer('MC_BKG2', 2018, '202410', 'GF', 'CMSSW_13_3_0', weights=use_weight)
            mcbkg1.readSnapshotSAMP('MC_BKG1', f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_{sfx}.root',
                                    sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            mcbkg2.readSnapshotSAMP('MC_BKG2', f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_20241119_{today_format}_{sfx}.root',
                                    sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            # an.readSnapshotSAMP('DATA', 'snapshot_DATA_2018_GF_v202406_20240802_WEIGHT.root')
            mcsig.readSnapshotSAMP('MC_SIG', f'snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_20241119_{today_format}_{sfx}.root')

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