from jpsicc import analyzer as an
from datetime import date

if __name__=='__main__':
    print('Choose preset (mcbkg, databkg, mcsig, allsnap, gen, allplots, noweight, cut, ancut): ', end='')
    preset = input()
    today = date.today()
    today_format = f'{today.year}{today.month:02}{today.day:02}'
    match preset:
        case 'mcbkg':
            mcbkg = an.JPsiCCLoader('MC_BKG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                                   filter_jpsi=False, filter_jets=False, filter_higgs=False)
            mcbkg.snapshotRDF(user_sfx='no_filter')
            mcbkg.readSnapshot(f'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_no_filter.root')
            mcbkg.makeHistos(genplots=False, cdf=True, save=True, draw=True, normSIG=False,
                             user_sfx='no_filter')

        case 'databkg':
            databkg = an.JPsiCCLoader('DATA_BKG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            databkg.createDataRDF('event_spec_data_bkg_skim.json')
            databkg.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                                   filter_jpsi=False, filter_jets=False, filter_higgs=False)
            databkg.snapshotRDF(user_sfx='no_filter')
            databkg.readSnapshot(f'snapshot_DATA_BKG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_no_filter.root')
            databkg.makeHistos(genplots=False, cdf=True, save=True, draw=True, normSIG=False,
                             user_sfx='no_filter')

        case 'mcsig':
            # verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)
            mcsig = an.JPsiCCLoader('MC_SIG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                                   filter_jpsi=False, filter_jets=False, filter_higgs=False)
            mcsig.snapshotRDF(user_sfx='no_filter')
            mcsig.readSnapshot(f'snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_{today_format}_WEIGHT_no_filter.root')
            mcsig.makeHistos(genplots=False, cdf=True, save=True, draw=True, normSIG=False,
                             user_sfx='no_filter')
        
        case 'allsnap':
            mcbkg = an.JPsiCCLoader('MC_BKG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF()
            mcbkg.snapshotRDF()
            mcbkg.makeCutFlowPlot()
            
            databkg = an.JPsiCCLoader('DATA_BKG', 2018, '202410', 'GF', 'CMSSW_13_3_0')
            databkg.createDataRDF('event_spec_data_bkg_skim.json')
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
            mcsig.readSnapshot(f'snapshot_MC_SIG_2018_GF_v202410_{today_format}_WEIGHT_no_filter.root')
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
            mcbkg1.readSnapshotSAMP('MC_BKG1', 'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_20241119_WEIGHT_no_filter.root',
                                    sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            mcbkg2.readSnapshotSAMP('MC_BKG2', 'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_20241119_WEIGHT_no_filter.root',
                                    sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            # an.readSnapshotSAMP('DATA_BKG', 'snapshot_DATA_BKG_2018_GF_v202406_20240802_WEIGHT.root')
            mcsig.readSnapshotSAMP('MC_SIG', 'snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_20241119_WEIGHT_no_filter.root')

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

        case 'ancut':
            mcsig = an.JPsiCCAnalyzer('MC_SIG', 2018, '202410', 'GF', 'CMSSW_13_3_0', weights=use_weight)
            mcbkg1 = an.JPsiCCAnalyzer('MC_BKG1', 2018, '202410', 'GF', 'CMSSW_13_3_0', weights=use_weight)
            mcbkg2 = an.JPsiCCAnalyzer('MC_BKG2', 2018, '202410', 'GF', 'CMSSW_13_3_0', weights=use_weight)
            mcbkg1.readSnapshotSAMP('MC_BKG1', 'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_20241119_WEIGHT_no_filter.root',
                                    sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            mcbkg2.readSnapshotSAMP('MC_BKG2', 'snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_20241119_WEIGHT_no_filter.root',
                                    sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            # an.readSnapshotSAMP('DATA_BKG', 'snapshot_DATA_BKG_2018_GF_v202406_20240802_WEIGHT.root')
            mcsig.readSnapshotSAMP('MC_SIG', 'snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_20241119_WEIGHT_no_filter.root')
            mcsig.addCutOptions('')