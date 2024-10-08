from jpsicc import analyzer as an

def makePlotsNanoJet(SAMP, YEAR, VERS, CAT, CMSSW, run_spec_json, event_spec_json, user_sfx):
    '''Make plots for the NanoAOD jet samples.

    Args:
        SAMP (str): Either one of the following options.
            'DATA_BKG', 'MC_BKG', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        event_spec_json (str): Name of the JSON file for the "Events" tree.

    Returns:
        (None)
    '''
    loader = an.JPsiCCLoader(SAMP, YEAR, VERS, CAT, CMSSW=CMSSW)
    loader.createWeightedRDF(run_spec_json, event_spec_json)
    loader.defineColumnsRDF(filter_trigger=False, filter_vertex=False, filter_muons=False,
                            filter_jpsi=False, filter_jets=False, filter_higgs=False)
    loader.snapshotRDF(user_sfx=user_sfx)
    loader.readSnapshot(loader._snapshotname)
    loader.makeHistos(user_sfx=user_sfx)
    return

def run():
    '''Run script.

    Args:
        (None)
    
    Returns:
        (None)
    '''
    SAMP, YEAR, VERS, CAT, CMSSW = 'MC_SIG', 2018, 'v202407', 'NANOAOD_JETS', 'CMSSW_10_6_30'
    json_files = [
        ('run_spec_nanoaod_jets_default.json', 'event_spec_nanoaod_jets_default.json', 'default'),
        ('run_spec_nanoaod_jets_experiment.json', 'event_spec_nanoaod_jets_experiment.json', 'experiment'),
    ]
    for jf in json_files:
        makePlotsNanoJet(SAMP=SAMP, YEAR=YEAR, VERS=VERS, CAT=CAT, CMSSW=CMSSW,
                         run_spec_json=jf[0], event_spec_json=jf[1], user_sfx=jf[2])
    return

if __name__=='__main__':
    run()