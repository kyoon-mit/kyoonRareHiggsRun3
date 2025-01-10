'''
analyzer
+++++++++++++++++++++++++++++++++++++++++++

JPsiCC analysis-specific analyzer
'''

from kytools import jsonreader, rootpdf
from jpsicc import rdfdefines
from datetime import date
import os, pickle
import matplotlib.pyplot as plt
import seaborn as sns
import ROOT
from time import process_time

ROOT.EnableImplicitMT()

class JPsiCCLoader:
    '''JPsiCC loader.

    Useful for loading one sample at a time, or for stacking histograms.

    Args:
        SAMP (str): Either one of the following options.
            'DATA', 'MC_BKG', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        weights (bool, optional): Whether to use weights.
            Defaults to True
    
    Raises:
        ValueError: If the string provided for SAMP is not among the options.
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
    '''
    def __init__(self, SAMP, YEAR, VERS, CAT, CMSSW, weights=True):
        match SAMP:
            case 'DATA': self._DATA, self._MODE = True, 'BKG'
            case 'MC_BKG': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG1': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG2': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG3': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG4': self._DATA, self._MODE = False, 'BKG'
            case 'MC_SIG': self._DATA, self._MODE = False, 'SIG'
            case _: raise ValueError(f'SAMP={SAMP} is not a valid option.')
        if not type(YEAR) is int: raise TypeError(f'YEAR must be an integer.')
        if not type(VERS) is str: raise TypeError(f'VERS must be a string.')
        if not type(CAT) is str: raise TypeError(f'CAT must be a string.')
        self.SAMPLE, self.YEAR, self.VERSION, self.CAT, self.CMSSW = SAMP, YEAR, VERS, CAT, CMSSW
        self.SAMPLENAME = f'{self.SAMPLE}' # Book sample name for plots
        self._weights = weights # bool

        self._branches = []   # TBranches to take RDF snapshot of
        self._rdf = None      # Placeholder for RDF
        self._rdfbranches = [] # TBranches in the RDF
        self._hists = dict()  # Dictionary of histogram objects
        self._models = dict() # Dictionary of histogram models
        self._cutflow = dict() # Dictionary of cut flow ({str: float})

        today = date.today()
        self._date = f'{today.year}{today.month:02}{today.day:02}'
        self._anpath = 'JPsiCC'
        self._plotsavedir = os.path.join(os.environ['HRARE_DIR'], self._anpath, 'plots', f'v{self.VERSION}', self._date, CAT)
        if not os.path.exists(self._plotsavedir): os.makedirs(self._plotsavedir)
        self._sfx = f'{SAMP}_{self.YEAR}_{self.CAT}_v{self.VERSION}_{self.CMSSW}_{self._date}_{"WEIGHT" if weights else "NOWEIGHT"}'

        # Color scheme: http://arxiv.org/pdf/2107.02270
        self._orange = ROOT.kOrange + 1
        self._blue = ROOT.kAzure - 4
        self._trueblue = ROOT.kBlue
        self._red = ROOT.kRed - 4
        self._purple = ROOT.kMagenta - 5
        self._gray = ROOT.kGray + 1
        self._violet = ROOT.kViolet + 2

        self._draw_option = 'P SAME' if self._DATA else 'HIST'
    
    def __draw_hist(self, histo1d, model1d, draw_option='HIST'):
        '''Internal method for drawing histogram.

        Args:
            histo1d (ROOT.TH1D): Histogram object.
            model1d (tuple(str)): Tuple of the histogram description.
                (name, title, nbins, xmin, xmax)
            draw_option (str, optional): Draw option to pass to TH1D.Draw() method.
                Defaults to 'HIST'.

        Returns:
            c (ROOT.TCanvas): ROOT.TCanvas object containing the histogram.
        '''
        c = ROOT.TCanvas()
        bin_width = (model1d[4]-model1d[3])/model1d[2]
        histo1d.Draw(draw_option)
        histo1d.SetTitle(model1d[1])
        histo1d.GetXaxis().SetTitle(model1d[1])
        histo1d.GetYaxis().SetTitle(f'Events / {bin_width:.3g}')
        c.Update()
        return c

    def __set_hist_style(self, histo1d, SAMP):
        '''Internal method for setting the histogram style.

        Args:
            histo1d (ROOT.TH1D): The histogram.
            SAMP (str):  Either one of the following three options.
                'DATA', 'MC_BKG1', 'MC_BKG2', 'MC_BKG3', MC_BKG4', 'MC_SIG'.

        Returns:
            histo1d (ROOT.TH1D): The histogram.
        '''
        match SAMP:
            case 'DATA':
                histo1d.SetMarkerStyle(ROOT.kFullSquare)
                histo1d.SetMarkerSize(0.5)
            case 'MC_BKG':
                histo1d.SetFillColorAlpha(ROOT.kWhite, 1.0)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
            case 'MC_BKG1':
                histo1d.SetFillColorAlpha(self._red, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
            case 'MC_BKG2':
                histo1d.SetFillColorAlpha(self._orange, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
            case 'MC_BKG3':
                histo1d.SetFillColorAlpha(self._gray, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
            case 'MC_BKG4':
                histo1d.SetFillColorAlpha(self._violet, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
            case 'MC_SIG':
                histo1d.SetFillColorAlpha(ROOT.kWhite, 0.0)
                histo1d.SetLineColorAlpha(self._blue, 1.0)
                histo1d.SetLineWidth(2)
            case _:
                raise ValueError(f'{SAMP} is not a valid option.')
        return histo1d

    def plot_hists(self, hist_dict, draw=True, normSIG=True, cdf=False, save=False, user_sfx=''):
        '''Plot histograms from dictionary of definitions.

        Args:
            hist_dict (dict): Dictionary of histogram definitions.
            draw (bool, optional): Whether to draw histogram.
                Defaults to True.
            normSIG (bool, optional): Whether to normalize the signal histogram to 1.
                This normalization applies to the events that are within the x axis range.
                Defaults to True.
            cdf (bool, optional): Whether to draw a CDF of each histogram.
                If True, additional histograms will be created with keys ending
                with the suffix, '_cdf'.
                Defaults to False.
            save (bool, optional): Whether to save the histograms to a ROOT file.
                Defaults to False.
            user_sfx (str, optional): Suffix to add to the names of the files.
                Defaults to ''.

        Returns:
            (None)
        '''
        for key, hdef in hist_dict.items():
            # Create Histo1D
            try:
                hname = f'{hdef["name"]}_{self._sfx}{"_" if user_sfx else ""}{user_sfx}'
                model1d = (hname, hdef['title'], hdef['bin'], hdef['xmin'], hdef['xmax'])
            except KeyError:
                print(f'Not drawing histogram for {key} because of KeyError.')
                continue
            if draw: ROOT.gStyle.SetOptStat('eimruo')
            if self._weights: histo1d = self._rdf.Histo1D(model1d, hdef['name'], 'w')
            else: histo1d = self._rdf.Histo1D(model1d, hdef['name'])

            # Set color
            histo1d = self.__set_hist_style(histo1d, self.SAMPLE)
            if self._MODE=='SIG' and normSIG:
                histo1d.Scale(1/histo1d.Integral())

            # Save histogram
            histo1d.SetDirectory(0)
            self._hists[key] = histo1d
            self._models[key] = model1d
            print(f'INFO: Made histogram {key} for {self.SAMPLE}.')

            # Add CDF
            if cdf:
                normed_histo1d = histo1d.Clone()
                normed_histo1d.SetDirectory(0)
                if normSIG: normed_histo1d.Scale(1/normed_histo1d.Integral())
                histo1d_cdf = normed_histo1d.GetCumulative()
                model1d_cdf = (f'{hname}_cdf', f'{hdef["title"]}_cdf', hdef['bin'], hdef['xmin'], hdef['xmax'])
                self._hists[f'{key}_cdf'] = histo1d_cdf
                self._models[f'{key}_cdf'] = model1d_cdf

            # Draw histogram
            if draw:
                c = self.__draw_hist(histo1d=histo1d, model1d=model1d, draw_option=self._draw_option)
                c.SaveAs(os.path.join(self._plotsavedir, f'{hname}.png'))
                c.Close()
                ROOT.gStyle.SetOptStat(0)
                if cdf:
                    c = self.__draw_hist(histo1d=histo1d_cdf, model1d=model1d_cdf, draw_option=self._draw_option)
                    c.SaveAs(os.path.join(self._plotsavedir, f'{hname}_cdf.png'))
        return
    
    def saveHistos(self, user_sfx=''):
        '''Save histograms to a ROOT file.

        Args:
            user_sfx (str, optional): Suffix to add to the names of the ROOT file.
                Defaults to ''.

        Returns:
            (None)
        '''
        fname = f'histos_{self._sfx}{"_" if user_sfx else ""}{user_sfx}.root'
        rootfile = ROOT.TFile.Open(fname, 'RECREATE')
        rootfile.cd()
        print('{}INFO: a ROOT file has been created. >> {} {}'.format('\033[1;33m', fname, '\033[0m'))
        for histo in self._hists.values():
            histo.Write()
        rootfile.Write()
        rootfile.Close()
        return

    def createDataRDF(self, event_spec_json):
        '''Open a sample using JSON spec and load it onto a RDF.

        Args:
            event_spec_json (str): Name of the JSON file for the "Events" tree.
            data (bool, optional): Whether to use data.

        Returns:
            (None)

        Raises:
            Exception: If attempted for not DATA.
        '''
        if not self._DATA:
            raise Exception('Please try again for data.')
        self._rdf = jsonreader.get_rdf_from_json_spec(self._anpath, event_spec_json)
        self._rdf, self._branches = rdfdefines.rdf_def_sample_meta(self._rdf, self._branches)
        self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_weights(None, self._rdf, self._branches, self._cutflow, data=self._DATA)
        return
    
    def createWeightedRDF(self, run_spec_json, event_spec_json):

        '''Open a sample using JSON spec and load it onto a RDF with proper weights.

        Args:
            run_spec_json (str): Name of the JSON file for the "Runs" tree.
            event_spec_json (str): Name of the JSON file for the "Events" tree.
            data (bool, optional): Whether to use data.

        Returns:
            (None)
        '''
        rdf_runs = jsonreader.get_rdf_from_json_spec(self._anpath, run_spec_json)
        self._rdf = jsonreader.get_rdf_from_json_spec(self._anpath, event_spec_json)
        rdf_runs, _ = rdfdefines.rdf_def_sample_meta(rdf_runs)
        self._rdf, self._branches = rdfdefines.rdf_def_sample_meta(self._rdf, self._branches)
        self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_weights(rdf_runs, self._rdf, self._branches, self._cutflow, data=self._DATA)
        return

    def defineColumnsRDF(self,
                         filter_trigger=True,
                         filter_vertex=True,
                         filter_muons=True,
                         filter_jpsi=True,
                         filter_jets=True,
                         filter_higgs=True):
        '''Define columns for the internal RDF.

        The columns will be different depending on the category.

        Args:
            rdf (ROOT.RDataFrame): Input ROOT.RDataFrame.
            filter_trigger (optional): Defaults to True.
            filter_vertex (optional): Defaults to True.
            filter_muons (optional): Defaults to True.
            filter_jpsi (optional): Defaults to True.
            filter_jets (optional): Defaults to True.
            filter_higgs (optional): Defaults to True.

        Returns:
            (None)
        '''
        start_time = process_time()
        match self.CAT:
            case 'GF':
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_triggers(self._rdf, self.CAT, self.YEAR, self.CMSSW, self._branches, self._cutflow, filter=filter_trigger)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_vertex(self._rdf, self.CMSSW, self._branches, self._cutflow, filter=filter_vertex)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_muons(self._rdf, self.CMSSW, self._branches, self._cutflow, filter=filter_muons)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_jpsi(self._rdf, self.CMSSW, self._branches, self._cutflow, filter=filter_jpsi)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_jets(self._rdf, self.CMSSW, self._branches, self._cutflow, data=self._DATA, filter=filter_jets)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_muon_jet_matching(self._rdf, self.CMSSW, self._branches, self._cutflow, data=self._DATA, filter=False)
                if self.SAMPLE=='MC_SIG':
                    self._rdf, self._branches = rdfdefines.rdf_def_genpart(self._rdf, self._branches)
                if not self._DATA:
                    self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_higgs(self._rdf, self.CMSSW, self._branches, self._cutflow, filter=filter_higgs)
            case 'NANOAOD_JETS':
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_jets(self._rdf, self.CMSSW, self._branches, self._cutflow, data=self._DATA, filter=filter_jets)
            case _: pass
        end_time = process_time()
        print('{}INFO: defineColumnsRDF process time is {:.8f}s.{}'.format('\033[1;32m', end_time - start_time, '\033[0m'))
        return
    
    def cut(self, cut, define_col='', define_exp='', weights=True):
        '''Filter the internal RDF using the cut expression.

        The internal self._cutflow object is updated when the cut is applied.

        Args:
            cut (str): String expression for the cut filter.
            define_col (str, optional): Create a new RDF column by this name.
                If '', no definition is created.
                Defaults to ''.
            define_exp (str, optional): The RDF column definition.
                Defaults to ''.
            weights (bool, optional): Whether to use weights in the computation
                of the cut flow and the number of events.
                In order to use this, there must be a 'w' column in the RDF.
                Defaults to True.

        Raises:
            Exception: If the cut cannot be applied.
            KeyError: If the weights option is specified but a column named 'w' does not exist.

        Returns:
            nevents (float): Number of events after the cut.
        '''
        if define_col:
            try: self._rdf = self._rdf.Define(str(define_col), str(define_exp))
            except Exception as err:
                print(err)
        try: self._rdf = self._rdf.Filter(cut)
        except Exception: raise Exception(f'The cut {cut} cannot be applied.')
        if weights:
            if 'w' not in self._rdf.GetColumnNames():
                raise KeyError('Weight column \'w\' does not exist in the RDF.')
            nevents = self._rdf.Sum('w').GetValue()
        else: nevents = self._rdf.Count().GetValue()
        self._cutflow = rdfdefines.add_cut_flow(self._cutflow, cut, nevents)
        return nevents
    
    def countNevents(self, cut='', weights=True):
        '''Count the number of events in the range.

        The cut which is applied does not affect the RDF stored in this class.

        Args:
            cut (str, optional): String expression for the cut filter, which applies the range.
                If '', no cut is applied.
                Defaults to ''.
            weights (bool, optional): Whether to use weights in this computation.
                In order to use this, there must be a 'w' column in the RDF.
                Defaults to True.

        Raises:
            Exception: If the cut cannot be applied.
            KeyError: If the weights option is specified but a column named 'w' does not exist.

        Returns:
            nevents (float): Number of events in the range.
        '''
        if cut:
            try: tmp_rdf = self._rdf.Filter(cut)
            except Exception: raise Exception(f'The cut {cut} cannot be applied.')
        else: tmp_rdf = self._rdf
        if weights:
            if 'w' not in tmp_rdf.GetColumnNames():
                raise KeyError('Weight column \'w\' does not exist in the RDF.')
            nevents = tmp_rdf.Sum('w').GetValue()
        else: nevents = tmp_rdf.Count().GetValue()
        return nevents
    
    def selectSampleRDF(self, sample_name):
        '''Select sample name in RDF.

        Args:
            rdf (ROOT.RDataFrame): Input ROOT.RDataFrame.
            sample_name (str): Name of the sample.

        Raises:
            ValueError: If sample_name provided is not str.
        
        Returns:
            (None)
        '''
        if type(sample_name) is not str: raise ValueError('sample_name must be str.')
        self._rdf = (self._rdf.Filter(f'sample=="{sample_name}"'))
        self.SAMPLENAME += f': {sample_name.split("_")[0]}'
        sum_weights = self._rdf.Sum('w').GetValue()
        print(f'Sum of weights for {sample_name}: {sum_weights:.6f}')
        return

    def retrieveRDF(self):
        '''Retrieve the internal RDF.

        Returns:
            rdf (ROOT.RDataFrame): The retrieved RDF.
        '''
        return self._rdf
    
    def retrieveCutFlow(self):
        '''Retrieve the internal cut flow dictionary.

        Returns:
            cutflow (dict): The retrieved cut flow dictionary.
        '''
        return self._cutflow

    def snapshotRDF(self, name='', make_pkl=True, user_sfx=''):
        '''Create a snapshot of the RDF.

        Outputs a file whose name is printed in the terminal.

        Args:
            name (str, optional): Name of the snapshot to create.
                Defaults to ''. If '', then the snapshot name will be determined automatically.
            make_pkl (bool, optional): Whether to save a pickle file.
                If True, a file with the extension .pkl and same name/suffix as the snapshot .root file
                will be created. This will contain the self._cutflow dictionary.
                Defaults to True.
            user_sfx (str, optional): Suffix to add to the names of the files.
                Only valid when name is ''.
                Defaults to ''.

        Returns:
            (None)
        '''
        start_time = process_time()
        if name=='': fname = f'snapshot_{self._sfx}{"_" if user_sfx else ""}{user_sfx}.root'
        else: fname = name
        self._rdf.Snapshot('Events', fname, self._branches)
        if not self._rdfbranches: self._rdfbranches = self._rdf.GetColumnNames()
        self._snapshotname = fname
        print('{}INFO: a snapshot has been created. >> {} {}'.format('\033[1;33m', fname, '\033[0m'))
        if make_pkl:
            pklname = f'{fname.removesuffix(".root")}.pkl'
            with open(pklname, 'wb') as f:
                pickle.dump(self._cutflow, f)
            self._picklename = pklname
            print('{}INFO: a pickle has been created. >> {} {}'.format('\033[1;33m', pklname, '\033[0m'))
        print(self._cutflow)
        end_time = process_time()
        print('{}INFO: snapshotRDF process time is {:.8f}s.{}'.format('\033[1;32m', end_time - start_time, '\033[0m'))
        return
    
    def readSnapshot(self, filename, treename='Events', read_pkl=True):
        '''Retrieve RDF from a saved snapshot ROOT file.

        Args:
            filename (str): Name of the snapshot ROOT file.
            treename (str, optional): Name of the TTree in the ROOT file.
                Defaults to 'Events'.
            read_pkl (bool, optional): Whether to read the corresponding pickle as well.

        Returns:
            (None)
        '''
        self._rdf = ROOT.RDataFrame(treename, filename)
        if not self._rdfbranches: self._rdfbranches = self._rdf.GetColumnNames()
        if read_pkl:
            pklname = f'{filename.removesuffix(".root")}.pkl'
            with open(pklname, 'rb') as f:
                self._cutflow = pickle.load(f)
        return

    def getVarInfo(self, varname):
        '''Get the information associated with a variable in the RDF.

        Args:
            varname (str): Name of the variable.

        Returns:
            stats (ROOT.TStatistic): TStatistic of the specified variable.
            vartype (str): The type of the variable.

        Raises:
            KeyError: If the variable is not in the RDF.
        '''
        if varname not in self._rdf.GetColumnNames():
            raise KeyError(f'The variable {varname} is not in the RDF.')
        stats, vartype = self._rdf.Stats(varname), self._rdf.GetColumnType(varname)
        return stats, vartype

    def makeHistos(self, plot=True, draw=True, genplots=False, keys=[], normSIG=True, cdf=False, save=False, user_sfx=''):
        '''Make histograms from the internal RDF.

        If 'plot' is True, outputs histograms in the 'plots' directory. Otherwise,
        returns a dictionary containing histogram definitions.

        Args:
            plot (bool, optional): If not true, function returns hist_dict.
                Defaults to True.
            draw (bool, optional): Whether to draw histograms.
                Defaults to True.
            genplots (bool, optional): Whether to add gen-level plots.
                Defaults to False.
            keys (list(str), optional): Specify the keys to draw. If the list is
                empty, it will draw every histogram.
                Defaults to [].
            normSIG (bool, optional): Whether to normalize the signal histogram to 1.
                This normalization applies to the events that are within the x axis range.
                Defaults to True.
            cdf (bool, optional): Whether to draw a CDF of each histogram.
                If True, additional histograms will be created with keys ending
                with the suffix, '_cdf'.
                Defaults to False.
            save (bool, optional): Whether to save the histograms to a ROOT file.
                Defaults to False.
            user_sfx (str, optional): Suffix to add to the names of the files.
                Defaults to ''.

        Returns:
            (None) or hist_dict (dict)
        '''
        hist_defs = jsonreader.get_object_from_json(anpath='JPsiCC',
                                                    jsonname=rdfdefines.select_json_file(CMSSW=self.CMSSW),
                                                    keys=['NANOAOD_to_RDF'])
        hist_dict = dict()
        match self.CAT:
            case 'GF':
                hist_dict.update(hist_defs['Jpsi'])
                hist_dict.update(hist_defs['muon'])
                hist_dict.update(hist_defs['vertex'])
                hist_dict.update(hist_defs['jet'])
                hist_dict.update(hist_defs['muon_jet_matching'])
                hist_dict.update(hist_defs['higgs'])
                if not self._DATA: hist_dict.update(hist_defs['jet_mconly'])
                if genplots: hist_dict.update(hist_defs['gen'])
            case 'NANOAOD_JETS':
                hist_dict.update(hist_defs['jet'])
            case _: pass
        if keys: hist_dict = {key: val for key, val in hist_dict.items() if key in keys}
        if plot: self.plot_hists(hist_dict, draw=draw, normSIG=normSIG, cdf=cdf, user_sfx=user_sfx)
        if save: self.saveHistos(user_sfx=user_sfx)
        else: return hist_dict
    
    def retrieveHisto(self, key):
        '''Retrieve a histogram object using a key.

        Args:
            key (str): The key to the histogram.

        Returns:
            histo1d (ROOT.TH1D)/(None): Histogram object if key exists. Otherwise,
                returns NoneType.
        '''
        if key in self._hists:
            return self._hists[key]
        return
    
    def makeCutFlowPlot(self):
        '''Make a plot showing the cut flow.

        Returns:
            (None)

        Raises:
            RuntimeError: If the internal self._cutflow object is empty.
        '''
        if len(self._cutflow)==0: raise RuntimeError('The self._cutflow object is empty.')
        # plt.figure(figsize=(12,6))
        ax = sns.barplot(y=list(self._cutflow.keys()), x=list(self._cutflow.values()), orient='y', errorbar=None)
        ax.bar_label(ax.containers[0], fontsize=9)
        ax.margins(x=0.2)
        plt.xlabel('sum of weights')
        plt.ylabel('cut')
        plt.title(f'cutflow_{self._sfx}')

        fname = os.path.join(self._plotsavedir, f'cutflow_{self._sfx}.png')
        plt.savefig(fname=fname, bbox_inches='tight', dpi=120)
        print('{}INFO: a cutflow figure has been created. >> {} {}'.format('\033[1;33m', fname, '\033[0m'))
        plt.close('all')
        return

    def stackHistos(self, analyzer):
        '''Stack and plot histograms.

        Args:
            analyzer (JPsiCCLoader): Another analyzer object.

        Returns:
            (None)
        '''
        thisSAMPLE, otherSAMPLE = self.SAMPLE, analyzer.SAMPLE
        for key, item in self._hists.items():
            other_hist = analyzer.retrieveHisto(key)
            if other_hist is not None:
                this_hist = self.__set_hist_style(item, thisSAMPLE)
                other_hist = self.__set_hist_style(other_hist, otherSAMPLE)
                # Create THStack
                hs = ROOT.THStack()
                hs.Add(this_hist.GetPtr(), self._draw_option)
                hs.Add(other_hist.GetPtr(), analyzer._draw_option)
                # Draw histograms
                model = self._models[key]
                c_hstack = self.__draw_hist(hs, model, draw_option='HIST NOSTACK') #TODO: fix model hname
                # Save TCanvas
                sfx = 'STACK_' + self._sfx.lstrip(f'{self.SAMPLE}_')
                c_hstack.SaveAs(os.path.join(self._plotsavedir, f'{key}_{sfx}.png'))
        return

################################################################################

       ######                 #####      #      #####    ##    #  #     #   
      #### ###                #    #    # #     #    #   # #   #   #   #
       #######                ######   #####    #####    #  #  #    ###
             #                #    #  #     #   #    #   #   # #     #
          ####                #####  #       #  #     #  #    ##     #
             ##
              ##
               ##
                ###
                 ####     
                  ########
                   #############
                  ############################
                 ########################################  
                 ##########################################
                 ###########################################
                  ###########################################
                   ##########################################
                    ########            #################### # 
                     ####                    ##############   # 
                     ###                         #########   #
                     ###                           #######  # 
                     ###                            ######   #
                     ####                            #####    #
                     ###                              ####   ##
                    ###                               ### 
                   ###                               ### 
                  ###                               ### 
                 ###                               ###
                 ###                               ### 

################################################################################

class JPsiCCAnalyzer:
    '''JPsiCC analyzer.

    Args:
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        weights (bool, optional): Whether to use weights.
            Defaults to True.
    
    Raises:
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
    '''
    def __init__(self, YEAR, VERS, CAT, CMSSW, weights=True):
        if not type(YEAR) is int: raise TypeError(f'YEAR must be an integer.')
        if not type(VERS) is str: raise TypeError(f'VERS must be a string.')
        if not type(CAT) is str: raise TypeError(f'CAT must be a string.')
        self.YEAR, self.VERSION, self.CAT, self.CMSSW = YEAR, VERS, CAT, CMSSW
        self._loaders = dict() # Placeholder for JPsiCCLoader
        self._weights = weights
        self._selections = []

        # Color scheme: http://arxiv.org/pdf/2107.02270
        self._orange = ROOT.kOrange + 1
        self._blue = ROOT.kAzure - 4
        self._trueblue = ROOT.kBlue
        self._red = ROOT.kRed - 4
        self._purple = ROOT.kMagenta - 5
        self._gray = ROOT.kGray + 1
        self._violet = ROOT.kViolet + 2

        today = date.today()
        self._date = f'{today.year}{today.month:02}{today.day:02}'
        self._anpath = 'JPsiCC'
        self._plotsavedir = os.path.join(os.environ['HRARE_DIR'], self._anpath, 'plots', f'v{self.VERSION}', self._date, CAT)
        self._sfx = f'{self.YEAR}_{self.CAT}_v{self.VERSION}_{self._date}'

    def __set_hist_style(self, histo1d, SAMP):
        '''Set the histogram style.

        Args:
            histo1d (ROOT.TH1D): The histogram.
            SAMP (str):  Either one of the following three options.
                'DATA', 'MC_BKG1', 'MC_BKG2', 'MC_BKG3', MC_BKG4', 'MC_SIG'.

        Returns:
            histo1d (ROOT.TH1D): The histogram.
        '''
        draw_option = ''
        match SAMP:
            case 'DATA':
                histo1d.SetMarkerStyle(ROOT.kFullSquare)
                histo1d.SetMarkerSize(0.5)
                draw_option = 'P SAME'
            case 'MC_BKG':
                histo1d.SetFillColorAlpha(ROOT.kWhite, 1.0)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
                draw_option = 'HIST'
            case 'MC_BKG1':
                histo1d.SetFillColorAlpha(self._red, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
                draw_option = 'HIST'
            case 'MC_BKG2':
                histo1d.SetFillColorAlpha(self._orange, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
                draw_option = 'HIST'
            case 'MC_BKG3':
                histo1d.SetFillColorAlpha(self._gray, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
                draw_option = 'HIST'
            case 'MC_BKG4':
                histo1d.SetFillColorAlpha(self._violet, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
                draw_option = 'HIST'
            case 'MC_SIG':
                histo1d.SetFillColorAlpha(ROOT.kWhite, 0.0)
                histo1d.SetLineColorAlpha(self._blue, 1.0)
                histo1d.SetLineWidth(2)
                draw_option = 'HIST SAME'
            case _:
                raise ValueError(f'{SAMP} is not a valid option.')
        return histo1d, draw_option

    def readSnapshotSAMP(self, SAMP, filename, treename='Events', sample_name='', read_pkl=True):
        '''Create JPsiCCLoader from a saved snapshot ROOT file.

        Args:
            SAMP (str): Either one of the following three options.
                'DATA', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
            filename (str): Name of the snapshot ROOT file.
            treename (str, optional): Name of the TTree in the ROOT file.
                Defaults to 'Events'.
            sample (str, optional): Name of the sample to filter on. If '' provided,
                it will use the entire events instead.
                Defaults to ''.
            read_pkl (bool, optional): Whether to read the corresponding pickle as well.
                Defaults to True.            

        Raises:
            ValueError: If the string provided for SAMP is not among the options.

        Returns:
            (None)
        '''
        match SAMP:
            case 'DATA': self._DATA, self._MODE = True, 'BKG'
            case 'MC_BKG': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG1': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG2': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG3': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG4': self._DATA, self._MODE = False, 'BKG'
            case 'MC_SIG': self._DATA, self._MODE = False, 'SIG'
            case _: raise ValueError(f'SAMP={SAMP} is not a valid option.')
        self._loaders[SAMP] = JPsiCCLoader(SAMP, self.YEAR, self.VERSION, self.CAT, self.CMSSW, self._weights)
        self._loaders[SAMP].readSnapshot(filename, treename, read_pkl=read_pkl)
        if sample_name: self._loaders[SAMP].selectSampleRDF(sample_name)

    def stackMultiHistos(self, draw_indiv=False, keys=[], user_sfx=''):
        '''Stack and plot mlutiple histograms.

        Args:
            draw_indiv (bool, optional): Whether to draw individual histograms.
                Defaults to False.
            keys (list(str), optional): Specify the keys to draw. If the list is
                empty, it will draw every histogram.
                Defaults to [].
            user_sfx (str, optional): Suffix to add to the names of the files.
                Defaults to ''.

        Raises:
            KeyError: If any one of MC_BKG, MC_SIG, or DATA is missing.

        Returns:
            (None)
        '''
        samples = list(self._loaders.keys())
        mcbkgs = [self._loaders[s] for s in samples if 'MC_BKG' in s]
        if not mcbkgs: raise KeyError('No MC_BKG sample is found. Please load first.')
        if 'MC_SIG' not in samples: raise KeyError('No MC_SIG sample is found. Please load first.')
        if 'DATA' not in samples: raise KeyError('No DATA sample is found. Please load first.')
        mcsig = self._loaders['MC_SIG']
        databkg = self._loaders['DATA']

        hist_dict = databkg.makeHistos(plot=False, draw=False, keys=keys)

        for samp in samples:
            self._loaders[samp].plot_hists(hist_dict, draw=draw_indiv, user_sfx=user_sfx)

        ROOT.gStyle.SetOptStat(0)
        for key, hdef in hist_dict.items():
            # Create canvas
            c_hstack = ROOT.TCanvas()

            # Set legend
            legend = ROOT.TLegend(0.64, 0.78, 0.95, 0.94)
            legend.SetBorderSize(1)
            legend.SetFillColorAlpha(ROOT.kWhite, 0.8)
            legend.SetTextSize(0.024)
            legend.SetMargin(0.2)

            # Create THStack
            hname = f'{hdef["name"]}_MULTI_STACK_{self.YEAR}_{self.CAT}'
            model1d = (hname, hdef['title'], hdef['bin'], hdef['xmin'], hdef['xmax'])
            hs = ROOT.THStack()

            # Stack mcbkg
            mcbkg_norm, mcbkg_max = 0., 0.
            for loader in mcbkgs:
                h = loader.retrieveHisto(key)
                h.SetStats(False)
                if h:
                    h, mcbkg_d_opt = self.__set_hist_style(h, loader.SAMPLE)
                    mcbkg_norm += h.Integral()
                    mcbkg_max += h.GetMaximum()
                    hs.Add(h.GetPtr(), mcbkg_d_opt)
                    legend.AddEntry(h.GetPtr(), f'{loader.SAMPLENAME} ({h.Integral():.2E})', 'F')
            hs.Draw()

            # Draw data
            data_hist = databkg.retrieveHisto(key)
            if data_hist:
                data_hist.SetStats(False)
                data_hist, databkg_d_opt = self.__set_hist_style(data_hist, databkg.SAMPLE)
                data_hist.Draw(databkg_d_opt)
                legend.AddEntry(data_hist.GetPtr(), f'{databkg.SAMPLENAME} ({data_hist.Integral():.2E})', 'P')

            # Draw signal
            sig_hist = mcsig.retrieveHisto(key)
            if sig_hist:
                sig_hist.SetStats(False)
                sig_hist, mcsig_d_opt = self.__set_hist_style(sig_hist, mcsig.SAMPLE)
                sig_hist.Scale(data_hist.Integral()/sig_hist.Integral()) # scale sig hist
                sig_hist.Draw(mcsig_d_opt)
                legend.AddEntry(sig_hist.GetPtr(), f'{mcsig.SAMPLENAME} ({sig_hist.Integral():.2E})', 'L')

            # Draw histograms
            hist_max_vals = [mcbkg_max, data_hist.GetMaximum(), sig_hist.GetMaximum()]
            hs.SetMaximum(round(max(hist_max_vals)*1.10))
            bin_width = (model1d[4]-model1d[3])/model1d[2]
            hs.SetTitle(key)
            hs.GetXaxis().SetTitle(model1d[1])
            hs.GetYaxis().SetTitle(f'Events / {bin_width:.3g}')
            legend.Draw('SAME')
            c_hstack.Update()

            # Save TCanvas
            sfx = f'MULTI_STACK_{self._sfx}{"_" if user_sfx else ""}{user_sfx}'
            c_hstack.SaveAs(os.path.join(self._plotsavedir, f'{key}_{sfx}.png'))

            # Print integrals
            print(f'MCBKG integral: {mcbkg_norm:.3f}, DATA integral: {data_hist.Integral():.3f}')
        return

    def countNeventsPerSample(self, cut, weights=True):
        '''Count the number of events in the range per sample.

        The cut which is applied does not affect the RDF stored in this class.

        Args:
            cut (str): String expression for the cut filter, which applies the range.
                If '', no cut is applied.
                Defaults to ''.
            weights (bool, optional): Whether to use weights in this computation.
                In order to use this, there must be a 'w' column in the RDF.
                Defaults to True.

        Returns:
            table_dict (float): Dictionary of the following format.
                {'selection': cut, <sample 1>: <nevents>, ...}
        '''
        table_dict = {'selection': cut}
        for sample, loader in self._loaders.items():
            nevents = loader.countNevents(cut, weights)
            table_dict[sample] = nevents
        return table_dict

    def applyCut(self, cut):
        '''Quickly apply a selection to all the loaded samples. 

        Args:
            cut (str): RDF expression for filtering based on selection.

        Returns:
            (None)
        '''
        for sample, loader in self._loaders.items():
            print(f'Applying cut {cut} to {sample}.')
            loader.cut(f'{cut}')
        return

    def addCutOptions(self, cuts):
        '''Does not apply the selections immediately, but saves them for later use.

        Args:
            cuts (str or list(str)): RDF expression or list of expressions.

        Raises:
            TypeError: If cuts is not str or list(str).
        
        Returns:
            (None)
        '''
        if isinstance(cuts, str): self._selections.append(cuts)
        elif isinstance(cuts, list) and all(isinstance(c, str) for c in cuts):
            self._selections += cuts
        else:
            raise TypeError('The value provided for cuts is not str or list(str)')

    def makeCutTable(self, modify=False, plot=True):
        '''Make a table of the number of events after applying each cut stored.

        Args:
            modify (bool, optional): Whether to modify the original RDFs.
                If False, each selection is applied independently from one another
                and have no effect the original RDFs.
                If True, the selections are applied sequentially.
            plot (bool, optional): Whether to plot the resulting table.
                Defaults to True.
        
        Returns:
            table_pd (pandas.DataFrame): Object containing the tabulated nevents after each cut.
        '''
        import pandas as pd # won't be using this outside the function
        from pandas.plotting import table
        table_dict = {'selection': 'None'}
        table_dict = self.countNeventsPerSample(cut='', weights=self._weights)
        table_df = pd.DataFrame(data=table_dict, index=[0])
        for cut in self._selections:
            if modify:
                self.applyCut(cut=cut)
            else:
                nevents = self.countNeventsPerSample(cut=cut, weights=self._weights)
                nevents['selection'] = cut
                table_df = pd.concat([table_df, pd.DataFrame(data=nevents, index=[0])], ignore_index=True)
        if modify:
            new_rows = {}
            for sample, loader in self._loaders.items(): # TODO: can't you be more efficient?
                cutflow = loader.retrieveCutFlow()
                if not 'selection' in new_rows:
                    new_rows['selection'] = list(cutflow)
                new_rows[sample] = list(cutflow.values())
            table_df = pd.concat([table_df, nevents], ignore_index=True)
        for sample in self._loaders:
            table_df[f'{sample}_ratio'] = table_df[sample]/table_df[sample].iloc[0]
        if plot:
            _, ax = plt.subplots()
            ax.axis('off')
            _ = table(ax, table_df, loc='center', cellLoc='center', colWidths=[.4]*len(table_df.index))
            plt_title = f'cut_table_{"sequential" if modify else "individual"}_{self._sfx}'
            plt_path = f'{os.path.join(self._plotsavedir, plt_title)}.png'
            plt.title(plt_title)
            plt.savefig(plt_path, bbox_inches='tight', dpi=120)
            print('{}INFO: a cut table figure has been created. >> {} {}'.format('\033[1;33m', plt_path, '\033[0m'))
            plt.close('all')
        return table_df

    def makeCutFlowPlotAll(self):
        '''Make a cut flow plot of all the loaders.
        
        Returns:
            (None)
        '''
        for loader in self._loaders.values():
            loader.makeCutFlowPlot()
        return

    def fitSignal(self, key_var, template):
        """Fit the signal RDF to a desired template.

        Args:
            key_var (str): Name of the variable to fit.
            template (str): Template to use for fitting.
        """
        fit = rootpdf.FittingTool(YEAR=self.YEAR, CAT=self.CAT, VARNAME=key_var)
        return