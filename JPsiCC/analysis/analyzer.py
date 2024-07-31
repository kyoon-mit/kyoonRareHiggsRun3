'''
analyzer
+++++++++++++++++++++++++++++++++++++++++++

JPsiCC analysis-specific analyzer
'''

from kytools import jsonreader, rootpdf
import rdfdefines
from datetime import date
import os, pickle
import matplotlib.pyplot as plt
import seaborn as sns
import ROOT

ROOT.EnableImplicitMT()

class JPsiCCLoader:
    '''JPsiCC loader.

    Useful for loading one sample at a time, or for stacking histograms.

    Args:
        SAMP (str): Either one of the following options.
            'DATA_BKG', 'MC_BKG', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        weights (bool): Whether to use weights.
    
    Raises:
        ValueError: If the string provided for SAMP is not among the options.
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
    '''
    def __init__(self, SAMP, YEAR, VERS, CAT, weights=True):
        match SAMP:
            case 'DATA_BKG': self._DATA, self._MODE = True, 'BKG'
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
        self.SAMPLE, self.YEAR, self.VERSION, self.CAT = SAMP, YEAR, VERS, CAT
        self.SAMPLENAME = f'{self.SAMPLE}'
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
        self._sfx = f'{SAMP}_{self.YEAR}_{self.CAT}_v{self.VERSION}_{self._date}'
        if not weights: self._sfx += '_NOWEIGHT'
        if not os.path.exists(self._plotsavedir): os.makedirs(self._plotsavedir)

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
                'DATA_BKG', 'MC_BKG1', 'MC_BKG2', 'MC_BKG3', MC_BKG4', 'MC_SIG'.

        Returns:
            histo1d (ROOT.TH1D): The histogram.
        '''
        match SAMP:
            case 'DATA_BKG':
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

    def plot_hists(self, hist_dict, draw=True, scaleSIG=1., user_sfx=''):
        '''Plot histograms from dictionary of definitions.

        Args:
            hist_dict (dict): Dictionary of histogram definitions.
            draw (bool, optional): Whether to draw histogram.
                Defaults to True.
            scaleSIG (float, optional): Scale factor for the signal histogram.
                Defaults to 1.
            user_sfx (str, optional): Suffix to add to the names of the files.
                Defaults to ''.

        Returns:
            (None)
        '''
        for key, hdef in hist_dict.items():
            # Create Histo1D
            hname = f'{hdef["name"]}_{self._sfx}{"_" if user_sfx else ""}{user_sfx}'
            model1d = (hname, hdef['title'], hdef['bin'], hdef['xmin'], hdef['xmax'])
            if draw: ROOT.gStyle.SetOptStat('eimruo')
            if self._weights: histo1d = self._rdf.Histo1D(model1d, hdef['name'], 'w')
            else: histo1d = self._rdf.Histo1D(model1d, hdef['name'])

            # Set color
            histo1d = self.__set_hist_style(histo1d, self.SAMPLE)
            if self._MODE=='SIG': histo1d.Scale(scaleSIG)

            # Save histogram
            histo1d.SetDirectory(0)
            self._hists[key] = histo1d
            self._models[key] = model1d
            print(f'INFO: Made histogram {key} for {self.SAMPLE}.')

            # Draw histogram
            if draw:
                c = self.__draw_hist(histo1d=histo1d, model1d=model1d, draw_option=self._draw_option)
                c.SaveAs(os.path.join(self._plotsavedir, f'{hname}.png'))
                ROOT.gStyle.SetOptStat(0)
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

    def defineColumnsRDF(self):
        '''Define columns for the internal RDF.

        The columns will be different depending on the category.

        Args:
            rdf (ROOT.RDataFrame): Input ROOT.RDataFrame.

        Returns:
            (None)
        '''
        match self.CAT:
            case 'GF':
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_triggers(self._rdf, self.CAT, self.YEAR, self._branches, self._cutflow)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_jpsi(self._rdf, self._branches, self._cutflow)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_muons(self._rdf, self._branches, self._cutflow)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_vertex(self._rdf, self._branches, self._cutflow)
                self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_jets(self._rdf, self.CAT, self.YEAR, self._branches, self._cutflow, data=self._DATA)
                if self.SAMPLE=='MC_SIG':
                    self._rdf, self._branches = rdfdefines.rdf_def_genpart(self._rdf, self._branches)
                if not self._DATA:
                    self._rdf, self._branches, self._cutflow = rdfdefines.rdf_def_higgs(self._rdf, self._branches, self._cutflow)
            case _: pass
        return
    
    def cut(self, cut):
        '''Filter the internal RDF using the cut expression.

        Args:
            cut (str): String expression for the cut filter.

        Raises:
            Exception: If the cut cannot be applied.

        Returns:
            (None)
        '''
        try: self._rdf = self._rdf.Filter(cut)
        except Exception: raise Exception(f'The cut {cut} cannot be applied.')
        self._cutflow = rdfdefines.add_cut_flow(self._cutflow, cut, self._rdf.Sum('w').GetValue())
    
    def countNevents(self, cut='', weights=True):
        '''Count the number of events in the range.

        The cut which is applied does not affect the RDF stored in this class.

        Args:
            cut (str): String expression for the cut filter, which applies the range.
                If '', no cut is applied.
                Defaults to ''.
            weights (bool, optional): Whether to use weights in this computation.
                In order to use this, there must be a 'w' column in the RDF.
                Defaults to True.

        Raises:
            Exception: If the cut cannot be applied.
            KeyError: If the weight column named 'w' does not exist.

        Returns:
            nevents (float): Number of events in the range.
        '''
        if 'w' not in self._rdf.GetColumnNames(): raise KeyError('Weight column \'w\' does not exist in the RDF.')
        try: tmp_rdf = self._rdf.Filter(cut)
        except Exception: raise Exception(f'The cut {cut} cannot be applied.')
        if weights: nevents = tmp_rdf.Sum('w').GetValue()
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

    def snapshotRDF(self, name='', make_pkl=True):
        '''Create a snapshot of the RDF.

        Outputs a file whose name is printed in the terminal.

        Args:
            name (str, optional): Name of the snapshot to create.
                Defaults to ''. If '', then the snapshot name will be determined automatically.
            make_pkl (bool, optional): Whether to save a pickle file.
                If True, a file with the extension .pkl and same name/suffix as the snapshot .root file
                will be created. This will contain the self._cutflow dictionary.
                Defaults to True.

        Returns:
            (None)
        '''
        if name=='': fname = f'snapshot_{self._sfx}.root'
        else: fname = name
        self._rdf.Snapshot('Events', fname, self._branches)
        if not self._rdfbranches: self._rdfbranches = self._rdf.GetColumnNames()
        print('{}INFO: a snapshot has been created. >> {} {}'.format('\033[1;33m', fname, '\033[0m'))
        if make_pkl:
            pklname = f'{fname.rstrip(".root")}.pkl'
            with open(pklname, 'wb') as f:
                pickle.dump(self._cutflow, f)
            print('{}INFO: a pickle has been created. >> {} {}'.format('\033[1;33m', pklname, '\033[0m'))
        print(self._cutflow)
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
            pklname = f'{filename.rstrip(".root")}.pkl'
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

    def makeHistos(self, plot=True, draw=True, genplots=False, keys=[], scaleSIG=1., user_sfx=''):
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
            scaleSIG (float, optional): Scale factor for the signal histogram.
                Defaults to 1.
            user_sfx (str, optional): Suffix to add to the names of the files.
                Defaults to ''.

        Returns:
            (None) or hist_dict (dict)
        '''
        hist_defs = jsonreader.get_object_from_json(anpath='JPsiCC',
                                                    jsonname='rdf_hists.json',
                                                    keys=['NANOAOD_to_RDF'])
        hist_dict = dict()
        match self.CAT:
            case 'GF':
                hist_dict.update(hist_defs['Jpsi'])
                hist_dict.update(hist_defs['muon'])
                hist_dict.update(hist_defs['vertex'])
                if not self._DATA: hist_dict.update(hist_defs['higgs'])
                if genplots: hist_dict.update(hist_defs['gen'])
                # self.plot_hists(hist_defs['jet'])
            case _: pass
        if keys: hist_dict = {key: val for key, val in hist_dict.items() if key in keys}
        if plot: self.plot_hists(hist_dict, draw=draw, scaleSIG=scaleSIG, user_sfx=user_sfx)
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

        fname=os.path.join(self._plotsavedir, f'cutflow_{self._sfx}.png')
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
        weights (bool): Whether to use weights.
    
    Raises:
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
    '''
    def __init__(self, YEAR, VERS, CAT, weights=True):
        if not type(YEAR) is int: raise TypeError(f'YEAR must be an integer.')
        if not type(VERS) is str: raise TypeError(f'VERS must be a string.')
        if not type(CAT) is str: raise TypeError(f'CAT must be a string.')
        self.YEAR, self.VERSION, self.CAT = YEAR, VERS, CAT
        self._loaders = dict() # Placeholder for JPsiCCLoader
        self._weights = weights

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
                'DATA_BKG', 'MC_BKG1', 'MC_BKG2', 'MC_BKG3', MC_BKG4', 'MC_SIG'.

        Returns:
            histo1d (ROOT.TH1D): The histogram.
        '''
        draw_option = ''
        match SAMP:
            case 'DATA_BKG':
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

    def readSnapshotSAMP(self, SAMP, filename, treename='Events', sample_name=''):
        '''Create JPsiCCLoader from a saved snapshot ROOT file.

        Args:
            SAMP (str): Either one of the following three options.
                'DATA_BKG', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
            filename (str): Name of the snapshot ROOT file.
            treename (str, optional): Name of the TTree in the ROOT file.
                Defaults to 'Events'.
            sample (str, optional): Name of the sample to filter on. If '' provided,
                it will use the entire events instead.
                Defaults to ''.

        Raises:
            ValueError: If the string provided for SAMP is not among the options.

        Returns:
            (None)
        '''
        match SAMP:
            case 'DATA_BKG': self._DATA, self._MODE = True, 'BKG'
            case 'MC_BKG': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG1': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG2': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG3': self._DATA, self._MODE = False, 'BKG'
            case 'MC_BKG4': self._DATA, self._MODE = False, 'BKG'
            case 'MC_SIG': self._DATA, self._MODE = False, 'SIG'
            case _: raise ValueError(f'SAMP={SAMP} is not a valid option.')
        self._loaders[SAMP] = JPsiCCLoader(SAMP, self.YEAR, self.VERSION, self.CAT)
        self._loaders[SAMP].readSnapshot(filename, treename, read_pkl=True)
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
            KeyError: If any one of MC_BKG, MC_SIG, or DATA_BKG is missing.

        Returns:
            (None)
        '''
        samples = list(self._loaders.keys())
        mcbkgs = [self._loaders[s] for s in samples if 'MC_BKG' in s]
        if not mcbkgs: raise KeyError('No MC_BKG sample is found. Please load first.')
        if 'MC_SIG' not in samples: raise KeyError('No MC_SIG sample is found. Please load first.')
        if 'DATA_BKG' not in samples: raise KeyError('No DATA_BKG sample is found. Please load first.')
        mcsig = self._loaders['MC_SIG']
        databkg = self._loaders['DATA_BKG']

        hist_dict = databkg.makeHistos(plot=False, draw=False, keys=keys)

        for samp in samples:
            self._loaders[samp].plot_hists(hist_dict, draw=draw_indiv, user_sfx=user_sfx)
        
        ROOT.gStyle.SetOptStat(0)
        for key, hdef in hist_dict.items():
            # Create canvas
            c_hstack = ROOT.TCanvas()

            # Set legend
            legend = ROOT.TLegend(0.64, 0.78, 0.95, 0.95)
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
            print(f'MCBKG integral: {mcbkg_norm:.3f}, DATABKG integral: {data_hist.Integral():.3f}')
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
            nevents (float): Number of events in the range.
        '''
        for sample, loader in self._loaders.items():
            nevents = loader.countNevets(cut, weights)
        return 

    def applyCut(self, cut):
        '''Apply a selection to all the loaded samples. 

        Args:
            cut (str): RDF expression for filtering based on selection.

        Returns:
            (None)
        '''
        for sample, loader in self._loaders.items():
            print(f'Applying cut {cut} to {sample}.')
            loader.cut(f'{cut}')

    def scanCut(self, var, scantype):
        '''Scan the cut over a variable.

        Args:
            var (str): Name of the variable to scan over.
            scantype (str): Type of scanning. Options are the following.
                'floor': applies a lower limit
                'ceiling': applies an upper limit
                'window': applies a window

        Raises:
            ValueError: If the given scantype is not one of the options.
        '''
        # stats, vartype = 
        pass
    
    def scanMultiCut(self, varlist):
        '''Scan a selection of variables.

        Args:
            varlist (list(str)): List of variables to scan over.
        '''
        # for var in varlist:
        #     for sample
        pass

    def fitSignal(self, key_var, template):
        """Fit the signal RDF to a desired template.

        Args:
            key_var (str): Name of the variable to fit.
            template (str): Template to use for fitting.
        """
        fit = rootpdf.FittingTool(YEAR=self.YEAR, CAT=self.CAT, VARNAME=key_var)
        return

if __name__=='__main__':
    print('Choose preset (mcbkg, databkg, mcsig, allsnap, gen, allplots, noweight, cut): ', end='')
    preset = input()
    match preset:
        case 'mcbkg':
            mcbkg = JPsiCCLoader('MC_BKG', 2018, '202406', 'GF')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF()
            mcbkg.snapshotRDF()
            mcbkg.makeCutFlowPlot()

        case 'databkg':
            databkg = JPsiCCLoader('DATA_BKG', 2018, '202406', 'GF')
            # databkg.createDataRDF('event_spec_data_bkg_skim.json')
            # databkg.defineColumnsRDF()
            # databkg.snapshotRDF()
            databkg.readSnapshot('snapshot_DATA_BKG_2018_GF_v202406_20240731.root')
            databkg.makeCutFlowPlot()

        case 'mcsig':
            mcsig = JPsiCCLoader('MC_SIG', 2018, '202407', 'GF')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF()
            mcsig.snapshotRDF()
            mcsig.makeCutFlowPlot()
        
        case 'allsnap':
            mcbkg = JPsiCCLoader('MC_BKG', 2018, '202406', 'GF')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF()
            mcbkg.snapshotRDF()
            mcbkg.makeCutFlowPlot()
            
            databkg = JPsiCCLoader('DATA_BKG', 2018, '202406', 'GF')
            databkg.createDataRDF('event_spec_data_bkg_skim.json')
            databkg.defineColumnsRDF()
            databkg.snapshotRDF()
            databkg.makeCutFlowPlot()
            
            mcsig = JPsiCCLoader('MC_SIG', 2018, '202407', 'GF')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF()
            mcsig.snapshotRDF()
            mcsig.makeCutFlowPlot()

        case 'gen':
            mcsig = JPsiCCLoader('MC_SIG', 2018, '202406', 'GF')
            mcsig.readSnapshot('snapshot_MC_SIG_2018_GF_v202406_20240628.root')
            mcsig.makeHistos(genplots=True)

        case 'allplots' | 'noweight' | 'cut':
            if preset=='noweight': use_weight, draw_indiv = False, False
            else: use_weight, draw_indiv = True, True
            
            an = JPsiCCAnalyzer(2018, '202407', 'GF', weights=use_weight)
            an.readSnapshotSAMP('MC_BKG1', 'snapshot_MC_BKG_2018_GF_v202406_20240731.root',
                                sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            an.readSnapshotSAMP('MC_BKG2', 'snapshot_MC_BKG_2018_GF_v202406_20240731.root',
                                sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            an.readSnapshotSAMP('DATA_BKG', 'snapshot_DATA_BKG_2018_GF_v202406_20240731.root')
            an.readSnapshotSAMP('MC_SIG', 'snapshot_MC_SIG_2018_GF_v202407_20240731.root')

            if preset != 'cut': an.stackMultiHistos(draw_indiv=draw_indiv)
            else:
                an.applyCut('trigger_user>0')
                an.stackMultiHistos(draw_indiv=draw_indiv, user_sfx='trigger_user>0')