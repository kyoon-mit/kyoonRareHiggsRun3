'''
analyzer
+++++++++++++++++++++++++++++++++++++++++++

JPsiCC analysis-specific analyzer
'''

from kytools import jsonreader
from datetime import date
import rdfdefines
import os, json
import ROOT

# Disable multithreading for now
# ROOT.DisableImplicitMT()
ROOT.EnableImplicitMT()

class JPsiCCLoader:
    '''JPsiCC loader.

    Useful for loading one sample at a time, or for stacking histograms.

    Args:
        SAMP (str): Either one of the following three options.
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
        self._weights = weights
        self.SAMPLENAME = f'{self.SAMPLE}'

        self._branches = []   # TBranches to take RDF snapshot of
        self._rdf = None      # Placeholder for RDF
        self._hists = dict()  # Dictionary of histogram objects
        self._models = dict() # Dictionary of histogram models

        today = date.today()
        self._date = f'{today.year}{today.month:02}{today.day:02}'
        self._anpath = 'JPsiCC'
        self._plotsavedir = os.path.join(os.environ['HRARE_DIR'], self._anpath, 'plots', f'v{self.VERSION}', self._date, CAT)
        self._sfx = f'{SAMP}_{self.YEAR}_{self.CAT}_v{self.VERSION}_{self._date}'
        if not weights: self._sfx += '_NOWEIGHT'

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

    def plot_hists(self, hist_dict, draw=True, scaleSIG=1.):
        '''Plot histograms from dictionary of definitions.

        Args:
            hist_dict (dict): Dictionary of histogram definitions.
            draw (bool, optional): Whether to draw histogram.
                Defaults to True.
            scaleSIG (float, optional): Scale factor for the signal histogram.
                Defaults to 1.

        Returns:
            (None)
        '''
        for key, hdef in hist_dict.items():
            # Create Histo1D
            hname = f'{hdef["name"]}_{self._sfx}'
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
        rdf_events = jsonreader.get_rdf_from_json_spec(self._anpath, event_spec_json)
        rdf_events, br1 = rdfdefines.rdf_defSAMPLE_meta(rdf_events)
        self._rdf, br2 = rdfdefines.rdf_def_weights(None, rdf_events, data=self._DATA)
        self._branches = self._branches + br1 + br2
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
        rdf_events = jsonreader.get_rdf_from_json_spec(self._anpath, event_spec_json)
        rdf_runs, _ = rdfdefines.rdf_defSAMPLE_meta(rdf_runs)
        rdf_events, br1 = rdfdefines.rdf_defSAMPLE_meta(rdf_events)
        self._rdf, br2 = rdfdefines.rdf_def_weights(rdf_runs, rdf_events, data=self._DATA)
        self._branches = self._branches + br1 + br2
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
                new_rdf, br1 = rdfdefines.rdf_filter_triggers(self._rdf, self.CAT, self.YEAR)
                new_rdf, br2 = rdfdefines.rdf_def_jpsi(new_rdf)
                new_rdf, br3 = rdfdefines.rdf_def_muons(new_rdf)
                new_rdf, br4 = rdfdefines.rdf_def_vertex(new_rdf)
                new_rdf, br5 = rdfdefines.rdf_def_jets(new_rdf, self.CAT, self.YEAR, data=self._DATA)
                if self.SAMPLE=='MC_SIG':
                    new_rdf, br_gen = rdfdefines.rdf_def_genpart(new_rdf)
                    self._branches += br_gen
                self._branches = self._branches + br1 + br2 + br3 + br4 + br5
            case _:
                new_rdf = self._rdf
        self._rdf = new_rdf
        return
    
    def simpleCut(self, cut):
        '''Filter the internal RDF using the cut expression.

        Args:
            cut (str): String expression for the cut filter.

        Returns:
            (None)
        '''
        self._rdf = self._rdf.Filter(cut)
    
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

    def snapshotRDF(self, name=''):
        '''Create a snapshot of the RDF.

        Outputs a file whose name is printed in the terminal.

        Args:
            name (str, optional): Name of the snapshot to create.
                Defaults to ''. If '', then the snapshot name will be determined automatically.

        Returns:
            (None)
        '''
        if name=='': fname = f'snapshot_{self._sfx}.root'
        else: fname = name
        self._rdf.Snapshot('Events', fname, self._branches)
        print('{}INFO: a snapshot has been created. >> {} {}'.format('\033[1;33m', fname, '\033[0m'))
        return
    
    def readSnapshot(self, filename, treename='Events'):
        '''Retrieve RDF from a saved snapshot ROOT file.

        Args:
            filename (str): Name of the snapshot ROOT file.
            treename (str, option): Name of the TTree in the ROOT file.

        Returns:
            (None)
        '''
        self._rdf = ROOT.RDataFrame(treename, filename)
        return

    def makeHistos(self, plot=True, draw=True, genplots=False, keys=[], scaleSIG=1.):
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

        Returns:
            (None) or hist_dict (dict)
        '''
        hist_defs = jsonreader.get_object_from_json(anpath='JPsiCC',
                                                    jsonname='rdf_hists.json',
                                                    keys=['NANOAOD_to_RDF'])
        if not os.path.exists(self._plotsavedir): os.makedirs(self._plotsavedir)
        hist_dict = dict()
        match self.CAT:
            case 'GF':
                hist_dict.update(hist_defs['Jpsi'])
                hist_dict.update(hist_defs['muon'])
                hist_dict.update(hist_defs['vertex'])
                if genplots: hist_dict.update(hist_defs['gen'])
                # self.plot_hists(hist_defs['jet'])
            case _: pass
        if keys: hist_dict = {key: val for key, val in hist_dict.items() if key in keys}
        if plot: self.plot_hists(hist_dict, draw=draw, scaleSIG=scaleSIG)
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
        self._loaders[SAMP].readSnapshot(filename, treename)
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

        hist_dict = mcsig.makeHistos(plot=False, draw=False, keys=keys)

        for samp in samples:
            self._loaders[samp].plot_hists(hist_dict, draw=draw_indiv)
        
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
                if h is None: continue
                h, mcbkg_d_opt = self.__set_hist_style(h, loader.SAMPLE)
                mcbkg_norm += h.Integral()
                mcbkg_max += h.GetMaximum()
                hs.Add(h.GetPtr(), mcbkg_d_opt)
                legend.AddEntry(h.GetPtr(), f'{loader.SAMPLENAME} ({h.Integral():.2E})', 'F')
            hs.Draw()

            # Draw data
            data_hist = databkg.retrieveHisto(key)
            if data_hist is None: continue
            data_hist, databkg_d_opt = self.__set_hist_style(data_hist, databkg.SAMPLE)
            data_hist.Draw(databkg_d_opt)
            legend.AddEntry(data_hist.GetPtr(), f'{databkg.SAMPLENAME} ({data_hist.Integral():.2E})', 'P')

            # Draw signal
            sig_hist = mcsig.retrieveHisto(key)
            if sig_hist is None: continue
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

    def select(self, cut, draw=False):
        '''

        Args:
            cut (str):
            draw (bool, optional): Whether to draw histograms.
                Defaults to False.
        '''
        for sample, loader in self._loaders.items():
            print(f'Applying cut {cut} to {sample}.')
            loader.simpleCut(cut)
        if draw:
            return
    
    def scanSelection(self, variable, scansize):
        pass

if __name__=='__main__':
    print('Choose preset (mcbkg, databkg, mcsig, allsnap, gen, allplots, noweight): ', end='')
    preset = input()
    match preset:
        case 'mcbkg':
            mcbkg = JPsiCCLoader('MC_BKG', 2018, '202406', 'GF')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF()
            mcbkg.snapshotRDF()

        case 'databkg':
            databkg = JPsiCCLoader('DATA_BKG', 2018, '202406', 'GF')
            databkg.createDataRDF('event_spec_data_bkg_skim.json')
            databkg.defineColumnsRDF()
            databkg.snapshotRDF()

        case 'mcsig':
            mcsig = JPsiCCLoader('MC_SIG', 2018, '202406', 'GF')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF()
            mcsig.snapshotRDF()
        
        case 'allsnap':
            mcbkg = JPsiCCLoader('MC_BKG', 2018, '202406', 'GF')
            mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
            mcbkg.defineColumnsRDF()
            mcbkg.snapshotRDF()
            
            databkg = JPsiCCLoader('DATA_BKG', 2018, '202406', 'GF')
            databkg.createDataRDF('event_spec_data_bkg_skim.json')
            databkg.defineColumnsRDF()
            databkg.snapshotRDF()
            
            mcsig = JPsiCCLoader('MC_SIG', 2018, '202406', 'GF')
            mcsig.createWeightedRDF('run_spec_mc_sig.json', 'event_spec_mc_sig.json')
            mcsig.defineColumnsRDF()
            mcsig.snapshotRDF()

        case 'gen':
            mcsig = JPsiCCLoader('MC_SIG', 2018, '202406', 'GF')
            mcsig.readSnapshot('snapshot_MC_SIG_2018_GF_v202406_20240628.root')
            mcsig.makeHistos(genplots=True)

        case 'allplots' | 'noweight':
            if preset=='noweight': use_weight, draw_indiv = False, False
            else: use_weight, draw_indiv = True, True
            
            an = JPsiCCAnalyzer(2018, '202406', 'GF', weights=use_weight)
            an.readSnapshotSAMP('MC_BKG1', 'snapshot_MC_BKG_2018_GF_v202406_20240712.root',
                                sample_name='BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM')
            an.readSnapshotSAMP('MC_BKG2', 'snapshot_MC_BKG_2018_GF_v202406_20240712.root',
                                sample_name='JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM')
            an.readSnapshotSAMP('DATA_BKG', 'snapshot_DATA_BKG_2018_GF_v202406_20240710.root')
            an.readSnapshotSAMP('MC_SIG', 'snapshot_MC_SIG_2018_GF_v202406_20240710.root')

            an.stackMultiHistos(draw_indiv=draw_indiv, keys=['Jpsi_kin_mass'])