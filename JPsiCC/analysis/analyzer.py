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
ROOT.DisableImplicitMT()

class JPsiCCAnalyzer:
    '''JPsiCC analyzer.

    Args:
        SAMP (str): Either one of the following three options.
            'DATA_BKG', 'MC_BKG', 'MC_SIG'
        YEAR (int): Year of data-taking. Provide any number if using MC.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
    
    Raises:
        ValueError: If the string provided for SAMP is not among the three options.
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
    '''
    def __init__(self, SAMP, YEAR, VERS, CAT):
        match SAMP:
            case 'DATA_BKG': self._DATA, self._MODE = True, 'BKG'
            case 'MC_BKG': self._DATA, self._MODE = False, 'BKG'
            case 'MC_SIG': self._DATA, self._MODE = False, 'SIG'
            case _: raise ValueError(f'SAMP={SAMP} is not a valid option.')
        if not type(YEAR) is int: raise TypeError(f'YEAR must be an integer.')
        if not type(VERS) is str: raise TypeError(f'VERS must be a string.')
        if not type(CAT) is str: raise TypeError(f'CAT must be a string.')
        self._SAMPLE, self._YEAR, self._VERSION, self._CAT = SAMP, YEAR, VERS, CAT

        self._branches = []   # TBranches to take RDF snapshot of
        self._rdf = None      # Placeholder for RDF
        self._hists = dict()  # Dictionary of histogram objects
        self._models = dict() # Dictionary of histogram models

        today = date.today()
        self._date = f'{today.year}{today.month}{today.day}'
        self._anpath = 'JPsiCC'
        self._plotsavedir = os.path.join(os.environ['HRARE_DIR'], self._anpath, 'plots', self._VERSION, CAT)
        self._sfx = f'{SAMP}_{self._YEAR}_{self._CAT}_{self._VERSION}'

        # Color scheme: http://arxiv.org/pdf/2107.02270
        self._orange = ROOT.kOrange + 1
        self._blue = ROOT.kAzure - 4
        self._red = ROOT.kRed - 4
        self._purple = ROOT.kMagenta - 5
        self._gray = ROOT.kGray + 1
        self._violet = ROOT.kViolet + 2

        self._draw_option = 'P' if self._DATA else 'HIST'
    
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
                'DATA_BKG', 'MC_BKG', 'MC_SIG'.

        Returns:
            histo1d (ROOT.TH1D): The histogram.
        '''
        match SAMP:
            case 'DATA_BKG':
                histo1d.SetMarkerStyle(ROOT.kFullSquare)
                histo1d.SetMarkerSize(0.5)
            case 'MC_BKG':
                histo1d.SetFillColorAlpha(self._orange, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
            case 'MC_SIG':
                histo1d.SetFillColorAlpha(self._blue, 0.6)
                histo1d.SetLineColorAlpha(ROOT.kBlack, 1.0)
        return histo1d

    def __plot_hists(self, hist_dict, scaleSIG=1.):
        '''Internal method to plot histograms.

        Args:
            hist_dict (dict): Dictionary of histogram definitions.
            scaleSIG (float): Scale factor for the signal histogram.

        Returns:
            (None)
        '''
        for key, hdef in hist_dict.items():
            # Create Histo1D
            hname = f'{hdef["name"]}_{self._SAMPLE}_{self._YEAR}_{self._CAT}'
            model1d = (hname, hdef['title'], hdef['bin'], hdef['xmin'], hdef['xmax'])
            histo1d = self._rdf.Histo1D(model1d, hdef['name'], 'w')
            # Set color
            histo1d = self.__set_hist_style(histo1d, self._SAMPLE)
            if self._MODE=='SIG': histo1d.Scale(scaleSIG)
            # Save histogram
            histo1d.SetDirectory(0)
            self._hists[key] = histo1d
            self._models[key] = model1d
            # Draw histogram
            c = self.__draw_hist(histo1d=histo1d, model1d=model1d, draw_option=self._draw_option)
            c.SaveAs(os.path.join(self._plotsavedir, f'{hname}.png'))
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
        rdf_runs, _ = rdfdefines.rdf_def_sample_meta(rdf_runs)
        rdf_events, br1 = rdfdefines.rdf_def_sample_meta(rdf_events)
        if self._DATA: sum_weights = 1.
        else: sum_weights = rdfdefines.compute_sum_weights(rdf_runs)[0]
        self._rdf, br2 = rdfdefines.rdf_def_weights(rdf_events, sum_weights, data=self._DATA)
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
        match self._CAT:
            case 'GF':
                new_rdf, br1 = rdfdefines.rdf_filter_triggers(self._rdf, self._CAT, self._YEAR)
                new_rdf, br2 = rdfdefines.rdf_def_goodjets(new_rdf, self._CAT, self._YEAR)
                new_rdf, br3 = rdfdefines.rdf_def_jpsi(new_rdf)
                self._branches = self._branches + br1 + br2 + br3
            case _:
                new_rdf = self._rdf
        self._rdf = new_rdf
        return

    def retrieveRDF(self):
        '''Retrieve the internal RDF.

        Returns:
            rdf (ROOT.RDataFrame): The retrieved RDF.
        '''
        return self._rdf

    def snapshotRDF(self):
        '''Create a snapshot of the RDF.

        Outputs a file whose name is printed in the terminal.

        Returns:
            (None)
        '''
        fname = f'snapshot_{self._sfx}.root'
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

    def makeHistos(self):
        '''Make histograms from the internal RDF.

        Outputs histograms in the 'plots' directory.

        Returns:
            (None)
        '''
        hist_defs = jsonreader.get_object_from_json(anpath='JPsiCC',
                                                    jsonname='rdf_hists.json',
                                                    keys=['NANOAOD_to_RDF'])
        if not os.path.exists(self._plotsavedir): os.makedirs(self._plotsavedir)
        match self._CAT:
            case 'GF':
                self.__plot_hists(hist_defs['Jpsi'])
                # self.__plot_hists(hist_defs['jet'])
            case _: pass
        return
    
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
            analyzer (JPsiCCAnalyzer): another analyzer object.

        Returns:
            (None)
        '''
        this_SAMPLE, other_SAMPLE = self._SAMPLE, analyzer._SAMPLE
        for key, item in self._hists.items():
            other_hist = analyzer.retrieveHisto(key)
            if other_hist is not None:
                this_hist = self.__set_hist_style(item, this_SAMPLE)
                other_hist = self.__set_hist_style(other_hist, other_SAMPLE)
                # Create THStack
                hs = ROOT.THStack()
                hs.Add(this_hist.GetPtr(), self._draw_option)
                hs.Add(other_hist.GetPtr(), analyzer._draw_option)
                # Draw histograms
                model = self._models[key]
                c_hstack = self.__draw_hist(hs, model, draw_option='HIST NOSTACK')
                # Save TCanvas
                sfx = 'STACK_' + self._sfx.lstrip(f'{self._SAMPLE}_')
                c_hstack.SaveAs(os.path.join(self._plotsavedir, f'{key}_{sfx}.png'))
        return
    
# Draw
# save_dir = os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'plots', 'jetstudies', 'control')
# for hist in my_hists:
#     c = ROOT.TCanvas()
#     hist.Draw()
#     c.SaveAs(os.path.join(save_dir, f'{hist.GetTitle()}.png'))
#     c.Close()

if __name__=='__main__':
    mcbkg = JPsiCCAnalyzer('MC_BKG', 2018, '202405', 'GF')
    mcbkg.createWeightedRDF('run_spec_mc_bkg_2018.json', 'event_spec_mc_bkg_2018.json')
    mcbkg.defineColumnsRDF()
    mcbkg.snapshotRDF()
    # mcbkg.readSnapshot('snapshot_MC_BKG_2018_GF_202405.root')
    # mcbkg.makeHistos()

    databkg = JPsiCCAnalyzer('DATA_BKG', 2018, '202405', 'GF')
    databkg.createWeightedRDF('run_spec_data_bkg_2018.json', 'event_spec_data_bkg_2018.json')
    databkg.defineColumnsRDF()
    databkg.snapshotRDF()
    # databkg.readSnapshot('snapshot_DATA_BKG_2018_GF_202405.root')
    # databkg.makeHistos()

    mcbkg.stackHistos(databkg)