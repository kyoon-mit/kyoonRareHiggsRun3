'''
analyzer
+++++++++++++++++++++++++++++++++++++++++++

JPsiCC analysis-specific analyzer
'''

from kytools import jsonreader
from rdfdefines import rdf_def_weights, rdf_def_generic
import os, json
import ROOT

# Disable multithreading for now
ROOT.DisableImplicitMT()

class JPsiCCAnalyzer:
    '''JPsi analyzer.

    Args:
        DATA (bool): Whether to use DATA
        YEAR (int): Year of data-taking. Provide any number if using MC.
        VERSION (str): Version of the files.
        LUMI (float): Luminosity
    '''
    def __init__(self, DATA, YEAR, VERSION, LUMI, CAT):
        self._DATA = DATA
        self._YEAR = YEAR
        self._VERSION = VERSION
        self._LUMI = LUMI
        self._CAT = CAT
        self._anpath = 'JPsiCC'
        self._savedir = '.'
        self._chainBKG = None
        self._chainSIG = None
        self._rdfBKG = None
        self._rdfSIG = None
        self._sfxBKG = ''
        self._sfxCOMB = ''
        if self._DATA:
            self._sfxBKG = f'BKG_DATA_{self._YEAR}_lumi{self._LUMI}_{self._CAT}_{self._VERSION}'
            self._sfxCOMB = f'COMB_DATA_{self._YEAR}_lumi{self._LUMI}_{self._CAT}_{self._VERSION}'
        else:
            self._sfxBKG = f'BKG_MC_{self._YEAR}_lumi{self._LUMI}_{self._CAT}_{self._VERSION}'
            self._sfxCOMB = f'COMB_MC_{self._YEAR}_lumi{self._LUMI}_{self._CAT}_{self._VERSION}'
        self._sfxSIG = f'SIG_MC_{self._YEAR}_lumi{self._LUMI}_{self._CAT}_{self._VERSION}'
        # Color scheme: http://arxiv.org/pdf/2107.02270
        self._orange = ROOT.kOrange + 1
        self._blue = ROOT.kAzure - 4
        self._red = ROOT.kRed - 4
        self._purple = ROOT.kMagenta - 5
        self._gray = ROOT.kGray + 1
        self._violet = ROOT.kViolet + 2


    def __getfilesBKG(self, treename, xrtd_proxy):
        '''Internal method.
        '''
        if self._DATA:
            meta_json_name = 'data_names.json'
        else:
            meta_json_name = 'MC_bkg_names.json'
        meta_info = jsonreader.get_object_from_json(self._anpath, meta_json_name, ['kraken', self._VERSION, self._YEAR])
        chain = ROOT.TChain(treename)
        for info in meta_info:
            events = jsonreader.get_chain_from_json_xrtd(anpath=self._anpath,
                                                         jsonname='NANOAOD.json',
                                                         keys=['kraken', self._VERSION, info['dataset']],
                                                         treename=treename,
                                                         xrtd_proxy=xrtd_proxy)
            chain.Add(events)
        return chain
    
    def __getfilesSIG(self, treename):
        '''Internal method.
        '''
        chain = ROOT.TChain(treename)
        events = jsonreader.get_chain_from_json(anpath=self._anpath,
                                                jsonname='NANOAOD.json',
                                                keys=['mariadlf', '202405', 'GEN-signal'],
                                                treename=treename)
        chain.Add(events)
        return chain
    
    def __createWeightedRDF(self, treename, nanoaodjson_key, dataset, xsec, xsec_sigma, xrtd_proxy=''):
        '''Internal method.

        Open a dataset and load it onto a RDF with proper weights.
        Each dataset will have different cross sections, and this is the base-level
        method to load a dataset with a single cross section.

        Args:
            treename (str): Name of the TTree in the target ROOT files.
            nanoaodjson_key (str): Top-level key in the NANOAOD.json.
                e.g. 'mariadlf', 'kraken', etc.
            dataset (str): Name of the dataset.
                For gen-level signal, provide 'GEN-signal.'
            xsec (float): The cross section of this dataset.
            xsec_sigma (float): The uncertainty on the cross section.
            xrtd_proxy (str, optional): If provided, it will open a remote file with the proxy.
                Otherwise, it will open a local file.
                Defaults to ''.
        '''
        filenames = jsonreader.get_filenames_from_json(anpath=self._anpath,
                                                       jsonname='NANOAOD.json',
                                                       keys=[nanoaodjson_key, self._VERSION, dataset],
                                                       xrtd_proxy=xrtd_proxy)
        rdf = ROOT.RDataFrame(treename, filenames)
        rdf = rdf_def_weights(rdf=rdf, lumi=self._LUMI, data=self._DATA, xsec=xsec, xsec_sigma=xsec_sigma)
        return rdf

    def __draw_hist(self, hist, model1d, draw_option='HIST'):
        '''Internal method.
        '''
        c = ROOT.TCanvas()
        bin_width = (model1d[4]-model1d[3])/model1d[2]
        hist.Draw(draw_option)
        hist.SetTitle(model1d[1])
        hist.GetXaxis().SetTitle(model1d[1])
        hist.GetYaxis().SetTitle(f'Events / {bin_width:.3g}')
        c.Update()
        return c

    def getfilesBKG(self, treename, xrtd_proxy):
        '''Retrieve files.

        Args:
            treename (str): Name of the TTree.
            xrtd_proxy (str): Name of the XRootD proxy.

        Returns:
            (None)
        '''
        self._chainBKG = self.__getfilesBKG(treename, xrtd_proxy)

    def getfilesSIG(self, treename):
        '''Retrieve files.

        Args:
            treename (str): Name of the TTree.

        Returns:
            (None)
        '''
        self._chainSIG = self.__getfilesSIG(treename)

    def getfiles_local(self, treename, filenames, signal=True):
        '''Retrieve files locally.

        Args:
            treename (str): Name of the TTree.
            filenames (list(str)): Names of the files.
            signal (bool): Whether to retrieve signal files.
                Defaults to True.

        Returns:
            (None)
        '''
        chain = ROOT.TChain(treename)
        for name in filenames:
            chain.Add(name)
        if signal: self._chainSIG = chain
        else: self._chainBKG = chain

    def doBKG(self):
        '''Do the background analysis.

        Returns:
           (None)
        '''
        if self._chainBKG is None:
            raise Exception("First open the files using one of the \'gefiles\' methods.")
        rdf = ROOT.RDataFrame(self._chainBKG)
        self._rdfBKG = rdf_def_generic(rdf, self._LUMI, self._DATA)

    def doSIG(self):
        '''Do the signal analysis.

        Args:
            scale (float): Scale the signal by this factor.

        Returns:
           (None)
        '''
        if self._chainSIG is None:
            raise Exception("First open the files using one of the \'getfiles\' methods.")
        rdf = ROOT.RDataFrame(self._chainSIG)
        self._rdfSIG = rdf_def_generic(rdf, self._LUMI, data=False)

    def snapshot(self, branches, outname, signal=True):
        '''Create snapshot of the RDF.

        Args:
            branches (list(str)): Name of the branches to create a snapshot.
            outname (str): Name of the output file.
            signal (bool): Whether to retrieve signal files.
                Defaults to True.

        Returns:
            (None)
        '''
        pass
    
    def makehists(self, scaleSIG=1.):
        '''Make histograms.

        Args:
            scaleSIG (float): Scale factor for the signal histogram.

        Returns:
            (None)
        '''
        if self._rdfBKG is None:
            raise Exception("First run \'doBKG\'.")
        hlist_BKG, hlist_SIG = dict(), dict()
        hist_defs = jsonreader.get_object_from_json(anpath='JPsiCC',
                                                    jsonname='rdf_hists.json',
                                                    keys=['NANOAOD_to_RDF'])
        savedir = os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'plots', self._savedir)
        if not os.path.exists(savedir): os.makedirs(savedir)
        for key, hdef in hist_defs['test'].items():
            # Create Histo1D
            model1d = (f'{hdef["name"]}_{self._YEAR}', hdef['title'], hdef['bin'], hdef['xmin'], hdef['xmax'])
            hbkg = self._rdfBKG.Histo1D(model1d, hdef['name'], 'w')
            hsig = self._rdfSIG.Histo1D(model1d, hdef['name'], 'w')

            # Set color
            hbkg.SetFillColorAlpha(self._orange, 0.6)
            hsig.SetFillColorAlpha(self._blue, 0.6)

            # Scale signal
            hsig.Scale(scaleSIG)

            # Create THStack
            hs = ROOT.THStack()
            hs.Add(hbkg.GetPtr())
            hs.Add(hsig.GetPtr())

            # Draw histograms
            c_bkg = self.__draw_hist(hbkg, model1d)
            c_sig = self.__draw_hist(hsig, model1d)
            c_hstack = self.__draw_hist(hs, model1d, draw_option='HIST NOSTACK')

            # Save TCanvas
            c_bkg.SaveAs(os.path.join(savedir, f'{key}_{self._sfxBKG}.png'))
            c_sig.SaveAs(os.path.join(savedir, f'{key}_{self._sfxSIG}.png'))
            c_hstack.SaveAs(os.path.join(savedir, f'{key}_stack_{self._sfxCOMB}.png'))
    
# Draw
# save_dir = os.path.join(os.environ['HRARE_DIR'], 'JPsiCC', 'plots', 'jetstudies', 'control')
# for hist in my_hists:
#     c = ROOT.TCanvas()
#     hist.Draw()
#     c.SaveAs(os.path.join(save_dir, f'{hist.GetTitle()}.png'))
#     c.Close()

if __name__=='__main__':
    analyzer = JPsiCCAnalyzer(DATA=False, YEAR=2018, VERSION='test', LUMI=100, CAT='generic')
    # analyzer.getfiles('Events', 'root://bost-cms-xcache01.lhcone.es.net:1094//') # xcache server down
    analyzer.getfilesBKG('Events', 'root://xrootd.cmsaf.mit.edu//')
    analyzer.getfilesSIG('Events')
    analyzer.doBKG()
    analyzer.doSIG()
    analyzer.makehists(scaleSIG=10e4)