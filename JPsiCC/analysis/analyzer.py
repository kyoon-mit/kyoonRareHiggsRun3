'''
analyzer
+++++++++++++++++++++++++++++++++++++++++++

JPsiCC analysis-specific analyzer
'''

from kytools import jsonreader
from rdfdefines import rdf_def_generic
import os, json
import ROOT

# Disable multithreading for now
ROOT.DisableImplicitMT()

class JPsiCCAnalyzer:
    """JPsi analyzer.

    Args:
        DATA (bool): Whether to use DATA
        YEAR (int): Year of data-taking. Provide any number if using MC.
        VERSION (str): Version of the files.
        LUMI (float): Luminosity
    """
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

    def __getfilesBKG(self, treename, xrtd_proxy):
        """Internal method.
        """
        if self._DATA:
            meta_json_name = 'data_names.json'
        else:
            meta_json_name = 'MC_bkg_names.json'
        meta_names = jsonreader.get_object_from_json(self._anpath, meta_json_name, ['kraken', self._VERSION])
        chain = ROOT.TChain(treename)
        for name in meta_names:
            events = jsonreader.get_chain_from_json_xrtd(anpath=self._anpath,
                                                         jsonname='NANOAOD.json',
                                                         keys=['kraken', self._VERSION, name],
                                                         treename=treename,
                                                         xrtd_proxy=xrtd_proxy)
            chain.Add(events)
        return chain
    
    def __getfilesSIG(self, treename):
        """Internal method.
        """
        chain = ROOT.TChain(treename)
        events = jsonreader.get_chain_from_json(anpath=self._anpath,
                                                jsonname='NANOAOD.json',
                                                keys=['mariadlf', '202405', 'GEN-signal'],
                                                treename=treename)
        chain.Add(events)
        return chain
    
    def getfilesBKG(self, treename, xrtd_proxy):
        """Retrieve files.

        Args:
            treename (str): Name of the TTree.
            xrtd_proxy (str): Name of the XRootD proxy.

        Returns:
            (None)
        """
        self._chainBKG = self.__getfilesBKG(treename, xrtd_proxy)

    def getfilesSIG(self, treename):
        """Retrieve files.

        Args:
            treename (str): Name of the TTree.

        Returns:
            (None)
        """
        self._chainSIG = self.__getfilesSIG(treename)

    def getfiles_local(self, treename, filenames, signal=True):
        """Retrieve files locally.

        Args:
            treename (str): Name of the TTree.
            filenames (list(str)): Names of the files.
            signal (bool): Whether to retrieve signal files.
                Defaults to True.

        Returns:
            (None)
        """
        chain = ROOT.TChain(treename)
        for name in filenames:
            chain.Add(name)
        if signal: self._chainSIG = chain
        else: self._chainBKG = chain

    def doBKG(self):
        """Do the background analysis.

        Returns:
           (None)
        """
        if self._chainBKG is None:
            raise Exception("First open the files using one of the \'gefiles\' methods.")
        rdf = ROOT.RDataFrame(self._chainBKG)
        self._rdfBKG = rdf_def_generic(rdf, self._LUMI, self._DATA)

    def doSIG(self):
        """Do the signal analysis.

        Returns:
           (None)
        """
        if self._chainSIG is None:
            raise Exception("First open the files using one of the \'gefiles\' methods.")
        rdf = ROOT.RDataFrame(self._chainSIG)
        self._rdfSIG = rdf_def_generic(rdf, self._LUMI, data=False)

    def snapshot(self, branches, outname, signal=True):
        """Create snapshot of the RDF.

        Args:
            branches (list(str)): Name of the branches to create a snapshot.
            outname (str): Name of the output file.
            signal (bool): Whether to retrieve signal files.
                Defaults to True.

        Returns:
            (None)
        """
        pass
    
    def makehistsBKG(self):
        """Make histograms.

        Returns:
            (None)
        """
        if self._rdfBKG is None:
            raise Exception("First run \'doBKG\'.")
        hist_list = dict()
        hist_defs = jsonreader.get_object_from_json(anpath='JPsiCC',
                                                    jsonname='rdf_hists.json',
                                                    keys=['NANOAOD_to_RDF'])
        for key, hdef in hist_defs['test'].items():
            model1d = (f'{hdef["name"]}_{self._YEAR}', hdef['title'], hdef['bin'], hdef['xmin'], hdef['xmax'])
            print(model1d)
            hist_list[key] = self._rdfBKG.Histo1D(model1d, hdef['name'], 'w')
            c = ROOT.TCanvas()
            hist_list[key].Draw()
            c.SaveAs(os.path.join(self._savedir, f'{key}.png'))
            c.Close()

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
    analyzer.doBKG()
    analyzer.makehistsBKG()