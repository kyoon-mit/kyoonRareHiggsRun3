'''
limits
+++++++++++++++++++++++++++++++++++++++++++

Library to compute limits
'''

import os
from datetime import date
import ROOT

ROOT.EnableImplicitMT()

class RooWorkspaceCreator:
    '''
    RooWorkspace creator.

    To create a RooWorkspace for a given sample.

    Args:
        SAMP (str): Either one of the following options.
            'DATA_BKG', 'MC_BKG', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        weights (bool): Whether weights were used.

    Raises:
        ValueError: If the string provided for SAMP is not among the options.
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
    '''
    def __init__(self, SAMP, YEAR, VERS, CAT, CMSSW, weights=True):
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
        self.SAMPLE, self.YEAR, self.VERSION, self.CAT, self.CMSSW = SAMP, YEAR, VERS, CAT, CMSSW
        self._weights = weights # bool

        today = date.today()
        self._date = f'{today.year}{today.month:02}{today.day:02}'
        self._anpath = 'JPsiCC'
        self._plotsavedir = os.path.join(os.environ['HRARE_DIR'], self._anpath, 'plots', f'v{self.VERSION}', self._date, CAT)
        self._sfx = f'{SAMP}_{self.YEAR}_{self.CAT}_v{self.VERSION}_{self.CMSSW}_{self._date}_{"WEIGHT" if weights else "NOWEIGHT"}'
        self._wname = f'w' # Name of the ROOT workspace
        self._wfilename = f'workspace_{self._sfx}.root' # Name of the ROOT workspace file
        w = ROOT.RooWorkspace(self._wname)
        w.writeToFile(self._wfilename)
        print('{}INFO: a workspace file has been created. >> {} {}'.format('\033[1;33m', self._wfilename, '\033[0m'))

    def __importToWorkspace(self, *args):
        '''Internal method for importing objects into the workspace.

        Args:
            *args (TObject): The objects to import are passed through here.

        Returns:
            (None)
        '''
        wfile = ROOT.TFile.Open(self._wfilename, 'UPDATE')
        w = wfile.Get(self._wname)
        for arg in args:
            w.Import(arg)
        w.writeToFile(self._wfilename, recreate=True)
        wfile.Close()
        return

    def __readVarFromWorkspace(self, var_name):
        '''Internal method for reading an variable from the workspace.

        Args:
            var_name (str): The name of the variable to read.

        Returns:
            var (ROOT.RooRealVar): The variable.
        '''
        wfile = ROOT.TFile.Open(self._wfilename, 'UPDATE')
        w = wfile.Get(self._wname)
        var = w.var(var_name)
        return var

    def addVar(self, var_name, var_title, value, var_min, var_max, unit=''):
        '''Add a variable to the workspace.

        Args:
            var_name (str): Variable name.
            var_title (str): Variable title.
            value (float): Nominal value.
            var_min (float): Lower bound.
            var_max (float): Upper bound.
            unit (str, optional): Unit.

        Returns:
            (None)

        Raises:
            ValueError: If var_min >= var_max.
        '''
        if not (var_min < var_max):
            raise ValueError('Variable lower bound should be smaller than upper bound.')
        var = ROOT.RooRealVar(var_name, var_title, value, var_min, var_max, unit)
        self.__importToWorkspace(var)
        return

    def addData(self, filename, treename, col_name, var_name, binned=True):
        '''Add data to the workspace.

        If binned, it will add a RooDataHist. Otherwise, it will add a RooDataSet.
        For reference, see https://root.cern/doc/v632/rf408__RDataFrameToRooFit_8py_source.html.

        Args:
            filename (str): Name of the file containing the data.
            treename (str): Name of the ROOT.TTree.
            col_name (str): Name of the column of the data.
            var_name (str): Name of the variable as appears in the workspace.
        
        Returns:
            (None)
        '''
        rdf = ROOT.RDataFrame(treename, filename)
        var = self.__readVarFromWorkspace(var_name)
        data_name = f'data{"hist" if binned else "set"}_{self.SAMPLE}_{self.YEAR}_{self.CAT}'
        data_title = data_name
        if binned:
            data_maker = ROOT.RooDataHistHelper(data_name, data_title, ROOT.RooArgSet(var))
        else:
            data_maker = ROOT.RooDataSetHelper(data_name, data_title, ROOT.RooArgSet(var))
        data_result = rdf.Book(ROOT.std.move(data_maker), [col_name])
        data = data_result.GetValue()
        self.__importToWorkspace(data)
        return

    def defineRegions(self, SR_low, SR_high, CR_low, CR_high):
        '''Define signal and control regions for the key observable.
        TODO: make regions in workspace.

        Args:
            SR_low (float): Lower bound of the signal region.
            SR_high (float): Upper bound of the signal region.
            CR_low (float): Lower bound of the control region.
            CR_high (float): Upper bound of the control region.
        
        Returns:
            (None)

        Raises:
            TypeError: If input arguments are not float.
            ValueError: If lower bound >= upper bound for any of the inputs.
            ValueError: If signal region is outside the control region.
        '''
        if not all(isinstance(arg, float) for arg in list(locals.values())):
            raise TypeError('Input arguments must be float.')
        elif not (SR_low < SR_high):
            raise ValueError('SR_low must be smaller than SR_high.')
        elif not (CR_low < CR_high):
            raise ValueError('CR_low must be smaller than CR_high.')
        elif not (CR_low < SR_low and SR_high < CR_high):
            raise ValueError('Signal region must lie outside the control region.')
        self._SR_low, self._SR_high, self._CR_low, self._CR_high = SR_low, SR_high, CR_low, CR_high
        return

    def specifyFit():
        '''Specify what fit to create in the workspace.

        Args:

        Returns:

        '''
        pass

if __name__=='__main__':
    verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)
    snapshot_dir = '/work/submit/kyoon/CMSSW_13_3_0/src/RareHiggsRun3/JPsiCC/analysis'
    snapshot_name = os.path.join(snapshot_dir, 'snapshot_MC_SIG_2018_GF_v202407_20241024_WEIGHT_no_filter.root')
    rwc = RooWorkspaceCreator('MC_BKG', 2018, '202407', 'GF', 'CMSSW_13_3_0')
    rwc.addVar('mH', 'm_{#mu#bar{#mu}jj}', 125, 30, 200, 'GeV/c^2')
    rwc.addData(snapshot_name, 'Events', 'higgs_mass_corr', 'mH')