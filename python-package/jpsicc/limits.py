'''
limits
+++++++++++++++++++++++++++++++++++++++++++

Module to compute limits.
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
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        weights (bool): Whether weights were used.
        suffix (str, optional): Suffix of the workspace file name.

    Raises:
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
        TypeError: If the value provided for suffix is not a string.
    '''
    def __init__(self, YEAR: int, VERS: str, CAT: str, CMSSW: str, weights=True, suffix=''):
        if not isinstance(YEAR, int): raise TypeError(f'YEAR must be an integer.')
        if not isinstance(VERS, str): raise TypeError(f'VERS must be a string.')
        if not isinstance(CAT, str): raise TypeError(f'CAT must be a string.')
        if not isinstance(suffix, str): raise TypeError(f'The value for suffix must be a string.')
        self.YEAR, self.VERSION, self.CAT, self.CMSSW = YEAR, VERS, CAT, CMSSW
        self._weights = weights # bool

        today = date.today()
        self._date = f'{today.year}{today.month:02}{today.day:02}'
        self._anpath = 'JPsiCC'
        self._plotsavedir = os.path.join(os.environ['HRARE_DIR'], self._anpath, 'plots', f'v{self.VERSION}', self._date, CAT)
        if not os.path.exists(self._plotsavedir): os.makedirs(self._plotsavedir)
        self._plot_configs = {}
        self._sfx = f'{self.YEAR}_{self.CAT}_v{self.VERSION}_{self.CMSSW}_{self._date}_{"WEIGHT" if weights else "NOWEIGHT"}'
        self._wname = f'w' # Name of the ROOT workspace
        self._wfilename = f'workspace_{self._sfx}{"_" if suffix else ""}{suffix}.root' # Name of the ROOT workspace file
        self._SR_low, self._SR_high, self._CR_low, self._CR_high = 120, 130, 110, 140 # default values
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

    def __readFromWorkspace(self, obj_name: str, obj_type: str):
        '''Internal method for reading an object from the workspace.

        Args:
            obj_name (str): The name of the object to read.
            obj_type (str): The type of the object to read. The options are:
                'var', 'pdf', or 'data'.

        Returns:
            object (ROOT.RooRealVar, ROOT.RooAbsPdf, or ROOT.RooAbsData): The object.

        Raises:
            ValueError: If obj_type is not one of the following:
                'var', 'pdf', or 'data'.
        '''
        wfile = ROOT.TFile.Open(self._wfilename, 'UPDATE')
        w = wfile.Get(self._wname)
        if obj_type=='var':
            obj = w.var(obj_name)
        elif obj_type=='pdf':
            obj = w.pdf(obj_name)
        elif obj_type=='data':
            obj = w.data(obj_name)
        else:
            raise ValueError('Invalid obj_type.')
        return obj

    def addVar(self, var_name: str, var_title: str, value: float, var_min: float, var_max: float, unit=''):
        '''Add a variable to the workspace.

        Args:
            var_name (str): Variable name.
            var_title (str): Variable title.
            value (float): Nominal value.
            var_min (int or float): Lower bound.
            var_max (int or float): Upper bound.
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

    def addDataSet(self, SAMP, filename: str, treename: str, var_name: str, col_name: str, weight_name=''):
        '''Add dataset to the workspace.

        For reference, see https://root.cern/doc/v632/rf408__RDataFrameToRooFit_8py_source.html.

        Args:
            SAMP (str): Name of the sample.
            filename (str): Name of the file containing the data.
            treename (str): Name of the ROOT.TTree.
            var_name (str): Name of the variable as appears in the workspace.
            col_name (str): Name of the column of the data.
            weight_name (str, optional): Name of the weight variable as appears in the workspace.
                If none specified, the data will be unweighted. Defaults to ''.
        
        Returns:
            (None)
        '''
        rdf = ROOT.RDataFrame(treename, filename)
        var = self.__readFromWorkspace(var_name, 'var')
        data_name = f'{SAMP}_{self.YEAR}_{self.CAT}_dataset'
        data_title = data_name
        if weight_name:
            rdf = rdf.Redefine(weight_name, f'float({weight_name})')
            var_weight = ROOT.RooRealVar(weight_name, weight_name, -1e6, 1e6)
            rooargset = ROOT.RooArgSet(var, var_weight)
            cols = [col_name, weight_name]
        else:
            rooargset = ROOT.RooArgSet(var)
            cols = [col_name]
        dataset_maker = ROOT.RooDataSetHelper(data_name, data_title, rooargset)
        dataset_result = rdf.Book(ROOT.std.move(dataset_maker), cols)
        data = dataset_result.GetValue() # TODO: needs to be fixed; want weighted dataset
        if weight_name:
            data = ROOT.RooDataSet(data, weight_name)
        self.__importToWorkspace(data)
        return
    
    def addDataHist(self, SAMP, filename: str, treename: str, var_name: str, col_name: str, nbins: int, weight_name=''):
        '''Add datahist to workspace.

        Args:
            SAMP (str): Name of the sample.
            filename (str): Name of the file containing the data.
            treename (str): Name of the ROOT.TTree.
            var_name (str): Name of the variable as appears in the workspace.
            col_name (str): Name of the column of the data.
            nbins (int): Number of bins.
            weight_name (str, optional): Name of the weight variable as appears in the workspace.
                If none specified, the data will be unweighted. Defaults to ''.

        Returns:
            (None)
        '''
        rdf = ROOT.RDataFrame(treename, filename)
        var = self.__readFromWorkspace(var_name, 'var')
        data_name = f'{SAMP}_{self.YEAR}_{self.CAT}_datahist'
        data_title = data_name
        histo1d = rdf.Histo1D((data_name, data_title, nbins, var.getBinning().lowBound(), var.getBinning().highBound()), col_name, weight_name)
        datahist = ROOT.RooDataHist(ROOT.RooStringView(data_name), ROOT.RooStringView(data_title), ROOT.RooArgList(var), histo1d.GetValue())
        self.__importToWorkspace(datahist)
        return

    def defineRegions(self, SR_low: float, SR_high: float, CR_low: float, CR_high: float):
        '''Define signal and control regions for the key observable.
        TODO: make regions in workspace.

        Args:
            SR_low (int or float): Lower bound of the signal region.
            SR_high (int or float): Upper bound of the signal region.
            CR_low (int or float): Lower bound of the control region.
            CR_high (int or float): Upper bound of the control region.
        
        Returns:
            (None)

        Raises:
            TypeError: If input arguments are not int or float.
            ValueError: If lower bound >= upper bound for any of the inputs.
            ValueError: If signal region is outside the control region.
        '''
        if not all(isinstance(arg, (int, float)) for arg in list(locals().values())[1:]):
            raise TypeError('Input arguments must be float.')
        elif not (SR_low < SR_high):
            raise ValueError('SR_low must be smaller than SR_high.')
        elif not (CR_low < CR_high):
            raise ValueError('CR_low must be smaller than CR_high.')
        elif not (CR_low <= SR_low and SR_high <= CR_high):
            raise ValueError('Signal region must lie outside or on the boundary of control region.')
        self._SR_low, self._SR_high, self._CR_low, self._CR_high = SR_low, SR_high, CR_low, CR_high
        return

    def addPDF(self, SAMP: str, pdf_type: str, var_name: str, max_tries=12, strategy=1, binned=True):
        '''Add PDF to the workspace.

        The PDF type is specified. The PDF is fitted to the existing data.

        Args:
            SAMP (str): Name of the sample.
            pdf_type (str): PDF type, e.g. 'gaussian', 'double_gaussian', 'crytal_ball', etc.
                See code body for exhaustive list of PDFs.
            var_name (str):
            max_tries (int, optional):
            strategy (int, optional): 
            binned (bool, optional): Whether the data is binned. Defaults to True

        Returns:

        Raises:
            ValueError: If pdf_type does not match an existing option.
            ImportError: If libHiggsAnalysisCombinedLimit.so does not exist.
        '''
        var = self.__readFromWorkspace(var_name, 'var')
        if pdf_type=='gaussian':
            mu = ROOT.RooRealVar('gauss_mu', 'guass_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gauss_sigma', 'guass_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            pdf = ROOT.RooGaussian(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, mu, sigma)
        elif pdf_type=='double_gaussian':
            pass
        elif pdf_type=='crystal_ball':
            mu = ROOT.RooRealVar('CB_mu', 'CB_mu', self._SR_low, self._SR_high)
            sigmaL = ROOT.RooRealVar('CB_sigmaL', 'CB_sigmaL', 3., 40.)
            sigmaR = ROOT.RooRealVar('CB_sigmaR', 'CB_sigmaR', 3., 20.)
            alphaL = ROOT.RooRealVar('CB_alphaL', 'CB_alphaL', .1, 7.)
            alphaR = ROOT.RooRealVar('CB_alphaR', 'CB_alphaR', .1, 7.)
            nL = ROOT.RooRealVar('CB_nL', 'CB_nL', .005, 6.)
            nR = ROOT.RooRealVar('CB_nR', 'CB_nR', .005, 6.)
            pdf = ROOT.RooCrystalBall(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}',
                                        var, mu, sigmaL, sigmaR, alphaL, alphaR, nL, nR)
        elif pdf_type=='single_double_sided_crystal_ball':
            mu = ROOT.RooRealVar('SDSCB_mu', 'SDSCB_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('SDSCB_sigma', 'SDSCB_sigma', 1., 30.)
            alpha = ROOT.RooRealVar('SDSCB_alpha', 'SDSCB_alpha', .1, 10.)
            n = ROOT.RooRealVar('SDSCB_n', 'SDSCB_n', .1, 10.)
            pdf = ROOT.RooCrystalBall(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}',
                                        var, mu, sigma, alpha, n)
        elif pdf_type=='double_sided_crystal_ball':
            print('{}kytools: Loading libHiggsAnalysisCombinedLimit.so.{}'.format('\033[0;33m', '\033[0m'))
            load = ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit.so')
            if load < 0: raise ImportError('Please compile the Combine package in the README.md instruction.')
            mu = ROOT.RooRealVar('DSCB_mu', 'DSCB_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('DSCB_sigma', 'DSCB_sigma', 1., 30.)
            alphaL = ROOT.RooRealVar('DSCB_alphaL', 'DSCB_alphaL', .1, 10.)
            alphaR = ROOT.RooRealVar('DSCB_alphaR', 'DSCB_alphaR', .1, 10.)
            nL = ROOT.RooRealVar('DSCB_nL', 'DSCB_nL', .1, 10.)
            nR = ROOT.RooRealVar('DSCB_nR', 'DSCB_nR', .1, 10.)
            pdf = ROOT.RooDoubleCBFast(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}',
                                        var, mu, sigma, alphaL, alphaR, nL, nR)
        elif pdf_type=='gaussian_plus_single_double_sided_crystal_ball':
            mu = ROOT.RooRealVar('gauss_mu', 'guass_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gauss_sigma', 'guass_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gauss_pdf = ROOT.RooGaussian(f'{SAMP}_gauss', f'{SAMP}_gauss', var, mu, sigma)
            mu = ROOT.RooRealVar('SDSCB_mu', 'SDSCB_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('SDSCB_sigma', 'SDSCB_sigma', 1., 30.)
            alpha = ROOT.RooRealVar('SDSCB_alpha', 'SDSCB_alpha', .1, 10.)
            n = ROOT.RooRealVar('SDSCB_n', 'SDSCB_n', .1, 10.)
            sdscb_pdf = ROOT.RooCrystalBall(f'{SAMP}_SDSCB', f'{SAMP}_SDSCB',
                                            var, mu, sigma, alpha, n)
            f_sdscb = ROOT.RooRealVar('f_sdscb', 'f_sdscb', 0., 1.) # fraction of the SDSCB
            pdf = ROOT.RooAddPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}',
                                 ROOT.RooArgList(sdscb_pdf, gauss_pdf),
                                 ROOT.RooArgList(f_sdscb),
                                 True)
        elif pdf_type=='exponential':
            exp1_pow = ROOT.RooRealVar('exp1_pow', 'exp1_pow', -0.1, -10., 0.)
            pdf = ROOT.RooExponential(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, exp1_pow)
        elif pdf_type=='double_exponential':
            exp1_pow = ROOT.RooRealVar('exp1_pow', 'exp1_pow', -0.1, -10., 0.)
            exp2_pow = ROOT.RooRealVar('exp2_pow', 'exp2_pow', -0.01, -10., 0.)
            exp1_pdf = ROOT.RooExponential(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, exp1_pow)
            exp2_pdf = ROOT.RooExponential(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, exp2_pow)
            c_exp1 = ROOT.RooRealVar('c_exp1', 'c_exp1', 0., 1.)
            c_exp2 = ROOT.RooRealVar('c_exp2', 'c_exp2', 0., 1.)
            pdf = ROOT.RooAddPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}',
                                 ROOT.RooArgList(exp1_pdf, exp2_pdf),
                                 ROOT.RooArgList(c_exp1, c_exp2))
        elif pdf_type=='gaussian_X_exponential': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            expX_pow = ROOT.RooRealVar('expX_pow', 'expX_pow', -0.1, -10., 0.)
            expX_pdf = ROOT.RooExponential(f'{SAMP}_expX', f'{SAMP}_expX', var, expX_pow)
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, expX_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_step_exponential': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            expX_pow = ROOT.RooRealVar('expX_pow', 'expX_pow', -0.1, -10., 0.)
            expX_pdf = ROOT.RooExponential(f'{SAMP}_expX', f'{SAMP}_expX', var, expX_pow)
            turnon = ROOT.RooRealVar('turnon', 'turnon', self._SR_low, self._SR_high)
            step_pdf = ROOT.RooStats.Heaviside(f'{SAMP}_step', f'{SAMP}_step', var, turnon)
            step_expX_pdf = ROOT.RooProdPdf(f'{SAMP}_step_expX', f'{SAMP}_step_expX', ROOT.RooArgList(step_pdf, expX_pdf))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, step_expX_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_double_exponential':
            exp1_pow = ROOT.RooRealVar('exp1_pow', 'exp1_pow', -0.1, -10., 0.)
            exp2_pow = ROOT.RooRealVar('exp2_pow', 'exp2_pow', -0.01, -10., 0.)
            exp1_pdf = ROOT.RooExponential(f'{SAMP}_exp1', f'{SAMP}_{pdf_type}', var, exp1_pow)
            exp2_pdf = ROOT.RooExponential(f'{SAMP}_exp2', f'{SAMP}_{pdf_type}', var, exp2_pow)
            c_exp1 = ROOT.RooRealVar('c_exp1', 'c_exp1', 0., 1.)
            c_exp2 = ROOT.RooRealVar('c_exp2', 'c_exp2', 0., 1.)
            double_expX_pdf = ROOT.RooAddPdf(f'{SAMP}_double_expX', f'{SAMP}_double_expX',
                                             ROOT.RooArgList(exp1_pdf, exp2_pdf),
                                             ROOT.RooArgList(c_exp1, c_exp2),
                                             False)
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, double_expX_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_step_double_exponential':
            exp1_pow = ROOT.RooRealVar('exp1_pow', 'exp1_pow', -0.1, -10., 0.)
            exp2_pow = ROOT.RooRealVar('exp2_pow', 'exp2_pow', -0.01, -10., 0.)
            exp1_pdf = ROOT.RooExponential(f'{SAMP}_exp1', f'{SAMP}_{pdf_type}', var, exp1_pow)
            exp2_pdf = ROOT.RooExponential(f'{SAMP}_exp2', f'{SAMP}_{pdf_type}', var, exp2_pow)
            c_exp1 = ROOT.RooRealVar('c_exp1', 'c_exp1', 0., 1.)
            c_exp2 = ROOT.RooRealVar('c_exp2', 'c_exp2', 0., 1.)
            double_expX_pdf = ROOT.RooAddPdf(f'{SAMP}_double_expX', f'{SAMP}_double_expX',
                                             ROOT.RooArgList(exp1_pdf, exp2_pdf),
                                             ROOT.RooArgList(c_exp1, c_exp2),
                                             False)
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            turnon = ROOT.RooRealVar('turnon', 'turnon', self._SR_low, self._SR_high)
            step_pdf = ROOT.RooStats.Heaviside(f'{SAMP}_step', f'{SAMP}_step', var, turnon)
            step_double_expX_pdf = ROOT.RooProdPdf(f'{SAMP}_step_double_expX', f'{SAMP}_step_double_expX', ROOT.RooArgList(step_pdf, double_expX_pdf))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, step_double_expX_pdf, gaussX_pdf)
        elif pdf_type=='bernstein_0th_order':
            bern0_c0 = ROOT.RooRealVar('bern0_c0', 'bern0_c0', 0.5, 0., 1.)
            pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                    ROOT.RooArgList(bern0_c0))
        elif pdf_type=='bernstein_1st_order':
            bern1_c0 = ROOT.RooRealVar('bern1_c0', 'bern1_c0', 0.5, 0., 1.)
            bern1_c1 = ROOT.RooRealVar('bern1_c1', 'bern1_c1', 0.1, 0., 1.)
            pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                    ROOT.RooArgList(bern1_c0, bern1_c1))
        elif pdf_type=='bernstein_2nd_order':
            bern2_c0 = ROOT.RooRealVar('bern2_c0', 'bern2_c0', 0.5, 0., 1.)
            bern2_c1 = ROOT.RooRealVar('bern2_c1', 'bern2_c1', 0.1, 0., 1.)
            bern2_c2 = ROOT.RooRealVar('bern2_c2', 'bern2_c2', 0.1, 0., 1.)
            pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                    ROOT.RooArgList(bern2_c0, bern2_c1, bern2_c2))
        elif pdf_type=='bernstein_3rd_order':
            bern3_c0 = ROOT.RooRealVar('bern3_c0', 'bern3_c0', 0.5, 0., 1.)
            bern3_c1 = ROOT.RooRealVar('bern3_c1', 'bern3_c1', 0.1, 0., 1.)
            bern3_c2 = ROOT.RooRealVar('bern3_c2', 'bern3_c2', 0.1, 0., 1.)
            bern3_c3 = ROOT.RooRealVar('bern3_c3', 'bern3_c3', 0.1, 0., 1.)
            pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                    ROOT.RooArgList(bern3_c0, bern3_c1, bern3_c2, bern3_c3))
        elif pdf_type=='bernstein_4th_order':
            bern4_c0 = ROOT.RooRealVar('bern4_c0', 'bern4_c0', 0.5, 0., 1.)
            bern4_c1 = ROOT.RooRealVar('bern4_c1', 'bern4_c1', 0.1, 0., 1.)
            bern4_c2 = ROOT.RooRealVar('bern4_c2', 'bern4_c2', 0.1, 0., 1.)
            bern4_c3 = ROOT.RooRealVar('bern4_c3', 'bern4_c3', 0.1, 0., 1.)
            bern4_c4 = ROOT.RooRealVar('bern4_c4', 'bern4_c4', 0.1, 0., 1.)
            pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                    ROOT.RooArgList(bern4_c0, bern4_c1, bern4_c2, bern4_c3, bern4_c4))
        elif pdf_type=='bernstein_5th_order':
            bern5_c0 = ROOT.RooRealVar('bern5_c0', 'bern5_c0', 0.5, 0., 1.)
            bern5_c1 = ROOT.RooRealVar('bern5_c1', 'bern5_c1', 0.1, 0., 1.)
            bern5_c2 = ROOT.RooRealVar('bern5_c2', 'bern5_c2', 0.1, 0., 1.)
            bern5_c3 = ROOT.RooRealVar('bern5_c3', 'bern5_c3', 0.1, 0., 1.)
            bern5_c4 = ROOT.RooRealVar('bern5_c4', 'bern5_c4', 0.1, 0., 1.)
            bern5_c5 = ROOT.RooRealVar('bern5_c5', 'bern5_c5', 0.1, 0., 1.)
            pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                    ROOT.RooArgList(bern5_c0, bern5_c1, bern5_c2, bern5_c3, bern5_c4, bern5_c5))
        elif pdf_type=='gaussian_X_bernstein_1st_order': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            bern1_c0 = ROOT.RooRealVar('bern1_c0', 'bern1_c0', 0.5, 0., 1.)
            bern1_c1 = ROOT.RooRealVar('bern1_c1', 'bern1_c1', 0.1, 0., 1.)
            bern1X_pdf = ROOT.RooBernstein(f'{SAMP}_bern1X', f'{SAMP}_bern1X', var,
                                            ROOT.RooArgList(bern1_c0, bern1_c1))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, bern1X_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_bernstein_2nd_order': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            bern2_c0 = ROOT.RooRealVar('bern2_c0', 'bern2_c0', 0.5, 0., 1.)
            bern2_c1 = ROOT.RooRealVar('bern2_c1', 'bern2_c1', 0.1, 0., 1.)
            bern2_c2 = ROOT.RooRealVar('bern2_c2', 'bern2_c2', 0.1, 0., 1.)
            bern2X_pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                            ROOT.RooArgList(bern2_c0, bern2_c1, bern2_c2))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, bern2X_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_bernstein_3rd_order': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            bern3_c0 = ROOT.RooRealVar('bern3_c0', 'bern3_c0', 0.5, 0., 1.)
            bern3_c1 = ROOT.RooRealVar('bern3_c1', 'bern3_c1', 0.1, 0., 1.)
            bern3_c2 = ROOT.RooRealVar('bern3_c2', 'bern3_c2', 0.1, 0., 1.)
            bern3_c3 = ROOT.RooRealVar('bern3_c3', 'bern3_c3', 0.1, 0., 1.)
            bern3X_pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                            ROOT.RooArgList(bern3_c0, bern3_c1, bern3_c2, bern3_c3))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, bern3X_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_step_double_exponential':
            exp1_pow = ROOT.RooRealVar('exp1_pow', 'exp1_pow', -0.1, -10., 0.)
            exp2_pow = ROOT.RooRealVar('exp2_pow', 'exp2_pow', -0.01, -10., 0.)
            exp1_pdf = ROOT.RooExponential(f'{SAMP}_exp1', f'{SAMP}_{pdf_type}', var, exp1_pow)
            exp2_pdf = ROOT.RooExponential(f'{SAMP}_exp2', f'{SAMP}_{pdf_type}', var, exp2_pow)
            c_exp1 = ROOT.RooRealVar('c_exp1', 'c_exp1', 0., 1.)
            c_exp2 = ROOT.RooRealVar('c_exp2', 'c_exp2', 0., 1.)
            double_expX_pdf = ROOT.RooAddPdf(f'{SAMP}_double_expX', f'{SAMP}_double_expX',
                                             ROOT.RooArgList(exp1_pdf, exp2_pdf),
                                             ROOT.RooArgList(c_exp1, c_exp2),
                                             False)
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            turnon = ROOT.RooRealVar('turnon', 'turnon', self._SR_low, self._SR_high)
            step_pdf = ROOT.RooStats.Heaviside(f'{SAMP}_step', f'{SAMP}_step', var, turnon)
            step_double_expX_pdf = ROOT.RooProdPdf(f'{SAMP}_step_double_expX', f'{SAMP}_step_double_expX',
                                                   ROOT.RooArgList(step_pdf, double_expX_pdf))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, step_double_expX_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_step_bernstein_3rd_order': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            bern3_c0 = ROOT.RooRealVar('bern3_c0', 'bern3_c0', 0.5, 0., 1.)
            bern3_c1 = ROOT.RooRealVar('bern3_c1', 'bern3_c1', 0.1, 0., 1.)
            bern3_c2 = ROOT.RooRealVar('bern3_c2', 'bern3_c2', 0.1, 0., 1.)
            bern3_c3 = ROOT.RooRealVar('bern3_c3', 'bern3_c3', 0.1, 0., 1.)
            bern3X_pdf = ROOT.RooBernstein(f'{SAMP}_bern3X', f'{SAMP}_bern3X', var,
                                            ROOT.RooArgList(bern3_c0, bern3_c1, bern3_c2, bern3_c3))
            turnon = ROOT.RooRealVar('turnon', 'turnon', self._SR_low, self._SR_high)
            step_pdf = ROOT.RooStats.Heaviside(f'{SAMP}_step', f'{SAMP}_step', var, turnon)
            step_bern3X_pdf = ROOT.RooProdPdf(f'{SAMP}_step_bern3X', f'{SAMP}_step_bern3X',
                                              ROOT.RooArgList(step_pdf, bern3X_pdf))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, step_bern3X_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_bernstein_4th_order': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            bern4_c0 = ROOT.RooRealVar('bern4_c0', 'bern4_c0', 0.5, 0., 1.)
            bern4_c1 = ROOT.RooRealVar('bern4_c1', 'bern4_c1', 0.1, 0., 1.)
            bern4_c2 = ROOT.RooRealVar('bern4_c2', 'bern4_c2', 0.1, 0., 1.)
            bern4_c3 = ROOT.RooRealVar('bern4_c3', 'bern4_c3', 0.1, 0., 1.)
            bern4_c4 = ROOT.RooRealVar('bern4_c4', 'bern4_c4', 0.1, 0., 1.)
            bern4X_pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                            ROOT.RooArgList(bern4_c0, bern4_c1, bern4_c2, bern4_c3, bern4_c4))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, bern4X_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_step_bernstein_4th_order': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            bern4_c0 = ROOT.RooRealVar('bern4_c0', 'bern4_c0', 0.5, 0., 1.)
            bern4_c1 = ROOT.RooRealVar('bern4_c1', 'bern4_c1', 0.1, 0., 1.)
            bern4_c2 = ROOT.RooRealVar('bern4_c2', 'bern4_c2', 0.1, 0., 1.)
            bern4_c3 = ROOT.RooRealVar('bern4_c3', 'bern4_c3', 0.1, 0., 1.)
            bern4_c4 = ROOT.RooRealVar('bern4_c4', 'bern4_c4', 0.1, 0., 1.)
            bern4X_pdf = ROOT.RooBernstein(f'{SAMP}_bern4X', f'{SAMP}_bern4X', var,
                                            ROOT.RooArgList(bern4_c0, bern4_c1, bern4_c2, bern4_c3, bern4_c4))
            turnon = ROOT.RooRealVar('turnon', 'turnon', self._SR_low, self._SR_high)
            step_pdf = ROOT.RooStats.Heaviside(f'{SAMP}_step', f'{SAMP}_step', var, turnon)
            step_bern4X_pdf = ROOT.RooProdPdf(f'{SAMP}_step_bern4X', f'{SAMP}_step_bern4X',
                                              ROOT.RooArgList(step_pdf, bern4X_pdf))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, step_bern4X_pdf, gaussX_pdf)
        elif pdf_type=='gaussian_X_bernstein_5th_order': # convolution
            mu = ROOT.RooRealVar('gaussX_mu', 'gaussX_mu', self._SR_low, self._SR_high)
            sigma = ROOT.RooRealVar('gaussX_sigma', 'gaussX_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
            gaussX_pdf = ROOT.RooGaussian(f'{SAMP}_gaussX', f'{SAMP}_gaussX', var, mu, sigma)
            bern5_c0 = ROOT.RooRealVar('bern5_c0', 'bern5_c0', 0.5, 0., 1.)
            bern5_c1 = ROOT.RooRealVar('bern5_c1', 'bern5_c1', 0.1, 0., 1.)
            bern5_c2 = ROOT.RooRealVar('bern5_c2', 'bern5_c2', 0.1, 0., 1.)
            bern5_c3 = ROOT.RooRealVar('bern5_c3', 'bern5_c3', 0.1, 0., 1.)
            bern5_c4 = ROOT.RooRealVar('bern5_c4', 'bern5_c4', 0.1, 0., 1.)
            bern5_c5 = ROOT.RooRealVar('bern5_c5', 'bern5_c5', 0.1, 0., 1.)
            bern5X_pdf = ROOT.RooBernstein(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var,
                                            ROOT.RooArgList(bern5_c0, bern5_c1, bern5_c2, bern5_c3, bern5_c4, bern5_c5))
            pdf = ROOT.RooFFTConvPdf(f'{SAMP}_{pdf_type}', f'{SAMP}_{pdf_type}', var, gaussX_pdf, bern5X_pdf)
        # TODO: add laurent
        else:
            raise ValueError('Invalid pdf_type.')
        params = pdf.getParameters(0)
        status = -1
        data_name = f'{SAMP}_{self.YEAR}_{self.CAT}_data{"hist" if binned else "set"}'
        data = self.__readFromWorkspace(data_name, 'data')
        print ('{}Performing likelihood fit of {} to {}.{}'.format('\033[1;36m', pdf.GetTitle(), data.GetTitle(), '\033[0m'))
        for ntries in range(1, max_tries+1):
            print ('{}kytools: Fit trial #{}.{}'.format('\033[0;36m', ntries, '\033[0m'))
            fit_result = pdf.fitTo(data,
                                   ROOT.RooFit.Save(True),
                                   ROOT.RooFit.Minimizer('Minuit2', 'minimize'),
                                   ROOT.RooFit.Strategy(strategy),
                                   ROOT.RooFit.PrintLevel(-1),
                                   ROOT.RooFit.Warnings(False),
                                   ROOT.RooFit.PrintEvalErrors(-1),
                                   ROOT.RooFit.SumW2Error(True))
            status = fit_result.status()
            if (status != 0):                                                                                                                              
                params.assignValueOnly(fit_result.randomizePars())                                                                                                                              
                ntries += 1
            else:
                break
        print ('{}kytools: Likelihood fit has exited with status {}.{}'.format('\033[0;36m', status, '\033[0m'))
        if status==0:
            print ('{}         Likelihood fit has converged.{}'.format('\033[0;36m', '\033[0m'))
        elif status==1:
            print ('{}         Covariance was made positive definite.{}'.format('\033[0;36m', '\033[0m'))
        elif status==2:
            print ('{}         Hessian is invalid.{}'.format('\033[0;36m', '\033[0m'))
        elif status==3:
            print ('{}         EDM is above max.{}'.format('\033[0;36m', '\033[0m'))
        elif status==4:
            print ('{}         Reached call limit.{}'.format('\033[0;36m', '\033[0m'))
        elif status==5:
            print ('{}         Please investigate.{}'.format('\033[0;36m', '\033[0m'))
        else:
            print ('{}         DISASTER!{}'.format('\033[0;36m', '\033[0m'))
        self.__importToWorkspace(pdf)
        return

    def configPlot(self, SAMP: str, preset: str, pdf_type='', binned=True):
        '''Configure plot for the given data or pdf.

        Args:
            SAMP (str): Name of the sample.
            preset (str): The pre-set plot configuration to apply. The options are
                'mc_sig', 'mc_bkg1', 'mc_bkg2', 'mc_bkg3', or 'data'.
            pdf_type (str, optional): PDF type. Defaults to ''.
                Argument must be provided if preset is pdf.
            binned (bool, optional): Whether the sample data is binned. Defaults to True.
                Irrelevant if preset is not 'data'.
        
        Returns:
            (None)
        
        Raises:
            ValueError (str): If preset is not one of the following:
                'mc_sig', 'mc_bkg1', 'mc_bkg2', 'mc_bkg3', or 'data'.
        '''
        if preset=='sig_pdf':
            self._plot_configs['sig_pdf'] = {
                'draw_type': 'line',
                'obj_name': f'{SAMP}_{pdf_type}',
                'obj_title': pdf_type,
                'obj_type': 'pdf',
                'SAMP': SAMP,
                'LineColor': ROOT.kRed,
                'LineStyle': ROOT.kSolid,
                'LineWidth': 2,
                'MarkerColor': ROOT.kWhite,
                'MarkerSize': 0,
                'MarkerStyle': ROOT.kDot,
                'FillColor': ROOT.kWhite,
                'FillStyle': 0, # hollow
            }
        elif preset=='bkg_pdf':
            self._plot_configs['sig_pdf'] = {
                'draw_type': 'line',
                'obj_name': f'{SAMP}_{pdf_type}',
                'obj_title': pdf_type,
                'obj_type': 'pdf',
                'SAMP': SAMP,
                'LineColor': ROOT.kBlue,
                'LineStyle': ROOT.kSolid,
                'LineWidth': 2,
                'MarkerColor': ROOT.kWhite,
                'MarkerSize': 0,
                'MarkerStyle': ROOT.kDot,
                'FillColor': ROOT.kWhite,
                'FillStyle': 0, # hollow
            }
        elif preset=='mc_bkg1':
            self._plot_configs['mc_bkg1'] = {
                'draw_type': 'solid',
                'obj_name': f'{SAMP}_{self.YEAR}_{self.CAT}_data{"hist" if binned else "set"}',
                'obj_title': f'{SAMP}_{self.YEAR}_{self.CAT}',
                'obj_type': 'data',
                'SAMP': SAMP,
                'LineColor': ROOT.kBlack,
                'LineStyle': ROOT.kSolid,
                'LineWidth': 0,
                'MarkerColor': ROOT.kWhite,
                'MarkerSize': 0,
                'MarkerStyle': ROOT.kDot,
                'FillColor': ROOT.kRed - 4,
                'FillStyle': 1001 # solid
            }
        elif preset=='mc_bkg2':
            self._plot_configs['mc_bkg1'] = {
                'draw_type': 'solid',
                'obj_name': f'{SAMP}_{self.YEAR}_{self.CAT}_data{"hist" if binned else "set"}',
                'obj_title': f'{SAMP}_{self.YEAR}_{self.CAT}',
                'obj_type': 'data',
                'SAMP': SAMP,
                'LineColor': ROOT.kBlack,
                'LineStyle': ROOT.kSolid,
                'LineWidth': 0,
                'MarkerColor': ROOT.kWhite,
                'MarkerSize': 0,
                'MarkerStyle': ROOT.kDot,
                'FillColor': ROOT.kOrange + 1,
                'FillStyle': 1001
            }
        elif preset=='mc_bkg3':
            pass
        elif preset=='data':
            self._plot_configs['data'] = {
                'draw_type': 'marker',
                'obj_name': f'{SAMP}_{self.YEAR}_{self.CAT}_data{"hist" if binned else "set"}',
                'obj_title': f'{SAMP}_{self.YEAR}_{self.CAT}',
                'obj_type': 'data',
                'SAMP': SAMP,
                'LineColor': ROOT.kBlack,
                'LineStyle': ROOT.kSolid,
                'LineWidth': 1,
                'MarkerColor': ROOT.kBlack,
                'MarkerSize': 0.8,
                'MarkerStyle': ROOT.kFullCircle,
                'FillColor': ROOT.kWhite,
                'FillStyle': 0, # hollow
            }
        else:
            raise ValueError('Invalid option for preset.')
        return
    
    def clearConfigPlot(self):
        '''Clear the plot configs.

        Args:
            (None)

        Returns:
            (None)
        '''
        self._plot_configs = {}
        return
    
    def makePlot(self, var_name: str, plot_name: str, plot_title: str, show_SR=True, resid=True, plot_width=1200, plot_height=800):
        '''Make plot.

        Args:
            var_name (str): Name of the variable to plot.
            plot_name (str): Plot name to use for the saved plot.
            plot_title (str): Plot title to display
            show_SR (bool, optional): Display SR as vertical lines on the plot.
                Defaults to True.
            resid (bool, optional): Show residuals.
                Defaults to True.
            plot_width (int, optional): Plot width. Defaults to 600.
            plot_height (int, optional): Plot height. Defaults to 400.

        Returns:
            (None)

        Raises:
            RuntimeError: IF the plot configurations are not specified.
        '''
        if not self._plot_configs:
            raise RuntimeError('Plot configurations are not specified. Call configPlot first.')
        var = self.__readFromWorkspace(var_name, 'var')
        xframe = var.frame(Title=plot_title)
        self._plot_configs = dict(sorted(self._plot_configs.items())) # data, mc_bkg, mc_sig
        c0 = ROOT.TCanvas('c0', plot_title, plot_width, plot_height) # Need to create canvas before legend
        legend = ROOT.TLegend(0.64, 0.84, 0.97, 0.95)
        legend.SetTextSize(.024)
        for config in self._plot_configs.values():
            obj = self.__readFromWorkspace(config['obj_name'], config['obj_type'])
            if not obj:
                print (f'{config["obj_name"]} is a null pointer.')
                continue
            if config['draw_type']=='line':
                obj.plotOn(xframe,
                            Name=config['obj_name'],
                            LineColor=config['LineColor'],
                            LineStyle=config['LineStyle'],
                            LineWidth=config['LineWidth'],
                            DrawOption='L')
                legend.AddEntry(xframe.findObject(config['obj_name']), config['obj_title'], 'LC')
            elif config['draw_type']=='solid':
                obj.plotOn(xframe,
                            Name=config['obj_name'],
                            FillColor=config['FillColor'],
                            FillStyle=config['FillStyle'],
                            LineWidth=config['LineWidth'],
                            MarkerSize=config['MarkerSize'],
                            DrawOption='B')
                legend.AddEntry(xframe.findObject(config['obj_name']), config['obj_title'], 'F')
            elif config['draw_type']=='marker':
                obj.plotOn(xframe,
                            Name=config['obj_name'],
                            MarkerColor=config['MarkerColor'],
                            MarkerSize=config['MarkerSize'],
                            MarkerStyle=config['MarkerStyle'])
                legend.AddEntry(xframe.findObject(config['obj_name']), config['obj_title'], 'EP')
            else:
                pass
        xframe.Draw()
        if show_SR:
            line_SR_low = ROOT.TLine(self._SR_low, xframe.GetMinimum(), self._SR_low, xframe.GetMaximum())
            line_SR_high = ROOT.TLine(self._SR_high, xframe.GetMinimum(), self._SR_high, xframe.GetMaximum())
            line_SR_low.SetLineWidth(1)
            line_SR_low.SetLineStyle(ROOT.kDashed)
            line_SR_high.SetLineWidth(1)
            line_SR_high.SetLineStyle(ROOT.kDashed)
            line_SR_low.Draw()
            line_SR_high.Draw()
        legend.Draw()
        c0.SaveAs(os.path.join(self._plotsavedir, f'{plot_name}.png'))
        c0.Close()
        if resid:
            c1 = ROOT.TCanvas('c1', plot_title, plot_width, int(plot_height*1.25))
            pad1 = ROOT.TPad('pad1', 'pad1', 0, .2, 1., 1.)
            pad2 = ROOT.TPad('pad1', 'pad1', 0, 0., 1., .2)
            pad1.Draw()
            pad2.Draw()
            pad1.cd()
            xframe.Draw()
            legend.Draw()
            pad2.cd()
            hresid = xframe.residHist()
            rframe = var.frame(Title='residuals')
            rframe.addPlotable(hresid, 'P')
            rframe.Draw()
            c1.Update()
            c1.SaveAs(os.path.join(self._plotsavedir, f'{plot_name}_resid.png'))
            c1.Close()
        return