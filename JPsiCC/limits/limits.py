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
        # SAMP (str): Either one of the following options.
        #     'DATA_BKG', 'MC_BKG', 'MC_BKG1','MC_BKG2', 'MC_BKG3', 'MC_BKG4', 'MC_SIG'
        YEAR (int): Year of data-taking.
        VERS (str): Version of the files.
        CAT (str): Category of the analysis.
        CMSSW (str): Version of the CMSSW.
        weights (bool): Whether weights were used.

    Raises:
        # ValueError: If the string provided for SAMP is not among the options.
        TypeError: If the value provided for YEAR is not an integer.
        TypeError: If the value provided for VERS is not a string.
        TypeError: If the value provided for CAT is not a string.
    '''
    def __init__(self, YEAR: int, VERS: str, CAT: str, CMSSW: str, weights=True):
        # match SAMP:
        #     case 'DATA_BKG': self._DATA, self._MODE = True, 'BKG'
        #     case 'MC_BKG': self._DATA, self._MODE = False, 'BKG'
        #     case 'MC_BKG1': self._DATA, self._MODE = False, 'BKG'
        #     case 'MC_BKG2': self._DATA, self._MODE = False, 'BKG'
        #     case 'MC_BKG3': self._DATA, self._MODE = False, 'BKG'
        #     case 'MC_BKG4': self._DATA, self._MODE = False, 'BKG'
        #     case 'MC_SIG': self._DATA, self._MODE = False, 'SIG'
        #     case _: raise ValueError(f'SAMP={SAMP} is not a valid option.')
        if not type(YEAR) is int: raise TypeError(f'YEAR must be an integer.')
        if not type(VERS) is str: raise TypeError(f'VERS must be a string.')
        if not type(CAT) is str: raise TypeError(f'CAT must be a string.')
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
        self._wfilename = f'workspace_{self._sfx}.root' # Name of the ROOT workspace file
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
        match obj_type:
            case 'var':
                obj = w.var(obj_name)
            case 'pdf':
                obj = w.pdf(obj_name)
            case 'data':
                obj = w.data(obj_name)
            case _:
                raise ValueError('Invalid obj_type.')
        return obj

    def addVar(self, var_name: str, var_title: str, value: float, var_min: float, var_max: float, unit=''):
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

    def addData(self, SAMP, filename: str, treename: str, col_name: str, var_name: str, binned=True):
        '''Add data to the workspace.

        If binned, it will add a RooDataHist. Otherwise, it will add a RooDataSet.
        For reference, see https://root.cern/doc/v632/rf408__RDataFrameToRooFit_8py_source.html.

        Args:
            SAMP (str): Name of the sample.
            filename (str): Name of the file containing the data.
            treename (str): Name of the ROOT.TTree.
            col_name (str): Name of the column of the data.
            var_name (str): Name of the variable as appears in the workspace.
            binned (bool, optional): Whether the data is binned. Defaults to True.
        
        Returns:
            (None)
        '''
        rdf = ROOT.RDataFrame(treename, filename)
        var = self.__readFromWorkspace(var_name, 'var')
        data_name = f'{SAMP}_{self.YEAR}_{self.CAT}_data{"hist" if binned else "set"}'
        data_title = data_name
        if binned:
            data_maker = ROOT.RooDataHistHelper(data_name, data_title, ROOT.RooArgSet(var))
        else:
            data_maker = ROOT.RooDataSetHelper(data_name, data_title, ROOT.RooArgSet(var))
        data_result = rdf.Book(ROOT.std.move(data_maker), [col_name])
        data = data_result.GetValue()
        self.__importToWorkspace(data)
        return

    def defineRegions(self, SR_low: float, SR_high: float, CR_low: float, CR_high: float):
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
        if not all(isinstance(arg, float) for arg in list(locals().values())[1:]):
            raise TypeError('Input arguments must be float.')
        elif not (SR_low < SR_high):
            raise ValueError('SR_low must be smaller than SR_high.')
        elif not (CR_low < CR_high):
            raise ValueError('CR_low must be smaller than CR_high.')
        elif not (CR_low < SR_low and SR_high < CR_high):
            raise ValueError('Signal region must lie outside the control region.')
        self._SR_low, self._SR_high, self._CR_low, self._CR_high = SR_low, SR_high, CR_low, CR_high
        return

    def addPDF(self, SAMP: str, pdf_type: str, var_name: str, max_tries=10, strategy=1, binned=True):
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
        '''
        var = self.__readFromWorkspace(var_name, 'var')
        match pdf_type:
            case 'gaussian':
                mu = ROOT.RooRealVar('gaussian_mu', 'guassian_mu', self._SR_low, self._SR_high)
                sigma = ROOT.RooRealVar('gaussian_sigma', 'guassian_sigma', 0.1, (self._SR_high-self._SR_low)*10.)
                pdf = ROOT.RooGaussian('gaussian', 'gaussian', var, mu, sigma)
            case 'double_gaussian':
                pass
            case 'crystal_ball':
                pass
            case _:
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
                                   ROOT.RooFit.PrintEvalErrors(-1))
            status = fit_result.status()
            if (status != 0):                                                                                                                              
                params.assignValueOnly(fit_result.randomizePars())                                                                                                                              
                ntries += 1
            else:
                break
        print ('{}kytools: Likelihood fit has exited with status {}.{}'.format('\033[0;36m', status, '\033[0m'))
        match status:
            case 0:
                print ('{}         Likelihood fit has converged.{}'.format('\033[0;36m', '\033[0m'))
            case 1:
                print ('{}         Covariance was made positive definite.{}'.format('\033[0;36m', '\033[0m'))
            case 2:
                print ('{}         Hessian is invalid.{}'.format('\033[0;36m', '\033[0m'))
            case 3:
                print ('{}         EDM is above max.{}'.format('\033[0;36m', '\033[0m'))
            case 4:
                print ('{}         Reached call limit.{}'.format('\033[0;36m', '\033[0m'))
            case 5:
                print ('{}         Please investigate.{}'.format('\033[0;36m', '\033[0m'))
            case _:
                print ('{}         DISASTER!{}'.format('\033[0;36m', '\033[0m'))
        self.__importToWorkspace(pdf)
        return

    def configPlot(self, preset: str, SAMP: str, pdf_type='', binned=True):
        '''Configure plot for the given data or pdf.

        Args:
            preset (str): The pre-set plot configuration to apply. The options are
                'mc_sig', 'mc_bkg1', 'mc_bkg2', 'mc_bkg3', or 'data'.
            SAMP (str): Name of the sample.
            pdf_type (str, optional): PDF type. Defaults to ''.
                Argument must be provided if preset is not 'data'.
            binned (bool, optional): Whether the sample data is binned. Defaults to True.
                Irrelevant if preset is not 'data'.
        
        Returns:
            (None)
        
        Raises:
            ValueError (str): If preset is not one of the following:
                'mc_sig', 'mc_bkg1', 'mc_bkg2', 'mc_bkg3', or 'data'.
        '''
        match preset:
            case 'mc_sig':
                self._plot_configs['mc_sig'] = {
                    'obj_name': pdf_type,
                    'obj_type': 'pdf',
                    'SAMP': SAMP,
                    'LineColor': ROOT.kRed,
                    'LineStyle': ROOT.kSolid,
                    'LineWidth': 1,
                    'MarkerColor': ROOT.kWhite,
                    'MarkerSize': 0,
                    'MarkerStyle': ROOT.kDot,
                    'FillColor': ROOT.kWhite,
                    'FillStyle': 0, # hollow
                }
            case 'mc_bkg1':
                pass
            case 'mc_bkg2':
                pass
            case 'mc_bkg3':
                pass
            case 'data':
                self._plot_configs['data'] = {
                    'obj_name': f'{SAMP}_{self.YEAR}_{self.CAT}_data{"hist" if binned else "set"}',
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
            case _:
                raise ValueError('Invalid option for preset.')
        return
    
    def makePlot(self, var_name: str, plot_title: str, plot_width=1200, plot_height=800):
        '''Make plot.

        Args:
            var_name (str): Name of the variable to plot.
            plot_title (str): Plot title.
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
        for config in self._plot_configs.values():
            obj = self.__readFromWorkspace(config['obj_name'], config['obj_type'])
            obj.plotOn(xframe, # TODO: make args different according to preset
                       LineColor=config['LineColor'],
                       LineStyle=config['LineStyle'],
                       LineWidth=config['LineWidth'],
                       MarkerColor=config['MarkerColor'],
                       MarkerSize=config['MarkerSize'],
                       MarkerStyle=config['MarkerStyle'],
                    #    FillColor=config['FillColor'],
                    #    FillStyle=config['FillStyle'],
                       XErrorSize=0)
        c = ROOT.TCanvas('c', plot_title, plot_width, plot_height)
        xframe.Draw()
        c.SaveAs(os.path.join(self._plotsavedir, 'roofitplot.png'))
        return

if __name__=='__main__':
    verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)
    snapshot_dir = '/work/submit/kyoon/CMSSW_13_3_0/src/RareHiggsRun3/JPsiCC/analysis'
    sig_data_name = os.path.join(snapshot_dir, 'snapshot_MC_SIG_2018_GF_v202407_20241024_WEIGHT_no_filter.root')
    rwc = RooWorkspaceCreator(2018, '202407', 'GF', 'CMSSW_13_3_0')
    rwc.addVar(var_name='mH', var_title='m_{#mu#bar{#mu}jj}', value=125, var_min=0., var_max=200., unit='GeV/c^2')
    rwc.defineRegions(SR_low=80., SR_high=140., CR_low=0., CR_high=200.)
    rwc.addData(SAMP='MC_SIG', filename=sig_data_name, treename='Events', col_name='higgs_mass_corr', var_name='mH')
    rwc.addPDF(SAMP='MC_SIG', pdf_type='gaussian', var_name='mH')
    rwc.configPlot(preset='mc_sig', SAMP='MC_SIG', pdf_type='gaussian')
    rwc.configPlot(preset='data', SAMP='MC_SIG')
    rwc.makePlot(var_name='mH', plot_title='mc sig + gaussian')