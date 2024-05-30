"""
rootpdf
=========================================

A module to create, fit, and store PDF models in ROOT.

"""
from ROOT import RooRealVar, RooArgList
from ROOT import RooGaussian, RooBernstein
from ROOT import RooWorkspace, RooDataHist, RooDataSet
from ROOT import RooFit
import os

class Models:
    """Creates PDF models.

    The __init__ method creates two internal lists for data storage, one for variables
    and another for PDFs. This ensures that the variable and PDF objects do not get
    delete during the lifetime of the class instance.
    """
    def __init__(self):
        print(f'{'\033[1;34m'}kytools: You have created an instance of rootpdf.Models.{'\033[0m'}')
        self._varlist = []
        self._pdflist = []

    def __var(self, name, title, value, low, high, unit=''):
        """Internal method for creating a RooRealVar.

        Args:
            name (str): Name of the variable.
            title (str): Display title.
            value (float): Nominal value.
            low (float): Lower limit.
            high (float): Upper limit.
            unit (str, optional): Units of this value.
        """
        var = RooRealVar(name, title, value, low, high, unit)
        self._varlist.append(var)
        print (f'{'\033[0;36m'}kytools: Created a variable with the name, {name}.{'\033[0m'}')
        return var

    def gaussian(self, x, mean, mean_low, mean_high,
                 sigma, sigma_low, sigma_high, suffix):
        """Creates and returns a RooGaussian.

        Args:
            x (RooRealVar): Key variable.
            mean (float): Nominal mean value.
            mean_low (float): Lower limit of the mean.
            mean_high (float): Upper limit of the mean.
            sigma (float): Nominal sigma value.
            sigma_low (float): Lower limit of the mean.
            sigma_high (float): Upper limit of the mean.
            suffix (str): Suffix for the name of the variables.
        """
        mu = self.__var('gauss_mu_{}'.format(suffix), 'gauss_mu', mean, mean_low, mean_high)
        sigma = self.__var('gauss_sigma_{}'.format(suffix), 'gauss_sigma', sigma, sigma_low, sigma_high)
        pdfname = 'gauss'
        pdf = RooGaussian(pdfname, pdfname, x, mu, sigma)
        self._pdflist.append(pdf)
        print (f'{'\033[0;32m'}kytools: Created a Gaussian with the name, {pdfname}.{'\033[0m'}')
        return pdf
    
    def crystalball(self, x, mean, mean_low, mean_high,
                    width=1., width_low=0., width_high=3.,
                    aL=1., aL_low=0., aL_high=5.,
                    aR=1., aR_low=0., aR_high=5.,
                    nL=5., nL_low=2., nL_high=50.,
                    nR=2., nR_low=0., nR_high=5.):
        """Creates and returns a RooCrystalBall

        """
        return
    
    def bernstein(self, x, degree, suffix, coeff_vals):
        """Creates and returns a Bernstein polynomial of degree n.

        Args:
            x (RooRealVar): Key variable.
            degree (int): Degree of the polynomial
            suffix (str): Suffix for the name of the variables.
            coeff_vals (list of tuple(float) or list of list(float)): For each
                Bernstein coefficient, provide the nominal value, lower limit, and
                upper limit as a tuple or list. For example, if I create a Bernstein
                polynomial of 2nd order, I would provide for this argument, e.g.
                    coeffs=[(.2, 0., .5), (.1, 0., 1.), (.1, 0., 1.)]
        """
        coefflist = RooArgList('bern{}_coeffs'.format(degree))
        if not len(coeff_vals) == (degree+1):
            raise Exception('Length of coeff_vals must be equal to the degrees of freedom.')
        for i in range(0, degree+1):
            if not len(coeff_vals[i]) == 3:
                raise Exception('Each item in coeff_vals must be a tuple or list of length 3.')
            coeff = RooRealVar('bern_c{}_{}'.format(i, suffix), 'bern_c{}'.format(i), coeff_vals[i][0], coeff_vals[i][1], coeff_vals[i][2])
            self._varlist.append(coeff)
            coefflist.add(coeff)
        pdfname = 'bern{}_{}'.format(degree, suffix)
        pdf = RooBernstein(pdfname, pdfname, x, coefflist)
        self._pdflist.append(pdf)
        print (f'{'\033[0;32m'}kytools: Created a Bernstein polynomial of degree {degree} (d.o.f. = {degree+1}) with the name, {pdfname}.{'\033[0m'}')
        return pdf

# ==============================================================================

class FitIt:
    """This class provides methods for fitting the model.

    TODO: __init__ description

    Args:
        x (RooRealVar): Key variable.
        wsname (str): Name of the RooWorkspace.

    Attributes: TODO
    """

    def __init__(self, x, ws_name):
        self._x = x
        self._workspace = RooWorkspace(ws_name, 'workspace') # RooWorkspace
        self._data = dict() # Dictionary for RooDataHist or RooDataSet
        self._pdf = dict() # Dictionary for RooAbsPdf

        print(f'{'\033[1;34m'}kytools: You have created an instance of rootpdf.FitIt.{'\033[0m'}')
        print(f'{'\033[0;32m'}kytools: Created a RooWorkspace with the name, {ws_name}.{'\033[0m'}')

    # def blind(self, range_low, range_high):
    #     """TODO

    #     Attributes:
    #         range_low (float): Low range of the blinded region.
    #         range_high (float): High range of the blinded region.
    #     """
    #     self._x.setRange('left', xlow, range_low)
    #     self._x.setRange('right', range_high, xhigh)

    def add_data(self, data, data_key):
        """Adds RooDataHist or RooDataSet data to this class instance.

        Attributes:
            data (RooDataHist or RooDataSet): The data.
            key (str): String key for storing this data in a dictionary.
        """
        self._data[data_key] = data

    def add_pdf(self, pdf, pdf_key):
        """Adds RooAbsPdf to this class instance.

        Attributes:
            pdf (RooAbsPdf): The pdf.
            key (str): String key for storing this pdf in a dictionary.
        """
        self._pdf[pdf_key] = pdf
    
    def fit(self, data_key, pdf_key, strategy=2, max_tries=10):
        """Fit the data to the model.

        Used here is RooFit's `FitTo` method, which is based on likelihood maximization.

        If the fit fails, the method will try re-fitting with randomized parameters
        until max tries are reached.

        Attributes:
            data_key (str): String key of the data to fit.
            pdf_key (str): String key of the PDF model to fit.
            strategy (int, optional): RooFit fitting strategy. See RooFit documentation.
                Defaults to 2.
            max_tries (int, optional): Number of max tries to reach. Defaults to 10.

        Returns:
            TODO
        """
        data = self._data[data_key]
        pdf = self._pdf[pdf_key]
        params = pdf.getParameters(0)
        status = -1
        print (f'{'\033[1;36m'}kytools: Performing likelihood fit of {pdf.GetTitle()} to {data.GetTitle()}.{'\033[0m'}')
        for ntries in range(1, max_tries+1):
            print (f'{'\033[0;36m'}kytools: Fit trial #{ntries}.{'\033[0m'}')
            fit_result = pdf.fitTo(data,
                                   RooFit.Save(True),
                                   RooFit.Minimizer('Minuit2', 'minimize'),
                                   RooFit.Strategy(strategy),
                                   RooFit.PrintLevel(-1),
                                   RooFit.Warnings(False),
                                   RooFit.PrintEvalErrors(-1))
            status = fit_result.status()
            if (status != 0):                                                                                                                              
                params.assignValueOnly(fit_result.randomizePars())                                                                                                                              
                ntries += 1
            else:
                break
        print (f'{'\033[0;36m'}kytools: Likelihood fit has exited with status {status}.{'\033[0m'}')
        match status:
            case 0:
                print (f'{'\033[0;36m'}         Likelihood fit has converged.{'\033[0m'}')
            case 1:
                print (f'{'\033[0;36m'}         Covariance was made positive definite.{'\033[0m'}')
            case 2:
                print (f'{'\033[0;36m'}         Hessian is invalid.{'\033[0m'}')
            case 3:
                print (f'{'\033[0;36m'}         EDM is above max.{'\033[0m'}')
            case 4:
                print (f'{'\033[0;36m'}         Reached call limit.{'\033[0m'}')
            case 5:
                print (f'{'\033[0;36m'}         Please investigate.{'\033[0m'}')
            case _:
                print (f'{'\033[0;36m'}         DISASTER!{'\033[0m'}')

    def add_norm(self):
        pass

    def draw(self):
        pass

    def wsave(self, fname, path='.'):
        """Saves all data and pdf to workspace.

        Attributes:
            fname (str): File name to store your workspace.
            path (str, optional): Directory name to save your workspace.
        """
        print (f'{'\033[1;36m'}kytools: Importing contents to workspace {self._workspace.GetName()}.{'\033[0m'}')
        for _, data in self._data.items():
            getattr(self._workspace, 'import')(data)
        for _, pdf in self._pdf.items():
            getattr(self._workspace, 'import')(pdf)
        os.makedirs(path, exist_ok=True)
        save_path = os.path.join(path, fname)
        self._workspace.writeToFile(save_path)
        print (f'{'\033[1;36m'}kytools: Saved workspace {self._workspace.GetName()} to {save_path}.{'\033[0m'}')