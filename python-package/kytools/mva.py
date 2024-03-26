"""
mva
=========================================

A module to handle MVA training and validation.

"""
import ROOT

# TODO: fix path to make it more general
ROOT.gInterpreter.Declare('#include "/work/submit/kyoon/HiggsRareDecay/python-package/kytools/tmva_helper_xml.h"')

class TMVAHelperXML():

    def __init__(self, model_input, model_name=""):
        model_tmp = ROOT.TMVA.Experimental.RReader(model_input)
        self.variables = [str(var) for var in model_tmp.GetVariableNames()]
        self.model_input = model_input
        self.model_name = model_name
        self.nthreads = ROOT.GetThreadPoolSize()

        self.tmva_helper = ROOT.tmva_helper_xml(self.model_input, self.nthreads)
        self.var_col = f"tmva_vars_{self.model_name}"

    def run_inference(self, rdf, col_name = "mva_score"):
        # check if columns exist in the dataframe
        cols = rdf.GetColumnNames()
        for var in self.variables:
            if not var in cols:
                raise Exception(f"Variable {var} not defined in dataframe.")

        vars_str = ', '.join(self.variables)
        rdf = rdf.Redefine(self.var_col, f"ROOT::VecOps::RVec<float>{{{vars_str}}}")
        rdf = rdf.RedefineSlot(col_name, self.tmva_helper, [self.var_col])
        return rdf
    
    def run_inference_deprecated(self, rdf, col_name = 'mva_score'):
         # check if columns exist in the dataframe
        cols = rdf.GetColumnNames()
        for var in self.variables:
            if not var in cols:
                raise Exception(f"Variable {var} not defined in dataframe.")
            
        process_string = f'''
        TMVA::Experimental::RReader model("{self.model_input}");
        computeModel = TMVA::Experimental::Compute<{len(self.variables)}, float>(model);
        '''
        print(process_string)
        ROOT.gInterpreter.ProcessLine(process_string)

        rdf = rdf.Redefine(col_name, ROOT.computeModel, self.variables)
        return rdf