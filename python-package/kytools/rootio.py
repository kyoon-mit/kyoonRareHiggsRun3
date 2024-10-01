import uproot
import pandas as pd
from ROOT import RDataFrame, std

def root_to_pandas(filename, treename, branches):
    """Converts ROOT NTuple to Pandas DataFrame format.

    Along with the pandas DataFrame, a Numpy array of labels is returned. The length
    of this array is the number of rows. The values are 1 if signal=True, and -1 if
    signal=False.

    Args:
        filename (str or list(str)): Name(s) of the ROOT file(s).
        treename (str): Name of the TTree.
        branches (list(str)): List of the branches to convert.

    Returns:
        df (pandas.DataFrame): Pandas DataFrame of the chosen variables.

    Raises:
        TypeError: If filename is neither str nor list(str).
    """
    df_list = []
    if type(filename) is list:
        for f_ in filename:
            with uproot.open('{}:{}'.format(f_, treename)) as events:
                df_list.append(events.arrays(branches, library='pd'))
    elif type(filename) is str:
        with uproot.open('{}:{}'.format(filename, treename)) as events:
            df_list.append(events.arrays(branches, library='pd'))
    else:
        raise TypeError('Argument provided for filename is neither str nor list(str)')
    df_ = pd.concat(df_list, ignore_index=True)
    return df_

def root_to_rdf(filename, treename, defines=[], redefines=[]):
    """Converts ROOT NTuple to ROOT RDataFrame format.

    If `defines` is provided, the RDataFrame will contain new column definitions
    along with the old.

    Args:
        filename (str or list(str)): Name(s) of the ROOT file(s).
        treename (str): Name of the TTree.
        defines (list((str, str)), optional): List containing new column definitions.
            Defaults to [].
        redefines (list((str, str)), optional): List containing column redefinitions.
            Useful for converting types. Defaults to [].

    Returns:
        rdf (ROOT.RDataFrame): ROOT RDataFrame object.

    Raises:
        TypeError: If filename is neither str nor list(str).
    """
    if type(filename) is list:
        filenames_C = std.vector('string')()
        print('rootio: adding the following files to RDataFrame.')
        for fname in filename:
            print(f'\t{fname}')
            filenames_C.push_back(fname)
        rdf_ = RDataFrame(treename, filenames_C)
    elif type(filename) is str:
        print('rootio: adding the following file to RDataFrame.')
        print(f'\t{filename}')
        rdf_ = RDataFrame(treename, filename)
    else:
        raise TypeError('Argument provided for filename is neither str nor list(str)')
    for def_ in defines:
        rdf_ = rdf_.Define(def_[0], def_[1])
    for redef_ in redefines:
        rdf_ = rdf_.Redefine(redef_[0], redef_[1])
    return rdf_