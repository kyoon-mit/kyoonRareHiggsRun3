from os import path, environ, makedirs
from kytools import rootio
from matplotlib import pyplot as plt
import numpy as np
# TODO: Parse arguments

def get_pandas_dataframes(sig_filename, bkg_filename, tree_name, branches, weight_var='', **kwargs):
    """Converts the given TTree in ROOT files into pandas.DataFrames.

    This method opens the ROOT files, reads the values in the specified TTree, and converts them into
    pandas.DataFrame objects. Two dataframes are returned: one for signal and
    another for background.

    Args:
        sig_filename (str): Name of the ROOT file containing the signal snapshot.
        bkg_filename (str): Name of the ROOT file containing the background snapshot.
        tree_name (str): Name of the TTree in the ROOT file.
        branches (list(str)): The MVA input variables to include in the pandas.DataFrame.
            It is okay to pass variables that are not used in training.
        weight_var (str, optional): If provided, the function will look for the
            weight column of that name in the pandas.DataFrame.
            Defaults to ''.
        **kwargs (optional): Other arguments that are passed but not used.
            This is to avoid getting an error in the case where **dict is passed and
            there are too many keys.
    
    Returns:
        sig_df (pandas.DataFrame): DataFrame of the signal dataset.
        bkg_df (pandas.DataFrame): DataFrame of the background dataset.
    """
    print('INFO: executing get_pandas_dataframes.')
    if weight_var: branches += weight_var
    sig_df, bkg_df =\
        rootio.root_to_pandas(sig_filename, tree_name, branches),\
        rootio.root_to_pandas(bkg_filename, tree_name, branches)
    return sig_df, bkg_df


def compute_correlation(sig_df, bkg_df, corr_var, save_path, matrix=True, scatter=True,
                        heatmap=True, terminal=False, bin_range=[[-1., 1.], [-1., 1.]],
                        weight_var='', sfx='',
                        **kwargs):
    """Computes numerical and graphical correlations from the given pandas.DataFrames.

    Each and every variable in the pandas.DataFrame object is compared against
    the corr_var variable. The results are shown for the background
    and signal dataset each.

    If `matrix` is True, it computes the correlation matrix. The matrix is saved
    as a text file in the specified directory and optionally printed in the terminal.

    If `scatter` is True, scatter plots are created and saved. For each column
    that is common in both bkg_df and sig_df, three scatter plots are shown:
    signal, background, and combined (in different colors).

    The option 'heatmap' is similar to 'scatter', except that a heatmap is shown.

    If `slices` is True, a histogram of the MVA discrimininator is created in binned ranges
    of each variable.

    Args:
        sig_df (pandas.DataFrame): pandas.DataFrame of the signal dataset.
        bkg_df (pandas.DataFrame): pandas.DataFrame of the background dataset.
        corr_var (str): Name of the correlation variable in the pandas.DataFrame
            objects to be compared against.
        save_path (str): Directory where the outputs of this method is saved.
        matrix (bool, optional): Whether to compute the correlation matrix.
            Defaults to True.
        scatter (bool, optional): Whether to plot the scatter plots.
            Defaults to True.
        heatmap (bool, optional): Wheter to plot the heatmaps.
            Defaults to True.
        terminal (bool, optional): Whether to show results in the terminal.
            Defaults to True.
        bin_range ([[float, float], [float, float]], optional): Used only if
            'heatmap' is True. Provides the boundaries of the outermost bins in
            the format of [[xmin, xmax], [ymin, ymax]].
            Defaults to [[-1., 1.], [-1., 1.]].
        weight_var (str, optional): If provided, the function will look for the
            weight column of that name in the pandas.DataFrame. Relevant only if
            'heatmap' is True.
            Defaults to ''.
        sfx (str, optional): Suffix of the file name.
            Defaults to ''.
        **kwargs (optional): Other arguments that are passed but not used.
            This is to avoid getting an error in the case where **dict is passed and
            there are too many keys.
    
    Returns:
        None
    """
    if not path.exists(save_path): makedirs(save_path)
    print('INFO: executing compute_correlation.')
    # Calculate correlation matrix
    if matrix:
        print('INFO: calculating correlation matrix.')
        sig_corr = sig_df.corr(method='pearson').to_string()
        bkg_corr = bkg_df.corr(method='pearson').to_string()
        text_file_name = path.join(save_path, f'CORRELATION_MATRIX{'_' if sfx else ''}{sfx}.txt')
        text_input = [f'<Signal>\n{sig_corr}',
                      '\n',
                      f'<Background>\n{bkg_corr}']
        with open(text_file_name, 'wt') as txtfile:
            txtfile.writelines(text_input)
        if terminal:
            print(*text_input, sep='\n')
    # Create scatter plots
    if scatter:
        print('INFO: creating scatter plot.')
        sig_columns, bkg_columns = list(sig_df), list(bkg_df)
        sig_columns.remove(corr_var)
        bkg_columns.remove(corr_var)
        for col in list(set(sig_columns).intersection(bkg_columns)):
            ax_sig = sig_df.plot.scatter(x=col, y=corr_var, c='Red', s=2.5, alpha=0.01)
            ax_sig.set_title(f'Signal: \n{corr_var} vs. {col}')
            plt.savefig(path.join(save_path, f'SCATTER_PLOT_sig_{corr_var}_VS_{col}{'_' if sfx else ''}{sfx}.png'),
                        dpi=400, pad_inches=0.05)
            plt.close()
            ax_bkg = bkg_df.plot.scatter(x=col, y=corr_var, c='Blue', s=2.5, alpha=0.002)
            ax_bkg.set_title(f'Background: \n{corr_var} vs. {col}')
            plt.savefig(path.join(save_path, f'SCATTER_PLOT_bkg_{corr_var}_VS_{col}{'_' if sfx else ''}{sfx}.png'),
                        dpi=400, pad_inches=0.05)
            plt.close()
            ax_comb = plt.subplot(111)
            ax_comb.scatter(x=bkg_df[col], y=bkg_df[corr_var], c='Blue', s=2.5, alpha=0.002, label='background')
            ax_comb.scatter(x=sig_df[col], y=sig_df[corr_var], c='Red', s=2.5, alpha=0.01, label='signal')
            ax_comb.set_title(f'Signal & Background: \n{corr_var} vs. {col}')
            ax_comb.set_xlabel(col)
            ax_comb.set_ylabel(corr_var)
            for lh in ax_comb.legend().legend_handles: lh.set_alpha(1)
            plt.savefig(path.join(save_path, f'SCATTER_PLOT_comb_{corr_var}_VS_{col}{'_' if sfx else ''}{sfx}.png'),
                        dpi=400, pad_inches=0.05)
            plt.close()
    # Create heatmaps
    if heatmap:
        sig_columns, bkg_columns = list(sig_df), list(bkg_df)
        sig_columns.remove(corr_var)
        bkg_columns.remove(corr_var)
        if weight_var:
            sig_columns.remove(weight_var)
            bkg_columns.remove(weight_var)
        for col in list(set(sig_columns).intersection(bkg_columns)):
            print('INFO: creating heatmap.')
            weights = bkg_df[weight_var] if weight_var else None
            heatmap, xedges, yedges = np.histogram2d(x=bkg_df[col], y=bkg_df[corr_var],
                                                    bins=40, range=bin_range, weights=weights)
            # heatmap = np.array([[1, 4], [2, 3]]) # for testing
            # xedges, yedges = [0., 0.5, 1.], [0., 0.5, 1.]
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            plt.clf()
            _, ax_bkg = plt.subplots()
            im = plt.imshow(heatmap.T, extent=extent, cmap='viridis', origin='lower')
            ax_bkg.set_title(f'Background: \n{corr_var} vs. {col}')
            ax_bkg.set_xlabel(col)
            ax_bkg.set_ylabel(corr_var)
            plt.colorbar(im)
            plt.savefig(path.join(save_path, f'HEATMAP_bkg_{corr_var}_VS_{col}{'_' if sfx else ''}{sfx}.png'),
                        dpi=400, pad_inches=0.05)
            plt.close()
            plt.clf()
            #
            weights = sig_df[weight_var] if weight_var else None
            heatmap, xedges, yedges = np.histogram2d(x=sig_df[col], y=sig_df[corr_var],
                                                    bins=40, range=bin_range, weights=weights)
            _, ax_sig = plt.subplots()
            im = plt.imshow(heatmap.T, extent=extent, cmap='viridis', origin='lower')
            ax_sig.set_title(f'Signal: \n{corr_var} vs. {col}')
            ax_sig.set_xlabel(col)
            ax_sig.set_ylabel(corr_var)
            plt.colorbar(im)
            plt.savefig(path.join(save_path, f'HEATMAP_sig_{corr_var}_VS_{col}{'_' if sfx else ''}{sfx}.png'),
                        dpi=400, pad_inches=0.05)
            plt.close()
    return
    

def sliced_hams(sig_df, bkg_df, x_col, y_col, xbins, save_path, **kwargs):
    """Creates histograms of the variable y in bins of the variable x and plots
    them on top of each other.

    This method is similar to projecting the scatter plot as a histogram to the 
    y axis, but instead of summing over the entire range of x, the data counts
    are summed over sub-ranges of x.

    Three plots are created: histograms that are stacked in steps, bar charts of
    different sample sizes shown in parallel, and normalized bar charts shown in
    parallel.

    The resultings plots are saved in the specified directory.

    Args:
        sig_df (pandas.DataFrame): pandas.DataFrame of the signal dataset.
        bkg_df (pandas.DataFrame): pandas.DataFrame of the background dataset.
        x_col (str): The column name of x, the variable to constrain the sub-ranges.
        y_col (str): The column name of y, or the variable to be plotted in the x-axis
            in the histograms
        xbins (list(float)): The boundaries of the x sub-ranges. The list must have
            at least two elements. For a list [x_1, ..., x_n], the sub-ranges are
            [x_1, x_2), ..., [x_(n-1), x_n).
        save_path (str): Directory where the outputs of this method is saved.
        **kwargs (optional): Other arguments that are passed but not used.
            This is to avoid getting an error in the case where **dict is passed and
            there are too many keys.

    Raises:
        ValueError: If the length of xbins is less than 2, or if the values are
            not strictly monotonically increasing.

    Returns:
        None
    """
    if not path.exists(save_path): makedirs(save_path)
    if len(xbins) < 2:
        raise ValueError('The list must have at least two elements.')
    sig_y_values, bkg_y_values = [], []
    sliced_labels = []
    for i in range(1, len(xbins)):
        xlow = xbins[i-1]
        xhigh = xbins[i]
        if xlow >= xhigh:
            raise ValueError('The list must be strictly monotonically increasing.')
        sliced_labels.append(f'[{xlow:.1f}, {xhigh:.1f})')
        sig_sliced_ham = sig_df.loc[(xlow <= sig_df[x_col]) & (sig_df[x_col] < xhigh)]
        bkg_sliced_ham = bkg_df.loc[(xlow <= bkg_df[x_col]) & (bkg_df[x_col] < xhigh)]
        sig_y_values.append(sig_sliced_ham[y_col])
        bkg_y_values.append(bkg_sliced_ham[y_col])
    # TODO: change the following lines to for loop
    # Signal unnormalized bar charts
    plt.hist(sig_y_values, density=False, histtype='bar', label=sliced_labels)
    plt.legend(prop={'size': 10})
    plt.title(f'Signal: \n{y_col} in bins of {x_col}')
    plt.xlabel(y_col)
    plt.ylabel('Events')
    plt.savefig(path.join(save_path, f'HISTOGRAM_UNNORMED_sig_{y_col}_BINS_{x_col}.png'),
                dpi=1000, pad_inches=0.05)
    plt.close()
    # Signal normalized bar charts
    plt.hist(sig_y_values, density=True, histtype='bar', label=sliced_labels)
    plt.legend(prop={'size': 10})
    plt.title(f'Signal: \nNormalized {y_col} in bins of {x_col}')
    plt.xlabel(y_col)
    plt.ylabel('Events')
    plt.savefig(path.join(save_path, f'HISTOGRAM_NORMED_sig_{y_col}_BINS_{x_col}.png'),
                dpi=1000, pad_inches=0.05)
    plt.close()
    # Signal stacked-stepped histograms
    plt.hist(sig_y_values, density=True, histtype='step', label=sliced_labels)
    plt.legend(prop={'size': 10})
    plt.title(f'Signal: \nNormalized {y_col} in bins of {x_col}')
    plt.xlabel(y_col)
    plt.ylabel('Events')
    plt.savefig(path.join(save_path, f'HISTOGRAM_STEP_NORMED_sig_{y_col}_BINS_{x_col}.png'),
                dpi=1000, pad_inches=0.05)
    plt.close()
    # Background unnormalized bar charts
    plt.hist(bkg_y_values, density=False, histtype='bar', label=sliced_labels)
    plt.legend(prop={'size': 10})
    plt.title(f'Background: \n{y_col} in bins of {x_col}')
    plt.xlabel(y_col)
    plt.ylabel('Events')
    plt.savefig(path.join(save_path, f'HISTOGRAM_UNNORMED_bkg_{y_col}_BINS_{x_col}.png'),
                dpi=1000, pad_inches=0.05)
    plt.close()
    # Background normalized bar charts
    plt.hist(bkg_y_values, density=True, histtype='bar', label=sliced_labels)
    plt.legend(prop={'size': 10})
    plt.title(f'Background: \nNormalized {y_col} in bins of {x_col}')
    plt.xlabel(y_col)
    plt.ylabel('Events')
    plt.savefig(path.join(save_path, f'HISTOGRAM_NORMED_bkg_{y_col}_BINS_{x_col}.png'),
                dpi=1000, pad_inches=0.05)
    plt.close()
    # Background stack-stepped bar charts
    plt.hist(bkg_y_values, density=True, histtype='step', label=sliced_labels)
    plt.legend(prop={'size': 10})
    plt.title(f'Background: \nNormalized {y_col} in bins of {x_col}')
    plt.xlabel(y_col)
    plt.ylabel('Events')
    plt.savefig(path.join(save_path, f'HISTOGRAM_STEP_NORMED_bkg_{y_col}_BINS_{x_col}.png'),
                dpi=1000, pad_inches=0.05)
    return

if __name__=='__main__':
    save_path = path.abspath(path.join(environ['HRARE_DIR'], 'JPsiCC', 'plots', 'v202410'))

    ### Parameters
    params = {
        'save_path': save_path,
        'sig_filename': '/work/submit/kyoon/CMSSW_13_3_0/src/RareHiggsRun3/JPsiCC/analysis/snapshot_MC_SIG_2018_GF_v202410_CMSSW_13_3_0_20241120_WEIGHT_no_filter.root',
        'bkg_filename': '/work/submit/kyoon/CMSSW_13_3_0/src/RareHiggsRun3/JPsiCC/analysis/snapshot_MC_BKG_2018_GF_v202410_CMSSW_13_3_0_20241120_WEIGHT_no_filter.root',
        'tree_name': 'Events',
        'branches': ['jetFar_CvL', 'jetFar_CvB', 'jetClose_CvL', 'jetClose_CvB'],
        'weight_var': 'w',
        'corr_var': 'jetFar_CvL'
    }

    ### Get pandas.DataFrames
    params['sig_df'], params['bkg_df'] = get_pandas_dataframes(**params)

    ### Compute correlation
    compute_correlation(scatter=False, terminal=False, bin_range=[[0., 1.], [0., 1.]], sfx='no_filter_jet_pt_0', **params)
    params['corr_var'] = 'jetClose_CvL'
    compute_correlation(scatter=False, terminal=False, bin_range=[[0., 1.], [0., 1.]], sfx='no_filter_jet_pt_0', **params)