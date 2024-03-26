import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap

def hist2d_heatmap(x_data, y_data, bins, save_name, title, x_title, y_title, aspect=1, dpi=300, range='default'):
    """Creates a heatmap on a 2d histogram.

    Args:
        x_data (array_like, shape (N,)): Array of x coordinates
        y_data (array_like, shape (N,)): Array of y coordinates
        bins (int or array_like or [int, int] or [array, array]):
            The bin specification:
            If int, the number of bins for the two dimensions (nx=ny=bins).
            If array_like, the bin edges for the two dimensions (x_edges=y_edges=bins).
            If [int, int], the number of bins in each dimension (nx, ny = bins).
            If [array, array], the bin edges in each dimension (x_edges, y_edges = bins).
            A combination [int, array] or [array, int], where int is the number
            of bins and array is the bin edges.
        save_name ('str'): Full path to save the plot, including name and extesion.
        title ('str'): Title of the plot.
        x_title ('str'): Title of the x axis.
        y_title ('str'): Title of the y axis.
        aspect (int, optional): Aspect ratio of the plot. Defaults to 1.
        dpi (int, optional): Resolution of the plot to save. Defaults to 300.
        range (str or [float, float, float, float], optional):
            The range of the plot:
            If str and 'default', the range is set to the min and max values of each data.
            If [float, float, float, float], the numbers represent min(x), max(x),
            min(y), max(y).
            Defaults to 'default'.

    Returns:
        H (numpy.ndarray): 2D array of shape(nx, ny) containing the coordinates.
    """
    H, xedges, yedges = np.histogram2d(x_data, y_data, bins=bins)
    if range=='default':
        range_values=[xedges[0], xedges[-1], yedges[0], yedges[-1]]
    elif len(range)==4 and all((isinstance(x, float) or isinstance(x, int)) for x in range):
        range_values=range
    else:
        raise Exception('Please provide a valid argument for `range`.')
    plt.clf()
    plt.imshow(H.T, extent=range_values, origin='lower', aspect=aspect, cmap=cmap.hot)
    plt.title(title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.colorbar()
    plt.show()
    plt.savefig(save_name)
    print(f'Created figure {save_name}')

    return H

def matrix_heatmap(H, range, bins, save_name, title, x_title, y_title, aspect=1, dpi=300):
    """Creates a heatmap from a 2d array.

    Args:
        H (numpy.ndarray): 2D array of shape(nx, ny) containing the coordinates.
        range ([float, float, float, float]):
            The range of the plot. The numbers represent min(x), max(x),
            min(y), max(y).
        bins (int or array_like or [int, int] or [array, array]):
            The bin specification:
            If int, the number of bins for the two dimensions (nx=ny=bins).
            If array_like, the bin edges for the two dimensions (x_edges=y_edges=bins).
            If [int, int], the number of bins in each dimension (nx, ny = bins).
            If [array, array], the bin edges in each dimension (x_edges, y_edges = bins).
            A combination [int, array] or [array, int], where int is the number
            of bins and array is the bin edges.
        save_name ('str'): Full path to save the plot, including name and extesion.
        title ('str'): Title of the plot.
        x_title ('str'): Title of the x axis.
        y_title ('str'): Title of the y axis.
        aspect (int, optional): Aspect ratio of the plot. Defaults to 1.
        dpi (int, optional): Resolution of the plot to save. Defaults to 300.

    Returns:
        H (numpy.ndarray): The input argument H.
    """
    if len(range)==4 and all((isinstance(x, float) or isinstance(x, int)) for x in range):
        range_values=range
    else:
        raise Exception('Please provide a valid argument for `range`.')
    plt.clf()
    plt.imshow(H.T, extent=range_values, origin='lower', aspect=aspect, cmap=cmap.hot)
    plt.title(title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.colorbar()
    plt.show()
    plt.savefig(save_name)
    print(f'Created figure {save_name}')

    return H