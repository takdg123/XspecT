import matplotlib.pyplot as plt
import numpy as np
from gbm.data.containers import RateHisto
from .globals import *


def histo(ratehistos, ax, **kwargs):
    """Plot a series of Ratehistos, either lightcurves or count spectra
    
    Parameters:
    -----------
    ratehistos: list of RateHistos
        The Rate Histogram objects to plot
    ax: matplotlib.axes
        The axis on which to plot
    **kwargs: optional
        Other plotting options
    
    Returns:
    -----------
    refs: list
        The reference to the plot objects
    """    
    refs = []
    for ratehisto in ratehistos:
        # Plot the data
        rate = ratehisto.rate
        rate = np.append(rate, ratehisto.rate[-1])
        p = ax.step(ratehisto.edges, rate, where='post', **kwargs)
        refs.append(p)
    return refs
    
def histo_errorbars(ratehistos, ax, **kwargs):
    """Plot errorbars for a series of RateHistos
    
    Parameters:
    -----------
    ratehistos: list of RateHistos
        The Rate Histogram objects to plot
    ax: matplotlib.axes
        The axis on which to plot
    
    Returns:
    -----------
    refs: list
        The reference to the errorbar objects
    """    
    refs = []
    for ratehisto in ratehistos:
        p = ax.errorbar(ratehisto.bin_centers(), ratehisto.rate, 
                        ratehisto.uncertainty, capsize=0, fmt='none', **kwargs)
        refs.append(p)
    return refs

def filled_histo(ratehistos, ax, color=DATA_SELECTED_COLOR, 
                 fill_alpha=DATA_SELECTED_ALPHA, **kwargs):
    """Plota filled histogram
    
    Parameters:
    -----------
    ratehistos: list of RateHistos
        The Rate Histogram objects to plot
    ax: matplotlib.axes
        The axis on which to plot
    color: str, optional
        The color of the filled histogram
    fill_alpha: float, optional
        The opacity of the fill
    
    Returns:
    -----------
    refs1:
        Reference to the histogram
    refs2:
        Reference to the errorbars
    refs3:
        Reference to the lower boundary
    refs4:
        Reference to the upper boundary
    refs5:
        Reference to the fill artist
    """     
    # plot the histogram
    refs1 = histo(ratehistos, ax, c=color, **kwargs, zorder=3)
    refs1 = [ref[0] for ref in refs1]
    
    # plot the errorbars
    refs2 = histo_errorbars(ratehistos, ax, ecolor=color, **kwargs, zorder=1)
    
    refs3 = []
    refs4 = []
    refs5 = []
    for ratehisto in ratehistos:
        zeros = np.zeros(ratehisto.edges.size)
        zeros[:] = -1e10
        rate = ratehisto.rate
        rate = np.append(rate, ratehisto.rate[-1])
        
        # Change the color of the boundaries of the selected data     
        p = ax.plot([ratehisto.edges[0], ratehisto.edges[0]], [zeros[0], rate[0]], 
                    color=color, zorder=3)
        refs3.append(p[0])
        p = ax.plot([ratehisto.edges[-1], ratehisto.edges[-1]], [zeros[-1], rate[-1]], 
                    color=color, zorder=3)
        refs4.append(p[0])
        
        # do the fill
        p = ax.fill_between(ratehisto.edges, rate, zeros, color=color, step='post',
                            alpha=fill_alpha, **kwargs, zorder=4)
        refs5.append(p)
    
    return (refs1, refs2, refs3, refs4, refs5)
    
def selection_line(xpos, ax, **kwargs):
    """Plot a selection line
    
    Parameters:
    -----------
    xpos: float
        The position of the selection line
    ax: matplotlib.axes
        The axis on which to plot
    
    Returns:
    -----------
    refs: list
        The reference to selection line
    """    
    ylim = ax.get_ylim()
    ref = ax.plot([xpos, xpos], ylim, **kwargs)
    return ref

def selections(bounds, ax, **kwargs):
    """Plot selection bounds
    
    Parameters:
    -----------
    bounds: list of tuples
        The list of selection bounds
    ax: matplotlib.axes
        The axis on which to plot
    
    Returns:
    -----------
    refs1:
        The reference to the lower selection line
    refs2:
        The reference to the upper selection line
    """        
    refs1 = []
    refs2 = []
    for bound in bounds:
        p = selection_line(bound[0], ax, **kwargs)
        refs1.append(p[0])
        p = selection_line(bound[1], ax, **kwargs)
        refs2.append(p[0])
    
    return (refs1, refs2)

def errorband(x, y_upper, y_lower, ax, **kwargs):
    """Plot an error band
    
    Parameters:
    -----------
    x: np.array
        The x values
    y_upper: np.array
        The upper y values of the error band
    y_lower: np.array
        The lower y values of the error band
    ax: matplotlib.axes
        The axis on which to plot
    
    Returns:
    -----------
    refs: 
        The reference to error band artist
    """        
    refs = ax.fill_between(x, y_upper, y_lower, **kwargs)
    return refs

def lightcurve_background(background_model, selection, ax, cent_alpha=None, 
                          err_alpha=None, **kwargs):
    """Plot a lightcurve background model with an error band
    
    Parameters:
    -----------
    background_model: Background
        A background fitted background model
    selection: Selection
        A lightcurve selection over which the background model will be plotted    
    ax: matplotlib.axes
        The axis on which to plot
    cent_alpha: float, optional
        The opacity of the centroid line
    err_alpha: float, optional
        The opacity of the error band
    
    Returns:
    -----------
    centroid:
        Reference to the plotted centroid
    sigma:
        Reference to the error band artist
    """            
    if 'zorder' not in kwargs:
        kwargs['zorder'] = 1
        
    b, u = background_model.integrate_energy(selection)
    times = []
    for lightcurves in selection.lightcurve:
        times.extend(lightcurves.bin_centers())
    # plot the centroid line
    centroid = ax.plot(times, b, alpha=cent_alpha, **kwargs)
    
    # plot the error band
    kwargs['linestyle'] = '-'
    kwargs['zorder'] -= 1
    sigma = errorband(times, b+u, b-u, ax, alpha=err_alpha, **kwargs)
    
    return (centroid, [sigma])

def spectrum_background(background_model, selection, ax, cent_alpha=None, 
                        err_alpha=None, **kwargs):
    """Plot a count spectrum background model with an error band
    
    Parameters:
    -----------
    background_model: Background
        A background fitted background model
    selection: Selection
        A count spectrum selection over which the background model will be plotted    
    ax: matplotlib.axes
        The axis on which to plot
    cent_alpha: float, optional
        The opacity of the centroid line
    err_alpha: float, optional
        The opacity of the error band
    
    Returns:
    -----------
    centroid:
        Reference to the plotted centroid
    sigma:
        Reference to the error band artist
    """        
    if 'zorder' not in kwargs:
        kwargs['zorder'] = 1

    bkgd_spectrum = background_model.source_spectrum_integrated(selection)
    # have to multiply by binwidth this is because Rob demanded that he be 
    # given counts/sec for xspec
    bkgd_spectrum.exposure = bkgd_spectrum.exposure * bkgd_spectrum.bin_widths()
    bkgd_spectrum.rate = RateHisto.calc_rate(bkgd_spectrum.values, 
                                             bkgd_spectrum.exposure)
    
    # plot the centroid line
    centroid = histo([bkgd_spectrum], ax, alpha=cent_alpha, **kwargs)
    
    # have to create the stepped functions so that we can fill between them
    t = np.array(bkgd_spectrum.bounds).T.flatten().tolist()
    r1 = [bkgd_spectrum.rate + bkgd_spectrum.uncertainty,
          bkgd_spectrum.rate + bkgd_spectrum.uncertainty]
    r1 = np.array(r1).T.flatten().tolist()
    r2 = [bkgd_spectrum.rate - bkgd_spectrum.uncertainty,
          bkgd_spectrum.rate - bkgd_spectrum.uncertainty]
    r2 = np.array(r2).T.flatten().tolist()
    
    # plot the error band
    kwargs['linestyle'] = '-'
    kwargs['zorder'] -= 1
    sigma = errorband(t, r1, r2, ax, alpha=err_alpha, **kwargs)
    
    return (centroid[0], [sigma])
