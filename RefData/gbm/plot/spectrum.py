import numpy as np
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as  plt

from .globals import *
from .gspecplot import GspecPlot
from .lib import *
from gbm.data.containers import RateHisto

class Spectrum(GspecPlot):
    """Class for plotting count spectra and count spectra paraphernalia
    """
    def __init__(self, gspec, dataname, canvas=None):
        """
        Parameters:
        --------------
        gspec: GspecManager
            Reference to an initialized GspecManager
        dataname: str
            The data filename to plot.  Must be registered in gspec
        canvas: Canvas Backend object, optional
            If interfacing with a backend (e.g. Tk), pass the relavant Canvas
            Backend object
        """
        super(Spectrum, self).__init__(gspec, dataname, canvas=canvas)
                
        # store the artists/collections for all of the plot items
        # as well as if they are hidden or displayed
        self.settings = {'data': {'artists': None, 'hidden': False},
                         'errorbars': {'artists': None, 'hidden': False},
                         'source': {'artists': None, 'hidden': False},
                         'background': {'artists': None, 'hidden': False},
                         'yscale': 'log', 'xscale': 'log',
                         'xlim': None, 'ylim': None}

    def clear(self, data=True, errorbars=True, source=True, background=True, 
              selections=True, text=True):
        """Clear the plot
        
        Parameters:
        --------------
        data: bool, optional
            Set to True to remove the data. Default is True.
        errorbars: bool, optional
            Set to True to remove the errorbars. Default is True.
        source: bool, optional
            Set to True to remove any sources. Default is True.
        background: bool, optional
            Set to True to remove the background model. Default is True.
        selections: bool, optional
            Set to True to remove any selections. Default is True.
        text: bool, optional
            Set to True to remove any text. Default is True.
        """
        if data:
            self.remove_data()
        if errorbars:
            self.remove_errorbars()
        if source:
            self.remove_source()
        if background:
            self.remove_background()
        if selections:
            self.remove_selection_lines()
        if text:
            self.remove_text()
    
    def plot(self, data=True, errorbars=True, source=True, background=True, 
             view=True, **kwargs):
        """Plot the data, selections, and background
        
        Parameters:
        --------------
        data: bool, optional
            Set to True to plot data. Default is True.
        errorbars: bool, optional
            Set to True to plot errorbars. Default is True.
        source: bool, optional
            Set to True to plot any source selections. Default is True.
        background: bool, optional
            Set to True to plot the background model. Default is True.
        view: bool, optional
            Set to True set the view window. Default is True
        """
        if view:
            self.display_lookup_view()

        if data:
            self.plot_data(**kwargs)
        
        if errorbars:
            self.plot_errorbars()
        
        if source:
            self.plot_source()
        
        if background:
            self.plot_background()
        
    def plot_data(self, yscaling=True):
        """Plot the count spectrum data

        Parameters:
        --------------
        yscaling: bool, optional
            Set to True to automatically scale the y-axis to the data range.
            Default is True.
        """
        # Extract the energy spectrum from the file
        ymin = []
        ymax = []
        spectra = self.gspec.data[self.dataname]['energyview'].spectrum
        for spectrum in spectra:
            ymin.append(np.min(spectrum.rate))
            ymax.append(np.max(spectrum.rate))
        
        plotrefs = histo(spectra, self.ax, c=DATA_COLOR, zorder=2)
        self.settings['data']['artists'] = plotrefs

        # set the plot range
        xlim = (spectra[0].bounds[0][0], spectra[-1].bounds[1][-1])
        ylim = (0.90*np.min(ymin), 1.10*np.max(ymax))
        if self.settings['xlim'] is None:
            self.ax.set_xlim(xlim)
            self.settings['xlim'] = xlim
        self.ax.set_ylim(ylim)
        self.settings['ylim'] = ylim

        # Set the x and y labels
        self.ax.set_xlabel('Energy (keV)', fontsize=PLOTFONTSIZE)
        self.ax.set_ylabel('Rate (counts/s-keV)', fontsize=PLOTFONTSIZE)

        # Set the tick label size
        self.ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self.ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        
        # Set the x and y scales
        self.ax.set_xscale(self.settings['xscale'])
        if yscaling:
            self.ax.set_yscale(self.settings['yscale'])
        
        if self.canvas is None:
            plt.draw()

    def plot_errorbars(self):
        """Plot the count spectrum errorbars
        """
        # Extract the energy spectrum from the file
        spectra = self.gspec.data[self.dataname]['energyview'].spectrum
        
        plotrefs = histo_errorbars(spectra, self.ax, ecolor=DATA_ERROR_COLOR, 
                                   alpha=DATA_ERROR_ALPHA, zorder=1)
        self.settings['errorbars']['artists'] = plotrefs

        if self.canvas is None:
            plt.draw()

    def plot_source(self):
        """Plot the source selections
        """
        # Display any existing user selections
        if self.gspec.data[self.dataname]['timeview'] is None:
            #print('No energy selection available for {0}'.format(self.dataname))
            return
        
        # Get the current axis limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        
        # plot the selections
        spectra = self.gspec.data[self.dataname]['timeview'].spectrum
        bounds = [(s.bounds[0][0], s.bounds[1][-1]) for s in spectra]
        plotrefs = selections(bounds, self.ax, c=DATA_SELECTED_COLOR, 
                              linestyle='--', linewidth=0.75) 
        
        # Set the axis limits to what they were before plotting the selection
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        
        self.settings['source']['artists'] = plotrefs#(plotrefs1, plotrefs2)
        
        if self.canvas is None:
            plt.draw()
    
    def plot_background(self):
        """Plot the background
        """
        # Display the background model, if defined
        if self.gspec.data[self.dataname]['bkgdmodel'] is None:
            #print('No background model available for {0}'.format(self.dataname))
            return
        
        # Display the background model, if defined
        background_model = self.gspec.data[self.dataname]['bkgdmodel']
        if self.gspec.data[self.dataname]['sourceview'] is not None:
            selection = self.gspec.data[self.dataname]['sourceview']
        else:
            selection = self.gspec.data[self.dataname]['timeview']
        plotrefs = spectrum_background(background_model, selection, self.ax,
                                       cent_alpha=BKGD_ALPHA, 
                                       err_alpha=BKGD_ERROR_ALPHA, zorder=11, 
                                       color=BKGD_COLOR, linewidth=BKGD_WIDTH)
        
        self.settings['background']['artists'] = plotrefs
        
        if self.canvas is None:
            plt.draw()

    def plot_source_selection_line(self, xpos):
        """Plot a source selection line
        
        Parameters:
        -----------
        xpos: float
            The position of the selection line
        """
        self.plot_selection_line(xpos, linestyle='--', c=DATA_SELECTED_COLOR, 
                                 linewidth=BKGD_WIDTH)    
    
    def plot_binning_selection_line(self, xpos):
        """Plot a bining selection line
        
        Parameters:
        -----------
        xpos: float
            The position of the binning selection line
        """
        self.plot_selection_line(xpos, linestyle='--', c=BINNING_SELECTED_COLOR, 
                                 linewidth=BKGD_WIDTH)  


    def display_lookup_view(self):
        """Set the view range based on the lookup
        """
        view = self.gspec.lookup[self.dataname].views.energy
        if view is None:
            #print('No view range in lookup for {0}'.format(self.dataname))
            return

        xlim = view.xrange()
        ylim = view.yrange()

        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.settings['xlim'] = xlim
        self.settings['ylim'] = ylim
        
        if self.canvas is None:
            plt.draw()
    
    def remove_data(self):
        """Remve the data display
        """
        self._remove(self.settings['data'])
        self.remove_errorbars()

    def remove_errorbars(self):
        """Remove the errorbars display
        """
        self._remove(self.settings['errorbars'])

    def remove_source(self):
        """Remove the source selections display
        """
        self._remove(self.settings['source'])

    def remove_background(self):
        """Remove the background display
        """
        self._remove(self.settings['background'])
   
    def toggle_data(self):
        """Toggle the data display
        """
        self._toggle(self.settings['data'])

    def toggle_errorbars(self):
        """Toggle the errorbars display
        """
        self._toggle(self.settings['errorbars'])
        
    def toggle_source(self):
        """Toggle the energy selections display
        """
        self._toggle(self.settings['source'])
            
    def toggle_background(self):
        """Toggle the background display
        """
        self._toggle(self.settings['background'])
    