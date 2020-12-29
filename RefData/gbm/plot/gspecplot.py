from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as  plt
from .lib import *

class GspecPlot(object):
    """Class for plotting lightcurves and lightcurve paraphernalia
    """
    
    def __init__(self, gspec, dataname, canvas=None):
        """
        Parameters:
        --------------
        gspec: GspecManager
            Reference to an initialized GspecManager
        dataname: str
            The data filename to plot.  Must be registered in gspec
        canvas: (Canvas Backend object, master window), optional
            If interfacing with a backend (e.g. Tk), pass the relavant Canvas
            Backend object and the master widget window
        *args:
            Options to be passed to the canvas
        """
        self.gspec = gspec
        self.dataname = dataname
        self.canvas = None
        
        # if using a backend canvas, like Tk       
        if canvas is not None:
            self.figure = Figure(figsize=(5, 4), dpi=100)
            # if tkagg, remember to assign to frame
            master = canvas[1]
            canvas = canvas[0]
            self.canvas = canvas(self.figure, master)
        else:
            # just use pyplot, but disable blocking
            self.figure = plt.figure(figsize=(7.7, 4.7), dpi=100)
            self.canvas = self.figure.canvas
            plt.ion()
            
        self.ax = self.figure.add_subplot(111)
        
        # Set the format of the coordinate readout
        self.ax.format_coord = lambda x, y: ""

        # Set the x minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.xaxis.set_minor_locator(minorLocator)

        # Set the y minor tick frequency
        minorLocator = AutoMinorLocator()
        self.ax.yaxis.set_minor_locator(minorLocator)
        
        self.annotations = {'selections': {'artists': [], 'hidden': False}, 
                            'selectfills': {'artists': [], 'hidden': False},
                            'text': {'artists': None, 'hidden': False}}
        self.settings = {'xscale': 'linear', 'yscale': 'linear', 
                         'xlim': None, 'ylim': None}
    
    def plot_selection_line(self, xpos, **kwargs):
        """Plot a selection line, and if we have a pair then fill between
        
        Parameters:
        -----------
        xpos: float
            The position of the selection line
        kwargs: dict
            The line properties
        """
        ylim = self.ax.get_ylim()
        plotref = selection_line(xpos, self.ax, **kwargs)
        self.ax.set_ylim(ylim)
        self.annotations['selections']['artists'].append(plotref)
        
        # if we have a pair of selections, then do a fill to show the range
        if len(self.annotations['selections']['artists']) % 2 == 0:
            sel1 = self.annotations['selections']['artists'][-2][0].get_data()
            sel2 = self.annotations['selections']['artists'][-1][0].get_data()
            color = None
            if 'c' in kwargs:
                color = kwargs['c']
            plotref = self.ax.axvspan(sel1[0][0], sel2[0][0], alpha=0.2, color=color)
            self.annotations['selectfills']['artists'].append([plotref])
            
    def plot_text(self, text, xypos, **kwargs):
        if self.annotations['text']['artists'] is not None:
            self.remove_text()
        artist = self.ax.annotate(text, xy=xypos, **kwargs)
        self.annotations['text']['artists'] = [[artist]]
        
    def toggle_selection_lines(self):
        """Toggle the selection lines
        """
        self._toggle(self.annotations['selections'])
        
    def toggle_text(self):
        """Toggle the text
        """
        self._toggle(self.annotations['text'])
    
    def remove_selection_lines(self):
        """Remove the selection lines
        """
        self._remove(self.annotations['selections'])
        self.annotations['selections']['artists'] = []
        self._remove(self.annotations['selectfills'])
        self.annotations['selectfills']['artists'] = []
    
    def remove_text(self):
        """Remove the text
        """
        self._remove(self.annotations['text'])
    
    def toggle_xscale(self):
        """Toggle the xscale between 'linear' and 'log'
        """
        if self.settings['xscale'] == 'linear':
              self.settings['xscale'] = 'log'
              xlim = self.ax.get_xlim()
              if xlim[0] <= 0.0:
                self.ax.set_xlim((0.01, xlim[1]))
        else:
            self.settings['xscale'] = 'linear'
        self.ax.set_xscale(self.settings['xscale'])

    def toggle_yscale(self):
        """Toggle the yscale between 'linear' and 'log'
        """
        if self.settings['yscale'] == 'linear':
              self.settings['yscale'] = 'log'
              ylim = self.ax.get_ylim()
              if ylim[0] <= 0.0:
                self.ax.set_ylim((0.1, ylim[1]))
        else:
            self.settings['yscale'] = 'linear'
        self.ax.set_yscale(self.settings['yscale'])
    
    def _remove(self, artist_dict):
        """Remove/erase an artist collection
        """        
        if artist_dict['artists'] is None:
            return

        for artist in artist_dict['artists']:
            if not isinstance(artist, list):
                artistlist = list(artist)
            else:
                artistlist = artist
            for element in artistlist:
                if element is not None:
                    if isinstance(element, tuple):
                        if len(element) > 0:
                            try:
                                element[0].remove()
                            except:
                                tuple(element)[-1][0].remove()
                    else:
                        element.remove()
        
        artist_dict['artists'] = None
        artist_dict['hidden'] = False
    
    def _toggle(self, artist_dict):
        """Toggle an artist collection on and off
        """        
        if artist_dict['artists'] is None:
            return
        if artist_dict['hidden']:
            visible = True
        else:
            visible = False
        
        # each artist collection is a collection of artists, and each artist
        # is a collection of elements.  Matplotlib isn't exactly consistent
        # on how each of these artists are organized for different artist 
        # classes, so we have to have some contingencies
        for artist in artist_dict['artists']:
            if not isinstance(artist, list):
                artistlist = list(artist)
            else:
                artistlist = artist
            for element in artistlist:
                try:
                    if isinstance(element, tuple):
                        element[0].set_visible(visible)
                    else:
                        element.set_visible(visible)
                except:
                    pass
        
        artist_dict['hidden'] = not artist_dict['hidden']
        
