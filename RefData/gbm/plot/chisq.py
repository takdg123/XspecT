import numpy as np
import matplotlib.pyplot as  plt

from .globals import *
from .gspecplot import GspecPlot

class BackgroundChisq(GspecPlot):
    """Class for plotting the background reduced chisq
    """
    def __init__(self, gspec, dataname, canvas=None):
        super(BackgroundChisq, self).__init__(gspec, dataname, canvas=canvas)
        self.settings = {'redchisq': {'artists': None, 'hidden': False}, 
                         'yscale': 'linear', 'xscale': 'log'}
    
    def plot_reduced_chisq(self):
        """Plot the reduced chisq as a fucnction of energy
        """
        chisq = self.gspec.data[self.dataname]['bkgdmodel'].chisq
        dof = self.gspec.data[self.dataname]['bkgdmodel'].dof
        bounds = self.gspec.data[self.dataname]['bkgdview']._get_energy_channel_bounds()
        edges = [edge for edge in zip(*bounds)][1]
        redchisq = chisq/dof

        hist = self.ax.step(edges, redchisq, where='post', c=DATA_COLOR, zorder=2)
        self.settings['redchisq']['artists'] = [hist]
        self.settings['redchisq']['hidden'] = False

        # set the plot range
        xlim = (bounds[0][0], bounds[-1][-1])
        ylim = (0.90*np.min(redchisq), 1.10*np.max(redchisq))
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        
        # Set the axis labels
        self.ax.set_ylabel('Chi-Square/DoF', fontsize=PLOTFONTSIZE)
        self.ax.set_xlabel('Energy (keV)', fontsize=PLOTFONTSIZE)
        self.ax.set_xscale(self.settings['xscale'])
        self.ax.set_yscale(self.settings['yscale'])
        
        # Set the tick label size
        self.ax.xaxis.set_tick_params(labelsize=PLOTFONTSIZE)
        self.ax.yaxis.set_tick_params(labelsize=PLOTFONTSIZE)
                
        self.plot_baseline(xlim)
        
        if self.canvas is None:
            plt.draw()
    
    def plot_baseline(self, xrange):
        """Plot the reduced chisq baseline
        """
        self.ax.plot(xrange, (1.0, 1.0), '--', c='#8b0000')
    
