from gbm.data.phaii import Cspec, Ctime, TTE, RateHisto, EventList
from gbm.selection import Selection
from gbm.lookup import RmfitLookup, GspecLookup
from .background.background import Background
from .background.binned import polynomial
from .binning.binned import rebin_by_edge_index, rebin_by_snr
from .binning.unbinned import bin_by_edges, bin_by_time
import numpy as np
import matplotlib.pyplot as plt

"""
All of the stuff in this comment is deprecated
#################################################

def test_cspec():
    # read file
    cspec = Cspec('./test/glg_cspec_b0_bn120415958_v00.pha')
    #read lookup
    lu = RmfitLookup('./test/glg_cspec_b0_bn120415958_v00.lu')
    time_ranges = lu.background_selection()
    time_ranges = [t.tolist() for t in time_ranges]
    # one energy range
    energy_range = lu.energy_selection()[0].tolist()
    # only want the xmin, xmax of the display
    view_range = lu.time_display_range[0:2]
    
    # background selection
    bkgd_select = Selection(cspec, time_selections=time_ranges, \
                            energy_selections=energy_range)

    # background fit
    bkgd = Background(bkgd_select)
    bkgd.fit(polynomial, lu.background_poly_order)
    
    # view selection
    view_select = Selection(cspec, time_selections=view_range, \
                            energy_selections=energy_range)
    # background interpolated over the lightcurve view and integrated
    # over energy channels       
    b, u = bkgd.integrate_energy(view_select)
    
    # the background spectrum integrated over the view selection
    bkgd_spectrum = bkgd.source_spectrum(view_select)

    # plot the lightcurve and the background fit
    times = view_select.lightcurve[0].bin_centers()
    plt.subplot(2,1,1)
    plt.fill_between(times, b+u, b-u, color='red', alpha=0.5, zorder=10)
    RateHisto.plot(view_select.lightcurve)
    plt.subplot(2,1,2)
    RateHisto.plot(view_select.spectrum, log=True, uncertainty=True)
    RateHisto.plot(bkgd_spectrum, log=True, color='red', alpha=0.9, credible=True)
    plt.show()

def test_ctime():
    # read file
    ctime = Ctime('./test/glg_ctime_nb_bn120415958_v00.pha')
    # read lookup
    lu = RmfitLookup('./test/glg_ctime_nb_bn120415958_v00.lu')
    time_ranges = lu.background_selection()
    time_ranges = [t.tolist() for t in time_ranges]
    # one energy range
    energy_range = lu.energy_selection()[0].tolist()
    # only want xmin, xmax of the display
    view_range = lu.time_display_range[0:2]
    # the time history is rebinned, get new edges
    time_bin_indices = lu.rebinned_time_edges()
    
    # set the default time binning for the data
    ctime.set_time_binning(rebin_by_edge_index, time_bin_indices)
    
    # background selection
    bkgd_select = Selection(ctime, time_selections=time_ranges, \
                            energy_selections=energy_range)
    # background fit
    bkgd = Background(bkgd_select)
    bkgd.fit(polynomial, lu.background_poly_order)
    
    # view selection
    view_select = Selection(ctime, time_selections=view_range, \
                            energy_selections=energy_range)
        
    # background interpolated over the lightcurve view and integrated
    # over energy channels       
    b, u = bkgd.integrate_energy(view_select)
        
    # the background spectrum integrated over the view selection
    bkgd_spectrum = bkgd.source_spectrum(view_select)

    # plot the lightcurve and the background fit
    times = []
    for lc in view_select.lightcurve:
        times.extend(lc.bin_centers())
    plt.subplot(2,1,1)
    plt.fill_between(times, b+u, b-u, color='red', alpha=0.5, zorder=10)
    RateHisto.plot(view_select.lightcurve, uncertainty=True)
    plt.subplot(2,1,2)
    RateHisto.plot(view_select.spectrum, log=True, uncertainty=True)
    RateHisto.plot(bkgd_spectrum, log=True, color='red', alpha=0.9, credible=True)
    plt.show()

def test_tte():
    # read file
    tte = TTE('./test/glg_tte_n9_bn090131090_v00.fit')
    # read lookup
    lu = RmfitLookup('./test/glg_tte_n9_bn090131090_v00.lu', \
                     ti_file='./test/glg_tte_n9_bn090131090_v00.ti')
    time_ranges = lu.background_selection()
    time_ranges = [t.tolist() for t in time_ranges]
    # one energy range
    energy_range = lu.energy_selection()[0].tolist()
    # only want xmin, xmax of the display
    view_range = lu.time_display_range[0:2]
    # the TTE time bin edges
    time_edges = lu.tte_edges()
    # the source selection
    source_range = lu.time_selection()[0].tolist()
    
                              
    # set the binning method for the data
    tte.set_time_binning(bin_by_edges, time_edges)
    
    # background selection
    bkgd_select = Selection(tte, time_selections=time_ranges, \
                            energy_selections=energy_range)
    # background fit
    bkgd = Background(bkgd_select)
    bkgd.fit(polynomial, lu.background_poly_order)
    
    # view selection
    view_select = Selection(tte, time_selections=view_range, \
                            energy_selections=energy_range)    
    
    # background interpolated over the lightcurve view and integrated
    # over energy channels       
    b, u = bkgd.integrate_energy(view_select)

    # source selection
    source_select = Selection(tte, time_selections=source_range, \
                              energy_selections=energy_range)
    b2, u2 = bkgd.integrate_energy(source_select)
    
    # the background spectrum integrated over the source selection
    bkgd_spectrum = bkgd.source_spectrum(source_select)
    
    # some ganky stuff to get the source selection to plot
    merged_source = RateHisto.merge(source_select.lightcurve)
    source_bounds = [merged_source.bounds[0].T, merged_source.bounds[1].T]
    source_bounds = np.array(source_bounds).T.flatten().tolist()
    source_rates = [merged_source.rate, merged_source.rate]
    source_rates = np.array(source_rates).T.flatten().tolist()
    b2 = np.array([b2,b2]).T.flatten().tolist()

    # plot the lightcurve and the background fit
    times = []
    for lc in view_select.lightcurve:
        times.extend(lc.bin_centers())
    plt.subplot(2,1,1)
    plt.fill_between(times, b+u, b-u, color='red', alpha=0.9, zorder=10)
    plt.fill_between(source_bounds, source_rates, b2, color='purple', alpha=0.5)
    RateHisto.plot(view_select.lightcurve, uncertainty=True)
    plt.subplot(2,1,2)
    RateHisto.plot(source_select.spectrum, log=True, uncertainty=True, zorder=10)
    RateHisto.plot(bkgd_spectrum, log=True, color='red', alpha=0.9, credible=True)
    plt.show()
"""

def test_gspec():
    
    # initialize GSpec
    from .gspec import GspecManager
    gspec = GspecManager()
    
    # read a rmfit lookup file
    lu = RmfitLookup('./gbm/test/glg_tte_n9_bn090131090_v00.lu', \
                     ti_file='./gbm/test/glg_tte_n9_bn090131090_v00.ti')
    
    # and load the data and lookup into gspec
    gspec.add_data('./gbm/test/glg_tte_n9_bn090131090_v00.fit', lookup=lu)
    
    # add a new detector to gspec.  this one does not have a lookup
    gspec.add_data('./gbm/test/glg_tte_na_bn090131090_v00.fit') 
    
    # import the time binning from n9 to na
    gspec.import_time_binning('n9', ['na'])
    
    # import the background selections and fitter from n9 to na
    gspec.import_background('n9', ['na'])  
    
    # import the source selection from n9 to na
    gspec.import_source_selection('n9', ['na'])
    
    # import the display ranges from n9 to na
    gspec.import_time_display_view('n9', ['na'])
    gspec.import_energy_display_view('n9', ['na'])
    
    # refresh na with all of the new specifications
    gspec.refresh('na')
    
    # display the detectors
    gspec.display_detector('n9')
    gspec.display_detector('na')
