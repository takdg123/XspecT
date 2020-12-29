import copy
import numpy as np
from gbm.data.containers import RateHisto, EventList
from gbm.data.phaii import Ctime, Cspec, TTE

class Selection(object):
    """Class which contains and links together the time and energy selections
    Attributes:
    -----------
    lightcurve: list
        A list of RateHisto lightcurve objects
    spectrum: list
        A list of RateHisto spectrum objects
    livetime: float
        The total livetime of the selection
    time_bound: (float, float)
        The convex hull of the time selection
    energy_bound: (float, float)
        The convex hull of the energy selection
    
    Public Methods:
    ---------------
    set_time_selections:
        Set the time selection
    set_time_binning:
        Set the default time binning method of the selection
    set_energy_selections:
        Set the energy selection
    set_energy_binning:
        Set the default energy binning method of the selection
    update_selection:
        Update the object with the current time and energy selections
    full_spectrum_for_each_bin:
        The full count spectrum for each time bin in the selection
    full_spectrum_integrated:
        The full count spectrum integrated in time
    slice:
        Return a new Selection that is a slice of the current selection   
    """    
    def __init__(self, data_obj, time_selections=None, energy_selections=None):
        """Initializes a Selection object with some data and optionally the 
        time and energy selections

        Parameters:
        -----------
        data_obj: data.PHAII
            A data object (Ctime, Cspec, TTE)
        time_selections: tuple or list of tuples, optional
            If the selection is contiguous, this parameter should be a tuple
            of the form (lo, hi).  If the selection is composed of disjoint
            selections, this parameter should be a list of tuples
        energy_selection:  tuple or list of tuples, optional
            If the selection is contiguous, this parameter should be a tuple
            of the form (lo, hi).  If the selection is composed of disjoint
            selections, this parameter should be a list of tuples
        """
        self._data_obj = data_obj
        self.type = data_obj.type
        self.lightcurve = None
        self.single_channels = None
        self.spectrum = None
        self.livetime = None
        self.time_bound = None
        self.energy_bound = None
        self._time_binning = None
        self._energy_binning = None
        
        if time_selections is None:
            self.set_time_selections(data_obj.time_range)
        else:
            self.set_time_selections(time_selections)
        if energy_selections is None:
            self.set_energy_selections(data_obj.energy_range)
        else:
            self.set_energy_selections(energy_selections)

        self.update_selection()
    
    def __str__(self):
        return 'Time Range: {0}\nEnergy Range: {1}'.format(self.time_bound, self.energy_bound)
    
    def set_time_selections(self, time_selections):
        """Set the time selection(s).  Can be contiguous or disjoint.

        Parameters:
        -----------
        time_selections: tuple or list of tuples
            If the selection is contiguous, this parameter should be a tuple
            of the form (lo, hi).  If the selection is composed of disjoint
            selections, this parameter should be a list of tuples
        """
        time_selections = self._assert_selections(time_selections)
        self._time_selections = time_selections
        self.time_bound = (time_selections[0][0], time_selections[-1][1])
        
    def set_energy_selections(self, energy_selections):
        """Set the energy selection(s).  Can be contiguous or disjoint.

        Parameters:
        -----------
        energy_selections: tuple or list of tuples
            If the selection is contiguous, this parameter should be a tuple
            of the form (lo, hi).  If the selection is composed of disjoint
            selections, this parameter should be a list of tuples
        """
        energy_selections = self._assert_selections(energy_selections)
        self._energy_selections = energy_selections
        self.energy_bound = (energy_selections[0][0], energy_selections[-1][1])
        
    def _assert_selections(self, selections):
        """Check to ensure the selections are of the correct form.

        Parameters:
        -----------
        selections: tuple or list of tuples
            The selection(s) to check
        
        Returns:
        --------
        selections: list
        """
        if (all(isinstance(selection, list) for selection in selections)) | \
           (all(isinstance(selection, tuple) for selection in selections)):
            if any(len(selection) != 2 for selection in selections):
                raise ValueError('Each range in selections must be of the '
                                 'form (lo, hi)')
            else:
                return selections
        else:
            if len(selections) != 2:
                raise ValueError('Selections must either be a range of '
                                   'the form (lo, hi) or a list of ranges')
            else:
                return [selections]
        
    def set_time_binning(self, method, *args, start=None, stop=None, **kwargs):
        """Set the default time (re)binning for the selection

        Parameters:
        -----------
        method: function
            The binning method to use
        *args:
            Other optional arguments to be passed to the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the 
            beginning of the selection
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the selection.
        **kwargs:
            Other optional keywords to be passed to the binning function
        """
        if start is None:
            start = self._data_obj.time_range[0]
        if stop is None:
            stop = self._data_obj.time_range[1]
        if method is None:
            self._time_binning = None
        else:
            self._time_binning = {'method': method, 'args': args, 'start': start, 
                                  'stop': stop, 'kwargs': kwargs}
        self._single_channel_lightcurves()
        self._integrate_over_energy()
    
    def set_energy_binning(self, method, *args, start=None, stop=None, **kwargs):
        """Set the default energy (re)binning for the selection

        Parameters:
        -----------
        method: function
            The binning method to use
        *args:
            Other optional arguments to be passed to the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the 
            beginning of the selection.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the selection.
        **kwargs:
            Other optional keywords to be passed to the binning function
        """
        if start is None:
            start = self._data_obj.energy_range[0]
        if stop is None:
            stop = self._data_obj.energy_range[1]
        if method is None:
            self._energy_binning = None
        else:
            self._energy_binning = {'method': method, 'args': args, 'start': start, 
                                    'stop': stop, 'kwargs': kwargs}
        self._integrate_over_time()
    
    def _rebin_binned_time(self, rates):
        """Rebins the binned lightcurve
    
        Parameters:
        -----------
        rates: list of RateHisto
            The lightcurve to be rebinned
       
        Returns:
        -----------
        rates: list of RateHisto
            The rebinned lightcurve
        """
        # if using a bin index list (like in rmfit), then this needs to be
        # split over SAA gaps
        if self._time_binning['method'].__name__ == 'rebin_by_edge_index':
            bounds = [rate.bounds for rate in rates]
            split_indices = self._data_obj._split_bin_index_on_saa(bounds, 
                                                    *self._time_binning['args'])
            rates = [rates[i].rebin(self._time_binning['method'], 
                                    split_indices[i]) for i in range(len(rates))]
        else:
            rates = [rate.rebin(self._time_binning['method'], 
                                *self._time_binning['args'], 
                                start=self._time_binning['start'], 
                                stop=self._time_binning['stop']) for rate in rates]
        
        return rates

    def _rebin_unbinned_time(self, eventlist):
        """Rebins the unbinned lightcurve
    
        Parameters:
        -----------
        eventlist:
            The unbinned data to be rebinned
       
        Returns:
        -----------
        rates: list of RateHisto
            The rebinned lightcurve
        """
        if len(eventlist) != 1:
            print('Warning: More than one Eventlist (Selection._rebin_unbinned_time)')
                
        rates = eventlist[0].bin(self._time_binning['method'], 
                                 *self._time_binning['args'],
                                 start=self._time_binning['start'], 
                                 stop=self._time_binning['stop'],
                                 **self._time_binning['kwargs'])

        # calculate the TTE exposure -- must be calculated over all channels
        exposure = self._data_obj.calculate_exposure(rates.edges)
        rates.exposure = exposure
        # split the data over SAA traversal - no data
        counts, bounds, exposure = self._data_obj._split_on_saa(rates.values, \
                                                                rates.bounds, \
                                                                rates.exposure)
        num_intervals = len(counts)
        rates = []
        for i in range(num_intervals):
            rates.append(RateHisto(counts[i], bounds[i], exposure=exposure[i]))
        
        return rates
    
    def _merge_unbinned_rebins(self, lightcurve, new_rebin):
        """Merge the rebinned segment with the previous binned lightcurve
    
        Parameters:
        -----------
        lightcurve: list of RateHistos
            The lightcurve with the old binning
        new_rebin: list of RateHistos
            The new rebinned segment to be merged
       
        Returns:
        -----------
        new_rates: list of RateHistos
            The new merged lightcurve
        """
        numseg = len(lightcurve)
        new_rates = []
        # cycle over the number of lightcurve segments (may be split over SAA)
        for iseg in range(numseg):
            # cycle over the new rebinned segment (may also be split over SAA)
            for segment in new_rebin:
                start = segment.edges[0]
                stop = segment.edges[-1]
                # if new rebinned segment is not contained in this lightcurve
                # segment, then we don't have to do anything
                if (start > lightcurve[iseg].edges[-1]) | (stop < lightcurve[iseg].edges[0]):
                    new_rates.append(lightcurve[iseg])
                    continue
                
                # split up the lightcurve segment into pre- and post rebinned
                # subsegments
                pre = None
                post = None
                if (start != lightcurve[iseg].edges[0]) & \
                   (stop == lightcurve[iseg].edges[-1]):
                    pre = lightcurve[iseg].slice((lightcurve[iseg].edges[0], start))
                elif (start == lightcurve[iseg].edges[0]) & \
                 (stop != lightcurve[iseg].edges[-1]):
                    post = lightcurve[iseg].slice((stop, lightcurve[iseg].edges[-1]))
                elif (start != lightcurve[iseg].edges[0]) & \
                 (stop != lightcurve[iseg].edges[-1]):
                    pre = lightcurve[iseg].slice((lightcurve[iseg].edges[0], start))
                    post = lightcurve[iseg].slice((stop, lightcurve[iseg].edges[-1]))
                else:
                    pass
                
                # now we can merge
                rates_to_merge = [i for i in (pre, segment, post) if i is not None]
                new_rates.append(RateHisto.merge(rates_to_merge))
        
        return new_rates
        
    def _rebin_energy(self, rates):
        """Rebins the binned count spectrum
    
        Parameters:
        -----------
        rates: list of RateHisto
            The count spectrum to be rebinned
       
        Returns:
        -----------
        rates: list of RateHisto
            The rebinned count spectrum
        """
        # if using a bin index list (like in rmfit), then this needs to be
        # split/shifted for any energy selections
        if self._energy_binning['method'].__name__ == 'rebin_by_edge_index':
            bounds = self._energy_selections
            split_indices = self._data_obj._split_bin_index_energy_bounds(bounds, 
                                                  *self._energy_binning['args'])
            rates = [rates[i].rebin(self._energy_binning['method'], 
                                    split_indices[i]) for i in range(len(rates))]
        else:
            rates = [rate.rebin(self._energy_binning['method'], 
                                *self._energy_binning['args'], 
                                start=self._energy_binning['start'], 
                                stop=self._energy_binning['stop']) for rate in rates]
        return rates
    
    def _get_energy_channel_indices(self):
        """Returns energy channel indices in the energy selection
        """
        energy_channels = []
        for energy_range in self._energy_selections:
            chan_lo, chan_hi = self._data_obj._assert_energy_range(energy_range)
            channels = np.arange(chan_lo, chan_hi+1)
            energy_channels.extend(channels)
        return energy_channels
    
    def _get_energy_channel_bounds(self):
        """Returns energy channel bounds in the energy selection
        """
        energy_ranges = []
        for energy_range in self._energy_selections:
            chan_lo, chan_hi = self._data_obj._assert_energy_range(energy_range)
            channels = np.arange(chan_lo, chan_hi+1)
            energy_ranges.extend([self._data_obj._channel_range_to_energy((chan, chan)) \
                                  for chan in channels])
        return energy_ranges
        
    def _single_channel_lightcurves(self):
        """Creates a lightcurve for each energy channel in the selection
        """
        energy_ranges = self._get_energy_channel_bounds()
        
        # get the data for each channel
        ichan = 0
        lightcurves = []
        for energy_range in energy_ranges:
            segments = []
            for time_range in self._time_selections:
                    view = self._data_obj.view(time_range=time_range, 
                                               energy_range=energy_range, 
                                               integration_axis='energy')
                    # do the binning
                    if self._time_binning is not None:
                        if self.type == 'binned':
                            view = self._rebin_binned_time(view)
                        elif self.type == 'unbinned':
                            view = self._rebin_unbinned_time(view)
                            # merge the new binning with previous binning
                            if type(self.single_channels[ichan][0]).__name__ != 'EventList':
                                view = self._merge_unbinned_rebins(self.single_channels[ichan], 
                                                                   view)
                            
                    segments.extend(view)
            #segments = np.array(segments).T.tolist()
            lightcurves.append(segments)
            ichan += 1
        self.single_channels = lightcurves
        
    def _integrate_over_energy(self):
        """Integrate the lightcurve over energy channels in the selecion
        """
        if self.single_channels is None:
            raise ValueError("Single channel lightcurves have not been calculated")
        
        numsegments = len(self.single_channels[0])
        lightcurve = []
        for iseg in range(numsegments):
            rates = [one_channel[iseg] for one_channel in self.single_channels]
            if (self.type == 'unbinned') & (self._time_binning is None):
                lightcurve.append(EventList.merge(rates))
            else:
                lightcurve.append(RateHisto.sum_energy(rates))
        self.lightcurve = lightcurve
        self._calculate_livetime()
    
    def _integrate_over_time(self):
        """Creates a count spectrum by integrating over the time selection
        """        
        # This sums over any disjoint selections in time and then creates
        # a list of spectra, splitting over disjoint selections in energy
        spectrum = []
        for energy_range in self._energy_selections:
            segments = []
            for time_range in self._time_selections:
                view = self._data_obj.view(time_range=time_range, 
                                           energy_range=energy_range, 
                                           integration_axis='time')
                # do the binning
                if self._energy_binning is not None:
                    view = self._rebin_energy(view)
                segments.append(view)
            segments = np.array(segments).T.tolist()
            segments = [RateHisto.sum_time(segment) for segment in segments]
            spectrum.extend(segments)
        self.spectrum = spectrum
    
    def update_selection(self):
        """Force an update of the selection.  Recalculates the integrated 
        lightcurve, integrated spectrum, and the lightcurve in each channel
        """
        self._single_channel_lightcurves()
        
        self._integrate_over_energy()
        
        self._integrate_over_time()
        
    def _calculate_livetime(self):
        """Total time exposure of the selection
        """
        self.livetime = 0
        for segment in self.lightcurve:
            self.livetime += np.sum(segment.exposure)
    
    def full_spectrum_for_each_bin(self):
        """Return the full spectrum for each time history bin in the selection.
        This is needed particularly for the creation of PHAII files for XSPEC
       
        Returns:
        -----------
        spectra: list of RateHisto objects
            The spectrum for each time bin
        """
        # list of time bins
        time_bins = []
        for lc in self.lightcurve:
            time_bins.extend(np.array(lc.bounds).T.tolist())
        
        # get the data for each time bin
        spectra = []
        for time_bin in time_bins:
            view = self._data_obj.view(time_range=time_bin, 
                                       integration_axis='time')
            
            #spectra.append(RateHisto.merge(view))
            # specially made for writeXSPEC
            view = RateHisto.merge(view)
            exposure = view.exposure / view.bin_widths()
            newview = RateHisto(view.values, view.bounds, exposure=exposure)
            spectra.append(newview)
            
        return spectra
    
    def full_spectrum_integrated(self):
        """Return the full spectrum integrated over time.
        This is needed particularly for the creation of PHAII files for XSPEC
        Note: This spectrum contains zeros for channels outside the source selection
        
        Returns:
        -----------
        integrated_spectrum: RateHisto
            The full time-integrated spectrum
        """
        numchans = self._data_obj.num_chans
        counts = np.zeros(numchans)
        bounds = (self._data_obj.ebounds['E_MIN'], self._data_obj.ebounds['E_MAX'])
        binwidths = bounds[1]-bounds[0]
        #exposure = self.livetime*binwidths
        exposure = float(self.livetime) # specially made for writeXSPEC
        
        spectrum = RateHisto.merge(self.spectrum)
        idx = self._get_energy_channel_indices()
        counts[idx] = spectrum.values
        integrated_spectrum = RateHisto(counts, bounds, exposure=exposure)
        
        return integrated_spectrum
    
    def slice(self, time_selections=None, energy_selections=None):
        """Make a slice from the current selection
    
        Parameters:
        -----------
        time_selections: tuple or list of tuples, optional
            If the selection is contiguous, this parameter should be a tuple
            of the form (lo, hi).  If the selection is composed of disjoint
            selections, this parameter should be a list of tuples
        energy_selection:  tuple or list of tuples, optional
            If the selection is contiguous, this parameter should be a tuple
            of the form (lo, hi).  If the selection is composed of disjoint
            selections, this parameter should be a list of tuples
       
        Returns:
        -----------
        new_select: Selection
            A new selection object
        """
        # make a copy of this object
        new_select = copy.deepcopy(self)
                
        # update energy selection
        if energy_selections is not None:
            new_select.set_energy_selections(energy_selections)
            
            # only keep the single channel lightcurves that correspond to the
            # new energy selection
            orig_chans = np.array(self._get_energy_channel_bounds())
            new_chans = np.array(new_select._get_energy_channel_bounds())
            numchans = len(new_chans)
            idx = [np.where(orig_chans[:,0] == new_chans[i,0])[0][0] for i in 
                   range(numchans)]
            single_channels = [self.single_channels[i] for i in idx]
            new_select.single_channels = single_channels
        
        # update spectrum
        spectrum = RateHisto.merge(self.spectrum)
        new_select.spectrum = [spectrum.slice(select) for select in 
                               new_select._energy_selections]
            
            
        # update time selection
        if time_selections is not None:
            new_select.set_time_selections(time_selections)
        
        # update the single channel lightcurves
        new_single_channels = []
        for one_channel in new_select.single_channels:
            if (new_select.type == 'unbinned') & (self._time_binning is None):
                lightcurve = EventList.merge(one_channel)
                new_single_channels.append([lightcurve.slice('TIME', select) for 
                                            select in new_select._time_selections])
            else:
                lightcurve = RateHisto.merge(one_channel)
                new_single_channels.append([lightcurve.slice(select) for select
                                            in new_select._time_selections])
        new_select.single_channels = new_single_channels    
        
        # update the integrated lightcurve
        new_select._integrate_over_energy()
        new_select._integrate_over_time()           
        
        return new_select

    
    