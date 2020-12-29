import numpy as np
from warnings import warn
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt


class RateHisto(object):
    """Class containing rate-based histogram data, e.g. lightcurves, spectra

    Attributes:
    -----------
    values: np.array
        count values in a bin, size n
    num_bins: int
        The number of bins
    bounds: (np.array, np.array)
        tuple holding lower and upper bounds of the bins, size (n,n)
    edges: np.array
        array of bin edges, size (n+1)
    rate: np.array
        rate per bin, size n
    exposure: np.array
        The exposure of each bin
    uncertainty: np.array
        The uncertainty in the rate for each bin
    
    Public Methods:
    ---------------
    bin_centers:
        Return the center of each bin in the RateHisto
    bin_widths:
        Return the width of each bin in the RateHisto
    slice:
        Return a slice of the RateHisto based on a time or energy range
    rebin:
        Rebin the RateHisto based on a provided binning algorithm and return
        a new RateHisto
        
    Static Methods:
    ---------------
    sum_time:
        Sum a list of time rate histograms
    sum_energy:
        Sum a list of energy rate histograms
    merge:
        Merge a list of rate histrograms
    plot:
        Plot the data contained in the RateHisto object.  This is currently a diagnostic
        plot and does not contain full functionality
    """
    def __init__(self, values, bounds, exposure=None):
        """Initializes the attributes and calculates the rate per bin

        Parameters:
        -----------
        values: np.array
            counts per bin
        bounds: (np.array, np.array)
            tuple of lower and higher bin edges
        exposure: float or np.array; optional
            The exposure of the bin.  If omitted, the bin width is used.
        """
        assert values is not None, "The input must be greater than length 0"
        assert values.size > 0, "The input must be greater than length 0"
        if np.any(values < 0.0):
            warn('Negative rate values present')

        self.values = values
        self.num_bins = values.size
        self.bounds = bounds
        self.exposure = self._calc_exposure(exposure)
        self.rate = RateHisto.calc_rate(self.values, self.exposure)
        self.uncertainty = RateHisto.calc_uncertainty(self.values, self.exposure)
        self.edges = self.generate_edges(bounds)

    def _calc_exposure(self, exposure):
        """Sets the exposure from either None, float, or np.array

        Parameters:
        -----------
        exposure: None, float, or np.array
            The exposure

        Returns:
        -----------
        exposure: np.array
            The calculated or applied exposure
        """
        if exposure is None:
            exposure = self.bin_widths()
        if type(exposure) == np.ndarray:
            assert exposure.size == self.num_bins, "The exposure, when input as " \
                   "an array, must be the same length as the number of bins"
        elif type(exposure) == float:
            exposure = np.array([exposure]*self.num_bins)
        else:
            raise TypeError("Exposure must either be a float or an array")

        mask = (exposure < 0.0)
        if np.sum(mask) > 0:
            print("An exposure is negative.  Setting to zero.")
            exposure[mask] = 0.0

        return exposure

    @staticmethod
    def generate_edges(bounds):
        """Calculate bin edges from input bin bounds

        Parameters:
        -----------
        bounds: (np.array, np.array)
            tuple of lower and higher bin edges

        Returns:
        -----------
        edges: np.array
            array of bin edges
        """
        edges = np.concatenate((bounds[0], bounds[1][-1:]))
        return edges

    @staticmethod
    def generate_bounds(edges):
        """Calculate bounds from input bin edges

        Parameters:
        -----------
        edges: np.array
            array of bin edges

        Returns:
        -----------
        bounds: (np.array, np.array)
            tuple of lower and higher bin edges

        """
        bounds = (edges[:-1], edges[1:])
        return bounds

    @staticmethod
    def calc_rate(values, exposure):
        """Calculate rate in each bin

        The rate is caculated as the counts divided by the exposure.

        Parameters:
        -----------
        values: np.array
            Counts per bin
        exposure: np.array
            The exposure
        Returns:
        -----------
        rate: np.array
            count rate per bin
        """
        rate = np.zeros_like(values, dtype=float)
        mask = (exposure > 0.0)
        rate[mask] = values[mask].astype(float)/exposure[mask]
        return rate

    @staticmethod
    def calc_uncertainty(values, exposure):
        """Calculate uncertainty in the rate in each bin

        Parameters:
        -----------
        values: np.array
            Counts per bin
        exposure: np.array
            The exposure
        Returns:
        -----------
        uncertainty: np.array
            The count rate uncertainty for each bin
        """
        uncertainty = np.zeros_like(values, dtype=float)
        mask = (exposure > 0.0)
        uncertainty[mask] = np.sqrt(values[mask].astype(float))/exposure[mask]
        return uncertainty

    def bin_centers(self):
        """Calculate bin centers

        Returns:
        -----------
        bin_centers: np.array
            The center of each bin
        """
        bin_centers = (self.bounds[0]+self.bounds[1])/2.0
        return bin_centers
    
    def bin_widths(self):
        """Calculate bin widths

        Returns:
        -----------
        bin_widths: np.array
            The width of each bin
        """
        return self.bounds[1]-self.bounds[0]
    
    def slice(self, param_range):
        """Slice the current RateHisto based on a given range

        Parameters:
        -----------
        range: tuple (low, high) 
            The range of time or energy over which to slice
        
        Returns:
        -----------
        sliced_rates: RateHisto
            The sliced rate histogram
        """
        mask = (self.bounds[0] < param_range[1]) & (self.bounds[1] > param_range[0])
        sliced_values = self.values[mask]
        # if slice results in no data
        if sliced_values.size == 0:
            return None
        sliced_bounds = (self.bounds[0][mask], self.bounds[1][mask])
        sliced_exposure = self.exposure[mask]
        sliced_rates = RateHisto(sliced_values, sliced_bounds, exposure=sliced_exposure)
        return sliced_rates

    def rebin(self, method, *args, start=None, stop=None):
        """Rebin the current RateHisto

        Parameters:
        -----------
        method: function
            The binning method to implement.  The function must take in arrays that
            represents the counts in each current bin, the exposure in each current bin, 
            and the current bin edges in that order.  The function must return the new
            counts in each bin, the new exposure in each bin, and the new bin edges in 
            that order.
        *args: 
            Other optional arguments to be passed to the binning function
        start: float, optional
            The start of the data range to be rebinned. The default is to start at the 
            beginning of the RateHisto segment.
        stop: float, optional
            The end of the data range to be rebinned. The default is to stop at
            the end of the RateHisto segment.
        
        Returns:
        -----------
        new_rates: RateHisto
            A new re-binned RateHisto
        """
        # set the start and stop of the rebinning segment
        if (start is None):
            start = self.edges[0]
        if (stop is None):
            stop = self.edges[-1]
        if (start < self.edges[0]):
            start = self.edges[0]
        if (stop > self.edges[-1]):
            stop = self.edges[-1]
        # split the RateHisto into pieces so that we only rebin the piece that
        # needs to be rebinned
        pre = None
        post = None
        if (start != self.edges[0]) & (stop == self.edges[-1]):
            pre = self.slice((self.edges[0], start))
            rate = self.slice((start, self.edges[-1]))
        elif (start == self.edges[0]) & (stop != self.edges[-1]):
            rate = self.slice((self.edges[0], stop))
            post = self.slice((stop, self.edges[-1]))
        elif (start != self.edges[0]) & (stop != self.edges[-1]):
            pre = self.slice((self.edges[0], start))
            rate = self.slice((start, stop))
            post = self.slice((stop, self.edges[-1]))
        else:
            rate = self
        
        # perform the rebinning and create a new RateHisto with the rebinned rates
        if rate is not None:
            new_counts, new_exposure, new_edges = method(rate.values, rate.exposure, 
                                                         rate.edges, *args)
            new_bounds = RateHisto.generate_bounds(new_edges)
            new_rates = RateHisto(new_counts, new_bounds, exposure=new_exposure)
        else:
            new_rates = None
        
        # now merge the split RateHisto back together again
        rates_to_merge = [i for i in (pre, new_rates, post) if i is not None]
        new_rates = RateHisto.merge(rates_to_merge)
        
        return new_rates

    @staticmethod
    def sum_time(ratehistos):
        """Sum count rate histograms with different exposures

        Parameters:
        -----------
        ratehistos: list of RateHisto
            The rate histograms to sum
        
        Returns:
        -----------
        sum_rates: RateHisto
            The summed rate histogram
        """
        values = np.zeros(ratehistos[0].num_bins)
        exposure = np.zeros(ratehistos[0].num_bins)
        for ratehisto in ratehistos:
            assert np.all(ratehisto.edges == ratehistos[0].edges), "The two histograms " \
                   "must have the same support"
            #assert np.all(ratehisto.exposure == ratehistos[0].exposure), "The two " \
            #       "histograms must have the same exposure" 
            values += ratehisto.values
            binwidths = ratehisto.bounds[1]-ratehisto.bounds[0]
            exposure += ratehisto.exposure/binwidths

        exposure *= binwidths
        sum_rates = RateHisto(values, ratehistos[0].bounds, 
                              exposure=ratehistos[0].exposure)
        return sum_rates

    @staticmethod
    def sum_energy(ratehistos):
        """Sum count rate histograms that have the same exposure

        Parameters:
        -----------
        ratehistos: list of RateHisto
            The rate histograms to sum
        
        Returns:
        -----------
        sum_rates: RateHisto
            The summed rate histogram
        """
        values = np.zeros(ratehistos[0].num_bins)
        exposure = ratehistos[0].exposure
        for ratehisto in ratehistos:
            assert np.all(ratehisto.edges == ratehistos[0].edges), "The two histograms " \
                   "must have the same support"

            values += ratehisto.values

        sum_rates = RateHisto(values, ratehistos[0].bounds, exposure=exposure)
        return sum_rates

    @staticmethod
    def merge(ratehistos):
        """Merge rate histograms

        Parameters:
        -----------
        ratehistos: list of RateHisto
            The rate histograms to concatenate
        
        Returns:
        -----------
        merged_rates: RateHisto
            The merged rate histogram
        """
        num = len(ratehistos)
        # sort by start times
        tstarts = np.array([ratehisto.edges[0] for ratehisto in ratehistos])
        idx = np.argsort(tstarts)

        # concatenate the RateHistos in order
        counts = ratehistos[idx[0]].values
        bounds = list(ratehistos[idx[0]].bounds)
        exposure = ratehistos[idx[0]].exposure
        for i in range(1,num):
            bin_starts = ratehistos[idx[i]].bounds[0]
            # make sure there is no overlap - remove if there is
            mask = (bin_starts >= bounds[1][-1])
            counts = np.concatenate((counts, ratehistos[idx[i]].values[mask]))
            bounds[0] = np.concatenate((bounds[0], 
                                        ratehistos[idx[i]].bounds[0][mask]))
            bounds[1] = np.concatenate((bounds[1], 
                                        ratehistos[idx[i]].bounds[1][mask]))
            exposure = np.concatenate((exposure, 
                                       ratehistos[idx[i]].exposure[mask]))
        # new RateHisto
        merged_rates = RateHisto(counts, bounds, exposure=exposure)
        return merged_rates

    @staticmethod
    def plot(ratehistos, uncertainty=False, credible=False, log=False, color='blue', **kwargs):
        """Basic plot method - for diagnostics now, may be expanded in the future

        Parameters:
        -----------
        ratehistos: list
            A list of RateHisto objects to plot
        uncertainty: bool, optional
            If True, then plot error bars. Defaults is False
        credible: bool, optional
            If True, then plot the credible region.
        log: bool, optional
            If True, then plot axes in log scale.  Default is False
        """
        for rates in ratehistos:
            plt.step(rates.bounds[0], rates.rate, where='post', c=color, **kwargs)
        if log:
            plt.gca().set_yscale('log')
            plt.gca().set_xscale('log')
        if uncertainty:
            for rates in ratehistos:
                bincent = rates.bin_centers()
                plt.errorbar(bincent, rates.rate, yerr=rates.uncertainty, 
                                   linewidth=0, capthick=0, elinewidth=0.5, color=color)
        elif credible:
            for rates in ratehistos:
                t = np.array(rates.bounds).T.flatten().tolist()
                r1 = [rates.rate-rates.uncertainty, rates.rate-rates.uncertainty]
                r1 = np.array(r1).T.flatten().tolist()
                r2 = [rates.rate+rates.uncertainty, rates.rate+rates.uncertainty]
                r2 = np.array(r2).T.flatten().tolist()
                plt.fill_between(t, r1, r2, color=color, alpha=0.9)


class EventList(object):
    """Class containing a list of discrete events, each with a set of attributes

    Attributes:
    -----------
    size: int
        The number of contained events
    
    Public Methods:
    slice:
        Slice the EventList based on a range of attribute values
    bin:
        Bin the EventList data based on a provided binning algorithm
    
    Static Methods:
    ---------------
    merge:
        Merge a list of EventLists into a single EventList
    """
    def __init__(self, events):
        """Initializes the EventList

        Parameters:
        -----------
        events: np.recarray
            The events to be stored, where each field represents an attribute
        """
        assert events is not None, "events must be an array"
        self._events = events
        self.size = self._events.size
        self.exposure = 0.0
        if self.size > 0:
            self._calc_exposure()

    def _calc_exposure(self):
        """Calculates the exposure based on the deadtime definitions for GBM
        """
        tstart = np.min(self._events['TIME'])
        tstop = np.max(self._events['TIME'])
        mask = (self._events['PHA'] == 127)
        deadtime = np.sum(~mask)*2.6e-6 + np.sum(mask)*7.4e-6
        self.exposure = (tstop-tstart) - deadtime

    def slice(self, axis, range):
        """Slice the EventList based on an attribute and return a new EventList

        Parameters:
        -----------
        axis: str
            The name of the attribute over which to perform the slice
            
        range: tuple
            The slice range (lo, hi)
        
        Returns:
        -----------
        new_eventlist: EventList
            The sliced event list
        """
        if axis not in self._events.dtype.names:
            raise ValueError("Axis name {0} does not exist in EventList".format(axis))
        assert range[1] >= range[0], "The slice range must be in ascending order"

        mask = (self._events[axis] >= range[0]) & (self._events[axis] <= range[1])
        new_eventlist = EventList(self._events[mask])
        return new_eventlist

    def bin(self, method, *args, start=None, stop=None, **kwargs):
        """Bin the current EventList

        Parameters:
        -----------
        method: function
            The binning method to implement. The function must take in a single
            array representing the time of each event.  The function must 
            return the counts in each bin and the new bin edges in that order.
            
        *args: 
            Other optional arguments to be passed to the binning function
        
        **kwargs:
            Other optional arguments to be passed to the binning function
        
        Returns:
        -----------
        rates: RateHisto
            The binned data
        """
        # set the start and stop of the rebinning segment
        if (start is None):
            start = self._events['TIME'][0]
        if (stop is None):
            stop = self._events['TIME'][-1]
        
        mask = (self._events['TIME'] >= start) & (self._events['TIME'] <= stop)
        kwargs['tstart']=start
        kwargs['tstop']=stop
        counts, edges = method(self._events['TIME'][mask], *args, **kwargs)
        
        # calculate deadtime
        himask = (self._events['PHA'][mask] == 127)
        overflow_counts = np.histogram(self._events['TIME'][mask][himask], edges)[0]
        deadtime = counts*2.6e-6 + overflow_counts*7.4e-6

        # calculate exposure
        exposure = (edges[1:]-edges[:-1])-deadtime

        # create a RateHisto
        rates = RateHisto(counts, (edges[:-1], edges[1:]), exposure=exposure)

        return rates

    def plot(self):
        p = plt.scatter(self._events['TIME'], self._events['PHA'], s=1)
        plt.show()

    @staticmethod
    def merge(eventlists, sort_axis=None):
        """Merge a list of EventLists into a single EventList

        Parameters:
        -----------
        eventlists: list
            The EventLists to be merged
            
        sort_axis: str, optional
            The name of the axis to sort over. If omitted, no sorting is done. 
        
        Returns:
        -----------
        new_eventlist: EventList
            The merged event list
        """
        # get the events recarray for each EventList and stack
        #new_events = np.hstack([eventlist._events for eventlist in eventlists])

        # Discovered an issue in the (supposed) way that concatenation is done 
        # for FITS_rec arrays.  FITS_rec inherits from numpy.recarray, but adds
        # some extra stuff.  Apparently when FITS_rec objects are concatenated, 
        # they are converted to pure numpy.recarray objects, which removes the
        # extra stuff in FITS_rec.  This is not a problem here, however, 
        # FITS_rec apparently applies the TZERO scaling *before* converting to
        # recarrrays, and then throws away the TZERO scaling, meaning that we
        # can't convert it back.  This has major implications for continuous
        # data with no reference point, because GBM (incorrectly) sets TZERO
        # to TSTART for continuous data.  This results is an array of METs for 
        # the FITS_rec, but the resulting array after the concatenation is 
        # relative to TSTART, causing everything that uses this to fail for
        # non-trigger files.
        new_events = []
        for eventlist in eventlists:
            if eventlist.size == 0:
                continue
            ref_time = eventlist._events['TIME'][0]
            temp_events = np.array(eventlist._events)
            if temp_events['TIME'][0] != ref_time:
                offset = ref_time-temp_events['TIME'][0]
                temp_events['TIME'] += offset
            new_events.append(temp_events)
        new_events = np.hstack(new_events)
                
        if sort_axis is not None:
            if sort_axis not in new_events.dtype.names:
                raise ValueError("Axis name {0} does not exist in "
                                 "recarray".format(sort_axis))
            # sort by a recarray field
            idx = np.argsort(new_events[sort_axis])
            new_events = new_events[idx]
        # create new EventList
        new_eventlist = EventList(new_events)
        return new_eventlist
