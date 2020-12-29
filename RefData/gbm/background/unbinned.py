import numpy as np
from scipy.interpolate import interp1d

class NaivePoisson(object):
    """ A class to estimate the background of un-binned data using the naive Poisson 
    maximum likelihood.
    
    This method is approximately equivalent to sliding a window of fixed length 
    through un-binned data and calculating the Poisson maximum likelihood for the rate.
    The rate estimate is applied to the center of the sliding window, therefore, the 
    amount of data equivalent to half of the sliding window at the beginning and half 
    of the window at the end of the data is a constant.

    This naive approach assumes there is either no 'strong' signal in the data, or the
    presence of a 'weaker' signal has a duration much less than the window width of the
    sliding window.    
    """
    
    def __init__(self, window_width):
        """
        Parameters:
        --------------        
        window_width: float
        The width of the sliding window in seconds. NOTE: if the range of the data is 
        shorter than window_width, then the window_width will automatically be shortened 
        to the range of the data.
        """
        self._min_dt = 2e-6
        self.window_width = window_width
        self._actual_widths = None
        self._rates = None
        self._times = None

    def fit(self, events):
        """Estimates the background of un-binned data using the naive Poisson maximum 
        likelihood
       
        Parameters:
        -----------
        events:  np.array
            an array of consecutive event times
    
        Returns:
        -----------
        brates:  np.array
            The Poisson rate calculated at the time of each event.
        uncert: np.array
            The Poisson rate uncertainty based on the duration of the sliding window
            at the time of each event.
        """
        num_events = len(events)
        if num_events == 0:
            return np.array([])
    
        # get extent of window at the start of the data
        end_idx = np.where(events <= events[0] + self.window_width)[0][-1]
        # from the middle of the sliding window
        mid_idx = np.where(events <= events[0] + self.window_width/2.0)[0][-1]
    
        # Instead of employing an actual sliding window, set the window width, but require 
        # the window to always contain the same number of counts.  This effectively slides
        # the window from event to event, but allows the width to vary a bit. This is as
        # accurate as a fixed width window, but much much more efficient. 
    
        # number of 'bins' containing a count (= number of counts in window)
        full_bins = float(end_idx+1)
        # number of 'bins' containing a count not in the window
        num_back_bins = num_events - end_idx + 1
        array_back = np.arange(num_back_bins)
    
        # number of empty 'bins' in the window, shifted along the entire event list
        empty_bins = np.ceil((events[end_idx+array_back-1] - 
                              events[array_back])/self._min_dt)-full_bins
    
        # actual window widths
        dts = (full_bins+empty_bins)*self._min_dt
    
        # Poisson maximum likelihood rate (just the fraction of total bins in the window that
        # contain a count!)
        back_rates = full_bins / dts
    
        # insert the background rates for the full window width into the time arrays    
        brates = np.zeros_like(events)
        brates[mid_idx+array_back-1] = back_rates    
        dt = np.zeros_like(events)
        dt[mid_idx+array_back-1] = dts
    
        # Now we want to estimate the rate at the ends. This means that the sliding window
        # will truncate at either end until it's the width of one event.  The main concern
        # is to calculate the width of the truncated window correctly so that the uncertainty
        # in the rate will properly increase as the window width decreases
    
        # indices into the pre- and post-full-window edges of the data
        pre_full_bins = np.arange(1, end_idx)
        post_full_bins = pre_full_bins[::-1]
    
        # the rate and window width at the start of the data
        pre_empty_bins = np.ceil((events[end_idx-post_full_bins+1]-events[0])/self._min_dt) \
                         - pre_full_bins
        pre_dts = (pre_full_bins + pre_empty_bins) * self._min_dt
        pre_back_rates = pre_full_bins / pre_dts
    
        # the rate and window width at the end of the data
        post_empty_bins = np.ceil((events[-1]-events[end_idx+array_back[-1]-post_full_bins])/\
                                   self._min_dt) - post_full_bins
        post_dts = (post_full_bins + post_empty_bins) * self._min_dt
        post_back_rates = post_full_bins / post_dts
    
        # insert the pre- and post-full-window rates into the times array
        brates[0:mid_idx-1] = pre_back_rates[0:mid_idx-1]
        dt[0:mid_idx-1] = pre_dts[0:mid_idx-1]
        n = len(brates[mid_idx+array_back[-1]:])
        # only do this if our data range is longer than our sliding window width
        if n > 0:
            brates[mid_idx+array_back[-1]:] = post_back_rates[-n:]
            dt[mid_idx+array_back[-1]:] = post_dts[-n:]
    
        # copy the second (next to last) rates and widths into the first (last) elements
        # otherwise the window widths are identically 0.
        brates[0] = brates[1]
        dt[0] = dt[1]
        brates[-1] = brates[-2]
        dt[-1] = dt[-2]
        
        self._actual_widths = dt
        self._rates = brates
        self._times = events
        uncert = np.sqrt(brates/dt)
        return brates, uncert
    
    def interpolate(self, times, *args):
        """Interpolate the background at the given times
        
        Parameters:
        --------------
        times: np.array
            The times at which to interpolate

        Returns:
        --------------
        np.array
            The interpolated background rate at each time
        np.array
            The background uncertainty at each time
        """
        rates_interp = interp1d(self._times, self._rates, fill_value='extrapolate')
        width_interp = interp1d(self._times, self._actual_widths, fill_value='extrapolate')
        rates = rates_interp(times)
        widths = width_interp(times)
        uncert = np.sqrt(rates/widths)
        return rates, uncert
        
        

