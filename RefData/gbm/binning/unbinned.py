import numpy as np

def bin_by_time(times, dt, tstart=None, tstop=None, time_ref=None):
    """Bins unbinned data to a specified temporal bin width.
        
    Parameters:
    -----------
    times: np.array
        The time of each event
    dt: float
        The requested temporal bin width in seconds
    tstart: float, optional
        The first bin edge time.  Will use the first event time if omitted.
    tstop: float, optional
        The last bin edge time.  Will use the last event time if omitted.
    time_ref: float, optional
        The reference time at which the binning will be based.  If the set, 
        The binning will proceed starting at time_ref and moving forward in time
        as well as starting at time_red and moving backward in time.  If not set,
        the binning will start at the beginning of the data.
       
    Returns:
    -----------
    np.array
        The counts in each bin
    np.array
        The time edges of the binned data
    """
    assert dt > 0.0, "Requested bin width must be > 0.0 s"
    if time_ref is not None:
        time_ref = float(time_ref)
    
    if tstart is None:
        tstart = np.min(times)
    if tstop is None:
        tstop = np.max(times)
    
    # if we are using a reference time
    if time_ref is not None:
        pre_edges = np.arange(time_ref, tstart-dt, -dt)[::-1]
        post_edges = np.arange(time_ref, tstop+dt, dt)
        edges = np.concatenate((pre_edges[:-1], post_edges))
    else:
        edges = np.arange(tstart, tstop+dt, dt)
    
    # bin the data
    counts = np.histogram(times, edges)[0]

    return counts, edges

def combine_into_one(times, tstart, tstop, **kwargs):
    """Bins unbinned data to a single bin.
        
    Parameters:
    -----------
    times: np.array
        The time of each event
    tstart: float
        The starting time
    tstop: float
        The stopping time
       
    Returns:
    -----------
    np.array
        The counts in each bin
    np.array
        The time edges of the binned data
    """
    assert tstart < tstop, "The time range must be set in " \
                                          "ascending order"
    
    # if no events, then the we have no counts
    if times.size == 0:
        return np.array([0]), np.array((tstart, tstop))
        
    if (np.min(times) > tstop) | (np.max(times) < tstart):
        raise ValueError("Requested time range is outside data range")
    
    # bin the data
    time_range = (tstart, tstop)
    counts = np.histogram(times, time_range)[0]
    return counts, np.array(time_range)

def combine_by_factor(times, old_edges, bin_factor, tstart=None, tstop=None):
    """Bins data to a multiple factor of bins given a set of bin edges
    
    Parameters:
    -----------
    times: np.array
        The time of each event
    old_edges: np.array
        The edges to be combined
    bin_factor: int
        The number of bins to be combined
    tstart: float, optional
        The first bin edge time.  Will use the first edge time if omitted.
    tstop: float, optional
        The last bin edge time.  Will use the last edge time if omitted.
       
    Returns:
    -----------
    np.array
        The counts in each bin
    np.array
        The time edges of the binned data
    """
    
    assert bin_factor >= 1, "bin_factor must be a positive integer"
    bin_factor = int(bin_factor)
    
    # mask the old edges for the desired data range
    if tstart is None:
        tstart = old_edges[0]
    if tstop is None:
        tstop = old_edges[-1]
    old_edges = old_edges[old_edges >= tstart]
    old_edges = old_edges[old_edges <= tstop]
    
    # create the new edges and counts within the bins
    new_edges = old_edges[::bin_factor]
    counts = np.histogram(times, new_edges)[0]
    return counts, new_edges

def bin_by_snr(times, back_rates, snr, **kwargs):
    """Bins unbinned data by SNR
        
    Parameters:
    -----------
    times: np.array
        The time of each event
    back_rates: np.array
        The background rate at the time of each event
    snr: float
        The signal-to-noise ratio
       
    Returns:
    -----------
    np.array
        The counts in each bin
    np.array
        The time edges of the binned data
    """
    # get the background rates and differential counts at each event time
    #back_rates, _ = background.interpolate(times)
    back_counts = np.zeros_like(back_rates)
    back_counts[:-1] = back_rates[:-1]*(times[1:]-times[:-1])
    back_counts[-1] = back_counts[-2]
    
    num_events = len(times)
    istart = 0
    edges = [istart]
    while (True):
        # cumulative sum of the counts, background counts, and snr
        countscum = np.arange(1, num_events+1-istart)
        backgroundscum = np.cumsum(back_counts[istart:])
        snrcum = (countscum-backgroundscum)/np.sqrt(backgroundscum)
        # determine where to make the cut
        below_thresh = np.sum(snrcum <= snr)
        if below_thresh == 0:
            iend = istart
        else:
            iend = istart + below_thresh-1
        edges.append(iend)
        if (iend >= num_events-1):
            break
        istart = iend+1
    
    # get the finalized edges and counts
    edges.append(num_events-1)
    edges = times[np.array(edges)]
    counts = np.histogram(times, edges)[0]
    return counts, edges

def time_to_spill(times, threshold, **kwargs):
    """Time-to-Spill Binning for an event list
    Bins an event list by accumulating counts until the set threshold is reached, 
    and then creating a new bin.
    
    Parameters
    -----------
    times: np.array
        a list of event times
    threshold: float
        the count threshold for the histogram
    
    Returns:
    -----------
    np.array
        The counts in each bin
    np.array
        The time edges of the binned data
    """
    threshold = int(threshold)
    assert threshold > 0, "Threshold must be positive"
    # this is easy: take only the nth time (set by threshold) as the bin edge
    edges = times[::threshold]
    counts = np.histogram(times, edges)[0]
    return counts, edges

def bin_by_edges(times, time_edges, **kwargs):
    """Bins unbinned data by pre-defined edges
        
    Parameters:
    -----------
    times: np.array
        The time of each event
    time_edges: np.array
        The pre-defined time edges
       
    Returns:
    -----------
    np.array
        The counts in each bin
    np.array
        The time edges of the binned data
    """
    counts = np.histogram(times, time_edges)[0]
    return counts, time_edges