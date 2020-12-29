import numpy as np

def combine_by_factor(counts, exposure, old_edges, bin_factor):
    """Rebins binned data to a multiple factor
    
    Parameters:
    -----------
    counts: np.array
        The counts in each bin
    exposure: np.array
        The exposure of each bin
    old_edges: np.array
        The time edges of each bin
    bin_factor: int
        The number of bins to be combined
       
    Returns:
    -----------
    np.array
        The counts in each new bin
    np.array
        The exposure of each new bin
    np.array
        The new time edges of the rebinned data
    """
    assert bin_factor >= 1, "bin_factor must be a positive integer"
    bin_factor = int(bin_factor)
    new_edges = old_edges[::bin_factor]
    # make sure the number of bins is a multiple of the bin_factor
    mod = counts.shape[0] % bin_factor
    # if the number of bins is not a multiple of the bin_factor, 
    # then remove the modulo bins
    if mod != 0:
        counts = counts[:-mod]
        exposure = exposure[:-mod]
    # reshape and combine counts and exposure
    new_counts = np.sum(counts.reshape(-1,bin_factor), axis=1)
    new_exposure = np.sum(exposure.reshape(-1,bin_factor), axis=1)
    return new_counts, new_exposure, new_edges

def rebin_by_time(counts, exposure, old_edges, dt):
    """Rebins binned data to a specified temporal bin width.
    
    If the requested bin width is smaller than some of the original bin widths, those bins
    will be left as is.
    
    If the requested bin width is an exact factor of all the current bin widths, the
    resulting bin width will be exactly as requested.  If the requested bin width is not
    an exact factor, then the resulting bin width will be approximately the requested bin
    width without exceeding the requested bin width.
    
    Parameters:
    -----------
    counts: np.array
        The counts in each bin
    exposure: np.array
        The exposure of each bin
    old_edges: np.array
        The time edges of each bin
    dt: float
        The requested temporal bin width in seconds
       
    Returns:
    -----------
    np.array
        The counts in each new bin
    np.array
        The exposure of each new bin
    np.array
        The new time edges of the rebinned data
    """
    assert dt > 0.0, "Requested bin width must be > 0.0 s"
    
    num_old_bins = counts.size
    dts = old_edges[1:]-old_edges[0:-1]
    
    # if the the requested bin width is a factor of the current bin width for all bins, 
    # call combine_by_factor.  this is an easier task
    if np.sum((dt % dts) == 0.0) == num_old_bins:
        return combine_by_factor(counts, exposure, old_edges, int(dt/dts[0]))
    
    #print('Temporal factor for rebinning is not a perfect multiple of all bin widths.')
    #print('Bin widths will be approximate to and will not exceed the temporal factor.')
    # cycle through current bins and add up bins until we have approximately reached,
    # but not exceeded, the requested binwidth
    dtcum = np.cumsum(dts)
    new_edges = [0]
    istart = 0
    while True:
        dtscum = np.cumsum(dts[istart:])
        iend = istart + np.sum(dtscum <= dt)
        if iend < istart:
            iend = istart
        new_edges.append(iend)
        if iend >= num_old_bins-1:
            break
        istart = iend+1
    
    # create the new rebinned arrays
    new_edges = np.array(new_edges)
    bounds_idx = np.array((new_edges[0:-1], new_edges[1:])).T
    new_exposure = np.array([np.sum(exposure[i[0]:i[1]+1]) for i in bounds_idx])
    new_counts = np.array([np.sum(counts[i[0]:i[1]+1]) for i in bounds_idx])
    new_edges = old_edges[new_edges]
    return new_counts, new_exposure, new_edges

def combine_into_one(counts, exposure, old_edges):
    """Combines binned data into a single bin

    Parameters:
    -----------
    counts: np.array
        The counts in each bin
    exposure: np.array
        The exposure of each bin
    old_edges: np.array
        The time edges of each bin
           
    Returns:
    -----------
    np.array
        The counts in each new bin
    np.array
        The exposure of each new bin
    np.array
        The new time edges of the rebinned data
    """
    new_counts = np.array([np.sum(counts)])
    new_exposure = np.array([np.sum(exposure)])
    new_edges = old_edges[[0,-1]]
    return new_counts, new_exposure, new_edges

def rebin_by_snr(counts, exposure, old_edges, background_counts, snr):
    """Rebins binned data such that each bin is above a minimum signal-to-noise ratio
    
    
    Parameters:
    -----------
    counts: np.array
        The counts in each bin
    exposure: np.array
        The exposure of each bin
    old_edges: np.array
        The time edges of each bin
    background_counts: np.array
        The background counts in each bin
    snr: float
        The minimum signal-to-ratio threshold
       
    Returns:
    -----------
    np.array
        The counts in each new bin
    np.array
        The exposure of each new bin
    np.array
        The new time edges of the rebinned data
    """
    num_old_bins = counts.size
    # make sure there is a non-zero background
    mask = (background_counts > 0.0)
    if np.sum(mask) == 0:
        raise ValueError('Background counts are all non-positive.  Cannot bin by SNR.')
    
    # cycle through current bins and combine bins until we have exceeded the requested snr
    new_edges = [0]
    istart = 0
    while True:
        countscum = np.cumsum(counts[istart:])
        backgroundscum = np.cumsum(background_counts[istart:])
        snrcum = np.cumsum((countscum-backgroundscum)/np.sqrt(backgroundscum))
        iend = istart + np.sum(snrcum <= snr)-1
        if iend < istart:
            iend = istart     
        new_edges.append(iend)
        if (iend >= num_old_bins-1):
            break
        istart = iend+1
    
    # create the new rebinned arrays
    new_edges = np.array(new_edges)
    bounds_idx = np.array((new_edges[0:-1], new_edges[1:])).T
    new_exposure = np.array([np.sum(exposure[i[0]:i[1]+1]) for i in bounds_idx])
    new_counts = np.array([np.sum(counts[i[0]:i[1]+1]) for i in bounds_idx])
    new_edges = old_edges[new_edges]
    return new_counts, new_exposure, new_edges

def rebin_by_edge_index(counts, exposure, old_edges, new_edge_index):
    """Rebins binned data based on an array of bin edge indices
    
    Parameters:
    -----------
    counts: np.array
        The counts in each bin
    exposure: np.array
        The exposure of each bin
    old_edges: np.array
        The time edges of each bin
    new_edge_index: np.array
        The edge indices of the for the new binned data
       
    Returns:
    -----------
    np.array
        The counts in each new bin
    np.array
        The exposure of each new bin
    np.array
        The new time edges of the rebinned data
    """
    # new edges
    num_bins = new_edge_index.size-1
    nei = new_edge_index.astype(dtype=int)
    new_edges = old_edges[nei]
    
    # combine the counts and exposure
    new_exposure = np.zeros(num_bins, dtype=float)
    new_counts = np.zeros(num_bins, dtype=float)
    for i in range(num_bins):
        start_idx = new_edge_index[i]
        end_idx = new_edge_index[i+1]
        new_counts[i] = np.sum(counts[start_idx:end_idx])
        new_exposure[i] = np.sum(exposure[start_idx:end_idx])
    
    return new_counts, new_exposure, new_edges

def rebin_by_edge(counts, exposure, old_edges, new_edges):
    """Rebins binned data based on an array of bin edge indices
    
    Parameters:
    -----------
    counts: np.array
        The counts in each bin
    exposure: np.array
        The exposure of each bin
    old_edges: np.array
        The time edges of each bin
    new_edges: np.array
        The new edges of the binned data
       
    Returns:
    -----------
    np.array
        The counts in each new bin
    np.array
        The exposure of each new bin
    np.array
        The new time edges of the rebinned data
    """
    # new edges
    num_bins = new_edges.size-1

    # combine the counts and exposure
    old_edges_list = old_edges.tolist()
    new_exposure = np.zeros(num_bins, dtype=float)
    new_counts = np.zeros(num_bins, dtype=float)
    for i in range(num_bins):
        start_idx = old_edges_list.index(new_edges[i])
        end_idx = old_edges_list.index(new_edges[i+1])
        new_counts[i] = np.sum(counts[start_idx:end_idx])
        new_exposure[i] = np.sum(exposure[start_idx:end_idx])
    
    return new_counts, new_exposure, new_edges