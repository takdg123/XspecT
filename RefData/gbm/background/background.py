import numpy as np
from . import binned
from . import unbinned
from gbm.data.containers import RateHisto, EventList

class Background(object):
    """Class containing a background model, either fitted to data or estimated
    from data.  This class is a wrapper around arbitrary background fitting
    and estimation algorithms, providing a standardized interface.

    Attributes:
    -----------
    rates: np.array
        The fitted background rates
    uncertainty: np.array
        The uncertainty in the fitted/estimated background
    livetime: float
        The background livetime
    chisq: np.array
        The chi-squared statistic for each energy channel.  Is null if no
        chisq variable exists in the fitting method
    dof: np.array
        The degrees-of-freedm of the fit for each energy channel.  Is null if 
        no dof variable exists in the fitting method
    
    Public Methods:
    fit:
        Fit the data with a background model
    interpolate:
        Interpolate the background at requested times
    integrate_energy:
        The background model integrated over energy
    source_spectrum_integrated:
        The background model integrated over time as a spectrum
    source_spectrum_for_each_bin:
        The background model spectrum for each time bin in the selection
    """     
    def __init__(self, selection):
        """Initializes the background model with the background selection 
        from which the background will be estimated

        Parameters:
        -----------
        selection: Selection
            The background selection object
        """
        self._datatype = type(selection._data_obj).__name__
        self._bkgd_selection = selection
        self._bkgds = []
        self.rates = []
        self.uncertainty = []
        self.livetime = selection.livetime
        self.chisq = None
        self.dof = None
        
    def __str__(self):
        return 'Background model for {0}'.format(self._bkgd_selection._data_obj.detector)
    
    def fit(self, method, *args, datatype='binned'):
        """Fit/estimate the background for each channel using a 
        specified algorithm

        Parameters:
        -----------
        method: class
            The class name of the background fitting algorithm to use.  The 
            background fitting algorithm should be written as a class with a
            public 'fit' method that takes the counts, exposure, and bin edges 
            and returns the background rates and fit uncertainty in that order.
            Additionally, the class can contain an optional public 'interpolate'
            method which takes in an array of times outputs the interpolated
            background and uncertainty.  If the class contains variables named 
            "chisq" and "dof", then the chisq and dof values will be compiled
            for each energy channel.
        *args: 
            Any extra arguments to be passed to the __init__ method of the 
            background fitting algorithm 
        datatype: str
            'binned' for binned data and 'unbinned' for unbinned data
        """
        for channel in self._bkgd_selection.single_channels:
            if datatype == 'binned':
                lc = RateHisto.merge(channel)
                bkgd = method(*args)
                vals, errs = bkgd.fit(lc.values, lc.exposure, lc.bounds)
            elif datatype == 'unbinned':
                lc = EventList.merge(channel)
                bkgd = method(*args)
                vals, errs = bkgd.fit(lc._events['TIME'])
            # temp fix for unbinned
            #if self._datatype == 'TTE':
            #    vals, errs = bkgd.fit(channel._events['TIME'])
            #else:
            self.rates.append(vals)
            self.uncertainty.append(errs)
            self._bkgds.append(bkgd)
        
        # check if method has chisq and dof variables.  if so then grab for each
        # energy channels
        if hasattr(self._bkgds[0], 'chisq'):
            self.chisq = np.array([bkgd.chisq for bkgd in self._bkgds])
        if hasattr(self._bkgds[0], 'dof'):
            self.dof = np.array([bkgd.dof for bkgd in self._bkgds])
    
    def integrate_energy(self, selection):
        """Interpolate the background model over the provided selection and 
        integrate over energy.
        
        Parameters:
        -----------
        selection: Selection
            The selection over which to interpolate (e.g. full time history
            selection)
        
        Returns:
        -----------
        rates: np.array
            The background rates
        uncert: np.array
            The background uncertainties
        """
        # temp fix for unbinned
        if type(selection.lightcurve[0]).__name__ == 'EventList':
            lc = EventList.merge(selection.lightcurve)
        else:    
            lc = RateHisto.merge(selection.lightcurve)
        # interpolate the lightcurve for each channel
        rates, uncerts = self._interpolate_each_channel(lc)

        # integrate over energy channels
        rate = np.array(rates[0])
        variance = np.array(uncerts[0]**2)
        for i in range(1, len(rates)):
            rate += np.array(rates[i])
            variance += np.array(uncerts[i])**2
        uncert = np.sqrt(variance)
        
        return rate, uncert
    
    def interpolate(self, bkgd, times, *args):
        """Interpolate the background at requested times
        
        Parameters:
        -----------
        bkgd: object
            The background model object to use
        times: np.array
            An array of times at which to interpolate the background
        *args: 
            Any extra arguments to be passed to the 'interpolate' method of the 
            background fitting algorithm 
        
        Returns:
        -----------
        rates: np.array
            The background rates
        errs: np.array
            The background uncertainties
        """
        rates, errs = bkgd.interpolate(times, *args)
        return rates, errs
    
    def source_spectrum_for_each_bin(self, selection):
        """Interpolate the background model over the provided selection and 
        return a background count rate spectrum for each time bin in the selection.
        Note: These spectra contain zeros for channels outside the source selection
        
        Parameters:
        -----------
        selection: Selection
            The selection over which to interpolate (e.g. full time history
            selection)
        
        Returns:
        -----------
        rates: list of RateHisto
            A list of spectra, representing the count rate spectrum in each time
            bin
        """        
        lc = RateHisto.merge(selection.lightcurve)
        
        # interpolate the lightcurve for each channel
        rates, uncerts = self._interpolate_each_channel(lc)
        rates = np.array(rates)
        uncerts = np.array(uncerts)
        widths = lc.bin_widths()
        counts = rates * widths[np.newaxis,:]
        
        # spectral bin exposures (full energy range)
        bounds = (selection._data_obj.ebounds['E_MIN'], \
                  selection._data_obj.ebounds['E_MAX'])
        binwidths = bounds[1]-bounds[0]
        #exposure = binwidths[:,np.newaxis] * lc.exposure
        exposure = np.ones((binwidths.size, lc.num_bins)) * \
                   lc.exposure[np.newaxis,:] # specially made for writeXSPEC
        
        # pad the unused channels with zeros (needed to create the xspec bak files
        all_counts = np.zeros((binwidths.size, lc.num_bins))
        all_uncerts = np.zeros((binwidths.size, lc.num_bins))
        idx = selection._get_energy_channel_indices()
        all_counts[idx,:] = counts
        all_uncerts[idx,:] = uncerts
        
        # create a spectrum RateHisto for each time bin
        spectra = []
        for i in range(lc.num_bins):
            r = RateHisto(all_counts[:,i], bounds, exposure=exposure[:,i])
            r.uncertainty = all_uncerts[:,i]
            spectra.append(r)
        
        return spectra
    
    def source_spectrum_integrated(self, selection):
        """Interpolate the background model over the provided selection and 
        integrate over time, creating a background count rate spectrum.
        Note: This spectrum contains zeros for channels outside the source selection
        
        Parameters:
        -----------
        selection: Selection
            The selection over which to interpolate (e.g. full time history
            selection)
        
        Returns:
        -----------
        integrated_spectrum: RateHisto
            The time-integrated count rate spectrum
        """        
        # get the spectrum for each time bin as well as the total exposure
        spectra = self.source_spectrum_for_each_bin(selection)
        binwidths = spectra[0].bin_widths()
        #exposure = selection.livetime*binwidths
        exposure = float(selection.livetime) # specially made for writeXSPEC
        
        # calculate the background uncertainty
        uncerts = np.zeros_like(spectra[0].uncertainty)
        for spectrum in spectra:
            livetime = (spectrum.exposure/binwidths)[0]
            uncerts += (spectrum.uncertainty * livetime)**2
        uncerts = np.sqrt(uncerts/selection.livetime)
        
        # construct the time-integrated spectrum
        temp_spectrum = RateHisto.sum_energy(spectra)
        integrated_spectrum = RateHisto(temp_spectrum.values, 
                                             temp_spectrum.bounds,
                                             exposure=exposure)
        integrated_spectrum.uncertainty = uncerts
        
        return integrated_spectrum
        
    def _interpolate_each_channel(self, lc):
        """Interpolate the background model for the lightcurve over each channel
        
        Parameters:
        -----------
        lc: RateHisto
            The lightcurve over which to interpolate
        
        Returns:
        -----------
        rates: np.array
            The background rates
        uncert: np.array
            The background uncertainties
        """
        rates = []
        uncerts = []
        for bkgd in self._bkgds:
            # temp fix for unbinned
            if type(lc).__name__ == 'EventList':
                rate, uncert = self.interpolate(bkgd, lc._events['TIME'][:-1], 
                                                lc._events['TIME'][1:])
            else:
                rate, uncert = self.interpolate(bkgd, lc.bounds[0], lc.bounds[1])
            rates.append(rate)
            uncerts.append(uncert)
        return rates, uncerts
    
