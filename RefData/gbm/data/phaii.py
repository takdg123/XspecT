import astropy.io.fits as fits
import numpy as np
from warnings import warn
from collections import OrderedDict


from .data import Data
from .containers import RateHisto, EventList
import gbm.data.headers as hdr
from gbm.time import Met

class PHAII(object):
    """PHAII container class for time history and spectra

    Attributes:
    -----------
    energy_range: (float, float)
        full energy range of the PHAII data
    header: OrderedDict
        The headers of the PHAII data
    num_chans: int
        number of energy channels
    time_range: (float, float)
        full time range of the PHAII data in MET
    trigtime: float
        The trigger time
    type: The type of data (binned or unbinned)
        
    Public Methods:
    ---------------
    view:
        Generates a view of PHAII data by making an optional time and energy slice and
        then integrating over either the time or energy axis.
    write:
        Write the PHAII data to a GBM FITS file
    
    Class Methods:
    ---------------
    create:
        Create a full PHAII object when providing the appropriate record arrays 

    """
    def __init__(self):
        
        self.header = OrderedDict()
        self.ebounds = None
        self.spectrum = None
        self.gti = None

        self.num_chans = None
        self.energy_range = None
        self.time_range = None
        self.trigtime = None
        self.type = None
        
        if self.filename is not None:
            self._from_file()
            self._set_properties()
    
    def __str__(self):
        return '{0} data from {1}'.format(self.data_type, self.detector)

    def _from_file(self):
        """Construct PHAII object from a GBM PHAII FITS file
        """
        with fits.open(self.filename) as hdulist:
            for hdu in hdulist:
                self.header.update({hdu.name: hdu.header})
            
            self.ebounds = hdulist['EBOUNDS'].data
            self.spectrum = hdulist['SPECTRUM'].data
            self.gti = hdulist['GTI'].data
        
    def _set_properties(self):
        """Set the higher level PHAII properties
        """
        self.num_chans = self.ebounds.size
        self.energy_range = self._get_energy_range(self.ebounds)
        self.time_range = self._get_time_range(self.header)
        self._assert_exposure()
        if 'TRIGTIME' in self.header['PRIMARY'].keys():
            self._relative_to_trigger()
      
    def _relative_to_trigger(self):
        """If there is a trigger time, scale the times in the file relative
        to the trigger time
        """
        self.trigtime = self.header['PRIMARY']['TRIGTIME']
        self.spectrum['TIME'] = self.spectrum['TIME']-self.trigtime
        self.spectrum['ENDTIME'] = self.spectrum['ENDTIME']-self.trigtime
        self.time_range = (self.header['PRIMARY']['TSTART']-self.trigtime, 
                           self.header['PRIMARY']['TSTOP']-self.trigtime)
        self.gti['START'] = self.gti['START']-self.trigtime
        self.gti['STOP'] = self.gti['STOP']-self.trigtime

    def _assert_exposure(self):
        """Sometimes the data files contain a negative exposure. That is very 
        very naughty, and we don't like dealing with naughty things.
        """
        if self.spectrum is not None:
            mask = self.spectrum['EXPOSURE'] > 0.0
            self.spectrum = self.spectrum[mask]

    def _get_energy_range(self, ebounds):
        return ebounds['E_MIN'][0], ebounds['E_MAX'][-1]

    def _get_time_range(self, header):
        return header['PRIMARY']['TSTART'], header['PRIMARY']['TSTOP']

    def _slice_time(self, counts_data, tstart, tstop):
        """Slices the PHAII data by time

        Parameters:
        -----------
        counts_data: np.recarray
            The time history data
        tstart: float
            The start of the time slice, in MET
        tstop
            The end of the time slice, in MET

        Returns:
        -----------
        sliced_counts_data: np.recarray
            The sliced time history data
        """
        mask = (counts_data['ENDTIME'] > tstart) & (counts_data['TIME'] < tstop)
        #sliced_counts_data = counts_data[mask]
        # there is a bug in some numpy versions that causes the above line to
        # do the wrong thing and/or throw a warning/error
        # So now we have to do this the hard way
        numrows = np.sum(mask)
        numchans = counts_data['COUNTS'].shape[1]
        dtypes = [('COUNTS', '>f4', (numchans,)), ('EXPOSURE', '>f4'), 
                  ('QUALITY', '>i2'), ('TIME', '>f8'), ('ENDTIME', '>f8')]
        sliced_counts_data = np.recarray((numrows,), dtype=dtypes)
        sliced_counts_data['COUNTS'] = counts_data['COUNTS'][mask,:]
        sliced_counts_data['EXPOSURE'] = counts_data['EXPOSURE'][mask]
        sliced_counts_data['QUALITY'] = counts_data['QUALITY'][mask]
        sliced_counts_data['TIME'] = counts_data['TIME'][mask]
        sliced_counts_data['ENDTIME'] = counts_data['ENDTIME'][mask]
        return sliced_counts_data

    def _slice_energy(self, counts_data, min_channel, max_channel):
        """Slices the PHAII data by energy
        
        Parameters:
        -----------
        counts_data: np.recarray
            The time history data
        min_channel: float
            The minimum channel energy
        max_channel: float
            The maximum channel energy

        Returns:
        -----------
        sliced_counts_data: np.recarray
            The sliced time history data
        """
        # Note: originally doing a resize on the counts_data['COUNTS']
        # so that the channel-sliced counts array can be placed in it.
        # However, it fails sometimes with a 'ValueError: resize only works on 
        # single-segment arrays' and it's not clear why.  So we have to do this
        # the hard way to make sure that it always works.
        mask = (self.ebounds['CHANNEL'] >= min_channel) & \
               (self.ebounds['CHANNEL'] <= max_channel)
        numrows = counts_data['TIME'].shape[0]
        numchans = np.sum(mask)
        dtypes = [('COUNTS', '>f4', (numchans,)), ('EXPOSURE', '>f4'), 
                  ('QUALITY', '>i2'), ('TIME', '>f8'), ('ENDTIME', '>f8')]
        sliced_counts_data = np.recarray((numrows,), dtype=dtypes)
        sliced_counts_data['COUNTS'] = counts_data['COUNTS'][:,mask]
        sliced_counts_data['EXPOSURE'] = counts_data['EXPOSURE']
        sliced_counts_data['QUALITY'] = counts_data['QUALITY']
        sliced_counts_data['TIME'] = counts_data['TIME']
        sliced_counts_data['ENDTIME'] = counts_data['ENDTIME']
        return sliced_counts_data

    def _integrate_time(self, counts_data):
        """Integrate the PHAII data over time

        Parameters:
        -----------
        counts_data: np.recarray
            The time history data

        Returns:
        -----------
        integrated_data: np.recarray
            The time-integrated data
        """
        num_chans = counts_data['COUNTS'].shape[1]
        dtypes = [('COUNTS', '>f4', (num_chans,)), ('EXPOSURE', '>f4'), 
                  ('QUALITY', '>i2'), ('TIME', '>f8'), ('ENDTIME', '>f8')]
        integrated_data = np.recarray((1,), dtype=dtypes)
        integrated_data['COUNTS'] = np.sum(counts_data['COUNTS'], axis=0)
        integrated_data['EXPOSURE'] = np.sum(counts_data['EXPOSURE'])
        integrated_data['QUALITY'] = 1
        integrated_data['TIME'] = counts_data['TIME'][0]
        integrated_data['ENDTIME'] = counts_data['ENDTIME'][-1]
        return integrated_data

    def _integrate_energy(self, counts_data):
        """Integrate the PHAII data over energy

        Parameters:
        -----------
        counts_data: np.recarray
            The time history data

        Returns:
        -----------
        integrated_data: np.recarray
            The energy-integrated data
        """
        num_bins = counts_data.size
        dtypes = [('COUNTS', '>f4'), ('EXPOSURE', '>f4'), ('QUALITY', '>i2'), 
                  ('TIME', '>f8'), ('ENDTIME', '>f8')]
        integrated_data = np.recarray((num_bins,), dtype=dtypes)
        integrated_data['COUNTS'] = np.sum(counts_data['COUNTS'], axis=1)
        integrated_data['EXPOSURE'] = counts_data['EXPOSURE']
        integrated_data['QUALITY'] = counts_data['QUALITY']
        integrated_data['TIME'] = counts_data['TIME']
        integrated_data['ENDTIME'] = counts_data['ENDTIME']
        return integrated_data

    def _assert_time_range(self, time_range):
        """Ensure the time range is of correct size and order and sanitized

        Parameters:
        -----------
        time_range: tuple
            The input time range

        Returns:
        -----------
        tstart: float
            The start time
        tstop: float
            The stop time
        """
        if time_range is None:
            time_range = self.time_range
        assert time_range[0] <= time_range[1]
        tstart = time_range[0]
        tstop = time_range[1]
        return (tstart, tstop)

    def _assert_energy_range(self, energy_range):
        """Ensure the energy range is of correct size and order and sanitized 
        and converts to the channel number

        Parameters:
        -----------
        energy_range: tuple
            The input energy range

        Returns:
        -----------
        min_chan: int
            The minimum channel number in the energy range
        max_chan: int
            The maximum channel number in the energy range
        """
        assert energy_range[0] <= energy_range[1]
        min_energy = energy_range[0]
        max_energy = energy_range[1]
        assert (min_energy > 0) & (max_energy > 0), \
               'Energies must be positive'
        if min_energy < self.energy_range[0]:
            warn("Requested min_energy is below the energy range of the data. "
                 "Will use {0} as min_energy".format(self.energy_range[0]))
            min_energy = self.energy_range[0]
        if max_energy > self.energy_range[1]:
            warn("Requested max_energy is above the energy range of the data. "
                "Will use {0} as max_energy".format(self.energy_range[0]))
            max_energy = self.energy_range[1]
        mask = (self.ebounds['E_MIN'] <= min_energy) & (self.ebounds['E_MAX'] > min_energy)
        min_chan = self.ebounds['CHANNEL'][mask][0]
        mask = (self.ebounds['E_MIN'] < max_energy) & (self.ebounds['E_MAX'] >= max_energy)
        max_chan = self.ebounds['CHANNEL'][mask][0]
        return (min_chan, max_chan)

    def _assert_channel_range(self, channel_range):
        """Ensure the channel range is of correct size and order and sanitized

        Parameters:
        -----------
        channel_range: tuple
            The input channel range

        Returns:
        -----------
        min_chan:
            The minimum energy channel
        max_chan:
            The maximum energy channel
        """
        assert channel_range[0] <= channel_range[1]
        min_chan = int(channel_range[0])
        max_chan = int(channel_range[1])
        assert (min_chan >= 0) & (max_chan >= 0), \
               'Energy channels numbers must be non-negative'
        if max_chan > self.num_chans-1:
            warn("Requested max_chan is too large. Will use {0} as max_chan".format(self.num_chans-1))
            max_chan = self.num_chans-1
        return (min_chan, max_chan)

    def _channel_range_to_energy(self, channel_range):
        """Convert channel range to energy range

        Parameters:
        -----------
        channel_range: tuple
            The input channel range

        Returns:
        -----------
        min_energy: float
            The minimum energy bound
        max_energy: float
            The maximum energy bound
        """
        min_channel = int(channel_range[0])
        max_channel = int(channel_range[1])
        min_energy = self.ebounds['E_MIN'][min_channel]
        max_energy = self.ebounds['E_MAX'][max_channel]
        return (min_energy, max_energy)

    def _split_on_saa(self, counts, bounds, exposure):
        """Splits the data over SAA boundaries

        Parameters:
        -----------
        counts: np.array
            The counts in each bin
        bounds: tuple
            The bin bounds
        exposure: np.array
            The bin exposures

        Returns:
        -----------
        split_counts: list
            The counts array, split over SAA
        split_bounds: list
            The bounds tuple, split over SAA
        split_exposure: list
            The exposute array, split over SAA
        """
        split_counts = []
        split_bounds = []
        split_exposure = []
        for gti in self.gti:
            mask = (bounds[0] < gti[1]) & (bounds[1] > gti[0])
            if np.sum(mask) > 0:
                split_counts.append(counts[mask])
                split_bounds.append((bounds[0][mask], bounds[1][mask]))
                split_exposure.append(exposure[mask])

        return (split_counts, split_bounds, split_exposure)

    def _split_bin_index_on_saa(self, bounds, index_list):
        """Splits the bin index list over SAA boundaries

        Parameters:
        -----------
        bounds: tuple
            The bin bounds
        index_list: np.array
            The bin index list

        Returns:
        -----------
        split_indices: list
            The index list, split over SAA
        """
        split_indices = []
        startidx = 0
        numgti = len(self.gti)
        for i in range(numgti):
            numbins = np.sum((bounds[i][0] >= self.gti[i][0]) & \
                             (bounds[i][1] <= self.gti[i][1]))
            mask = (index_list >= startidx) * (index_list < startidx+numbins)
            split_indices.append(index_list[mask]-startidx)
            startidx = startidx+numbins
        return split_indices
    
    def _split_bin_index_energy_bounds(self, bounds, index_list):
        """Splits the bin index list over energy selection bounds

        Parameters:
        -----------
        bounds: list
            The energy selection bounds
        index_list: np.array
            The bin index list

        Returns:
        -----------
        split_indices: list
            The index list, split over the energy selection bounds
        """
        index_list = np.asarray(index_list)
        split_indices = []
        startidx = 0
        num_bounds = len(bounds)
        for i in range(num_bounds):
            numbins = np.sum((self.ebounds['E_MIN'] >= bounds[i][0]) & \
                             (self.ebounds['E_MAX'] <= bounds[i][1]))
            mask = (index_list >= startidx) * (index_list < startidx+numbins)
            split_indices.append(index_list[mask]-startidx)
            startidx += numbins
        return split_indices
    
    def view(self, time_range=None, energy_range=None, channel_range=None, 
             integration_axis='energy'):
        """Generate a view of the PHAII data by making an optional time and energy slice 
        and then integrating over either the time or energy axis.

        Parameters:
        -----------
        time_range: (float, float), optional
            The MET time range of the view: (tstart, tstop).  If not set, then uses the 
            full time range of the data.
        energy_range: (float, float), optional
            The energy range of the view: (energy_lo, energy_hi).  If not set, then uses
            the full energy range of the data.  This is ignored if channel_range is set.
        channel_range: (int, int), optional
            The channel range of the view: (channel_lo, channel_hi).  If not set, then uses
            the full channel range of the data.
        integration_axis: str, optional
            The axis to integrate over.  Options are 'time' or 'energy'.  Default is
            'energy'.
        
        Returns:
        -----------
        rates: list
            A list of RateHisto objects containing the PHAII view.
        """
        # make sure integration axis is valid
        if (integration_axis != 'time') & (integration_axis != 'energy'):
            raise ValueError("integration_axis keyword must either be 'time' or 'energy'")

        # make sure time range is valid
        tstart, tstop = self._assert_time_range(time_range)

        # make sure channel range is valid...
        if (energy_range is None) & (channel_range is None):
            channel_range = [0, self.num_chans-1]
        if channel_range is not None:
            min_chan, max_chan = self._assert_channel_range(channel_range)
            min_energy, max_energy = self._channel_range_to_energy(channel_range)

        # ...or make sure energy range is valid and convert to channel range
        if energy_range is not None:
            min_chan, max_chan = self._assert_energy_range(energy_range)
            min_energy, max_energy = energy_range

        # Here is the following process:
        # 1) Integrate over an axis and create a RateHisto from the full data
        # 2) Make the view slice from the RateHisto

        counts_data = self.spectrum.copy()
        # integrate over time?
        if integration_axis == 'time':
            sliced_data = self._slice_time(counts_data, tstart, tstop)
            sliced_data = self._integrate_time(sliced_data)
            bounds = (self.ebounds['E_MIN'], self.ebounds['E_MAX'])
            binwidths = bounds[1]-bounds[0]
            exposure = sliced_data['EXPOSURE'][0]*binwidths
            rates = [RateHisto(sliced_data['COUNTS'].squeeze(), bounds, 
                               exposure=exposure)]
        # integrate over energy?
        elif integration_axis == 'energy':
            counts_data = self._slice_energy(counts_data, min_chan, max_chan)
            counts_data = self._integrate_energy(counts_data)
            bounds = (counts_data['TIME'], counts_data['ENDTIME'])
            # split the data over SAA traversal - no data
            counts, bounds, exposure = self._split_on_saa(counts_data['COUNTS'],
                                                          bounds, 
                                                          exposure=counts_data['EXPOSURE'])
            num_intervals = len(counts)
            rates = []
            for i in range(num_intervals):
                rates.append(RateHisto(counts[i], bounds[i], exposure=exposure[i]))

        # finally do the view slices                     
        if integration_axis == 'energy':
            rates = [rate.slice((tstart, tstop)) for rate in rates]
        elif integration_axis == 'time':
            rates = [rate.slice((min_energy, max_energy)) for rate in rates]
        # filter out null rates (things that were outside the slice)
        rates = [rate for rate in rates if rate is not None]

        return rates
    
    @classmethod
    def create(cls, events=None, ebounds=None, spectrum=None, gti=None, 
               trigtime=None, tstart=None, tstop=None, detector=None, object=None, 
               ra_obj=None, dec_obj=None, err_rad=None, detchans=None):
        """Create a PHAII object from the given record arrays

        Parameters:
        -----------
        events: np.recarray, optional
            A TTE events record array
        ebounds: np.recarray, optional
            An Ebounds record array
        spectrum: np.recarray, optional
            A counts/channel record array
        gti: np.recarray, optional
            A GTI record array
        trigtime: float, optional
            Trigger time
        tstart: float, optional
            The starting time of the file
        tstop: float, optional
            The stopping time of the file
        detector: str, optional
            The detector to which the data belongs
        object: str, optional
            The object class for which the data describes
        ra_obj: float, optional
            The RA of the object
        dec_obj: float, optional
            The Dec of the object
        err_rad: float, optional
            The circular equivalent uncertainty in the location of the object
        detchans: int, optional
            The number of detector channels
        
        Returns:
        -----------
        obj: PHAII
            The created PHAII object
        """
        obj = cls(None)
        
        # create the primary extension
        filetype='PHAII'
        primary_header = hdr.primary(detnam=detector, filetype=filetype, tstart=tstart,
                                     tstop=tstop, trigtime=trigtime, object=object, 
                                     ra_obj=ra_obj, dec_obj=dec_obj, err_rad=err_rad)
        headers = [primary_header]
        header_names = ['PRIMARY']
        
        # ebounds extension
        if ebounds is not None:
            ebounds_header = hdr.ebounds(detnam=detector, tstart=tstart, tstop=tstop, 
                                         trigtime=trigtime, object=object, 
                                         ra_obj=ra_obj, dec_obj=dec_obj, 
                                         err_rad=err_rad, detchans=detchans)
            headers.append(ebounds_header)
            header_names.append('EBOUNDS')
            obj.ebounds = ebounds
        
        # spectrum extension
        if spectrum is not None:
            spectrum_header = hdr.spectrum(detnam=detector, tstart=tstart, tstop=tstop, 
                                           trigtime=trigtime, object=object, 
                                           ra_obj=ra_obj, dec_obj=dec_obj, 
                                           err_rad=err_rad, detchans=detchans)
            headers.append(spectrum_header)
            header_names.append('SPECTRUM')
            obj.spectrum = spectrum
        
        # events extension
        if events is not None:
            events_header = hdr.events(detnam=detector, tstart=tstart, tstop=tstop, 
                                           trigtime=trigtime, object=object, 
                                           ra_obj=ra_obj, dec_obj=dec_obj, 
                                           err_rad=err_rad, detchans=detchans)
            headers.append(events_header)
            header_names.append('EVENTS')
            obj.events = events

        # gti extension
        if gti is not None:
            gti_header = hdr.gti(detnam=detector, tstart=tstart, tstop=tstop, 
                                 trigtime=trigtime, object=object, ra_obj=ra_obj, 
                                 dec_obj=dec_obj, err_rad=err_rad)
            headers.append(gti_header)
            header_names.append('GTI')
            obj.gti = gti
        
        # store headers and set data properties 
        obj.header = {name:header for name, header in zip(header_names, headers)}           
        obj._set_properties()
        
        # set file properties
        obj.is_gbm_file = True
        obj.detector = detector
        if trigtime is not None:
            obj.is_trigger = True
            met = Met(trigtime)
            obj.id = met.bn
        else:
            obj.is_trigger = False
            met = Met(tstart)
            obj.id = met.ymd_h
        return obj
    
    def write(self, directory, filename=None):
        """Write the PHAII object to a GBM FITS file

        Parameters:
        -----------
        directory: str
            The directory to which the file is written
        filename: str, optional
            The filename.  If omitted, will construct a standardized filename
        """
        if (self.filename is None) and (filename is None):
            raise NameError('Filename not set')
        if filename is None:
            filename = self.filename_obj.basename()
        self.set_filename(filename, directory=directory)
        self.header['PRIMARY']['FILENAME'] = filename
        
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=self.header['PRIMARY'])
        hdulist.append(primary_hdu)
        
        try:
            ebounds_hdu = fits.BinTableHDU(data=self.ebounds, name='EBOUNDS',
                                           header=self.header['EBOUNDS'])
            hdulist.append(ebounds_hdu)
        except:
            pass
        
        try:
            spectrum_hdu = fits.BinTableHDU(data=self.spectrum, name='SPECTRUM', 
                                            header=self.header['SPECTRUM'])
            hdulist.append(spectrum_hdu)
        except:
            pass

        try:
            events_hdu = fits.BinTableHDU(data=self.events, name='EVENTS',
                                          header=self.header['EVENTS'])
            hdulist.append(spectrum_hdu)
        except:
            pass
        
        try:
            gti_hdu = fits.BinTableHDU(data=self.gti, header=self.header['GTI'], 
                                       name='GTI')
            hdulist.append(gti_hdu)
        except:
            pass
        
        hdulist.writeto(self.filename, checksum=True)

class Ctime(Data, PHAII):
    """Container class for a CTIME PHAII file
    """
    def __init__(self, filename):
        Data.__init__(self, filename)
        PHAII.__init__(self)
        self.type = 'binned'

    @classmethod
    def create(cls, ebounds=None, spectrum=None, gti=None, trigtime=None, 
               tstart=None, tstop=None, detector=None, object=None, 
               ra_obj=None, dec_obj=None, err_rad=None):
        
        if tstart is None:
            if spectrum is not None:
                tstart = spectrum['TIME'][0]
        if tstop is None:
            if spectrum is not None:
                tstop = spectrum['ENDTIME'][-1]
        obj = super(Ctime, cls).create(ebounds=ebounds, spectrum=spectrum, gti=gti, 
                                 trigtime=trigtime, tstart=tstart, tstop=tstop, 
                                 detector=detector, object=object, ra_obj=ra_obj, 
                                 dec_obj=dec_obj, err_rad=err_rad, detchans=8)
        obj.data_type='ctime'
        obj.create_filename(extension='pha')
        return obj

class Cspec(Data, PHAII):
    """Container class for a CSPEC PHAII file
    """
    def __init__(self, filename):
        Data.__init__(self, filename)
        PHAII.__init__(self)
        self.type = 'binned'

    @classmethod
    def create(cls, ebounds=None, spectrum=None, gti=None, trigtime=None, 
               tstart=None, tstop=None, detector=None, object=None, 
               ra_obj=None, dec_obj=None, err_rad=None):
        
        if tstart is None:
            if spectrum is not None:
                tstart = spectrum['TIME'][0]
        if tstop is None:
            if spectrum is not None:
                tstop = spectrum['ENDTIME'][-1]
        obj = super(Cspec, cls).create(ebounds=ebounds, spectrum=spectrum, gti=gti, 
                                 trigtime=trigtime, tstart=tstart, tstop=tstop, 
                                 detector=detector, object=object, ra_obj=ra_obj, 
                                 dec_obj=dec_obj, err_rad=err_rad, detchans=128)
        obj.data_type='cspec'
        obj.create_filename(extension='pha')
        return obj

class TTE(Data, PHAII):
    """Container class for a TTE FITS file
    """
    def __init__(self, filename):
        Data.__init__(self, filename)
        self.events = None
        PHAII.__init__(self)
        #self.spectrum = None
        self.type = 'unbinned'

    def _from_file(self):
        with fits.open(self.filename) as hdulist:
            for hdu in hdulist:
                self.header.update({hdu.name: hdu.header})
            
            self.events = hdulist['EVENTS'].data
            self.ebounds = hdulist['EBOUNDS'].data
            self.gti = hdulist['GTI'].data

    def _relative_to_trigger(self):
        """If there is a trigger time, scale the times in the file relative
        to the trigger time
        """
        self.trigtime = self.header['PRIMARY']['TRIGTIME']
        self.events['TIME'] = self.events['TIME']-self.trigtime
        self.time_range = (self.header['PRIMARY']['TSTART']-self.trigtime, 
                           self.header['PRIMARY']['TSTOP']-self.trigtime)
        self.gti['START'] = self.gti['START']-self.trigtime
        self.gti['STOP'] = self.gti['STOP']-self.trigtime

    def _integrate_time(self, events_data):
        """Integrate the TTE data over time

        Parameters:
        -----------
        events_data: np.recarray
            The TTE event data

        Returns:
        -----------
        integrated_data: np.recarray
            The time-integrated data
        """
        num_chans = self.num_chans
        dtypes = [('COUNTS', '>f4', (num_chans,)), ('EXPOSURE', '>f4'), 
                  ('QUALITY', '>i2'), ('TIME', '>f8'), ('ENDTIME', '>f8')]
        integrated_data = np.recarray((1,), dtype=dtypes)

        dt = np.max(events_data['TIME'])-np.min(events_data['TIME'])
        channel_edges = np.arange(num_chans+1)
        counts, _ = np.histogram(events_data['PHA'], channel_edges)

        channels = channel_edges[:-1]
        #channels, counts = np.unique(events_data['PHA'], return_counts=True)
        # 2.6 microseconds deadtime for each count...
        deadtime = 0.0
        mask = (channels != 127)
        deadtime += np.sum(counts[mask])*2.6e-6
        # ...except in the last channel which is 10 microseconds for each count
        deadtime += np.sum(counts[~mask])*1e-5

        # fill up counts data recarray
        integrated_data['COUNTS'] = counts
        integrated_data['EXPOSURE'] = dt-deadtime
        integrated_data['QUALITY'] = 1
        integrated_data['TIME'] = np.min(events_data['TIME'])
        integrated_data['ENDTIME'] = np.max(events_data['TIME'])
        return integrated_data

    def calculate_exposure(self, bin_edges):
        """Calculate the exposure given bin edges

        Parameters:
        -----------
        bin_edges: np.array
            The time bin edges

        Returns:
        -----------
        exposure: np.array
            The exposure corresponding to the bin edges
        """
        
        events = EventList(self.events)
        events = events.slice('TIME', (bin_edges[0], bin_edges[-1]))
        
        # all counts
        counts = np.histogram(events._events['TIME'], bin_edges)[0]
        
        # all overflow counts
        mask = (events._events['PHA'] == 127)
        overflow_counts = np.histogram(events._events['TIME'][mask], bin_edges)[0]
        
        # deadtime and livetime
        deadtime = counts*2.6e-6 + overflow_counts*7.4e-6
        exposure = (bin_edges[1:]-bin_edges[:-1]) - deadtime
        
        return exposure
    
    def view(self, time_range=None, energy_range=None, channel_range=None, 
             integration_axis='energy'):
        """Generate a view of the TTE data by making an optional time and energy slice 
        and then integrating over either the time or energy axis.
        
        This method overrides the PHAII view method

        Parameters:
        -----------
        time_range: (float, float), optional
            The MET time range of the view: (tstart, tstop).  If not set, then uses the 
            full time range of the data.
        energy_range: (float, float), optional
            The energy range of the view: (energy_lo, energy_hi).  If not set, then uses
            the full energy range of the data.  This is ignored if channel_range is set.
        channel_range: (int, int), optional
            The channel range of the view: (channel_lo, channel_hi).  If not set, then uses
            the full channel range of the data.
        integration_axis: str, optional
            The axis to integrate over.  Options are 'time' or 'energy'.  Default is
            'energy'.
        
        Returns:
        -----------
        rates: list
            A list of RateHisto objects containing the PHAII view.  
            Alternatively, if a time binning method is not set, will return a 
            list of EventLists
        """
        # make sure integration axis is valid
        if (integration_axis != 'time') & (integration_axis != 'energy'):
            raise ValueError("integration_axis keyword must either be 'time' or 'energy'")

        # make sure time range is valid
        tstart, tstop = self._assert_time_range(time_range)

        # make sure channel range is valid...
        if (energy_range is None) & (channel_range is None):
            channel_range = [0, self.num_chans-1]
        if channel_range is not None:
            min_chan, max_chan = self._assert_channel_range(channel_range)
            min_energy, max_energy = self._channel_range_to_energy(channel_range)

        # ...or make sure energy range is valid and convert to channel range
        if energy_range is not None:
            min_chan, max_chan = self._assert_energy_range(energy_range)
            min_energy, max_energy = energy_range

        # Here is the following process:
        # 1) Integrate over an axis and create a RateHisto from the full data
        # 2) Make the view slice from the RateHisto

        eventlist = EventList(self.events)
        # integrate over time?
        if integration_axis == 'time':
            eventlist = eventlist.slice('TIME', (tstart, tstop))
            counts_data = self._integrate_time(eventlist._events)
            bounds = (self.ebounds['E_MIN'], self.ebounds['E_MAX'])
            binwidths = bounds[1]-bounds[0]
            exposure = counts_data['EXPOSURE'][0]*binwidths
            rates = [RateHisto(counts_data['COUNTS'].squeeze(), bounds, 
                               exposure=exposure)]
        # integrate over energy?
        elif integration_axis == 'energy':
            eventlist = eventlist.slice('PHA', (min_chan, max_chan))

        # finally do the view slices                     
        if integration_axis == 'energy':
            rates = [eventlist.slice('TIME', (tstart, tstop))]
        elif integration_axis == 'time':
            rates = [rate.slice((min_energy, max_energy)) for rate in rates]
        # filter out null rates (things that were outside the slice)
        rates = [rate for rate in rates if rate is not None]

        return rates

    @classmethod
    def create(cls, events=None, ebounds=None, gti=None, trigtime=None, 
               tstart=None, tstop=None, detector=None, object=None, 
               ra_obj=None, dec_obj=None, err_rad=None, detchans=None):
        
        if tstart is None:
            if events is not None:
                tstart = events['TIME'][0]
        if tstop is None:
            if events is not None:
                tstop = events['TIME'][-1]
        obj = super(TTE, cls).create(ebounds=ebounds, events=events, gti=gti, 
                                     trigtime=trigtime, tstart=tstart, tstop=tstop, 
                                     detector=detector, object=object, ra_obj=ra_obj, 
                                     dec_obj=dec_obj, err_rad=err_rad, detchans=128)
        obj.data_type='tte'
        obj.create_filename(extension='fit')
        return obj
