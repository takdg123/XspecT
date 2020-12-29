import astropy.io.fits as fits
import numpy as np
from collections import OrderedDict

from gbm.detectors import Detector
from .data import Data
from .phaii import PHAII
from .poshist import PosHist

classifications = {0: 'ERROR', 1: 'UNRELOC', 2: 'LOCLPAR', 3: 'BELOWHZ', 
                   4: 'GRB', 5: 'SGR', 6: 'TRANSNT', 7: 'DISTPAR', 
                   8: 'SFL', 9: 'CYGX1', 10: 'SGR1806', 11: 'GROJ422',
                   19: 'TGF', 20: 'UNCERT', 21: 'GALBIN'}
# TODO: Need to expose the other extensions of the trigdat file     

class Trigdat(Data, PHAII):
    """Trigdat container class for trigger data

    Attributes:
    -----------
    trigtime: float
        The trigger MET
    header: list
    	The header information
	poshist: PosHist
		The position history for the trigdat
	        
    Public Methods:
    ---------------
    to_phaii:
    	Return a detector as a PHAII object
    """    
    def __init__(self, filename):
        self.header = OrderedDict()
        self.detectors = [det.short_name for det in Detector]
        # trigdat has an assumed (hard-coded) ebounds
        self.ebounds = np.rec.array([(0, 3.4, 10.0), 
                                     (1, 10.0, 22.0), 
                                     (2, 22.0, 44.0),
                                     (3, 44.0, 95.0),
                                     (4, 95.0, 300.0),
                                     (5, 300.0, 500.0),
                                     (6, 500.0, 800.0),
                                     (7, 800.0, 2000.0)], 
                                     formats='>i2,>f4,>f4',
                                     names='CHANNEL,E_MIN,E_MAX')
        self.ebounds = fits.BinTableHDU(data=self.ebounds).data
        
        Data.__init__(self, filename)
        if self.filename is not None:
            self._from_file()
            
    def _from_file(self):
        """Strips data from the file, and stores the header.

        Extracts the evntrate extensions and store the contained info as attributes
        """
        with fits.open(self.filename) as hdulist:
            for hdu in hdulist:
                self.header.update({hdu.name: hdu.header})
        
            self.trigrate = hdulist['TRIGRATE'].data
            self.backrate = hdulist['BCKRATES'].data
            self.onboard_calc = hdulist['OB_CALC'].data
            self.maxrates = hdulist['MAXRATES'].data
            self.lightcurve = hdulist['EVNTRATE'].data
        
        self.lightcurve.sort(order='TIME')
        
        # If read naively, the rates will be incorrect - need to be reshaped first
        self.rates = self.lightcurve['RATE'].reshape(-1, 14, 8)
        
        # create the position history
        time_idx, _ = self._time_indices(64)
        self.poshist = PosHist.from_trigdat(self.lightcurve['TIME'][time_idx], 
                                            self.lightcurve['SCATTITD'][time_idx,:].T, 
                                            self.lightcurve['EIC'][time_idx,:])
        
        self.trigtime = self.header['PRIMARY']['TRIGTIME']
        self.poshist.time -= self.trigtime
        return
    
    def to_phaii(self, detector, timescale=1024):
        """Return a detector as a PHAII object
        
        Parameters:
        -----------
		detector: string
			string name of the GBM detector
		timescale: int, optional
			Include data down to this timescale in milliseconds.  
			Possible values are 1024 (default), 256, and 64.  The 8192 ms timescale
			is always included.
			
		Returns:
		-----------
		phaii: PHAII
			The trigdat data for the detector in a PHAII object
        """
        detector = detector.lower()
        if detector not in self.detectors:
            raise ValueError('Illegal detector name')
        
        if (timescale != 1024) and (timescale != 256) and (timescale != 64):
            raise ValueError('Illegal Trigdat resolution. Available resolutions: \
                              1024, 256, 64')        

        # grab the correct detector rates (stored as 14 dets x 8 channels)
        det_idx = self.detectors.index(detector)

        # return the requested timescales
        time_idx, dt = self._time_indices(timescale)
        
        # calculate counts and exposure
        counts = self.rates[time_idx, det_idx,:] * dt[:,np.newaxis]
        exposure = self._calc_exposure(counts, dt)

        # the quality flags
        quality = np.zeros_like(exposure)
        quality[exposure <= 0.01] = 1
        quality[exposure > 33.0] = 1
        
        # create the spectrum record array
        numtimes = dt.size
        record = [(counts[i,:], exposure[i], quality[i], 
                   self.lightcurve['TIME'][time_idx[i]], 
                   self.lightcurve['ENDTIME'][time_idx[i]]) for i in range(numtimes)]
        spectrum = np.rec.array(record, dtype=[('COUNTS','(8,)>f4'), ('EXPOSURE','>f4'),
                                               ('QUALITY','>i2'), ('TIME','>f8'),
                                               ('ENDTIME','>f8')])
        spectrum = fits.BinTableHDU(data=spectrum).data
        
        # create the gti record array
        # Note: We don't actually have any way of determining the GTI
        tstart = self.header['PRIMARY']['TSTART']
        tstop = self.header['PRIMARY']['TSTOP']
        gti = np.rec.array([(tstart, tstop)], dtype=[('START', '>f8'), ('STOP', '>f8')])
        gti = fits.BinTableHDU(data=gti).data
        
        # create the PHAII object
        object = self.header['PRIMARY']['OBJECT']
        ra_obj = self.header['PRIMARY']['RA_OBJ']
        dec_obj = self.header['PRIMARY']['DEC_OBJ']
        err_rad = self.header['PRIMARY']['ERR_RAD']
        phaii = super(Trigdat, self).create(ebounds=self.ebounds, spectrum=spectrum, gti=gti, 
                             trigtime=self.trigtime, tstart=tstart, tstop=tstop,
                             detector=detector, object=object, ra_obj=ra_obj,
                             dec_obj=dec_obj, err_rad=err_rad, detchans=8)
        phaii.data_type='trigdat'
        phaii.create_filename(extension='fit')
        return phaii
            
    def _calc_exposure(self, counts, dt):
        """Calculate the exposure
            
		Parameters:
		-----------
		counts: np.array
		    The observed counts in each bin
		dt: np.array
		    The time bin widths 
			
		Returns:
		-----------
		exposure: np.array
		    The exposure of each bin
        """
        deadtime = np.copy(counts)
        deadtime[:,:7] *= 2.6e-6
        deadtime[:,7] *= 1e-5
        total_deadtime = np.sum(deadtime, axis=1)
        exposure = (1.0-total_deadtime)*dt
        return exposure
       
    def _time_indices(self, time_res):
        """Indices into the Trigdat arrays corresponding to the desired time resolution(s) 
            
		Parameters:
		-----------
		time_res: int
			string name of the GBM detector
		timescale: int, optional
			Include data down to this timescale in milliseconds.  
			
		Returns:
		-----------
		np.array
			indices into Trigdat arrays
		np.array
			the bin widths in seconds
        """
        # bin widths
        dt = np.round((self.lightcurve['ENDTIME']-self.lightcurve['TIME'])*1000)
        # background bins
        back_idx = np.where(dt == 8192)[0]
        # 1 s scale bins - this is the minimum amount returned
        idx = np.where(dt == 1024)[0]
        cnt = len(idx)
        # reconcile 8 s and 1 s data
        idx = self._reconcile_timescales(back_idx, idx)

        # reconcile 8 s + 1 s and 256 ms data
        if time_res <= 256:
            tidx = np.where(dt == 256)[0]
            idx = self._reconcile_timescales(idx, tidx)
        
        # reconcile 8 s + 1 s + 256 ms and 64 ms data
        if time_res == 64:
            tidx = np.where(dt == 64)[0]
            idx = self._reconcile_timescales(idx, tidx)
        
        # return reconciled indices
        return idx, np.reshape(dt[idx]/1000.0, len(idx))
    
    def _reconcile_timescales(self, idx1, idx2):
        """Reconcile indices representing different timescales and glue them 
           together to form a complete (mostly) continuous set of indices
            
		Parameters:
		-----------
		idx1: np.array
			indices of the "bracketing" timescale
		idx2: np.array
			indices of the "inserted" timescale 
			
		Returns:
		-----------
		idx: np.array
			indices of idx2 spliced into idx1
		"""
        # bin edges for both selections
        start_times1 = self.lightcurve['TIME'][idx1]
        end_times1 = self.lightcurve['ENDTIME'][idx1]
        start_times2 = self.lightcurve['TIME'][idx2]
        end_times2 = self.lightcurve['ENDTIME'][idx2]
        
        # find where bracketing timescale ends and inserted timescale begins
        start_idx = (np.where(end_times1 >= start_times2[0]))[0][0]
        idx = np.concatenate((idx1[0:start_idx], idx2))
        
        # find wehere inserted timescale ends and bracketing timescale begins again
        end_idx = (np.where(start_times1 >= end_times2[-1]))[0][0]
        idx = np.concatenate((idx, idx1[end_idx:]))
        
        return idx
    