import numpy as np
import astropy.io.fits as fits
import gbm.data.headers as hdr
from gbm.time import Met
from .data import Data

class Parameter(object):
    """A fit parameter class
    
    Attributes:
    -----------
    name: str
        The name of the parameter
    value: float
        The parameter fit centroid
    uncertainty: tuple
        A one or two element tuple corresponding to a symmetric (Gaussian) errors
        or asymmetric errors, respectively.  In the case of asymmetric errors,
        The order of the tuple is assumed to be (upper, lower).
    valid_range: tuple
        A two-element tuple representing the valid range for the parameter.
        Default is (-np.inf, np.inf)
    
    Public Methods:
    ---------------
    one_sigma_range:
        Return the 1 sigma range of the parameter fit
    to_fits_value:
        Return as a tuple to be used for a FITS file
    valid_value:
        Check if the parameter value is within the allowed parameter range
    
    Class Methods:
    --------------
    create:
        Construct a Parameter object from keywords or dictionary
    """    
    def __init__(self):
        self.name = ''
        self._value = None
        self._uncertainty = (None,)
        self.valid_range = (-np.inf, np.inf)

    def __str__(self):
        value, uncertainty = self._str_format()
        
        if len(uncertainty) == 1:
            s = '+/- {0}'.format(uncertainty[0])
        else:
            if uncertainty[0] == uncertainty[1]:
                s = '+/- {0}'.format(uncertainty[0])
            else:
                s = '+{0}/-{1}'.format(uncertainty[0], uncertainty[1])
        return '{0}: {1} {2}'.format(self.name, value, s)
    
    def _str_format(self):
        if (self.value > 0.005) and (self.uncertainty[0] > 0.005):
            value = '{0:.2f}'.format(self.value)
            uncertainty = tuple(['{0:.2f}'.format(u) for u in self.uncertainty])
        else:
            value = '{0:.2e}'.format(self.value)
            val_coeff, val_exp = value.split('e')
            val_exp = int(val_exp)
            uncertainty = ['{0:.2e}'.format(u) for u in self.uncertainty]
            uncert_coeff = []
            uncert_exp = []
            for uncert in uncertainty:
                uncert_coeff.append(uncert.split('e')[0])
                uncert_exp.append(int(uncert.split('e')[1]))
        return (value, uncertainty)
      
    def valid_value(self):
        """Check if the parameter value is within the allowed parameter range
        """
        if (self.value >= self.valid_range[0]) and \
           (self.value <= self.valid_range[1]):
            return True
        else:
            return False

    def one_sigma_range(self):
        """Return the 1 sigma range of the parameter fit
        """
        if len(self._uncertainty) == 1:
            return (self._value-self._uncertainty[0], self._value+self._uncertainty[0])
        else:
            return (self._value-self._uncertainty[1], self._value+self._uncertainty[0])
    
    def to_fits_value(self, asymmetric=False):
        """Return as a tuple to be used for a FITS file
        
        Parameters:
        -----------
        asymmetric: bool
            If True, then even if the uncertainty is symmetric, it returns
            it in the form of an asymmetric uncertainty. Default is False
        
        Returns:
        -----------
        tuple
            2-value tuple (value, uncertainty) or 3-value tuple 
            (value, +uncertainty, -uncertainty)   
        """
        if asymmetric and len(self.uncertainty) == 1:
            return (self.value, self.uncertainty[0], self.uncertainty[0])
        
        if len(self.uncertainty) == 2:
            return (self.value, self.uncertainty[0], self.uncertainty[1])
        else:
            return (self.value, self.uncertainty[0])
    
    def to_table_value(self):
        s = self.__str__()
        s = s.split(':')[1].strip()
        return s
     
    @classmethod
    def create(cls, **kwargs):
        obj = cls()
        obj._init_by_dict(kwargs)
        return obj  

    @property
    def value(self):
        return self._value
    @value.setter
    def value(self, val):
        try:
            val = float(val)
        except:
            raise ValueError('Parameter.value must be a real number')
        self._value = val
    
    @property
    def uncertainty(self):
        return self._uncertainty
    @uncertainty.setter
    def uncertainty(self, value):
        try:
            len(value)
        except:
            value = (value,)
        if len(value) > 2:
            raise ValueError('Parameter.uncertainty must be no more than two elements')
        for v in value:
            try:
                float(v)
            except:
                raise ValueError('Parameter.uncertainty must be a real number')
        self._uncertainty = value

    def _init_by_dict(self, values):
        for key, val in values.items():
            try:
                p = getattr(self, key)
                if isinstance(p, property):
                    p.__set__(self, val)
                else:
                    self.__setattr__(key, val)
            except AttributeError:
                raise ValueError("{} is not a valid attribute".format(key))

class ModelFit(object):
    """A model fit class
    
    Attributes:
    -----------
    name: str
        The name of the model
    time_range: (float, float)
        The time range of the fit
    parameters: list
        A list of Parameters
    photon_flux: Parameter
    energy_flux: Parameter
    photon_fluence: Parameter
    energy_fluence: Parameter
    photon_fluence_50_300: Parameter
    energy_fluence_50_300: Parameter
    duration_fluence: Parameter
    stat_name: str
        The name of the fit statistic
    stat_value: float
        The fit statistic
    dof: int
        The fit degrees-of-freedom
    covariance: np.array
        A (nparam, nparam) covariance matrix
    
    Public Methods:
    ---------------
    parameter_list:
        List the parameters in the model
    to_fits_row:
        Return the contained data as a FITS table row
    
    Class Methods:
    --------------
    create:
        Construct a ModelFit object from keywords or dictionary
    from_xspec_data_dict:
        Construct a ModelFit object from the XSPEC data dictionary
    """    
    def __init__(self):
        self.name = ''
        self.flux_energy_range = (None, None)
        self.time_range = None
        self.parameters = []
        self.photon_flux = Parameter.create(name='photon flux')
        self.energy_flux = Parameter.create(name='energy flux')
        self.photon_fluence = Parameter.create(name='photon fluence')
        self.energy_fluence = Parameter.create(name='energy fluence')
        self.photon_flux_50_300 = Parameter()
        self.energy_fluence_50_300 = Parameter()
        self.duration_fluence = Parameter()
        self.stat_name = None
        self.stat_value = None
        self.dof = None
        self._covariance = None

    def __str__(self):
        param_str = '\n   '.join([str(param) for param in self.parameters])
        return '{0}\n   {1}'.format(self.name, param_str) 
        
    def parameter_list(self):
        return [param.name for param in self.parameters]
    
    def to_fits_row(self):
        """Return the contained data as a FITS table row
        
        Returns:
        -----------
        hdu: astropy.io.fits.BinTableHDU
            The FITS table
        """       
        numparams = len(self.parameters)
        cols = []
        cols.append(fits.Column(name='TIMEBIN', format='2D', array=[self.time_range]))
        i = 0
        for param in self.parameters:
            col = fits.Column(name='PARAM{0}'.format(i), format='3E', 
                              array=[param.to_fits_value(asymmetric=True)])
            cols.append(col)
            i+=1
        
        cols.append(fits.Column(name='PHTFLUX', format='3E',
                    array=[self.photon_flux.to_fits_value()]))
        cols.append(fits.Column(name='PHTFLNC', format='3E',
                    array=[self.photon_fluence.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLUX', format='3E',
                    array=[self.energy_flux.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLNC', format='3E',
                    array=[self.energy_fluence.to_fits_value()]))
        cols.append(fits.Column(name='REDCHSQ', format='2E',
                    array=[[self.stat_value/self.dof]*2]))
        cols.append(fits.Column(name='CHSQDOF', format='1I', array=[self.dof]))
        cols.append(fits.Column(name='PHTFLUXB', format='3E',
                    array=[self.photon_flux_50_300.to_fits_value()]))
        cols.append(fits.Column(name='DURFLNC', format='3E',
                    array=[self.duration_fluence.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLNCB', format='3E',
                    array=[self.energy_fluence_50_300.to_fits_value()]))
        cols.append(fits.Column(name='COVARMAT', format='{0}E'.format(numparams*numparams), 
                    dim='({0},{0})'.format(numparams), 
                    array=[self.covariance]))

        hdu = fits.BinTableHDU.from_columns(cols, name='FIT PARAMS')
        return hdu
    
    @classmethod
    def from_xspec_data_dict(cls, data_dict):
        """Create a ModelFit object from an XSPEC SCAT data dictionary
        
        Returns:
        -----------
        obj: ModelFit
        """       
        # extract the parameters
        numparams = len(data_dict['parameter_names'])
        params = []
        for i in range(numparams):
            p = Parameter.create(name=data_dict['parameter_names'][i], 
                                 value=data_dict['parameter_values'][i],
                                 uncertainty=data_dict['parameter_sigmas'][i])
            params.append(p)
        
        # extract the fluxes/fluences
        dt = np.abs(data_dict['time_range'][1]-data_dict['time_range'][0])
        eflux = data_dict['energy_flux']
        efluence = eflux*dt
        pflux = data_dict['photon_flux']
        pfluence = pflux*dt
        eflux_50_300 = data_dict['energy_flux_50_300']
        efluence_50_300 = eflux_50_300*dt
        pflux_50_300 = data_dict['photon_flux_50_300']
        
        eflux = Parameter.create(name='energy flux (10-1000 keV)', value=eflux[0], 
                                 uncertainty=eflux[::-1][:2], valid_range=(0., np.inf))
        efluence = Parameter.create(name='energy fluence (10-1000 keV)', value=efluence[0], 
                                 uncertainty=efluence[::-1][:2], valid_range=(0., np.inf))
        pflux = Parameter.create(name='photon flux (10-1000 keV)', value=pflux[0], 
                                 uncertainty=pflux[::-1][:2], valid_range=(0., np.inf))
        pfluence = Parameter.create(name='photon fluence (10-1000 keV)', value=pfluence[0], 
                                 uncertainty=pfluence[::-1][:2], valid_range=(0., np.inf))
        pflux50300 = Parameter.create(name='photon flux (50-300 keV)', value=pflux_50_300[0], 
                                 uncertainty=pflux_50_300[::-1][:2], valid_range=(0., np.inf))
        efluence50300 = Parameter.create(name='energy fluence (50-300 keV)', value=efluence_50_300[0], 
                                 uncertainty=efluence_50_300[::-1][:2], valid_range=(0., np.inf))
        
        # create the object
        obj = cls.create(parameters=params, stat_name=data_dict['stat_name'], 
                         stat_value=data_dict['stat_value'], dof=data_dict['dof'],
                         energy_flux=eflux, photon_flux=pflux, energy_fluence=efluence,
                         photon_fluence=pfluence, photon_flux_50_300=pflux50300,
                         energy_fluence_50_300=efluence50300, 
                         time_range=data_dict['time_range'])
        obj.covariance=data_dict['covar']
        return obj
    
    @classmethod
    def create(cls, **kwargs):
        obj = cls()
        obj._init_by_dict(kwargs)
        return obj
    
    @property
    def covariance(self):
        return self._covariance
    @covariance.setter
    def covariance(self, value):
        try:
            nrow, ncol = value.shape
            if nrow != ncol:
                raise ValueError('ModelFit.covariance must be a square array')
            if nrow != len(self.parameters):
                raise ValueError('ModelFit.covariance wrong size')
        except:
            raise ValueError('ModelFit.covariance must be an array')
        self._covariance = value
    
    def _init_by_dict(self, values):
        for key, val in values.items():
            try:
                p = getattr(self, key)
                if isinstance(p, property):
                    p.__set__(self, val)
                else:
                    self.__setattr__(key, val)
            except AttributeError:
                raise ValueError("{} is not a valid attribute".format(key))
    
class DetectorData(object):
    """A Detector Data class
    
    Attributes:
    -----------
    instrument: str
        The name of the instrument
    detector: str
        The name of the detector
    datatatyp: str
        The name of the datatype
    filename: str
        The filename of the data file
    response: str
        The filename of the detector response
    time_range: (float, float)
        The time range of the data used
    energy_range: (float, float)
        The energy range of the data used
    channel_range = (int, int)
        The energy channel range of the data used
    energy_edges: np.array
        The edges of the energy channels
    photon_counts: np.array
        The deconvolved photon counts for the detector
    photon_model: np.array
        The photon model for the detector
    photon_errors: np.array
        The deconvolved photon count errors for the detector
    
    Public Methods:
    ---------------
    to_fits_row:
        Return the contained data as a FITS table row
    
    Class Methods:
    --------------
    create:
        Construct a ModelFit object from keywords or dictionary
    """
    def __init__(self):
        self.instrument = ''
        self.detector = ''
        self.datatype = ''
        self.filename = ''
        self.response = ''
        self.time_range = (None, None)
        self.energy_range = (None, None)
        self.channel_range = (None, None)
        self.energy_edges = None
        self.photon_counts = None
        self.photon_model = None
        self.photon_errors = None
    
    @classmethod
    def create(cls, **kwargs):
        obj = cls()
        obj._init_by_dict(kwargs)
        return obj
    
    def to_fits_row(self):
        """Return the contained data as a FITS table row
        
        Returns:
        -----------
        hdu: astropy.io.fits.BinTableHDU
            The FITS table
        """
        numchans = len(self.energy_edges)
        e_dim = str(numchans)+'E'
        p_dim = str(numchans-1)+'E'
        p_unit = 'Photon cm^-2 s^-1 keV^-1'
        fit_int = '{0}: {1} s, '.format(self.time_range[0], self.time_range[1])
        fit_int+= '{0}: {1} keV, '.format(self.energy_range[0], self.energy_range[1])
        fit_int+= 'channels {0}: {1}'.format(self.channel_range[0], self.channel_range[1])
        
        col1 = fits.Column(name='INSTRUME', format='20A', array=[self.instrument])
        col2 = fits.Column(name='DETNAM', format='20A', array=[self.detector])
        col3 = fits.Column(name='DATATYPE', format='20A', array=[self.datatype])
        col4 = fits.Column(name='DETSTAT', format='20A', array=['INCLUDED'])
        col5 = fits.Column(name='DATAFILE', format='60A', array=[self.filename])
        col6 = fits.Column(name='RSPFILE', format='60A', array=[self.response])
        col7 = fits.Column(name='FIT_INT', format='60A', array=[fit_int])
        col8 = fits.Column(name='CHANNUM', format='1J', array=[numchans-1])
        col9 = fits.Column(name='FITCHAN', format='2J', array=[self.channel_range])
        col10= fits.Column(name='E_EDGES', format=e_dim, unit='keV', 
                           array=[self.energy_edges])
        col11= fits.Column(name='PHTCNTS', format=p_dim, unit=p_unit, 
                           array=[self.photon_counts])
        col12= fits.Column(name='PHTMODL', format=p_dim, unit=p_unit, 
                           array=[self.photon_model])
        col13= fits.Column(name='PHTERRS', format=p_dim, unit=p_unit, 
                           array=[self.photon_errors])
        hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, 
                                             col7, col8, col9, col10, col11, 
                                             col12, col13], name='DETECTOR DATA')
        return hdu
    
    def _init_by_dict(self, values):
        for key, val in values.items():
            try:
                p = getattr(self, key)
                if isinstance(p, property):
                    p.__set__(self, val)
                else:
                    self.__setattr__(key, val)
            except AttributeError:
                raise ValueError("{} is not a valid attribute".format(key))

class SCat(Data):
    """A container class for the data in an SCAT file
    
    Attributes:
    -----------
    detectors: list
        The DetectorData objects used in the analysis
    model_fits: list
        The ModelFit objects, one for each model fit
    trigtime: float
        The trigger time of the data used
    tstart: float
        The start time of the data used
    tstop: float
        The stop time of the data used
    
    Public Methods:
    ---------------
    add_detector_data:
        Add a new DetectorData object to the SCat
    add_model_fit:
        Add a new fit object to the SCat
    write:
        Write the SCat object to a FITS file
    
    Class Methods:
    --------------
    from_file:
        Create an SCat object from an SCAT FITS file
    """
    def __init__(self, tstart=None, tstop=None, trigtime=None):
        self.detectors = []
        self.model_fits = []
        self.tstart = tstart
        self.tstop = tstop
        self.trigtime = trigtime

        super(SCat, self).__init__(None)
        self.is_gbm_file = True
        self.data_type = 'scat'
        if trigtime is not None:
            self.is_trigger = True
            met = Met(trigtime)
            self.id = met.bn
        else:
            self.is_trigger = False
            if tstart is not None:
                met = Met(tstart)
                self.id = met.ymd_h
        self.create_filename(extension='fit')
    
    def add_detector_data(self, detector_data):
        self.detectors.append(detector_data)
    
    def add_model_fit(self, model_fit):
        self.model_fits.append(model_fit)
    
    def to_table(self):
        used_dets = '+'.join([det.detector for det in self.detectors])
        header = ['Detectors', 'Dt (s)', 'Model']
        params = []
        for model in self.model_fits:
            new_params = model.parameter_list()
            for p in new_params:
                if p not in params:
                    params.append(p)
        header.extend(params)
        header.extend(['Statistic/DoF', 'Ph. Flux', 'En. Flux', 'En. Fluence'])
        numcols = len(header)
        
        table = []
        for model in self.model_fits:
            row = ['']*11
            row[0] = used_dets
            row[1] = ':'.join(model.time_range.astype(str))
            row[2] = model.name
            row[-4] = '{0:.2f}'.format(model.stat_value) + '/' + str(model.dof)
            row[-3] = model.photon_flux.to_table_value()
            row[-2] = model.energy_flux.to_table_value()
            row[-1] = model.energy_fluence.to_table_value()
            
            model_params = model.parameter_list()
            for i in range(len(model_params)):
                idx = header.index(model_params[i])
                row[idx] = model.parameters[i].to_table_value()
            table.append(row)
            
        
        sys.exit()
    
    @classmethod
    def from_file(cls, filename):
        """Create an SCat object from an SCAT FITS file
        
        Returns:
        -----------
        obj: SCat
            The SCat object
        """
        obj = cls()
        Data.__init__(obj,filename)
        
        with fits.open(filename) as hdulist:
            obj.tstart = hdulist[0].header['TSTART']
            obj.tstop = hdulist[0].header['TSTOP']
            obj.trigtime = hdulist[0].header['TRIGTIME']
            
            # read the detector data HDU
            det_data = hdulist['DETECTOR DATA'].data
            obj._from_detector_hdu(det_data)
            
            # read the fit params HDU
            fit_hdu = hdulist['FIT PARAMS']
            obj._from_fitparam_hdu(fit_hdu)
                
        return obj
    
    def write(self, directory, filename=None):
        """Write an SCat object to a FITS file
        
        Parameters:
        -----------
        directory: str
            The directory where the file is to be written
        filename: str, optional
            The filename.  If not set, a default filename will be chosen.
        """
        if (self.filename is None) and (filename is None):
            raise NameError('Filename not set')
        if filename is None:
            filename = self.filename_obj.basename()
        self.set_filename(filename, directory=directory)
        
        # initialize the FITS file
        hdulist = fits.HDUList()
        prihdr = self._primary_header()
        primary_hdu = fits.PrimaryHDU(header=prihdr)
        primary_hdu.add_checksum()
        hdulist.append(primary_hdu)
        
        # construct the detector data extension
        det_hdu = None
        for detector in self.detectors:
            if det_hdu is None:
                det_hdu = detector.to_fits_row()
            else:
                det_hdu.data = np.concatenate((det_hdu.data, detector.to_fits_row().data))
        det_hdu.header = self._update_detector_hdr(det_hdu.header)
        det_hdu.add_checksum()
        hdulist.append(det_hdu)
        
        # construct the fit extension
        fit_hdu = None
        for fit in self.model_fits:
            if fit_hdu is None:
                fit_hdu = fit.to_fits_row()
            else:
                fit_hdu.data = np.concatenate((fit_hdu.data, fit.to_fits_row().data))
        fit_hdu.header = self._update_fitparam_hdr(fit_hdu.header)
        fit_hdu.add_checksum()
        hdulist.append(fit_hdu)
        
        # write out the file
        filename = directory+filename
        hdulist.writeto(self.filename, checksum=True, clobber=True)
    
    def _from_detector_hdu(self, det_hdu):
        for row in det_hdu:
            # have to split out the time range and energy range from the
            # FIT_INT column
            t_range, e_range, _ = row['FIT_INT'].split(',')
            t_range = t_range.split(':')
            t_range = (float(t_range[0]), float(t_range[1].split()[0]))
            e_range = e_range.split(':')
            e_range = (float(e_range[0]), float(e_range[1].split()[0]))
            c_range = tuple(row['FITCHAN'])
            # create a Detector Data object
            det = DetectorData.create(instrument=row['INSTRUME'], detector=row['DETNAM'],
                               datatype=row['DATATYPE'], filename=row['DATAFILE'],
                               response=row['RSPFILE'], time_range=t_range,
                               energy_range=e_range, channel_range=c_range,
                               energy_edges=row['E_EDGES'], photon_counts=row['PHTCNTS'],
                               photon_model=row['PHTMODL'], photon_errors=row['PHTERRS'])
            self.add_detector_data(det)
        
    def _from_fitparam_hdu(self, fit_hdu):
            
            fit_data = fit_hdu.data
            fit_hdr = fit_hdu.header
            e_range = (fit_hdr['FLU_LOW'], fit_hdr['FLU_HIGH'])
            nparams = fit_hdr['N_PARAM']
            stat_name = fit_hdr['STATISTC']
            # find unique models in the event there are multiple components
            models = [fit_hdr.comments['TTYPE'+str(2+i)].split(':')[0] 
                      for i in range(nparams)]
            model = '+'.join(list(set(models)))
            # the parameter names
            param_names = [fit_hdr.comments['TTYPE'+str(2+i)].split(':')[1].strip() 
                           for i in range(nparams)]
            
            for row in fit_data:
                # the array of parameters
                params = []
                for i in range(nparams):
                    param = row['PARAM'+str(i)]
                    params.append(Parameter.create(name=param_names[i], 
                                                   value=param[0], 
                                                   uncertainty=param[1:]))
                # fluxes and fluences
                pflux = Parameter.create(name='photon flux', value=row['PHTFLUX'][0],
                                         uncertainty=row['PHTFLUX'][1], 
                                         valid_range=(0.0, np.inf))
                pflnc = Parameter.create(name='photon fluence', value=row['PHTFLNC'][0],
                                         uncertainty=row['PHTFLNC'][1], 
                                         valid_range=(0.0, np.inf))
                eflux = Parameter.create(name='energy flux', value=row['NRGFLUX'][0],
                                         uncertainty=row['NRGFLUX'][1], 
                                         valid_range=(0.0, np.inf))
                eflnc = Parameter.create(name='energy fluence', value=row['NRGFLNC'][0],
                                         uncertainty=row['NRGFLNC'][1], 
                                         valid_range=(0.0, np.inf))
                pfluxb = Parameter.create(name='photon flux', value=row['PHTFLUXB'][0],
                                         uncertainty=row['PHTFLUXB'][1], 
                                         valid_range=(0.0, np.inf))
                durflnc = Parameter.create(name='photon fluence', value=row['DURFLNC'][0],
                                         uncertainty=row['DURFLNC'][1], 
                                         valid_range=(0.0, np.inf))
                eflncb = Parameter.create(name='energy fluence', value=row['NRGFLNCB'][0],
                                         uncertainty=row['NRGFLNCB'][1], 
                                         valid_range=(0.0, np.inf))
                # create a ModelFit object
                fit = ModelFit.create(name=model, flux_energy_range=e_range,
                                      time_range=row['TIMEBIN'], parameters=params,
                                      photon_flux=pflux, energy_flux=eflux,
                                      photon_fluence=pflnc, energy_fluence=eflnc,
                                      photon_flux_50_300=pfluxb, duration_fluence=durflnc,
                                      energy_fluence_50_300=eflncb, stat_name=stat_name,
                                      stat_value=row['REDCHSQ'][1]*row['CHSQDOF'],
                                      dof=row['CHSQDOF'])
                fit.covariance=row['COVARMAT']
                self.add_model_fit(fit)
    
    def _update_detector_hdr(self, det_hdr):
        det_hdr.comments['TTYPE1'] = 'Instrument name for this detector'
        det_hdr.comments['TTYPE2'] = 'Detector number; if one of several available'
        det_hdr.comments['TTYPE3'] = 'Data type used for this analysis'
        det_hdr.comments['TTYPE4'] = 'Was this detector INCLUDED or OMITTED'
        det_hdr.comments['TTYPE5'] = 'Data file name for this dataset'
        det_hdr.comments['TTYPE6'] = 'Response file name for this dataset'
        det_hdr.comments['TTYPE7'] = 'Fit intervals'
        det_hdr.comments['TTYPE8'] = 'Total number of energy channels for this detector'
        det_hdr.comments['TTYPE9'] = 'Channels selected in fitting this detector'
        det_hdr.comments['TTYPE10'] = 'Energy edges for each selected detector'
        det_hdr.comments['TTYPE11'] = 'Array of photon counts data'
        det_hdr.comments['TTYPE12'] = 'Array of photon model data'
        det_hdr.comments['TTYPE13'] = 'Array of errors in photon counts data'
        det_hdr['NUMFITS'] = (len(self.model_fits), 'Number of spectral fits in the data')
        prihdu = self._primary_header()
        keys = ['ORIGIN', 'TELESCOP', 'INSTRUME', 'OBSERVER', 'MJDREFI', 'MJDREFF',
                'TIMESYS', 'TIMEUNIT', 'DATE-OBS', 'DATE-END', 'TSTART', 'TSTOP', 
                'TRIGTIME']
        for key in keys:
            det_hdr[key] = (prihdu[key], prihdu.comments[key])
        return det_hdr
    
    def _update_fitparam_hdr(self, fit_hdr):
        e_range = self.model_fits[0].flux_energy_range
        model_name = self.model_fits[0].name
        param_names = self.model_fits[0].parameter_list()
        numparams = len(param_names)
        statistic = self.model_fits[0].stat_name
        
        g_range = '({0}-{1} keV)'.format(e_range[0], e_range[1])
        b_range = '(50-300 keV)'
        
        fit_hdr.comments['TTYPE1'] = 'Start and stop times relative to trigger'
        for i in range(numparams):
            colname = 'TTYPE{0}'.format(i+2)
            fit_hdr.comments[colname] = '{0}: {1}'.format(model_name, param_names[i])
        ttypes = ['TTYPE'+str(numparams+2+i) for i in range(10)]
        fit_hdr.comments[ttypes[0]] = 'Photon Flux (ph/s-cm^2) std energy '+g_range
        fit_hdr.comments[ttypes[1]] = 'Photon Fluence (ph/cm^2) std energy '+g_range
        fit_hdr.comments[ttypes[2]] = 'Energy Flux (erg/s-cm^2) std energy '+g_range
        fit_hdr.comments[ttypes[3]] = 'Energy Fluence (erg/cm^2) std energy '+g_range
        fit_hdr.comments[ttypes[4]] = 'Reduced Chi^2 (1) and fitting statistic (2)'
        fit_hdr.comments[ttypes[5]] = 'Degrees of Freedom'
        fit_hdr.comments[ttypes[6]] = 'Photon Flux (ph/s-cm^2) BATSE energy '+b_range
        fit_hdr.comments[ttypes[7]] = 'Photon Fluence (ph/cm^2) for durations (user)'
        fit_hdr.comments[ttypes[8]] = 'Energy Fluence (erg/cm^2) BATSE energy '+b_range
        fit_hdr.comments[ttypes[9]] = 'Covariance matrix for the fir (N_PARAM^2)'
        fit_hdr['N_PARAM'] = (numparams, 'Total number of fit parameters (PARAMn)')
        fit_hdr['FLU_LOW'] = (e_range[0], 'Lower limit of flux/fluence integration (keV)')
        fit_hdr['FLU_HIGH'] = (e_range[1], 'Upeer limit of flux/fluence integration (keV)')
        fit_hdr['STATISTC'] = (statistic, 'Indicates merit function used for fitting')
        fit_hdr['NUMFITS'] = (len(self.model_fits), 'Number of spectral fits in the data')
        
        prihdu = self._primary_header()
        keys = ['ORIGIN', 'TELESCOP', 'INSTRUME', 'OBSERVER', 'MJDREFI', 'MJDREFF',
                'TIMESYS', 'TIMEUNIT', 'DATE-OBS', 'DATE-END', 'TSTART', 'TSTOP', 
                'TRIGTIME']
        for key in keys:
            fit_hdr[key] = (prihdu[key], prihdu.comments[key])
        return fit_hdr
    
    def _primary_header(self):
        # create standard GBM primary header
        filetype = 'SPECTRAL FITS'
        prihdr = hdr.primary(filetype=filetype, tstart=self.tstart, 
                             filename=self.filename_obj.basename(), 
                             tstop=self.tstop, trigtime=self.trigtime)
        # remove the keywords we don't need
        del prihdr['DETNAM'], prihdr['OBJECT'], prihdr['RA_OBJ'], \
            prihdr['DEC_OBJ'], prihdr['ERR_RAD'], prihdr['INFILE01']
        
        return prihdr