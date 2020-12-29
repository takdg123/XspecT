import astropy.io.fits as fits
import numpy as np
from collections import OrderedDict
from warnings import warn
from gbm import coords
from .data import Data

class PosHist(Data):
    """Container class for position history data

    Attributes:
    -----------
    angular_velocity: np.array
    	The angular velocity of the spacecraft
    in_saa: np.array
        State flag for Fermi in/out of SAA
    lat_long: np.array
    	The latitude and longitude position of the spacecraft
    position: np.array
    	The Geocentric position of the spacecraft
    quaternion: np.array
    	The spacecraft attitude quaternions
    solar_ny: np.array
    	The negative y solar panel pointings
    solar_py: np.array
    	The positive y solar panel pointings
    sun_visible: np.array
        State flag for sun visibility
    time: np.array
    	The times at which the position history is sampled
    velocity: np.array
    	The velocity of the spacecraft

    Public Methods:
    ---------------
    angular_velocity_from_quaternion:
        Estimate the angular velocity vector from the change in quaternions
    lat_long_from_geocentric:
        Estimate the latitude/longitude from the spacecraft position and time
    velocity_from_scpos:
        Estimate the velocity vector from the change in spacecraft position


    Class Methods:
    --------------
    from_trigdat:
        Create a PosHist object from data contained in a TrigDat
    """    
    def __init__(self, filename):
        
        self.header = OrderedDict()
        self.angular_velocity = None
        self.in_saa = None
        self.lat_long = None
        self.position = None
        self.quaternion = None
        self.solar_ny = None
        self.solar_py = None
        self.sun_visible = None
        self.time = None
        self.velocity = None
        
        super(PosHist, self).__init__(filename)
        if self.filename is not None:
            self._from_file()
    		
    def _from_file(self):
        """Strips data from the file, and stores the header.

        Extracts the GLAST POS HIST extension and stores the contained info as attributes
        """
        with fits.open(self.filename) as hdulist:
            for hdu in hdulist:
                self.header.update({hdu.name: hdu.header})
            data = hdulist['GLAST POS HIST'].data
        
        self.time = data['SCLK_UTC']
        self.quaternion = np.array((data['QSJ_1'], data['QSJ_2'], data['QSJ_3'],
                                     data['QSJ_4']))
        self.angular_velocity = np.array((data['WSJ_1'], data['WSJ_2'], data['WSJ_3']))
        self.position = np.array((data['POS_X'], data['POS_Y'], data['POS_Z']))
        self.velocity = np.array((data['VEL_X'], data['VEL_Y'], data['VEL_Z']))
        try:
            self.lat_long = np.array((data['SC_LAT'], data['SC_LON']))
        except:
            warn('No Lat/Lon Information')
        try:
            self.solar_py = data['SADA_PY']
            self.solar_ny = data['SADA_NY']
        except:
            warn('No Solar Panel Pointing Information')
        
        flags = data['FLAGS'] 
        self.sun_visible = flags & 0x01
        self.in_saa = flags & 0x02
    
    @classmethod
    def from_trigdat(cls, time, quaternion, scpos):
        """Constructs poshist container from trigdat position information
		
		Parameters:
		-----------
		time: np.array
			time in MET of the sc positions
		quaternion: np.array
			quaternions; one for each time
		scpos: np.array
			spacecraft positions; one for each time

		Returns:
		-----------
        obj: PosHist
            The created poshist object
        """
        obj = cls(None)
        obj.time = time
        obj.quaternion = quaternion.T
        obj.position = np.transpose(scpos) # in km
        # calculate angular velocity vectors from time and quaternions
        obj.angular_velocity = obj.angular_velocity_from_quaternion(time, 
                                                                    quaternion)
        # calculate velocity vectors from positions and time
        obj.velocity = obj.velocity_from_scpos(obj.time, obj.position)
        
        # calculate the latitude and longitude position of the spacecraft
        lat_long = []
        for i in range(len(time)):
            lat, long, _ = obj.lat_long_from_geocentric(obj.position[:,i], time[i])
            lat_long.append([lat, long])
        obj.lat_long = np.array(lat_long).T
        
        return obj
     
    def angular_velocity_from_quaternion(self, times, quaternions):
        """Calculate angular velocity from changes in quaternion over time.
        This is accurate for small rotations over very short times.
		
		Parameters:
		-----------
		times: np.array
			times in MET of the quaternions
		quaternions: np.array
			quaternions; one for each time
		
		Returns:
		-----------
		angular_velocity: np.array
			The angular velocity
        """
        dt = times[1:]-times[0:-1]        
        dquat = quaternions[:,1:]-quaternions[:,:-1]
        prod_quat = dquat
        for i in range(len(dt)):
            conj_quat = coords.quaternion_conj(quaternions[:,i])
            prod_quat[:,i] = coords.quaternion_prod(conj_quat, dquat[:,i])
        angular_velocity = 2.0/dt[np.newaxis,:]*prod_quat[0:3,:]
        angular_velocity=np.append(angular_velocity, angular_velocity[:,-1:], axis=1)
        return angular_velocity

    def velocity_from_scpos(self, times, positions):
        """Calculate velocity from changes in position over time.
		
		Parameters:
		-----------
		times: np.array
			times in MET of the quaternions
		positions: np.array
			Geocentric positions of the spacecraft
		
		Returns:
		-----------
		velocity: np.array
			The velocity
        """
        dt = times[1:]-times[0:-1]
        dpos = positions[:,1:]-positions[:,0:-1]
        velocity = dpos/dt
        # copy last entry
        velocity = np.append(velocity, velocity[:,-1:], axis=1)
        return velocity
    
    def lat_long_from_geocentric(self, iposition, met, complex=False):
        """Convert Geocentric Cartesian coordinates of the spacecraft to latitude, 
		East longitude, and altitude.
		
		Parameters:
		-----------
		iposition: np.array
			Geocentric cartesian coordinates of the spacecraft location in km
		met: float
			The Fermi MET
		complex: bool, optional
			If True, then uses the non-spherical Earth model and accurate
			Earth rotation model.  Default is False; Earth is a sphere and 
			the Earth rotation model is less accurate.  GBM FSW uses the 
			latter option.
		Returns:
		-----------
		float
			latitude
		float
			East longitude
		float
			altitude
        """
        if complex is False:
            latitude, altitude = coords.latitude_from_geocentric_coords_simple(iposition)
        else:
            latitude, altitude = coords.latitude_from_geocentric_coords_complex(iposition)
        
        longitude = coords.longitude_from_geocentric_coords(iposition, met, ut1=complex)
        
        return (latitude, longitude, altitude)
