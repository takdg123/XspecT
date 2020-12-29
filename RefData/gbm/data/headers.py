import astropy.io.fits as fits
from gbm.time import Met
import gbm

class GBMDefinitions(object):
    def __init__(self):
        self.telescope = ('telescop', 'GLAST', 'Name of mission/satellite')
        self.instrument = ('instrume', 'GBM', 'Specific instrument used for observation')
        self.observer = ('observer', 'Meegan', 'GLAST Burst Monitor P.I.')
        self.origin = ('origin', 'GIOC', 'Name of organization making file')
        self.timesys = ('timesys', 'TT', 'Time system used in time keywords')
        self.timeunit = ('timeunit', 's', 'Time since MJDREF, used in TSTART and TSTOP')
        self.mjdrefi = ('mjdrefi', 51910, 'MJD of GLAST reference epoch, integer part')
        self.mjdreff = ('mjdreff', '7.428703703703703e-4', 'MJD of GLAST reference epoch, fractional part')
        self.radecsys = ('radecsys', 'FK5', 'Stellar reference frame')
        self.equinox = ('equinox', 2000.0, 'Equinox for RA and Dec')

def primary(detnam=None, filetype=None, tstart=None, tstop=None, filename=None, 
            trigtime=None, object=None, ra_obj=None, dec_obj=None, err_rad=None,
            infile=None):
    
    # enforce strings
    if detnam is not None:
        detnam = str(detnam)
    if filetype is not None:
        filetype = str(filetype)
    if filename is not None:
        filename = str(filename)
    if object is not None:
        object = str(object)
    if infile is not None:
        infile = str(infile)
    
    # enforce floats
    if tstart is not None:
        tstart = float(tstart)
    if tstop is not None:
        tstop = float(tstop)
    if trigtime is not None:
        trigtime = float(trigtime)
    if ra_obj is not None:
        ra_obj = float(ra_obj)
    if dec_obj is not None:
        dec_obj = float(dec_obj)
    if err_rad is not None:
        err_rad = float(err_rad)
    
    # do MET -> UTC time conversion
    date_obs = None
    if tstart is not None:
        met = Met(tstart)
        date_obs = met.iso()
    date_end = None
    if tstop is not None:
        met = Met(tstop)
        date_end = met.iso()        
    
    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(('creator', 'GBM Data Tools {:} Software and version creating file'.format(gbm.__version__)))
    header.append(('filetype', filetype, 'Name for this type of FITS file'))
    header.append(('file-ver', '1.0.0', 'Version of the format for this filetype'))    
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(('detnam', detnam, 'Individual detector name'))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(('date', current_time, 'file creation date (YYYY-MM-DDThh:mm:ss UT)'))
    header.append(('date-obs', date_obs, 'Date of start of observation')) 
    header.append(('date-end', date_end, 'Date of end of observation'))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(('tstart', tstart, '[GLAST MET] Observation start time'))
    header.append(('tstop', tstop, '[GLAST MET] Observation stop time'))
    header.append(('filename', filename, 'Name of this file'))
    header.append(('trigtime', trigtime, 'Trigger time relative to MJDREF, double precision'))
    header.append(('object', object, 'Burst name in standard format, yymmddfff'))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(('ra_obj', ra_obj, 'Calculated RA of burst'))
    header.append(('dec_obj', dec_obj, 'Calculated Dec of burst'))
    header.append(('err_rad', err_rad, 'Calculated Location Error Radius'))
    header.append(('infile01', infile, 'Level 1 input lookup table file'))
    
    return header


def ebounds(detnam=None, tstart=None, tstop=None, trigtime=None, object=None, 
            ra_obj=None, dec_obj=None, err_rad=None, detchans=None, ch2e_ver=None,
            gain=None, infile=None):

    # enforce strings
    if detnam is not None:
        detnam = str(detnam)
    if object is not None:
        object = str(object)
    if ch2e_ver is not None:
        ch2e_ver = str(ch2e_ver)
    if infile is not None:
        infile = str(infile)
        
    # enforce floats
    if tstart is not None:
        tstart = float(tstart)
    if tstop is not None:
        tstop = float(tstop)
    if trigtime is not None:
        trigtime = float(trigtime)
    if ra_obj is not None:
        ra_obj = float(ra_obj)
    if dec_obj is not None:
        dec_obj = float(dec_obj)
    if err_rad is not None:
        err_rad = float(err_rad)
    if gain is not None:
        gain = float(gain)
    
    if detchans is not None:
        detchans = int(detchans)
    
    # do MET -> UTC time conversion
    date_obs = None
    if tstart is not None:
        met = Met(tstart)
        date_obs = met.iso()
    date_end = None
    if tstop is not None:
        met = Met(tstop)
        date_end = met.iso()        

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append('extname', 'EBOUNDS', 'name of this binary table extension')
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(('detnam', detnam, 'Individual detector name'))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(('date', current_time, 'file creation date (YYYY-MM-DDThh:mm:ss UT)'))
    header.append(('date-obs', date_obs, 'Date of start of observation')) 
    header.append(('date-end', date_end, 'Date of end of observation'))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(('tstart', tstart, '[GLAST MET] Observation start time'))
    header.append(('tstop', tstop, '[GLAST MET] Observation stop time'))
    header.append(('trigtime', trigtime, 'Trigger time relative to MJDREF, double precision'))
    header.append(('object', object, 'Burst name in standard format, yymmddfff'))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(('ra_obj', ra_obj, 'Calculated RA of burst'))
    header.append(('dec_obj', dec_obj, 'Calculated Dec of burst'))
    header.append(('err_rad', err_rad, 'Calculated Location Error Radius'))
    header.append(('hduclass', 'OGIP', 'Conforms to OGIP standard indicated in HDUCLAS1'))
    header.append(('hduclas1', 'RESPONSE', 'These are typically found in RMF files'))
    header.append(('hduclas2', 'EBOUNDS', 'From CAL/GEN/92-002'))
    header.append(('hduvers', '1.2.0', 'Version of HDUCLAS1 format in use'))
    header.append(('chantype', 'PHA', 'No corrections have been applied'))
    header.append(('filter', 'None', 'The instrument filter in use (if any)'))
    header.append(('detchans', detchans, 'Total number of channels in each rate'))
    header.append(('extver', 1, 'Version of this extension format'))
    header.append(('ch2e_ver', ch2e_ver, 'Channel to energy conversion scheme used'))
    header.append(('gain_cor', gain, 'Gain correction factor applied to energy edges'))
    header.append(('infile01', infile, 'Level 1 input lookup table file'))
    
    return header

def spectrum(detnam=None, tstart=None, tstop=None, trigtime=None, object=None, 
            ra_obj=None, dec_obj=None, err_rad=None, detchans=None, poisserr=True):

    # enforce strings
    if detnam is not None:
        detnam = str(detnam)
    if object is not None:
        object = str(object)
        
    # enforce floats
    if tstart is not None:
        tstart = float(tstart)
    if tstop is not None:
        tstop = float(tstop)
    if trigtime is not None:
        trigtime = float(trigtime)
    if ra_obj is not None:
        ra_obj = float(ra_obj)
    if dec_obj is not None:
        dec_obj = float(dec_obj)
    if err_rad is not None:
        err_rad = float(err_rad)
    
    if detchans is not None:
        detchans = int(detchans)
    
    # do MET -> UTC time conversion
    date_obs = None
    if tstart is not None:
        met = Met(tstart)
        date_obs = met.iso()
    date_end = None
    if tstop is not None:
        met = Met(tstop)
        date_end = met.iso()        

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(('tzero4', trigtime, 'Offset, equal to TRIGTIME'))
    header.append(('tzero5', trigtime, 'Offset, equal to TRIGTIME'))
    header.append('extname', 'SPECTRUM', 'name of this binary table extension')
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(('detnam', detnam, 'Individual detector name'))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(('date', current_time, 'file creation date (YYYY-MM-DDThh:mm:ss UT)'))
    header.append(('date-obs', date_obs, 'Date of start of observation')) 
    header.append(('date-end', date_end, 'Date of end of observation'))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(('tstart', tstart, '[GLAST MET] Observation start time'))
    header.append(('tstop', tstop, '[GLAST MET] Observation stop time'))
    header.append(('trigtime', trigtime, 'Trigger time relative to MJDREF, double precision'))
    header.append(('object', object, 'Burst name in standard format, yymmddfff'))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(('ra_obj', ra_obj, 'Calculated RA of burst'))
    header.append(('dec_obj', dec_obj, 'Calculated Dec of burst'))
    header.append(('err_rad', err_rad, 'Calculated Location Error Radius'))
    header.append(('filter', 'None', 'The instrument filter in use (if any)'))
    header.append(('areascal', 1., 'No special scaling of effective area by channel'))
    header.append(('backfil', 'none', 'Name of corresponding background file (if any)'))
    header.append(('backscal', 1., 'No scaling of background'))
    header.append(('corrfile', 'none', 'Name of corresponding correction file (if any)'))
    header.append(('corrscal', 1., 'Correction scaling file'))
    header.append(('respfile', 'none', 'Name of corresponding RMF file (if any)'))
    header.append(('ancrfile', 'none', 'Name of corresponding ARF file (if any)'))
    header.append(('sys_err', 0., 'No systematic errors'))
    header.append(('poisserr', poisserr, 'Assume Poisson Errors'))
    header.append(('grouping', 0, 'No special grouping has been applied'))    
    header.append(('hduclass', 'OGIP', 'Format conforms to OGIP standard'))
    header.append(('hduclas1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'))
    header.append(('hduclas2', 'TOTAL', 'Indicates gross data (source + background)'))
    header.append(('hduclas3', 'COUNT', 'Indicates data stored as counts'))
    header.append(('hduclas4', 'TYPEII', 'Indicates PHA Type II file format'))    
    header.append(('hduvers', '1.2.1', 'Version of HDUCLAS1 format in use'))
    header.append(('chantype', 'PHA', 'No corrections have been applied'))
    header.append(('detchans', detchans, 'Total number of channels in each rate'))
    header.append(('extver', 1, 'Version of this extension format'))
    
    return header

def events(detnam=None, tstart=None, tstop=None, trigtime=None, object=None, 
            ra_obj=None, dec_obj=None, err_rad=None, detchans=None):

    # enforce strings
    if detnam is not None:
        detnam = str(detnam)
    if object is not None:
        object = str(object)
        
    # enforce floats
    if tstart is not None:
        tstart = float(tstart)
    if tstop is not None:
        tstop = float(tstop)
    if trigtime is not None:
        trigtime = float(trigtime)
    if ra_obj is not None:
        ra_obj = float(ra_obj)
    if dec_obj is not None:
        dec_obj = float(dec_obj)
    if err_rad is not None:
        err_rad = float(err_rad)
    
    if detchans is not None:
        detchans = int(detchans)
    
    # do MET -> UTC time conversion
    date_obs = None
    if tstart is not None:
        met = Met(tstart)
        date_obs = met.iso()
    date_end = None
    if tstop is not None:
        met = Met(tstop)
        date_end = met.iso()        

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(('tzero1', trigtime, 'Offset, equal to TRIGTIME'))
    header.append('extname', 'EVENTS', 'name of this binary table extension')
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(('detnam', detnam, 'Individual detector name'))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(('date', current_time, 'file creation date (YYYY-MM-DDThh:mm:ss UT)'))
    header.append(('date-obs', date_obs, 'Date of start of observation')) 
    header.append(('date-end', date_end, 'Date of end of observation'))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(('tstart', tstart, '[GLAST MET] Observation start time'))
    header.append(('tstop', tstop, '[GLAST MET] Observation stop time'))
    header.append(('trigtime', trigtime, 'Trigger time relative to MJDREF, double precision'))
    header.append(('object', object, 'Burst name in standard format, yymmddfff'))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(('ra_obj', ra_obj, 'Calculated RA of burst'))
    header.append(('dec_obj', dec_obj, 'Calculated Dec of burst'))
    header.append(('err_rad', err_rad, 'Calculated Location Error Radius'))
    header.append(('respfile', 'none', 'Name of corresponding RMF file (if any)'))
    header.append(('evt_dead', 2.6e-6, 'Deadtime per event (s)'))
    header.append(('detchans', detchans, 'Total number of channels in each rate'))
    header.append(('hduclass', 'OGIP', 'Format conforms to OGIP standard'))
    header.append(('hduclas1', 'EVENTS', 'Extension contains Events'))
    header.append(('extver', 1, 'Version of this extension format'))
        
    return header

def gti(detnam=None, tstart=None, tstop=None, trigtime=None, object=None, 
            ra_obj=None, dec_obj=None, err_rad=None):

    # enforce strings
    if detnam is not None:
        detnam = str(detnam)
    if object is not None:
        object = str(object)
        
    # enforce floats
    if tstart is not None:
        tstart = float(tstart)
    if tstop is not None:
        tstop = float(tstop)
    if trigtime is not None:
        trigtime = float(trigtime)
    if ra_obj is not None:
        ra_obj = float(ra_obj)
    if dec_obj is not None:
        dec_obj = float(dec_obj)
    if err_rad is not None:
        err_rad = float(err_rad)
        
    # do MET -> UTC time conversion
    date_obs = None
    if tstart is not None:
        met = Met(tstart)
        date_obs = met.iso()
    date_end = None
    if tstop is not None:
        met = Met(tstop)
        date_end = met.iso()        

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append(('tzero1', trigtime, 'Offset, equal to TRIGTIME'))
    header.append(('tzero2', trigtime, 'Offset, equal to TRIGTIME'))
    header.append('extname', 'GTI', 'name of this binary table extension')
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(('detnam', detnam, 'Individual detector name'))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(('date', current_time, 'file creation date (YYYY-MM-DDThh:mm:ss UT)'))
    header.append(('date-obs', date_obs, 'Date of start of observation')) 
    header.append(('date-end', date_end, 'Date of end of observation'))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(('tstart', tstart, '[GLAST MET] Observation start time'))
    header.append(('tstop', tstop, '[GLAST MET] Observation stop time'))
    header.append(('hduclass', 'OGIP', 'Format conforms to OGIP standard'))
    header.append(('hduclas1', 'GTI', 'Indicates good time intervals'))
    header.append(('hduvers', '1.2.0', 'Version of HDUCLAS1 format in use'))
    header.append(('extver', 1, 'Version of this extension format'))
    header.append(('trigtime', trigtime, 'Trigger time relative to MJDREF, double precision'))
    header.append(('object', object, 'Burst name in standard format, yymmddfff'))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(('ra_obj', ra_obj, 'Calculated RA of burst'))
    header.append(('dec_obj', dec_obj, 'Calculated Dec of burst'))
    header.append(('err_rad', err_rad, 'Calculated Location Error Radius'))
        
    return header

def specresp(detnam=None, tstart=None, tstop=None, trigtime=None, object=None, 
            ra_obj=None, dec_obj=None, mat_type=None, rsp_num=None, src_az=None, 
            src_el=None, geo_az=None, geo_el=None, det_ang=None, geo_ang=None, 
            numebins=None, detchans=None, infiles=None, atscat=None):

    # enforce strings
    if detnam is not None:
        detnam = str(detnam)
    if object is not None:
        object = str(object)
    if mat_type is not None:
        mat_type = str(mat_type)
    if infiles is None:
        infiles = [None]*3
    else:
        infiles = [str(infile) for infile in infiles]
    if atscat is not None:
        atscat = str(atscat)
        
    # enforce floats
    if tstart is not None:
        tstart = float(tstart)
    if tstop is not None:
        tstop = float(tstop)
    if trigtime is not None:
        trigtime = float(trigtime)
    if ra_obj is not None:
        ra_obj = float(ra_obj)
    if dec_obj is not None:
        dec_obj = float(dec_obj)
    if src_az is not None:
        src_az = float(src_az)
    if src_el is not None:
        src_el = float(src_el)
    if geo_az is not None:
        geo_az = float(geo_az)
    if geo_el is not None:
        geo_el = float(geo_el)
    if det_ang is not None:
        det_ang = float(det_ang)
    if geo_ang is not None:
        geo_ang = float(geo_ang)
    
    # enforce integers
    if detchans is not None:
        detchans = int(detchans)
    if rsp_num is not None:
        rsp_num = int(rsp_num)
    if numebins is not None:
        numebins = int(numebins)
    
    # do MET -> UTC time conversion
    date_obs = None
    if tstart is not None:
        met = Met(tstart)
        date_obs = met.iso()
    date_end = None
    if tstop is not None:
        met = Met(tstop)
        date_end = met.iso()        

    header = fits.Header()
    gbmdefs = GBMDefinitions()
    current_time = Met.now().iso()
    header.append('extname', 'SPECRESP MATRIX', 'name of this binary table extension')
    header.append(('extver', 1, 'Version of this extension format'))
    header.append(('date', current_time, 'file creation date (YYYY-MM-DDThh:mm:ss UT)'))
    header.append(('date-obs', date_obs, 'Date of start of observation')) 
    header.append(('date-end', date_end, 'Date of end of observation'))
    header.append(gbmdefs.mjdrefi)
    header.append(gbmdefs.mjdreff)
    header.append(('tstart', tstart, '[GLAST MET] Observation start time'))
    header.append(('tstop', tstop, '[GLAST MET] Observation stop time'))
    header.append(('trigtime', trigtime, 'Trigger time relative to MJDREF, double precision'))
    header.append(gbmdefs.timesys)
    header.append(gbmdefs.timeunit)
    header.append(gbmdefs.telescope)
    header.append(gbmdefs.instrument)
    header.append(('detnam', detnam, 'Individual detector name'))
    header.append(gbmdefs.observer)
    header.append(gbmdefs.origin)
    header.append(('mat_type', mat_type, 'Response Matrix Type'))
    header.append(('rsp_num', rsp_num, 'Response matrix index number'))
    header.append(('object', object, 'Burst name in standard format, yymmddfff'))
    header.append(gbmdefs.radecsys)
    header.append(gbmdefs.equinox)
    header.append(('ra_obj', ra_obj, 'Calculated RA of burst'))
    header.append(('dec_obj', dec_obj, 'Calculated Dec of burst'))
    header.append(('src_az', src_az, 'Azimuth of source in spacecraft coordinates'))
    header.append(('src_el', src_el, 'Elevation of source in spacecraft coordinates'))
    header.append(('geo_az', geo_az, 'Azimuth of geocenter in spacecraft coordinates'))
    header.append(('geo_el', geo_el, 'Elevation of geocenter in spacecraft coordinates'))
    header.append(('det_ang', det_ang, 'Angle between source and detector normal'))
    header.append(('geo_ang', geo_ang, 'Angle between geocenter and detector normal'))
    header.append(('filter', 'none', 'Not applicable to GBM'))
    header.append(('chantype', 'PHA', 'No corrections have been applied'))
    header.append(('numebins', numebins, 'Number of true energy bins of the MATRIX'))
    header.append(('detchans', detchans, 'Total number of channels in each rate'))
    header.append(('infile01', infiles[0], 'Detector response database in'))
    header.append(('infile02', infiles[1], 'Detector response database in'))
    header.append(('infile03', infiles[2], 'Detector response database in'))
    header.append(('infile04', atscat, 'Atmospheric scattering datab'))
    header.append(('hduclass', 'OGIP', 'Format conforms to OGIP standard'))
    header.append(('hduvers', '1.3.0', 'Version of HDUCLAS1 format in use'))
    header.append(('hduclas1', 'Typically found in RMF files'))
    header.append(('hduclas2', 'RSP_MATRIX', 'From CAL/GEN/92-002'))

    return header
