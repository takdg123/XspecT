from gbm.data.containers import RateHisto, EventList
from gbm.data.phaii import Cspec, Ctime, TTE
from gbm.data.drm import RSP
from .selection import Selection
from .lookup_old import RmfitLookup, GspecLookup
from .background.background import Background
from .background.binned import polynomial
from .binning.binned import rebin_by_edge_index
from .binning.unbinned import bin_by_edges
import astropy.io.fits as fits
import numpy as np
import os
import sys
from gbm.detectors import Detector
from gbm.time import Met
from shutil import copyfile

def updateHDU(myHDU, filename, FITSfile, tstart, tstop, startStr, stopStr):
    """Single function that will update the DATE/TIME keywords in an HDU.
    """
    todaysDate = Met.now().iso()
    myHDU.header.set("DATE", todaysDate)
    myHDU.header.set("TSTART", tstart)
    myHDU.header.set("TSTOP", tstop)
    myHDU.header.set("DATE-OBS", startStr)
    myHDU.header.set("DATE-END", stopStr)
    myHDU.header.set("FILENAME", os.path.basename(filename))
    myHDU.header.set("HISTORY", "Created by writeXSPEC.py from GSPEC")
    myHDU.header.set("HISTORY", "from original file: %s" % os.path.basename(FITSfile))
    # The locations in the production PHA files were never correct! Time for them to go:
    try:
        del myHDU.header['RA_OBJ']
        del myHDU.header['DEC_OBJ']
        del myHDU.header['ERR_RAD']
    except:
        pass

def writePHAII(filename, FITSfile, tint, history, livetime, quality = None, 
               BACKGROUND = None, RATES = None, BACKFILE = None, RESPFILE = None):
    """Code that will output PHA data in the proper FITS format.
    
    INPUTS:
       filename   : STRING containing the output FITS filename
       FITSfile   : STRING containing a data FITS filename
       tint       : 2 X numTimes FLOAT ARRAY of the accumulation start and stop times
       history    : numTimes X numbins element LONG ARRAY containing the PHA counts data
       livetime   : numTimes element FLOAT ARRAY containing the livetime per accumulation
       
    KEYWORDS:
       quality    : If set, contains a numbins element array of 1s and 0s, to indicate the 
                    channels for XSPEC to ignore and notice, respectively. Should be based 
                    upon the energy selection from the lookup file. If left blank, sets 
                    the QUALITY keyword to 0 in the header.
       RATES      : If set, will write SPECTRUM extension data as rates: HDUCLAS3 = RATE.
                    In this case, the data maps as: history -> RATE and livetime -> STAT_ERR.
       BACKGROUND : If set, will change to keywords appropriate for background data. In 
                    particular, HDUCLAS2 = BKG is set. Works with RATES keyword, assuming 
                    Poisson errors in STAT_ERR.
       BACKFILE   : A STRING that names the background file to be associated with the 
                    file to be created.
       RESPFILE   : A STRING array constructed out of a response matrix filename and an 
                    index into the file that together specify the extension number of 
                    the detector redistribution matrix to be associated with each spectrum 
                    in the PHAII data file to be created. Needless to say, this should 
                    be an RSP2 file, with more than a single response matrix. If single-
                    valued, the solo entry is assumed to be applicable to every spectrum 
                    (and is written to a header keyword, rather than a column in the table).
    """
    
    numTimes = tint.shape[1]
    assert numTimes == history.shape[0], "The number of time bins must match number of spectra"

    numbins = history.shape[1]
    
    baseFile = fits.open(FITSfile)
    primary = baseFile[0].copy()
    trigtime = 0.0

    if 'TRIGTIME' in primary.header.keys():
        trigtime = primary.header['TRIGTIME']

    tstart = tint[0, 0] + trigtime
    startStr = Met(tstart).iso()
    tstop = tint[1, -1] + trigtime
    stopStr = Met(tstop).iso()
    
    updateHDU(primary, filename, FITSfile, tstart, tstop, startStr, stopStr)

    # If not set already (could happen with TTE data file)
    primary.header.set('FILETYPE', 'PHAII')

    ebounds = baseFile['EBOUNDS'].copy()
    assert numbins == ebounds.data.shape[0], "The number of channels must match number of energy bins"

    updateHDU(ebounds, filename, FITSfile, tstart, tstop, startStr, stopStr)
    
    # Our data file might be TTE, in which case, we don't want to copy the data HDU; 
    # it is constructed differently. Create the HDU from scratch:
    dataHDR = primary.header.copy()
    # Get rid of superfluous keywords (CHECKSUM and DATASUM will be added at the end):

    del dataHDR['CHECKSUM']
    del dataHDR['DATASUM']

    infiles = []
    for key in dataHDR.keys():
        if 'INFILE'in key:
            infiles.append(key)

    for key in infiles:
        del dataHDR[key]
    
    # Now build up the data header structure:
    dataHDR.set('EXTNAME', 'SPECTRUM', 'name of this binary table extension')
    bfilename = 'none'

    if BACKFILE is not None:
        bfilename = os.path.basename(BACKFILE)

    dataHDR.set('BACKFILE', bfilename, 'associated background filename')
    dataHDR.set('BACKSCAL', 1.0, 'background file scaling factor')
    dataHDR.set('CORRFILE', 'none', 'associated correction filename')
    dataHDR.set('CORRSCAL', 1.0, 'correction file scaling factor')
    rfilename = 'none'
    rspIsArray = False

    if RESPFILE is not None:
        if len(RESPFILE) == 1:
            rfilename = os.path.basename(RESPFILE[0])
            dataHDR.set('RESPFILE', rfilename, 'associated redistribution matrix filename')
        else:
            rspCol = fits.Column(name='RESPFILE', format='40A', array = RESPFILE)
            rspIsArray = True
    else:
        dataHDR.set('RESPFILE', rfilename, 'associated redistribution matrix filename')

    dataHDR.set('ANCRFILE', 'none', 'associated ancilary response filename')
    dataHDR.set('AREASCAL', 1.0, 'Effective area scaling factor')
    dataHDR.set('SYS_ERR', 0, 'no systematic error specified')
    #dataHDR.set('QUALITY', 0, 'no systematic error specified')
    dataHDR.set('GROUPING', 0, 'Grouping flag')
    
    dataHDR.set('HDUCLASS', 'OGIP', 'format conforms to OGIP standard')
    dataHDR.set('HDUCLAS1', 'SPECTRUM', 'PHA dataset (OGIP-92-007)')
    hduClas2 = 'TOTAL'
    clas2comm  = 'Indicates gross data (source + background)'
    myPoissErr = True

    if BACKGROUND is not None:
        hduClas2 = 'BKG'
        clas2comm  = 'Background PHA Spectrum'
        myPoissErr = False

    dataHDR.set('HDUCLAS2', hduClas2, clas2comm)
    hduClas3  = 'COUNT'
    clas3comm = 'PHA data stored as counts'

    if RATES is not None:
        hduClas3 = 'RATE'
        clas3comm = 'PHA data stored as rates'

    dataHDR.set('HDUCLAS3', hduClas3, clas3comm)
    dataHDR.set('HDUCLAS4', 'TYPEII', 'Multiple PHA dataset')
    dataHDR.set('HDUVERS', '1.2.1', 'Format version number')
    dataHDR.set('POISSERR', myPoissErr, 'Whether Poisson errors are applicable')
    dataHDR.set('CHANTYPE', 'PHA', 'channel type (PHA, PI etc.)')
    dataHDR.set('DETCHANS', numbins, 'total number possible channels')

    if quality is None:
        quality = np.zeros(numTimes)
        qualCol = fits.Column(name='QUALITY', format='1I', array = quality)
        qualCom = 'Quality (no non-zero bits are set)'
    else:
        qual_arr = np.tile(quality, (numTimes,1))
        qualCol = fits.Column(name='QUALITY', format='%iI' % numbins, array = qual_arr)
        qualCom = 'Quality (non-zero bits indicate bad channels)'

    timeCol = fits.Column(name='TIME', unit = 's', format='1D', bzero = trigtime, array = tint[0, :])
    timeCom = 'Start time of the accumulation'
    endtCol = fits.Column(name='ENDTIME', unit = 's', format='1D', bzero = trigtime, array = tint[1, :])
    endtCom = 'End time of the accumulation'

    if RATES is not None:
        rateCol = fits.Column(name='RATE', unit = 'count/s', format='%iD' % numbins, array = history)
        col1Com = 'Rate per Pulse Height Analyzer (PHA) Channel'
        statCol = fits.Column(name='STAT_ERR', unit = 'count/s', format='%iD' % numbins, array = livetime)
        col2Com = 'Statistical errors in RATE (assumes Poisson)'
        dataColArr = [rateCol, statCol, qualCol, timeCol, endtCol]

        if BACKGROUND is not None:
            expCol = fits.Column(name='EXPOSURE', unit = 's', format='1E', array = tint[1, :] - tint[0, :])
            expCom = 'Exposure time = Accumulation time'
            dataColArr.append(expCol)
    else:
        countsCol = fits.Column(name='COUNTS', unit = 'counts', format='%iJ' % numbins, array = history)
        col1Com = 'Counts per Pulse Height Analyzer (PHA) Channel'
        expCol = fits.Column(name='EXPOSURE', unit = 's', format='1E', array = livetime)
        col2Com = 'Accumulation time - deadtime'
        dataColArr = [countsCol, expCol, qualCol, timeCol, endtCol]
        #dataCols = fits.ColDefs([countsCol, expCol, qualCol, timeCol, endtCol])
    
    if rspIsArray is True:
        dataColArr.append(rspCol)

    dataCols = fits.ColDefs(dataColArr)
    dataHDU = fits.BinTableHDU.from_columns(dataCols, header = dataHDR)

    # Now that we have an HDU, we can update the comments for each of the columns:
    dataHDU.header.comments['TTYPE1'] = col1Com
    dataHDU.header.comments['TTYPE2'] = col2Com
    dataHDU.header.comments['TTYPE3'] = qualCom
    dataHDU.header.comments['TTYPE4'] = timeCom
    dataHDU.header.comments['TTYPE5'] = endtCom

    if BACKGROUND is not None:
        dataHDU.header.comments['TTYPE6'] = expCom
    
    # Construct the GTI extension from the original, but with the new bounds:
    gtiHDR = baseFile['GTI'].header.copy()

    startCol = fits.Column(name='START', unit = 's', format='1D', bzero = trigtime, array = (tint[0, 0],))
    stopCol = fits.Column(name='STOP', unit = 's', format='1D', bzero = trigtime, array = (tint[1, -1],))
    gtiCols = fits.ColDefs([startCol, stopCol])
    gtiHDU = fits.BinTableHDU.from_columns(gtiCols, header = gtiHDR)

    updateHDU(gtiHDU, filename, FITSfile, tstart, tstop, startStr, stopStr)

    # Now that we have an HDU, we can update the comments for each of the columns:
    gtiHDU.header.comments['TTYPE1'] = 'Beginning time of accumulation'
    gtiHDU.header.comments['TZERO1'] = '[s] Offset, equal to TRIGTIME'
    gtiHDU.header.comments['TTYPE2'] = 'Ending time of accumulation'
    gtiHDU.header.comments['TZERO2'] = '[s] Offset, equal to TRIGTIME'
    
    # Construct and write the file:
    myHDUlist = fits.HDUList([primary, ebounds, dataHDU, gtiHDU])
    myHDUlist.writeto(filename, checksum = True)
    

def writePHAI(filename, FITSfile, tint, history, livetime, quality = None, 
               BACKGROUND = None, RATES = None, BACKFILE = None, RESPFILE = None):
    """Code that will output PHA data in the proper FITS format.
    
    INPUTS:
       filename   : STRING containing the output FITS filename
       FITSfile   : STRING containing a data FITS filename
       tint       : 2 X 1 FLOAT ARRAY of the accumulation start and stop times
       history    : numbins element LONG ARRAY containing the PHA counts data
       livetime   : FLOAT containing the livetime for the accumulation
       
    KEYWORDS:
       quality    : If set, contains a numbins element array of 1s and 0s, to indicate the 
                    channels for XSPEC to ignore and notice, respectively. Should be based 
                    upon the energy selection from the lookup file. If left blank, sets 
                    the QUALITY keyword to 0 in the header.
       RATES      : If set, will write SPECTRUM extension data as rates: HDUCLAS3 = RATE.
                    In this case, the data maps as: history -> RATE and livetime -> STAT_ERR.
       BACKGROUND : If set, will change to keywords appropriate for background data. In 
                    particular, HDUCLAS2 = BKG is set. Works with RATES keyword, assuming 
                    Poisson errors in STAT_ERR.
       BACKFILE   : A STRING that names any background file to be associated with the 
                    file to be created.
       RESPFILE   : A STRING that names any detector redistribution matrix file to be  
                    associated with the file to be created.
    """
    
    numbins = history.shape[0]
    
    baseFile = fits.open(FITSfile)
    primary = baseFile[0].copy()
    trigtime = 0.0
    if 'TRIGTIME' in primary.header.keys():
        trigtime = primary.header['TRIGTIME']
    tstart = tint[0] + trigtime
    startStr = Met(tstart).iso()
    tstop = tint[-1] + trigtime
    stopStr = Met(tstop).iso()
    todaysDate = Met.now().iso()
    
    updateHDU(primary, filename, FITSfile, tstart, tstop, startStr, stopStr)
    if BACKGROUND is not None:
        filetype = 'GBM BACK'
    else:
        filetype = 'SPECTRUM'
    primary.header.set('FILETYPE', filetype)

    ebounds = baseFile['EBOUNDS'].copy()
    assert numbins == ebounds.data.shape[0], "The number of channels must match number of energy bins"
    updateHDU(ebounds, filename, FITSfile, tstart, tstop, startStr, stopStr)
    
    # Our data file might be TTE, in which case, we don't want to copy the data HDU; 
    # it is constructed differently. Create the HDU from scratch:
    dataHDR = primary.header.copy()
    # Get rid of superfluous keywords (CHECKSUM and DATASUM will be added at the end):
    del dataHDR['CHECKSUM']
    del dataHDR['DATASUM']
    infiles = []
    for key in dataHDR.keys():
        if 'INFILE'in key:
            infiles.append(key)
    for key in infiles:
        del dataHDR[key]
    
    # Now build up the data header structure:
    dataHDR.set('EXTNAME', 'SPECTRUM', 'name of this binary table extension')
    bfilename = 'none'
    if BACKFILE is not None:
        bfilename = os.path.basename(BACKFILE)
    dataHDR.set('BACKFILE', bfilename, 'associated background filename')
    dataHDR.set('BACKSCAL', 1.0, 'background file scaling factor')
    dataHDR.set('CORRFILE', 'none', 'associated correction filename')
    dataHDR.set('CORRSCAL', 1.0, 'correction file scaling factor')
    rfilename = 'none'
    if RESPFILE is not None:
        rfilename = os.path.basename(RESPFILE)
    dataHDR.set('RESPFILE', rfilename, 'associated redistribution matrix filename')
    dataHDR.set('ANCRFILE', 'none', 'associated ancilary response filename')
    dataHDR.set('AREASCAL', 1.0, 'Effective area scaling factor')
    dataHDR.set('SYS_ERR', 0, 'no systematic error specified')
    if quality is None:
        dataHDR.set('QUALITY', 0, 'no bad channels to be excluded')
    dataHDR.set('GROUPING', 0, 'no energy channel grouping in data')
    
    dataHDR.set('HDUCLASS', 'OGIP', 'format conforms to OGIP standard')
    dataHDR.set('HDUCLAS1', 'SPECTRUM', 'PHA dataset (OGIP-92-007)')
    hduClas2 = 'TOTAL'
    clas2comm  = 'Indicates gross data (source + background)'
    myPoissErr = True
    if BACKGROUND is not None:
        hduClas2 = 'BKG'
        clas2comm  = 'Background PHA Spectrum'
        myPoissErr = False
    dataHDR.set('HDUCLAS2', hduClas2, clas2comm)
    hduClas3  = 'COUNT'
    clas3comm = 'PHA data stored as counts'
    if RATES is not None:
        hduClas3 = 'RATE'
        clas3comm = 'PHA data stored as rates'
    dataHDR.set('HDUCLAS3', hduClas3, clas3comm)
    dataHDR.set('HDUCLAS4', 'TYPEI', 'Single PHA dataset')
    dataHDR.set('HDUVERS', '1.2.1', 'Format version number')
    dataHDR.set('POISSERR', myPoissErr, 'Whether Poisson errors are applicable')
    dataHDR.set('CHANTYPE', 'PHA', 'channel type (PHA, PI etc.)')
    dataHDR.set('DETCHANS', numbins, 'total number possible channels')
    if BACKGROUND is not None:
        myExposure = tint[-1] - tint[0]
        dataHDR.set('EXPOSURE', myExposure, 'Accumulation time')
    else:
        dataHDR.set('EXPOSURE', livetime, 'Accumulation time - deadtime')

    channels = np.arange(numbins) + 1
    chanCol = fits.Column(name='CHANNEL', format='1I', array = channels)
    chanCom = 'Pulse Height Analyzer (PHA) Channel'
    colList = [chanCol]
    if RATES is not None:
        rateCol = fits.Column(name='RATE', unit = 'count/s', format='1D', array = history)
        rateCom = 'Rate per Pulse Height Analyzer (PHA) Channel'
        colList.append(rateCol)
        statCol = fits.Column(name='STAT_ERR', unit = 'count/s', format='1D', array = livetime)
        statCom = 'Statistical errors in RATE (assumes Poisson)'
        colList.append(statCol)
    else:
        countsCol = fits.Column(name='COUNTS', unit = 'counts', format='1J', array = history)
        countsCom = 'Counts per Pulse Height Analyzer (PHA) Channel'
        colList.append(countsCol)
    
    if quality is not None:
        qualCol = fits.Column(name='QUALITY', format='1I', array = quality)
        qualCom = 'indicates bad channels to be excluded'
        colList.append(qualCol)
    dataCols = fits.ColDefs(colList)
    dataHDU = fits.BinTableHDU.from_columns(dataCols, header = dataHDR)
    # Now that we have an HDU, we can update the comments for each of the columns:
    dataHDU.header.comments['TTYPE1'] = chanCom
    if RATES is not None:
        dataHDU.header.comments['TTYPE2'] = rateCom
        dataHDU.header.comments['TTYPE3'] = statCom
    else:
        dataHDU.header.comments['TTYPE2'] = countsCom
    
    # Construct the GTI extension from the original, but with the new bounds:
    gtiHDR = baseFile['GTI'].header.copy()

    startCol = fits.Column(name='START', unit = 's', format='1D', bzero = trigtime, array = (tint[0],))
    stopCol = fits.Column(name='STOP', unit = 's', format='1D', bzero = trigtime, array = (tint[-1],))
    gtiCols = fits.ColDefs([startCol, stopCol])
    gtiHDU = fits.BinTableHDU.from_columns(gtiCols, header = gtiHDR)
    updateHDU(gtiHDU, filename, FITSfile, tstart, tstop, startStr, stopStr)
    # Now that we have an HDU, we can update the comments for each of the columns:
    gtiHDU.header.comments['TTYPE1'] = 'Beginning time of accumulation'
    gtiHDU.header.comments['TZERO1'] = '[s] Offset, equal to TRIGTIME'
    gtiHDU.header.comments['TTYPE2'] = 'Ending time of accumulation'
    gtiHDU.header.comments['TZERO2'] = '[s] Offset, equal to TRIGTIME'
    
    # Construct and write the file:
    myHDUlist = fits.HDUList([primary, ebounds, dataHDU, gtiHDU])
    myHDUlist.writeto(filename, checksum = True)


def writeXSPECFiles(source_select, bkgd, RSPobj, outPHAfile, outBAKfile, outRSPfile, MULTI=False):
    """Function that gathers all the info required to write out XSPEC PHA and BAK files
    
    Parameters:
    -----------
    source_select: an instance of Selection; made from Ctime, Cspec or TTE data objects
        Encapsulates the user's time and energy selections and binning.
        The time selection is used to extract the counts data. What this means will depend 
        on the MULTI keyword. If it is False, the time selection will be integrated.
    bkgd: an instance of Background.
        The caller must perform the fit before handing off the background to this routine.
    RSPobj: a CTIME or CSPEC multiple-matrix reponse file object RSP
    outPHAfile: string
        (Fully-qualified) name of the PHA file to be written
    outBAKfile: string
        (Fully-qualified) name of the BAK file to be written
    outRSPfile: string
        (Fully-qualified) name of the RSP file to be written
    MULTI: Boolean
        If True, will attempt to write PHAII source and background files (one each) 
        for batch operations with XSPEC
    """
    
    FITSfileName = source_select._data_obj.filename
    
    energy_sel = source_select._energy_selections # energy_bound
    spectrum = source_select.full_spectrum_integrated()
    src_spectrum = np.array(spectrum.values)
    livetime = source_select.livetime

    numBins = len(src_spectrum)
    luQuality = np.zeros(numBins) + 1.

    for energy_intvl in energy_sel:
        mask = (spectrum.edges[1:] > energy_intvl[0]) & (spectrum.edges[0:-1] < energy_intvl[1])
        luQuality[mask] = 0.
    
    # Now branch to write either single or multiple PHA files for XSPEC
    if MULTI == False:

        #print("Creating PHAI files...")

        # the background spectrum integrated over the source selection
        bkgd_spectrum = bkgd.source_spectrum_integrated(source_select)
        bak_spectrum = np.array(bkgd_spectrum.values)/livetime
        bak_errors = np.array(bkgd_spectrum.uncertainty)/livetime
    
        # get the lightcurve data and rebin
        rate_hist = RateHisto.merge(source_select.lightcurve)
        full_rates = np.array(rate_hist.values)
        norm = sum(full_rates)
        norm_hist = full_rates / norm
        time_bins = np.array(rate_hist.bounds)
        time_bounds = source_select.time_bound

        # write out the resulting files for XSPEC:
        
        writePHAI(outPHAfile, FITSfileName, time_bounds, src_spectrum, livetime, 
                  quality = luQuality, BACKFILE = outBAKfile, RESPFILE = outRSPfile)
        writePHAI(outBAKfile, FITSfileName, time_bounds, bak_spectrum, bak_errors, 
                  BACKGROUND = True, RATES = True)
        RSPobj.writeResponse(time_bins, norm_hist, rspFileName = outRSPfile)

    else:

        #print("Creating PHAII files...")

        # source selection for each time bin
        merged_lightcurve = RateHisto.merge(source_select.lightcurve)
        tb = np.array(merged_lightcurve.bounds)
        time_bins = np.array([merged_lightcurve.bounds[0],merged_lightcurve.bounds[1]]).T
        livetimes = merged_lightcurve.exposure
        all_spectra = source_select.full_spectrum_for_each_bin()
        src_spectra = np.array([spectrum.values for spectrum in all_spectra])
        
        # the background spectrum integrated over the time bin
        bkgd_spectra = bkgd.source_spectrum_for_each_bin(source_select)
        bak_spectra = np.array([spectrum.rate for spectrum in bkgd_spectra])
        bak_uncerts = np.array([spectrum.uncertainty for spectrum in bkgd_spectra])

        assert time_bins.shape[0] > 1, 'PHAII file requested but found only a single spectrum in the selection!\n'

        # Extract the path to the output pha file        
        outPHAfile_basename = os.path.basename(outPHAfile)
        outDirectory = outPHAfile.strip(outPHAfile_basename)

        rspIdxArr = []

        for ti in time_bins:
            # Build up the indexed list of rsp2 matrices corresponding with the time bins:
            rspIdxName = RSPobj.checkResponse(ti)
            rspIdxArr.append(rspIdxName)
            
        # write out the resulting files for XSPEC:
        writePHAII(outPHAfile, FITSfileName, tb, src_spectra, livetimes, \
                  quality = luQuality, BACKFILE = outBAKfile, RESPFILE = rspIdxArr)
        writePHAII(outBAKfile, FITSfileName, tb, bak_spectra, bak_uncerts, \
                  BACKGROUND = True, RATES = True)

        # Generate the destination rsp2 filename
        inRSPFile = RSPobj.filename
        inRSPFile_basename = os.path.basename(inRSPFile)
        outRSPfile = outDirectory + inRSPFile_basename

        # print("Original RSP2 file:")
        # print(inRSPFile)

        # print("Destination RSP2 file:")
        # print(outRSPfile)

        # Copy the rsp2 file
        copyfile(inRSPFile, outRSPfile)
    


if __name__ == '__main__':
    FITSfile = 'test/glg_cspec_n6_bn090926181_v00.pha'
    PHAobj = Cspec(FITSfile)  # TTE(FITSfile)
    LUfile = 'test/glg_cspec_n6_bn090926181_v04.lu'
    LUobj = RmfitLookup(LUfile)  #, ti_file='test/glg_tte_n9_bn090131090_v01.ti')
    # Background time selections
    time_ranges = LUobj.background_selection()
    time_ranges = [t.tolist() for t in time_ranges]
    # one energy range
    energy_range = [PHAobj.ebounds[0][1],PHAobj.ebounds[-1][2]] 
    energy_sel = [selection.tolist() for selection in LUobj.energy_selection()] #LUobj.energy_selection() #[0].tolist()
    # background selection
    bkgd_select = Selection(PHAobj, time_selections=time_ranges, 
                            energy_selections=energy_sel)

    # background fit
    bkgd = Background(bkgd_select)
    bkgd.fit(polynomial, LUobj.background_poly_order)
        
    # the source selection
    source_range = [selection.tolist() for selection in LUobj.time_selection()] #[0].tolist()
    # energy_sel = LUobj.energy_selection()[0].tolist()
    source_select = Selection(PHAobj, time_selections=source_range, 
                              energy_selections=energy_sel)
    # check rebinning and adjust edges if necessary:
    if LUobj.is_time_rebinned == True:
        if isinstance(PHAobj, TTE):
            time_edges = LUobj.tte_edges()
            source_select.set_time_binning(bin_by_edges, time_edges)
        else:
            time_bin_indices = LUobj.rebinned_time_edges()
            source_select.set_time_binning(rebin_by_edge_index, time_bin_indices)
        
    RSPfile = 'test/glg_cspec_n6_bn090926181_v00.rsp2'
    RSP2 = RSP(RSPfile)
    # These output files should be compared with the v00 ones as follows:
    # cd test
    # fdiff phaI_n6_bn090926181_v00.pha phaI_n6_bn090926181_v01.pha exclude="FILENAME,BACKFILE,RESPFILE"
    outPHAfile = 'test/phaI_n6_bn090926181_v01.pha'
    outBAKfile = 'test/phaI_n6_bn090926181_v01.bak'
    outRSPfile = 'test/phaI_n6_bn090926181_v01.rsp'
    outPHAIIfile = 'test/phaII_n6_bn090926181_v01.pha'
    outBAKIIfile = 'test/phaII_n6_bn090926181_v01.bak'
    outRSPIIfile = 'test/glg_cspec_n6_bn090926181_v00.rsp2'  #(Ignored)
    writeXSPECFiles(source_select, bkgd, RSP2, outPHAfile, outBAKfile, outRSPfile)
    writeXSPECFiles(source_select, bkgd, RSP2, outPHAIIfile, outBAKIIfile, outRSPIIfile, MULTI = True)
    