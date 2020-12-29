import astropy.io.fits as fits
import numpy as np
from os.path import basename
from collections import OrderedDict

from gbm.time import Met
from .data import Data

class RSP(Data):
    """Container class for a RSP or RSP2 file
    Attributes:
    -----------
    drm_list: list
        The list of DRMs contained in the file
    num_drms: int
        The number of DRMs contained in the file
    trigtime: float
        The trigger time
    times: np.array
        The times for which the DRM(s) were generated

    Public Methods:
    ---------------
    checkResponse: 
        Check to see if a time interval is covered by a single matrix
    weightResponse: 
        Constructs a single weighted response matrix
    writeResponse: 
        Write a response file
    """
    def __init__(self, filename):
        self.header = OrderedDict()
        self.ebounds = None
        self.pbounds = None
        self.drm_list = []
        
        self.num_drms = 0
        self.trigtime = None
        self.times = None
        
        super(RSP, self).__init__(filename)
        if self.filename is not None:
            self._from_file()
            self._set_properties()
        
    def _from_file(self):

        with fits.open(self.filename) as hdulist:
            ispec = 1
            for hdu in hdulist:
                if hdu.name == 'SPECRESP MATRIX':
                    name = hdu.name+str(ispec)
                    ispec += 1
                else:
                    name = hdu.name
                self.header.update({name: hdu.header})
            self.ebounds = hdulist['EBOUNDS'].data

            self.num_drms = self.header['PRIMARY']['DRM_NUM']
            for idrm in range(1,self.num_drms+1):
                drm_data = hdulist['SPECRESP MATRIX', idrm].data
                self.pbounds = (drm_data['ENERG_LO'], drm_data['ENERG_HI'])
                self.drm_list.append(self._decompress_drm(drm_data, idrm))
    
    def _set_properties(self):
        self.trigtime = self.header['PRIMARY']['TRIGTIME']
        
        rsp_times = []
        for ext, header in self.header.items():
            if 'SPECRESP MATRIX' not in ext:
                continue
            tstart = header['TSTART']
            tstop = header['TSTOP']
            rsp_times.append(tstart-self.trigtime)
        rsp_times.append(tstop-self.trigtime)
        self.times = rsp_times
    
    def _decompress_drm(self, drm_data, spec_index):
        """The GBM DRM is not stored as 2D matrix and may be compressed, so this
        constructs and decompresses the DRM.

        Parameters:
        -----------
        drm_data: np.rec_array
            The DRM data
        spec_index: int
            The index of the DRM in the event of multiple DRMS

        Returns:
        -----------
        drm: np.array
            2D matrix containing the unpacked DRM
        """

        num_photon_bins = self.header['SPECRESP MATRIX'+str(spec_index)]['NAXIS2']
        first_chans = drm_data['F_CHAN']

        matrix = drm_data['MATRIX']
        energy_channel_max = self.ebounds['E_MAX']
        num_channels = energy_channel_max.size

        drm = np.zeros((num_photon_bins, num_channels))
        for ibin, effective_area in enumerate(matrix):
            drm[ibin, first_chans[ibin][0]-1:] = effective_area

        return drm

    def checkResponse(self, tint):
        """Given a time interval, checks to see if it is covered by a single matrix in an 
        RSP2 file. Constructs a string corresponding to the XSPEC response command param.
        
        Parameters:
        -----------
        tint: If not empty, contains the start and end times for requested response. 

        Returns:
        -----------
        responseStr: If tint is covered by a single matrix, the XSPEC response command 
        string is constructed: self.filename{i}, where 'i' is the matrix index.
        """
        
        responseStr = ''
        myStart = tint[0]
        myStop = tint[1]
        tiList = np.array(self.times[:-1])
        myFilename = basename(self.filename)
        minInd = np.argmin(np.absolute(tiList - myStart))
        maxInd = np.argmin(np.absolute(tiList - myStop))
        if minInd == maxInd:
            #Python arrays being zero-indexed and the RSP file responses not:
            responseStr = myFilename +'{' + str(minInd+1) + '}'
        else:
            # Find the matrix whose time span contains the center of our time bin:
            aveTi = (myStart + myStop) / 2.0
            myInd = np.argmin(np.absolute(tiList - aveTi))
            if (tiList[myInd] - aveTi) > 0.0:
                responseStr = myFilename +'{' + str(myInd) + '}'
            else:
                responseStr = myFilename +'{' + str(myInd+1) + '}'
            
        return responseStr

    def weightResponse(self, tint = None, norm_hist = None, time_bins = None):
        """Constructs a single weighted response matrix, given a time interval or set of 
        normalized weights over several time bins.
        
        Parameters:
        -----------
        tint: If not empty, contains the start and end times for requested response. These 
            will be averaged and the response closest to the averaged time returned.
        norm_hist and time_bins: For longer spectral accumulations, the time series of 
            responses needs to be weighted by the counts history. norm__hist provides the 
            weights, and must sum to 1.0. time_bins holds the start and end times of the 
            corresponding spectral acumulations and must be equal in length to norm_hist. 
            It has the shape 2 X n.

        Returns:
        -----------
        response: The weighted response matrix.
        """

        if self.num_drms > 1:
            # Note that there is an extra time in the array and we don't want to
            # generate an index error due to an edge case:
            tiList = np.array(self.times[:-1])

            if norm_hist is not None:
                lenNorm = len(norm_hist)
                myTimes = np.array(time_bins).sum(0) / 2.0
                if myTimes.shape[0] != lenNorm:
                    raise ValueError("norm_hist and time_bins are not equal length")
                #== We want to weight the matrices by counts:
                weightV = []
                weightS = 0.0
                weightD = []
                currDRM = -1
                for hi in range(lenNorm):
                    myTime = myTimes[hi]
                    minInd = np.argmin(np.absolute(tiList - myTime))
                    if currDRM == -1:
                        currDRM = minInd
                        weightD.append(currDRM)

                    if currDRM == minInd:
                        weightS += norm_hist[hi]
                    else:
                        currDRM = minInd
                        weightV.append(weightS)
                        weightD.append(currDRM)
                        weightS = norm_hist[hi]

                weightV.append(weightS)
                num_photon_bins = len(self.pbounds[0])
                num_channels = len(self.ebounds)
                response = np.zeros((num_photon_bins, num_channels))
                for wi, myWeight in zip(weightD, weightV):
                    response = response + self.drm_list[wi] * myWeight

            else:
                # Start simply - average the time interval:
                targetT = sum(tint) / 2.0
                # Find the matrix closest to the trigger:
                minInd = np.argmin(np.absolute(tiList - targetT))
                # TODO: Inspection gave unexpected type error. minInd may be an integer or an array slice (check this).
                response = self.drm_list[minInd]

        else:
            # Only one DRM in the file; return it:
            response = self.drm_list[0]

        return response

    def writeResponse(self, time_bins, norm_hist, rspFileName = 'XSPEC_001.rsp'):
        """Given a set of time intervals, construct a weighted single matrix RSP file.
        
        Parameters:
        -----------
        time_bins: Provides the start and end times of the counts history time bins 
            which must be equal in length to norm_hist. It has the shape 2 X n.
        norm_hist: A vector of n source counts history weights, which must sum to 1.0. 
            As the weights are derived from the source counts, any time bins that cover 
            background-only periods should be zero. The weights should never be negative!
        rspFileName: The name of the output response matrix file. Detector, datatype and 
            temporal info should be included to make it unique. 

        Returns:
        -----------
        rspFileName upon success.
        
        Side Effect:
        -----------
        Writes out a single weighted matrix response file using the base name in 
        rspFileName. The intent is that this will be a temporary file, living long 
        enough to support the analysis of a spectrum based upon the specific time 
        interval accumulation requested. 
        """

        myFilename = self.filename
        rsp2File = fits.open(myFilename)
        primary = rsp2File[0].copy()
        primary.header.set("DRM_NUM", 1)
        trigtime = 0.0
        if self.trigtime is not None:
            trigtime = self.trigtime
        tstart = time_bins[0, 0] + trigtime
        startStr = Met(tstart).iso()
        tstop = time_bins[1, -1] + trigtime
        stopStr = Met(tstop).iso()
        todaysDate = Met.now().iso()
        primary.header.set("DATE", todaysDate)
        primary.header.set("TSTART", tstart)
        primary.header.set("TSTOP", tstop)
        primary.header.set("DATE-OBS", startStr)
        primary.header.set("DATE-END", stopStr)
        primary.header.set("FILENAME", rspFileName)
        primary.header.set("HISTORY", "Created by writeResponse() in data.py from GSPEC")
        primary.header.set("HISTORY", "from original file: %s" % myFilename)
        primary.header.set("COMMENT", "The matrix has been weighted by the count history")

        ebounds = rsp2File['EBOUNDS'].copy()
        ebounds.header.set("TSTART", tstart)
        ebounds.header.set("TSTOP", tstop)
        ebounds.header.set("DATE", todaysDate)
        ebounds.header.set("DATE-OBS", startStr)
        ebounds.header.set("DATE-END", stopStr)

        matrix = self.weightResponse(norm_hist = norm_hist, time_bins = time_bins)

        matrixHDU = rsp2File[2]
        newcols = []
        for col in matrixHDU.columns:
            if col.format.find("P") == 0:
                #Variable-length
                colF = col.format.split("(")
                coltype = colF[0].replace("P","")
                length = colF[1].replace(")","")
                newFormat = '%s%s' %(length,coltype)
                if col.name != 'MATRIX':
                    newMatrix = matrixHDU.data.field(col.name)
                else:
                    newMatrix = matrix
                newcols.append(fits.Column(col.name, newFormat, col.unit, col.null,
                        col.bscale, col.bzero, col.disp, col.start, col.dim, newMatrix))
            else:
                newcols.append(fits.Column(col.name, col.format, col.unit, col.null,
                        col.bscale, col.bzero, col.disp, col.start, col.dim,
                        matrixHDU.data.field(col.name)))

        response = fits.BinTableHDU.from_columns(newcols,header=matrixHDU.header)
        response.header.set("TSTART", tstart)
        response.header.set("TSTOP", tstop)
        response.header.set("DATE", todaysDate)
        response.header.set("DATE-OBS", startStr)
        response.header.set("DATE-END", stopStr)
        HDUlist = fits.HDUList([primary, ebounds, response])
        HDUlist.writeto(rspFileName, output_verify='fix', checksum=True)

        return rspFileName
