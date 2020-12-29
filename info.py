import os
import numpy as np
from astropy.io import fits
from XspecT.Preparation.downinfo import DownloadGBMcat

class EventInfo(object):
    _GBM_dtrs = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']
    _GBM_download_address = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts"
    _LAT_downlad_address = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/triggers"
    
    def __init__(self, grb, tRan = [None, None], eRan = [100, 10000], GBM_catalog=None, LAT_catalog=None, **kwargs):

        self.__dict__.update(**kwargs)

        grb = str(grb)
        self.event_name = grb[:6]
        self.yr = grb[:2]
        self.full_name = grb

        if GBM_catalog==None:
            SCRIPT_DIR=os.environ.get('XspecT')
            GBM_catalog = SCRIPT_DIR+'/RefData/GBMcatalog.fits'
        
        if LAT_catalog==None:
            SCRIPT_DIR=os.environ.get('XspecT')
            LAT_catalog = SCRIPT_DIR+'/RefData/LAT2CATALOG-v16-LTF.fits'

        self.event_info = fits.open(GBM_catalog)
        self.event_info_LAT = fits.open(LAT_catalog)

        self.usedGBM = []
        self._info_num = self.event_info[1].data["TRIGGER_NAME"]=='BN'+self.full_name
        self._LAT_info_num = self.event_info_LAT[1].data["GRBNAME"]==self.full_name
        
        self._address_LAT = self.full_name+'/LAT/'+self.event_name
        self._address_Xspec = self.full_name+'/Xspec/'+self.event_name

        
        if sum(self._info_num) == 1:
            detNm = self.event_info[1].data["SCAT_DETECTOR_MASK"][self._info_num]
            
            for det, i in zip(detNm[0], range(len(detNm[0]))):
                if int(det):
                    self.usedGBM.append(self._GBM_dtrs[i])

            try:
                self._bcat=fits.open('./{}/GBM/{}-bcat.fit'.format(self.full_name,self.event_name))
                self._scat=fits.open('./{}/GBM/{}-scat.fit'.format(self.full_name,self.event_name))
            except:
                DownloadGBMcat(self.full_name)
                self._bcat=fits.open('./{}/GBM/{}-bcat.fit'.format(self.full_name,self.event_name))
                self._scat=fits.open('./{}/GBM/{}-scat.fit'.format(self.full_name,self.event_name))
            

            if len(self.usedGBM) == 0:
                for dtr in self._scat[1].data["DETNAM"]:
                    dtr=dtr.replace("NAI_0", "n")
                    dtr=dtr.replace("NAI_10", "na")
                    dtr=dtr.replace("NAI_11", "nb")
                    dtr=dtr.replace("BGO_0", "b")
                    self.usedGBM.append(dtr)

            # Check consistancy
            self.trigger = self._bcat[0].header["TRIGTIME"]
            if self._bcat[0].header["T90START"] - self.event_info[1].data["T90_START"][self._info_num ][0] > 1e-2:
                print("T90_start Error")
                print(self._bcat[0].header["T90START"] - self.event_info[1].data["T90_START"][self._info_num])
            if self._bcat[0].header["T90"] - self.event_info[1].data["T90"][self._info_num][0] > 1e-2:
                print("T90 Error")
                print(self._bcat[0].header["T90"], self.event_info[1].data["T90"][self._info_num])
            

        elif sum(self._info_num) == 0:
            DownloadGBMcat(self.full_name)
            self._bcat=fits.open('./{}/GBM/{}-bcat.fit'.format(self.full_name,self.event_name))
            self._scat=fits.open('./{}/GBM/{}-scat.fit'.format(self.full_name,self.event_name))

            for dtr in self._scat[1].data["DETNAM"]:
                dtr=dtr.replace("NAI_0", "n")
                dtr=dtr.replace("NAI_10", "na")
                dtr=dtr.replace("NAI_11", "nb")
                dtr=dtr.replace("BGO_0", "b")
                self.usedGBM.append(dtr)

            self._bcat=fits.open('./{}/GBM/{}-bcat.fit'.format(self.full_name,self.event_name))
            self._scat=fits.open('./{}/GBM/{}-scat.fit'.format(self.full_name,self.event_name))


        self.t90 = float(self._bcat[0].header["T90"])
        self.t05 = float(self._bcat[0].header["T90START"])
        self.t95 = self._bcat[0].header["T90START"]+self._bcat[0].header["T90"]

        if sum(self._LAT_info_num) == 0:
            print("This GRB is not in the LAT catalog. Please update the catalog file.")
            print("You need to define following parameters.")
            print("refRa, refDec")

        elif sum(self._LAT_info_num) == 1:
            self.t05 = self.event_info_LAT[1].data["GBMt05"][self._LAT_info_num][0]
            self.t95 = self.event_info_LAT[1].data["GBMT95"][self._LAT_info_num][0]
            
            self.refRa = self.event_info_LAT[1].data["RA"][self._LAT_info_num][0]
            self.refDec = self.event_info_LAT[1].data["DEC"][self._LAT_info_num][0]

        if tRan[0] !=None:
            if tRan[0] == "T0": _time_start = 0.
            elif tRan[0] == "T05": _time_start = self.t05
            elif tRan[0] == "T95": _time_start = self.t95
            elif tRan[0] == "LATT0" and sum(self._LAT_info_num) == 1:
                _time_start = self.event_info_LAT[1].data["TL0"][self._LAT_info_num][0]         
            elif tRan[0] == "LATT1" and sum(self._LAT_info_num) == 1:
                _time_start = self.event_info_LAT[1].data["TL100"][self._LAT_info_num][0]
            elif isinstance(tRan[0], int) or isinstance(tRan[0], float):
                _time_start = tRan[0]
            else:
                raise("Check your start time (tRan[0]): e.g., T05, T95, LATT0, LATT1, float, int, ...")
        else:
            _time_start = self.t05


        if tRan[1] !=None:
            if tRan[1] == "T0": _time_end = 0.
            elif tRan[1] == "T05": _time_end = self.t05
            elif tRan[1] == "T95": _time_end = self.event_info_LAT[1].data["GBMT95"][self._LAT_info_num][0]         
            elif tRan[1] == "LATT0" and sum(self._LAT_info_num) == 1:
                _time_end = self.event_info_LAT[1].data["TL0"][self._LAT_info_num][0]         
            elif tRan[1] == "LATT1" and sum(self._LAT_info_num) == 1:
                _time_end = self.event_info_LAT[1].data["TL100"][self._LAT_info_num][0]
            elif isinstance(tRan[1], int) or isinstance(tRan[1], float):
                _time_end = tRan[1]
            else:
                raise("Check your end time (tRan[1]): e.g., T05, T95, LATT0, LATT1, float, int, ...")
        else:
            _time_end = self.t95

        if _time_end < _time_start:
            _time_start, _time_end = _time_end, _time_start 

        self.toi = [_time_start, _time_end, self.trigger]
        if not(os.path.isfile(self._address_Xspec+"-ti.npy")):
            self.time_intervals = np.asarray([self.toi])
            self._ti_num = 0
            np.save("./"+self._address_Xspec+"-ti", self.time_intervals)
        else:
            self.time_intervals = np.load("./"+self._address_Xspec+"-ti.npy")
            
            self._ti_num = None
            
            for i in range(len(self.time_intervals)):
                if sum(self.toi == self.time_intervals[i]) == 3:
                    self._ti_num = i
            
            if self._ti_num == None:
                self.time_intervals = self.time_intervals.tolist()
                self.time_intervals.append(self.toi)
                np.save("./"+self._address_Xspec+"-ti", self.time_intervals)
                self.time_intervals = np.asarray(self.time_intervals)
                self._ti_num = len(self.time_intervals)-1                

        self._m_num = 0
        self.model=""

        self._energy_start = eRan[0]
        self._energy_end = eRan[1]
        
        self._irf_class = 'P8R2_TRANSIENT020_V6'
        self._evt_class = 16
        self._zmax = 105
        self._rad = 12

