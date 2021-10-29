import os
import numpy as np
from astropy.io import fits
from XspecT.Preparation.downinfo import DownloadGBMcat

class EventInfo(object):
    _GBM_dtrs = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']
    _GBM_download_address = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts"
    _LAT_download_address = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/triggers"
    
    def __init__(self, grb, tRan = [None, None], eRan = [100, 10000], GBM_catalog=None, LAT_catalog=None, EvtInfo = {},  **kwargs):

        self.__dict__.update(**kwargs)

        if type(grb) == str or type(grb) == int or type(grb) == np.str_:
            self.importConfig(EvtInfo)
            
        elif type(grb) == dict:
            info = {**grb, **EvtInfo}
            self.importConfig(info)
            grb = self.full_name

        grb = str(grb)
        self.event_name = grb[:6]
        self.yr = grb[:2]
        self.full_name = grb
        self.current_dir = os.path.abspath('.')

        if not(os.path.isdir(self.current_dir+"/{}".format(grb))):
            self.makeDir()

        if GBM_catalog==None:
            SCRIPT_DIR=os.environ.get('XspecT')
            GBM_catalog = SCRIPT_DIR+'/RefData/GBMcatalog.fits'
        
        if LAT_catalog==None:
            SCRIPT_DIR=os.environ.get('XspecT')
            LAT_catalog = SCRIPT_DIR+'/RefData/LAT2CATALOG-v16-LTF.fits'


        self.event_info = fits.open(GBM_catalog)
        self.event_info_LAT = fits.open(LAT_catalog)

        self._info_num = self.event_info[1].data["TRIGGER_NAME"]=='BN'+self.full_name
        self._LAT_info_num = self.event_info_LAT[1].data["GRBNAME"]==self.full_name
        
        self._address_LAT = self.full_name+'/LAT/'+self.event_name
        self._address_Xspec = self.full_name+'/Xspec/'+self.event_name

        # bcat file
        if not(os.path.isfile(self.current_dir+'/{}/GBM/{}-bcat.fit'.format(self.full_name,self.event_name))):
            DownloadGBMcat(self.full_name)

        self._bcat=fits.open(self.current_dir+'/{}/GBM/{}-bcat.fit'.format(self.full_name,self.event_name))
        try:
            self._scat=fits.open(self.current_dir+'/{}/GBM/{}-scat.fit'.format(self.full_name,self.event_name))
        except:
            self._scat = None

        # GBM detectors
        if not(hasattr(self, "usedGBM")):
            
            if sum(self._info_num) == 1:
                detNm = self.event_info[1].data["SCAT_DETECTOR_MASK"][self._info_num]
                
                for det, i in zip(detNm[0], range(len(detNm[0]))):
                    if int(det):
                        self.usedGBM.append(self._GBM_dtrs[i])

            if len(self.usedGBM) == 0:
                if self._scat!=None:
                    for dtr in self._scat[1].data["DETNAM"]:
                        dtr=dtr.replace("NAI_0", "n")
                        dtr=dtr.replace("NAI_10", "na")
                        dtr=dtr.replace("NAI_11", "nb")
                        dtr=dtr.replace("BGO_0", "b")
                        self.usedGBM.append(dtr)
                else:
                    print('scat file does not exist. "usedGBM" should be specified.')
        
        # Trigger time
        if not(hasattr(self, "trigger")):    
            self.trigger = self._bcat[0].header["TRIGTIME"]
            if self._bcat[0].header["T90START"] - self.event_info[1].data["T90_START"][self._info_num ][0] > 1e-2:
                print("T90_start Error")
                print(self._bcat[0].header["T90START"] - self.event_info[1].data["T90_START"][self._info_num])
            if self._bcat[0].header["T90"] - self.event_info[1].data["T90"][self._info_num][0] > 1e-2:
                print("T90 Error")
                print(self._bcat[0].header["T90"], self.event_info[1].data["T90"][self._info_num])
        
        # T90
        if not(hasattr(self, "t90")): 
            self.t90 = float(self._bcat[0].header["T90"])
            self.t05 = float(self._bcat[0].header["T90START"])
            self.t95 = self._bcat[0].header["T90START"]+self._bcat[0].header["T90"]

        # Localization
        if not(hasattr(self, "refRa")) and not(hasattr(self, "refDec")):
            if sum(self._LAT_info_num) == 0:
                print("This GRB is not in the LAT catalog. Please update the catalog file.")
                print("You need to define following parameters.")
                print("refRa, refDec")

            elif sum(self._LAT_info_num) == 1:
                self.t05 = self.event_info_LAT[1].data["GBMT05"][self._LAT_info_num][0]
                self.t95 = self.event_info_LAT[1].data["GBMT95"][self._LAT_info_num][0]
                self._tl0 = self.event_info_LAT[1].data["TL0"][self._LAT_info_num][0]
                self._tl1 = self.event_info_LAT[1].data["TL1"][self._LAT_info_num][0]
                self.refRa = self.event_info_LAT[1].data["RA"][self._LAT_info_num][0]
                self.refDec = self.event_info_LAT[1].data["DEC"][self._LAT_info_num][0]

        # Time range
        if tRan[0] !=None:
            if tRan[0] == "T0": _time_start = 0.
            elif tRan[0] == "T05": _time_start = self.t05
            elif tRan[0] == "T95": _time_start = self.t95
            elif tRan[0] == "TL0" and sum(self._LAT_info_num) == 1:
                _time_start = self._tl0
            elif tRan[0] == "TL1" and sum(self._LAT_info_num) == 1:
                _time_start = self._tl1
            elif isinstance(tRan[0], int) or isinstance(tRan[0], float):
                _time_start = tRan[0]
            else:
                raise("Check your start time (tRan[0]): e.g., T05, T95, TL0, TL1, float, int, ...")
        else:
            _time_start = self.t05


        if tRan[1] !=None:
            if tRan[1] == "T0": _time_end = 0.
            elif tRan[1] == "T05": _time_end = self.t05
            elif tRan[1] == "T95": _time_end = self.event_info_LAT[1].data["GBMT95"][self._LAT_info_num][0]         
            elif tRan[1] == "TL0" and sum(self._LAT_info_num) == 1:
                _time_end = self._tl0
            elif tRan[1] == "TL1" and sum(self._LAT_info_num) == 1:
                _time_end = self._tl1
            elif isinstance(tRan[1], int) or isinstance(tRan[1], float):
                _time_end = tRan[1]
            else:
                raise("Check your end time (tRan[1]): e.g., T05, T95, TL0, TL1, float, int, ...")
        else:
            _time_end = self.t95

        if _time_end < _time_start:
            _time_start, _time_end = _time_end, _time_start 

        self.toi = [_time_start, _time_end, self.trigger]
        
        # Setting time intervals
        if hasattr(self, "time_intervals"):
            self._ti_num = 0
            self.toi = [self.time_intervals[0][0], self.time_intervals[0][1]]
            np.save(self.current_dir+"/"+self._address_Xspec+"-ti", self.time_intervals)
        else:
            if not(os.path.isfile(self._address_Xspec+"-ti.npy")):
                self.time_intervals = np.asarray([self.toi])
                self._ti_num = 0
                np.save(self.current_dir+"/"+self._address_Xspec+"-ti", self.time_intervals)
            else:
                self.time_intervals = np.load(self.current_dir+"/"+self._address_Xspec+"-ti.npy")
                
                self._ti_num = None
                
                for i in range(len(self.time_intervals)):
                    if sum(self.toi == self.time_intervals[i]) == 3:
                        self._ti_num = i
                
                if self._ti_num == None:
                    self.time_intervals = self.time_intervals.tolist()
                    self.time_intervals.append(self.toi)
                    np.save(self.current_dir+"/"+self._address_Xspec+"-ti", self.time_intervals)
                    self.time_intervals = np.asarray(self.time_intervals)
                    self._ti_num = len(self.time_intervals)-1                

        self._m_num = 0
        self.model=""

        self._energy_start = eRan[0]
        self._energy_end = eRan[1]
        
        self._irf_class = 'P8R3_TRANSIENT020E_V3'
        self._evt_class = 8
        self._zmax = 105
        self._rad = 12


    def importConfig(self, config):
        keys = config.keys()
        if "FullName" in keys:      self.full_name = config["FullName"]
        if "Name" in keys:          self.event_name = config['Name']
        if "Trigger" in keys:       self.trigger = config['Trigger']
        if "T90" in keys:           
            self.t05 = config['T90'][0]
            self.t95 = config['T90'][1]
            self.t90 = config['T90'][1] - config['T90'][0]
        if "RA" in keys:            self.refRa = config['RA']
        if "DEC" in keys:           self.refDec = config['DEC']
        if "usedGBM" in keys:       self.usedGBM = config['usedGBM']
    
        if 'redshift' in keys:      self.redshift = config['redshift']
        if 'nH' in keys:            self._nH = config['nH']
        
        if 'only' in keys:          self._only = config['only']

        if ("t_start" in keys) and ("t_end" in keys):
            self.time_intervals = []
            for ti, tf in zip(config["t_start"],config["t_end"]) :
                self.time_intervals.append([ti, tf])
            self.time_intervals = np.asarray(self.time_intervals)

        if "TimeBins" in keys:
            self.time_interal = [config["TimeBins"][:-1], config["TimeBins"][1:]]
            self.time_interal = np.asarray(self.time_interal).T


    def makeDir(self):
        os.system("mkdir {}".format(self.full_name))
        os.system("mkdir {}/GBM".format(self.full_name))
        os.system("mkdir {}/Xspec".format(self.full_name))
        os.system("mkdir {}/LAT".format(self.full_name))