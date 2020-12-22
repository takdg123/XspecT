import os
from astropy.io import fits

class EventInfo(object):
    _gbmDet = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']
    _GBMdownAdd = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts"
    _LATdownAdd = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/triggers"
    
    def __init__(self, grb, tRan = [None, None], GBM_catalog=None, LAT_catalog=None, **kwargs):

        if hasattr(grb, 'dtype'):
            grb = grb.decode('UTF-8')

        self.__dict__.update(**kwargs)


        if GBM_catalog==None:
            SCRIPT_DIR=os.environ.get('SCRIPT_DIR')
            GBM_catalog = SCRIPT_DIR+'/RefData/GBMcatalog.fits'
        
        if LAT_catalog==None:
            SCRIPT_DIR=os.environ.get('SCRIPT_DIR')
            LAT_catalog = SCRIPT_DIR+'/RefData/LAT2CATALOG-v16-LTF.fits'
        
        self.GRBinfo = fits.open(GBM_catalog)
        self.LATGRBinfo = fits.open(LAT_catalog)
        
        self.GRB_name = grb[:6]
        self.yr = grb[:2]
        self.fullname = grb
        
        self.usedGBM = []
        self.infoNum = self.GRBinfo[1].data["TRIGGER_NAME"]=='BN'+self.fullname
        self.LATinfoNum = self.LATGRBinfo[1].data["GRBNAME"]==self.fullname
        
        detNm = self.GRBinfo[1].data["SCAT_DETECTOR_MASK"][self.infoNum]
        
        for det, i in zip(detNm[0], range(len(detNm[0]))):
            if int(det):
                self.usedGBM.append(self._gbmDet[i])
        
        try:
            self.info=fits.open('./{}/GBM/{}-bcat.fit'.format(self.fullname,self.GRB_name))
        except:
            self.makeDir()
            self.download(initial=True)
            self.info=fits.open('./{}/GBM/{}-bcat.fit'.format(self.fullname,self.GRB_name))
            
        self.Trigger = self.info[0].header["TRIGTIME"]
        if self.info[0].header["T90START"] - self.GRBinfo[1].data["T90_START"][self.infoNum ][0] > 1e-2:
            print("T90_start Error")
            print(self.info[0].header["T90START"] - self.GRBinfo[1].data["T90_START"][self.infoNum])
        if self.info[0].header["T90"] - self.GRBinfo[1].data["T90"][self.infoNum][0] > 1e-2:
            print("T90 Error")
            print(self.info[0].header["T90"], self.GRBinfo[1].data["T90"][self.infoNum])
        
        self.T90 = [self.info[0].header["T90START"],self.info[0].header["T90START"]+self.info[0].header["T90"],self.info[0].header["T90"]]
        
        self._refRa = self.LATGRBinfo[1].data["RA"][self.LATinfoNum][0]
        self._refDec = self.LATGRBinfo[1].data["DEC"][self.LATinfoNum][0]

        self._address_GBM = self.fullname+'/GBM/'+self.GRB_name
        self._address_Xspec = self.fullname+'/Xspec/'+self.GRB_name
        self._address_LAT = self.fullname+'/LAT/'+self.GRB_name

        if tRan[0] !=None:
            if tRan[0] == "T0":
                self._tStart = 0.
            elif tRan[0] == "T05":
                self._tStart = self.LATGRBinfo[1].data["GBMT05"][self.LATinfoNum][0]
            elif tRan[0] == "T95":
                self._tStart = self.LATGRBinfo[1].data["GBMT95"][self.LATinfoNum][0]         
            elif tRan[0] == "LATT0":
                self._tStart = self.LATGRBinfo[1].data["TL0"][self.LATinfoNum][0]         
            elif tRan[0] == "LATT1":
                self._tStart = self.LATGRBinfo[1].data["TL100"][self.LATinfoNum][0]
            elif isinstance(tRan[0], int) or isinstance(tRan[0], float):
                self._tStart = tRan[0]
            else:
                raise("Check your start time (tRan[0])")
        else:
            self._tStart = self.LATGRBinfo[1].data["GBMT05"][self.LATinfoNum][0]         


        if tRan[1] !=None:
            if tRan[1] == "T0":
                self._tEnd = 0.
            elif tRan[1] == "T05":
                self._tEnd = self.LATGRBinfo[1].data["GBMT05"][self.LATinfoNum][0]
            elif tRan[1] == "T95":
                self._tEnd = self.LATGRBinfo[1].data["GBMT95"][self.LATinfoNum][0]         
            elif tRan[1] == "LATT0":
                self._tEnd = self.LATGRBinfo[1].data["TL0"][self.LATinfoNum][0]         
            elif tRan[1] == "LATT1":
                self._tEnd = self.LATGRBinfo[1].data["TL100"][self.LATinfoNum][0]
            elif isinstance(tRan[1], int) or isinstance(tRan[1], float):
                self._tEnd = tRan[1]
            else:
                raise("Check your end time (tRan[1])")
        else:
            self._tEnd = self.LATGRBinfo[1].data["GBMT95"][self.LATinfoNum][0]

        if self._tEnd < self._tStart:
            self._tStart, self._tEnd = self._tEnd, self._tStart 

        self._tStart = self.LATGRBinfo[1].data["GBMT05"][self.LATinfoNum][0]
        self._tEnd = self.LATGRBinfo[1].data["GBMT95"][self.LATinfoNum][0]
        self._eStart = 100
        self._eEnd = 10000
        self._irfC = 'P8R2_TRANSIENT020_V6'
        self._eC = 16
        self._zmax = 105
        self._rad = 12