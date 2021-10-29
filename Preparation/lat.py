import os, sys

from astropy.io import fits
from astropy.stats import *

from scipy import stats

import gt_apps
from gt_apps import *
from GtApp import GtApp
bindef = GtApp('gtbindef','Likelihood')
gtbkg = GtApp('gtbkg')
gtsrcprob = GtApp('gtsrcprob')
gtlike = GtApp('gtlike')

from .make4FGLxml import *
from ..info import EventInfo



class PrepLAT(EventInfo):
    
    
    def __init__(self, grb, tRan = [None, None], eRan = [100, 10000], zmax = 105, rad = 12, allBins=False, verbose=False, forcedGen=False, prep=True, **kwargs):
        
        super().__init__(grb, tRan = tRan, eRan=eRan, **kwargs)

        self._zmax = zmax
        self._rad = rad
        self.verbose=verbose
        self.forcedGen = forcedGen
        self.current_dir = os.path.abspath('.')
        path = self.current_dir+"/"+self.full_name+'/LAT/'
        files = os.listdir(path)
        self.scfile = path+[f for f in files if "SC" in f][0]

        if prep: 
            if allBins and hasattr(self, "time_intervals"):

                for i in range(len(self.time_intervals)):
                    self._ti_num = i
                    self.toi = [self.time_intervals[self._ti_num][0], self.time_intervals[self._ti_num][1]]
                
                    self.Preparation()
            else:
                self.Preparation()

    @property
    def ti_num(self):
        return self._ti_num
        
    def Preparation(self):
        self.__Files__()
        file_flag = os.path.isfile(self.current_dir+'/{}-LAT-{}.bak'.format(self._address_Xspec, self._ti_num))

        if not(self.forcedGen) and file_flag:
            print("LAT files (name: {}, time interval: {:.3f}-{:.3f}, ti_num: {}) already exist.".format(self.full_name, self.toi[0], self.toi[1], self._ti_num))
            return

        self.__FilterData__()
        try:
            self.__MaketimeData__()
        except:
            self._zmax = 110
            self._rad = 8
            self.__FilterData__()
            self.__MaketimeData__()
        status = self.__CheckFile__()
        if status == 0:
            self.__FilterData__()
            self.__MaketimeData__()
        elif status == 2:
            return
        self.__EvtbinData__()
        self.__RspgenData__()
        self.__ExpCubeData__()
        self.__GenerateXML__()
        self.__ExpMapData__()

        self.__BkgData__()
    
    def SrcProb(self, tran=[0, 1000]):
        self.__Files__()
        self.__GenerateXML__()
        self.__SrcProb__(tran=tran)
        
    def __Files__(self):
        os.system("ls {}/{}/LAT/*EV* > {}/{}/LAT/{}-events.txt".format(self.current_dir, self.full_name, self.current_dir, self.full_name,self.event_name))

    def __CheckFile__(self, file=None):
        if file == None:
            Nevts = len(fits.open(self.current_dir+"/{}-mktime-{}.fits".format(self._address_LAT, self.ti_num))[1].data)
        else:
            Nevts = len(fits.open(file)[1].data)

        if Nevts == 0:
            if self._zmax!=110 and self._rad!=8:
                if self.verbose:
                    print("With the current selection (zmax={} and roi={}), there is no LAT event.")
                    print("We increase the selection; i.e., a losse cut, zmax = 110 and roi = 8.")
                self._zmax = 110
                self._rad = 8
                
                return 0
            else:
                if self.verbose:
                    print("With the loose selection, there is no LAT event.")
                    print("We may need a caution on analyzing this burst.")
                return 2
        else:
            return 1

    def __FilterData__(self, tRan = None, srcProb=False): 
        filter['infile'] = self.current_dir+'/{}-events.txt'.format(self._address_LAT)
        filter['ra'] = self.refRa
        filter['dec'] = self.refDec
        filter['rad'] = self._rad
        if srcProb and tRan!=None:
            filter['tmin'] = self.trigger+tRan[0]
            filter['tmax'] = self.trigger+tRan[1]
            filter['outfile'] = self.current_dir+'/{}-filtered_temp.fits'.format(self._address_LAT)
        else:
            filter['tmin'] = self.toi[0]+self.trigger
            filter['tmax'] = self.toi[1]+self.trigger
            filter['outfile'] = self.current_dir+'/{}-filtered-{}.fits'.format(self._address_LAT, self.ti_num)
        filter['emin'] = self._energy_start
        filter['emax'] = self._energy_end
        filter['zmax'] = self._zmax
        filter['evclass'] = self._evt_class
        filter['evtype'] = 3
        filter['chatter'] = int(self.verbose)
        filter.run(print_command=self.verbose)

    def __MaketimeData__(self, srcProb = False):
        if srcProb:
            maketime['evfile'] = self.current_dir+'/{}-filtered_temp.fits'.format(self._address_LAT)
            maketime['outfile'] = self.current_dir+'/{}-mktime_temp.fits'.format(self._address_LAT)
        else:
            maketime['evfile'] = self.current_dir+'/{}-filtered-{}.fits'.format(self._address_LAT, self.ti_num)
            maketime['outfile'] = self.current_dir+'/{}-mktime-{}.fits'.format(self._address_LAT, self.ti_num)
        maketime['scfile'] = self.scfile
        maketime['filter'] = '(DATA_QUAL>0||DATA_QUAL==-1||DATA_QUAL==1)&&(LAT_CONFIG==1)'
        maketime['apply_filter'] = 'yes'
        maketime['roicut'] = 'yes'
        maketime['chatter'] = int(self.verbose)
        maketime.run(print_command=self.verbose)
        

    def __EvtbinData__(self):
        evtbin['evfile'] = self.current_dir+'/{}-mktime-{}.fits'.format(self._address_LAT, self.ti_num)
        evtbin['scfile'] = self.scfile
        evtbin['outfile'] = self.current_dir+'/{}-LAT-{}.pha'.format(self._address_Xspec, self.ti_num)
        evtbin['algorithm'] = 'PHA1'
        evtbin['ebinalg'] = 'LOG'
        evtbin['emin'] = self._energy_start
        evtbin['emax'] = self._energy_end
        evtbin['enumbins'] = 10
        evtbin['tbinalg'] = 'LIN'
        evtbin['tstart'] = self.toi[0]+self.trigger
        evtbin['tstop'] = self.toi[1]+self.trigger
        evtbin['dtime'] = self.toi[1]-self.toi[0]
        evtbin['chatter'] = int(self.verbose)
        evtbin['snratio'] = 0.
        evtbin['lcemin'] = 0.
        evtbin['lcemax'] = 0.
        evtbin['nxpix'] = 240
        evtbin['nypix'] = 240
        evtbin['binsz'] = 0.1
        evtbin['coordsys'] = 'CEL'
        evtbin['xref'] = self.refRa
        evtbin['yref'] = self.refDec
        evtbin.run(print_command=self.verbose)
        
    def __RspgenData__(self):
        rspgen['respalg'] = 'PS'
        rspgen['specfile'] = self.current_dir+'/{}-LAT-{}.pha'.format(self._address_Xspec, self.ti_num)
        rspgen['scfile'] = self.scfile
        rspgen['outfile'] = self.current_dir+'/{}-LAT-{}.rsp'.format(self._address_Xspec, self.ti_num)
        rspgen['thetacut'] = 90
        rspgen['irfs'] = self._irf_class
        rspgen['emin'] = self._energy_start
        rspgen['emax'] = self._energy_end
        rspgen['enumbins'] = 10
        rspgen['chatter'] = int(self.verbose)
        rspgen.run(print_command=self.verbose) 
            
    def __ExpCubeData__(self):
        expCube['evfile'] = self.current_dir+'/{}-mktime-{}.fits'.format(self._address_LAT, self.ti_num)
        expCube['outfile'] = self.current_dir+'/{}-ltcube-{}.fits'.format(self._address_LAT, self.ti_num)
        expCube['scfile'] = self.scfile
        expCube['dcostheta'] = 0.025
        expCube['binsz'] = 1
        expCube['zmax'] = self._zmax
        expCube['chatter'] = int(self.verbose)
        expCube['clobber'] = "yes"
        expCube.run(print_command=self.verbose)
 
    def __ExpMapData__(self):
        expMap['evfile'] = self.current_dir+'/{}-mktime-{}.fits'.format(self._address_LAT, self.ti_num)
        expMap['outfile'] = self.current_dir+'/{}-expmap-{}.fits'.format(self._address_LAT, self.ti_num)
        expMap['scfile'] = self.scfile
        expMap['expcube'] = self.current_dir+'/{}-ltcube-{}.fits'.format(self._address_LAT, self.ti_num)
        expMap['irfs'] = self._irf_class
        expMap['srcrad'] = 30
        expMap['nlong'] = 120
        expMap['nlat'] = 120
        expMap['nenergies'] = 10
        expMap['chatter'] = int(self.verbose)
        expMap.run(print_command=self.verbose)

    def __DiffRespsData__(self, srcProb=False):
        if srcProb:
            diffResps['evfile'] = self.current_dir+'/{}-mktime_temp.fits'.format(self._address_LAT)
        else:
            diffResps['evfile'] = self.current_dir+'/{}-mktime-{}.fits'.format(self._address_LAT, self.ti_num)
        diffResps['scfile'] = self.scfile
        diffResps['srcmdl'] = self.current_dir+'/{}-model.xml'.format(self._address_LAT)
        diffResps['irfs'] = self._irf_class
        diffResps['evclass'] = self._evt_class
        diffResps['evtype'] = 3
        diffResps['chatter'] = int(self.verbose)
        diffResps.run(print_command=self.verbose)

    def __GenerateXML__(self):
        FERMI_DIR = os.environ.get('FERMI_DIR')
        SCRIPT_DIR = os.environ.get('XspecT')
        mymodel = srcList('{}/RefData/gll_psc_v26.fit'.format(SCRIPT_DIR), ft1=self.current_dir+'/{}-mktime-{}.fits'.format(self._address_LAT, self.ti_num), out=self.current_dir+'/{}-model.xml'.format(self._address_LAT))
        mymodel.makeModel(radLim=0.001, GDname = "Galactic diffuse", GDfile = '{}/refdata/fermi/galdiffuse/gll_iem_v07.fits'.format(FERMI_DIR),  ISOname = "Isotropic diffuse", ISOfile = '{}/refdata/fermi/galdiffuse/iso_{}_v1.txt'.format(FERMI_DIR, self._irf_class), ExtraRad=0, varFree=False)

    def __BkgData__(self):
        self.__DiffRespsData__()
        gtbkg['phafile'] = self.current_dir+'/{}-LAT-{}.pha'.format(self._address_Xspec, self.ti_num)
        gtbkg['scfile'] = self.scfile
        gtbkg['outfile'] = self.current_dir+'/{}-LAT-{}.bak'.format(self._address_Xspec, self.ti_num)
        gtbkg['irfs'] = self._irf_class
        gtbkg['expcube'] = self.current_dir+'/{}-ltcube-{}.fits'.format(self._address_LAT, self.ti_num)
        gtbkg['expmap'] = self.current_dir+'/{}-expmap-{}.fits'.format(self._address_LAT, self.ti_num)
        gtbkg['srcmdl'] = self.current_dir+'/{}-model.xml'.format(self._address_LAT)
        gtbkg['target'] = 'GRB'
        gtbkg['debug'] = True
        gtbkg['chatter'] = int(self.verbose)
        gtbkg.run(print_command=self.verbose) 

    def __Likelihood__(self):
        gtlike['evfile'] = self.current_dir+'/{}-mktime-{}.fits'.format(self._address_LAT, self.ti_num)
        gtlike['scfile'] = self.scfile
        gtlike['expcube'] = self.current_dir+'/{}-ltcube-{}.fits'.format(self._address_LAT, self.ti_num)
        gtlike['expmap'] = self.current_dir+'/{}-expmap-{}.fits'.format(self._address_LAT, self.ti_num)
        gtlike['irfs'] = self._irf_class
        gtlike['srcmdl'] = self.current_dir+'/{}-model.xml'.format(self._address_LAT)
        gtlike['optimizer'] = "NEWMINUIT"
        gtlike['chatter'] = int(self.verbose)
        gtlike.run(print_command=self.verbose) 

    def __SrcProb__(self, tran = None):
        self._zmax = 105
        self._rad = 12

        if tran == None:
            if hasattr(self, "_tl0"):
                tRan = [self._tl0-1, min(self._tl1+1, 1000)]
            else:
                tRan = [-100, 1000]
        else:
            tRan=tran
        
        self.__FilterData__(tRan = tRan, srcProb=True)
        try:
            self.__MaketimeData__(srcProb=True)
        except:
            self._zmax = 110
            self._rad = 8
            self.__FilterData__(tRan = tRan, srcProb=True)
            self.__MaketimeData__(srcProb=True)

        status = self.__CheckFile__(file=self.current_dir+'/{}-mktime_temp.fits'.format(self._address_LAT))
        if status == 0:
            self.__FilterData__(tRan = tRan, srcProb=True)
            self.__MaketimeData__(srcProb=True)            

        self.__DiffRespsData__(srcProb=True)
        
        
        gtsrcprob["evfile"] = self.current_dir+'/{}-mktime_temp.fits'.format(self._address_LAT)
        gtsrcprob["scfile"] = self.scfile
        gtsrcprob["outfile"] = self.current_dir+'/{}-srcprob.fits'.format(self._address_LAT)
        gtsrcprob["srcmdl"] = self.current_dir+'/{}-model.xml'.format(self._address_LAT)
        gtsrcprob["irfs"] = self._irf_class
        gtsrcprob['chatter'] = int(self.verbose)
        gtsrcprob.run(print_command=self.verbose) 

        os.system("rm self.current_dir+'/{}-filtered_temp.fits'".format(self._address_LAT))
        os.system("rm self.current_dir+'/{}-mktime_temp.fits'".format(self._address_LAT))