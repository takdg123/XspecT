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
    
    
    def __init__(self, grb, tRan = [None, None], eRan = [100, 10000], zmax = 105, rad = 12, verbose=False, **kwargs):
        
        super().__init__(grb, tRan = tRan, eRan=eRan, **kwargs)

        self._zmax = zmax
        self._rad = rad
        self.verbose=verbose

        self.Preparation()
        
    def Preparation(self):
        self.__Files__()
        self.__FilterData__()
        self.__MaketimeData__()
        self.__CheckFile__()
        self.__EvtbinData__()
        self.__RspgenData__()
        self.__ExpCubeData__()
        self.__GenerateXML__()
        self.__ExpMapData__()
        self.__SrcProb__()
        self.__BkgData__()
        
        
    def __Files__(self):
        os.system("ls ./{}/LAT/*EV* > ./{}/LAT/{}-events.txt".format(self.full_name, self.full_name,self.event_name))

    def __CheckFile__(self):
        Nevts = len(fits.open("./090720710/LAT/090720-mktime.fits")[1].data)
        if Nevts == 0:
            if self._zmax!=110 and self._rad!=8:
                if self.verbose:
                    print("With the current selection (zmax={} and roi={}), there is no LAT event.")
                    print("We increase the selection; i.e., a losse cut, zmax = 110 and roi = 8.")
                self._zmax = 110
                self._rad = 8
                self.__FilterData__()
                self.__MaketimeData__()
            else:
                if self.verbose:
                    print("With the loose selection, there is no LAT event.")
                    print("We may need a caution on analyzing this burst.")


    def __FilterData__(self, tRan = None, srcProb=False): 
        filter['infile'] = './{}-events.txt'.format(self._address_LAT)
        filter['ra'] = self.refRa
        filter['dec'] = self.refDec
        filter['rad'] = self._rad
        if srcProb and tRan!=None:
            filter['tmin'] = self.trigger+tRan[0]
            filter['tmax'] = self.trigger+tRan[1]
            filter['outfile'] = './{}-filtered_temp.fits'.format(self._address_LAT)
        else:
            filter['tmin'] = self._time_start+self.trigger
            filter['tmax'] = self._time_end+self.trigger
            filter['outfile'] = './{}-filtered.fits'.format(self._address_LAT)
        filter['emin'] = self._energy_start
        filter['emax'] = self._energy_end
        filter['zmax'] = self._zmax
        filter['evclass'] = self._evt_class
        filter['evtype'] = 3
        filter['chatter'] = int(self.verbose)
        filter.run(print_command=self.verbose)

    def __MaketimeData__(self, srcProb = False):
        if srcProb:
            maketime['evfile'] = './{}-filtered_temp.fits'.format(self._address_LAT)
            maketime['outfile'] = './{}-mktime_temp.fits'.format(self._address_LAT)
        else:
            maketime['evfile'] = './{}-filtered.fits'.format(self._address_LAT)
            maketime['outfile'] = './{}-mktime.fits'.format(self._address_LAT)
        maketime['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        maketime['filter'] = '(DATA_QUAL>0||DATA_QUAL==-1||DATA_QUAL==1)&&(LAT_CONFIG==1)'
        maketime['apply_filter'] = 'yes'
        maketime['roicut'] = 'yes'
        maketime['chatter'] = int(self.verbose)
        maketime.run(print_command=self.verbose)
        

    def __EvtbinData__(self):
        evtbin['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        evtbin['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        evtbin['outfile'] = './{}-LAT.pha'.format(self._address_Xspec)
        evtbin['algorithm'] = 'PHA1'
        evtbin['ebinalg'] = 'LOG'
        evtbin['emin'] = self._energy_start
        evtbin['emax'] = self._energy_end
        evtbin['enumbins'] = 10
        evtbin['tbinalg'] = 'LIN'
        evtbin['tstart'] = self._time_start+self.trigger
        evtbin['tstop'] = self._time_end+self.trigger
        evtbin['dtime'] = self._time_end-self._time_start
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
        rspgen['specfile'] = './{}-LAT.pha'.format(self._address_Xspec)
        rspgen['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        rspgen['outfile'] = './{}-LAT.rsp'.format(self._address_Xspec)
        rspgen['thetacut'] = 90
        rspgen['irfs'] = self._irf_class
        rspgen['emin'] = self._energy_start
        rspgen['emax'] = self._energy_end
        rspgen['enumbins'] = 10
        rspgen['chatter'] = int(self.verbose)
        rspgen.run(print_command=self.verbose) 
            
    def __ExpCubeData__(self):
        expCube['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        expCube['outfile'] = './{}-ltcube.fits'.format(self._address_LAT)
        expCube['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        expCube['dcostheta'] = 0.025
        expCube['binsz'] = 1
        expCube['zmax'] = self._zmax
        expCube['chatter'] = int(self.verbose)
        expCube.run(print_command=self.verbose)
 
    def __ExpMapData__(self):
        expMap['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        expMap['outfile'] = './{}-expmap.fits'.format(self._address_LAT)
        expMap['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        expMap['expcube'] = './{}-ltcube.fits'.format(self._address_LAT)
        expMap['irfs'] = self._irf_class
        expMap['srcrad'] = 30
        expMap['nlong'] = 120
        expMap['nlat'] = 120
        expMap['nenergies'] = 10
        expMap['chatter'] = int(self.verbose)
        expMap.run(print_command=self.verbose)

    def __DiffRespsData__(self, srcProb=False):
        if srcProb:
            diffResps['evfile'] = './{}-mktime_temp.fits'.format(self._address_LAT)
        else:
            diffResps['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        diffResps['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        diffResps['srcmdl'] = './{}-model.xml'.format(self._address_LAT)
        diffResps['irfs'] = self._irf_class
        diffResps['evclass'] = 16
        diffResps['evtype'] = 3
        diffResps['chatter'] = int(self.verbose)
        diffResps.run(print_command=self.verbose)

    def __GenerateXML__(self):
        FERMI_DIR = os.environ.get('FERMI_DIR')
        SCRIPT_DIR = os.environ.get('XspecT')
        mymodel = srcList('{}/RefData/gll_psc_v26.fit'.format(SCRIPT_DIR), ft1='./{}-mktime.fits'.format(self._address_LAT), out='./{}-model.xml'.format(self._address_LAT))
        mymodel.makeModel(radLim=0.001, GRB=True, GDname = "Galactic diffuse", GDfile = '{}/refdata/fermi/galdiffuse/gll_iem_v07.fits'.format(FERMI_DIR),  ISOname = "Isotropic diffuse", ISOfile = '{}/refdata/fermi/galdiffuse/iso_{}_v06.txt'.format(FERMI_DIR, self._irf_class), varFree=False)

    def __BkgData__(self):
        self.__DiffRespsData__()
        gtbkg['phafile'] = './{}-LAT.pha'.format(self._address_Xspec)
        gtbkg['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        gtbkg['outfile'] = './{}-LAT.bak'.format(self._address_Xspec)
        gtbkg['irfs'] = self._irf_class
        gtbkg['expcube'] = './{}-ltcube.fits'.format(self._address_LAT)
        gtbkg['expmap'] = './{}-expmap.fits'.format(self._address_LAT)
        gtbkg['srcmdl'] = './{}-model.xml'.format(self._address_LAT)
        gtbkg['target'] = 'GRB'
        gtbkg['debug'] = True
        gtbkg['chatter'] = int(self.verbose)
        gtbkg.run(print_command=self.verbose) 

    def __Likelihood__(self):
        gtlike['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        gtlike['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        gtlike['expcube'] = './{}-ltcube.fits'.format(self._address_LAT)
        gtlike['expmap'] = './{}-expmap.fits'.format(self._address_LAT)
        gtlike['irfs'] = self._irf_class
        gtlike['srcmdl'] = './{}-model.xml'.format(self._address_LAT)
        gtlike['optimizer'] = "NEWMINUIT"
        gtlike['chatter'] = int(self.verbose)
        gtlike.run(print_command=self.verbose) 

    def __SrcProb__(self):
        self.__FilterData__(tRan = [-100, 1000], srcProb=True)
        self.__MaketimeData__(srcProb=True)
        self.__DiffRespsData__(srcProb=True)
        
        gtsrcprob["evfile"] = './{}-mktime_temp.fits'.format(self._address_LAT)
        gtsrcprob["scfile"] = './{}-SC00.fits'.format(self._address_LAT)
        gtsrcprob["outfile"] = './{}-srcprob.fits'.format(self._address_LAT)
        gtsrcprob["srcmdl"] = './{}-model.xml'.format(self._address_LAT)
        gtsrcprob["irfs"] = self._irf_class
        gtsrcprob['chatter'] = int(self.verbose)
        gtsrcprob.run(print_command=self.verbose) 

        os.system("rm './{}-filtered_temp.fits'".format(self._address_LAT))
        os.system("rm './{}-mktime_temp.fits'".format(self._address_LAT))