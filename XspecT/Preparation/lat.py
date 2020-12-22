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
gtltcube = GtApp('gtltcube')
gtlike = GtApp('gtlike')

from .make4FGLxml import *
from XspecT.info import EventInfo


class PrepLAT(EventInfo):
    
    
    def __init__(self, grb, tRan = [None, None], zmax = 105, rad = 12, verbose=False, **kwargs):
        
        if hasattr(grb, 'dtype'):
            grb = grb.decode('UTF-8')
        super().__init__(grb, tRan = tRan, **kwargs)

        self._zmax = zmax
        self._rad = rad
        self.verbose=verbose

        self.Preparation()
        
    def Preparation(self):
        self.__EVSCfiles__()
        self.__filter_data__()
        self.__maketime_data__()
        self.__check_file__()
        self.__evtbin_LAT__()
        self.__rspgen_data__()
        #self.__expCube_data__()
        self.__ltCube_data__()
        self.__generateXML__()
        self.__expMap_data__()
        
        self.__srcProb__()
        self.__bkGenLAT__()
        
        
    def __EVSCfiles__(self):
        os.system("ls ./{}/LAT/*EV* > ./{}/LAT/{}-events.txt".format(self.fullname, self.fullname,self.GRB_name))

    def __check_file__(self):
        Nevts = len(fits.open("./090720710/LAT/090720-mktime.fits")[1].data)
        if Nevts == 0:
            if self._zmax!=110 and self._rad!=8:
                if self.verbose:
                    print("With the current selection (zmax={} and roi={}), there is no LAT event.")
                    print("We increase the selection; i.e., a losse cut, zmax = 110 and roi = 8.")
                self._zmax = 110
                self._rad = 8
                self.__filter_data__()
                self.__maketime_data__()
            else:
                if self.verbose:
                    print("With the loose selection, there is no LAT event.")
                    print("We may need a caution on analyzing this burst.")


    def __filter_data__(self, tRan = None, srcProb=False): 
        filter['infile'] = './{}-events.txt'.format(self._address_LAT)
        filter['ra'] = self._refRa
        filter['dec'] = self._refDec
        filter['rad'] = self._rad
        if srcProb and tRan!=None:
            filter['tmin'] = self.Trigger+tRan[0]
            filter['tmax'] = self.Trigger+tRan[1]
            filter['outfile'] = './{}-filtered_temp.fits'.format(self._address_LAT)
        else:
            filter['tmin'] = self._tStart+self.Trigger
            filter['tmax'] = self._tEnd+self.Trigger
            filter['outfile'] = './{}-filtered.fits'.format(self._address_LAT)
        filter['emin'] = self._eStart
        filter['emax'] = self._eEnd
        filter['zmax'] = self._zmax
        filter['evclass'] = self._eC
        filter['evtype'] = 3
        filter['chatter'] = int(self.verbose)
        filter.run(print_command=self.verbose)

    def __maketime_data__(self, srcProb = False):
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
        #maketime['tstart'] = self._tStart+self.Trigger
        #maketime['tstop'] = self._tEnd+self.Trigger
        maketime['chatter'] = int(self.verbose)
        maketime.run(print_command=self.verbose)
        

    def __evtbin_LAT__(self):
        evtbin['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        evtbin['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        evtbin['outfile'] = './{}-LAT.pha'.format(self._address_Xspec)
        evtbin['algorithm'] = 'PHA1'
        evtbin['ebinalg'] = 'LOG'
        evtbin['emin'] = self._eStart
        evtbin['emax'] = self._eEnd
        evtbin['enumbins'] = 10
        evtbin['tbinalg'] = 'LIN'
        evtbin['tstart'] = self._tStart+self.Trigger
        evtbin['tstop'] = self._tEnd+self.Trigger
        evtbin['dtime'] = self._tEnd-self._tStart
        evtbin['chatter'] = int(self.verbose)
        evtbin['snratio'] = 0.
        evtbin['lcemin'] = 0.
        evtbin['lcemax'] = 0.
        evtbin['nxpix'] = 240
        evtbin['nypix'] = 240
        evtbin['binsz'] = 0.1
        evtbin['coordsys'] = 'CEL'
        evtbin['xref'] = self._refRa
        evtbin['yref'] = self._refDec
        evtbin.run(print_command=self.verbose)
        
    def __rspgen_data__(self):
        rspgen['respalg'] = 'PS'
        rspgen['specfile'] = './{}-LAT.pha'.format(self._address_Xspec)
        rspgen['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        rspgen['outfile'] = './{}-LAT.rsp'.format(self._address_Xspec)
        rspgen['thetacut'] = 90
        rspgen['irfs'] = self._irfC
        rspgen['emin'] = self._eStart
        rspgen['emax'] = self._eEnd
        rspgen['enumbins'] = 10
        rspgen['chatter'] = int(self.verbose)
        rspgen.run(print_command=self.verbose) 
            
    def __expCube_data__(self):
        expCube['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        expCube['outfile'] = './{}-ltcube.fits'.format(self._address_LAT)
        expCube['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        expCube['dcostheta'] = 0.025
        expCube['binsz'] = 1
        expCube['zmax'] = self._zmax
        expCube['chatter'] = int(self.verbose)
        expCube.run(print_command=self.verbose)
    
    def __ltCube_data__(self):
        gtltcube['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        gtltcube['outfile'] = './{}-ltcube.fits'.format(self._address_LAT)
        gtltcube['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        gtltcube['dcostheta'] = 0.025
        gtltcube['binsz'] = 1
        gtltcube['tmin'] = self._tStart+self.Trigger
        gtltcube['tmax'] = self._tEnd+self.Trigger
        gtltcube['zmax'] = self._zmax
        gtltcube['chatter'] = int(self.verbose)
        gtltcube.run(print_command=self.verbose)

    def __expMap_data__(self):
        expMap['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        expMap['outfile'] = './{}-expmap.fits'.format(self._address_LAT)
        expMap['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        expMap['expcube'] = './{}-ltcube.fits'.format(self._address_LAT)
        expMap['irfs'] = self._irfC
        expMap['srcrad'] = 30
        expMap['nlong'] = 120
        expMap['nlat'] = 120
        expMap['nenergies'] = 10
        expMap['chatter'] = int(self.verbose)
        expMap.run(print_command=self.verbose)

    def __diffResps_data__(self, srcProb=False):
        if srcProb:
            diffResps['evfile'] = './{}-mktime_temp.fits'.format(self._address_LAT)
        else:
            diffResps['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        diffResps['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        diffResps['srcmdl'] = './{}-model.xml'.format(self._address_LAT)
        diffResps['irfs'] = self._irfC
        diffResps['evclass'] = 16
        diffResps['evtype'] = 3
        diffResps['chatter'] = int(self.verbose)
        diffResps.run(print_command=self.verbose)

    def __generateXML__(self):
        FERMI_DIR = os.environ.get('FERMI_DIR')
        SCRIPT_DIR = os.environ.get('SCRIPT_DIR')
        mymodel = srcList('{}/RefData/gll_psc_v26.fit'.format(SCRIPT_DIR), ft1='./{}-mktime.fits'.format(self._address_LAT), out='./{}-model.xml'.format(self._address_LAT))
        mymodel.makeModel(radLim=0.001, GRB=True, GDname = "Galactic diffuse", GDfile = '{}/refdata/fermi/galdiffuse/gll_iem_v07.fits'.format(FERMI_DIR),  ISOname = "Isotropic diffuse", ISOfile = '{}/refdata/fermi/galdiffuse/iso_{}_v06.txt'.format(FERMI_DIR, self._irfC), varFree=False)

    def __bkGenLAT__(self):
        self.__diffResps_data__()
        gtbkg['phafile'] = './{}-LAT.pha'.format(self._address_Xspec)
        gtbkg['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        gtbkg['outfile'] = './{}-LAT.bak'.format(self._address_Xspec)
        gtbkg['irfs'] = self._irfC
        gtbkg['expcube'] = './{}-ltcube.fits'.format(self._address_LAT)
        gtbkg['expmap'] = './{}-expmap.fits'.format(self._address_LAT)
        gtbkg['srcmdl'] = './{}-model.xml'.format(self._address_LAT)
        gtbkg['target'] = 'GRB'
        gtbkg['debug'] = True
        gtbkg['chatter'] = int(self.verbose)
        gtbkg.run(print_command=self.verbose) 

    def __likelihood__(self):
        gtlike['evfile'] = './{}-mktime.fits'.format(self._address_LAT)
        gtlike['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        gtlike['expcube'] = './{}-ltcube.fits'.format(self._address_LAT)
        gtlike['expmap'] = './{}-expmap.fits'.format(self._address_LAT)
        gtlike['irfs'] = self._irfC
        gtlike['srcmdl'] = './{}-model.xml'.format(self._address_LAT)
        gtlike['optimizer'] = "NEWMINUIT"
        gtlike['chatter'] = int(self.verbose)
        gtlike.run(print_command=self.verbose) 

    def __srcProb__(self):
        self.__filter_data__(tRan = [-100, 1000], srcProb=True)
        self.__maketime_data__(srcProb=True)
        self.__diffResps_data__(srcProb=True)
        
        gtsrcprob["evfile"] = './{}-mktime_temp.fits'.format(self._address_LAT)
        gtsrcprob["scfile"] = './{}-SC00.fits'.format(self._address_LAT)
        gtsrcprob["outfile"] = './{}-srcprob.fits'.format(self._address_LAT)
        gtsrcprob["srcmdl"] = './{}-model.xml'.format(self._address_LAT)
        gtsrcprob["irfs"] = self._irfC
        gtsrcprob['chatter'] = int(self.verbose)
        gtsrcprob.run(print_command=self.verbose) 

        os.system("rm './{}-filtered_temp.fits'".format(self._address_LAT))
        os.system("rm './{}-mktime_temp.fits'".format(self._address_LAT))