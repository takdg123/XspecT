import numpy as np
import os, sys
import urllib
import time

from astropy.io import fits
from GtBurst import html2text

import re

from .downinfo import DownloadGBMcat
from ..info import EventInfo

class DownloadData(EventInfo):
    
    def __init__(self, grb, verbose=False, **kwargs):
        super().__init__(grb, **kwargs)
        self.verbose=verbose

    def DownloadAll(self, gbm=True, lle=False, lat=True, forcedDwn=False):

        self.forcedDwn=forcedDwn
        if gbm:
            self.DownloadGBM()
        if lle:
            self.DownloadLLE()
        if lat:
            self.DownloadLAT()

    def DownloadGBM(self):
        dwDet = self.usedGBM

        urllib.request.urlretrieve("{}/20{}/bn{}/quicklook/glg_lc_medres34_bn{}.gif".format(self._GBM_download_address,self.yr,self.full_name,self.full_name), "./{}/GBM/{}-lc.gif".format(self.full_name,self.event_name))

        if self.verbose: print("Checking tte file version...")
        for version in range(10):
            version=9-version
            try:
                urllib.request.urlopen(urllib.request.Request("{}/20{}/bn{}/current/glg_tte_{}_bn{}_v0{}.fit".format(self._GBM_download_address,self.yr,self.full_name,dwDet[0],self.full_name,version))) 
                break
            except:
                continue
        tteV = version

        if self.verbose: print("Checking pha file version...")              
        for version in range(10):
            version=9-version
            try:
                urllib.request.urlopen(urllib.request.Request("{}/20{}/bn{}/current/glg_cspec_{}_bn{}_v0{}.pha".format(self._GBM_download_address,self.yr,self.full_name,dwDet[0],self.full_name,version))) 
                break
            except:
                continue
        cspecV = version
        
        if self.verbose: print("Checking rsp file version...")
        rsp2Flag=False
        for version in range(10):
            version=9-version
            try:
                urllib.request.urlopen(urllib.request.Request("{}/20{}/bn{}/current/glg_cspec_{}_bn{}_v0{}.rsp2".format(self._GBM_download_address,self.yr,self.full_name,dwDet[0],self.full_name,version))) 
                rsp2Flag=True
                break
            except:
                try:
                    urllib.request.urlopen(urllib.request.Request("{}/20{}/bn{}/current/glg_cspec_{}_bn{}_v0{}.rsp".format(self._GBM_download_address,self.yr,self.full_name,dwDet[0],self.full_name,version)))
                    break
                except:
                    continue
        rspV = version
        
        if self.verbose:
            print("GRB{} data version: TTE ({}), CSPEC ({}), RSP ({})".format(self.full_name, tteV, cspecV, rspV))
            if rsp2Flag: print("                    Use rsp2 file")

        for gbm in dwDet:
            
            # tte files
            dwUrl = "{}/20{}/bn{}/current/glg_tte_{}_bn{}_v0{}.fit".format(self._GBM_download_address,self.yr,self.full_name,gbm,self.full_name, tteV)
            dwFile = "./{}/GBM/glg_tte_{}_bn{}_v00.fit".format(self.full_name,gbm,self.full_name)
            self.__DownloadProcess__(dwUrl, dwFile)

            # pha files
            dwUrl = "{}/20{}/bn{}/current/glg_cspec_{}_bn{}_v0{}.pha".format(self._GBM_download_address,self.yr,self.full_name,gbm,self.full_name, cspecV)
            dwFile = "./{}/GBM/glg_cspec_{}_bn{}_v00.pha".format(self.full_name,gbm,self.full_name)
            self.__DownloadProcess__(dwUrl, dwFile)
        
            # rsp file           
            if rsp2Flag:
                dwUrl = "{}/20{}/bn{}/current/glg_cspec_{}_bn{}_v0{}.rsp2".format(self._GBM_download_address,self.yr,self.full_name,gbm,self.full_name, rspV)
                dwFile = "./{}/GBM/glg_cspec_{}_bn{}_v00.rsp2".format(self.full_name,gbm,self.full_name)
                self.__DownloadProcess__(dwUrl, dwFile)
            else:
                dwUrl = "{}/20{}/bn{}/current/glg_cspec_{}_bn{}_v0{}.rsp".format(self._GBM_download_address,self.yr,self.full_name,gbm,self.full_name, rspV)
                dwFile = "./{}/GBM/glg_cspec_{}_bn{}_v00.rsp".format(self.full_name,gbm,self.full_name)
                self.__DownloadProcess__(dwUrl, dwFile)
        
        if self.verbose: print("Completed (GBM)")

                
    def DownloadLLE(self):

        try:
            urllib.request.urlretrieve("{}/20{}/bn{}/current/".format(self._LAT_download_address,self.yr,self.full_name))
        except:
            print("LLE data folder does not exist. Check {}/20{}/bn{}/current/".format(self._LAT_download_address,self.yr,self.full_name))
            return
            
        for version in range(12):
            version=11-version

            if version < 10:
                version = "0"+str(version)

            try:
                urllib.request.urlopen(urllib.request.Request("{}/20{}/bn{}/current/gll_selected_bn{}_v{}.fit".format(self._LAT_download_address,self.yr,self.full_name,self.full_name,version)))
                break
            except:
                continue

        # tte file
        dwUrl = "{}/20{}/bn{}/current/gll_selected_bn{}_v{}.fit".format(self._LAT_download_address,self.yr,self.full_name,self.full_name,version)
        dwFile = "./{}/LAT/gll_lle_bn{}_v00.fit".format(self.full_name,self.full_name)
        self.__DownloadProcess__(dwUrl, dwFile)
        
        # pha file
        dwUrl = "{}/20{}/bn{}/current/gll_cspec_bn{}_v{}.pha".format(self._LAT_download_address,self.yr,self.full_name,self.full_name,version)
        dwFile = "./{}/LAT/glg_cspec_ll_bn{}_v00.pha".format(self.full_name,self.full_name)
        self.__DownloadProcess__(dwUrl, dwFile)

        # rsp file
        dwUrl = "{}/20{}/bn{}/current/gll_cspec_bn{}_v{}.rsp".format(self._LAT_download_address,self.yr,self.full_name,self.full_name,version)
        dwFile = "./{}/LAT/glg_cspec_ll_bn{}_v00.rsp".format(self.full_name,self.full_name)
        self.__DownloadProcess__(dwUrl, dwFile)

        if self.verbose: print("Completed (LLE)")

        
    def DownloadLAT(self):
        if not(self.forcedDwn):
            if os.path.isfile("./{}/LAT/{}-EV00.fits".format(self.full_name,self.event_name)):
                return 
        url                         = "https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi"
        parameters                  = {}
        parameters['coordfield']    = "%s,%s" %(self.refRa,self.refDec)
        parameters['coordsystem']   = "J2000"
        parameters['shapefield']    = "%s" %(15)
        parameters['timefield']     = "%s,%s" %(self.trigger-100.,self.trigger+2000)
        parameters['timetype']      = "%s" %("MET")
        parameters['energyfield']   = "100,10000"
        parameters['photonOrExtendedOrNone'] = "Extended"
        parameters['destination']   = 'query'
        parameters['spacecraft']    = 'checked'

        if self.verbose:
            print("Query parameters:")
            for k,v in parameters.items():
                print("%30s = %s" %(k,v))

        postData                    = urllib.parse.urlencode(parameters).encode("utf-8")
        temporaryFileName           = "__temp_query_result.html"
        try:
            os.remove(temporaryFileName)
        except:
            pass
        pass

        urllib.request.urlcleanup()

        urllib.request.urlretrieve(url, temporaryFileName, lambda x,y,z:0, postData)

        htmlFile = open(temporaryFileName)
        lines = []
        for line in htmlFile:
            lines.append(line.encode('utf-8'))
        pass
        html = "".join(str(lines)).strip()
        htmlFile.close()
        if self.verbose: print("\nAnswer from the LAT data server:\n")
        text = html2text.html2text(html.strip()).split("\n")

        text = list(filter(lambda x:x.find("[") < 0 and  x.find("]") < 0 and x.find("#") < 0 and x.find("* ") < 0 and
                        x.find("+") < 0 and x.find("Skip navigation")<0,text))
        text = list(filter(lambda x:len(x.replace(" ",""))>1,text))

        text = [t for t in text if t[0] != '\\']

        if self.verbose: 
            for t in text: print(t)
        
        os.remove(temporaryFileName)
        estimatedTimeForTheQuery = re.findall("The estimated time for your query to complete is ([0-9]+) seconds",text[11])[0]
        httpAddress = text[12].split()[-1][1:] + text[13][:-2]

        startTime = time.time()
        timeout = 2.*max(5.0,float(estimatedTimeForTheQuery))
        refreshTime  = 3.0  
        endString = "The state of your query is 2 (Query complete)"
        regexpr = re.compile("wget (.*.fits)")

        links = None
        fakeName = "__temp__query__result.html"

        time.sleep(refreshTime)
        
        while(time.time() <= startTime+timeout):
            try:
                (filename, header) = urllib.request.urlretrieve(httpAddress,fakeName)
            except:
                urllib.request.urlcleanup()
                continue
            f = open(fakeName)
            html = " ".join(f.readlines())
            try:
                status = re.findall("The state of your query is ([0-9]+)",html)[0]
            except:
                status = '0'
                pass

            if(status=='2'):
                links = regexpr.findall(html)
                break
            f.close()
            os.remove(fakeName)
            urllib.request.urlcleanup()
        
        np.save("./{}/LAT/link.npy".format(self.full_name), links)
        for lk in links:
            print("Downloading ... ", lk)
            urllib.request.urlretrieve(lk, "./{}/LAT/{}-{}.fits".format(self.full_name,self.event_name, lk[-9:-5]))

        if self.verbose: print("Completed (LAT)")
 
    def __DownloadProcess__(self, dwUrl, dwFile):
        if self.forcedDwn or not(os.path.isfile(dwFile)):
            if self.verbose: print("Downloading... ", dwUrl)
            urllib.request.urlretrieve(dwUrl, dwFile)

