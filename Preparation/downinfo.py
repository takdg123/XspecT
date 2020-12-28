import os
import urllib

def DownloadGBMcat(grb, verbose=True, scat=True, bcat=True):
    fullname = grb
    name = grb[:6]
    yr = grb[:2]
    address="https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts/"
    
    if bcat:
        if not(os.path.isfile("./{}/GBM/{}-bcat.fit".format(fullname,name))):
            for version in range(10):
                version=9-version
                try:
                    urllib.request.urlopen(urllib.request.Request("{}/20{}/bn{}/current/glg_bcat_all_bn{}_v0{}.fit".format(address,yr,fullname,fullname,version))) 
                    break
                except:
                    continue

            urllib.request.urlretrieve("{}/20{}/bn{}/current/glg_bcat_all_bn{}_v0{}.fit".format(address,yr,fullname,fullname,version), "./{}/GBM/{}-bcat.fit".format(fullname,name))
            if verbose: print("Downloading bcat is completed.")
    
    if scat:
        if not(os.path.isfile("./{}/GBM/{}-scat.fit".format(fullname,name))):
            for version in range(10):
                version=9-version
                try:
                    urllib.request.urlopen(urllib.request.Request("{}/20{}/bn{}/current/glg_scat_all_bn{}_flnc_plaw_v0{}.fit".format(address,yr,fullname,fullname,version))) 
                    break
                except:
                    continue

            try:
                urllib.request.urlretrieve("{}/20{}/bn{}/current/glg_scat_all_bn{}_flnc_plaw_v0{}.fit".format(address,yr,fullname,fullname,version), "./{}/GBM/{}-scat.fit".format(fullname,name))
            except:
                print("scat file does not exist.")
            if verbose: print("Downloading scat is completed.")

