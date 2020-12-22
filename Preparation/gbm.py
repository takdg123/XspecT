import gbm
from gbm import gspec
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from astropy.io import fits

from XspecT.info import EventInfo


class PrepGBM(EventInfo):
    
    def __init__(self, grb, tRan = [None, None], forcedGen = False, verbose=False, version=0, **kwargs):
        
        if hasattr(grb, 'dtype'):
            grb = grb.decode('UTF-8')
        else:
            try:
                grb = grb.astype(str)
            except:
                pass

        super().__init__(grb, tRan = tRan, version = version, **kwargs)
        self.verbose=verbose
        self.forcedGen=forcedGen

        self.Preparation()
        
    def Preparation(self):
        for det in self.usedGBM:

            bkFileFlag = os.path.isfile('./{}/Xspec/{}-{}-.bak'.format(self.fullname, self.GRB_name, det))
            if self.forcedGen:
                if bkFileFlag:
                    os.system("rm ./{}/Xspec/glg_cspec_{}_bn{}_xspec_v00.bak".format(self.fullname, det, self.fullname))
            elif bkFileFlag:
                return

            data = gspec.GspecManager()

            forder, maxt = self.polyFit(det)

            detName = data.add_data("./{}/GBM/glg_cspec_{}_bn{}_v00.pha".format(self.fullname, det, self.fullname))
            try:
                data.add_response(detName,"./{}/GBM/glg_cspec_{}_bn{}_v00.rsp2".format(self.fullname, det, self.fullname))
            except:
                data.add_response(detName,"./{}/GBM/glg_cspec_{}_bn{}_v00.rsp".format(self.fullname, det, self.fullname))
            if 'n' in det:
                data.energy_selection(detName, data.cspec_nai_default)
            elif 'b' in det:
                data.energy_selection(detName, data.cspec_bgo_default)
            else:
                continue
                
            data.clear_time_binning(detName)
            

            data.source_selection(detName, [self._tStart, self._tEnd])
            order = 1
            data.add_background(detName, [[-200, -50], [self._tEnd+50, self._tEnd+maxt]],'Polynomial', forder)
            data.export_to_xspec(detName, directory='./{}/Xspec/'.format(self.fullname))
            
            os.system('mv ./{}/Xspec/glg_cspec_{}_bn{}_xspec_v00.bak ./{}/Xspec/{}-{}.bak'.format(self.fullname, det, self.fullname, self.fullname, self.GRB_name, det))
            os.system('mv ./{}/Xspec/glg_cspec_{}_bn{}_xspec_v00.pha ./{}/Xspec/{}-{}.pha'.format(self.fullname, det, self.fullname, self.fullname, self.GRB_name, det))
            os.system('mv ./{}/Xspec/glg_cspec_{}_bn{}_xspec_v00.rsp ./{}/Xspec/{}-{}.rsp'.format(self.fullname, det, self.fullname, self.fullname, self.GRB_name, det))
            
            with fits.open('./{}/Xspec/{}-{}.pha'.format(self.fullname, self.GRB_name, det), mode='update') as phafile:
                phafile[0].header['FILENAME'] = '{}-{}.pha'.format(self.GRB_name, det)
                phafile[1].header['FILENAME'] = '{}-{}.pha'.format(self.GRB_name, det)
                phafile[2].header['FILENAME'] = '{}-{}.pha'.format(self.GRB_name, det)
                phafile[3].header['FILENAME'] = '{}-{}.pha'.format(self.GRB_name, det)
                phafile[2].header['RESPFILE'] = './{}/Xspec/{}-{}.rsp'.format(self.fullname, self.GRB_name, det)
                phafile[2].header['BACKFILE'] = './{}/Xspec/{}-{}.bak'.format(self.fullname, self.GRB_name, det)

            with fits.open('./{}/Xspec/{}-{}.bak'.format(self.fullname, self.GRB_name, det), mode='update') as backfile:
                backfile[0].header['FILENAME'] = '{}-{}.bak'.format(self.GRB_name, det)
                backfile[1].header['FILENAME'] = '{}-{}.bak'.format(self.GRB_name, det)
                backfile[2].header['FILENAME'] = '{}-{}.bak'.format(self.GRB_name, det)
                backfile[3].header['FILENAME'] = '{}-{}.bak'.format(self.GRB_name, det)
                
            with fits.open('./{}/Xspec/{}-{}.rsp'.format(self.fullname, self.GRB_name, det), mode='update') as rspfile:
                rspfile[0].header['FILENAME'] = '{}-{}.rsp'.format(self.GRB_name, det)
                rspfile[1].header['FILENAME'] = '{}-{}.rsp'.format(self.GRB_name, det)
                rspfile[2].header['FILENAME'] = '{}-{}.rsp'.format(self.GRB_name, det)

    def polyFit(self, det):
        for order in [1,2,3,4]:
            chidist = []
            ytot = 0
            mytot = 0
            for i in range(16, 120):
                
                y = fits.open("./{}/GBM/glg_cspec_{}_bn{}_v00.pha".format(self.fullname, det, self.fullname))[2].data['Counts'][:,i]/fits.open("./{}/GBM/glg_cspec_{}_bn{}_v00.pha".format(self.fullname, det, self.fullname))[2].data['Exposure']
                x = fits.open("./{}/GBM/glg_cspec_{}_bn{}_v00.pha".format(self.fullname, det, self.fullname))[2].data['Time']-self.Trigger
                ytot +=y
                bakRan = ((x>-200) * (x<-50))  + ((x>(self._tEnd+50))  * (x<(self._tEnd+300)))
                    
                xbak = x[bakRan]
                ybak = y[bakRan]
                
                if i==16: 
                    xlim = max(xbak[ybak!=0])-10
                
                ybak = ybak[xbak<xlim]
                xbak = xbak[xbak<xlim]
                
                z = np.polyfit(xbak,ybak,order)
                my = np.poly1d(z)
                mytot += my(x)
                
                chi = sum((ybak-my(xbak))**2./my(xbak))
                chi_nu = chi/(len(ybak)-order-1)
                chidist.append(chi_nu)
            
            if max(chidist) < 1.2:
                forder = order
                if self.verbose: 

                    print("Background polynomial fit for {} (order = {}): Reduced chi^2 = {:.2f}".format(forder, det, np.average(chidist), det, forder))
                break
            else:
                forder = 1
        
        if self.verbose:
            f, ax = plt.subplots(1,2, figsize=(12, 5))
            ax[0].plot(x, ytot)
            ax[0].plot(x, mytot)
            ax[0].axvline(0, ls=":", color='r')
            ax[0].axvspan(self._tStart, self._tEnd, facecolor='orange', edgecolor="k", alpha=0.2, lw=1, zorder=-1)
            ax[0].set_xlim(-200, max(xbak))
            ax[0].set_ylim(-10,)
            ax[0].set_title("GRB{}".format(self.fullname))
            ax[0].set_xlabel(r"Time Since Trigger ($T_0$)")
            ax[0].set_ylabel("Counts")

            ax[1].plot(range(len(chidist)), chidist)
            ax[1].set_xlabel("Energy channel")
            ax[1].set_ylabel(r"$\chi^2$")
            ax[1].axhline(1, ls=':', color='gray')
            ax[1].set_title("{}, {}-order polynomial fit".format(det, forder))
            plt.savefig("./{}/GBM/{}_bk_chisq.png".format(self.fullname, det))
            plt.show(block=False)
        
        
        return forder, max(xbak)
