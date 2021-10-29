import gbm
from gbm import gspec
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from astropy.io import fits
from gt_apps import *

from ..info import EventInfo


class PrepLLE(EventInfo):
    
    def __init__(self, grb, tRan = [None, None], forcedGen = False, verbose=False, plotting=False, **kwargs):
 

        super().__init__(grb, tRan = tRan, **kwargs)
        self.verbose=verbose
        self.forcedGen=forcedGen

        self.PhotonEvent()

        self.Background(plotting=plotting)

    def PhotonEvent(self):
        evtbin['evfile'] = "./{}/LAT/gll_lle_bn{}_v00.fit".format(self.full_name, self.full_name)
        evtbin['scfile'] = './{}-SC00.fits'.format(self._address_LAT)
        evtbin['outfile'] = './{}-{}-{}.pha'.format(self._Address, det_spec, timeInt)
        evtbin['algorithm'] = 'PHA2'
        evtbin['tbinalg'] = 'LIN'
        evtbin['tstart'] = self.toi[0] + self.trigger
        evtbin['tstop'] = self.toi[1] + self.trigger
        evtbin['dtime'] = self.toi[1] - self.toi[0]
        evtbin['ebinalg'] = 'LOG'
        evtbin['emin'] = 10
        evtbin['emax'] = 1e4
        evtbin['enumbins'] = 50
        evtbin['chatter'] = self.verbose
        evtbin.run(print_command=self.verbose)

        with fits.open('./Xspec/{}-LLE-{}.pha'.format(self._address_Xspec, self._ti_num), mode='update') as phafile:
                phafile[0].header['FILENAME'] = '{}-LLE-{}.pha'.format(self.event_name,  self._ti_num)
                phafile[1].header['FILENAME'] = '{}-LLE-{}.pha'.format(self.event_name,  self._ti_num)
                phafile[2].header['FILENAME'] = '{}-LLE-{}.pha'.format(self.event_name,  self._ti_num)
                phafile[3].header['FILENAME'] = '{}-LLE-{}.pha'.format(self.event_name,  self._ti_num)
                phafile[2].header['RESPFILE'] = '{}-LLE-{}.rsp'.format(self.event_name,  self._ti_num)
                phafile[2].header['BACKFILE'] = '{}-LLE-{}.bak'.format(self.event_name,  self._ti_num)

    def Background(self, plotting=False):
    
        data = gspec.GspecManager()

        forder, maxt = self.BackgroundFit(det, plotting=plotting)

        detName = data.add_data("./{}/LAT/glg_cspec_ll_bn{}_v00.pha".format(self.full_name, self.full_name))
        data.add_response(detName,"./{}/LAT/glg_cspec_ll_bn{}_v00.rsp".format(self.full_name, det, self.full_name))
        ata.energy_selection(detName, data.cspec_lle_default)
        data.clear_time_binning(detName)
        
        data.source_selection(detName, [self.toi[0], self.toi[1]])
        data.add_background(detName, [[-300, -50], [self._tl100, self._tl100+200]],'Polynomial', forder)
        data.export_to_xspec(detName, directory='./{}/Xspec/'.format(self.full_name))
        
        os.system('rm ./{}/Xspec/glg_cspec_ll_bn{}_xspec_v00.pha'.format(self.full_name, self.full_name))
        os.system('mv ./{}/Xspec/glg_cspec_ll_bn{}_xspec_v00.bak ./{}-LLE-{}.bak'.format(self.full_name, self.full_name, self._address_Xspec, self._ti_num))
        os.system('cp ./{}/LAT/glg_cspec_ll_bn{}_v00.rsp ./{}-LLE-{}.rsp'.format(self.full_name, self.full_name, self._address_Xspec, self._ti_num))
        
        if self.verbose:
            print("A background file is created: ./{}-LLE-{}.bak".format(self._address_Xspec, self._ti_num))
            print("A response file is copied: ./{}-LLE-{}.rsp".format(self._address_Xspec, self._ti_num))

        with fits.open('./{}-LLE-{}.bak'.format(self._address_Xspec, self._ti_num), mode='update') as backfile:
            backfile[0].header['FILENAME'] = '{}-LLE-{}.bak'.format(self.event_name, self._ti_num)
            backfile[1].header['FILENAME'] = '{}-LLE-{}.bak'.format(self.event_name, self._ti_num)
            backfile[2].header['CHANTYPE'] = 'PI'
            backfile[2].header['FILENAME'] = '{}-LLE-{}.bak'.format(self.event_name, self._ti_num)
            backfile[3].header['FILENAME'] = '{}-LLE-{}.bak'.format(self.event_name, self._ti_num)
        
    def Background(self, plotting=False):
    
        data = gspec.GspecManager()

        forder, maxt = self.BackgroundFit(det, plotting=plotting)

        detName = data.add_data("./{}/LAT/glg_cspec_ll_bn{}_v00.pha".format(self.full_name, self.full_name))
        data.add_response(detName,"./{}/LAT/glg_cspec_ll_bn{}_v00.rsp".format(self.full_name, det, self.full_name))
        ata.energy_selection(detName, data.cspec_lle_default)
        data.clear_time_binning(detName)
        
        data.source_selection(detName, [self.toi[0], self.toi[1]])
        data.add_background(detName, [[-300, -50], [self._tl100, self._tl100+200]],'Polynomial', forder)
        data.export_to_xspec(detName, directory='./{}/Xspec/'.format(self.full_name))
        
        os.system('rm ./Xspec/{}_xspec_v00.pha'.format(detName[:-8]))
        os.system('mv ./Xspec/{}_xspec_v00.bak ./Xspec/{}-{}-{}.bak'.format(detName[:-8], GRBname, det, ti))
        os.system('rm ./Xspec/{}_xspec_v00.rsp'.format(detName[:-8]))
        if self.verbose:
            print("A photon-event file is created: ./{}-{}-{}.pha".format(self._address_Xspec, det, self._ti_num))
            print("A background file is created: ./{}-{}-{}.bak".format(self._address_Xspec, det, self._ti_num))
            print("A response file is created: ./{}-{}-{}.rsp".format(self._address_Xspec, det, self._ti_num))

        with fits.open('./Xspec/{}-{}-{}.bak'.format(GRBname, det, ti), mode='update') as backfile:
            backfile[0].header['FILENAME'] = '{}-{}-{}.bak'.format(GRBname, det, ti)
            backfile[1].header['FILENAME'] = '{}-{}-{}.bak'.format(GRBname, det, ti)
            backfile[2].header['CHANTYPE'] = 'PI'
            backfile[2].header['FILENAME'] = '{}-{}-{}.bak'.format(GRBname, det, ti)
            backfile[3].header['FILENAME'] = '{}-{}-{}.bak'.format(GRBname, det, ti)


    def BackgroundFit(self, plotting=False):

        for order in [1,2,3,4]:
            chidist = []
            ytot = 0
            mytot = 0

            events = fits.open("./{}/LAT/gll_lle_bn{}_v00.fit".format(self.full_name, self.full_name))[1].data.field("Time")
            events = events - self.trigger
            #plt.hist(events, np.arange(min(events), max(events), step=1), histtype='step')
            y1, x1, etc = histogram(events, np.arange(max(min(events)-50, self.trigger-200), self.trigger-50, step=1))
            y2, x2, etc = histogram(events, np.arange(self._tl100, min(self._tl100+200, max(events)), step=1))
            x1 = ((x1[1:]+x1[:-1])/2.).tolist()
            x2 = ((x2[1:]+x2[:-1])/2.).tolist()
            xbak = x1+x2
            ybak = y1.tolist()+y2.tolist()
            z = np.polyfit(x,y,order)
            my = np.poly1d(z)
            mytot += my(x)
            chi = sum((ybak-my(xbak))**2./my(xbak))
            chi_nu = chi/(len(ybak)-order-1)
            chidist.append(chi_nu)
            
            if max(chidist) < 1.2:
                forder = order
                if self.verbose: 
                    print("Background polynomial fit for {} (order = {}): Reduced chi^2 = {:.2f}".format(det, forder, np.average(chidist)))
                break
            else:
                forder = 1
        
        if plotting or self.verbose > 1:
            f, ax = plt.subplots(1,2, figsize=(12, 5))
            ax[0].plot(x, ytot)
            ax[0].plot(x, mytot)
            ax[0].axvline(0, ls=":", color='r')
            ax[0].axvspan(self.toi[0], self.toi[1], facecolor='orange', edgecolor="k", alpha=0.2, lw=1, zorder=-1)
            ax[0].set_xlim(-200, max(xbak))
            ax[0].set_ylim(-10,)
            ax[0].set_title("GRB{}".format(self.full_name))
            ax[0].set_xlabel(r"Time Since Trigger ($T_0$)")
            ax[0].set_ylabel("Counts")

            ax[1].plot(range(len(chidist)), chidist)
            ax[1].set_xlabel("Energy channel")
            ax[1].set_ylabel(r"$\chi^2$")
            ax[1].axhline(1, ls=':', color='gray')
            ax[1].set_title("lle, {}-order polynomial fit".format(forder))
            plt.savefig("./{}/LAT/lle_bk_chisq.png".format(self.full_name))
            plt.show(block=False)
        
        
        return forder, max(xbak)
