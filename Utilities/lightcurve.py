import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

from .background import BackgroundFit
from . import utils 
from ..LATanalysis import UnbinnedLAT
from ..external.LinearFit import *

from astropy.table import Table

def GBMlcData(info, eCut=[50, 300]):
    rawData = []
    for i in range(np.size(eCut)-1):
        rawData.append([])

    for dtr, i in zip(info.usedGBM, range(np.size(info.usedGBM))):

        with fits.open("./{}/GBM/glg_tte_{}_bn{}_v00.fit".format(info.full_name,dtr,info.full_name)) as sourceFile:

            if 'n' in dtr:
                cut = (sourceFile[1].data.field('E_MIN')>10)*(sourceFile[1].data.field('E_MAX')<1000)
            elif 'b' in dtr:
                cut = (sourceFile[1].data.field('E_MIN')>200)*(sourceFile[1].data.field('E_MAX')<40000)

            ecut_temp=[]    
            for elim in eCut:
                ecut_temp.append(128-sum(sourceFile[1].data.field('E_MIN')>elim))

            dataT = sourceFile[2].data.field('TIME')-info.trigger
            dataChannel = sourceFile[2].data.field('PHA')

            for j in range(np.size(ecut_temp)-1): 
                lcData_temp=np.asarray(dataT)*np.asarray(dataChannel>=ecut_temp[j])*np.asarray(dataChannel<=ecut_temp[j+1])*np.asarray(dataT>=-300)*np.asarray(dataT<=300)
                lcData_temp = lcData_temp[lcData_temp != 0]
                rawData[j]=rawData[j]+lcData_temp.tolist()

    return np.asarray(rawData)

def LATlcData(info, prob=0.9, emin=100, target = "GRB", return_arrray=False):
    if os.path.isfile('./{}-srcprob.fits'.format(info._address_LAT)):
        rawData = []    
        with fits.open('./{}-srcprob.fits'.format(info._address_LAT)) as sourceFile:

            for t, e, p in zip(sourceFile[1].data["TIME"], sourceFile[1].data["ENERGY"], sourceFile[1].data[target]):
                if p >= prob and e>=emin:
                    rawData.append([t - info.trigger, e, p])
            rawData.sort()
        rawData = np.asarray(rawData)
        if return_arrray:
            return rawData
        else:
            tab = Table(rawData, names=["Time", "Energy", "Probability"])
            return tab
    else:
        return []

def LightCurve(info, ext_ax=None, sed_fig=False, return_ax=True, toi = [None, None], eRan = [50, 300], tb=0.064, ts = 0.064, ti_flag = False, savefig=False):

    rawData = GBMlcData(info, eCut=eRan)[0]

    bk, err = BackgroundFit(rawData, t1=-10, t2=info.t95+50)
    lcEdge = np.arange(-info.t90, info.t90*2, step=ts)
    cnts, t = np.histogram(rawData, bins=lcEdge)
    ts = np.asarray(t[1:]-t[:-1])
    t = utils.center_pt(t)
    B = np.asarray(bk(t))
    
    if ext_ax == None:
        f, ax = plt.subplots(1, 1, figsize=(8, 3))
    else:
        ax = ext_ax

    ax.step(lcEdge, [((cnts/ts-B/tb)/1000.)[0]]+((cnts/ts-B/tb)/1000.).tolist(), color='k', label="GBM Events")
    ax.set_xlabel("Time since triggered [s]", fontsize=12)
    ax.set_ylabel(r"Counts/s", fontsize=12)
    ax.text(-0.03, 1.03, r"$\times$ 10$^3$", transform=ax.transAxes)
    ax.text(1, 1.03, "{:.0f} - {:.0f} keV".format(eRan[0], eRan[1]), transform=ax.transAxes, ha="right")
    ax.axhline(0, ls=":", color='k')
    ax.axvline(0, lw=0.5, color='r')
    ax.set_xlim(min(-info.t90*0.1, info.t05-info.t90*0.1), max(info.t90*1.2, info.t95+info.t90*0.2))

    ax.set_title(info.full_name)
    ax.legend(loc=2)
    ax.set_ylim(-0.2, )
    rawData_LAT_9 = LATlcData(info, return_arrray=True)
    rawData_LAT = LATlcData(info, prob=0, return_arrray=True)

    if len(rawData_LAT) !=0:
        ax2=ax.twinx()
        ax2.scatter(rawData_LAT[:,0], rawData_LAT[:,1], color="w", edgecolor="k")
        if len(rawData_LAT_9) != 0:
            ax2.scatter(rawData_LAT_9[:,0], rawData_LAT_9[:,1], color='r', label="LAT Events")  
        ax2.set_ylabel(r"Energy [MeV]", fontsize=12)
        ax2.set_yscale("log")
        if max(rawData_LAT_9[:,1]) >= 1e4:
            ax2.set_ylim(45, 1e5)
        else:
            ax2.set_ylim(60, 1e4)

        ax2.legend()

    if sed_fig:
        ax.axvspan(toi[0],toi[1], facecolor='orange', edgecolor="k", alpha=0.5, lw=1, zorder=-1)
    else:
        if ti_flag:
            h=0
            for toi in info.time_intervals:
                ax.axvspan(toi[0], toi[1], ymin = 0.9-0.1*(h+1), ymax = 0.9-0.1*h, facecolor='orange', edgecolor="k", alpha=0.2, lw=1, zorder=-1)   		
                h+=1
        else:
            if toi[0] == None and toi[1] == None:
                ax.axvspan(info.t05,info.t95, facecolor='orange', edgecolor="k", alpha=0.2, lw=1, zorder=-1)
            else:
                ax.axvspan(toi[0],toi[1], facecolor='orange', edgecolor="k", alpha=0.2, lw=1, zorder=-1)

    if savefig: plt.savefig("./{}/GBM/lightcurve.png".format(info.full_name))
    
    if ext_ax==None and return_ax==False:
        plt.show(block=False)
    else:
        return ax

def LightcurveLAT(tab, tran=[0.1, 1000], fran=[1e-9, None]):
    
    f, ax = plt.subplots(2,1, figsize=(9,9), gridspec_kw = {'height_ratios':[7, 2]})

    lfit = Linear_fit(tab["Tc"], tab["eFlux"], y_err = tab["eFlux_err_high"], logy=True, logx=True)
    lfit.runFit()

    ax[0].errorbar(tab["Tc"], tab["eFlux"], xerr=[tab["Tc"]-tab["From"], tab["To"]-tab["Tc"]], yerr=[tab["eFlux_err_low"], tab["eFlux_err_high"]], ls="",label="Energy flux")
    ax[0].plot(np.linspace(tab["From"][0], tab["To"][-1], 10), lfit.linear(np.linspace(tab["From"][0], tab["To"][-1], 10)), lw=2, label="Index={:.1f}+/-{:.1f}".format(lfit.p[0], lfit.perr[0]))
    
    ax[1].errorbar(tab["Tc"], tab["Index"], xerr=[tab["Tc"]-tab["From"], tab["To"]-tab["Tc"]], yerr=tab["Index_err"], ls="")
    
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
    ax[0].set_yscale("log")

    ax[0].set_xlim(tran[0], tran[1])
    ax[1].set_xlim(tran[0], tran[1])
    ax[0].set_ylim(fran[0], fran[1])

    ax[0].set_ylabel(r"Flux (0.1 - 10 GeV) [erg cm$^{-2}$ s$^{-1}$]", fontsize=12)
    ax[1].set_ylabel("Index", fontsize=12)
    ax[1].set_xlabel(r"Time since T$_0$ [s]", fontsize=12)
    ax[0].grid()
    ax[1].grid()
    ax[1].axhline(-2, ls=":", color="k")
    ax[0].legend(frameon=False, fontsize=12)

    return tab, ax