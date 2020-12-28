import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

from .background import BackgroundFit
from . import utils 

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

def LATlcData(info):
    if os.path.isfile('./{}-srcprob.fits'.format(info._address_LAT)):
        rawData = []    
    
        with fits.open('./{}-srcprob.fits'.format(info._address_LAT)) as sourceFile:

            for t, e, p in zip(sourceFile[1].data["TIME"], sourceFile[1].data["ENERGY"], sourceFile[1].data["GRB"]):
                if p >= 0.9:
                    rawData.append([t - info.trigger, e])
            rawData.sort()

        return np.asarray(rawData)
    else:
        return []

def LightCurve(info, ax=None, toi = [None, None], eRan = [50, 300], tb=0.064, ts = 0.064, savefig=False):

    rawData = GBMlcData(info, eCut=eRan)[0]

    bk, err = BackgroundFit(rawData, t1=-10, t2=info.t95+50)
    lcEdge = np.arange(-info.t90, info.t90*2, step=ts)
    cnts, t = np.histogram(rawData, bins=lcEdge)
    ts = np.asarray(t[1:]-t[:-1])
    t = utils.center_pt(t)
    B = np.asarray(bk(t))
    
    if ax == None:
        f, ax = plt.subplots(1, 1, figsize=(8, 3))

    ax.step(lcEdge, [((cnts/ts-B/tb)/1000.)[0]]+((cnts/ts-B/tb)/1000.).tolist(), color='k', label="GBM Events")
    ax.set_xlabel("Time since triggered [s]", fontsize=12)
    ax.set_ylabel(r"Counts/s", fontsize=12)
    ax.text(-0.03, 1.03, r"$\times$ 10$^3$", transform=ax.transAxes)
    ax.text(1, 1.03, "{:.0f} - {:.0f} keV".format(eRan[0], eRan[1]), transform=ax.transAxes, ha="right")
    ax.axhline(0, ls=":", color='k')
    ax.axvline(0, lw=0.5, color='r')
    ax.set_xlim(min(-info.t90*0.1, info.t05-info.t90*0.1), max(info.t90*1.2, info.t95+info.t90*0.2))

    if toi[0] == None and toi[1] == None:
        ax.axvspan(info.t05,info.t95, facecolor='orange', edgecolor="k", alpha=0.2, lw=1, zorder=-1)
    else:
        ax.axvspan(toi[0],toi[1], facecolor='orange', edgecolor="k", alpha=0.2, lw=1, zorder=-1)

    ax.set_title(info.full_name)
    ax.legend(loc=2)
    ax.set_ylim(-0.2, )
    rawData_LAT = LATlcData(info)

    if len(rawData_LAT) !=0:
        ax2=ax.twinx()
        ax2.scatter(rawData_LAT[:,0], rawData_LAT[:,1], color='r', label="LAT Events")    
        ax2.set_ylabel(r"Energy [MeV]", fontsize=12)
        ax2.set_yscale("log")
        if max(rawData_LAT[:,1]) >= 1e4:
            ax2.set_ylim(45, 1e5)
        else:
            ax2.set_ylim(60, 1e4)

        ax2.legend()

    if savefig: plt.savefig("./{}/GBM/lightcurve.png".format(info.full_name))
    
    if ax==None:
        plt.show(block=False)
    else:
        return ax
