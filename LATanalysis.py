from UnbinnedAnalysis import *
from astropy.table import Table
import numpy as np
import os

def AnalysisLAT(info, time_intervals, emin=100, emax=10000, verbose=True):
    tab=None

    for i in time_intervals:
        if verbose:
            print("Analysis on the time interval {} starts.".format(i), end="\r")
        tab = UnbinnedLAT(info, i, tab=tab, emin=emin, emax=emax)
        if verbose:
            print("Analysis on the time interval {} is completed.".format(i), end="\r")
    return tab

def UnbinnedLAT(info, ti_num, verbose=False, full_output=False, tab=None, emin=100, emax=10000):
    path = info.current_dir+"/"+info.full_name+'/LAT/'
    files = os.listdir(path)
    scfile = path+[f for f in files if "SC" in f][0]

    obs = UnbinnedObs('./{}-mktime-{}.fits'.format(info._address_LAT, ti_num), scfile,
                    expMap='./{}-expmap-{}.fits'.format(info._address_LAT, ti_num), 
                    expCube='./{}-ltcube-{}.fits'.format(info._address_LAT, ti_num),
                    irfs=info._irf_class)
    like = UnbinnedAnalysis(obs,'./{}-model.xml'.format(info._address_LAT),optimizer='NewMinuit')
    likeobj = pyLike.NewMinuit(like.logLike)
    like.tol = 0.1
    like.tolType = 1

    logLikelihood = like.fit(verbosity=1,covar=True,optObject=likeobj)
    like.tol = 0.001
    try:
        logLikelihood = like.fit(verbosity=0,covar=True,optObject=likeobj)
    except:
        pass
    
    errCode = likeobj.getRetCode()

    if errCode != 0:

        for source in like.sourceNames():
            TS = like.Ts(source)
            if source not in ['GRB','iso_P8R2_TRANSIENT020E_V6_v06','gll_iem_v06'] and TS < 5:
                like.deleteSource(source)

        logLikelihood = like.fit(verbosity=0,covar=True,optObject=likeobj)


    TS = like.Ts('GRB')
    
    EnergyFlux = like.energyFlux('GRB', emin=emin, emax=emax)*1.6022*10**-6
    EnergyFlux_err = like.energyFluxError('GRB', emin=emin, emax=emax)*1.6022*10**-6
    N = like.model['GRB'].funcs['Spectrum'].getParam('Integral').value()*1e-06
    N_err = like.model['GRB'].funcs['Spectrum'].getParam('Integral').error()*1e-06
    Index = like.model['GRB'].funcs['Spectrum'].getParam('Index').value()
    Index_err = like.model['GRB'].funcs['Spectrum'].getParam('Index').error()
    
    eFlux = lambda E1, E2, N0, gamma: N0*(gamma+1)/(gamma+2)*(E1**(gamma+2)-E2**(gamma+2))/(10000**(gamma+1)-100**(gamma+1))

    isamp = np.random.normal(Index,Index_err,10000) 
    fsamp = np.random.normal(N,N_err,10000)

    eflx = eFlux(emax, emin, fsamp, isamp)
    eflx_alt = np.asarray([np.percentile(eflx,50), np.percentile(eflx,50)-np.percentile(eflx,16), np.percentile(eflx,84)-np.percentile(eflx,50)])
    eflx_alt = eflx_alt*1.6022*10**-6

    properties=[[info.time_intervals[ti_num][0]],
                [info.time_intervals[ti_num][1]],
                [10**((np.log10(info.time_intervals[ti_num][0])+np.log10(info.time_intervals[ti_num][1]))/2.)],
                [TS], 
                [eflx_alt[0]], [eflx_alt[1]], [eflx_alt[2]], 
                [N], [N_err],
                [Index], [Index_err]]
    
    if tab == None:
        prop_table = Table(properties, names=["From", "To", "Tc", "TS", "eFlux", "eFlux_err_low","eFlux_err_high", "N", "N_err", "Index", "Index_err"])
    else:
        prop_table = tab
        prop_table.add_row(properties)
        
    if errCode != 0:
        print('NOT CONVERGE({}) '.format(likeobj.getRetCode()))
    else:
        if full_output:
            return prop_table, like, likeobj, errCode
        else:
            return prop_table
