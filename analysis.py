import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import glob

from astropy.io import fits

import gt_apps
from gt_apps import *
from GtApp import GtApp
bindef = GtApp('gtbindef','Likelihood')

from scipy.integrate import quad
from datetime import datetime

from xspec import *
from . import Utilities
from .Utilities.functions import *
from .info import EventInfo

class Analysis(EventInfo):

    def __init__(self, info=None, EvtInfo={}, verbose = False, only='GBMLAT', **kwargs):
        self.verbose = verbose
        if self.verbose == 2 or info==None:
            print("config = {}")
            print("config['Name']: e.g., '131108A'")
            print("config['Trigger']: e.g., 405636118")
            print("config['RA']: e.g., 156.50187")
            print("config['DEC']: e.g., 9.66225")
            print("config['T90']: [T05, T95]; e.g., [0.319986, 18.496286]")
            print("config['usedGBM']: ['n0','n3','n6','n7','b0','b1']")
            print("config['TimeEdge']: Time edges; e.g., [0, 1, 2, 3, 4, 5]")
            print("config['TimeInterval']: Time intervals; e.g., [[0, 1], [1, 2], [2, 3]]")
            print("config['TimeStart']: Start times; e.g., [0, 1, 2, 3, 4] ")
            print("config['TimeEnd']: End times; e.g., [1, 2, 3, 4, 5] ")
            print("Fitting = Analysis(grb, EvtInfo=config) or Fitting.Configuration(config)")
            print("Fitting.RunFit(1,1)")
            print("Fitting.PrintSED()")
            print("Fitting.PrintResult()")
            sys.exit()

        if type(info) == str or type(info) == int:
            super().__init__(str(info),  **kwargs)
            self.Configuration(EvtInfo, only=only)
        elif type(info) == dict:
            info = {**info, **EvtInfo}
            self.Configuration(info, only=only)

    def Configuration(self, config, only='GBMLAT'):
        keys = config.keys()
        if "Name" in keys:          self.event_name = config['Name']
        if "Trigger" in keys:       self.trigger = config['Trigger']
        if "T90" in keys:           
            self.t05 = config['T90'][0]
            self.t95 = config['T90'][1]
            self.t90 = config['T90'][1] - config['T90'][0]
        if "RA" in keys:            self.refRa = config['RA']
        if "DEC" in keys:           self.refDec = config['DEC']
        if "usedGBM" in keys:       self.usedGBM = config['usedGBM']
    
        if 'redshift' in keys:      self.redshift = config['redshift']
        if 'nH' in keys:            self._nH = config['nH']
        
        self._only = only
        self._usedData = self.__SelectedData__(self._only)        
        self.__Reset__()
        if self.verbose:
            self.PrintTimeItv()


    def SetInit(self, m_num, init = {}, verbose=True):
        try:
            AllModels.addPyMod(bkn2power, bkn2powParInfo, 'add')
            #AllModels.addPyMod(doublegrbm, doublegrbmInfo, 'add')
            #AllModels.addPyMod(sbpl_gbm, sbpl_gbm_Info, 'add')
            #AllModels.addPyMod(sbpl_smp, sbpl_smp_Info, 'add')
            #AllModels.addPyMod(sbpl_granot, sbpl_granot_Info, 'add')
            #AllModels.addPyMod(PowerLaw_piv, pl_pivInfo, 'add')
            AllModels.addPyMod(multiBB, mBB_Info, 'add')
        except:
            pass
        models = self.__RearrModel__(m_num)
        Model_Xspec = Model(models)
        parNum = 1
        for comp in Model_Xspec.componentNames:
            parameters = Model_Xspec.__getattribute__(comp)
            if len(init)==0: 
                print(comp, parameters.parameterNames)
                print(self._setP[m_num])
            else:
                for par in parameters.parameterNames:
                    try:
                        if verbose: print(comp, par, init[parNum])
                    except:
                        if verbose: print(comp, par)
                    parNum+=1

        changedInitP = np.asarray(self._newP)
        if len(self._newP)>0:
            if m_num not in changedInitP[:,0] and init!={}:
                self._newP.append([m_num, init])
            else:
                changedInitP[np.where(changedInitP[:,0] == m_num)[0][0]][1] = init
                self._newP = changedInitP.tolist()
        else:
            self._newP.append([m_num, init])

    @classmethod
    def PrintModels(self):
        makeMtable = []
        i = 0
        j = 0
        temp = []
        for mod in Utilities.models.model_list:
            temp.append(j)
            temp.append(mod)
            i+=1
            if i >4:
                makeMtable.append(temp)
                temp = []
                i=0
            j+=1
        makeMtable.append(temp)

        return Utilities.table.ModelTable(makeMtable)

    def PrintTimeItv(self):
        makeTtable = []
        i = 0
        j = 0
        temp = []
        for toi in self.time_intervals:
            temp.append(j)
            temp.append("{:.2f} - {:.2f}".format(toi[0], toi[1]))
            i+=1
            if i >4:
                makeTtable.append(temp)
                temp = []
                i=0
            j+=1
        makeTtable.append(temp)
        return Utilities.table.IntvTable(makeTtable)

    def SaveFit(self, filename):
        np.save(filename, self.output)
        np.save(filename+"_plot", self.data_points)

    def LoadFit(self, filename):
        self.output = np.load(filename+".npy")
        self.data_points = np.load(filename+"_plot.npy")

    def RunXspec(self, m_num = 0, ti_num = 0, only=None, newInit = False, paired = False, dataType = 1, rspType = True, BAT = False, XRT = False, debug = False, **kwargs):
        if self.verbose == 2:
            print("m_num: check PrintModels")
            print("ti_num: check PrintTimeItv")
            print("paired (boolen): If True, set [[m_num, ti_num], ...] = [[m, t], ...].")
            print("newInit: update the initial parameters")
            print("rspType (boolen): use rsp2 file.")
            print("dataType (boolen)")
            sys.exit()

        self.__dict__.update(kwargs)
        
        if only !=None:
            self._only = only

        if paired:
            if np.size(m_num)>2:
                self._m_set = np.asarray(m_num)[:,0]
                self._ti_set = np.asarray(m_num)[:,1]
            else:
                self._m_set = [m_num[0]]
                self._ti_set = [m_num[1]]

        elif np.size(ti_num) == np.size(m_num):
            if np.size(m_num) == 1:
                self._ti_set = [ti_num]
                self._m_set = [m_num]
            else:
                self._ti_set = ti_num
                self._m_set = m_num

        elif ti_num == -1:
            self._ti_set = (np.zeros(np.size(m_num)) -1).astype(int)
            self._m_set = m_num

        elif np.size(ti_num) == 1 and np.size(m_num) != 1:
            self._ti_set = (np.zeros(np.size(m_num)) + ti_num).astype(int)
            self._m_set = m_num

        elif np.size(ti_num) >= 1 and np.size(m_num) >= 1:
            self._ti_set = []
            self._m_set = []
            for t in ti_num:
                for m in m_num:
                    self._ti_set.append(t)
                    self._m_set.append(m)
        else:
            sys.exit()
        
        self.output = []
        self.data_points = []
        if self.verbose: 
            print("="*75)
            print("{:^25} | {:^15} | {:^27}".format("Model","Time Intv","Status"))
            print("-"*75)

        for m_num, ti_num in zip(self._m_set, self._ti_set):
            self.toi = self.time_intervals[ti_num]
            self.model = self._model_names[m_num]
            self._m_num = m_num
            self._ti_num = ti_num
            if self.verbose and ti_num!=-1: 
                print("{:^25} | {:>6.2f} - {:<6.2f} | {:^27}".format(self.model, self.toi[0], self.toi[1], 'Working'), end='\r')
            output, data_points= self.__Fitting__(newInit, debug=debug, BAT=BAT, XRT=XRT, dataType=dataType, rspType=rspType)
            self.output.append(output)
            self.data_points.append(data_points)
            
            if self.verbose: 
                print("{:^25} | {:>6.2f} - {:<6.2f} | {:^27}".format(self.model, self.toi[0], self.toi[1], 'Done'))
        self._current_output = self.output[-1]
        self._current_data_points = self.data_points[-1]
        
        if self.verbose: print("="*75)

    def PrintResult(self, output_num = -1, clear=False, latex=False, norm=True, Epiv=100):
        if np.size(self.output) == 1:
            self._current_output = self.output[0]
            self.__MakeTable__(clear=clear, norm=norm, latex=latex, Epiv=Epiv)
            
        elif output_num!=-1:
            self._current_output = self.output[output_num]
            self._m_num = self._m_set[output_num]
            self._ti_num =  self._ti_set[output_num]
            self.__MakeTable__(clear=clear, norm=norm, latex=latex, Epiv=Epiv)
            
        else:
            self.__ClearTable__()
            for output, ti_num in zip(self.output, self._ti_set):
                self._current_output = output
                self.__MakeTable__(clear=False, norm=norm, latex=latex, Epiv=Epiv)
   
        if latex:
            FT = Utilities.table.FitTable(self.__TableList__(latex=True), parlist=self._FparN, LaTex=latex)
            print(FT._repr_latex_())
        else:
            return Utilities.table.FitTable(self.__TableList__(), parlist=self._FparN)

    
    def PrintSED(self, output_num = -1, lc_flag = False, save_flag = False, overlap=True, time_range = [], energy_range = [50, 300], verbose=False): 
        if np.size(self.output) == 1:
            self._current_output = self.output[0]
            self.__MakeSED__(lc_flag = lc_flag, save_flag = save_flag, time_range=time_range, energy_range=energy_range, overlap=overlap, verbose=verbose)
        elif output_num!=-1:
            self._current_output = self.output[output_num]
            self._m_num = self._m_set[output_num]
            self._ti_num =  self._ti_set[output_num]
            self.__MakeSED__(lc_flag = lc_flag, save_flag = save_flag, time_range=time_range, energy_range=energy_range, overlap=overlap, verbose=verbose)
        else:
            for output, m_num, ti_num in zip(self.output, self._m_set, self._ti_set):
                self._current_output = output
                self._ti_num = ti_num
                self._m_num = m_num
                self.__MakeSED__(lc_flag = lc_flag, save_flag = save_flag, time_range=time_range, energy_range=energy_range, overlap=overlap, verbose=verbose)
                '''
    def printSED_tot(self, save_flag = False, verbose=False):
        fig, ax = plt.subplots(len(self.output)+1,1,figsize=(10, 7), gridspec_kw = {'height_ratios':[7]+(np.zeros(len(self.output))+1.5).tolist()}) 
        plt.subplots_adjust(hspace=0.01)

        for i in range(len(self.output)):

            self._current_output = self.output[i]
            emax = self.__MaxEnergy__(i)
            E = np.logspace(1, emax, 100)
            Label = r"Episode {} ({:.1f}s - {:.1f}s): {}".format(i+1, self._time_start_set[i]-self.trigger, self._time_end_set[i]-self.trigger, self._model_names[self._m_set[i]])
            
            modelF, model_compF, modelF_lu=self.__ModelFlux__(E, butterfly=True, erg=True, verbose=verbose)
            PLT, = ax[0].loglog(E, E**2*modelF, label=Label)
            ax[0].fill_between(E, E**2*modelF_lu[:,0], E**2*modelF_lu[:,1], facecolor='none', alpha=0.3, edgecolor=PLT.get_c(), linewidth=1, antialiased=True)

            x = self._current_output["PlotX"]
            y = self._current_output["PlotY"]
            newArr = x[:,0].argsort()
            y = y[newArr]
            x = x[newArr]

            modelF, model_compF, etc = self.__ModelFlux__(x[:,0], verbose=verbose)
            try:
                ratio = np.sign(y[:,0]-x[:,0]**2.*modelF)*abs(np.asarray(y[:,0]-x[:,0]**2.*modelF)/abs(y[:,1]))
            except:
                ratio = 0
            ratio_err = np.zeros(len(x[:,0]))+1

            ax[i+1].errorbar(x[:,0], ratio, yerr=ratio_err, xerr=x[:,1], fmt='', ls='', c=PLT.get_c(), elinewidth=1)
            ax[i+1].text(1e6, 0,"Episode {}".format(i+1), fontsize=12)
            ax[i+1].axhline(0, xmin=0, xmax=1, color='k', lw=0.5, ls=':')
            ax[i+1].set_xlim(5, 5e6)
            ax[i+1].set_ylim(-2.5, 2.5)
            ax[i+1].set_xscale("log")
            ax[i+1].xaxis.grid()
            ax[i+1].yaxis.set_label_coords(-0.055, 0.5)  
            ax[i+1].set_yticks([-1.5, 0, 1.5])
            ax[i].set_xticklabels([])

        ax[0].set_ylabel(r'$\nu$  F$_\nu$ [erg s$^{-1}$ cm$^{-2}$]', fontsize=12)
        ax[0].legend(fontsize=12)
        ax[0].set_xlim(5, 5e6)
        ax[0].set_ylim(1e-8, 1e-5)
        ax[0].tick_params(labelsize=12)
        ax[0].grid()
        ax[-1].set_ylabel("Residual [$\sigma$]", fontsize=12)
        ax[(len(self.output)+1)/2].yaxis.set_label_coords(-0.07, 0.5)
        ax[-1].set_xlabel('Energy [keV]', fontsize=12)
        ax[-1].tick_params(labelsize=12, axis='x')

        if save_flag:
            plt.save_flag("SED.pdf")
        plt.show(block=False)
        '''

    def CalculateFlux(self, energy_range = [10, 1e7], runs = 100000, erg=False, filename='', output=False, phoflx=False):
        TotalF = []
        CompF = []
        ergCov = 1.6022*10**-9

        print("="*75)
        print("{:^25} | {:^15} | {:^27}".format("Model","Time Intv","Status"))
        print("-"*75)
        for i in range(np.size(self.output)):
            self._current_output = self.output[i]
            parSet = self.__GeneratePars__(runs = runs, sim=True)
            synF = []
            synF_comp = []
            for parS, k in zip(parSet, range(runs)):
                parNum = 0
                totF = 0
                cF = []
                for fresult in self._current_output['Model']:

                    model_name = fresult['name']
                    if model_name =='TBabs' or model_name == 'zTBabs': continue
                    
                    if phoflx:
                        func = globals()[model_name.upper()]
                    else:
                        func = globals()["e"+model_name.upper()]
                    thisParNum = func.__code__.co_argcount-1
                    flux = quad(func, energy_range[0], energy_range[1], args=tuple(parS[np.arange(thisParNum)+parNum]))[0]
                    cF.append(flux)
                    totF+=flux
                    parNum += thisParNum
                synF.append(totF)
                synF_comp.append(cF)
                k = int((k+1)*20./runs)
                print("{:^25} | {:>6.2f} - {:<6.2f} | [{:<20}] {:.0f}%".format(self.model, self.toi[0], self.toi[1], '='*k, 5*k), end='\r')
            
            print("{:^25} | {:>6.2f} - {:<6.2f} | {:^27}".format(self.model, self.toi[0], self.toi[1], 'Done'))
            print("="*75)

            if erg:
                TotalF.append([synF[-1]*ergCov, np.percentile(synF,16)*ergCov, np.percentile(synF, 84)*ergCov, self.toi[0], self.toi[1]])
            else:
                TotalF.append([synF[-1], np.percentile(synF,16), np.percentile(synF, 84), self.toi[0], self.toi[1]])
            
            CompF_temp = []

            if np.size(synF_comp, axis=1) > 1:
                synF_comp = np.asarray(synF_comp)
                CompF_temp = []
                for j in range(np.size(synF_comp, axis=1)):
                    if erg:
                        CompF_temp.append([synF_comp[:,j][-1]*ergCov, np.percentile(synF_comp[:,j],16)*ergCov, np.percentile(synF_comp[:,j], 84)*ergCov, self.toi[0], self.toi[1]])
                    else:
                        CompF_temp.append([synF_comp[:,j][-1], np.percentile(synF_comp[:,j],16), np.percentile(synF_comp[:,j], 84)*ergCov, self.toi[0], self.toi[1]])
                CompF.append(CompF_temp)
            else:
                if erg:
                    CompF.append([synF[-1]*ergCov, np.percentile(synF,16)*ergCov, np.percentile(synF, 84)*ergCov, self.toi[0], self.toi[1]])
                else:
                    CompF.append([synF[-1], np.percentile(synF,16), np.percentile(synF, 84), self.toi[0], self.toi[1]])
        
        if filename == '':
            now = datetime.now()
            filename = str(now.month)+'-'+str(now.day)
            extra = 0
            while True:
                filename = filename+'-'+str(extra)
                if os.path.isfile(filename+'_tot.npy'):
                    extra+=1
                else:
                    break

        np.save(filename+"_tot", TotalF)
        np.save(filename+"_comp", CompF)

        print("The result is automatically saved in "+filename+"_tot.npy")
        print("The result is automatically saved in "+filename+"_comp.npy")
        
        if output:
            return TotalF, CompF


    def __Reset__(self, newP = []):

        self._model_names = Utilities.models.model_list

        self._modelset = Utilities.models.model_xspec_st

        self._setP = Utilities.models.pars_xspec_st

        self._newP = newP


    def __Fitting__(self, newInit, debug=False, BAT=False, XRT=False, dataType=1, rspType=True):    
        from xspec import AllData, AllModels, Xset, Model, Fit, Plot
        
        AllData.clear()
        AllModels.clear()
        Xset.allowPrompting=False

        self.__Reset__(newP = self._newP)
        self._AllData = AllData
        self._AllModels = AllModels
        self._data_num = 1
        self._group_num = []
        self._dof = 0
        
        AllModels.addPyMod(Utilities.functions.bkn2power, Utilities.functions.bkn2pow_Info, 'add')
        AllModels.addPyMod(Utilities.functions.multiBB, Utilities.functions.multiBB_Info, 'add')
         
        self._usedData = self.__SelectedData__(self._only)
        
        if self._only != 'Swiftonly': self.__PrepGBMLLE__(dataType=dataType, rspType = rspType)

        if 'LAT' in self._usedData: self.__PrepLAT__(dataType=dataType)
        
        if BAT: self.__PrepBAT__()

        if XRT: self.__PrepXRT__()


        self._usedData = self.__SelectedData__(self._only, BAT=BAT, XRT=XRT)

        # Set modelname
        models = self.__RearrModel__(self._m_num)

        if XRT: 
            models = models+'* TBabs * zTBabs'

        # Input the model to Xspec
        if newInit and len(self._newP)>0:
            changedInitP = np.asarray(self._newP)
            if m_num in changedInitP[:,0]:
                paramSET = changedInitP[:,1][changedInitP[:,0] == self._m_num][0]
                
            else:
                paramSET = self._setP[self._m_num]
        else:
            paramSET = self._setP[self._m_num]

        Mmanager=Model(models, setPars=paramSET) 

        self._group_num = np.asarray(self._group_num)

        if XRT:
            paramSET.update({Mmanager.nParameters: '{}, 0.'.format(self.redshift)})
            paramSET.update({Mmanager.nParameters-2: '{}, 0.'.format(self._nH)})
            Mmanager=Model(models, setPars=paramSET) 
        
        # Constant factor
        if hasattr(self, 'const'):
            if self.const:
                Mmanager=Model(models+'* const', setPars=paramSET)
                Mmanager.constant.factor.frozen=True
                #if XRT:
                #    if self._only != "Swiftonly":
                #        AllModels(np.int(self._group_num[self._group_num[:,0] == 'XRT'][0][1])).setPars({Mmanager.nParameters:'1.0, 0.1, 0.67, 0.67, 1.5, 1.5'})
                #    else:
                #        Mmanager.constant.factor.frozen=True
                if BAT:
                    AllModels(np.int(self._group_num[self._group_num[:,0] == 'BAT'][0][1])).setPars({Mmanager.nParameters:'1.0, 0.1, 0.67, 0.8, 1.25, 1.5'})
            else:
                Mmanager=Model(models, setPars=paramSET)

        

        if hasattr(self, 'sync'):
            if self.sync:
                if 'bknpower' in models:
                    Mmanager.bknpower.PhoIndx2.link = "-0.5 + 1"
                elif 'grbm' in models:
                    Mmanager.grbm.beta.link = "-0.5 + 1"
                
            else:
                if 'bknpower' in models:
                    Mmanager.bknpower.PhoIndx2 = ""
                elif 'grbm' in models:
                    Mmanager.grbm.beta.link = ""

        # Fitting
        if XRT:
            Mmanager.zTBabs.Redshift.frozen=True
            Mmanager.TBabs.nH.frozen=True

            Xset.abund = "wilm"
            Xset.xsect = "vern"
        
        Fit.statMethod = 'cstat'
        
        if self._only == 'Swiftonly':
            if XRT and BAT:
                Fit.statMethod = "chi 1"
                Fit.statMethod = "cstat 2"
            elif BAT:
                Fit.statMethod = "chi"
            else:
                Fit.statMethod = "cstat"
        else:
            if XRT and BAT:
                Fit.statMethod = "pgstat"
                Fit.statMethod = "chi {}".format(int(self._data_num-2))
                Fit.statMethod = "cstat {}".format(int(self._data_num-1))
            elif BAT:
                Fit.statMethod = "pgstat"
                Fit.statMethod = "chi {}".format(int(self._data_num-1))
            elif XRT:
                Fit.statMethod = "pgstat"
                Fit.statMethod = "cstat {}".format(int(self._data_num-1))
            else:
                Fit.statMethod = "pgstat"
        
        Fit.nIterations = 10000

        Fit.perform()
    
        Xset.parallel.error = Mmanager.nParameters
        
        errPar = ''
        for i in range(Mmanager.nParameters):
            errPar+='{} '.format(i+1)
        Fit.error('{}'.format(errPar))
    
        for i in range(Mmanager.nParameters):
            Fit.error('1.0 {}'.format(i+1))

        x, y, m = self.__PlotSetting__(debug = debug, BAT=BAT, XRT=XRT)
        
        self._Mmanager = Mmanager
        self._Fit = Fit

        n = self._dof
        k = self._dof-Fit.dof
    
        bic= Fit.statistic + k*np.log(n)

        Statistic = [Fit.statistic, Fit.dof, bic]

        output = {"Model":self.__FitResult__(), "TimeInterval":self.toi, "Statistics":Statistic, "Covariance":Fit.covariance}
        dataPoints = {"TimeInterval":self.toi, "PlotX":x, "PlotY":y, "PlotModel": m}

        return output, dataPoints
            
    def __PrepGBMLLE__(self, dataType=1, rspType = True):

        # Import GBM and LLE dataset

        for det in self._usedData:

            if det == 'LAT': break

            if self._ti_num == -1:
                pha_file = './{}-{}-0.pha'.format(self._address_Xspec, det)
                bak_file = './{}-{}-0.bak'.format(self._address_Xspec, det)
            else:
                pha_file = './{}-{}-{}.pha'.format(self._address_Xspec, det, self._ti_num)
                bak_file = './{}-{}-{}.bak'.format(self._address_Xspec, det, self._ti_num)
            
            if self._data_num == 1:
                if dataType == 1:
                    self._AllData(pha_file)
                elif dataType == 2:
                    self._AllData(pha_file+ '{1}')
            else:
                if dataType == 1:
                    self._AllData('{}:{} '.format(self._data_num, self._data_num)+pha_file)
                elif dataType== 2:
                    self._AllData('{}:{} '.format(self._data_num, self._data_num)+pha_file+ '{1}')
            
            with fits.open(pha_file, mode='update') as PHA:
                PHA['SPECTRUM'].header['POISSERR'] = True
            
            if det !='LLE':
                if rspType:
                    self._AllData(self._data_num).response = './{}-{}-{}.rsp2'.format(self._address_Xspec, det, self._ti_num)
                else:
                    self._AllData(self._data_num).response = './{}-{}-{}.rsp'.format(self._address_Xspec, det, self._ti_num)
            else:
                if rspType:
                    self._AllData(self._data_num).response = './{}-{}-{}.rsp2'.format(self._address_Xspec, det, self._ti_num)
                else:
                    self._AllData(self._data_num).response = './{}-{}-{}.rsp'.format(self._address_Xspec, det, self._ti_num)

            try:
                with fits.open(bak_file, mode='update') as BAK:
                    BAK['SPECTRUM'].header['POISSERR'] = False
                self._AllData(self._data_num).background = './{}-{}-{}.bak'.format(self._address_Xspec, det, self._ti_num)
            except:
                self._AllData(self._data_num).background = './{}-{}-{}.bak'.format(self._address_Xspec, det, self._ti_num)

            if det == 'b0' or det == 'b1': 
                self.__IgnoreChan__(det)
                
            elif det == 'LLE': 
                self.__IgnoreChan__(det)
                
            else: 
                self.__IgnoreChan__(det)
                self._AllData(self._data_num).ignore("30.-40.")

            
            self._dof+=len(self._AllData(self._data_num).energies)
            self._group_num.append([det, self._data_num])
            self._data_num+=1
            
    
    def __PrepLAT__(self, dataType=1):
        
        # Import LAT dataset

        if self._ti_num == -1:
            pha_file = './{}-LAT-0.pha'.format(self._address_Xspec)
            rsp_file = './{}-LAT-0.rsp'.format(self._address_Xspec)
            bak_file = './{}-LAT-0.bak'.format(self._address_Xspec)
        else:
            pha_file = './{}-LAT-{}.pha'.format(self._address_Xspec, self._ti_num)
            rsp_file = './{}-LAT-{}.rsp'.format(self._address_Xspec, self._ti_num)
            bak_file = './{}-LAT-{}.bak'.format(self._address_Xspec, self._ti_num)
            

        with fits.open(pha_file, mode='update') as PHA:
            PHA['SPECTRUM'].header['POISSERR'] = True

        # Source file(Neglect background)
        if self._data_num == 1:
            self._AllData(pha_file)
        else:
            self._AllData('{}:{} '.format(self._data_num, self._data_num)+pha_file)

        # Response file
        self._AllData(self._data_num).response = rsp_file

        # Background file
        if dataType==1:
            with fits.open(bak_file, mode='update') as BAK:
                BAK['SPECTRUM'].header['POISSERR'] = False
            self._AllData(self._data_num).background = bak_file

        if hasattr(self, 'energy_range'):
            energy_range = self.energy_range
            self._AllData(self._data_num).ignore("**-{:.1f}".format(energy_range[0]))
            self._AllData(self._data_num).ignore("{:.1f}-**".format(energy_range[1]))
        else:
            energy_range = [100000.,100000000.]

        self._dof+=len(self._AllData(self._data_num).energies)
        self._group_num.append(["LAT", self._data_num])
        self._data_num+=1
    '''
    def __PrepXRT__(self):

        # Import XRT dataset

        os.chdir("./XRT/")
        if ti_num == 7: txrt = 1
        elif ti_num == 8: txrt = 2
        elif ti_num == 9: txrt = 3
        elif ti_num == 10: txrt = 4
        elif ti_num == 11: txrt = 'Ma5'
        elif ti_num == 12: txrt = '1_2'

        if ti_num == 11:
            self._AllData('{}:{} 77to180wt.pi'.format(self._data_num, self._data_num))
        else:
            if self._data_num == 1:
                self._AllData('sw00883832000xwtw2po_gti{}_g0_gr1.pi'.format(txrt))
            else:
                self._AllData('{}:{} sw00883832000xwtw2po_gti{}_g0_gr1.pi'.format(self._data_num, self._data_num, txrt))

        if hasattr(self, 'energy_range'):
            energy_range = self.energy_range
            self._AllData(self._data_num).ignore("**-{:.1f}".format(energy_range[0]))
            self._AllData(self._data_num).ignore("{:.1f}-**".format(energy_range[1]))
            #self._AllData(self._data_num).ignore("1-{} {}-{}".format(len(np.asarray(AllData(self._data_num).energies)[:,0])-sum(np.asarray(AllData(self._data_num).energies)[:,0]>energy_range[0]), sum(np.asarray(AllData(self._data_num).energies)[:,0]<energy_range[1])+1, len(np.asarray(AllData(self._data_num).energies)[:,0])))
        else:
            self._AllData(self._data_num).ignore("**-1.")
            self._AllData(self._data_num).ignore("10.-**")

        os.chdir("../")
        self._dof+=len(self._AllData(self._data_num).energies)
        self._group_num.append(["XRT", self._data_num])
        self._data_num+=1


    def __PrepBAT__(self, ti_num):

        # Import BAT dataset

        if self._data_num == 1:
            self._AllData('./BAT/{}-{}-{}.pha'.format(self.event_name, 'BAT', ti_num))
        else:
            self._AllData('{}:{} ./BAT/{}-{}-{}.pha'.format(self._data_num, self._data_num,self.event_name, 'BAT', ti_num))
        self._AllData(self._data_num).response = './BAT/{}-{}-{}.rsp'.format(self.event_name, 'BAT', ti_num)
        
        if hasattr(self, 'energy_range'):
            energy_range = self.energy_range
        else:
            energy_range = [15.,150.]
        self._AllData(self._data_num).ignore("**-{:.1f}".format(energy_range[0]))
        self._AllData(self._data_num).ignore("{:.1f}-**".format(energy_range[1]))

        self._dof+=len(self._AllData(self._data_num).energies)
        self._group_num.append(["BAT", self._data_num])
        self._data_num+=1
    '''    
    def __IgnoreChan__(self, det):

        if hasattr(self, 'energy_range'):
            if det in ['b0', 'b1']:
                if self.energy_range[0]<250.0:
                    energy_range =  [250.0, self.energy_range[1]]
                else:
                    energy_range = self.energy_range
            elif det == 'LLE':
                energy_range = self.energy_range
                #energy_range = [30000.0, 100000.0]
            else:
                if self.energy_range[1]>900.0:
                    energy_range = [self.energy_range[0], 900.0]
                else:
                    energy_range = self.energy_range
        else:
            if det in ['b0', 'b1']:
                energy_range =  [260.0, 38000.0]
            elif det == 'LLE':
                energy_range = [30000.0, 100000.0]
            else:
                energy_range = [10.0, 900.0]

        self._AllData(self._data_num).ignore("**-{:.1f}".format(energy_range[0]))
        self._AllData(self._data_num).ignore("{:.1f}-**".format(energy_range[1]))
 
    
    def __RearrModel__(self, m_num):

        # Input: ["Model_1", "Model_2"]
        # Output: "Model_1+Model_2"
        
        models = ''
        for mn in self._modelset[m_num]:
            models=models+'+'+mn
        models=models[1:]
        return models
    
    def __SelectedData__(self, only, BAT=False, XRT=False):
        if only == 'LATonly': Detectors=['LAT']
        elif only == 'GBMonly': Detectors=self.usedGBM
        elif only == 'LLEonly': Detectors=['LLE']    
        elif only == 'LATLLE': Detectors=['LLE','LAT']
        elif only == 'GBMLLE': Detectors=self.usedGBM + ['LLE']  
        elif only == 'GBMLAT': Detectors=self.usedGBM + ['LAT']  
        elif only == 'Swiftonly': Detectors=[]
        else: Detectors = self.usedGBM+['LLE','LAT'] 
        if BAT : Detectors = ['BAT'] + Detectors
        if XRT: Detectors = Detectors+['XRT']

        return Detectors

    def __PlotSetting__(self, debug = False, BAT=False, XRT=False):

        Plot.commands = ()
        Plot.xLog = True
        Plot.yLog = True
        Plot.device = "/xs"
        Plot.xAxis="keV"
        #Plot('data', "eemodel")
        
        for det, groupN in self._group_num:
            if 'n' in det or 'b' in det:
                Plot.setRebin(1, 1, groupNum="{}".format(groupN))
            elif det == 'XRT':
                Plot.setRebin(10, 100, groupNum="{}".format(groupN))
            elif det == 'BAT':
                Plot.setRebin(5, 10, groupNum="{}".format(groupN))
            #elif det == 'LAT':
            #    Plot.setRebin(1.5, 10, groupNum="{}".format(groupN))
        

        if not(BAT) and not(XRT):
            Plot('eeufspec', "delchi")
            #Plot.addCommand("la y \gn F\d\gn\u (keV cm\u-2\d s\u-1\d)")
            Plot.addCommand("wind 2")
            #Plot.addCommand("la y ratio")
            Plot.addCommand("rescale y -5 5")
            Plot.addCommand("LIne OFf")
            Plot.addCommand("Uncertainty y ON")
            Plot.addCommand("Error y ON")
            Plot.iplot()
        else:
            Plot('data', 'residuals')
            Plot.iplot()
        
        self._Plot = Plot
        x = []
        y = []
        m = []
        
        #if not(BAT):
        if True:
            
            for i in range(1,np.size(self._usedData)+1):
                x += [[median, sigma] for median, sigma in zip(Plot.x(plotGroup=i), Plot.xErr(plotGroup=i))]
                y += [[median, sigma] for median, sigma in zip(Plot.y(plotGroup=i), Plot.yErr(plotGroup=i))]
                m += [value for value in Plot.model(plotGroup=i)]

            x = np.asarray(x)
            y = np.asarray(y)
            m = np.asarray(m)
        
        return x, y, m
    
    def __MakeTable__(self, clear=False, norm=True, latex=False, Epiv=100):
        multi = False

        if clear or not(hasattr(self, "_Ftable")):
            self.__ClearTable__()
            newT = True
        else:
            newT = False
            

        for fresult in self._current_output["Model"]:
            self._FparN +=fresult['pars']
        self._FparN = list(set(self._FparN))
        self._FparN.sort()

        if not(norm): 
            self._FparN.remove('norm')

        for fresult in self._current_output["Model"]:
            for par in self._FparN:
                if par not in self._Ftable.keys():
                    self._Ftable.update({par:[]})
                    if not(newT):
                        for i in range(len(self._Ftable['Model'])):
                            self._Ftable[par].append('')
                try:
                    #if par == 'norm' and fresult['name'] == 'powerlaw':
                    #    self._Ftable[par].append(r"{:.3f}$^{{+{:.3f}}}_{{-{:.3f}}}$".format(Ampli[0], Ampli[1]-Ampli[0], Ampli[0]-Ampli[2]))
                    if par == 'ep2' and fresult['name'] == 'grbm*highcut':
                        Ep2 = self.__FoldEtoEp__()
                        self._Ftable[par].append(r"{:.3f}$^{{+{:.3f}}}_{{-{:.3f}}}$".format(Ep2[0], Ep2[1]-Ep2[0], Ep2[0]-Ep2[2]))
                    else:
                        if (fresult[par+'_err'][1] == 0.) and (fresult[par+'_err'][0] == 0.) :
                            self._Ftable[par].append(r"{:.3f}$_{{fixed}}$".format(fresult[par]))
                        else:
                            self._Ftable[par].append(r"{:.3f}$^{{+{:.3f}}}_{{-{:.3f}}}$".format(fresult[par], fresult[par+'_err'][1]-fresult[par], fresult[par]-fresult[par+'_err'][0]))
                except:
                    self._Ftable[par].append("")

            if fresult['name'] == 'step':
                modelName = 'cutoffpl'
            else:
                modelName = fresult['name']

            if multi:
                self._Ftable['Model'].append(modelName)
                self._Ftable['Time'].append(['',''])
                self._Ftable['Stats'].append(['','',''])
                self._Ftable['T_num'].append('')
            else:
                self._Ftable['Model'].append(modelName)
                self._Ftable['Time'].append(["{:.2f}".format(self.toi[0]), "{:.2f}".format(self.toi[1])])
                self._Ftable['Stats'].append(["{:.1f}".format(self._current_output["Statistics"][0]), self._current_output["Statistics"][1], "{:.1f}".format(self._current_output["Statistics"][2])])
                self._Ftable['T_num'].append(self._ti_num)

            multi = True

    def __ClearTable__(self):
        self._Ftable = {}
        self._FparN = []
        self._Ftable.update({"Model":[]})
        self._Ftable.update({"Stats":[]})
        self._Ftable.update({"Time":[]})
        self._Ftable.update({"T_num":[]})

    def __TableList__(self, latex=False):
        tList = []
        for i in range(len(self._Ftable['Model'])):
            temp = []
            if latex:
                parlist = ['Time','Model']+self._FparN+['Stats']
            else:
                parlist = ['Time','Model']+self._FparN+['Stats', 'T_num']
            for par in parlist:
                if np.size(self._Ftable[par][i])==1:
                    temp+=[self._Ftable[par][i]]
                else:
                    temp+=self._Ftable[par][i]
            tList.append(temp)
        return tList

    def __FoldEtoEp__(self):
        # Input: E_f
        # Output: Ep = (alpha + 2)*E_f
        parSet = self.__GeneratePars__(sim=10000)
        Ep2 = (parSet[:,1]+2.)*parSet[:,5]
        return [np.percentile(Ep2, 50), np.percentile(Ep2, 84),np.percentile(Ep2, 16)]


    def __MakeSED__(self, lc_flag = False, save_flag = False, time_range = [], energy_range = [50, 300], overlap=True, verbose=False):
        if not(lc_flag):
            fig = plt.figure(figsize=(5,6))
            gs = plt.GridSpec(2, 1, height_ratios=[3, 12])
            gs.update(hspace=0.32)
            topax = fig.add_subplot(gs[0])

            topax = Utilities.lightcurve.LightCurve(self, ax = topax, toi=self.toi+self.trigger)
            
            gs_base = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], height_ratios=[9, 2], hspace=0.05)
            ax = [fig.add_subplot(gs_base[0,:]),fig.add_subplot(gs_base[1,:])]
        else:
            f, ax = plt.subplots(2, 1, figsize=(5,5), gridspec_kw = {'height_ratios':[7, 2]})
        x = self._current_data_points["PlotX"]
        y = self._current_data_points["PlotY"]
        m = self._current_data_points["PlotModel"]
        newArr = x[:,0].argsort()
        x = x[newArr]
        y = y[newArr]
        m = m[newArr]

        #try:
        modelF, model_compF, etc = self.__ModelFlux__(x[:,0], verbose=verbose)
        ratio = np.sign(y[:,0]-x[:,0]**2.*modelF)*abs(np.asarray(y[:,0]-x[:,0]**2.*modelF)/abs(y[:,1]))
        ratio_err = np.zeros(len(x[:,0]))+1
        #except:
        #    print("At least one model is not in the list")
        #    modelF = m
        #    ratio = np.sign(y[:,0]-m)*abs(np.asarray(y[:,0]-m)/abs(y[:,1]))
        #    ratio_err = np.zeros(len(x[:,0]))+1

        if overlap: ax[0].errorbar(x[:,0], y[:,0], xerr=x[:,1], yerr=y[:,1], fmt='', ls='', c='gray', elinewidth=1)
        ax[1].errorbar(x[:,0], ratio, yerr=ratio_err, xerr=x[:,1], fmt='', ls='', c='gray', elinewidth=1)

        if 'XRT' in self._usedData:
            E=np.logspace(-1, 7, 150)
        else:
            E=np.logspace(1, 7, 100)

        try:
            modelF, model_compF, etc = self.__ModelFlux__(E, verbose=verbose)
            ax[0].loglog(E, E**2*modelF, color='k', label=Utilities.models.model_list[self._m_num])
            for modelF in model_compF:
                ax[0].loglog(E, E**2*modelF, color='k', ls=":", lw=0.7, )
        except:
            ax[0].loglog(x[:,0], m, color='k', label=Utilities.models.model_list[self._m_num])
        
        if 'XRT' in self._usedData:
            ax[0].set_xlim(0.8, 1e8)
            ax[1].set_xlim(0.8, 1e8)
            ax[0].set_ylim(1e-2, 5e4)
        else:
            ax[0].set_xlim(10, 1e8)
            ax[1].set_xlim(10, 1e8)
            ax[0].set_ylim(2, 5e4)

        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        ax[1].set_xscale("log")
        ax[0].set_xticklabels([])
        ax[0].set_ylabel(r"$\nu$ F$\nu$ [keV cm$^{-2}$ s$^{-1}$]", fontsize=12)
        ax[1].set_xlabel("Energy [keV]", fontsize=12)
        ax[1].set_ylabel("Residual", fontsize=12)
        ax[1].set_ylim(-3, 3)
        ax[1].axhline(0, color='k', ls=":")
        ax[0].grid()
        ax[1].grid()
        ax[0].legend()

        if save_flag:
            plt.save_flag("{}-{}.pdf".format(self._ti_num,self._m_num), bbox_inches="tight")
        
        plt.show(block=False)

    def __ModelFlux__(self, E, erg=False, butterfly = False, verbose=False, XRT=False):
        parSet = self.__GeneratePars__(sim=butterfly)

        synF = []
        
        if XRT:
            parSet = [[x, y] for x, y in zip(parSet[:,0], parSet[:,1])]
        
        for parS in parSet:
            CompF = []
            TotalF = np.zeros(len(E))
            parNum = 0
            for fresult in self._current_output['Model']:
                model_name = fresult['name']
                if model_name =='TBabs' or model_name == 'zTBabs': continue

                func = globals()[model_name.upper()]
                thisParNum = func.__code__.co_argcount-1

                nufnu = func(E, *parS[np.arange(thisParNum)+parNum])
                
                CompF.append(nufnu)
                TotalF+=np.asarray(nufnu)
                parNum += thisParNum

            synF.append(TotalF)
            if verbose :
                print("{:.1f} %".format(len(synF)*100./len(parSet)), end="\r")
        
        TotalF_lu = []
        if butterfly:
            synF = np.asarray(synF)

            for i in range(len(E)):
                TotalF_lu.append([np.percentile(synF[:,i],16), np.percentile(synF[:,i], 84)])
        if erg:
            TotalF = np.asarray(TotalF)*1.6022*10**-9
            CompF = np.asarray(CompF)*1.6022*10**-9
            TotalF_lu = np.asarray(TotalF_lu)*1.6022*10**-9
        return TotalF, np.asarray(CompF), np.asarray(TotalF_lu)

    def __FitResult__(self):
        fitPar = self._Mmanager
        FitOutput = []
        for comp in fitPar.componentNames:
            parameters = fitPar.__getattribute__(comp)
            if comp == 'highecut':
                temp = FitOutput[-1]
                temp['name'] = temp['name']+'*highcut'
            else:
                temp = {}
                temp['name'] = comp
                temp['pars'] = []
                temp['frozen'] = []
            for parN in parameters.parameterNames:
                par = parameters.__getattribute__(parN)
                if parN == 'PhoIndex' and len(parameters.parameterNames) == 2:
                    parN = 'gamma'
                elif parN == 'PhoIndex' or parN == 'PhoIndx1':
                    parN = 'alpha'
                elif parN == 'PhoIndx2':
                    parN = 'beta'
                elif parN == 'PhoIndx3':
                    parN = 'beta2'
                elif parN == 'Sigma' or parN == 'tem' or parN == 'BreakE' or parN == 'BreakE1' or parN == 'HighECut':
                    parN = 'ep'
                elif parN == 'cutoffE':
                    parN = 'cutoffE'
                elif parN == 'BreakE2':
                    parN = 'ep2'
                if parN == 'Energy':
                    parN = 'alpha'
                    temp['pars'].append(parN)
                    temp[parN] = -par.values[0]
                    temp[parN+'_err'] = (-np.asarray(par.error[:-1]))[::-1]
                elif parN == 'factor':
                    temp['pars'].append(parN)
                    temp[parN] = self._AllModels(self._data_num-1).constant.factor.values[0]
                    temp[parN+'_err'] = [self._AllModels(self._data_num-1).constant.factor.values[0]-abs(self._AllModels(self._data_num-1).constant.factor.sigma),self._AllModels(self._data_num-1).constant.factor.values[0]+abs(self._AllModels(self._data_num-1).constant.factor.sigma)]
                else:
                    temp['pars'].append(parN)
                    temp[parN] = par.values[0]
                    temp[parN+'_err'] = np.asarray(par.error[:-1])

                if par.frozen and parN!='factor':
                    temp['frozen'].append(parN)
            
            if comp != 'highecut':
                FitOutput.append(temp)
        return FitOutput

    def __MaxEnergy__(self, fileN):
        DATA = fits.open("./Xspec/{}-mktime-{}.fits".format(self.event_name,fileN))[1].data
        maxE = DATA.field("Energy")
        if len(maxE) == 0:
            Chan = fits.open("./Xspec/{}-LLE-{}.pha".format(self.event_name,fileN))[2].data
            DATA = fits.open("./Xspec/{}-LLE-{}.pha".format(self.event_name,fileN))[1].data.field('Counts')[0]
            for i, c in zip(DATA[::-1], Chan[::-1]):
                if i>0:
                    maxE = c[2]
                    break
            if np.size(maxE) == 0:
                maxE = 40000
            elif maxE > 1e6:
                maxE = 1e6
        else:
            maxE = max(maxE)*1e3
        return np.log10(maxE)

    def __CovMatrix__(self):
        temp=0
        for matSize in range(20):
            temp+=matSize
            if np.size(self._current_output['Covariance']) == temp:
                break

        covReshape = np.zeros(shape=(matSize,matSize))

        k=0
        for j in range(matSize):
            for i in range(matSize):
                if i <= j and covReshape[i][j] == 0: 
                    covReshape[i][j]=self._current_output['Covariance'][k]
                    covReshape[j][i]=self._current_output['Covariance'][k]
                    k+=1
                    continue
                continue

        return covReshape

    def __GeneratePars__(self, runs = 10000, sim=False):
        mean = []
        fixed_mean = []
        for pars in self._current_output['Model']:
            for par in pars['pars']:
                if hasattr(self, 'sync'):
                    if self.sync:
                        if par == 'beta' or par == 'PhoIndx2' or par in pars['frozen']:
                            fixed_mean.append(False)
                        else:
                            fixed_mean.append(True)
                else:
                    if par in pars['frozen']:
                        fixed_mean.append(False)
                    else:
                        fixed_mean.append(True)
                mean.append(pars[par])

        if sim:
            cov = self.__CovMatrix__()
            parSet = np.random.multivariate_normal(np.asarray(mean)[fixed_mean], cov, runs-1)

            if fixed_mean[0]:
                new_parSet = parSet[:,0]
                i=1
            else:
                new_parSet = np.zeros(runs-1)+mean[0]
                i=0

            for fm, j in zip(fixed_mean[1:], range(1,len(fixed_mean))):
                if fm:
                    new_parSet = np.vstack((new_parSet, parSet[:,i]))
                    i+=1
                else:
                    new_parSet = np.vstack((new_parSet,  np.zeros(runs-1)+mean[j]))
            parSet = np.asarray(new_parSet.transpose().tolist()+[mean])
        else:
            parSet = np.asarray([mean])

        return parSet

