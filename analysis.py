import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import glob

from xspec import *
import XspecModels
import XspecFunctions
from XspecFunctions import *
import XspecBackground
from XspecTable import *
from XspecUtils import *

from astropy.io import fits

import gt_apps
from gt_apps import *
from GtApp import GtApp

from scipy.integrate import quad
from datetime import datetime
bindef = GtApp('gtbindef','Likelihood')

from IPython.display import clear_output



class Xspec_T:

    def __init__(self, EvtInfo={}, verbose=False, only='', folder = ''):
        if verbose:
            print("config = {}")
            print("config['Name'] = '131108A'")
            print("config['Trigger'] = 405636118")
            print("config['RA'] = 156.50187")
            print("config['DEC'] = 9.66225")
            print("config['T90'] = [0.319986, 18.496286]")
            print("config['usedGBM'] = ['n0','n3','n6','n7','b0','b1']")
            print("config['timeintervals'] = cuts")
            print("Fitting = Xspec_T(EvtInfo=config) or Fitting.configuration(config)")
            print("Fitting.runFit(1,1)")
            print("Fitting.printSED()")
            print("Fitting.printTable()")

        if len(EvtInfo)!=0:
            self.configuration(EvtInfo, only=only, folder=folder)

    def configuration(self, config, only='', folder = ''):
        self._GRBname = config['Name']
        self._T0 = config['Trigger']
        self._T90 = config['T90']
        self._refRa = config['RA']
        self._refDec = config['DEC']
        self._usedGBM = config['usedGBM']
        self._timeInt = np.asarray(config['timeintervals'])
        if 'redshift' in config.keys():
            self._redshift = config['redshift']
        if 'nH' in config.keys():
            self._nH = config['nH']
        self._only = only
        try:
            self._tS_set = self._T0+np.asarray(config['t_start'])
            self._tE_set = self._T0+np.asarray(config['t_end'])
        except:
            self._tS_set = self._T0+np.asarray(self._timeInt)[:-1]
            self._tE_set = self._T0+np.asarray(self._timeInt)[1:]
            
        self._usedData = self.__SelectedDet__(self._only)
        self._address = folder
        self._addressToData = self._address+"/Xspec/"+self._GRBname+''
        self.__reset__()


    def setInit(self, mod_num, init = {}, verbose=True):
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
        models = self.__rearrModel__(mod_num)
        Model_Xspec = Model(models)
        parN = 1
        for comp in Model_Xspec.componentNames:
            parameters = Model_Xspec.__getattribute__(comp)
            if len(init)==0: 
                print(comp, parameters.parameterNames)
                print(self._setP[mod_num])
            else:
                for par in parameters.parameterNames:
                    try:
                        if verbose: print(comp, par, init[parN])
                    except:
                        if verbose: print(comp, par)
                    parN+=1

        changedInitP = np.asarray(self._newP)
        if len(self._newP)>0:
            if mod_num not in changedInitP[:,0] and init!={}:
                self._newP.append([mod_num, init])
            else:
                changedInitP[np.where(changedInitP[:,0] == mod_num)[0][0]][1] = init
                self._newP = changedInitP.tolist()
        else:
            self._newP.append([mod_num, init])

    @classmethod
    def printModels(self):
        makeMtable = []
        i = 0
        j = 0
        temp = []
        for mod in XspecModels.ModelList:
            temp.append(j)
            temp.append(mod)
            i+=1
            if i >4:
                makeMtable.append(temp)
                temp = []
                i=0
            j+=1
        makeMtable.append(temp)

        return ModelTable(makeMtable)

    def printTimeItv(self):
        makeTtable = []
        i = 0
        j = 0
        temp = []
        for ti, te in zip(self._tS_set, self._tE_set):
            temp.append(j)
            temp.append("{:.2f} - {:.2f}".format(ti-self._T0, te-self._T0))
            i+=1
            if i >4:
                makeTtable.append(temp)
                temp = []
                i=0
            j+=1
        makeTtable.append(temp)
        return IntvTable(makeTtable)

    def saveFit(self, filename):
        self.OUTPUT = np.save(filename, self.OUTPUT)

    def loadFit(self, filename):
        self.OUTPUT = np.load(filename+'.npy')

    def runFit(self, mod_num = 0, ti_num = 0, newInit = False, debug = False, paired = False, BAT=False, XRT=False, dtType=1, rspT = True,verbose=True,  **kwargs):
        self.__dict__.update(kwargs)
        try:
            self._GRBname
        except:
            print("Initiate configuration first")
            self.__init__(verbose=True)
            sys.exit() 

        if paired:
            if np.size(mod_num)>2:
                self._mod_num = np.asarray(mod_num)[:,0]
                self._ti_num = np.asarray(mod_num)[:,1]
            else:
                self._mod_num = [mod_num[0]]
                self._ti_num = [mod_num[1]]
        elif np.size(ti_num) == np.size(mod_num):
            if np.size(mod_num) == 1:
                self._ti_num = [ti_num]
                self._mod_num = [mod_num]
            else:
                self._ti_num = ti_num
                self._mod_num = mod_num
        elif ti_num == -1:
            self._ti_num = (np.zeros(np.size(mod_num)) -1).astype(int)
            self._mod_num = mod_num
        elif np.size(ti_num) == 1 and np.size(mod_num) != 1:
            self._ti_num = (np.zeros(np.size(mod_num)) + ti_num).astype(int)
            self._mod_num = mod_num
        elif np.size(ti_num) >= 1 and np.size(mod_num) >= 1:
            self._ti_num = []
            self._mod_num = []
            for t in ti_num:
                for m in mod_num:
                    self._ti_num.append(t)
                    self._mod_num.append(m)
        else:
            sys.exit()
        
        self.OUTPUT = []
        if verbose and ti_num!=-1: print("{:^25} | {:^15} | {:^27}".format("Model","Time Intv","Status"))
        sys.stdout.flush()

        for mod_num, ti_num in zip(self._mod_num, self._ti_num):
            if verbose and ti_num!=-1: print("{:^25} | {:>6.2f} - {:<6.2f} | {:^27} \r".format(self._modNames[mod_num], self._tS_set[ti_num]-self._T0, self._tE_set[ti_num]-self._T0, 'Working'))
            output = self.__Xspec__(mod_num, ti_num, newInit, debug=debug, BAT=BAT, XRT=XRT, dtType=dtType, rspT=rspT)
            self.OUTPUT.append(output)
            
            if verbose: print("{:^25} | {:>6.2f} - {:<6.2f} | {:^27}".format(self._modNames[mod_num], self._tS_set[ti_num]-self._T0, self._tE_set[ti_num]-self._T0, 'Done'))
        self._output = self.OUTPUT[-1]

        
    def printTable(self, outputN = -1, clear=False, latex=False, norm=True, Epiv=100):
        if np.size(self.OUTPUT) == 1:
            self._output = self.OUTPUT[0]
            self.__makeTable__(self._ti_num[0], clear=clear, norm=norm, latex=latex, Epiv=Epiv)
            
        elif outputN!=-1:
            self._output = self.OUTPUT[outputN]
            self.__makeTable__(self._ti_num[outputN], clear=clear, norm=norm, latex=latex, Epiv=Epiv)
            
        else:
            self.__clearTable__()
            for output, ti_num in zip(self.OUTPUT, self._ti_num):
                self._output = output
                self.__makeTable__(ti_num, clear=False, norm=norm, latex=latex, Epiv=Epiv)
   
        if latex:
            FT = FitTable(self.__tableList__(latex=True), parlist=self._FparN, LaTex=latex)
            print(FT._repr_latex_())
        else:
            return FitTable(self.__tableList__(), parlist=self._FparN)

    
    def printSED(self, outputN = -1, noLC = False, saveFig = False, overlap=True, tRan = [], eRan = [50, 300], verbose=False): 
        if np.size(self.OUTPUT) == 1:
            self._output = self.OUTPUT[0]
            self.__makeSED__(self._mod_num, self._ti_num, noLC = noLC, saveFig = saveFig, tRan=tRan, eRan=eRan, overlap=overlap, verbose=verbose)
        elif outputN!=-1:
            self._output = self.OUTPUT[outputN]
            self.__makeSED__(self._mod_num[outputN], self._ti_num[outputN], noLC = noLC, saveFig = saveFig, tRan=tRan, eRan=eRan, overlap=overlap, verbose=verbose)
        else:
            for output, mod_num, ti_num in zip(self.OUTPUT, self._mod_num, self._ti_num):
                self._output = output
                self.__makeSED__(mod_num, ti_num, noLC = noLC, saveFig = saveFig, tRan=tRan, eRan=eRan, overlap=overlap, verbose=verbose)

    def printSED_tot(self, saveFig = False, verbose=False):
        fig, ax = plt.subplots(len(self.OUTPUT)+1,1,figsize=(10, 7), gridspec_kw = {'height_ratios':[7]+(np.zeros(len(self.OUTPUT))+1.5).tolist()}) 
        plt.subplots_adjust(hspace=0.01)

        for i in range(len(self.OUTPUT)):

            self._output = self.OUTPUT[i]
            emax = self.__MaxEnergy__(i)
            E = np.logspace(1, emax, 100)
            Label = r"Episode {} ({:.1f}s - {:.1f}s): {}".format(i+1, self._tS_set[i]-self._T0, self._tE_set[i]-self._T0, self._modNames[self._mod_num[i]])
            
            modelF, model_compF, modelF_lu=self.__modelFlux__(E, butterfly=True, erg=True, verbose=verbose)
            PLT, = ax[0].loglog(E, E**2*modelF, label=Label)
            ax[0].fill_between(E, E**2*modelF_lu[:,0], E**2*modelF_lu[:,1], facecolor='none', alpha=0.3, edgecolor=PLT.get_c(), linewidth=1, antialiased=True)

            x = self._output["PlotX"]
            y = self._output["PlotY"]
            newArr = x[:,0].argsort()
            y = y[newArr]
            x = x[newArr]

            modelF, model_compF, etc = self.__modelFlux__(x[:,0], verbose=verbose)
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
        ax[(len(self.OUTPUT)+1)/2].yaxis.set_label_coords(-0.07, 0.5)
        ax[-1].set_xlabel('Energy [keV]', fontsize=12)
        ax[-1].tick_params(labelsize=12, axis='x')

        if saveFig:
            plt.savefig("SED.pdf")
        plt.show(block=False)

    def calculateFlux(self, eRan = [10, 1e7], Runs = 100000, erg=False, filename='', output=False, phoflx=False):
        TotalF = []
        CompF = []
        ergCov = 1.6022*10**-9
        print("{:^25} | {:^15} | {:^27}".format("Model","Time Intv","Status"))
        for i in range(np.size(self.OUTPUT)):
            self._output = self.OUTPUT[i]
            parSet = self.__genPars__(runs = Runs, sim=True)
            synF = []
            synF_comp = []
            for parS, k in zip(parSet, range(Runs)):
                parN = 0
                totF = 0
                cF = []
                for fresult in self._output['Model']:
                    model_name = fresult['name']
                    if phoflx:
                        model = ModelCov[model_name]
                    else:
                        model = eModelCov[model_name]
                    flux = quad(model, eRan[0], eRan[1], args=tuple(parS[ParNCov[model_name]+parN]))[0]
                    cF.append(flux)
                    totF+=flux
                    parN = max(ParNCov[model_name])+1
                synF.append(totF)
                synF_comp.append(cF)
                k = int((k+1)*20./Runs)
                print("{:^25} | {:>6.2f} - {:<6.2f} | [{:<20}] {:.0f}%\r".format(self._modNames[self._mod_num[i]], self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0, '='*k, 5*k))
                sys.stdout.flush()
            
            print("{:^25} | {:>6.2f} - {:<6.2f} | {:^27}".format(self._modNames[self._mod_num[i]], self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0, 'Done'))

            if erg:
                TotalF.append([synF[-1]*ergCov, np.percentile(synF,16)*ergCov, np.percentile(synF, 84)*ergCov, self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0])
            else:
                TotalF.append([synF[-1], np.percentile(synF,16), np.percentile(synF, 84), self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0])
            
            CompF_temp = []

            if np.size(synF_comp, axis=1) > 1:
                synF_comp = np.asarray(synF_comp)
                CompF_temp = []
                for j in range(np.size(synF_comp, axis=1)):
                    if erg:
                        CompF_temp.append([synF_comp[:,j][-1]*ergCov, np.percentile(synF_comp[:,j],16)*ergCov, np.percentile(synF_comp[:,j], 84)*ergCov, self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0])
                    else:
                        CompF_temp.append([synF_comp[:,j][-1], np.percentile(synF_comp[:,j],16), np.percentile(synF_comp[:,j], 84)*ergCov, self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0])
                CompF.append(CompF_temp)
            else:
                if erg:
                    CompF.append([synF[-1]*ergCov, np.percentile(synF,16)*ergCov, np.percentile(synF, 84)*ergCov, self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0])
                else:
                    CompF.append([synF[-1], np.percentile(synF,16), np.percentile(synF, 84), self._tS_set[self._ti_num[i]]-self._T0, self._tE_set[self._ti_num[i]]-self._T0])
        
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


    def __reset__(self, newP = []):

        self._modNames = XspecModels.ModelList

        self._modelset = XspecModels.model_xspec_st

        self._setP = XspecModels.pars_xspec_st

        self._newP = newP
            
    def __prepGBMLLE__(self, ti_num, dtType=1, rspT = True):

        # Import GBM and LLE dataset

        for det in self._usedData:

            if det == 'LAT': break

            if ti_num == -1:
                fileAdd = './{}-{}.pha'.format(self._addressToData, det)
                fileAddB = './{}-{}.bak'.format(self._addressToData, det)
            else:
                fileAdd = './{}-{}-{}.pha'.format(self._addressToData, det, ti_num)
                fileAddB = './{}-{}-{}.bak'.format(self._addressToData, det, ti_num)
            
            if self._dataN == 1:
                if dtType == 1:
                    self._AllData(fileAdd)
                elif dtType == 2:
                    self._AllData(fileAdd+ '{1}')
            else:
                if dtType == 1:
                    self._AllData('{}:{} '.format(self._dataN, self._dataN)+fileAdd)
                elif dtType== 2:
                    self._AllData('{}:{} '.format(self._dataN, self._dataN)+fileAdd+ '{1}')
            
            with fits.open(fileAdd, mode='update') as PHA:
                PHA['SPECTRUM'].header['POISSERR'] = True
            
            if det !='LLE':
                if rspT:
                    self._AllData(self._dataN).response = './{}-{}.rsp2'.format(self._addressToData, det)
                else:
                    self._AllData(self._dataN).response = './{}-{}.rsp'.format(self._addressToData, det)
            else:
                if rspT:
                    self._AllData(self._dataN).response = './{}-{}.rsp2'.format(self._addressToData, det)
                else:
                    self._AllData(self._dataN).response = './{}-{}.rsp'.format(self._addressToData, det)

            try:
                with fits.open(fileAddB, mode='update') as BAK:
                    BAK['SPECTRUM'].header['POISSERR'] = False
                self._AllData(self._dataN).background = './{}-{}-{}.bak'.format(self._addressToData, det, ti_num)
            except:
                self._AllData(self._dataN).background = './{}-{}.bak'.format(self._addressToData, det)

            if det == 'b0' or det == 'b1': 
                self.__igChannel__(ti_num, det)
                
            elif det == 'LLE': 
                self.__igChannel__(ti_num, det)
                
            else: 
                self.__igChannel__(ti_num, det)
                self._AllData(self._dataN).ignore("30.-40.")

            
            self._dofData+=len(self._AllData(self._dataN).energies)
            self._groupN.append([det, self._dataN])
            self._dataN+=1
            
    
    def __prepLAT__(self, ti_num, dtType=1):
        
        # Import LAT dataset

        if ti_num == -1:
            fileAdd = './{}-LAT.pha'.format(self._addressToData)
            fileAddR = './{}-LAT.rsp'.format(self._addressToData)
            fileAddB = './{}-LAT.bak'.format(self._addressToData)
        else:
            fileAdd = './{}-LAT-{}.pha'.format(self._addressToData, ti_num)
            fileAddR = './{}-LAT-{}.rsp'.format(self._addressToData, ti_num)
            fileAddB = './{}-LAT-{}.bak'.format(self._addressToData, ti_num)
            

        with fits.open(fileAdd, mode='update') as PHA:
            PHA['SPECTRUM'].header['POISSERR'] = True

        # Source file(Neglect background)
        if self._dataN == 1:
            self._AllData(fileAdd)
        else:
            self._AllData('{}:{} '.format(self._dataN, self._dataN)+fileAdd)

        # Response file
        self._AllData(self._dataN).response = fileAddR

        # Background file
        if dtType==1:
            with fits.open(fileAddB, mode='update') as BAK:
                BAK['SPECTRUM'].header['POISSERR'] = False
            self._AllData(self._dataN).background = fileAddB

        if hasattr(self, 'eran'):
            eran = self.eran
            self._AllData(self._dataN).ignore("**-{:.1f}".format(eran[0]))
            self._AllData(self._dataN).ignore("{:.1f}-**".format(eran[1]))
        else:
            eran = [100000.,100000000.]

        self._dofData+=len(self._AllData(self._dataN).energies)
        self._groupN.append(["LAT", self._dataN])
        self._dataN+=1

    def __prepXRT__(self, ti_num):

        # Import XRT dataset

        os.chdir("./XRT/")
        if ti_num == 7: txrt = 1
        elif ti_num == 8: txrt = 2
        elif ti_num == 9: txrt = 3
        elif ti_num == 10: txrt = 4
        elif ti_num == 11: txrt = 'Ma5'
        elif ti_num == 12: txrt = '1_2'

        if ti_num == 11:
            self._AllData('{}:{} 77to180wt.pi'.format(self._dataN, self._dataN))
        else:
            if self._dataN == 1:
                self._AllData('sw00883832000xwtw2po_gti{}_g0_gr1.pi'.format(txrt))
            else:
                self._AllData('{}:{} sw00883832000xwtw2po_gti{}_g0_gr1.pi'.format(self._dataN, self._dataN, txrt))

        if hasattr(self, 'eran'):
            eran = self.eran
            self._AllData(self._dataN).ignore("**-{:.1f}".format(eran[0]))
            self._AllData(self._dataN).ignore("{:.1f}-**".format(eran[1]))
            #self._AllData(self._dataN).ignore("1-{} {}-{}".format(len(np.asarray(AllData(self._dataN).energies)[:,0])-sum(np.asarray(AllData(self._dataN).energies)[:,0]>eran[0]), sum(np.asarray(AllData(self._dataN).energies)[:,0]<eran[1])+1, len(np.asarray(AllData(self._dataN).energies)[:,0])))
        else:
            self._AllData(self._dataN).ignore("**-1.")
            self._AllData(self._dataN).ignore("10.-**")

        os.chdir("../")
        self._dofData+=len(self._AllData(self._dataN).energies)
        self._groupN.append(["XRT", self._dataN])
        self._dataN+=1


    def __prepBAT__(self, ti_num):

        # Import BAT dataset

        if self._dataN == 1:
            self._AllData('./BAT/{}-{}-{}.pha'.format(self._GRBname, 'BAT', ti_num))
        else:
            self._AllData('{}:{} ./BAT/{}-{}-{}.pha'.format(self._dataN, self._dataN,self._GRBname, 'BAT', ti_num))
        self._AllData(self._dataN).response = './BAT/{}-{}-{}.rsp'.format(self._GRBname, 'BAT', ti_num)
        
        if hasattr(self, 'eran'):
            eran = self.eran
        else:
            eran = [15.,150.]
        self._AllData(self._dataN).ignore("**-{:.1f}".format(eran[0]))
        self._AllData(self._dataN).ignore("{:.1f}-**".format(eran[1]))

        self._dofData+=len(self._AllData(self._dataN).energies)
        self._groupN.append(["BAT", self._dataN])
        self._dataN+=1
    
    def __igChannel__(self, ti_num, det):

        if hasattr(self, 'eran'):
            if det in ['b0', 'b1']:
                if self.eran[0]<250.0:
                    eran =  [250.0, self.eran[1]]
                else:
                    eran = self.eran
            elif det == 'LLE':
                eran = self.eran
                #eran = [30000.0, 100000.0]
            else:
                if self.eran[1]>900.0:
                    eran = [self.eran[0], 900.0]
                else:
                    eran = self.eran
        else:
            if det in ['b0', 'b1']:
                eran =  [260.0, 38000.0]
            elif det == 'LLE':
                eran = [30000.0, 100000.0]
            else:
                eran = [10.0, 900.0]

        self._AllData(self._dataN).ignore("**-{:.1f}".format(eran[0]))
        self._AllData(self._dataN).ignore("{:.1f}-**".format(eran[1]))
 
    
    def __rearrModel__(self, mod_num):

        # Input: ["Model_1", "Model_2"]
        # Output: "Model_1+Model_2"
        
        models = ''
        for mn in self._modelset[mod_num]:
            models=models+'+'+mn
        models=models[1:]
        return models
    
    def __SelectedDet__(self, only, BAT=False, XRT=False):
        if only == 'LATonly': Detectors=['LAT']
        elif only == 'GBMonly': Detectors=self._usedGBM
        elif only == 'LLEonly': Detectors=['LLE']    
        elif only == 'LATLLE': Detectors=['LLE','LAT']
        elif only == 'GBMLLE': Detectors=self._usedGBM + ['LLE']  
        elif only == 'GBMLAT': Detectors=self._usedGBM + ['LAT']  
        elif only == 'Swiftonly': Detectors=[]
        else: Detectors = self._usedGBM+['LLE','LAT'] 
        if BAT : Detectors = ['BAT'] + Detectors
        if XRT: Detectors = Detectors+['XRT']

        return Detectors

    def __Xspec__(self, mod_num, ti_num, newInit, debug=False, BAT=False, XRT=False, dtType=1, rspT=True):    
        from xspec import AllData, AllModels, Xset, Model, Fit, Plot
        
        AllData.clear()
        AllModels.clear()
        Xset.allowPrompting=False

        self.__reset__(newP = self._newP)
        self._AllData = AllData
        self._AllModels = AllModels
        self._dataN = 1
        self._groupN = []
        self._dofData = 0
        
        AllModels.addPyMod(XspecFunctions.bkn2power, XspecFunctions.bkn2pow_Info, 'add')
        AllModels.addPyMod(XspecFunctions.multiBB, XspecFunctions.multiBB_Info, 'add')
         
        self._usedData = self.__SelectedDet__(self._only)
        
        if self._only != 'Swiftonly': self.__prepGBMLLE__(ti_num, dtType=dtType, rspT = rspT)

        if 'LAT' in self._usedData: self.__prepLAT__(ti_num, dtType=dtType)
        
        if BAT: self.__prepBAT__(ti_num)

        if XRT: self.__prepXRT__(ti_num)


        self._usedData = self.__SelectedDet__(self._only, BAT=BAT, XRT=XRT)

        # Set modelname
        models = self.__rearrModel__(mod_num)

        if XRT: 
            models = models+'* TBabs * zTBabs'

        # Input the model to Xspec
        if newInit and len(self._newP)>0:
            changedInitP = np.asarray(self._newP)
            if mod_num in changedInitP[:,0]:
                paramSET = changedInitP[:,1][changedInitP[:,0] == mod_num][0]
                
            else:
                paramSET = self._setP[mod_num]
        else:
            paramSET = self._setP[mod_num]

        Mmanager=Model(models, setPars=paramSET) 

        self._groupN = np.asarray(self._groupN)

        if XRT:
            paramSET.update({Mmanager.nParameters: '{}, 0.'.format(self._redshift)})
            paramSET.update({Mmanager.nParameters-2: '{}, 0.'.format(self._nH)})
            Mmanager=Model(models, setPars=paramSET) 
        
        # Constant factor
        if hasattr(self, 'const'):
            if self.const:
                Mmanager=Model(models+'* const', setPars=paramSET)
                Mmanager.constant.factor.frozen=True
                #if XRT:
                #    if self._only != "Swiftonly":
                #        AllModels(np.int(self._groupN[self._groupN[:,0] == 'XRT'][0][1])).setPars({Mmanager.nParameters:'1.0, 0.1, 0.67, 0.67, 1.5, 1.5'})
                #    else:
                #        Mmanager.constant.factor.frozen=True
                if BAT:
                    AllModels(np.int(self._groupN[self._groupN[:,0] == 'BAT'][0][1])).setPars({Mmanager.nParameters:'1.0, 0.1, 0.67, 0.8, 1.25, 1.5'})
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
                Fit.statMethod = "chi {}".format(int(self._dataN-2))
                Fit.statMethod = "cstat {}".format(int(self._dataN-1))
            elif BAT:
                Fit.statMethod = "pgstat"
                Fit.statMethod = "chi {}".format(int(self._dataN-1))
            elif XRT:
                Fit.statMethod = "pgstat"
                Fit.statMethod = "cstat {}".format(int(self._dataN-1))
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

        x, y, m = self.__plotSetting__(debug = debug, BAT=BAT, XRT=XRT)
        
        self._Mmanager = Mmanager
        self._Fit = Fit

        n = self._dofData
        k = self._dofData-Fit.dof
    
        bic= Fit.statistic + k*np.log(n)

        Statistic = [Fit.statistic, Fit.dof, bic]

        output = {"Model":self.__FitResult__(), "Statistics":Statistic, "PlotX":x, "PlotY":y, "Covariance":Fit.covariance, "PlotModel": m, "TimeInterval":[self._tS_set[ti_num]-self._T0,self._tE_set[ti_num]-self._T0]}
        return output
    
    def __plotSetting__(self, debug = False, BAT=False, XRT=False):

        Plot.commands = ()
        Plot.xLog = True
        Plot.yLog = True
        Plot.device = "/xs"
        Plot.xAxis="keV"
        #Plot('data', "eemodel")
        
        for det, groupN in self._groupN:
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
    
    def __makeTable__(self, ti_num, clear=False, norm=True, latex=False, Epiv=100):
        tS = self._tS_set[ti_num]-self._T0
        tE = self._tE_set[ti_num]-self._T0
        multi = False

        if clear or not(hasattr(self, "_Ftable")):
            self.__clearTable__()
            newT = True
        else:
            newT = False
            

        for fresult in self._output["Model"]:
            self._FparN +=fresult['pars']
        self._FparN = list(set(self._FparN))
        self._FparN.sort()

        if not(norm): 
            self._FparN.remove('norm')

        for fresult in self._output["Model"]:
            for par in self._FparN:
                if not(self._Ftable.has_key(par)):
                    self._Ftable.update({par:[]})
                    if not(newT):
                        for i in range(len(self._Ftable['Model'])):
                            self._Ftable[par].append('')
                try:
                    #if par == 'norm' and fresult['name'] == 'powerlaw':
                    #    self._Ftable[par].append(r"{:.3f}$^{{+{:.3f}}}_{{-{:.3f}}}$".format(Ampli[0], Ampli[1]-Ampli[0], Ampli[0]-Ampli[2]))
                    if par == 'ep2' and fresult['name'] == 'grbm*highcut':
                        Ep2 = self.__foldEtoEp__()
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
                self._Ftable['Time'].append(["{:.2f}".format(tS), "{:.2f}".format(tE)])
                self._Ftable['Stats'].append(["{:.1f}".format(self._output["Statistics"][0]), self._output["Statistics"][1], "{:.1f}".format(self._output["Statistics"][2])])
                self._Ftable['T_num'].append(ti_num)

            multi = True

    def __clearTable__(self):
        self._Ftable = {}
        self._FparN = []
        self._Ftable.update({"Model":[]})
        self._Ftable.update({"Stats":[]})
        self._Ftable.update({"Time":[]})
        self._Ftable.update({"T_num":[]})

    def __tableList__(self, latex=False):
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

    def __calNorm__(self, Epiv):
        # Input: PL = N*E^(-gamma)
        # Output: PL = N*(E/100)^(-gamma)
        parSet = self.__genPars__(sim=10000)
        Comp = []
        for parS in parSet:
            parN = 0
            for fresult in self._output['Model']:
                model_name = fresult['name']
                if model_name == 'powerlaw':
                    nufnu = ModelCov[model_name](Epiv, *parS[ParNCov[model_name]+parN])
                    Comp.append(nufnu)
                parN += len(fresult['pars'])
        Comp = np.asarray(Comp)
        return [np.percentile(Comp, 50), np.percentile(Comp, 84),np.percentile(Comp, 16)]
    
    def __foldEtoEp__(self):
        # Input: E_f
        # Output: Ep = (alpha + 2)*E_f
        parSet = self.__genPars__(sim=10000)
        Ep2 = (parSet[:,1]+2.)*parSet[:,5]
        return [np.percentile(Ep2, 50), np.percentile(Ep2, 84),np.percentile(Ep2, 16)]


    def __makeSED__(self, mod_num, ti_num, noLC = False, saveFig = False, tRan = [], eRan = [50, 300], overlap=True, verbose=False):
        if not(noLC):
            fig = plt.figure(figsize=(5,6))
            gs = plt.GridSpec(2, 1, height_ratios=[3, 12])
            gs.update(hspace=0.32)
            topax = fig.add_subplot(gs[0])

            lcData = self.__refinedData__(ecut = eRan)[0]

            self.lightCurve = lcData

            binsDef = [-10]+(self._tS_set-self._T0).tolist()+[self._tE_set[-1]-self._T0]+[self._tE_set[-1]-self._T0+100]

            bk, err = XspecBackground.PolyFit(lcData, t1=self._T90[0], t2=self._T90[1])
            tb = 0.064
            if not(binsDef[1:]>=binsDef[:-1]):

                cnts, t = np.histogram(lcData, bins=binsDef)
                ts = np.asarray(t[1:]-t[:-1])
                t = centerPt(t)
                B = np.asarray(bk(t))
                topax.step(binsDef, [((cnts/ts-B/tb)/1000.)[0]]+((cnts/ts-B/tb)/1000.).tolist(), color='k')
                

            lcBins = np.arange(self._tS_set[0]-self._T0-10, self._tE_set[-1]-self._T0+100, step=0.064)
            cnts, t = np.histogram(lcData, bins=lcBins)
            ts = np.asarray(t[1:]-t[:-1])
            t = centerPt(t)
            B = np.asarray(bk(t))

            topax.step(lcBins, [((cnts/ts-B/tb)/1000.)[0]]+((cnts/ts-B/tb)/1000.).tolist(), color='k', alpha=0.5)

            topax.set_xlabel("Time since triggered [s]", fontsize=12)
            topax.set_ylabel(r"Counts/s", fontsize=12)
            topax.text(-0.03, 1.05, r"$\times$ 10$^3$", transform=topax.transAxes)
            topax.text(1, 1.05, "50 - 300 keV", transform=topax.transAxes, ha="right")

            if len(tRan) == 0:
                topax.set_xlim(self._T90[0]-5, self._T90[1]+10)
            else:
                topax.set_xlim(tRan[0], tRan[1])

            topax.axhline(0, ls=":", color='k')
            topax.axvspan(self._tS_set[ti_num]-self._T0,self._tE_set[ti_num]-self._T0, facecolor='yellow', edgecolor="k", alpha=0.2, lw=1, zorder=-1)
            topax.axvline(self._T90[0], color='red', lw=1)
            topax.axvline(self._T90[1], color='red', lw=1)
            gs_base = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], height_ratios=[9, 2], hspace=0.05)
            ax = [fig.add_subplot(gs_base[0,:]),fig.add_subplot(gs_base[1,:])]
        else:
            f, ax = plt.subplots(2, 1, figsize=(5,5), gridspec_kw = {'height_ratios':[7, 2]})
        x = self._output["PlotX"]
        y = self._output["PlotY"]
        m = self._output["PlotModel"]
        newArr = x[:,0].argsort()
        x = x[newArr]
        y = y[newArr]
        m = m[newArr]

        #try:
        modelF, model_compF, etc = self.__modelFlux__(x[:,0], verbose=verbose)
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
            modelF, model_compF, etc = self.__modelFlux__(E, verbose=verbose)
            ax[0].loglog(E, E**2*modelF, color='k')
            for modelF in model_compF:
                ax[0].loglog(E, E**2*modelF, color='k', ls=":", lw=0.7)
        except:
            ax[0].loglog(x[:,0], m, color='k')
        
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

        if saveFig:
            plt.savefig("{}-{}.pdf".format(ti_num,mod_num), bbox_inches="tight")
        
        plt.show(block=False)

    def __modelFlux__(self, E, erg=False, butterfly = False, verbose=False, XRT=False):
        parSet = self.__genPars__(sim=butterfly)

        synF = []
        
        if XRT:
            parSet = [[x, y] for x, y in zip(parSet[:,0], parSet[:,1])]
        
        for parS in parSet:
            CompF = []
            TotalF = np.zeros(len(E))
            parNum = 0
            for fresult in self._output['Model']:
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
                print("{:.1f} %".format(len(synF)*100./len(parSet)))
                clear_output(wait=True)
        #if verbose: 
        
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
                    temp[parN] = self._AllModels(self._dataN-1).constant.factor.values[0]
                    temp[parN+'_err'] = [self._AllModels(self._dataN-1).constant.factor.values[0]-abs(self._AllModels(self._dataN-1).constant.factor.sigma),self._AllModels(self._dataN-1).constant.factor.values[0]+abs(self._AllModels(self._dataN-1).constant.factor.sigma)]
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
        DATA = fits.open("./Xspec/{}-mktime-{}.fits".format(self._GRBname,fileN))[1].data
        maxE = DATA.field("Energy")
        if len(maxE) == 0:
            Chan = fits.open("./Xspec/{}-LLE-{}.PHA".format(self._GRBname,fileN))[2].data
            DATA = fits.open("./Xspec/{}-LLE-{}.PHA".format(self._GRBname,fileN))[1].data.field('Counts')[0]
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
            if np.size(self._output['Covariance']) == temp:
                break

        covReshape = np.zeros(shape=(matSize,matSize))

        k=0
        for j in range(matSize):
            for i in range(matSize):
                if i <= j and covReshape[i][j] == 0: 
                    covReshape[i][j]=self._output['Covariance'][k]
                    covReshape[j][i]=self._output['Covariance'][k]
                    k+=1
                    continue
                continue

        return covReshape

    def __genPars__(self, runs = 10000, sim=False):
        mean = []
        fixed_mean = []
        for pars in self._output['Model']:
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


    def __refinedData__(self, ecut=[50, 300]):
            
        lcData = []
        for i in range(np.size(ecut)-1):
            lcData.append([])

        for uG, i in zip(self._usedGBM, range(np.size(self._usedGBM))):

            with fits.open(glob.glob(self._address+"/GBM/*_"+uG+"_*.fit")[0]) as sourceFile:

                if 'n' in uG:
                    cut = (sourceFile[1].data.field('E_MIN')>10)*(sourceFile[1].data.field('E_MAX')<1000)
                elif 'b' in uG:
                    cut = (sourceFile[1].data.field('E_MIN')>200)*(sourceFile[1].data.field('E_MAX')<40000)

                ecut_temp=[]    
                for elim in ecut:
                    ecut_temp.append(128-sum(sourceFile[1].data.field('E_MIN')>elim))

                dataT = sourceFile[2].data.field('TIME')-self._T0
                dataChannel = sourceFile[2].data.field('PHA')

                for j in range(np.size(ecut_temp)-1): 
                    lcData_temp=np.asarray(dataT)*np.asarray(dataChannel>=ecut_temp[j])*np.asarray(dataChannel<=ecut_temp[j+1])*np.asarray(dataT>=-300)*np.asarray(dataT<=300)
                    lcData_temp = lcData_temp[lcData_temp != 0]
                    lcData[j]=lcData[j]+lcData_temp.tolist()

        for k in range(np.size(ecut)-1):
            lcData[k]=np.asarray(lcData[k])
            
        return np.asarray(lcData)

