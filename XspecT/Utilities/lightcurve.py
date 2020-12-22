import numpy as np
from XspecT.info import EventInfo

def LightCurve(grb, energy_cut=[50, 300], time_cut=[-10, 50]):
    GRBinfo = EventInfo(grb)

    lcData = []
    for i in range(np.size(energy_cut)-1):
        lcData.append([])

    for uG, i in zip(GRBinfo.usedGBM, range(np.size(GRBinfo.usedGBM))):

        with fits.open(glob.glob(GRBinfo._address+"/GBM/*_"+uG+"_*.fit")[0]) as sourceFile:

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