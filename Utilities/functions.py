import numpy as np
from numpy import *

POWERLAW = lambda E, gamma, N0: N0*(E/100.)**(gamma)
CUTOFFPL = lambda E, gamma, Ep, N0: N0*(E/100.)**(gamma)*np.exp(-E*(2.+gamma)/Ep)
STEP = CUTOFFPL
BB = lambda E, T, N0: N0*E**2/(np.exp(E/T)-1)/1e8

ePOWERLAW = lambda E, gamma, N0: N0*(E/100.)**(gamma)*E
eCUTOFFPL = lambda E, gamma, Ep, N0: N0*(E/100.)**(gamma)*np.exp(-E*(2.+gamma)/Ep)*E
eSTEP = eCUTOFFPL

def BBODY(E, T, N0):
    if np.size(E)!=1:
        BF = []
        for eng in E:
            if eng/T < 200: 
                BF.append(BB(eng, N0, T))
            else: 
                BF.append(10**-20*eng) 
    else:
        if E/T < 200: 
            BF = BB(E, N0, T)
        else: 
            BF = 10**-20*E
    return BF

def eBBODY(E, T, N0):
    if np.size(E)!=1:
        BF = []
        for eng in E:
            if eng/T < 200: 
                BF.append(BB(eng, N0, T)*eng)
            else: 
                BF.append(10**-20*eng*eng) 
    else:
        if E/T < 200: 
            BF = BB(E, N0, T)*E
        else: 
            BF = 10**-20*E*E
    return BF

def GRBM(E, alpha, beta, Ep, N0):
    BAND1 = lambda E, N0, alpha, beta, Ep: N0*(E/100.)**alpha*np.exp(-E*(2.+alpha)/Ep)
    BAND2 = lambda E, N0, alpha, beta, Ep: N0*(E/100.)**beta*((Ep/100.)*(alpha-beta)/(2+alpha))**(alpha-beta)*np.exp(-(alpha-beta))

    cvt = 2.+alpha
    if cvt <= 0: cvt = 0.0001
    amb = (alpha-beta)
    
    eb = Ep/cvt*amb
    if amb < 1./eb: amb = 1./eb
    
    ec = Ep/cvt

    if np.size(E) == 1:
        if E < eb:
            val = BAND1(E, N0, alpha, beta, Ep)
        else:
            val = BAND2(E, N0, alpha, beta, Ep)
    else:
        engs = np.asarray(E)
        
        if eb <= 0:
            val1 = N0*(engs/100.)**alpha*np.exp(-E/ec)
            val2 = np.asarray([])
        else:
            val1 = N0*(engs[engs<eb]/100.)**alpha*np.exp(-engs[engs<eb]/ec)
            val2 = N0*(engs[engs>=eb]/100.)**beta*((ec/100.)*amb)**(amb)*np.exp(-(amb))
        val = val1.tolist() + val2.tolist()
    
    return np.asarray(val)

def eGRBM(E, alpha, beta, Ep, N0):
    BAND1 = lambda E, N0, alpha, beta, Ep: N0*(E/100.)**alpha*np.exp(-E*(2.+alpha)/Ep)
    BAND2 = lambda E, N0, alpha, beta, Ep: N0*(E/100.)**beta*((Ep/100.)*(alpha-beta)/(2+alpha))**(alpha-beta)*np.exp(-(alpha-beta))

    if Ep<0: Ep = 10

    cvt = 2.+alpha
    if cvt <= 0: cvt = 0.0001
    amb = (alpha-beta)

    eb = Ep/cvt*amb

    if amb < 1./eb: amb = 1./eb
    
    ec = Ep/cvt

    if np.size(E) == 1:
        if eb<0:
            val = E*N0*(E/100.)**beta*((ec/100.)*amb)**(amb)*np.exp(-amb)
        elif E < eb:
            val = E*N0*(E/100.)**alpha*np.exp(-E/ec)
        else:
            val = E*N0*(E/100.)**beta*((ec/100.)*amb)**(amb)*np.exp(-amb)
    else:
        engs = np.asarray(E)
        vals = []
        
        if eb <= 0:
            val1 = engs*N0*(engs/100.)**alpha*np.exp(-engs/ec)
            val2 = np.asarray([])
        else:
            val1 = engs*N0*(engs[engs<eb]/100.)**alpha*np.exp(-engs[engs<eb]/ec)
            val2 = engs*N0*(engs[engs>=eb]/100.)**beta*((ec/100.)*amb)**(amb)*np.exp(-amb)
        val = val1.tolist() + val2.tolist()
    
    return np.asarray(val)

def BKNPOWER(E, alpha, Eb, beta, N0):
    if np.size(E)==1:
        if E < Eb:
            val = POWERLAW(E, alpha, N0)
        else:
            val = POWERLAW(E, beta, N0)*Eb**(-beta+alpha)
    else:
        engs = np.asarray(E)
        cutoff = -1
        for i in range(len(engs)):
            if engs[i] >= Eb: 
                cutoff = i
                break
        if min(E)>Eb:
            val1 = POWERLAW(engs, beta, N0)*Eb**(-beta+alpha)
            val2 = np.asarray([])
        elif max(E)<Eb:
            val1 = np.asarray([])
            val2 = POWERLAW(engs, alpha, N0)
        else:
            val1 = POWERLAW(engs[:cutoff], alpha, N0)
            val2 = POWERLAW(engs[cutoff:], beta, N0)*Eb**(-beta+alpha)
        val = val1.tolist() + val2.tolist()
        
    return np.asarray(val)
    
def GRBM_HIGHECUT(E, alpha, beta, Ep, N0, Ecut, Ef):
    BAND1 = lambda E, N0, alpha, beta, Ep: N0*(E/100.)**alpha*np.exp(-E*(2.+alpha)/Ep)
    BAND2 = lambda E, N0, alpha, beta, Ep: N0*(E/100.)**beta*((Ep/100.)*(alpha-beta)/(2+alpha))**(alpha-beta)*np.exp(-(alpha-beta))
    BAND3 = lambda E, N0, alpha, beta, Ep, Ecut, Ef: (N0*(E/100.)**beta*((Ep/100.)*(alpha-beta)/(2.+alpha))**(alpha-beta)*np.exp(-(alpha-beta)))*np.exp((Ecut-E)/Ef)

    engs = np.asarray(E)
    cvt = 2.+alpha
    if cvt <= 0: cvt = 0.0001
    amb = (alpha-beta)
    
    eb = Ep/cvt*amb
    if amb < 1./eb: amb = 1./eb
    
    ec = Ep/cvt

    cutoff = -1
    for i in range(len(engs)):
        if engs[i] >= eb: 
            cutoff = i
            break
    highcutoff = -1
    for j in range(len(engs)):
        if engs[j] >= Ecut and Ecut>=0 and Ef>=0 and Ecut>=Ep: 
            highcutoff = j
            break
    if cutoff == -1:
        val1 = BAND1(engs, N0, alpha, beta, Ep)
        val = val1.tolist()
    elif cutoff != -1 and highcutoff == -1:
        val1 = BAND1(engs[:cutoff], N0, alpha, beta, Ep)
        val2 = BAND2(engs[cutoff:], N0, alpha, beta, Ep)
        val = val1.tolist() + val2.tolist()
    else:
        val1 = BAND1(engs[:cutoff], N0, alpha, beta, Ep)
        val2 = BAND2(engs[cutoff:highcutoff], N0, alpha, beta, Ep)
        val3 = BAND3(engs[highcutoff:], N0, alpha, beta, Ep, Ecut, Ef)
        val = val1.tolist() + val2.tolist() + val3.tolist()
        
    return np.asarray(val)

def BKNPOWER_HIGHECUT(E, alpha, Eb, beta, N0, Ecut, Ef):
    engs = np.asarray(E)
    
    highcutoff = -1
    for j in range(len(engs)):
        if engs[j] >= Ecut and Ecut>=0 and Ef>=0 and Ecut>=Ep: 
            highcutoff = j
            break
    
    if highcutoff == -1:
        val = BKNPL(E, N0, alpha, beta, Eb)    
    else:
        val1 = BKNPL(engs[:highcutoff], N0, alpha, beta, Eb)
        val2 = BKNPL(engs[:highcutoff], N0, alpha, beta, Eb)
        val = val1.tolist() + val2.tolist()
        
    return np.asarray(val)


def bkn2power(engs, params, flux):
    E = np.asarray(engs)
    cutoff = -1
    idx1 = params[0]
    idx2 = params[1]
    idx3 = params[2]
    Eb1 = params[3]
    Eb2 = params[4]
    N0 = 1
    norm1 = Eb1**(-idx2+idx1)
    norm2 = norm1*Eb2**(-idx3+idx2)
    for i in range(len(E)-1):
        if E[i] >= Eb2: 
            val = N0*(E[i+1]**(idx3+1.)-E[i]**(idx3+1.))/(idx3+1.)
            if np.isnan(val):
                flux[i] = 1e-20
            else:
                flux[i] = val*norm2
        elif E[i] >= Eb1:
            val = N0*(E[i+1]**(idx2+1.)-E[i]**(idx2+1.))/(idx2+1.)
            if np.isnan(val):
                flux[i] = 1e-20
            else:
                flux[i] = val*norm1
        else:
            val = N0*(E[i+1]**(idx1+1.)-E[i]**(idx1+1.))/(idx1+1.)
            if np.isnan(val):
                flux[i] = 1e-20
            else:
                flux[i] = val
                

def BKN2POWER(E, idx1, idx2, idx3, Eb1, Eb2, N0):
    engs = np.asarray(E)
    vals = []
    cutoff = -1
    for i in range(len(engs)):
        if engs[i] >= Eb1: 
            cutoff1 = i
            break
    for i in range(len(engs)):
        if engs[i] >= Eb2: 
            cutoff2 = i
            break
    if min(E)>Eb1:
        val1 = np.asarray([])
        val2 = BKNPL(engs, N0, idx2, idx3, Eb2)
        val3 = np.asarray([])
    elif max(E)<Eb2:
        val1 = BKNPL(engs, N0, idx1, idx2, Eb1)
        val2 = np.asarray([])
        val3 = np.asarray([])
    else:
        val1 = POWERLAW(engs[:cutoff1], idx1, N0)
        val2 = POWERLAW(engs[cutoff1:cutoff2], idx2, N0)*Eb1**(-idx2+idx1)
        val3 = POWERLAW(engs[cutoff2:], idx3, N0)*Eb1**(-idx2+idx1)*Eb2**(-idx3+idx2)
    val = val1.tolist() + val2.tolist() + val3.tolist()
    return np.asarray(val)


def multiBB(engs, params, flux):
    BBmodel = lambda x, m: x**(2.-m)/(np.exp(x)-1.)
    E = np.asarray(engs)
    m = params[0]
    kT_min = params[1]
    kT_max = params[2]
    Tratio = ((kT_max/kT_min)**(m+1.)-1.)
    IE = np.zeros(len(E))
    
    for i in range(len(E)):
        IE[i] = quad(BBmodel, E[i]/kT_max, E[i]/kT_min, args=(m))[0]

    N_edp = (E/kT_min)**(m-1.)
    N = 8.0525 * (m + 1.)/Tratio * (kT_min)**(-2.)
    fdensity = N*N_edp*IE
    for i in range(len(E)-1):
        val = (fdensity[i+1]+fdensity[i])*(E[i+1]-E[i])/2.
        if np.isnan(val):
            flux[i] = 1e-20
        else:
            flux[i] = val
            
def MULTIBB(E, m, kT_min, kT_max, N0):
    BBmodel = lambda x, m: x**(2.-m)/(np.exp(x)-1.)
    E = np.asarray(E)
    
    Tratio = ((kT_max/kT_min)**(m+1.)-1.)

    N = N0*8.0525 * (m + 1.)/Tratio * (kT_min)**(-2.)
    N_edp = (E/kT_min)**(m-1.)
    IE = np.zeros(len(E))
    
    for i in range(len(E)):
        IE[i] = quad(BBmodel, E[i]/kT_max, E[i]/kT_min, args=(m))[0]

    return N*N_edp*IE



bkn2pow_Info = ("PhoIndx1   \"\"  -0.7  -10.0  -9.0  9.0  10.0  0.01",
                  "PhoIndx2   \"\"  -1.5  -10.0  -9.0  9.0  10.0  0.01",
                  "PhoIndx3   \"\"  -2.2  -10.0  -9.0  9.0  10.0  0.01",
                  "BreakE1  keV  200 10  10  10000. 200000. 1" ,
                  "BreakE2  keV  2000 10  10  10000000. 20000000. 1")

multiBB_Info = ("m   \"\"  0.0  -2.0  -2.0  2.0  2.0  0.01",
              "kT_min  keV  50.  1.0  1.0  1000.0  1000.0  0.1",
              "kT_max  keV  300.  1.0  1.0  2000.0  2000.0  0.1" )
