import numpy as np
import matplotlib.pyplot as plt

def center_pt(t):
    return (t[1:]+t[:-1])/2.

def linear(T, alpha, T0):
    return alpha*T+T0
    
def gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.np.exp(-(x-mu)**2/(2.*sigma**2))

def text_to_npy(fname):
    textfile = open("./{}.txt".format(fname), "r")
    fdata=[]
    data = textfile.readlines()
    for line in data:
        pars = line.split()
        ffdata=[]
        for par in pars:
            try: ffdata.append(float(par))
            except: ffdata.append(par)
        if np.size(ffdata)==1: fdata.append(ffdata[0])
        else: fdata.append(ffdata)
    textfile.close()
    np.save(fname,np.asarray(fdata))
