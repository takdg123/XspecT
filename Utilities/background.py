import numpy as np

def PolyFit(data, tbin_size=0.064, t1=0, t2=100):
    bk1 = np.histogram(data[(data<t1) + (data>t2)], bins=np.arange(min(data), max(data), step=tbin_size))
    bk2 = np.histogram(data[(data<t1) + (data>t2)], bins=np.arange(min(data), max(data), step=tbin_size))
    bk = np.asarray((bk1[0].tolist()+bk2[0].tolist()))
    bk_t = np.asarray(((bk1[1][1:]+bk1[1][:-1])/2.).tolist()+((bk2[1][1:]+bk2[1][:-1])/2.).tolist())
    bk_t = bk_t[bk>0][1:-2]
    bk = bk[bk>0][1:-2]

    chi_min=1e8
    for i in range(5):
        p_bk = np.polyfit(bk_t[bk>0], bk[bk>0], i)
        bk_temp=np.poly1d(p_bk)
        chi = sum((bk-bk_temp(bk_t))**2./abs(bk))
        if (chi_min - chi) > 9:
            chi_min=chi
            bk_min=bk_temp
        else:
            break
    return bk_min, np.average(bk-bk_min(bk_t))
