import numpy as np
from numba import vectorize, njit, float64
from scipy.special import kv
from scipy.interpolate import interp1d

# cgs
me = 9.1e-28
e = 4.8e-10
c = 3e10
mec2 = me * c**2
h = 6.62e-27
hbar = h/2/np.pi
alphaF = 1/137
sigma_T = 6.65245854533e-25
Bc = me**2 * c**3 / e / hbar

eV2erg = 1.60218e-12
erg2eV = 1/eV2erg


def Eph_syn(B, Eel):
    return B/Bc * Eel**2/mec2


def Eph_IC_max(Eel, Eph):
    b = 4 * Eel * Eph / mec2**2
    return Eel * b/(1+b)


def dN_dE_el(Eel, Eelmin, Eelb, Eelmax, p, s, chimin, chimax, dp=1):
    return (
        (Eel/Eelmin)**(-p) * (1 + (Eel/Eelb)**(dp*s))**(-1/s) *
        np.exp(-(Eel/Eelmax)**chimax) * np.exp(-(Eelmin/Eel)**chimin)
    )


def Kernel_Syn(E_electron, E_photon, B):
    '''pitch angle averaged spectral synchrotron term'''
    xc = mec2 * E_photon/(3 * B/Bc * E_electron**2)
    constants = 2*np.sqrt(3)/9 * alphaF
    norm = mec2**4 / h
    return constants * norm * Bc/B * E_photon/E_electron**4 * (
        (kv(4/3, xc)*kv(1/3, xc) - 3/5*xc*(kv(4/3, xc)**2 - kv(1/3, xc)**2))
    )

# full inverse Compton cross section


@vectorize([float64(float64, float64, float64)])
def s_diff_full(e, z, g):
    '''
    e = incoming photon
    z = outgoing photon
    g = electron
    all normalised to mec2
    '''
    L = 2.0 * g * e
    K = z / g
    K_ = 1.0 - K
    R = K / L

    a1 = (R - 2.0 + 2.0 * K) * (R + 1.0 + 0.5 * K * K - K)
    a2 = -2.0 * K_ * R * np.log(0.5 * R / K_)
    a3 = -0.75 / g / L / K_ / K_

    Kmax = 2.0 * L / (1.0 + 2.0 * L)

    if (K < Kmax):
        return a3 * (a1 + a2)
    else:
        return 0.0

# approx for highest bins


@njit
def s_diff_KN_L(e, g, D_X):
    a1 = 0.375 / D_X / g / g / e
    a2 = 1.0 / (1.0 - D_X)
    return a1 * a2


@njit
def s_diff_KN_R(e, g, D_X):
    a1 = 0.375 / D_X / g / g / e
    a2 = np.log(D_X * g * e) + 2.0 * np.log(2.0) - 2.0
    return a1 * a2


@vectorize([float64(float64, float64, float64, float64)])
def Kernel_IC(Eph, Ephout, Eel, D_X):
    # units
    Eph *= 1/mec2
    Ephout *= 1/mec2
    Eel *= 1/mec2
    # factor for grid for averaging
    Nm = 100
    m = np.exp(np.linspace(-D_X/2, D_X/2, Nm))
    conv = c * sigma_T / mec2  # unit conversion
    # determine regime
    L = 2 * Eel * Eph
    if Eph <= Ephout:
        if L < 10:
            # full expression
            return s_diff_full(Eph, Ephout, Eel)*conv
        elif L < 1000:
            if np.log(Ephout) < np.log(Eel) - 5 * D_X:
                return s_diff_full(Eph, Ephout, Eel)*conv
            else:
                # average
                return np.sum(s_diff_full(Eph, Ephout*m, Eel))/Nm*conv
        else:
            if np.log(Eel)-D_X > np.log(Ephout) + D_X/2:
                return s_diff_full(Eph, Ephout, Eel)*conv
            if (np.log(Ephout) - D_X/2 < np.log(Eel)-D_X) and (np.log(Eel)-D_X <= np.log(Ephout) + D_X/2):
                # bin before electron energy
                return conv * s_diff_KN_L(Eph, Eel, D_X)
            if (np.log(Ephout) - D_X/2 < np.log(Eel)) and (np.log(Eel) <= np.log(Ephout) + D_X/2):
                # bin with electron energy
                return conv * s_diff_KN_R(Eph, Eel, D_X)
            else:
                return 0
    else:
        return 0.


class SSCfast:

    def __init__(self, Eel_min=1e6 * eV2erg, Eel_max=1e14 * eV2erg,
                 Eph_min=1e-3 * eV2erg, Eph_max=1e14 * eV2erg,
                 dlgEel=0.1, dlgEph=0.1, z=0):
        # grid specifications
        self.Eel_min = Eel_min
        self.Eel_max = Eel_max
        self.Eph_min = Eph_min
        self.Eph_max = Eph_max
        self.dlgEel = dlgEel
        self.dlgEph = dlgEph
        # build grids
        self.Eelgrid = np.exp(np.arange(
            np.log(Eel_min),
            np.log(Eel_max) + dlgEel,
            step=dlgEel
        ))
        self.Ephgrid = np.exp(np.arange(
            np.log(Eph_min),
            np.log(Eph_max) + dlgEph,
            step=dlgEph
        ))
        self.z = z
        self.E0obs = 1e3 * eV2erg 
        print("number of bins:")
        print("N_el:", len(self.Eelgrid))
        print("N_ph:", len(self.Ephgrid))
        print("E0_obs (norm):", self.E0obs*erg2eV, "eV")
        # precalculate IC kernel
        self.calcICKernel()

    def calcICKernel(self):
        self.KIC = Kernel_IC(
            self.Ephgrid[None, :, None], self.Ephgrid[:, None, None],
            self.Eelgrid[None, None, :], self.dlgEph
        )

    def evaluate(
            self, Eobs, Eelmin, Eelb, Eelmax, p, s, chimin, chimax, B, N_IC,
            Norm, Gamma):
        
        # doppler factor
        beta = np.sqrt(1-1/Gamma**2)
        doppler = 1 / Gamma / (1 - beta)
        red = (1+self.z)
        
        # calc electron spectrum
        self.dN_dE_el = dN_dE_el(
            self.Eelgrid, Eelmin, Eelb, Eelmax, p, s, chimin, chimax, dp=1)
        # calc synchrotron spectrum
        self.dN_dE_ph = np.sum(
            self.dlgEel*self.Eelgrid[:, None] * self.dN_dE_el[:, None] *
            Kernel_Syn(
                self.Eelgrid[:, None], self.Ephgrid[None, :],
                B),
            axis=0
        )
        # calc IC spectrum
        integrand = (
            self.dlgEel * self.dlgEph *
            self.Eelgrid[None, None, :] * self.dN_dE_el[None, None, :] *
            self.Ephgrid[None, :, None] * self.dN_dE_ph[None, :, None] *
            self.KIC
        )
        dN_dE_IC = np.sum(
            np.sum(
                integrand,
                axis=1
            ), axis=1
        )
        # constant factors of eV2erg and doppler cancel out
        dN_dE_IC *= np.max(self.Ephgrid**2 * self.dN_dE_ph)/np.max(self.Ephgrid**2 * dN_dE_IC)
        self.dN_dE_ph += N_IC * dN_dE_IC
        
        
        # interpolate to grid
        normFac = self.E0obs**2 * np.exp(
            np.interp(
                np.log(self.E0obs/doppler*red),
                np.log(self.Ephgrid),
                np.log(self.dN_dE_ph),
                left=-1e100, right=-1e100)
        )
        return Norm/normFac * np.exp(
            np.interp(
                np.log(Eobs/doppler*red),
                np.log(self.Ephgrid),
                np.log(self.dN_dE_ph),
                left=-1e100, right=-1e100))
