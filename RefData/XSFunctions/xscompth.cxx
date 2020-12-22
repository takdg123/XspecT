//c    written by Paolo  Coppi 
//c    changed by Juri to increase the speed (Oct 28, 1997) 
//c    relativistic reflection is added Oct 29, 1997
//c    parameters changed by MG on Nov 5, 1997
//c    diskbb replaced by diskpn by MG, Jun 28, 1998
//c    fits data load routines added by PC, Jul 23, 2001
//
//      SUBROUTINE XSCOMPTH(Ear,Ne,Param,Ifl,Photar,Photer)
//
//      IMPLICIT NONE
//
//      INTEGER Ne , Ifl
//      REAL Param(20) , Ear(0:Ne) , Photar(Ne), Photer(Ne)
//
//c
//c     driver for the nonthermal/thermal pair program
//c
//c     number of model parameters: 18
//c     parameters are given below in the order of decreasing importance:
//c     1: hard-to-soft compactness ratio
//c     2: blackbody compactness (if <0, uncomptonized black body component only)
//c     3: blackbody temperature in eV (if <0, use diskpn, T_{max}=abs(Par(3)) )
//c     4: non-thermal-to-hard compactness ratio (0 for pure thermal plasma)
//c     5: Thomson optical depth of ionization electrons (e.g., 0)
//c     6: source region radius in cm (for Coulomb/bremsstrahlung only)
//c     7: minimum Lorentz factor for power law injection (not used for mono)
//c     8: maximum Lorentz factor for power law injection 
//c     9: power law injection index (if <0, inject monoenergetically at #8)
//c     10 choose primary pair or electron injection (0=electron, 1=pair)
//c     11: cosIncl, cosine of inclination angle
//c     12, R, amount of reflection (if < 0 then only reflection component)
//c     13: FeAb, iron abundance in units of solar
//c     14: abundance of elements heavier than He relative to solar
//c     15: disk temperature (used for reflection)
//c     16: xi, ionization parameter
//c     17: beta, reflection emissivity law (r^beta) IF=-10 then 
//c         non-rotating disk, IF=10 then  (1.-sqrt(6./rg))/rg**3
//c     18: radin, inner radius of the disk in M
//c     19: radout, outer radius of the disk in M
//c     20: redshift
//c
//c    Note: if #3 OR #12 < 0, only show BB and/or Reflected components
//c
//c     algorithm:
//c           a(x)=(non-thermal/thermal pair spectrum)+reflection
//c

#include <xsTypes.h>
#include <functionMap.h>

// declaration for main wrap routine for eqpair, eqtherm and compth

void eqpwrap(const RealArray& energyArray, const RealArray& params,
             RealArray& flux, RealArray& fluxerr, int modelNumber);



void xscompth(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
  eqpwrap(energyArray, params, flux, fluxErr, 3);
}


