#include "MZCompRefl.h"

// Prototype for the routine which actually does all the work (included in pexriv.cxx)

void doIonizedReflection(string ModelName, const RealArray& energyArray, 
			 const RealArray& params, RealArray& flux);


// Main model routine

void xsbexriv(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{
   //  Driver for angle-dependent reflection
   //
   //  See Magdziarz & Zdziarski, 1995, MNRAS.
   //  See Zdziarski et al., 1995, ApJL 438, L63 for description of
   //  calculation of ionization (based on Done et al. 1992, ApJ 395, 275).
   //  The abundances are defined by the command abund
  //
  //     number of model parameters:11
  //     0: Gamma1, first power law photon index
  //     1: E_break, break energy
  //     2: Gamma2, second power law photon index
  //     3: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
  //        one needs to change the lower limit for that)
  //     4: scale, scaling factor for reflection; if <0, no direct component
  //        (scale=1 for isotropic source above disk)
  //     5: redshift, z
  //     6: abundance of elements heavier than He relative to
  //        the solar abundances
  //     7: iron abundance relative to the solar abundances
  //     8: cosine of inclination angle
  //     9: disk temperature in K
  //     10: disk ionization parameter = L/nR^2
  //     algorithm:
  //          a[E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
  //     Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
  //     of the cutoff power law only (without reflection)
  //     and in the earth frame.
  //
  //
  // also uses the following values which can be set using xset
  //   BEXRIV_PRECISION  fractional precision for Greens' fn. adaptive integration
  //  Based on ireflct.cxx

   static bool isFirst = true;
   if (isFirst)
   {
      string msg("Compton reflection from ionized medium.");
      msg += "\nSee help for details.\nIf you use results of this model in a paper,";
      msg += "\nplease refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837";
      xs_write(const_cast<char*>(msg.c_str()), 5);
      isFirst = false;
   }

   doIonizedReflection("BEXRIV", energyArray, params, flux);

}

