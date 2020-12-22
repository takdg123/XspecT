#include "MZCompRefl.h"

// Prototype for the routine which actually does all the work and routine which is
// required for the broken power-law case (bexrav model).

void doIonizedReflection(string ModelName, const RealArray& energyArray, 
			 const RealArray& params, RealArray& flux);

// Main model routine

void xspexrav(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{

   //    Driver for angle-dependent reflection from an exponentially-cutoff
   //    power law and neutral medium.
   //    See Magdziarz & Zdziarski, 1995, MNRAS.

   //    The output spectrum is the sum of the cutoff power law and the
   //    reflection component.
   //    The reflection component alone can be obtained
   //    for scale (see below) = rel_refl < 0. Then the actual
   //    reflection normalization is |scale|. Note that you need to
   //    change then the limits of rel_refl. The range of rel_refl in that case
   //    should exclude zero (as then the direct component appears).

   //    If E_c=0, there is no cutoff in the power law.

   //    Version with variable iron abundance and new opacities of Balucinska
   //    & McCammon (1992, and 1994, private communication). As expected in AGNs,
   //    H and He are assumed to be fully ionized.

   //    This version allows for changes of the vector 'ear' between subsequent
   //      calls.


   //    number of model parameters: 7
   //    1: Gamma, power law photon index, N_E prop. to E^{-Gamma}
   //    2: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
   //       one needs to change the lower limit for that)
   //    3: scale, scaling factor for reflection; if <0, no direct component
   //       (scale=1 for isotropic source above disk)
   //    4: redshift, z
   //    5: abundance of elements heavier than He relative to
   //       the solar abundances
   //    6: iron abundance relative to the solar iron abundance
   //    7: cosine of inclination angle
   //    algorithm:
   //         a[E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
   //    Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
   //    of the cutoff power law only (without reflection)
   //    and in the earth frame.
   //
   // also uses the following values which can be set using xset
   //   PEXRAV_PRECISION  fractional precision for Greens' fn. adaptive integration

   static bool isFirst = true;
   if (isFirst)
   {
      string msg("Compton reflection from neutral medium.");
      msg += "\nSee help for details.\nIf you use results of this model in a paper,";
      msg += "\nplease refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837";
      xs_write(const_cast<char*>(msg.c_str()), 5);
      isFirst = false;
   }

   // add parameter entries for Temp and Xi set to zero which will force the
   // use of NeutralOpacity instead of IonizedOpacity.

   RealArray passedParams(9);
   for (size_t i=0; i<7; i++) passedParams[i] = params[i];
   passedParams[7] = 0.0;
   passedParams[8] = 0.0;

   doIonizedReflection("PEXRAV", energyArray, passedParams, flux);

}
