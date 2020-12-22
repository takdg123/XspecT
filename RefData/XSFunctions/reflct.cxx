#include "MZCompRefl.h"

// Main model routine

void reflct(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{

   //     Driver for angle-dependent reflection from a neutral medium.
   //     See Magdziarz & Zdziarski, 1995, MNRAS.

   //     The output spectrum is the sum of the input (the original contents
   //     of Photar) and the reflection component.
   //     The reflection component alone can be obtained
   //     for scale (see below) = rel_refl < 0. Then the actual
   //     reflection normalization is |scale|. Note that you need to
   //     change then the limits of rel_refl. The range of rel_refl in that case
   //     should exclude zero (as then the direct component appears).

   //     Version with variable iron abundance and new opacities of Balucinska
   //     & McCammon (1992, and 1994, private communication). As expected in AGNs,
   //     H and He are assumed to be fully ionized.

   //     This version allows for changes of the vector 'ear' between subsequent
   //       calls.


   //     number of model parameters: 5
   //     1: scale, scaling factor for reflection; if <0, no direct component
   //        (scale=1 for isotropic source above disk)
   //     2: redshift, z
   //     3: abundance of elements heavier than He relative to
   //        the solar abundances
   //     4: iron abundance relative to the solar iron abundance
   //     5: cosine of inclination angle
   //     algorithm:
   //          a[E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
   //     Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
   //     of the cutoff power law only (without reflection)
   //     and in the earth frame.

   // also uses the following values which can be set using xset
   //   REFLECT_MAX_E      maximum energy for which to calculate the output spectrum
   //   REFLECT_PRECISION  fractional precision for Greens' fn. adaptive integration


  RealArray ireflPars(7);
  for (size_t i=0; i<5; i++) ireflPars[i] = params[i];
  ireflPars[5] = 0.0;
  ireflPars[6] = 0.0;

  string modelName = "REFLECT";
  if (initString.find_first_not_of(" ") != string::npos) modelName = initString;

  ireflct(energyArray, ireflPars, spectrumNumber, flux, fluxErr, modelName);

}
