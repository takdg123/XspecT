#include "MZCompRefl.h"

// Prototype for the routine which actually does all the work (included in pexriv.cxx)

void doIonizedReflection(string ModelName, const RealArray& energyArray, 
			 const RealArray& params, RealArray& flux);


// Main model routine

void xsbexrav(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{

  //  Driver for angle-dependent reflection from an exponentially-cutoff
  //  BROKEN power law and neutral medium.

  //  number of parameters: 9
  //  0: Gamma1, first power law photon index
  //  1: E_break, break energy
  //  2: Gamma2, second power law photon index
  //  3: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
  //     one needs to change the lower limit for that)
  //  4: scaling factor for reflection (1 for isotropic source above disk) 
  //  5: cosine of inclination angle 
  //  6: abundance of elements heavier than He relative to
  //     the solar abundances
  //  7: iron abundance relative to the solar iron abundance (7.67)
  //  8: redshift 

  // also uses the following values which can be set using xset
  //   BEXRAV_PRECISION  fractional precision for Greens' fn. adaptive integration

   static bool isFirst = true;
   if (isFirst)
   {
      string msg("Compton reflection from ionized medium.");
      msg += "\nSee help for details.\nIf you use results of this model in a paper,";
      msg += "\nplease refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837";
      xs_write(const_cast<char*>(msg.c_str()), 5);
      isFirst = false;
   }

   // Create standard passed parameter array - for some bizarre reason the
   // parameters in bexrav are in a different order from bexriv

   RealArray passedParams(11);
   for (size_t i=0; i<5; i++) passedParams[i] = params[i];
   passedParams[5] = params[8];
   passedParams[6] = params[6];
   passedParams[7] = params[7];   
   passedParams[8] = params[5];
   passedParams[9] = 0.0;
   passedParams[10] = 0.0;

   doIonizedReflection("BEXRAV", energyArray, passedParams, flux);

}

