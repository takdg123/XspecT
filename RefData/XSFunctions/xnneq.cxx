#include <xsTypes.h>
#include <functionMap.h>


void xnneq(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   //    Feb 98    Nonequlibrium ionization models with two electron
   //              temperatures (final and ionization-timescale averaged)
   //		   additive model for XSPEC to use with general models.
   //		   parameters (param):
   //              (1):    temperature t(keV)
   //		   (2):    heavy metal abundance (with respect to solar)
   //              (3):    ionization parameter (cm^-3 s^1)
   //              (4):    tau-averaged temperature tav(keV)
   //		   (5):	   redshift z

   // Translated from xnneq.f 03/09

   RealArray vparam(17);
   vparam[0] = params[0];
   vparam[1] = 1.0;
   vparam[2] = 1.0;
   for (size_t i=3; i<14; ++i)
      vparam[i] = params[1];
   vparam[14] = params[2];
   vparam[15] = params[3];
   vparam[16] = params[4];

   // Call the xsneq subroutine for variable abundances.

   xsnneq(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);
}
