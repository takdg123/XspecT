#include <xsTypes.h>
#include <functionMap.h>

void xneq(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   //
   //   Jan 97  Nonequlibrium ionization models with constant electron
   //           temperature
   //	        additive model for XSPEC to use with general models.
   //		parameters (param):
   //           (1):    temperature t(keV)
   //		(2):    heavy metal abundance (with respect to solar)
   //           (3):    ionization parameter (cm^-3 s^-1)
   //		(4):	redshift z
   //
   //

   RealArray vparam(16);
   vparam[0] = params[0];
   vparam[1] = 1.0;
   vparam[2] = 1.0;
   for (size_t i=3; i<14; ++i)
      vparam[i] = params[1];
   vparam[14] = params[2];
   vparam[15] = params[3];

   xsneq(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);
}
