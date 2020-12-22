#include <xsTypes.h>
#include <functionMap.h>

void nei(const RealArray& energyArray, const RealArray& params,
	 int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	 const string& initString)
{
   //
   //   Jan 97  Nonequlibrium ionization models with constant electron
   //           temperature
   //	        additive model for XSPEC to use with general models.
   //		parameters (param):
   //           (0):    temperature t(keV)
   //		(1):    heavy metal abundance (with respect to solar)
   //           (2):    ionization parameter (cm^-3 s^-1)
   //		(3):	redshift z
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

   vnei(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);
}
