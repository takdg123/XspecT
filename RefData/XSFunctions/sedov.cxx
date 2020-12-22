//    Dec 98    Sedov model
//		   additive model for XSPEC to use with general models.
//		   parameters (param):
//              (1):    postshock temperature t(keV)
//              (2):    postshock electron temperature t_e(keV)
//		(3):    heavy metal abundance (with respect to solar)
//              (4):    ionization parameter (cm^-3 s)
//		(5):	   redshift z

#include <xsTypes.h>
#include <functionMap.h>

void sedov(const RealArray& energyArray, const RealArray& params,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
  RealArray vparam(17);
  vparam[0] = params[0];
  vparam[1] = params[1];
  vparam[2] = 1.0;
  vparam[3] = 1.0;
  for (size_t i=4; i<=14; ++i) vparam[i] = params[2];
  vparam[15] = params[3];
  vparam[16] = params[4];

  // Call the vsedov subroutine for variable abundances.

  vsedov(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);
  return;
}

