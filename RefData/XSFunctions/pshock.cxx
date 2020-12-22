//    Aug 9,99  Nonequlibrium ionization models with constant electron
//              temperature, with a distribution of ionization timescales
//              appropriate for a plane parallel shock or its section.
//              additive model for XSPEC to use with general models.
//              parameters (param):
//              (1):    temperature t(keV)
//              (2):    heavy metal abundance (with respect to solar)
//                    ionization range bounded by:
//              (3):    ionization parameter taul(cm^-3 s)
//              (4):    ionization parameter tauu(cm^-3 s)
//              (5):    redshift z

#include <xsTypes.h>
#include <functionMap.h>

void pshock(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	    const string& initString)
{
  RealArray vparam(17);
  vparam[0] = params[0];
  vparam[1] = 1.0;
  vparam[2] = 1.0;
  for (size_t i=3; i<14; ++i) vparam[i] = params[1];
  vparam[14] = params[2];
  vparam[15] = params[3];
  vparam[16] = params[4];

  // Call the subroutine for variable abundances.

  vpshock(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);
  return;
}

