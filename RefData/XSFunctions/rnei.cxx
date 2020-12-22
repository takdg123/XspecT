//    Rnei model
//	Recombining plasma model
//		   parameters (param):
//              (0):    current temperature t(keV)
//              (1):    initial temperature t_i(keV)
//		(2):    heavy metal abundance (with respect to solar)
//              (3):    ionization parameter (cm^-3 s)
//		(4):	redshift z

#include <xsTypes.h>
#include <functionMap.h>

void rnei(const RealArray& energyArray, const RealArray& params,
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

  // Call the vrnei subroutine for variable abundances.

  vrnei(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);
  return;
}

