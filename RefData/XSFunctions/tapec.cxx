// Runs the APED interpolation
//      parameters (params):
//           0:      kT temperature in keV
//           1:      kT ion temperature in keV
//           2:      Metal abundances (H and He fixed at cosmic)
//           3:      redshift z
//           Norm = (4 * pi * 1e14)^-1 * Int(n^2)dV / D^2 where n is
//           in cm^-3 and D is the distance in cm.

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

void tapec(const RealArray& energyArray, const RealArray& params, 
	   int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	   const string& initString)
{

  // set up parameter array for vtapec model and call it

  RealArray vparams(16);
  vparams[0] = params[0];
  vparams[1] = params[1];
  vparams[2] = 1.0;
  for (size_t i=3; i<15; i++) vparams[i] = params[2];
  vparams[15] = params[3];

  vtapec(energyArray, vparams, spectrumNumber, fluxArray, fluxErrArray, initString);

  return;
}
