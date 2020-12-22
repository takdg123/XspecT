// Runs the APED interpolation with thermal and velocity broadening
//      parameters (params):
//           0:      kT temperature in keV
//           1:      kTi temperature in keV
//           2:      Metal abundances (H and He fixed at cosmic)
//           3:      redshift z
//           4:      gaussian velocity broadening
//           Norm = (4 * pi * 1e14)^-1 * Int(n^2)dV / D^2 where n is
//           in cm^-3 and D is the distance in cm.

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

void btapec(const RealArray& energyArray, const RealArray& params, 
	   int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	   const string& initString)
{

  // set up parameter array for bvtapec model and call it

  RealArray vparams(17);
  vparams[0] = params[0];
  vparams[1] = params[1];
  vparams[2] = 1.0;
  for (size_t i=2; i<15; i++) vparams[i] = params[2];
  vparams[15] = params[3];
  vparams[16] = params[4];

  bvtapec(energyArray, vparams, spectrumNumber, fluxArray, fluxErrArray, initString);

  return;
}
