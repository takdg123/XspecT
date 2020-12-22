// Runs the APED interpolation with thermal and velocity broadening
//      parameters (params):
//           0:      kT temperature in keV
//           1:      Metal abundances (H and He fixed at cosmic)
//           2:      redshift z
//           3:      gaussian velocity broadening
//           Norm = (4 * pi * 1e14)^-1 * Int(n^2)dV / D^2 where n is
//           in cm^-3 and D is the distance in cm.

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

void bapec(const RealArray& energyArray, const RealArray& params, 
	   int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	   const string& initString)
{

  // set up parameter array for vbapec model and call it

  RealArray vparams(16);
  vparams[0] = params[0];
  vparams[1] = 1.0;
  for (size_t i=2; i<14; i++) vparams[i] = params[1];
  vparams[14] = params[2];
  vparams[15] = params[3];

  bvapec(energyArray, vparams, spectrumNumber, fluxArray, fluxErrArray, initString);

  return;
}
