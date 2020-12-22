// XSPEC subroutine to calculate the new Mewe-Gronenschild
// plasma emission spectrum.

// parameters :
//    0.........kT (keV)
//    1.........nH (cm^-3)           fixed at 1 for most applications
//    2.........Heavy metal abundance
//    3.........Redshift


#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

void meka(const RealArray& energyArray, const RealArray& params, 
	 int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	 const string& initString)
{

  // set up parameter array for vmeka model and call it

  RealArray vparams(17);
  vparams[0] = params[0];
  vparams[1] = params[1];
  vparams[2] = 1.0;
  for (size_t i=3; i<16; i++) vparams[i] = params[2];
  vparams[16] = params[3];

  vmeka(energyArray, vparams, spectrumNumber, fluxArray, fluxErrArray, initString);

  return;
}
