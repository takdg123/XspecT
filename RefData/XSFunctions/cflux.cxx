// Convolution model that can be used to create a parameter which is the flux
// in a particular energy range
// Parameters are     energ_lo     Low energy over which to calculate flux
//                    energ_hi     High energy over which to calculate flux
//                    logflux      Log of flux (cgs)

#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

void cflux (const RealArray& energyArray, 
	    const RealArray& params, 
	    int spectrumNumber,
	    RealArray& fluxArray, 
	    RealArray& fluxErrArray,
	    const string& initString)
{

  using namespace std;
  using namespace Numerics;

  Real eMin (params[0]);
  Real eMax (params[1]);
  Real flux (pow(10.0, params[2]));

  // Integrate the input array between elow and ehi into fluxsum.

  pair<Real,Real> fluxsum (integrationKernel(energyArray,fluxArray,eMin,eMax));

  // Scale by flux/fluxsum

  fluxArray *= flux/fluxsum.second;

  return;

}
