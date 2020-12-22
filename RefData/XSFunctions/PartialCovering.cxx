// Convolution model for partial covering. This should be used to modify a multiplicative
// model using a syntax (PM)A where P is this model, M is a multiplicative model and A
// an additive component. Note that this is not the same as P(MA).
//
// Parameters are     fract        Covering fraction

#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

void PartialCovering (const RealArray& energyArray, 
		      const RealArray& params, 
		      int spectrumNumber,
		      RealArray& fluxArray, 
		      RealArray& fluxErrArray,
		      const string& initString)
{

  using namespace std;

  Real fract (params[0]);
  Real fractc (1.0-fract);

  size_t Nbins (fluxArray.size());

  for (size_t i=0; i<Nbins; i++) {
    fluxArray[i] = fractc + fract * fluxArray[i];
  }

  if ( fluxErrArray.size() > 0 ) {
    for (size_t i=0; i<Nbins; i++) {
      fluxErrArray[i] = fractc + fract * fluxErrArray[i];
    }
  }

  return;

}
