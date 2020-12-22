// Convolution model to redshift a multiplicative model. Just does the (1+z)
// energy correction.
// kaa  3/7/11

// Parameter is     z     Redshift

#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

void zmshift (const RealArray& energyArray, 
	      const RealArray& params, 
	      int spectrumNumber,
	      RealArray& fluxArray, 
	      RealArray& fluxErrArray,
	      const string& initString)
{

  using namespace std;
  using namespace Numerics;
  using namespace Rebin;

  Real zfactor (1.0/(1.0+params[0]));

  // set up array of shifted energies for the redshift

  RealArray shiftedEnergy(energyArray*zfactor);

  // set up bin matching to go from rest to observed frame

  size_t inputBin;
  size_t outputBin;

  IntegerVector startBin(fluxArray.size());
  IntegerVector endBin(fluxArray.size());
  RealArray startWeight(fluxArray.size());
  RealArray endWeight(fluxArray.size());

  const Real FUZZY = 1.0e-6;

  findFirstBins(shiftedEnergy, energyArray, FUZZY, inputBin, outputBin);
  initializeBins(shiftedEnergy, energyArray, FUZZY, inputBin, outputBin, 
		 startBin, endBin, startWeight, endWeight);

  // interpolate onto output temporary arrays

  RealArray TempFlux (fluxArray);

  interpolate(fluxArray, startBin, endBin, startWeight, endWeight, 
	      TempFlux, false);

  fluxArray = TempFlux;

  if ( fluxErrArray.size() > 0 ) {
    RealArray TempFluxErr (fluxErrArray);
    interpolate(fluxErrArray, startBin, endBin, startWeight, endWeight, 
		TempFluxErr, false);
    fluxErrArray = TempFluxErr;
  }

  return;

}
