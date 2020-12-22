// Convolution model to redshift an additive model. Includes the (1+z)
// energy correction and (1+z) time correction.
// kaa  3/7/11

// Parameter is     z     Redshift

#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

void zashift (const RealArray& energyArray, 
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

  // rebin into output temporary arrays then copy to output arrays
  // including the time redshift factor

  RealArray TempFlux (fluxArray);

  rebin(fluxArray, startBin, endBin, startWeight, endWeight, TempFlux);

  fluxArray = TempFlux*zfactor;

  if ( fluxErrArray.size() > 0 ) {
    RealArray TempFluxErr (fluxErrArray);
    rebin(fluxErrArray, startBin, endBin, startWeight, endWeight, TempFluxErr);
    fluxErrArray = TempFluxErr*zfactor;
  }

  return;

}
