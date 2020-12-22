// Convolution model to velocity shift an additive model. The parameter is
// the velocity, v, so the shift is dE = -Ev/c.
// kaa  2/16/16

// Parameter is     v     Velocity in km/s

#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

void vashift (const RealArray& energyArray, 
	      const RealArray& params, 
	      int spectrumNumber,
	      RealArray& fluxArray, 
	      RealArray& fluxErrArray,
	      const string& initString)
{

  using namespace std;
  using namespace Numerics;
  using namespace Rebin;

  if ( params[0] == 0.0 ) return;

  Real shiftFactor (1.0 - params[0]/LIGHTSPEED);

  // set up array of shifted energies

  RealArray shiftedEnergy(energyArray*shiftFactor);

  // set up bin matching to shift the flux

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

  RealArray TempFlux (fluxArray);

  rebin(fluxArray, startBin, endBin, startWeight, endWeight, TempFlux);

  fluxArray = TempFlux;

  if ( fluxErrArray.size() > 0 ) {
    RealArray TempFluxErr (fluxErrArray);
    rebin(fluxErrArray, startBin, endBin, startWeight, endWeight, TempFluxErr);
    fluxErrArray = TempFluxErr;
  }

  return;

}
