// Convolution model with a lorentzian kernel whose width can vary with energy
// Parameters:
//     0    Lorentzian sigma at 6 keV
//     1    Index for energy dependence of sigma

#include <xsTypes.h>
#include <functionMap.h>


// Uses the calcManyLines function from calcLines.cxx

void calcManyLines(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevel,
		   const int lineShape, const bool qspeedy, RealArray& fluxArray);

void lsmooth(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
  size_t nE = energyArray.size();

  // construct the array of lorentzian sigmas
  std::vector<RealArray> width(nE-1);
  RealArray centerEnergy(nE-1);
  for (size_t i=0; i<nE-1; i++) {
    width[i].resize(1);
    centerEnergy[i] = 0.5*(energyArray[i]+energyArray[i+1]);
    width[i][0] = params[0]*pow((centerEnergy[i]/6.0),params[1]);
  }

  // do the smoothing by just treating the input flux as the sum of a lot of
  // lorentzians with the center in each energy bin and the normalization based
  // on the flux in that bin

  RealArray inputFlux(flux);
  flux = 0.0;
  calcManyLines(energyArray, centerEnergy, width, inputFlux, (Real)1.0e-6, 1,
		true, flux);
}
