// gaussian model in wavelength space using Angstrom units
// just calls calcGaussianLine

#include <xsTypes.h>
#include <functionMap.h>
#include "Numerics.h"

void calcLine(const RealArray& energyArray, const Real ecenter, 
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel, 
	      const int lineShape, const bool qspeedy, RealArray& fluxArray);

void agauss(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& fluxArray, 
	    RealArray& fluxErrArray, const string& initString)
{
  // convert the energy from keV to Angstroms. Need to make sure that we are
  // in increasing order of wavelength to be on the safe side.
  size_t nE = energyArray.size();
  RealArray angstromArray(nE);
  for (size_t i=0; i<nE; i++) angstromArray[i] = Numerics::KEVTOA/energyArray[nE-1-i];

  const Real crtLevel = 1.0e-6;
  RealArray lineParams(params[1],1);

  RealArray angstromFluxArray(nE-1);
  calcLine(angstromArray, params[0], lineParams, (Real)1.0, crtLevel, 0,
	   false, angstromFluxArray);

  fluxArray.resize(energyArray.size()-1);
  for (size_t i=0; i<nE-1; i++) fluxArray[i] = angstromFluxArray[nE-2-i];
  fluxErrArray.resize(0);
}
