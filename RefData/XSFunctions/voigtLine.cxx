// Voigt profile emission line
//  Parameters are :
//    0:   line center (keV)
//    1:   sigma from Gaussian  (keV)
//    2:   gamma from Lorentzian  (keV)

#include <XSFunctions/functionMap.h>
#include <xsTypes.h>
#include <XSUtil/Numerics/BinarySearch.h>
#include <cmath>

void calcLine(const RealArray& energyArray, const Real ecenter, 
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel, 
	      const int lineShape, const bool qspeedy, RealArray& fluxArray);

void voigtLine(const RealArray& energyArray, const RealArray& params, 
	       int spectrumNumber, RealArray& fluxArray, 
	       RealArray& fluxErrArray, const string& initString)
{
  const Real crtLevel = 1.0e-6;
  RealArray lineParams(1);
  lineParams[0] = params[2];

  fluxArray.resize(energyArray.size()-1);
  fluxArray = 0.0;

  // first calculate the lorentzian line
  calcLine(energyArray, params[0], lineParams, (Real)1.0, crtLevel,
	   1, false, fluxArray);
  fluxErrArray.resize(0);

  // it is also useful to calculate a gaussian centered on the same energy
  lineParams[0] = params[1];
  RealArray gaussFluxArray(fluxArray.size());
  calcLine(energyArray, params[0], lineParams, (Real)1.0, crtLevel,
	   0, false, gaussFluxArray);

  // now find the flux bins for which either the gaussian or the lorentzian is
  // above zero. The calcLine routine will have zeroed the flux when the integrated
  // flux under the line is within crtLevel of the total integrated flux.
  int ie0 = Numerics::BinarySearch(energyArray, params[0]);
  int iLow = ie0;
  while ( (iLow > 0) && (fluxArray[iLow] > 0.0) && 
	  (gaussFluxArray[iLow] > 0.0) ) iLow--;
  int iHigh = ie0;
  while ( (iHigh < (int)fluxArray.size()-1) && (fluxArray[iHigh] > 0.0) && 
	  (gaussFluxArray[iHigh] > 0.0) ) iHigh++;

  // put the lorentzian line in a temporary array
  RealArray energies(iHigh-iLow+2);
  RealArray fluxes(iHigh-iLow+1);
  for (size_t i=0; i<energies.size(); i++) energies[i] = energyArray[i+iLow];
  for (size_t i=0; i<fluxes.size(); i++) fluxes[i] = fluxArray[i+iLow];

  // now smooth this temporary array with a gaussian
  RealArray gsparams(2);
  gsparams[0] = params[1];
  gsparams[1] = 0.0;
  gsmooth(energies, gsparams, spectrumNumber, fluxes, fluxErrArray, " ");

  // insert the smoothed array back into the output fluxArray
  for (size_t i=0; i<fluxes.size(); i++) fluxArray[i+iLow] = fluxes[i];

  return;
}


