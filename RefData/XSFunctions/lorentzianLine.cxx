// Based on gaussianLine.cxx but for a Lorentzian

#include <XSFunctions/functionMap.h>
#include <xsTypes.h>
#include <XSstreams.h>
#include <cmath>

void calcLine(const RealArray& energyArray, const Real ecenter, 
	      const RealArray& lineParams, const Real lineflux, const Real crtLevel, 
	      const int lineShape, const bool qspeedy, RealArray& fluxArray);

void lorentzianLine(const RealArray& energyArray, const RealArray& params, 
		    int spectrumNumber, RealArray& fluxArray, 
		    RealArray& fluxErrArray, const string& initString)
{
  const Real crtLevel = 1.0e-6;
  RealArray lineParams(params[1],1);

  fluxArray.resize(energyArray.size()-1);
  fluxArray = 0.0;

  calcLine(energyArray, params[0], lineParams, (Real)1.0, crtLevel,
	   1, false, fluxArray);
  fluxErrArray.resize(0);

  return;
}


