// logpar model
//   Parameters are :
//      1    alpha
//      2    beta
//      3    pivot energy (E0: emitted frame: keV)
//
//  N(E) = (E/E0)**(-alpha-beta*log(E/E0))
//
// redshifted version
//   Parameters are :
//      1    alpha
//      2    beta
//      3    pivot energy (E0: emitted frame: keV)
//      4    redshift
//
//  N(E) = ([E(1+z)]/E0)**(-alpha-beta*log([E(1+z)]/E0))

#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

// prototype in this file

void calcLogpar(const RealArray& energyArray, const RealArray& params,
		RealArray& fluxArray);
// ***************************************************************

void logpar (const RealArray& energyArray, 
	      const RealArray& params, 
	      int spectrumNumber,
	      RealArray& fluxArray, 
	      RealArray& fluxErrArray,
	      const string& initString)
{
  fluxErrArray.resize(0);

  calcLogpar(energyArray, params, fluxArray);

  return;
}

// ***************************************************************

void zLogpar (const RealArray& energyArray, 
	      const RealArray& params, 
	      int spectrumNumber,
	      RealArray& fluxArray, 
	      RealArray& fluxErrArray,
	      const string& initString)
{
  fluxErrArray.resize(0);

  Real zfactor = (1.0+params[3]);
  RealArray bparams(3);
  for (size_t i=0; i<3; i++) bparams[i] = params[i];

  const RealArray energy(energyArray*zfactor);
  
  calcLogpar(energy, bparams, fluxArray);
  fluxArray /= zfactor;

  return;
}

// ***************************************************************

void calcLogpar(const RealArray& energyArray, const RealArray& params,
		RealArray& fluxArray)
{
  fluxArray.resize(energyArray.size()-1);
  
  const Real h(0.5/3.0);      // constant for Simpson's rule

  const Real alpha (params[0]);
  const Real beta (params[1]);
  const Real e0 (params[2]);

  Real e = energyArray[0]/e0;
  Real a = pow(e, -alpha-beta*log10(e));
  for (size_t i=0; i<fluxArray.size(); i++) {
    Real ea = 0.5*(energyArray[i]+energyArray[i+1])/e0;
    Real b = pow(ea, -alpha-beta*log10(ea));
    e = energyArray[i+1]/e0;
    Real c = pow(e, -alpha-beta*log10(e));
    // Simpson's rule integration
    fluxArray[i] = h*(energyArray[i+1]-energyArray[i])*(a+c+4*b);
    a = c;
  }

  return;
}
