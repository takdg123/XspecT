// Multiplicative super-exponential cut-off model often used for gamma-ray
// observations of pulsars
//
//   M(E) = exp(-(E/Ec)^alpha)
//
// Parameters are     cutoff       Cutoff-energy scale
//                    alpha        Exponent on exponential

#include "xsTypes.h"

extern "C" void superExpCutoff (const RealArray& energyArray, 
		       const RealArray& params, 
		       int spectrumNumber,
		       RealArray& fluxArray, 
		       RealArray& fluxErrArray,
		       const string& initString)
{

  Real Ec (params[0]);
  Real alpha (params[1]);

  Real a, b;

  fluxArray.resize(energyArray.size()-1);
  fluxErrArray.resize(0);

  a = exp(-pow(energyArray[0]/Ec,alpha));
  for (size_t i=0; i<fluxArray.size(); i++) {
    b = exp(-pow(energyArray[i+1]/Ec,alpha));
    fluxArray[i] = 0.5 * (a + b);
    a = b;
  }

  return;
}
