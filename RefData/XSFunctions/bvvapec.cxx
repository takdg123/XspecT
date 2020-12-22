// Runs the APED interpolation with velocity and thermal broadening and all abundances
// free to be fit
//      parameters (params):
//           0:      kT temperature in keV
//       1..30:      Metal abundances
//          31:      redshift z
//          32:      gaussian velocity width
//           Norm = (4 * pi * 1e14)^-1 * Int(n^2)dV / D^2 where n is
//           in cm^-3 and D is the distance in cm.

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

// from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real DEM, const int ifl, const bool qtherm, 
			const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray);


void bvvapec(const RealArray& energyArray, const RealArray& params, 
	    int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	    const string& initString)
{

  RealArray abun(params.size()-3);
  IntegerVector Zarray(params.size()-3);
  for (size_t i=0; i<abun.size(); i++) {
    abun[i] = params[i+1];
    Zarray[i] = i+1;
  }

  Real T = params[0];
  Real z = params[31];
  Real v = params[32];
  calcMultiTempPlasma(energyArray, 6, Zarray, abun, (Real)1.0, z, T, 
		      (Real)1.0, spectrumNumber, true, v, fluxArray, 
		      fluxErrArray);
  return;
}
