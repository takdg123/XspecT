// Runs the APED interpolation with all 30 elements allowed to be free parameters
//      parameters (params):
//           0:      kT temperature in keV
//           1:      kTi ion temperature in keV
//       2..31:      Metal abundances
//          32:      redshift z
//           Norm = (4 * pi * 1e14)^-1 * Int(n^2)dV / D^2 where n is
//           in cm^-3 and D is the distance in cm.

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

// from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real DEM, const Real Tb, const int ifl, 
			const bool qtherm, const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray);


void vvtapec(const RealArray& energyArray, const RealArray& params, 
	     int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	     const string& initString)
{

  RealArray abun(params.size()-3);
  IntegerVector Zarray(params.size()-3);
  for (size_t i=0; i<abun.size(); i++) {
    abun[i] = params[i+2];
    Zarray[i] = i+1;
  }

  Real T = params[0];
  Real Tb = params[1];
  Real z = params[32];
  calcMultiTempPlasma(energyArray, 6, Zarray, abun, (Real)1.0, z, T,
		      Tb, (Real)1.0, spectrumNumber, false, (Real)0.0, 
		      fluxArray, fluxErrArray);
  return;
}
