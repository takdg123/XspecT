// raymond smith model with all abundances free using data from model file
// for this model. The file is in the general format for table model files.
//   parameters (param):
//     0:      kT temperature in keV
// 1..12:      Abundances for He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni 
//             wrt Solar (defined by the abund command)
//    13:      redshift z
//     Norm = (4 * pi * 1e14)^-1 * Int(n^2)dV / D^2 where n is
//             in cm^-3 and D is the distance in cm.


#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

// from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real DEM, const int ifl, const bool qtherm, 
			const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray);


void vraysmith(const RealArray& energyArray, const RealArray& params, 
	       int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	       const string& initString)
{
  size_t nRS(13);
  const int rselt[] = {1, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28};
  bool newVersion(false);

  if ( !newVersion ) nRS--;

  RealArray abun(nRS);
  IntegerVector Zarray(nRS);

  if ( newVersion ) {
    for(size_t i=0; i<nRS; i++) Zarray[i] = rselt[i];
    // Add H
    abun[0] = 1.0;
    // Set the other abundances
    for(size_t i=1; i<nRS; i++) abun[i] = params[i];
  } else {
    for(size_t i=0; i<nRS; i++) Zarray[i] = rselt[i+1];
    for(size_t i=0; i<nRS; i++) abun[i] = params[i+1];
  }


  Real z = params[13];
  Real T = params[0];
  // uses plasmaType=2 which is new the APED-style RS files. May want to add an option
  // here to get back the old plasmaType=1 using raysmith.mod. Note that if we do use
  // the old version then we cannot have H in the abun and Zarray arrays.

  if ( newVersion ) {
    calcMultiTempPlasma(energyArray, 2, Zarray, abun, (Real)1.0, z, T, 
			(Real)1.0, spectrumNumber, false, (Real)0.0, 
			fluxArray, fluxErrArray);
  } else {
    calcMultiTempPlasma(energyArray, 1, Zarray, abun, (Real)1.0, z, T, 
			(Real)1.0, spectrumNumber, false, (Real)0.0, 
			fluxArray, fluxErrArray);
  }

  return;
}
