// XSPEC subroutine to calculate the new Mewe-Gronenschild
// plasma emission spectrum.

// parameters :
//    0.........kT (keV)
//    1.........nH (cm^-3)           fixed at 1 for most applications
//    2.........He abundance
//    3.........C     "
//    4.........N     "
//    5.........O     "
//    6.........Ne    "
//    7.........Na    "
//    8.........Mg    "
//    9.........Al    "
//   10.........Si    "
//   11.........S     "
//   12.........Ar    "
//   13.........Ca    "
//   14.........Fe    "
//   15.........Ni    "
//   16.........Redshift


#include "xsTypes.h"
#include <XSFunctions/functionMap.h>

// from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const Real T, 
			const Real DEM, const int ifl, const bool qtherm, 
			const Real velocity, RealArray& fluxArray, 
			RealArray& fluxErrArray);


void vmeka(const RealArray& energyArray, const RealArray& params, 
 	   int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	   const string& initString)
{
  const size_t nMK(15);
  const int mkelt[] = {1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};

  RealArray abun(nMK);
  IntegerVector Zarray(nMK);
  for(size_t i=0; i<nMK; i++) Zarray[i] = mkelt[i];
  abun[0] = 1.0;
  for(size_t i=1; i<nMK; i++) abun[i] = params[i+1];

  Real z = params[16];
  Real T = params[0];
  calcMultiTempPlasma(energyArray, 5, Zarray, abun, params[1], z, T, 
		      (Real)1.0, spectrumNumber, false, (Real)0.0,
		      fluxArray, fluxErrArray);
  return;
}
