// This file is here for the benefit of any other routines still calling sumdem.
// sumdem is now just a wrapper for calcMultiTempPlasma.

#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <sstream>
#include <cfortran.h>

int calcMultiTempPlasma(const RealArray& energyArray, const int PlasmaType, 
			const IntegerVector& Zarray, const RealArray& abun, 
			const Real dens, const Real z, const RealArray& Tarr, 
			const RealArray& DEMarr, const int ifl, 
			const bool qtherm, const Real velocity, 
			RealArray& fluxArray, RealArray& fluxErrArray);


void sumdem(int itype, int swtch, float* ear, int ne, float* abun,
            float dens, float z, int ninputt, float* inputt, float* dem,
            int ifl, bool qtherm, float velocity, float* photar, float* photer,
	    int* status);
FCALLSCSUB16(sumdem,SUMDEM,sumdem,INT,INT,FLOATV,INT,FLOATV,FLOAT,FLOAT,INT,FLOATV,FLOATV,INT,LOGICAL,FLOAT,FLOATV,FLOATV,PINT)




void sumdem(int itype, int swtch, float* ear, int ne, float* abun,
            float dens, float z, int ninputt, float* inputt, float* dem,
            int ifl, bool qtherm, float velocity, float* photar, float* photer,
	    int* status)
{

// Wrapper for calcMultiTempPlasma allowing it to be called from Fortran
// Arguments :
//      itype   I        i: type of plasma emission file
//                           1 = R-S
//                           2 = Mekal
//                           3 = Meka
//                           4 = APEC
//                           5 = APEC with all elements given abundances
//      swtch  I         i: 0==calculate, 1==interpolate, 2=apec interpolate
//      ear     R        i: model energy ranges
//      ne      I        i: number of model energies
//      abun    R        i: abundances
//      dens    R        i: density (cm^-3)
//      z       R        i: redshift
//      ninputt I        i: number of temperatures
//      inputt  R        i: temperatures
//      dem     R        i: emission measures for input temperatures
//      ifl     I        i: dataset number (unused at present)
//      qtherm  R        i: apply thermal broadening (APEC models only)
//      velocityR        i: gaussian velocity broadening (APEC models only)
//      photar  R        r: spectrum
//      status  I        r: 0==OK
// 

  const int TOTEL=30, MKEL=14, RSEL=12, APEL=13;
  const int rselt[] = {2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28};
  const int apelt[] = {2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28};
  const int mkelt[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};

  // set up input energy array for calcMultiTempPlasma

  RealArray energyArray(ne+1);
  for (size_t ie=0; ie<energyArray.size(); ie++) energyArray[ie] = ear[ie];

  // set up abundance and atomic number arrays

  size_t nAbun;
  if ( itype == 1 ) nAbun = RSEL;
  if ( itype == 2 || itype == 3 ) nAbun = MKEL;
  if ( itype == 4 ) nAbun = APEL;
  if ( itype == 5 ) nAbun = TOTEL;
  RealArray abundance(nAbun);
  for (size_t i=0; i<nAbun; i++) abundance[i] = abun[i];
  IntegerVector Zarray(nAbun);
  if ( itype == 1 ) {
    for (size_t i=0; i<nAbun; i++) Zarray[i] = rselt[i];
  } else if ( itype == 2 || itype == 3 ) {
    for (size_t i=0; i<nAbun; i++) Zarray[i] = mkelt[i];
  } else if ( itype == 4 ) {
    for (size_t i=0; i<nAbun; i++) Zarray[i] = apelt[i];
  } else if ( itype == 5 ) {
    for (size_t i=0; i<nAbun; i++) Zarray[i] = i+1;
  }

  // define the plasmaType variable

  int plasmaType;
  if ( itype == 1 ) plasmaType = 1;
  if ( itype == 2 ) {
    if ( swtch == 0 ) plasmaType = 3;
    if ( swtch == 1 ) plasmaType = 4;
    if ( swtch == 2 ) plasmaType = 6;
  }
  if ( itype == 3 ) plasmaType = 5;
  if ( itype == 4 || itype == 5 ) plasmaType = 6;

  // set up the temperature and DEM arrays

  RealArray Tarr(ninputt);
  RealArray DEMarr(ninputt);
  for (size_t i=0; i<(size_t)ninputt; i++) {
    Tarr[i] = inputt[i];
    DEMarr[i] = dem[i];
  }

  // calculate the model spectrum

  RealArray fluxArray, fluxErrArray;

  *status = calcMultiTempPlasma(energyArray, plasmaType, Zarray, abundance, (Real)dens, 
				(Real)z, Tarr, DEMarr, ifl, qtherm, (Real)velocity,
				fluxArray, fluxErrArray);

  // transfer the model spectrum to photar
  for (size_t ie=0; ie<fluxArray.size(); ie++) photar[ie] = fluxArray[ie];

  return;

}
