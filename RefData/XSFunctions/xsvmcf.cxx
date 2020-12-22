#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <iostream>
#include <sstream>


// function definition from calcCoolingFlow.cxx
void calcCoolingFlow(const RealArray& energyArray, const Real tlow, 
		     const Real thigh, const Real slope,
		     const IntegerVector& Zarray, const RealArray& abun,
		     const Real z, const int plasmaType, const int ifl,
		     RealArray& flux, RealArray& fluxErr, const Real tpeak = -1.0);



void xsvmcf(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{

  //  XSPEC model subroutine to calculate cooling flow spectrum
  //  from sum of MEKAL spectra. Variable abundances.
  //  Parameters :
  //        1..................low temperature
  //        2..................high temperature
  //        3..................He abundance
  //        4..................C     "
  //        5..................N     "
  //        6..................O     "
  //        7..................Ne    "
  //        8..................Na    "
  //        9..................Mg    "
  //       10..................Al    "
  //       11..................Si    "
  //       12..................S     "
  //       13..................Ar    "
  //       14..................Ca    "
  //       15..................Fe    "
  //       16..................Ni    "
  //       17..................redshift
  //       18..................switch(0=calculate MEKAL model, 
  //                                  1=interpolate MEKAL model
  //                                  2=APEC model)

  //  Norm is mass accretion rate in units of Msun/yr
  
  //  kaa 8/3/93      based on XSCFLW.
  //      10/4/96     added in slow and fast options

  //  cg    2/6/09    Front end and dynamic memory allocation translated
  //                  to C++.
  //  kaa 12/19/17    Modified to use calcMultiTempPlasma and associated functions


  const int elements[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};
  IntegerVector Zarray(14);
  for (size_t i=0; i<14; i++) Zarray[i] = elements[i];

  // Setting variables with parameter values

  Real tlow = params[0];
  Real thigh = params[1];
  Real slope = 0.0;
  RealArray abun(14);
  for (size_t i=0; i<14; i++) abun[i] = params[i+2];
  Real z = params[16];
  int switchPar = static_cast<int>(params[17]);

  if ( z <= 0.0 ) {
    FunctionUtility::xsWrite("\n XSVMCF: Require z > 0 for cooling flow models",10);
    return;
  }

  int plasmaType;
  if ( switchPar == 0 || switchPar == 1) {
    plasmaType = switchPar + 3;
  } else if ( switchPar == 2 ) {
    plasmaType = switchPar + 4;
  } else {
    FunctionUtility::xsWrite("\n XSVMCF: Invalid switch parameter value",2);
    FunctionUtility::xsWrite("         Must be 0, 1, or 2",2);
    return;
  }

  calcCoolingFlow(energyArray, tlow, thigh, slope, Zarray, abun, z, plasmaType,
		  spectrumNumber, flux, fluxErr);

}
