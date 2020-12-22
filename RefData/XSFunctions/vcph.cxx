#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <iostream>
#include <sstream>

// function definition in calcCoolingFlow.cxx
void calcCoolingFlow(const RealArray& energyArray, const Real tlow, 
		     const Real thigh, const Real slope,
		     const IntegerVector& Zarray, const RealArray& abun,
		     const Real z, const int plasmaType, const int ifl,
		     RealArray& flux, RealArray& fluxErr, const Real tpeak);


void vcph(const RealArray& energyArray, const RealArray& params,
	  int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	  const string& initString)
{

  //  XSPEC model subroutine to calculate modified cooling flow spectrum
  //  from sum of spectra. Variable abundances.
  //  params array :
  //        0..................peak temperature
  //        1..................He abundance
  //        2..................C     "
  //        3..................N     "
  //        4..................O     "
  //        5..................Ne    "
  //        6..................Na    "
  //        7..................Mg    "
  //        8..................Al    "
  //        9..................Si    "
  //       10..................S     "
  //       11..................Ar    "
  //       12..................Ca    "
  //       13..................Fe    "
  //       14..................Ni    "
  //       15..................redshift
  //       16..................switch(0=calculate MEKAL model, 
  //                                  1=interpolate MEKAL model
  //                                  2=APEC model)

  //  Norm is mass accretion rate in units of Msun/yr

  const int elements[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};
  IntegerArray Zarray(14);
  for (size_t i=0; i<14; i++) Zarray[i] = elements[i];

  // integrate over the whole temperature range
  
  Real tlow = 0.01;
  Real thigh = 50.0;
  
  Real slope = 0.0;
  RealArray abun(14);
  for (size_t i=0; i<14; i++) abun[i] = params[i+1];
  Real z = params[15];
  Real tpeak = params[0];

  int switchPar = static_cast<int>(params[16]);
  int plasmaType;
  if ( switchPar == 0 || switchPar == 1) {
    plasmaType = switchPar + 3;
  } else if ( switchPar == 2 ) {
    plasmaType = switchPar + 4;
  } else {
    FunctionUtility::xsWrite("\n vcph: Invalid switch parameter value",2);
    FunctionUtility::xsWrite("         Must be 0, 1, or 2",2);
    return;
  }

  calcCoolingFlow(energyArray, tlow, thigh, slope, Zarray, abun, z, plasmaType,
		  spectrumNumber, flux, fluxErr, tpeak);

}
