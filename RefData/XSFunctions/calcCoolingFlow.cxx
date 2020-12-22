// Common code to calculate all the cooling flow model options

#include <functionMap.h>
#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <iostream>
#include <sstream>

#define keVtoK 1.1604505e7

// function definitions calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
                        const IntegerVector& Zarray, const RealArray& abun, 
                        const Real dens, const Real z, const RealArray& Tarr, 
                        const RealArray& DEMarr, const int ifl, const bool qtherm, 
                        const Real velocity, RealArray& fluxArray, 
                        RealArray& fluxErrArray);
int plasmaFileInfo(const int plasmaType, RealArray& energyArray, RealArray& Tarray);
int plasmaBolo(const int plasmaType, const IntegerVector& Zarray, 
	       const RealArray& abun, const Real dens, const Real z,
	       const RealArray& Tarray, const int ifl, const bool qtherm,
	       const Real velocity, RealArray& boloArray);

// function definition for old version (CFLOW_VERSION 1) included later in the file
void calcOldCoolingFlow(const RealArray& energyArray, const Real tlow, 
			const Real thigh, const Real slope, const IntegerVector& Zarray,
			const RealArray& abun, const Real z, const int plasmaType, 
			const int ifl, RealArray& flux, RealArray& fluxErr);



void calcCoolingFlow(const RealArray& energyArray, const Real tlow, 
		     const Real thigh, const Real slope, const IntegerVector& Zarray,
		     const RealArray& abun, const Real z, const int plasmaType, 
		     const int ifl, RealArray& flux, RealArray& fluxErr,
		     const Real tpeak)
{

  //  Routine arguments
  //       energyArray      standard XSPEC energy array
  //       tlow             minimum temperature in cooling flow
  //       thigh            maximum temperature in cooling flow
  //       slope            slope of emissivity relation wrt bolometric lum.
  //       Zarray           Atomic numbers of elements
  //       abun             Abundances of elements in Zarray
  //       z                redshift (must be >0)
  //       plasmaType       1 = R-S  table model (old version)
  //                        2 = R-S  APED-style model
  //                        3 = Mekal  calculate
  //                        4 = Mekal  interpolate
  //                        5 = Meka
  //                        6 = APEC
  //       ifl              file number (probably not required)
  //       flux             returned flux
  //       fluxErr          error on flux
  //       tpeak            if >0 add heating term (for cph model)


  //  Norm is mass accretion rate in units of Msun/yr
  
  //  kaa 8/3/93      based on XSCFLW.
  //      10/4/96     added in slow and fast options

  //  cg    2/6/09    Front end and dynamic memory allocation translated
  //                  to C++.
  //  kaa 12/19/17    Modified to use calcMultiTempPlasma and associated functions
  //  kaa 12/20/17    Broken out from xsvmcf


  const size_t nE = energyArray.size()-1;
  flux.resize(nE);
  fluxErr.resize(0);

  Real dens = 1.0;

  if ( z <= 0.0 ) {
    FunctionUtility::xsWrite("\n calcCoolingFlow: Require z > 0 for cooling flow models",10);
    return;
  }

  if ( plasmaType < 0 || plasmaType > 6) {
    FunctionUtility::xsWrite("\n calcCoolingFlow: Invalid plasmaType value",2);
    FunctionUtility::xsWrite("                    Must be between 1 and 6",2);
    return;
  }

  // Check for a version number

  int cflver = 2;
  const string& verStr(FunctionUtility::getModelString("CFLOW_VERSION"));
  if (verStr.length() && verStr != FunctionUtility::NOT_A_KEY()) {
    std::istringstream iss(verStr);
    int test=0;
    if (!(iss >> test) || !iss.eof()) {
      std::ostringstream err;
      err << "Invalid CFLOW_VERSION value: " << verStr
	  <<"\n will use version = " << cflver;
      FunctionUtility::xsWrite(err.str(), 10);
    } else {
      cflver = test;
    }
  }
  if ( cflver == 1 || cflver == 2 ) {
    std::ostringstream oss;
    oss << "Cooling flow model version " << cflver;
    FunctionUtility::xsWrite(oss.str(), 25);
  } else {
    std::ostringstream err;
    err << "Invalid CFLOW version number: " << cflver;
    FunctionUtility::xsWrite(err.str(), 10);
    return;
  }

  // If we want the old version then just call the routine and return

  if ( cflver == 1 ) {
    calcOldCoolingFlow(energyArray, tlow, thigh, slope, Zarray, abun, z,
		       plasmaType, ifl, flux, fluxErr);
    return;
  }
  
  // Get the number of temperature steps. For the new version use 
  // either 10 or the value set by CFLOW_NTEMPS

  size_t nsteps = 0;
  RealArray ttab;
  nsteps = 10;
  const string& nTempsStr(FunctionUtility::getModelString("CFLOW_NTEMPS"));
  if (nTempsStr.length() && nTempsStr != FunctionUtility::NOT_A_KEY()) {
    std::istringstream iss(nTempsStr);
    int test=0;
    if (!(iss >> test) || !iss.eof()) {
      std::ostringstream err;
      err << "Invalid CFLOW_NTEMPS value: " << nTempsStr
	  <<"\n will use ntemps = " << nsteps;
      FunctionUtility::xsWrite(err.str(), 10);
    } else {
      nsteps = test;
    }
  }

  // Set up the temperature array

  Real altlow = log(tlow);
  Real althigh = log(thigh);

  RealArray tval(nsteps);

  for (size_t i=0; i<nsteps; i++) {
    tval[i] = exp(altlow + i*(althigh-altlow)/(nsteps-1));
  }

  // Get the cosmology parameters

  Real q0 = FunctionUtility::getq0();
  Real h0 = FunctionUtility::getH0();
  Real Lambda0 = FunctionUtility::getlambda0();

  // Fix up the norm. The numerical constant assumes distance linearly depends
  // on redshift with H0=50 hence cosmology factors correct this.

  Numerics::FZSQ fzsq;
  Real norm = 3.16e-15 * (h0/50.)*(h0/50.) / fzsq(z, q0, Lambda0);

  //  norm is C/4piD^2 in Mushotzky & Szymkowiak - uses :
  //	L = (5/2) (Mdot/mu/m(H)) k dT
  //  and :
  //   3.16e-15 = (5/2) (1.989e33/3.15e7) (1/0.6) (1/1.66e-24)
  //              (1/4pi) (6000 3.09e24)^-2 (1.38e-16/1.60e-9)


  // Calculate the bolometric luminosities for each temperature.

  RealArray bolo(nsteps);
  plasmaBolo(plasmaType, Zarray, abun, dens, 0.0, tval, ifl, false, 0.0, bolo);

  // In the heating case we need the bolometric luminosity at the peak temperature

  Real sqrtTpeak5K, tpeakBolo;
  if ( tpeak > 0.0 ) {
    RealArray tbolo(1);
    plasmaBolo(plasmaType, Zarray, abun, dens, 0.0, RealArray(1,tpeak), ifl,
	       false, 0.0, tbolo);
    tpeakBolo = tbolo[0];
    sqrtTpeak5K = pow(keVtoK*tpeak,2.5);
  }

  //  Now calculate the DEMs. The integral is an extended trapezium rule 
  //  over the temperatures. The integration is performed in log T space.

  RealArray dem(nsteps);
  for (size_t i=0; i<nsteps; i++) dem[i] = 0.0;

  //  Now loop over the temperatures.

  for ( size_t i=0; i<nsteps; i++) {

    Real factor = (althigh-altlow)/(nsteps-1);
    if ( i == 0 || i == nsteps-1 ) factor = factor/2.;

    //  calculate the emission-weighting for this temperature. This
    //  includes the division by the bolometric emissivity and an
    //  extra factor of T since we are integrating in log space.

    Real tstep = keVtoK*tval[i];
    dem[i] = pow((tstep/1.e7),slope) * tstep * factor * norm / bolo[i];

    // if the heating term is included for tpeak>0

    if ( tpeak > 0.0 ) {
      Real sqrtTstep5 = pow(tstep,2.5);
      dem[i] *= 1.0/fabs(1-sqrtTstep5/bolo[i]/(sqrtTpeak5K/tpeakBolo));
    }

  }

  // Now call the routine that does the actual calculation of the
  // spectrum

  calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, dens, z,
		      tval, dem, ifl, false, 0.0, flux, fluxErr);

}


// Old version (CFLOW_VERSION 1)

void calcOldCoolingFlow(const RealArray& energyArray, const Real tlow, 
			const Real thigh, const Real slope,
			const IntegerVector& Zarray, const RealArray& abun,
			const Real z, const int plasmaType, 
			const int ifl, RealArray& flux, RealArray& fluxErr)
{

  //  Routine arguments
  //       energyArray      standard XSPEC energy array
  //       tlow             minimum temperature in cooling flow
  //       thigh            maximum temperature in cooling flow
  //       slope            slope of emissivity relation wrt bolometric lum.
  //       Zarray           Atomic numbers of elements
  //       abun             Abundances of elements in Zarray
  //       z                redshift (must be >0)
  //       plasmaType       1 = R-S  table model (old version)
  //                        2 = R-S  APED-style model
  //                        3 = Mekal  calculate
  //                        4 = Mekal  interpolate
  //                        5 = Meka
  //                        6 = APEC
  //       ifl              file number (probably not required)
  //       flux             returned flux
  //       fluxErr          error on flux


  //  Norm is mass accretion rate in units of Msun/yr
  
  //  kaa 8/3/93      based on XSCFLW.
  //      10/4/96     added in slow and fast options

  //  cg    2/6/09    Front end and dynamic memory allocation translated
  //                  to C++.
  //  kaa 12/19/17    Modified to use calcMultiTempPlasma and associated functions
  //  kaa 12/20/17    Broken out from xsvmcf

  const size_t nE = energyArray.size()-1;
  flux.resize(nE);
  fluxErr.resize(0);

  Real dens = 1.0;

  if ( z <= 0.0 ) {
    FunctionUtility::xsWrite("\n calcCoolingFlow: Require z > 0 for cooling flow models",10);
    return;
  }

  if ( plasmaType < 0 || plasmaType > 6) {
    FunctionUtility::xsWrite("\n calcCoolingFlow: Invalid plasmaType value",2);
    FunctionUtility::xsWrite("                    Must be between 1 and 6",2);
    return;
  }

  // Get the number of temperature steps.

  size_t nsteps = 0;
  RealArray ttab;
  RealArray etab;
  plasmaFileInfo(plasmaType, etab, ttab);
  nsteps = ttab.size();

  // Set up the temperature array

  Real altlow = log(tlow);
  Real althigh = log(thigh);

  RealArray tval(nsteps);

  for (size_t i=0; i<nsteps; i++) tval[i] = ttab[i];

  // Get the cosmology parameters

  Real q0 = FunctionUtility::getq0();
  Real h0 = FunctionUtility::getH0();
  Real Lambda0 = FunctionUtility::getlambda0();

  // Fix up the norm. The numerical constant assumes distance linearly depends
  // on redshift with H0=50 hence cosmology factors correct this.

  Numerics::FZSQ fzsq;
  Real norm = 3.16e-15 * (h0/50.)*(h0/50.) / fzsq(z, q0, Lambda0);

  //  norm is C/4piD^2 in Mushotzky & Szymkowiak - uses :
  //	L = (5/2) (Mdot/mu/m(H)) k dT
  //  and :
  //   3.16e-15 = (5/2) (1.989e33/3.15e7) (1/0.6) (1/1.66e-24)
  //              (1/4pi) (6000 3.09e24)^-2 (1.38e-16/1.60e-9)


  // Calculate the bolometric luminosities for each temperature.

  RealArray bolo(nsteps);
  plasmaBolo(plasmaType, Zarray, abun, dens, 0.0, tval, ifl, false, 0.0, bolo);

  //  Now calculate the DEMs. The integral is an extended trapezium rule 
  //  over the temperatures. The integration is performed in log T space.

  RealArray dem(nsteps);
  for (size_t i=0; i<nsteps; i++) dem[i] = 0.0;

  //  Now loop over the temperatures. In the old version we perform the
  //  function evaluations at tabulated temperatures. The integral
  //  is an extended trapezium rule over the tabulated temperatures
  //  within the range then simple trapezium rules to handle the bits
  //  left over at the end. The integration is performed in log T space.

  //  Find the index for the tabulated temperatures immediately above the 
  //  minimum and below the maximum. tmin and tmax are the log T corresponding
  //  to these. tbot is the log of the lowest tabulated T.

  //  NB this assumes that the temperatures are logarithmically evenly 
  //  distributed

  Real deltat = (log(tval[nsteps-1])-log(tval[0]))/(nsteps-1);
  Real tbot = log(tval[0]);
  size_t itmin = (size_t)((altlow-tbot)/deltat + 1);
  if ( itmin > 1 ) itmin = 1;
  size_t itmax = (size_t)((althigh-tbot)/deltat);
  if ( itmax > nsteps - 2 ) itmax = nsteps - 2;
  Real tmin = tbot + deltat*(itmin);
  Real tmax = tbot + deltat*(itmax);

  for (size_t i=(itmin-1); i<=(itmax+1); i++) {

    //  Calculate the weighting factor in the integral (this is
    //  a bit messy !).

    Real factor = deltat;
    if (i == itmin-1) {
      factor *= (1.+(tmin-altlow)/deltat)/2;
    } else if (i == itmin) {
      factor *= (1.+(altlow-tmin+deltat)/deltat)/2;
    } else if (i == itmax) {
      factor *= (1.+(tmax+deltat-althigh)/deltat)/2;
    } else if (i == itmax+1) {
      factor *= (1.+(althigh-tmax)/deltat)/2;
    }
    
    //  calculate the emission-weighting for this temperature. This
    //  includes the division by the bolometric emissivity and an
    //  extra factor of T since we are integrating in log space.

    Real tstep = 1.16e7*tval[i];
    dem[i] = pow((tstep/1.e7),slope) * tstep * factor * norm / bolo[i];

  }

  // Now call the routine that does the actual calculation of the
  // spectrum

  calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, dens, z,
		      tval, dem, ifl, false, 0.0, flux, fluxErr);

}
