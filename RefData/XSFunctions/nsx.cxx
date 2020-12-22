//-----------------------------------------------------------------------------
//    Model neutron star spectrum for unmagnetized, partially ionized atmosphere
//    (see Ho, WCG, Heinke, CO 2009, Nature, 462, 71)
//     Converted to C++ and to use a single input FITS file by kaa on 5/17/17
//
//    parameter[0] = log(unredshifted effective temperature, in K)
//    parameter[1] = neutron star mass, in solar masses
//    parameter[2] = neutron star radius, in km
//    parameter[3] = distance to neutron star, in kpc
//    parameter[4] = model to use, see MODELINFO extension of nsx.fits
//-----------------------------------------------------------------------------

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>
#include "XSFunctions/Utilities/FunctionUtility.h"
#include "BinarySearch.h"
#include <CCfits/CCfits>

using namespace std;

// routine to do the interpolation over temperature and surface gravity
// returns tableEnergy and tableFlux. code is in nsmaxg.cxx.
void calcHoNSatmos(Real& temperature, Real& gg, const string filename,
		   const int imodel, RealArray& tableEnergy,
		   RealArray& tableFlux);

// spline routines (should these be in XSUtil?)
void DoSpline(const RealArray& x, const RealArray& y, Real& firstDeriv, Real& lastDeriv,
	      RealArray& y2);
void DoSplint(const RealArray& x, const RealArray& y, const RealArray& y2, 
	      const Real xtarg, Real& yout);


void nsx(const RealArray& energy, const RealArray& parameter, int spectrum, 
	 RealArray& flux, RealArray& fluxError, const string& init)
{

  Real temperature(parameter[0]);
  const Real redshift(1.0/sqrt(1.0-2.95316*(parameter[1]/parameter[2])));
  Real gg(log10(1.3271e16*((parameter[1]/(parameter[2]*parameter[2]))*redshift)));
  const int imodel(parameter[4]);

  // check for NSX_DIR and set the filename
  string directory = FunctionUtility::getModelString("NSX_DIR");
  if (directory.compare("$$NOT$$") == 0) {
    directory = FunctionUtility::modelDataPath();
  }
  const string filename(directory+"nsx.fits");

  RealArray tableEnergy, tableFlux;
  calcHoNSatmos(temperature, gg, filename, imodel, tableEnergy, tableFlux);

  // interpolate onto the output flux array using spline and 5-pt Simpson formula.

  RealArray splineFlux;
  Real deriv1 = 1.0e32;
  Real derivN = 1.0e32;
  DoSpline(tableEnergy, tableFlux, deriv1, derivN, splineFlux);

  flux.resize(energy.size()-1);
  int Nenergies(tableEnergy.size());
  for (size_t ie=0; ie<energy.size()-1; ie++) {
    RealArray ce(5), zE(5);
    for (size_t i=0; i<5; i++) {
      zE[i] = ((1.0-i/4.0)*energy[ie] + (i/4.0)*energy[ie+1])*redshift;
      if ( zE[i] < tableEnergy[0] || zE[i] > tableEnergy[Nenergies-1] ) {
	ce[i] = 0.0;
      } else {
	DoSplint(tableEnergy, tableFlux, splineFlux, zE[i], ce[i]);
      }
    }
    // Simpson five-point formula for integral within the bin
    flux[ie] = ( 14.0*(ce[0]+ce[4])+64.0*(ce[1]+ce[3])+24.0*(ce[2]) )/180.0;

    // Scale by redshift
    flux[ie] -= log10(redshift);
    // Scale by (radius/distance)^2, where radius (in km) and distance (in kpc)
    // PC=3.0856775807e18 cm
    flux[ie] += -32.978701090 + 2.0*log10(parameter[2]/parameter[3]);
    // Convert from ergs/(s cm^2 Hz) to counts/(s cm^2 keV)
    // HH=6.6260693e-27 ergs s
    flux[ie] += 26.178744 - log10(energy[ie+1]);
    flux[ie] = pow(10.0, flux[ie]);
    // Convert counts/(s cm^2 bin)
    flux[ie] *= (energy[ie+1] - energy[ie]);

  }

  return;
}

void DoSpline(const RealArray& x, const RealArray& y, Real& firstDeriv, Real& lastDeriv,
	      RealArray& y2)
{
  size_t n(x.size());
  y2.resize(n);
  RealArray u(n);

  if ( firstDeriv > 0.99e30 ) {
    y2[0] = 0.0;
    u[0] = 0.0;
  } else {
    y2[0] = -0.5;
    u[0] = (3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-firstDeriv);
  }
  for (size_t i=1; i<n-1; i++) {
    Real sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    Real p = sig * y2[i-1] + 2.0;
    y2[i] = (sig - 1.0)/p;
    u[i] = ( 6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/
	     (x[i+1]-x[i-1])-sig*u[i-1] )/p;
  }
  Real qn(0.0), un(0.0);
  if ( lastDeriv <= 0.99e30 ) {
    qn = 0.5;
    un = (3./(x[n-1]-x[n-2]))*(lastDeriv-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.);

  for (int i=n-2; i>=0; i--) {
    y2[i] = y2[i]*y2[i+1] + u[i];
  }

  return;

}


void DoSplint(const RealArray& x, const RealArray& y, const RealArray& y2, 
	      const Real xtarg, Real& yout)
{
  int n(x.size());
  int klo(0), khi(n-1);

  while ( (khi-klo) > 1 ) {
    int k = (khi+klo)/2;
    if ( x[k] > xtarg ) {
      khi = k;
    } else {
      klo = k;
    }
  }

  Real h = x[khi] - x[klo];
  if ( h != 0.0 ) {
    Real a = (x[khi]-xtarg)/h;
    Real b = (xtarg-x[klo])/h;
    yout = a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.;
  } else {
    yout = 0.0;
  }

  return;

}

