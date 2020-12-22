// Converted from Fortran in xsgabs.f and modified to average exp(-tau) over
// the energy bin, not tau.
// number of model parameters:3
//       0       EL      line energy (in energy units, e.g. keV)
//       1       W       line width (sigma) (in energy units)
//       2       Tau     line optical depth
// intrinsic energy range: none
// algorithm:
//       M(E) = exp(-(tau/sqrt(2*PI)/W) * exp(-0.5*((E-EL)/W)**2)))
//       If W is less than or equal to zero, the line is treated as a delta
//            function.
//       N.B. when the energy spacing is much greater than W and the stepsize
//            for EL then the partial derivative determinations can lead
//            to fit error conditions.  The solution is to increase the
//            stepsize for EL


#include <XSFunctions/functionMap.h>
#include <xsTypes.h>
#include <XSstreams.h>
#include <cmath>
#include <XSUtil/Numerics/AdaptiveIntegrate.h>
#include <XSUtil/Numerics/BinarySearch.h>

using namespace Numerics;

RealArray gaussAbsIntegrand(const RealArray& x, void *p);

void gaussianAbsorptionLine(const RealArray& energyArray, const RealArray& params, 
			    int spectrumNumber, RealArray& fluxArray, 
			    RealArray& fluxErrArray, const string& initString)
{
  int nE = energyArray.size();
  fluxArray.resize(nE-1);
  fluxErrArray.resize(0);
  //fluxErrArray.resize(nE-1);
  fluxArray = 1.0;
  //fluxErrArray = 0.0;

  Real efirst = energyArray[0];
  Real elast = energyArray[nE-1];

  // Set the line energy, width and the energy range over which the
  // line will be calculated (essentially +/- crtsig sigma)

  Real crtsig = 6.0;
  Real emin = params[0] - crtsig*params[1];
  Real emax = params[0] + crtsig*params[1];

  // If the line calculation energies lie outside the energy range
  // then do not have to do anything

  if ( emax < efirst || emin > elast ) return;

  // find the first and last bins on which we need to calculate the absorption line

  int ifirst(0), ilast(nE-2);
  if ( emin >= efirst ) {
    ifirst = BinarySearch(energyArray, emin);
  }
  if ( emax <= elast ) {
    ilast = BinarySearch(energyArray, emax);
  }

  // if the width is zero so the entire line is in one bin

  if ( params[1] == 0.0 ) {
    fluxArray[ifirst] = 1.0;
    return;
  }

  Real Precision = 1.0e-6;
  RealArray p(params);

  int ModelEvaluations(0);
  for (int i=ifirst; i<=ilast; i++) {

    Real Integral, IntegralError;
    ModelEvaluations += 
      AdaptiveIntegrate<gaussAbsIntegrand>(energyArray[i], energyArray[i+1],
					   &p, Precision, Integral, IntegralError);
    fluxArray[i] = Integral;
    //fluxErrArray[i] = IntegralError;
    Real ebin = energyArray[i+1]-energyArray[i];
    if ( ebin > 0.0 ) {
      fluxArray[i] /= ebin;
     // fluxErrArray[i] /= ebin;
    }

  }

  return;
}

// evaluate the gaussian absorption line integrand with p used to pass in the
// parameter values

RealArray gaussAbsIntegrand(const RealArray& x, void *p)
{
  RealArray OutArray(x.size());

  RealArray params(3);
  params = *(static_cast<RealArray*>(p));

  Real center = params[0];
  Real width  = params[1];
  Real tau    = params[2];

  Real gnorm = 1.0/width/sqrt(2*M_PI);

  for (size_t i=0; i<x.size(); i++) {
    Real gauss = gnorm * exp(-0.5*(x[i]-center)*(x[i]-center)/width/width);
    OutArray[i] = exp(-tau * gauss);
  }

  return OutArray;
}
