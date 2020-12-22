// Convolution model that can be used to create a parameter which is the 
// luminosity in a particular energy range for the input redshift
// Parameters are     energ_lo     Source frame low energy over which to 
//                                 calculate luminosity
//                    energ_hi     Source frame high energy over which to 
//                                 calculate luminosity
//                    z            Source redshift
//                    lumin        Luminosity (10^44 erg/s)

#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/CosmologyFunction.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

void clumin (const RealArray& energyArray, 
	     const RealArray& params, 
	     int spectrumNumber,
	     RealArray& fluxArray, 
	     RealArray& fluxErrArray,
	     const string& initString)
{

  using namespace std;
  using namespace Numerics;

  static const Real LUMCON (1.07057e61);
  Numerics::FZSQ fzsq;
  Real H0 = FunctionUtility::getH0();
  Real q0 = FunctionUtility::getq0();
  Real lambda0 = FunctionUtility::getlambda0();

  Real redshift (params[2]);
  Real eMin (params[0]/(1.0+redshift));
  Real eMax (params[1]/(1.0+redshift));
  Real lumin (pow(10.0, params[3]));

  // Integrate the input array between elow and ehi into fluxsum.

  pair<Real,Real> fluxsum (integrationKernel(energyArray,fluxArray,eMin,eMax));

  Real luminsum = fluxsum.second * LUMCON * fzsq(redshift,q0,lambda0)/H0/H0;

  // Scale by lumin/luminsum

  fluxArray *= lumin/luminsum;
  fluxErrArray *= lumin/luminsum;

  return;

}
