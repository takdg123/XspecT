// Subroutine to smooth the model spectrum by relativistic effects from a
// disk using the Brenneman & Reynolds spin model. 
// ACF and RMJ May/June 1998
//
// Updated by LB, Feb. 2006
//  Updated Mar 2009 - renamed from spinconv to dospinconv, tmp array
//       is now allocated in outer c++ wrapper spinconv.cxx (CG)
// Updated to use blurring and remove fixed 20 eV resolution.
// Updated to use FFTs  kaa 8/19
//
// Parameters
//       0        power law index for emissivity for inner disk
//       1        power law index for emissivity for outer disk
//       2        break radius (GM/c**2)
//       3        BH spin (dimensionless)
//       4        inclination  (degrees)
//       5        inner radius of disk (units of marginal stability)
//       6        outer radius of disk (units of marginal stability)

#include <xsTypes.h>
#include <functionMap.h>
#include <XSUtil/Numerics/Convolution.h>

void spinconv(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{

  RealArray spinParams(9);

  // set the line energy to the midpoint of the energy range (this is abitrary)

  spinParams[0] = (energyArray[energyArray.size()-1]+energyArray[0])/2.0;
  for (size_t i=0; i<7; i++) spinParams[i+1] = params[i];
  spinParams[8] = 0.0;

  Numerics::ConvolutionInLnSpace<spin>(energyArray, spinParams, spinParams[0], 
				       spectrumNumber, "", flux, fluxErr);
  return;
}
