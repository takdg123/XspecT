// Subroutine to smooth the model spectrum by relativistic effects from a
// disk in the presence of a non-spinning black hole - uses laor.
// ACF and RMJ May/June 1998
// kaa converted to C++ Aug 2018
// kaa changed to use FFTs Aug 2019

//  parameters :
//       0        power law index for emissivity (10 for disk)
//       1        inner radius (GM/c**2)
//       2        outer radius (GM/c**2)
//       3        inclination  (degrees)

#include <xsTypes.h>
#include <functionMap.h>
#include <XSUtil/Numerics/Convolution.h>


void kdblur(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	    const string& initString)
{

  RealArray laorParams(5);

  // set the laor energy to the midpoint of the energy range (this is arbitrary)
 
  laorParams[0] = (energyArray[energyArray.size()-1]+energyArray[0])/2.0;
  for (size_t i=0; i<4; i++) laorParams[i+1] = params[i];

  Numerics::ConvolutionInLnSpace<laor>(energyArray, laorParams, laorParams[0], 
				       spectrumNumber, "", flux, fluxErr);

  return;
  
}
