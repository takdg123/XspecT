// gaussian model with redshift in wavelength space using Angstrom units
// just calls agauss

#include <xsTypes.h>
#include <functionMap.h>

void zagauss(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& fluxArray, 
	     RealArray& fluxErrArray, const string& initString)
{
  // create energy array in source frame
  Real zfac = 1.0 + params[2];
  RealArray zenergyArray = energyArray*zfac;

  // and parameter array to pass to agauss
  RealArray zparams(params.size()-1);
  for (size_t i=0; i<zparams.size(); i++) zparams[i] = params[i];

  agauss(zenergyArray, zparams, spectrumNumber, fluxArray, fluxErrArray, initString);

  // include time dilation factor

  fluxArray /= zfac;

}
