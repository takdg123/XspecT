// This a version of the simple NEI model using the Aped class.

#include <xsTypes.h>
#include <functionMap.h>

void vvnei(const RealArray& energyArray, const RealArray& params,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
   //
   //  aug 95   Non-equlibrium ionization models with constant electron
   //           temperature
   //		additive model for XSPEC to use with general models.
   //		parameters (params):
   //           (0):    temperature t(keV)
   //           (1):    hydrogen abundance (switches on and off free-free
   //                   continuum
   //		(2..30): abundances of all elements from Z=2 to Z=30
   //                                  with respect to solar values
   //           (31):   ionization timescale tau (cm^-3 s^-1)
   //		(32):	redshift z
   //
   //           Parameter (1) should be usually set to one. Its only use
   //           in the code is to set hydrogen abundance to zero; in this
   //           case this parameter should be set to zero. It is not the
   //           density parameter.
   //


  // vvnei is just a special case of vvgnei with the final and ionization-timescale
  // averaged temperatures equal.

  RealArray tparams(34);
  for (size_t i=0; i<32; i++) tparams[i] = params[i];
  tparams[32] = params[0];
  tparams[33] = params[32];

  vvgnei(energyArray, tparams, spectrumNumber, flux, fluxErr, initString);

  return;

}
