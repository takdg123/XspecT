// This a version of the simple NEI model using the Aped
// class.

#include <xsTypes.h>
#include <functionMap.h>

void vnei(const RealArray& energyArray, const RealArray& params,
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
   //		(2..13): abundances of He,C,N,O,Ne,Mg,Si,S,Ar,Ca,Fe,Ni
   //                                  with respect to solar values
   //           (14):   ionization timescale tau (cm^-3 s^-1)
   //		(15):	redshift z
   //
   //           Parameter (1) should be usually set to one. Its only use
   //           in the code is to set hydrogen abundance to zero; in this
   //           case this parameter should be set to zero. It is not the
   //           density parameter.
   //


  // vnei is just a special case of vgnei with the final and ionization-timescale
  // averaged temperatures equal.

  RealArray tparams(17);
  for (size_t i=0; i<15; i++) tparams[i] = params[i];
  tparams[15] = params[0];
  tparams[16] = params[15];

  vgnei(energyArray, tparams, spectrumNumber, flux, fluxErr, initString);

  return;

}
