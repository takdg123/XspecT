#include <xsTypes.h>
#include <functionMap.h>


void equil(const RealArray& energyArray, const RealArray& params,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
   //    March 95  Equlibrium ionization models 
   //		   additive model for XSPEC to use with general models.
   //		   parameters (param):
   //              [0]:    temperature t(keV)
   //		   [1]:    heavy metal abundance
   //		   [2]:	   redshift z

   // Translated from xeq.f 03/09

   RealArray vparam(14);
   vparam[0] = params[0];
   vparam[1] = 1.0;
   for (size_t i=2; i<13; ++i)
      vparam[i] = params[1];
   vparam[13] = params[2];

   // Call the vequil subroutine for variable abundances.

   vequil(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);
}
