#include <xsTypes.h>
#include <functionMap.h>

// routine from tbvabs.cxx file
void tbdefaults(RealArray& params);

void tbpcf(const RealArray& energyArray, const RealArray& par,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
   //    Compute X-ray absorptivity of cold gas using the formalism given
   //    by Wilms, Allen, and McCray, ApJ, 2000, in press

   //    Partial covering version with all other parameters set to defaults

   //    The parameters are

   //     par[0]   hydrogen column (1E22/cm^2)
   //     par[1]   covering fraction
   //     par[2]   redshift

   // Preset other values with their defaults

   RealArray param(.0, 42);
   tbdefaults(param);

   //    H column
   param[0]=par[0];
   // Redshift
   param[41] = par[2];

   // Compute the absorption
   tbvabs(energyArray, param, spectrumNumber, flux, fluxErr, initString);

   // Partial covering
   for (size_t i=0; i<flux.size(); i++) flux[i] = 1.0 + par[1]*(flux[i]-1.0);

}
