#include <xsTypes.h>
#include <functionMap.h>

// routine from tbvabs.cxx file
void tbdefaults(RealArray& params);

void tbgas(const RealArray& energyArray, const RealArray& par,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
   //    Compute X-ray absorptivity of cold gas using the formalism given
   //    by Wilms, Allen, and McCray, ApJ, 2000, in press

   //    Pure gas absorption (no H2, no dust)

   //    The parameters are

   //     par[0]   hydrogen column (1E22/cm^2)
   //     par[1]   redshift

   // Preset other values with their defaults

   RealArray param(.0, 42);
   tbdefaults(param);

   //    H column
   param[0]=par[0];
   // Redshift
   param[41] = par[1];


   // turn off H2 and gas
   param[18] = 0.0;
   param[19] = 0.0;

   // Compute the absorption
   tbvabs(energyArray, param, spectrumNumber, flux, fluxErr, initString);

}
