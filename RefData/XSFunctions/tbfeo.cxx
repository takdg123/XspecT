#include <xsTypes.h>
#include <functionMap.h>

// routine from tbvabs.cxx file
void tbdefaults(RealArray& params);

void tbfeo(const RealArray& energyArray, const RealArray& par,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
   //    Compute X-ray absorptivity of cold gas using the formalism given
   //    by Wilms, Allen, and McCray, ApJ, 2000, in press

   //    This version allows fitting the nH, and the O and Fe abundances.
   //    All other parameters are left at their default values

   //    The parameters are

   //     par[0]   hydrogen column (1E22/cm^2)
   //     par[1]   O abundance relative to Solar
   //     par[2]   Fe abundance relative to Solar
   //     par[3]   redshift

   //     Depletion et al. are still at their default values

   // Preset other values with their defaults

   RealArray param(.0, 42);
   tbdefaults(param);

   //    H column
   param[0]=par[0];
   //    oxygen      [7]
   param[4]=par[1];
   //    iron        [25]
   param[15]=par[2];
   // Redshift
   param[41] = par[3];

   // Compute the absorption
   tbvabs(energyArray, param, spectrumNumber, flux, fluxErr, initString);

}
