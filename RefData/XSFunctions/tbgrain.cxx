#include <xsTypes.h>
#include <functionMap.h>

// routine from tbvabs.cxx file
void tbdefaults(RealArray& params);

void tbgrain(const RealArray& energyArray, const RealArray& par,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   //    Compute X-ray absorptivity of cold gas using the formalism given
   //    by Wilms, Allen, and McCray, ApJ, 2000, in press

   //    This version allows to vary the grain and h2 parameters

   //    The parameters are

   //     par[0]   hydrogen column (1E22/cm^2)
   //     par[1]   fraction of h in h2 (default: 0.2)
   //     par[2]   density of dust (default: 1 g/cm^3)
   //     par[3]   minimum thickness, mum (0.025)
   //     par[4]   maximum thickness, mum (0.25)
   //     par[5]   PL index with implicit minus sign (3.5)

   //     Depletion et al. are still at their default values


   RealArray param(.0, 42);
   tbdefaults(param);

   param[0] = par[0];
   for (size_t i=1; i<6; i++) param[i+17] = par[i];

   // Compute the absorption
   tbvabs(energyArray, param, spectrumNumber, flux, fluxErr, initString);

}
