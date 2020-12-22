#include <xsTypes.h>
#include <functionMap.h>

// routine from tbvabs.cxx file
void tbdefaults(RealArray& params);

void ztbabs(const RealArray& energyArray, const RealArray& par,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   //    Compute X-ray absorptivity of cold gas using the formalism given
   //    by Wilms, Allen, and McCray, ApJ, 2000, in press

   //    This version is redshift-dependent and includes no grains

   //    The parameters are

   //    par[0]   hydrogen column (1E22/cm^2)
   //    par[1]   redshift

   //    This model assumes that 20% of the hydrogen is molecular
   //    and that there is NO MATERIAL IN GRAINS.

   //    Version 1.0, 2000/06/16, Joern Wilms, wilms@astro.uni-tuebingen.de

   RealArray param(0.0, 42);
   tbdefaults(param);

   param[0] = par[0];
   param[41] = par[1];

   // turn off grains
   param[19]=0.0;

   // Compute the absorption
   tbvabs(energyArray, param, spectrumNumber, flux, fluxErr, initString);

}
