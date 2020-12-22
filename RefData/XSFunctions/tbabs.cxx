#include <xsTypes.h>
#include <functionMap.h>

// routine from tbvabs.cxx file
void tbdefaults(RealArray& params);

void tbabs(const RealArray& energyArray, const RealArray& par,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   //    Compute X-ray absorptivity of cold gas using the formalism given
   //    by Wilms, Allen, and McCray, ApJ, 2000, in press

   //    This version has all parameters fixed at their default values

   //    The parameter is

   //    par[0]   hydrogen column (1E22/cm^2)

   //    This model assumes that 20% of the hydrogen is molecular
   //    and that grains with the depletions given in the above paper
   //    are present (see subroutine fggrain).

   RealArray param(0.0, 42);
   tbdefaults(param);
   param[0] = par[0];

   // Compute the absorption
   tbvabs(energyArray, param, spectrumNumber, flux, fluxErr, initString);


}
