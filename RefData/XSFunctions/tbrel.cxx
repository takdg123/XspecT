#include <xsTypes.h>
#include <functionMap.h>

// routine from tbvabs.cxx file
void tbdefaults(RealArray& params);

void tbrel(const RealArray& energyArray, const RealArray& par,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
   //    Compute X-ray absorptivity of cold gas using the formalism given
   //    by Wilms, Allen, and McCray, ApJ, 2000, in press

   //    Allows a relative column. Parameters are same as for tbvabs
   //    except that if the H column is < 0 then the output is 1/absorption
   //    instead of absorption.

   RealArray param(0.0,42);

   int signh = (par[0]<0)? -1: +1;
   param[0] = fabs(par[0]);
   for (size_t i=1; i<param.size(); i++) param[i] = par[i];

   // Compute the absorption
   tbvabs(energyArray, param, spectrumNumber, flux, fluxErr, initString);

   if ( signh < 0 ) {
     for (size_t i=0; i<flux.size(); i++) flux[i] = 1.0/flux[i];
   }

}
