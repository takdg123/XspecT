#include <functionMap.h>
#include <xsTypes.h>

void xsvmcf(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString);


void xsmkcf(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{

   //  XSPEC model subroutine to calculate cooling flow spectrum
   //  from sum of MEKAL spectra.
   //  Parameters :
   //        1..................low temperature
   //        2..................high temperature
   //        3..................abundance
   //        4..................redshift
   //        5..................switch(0=calculate MEKAL model, 
   //                                  1=interpolate MEKAL model)

   //  Norm is mass accretion rate in units of Msun/yr

   //  kaa 8/3/93      based on XSCFLW.

   // cg 2/9/09   Translated to C++

   RealArray vparam(18);

   // Set the parameters for the call to the subroutine with variable
   // abundances (fix He to Solar).

   vparam[0] = params[0];
   vparam[1] = params[1];
   vparam[2] = 1.0;
   for (size_t i=3; i<16; ++i)
      vparam[i] = params[2];
   vparam[16] = params[3];
   vparam[17] = params[4];

   // Call the MEKA cooling flow routine with variable abundances

   xsvmcf(energyArray, vparam, spectrumNumber, flux, fluxErr, initString);

}
