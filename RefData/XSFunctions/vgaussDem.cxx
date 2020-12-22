#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <xsTypes.h>
#include <cmath>
#include <iostream>
#include <sstream>

// function from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
                        const IntegerVector& Zarray, const RealArray& abun, 
                        const Real dens, const Real z, const RealArray& Tarr, 
                        const RealArray& DEMarr, const int ifl, const bool qtherm, 
                        const Real velocity, RealArray& fluxArray, 
                        RealArray& fluxErrArray);



// XSPEC model subroutine to calculate collisional plasma with a gaussian DEM
// 
// Parameters:
//    param(1) = Temperature mean
//    param(2) = Temperature sigma
//    param(3) = nH (cm^-3)  Fixed at 1 for most applications
//    param(4) = He abundance
//    param(5) = C   "
//    param(6) = N   "
//    param(7) = O   "
//    param(8) = Ne  "
//    param(9) = Na  "
//    param(10)= Mg  "
//    param(11)= Al  "
//    param(12)= Si  "
//    param(13)= S   "
//    param(14)= Ar  "
//    param(15)= Ca  "
//    param(16)= Fe  "
//    param(17)= Ni  " 
//    param(18)= redshift
//    param(19) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model,
//                       2=AtomDB model)


void vgaussDem(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString)
{


   using namespace XSutility;
   using namespace Numerics;

   const Real Tmean = params[0];
   const Real Tsigma = params[1];

   // *******************************************************************
   // set up arrays of temperature and DEM values
   // use nT temperature steps running from -nSig sigma to +nSig sigma

   int nT = 21;
   int nSig = 3.0;

   RealArray Tarray(nT);
   RealArray demarray(nT);

   Real Tmin = Tmean - nSig*Tsigma;
   if ( Tmin < 0.0 ) Tmin = 0.0;
   Real Tdelta = 2 * nSig * Tsigma / nT;

   for (int i=0; i<nT; i++) {
     Real T1 = Tmin + Tdelta * i;
     Real T2 = T1 + Tdelta;
     Tarray[i] = 0.5 * (T1 + T2);
     demarray[i] = erf((T2-Tmean)/Tsigma) - erf((T1-Tmean)/Tsigma);
   }

   // end of set up arrays of temperature and DEM values
   // *******************************************************************


   // set up all the variables to pass to calcMultiTempPlasma

   int swtch = static_cast<int>(params[18]);
   int plasmaType(6);
   if ( swtch == 0 ) {
     plasmaType = 3;
   } else if ( swtch == 1 ) {
     plasmaType = 4;
   } else if ( swtch == 2 ) {
     plasmaType = 6;
   }

   const Real density = params[2];
   const Real redshift = params[17];

   RealArray abun(14);
   for (size_t i=0; i<abun.size(); i++) abun[i] = params[i+3];
   const int elements[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};
   IntegerVector Zarray(14);
   for (size_t i=0; i<14; i++) Zarray[i] = elements[i];

   const bool qtherm = false;
   const Real velocity = 0.0;

   int status=0;
   status = calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, density,
                                redshift, Tarray, demarray, spectrumNumber, 
				qtherm, velocity, flux, fluxErr);

   if (status != 0) {
     std::ostringstream msg;
     msg << "vgaussDem: error status " << status << " returned from calcMultiTempPlasma";
     FunctionUtility::xsWrite(msg.str(),5);
   }

}
