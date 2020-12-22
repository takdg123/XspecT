#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <xsTypes.h>
#include <cmath>

// function from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
                        const IntegerVector& Zarray, const RealArray& abun, 
                        const Real dens, const Real z, const RealArray& Tarr, 
                        const RealArray& DEMarr, const int ifl, const bool qtherm, 
                        const Real velocity, RealArray& fluxArray, 
                        RealArray& fluxErrArray);



void cemVMekal(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString)
{
   /*  
     c
     c XSPEC model subroutine to calculate Continous Emission
     c Measure of the form :
     c   Q(T) = Norm*(T/Tmax)^alpha  
     c See, for example, Schmitt et al. ApJ 365, 704 (1990), 
     c but note that this program yields Schmitt's alpha - 1.0.
     c
     c This program calls 'vmeka' thus allowing one to vary the elemental
     c abundances.
     c 
     c Parameters:
     c    param(1) = slope of CEM, alpha
     c    param(2) = maximum temperature, tmax
     c    param(3) = nH (cm^-3)  Fixed at 1 for most applications
     c    param(4) = He abundance
     c    param(5) = C   "
     c    param(6) = N   "
     c    param(7) = O   "
     c    param(8) = Ne  "
     c    param(9) = Na  "
     c    param(10)= Mg  "
     c    param(11)= Al  "
     c    param(12)= Si  "
     c    param(13)= S   "
     c    param(14)= Ar  "
     c    param(15)= Ca  "
     c    param(16)= Fe  "
     c    param(17)= Ni  " 
     c    param(18) = redshift used in the Mewe-Kaastra plasma model
     c    param(19) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model)
     c
     c K. P. Singh    April 22, 1994
     c
     c Disclaimer: Any resemblance to a real program is purely
     c             coincidental
     c
       Translated from Fortran cevmkl.f by C. Gordon   Dec. 2007
       Switched to using calcMultiTempPlasma   kaa 12/21/17
   */

   using namespace XSutility;

   const Real k = 8.6171e-8;
   const Real alpha = params[0];
   const Real tMax = params[1];
   const Real max = log10(tMax/k);
   int nt = static_cast<int>((max - 5.5)*10.0 + 1.0);

   RealArray Tarray(nt);
   RealArray demarray(nt);

   /*
      c Integrate contributions in form:
      c    f = (sum over i) Wi * Fi * logdeltT
      c where
      c   Wi= (Ti/Tmax)^alpha and Fi is F(Ti) from meka:
      c
      c it is important to do this in uniform Log(T) steps.
      c
      c I am shuffling the definitions of Fi and photar
      c temporarily for ease.
      c
      c Note: Doing the stepping from logT=log(Tmax) in nt 0.1 steps.
      c       This suggestion from Chris Done to ensure that changing
      c       Tmax by less than 0.1 gives a change in chi-squared.
   */

   for (int i=nt-1; i>=0; --i)
   {
      Real logTemp = max - 0.1*i;
      Tarray[i] = k*pow(10.0, logTemp);
      demarray[i] = 0.1*pow(Tarray[i]/tMax, alpha);      
   }

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
     msg << "cemVMekal: error status " << status << " returned from calcMultiTempPlasma";
     FunctionUtility::xsWrite(msg.str(),5);
   }
}
