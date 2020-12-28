#include <xsTypes.h>
#include <XSUtil/FunctionUtils/funcType.h>
#include <XSFunctions/functionMap.h>
#include <XSUtil/Numerics/IncGamma.h>

// prototype for routine which does the work

Real calcCutoffPowerLaw(const RealArray& energyArray, const Real& photIdx, const Real& cutoff, bool isRenorm, RealArray& fluxArray); 

// prototypes for routine from powerLaw.cxx

Real calcPowerLaw(const RealArray& energyArray, const Real& index, bool isRenorm, RealArray& fluxArray); 
Real pegRenormalize(const RealArray& energyArray, RealArray& fluxArray);


void 
cutoffPowerLaw (const RealArray& energyArray, 
                const RealArray& params, 
                int spectrumNumber,
                RealArray& fluxArray, 
                RealArray& fluxErrArray,
                const string& initString)
{
   // Power law with high energy exponential cutoff.  
   // Number of model parameters: 2
   //   1       photIdx         powerlaw photon index
   //   2       cutoff          energy of exponential cutoff (in
   //                           energy units, e.g. keV). if <= 0
   //                           then no cutoff applied.
   // Intrinsic energy range:
   //   Emin = epsilon(>0), Emax = infinity
   //
   // algorithm:
   //   n(E)= E**(-photIdx) * exp(-E/cutoff) dE
   //   This relies on an approximate incomplete gamma function
   //   calculation to perform the integral over n(E).  
   //   WARNING: The approximation loses accuracy in the region,
   //   10^-6 > (1+photIdx) > 0.

   const Real& photIdx = params[0];
   const Real& epeak = params[1];

   calcCutoffPowerLaw(energyArray, photIdx, epeak, true, fluxArray); 
   fluxErrArray.resize(0);

}


// ***************************************************************************
// routine which does the actual calculation

Real calcCutoffPowerLaw(const RealArray& energyArray, const Real& photIdx, const Real& epeak, bool isRenorm, RealArray& fluxArray)
{

   const size_t nBins = energyArray.size()-1;
   fluxArray.resize(nBins);

   Real cutoff = epeak/(2.0+photIdx);
   
   if ( cutoff <= 0.0 ) {

     calcPowerLaw(energyArray, photIdx, false, fluxArray)*pow(100, -photIdx);

   } else {

     Numerics::IncGamma incGamma;
     Real a = 1.0 + photIdx;
     Real x = energyArray[0]/cutoff;
     const Real multiplier = pow(cutoff, a);
     Real lowIntegral = incGamma(a,x);
     for (size_t i=1; i<=nBins; ++i)
       {
	 x = energyArray[i]/cutoff;
	 Real highIntegral = incGamma(a,x);
	 fluxArray[i-1] = multiplier*(lowIntegral - highIntegral)*pow(100, -photIdx);
	 lowIntegral = highIntegral;
       }

   }

   if ( isRenorm ) {

     return(pegRenormalize(energyArray, fluxArray));

   } else {

     return(1.0);

   }

}

