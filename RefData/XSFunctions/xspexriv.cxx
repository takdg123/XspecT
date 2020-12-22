#include "MZCompRefl.h"

// Prototype for the routine which actually does all the work and routine which is
// required for the broken power-law case (bexriv model).

void doIonizedReflection(string ModelName, const RealArray& energyArray, 
			 const RealArray& params, RealArray& flux);
void calcCutoffBrokenPowerLaw(const RealArray& energyArray, const RealArray& params,
			      RealArray& flux);

// Prototype from cutoffPowerLaw.cxx

Real calcCutoffPowerLaw(const RealArray& energyArray, const Real& photIdx, const Real& cutoff, bool isRenorm, RealArray& flux); 

// Prototype from brokenPowerLaw.cxx

Real calcBrokenPowerLaw(const RealArray& energyArray, const RealArray& params, 
                        bool isRenorm, RealArray& flux); 


// Main model routine

void xspexriv(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{
   //  Driver for angle-dependent reflection
   //
   //  See Magdziarz & Zdziarski, 1995, MNRAS.
   //  See Zdziarski et al., 1995, ApJL 438, L63 for description of
   //  calculation of ionization (based on Done et al. 1992, ApJ 395, 275).
   //  The abundances are defined by the command abund
   //
   //  The output spectrum is the sum of the cutoff power law and the
   //  reflection component.
   //  The reflection component alone can be obtained
   //  for scale (see below) = rel_refl < 0. Then the actual
   //  reflection normalization is |scale|. Note that you need to
   //  change then the limits of rel_refl. The range of rel_refl in that case
   //  should exclude zero (as then the direct component appears).
   //
   //  If E_c=0, there is no cutoff in the power law.
   //
   //  number of model parameters:9
   //  1: Gamma, power law photon index, N_E prop. to E^{-Gamma}
   //  2: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
   //     one needs to change the lower limit for that)
   //  3: scale, scaling factor for reflection; if <0, no direct component
   //     (scale=1 for isotropic source above disk)
   //  4: redshift, z
   //  5: abundance of elements heavier than He relative to
   //     the solar abundances
   //  6: iron abundance relative to the solar abundances
   //  7: cosine of inclination angle
   //  8: disk temperature in K
   //  9: disk ionization parameter = L/nR^2
   //  algorithm:
   //       a[E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
   //  Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
   //  of the cutoff power law only (without reflection)
   //  and in the earth frame.
   //
   // also uses the following values which can be set using xset
   //   PEXRIV_PRECISION  fractional precision for Greens' fn. integration
   //  Based on ireflct.cxx

   static bool isFirst = true;
   if (isFirst)
   {
      string msg("Compton reflection from ionized medium.");
      msg += "\nSee help for details.\nIf you use results of this model in a paper,";
      msg += "\nplease refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837";
      xs_write(const_cast<char*>(msg.c_str()), 5);
      isFirst = false;
   }

   doIonizedReflection("PEXRIV", energyArray, params, flux);

}


// *********************************************************************************
// Routine which does all the work for pexriv and bexriv models

void doIonizedReflection(string ModelName, const RealArray& energyArray, const RealArray& params, RealArray& flux)
{
   using namespace XSutility;

   Real Gamma(0.0), Gamma1(0.0), Ebreak(0.0), Gamma2(0.0);

   static int nesave = -1;
   const int Ne = static_cast<int>(energyArray.size()) - 1;
   static RealArray Sptot;
   static RealArray Spref;
   static RealArray Spinc;
   static RealArray X;
   static RealArray zeArray;

   // parameters

   int ipar(0);

   if (ModelName == "PEXRIV" || ModelName == "PEXRAV") {

     Gamma = params[ipar++];

   } else if ( ModelName == "BEXRIV" || ModelName == "BEXRAV") {

     Gamma1 = params[ipar++];
     Ebreak = params[ipar++];
     Gamma2 = params[ipar++];
   }

   Real Ecut = params[ipar++];
   Real Scale = params[ipar++];
   Real zshift = 1.0 + params[ipar++];
   Real Abund = log10(params[ipar++]);
   Real FeAbund  = log10(params[ipar++]);
   Real cosIncl = params[ipar++];
   if ( cosIncl < 0.05 ) cosIncl = 0.05;
   if ( cosIncl > 0.95 ) cosIncl = 0.95;
   Real Temp = params[ipar++];
   Real Xi   = params[ipar++];

   // In order for the absorption and reflection calculations to be performed
   // correctly the energies on which the input spectrum are calculated should
   // run from 0.005 keV to 511.0/YMIN keV.

   int NewBinsHigh = 0;
   int NewBinsLow = 0;
   Real dellogEHigh = 0.01;
   Real dellogELow = 0.01;

   if ( energyArray[0]*zshift > 0.005 ) NewBinsLow = (int) ((log10(energyArray[0]*zshift) - log10(0.005)) * 100.0);

   if ( energyArray[Ne]*zshift < 511.0/YMIN ) NewBinsHigh = (int) ((log10(511.0/YMIN)-log10(energyArray[Ne]*zshift)) * 100.0);

   int NewBins = NewBinsLow + NewBinsHigh;

   // Resize arrays

   if ( Ne+NewBins != nesave ) {
     zeArray.resize(Ne+1+NewBins);
     flux.resize(Ne+NewBins);
     X.resize(Ne+NewBins);
     Sptot.resize(Ne+NewBins);
     Spref.resize(Ne+NewBins);
     Spinc.resize(Ne+NewBins);
     nesave = Ne + NewBins;
   }

   // and set extended redshifted energy array

   for (int i=0; i<NewBinsLow; i++) zeArray[i] = pow(10.0,log10(0.005)+i*dellogELow);
   for (int i=NewBinsLow; i<=Ne+NewBinsLow; i++) zeArray[i] = energyArray[i-NewBinsLow] * zshift;
   for (int i=Ne+NewBinsLow+1; i<(int)zeArray.size(); i++) 
     zeArray[i] = pow(10.0,log10(zeArray[Ne+NewBinsLow])+(i-Ne-NewBinsLow)*dellogEHigh);


   Real factor = 1.0/(2.0*511.0);
   for (int i=0; i<(int)X.size(); i++) X[i] = (zeArray[i+1]+zeArray[i]) * factor;

   // Generate the input spectrum and the normalization factor we will need later

   Real NormFac(1.0);
   RealArray Etmp(2), Ftmp(1);
   Etmp[0] = zshift - 0.01;
   Etmp[1] = zshift + 0.01;

   if ( ModelName == "PEXRIV" || ModelName == "PEXRAV") {

     RealArray inpar(2);
     inpar[0] = Gamma;
     inpar[1] = Ecut;

     calcCutoffPowerLaw(zeArray, Gamma, Ecut, false, flux);

     calcCutoffPowerLaw(Etmp, Gamma, Ecut, false, Ftmp);
     NormFac = Ftmp[0] * zshift*zshift / 0.02;

   } else if ( ModelName == "BEXRIV" || ModelName == "BEXRAV") {

     RealArray inpar(4);
     inpar[0] = Gamma1;
     inpar[1] = Gamma2;
     inpar[2] = Ebreak;
     inpar[3] = Ecut;

     calcCutoffBrokenPowerLaw (zeArray, inpar, flux);

     calcCutoffBrokenPowerLaw (Etmp, inpar, Ftmp);
     NormFac = Ftmp[0]  * zshift*zshift / 0.02;

   }

   // Generate the spinc array from the flux array.

   for (int i=0; i<(int)X.size(); ++i) {
     if (zeArray[i+1] != zeArray[i]) {
       Real avgE = (zeArray[i+1]+zeArray[i])/2.;
       Spinc[i] = flux[i]*avgE*avgE/(zeArray[i+1]-zeArray[i]);
     }
   }

   //     X is the source frame energy array (units m_e c^2)
   //     spinc is the input spectrum array (E F_E)
   //     spref is the reflected spectrum array (E F_E)
   //     sptot is the total spectrum array (E F_E), = spinc if no reflection
   //     all dimensions = Ne

   Real Xmax(zeArray[Ne+NewBinsLow]/511.0);

   // calculate the total flux including direct and Compton reflection (non-relativistic
   // and relativistic).

   calcCompReflTotalFlux(ModelName, Scale, cosIncl, Abund, FeAbund, Xi, Temp, Xmax,
			 X, Spinc, Sptot);

   // Renormalized based on the input spectrum

   Sptot /= NormFac;

   // convert back to xspec internal units

   flux.resize(Ne);
   for (int i=0; i<Ne; ++i)
   {
     Real avgE = (energyArray[i+1]+energyArray[i])/2.;
     flux[i] = Sptot[i+NewBinsLow]*(energyArray[i+1]-energyArray[i])/(avgE*avgE);
   }

  return;

}

#include <XSUtil/Numerics/IncGamma.h>

// **************************************************************************
void calcCutoffBrokenPowerLaw(const RealArray& energyArray, const RealArray& params,
			      RealArray& fluxArray)
{
  // Broken power law with high energy exponential cutoff.  
  // Number of model parameters: 2
  //   1       photIdx1        low E powerlaw photon index
  //   2       photIdx2        high E powerlaw photon index
  //   3       breakE          energy of break (keV)
  //   4       cutoff          energy of exponential cutoff (in
  //                           energy units, e.g. keV)
  // Intrinsic energy range:
  //   Emin = epsilon(>0), Emax = infinity
  //
  // algorithm:
  //   n(E)= E^(-photIdx1) * exp(-E/cutoff) dE       for E < breakE
  //   n(E)= breakE^(photIdx2-photIdx1)*E^(-photIdx2) * exp(-E/cutoff) dE   for E < breakE
  //   This relies on an approximate incomplete gamma function
  //   calculation to perform the integral over n(E).  
  //   WARNING: The approximation loses accuracy in the region,
  //   10^-6 > (1-photIdx) > 0.

  const size_t nBins = energyArray.size()-1;
  fluxArray.resize(nBins);

  const Real& photIdx1 = params[0];
  const Real& photIdx2 = params[1];
  const Real& breakE = params[2];
  const Real& cutoff = params[3];

  fluxArray = 0.0;

  // trap out case of cutoff <= 0. in this case just do a simple power-law

  if ( cutoff <= 0.0 ) {
    calcBrokenPowerLaw(energyArray, params, false, fluxArray);
    return;
  }

  // continue with non-zero cutoff...

  Numerics::IncGamma incGamma;
  Real a;
  if ( energyArray[0] < breakE ) {
    a = 1.0 - photIdx1;
  } else {
    a = 1.0 - photIdx2;
  }
  Real x = energyArray[0]/cutoff;
  Real multiplier = pow(cutoff, a);
  Real lowIntegral = incGamma(a,x);
  for (size_t i=0; i<nBins; ++i) {
    if ( energyArray[i] < breakE && energyArray[i+1] >= breakE ) {
      x = breakE/cutoff;
      Real highIntegral = incGamma(a,x);
      fluxArray[i] = multiplier*(lowIntegral - highIntegral);
      a = 1.0 - photIdx2;
      multiplier = pow(cutoff, a) * pow(breakE,(photIdx2-photIdx1));
      lowIntegral = incGamma(a,x);
    }
    x = energyArray[i+1]/cutoff;
    Real highIntegral = incGamma(a,x);
    fluxArray[i] += multiplier*(lowIntegral - highIntegral);
    lowIntegral = highIntegral;
  }

  return;
}
