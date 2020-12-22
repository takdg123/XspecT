#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <xsTypes.h>
#include <cmath>

// function from calcMultiTempPlasma.cxx
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType, 
                        const IntegerVector& Zarray, const RealArray& abun, 
                        const Real dens, const Real z, const RealArray& Tarr, 
                        const RealArray& DEMarr, const int ifl, const bool qtherm, 
                        const Real velocity, RealArray& fluxArray, 
                        RealArray& fluxErrArray);


/*
  c XSPEC model subroutine to calculate Differential Emission
  c Measure of the form:
  c   Q(T) = Norm*(exp(w(T))), where w(T) = Sum over 6 orders of Chebyshev
  c                                        polynomials with 6 coeffs. 
  C CAUTION : The DEM here is constrained to be positive throughout,
  C by taking the exponential of the sum of the 
  C polynomial and which follows Lemen et al. ApJ 341, 474 (1989), 
  c 
  c Parameters:
  c    param(1) = coeff a1 of Chebyshev polynomial order 1 
  c    param(2) = coeff a2 of Chebyshev polynomial order 2 
  c    param(3) = coeff a3 of Chebyshev polynomial order 3 
  c    param(4) = coeff a4 of Chebyshev polynomial order 4 
  c    param(5) = coeff a5 of Chebyshev polynomial order 5 
  c    param(6) = coeff a6 of Chebyshev polynomial order 6 
  c    param(7) = nH (cm^-3)  Fixed at 1 for most applications
  c    param(8) = He abundance
  c    param(9) = C   "
  c    param(10) = N   "
  c    param(11) = O   "
  c    param(12) = Ne  "
  c    param(13) = Na  "
  c    param(14)= Mg  "
  c    param(15)= Al  "
  c    param(16)= Si  "
  c    param(17)= S   "
  c    param(18)= Ar  "
  c    param(19)= Ca  "
  c    param(20)= Fe  "
  c    param(21)= Ni  "
  c    param(22) = redshift used in the Mewe-Kaastra plasma model
  c    param(23) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model)
  c
  c K. P. Singh    September 15, 1995 
  c                Jan. 13, 1996  Modified for XSPEC V 9.0
  c                               Modified for orthogonality condition
  c kaa   12/21/17 converted to C++
  c
  c Disclaimer: Any resemblance to a real program is purely
  c             coincidental
*/

void c6pvmkl(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)  
{

  using namespace XSutility;
  using namespace Numerics;

  // Initialize values:

  const Real k(8.6171e-8);
  RealArray ak(6);
  for (size_t i=0; i<6; i++) ak[i] = params[i];

  // Integrate contributions in form:
  //    f = (sum over j) Q * Fj * logdeltT
  // where
  //   Q= (sum over k) ak(Pk)
  // where Pk is Chebyshev polynomial of order k,
  // and Fj is F(Tj) from meka for diff. energies or wavelengths
  //
  // it is important to do sum over j in uniform Log(T) steps.
  //
  // Note: Doing the stepping from logT=5.5 to 8.0 in 0.1 steps.
  //
  //    tarr is in units of keV
  //    temp is the temperature scaled with a constant=1E-6

  size_t Ntemp(25);
  RealArray Tarr(Ntemp);
  RealArray demarr(Ntemp);

  for (size_t it=0; it<Ntemp; it++) {

    Real logtemp = 5.4 + 0.1*(it+1);
    Tarr[it] = pow(10,logtemp) * k;

    Real temp  = (1E-6)*pow(10,logtemp);
    Real X = ((logtemp - 5.5)*0.8 - 1.0);
    RealArray pk(6);
    pk[0] = X;
    pk[1] = 2*X*X - 1;
    for(size_t k=2; k<6; k++) pk[k] = 2*X*pk[k-1] - pk[k-2];
    Real Q(0.0);
    for(size_t k=0; k<6; k++) Q += ak[k]*pk[k];
    demarr[it] = exp(Q) * temp * 0.1;

  }

   // set up all the variables to pass to calcMultiTempPlasma

   int swtch = static_cast<int>(params[22]);
   int plasmaType(6);
   if ( swtch == 0 ) {
     plasmaType = 3;
   } else if ( swtch == 1 ) {
     plasmaType = 4;
   } else if ( swtch == 2 ) {
     plasmaType = 6;
   }

   const Real density = params[6];
   const Real redshift = params[21];

   RealArray abun(14);
   for (size_t i=0; i<abun.size(); i++) abun[i] = params[i+7];
   const int elements[] = {2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 26, 28};
   IntegerVector Zarray(14);
   for (size_t i=0; i<14; i++) Zarray[i] = elements[i];

   const bool qtherm = false;
   const Real velocity = 0.0;

   calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun, density,
		       redshift, Tarr, demarr, spectrumNumber, 
		       qtherm, velocity, flux, fluxErr);

}
