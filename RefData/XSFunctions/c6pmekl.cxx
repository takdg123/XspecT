#include <functionMap.h>
#include <xsTypes.h>

/*
  c XSPEC model subroutine to calculate Differential Emission
  c Measure of the form:
  c   Q(T) = Norm*(exp(w(T))), where w(T) = Sum over 6 orders of Chebyshev
  c                                        polynomials with 6 coeffs. 
  c See, for example, Lemen et al. ApJ 341, 474 (1989), 
  c The DEM is constrained to be positive here, as the exponential 
  c of the Chebyshev Polynomial is being used as in Lemen et al. 
  c Parameters:
  c    param(1) = coeff a1 of Chebyshev polynomial order 1 
  c    param(2) = coeff a2 of Chebyshev polynomial order 2 
  c    param(3) = coeff a3 of Chebyshev polynomial order 3 
  c    param(4) = coeff a4 of Chebyshev polynomial order 4 
  c    param(5) = coeff a5 of Chebyshev polynomial order 5 
  c    param(6) = coeff a6 of Chebyshev polynomial order 6 
  c    param(7) = nH (cm^-3)  Fixed at 1 for most applications
  c    param(8) = metallicity used in the Mewe-Kaastra plasma model
  c    param(9) = redshift used in the Mewe-Kaastra plasma model
  c    param(10) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model)
  c
  c K. P. Singh    March 23, 1995 
  c                June 1, 1995 Changed for XSPEC 8.70
  c                Jan 13, 1996 CP6 defined properly within the limits
  c                             of orthogonality and Changed for XSPEC 9.0
  c kaa  12/21/17  Converted to C++
  c Disclaimer: Any resemblance to a real program is purely
  c             coincidental
*/

void c6pmekl(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)  
{
  RealArray pparams(23);

  for (size_t i=0; i<7; i++) pparams[i] = params[i];
  pparams[7] = 1.0;
  for (size_t i=8; i<21; i++) pparams[i] = params[7];
  pparams[21] = params[8];
  pparams[22] = params[9];

  c6pvmkl(energyArray, pparams, spectrumNumber, flux, fluxErr, initString);
}

