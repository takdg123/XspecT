#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <sstream>

// function definition from calcCoolingFlow.cxx
void calcCoolingFlow(const RealArray& energyArray, const Real tlow, 
		     const Real thigh, const Real slope, const IntegerVector& Zarray,
		     const RealArray& abun, const Real z, const int plasmaType, const int ifl,
		     RealArray& flux, RealArray& fluxErr, const Real tpeak = -1.0);



void xscflw(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{

  //  XSPEC model subroutine to calculate cooling flow spectrum
  //  Parameters :
  //        1..................emission measure distribution
  //        2..................low temperature
  //        3..................high temperature
  //        4..................abundance
  //        5..................redshift

  //  Norm is mass accretion rate in units of Msun/yr

  //  kaa  3/24/88
  //  kaa 12/21/92    attempt to speed up and fix intermittent errors.
  //                  NB the workphot array is used to store the R-S data
  //                  as read in and not after rebinning to the ear array
  //                  as in other routines.
  //                  Following a couple of bug fixes from daw this now gives
  //                  results for slope=0 consistent with rmj's model.
  //  kaa   9/3/94    modified for FITS format files
  //  kaa   9/6/96    use of dynamic memory - shares loadrs routine with R-S model
  //  kaa  10/3/96    replaced loadrs by more general ldpfil call
  //  kaa  10/4/96    now just does a call to XSVCFL
  //  cg    2/6/09    Front end and dynamic memory allocation translated
  //                  to C++.
  //  kaa  12/20/17   Modified to call calcCoolingFlow

   // Set abundances - assume He to be Solar
   const int elements[] = {2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28};
   IntegerVector Zarray(12);
   for (size_t i=0; i<12; i++) Zarray[i] = elements[i];

   Real slope = params[0];
   Real tlow = params[1];
   Real thigh = params[2];

   Real z = params[4];
   if ( z <= 0.0 ) {
     FunctionUtility::xsWrite("\n XSCFLW: Require z > 0 for cooling flow models",10);
     return;
   }

   RealArray abun(12);
   abun[0] = 1.0;
   for (int i=1; i<12; ++i) abun[i] = params[3];

   int plasmaType = 1;

   calcCoolingFlow(energyArray, tlow, thigh, slope, Zarray, abun, z, plasmaType,
		   spectrumNumber, flux, fluxErr);

}
