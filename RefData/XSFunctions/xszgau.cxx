#include <xsTypes.h>
#include <functionMap.h>
#include <stlToCArrays.h>
#include <XSUtil/Utils/XSutility.h>
#include <cfortran.h>

void xszgau(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // ---
   //  XSPEC model subroutine
   //  Simple gaussian line profile
   // ---
   //  see ADDMOD for parameter descriptions
   //  number of model parameters:3
   //        1       EL      line energy (in energy units, e.g. keV)
   //        2       W	line wiflh (sigma) (in energy units)
   //        3 	     redshift
   //  intrinsic energy range: none
   //        N.B. the line may have significant flux at unphysical (i.e. negative)
   //             energies when EL<~W.
   //  algorithm:
   //        n(E) = (1./(sqrt(2pi*W*W)) * exp(-0.5*((E-EL)/W)**2)) dE
   //        N.B. the norm is equal to the total integrated number of counts in
   //             the line.
   //        If W is less than or equal to zero, the line is treated as a delta
   //             function.
   //        N.B. when the energy spacing is much greater than W and the stepsize
   //             for EL then the partial derivative determinations can lead
   //             to fit error conditions.  The solution is to increase the
   //             stepsize for EL
   // ---
   //  Just calls gaussian line model with redshifted energies.

   Real zfac = 1.0 + params[2];
   RealArray Ez(energyArray*zfac);
   RealArray gparams(2);
   gparams[0] = params[0];
   gparams[1] = params[1];

   gaussianLine(Ez, gparams, spectrumNumber, flux, fluxErr, initString);

   // Include time dilation factor

   flux /= zfac;       
}
