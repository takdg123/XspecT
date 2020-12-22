
      SUBROUTINE xsgaul(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C Simple gaussian line profile
C---
C see ADDMOD for parameter descriptions
C number of model parameters:2
C       1       EL      line energy (in energy units, e.g. keV)
C       2       W	line width (sigma) (in energy units)
C intrinsic energy range: none
C       N.B. the line may have significant flux at unphysical (i.e. negative)
C            energies when EL<~W.
C algorithm:
C       n(E) = (1./(sqrt(2pi*W*W)) * exp(-0.5*((E-EL)/W)**2)) dE
C       N.B. the norm is equal to the total integrated number of counts in
C            the line.
C       If W is less than or equal to zero, the line is treated as a delta
C            function.
C       N.B. when the energy spacing is much greater than W and the stepsize
C            for EL then the partial derivative determinations can lead
C            to fit error conditions.  The solution is to increase the
C            stepsize for EL
C---
C 17 Dec 1983 - rashafer
C 1.1: An attempt to use simple linear interpolation for the
C      response to better handle very narrow lines.  This
C      is not strictly correct, and effectively the lines will have
C      widths of order of the ear spacing at that point.
C 4/19/98  kaa  Reverted to putting whole line in a single bin if the
C               width is small. This means for lines with width < ear
C               spacing the energy can only be determined to the ear
C               spacing. This change is necessary to determine line widths
C               correctly using high resolution spectroscopy.
C---

      call gaussianline(ear, ne, param, ifl, photar, photer)

      return
      end

