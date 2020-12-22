
      SUBROUTINE xsxpdec(ear, ne, param, ifl, photar, photer)

      IMPLICIT NONE

      INTEGER ne, ifl
      REAL ear(0:ne), param(1), photar(ne), photer(ne)

C---
C additive exponential decay model for XSPEC  exp(-param(1)*E)
C---
C see ADDMOD for parameter descriptions
C number of parameters: 1
C      1      exponential factor
C no intrinsic energy range limit
C---
C  01/05/03   kaa
C---
      REAL factor, explow, exphi
      INTEGER i
C---

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      DO i = 1, ne
         photar(i) = 0.
      ENDDO

c  loop over energies

      IF ( param(1) .EQ. 0.0 ) RETURN

      factor = 1./param(1)

      explow = exp(-param(1)*ear(0))

      DO i = 1, ne

         exphi = exp(-param(1)*ear(i))
         photar(i) = (explow - exphi) * factor
         explow = exphi

      ENDDO

      RETURN
      END
