
      SUBROUTINE xsabsc(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(1), photar(ne), photer(ne)

C---
C XSPEC model subroutine:
C simple exponential low energy cutoff for absorption
C---
C            see MULMOD for parameter descriptions
C      number of model parameters: 1
C            1      E0: The e-folding energy of the absorption
C                  nominal values:
C      no intrinsic energy range limit
C      model form A(E) = exp(-E0/E)
C---
C 10 nov 1983 - rashafer
C---
      REAL e0n, a, b
      INTEGER i
C---

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      e0n = -param(1)
      a = exp(e0n/ear(0))
      DO i = 1, ne
         b = exp(e0n/ear(i))
         photar(i) = 0.5*(a+b)
         a = b
      ENDDO
      RETURN
      END
