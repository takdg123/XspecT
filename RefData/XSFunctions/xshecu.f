      SUBROUTINE xshecu(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C XSPEC model subroutine:
C simple exponential high energy cutoff
C---
C see MULMOD for parameter descriptions
C number of model parameters: 2
C      1      E0: The cutoff energy
C      2      E1: normalization of the exponent
C model form A(E) = exp((E0-E)/E1)
C---
C 25 june 1986 - fwj haberl
C---
      REAL e0n, e1n, a, b
      INTEGER i
C---
c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      e0n = param(1)
      e1n = param(2)
      a = exp(-amax1((ear(0)-e0n),.001)/e1n)
      DO i = 1, ne
         b = exp(-amax1((ear(i)-e0n),.001)/e1n)
         photar(i) = 0.5*(a+b)
         a = b
      ENDDO
      RETURN
      END
