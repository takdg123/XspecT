      SUBROUTINE xscnst(ear,ne,param,ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne),param(1),photar(ne), photer(ne)

C---
C multiplicative model for XSPEC  -  multiplies by an energy independent factor
C---
C see MULMOD for parameter descriptions
C number of parameters: 1
C      1      multiplicative factor
C no intrinsic energy range limit
C---
C  Feb 24, 1994  kaa
C---

      INTEGER i

c suppress warning messages from the compiler
      i = ifl
      i = INT(ear(0))


      DO i = 1, ne
         photar(i) = param(1)
         photer(i) = 0.0
      ENDDO

      RETURN
      END
