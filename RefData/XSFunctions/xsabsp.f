      SUBROUTINE xsabsp(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C XSPEC model subroutine for "partial covering absorption".
C Assumes that part of the emitter is covered by the given
C absorption and the rest of the emitter is unobscured.
C---
C number of model parameters: 2
C      1      ANH    Hydrogen column density (in units of 10**22
C                    atoms per square centimeter
C      2      FRACT  Covering fraction (0 implies no absorption,
C                    1 implies the emitter is all absorbed with
C                    the indicated column ANH.
C---
C Arguments :
C    ear      r         i: energy ranges
c    ne       i         i: number of energy ranges
c    param    r         i: model parameter values
c    ifl      i         i: the dataset
c    photar   r         r: the transmission fraction
C---

      REAL fract, fractc
      INTEGER ie

C calculate the standard photo-electric absorption

      CALL xsphab(ear, ne, param(1), ifl, photar, photer)

C now modify for the partial covering

      fract = param(2)
      fractc = 1. - fract

      DO ie = 1, ne
         photar(ie) = fractc + fract * photar(ie)
      ENDDO

      RETURN
      END
