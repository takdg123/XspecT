      SUBROUTINE xswnab(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C XSPEC model subroutine for "window absorption".
C There is no absorption below the energy specified in the
C model parameter.
C---
C number of model parameters: 2
C      1      ANH      Hydrogen column density (in units of 10**22
C                      atoms per square centimeter
C      2      WINDOWE  Low energy boundary of absorption
C---
C Arguments :
C    ear      r         i: energy ranges
c    ne       i         i: number of energy ranges
c    param    r         i: model parameter values
c    ifl      i         i: the dataset
c    photar   r         r: the transmission fraction
C---

      REAL e0
      INTEGER ie

C calculate the standard photo-electric absorption

      CALL xsphab(ear, ne, param(1), ifl, photar, photer)

C now modify to remove any absorption below the window energy

      e0 = param(2)

      DO ie = 1, ne
         IF ( ear(ie) .LE. e0 ) THEN
            photar(ie) = 1.
         ELSEIF ( (ear(ie) .GT. e0) .AND. (ear(ie-1) .LT. e0) ) THEN
            photar(ie) = ( (e0-ear(ie-1)) + (ear(ie)-e0)*photar(ie) ) /
     &                   (ear(ie)-ear(ie-1))
         ENDIF
      ENDDO

      RETURN
      END
