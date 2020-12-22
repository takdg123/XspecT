      SUBROUTINE xsntch(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(3), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C Notch line absorption, equivalent to a very saturated line.
C---
C see ADDMOD for parameter descriptions
C number of model parameters:3
C      1      EL      line energy (in energy units, e.g. keV)
C      2      EQW     line equivalent wiflh (in energy units)
C      3      CF      covering fraction
C intrinsic energy range: none
C      N.B. the line may have significant effects at unphysical(i.e. negative)
C      energies when EL<~W.
C algorithm:
C      The multiplicative factor f = (1-CF) over the range from EL +/-
C      EQW/2.  For energy ranges that are only part covered,
C      by the range of the line, the absorption is scaled appropriately.
C N.B. when the energy spacing is much greater than EQW and the stepsize
C      for EL then the partial derivative determinations can lead
C      to fit error conditions.  The solution is to increase the
C      stepsize for EL
C---
C 17 Dec 1983 - rashafer
C---
      REAL el, eqw, cf, cfc, emin, emax, del
      INTEGER ie, ielow

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO


C---
      el = param(1)
      eqw = max(0., param(2))
      cf = param(3)
      cfc = 1. - cf
      emin = el - eqw*.5
      emax = el + eqw*.5
      DO ie = 1, ne
         photar(ie) = 1.
      ENDDO
      IF ((emax.GT.ear(0)) .AND. (emin.LT.ear(ne))) THEN
         ielow = 0
         DO WHILE ((ear(ielow).LE.emin))
            ielow = ielow + 1
         ENDDO
         IF ((ielow.NE.0)) THEN
            del = (emin-ear(ielow-1))/(ear(ielow)-ear(ielow-1))
            photar(ielow) = del + cfc*(1.-del)
         ELSE
            photar(1) = cfc
            ielow = 1
         ENDIF
         IF (ear(ielow).LT.emax) THEN
            ielow = ielow + 1
            DO WHILE ((ielow.LE.ne).AND.(ear(ielow).LT.emax))
               photar(ielow) = cfc
               ielow = ielow + 1
            ENDDO
            IF (ielow.LE.ne) photar(ielow) = cfc
         ENDIF
         IF (ielow.LE.ne) THEN
C **        ** restore the part of the notch from emin to the top of the bin
            del = (ear(ielow)-emax)/(ear(ielow)-ear(ielow-1))
            photar(ielow) = photar(ielow) + cf*del
         ENDIF
      ENDIF
      RETURN
      END
