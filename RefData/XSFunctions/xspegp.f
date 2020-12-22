
      SUBROUTINE xspegp(ear,ne,param,ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne),param(3),photar(ne), photer(ne)

C---
C simple single power law spectral model for XSPEC programs
C modified form XSPWLW, so that the norm is adjusted to be
C proportional to the flux of the model over a specified range.
C---
C see ADDMOD for parameter descriptions
C number of parameters: 3
C      1      PHOTON power law index: 'index', no range limit
C             nominal values 1.5, 0.0, 5.0, 0.001
C      2      EMIN lower peg energy range
C      3      EMAX upper peg energy range
C             N.B., EMIN and EMAX should NOT be fitted parameters
C no intrinsic energy range limit
C The norm is equal to the flux, in units of 1.e-12 ergs/cmsq/sec over
C the indicated range, OR if EMIN=EMAX, it is in micro-Janskys at the
C point EMIN.  If EMIN or EMAX is <= 0, then the norm is just the
C flux at 1 keV in photons/cmsq/sec/keV (as in the standard POWERLAW
C model).
C model form N(E) = E**(-index)
C---
C 14 aug 1983 - rashafer
C---

      REAL alphan, emin, emax, anorm, alphnp, b, a, alphni
      INTEGER i

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      alphan=1.-param(1)
      emin=param(2)
      emax=param(3)

c If either EMIN or EMAX is <= 0 then set norm to be flux at 1 keV

      IF((emin.LE.0.).OR.(emax.LE.0.)) THEN

         anorm=1.

c else if EMIN=EMAX then make the norm proportional to the flux in 
c micro-Janskys at the peg point.

      ELSEIF(emin.EQ.emax) THEN

            anorm=1.509E-3/(emin**alphan)

c else make the norm proportional to the flux (in units of 1.e-12 ergs) 
c over the indicated range

      ELSE

         alphnp=alphan+1.
         IF(alphnp.EQ.0.) THEN
            anorm=1./((log(emax)-log(emin))*1.602e3)
         ELSE
            anorm=alphnp/(((emax**alphnp)-(emin**alphnp))*1.602e3)
         ENDIF

      ENDIF

c Now do the power-law integration. Note the special case if index=1.

      IF (alphan.EQ.0) THEN
         b=alog(ear(0))
         DO i=1,ne
            a=alog(ear(i))
            photar(i)=anorm*(a-b)
            b=a
         ENDDO
      ELSE
         alphni=anorm/alphan
         b=alphni*(ear(0)**alphan)
         DO i=1,ne
            a=alphni*(ear(i)**alphan)
            photar(i)=a-b
            b=a
         ENDDO
      ENDIF

      RETURN
      END
