
      SUBROUTINE xszbod(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C      black body, uses simple
C      2 point approximation to the integral
C---
C see ADDMOD for parameter descriptions
C number of model parameters: 2
C      1      kT      (keV), must be > 0.
C      2      redshift
C intrinsic energy range:
C      Emin=epsilon(>0),Emax=infinity
C algorithm:
C      n(E) = K 8.0525  E**2 dE / ((kT)**4 (exp(E/kt)-1))
C normalization:
C      The above normalization is such that K = L39 / (D10)**2
C      that is, the norm is the source luminosity (in units of 10^39 ergs/sec)
C      divided by the distance to the source (in units of 10 kpc) squared.
C---
C 18 Feb 1984 - rashafer
C---

      REAL t, anorm, elow, x, tinv, anormh, alow, ehi, ahi,zfac
      INTEGER ie, je

C---

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

      zfac = 1.0+param(2)

      t = param(1)
      tinv = 1./t
      anorm = 8.0525*tinv*tinv*tinv*tinv
      anormh = 0.5*anorm
      elow = ear(0)*zfac
      x = elow*tinv
      IF (x.LE.1.e-4) THEN
         alow = elow*t
      ELSEIF (x.GT.60.) THEN
         DO ie = 1, ne
            photar(ie) = 0.0
         ENDDO
         RETURN
      ELSE
         alow = elow*elow/sngl(dble(exp(x))-1.D0)
      ENDIF
      DO ie = 1, ne
         ehi = ear(ie)*zfac
         x = ehi*tinv
         IF (x.LE.1.e-4) THEN
            ahi = ehi*t
         ELSEIF (x.GT.60.) THEN
            DO je = ie, ne
               photar(je) = 0.0
            ENDDO
            RETURN
         ELSE
            ahi = ehi*ehi/sngl(dble(exp(x))-1.D0)
         ENDIF
         photar(ie) = anormh*(alow+ahi)*(ehi-elow)/zfac
         elow = ehi
         alow = ahi
      ENDDO

      RETURN
      END
