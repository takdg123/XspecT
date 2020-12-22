      SUBROUTINE xsbrmv(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C simple thermal bremsstrahlung, based on GSFC routine
C FBG, as derived by GSFC for low temperature (<20 keV) plasmas
C the number of photons in the bin is estimated by a simple 2
C point approximation to the integral. Includes variable He/H.
C---
C see ADDMOD for parameter descriptions
C number of model parameters: 1
C      1	T(keV)	Plasma temperature in keV
C      2	n(He)/n(H)  Helium to Hydrogen ratio
C intrinsic energy range:
C      Emine=epsilon, Emax=infinity
C algorithm:
C      Uses the function Gaunt, which returns the photon spectrum,
C      based on the GSFC routine FBG
C---
C 28 July 1988 - kaa
C 16 Dec 1994 - kaa         switched to Gaunt routine and included
C                           1/sqrt(T) factor to make the normalization
C                           independent of T. Normalization is now :
C                             K = (3.02e-15/4/pi/R^2) Int n_e n_H dV
C                                 
C---
      REAL t, ab, elow, alow, ehi, ahi, tfac
      INTEGER ie

      REAL gaunt
      EXTERNAL gaunt

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

      t = param(1)
      tfac = 1./SQRT(t)

      ab = param(2)

      elow = ear(0)

      alow = 0.
      IF ( (t.GT.0) .AND. (elow.LT.50*t) .AND. (elow.GT.0.) )
     &    alow = (gaunt(elow, t, 1.) + 4*ab*gaunt(elow, t, 2.))
     &          *tfac*EXP(-elow/t)/elow

      DO ie = 1, ne

         ehi = ear(ie)

         ahi = 0.
         IF ( (t.GT.0) .AND. (ehi.LT.50*t) .AND. (ehi.GT.0.) )
     &      ahi = (gaunt(ehi, t, 1.) + 4*ab*gaunt(ehi, t, 2.))
     &             *tfac*EXP(-ehi/t)/ehi

         photar(ie) = 0.5*(alow+ahi)*(ehi-elow)

         elow = ehi
         alow = ahi

      ENDDO

      RETURN
      END
