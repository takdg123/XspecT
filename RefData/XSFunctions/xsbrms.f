      SUBROUTINE xsbrms(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(1), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C simple thermal bremsstrahlung, based on GSFC routine
C FBG, as derived by GSFC for low temperature (<20 keV) plasmas
C the number of photons in the bin is estimated by a simple 2
C point approximation to the integral. assumes that n(He)/n(H)=0.085.
C---
C see ADDMOD for parameter descriptions
C number of model parameters: 1
C      1	T(keV)	Plasma temperature in keV
C intrinsic energy range:
C      Emine=epsilon, Emax=infinity
C---
C 18 May 1984 - rashafer
C 16 Dec 1994 - kaa        changed to call XSBRMV.       
C---

      REAL    par(2)

      par(1) = param(1)
      par(2) = 0.085

      CALL xsbrmv(ear, ne, par, ifl, photar, photer)

      RETURN
      END
