
      SUBROUTINE xszbrm(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C XSPEC model subroutine
C simple thermal bremsstrahlung, based on GSFC routine
C FBG, as derived by GSFC for low temperature (<20 keV) plasmas
C the number of photons in the bin is estimated by a simple 2
C point approximation to the integral
C---
C see ADDMOD for parameter descriptions
C number of model parameters: 2
C      1      T(keV)      Plasma temperature in keV
C      2      redshift
C intrinsic energy range:
C      Emine=epsilon, Emax=infinity
C---
C 18 May 1984 - rashafer
C 16 Dec 1994 - kaa           simply calls XSBRMV
C---

      REAL par(2), zfac
      INTEGER ie


c Shift the energy array

      zfac = 1.00+param(2)
      DO ie = 0, ne
         ear(ie) = ear(ie)*zfac
      ENDDO

c Calculate the model

      par(1) = param(1)
      par(2) = 0.085

      CALL xsbrmv(ear, ne, par, ifl, photar, photer)

c Shift the energies back and correct the model

      DO ie = 0, ne
         ear(ie) = ear(ie)/zfac
      ENDDO

      DO ie = 1, ne
         photar(ie) = photar(ie)/zfac
      ENDDO

      RETURN
      END



