      SUBROUTINE xsdust(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(*), photar(ne), photer(ne)

C---
C Dust halo scattering model. Assumes a 1/E**2 dependence
C for scattering fraction and a halo of uniform surface
C brightness with a 1/E dependence.
C Parameters are :
C      1      FRAC      Scattering fraction at 1 keV
C      2      HALOSZ    Size of halo at 1 keV in detector beamsizes
C---
C EAR     I     The array of model energies
C NE      I     The number of model energies
C PARAM   I     The model parameters
C PHOTAR  I     Multiplicative correction to input spectrum
C               for scattering out of beam
C---

      REAL e, invsze
      INTEGER ie

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors

      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO


C Loop round energies

      DO ie = 1, ne
         e = (ear(ie-1)+ear(ie))/2.
         invsze = e/param(2)
         IF (invsze .GE. 1.) THEN
            photar(ie) = 1.
         ELSE
            photar(ie) = 1. - param(1)*(1-invsze**2)/e**2
         ENDIF
      ENDDO

      RETURN
      END


