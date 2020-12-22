      SUBROUTINE xscabs(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(1), photar(ne), photer(ne)

C---
C XSPEC model subroutine to calculate absorption factor due
C to Compton scattering. Includes Klein-Nishina X-sections
C for higher energies but not energy downshifting
C---
C number of model parameters: 1
C      1      ANH      Hydrogen column density (in units of 10**22
C                      atoms per square centimeter
C---
C Arguments :
C    ear      r         i: energy ranges
c    ne       i         i: number of energy ranges
c    param    r         i: model parameter values
c    ifl      i         i: the dataset
c    photar   r         r: the transmission fraction
C---

      REAL energy

      INTEGER ie

      REAL sigmakn
      EXTERNAL sigmakn

c suppress a warning message from the compiler
      ie = ifl

      DO ie = 1, ne
         energy = 0.5*(ear(ie-1)+ear(ie))
         photar(ie) = EXP(-param(1)*1.21*0.665e-2*sigmakn(energy))
         photer(ie) = 0.0
      ENDDO

      RETURN
      END


      REAL FUNCTION sigmakn(x)

      REAL erest, x
      DOUBLE PRECISION w, z1, z2, z3, z4, z5, z6

c erest is 511 keV, x is the photon energy in keV, sigmakn is
c the Klein-Nishina cross section in units of the Thomson one

      erest = 510.999
      w = x/erest

      IF ( w .LE. 0.025 ) THEN

c use the asymptotic limit for w less than 0.025
c the relative error at w<0.025 is <5E-8

         sigmakn = SNGL( 1.0d0
     &        +w*(-140+w*(364+w*(-931+w*(2288+w*(-5440)))))/70.0d0 )

      ELSE

         z1 = (1+w)/w**3
         z2 = 1+2*w
         z3 = LOG(z2)
         z4 = 2*w*(1+w)/z2
         z5 = z3/2/w
         z6 = (1+3*w)/z2/z2
         sigmakn = SNGL(0.75d0*(z1*(z4-z3)+z5-z6))

      ENDIF

      END
