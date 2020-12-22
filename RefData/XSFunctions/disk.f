      SUBROUTINE DISK(ear, ne, param, ifl, photar, photer)

      IMPLICIT NONE

      INTEGER ne, ifl
      REAL ear(0:ne), param(3), photar(ne), photer(ne)

c **	** subroutine modified for XSPEC
c **	** fwj haberl 10 july 1986
c
c **	** STANDARD OPTICALLY THICK DISK
c
c **	** see ADDMOD for parameter descriptions
c	   number of parameters: 3
c		1	accretion rate
c		2	neutron star mass in solar mass units
c		3	inner disk radius in 3*Schwarzschild radius units

c 8/24/95  - factor of two bug in normalization fixed (new)

      REAL amdot, amass, rin, dist, routl, rinl, tcon, acons, tcons
      REAL dr, r1, r1l, adisk, e, r2, r2l, r, dx, anorm, t, edelt
      INTEGER nnn, i

      REAL bbhalm
      EXTERNAL bbhalm

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors

      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      amdot = param(1)
      AMASS = param(2)
      RIN = param(3)

c with the distance d fixed to 10 kpc the norm is then cos(i)/d**2
c with d in units of 10kpc

      DIST = 10.0

c normalization to a distance of 10kpc

      DIST = DIST/10.0
      TCON = 2.54/(AMASS**0.25)
      ACONS = 0.088*AMASS*AMASS/(DIST*DIST)

      ROUTL = 2.5
      RINL = ALOG10(RIN)
      NNN = 18
      DR = (ROUTL-RINL)/NNN

      TCONS = TCON*AMDOT**0.25

      DO i = 1, ne
         R1L = RINL - DR*0.5
         R1 = 10.0**R1L
         adisk = 0.
         DO WHILE (r1l.LE.routl)
            e = (ear(i)+ear(i-1))/2.
            R2L = R1L + DR
            R2 = 10.0**R2L
            DX = R2 - R1
            R = (R2+R1)*0.5
            IF (R.GT.1.01) THEN
               T = TCONS/R**0.75
               T = T*(1.0-1.0/R**0.5)**0.25
               ANORM = T*T*T
               aDISK = aDISK + BBHALM(E, T)*R*DX*ANORM
            ENDIF
            R1L = R2L
            R1 = R2
         ENDDO
         aDISK = aDISK*ACONS
         edelt = ear(i) - ear(i-1)
         photar(i) = edelt*adisk
      ENDDO

      RETURN
      END

c*********************************
      REAL FUNCTION BBHALM(E, T)
      
C **  VERSION 12.11.99, I. HALM
      
      REAL E, T, X, EX
      
      X = E/T
      BBHALM = X*X*X
      EX = EXP(-X)
      BBHALM = BBHALM*EX/(1-EX)
      BBHALM = BBHALM/E
      RETURN
      END
