
      SUBROUTINE DISKM(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(4), photar(ne), photer(ne)

c
c **	** subroutine for XSPEC
c **	** modified for use on the VAX fwj haberl 15 july 1986
c
c **	**LUIGI TYPE DISK WHERE ALPHA SCALES AS ONLY GAS PRESSURE
c
c **	** see ADDMOD for parameter description
c	   number of parameters: 4
c		1	accretion rate in Eddington units
c		2	neutron star mass in solar mass units
c		3	inner disk radius in Schwarzschild radius units
c		4	alpha viscosity
c
c  The normalization outside this routine gives the distance to the
c  source, but also depends on inclination
c           [d/cosi]**2=1.0/norm where d is in units of 10 kpc.
c
      REAL amdot, amass, rin, routl, rinl, tcon, acons, tcons
      REAL dr, r1, r1l, e, r2, r2l, r, dx, anorm, t, edelt
      REAL alpha, acon, adiskm, rp5, bbm
      INTEGER nnn, i
C---
c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors

      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      amdot = param(1)
      AMASS = param(2)
      RIN = param(3)
      ALPHA = param(4)
c
      ACON = 0.77*AMASS*AMASS/(AMASS**0.73*ALPHA**0.178)
      TCON = 9.91*ALPHA**0.178/AMASS**0.2667
      ROUTL = 2.5
      RINL = ALOG10(RIN)
      NNN = 20
      DR = (ROUTL-RINL)/NNN
      ACONS = ACON*AMDOT**0.470
      TCONS = TCON*AMDOT**0.533
      DO i = 1, ne
         R1L = RINL - DR*0.5
         R1 = 10**R1L
         adiskm = 0.
         DO WHILE (r1l.LE.routl)
            e = (ear(i)+ear(i-1))/2.
            R2L = R1L + DR
            R2 = 10.0**R2L
            DX = R2 - R1
            R = (R1+R2)*0.5
            IF (R.GT.1.01) THEN
               ANORM = 1.0/R**1.80
               RP5 = 1.0 - 1.0/R**0.5
               ANORM = ANORM*RP5**0.467
               T = RP5**0.533/R**1.2
               T = T*TCONS
               aDISKM = aDISKM + BBM(E, T)*R*DX*ANORM
            ENDIF
            R1 = R2
            R1L = R2L
         ENDDO
         aDISKM = aDISKM*2.0*ACONS
         edelt = ear(i) - ear(i-1)
         photar(i) = edelt*adiskm
      ENDDO
      RETURN
      END
