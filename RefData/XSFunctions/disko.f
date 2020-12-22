      SUBROUTINE DISKO(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(4), photar(ne), photer(ne)

c
c **	** subroutine for XSPEC
c **	** modified for use on VAX fwj haberl 11 july 1986
c
c **	** DISK MODEL USING S&S DISK +MBB
c
c **	** see ADDMOD for parameter description
c	   number of parameters: 4
c		1	accretion rate in eddington units
c		2	neutron star mass in solar mass units
c		3	inner disk radius in Schwarzschild radius units
c		4	alpha viscosity
c
      REAL amdot, amass, rin, routl, rinl, tcon, acons, tcons
      REAL dr, r1, r1l, e, r2, r2l, r, dx, anorm, t, edelt
      REAL alpha, acon, adisko, rp5, bbm
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
      IF (ALPHA.EQ.2.) ALPHA = 1.0/10.0**5.0
C
      ACON = 0.14*AMASS*AMASS/(AMASS**0.78*ALPHA**0.22)
      TCON = 127.0*ALPHA**0.22/AMASS**0.22
C
      ROUTL = 2.5
      RINL = ALOG10(RIN)
      NNN = 15
      DR = (ROUTL-RINL)/NNN
      ACONS = ACON*AMDOT**0.11
      TCONS = TCON*AMDOT**0.889
      DO i = 1, ne
         R1L = RINL - DR*0.5
         R1 = 10.0**R1L
         adisko = 0.
         DO WHILE (r1l.LE.routl)
            e = (ear(i)+ear(i-1))/2.
            R2L = R1L + DR
            R2 = 10.0**R2L
            DX = R2 - R1
            R = (R2+R1)/2.0
            IF (R.GT.1.001) THEN
               ANORM = 1.0/R**1.33
               RP5 = 1.0 - 1.0/R**0.5
               ANORM = ANORM*RP5**0.111
               T = TCONS*RP5**0.889/R**1.667
               aDISKO = aDISKO + BBM(E, T)*R*DX*ANORM
            ENDIF
            R1L = R2L
            R1 = R2
         ENDDO
         aDISKO = aDISKO*ACONS
         edelt = ear(i) - ear(i-1)
         photar(i) = edelt*adisko
      ENDDO
      RETURN
      END
