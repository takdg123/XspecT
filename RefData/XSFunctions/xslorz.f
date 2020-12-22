
      SUBROUTINE xslorz (ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

c XSPEC model subroutine to calculate a Lorentzian line profile.
c Parameters are :
c         1        Line energy in keV
c         2        Line width in keV.
c
c The Lorentzian is given by :
c
c       L(E) = A (W/2) / ( (E-E0)**2 + (W/2)**2 )
c
c where E0 is the line energy and W the line width. This integrates
c analytically to give the contribution in the energy range (E1,E2) as
c
c       L(E1,E2) = A ( arctan(2(E2-E0)/W) - arctan(2(E1-E0)/W) )
c
c Note that if W=0 then L=0. The factor A is set so that L(0,inf) = 1
c since the energy cannot be negative =>
c       A = 1 / ((pi/2) - arctan(-2E0/W))

      REAL f1, f2, a
      INTEGER ie

c suppress a warning message from the compiler
      ie = ifl

c this model has no errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

c Trap the case of W=0

      IF ( param(2) .LE. 0 ) THEN
         DO ie = 1, ne
            photar(ie) = 0.
         ENDDO
         RETURN
      ENDIF

c Set the normalization factor

      a = 1.0 / ( asin(1.0) - atan(-2.0*param(1)/param(2)) )

c Loop round the energies

      f2 = atan(2*(ear(0)-param(1))/param(2))
      sum = 0.0

      DO ie = 1, ne

         f1 = f2
         f2 = atan(2*(ear(ie)-param(1))/param(2))

         photar(ie) = (f2 - f1) * a

      ENDDO

      RETURN
      END
