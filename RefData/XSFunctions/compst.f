      SUBROUTINE compst(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

c
c **	** subroutine modified for XSPEC
c **	** fwj haberl 28 july 1986
c
C-    TO DERIVE THE FLUX EXPECTED FROM A SUNYAEV-TITARCHUK COMPTONISATIO
C-    SPECTRA BETWEEN E1 AND E2 KEV NORMALISED TO THAT AT EREF KEV
c
c **	** see ADDMOD for parameter descriptions
c	   number of parameters: 2
c		1	temperature
c		2	optical depth
c

      DOUBLE PRECISION PI, RME
      PARAMETER (PI=3.141592654D0, RME=1681.122615D0)

      DOUBLE PRECISION X(3), COM(3)
      DOUBLE PRECISION G1, G2, G3, G4, Y, X1, X2, Z, v, vold
      DOUBLE PRECISION temp, depth, eref, elow, ehi
      INTEGER i, j

      EXTERNAL TAMMA

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors

      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      temp = DBLE(param(1))
      depth = DBLE(param(2))
      eref = 1.d0
      vold = 0.d0

      DO i = 1, ne

         elow = DBLE(ear(i-1))
         ehi = DBLE(ear(i))

         V = temp*depth
         IF ( V. NE. VOLD ) THEN

            Z = -1.5D0 + SQRT(2.25D0+RME/temp/(depth+0.66667D0)
     &          /(depth+0.66667D0))
            CALL TAMMA(Z, G1)
            CALL TAMMA(-3.D0-Z, G2)
            CALL TAMMA(2.D0*Z+4.D0, G3)
            CALL TAMMA(-2.D0*Z-2.D0, G4)

         ENDIF

         X(1) = elow
         X(2) = ehi
         X(3) = eref
         DO j = 1, 3
            Y = X(j)/temp
            IF ( Y .GT. 100 ) THEN
               COM(j) = 0.
            ELSE
               CALL KUMMR(Z, 2.D0*Z+4.D0, X(j), temp, X1)
               CALL KUMMR(-3.D0-Z, -2.D0*Z-2.D0, X(j), temp, X2)
               COM(j) = Y**(2.D0+Z)*DEXP(-Y)*G1*PI/DSIN(2.D0*PI*Z)
     &                  *(X1/G2/G3-Y**(-2.D0*Z-3.D0)*X2/G1/G4)
            ENDIF
         ENDDO
     
         photar(i) = SNGL( (ehi-elow)*(COM(1)+COM(2))/2.D0/COM(3) )

         VOLD = V

      ENDDO

      RETURN
      END
