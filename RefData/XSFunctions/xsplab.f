      SUBROUTINE xsplab(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

C---
C simple absorption as a power-law in energy - useful for things
C like dust. calculates mean transmission over each energy interval
C---
C number of parameters: 2
C      1     power law index: 'index', no range limit
C            nominal values 2.0, 0.0, 5.0, 0.001
C      2     coefficient
C no intrinsic energy range limit
C model form A(E) = coef * E**(-index)
C---
C kaa  2/1/95
C---
      DOUBLE PRECISION a, b, alphan, alphni, ae, be
      INTEGER i
C---
c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      alphan = 1.d0 - DBLE(param(1))

C  special case of index=1, then must use logarithmic integral

      IF (alphan .EQ. 0.d0) THEN

         be = DBLE(ear(0))
         b = LOG(be)
         DO i = 1, ne
            ae = DBLE(ear(i))
            a = LOG(ae)
            IF ( ae .NE. be ) THEN
               photar(i) = param(2)*SNGL((a - b)/(ae - be))
            ELSE
               photar(i) = 1.
            ENDIF
            b = a
            be = ae
         ENDDO

C  boring old standard integral

      ELSE

         alphni = 1.d0/alphan
         be = DBLE(ear(0))
         b = alphni*(be**alphan)
         DO i = 1, ne
            ae = DBLE(ear(i))
            a = alphni*(ae**alphan)
            IF ( ae .NE. be ) THEN
               photar(i) = param(2)*SNGL((a - b)/(ae-be))
            ELSE
               photar(i) = 1.
            ENDIF
            be = ae
            b = a
         ENDDO

      ENDIF

      RETURN
      END



