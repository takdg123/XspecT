      subroutine XSSTEP(ear, ne, param, ifl, photar, photer)
C
      IMPLICIT NONE
      INTEGER          ne, ifl, ie
      REAL             ear(0:ne), param(2), photar(ne), photer(ne)
      DOUBLE PRECISION alfa, cvt, e_p,
     +                 tem, temi, 
     +                 elo, vlo, ehi, vhi, output
      PARAMETER        (e_p=1.d+2)
C
C--> Cutoff powerlaw
C

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

C-->  Initialize variables...
      alfa = -param(1)
      cvt  = alfa + 2.d0
      if (cvt .le. 0) cvt = 0.0001
      tem  = param(2)/cvt
C
C-->  Define some secondary variables...
      temi = 1.d0/tem
C
C-->  Get first data point and consider the cases
C     when E is >/< (alfa-beta)tem; finally, also   
C     consider also the special case of beta=-1.
      elo = ear(0)
      vlo = (elo/e_p)**alfa*exp(-elo*temi)
C
C-->  Loop over all energy bins, considering again
C     the cases when E is >/< (alfa-beta)tem and, again,
C     also consider the special case of beta=-1.
C     Note that if both elo and ehi are >(alfa-beta)tem,
C     we can do the exact integral; otherwise, we do a
C     simple two-point approximation.
      do ie = 1, ne
         ehi = ear(ie)
         vhi = (ehi/e_p)**alfa*exp(-ehi*temi)
C-->     Fill output array... 
         output = (vhi+vlo)*(ehi-elo)/2.d0
         IF ( output .LT. 2.0e-38 ) THEN
            photar(ie) = 0.0
         ELSE IF ( output .GT. 2.0e38 ) THEN
            photar(ie) = 2.0e38
         ELSE
            photar(ie) = SNGL(output)
         ENDIF

C-->     Set up for next loop...
         elo = ehi
         vlo = vhi
      enddo
C
      return
      end
