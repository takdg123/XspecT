      subroutine XSGRBM(ear, ne, param, ifl, photar, photer)
C
      IMPLICIT NONE
      INTEGER          ne, ifl, ie
      REAL             ear(0:ne), param(3), photar(ne), photer(ne)
      DOUBLE PRECISION nprime, alfa, beta, amb, cvt,
     +                 tem, temi, e_b, e_p, betap1, 
     +                 elo, vlo, ehi, vhi, foo, output
      PARAMETER        (e_p=1.d+2)
C
C-->  Approximation to the Gamma-Ray Burst model by D. Band, et al.,
C     1993 (ApJ 413, 281) for XSPEC programs:
C     N(E) = n  * (E/E_p)**alfa * exp(-E/tem), if E < (alfa-beta)tem
C          = n' * (E/E_p)**beta              , if E > (alfa-beta)tem
C
C     The parameters in order are n, alfa, beta, and tem.
C     Suggested default values and increments are (note that
C     the normalization "n" is handled outside this function):
C         n    = 0.01        Increment = 1e-4
C         alfa = -1.0        Increment = 1e-2
C         beta = -2.0        Increment = 1e-2
C         tem  = 300.        Increment = 10.
C
C-->  Revisions:
C     Written and revised for BSAS by D. Band and S. Bansal (1/27/92; 
C     7/06/92; 3/02/93);  Adapted for TDAS by S. Bansal (4/25/94); 
C     Modified for XSPEC by H. Seifert (7/23/96);
C

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

C-->  Initialize variables...
      alfa = param(1)
      beta = param(2)
      cvt  = alfa + 2.d0
      if (cvt .le. 0) cvt = 0.0001
      tem  = param(3)/cvt
C
C-->  Define some secondary variables...
      temi = 1.d0/tem
C
      amb = alfa-beta
      if (amb .lt. temi) amb = temi
      e_b = amb*tem
C
      nprime = ((e_b/e_p)**amb)*exp(-amb)
      betap1 = beta+1.d0

C
C-->  Get first data point and consider the cases
C     when E is >/< (alfa-beta)tem; finally, also   
C     consider also the special case of beta=-1.
      elo = ear(0)
      vlo = 0.0d0
      if (elo .le. e_b) then
         vlo = ((elo/e_p)**alfa)*exp(-elo*temi)
      else
         if (beta .eq. -1.d0) then
            vlo = nprime*e_p*LOG(elo)
         else
            vlo = nprime*e_p*((elo/e_p)**betap1)/betap1
         endif
      endif

C
C-->  Loop over all energy bins, considering again
C     the cases when E is >/< (alfa-beta)tem and, again,
C     also consider the special case of beta=-1.
C     Note that if both elo and ehi are >(alfa-beta)tem,
C     we can do the exact integral; otherwise, we do a
C     simple two-point approximation.
      do ie = 1, ne
         ehi = ear(ie)
         if (ehi .le. e_b) then
            vhi = ((ehi/e_p)**alfa)*exp(-ehi*temi)
         else
            if (elo .gt. e_b) then
               if (beta .eq. -1.d0) then
                  vhi = nprime*e_p*LOG(ehi)
               else
                  vhi = nprime*e_p*((ehi/e_p)**betap1)/betap1
               endif
            else
               vhi = nprime*(ehi/e_p)**beta
               foo = nprime*e_p*((ehi/e_p)**betap1)/betap1
            endif
         endif
C-->     Fill output array... 
         if (elo .gt. e_b) then
            output = vhi-vlo  
         else
            output = (vhi+vlo)*(ehi-elo)/2.d0
         endif

         IF ( output .LT. 2.0e-38 ) THEN
            photar(ie) = 0.0
         ELSE IF ( output .GT. 2.0e38 ) THEN
            photar(ie) = 2.0e38
         ELSE
            photar(ie) = SNGL(output)
         ENDIF

C-->     Set up for next loop...
         if ((ehi .gt. e_b) .and. (elo .le. e_b)) then
            vlo = foo
         else
            vlo = vhi
         endif
         elo = ehi
      enddo
C
      return
      end

