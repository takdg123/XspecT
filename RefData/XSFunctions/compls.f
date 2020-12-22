      SUBROUTINE compls(ear,ne,param,ifl, photar, photer)
      IMPLICIT NONE

      INTEGER ne, ifl
      REAL ear(0:ne),param(2),photar(ne), photer(ne)

c
c **      ** subroutine modified for XSPEC
c **      ** fwj haberl 28 july 1986
c
c **      ** CALCULATES A COMPTONIZATION SPECTRUM AFTER LAMB ANS SANFORD 1983
c **      ** see ADDMOD for parameter descriptions
c         number of parameters: 2
c            1      temperature
c            2      optical depth
C---
      REAL tkt, tau, at, emid, ek2, ek, ek3, pk, qk, eat
      REAL fac, bre, thermbr, acompls, edelt
      INTEGER i

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors

      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      tkt=param(1)
      tau=param(2)

      AT=2*TKT/511.*TAU*TAU

      DO i=1,ne
         EMID=(ear(i-1)+ear(i))/2.
         EK=EMID/TKT
         EK2=EK*EK
         EK3=EK2*EK
         PK=1.5*(EK3-EK2-EK)/AT
         QK=1.5*EK2/AT-0.75*EK3/AT
         EAT=0
         IF(AT.LT.170.)EAT=EXP(-AT)
         FAC=(1.+0.375*EK3*AT-PK*(EAT-1.+AT)+QK*(EAT*(AT+1.)-1.))
         BRE=THERMBR(EMID,TKT)
         aCOMPLS=FAC*BRE
         edelt=ear(i)-ear(i-1)
         photar(i)=edelt*acompls
      ENDDO

      RETURN
      END

C ******************************************************
      REAL FUNCTION THERMBR(E, T)
      IMPLICIT NONE
      REAL e, t, fbg, cfac, alf, gf
      IF (T.GT.20.) THEN
         GOTO 100
      ENDIF
      THERMBR = FBG(E, T, 1.) + .34*FBG(E, T, 2.)
      RETURN
 100  CFAC = 1./(.31968*EXP(-.3485*E)+.7840)
      ALF = .37*(30./T)**0.15
      GF = 0.90*(T/E)**ALF
      THERMBR = CFAC*GF/E*EXP(-E/T)
      RETURN
      END

C ******************************************************
      REAL FUNCTION FBG(E, Q, Z)
      IMPLICIT NONE
      REAL A(6, 7, 3), A1(42), A2(42), A3(42), GAM2(6), GAM3(6)
      REAL e, q, z, gam1, gam, u, u2, u1, ai, t, ak, u12
      REAL g, born, g1, g2, p
      INTEGER n, m, m1
      EQUIVALENCE (A1(1), A(1,1,1)), (A2(1), A(1,1,2)), 
     &            (A3(1), A(1,1,3))
      DATA GAM2/.7783, 1.2217, 2.6234, 4.3766, 20., 70./
      DATA GAM3/1., 1.7783, 3., 5.6234, 10., 30./
      DATA A1/1.001, 1.004, 1.017, 1.036, 1.056, 1.121, 1.001, 1.005,
     &     1.017, 1.046, 1.073, 1.115, .9991, 1.005, 1.03, 1.055, 1.102,
     &     1.176, .997, 1.005, 1.035, 1.069, 1.134, 1.186, .9962, 1.004,
     &     1.042, 1.1, 1.193, 1.306, .9874, .9962, 1.047, 1.156, 1.327,
     &     1.485, .9681, .9755, .8363, 1.208, 1.525, 1.965/
      DATA A2/.3029, .1616, .04757, .013, .0049, -.0032, .4905, .2155,
     &     .08357, .02041, .00739, .00029, .654, .2833, .08057, .03257,
     &     .00759, -.00151, 1.029, .391, .1266, .05149, .01274, .00324,
     &     .9569, .4891, .1764, .05914, .01407, -.00024, 1.236, .7579,
     &     .326, .1077, .028, .00548, 1.327, 1.017, 1.398, .205, .0605,
     &     .00187/
      DATA A3/ - 1.323, -.254, -.01571, -.001, -.000184, .00008, -4.762,
     &     -.3386, -.03571, -.001786, -.0003, .00001, -6.349, -.4206,
     &     -.02571, -.003429, -.000234, .00005, -13.231, -.59, -.04571,
     &     -.005714, -.000445, -.00004, -7.672, -.6852, -.0643,
     &     -.005857, -.00042, .00004, -7.143, -.9947, -.12, -.01007,
     &     -.000851, -.00004, 3.175, -1.116, -.8414, -.01821, -.001729,
     &     .00023/
      FBG = 0.0
      IF (Q.EQ.0.) RETURN
      GAM = Z*Z*.01358/Q
      GAM1 = GAM*1000.
      U = E/Q
      IF (U.GE.50.) RETURN
      G = 1.0
      IF (GAM1.GT.100.) GOTO 700
      U2 = U*U
      U1 = U/2.
      T = U1/3.75
      IF (U1.GT.2.) THEN
         GOTO 299
      ENDIF
      AI = 1. + 3.5156229*T*T + 3.089942*T**4 + 1.2067492*T**6 +
     &     .2659732*T**8 + .0360768*T**10 + .0045813*T**12
      U12 = U1/2.
      AK = -ALOG(U12)*AI - .57721566 + .4227842*U12**2 +
     &     .23069756*U12**4 + .0348859*U12**6 + .00262698*U12**8 +
     &     .0001075*U12**10 + .0000074*U12**12
      GOTO 297
 299  U12 = 2./U1
      AK = 1.25331414 - .07832358*U12 + .02189568*U12**2 -
     &     .01062446*U12**3 + .00587872*U12**4 - .0025154*U12**5 +
     &     .00053208*U12**6
      AK = AK/(EXP(U1)*SQRT(U1))
 297  BORN = .5513*EXP(U1)*AK
      IF (GAM1.LT.1. .OR. U.LT..003) THEN
         G = BORN
         GOTO 700
      ENDIF
      IF (U.LE..03) N = 1
      IF (U.LE..30 .AND. U.GT..03) N = 2
      IF (U.LE.1. .AND. U.GT..30) N = 3
      IF (U.LE.5. .AND. U.GT.1.0) N = 4
      IF (U.LE.15. .AND. U.GT.5.) N = 5
      IF (U.GT.15.) N = 6
      IF (GAM1.LE.1.7783) M = 1
      IF (GAM1.LE.3. .AND. GAM1.GT.1.7783) M = 2
      IF (GAM1.LE.5.6234 .AND. GAM1.GT.3.) M = 3
      IF (GAM1.LE.10. .AND. GAM1.GT.5.6234) M = 4
      IF (GAM1.LE.30. .AND. GAM1.GT.10.) M = 5
      IF (GAM1.LE.100. .AND. GAM1.GT.30.) M = 6
      M1 = M + 1
      G1 = (A(N,M,1)+A(N,M,2)*U+A(N,M,3)*U2)*BORN
      G2 = (A(N,M1,1)+A(N,M1,2)*U+A(N,M1,3)*U2)*BORN
      P = (GAM1-GAM3(M))/GAM2(M)
      G = (1.-P)*G1 + P*G2
 700  FBG = G*EXP(-E/Q)/E
      RETURN
      END

