      FUNCTION gaunt(E, Q, Z)

      REAL  gaunt, E, Q, Z

c Calculation of the Gaunt factor based on the program given in
c Kellogg, Baldwin, & Koch (ApJ 199, 299). Fixes supplied by Larry
c Molnar and Jack Hughes for inaccuracies in the polynomial fits
c to the Karzas & Latter (1961) numerical values. The routine by
c Kurucz for the low T case was supplied by Larry Molnar.

c arguments
c		e	r	i:energy (keV)
c		q	r	i: temperature (keV)
c		z	r	i: ion charge used for correction
c	returned value
c		gaunt	r	r: the gaunt factor

      REAL A(6, 7, 3), A1(42), A2(42), A3(42), GAM2(6), GAM3(6)
      REAL gam1, gam, u, u2, u1, ai, t, ak, u12
      REAL g, born, g1, g2, p, t2, t4, t8, u122, u12im
      INTEGER n, m, m1

      EQUIVALENCE (A1(1), A(1,1,1)), (A2(1), A(1,1,2)), 
     &            (A3(1), A(1,1,3))

      DATA GAM2/.7783, 1.2217, 2.6234, 4.3766, 20., 70./
      DATA GAM3/1., 1.7783, 3., 5.6234, 10., 30./

      DATA A1/1.001, 1.004, 1.017, 1.036, 1.056, 1.121, 1.001, 1.005,
     &     1.017, 1.046, 1.073, 1.115, .9991, 1.005, 1.03, 1.055, 1.102,
     &     1.176, .997, 1.005, 1.035, 1.069, 1.134, 1.186, .9962, 1.004,
     &     1.042, 1.1, 1.193, 1.306, .9874, .9962, 1.047, 1.156, 1.327,
     &     1.485, .9681, .9755, 1.020, 1.208, 1.525, 1.965/

      DATA A2/.3029, .1616, .04757, .013, .0049, -.0032, .4905, .2155,
     &     .08357, .02041, .00739, .00029, .654, .2833, .08057, .03257,
     &     .00759, -.00151, 1.029, .391, .1266, .05149, .01274, .00324,
     &     .9569, .4891, .1764, .05914, .01407, -.00024, 1.236, .7579,
     &     .326, .1077, .028, .00548, 1.327, 1.017, .6017, .205, .0605,
     &     .00187/

      DATA A3/ - 1.323, -.254, -.01571, -.001, -.000184, .00008, -4.762,
     &     -.3386, -.03571, -.001786, -.0003, .00001, -6.349, -.4206,
     &     -.02571, -.003429, -.000234, .00005, -13.231, -.59, -.04571,
     &     -.005714, -.000445, -.00004, -7.672, -.6852, -.0643,
     &     -.005857, -.00042, .00004, -7.143, -.9947, -.12, -.01007,
     &     -.000851, -.00004, -3.175, -1.116, -.2270, -.01821, -.001729,
     &     .00023/

c if temperature or energy is zero or if energy/temperature is
c too high then give up.

      IF ((Q.EQ.0.) .OR. (e.GT.50.*q) .OR. (E.EQ.0)) THEN
         gaunt = 0.
         RETURN
      ENDIF

c Convert the energy and temperature to the units used by Karzas & Latter

      U = E/Q
      GAM = Z*Z*.01358/Q
      GAM1 = MIN(GAM*1000., 100.)

c If we are in the low T regime then use Kurucz's algorithm

      IF (GAM.GT.0.1) THEN
         CALL kurucz(u, gam, g)
         gaunt = g
         RETURN
      ENDIF

c Calculate the Born approximation

      g = 1.
      U2 = U*U
      U1 = U*0.5
      T = U1/3.75
      u12 = u1*0.5
      IF (u12.LE.1.) THEN
         t2 = t*t
         t4 = t2*t2
         t8 = t4*t4
         AI = 1. + 3.5156229*T2 + 3.089942*T4 + 1.2067492*T2*t4 +
     &        .2659732*T8 + .0360768*T8*t2 + .0045813*T8*t4
         u122 = u12*u12
         AK = -ALOG(U12)*AI - .57721566 +
     &        u122*(.4227842+U122*(.23069756+
     &        U122*(.0348859+u122*(.00262698+
     &        U122*(.0001075+U122*.0000074)))))
      ELSE
         u12im = -1./u12
         AK = 1.25331414 +
     &        u12im*(.07832358+u12im*(.02189568+u12im*(.01062446+
     &        u12im*(.00587872+u12im*(.0025154+u12im*.00053208)))))
         AK = AK/(EXP(U1)*SQRT(U1))
      ENDIF
      BORN = .5513*EXP(U1)*AK

c If Born approximation is valid then go with it

      IF (gam1.LT.1.) THEN
         gaunt = born
         RETURN
      ENDIF

c All that is left to do is the polynomial approximation

      IF (u.LT.0.003) THEN
         u = .003
         u2 = u*u
      ENDIF
      IF (u.LE.0.03) THEN
         n = 1
      ELSEIF (u.LE.0.30) THEN
         n = 2
      ELSEIF (u.LE.1.0) THEN
         n = 3
      ELSEIF (u.LE.5.0) THEN
         n = 4
      ELSEIF (u.LE.15.0) THEN
         n = 5
      ELSE
         n = 6
      ENDIF
      IF (gam1.LE.1.773) THEN
         m = 1
      ELSEIF (gam1.LE.3.0) THEN
         m = 2
      ELSEIF (gam1.LE.5.6234) THEN
         m = 3
      ELSEIF (gam1.LE.10.) THEN
         m = 4
      ELSEIF (gam1.LE.30.) THEN
         m = 5
      ELSE
         m = 6
      ENDIF

      M1 = M + 1
      G1 = (A(N,M,1)+A(N,M,2)*U+A(N,M,3)*U2)*BORN
      G2 = (A(N,M1,1)+A(N,M1,2)*U+A(N,M1,3)*U2)*BORN
      P = (GAM1-GAM3(M))/GAM2(M)

      gaunt = (1.-P)*G1 + P*G2

      RETURN
      END
c
c  routine for low T case (from Larry Molnar)
c
      SUBROUTINE KURUCZ(UIN, GAM, GAUNT)
      INTEGER J, K
      REAL GAM, GAUNT, RJ, RK, T, U, UIN, YA(7, 12)
      DATA YA/5.40, 5.25, 5.00, 4.69, 4.48, 4.16, 3.85, 4.77, 4.63,
     &     4.40, 4.13, 3.87, 3.52, 3.27, 4.15, 4.02, 3.80, 3.57, 3.27,
     &     2.98, 2.70, 3.54, 3.41, 3.22, 2.97, 2.70, 2.45, 2.20, 2.94,
     &     2.81, 2.65, 2.44, 2.21, 2.01, 1.81, 2.41, 2.32, 2.19, 2.02,
     &     1.84, 1.67, 1.50, 1.95, 1.90, 1.80, 1.68, 1.52, 1.41, 1.30,
     &     1.55, 1.56, 1.51, 1.42, 1.33, 1.25, 1.17, 1.17, 1.30, 1.32,
     &     1.30, 1.20, 1.15, 1.11, 0.86, 1.00, 1.15, 1.18, 1.15, 1.11,
     &     1.08, 0.59, 0.76, 0.97, 1.09, 1.13, 1.10, 1.08, 0.38, 0.53,
     &     0.76, 0.96, 1.08, 1.09, 1.09/
      RJ = 2.*LOG10(GAM) + 3.
      J = INT(RJ)
      RJ = J
      RK = 2.*LOG10(UIN) + 9.
      K = INT(RK)
      IF (K.LT.1) K = 1
      RK = K
      IF (J.GE.7 .OR. J.LT.1 .OR. K.GE.12) THEN
         GAUNT = 0.
      ELSE
         T = (LOG10(GAM)-(RJ-3.)/2.)/0.5
         U = (LOG10(UIN)-(RK-9.)/2.)/0.5
         GAUNT = (1.-T)*(1.-U)*YA(J, K) + T*(1.-U)*YA(J+1, K)
     &            + T*U*YA(J+1, K+1) + (1.-T)*U*YA(J, K+1)
      ENDIF
      RETURN
      END
