      SUBROUTINE DISKPBB(EAR,NE,PARAM,IFL,PHOTAR,PHOTER)
      implicit none
      INTEGER NE, IFL
      REAL    EAR(0:NE), PARAM(2), PHOTAR(NE), PHOTER(NE)
      
c     Multicolour disk blackbody model used in ISAS, Japan.
c     See Mitsuda et al. PASJ, 36, 741 (1984)
c     & Makishima et al. ApJ 308, 635 1986)
c     Ken Ebisawa 1992/12/22
      
C     Modified to use double precision and to make numerical
C     integration faster.
C     Ken Ebisawa 1993/07/29
      
C     The exponent p (T(r) ~ r^-p) is made a free parameter.
C     Ken Ebisawa 2000-3-1

C     The bug was fixed that the original normalization was wrong by
C     factor (0.75/p).  The true normalization is (0.75/p) times the original
C     normalization.  
C     Also, the numerical integration in temperature (disk radius)
C     is made more precise to lower energies (to outer
C     radis) so that the model gives more precise fluxes at E<<0.1 keV.
C     Ken Ebisawa and Aya Kubota 2006-10-26

      INTEGER I, J
      double precision XN, XH
      
      double precision TIN, E, photon
      real p

C     These coefficients are taken from the spectral fitting program SPFD
C     in ISAS.
      DOUBLE PRECISION  GAUSS(5,2)
      DATA GAUSS /
     $     0.236926885,  0.478628670, 0.568888888, 0.478628670,
     $     0.236926885, -0.906179846, -0.538469310, 0.0,  0.538469310,
     $     0.906179846 /   

c suppress a warning message from the compiler
      i = ifl

c this model has no errors
      DO I = 1, NE
         PHOTER(I) = 0.0
      ENDDO
      
      TIN = DBLE(PARAM(1))
      p   = PARAM(2)

      DO 100 I = 1, NE
         XN = (dble(EAR(I))-dble(EAR(I-1)))/2.0D0
         PHOTAR(I) = 0.0
         XH = XN + dble(EAR(I-1))
         DO 200 J = 1, 5
            E = XN * GAUSS(J,2) + XH
            CALL DKBFLX(TIN, p, E, PHOTON)
            PHOTAR(I) = PHOTAR(I) + REAL(GAUSS(J,1) * PHOTON)
 200     CONTINUE
         PHOTAR(I) = PHOTAR(I) * SNGL(XN)
 100  CONTINUE
      END
       
      SUBROUTINE DKBFLX(TIN,p,E,PHOTON)
      implicit none
      integer i, j, nnn
      double precision TIN,E,PHOTON, DKFUNC
      double precision X, XN, XH
      real p
      
C     DISK BLACKBODY SPECTRUM USED IN ISAS. KEN EBISAWA 1991/01/14
C     SEE MITSUDA ET AL. 1984 PASJ 36, 741.
C     & MAKISHIMA ET AL. APJ 308, 635 1986)
C     NORMALIZATION={RIN(KM)/(D/10KPC)}^2*COS(THETA)
C     TIN=     INNER TEMPERATURE OF THE DISK
C     E  =     ENERGY
C     PHOTON = PHOTONS/S/KEV/CM2
      
C     modified to use double precision variables and to make numerical
C     integration faster 
C     Ken Ebisawa 93/07/29
      
C     These coefficients are taken from the spectral fitting program SPFD
C     in ISAS.

C     New paramer p is introduced for the exponent of the
C     temperature radial dependence (T(r) ~ r^-p).
C     Ken Ebisawa 2000-3-1

      DOUBLE PRECISION  GAUSS(5,2)
      DATA GAUSS /
     $     0.236926885,  0.478628670, 0.568888888, 0.478628670,
     $     0.236926885, -0.906179846, -0.538469310, 0.0,  0.538469310,
     $     0.906179846 /   
      
C     NNN determines the accuracy of the numerical integration and the time
C     of the calculation.
      IF(E.GT.TIN) THEN
         NNN=5
      ELSEIF(E.GT.0.2*TIN)THEN
         NNN = 10
      ELSEIF(E.GT.0.1*TIN)THEN
         NNN = 100
      ELSEIF(E.GT.0.01*TIN)THEN
         NNN = 500
      ELSEIF(E.GT.0.001*TIN)THEN
         NNN = 1000
      ELSE
         NNN = 10000
      ENDIF
      
      XN = 1.0D0/NNN/2.0D0      
      
      PHOTON=0.0D0
      
      DO 100 I = 1, NNN
         XH = XN * (2*I - 1)
         DO 200 J = 1, 5
            X = XN * GAUSS(J,2) + XH
            PHOTON = PHOTON + GAUSS(J,1)*DKFUNC(E,p,TIN, X)*XN
 200     CONTINUE
 100  CONTINUE
      
c      PHOTON=2.78D-3*E*E*PHOTON
C     We have to take into account the "p" value here.
C     2006-10-26 by Ken Ebisawa
      PHOTON=2.78D-3*E*E*PHOTON*(0.75D0/p)
      
      END
      
      FUNCTION DKFUNC(E,p,TIN,X)
      implicit none
      double precision DKFUNC,E, TIN,X
      real p

      IF(E/TIN/X.GE.170.0D0) THEN
c         DKFUNC=X**(-11.0D0/3.0D0)*EXP(-E/TIN/X)
         DKFUNC=X**(-2.0D0/DBLE(p)-1.0D0)*EXP(-E/TIN/X)
      ELSE
c         DKFUNC=X**(-11.0D0/3.0D0)/( EXP(E/TIN/X)-1.0D0 )
         DKFUNC=X**(-2.0D0/DBLE(p)-1.0D0)/( EXP(E/TIN/X)-1.0D0 )
      ENDIF
      END



