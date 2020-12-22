C*******************************************************************
C HEI CALCULATES THE VOIGT ABSORPTION PROFILES FOR THE HE I SERIES *
C SERIES WLIM = 504A.  ABSORPTION PROFILES (CM) AS A FUNCTION OF   *
C WAVELENGTH ARE RETURNED FOR A DOPPLER VELOCITY B(KM/S) BETWEEN   *
C500 & 520 A.           3/8/96                                     *
C SUBROUTINES: VOIGTUTL ---> FUNCTIONS: DAWSON                        *
C*******************************************************************
      subroutine xsphei(ear,ne,param,ifl,photar,photer)
c               see MULMOD for parameter descriptions
c       number of model parameters: 5
c               1       HeI - column density (Units of 10**22)
c               2       b - doppler velocity (km/s)
c               3       Z - redshift
      IMPLICIT NONE
      INTEGER ne,n,jmax,j,ifl
      REAL param(3),ear(0:ne)
      REAL photar(ne),photer(ne)
      REAL hei,atokev,z,b,g1,g2,wmin,wmax
      REAL wn,tau,pi,sqrtpi,c,s,ang,kmcm
      REAL asum,beta,alpha
      REAL WLAM,ALAM,U
      REAL W(20),A(20),F(20),V0(20)

      REAL voigtutl
      EXTERNAL voigtutl

C****  WAVELENGTH & F VALUES FROM THEODOSIOU (1987) ************** 

      DATA W/ 584.3342, 537.0297, 522.2130, 515.6166, 512.0983,
     +        509.9998, 508.6431, 507.7179, 507.0578, 506.5704, 
     +        506.2002, 505.9124, 505.6840, 505.5002, 505.3496, 
     +        505.2250, 505.1205, 505.0322, 504.9568, 504.8918/      
     
      DATA F/0.2764300,  0.0733360,  0.0298040,  0.0149950,  0.0086030,       
     +       0.0053897,  0.0035999,  0.0025233,  0.0018356,  0.0013779,
     +       0.0010607,  0.00083368, 0.00066728, 0.00054237, 0.00044684,
     +       0.00037251, 0.00031372, 0.00026894, 0.00022886, 0.00019803/
C******************* DEFINITIONS  ********************************
       
C    C = Speed of Light (cm/s)
C    S = Line Strength (pi*e^2/mc^2)  cgs (cm)
C    B = Doppler Parameter b (km/s)
C BETA = B/C
C  ANG = Angstrom Conversion Factor (cm/A)
C    F = Absorption Oscillator Strength f
C   G1 = Upper g Value
C   G2 = Lower g Value    
C    A = Voigt Parameter (natural width to doppler width)
C   V0 = Voigt Coefficient (cm^2)
C    U = Off-Resonance Distance in Doppler Units
C WLAM = Output Wavelength Vector (A)
C ALAM = Output Absorption Coeff. Vector (cm)

c suppress a warning message from the compiler
      j = ifl

c this model does not calculate errors
      DO j = 1, ne
         photer(j) = 0.0
      ENDDO


C********************  CONSTANTS  ******************************

      PI   = 4.0*ATAN(1.0)
      SQRTPI = SQRT(PI)
      C    = 2.99792458E10
      S    = 8.85282254E-13
      ANG  = 1.0E-8 
      KMCM = 1.0E5
      G1   = 1.0
      G2   = 3.0
      WMIN = 500.0
      WMAX = 595.0
      JMAX = 20
      atokev = 12.39854


      hei=param(1)
      b=param(2)
      if(b.le.0.1)b=0.1
      z=param(3)


c **********   compute absorption profiles *****************************




C****************** DOPPLER CONSTANTS   *************************

      B = B*KMCM
      BETA  = B/C
      ALPHA = (1.0/SQRTPI)*S/BETA

C*******   CALCULATE INDIVIDUAL VOIGT PARAMETERS  ***************

      DO J = 1,JMAX                          
         V0(J) = ALPHA*F(J)*W(J)*ANG                 
         A(J) = 2.0*SQRTPI*ALPHA*(G1/G2)*F(J)/(W(J)*ANG)       
      ENDDO

C*******  CALCULATE VOIGT PROFILES  *****************************
        

      do n=1,ne
         wn = atokev/ear(n)
         WLAM = wn-z*wn
         if(WLAM.lt.WMIN.or.WLAM.gt.WMAX)then
            ALAM=0.0
         ELSE
            ASUM = 0.0
            DO J = 1, JMAX 
               U = ABS((1.0 - W(J)/WLAM)/BETA)
               ASUM = V0(J)*VOIGTUTL(A(J),U) + ASUM
            ENDDO
            ALAM = ASUM
         endif

         TAU = 1.0e22*hei*ALAM
         IF (TAU. GT. 35.0)  TAU = 35.0
C --(this traps large negative exponents which could cause underflow)
         photar(n) = EXP(-TAU)

      enddo
c       

      RETURN
      END


      REAL FUNCTION VOIGTUTL(A,U)
C****************************************************************
C* VOIGT FUNCTION - Using the expansion of Harris (1949) for    *
C* small 'A'.  CALLS DAWSON'S FUNCTION                          *
C*    SINGLE PRECISION  3/8/96                                  *
C****************************************************************
      IMPLICIT NONE
      REAL H(6),c1,u,u2,u4,u6,a,h0
      REAL D
      INTEGER I

      DOUBLE PRECISION DAWSON
      EXTERNAL DAWSON

      C1 = -1.0/SQRT(ATAN(1.0))
      U2 = U*U
      IF(U2. GT. 72.) U2 = 72.
      U4 = U2*U2
      U6 = U4*U2
      D  = REAL(DAWSON(U))
      H0   = EXP(-U2)
      H(2) = (1.0 - 2.0*U2)*H0
      H(4) = (1./6.0)*(3.0 - 12.0*U2 + 4.0*U4)*H0
      H(6) = (1./90.)*(15.0 - 90.0*U2 + 60.0*U4 - 8.0*U6)*H0
      H(1) = C1*(1.0 - 2.0*U*D)
      H(3) = C1*(2./3.)* (1.0 - U2 - U*(3.0 - 2.0*U2)*D)
      H(5) = C1*(1./30.)*(8.0 - 18.0*U2 + 4.0*U4 - U*
     +     (30.0 - 40.0*U2 + 8.0*U4)*D)

C     SUM TERMS

      VOIGTUTL  =  0.0
      DO I = 6,1,-1
         VOIGTUTL = H(I) + VOIGTUTL*A
      ENDDO
      VOIGTUTL = H0   + VOIGTUTL*A
      RETURN
      END

      DOUBLE PRECISION FUNCTION DAWSON(U) 
C************************************************************
C* DAWSON'S FUNCTION  USES CHEBYSHEV COEFFICIENTS OF HUMMER *
C* (1964), EXPANSION TERIMINATED AT N=22, ERROR < 1E-7      *
C*        -  DOUBLE PRECISION 3/14/96                       * 
C************************************************************
      IMPLICIT NONE
      REAL U
      INTEGER I,N
      DOUBLE PRECISION C(22), T(43)
      DOUBLE PRECISION X,D,TERM
      DATA C / 0.1999999999972224D0, -0.1840000000029998D0,
     +     0.1558399999965025D0, -0.1216640000043988D0,
     +     0.0877081599940391D0, -0.0585141248086907D0,
     +     0.0362157301623914D0, -0.0208497654398036D0,
     +     0.0111960116346270D0, -0.0056231896167109D0,
     +     0.0026487634172265D0, -0.0011732670757704D0,
     +     0.0004899519978088D0, -0.0001933630801528D0,
     +     0.0000722877446788D0, -0.0000256555124979D0,
     +     0.0000086620736841D0, -0.0000027876379719D0,
     +     0.0000008566873627D0, -0.0000002518433784D0, 
     +     0.0000000709360221D0, -0.0000000191732257D0/
C     +          0.0000000049801256D0, -0.0000000012447734D0,
C     +          0.0000000002997777D0, -0.0000000000696450D0,
C     +          0.0000000000156262D0, -0.0000000000033897D0,
C     +          0.0000000000007116D0, -0.0000000000001447D0,
C     +          0.0000000000000285D0, -0.0000000000000055D0,
C     +          0.0000000000000010D0, -0.0000000000000002D0 /

      X =DBLE(U)
      IF(X.LT.5.0D0) THEN      
      
         N = 20

C *************** COMPUTE CHEBYSHEV POLYNOMIALS ****************** 
      
         X = 0.2D0*X
         T(1) = X
         T(2) = 2.0D0*X*X - 1.0D0
      
         DO I = 2,2*N+1
            T(I+1)   = 2.0D0*X*T(I) - T(I-1)
         ENDDO

         C(1) = 0.1999999999972224D0
         D = C(1)*T(1)
         
         DO I=1,N
            D = C(I+1)*T(2*I+1) + D
         ENDDO

         DAWSON = D
      
      ELSE
      
         TERM = 0.5D0/X
         D = 0.0

         DO I = 1,11
            D = TERM + D
            TERM = (2*I-1)*TERM/(2.D0*X*X)
         ENDDO
      
         DAWSON = D
      
      ENDIF

      RETURN
      END












