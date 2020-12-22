C*******************************************************************
C HLYMAN CALCULATES THE VOIGT ABSORPTION PROFILES FOR THE H I AND  *
C HE II LYMAN SERIES.  ABSORPTION PROFILES (CM) AS A FUNCTION OF   *
C WAVELENGTH ARE RETURNED FOR A DOPPLER VELOCITY B(KM/S) BETWEEN   *
C NELEM = 1 (HYDROGEN ATOM) AND NELEM = 2 (HELIUM).    3/10/96     *
C SUBROUTINES: VOIGTUTL ---> FUNCTIONS: DAWSON                        *
C*******************************************************************
      subroutine xslyman(ear,ne,param,ifl,photar,photer)
c               see MULMOD for parameter descriptions
c       number of model parameters: 5
c               1       N - column density (Units of 10**22)
c               2       b - doppler velocity (km/s)
c               3       Z - redshift
c               4       ZA - atomic number
      IMPLICIT NONE
      INTEGER ne,n,jmax,j,ifl
      REAL param(4),ear(0:ne)
      REAL photar(ne),photer(ne)
      REAL heii,atokev,z,b,g1,g2,wmin,wmax
      REAL wn,tau,pi,sqrtpi,c,s,ang,kmcm
      REAL asum,beta,alpha,voigtutl,za
      REAL WLAM,ALAM,U
      REAL W(39),A(39),F(39),V0(39),W1(39)

C***  WAVELENGTH & F VALUES FROM WISE, SMITH & GLENNON (1966) ***** 

      DATA W/1215.6701, 1025.7223,  972.5368,  949.7431,  937.8035,
     +        930.7483,  926.2257,  923.1504,  920.9631,  919.3514,
     +        918.1294,  917.1806,  916.4290,  915.824 ,  915.329 ,
     +        914.919 ,  914.576 ,  914.286 ,  914.039 ,  913.826 ,
     +        913.641 ,  913.48  ,  913.339 ,  913.215 ,  913.104 ,
     +        913.006 ,  912.918 ,  912.839 ,  912.768 ,  912.703 ,
     +        912.645 ,  912.592 ,  912.543 ,  912.499 ,  912.458,
     +         912.42 ,  912.385 ,  912.353 ,  912.324/      
     
      DATA F/0.4164   ,  0.07912 ,  0.02900 ,  0.01394 ,  0.007804,
     +       0.004816 ,  0.003183,  0.002216,  0.001605,  0.001200,
     +       0.0009213, 0.0007226, 0.0005774, 0.0004686, 0.0003856,
     +       0.0003211, 0.0002702, 0.0002296, 0.0001967, 0.0001698,
     +       0.0001476, 0.0001291, 0.0001136, 0.0001005, 8.928e-5 ,
     +       7.97e-5  , 7.144e-5 , 6.429e-5 , 5.806e-5 , 5.261e-5 ,
     +       4.782e-5 , 4.360e-5 , 3.986e-5 , 3.653e-5 , 3.357e-5 ,
     +       3.092e-5 , 2.854e-5 , 2.640e-5 , 2.446e-5/
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
      n = ifl

c this model does not calculate errors
      DO n = 1, ne
         photer(n) = 0.0
      ENDDO

C********************  CONSTANTS  ******************************
      atokev = 12.39854
      PI   = 4.0*ATAN(1.0)
      SQRTPI = SQRT(PI)
      C    = 2.99792458E10
      S    = 8.85282254E-13
      ANG  = 1.0E-8 
      KMCM = 1.0E5
      G1   = 2.0
      G2   = 6.0
      JMAX = 39



      heii=param(1)
      b=param(2)
      if(b.le.0.1)b=0.1
      z=param(3)
      za=param(4)

      WMIN = 912/za**2
      WMAX = 1220/za**2

C****************** DOPPLER CONSTANTS   *************************

      B = B*KMCM
      BETA  = B/C
      ALPHA = (1.0/SQRTPI)*S/BETA

C*******   CALCULATE INDIVIDUAL VOIGT PARAMETERS  ***************

      DO J = 1,JMAX
         W1(J) = W(J)/za**2                          
         V0(J) = ALPHA*F(J)*W1(J)*ANG                 
         A(J) = 2.0*SQRTPI*ALPHA*(G1/G2)*F(J)/(W1(J)*ANG)       
      ENDDO

C*******  CALCULATE VOIGT PROFILES  *****************************
        

      DO n=1,ne
         wn = atokev/ear(n)
         WLAM = wn-z*wn
         if(WLAM.lt.WMIN.or.WLAM.gt.WMAX)then
            ALAM=0.0
         else
            ASUM = 0.0
            DO J = 1, JMAX 
               U = ABS((1.0 - W1(J)/WLAM)/BETA)
               ASUM = V0(J)*VOIGTUTL(A(J),U) + ASUM
            ENDDO
            ALAM = ASUM
         endif

         TAU = 1.0e22*heii*ALAM
         IF (TAU. GT. 35.0)  TAU = 35.0
C --(this traps large negative exponents which could cause underflow)
         photar(n) = EXP(-TAU)
         
      enddo
c
      RETURN
      END













