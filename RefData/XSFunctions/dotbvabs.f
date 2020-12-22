
      subroutine dotbvabs(ear, ne, param, ifl, photar, photer, siggas,
     &                    siggrains, sigmol, sigavg)
      INTEGER ne, ifl
      REAL ear(0:ne), param(42), photar(ne), photer(ne)
      REAL siggas(ne),siggrains(ne),sigmol(ne),sigavg(ne)

c
c     Compute X-ray absorptivity of cold gas using the formalism given
c     by Wilms, Allen, and McCray, ApJ, 2000, in press
c
c     This is the main subroutine -- beware!
c
c     The parameters are
c
c     param( 1)   hydrogen column (1E22/cm^2)
c
c     param(2..30): relative abundance of the element wrt the currently
c          set abundance table
c     param( 2)   helium      (2)
c     param( 3)   carbon      (6)
c     param( 4)   nitrogen    (7)
c     param( 5)   oxygen      (8)
c     param( 6)   neon        (10)
c     param( 7)   sodium      (11)
c     param( 8)   magnesium   (12)
c     param( 9)   aluminum    (13)
c     param(10)   silicon     (14)
c     param(11)   sulphur     (16)
c     param(12)   chlorine    (17)
c     param(13)   argon       (18)
c     param(14)   calcium     (20)
c     param(15)   chromium    (24)
c     param(16)   iron        (26)
c     param(17)   cobalt      (27)
c     param(18)   nickel      (28)
c     param(19): fraction of hydrogen existing in molecular form
c     param(20): density of dust in g/cm^3
c        if param(20) <= 0 then no dust is included in the modeling
c        (i.e., the depletion factor is also not taken into account)
c     param(21): minimum thickness of dust (mum)
c     param(22): maximum thickness of dust (mum)
c     param(23): PL index for dust size distribution (implicit minus-sign!)
c        if param(21) = param(22): use single size grain
c
c     param(24..42): depletion factor (ratio of GAS to total ISM abundance)
c     param(24)   hydrogen    (1)
c     param(25)   helium      (2)
c     param(26)   carbon      (6)
c     param(27)   nitrogen    (7)
c     param(28)   oxygen      (8)
c     param(29)   neon        (10)
c     param(30)   sodium      (11)
c     param(31)   magnesium   (12)
c     param(32)   aluminium   (13)
c     param(33)   silicon     (14)
c     param(34)   sulphur     (16)
c     param(35)   chlorine    (17)
c     param(36)   argon       (18)
c     param(37)   calcium     (20)
c     param(38)   chromium    (24)
c     param(39)   iron        (26)
c     param(40)   cobalt      (27)
c     param(41)   nickel      (28)
c
c     param(42):  Redshift

c     ... Ratio of H2 to H I
      REAL h2

c     ... Dust parameters: min. and max. thickness (cm!!), PL index
c     ...   density
      REAL amin,amax,pl,rho
      REAL massgrain,a,ngrains,nbar,fac,k,sigmin,sigmax,sigdummy

c     ... External references to cross section function
      REAL photo, gphoto, fgabnd
      CHARACTER(4) fgxsct
      EXTERNAL photo, gphoto, fgabnd, fgxsct

c     ... External references to the incomplete gamma function
c     ... from Slatec library
      DOUBLE PRECISION dgami
      EXTERNAL dgami

c     ... The H, He, and H2 cross sections
      REAL hphoto, hephoto, h2photo
      EXTERNAL hphoto, hephoto, h2photo

c     ... Other variables
      REAL zfac,ab,elow,ehi,sig,meanmol,normal

c     ... Atomic mass unit in grams
      REAL amu
      PARAMETER(amu=1.6605655E-24)

c     ... a number :-)
      REAL PI
      PARAMETER(PI=3.1415926535897932384626433)

      INTEGER NELE
      PARAMETER (NELE=18)

      INTEGER atomic(NELE)
      INTEGER ie, status, i,z

      character(2) celts(NELE)
      REAL weight(NELE)
      REAL abun(NELE),beta(NELE)
      LOGICAL startup

      DATA atomic /1,    2,    6,    7,    8,    10,   11,   12,
     &             13,   14,   16,   17,   18,   20,   24,   26,
     &             27,   28/

c
c     Missing compared to the paper
c     ...P
c     ...Ti
c     ...Mn

c
      DATA celts /'H ', 'He', 'C ', 'N ', 'O ', 'Ne', 'Na', 'Mg',
     &            'Al', 'Si', 'S ', 'Cl', 'Ar', 'Ca', 'Cr', 'Fe',
     &            'Co', 'Ni'/

      DATA weight /1.,  4.,  12.,  14.,  16.,   20.,  23.,  24.,
     &            27.,  28.,  32.,  35.,  40.,  40.,  52.,  56.,
     &            59.,  59./
      
      DATA startup/.true./ 

      IF (startup) THEN
c        
c        Copyright (only once, and only to stdout)
c
         call xwrite('tbvabs Version 1.0',1) 
         call xwrite('Cosmic absorption with grains and H2',1)
         call xwrite('Wilms, Allen, & McCray 2000',1) 
         call xwrite('Questions: Adrienne Allen or Joern Wilms',1)
         call xwrite('allenu@super.colorado.edu',1)
         call xwrite('wilms@astro.uni-tuebingen.de',1)
         startup=.false.
      ENDIF 

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      status=0

c     ... Redshift
      zfac=1.+param(42)

c
c     ... Put parameters in more meaningful variables
c
      h2=param(19)

c     ... grains
      rho=param(20)
      amin=param(21)*1E-4
      amax=param(22)*1E-4
      pl=param(23)
c     ... Set the element abundances
      abun(1)=1.
      DO i = 2, NELE
         abun(i)=param(i)*fgabnd(celts(i))
      ENDDO

c     ... Fraction of atoms in the grain phase
      IF (rho.gt.0.) THEN 
         DO i= 1, NELE
            beta(i)=1.0-param(23+i)
         ENDDO
      ELSE 
         DO i= 1, NELE
            beta(i)=0.0
         ENDDO
      ENDIF      

c     ... Reset the arrays
      DO ie=1,ne
         siggas(ie)=0.
         siggrains(ie)=0.
         sigmol(ie)=0.
         sigavg(ie)=0.
      ENDDO 

c
c     ... computation for the gas phase
c
      sig = 0.0
      DO i=1,NELE 
         z=atomic(i)
         ab=abun(i)*(1.-beta(i))

c        ... take care of molecular H
         IF (z.eq.1) THEN
            ab=ab*(1.-h2)
         ENDIF

c        ... sum total cross section for the gas phase
         IF (ab.ne.0.) THEN
            elow = ear(0) * zfac
            DO ie=1,ne
               ehi=ear(ie)*zfac
c              ... compute below 100keV only
               IF (elow.le.100.) THEN
                  IF (z.eq.1) sig=hphoto(elow,ehi)*1E18
                  IF (z.eq.2) sig=hephoto(elow,ehi)*1E18
                  IF (z.gt.2) sig=gphoto(elow,ehi,z,status)*1E18
                  siggas(ie)=siggas(ie)+ab*sig
                  elow=ehi
               ELSE
                  siggas(ie)=0.
               ENDIF
            ENDDO
         ENDIF
      ENDDO

c
c     ... computation for the grain phase
c
      normal=0.
      IF (rho.gt.0.) THEN
c
c        ... Mean Weights
c
         massgrain=0.
         DO i=1,NELE
            massgrain=massgrain+abun(i)*beta(i)*weight(i)
            normal=normal+abun(i)*beta(i)
         ENDDO
      ENDIF 

c     ... Sanity check to see whether beta tells us that grain exists
      IF (normal.gt.0.) THEN 
c        ... Mass in grain (in g per H atom)
         massgrain=massgrain*amu
       
c        ... Mean weight per atom in g
         meanmol=massgrain/normal
c
c        ... Compute average cross-section per atom, sigavg
c
         DO i=1,NELE
            z=atomic(i)
            ab=abun(i)*beta(i)
            IF (ab.ne.0.) THEN
               elow = ear(0) * zfac
               DO ie=1,ne
                  ehi=ear(ie)*zfac
c                 ... compute below 100keV only
                  IF (elow.le.100.) THEN
                     IF (z.eq.1) sig=hphoto(elow,ehi)*1E18
                     IF (z.eq.2) sig=hephoto(elow,ehi)*1E18
                     IF (z.gt.2) sig=gphoto(elow,ehi,z,status)*1E18
                     sigavg(ie)=sigavg(ie)+ab*sig
                     elow=ehi
                  ELSE
                     sigavg(ie)=0.
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
c        ... normalize
         DO ie=1,ne
            sigavg(ie)=sigavg(ie)*1E-18/normal
         ENDDO

c
c        ... Computation for the grains
c
         IF (amin.eq.amax) THEN
c
c           ... Case of grains of single size
c
            a=amin

c           ... Number of grains per H atom
            ngrains=massgrain/(rho*4./3.*PI*a**3.)

c           ... Column of grain (atoms/cm^2)
            nbar=(4./3.)*rho*a/meanmol

c           ... grain cross-section per H atom 
            fac=ngrains*PI*a*a
            DO ie=1,ne
               siggrains(ie)=fac*(1.-exp(-sigavg(ie)*nbar))*1E18
            ENDDO
         ELSE
c
c           ... Case of grains with a power-law distribution
c

c           ... Normalization of PL Distribution
            IF (pl.ne.1.) THEN
               k=(1.-pl) / ( amax**(1.-pl) - amin**(1.-pl) )
            ELSE 
               k=1./alog(amax/amin)
            ENDIF

c           ... Number of grains per H atom
            ngrains=massgrain/(rho*4./3.*PI*k)
            IF (pl.ne.4.) THEN 
               ngrains=ngrains*(4.-pl)/( amax**(4.-pl) - amin**(4.-pl) )
            ELSE 
               ngrains=ngrains/alog(amax/amin)
            ENDIF

c
c           ... Formulae from the appendix of the Wilms et al. paper
c           ... incl. conversion to Mbarns
c
            DO ie=1,ne
               sigdummy=(4./3.)*sigavg(ie)*rho/meanmol
               sigmax=sigdummy*amax
               sigmin=sigdummy*amin
               siggrains(ie)=SNGL(ngrains* PI *k/(pl-3.) 
     &              * (amin**(3.-pl)*(1.-exp(-sigmin)) - 
     &                 amax**(3.-pl)*(1.-exp(-sigmax)) + 
     &                 sigdummy**(pl-3.) * ( 
     &                   dgami(dble(4.-pl),dble(sigmax)) - 
     &                   dgami(dble(4.-pl),dble(sigmin))
     &                                     ) 
     &                 )*1E18)
            ENDDO
         ENDIF
      ENDIF

c
c     ... Molecular phase
c     ... (everything for H2 only, but other molecules can be added)
c
      IF (h2.ne.0.) THEN
c        ... not that H will ever be depleted (ice??)
         ab=abun(1)*(1.-beta(1))
c        ... abundance of H2 by number (two H go into one H2)
         ab=ab*h2*0.5
         elow = ear(0) * zfac
         DO ie=1,ne
            ehi=ear(ie)*zfac
            sigmol(ie)=h2photo(elow,ehi)*ab*1E18
            elow=ehi
         ENDDO
      ENDIF
c        
c     ... add everything up and convert to cm^2
c     
      DO ie=1,ne
         photar(ie)=exp(-param(1)*1E22*
     1        (siggas(ie)+siggrains(ie)+sigmol(ie))*1E-18)
      ENDDO
      return
      END


      FUNCTION h2photo(keV1,keV2)
c
c     Wrapper to h2abs to return the photoelectric absorption cross
c     section in cm**2 of the hydrogen molecule
c     
c     calls h2abs below
c
      implicit none
      REAL h2photo
      REAL keV1,keV2,h2abs
      EXTERNAL h2abs
      h2photo=(h2abs(keV1*1000.)+h2abs(keV2*1000.))/2. * 1E-18
      return
      END

      FUNCTION h2abs(e)
c   
c     Cross section for the H2 molecule according to
c     Yan et al., 1998, ApJ 496, 1044-1050
c
c     e: energy in eV
c     sig: cross-section in Mbarns
c
c     Version 1.0, 1999/05/05 Joern Wilms,
c     wilms@astro.uni-tuebingen.de
c
c     Version 1.1, 1999/07/06 Joern Wilms
c     Updated to cross section of Wilms, Allen, McCray,1999,ApJ, submitted
c     for energies between 30 and 85eV
c
      implicit none
      real h2abs
      real x,e,xhlp,sig
      integer i

      real coeff(6)
      DATA coeff/  0.664E6, -11.768E6, 78.118E6, -231.339E6,
     &           368.053E6,-189.953E6/

      IF (e.LE.15.4) THEN
         h2abs=0.
         return
      ENDIF

      x=e/15.4E0
      xhlp=1./x
      sig = 0.0
      
      IF (e.GT.15.4.AND.e.LE.18.) THEN 
         sig=1E8*(((25.4*x-87.227)*x+99.723)*x-37.895)
      ENDIF 

c     ... original in yan: from 18 to 85eV
      IF (E.GT.18..AND.E.LE.30.) THEN
         sig=(((-0.692*xhlp+1.977)*xhlp-0.673)*xhlp+0.071) 
     &        *x**(-0.252)*2E7
      ENDIF 

c     ... Wilms, Allen, McCray fit
      IF (e.GT.30..AND.E.LE.85.) THEN 
         DO i=6,1,-1 
            sig=sig*xhlp+coeff(i)
         ENDDO
      ENDIF 

c     jw: ... might be optimized with polys....
      IF (e.GT.85.) THEN 
         sig=45.57*(1.-2.003*x**(-0.5)-4.806*x**(-1.)+ 
     &              50.577*x**(-1.5)-171.044*x**(-2.)+ 
     &              231.608*x**(-2.5)-81.885*x**(-3.))
     &        /(e/1000.)**(3.5)
      ENDIF
   
      h2abs=sig*1E-6
      return 
      END 
      FUNCTION hephoto(keV1,keV2)
c
c     Wrapper to heabs to return the photoelectric absorption cross section
c     in cm**2 of the helium atom
c
c     calls habs below
c
c
      REAL hephoto
      REAL keV1,keV2
      REAL heabs
      EXTERNAL heabs
      hephoto=(heabs(keV1*1000.)+heabs(keV2*1000.))/2. * 1E-18
      return
      END

      FUNCTION heabs(ener)
c
c     Photoionization cross-section for helium after
c     Yan et al., 1998, ApJ 496, 1044-1050
c
c     ener: energy in eV
c     sig: cross-section in Mbarns
c
c     Version 1.0
c     J. Wilms, wilms@astro.uni-tuebingen.de, 1999/04/22
c
      REAL heabs
      REAL ener
      DOUBLE PRECISION sig,sq
      INTEGER i

      external jwfano
      REAL q(4),nu(4),gamma(4),jwfano

c     Fitting coefficients (Yan et al., Tab. 4)
      DOUBLE PRECISION a(7)
      DATA a/1.D0,-4.7416D0,14.8200D0,-30.8678D0,37.3584D0,
     1     -23.4585D0,5.9133D0/

C     ... parameters Q for resonances (Fernley et al. 1987) 
      DATA q /2.81E0, 2.51E0, 2.45E0, 2.44E0/

C     ... parameters NU for resonances (Oza 1986) 
      DATA nu /1.610E0, 2.795E0, 3.817E0, 4.824E0/

C     ... parameters GAMMA for resonances (Oza 1986) 
      DATA gamma /2.64061E-3, 6.20116E-4, 2.56061E-4,
     +             1.320159E-4/

c     Threshold
      IF (ener.lt.24.58) THEN
         heabs=0.
         return
      ENDIF

c     1./x**(1/2) in Yan notation
      sq=sqrt(24.58D0/ener)

c     equation (14)
      sig=a(7)
      DO i=6,1,-1 
         sig=sig*sq+a(i)
      ENDDO
      sig=sig*733D-6/(ener/1000.)**3.5D0

C     Add in autoionization resonances
      DO i=1,4
         sig =sig*jwfano(q(i),nu(i),gamma(i),ener)
      ENDDO

      heabs=real(sig)
      return
      END 

      FUNCTION jwfano(Q,NU,GAMMA,ENER)
C
C
C    Source :
C             Fernley, J. A., Taylor, K. T., and Seaton, M. J., (1987),
C                J. Phys. B., 20, 6457.
C             Oza, D. H., (1986), Phys. Rev. A, 33,  824.
C
C    Description :
C     returns a Fano line profile for use by HELIUM. The form of
C     this line shape is taken from Fernley, Taylor and Seaton (above).
C
C    Deficiencies :
C     works in the energy range from 30 eV to 10,000 eV
C
C    Bugs :
C    if any are found please report to the authors
C
C
C    History :
C     04 jan 93 : original, based on Jelinsky's program written in
C                 C (EUVE Archive) - (19775::MBC)
C
C     23 feb 93 : comments added and modified to remove VAX
C                fortran 77 extensions (19775::MC)
C
C     21 sep 93 : further remnants of VAX fortran 77 extensions
C                  removed (19775::MBC)
C
C     22 sep 93 : bug fixed in calculations of EPSI - divided by 1.807317
C                 not added (19775::MBC)
C
C     16 jul 2000: modified by Joern Wilms (wilms@astro.uni-tuebingen.de)
c                 fourth argument is now energy in eV
c                 renamed to jwfano
C
C----------------------------------------------------------------------------
C
c     ... energy in Rydbergs
      REAL EPS
      REAL EPSI
c     ... log_10 of wavelength in Angstroms
      REAL X
C     Import :
c     ... Q coefficient: Fernley et al. 1987
c     ... nu, gamma coefficients: Oza 1986
      REAL Q,NU,GAMMA
c     ... energy in eV
      REAL ener
C     Export :
      REAL jwfano
c     ....
      EPS=ENER/13.598
      EPSI=3.0 - 1./(NU*NU) + 1.807317
      X=2.0*(EPS - EPSI)/GAMMA
      jwfano=(X-Q)*(X-Q)/(1.0+X*X)
c      write(*,*) (x-q)**2./(1.+X*X),(X-Q)*(X-Q)/(1.0+X*X)
      END
C-----------------------------------------------------------------------
      FUNCTION hphoto(keV1,keV2)
c
c     Wrapper to habs to return the photoelectric absorption cross section
c     in cm**2 of the hydrogen atom
c
c     calls habs below
c
c
      REAL hphoto,keV1,keV2
      REAL tbhabs
      EXTERNAL tbhabs

      hphoto=(tbhabs(1,keV1*1000.)+tbhabs(1,keV2*1000.))/2. * 1E-18
      return
      END

      FUNCTION tbhabs(z,e)
c
c     Photoionization cross-section for hydrogenic ions
c     after Band et al., 1990, A&A 237, 267
c
c     e: energy in eV
c     sig: cross-section in Mbarns
c
c     J. Wilms, wilms@astro.uni-tuebingen.de, 1999/04/21
c
      REAL tbhabs
      REAL e,e1s,y,denom
      INTEGER z
      
      e1s=13.599*z*z
      IF (e.lt.e1s) THEN
         tbhabs=0.
         return
      ENDIF 
      y=2.*e/e1s
      denom=z*z*y**(1.5)*(1.+sqrt(y))**4.
      tbhabs=605./denom
      
      return
      END 
c
c     simplified version of the SLATEC routines needed to compute the
c     gamma-function and incomplete gamma function (these special functions 
c     are needed, e.g., in the absorption models that contain grains).
c     The original subroutines can be obtained from netlib. The major
c     change has been to use fixed (i.e., non-machine-dependent) machine 
c     constants in function D1MACH. These values should be appropriate
c     for all machines on which XSPEC is supposed to run. Further
c     changes have been to REALLY simplify the error handling (which
c     is not yet adapted to the XSPEC way of doing things --> ask Keith!).
c
c     My comments are designated by "cjw".
c
c     Questions, Comments, Bug Reports to:
c        Joern Wilms, IAA Tuebingen, Astronomie, wilms@astro.uni-tuebingen.de
c
c


      DOUBLE PRECISION FUNCTION D1MACH (I)
cjw   customized 2000 March by Joern Wilms (wilms@astro.uni-tuebingen.de)
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
      DATA DMACH / 2.23D-308,1.79D+308,1.11D-16,2.22D-16,
     &     0.301029995663981195D0 /
***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH',
     +   'I OUT OF BOUNDS', 1, 2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
*DECK XERCLR
      SUBROUTINE XERCLR
C***BEGIN PROLOGUE  XERCLR
C***PURPOSE  Reset current error number to zero.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCLR-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        This routine simply resets the current error number to zero.
C        This may be necessary in order to determine that a certain
C        error has occurred again since the last time NUMXER was
C        referenced.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCLR
C***FIRST EXECUTABLE STATEMENT  XERCLR
      JUNK = J4SAVE(1,0,.TRUE.)
      RETURN
      END
*DECK D9GMIT
      DOUBLE PRECISION FUNCTION D9GMIT (A, X, ALGAP1, SGNGAM)
C***BEGIN PROLOGUE  D9GMIT
C***SUBSIDIARY
C***PURPOSE  Compute Tricomi's incomplete Gamma function for small
C            arguments.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      DOUBLE PRECISION (R9GMIT-S, D9GMIT-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
C             SPECIAL FUNCTIONS, TRICOMI
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute Tricomi's incomplete gamma function for small X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9GMIT
      DOUBLE PRECISION A, X, ALGAP1, SGNGAM, AE, AEPS, ALGS, ALG2,
     1  BOT, EPS, FK, S, SGNG2, T, TE, D1MACH, DLNGAM
      LOGICAL FIRST
      SAVE EPS, BOT, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9GMIT
      IF (FIRST) THEN
         EPS = 0.5D0*D1MACH(3)
         BOT = LOG (D1MACH(1))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'D9GMIT',
     +   'X SHOULD BE GT 0', 1, 2)
C
      MA = INT(A + 0.5D0)
      IF (A.LT.0.D0) MA = INT(A - 0.5D0)
      AEPS = A - MA
C
      AE = A
      IF (A.LT.(-0.5D0)) AE = AEPS
C
      T = 1.D0
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'D9GMIT',
     +   'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES', 2, 2)
C
 30   IF (A.GE.(-0.5D0)) THEN
         ALGS = -ALGAP1 + LOG(S)
         GO TO 60
      ENDIF
C
      ALGS = -DLNGAM(1.D0+AEPS) + LOG(S)
      S = 1.0D0
      M = -MA - 1
      IF (M.EQ.0) GO TO 50
      T = 1.0D0
      DO 40 K=1,M
        T = X*T/(AEPS-(M+1-K))
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
C
 50   D9GMIT = 0.0D0
      ALGS = -MA*LOG(X) + ALGS
      IF (S.EQ.0.D0 .OR. AEPS.EQ.0.D0) GO TO 60
C
      SGNG2 = SGNGAM * SIGN (1.0D0, S)
      ALG2 = -X - ALGAP1 + LOG(ABS(S))
C
      IF (ALG2.GT.BOT) D9GMIT = SGNG2 * EXP(ALG2)
      IF (ALGS.GT.BOT) D9GMIT = D9GMIT + EXP(ALGS)
      RETURN
C
 60   D9GMIT = EXP (ALGS)
      RETURN
C
      END
*DECK D9LGIC
      DOUBLE PRECISION FUNCTION D9LGIC (A, X, ALX)
C***BEGIN PROLOGUE  D9LGIC
C***SUBSIDIARY
C***PURPOSE  Compute the log complementary incomplete Gamma function
C            for large X and for A .LE. X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      DOUBLE PRECISION (R9LGIC-S, D9LGIC-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
C             LOGARITHM, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log complementary incomplete gamma function for large X
C and for A .LE. X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9LGIC
      DOUBLE PRECISION A, X, ALX, EPS, FK, P, R, S, T, XMA, XPA, D1MACH
      SAVE EPS
      DATA EPS / 0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9LGIC
      IF (EPS.EQ.0.D0) EPS = 0.5D0*D1MACH(3)
C
      XPA = X + 1.0D0 - A
      XMA = X - 1.D0 - A
C
      R = 0.D0
      P = 1.D0
      S = P
      DO 10 K=1,300
        FK = K
        T = FK*(A-FK)*(1.D0+R)
        R = -T/((XMA+2.D0*FK)*(XPA+2.D0*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'D9LGIC',
     +   'NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION', 1, 2)
C
 20   D9LGIC = A*ALX - X + LOG(S/XPA)
C
      RETURN
      END
*DECK D9LGIT
      DOUBLE PRECISION FUNCTION D9LGIT (A, X, ALGAP1)
C***BEGIN PROLOGUE  D9LGIT
C***SUBSIDIARY
C***PURPOSE  Compute the logarithm of Tricomi's incomplete Gamma
C            function with Perron's continued fraction for large X and
C            A .GE. X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      DOUBLE PRECISION (R9LGIT-S, D9LGIT-D)
C***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
C             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log of Tricomi's incomplete gamma function with Perron's
C continued fraction for large X and for A .GE. X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9LGIT
      DOUBLE PRECISION A, X, ALGAP1, AX, A1X, EPS, FK, HSTAR, P, R, S,
     1  SQEPS, T, D1MACH
      LOGICAL FIRST
      SAVE EPS, SQEPS, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9LGIT
      IF (FIRST) THEN
         EPS = 0.5D0*D1MACH(3)
         SQEPS = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0 .OR. A .LT. X) CALL XERMSG ('SLATEC', 'D9LGIT',
     +   'X SHOULD BE GT 0.0 AND LE A', 2, 2)
C
      AX = A + X
      A1X = AX + 1.0D0
      R = 0.D0
      P = 1.D0
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.D0+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'D9LGIT',
     +   'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 3, 2)
C
 30   HSTAR = 1.0D0 - X*S/A1X
      IF (HSTAR .LT. SQEPS) CALL XERMSG ('SLATEC', 'D9LGIT',
     +   'RESULT LESS THAN HALF PRECISION', 1, 1)
C
      D9LGIT = -X - ALGAP1 - LOG(HSTAR)
      RETURN
C
      END
*DECK D9LGMC
      DOUBLE PRECISION FUNCTION D9LGMC (X)
C***BEGIN PROLOGUE  D9LGMC
C***SUBSIDIARY
C***PURPOSE  Compute the log Gamma correction factor so that
C            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
C            + D9LGMC(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
C             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10. so that
C LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000E-02
C                                        with weighted error   1.28E-31
C                                         log weighted error  30.89
C                               significant figures required  29.81
C                                    decimal places required  31.48
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9LGMC
      DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS(  1) / +.1666389480 4518632472 0572965082 2 D+0      /
      DATA ALGMCS(  2) / -.1384948176 0675638407 3298605913 5 D-4      /
      DATA ALGMCS(  3) / +.9810825646 9247294261 5717154748 7 D-8      /
      DATA ALGMCS(  4) / -.1809129475 5724941942 6330626671 9 D-10     /
      DATA ALGMCS(  5) / +.6221098041 8926052271 2601554341 6 D-13     /
      DATA ALGMCS(  6) / -.3399615005 4177219443 0333059966 6 D-15     /
      DATA ALGMCS(  7) / +.2683181998 4826987489 5753884666 6 D-17     /
      DATA ALGMCS(  8) / -.2868042435 3346432841 4462239999 9 D-19     /
      DATA ALGMCS(  9) / +.3962837061 0464348036 7930666666 6 D-21     /
      DATA ALGMCS( 10) / -.6831888753 9857668701 1199999999 9 D-23     /
      DATA ALGMCS( 11) / +.1429227355 9424981475 7333333333 3 D-24     /
      DATA ALGMCS( 12) / -.3547598158 1010705471 9999999999 9 D-26     /
      DATA ALGMCS( 13) / +.1025680058 0104709120 0000000000 0 D-27     /
      DATA ALGMCS( 14) / -.3401102254 3167487999 9999999999 9 D-29     /
      DATA ALGMCS( 15) / +.1276642195 6300629333 3333333333 3 D-30     /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (FIRST) THEN
         NALGM = INITDS (ALGMCS, 15, REAL(D1MACH(3)) )
         XBIG = 1.0D0/SQRT(D1MACH(3))
         XMAX = EXP (MIN(LOG(D1MACH(2)/12.D0), -LOG(12.D0*D1MACH(1))))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 10.D0) CALL XERMSG ('SLATEC', 'D9LGMC',
     +   'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
C
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS,
     1  NALGM) / X
      RETURN
C
 20   D9LGMC = 0.D0
      CALL XERMSG ('SLATEC', 'D9LGMC', 'X SO BIG D9LGMC UNDERFLOWS', 2,
     +   1)
      RETURN
C
      END
*DECK DCSEVL
      DOUBLE PRECISION FUNCTION DCSEVL (X, CS, N)
C***BEGIN PROLOGUE  DCSEVL
C***PURPOSE  Evaluate a Chebyshev series.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
C***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Evaluate the N-term Chebyshev series CS at X.  Adapted from
C  a method presented in the paper by Broucke referenced below.
C
C       Input Arguments --
C  X    value at which the series is to be evaluated.
C  CS   array of N terms of a Chebyshev series.  In evaluating
C       CS, only half the first coefficient is summed.
C  N    number of terms in array CS.
C
C***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
C                 Chebyshev series, Algorithm 446, Communications of
C                 the A.C.M. 16, (1973) pp. 254-256.
C               L. Fox and I. B. Parker, Chebyshev Polynomials in
C                 Numerical Analysis, Oxford University Press, 1968,
C                 page 56.
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900329  Prologued revised extensively and code rewritten to allow
C           X to be slightly outside interval (-1,+1).  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCSEVL
      DOUBLE PRECISION B0, B1, B2, CS(*), ONEPL, TWOX, X, D1MACH
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DCSEVL
      IF (FIRST) ONEPL = 1.0D0 + D1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
C
      B1 = 0.0D0
      B0 = 0.0D0
      B2 = 0.0D0
      TWOX = 2.0D0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
C
      DCSEVL = 0.5D0*(B0-B2)
C
      RETURN
      END
*DECK DGAMI
      DOUBLE PRECISION FUNCTION DGAMI (A, X)
C***BEGIN PROLOGUE  DGAMI
C***PURPOSE  Evaluate the incomplete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      DOUBLE PRECISION (GAMI-S, DGAMI-D)
C***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate the incomplete gamma function defined by
C
C DGAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
C
C DGAMI is evaluated for positive values of A and non-negative values
C of X.  A slight deterioration of 2 or 3 digits accuracy will occur
C when DGAMI is very large or very small, because logarithmic variables
C are used.  The function and both arguments are double precision.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DGAMIT, DLNGAM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DGAMI
      DOUBLE PRECISION A, X, FACTOR, DLNGAM, DGAMIT
C***FIRST EXECUTABLE STATEMENT  DGAMI
      IF (A .LE. 0.D0) CALL XERMSG ('SLATEC', 'DGAMI',
     +   'A MUST BE GT ZERO', 1, 2)
      IF (X .LT. 0.D0) CALL XERMSG ('SLATEC', 'DGAMI',
     +   'X MUST BE GE ZERO', 2, 2)
C
      DGAMI = 0.D0
      IF (X.EQ.0.0D0) RETURN
C
C THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
      FACTOR = EXP (DLNGAM(A) + A*LOG(X))
C
      DGAMI = FACTOR * DGAMIT (A, X)
C
      RETURN
      END
*DECK DGAMIT
      DOUBLE PRECISION FUNCTION DGAMIT (A, X)
C***BEGIN PROLOGUE  DGAMIT
C***PURPOSE  Calculate Tricomi's form of the incomplete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      DOUBLE PRECISION (GAMIT-S, DGAMIT-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
C             SPECIAL FUNCTIONS, TRICOMI
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C   Evaluate Tricomi's incomplete Gamma function defined by
C
C   DGAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
C              T**(A-1.)
C
C   for A .GT. 0.0 and by analytic continuation for A .LE. 0.0.
C   GAMMA(X) is the complete gamma function of X.
C
C   DGAMIT is evaluated for arbitrary real values of A and for non-
C   negative values of X (even though DGAMIT is defined for X .LT.
C   0.0), except that for X = 0 and A .LE. 0.0, DGAMIT is infinite,
C   which is a fatal error.
C
C   The function and both arguments are DOUBLE PRECISION.
C
C   A slight deterioration of 2 or 3 digits accuracy will occur when
C   DGAMIT is very large or very small in absolute value, because log-
C   arithmic variables are used.  Also, if the parameter  A  is very
C   close to a negative integer (but not a negative integer), there is
C   a loss of accuracy, which is reported if the result is less than
C   half machine precision.
C
C***REFERENCES  W. Gautschi, A computational procedure for incomplete
C                 gamma functions, ACM Transactions on Mathematical
C                 Software 5, 4 (December 1979), pp. 466-481.
C               W. Gautschi, Incomplete gamma functions, Algorithm 542,
C                 ACM Transactions on Mathematical Software 5, 4
C                 (December 1979), pp. 482-489.
C***ROUTINES CALLED  D1MACH, D9GMIT, D9LGIC, D9LGIT, DGAMR, DLGAMS,
C                    DLNGAM, XERCLR, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
C***END PROLOGUE  DGAMIT
      DOUBLE PRECISION A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNG, ALX,
     1  BOT, H, SGA, SGNGAM, SQEPS, T, D1MACH, DGAMR, D9GMIT, D9LGIT,
     2  DLNGAM, D9LGIC
      LOGICAL FIRST
      SAVE ALNEPS, SQEPS, BOT, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DGAMIT
      IF (FIRST) THEN
         ALNEPS = -LOG (D1MACH(3))
         SQEPS = SQRT(D1MACH(4))
         BOT = LOG (D1MACH(1))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 0.D0) CALL XERMSG ('SLATEC', 'DGAMIT', 'X IS NEGATIVE'
     +   , 2, 2)
C
      IF (X.NE.0.D0) ALX = LOG (X)
      SGA = 1.0D0
      IF (A.NE.0.D0) SGA = SIGN (1.0D0, A)
      AINTA = AINT (A + 0.5D0*SGA)
      AEPS = A - AINTA
C
      IF (X.GT.0.D0) GO TO 20
      DGAMIT = 0.0D0
      IF (AINTA.GT.0.D0 .OR. AEPS.NE.0.D0) DGAMIT = DGAMR(A+1.0D0)
      RETURN
C
 20   IF (X.GT.1.D0) GO TO 30
      IF (A.GE.(-0.5D0) .OR. AEPS.NE.0.D0) CALL DLGAMS (A+1.0D0, ALGAP1,
     1  SGNGAM)
      DGAMIT = D9GMIT (A, X, ALGAP1, SGNGAM)
      RETURN
C
 30   IF (A.LT.X) GO TO 40
      T = D9LGIT (A, X, DLNGAM(A+1.0D0))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = EXP (T)
      RETURN
C
 40   ALNG = D9LGIC (A, X, ALX)
C
C EVALUATE DGAMIT IN TERMS OF LOG (DGAMIC (A, X))
C
      H = 1.0D0
      IF (AEPS.EQ.0.D0 .AND. AINTA.LE.0.D0) GO TO 50
C
      CALL DLGAMS (A+1.0D0, ALGAP1, SGNGAM)
      T = LOG (ABS(A)) + ALNG - ALGAP1
      IF (T.GT.ALNEPS) GO TO 60
C
      IF (T.GT.(-ALNEPS)) H = 1.0D0 - SGA * SGNGAM * EXP(T)
      IF (ABS(H).GT.SQEPS) GO TO 50
C
      CALL XERCLR
      CALL XERMSG ('SLATEC', 'DGAMIT', 'RESULT LT HALF PRECISION', 1,
     +   1)
C
 50   T = -A*ALX + LOG(ABS(H))
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = SIGN (EXP(T), H)
      RETURN
C
 60   T = T - A*ALX
      IF (T.LT.BOT) CALL XERCLR
      DGAMIT = -SGA * SGNGAM * EXP(T)
      RETURN
C
      END
*DECK DGAMLM
      SUBROUTINE DGAMLM (XMIN, XMAX)
C***BEGIN PROLOGUE  DGAMLM
C***PURPOSE  Compute the minimum and maximum bounds for the argument in
C            the Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A, R2
C***TYPE      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Calculate the minimum and maximum legal bounds for X in gamma(X).
C XMIN and XMAX are not the only bounds, but they are the only non-
C trivial ones to calculate.
C
C             Output Arguments --
C XMIN   double precision minimum legal value of X in gamma(X).  Any
C        smaller value of X might result in underflow.
C XMAX   double precision maximum legal value of X in gamma(X).  Any
C        larger value of X might cause overflow.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DGAMLM
      DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
C***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = LOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML)
     1    / (XMIN*XLN+0.5D0)
        IF (ABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMIN', 1, 2)
C
 20   XMIN = -XMIN + 0.01D0
C
      ALNBIG = LOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG)
     1    / (XMAX*XLN-0.5D0)
        IF (ABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMAX', 2, 2)
C
 40   XMAX = XMAX - 0.01D0
      XMIN = MAX (XMIN, -XMAX+1.D0)
C
      RETURN
      END
*DECK SLADGAMMA
      DOUBLE PRECISION FUNCTION SLADGAMMA (X)
C***BEGIN PROLOGUE  DGAMMA
C***PURPOSE  Compute the complete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      DOUBLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DGAMMA(X) calculates the double precision complete Gamma function
C for double precision argument X.
C
C Series for GAM        on the interval  0.          to  1.00000E+00
C                                        with weighted error   5.79E-32
C                                         log weighted error  31.24
C                               significant figures required  30.00
C                                    decimal places required  32.05
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9LGMC, DCSEVL, DGAMLM, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  DGAMMA
      DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX,
     1  XMIN, Y, D9LGMC, DCSEVL, D1MACH
      LOGICAL FIRST
C
      SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
      DATA GAMCS(  1) / +.8571195590 9893314219 2006239994 2 D-2      /
      DATA GAMCS(  2) / +.4415381324 8410067571 9131577165 2 D-2      /
      DATA GAMCS(  3) / +.5685043681 5993633786 3266458878 9 D-1      /
      DATA GAMCS(  4) / -.4219835396 4185605010 1250018662 4 D-2      /
      DATA GAMCS(  5) / +.1326808181 2124602205 8400679635 2 D-2      /
      DATA GAMCS(  6) / -.1893024529 7988804325 2394702388 6 D-3      /
      DATA GAMCS(  7) / +.3606925327 4412452565 7808221722 5 D-4      /
      DATA GAMCS(  8) / -.6056761904 4608642184 8554829036 5 D-5      /
      DATA GAMCS(  9) / +.1055829546 3022833447 3182350909 3 D-5      /
      DATA GAMCS( 10) / -.1811967365 5423840482 9185589116 6 D-6      /
      DATA GAMCS( 11) / +.3117724964 7153222777 9025459316 9 D-7      /
      DATA GAMCS( 12) / -.5354219639 0196871408 7408102434 7 D-8      /
      DATA GAMCS( 13) / +.9193275519 8595889468 8778682594 0 D-9      /
      DATA GAMCS( 14) / -.1577941280 2883397617 6742327395 3 D-9      /
      DATA GAMCS( 15) / +.2707980622 9349545432 6654043308 9 D-10     /
      DATA GAMCS( 16) / -.4646818653 8257301440 8166105893 3 D-11     /
      DATA GAMCS( 17) / +.7973350192 0074196564 6076717535 9 D-12     /
      DATA GAMCS( 18) / -.1368078209 8309160257 9949917230 9 D-12     /
      DATA GAMCS( 19) / +.2347319486 5638006572 3347177168 8 D-13     /
      DATA GAMCS( 20) / -.4027432614 9490669327 6657053469 9 D-14     /
      DATA GAMCS( 21) / +.6910051747 3721009121 3833697525 7 D-15     /
      DATA GAMCS( 22) / -.1185584500 2219929070 5238712619 2 D-15     /
      DATA GAMCS( 23) / +.2034148542 4963739552 0102605193 2 D-16     /
      DATA GAMCS( 24) / -.3490054341 7174058492 7401294910 8 D-17     /
      DATA GAMCS( 25) / +.5987993856 4853055671 3505106602 6 D-18     /
      DATA GAMCS( 26) / -.1027378057 8722280744 9006977843 1 D-18     /
      DATA GAMCS( 27) / +.1762702816 0605298249 4275966074 8 D-19     /
      DATA GAMCS( 28) / -.3024320653 7353062609 5877211204 2 D-20     /
      DATA GAMCS( 29) / +.5188914660 2183978397 1783355050 6 D-21     /
      DATA GAMCS( 30) / -.8902770842 4565766924 4925160106 6 D-22     /
      DATA GAMCS( 31) / +.1527474068 4933426022 7459689130 6 D-22     /
      DATA GAMCS( 32) / -.2620731256 1873629002 5732833279 9 D-23     /
      DATA GAMCS( 33) / +.4496464047 8305386703 3104657066 6 D-24     /
      DATA GAMCS( 34) / -.7714712731 3368779117 0390152533 3 D-25     /
      DATA GAMCS( 35) / +.1323635453 1260440364 8657271466 6 D-25     /
      DATA GAMCS( 36) / -.2270999412 9429288167 0231381333 3 D-26     /
      DATA GAMCS( 37) / +.3896418998 0039914493 2081663999 9 D-27     /
      DATA GAMCS( 38) / -.6685198115 1259533277 9212799999 9 D-28     /
      DATA GAMCS( 39) / +.1146998663 1400243843 4761386666 6 D-28     /
      DATA GAMCS( 40) / -.1967938586 3451346772 9510399999 9 D-29     /
      DATA GAMCS( 41) / +.3376448816 5853380903 3489066666 6 D-30     /
      DATA GAMCS( 42) / -.5793070335 7821357846 2549333333 3 D-31     /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DGAMMA
      IF (FIRST) THEN
         NGAM = INITDS (GAMCS, 42, 0.1*REAL(D1MACH(3)) )
C
         CALL DGAMLM (XMIN, XMAX)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.10.D0) GO TO 50
C
C COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
C GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
C
      N = INT(X)
      IF (X.LT.0.D0) N = N - 1
      Y = X - N
      N = N - 1
      SLADGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
C COMPUTE GAMMA(X) FOR X .LT. 1.0
C
      N = -N
      IF (X .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', 'X IS 0', 4, 2)
      IF (X .LT. 0.0 .AND. X+N-2 .EQ. 0.D0) CALL XERMSG ('SLATEC',
     +   'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
      IF (X .LT. (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X) .LT. DXREL)
     +   CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',
     +   1, 1)
C
      DO 20 I=1,N
        SLADGAMMA = SLADGAMMA/(X+I-1 )
 20   CONTINUE
      RETURN
C
C GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
C
 30   DO 40 I=1,N
        SLADGAMMA = (Y+I) * SLADGAMMA
 40   CONTINUE
      RETURN
C
C GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
C
 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'X SO BIG GAMMA OVERFLOWS', 3, 2)
C
      SLADGAMMA = 0.D0
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
      IF (X.LT.XMIN) RETURN
C
      SLADGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
C
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
     +   'DGAMMA',
     +   'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
C
      SINPIY = SIN (PI*Y)
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'X IS A NEGATIVE INTEGER', 4, 2)
C
      SLADGAMMA = -PI/(Y*SINPIY*SLADGAMMA)
C
      RETURN
      END
*DECK DGAMR
      DOUBLE PRECISION FUNCTION DGAMR (X)
C***BEGIN PROLOGUE  DGAMR
C***PURPOSE  Compute the reciprocal of the Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      DOUBLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
C***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DGAMR(X) calculates the double precision reciprocal of the
C complete Gamma function for double precision argument X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DGAMMA, DLGAMS, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DGAMR
      DOUBLE PRECISION X, ALNGX, SGNGX, SLADGAMMA
      EXTERNAL SLADGAMMA
C***FIRST EXECUTABLE STATEMENT  DGAMR
      DGAMR = 0.0D0
      IF (X.LE.0.0D0 .AND. AINT(X).EQ.X) RETURN
C
      CALL XGETF (IROLD)
      CALL XSETF (1)
      IF (ABS(X).GT.10.0D0) GO TO 10
      DGAMR = 1.0D0/SLADGAMMA(X)
      CALL XERCLR
      CALL XSETF (IROLD)
      RETURN
C
 10   CALL DLGAMS (X, ALNGX, SGNGX)
      CALL XERCLR
      CALL XSETF (IROLD)
      DGAMR = SGNGX * EXP(-ALNGX)
      RETURN
C
      END
*DECK DLGAMS
      SUBROUTINE DLGAMS (X, DLGAM, SGNGAM)
C***BEGIN PROLOGUE  DLGAMS
C***PURPOSE  Compute the logarithm of the absolute value of the Gamma
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      DOUBLE PRECISION (ALGAMS-S, DLGAMS-D)
C***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
C             FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DLGAMS(X,DLGAM,SGNGAM) calculates the double precision natural
C logarithm of the absolute value of the Gamma function for
C double precision argument X and stores the result in double
C precision argument DLGAM.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DLNGAM
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DLGAMS
      DOUBLE PRECISION X, DLGAM, SGNGAM, DLNGAM
C***FIRST EXECUTABLE STATEMENT  DLGAMS
      DLGAM = DLNGAM(X)
      SGNGAM = 1.0D0
      IF (X.GT.0.D0) RETURN
C
      IF (INT(MOD (-AINT(X), 2.0D0) + 0.1D0).EQ.0) SGNGAM = -1.0D0
C
      RETURN
      END
*DECK DLNGAM
      DOUBLE PRECISION FUNCTION DLNGAM (X)
C***BEGIN PROLOGUE  DLNGAM
C***PURPOSE  Compute the logarithm of the absolute value of the Gamma
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      DOUBLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
C***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DLNGAM(X) calculates the double precision logarithm of the
C absolute value of the Gamma function for double precision
C argument X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9LGMC, DGAMMA, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DLNGAM
      DOUBLE PRECISION X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX,
     1  Y, SLADGAMMA, D9LGMC, D1MACH, TEMP
      LOGICAL FIRST
      EXTERNAL SLADGAMMA
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA SQPI2L / +.2257913526 4472743236 3097614947 441 D+0    /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DLNGAM
      IF (FIRST) THEN
         TEMP = 1.D0/LOG(D1MACH(2))
         XMAX = TEMP*D1MACH(2)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS (X)
      IF (Y.GT.10.D0) GO TO 20
C
C LOG (ABS (DGAMMA(X)) ) FOR ABS(X) .LE. 10.0
C
      DLNGAM = LOG (ABS (SLADGAMMA(X)) )
      RETURN
C
C LOG ( ABS (DGAMMA(X)) ) FOR ABS(X) .GT. 10.0
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DLNGAM',
     +   'ABS(X) SO BIG DLNGAM OVERFLOWS', 2, 2)
C
      IF (X.GT.0.D0) THEN
         DLNGAM = SQ2PIL + (X-0.5D0)*LOG(X) - X + D9LGMC(Y)
         RETURN
      ENDIF
C
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DLNGAM',
     +   'X IS A NEGATIVE INTEGER', 3, 2)
C
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
     +   'DLNGAM',
     +   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',
     +   1, 1)
C
      DLNGAM = SQPI2L + (X-0.5D0)*LOG(Y) - X - LOG(SINPIY) - D9LGMC(Y)
      RETURN
C
      END
*DECK INITDS
      FUNCTION INITDS (OS, NOS, ETA)
C***BEGIN PROLOGUE  INITDS
C***PURPOSE  Determine the number of terms needed in an orthogonal
C            polynomial series so that it meets a specified accuracy.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
C***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
C             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Initialize the orthogonal series, represented by the array OS, so
C  that INITDS is the number of terms needed to insure the error is no
C  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
C  machine precision.
C
C             Input Arguments --
C   OS     double precision array of NOS coefficients in an orthogonal
C          series.
C   NOS    number of coefficients in OS.
C   ETA    single precision scalar containing requested accuracy of
C          series.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891115  Modified error message.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  INITDS
      DOUBLE PRECISION OS(*)
C***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'INITDS',
     +   'Number of coefficients is less than 1', 2, 1)
C
      ERR = 0.
      I = NOS
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
C
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITDS',
     +   'Chebyshev series too short for specified accuracy', 1, 1)
      INITDS = I
C
      RETURN
      END
*DECK XGETF
      SUBROUTINE XGETF (KONTRL)
C***BEGIN PROLOGUE  XGETF
C***PURPOSE  Return the current value of the error control flag.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETF-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C   Abstract
C        XGETF returns the current value of the error control flag
C        in KONTRL.  See subroutine XSETF for flag value meanings.
C        (KONTRL is an output parameter only.)
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XGETF
C***FIRST EXECUTABLE STATEMENT  XGETF
      KONTRL = J4SAVE(2,0,.FALSE.)
      RETURN
      END
*DECK XSETF
      SUBROUTINE XSETF (KONTRL)
C***BEGIN PROLOGUE  XSETF
C***PURPOSE  Set the error control flag.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3A
C***TYPE      ALL (XSETF-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XSETF sets the error control flag value to KONTRL.
C        (KONTRL is an input parameter only.)
C        The following table shows how each message is treated,
C        depending on the values of KONTRL and LEVEL.  (See XERMSG
C        for description of LEVEL.)
C
C        If KONTRL is zero or negative, no information other than the
C        message itself (including numeric values, if any) will be
C        printed.  If KONTRL is positive, introductory messages,
C        trace-backs, etc., will be printed in addition to the message.
C
C              ABS(KONTRL)
C        LEVEL        0              1              2
C        value
C          2        fatal          fatal          fatal
C
C          1     not printed      printed         fatal
C
C          0     not printed      printed        printed
C
C         -1     not printed      printed        printed
C                                  only           only
C                                  once           once
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Change call to XERRWV to XERMSG.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XSETF
      CHARACTER(8) XERN1
C***FIRST EXECUTABLE STATEMENT  XSETF
      IF (ABS(KONTRL) .GT. 2) THEN
         WRITE (XERN1, '(I8)') KONTRL
         CALL XERMSG ('SLATEC', 'XSETF',
     *      'INVALID ARGUMENT = ' // XERN1, 1, 2)
         RETURN
      ENDIF
C
      JUNK = J4SAVE(2,KONTRL,.TRUE.)
      RETURN
      END
