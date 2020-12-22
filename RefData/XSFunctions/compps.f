c Fortran routines used by the compps model (and called from xscompps.cxx)

c*********************************************************************
c********************************************************************* 
      SUBROUTINE ISMCO(PARM,OBJ,PHCON,PHBLB,PHREF,NE,PHNORM,FLAGFRE)
c********************************************************************* 
c this program calculates the photon flux 
c for the Compton scattered component (PHCON, normalized to 1 at 1 keV) 
c Reflected component (PHREF), and Black body (PHBLB)  
c for the disk-corona model neglecting energy balance 
c for the given parameters of the problem: 
c A: Maxwellian electrons 
c  TE (keV)  - electron temperature of the corona
c B: Power-law electrons
c  PLIND     - index of the power-law 
c  GMIN      - minimum Lorentz factor  
c  GMAX      - maximum Lorentz factor  
c C: Maxwellian + power-law tail
c  TE (keV)  - electron temperature of the corona
c  PLIND     - index of the power-law
c  GMIN      - Lorentz factor   of the break 
c  GMAX      - maximum Lorentz factor  
c 
c  TPH (keV) - temperature of the initial radiation (black-body) 
c  TAU       - vertical optical depth of  corona 
c  ETA       - cosine of the inclination angle 
c  FILFAC    - covering factor of cold material in the central plane
c  REDSHIFT  - redshift 
c  t0taug    - ratio of the height to the radius of the cylinder 
c********************************************************************* 
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      integer NE

      DOUBLE PRECISION DSPXX(MAXFRE),F2(MAXFRE),COMPTSP(MAXFRE)
      DOUBLE PRECISION BLBAVE(MAXFRE),REFLSP(MAXFRE)
      DOUBLE PRECISION PARM(11),FINT(MAXANG),pa0(11)
      DOUBLE PRECISION OBJ(NE)
      DOUBLE PRECISION PHCON(NE),PHBLB(NE),PHREF(NE)

      DOUBLE PRECISION xlmin, xlmax, eps, epsref, cmin, ckev, pi, phnorm
      DOUBLE PRECISION aj, cons, dely, eta, filf1, filfac, fly
      DOUBLE PRECISION hx, redshift, sumbb, sumdo, sumup
      DOUBLE PRECISION t0taug, tau, te, tph, yyy

      INTEGER neold, igeold, inpold
      INTEGER i, iang, j, jx, k, l, n, nl

      logical firstcall,flagfre,flagcs 


      save    pa0, firstcall, neold, igeold
      save    comptsp, blbave, reflsp

      data    pa0/11*9999./
      data    firstcall/.true./
      data    flagcs/.true./
      data    neold/0/,igeold/10/,inpold/10/  
      DATA CKEV/0.0019569472D0/,PI/3.141592653589792384D0/,EPS/1D-14/
      DATA CMIN/1d-30/
C XLMIN - decimal logarithms of low  photon energy
C               cutoff (in units of electron rest mass)
c XLMAX - maximum of the spectral range
c          0.1 eV         50 MeV 
      DATA XLMIN/-6.7D0/,XLMAX/2.D0/
      DATA EPSREF/1D-6/  
C ************************************
C ************************************

c EPSMIN accuracy of the spectrum computed 
c XEMAX  maximum energy up to what accuracy is checked 
      EPSMIN=3.D-3
      XEMAX=2.D0 

c number of scatterings 
cjp120504      MAXSC=50 
      MAXSC=50+INT(4d0*TAU**2) 
C
C other parameters: 
C ICHAND=0 isotropic distribution of initial radiation 
C ICHAND=1 Chandrasekhar-Sobolev distribution of intensity (1+2*mu)
C ICHAND=2 Isotropic flux I=1/mu 
      ICHAND=0 
C

      ETA=DABS(PARM(9))
      FILFAC=PARM(10)
      REDSHIFT=PARM(11)

C*****************************************
C IDISTF=0 Maxwellian electron distribution 
C IDISTF=1 monoenergetic electron distribution 
C IDISTF=2 power-law  electron distribution 
C IDISTF=3 cutoff Maxwellian 
C IDISTF=4 Maxwellian + power-law tail

C default is Maxellian
      IDISTF = 0
      if (PARM(3) .GE. 1.0 .AND. PARM(4) .GE. 1.0) THEN
c cut-off Maxwellian
         if (PARM(3) .GE. PARM(4) ) THEN
            IDISTF = 3
         else
            if ( ABS(PARM(1)) .LE. 1.0 ) THEN
c power-law electrons
               IDISTF = 2
            else
c hybrid
               IDISTF = 4
            endif
         endif
      endif

      if ( IDISTF .EQ. 0 .OR. IDISTF .EQ. 3 ) THEN
         if ( ABS(PARM(1)) .LE. 1.0 ) THEN
            CALL xwrite("COMPPS: kTe too small",5)
            RETURN
         endif
      endif

C***************************************

c set up electron distribution internal variables
c This line added to suppress a compiler warning
      TE=0.0

      IF(IDISTF.eq.0.or.IDISTF.eq.3.or.IDISTF.eq.4) THEN
         TE=DABS(PARM(1)/511.)
         YE=1D0/TE 
      ENDIF

      IF(IDISTF.eq.2) THEN 
         PLIND=DABS(PARM(2))
         GMIN=PARM(3)
         GMAX=PARM(4)
      ELSEIF(IDISTF.eq.3) THEN
         PLIND=PARM(2) 
         GMIN=PARM(3)
      ELSEIF(IDISTF.eq.4) THEN
         PLIND=PARM(2) 
         GMIN=PARM(3)
         GMAX=PARM(4) 
      ENDIF 

      TPH=DABS(PARM(5)) 
      IF(PARM(6).gt.0D0.or.IDISTF.eq.2) THEN
         TAU=DABS(PARM(6)) 
      ELSE
c tau=y/(4*Theta) 
         TAU=DABS(PARM(6)/(4D0*TE)) 
      ENDIF 
C *****************************************       
      YPH=511D0/TPH
C *****************************************       
c**********************************
c normalization using bbody emission 
c 1.0344E-3 is taken from xsbbrd.f 
c phnorm=1/[Tbb (keV)**3 1.0344E-3 / 511] since we divide by PHNORM
      PHNORM=511D0/1.0344D-3/TPH**3
c**********************************

C IGEOM=0 escape probability (sphere) 
C IGEOM=1 SLAB geometry
C IGEOM=2 CYLINDER geometry (with T0TAUG=H/R)
C IGEOM=3 HEMISPHERE geometry  (TOTAUG=1, H=R) 
C IGEOM=4,5 SPHERICAL  geometry 

c IGEOM.gt.0 sources at the bottom (or central plane) 
c IGEOM.lt.0 homogeneous isotropic sources 
c IGEOM=-5    source of photons distributed
c       according to the eigen function of the diffusion equation
c       f(tau')=sim(pi*tau'/tau)/(pi*tau'/tau)
c       where tau' varies between 0 and tau.

      IGEOM = NINT(PARM(7))
      IF ( IGEOM .GT. 5 ) THEN
         CALL xwrite("COMPPS: geometry parameter out of range",5)
         RETURN
      ENDIF

c set up geometry internal variables

c sources at the bottom (or central plane) or center of the sphere
      if(IGEOM.gt.0) then
         INPUT=0
      else
c sources  are homogeneous and isotropic or sim(pi*tau'/tau)/(pi*tau'/tau)
         INPUT=1      
      endif
      IGEOM = ABS(IGEOM)

      if(IGEOM.eq.4.or.IGEOM.eq.5) FILFAC=1D0 
      FILF1=1D0-FILFAC 
      if(DABS(FILF1).lt.EPSREF) THEN 
         IGROUND=0
      else 
         IGROUND=2
      endif 
      
c T0TAUG=H/R 
c for slab T0TAUG=0; for cylinder =parm(10),  for hemisphere =1
      IF(IGEOM.eq.1) T0TAUG=0.
      IF(IGEOM.eq.2) T0TAUG=PARM(8) 
      IF(IGEOM.eq.3) T0TAUG=1.
C IREDF =0 exact redistribution function    
C IREDF =1 isotropic approximation in the electron rest frame 
c NOW it is always IREDF=0 
c      IF(PARM(1).gt.0.) THEN 
      IREDF=0
c      ELSE 
c      IREDF=1
c      ENDIF  
C NISO=0 - exact scattering process
C NISO=1 - anisotropy for first 1,..ISOTR scatterings,
C          for all other, isotropic TAU dependent source function
      NISO=1
cjp000110
c      IF(PARM(6).gt.0.) THEN 
c      ISOTR=5 
cjp 
c      ELSE 
      ISOTR=20 
c      ENDIF  
      if(INPUT.eq.0.and.ISOTR.lt.1) ISOTR=1
      if(IGEOM.eq.0) then 
         NISO=1
         ISOTR=0 
      endif 

c Compton spectral parameters have changed  
      flagcs = .false.
      do i = 1, 10
         if ( pa0(i) .NE. parm(i) ) flagcs = .true.
      enddo

      if(igeom.ne.igeold.or.input.ne.inpold) then 
c for the central injection - Gauss-Chebyshev angular quadrature 
         if(INPUT.eq.0.and.(IGEOM.eq.4.or.IGEOM.eq.5)) then 
            CALL QCHEB(UANG,AANG,MAXANG)
         else 
c for any  other  injection - Gauss angular quadrature 
            CALL QDRGSDO(UANG,AANG,MAXANG)  
            DO 20 L=1,MAXANG
               UANG(L)=0.5D0*(1D0+UANG(L))                         
               AANG(L)=AANG(L)*5D-1                                
 20         CONTINUE 
         endif 
         flagcs = .true.
      endif 

c set up internal arrays on first invocation of function

      if(firstcall) then 
C NL    -  number of points in Gauss-Laguerre quadrature  
         NL=10
         CALL AGUEGALA(NL,EPS) 
         CALL QDRGSDO(UG,AG,NG)  
         DO 23 L=1,NG
            UG(L)=0.5D0*(1D0+UG(L))                         
            AG(L)=AG(L)*5D-1                                
 23      CONTINUE 
         CALL QDRGSDO(UGS,AGS,NS)
         DO 26 L=1,NS
            UGS(L)=0.5D0*(1D0+UGS(L))
            AGS(L)=AGS(L)*5D-1
 26      CONTINUE
         CALL QDRGSDO(UGSPH,AGSPH,NSPH)
         DO 28 L=1,NSPH
            UGSPH(L)=0.5D0*(1D0+UGSPH(L))
            AGSPH(L)=AGSPH(L)*5D-1
 28      CONTINUE   
c frequency step 
         HX=(XLMAX-XLMIN)/DFLOAT(MAXFRE-1)
         CONS=DLOG(1D1)*HX 
C calculation of frequencies
         CALL FREQPOIN(XLMIN,HX)
C weights of frequency integration 
         DO 27 JX=1,MAXFRE
            AJ=1D0
            IF(JX .eq. 1 .or. JX .eq. MAXFRE) AJ=5D-1 
            A(JX)=AJ*CONS
 27      CONTINUE
c weights of frequency-angular integration 
         DO N=1,ND
            DO I=1,MAXFRE
               DO J=1,MAXANG
                  K=J+(I-1)*MAXANG 
                  L=(N-1)*KK+K
c     AINT and CINT are FREQUENCY and ANGLE  supervectors.
                  AINT(L)=A(I)
                  CINT(L)=AANG(J)
                  AC(L)=A(I)*AANG(J) 
               ENDDO
            ENDDO
         ENDDO
         firstcall=.false.
      endif 


C******************************************************
C   if necessary recalculate Compton scattering       *
C******************************************************
      if(pa0(1).ne.parm(1).or.pa0(2).ne.parm(2).or.
     &     pa0(3).ne.parm(3).or.pa0(4).ne.parm(4))then
         CALL COMWY(EPS) 
C Compton cross-section 
         CALL CROSCOMP 
c COMPTON scattering redistribution matrix 
         CALL CSRF 
      endif 
      
C******************************************************
C   if necessary recalculate blackbody                *
C******************************************************
C blackbody intrinsic radiation (isotropic incident radiation) 
      if(pa0(5).ne.parm(5)) then
         if(PARM(5).ge.0d0) CALL BLACKB(XEN,DINTEN,MAXFRE) 
         if(PARM(5).lt.0d0) CALL BBMULTI(XEN,DINTEN,MAXFRE) 
      endif 

C******************************************************
C   if necessary recalculate comptonized spectrum     *
c******************************************************
C calculation of the comptonized spectrum
      if(flagcs) then
         CALL COMPSPEC(TAU,T0TAUG,FILF1)
c flux averaged over angles
         DO 212 I=1,MAXFRE
            sumdo=0d0
            sumup=0d0
            sumbb=0d0
            DO 210 IANG=1,MAXANG
               K=IANG+(I-1)*MAXANG     
c flux down
               sumdo=sumdo+aang(iang)*SUIMI(K)*uang(iang)
c flux up 
               sumup=sumup+aang(iang)*SUIPL(K)*uang(iang)
               sumbb=sumbb+aang(iang)*DINTPL(K)*uang(iang)
 210        CONTINUE
            if(igeom.eq.1.or.igeom.eq.2.or.igeom.eq.3) then 
c slab + hemisphere + cylinder : flux down 
               REFLSP(I)=sumdo  
            else 
c sphere : flux up  
               REFLSP(I)=sumup 
            endif 
            BLBAVE(I)=sumbb
 212     CONTINUE

      endif 


C ****************************************************
C case of slab, hemisphere or cylinder
      if (flagcs) then
         if (igeom.eq.1.or.igeom.eq.2.or.igeom.eq.3) then 

c calculation of COMPTON SCATTERED FLUX at obs. grid 
c interpolation between angles to get N_nu=I_nu/nu at ETA

            DO I=1,MAXFRE
               DO IANG=1,MAXANG 
                  K=IANG+(I-1)*MAXANG 
                  FINT(IANG)=DLOG(SUIPL(K)*UANG(IANG)+CMIN)
               ENDDO
               CALL POLINTQ(UANG,FINT,MAXANG,ETA,FLY,DELY)
               COMPTSP(I)=DEXP(FLY)
            ENDDO

c calculation of BLACK BODY FLUX at obs. grid 
c interpolation between angles to get N_nu=I_nu/nu at ETA

            DO I=1,MAXFRE
               DO IANG=1,MAXANG 
                  K=IANG+(I-1)*MAXANG 
                  FINT(IANG)=DLOG(DINTPL(K)*UANG(IANG)+CMIN)  
               ENDDO
               CALL POLINTQ(UANG,FINT,MAXANG,ETA,FLY,DELY)
               BLBAVE(I)=DEXP(FLY)
            ENDDO
            
         ELSE

c otherwise just copy REFLSP into COMPTSP

            DO I=1, MAXFRE
               COMPTSP(I) = REFLSP(I)
            ENDDO

         ENDIF

      ENDIF

C******************************************************
C   remap from internal energies                      *
C******************************************************

      if(flagcs.or.flagfre) then
c calculation of COMPTON SCATTERED FLUX at obs. grid 
         DO 170 I=1,MAXFRE
            DSPXX(I)=DLOG(COMPTSP(I)+CMIN)  
 170     CONTINUE 
         CALL QSPLINE(XLOG,DSPXX,MAXFRE,1d30,1d30,F2) 
         DO 172 I=1,NE 
            IF(OBJ(I).lt.XEN(1).or.OBJ(I).ge.XEN(MAXFRE)) THEN 
               PHCON(I)=CMIN 
            ELSE 
               CALL QSPLINT(XLOG,DSPXX,F2,MAXFRE,DLOG10(OBJ(I)),YYY)
               PHCON(I)=DEXP(YYY)/OBJ(I) 
            ENDIF 
 172     CONTINUE 
c calculation of input for reflection at obs. grid 
         DO I=1,MAXFRE
            DSPXX(I)=DLOG(REFLSP(I)+CMIN)  
         ENDDO
         CALL QSPLINE(XLOG,DSPXX,MAXFRE,1d30,1d30,F2) 
         DO I=1,NE 
            IF(OBJ(I).lt.XEN(1).or.OBJ(I).ge.XEN(MAXFRE)) THEN 
               PHREF(I)=CMIN 
            ELSE 
               CALL QSPLINT(XLOG,DSPXX,F2,MAXFRE,DLOG10(OBJ(I)),YYY)
               PHREF(I)=DEXP(YYY)/OBJ(I) 
            ENDIF 
         ENDDO
c calculation of BLACK BODY FLUX at obs. grid 
         DO 174 I=1,MAXFRE 
            DSPXX(I)=DLOG(BLBAVE(I)+CMIN)  
 174     CONTINUE 
         CALL QSPLINE(XLOG,DSPXX,MAXFRE,1d30,1d30,F2) 
         DO 176 I=1,NE 
            IF(OBJ(I).lt.XEN(1).or.OBJ(I).ge.XEN(MAXFRE)) THEN 
               PHBLB(I)=CMIN 
            ELSE 
               CALL QSPLINT(XLOG,DSPXX,F2,MAXFRE,DLOG10(OBJ(I)),YYY)
               PHBLB(I)=DEXP(YYY)/OBJ(I)  
            ENDIF 
 176     CONTINUE          
      endif

      igeold=igeom 
      inpold=input

      do n=1,11
         pa0(n)=parm(n)
      enddo
      RETURN 
      END                  

C**********************************************************************
C
C IGEOM =  0 escape probability (sphere) 
C IGEOM =  1 SLAB     
C IGEOM =  2 CYLINDER 
C IGEOM =  3 HEMISPHERE 
C IGEOM =  4,5 SPHERE
C integration over the optical depth
C iteration prosedure:  Neuman series 
C                       (separation of scattering orders) 
C 
      SUBROUTINE COMPSPEC(TAU,T0TAUG,FILF1)
      IMPLICIT NONE

      DOUBLE PRECISION TAU, T0TAUG, FILF1

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION SOUPL(MAXTAU,LL),SOUMI(MAXTAU,LL), 
     2 BGPL1(MAXTAU,MAXTAU,KK),BGMI1(MAXTAU,MAXTAU,KK),   
     3 BGPL2(MAXTAU,MAXTAU,KK),BGMI2(MAXTAU,MAXTAU,KK), 
     4 BGMIN(MAXTAU,MAXTAU,KK),BGMEX(MAXTAU,MAXTAU,KK), 
     5 BGPLA(MAXTAU,MAXTAU,MAXFRE),BGMIA(MAXTAU,MAXTAU,MAXFRE),
     6 BGMINA(MAXTAU,MAXTAU,MAXFRE)
      DOUBLE PRECISION TAUVEC(MAXTAU),SOUPL0(MAXTAU,LL),
     1 SOUMI0(MAXTAU,LL), DITAPL0(LL,MAXTAU),DITAMI0(LL,MAXTAU) 
      DOUBLE PRECISION SOISPL(MAXTAU,MAXFRE),SOISMI(MAXTAU,MAXFRE),
     1 DIAVER(MAXFRE,MAXTAU), DIFLUX(MAXFRE,MAXTAU) 
      DOUBLE PRECISION DINUP(LL),DINDO(LL)
      DOUBLE PRECISION CINP(MAXTAU,MAXTAU,MAXANG),
     1         CINM(MAXTAU,MAXTAU,MAXANG),CEMP(MAXTAU,MAXTAU,MAXANG),
     2         CEMM(MAXTAU,MAXTAU,MAXANG),CGRIN(MAXTAU,MAXTAU,MAXANG),
     3         CGREX(MAXTAU,MAXTAU,MAXANG)

      DOUBLE PRECISION T0TAUGOLD,TAUOLD,YOLD,FIL1OLD,YPHOLD,PI
      DOUBLE PRECISION ATAU, DELTAU, DIFMAX, DINSTA, DINT, DINTAU
      DOUBLE PRECISION DLN, DNORM, DSTA, DTT, FLPL, SINPI, DNTR
      DOUBLE PRECISION STA, sumneg, sumpos, T0TA2, TAUPI, TAURE
      DOUBLE PRECISION TAUSC, XPI

      INTEGER IGEOLD,IGROLD
      INTEGER IANG, ISC, ITAU, IX, JX, K, L

      DOUBLE PRECISION SINXX
      EXTERNAL SINXX

      save T0TAUGOLD,TAUOLD,YOLD,FIL1OLD,IGEOLD,IGROLD,YPHOLD
      save BGPL1,BGMI1,BGPL2,BGMI2,BGMIN,BGMEX
      save BGPLA,BGMIA,BGMINA 
      save CINP,CINM,CEMP,CEMM,CGRIN,CGREX
     
      DATA PI/3.141592653589792384D0/
      data TAUOLD/9D9/,YOLD/-9D9/,FIL1OLD/9D9/,YPHOLD/-9D9/ 
      data T0TAUGOLD/9D9/,IGEOLD/10/,IGROLD/10/

      DO 377 K=1,KK 
         SUIPL(K)=0D0 
         SUIMI(K)=0D0 
         DINTPL(K)=0D0
         DINTMI(K)=0D0
 377  CONTINUE 

      T0TA2=T0TAUG/2D0
      DNTR=1D0/DFLOAT(MAXTAU-1)
      DELTAU=TAU*DNTR 
      DLN=DLOG(1D1)

      DTT=0D0
      DO 4 ITAU=1,MAXTAU
         TAUVEC(ITAU)=DTT
         DTT=DTT+DNTR
    4 CONTINUE 

c if not a SPHERE then
c compute geometrical factors 
c      if(IGEOM.ne.4) then
      if(IGEOM.eq.1.or.IGEOM.eq.2.or.IGEOM.eq.3) then
         if(IGEOM.ne.IGEOLD.or.T0TAUG.ne.T0TAUGOLD) then 
            CALL GEOMFAC(TAUVEC,T0TAUG,CINP,CINM,CEMP,CEMM)
         endif
c for disk with holes            
         if(IGROUND.eq.2.and.(IGEOM.ne.IGEOLD.or.IGROUND.ne.IGROLD 
     &        .or.T0TAUG.ne.T0TAUGOLD)) then 
            CALL GEOMGR(TAUVEC,T0TAUG,CGRIN,CGREX)
         endif

c computation of energy and optical depth dependent factors
         if(IGEOM.ne.IGEOLD.or.YOLD.ne.YE.or.TAUOLD.ne.TAU
     &        .or.T0TAUG.ne.T0TAUGOLD) then 
            CALL  ENERFAC(TAUVEC,TAU,T0TAUG,CINP,CINM,CEMP,CEMM,
     &           BGPL1,BGMI1,BGPL2,BGMI2,BGPLA,BGMIA)
         endif
c for disk with holes      
         if(IGROUND.eq.2.and.(IGEOM.ne.IGEOLD.or.
     &        IGROUND.ne.IGROLD.or.YOLD.ne.YE.or.TAUOLD.ne.TAU
     &        .or.T0TAUG.ne.T0TAUGOLD)) then 
            CALL  ENERGR(TAUVEC,TAU,CGRIN,CGREX,
     &           BGMIN,BGMEX,BGMINA)
         endif
c end IGEOM.ne.4,5      
      endif

      ISC=0
         
C SOURCE FUNCTION FOR UNSCATTERED RADIATION
      IF(INPUT.eq.0) THEN
c S_0(tau)=delta(tau)      
         if(IGEOM.eq.4.or.IGEOM.eq.5) then
c sphere
            DO IX=1,MAXFRE
               XPI=XEN(IX)*PI
               sta=TAU*S0(IX)
               dsta=0d0
               if(sta.lt.7d1) dsta=dexp(-sta)*DINTEN(IX)
               DO IANG=1,MAXANG
                  K=IANG+(IX-1)*MAXANG 
c emergent black body radiation
                  DINTPL(K)=dsta
c      DO 10 ITAU=1,NT
                  SOUPL(1,K)=0d0
                  SOUMI(1,K)=0d0
                  SOUPL0(1,K)=0d0
                  SOUMI0(1,K)=0d0
                  DO ITAU=2,MAXTAU
                     TAUSC=XPI/(TAUVEC(ITAU)*TAUVEC(ITAU))
                     TAURE=TAU*TAUVEC(ITAU)
c integration over frequency  
                     sumneg=0d0    
                     sumpos=0d0  
                     do JX=1,MAXFRE
                        sta=TAURE*S0(JX)
                        dinsta=0d0
                        if(sta.lt.7d1) dinsta=dexp(-sta)*DINTEN(JX)
                        sumneg=sumneg+FRDMU(JX,IX,MAXANG-IANG+1)
     &                       *dinsta*A(JX)
                        sumpos=sumpos+FRDMU(JX,IX,IANG+MAXANG)
     &                       *dinsta*A(JX)
                     enddo
                     SOUPL(ITAU,K)=sumpos*TAUSC
                     SOUMI(ITAU,K)=sumneg*TAUSC
                     SOUPL0(ITAU,K)=SOUPL(ITAU,K)
                     SOUMI0(ITAU,K)=SOUMI(ITAU,K)
                  ENDDO
               ENDDO
            ENDDO
      
c not a sphere  
         else
            DO IX=1,MAXFRE
               DINT=DINTEN(IX)
               DO IANG=1,MAXANG
                  K=IANG+(IX-1)*MAXANG 
                  DO 22 ITAU=1,MAXTAU
                     SOUPL(ITAU,K)=0D0
                     SOUMI(ITAU,K)=0D0
                     SOUPL0(ITAU,K)=0D0
                     SOUMI0(ITAU,K)=0D0
 22               CONTINUE 
                  SOUPL(1,K)=0D0
                  SOUPL0(1,K)=UANG(IANG)*DINT*2D0/DELTAU
c emergent black body radiation for delta(tau) injection (not sphere) 
                  FLPL=0D0 
                  DO 24 ITAU=1,MAXTAU
                     ATAU=1D0
                     IF(ITAU.eq.1.or.ITAU.eq.MAXTAU) ATAU=5D-1  
                     FLPL=FLPL+ATAU*BGPL2(1,ITAU,K)
 24               CONTINUE 
                  DINTPL(K)=FLPL*DNTR/UANG(IANG)*DINT 
               ENDDO
            ENDDO
         endif
      ELSE
c INPUT=1 
         if(IGEOM.ne.5) then
c S_0(tau) homogeneous
            DO IX=1,MAXFRE
               DINT=DINTEN(IX)*S0(IX) 
               if(IGEOM.eq.0) DINT=DINTEN(IX)
               DO ITAU=1,MAXTAU
                  SOISPL(ITAU,IX)=DINT
                  SOISMI(ITAU,IX)=DINT
                  DO IANG=1,MAXANG
                     K=IANG+(IX-1)*MAXANG       
                     SOUPL(ITAU,K)=0D0
                     SOUMI(ITAU,K)=0D0
                     SOUPL0(ITAU,K)=DINT
                     SOUMI0(ITAU,K)=DINT
                  ENDDO
               ENDDO
            ENDDO
            IF(IGEOM.eq.2.or.IGEOM.eq.3) THEN
               CALL SIDEESC(SOUPL0,SOUMI0,DELTAU,DNTR,FILF1,BGMEX,BGPL2,
     &              BGMI2,DINUP,DINDO)
               DO 31 L=1,LL
                  DINTPL(L)=DINTPL(L)+DINUP(L)
 31            CONTINUE
            ENDIF
         else
c S_0(tau) propto sin(pi*tau'/tau)/(pi*tau'/tau)  
            DO IX=1,MAXFRE
               DINT=DINTEN(IX)*S0(IX)
               DO ITAU=1,MAXTAU
c tp=pi*tau'/tau 
                  TAUPI=PI*TAUVEC(ITAU) 
c sin(tp)/tp 
                  SINPI=SINXX(TAUPI) 
                  DINTAU=DINT*SINPI 
                  SOISPL(ITAU,IX)=DINTAU
                  SOISMI(ITAU,IX)=DINTAU 
                  DO IANG=1,MAXANG
                     K=IANG+(IX-1)*MAXANG
                     SOUPL(ITAU,K)=0D0
                     SOUMI(ITAU,K)=0D0
                     SOUPL0(ITAU,K)=DINTAU 
                     SOUMI0(ITAU,K)=DINTAU 
                  ENDDO
               ENDDO
            ENDDO
         endif 
      ENDIF
C          
  
      ISC=1 
      DIFMAX=1D0 
C MULTIPLY SCATTERING: ITERATION PROCEDURE 
      DO WHILE(ISC .lt. MAXSC .and. DIFMAX .gt. EPSMIN) 
         IF(ISC.gt.1) DIFMAX=0D0 

c not escape probability 
         if(IGEOM.ne.0) then
c INTERNAL RADIATION FIELD 
c EXACT INTENSITIES
            IF(NISO.eq.0.or.ISC.le.ISOTR) THEN 
               if(IGEOM.eq.4.or.IGEOM.eq.5)  then 
c sphere 
                  CALL EXASPH(SOUPL0,SOUMI0,DITAPL0,DITAMI0,DIAVER,
     &                 DIFLUX,TAUVEC,TAU,ISC)
               else
                  CALL EXAINT(SOUPL0,SOUMI0,DITAPL0,DITAMI0,DIAVER,
     &                 DIFLUX,BGPL1,BGMI1,BGMIN,DELTAU,FILF1,ISC)
               endif
            ELSE 
c ISOTROPIC INTENSITIES    
               if(IGEOM.eq.4.or.IGEOM.eq.5)  then
c sphere 
                  CALL ISOSPH(SOISPL,SOISMI,DIAVER,DIFLUX,TAUVEC,
     &                        TAU,ISC)
               else 
                  CALL ISOINT(SOISPL,SOISMI,DIAVER,DIFLUX, 
     &                 BGPL1,BGMI1,BGMIN,BGMINA,BGPLA,BGMIA,DELTAU,
     &                 FILF1,ISC)
               endif      
            ENDIF 

c EXACT SOURCE FUNCTIONS
            IF(NISO.eq.0.or.ISC.lt.ISOTR) THEN 
               CALL EXASOUR(SOUPL,SOUMI,DITAPL0,DITAMI0,SOUPL0,SOUMI0,
     &                      DIFMAX)     
            ELSE
c ISOTROPIC SOURCE FUNCTIONS dependent on TAU
               CALL ISOSOUR(SOUPL,SOUMI,DIAVER,DIFLUX,SOISPL,SOISMI,
     &                      DIFMAX)
            ENDIF

c end of not escape probability 
         else
c escape probability (for sphere) 
            CALL ESCINT(SOISPL,DIAVER,TAU,ISC)
            CALL ESCSOUR(SOUPL,DIAVER,SOISPL,DIFMAX) 
         endif 

         ISC=ISC+1 
      
      ENDDO       
c end of iteration procedure 
c***************************
      
c emergent Comptonized radiation field through sides 
      IF(IGEOM.eq.2.or.IGEOM.eq.3) THEN 
         CALL SIDEESC(SOUPL,SOUMI,DELTAU,DNTR,FILF1,BGMEX,BGPL2,
     &                BGMI2,DINUP,DINDO)
         DO 380 L=1,LL 
            SUIPL(L)=SUIPL(L)+DINUP(L) 
            SUIMI(L)=SUIMI(L)+DINDO(L) 
 380     CONTINUE
      ENDIF 


      DNORM=1d0 
      if(INPUT.eq.1) then 
         if(IGEOM.eq.1) DNORM=1D0/(TAU*(1D0+FILF1)) 
         if(IGEOM.eq.2) DNORM=1D0/(TAU*(1D0+FILF1)) 
         if(IGEOM.eq.3) DNORM=0.5*PI/(TAU*(1D0+FILF1)) 
         if(IGEOM.eq.0) DNORM=2D0 
         if(IGEOM.eq.4) DNORM=0.75D0/TAU*2d0  
         if(IGEOM.eq.5) DNORM=PI*PI*0.25D0/TAU*2d0  
      endif 
      if(INPUT.eq.0) DNORM=2d0 
      do k=1,kk 
         DINTPL(K)=DINTPL(K)*DNORM
         SUIPL(K)=SUIPL(K)*DNORM 
         SUIMI(K)=SUIMI(K)*DNORM 
      enddo


      T0TAUGOLD=T0TAUG 
      TAUOLD=TAU 
      YPHOLD=YPH 
      YOLD=YE
      FIL1OLD=FILF1 
      IGEOLD=IGEOM 
      IGROLD=IGROUND
      RETURN 
      END

c***********************************************************************
c ESCAPE through the sides of hemisphere or cylinder 
      SUBROUTINE SIDEESC(SOUPL,SOUMI,DELTAU,DNTR,FILF1,BGMEX,BGPL2,
     &                   BGMI2,DINUP,DINDO) 
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION SOUPL(MAXTAU,LL),SOUMI(MAXTAU,LL),DINUP(LL)
      DOUBLE PRECISION DINDO(LL),BGMEX(MAXTAU,MAXTAU,KK) 
      DOUBLE PRECISION BGPL2(MAXTAU,MAXTAU,KK),BGMI2(MAXTAU,MAXTAU,KK)
      DOUBLE PRECISION DELTAU,DNTR,FILF1

      DOUBLE PRECISION AT, ATAU, EPTA, ETA, FLMI, FLPL, SUMI2, SUPL2

      INTEGER IANG, IT, ITAU, IX, K, L, LTAU

      DO IX=1,MAXFRE
         DO IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            L=K
            ETA=UANG(IANG)
            EPTA=DELTAU/ETA
            FLPL=0D0
            FLMI=0D0
            DO 360 ITAU=1,MAXTAU
               LTAU=MAXTAU-ITAU+1
               ATAU=1D0
               IF(ITAU.eq.1.or.ITAU.eq.MAXTAU) ATAU=5D-1

               SUPL2=0D0
               IF(IGROUND.eq.2) THEN
                  DO 362 IT=1,MAXTAU
                     AT=1D0
                     IF(IT.eq.1.or.IT.eq.MAXTAU) AT=5D-1
                     SUPL2=SUPL2+AT*SOUMI(IT,L)*BGMEX(IT,ITAU,K)
 362              CONTINUE
                  SUPL2=SUPL2*FILF1
               ENDIF
               DO 364 IT=1,ITAU
                  AT=1D0
                  IF(IT.eq.1.or.IT.eq.ITAU) AT=5D-1
                  SUPL2=SUPL2+AT*SOUPL(IT,L)*BGPL2(IT,ITAU,K)
 364           CONTINUE
               SUPL2=SUPL2*EPTA
               SUMI2=0D0
               DO 366 IT=LTAU,MAXTAU
                  AT=1D0
                  IF(IT.eq.LTAU.or.IT.eq.MAXTAU) AT=5D-1
                  SUMI2=SUMI2+AT*SOUMI(IT,L)*BGMI2(LTAU,IT,K)
 366           CONTINUE
               SUMI2=SUMI2*EPTA
c emergent side
               FLPL=FLPL+ATAU*SUPL2
               FLMI=FLMI+ATAU*SUMI2
 360        CONTINUE
            FLPL=FLPL*DNTR/ETA
            FLMI=FLMI*DNTR/ETA
            DINUP(L)=FLPL
            DINDO(L)=FLMI
         ENDDO
      ENDDO

      RETURN
      END
c***********************************************************************
c for the escape probability (sphere) 
c for angle averaged radiation field 
c emergent up-down and internal angle averaged radiation  
      SUBROUTINE ESCINT(SOISTR,DISOTR,TAU,ISC)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION TAU
      INTEGER ISC

      DOUBLE PRECISION SOISTR(MAXTAU,MAXFRE),DISOTR(MAXFRE,MAXTAU)
      DOUBLE PRECISION ESCAPE(MAXFRE)

      DOUBLE PRECISION S00, SUESC, TAUX

      INTEGER IANG, ix, K

      DOUBLE PRECISION escsphere
      EXTERNAL escsphere

      save ESCAPE

      if(ISC.eq.1) then
         do ix=1,maxfre
            TAUX=S0(IX)*TAU
c escsphere computes 1-escape probability from a sphere
            ESCAPE(IX)=escsphere(TAUX)
         enddo
      endif                 

c for angle averaged radiation field 
c emergent radiation field up-down       
      DO IX=1,MAXFRE 
         S00=SOISTR(1,IX)
c internal radiation 
         DISOTR(IX,1)=S00*ESCAPE(IX)/S0(IX)
c escape radiation 
         SUESC=S00*(1d0-ESCAPE(IX)) 
         DO IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            if(ISC.eq.1) then 
               DINTPL(K)=SUESC
               SUIMI(K)=SUESC
            else
               SUIPL(K)=SUIPL(K)+SUESC
               SUIMI(K)=SUIPL(K)
            endif
         ENDDO
      ENDDO
      RETURN  
      END     
c***********************************************************************
c calculations of isotropic source function 
c for escape probability
      SUBROUTINE ESCSOUR(SOUPL,DISOTR,SOISTR,DIFMAX)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION SOISTR(MAXTAU,MAXFRE),DISOTR(MAXFRE,MAXTAU)
      DOUBLE PRECISION SOUPL(MAXTAU,LL),DIFMAX

      DOUBLE PRECISION DIF, SUMI, X, XXPM
      INTEGER IX, JX, L

      DO 445 IX=1,MAXFRE
         X=XEN(IX)
         SUMI=0D0
         DO 446 JX=1,MAXFRE
            SUMI=SUMI+FRDX(JX,IX)*DISOTR(JX,1)*A(JX)
 446     CONTINUE
         SOISTR(1,IX)=SUMI*X
 445  CONTINUE

      DO 450 IX=1,MAXFRE
         X=XEN(IX)
         XXPM=SOISTR(1,IX)
         L=IX*MAXANG
         SOUPL(1,L)=XXPM+SOUPL(1,L)
         IF(XXPM.gt.1D-30.and.X.lt.XEMAX) THEN
            DIF=XXPM/SOUPL(1,L)
            DIFMAX=DMAX1(DIF,DIFMAX)
         ENDIF
 450  CONTINUE

      RETURN
      END     
c***********************************************************************
c computes 1-p where 
c p=escape probability from a homogeneous sphere
c p(tau)=(3/4*tau){1-1/(2*tau^2)+[1/tau+1/(2*tau^2)]*exp(-2*tau)}
c for small tau
c p(tau)=1-3/4*tau+12*tau^2*\sum_(n=0) (-2*tau)^n * (n+4)/(n+5)!
c
      DOUBLE PRECISION FUNCTION escsphere(X)

      IMPLICIT NONE

      DOUBLE PRECISION x,x2,xr,xr2,accur,eps,del,dn,sum

      DATA accur/1d-14/,eps/3d-1/

      x2=x*2d0
      if(x.gt.eps) then
         xr2=0.5d0/x**2
         xr=1d0/x
         if(x.gt.40d0) then
            escsphere=1d0-0.75d0*xr*(1d0-xr2)
         else
            escsphere=1d0-0.75d0*xr*(1d0-xr2+(xr+xr2)*dexp(-x2))
         endif
      else
         sum=1d0/30d0
         dn=6d0
         del=1d0/120d0
         do while(dabs(del).gt.accur)
            del=-del*x2/dn
            sum=sum+del*(dn-1d0)
            dn=dn+1d0
         enddo
         escsphere=0.75d0*x-12d0*x**2*sum
      endif
      
      RETURN
      END                 
c***********************************************************************
c exact intensities for a spherical geometry
      SUBROUTINE EXASPH(SOUPL0,SOUMI0,DITAPL0,DITAMI0,
     &                  DIAVER,DIFLUX,TAUVEC,TAU,ISC)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION SOUPL0(MAXTAU,LL),SOUMI0(MAXTAU,LL), 
     1 DIAVER(MAXFRE,MAXTAU),DIFLUX(MAXFRE,MAXTAU), 
     2 DITAPL0(LL,MAXTAU),DITAMI0(LL,MAXTAU),TAUVEC(MAXTAU) 

      DOUBLE PRECISION TAU, ARX1, ARX2, ARY1, ARY2, AS, DZETA, DZETABS
      DOUBLE PRECISION EGMI, EGPL, ETA, F1, F2, F3, F4, GMI, GPL
      DOUBLE PRECISION R02, S0XT, SMAX, SMIN, SMITI, SOURINT, SPLTI
      DOUBLE PRECISION SPR, SS, SUMI1, SUPL1, TAUPR, TCUR, TT, UU

      INTEGER ISC, IANG, IS, ITAU, IX, JA, JTA, K, KA

      DOUBLE PRECISION CMIN
      DATA CMIN/1D-30/
        


      DO IX=1,MAXFRE 
         S0XT=S0(IX)*TAU
         DO IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            ETA=UANG(IANG)
            DO 260 ITAU=1,MAXTAU
               TCUR=TAUVEC(ITAU)
               R02=TCUR**2*(1D0-ETA**2)
               SS=TCUR*ETA
      
c positive angles    
               SMIN=-DSQRT(1D0-R02)

               SUPL1=0D0 
          
               DO 250 IS=1,NSPH
                  AS=AGSPH(IS) 
                  SPR=(SS-SMIN)*UGSPH(IS)+SMIN
                  TAUPR=DSQRT(R02+SPR**2)                          
                  DZETA=SPR/TAUPR
                  DZETABS=DABS(DZETA) 
                  CALL QHUNT(TAUVEC,MAXTAU,TAUPR,JTA)
                  if(JTA.eq.1.and.ISC.eq.1.and.INPUT.eq.0) JTA=JTA+1 
                  ARX1=TAUVEC(JTA)
                  ARX2=TAUVEC(JTA+1)

                  CALL DLOCATE(UANG,MAXANG,DZETABS,JA)
                  if(JA.eq.MAXANG) JA=MAXANG-1
                  if(JA.ne.0) then
                     ARY1=UANG(JA)
                     ARY2=UANG(JA+1)
                  else
                     ARY1=-UANG(1)
                     ARY2=UANG(1)
                  endif
                  KA=JA+(IX-1)*MAXANG
                  if(JA.ne.0) then
                     if(DZETA.ge.0d0) then
                        F1=DLOG(SOUPL0(JTA,KA)+CMIN)
                        F2=DLOG(SOUPL0(JTA+1,KA)+CMIN)
                        F3=DLOG(SOUPL0(JTA+1,KA+1)+CMIN)
                        F4=DLOG(SOUPL0(JTA,KA+1)+CMIN)
                     else
                        F1=DLOG(SOUMI0(JTA,KA)+CMIN)
                        F2=DLOG(SOUMI0(JTA+1,KA)+CMIN)
                        F3=DLOG(SOUMI0(JTA+1,KA+1)+CMIN)
                        F4=DLOG(SOUMI0(JTA,KA+1)+CMIN)
                     endif
                  else
c JA=0
                     F1=DLOG(SOUMI0(JTA,KA+1)+CMIN)
                     F2=DLOG(SOUMI0(JTA+1,KA+1)+CMIN)
                     F3=DLOG(SOUPL0(JTA+1,KA+1)+CMIN)
                     F4=DLOG(SOUPL0(JTA,KA+1)+CMIN)
                  endif
c interpolate to get source at an angle DZETA
                  TT=(TAUPR-ARX1)/(ARX2-ARX1)
                  if(JA.ne.0) then
                     UU=(DZETABS-ARY1)/(ARY2-ARY1)
                  else
                     UU=(DZETA-ARY1)/(ARY2-ARY1)
                  endif
                  SOURINT=(1D0-UU)*((1D0-TT)*F1+TT*F2)+
     &                 UU*(TT*F3+(1D0-TT)*F4)
                  SOURINT=DEXP(SOURINT)

                  GPL=(SS-SPR)*S0XT
                  EGPL=0D0
                  if(GPL.lt.7d1) EGPL=DEXP(-GPL)      
                  SUPL1=SUPL1+AS*SOURINT*EGPL
 250           CONTINUE 
               SUPL1=SUPL1*(SS-SMIN)*TAU

c negative angles
               SMAX=-SMIN
      
               SUMI1=0D0 
               if(ITAU.ne.MAXTAU) then
                  DO 255 IS=1,NSPH
                     AS=AGSPH(IS)                              
                     SPR=(SMAX-SS)*UGSPH(IS)+SS   
                     TAUPR=DSQRT(R02+SPR**2)                          
                     DZETA=SPR/TAUPR
                     DZETABS=DABS(DZETA) 
                     CALL QHUNT(TAUVEC,MAXTAU,TAUPR,JTA)        
                     if(JTA.eq.1.and.ISC.eq.1.and.INPUT.eq.0) JTA=JTA+1 
                     ARX1=TAUVEC(JTA)
                     ARX2=TAUVEC(JTA+1)
        
                     CALL DLOCATE(UANG,MAXANG,DZETA,JA)
                     if(JA.eq.MAXANG) JA=MAXANG-1
                     if(JA.ne.0) then
                        ARY1=UANG(JA)
                        ARY2=UANG(JA+1)
                     else
                        ARY1=-UANG(1)
                        ARY2=UANG(1)
                     endif
                     KA=JA+(IX-1)*MAXANG
                     if(JA.ne.0) then
                        F1=DLOG(SOUMI0(JTA,KA)+CMIN)
                        F2=DLOG(SOUMI0(JTA+1,KA)+CMIN)
                        F3=DLOG(SOUMI0(JTA+1,KA+1)+CMIN)
                        F4=DLOG(SOUMI0(JTA,KA+1)+CMIN)
                     else
c JA=0
                        F1=DLOG(SOUPL0(JTA,KA+1)+CMIN)
                        F2=DLOG(SOUPL0(JTA+1,KA+1)+CMIN)
                        F3=DLOG(SOUMI0(JTA+1,KA+1)+CMIN)
                        F4=DLOG(SOUMI0(JTA,KA+1)+CMIN)
                     endif
                     TT=(TAUPR-ARX1)/(ARX2-ARX1)
                     UU=(DZETA-ARY1)/(ARY2-ARY1)
                     SOURINT=(1D0-UU)*((1D0-TT)*F1+TT*F2)+
     &                    UU*(TT*F3+(1D0-TT)*F4)
                     SOURINT=DEXP(SOURINT)

                     GMI=(SPR-SS)*S0XT
                     EGMI=0D0
                     if(GMI.lt.7d1) EGMI=DEXP(-GMI)
                     SUMI1=SUMI1+AS*SOURINT*EGMI
 255              CONTINUE                  
                  SUMI1=SUMI1*(SMAX-SS)*TAU
               endif 
                                
               DITAPL0(K,ITAU)=SUPL1
               DITAMI0(K,ITAU)=SUMI1 
 260        CONTINUE 
            if(INPUT.eq.1.and.ISC.eq.1) then 
               DINTPL(K)=DITAPL0(K,MAXTAU)
            else
               SUIPL(K)=SUIPL(K)+DITAPL0(K,MAXTAU)
               SUIMI(K)=SUIPL(K)
            endif
            
         ENDDO
      ENDDO
c averaged intensities  
      IF(NISO.eq.1.and.ISC.eq.ISOTR) THEN
         DO ITAU=1,MAXTAU
            DO IX=1,MAXFRE
               SPLTI=0D0
               SMITI=0D0
               DO 419 IANG=1,MAXANG
                  K=IANG+(IX-1)*MAXANG
                  SPLTI=SPLTI+AANG(IANG)*DITAPL0(K,ITAU) 
                  SMITI=SMITI+AANG(IANG)*DITAMI0(K,ITAU)
 419           CONTINUE
               DIAVER(IX,ITAU)=DMAX1(CMIN,(SPLTI+SMITI)/2D0)  
               DIFLUX(IX,ITAU)=(SPLTI-SMITI)/2D0
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END
c***********************************************************************
c for spherical geometry 
c for angle averaged radiation field 
c emergent up-down and internal angle averaged radiation  
      SUBROUTINE ISOSPH(SOISPL,SOISMI,DIAVER,DIFLUX,TAUVEC,TAU,ISC)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION SOISPL(MAXTAU,MAXFRE),SOISMI(MAXTAU,MAXFRE), 
     1 DIAVER(MAXFRE,MAXTAU),DIFLUX(MAXFRE,MAXTAU),TAUVEC(MAXTAU) 

      DOUBLE PRECISION TAU, AS, EGMI, EGPL, ETA, GMI, GPL
      DOUBLE PRECISION R02, S0XT, SMAX, SMIN, SOURINT
      DOUBLE PRECISION SPR, SS, SUMI1, SUPL1, TAUPR, TCUR
      DOUBLE PRECISION SUMMI1, SUMPL1

      INTEGER ISC, IANG, IS, ITAU, IX, JTA, K

      DOUBLE PRECISION CMIN
      DATA CMIN/1D-30/

c for angle averaged radiation field 
c emergent radiation field up-down       
        
      DO IX=1,MAXFRE 
         S0XT=S0(IX)*TAU
         DO IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            ETA=UANG(IANG)     
            R02=1D0-ETA**2
            SS=ETA 
            SMIN=-ETA
            SUPL1=0D0 
         
            DO 150 IS=1,NSPH
               AS=AGSPH(IS) 
               SPR=(SS-SMIN)*UGSPH(IS)+SMIN
               TAUPR=DSQRT(R02+SPR**2)                          
               CALL QHUNT(TAUVEC,MAXTAU,TAUPR,JTA)
c interpolate in SOISPL or SOISMI to get SOURINT at a given TAUPR
               if(SPR.ge.0d0) then 
                  CALL DLINEAR(TAUVEC(JTA),DLOG(SOISPL(JTA,IX)+CMIN),
     *                 TAUVEC(JTA+1),DLOG(SOISPL(JTA+1,IX)+CMIN),TAUPR,
     +                 SOURINT)
               else
                  CALL DLINEAR(TAUVEC(JTA),DLOG(SOISMI(JTA,IX)+CMIN),
     *                 TAUVEC(JTA+1),DLOG(SOISMI(JTA+1,IX)+CMIN),TAUPR,
     +                 SOURINT)
               endif 
               GPL=(SS-SPR)*S0XT
               EGPL=0D0
               if(GPL.lt.7d1) EGPL=DEXP(-GPL)      
               SUPL1=SUPL1+AS*DEXP(SOURINT)*EGPL
 150        CONTINUE 
            SUPL1=SUPL1*(SS-SMIN)*TAU


            if(INPUT.eq.1.and.ISC.eq.1) then 
               DINTPL(K)=SUPL1
               SUIMI(K)=SUPL1
            else
               SUIPL(K)=SUIPL(K)+SUPL1
               SUIMI(K)=SUIPL(K)
            endif
c no emergent inside  SUIMI(K)=0d0
         ENDDO
      ENDDO
  
c internal angle averaged radiation field          
      DO IX=1,MAXFRE 
         S0XT=S0(IX)*TAU
         DO ITAU=1,MAXTAU
            TCUR=TAUVEC(ITAU)
      
c integral over angles      
            SUMPL1=0D0
            SUMMI1=0D0
            DO 260 IANG=1,MAXANG
               ETA=UANG(IANG)
               R02=TCUR**2*(1D0-ETA**2)
               SS=TCUR*ETA
      
c positive angles    
               SMIN=-DSQRT(1D0-R02)

               SUPL1=0D0 
          
               DO 250 IS=1,NSPH
                  AS=AGSPH(IS) 
                  SPR=(SS-SMIN)*UGSPH(IS)+SMIN
                  TAUPR=DSQRT(R02+SPR**2)                          
                  CALL QHUNT(TAUVEC,MAXTAU,TAUPR,JTA)
c interpolate in SOISTA to get SOURINT at a given TAUPR
                  if(SPR.ge.0d0) then 
                     CALL DLINEAR(TAUVEC(JTA),DLOG(SOISPL(JTA,IX)+CMIN),
     *                    TAUVEC(JTA+1),DLOG(SOISPL(JTA+1,IX)+CMIN),
     +                    TAUPR,SOURINT)
                  else
                     CALL DLINEAR(TAUVEC(JTA),DLOG(SOISMI(JTA,IX)+CMIN),
     *                    TAUVEC(JTA+1),DLOG(SOISMI(JTA+1,IX)+CMIN),
     +                    TAUPR,SOURINT)
                  endif 
                  GPL=(SS-SPR)*S0XT
                  EGPL=0D0
                  if(GPL.lt.7d1) EGPL=DEXP(-GPL)      
                  SUPL1=SUPL1+AS*DEXP(SOURINT)*EGPL
 250           CONTINUE 
               SUPL1=SUPL1*(SS-SMIN)

c negative angles
               SMAX=-SMIN      
               SUMI1=0D0 
               DO 255 IS=1,NSPH
                  AS=AGSPH(IS)                              
                  SPR=(SMAX-SS)*UGSPH(IS)+SS   
                  TAUPR=DSQRT(R02+SPR**2)                          
                  CALL QHUNT(TAUVEC,MAXTAU,TAUPR,JTA)        
c interpolate in SOISMI to get SOURINT at a given TAUPR
                  CALL DLINEAR(TAUVEC(JTA),DLOG(SOISMI(JTA,IX)+CMIN),
     *                 TAUVEC(JTA+1),DLOG(SOISMI(JTA+1,IX)+CMIN),
     +                 TAUPR,SOURINT)
                  GMI=(SPR-SS)*S0XT
                  EGMI=0D0
                  if(GMI.lt.7d1) EGMI=DEXP(-GMI)
                  SUMI1=SUMI1+AS*DEXP(SOURINT)*EGMI
 255           CONTINUE                  
               SUMI1=SUMI1*(SMAX-SS)
               SUMPL1=SUMPL1+AANG(IANG)*SUPL1
               SUMMI1=SUMMI1+AANG(IANG)*SUMI1
 260        CONTINUE 
            DIAVER(IX,ITAU)=DMAX1(CMIN,(SUMPL1+SUMMI1)*TAU/2D0)  
            DIFLUX(IX,ITAU)=(SUMPL1-SUMMI1)*TAU/2D0  
         ENDDO
      ENDDO
      RETURN  
      END     
c***********************************************************************     
c for angle averaged radiation field 
c emergent up-down and internal angle averaged radiation
      SUBROUTINE ISOINT(SOISPL,SOISMI,DIAVER,DIFLUX, 
     &      BGPL1,BGMI1,BGMIN,BGMINA,BGPLA,BGMIA,DELTAU,FILF1,ISC)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION BGPL1(MAXTAU,MAXTAU,KK),BGMI1(MAXTAU,MAXTAU,KK),
     1 BGMIN(MAXTAU,MAXTAU,KK),BGMINA(MAXTAU,MAXTAU,MAXFRE),   
     2 BGPLA(MAXTAU,MAXTAU,MAXFRE),BGMIA(MAXTAU,MAXTAU,MAXFRE),
     3 SOISPL(MAXTAU,MAXFRE),SOISMI(MAXTAU,MAXFRE), 
     4 DIAVER(MAXFRE,MAXTAU),DIFLUX(MAXFRE,MAXTAU)

      DOUBLE PRECISION DELTAU, FILF1, AT, EPTA, ETA, SUMI1, SUPL1, SUPL2

      INTEGER ISC, IANG, IT, ITAU, IX, K, L

      DOUBLE PRECISION CMIN
      DATA CMIN/1D-30/

c for angle averaged radiation field 
c emergent radiation field up-down
      DO 301 IX=1,MAXFRE
         DO 303 IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            L=K
            ETA=UANG(IANG)
            EPTA=DELTAU/ETA
            SUPL1=0D0
            SUPL2=0D0
            IF(IGEOM.ne.3)THEN 
               DO 305 IT=1,MAXTAU
                  AT=1D0
                  IF(IT.eq.1.or.IT.eq.MAXTAU) AT=5D-1
                  IF(IGROUND.eq.2) SUPL2=SUPL2+AT*SOISMI(IT,IX)
     +                                   *BGMIN(IT,MAXTAU,K)
                  SUPL1=SUPL1+AT*SOISPL(IT,IX)*BGPL1(IT,MAXTAU,K)
 305           CONTINUE
               SUPL1=(SUPL1+SUPL2*FILF1)*EPTA  
            ENDIF 
            SUMI1=0D0
            DO 307 IT=1,MAXTAU
               AT=1D0
               IF(IT.eq.1.or.IT.eq.MAXTAU) AT=5D-1
               SUMI1=SUMI1+AT*SOISMI(IT,IX)*BGMI1(1,IT,K)
 307        CONTINUE
            SUMI1=SUMI1*EPTA

            IF(IGEOM.ne.3) then 
               if(ISC.eq.1) then 
                  DINTPL(K)=DINTPL(K)+SUPL1
               else
                  SUIPL(L)=SUIPL(L)+SUPL1
               endif
            ENDIF 
            SUIMI(L)=SUIMI(L)+SUMI1
 303     CONTINUE
 301  CONTINUE
c internal angle averaged radiation field 
      DO 311 IX=1,MAXFRE
         DO 313 ITAU=1,MAXTAU
            SUPL1=0D0
            IF(IGROUND.eq.2) THEN
               DO 315 IT=1,MAXTAU
                  AT=1D0
                  IF(IT.eq.1.or.IT.eq.MAXTAU) AT=5D-1
                  SUPL1=SUPL1+AT*SOISMI(IT,IX)*BGMINA(IT,ITAU,IX)
 315           CONTINUE
               SUPL1=SUPL1*FILF1
            ENDIF
            IF(ITAU.ne.1) THEN
               DO 317 IT=1,ITAU
                  AT=1D0
                  IF(IT.eq.1.or.IT.eq.ITAU) AT=5D-1
                  SUPL1=SUPL1+AT*SOISPL(IT,IX)*BGPLA(IT,ITAU,IX)
 317           CONTINUE
            ENDIF
            SUMI1=0D0
            IF(ITAU.ne.MAXTAU) THEN 
               DO 319 IT=ITAU,MAXTAU
                  AT=1D0
                  IF(IT.eq.ITAU.or.IT.eq.MAXTAU) AT=5D-1
                  SUMI1=SUMI1+AT*SOISMI(IT,IX)*BGMIA(ITAU,IT,IX)
 319           CONTINUE
            ENDIF 
            DIAVER(IX,ITAU)=DMAX1(CMIN,(SUPL1+SUMI1)*DELTAU/2D0) 
            DIFLUX(IX,ITAU)=(SUPL1-SUMI1)*DELTAU/2D0 
 313     CONTINUE
 311  CONTINUE


      RETURN
      END
c***********************************************************************          
      SUBROUTINE EXAINT(SOUPL0,SOUMI0,DITAPL0,DITAMI0,DIAVER,DIFLUX,
     &      BGPL1,BGMI1,BGMIN,DELTAU,FILF1,ISC)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION BGPL1(MAXTAU,MAXTAU,KK),BGMI1(MAXTAU,MAXTAU,KK),
     1 BGMIN(MAXTAU,MAXTAU,KK),SOUPL0(MAXTAU,LL),SOUMI0(MAXTAU,LL), 
     2 DITAPL0(LL,MAXTAU),DITAMI0(LL,MAXTAU), 
     3 DIAVER(MAXFRE,MAXTAU),DIFLUX(MAXFRE,MAXTAU) 

      DOUBLE PRECISION DELTAU, FILF1, AT, EPTA, ETA, SUMI1, SUPL1
      DOUBLE PRECISION SMITI, SPLTI

      INTEGER ISC, IANG, IT, ITAU, IX, K, L, LTAU

      DOUBLE PRECISION CMIN
      DATA CMIN/1D-30/
        
      DO IX=1,MAXFRE 
         DO IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG 
            L=K 
            ETA=UANG(IANG)
            EPTA=DELTAU/ETA 
            DO 260 ITAU=1,MAXTAU
               LTAU=MAXTAU-ITAU+1 
               SUPL1=0D0 
               IF(IGROUND.eq.2) THEN
                  DO 251 IT=1,MAXTAU
                     AT=1D0
                     IF(IT.eq.1.or.IT.eq.MAXTAU) AT=5D-1
                     SUPL1=SUPL1+AT*SOUMI0(IT,L)*BGMIN(IT,ITAU,K)
 251              CONTINUE
                  SUPL1=SUPL1*FILF1 
               ENDIF
               IF(ITAU.ne.1) THEN 
                  DO 250 IT=1,ITAU 
                     AT=1D0 
                     IF(IT.eq.1.or.IT.eq.ITAU) AT=5D-1 
                     SUPL1=SUPL1+AT*SOUPL0(IT,L)*BGPL1(IT,ITAU,K)
 250              CONTINUE 
               ENDIF 
               SUPL1=SUPL1*EPTA 
               SUMI1=0D0 
               IF(LTAU.ne.MAXTAU) THEN 
                  DO 240 IT=LTAU,MAXTAU 
                     AT=1D0 
                     IF(IT.eq.LTAU.or.IT.eq.MAXTAU) AT=5D-1 
                     SUMI1=SUMI1+AT*SOUMI0(IT,L)*BGMI1(LTAU,IT,K)
 240              CONTINUE 
                  SUMI1=SUMI1*EPTA 
               ENDIF
               DITAPL0(L,ITAU)=SUPL1
               DITAMI0(L,LTAU)=SUMI1 
 260        CONTINUE 
c emergent radiation field up-down 


            IF(IGEOM.ne.3) then 
               if(ISC.eq.1) then
                  DINTPL(K)=DINTPL(K)+DITAPL0(L,MAXTAU)
               else
                  SUIPL(L)=SUIPL(L)+DITAPL0(L,MAXTAU)
               endif
            ENDIF 
            SUIMI(L)=SUIMI(L)+DITAMI0(L,1)
         ENDDO
      ENDDO

      IF(NISO.eq.1.and.ISC.eq.ISOTR) THEN
         DO ITAU=1,MAXTAU
            DO IX=1,MAXFRE
               SPLTI=0D0
               SMITI=0D0
               DO 419 IANG=1,MAXANG
                  K=IANG+(IX-1)*MAXANG
                  SPLTI=SPLTI+AANG(IANG)*DITAPL0(K,ITAU) 
                  SMITI=SMITI+AANG(IANG)*DITAMI0(K,ITAU)
 419           CONTINUE
               DIAVER(IX,ITAU)=DMAX1(CMIN,(SPLTI+SMITI)/2D0)  
               DIFLUX(IX,ITAU)=(SPLTI-SMITI)/2D0
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END

c***********************************************************************
c calculations of isotropic source function 
      SUBROUTINE ISOSOUR(SOUPL,SOUMI,DIAVER,DIFLUX,SOISPL,SOISMI,DIFMAX)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION SOUPL(MAXTAU,LL),SOUMI(MAXTAU,LL), 
     1 SOISPL(MAXTAU,MAXFRE),SOISMI(MAXTAU,MAXFRE),
     2 DIAVER(MAXFRE,MAXTAU),DIFLUX(MAXFRE,MAXTAU)  

      DOUBLE PRECISION DIFMAX, DIF, SUM1, SUM2, X, XXPM

      INTEGER IANG, ITAU, IX, JX, L



      DO ITAU=1,MAXTAU
         DO IX=1,MAXFRE
            X=XEN(IX)
            SUM1=0D0
            SUM2=0D0
            DO 445 JX=1,MAXFRE
               SUM1=SUM1+FRDX(JX,IX)*DIAVER(JX,ITAU)*A(JX)
               SUM2=SUM2+FRDS(JX,IX)*DIFLUX(JX,ITAU)*A(JX)
 445        CONTINUE
c     SOISTA(ITAU,IX)=SUMI*X/2D0 
            SOISPL(ITAU,IX)=(SUM1+SUM2)*X
            SOISMI(ITAU,IX)=(SUM1-SUM2)*X
         ENDDO
      ENDDO

      DO IX=1,MAXFRE
         X=XEN(IX) 
         DO IANG=1,MAXANG
            L=IANG+(IX-1)*MAXANG 
            DO ITAU=1,MAXTAU
               XXPM=SOISPL(ITAU,IX)
               SOUPL(ITAU,L)=XXPM+SOUPL(ITAU,L)
               SOUMI(ITAU,L)=XXPM+SOUMI(ITAU,L)
               IF(ITAU.eq.MAXTAU.and.XXPM.gt.1D-30.and.X.lt.XEMAX) THEN
                  DIF=XXPM/SOUPL(MAXTAU,L)
                  DIFMAX=DMAX1(DIF,DIFMAX)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END     
c***********************************************************************
c calculations of exact source function 
      SUBROUTINE EXASOUR(SOUPL,SOUMI,DITAPL0,DITAMI0,SOUPL0,SOUMI0,
     &                   DIFMAX) 

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION SOUPL(MAXTAU,LL),SOUMI(MAXTAU,LL), 
     1 SOUPL0(MAXTAU,LL),SOUMI0(MAXTAU,LL), 
     2 DITAPL0(LL,MAXTAU),DITAMI0(LL,MAXTAU)    

      DOUBLE PRECISION DIFMAX, DIACMI, DIACPL, DIF, FRDPMEL, FRDPPEL
      DOUBLE PRECISION X, XIMI, XIPL, XXMI, XXPL

      INTEGER IANG, ITAU, IX, JANG, JX, L, L1


      DO IX=1,MAXFRE 
         X=XEN(IX) 
         DO IANG=1,MAXANG
            L=IANG+(IX-1)*MAXANG 
            DO 185 ITAU=1,MAXTAU
               XIPL=0D0 
               XIMI=0D0 
               DO JX=1,MAXFRE 
                  DO JANG=1,MAXANG
                     L1=JANG+(JX-1)*MAXANG 
                     FRDPPEL=FRDPP(L1,L)
                     FRDPMEL=FRDPM(L1,L)
                     DIACPL=DITAPL0(L1,ITAU)*AC(L1)
                     DIACMI=DITAMI0(L1,ITAU)*AC(L1)
                     XIPL=XIPL+FRDPPEL*DIACPL+FRDPMEL*DIACMI
                     XIMI=XIMI+FRDPMEL*DIACPL+FRDPPEL*DIACMI 
                  ENDDO
               ENDDO
               XXPL=XIPL*X 
               XXMI=XIMI*X  
               SOUPL0(ITAU,L)=XXPL
               SOUMI0(ITAU,L)=XXMI
               SOUPL(ITAU,L)=XXPL+SOUPL(ITAU,L) 
               SOUMI(ITAU,L)=XXMI+SOUMI(ITAU,L)  
 185        CONTINUE 
            IF(XXPL.gt.1D-30.and.X.lt.XEMAX) THEN
               DIF=XXPL/SOUPL(MAXTAU,L)
               DIFMAX=DMAX1(DIF,DIFMAX)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END     
c***********************************************************************
c geometrical factors depending of IGEOM 
      SUBROUTINE GEOMFAC(TAUVEC,T0TAUG,CINP,CINM,CEMP,CEMM)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION CINP(MAXTAU,MAXTAU,MAXANG),
     1            CINM(MAXTAU,MAXTAU,MAXANG),CEMP(MAXTAU,MAXTAU,MAXANG),
     1            CEMM(MAXTAU,MAXTAU,MAXANG),TAUVEC(MAXTAU)

      DOUBLE PRECISION T0TAUG, CEMMI, CEMPL, CINMI, CINPL, CPL, DQRR
      DOUBLE PRECISION DTI, DTL, ETA, ETADTI, PHI0, PHIDU, PHIST, RI
      DOUBLE PRECISION RI2, RL, RL2, SIET, SIETA, SINPHI0, SINSTI
      DOUBLE PRECISION SINSTL, TSI, TT0, TT0RI2, TT0RL2, TT2
      DOUBLE PRECISION TTAG, T0TA2

      INTEGER IANG, LTAU, ITAU, ITAU1

c PI2=2/pi    
      DOUBLE PRECISION PI, PI2
      DATA PI/3.141592653589792384D0/,PI2/0.63661977236D0/ 


      T0TA2=T0TAUG/2D0
      DO 9 IANG=1,MAXANG
         ETA=UANG(IANG)
         SIET=DSQRT(1D0-ETA*ETA)
         SIETA=SIET/ETA
         TSI=PI2*T0TAUG*SIET
         DO 7 ITAU=2,MAXTAU
            ITAU1=ITAU-1
            DTI=TAUVEC(ITAU)
            ETADTI=ETA*DTI
            DO 5 LTAU=1,ITAU1
               DTL=TAUVEC(LTAU)
               CINPL=0D0
               CINMI=0D0
               CEMPL=0D0
               CEMMI=0D0
               TT0=DABS(DTL-DTI)*SIETA
c slab and cylinder      
               IF(IGEOM.eq.1 .or. IGEOM.eq.2) THEN
                  TTAG=TT0*T0TA2
                  IF(TTAG .lt. 1D0) THEN
                     CEMPL=DSQRT(1D0-TTAG**2)
                     CINPL=(DACOS(TTAG)-TTAG*CEMPL)*PI2
                     CEMPL=CEMPL*TSI
                     CINMI=CINPL
                     CEMMI=CEMPL
                  ENDIF
               ENDIF
c end of slab and cylinder 

c hemisphere
               IF(IGEOM.eq.3) THEN
                  TT2=TT0*TT0
                  IF(ITAU.eq.MAXTAU) THEN
                     RI=0D0
                     RI2=0D0
                  ELSE
                     RI2=1D0-DTI*DTI
                     RI=DSQRT(RI2)
                  ENDIF
                  RL2=1D0-DTL*DTL
                  RL=DSQRT(RL2)
                  IF(TT0.le.(RL-RI)) THEN
                     CINPL=1D0
                     CINMI=RI2/RL2
                     CEMPL=2D0*ETADTI
                     CEMMI=0D0
                  ENDIF
                  IF(TT0.gt.(RL-RI).and.TT0.lt.(RI+RL)) THEN
                     TT0RI2=2D0*TT0*RI
                     TT0RL2=2D0*TT0*RL
                     PHIST=DACOS((RI2-RL2+TT2)/TT0RI2)
                     PHIDU=DACOS((RL2-RI2+TT2)/TT0RL2)
                     DQRR=DSQRT(4D0*RI2*RL2-(TT2-RI2-RL2)**2)
                     CPL=(RI2*PHIST+RL2*PHIDU-0.5D0*DQRR)/PI
                     CINPL=CPL/RI2
                     CINMI=CPL/RL2
                     SINSTI=DQRR/TT0RI2
c      IF(ETA.ge.ETATANI) THEN
                     PHI0=PHIST
                     SINPHI0=SINSTI
c      ELSE
c      PHI0=DMIN1(PHIST,PI/2D0)
c      SINPHI0=DMIN1(1D0,SINSTI)
c      ENDIF
                     CEMPL=PI2*(ETADTI*PHI0+SIET*RI*SINPHI0)
                     SINSTL=DQRR/TT0RL2
                     CEMMI=PI2*(-ETA*DTL*PHIDU+SIET*RL*SINSTL)
                  ENDIF
               ENDIF
c     end of hemisphere

               CINP(LTAU,ITAU,IANG)=CINPL
               CINM(LTAU,ITAU,IANG)=CINMI
               CEMP(LTAU,ITAU,IANG)=CEMPL
               CEMM(LTAU,ITAU,IANG)=CEMMI
 5          CONTINUE
 7       CONTINUE
    9 CONTINUE

      RETURN
      END
c***********************************************************************
c geometrical factors depending of IGEOM for DISK with HOLE geometry
      SUBROUTINE GEOMGR(TAUVEC,T0TAUG,CGRIN,CGREX) 

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION CGRIN(MAXTAU,MAXTAU,MAXANG),
     2           CGREX(MAXTAU,MAXTAU,MAXANG),TAUVEC(MAXTAU)

      DOUBLE PRECISION T0TAUG, CEMPL, CINPL, CPL, DQRR, DTI
      DOUBLE PRECISION DTL, ETA, ETADTI, PHIDU, PHIST, RI, RI2, RL
      DOUBLE PRECISION RL2, SIET, SIETA, SINSTI, T0TA2, TT0, TT0RI2
      DOUBLE PRECISION TT0RL2, TT2, TTAG

      INTEGER IANG, ITAU, LTAU

c PI2=2/pi    
      DOUBLE PRECISION PI, PI2
      DATA PI/3.141592653589792384D0/,PI2/0.63661977236D0/ 

      T0TA2=T0TAUG/2D0       
      DO 109 IANG=1,MAXANG
         ETA=UANG(IANG)
         SIET=DSQRT(1D0-ETA*ETA)
         SIETA=SIET/ETA
         DO 107 ITAU=1,MAXTAU
            DTI=TAUVEC(ITAU)
            ETADTI=ETA*DTI
            DO 105 LTAU=1,MAXTAU
               DTL=TAUVEC(LTAU)
               CINPL=0D0
               CEMPL=0D0
               TT0=(DTL+DTI)*SIETA

c SLAB AND CYLINDER
               IF(IGEOM.eq.1 .or. IGEOM.eq.2) THEN
                  TTAG=TT0*T0TA2
                  IF(TTAG .lt. 1D0) THEN
                     CEMPL=DSQRT(1D0-TTAG**2)
                     CINPL=(DACOS(TTAG)-TTAG*CEMPL)*PI2
                     CEMPL=CEMPL*PI2*T0TAUG*SIET
                  ENDIF
               ENDIF
c end of SLAB AND CYLINDER

c HEMISPHERE
               IF(IGEOM.eq.3) THEN
                  TT2=TT0*TT0
                  IF(ITAU.eq.MAXTAU) THEN
                     RI2=0D0
                  ELSE
                     RI2=1D0-DTI*DTI
                  ENDIF
                  IF(LTAU.eq.MAXTAU) THEN
                     RL2=0D0
                  ELSE
                     RL2=1D0-DTL*DTL
                  ENDIF
                  RI=DSQRT(RI2)
                  RL=DSQRT(RL2)
      
                  IF(TT0.le.DABS(RL-RI)) THEN
                     IF(RL.ge.RI) THEN
                        CINPL=1D0
                        CEMPL=2D0*ETADTI
                     ELSE
                        CINPL=RL2/RI2
                        CEMPL=0D0
                     ENDIF
                  ENDIF
                  IF(TT0.gt.DABS(RL-RI).and.TT0.lt.(RI+RL)) THEN
                     TT0RI2=2D0*TT0*RI
                     TT0RL2=2D0*TT0*RL
                     PHIST=DACOS((RI2-RL2+TT2)/TT0RI2)
                     PHIDU=DACOS((RL2-RI2+TT2)/TT0RL2)
                     DQRR=DSQRT(4D0*RI2*RL2-(TT2-RI2-RL2)**2)
                     CPL=(RI2*PHIST+RL2*PHIDU-0.5D0*DQRR)/PI
                     SINSTI=DQRR/TT0RI2
                     CINPL=CPL/RI2
                     CEMPL=PI2*(ETADTI*PHIST+SIET*RI*SINSTI)
                  ENDIF
               ENDIF
c end of HEMISPHERE      

               CGRIN(LTAU,ITAU,IANG)=CINPL
               CGREX(LTAU,ITAU,IANG)=CEMPL
 105        CONTINUE
 107     CONTINUE
 109  CONTINUE

            
      RETURN
      END
c***********************************************************************      
      SUBROUTINE ENERFAC(TAUVEC,TAU,T0TAUG,CINP,CINM,CEMP,CEMM,
     &      BGPL1,BGMI1,BGPL2,BGMI2,BGPLA,BGMIA)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION CINP(MAXTAU,MAXTAU,MAXANG),
     1 CINM(MAXTAU,MAXTAU,MAXANG),CEMP(MAXTAU,MAXTAU,MAXANG),
     2 CEMM(MAXTAU,MAXTAU,MAXANG),BGPL1(MAXTAU,MAXTAU,KK),
     3 BGMI1(MAXTAU,MAXTAU,KK),BGPL2(MAXTAU,MAXTAU,KK),
     4 BGMI2(MAXTAU,MAXTAU,KK),BGPLA(MAXTAU,MAXTAU,MAXFRE),
     5 BGMIA(MAXTAU,MAXTAU,MAXFRE),TAUVEC(MAXTAU)

      DOUBLE PRECISION TAU, T0TAUG, AZ, DTI, DTL, EGPL, ETA, ETADTI
      DOUBLE PRECISION S0IX, SIET, SUM1, SUM2, SUPL, TAUETA, TDDE, TSI

      INTEGER IANG, ITAU, ITAU1, IX, K, LTAU


c PI2=2/pi    
      DOUBLE PRECISION PI2
      DATA PI2/0.63661977236D0/ 
           
      DO IX=1,MAXFRE
         DO IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            ETA=UANG(IANG)
            SIET=DSQRT(1D0-ETA*ETA)
            TSI=PI2*T0TAUG*SIET
            DO ITAU=1,MAXTAU
               BGPL1(ITAU,ITAU,K)=1D0
               BGMI1(ITAU,ITAU,K)=1D0
               BGPL2(ITAU,ITAU,K)=0D0
               BGMI2(ITAU,ITAU,K)=0D0
            ENDDO
         ENDDO
      ENDDO

      DO IX=1,MAXFRE
         S0IX=S0(IX)
         DO IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            ETA=UANG(IANG)
            TAUETA=TAU/ETA
            DO ITAU=2,MAXTAU
               ITAU1=ITAU-1
               DTI=TAUVEC(ITAU)
               ETADTI=ETA*DTI
               DO LTAU=1,ITAU1
                  DTL=TAUVEC(LTAU)
                  TDDE=TAUETA*(DTI-DTL)
                  EGPL=0D0
                  SUPL=S0IX*TDDE
                  IF(SUPL.lt.7d1) EGPL=DEXP(-SUPL)
c for radiation field inside
                  BGPL1(LTAU,ITAU,K)=EGPL*CINP(LTAU,ITAU,IANG)
                  BGMI1(LTAU,ITAU,K)=EGPL*CINM(LTAU,ITAU,IANG)
c     for the emergent radiation field
                  BGPL2(LTAU,ITAU,K)=EGPL*CEMP(LTAU,ITAU,IANG)
                  BGMI2(LTAU,ITAU,K)=EGPL*CEMM(LTAU,ITAU,IANG)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c angle averaged
      DO ITAU=1,MAXTAU
         DO LTAU=1,ITAU
            DO IX=1,MAXFRE
               SUM1=0D0
               SUM2=0D0
               DO 23 IANG=1,MAXANG
                  AZ=AANG(IANG)/UANG(IANG)
                  K=IANG+(IX-1)*MAXANG
                  SUM1=SUM1+AZ*BGPL1(LTAU,ITAU,K)
                  SUM2=SUM2+AZ*BGMI1(LTAU,ITAU,K)
 23            CONTINUE
               BGPLA(LTAU,ITAU,IX)=SUM1
               BGMIA(LTAU,ITAU,IX)=SUM2
            ENDDO
         ENDDO
      ENDDO
            
      RETURN
      END
c***********************************************************************
c FOR THE DISK WITH HOLES
      SUBROUTINE ENERGR(TAUVEC,TAU,CGRIN,CGREX,
     &      BGMIN,BGMEX,BGMINA)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION CGRIN(MAXTAU,MAXTAU,MAXANG),
     2 CGREX(MAXTAU,MAXTAU,MAXANG),BGMIN(MAXTAU,MAXTAU,*),
     3 BGMEX(MAXTAU,MAXTAU,*),BGMINA(MAXTAU,MAXTAU,*),TAUVEC(MAXTAU)

      DOUBLE PRECISION TAU, DTI, DTL, EGPL, ETA, CINPL
      DOUBLE PRECISION S0IX, SUM1, SUPL, TAUETA, TDDE

      INTEGER IANG, ITAU, IX, K, LTAU


      
      DO 119 IX=1,MAXFRE
         S0IX=S0(IX) 
         DO 117 IANG=1,MAXANG
            K=IANG+(IX-1)*MAXANG
            ETA=UANG(IANG)
            TAUETA=TAU/ETA
            DO 115 ITAU=1,MAXTAU
               DTI=TAUVEC(ITAU)
               DO 113 LTAU=1,MAXTAU
                  DTL=TAUVEC(LTAU)
                  TDDE=TAUETA*(DTI+DTL) 
                  CINPL=0D0
                  SUPL=S0IX*TDDE 
                  EGPL=0D0
                  IF(SUPL.lt.7d1) EGPL=DEXP(-SUPL)
c for radiation field inside
                  BGMIN(LTAU,ITAU,K)=EGPL*CGRIN(LTAU,ITAU,IANG)
c     for the emergent radiation field
                  BGMEX(LTAU,ITAU,K)=EGPL*CGREX(LTAU,ITAU,IANG)
 113           CONTINUE
 115        CONTINUE
 117     CONTINUE
 119  CONTINUE
c angle averaged 
      DO ITAU=1,MAXTAU
         DO LTAU=1,MAXTAU
            DO IX=1,MAXFRE 
               SUM1=0D0
               DO 123 IANG=1,MAXANG 
                  ETA=UANG(IANG)
                  K=IANG+(IX-1)*MAXANG
                  SUM1=SUM1+AANG(IANG)/ETA*BGMIN(LTAU,ITAU,K) 
 123           CONTINUE 
               BGMINA(LTAU,ITAU,IX)=SUM1 
            ENDDO
         ENDDO
      ENDDO
     
      RETURN
      END
c***********************************************************************
c        REDISTRIBUTION MATRICES FOR COMPTON SCATTERING 
c***********************************************************************
      SUBROUTINE CSRF  
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION COMPHI11(MAXFRE),REDXX(NMU),DMUA(NMU),F2(NMU)   

      DOUBLE PRECISION CFI, COKL, CONS, DMUPM, DMUPP, ETA, ETA1
      DOUBLE PRECISION ETA12, ETA2, EXXJS, H3, DLN, PHI, PHKL, QMU, SI
      DOUBLE PRECISION SI1, SI12, SI2, SIKL, SM, SMU, SP, SUM, SUM1
      DOUBLE PRECISION SUM11M, SUM11P, SUM2, X, X1, XXR

      INTEGER I, I1, IAU, ID, IX, J, J1, JSX, JX, JXMIN, K, K1, KE
      INTEGER KE1, KETA, L, L1, LETA, NX

      DOUBLE PRECISION PI, EPS
      DATA PI/3.141592653589792384D0/,EPS/1D-14/ 

      NX=MAXFRE 
c-------------------------------------------------------------------------
C CALCULATIONS OF ANGLE DEPENDENT REDISTRIBUTION FUNCTIONS EXPLICITELY FOR
c (CSM+RT) CODE. 
c-------------------------------------------------------------------------
      DLN=DLOG(1D1) 
      H3=PI*2D0 
      DO 11 IAU=1,NMU
         if(IAU.le.MAXANG) DMUA(IAU)=-UANG(MAXANG-IAU+1) 
         if(IAU.gt.MAXANG) DMUA(IAU)=UANG(IAU-MAXANG) 
 11   CONTINUE 
c------------------------------------------------------------------
c begin  STAGE I 
C mainly integrals over azimuth angle:
      DO L1=1,LL
         DO L=1,LL
            FRDPP(L,L1)=0.D0
            FRDPM(L,L1)=0.D0
         ENDDO
      ENDDO
      JXMIN=1
      DO 400 IX=1,NX
         X=XEN(IX)
         IF(IDISTF.eq.0) JXMIN=IX
         DO 300 JX=JXMIN,NX
            X1=XEN(JX)
c here computation of R(x,x1,mu) 
            DO 350 IAU=1,NMU 
               QMU=1D0-DMUA(IAU) 
               IF(QMU.gt.EPS) THEN 
                  CALL REDFUNC(X,X1,QMU,SMU) 
               ELSE 
                  SMU=0D0
               ENDIF 
               REDXX(IAU)=SMU  
               FRDMU(JX,IX,IAU)=SMU 
 350        CONTINUE 
            CALL QSPLINE(DMUA,REDXX,NMU,1d30,1d30,F2)

            DO 200 KETA=1,MAXANG  
               K=KETA+(IX-1)*MAXANG   
               ETA=UANG(KETA) 
               ETA2=ETA*ETA
               SI2=1D0-ETA2
               SI=DSQRT(SI2) 
               DO 100 LETA=1,KETA
c              DO 100 LETA=1,MAXANG 
                  K1=LETA+(JX-1)*MAXANG   
                  ETA1=UANG(LETA)
                  ETA12=ETA1*ETA1
                  SI12=1D0-ETA12
                  SI1=DSQRT(SI12)  
                  SIKL=SI*SI1
                  COKL=ETA*ETA1
                  SUM11P=0D0 
                  SUM11M=0D0
c BEGIN DELTA (azimuth) INTEGRATION LOOP
                  DO 50 ID=1,NS
                     PHI=UGS(ID)*PI 
                     CFI=DCOS(PHI)
                     PHKL=CFI*SIKL
c BEGIN CALN. OF  S+- (+mu,-mu') = S-+ (-mu,+mu') BLOCK MATRIX:
                     DMUPM=-COKL+PHKL
                     CALL QSPLINT(DMUA,REDXX,F2,NMU,DMUPM,SM)
c END CALN. OF  S+- (+mu,-mu') = S-+ (-mu,+mu') BLOCK MATRIX.
c BEGIN CALN. OF  S++ (+mu,+mu') = S-- (-mu,-mu') BLOCK MATRIX:
                     DMUPP=COKL+PHKL
                     CALL QSPLINT(DMUA,REDXX,F2,NMU,DMUPP,SP)
c END CALN. OF  S++ (+mu,+mu') = S-- (-mu,-mu') BLOCK MATRIX.
                     SUM11P=SUM11P+SP*AGS(ID)
                     SUM11M=SUM11M+SM*AGS(ID)
 50               CONTINUE              
c END azimuth INTEGRATION LOOP 
c 
c--------------------------------------------------------------------------
                  FRDPP(K1,K)=DMAX1(H3*SUM11P,0D0)
                  FRDPM(K1,K)=DMAX1(H3*SUM11M,0D0)
c CONSTRUCTION OF A PART OF THE REDISTRIBUTION MATRICES FRD's IS OVER.
 100           CONTINUE 
 200        CONTINUE 
 300     CONTINUE 
 400  CONTINUE 
c end of STAGE I
c---------------
c------------------------------------------------------------------
c begin  STAGE II 
c CONSTRUCTION OF FULL STRUCTURE OF FRD's IS DONE BELOW BY EXPLOITING THE
c SYMMETRY PROPERTIES OF THE CSM.
c ATTENTION! DO NOT CHANGE THE ORDER of LOOPS !!! 
c power-law electrons 
      IF(IDISTF.NE.0) THEN
         DO 180 I=1,NX
            X=XEN(I)
            DO 170 J=1,MAXANG
               K=J+(I-1)*MAXANG
               DO 160 I1=1,NX
                  X1=XEN(I1)
                  XXR=X/X1
                  KE1=J+(I1-1)*MAXANG
                  DO 150 J1=1,MAXANG
                     KE=J1+(I-1)*MAXANG
                     K1=J1+(I1-1)*MAXANG
                     IF(J1.gt.J) THEN
                        FRDPP(K1,K)=FRDPP(KE1,KE)*XXR
                        FRDPM(K1,K)=FRDPM(KE1,KE)*XXR
                     ELSE
                        FRDPP(K1,K)=FRDPP(K1,K)*XXR
                        FRDPM(K1,K)=FRDPM(K1,K)*XXR
                     ENDIF
 150              CONTINUE
c     END OF mu' LOOP
 160           CONTINUE
c     END OF X' LOOP
 170        CONTINUE
c     END OF mu LOOP
 180     CONTINUE
c     END OF X LOOP
         do I=1,NX
            do I1=1,NX
               do IAU=1,NMU 
                  FRDMU(I1,I,IAU)=FRDMU(I1,I,IAU)*XEN(I)/XEN(I1)
               enddo
            enddo
         enddo
c--------------------------------------------------------------------------
      ELSE 
c Maxwellian electrons 
c This line included to suppress a compiler warning
         EXXJS=0.0
         DO 280 I1=1,NX
            X1=XEN(I1)
            DO 270 I=1,NX
               X=XEN(I)
               XXR=X/X1 
               JSX=I1+I*(I-1)/2
               IF(I.GT.I1)  EXXJS=EXXY(JSX)*XXR 
               DO 260 J=1,MAXANG
                  K=J+(I-1)*MAXANG
                  KE1=J+(I1-1)*MAXANG
                  DO 250 J1=1,MAXANG
                     KE=J1+(I-1)*MAXANG
                     K1=J1+(I1-1)*MAXANG
                     IF(I.LE.I1) THEN
                        IF(J1.gt.J) THEN
                           FRDPP(K1,K)=FRDPP(KE1,KE)*XXR 
                           FRDPM(K1,K)=FRDPM(KE1,KE)*XXR 
                        ELSE
                           FRDPP(K1,K)=FRDPP(K1,K)*XXR
                           FRDPM(K1,K)=FRDPM(K1,K)*XXR
                        ENDIF
                     ELSE 
                        IF(J1.lt.J) THEN
                           FRDPP(K1,K)=FRDPP(KE,KE1)*EXXJS
                           FRDPM(K1,K)=FRDPM(KE,KE1)*EXXJS
                        ELSE 
                           FRDPP(K1,K)=FRDPP(K,K1)*EXXJS
                           FRDPM(K1,K)=FRDPM(K,K1)*EXXJS
                        ENDIF 
                     ENDIF 
c--------------------------------------------------------------------------
c--------------------------------------------------------------------------
 250              CONTINUE
c     END OF mu' LOOP
 260           CONTINUE
c     END OF mu LOOP
 270        CONTINUE
c     END OF X LOOP
 280     CONTINUE
c     END OF X' LOOP
         DO I1=1,NX
            X1=XEN(I1)
            DO I=1,NX
               X=XEN(I)
               XXR=X/X1
               JSX=I1+I*(I-1)/2
               IF(I.GT.I1) THEN
                  EXXJS=EXXY(JSX)*XXR
                  DO IAU=1,NMU 
                     FRDMU(I1,I,IAU)=FRDMU(I,I1,IAU)*EXXJS 
                  ENDDO 
               ELSE
                  DO IAU=1,NMU 
                     FRDMU(I1,I,IAU)=FRDMU(I1,I,IAU)*XXR 
                  ENDDO 
               ENDIF
            ENDDO 
         ENDDO 

      ENDIF 
c NOW WE SPECIFICALLY USE THE SYMMETRY RELATIONS: FRDMP=FRDPM; FRDMM=FRDPP:
c end of STAGE II
c begin  STAGE III
c CHECKING THE NORMALIZATION OF THE ANGLE DEPENDENT GENERAL  CSM 
c REDISTRIBUTION MATRIX:
c FIRST DEVELOP THE ANGLE AVERAGED REDISTRIBUTION MATRICES FOR WORK 
      DO I=1,NX
         DO I1=1,NX
c     sum=0d0 
            sum1=0d0 
            sum2=0d0 
            DO J=1,MAXANG
               K=J+(I-1)*MAXANG
               DO J1=1,MAXANG
                  K1=J1+(I1-1)*MAXANG
                  sum1=sum1+ CINT(K)*FRDPP(K1,K)*CINT(K1)
                  sum2=sum2+ CINT(K)*FRDPM(K1,K)*CINT(K1)
               ENDDO
            ENDDO
            FRDX(I1,I)=sum1+sum2  
            FRDS(I1,I)=sum1-sum2 
         ENDDO
      ENDDO
c     FRDX GIVEN ABOVE IS AN ANGLE AVERAGED REDISTRIBUTION FUNCTION. 
      DO 840 I=1,NX
         X=XEN(I)
         SUM=0.D0
         DO 841 I1=1,NX
            SUM=SUM + FRDX(I,I1)*XEN(I1)*A(I1)
 841     CONTINUE
         COMPHI11(I)=SUM 
 840  CONTINUE
c-----------------------------------------------------------------------------
c COMPHI(NX) GIVEN ABOVE IS <computed CS profile function>. IN PRINCIPLE,
c IT SHOULD BE SAME AS  S0(NX)  EVALUATED DIRECTLY IN THE CSM PART OF THE CODE
c (KN CROSS SECTION BASICALLY). BUT IN PRACTICE IT IS <not>, DUE TO THE FINITE
c QUADRATURE APPROXIMATIONS IN PERFORMING THE NUMERICAL INTEGRATIONS. HENCE 
c THERE IS A NEED TO RENORMALIZE THE REDISTRIBUTION FUNCTION.
c-----------------------------------------------------------------------------
      DO I=1,NX 
         DO I1=1,NX
            CONS=S0(I1)/COMPHI11(I1)
c NOW THE RENORMALIZATION OF THE ANGLE AVERAGED RF 
            FRDX(I1,I)=FRDX(I1,I) * CONS
            FRDS(I1,I)=FRDS(I1,I) * CONS
         ENDDO
      ENDDO
      DO I=1,NX 
         DO J=1,MAXANG 
            K=J+(I-1)*MAXANG
            DO I1=1,NX
               CONS=S0(I1)/COMPHI11(I1)
               DO J1=1,MAXANG
                  K1=J1+(I1-1)*MAXANG
c NOW THE RENORMALIZATION OF THE FULL ANGLE DEPENDENT  <CSM>
c REDISTRIBUTION MATRICES  FRDPP,FRDPM 
                  FRDPP(K1,K)=FRDPP(K1,K) * CONS
                  FRDPM(K1,K)=FRDPM(K1,K) * CONS
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c END  OF THE OF CALCULATION OF NORMALIZED CSM REDISTRIBUTION MATRICES 
c [vector S_{ij} (++)] = [vector S_{ij} (--)]   AND
c [vector S_{ij} (+-)] = [vector S_{ij} (-+)],  i,j = 1,2.
c WE ARE NOW READY TO EMPLOY THESE REDISTRIBUTION MATRICES IN (CSM+RT).
c end of STAGE III
c----------------
      RETURN
      END
c***********************************************************************
c multi-colour disk 
      SUBROUTINE BBMULTI(ENER,DIMULT,NN) 

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION ENER(1),DIMULT(1)

      DOUBLE PRECISION CONST, EY, EYT, SUM, TT, X, YX

      INTEGER NN, IA, IX

      DOUBLE PRECISION CMAX, C83, C53, C13
      DATA CMAX/7D1/
      DATA C83/2.6666666667D0/,C53/1.666666667D0/,C13/0.333333333D0/ 

      CONST=C83
      DO 10 IX=1,NN 
         X=ENER(IX) 
         YX=X*YPH 
         EYT=0D0
         IF(YX.lt.CMAX) EYT=DEXP(-YX) 
         SUM=0D0 
         DO 20 IA=1,N0  
            TT=UL(IA)+YX  
            EY=1D0 
            IF(TT.lt.CMAX) EY=1D0/(1D0-DEXP(-TT))    
            SUM=SUM+AL(IA)*EY*TT**C53 
 20      CONTINUE 
         DIMULT(IX)=SUM*CONST*YX**C13*EYT  
 10   CONTINUE 
      RETURN
      END

c***********************************************************************
C Incident black body spectrum:
      SUBROUTINE BLACKB(ENER,BLINT,N) 

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION ENER(1),BLINT(1) 

      DOUBLE PRECISION X, EYPH, YX

      INTEGER N, IX

      DO 10 IX=1,N 
         X=ENER(IX)
         EYPH=0D0
         YX=YPH*X 
         IF(YX.LT.7D1) EYPH=DEXP(-YX) 
         BLINT(IX)=YX**3*EYPH/(1D0-EYPH)
 10   CONTINUE 
      RETURN
      END 

c***********************************************************************
C  calculation of  frequencies 
      SUBROUTINE FREQPOIN(XLMIN,HX)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION XLMIN, HX, DLN, XL, X, XEV

      INTEGER IX

      DOUBLE PRECISION CKEV
      DATA CKEV/0.0019569472D0/
      DLN=DLOG(1D1) 
      XL=XLMIN  
      DO 15 IX=1,MAXFRE
         XLOG(IX)=XL 
         X=DEXP(XL*DLN)
         XEN(IX)=X                
         XEV=X/CKEV             
         XKEV(IX)=XEV          
         XL=XL+HX             
 15   CONTINUE            
      RETURN
      END  
c***********************************************************************
c***********************************************************************
C  calculation of cross-section S0(x) and auxiliary quantities
      SUBROUTINE CROSCOMP

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION XL, X, YXX, AB, X1

      INTEGER IX, JX, ISX

      DOUBLE PRECISION EPS
      DATA EPS/1D-14/

      DO 15 IX=1,MAXFRE
         XL=XLOG(IX)
         X=XEN(IX) 
         IF(IDISTF.eq.0) THEN 
c Maxwellian electrons
            CALL MXWELL(X,AB,EPS) 
            DO 30 JX=1,IX 
               X1=XEN(JX)
               ISX=JX+IX*(IX-1)/2 
               YXX=YE*(X-X1) 
               EXXY(ISX)=0D0
               IF(YXX.LT.7D1) EXXY(ISX)=DEXP(-YXX) 
 30         CONTINUE
         ELSEIF(IDISTF.eq.2) THEN
c Power-law electrons
            CALL PLCROS(X,AB,EPS) 
         ELSEIF(IDISTF.eq.3.or.IDISTF.eq.4) THEN
c     Maxwellian electrons + Power-law tail
            CALL PMCROS(X,AB,EPS) 
         ENDIF 
         S0(IX)=AB
 15   CONTINUE    
      RETURN
      END   
c***********************************************************************

c***********************************************************************
C calculation of cross-section for Power-law electron distribution 
c eq. (3.4.1) in Nagirner & Poutanen 1994 
      SUBROUTINE PLCROS(X,AB,EPS)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, AB, EPS, DL21, SUM, DLU, DLG, GP, GPZ

      INTEGER I

      DOUBLE PRECISION PS0JP
      EXTERNAL PS0JP

      DOUBLE PRECISION PI
      DATA PI/3.141592653589792384D0/ 

      DL21=DLOG(GMAX/GMIN)
c integral over electron distribution 
c loop over gamma 
      SUM=0D0 
      DO 10 I=1,NG
         DLU=UG(I)*DL21
         DLG=DLU+DLOG(GMIN)
         G=GMIN*DEXP(DLU)
         CALL COMWG
         GP=DEXP(-PLIND*DLG)
         GPZ=GP/Z*AG(I) 
c function of x and gamma (3.3.3) 
         SUM=SUM+PS0JP(X,EPS)*GPZ  
 10   CONTINUE 
      AB=SUM*DNORPOW*DL21*4D0*PI  
      RETURN
      END
c***********************************************************************
C calculation of cross-section for Maxwellian electrons + Power-law tail
c eq. (3.4.1) in Nagirner & Poutanen 1994 
      SUBROUTINE PMCROS(X,AB,EPS)

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, AB, EPS
      DOUBLE PRECISION DL21, DLGM, SUM1, DLU, DLG, GP, GPZ, S, DGY, EGA
      DOUBLE PRECISION DGZ

      INTEGER I

      DOUBLE PRECISION PS0JP
      EXTERNAL PS0JP

      DOUBLE PRECISION PI
      DATA PI/3.141592653589792384D0/ 

      DL21=DLOG(GMAX/GMIN)
      DLGM=DLOG(GMIN) 
      SUM1=0D0 
      if(IDISTF.eq.4) then 
c integral over power-law part of the electron distribution 
      SUM1=0D0 
      DO 10 I=1,NG
         DLU=UG(I)*DL21
         DLG=DLU+DLGM
         G=GMIN*DEXP(DLU)
         CALL COMWG
         GP=DEXP(-PLIND*DLG)
         GPZ=GP/Z*AG(I) 
c function of x and gamma (3.3.3) 
         SUM1=SUM1+PS0JP(X,EPS)*GPZ  
 10   CONTINUE 
      SUM1=SUM1*DNORPOW*DL21*4D0*PI  
      endif
c integral over Maxwellian part
      S=0D0
      DO 20 I=1,NG
         DGZ=UG(I)*DLGM
         G=DEXP(DGZ)
         DGY=(G-1D0)*YE
         EGA=0D0
         IF(DGY.lt.7d1) EGA=DEXP(-DGY)
         CALL COMWG
         S=S+PS0JP(X,EPS)*AG(I)*EGA*G
 20   CONTINUE
      S=S*DNORMAX*DLGM*YKY*YE
      AB=SUM1+S
      RETURN
      END
c***********************************************************************
C calculation of cross-section, mean energy of scattered photon, 
C dispersion of energy of scattered photon, and radiative pressure for   
C Maxwellian electron distribution. 
C THIS VERSION CALCULATES ONLY THE TOTAL  
C C S CROSS-SECTION:   AB (in units of Thomson 
C                                         cross-section);  
C INPUT: 
C 1. X - PHOTON ENERGY IN UNITS OF ELECTRON REST MASS m_ec^2 ; 
C 2. E - RELATIVE ERROR; 
C OUTPUT: 
C 3. AB - TOTAL CS CROSS-SECTION (in units of Thomson 
C         cross-section;  
C  
      SUBROUTINE MXWELL(X,AB,E)                              

      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION X, AB, E
      DOUBLE PRECISION U, AEL, AGL, G1, GY1, YG1, DZ1, DR1, Z1, GP1, ZP1
      DOUBLE PRECISION GM1, ZM1, AB1, AB2, DR2, DZ2, G2, GM2, GP2, GY2
      DOUBLE PRECISION SSM, SSP, YG2, Z2, ZM2, ZP2

      INTEGER L

      DOUBLE PRECISION DE1
      DATA DE1/3.67879441171442D-01/

c      X2=X*X      
      AB=0D0                                   
C      XM=0D0                             
C      QM=0D0                                   
C      PM=0D0                                   
      DO 600 L=1,N0                            
         U=UL(L)+1D0                           
         AEL=AE(L)                             
         AGL=AL(L)                             
         G1=1D0+AU(L)*AU(L)*Y1                         
         GY1=G1+Y1                             
         YG1=G1*(G1+YR1)+YR2                   
         DZ1=DSQRT(1D0+G1)                     
         DR1=1D0/DZ1                           
         Z1=AU(L)*DZ1*DY                           
         GP1=G1+Z1                            
         ZP1=GP1*GP1                          
         GM1=1D0/GP1                          
         ZM1=GM1*GM1                          
c         CALL QSN(X*GP1,0,SSP,SSP1,SSP2,SP1,SP2,SP3,SP4,SP5,SP6,SP7,E) 
c         CALL QSN(X*GM1,0,SSM,SSM1,SSM2,SM1,SM2,SM3,SM4,SM5,SM6,SM7,E) 
         CALL QSN(X*GP1,0,SSP,E) 
         CALL QSN(X*GM1,0,SSM,E) 
         AB1=(ZP1*SSP+ZM1*SSM)*DR1                                    
C         XM1=((GP1*ZP1*(GY1*SP1+X*SP2)+GM1*ZM1*(GY1*SM1+X*SM2)))*DR1  
C         QM1=(ZP1*(0.5D0*SP6-GP1*(GY1*SP5+GP1*(SP7-YG1*SP4)))+        
C     +        ZM1*(0.5D0*SM6-GM1*(GY1*SM5+GM1*(SM7-YG1*SM4))))*DR1    
C         PM1=(ZP1*ZP1*SP1+ZM1*ZM1*SM1)*DR1                       
         G2=1D0+U*Y1                                             
         GY2=G2+Y1                                              
         YG2=G2*(G2+YR1)+YR2                                     
         DZ2=DSQRT(U*(1D0+G2))                                   
         DR2=1D0/DZ2                                             
         Z2=DZ2*DY                                               
         GP2=G2+Z2                                               
         ZP2=GP2*GP2                                             
         GM2=1D0/GP2                                             
         ZM2=GM2*GM2                                             
c         CALL QSN(X*GP2,0,SSP,SSP1,SSP2,SP1,SP2,SP3,SP4,SP5,SP6,SP7,E)
c         CALL QSN(X*GM2,0,SSM,SSM1,SSM2,SM1,SM2,SM3,SM4,SM5,SM6,SM7,E)
         CALL QSN(X*GP2,0,SSP,E) 
         CALL QSN(X*GM2,0,SSM,E) 
         AB2=(ZP2*SSP+ZM2*SSM)*DR2                                   
C         XM2=(GP2*ZP2*(GY2*SP1+X*SP2)+GM2*ZM2*(GY2*SM1+X*SM2))*DR2   
C         QM2=(ZP2*(0.5D0*SP6-GP2*(GY2*SP5+GP2*(SP7-YG2*SP4)))+       
C     +        ZM2*(0.5D0*SM6-GM2*(GY2*SM5+GM2*(SM7-YG2*SM4))))*DR2   
C         PM2=(ZP2*ZP2*SP1+ZM2*ZM2*SM1)*DR2 
         AB=AB+AGL*(AEL*AB1+DE1*AB2) 
C         XM=XM+AGL*(AEL*XM1+DE1*XM2) 
C         QM=QM+AGL*(AEL*QM1+DE1*QM2) 
C         PM=PM+AGL*(AEL*PM1+DE1*PM2) 
  600 CONTINUE  
      AB=AB*DK2Y                     
C      XM=XM*DK2Y                     
C      QM=QM*DK2Y                     
C      PM=PM*DK2Y                     
      RETURN                         
      END                            


C  S- functions. if K=0, then only  SS -- mean cross-section; 
C if  K=1, then  SS1 &  S1 -- mean frequency;
C if K=2, all functions -- and the square of frequency
C x -- small
c      SUBROUTINE QSN(X,K,SS,SS1,SS2,S1,S2,S3,S4,S5,S6,S7,E) 
      SUBROUTINE QSN(X,K,SS,E) 

      IMPLICIT NONE

      DOUBLE PRECISION X, SS, E
      DOUBLE PRECISION A, B, D1, D2, D3, D4, D5, DLX0, DN, DN1
      DOUBLE PRECISION DN124, DN13, DN134, DN14, DN2, DN23, DN234
      DOUBLE PRECISION DN3, DN4, S1, S2, S3, S4, S5, S6, S7, SS1
      DOUBLE PRECISION SS2, X0, X02, X2, X2R, XR

      INTEGER K

      DOUBLE PRECISION XT
      DATA XT /4D-1/
                                  
      IF (X.LT.XT) GO TO 110                               
      X0=1D0/(1D0+2D0*X)                                   
      X2=X*X                                               
      XR=1D0/X                                             
      DLX0=-DLOG(X0)                                       
      X02=X0*X0                                            
      X2R=XR*XR                                            
      SS=0.375D0*(4D0-(2D0*XR+2D0-X)*DLX0+2D0*X2*(1D0+X)*X02)*X2R    
      IF(K.EQ.0) RETURN                                              
      SS1=0.125D0*X2R*XR*(3D0*DLX0+4D0*X2-4.5D0*X-X*X0*(1.5D0+X*X02))
      S1=(SS-SS1)*XR                             
      S2=(SS1-S1)*XR                             
      IF(K.EQ.1) RETURN                          
      SS2=(1D0+X*(4D0+X*(7D0+4.5D0*X)))*X02*X02
      S3=(SS1-SS2)*XR                                
      S4=(S1-S3)*XR                            
      GO TO 160                                      
  110 A=1D0                                          
      B=-2D0*X                                       
      DN=0D0                                         
      SS=0D0                                         
      IF (K.EQ.0) GO TO 120                          
      SS1=0D0                                        
      S1=0D0                                         
      S2=0D0                                         
      IF (K.EQ.1) GO TO 120                          
      SS2=0D0                                        
      S3=0D0                                         
      S4=0D0                                         
  120 DN1=DN+1D0                                     
      DN2=DN+2D0                                     
      DN3=DN+3D0                                     
      D1=1D0/DN1                                     
      D2=1D0/DN2                                     
      D3=1D0/DN3                                     
      DN4=DN+4D0                                     
      DN23=DN2*DN3                                   
      DN234=0.25D0*DN23*DN4                          
      SS=SS+A*(DN2+2D0*D1+8D0*D2-16D0*D3)            
      IF (K.EQ.0) GO TO 140                          
      D4=1D0/DN4                               
      D5=1D0/(DN4+1D0)                               
      DN13=DN1*DN3                             
      DN14=DN1*DN4                       
      SS1=SS1+A*(DN23-6D0+24D0*D3)                   
      S1=S1+A*(DN13-6D0-6D0*D2-24D0*D3+72D0*D4)      
      S2=S2+A*(DN14-6D0-12D0*D3-72D0*D4+144D0*D5)    
      IF (K.EQ.1) GO TO 140                          
      DN134=0.25D0*DN13*DN4                          
      DN124=0.25D0*DN2*DN14                          
      SS2=SS2+A*(DN234+2D0-DN)                       
      S3=S3+A*(DN134+7D0-DN-24D0*D4)                 
      S4=S4+A*(DN124+12D0-DN+6D0*D3+24D0*D4-96D0*D5) 
  140 IF (DABS(4D0*A)*DN234.LT.E) GO TO 150          
      DN=DN+1D0                                      
      A=A*B                                          
      GO TO 120                                      
  150 SS=SS*0.375D0                                  
      IF (K.EQ.0) RETURN                             
      SS1=SS1*0.125D0                                
      S1=S1*0.25D0                                   
      S2=S2*0.25D0                                   
      IF (K.EQ.1) RETURN                             
      SS2=SS2*0.125D0                          
      S3=S3*0.25D0                                   
      S4=S4*0.5D0                                    
  160 S5=3D0*S4-4D0*S3                               
      S7=S3-0.5D0*S4                                 
      S6=2D0*SS2-6D0*S7                              
      RETURN                                         
      END                                            


      SUBROUTINE COMWY(EPS)                                           
C*** CALCULATES SOME FUNCTION OF ELECTRON TEMPERATURE. 
C* * INPUT: 
C* * 1. Y - INVERSE TEMPERATURE (m_c^2/kT_e) .
C* * OUTPUT: 
C* * 2. YR= 1/Y, YR1=2/Y, YR2=2/(Y*Y), DY=1./SQRT(Y)  
C* * 3. YKY=EXP(-Y)/K_2(Y) 
C* * 4. DK2Y=YKY/(2.*SQRT(Y))
C* * 5. DK02=K_0(Y)/K_2(Y) 
C* * 6. DK12=K_1(Y)/K_2(Y) 
C* * 7. GAMMEAN=3D0/Y + K_1(y)/K_2(y) 
C* * 8. GAMSQUA=1D0+3/y*[K_3(y)/K_2(y)]  
C***
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION DK(20)  

      DOUBLE PRECISION EPS, C21, DGPMI, DGY, DGZ, DK2, DK2R, DLGM
      DOUBLE PRECISION DLGMAMI, EDY, EGA, EGPMI, GAM1, GMID, GMY
      DOUBLE PRECISION P1, PI2, PI4, SUM1, YKYEGA, YPI, ZMIN

      INTEGER I

      DOUBLE PRECISION PI, PIQR
      DATA PI/3.141592653589792384D0/
C    PIQR=1./SQRT(2.*PI)
      DATA PIQR/0.39894228040143267794D0/

      IF(IDISTF.eq.0.or.IDISTF.eq.3.or.IDISTF.eq.4) THEN 
         PI2=2D0/PI                                                    
         CALL QDBESK(YE,4,DK,8D0)
         DK2=DK(3) 
         Y1=1D0/YE   
         G0=-Y1*DLOG(EPS)        
         IF(YE.LT.1D0) G0=-Y1*DLOG(EPS*G0) 
         DK02=DK(1)/DK2
         DK12=DK(2)/DK2
         GAMMEAN=3D0*Y1 + DK12 
         ZMEAN=DSQRT(GAMMEAN*GAMMEAN-1D0)
         GAM1=GAMMEAN-1D0 
         GAMSQUA=1D0+3D0/YE*DK(4)/DK2 
         YR1=2D0*Y1             
         YR2=YR1*Y1             
         DY=DSQRT(Y1)           
         DK2R=1D0/DK2            
         IF(YE.GT.8D0) THEN 
            YPI=DSQRT(PI2*YE)        
            YKY=YPI/DK2           
            DK2Y=PIQR*DK2R          
         ELSE 
            EDY=DEXP(-YE)            
            YKY=EDY/DK2      
            DK2Y=0.5D0*EDY*DY*DK2R  
         ENDIF 
      ENDIF 

      IF(IDISTF.eq.2.or.IDISTF.eq.4.and.GMIN.gt.GMAX) THEN 
         GMID=GMIN
         GMIN=DMIN1(GMIN,GMAX)
         GMAX=DMAX1(GMID,GMAX)
      ENDIF 
      IF(IDISTF.eq.4.and.GMIN.eq.GMAX) IDISTF=3 

      SUM1=0D0
      IF(IDISTF.eq.3.or.IDISTF.eq.4) THEN 
         DLGM=DLOG(GMIN)
         DO 10 I=1,NG
            DGZ=UG(I)*DLGM
            G=DEXP(DGZ)
            DGY=(G-1D0)*YE
            EGA=0D0
            IF(DGY.lt.7d1) EGA=DEXP(-DGY)
            CALL COMWG
            SUM1=SUM1+Z*G2*AG(I)*EGA
 10      CONTINUE
         SUM1=SUM1*DLGM*YE*YKY
      ENDIF

c normalization constant of power-law part 
      IF(IDISTF.eq.2.or.IDISTF.eq.4) THEN 
         DLGMAMI=DLOG(GMAX/GMIN)
         P1=PLIND-1D0
         PI4=4D0*PI
         DNORPOW=1D0/(PI4*DLGMAMI)
         IF(DABS(P1).GT.1D-5) DNORPOW=P1/PI4*DEXP(P1*DLOG(GMIN))/
     *        (1D0-DEXP(-P1*DLGMAMI))
      ENDIF
 
c normalization constant for maxwellian part  
      IF(IDISTF.eq.3) THEN 
         DNORMAX=1D0/SUM1 
      ENDIF

      IF(IDISTF.eq.4) THEN 
         DLGMAMI=DLOG(GMAX/GMIN)
         ZMIN=DSQRT(GMIN*GMIN-1D0) 
         GMY=YE*(GMIN-1D0) 
         EGA=0D0 
         IF(GMY.lt.7d1) EGA=DEXP(-GMY) 
         YKYEGA=YKY*EGA
         DGPMI=(1D0+PLIND)*DLOG(GMIN) 
         EGPMI=0D0
         IF(DGPMI.lt.5D1) EGPMI=DEXP(DGPMI)
c ratio of normalization factors: Power-law/Maxw 
         C21=YE*YKYEGA*ZMIN*EGPMI/(PI*4D0)
         DNORMAX=1D0/(SUM1+C21/DNORPOW)
         DNORPOW=DNORMAX*C21 
      ENDIF

      RETURN   
      END      


C
C AUXILARY QUANTITIES DEPENDING ON X,X1 AND MU ; 
C FIRST CALL COMX(X,X1) SUBROUTINE 
C 
      SUBROUTINE COMXM(QM)                                 

      IMPLICIT NONE

      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION QM

      DOUBLE PRECISION EPS
      DATA EPS/1D-10/ 

      QO=T*QM    
      Q2=DX2+2D0*QO                   
      Q=DSQRT(Q2)                     
      IF(QM.gt.EPS) THEN 
         QR=1D0/QO                       
         QOT=2D0*QR                      
         QT=2D0/Q                        
         QMT=2D0/QM                      
         PM=2D0-QM                       
         RS=PM/QM                        
         QS=T*PM                         
         GST=(DX+Q*DSQRT(1D0+QOT))*5D-1
         GST1=GST-1D0
      ELSE 
         GST=1D0 
         GST1=0D0
         QT=0D0 
      ENDIF 
      RETURN                                          
      END                           


C 
C AUXILARY QUANTITIES DEPENDING ON GAMMA 
C 
      SUBROUTINE COMWG            
      IMPLICIT NONE

      INCLUDE 'comppsj.inc'

      G1=G-1D0                     
      G2=G*G                       
      GT=2D0*G                     
      Z2=G1*(G1+2D0)               
      Z=DSQRT(Z2)                  
      RETURN                       
      END                          

      SUBROUTINE COMWZ
      IMPLICIT NONE

      INCLUDE 'comppsj.inc'

      Z2=Z*Z
      G2=Z2+1D0
      G=DSQRT(G2)
      G1=G-1D0
      GT=2D0*G
      RETURN
      END

C 
C AUXILARY QUANTITIES DEPENDING ON  X AND X1
C 
      SUBROUTINE COMX(X,X1)                                 
      IMPLICIT NONE

      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, X1, GA, TOO, TT, TTT, XMAX, XMIN, XXX

      DOUBLE PRECISION C13
      DATA C13/0.33333333333333D0/

      XS=1D0+2D0*X                                          
      XA=X/XS                                               
      XS1=1D0+2D0*X1                                        
      XA1=X1/XS1                                            
      XX=2D0*XA*X     
      XX1=-2D0*XA1*X1 
      SX=X+X1         
      DX=X-X1         
      T=X*X1          
      TT=T*T          
      TR=1D0/T        
      XM2=SX*DX       
      DX2=DX*DX       
      TR2=TR*TR       
      TRT=-2D0*TR     
      TRF=4D0*TR      
      TRR=TR+1D0      
      TRS=5D-1*TR2                    
      TTT=DX2*TR2-TRF                 
      TOO=1D0-TRT                     
      XMAX=X                          
      IF(X.LT.X1) XMAX=X1             
      XMIN=X                          
      IF(X.GT.X1) XMIN=X1             
      RB=4D0/XMAX                     
      RWB=RB-4D0*C13*XMIN/(XMAX*XMAX) 
      GA=0.5D0*(DX+SX*DSQRT(TRR))     
      XXA=X-XA1                       
      XXA1=X1-XA                      
      XXX=(1D0-2D0*X)/XS              
      GA1=GA-1D0                      
      IF(XXA.GT.-0.1D0)      GA1=0.25D0*XS1*XS1*TR*XXA*XXA/(GA+1D0-DX)  
      IF(XXA1.GT.-0.1D0*XXX) GA1=0.25D0*XS*XS*TR*XXA1*XXA1/(GA+1D0)+DX  
      IF(XXA.LE.0D0) NR=1                 
      IF(XXA1.LE.0D0) NR=2                
      IF((XXA.GT.0D0).AND.(X.LE.X1)) NR=3 
      IF((XXA1.GT.0D0).AND.(X1.LT.X)) NR=4
      IF(NR.LT.0) THEN
         GM1=GA1
      ELSEIF(NR.EQ.0) THEN
         GM1=0D0                             
      ELSE
         GM1=DX                             
      ENDIF
      RETURN                             
      END                                


C 
C REDISTRIBUTION MATRIX : S(X,X1,MU)
C INPUT: 
C 1. X, X1 - ENERGIES OF SCATTERED AND  INITIAL PHOTONS, RESPECTIVELY; 
C 2. QM=1-MU , AND MU - THE COSINE OF THE SCATTERING ANGLE; 
C OUTPUT: 
C 3. S, SI, SQ, SU, SV, SA=SQ-SU  (X,X1,MU) 
C 
      SUBROUTINE REDFUNC(X,X1,QM,S) 
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, X1, QM, S, SM, SP

      DOUBLE PRECISION xold, x1old
      DATA xold/0./,x1old/0./ 

      if(x.ne.xold.or.x1.ne.x1old) then 
         CALL COMX(X,X1) 
         xold=x 
         x1old=x1  
      endif  

      CALL COMXM(QM) 

c Maxwellian electron distribution 
      IF(IDISTF.eq.0) THEN 
         CALL MAXGALA(X,X1,QM,S) 
      ENDIF
c Maxwellian electron distribution with cut-off 
c and Maxwellian part of RF for Maxw + power-law
      IF(IDISTF.eq.3.or.IDISTF.eq.4) THEN 
         IF(GMIN.gt.GST) THEN 
            CALL MAXCUT(X,X1,QM,SM)  
         ELSE 
            SM=0D0 
         ENDIF
      ENDIF

      IF(IDISTF.eq.3) S=SM 
c power-law part of the hybrid electron distribution 
      IF(IDISTF.eq.4) THEN 
         IF(GMAX.gt.GST) THEN 
            CALL POWMATR(X,X1,QM,SP)
            S=SM+SP 
         ELSE 
            S=0D0 
         ENDIF 
      ENDIF                                          
c power-law electron distribution 
      IF(IDISTF.eq.2) CALL POWMATR(X,X1,QM,S) 
      RETURN                                               
      END                                                  


C 
C SCATTERING MATRIX FOR MAXWELLIAN ELECTRONS 
C 
      SUBROUTINE MAXGALA(X,X1,QM,S) 
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, X1, QM, S
      DOUBLE PRECISION CONST, DGZ, EXGST, R, YGST

      INTEGER I

C PI332=3/(32*PI)
      DOUBLE PRECISION PI332
      DATA PI332/0.2984155182973039D-01/

      EXGST=0D0                                            
      YGST=YE*GST1                                    
      IF(YGST.LT.7D1) EXGST=DEXP(-YGST)     
      CONST=PI332*YKY*EXGST
      S=0D0
      DO 10 I=1,N0
         DGZ=UL(I)*Y1 
         G=DGZ+GST
         CALL COMWG
         CALL RDSTRB(X,X1,QM,R) 
         S=S+R*AL(I)                                       
 10   CONTINUE 
      S=S*CONST 
      RETURN                                               
      END                                                  

C SCATTERING MATRIX FOR MAXWELLIAN ELECTRONS with cut-off at GMIN 
C 
      SUBROUTINE MAXCUT(X,X1,QM,S)
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, X1, QM, S
      DOUBLE PRECISION CONST, DGMS, DGST, DGY, DGZ, EGA, EXGST, GSTY
      DOUBLE PRECISION R, SUM1, YGST

      INTEGER I

C PI332=3/(32*PI)
      DOUBLE PRECISION PI332
      DATA PI332/0.2984155182973039D-01/
   
      S=0D0
      IF(GMIN.gt.GST) THEN 

         EXGST=0D0
         GSTY=GST+Y1
         YGST=YE*GST1
         IF(YGST.LT.7D1) EXGST=DEXP(-YGST)
         CONST=PI332*YKY*EXGST*DNORMAX  
         SUM1=0D0

         DGMS=DLOG(GMIN/GST) 
         DGST=DLOG(GST) 
         DO 10 I=1,NG 
            DGZ=UG(I)*DGMS+DGST 
            G=DEXP(DGZ) 
            DGY=(G-GST)*YE 
            EGA=0D0 
            IF(DGY.lt.7d1) EGA=DEXP(-DGY) 
            CALL COMWG 
            CALL RDSTRB(X,X1,QM,R) 
            SUM1=SUM1+R*AG(I)*EGA*G 
 10      CONTINUE 
         SUM1=SUM1*DGMS*YE 
         S=SUM1*CONST 
      ENDIF 

      RETURN                                               
      END 
                                                 

C Procedure calculates the S(x,x1,mu) matrix for the power-law electron 
C distribution:
C P- power-law index; GMIN,GMAX- limits of distribution; 
C X,X1 - frequencies in units of m_ec2; QM=1-mu; 
C GST=gamma_*; 
C need in advance calling the COMXX(X,X1) prosedure
      SUBROUTINE POWMATR(X,X1,QM,S)
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'
      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, X1, QM, S
      DOUBLE PRECISION AGI, DL1, DL21, DLG, DLU, GM, GMM, GP, GP1, GPZ
      DOUBLE PRECISION R, SUM1, SUM2, ZM, ZMM

      INTEGER I

      DOUBLE PRECISION C38, C23, GCR
      DATA C38/0.375D0/
      DATA C23/0.6666666666666667D0/
      DATA GCR/2D0/ 

      S=1D-40
      IF(GST.GE.GMAX) RETURN
      IF(IREDF.eq.1) R=C23*QT
      GM=DMAX1(GST,GMIN)

      SUM1=0D0
      if(GCR.lt.GMAX) then
c int d(ln g)/z g^-p R(g)
         GMM=DMAX1(GM,GCR)
         DL21=DLOG(GMAX/GMM)
         DL1=DLOG(GMM)
         DO 10 I=1,NG
            AGI=AG(I)
            DLU=UG(I)*DL21
            DLG=UG(I)*DL21+DL1
            G=DEXP(DLG)
            CALL COMWG
            GP=DEXP(-PLIND*DLG)
            GPZ=GP/Z*AGI
            IF(IREDF.eq.0) CALL RDSTRB(X,X1,QM,R)
            SUM1=SUM1+R*GPZ
 10      CONTINUE
         SUM1=SUM1*DL21
      endif

      SUM2=0D0
      if(GCR.gt.GM) then
c int dz/g^2 g^-p R(g)
         GMM=DMIN1(GCR,GMAX)
         ZMM=DSQRT(GMM*GMM-1D0)
         ZM=DSQRT(GM*GM-1D0)
         DL21=ZMM-ZM
         DO 20 I=1,NG
            AGI=AG(I)
            Z=UG(I)*DL21+ZM
            CALL COMWZ
            DLG=DLOG(G)
            GP1=DEXP(-(PLIND+2D0)*DLG)
            IF(IREDF.eq.0) CALL RDSTRB(X,X1,QM,R)
            SUM2=SUM2+R*GP1*AGI
 20      CONTINUE
         SUM2=SUM2*DL21
      endif

      S=SUM1+SUM2
      S=S*C38*DNORPOW
      RETURN  
      END    
c ****end of POWMATR 

C 
C SCATTERING MATRIX R(X,X1,MU,GAMMA) 
C needed first to call subroutines: 
C COMX(X,X1); COMXM(QM); COMWG 
C  INPUT: 
C  1. X,X1 - ENERGIES OF  SCATTERED AND INITIAL PHOTONS;  
C  2. QM=1-MU (MU- COSINE OF SCATTERING ANGLE); 
C  3.  GAMMA=G (in COMMON/WZ/) 
C  OUTPUT: R,RI,RQ,RU,RV, RF=RQ-RU 
C 
      SUBROUTINE RDSTRB(X,X1,QM,R) 
      IMPLICIT NONE

      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, X1, QM, R
      DOUBLE PRECISION A, A1, C, CD, D, DFG, DG, GX, QQQ, QSR, QU
      DOUBLE PRECISION RF, RFO, RG, RGO, RHO, U, U2, UQ, UQQ, UV
      DOUBLE PRECISION UVQ, V, VQ, VR

      DOUBLE PRECISION EPS
      DATA EPS/1D-10/
                   
      IF(Z.LE.0D0)  stop 
      IF(QM.GT.2D0) stop 
      IF(X.LE.0D0)  stop
      IF(X1.LE.0D0) stop
      IF(QM.lt.EPS) THEN 
         R=0D0 
         RETURN  
      ELSE 
         A=DSQRT(Z2-(GT-X)*X+QMT)        
         A1=DSQRT(Z2+(GT+X1)*X1+QMT)     
         U=SX*(GT-DX)/(A+A1)             
         U2=U*U                                  
         V=A*A1                                  
         VR=1D0/V                                
         VQ=QR*VR                                
         UV=U*VR                                 
         UVQ=UV*QR                               
         QU=U+Q                                  
         UQ=U-Q                                   
         IF(UQ.LT.3D-1) GO TO 500               
         QSR=1D0/QS
         UQQ=U2-Q2                               
         RF=0.5D0*UV*VQ*VQ*(U2*UQQ+V*(5D0*U2-3D0*Q2)) 
         RG=QT+UV*(1D0-QOT)                      
         R=RF+RG                                 
         RETURN                                  
 500     GX=G-DX                                 
         DG=G-GST
         C=2D0/(G*GX+RS+T-QO+V)                  
         D=2D0*(GX+GST)*DG                        
         CD=C*D                                  
         UQQ=QS*CD                               
         QQQ=CD/QU                               
         RFO=0.5D0*UV*VQ*VQ*UQQ*(U2+5D0*V)       
         RGO=QT+UV                                
         RHO=UVQ*CD                              
         DFG=UV*VQ*((1D0-QS*C)*D+RS*(2D0+QO))    
         R=RFO+RGO-DFG                           
      ENDIF 
      RETURN                                 
      END                                                              


c***********************************************************************
c**************************************************
c function Psi_0 (x)*G*Z used in calculation of the total
c cross-section for Compton scattering
c eq. (3.3.3) in Nagirner & Poutanen 1994 
      DOUBLE PRECISION FUNCTION PS0JP(X,E) 
      IMPLICIT NONE

      INCLUDE 'comppsj.inc'

      DOUBLE PRECISION X, E, ARG1, ARG2, GZP, GZP2

      DOUBLE PRECISION PSIJP
      EXTERNAL PSIJP

      GZP=G+Z 
      GZP2=GZP**2  
      ARG1=X*GZP 
      ARG2=X/GZP 
      PS0JP=(PSIJP(ARG1,E)*GZP2-PSIJP(ARG2,E)/GZP2)/4D0  
      RETURN
      END
c**************************************************
c function psi_10 (x) used in calculation of the total 
c cross-section for Compton scattering 
c eq. (3.3.10) in Nagirner & Poutanen 1994 
      DOUBLE PRECISION FUNCTION PSIJP(X,E) 
      IMPLICIT NONE

      DOUBLE PRECISION X, E
      DOUBLE PRECISION A, B, D1, D2, D3, D33, DLX0, DN, DN1, DN2, DN3
      DOUBLE PRECISION DNQ, G0, P10, X0, X02, X03, X2, XR, XR2

      DOUBLE PRECISION GINT
      EXTERNAL GINT

      DOUBLE PRECISION XT
      DATA XT /4D-1/

      IF(X.LT.XT) GO TO 210 
      X0=1D0/(1D0+2D0*X)   
      X02=X0*X0           
      X03=X0*X02         
      X2=X*X            
      XR=1D0/X         
      DLX0=-DLOG(X0)  
      G0=GINT(X,E)      
      XR2=XR*XR     
      P10=0.75D0*XR2*((X+4.5D0+2D0*XR)*DLX0-4D0-X+X2*X0-2D0*G0)       
      GO TO 240 
  210 DN=0D0    
      A=1D0     
      B=-2D0*X  
      P10=0D0      
  220 DN1=DN+1D0 
      DN2=DN+2D0
      DN3=DN+3D0  
      DNQ=DN*DN    
      D1=1D0/DN1     
      D2=1D0/DN2    
      D3=1D0/DN3      
      D33=D3*D3         
      P10=P10+A*(1D0+D2*(2D0*D1+8D0*D2-16D0*D3)) 
      IF(DABS(A)*(DNQ+1D0).LT.E) GO TO 230  
      A=B*A                                
      DN=DN+1D0     
      GO TO 220    
  230 CONTINUE 
      P10=P10*0.75D0
  240 PSIJP=P10 
      RETURN  
      END  

      DOUBLE PRECISION FUNCTION GINT(X,E)    
      IMPLICIT NONE

      DOUBLE PRECISION X, E, A, B, DK, DN, U, V, Y, Z

      DOUBLE PRECISION XT, XTR, DL, PIH, PITW
      DATA XT,XTR /0.4D0,0.625D0/ 
      DATA DL/0.6931471805599453094D0/
      DATA PIH,PITW /1.64493406684823D0,0.822467033424113D0/   

      IF (X.LT.XT) GO TO 40                              
      IF (X.GT.XTR) GO TO 40  
      GINT=PITW                  
      IF (X.EQ.0.5D0) RETURN  
      Y=2D0*X                 
      Z=DLOG(Y)               
      GINT=GINT+DL*Z                
      IF (X.GT.0.5D0) GINT=GINT+0.5D0*Z*Z
      IF (X.GT.0.5D0) Y=1D0/Y  
      Y=1D0-Y   
      B=1D0    
      U=0D0     
      DK=1D0    
      V=0D0     
      A=1D0     
   30 A=0.5D0*A 
      V=V+A/DK  
      DK=DK+1D0 
      U=U+B*V/DK  
      B=B*Y       
      IF (B.GT.E) GO TO 30   
      U=U*Y*Y                
      IF (X.GT.0.5D0) GINT=GINT-U  
      IF (X.LT.0.5D0) GINT=GINT+U  
      RETURN                 
   40 IF (X.LT.XT) Y=2D0*X   
      IF (X.GT.XTR) Y=0.5D0/X
      GINT=1D0                 
      A=-Y                 
      DN=1D0     
      B=1D0      
   50 DN=DN+1D0  
      B=A*B      
      GINT=GINT+B/(DN*DN)  
      IF (DABS(B).GT.E) GO TO 50 
      GINT=Y*GINT          
      IF (X.LT.XT) RETURN  
      Z=DLOG(Y)     
      GINT=PIH+0.5D0*Z*Z-GINT     
      RETURN       
      END         


c**************************************************
C linear interpolation:
C abciss:  XX
C ordinate: FUNLOG
C A,B: coefficients of the linear function Ax+B
      SUBROUTINE DLININT(FUNLOG,XLOG,A,B,NMIN,NMAX,MAXFRE)
      IMPLICIT NONE

      INTEGER MAXFRE
      DOUBLE PRECISION FUNLOG(MAXFRE),XLOG(MAXFRE)
      DOUBLE PRECISION A, B, DELTA, DI, S, SX, SXX, SXY, SY, X

      INTEGER NMIN, NMAX, IX

      DOUBLE PRECISION ERR
      DATA ERR/1D30/

      IF(NMAX.le.NMIN) THEN
         A=ERR
         B=ERR
      ELSE
         S=NMAX-NMIN+1
         SX=0. 
         SY=0. 
         SXX=0. 
         SXY=0. 
         DO 15 IX=NMIN,NMAX
            X=XLOG(IX)
            DI=FUNLOG(IX)
            SX=SX+X
            SY=SY+DI
            SXX=SXX+X*X
            SXY=SXY+DI*X
 15      CONTINUE
         DELTA=S*SXX-SX*SX
         B=(SXX*SY-SX*SXY)/DELTA
         A=(S*SXY-SX*SY)/DELTA
      ENDIF
      RETURN
      END
c***********************************************************************
      SUBROUTINE QSPLINE(X,Y,N,YP1,YPN,Y2)
      IMPLICIT NONE

      INTEGER NMAX
      PARAMETER (NMAX=500)

      INTEGER N
      DOUBLE PRECISION X(N),Y(N),Y2(N),U(NMAX)
      DOUBLE PRECISION YP1, YPN, P, QN, SIG, UN

      INTEGER I, K

      IF (YP1.GT..99D30) THEN
         Y2(1)=0.
         U(1)=0.
      ELSE
         Y2(1)=-0.5
         U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.
         Y2(I)=(SIG-1.)/P
         U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
 11   CONTINUE
      IF (YPN.GT..99D30) THEN
         QN=0.
         UN=0.
      ELSE
         QN=0.5
         UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
 12   CONTINUE
      RETURN
      END


      SUBROUTINE QSPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION XA(N),YA(N),Y2A(N)
      DOUBLE PRECISION X, Y, A, B, H

      INTEGER K, KHI, KLO

      KLO=1
      KHI=N
 1    IF (KHI-KLO.GT.1) THEN
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
            KHI=K
         ELSE
            KLO=K
         ENDIF
         GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) CALL XWRITE('Bad XA input.',5)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *     ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
c*********************************************
c polinomial interpolation and extrapolation    
      SUBROUTINE POLINTQ(XA,YA,N,X,Y,DY)
      IMPLICIT NONE

      INTEGER NMAX
      PARAMETER (NMAX=10) 

      INTEGER N
      DOUBLE PRECISION XA(N),YA(N),C(NMAX),D(NMAX)
      DOUBLE PRECISION X, Y, DY, DEN, DIF, DIFT, HO, HP, W

      INTEGER I, M, NS

      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
         DIFT=ABS(X-XA(I))
         IF (DIFT.LT.DIF) THEN
            NS=I
            DIF=DIFT
         ENDIF
         C(I)=YA(I)
         D(I)=YA(I)
 11   CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
         DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.D0)CALL XWRITE('In POLINTQ',5)
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
 12    CONTINUE
       IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
       ELSE
          DY=D(NS)
          NS=NS-1
       ENDIF
       Y=Y+DY
 13   CONTINUE
      RETURN
      END

c**************************************************
c linear interpolation 
      SUBROUTINE DLINEAR(X1,Y1,X2,Y2,X,Y) 
      IMPLICIT NONE

      DOUBLE PRECISION X1, Y1, X2, Y2, X, Y
      DOUBLE PRECISION DX, T

      DX=X2-X1 
      IF(DX.eq.0D0) CALL XWRITE('In DLINEAR',5)
      T=(Y2-Y1)/DX 
      Y=Y1+T*(X-X1) 
      RETURN
      END

c**********************************************************************
      SUBROUTINE QHUNT(XIN,N,X,JLO)
      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION XIN(N)
      DOUBLE PRECISION X

      INTEGER JLO, INC, JHI, JM

      LOGICAL ASCND

      ASCND=XIN(N).GT.XIN(1)
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
         JLO=0
         JHI=N+1
         GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XIN(JLO).EQV.ASCND)THEN
 1       JHI=JLO+INC
         IF(JHI.GT.N)THEN
            JHI=N+1
         ELSE IF(X.GE.XIN(JHI).EQV.ASCND)THEN
            JLO=JHI
            INC=INC+INC
            GO TO 1
         ENDIF
      ELSE
         JHI=JLO
 2       JLO=JHI-INC
         IF(JLO.LT.1)THEN
            JLO=0
         ELSE IF(X.LT.XIN(JLO).EQV.ASCND)THEN
            JHI=JLO
            INC=INC+INC
            GO TO 2
         ENDIF
      ENDIF
 3    IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XIN(JM).EQV.ASCND)THEN
         JLO=JM
      ELSE
         JHI=JM
      ENDIF
      GO TO 3
      END
c**************************************************
      SUBROUTINE DLOCATE(XIN,N,X,J)
      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION XIN(N)
      DOUBLE PRECISION X

      INTEGER J, JL, JU, JM

      JL=0
      JU=N+1
 10   IF(JU-JL.GT.1)THEN
         JM=(JU+JL)/2
         IF((XIN(N).GT.XIN(1)).EQV.(X.GT.XIN(JM)))THEN
            JL=JM
         ELSE
            JU=JM
         ENDIF
         GO TO 10
      ENDIF
      J=JL
      RETURN
      END
c**************************************************
      SUBROUTINE AGUEGALA(NL,E)             
      IMPLICIT NONE

      INCLUDE 'comppsq.inc'

      DOUBLE PRECISION E, DLE, EXPUL

      INTEGER NL, L

      DLE=-DLOG(E)                           
      CALL QFQLAG0(NL,E,UL,AL)               
      DO 200 L=1,NL                          
C     IF(DABS(AL(L)).LT.E) RETURN        
         EXPUL=DEXP(-UL(L))                      
         AU(L)=EXPUL                             
         AE(L)=DEXP(-EXPUL*EXPUL)*2D0                
C     N0=L                                
 200  CONTINUE                               
      N0=NL
      RETURN                                 
      END                                    
c**************************************************
c sin(x)/x  
      DOUBLE PRECISION FUNCTION SINXX(X) 
      IMPLICIT NONE

      DOUBLE PRECISION X, SUM, DN, X2, DEL

      DOUBLE PRECISION EPS, ACCUR
      DATA EPS/1D-2/,ACCUR/1D-14/
 
      if(X.gt.EPS) then 
         SINXX=DSIN(X)/X 
      else
         SUM=1D0 
         DN=2D0 
         X2=X*X 
         DEL=1D0 
         DO WHILE(DABS(DEL).gt.ACCUR) 
            DEL=-DEL*X2/(DN*(DN+1D0)) 
            SUM=SUM+DEL 
            DN=DN+2D0  
         ENDDO 
         SINXX=SUM 
      endif 

      RETURN
      END

      SUBROUTINE QDRGSDO(Z,W,NUM)
      IMPLICIT NONE

      DOUBLE PRECISION Z(1),W(1)
      DOUBLE PRECISION DIV, DP, DQ, FJ, FN, FR, HP, HP2
      DOUBLE PRECISION P, Q, TY, X, Y

      INTEGER NUM, N, I, IND, ITER, J, K, KIND, M, NP

      INTEGER ITRMAX
      DATA ITRMAX/10/

      N=NUM
      IF (N.GT.1) GO TO 20
      IF (N.EQ.1) GO TO 10
      N=5
      GO TO 20
 10   Z(1)=0.D0
      W(1)=2.D0
      RETURN
 20   IND=MOD(N,2)
      K=N/2
      IF(IND.EQ.0) GOTO 40
      HP=1.D0
      FJ=1.D0
      DO 30 J=3,N,2
         FJ=FJ+2.D0
         HP=(FJ-1.D0)*HP/FJ
 30   CONTINUE
      HP2=HP*HP
      KIND=K+IND
      Z(KIND)=0.D0
      W(KIND)=HP2+HP2
 40   FN=DFLOAT(N)
      NP=N+1
      X=1.D0-3.D0/(FN+2.D0)/FN
      DO 100 I=1,K
         M=NP-I
         IF(I.EQ.2) X=(X-1.D0)*5.25D0+1.D0
         IF(I.EQ.3) X=(X-Z(N))*1.6D0+X
         IF(I.GT.3) X=(X-Z(M+2))*3.D0+Z(M+3)
         ITER=0
 50      Y=X-1.D0
         ITER=ITER+1
         TY=2.D0*Y
         FJ=1.D0
         DP=Y
         P=X
         Q=1.D0
         DQ=1.D-0
         DO 60 J=2,N
            DQ=TY*Q+DQ+P
            Q=DQ+Q
            FJ=FJ+1.D0
            DP=(2.D0*FJ-1.D0)*Y*P+DP
            P=DP/FJ+P
 60      CONTINUE
         DIV=P/Q
         X=X-DIV
         FR=DIV/X
         IF(ITER.GT.ITRMAX) GOTO 190
         IF(DABS(FR).GT.1.D-16) GOTO 50
 190     CONTINUE
         Z(M)=X
         W(M)=2.D0/(1.D0-X*X)/Q/Q
         Z(I)=-X
         W(I)=W(M)
 100  CONTINUE
      RETURN
      END

c**********************************************************************
c Gauss-Chebyshev quadrature for
c \int_0^1 d\mu/\sqrt(1-\mu^2) * \sqrt(1-\mu^2) f(\mu)
      SUBROUTINE QCHEB(Z,W,N)
      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION Z(N),W(N)
      DOUBLE PRECISION PIN, ARG

      INTEGER J

      DOUBLE PRECISION PI
      DATA PI/3.141592653589792384D0/

      PIN=PI/2D0/DFLOAT(N)
      ARG=PIN/2D0
      DO 60 J=1,N
         Z(J)=DSIN(ARG)
         W(J)=PIN*DCOS(ARG)
         ARG=ARG+PIN
 60   CONTINUE
      RETURN
      END
c**********************************************************************

C     SUBROUTINE FQLAG0(K,EPS,T,A) CALCULATE WEIGHTS AND POINTS 
C     K-POINTS GAUSSIAN QUADRATURE  
C     FROM INTERVAL  [0,INFINITY) WITH WEIGHT  EXP(-X).
C     EPS- RELATIVE ERROR OF POINTS CALCULATION 
C     T(1:K),A(1:K)-VECTORS OF POINTS AND WEIGHTS 
C     LIT.: A.K.PONOMARENKO,    LENINGRAD UNIVERSITY,1974 
C 
      SUBROUTINE QFQLAG0 (K,EPS,T,A)    
      IMPLICIT NONE

      INTEGER K
      DOUBLE PRECISION T(K),A(K)              
      DOUBLE PRECISION EPS, DP, X, DJ, E2, E3, E1, DJ1

      INTEGER I, J

      I=0                              
      DP=K                             
      X=3.0D0/(1.0D0 + 2.4D0 * DP)     
      GO TO 30                         
   10 X=T(I) + 6.0D0 / (0.4D0 + DP)    
      GO TO 30                         
   20 J=I-1                            
      DJ=J                             
      X=T(I) + (T(I) - T(J)) * (1.0D0 + 2.55D0 * DJ)/1.9D0/DJ
   30 E2=1.0D0                         
      E3=1.0D0 - X                         
      IF (K.EQ.1) GO TO 41                 
      DO 40 J=2,K                          
         E1 = E2                             
         E2 = E3                              
         DJ = J                               
         DJ1= DJ - 1.0D0                      
         E3 = ((DJ + DJ1 - X)*E2 - DJ1 * E1)/DJ
 40   CONTINUE
 41   E1=DP*(E3-E2)                         
      E2 = E3/E1                            
      X  = X*(1.0D0 - E2)                   
      IF(DABS(E2).GE.EPS) GOTO 30
      I=I+1                                 
      T(I)=X                                
      A(I)= X/E1/E1                         
      IF (I.EQ.1) GOTO 10
      IF (I.GE.1 .AND. I.LT.K) GOTO 20
      RETURN                                
      END                                   


C SUBROUTINE  DBESK (X,N,BK,X1) CALCULATE 
C FUNCTION CONNECTED WITH K(N,X) -- BESSEL FUNCTIONS OF ORDER 
C FROM  0 TO N, WHICH PUT TO THE FIRST N+1 ELEMENTS OF MASSIVE BK.
C C   X1 - BOUNDARY OF DIFFERENT  TYPES OF CALCULATIONS
C        (5<=X1<=8). 
C IF  0<=X<=X1  IS CALCULATED  K(N,X), IF   X>X1
C IS CALCULATED FUNCTION:
C
C    K(N,X)*EXP(X)*SQRT(2*X/PI), WHERE PI=3.14... 
C
C USE SUBROUTINES:
C
C    DBI0LE(X)      -8<=X<=8  
C    DBI1LO(X)      -8<=X<=8      DBK0GS(X)    X>=5 
C    DBK0LE(X)       0< X<=8      DBK1GS(X)    X>=5 
C    DBK1LO(X)       0< X<=8  
C
      SUBROUTINE QDBESK(X,N,BK,X1)        
      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION BK(N) 
      DOUBLE PRECISION X, X1, ED

      INTEGER I

      DOUBLE PRECISION DBK0LE, DBI0LE, DBK0GS, DBK1LO, DBI1LO, DBK1GS
      EXTERNAL DBK0LE, DBI0LE, DBK0GS, DBK1LO, DBI1LO, DBK1GS

      DOUBLE PRECISION EILER, C1
      DATA EILER /0.57721 56649 01532 86061D0/
      DATA C1 /1D0/         
                  
C CONTROL OF ORDER:
      IF (N.LT.1) THEN
C REACTION ON THE NEGATIVE ORDER:
         CALL XWRITE('DBESK: N < 1', 5)
         RETURN 
      ENDIF
C CONTROL OF NEGATIVE ARGUMENT (X):
      IF (X.LT.0.0D0) THEN
C REACTION ON THE NEGATIVE ARGUMENT 
         CALL XWRITE('DBESK: X < 0.0D0', 5)
         RETURN 
      ELSEIF (X.EQ.0.0D0) THEN
C CALCULATION IF X=0
         BK(1)=1D0     
         DO 17 I=2,N   
            BK(I)=0D0 
 17      CONTINUE
         RETURN        
      ENDIF
C       CALCULATION OF K(0,X)
      IF (X.LE.X1) BK(1)=DBK0LE(X) 
      IF (X.LE.X1) BK(1)=BK(1)-(EILER+DLOG(X/2D0))*DBI0LE(X)
      IF (X.GT.X1) BK(1)=DBK0GS(X)                          
      IF (N.EQ.1) RETURN                                    
C CALCULATION OF K(1,X)
      IF (X.LE.X1) THEN
         BK(2)=DBK1LO(X) 
         ED=EILER+DLOG(X/2D0) 
         BK(2)=ED*DBI1LO(X)+C1/X-BK(2) 
      ENDIF
      IF (X.GT.X1) BK(2)=DBK1GS(X)               
      IF (N.EQ.2) RETURN                         
C CALCULATION OF K(N,X)
      DO 30 I=3,N      
         BK(I)=BK(I-2)+2*(I-2)*BK(I-1)/X  
   30 CONTINUE                             
      RETURN                               
      END                                        


C SUBROUTINE   DCHEB0 (A,X,N,KEY)  CALCULATE 
C THE VALUES OF COEFFICIENTS OF EXPANSION OVER THE 
C CHEBYSHEV POLYNOMS
C IF ARE KNOWN FIRST  N+1 ELEMENTS OF VESTOR  A 
C 
C           F(N,X)=[K=0,N] A(K) T(L(KEY),X),     
C                                                
C IF        L        KEY                         
C         ______________                         
C                                                
C           N         0                          
C         2N+1        1                                    
C          2N         2                           
C           N (DISPL) 3,                          
C                                                 
C ALGORITHM KLENSHOW
C LITERATURE:Y.LUKE "SPECIAL MATHEMATICAL FUNCTIONS 
C                   AND ITS APPROXIMATIONS"
C                "MIR",MOSCOW  1980 608 P.
C  P.511
      DOUBLE PRECISION FUNCTION QDCHEB0(A,X,N,KEY)   
      IMPLICIT NONE

      DOUBLE PRECISION A(*)
      DOUBLE PRECISION X, B0, B1, B2, Z

      INTEGER N, KEY, I, IER, K, N1

      DOUBLE PRECISION CMAX, C0, C2, C4
      DATA CMAX /1D76/                  
      DATA C0,C2,C4 /0D0,2D0,4D0/       
C
C
C CONTROL OF NUMBER OF ITEMS
C
      IER=0  
      IF (N.LT.0) THEN
C        REACTION ON THE NEGATIVE NUMBER
         QDCHEB0=CMAX    
         IER=1           
      ENDIF

      IF (N.LE.0 .AND. ((KEY.LT.0).OR.(KEY.GT.3)) ) THEN
C   CONTROL OF KEY
C  REACTION ON THE NONPOSITIVE VALUE OF KEY
         QDCHEB0=1D76     
         IER=1
      ENDIF

      IF (IER.NE.0) RETURN                 

      N1=N+1
      Z=0
      IF (KEY.EQ.0) Z=C2*X                 
      IF (KEY.EQ.1) Z=C4*X*X-C2            
      IF (KEY.EQ.2) Z=C4*X*X-C2            
      IF (KEY.EQ.3) Z=C4*X-C2              
      B0=C0                                
      B1=C0                                
      DO 10 I=1,N1                         
         K=N1-I                             
         B2=B1                              
         B1=B0                              
         B0=Z*B1-B2+A(K+1)                  
 10   CONTINUE                             
      IF (KEY.EQ.0) QDCHEB0=B0-X*B1         
      IF (KEY.EQ.1) QDCHEB0=X*(B0-B1)       
      IF (KEY.EQ.2) QDCHEB0=B0-B1*Z/C2      
      IF (KEY.EQ.3) QDCHEB0=B0-B1*Z/C2      
      RETURN                               
      END                                       



C FUNCTION  DBK0LE(X) CALCULATE 
C BESSEL FUNCTION   K(0,X) + I(0,X)*[ EILER+DLOG(X/2) ] 
C IF  X IS FROM INTERVAL [0,+8]
C  ( EILER=0.57721 56649 01532 86061 ).
C THROUGH EXPANSIONS OVER THE CHEBYSHEV POLYNOMS OF EVEN ORDER
C LITERATURE:Y.LUKE "SPECIAL MATHEMATICAL FUNCTIONS 
C                   AND ITS APPROXIMATIONS"
C                "MIR",MOSCOW  1980 608 P.
C TABLE  9.5  ; P.350
C COEFFICIENTS C
      DOUBLE PRECISION FUNCTION DBK0LE(X)    
      IMPLICIT NONE

      DOUBLE PRECISION X, W

      DOUBLE PRECISION QDCHEB0
      EXTERNAL QDCHEB0

      DOUBLE PRECISION B(19)             
      DATA B / 240.27705 96407 20389 10102D0, 
     1         369.47407 39728 67282 63764D0, 
     2         169.97341 16984 01148 04378D0, 
     3          49.02046 37772 63439 39371D0, 
     4           9.38849 73252 68442 32756D0, 
     5           1.25947 97636 67703 58618D0, 
     6           0.12377 69641 14924 54118D0, 
     7           0.00924 43098 62866 90621D0, 
     8           0.00054 06238 96492 55807D0, 
     9           0.00002 53737 96028 08704D0, 
     *           0.00000 09754 78302 83898D0, 
     1           0.00000 00312 49571 77932D0, 
     2           0.00000 00008 46434 70610D0, 
     3           0.00000 00000 19628 88451D0, 
     4           0.00000 00000 00393 96098D0, 
     5           0.00000 00000 00006 90835D0, 
     6           0.10673D-15, 0.146D-17, 0.2D-19  /
C                                                  
      W=X*0.125D0                                  
      DBK0LE=QDCHEB0(B,W,18,2)                 
      RETURN                                       
C                                                  
      END                                          



C FUNCTION  DBK1LO(X) CALCULATE 
C BESSEL FUNCTION 
C            1/X + I(1,X)*[EILER+LN(X/2)]- K(1,X)       
C  IF X IS FROM INTERVAL [0,+8]
C ( EILER=0.57721 56649 01532 86061 ).
C THROUGH EXPANSIONS OVER THE CHEBYSHEV POLYNOMS OF ODD  ORDER
C LITERATURE:Y.LUKE "SPECIAL MATHEMATICAL FUNCTIONS 
C                   AND ITS APPROXIMATIONS"
C                "MIR",MOSCOW  1980 608 P.
C TABLE  9.6  ; P.352
C COEFFICIENTS C 
      DOUBLE PRECISION FUNCTION DBK1LO(X)                      
      IMPLICIT NONE

      DOUBLE PRECISION X, W

      DOUBLE PRECISION QDCHEB0
      EXTERNAL QDCHEB0
      
      DOUBLE PRECISION C(18)                         
      DATA C / 418.88944 61663 96890 97522D0, 
     1         249.89554 90428 68080 38961D0, 
     2          91.18031 93387 41787 75763D0, 
     3          21.44499 50539 62240 43921D0, 
     4           3.43841 53928 80464 59793D0, 
     5           0.39484 60929 40938 23432D0, 
     6           0.03382 87455 26884 19281D0, 
     7           0.00223 57203 34170 88760D0, 
     8           0.00011 71310 22460 84561D0, 
     9           0.00000 49754 27122 13645D0, 
     *           0.00000 01746 04931 76984D0, 
     1           0.00000 00051 43294 11806D0, 
     2           0.00000 00001 28903 39664D0, 
     3           0.00000 00000 02780 94119D0, 
     4           0.00000 00000 00052 17097D0, 
     5           0.00000 00000 00000 85869D0, 
     6           0.00000 00000 00000 01250D0, 
     7           0.00000 00000 00000 00016D0 /
      W=X*0.125D0                             
      DBK1LO=QDCHEB0(C,W,17,1)            
      RETURN                                  
C                                             
      END                                     



C FUNCTION  DBK0GS(X) CALCULATE 
C BESSEL FUNCTION K(0,X)*SQRT(2*X/PI)*EXP(X) IF X>=5
C THROUGH EXPANSIONS OVER THE SHIFTED(DISPLACED) CHEBYSHEV POLYNOMS 
C LITERATURE:Y.LUKE "SPECIAL MATHEMATICAL FUNCTIONS 
C                   AND ITS APPROXIMATIONS"
C                "MIR",MOSCOW  1980 608 P.
C TABLE  9.5  ; P.351
C COEFFICIENTS D
      DOUBLE PRECISION FUNCTION DBK0GS(X)      
      IMPLICIT NONE

      DOUBLE PRECISION X, W

      DOUBLE PRECISION QDCHEB0
      EXTERNAL QDCHEB0

      DOUBLE PRECISION D(21)            
      DATA D / 0.98840 81742 30825 80035D0, 
     1        -0.01131 05046 46928 28069D0, 
     2         0.00026 95326 12762 72369D0, 
     3        -0.00001 11066 85196 66535D0, 
     4         0.00000 06325 75108 50049D0, 
     5        -0.00000 00450 47337 64110D0, 
     6         0.00000 00037 92996 45568D0, 
     7        -0.00000 00003 64547 17921D0, 
     8         0.00000 00000 39043 75576D0, 
     9        -0.00000 00000 04579 93622D0, 
     *         0.00000 00000 00580 81063D0, 
     1        -0.00000 00000 00078 83236D0, 
     2         0.00000 00000 00011 36042D0, 
     3        -0.00000 00000 00001 72697D0, 
     4  0.27545D-15,   -0.4589D-16,    0.796D-17,  -0.143D-17,
     5  0.27D-18,      -0.5D-19,       0.1D-19  /             
C
C
      W=5D0/X                                                 
      DBK0GS=QDCHEB0(D,W,20,3)                           
      RETURN                                                 
C
      END                                                    



C FUNCTION  DBK1GS(X) CALCULATE 
C BESSEL FUNCTION K(1,X)*SQRT(2*X/PI)*EXP(X) IF X>=5
C THROUGH EXPANSIONS OVER THE SHIFTED(DISPLACED) CHEBYSHEV POLYNOMS 
C LITERATURE:Y.LUKE "SPECIAL MATHEMATICAL FUNCTIONS 
C                   AND ITS APPROXIMATIONS"
C                "MIR",MOSCOW  1980 608 P.
C TABLE  9.6  ; P.351
C     COEFFICIENTS   E  
      DOUBLE PRECISION FUNCTION DBK1GS(X)                      
      IMPLICIT NONE

      DOUBLE PRECISION X, W

      DOUBLE PRECISION QDCHEB0
      EXTERNAL QDCHEB0

      DOUBLE PRECISION D(21)                         
      DATA D / 1.03595 08587 72358 33071D0,   
     1         0.03546 52912 43331 11380D0,   
     2        -0.00046 84750 28166 88856D0,   
     3         0.00001 61850 63810 05343D0,   
     4        -0.00000 08451 72048 12368D0,   
     5         0.00000 00571 32218 10284D0,   
     6        -0.00000 00046 45554 60661D0,   
     7         0.00000 00004 35417 33857D0,   
     8        -0.00000 00000 45757 29704D0,   
     9         0.00000 00000 05288 13281D0,   
     *        -0.00000 00000 00662 61293D0,   
     1         0.00000 00000 00089 04792D0,   
     2        -0.00000 00000 00012 72607D0,   
     3         0.00000 00000 00001 92086D0,   
     4 -0.30450D-15,    0.5046D-16,   -0.871D-17,   0.156D-17, 
     5 -0.29D-18,       0.6D-19,      -0.1D-19  / 
C
C
      W=5D0/X                                     
      DBK1GS=QDCHEB0(D,W,20,3)                
      RETURN                                      
C                                                 
      END                                         



C FUNCTION  DBI0LE(X) CALCULATE 
C BESSEL FUNCTION I(0,X) IF X IS FROM INTERVAL [-8,+8]
C THROUGH EXPANSIONS OVER THE CHEBYSHEV POLYNOMS OF EVEN ORDER
C LITERATURE:Y.LUKE "SPECIAL MATHEMATICAL FUNCTIONS 
C                   AND ITS APPROXIMATIONS"
C                "MIR",MOSCOW  1980 608 P.
C TABLE  9.5  ; P.350
      DOUBLE PRECISION FUNCTION DBI0LE(X)                      
      IMPLICIT NONE

      DOUBLE PRECISION X, W

      DOUBLE PRECISION QDCHEB0
      EXTERNAL QDCHEB0

      DOUBLE PRECISION A(19)                         
      DATA A / 127.73343 98121 81083 56301D0, 
     1         190.49432 01727 42844 19322D0, 
     2          82.48903 27440 24099 61321D0, 
     3          22.27481 92424 62230 87742D0, 
     4           4.01167 37601 79348 53351D0, 
     5           0.50949 33654 39982 87079D0, 
     6           0.04771 87487 98174 13524D0, 
     7           0.00341 63317 66012 34095D0, 
     8           0.00019 24693 59688 11366D0, 
     9           0.00000 87383 15496 62236D0, 
     *           0.00000 03260 91050 57896D0, 
     1           0.00000 00101 69726 72769D0, 
     2           0.00000 00002 68828 12895D0, 
     3           0.00000 00000 06096 89280D0, 
     4           0.00000 00000 00119 89083D0, 
     5           0.00000 00000 00002 06305D0, 
     6           0.3132D-16, 0.42D-18, 0.1D-19  /
C                                                
      W=X*0.125D0                                
      DBI0LE=QDCHEB0(A,W,18,2)               
      RETURN                                     
C                                                
      END                                        



C FUNCTION  DBI1LO(X) CALCULATE 
C BESSEL FUNCTION I(1,X) IF X IS FROM INTERVAL [-8,+8]
C THROUGH EXPANSIONS OVER THE CHEBYSHEV POLYNOMS OF ODD  ORDER
C LITERATURE:Y.LUKE "SPECIAL MATHEMATICAL FUNCTIONS 
C                   AND ITS APPROXIMATIONS"
C                "MIR",MOSCOW  1980 608 P.
C TABLE  9.6  ; P.352
      DOUBLE PRECISION FUNCTION DBI1LO(X)                      
      IMPLICIT NONE

      DOUBLE PRECISION X, W

      DOUBLE PRECISION QDCHEB0
      EXTERNAL QDCHEB0

      DOUBLE PRECISION A(18)                         
      DATA A / 220.60142 69235 23778 56112D0, 
     1         125.35426 68371 52356 46451D0, 
     2          42.86523 40931 28256 85130D0, 
     3           9.45300 52294 34910 53517D0, 
     4           1.42965 77090 76213 46814D0, 
     5           0.15592 42954 76256 29116D0, 
     6           0.01276 80490 81733 88545D0, 
     7           0.00081 08879 00690 69214D0, 
     8           0.00004 10104 61938 23750D0, 
     9           0.00000 16880 42203 43687D0, 
     *           0.00000 00575 86950 54206D0, 
     1           0.00000 00016 53453 53976D0, 
     2           0.00000 00000 40484 76606D0, 
     3           0.00000 00000 00854 96289D0, 
     4           0.00000 00000 00015 72708D0, 
     5           0.00000 00000 00000 25419D0, 
     6           0.00000 00000 00000 00364D0, 
     7           0.00000 00000 00000 00005D0 /
      W=X*0.125D0                             
      DBI1LO=QDCHEB0(A,W,17,1)  
      RETURN                                  
C                                             
      END                                     
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------
