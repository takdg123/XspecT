c Fortran subroutines for the ntee model. NONTH is called from xsnteea.cxx

c ------------------------------------------------------------------- c
 
      SUBROUTINE NONTH(Param,X,Jmax,Spnth,Sphth)
c  version: June 93
c  revised for reflection 1995
c
c
C  NONTHERMAL MODEL OF A COMPACT SOURCE
C  INITIAL POWER LAW ELECTRON INJECTION BETWEEN GMIN AND GMAX
C  COOLING ON MAGNETIC FIELD AND ON SYNCHROTRON PHOTONS
c  PAIR PRODUCTION  -> A STEADY STATE ELECTRON DISTRIBUTION
c
c Units: photon density  1/(sigma_T*R)
c        injection rate  c/(sigma_T*R^2)
c        gamma dot       c/R
c        radiated spectrum  dl/dln x
c
      IMPLICIT NONE

      INTEGER MAXENE
      PARAMETER (MAXENE=900)

      REAL besc , scale , ctaug , xplasm , delta0 , delta , xmin , b0 , 
     &     syncon , xmax ,  planck , 
     &     emagn , aa , cn , cndis , brcth , tautom , xl , sum , ga , 
     &     clpair , coulht , brcnth , gth , annlum , annht , thetaan , 
     &     dlt , xanlo , xanup , deltal , xnr , xr , taukn , arg , flz , 
     &     xs , xlth , taup , pterm
c     real functions
      REAL COMPTT1 , ANNSP
      REAL xm , ab0
      REAL Param(15) , g(MAXENE) , del(MAXENE) , qint0(MAXENE) , 
     &     qint(MAXENE) , taug(MAXENE) , dphth(MAXENE) , dwork(MAXENE) ,
     &     qintp(MAXENE) , eint(MAXENE), tauc(MAXENE) , gcoul(MAXENE) , 
     &     brnth(MAXENE) , gbrem(MAXENE) , dsyn(MAXENE) , brth(MAXENE) , 
     &     X(MAXENE) , DPH(MAXENE) , DPHesc(MAXENE) , DPHdot(MAXENE) , 
     &     TAUabs(MAXENE) , REL(MAXENE) , BET(MAXENE) , C1(MAXENE) , 
     &     C2(MAXENE) , Spnth(MAXENE) , Sphth(MAXENE)
      INTEGER jmin2(21) , jmax2(21)
      INTEGER ith , ith2 , n34 , i , j , jsyn , jmc22 , ju , j34 , 
     &        imaxth , iminth , ijmax , k1 , k2 , k3 , k , n , jlold , 
     &        jl , il , iu , kl , ku , janlo , janup , jnr , jrel , kk
      DOUBLE PRECISION q1 , q2 , q3 , q4 , q5 , q6 
      DOUBLE PRECISION qa , qb , qc , w , w1
      DOUBLE PRECISION z1 , z2 , z3 , z4 , z5 , z6
c  physical constants
      REAL A0 , SIGmat , VLIght , BCRit , XLAmbd , XME , PI , XRS
c  general parameters:
      REAL GMIn , GMAx , EXP0 , BFIeld , TEMpbb , OMEga , RADius , 
     &     XLEl , XLBb , GAMma0 , GT , theta
      INTEGER IMIns , IMAx , Jmax , JMC2 , ITEr
c
      COMMON /PHYSNT/ A0 , SIGmat , VLIght , BCRit , XLAmbd , XME , PI , 
     &                XRS
      COMMON /PARNT / GMIn , GMAx , EXP0 , BFIeld , TEMpbb , OMEga , 
     &                RADius , XLEl , XLBb , GAMma0 , GT
      COMMON /INTNT / IMIns , IMAx , JMC2 , ITEr
      COMMON /NTEE_THERM/ DPH , DPHesc , DPHdot , TAUabs , REL , BET , 
     &                    C1 , C2
c
C  INPUT:
c   the minimum and maximum Lorentz factor for the electron injection
c   gmin used only for power law injection
      GMIn = Param(9)
      GMAx = Param(5)
C  EXP0 - INDEX OF POWER LAW INJECTION DISTRIBUTION
C   EXP0=0 - DELTA-FUNCTION ELECTRON INJECTION
      EXP0 = Param(8)
C  xlel - DIMENSIONLESS LUMINOSITY IN INJECTED ELECTRON SPECTRUM
      XLEl = Param(1)
c  xlbb - dimensionless blackbody luminosity
      XLBb = Param(2)
c     the thermal compactness:
c     For pure nonthermal plasmas, xlth=0
      xlth = Param(6)
C  TEMPBB - TEMPERATURE OF EXTERNAL PHOTONS IN MC**2
      TEMpbb = Param(4)/5.11E5
c    taup is the optical depth of ionization electrons
      taup = Param(7)
c  besc - the beta escape parameters of pairs, defined as in Zdziarski (1985)
      besc = Param(12)
C  RADIUS - RADIUS OF THE SOURCE in cm
      RADius = Param(11)
C  GAMMA0 - MINIMUM LORENTZ FACTOR (E.G., =1.1) used for nonthermal
c  reprocessing; gmin.geq.gamma0; gamma0>1+3*Theta
      GAMma0 = Param(10)
c     the scaling factor for the reflected spectrum;
c     1 corresponds to seeing equal
c     contributions from the reflected and direct spectra
      scale = Param(3)
c  inclination angle
      xm = Param(13)
c  iron abundance
      ab0 = ALOG10(Param(14))
c
C  ITER - NUMBER OF ITERATIONS
c  ith  - number of iterations with 'thermal'
c  ith2 - number of calls to 'thermal' within one iteration
      ITEr = 16
      ith = 11
      ith2 = 4
C  BFIELD - MAGNETIC FIELD STRENGTH
C     IF(BFIELD.LE.0.) BFIELD IS IN EQUIPARTITION WITH the TOTAL PHOTON ENERGY
      BFIeld = 0
C  GT - VALUE OF THE LORENTZ FACTOR CORRESPONDING TO THE SELF-ABSORPTION
C    FREQUENCY  -- if bfield=0, no self-absorption
      GT = 50
 
c  clear arrays (important for repeated calls)
 
      DO 100 j = 1 , MAXENE
         dphth(j) = 0
         dsyn(j) = 0
         taug(j) = 0
         tauc(j) = 0
         g(j) = 0
         del(j) = 0
         qint0(j) = 0
         qint(j) = 0
         dwork(j) = 0
         qintp(j) = 0
         eint(j) = 0
         gcoul(j) = 0
         brnth(j) = 0
         gbrem(j) = 0
         brth(j) = 0
         X(j) = 0
         DPH(j) = 0
         DPHesc(j) = 0
         DPHdot(j) = 0
         TAUabs(j) = 0
         REL(j) = 0
         BET(j) = 0
         C1(j) = 0
         C2(j) = 0
 100  CONTINUE
 
C
C  PHYSICAL CONSTANTS:
C
      BCRit = 4.414E13
      VLIght = 3E10
      SIGmat = 6.65E-25
      XLAmbd = 2.426E-10
      XME = 9.1095E-28
      A0 = 5.29E-9
      PI = 3.14159
c     dimensionless volume
      XRS = 4*PI/3
c     pair optical depth constant (Svensson 1987)
      ctaug = 0.2
c     dimensionless plasma frequency for unit Thomson depth,
c     should be multiplied by sqrt(tautom):
      xplasm = SQRT(3*A0/2/RADius)
c     n34 - number of bins between 3/4 and 1
      n34 = 4
c
c
c
C  ARRAY FOR LORENTZ FACTORS:
C
      delta0 = LOG10(4./3)/n34
      IMAx = INT(LOG10(GMAx/GAMma0)/delta0 + 1)
c     delta .ge. delta0
c     delta is the 10-log interval between the energy points
      delta = LOG10(GMAx/GAMma0)/(IMAx-1)
c  deltal - the e-log interval between the energy points
      deltal = delta*LOG(10.)
      DO 200 i = 1 , IMAx
         g(i) = 10.**((i-1)*delta)*GAMma0
 200  CONTINUE
      IF ( ABS(EXP0).GT.1.E-15 ) IMAx = IMAx - 1
C
C SYNCHROTRON -----
C
C       DIMENSIONLESS CYCLOTRON FREQUENCY:
      b0 = BFIeld/BCRit
c       constant for synchrotron radiation (depending on radius)
      syncon = RADius/A0/8
 
      IF ( ABS(BFIeld).LT.1.E-15 ) THEN
         GT = GAMma0
         IMIns = 0
c       the minimum blackbody energy:
         xmin = 8.E-2*TEMpbb
      ELSE
C       THE EQUIPARTITION VALUE OF BFIELD - with synchrotron photons
c       only - assuming Lsyn=Lcom=1/2 Ltotal
c       the mean photon time inside the source =R/c (no 3/4, taut<1)
         IF ( BFIeld.LT.0. ) BFIeld = VLIght*SQRT(3*XLEl*XME/SIGmat/
     &                                RADius)
C       del(imins) corresponds to the synchrotron turnover frequency
         IMIns = INT(LOG10(GT/GAMma0)/delta + 1)
         xmin = 1.3334*b0*g(IMIns+1)**2
      ENDIF
 
c       maximum synchrotron energy (reduced since 'synchr' interpolates
c       forward):
      jsyn = 2*(IMAx-IMIns) - 1
C
C -----------------
 
C
C  JMAX - # OF PHOTON ENERGIES
C  X(JMAX) < GMAX-1
C
      Jmax = INT(LOG10((g(IMAx)-1)/xmin)/delta + 1)
      IF ( Jmax .GT. MAXENE ) THEN
         CALL xwrite(
     &      'ntee: too many internal energy bins: increase MAXENE', 5)
      ENDIF
C  X(JMC2) = ME C**2
      JMC2 = INT(LOG10(1/xmin)/delta + 1)
C
C  X - ARRAY FOR PHOTON ENERGIES
C
      DO 300 j = 1 , Jmax + 1
         X(j) = 10.**((j-JMC2)*delta)
 300  CONTINUE
      xmin = X(1)
      xmax = X(Jmax)
c     used in calls to comptt
      jmc22 = INT(LOG10(0.5/xmin)/delta + 2)
c     minimum upper limit in the tauc integral
      ju = INT(LOG10(0.75/g(IMAx)/xmin)/delta + 1)
c     x(j34)=3/4/gamma0
      j34 = INT(JMC2 + LOG10(0.75/GAMma0)/delta)
c
c compute c1(x), c2(x), and rel(x) arrays
c c1(x) and c2(x) are the relativistic corrections to K-equation
c rel(x) is Klein-Nishina cross section divided by thomson crossection
      DO 400 j = 1 , Jmax
         w = X(j)
c  c1, c2 is the Cooper's coefficient calculated at w and w1, resp.
         C1(j) = SNGL(w**4/(1+4.6*w+1.1*w*w))
c   the commented lines are for 0-Theta Fokker-Planck eq. (Spring 1987)
c   use the file coeff.f
c d1(x)=c1(x), d2 is derivative of d1, used to determine theta
c       d1(j)=coeff1(w)
c   coeff3 is derivative of x**4*coeff1
c       if (w.gt.5.e-3) then
c         d2(j)=coeff3(w)/w**4-d1(j)*4/w
c       else
c         d2(j)=-161./5*w+2304./7*w*w
c       end if
         w1 = SQRT(X(j)*X(j+1))
         C2(j) = SNGL(w1**4/(1+4.6*w1+1.1*w1*w1))
c       w1 is x(j+1/2) (x(i) defined up to jmax+1)
c       c1(j)=coeff1(w1)
c       c2(j)=coeff2(w1)
         IF ( w.LE.0.05 ) THEN
c         use asymptotic limit for rel(x) for x less than 0.05
            REL(j) = SNGL(1. - 2.*w + 26.*w*w/5.)
         ELSE
            z1 = (1.+w)/w**3
            z2 = 1. + 2.*w
            z3 = LOG(z2)
            z4 = 2.*w*(1.+w)/z2
            z5 = z3/2./w
            z6 = (1.+3.*w)/z2/z2
            REL(j) = SNGL(0.75*(z1*(z4-z3)+z5-z6))
         ENDIF
 400  CONTINUE
 
c  the thermal emision spectrum 'dphth'
      IF ( ABS(XLBb).GT.1.E-15 ) THEN
c       blackbody dilution factor
         OMEga = 45*XLAmbd**3/32/PI/PI/SIGmat/RADius/(PI*TEMpbb)**4*XLBb
         imaxth = INT(LOG10(50*TEMpbb/xmin)/delta)
         iminth = INT(MAX(LOG10(8.E-2*TEMpbb/xmin)/delta+1,1.))
         planck = 15/(PI*TEMpbb)**4/XRS*XLBb
         DO 450 j = iminth , imaxth
            dphth(j) = planck*X(j)**2/(EXP(X(j)/TEMpbb)-1)
 450     CONTINUE
      ENDIF
C
C  2*G(IJMAX) < X(JMAX) -
C   THE CONDITION FOR THE MAXIMUM ENERGY OF THE PRODUCED PAIRS
C
      ijmax = INT(LOG10(xmax/2/GAMma0)/delta + 1)
C
C  PHOTON ENERGY RANGES .4-2 AND 2-10 keV
C  (TO COMPUTE LUMINOSITIES & SPECTRAL INDICES IN THAT REGIONS)
C
      k1 = INT(LOG10(7.83E-4/xmin)/delta + 2)
      k2 = INT(LOG10(3.91E-3/xmin)/delta + 1)
      k3 = INT(LOG10(1.96E-2/xmin)/delta + 1)
C
C  X(JMIN2(K)) - MINIMUM PHOTON ENERGY IN THE K-TH SCATTERING
C  X(JMAX2(K)) - MAXIMUM PHOTON ENERGY IN THE K-TH SCATTERING
C
C  X(JMAX2(1)) < 1.333*GMAX**2*X0
C
      jmin2(1) = 1
c  the maximum energy in the 0-order photons:
c     I assume that the non-zero synchrotron spectrum always extends to a
c     larger energy then the thermal spectrum
      IF ( ABS(BFIeld).GT.1.E-15 ) THEN
         jmax2(1) = imaxth
      ELSE
         jmax2(1) = IMAx
      ENDIF
      DO 500 k = 2 , ITEr + 1
         jmin2(k) = MIN(Jmax+1,INT(LOG10(4./3*GAMma0**2*X(jmin2(k-1))/
     &              xmin)/delta+2))
         jmax2(k) = MIN(Jmax,INT(LOG10(4./3*g(IMAx)**2*X(jmax2(k-1))/
     &              xmin)/delta+1))
 500  CONTINUE
C
c     set up the initial (nonzero) theta and tautom
      theta = 2.E-2
      tautom = MAX(1.E-4,taup)
c
c  determine the initial electron distribution constants:
c     constants:
      CALL ZEROCN(emagn,aa,cn,cndis)
C
c     the electron injection and initial electron distribution:
      CALL ELDIS(g,aa,cn,cndis,qint0,del)
c     thermal bremsstrahlung emission
      CALL THBREM(X,brth,brcth,theta,tautom,taup,deltal,Jmax)
c     the synchrotron component
      IF ( ABS(BFIeld).GT.1.E-15 )
     &     CALL SYNCHR(dsyn,del,g,X,jsyn,b0,syncon,deltal)
 
c     the initial photon spectrum:
      DO 600 j = 1 , Jmax
         DPH(j) = dphth(j) + dsyn(j)
 600  CONTINUE
c
C*****  iteration loop
c      ipl is the index used to call 'table', to print out results
c       ipl=1
      DO 1900 n = 1 , ITEr
C
C       PAIR INJECTION DISTRIBUTION:
         DO 650 j = JMC2 , Jmax
            dwork(j) = DPH(j)*taug(j)*X(j)*deltal
 650     CONTINUE
c       integral over dwork; see Lighman and Zdziarski (1987)
c       ijmax assumed ge. 2
c       qintp is equal 2 times the integral, so no /2
         qintp(ijmax) = dwork(Jmax)
         qintp(ijmax-1) = dwork(Jmax) + dwork(Jmax-1)
         jlold = Jmax - 1
c       trapezoidal integration; jlold is previous lower limit;
c       loop from high to low values of g(i)
         DO 700 i = ijmax - 2 , 1 , -1
c         the lower limit
            xl = g(i) + SQRT(g(i)*g(i)-1)
            jl = INT(LOG10(xl/xmin)/delta + 2)
            sum = qintp(i+1)
            DO 660 j = jl , jlold - 1
               sum = sum + dwork(j) + dwork(j+1)
 660        CONTINUE
            qintp(i) = sum
            jlold = jl
 700     CONTINUE
c       the sum of electron and pair integrals
         DO 750 i = 1 , IMAx
            qint(i) = qintp(i) + qint0(i)
 750     CONTINUE
C
C       EQUILIBRIUM THOMSON THICKNESS AND PAIR LUMINOSITY:
c   relativistic correction
         ga = 1/(1+2*theta**2/LOG(1.12*theta+1.3))
c        Thomson depth: acceleration has zero net effect; only secondary
c        pairs contribute; ionization electrons with the optical depth taup
c        included
         pterm = qintp(1)
c         assume tautom larger than zero, to calculate 'bet'
         IF ( pterm.NE.0 ) tautom = (0.375*ga*taup**2+2*(besc*taup+pterm
     &                              ))
     &                              /(besc+SQRT(besc**2+(0.375*ga*taup)
     &                              **2+0.75*ga*(besc*taup+pterm)))
c       pair luminosity (in the pair rest mass)
         clpair = pterm*XRS
c
c
c       the Coulomb and bremsstrahlung cooling rates
         DO 800 i = 1 , IMAx
            gcoul(i) = 0.75*tautom*LOG(g(i)/xplasm**2/tautom)
            gbrem(i) = MAX(0.,3.485E-3*tautom*(LOG(2*(g(i)-0.999))-1/3.)
     &                 )
 800     CONTINUE
C
C       COMPTON COOLING INTEGRAL
         DO 850 j = 1 , JMC2
            dwork(j) = DPH(j)*X(j)**2*deltal
 850     CONTINUE
c       calculate first integral from xmin to 3/4/g(imax)
         sum = (dwork(1)+dwork(ju))/2
         DO 900 j = 2 , ju - 1
            sum = sum + dwork(j)
 900     CONTINUE
         eint(IMAx) = sum
c       add imax-1 further terms to calculate eint(i)
         DO 950 i = IMAx - 1 , 1 , -1
            eint(i) = eint(i+1) + (dwork(ju+IMAx-i-1)+dwork(ju+IMAx-i))
     &                /2
 950     CONTINUE
C       CALCULATE THE ELECTRON DISTRIBUTION
         DO 1000 i = 1 , IMIns
            del(i) = qint(i)
     &               /(eint(i)*(1.333*g(i)**2-1)+gcoul(i)+gbrem(i))
 1000    CONTINUE
         DO 1050 i = IMIns + 1 , IMAx
            del(i) = qint(i)/(eint(i)*(1.333*g(i)**2-1)+emagn*1.333*g(i)
     &               **2+gcoul(i)+gbrem(i))
 1050    CONTINUE
c
c       calculate the Coulomb heating rate
c       (the integral over gcoul(i)*del(i) )
         sum = (gcoul(1)*del(1)*g(1)+gcoul(IMAx)*del(IMAx)*g(IMAx))/2
         DO 1100 i = 2 , IMAx - 1
            sum = sum + gcoul(i)*del(i)*g(i)
 1100    CONTINUE
         coulht = sum*deltal
C
C       THE NEW LORENTZ FACTOR CORRESPONDING TO THE TURNOVER FREQUENCY
c        if(imins.ge.2) then
c          pindex=-(log(del(imins))-log(del(imins-1)))/deltal
c         THE NEW LORENTZ FACTOR CORRESPONDING TO THE TURNOVER FREQUENCY
c          gtnew=gt0(del(imins)*g(imins)**pindex,pindex,radius,b0)
c        end if
c
c       Optical depth to Compton scatterings with relativistic electrons --
c       integral over the electron distribution
         DO 1150 i = 1 , IMAx
            dwork(i) = del(i)*g(i)*deltal
 1150    CONTINUE
         tauc(j34) = dwork(1)/2
         tauc(j34-1) = (dwork(1)+dwork(2))/2
c       continuous continuation up to 1:
         DO 1200 j = j34 + 1 , JMC2 - 1
            tauc(j) = tauc(j34)*LOG(X(j))/LOG(X(j34))
 1200    CONTINUE
c       trapezoidal integral
         DO 1250 j = j34 - 2 , j34 - IMAx + 1 , -1
            tauc(j) = tauc(j+1) + (dwork(j34-j)+dwork(j34-j+1))/2
 1250    CONTINUE
c       constant tauc:
         DO 1300 j = 1 , j34 - IMAx
            tauc(j) = tauc(j34-IMAx+1)
 1300    CONTINUE
c       synchrotron spectrum:
         IF ( ABS(BFIeld).GT.1.E-15 )
     &        CALL SYNCHR(dsyn,del,g,X,jsyn,b0,syncon,deltal)
c       nonthermal and thermal bremsstrahlung spectra
         CALL NTHBREM(brnth,brcnth,g,del,X,deltal,IMAx,Jmax,tauc(1),
     &                tautom)
         CALL THBREM(X,brth,brcth,theta,tautom,taup,deltal,Jmax)
c       annihilation heating and cooling
c***    IMPORTANT!
c**     for external pair injection, replace pterm by qint(1)
         CALL ANNIH(theta,tautom,taup,pterm,gth,annlum,annht)
c
c       compute Compton photon production rate
         DO 1350 j = jmin2(2) , Jmax
c         for each x, we will integrate over gamma
c         limits of integration:
            il = MAX(1,INT(LOG10(X(j)/GAMma0)/delta+1))
            iu = MIN(IMAx,INT(LOG10(0.75*X(j)/xmin/GAMma0**2)/2/delta+1)
     &           )
c         the boundary elements of the integral (kl.ge.ku)
            kl = INT(LOG10(0.75*X(j)/g(il)**2/xmin)/delta + 1)
            ku = INT(LOG10(0.75*X(j)/g(iu)**2/xmin)/delta + 1)
            sum = (DPH(kl)*del(il)/g(il)+DPH(ku)*del(iu)/g(iu))/2
            DO 1320 i = il + 1 , iu - 1
               k = INT(LOG10(0.75*X(j)/g(i)**2/xmin)/delta + 1)
               sum = sum + DPH(k)*del(i)/g(i)
 1320       CONTINUE
c         add contribution to dphdot from blackbody
c         synchrotron, and bremsstrahlung photons
            DPHdot(j) = 0.75*sum*deltal + dphth(j) + dsyn(j) + brnth(j)
     &                  + brth(j)
 1350    CONTINUE
c       define dphdot for lowest x's
         DO 1400 j = 1 , jmin2(2) - 1
            DPHdot(j) = dphth(j) + dsyn(j) + brnth(j) + brth(j)
 1400    CONTINUE
c
c       for n=1, tabulate the photon production rate w/o any other effects
c        if (n.eq.1) call table(dphdot,x,jmax,del,g,ipl,xrs)
c
c    Include in dphdot those annihilation photons
c    that get either scattered or absorbed
c    The annihilation spectrum is based on temperature from the previous
c    iteration
 
c  Correct dphdot for the case when theta is so low that the line width
c  becomes less than the bin width
         IF ( SQRT(PI*theta).LT.deltal ) THEN
            thetaan = deltal**2/PI
         ELSE
            thetaan = theta
         ENDIF
 
c     Maximum and minimum energies of annihilation photons
c     (see Svensson 1983; exp(-10) maximum decrease assumed)
         dlt = SQRT((10*thetaan)**2+40*thetaan)
         xanlo = 1 + 5*thetaan - dlt/2
         xanup = 1 + 5*thetaan + dlt/2
         janlo = INT(LOG(xanlo/xmin)/deltal + 2)
         janup = INT(LOG(xanup/xmin)/deltal + 1)
         DO 1450 j = janlo , janup
            DPHdot(j) = DPHdot(j) + ANNSP(X(j),thetaan,tautom,taup)
     &                  *(TAUabs(j)+tautom*REL(j))
     &                  /(1+tautom*REL(j)+TAUabs(j))
 1450    CONTINUE
c
c compute beta array, the probalility of escape per Thomson time
c bet evaluated for spherical geometry and nearly uniform sources.
c Between x=0.1 and 1.0, a function flz modifies beta to allow
c the increasingly large energy change per scattering to gradually
c eliminate spatial diffusion
         jnr = INT(JMC2 + LOG10(0.1)/delta)
         jrel = INT(JMC2 + LOG10(1.)/delta)
         xnr = X(jnr)
         xr = X(jrel)
         DO 1500 j = 1 , jnr - 1
            taukn = tautom*REL(j)
            BET(j) = 1./tautom/(1.+taukn/3.)
 1500    CONTINUE
         DO 1550 j = jnr , jrel
            taukn = tautom*REL(j)
            arg = (X(j)-xnr)/(xr-xnr)
            flz = 1. - arg
            BET(j) = 1./tautom/(1.+taukn/3.*flz)
 1550    CONTINUE
         DO 1600 j = jrel + 1 , Jmax
            BET(j) = 1./tautom
 1600    CONTINUE
c
c skip over solution of dph and taug if thermal already called
         IF ( tautom.LT.0.3 .OR. n.LE.ITEr-ith+1 ) THEN
c
c       compute dph and taug self-consistently
c
            DO 1620 j = 2*JMC2 - Jmax , Jmax
               q1 = DPHdot(j)
               q2 = BET(j)*tautom + tauc(j)
               q3 = ctaug/X(j)
               k = 2*JMC2 - j
               q4 = DPHdot(k)
               q5 = BET(k)*tautom + tauc(k)
               q6 = ctaug/X(k)
               qa = q6*q2
               qb = q2*q5 + q3*q4 - q1*q6
               qc = -q1*q5
c       quadratic equations
               IF ( qb.GE.0. ) THEN
                  DPH(j) = SNGL(-2*qc/(qb+SQRT(qb*qb-4*qa*qc)))
               ELSE
                  DPH(j) = SNGL((-qb+SQRT(qb*qb-4*qa*qc))/2/qa)
               ENDIF
 1620       CONTINUE
            DO 1640 j = 2*JMC2 - Jmax , Jmax
               k = 2*JMC2 - j
               taug(j) = ctaug*DPH(k)/X(j)
 1640       CONTINUE
c       define tauabs(x)
            DO 1660 j = 1 , Jmax
               TAUabs(j) = taug(j) + tauc(j)
 1660       CONTINUE
c       define dph(x) in the region not covered above
            DO 1680 j = 1 , 2*JMC2 - Jmax - 1
               DPH(j) = DPHdot(j)/(BET(j)*tautom+TAUabs(j))
 1680       CONTINUE
 
c  compute theta (x(jmc22)=0.5)
c  XLTH is the direct thermal heating in addition to the radiative and
c  Coulomb ones. It appears in the command below and in the call to THERML
c  a few lines below.
            theta = COMPTT1(DPH,jmc22,Jmax,X,deltal,REL,
     &              (xlth/XRS+coulht+annht-brcth)/tautom,theta)
c compute escaping photon density
            DO 1700 j = 1 , Jmax
               DPHesc(j) = DPH(j)*BET(j)*tautom
 1700       CONTINUE
         ENDIF
 
         IF ( n.GT.ITEr-ith .AND. tautom.GE.0.3 ) THEN
            DO 1720 kk = 1 , ith2
               CALL THERML(tautom,theta,deltal,XRS,xlth/XRS+coulht+
     &                     annht-brcth,X,Jmax)
c     tauabs for small indices not affected
               DO 1710 j = 2*JMC2 - Jmax , Jmax
                  k = 2*JMC2 - j
                  taug(j) = ctaug*DPH(k)/X(j)
                  TAUabs(j) = taug(j) + tauc(j)
 1710          CONTINUE
 1720       CONTINUE
         ENDIF
 
c compute xstar
         DO 1750 j = 1 , Jmax
            IF ( taug(j).GT.1. ) GOTO 1800
 1750    CONTINUE
 1800    xs = 0.5*(X(j)+X(j-1))
 
c  add freely escaping annihilation photons
         DO 1850 j = janlo , janup
            DPHesc(j) = DPHesc(j) + ANNSP(X(j),theta,tautom,taup)
     &                  /(1+tautom*REL(j)+TAUabs(j))
 1850    CONTINUE
c
c       store the distributions
c        if(n.eq.iter.or.n.eq.iter-ith) then
c          ipl=ipl+1
c          call table(dphesc,x,jmax,del,g,ipl,xrs)
c        end if
c
c     end of the main loop
c
 1900 CONTINUE
c
c     the nonthermal component in E F_E
      DO 2000 j = 1 , Jmax
         Spnth(j) = DPHesc(j)*X(j)**2*XRS
 2000 CONTINUE

c     the thermal component in E F_E

      DO j = 1 , Jmax
         Sphth(j) = Dphth(j)*X(j)**2
      ENDDO

      RETURN
c
C****  WRITE MODEL PARAMETERS:
C
c      write(6,41)GMIN,GMAX,EXP0,RADIUS,GAMMA0,XLEL,XLBB,
c     >  imins,IMAX,JMAX,IJMAX,jmc2,ju,j34,K1,K2,K3,n34,ith,ith2,
c     >  CN,TEMPBB,OMEGA
c99001 FORMAT (//,11X,'PARAMETERS:',/,' GMIN = ',1PE9.2,' GMAX = ',E9.2,
c     &        ' EXP0 = ',0PF7.3,/,' RADIUS [CM] = ',1PE9.2,' GAMMA0 = ',
c     &        E9.2,' XLEL = ',E9.2,/,' XLBB=',E9.2,' imins=',i3,
c     &        ' IMAX= ',I3,' JMAX= ',I3,/,' IJMAX= ',I3,' jmc2=',i3,
c     &        ' ju=',i3,' j34=',i3,' K1= ',I3,' K2= ',I3,' K3= ',I3,/,
c     &        ' n34=',i3,' ith=',i3,' ith2=',i3,' CN = ',E9.2,/,
c     &        ' TEMPBB = ',E9.2,' OMEGA = ',E9.2)
c
c       write some output
c        write(6,60) n
c99002 FORMAT (/,1x,'iteration # ',i3)
c       write the results for electrons:
C
c        write(6,42)tautom,clpair
c99003 FORMAT (/,' tau_T=',1pe9.2,' L(pair rest mass)=',e9.2)
c       write the bremsstrahlung and Coulomb luminosities
c        write(6,66) brcth*xrs,brcnth*xrs,coulht*xrs
c99004 FORMAT (1x,'thbr lum.=',1pe10.3,' nthbr lum.=',e10.3,
c     &        ' Coul. lum.=',e10.3)
c        write(6,67) annlum*xrs,annht*xrs,gth
c99005 FORMAT (1x,'ann. lum.=',1pe9.2,' ann. heat.=',e9.2,' gth=',e9.2)
 
 
c       print the results for photons
c        write(6,777) xs, tauc(1), theta
c99006 FORMAT (' E(tau_gg=1) = ',1pe9.2,' tau_nth = ',e9.2,' theta = ',
c     &        e9.2)
      END
**==therml.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
      SUBROUTINE THERML(Tautom,Theta,Deltal,Xrs,Extht,X,Jmax)
c this program computes the effects of Comptonization by
c nonrelativistic thermal electrons in a sphere of radius R,
c including absorption, escape, and
c relativistic corrections up to photon energies of 1 MeV.
c the dimensionless photon energy is x=hv/(m*c*c)
c
c the input parameters and functions are:
c dphdot(x), the photon production rate (photons/volume/time/energy)
c tautom, the Thomson scattering depth
c tauabs(x), the absorption depth
c theta, the initial value of temperature in units of m*c*c
c c1(x), c2(x), and bet(x), the coefficients in the K-equation and the
c   probability of photon escape per Thomson time, respectively,
c   including Klein-Nishina corrections
c rel(x) is the ratio of k-n crossection to thomson crossection
c the output parameters and functions are:
c dph(x), the interior photon density, in units of R/c
c dphesc(x), the escaping photon density, in units of R/c
c theta, the equilibrium electron temperature
c exheat is the rate of nonradiative (e.g., Coulomb) heating.
c If exheat=0, the equilibrium temperature is the Compton temperature.
 
      IMPLICIT NONE

      INTEGER MAXENE
      PARAMETER (MAXENE=900)

      REAL aa , BET , C1 , C2 , c20 , CLUM , COMPTT1 , Deltal , DPH , 
     &     DPHdot , DPHesc , dy , ein , eout , eta , Extht , REL , t1 , 
     &     t2 , t3
      REAL tau , TAUabs , Tautom , Theta , w , w1 , w2 , x12 , x23 , 
     &     x32 , Xrs , ynew , yp
      INTEGER IMAx , IMIns , ITEr , j , jj , jj1 , Jmax , JMC2 , jx1 , 
     &        jx2
      COMMON /NTEE_THERM/ DPH(MAXENE) , DPHesc(MAXENE) , DPHdot(MAXENE),
     &                    TAUabs(MAXENE) , REL(MAXENE) , BET(MAXENE) ,
     &                    C1(MAXENE) , C2(MAXENE)
      COMMON /INTNT / IMIns , IMAx , JMC2 , ITEr
      REAL a(MAXENE) , b(MAXENE) , c(MAXENE) , f1(MAXENE) , X(MAXENE)
      REAL d(MAXENE) , alp(MAXENE) , g(MAXENE) , gam(MAXENE) , u(MAXENE)
c u(x) is the dimensionless photon occupation number
c jmax is the maximum electron photon energy
c jmc2 is defined by x(jmc2)=1.0
 
      c20 = Tautom/Deltal
      x12 = 1.2
c x12 is the value of x separating region 1 from
c region 2. For larger x, photons are simply absorbed in
c relativistic scatterings and the K-equation is not used
      jx1 = INT(JMC2 + LOG(x12)/Deltal + 1)
c x23 separates regions where the first and the second order equations are
c used; it is not used when a 0-theta dispersion term is added to temperature
      x23 = MIN(1.2,0.2+80*Theta)
      jx2 = INT(JMC2 + LOG(x23)/Deltal + 1)
c compute input energy, 'ein'
      DO 100 j = 1 , Jmax
         f1(j) = (DPHdot(j)-DPH(j)*TAUabs(j))*X(j)**2
 100  CONTINUE
      ein = CLUM(f1,1,Jmax,Deltal*Xrs)
 
c determine u in region 2, j.ge.jx1
      DO 200 j = jx1 , Jmax
         tau = TAUabs(j) + Tautom*REL(j)
         u(j) = DPHdot(j)/X(j)/X(j)/(tau+BET(j)*Tautom)
 200  CONTINUE
 
c  the first-order equation:
c  the boundary value is at the point next to x(jx1)
c  set jx2=jx1 when removing this region
      DO 300 j = 1 , jx1 - jx2
         jj = jx1 - j
         jj1 = jj + 1
         w = X(jj1)
         eta = -w*DPHdot(jj1)
     &         /Tautom + w**3*(TAUabs(jj1)+BET(jj1)*Tautom)*u(jj1)
     &         /Tautom
         dy = -eta*Deltal
         yp = u(jj1)*C1(jj1) + dy
         ynew = C1(jj1)*u(jj1)
     &          + 0.5*(dy-Deltal*(-X(jj)*DPHdot(jj)/Tautom+X(jj)
     &          **3*(TAUabs(jj)+BET(jj)*Tautom)*yp/C1(jj)/Tautom))
         u(jj) = ynew/C1(jj)
 300  CONTINUE
c determine u in region 1, 1.ge.j.lt.jx2
c define coefficients going into equation
c a(j)*u(j+1)+b(j)*u(j)+c(j)*u(j-1)=d(j)
      DO 400 j = 2 , jx2 - 1
         w1 = SQRT(X(j)*X(j+1))
         w2 = SQRT(X(j-1)*X(j))
c  w1 is x(j+1/2)
c  w2 is x(j-1/2)
         a(j) = -c20*C2(j)*(Theta/Deltal/w1+0.5)
         t1 = -c20*C2(j)*(0.5-Theta/Deltal/w1)
         t2 = c20*C2(j-1)*(Theta/Deltal/w2+0.5)
         t3 = X(j)**3*(TAUabs(j)+Tautom*BET(j))
         b(j) = t1 + t2 + t3
         c(j) = c20*C2(j-1)*(0.5-Theta/Deltal/w2)
         d(j) = X(j)*DPHdot(j)
c  the 0-theta F-P eq.
c       a(j)=-c20*w1**4*((c2(j)+theta)/deltal/w1+0.5*c1(j))
c       t1=-c20*w1**4*(0.5*c1(j)-(c2(j)+theta)/deltal/w1)
c       t2=c20*w2**4*((c2(j-1)+theta)/deltal/w2+0.5*c1(j-1))
c       c(j)=c20*w2**4*(0.5*c1(j-1)-(c2(j-1)+theta)/deltal/w2)
 400  CONTINUE
 
c define constants going into boundary terms
c  u(1)=aa*u(2) (zero flux at lowest energy)
c  u(jx2) given from region 2 above
      x32 = SQRT(X(1)*X(2))
      aa = (Theta/Deltal/x32+0.5)/(Theta/Deltal/x32-0.5)
 
c invert tridiagonal matrix
      alp(2) = b(2) + c(2)*aa
      gam(2) = a(2)/alp(2)
      DO 500 j = 3 , jx2 - 1
         alp(j) = b(j) - c(j)*gam(j-1)
         gam(j) = a(j)/alp(j)
 500  CONTINUE
      g(2) = d(2)/alp(2)
      DO 600 j = 3 , jx2 - 2
         g(j) = (d(j)-c(j)*g(j-1))/alp(j)
 600  CONTINUE
      g(jx2-1) = (d(jx2-1)-a(jx2-1)*u(jx2)-c(jx2-1)*g(jx2-2))/alp(jx2-1)
      u(jx2-1) = g(jx2-1)
      DO 700 j = 3 , jx2 - 1
         jj = jx2 + 1 - j
         u(jj) = g(jj) - gam(jj)*u(jj+1)
 700  CONTINUE
      u(1) = aa*u(2)
 
c compute new value of dph(x) and new value of dphesc(x)
      DO 800 j = 1 , Jmax
         DPH(j) = X(j)*X(j)*u(j)
         DPHesc(j) = DPH(j)*BET(j)*Tautom
         f1(j) = DPHesc(j)*X(j)**2
 800  CONTINUE
c  compute output energy
      eout = CLUM(f1,1,Jmax,Deltal*Xrs)
c  the final output energy is larger by the energy of freely escaping
c  annihilation photons
 
c determine revised temperature based on new value of u
      Theta = COMPTT1(DPH,jx1,Jmax,X,Deltal,REL,Extht/Tautom,Theta)
      RETURN
c      write(6,29) theta, ein, eout
c99001 FORMAT (10x,'Theta= ',1pe9.2,5x,'P_in=',e9.2,5x,'P_out=',e9.2)
      END
**==comptt.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
      FUNCTION COMPTT1(Dph,Jx1,Jmax,X,Deltal,Rel,Z8,Theta)
c     calculates the Compton temperature based on the internal photon
c     density dph
      IMPLICIT NONE
      INTEGER MAXENE
      PARAMETER (MAXENE=900)
      REAL COMPTT1 , Deltal , sum , Theta , w , z , z7 , Z8 , zjx1
      INTEGER j , Jmax , Jx1
      REAL X(Jmax) , Dph(Jmax) , Rel(Jmax) , f(MAXENE) , f1(MAXENE)
      REAL f2(MAXENE)
      DOUBLE PRECISION z1 , z2 , z3 , z4 , z5 , z6
c     for the 0-theta F-P eq.
c     common /fp/ d1(900),d2(900)
      DO 100 j = 1 , Jx1
         w = X(j)
c     for the 0-theta F-P eq.
c       f(j)=dph(j)*w**2*d1(j)
c       f1(j)=dph(j)*w**3*d2(j)
         z = 1/(1+4.6*w+1.1*w*w)
         f(j) = Dph(j)*w**2*z
         f1(j) = -Dph(j)*w**3*z*z*(4.6+2.2*w)
 100  CONTINUE
      zjx1 = 1/(1+4.6*X(Jx1)+1.1*X(Jx1)*X(Jx1))
      DO 200 j = Jx1 , Jmax
         w = X(j)
         f2(j) = Dph(j)*w**2*Rel(j)
 200  CONTINUE
      sum = X(1)*f(1) + X(Jx1)*f(Jx1)
      DO 300 j = 2 , Jx1 - 1
         sum = sum + 2.*X(j)*f(j)
 300  CONTINUE
      z1 = sum*Deltal/2.
      sum = f(1) + f(Jx1)
      DO 400 j = 2 , Jx1 - 1
         sum = sum + 2.*f(j)
 400  CONTINUE
      z2 = sum*4.*Deltal/2.
      sum = f1(1) + f1(Jx1)
      DO 500 j = 2 , Jx1 - 1
         sum = sum + 2.*f1(j)
 500  CONTINUE
      z3 = sum*Deltal/2.
      sum = f2(Jx1) + f2(Jmax)
      DO 600 j = Jx1 + 1 , Jmax - 1
         sum = sum + 2.*f2(j)
 600  CONTINUE
      z4 = sum*Deltal/2.
      z5 = zjx1*X(Jx1)**2*Dph(Jx1)
      z6 = z5*X(Jx1)
c  the derivative term (neglected in LZ87)
      z7 = (Dph(Jx1)/X(Jx1)**2-Dph(Jx1-1)/X(Jx1-1)**2)/(X(Jx1)-X(Jx1-1))
     &     *X(Jx1)**5*zjx1
c  comptt1 is the Compton temperature
      COMPTT1 = SNGL((z1-z6+z4+Z8)/(z2+z3-z5+z7))
c  the average of the previous and the new temperature is used to avoid
c   oscillations
      COMPTT1 = (COMPTT1+Theta)/2
c  theta>0.15 causes instabilities; the method is not applicable then
c        write(6,*) 'Warning: comptt1=',comptt1
      IF ( COMPTT1.GT.0.15 .OR. COMPTT1.LT.1.E-10 ) COMPTT1 = 0.15
      RETURN
      END
**==annsp.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
      FUNCTION ANNSP(X,Temp,Tautom,Taup)
c     This program calculates the annihilation spectrum from a thermal
c     pure pair plasma of the total Thomson depth of tautom, which includes
c     the ionization electrons of the depth taup. (tautom=taup for pure
c     pair plasma.)
c     Using the formulae by Svensson 1983.
c     x and temp are in units of m c**2. The result gives photon number
c     spectrum in units of sigma(Thomson) c m c**2/ m c**2
      IMPLICIT NONE
      REAL a1 , a2 , ANNSP , cc , pi , Taup , Tautom , Temp , X , y , y4
      pi = 3.14159
      a1 = (X/Temp)**2/2/(1+2.0049/Temp+1.4774/Temp/Temp+pi/(2*Temp)**3)
      y = X*Temp
      IF ( y.LE.4. ) THEN
         a2 = pi/2*SQRT(pi/y)*EXP((2-X-1/X)/Temp)
         IF ( y.LE.0.25 ) THEN
            y4 = 4*y
            cc = 1 + 0.5288*y4 - 0.4483*y4**2 + 0.1643*y4**3
         ELSEIF ( y.LE.1. ) THEN
            cc = 1.125 + 0.6600*y - 0.7972*y**2 + 0.3060*y**3
         ELSE
            cc = 0.9049 + 1.1613/y - 1.2487/y**2 + 0.4770/y**3
         ENDIF
      ELSE
         a2 = pi/y*EXP((2-X)/Temp)*(LOG(4*0.56146*y)-1)
         y4 = 4/y
         cc = 1 + 0.5289*y4 - 0.8254*y4**2 + 0.9811*y4**3 - 0.3895*y4**4
      ENDIF
c  3/32/pi follow from assuming tautom=R sigma_T (n+ + n-)
      ANNSP = a1*a2*cc/Temp*3/32./pi*(Tautom**2-Taup**2)
      RETURN
      END
**==zerocn.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
c      FUNCTION GT0(CONST,P,RADIUS,B0)
c      implicit real (a-h,o-z)
c      implicit integer (i-n)
C  physical constants
c      COMMON /PHYSNT/ A0,SIGMAT,VLIGHT,bcrit,xlambd,xme,pi,xrs
C
C  THE NEW TURNOVER FREQUENCY
C
c      XTURN=(CONST*3.**((3+P)/2)/16.*137.*B0**(1+P/2))
c     >**(1/(2+P/2))
C  AND THE CORRESPONDING LORENTZ FACTOR
c      GT0=sqrt(3*XTURN/4/B0)
c      RETURN
c      END
 
      SUBROUTINE ZEROCN(Emagn,Aa,Cn,Cndis)
c    computes the initial parameters of the electron distribution
c  emagn -  the magnetic energy density
c  aa - the proportionality constant in the initial electron injection
c  cn - proportionality constant in the initial electron distribution
c      above gt
c  cndis - the same below gt
      IMPLICIT NONE
      REAL A0 , Aa , BCRit , BFIeld , Cn , Cndis , CONSTK , edens , 
     &     Emagn , EXP0 , fi , GAMma0 , GMAx , GMIn , GT , OMEga , PI , 
     &     RADius , SIGmat , TEMpbb
      REAL VLIght , XLAmbd , XLBb , XLEl , XME , XRS
C  physical constants
      COMMON /PHYSNT/ A0 , SIGmat , VLIght , BCRit , XLAmbd , XME , PI , 
     &                XRS
c  general parameters:
      COMMON /PARNT / GMIn , GMAx , EXP0 , BFIeld , TEMpbb , OMEga , 
     &                RADius , XLEl , XLBb , GAMma0 , GT
C
C  AA - ELECTRON INJECTION COEFFICIENT:
C    Q(GAMMA)=AA*GAMMA**(-EXP0)   CM**(-3)/SEC
c     take gamma0 as the last argument for electron injection, 0 for pair
      Aa = XLEl/XRS/CONSTK(EXP0,GMIn,GMAx,GAMma0)
      Emagn = BFIeld**2/8/PI/XME/VLIght**2*SIGmat*RADius
C     EDENS - DIMENSIONLESS ENERGY DENSITY OF MAGNETIC FIELD AND SOFT PHOTONS
c     for models with bremsstrahlung only, the initial energy density is
c     increased by alpha(fine)
      edens = Emagn + XLBb/XRS + 1/137.
C
      IF ( ABS(BFIeld).GT.1.E-15 ) THEN
C       FI - COEFFICIENT IN THE BINOMIAL SOLUTION FOR
c       INITIAL ELECTRON DISTRIBUTION (Zdziarski 1985, preprint)
         fi = 3*Emagn/edens**2
c       cn,cndis - proportionality constants in the electron distributions
         Cn = (SQRT(1+fi)-1)*edens/Emagn/2/CONSTK(EXP0,GMIn,GMAx,GT)/Aa
         Cndis = 1/edens
      ELSE
         Cn = 1/edens
      ENDIF
      RETURN
      END
**==eldis.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
      SUBROUTINE ELDIS(G,Aa,Cn,Cndis,Qint0,Del)
      IMPLICIT NONE
      REAL Aa , BFIeld , Cn , Cndis , EXP0 , GAMma0 , GMAx , GMIn , GT , 
     &     OMEga , RADius , TEMpbb , XLBb , XLEl
      INTEGER i , IMAx , IMIns , ITEr , JMC2
      REAL G(IMAx) , Qint0(IMAx) , Del(IMAx)
      COMMON /PARNT / GMIn , GMAx , EXP0 , BFIeld , TEMpbb , OMEga , 
     &                RADius , XLEl , XLBb , GAMma0 , GT
      COMMON /INTNT / IMIns , IMAx , JMC2 , ITEr
c  computes an INTEGRAL OVER THE the electron injection rate qint0
C  and the initial electron distribution del0, del1
C
      DO 100 i = 1 , IMAx
         IF ( EXP0.NE.0 ) THEN
            Qint0(i) = Aa*(MAX(G(i),GMIn)**(-EXP0+1)-GMAx**(-EXP0+1))
     &                 /(EXP0-1)
         ELSE
            Qint0(i) = Aa
         ENDIF
C       NORMALIZATION FOR COOLING:
         IF ( i.LE.IMIns ) THEN
            Del(i) = Qint0(i)*Cndis/(1.333*G(i)**2-1)
         ELSE
            Del(i) = Qint0(i)*Cn/(1.333*G(i)**2-1)
         ENDIF
 100  CONTINUE
      RETURN
      END
**==synchr.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
      SUBROUTINE SYNCHR(Dsyn,Del,G,X,Jsyn,B0,Syncon,Deltal)
      IMPLICIT NONE
      INTEGER MAXENE
      PARAMETER (MAXENE=900)
      REAL B0 , BFIeld , Deltal , el , EXP0 , gamma , GAMma0 , GMAx , 
     &     GMIn , GT , OMEga , RADius , Syncon , TEMpbb , XLBb , XLEl
      INTEGER i , j , Jsyn
      REAL Dsyn(MAXENE) , Del(MAXENE) , G(MAXENE) , X(Jsyn)
c  general parameters:
      COMMON /PARNT / GMIn , GMAx , EXP0 , BFIeld , TEMpbb , OMEga , 
     &                RADius , XLEl , XLBb , GAMma0 , GT
C  SYNCHROTRON EMISSION:
      DO 100 j = 1 , Jsyn
         gamma = SQRT(0.75*X(j)/B0)
         i = INT(LOG(gamma/GAMma0)/Deltal + 1)
c       the interpolated value for the electron distribution
         el = EXP(LOG(Del(i))+LOG(Del(i+1)/Del(i))
     &        /Deltal*LOG(gamma/G(i)))
         Dsyn(j) = Syncon*el/gamma
 100  CONTINUE
      RETURN
      END
**==clum.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
      FUNCTION CLUM(Dph,J1,J2,Const)
      IMPLICIT NONE
      INTEGER MAXENE
      PARAMETER (MAXENE=900)
      REAL Dph(MAXENE) , Const , sum , CLUM
      INTEGER j , J1 , J2
c     calculates a trapezoidal integral to obtain luminosity
c     dph is E F_E
      sum = (Dph(J1)+Dph(J2))/2
      DO 100 j = J1 + 1 , J2 - 1
         sum = sum + Dph(j)
 100  CONTINUE
      CLUM = sum*Const
      RETURN
      END
**==constk.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
      FUNCTION CONSTK(Exp0,Gmin,Gmax,Gt)
      IMPLICIT NONE
      REAL CONSTK , Exp0 , glower , Gmax , Gmin , Gt
C  SUBROUTINE COMPUTES AN INTEGRAL OVER THE INITIAL ELECTRON
C  INJECTION DISTRIBUTION
C  EXP0=0 - DELTA-FUNCTION ELECTRON INJECTION
C  OTHERWISE EXP0>1
C  EXP0=2 - INTEGRAL OF 1/X
C
      IF ( Gmax.LE.Gt ) CALL xwrite('gt > gmax',5)
      IF ( Exp0.EQ.0 ) THEN
         CONSTK = Gmax - Gt
      ELSE
         glower = MAX(Gt,Gmin)
         IF ( Exp0.EQ.2 ) THEN
            CONSTK = LOG(Gmax/glower) - Gt*(1/glower-1/Gmax)
         ELSE
            CONSTK = (Gmax**(-Exp0+2)-glower**(-Exp0+2))/(-Exp0+2)
     &               - Gt*(Gmax**(-Exp0+1)-glower**(-Exp0+1))/(-Exp0+1)
         ENDIF
      ENDIF
      RETURN
      END
**==table.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
      SUBROUTINE TABLE(Dph,X,Jmax,Del,G,N,Phcon)
c     tabulates the results
      IMPLICIT NONE
      INTEGER MAXENE
      PARAMETER (MAXENE=900)
      INTEGER i , IMAx , IMIns , ITEr , j , Jmax , JMC2 , N
      REAL Phcon
      REAL Dph(Jmax) , Del(IMAx) , X(Jmax) , G(IMAx)
      REAL ELDst(3,MAXENE) , PHDst(3,MAXENE)
      COMMON /INTNT / IMIns , IMAx , JMC2 , ITEr
      COMMON /NTHPLT/ ELDst , PHDst
      DO 100 j = 1 , Jmax
         PHDst(N,j) = Dph(j)*X(j)**2*Phcon
 100  CONTINUE
      DO 200 i = 1 , IMAx
         ELDst(N,i) = Del(i)*G(i)**2
 200  CONTINUE
      RETURN
      END
**==thbrem.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
      SUBROUTINE THBREM(X,Brth,Brcth,Theta,Tautom,Taup,Deltal,Jmax)
c     bremsstrahlung from thermal electrons
c     tautom is the total Thomson depth, taup is the ionization depth

      use fgsl
      IMPLICIT NONE
      REAL besexp , Brcth , Brth , Deltal , dp , Taup , Tautom , 
     &     THBR , Theta , X
      INTEGER j , Jmax , jmaxbr
      DIMENSION X(Jmax) , Brth(Jmax)
c     dp is the ratio of positrons to electrons+positrons, n_+/(n_+ + n_-)
      dp = 0.5 - Taup/2/Tautom
c     the maximum x for which the thermal bremsstrahlung is calculated
      jmaxbr = INT(LOG(50*Theta/X(1))/Deltal + 1)
c     besexp is K_2(1/theta)*exp(1/theta)
      IF ( Theta.GT.0.99 ) then
         besexp = SNGL(fgsl_sf_bessel_kcn(2, DBLE(1/Theta)))*
     &        EXP(1/Theta)
      endif
      DO 100 j = 1 , jmaxbr
         Brth(j) = THBR(X(j),Theta,dp,Tautom,besexp)
 100  CONTINUE
      DO 200 j = jmaxbr + 1 , Jmax
         Brth(j) = 0
 200  CONTINUE
c     integral thermal cooling rate
      Brcth = (Brth(1)*X(1)**2+Brth(jmaxbr)*X(jmaxbr)**2)/2
      DO 300 j = 2 , jmaxbr - 1
         Brcth = Brcth + Brth(j)*X(j)**2
 300  CONTINUE
      Brcth = Brcth*Deltal
      RETURN
      END
**==thbr.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
      FUNCTION THBR(E,Theta,Dp,Tautom,Besexp)
c     bremsstrahlung photon number spectrum from thermal electrons
c     e is the photon energy and theta the temperature, in units of m_e c^2
c     dp=n+/(n+ + n-)

      use fgsl
      IMPLICIT NONE
      REAL bes0 , bes1 , Besexp , Dp , E , e1x , ee0 , 
     &     ee1 , ee9 , ep0 , ep1 , ep2 , ep3 , ep4 , ep5 , ep9 , epm , 
     &     EXP1
      REAL f , pol1 , pol2 , Tautom , THBR , Theta , x , x2 , y
      x = E/Theta
C***  CHOICE OF THE APPRIOPRIATE FORMULA
      IF ( Theta.GE.1. ) THEN
C    THE ULTRARELATIVISTIC FORMULAS (THETA.GE.1.)
         y = EXP(-x)
         f = 1/Theta
c      EF=EXP(F)
         e1x = EXP1(x)
         pol1 = 6*x*x + 8*x + 16
         pol2 = pol1 - 16*x
C    ELECTRON-PROTON BREMSSTRAHLUNG
         ep0 = y*2./3./Besexp
         ep1 = LOG(2.)*(pol1*Theta*Theta+(8*x+16)*Theta+8)
         ep2 = e1x*(pol2*Theta*Theta+(-8*x+16)*Theta+8)
         ep3 = EXP1(f)*pol1*Theta*Theta
         ep4 = (-3*x*x+4*x+40)*Theta*Theta + (-4*x+16)*Theta
         ep5 = 3*x*x/(f+x) + EXP1(f+x)*(-3*x*x-4*x)
         ep9 = ep0*(ep1+ep2+ep3+ep4+ep5)
C    ELECTRON-ELECTRON BREMSSTRAHLUNG
         ee9 = y*2./3.*(pol2*e1x+2*(LOG(2*Theta)-.5772)
     &         *pol1+(3*x*x+12*x+56))
C    ELECTRON-POSITRON BREMSSTRAHLUNG
         epm = ee9
      ELSE
C     THE NONRELATIVISTIC FORMULAS (THETA.LE.1.)
         x2 = x/2
C     BESSKi CALCULATES VALUES OF THE K FUNCTION
         bes0 = SNGL(fgsl_sf_bessel_kc0(DBLE(x2)))
         bes1 = SNGL(fgsl_sf_bessel_kc1(DBLE(x2)))
C    ELECTRON-PROTON BREMSSTRAHLUNG
         ep0 = 4.2554/SQRT(Theta)*EXP(-x2)*bes0
         ep1 = ep0*(1+2*Theta+2*Theta**2+E*E*(.51875+.44018*Theta)
     &         +E*bes1/bes0*(.25+1.975*Theta+.073884*E*E))
     &         /(1+1.875*Theta+.82031*Theta*Theta-.30762*Theta**3)
         ep9 = ep1/(1+4.56*E/(168+E))
C    ELECTRON-ELECTRON BREMSSTRAHLUNG
         IF ( E.GE.1. ) THEN
            ee0 = ep9
            ee9 = ee0*(1+3*Theta)
         ELSE
            ee0 = ep0*Theta
            ee1 = ee0*2.2627*(.75+x2*bes1/bes0+.91438/bes0/EXP(x2))
            ee9 = ee1*(1+Theta)
         ENDIF
C    ELECTRON-POSITRON BREMSSTRAHLUNG
         epm = ep9*(1+3*Theta+(2.8284*(1-Theta)-1)/(1+E))
      ENDIF
C    THE FULL FORMULA; THE CONSTANT IS 3*ALPHA/8/PI
      THBR = 8.712E-4*(ep9*(1-2*Dp)+ee9*(Dp*Dp+(1-Dp)*(1-Dp))
     &       /2+epm*Dp*(1-Dp))*Tautom**2/E
      RETURN
      END
**==exp1.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
      FUNCTION EXP1(X)
C*** THE EXPONENTIAL INTEGRAL APPROXIMATION
c    TIMES EXP(X) TO AVOID UNDERFLOWS
      IMPLICIT NONE
      REAL EXP1 , f , s , X , xn , z
      INTEGER n
      IF ( X.LT.1. ) THEN
         f = 1
         s = 0.
         xn = 1
         DO 50 n = 1 , 5
            f = f*n
            xn = -xn*X
            s = s + xn/n/f
 50      CONTINUE
         EXP1 = -.5772 - LOG(X) - s
         EXP1 = EXP1*EXP(X)
      ELSE
         z = (X*X+2.334733*X+.250621)/(X*X+3.330657*X+1.681534)
         EXP1 = z/X
      ENDIF
      RETURN
      END
**==nthbrem.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
      SUBROUTINE NTHBREM(Brnth,Brcnth,G,Del,X,Deltal,Imax,Jmax,Taunth,
     &                   Tautom)
c     Gives bremsstrahlung emissivity from nonthermal electrons
c     interacting with cold plasma.
      IMPLICIT NONE
      REAL Brcnth , BRSIG , Deltal , Taunth , Tautom
      INTEGER i , Imax , iminbr , j , Jmax
      REAL Brnth(Jmax) , G(Imax) , X(Jmax) , Del(Imax)
      DO 100 j = 1 , Jmax
c     it used to be min, which caused problems ???
c     obviously it should be max
         iminbr = MAX(1,INT(LOG((1+X(j))/G(1))/Deltal+1))
c        iold=min(1,int(log((1+x(j))/g(1))/deltal+2))
c        write(6,*) j,jmax,iold,iminbr,imax
         Brnth(j) = (BRSIG(G(iminbr),X(j))*G(iminbr)*Del(iminbr)
     &              +BRSIG(G(Imax),X(j))*G(Imax)*Del(Imax))/2
         DO 50 i = iminbr + 1 , Imax - 1
            Brnth(j) = Brnth(j) + BRSIG(G(i),X(j))*Del(i)*G(i)
 50      CONTINUE
         Brnth(j) = Brnth(j)*Deltal*(Tautom+2*Taunth)
 100  CONTINUE
c     the nonthermal bremsstrahlung integral emissivity
      Brcnth = (Brnth(1)*X(1)**2+Brnth(Jmax)*X(Jmax)**2)/2
      DO 200 j = 2 , Jmax - 1
         Brcnth = Brcnth + Brnth(j)*X(j)**2
 200  CONTINUE
      Brcnth = Brcnth*Deltal
      RETURN
      END
**==brsig.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
      FUNCTION BRSIG(G,X)
c     cross section for relativistic bremsstrahlung
c     in units of sigma_T, 3.485e-3 is 3/2/pi*alpha
c     the electron energy is taken to be the kinetic one only
      IMPLICIT NONE
      REAL BRSIG , e , e2 , G , X
      e = G - 1
      e2 = G - 1 - X
      IF ( e2.GT.0. ) THEN
         BRSIG = 3.485E-3*e2/e*(e2/e+e/e2-.6667)*(LOG(2*e*e2/X)-0.5)/X
      ELSE
         BRSIG = 0
      ENDIF
      IF ( BRSIG.LT.0. ) BRSIG = 0
      RETURN
      END
**==annih.spg  processed by SPAG 4.50J  at 11:13 on 27 Feb 1996
 
      SUBROUTINE ANNIH(Theta,Tautom,Taup,Pprate,Gth,Annlum,Annht)
c     annihilation luminosity, and
c     residual annihilation heating
c     Calculate first the average thermal Lorentz factor
      use fgsl
      IMPLICIT NONE
      REAL Annht , Annlum , Gth , Pprate , Taup , 
     &     Tautom , Theta
      IF ( Theta.GT.0.014 ) THEN
         Gth = 3*Theta + SNGL(fgsl_sf_bessel_kc1(DBLE(1/Theta))/
     &             fgsl_sf_bessel_kcn(2,DBLE(1/Theta)))
      ELSE
         Gth = 1 + 1.5*Theta
      ENDIF
c     annihilation luminosity - a fit by Svensson (1982a)
      Annlum = 3/16./(1/(1+6*Theta)+Theta/(LOG(1.12*Theta+1)+0.25))
     &         *(Tautom**2-Taup**2)
c     residual heating
c  ********
c     replace pprate by the actual annihilation rate, which is important
c     for nonzero escape rate (e.g., use a fit by Svensson)
c     THIS IS NOT DONE YET
      Annht = Gth*Pprate - Annlum
      RETURN
      END
