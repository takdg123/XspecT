**==xspexrav.spg  processed by SPAG 4.50J  at 11:22 on 27 Feb 1996
      SUBROUTINE pexmon(Ear,Ne,Param,Ifl,Photar,Photer)
 
      IMPLICIT NONE
      INTEGER Ne , Ifl
      REAL Ear(0:Ne) , Param(7) , Photar(Ne), Photer(Ne)

c     Neutral Compton reflection with associated Fe and Ni lines. Built on top
c     of the pexrav model. Nandra et al. 2007, MNRAS 382, 194.

c     Driver for angle-dependent reflection from an exponentially-cutoff
c     power law and neutral medium. See Magdziarz & Zdziarski, 1995, MNRAS.
c
c     The output spectrum is the sum of the cutoff power law and the
c     reflection component.
c     The reflection component alone can be obtained
c     for scale (see below) = rel_refl < 0. Then the actual
c     reflection normalization is |scale|. Note that you need to
c     change then the limits of rel_refl. The range of rel_refl in that case
c     should exclude zero (as then the direct component appears).
c
c     If E_c=0, there is no cutoff in the power law.
c
c     Version with variable iron abundance and new opacities of Balucinska
c     & McCammon (1992, and 1994, private communication). As expected in AGNs,
c     H and He are assumed to be fully ionized.
c
c     Fe Kalpha, Fe Kbeta, Ni Kalpha, Fe Kalpha Compton shoulder included. The
c     line strength is based on Monte Carlo calculations by George & Fabian (1991)
c     parametrized by
c        EW = 9.66 EW_0 (Gamma^-2.8 - 0.56)    for 1.1 < Gamma < 2.5
c     with inclination dependence
c        EW = EW_0 (2.20 cos i - 1.749(cos i)^2 + 0.541(cos i)^3)  for i < 85 deg
c     and abundance dependence
c        log EW = log EW_0 (0.641 log A_Fe - 0.172 (log A_Fe)^2)
c     The Fe Kbeta and Ni Kalpha lines are given fluxes of 11.3% and 5%, 
c     respectively, of Fe Kalpha. The Fe Kalpha Compton shoulder
c     is approximated by a gaussian with E = 6.315 keV and sigma = 0.035 keV. The
c     inclination dependence is taken from Matt (2002) such that
c        EW_shoulder = EW_FeKalpha (0.1 + 0.1 cos i)
c     
c     number of model parameters: 7
c     1: Gamma, power law photon index, N_E prop. to E^{-Gamma}
c     2: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
c        one needs to change the lower limit for that)
c     3: scale, scaling factor for reflection; if <0, no direct component
c        (scale=1 for isotropic source above disk)
c     4: redshift, z
c     5: abundance of elements heavier than He relative to
c        the solar abundances
c     6: iron abundance relative to the solar iron abundance
c     7: inclination angle (deg)
c     algorithm:
c          a[E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
c     Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
c     of the cutoff power law only (without reflection)
c     and in the earth frame.
c
c

* KN
      INTEGER ie, ierr
      real reln_fe, reln_feb, reln_nia, reln_cs, zfac
      real ab_log, ab_qfit, inc_cfit

      REAL gparam(2)
      REAL gphotar(Ne), gphoter(Ne)

      LOGICAL firstcall

      SAVE firstcall
* KN 

      ierr = 0
 
      IF ( firstcall ) THEN
         CALL xwrite('Neutral Compton reflection with Fe emission',5)
         CALL xwrite('Fe line: George & Fabian 1991, 249, 352',5)
         CALL xwrite('Reflection: Magdziarz & Zdziarski 1995 MNRAS',5)
         CALL xwrite('Full description in: Nandra et al. 2007, MNRAS',5)
         firstcall = .FALSE.
      ENDIF

c *KN Change cos parameter
      param(7)=cos(param(7)*3.14159265/180.0)

c *KN Calculate relative norm from power law  (from FGF+ page)    
c Gamma (FGF+ parameterization)
      if (param(1).gt.1.1) then
         if (param(1).lt.2.5) then
            reln_fe = 4.75E-3*(9.66 * (param(1)**(-2.80)) - 0.56)
         else
            CALL xwrite ('*** pexmon: Gamma >2.5 - model invalid',5)
            reln_fe = 4.75E-3*0.182
         end if
      else
            CALL xwrite ('*** pexmon: Gamma <1.1 - model invalid',5)
            reln_fe = 4.75E-3*6.83
      end if   
c Inclination - cubic fit
      inc_cfit = 2.210 * param(7) - 1.749 *param(7) * param(7) 
     &     + 0.541 * param(7) * param(7) * param(7)
      reln_fe = reln_fe * inc_cfit 
c Abundance - quadratic fit to FGF Figure 17 
      if (param(6).lt.100.0) then
         if (param(6).lt.1.e-7) then
            reln_fe = 0. 
         else
            ab_log = log10(param(6))
            ab_qfit = ab_log * 0.641 - ab_log * ab_log * 0.172 
            reln_fe = reln_fe * 10.0 ** ab_qfit
         end if
      else
         CALL xwrite(' *** pexmon: Fe abun > 100 solar - model invaid',
     &               5)
      end if
c Special case for no line
      if (param(2).eq.999999.9) then
         reln_fe=0.0
      end if

c Rescale for Reflection fraction
      reln_fe = reln_fe * abs(param(3)) 

c K-beta and Nickel
      reln_feb = reln_fe * 17.0/ 150.0 
      reln_nia = reln_fe * 0.05

c Compton shoulder has approximate inclination dependence, from Matt (2002)
c Currently using gaussian (poor approximateion)
      reln_cs = reln_fe * (0.1 + param(7) * 0.1)

c Power-law reflection model

      call xspexrav(Ear,Ne,Param,Ifl,Photar,Photer)

c Redshift the energy array

      zfac = 1.0 + param(4)
      DO ie = 0, ne
         ear(ie) = ear(ie) * zfac
      ENDDO

c Add in the lines
      
      gparam(1)=6.4
      gparam(2)=0.005
      DO ie = 1, ne
         gphotar(ie) = 0.0
      ENDDO

      call gaussianline(ear, ne, gparam, ifl, gphotar, gphoter)
      DO ie = 1, ne
         photar(ie) = photar(ie) + gphotar(ie)* reln_fe
         gphotar(ie) = 0.0
      ENDDO


      gparam(1)=7.05
      call gaussianline(ear, ne, gparam, ifl, gphotar, gphoter)
      DO ie = 1, ne
         photar(ie) = photar(ie) + gphotar(ie)* reln_feb
         gphotar(ie) = 0.0
      ENDDO

      gparam(1)=7.47
      call gaussianline(ear, ne, gparam, ifl, gphotar, gphoter)
      DO ie = 1, ne
         photar(ie) = photar(ie) + gphotar(ie)* reln_nia
         gphotar(ie) = 0.0
      ENDDO

      gparam(1)=6.315
      gparam(2)=0.035
      call gaussianline(ear, ne, gparam, ifl, gphotar, gphoter)
      DO ie = 1, ne
         photar(ie) = photar(ie) + gphotar(ie)* reln_cs
         gphotar(ie) = 0.0
      ENDDO

c Unredshift the energy array

      DO ie = 0, ne
         ear(ie) = ear(ie) / zfac
      ENDDO

      RETURN
      END
