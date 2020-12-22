      subroutine zigm(ear, ne, param, ifl, photar, photer)

* --------------------------------------------------------------------
* computes the stochastic Hydrogen attenuation due to intervening clouds
* User can select either the prescription of 
*    Avery Meiksin, 2006, MNRAS 365,805.
* or Madau 1995, ApJ, 441, 18.

* Code started by Martin Still Feb-26-2010 (NASA Ames)
* Updated by Frank Marshall March-29-2010 (NASA/GSFC)
* questions and comments: frank.marshall@nasa.gov                
*
* param(1) = redshift (z)
* param(2) = Madau (0) or Meiksin (1)    
* param(3) = Lyman limit on/off
* --------------------------------------------------------------------

      implicit none

c number of energy bins from xspec
      integer ne 
c spectrum number from xspec (apparently not used)         
      integer ifl
c energy loop
      integer ie          

c energy (keV) from xspec
      real ear(0:ne)
c output spectrum to xspec      
      real photar(ne)
c output error spectrum to xspec     
      real photer(ne)
c fit parameters from xspec     
      real param(3)
c redshift       
      real zfac           

* convert xspec parameters to subroutine parameters

      zfac = 1.0 + param(1)
      
* shift the energy array to the emitted frame

      do ie = 0, ne
         ear(ie) = ear(ie) * zfac
      enddo

* calculate the transmission

      call fmabsspecam(ear, ne, param, ifl, photar, photer)

* Now shift the energy array back again

      do ie = 0, ne
         ear(ie) = ear(ie) / zfac
      enddo

      return
      end

* --------------------------------------------------------------------
* absorb the source spectrum

      subroutine fmabsspecam(e, ne, param, ifl, photar, photer)

      implicit none
      integer ne,ifl
      real e(0:ne),param(*),photar(ne), photer(ne)

      integer i
      real xl,xh,sl,sh,fl,fh
      real zfac
      real fmmeiksin
      real fmmadau
      integer lylim
      integer metals
      real lls_fact

* conversion constant in keV*\AA

      real hc
      parameter (hc=12.3982362)

* photon index of the base power law function

      real index
      parameter (index=2.)

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

* convert xspec parameters to subroutine parameters

      zfac = 1.0 + param(1)
      lylim = int(param(3) + 0.5)
* code allows for a multiplicative factor to be applied
* to LLS in Meiksin, but it is not currently used.
      lls_fact = 1.0
* code allows for metals not to be used,
* but currently absorption due to metals is always used by Madau
      metals = 1
      
* calculate absorption at lambda
* average the values at the lower and upper energy bounds

      xl = hc/e(0)
      sl = xl**(-index)
      if (param(2) .gt. 0.5) then
        fl = fmmeiksin(xl,zfac,lylim,lls_fact)
      else
        fl = fmmadau(xl,zfac,lylim,metals)
      endif
      do i = 1, ne
         xh = hc/e(i)
         sh = xh**(-index)
       if (param(2) .gt. 0.5) then
           fh = fmmeiksin(xh,zfac,lylim,lls_fact)
       else
           fh = fmmadau(xh,zfac,lylim,metals)
         endif
         photar(i) = (sl*fl+sh*fh)/(sl+sh)
         xl = xh
         sl = sh
         fl = fh
      enddo

      return
      end

* --------------------------------------------------------------------
* Meiksin model

      function fmmeiksin(lambda,zfac,lylim,lls_fact)

* lambda is wavelength in the emitted frame
* zfac is 1+z
* lylim is flag whether to include photo-electric absorption
* lls_flag is factor to increase LLS optical depth

      implicit none

      real fmmeiksin    ! attenuation factor
      real lambda       ! wavelength (Angstroms)
      real zfac         ! redshift factor (z + 1)
      real lambda_lim   ! wavelength of Lyman limit (Angstroms)
      real lamn            ! wavelength for line n
      real tau_eff      ! optical depth
      real tau2            ! tau_alpha
      real tau_igm      ! tau due to optically thin clouds
      real tau_lls      ! tau due to LLS
      real xc           ! observed wavelength
      real xc25            ! xc**2.5
      real term1      ! part of LLS sum
      real term2      ! part of LLS sum
      real term3      ! part of LLS sum
      real term4      ! part of LLS sum
      real zn125      ! temporary value      
      integer lylim     ! switch for treating Ly limit (0 = off, 1 = on)
      real lls_fact     ! factor to increase optical depth due to LLS      
      integer n         ! counts Lyman series
      real an            ! =n
      integer m         ! counts LLS summation
      real am            ! =m
      integer nmax      ! maximum value for Lyman series
      integer mmax      ! maximum value for LLS summation
      real fact_n       ! keeps track of factorial(n)*(-1)**n
      real lambda_obs   ! observed wavelength
      real zn1          ! observed_lambda / lambda of transition n
      
      real fac(9)      ! relative Lyman Series coefficients

      parameter (nmax=31)
      parameter (mmax=10)
      parameter (lambda_lim=912.)
      
      data fac/0.,0.,0.348,0.179,0.109,0.0722,0.0508,0.0373,0.0283/
* fac contains the ratio tau_n to tau_alpha (n=2) for n=3 to 9

      tau_eff = 0.
      tau_igm=0.0
      tau_lls=0.0

* observed wavelength in Angstroms
      lambda_obs = lambda * zfac
* limit model to wavelengths more than 900 Angstroms (<13.8 eV)
      if (lambda_obs .gt. 900.0) then

* Attenuation due to Lyman Series lines
      lamn = lambda_lim / 0.75
      if (lambda .le. lamn) then
       zn1 = lambda_obs / lamn
       if (zn1 .le. 5.0) then
* Meiksin Eq. 2
          tau_eff = 0.00211*((zn1)**3.7)
       else
* Meiksin Eq. 3
          tau_eff = 0.00058*((zn1)**4.5)
       endif       
       do n = 3, nmax
        an=n
        lamn = lambda_lim / (1.0 - 1.0/an/an)
      if (lambda .gt. lamn) then
        exit
      else
        zn1 = lambda_obs / lamn
        zn125 = zn1*0.25
* I interpret Meiksin as using tau_alpha at zn1 in his Table 1
* calculate tau2 at zn1
        if (zn1 .le. 5.0) then
          tau2 = 0.00211*((zn1)**3.7)
          else
          tau2 = 0.00058*((zn1)**4.5)
          endif
        if (n .lt. 6) then
          if (zn1 .le. 4.0) then
            tau_eff = tau_eff + tau2*fac(n)*zn125**(1.0/3.0)
          else
            tau_eff = tau_eff + tau2*fac(n)*zn125**(1.0/6.0)
          endif
        else
          if (n .lt. 10) then
            tau_eff = tau_eff + tau2*fac(n)*zn125**(1.0/3.0)
          else
            tau_eff = tau_eff + tau2*fac(9)*zn125**(1.0/3.0) 
     &                * 720.0/an/(an*an -1.0)
          endif
        endif 
       endif 
       enddo
      endif
      
* Lyman continuum attenuation 
      if (lambda.lt.lambda_lim .and. lylim.eq.1) then      
      xc = lambda * zfac / lambda_lim
      xc25 = xc**2.5

* do IGM -- Meiksin Eq. 5
      tau_igm=0.805*xc**3 * (1.0/xc - 1.0/zfac)
      
* do LLS -- Meiksin Eq. 7
* term1 is gamma(0.5,1) - 1/e
        term1 = 0.2788 - 0.3679
        term2 = 0.0
      fact_n = 1.0
        do m = 0, mmax-1
       am = m
       if (m .gt. 0) then
         fact_n = -fact_n * am
       endif
         term2 = term2 + 1.0 / (fact_n * (2.0*am - 1.0))
        enddo
        term3 = zfac * xc**1.5 - xc25
        term4 = 0.0
      fact_n = 1.0
        do m = 1, mmax
       am = m
       fact_n = -fact_n * am
         term4 = term4 + (1.0 / 
     &            (fact_n * (2.0*am - 1.0) * (6.0*am - 5.0)) * 
     &            (zfac ** (2.5 - (3.0*am)) * xc**(3 * m) - xc25))
        enddo
        tau_lls = 0.25*((term1 - term2) * term3 - 2.0*term4) 
* multiply by specified factor to be able to account for variations across the sky
      tau_lls = tau_lls * lls_fact
      endif
      
      tau_eff=tau_eff + tau_igm + tau_lls
      endif

* check for reasonable values of tau_eff
      if (tau_eff .gt. 0.) then
        if (tau_eff .lt. 100.) then
          fmmeiksin = exp(-tau_eff)
      else
        fmmeiksin = 0.0
      endif
      else
        fmmeiksin = 1.0
      endif

      return
      end

* --------------------------------------------------------------------
* Madau model

      function fmmadau(lambda,zfac,lylim,metals)

* The coefficients for the Lyman series were provided by e-mail from P.Madau
* They differ slightly from Madau (1995)

      implicit none

      real fmmadau      ! attenuation factor
      real lambda       ! wavelength (Angstroms)
      real zfac         ! redshift factor (z + 1)
      real lambda_lim   ! wavelength of Lyman limit (Angstroms)
      real a_metal      ! optical depth coefficient of metals
      real tau_eff      ! optical depth
      real xc           ! observed wavelength
      real xc0            ! temporary variable
      
      real ly_coef(17)      ! Lyman Series coefficients
      real ly_wave(17)      ! Lyman Series wavelengths

      integer lylim     ! switch for treating Ly limit (0 = off, 1 = on)
      integer metals    ! switch for treating metal blanketing (0 = off, 1 = on)
      integer m         ! do loop index

      parameter (lambda_lim=911.75)
      parameter (a_metal=0.0017)
      
*2345678901234567890123456789012345678901234567890123456789012345678901    
      data ly_coef/0.0036,   0.0017,   0.0011846,0.0009410,0.0007960,
     &             0.0006967,0.0006236,0.0005665,0.0005200,0.0004817,
     &             0.0004487,0.0004200,0.0003947,0.000372 ,0.000352 ,
     &               0.0003334,0.00031644/
      data ly_wave/1215.67 ,1025.72 ,972.537,949.743,937.803,
     &              930.748, 926.226,923.15 ,920.963,919.352,
     &              918.129, 917.181,916.429,915.824,915.329,
     &              914.919, 914.576/

      tau_eff = 0.

* observed wavelength

      xc0= lambda * zfac
      xc = xc0 / lambda_lim
      
      if (xc0 .gt. 900.0) then

* Lyman alpha attenuation
      do m = 1, 17
        if (lambda.lt.ly_wave(m)) then
          tau_eff = tau_eff + ly_coef(m) * 
     &        (xc0 / ly_wave(m))**3.46
          if (metals.eq.1 .and. m.eq.1) then
          tau_eff = tau_eff + a_metal * 
     &        (xc0 / ly_wave(1))**1.68
          endif
      else
        exit
      endif 
      enddo
            
* Lyman continuum attenuation
* This uses the approximation given in footnote 3 to the
* integral in Eq. 16 of Madau (1995).
* It appears to be a poor approximation for observed wavelengths
* less than 912 Angstroms (xc < 1).
      if (lambda.lt.lambda_lim.and.lylim.eq.1) then
         tau_eff = tau_eff 
     &      + (0.25 * xc**3.0 * ((zfac**0.46) - (xc**0.46)))
     &      + (9.4 * xc**1.5 * ((zfac**0.18) - (xc**0.18)))
     &      - (0.7 * xc**3 * ((xc**(-1.32)) - (zfac**(-1.32))))
     &      - (0.023 * ((zfac**1.68) - (xc**1.68)))
      endif

* close if (xc0 .gt. 900     
      endif
 
* check for reasonable values of tau_eff
      if (tau_eff .gt. 0.) then
        if (tau_eff .lt. 100.) then
          fmmadau = exp(-tau_eff)
      else
        fmmadau = 0.0
      endif
      else
        fmmadau = 1.0
      endif

      return
      end
