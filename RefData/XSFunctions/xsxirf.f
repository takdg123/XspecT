c updated version of the program. May 24 2001. 
**==xsxionref.spg  processed by SPAG 4.50J  at 18:52 on  3 Oct 2000
 
      SUBROUTINE XSXIRF(Ear1,Ne1,Param,Ifl,Photar1,Photer1)
      IMPLICIT NONE
 
      INTEGER Ne1 , Ifl
      REAL Ear1(0:Ne1) , Param(13) , Photar1(Ne1), Photer1(Ne1)

c Driver routine for Sergei Nayakshin's reflection code. Now doesn't have
c to do anything but run xion_ref.

      INTEGER i

c this model has no errors
      DO i = 1, Ne1
         Photer1(i) = 0.0
      ENDDO

      CALL XION_REF(Ear1,Ne1,Param,Ifl,Photar1)

      RETURN
      END
**==xion_ref.spg  processed by SPAG 4.50J  at 18:52 on  3 Oct 2000
 

c----------------------------------------------------------------
c
c     the SUBROUTINE for XSPEC 
C
C     
      subroutine xion_ref(ear1, ne1, param, ifl, photar1)

      implicit none

      INTEGER ifl, ne1
      real photar1(ne1), ear1(0:ne1)
      real param(13) 

c     this is number of bins in my tables (between ~0.1 and 200 keV)
      INTEGER nmaxp
      parameter (nmaxp = 582)   ! number of bins in phot. energy

      INTEGER ne, nr, ne2
      parameter (ne = nmaxp-2, ne2 = 80)
      parameter (nr = 30)! number of radial points
      REAL wp(nmaxp)  !photon energy grid in m_e c^2
      REAL spectr0(nmaxp), gamt(4)
      character(255) slovo



      REAL ear(0:ne), photar(ne), phout(ne), fxinc0(ne)
      real fxinc(ne), xnorm(ne), photnor(ne), eflux(nr), emis(nr),
     $ phnorm(ne), rd(nr), rds(nr), fx(nr), spectr(nr, nmaxp), 
     $     dgam(4), pwl(4), fxtemp(4, nmaxp), renorm,
     $     sptemp(4, nr, nmaxp), eart(ne2), phottmp(ne2), enel
      real mdot, param0(6), earmin, earmax, etgrow, eelow, eeup
      REAL h_x, eta_x1, cosv, gamma0, z, airon, ecut1, reftype, ecut0
      REAL elow0, tx0, tx1, ac, ap_cor, dera, ri, rilog10, dlogr
      REAL rilog1, extot, einc, drads, rads, ro, dvolume, phia, dphia
      REAL thea, dthea, ddist, defl, radius, omega, draa, rnorm
      REAL dgamto, a, b, gamma, eflux1
      REAL atmp(1),btmp(1)
      REAL phinc, powl, rinit
      REAL vphin, grin, cosfake

      INTEGER lenn, ilun, ierr, iilow, iiup, nshag, nkount
      INTEGER icall, n, kimax, kimax1, lo, ki, ji, igmin, igmax, ig
      INTEGER np, nread, ir, lnm, iener, len, block

c     chris' definitions
      logical up,low,qanyf
      real emax,emin,grad,y1,y2
      integer i,j,jup,jlow

      CHARACTER(255) contxt

      LOGICAL firstcall
      SAVE firstcall, wp

      INTEGER lenact
      CHARACTER(128) fgmodf
      EXTERNAL fgmodf, lenact

      DATA dera /57.295779/
c     this is normalization of the incident power law for the first point
      data pwl /4.8112E+05, 1.3158E+06, 2.9348E+06, 5.1738E+06/
      DATA firstcall/.TRUE./
      data gamt/1.6, 1.8, 2.0, 2.2/

c suppress a warning message from the compiler
      i = ifl

 
      IF ( firstcall ) THEN
         icall = 0

c         CALL xwrite('Photo-ionized X-ray reflection',5)
c         CALL xwrite('for incident power law with no exp. cutoff.',5)
c         CALL xwrite('Questions: Sergei Nayakshin',5)
         firstcall = .FALSE.
      else
         icall = 1
      ENDIF

c     ******************************
      h_x = param(1)            !in units of R_S

      eta_x1 = param(2)         !RATIO OF L_X/L_DISK

      mdot = param(3)           !ACCRETION RATE THROUGH THE DISK

      cosv = param(4)   !cosine of viewing angle

c inner and outer disk radii. Note that if param(13) is between 2 and 3,
c then R_in = param(1) by default

      if (param(13).ge.2.and.param(13).lt.2.999) then
         param0(3) = param(1) * 2.
      else
         param0(3) = param(5)*2. !IN UNTIS OF R_G = R_S/2!!!
      endif
      param0(4) = param(6)*2.

      gamma0 = param(7)

      z =  param(8) !redshift
c      print*, ' z = ', real(z)

      airon = param(9)!iron abundance

      ecut1 = param(10)         ! cutoff energy in keV

      reftype = param(11)
      
      ecut0 = 150.
      if (ecut1.gt.ecut0) ecut1 = ecut0
      elow0 = 0.1
      if (abs(gamma0-2.).le.3.e-3) then
         tx0 = log(ecut0/elow0)
         tx1 = log(ecut1/elow0)
         tx0 = tx0/(elow0**(1.-gamma0) - ecut0**(1.-gamma0))
         tx1 = tx1/(elow0**(1.-gamma0) - ecut1**(1.-gamma0))
      else
         tx0 = (ecut0**(-gamma0+2.) - elow0**(-gamma0+2.))/
     1        (-gamma0+2.)
         tx1 = (ecut1**(-gamma0+2.) - elow0**(-gamma0+2.))/
     1        (-gamma0+2.)
         tx0 = tx0/(elow0**(1.-gamma0) - ecut0**(1.-gamma0))
         tx1 = tx1/(elow0**(1.-gamma0) - ecut1**(1.-gamma0))
      endif

      ac = 2.
      eta_x1 = 1./((tx0/tx1) * (1. + 1./(ac * eta_x1)) - 1.)/ac
      ap_cor = eta_x1/param(2)  !this tells you by how much the
c     X-ray flux was requested to increase. This keeps the Compton
c     temperature at the approximately ``right'' value. However, it
c     artificially changes the A-parameter (because L_disk is really
c     what is kept constant by setting mdot). To avoid doing that, one
c     has to correct the A-parameter by multiplying it by ap_cor
      
c     ******************************

      param0(5) = acos(cosv) * dera
      param0(6) = param(12)



c     calculate radii. 
      ri = param0(3)
      rilog10 = alog10(ri)
      ro = param0(4)

         
      dlogr = (alog10(ro)-rilog10)/float(nr-1)
      rd(1) = ri
      DO n = 2, nr
         rd(n) = 10.**(rilog10+float(n-1)*dlogr)
      ENDDO

      if (param(12).ge.2.0) then
         rinit = param(5)
      else
         rinit = 3.
      endif

c     compute the disk illumination law for sphere
      if (param(13).ge.2.and.param(13).lt.2.999) then
         
         powl = 10. * (param(13)-2.)!powl is the index in the radial
c     power-law dependence of the X-ray emissivity within sphere
         rds(1) = rinit
         rilog1 = alog10(rds(1))
         dlogr = (alog10(param(1))-rilog1)/float(nr-1)
         DO n = 2, nr
            rds(n) = 10.**(rilog1+float(n-1)*dlogr)
         ENDDO
         extot = 0.
         einc = 0.
         kimax = 20
         kimax1 = kimax/2
         do n = 1, nr-1
            drads = rds(n+1) - rds(n)
            rads = (rds(n+1) + rds(n))/2.
c            do ki = 1, kimax
c               phia = float(ki) * 180./(float(kimax)*dera)
c     2 below comes from symmetry in phi-angle
c               dphia = 2.*  180./(float(kimax)*dera)
c               do ji = 1, kimax
c                  thea = float(ji) *90./(float(kimax)*dera)
c                  dthea = 90./(float(kimax1)*dera)
c                  dvolume = rads**2 * drads * sin(thea) 
c     $                 * dphia * dthea
c                  if (ji.eq.1.and.n.eq.nr-1) then
c                     emis(n) = 1.
c                  else
c                     emis(n) = 0.
c                  endif
c            emis(n) = (1.-sqrt(3./rads))/rads**4
c            emis(n) = 1./rads**4
            emis(n) = 1./rads**powl
c            print*, rads, emis(n)
                  dvolume = 2* 3.1415 * rads**2 * drads
                  extot = extot + emis(n) * dvolume
c               enddo
c            enddo
         enddo
         do 27 lo = 1, nr-1
            eflux(lo) = 0.
            do n = 1, nr-1
               drads = rds(n+1) - rds(n)
               rads = (rds(n+1) + rds(n))/2.
               do ki = 1, kimax
                  phia = float(ki) * 180./(float(kimax)*dera)
                  dphia =   180./(float(kimax)*dera)
                  do ji = 1, kimax1
                     thea = float(ji) *90./(float(kimax1)*dera)
                     dthea = 90./(float(kimax1)*dera)
                     dvolume = rads**2 * drads * sin(thea) 
     $                    * dphia * dthea
                     ddist = sqrt(rads**2 * sin(phia)**2 * 
     1                    sin(thea)**2 + (rads* cos(thea))**2 + 
     2                    (0.25*(rd(lo)+rd(lo+1))
     $                    - rads * sin(thea) * cos(phia))**2)
                     defl = rads * cos(thea)* dvolume/(ddist**3)
                     eflux(lo) = eflux(lo) + defl * emis(n)
                  enddo
               enddo
            enddo
                                !2 comes because of simetry in phi-angle
            eflux(lo) = 2. * eflux(lo)/(4. * 3.1415 * extot)
            einc = einc + 2.* 3.1415 * (rd(lo) + rd(lo+1)) *
     1           (rd(lo+1) - rd(lo)) * eflux(lo) /8.
 27      enddo
         
      endif

      do n = 1, nr-1
         radius = (rd(n) + rd(n+1))/4.
         if (param(13).ge.2.and.param(13).lt.2.999) then
            fx(n) = 2. * 3.1415 * eflux(n)/2
         endif
         if (abs(param(13)-1.).lt.1.e-6) then
c     note: fx(n) is the dimensionless incident flux per dr where
c     r is in R_G, not R_S (even though radius is in R_S), that's where
c     0.5 comes from.
            fx(n) = 0.5 * h_x/(h_x**2 + radius**2)**1.5 
         endif
c     for magnetic flares: similar to Shakura-Sunyaev
         if (abs(param(13)-3.).le.1.e-5) then
            fx(n) = (1. - sqrt(rinit/radius)) * (1./radius)**3
c            fx(n) = (1./radius)**3
         endif
C     THIS IS FOR A SPHERE ABOVE DISK -- WHEN R < R_SPHERE,
C     F_X = CONST
         if (abs(param(13)-4.).lt.1.e-5) then
            if (radius.ge.h_x) then
               fx(n) = (1. - sqrt(rinit/radius)) * (1./radius)**3
            else
               fx(n) = (1. - sqrt(rinit/h_x)) * (1./h_x)**3               
            endif
         endif
      enddo
      omega = 0. !angle covered by cold disk in Omega/2pi
      do n = 1, nr-1
         radius = (rd(n) + rd(n+1))/4.
         draa = rd(n+1) - rd(n)
         omega = omega +  fx(n) * draa * radius
      enddo
      rnorm = h_x/sqrt(h_x**2 + param(5)**2) - h_x/sqrt(h_x**2 + param(6
     $     )**2)
c      print*, rnorm, omega, '!'
      
      if (param(13).gt.2.999) then
         do n= 1, nr-1
            fx(n) = fx(n)/omega
         enddo
      endif
      if (param(13).gt.2.999) omega = 1.


c     CALCULATE SMEARING FUNCTION IF APPROXIMATE SOLUTION IS SOUGHT
      if (param(12).ge.3.0) then
c     this finds maximum and minimum fractional energy shift
c     and picks the grid
         grin = sqrt((rd(1)-2.)/(rd(1)-3.))
         vphin = 1./sqrt(rd(1)-2.)
         cosfake = 1. - param(4)
         earmin = 0.9 * sqrt(1.-2/rd(1))/(grin*(1.+ vphin))
         earmax = 1.1 * sqrt(1.-2/rd(1))/(grin*(1.- vphin))

         etgrow = (earmax/earmin)**(1./float(ne2-1))
         do i = 1, ne2
            eart(i) = earmin * etgrow**float(i-1)
            phottmp(i) = 0.
         enddo
         param0(1) = 1.

         call schw(eart, ne2, param0, nr, rd, fx, phottmp)

      endif 

      if (gamma0.lt.gamt(1)) then
         igmin =1
         igmax = 1
      endif

      if (gamma0.gt.gamt(4)) then
         igmin =4
         igmax = 4
      endif

      do ig = 1,3
         if (gamma0.gt.gamt(ig).and.gamma0.lt.gamt(ig+1)) then
            igmin = ig
            igmax = ig +1
         endif
      enddo

      do ig = 1,4
         if (gamma0.eq.gamt(ig)) then
            igmin = ig
            igmax = ig
         endif
      enddo

      dgam(igmin) = abs(gamma0 - gamt(igmin))
      dgam(igmax) = abs(gamt(igmax) - gamma0)
      dgamto = dgam(igmin) + dgam(igmax)

      if (icall.eq.0) then
         slovo = fgmodf()
         lenn = lenact(slovo)
         slovo = slovo(1:lenn)//'xiondata.fits'
         CALL getlun(ilun)
         call ftopen(ilun, slovo, 0, block, ierr)
         IF ( ierr .NE. 0 ) THEN
            contxt = 'XION_REF: Failed to open '//slovo(:lenact(slovo))
            CALL xwrite(contxt, 5)
            WRITE(contxt, '(a,i4)') '           iostat = ', ierr
            CALL xwrite(contxt, 5)
            RETURN
         ENDIF
         call ftmnhd(ilun, 2, 'XSPECTR16', 1, ierr)
         np = 0
         do nread = 1, 699
            call ftgcve(ilun, 1, nread, 1, 1, 0.0, atmp, qanyf, ierr)
            a = atmp(1)
            call ftgcve(ilun, 2, nread, 1, 1, 0.0, btmp, qanyf, ierr)
            b = btmp(1)
            if (a.gt.100.) then
               np = np + 1
               wp(np) = a/5.11e5
c               fxtemp(ig, np) = b/a
c               print*, a, real(fxtemp(ig, np))
            endif
         enddo
         do np = np, nmaxp-1
            wp(np+1) = wp(np) * (wp(np)/wp(np-1))
         enddo
         call ftclos(ilun, ierr)
         CALL frelun(ilun)
      endif
         

      iilow = 1
      iiup = ne
      do 1050 ig = igmin, igmax, 1
         

c     energy grid and illuminating power-law
c         if (abs(gamt(ig)-1.6).lt.1.e-5) then
c     slovo = '/homeA/lhea3/serg/lmodel/xspectr1.6'
c            slovo = 's1.6'
c         endif
c         if (abs(gamt(ig)-1.8).lt.1.e-5) then
c            slovo = 'xspectr1.8'
c         endif
c         if (abs(gamt(ig)-2.0).lt.1.e-5) then
c            slovo = 'xspectr2.0'
c            slovo = 'new-runs/xspectr2.0'
c         endif
c         if (abs(gamt(ig)-2.2).lt.1.e-5) then
c            slovo = 'xspectr2.2'
c         endif
         gamma = gamt(ig)

c         np = 0
c         do nread = 1, 699
c            read(11,*) a, b
c            if (a.gt.100.) then
c               np = np + 1
c               wp(np) = a/5.11e5
c               fxtemp(ig, np) = b/a
c            endif
c         enddo
c         do np = np, nmaxp-1
c            wp(np+1) = wp(np) * (wp(np)/wp(np-1))
c         enddo


         do i = 1, ne
            ear(i) = wp(i) * 5.11e5/(1. + z)
            spectr0(i) = 0.
            fxtemp(ig,i) = (pwl(ig)) * (wp(1)/wp(i))**
     1           (gamma)
         enddo

c     find the minimum and maximum energies to be treated

         iilow = 1
         iiup = ne
         if (param(12).ge.3.0) then
            eelow = ear1(0)*1.e3
            eeup = ear1(ne1)*1.e3

            do i = 1, ne
               if (ear(i).le.0.5*eelow) iilow = i
               if (ear(i).le.2.*eeup) iiup = i
            enddo
            if (iiup.eq.ne) iiup = ne -1
            if (iilow.le.1) iilow = 1
         endif

c     cccccccccccccc calculating the reflected spectra for all radii
         nshag = 5
         nkount = nshag-1
         if (param(12).ge.3.0) then
            do 2069 ir = 1, nr -1, 1
               nkount = nkount + 1
               if ((ir+nshag-1).le.nr-1) then
                  radius = (rd(ir) + rd(ir+nshag-1))/4.
                  eflux1 = (fx(ir) + fx(ir+nshag-1))/2.
               else
                  radius = (rd(ir) + rd(ir+1))/4.
                  eflux1 = (fx(ir) + fx(ir+1))/2.
               endif

               if (nkount.eq.nshag) then
c                  print*, ' nkount = ', nkount
                  nkount = 0

                  call x_illumination(h_x, radius, eta_x1,
     $                 ecut0, mdot, gamma, cosv, airon, icall, 
     $                 spectr0, ap_cor, eflux1, 
     $                 ig, rinit)
                  icall = icall +1
               endif
               do np = 2, ne
                  sptemp(ig, ir, np) = spectr0(np) *
     1                 (wp(np+1) - wp(np))*5.11e5
               enddo
                  
 2069       continue
               

         else
            do 2070 ir = 1, nr -1
               radius = (rd(ir) + rd(ir+1))/4.
               eflux1 = fx(ir)
               call x_illumination(h_x, radius, eta_x1, ecut0, mdot, 
     $              gamma, cosv, airon, icall, spectr0, ap_cor, eflux1, 
     $              ig, rinit)
               icall = icall +1
               
               do np = 2, ne
                  sptemp(ig, ir, np) = spectr0(np) *
     1                 (wp(np+1) - wp(np))*5.11e5
               enddo
 2070       enddo
         endif

         do i = 1, ne
            phout(i) = 0.
            phnorm(i) = 0.
            photar(i) = 0.
            photnor(i) = 0.
         enddo

c      close(11)

 1050 continue
      if (igmin.ne.igmax) then
         do np = 1, ne
            fxinc(np) = fxtemp(igmin, np)**(dgam(igmax)/dgamto)
     1           * fxtemp(igmax, np) ** (dgam(igmin)/dgamto)
         enddo
      else
         do np = 1, ne
            fxinc(np) = fxtemp(igmin, np)
         enddo
      endif

c     spectr(ir, np) is the reflected un-smeared spectrum as a function
c     of radial bin number (ir)

      do ir = 1, nr-1
         if (igmin.ne.igmax) then
            do np = 1, ne
               spectr(ir, np) = sptemp(igmin, ir, np)**
     1              (dgam(igmax)/dgamto) 
     2              * sptemp(igmax, ir, np)**(dgam(igmin)/dgamto)       
            enddo
         else
            do np = 1, ne
               spectr(ir, np) = sptemp(igmin, ir, np)
            enddo
         endif
      enddo

c      do lnm = 1, ne-1
      do lnm = iilow, iiup-1
         fxinc(lnm) = fxinc(lnm) *(ear(lnm+1) - ear(lnm))
         xnorm(lnm) = (ear(lnm+1) - ear(lnm))
      enddo

      
c     ccccccccccccccc smearing and integration over the whole disk
c     surface
      
c      do 5700 iener = 1, ne-1
      do 5700 iener = iilow, iiup
         len = iener
         param0(1) = ear(len)

         call xsd(len, ear, param0, rd, nr, photar, len,
     1            spectr, xnorm, photnor, fx)
         
         do i = 1, ne-1
            phout(i) = phout(i) + photar(i)
            phnorm(i) = phnorm(i) + photnor(i)
         enddo
         
 5700 continue


c      do i = 1, ne-1
      do i = iilow, iiup-1
         phnorm(i) = phnorm(i)/(ear(i+1)-ear(i))
         phout(i) = phout(i)/(phnorm(i)+1.e-10)
         spectr0(i)  = omega*phout(i)/(1.e-15 + fxinc(i))
      enddo

      if (param(12).ge.3.0) then      
c         do i = 1, ne
         do i = iilow, iiup
            spectr0(i) = 0
            fxinc0(i)=0.
         enddo
         
c         do 501 len = 1, ne-1
         do 501 len = iilow, iiup-1
c     len = 40
            enel = ear(len)
            call schw_tot(ear, ne, eart, phottmp, ne2, enel, photar, 
     1           renorm, len)
c            do i = 1, ne-1
            do i = iilow, iiup-1
c     spectr0(i) = spectr0(i) + spectr(1, len)
c     1              * photar(i)*(ear(i+1)-ear(i))/renorm
               spectr0(i) = spectr0(i) + phout(len)
     1              * photar(i)*(ear(i+1)-ear(i))/renorm
               fxinc0(i) = fxinc0(i) + fxinc(len)
     1              * photar(i)*(ear(i+1)-ear(i))/renorm
            enddo
 501     continue

c         do 502 i = 1, ne-1
         do 502 i = iilow, iiup-1
            spectr0(i) = omega*spectr0(i)/(fxinc0(i))
 502     enddo
      endif

c     sort out the photon spectrum on the initial grid - in eV
c      do i = 1, ne
      do i = iilow, iiup
         phinc = ear(i)**(-gamma0)
     &        * exp(-ear(i)/(ecut0*1.0e3))
         if (reftype.eq.1.) then
            phout(i) = phinc*(1. + cosv * spectr0(i))
         endif
         if (reftype.eq.2.) then
            phout(i) = phinc*(1. + cosv * spectr0(i)/omega)
         endif
         if (reftype.eq.3.) then
            phout(i) = phinc * spectr0(i)/omega
         endif
      enddo

c     check arrays for final energy grid are zero'd and put energy in eV
      ear1(0) = ear1(0)*1000.
      do i = 1, ne1
         photar1(i) = 0.
         ear1(i)=ear1(i)*1000.
      enddo

c     rebin onto output energy grid
      jlow=1
      jup=1
      do i=1,ne1,1
         low=.true.
         up=.true.
         emin=ear1(i-1)
         emax=ear1(i)
         do j=1,ne,1
            if (ear(j).gt.emin.and.low) then
               jlow=j
               low=.false.
            end if
            if (ear(j).gt.emax.and.up) then
               jup=j-1
               up=.false.
            end if
         end do          

         if (jup.eq.jlow-1) then
            grad=(phout(jup)-phout(jlow))/
     &           (ear(jup)-ear(jlow))
            y1=phout(jup)+grad*(emin-ear(jup))
            y2=phout(jup)+grad*(emax-ear(jup))
            photar1(i)=0.5*(y2+y1)*(emax-emin)
         end if
         if (jup.eq.jlow) then
            grad=(phout(jlow-1)-phout(jlow))
     &           /(ear(jlow-1)-ear(jlow))
            y1=phout(jlow-1)+grad*(emin-ear(jlow-1))
            photar1(i)=0.5*(phout(jlow)+y1)*(ear(jlow)-emin)
            grad=(phout(jlow)-phout(jlow+1))
     &           /(ear(jlow)-ear(jlow+1))
            y2=phout(jlow)+grad*(emax-ear(jlow))
            photar1(i)=photar1(i)+0.5*(emax-ear(jlow))*
     &           (phout(jlow)+y2)
         end if
         if (jup.gt.jlow) then
            grad=(phout(jlow-1)-phout(jlow))
     &           /(ear(jlow-1)-ear(jlow))
            y1=phout(jlow-1)+grad*(emin-ear(jlow-1))
            photar1(i)=0.5*(phout(jlow)+y1)*(ear(jlow)-emin)
            do j=jlow,jup-1,1
               photar1(i)=photar1(i)+0.5*(ear(j+1)-ear(j))
     &              *(phout(j+1)+phout(j))
            end do
            grad=(phout(jup+1)-phout(jup))
     &           /(ear(jup+1)-ear(jup))
            y2=phout(jup)+grad*(emax-ear(jup))
            photar1(i)=photar1(i)+0.5*(emax-ear(jup))
     &           *(y2+phout(jup))
         end if
       end do   
       

c     make the output to be the ratio of total/direct
       do i=1, ne1
          photar1(i) = photar1(i) * ear1(i)**(gamma0) 
     &         * exp(-ear1(i)/(ecut0*1000.))
     &         / (ear1(i) - ear1(i-1))
       enddo

c     then redo the output energy grid
       DO i = 0, ne1
          ear1(i) = ear1(i)/1000.
       ENDDO

       do i = 1, ne
          ear(i) = ear(i) * (1. + z)
       enddo
       
       return
       end



c --------------------------------------------------------      
      
      SUBROUTINE xsd(len, ear, param, rd, nr, photar,
     $     iclose, spectr, xnorm, photnor, fx)

      implicit none

      INTEGER ne, nr, iclose
      parameter (ne = 580)   ! number of bins in phot. energy

      REAL ear(0:ne), param(6), photar(ne), photar0(ne),
     1     xnorm(ne), photnor(ne), fx(nr), spectr(nr, ne)


c  model to calculate the line shape for a rotating accretion
c  disk. does not include GR effects. note that if param(2) is
c  set to 10 then do the special case of an accretion disk
c  emissivity.

c  parameters :
c	1        line energy
c       2        power law index for emissivity (10 for disk)
c       3        inner radius (GM/c**2)
c       4        outer radius (GM/c**2)
c       5        inclination  (degrees)


      REAL dera, 
     &     tgcsi2, cosdel, fra, ra, dra,
     &     rd(nr), enel
      REAL beta, betal, betah, zpo, zpol, zpoh
      REAL sincl, cincl, rafact, flux, total
      REAL enobl, enobh
      REAL radi, f_x, fra1, fluxnor
      INTEGER n, k, ie, len, kdeloc, ndeloc, nemin, nemax

      INTEGER ndel, kdel
      parameter (ndel = 60, kdel = 4)
      DATA dera /57.295779/

      rd(1) = param(3)


c convert input parameters - inclination set to radians.

      enel = param(1)
      sincl = sin(param(5)/dera)
      cincl = cos(param(5)/dera)

c set spectrum to 0

      if (param(6).gt.1.1) then
         photar(len) = 0.
         photnor(len) = 0.         
         DO n = 1, nr - 1
            ra = (rd(n)+rd(n+1))/2.
            dra = rd(n+1) - rd(n)
            f_x = fx(n) * ra * dra/2.
            photar(len) = photar(len) + f_x * spectr(n, len)
            photnor(len) = photnor(len) + f_x * xnorm(len)
c     endif
c     ENDDO
         enddo
         go to 2000
      endif

      DO n = 1, ne
         photar(n) = 0.
         photnor(n) = 0.         
         photar0(n) = 0.
      ENDDO

c trap case where inner radius is greater than outer

c big loop for radii. note that the radius cannot reach 3 else
c     the metric goes singular

      DO n = 1, nr - 1

         ra = (rd(n)+rd(n+1))/2.

         dra = rd(n+1) - rd(n)

         rafact = sqrt(1.-3./ra)

         radi = ra/2.
         f_x = fx(n)
         fra = f_x * spectr(n, len)
         fra1 = f_x * xnorm(len)


c loop over azimuthal angles

         kdeloc = INT(kdel * sqrt(ra/6.))
         if (kdeloc.gt.18) kdeloc = 18
c         kdeloc = 2
         DO k = 1, 179, kdeloc


            beta = float(k)/dera
            betal = float(k-1)/dera
            betah = float(k+1)/dera

c calculate mean redshift (1+z = zpo) for the bin

            tgcsi2 = (sincl*sin(beta))
     &               **2/(1.-(sincl*sin(beta))**2)
            cosdel = sincl*cos(beta)
     &               /sqrt((sincl*cos(beta))**2+cincl**2)
            zpo = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact

            if (param(6).eq.2.) zpo = 1.

c     and the low and high redshifts for the bin. note the traps for
c     the case of an inclination of 90 degrees.
            
            
            IF ( param(5) .GT. 89.9 .AND. k .EQ. 91 ) THEN

               zpol = 1./rafact
               if (param(6).eq.2.) zpol = 1.               

            ELSE

               tgcsi2 = (sincl*sin(betal))
     &                  **2/(1.-(sincl*sin(betal))**2)
               cosdel = sincl*cos(betal)
     &                  /sqrt((sincl*cos(betal))**2+cincl**2)
               zpol = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact
               if (param(6).eq.2.) zpol = 1.                              

            ENDIF

            IF ( param(5) .GT. 89.9 .AND. k .EQ. 89 ) THEN

               zpoh = 1./rafact
               if (param(6).eq.2.) zpoh = 1.

            ELSE

               tgcsi2 = (sincl*sin(betah))
     &                  **2/(1.-(sincl*sin(betah))**2)
               cosdel = sincl*cos(betah)
     &                  /sqrt((sincl*cos(betah))**2+cincl**2)
               zpoh = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact
               if (param(6).eq.2.) zpoh = 1.               

            ENDIF

c  enobl and enobh are the lower and upper observed energy from
c  this azimuthal and radial bin

            enobl = min(enel/zpol, enel/zpoh)
            enobh = max(enel/zpol, enel/zpoh)
            total = enobh - enobl

c  calculate total emission from this bin

            flux = 0.5*float(kdeloc) * fra*ra*dra/(179. * zpo**3)
            fluxnor = 0.5*float(kdeloc) * fra1*ra*dra/(179. * zpo**3)


c  find fractions of emission from this bin to place in each energy
c  range.

            IF (enobh .LT. ear(0) .OR. enobl .GT. ear(ne) ) GOTO 10

            ndeloc = INT(ndel * (0.1 + sqrt(6./ra)))
c            ndeloc = ndel
c            nemin = 1
c            nemax = ne
            nemin = max(1, iclose - ndeloc)
            nemax = min( ne, iclose + ndeloc)
            DO ie = nemin, nemax
               IF ( ear(ie) .GE. enobh ) THEN
                  IF ( ear(ie-1) .LE. enobl ) THEN
                     photar(ie) = photar(ie) + flux
                     photnor(ie) = photnor(ie) + fluxnor                     
c                     photar0(ie) = photar0(ie) + flux0                     
                     GOTO 10
                  ELSEIF ( ear(ie-1) .GE. enobh ) THEN
                     GOTO 10
                  ELSE
                     IF ( total .GT. 0. ) THEN
                        photar(ie) = photar(ie) 
     &                       + flux*(enobh-ear(ie-1))/total
                        photnor(ie) = photnor(ie) 
     &                       + fluxnor*(enobh-ear(ie-1))/total
c                        photar0(ie) = photar0(ie) 
c     &                       + flux0*(enobh-ear(ie-1))/total
                     ENDIF
                  ENDIF
               ELSEIF ( ear(ie) .GE. enobl ) THEN
                  IF ( ear(ie-1) .GE. enobl ) THEN
                     IF ( total .GT. 0. ) THEN
                        photar(ie) = photar(ie) 
     &                       + flux*(ear(ie)-ear(ie-1))/total
                        photnor(ie) = photnor(ie) 
     &                       + fluxnor*(ear(ie)-ear(ie-1))/total

c                        photar0(ie) = photar0(ie) 
c     &                     + flux0*(ear(ie)-ear(ie-1))/total                        
                     ENDIF
                  ELSE
                     IF ( total .GT. 0. ) THEN
                        photar(ie) = photar(ie) 
     &                       + flux*(ear(ie)-enobl)/total
                        photnor(ie) = photnor(ie) + fluxnor*(ear(ie)
     $                       -enobl)/total 
c                        photar0(ie) = photar0(ie) 
c     &                       + flux0*(ear(ie)-enobl)/total
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

 10         CONTINUE

         ENDDO
      ENDDO

 2000 RETURN
      END


c     -----------------------------------------
c     program that solves the illumination problem approximately.

      subroutine x_illumination(h_x, radius, eta_x1, ecut0, mdot, 
     $     gamma, cosv, airon, icall, spectr, ap_cor, eflux1, ig, 
     $     rinit)

      implicit none

      INTEGER nmaxp, mumax_v, nbegin, nend
      parameter (nmaxp = 582)   ! number of bins in phot. energy
      parameter (mumax_v = 5)  ! same but for output spectra

      INTEGER ntest, ncmax
      parameter (ntest = 2882, ncmax = 16)
c      parameter (ntest = 1441, ncmax = 16)      

      REAL int_out(ntest, nmaxp),
     2     mu_v(mumax_v), dmu_v(mumax_v), spectr(nmaxp),
     3     gam1(ntest), ecut1(ntest), fxd1(ntest), 
     4     tauex(ntest), d(ncmax), xdrat(9), 
     5     weight(ncmax), afe(3)
      REAL grav(ntest), abund(ntest)

      integer iang(ntest), nc(ncmax), ife, ig

      REAL ergsev, sigma_t, taumin1, eta_x, eta_x1, ampl1, gioniz
      REAL gamma, radius, fdfx, mdot, apar, h_x, eflux1, ap_cor
      REAL fxfda, ecut0, cosv, airon, fxdtem, costem
      REAL afetemp, dgrav, dtot, sum, devtot, dfxd, diang, dabun
      REAL v1, v2, v3, v4, rinit
      INTEGER lenn, ilun, ierr, jh, m, icall, i, np, ntmax
      INTEGER jhmin, jhmax, jmmin, jmmax, ifemin, ifemax, j
      INTEGER ncexc, nclose, io, jcount, block
      CHARACTER(255) filenm
      LOGICAL qanyf

      INTEGER lenact
      CHARACTER(128) fgmodf
      EXTERNAL lenact, fgmodf

      SAVE abund, gam1, ecut1, fxd1, tauex, iang, grav, int_out, ntmax

      ergsev = 1.602197e-12
      sigma_t = 6.65e-25

c constants: m_e c^3/F_x = 2.45 * 10^{-13} /F_17 cm^3, if F_17 
c is flux/10^17 erg/sec/cm^2

      taumin1 = 0.001
      eta_x = eta_x1
      ampl1 = 1.
      gioniz = gamma

      call approximate(radius, fdfx, mdot, eta_x, apar, h_x, eflux1, 
     1     rinit)

      apar = apar * ap_cor

      fxfda = 1./fdfx

      ECUT0 = 20000.
      do jh = 1, 9
         xdrat (jh) = 2.5**float(jh-4)
      enddo

      do m = 1, mumax_v
c         mu_v(m) = (m-0.5)/(float(mumax_v))
         dmu_v(m) = 1./float(mumax_v)
         mu_v(m) = 0.15 + 0.2 * float(m-1)
      enddo

      afe(1) = 0.2
      afe(2) = 1.
      afe(3) = 4.
         

      if (icall.gt.0) go to 130

      i = 1
      filenm = fgmodf()
      lenn = lenact(filenm)
      filenm = filenm(1:lenn)//'xiondata.fits'      
      CALL getlun(ilun)
      CALL ftopen(ilun, filenm, 0, block, ierr)
      CALL ftmnhd(ilun, 2, 'NEWTABLE12', 1, ierr)
      CALL ftgnrw(ilun, ntmax, ierr)
      ntmax = min(ntmax, ntest)
      CALL ftgcve(ilun,1,1,1,ntmax,0.0,abund,qanyf,ierr)
      CALL ftgcve(ilun,2,1,1,ntmax,0.0,gam1,qanyf,ierr)
      CALL ftgcve(ilun,3,1,1,ntmax,0.0,ecut1,qanyf,ierr)
      CALL ftgcve(ilun,4,1,1,ntmax,0.0,fxd1,qanyf,ierr)
      CALL ftgcve(ilun,5,1,1,ntmax,0.0,tauex,qanyf,ierr)
      CALL ftgcvj(ilun,6,1,1,ntmax,0.0,iang,qanyf,ierr)
      CALL ftgcve(ilun,7,1,1,ntmax,0.0,grav,qanyf,ierr)

      DO i = 1, ntmax
         iang(i) = iang(i)/2
         CALL ftgcve(ilun,10,i,1,nmaxp,0.0,spectr,qanyf,ierr)
         DO j = 1, nmaxp-1
            int_out(i,j+1) = spectr(j)
         ENDDO
      ENDDO

c      print*, 'number of models in the table: ', ntmax-1

      CALL ftclos(ilun,ierr)
      CALL frelun(ilun)

 130  continue

c     find the closest gridpoint in FX/Fd ratio
c     

      jhmin = 1
      jhmax = 1
      if (fxfda.lt.xdrat(1)) go to 135

      if (fxfda.gt.xdrat(9)) then
         jhmin = 9
         jhmax = 9
         go to 135
      endif


      do jh = 1, 8
         if (xdrat(jh).lt.fxfda.and.xdrat(jh+1).gt.fxfda) then
            jhmin = jh
            jhmax = jh+1
         endif
      enddo

      do jh = 1, 9
c         if (xdrat(jh).eq.fxfda) then
         if (abs(xdrat(jh)-fxfda).le.0.0001*(fxfda+xdrat(jh))) then
            jhmin = jh
            jhmax = jh
         endif
      enddo

 135  continue

c     find the closest gridpoint in angle
c     

      jmmin = 1
      jmmax = 1
      if (cosv.le.mu_v(1)) go to 140

      if (cosv.ge.mu_v(mumax_v)) then
         jmmin = mumax_v
         jmmax = mumax_v
         go to 140
      endif

      do m = 1, mumax_v-1
         if (mu_v(m).le.cosv.and.mu_v(m+1).ge.cosv) then
            jmmin = m
            jmmax = m +1
         endif
      enddo

      do m = 1, mumax_v
         if (abs(mu_v(m)-cosv).le.0.0001) then
            jmmin = m
            jmmax = m 
         endif
      enddo

 140  continue

c find the right abundances. 


      ifemin = 1
      ifemax = 1
      do ife = 1, 3
         if (abs(airon-afe(ife)).le.0.01 * afe(ife)) then
            ifemin = ife
            ifemax = ife
            go to 145
         else
            if (airon.le.afe(1)) then
               ifemin = 1
               ifemax = 1
               go to 145
            endif
            if (airon.ge.afe(3)) then
               ifemin = 3
               ifemax = 3
               go to 145
            endif
            if (airon.gt.afe(1).and.airon.lt.afe(2)) then
               ifemin = 1
               ifemax = 2
               go to 145
            endif
            if (airon.gt.afe(2).and.airon.lt.afe(3)) then
               ifemin = 2
               ifemax = 3
               go to 145
            endif
         endif
      enddo

 145  continue

      nbegin = 1
      nend = ntmax-1
c     search only spectra with the right abundances in the table
      if (ifemin.eq.2.and.ifemax.eq.2) then
         nbegin = 1
         nend = 1440

         if (ig.eq.1) then
            nend = nbegin + nend/4
         endif

         if (ig.eq.2) then
            nbegin = nbegin + nend/4
            nend = nbegin + nend/4
         endif

         if (ig.eq.3) then
            nbegin = nbegin + nend/2
            nend = nbegin + nend/4
         endif

         if (ig.eq.4) then
            nbegin = nend - nend/4
         endif

      endif
      if (ifemin.eq.3.and.ifemax.eq.3) then
         nbegin = 1441
         nend = 2880
      endif



      do j = 1, ncmax
         d(j) = 1.e6
         nc(j) = 1000000
      enddo

      ncexc = 1000000
      fxdtem = xdrat(jhmin)
      costem = mu_v(jmmin)
      afetemp = afe(ifemin)
      do j = 1, ncmax
         if (j.le.2) then
            fxdtem = xdrat(jhmin)
            costem = mu_v(jmmin)
            afetemp = afe(ifemin)
            if (j.eq.2) ncexc = nc(j-1)
         endif
         if (j.gt.2.and.j.le.4) then
            fxdtem = xdrat(jhmin)
            costem = mu_v(jmmin)
            afetemp = afe(ifemax)
            if (j.eq.4) ncexc = nc(j-1)
         endif
         if (j.eq.5.or.j.eq.6) then
            fxdtem = xdrat(jhmin)
            costem = mu_v(jmmax)
            afetemp = afe(ifemax)
            if (j.eq.6) ncexc = nc(j-1)
         endif
         if (j.eq.7.or.j.eq.8) then
            fxdtem = xdrat(jhmax)
            costem = mu_v(jmmin)
            afetemp = afe(ifemax)
            if (j.eq.8) ncexc = nc(j-1)
         endif
         if (j.eq.9.or.j.eq.10) then
            fxdtem = xdrat(jhmax)
            costem = mu_v(jmmax)
            afetemp = afe(ifemax)
            if (j.eq.10) ncexc = nc(j-1)
         endif
         if (j.eq.11.or.j.eq.12) then
            fxdtem = xdrat(jhmax)
            costem = mu_v(jmmin)
            afetemp = afe(ifemin)
            if (j.eq.12) ncexc = nc(j-1)
         endif
         if (j.eq.13.or.j.eq.14) then
            fxdtem = xdrat(jhmax)
            costem = mu_v(jmmax)
            afetemp = afe(ifemin)
            if (j.eq.14) ncexc = nc(j-1)
         endif
         if (j.eq.15.or.j.eq.16) then
            fxdtem = xdrat(jhmin)
            costem = mu_v(jmmax)
            afetemp = afe(ifemin)
            if (j.eq.16) ncexc = nc(j-1)
         endif

c         do i = 1, ntmax-1
         do i = nbegin, nend, 1

c     the poor acuracy of the comparison below is to allow
c     for roundoff errors on systems such as linux.... In reality it
c     does not matter because the grid points are separated by more than
c     the error below

            if (abs(gamma-gam1(i)).gt.1.e-4.or.
     $           abs(afetemp-abund(i)).gt.0.01*(afetemp+abund(i))) then
               go to 200
            endif
            if (abs(fxd1(i)-fxdtem).gt.(0.02*fxdtem).or.
     $           abs(mu_v(iang(i))-costem).gt.1.e-2) then
               go to 200
            endif
            
            if (i.eq.ncexc) go to 200

            if (grav(i)/apar.ge.1.) then
               dgrav = grav(i)/apar - 1.
            else
               dgrav = apar/grav(i) - 1.
            endif

            dtot = abs(dgrav)

            if (dtot.le.d(j)) then
               d(j) = dtot
               nc(j) = i
            endif
 200        continue
         enddo

      enddo

      nclose = nc(1)
      io = nclose

      do np = 1, nmaxp-1
         spectr(np) = 0.
      enddo

      sum = 1.


      devtot = 0.
      jcount = 0
      do j = 1, ncmax
        i = nc(j)
        jcount = jcount + 1


        if (fxfda.le.fxd1(i)) then
           dfxd = fxd1(i)/fxfda - 1.
        else
           dfxd = fxfda/fxd1(i) - 1.
        endif

        diang =(cosv - mu_v(iang(i)))/0.1

        if (grav(i)/apar.ge.1.) then
           dgrav = grav(i)/apar - 1.
        else
           dgrav = apar/grav(i) - 1.
        endif


        if (airon.ge.abund(i)) then
           dabun = airon/abund(i) - 1.
        else
           dabun = abund(i)/airon- 1.
        endif

        v1 = abs(dfxd)
        v2 = abs(diang)
        v3 = abs(dgrav)
        v4 = abs(dabun)

c     if difference between the variable and the grid point
c     is zero, I can just use weight 1 because it means that
c     the value at the grid point will be taken

         if (v1.le.1.e-3) v1 = 1.e-3
         if (v2.le.1.e-3) v2 = 1.e-3
         if (v3.le.1.e-3) v3 = 1.e-3
         if (v4.le.1.e-3) v4 = 1.e-3

         weight(j) = 1./(v1 * v2 * v3 * v4)

         if (jcount.eq.2) then
            jcount = 0
         endif

         devtot = devtot + weight(j)

      enddo

      do j = 1, ncmax
         do np = 2, nmaxp-1
            spectr(np) =  spectr(np) + int_out(nc(j), np)
     1          * weight(j)/devtot
c            spectr(np) = int_out(nclose, np)
         enddo
      enddo

      return
      end


c     ------------------------------------------------------------------
c     program tau_hot
      subroutine approximate(radi, fdfx, mdot, eta_x, apar, h_x, eflux1,
     1     rinit)

      implicit none

      real m8, mdot, mdot0, fcor, alphaf, rinit

c     parameters: m8 = mass in 10^8 solar masses; fcor = fraction of
c     disk energy chanelled into the corona; alphaf = viscosity parameter
c     h_x = height of the source above the disk in R_s
      parameter (m8 = 3.e0, fcor = 0.001, alphaf = 1.e-2)

      real loclum, speedl, flconst, sigmat, hratio, pdisk, tauss
      REAL thetass, prpg, radi, h_1, h_x, f_x, eta_x, eflux1
      REAL f_d, fdfx, apar

      REAL xsxfd
      external xsxfd
      
      speedl = 3.e10

      flconst = 2.26e18/m8      !m_pc^3/sigma_t R_s
      sigmat = 6.65e-25

      mdot0 = mdot

      call shakura_sunyaev(hratio, loclum, mdot0, pdisk,
     1     tauss, thetass, prpg, radi, rinit)
            
      h_1 = h_x
      f_d = mdot0 * xsxfd(radi)
      f_x = eta_x * mdot0 * eflux1      

      apar = (1./(2* radi**2)) * hratio/f_x
      apar = apar * sqrt(radi * 8/(1836.*511.)) /hratio

      fdfx = f_d/f_x

         
      return
      end
      
c-----------------------------------------------------------------
c solves for the disk structure in the standard Shakura-Sunyaev model
c for a given radius and accretion rate WITH ARBITRARY VISCOSITY PRESCRIPTION
c used to give initial values of the variables

      subroutine shakura_sunyaev(hratio, loclum, mdot0, pdisk,
     1     tauss, thetass, prpg, radl, rinit)

      implicit none

c     on the input:
c     rad = radius in R_g
c     alpha = viscosity parameter
c     m8 = mass in solar masses
c     mdot0 = accretion rate 
c     fcor = fraction of power channeled to the corona

c on the output: 
c prpg = radiation to gas pressure ratio
c hratio is H/R
c thetass is the gas temperature in the disk midplane
c tauss is the Total Thomson depth of the disk.

      real m8, mdot0, fcor, alphaf


c parameters: m8 = mass in 10^8 solar masses; fcor = fraction of
c disk energy chanelled into the corona; alphaf = viscosity parameter
      parameter (m8 = 3.e0, fcor = 0.001, alphaf = 1.e-2)


      REAL lumin, loclum, eta, aconst, sigma_t, ergsev, c, x_1, xrmax
      REAL xrmin, funj, radl, const, prpg, theta_1, tau_1, h_1, pr
      REAL pg, hratio, tauss, thetass, pdisk , rinit
      integer j, itermax
      parameter (itermax = 150, eta = 0.06)
      real xr(itermax + 1), m1
      character(72) contxt

      m1 = 1e7 * m8
      aconst = 4.13e9 * m1
c      radl = radi
      sigma_t = 6.65e-25
      ergsev = 1.602197e-12


      c = 0.8
      x_1 = 1.0
      
      xrmax = 0.1
      xrmin = 1.0e-10

      funj = 1. - sqrt(rinit/radl)

      const = 9./sqrt(8.)
      xr(1) = 0.0001
      prpg = 0.1

      do 25 j = 2, itermax
         loclum = (3./8.) * mdot0 * funj/(eta * radl**3)
         loclum = loclum * (1 - fcor)
c pick a temperature            
         theta_1 = xr(j-1)

c find corresponding optical depth
         tau_1 = 2. * aconst * (theta_1**4)/loclum
c find h
         h_1 = 2. * aconst * (theta_1**4) * (radl**3)/tau_1

         h_1 = SNGL(h_1 + dsqrt(8.d0 * theta_1 * (radl**3)/1836.d0 +
     1        h_1**2))

         pr = 4. * aconst * theta_1 ** 4
         pg = 8. * tau_1 * theta_1/(1836. * h_1)
         prpg = pr/pg
         hratio = h_1/radl
            
         lumin = 0.5 * const * radl ** (-4.5) * 
     1        h_1 **2 * tau_1 * alphaf
     2        * (1. - fcor)
c this would mean the lightman-eardley prescription.
c /(1. + prpg)

         if (lumin.lt.loclum) then
            xr(j) = xr(j-1)**c * xrmax ** (1-c)
            xrmin = xr(j-1)
         else
            xr(j) = xr(j-1) **c * xrmin ** (1-c)
            xrmax = xr(j-1)
         endif

 25   enddo
      if (abs(lumin-loclum).gt.0.01 * loclum) then
         WRITE(contxt,'(a)') 'Shakura-Sunyaev model solution failed'
         CALL xwrite(contxt, 5)
         WRITE(contxt,*) real((lumin-loclum)/loclum), prpg
         CALL xwrite(contxt, 5)
      endif
      tauss = tau_1
      thetass = xr(itermax)

       pdisk = 0.25 * (1836 * 5.11e5 * 1.6e-12/(6.65e-25 *
     1     3e13 * m8)) * tauss * hratio/(radl**2)

      return
      end



c     ----------------------------------------------
c     disk intrinsic flux F_d in units of m_pc^3/\sigma_T R_s
      function xsxfd(radi)

      implicit none
      real xsxfd, radi

      xsxfd = (3./8.) * (1. - sqrt(3./radi))/
     $     (0.06 * radi**3)

      return
      end

       

c --------------------------------------------------------      
      
      SUBROUTINE schw(ear, ne, param, nr, rd, fra, photar)

      implicit none
      INTEGER ne
      REAL ear(ne), param(5), photar(ne)

c  model to calculate the line shape for a rotating accretion
c  disk. does not include GR effects. note that if param(2) is
c  set to 10 then do the special case of an accretion disk
c  emissivity.

c  parameters :
c	1        line energy
c       2        power law index for emissivity (10 for disk)
c       3        inner radius (GM/c**2)
c       4        outer radius (GM/c**2)
c       5        inclination  (degrees)

      INTEGER nr
c      PARAMETER (nr=100)

      INTEGER n, k, ie
      REAL dera,
     &     tgcsi2, cosdel, fra(nr), ra, dra,
     &     rd(nr), enel, spm
      REAL beta, betal, betah, zpo, zpol, zpoh
      REAL sincl, cincl, rafact, flux, total
      REAL enobl, enobh

      DATA dera /57.295779/

c convert input parameters - inclination set to radians.

      enel = param(1)

      sincl = sin(param(5)/dera)
      cincl = cos(param(5)/dera)

c set spectrum to 0

      DO n = 1, ne
         photar(n) = 0.
      ENDDO

c calculate log step in radius

c big loop for radii. note that the radius cannot reach 3 else
c the metric goes singular

      DO n = 1, nr - 1

         ra = (rd(n)+rd(n+1))/2.

         dra = rd(n+1) - rd(n)

         rafact = sqrt(1.-3./ra)

c if power-law index is less than ten use to calculate emissivity

c loop over azimuthal angles

         DO k = 1, 179, 4

            beta = float(k)/dera
            betal = float(k-1)/dera
            betah = float(k+1)/dera

c calculate mean redshift (1+z = zpo) for the bin

            tgcsi2 = (sincl*sin(beta))
     &               **2/(1.-(sincl*sin(beta))**2)
            cosdel = sincl*cos(beta)
     &               /sqrt((sincl*cos(beta))**2+cincl**2)
            zpo = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact

c and the low and high redshifts for the bin. note the traps for
c the case of an inclination of 90 degrees.


            IF ( param(5) .GT. 89.9 .AND. k .EQ. 91 ) THEN

               zpol = 1./rafact

            ELSE

               tgcsi2 = (sincl*sin(betal))
     &                  **2/(1.-(sincl*sin(betal))**2)
               cosdel = sincl*cos(betal)
     &                  /sqrt((sincl*cos(betal))**2+cincl**2)
               zpol = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact

            ENDIF

            IF ( param(5) .GT. 89.9 .AND. k .EQ. 89 ) THEN

               zpoh = 1./rafact

            ELSE

               tgcsi2 = (sincl*sin(betah))
     &                  **2/(1.-(sincl*sin(betah))**2)
               cosdel = sincl*cos(betah)
     &                  /sqrt((sincl*cos(betah))**2+cincl**2)
               zpoh = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact

            ENDIF

c  enobl and enobh are the lower and upper observed energy from
c  this azimuthal and radial bin

            enobl = min(enel/zpol, enel/zpoh)
            enobh = max(enel/zpol, enel/zpoh)
            total = enobh - enobl

c  calculate total emission from this bin

            flux = fra(n)*ra*dra*4./dera/zpo**3

c  find fractions of emission from this bin to place in each energy
c  range.

            IF (enobh .LT. ear(1) .OR. enobl .GT. ear(ne) ) GOTO 10

            DO ie = 1, ne
               IF ( ear(ie) .GE. enobh ) THEN
                  IF ( ear(ie-1) .LE. enobl ) THEN
                     photar(ie) = photar(ie) + flux
                     GOTO 10
                  ELSEIF ( ear(ie-1) .GE. enobh ) THEN
                     GOTO 10
                  ELSE
                     IF ( total .GT. 0. ) THEN
                        photar(ie) = photar(ie) 
     &                     + flux*(enobh-ear(ie-1))/total
                     ENDIF
                  ENDIF
               ELSEIF ( ear(ie) .GE. enobl ) THEN
                  IF ( ear(ie-1) .GE. enobl ) THEN
                     IF ( total .GT. 0. ) THEN
                        photar(ie) = photar(ie) 
     &                     + flux*(ear(ie)-ear(ie-1))/total
                     ENDIF
                  ELSE
                     IF ( total .GT. 0. ) THEN
                        photar(ie) = photar(ie) 
     &                     + flux*(ear(ie)-enobl)/total
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

 10         CONTINUE

         ENDDO
      ENDDO

c normalise values to total

      spm = 0.
      DO n = 1, ne-1
         spm = spm + photar(n) * (ear(n+1) - ear(n))
      ENDDO
c      spm = spm * (ear(10) - ear(9))

c write values

      IF ( spm .NE. 0. ) THEN
         DO n = 1, ne
            photar(n) = photar(n)/spm
         ENDDO
      ENDIF


      RETURN
      END


c --------------------------------------------------------      
      
      SUBROUTINE schw_tot(ear, ne, eart, phottmp, ne2, enel, photar, 
     1     norm, icl)

      implicit none
      INTEGER ne, ne2
      REAL ear(0:ne), photar(ne), eart(ne2), phottmp(ne2)

      INTEGER n, k, ie, imin, imax, icl
      REAL enel
      REAL flux, total
      REAL enobl, enobh, norm

c set spectrum to 0
      DO n = 1, ne
         photar(n) = 0.
      ENDDO

      DO k = 1, ne2-1

c  enobl and enobh are the lower and upper observed energy from
c  this azimuthal and radial bin

         enobl = enel*eart(k)
         enobh = enel*eart(k+1)
         total = enobh - enobl
         
c     calculate total emission from this bin

         flux = 0.5*(phottmp(k) + phottmp(k+1))

c     find fractions of emission from this bin to place in each energy
c     range.
         IF (enobh .LT. ear(1) .OR. enobl .GT. ear(ne) ) GOTO 50

         imin = max(1, icl-80)
         imax = min(ne, icl+60)

         DO ie = imin, imax
            IF ( ear(ie) .GE. enobh ) THEN
               IF ( ear(ie-1) .LE. enobl ) THEN
                  photar(ie) = photar(ie) + flux
                  GOTO 50
               ELSEIF ( ear(ie-1) .GE. enobh ) THEN
                  GOTO 50
               ELSE
                  IF ( total .GT. 0. ) THEN
                     photar(ie) = photar(ie) 
     &                    + flux*(enobh-ear(ie-1))/total
                  ENDIF
               ENDIF
            ELSEIF ( ear(ie) .GE. enobl ) THEN
               IF ( ear(ie-1) .GE. enobl ) THEN
                  IF ( total .GT. 0. ) THEN
                     photar(ie) = photar(ie) 
     &                    + flux*(ear(ie)-ear(ie-1))/total
                  ENDIF
               ELSE
                  IF ( total .GT. 0. ) THEN
                     photar(ie) = photar(ie) 
     &                    + flux*(ear(ie)-enobl)/total
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         
 50      CONTINUE
         
      ENDDO

      norm = 0.
      DO ie = 1, ne-1
         norm = norm + photar(ie) * (ear(ie+1)-ear(ie))
      enddo

      RETURN
      END
