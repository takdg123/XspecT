!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine jet(ear,ne,param,ifl,photar)

      implicit none
      integer ne,nn,n,i,ix,ifl
      parameter (nn=3000)
      real ear(0:ne),photar(ne),param(*)
      real e1,e2,emin,emax,zfac
      real z4,M4,rco4,mdot4,a_star4,thetaobs4,BLF4,Rdiss4,phi4,B4,gmax4
      real logPrel4,gmin_inj4,gb4,s14,s24,eobs4(nn),Feobs4(nn),logPj4
      save eobs4,Feobs4

c unused variables
c      real alphax4,exc4,logL_BLR4,logL_IR4,logL_X4,RBLR4,Rdin4
c      real Rdout4,RIR4,Rx4,T_IR4


c     param(1) mass in solar
c     param(2) comoving distance in Mpc
c     param(3) log mass accretion rate in L/Ledd
c     param(4) inclination in degrees
c     param(5) BLF
c     param(6) phi - opening angle in radians of jet rdiss=phi*zdiss
c     param(7) zdiss in Rg 
c     param(8) B field in G
c     param(9) log(Prel) in erg s-1
c     param(10) gmin_inj
c     param(11) gamma_b
c     param(12) gmax_inj
c     param(13) s1 - index below the break
c     param(14) s2 - index above the break
c     param(15) z

c     Extra files written out (write statements currently commented out): 
c     fort.2 contains luminosity of individual components
c     (log(v)  log(vL_sync(v))  log(vL_ssc(v))  log(vL_exc(v))  log(vL_tot(v)))
c     fort.3 contains seed photon energy densities
c     (log(v)  log(vU_disc(v))  log(vU_corona(v))  log(vU_BLR(v))  log(vU_reflcorona(v))
c      log(vU_torus(v))  log(vU_externaltotal(v))  log(vU_sync(v)))

c     silly line to suppress a compiler warning
      i = ifl
      
      z4=param(15) 
      m4=param(1) 
      rco4=param(2)
      mdot4=10.0**(param(3)) 
!      a_star4=param(4) 
      a_star4=0.0
      thetaobs4=param(4)
      blf4=param(5) 
      rdiss4=param(6)*param(7)*(1.48222e5*param(1)) !Get rdiss from zdiss
      phi4=param(6) 
      b4=param(8) 
      logprel4=param(9) 
      gmin_inj4=param(10) 
      gmax4=param(12) 
      gb4=param(11) 
      s14=param(13) 
      s24=param(14)
      zfac=1.0+param(15)

      call jetspec(z4,m4,rco4,mdot4,a_star4,thetaobs4,blf4,rdiss4,phi4,
     &                         b4,logprel4,gmin_inj4,gmax4,gb4,s14,s24,
     &                                             eobs4,Feobs4,logPj4)

c     now redshift the energy bins back
      do n=1,nn,1
         eobs4(n)=eobs4(n)/zfac
      end do

c     and now correct the flux
      do n=1,nn,1
         Feobs4(n)=Feobs4(n)/zfac
      end do

c     rebin onto original grid
      ix = 1
      do i = 1, ne
         photar(i) = 0.0
         emin = ear(i - 1)
         emax = ear(i)
         do while(eobs4(ix) < emin)
           ix = ix + 1
         end do
         do while(eobs4(ix-1) < emax)
           e1 = max(eobs4(ix - 1), emin)
           e2 = min(eobs4(ix), emax)
           photar(i) = photar(i) + 
     &                      Feobs4(ix)*(e2-e1)/(eobs4(ix)-eobs4(ix-1))
           ix = ix + 1
         end do
         ix = ix - 1
      end do

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine jetspec(z4,M4,rco4,mdot4,a_star4,thetaobs4,BLF4,Rdiss4,
     &                     phi4,B4,logPrel4,gmin_inj4,gmax4,gb4,s14,s24,
     &                                              eobs4,Feobs4,logPj4)
c     Calculates jet spectrum in jet frame and then converts to 
c     observer frame.
c     Includes synchrotron and synchrotron-self-Compton.
c     Includes inverse-Compton from external seed photons if mdot>=0.01. 
c     Calculates self-consistent electron distribution. 
c     go = 1 
c     gmaxo = 1.0e6 (!gmax must be smaller than this!)
c     gmin = minimum gamma of steady state electron distribution
c     gmin_inj = minimum gamma of injected electron distribution
c     gmax = maximum gamma of injected, and therefore steady state,
c                                            electron distributions.
c     Includes cooling.
c     Assumes thetaobs is in radians.

      implicit none
      integer i,j,nv,ng,imin_sync,imax_sync,gcoolindex,loop
      integer gmin_injindex,gmin_index,gmax_index
      parameter (nv=3000)
      parameter (ng=1716)
      double precision pi,c,me,sigmat,BLF,B,Rdiss,s1,s2,gb,thetaobs
      double precision deltaobs,K,gmin,gmax,dlogg,g(ng),n(ng),vssa
      double precision nn(ng),dn,gmin_inj,go,A,z,logvmax,logvmin,dlogv
      double precision v(nv),L(nv),vobs(nv),Lobs(nv),dv(nv),js(nv)
      double precision jc(nv),Ub,Uvexseed_j(nv),M,mdot,Pb,Prel,gcool
      double precision kevhz,sum,dg(ng),q(ng),Qo,Uvsync(nv),Useed,vb
      double precision Pe,Pp,Pj,Pr,lprime,zdiss,phi,Rg,Lsobs(nv)
      double precision Lexcobs(nv),Mpc,Fobs(nv),integral,RBLR,Ldisc,RIR
      double precision a_star,vpeak,vfvpeak,gmaxo,UvsyncKN(nv),mp,h
      double precision dvobs(nv),KNindex(nv),UsyncKN,UexKN,Uv(nv,9)
      double precision Lsscobs(nv)
      logical repeat
      real z4,M4,mdot4,a_star4,thetaobs4,BLF4,Rdiss4,phi4,B4,gmax4
      real logPrel4,gmin_inj4,gb4,s14,s24,eobs4(nv),Feobs4(nv),logPj4
      real rco4
      double precision jssc(nv),jexc(nv),Uzeros(nv)
      double precision Lcomp,Rco
      double precision jexcd(nv),jexcc(nv),jexcb(nv)
      double precision jexcr(nv),jexct(nv),Lexcdobs(nv),Lexccobs(nv)
      double precision Lexcbobs(nv),Lexcrobs(nv),Lexctobs(nv)
      double precision denom

c unused variables
c      integer num
c      real alphax4,exc4,logL_BLR4,logL_IR4,logL_X4,RBLR4,,Rdin4
c      real Rdout4,RIR4,Rx4,T_IR4
c      double precision alphax,exc,dve(nv),L_BLR,L_IR,L_X,Rdisc,Rin,Rx
c      double precision ve(nv),T_IR

      pi = 6.0*asin(0.5)

c     Convert from reals to double precision
      z = dble(z4)
      M = dble(m4)
      Rco = dble(rco4)
      mdot = dble(mdot4) 
      a_star = dble(a_star4)
      thetaobs = dble(thetaobs4*pi/180.0) !convert to radians
      BLF = dble(BLF4) 
      Rdiss = dble(Rdiss4)
      phi = dble(phi4)
      B = dble(B4)
      Prel = 10.0**(dble(logprel4)) 
      gmin_inj = dble(gmin_inj4)
      gmax = dble(gmax4)
      gb = dble(gb4)
      s1 = dble(s14)
      s2 = dble(s24)
      !Rin = dble(Rdin4)
      !Rdisc = dble(Rdout4)
      !Rx = dble(Rx4)
      !L_X = 10.0**dble(logL_X4)
      !alphax = dble(alphax4)
      !exc = dble(exc4)
      !L_BLR = 10.0**dble(logL_BLR4)
      !RBLR = dble(RBLR4)
      !L_IR = 10.0**dble(logL_IR4)
      !RIR = dble(RIR4)
      !T_IR = dble(T_IR4)

      c = 3.0e10
      h = 6.63e-27
      me = 9.1e-28
      mp = 1.67e-24
      sigmat = 6.65e-25
      kevhz = 1000.0*1.6e-19/6.63e-34
      Rg = (6.67e-8*M*2.0e33)/(3.0e10)**2
      Mpc = 1.0e6*3.0e18 !in cm

      deltaobs = 1.0/(BLF-cos(thetaobs)*sqrt(BLF**2-1.0))
      Ub = B**2.0/(8.0*pi)
      vb = 2.8e6*B
      write(*,*) 'delta',deltaobs

c     Make grid of frequencies
      logvmax = 28.0
      logvmin = 7.0
      dlogv = (logvmax-logvmin)/real(nv)
      do i=1,nv,1
         vobs(i) = 10.0**(logvmin+real(i)*dlogv-dlogv/2.0)
         dvobs(i) = 10.0**(log10(vobs(i))+dlogv/2.0)-
     &                    10.0**(log10(vobs(i))-dlogv/2.0)
         v(i) = vobs(i)/deltaobs
         dv(i) = 10.0**(log10(v(i))+dlogv/2.0)-
     &              10.0**(log10(v(i))-dlogv/2.0)
c         ve(i) = vobs(i)/deltaobs**2.0
c         dve(i) = 10.0**(log10(ve(i))+dlogv/2.0)-
c     &              10.0**(log10(ve(i))-dlogv/2.0)
         Uzeros(i) = 0.0
      end do

c     Make grid of gammas
      go = 1.0
      gmaxo = 1.0e6
      dlogg = (log10(gmaxo)-log10(go))/real(ng)
      do i=1,ng,1
         g(i) = 10.0**(log10(go)+real(i)*dlogg-dlogg/2.0)
         dg(i) = 10.0**(log10(g(i))+dlogg/2.0)-
     &              10.0**(log10(g(i))-dlogg/2.0)
      end do

c     Make grid of injected electron distribution
c     q(gamma) = (gamma/gammab)^(-s1)/(1+(gamma/gammab)^(-s1+s2))
      sum = 0.0
      gmin_injindex = int((log10(gmin_inj)-log10(go))/dlogg)+2
      if (gmin_inj.eq.go) then
       gmin_injindex = 1
      end if
      gmin_index = gmin_injindex
      gmax_index = int((log10(gmax)-log10(go))/dlogg)+2
      if (gmax_index.gt.ng) then
        gmax_index = ng
      end if
      do i=1,ng,1
       if (i.lt.gmin_injindex) then
         q(i) = 0.0
       end if
       if (i.eq.gmin_injindex) then
         q(i) = (g(i)/gb)**(-s1)/(1.0+(g(i)/gb)**(-s1+s2))
         sum = sum+g(i)*q(i)*
     &          (10.0**(log10(g(i))+dlogg/2.0)-gmin_inj)
       end if
       if ((i.gt.gmin_injindex).and.(i.lt.gmax_index)) then
         q(i) = (g(i)/gb)**(-s1)/(1.0+(g(i)/gb)**(-s1+s2))
         sum = sum+g(i)*q(i)*dg(i)
       end if
       if (i.eq.gmax_index) then
         q(i) = (g(i)/gb)**(-s1)/(1.0+(g(i)/gb)**(-s1+s2))
         sum = sum+g(i)*q(i)*
     &          (10.0**(log10(g(i))+dlogg/2.0)-gmax)
       end if
       if (i.gt.gmax_index) then
         q(i) = 0.0
       end if
      end do
c     Calculate Qo
      Qo = Prel/(me*c**2.0*4.0/3.0*pi*Rdiss**3.0*sum)
c     Normalise Q(gamma) = Qo q(gamma) and guess at N(gamma)
      do i=1,ng,1
         Q(i) = Qo*q(i)
         n(i) = Q(i)*(Rdiss/c)
      end do
      K = Qo*(Rdiss/c)
      gmin = gmin_inj
      Pb = 4.0*pi*Rdiss**2.0*c*Ub

c     Get external seed photons
      zdiss = Rdiss/(phi*Rg)
      Ldisc = mdot*1.38e38*M
      RBLR = 1.0e17*sqrt(Ldisc/1.0d45)/Rg !Default value
      RIR = 2.5e18*sqrt(Ldisc/1.0d45)/Rg !Default value
      if (mdot.ge.0.01) then
        write(*,*) 'Ldisc',Ldisc,'erg s^-1'
        write(*,*) 'RBLR',RBLR,'Rg',RBLR*Rg,'cm'
        write(*,*) 'RIR',RIR,'Rg',RIR*Rg,'cm'
        write(*,*) 'mdot>=0.01, so jet = SSC + EC'
        call SimpleSeedPhotons(M,a_star,BLF,zdiss,Ldisc,
     &                           RBLR,RIR,v,dv,Uvexseed_j,Uv)
      else
        write(*,*) 'mdot<0.01, so jet = SSC'
        do j=1,nv,1
          Uvexseed_j(j) = 0.0
        end do
      end if
      !write(7,*) log10(zdiss),log10(K),log10(B)

c     Check compactness for pair production
      do i=1,nv,1
        if (v(i).le.511.0*kevhz) then
          lprime = lprime+Uvexseed_j(i)*dv(i)
        end if
c        write(4,*) log10(v(i)),
c     &         log10(v(i)*4.0*pi*sigmat*Rdiss*Uvexseed_j(i)/(me*c**2.0))
      end do
      lprime = 4.0*pi*sigmat*Rdiss*lprime/(me*c**2.0)
      !write(*,*) lprime

c     Loop for calculating steady state electron distribution and spectrum
      repeat = .true.
      loop = 1
      do while(repeat)

        !Zero spectra and energy densities
        do i=1,nv,1
          js(i) = 0.0
          jc(i) = 0.0
          Uvsync(i) = 0.0
          UvsyncKN(i) = 0.0
        end do

!       Calculate vssa - this is hardwired for alpha0=1
!       Ghisellini maraschi & Treves 1985A&A...146..204
        vssa = (2.0*2.32e14*K*B**2.5*Rdiss/0.7)**(2.0/7.0)

!        write(*,*) vssa,k,b,rdiss       
!        write(*,*) 'dermer ',2.65e7*(rdiss/1e15 *K * B)**(1.0/3.0)

        !Find imin_sync and approx imax_sync
        imin_sync = int((log10(vssa)-log10(v(1)))/dlogv)+2
        if (imin_sync.lt.0) then
          imin_sync = 1
        end if
!        write(*,*) vssa,v(imin_sync)
        do i=1,nv,1
          if (v(i).gt.(4.0/3.0*vb*gmax**2.0)) then
            imax_sync = i-1
            exit
          end if
        end do

        !Do synchrotron from vssa to vmax_sync with Delta Fn Approximation
        call DeltaSync(ng,g,n,dlogg,gmin_index,gmax_index,v,vb,Ub,
     &                              imin_sync,imax_sync,Rdiss,js,Uvsync)

        !Do compton from each synchrotron vseed
        call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                            Uvsync,Uvexseed_j,jc,UvsyncKN,KNindex)

        !Calculate radiation energy density in jet and gcool
        Useed = Ub
        do i=1,nv,1
          Useed = Useed+UvsyncKN(i)*dv(i)
        end do 
        gcool = 3.0*me*c**2.0/(4.0*sigmat*Rdiss*Useed)
        gcoolindex = int((log10(gcool)-log10(go))/dlogg)+2
        if (gcoolindex.lt.1) then
          gcoolindex = 1
        end if
        if (gcool.lt.1.0) then
          gcool = 1.0
        end if

        if (gcoolindex.gt.gmax_index) then
          !No cooling so steady state distribution same as injected
          !write(*,*) 'gmax < gcool'
          gmin = gmin_inj
          gmin_index = gmin_injindex
          do i=1,gmin_injindex-1,1
            nn(i) = 0.0
          end do
          do i=gmin_injindex,gmax_index,1
            nn(i) = n(i)
          end do
          if (gmax_index.lt.ng) then
            do i=gmax_index+1,ng,1
              nn(i) = 0.0
            end do
          end if
          integral = 0.0
        end if

        if ((gcoolindex.le.gmax_index).and.
     &                              (gcoolindex.ge.gmin_injindex)) then
          !Partial cooling
          !write(*,*) 'gmin_inj < gcool < gmax'
          if (gcoolindex.gt.gmax_index) then
            gcoolindex = gmax_index
          end if
          !Recalculate K
          sum = 0.0
          do i=gcoolindex,gmax_index,1
            sum = sum+Q(i)*dg(i)
          end do
          A = 3.0*me*c/(4.0*sigmat*Useed)*sum/
     &                               (g(gcoolindex)**2.0*Q(gcoolindex))
          K = A*Qo
          gmin = gmin_inj
          gmin_index = gmin_injindex
          integral = sum/(Q(gcoolindex))
          !Recalculate electron distribution
          do i=1,gmin_injindex-1,1
             nn(i) = 0.0
          end do
          do i=gmin_injindex,gmax_index,1
           if (i.le.gcoolindex) then
             nn(i) = A*Q(i)
           else if (i.gt.gcoolindex) then
             nn(i) = 0.0
             do j=i,gmax_index,1
               nn(i) = nn(i)+Q(j)*dg(j)
             end do
             nn(i) = nn(i)*3.0*me*c**2.0/(4.0*sigmat*c*Useed*g(i)**2.0)
           end if
          end do
          if (gmax_index.lt.ng) then
            do i=gmax_index+1,ng,1
               nn(i) = 0.0
            end do
          end if
        end if

        if (gcoolindex.lt.gmin_injindex) then
          !Complete cooling
          !write(*,*) 'gcool < gmin_inj'
          !Recalculate K
          sum = 0.0
          do i=gmin_injindex,gmax_index,1
            sum = sum+Q(i)*dg(i)
          end do
          K = 3.0*me*c/(4.0*sigmat*Useed)*sum
          gmin = gcool
          gmin_index = gcoolindex
          integral = 0.0
          !Recalculate electron distribution
          do i=1,ng,1
           if (i.lt.gcoolindex) then
             nn(i) = 0.0
           end if
           if ((i.ge.gcoolindex).and.(i.le.gmax_index)) then
             nn(i) = 0.0
             do j=i,gmax_index,1
               nn(i) = nn(i)+Q(j)*dg(j)
             end do
             nn(i) = nn(i)*3.0*me*c**2.0/(4.0*sigmat*c*Useed*g(i)**2.0)
           end if
           if (i.gt.gmax_index) then
             nn(i) = 0.0
           end if
          end do
        end if

        !Check how much electron distribution has changed and loop again
        !or finish
        do i=1,gmax_index,1
          repeat = .false.
          dn = abs((nn(i)-n(i))/n(i))
          if ((i.ge.gmin_index).and.(dn.gt.0.1)) then
            repeat = .true.
            do j=1,ng,1
             if ((j.lt.gmin_index).or.(j.gt.gmax_index)) then
              n(j) = 0.0
             else
              n(j) = (nn(j)+n(j))/2.0
             end if
            end do
            exit
          end if
        end do
        !write(*,*) loop,repeat,sngl(gcool)
        loop = loop+1
      end do
      write(*,*) 'gcool',gcool

      UsyncKN = 0.0
      UexKN = 0.0
      !write(3,*) 'log(v) log(vUvsync)'
      do i=1,nv,1
         Uv(i,5) = Uvsync(i)
         !write(num,*) v(i),(v(i)*Uv(i,j),j=1,5)
         if (Uvsync(i).gt.0.0) then
           !write(3,*) log10(v(i)),log10(v(i)*Uvsync(i))
         end if
         if (KNindex(i).gt.0.0) then
            UsyncKN = UsyncKN+Uvsync(i)*dv(i)
            UexKN = UexKN+Uvexseed_j(i)*dv(i)
         end if 
      end do
      !write(8,*) log10(M),log10(Ub),log10(UsyncKN),log10(UexKN)
      !write(3,*) 'no no'

c     Get just EC
      if (mdot.ge.0.01) then
        call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                          Uzeros,Uvexseed_j,jexc,UvsyncKN,KNindex)
      else
        do i=1,nv,1
          jexc(i) = 0.0
        end do
      end if

c     Get just disc EC
      do i=1,nv,1
        Uvexseed_j(i) = Uv(i,6)
      end do
      if (mdot.ge.0.01) then
        call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                         Uzeros,Uvexseed_j,jexcd,UvsyncKN,KNindex)
      else
        do i=1,nv,1
          jexc(i) = 0.0
        end do
      end if

c     Get just corona EC
      do i=1,nv,1
        Uvexseed_j(i) = Uv(i,7)
      end do
      if (mdot.ge.0.01) then
        call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                         Uzeros,Uvexseed_j,jexcc,UvsyncKN,KNindex)
      else
        do i=1,nv,1
          jexc(i) = 0.0
        end do
      end if

c     Get just BLR EC
      do i=1,nv,1
        Uvexseed_j(i) = Uv(i,8)
      end do
      if (mdot.ge.0.01) then
        call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                         Uzeros,Uvexseed_j,jexcb,UvsyncKN,KNindex)
      else
        do i=1,nv,1
          jexc(i) = 0.0
        end do
      end if

c     Get just reflection of coronal X-rays EC
      do i=1,nv,1
        Uvexseed_j(i) = Uv(i,9)
      end do
      if (mdot.ge.0.01) then
        call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                         Uzeros,Uvexseed_j,jexcr,UvsyncKN,KNindex)
      else
        do i=1,nv,1
          jexc(i) = 0.0
        end do
      end if

c     Get just torus EC
      do i=1,nv,1
        Uvexseed_j(i) = Uv(i,3)
      end do
      if (mdot.ge.0.01) then
        call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                         Uzeros,Uvexseed_j,jexct,UvsyncKN,KNindex)
      else
        do i=1,nv,1
          jexc(i) = 0.0
        end do
      end if

c     Get just SSC
      do i=1,nv,1
        Uvexseed_j(i) = 0.0
      end do
      call Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                          Uvsync,Uvexseed_j,jssc,UvsyncKN,KNindex)

c     Sum up
      Pr = 0.0
      Lcomp = 0.0
      denom = 3.0*pi*Rdiss**3.0*deltaobs**3.0
      do i=1,nv,1
         L(i) = 4.0*pi*(js(i)+jc(i))*4.0/3.0*pi*Rdiss**3.0
         Lobs(i) = L(i)*deltaobs**3.0
         Fobs(i) = (js(i)*deltaobs**3.0+jssc(i)*deltaobs**3.0+
     &                         jexc(i)*deltaobs**3.0)
     &                     *4.0/3.0*pi*Rdiss**3.0/(Rco*Mpc)**2.0
         Lsobs(i) = 4.0*pi*js(i)*4.0/denom
         Lsscobs(i) = 4.0*pi*jssc(i)*4.0/denom
         Lexcobs(i) = 4.0*pi*jexc(i)*4.0/denom
         Lexcdobs(i) = 4.0*pi*jexcd(i)*4.0/denom
         Lexccobs(i) = 4.0*pi*jexcc(i)*4.0/denom
         Lexcbobs(i) = 4.0*pi*jexcb(i)*4.0/denom
         Lexcrobs(i) = 4.0*pi*jexcr(i)*4.0/denom
         Lexctobs(i) = 4.0*pi*jexct(i)*4.0/denom
         Pr = Pr+Lobs(i)*dvobs(i)/(4.0*deltaobs**2.0)
         Lcomp = Lcomp+Lexcobs(i)*dvobs(i)
         !Uncomment this to write out the individual jet spectral components  
c         write(2,*) log10(vobs(i)),
c     &              log10(vobs(i)*Lsobs(i)),log10(vobs(i)*Lsscobs(i)),
c     &            log10(vobs(i)*Lexcobs(i)),log10(vobs(i)*Lexcdobs(i)),
c     &            log10(vobs(i)*Lexccobs(i)),log10(vobs(i)*Lexcbobs(i)),
c     &            log10(vobs(i)*Lexcrobs(i)),log10(vobs(i)*Lexctobs(i)),
c     &            log10(vobs(i)*Lobs(i))
      end do

c     Find v_sync,peak
      vpeak = vobs(imin_sync)
      vfvpeak = vobs(imin_sync)*js(imin_sync)
      do i=imin_sync+1,imax_sync,1
        if (vobs(i)*js(i).gt.vfvpeak) then
          vfvpeak = vobs(i)*js(i)
          vpeak = vobs(i)
        end if
      end do
      !write(*,*) 'vpeak',log10(vpeak),log10(vfvpeak)

c     Calculate jet power
      Pe = 0.0
      Pp = 0.0
      do i=1,ng,1
        Pe = Pe+g(i)*N(i)*dg(i)
        Pp = Pp+N(i)*dg(i)
c       write(908,*) g(i),n(i)
      end do
c      write(908,*) 'no no'
      Pe = pi*Rdiss**2.0*BLF**2.0*sqrt(1.0-1.0/BLF**2.0)*c*Pe*me*c**2.0
      Pp = pi*Rdiss**2.0*BLF**2.0*sqrt(1.0-1.0/BLF**2.0)*c*Pp*mp*c**2.0
      Pb = pi*Rdiss**2.0*BLF**2.0*sqrt(1.0-1.0/BLF**2.0)*c*Ub
      Pj = Pr+Pb+Pp+Pe
      write(*,*) 'log Pr',sngl(log10(Pr)),'log Pb',sngl(log10(Pb))
      write(*,*) 'log Pe',sngl(log10(Pe)),'log Pp',sngl(log10(Pp))
      write(*,*) 'log Pj',sngl(log10(Pj))

c     Convert from double precision to reals and xspec units
      do i=1,nv,1
        eobs4(i) = sngl(vobs(i)*h*511.0/(me*c**2.0))
        Feobs4(i) = sngl(Fobs(i)/(h*vobs(i))*dvobs(i))
      end do
      logPj4 = sngl(log10(Pj))

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SimpleSeedPhotons(M,a,BLF,zdiss,Ldisc,RBLR,RIR,
     &                                              v,dv,Uvexseed_j,Uv)

      implicit none
      integer i,j,nvtot,nphi,nv,nmu
      parameter (nvtot=3000)
      parameter (nv=1000)
      parameter (nphi=500)
      parameter (nmu=1000)
      double precision v(nvtot),dv(nvtot),dlogv,pi,c,Rg,M,zdiss
      double precision BLF,UBLR,LBLR,RBLR
      double precision f_BLR,norm_BLR
      double precision Ldisc
      double precision v_lyalpha,Lvexseed_j(nvtot),norm_IR
      double precision Uvexseed_j(nvtot),mu1,mu2,beta,f_IR,UIR,mu,v_IR
      double precision LIR,RIR,Rd,mud(nmu),mumin,mumax,Rdisc,Ud(nvtot)
      double precision b,dmu,Id,h,k,a,T,alphax,fx,Lx,norm_x,mux,Rx
      double precision vxc,Ux,vdpeak,fXBLR,UXBLR,LXBLR,norm_XBLR,Udtot
      double precision const,grad,U1,U2,isco,Rin,Uv(nvtot,9),Tmax

c unused variables
c      integer index
c      double precision d,dA,deltaBLR,dlogv_BLR,dphi,dv_lyalpha
c      double precision dvBLR(nv),exc,fL,fLtot,L_BLR,L_IR,L_X
c      double precision logvmax_BLR,logvmin_BLR,LvBLR(nv)
c      double precision phi(nphi),phimax,phimaxlos,phimin,r,T_IR
c      double precision thetaBLR(nphi),Uex(nv),vBLR(nv)

      pi = 6.0*asin(0.5)
      c = 3.0e10
      h = 6.63e-27
      k = 1.36e-16
      Rg = (6.67e-8*M*2.0e33)/(3.0e10)**2
      beta = sqrt(1.0-1.0/BLF**2.0)

      do i=1,nvtot,1
       Lvexseed_j(i) = 0.0
       Uvexseed_j(i) = 0.0
       do j=1,9,1
         Uv(i,j) = 0.0
       end do
      end do
      dlogv = log10(v(2))-log10(v(1))
      !write(3,*) 'skip on'

c     Do disc spectrum in jet frame
      Rdisc = 2.0*500.0 !Default value
      Rin = isco(a) !Default value
      mumin = (1.0+(Rdisc/zdiss)**2.0)**(-1.0/2.0)
      mumax = 1.0 
      do j=1,nmu,1
       mud(j) = mumin+real(j)*(mumax-mumin)/real(nmu-1)
      end do
      !write(3,*) 'log(v) log(vUd)'
      Udtot = 0.0
      Tmax = 0.0
      do i=1,nvtot,1
       Ud(i) = 0.0
       do j=1,nmu-1,1
        dmu = mud(j+1)-mud(j)
        b = BLF*(1.0-beta*mud(j))
        Rd = zdiss*tan(acos(mud(j)))
        if (Rd.gt.Rin) then
         T = (6.0*Rg*Ldisc/(16.0*pi*0.08*5.67e-5*(Rd*Rg)**3.0)*
     &                               (1.0-sqrt(6.0/Rd)))**(1.0/4.0)
         Id = 2.0*h/c**2.0*(v(i)/b)**3.0/(exp(h*(v(i)/b)/(k*T))-1.0)
         Ud(i) = Ud(i)+2.0*pi/c*b*Id*dmu
         if (T.gt.Tmax) then
          Tmax = T
         end if
        end if
       end do
       Udtot = Udtot+Ud(i)*dv(i)
       Uvexseed_j(i) = Uvexseed_j(i)+Ud(i)
       Uv(i,6) = Uv(i,6)+Ud(i)
       if (Ud(i).gt.0.0) then
        !write(3,*) log10(v(i)),log10(v(i)*Ud(i))
       end if
      end do
      !write(3,*) 'no no'
      vdpeak = 4.0*k*Tmax/h

c     Do corona spectrum in jet frame
      Rx = 2.0*30.0 !Default value
      fx = 0.1 !Default value
      mux = (1.0+(Rx/zdiss)**2.0)**(-1.0/2.0)
      b = BLF*(1.0-beta*mux)
      Ux = fx*Ldisc*BLF**2.0/(pi*c*(Rx*Rg)**2.0)*
     &       (1.0-mux-beta*(1.0-mux**2.0)+beta**2.0/3.0*(1.0-mux**3.0))
      vxc = 150.0e3*(1.6e-12)/h !Default value
      alphax = 1.0 !Default value
      !write(3,*) 'log(v) log(vUx)'
      Lx = 0.0
      do i=1,nvtot,1
       Lvexseed_j(i) = 0.0
       if (v(i).gt.vdpeak*b) then
        Lvexseed_j(i) = v(i)**(-alphax)*exp(-v(i)/(b*vxc))
        Lx = Lx+Lvexseed_j(i)*dv(i)
       end if
      end do
      norm_x = Ux/Lx
      do i=1,nvtot,1
       if (v(i).gt.vdpeak*b) then
        Uvexseed_j(i) = Uvexseed_j(i)+norm_x*Lvexseed_j(i)
        Uv(i,7) = norm_x*Lvexseed_j(i)
        if (Lvexseed_j(i).gt.0.0) then
         !write(3,*) log10(v(i)),log10(v(i)*norm_x*Lvexseed_j(i))
        end if
       end if
       Uv(i,1) = Uvexseed_j(i)
      end do
      !write(3,*) 'no no'

c     Do Broad Line Region Spectrum in jet frame
      f_BLR = 0.1 !Default value
      !write(3,*) 'log(v) log(vUblr)'
      if (zdiss.lt.RBLR) then
       UBLR = (17.0*BLF**2.0/12.0)*f_BLR*Ldisc/(4.0*pi*c*(RBLR*Rg)**2.0)
       mu1 = (1.0+(RBLR/zdiss)**2.0)**(-1.0/2.0)
       b = BLF 
       !write(*,*) mu1,beta,BLF,b
      end if
      if ((RBLR.le.zdiss).and.(zdiss.le.3.0*RBLR)) then
       write(*,*) '! zdiss > R_BLR !'
       U1 = (17.0*BLF**2.0/12.0)*f_BLR*Ldisc/(4.0*pi*c*(RBLR*Rg)**2.0)
       mu1 = (1.0+(1.0/3.0)**2.0)**(-1.0/2.0)
       mu2 = (1.0-(1.0/3.0)**2.0)**(1.0/2.0)
       U2 = f_BLR*Ldisc/(4.0*pi*c*(RBLR*Rg)**2.0)*BLF**2.0/(3.0*beta)*
     &     (2.0*(1.0-beta*mu1)**3.0-(1.0-beta*mu2)**3.0-(1.0-beta)**3.0)
       grad = (log10(U2)-log10(U1))/(log10(3.0*RBLR)-log10(RBLR))
       const = log10(U1)-grad*log10(RBLR)
       UBLR = 10.0**(grad*log10(zdiss)+const)
       b = BLF*(1.0-beta*mu1) 
      end if
      if (zdiss.gt.3.0*RBLR) then
       write(*,*) '! zdiss > 3R_BLR !'
       mu1 = (1.0+(RBLR/zdiss)**2.0)**(-1.0/2.0)
       mu2 = (1.0-(RBLR/zdiss)**2.0)**(1.0/2.0)
       UBLR = f_BLR*Ldisc/(4.0*pi*c*(RBLR*Rg)**2.0)*BLF**2.0/(3.0*beta)*
     &     (2.0*(1.0-beta*mu1)**3.0-(1.0-beta*mu2)**3.0-(1.0-beta)**3.0)
       b = BLF*(1.0-beta*mu1) 
      end if
      if (RBLR.ne.0.0) then
       v_lyalpha = (c/1216.0e-8)*b/4.0
       LBLR = 0.0
       do i=1,nvtot,1
         Lvexseed_j(i) = v(i)**3.0/(exp(v(i)/v_lyalpha)-1.0)
         LBLR = LBLR+Lvexseed_j(i)*dv(i)
       end do
       norm_BLR = (4.0*pi*c*(RBLR*Rg)**2.0*UBLR)/LBLR
       do i=1,nvtot,1
         Lvexseed_j(i) = norm_BLR*Lvexseed_j(i)
         Uvexseed_j(i) = Uvexseed_j(i)+
     &                       Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0)
         Uv(i,2) = Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0)
         Uv(i,8) = Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0)
         if (Lvexseed_j(i).gt.0.0) then
c          write(3,*) log10(v(i)),
c     &           log10(v(i)*Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0))
         end if
       end do
      end if
      !write(3,*) 'no no'

c     Do reflection of corona X-rays by BLR
      fXBLR = 0.01
      !write(3,*) 'log(v) log(vUxblr)'
      if (zdiss.lt.RBLR) then
       UXBLR = (17.0*BLF**2.0/12.0)*fXBLR*fX*Ldisc
     &                                   /(4.0*pi*c*(RBLR*Rg)**2.0)
       mu1 = (1.0+(RBLR/zdiss)**2.0)**(-1.0/2.0)
       b = BLF 
      end if
      if ((RBLR.le.zdiss).and.(zdiss.le.3.0*RBLR)) then
       U1 = (17.0*BLF**2.0/12.0)*fXBLR*fX*Ldisc
     &                                   /(4.0*pi*c*(RBLR*Rg)**2.0)
       mu1 = (1.0+(1.0/3.0)**2.0)**(-1.0/2.0)
       mu2 = (1.0-(1.0/3.0)**2.0)**(1.0/2.0)
       U2 = fXBLR*fX*Ldisc/(4.0*pi*c*(RBLR*Rg)**2.0)
     &          *BLF**2.0/(3.0*beta)*
     &     (2.0*(1.0-beta*mu1)**3.0-(1.0-beta*mu2)**3.0-(1.0-beta)**3.0)
       grad = (log10(U2)-log10(U1))/(log10(3.0*RBLR)-log10(RBLR))
       const = log10(U1)-grad*log10(RBLR)
       UXBLR = 10.0**(grad*log10(zdiss)+const)
       b = BLF*(1.0-beta*mu1) 
      end if
      if (zdiss.gt.3.0*RBLR) then
       mu1 = (1.0+(RBLR/zdiss)**2.0)**(-1.0/2.0)
       mu2 = (1.0-(RBLR/zdiss)**2.0)**(1.0/2.0)
       UXBLR = fXBLR*fX*Ldisc/(4.0*pi*c*(RBLR*Rg)**2.0)
     &          *BLF**2.0/(3.0*beta)*
     &     (2.0*(1.0-beta*mu1)**3.0-(1.0-beta*mu2)**3.0-(1.0-beta)**3.0)
       b = BLF*(1.0-beta*mu1) 
      end if
      if (RBLR.ne.0.0) then
       LXBLR = 0.0
       do i=1,nvtot,1
         Lvexseed_j(i) = 0.0
         if (v(i).gt.vdpeak*b) then
           Lvexseed_j(i) = v(i)**(-alphax)*exp(-v(i)/(b*vxc))
           LXBLR = LXBLR+Lvexseed_j(i)*dv(i)
         end if
       end do
       norm_XBLR = (4.0*pi*c*(RBLR*Rg)**2.0*UXBLR)/LXBLR
       do i=1,nvtot,1
        if (v(i).gt.vdpeak*b) then
         Lvexseed_j(i) = norm_XBLR*Lvexseed_j(i)
         Uvexseed_j(i) = Uvexseed_j(i)+
     &                       Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0)
         Uv(i,2) = Uv(i,2)+Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0)
         Uv(i,9) = Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0)
         if (Lvexseed_j(i).gt.0.0) then
c          write(3,*) log10(v(i)),
c     &           log10(v(i)*Lvexseed_j(i)/(4.0*pi*c*(RBLR*Rg)**2.0))
         end if
        end if
       end do
      end if
      !write(3,*) 'no no'

c     Do Torus Spectrum in jet frame
      f_IR = 0.3 !Default value
      !write(3,*) 'log(v) log(vUir)'
      if (zdiss.lt.RIR) then
       UIR = f_IR*Ldisc*BLF**2.0/(4.0*pi*c*(RIR*Rg)**2.0)
       mu1 = (1.0+(RIR/zdiss)**2.0)**(-1.0/2.0)
      end if
      if ((RIR.le.zdiss).and.(zdiss.le.3.0*RIR)) then
       write(*,*) '! zdiss > R_IR !'
       U1 = f_IR*Ldisc*BLF**2.0/(4.0*pi*c*(RIR*Rg)**2.0)
       mu1 = (1.0+(1.0/3.0)**2.0)**(-1.0/2.0)
       mu2 = (1.0-(1.0/3.0)**2.0)**(1.0/2.0)
       U2 = f_IR*Ldisc/(4.0*pi*c*(RIR*Rg)**2.0)*BLF**2.0/(3.0*beta)*
     &     (2.0*(1.0-beta*mu1)**3.0-(1.0-beta*mu2)**3.0-(1.0-beta)**3.0)
       grad = (log10(U2)-log10(U1))/(log10(3.0*RIR)-log10(RIR))
       const = log10(U1)-grad*log10(RIR)
       UIR = 10.0**(grad*log10(zdiss)+const)
      end if
      if (zdiss.gt.3.0*RIR) then
       write(*,*) '! zdiss > 3R_IR !'
       mu1 = (1.0+(RIR/zdiss)**2.0)**(-1.0/2.0)
       mu2 = (1.0-(RIR/zdiss)**2.0)**(1.0/2.0)
       UIR = f_IR*Ldisc/(4.0*pi*c*(RIR*Rg)**2.0)*BLF**2.0/(3.0*beta)*
     &     (2.0*(1.0-beta*mu1)**3.0-(1.0-beta*mu2)**3.0-(1.0-beta)**3.0)
      end if
      if (RIR.ne.0.0) then
       mu = cos(atan(RIR/zdiss))
       b = BLF*(1.0-beta*mu)
       v_IR = k/h*370.0*b !Default value
       LIR = 0.0
       do i=1,nvtot,1
         Lvexseed_j(i) = 0.0
         Lvexseed_j(i) = v(i)**3.0/(exp(v(i)/v_IR)-1.0)
         LIR = LIR+Lvexseed_j(i)*dv(i)
       end do
       norm_IR = (4.0*pi*c*(RIR*Rg)**2.0*UIR)/LIR
       do i=1,nvtot,1
         Lvexseed_j(i) = norm_IR*Lvexseed_j(i)
         Uvexseed_j(i) = Uvexseed_j(i)+
     &                       Lvexseed_j(i)/(4.0*pi*c*(RIR*Rg)**2.0)
         Uv(i,3) = Lvexseed_j(i)/(4.0*pi*c*(RIR*Rg)**2.0)
         if (Lvexseed_j(i).gt.0.0) then
c          write(3,*) log10(v(i)),
c     &           log10(v(i)*Lvexseed_j(i)/(4.0*pi*c*(RIR*Rg)**2.0))
         end if
       end do
      end if
      !write(3,*) 'no no'

      !write(3,*) 'log(v) log(vUtot)'
      do i=1,nvtot,1
        Uv(i,4) = Uvexseed_j(i)
        if (Uvexseed_j(i).gt.0.0) then
         !write(3,*) log10(v(i)),log10(v(i)*Uvexseed_j(i))
        end if
      end do
      !write(3,*) 'no no'

c      write(7,*) log10(M),log10(Udtot),log10(Ux),log10(UBLR+UXBLR),
c     &                                                       log10(UIR)

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DeltaSync(ng,g,n,dlogg,gmin_index,gmax_index,v,vb,Ub,
     &                           imin_sync,imax_sync,Rdiss,jsync,Uvsync)
c     Does synchrotron with delta function approximation. 

      implicit none
      integer i,imin_sync,imax_sync,nv,ng,gindex,gmin_index,gmax_index
      parameter (nv=3000)
      double precision pi,c,sigmat,Rdiss,dlogg
      double precision vb,v(nv),g(ng),n(ng),jsync(nv),Uvsync(nv),gi,Ub

      pi = 6.0*asin(0.5)
      c = 3.0e10
      sigmat = 6.65e-25

      do i=1,nv,1
         jsync(i) = 0.0
         Uvsync(i) = 0.0
      end do

c     Do synchrotron from vssa to vmax_sync
      do i=imin_sync,imax_sync,1
         gi = sqrt(3.0/4.0*(v(i)/vb))
         !Find gamma index
         gindex = int((log10(gi)-log10(g(1)))/dlogg)+2
         if ((gindex.ge.gmin_index).and.(gindex.le.gmax_index)) then
           jsync(i) = (sigmat*c*Ub)/(6.0*pi*vb)*g(gindex)*n(gindex)
           Uvsync(i) = (4.0*pi*Rdiss)/(3.0*c)*jsync(i)
         end if
      end do

c     Do self-absorbed spectrum below vssa with j prop to v^(5/2)
      if (imin_sync.gt.1) then
         do i=1,imin_sync-1,1
            jsync(i) = (v(i)/v(imin_sync))**(5.0/2.0)*jsync(imin_sync)
         end do
      end if

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Comp(ng,g,n,gmin_index,gmax_index,v,dv,dlogv,imin_sync,
     &                                   Uvsync,Uvex,jcomp,UvKN,KNindex)
c     Does comptonisation of synchrotron and external seed photons.
c     Calculates actual Usync used for Compton scattering leaving out 
c     photons cut off by Klein-Nishina limit. 

      implicit none
      integer i,j,nv,ng,imin_sync,vcindex,gmin_index,gmax_index,vc_old
      parameter (nv=3000)
      double precision pi,c,me,h,sigmat,Uvseed(nv)
      double precision g(ng),n(ng),v(nv),dv(nv),Uvsync(nv),Uvex(nv)
      double precision jcomp(nv),UvKN(nv),KNindex(nv),vc,dlogv

      pi = 6.0*asin(0.5)
      c = 3.0e10
      me = 9.1e-28
      h = 6.67e-27
      sigmat = 6.65e-25

      do i=1,nv,1
        jcomp(i) = 0.0
        UvKN(i) = 0.0
        KNindex(i) = 0.0
        Uvseed(i) = Uvsync(i)+Uvex(i)
      end do

      do j=imin_sync,nv,1
        vc_old = 0
        do i=gmin_index,gmax_index,1
          vc = 4.0/3.0*g(i)**2.0*v(j)
          !Find vc index 
          vcindex = int((log10(vc)-log10(v(1)))/dlogv)+2
          if ((vcindex.ge.1).and.(vcindex.le.nv)
     &                                 .and.(vcindex.ne.vc_old)) then
            if (vc.lt.g(i)*me*c**2.0/h) then !Klein-Nishina Limit
              jcomp(vcindex) = jcomp(vcindex)+
     &                  (sigmat*c*Uvseed(j)*dv(j))/(6.0*pi*v(j))*
     &                            g(i)*n(i)
              KNindex(j) = 1.0
            end if
            vc_old = vcindex
          end if
        end do
      end do

      do i=1,nv,1
         if (KNindex(i).gt.0.0) then
            UvKN(i) = Uvseed(i)
         end if
      end do

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine GetTemp(M,smdot,a,R,T)

      implicit none
      double precision M,R,T,y,y1,y2,y3,yms,a,pi,rms,isco,A1,B1,B2,B3,B
      double precision Mdot,G,Msun,sigma,Rg,c,eta,D,smdot

      pi = 6.0*asin(0.5)
      G = 6.67e-8
      sigma = 5.67e-5 !erg cm^-2 s^-1 K^-4
      c = 3.0e10
      Msun = 2.0e33
      Rg = G*M*Msun/c**2.0

      rms = isco(a)
      eta = 1.0-sqrt(1.0-2.0/(3.0*rms))
      Mdot = smdot*1.39e38*M/(eta*c**2.0)

      y = sqrt(R)
      y1 = 2.0*cos((acos(a)-pi)/3.0)
      y2 = 2.0*cos((acos(a)+pi)/3.0)
      y3 = -2.0*cos(acos(a)/3)
      yms = sqrt(rms)
      A1 = 1.0-(yms/y)-((3.0*a*log(y/yms))/(2.0*y))
      B1 =(3.0*(y1-a)**2*log((y-y1)/(yms-y1)))/(y*y1*(y1-y2)*(y1-y3))
      B2 =(3.0*(y2-a)**2*log((y-y2)/(yms-y2)))/(y*y2*(y2-y1)*(y2-y3))
      B3 =(3.0*(y3-a)**2*log((y-y3)/(yms-y3)))/(y*y3*(y3-y1)*(y3-y2))
      B = B1+B2+B3
      D = 1.0-(3.0/R)+((2.0*a)/R**(3.0/2.0))
      T =(3.0*G*M*Msun*Mdot*(A1-B))/(8.0*pi*sigma*R**3*Rg**3*D)
      T = T**(1.0/4.0)

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function isco(a)

      implicit none
      double precision a,z1,z2,rms,isco

      z1 = 1.0+(1.0-a**2)**(1.0/3.0)*((1.0+a)**(1.0/3.0)+(1.0-a)**
     &                                                       (1.0/3.0))
      z2 = sqrt(3*a**2+z1**2)
      rms = 3.0+z2-sqrt((3.0-z1)*(3.0+z1+2.0*z2))
      !rms = 3.0+z2+sqrt((3.0-z1)*(3.0+z1+2.0*z2))
      
      isco = rms

      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
