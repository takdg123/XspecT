c Subroutine tree starting at skclt which is called from tcaf
c
c  skclt     den
c            dsk2sok       dsksok
c            initial       den
c            zgamln
c            integre       den
c                          rk4         derivs       den
c            alfabeta      d0          zgamln
c            sok2dsk       sokdsk
c                          thialfhi
c                          thialflo
c                          tloalfhi
c                          tloalflo
c                          storall
c            ctcompspec    compt       zgammi       zgamln
c                                                   zgamln
c                                      zgamln
c                                      zyyit2       zgamln
c                          d00         zgamln
c            bbspec


      SUBROUTINE tcaf(ear, ne, param,ifil,photar,photer)
      implicit double precision(a-h,o-z)


      INTEGER ne, ifil

      REAL param(5), ear(0:ne), photar(ne),photer(ne)
      REAL start(ne), end(ne), fstart(ne), fend(ne)
        
C Subroutine to do the actual calculation of the tcaf model
C Parameters :
C       ear      r           i: energies for output
C       ne       i           i: number of energies
C       param    r           i: model parameters
C       start    r           w: workspace array used for rebinning
C       end      r           w: workspace array used for rebinning
C       fstart   r           w: workspace array used for rebinning
C       fend     r           w: workspace array used for rebinning
C       photar   r           r: final model spectrum
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

cc        SUBROUTINE ct95flux(emdotdsk,emdot0,ebh,xs,cr, spec,beng,eng) 
cc   ! ebh in solar mass unit,              ADDED by DD on 14/01/2012
cc   ! Here 4 variable parameters (emdotdsk,emdot0,ebh,xs) from from program
cc   ! and other 3 parameters (spectra,beng,eng) are output.
cc   ! spec are in photons/cm^2/s & beng are corresponding energy in keV,
cc   ! eng are the actual energy in keV, from where beng are calculated.

      include 'ct.inc'

c        implicit double precision(a-h,o-z)
      real spec(NPTSm1),beng(NPTSm1)
      dimension eng(NPTS),totflux(NPTS)
c      common /presok/ hdd,prex(1000),pret(1000),ipresok  
c      common/resu/ prex1(200),tbb1(200),ts1(200),coskhi1(200),
c     &     replum1(200),frac1(200),frfr1(200)


c          character*100 abcd   !Line Added by DD
!        do 10 i=1,1
c        read (11,*) emdotdsk,emdot0,embh,xout,tout,dr0,xs,xin
c          close(11)       !Line Added by DD
c        read (11,92) emdotdsk,emdot0,embh,abcd  !Line Added by DD
c92             format(f4.2,1x,f4.2,1x,f10.9,a100) !Line Added by DD

cc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Print*,'Infile: fort.12[5para:emdotdsk,emdot0,ebh(in M_Sun),xs,R]'
!        read (12,*) emdotdsk,emdot0,ebh,xs,cr
!        print*, emdotdsk,emdot0,ebh,xs,cr
cc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      emdotdsk=param(1)
      emdot0=param(2)
      ebh=param(3)
      xs=param(4)
      cr=param(5)

      embh=ebh*1.d-8
      xout=500.0
!!        tout=1.6e+9
      emu=11./24.
      gma=5./3.
      tout=emu*EMP*CLIGHT**2/(2.*gma*BOLTZ*xout) !!! CHANGED in v.0
      dr0=-0.001
      xin=3.0
      write(13,93) emdotdsk,emdot0,embh,xout,tout,dr0,xs,xin,cr 
 93   format(2f8.3,e12.3,f8.2,e10.2,f7.3,3f8.2)

!       Print*,'emdotdsk, emdot0, embh, xout, tout, dr0, xs, xin, R:'
!       Print*, emdotdsk,emdot0,embh,xout,tout,dr0,xs,xin,cr
c                pause
      call skclt(i,emdotdsk,emdot0,embh,xout,tout,dr0,xs,xin,cr,
     &     totflux,eng)

      distconst = 1.0d0/9.523d44
      do l = 1, NPTS
         totflux(l) = totflux(l)*distconst
      enddo

cc ! Calculate integrated flux (spectral flux) in eng(l) to eng(l+1)
cc ! and write it to spec(i). Unit of spec is photon/cm^2/s
cc ! In 'fort.23' avg. bin energy spectra

      amaxv=1.d-30
       
      do 40 l=2,NPTS
         x = (eng(l)+eng(l-1))/2.0d0
         y = totflux(l)
         ym1 = totflux(l-1)

         dx = abs(eng(l)-eng(l-1))
         dy = abs(y-ym1)
         if (y.lt.ym1) then
            spec(l-1)=SNGL(dx*y + 0.5*dx*dy)
            beng(l-1)=SNGL(x)
            write(23,*)beng(l-1),spec(l-1)
         else
            spec(l-1)=SNGL(dx*ym1 + 0.5*dx*dy)
            beng(l-1)=SNGL(x)
            write(23,*)beng(l-1),spec(l-1)
         endif
c            write(*,*)beng(l-1),spec(l-1)
         if(spec(l-1).gt.amaxv) amaxv=spec(l-1)
 40   continue
!c                 write(23,*)
!          print*, 'Maximum value of spectral flux: ', amaxv
!           
!c        do 11 k=ipresok,2,-1
!c        write(15,51)prex1(k),tbb1(k),ts1(k),coskhi1(k), !Commented by DD
!c     &  replum1(k),frac1(k),frfr1(k)                    !Commented by DD
!c51      format(7e11.3)                                  !Commented by DD
!c11      continue
!10      continue
!
!        close(all)

c  Rebin onto passed energies

      CALL inibin(549, beng, ne, ear, start, end,
     &            fstart, fend, 0)
      CALL erebin(549, spec, ne, start, end, fstart, fend, photar)

      RETURN                    ! ADDED FOR XSPEC
      end

      subroutine skclt(it,emdodi,emdodi0,embhh,xxout,ttout,
     &     ddr0,xxs,xxin,ccr, totflux,eng)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      dimension totflux(NPTS),eng(NPTS)
      dimension comflux(NPTS),bbflux(NPTS)
      dimension xx(7000),tt(7000),tel(7000)
      common /dat/ yy(2),dydt(2)
      common /bbneed/ bb(NPTS),xxdd(MXANUL),tbb(MXANUL),ts(3000)
      common /comneed/com(NPTS),eemm,tin,tell,tao,teff,xshock,schi,xinn
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /presok/ hdd,prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /alfabet/ x0,alfax0,rdel,temp,t5
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr 
      character(20) csumtel
c     print *, 'm_dotdsk, m_dothalo, mass of black hole(10^8)'

      emdotdsk=emdodi
      emdot0=emdodi0
      embh=embhh
      xout=xxout
      tout=ttout
      dr0=ddr0
      xs=xxs
      xin=xxin
      cr=ccr
c     print *, emdotdsk, emdot0,embh
c     print *, xout,tout,dr0,xs,xin
      fracvec=1.0
      eemm=embh
      enfac=20.0
c  !     This is the emdot of the halo
      emdot=0.2*2.e+33/3.15e+7*emdot0*embh
      emdotall=0.2*2.e+33/3.15e+7*(cr*emdot0+emdotdsk)*embh !!! CHANGED in v.0
      embhh=embh*1.e+8
c  ! temp 1.56e+11 for xs=10.0 when vertical equilibrium is used
c  ! temp 1.e+11 for xs=10.0 when thin disk assumption is used

!************************************************
cc        tshock=1.56e+11*(10.0/xs)
      tshock_skc=1.56e+11*(10.0/xs) 
!        emu=0.5   
cc new emu by considering 75% H (1e, 1p) & 25% He (2e, 1nucleus) =>(3/4*1/2 + 1/4*1/3)
      emu=11./24.               !!! CHANGED in v.0  
!        print*,emp,c,boltz,xs
      tshock=emu*EMP*(cr-1.)*CLIGHT**2/(cr**2*BOLTZ*(xs-1.)) !! tshock interms of compression ratio

!c        print *,'tshocks (new & sir):',tshock,tshock_skc
!************************************************

      iter=0
      rho0=emdot/FORPI/CLIGHT
      rhoall=emdotall/FORPI/CLIGHT
      sigma=7.56e-15*CLIGHT/4.0 !s-b const
      tercep=0.25
      rs=3.e+5*embh*1.e+8
      schi=rs
      surface=PI*((xout-xs)*rs)**2 !!! CHANGED in v.0
c    !   shock height is obtained from a*sqrt(x)*(x-1.0)
c     
c        xs=9.9999

!************************************************
cc        sokheit=1./4.*sqrt(5./xs)*xs**1.5       !Shapiro p.441
      sokheit_skc=1./4.*sqrt(5./xs)*xs**1.5
      pol=5./3.
      sokheit=sqrt(pol*(cr-1.)*xs*(xs)/cr**2) !! sokheit in terms of compression ratio

!c        print *,'sokheits (new & sir):',sokheit,sokheit_skc
!************************************************
        
c    !    diskheit is assumed=1.0 here.
      hdd=0.0
      sokarea=TWOPI*xs*(sokheit-hdd)*rs**2
      opacity=0.4
      ogamda=20.0
      den0=rho0/rs**2
      denall=rhoall/rs**2
      gamhi=5./3.
      gamlo=4./3.
      te=EME*CLIGHT**2/BOLTZ
      ej=1.44e-27               !j
      capc=9.3e+12
      conp=4.*EMP*rs**3*PI/BOLTZ/emdot
      conep=1.5*BOLTZ*ogamda/(capc*EMP*EME)*EME**1.5
      conpro=BOLTZ/EME/CLIGHT**2
c   !    This is the emdot of the disk
      emdo=0.2*2.e+33/3.15e+7*emdotdsk*embh
      emdot17=emdo/1.e+17
c  !  First we compute the flux intercepted by the shock from the disk.
c  ! factor of 8 is because r in Shapiro book (p.441) is in unit of GM/c^2
c          print *, 'emdot17',emdot17
      sscons=fracvec*5.e+26/embhh**2*emdot17/8.0
c    One dummy run is done here to compute tau00
      tau0=0.0
      do 82 k=1,100000
         rr=xs-0.01*(k-1)
         if(rr.le.xin)goto67
         tau0=tau0+0.4*den(rr)*0.01*rs
 82   continue
 67   continue
      alpha=0.5
c     tau00=0.92*emdot0
      tau00=tau0
c     print *, 'tau00=',tau00
      funtau=1.5*exp(-sqrt(tau00+2.))
      c1=1.0

c     loop point 
 500  continue
      call dsk2sok (bbtens,cept,sumtem,sumthe,cept2)
      ind=1
c   !    effective temperature of the disk
      diskarea=PI*(xout**2-xs**2)*rs**2
cc        ssflux=sigma*tss**4*emdotdsk
cc        ucons=ssflux*surface*tercep/sokarea*sokheit/(0.5*xout)
c          print *, 'bbtens',bbtens
c           write(10,*) iter,bbtens
      ssflux=bbtens/diskarea
c   ! This would have been ucons if monotemperature disk was present
      ucons=cept/sokarea
c   ! stupid sum
      tssbb=(ssflux/sigma)**0.25
c   ! better sum
      tsseff=(cept2/(sumtem*sigma))**0.25
c         print *, 'effective disk temp', tsseff,tssbb
      tssold=tss
      tss=tsseff
c   !  effective inclination theta of the shock to the disk
c   !  stupid sum
      tantheta=sokheit/(xout-xs)
      theta=atan(tantheta)
c         print *, 'bad theta', theta*180.0/pi
c  !  to compute the average flux emitted, we need to
c ! multiply by local cos(xi) on the sok surface.
      clina=1./cos(theta)
c ! better sum
      tantheta=sumthe/sumtem
      theta=atan(tantheta)
c         print *, 'better theta', theta*180.0/pi
      clina=1./cos(theta)
c
      conheit=sqrt(1.666*BOLTZ/EMP)
      call initial
c ! This factor is to take care of the recoil effect. For
c ! alpha>1, facalpha=1
      if(alpha.lt.1) then
         facalpha=1.-(1.+alpha/4.)*(1.-alpha)
         gamterm=exp(log(alpha)+log(alpha+3.0)+zgamln(alpha+4.0)
     &        +zgamln(alpha)+zgamln(1.-alpha)-zgamln(2.*alpha+4.))
      else
         facalpha=1.0
      endif
      d0alf=(3.*((alpha+3.0)*alpha+4.0)*exp(zgamln(2.*alpha+2.)))/
     &     ((alpha+3.)*(alpha+2.)**2)
c        print *,'gamterm', gamterm,d0alf
cc      iii=20000
      iii=20
      dr=dr0
      n=iii*50000
      in=1
      inp=1
      tauhori=0.0

      do 10 i=1,n
c        if(iter.eq.0) then
         if(r.ge.xs)then
            if(mod(i,1000).eq.0) then
               prex(inp)=r
               pret(inp)=tauhori
               inp=inp+1
cx        write(*,47)r,yy(1),yy(2),tauhori
c      endif
            endif
         endif
         if(r.lt.xin)goto99
         if(r.lt.xs) then 
            if(ind.eq.1) then
               yy(1)=tshock
               ind=ind+1
            endif
c         dr=-0.0000001
            dr=-0.0001
            if(mod(i,1000).eq.0) then
               xx(in)=r
               tt(in)=tauhori
               tel(in)=yy(2)
               in=in+1
cx        write(*,47)r,yy(1),yy(2),tauhori
               if(yy(1).lt.1.0e+7)goto99
            endif
         endif
         call integre(i)
 10   continue
 99   continue  

c 47   format(e14.5,3e15.5)
c   !  calculation of average T_e assuming a sphere
      num=in-1
      ipresok=inp-1
c     print *,num,ipresok
      tau00=tt(num)
      funtau=1.5*exp(-sqrt(tau00+2.))
      sumet=0.0
      sumtel=0.0
      sumcon=0.0
      do 30 j=1,num
         ttt=tt(j)
c  ! geometric factor
         geo=(tau00-ttt)**2     !2clcult avg elc temp
c ! the source factor.
         et=((1.-funtau)*cos(PIBI2*(1.-ttt/tau00))+funtau) !gtau
         et2=geo*et**2
         sumet=sumet+et2                        !deno f avg elc tmp
         sumcon=sumcon+geo*et*exp(-ttt*clina)
         sumtel=sumtel+et2*tel(j)               !numtorf avg elc tmp
 30   continue
      c1=sumcon/sumet
c       mean electron temperature
      telmean=sumtel/sumet      !avg (CT95 p.631)
        
!        print *,'sumtel,sumet',sumtel,sumet
ccccccccccccccccccccccccccccccc  NAN cccccccccc
cc        if (sumtel.lt.1.e-30) stop
      open(40)
      write(40,555)sumtel
cc        write(*,555)sumtel
 555  format   (d15.7)
      close(40)
      read(40,556)csumtel
 556  format(a15)
      close(40)
      if(csumtel.eq.'            NAN'.or.csumtel.eq.'           -NAN'
     & .or.csumtel.eq.'            NaN'.or.csumtel.eq.'           -NaN'
     & .or.csumtel.eq.' NaN'.or.csumtel.eq.'-NaN'.or.csumtel.eq.' NAN'
     & .or.csumtel.eq.'-NAN'.or.sumtel.lt.1.e-30)then
         write(*,*)'NAN / zero obtained'
         a=0.d0
         do k111=1,549
            write(23,*)a,a
            write(21,*)a,a
            write(22,*)a,a
         enddo
         write(15,52)emdotdsk,emdot0,a,a,a,a,a
 52      format(2f8.3, f8.5,4e14.5)
         write(21,*)a,a
         write(22,*)a,a
!                   close(all)
         close(21)
         close(22)
         close(23)
         stop
      endif
cccccccccccccccccccccccccccccccccccccccccccccc
c         print *,'c1', c1
c       in keV
      telmean0=telnew
      telmean=telmean*BOLTZ/KEVTOERG
      telnew=telmean
      if(iter.gt.2)then
         telmean=0.5*(telmean+telmean0)
         tss=0.5*(tss+tssold)
      endif
      tinject=BOLTZ*tss/KEVTOERG
c        print *, 'telmean,tinject', telmean,tinject
      xnotnew=3.*tinject/telmean
      alphaol=alpha
      call alfabeta
c        if(alfax0.gt.0.95) alfax0=0.95
c        if(alfax0.lt.0.35) alfax0=0.35
      if(iter.gt.2)then
         alpha=0.5*(alfax0+alphaol)
      else
         alpha=alfax0
      endif
      if(iter.gt.3) then
         if(alpha.gt.0.95) then
            if(alpha.lt.1.05) then
               return
            endif
         endif
      endif
c        alpha=(alfax0+alpha0)*0.5
c            endif
c ! The albido calculation from Kallman+Titarchuk stuff
c ! x_star=10/500=0.02 (kallman cutoff at 10 Kev)
c ! x_hat=2.0(in units of kt_e/m_ec^2 , Marachi papers
c ! the fraction absorbed formula 1-A= xstar**(1-alpha)...
c ! has come from the recoil effect and the photoionization effects
c     absorb=1-A
      if (alpha.lt.1.0) then
         alp1=1.-alpha
         alp3=1.5-alpha
         absorbhi=(sqrt(telmean/500)*
     &        (2.**alp3-0.02**alp3)*alp1/alp3)/2.**alp1
         absorblo=(0.02)**alp1/2.**alp1 + absorbhi
         facalpha=1.-(1.+alpha/4.)*(1.-alpha)
         gamterm=exp(log(alpha)+log(alpha+3.0)+zgamln(alpha+4.0)
     &        +zgamln(alpha)+zgamln(1.-alpha)-zgamln(2.*alpha+4.))
      else
         absorbhi=0.5
         absorblo=0.5
         facalpha=1.0
      endif
      d0alf=(3.*((alpha+3.0)*alpha+4.0)*exp(zgamln(2.*alpha+2.)))/
     &     ((alpha+3.)*(alpha+2.)**2)
c         print *,'d0alf', d0alf
c       enhancement factor
      alm1=alpha-1.
      enfacold=enfac
      if(alpha.lt.1.0)then
         enfac=gamterm*(xnotnew**alm1 -1.) !E_comp
      else
         enfac=(alpha*(alpha+3.)*(1.-(alpha+4)/(2.*alpha+3.)*
     &        xnotnew**alm1))/((alpha+4.0)*alm1)
      endif
      enhi=enfacold*10.0
      enlo=enfacold/10.0
      if(enfac.gt.enhi)enfac=enhi
      if(enfac.lt.enlo)enfac=enlo
      emitflux=enfac*ucons
c        print *, 'fac,cons,emitflux', enfac,ucons,emitflux
c        print *, 'absorb_hi, absorb_lo',absorbhi,absorblo
c          print *,'results',alpha,telmean,tinject,tss
      call sok2dsk(emitflux)
      if(iter.gt.3) then
         ap=1.05*alphaol
         am=0.95*alphaol
c       if((alpha.gt.am).and.(alpha.lt.ap)) then  !Commented by DD
         write(15,51)emdotdsk,emdot0,alpha,telmean,tinject,tss,tau00
 51      format(2f8.3, f8.5,4e14.5)
!      print *,'final parameters:alpha=',alpha
c      print *,'final parameters telmean=',telmean
c      print *,'final parameters injection temp=',tinject,tss
c        print *,'final tau00=',tau00
c        print *, 'alfa final',alpha,'telmean',telmean
         goto 300
c        endif          !Commented by DD
      endif            
      iter=iter+1
      if (iter.lt.30) goto 500
c     ! After we are satisfied with alpha calculation and the
c     ! reprocessed flux, the outgoing Compton spectrum has to be calculated.
 300  continue
      tin=tinject
      tell=telmean
      tao=tau00
      teff=tss
      xshock=xs
      schi=rs
      xinn=xin
      call ctcompspec
      call bbspec(embh,idisk)
      flo=13.5
      fhi=21.0

      dfreq=(fhi-flo)/DFLOAT(NPTS)

      hkev = HPLANCK/KEVTOERG
      do 20 ii=1,NPTS
         fo=flo+dfreq*(ii-1)  
         f=10.**fo              ! in Hz
         hf=HPLANCK*f

         eng(ii)=hkev*f     ! bin in keV
         totflux(ii)=((com(ii)+bb(ii))/f/(hf))/hkev
         comflux(ii)=(com(ii)/f/(hf))/hkev 
         bbflux(ii)=(bb(ii)/f/(hf))/hkev

c      !F_nu = (com(ii)+bb(ii))/f           in erg/cm^2/s/Hz
c      !     = ... /(h*f)                   in photons/cm^2/s/Hz
c      !     = ... *(1.6e-9/6.625e-27)      in photons/cm^2/s/keV
c      ! since /Hz can be converted to /keV by multiplying Hz/keV factor
c      ! now 1/h in egs unit Hz/erg=Hz/1.6e-9keV, implies 
c      ! Hz/keV = 1.6e-9/6.625e-27, since in cgs h = 6.625e-27 erg.s
c
cc     ! for BH distance (D) luminocity should decrease by 1/(4*pi*D^2)
cc     ! for e.g. 3kpc distant source, totflux should be divided by a 
cc     ! factor ~ 1.e+45, we should include this as a Normalization
cc     ! totflux(ii)=((com(ii)+bb(ii))/f/(h*f))*(1.6e-9/h)/1.e+45

!       write(22,*)eng(ii),totflux(ii)
!       write(21,*)f,(com(ii)+bb(ii))

         write(22,72)eng(ii),totflux(ii),comflux(ii),bbflux(ii)
         write(21,71)f,(com(ii)+bb(ii)),com(ii),bb(ii)
c       write(21,*)fo,log10(com(ii)+bb(ii)),log10(com(ii)),log10(bb(ii))

 20   continue
 71   format(5e14.6)
 72   format(f14.8,4e14.6)

      return 
      end

        
c
c     !  initial data
      subroutine initial
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      common /dat/ yy(2),dydt(2)
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
      yy(1)=tout
      yy(2)=tout/100.0
      r=xout
      gamap=0.0
c    !    energy exchange between electron and proton due to 
c    !  Coulomb scattering. 
c    !   Here r is in units of Schwarzschild radius.
      enr=den(r)/EMP
c         gamaep=1.5*enr**2*appa*(yy(1)-yy(2))*ogamda   !yy(1)=T_p,2=T_e,enr=n
c     &/(c*emp*eme)*(eme/yy(2))**1.5
      gamaep=conep*enr**2*(yy(1)-yy(2))/yy(2)**1.5
      amdaib=ej*enr**2*(EME/EMP*yy(1))**0.5
      amdap=gamaep+amdaib
      amdab=ej*enr**2*sqrt(yy(2))*(1.+4.4e-10*yy(2))
c         ignore now comptonization and cyclo-synchrotron emissions           
      dydt(1)=-yy(1)/r-conp*2./3.*r*r*(-amdap)
      if(yy(2).lt.te) then
         gamael=gamhi
      else
         gamael=gamlo
      endif
      dydt(2)=-yy(2)/r*1.5*(gamael-1.)-conp*(gamael-1.)*r**2*
     &     (gamaep-amdab)
c      print *,'initial done', dydt(1),dydt(2)
      return
      end

      subroutine integre(i)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      common /dat/ yy(2),dydt(2)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
      if(r.lt.xs) then 
         tauold=tauhori
         tauhori=tauhori+0.4*den(r)*abs(dr)*rs !den(r)
      endif
      call rk4(yy,dydt,2,r,dr)
      r=r+dr
      return
      end

      subroutine derivs(tt,ym,dymdt)
      implicit double precision (a-h,o-z)
      include 'ct.inc'
      dimension ym(2),dymdt(2)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
c      cooling and heating terms are from Marachi, Colpi and Treves
c  proton heating is ignored.
      gamap=0.0
c        energy exchange between electron and proton due to 
c      Coulomb scattering. 
c       Here r is in units of Schwarzschild radius.
      enr=den(r)/EMP
      if(emdot0.gt.3.0) then
         if(r.lt.xs) then
            gamaep=conep*enr**2*(ym(1)-ym(2))/ym(2)**1.5
            amdaib=ej*enr**2*(EME/EMP*ym(1))**0.5
            amdap=gamaep+amdaib
            amdab=ej*enr**2*sqrt(ym(2))*(1.+4.4e-10*ym(2))
         else
            gamaep=0.0
            amdaib=0.0
            amdap=0.0
            amdab=0.0
         endif
      else
         gamaep=conep*enr**2*(ym(1)-ym(2))/ym(2)**1.5
         amdaib=ej*enr**2*(EME/EMP*ym(1))**0.5
         amdap=gamaep+amdaib
         amdab=ej*enr**2*sqrt(ym(2))*(1.+4.4e-10*ym(2))
      endif               
c   ignore the effect of cyclo-synchrotron emission from the disk
c  use the shakura-sunyaeve photons for the comptonization
      if(r.le.xs)then
c     if(ym(2).lt.1.e+7) then
c        amdacomp=0.0
c        else
c     Here we use more corrected formula from Titarchuk and Lubierski
         coym=conpro*ym(2)
c       tlterm=(alpha+3.)*coym/(1.+coym)+4.*d0alf**(1./alpha)*coym**2
         tlterm=1.0
         amdacons=SIGMAT*enr*facalpha*tlterm
c       amdacons=4.*sigmat*enr*conpro*ym(2)*facalpha
c     &*(1.+4.*conpro*ym(2))
c        do 95 j=1,idisk
c   for a better calculation one needs to use separate annulus to be
c seperate sourcess
         xnot=3.0*tss/ym(2)
         alm1=alpha-1.0
         if(alpha.lt.1) then
            elbyel0=(xnot**alm1 -1.0)*gamterm !E_comp
         else
            elbyel0=((alpha+3.)*alpha*(1.-xnot**alm1))/(alpha+4.)/alm1
         endif
c     taueff=tauhori*clina
         et= ((1.-funtau)*cos(PIBI2*(1.-tauhori/tau00))+funtau) !g(tau)
         ubar=ucons*c1*et       !ucons?
         amdacomp=amdacons*ubar*elbyel0
c        amdacomp=amdacons*ubar
c         write(6,62) r,gamaep,amdacomp,elbyel0
c 62      format(4e13.5)
c           endif
c 95      continue
      else
c        heit=conheit*sqrt(ym(1)*r**3)*rs/c
c        tauvert=0.4*den(r)*heit
c        dsk=sscons/r**3*sqrt(1.-3.0/r)
c        ubar=dsk*exp(-tauvert)
c        amdacomp=4.*sigmat*ubar*enr*conpro*ym(2)*facalpha*elbyel0
c     &*(1.+4.*conpro*ym(2))
         amdacomp=0.0
      endif
      amdae=amdab+amdacomp
c        write(76,67) heit/rs, tauvert,tauhori,
c     &  amdacomp,amdab,gamaep,amdaib
c 67   format (7e11.3)
      acomp=amdacomp
      dymdt(1)=-ym(1)/r-conp*2./3.*r*r*(-amdap)
      if(ym(2).lt.te) then
         gamael=gamhi
      else
         gamael=gamlo
      endif
      gamnet=gamaep-amdae
      dymdt(2)=-ym(2)/r*1.5*(gamael-1.)-conp*(gamael-1.)*r**2*gamnet
c     & (gamaep-amdae)
c       write(*,22)r,amdacomp,enr
c         if (r.lt.xs) then
c        print *, 'ubar,amdacomp,amdab',ubar,amdacomp,amdab,
c     & gamaep, amdaib,amdap
c        write(*,22) r,ym(1),ym(2),dymdt(1),dymdt(2)
c          endif
c        write(*,22)r,tauhori,acomp,enr
c 22   format(5e11.3)
      return
      end
c
      SUBROUTINE RK4(Y,DYDX,N,X,H1)
      IMPLICIT NONE
      DOUBLE PRECISION Y(N),DYDX(N)
      DOUBLE PRECISION X, H1
      INTEGER N

      INTEGER NMAX
      PARAMETER (NMAX=8)
      DOUBLE PRECISION YT(NMAX),DYT(NMAX),DYM(NMAX)
      DOUBLE PRECISION HH, H6, XH
      INTEGER I

      HH=H1*0.5
      H6=H1/6.
      XH=X+HH
      DO 11 I=1,N
         YT(I)=Y(I)+HH*DYDX(I)
 11   CONTINUE
      CALL DERIVS(XH,YT,DYT)
      DO 12 I=1,N
         YT(I)=Y(I)+HH*DYT(I)
 12   CONTINUE
      CALL DERIVS(XH,YT,DYM)
      DO 13 I=1,N
         YT(I)=Y(I)+H1*DYM(I)
         DYM(I)=DYT(I)+DYM(I)
 13   CONTINUE
      CALL DERIVS(X+H1,YT,DYT)
      DO 14 I=1,N
         Y(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
 14   CONTINUE
      CALL DERIVS(X+H1,Y,DYDX)
      RETURN
      END

      function  den(x)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
c     density in cgs unit
      if(x.gt.xs)then
         den=den0/x**1.5
      else
         den=denall/x**1.5
      endif
      return
      end

c      double precision function  velo(x)
c      implicit none
c      double precision x
c      include 'ct.inc'
c      velocity in cgs units
c      velo=CLIGHT/sqrt(x)
c      return
c      end
      
      subroutine  alfabeta
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      PARAMETER (MXROW=2000, MXCOL=20, MXVEC=20, MXCMD=30)
      parameter (MLAG = 15)
      dimension param(5)
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /alfabet/ x0,alfax0,rdel,temp,t5
      common /paraIn/ cr
      param(1)=0.0
      param(2)=tinject
      param(3)=telmean
      param(4)=tau00
      param(5)=2.0
c=====>
C 1. Z    ** REDSHIFT
C 2. T0   ** WIEN TEMPERATURE (keV)
C 3. TEMP ** PLASMA TEMPERATURE (KEV)
C 4. TAUP ** PLASMA OPTICAL DEPTH
C 5. APPRX * <= 1.0 DISK
C           >  1.0 SPHERE
C PHOT ** PHOTON SPECTRUM IN PHOTONS/CM/CM/S/KEV
      ro=1.
      bol2i=0.0
      ZFAC = 1.+PARAM(1)
      T0 = PARAM(2)
      TEMP = PARAM(3)
      TAUP = PARAM(4)
      apprx=param(5)
      t5 = temp/511.
c dimensionless soft photon energy
      x0 = 3.*t0/temp
      t23=taup+(2./3.)
      taue=taup/(1.0+(temp/39.2)**0.86) !roseland mean(pac 83,HT 95)
c        print *,'taue',taue
      if (apprx.gt.1.0) then 
c BETA FOR SPHERE
         TB1=PI*PI*(1.-exp(-0.7*taup))/3./t23/t23
         TB2= dexp(-1.4*taup)*dlog(4./(3.*taup))
         beta = tb1+tb2
c              write(*,*) 'beta = ',beta
         if (taue.lt.0.01) then 
            rdel = taue/2.      !reduced optcl dpth/prbalty f phtn 2 b sctrd
         else
            rdel=1.0-3./taue*(1.-2./taue+2./taue/taue*(1-exp(-taue))) !HT95
         endif
      elseif (apprx.le.1.0) then
c BETA FOR DISK      
         tb1 = PI*PI*(1.-exp(-1.35*taup))/12./t23/t23 !CT95
         tb2 = 0.45*dexp(-3.7*taup)*dlog(10./(3.*taup))
         beta = tb1+tb2
c     write(*,*) 'beta = ',beta
         if (taue.lt.0.01) then
            rdel = taue/2.
         else
            rdel = 1.0-(1.0-exp(-taue))/taue !HT95
         endif
      endif
      f0theta=2.5*t5 + 1.875*t5*t5*(1.0-t5) !t5=men elc tmp/511
      gam0=beta/t5
c        write(*,*) 'gam0 = ',gam0
c compute the power-law index
      tal0 = dsqrt(2.25+gam0)-1.5
 500  bola = 1. + (tal0+3.)*t5/(1.+t5)+4.*d0(tal0)**(1./tal0)*t5*t5
      tal = beta/dlog(bola)
c      write(*,*) 'tal0 tal ',tal0,tal
      if (dabs(tal-tal0).le.1.0d-5) then
         alfax0 = tal
      else
         tal0=tal
         go to 500
      endif
c      print *, 'alfax0 = ',alfax0
      return
      end

      double precision function d0(x)
      implicit none
      double precision x
      double precision x1, x2, x3, bol, argg, db
      double precision zgamln
      x1=x+1.d0
      x2=x+2.d0
      x3=x+3.d0
      bol = x3*x + 4.d0
      argg = 2.d0*x1
      db = zgamln(argg)
      if (db.lt.-70.) then 
         d0=0.0
      elseif (db.gt.70.) then 
         d0=1.e34
      else
         d0 = 3.d0*bol*dexp(db)/x3/x2/x2
      endif
      return
      end

      DOUBLE PRECISION FUNCTION ZGAMLN(az)
      implicit none
      double precision az
      include 'ct.inc'
      double precision a(7)
      double precision z, s
      integer i
      A(1) = 1./12.
      A(2) = 1./30.
      A(3) = 53./210.
      A(4) = 195./371.
      A(5) = 22999./22737.
      A(6) = 29944523./19733142.
      A(7) = 109535241009./48264275462.
      z=0.d0
      z=az
c      write(*,*) 'zgamln : az  z = ',az,z
      s = z

      DO 1, I=1,6
         s = z + A(8-I)/s
 1    CONTINUE

      s = A(1) / s
      s = s - z + (z - 0.5)*DLOG(z) + 0.5*DLOG(TWOPI)
      ZGAMLN = s
      RETURN
      end

      subroutine sok2dsk(emitflux)         
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      common /bbneed/ bb(NPTS),xxdd(MXANUL),tbb(MXANUL),ts(3000)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /presok/hdd, prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
      stephan=7.56e-15
      appaes=0.4
      sigm=5.67e-5
c     height of the sok
      hs= sokheit*rs
c     height of the disk
      hd=hdd*rs
      sumn=0.0
      tensn=0.0
      sokarea=TWOPI*xs*(sokheit-hdd)*rs**2
c           print *,'sokarea',sokarea
c         print *, 'sokheit', sokheit
      cept=0.0
      xxdd(1) = xs
      do 31 k=ipresok,2,-1
         kk=ipresok-k+1
c       call the subroutine to sum over all the contributions
c from the shock, at the reprocessing site.
c           print *, 'emitflux before',emitflux
         call sokdsk(k,emitflux,rep,replum,coskhi)
c        cosss=(prex(k)+1.0-xs)/sqrt((prex(k)+1.0-xs)**2+(sokheit-hdd)**2)
         sinss=sqrt(1.-coskhi**2)
c the temperature at the reprocessing site, assuming the
c black body approximation only
c         tsnew=(rep/sigm)**0.25

!!****************************************************************
         tbb(k)= ((dskflux0(ipresok-k+1))/sigm)**0.25
!        if(emdotdsk.ge.1.0)then
!                 tbb(k)= (((dskflux0(ipresok-k+1))/sigm)**0.25)*1.8 !! spectral hardening factor (Shimura & Takahara 1995)
!         else
!                 tbb(k)= (((dskflux0(ipresok-k+1))/sigm)**0.25)
!         endif        
!!***************************************************************

c       In the next iteration, dskflux(ipresok-k+1) should
c be same as rep+dskflux(ipresok-k+1)
c           write (7,*)prex(k),dskflux0(kk)               
         dumtem=((dskflux0(kk)+rep)/sigm)**0.25
         if(dumtem.gt.5.e+5) then
            if(alpha.gt.1.0) then
               call thialfhi(sinss,absorb)
            else
               call thialflo(sinss,absorb)
c     absorb=0.1
            endif
         else
            if(alpha.ge.1.0)then
               call tloalfhi(sinss,absorb)
            else
               call tloalflo(sinss,absorb)
            endif
         endif
         dskflux(kk)=coskhi*fracvec*rep*absorb+dskflux0(kk)
c           dskflux(kk)=dskflux(kk)
c           write (7,*)prex(k),dskflux(kk)
c       ts(k) is the effective temperature of the annulus
c      when the reflected and SSdisk fluxs are considered
c      together. This  goes in the common statement
c      for iteration.
         ts(k)=((dskflux(kk))/sigm)**0.25
         xxdd(k)=prex(k)
c          print *,'replum inside',replum
c           print *,'sokarea',sokarea
         frfr=replum/(emitflux*sokarea)
         cept=cept+replum
c           print *, 'cept', cept,'emitflux after ',emitflux
         frac=cept/(emitflux*sokarea)
c      if(kk.lt.3)then
c      write(11,18)kk,dskflux0(kk),dskflux(kk),xd/rs,rep,replum,coskhi
c18    format(i4,5e13.4,f10.4)
c      endif
         call storall(k,prex(k),tbb(k),ts(k),coskhi,replum,frac,
     & frfr,ipresok)
c         write(6,66)prex(k),tbb(k),ts(k),coskhi,replum,frac,frfr
c 66      format(7e11.3)
 31   continue
      tensok=0.0
c   USING this temperature do the following reprocessing flux
c   calculation
      ttt=0.0
      tens=0.0
      sumsok=0.0
      sum=0.0
cc       f=10**fo
c     
c     frequency loop
cc       do 20 j=1,150
cc       olf=f
cc         f=fo+(j-1)*0.015
cc         f=10.0**f
cc          diff=f-olf
c  summing on all disk for each frequency
cc     do 30 k=npost+2,imax
c         take a mean tempeirature within an annulus
cc                 y(k)=ts(k)
c             t=(y(k)+y(k-1))*0.5
cc                 t=y(k)
cc                dx=x(k)-x(k-1)
cc             hmean=(hh(k)+hh(k-1))*0.5
cc             xmean=(x(k)+x(k-1))*0.5
cc             dh=hh(k)-hh(k-1)
cc               da=sqrt(dx**2+dh**2)
cc             z=c*f/t
c    USE BLACK BODY 
cc                 if(z.gt.80.0) then 
cc                 bb=0.0
cc                 else
cc             bb=2.*pi*h*(f/cv)**2*f/(exp(z)-1.0)
cc                 endif
cc              dx=x(k)-x(k-1)
c THE FOLLOWING LINE IS TO BE REMOVED IF NO COLD DISK IS ADDED
c                 if (k.eq.imax)da=dx
cc       s(k)=2.0*pi*xmean*da*bb
c         tens=tens+f*s(k)*(1.-convx/xmean)**2
c         tens=tens+f*s(k)*(1.-convx/xmean)**2
cc         ttt=ttt+s(k)*(1.-convx/xmean)**1.5
cc 30      continue
cc       tttall=ttt
cc      write(20,22)log10(f),log10(tttall)
cc22      format(2e15.5)
cc        tens=0.0
cc          ttt=0.0
cc 20   continue        
      end

      subroutine sokdsk(i,emitflux,rep,replum,coskhi)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
c        dimension ts(3000)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &  sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /presok/hdd,prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
c     the reprocessing site is actually the grid centre
      xss=xs*rs
      dx=-(prex(i)-prex(i-1))*rs
      xd=rs*(prex(i)+prex(i-1))*0.5
c         print *,'dx,xd',dx/rs,xd/rs
c                     hd=(hh(i)+hh(i-1))*0.5
c                     dh=hh(i)-hh(i-1)
      da=dx
      ratha=0.0
      ratxa=1.0
c     emitflux is not the flux yet. One has to multiply
c   by cos(khi) where khi is the local angle with normal 
c on the shock surface.
c In reality one should use different coskhi for different annulus
c this is used (not fixed coskhi) in the sok2dsk file. This
c should also change absorb with cosine angle. 
      emit=emitflux/(PI)
      dh=hs-hd
      coskhi=(xd-xss)/sqrt((xd-xss)**2+0.25*dh**2)
      ddh=dh/20.0
      rep=0
      replum=0.0
      sokar=0.0
c       the shock height is divided into 100 grids
c here phi0 is that phi for which psi=pi/2
      phi0=acos(xss/xd)
      dphi=phi0/20.0
      f1=dphi*xss*emit
      do 101 j=1,20
         hss=hd+ddh*(j-0.5)
         f2=f1*ddh
c      phi_0 is the maximum angle over which the
c      solved from cospsi(loospage)formula with cospsi=pi/2
ccc        the0=asin(xss/sqrt(xd**2+(hd-hss)**2))
c         print *,'the0',the0
ccc        phi0=pibi2-the0
c       phi0 is calculated differently
         do 100 k=1,41
c      phi coordinate at the shock
            phis=-float(k-21)*dphi
            if(((k.eq.1).or.(k.eq.41)).or.((j.eq.1).or.(j.eq.20))) then
               facto=0.5
            else
               facto=1.0
            endif
c  arbitrary geometry factor of 1.0
            facto=facto*0.5
            cosphi=cos(phis)
            sinphi=sin(phis)
            top=(xd-xss*cosphi)*cosphi-xss*sinphi**2
            bottom=sqrt((xd-xss*cosphi)**2+(hd-hss)**2+(xss*sinphi)**2)
c       print *,bottom
c      the angle between the local normal and the direction
c of reprocessing site from the emitting site
            cospsi=abs(top/bottom)
            top2=-(xd-xss*cosphi)*ratha+(hd-hss)*ratxa
c      the angle between the local normal at the reprocessing site
c with the direction of the ray from the emitting site.
            cosxi=abs(top2/bottom)
            f3=f2*facto/bottom**2*cospsi*cosxi
c   flux falling on the disk which we need to use for
c calculation of the fractional interception at each
c radius
c       replum=replum+facto*emit*ddh*dphi*xss/bottom**2*cospsi*cosxi
            replum=replum+f3
c        sokar=sokar+ddh*dphi*xss
c        print *,sokar
c       reprocessed flux 
            rep=replum
 100     continue
 101  continue
c       print *,'replum before',replum
      anularea=2.*pi*xd*da
c       print *,'anularea',anularea
      replum=replum*anularea
c     print *,'replum after',replum
      return
      end
C ====>

      subroutine dsk2sok(bbtens,cept,sumtem,sumthe,cept2)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /presok/hdd, prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
      stephan=7.56e-15
      sumflux=0.0
      bbtens=0.0
      appaes=0.4
      sigm=5.67e-5
c     height of the sok
      hs= sokheit*rs
c     height of the disk
      hd=hdd*rs
c      CALCULATION OF TS FOR reprocessing AT EACH RADIUS in the 
c      preshock flow
      sumn=0.0
      tensn=0.0
      sokarea=2.*pi*xs*(sokheit-hdd)*rs**2
      idisk=int(xout-xs)
c     print *,'idisk', idisk
      do 34 i=1,idisk+1
         rdisk(i)=(xs+1.0*(i-1))
 34   continue         
      xss=xs*rs
      dh=hs-hd
      ddh=dh/10.0
      cept=0.0
      cept2=0.0
      sumtem=0.0
      sumthe=0.0
      do 101 j=1,idisk
         dx=(rdisk(j+1)-rdisk(j))*rs
         xdd=(rdisk(j+1)+rdisk(j))*0.5
         xd=rs*xdd
c here phi0 is that phi for which psi=pi/2
         phi0=acos(xss/xd)
c       print *,'phi0',phi0
         dphi=phi0/20.0
         if(iter.eq.0) then
            dskflux0(j)=sscons/xdd**3*(1.-sqrt(3.0/xdd))
            dskflux(j)=dskflux0(j)
         endif
cc       tsbb(j)=(dskflux(j)/sigm/emdotdsk)**0.25
         tsbb(j)=(dskflux(j)/sigm)**0.25
         sumflux=sumflux+dskflux(j)
c       print *,dskflux(j),xd,dx
         ssemit=dskflux(j)/pi
         darea=2.*pi*xd*dx
         tens=dskflux(j)*darea
         tanthe=0.5*sokheit/(xdd-xs) 
         if (j.lt.3) then
cx        write(9,18)xdd,dskflux(j),sumflux,ssemit,tens,tanthe
cx 18         format (6e13.4)
         endif
c       costhe=cos(atan(tanthe))
c       call the subroutine to sum over all the contributions
c from the disk, at the reprocessing site on the shock.
         call dsksok(j,dx,xss,ssemit,sumflux,bbtens,repsok,
     &        repint,phi0,dphi)
         d2ss=repsok/tens
c          if(j.gt.3) then
         bbtens=bbtens+tens
         sumthe=sumthe+tanthe*darea*d2ss
         sumtem=sumtem+d2ss*darea
c   here the intensity has to computed so the cos(xi) will not
c be present. The photons are scattered within the
c medium. When rad. transfer is absent, cos(xi) is necessary.
         cept=cept+repint
c  cept2 is to find the fractional interception by the shock
c from the disk.
         cept2=cept2+repsok
c           endif
c          d2s(j)=repsok/sokarea
         frac=cept2/bbtens
cx         write(6,66)xdd,d2ss,tsbb(j),cept,cept2,frac
cx 66      format(6e13.5)
 101  continue
      tensok=0.0
c   USING this temperature do the following reprocessing flux
c   calculation
      ttt=0.0
      tens=0.0
      sumsok=0.0
      sum=0.0
cc       f=10**fo
c     
c     frequency loop
cc       do 20 j=1,150
cc       olf=f
cc         f=fo+(j-1)*0.015
cc         f=10.0**f
cc          diff=f-olf
c  summing on all disk for each frequency
cc     do 30 k=npost+2,imax
c         take a mean tempeirature within an annulus
cc                 y(k)=ts(k)
c             t=(y(k)+y(k-1))*0.5
cc                 t=y(k)
cc                dx=x(k)-x(k-1)
cc             hmean=(hh(k)+hh(k-1))*0.5
cc             xmean=(x(k)+x(k-1))*0.5
cc             dh=hh(k)-hh(k-1)
cc               da=sqrt(dx**2+dh**2)
cc             z=c*f/t
c    USE BLACK BODY 
cc                 if(z.gt.80.0) then 
cc                 bb=0.0
cc                 else
cc             bb=2.*pi*h*(f/cv)**2*f/(exp(z)-1.0)
cc                 endif
cc              dx=x(k)-x(k-1)
c THE FOLLOWING LINE IS TO BE REMOVED IF NO COLD DISK IS ADDED
c                 if (k.eq.imax)da=dx
cc       s(k)=2.0*pi*xmean*da*bb
c         tens=tens+f*s(k)*(1.-convx/xmean)**2
c         tens=tens+f*s(k)*(1.-convx/xmean)**2
cc         ttt=ttt+s(k)*(1.-convx/xmean)**1.5
cc 30      continue
cc       tttall=ttt
cc      write(20,22)log10(f),log10(tttall)
cc22      format(2e15.5)
cc        tens=0.0
cc          ttt=0.0
cc 20   continue        
      end

      subroutine dsksok(i,dx,xss,
     &     ssemit,sumflux,bbtens,repsok,repint,phi0,dphi)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
c     the reprocessing site is actually the grid centre
c        anularea=2.*pi*xss*ddh
      repsok=0.0
      repint=0.0
      f1=dphi*dx*ssemit*xd
      do 31 j=1,10
         hss=hd+ddh*(j-0.5)
c     print *,'hss',hss/rs
         do 100 k=1,41
            phid=-float(k-21)*dphi
            if(((k.eq.1).or.(k.eq.41)).or.((j.eq.1).or.(j.eq.10))) then
               facto=0.5
            else
               facto=1.0
            endif
            cosphi=cos(phid)
            sinphi=sin(phid)
            top=hss-hd
            bottom=sqrt((xss-xd*cosphi)**2+top**2+(xd*sinphi)**2)
c      the angle between the local normal and the direction
c of reprocessing site from the emitting site
            cospsi=abs(top/bottom)
            top2=(xss-xd*cosphi)
c      the angle between the local normal at the reprocessing site
c with the direction of the ray from the emitting site.
            cosxi=abs(top2/bottom)
c   flux falling on the disk
c       replum=replum+facto*emit*ddh*dphi*xss/bottom**2*cospsi*cosxi
            f4=f1*facto/bottom**2*cospsi
            f3=f4*cosxi
            repsok=repsok+f3
c  repint is used to compute incident intensity on the shock.
            repint=repint+f4
c        sokar=sokar+ddh*dphi*xss
c        print *,sokar
c       reprocessed flux (with the cos xi factor)
c       rep=rep+f3*cosxi
 100     continue
 31   continue
      anularea=TWOPI*xss*ddh
c    arbitrary geometric factor of 1.0
      repsok=1.0*repsok*anularea
      repint=1.0*repint*anularea
      return
      end
      subroutine thialflo(sinss,absorb)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
c      dimension xx(7000),tt(7000),tel(7000)
      common /dat/ yy(2),dydt(2)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /presok/ hdd,prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /alfabet/ x0,alfax0,rdel,temp,t5
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
c     alpha=0.2
c            coskhi=0.17
c            telmean=200.0
      tekelvin=telmean*KEVTOERG/BOLTZ
      theta=BOLTZ*tekelvin/EME/CLIGHT**2
      dcapx=2./50.0
      amu0=sinss
c     phi0=1.+2.*amu0
      zstar=7.81/511.0
      deno=1./(1.-alpha)*(2.*theta)**(1.-alpha)
      sum=0.0
      dzzz=theta*dcapx
      do 10 i=1,50
         capx=dcapx*i
         zzz=theta*capx
         rat=(zstar/zzz)**3
         delnu=rat/(1.+rat)
c       if(T.gt.5e+5) only recoil effect is important
         photo=0.0
         argu=sqrt(delnu/zzz)
         if(argu.gt.4.0) then
            asym=0.0
            sq1mlamda=zzz/sqrt(PI*delnu)+photo*sqrt(delnu)
         else
            erfx=1.-1./(1.+0.0705230784*argu+0.0422820123*argu**2+
     &           0.0092705272*argu**3+0.0001520143*argu**4+0.0002765672*
     &           argu**5+0.0000430638*argu**6)**16
            erfun=1.0-2./sqrt(PI)*erfx
            if(erfun.le.0.0) erfun=1.e-6
            sq1mlamda=(sqrt(zzz)*exp(delnu/zzz)*erfun+photo*sqrt(delnu)) !delta
         endif
         amda=1.-sq1mlamda**2
         hhh=(1.+sqrt(3.)*amu0)*(1.-amda*0.25*(1.+amda**2)*amu0*
     &        (log(amu0)+1.33-1.458*amu0**0.62))
     &        /(1.+sqrt(3.)*sq1mlamda*amu0) !phi(mu0)
c         print *,'hhh', hhh
c            anu=1.-hhh*(sqrt(zzz)*exp(delnu/zzz)*erfun+sqrt(delnu))
         anu=1.-hhh*sq1mlamda   !albedo
         sum=sum+(1.-anu)*zzz**(-alpha)*dzzz
c        print *,'delnu,zzz,argu,sum,anu', delnu,zzz,argu,sum,anu
 10   continue
      absorb=sum/deno
c            print *, 'High T absorb, alpha <1',absorb
      end

c      function fact(j)
c      implicit double precision(a-h,o-z)
c      fac=1
c      if(j.eq.0) then
c         fac=1
c      else
c         do 10 i=1,j
c            fac=i*fac
c 10      continue
c      endif
c      fact=fac
c      return
c      end

c      function prod(m)
c      implicit double precision(a-h,o-z)
c      pro=1.0
c      do 10 j=1,m
c         fac=(2.*j-1.)
c         pro=pro*fac
c 10   continue
c      prod=pro
c      return
c      end

      subroutine tloalflo(sinss,absorb)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
c      dimension xx(7000),tt(7000),tel(7000)
      common /dat/ yy(2),dydt(2)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /presok/ hdd,prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /alfabet/ x0,alfax0,rdel,temp,t5
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
c        alpha=0.2
c            coskhi=0.17
c            telmean=200.0
      tekelvin=telmean*KEVTOERG/BOLTZ
      theta=BOLTZ*tekelvin/EME/CLIGHT**2
      dcapx=2./50.0
      amu0=sinss
c            phi0=1.+2.*amu0
      zstar=7.81/511.0
      deno=1./(1.-alpha)*(2.*theta)**(1.-alpha)
      sum=0.0
      dzzz=theta*dcapx
      do 10 i=1,50
         capx=dcapx*i
         zzz=theta*capx
         rat=(zstar/zzz)**3
         delnu=rat/(1.+rat)
         photo=1.0
         argu=sqrt(delnu/zzz)
         if(argu.gt.4.0) then
            asym=0.0
            sq1mlamda=zzz/sqrt(PI*delnu)+photo*sqrt(delnu)
         else
            erfx=1.-1./(1.+0.0705230784*argu+0.0422820123*argu**2+
     &           0.0092705272*argu**3+0.0001520143*argu**4+0.0002765672*
     &           argu**5+0.0000430638*argu**6)**16
            erfun=1.0-2./sqrt(PI)*erfx
            if(erfun.le.0.0) erfun=1.e-6
            sq1mlamda=(sqrt(zzz)*exp(delnu/zzz)*erfun+photo*sqrt(delnu))
         endif
         amda=1.-sq1mlamda**2
         hhh=(1.+sqrt(3.)*amu0)*(1.-amda*0.25*(1.+amda**2)*amu0*
     &        (log(amu0)+1.33-1.458*amu0**0.62))
     &        /(1.+sqrt(3.)*sq1mlamda*amu0)
c         print *,'hhh', hhh
c            anu=1.-hhh*(sqrt(zzz)*exp(delnu/zzz)*erfun+sqrt(delnu))
         anu=1.-hhh*sq1mlamda
         sum=sum+(1.-anu)*zzz**(-alpha)*dzzz
c        print *,'delnu,zzz,argu,sum,anu', delnu,zzz,argu,sum,anu
 10   continue
      absorb=sum/deno
c            print *, 'low T absorb alpha < 1',absorb
      end
      subroutine tloalfhi(sinss,absorb)
c      This is for albedo at low temperature when alfa is higher than 1         
      implicit double precision(a-h,o-z)
      include 'ct.inc'
c     dimension xx(7000),tt(7000),tel(7000)
      common /dat/ yy(2),dydt(2)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /presok/ hdd,prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /alfabet/ x0,alfax0,rdel,temp,t5
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
c        alpha=0.2
c            coskhi=0.17
c            telmean=200.0
      tekelvin=telmean*KEVTOERG/BOLTZ
      theta=BOLTZ*tekelvin/EME/CLIGHT**2
      dcapx=20.0/500.0
      amu0=sinss
      deno=PI**4/15.0
      sum=0.0
      zstar=7.81/511.0
      do 10 i=1,500
         capx=dcapx*i
         zzz=theta*capx
         rat=(zstar/zzz)**3
         delnu=rat/(1.+rat)
         photo=1.0
         argu=sqrt(delnu/zzz)
         if(argu.gt.4.0) then
            asym=0.0
            sq1mlamda=zzz/sqrt(PI*delnu)+photo*sqrt(delnu)
         else
            erfx=1.-1./(1.+0.0705230784*argu+0.0422820123*argu**2+
     &           0.0092705272*argu**3+0.0001520143*argu**4+0.0002765672*
     &           argu**5+0.0000430638*argu**6)**16
            erfun=1.0-2./sqrt(PI)*erfx
            if(erfun.le.0.0) erfun=1.e-6
            sq1mlamda=(sqrt(zzz)*exp(delnu/zzz)*erfun+photo*sqrt(delnu))
         endif
         amda=1.-sq1mlamda**2
         hhh=(1.+sqrt(3.)*amu0)*(1.-amda*0.25*(1.+amda**2)*amu0*
     &        (log(amu0)+1.33-1.458*amu0**0.62))
     &        /(1.+sqrt(3.)*sq1mlamda*amu0)
c         print *,'hhh', hhh
         anu=1.-hhh*sq1mlamda
         sum=sum+(1.-anu)*capx**3/(exp(capx)-1.)*dcapx
c        print *,'delnu,zzz,argu,sum,anu', delnu,zzz,argu,sum,anu
 10   continue
      absorb=sum/deno
c            print *, 'low T absorb alpha>1',absorb
      end
      subroutine thialfhi(sinss,absorb)
c      This is for albedo at high temperature when alfa is higher than 1         
      implicit double precision(a-h,o-z)
      include 'ct.inc'
c      dimension xx(7000),tt(7000),tel(7000)
      common /dat/ yy(2),dydt(2)
      common /diskrad/ rdisk(1000),dskflux(1000),tsbb(1000),d2s(1000),
     &     sokarea,sokheit,sscons,sigm,emdotdsk,emdot0,idisk,iter
      common /cons/ogamda,conp,ej,te,rs,den0,tauhori,xd
      common /presok/ hdd,prex(1000),pret(1000),ipresok
      common /aa/ r,dr,tout,xout,conep,conpro,clina,ucons,conheit,denall
      common /ab/ rho0,emdot,gamhi,gamlo,embh,enr
      common /compton/ acomp,facalpha,d0alf,ssflux,elbyel0
      common /amdas/ amdacomp,amdab,amdaib,gamaep,tauvert,gamterm
      common /betalfa/ tau00,telmean,tinject,tss,xs,hd,hs,ddh,dh,alpha
      common /alfabet/ x0,alfax0,rdel,temp,t5
      common /tauu/ funtau,c1,absorbhi,absorblo,dskflux0(1000),fracvec
      common /paraIn/ cr
c        alpha=0.2
c            coskhi=0.17
c            telmean=200.0
      tekelvin=telmean*KEVTOERG/BOLTZ
      theta=BOLTZ*tekelvin/EME/CLIGHT**2
      dcapx=20.0/500.0
      amu0=sinss
      deno=PI**4/15.0
      sum=0.0
      zstar=7.81/511.0
      do 10 i=1,500
         capx=dcapx*i
         zzz=theta*capx
         rat=(zstar/zzz)**3
         delnu=rat/(1.+rat)
c       just the recoil effect at high temperature
         photo=0.0
         argu=sqrt(delnu/zzz)
         if(argu.gt.4.0) then
            asym=0.0
            sq1mlamda=zzz/sqrt(PI*delnu)+photo*sqrt(delnu)
         else
            erfx=1.-1./(1.+0.0705230784*argu+0.0422820123*argu**2+
     &           0.0092705272*argu**3+0.0001520143*argu**4+0.0002765672*
     &           argu**5+0.0000430638*argu**6)**16
            erfun=1.0-2./sqrt(PI)*erfx
            if(erfun.le.0.0) erfun=1.e-6
            sq1mlamda=(sqrt(zzz)*exp(delnu/zzz)*erfun+photo*sqrt(delnu))
         endif
         amda=1.-sq1mlamda**2
         hhh=(1.+sqrt(3.)*amu0)*(1.-amda*0.25*(1.+amda**2)*amu0*
     &        (log(amu0)+1.33-1.458*amu0**0.62))
     &        /(1.+sqrt(3.)*sq1mlamda*amu0)
c         print *,'hhh', hhh
         anu=1.-hhh*sq1mlamda
         sum=sum+(1.-anu)*capx**3/(exp(capx)-1.)*dcapx
c        print *,'delnu,zzz,argu,sum,anu', delnu,zzz,argu,sum,anu
 10   continue
      absorb=sum/deno
c            print *, 'high T absorb alpha > 1',absorb
      end            

      subroutine ctcompspec
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      PARAMETER (MXROW=2000, MXCOL=20, MXVEC=20, MXCMD=30)
      parameter (MLAG = 15)
      dimension param(5)
      common /comab/ rs,xs,tss,ro
      common /combet/ tau00,telmean,tinject,hd,hs,ddh,dh,alpha
      common /comalf/ x0,alfax0,rdel,temp,t5
      common /comneed/com(NPTS),eemm,tin,tell,tao,teff,xshock,schi,xinn
      common /paraIn/ cr
      embh=eemm
      tinject=tin
      telmean=tell
      tau00=tao
      tss=teff
      xs=xshock
      rs=schi
      xin=xinn
c         print *,'give embh(1.e+8),tinject,telmean,tau1,tss'
c         read *, embh,tinject,telmean,tau1,tss
c        tau00=tau1*2.0*0.8724
c          xs=5.0
c        rs=3.e+5*embh*1.e+8
      param(1)=0.0
      param(2)=tinject
      param(3)=telmean
      param(4)=tau00
      param(5)=2.0
c=====>
C 1. Z    ** REDSHIFT
C 2. T0   ** WIEN TEMPERATURE (keV)
C 3. TEMP ** PLASMA TEMPERATURE (KEV)
C 4. TAUP ** PLASMA OPTICAL DEPTH
C 5. APPRX * <= 1.0 DISK
C           >  1.0 SPHERE
C PHOT ** PHOTON SPECTRUM IN PHOTONS/CM/CM/S/KEV
      ro=1.
      bol2i=0.0
      ZFAC = 1.+PARAM(1)
      T0 = PARAM(2)
      TEMP = PARAM(3)
      TAUP = PARAM(4)
      apprx=param(5)
      t5 = temp/511.
c   because bb feed (if wein use 3.0 instead of 2.7)
c dimensionless soft photon energy is
      x0 = 2.7*t0/temp
      t23=taup+(2./3.)
      taue=taup/(1.0+(temp/39.2)**0.86)
c        print *,'taue',taue
      if (apprx.gt.1.0) then 
c BETA FOR SPHERE
         TB1=PI*PI*(1.-exp(-0.7*taup))/3./t23/t23
         TB2= dexp(-1.4*taup)*dlog(4./(3.*taup))
         beta = tb1+tb2
c              write(*,*) 'beta = ',beta
         if (taue.lt.0.01) then 
            rdel = taue/2.
         else
            rdel=1.0-3./taue*(1.-2./taue+2./taue/taue*(1-exp(-taue))) 
         endif
      elseif (apprx.le.1.0) then
c BETA FOR DISK      
         tb1 = PI*PI*(1.-exp(-1.35*taup))/12./t23/t23
         tb2 = 0.45*dexp(-3.7*taup)*dlog(10./(3.*taup))
         beta = tb1+tb2
c       write(*,*) 'beta = ',beta
         if (taue.lt.0.01) then
            rdel = taue/2.
         else
            rdel = 1.0-(1.0-exp(-taue))/taue !hua&tit 95
         endif
      endif
      f0theta=2.5*t5 + 1.875*t5*t5*(1.0-t5)
      gam0=beta/t5
c        write(*,*) 'gam0 = ',gam0
c compute the power-law index
      tal0 = dsqrt(2.25+gam0)-1.5
 500  bola = 1. + (tal0+3.)*t5/(1.+t5)+4.*d00(tal0)**(1./tal0)*t5*t5
      tal = beta/dlog(bola)
c     write(*,*) 'tal0 tal ',tal0,tal
      if (dabs(tal-tal0).le.1.0d-5) then
         alfax0 = tal
      else
         tal0=tal
         go to 500
      endif
c      print *, 'alfax0 inside comspec= ',alfax0
      call compt
      return
      end

      function d00(x)
      implicit double precision(a-h,o-z)
      x1=x+1.d0
      x2=x+2.d0
      x3=x+3.d0
      bol = x3*x + 4.d0
      argg = 2.d0*x1
      db = zgamln(argg)
      if (db.lt.-70.) then 
         d00=0.0
      elseif (db.gt.70.) then 
         d00=1.e34
      else
         d00 = 3.d0*bol*dexp(db)/x3/x2/x2
      endif
      return
      end

      subroutine compt
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      PARAMETER (MXROW=2001, MXCOL=20, MXVEC=20, MXCMD=30)
      parameter (MLAG = 15)
c      dimension  IERY(MXVEC)
c      dimension ens(mxrow), phot(mxrow), ear(2001)
      dimension y(mxrow,mxcol) 
c      dimension param(4), photar(2001)
      DIMENSION Xt(2000)
c      DIMENSION Fr(2000), Al(2000), Ak(2000), Aj(2000)
c      DIMENSION m1m(2000), mu(2000)
c      DIMENSION F0(2000),ft(2000),gamx(2000),BA(5)
c      dimension taut(5)
      common /comab/ rs,xs,tss,ro
      common /combet/ tau00,telmean,tinject,hd,hs,ddh,dh,alpha
      common /comalf/ x0,alfax0,rdel,temp,t5
      common /comneed/com(NPTS),eemm,tin,tell,tao,teff,xshock,schi,xinn
      common /paraIn/ cr
c
      aa=3./x0
      argg = 3.+ alfax0 
      arg1=argg
      xrfac = rdel**(1./argg)  
      algx0 =2.0*alfax0 + 3.0
      al3x0 = alfax0*(alfax0+3.)/2./algx0
c     at0 corresponds to 10**13.5Hz (DEPENDS UPON telmean)
c      atn corresponds to 10**21Hz (depends upon telmean)
c        at0=-14.0
c        atn=3.0
      flo=13.5                  !freqnc lw
      at0=log(10**13.5*HPLANCK/KEVTOERG/telmean)
      atn=log(10**21.0*HPLANCK/KEVTOERG/telmean)
      hhh=(atn-at0)/npts
      do 190 i = 1,npts
c      ens(i) = ear(i)*zfac
         at=at0+(i-1)*hhh
         xt(i)=dexp(at) 
c            if (ens(i).gt.0.0) then
c            xt(i) = ens(i)/temp
         x=xt(i)
         y(i,1)=xt(i)*temp
         z = xt(i)*t5
         xr=x*xrfac     
         bol3 = x*aa 
         bol3g=bol3
c see if fx is going to overflow
         tfx = (+1.-alfax0)*dlog10(bol3)-x*0.4343+dlog10(algx0)
         if (dabs(tfx).lt.27.) then
            fx = bol3**(-alfax0+1.)*dexp(-x)*algx0
         else
            fx = 0.0
         endif
         bol2i = zgammi(arg1,bol3g)*exp(zgamln(arg1))
         bol7i=zyyit2(xr,alfax0,ro)*fx*bol2i
         bol6= bol3**2.*dexp(-bol3)/(alfax0+3.d0)
c           y(i,3) = al3x0*(bol7i+bol6)*aa
         y(i,8) = al3x0*(bol7i+bol6)
c         phot(i) = y(i,3)
c       else
c            phot(i) = 0.0
c       endif
 190  CONTINUE
c       anueff=1.38e-16*tss/6.625e-27
      tel=telmean*KEVTOERG/BOLTZ
c        cept=0.042
c        factor=6.625e-27*5.7e-5*tss**4/1.38e-16/tel*pi
c     &*(xs*rs)**2*anueff*(tel/tss)
      anueff=2.7*BOLTZ*tss/HPLANCK

!!!    Condition for R < 1.02 ... pure BB

      if(cr.gt.1.02)then
         factor=4.0*BOLTZ**3*tss**3/CLIGHT**2/HPLANCK**2*anueff*PI
     &        *((xs*rs)**2 -(xinn*rs)**2)
      else
         factor=0.d0
      endif

!        print*,'cr & factor :',cr,factor

      facx=KEVTOERG/HPLANCK
c       print *, 'factor,facx,tss,anueff'
c      print *, factor,facx,tss,anueff
      do 291 i=1,npts
c     photar(i) = 0.5*(phot(i)+phot(i-1))
         com(i)=y(i,8)*factor
c          write(8,88) log10(y(i,1)*facx),log10(com(i))
c 88      format(2e15.5)
 291  CONTINUE
c====>
      END
c
c************************************************************************     
C ====>
      FUNCTION ZYYIT2(x,alfa,ro)
      implicit double precision(a-h,o-z)
      DIMENSION w(10) , z(10)
      DATA z/.1377934705 , .7294545495 , 1.8083429017 , 3.4014336978 , 
     &     5.55249614 , 8.3301527467 , 11.8437858379 , 16.2792578313 , 
     &     21.9965858119 , 29.9206970122/
      DATA w/3.0844111576D-01 , 4.0111992915D-01 , 2.180682876E-01 , 
     &     6.208745609E-02 , 9.501516975E-03 , 7.530083885E-04 , 
     &     2.825923349E-05 , 4.249313984E-07 , 1.839564823E-09 , 
     &     9.911827219E-13/
      a2 = alfa + 3.
      a2n = alfa + 2.
      a3 = alfa -1.
      argg=a2+a3+2.
      db = zgamln(argg)
      ZYYIT2 = 0.d0
      DO 50 i = 1 , 10
         v = w(i)*DEXP(a2n*DLOG(ro*x+z(i))+alfa*Dlog(z(i))-db)
         v1 = v/alfa*((x + z(i)) - a2)
         ZYYIT2 = ZYYIT2 + v1
 50   CONTINUE
      ZYYIT2 = ZYYIT2
      RETURN
      END

      DOUBLE PRECISION FUNCTION ZGAMMI(A,X)
      use fgsl
      implicit none
      DOUBLE PRECISION A, X
      IF (X.LT.0.0.OR.A.LE.0.0) pause
      ZGAMMI=fgsl_sf_gamma_inc_p(A,X)
      RETURN
      END


      subroutine bbspec(eemm,numanul)
      implicit double precision (a-h,o-z) 
      include 'ct.inc'
      common /bbneed/ bb(NPTS),xd(MXANUL),tb(MXANUL),ta(3000)
      common /paraIn/ cr
c        character*80 dumy
c        print *,'embh(1.e+8)'
c        read *, embh

      embh=eemm 
      ta(1)=ta(2)
      idisk=numanul
      rs=3.0e+5*embh*1.e+8
      hdd=0.0
      opacity=0.4
      ogamda=20.0
      gamhi=5./3.
      gamlo=4./3.
      te=EME*CLIGHT**2/BOLTZ
      ej=1.44e-27
      capc=9.3e+12
      c=HPLANCK/BOLTZ
      flo=13.5
      f=10**flo
      fhi=21.0
      dfreq=(fhi-flo)/DFLOAT(NPTS)
c        idisk=95
      do 29 k=1,idisk-1
c        read(3,89)dumy
c        read(3,*,end=99) xd(k),tb(k),ta(k),g,g,g,g
c         print *, xd(k),ta(k)
c89      format(a80)
         xd(k)=(xd(k)+0.5)*rs
 29   continue
      do 64 mm=1,NPTS
         fo=flo+dfreq*(mm-1)
         olf=f
         f=10**fo
         dx=rs
         tendisk=0.0
         do 30 k=1,idisk-1
            z=c*f/ta(k)
c    USE BLACK BODY 
            if(z.gt.80.0) then 
               bbb=0.0
            else
               bbb=TWOPI*HPLANCK*(f/CLIGHT)**2*f/(exp(z)-1.0)
            endif
            surfac=TWOPI*xd(k)*dx*bbb
c        write(21,12) xd(k),dx,z,f,ta(k)
c 12         format(5e15.5)
            tendisk=tendisk+f*surfac
 30      continue
         if(tendisk.le.0.0) tendisk=1.d-20
         bb(mm)=tendisk 
c      write(20,22)log10(f),log10(tendisk)
c 22      format(2e15.5)
 64   continue
      return
      end

      subroutine storall(k,prex,tbb,ts,coskhi,replum,frac,frfr,ipresok)
      implicit double precision(a-h,o-z)
      include 'ct.inc'
      common/resu/ prex1(1000),tbb1(MXANUL),ts1(3000),coskhi1(550),
     & replum1(550),frac1(550),frfr1(550)
      common /paraIn/ cr
c     do 10 k=ipresok,2,-1
      prex1(k)=prex
      tbb1(k)=tbb
      ts1(k)=ts
      coskhi1(k)=coskhi
      replum1(k)=replum
      frac1(k)=frac
      frfr1(k)=frfr
c10    continue
      end
