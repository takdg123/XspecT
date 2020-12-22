      Subroutine optxagnf(ear,ne,param,ifl,photar,photer)

      implicit none

      integer NPAR
      parameter(NPAR=11)

      integer ne, ifl
      real ear(0:ne), photar(ne), param(NPAR), photer(ne)

!------------------------------------------------------------------------------------ 
! This subroutine checks whether any of the parameters have changed, if they have, 
! then it creates a new energy binning scheme with 5000 bins increasing exponentially
! in size.  The input parameters and new binning scheme are used in a call to the 
! subroutine 'asdiskf'.  The returned fluxes are then corrected for redshift and
! shifted back into the original energy binning scheme.
!------------------------------------------------------------------------------------ 

!     ear = Energy array
!     ne = size of energy array (1 more than flux array as gives boundaries)
!     param = array containing model parameters
!     ifl = spectrum number
!     photar = output flux array (ne-1 parameters)

      integer NN
      parameter(NN=5000)

      integer n,i
      real e(0:NN),ph(NN)
      real dloge,zfac,oldpar(NPAR),newemin,newemax
      logical parchange, echange

      REAL fstart(ne), fend(ne)
      integer istart(ne), iend(ne)

      save oldpar,ph,e

c     param(1) mass in solar
c     param(2) luminosity distance
c     param(3) log mass accretion rate in L/LEdd
c     param(4) astar
c     param(5) rcorona/rg             if -ve plot only disc
c     param(6) log10 rout/rg
c     param(7) opt thick kte
c     param(8) opt thick tau              if -ve plot only comp
c     param(9) power law gamma
c     param(10) frac of coronal power in power law. if -ve plot only pl 
c     param(11) redshift

c this model has no errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      zfac = 1.0 + param(11)

!     test whether any parameters have changed
      parchange=.false.
      do i=1,NPAR,1
         if (param(i).ne.oldpar(i)) parchange=.true.
      end do

!     test whether own energy grid needs to be changed
      if (ear(0).eq.0.0) then 
         newemin=ear(1)-(ear(2)-ear(1))/10.0
      else
         newemin=min(1.0e-5,ear(0)*zfac)
      end if
      newemax=max(1.0e3,ear(ne)*zfac)

      if ((e(0).ne.newemin).or.(e(NN).ne.newemax)) then
         echange=.true.
      end if

!     if necessary recalculate the energy grid

      if (echange) then

         dloge=log10(newemax/newemin)/float(NN)
         e(0) = newemin
         do n=1,NN,1
!           populate the new energy array and zero flux
            e(n)=10**(log10(e(0))+dloge*float(n))
         end do

      endif
         

!     call the irradiated disc if parameters or energy grid have changed
      if (parchange .or. echange) then 

         call mydiskf(e,NN,param,ifl,ph)

!        now redshift the energy bins back
         do n=0,NN,1
            e(n)=e(n)/zfac
         end do

!        and now correct the flux for redshift
         do n=1,NN,1
            ph(n)=ph(n)/zfac
         end do

!        save the new parameters into the oldpar array
         do i=1,NPAR,1
            oldpar(i)=param(i)
         end do

      end if

!     rebin the calculated fluxes back onto original energy grid
      CALL inibin(NN, e, ne, ear, istart, iend, fstart, fend, 0)
      CALL erebin(NN, ph, ne, istart, iend, fstart, fend, photar)
     
      return
      end



      subroutine mydiskf(ear,ne,param,ifl,photar)


c     program to integrate the disk equations from shakura-sunyaev disk
c     as given by Novikov and Thorne

      implicit none
      double precision m,mdot,rsg
      double precision astar,z1,z2,rms,fcor,fpl
      double precision corona,cor,pow,rcor,mytemp,fcol
      double precision rgcm,pi,t,t0,r,dr,dlogr,dflux
      double precision en,tot,kkev,h,kevhz,d0,d,logrout
      integer i,iin,imax,n,ne,ifl
      logical first
      double precision flux(ne),disc,discu,mdotedd,alpha,eff

      real ear(0:ne),photar(ne),param(*)
      real lpar(5),lphot(ne),lphote(ne),hpar(5),hphot(ne),hphote(ne)
      


c     constants
      h=6.62617d-27
      kkev=1.16048d7
      kevhz=2.417965d17
      pi=4.0*atan(1.0)      

c     system parameters
      m=dble(param(1))         !in solar units
      d0=dble(param(2))        !in Mpc
      mdotedd=dble(10**(param(3)))      !in L/Ledd if -ve plot disc
      astar=dble(param(4))     
      alpha=0.1

c     get d and rg in cm
      d=d0*1d6*3d18
      rgcm=1.5d5*m
      pi=4.0*atan(1.0)

c     get rms
      z1=((1-astar**2)**(1./3.))
      z1=z1*(((1+astar)**(1./3.))+((1-astar)**(1./3.)))
      z1=1+z1
      z2=sqrt(3*astar*astar+z1*z1)
      rms=3.+z2-sqrt((3.-z1)*(3.+z1+2.0*z2))
      eff=1.0-sqrt(1.0-2.0/(3.0*rms))
c      write(*,*) 'eff= ',eff

c     mdot in g/s Ledd=1.39e38 m and c2=9e20 
      mdot=mdotedd*1.39d18*m/(9.0*eff)

c     disc parameters
      rcor=dble(abs(param(5)))       !in Rin/Rg
      logrout=dble(param(6))   !log10 Rout/Rg

c     iv -ve then calculate rout=rsg frmo laor & netzer 1989
      if (param(6).lt.0.0) then
         logrout=((m/1e9)**(-2.0/9.0)) * ((mdotedd)**(4.0/9.0))
         logrout=2150*logrout* (alpha**(2.0/9.0))
         write(*,*) 'rout= ',logrout
         logrout=log10(logrout)
      end if

      rsg=10**logrout

      fpl=abs(dble(abs(param(10))))  !if -ve plot pl

c     initialise 
      do n=1,ne,1
         photar(n)=0.0
         flux(n)=0.0
      end do
      
      if (rcor.gt.rms) then
         iin=50
      else
         iin=0
         rcor=rms
      end if
      imax=1000
      first=.true.
      corona=0.0
      disc=0.0d0
      discu=0.0d0
      t0=0.0d0

      do i=1,iin+imax,1
         if (i.le.iin) then
            dlogr=log10(rcor/rms)/float(iin-1)
            r=10**(log10(rms)+float(i-1)*dlogr+dlogr/2.0)
            else
            dlogr=log10(rsg/rcor)/float(imax-1)
            r=10**(log10(rcor)+float(i-iin-1)*dlogr+dlogr/2.0)
         end if
         dr=10**(log10(r)+dlogr/2.0) - 10**(log10(r)-dlogr/2.0)



         t=mytemp(m,astar,mdot,rms,r)
         if (t.gt.1.0e5) then
            fcol=(72.0/(t/kkev))
            fcol=fcol**(1.0/9.0)
            else 
            if (t.gt.3.0e4) then
c              linear from 1 to 2.7 between 4e4 and 1e5
               fcol=(t/3.0e4)**0.82
               else
               fcol=1.0
            end if
         end if

         if (i.le.iin) then
            discu=discu+2.0*2.0*pi*r*dr*rgcm*rgcm*5.67e-5*(t**4)
            else
            if (first) then
            	t0=t*fcol/kkev
c                write(*,*) 'seed photon',t,t0*kkev
            	first=.false.
            end if
         end if

         disc=disc+2.0*2.0*pi*r*dr*rgcm*rgcm*5.67e-5*(t**4)
         t=t/kkev

c        go over each photon energy - midpoint of bin
         do n=1,ne,1
            en=dble(log10(ear(n))+log10(ear(n-1)))
            en=en/2.0
            en=10**en

c           do blackbody spectrum
            if ((en.lt.30.0*t*fcol).and.(r.gt.rcor)) then
               dflux=pi*2.0*h*((en*kevhz)**3)/9.0d20
               dflux=dflux*4.0*pi*r*dr*rgcm*rgcm/(exp(en/(t*fcol))-1.0)
               dflux=dflux/(fcol**4)
               else
               dflux=0.0d0
            end if

            flux(n)=flux(n)+(dflux)

          end do 
      end do


      tot=0.0
      do n=1,ne,1
         tot=tot+flux(n)*(ear(n)-ear(n-1))*kevhz
         
c        this is ergs cm^2 s-1 Hz^-1 
         flux(n)=flux(n)/(4.0*pi*d*d)

c        photons is flux/hv - photons cm^2 s-1 Hz^-1
         flux(n)=flux(n)/(h*kevhz*ear(n))

c        now multiply by energy band in Hz
         photar(n) = sngl(flux(n)*kevhz) * (ear(n)-ear(n-1))

      end do


      fcor=discu/(disc-discu)
c      write(*,*) tot,disc,discu,fcort,t0*1.1e5

c     now add low temp comptonised emission with comptt
      lpar(1)=0.0
      lpar(2)=SNGL(t0)
      lpar(3)=param(7)
      lpar(4)=abs(param(8))   !if -ve plot comp
      lpar(5)=1

      call xstitg(ear,ne,lpar,ifl,lphot,lphote)

c     now add high temp compton emission
      hpar(1)=param(9)
      hpar(2)=100.0
      hpar(3)=SNGL(t0)
      hpar(4)=0.0
      hpar(5)=0.0
      call donthcomp(ear,ne,hpar,ifl,hphot,hphote)

      cor=0.0
      pow=0.0
      tot=0.0
      do n=1,ne,1
         cor=cor+lphot(n)*ear(n)
         pow=pow+hphot(n)*ear(n)
         tot=tot+photar(n)*ear(n)
      end do

c     now add in the rest of the components
      do n=1,ne,1
         if ((param(5).lt.0.0).or.(param(8).lt.0.0).or.
     &        (param(10).lt.0.0)) then
            if (param(10).lt.0.0) photar(n)=SNGL(hphot(n)*(tot/pow)
     &                                        *fcor*fpl)
            if (param(8).lt.0.0) 
     &            photar(n)=SNGL(lphot(n)*(tot/cor)*fcor*(1.0-fpl))
            if (param(5).lt.0.0) photar(n)=photar(n)
         else
            photar(n)=SNGL(photar(n)+lphot(n)*(tot/cor)*fcor*(1.0-fpl)
     &           +hphot(n)*(tot/pow)*fcor*fpl)
         end if   
      end do


      return
      end


