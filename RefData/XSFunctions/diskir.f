      subroutine diskir(ear,ne,param,ifl,photar,photer)

      implicit none
      INTEGER NEDISK
      PARAMETER(NEDISK=5000)
      integer ne,ifl
      real ear(0:ne),photar(ne),photer(ne),param(8)
      integer i
      external xwrite
      real e(0:NEDISK),dph(NEDISK),dphe(NEDISK),tph(NEDISK)
      real lx,ld,tpar(5)
      real tin,lcld,rirr,fin,fout
      
c     parameter 1 is intrinsic disc temperature 
      tin=param(1)
      e(0)=1.0e-3
      do i=1,NEDISK,1
         e(i)=e(0)*( 10**(float(i)*6.0/float(NEDISK)) )
         dph(i)=0.0
         dphe(i)=0.0
         tph(i)=0.0
      end do
      call xsdskb(e,NEDISK,tin,ifl,dph,dphe)
      
c     parameter 2 and 3 are gamma and kTe for thcomp
      tpar(1)=param(2)
      tpar(2)=param(3)

      lcld=param(4)
      fin=param(5)
      rirr=param(6)
      if (rirr .le. 1.0) then
        call xwrite('***DISKIR Error: Cannot calculate for rirr <= 1.0'
     +                ,5)
        goto 999
      end if
      fout=param(7)

c     new inner disc temperature
      tpar(3)=param(1) *(  (1.0 + fout*(1.0+lcld*(1.0+fin) )
     &     + 2.0*fin*lcld/(rirr**2-1.0) ) **0.25)

c     blackbody at each radius
      tpar(4)=0.0

c     redshift
      tpar(5)=0.0

      call donthcomp(e,NEDISK,tpar,ifl,tph,dphe)

      lx=0.0
      ld=0.0
      do i=1,NEDISK,1
         lx=lx+tph(i)*e(i)
         ld=ld+dph(i)*e(i)
      end do

      call dthircore(ear,ne,param,ifl,lx,ld,photar,photer)

 999  return 
      end 


      subroutine dthircore(ear,ne,param,ifl,lx,ld,photar,photer)
c     code to do diskbb for an irradiated disc where irradiation is
c     from diskbb+thcompml

      implicit none
      integer ne,ifl,n,i,imax,iin
      double precision r,logr,dlogr,dr,t,tin,dflux,dnphot,en
      real param(8),ear(0:ne),photar(ne),photer(ne)
      real t1,tpar(5),thphotar(ne),lx,ld
      double precision fin,fout,lcld,logrout,rlow,rhigh,rirr
      
      do n=1,ne,1
         photar(n)=0.0
      end do

c     this is tin it would have without irradiation
      tin=dble(param(1))

      lcld=dble(param(4))
      fin=dble(param(5))
      rirr=dble(param(6))
      fout=dble(param(7))
      logrout=dble(param(8))

c     thcomp parameters
      tpar(1)=param(2)
      tpar(2)=param(3)

c     blackbody at each radius
      tpar(4)=0.0

c     redshift
      tpar(5)=0.0

      imax=1000
      iin=50
      logr=0.0
      rlow=1.0d0
      t1 = 0.0
      do i=1,iin+imax,1
         if (i.le.iin) then
            dlogr=log10(rirr)/dble(iin)
            logr=dble(i-1)*dlogr+dlogr/2.0
            else
            dlogr=(logrout-log10(rirr))/dble(imax)
            logr=log10(rirr)+dble(i-iin-1)*dlogr+dlogr/2.0
         end if 
         r=10**logr
         rhigh=10**(logr+dlogr/2.0)
         dr=rhigh-rlow

         if (i.le.iin) then
            t=(tin*(r**(-0.75)))*(  (1.0 + fout*r*(1.0+(1.0+fin)*lcld) 
     &           + 2.0*(r**3)*lcld*fin/(rirr**2-1.0) )**0.25)
            else
            t=(tin*(r**(-0.75)))*
     &              (  (1.0 + fout*r*(1.0+(1.0+fin)*lcld)  )**0.25)
         end if
         if (i.eq.1) t1=sngl(t)

c        go over each photon energy
         do n=1,ne,1
c           this is in ergs cm^-2 s^-1 keV^-1 ie at 10kpc
            en=dble(ear(n))
            if (en.lt.20.0*t) then
               dflux=3.33d-22*(en**3)*r*dr/(exp(en/t)-1.0)
               else
               dflux=0.0d0
            end if

c           divide by energy to get photon spectrum, by need to change keV
c           to ergs. then multiply by de to get photons in the bin
            dnphot=dflux*(en-dble(ear(n-1)))/(1.602177e-9*en)
            
c           rin in km r dr = 1e10 factor
            photar(n)=photar(n)+1.e10*sngl(dnphot)

         end do 
         rlow=rhigh

      end do

c     set inner thcomp temperature to new inner disc temp.
      tpar(3)=t1

      call donthcomp(ear,ne,tpar,ifl,thphotar,photer)
      
      do n=1,ne,1
         photar(n)=photar(n)+sngl(lcld)*(ld/lx)*thphotar(n)
      end do

      return
      end
