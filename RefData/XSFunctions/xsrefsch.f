      subroutine xsrefsch(ear,ne,param,ifl,photar,photer)

      implicit none
      integer ne,ifl
      real ear(0:ne),param(*),photar(ne),photer(ne)

c     Driver for angle-dependent reflection from an exponentially-cutoff
c     power law and ionized medium.
c     It needs to be compiled together with xspexrav.f (neutral reflection).
c
c     See Magdziarz & Zdziarski, 1995, MNRAS.
c     See Zdziarski et al., 1995, ApJL January 10 for description of 
c     calculation of ionization (based on Done et al. 1992).
c     The abundances are of Morrison & McCammon.
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
c     This version allows for changes of the vector 'ear' between subsequent
c       calls.
c     and then calculates the broanening in a schwartzchild metric
c
c     number of model parameters: 13
c      1: Gamma, power law photon index, N_E prop. to E^{-Gamma}
c      2: E_c, the cutoff energy in keV (if E_c=0 there is no cutoff;
c         one needs to change the lower limit for that)
c      3: scale, scaling factor for reflection; if <0, no direct component
c         (scale=1 for isotropic source above disk)
c      4: redshift, z
c      5: abundance of elements heavier than He relative to 
c         the solar abundances 
c      6: iron abundance relative to the above  
c      7: inclination angle
c      8: disk temperature in K
c      9: disk ionization parameter = 4 pi F/n
c     10: power law index for emissivity
c     11: inner radius
c     12: outer radius
c     13: internal accurancy: points per decade


c parameters number for the model used

      integer nparam,nref,ndisk
      parameter(nparam=12,nref=9,ndisk=6)

c internal integration parameters
c   xmin-xmax            -- energy range (spectrum array)
c   nx_max               -- max number of data points
c   disk_xmin, disk_xmax -- relative energy range for line convolution
c   n_disk               -- max number of line points

      real xmin,xmax      
      real disk_xmin, disk_xmax
      integer nx_max, n_disk
      parameter(xmin=.01,xmax=2000.,nx_max=10000)
      parameter(disk_xmin = 0.5, disk_xmax = 1.5, n_disk = 400)

      logical firstflag,fl0,fl1
      integer i,j,k,nx,jlo,jhi,i1kev
      real x(0:nx_max),x_disk(0:n_disk),dx
      real phot_ref(nx_max),phot_err(nx_max),phot_rs(nx_max)
      real phot_disk(n_disk),phot_derr(n_disk)
      real logxmin, logxmax, logdmin, logdmax
      real pp,pn,r,xlo,xhi,norm
      real in_ref(nref),in_disk(ndisk)
      real cparam(13)
      real sppsc 
      save firstflag,x,nx,in_ref,phot_ref,phot_disk,in_disk,
     >     phot_rs,x_disk,i1kev
      data firstflag/.true./
      data nx/nx_max/
      data in_ref/nref*9999./,in_disk/ndisk*9999./

c this model has no errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      do i=1,13
        cparam(i) = param(i)
      end do
      cparam(7) = cos(0.0174532925*param(7))

      fl0 = .FALSE.
      fl1 = .FALSE.
      
c  nx - actual number of data points in energy array

      pp=min(float(nx_max),cparam(13)*log(xmax/xmin))
      if(nx.ne.int(pp))then
         nx=int(pp)
         firstflag=.true.
         fl0=.true.
         fl1=.true.
      endif
      
      logxmin = log10(xmin)
      logxmax = log10(xmax)
      logdmin = log10(disk_xmin)
      logdmax = log10(disk_xmax)
 
      if(firstflag)then
         i1kev = 1
         dx=log10(xmax/xmin)/nx
         x(0)=xmin
         do i=1,nx
            x(i)=xmin*10**(i*dx)
            if(x(i) .le. 1.0) i1kev = i
         enddo
         dx = log10(disk_xmax / disk_xmin) / n_disk
         x_disk(0) = disk_xmin
         do i = 1, n_disk
           x_disk(i) = disk_xmin * 10**(i * dx)
         enddo
      endif

      do i=1,nref
         if(in_ref(i).ne.cparam(i))then
            fl0=.true.
            in_ref(i)=cparam(i)
         endif 
      end do

      if(fl0)then
         call xspexriv(x,nx,in_ref,ifl,phot_ref, phot_err)
         do i=1,nx
            phot_ref(i) = phot_ref(i) / (x(i) - x(i-1))
         enddo
      endif

      if(in_disk(2).ne.cparam(10))then 
         in_disk(2)=cparam(10)
         fl1=.true.
      endif
      if(in_disk(3).ne.cparam(11))then
         in_disk(3)=cparam(11)
         fl1=.true.
      endif
      if(in_disk(4).ne.cparam(12))then
         in_disk(4)=cparam(12)
         fl1=.true.
      endif
      pp=180./3.1415*acos(cparam(7))
      if(in_disk(5).ne.pp)then 
         in_disk(5)=pp
         fl1=.true.
      endif
      if(in_disk(6).ne.cparam(4))then
         in_disk(6)=cparam(4)
         fl1=.true.
      endif
      
      if(fl1) then
        in_disk(1) = 1.0
        call xszdili(x_disk, n_disk, in_disk, ifl, phot_disk,
     &               phot_derr)
         do i=1, n_disk
           if(phot_disk(i) .ge. 1.0e-8) then
             phot_disk(i) = phot_disk(i) / (x_disk(i) - x_disk(i-1))
           else
             phot_disk(i) = 0.0
           endif
         enddo
      endif
      
      if(fl0.or.fl1) then
        do i=1,nx
           phot_rs(i)=0.
        enddo
        do i = 1, nx
          xlo = disk_xmin * x(i)
          xhi = disk_xmax * x(i)
          jlo = int(nx * (log10(xlo) - logxmin) / (logxmax-logxmin))
          jhi = int(nx * (log10(xhi) - logxmin) / (logxmax-logxmin))
          jlo = max(jlo, 1)
          jhi = min(jhi, nx)
          do j = jlo, jhi
            r = x(i) / x(j)
            k = int(n_disk * (log10(r) - logdmin) / (logdmax-logdmin))
            if((k.ge.1) .and. (k.le.n_disk)) then
              phot_rs(j) = phot_rs(j) + phot_ref(i) * phot_disk(k)
            endif
          enddo
        enddo
        norm = phot_rs(i1kev)
     >       + (phot_rs(i1kev + 1) - phot_rs(i1kev)) 
     >       * (1.0 - x(i1kev)) / (x(i1kev + 1) - x(i1kev))
        do i = 1, nx
          phot_rs(i) = phot_rs(i) / norm
        enddo
      endif

      fl0=.false.
      fl1=.false.

      pp=sppsc(ear(0),x,nx,phot_rs)
      do i=1,ne
         pn=sppsc(ear(i),x,nx,phot_rs)
         photar(i)=0.5*(pn+pp)*(ear(i)-ear(i-1))
         pp=pn
      enddo

      return
      end

c ------------------------------------------------------------------- c
 
      real function sppsc(xx,xt,jt,spt)
 
      implicit   none
      integer    jt
      real       xx,xt(0:jt),spt(jt)
 
      integer    il,ih
      real       sp
      save       ih
      data       ih/2/

      if((xx.gt.xt(1)).and.(xx.lt.xt(jt)))then
         if(xx.lt.xt(ih))then
            ih=2
         endif
 10      if(xx.gt.xt(ih))then
            ih=ih+1
            goto 10
         else
            il=ih-1
c            sp=spt(il)+
c     &         (spt(ih)-spt(il))*
c     &         (xx-xt(il))/(xt(ih)-xt(il))
        sp=spt(il)*exp(log(spt(ih)/spt(il))*
     *     log(xx/xt(il))/log(xt(ih)/xt(il))) 

         endif
      else
         sp=0
      endif
 
      sppsc=sp
 
      return
      end
 
c ------------------------------------------------------------------- c
 

      SUBROUTINE xszdili(earp,ne,param,ifl,photar,photer)

      implicit none
      INTEGER ne,ifl
      REAL earp(0:ne),param(*),photar(ne),photer(ne)
      real ear(0:10000)

c  model to calculate the line shape for a rotating accretion
c  disk. does not include GR effects. note that if param(2) is
c  set to 10 then do the special case of an accretion disk
c  emissivity.
c  modified 11/30/95 by P.Zycki to include object's redshift
c  corrected 20/02/96 by P.Magdziarz 

c  parameters :
c	1        line energy
c       2        power law index for emissivity (10 for disk)
c       3        inner radius (GM/c**2)
c       4        outer radius (GM/c**2)
c       5        inclination  (degrees)
c       6        redshift

      INTEGER nr
      PARAMETER (nr=100)

      INTEGER n, k, ie
      REAL ri, ro, dlogr, dera, alp,
     &     tgcsi2, cosdel, fra, ra, dra,
     &     rd(nr), enel, spm, rilog10, zfac  
      REAL beta, betal, betah, zpo, zpol, zpoh
      REAL sincl, cincl, rafact, flux, total
      REAL enobl, enobh

      DATA dera /57.295779/

c suppress a warning message from the compiler
      ie = ifl

c this model has no errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

c convert input parameters - inclination set to radians.

      enel = param(1)
      alp = param(2)
      ri = param(3)
      rilog10 = alog10(ri)
      ro = param(4)
      zfac = 1+param(6)

      sincl = sin(param(5)/dera)
      cincl = cos(param(5)/dera)

      do n=1,ne
         ear(n) = zfac*earp(n)
      end do

c set spectrum to 0

      DO n = 1, ne
         photar(n) = 0.
      ENDDO

c trap case where inner radius is greater than outer

      IF ( ri .GE. ro ) THEN
         CALL xwrite(
     &     'Inner radius > outer radius  -  model is invalid', 5)
         RETURN
      ENDIF

c trap case of inner radius being less than the last stable orbit.

      IF  (ri .LT. 6.) THEN
         CALL xwrite('Inner radius < 6 R_g  -  model is invalid', 5)
         RETURN
      ENDIF

c calculate log step in radius

      dlogr = (alog10(ro)-rilog10)/float(nr-1)

c calculate radii

      rd(1) = ri
      DO n = 2, nr
         rd(n) = 10.**(rilog10+float(n-1)*dlogr)
      ENDDO

c big loop for radii. note that the radius cannot reach 3 else
c the metric goes singular

      DO n = 1, nr - 1

         ra = (rd(n)+rd(n+1))/2.

         dra = rd(n+1) - rd(n)

         rafact = sqrt(1.-3./ra)

c if power-law index is less than ten use to calculate emissivity

         IF (alp.LT.9.9) THEN

            fra = ra**(alp)

c else use the accretion disk emissivity law.

         ELSE

            fra = (1.-sqrt(6./ra))/ra**3

         ENDIF

c loop over azimuthal angles

         DO k = 1, 179, 2

            beta = float(k)/dera
            betal = float(k-1)/dera
            betah = float(k+1)/dera

c calculate mean redshift (1+z = zpo) for the bin

            tgcsi2 = (sincl*sin(beta))
     &               **2/(1.-(sincl*sin(beta))**2)
            cosdel = sincl*cos(beta)
     &               /sqrt((sincl*cos(beta))**2+cincl**2)
            zpo = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact

c and the low and high redshifts for the bin

            tgcsi2 = (sincl*sin(betal))
     &               **2/(1.-(sincl*sin(betal))**2)
            cosdel = sincl*cos(betal)
     &               /sqrt((sincl*cos(betal))**2+cincl**2)
            zpol = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact

            tgcsi2 = (sincl*sin(betah))
     &               **2/(1.-(sincl*sin(betah))**2)
            cosdel = sincl*cos(betah)
     &               /sqrt((sincl*cos(betah))**2+cincl**2)
            zpoh = (1.+cosdel/sqrt(ra*(1.+tgcsi2)-2.))/rafact


c  enobl and enobh are the lower and upper observed energy from
c  this azimuthal and radial bin

            enobl = min(enel/zpol, enel/zpoh)
            enobh = max(enel/zpol, enel/zpoh)
            total = enobh - enobl

c  calculate total emission from this bin

            flux = fra*ra*dra*2./dera/zpo**3

c  find fractions of emission from this bin to place in each energy
c  range.

            IF (enobh .LT. ear(0) .OR. enobl .GT. ear(ne) ) GOTO 10

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
      DO n = 1, ne
         spm = spm + photar(n)
      ENDDO

c write values

      IF ( spm .NE. 0 ) THEN
         DO n = 1, ne
            photar(n) = photar(n)/spm
         ENDDO
      ENDIF

      RETURN
      END
