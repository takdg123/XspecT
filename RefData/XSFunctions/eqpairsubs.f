c This file contains fortran subroutines for the eqpair, eqtherm, and
c compth models.g

c ===================================


      SUBROUTINE pairiter(param,np,imodel,compspect,
     &                    blbdspect,wkev)

      IMPLICIT NONE

      INCLUDE "common.eq"

      INTEGER np, imodel
      real param(np)
      DOUBLE PRECISION bbspect(NPHENERG)
      DOUBLE PRECISION compspect(NPHENERG), blbdspect(NPHENERG)
      DOUBLE PRECISION wkev(NPHENERG), fnorm, radnorm

      INTEGER i

c set model parameters and initialize model

      call setup(param,np,imodel,radnorm)

c double check for bad input parameters
      if((rle .lt. 1.d-40).and.(rlth .lt. 1.d-40)
     &     .and.(taup.lt.1d-40)) then
c     do nothing!!  (no pairs!, gdist(i) = phinjfn(i) in setup)
         do i=1,NPHENERG
            esc(i) = gdist(i)
         enddo
         GOTO 999
      endif
      if((rle .lt. 1.d-40).and.(rlth .gt. 1.d-40)
     &     .and.(taup.lt.1d-40)) then
         print *,'ILLEGAL INPUT COMBINATION:  l_th is not zero, yet you'
         print *,'have no background electrons/positrons since tau_p=0'
         print *,'and l_nth=0!!     '
         do i=1,NPHENERG
            esc(i) = gdist(i)
         enddo
         GOTO 999
      endif

c get model spectrum

      IF ( imodel .EQ. 1 ) THEN
         call solve(.FALSE.)
      ELSE IF ( imodel .EQ. 2 ) THEN
         CALL solve(.TRUE.)
      ELSE IF ( imodel .EQ. 3 ) THEN
         CALL thcsolve
      ENDIF

c get uncomptonized blackbody spectrum

      call bbesc(bbspect)

c set up output arrays. compspect is output spectrum minus the black body

 999  CONTINUE

      fnorm = radnorm/511.0
      DO i = 1, NPHENERG
         esc(i) = esc(i)/gwidth(i)
         bbspect(i) = max(1.d-30,bbspect(i)/gwidth(i))
         compspect(i) = max(1.d-30,(esc(i) - bbspect(i))*fnorm)
         blbdspect(i) = bbspect(i)*fnorm
         wkev(i) = wgam(i)*511.0
      ENDDO

      return
      end

c ===================================================================
c This routine calculates the photon distribution given the electron
c distribution.
c Global variables:
c    i: ctab, compr, ictab, ppr, ppf, ipp, tesc, edist, posdist
c    i/r: gdist
c    r: gder

c Calls routines dotherm, dobrem, mkscat

      SUBROUTINE dopcomp(isthermonly)

      IMPLICIT NONE
      LOGICAL isthermonly
      INCLUDE "common.eq" 

c Local variables

      DOUBLE PRECISION scat(NPHENERG,NPHENERG)
      DOUBLE PRECISION ntscat(NPHENERG,NPHENERG)
      DOUBLE PRECISION pin(NPHENERG), pout(NPHENERG)
      DOUBLE PRECISION pin0(NPHENERG), pout0(NPHENERG)
      DOUBLE PRECISION thresh, dmax, erate, odmax, rate, tau, temp

      INTEGER iscat(NPHENERG,2)

      INTEGER i, idmax, ii, index, j, k

      data thresh /0.1/

c set convergence threshold for photon distribution

      if ((diff .lt. 0.3d0).and.(niter.gt.5)) then 
         thresh = 0.1d0
      else
         thresh = 2.d0
      endif
c      print *,'thresh',thresh,diff

c clear pout and set pin to the injected photon distribution defined in setup

      DO i = 1,NPHENERG
         pout(i) = 0.
         pin(i) = phinjfn(i)
      ENDDO

c assume e- dist fixed

      IF (.NOT.isthermonly ) THEN

         do i=2,NLORGAM
            tau = edist(i) + posdist(i)
            do j=1,NPHENERG
               if (compr(i,j).gt.0) then
                  pout(j) = pout(j) + tau*compr(i,j)
               endif
            enddo
         enddo

      ENDIF

c call annihilation + brem + anything else that makes photons

      do i=1,NPHENERG
         gder(i) = 0.d0
      enddo
      call dotherm(isthermonly)
      IF ( .NOT.isthermonly ) call dobrem     

      do i=1,NPHENERG
         pout(i) = pout(i) + tesc(i)
         pin(i) = pin(i) + gder(i)
         pin0(i) = pin(i)
      enddo

      IF ( .NOT.isthermonly ) THEN

c get nonthermal scattering matrix
c throw out row one
c can optimize this out if you want later

         do i=1,NPHENERG
            do j=1,NPHENERG
               ntscat(i,j) = 0.
            enddo
         enddo

         index = 0
         i=1
         do j=1,NPHENERG
            do k=ictab(i,j,1),ictab(i,j,2)
               index=index+1
            enddo
         enddo
 
         do i=2,NLORGAM
            tau = posdist(i) + edist(i)
            do j=1,NPHENERG
               rate = tau
               do k=ictab(i,j,1),ictab(i,j,2)
                  index=index+1
                  ntscat(k,j) =ntscat(k,j)+ctab(index)*rate
               enddo
            enddo
         enddo

         do i=1,NPHENERG
            pout(i) = pout(i) - ntscat(i,i)
            ntscat(i,i) = 0
         enddo

      ENDIF

      if ( dothermpairs ) then
         call mkscat(pout, scat, iscat)
      endif


      do i=1,NPHENERG
         pout0(i) = pout(i)
      enddo

c  now iterate over pin and pout till changes in gdist = pin/pout 
c  satisfy the convergence criterion

      dmax = thresh + 1.0

      DO WHILE ( dmax .GT. thresh ) 

         do i=1,NPHENERG
            pout(i) = pout0(i)
            pin(i) = pin0(i)
         enddo

         IF ( .NOT.isthermonly ) THEN

            do 10 i = 141,NPHENERG
               do 20 j = 1,i
                  ii = i-140
                  if(ii.gt.NPHENERG) CALL xwrite('hi',10)
                  if(ii.lt.1) CALL xwrite('lo',10)
                  if (j.gt.NPHENERG) CALL xwrite('hj',10)
                  if(j.lt.1) CALL xwrite('lj',10)
                  if (ppr(ii,j).lt. .005) goto 20
                  if (ipp(ii,j) .lt. 1) goto 20
                  erate = gdist(i)*gdist(j)*ppr(ii,j)
                  pout(i) = pout(i) + ppr(ii,j)*gdist(j)
                  pout(j) = pout(j) + ppr(ii,j)*gdist(i)
 20            continue
 10         continue


            do i=1,NPHENERG
               do j=1,NPHENERG
                  pin(i) = pin(i) + gdist(j)*ntscat(i,j)
               enddo
            enddo

         ENDIF

         if ( dothermpairs ) then
            do i=1,NPHENERG-20
               do j=iscat(i,1),iscat(i,2)
                  pin(j) = pin(j) + scat(i,j)*gdist(i)
               enddo
            enddo
         endif


         dmax = 0
         do i=1,NPHENERG
            temp = pin(i)/pout(i)
            odmax = dmax
            dmax = max(abs(temp-gdist(i))/(gdist(i)+1.d-40),dmax)
            if (dmax.ne.odmax) idmax = i
            gdist(i) = temp
         enddo

      ENDDO

      return
      end

c =========================================================
c dump out particle distributions

      SUBROUTINE info
      IMPLICIT NONE
      INCLUDE "common.eq"

      DOUBLE PRECISION enesc, etot, trnpart, rnpart, sum, xdiff

      INTEGER i, noiter
c      INTEGER nz

c      backspace(1)
c      nz = nz + 1
      write(1,*)'------------------------------------------------'
      write(1,*)'''ITERATION# ',iter,' le =',rle,' ls =',rls,''''
      write(1,*) '''rlth =',rlth,'taup=',taup,''''
      write(1,*)'''-- iterations since last report: ',iter-niter,''''
      write(1,*) '''Source Rad(m.) '''
      write(1,*) radius
c      write(1,*) '''B0,  Source Rad(m.) '''
c      write(1,*) b0,radius
      noiter = iter
      etot = 0.
      trnpart = 0.
      do 30 i = 2,NLORGAM
         etot = etot + een(i)*(edist(i)+posdist(i))
         trnpart = trnpart+(edist(i)+posdist(i))
 30   continue
      etot = etot+ 2.*eth
      trnpart = trnpart + (elecdens+posdens)
      enesc = 0.
      do 31 i = 1,NPHENERG
         etot = etot + wgam(i)*gdist(i)
         trnpart = trnpart+gdist(i)
         enesc = GEOZ * esc(i)*wgam(i) + enesc
c     print *,i,esc(i),gdist(i)
 31   continue
      write(1,*) '''theta = '''
      write(1,*) theta
      xdiff = abs(trnpart - rnpart)
      rnpart = trnpart
      write(6,*) 'Etotal = ',etot
      write(6,*) 'rnpart = ',rnpart
      write(6,*) 'idump = ',idump
      write(6,*) ' dt = ',dt
      write(1,*)'''Etotal = ''',etot
      write(1,*)'''Total # of particles:''',rnpart
      write(1,*)
      write(1,*) '''ne- and ne+ ''', elecdens, posdens
      do 10 i = 1,NLORGAM
c      write(1,*)een(i),edist(i),gnder(i),gcder(i)
      write(1,299)een(i),edist(i),posdist(i)
c299   format(2x,g12.6,4x,g13.7,4x,g13.7)
299   format(2x,g12.6,4x,g13.7,4x,g13.7)
10    continue
      write(1,*)
      write(1,*)'''photon bin values:'''
      sum = 0.
      do 20 i = 1,NPHENERG
         sum = sum + wgam(i)*esc(i)*GEOZ
         write(1,399) wgam(i),GEOZ*phinjfn(i),GEOZ*esc(i),sum
 20   continue
 399  format(2x,g12.6,4x,g13.7,4x,g13.7,4x,g13.7)
      write (1,*) '''Escaping photon luminosity: ''',enesc
      write (1,*) '''Escaping pair   luminosity: ''',.3*eth*2.*GEOZ
c      print *,nz
c      write (1,*) 'Pair yield: ',nz
c     -- For Piero
c      write(1,*) '''Tau_comp(x): '''
c      do 88 i = 1,NPHENERG
c        tau = 0.
c        do 89 j = 2,NLORGAM
c          tau = tau+2.*ecomp(j,i,1)*edist(j)
c89      continue
c        write(1,499) wgam(i),tau
c499     format(2x,g12.6,4x,g13.7)
c88    continue
c      endfile(1)
      call flush(1)
c        do i=1,NPHENERG
c          gder(i) = 0.
c        enddo
c        call dotherm(isthermonly)
c        
c        do i=1,NPHENERG
c        write(2,*) wgam(i),GEOZ*gder(i)*tesc(i)
c        enddo
      return
      end

c==============================================================
c Finds electron temperature
c Global variables:
c   i: tesc, otesc, taup
c   i/r: gdist, theta
c   r: elecdens, posdens

c Calls routines getethd, dotherm, docoul, dobrem

      SUBROUTINE getth(isthermonly, ein1, posin1, psum)
      IMPLICIT NONE
      INCLUDE "common.eq" 

      LOGICAL isthermonly
      DOUBLE PRECISION ein1, posin1, psum

c Local variables

      DOUBLE PRECISION coulethd, ctemp, dn, dnth, drdth, dth, ethder0
      DOUBLE PRECISION oethder, otheta, realotheta, sarg, temp
      DOUBLE PRECISION eout2

      INTEGER i, ith, ithmax

      LOGICAL iboink

c suppressing an unnecessary compiler warning
      coulethd = 0.d0

c  ------   Get Optical Depth in equilibrium
c for pair inject
      if (qpairinj) then
         ein1 = psum + erninj
         posin1 = psum + erninj
      else
         ein1 = psum
         posin1 = psum
      endif

c new experiment -- let's compute tau dynamically again and see what 
c happens
c      if (iread.gt.0) then
      call getethd(isthermonly, eout2)
      IF ( isthermonly ) THEN
         elecdens = taup
         posdens = 0.d0
      ELSE
         sarg = sqrt(taup*taup+4.d0*ein1/eout2)/2.d0
         elecdens = sarg+0.5d0*taup
         posdens = sarg-0.5d0*taup
      ENDIF
      
      call dotherm(isthermonly)
 
   
      do i=1,NPHENERG
         gdist(i) = gdist(i)*otesc(i)/tesc(i)
c     muy important 
         otesc(i)  = tesc(i)
      enddo
c      endif


c if we are not calculating thermal particle distribution then return

      if ( .NOT.dothermpairs ) RETURN

c otherwise find the temperature of the thermal particles

      realotheta = theta

      dnth = 2.d0*psum
c pair injection
      if (qpairinj) then 
         ethder0 = dnth*een(1) + (rlth/GEOZ) + 2.d0*erninj*een(1)
      else
c electron injection
         ethder0 = dnth*een(1) + (rlth/GEOZ)
     &        +(een(1)-1.d0)*(erninj+prninj)
      endif
         
c try a trial temperature and figure out which the gradient
c of ethder is pointing
      temp = theta
      if (ethder .gt. 0.) then
         theta = theta*1.001
      else
         theta = theta*.999
      endif
      ethder = ethder0
      call getethd(isthermonly, eout2)
c ** don't forget coulomb heating & brems. note that coulomb heating
c    term is saved for use in loop.
      ctemp = ethder
      IF ( .NOT.isthermonly ) THEN
         call docoul
         coulethd = ethder-ctemp
         call dobrem
      ENDIF
      oethder = ethder
      otheta = theta
      theta = temp
      iboink = .false.
      dn = 1.0 

cdebug
cc      WRITE(*,'(a,4(1x,1pg13.6))') 'getth initial: ', theta, ethder, 
cc     &   elecdens, posdens
cdebug

c loop on temperature search

      IF ( isthermonly ) THEN
cdebug         ithmax = 5000
         ithmax = 15
      ELSE
         ithmax = 15
      ENDIF
      ith= 0
      dth = 1.0

      DO WHILE ( ith .LE. ithmax .AND. ABS(dth) .GT. 1.0d-6 )

         ethder = ethder0
         call getethd(isthermonly, eout2)
         IF ( .NOT.isthermonly ) THEN
c ** don't forget coulomb heating & brems
c           call docoul
c BUT don't recompute coul heating during temperature search
c too unstable! Just add in saved coulomb heating term.
            ethder=ethder + coulethd
            call dobrem
         ENDIF
         if (ethder.eq.oethder) then
            print *,'temp ethder mucked!'
         endif
         if (theta.ne.otheta) then
            drdth=(ethder-oethder)/(theta - otheta)
         else
            drdth = 0.0
            print *, 'theta = otheta = ', theta, 
     &              '  this should not happen'
         endif
ccdebug
cc         write(*,'(a,i5,1x,3(1pg13.6,2x))') 
cc     &         'ith,theta,ethder,drdth : ',ith,theta,ethder,drdth
ccdebug
         otheta = theta
         oethder = ethder
         if(drdth.ne.0) theta = theta - ethder/(drdth)*dn
         if (theta.lt.1.d-4) then
            dn = (otheta-1.d-4)/ethder*drdth
            theta = otheta
            theta = 0.5*otheta
c     print *,'theta boink!... ouch'
            iboink  = .true.
         endif
         if (theta.gt.50.d0) then
            theta = min(50.d0,1.07*otheta)
c     print *,'theta bink!...eech'
         endif

         dth = abs(theta - otheta)/(abs(theta)+1.d-21)
c     print *,'th,oth,thd,othd',theta,otheta,ethder,oethder
         ith = ith + 1

      ENDDO

      if (ith.gt.ithmax) then
         WRITE(*,*) 'ith = ', ith, ' ZZZZZZZZZZZZZZZZZZZZZZZZZZZZAP'
         WRITE(*,'(a,4(1x,1pg13.6))') 'th,oth,thd,othd',theta,otheta,
     &      ethder,oethder
         if (iboink) then
            theta = max(0.8*theta,1.d-4)
         else
            theta = realotheta*1.2
         endif
      endif

      dn = min(.5d0,1.d0/(1.d0+3.d0*realotheta))
      theta = dn*theta + (1.d0-dn)*realotheta
   
c       print *,'final theta = ',theta

      return
      end


c=====================

      SUBROUTINE  load
      IMPLICIT NONE
      INCLUDE "common.eq"

      INTEGER ilun
      INTEGER lenact, blocksize, rwmode, status

      CHARACTER(512) ctmp
      character(255) FileName, datadir, fgmodf
      EXTERNAL fgmodf

c get free unit number from xspec
      call getlun(ilun)

      CALL xwrite('Reading data files...',15)

      rwmode = 0
      blocksize = 1

      datadir = fgmodf()

      CALL xwrite('ltcz',15)
      FileName = datadir(:lenact(datadir))//'ltcz.fits'
      status = 0
      call ftopen(ilun, FileName, rwmode,  blocksize, status)

      IF ( status .NE. 0 ) THEN
         ctmp = 'Failed to open '//Filename(:lenact(Filename))
         CALL xwrite(ctmp,10)
         RETURN
      ENDIF
c clav
      call fitsread(ilun,1,clav,2,0,61,41,status)
c ig
      call ifitsread(ilun,2,ig,3,61,41,2,status)
c rldisp
      call fitsread(ilun,3,rldisp,3,61,41,2,status)
      call ftclos(ilun,status)

      CALL xwrite('comp',15)
      FileName = datadir(:lenact(datadir))//'comp.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c ecomp
      call fitsread(ilun,1,ecomp,3,61,161,3,status)
c gcomp
      call fitsread(ilun,2,gcomp,3,61,161,3,status)
c iecomp
      call ifitsread(ilun,3,iecomp,3,61,161,3,status)
c igcomp
      call ifitsread(ilun,4,igcomp,3,61,161,3,status)
      call ftclos(ilun,status)

      CALL xwrite('azntc',15)
      FileName = datadir(:lenact(datadir))//'azntc.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c av(NLORGAM,NPHENERG)
      call fitsread(ilun,1,av,2,0,NLORGAM,NPHENERG,status)
c compr(NLORGAM,NPHENERG)
      call fitsread(ilun,2,compr,2,0,NLORGAM,NPHENERG,status)
      call ftclos(ilun,status)


      CALL xwrite('clt2',15)
      FileName = datadir(:lenact(datadir))//'clt2.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c cltr(141,241)
      call fitsread(ilun,1,cltr,2,0,141,241,status)
c clta(141,241)
      call fitsread(ilun,2,clta,2,0,141,241,status)
c cltd(141,241)
      call fitsread(ilun,3,cltd,2,0,141,241,status)
      call ftclos(ilun,status)


      CALL xwrite('azcomp',15)
      FileName = datadir(:lenact(datadir))//'azcomp.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c compr(NLORGAM,NPHENERG)
      call fitsread(ilun,1,compr,2,0,NLORGAM,NPHENERG,status)
c ictab(NLORGAM,NPHENERG,2)
      call ifitsread(ilun,3,ictab,3,NLORGAM,NPHENERG,2,status)
c ctab(520483)
      call fitsread(ilun,5,ctab,1,0,0,520483,status)
      call ftclos(ilun,status)

      CALL xwrite('azecomp',15)
      FileName = datadir(:lenact(datadir))//'azecomp.fits'
      call ftopen(ilun, filename, rwmode,  blocksize, status)
c iectab(NLORGAM,NPHENERG)
      call ifitsread(ilun,1,iectab,3,NLORGAM,NPHENERG,2,status)
      call fitsread(ilun,1,clav,2,0,61,41,status)
c ectab(218997)
      call fitsread(ilun,2,ectab,1,0,0,218997,status)
      call ftclos(ilun,status)


      CALL xwrite('azpp',15)
      FileName = datadir(:lenact(datadir))//'azpp.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c ppr(NLORGAM,NPHENERG)
      call fitsread(ilun,1,ppr,2,0,NLORGAM,NPHENERG,status)
c ipp(NLORGAM,NPHENERG)
      call ifitsread(ilun,2,ipp,2,0,NLORGAM,NPHENERG,status)
c ppf(NLORGAM,NPHENERG)
      call fitsread(ilun,3,ppf,2,0,NLORGAM,NPHENERG,status)
      call ftclos(ilun,status)


      CALL xwrite('thrm',15)
      FileName = datadir(:lenact(datadir))//'thrm.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c thann(3,2)
      call fitsread(ilun,1,thann,2,0,3,2,status)
c cr(162,2)
      call fitsread(ilun,3,cr,2,0,162,2,status)
      call ftclos(ilun,status)


      CALL xwrite('htc',15)
      FileName = datadir(:lenact(datadir))//'htc.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c htcr(161,4)
      call fitsread(ilun,1,htcr,2,0,161,4,status)
c htcav(161,4)
      call fitsread(ilun,2,htcav,2,0,161,4,status)
c htcd(161,4)
      call fitsread(ilun,3,htcd,2,0,161,4,status)
c htcavl(NLORGAM,4)
      call fitsread(ilun,4,htcavl,2,0,NLORGAM,4,status)
c htcdl(NLORGAM,4)
      call fitsread(ilun,5,htcdl,2,0,NLORGAM,4,status)
      call ftclos(ilun,status)

      CALL xwrite('brem',15)
      FileName = datadir(:lenact(datadir))//'brem.fits'
      call ftopen(ilun, FileName, rwmode,  blocksize, status)
c brem(61)
      call fitsread(ilun,1,brem,1,0,0,61,status)
c brema(7710)
      call fitsread(ilun,2,brema,1,0,0,7710,status)
c ibrem(61)
      call ifitsread(ilun,3,ibrem,1,0,0,61,status)
      call ftclos(ilun,status)

      CALL xwrite('Done.',15)

      call frelun(ilun)
      return
      end


      SUBROUTINE ifitsread(lun,hdu,dummy,naxis,nax3,nax2,nax1,status)
      IMPLICIT NONE

      INTEGER dummy(*), bitpix, naxis, naxes(4), nelements, group
      INTEGER firstpixel, lun , status, nax3, nax2, nax1, nullval
      INTEGER hdutype, hdu

      LOGICAL anyf

      call ftmahd(lun, hdu, hdutype, status)

      naxes(3) = nax3
      naxes(2) = nax2
      naxes(1) = nax1
      bitpix = 32
      firstpixel=1
      group=1
      status = 0
      nullval = 0
      nelements = naxes(1)
      if (naxis .gt. 1) nelements=nelements*naxes(2)
      if (naxis .gt. 2) nelements=nelements*naxes(3)
c     print *,'nel i',nelements, nelements*4

      call ftgpvj(lun, group, firstpixel, nelements, nullval, 
     &     dummy, anyf, status)
      if (status .ne.0) CALL xwrite('failed ifitsread',10)
      return
      end

      SUBROUTINE fitsread(lun,hdu,dummy,naxis,nax3,nax2,nax1,status)
      IMPLICIT NONE

      DOUBLE PRECISION dummy(*), nullval

      INTEGER bitpix, naxis, naxes(4), nelements, group, hdu
      INTEGER firstpixel, lun , status, nax3, nax2, nax1, hdutype

      LOGICAl anyf

      call ftmahd(lun, hdu, hdutype, status)

      naxes(3) = nax3
      naxes(2) = nax2
      naxes(1) = nax1
      bitpix = -64
      firstpixel=1
      group=1
      status = 0
      nullval = 0.d0
      nelements = naxes(1)
      if (naxis .gt. 1) nelements=nelements*naxes(2)
      if (naxis .gt. 2) nelements=nelements*naxes(3)
c     print *,'nel d',nelements,nelements*8

      call ftgpvd(lun, group, firstpixel, nelements, nullval, 
     &     dummy, anyf, status)

      if (status .ne.0) CALL xwrite('failed fitsread',10)
      return
      end

c=================================
c not currently used
c
c      SUBROUTINE  pprod
c      IMPLICIT NONE
c      INCLUDE "common.eq"
c
c      DOUBLE PRECISION erate
c      INTEGER i, j, ii
c
c      psum = 0.
c      do 10 i = 141,NPHENERG
c         do 20 j = 1,i
c            ii = i-140
c            if (ppr(ii,j).lt. .005) goto 20
c            if (ipp(ii,j) .lt. 1) goto 20
c            erate = gdist(i)*gdist(j)*ppr(ii,j)
c            psum = psum + erate
cc        if ((erate .lt. rmin) .or. (erate .lt. 1.e-9)) goto 99
c            pout(i) = pout(i) + ppr(ii,j)*gdist(j)
c            pout(j) = pout(j) + ppr(ii,j)*gdist(i)
c 20      continue
cc       if (i.gt.181) print *,i,pout(i)
c 10   continue
c      return
c      end

c=================================
      SUBROUTINE setup(param,np,imodel,radnorm)
      IMPLICIT NONE
      INCLUDE "common.eq"

      INTEGER np, imodel
      DOUBLE PRECISION radnorm
      REAL param(np)

      DOUBLE PRECISION ewidth(NLORGAM)

      DOUBLE PRECISION b0, Rin, rntot, tbb, temp
      DOUBLE PRECISION texp, ue, urad
      DOUBLE PRECISION upeenl, upeenh, plawin

      REAL tparam(2),tphotar(NPHENERG),tear(0:NPHENERG)
      REAL tphoter(NPHENERG)

      INTEGER i, igmin, igmax

      CHARACTER(255) comment

      LOGICAL qfirst
      SAVE qfirst

      DATA qfirst /.TRUE./


c      print *,'in setup' 
      IF ( imodel .EQ. 1 ) THEN
         rle = param(1)*param(2)*param(4)
         rls = param(2)
         rlth = param(1)*param(2)*(1-param(4))
      ELSE IF ( imodel .EQ. 2 ) THEN
         rle = 0.0d0
         rls = param(2)
         rlth = param(1)*param(2)*(1-param(4))
      ELSE IF ( imodel .EQ. 3 ) THEN
         rle = 0.0d0
         rls = 1.0d0
c electron temperature in m_ec^2
         theta = param(1)/5.11d2
         rlth = param(4)
      ENDIF
c blackbody temperature in m_ec^2
      bbtemp = param(3)/5.11d5
      taup = param(5)
      radius  = param(6)

      Rin = param(18)

      igmax = max(min(int(log10(param(8))/.05d0+1.5d0),NLORGAM),2)
      plawin = param(9)
      igmin = 1
      if (plawin .gt. 0.d0) then
         igmin = max(min(int(log10(param(7))/.05d0+1.5d0),NLORGAM),2)
      endif
      if (param(10).gt.0) then
         qpairinj = .true.
      else
         qpairinj = .false.
      endif
      b0 = 0.d0

      WRITE(comment, '(3(a,1pg13.6))') 'l_bb = ', rls, 
     &   '   l_nth = ', rle, '   l_th = ', rlth
      CALL xwrite(comment, 20)

c      print *,'db: rle,rls,theta,radius',rle,rls,theta,radius
c      print *,'igmin,igmax,plawin',igmin,igmax,plawin
c      print *,'MAXITER',MAXITER

c initial idump = bin# e- cascade into when using stationary approx
      idump = NLORGAM
c initialize energy bin values
      een(1) = (10.**.05d0 - 1.)/2./(10.**.025d0 - 1.)
      upeenl = 10.**.025d0
      ewidth(1) = upeenl - 1.d0
      do 58 i = 2,NLORGAM
         temp = 10.**(float(i-1)*.05d0)
         een(i) = .5d0*temp*(10.**.025d0 + 10.**(-.025d0))
         upeenh = 10.**.025d0*temp
         ewidth(i) = upeenh - upeenl
         upeenl = upeenh
 58   continue

c initialize energies and bin widths for photon arrays
      do 62 i = 1,NPHENERG
         temp = 10.**(float(i-101-IOFF)*.05d0)
         wgam(i) = .5d0*temp*(10.**.025d0 + 10.**(-.025d0))
         upwgam(i)= 10.**.025d0*temp
c        print *,i,wgam(i)
c        pause
 62   continue
      do 60 i = 2,NPHENERG
         gwidth(i) = upwgam(i) - upwgam(i-1)
 60   continue
      gwidth(1) = upwgam(1)-10.**(-7.025d0)
c
c  --- get model parameters
c      write(6,*) 'Hello...'
c      write(6,*) 'Model parameters are: '
c      write(6,*) 'bbtemp = ',bbtemp
c define source (injection) functions
c pairinjfn(i) - PAIR inject funct, phinjfn(i) - photon inject. func.
c power law(?) injection
      do 45 i = 1,NLORGAM
         pairinjfn(i) = 0.
         epinjfn(i) = 0.
 45   continue

c Take care of injection:
c ***

      urad= 0.0
      ue = 0.0
c     *factor of -1 is to agree with big Z *
      if (rle .ne. 0.) then
         if (plawin .lt. 0.) then
            write(6,*) '** Monoergetic injection at gamma =',een(igmax)
            pairinjfn(igmax) = 1.
         else
c          write(6,*) '** Power law injection with index: ',plawin
c          write(6,*) '   from gmin =',een(igmin),'to gmax =',een(igmax)
            do 422 i = igmin,igmax
               pairinjfn(i) = een(i)**(-plawin)*ewidth(i)
 422        continue
         endif
         rntot = 0.
         do 47 i = 1,NLORGAM
            rntot = rntot + pairinjfn(i)
            ue = pairinjfn(i)*(een(i)) + ue
 47      continue
c remove rest mass since it is reaccelerated away
c if have pair injection then things more complicated because
c of pair annihilation -- in practice we won't bother with it
         if (qpairinj) then 
c           print *,'PAIR INJECTION'
            ue = 2*ue
         else
            ue = ue - rntot
c           print *,'ELECTRON INJECTION'
         endif


         do 67 i = 1,NLORGAM
            pairinjfn(i) = pairinjfn(i)/ue*rle / GEOZ
 67      continue
         ue = 0.
         do 68 i = 1,NLORGAM
            ue = 2.*pairinjfn(i)*(een(i))*GEOZ + ue
 68      continue
      endif
c  black body photon injection
      if (rls .ne. 0.) then
         if (bbtemp .lt. 0.) then 
c         use diskpn model
            tparam(1) = SNGL(-bbtemp*511.0)
            tparam(2) = SNGL(Rin)
            do i=1,NPHENERG
               tear(i) = SNGL(upwgam(i)*511.0)
            enddo
            tear(0) = SNGL(upwgam(1)*10.**(-0.05d0)*511.0)
            call xsdiskpn(tear,NPHENERG,tparam,1,tphotar,tphoter)
            do i=1,NPHENERG
               phinjfn(i) = tphotar(i)
            enddo
c        ^^ may be a constant of 511 floating around, but gets normalized out
c           I think
         else 
c         use standard blackbody
            do 66 i = 1,NPHENERG
               tbb = 1.e18*gwidth(i)*wgam(i)*wgam(i)
               texp = wgam(i)/bbtemp
               texp = .5*(texp-1.+sign(texp+1.,60.-texp))
               tbb = tbb/(exp(texp)-1.)
               phinjfn(i) = dmax1(0.d0,tbb)
 66         continue
         endif
         do 48 i = 1,NPHENERG
            urad = phinjfn(i)*wgam(i)+ urad
 48      continue
         do 49 i = 1,NPHENERG
            phinjfn(i) = phinjfn(i)/urad*rls/GEOZ
 49      continue

c Normalization of disk component, will be used to renormalize total spectrum
c to make the XSPEC norm of the total spectrum equal norm of the disc spectrum.
c Thus, the inner disk radius can be computed.

         radnorm = urad * GEOZ / rls

         urad = 0.0
         do 50 i = 1,NPHENERG
            urad = phinjfn(i)*wgam(i)*GEOZ+ urad
 50      continue
      else
         do 102 i = 1,NPHENERG
            phinjfn(i) = 0.
 102     continue
      endif

c      write(6,*) 'urad-dot: ',urad
c      write(6,*) 'ue-dot: ',ue
c      write(6,*) 'pairinjfn(61) = ',pairinjfn(61)

      IF ( imodel .NE. 3 ) THEN

c      if don't read from file, put
c        initial non-zero starting distribution into bins
c !! kaa  added posdens and posdist initialization
         elecdens = 1.0e-7
         posdens = 1.0e-7
         do 87 i = 1,NLORGAM
            edist(i) = 1.e-7
            posdist(i) = 1.e-7
 87      continue
         eth = elecdens*een(1)
c fudge --
         theta = 2./3.*(eth/elecdens-1.)

      ENDIF

c  initialize tables
      if ( qfirst ) then 
         call load
         qfirst = .FALSE.
c     write(6,*) 'file tables loaded in...?'
      endif


c  ** Initial guesses for various parameters
c!! kaa  added posdens and posdist initialization
      elecdens = 0.0
      posdens = 0.0
      do 766 i = 1,NLORGAM
         edist(i) = 0.
         posdist(i) = 0.
 766  continue

      IF ( imodel .EQ. 3 ) THEN
         elecdens = taup  
         posdens = 0.d0
      ENDIF

      do i=1,NPHENERG
         gdist(i) = phinjfn(i)
      enddo

c      write(6,*) 'o.k. here goes nothing...'
      return
      end

c **************************************************************************
c Main iteration routine to derive photon and electron distributions. 
c Global variables:
c   i: pairinjfn, epinjfn, phinjfn, wgam, tesc, gdist, taup. 
c   r: esc, otesc, elecdens, posdens.

c Calls routines doelec, dopcomp

      SUBROUTINE  solve(isthermonly)
      IMPLICIT NONE

      LOGICAL isthermonly

      INCLUDE  "common.eq"

c Local variables. ECONV is convergence criterion for energy balancing
c (escape = inject). MAXITER is maximum number of iterations to try.

      DOUBLE PRECISION ECONV
      PARAMETER(ECONV=2.0d-2)
      INTEGER MAXITER
      PARAMETER(MAXITER=10000)

      DOUBLE PRECISION ogdist(NPHENERG)
      DOUBLE PRECISION tenesc
      DOUBLE PRECISION ediff, enesc
      DOUBLE PRECISION tdiff1, temp

      INTEGER i, id

      CHARACTER(255) comment

c total energy escaping must equal that being injected.

      tenesc = (rle+rls+rlth)/GEOZ

      dothermpairs = .FALSE.

      elecdens = 0.d0 + taup
      IF ( isthermonly ) posdens = 0.d0
c      print *,'>> initial guess for tau:',elecdens

c Total energy in non-thermal particles

      erninj = 0.
      prninj = 0.

      IF ( .NOT. isthermonly ) THEN
         DO i=2,NLORGAM
            erninj = erninj + pairinjfn(i)
            prninj = prninj + epinjfn(i)
         ENDDO
         erninj = max(erninj - 1.d-20,0.d0)
c     print *,'erninj',erninj
         prninj = max(prninj -1.d-20,0.d0)
      ENDIF

      idump = 1

c iteration difference to stop at before start thermal
      tdiff1 = 1.d0

c  --- start main loop
c  --- (assume theta, eth, elecdens initialized in setup! Ja?)

      niter = 0
      do i=1,NPHENERG
         otesc(i) = 1.
         ogdist(i) = 0.
      enddo

c Loop is over photon distribution. Convergence is when the largest change in any
c bin of the photon distribution is less than 1% and the total energy escaping is
c within 2% of that expected.

      DO WHILE ( niter .LE. MAXITER )

c get electron distribution for the current photon distribution.

         call doelec(isthermonly)

c now do photon distribution assuming the electron distribution just calculated.

         call dopcomp(isthermonly)

c set the distribution of escaping photons

         do i=1,NPHENERG
            esc(i) = tesc(i)*gdist(i)
            otesc(i) = tesc(i)
         enddo

c ---- Overall convergence check (look at photon distribution)
         diff = 0.
         temp = 0.0
         DO i = 1,NPHENERG
c modification from Joern Wilms
            if (gdist(i).gt.1.d-38) then
               temp=abs(ogdist(i)-gdist(i))/gdist(i)
               if (temp.gt.diff) then 
                  id = i
                  diff = temp
               endif
            endif
         ENDDO

c don't forget to update the array which stores the last calculation of
c the photon distribution

         DO i =1,NPHENERG
            ogdist(i) = gdist(i)
         ENDDO

         enesc = 0.0
         DO i=1,NPHENERG
            enesc = enesc + esc(i)*wgam(i)
         ENDDO
c      print *,'energy escaping',enesc,'vs. expected',tenesc

cdebug
c         WRITE(*,'(i4,1x,1pg14.7,1x,i6,3(1x,1pg14.7))') niter, theta, 
c     &          id, gdist(100), edist(50), enesc
cdebug
         
      

         if ( diff .lt. tdiff1 ) then
            if (.NOT.dothermpairs) then
               dothermpairs = .TRUE.
               tdiff1=1.d-2
            else
               ediff = abs(enesc-tenesc)/(tenesc+1.d-12)
               if( ediff .lt. ECONV ) then

                  IF ( isthermonly ) THEN
                     WRITE(comment, '(2(a,1pg14.7))') 
     &                 '***eqtherm: theta = ', theta*511.d0,
     &                 ' keV,    tau_T = ',elecdens+posdens
                     CALL PDBVAL("eqtherm_theta", theta*511.d0)
                     CALL PDBVAL("eqtherm_tau_T", elecdens*posdens)
                  ELSE
                     WRITE(comment, '(2(a,1pg14.7))') 
     &                 '***eqpair: theta = ', theta*511.d0,
     &                 ' keV,    tau_T = ',elecdens+posdens
                     CALL PDBVAL("eqpair_theta", theta*511.d0)
                     CALL PDBVAL("eqpair_tau_T", elecdens*posdens)
                  ENDIF
                  CALL xwrite(comment, 15)

                  return
               endif
            endif
         endif

         niter = niter + 1

c -- go back for another round of convergence talks
      ENDDO

      IF (niter .gt. MAXITER) THEN
         CALL xwrite(' >> Run ended prematurely: MAXITER hit', 10)
         RETURN
      ENDIF

      end 


c ================================================================
c Calculate brems spectrum from particle distribution
c Global variables:
c   i: wgam, gwidth, elecdens, posdens, edist, posdist
c   i/r: gdot

      SUBROUTINE dobrem
      IMPLICIT NONE
      INCLUDE "common.eq"

      DOUBLE PRECISION coeff, delta, e1, eloss, eugam, otheta, pee
      DOUBLE PRECISION rba, rbb, rbc, rbd, rga, rgb, rgc, rgd
      DOUBLE PRECISION rndot, rnprod, temp, ttheta, v1, v2, x, x1
      DOUBLE PRECISION x2

      INTEGER i, imax, imin, index, iswitch, it, iu, j, jj, jshift

c *** calculate brems spectra. Careful with factors of two.!
c (ethder is TOTAL energy now)
c     common /bremd/ spec(241)
c coeffs for Stepney fit:
      DOUBLE PRECISION ba(13),bb(13),bc(13),bd(13),ga(13),gb(13),gc(13)
      DOUBLE PRECISION gd(13),fiten(13)
      DOUBLE PRECISION ca(11),cb(11),cc(11),cd(11),ce(11)
      DOUBLE PRECISION ha(11),hb(11),hc(11),hd(11),he(11),bfiten(11)
      data ba / 1.584,1.357,1.197,1.023,
     *          .883, .700 , .572,.484,
     *          .417, .361, .322 , .286,
     *          .259 /
      data bb / .578, .437,.291,.204,
     *          .0835,-0.0494,-.139,-.181,
     *          -.209,-.240, -.244,-.257,
     *          -.258 /
      data bc / 4.565,3.842,3.506,3.036,
     *          2.831,2.545,2.352,2.175,
     *          2.028,1.914,1.795,1.705,
     *          1.617 /
      data bd / 2.091,1.855,1.672,1.593,
     *          1.487,1.364,1.254,1.179,
     *          1.108,1.030,0.982,0.923,
     *          0.879 /
      data ga / .0387,.0633,.0862,.128,
     *          .159, .208, .234, .245,
     *          .248, .247, .243, .239,
     *          .235 /
      data gb / .523,.540,.569,.596,
     *          .658,.633,.643,.695,
     *          .729,.756,.763,.755,
     *          .735 /
      data gc / 5.319,4.412,3.897,3.383,
     *          2.974,2.738,2.424,2.025,
     *          1.716,1.457,1.271,1.140,
     *          1.060 /
      data gd / .782,.689,.633,.523,
     *          .532,.326,.302,.394,
     *          .453,.500,.515,.508,
     *          .478 /
      data fiten / 0.0978,.1468,.1957,.2935,
     *             .39139,.5871,.7828,.9785,
     *             1.1742,1.3699,1.5656,1.7613,
     *             1.956947 /
c  coeffs for Haug fit:
      data ca / 3808.,657.,420.3,130.0,67.55,
     *          40.72, 20.28,12.91,10.52,9.61,9.22 /
      data cb / -2560.,-4427.,250.9,-0.8,-18.79,
     *          -28.02, -29.46,-22.30,-12.14,-3.04,4.52 /
      data cc / 1036.,-114.,120.8,28.1,12.61,
     *          3.26,-1.27, -3.56, -2.97, -2.52, -2.16 /
      data cd / 7237.,7209.,0.,134.0,111.33,
     *          107.21, 91.86, 75.21, 56.00, 39.96, 27.09 /
      data ce / -539.,-4927.,566.9,91.9,36.47,
     *          8.59,-4.59,-1.88,8.54,18.28,26.31 /
      data ha / 17.619,5.124,2.633,2.854,3.315,
     *          3.629,3.991,4.262,4.366,4.383,4.374 /
      data hb / -388.0,-98.27,-10.22,8.583,11.86,
     *          13.33, 14.56, 13.26, 12.56, 12.07, 11.29 /
      data hc / 4015.7,1132.4,242.9,99.08,76.26,
     *          63.96, 47.92, 44.57, 39.36, 34.58, 32.34 /
      data hd / 3587.6,1649.4,454.4,155.6,70.38,
     *          42.07, 26.37, 12.06, 9.05, 8.59, 5.86 /
      data he / -535.0,-406.0,-122.8,-40.24,-14.01,
     *          -7.19, -5.57, -0.82,-0.80,-1.38,-0.46 /
      data bfiten / 0.01957,.03913,0.0978,.1957,.2935,
     *              .39139,.5871,.7828,.9875,1.1742,1.3699 /

      otheta = theta
      theta = max(1.d-4,theta)

      if ( dothermpairs ) then 

         if (theta .le. bfiten(11)) then
c  use Haug/Stepney approx

c take care of  thermal bremsstrahlung.
c  Use Stepney for e-e-/e+e+  - typically accurate to w/in 10% over all range.
            if (theta .le. fiten(1)) then
c   use fudged up Haug 1975b non-rel approx; not so great, but o.k.
               iswitch = min(int(log10(theta*1.05)/.05+RIOFF1),NPHENERG)
               imax = min(int(log10(theta*10.)/.05+RIOFF1),NPHENERG)
               coeff = (elecdens*elecdens+posdens*posdens)
     &              /137.*sqrt(fiten(1))/sqrt(theta)
               do 161 i = 1,iswitch
                  x = wgam(i)/theta
                  rndot= gwidth(i)*coeff/x*exp(-x)*((-ba(1)-bb(1)*x)*
     &                 log(x)+bc(1)+bd(1)*x)
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 161           continue
               do 171 i = iswitch+1,imax
                  x = wgam(i)/theta
                  rndot = gwidth(i)*coeff/x*exp(-x)*(ga(1)*x*x+gb(1)*x+
     &                 gc(1)+gd(1)/x)
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 171           continue
            else
               ttheta = theta
               if (theta .gt. fiten(13)) then
                  CALL xwrite('TH EE BREM: THETA TOO HIGH !!',10)
                  ttheta = fiten(13)
               endif
               it = 1
c bozito!! the energies aren't evenly spaced!!
               do 96 i = 1,13
                  if (ttheta .ge. fiten(i)) it = i
 96            continue
               it = min(it,12)
               delta = (ttheta-fiten(it))/(fiten(it+1)-fiten(it))
               rba = ba(it) + (ba(it+1)-ba(it))*delta
               rbb = bb(it) + (bb(it+1)-bb(it))*delta
               rbc = bc(it) + (bc(it+1)-bc(it))*delta
               rbd = bd(it) + (bd(it+1)-bd(it))*delta
               rga = ga(it) + (ga(it+1)-ga(it))*delta
               rgb = gb(it) + (gb(it+1)-gb(it))*delta
               rgc = gc(it) + (gc(it+1)-gc(it))*delta
               rgd = gd(it) + (gd(it+1)-gd(it))*delta
               imin = max(int(log10(theta*.05)/.05+RIOFF1),0)
               temp = gder(imin)
               iswitch = min(int(log10(ttheta*1.05)/.05+RIOFF1),
     &                       NPHENERG)
               imax = min(int(log10(ttheta*10.)/.05+RIOFF1),NPHENERG)
c ***** DON'T forget the FACTOR of two -- e-e- AND e+e+ brem
c   HOWEVER, ethder is one-half total energy radiatied
c   since are only keeping track of the e- distribution.
c   (Check: energy conservation o.k.)
               coeff = (elecdens*elecdens+posdens*posdens)/137.
c        do 160 i = imin+1,iswitch
               do 160 i = 1,iswitch
                  x = wgam(i)/ttheta
                  rndot= gwidth(i)*coeff/x*exp(-x)*((-rba-rbb*x)*log(x)+
     &                 rbc+rbd*x)
c        trndot=8.*3.14/3.*exp(-x)*((rba+rbb*x)*log(1./x)+rbc+rbd*x)/x
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 160           continue
               do 170 i = iswitch+1,imax
                  x = wgam(i)/ttheta
                  rndot = gwidth(i)*coeff/x*exp(-x)*
     &                 (rga*x*x+rgb*x+rgc+rgd/x)
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 170           continue
c Extend Stepney results down below .05*theta (sie werden negatif?)
c Guess low energy behavior
c        temp = gder(imin+1)-temp
c        do 172 i = 1,imin
c        rndot = temp*(1.-.28*log(wgam(i)/theta))/
c     &    (1.-.28*log(wgam(imin+1)/theta))
c         gder(i) = gder(i) + rndot
c         ethder = ethder - wgam(i)*rndot/2.0
c172   continue
            endif

c Now do e+/e- case.
c Use Haug (1987) "simple(?) analytical fits"
c Note: now there is NO factor of two in the emission coeff.
c  integral takes into account all possible combinations.
c  However, ethder must be divided by two since only half of
c  the energy radiated comes from either the pos. or e- dist.
c (but NOT in this version)
            if (theta .le. bfiten(1)) then
c   use fudge to non-rel approx. (calculating the bessel
c   function K0(x) is too costly)
               imax = min(int(log10(theta*10.)/.05+RIOFF1),NPHENERG)
               iswitch = min(int(log10(theta*1.)/.05+RIOFF1),NPHENERG)
               coeff = elecdens*posdens*2.385225686e-6/theta**1.5
               do 61 i = 1,iswitch
                  x = wgam(i)/theta
                  rndot= gwidth(i)*coeff/x*exp(-x)
                  rndot=rndot*(cd(1)+x*ce(1)-(ca(1)+cb(1)*sqrt(x)
     *                 +cc(1)*x)*log(x))
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 61            continue
               do 71 i = iswitch+1,imax
                  x = wgam(i)/theta
                  rndot = gwidth(i)*coeff/x*exp(-x)
                  rndot =rndot*(ha(1)*x*x+hb(1)*x+
     &                 hc(1)+hd(1)/x+he(1)/x/x)
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 71            continue
            else
               ttheta = theta
               if (theta .gt. bfiten(11)) then
                  CALL xwrite('TH E+E- BREM: THETA TOO HIGH !!',10)
                  ttheta = bfiten(11)
               endif
               it = 1
               do 97 i = 1,11
                  if (ttheta .ge. bfiten(i)) it = i
 97            continue
               it = min(it,10)
c ** NOTE: Unlike Stepney fits, here have to interpolate VALUES not Coeffs
               iu = it+1
               delta = (ttheta-bfiten(it))/(bfiten(iu)-bfiten(it))
               iswitch = min(int(log10(ttheta*1.0)/.05+RIOFF1),NPHENERG)
               imax = min(int(log10(ttheta*10.)/.05+RIOFF1),NPHENERG)
               coeff = elecdens*posdens*8.7128618e-4
               do 60 i = 1,iswitch
                  x1 = wgam(i)/bfiten(it)
                  x2 = wgam(i)/bfiten(iu)
                  v1 = exp(-x1)/x1*(cd(it)+x1*ce(it)-(ca(it)+
     &                 cb(it)*sqrt(x1)+cc(it)*x1)*log(x1))
                  v2 = exp(-x2)/x2*(cd(iu)+x2*ce(iu)-(ca(iu)+
     &                 cb(iu)*sqrt(x2)+cc(iu)*x2)*log(x2))
                  rndot= gwidth(i)*coeff*(delta*(v2-v1)+v1)
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 60            continue
               do 70 i = iswitch+1,imax
                  x1 = wgam(i)/bfiten(it)
                  x2 = wgam(i)/bfiten(iu)
                  v1 = exp(-x1)/x1*((ha(it)*x1+hb(it))*x1+hc(it)+(hd(it)
     &                 + he(it)/x1)/x1)
                  v2 = exp(-x2)/x2*((ha(iu)*x1+hb(iu))*x2+hc(iu)+(hd(iu)
     &                 + he(iu)/x2)/x2)
                  rndot = gwidth(i)*coeff*(delta*(v2-v1)+v1)
                  gder(i) = gder(i) + rndot
                  ethder = ethder - wgam(i)*rndot
 70            continue
            endif

c *** end of non-mildly rel. thermal section
         else
c ***  use  Alexanian approx, pretty good for theta > 1
c Vergessen nicht faktor von 4! (2xe-e + e+e-=2xe-e)
            coeff = 2.*.001742572369619*(0.5d0*elecdens*elecdens+
     &           0.5d0*posdens*posdens + elecdens*posdens)
            eugam = .5772156649
            imax = min(int(log10(theta*10.)/.05+RIOFF1),NPHENERG)
            iswitch = min(int(log10(theta*1.05)/.05+RIOFF1),NPHENERG)
            do 92 i = 1,iswitch
               x = wgam(i) /theta
c Achtung mit dem Sign: Ei(-x) as def'd in Alex. = -!! E_1(x)  (Ab+Steg)
c ESTUPIDO !!! no es log(2theta/eugam) aber log(2theta)-eugam!!(READ next time)
               e1 = -.57721566 + x*(.99999193 + x*(-.24991055+
     &              x*(.05519968 +x*(-.00976004 + .00107857*x))))
               e1 = e1 - log(x)
               pee=28./3.+2.*x+.5*x*x+2.*(8./3.+4./3.*x+x*x)
     &              *(log(2.*theta)-eugam)+exp(x)*e1*(8./3.-4./3.*x+x*x)
               pee=coeff*exp(-x)/wgam(i)*gwidth(i)*pee
               gder(i) = gder(i) + pee
               ethder = ethder - wgam(i)*pee
 92         continue
            do 93 i = iswitch+1,imax
               x = wgam(i) /theta
               e1=(x*x+2.334733*x+.250621)/(x*x+3.330657*x + 1.681534)
               e1 = e1/x/exp(x)
               pee=28./3.+2.*x+.5*x*x+2.*(8./3.+4./3.*x+x*x)
     &              *(log(2.*theta)-eugam)+exp(x)*e1*(8./3.-4./3.*x+x*x)
               pee=coeff*exp(-x)/wgam(i)*gwidth(i)*pee
               gder(i) = gder(i) + pee
               ethder = ethder - wgam(i)*pee
 93         continue
         endif
c ** end thermal bremss.

      endif
c endif for dothermpairs

c finally take care of non-thermal bremsstrahlung
      index = 1
c*******************  NLORGAM = 61!! careful-- should really compute
c up to NLORGAM=81
      do 10 i = 2,61
         rnprod = elecdens*(edist(i)+posdist(i))/4.d0
     &        + posdens*(edist(i)+posdist(i))/4.d0
         eloss = brem(i)*(elecdens+posdens)/2.d0
c       write(3,*)een(i),eloss
         gdot(i) = gdot(i)- eloss
c kaa  gcder never initialized or used elsewhere
c        gcder(i) = gcder(i) + eloss
c       gnder(i) = gnder(i) + dmax1(eloss,0.d0)
        do 20 j = 1,ibrem(i)
           jshift  = j+IOFF
           jj = j + index-1
           gder(jshift) = gder(jshift)+rnprod*brema(jj)
 20     continue
        index = index + ibrem(i)
 10   continue

      theta = otheta 
      return
      end

c========================================================================
c Calculate Coulomb cooling of particles
c Global variables:
c     i: elecdens, posdens, edist, posdist, radius
c     i/r: gdot, ethder

      SUBROUTINE docoul
      IMPLICIT NONE
      INCLUDE "common.eq"

      DOUBLE PRECISION ap, betas, corr, e1, echop, es, p1, ps, rlam
      DOUBLE PRECISION rn, tethder, tmin, trate, trechop, wp, x, zerox

      INTEGER i


      if (elecdens.eq.0) return

      tethder = 0.d0

      zerox=1.0 + theta/(theta+1.0)
c I guess -- I'm putting all tau_thom in here
      rn = (elecdens+posdens)/radius/6.6d-25
      wp=4*3.1415926535*4.803d-10*4.803d-10/9.1d-28*rn
c relativistic correction (approx)
      wp = sqrt(wp/(1.+3.*theta))
c      imin = max(int(log10(1.d0+4.d0*theta)/.05d0+1.5d0),2)
c to prevent bouncing around in code, need a smooth formula that
c makes gdot grow fast for  e< (1+8)theta
      echop = 1.d0+4.d0*theta
      x=(echop-1.0)/theta
      e1 = echop
      p1=sqrt(e1*e1 - 1.0)
      ps=sqrt((e1-1)/2.0)
      es=sqrt((e1+1)/2.0)
      betas=ps/es
      tmin=1.05e-27/(betas/wp)/9.1e-28/9e20/ps/es/es*(es*es+ps*ps)
      rlam=0.5-log(sqrt(2.0)*tmin)
c         print *,'rlam',rlam,log(sqrt(e1)*9.1d-28*9d20/1.05d-27/wp)
      corr= (e1-1.0)*(e1-1.0)/4.0/e1/e1*(2.0*log(2.0)+0.25)
      ap=3.0/2.0*e1/p1*(rlam+corr)
c ******************************************
c         ERROR (Moeglich) - das sollte (rlam >>+<< corr) sein
c ****************************************
c according to haug.f should be "+" so I'm changing it! PSC 8/14/97
      ap=ap/(1.0/(1.0+theta)+2.0*theta)
c    { add correction for small x}
c         if ( x < 20.0) then
      ap=ap*(x -zerox)*sqrt(x)/(x*sqrt(x)+1.0)
c  FUDGE - can go negative, because sometimes x< theta (imin approx)
c change for split e+/e- dist., Moeller and Babha scat approximately
c the same, so should be o.k.
      trechop = abs(ap*(posdens+elecdens))
      do 14 i = NLORGAM,2,-1
         x=(een(i)-1.0)/theta
         e1 = een(i)
         p1=sqrt(e1*e1 - 1.0)
         ps=sqrt((e1-1)/2.0)
         es=sqrt((e1+1)/2.0)
         betas=ps/es
         tmin=1.05e-27/(betas/wp)/9.1e-28/9e20/ps/es/es*(es*es+ps*ps)
         rlam=0.5-log(sqrt(2.0)*tmin)
c         print *,'rlam',rlam,log(sqrt(e1)*9.1d-28*9d20/1.05d-27/wp)
         corr= (e1-1.0)*(e1-1.0)/4.0/e1/e1*(2.0*log(2.0)+0.25)
         ap=3.0/2.0*e1/p1*(rlam+corr)
c ******************************************
c         ERROR (Moeglich) - das sollte (rlam >>+<< corr) sein
c ****************************************
c according to haug.f should be "+" so I'm changing it! PSC 8/14/97
         ap=ap/(1.0/(1.0+theta)+2.0*theta)
c    { add correction for small x}
c         if ( x < 20.0) then
         ap=ap*(x -zerox)*sqrt(x)/(x*sqrt(x)+1.0)
c  FUDGE - can go negative, because sometimes x< theta (imin approx)
c change for split e+/e- dist., Moeller and Babha scat approximately
c the same, so should be o.k.
c***************  CHEcK THIS
         if (een(i).lt.echop) then 
            trate = trechop*((echop/een(i))**3.6)
         else
            trate = abs(ap*(posdens+elecdens))
         endif
c  smooth  turnon!
         trate = CTURNON * trate
         gdot(i) = gdot(i) - trate
c        print *,i,een(i),ap,gdot(i),gcder(i)
c remember: TOTAL ethder
         tethder=tethder+trate*(edist(i)+posdist(i))
 14   continue
      ethder = ethder + tethder
     
      return
      end

c===================================================================================
c Global variables:
c   i: wgam, upwgam, gwidth
c   i/r: gder

      SUBROUTINE  thscat(rate,eav,disp)
      IMPLICIT NONE
      DOUBLE PRECISION rate, eav, disp

      INCLUDE "common.eq"

      DOUBLE PRECISION c, de2, dw, emax, emin, f1, f2, frac
      DOUBLE PRECISION targ, tdisp

      INTEGER i1, it, j1, k


c      dimension store(241)
c      if(eav .lt.wgam(1)) then
c        write(1,*) '!tscat: eav out of bounds--',wgam(1),eav
c        eav = wgam(1)
c      endif
      tdisp = sqrt(3.)*disp
      emin = dmax1(WMIN,eav-tdisp)
      emax = dmin1(WMAX,eav+tdisp)
      i1 = max0(int(log10(emin)/.05 + RIOFF1),1)
      j1 = min0(int(log10(emax)/.05 + RIOFF1),NPHENERG)
c      do 777 i =1,NPHENERG
c777     store(i) = gder(i)
c  o.k. let's figure out where the peak goes    
c **** Wir haben 3 an zwei geandert.
      if (abs(i1-j1) .lt. 2) then
         it = max0(int(log10(eav)/.05 + RIOFF1),1)
         if (eav .eq. wgam(it)) then
            gder(it) = gder(it) + rate
         elseif (eav .lt. wgam(it)) then
            frac =(eav-wgam(it-1))/(wgam(it-1)-wgam(it))
            gder(it) = gder(it)-frac*rate
            gder(it-1)=gder(it-1)+(1.+frac)*rate
         else
            frac = (eav-wgam(it+1))/(wgam(it)-wgam(it+1))
            gder(it) = gder(it)+frac*rate
            gder(it+1)=gder(it+1)+(1.-frac)*rate
         endif
      elseif (eav .lt. wgam(i1)) then
         it = max0(int(log10(eav)/.05 + RIOFF1),1)
         frac =(eav-wgam(it-1))/(wgam(it-1)-wgam(it))
         gder(it) = gder(it)-frac*rate
         gder(it-1)=gder(it-1)+(1.+frac)*rate
      elseif (eav .gt. wgam(j1)) then
         it = max0(int(log10(eav)/.05 + RIOFF1),1)
         frac = (eav-wgam(it+1))/(wgam(it)-wgam(it+1))
         gder(it) = gder(it)+frac*rate
         gder(it+1)=gder(it+1)+(1.-frac)*rate
      else
         dw = 1./(emax - emin)
         f1 = (upwgam(i1)-emin)*dw
         f2 = (emax - upwgam(j1-1))*dw
         de2 = (upwgam(j1-1)-upwgam(i1))*dw
         de2 = de2*.5e0*(upwgam(j1-1)+upwgam(i1))
         de2 = wgam(i1)*f1+wgam(j1)*f2+de2
         if (eav .ne. de2) then
c begin
c  stick in spike to correct things. de2 is error due to
c  discrete bin energies and edge effects
            if (eav .gt. de2) then
               targ = (emax +eav)/2.
            else
               targ = (emin + eav)/2.
            endif
            it = max0(int(log10(targ)/.05 + RIOFF1),1)
c   ** c is fraction of # scat. part that goes into "spike"
            if (targ .eq. wgam(it)) then
               c  = (eav - de2)/(targ - de2)
               gder(it) = gder(it) + c*rate
            elseif (targ .lt. wgam(it)) then
               if ((it-1) .lt. i1) then
                  c  = (eav - de2)/(wgam(i1) - de2)
                  gder(i1) = gder(i1) +c*rate
               else
                  c  = (eav - de2)/(targ - de2)
                  frac = (targ-wgam(it-1))/(wgam(it-1)-wgam(it))
                  gder(it) = gder(it)-c*frac*rate
                  gder(it-1)=gder(it-1)+c*(1.+frac)*rate
               endif
            else
               if ((it+1) .gt. j1) then
                  c  = (eav - de2)/(wgam(j1) - de2)
                  gder(j1) = gder(j1) +c*rate
               else
                  c  = (eav - de2)/(targ - de2)
                  frac = (targ-wgam(it+1))/(wgam(it)-wgam(it+1))
                  gder(it) = gder(it)+c*frac*rate
                  gder(it+1)=gder(it+1)+c*(1.-frac)*rate
               endif
            endif
c end -- now distribute the rest of the particles uniformly
         else
            c = 0.
         endif
c (if de2 = eav, we have nothing do to of course.)
         c= (1.-c)*dw*rate
         gder(i1) =gder(i1)+(upwgam(i1)-emin)*c
         gder(j1) = gder(j1)+(emax-upwgam(j1-1))*c
         do 32 k = (i1+1),(j1-1)
            gder(k)=gder(k) + gwidth(k)*c
 32      continue
      endif
c      ave = 0.
c      rpart = 0.
c      do 787 kk = 1,NPHENERG
c        rpart = rpart + gder(kk) - store(kk)
c787     ave = ave + (gder(kk)- store(kk))*wgam(kk)
c      if (abs(rpart -rate)/rate .gt. 1.e-5) then
c        write(6,*) 'another horrid problem!'
c        write(6,*) i,j,'supposed',rate,'actual',rpart
c        write(6,*) 'i,j: ',i1,j1,eav,de2,targ,c,frac
c        pause 'LUnch!'
c      endif
c      if (rpart .ne. 0.) then
c        ave = ave / rpart
c        if (abs(ave -eav)/ave .gt. 1.e-5) then
c          write(6,*) 'we have a big problem!!!'
c          write(6,*) i,j,i1,j1,'ave',ave,'eav',eav
c          pause 'CRUNCH!'
c        endif
c      endif
      return
      end



      SUBROUTINE getfname(stuff,n)
      IMPLICIT NONE
      INTEGER n
      character(8) stuff

      INTEGER len, i, id
      character(1) digit(0:9)

      stuff = 'mc'
      digit(0) = '0'
      digit(1) ='1'
      digit(2) = '2'
      digit(3) = '3'
      digit(4) = '4'
      digit(5) = '5'
      digit(6) = '6'
      digit(7)  = '7'
      digit(8) = '8'
      digit(9) = '9'
      if (n .eq. 0) then 
         len = 1
      else
         len = int(log10(float(n))+1)
      endif
      do 10 i = 0,len-1
         id = n/10**(len-1-i) 
         id = mod(id,10)
         stuff(i+3:i+3) = digit(id)
 10   continue
      stuff(3+len:3+len) = '.'
      stuff(4+len:4+len) = 'd'
      return
      end

c ----------------------------------------------
c  compute part of escaping spectrum made of black body photons that
c  do NOT scatter (presumably they are emitted from disk and
c  don't reflect off disk).

c  Global variables:
c    i: av, thann, cr, phinjfn, wgam, elecdens, posdens, edist, posdist

      SUBROUTINE bbesc(bbspect)
      IMPLICIT NONE
      INCLUDE "common.eq"

      DOUBLE PRECISION bbspect(NPHENERG)

c Local variables

      DOUBLE PRECISION tauc(NPHENERG)
      DOUBLE PRECISION qrate, tau, trate, x, x1, x2

      INTEGER i, j, ii

c probably unnecessary but removes potential compiler warning

      trate = 0

c INCLUDE CONTRIBUTION to tau from NON-THERMAL scattering too

      do i=1,NPHENERG
         tauc(i) = 0.0
         do j=2,NLORGAM
            tauc(i) = tauc(i) + (edist(j)+posdist(j))*compr(j,i)
         enddo
      enddo



      x = log10(theta)


c ** do escape!!
c ** Remember N for Tau_thomp = (ne- + ne+)=2!!!*elecdens
      tau = 1./3.*(elecdens+posdens)
      do 303 i = 1,IOFF
         bbspect(i) = phinjfn(i)/(1.d0+tau+tauc(i))
 303  continue
      do 30 i = 1+IOFF,NLORGAM+IOFF
         ii = i-IOFF
         if (theta .lt. .15) then
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
            x1 = log10(.15)
            x2 = log10(.16)
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            qrate=10.**qrate
            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
            trate=dmin1(1.d0,trate)
c         elseif (theta .ge. .16d0) then
         else
            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            trate=10.**trate
            trate=dmin1(trate,1.d0)
         endif
         tau = 1./3.*trate*(elecdens+posdens)
         bbspect(i) = phinjfn(i)/(1. + tau+tauc(i))
 30   continue
      do 309 i =82+IOFF,100+IOFF
         ii = i - IOFF
         if (theta .lt. .15) then
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
            x1 = log10(.15)
            x2 = log10(.16)
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            qrate=10.**qrate
            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
            trate=dmin1(1.d0,trate)
c         elseif (theta .ge. .16) then
         else
            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            trate=10.**trate
            trate=dmin1(trate,1.d0)
         endif
         tau = 1./2.7*trate*(elecdens+posdens)*(1.-wgam(i))
         bbspect(i) = phinjfn(i)/(1. + tau+tauc(i))
 309  continue
      do 399 i =101+IOFF,181+IOFF
         bbspect(i)= phinjfn(i)/(1.d0 + tauc(i))
 399  continue
      return
      end

c ----------------------------------------
c Note that this routine is not used.
c
c      SUBROUTINE doesc
c      IMPLICIT NONE
c      INCLUDE "common.eq"
c
c      DOUBLE PRECISION qrate, tau, trate, x, x1, x2
c
c      INTEGER i, ii
c
cc probably unnecessary but removes potential compiler warning
c
c      trate = 0
c       
c      x = log10(theta)
c      print *,'ello theta  ',theta,elecdens,posdens
c
c
cc ** do escape!!
cc ** Remember N for Tau_thomp = (ne- + ne+)=2!!!*elecdens
c      tau = 1./3.*(elecdens+posdens)
c      do 303 i = 1,IOFF
c         pout(i) = pout(i)+1./(1.+tau)
c         esc(i) = gdist(i)/(1.+tau)
c 303  continue
c      do 30 i = 1+IOFF,NLORGAM+IOFF
c         ii = i-IOFF
c         if (theta .lt. .15) then
c            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
c         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
c            x1 = log10(.15)
c            x2 = log10(.16)
c            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
c            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
c            qrate=10.**qrate
c            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
c            trate=dmin1(1.d0,trate)
cc         elseif (theta .ge. .16d0) then
c         else
c            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
c            trate=10.**trate
c            trate=dmin1(trate,1.d0)
c         endif
c         tau = 1./3.*trate*(elecdens+posdens)
c         pout(i) = pout(i)+1./(1.+tau)
c         esc(i) = gdist(i)/(1. + tau)
c 30   continue
c      do 309 i =82+IOFF,100+IOFF
c         ii = i - IOFF
c         if (theta .lt. .15) then
c            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
c         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
c            x1 = log10(.15)
c            x2 = log10(.16)
c            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
c            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
c            qrate=10.**qrate
c            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
c            trate=dmin1(1.d0,trate)
cc         elseif (theta .ge. .16) then
c         else
c            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
c            trate=10.**trate
c            trate=dmin1(trate,1.d0)
c         endif
c         tau = 1./2.7*trate*(elecdens+posdens)*(1.-wgam(i))
c         pout(i) = pout(i)+1./(1.+tau)
c         esc(i) = gdist(i)/(1. + tau)
c 309  continue
c      do 399 i =101+IOFF,181+IOFF
c         pout(i) = pout(i)+1.
c         esc(i)= gdist(i)
c 399  continue
c      return
c      end

c============================================
c Global variables:
c   i: thann, cr, wgam, upwgam, gwidth, elecdens, posdens
c   r: tesc, gder

c Calls routine thscat.

      SUBROUTINE dotherm(isthermonly)
      IMPLICIT NONE
      LOGICAL isthermonly
      INCLUDE "common.eq"

      DOUBLE PRECISION disp, eav, esum, ga, qrate, rate, rl, rlb
      DOUBLE PRECISION ru, tau, trate, ttau, wtarg, x, x1, x2, xl, xu
      DOUBLE PRECISION annspme

      INTEGER i, ii

c probably unnecessary but removes potential compiler warning

      trate = 0

c -- START: check temperature first
c      if (theta .gt. .15) write(6,*) 'Warning theta too high!',theta
c      if (theta .lt. .001) write(6,*) 'Warning theta too low!',theta

      x = log10(theta)
c*********      IDIOTSKI DON'T FORGET
      eth=(1.d0+27.d0/8.d0*theta+45.d0/8.d0*theta*theta)
     &   /(1.d0+15.d0/8.d0*theta)*(elecdens+posdens)

c --- annihilation
c   (first do annihilation of thermal bin with itself)
c use Svennson approx.
c
c zero out gder(i) so we can still thscat routine verbatim
      do 233 i=1,NPHENERG
         gder(i) = 0.
 233  continue

      IF ( .NOT.isthermonly ) THEN

         if (theta .lt. 0.005) then
            ga = 1./(1.+2.*theta*theta/log(1.120*theta+1.3))
            trate = ga*.375e0
            rate = trate*elecdens*posdens
            if (elecdens .gt. 1.e-12) then
               eav = eth/(elecdens+posdens)
            else
               eav = 1+1.5*theta
            endif
            disp = thann(3,1)*theta+thann(3,2)
            call thscat(2.d0*rate,eav,disp)
         else
c use Svensson (83)
            xl = (1.0-sqrt(3.0)/3.0)/2.0
            xu = (1.0+sqrt(3.0)/3.0)/2.0
            ttau = elecdens*posdens
            do 510 i = 80+IOFF,142+IOFF
               rlb = upwgam(i-1)
               wtarg = rlb + gwidth(i)*xl
               rl =  annspme(wtarg,theta)
               wtarg = rlb + gwidth(i)*xu
               ru = annspme(wtarg,theta)
c    AZ annspme gives answer for unit opt. depth     
c    if elecdens = 1, => tau=2. AZ is factor 4 too small
               trate =  (ru+rl)*gwidth(i)
               rate =  trate*ttau*2.0
               gder(i) = gder(i) + rate
 510        continue
         endif

      ENDIF
      
      esum = 0
      do i=1,NPHENERG
         esum = esum + gder(i)*wgam(i)
      enddo


c ---- PAIR Escape!!!

c ** initialize escape array
c ** Remember N for Tau_thomp = (ne- + ne+)=2!!!*elecdens
      tau = 1./3.*(elecdens+posdens)
      do 303 i = 1,IOFF
         tesc(i) = 1.d0/(1.d0+tau)
 303  continue
      do 30 i = 1+IOFF,NLORGAM+IOFF
         ii = i-IOFF
         if (theta .lt. .15) then
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
            x1 = log10(.15)
            x2 = log10(.16)
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            qrate=10.**qrate
            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
            trate=dmin1(1.d0,trate)
c         elseif (theta .ge. .16d0) then
         else
            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            trate=10.**trate
            trate=dmin1(trate,1.d0)
         endif
         tau = 1./3.*trate*(elecdens+posdens)
         tesc(i) = 1./(1.+tau)
 30   continue
      do 309 i =82+IOFF,100+IOFF
         ii = i - IOFF
         if (theta .lt. .15) then
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
            x1 = log10(.15)
            x2 = log10(.16)
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            qrate=10.**qrate
            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
            trate=dmin1(1.d0,trate)
c         elseif (theta .ge. .16) then
         else
            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            trate=10.**trate
            trate=dmin1(trate,1.d0)
         endif
         tau = 1./2.7*trate*(elecdens+posdens)*(1.-wgam(i))
         tesc(i) = 1./(1.+tau)
 309  continue
      do 399 i =101+IOFF,181+IOFF
         tesc(i)= 1.d0
 399  continue
      return
      end

c ********************************************************************************
c Calculate photon scattering matrix and return in scat.
c Global variables:
c    i: ecomp, gcomp, iecomp, igcomp, cltr, clta, cltd, wgam, upwgam, gwidth,
c       elecdens, posdens

      SUBROUTINE mkscat(pout, scat, iscat)
      IMPLICIT NONE
      INCLUDE "common.eq"

      DOUBLE PRECISION pout(NPHENERG)
      DOUBLE PRECISION scat(NPHENERG,NPHENERG)
      INTEGER iscat(NPHENERG,2)

c Local variables

      DOUBLE PRECISION exn(141)
      DOUBLE PRECISION wi(15),a(6),b(6),xt(6)
      DOUBLE PRECISION dn1, dw, eav, eavp, emax, emin, frac, grate
      DOUBLE PRECISION ratdelt, rate, rm1, rm2, rnp, tdisp, tdn1
      DOUBLE PRECISION thint, trate, trnp, weight, xp
      
      INTEGER ien(15)
      INTEGER i, i1, ii, ip1, j, j1, jc, jj, k, lm

      LOGICAL qfirst
      
      SAVE exn, a, b, qfirst

c * xi values
      data  a /0.222846604179,
     *     1.188932101673,
     *     2.992736326059,
     *     5.775143569105,
     *     9.837467418383,
     *     15.982873980602/
c * wi values   
      data b /4.58964673950e-1,
     *     4.17000830772e-1,
     *     1.13373382074e-1,
     *     1.03991974531e-2,
     *     2.61017202815e-4,
     *     8.98547906430e-7/

      data qfirst /.TRUE./

c Initialize energy array used to interpolate on clt* arrays

      IF ( qfirst ) THEN
         DO i = 1, 141
            exn(i) =1.0+10.**(float(i-101)*.05)
         ENDDO
         qfirst = .FALSE.
      ENDIF
    
c
c --- compton scattering
c(Assume e- aren't scattered out of thermal dist.)
c --- Also do photon escape (only count thermal e+/e- in
c      determining escape rate. Escape rate formula is one
c      used by big Z.)

      do 43 i=1,NPHENERG-20
         iscat(i,1) = NPHENERG
         iscat(i,2) = 1
         do 44 j=1,NPHENERG 
            scat(i,j) = 0.d0
 44      continue
 43   continue
      do i=NPHENERG-20,NPHENERG
         do j=1,NPHENERG
            scat(i,j)  = 0.d0
         enddo
      enddo


c *** Use SLEDGE hammer approach for everything, ie. the laguerre
c   integration fudge. (may be a bit of overkill, but we must get
c   rid of glitches und have Ordnung hier!)
c*** sledge hammer approach
      thint = 0.0
      do 189 lm = 1,6
         xp = a(lm)*theta + 1.0
         xt(lm) = xp
         if (xp .ge. 101.) then 
            ien(lm) = max(int(log10(xp)/.05+1.5),1)
         else
            ien(lm)=-max(int(log10(xp-1.0)/.05+101),1)
         endif
         wi(lm) = b(lm)*xp*sqrt(xp*xp-1.0)
         thint = thint + wi(lm)
 189  continue
      do 190 lm=1,6
         wi(lm)=wi(lm)/thint*(elecdens+posdens)
190   continue
      do 310 ii = 1,6
         i = ien(ii)     
         weight = wi(ii)
         if (i .gt. 0) then 
c ====     use non-thermal comp. data
c ignore photons > x=10^3
            do 320 j = 1,NPHENERG-IOFF-20
               jj = j + IOFF
               grate = ecomp(i,j,1)*weight
               pout(jj) = pout(jj) + grate
c Take care of scattering of PHOTONS into other bins
               i1 = igcomp(i,j,1)+IOFF
               j1 = igcomp(i,j,2)+IOFF
               if (i1 .lt. iscat(jj,1)) iscat(jj,1) = i1
               if (j1 .gt. iscat(jj,2)) iscat(jj,2) = j1
               dn1 = gcomp(i,j,2)
               rnp = gcomp(i,j,3)
               if (i1 .eq. j1) then
                  scat(jj,i1) = scat(jj,i1) + grate
               elseif (abs(j1-i1) .eq. 1) then
                  scat(jj,i1)= scat(jj,i1) + grate*dn1
                  scat(jj,j1) = scat(jj,j1) + grate*(1. - dn1)
               else
                  trnp = grate*rnp
                  tdn1 = grate*dn1
                  dw = trnp/(upwgam(j1-1) - upwgam(i1))
                  do 340 k = (i1+1),(j1-1)
                     scat(jj,k) = scat(jj,k) + dw*gwidth(k)
 340              continue
                  scat(jj,i1) = scat(jj,i1) + tdn1
                  scat(jj,j1) = scat(jj,j1) + (grate  - tdn1) - trnp
               endif
 320        continue

c  go do low energy comp. the same way
c  NB: All rates are assumed to be=sigma_t, since should be in
c       Thompson limit.
c  *** 80 should be ibot !!!! (in program with lcinit)
c achtung, yup, but knocked off two decades
            do 420 j = 1,40
               rate = weight
               dn1 = rldisp(i,j,1)
               rnp = rldisp(i,j,2)
               pout(j) = pout(j) + weight
               i1 = ig(i,j,1)
               j1 = ig(i,j,2)
               if (i1.lt.iscat(j,1)) iscat(j,1) = i1
               if (j1.gt.iscat(j,2)) iscat(j,2) = j1
               if (i1 .eq. j1) then
                  scat(j,i1) = scat(j,i1) + rate
               elseif (abs(j1-i1) .eq. 1) then
                  scat(j,i1)= scat(j,i1) + rate*dn1
                  scat(j,j1) = scat(j,j1) + rate*(1. - dn1)
               else
                  trnp = rate*rnp
                  tdn1 = rate*dn1
                  dw = trnp/(upwgam(j1-1) - upwgam(i1))
                  do 490 k = (i1+1),(j1-1)
                     scat(j,k) = scat(j,k) + dw*gwidth(k)
 490              continue
                  scat(j,i1) = scat(j,i1) + tdn1
                  scat(j,j1) = scat(j,j1) + (rate  - tdn1) - trnp
               endif
 420        continue

c == end of non-therm. data use
         else
c == use extended, low gamma data
            i = -i
            ip1 = min(i+1,141)
            if (i .eq. ip1) then
               frac = 0
            else
c BUGUGUGU!!!!!!!!!!!!!!!!!!!!!!!!!????????
c         frac=(xt(ii)-exn(ip1))/(exn(ip1)-exn(i))
c YES, idiotski!!! DOUBLE CHECK  
               frac=(xt(ii)-exn(i))/(exn(ip1)-exn(i))
            endif
            do 321 j=1,NPHENERG -20
               jc = j + 40
               trate=cltr(i,jc)+frac*(cltr(ip1,jc)-cltr(i,jc))
               eav  =clta(i,jc)+frac*(clta(ip1,jc)-clta(i,jc))
               tdisp=cltd(i,jc)+frac*(cltd(ip1,jc)-cltd(i,jc))
c         print *,j,trate,eav,tdisp
c         b2 = (xt(ii)*xt(ii)-1)/xt(ii)/xt(ii)
c         print *,'ap',wgam(j)*(1.+4./3.*b2-wgam(j))
c         test=test+eav*trate*weight/(elecdens+posdens)
               grate = trate*weight
               pout(j) = pout(j) + grate
c go scatter!  ============================
               emin = eav-tdisp
               emax = eav+tdisp
               if (emax .lt. emin) then
                  emin = eav
                  emax = eav
               endif
               if (emin .lt. WMIN) then 
                  tdisp = abs(eav - WMIN - 1.e-10)
                  emax = eav+tdisp
                  emin = eav - tdisp
               endif
               if (emax .gt. WMAX) then
                  tdisp = abs(WMAX - eav - 1.e-10)
               endif
               emax = eav+tdisp
               emin = eav - tdisp
               i1 = max0(int(log10(emin)/.05 + RIOFF1),1)
               j1 = min0(int(log10(emax)/.05 + RIOFF1),NPHENERG)
c      do 777 klm =1,NPHENERG
c777     store2(klm) = gder(klm)

c -- changed from ge
               if (wgam(i1) .ge. eav) then
                  i1 = i1 -1  
                  if (i1 .lt. 1) then
                     scat(j,1) = scat(j,1) + grate
                     print *,'ERROR in thscat - eav < wgam(1)',eav
                     goto 321
                  endif
               endif
c -- changed from le
               if (wgam(j1) .le. eav) then
                  j1 = j1 + 1
                  if (j1 .gt. NPHENERG) then
                     scat(j,NPHENERG)=scat(j,NPHENERG)+grate
                     print *,'ERROR in thscat - eav > wgam(162)',eav
                     goto 321
                  endif
               endif
c      print *,'i1 = ',i1
c      print *,'j1 = ',j1
c      if (i1 .eq. j1) pause 'rat Foo!'
               if (i1.lt.iscat(j,1)) iscat(j,1) = i1
               if (j1.gt.iscat(j,2)) iscat(j,2) = j1
               if (abs(j1-i1) .eq. 1) then
                  dn1 = grate*(eav - wgam(j1))/(wgam(i1) - wgam(j1))
                  scat(j,i1)= scat(j,i1) + dn1
                  scat(j,j1) = scat(j,j1) + grate - dn1
               else
                  eavp = .5*(upwgam(i1)+ upwgam(j1-1))
                  dw = 1./(emax-emin) 
c        if (eavp .eq. wgam(j1)) pause 'Here it is!'
                  rm1 = (eav - wgam(j1))/(eavp - wgam(j1)) 
c     if (eavp .eq. wgam(i1)) pause 'no its here!'
                  rm2 = (eav - wgam(i1))/(eavp - wgam(i1))
                  ratdelt = upwgam(j1-1) - upwgam(i1)
                  dw = min(dw*ratdelt,rm1,rm2)/ratdelt
                  rnp = grate*ratdelt*dw
                  do 30 k = (i1+1),(j1-1)
                     scat(j,k) = scat(j,k) + 
     &                    grate*dw*(upwgam(k) - upwgam(k-1))
 30               continue
                  dn1 = grate*eav - rnp*eavp + wgam(j1)*(rnp - grate)
                  dn1 = dn1/(wgam(i1)-wgam(j1))
                  scat(j,i1) = scat(j,i1) + dn1
                  scat(j,j1) = scat(j,j1) + (grate  - dn1) - rnp
               endif
c      ave = 0.
c      rpart = 0.
c      do 787 kk = 1,NPHENERG
c        rpart = rpart + gder(kk) - store2(kk)
c787     ave = ave + (gder(kk)- store2(kk))*wgam(kk)
c      print *,rpart,grate
c      print *,ave/rpart,eav
c      if (abs(rpart -grate)/grate .gt. 1.e-5) then
c        write(6,*) 'another horrid problem!'
c        write(6,*) i,j,'supposed',grate,'actual',rpart
c        write(6,*) 'i,j: ',i1,j1,eav,de2,targ,c,frac
c        pause 'LUnch!'
c      endif
c      if (rpart .ne. 0.) then
c        ave = ave / rpart
c        if (abs(ave -eav)/ave .gt. 1.e-5) then
c          write(6,*) 'we have a big problem!!!'
c          write(6,*) i,j,i1,j1,'ave',ave,'eav',eav
c          pause 'CRUNCH!'
c        endif
c      endif
 321        continue
         endif
 310  continue

      do 343 i=1,NPHENERG-20
         pout(i) = pout(i) - scat(i,i)
         scat(i,i) =  0.0
 343  continue

      return
      end

c=============================================================================
c Calculate change in thermal energy
c Global variables:
c   i: rldisp, ig, ecomp, gcomp, iecomp, igcomp, clav, ectab, iectab, wgam,
c      upwgam, gwidth, gdist, elecdens, posdens.
c   i/r: ethder

      SUBROUTINE getethd(isthermonly, eout2)
      IMPLICIT NONE
      INCLUDE "common.eq"

      LOGICAL isthermonly
      DOUBLE PRECISION eout2

c Local variables

      DOUBLE PRECISION exn(141)
      DOUBLE PRECISION wi(15),a(6),b(6),xt(6)
      DOUBLE PRECISION eav, frac, ga, grate, rate, rl, rlb, ru
      DOUBLE PRECISION thint, trate, ttau, weight, wtarg
      DOUBLE PRECISION x, xl, xp, xu, zrate, annspme

      INTEGER ien(15)
      INTEGER i, ii, ip1, j, jc, jj, lm

      LOGICAL qfirst

      SAVE exn, a, b, qfirst

c * xi values
      data  a /0.222846604179,
     *     1.188932101673,
     *     2.992736326059,
     *     5.775143569105,
     *     9.837467418383,
     *     15.982873980602/
c * wi values   
      data b /4.58964673950e-1,
     *     4.17000830772e-1,
     *     1.13373382074e-1,
     *     1.03991974531e-2,
     *     2.61017202815e-4,
     *     8.98547906430e-7/

      data qfirst /.TRUE./

c Initialize energy array used to interpolate on clt* arrays

      IF ( qfirst ) THEN
         DO i = 1, 141
            exn(i) =1.0+10.**(float(i-101)*.05)
         ENDDO
         qfirst = .FALSE.
      ENDIF


      x = log10(theta)


c*********         IDIOTSKI DON'T FORGET
      eth=(1.d0+27.d0/8.d0*theta+45.d0/8.d0*theta*theta)
     &   /(1.d0+15.d0/8.d0*theta)*(elecdens+posdens)

      IF ( .NOT.isthermonly ) THEN

c --- annihilation
c   (first do annihilation of thermal bin with itself)
c use Svennson approx.
c
c *** remember ethder is TOTAL change in thermal energy 
         eout2 = 0.
         if (theta .lt. 0.005) then
            ga = 1./(1.+2.*theta*theta/log(1.120*theta+1.3))
            trate = ga*.375e0
            rate = trate*elecdens*posdens
            if (elecdens .gt. 1.e-12) then
               eav = eth/(elecdens+posdens)
            else
               eav = 1+1.5*theta
            endif
            eout2 = eout2 + trate
            ethder = ethder - rate*eav*2.
         else
c     use Svensson (83)
            xl = (1.0-sqrt(3.0)/3.0)/2.0
            xu = (1.0+sqrt(3.0)/3.0)/2.0
            ttau = elecdens*posdens
            do 510 i = 80+IOFF,142+IOFF
               rlb = upwgam(i-1)
               wtarg = rlb + gwidth(i)*xl
               rl =  annspme(wtarg,theta)
               wtarg = rlb + gwidth(i)*xu
               ru = annspme(wtarg,theta)
c    AZ annspme gives answer for unit opt. depth     
c    if elecdens = 1, => tau=2. AZ is factor 4 too small
               trate =  (ru+rl)*gwidth(i)
               rate =  trate*ttau*2.0
               eout2 = eout2 + trate
c  ATTENTO! ACHTUNG -- Gefahr
               ethder = ethder - rate*wgam(i)
 510        continue
         endif
c debug
         ga = 1./(1.+2.*theta*theta/log(1.120*theta+1.3))
         zrate = ga*.375e0
         zrate = zrate

      ENDIF

c
c --- compton scattering
c(Assume e- aren't scattered out of thermal dist.)
c --- Also do photon escape (only count thermal e+/e- in
c      determining escape rate. Escape rate formula is one
c      used by big Z.)

c *** Use SLEDGE hammer approach for everything, ie. the laguerre
c   integration fudge. (may be a bit of overkill, but we must get
c   rid of glitches und have Ordnung hier!)
c*** sledge hammer approach
      thint = 0.0
      do 189 lm = 1,6
         xp = a(lm)*theta + 1.0
         xt(lm) = xp
         if (xp .ge. 101.) then 
            ien(lm) = max(int(log10(xp)/.05+1.5),1)
         else
            ien(lm)=-max(int(log10(xp-1.0)/.05+101),1)
         endif
         wi(lm) = b(lm)*xp*sqrt(xp*xp-1.0)
         thint = thint + wi(lm)
 189  continue
      do 190 lm=1,6
         wi(lm)=wi(lm)/thint*(elecdens+posdens)
 190  continue
      do 310 ii = 1,6
         i = ien(ii)     
         weight = wi(ii)
         if (i .gt. 0) then 
c ====     use non-thermal comp. data
c ignore photons > x=10^3
            do 320 j = 1,NPHENERG-IOFF-20
               jj = j + IOFF
               grate = ecomp(i,j,1)*weight*gdist(jj)
               ethder = ethder + (wgam(jj)-av(i,jj))*grate
 320        continue
c  go do low energy comp. the same way
c  NB: All rates are assumed to be=sigma_t, since should be in
c       Thompson limit.
c  *** 80 should be ibot !!!! (in program with lcinit)
c achtung, yup, but knocked off two decades
            do 420 j = 1,40
               rate = weight*gdist(j)
c surely this should be rate not grate ? - kaa
c               ethder = ethder + (wgam(j)-clav(i,j))*grate
               ethder = ethder + (wgam(j)-clav(i,j))*rate
 420        continue
c == end of non-therm. data use
         else
c     == use extended, low gamma data
            i = -i
            ip1 = min(i+1,141)
            if (i .eq. ip1) then
               frac = 0
            else
c BUGUGUGU!!!!!!!!!!!!!!!!!!!!!!!!!????????
c         frac=(xt(ii)-exn(ip1))/(exn(ip1)-exn(i))
c YES, idiotski!!! DOUBLE CHECK
               frac=(xt(ii)-exn(i))/(exn(ip1)-exn(i))
            endif
            do 321 j=1,NPHENERG -20
               jc = j + 40
               trate=cltr(i,jc)+frac*(cltr(ip1,jc)-cltr(i,jc))
               eav  =clta(i,jc)+frac*(clta(ip1,jc)-clta(i,jc))
               grate = trate*weight*gdist(j)
               ethder = ethder + (wgam(j)-eav)*grate
 321        continue
         endif
 310  continue

      return
      end

c========================================================================
c Derive the electron distribution assuming the current photon distribution 
c (in gdist).
c Global variables:
c    i: av, clav, ectab, iectab, ppr, ppf, ipp, gdist
c    i/r: edist, posdist
c    r: gdot, ethder

      SUBROUTINE doelec(isthermonly)
      IMPLICIT NONE
      INCLUDE "common.eq" 

      LOGICAL isthermonly

c Local variables

      DOUBLE PRECISION ein(NLORGAM), posin(NLORGAM)
      DOUBLE PRECISION eout(NLORGAM), posout(NLORGAM)
      DOUBLE PRECISION ein0(NLORGAM), posin0(NLORGAM)
      DOUBLE PRECISION eout0(NLORGAM), posout0(NLORGAM)
      DOUBLE PRECISION deltag(NLORGAM)
      DOUBLE PRECISION escat(NLORGAM,NLORGAM)
      DOUBLE PRECISION tescat(NLORGAM,NLORGAM)
      DOUBLE PRECISION dmax, dth, erate, frac, superoth, temp, tfr
      DOUBLE PRECISION psum

      INTEGER i, j, ii, i1, i2, index, isth, k

      LOGICAL qdone

c initialize ein and posin particle distribution(s) to
c injected distributions

      ethder = 0.

      DO i=1,NLORGAM
         gdot(i) = 0.
         ein(i) = pairinjfn(i)
c     for pair injection
         if (qpairinj) then 
            posin(i) = pairinjfn(i)
c     for electron and/or positron injection
         else
            posin(i) = epinjfn(i)
         endif
      ENDDO

c Initialize work arrays

      DO i = 1, NLORGAM
         eout(i) = 0.0
         posout(i) = 0.0
      ENDDO

c -- set handy array of bin widths

      DO i=2,NLORGAM
         deltag(i) = een(i)-een(i-1)
      ENDDO


c get sources
c  -----------  pair production
c modify ein and posin arrays by adding pairs created from photons
c in gdist array

      superoth = theta
      isth = 0

      psum = 0.
      IF ( .NOT.isthermonly ) THEN
         DO i = 141,NPHENERG
            DO j = 1,i
               ii = i-140
               IF (ppr(ii,j).GE. .005 .AND. ipp(ii,j) .GE. 1) THEN
                  erate = gdist(i)*gdist(j)*ppr(ii,j)
                  i1 = ipp(ii,j)
                  i2 = i1 + 1
                  frac = ppf(ii,j)
                  tfr = frac*erate
                  ein(i1) = ein(i1) + tfr
                  ein(i2) = ein(i2) + erate - tfr
                  posin(i1) = posin(i1) + tfr
                  posin(i2) = posin(i2) + erate-tfr
                  psum = psum + erate
               ENDIF
            ENDDO
         ENDDO
      ENDIF

c get thomson optical depth
      call getth(isthermonly, ein(1), posin(1), psum)

      IF ( isthermonly ) THEN

         do i=1,1
            ein0(i) = ein(i)
            posin0(i) = posin(i)
            eout0(i) = eout(i)
            posout0(i) = posout(i)
         enddo

      ELSE

c     setup scattering matrix
         index = 0
         do i=1,NLORGAM
            do j=1,NLORGAM
               escat(i,j) = 0.0
            enddo
         enddo

c for my scheme
         i=1
         do j=1,NPHENERG
            do k=iectab(i,j,1),iectab(i,j,2)
               index=index+1
            enddo
         enddo
      
         do i=2,NLORGAM
            do j=1,NPHENERG
               do k=iectab(i,j,1),iectab(i,j,2)
                  index=index+1
                  escat(k,i) =escat(k,i)+ectab(index)*gdist(j)
               enddo
            enddo
         enddo

         do i=2,NLORGAM
            do j=1,NPHENERG
               if(compr(i,j).gt.0) then
                  eout(i) = eout(i) + compr(i,j)*gdist(j)
               endif
            enddo
            eout(i) = eout(i) - escat(i,i)
            escat(i,i) = 0.0
            posout(i) = eout(i)
         enddo

         do i=1,NLORGAM
            ein0(i) = ein(i)
            posin0(i) = posin(i)
            eout0(i) = eout(i)
            posout0(i) = posout(i)
         enddo
         do  i = 2,NLORGAM
            tescat(i-1,i) = escat(i-1,i)
         enddo

      ENDIF

c Start of loop searching for theta.

 1    continue

      IF ( isthermonly ) THEN

         do i=1,1
            ein(i) = ein0(i)
            posin(i) = posin0(i)
         enddo

         dmax = 0

      ELSE

c always call bremsstrahlung, and if are doing thermal distribution,docoul
c REMEMBER to zero out GDOT!
         do i=1,NLORGAM
            gdot(i) = 0.
         enddo
         call dobrem
         if ( dothermpairs ) then
c turn on cooling smoothly or else, boom!
            call docoul
         endif

c  dumb, first order advection scheme which explicitly conserves
c  energy.
         DO i = 2,NLORGAM
            if (gdot(i) .le. 0.) then
               frac = gdot(i)/deltag(i)
               eout(i) =  eout(i) -frac
               posout(i) = posout(i) - frac
               escat(i-1,i) = escat(i-1,i) - frac
            elseif (i .lt. NLORGAM) then
               CALL xwrite('positive gdot?',10)
               frac = gdot(i)/deltag(i+1)
               eout(i) = eout(i) + frac
               escat(i+1,i)=escat(i+1,i) + frac
            endif
         ENDDO

         qdone = .FALSE.
         DO WHILE (.NOT.qdone)
 
            do i=1,NLORGAM
               ein(i) = ein0(i)
               posin(i) = posin0(i)
               do j=2,NLORGAM
                  ein(i) = ein(i)+escat(i,j)*edist(j)
                  posin(i) = posin(i) + escat(i,j)*posdist(j)
               enddo
            enddo
 
            dmax = 0
      
            do i=2,NLORGAM
               if (eout(i).ne.0.d0) then 
                  temp = ein(i)/eout(i)
                  dmax = max(abs(temp-edist(i))/(edist(i)+1.d-40),dmax)
                  edist(i) = temp
               endif
               if (posout(i).ne.0.d0) then
                  temp = posin(i)/posout(i)
                  dmax = max(abs(temp-posdist(i))/(posdist(i)+1.d-40),
     &                       dmax)
                  posdist(i) = temp
               endif
            enddo

            if (dmax.LE. 1.d-4) qdone = .TRUE.

         ENDDO

      ENDIF

      call getth(isthermonly, ein(1), posin(1), psum)

      if ( superoth .EQ. 0.0 ) GOTO 928
      dth = abs(theta-superoth)/superoth
      superoth = theta
ccdebug
cc      write(*,'(2(a,1pg13.6))') 'In doelec loop: theta = ', 
cc     &                          theta, ' dth = ', dth
ccdebug

      isth=isth + 1
      if ((dth .lt. 1.d-2)) goto 928

      IF ( isthermonly ) THEN

         do i=1,1
            ein(i) = ein0(i)
            posin(i) = posin0(i)
            eout(i) = eout0(i)
            posout(i) = posout0(i)
         enddo

      ELSE

         do i=1,NLORGAM
            ein(i) = ein0(i)
            posin(i) = posin0(i)
            eout(i) = eout0(i)
            posout(i) = posout0(i)
         enddo
         do  i = 2,NLORGAM
            escat(i-1,i) = tescat(i-1,i)
         enddo

      ENDIF

c End of loop searching for theta.

      goto 1

928   continue

c      stop = dtime(time)
c      print *,'e- time=',time(1)

      return
      end

      double precision function annspme(x, temp)
      implicit double precision (a-h,o-z)
c     This program calculates the annihilation spectrum from a thermal
c     pure pair plasma of unit Thomson depth using the formulae by
c     Svensson 1983.
c     x and temp are in units of m c**2. The result gives photon number
c     spectrum in units of sigma(Thomson) c m c**2/ m c**2
      pi=3.14159
      a1=(x/temp)**2 /2 /(1+2.0049/temp +1.4774/temp/temp
     >     +pi/(2*temp)**3)
      y=x*temp
      if(y.le.4.) then
         a2=pi/2 *sqrt(pi/y) *exp((2-x-1/x)/temp)
         if(y.le.0.25) then
            y4=4*y
            cc=1+0.5288*y4-0.4483*y4**2+0.1643*y4**3
         else if(y.le.1.) then
            cc=1.125+0.6600*y-0.7972*y**2+0.3060*y**3
         else
            cc=0.9049+1.1613/y-1.2487/y**2+0.4770/y**3
         end if
      else
         a2=pi/y*exp((2-x)/temp)*(log(4*0.56146*y)-1)
         y4=4/y
         cc=1+0.5289*y4-0.8254*y4**2+0.9811*y4**3-0.3895*y4**4
      end if
c  3/32/pi follow from assuming n+=n- and using sigma(Thomson)
      annspme=a1*a2*cc/temp*3/32./pi
      return
      end


c ===================================

      SUBROUTINE thcdopcomp

      IMPLICIT NONE

      INCLUDE "common.eq" 

      DOUBLE PRECISION scat(NPHENERG,NPHENERG)
      DOUBLE PRECISION pin(NPHENERG), pout(NPHENERG)
      DOUBLE PRECISION pin0(NPHENERG), pout0(NPHENERG)
      DOUBLE PRECISION dmax, odmax, temp, thresh

      INTEGER iscat(NPHENERG,2)

      INTEGER i, j, idmax

c      print *,'o.k. bosses, theta,elecdens',theta,elecdens

      thresh = 1.d-5

c clear pout and set pin to phinjfn(i)

      DO i = 1,NPHENERG
         pout(i) = 0.
         pin(i) = phinjfn(i)
      ENDDO

c assume e- dist fixed


c call annihilation + brem + anything else that makes photons
      do i=1,NPHENERG
         gder(i) = 0.d0
      enddo
      call thcdotherm

      do i=1,NPHENERG
         pout(i) = pout(i) + tesc(i)
         pin(i) = pin(i) + gder(i)
         pin0(i) = pin(i)
      enddo

c       print *,'called dotherm'

      call mkscat(pout, scat, iscat)
c       print *,'made scattering matrix'


      do i=1,NPHENERG
         pout0(i) = pout(i)
      enddo

c  now iterate


 100  continue

       do i=1,NPHENERG-20
          do j=iscat(i,1),iscat(i,2)
             pin(j) = pin(j) + scat(i,j)*gdist(i)
          enddo
       enddo


      dmax = 0
      do i=1,NPHENERG
         temp = pin(i)/pout(i)
         odmax = dmax
         dmax = max(abs(temp-gdist(i))/(gdist(i)+1.d-40),dmax)
         if (dmax.ne.odmax) idmax = i
         gdist(i) = temp
      enddo

c      print *,'thresh?,dmax',thresh,dmax


      if (dmax .gt. thresh) Then
         do i=1,NPHENERG
            pout(i) = pout0(i)
            pin(i) = pin0(i)
         enddo
         goto 100
      endif



c      print *,'RETURNING!!'
      return
      end


      subroutine thcsolve

      IMPLICIT NONE

      include "common.eq"

      INTEGER i

      dothermpairs = .TRUE.

c  --- (assume theta, eth, elecdens initialized in setup! Ja?)

      do i=1,NPHENERG
         otesc(i) = 1.
      enddo

c   -- now do photon distribution
      call thcdopcomp

      do i=1,NPHENERG
         esc(i) = tesc(i)*gdist(i)
         otesc(i) = tesc(i)
      enddo

      RETURN
      end 



c ----------------------------------------

      SUBROUTINE thcdotherm

      IMPLICIT NONE

      INCLUDE "common.eq"

      DOUBLE PRECISION trate, qrate, tau, x, x1, x2

      INTEGER i, ii


c -- START: check temperature first
c      if (theta .gt. .15) write(6,*) 'Warning theta too high!',theta
c      if (theta .lt. .001) write(6,*) 'Warning theta too low!',theta

      x = log10(theta)
c*********      IDIOTSKI DON'T FORGET
      eth=(1.d0+27.d0/8.d0*theta+45.d0/8.d0*theta*theta)
     &   /(1.d0+15.d0/8.d0*theta)*(elecdens+posdens)


c ** initialize escape array
c ** Remember N for Tau_thomp = (ne- + ne+)=2!!!*elecdens
      tau = 1./3.*(elecdens+posdens)
      do 303 i = 1,IOFF
         tesc(i) = 1.d0/(1.d0+tau)
 303  continue
      do 30 i = 1+IOFF,NLORGAM+IOFF
         ii = i-IOFF
         if (theta .lt. .15) then
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
            x1 = log10(.15)
            x2 = log10(.16)
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            qrate=10.**qrate
            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
            trate=dmin1(1.d0,trate)
c         elseif (theta .ge. .16d0) then
         else
            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            trate=10.**trate
            trate=dmin1(trate,1.d0)
         endif
         tau = 1./3.*trate*(elecdens+posdens)
         tesc(i) = 1./(1.+tau)
 30   continue
      do 309 i =82+IOFF,100+IOFF
         ii = i - IOFF
         if (theta .lt. .15) then
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
         elseif ((theta .ge. .15).and.(theta .lt. .16)) then
            x1 = log10(.15)
            x2 = log10(.16)
            trate =dmin1(cr(ii,1)*theta+cr(ii,2),1.d0)
            qrate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            qrate=10.**qrate
            trate = trate+(qrate-trate)/(x2-x1)*(x-x1)
            trate=dmin1(1.d0,trate)
c         elseif (theta .ge. .16) then
         else
            trate=htcr(ii,1)+x*(htcr(ii,2)+x*(htcr(ii,3)+x*htcr(ii,4)))
            trate=10.**trate
            trate=dmin1(trate,1.d0)
         endif
         tau = 1./2.7*trate*(elecdens+posdens)*(1.-wgam(i))
         tesc(i) = 1./(1.+tau)
 309  continue
      do 399 i =101+IOFF,181+IOFF
         tesc(i)= 1.d0
 399  continue
      return
      end



