      subroutine dospin(EAR,NE,PARAM,IFL,PHOTAR,PHOTER)
c
c     XSPECv11 Subroutine to produce an emission line profile from a thin
c     keplerian disk around a Kerr black hole with arbitrary spin
c     parameter.   Method relies on the transfer function approach
c     of Cunningham (1975); the line profile can be computed via a
c     simple integral with a analytic integrand apart from the
c     effects of gravitational light bending/lensing.  This are
c     incorporated into a slowly varying transfer function computed
c     via the code of Speith et al. (1995) and sparsely tabulated in the
c     accompanying table "kerrtable.fits".   High quality results are
c     obtained despite the sparse sampling of the TF through
c     linear interpolation of this slowing varying function.
c
c     This code also generates the basic relativistic disk kernel
c     used for the convolution model kerrconv.
c
c     Two features have been hardwired into this code but are simple
c     to change.  Firstly, we assume that the line emissivity is described
c     with a broken powerlaw (see comment starting "EMISSIVITY PROFILE").
c     It is trivial to include any functional form.  Secondly, we have
c     used the same limb-darkening rule as used in Laor (1991; see
c     comment starting "LIMB DARKENING").  This is readily changed to any
c     functional form.
c
c     At the current time, kerrdisk and kerrconv do not allow for emission
c     within the innermost stable circular orbit
c
c     For other details of this code see Brenneman & Reynolds, 2006,
c     ApJ, volume 652, pp.1028.  Please reference this paper if you
c     publish results dereived from this model.
c
c     Developed by Laura Brenneman and Chris Reynolds
c     Dept. of Astronomy, University of Maryland, College Park
c
c
c     Modified Mar 09 for xspec12: calls fgmodf instead of using
c        xspec11's LMODDIR setting.  (CG)
c
c     Modified Jan 20: Replaced NR-based spingauleg routine with calls 
c        to fgsl library
c       
      use fgsl     
      IMPLICIT NONE
      INTEGER IFL, NE
      REAL EAR(0:NE), PARAM(9), PHOTAR(NE), PHOTER(NE)

c
c---------Initialization--------------
c
      integer i,j,k,ii,jj,ii1,ii2,igstar2,irad,ilgrad,ilgrad_max
      integer nradii,ng,ia,imu0,abins,mu0bins,ilun,ios,block,hdutyp
      integer ir
      parameter(nradii=50,ng=20,abins=20,mu0bins=20)
      DOUBLE PRECISION sumspec
      DOUBLE PRECISION a,theta0,mu0,gstar(ng)
      DOUBLE PRECISION a_tab(abins),mu0_tab(mu0bins)
      DOUBLE PRECISION aintfac,mu0intfac
      DOUBLE PRECISION trff_tab(nradii,ng,2,abins,mu0bins)
      DOUBLE PRECISION cosne_tab(nradii,ng,2,abins,mu0bins)
      DOUBLE PRECISION gmin_tab(nradii,abins,mu0bins)
      DOUBLE PRECISION gmax_tab(nradii,abins,mu0bins)
      DOUBLE PRECISION re(nradii),gmin(nradii)
      DOUBLE PRECISION gmax(nradii),trff(nradii,ng,2)
      DOUBLE PRECISION cosne(nradii,ng,2)
      DOUBLE PRECISION ispec,lspec(ne),normspec(ne),eeo(0:ne)
      DOUBLE PRECISION lspecfine(4*ne),eeofine(4*ne)
      DOUBLE PRECISION intgmin,intgmax,inttf(ng,2),intmu(ng,2),rad
      DOUBLE PRECISION intgs,intfac
      DOUBLE PRECISION rms,marginal,gee,gstar2,trf,mu,eem1,eem2
      DOUBLE PRECISION rmin,rmax,alp1,alp2,rbreak,lineE,z
      DOUBLE PRECISION rmin_grid,rmax_grid
      DOUBLE PRECISION pi,r1,r2,r(nradii),wr(nradii),g1,g2,wg(ng)
c      
c   Variables specifically for fgsl library calls
c
      DOUBLE PRECISION r_i, wr_i
      type(fgsl_integration_glfixed_table) integr_table;
      integer(fgsl_size_t) ngl, igl
      integer fgslstat

      character(255) datafile, contxt

      LOGICAL qfirst, qanyf

      integer lenact
      character(255) fgmodf
      external lenact, fgmodf

      SAVE qfirst
      save a_tab,trff_tab,cosne_tab,gmin_tab,gmax_tab,mu0_tab

      DATA qfirst / .true. /

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

c
c     Get parameters of model from xspec
c
      lineE=param(1)
      alp1=param(2)
      alp2=param(3)
      rbreak=param(4)
      a=param(5)
      theta0=param(6)
      rmin=param(7)
      rmax=param(8)
      z=param(9)

      rms=marginal(a)

      rmin=rmin*rms
      rmax=rmax*rms

c
c     Set pi and convert angle variables
c
      pi = 4.d0*atan(1.d0)
      theta0=theta0*pi/180.d0
      mu0=cos(theta0)
c
c     Read in transfer function from the kerrtable.fits file.
c
      if ( qfirst ) then
         call getlun(ilun)
         datafile = fgmodf()
         datafile=datafile(1:lenact(datafile))//'/kerrtable.fits'
         ios = 0
         call ftopen(ilun,datafile,0,block,ios)
         if ( ios .NE. 0 ) THEN
            contxt = 'KERRDISK: Failed to open '
     &           //datafile(:lenact(datafile))
            call xwrite(contxt, 5)
            write(contxt, '(a,i4)') 'Status = ', ios
            call xwrite(contxt, 5)
            call frelun(ilun)
            return
         endif
         CALL ftmrhd(ilun, 1, hdutyp, ios)
         CALL ftgcvd(ilun, 1, 1, 1, abins, 0.0, a_tab, qanyf, ios)
         CALL ftmrhd(ilun, 1, hdutyp, ios)
         CALL ftgcvd(ilun, 1, 1, 1, mu0bins, 0.0, mu0_tab, qanyf, ios)
         CALL ftmrhd(ilun, 1, hdutyp, ios)
         do i=1,abins
            do j=1,mu0bins
               ir = (i-1)*abins + j
               call ftgcvd(ilun,1,ir,1,nradii,0.0,gmin,qanyf,ios)
               call ftgcvd(ilun,2,ir,1,nradii,0.0,gmax,qanyf,ios)
               call ftgcvd(ilun,3,ir,1,nradii*ng*2,0.0,trff,qanyf,ios)
               call ftgcvd(ilun,4,ir,1,nradii*ng*2,0.0,cosne,qanyf,ios)
               do ii=1,nradii
                  gmin_tab(ii,i,j)=gmin(ii)
                  gmax_tab(ii,i,j)=gmax(ii)
                  do jj=1,ng
                     do k=1,2
                        trff_tab(ii,jj,k,i,j)=trff(ii,jj,k)
                        cosne_tab(ii,jj,k,i,j)=cosne(ii,jj,k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         call ftclos(ilun,ios)
         call frelun(ilun)
         qfirst = .false.
      endif
c
c     Work out interpolation factors in spin and mu0 directions
c
      ia=1
      do i=1,abins-1
         if (a_tab(i) .lt. a) ia=i
      enddo
      aintfac=(a-a_tab(ia))/(a_tab(ia+1)-a_tab(ia))
c
      imu0=1
      do i=1,mu0bins
         if (mu0_tab(i) .lt. mu0) imu0=i
      enddo
      mu0intfac=(mu0-mu0_tab(imu0))/(mu0_tab(imu0+1)-mu0_tab(imu0))
c     
c     attempt to fix low-inclination problem
c
      if (mu0_tab(mu0bins) .lt. mu0) then
         imu0=mu0bins-1
         mu0intfac=1
      endif

c
c     Interpolate transfer function in a and mu0 plane
c
      do i=1,nradii
         gmin(i)=(1.0-aintfac)*(1.0-mu0intfac) * 
     $        gmin_tab(i,ia,imu0)
     $        +aintfac*(1.0-mu0intfac)*
     $        gmin_tab(i,ia+1,imu0)
     $        +(1.0-aintfac)*mu0intfac*
     $        gmin_tab(i,ia,imu0+1)
     $        +aintfac*mu0intfac*
     $        gmin_tab(i,ia+1,imu0+1)
c
         gmax(i)=(1.0-aintfac)*(1.0-mu0intfac) * 
     $        gmax_tab(i,ia,imu0)
     $        +aintfac*(1.0-mu0intfac)*
     $        gmax_tab(i,ia+1,imu0)
     $        +(1.0-aintfac)*mu0intfac*
     $        gmax_tab(i,ia,imu0+1)
     $        +aintfac*mu0intfac*
     $        gmax_tab(i,ia+1,imu0+1)
         do j=1,ng
            do k=1,2
               trff(i,j,k)=(1.0-aintfac)*(1.0-mu0intfac) * 
     $              trff_tab(i,j,k,ia,imu0)
     $              +aintfac*(1.0-mu0intfac)*
     $              trff_tab(i,j,k,ia+1,imu0)
     $              +(1.0-aintfac)*mu0intfac*
     $              trff_tab(i,j,k,ia,imu0+1)
     $              +aintfac*mu0intfac*
     $              trff_tab(i,j,k,ia+1,imu0+1)
c
               cosne(i,j,k)=(1.0-aintfac)*(1.0-mu0intfac) * 
     $              cosne_tab(i,j,k,ia,imu0)
     $              +aintfac*(1.0-mu0intfac)*
     $              cosne_tab(i,j,k,ia+1,imu0)
     $              +(1.0-aintfac)*mu0intfac*
     $              cosne_tab(i,j,k,ia,imu0+1)
     $              +aintfac*mu0intfac*
     $              cosne_tab(i,j,k,ia+1,imu0+1)
            enddo
         enddo
      enddo
c
c     Set up a radial grid: inversely spaced, and define integration values
c     for g* as in radial case 
c
      rmin_grid=rms
      rmax_grid=500*rms
      r1=1.d0/sqrt(rmax_grid)
      r2=1.d0/sqrt(rmin_grid)
      
      ngl = nradii
      integr_table = fgsl_integration_glfixed_table_alloc(ngl)      
      do i=1,nradii
         igl = i-1
         fgslstat = fgsl_integration_glfixed_point(r1,r2,igl,r_i,wr_i,
     &             integr_table)
         r(i) = r_i
         wr(i) = wr_i
         re(i)=1.d0/(r(i)**2.0)
      enddo
      call fgsl_integration_glfixed_table_free(integr_table)
      
      g1=0.d0
      g2=1.d0
      ngl = ng
      integr_table = fgsl_integration_glfixed_table_alloc(ngl)      
      do i=1,ng
         igl = i-1
         fgslstat = fgsl_integration_glfixed_point(g1,g2,igl,r_i,wr_i,
     &            integr_table)
         gstar(i) = r_i
         wg(i) = wr_i
      enddo
      call fgsl_integration_glfixed_table_free(integr_table)
      
c
c---------Integrate the line profile-------------------------
c

c
c     Generate energy grid and finer grid within it (4 x finer) 
c     to effectively get greater resolution than before without 
c      smoothing.  We will linearly interpolate between grid point 
c     values later 
c
      eeo(0)=dble(ear(0))
      do ii=1,ne
         lspec(ii)=0.d0
         eeo(ii)=dble(ear(ii))
      enddo
      do ii=1,ne
         do j=1,4
            lspecfine((ii-1)*4+j)=0.d0
            intfac=float(j)/4.0
            eeofine((ii-1)*4+j) = 
     &           intfac*eeo(ii)+(1.0-intfac)*eeo(ii-1)
         enddo
      enddo
c
c     Latest generation of line integration.  Integrates over a
c     large number of radii, using linear radial interpolation 
c     of the TF as well as gmin and gmax 
c
      irad=nradii-1
      do ii=1,nradii
         if (rmin .lt. re(ii)) irad=ii
      enddo
c
c     work out the number of zones for the final radial integral
c     ... currently set to ne per decade of radius
c
      ilgrad_max=3*int(float(ne)*dlog10(re(1)/re(nradii)))
c
c     Perform radial integral...
c
      do ilgrad=1,ilgrad_max
c
c        work out corresponding rrelevant 
c
         rad=rmin * 10 ** ( dfloat(ilgrad-1)*dlog(re(1)/re(nradii))/
     &     dfloat(ilgrad_max-1))

         if ((rad .gt. rmin) .and. (rad .lt. rmax)) then

         if (rad .gt. re(irad)) irad=irad-1
         intfac=(rad-re(irad+1))/(re(irad)-re(irad+1))
         do j=1,ng
            do k=1,2
               inttf(j,k)=intfac*trff(irad,j,k) + 
     &              (1.0-intfac)*trff(irad+1,j,k)
               intmu(j,k)=intfac*cosne(irad,j,k) + 
     &              (1.0-intfac)*cosne(irad+1,j,k)
            enddo
         enddo
         intgmin=intfac*gmin(irad)+(1.0-intfac)*gmin(irad+1)
         intgmax=intfac*gmax(irad)+(1.0-intfac)*gmax(irad+1)
c
c     EMISSIVITY PROFILE is hardwired in here.  Currently a broken
c     powerlaw
c
         if (rad .lt. rbreak) then
            ispec=(rad/rbreak)**(-alp1)
         else 
            ispec=(rad/rbreak)**(-alp2)
         endif
c
         eem1=lineE*intgmin/(1.0d0+z)
         eem2=lineE*intgmax/(1.0d0+z)
         ii1=1
         ii2=1

         do ii=1,4*ne
            if (eeofine(ii) .lt. eem1) ii1=ii
         enddo
         do ii=ii1-1,4*ne
            if (ii .gt. 0) then
               if (eeofine(ii) .lt. eem2) ii2=ii
            endif
         enddo

         trf = 0.0
         mu = 0.0
         if ((ii1 .gt. 1) .and. (ii2 .gt. 1)) then
            do ii=ii1+1,ii2
            gee=eeofine(ii)/(lineE/(1.0d0+z))
            gstar2=(gee-intgmin)/(intgmax-intgmin)
            do k=1,2
               if (gstar2 .le. gstar(1)) then
                  trf=inttf(1,k)
                  mu=intmu(1,k)
               endif
               if (gstar2 .ge. gstar(ng)) then
                  trf=inttf(ng,k)
                  mu=intmu(ng,k)
               endif
               if ((gstar2 .lt. gstar(ng)) .and. 
     &              (gstar2 .gt. gstar(1))) then
                  igstar2=1
                  do j=1,ng-1
                     if (gstar(j) .lt. gstar2) igstar2=j
                  enddo
                  intgs=(gstar2-gstar(igstar2))/
     &                 (gstar(igstar2+1)-gstar(igstar2))
                  trf=intgs*inttf(igstar2,k) + 
     &                 (1.0-intgs)*inttf(igstar2+1,k)
                  mu=intgs*intmu(igstar2,k) + 
     &                 (1.0-intgs)*intmu(igstar2+1,k)
               endif
c
c     Next line actually is the guts of the integral. The
c     LIMB DARKENING is hardwired in here.  It is currently set to
c     mu*(1+2.06*mu)
c
               lspecfine(ii)=lspecfine(ii) +
     &              rad*gee*(2.0*pi*gee)**2.0
     &              *trf*ispec*(1.0+2.06*mu)*rad/
     &              (sqrt(gstar2-gstar2**2.0)*(intgmax-intgmin))
            enddo
         enddo
         endif
         endif
      enddo
c
c     Bin up lspecfine to give to lspec 
c
      do i=1,ne
         do j=1,4
            lspec(i)=lspec(i)+lspecfine((i-1)*4+j)
         enddo
      enddo

c
c---------Output of spectral luminosity-----------------
c
c     Divide each lspec(ii) by the observed energy of that 
c     gridpoint --> ph/cm^2/s units 
c
      do ii=1,ne
         lspec(ii)=lspec(ii)/eeo(ii)
      enddo
      sumspec=0.d0
      do ii=1,ne
         if (lspec(ii) .gt. 0.d0) then
c
c     Sumspec weighted by the energy bin size 
c
            sumspec=sumspec+lspec(ii)*(eeo(ii)-eeo(ii-1))
         endif
      enddo
c      open(7,file='lspec.dat',status='unknown')
      do ii=1,ne
         if (lspec(ii) .gt. 0.d0) then
c
c     Normspec weighted by the energy bin size 
c
            normspec(ii)=lspec(ii)*(eeo(ii)-eeo(ii-1))/sumspec
         else 
            normspec(ii)=0.d0
         endif
         photar(ii)=SNGL(normspec(ii))
c         write(7,*) eeo(ii),normspec(ii)
      enddo
c      close(unit=7)

      end
c
c---------Subroutines------------------------
c
      function marginal(a)
      implicit none
      DOUBLE PRECISION marginal,a,Z1,Z2
      Z1=1.0+(1.0-a**2.0)**0.33*((1.0+a)**0.33+(1.0-a)**0.33)
      Z2=((3.0*a**2.0)+(Z1**2.0))**0.5
      marginal=3.0+Z2-((3.0-Z1)*(3.0+Z1+(2*Z2)))**0.5
      return

      end 


C========================================================================


