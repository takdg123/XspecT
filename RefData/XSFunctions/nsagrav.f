      subroutine nsagrav(ear,ne,param,ifl,photar,photer)
      implicit real (a-h,o-z)
      implicit integer (i-n)

c------------------------------------------------------------------------------
c
c      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c      a SHORT DESCRIPTION:
c
c       Spectrum of X-ray radiation from neutron
c       star  atmosphere.
c       CURRENTLY ONLY NONMAGNETIC HYDROGEN MODELS AVAILABLE
c       with account for the Comptonization effect
c       (see Zavlin et al. 1996, A&A, 315, 141 and
c       Pavlov et al. 1992 MNRAS, 253, 193) 
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c      INPUT PARAMETERS:
c
c   param(1) - Log of the effective (UNREDSHIFTED) temperature of the 
c              neutron star surface (in K);
c              Log T=5.5-6.5
c   param(2) - neutron star gravitational mass (in solar mass)
c   param(3) - neutron star radius (in km)
c
c-----------------------------------------------------------------------------

      real ear(0:ne),param(3),photar(ne),photer(ne)
c

      INTEGER NINP, MINP, NMOD
      PARAMETER(NINP=1000, MINP=11, NMOD=5)

      real temp(MINP),ene(NINP),flux1(MINP,NINP),flux2(MINP,NINP)
      real flux(MINP,NINP,NMOD), fluxin(MINP)
      real t1, t2
      
      INTEGER ilun, lenn, ios, i, index
      CHARACTER(255) nsadir, filenm, contxt, pname
      character(9) modnam(NMOD)

      LOGICAL qanyf

      INTEGER lenact
      CHARACTER(255) fgmodf, fgmstr
      EXTERNAL lenact, fgmodf, fgmstr
      
      save t1, t2, temp, ene, flux

      DATA t1 / 0. /

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      t=param(1)
      rms=param(2)
      rs=param(3)

      rgr=1.-2.952*rms/rs
      gr=sqrt(1./3.)
      if(rgr.ge.(1./3.)) gr=sqrt(rgr)

      gs=1.33e2*rms/rs/rs/gr
      gsl=alog10(gs)
      if(gsl.lt.-1e0) gsl=-1e0
      if(gsl.gt.1e0) gsl=1e0

      sa=(rs/3.086e13)**2

c If a temperature has been set then we have already read in the model
c data files so jump straight to the model calculation      

      IF ( t1 .EQ. 0.) THEN

c Find the directory for the model data files. First check for an 
c NSAGRAV_DIR model string and if it is present use that. If it is not 
c then look for the file in the standard directory.

         pname = 'NSAGRAV_DIR'
         nsadir = fgmstr(pname)
         lenn = lenact(nsadir)
         IF ( lenn .EQ. 0 ) nsadir = fgmodf()
         lenn = lenact(nsadir)

c Open the file

         filenm = nsadir(:lenn)//'nsagravdata.fits'
         CALL getlun(ilun)
         ios = 0
         CALL ftopen(ilun, filenm, 0, index, ios)
         IF ( ios .NE. 0 ) THEN
            contxt = 'NSAGRAV: Failed to open nsagrav.fits'
            CALL xwrite(contxt, 5)
            WRITE(contxt, '(a,i4)') 'Status = ', ios
            CALL xwrite(contxt, 5)
            CALL frelun(ilun)
            RETURN
         ENDIF

c Read the MODELINFO extension. Note that the temperatures are assumed
c to be the same for every model so we just read the first row

         CALL ftmnhd(ilun, 2, 'MODELINFO', 1, ios)
         CALL ftgcvs(ilun, 1, 1, 1, NMOD, " ", modnam, qanyf, ios)
         CALL ftgcve(ilun, 3, 1, 1, MINP, 0.0, temp, qanyf, ios)
         IF ( ios .NE. 0 ) THEN
            contxt = 'NSAGRAV: Failed to read MODELINFO extension'
            CALL xwrite(contxt, 5)
            WRITE(contxt, '(a,i4)') 'Status = ', ios
            CALL xwrite(contxt, 5)
            ios = 0
            CALL ftclos(ilun, ios)
            CALL frelun(ilun)
            RETURN
         ENDIF

c Read the other extensions to get the fluxes

         do i=1,NMOD

            CALL ftmnhd(ilun, 2, modnam(i), 1, ios)
            CALL ftgcve(ilun, 1, 1, 1, NINP, 0.0, ene, qanyf, ios)

            do k=1,NINP
               CALL ftgcve(ilun, 2, k, 1, MINP, 0.0, fluxin, qanyf, ios)
               do j=1,MINP
                  if(fluxin(j).gt.0.0) then
                     flux(j,k,i)=alog10(fluxin(j))
                  else
                     flux(j,k,i)=flux(j,k-1,i)
                  endif
               enddo
               ene(k)=alog10(ene(k))
            enddo

         enddo

         CALL ftclos(ilun, ios)
         CALL frelun(ilun)

         t1=temp(1)
         t2=temp(MINP)
         if(t.lt.t1) t=t1
         if(t.gt.t2) t=t2

      ENDIF

      if(gsl.le.-0.5e0) then
         gsl1=-1e0
         gsl2=-0.5e0
         ig=1
      endif
      if(gsl.le.0.0e0.and.gsl.gt.-0.5e0) then
         gsl1=-0.5e0
         gsl2=0e0
         ig=2
      endif
      if(gsl.le.0.5e0.and.gsl.gt.0e0) then
         gsl1=0e0
         gsl2=0.5e0
         ig=3
      endif
      if(gsl.le.1.0e0.and.gsl.gt.0.5e0) then
         gsl1=0.5e0
         gsl2=1e0
         ig=4
      endif
      dg=(gsl-gsl1)/(gsl2-gsl1)


      do i=1,NINP
         do j=1,MINP
            flux1(j,i)=flux(j,i,ig)
            flux2(j,i)=flux(j,i,ig+1)
         enddo
      enddo

      do jt=2,MINP
         if(temp(jt).ge.t) go to 2
      enddo
      jt=MINP
 2    dt=(t-temp(jt-1))/(temp(jt)-temp(jt-1))

      kk=2

      do i=0,ne
         e=alog10(ear(i)/gr)
         if(e.lt.ene(1)) e=ene(1)
         if(e.gt.ene(NINP)) go to 4
         
         do k=kk,NINP
            if(ene(k).ge.e) go to 3
         enddo

 3       de=(e-ene(k-1))/(ene(k)-ene(k-1))
         f1=flux1(jt-1,k-1)+de*(flux1(jt-1,k)-flux1(jt-1,k-1))
         f2=flux1(jt,k-1)+de*(flux1(jt,k)-flux1(jt,k-1))
         f3=f1+dt*(f2-f1)
         f4=flux2(jt-1,k-1)+de*(flux2(jt-1,k)-flux2(jt-1,k-1))
         f5=flux2(jt,k-1)+de*(flux2(jt,k)-flux2(jt,k-1))
         f6=f4+dt*(f5-f4)
         f=f3+dg*(f6-f3)
         f=sa*10**f
         go to 5
 4       photar(i)=photar(i-1)
         go to 6
 5       if(i.eq.0) go to 7      
         photar(i)=(f+ff)/2.*(ear(i)-ear(i-1))
 7       ff=f
         kk=k
 6    enddo


      end
