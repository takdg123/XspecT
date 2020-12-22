
      SUBROUTINE nsa(ear,ne,param,ifl,photar,photer)

      IMPLICIT NONE

      INTEGER ne, ifl
      REAL ear(0:ne),param(4),photar(ne),photer(ne)
      

c------------------------------------------------------------------------------
c
c      ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c      a SHORT DESCRIPTION:
c
c       Spectrum of X-ray radiation from neutron
c       star  atmosphere.
c       with account for the Comptonization effect
c       (see Zavlin et al. 1996, A&A, 315, 141 and
c       Pavlov et al. 1992 MNRAS, 253, 193) 
c       Models are available for three magnetic field strengths 0, 1e12, 1e13.
c       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c      INPUT PARAMETERS:
c
c   param(1) - Log of the effective (UNREDSHIFTED) temperature of the 
c              neutron star surface (in K);
c              Log T=5.5-7.0
c   param(2) - neutron star gravitational mass (in solar mass)
c   param(3) - neutron star radius (in km)
c   param(4) - magnetic field strength
c
c-----------------------------------------------------------------------------

      REAL temp(21), ene(1000)
      REAL t, rms, rs, gr, sa, t1, t2, dt, e
      REAL de, f1, f2, f, ff, magfld, magsve
    
      DOUBLE PRECISION flux(21,1000)
      DOUBLE PRECISION fluxin(21)

      INTEGER ilun, lfil, ios, ninp, minp
      INTEGER i, j, jt, kk, k, index, irow
      CHARACTER(255) filenm, contxt
      CHARACTER(128) pname
      CHARACTER(8) extnm, modnam(3)

      LOGICAL qanyf

      CHARACTER(255) fgmstr, fgmodf
      INTEGER lenact
      EXTERNAL fgmstr, lenact, fgmodf

      SAVE ninp, minp, temp, ene, flux, t1, t2, magsve

      DATA t1, magsve / 0., -1./

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO
      
      t = param(1)
      rms = param(2)
      rs = param(3)
      magfld = param(4)


      gr = sqrt(1.-2.952*rms/rs)
      sa = (rs/3.086e13)**2
      

c find out whether we need to read the data files

      if ( t1 .gt. 0. .AND. magsve .EQ. magfld ) go to 1

      magsve = magfld

      CALL getlun(ilun)

c Open the model file. First check for an NSA_FILE model string and if it is present 
c use that. If it is not then look for the appropriate file in the standard manager 
c directory.

      pname = 'NSA_FILE'
      filenm = fgmstr(pname)
      lfil = lenact(filenm)

      IF ( lfil .EQ. 0 ) THEN
         filenm = fgmodf()
         lfil = lenact(filenm)
         filenm = filenm(1:lfil)//'nsadata.fits'
      ENDIF

      IF ( ABS(magfld) .LT. 1e9 ) THEN
         extnm = 'LOWFIELD'
      ELSEIF ( ABS(magfld-1e12) .LT. 1e9 ) THEN
         extnm = '12GFIELD'
      ELSEIF ( ABS(magfld-1e13) .LT. 1e9 ) THEN
         extnm = '13GFIELD'
      ELSE
         CALL xwrite(
     &        'The magnetic field must be one of 0, 1e12, or 1e13 G', 5)
         RETURN
      ENDIF

      contxt = 'Using '//filenm(:lenact(filenm))
      CALL xwrite(contxt, 25)

      CALL FTOPEN(ilun,filenm,0,index,ios)
      IF ( ios .NE. 0 ) THEN
         contxt = 'NSA: Failed to open '//filenm(:lenact(filenm))
         CALL xwrite(contxt, 5)
         WRITE(contxt, '(a,i4)') 'Status = ', ios
         CALL xwrite(contxt, 5)
         CALL frelun(ilun)
         RETURN
      ENDIF

c Read the MODELINFO extension.

      CALL ftmnhd(ilun, 2, 'MODELINFO', 1, ios)
      CALL ftgcvs(ilun, 1, 1, 1, 3, " ", modnam, qanyf, ios)

c Find the row containing information about the magnetic field case in use

      irow = 0
      DO i = 1, 3
         IF ( extnm .EQ. modnam(i) ) irow = i
      ENDDO
      if ( irow .EQ. 0 ) then
         contxt = 'NSA: Failed to find extension '//extnm
         CALL xwrite(contxt, 5)
         CALL ftclos(ilun, ios)
         CALL frelun(ilun)
         RETURN
      ENDIF

      CALL ftgcvj(ilun, 3, irow, 1, 1, 0.0, minp, qanyf, ios)
      CALL ftgcve(ilun, 4, irow, 1, minp, 0.0, temp, qanyf, ios)
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

c Go to the required extension for the energies and fluxes

      CALL FTMNHD(ilun, 2, extnm, 1, ios)
      CALL FTGNRW(ilun, ninp, ios)
      CALL FTGCVE(ilun, 1, 1, 1, ninp, 0.0, ene, qanyf, ios)

      do i=1,ninp
         ene(i)=alog10(ene(i))
         CALL FTGCVD(ilun, 2, i, 1, minp, 0.0, fluxin, qanyf, ios)
         do j=1,minp
            if(fluxin(j).gt.0.0) then
               flux(j,i)=log10(fluxin(j))
            else
               flux(j,i)=flux(j,i-1)
            endif
         enddo
      enddo
      CALL ftclos(ilun, ios)
      CALL frelun(ilun)

      t1=temp(1)
      t2=temp(minp)

c jump to here if we did not need to read a data file

1     continue
      
      do jt=2,minp
         if(temp(jt).ge.t) go to 2
      enddo
      jt=minp
2     dt=(t-temp(jt-1))/(temp(jt)-temp(jt-1))

      kk=2

      do i=0,ne
         e=alog10(ear(i)/gr)
         if(e.lt.ene(1)) e=ene(1)
         if(e.gt.ene(ninp)) go to 4

         do k=kk,ninp
            if(ene(k).ge.e) go to 3
         enddo

3        de=(e-ene(k-1))/(ene(k)-ene(k-1))
         f1=REAL(flux(jt-1,k-1))+de*(REAL(flux(jt-1,k)-flux(jt-1,k-1)))
         f2=REAL(flux(jt,k-1))+de*(REAL(flux(jt,k)-flux(jt,k-1)))
         f=f1+dt*(f2-f1)
         f=10**f*sa
         go to 5
4        photar(i)=photar(i-1)
         go to 6
5        if(i.eq.0) go to 7      
         photar(i)=(f+ff)/2.*(ear(i)-ear(i-1))
7        ff=f
         kk=k
6     enddo

      RETURN
      end

