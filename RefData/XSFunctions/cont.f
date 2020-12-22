*
*        Last modified Aug 23, 1999 (kjb) - subroutine hunt added.
*
      integer function ibine(e,ear,ne)
c   returns # of bin which contains energy e
      integer ne
      real e,ear(*)
      
      do ibine=0,ne-1
         if (e .lt. ear(ibine+1))return
      enddo
      ibine=0
      return
      end
*
*
*
      subroutine zrebin(specin,specout,numoutbins,numcmbins,
     $                  indexin,wl,wr,zstart,ein,area)

c  calculates area under continuum for output bins
c   requires znibin be called first to initialize weight and index arrays

      real ein(*), area(0:*)
      real specin(*),specout(*),wl(*),wr(*)
      INTEGER indexin(*),zstart(*)

      INTEGER icmb,iout,numcmbins,numoutbins,imach

      real small, r1mach
      real s1,s2,wrr,wll
      real log
      
      imach = 1
      small = r1mach(imach)
      area(0)=0.0
c   first calculate the area in each interval
      do icmb=1,numcmbins
*            Interpolation for exponential function.
         s1=specin(indexin(icmb))
         s2=specin(indexin(icmb)+1)
         wrr=wr(icmb)
         wll=wl(icmb)
*           Should be rewritten for the case of s2=s1
         if( s2 .gt. small .and. log(s2/s1) .ne. 0.0 )then
         area(icmb)= ( s1**(1.-wrr)*s2**wrr-s1**(1.-wll)*s2**wll )/
     .    log(s2/s1)*(ein(indexin(icmb)+1)-ein(indexin(icmb)))
         else
            area(icmb)=0.0
         endif
*         area(icmb)=wl(icmb)*specin(indexin(icmb))+
*     $        wr(icmb)*specin(indexin(icmb)+1)
      enddo
c   now sum up the intervals appropriate for each output bin
      do iout=1,numoutbins
         specout(iout)=0.0
         do icmb=zstart(iout), zstart(iout+1)-1
            specout(iout)=specout(iout)+area(icmb)
         enddo
      enddo

      return
      end

      subroutine znibin(numinbins,ein,numoutbins,eout,z,
     $     numcmbins,indexin,zstart,wl,wr,ecmb)

      real ein(*),eout(*),wl(*),wr(*),ecmb(*)
      INTEGER indexin(*),zstart(*)
      INTEGER numinbins,numoutbins, numcmbins

c   sets up arrays needed to calculate counts from continua in model bins
c    note that this assumes that the inputs are the continuum heights
c    sampled at energies, ein, (rather than integrated over ein bins)
c Parameters
c    numinbins       i        i: Number of bins in the input array
c    ein             r        i: Energies for the input array
c    numoutbins      i        i: Number of bins in the output array
c    eout            r        i: Energies for the input array
c    z               r        i: Redshift
c    numcmbins       i        r: Number of combined bins
c    indexin         i        r: input bin for this combined bin
c    zstart          i        r: start combined bin for this output bin
c    wl              r        r: weight for low end of this combined bin
c    wr              r        r: weight for high end of this combined bin
c    ecmb            r        r: start energy for this combined bin


      INTEGER istopout,lastin
      INTEGER iin,iout,icmb
      REAL z,zf
      REAL elresp, ehresp, elmod, ehmod

      INTEGER type_last

      LOGICAL qlow, qhigh

      CHARACTER(72) contxt

      SAVE qlow, qhigh
      SAVE elresp, ehresp, elmod, ehmod

      DATA qlow, qhigh /.FALSE., .FALSE./


      if (z .ne. -1.0) then
         zf=1.0/(1.0+z)
      else
         call xwrite(' *ERROR: ZNIBIN: Infinite negative redshift ', 2)
         zf=0.0
      endif

      icmb=1
      iin=1
      iout=1
      lastin=1
c   first, point iin at first relevant input energy
      if (zf*ein(1) .gt. eout(numoutbins)) then
         call xwrite(
     $    ' *ERROR*: ZNIBIN: input range all above output range',
     $        5)
         return
      endif

      if (eout(1) .lt. zf*ein(1)) goto 1100
 1000 continue
      if (eout(1) .lt. zf*ein(iin+1)) goto 1200
      iin=iin+1
      if (iin .gt. numinbins) then
         call xwrite(
     $    ' *ERROR*: ZNIBIN: input range all below output range',
     $        5)
         return
      endif
      goto 1000

c comes here if output range starts before beginning of input range
 1100 continue

      IF ( .NOT.qlow .OR. ABS(elresp-eout(1)) .GT. 1.e-6 .OR.
     &     ABS(elmod-zf*ein(1)) .GT. 1.e-6 ) THEN
         call xwrite(
     $  ' *WARNING*: ZNIBIN: some low energy bins before stored model',
     $     2)
         WRITE(contxt,'(a,1pg9.2,a,1pg9.2)') ' Response starts at ', 
     *     eout(1), ' and stored model at ', zf*ein(1)
         CALL xwrite(contxt, 2)
         qlow = .TRUE.
         elresp = eout(1)
         elmod = zf*ein(1)
      ENDIF
 1110 continue
      if (eout(iout) .ge. zf*ein(1)) goto 1200
      zstart(iout)=0
      iout=iout+1
      goto 1110

c  starting places fixed up - now work on ending
 1200 continue
      istopout=numoutbins
      if (eout(istopout) .lt. zf*ein(numinbins)) goto 2000
      IF ( .NOT.qhigh .OR. ABS(ehresp-eout(istopout)) .GT. 1.e-6 .OR.
     &     ABS(ehmod-zf*ein(numinbins)) .GT. 1.e-6 ) THEN
         call xwrite(
     $  ' *WARNING*: ZNIBIN: some high energy bins after stored model',
     $     2)
         WRITE(contxt,'(a,1pg9.2,a,1pg9.2)') ' Response ends at ', 
     *     eout(istopout), ' and stored model at ', zf*ein(numinbins)
         CALL xwrite(contxt, 2)
         qhigh = .TRUE.
      ENDIF
 1210 continue
      if (eout(istopout) .lt. zf*ein(numinbins)) goto 2000
      zstart(istopout)=0
      istopout=istopout-1
      goto 1210


c  ends of ranges fixed up, now generate combined info
 2000 continue
 2010 continue
      type_last=1
      if (zf*ein(iin) .LT. eout(iout) ) GOTO 2100
      if (zf*ein(iin) .EQ. eout(iout) ) GOTO 2200
      if (zf*ein(iin) .GT. eout(iout) ) GOTO 2300

c   eout > ein - next combined bin is ein
 2100 continue
      ecmb(icmb)=zf*ein(iin)
c     type(icmb)=1
      type_last=1
      lastin=iin
      indexin(icmb)=lastin
      icmb=icmb+1
      iin=iin+1
      goto 2010

c   eout = ein - next combined bin is either - advance both
 2200 continue
      ecmb(icmb)=zf*ein(iin)
c     type(icmb)=1
      if (type_last .eq. 1) zstart(iout)=icmb
      type_last=2
      lastin=iin
      indexin(icmb)=lastin
      icmb=icmb+1
      iin=iin+1
      iout=iout+1
      if (iout .gt. istopout) goto 3000
      goto 2010

c   ein > eout - next combined bin is eout
 2300 continue
      ecmb(icmb)=eout(iout)
      if (type_last .eq. 1) zstart(iout)=icmb
c     type(icmb)=2
      type_last=2
      indexin(icmb)=lastin
      icmb=icmb+1
      iout=iout+1
      if (iout .gt. istopout) goto 3000
      goto 2010

 3000 continue
      numcmbins=icmb-2
c   now calculate the weighting factors for the endpoints
      do icmb=1,numcmbins
         call wghts(zf*ein(indexin(icmb)),
     $        ecmb(icmb),
     $        ecmb(icmb+1),
     $        zf*ein(indexin(icmb)+1),
     $        wl(icmb),wr(icmb))
      enddo
      return
      end

      subroutine wghts(e1,el,er,e2,wl,wr)
c   calculates weights to multiply continuum samples (at e1 and e2)
c    in order to determine area under continuum between el and er

      real e1,el,er,e2,wl,wr

*      wr=(0.5/(e2-e1))*(er**2-2.0*e1*(er-el)-el**2)
*      wl=(er-el)-wr
*        Weights for exponential interpolation.
      wr=(er-e1)/(e2-e1)
      wl=(el-e1)/(e2-e1)
      return
      end
