      subroutine felines(t,fei,elfeo,xlinefeo,nfe)
*
*        Error with array elfeo fixed on Sep 28, 99.
*        Last modification - Aug 23, 1999 (kjb) - line sorting added.
*        Fe L-shell lines.
*        Input:  t        -  temperature (K)
*                fei      -  Fe ionic concentrations.
*        Output: elfeo    -  line energies in keV
*                xlinefeo -  line intensities 
*                nfe      -  number of lines
*
**************  Ne- to Li-like species ***************
*
      real keV
      parameter ( keV = 1.60219e-9 )
      real t, fei(*), elfeo(*), xlinefeo(*)
      integer nmaxfe
      parameter ( nmaxfe = 146000 )
      real elfe(nmaxfe), xlinefe(nmaxfe)
      integer index(nmaxfe)
      integer nionmn, nionmx, ion, ind
      integer istart, iend, iline
      parameter ( nionmn = 17, nionmx = 24 )
      integer nl(nionmn:nionmx)
      integer nfe
      logical qini
      save qini, index
      data qini /.true./
      if( qini )then
         call feread
      endif
      ion = nionmx
      call felin(ion,t,elfe,xlinefe,nl(ion))
      istart = 1
      iend = nl(ion)
      do 200 iline=istart,iend
         xlinefe(iline) = fei(ion)*xlinefe(iline)
 200  continue
      istart = iend + 1
      ind = 1
      do 220 ion = nionmx-1,nionmn,-1
         ind = ind + nl(ion+1)
         call felin(ion,t,elfe(ind),xlinefe(ind),nl(ion))
         iend = iend + nl(ion)
         do 240 iline=istart,iend
            xlinefe(iline) = fei(ion)*xlinefe(iline)
 240     continue
         istart = iend + 1
 220  continue
      nfe = iend
      if( qini )then
         call sortrx(nfe,elfe,index)
         qini = .false.
      endif
      do 310 iline=1,nfe
         elfeo(iline) = elfe(index(iline))
         xlinefeo(iline) = xlinefe(index(iline))
 310  continue
      return
      end
*
*
*
      subroutine feread
*
*        Reads Fe atomic data.
*                weight0 --   statistical weight of the ground level.
*                k       --   collisional transition type (0 - no data,
*                             1 - allowed transitions, 2 - forbidden
*                             transitions). All transitions are from the
*                             ground state to higher levels, ordered with
*                             increasing energy.
*                eij     --   excitation energy in Rydbergs.
*                c       --   scaling constant for spline fits to collision
*                             strengths data.
*                p1m's   --   collision strength spline coefficients.
*                nl      --   number of levels.
*                ntrpnti --   lower level of radiative transition.
*                ntrpntj --   upper level of radiative transition.
*                einsta  --   Einstein A transition value.
*                eln    --   line energy in keV.
*                ntr     --   number of transitions.
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      parameter ( nionmn = 17, nionmx = 24 )
      parameter ( ntrmx = 260 + 12034 + 17829 + 15719 + 41967 + 36029 +
     .                    13095 + 6999 )
      parameter ( nlmx = 42 + 302 + 357 + 848 + 1087 + 994 + 501 + 293 )
      dimension eln(ntrmx),ntrpnti(ntrmx),ntrpntj(ntrmx),einsta(ntrmx)
      dimension nl(nionmn:nionmx), ntr(nionmn:nionmx)
      dimension weight0(nionmn:nionmx)
      dimension k(nlmx),eij(nlmx),c(nlmx),
     .  p1m(nlmx),p2m(nlmx),p3m(nlmx),p4m(nlmx),p5m(nlmx)
      dimension ib(nlmx), iwork(126,nlmx), suma(nlmx)
      common /combr/ ib,iwork,suma
      common /comlev/ nl,k,eij,c,p1m,p2m,p3m,p4m,p5m,weight0
      common /comtr/ ntr,ntrpnti,ntrpntj,eln,einsta
c      character(1) strdum
c      character(18) strdum1
      character(255) pathd, filenm, contxt
      logical qanyf
      dimension rlambd(ntrmx)
      equivalence (eln,rlambd)
      character(255) fgmodf
      external lenact, fgmodf

c      data nl /293, 501, 994, 1087, 848, 357, 302, 42/
c      data ntr /6999, 13095, 36029, 41967, 15719, 17829, 12034, 260/
c initialization - to use a data statement here would require a BLOCK DATA module

      nl(nionmn)   = 293
      nl(nionmn+1) = 501
      nl(nionmn+2) = 994
      nl(nionmn+3) = 1087
      nl(nionmn+4) = 848
      nl(nionmn+5) = 357
      nl(nionmn+6) = 302
      nl(nionmn+7) = 42

      ntr(nionmn)   = 6999
      ntr(nionmn+1) = 13095
      ntr(nionmn+2) = 36029
      ntr(nionmn+3) = 41967
      ntr(nionmn+4) = 15719
      ntr(nionmn+5) = 17829
      ntr(nionmn+6) = 12034
      ntr(nionmn+7) = 260


      pathd = fgmodf()
      lenn = lenact(pathd)
      ierr = 0

c read the outfe.fits file - each ionization stage is in its own extension

      filenm = pathd(:lenn) // 'outfe.fits'
      CALL getlun(lun)
      CALL ftopen(lun, filenm, 0, iblock, ierr)
      contxt = 'Failed to open '// filenm(:lenact(filenm))
      IF ( ierr .NE. 0 ) GOTO 999

      nbl = 0
      DO ion = 1, 8

         CALL ftmrhd(lun, 1, hdutyp, ierr)
         contxt = 'Failed to move to next HDU in '// 
     &            filenm(:lenact(filenm))
         IF ( ierr .NE. 0 ) GOTO 999

         CALL ftgcvj(lun, 3, 1, 1, nl(nionmx-ion+1)-1, 0, k(2+nbl), 
     &               qanyf, ierr)
         CALL ftgcve(lun, 4, 1, 1, nl(nionmx-ion+1)-1, 0.0, eij(2+nbl), 
     &               qanyf, ierr)
         CALL ftgcve(lun, 5, 1, 1, nl(nionmx-ion+1)-1, 0.0, c(2+nbl), 
     &               qanyf, ierr)
         CALL ftgcve(lun, 6, 1, 1, nl(nionmx-ion+1)-1, 0.0, p1m(2+nbl), 
     &               qanyf, ierr)
         CALL ftgcve(lun, 7, 1, 1, nl(nionmx-ion+1)-1, 0.0, p2m(2+nbl), 
     &               qanyf, ierr)
         CALL ftgcve(lun, 8, 1, 1, nl(nionmx-ion+1)-1, 0.0, p3m(2+nbl), 
     &               qanyf, ierr)
         CALL ftgcve(lun, 9, 1, 1, nl(nionmx-ion+1)-1, 0.0, p4m(2+nbl), 
     &               qanyf, ierr)
         CALL ftgcve(lun, 10, 1, 1, nl(nionmx-ion+1)-1, 0.0, p5m(2+nbl), 
     &               qanyf, ierr)

         contxt = 'Failed to read data from '//filenm(:lenact(filenm))
         IF ( ierr .NE. 0 ) GOTO 999

         nbl = nbl + nl(nionmx-ion+1)

      ENDDO

      CALL ftclos(lun, ierr)
      CALL frelun(lun)

c now repeat for the ferad.fits file

      filenm = pathd(:lenn) // 'ferad.fits'
      CALL getlun(lun)
      CALL ftopen(lun, filenm, 0, iblock, ierr)
      contxt = 'Failed to open '// filenm(:lenact(filenm))
      IF ( ierr .NE. 0 ) GOTO 999

      nbt = 0
      DO ion = 1, 8

         CALL ftmrhd(lun, 1, hdutyp, ierr)
         contxt = 'Failed to move to next HDU in '// 
     &            filenm(:lenact(filenm))
         IF ( ierr .NE. 0 ) GOTO 999

         CALL ftgcvj(lun, 3, 1, 1, ntr(nionmx-ion+1), 0, 
     &               ntrpnti(1+nbt), qanyf, ierr)
         CALL ftgcvj(lun, 5, 1, 1, ntr(nionmx-ion+1), 0, 
     &               ntrpntj(1+nbt), qanyf, ierr)
         CALL ftgcve(lun, 7, 1, 1, ntr(nionmx-ion+1), 0.0, 
     &               rlambd(1+nbt), qanyf, ierr)
         CALL ftgcve(lun, 9, 1, 1, ntr(nionmx-ion+1), 0.0, 
     &               einsta(1+nbt), qanyf, ierr)

         contxt = 'Failed to read data from '//filenm(:lenact(filenm))
         IF ( ierr .NE. 0 ) GOTO 999

         nbt = nbt + ntr(nionmx-ion+1)

      ENDDO

      CALL ftclos(lun, ierr)
      CALL frelun(lun)




cc      open(67,file=pathd(:lenn) // 'outlifem',status='old')
cc      nl(nionmx) = 42
cc      read(67,915) ( il,jl,k(j),eij(j),c(j),
cc     .         p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j=2,nl(24) )
cc 915  format(3I4,F13.6,F10.5,5E12.5)
cc      close(67)
cc      open(68,file=pathd(:lenn) // 'lifexxiv.rad', status='old')
cc      ntr(nionmx) = 260
cc      read(68,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     .       oscstr,rlambd(i),strdum1,einsta(i),i=1,ntr(nionmx) )
cc 900  format(A1,I4,I5,I3,I5,E16.6,E15.6,A18,E12.5)
cc      close(68)
cc      open(69,file=pathd(:lenn) // 'outbefem',status='old')
cc      nl(nionmx-1) = 302
cc      nbl=nl(nionmx)
cc      read(69,915) ( il,jl,k(j),eij(j),c(j),
cc     . p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j = 2+nbl,nl(nionmx-1)+nbl )
cc      close(69)
cc      open(70,file=pathd(:lenn) // 'befexxiii.rad', status='old')
cc      ntr(nionmx-1) = 12034
cc      nbt=ntr(nionmx)
cc      read(70,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     . oscstr,rlambd(i),strdum1,einsta(i),i=1+nbt,ntr(nionmx-1)+nbt )
cc      close(70)
cc      open(71,file=pathd(:lenn) // 'outbfem',status='old')
cc      nl(nionmx-2) = 357
cc      nbl=nbl+nl(nionmx-1)
cc      read(71,915) ( il,jl,k(j),eij(j),c(j),
cc     . p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j = 2+nbl,nl(nionmx-2)+nbl )
cc      close(71)
cc      open(72,file=pathd(:lenn) // 'bfexxii.rad', status='old')
cc      ntr(nionmx-2) = 17829
cc      nbt=nbt+ntr(nionmx-1)
cc      read(72,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     . oscstr,rlambd(i),strdum1,einsta(i),i=1+nbt,ntr(nionmx-2)+nbt )
cc      close(72)
cc      open(73,file=pathd(:lenn) // 'outcfem',status='old')
cc      nl(nionmx-3) = 848
cc      nbl=nbl+nl(nionmx-2)
cc      read(73,915) ( il,jl,k(j),eij(j),c(j),
cc     . p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j = 2+nbl,nl(nionmx-3)+nbl )
cc      close(73)
cc      open(74,file=pathd(:lenn) // 'cfexxi.rad', status='old')
cc      ntr(nionmx-3) = 15719
cc      nbt=nbt+ntr(nionmx-2)
cc      read(74,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     . oscstr,rlambd(i),strdum1,einsta(i),i=1+nbt,ntr(nionmx-3)+nbt )
cc      close(74)
cc      open(75,file=pathd(:lenn) // 'outnfem',status='old')
cc      nl(nionmx-4) = 1087
cc      nbl=nbl+nl(nionmx-3)
cc      read(75,915) ( il,jl,k(j),eij(j),c(j),
cc     . p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j = 2+nbl,nl(nionmx-4)+nbl )
cc      close(75)
cc      open(76,file=pathd(:lenn) // 'nfexx.rad', status='old')
cc      ntr(nionmx-4) = 41967
cc      nbt=nbt+ntr(nionmx-3)
cc      read(76,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     . oscstr,rlambd(i),strdum1,einsta(i),i=1+nbt,ntr(nionmx-4)+nbt )
cc      close(76)
cc      open(77,file=pathd(:lenn) // 'outofem',status='old')
cc      nl(nionmx-5) = 994
cc      nbl=nbl+nl(nionmx-4)
cc      read(77,915) ( il,jl,k(j),eij(j),c(j),
cc     . p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j = 2+nbl,nl(nionmx-5)+nbl )
cc      close(77)
cc      open(78,file=pathd(:lenn) // 'ofexix.rad', status='old')
cc      ntr(nionmx-5) = 36029
cc      nbt=nbt+ntr(nionmx-4)
cc      read(78,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     . oscstr,rlambd(i),strdum1,einsta(i),i=1+nbt,ntr(nionmx-5)+nbt )
cc      close(78)
cc      open(79,file=pathd(:lenn) // 'outffem',status='old')
cc      nl(nionmx-6) = 501
cc      nbl=nbl+nl(nionmx-5)
cc      read(79,915) ( il,jl,k(j),eij(j),c(j),
cc     . p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j = 2+nbl,nl(nionmx-6)+nbl )
cc      close(79)
cc      open(80,file=pathd(:lenn) // 'ffexviii.rad', status='old')
cc      ntr(nionmx-6) = 13095
cc      nbt=nbt+ntr(nionmx-5)
cc      read(80,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     . oscstr,rlambd(i),strdum1,einsta(i),i=1+nbt,ntr(nionmx-6)+nbt )
cc      close(80)
cc      open(81,file=pathd(:lenn) // 'outnefem', status='old')
cc      nl(nionmx-7) = 293
cc      nbl=nbl+nl(nionmx-6)
cc      read(81,915) ( il,jl,k(j),eij(j),c(j),
cc     . p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),j = 2+nbl,nl(nionmx-7)+nbl )
cc      close(81)
cc      open(82,file=pathd(:lenn) // 'nefexvii.rad', status='old')
cc      ntr(nionmx-7) = 6999
cc      nbt=nbt+ntr(nionmx-6)
cc      read(82,900) ( strdum,ntrdum1,ntrpnti(i),ntrdum2,ntrpntj(i),
cc     . oscstr,rlambd(i),strdum1,einsta(i),i=1+nbt,ntr(nionmx-7)+nbt )
cc      close(82)


      weight0(nionmx) = 2.0
      weight0(nionmx-1) = 1.0
      weight0(nionmx-2) = 2.0
      weight0(nionmx-3) = 1.0
      weight0(nionmx-4) = 4.0
      weight0(nionmx-5) = 5.0
      weight0(nionmx-6) = 4.0
      weight0(nionmx-7) = 1.0
*        Calculate line energies in keV.
      ntot = 0
      do 10 i = nionmn, nionmx
         ntot = ntot + ntr(i)
 10   continue
      do 60 i=1,ntot
         eln(i) = 12.3986/rlambd(i)
 60   continue
*
      nbl = 0
      nbt = 0
      do 40 i = nionmx, nionmn, -1
         do 20 j = 1,nl(i)
            ib(j+nbl) = 0
            suma(j+nbl) = 0.0
            do 30 l=1+nbt,ntr(i)+nbt
               if( ntrpntj(l) .eq. j )then
                  ib(j+nbl) = ib(j+nbl) + 1
                  iwork(ib(j+nbl),j+nbl) = l-nbt
                  suma(j+nbl) = suma(j+nbl) + einsta(l)
               endif
 30         continue
 20      continue
         nbl = nbl + nl(i)
         nbt = nbt + ntr(i)
 40   continue


 999  CONTINUE
      IF ( ierr .NE. 0 ) THEN
         CALL xwrite(contxt, 5)
         WRITE(contxt, '(a, i5)' ) 'FITSIO status = ', ierr
         CALL xwrite(contxt, 5)
      ENDIF


      return
      end
*
*
*
      subroutine felin(ion,temp,eline,flx,nline)
*
*        Calculates line intensities for Fe ions.
*        Collisional excitation from ground state only, recombination not
*        included, coronal approximation (no collisional mixing between
*        excited states).
*        Input:  ion     --   Fe ion (23 for Be-like Fe XXIII, 24 for
*                             Li-like Fe XXIV, etc)
*                temp    --   temperature (K).
*        Output: eline   --   line energies (keV).
*                flx     --   line emissivities (ergs/s per ion per electron).
*                             Must multiply by electron and ion densities to
*                             get fluxes in ergs/cm**3/s.
*                nline   --   number of lines.
*                
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      parameter ( ntrmxl = 41967 )
      parameter ( nionmn = 17, nionmx = 24 )
      parameter ( ntrmx = 260 + 12034 + 17829 + 15719 + 41967 + 36029 +
     .                    13095 + 6999 )
      parameter ( nlmx = 42 + 302 + 357 + 848 + 1087 + 994 + 501 + 293 )
      dimension eln(ntrmx),ntrpnti(ntrmx),ntrpntj(ntrmx),einsta(ntrmx)
      dimension nl(nionmn:nionmx), ntr(nionmn:nionmx)
      dimension weight0(nionmn:nionmx)
      dimension k(nlmx),eij(nlmx),c(nlmx),
     .  p1m(nlmx),p2m(nlmx),p3m(nlmx),p4m(nlmx),p5m(nlmx)
      dimension ib(nlmx), iwork(126,nlmx), suma(nlmx)
      common /combr/ ib,iwork,suma
      common /comlev/ nl,k,eij,c,p1m,p2m,p3m,p4m,p5m,weight0
      common /comtr/ ntr,ntrpnti,ntrpntj,eln,einsta
      dimension eline(*), flx(*)
      dimension br(ntrmxl), rates(ntrmxl)
*        Calculate direct excitation rates.
      nbl = 0
      do 5 i = ion, nionmx-1
         nbl = nbl + nl(i+1)
 5    continue
      nbt = 0
      do 6 i = ion, nionmx-1
         nbt = nbt + ntr(i+1)
 6    continue
      weight = weight0(ion)
      rates(1) = 0.0
      do 10 j = 2+nbl,nl(ion)+nbl
         tij = 1.57888e5*eij(j)
         effstr = upsil(k(j),eij(j),c(j),
     .            p1m(j),p2m(j),p3m(j),p4m(j),p5m(j),temp)
         rates(j-nbl) = excdex(temp,tij,weight,effstr)
 10   continue
*        Reset branching ratios to zero
      do 15 i = 1, ntr(ion)
         br(i) = 0.0
 15   continue
*        Add cascades from higher levels.
      do 20 j = nl(ion),2,-1
         jnbl = j + nbl
         ibb = ib(jnbl)
         sum = suma(jnbl)
         if( ibb .gt. 0 ) sum = 1./sum
         do 40 i = 1, ibb
            l = iwork(i,jnbl)
            br(l) = einsta(l+nbt)*sum
            ill = ntrpnti(l+nbt)
            rates(ill) = rates(ill) + br(l)*rates(j)
 40      continue
 20   continue
*        Calculate line intensities (ergs/s)
      i = 0
      do 50 l=1,ntr(ion)
         el = eln(l+nbt)
*         if( el .gt. 0.09999 )then
            i = i + 1
            flx(i)=1.60219e-9*eln(l+nbt)*br(l)*rates( ntrpntj(l+nbt) )
            eline(i) = eln(l+nbt)
*         endif
 50   continue
      nline = i
*      nline = ntr(ion)
*      noutpt = 6
*      write(noutpt,910) ion, nline
* 910  format(I5,I8)
*      write(noutpt,920) ( ntrpnti(l+nbt),ntrpntj(l+nbt),
*     .                eline(l),flx(l),l=1,nline )
* 920  format(2I5,2E14.5)
      return
      end
*
*
*
      REAL FUNCTION EXCDEX( TEMP, TIJ, WEIGHT, EFFSTR )
*
*        Calculates (de)excitation rates per ion per electron.
*        TEMP   - temperature.
*        TIJ    - excitation potential in K (0 for deexcitation).
*        WEIGHT - statistical weight of the lower (excitation) or
*                 upper level (deexcitation).
*        EFFSTR - effective collision strength.
*
      IMPLICIT REAL (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( CONST = 8.629E-6 )
      IF( TIJ .EQ. .0 )THEN
         A = 1.
       ELSE
         A = EXP( -TIJ/TEMP )
      ENDIF
      EXCDEX = CONST/( WEIGHT*SQRT(TEMP) )*A*EFFSTR
      END
*
*
*
C From Leonard J. Moss of SLAC:

C Here's a hybrid QuickSort I wrote a number of years ago.  It's
C based on suggestions in Knuth, Volume 3, and performs much better
C than a pure QuickSort on short or partially ordered input arrays.  

      SUBROUTINE SORTRX(N,DATA,INDEX)
C===================================================================
C
C     SORTRX -- SORT, Real input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(*)
      REAL      DATA(*)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      REAL      DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
C===================================================================
C
C     All done
 
      END
