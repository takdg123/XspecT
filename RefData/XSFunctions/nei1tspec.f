c Routine to return the spectrum from a single-temperature NEI plasma. 
c Simply a wrap-up of Kazik's code.

c Arguments
c    temp         r        i: temperature (keV)
c    nion         i        i: size of ion arrays (may be > actual no. of ions)
c    ionfrac      r        i: fraction in ion
c    ne           i        i: number of response energies
c    ear          r        i: response energies
c    z            r        i: redshift
c    nelt         i        i: number of elements
c    outphot      r        r: output spectrum for each element

      SUBROUTINE nei1tspec(temp, nion, ionfrac, ne, ear, z, nelt, 
     &                     outphot)

      INTEGER ne, nion, nelt

      REAL temp, ionfrac(nion), ear(0:ne), z
      REAL outphot(ne, nelt)

      real t
      real zf
*        Conversion factors.
      real keV, boltz, tconv
      parameter ( keV = 1.60219e-9, boltz = 1.38062e-16 )
      parameter ( tconv = keV/boltz )

*
*          nelem         -  number of heavy elements
*          ielem         -  curent index of elements (2 for He^++,
*                           3-12 for heavy elements
*          istart, iend  -  dummy integer indexes
      integer nelem, nelemp
      parameter (nelem=10,nelemp=nelem+2)
      integer nlines(3:nelemp), ielem, istart, iend
*        nline - number of lines, el(iline) - energy of line iline in keV
      integer nline, iline
      parameter (nline=566)
      real el(nline)
c
c                  xline(iline) = line luminosity of line iline
c      xline(iline) = lines
c
      real xline(nline)
      integer nmaxfe
      parameter ( nmaxfe = 146000 )
      real elfe(nmaxfe), xlinefe(nmaxfe)
      integer nfe
      integer ifreq

      save nlines, xline, elfe, xlinefe, nfe

      data nlines / 15, 15, 15, 20, 22, 48, 70, 71, 157, 133 /

      t = tconv*temp
      zf = 1./(1. + z)

* Line emission.
      call emisl(t,ionfrac,xline,el)
* Ne- to Li-like Fe ions
      call felines(t,ionfrac(104),elfe,xlinefe,nfe)

*           Continua.
      DO ielem = 1, nelemp
         DO ifreq = 1, ne
            outphot(ifreq,ielem) = 0.0
         ENDDO
      ENDDO
      call continua(t,ionfrac,ear,ne,z,outphot)

*           Add lines into energy bins.
      istart = 1
      iend = 0
      do 2000 ielem=3,nelemp
         iend = iend + nlines(ielem)
         do 1200 iline=istart,iend           
            call nei1hunt(zf*el(iline),ear,ne,ifreq)
            if( ifreq .ne. 0) outphot(ifreq,ielem) =
     .           outphot(ifreq,ielem) + xline(iline)/(keV*el(iline))
 1200    continue
         istart = iend + 1
 2000 continue

*           Add Fe L-shell lines.
      do 3200 iline=1,nfe
         call nei1hunt(zf*elfe(iline),ear,ne,ifreq)
         if( ifreq .ne. 0) outphot(ifreq,11) = 
     &        outphot(ifreq,11) + xlinefe(iline)/(keV*elfe(iline))
 3200 continue

      RETURN
      END

c********************************************************************
      subroutine nei1hunt(e,ear,ne,jlo)
c   returns # of bin which contains energy e
      integer ne, jlo
      real e,ear(*)
      integer inc,jhi,jm
      if( e .lt. ear(1) .or. e .ge. ear(ne) )then
         jlo = 0
         return
      endif
      if( jlo .le. 0 .or. jlo .gt. ne )then
         jlo = 0
         jhi = ne + 1
         go to 3
      endif
      inc = 1
      if( e .ge. ear(jlo) )then
 1       jhi = jlo + inc
         if( jhi .gt. ne )then
            jhi = ne + 1
          elseif( e .ge. ear(jhi) )then
            jlo=jhi
            inc = inc + inc
            go to 1
         endif
       else
         jhi = jlo
 2       jlo = jhi - inc
         if( jlo .lt. 1 )then
            jlo = 0
          elseif( e .lt. ear(jlo) )then
            jhi = jlo
            inc = inc + inc
            go to 2
         endif
      endif
 3    if( jhi-jlo .eq. 1 )then
         if( e .eq. ear(1) ) jlo=1
         return
      endif
      jm = ( jhi + jlo )/2
      if( e .ge. ear(jm) )then
         jlo = jm
       else
         jhi = jm
      endif
      go to 3
      end

