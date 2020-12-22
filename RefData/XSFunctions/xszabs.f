
      SUBROUTINE xszabs(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

c	xszabs	rashafer 27 nov 1983
c	includes redshift T. Yaqoob 14 Dec 1992
c		XSPEC model subroutine:
c		Standard absorption, using the 'wisconsin' cross sections
c		and abundances (Morrison and McCammon, 1983 ApJ 270, p119-122)
c		for the case of no grains.
c		For energies above 10 keV, the continuation of the higher
c		energy range is used, assuming a pure E**-3 form.
c		For energies below 0.03 keV, the
c		cross section is assumed to be purely proportional to E**-3.
c             fwj haberl   15-FEB-1990
c               no absorption below 0.0136 keV

c		see MULMOD for parameter descriptions
c	number of model parameters: 2
c		1	ANH	Hydrogen column density (in units of 10**22
c				atoms per square centimeter
c		2		redshift
c	intrinsic energy range.
c		Emin=0.03, Emax=infinity (although an attempted continuation
c			to lower energies is made)
c	algorithm:
c		a = exp(-sigma(E)*Nh), where sigma(E) is derived from
c		an analytic fit to cross section data

      INTEGER NRANGE
      PARAMETER (NRANGE=16)
      INTEGER ir, ie
      REAL edge(NRANGE), c0(NRANGE), c1(NRANGE), c2(NRANGE)
      REAL anh, cedge, cc2, cc1, cc0, cear, e, zfac
      REAL facum, sigma, einv, cf, eold, oldf, oldear
      LOGICAL qlaste, qedge

      DATA edge/0.03, 0.1, .284, .4, .532, .707, .867, 1.303, 1.840,
     &     2.471, 3.210, 4.038, 7.111, 8.331, 10., 1.e10/
      DATA c0/.336, .173, .346, .781, .714, .955, 3.089, 1.206, 1.413,
     &     2.027, 3.427, 3.522, 4.339, 6.29, 7.012, 9.532/
      DATA c1/0., 6.081, 2.679, .188, .668, 1.458, -3.806, 1.693, 1.468,
     &     1.047, .187, .187, -.024, .309, .252, 0./
      DATA c2/0., -21.5, -4.761, 0.043, -.514, -.611, 2.94, -.477,
     &     -.315, -.17, 2*0., .0075, 3*0./

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

      anh = param(1)
      zfac = 1.0+param(2)

c find edge energy just greater than initial array energy

      ir = 1
      DO WHILE (edge(ir).LE.(ear(0)*zfac))
         ir = ir + 1
      ENDDO
      cedge = edge(ir)
      cc2 = c2(ir)
      cc1 = c1(ir)
      cc0 = c0(ir)
      oldf = 0.0
      eold = 0.0
      DO ie = 0, ne
         cear = ear(ie)*zfac
         e = min(cear, cedge)
         qlaste = .FALSE.
         facum = 0.
         DO WHILE (.NOT.qlaste)
            qedge = (e.EQ.cedge)
            qlaste = (.NOT.qedge) .AND. (e.EQ.cear)

c evaluate the absorption at e

            IF (e.LT.0.0136) THEN
               cf = 1.0
            ELSE
               einv = 1./e
               sigma = einv*(cc2+(einv*(cc1+(einv*cc0))))
               cf = exp(-anh*sigma)
            ENDIF
            IF (qedge .OR. qlaste) THEN
               facum = facum + (e-eold)*(oldf+cf)
               eold = e
               IF (qedge) THEN
                  ir = ir + 1
                  cedge = edge(ir)
                  cc2 = c2(ir)
                  cc1 = c1(ir)
                  cc0 = c0(ir)
                  qedge = .FALSE.
               ELSE

c this can be reached only if qlaste is true

                  oldf = cf
               ENDIF
            ELSE

c since it is not the end of an edge OR the end of the
c ear range, then this must be the initial part of the
c energy range for a new edge

               eold = e
               oldf = cf
               e = min(cedge, cear)
            ENDIF
         ENDDO

         IF (ie.NE.0) THEN
            IF (cear.NE.oldear) THEN
               photar(ie) = facum*0.5/(cear-oldear)
            ELSE
               photar(ie) = 0.0
            ENDIF
         ENDIF
         oldear = cear
      ENDDO

      RETURN
      END
