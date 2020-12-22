c  Routine to return the photoelectric absorption cross section in cm**2
c  at energy keV for element Z.

      FUNCTION gphoto(keV1, keV2, Z, status)

      REAL gphoto, keV1, keV2
      INTEGER Z, status

c Calls Dima Verner's routine

c  kaa  8/16/97

c  Arguments :
c	keV1	r	i: Lower energy of bin in keV.
c	keV2	r	i: Upper energy of bin in keV.
c       Z       i	i: Atomic number of element
c	status	i	r: 0 = OK
c			   1 = No data for this element
c	photo	r	r: Cross-section in cm**2

      REAL xsect
      INTEGER i

      REAL helxsc
      EXTERNAL helxsc

      IF ( Z .LT. 1 .OR. Z .GT. 30 ) THEN
         status = 1
         gphoto = 0.
         RETURN
      ENDIF
      status = 0

c Loop round shell numbers

      gphoto = 0.
      DO i = 1, 7

         CALL phfit2(Z, Z, i, keV1*1000., xsect)
         gphoto = gphoto + 0.5*xsect
         CALL phfit2(Z, Z, i, keV2*1000., xsect)
         gphoto = gphoto + 0.5*xsect

      ENDDO

      gphoto = gphoto*1.e-18

      END


