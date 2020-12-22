
      SUBROUTINE xszvab(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(19), photar(ne), photer(ne)

c  works out the redshifted transmission for material whose abundances 
c  are specified for the following 18 elements :
c    1   hydrogen    (1)
c    2   helium      (2)
c    3   carbon      (6)
c    4   nitrogen    (7)
c    5   oxygen      (8)
c    6   neon        (10)
c    7   sodium      (11)
c    8   magnesium   (12)
c    9   aluminium   (13)
c   10   silicon     (14)
c   11   sulphur     (16)
c   12   chlorine    (17)
c   13   argon       (18)
c   14   calcium     (20)
c   15   chromium    (24)
c   16   iron        (26)
c   17   cobalt      (27)
c   18   nickel      (28)
c  The parameters are the column densities of the 18 elements in units of
c  the column of each element in a solar abundance column of equivalent
c  hydrogen column density of 1e22 /cm/cm. Parameter 19 is the redshift.

c Arguments :
c      ear     r        i: energy ranges
c      ne      i        i: number of energies
c      param   r        i: model parameters
c      ifl     i        i: file number
c      photar  r        o: transmitted fraction

      INTEGER NELTS
      PARAMETER (NELTS=18)
      REAL LYLIMIT
      PARAMETER (LYLIMIT=0.0135984)

      INTEGER atomic(NELTS)
      INTEGER ie, status, i

      REAL column(NELTS)
      REAL zfac, elow, ehi, xsect, pret

      character(2) celts(NELTS)

c External references to cross section function

      REAL photo, gphoto, fgabnd
      CHARACTER(4) fgxsct
      EXTERNAL photo, gphoto, fgabnd, fgxsct

      DATA atomic /1,    2,    6,    7,    8,    10,   11,   12,
     &             13,   14,   16,   17,   18,   20,   24,   26,
     &             27,   28/

      DATA celts /'H ', 'He', 'C ', 'N ', 'O ', 'Ne', 'Na', 'Mg',
     &            'Al', 'Si', 'S ', 'Cl', 'Ar', 'Ca', 'Cr', 'Fe',
     &            'Co', 'Ni'/


c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

      status = 0

c Set the element columns

      DO i = 1, NELTS
         column(i) = param(i) * 1.e22 * fgabnd(celts(i))
      ENDDO

      zfac = 1.0 + param(19)
      elow = ear(0) * zfac

      DO ie = 1, ne

         ehi = ear(ie) * zfac
         xsect = 0.

c Set transmission to unity above 100 keV and below the Lyman limit

         IF ( elow .GT. 100. .OR. elow .LT. LYLIMIT ) THEN

            photar(ie) = 1.0

c else calculate the cross-section

         ELSE

            DO i = 1, NELTS
               pret = 0.0
               IF ( fgxsct() .EQ. 'vern' ) THEN
                  pret = gphoto(elow, ehi, atomic(i), status) 
               ELSEIF ( fgxsct() .EQ. 'bcmc' ) THEN
                  pret = photo(elow, ehi, atomic(i), 3, status)
               ELSEIF ( fgxsct() .EQ. 'obcm' ) THEN
                  pret = photo(elow, ehi, atomic(i), 2, status)
               ENDIF
               IF ( column(i) .GT. 0. ) THEN
                  IF ( pret .LT. 1e30/column(i) ) THEN
                     xsect = xsect + pret * column(i)
                  ELSE
                     xsect = 1e30
                  ENDIF
               ENDIF

            ENDDO

            photar(ie) = EXP(-xsect)

         ENDIF

         elow = ehi

      ENDDO

      RETURN
      END



