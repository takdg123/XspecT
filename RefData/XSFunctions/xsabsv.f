
      SUBROUTINE xsabsv(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(18), photar(ne), photer(ne)

c  works out the attenuation for material whose abundances are specified for
c  the following 18 elements :
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
c  hydrogen column density of 1e22 /cm/cm.

c Arguments :
c      ear     r        i: energy ranges
c      ne      i        i: number of energies
c      param   r        i: model parameters
c      ifl     i        i: file number
c      photar  r        o: transmitted fraction

c Author :
c  Andy Pollock
c History :
c  25 January 1988 : original
c   4 April   1992 : new photo code used

      INTEGER NPARM
      PARAMETER (NPARM=19)

      REAL pparam(NPARM)

      INTEGER i

      DO i = 1, NPARM-1
         pparam(i) = param(i)
      ENDDO
      pparam(NPARM) = 0.

      CALL xszvab(ear, ne, pparam, ifl, photar, photer)

      RETURN
      END









