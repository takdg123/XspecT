
      SUBROUTINE xsphab(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(1), photar(ne), photer(ne)

c "Wisconsin" absorption using the cross-sections of Balucinska-Church
c and McCammon, 1992, ApJ 400, 599. The elemental abundances are those
c set by the abund command.

c Arguments :
c     ear       r        i: the energy ranges on which to calculate the model
c     ne        i        i: the number of energy ranges
c     param     r        i: the H column density in 10^22 cm^-2
c     ifl       i        i: the dataset
c     photar    r        r: fractional transmission

      INTEGER NPARM
      PARAMETER (NPARM=19)

      REAL vparam(NPARM)

      INTEGER i


      vparam(1) = param(1)
      DO i = 2, NPARM-1
         vparam(i) = 1.
      ENDDO
      vparam(NPARM) = 0.

      CALL xszvph(ear, ne, vparam, ifl, photar, photer)

      RETURN
      END
