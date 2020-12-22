
      SUBROUTINE xszphb(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

c Redshifted photoelectric absorption. Uses the XSPHAB to calculate
c cross-sections

c Arguments :
c     ear       r        i: the energy ranges on which to calculate the model
c     ne        i        i: the number of energy ranges
c     param     r        i: the H column density in 10^22 cm^-2 and redshift
c     ifl       i        i: the dataset
c     photar    r        r: fractional transmission

      REAL zfac
      INTEGER ie

c Shift the energy array to the emitted frame

      zfac = 1.0 + param(2)

      DO ie = 0, ne
         ear(ie) = ear(ie) * zfac
      ENDDO

c Calculate the transmission

      CALL xsphab(ear, ne, param(1), ifl, photar, photer)

c Now shift the energy array back again

      DO ie = 0, ne
         ear(ie) = ear(ie) / zfac
      ENDDO

      RETURN
      END
