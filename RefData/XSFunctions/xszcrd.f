c ------------------------------------------------------------------- c
c IR/optical/UV extinction from Cardelli et al. (1989, ApJ, 345, 245) c
c contains a redshift correction to the redden model (xscred.f)       c
c                                                                     c
c number of parameters = 1                                            c
c param(1) = E(B-V)                                                   c
c param(2) = redshift                                                 c
c                                                                     c
c     created: 22-Mar-2005 by M. Still                                c
c questions and comments: Martin.Still@gsfc.nasa.gov                  c
c ------------------------------------------------------------------- c

      SUBROUTINE xszcrd(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

c Redshifted reddening law. Uses XSCRED.

c Arguments :
c     ear       r        i: the energy ranges on which to calculate the model
c     ne        i        i: the number of energy ranges
c     param     r        i: Reddening E(B-V) and redshift
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

      CALL xscred(ear, ne, param(1), ifl, photar, photer)

c Now shift the energy array back again

      DO ie = 0, ne
         ear(ie) = ear(ie) / zfac
      ENDDO

      RETURN
      END

