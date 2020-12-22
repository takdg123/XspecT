
      SUBROUTINE xssssi(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(1), photar(ne), photer(ne)

*- xssssi - Einstein sss ice absorption
* Description :
*  works out the attenuation for ice infront of the Einstein sss for the
*  single parameter of the ice thickness
* Author :
*  Andy Pollock
* History :
*  1 February 1988 : original
*  15 September 1988 : revised clumpy model from Nick White

c     ne       i       i: no of energy bins
c     ifl      i       i: dataset
c     ear      r       i: energy bounds
c     param    r       i: sss ice thickness
c     photar   r       r: ice transmission

* Local variables :
      REAL clumps
c						  ..param(1)
      REAL keV
c						  energy
      INTEGER ie
c						  energy index
* External reference :
      REAL sssica
* Status :
      INTEGER status
*-

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO

      status = 0

      clumps = param(1)
      DO ie = 1, ne
         keV = (ear(ie-1)+ear(ie))/2.
         photar(ie) = sssica(clumps, keV, status)
      ENDDO

      END

*- ice_abs - Einstein sss ice absorbtion model

      REAL FUNCTION sssica(clumps, keV, status)

      REAL clumps, keV
      INTEGER status

* Description :
*  code lifted from GSFC software for use in sss spectrum modelling of
*  the absorption caused by ice condensed on the front of the instrument
* Author :
*  Andy Pollock (EXOSAT::ANDY)
* History :
*  14 December 1987 : original
*  1 September 1988 : clumpy model devised by Nick White
c      implicit none
* Global constants :

c   clumps     r       i: no of clumps of thick ice
c   keV        r       i: energy
c   status     i     i/r: status

* Local constant :
      REAL oxygen_edge
      PARAMETER (oxygen_edge=0.5317)

* Local variables :

      REAL clump2
      REAL t1, t2, v1, v2, z

*-

      sssica = 1.

      IF (status .NE. 0) RETURN

      IF (keV.LE.oxygen_edge) THEN
         z = 213.*keV**(-2.84)
      ELSE
         z = 4860.*keV**(-2.84)
      ENDIF

c new parameters from Damian
c
c      t1 = (0.82+1.35*clumps)*1.0E-04
c      t2 = t1/16.6
c      clump2 = clumps*1.27
c      v1 = z*t1
c      v2 = z*t2

      t1=(1.165+1.944*clumps)*1.0e-04
      t2=t1/26.68
      clump2=clumps*3.02
      v1=z*t1
      v2=z*t2

      sssica = exp(-clumps*(1.-exp(-v1)))
      sssica = sssica*exp(-clump2*(1.-exp(-v2)))

      RETURN
      END
