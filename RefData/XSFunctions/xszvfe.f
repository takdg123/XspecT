
      SUBROUTINE xszvfe(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(5), photar(ne), photer(ne)

c  the redshifted absorption for material with a variable Fe abundance
c  and variable energy Fe K edge.

c T. Yaqoob June 1994 

c Arguments :
c     ear          r       i: energy ranges
c     ne           i       i: number of energies
c     param        r       i: model parameters
c     ifl          i       i: file number
c     photar       r       o: transmitted fraction

c Model parameters are :
c      1    equivalent H column in units of 10^22 cm^-2
c      2    metal abundance (except Fe) relative to Solar
c      3    Fe abundance relative to Solar
c      4    Fe K edge
c      5    Redshift


      REAL rparm(19), zfac, afe
      INTEGER ie, i

      REAL feabs, fgabnd
      EXTERNAL feabs, fgabnd

c Set all the abundances to the input value but zero out Fe.

      rparm(1) = param(1)
      DO i = 2, 18
         rparm(i) = param(2)
      ENDDO
      rparm(16) = 0.
      rparm(19) = param(5)

c Calculate the absorption from all the elements except iron

      CALL xszvph(ear, ne, rparm, ifl, photar, photer)

c Now multiply in the Fe absorption. The feabs function (found
c in xshrfl.f) assumes an [Fe/H] = 3.31e-5 so correct
c this to the current Solar abundance in use.

      afe = param(3) * fgabnd('Fe') / 3.31e-5
      zfac = 1.0 + param(5)
      DO ie = 1, ne
         photar(ie) = photar(ie) * EXP(-0.5*(
     &      feabs(ear(ie-1)*zfac, param(1)*10., afe, param(4))
     &     +feabs(ear(ie)*zfac  , param(1)*10., afe, param(4)) ) ) 
      ENDDO

      RETURN
      END

