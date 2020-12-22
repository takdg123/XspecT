* -------------------------------------------------------------- 
* Powerlaw extinction from Savaglio & Fall (2004, ApJ, 614, 293)
*
*                                                                     
* number of parameters = 4   
* param(1) = E(B-V)                                                   
* param(2) = Rv
* param(3) = gamma, powerlaw index 
* param(4) = redshift (z)                                                   
*                                                                     
*     created: 30-Jul-2005 by M. Still
*                               
* questions and comments: Martin.Still@gsfc.nasa.gov                  
* --------------------------------------------------------------


      SUBROUTINE msldst(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(4), photar(ne), photer(ne)

* Redshifted reddening law (lambda^gamma).

* Arguments :
*     ear       r        i: the energy ranges on which to calculate the model
*     ne        i        i: the number of energy ranges
*     param     r        i: E(B-V), Rv, index and redshift
*     ifl       i        i: the dataset
*     photar    r        r: fractional transmission

      real zfac
      integer ie

* shift the energy array to the emitted frame

      zfac = 1.0 + param(4)
 
      do ie = 0, ne
         ear(ie) = ear(ie) * zfac
      enddo

* calculate the transmission

      call msldust(ear, ne, param, ifl, photar, photer)

* Now shift the energy array back again

      do ie = 0, ne
         ear(ie) = ear(ie) / zfac
      enddo

      return
      end


* ------------------------------------------------------------------- c


      subroutine msldust(e,ne,param,ifl, photar, photer)

      implicit none
      integer ne,ifl
      real e(0:ne),param(*),photar(ne), photer(ne)

      integer i
      real xl,xh,sl,sh,fl,fh,a_v
      real savaglio

* conversion constant in keV*\AA

      real hc
      parameter (hc=12.3963)

* photon index of the base power law function

      real index
      parameter (index=2)

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

* extinction at V (a_v)

      a_v = param(1) * param(3)

* calculate extinction at lambda (a_lambda)

      xl = hc/e(0)
      sl = xl**(-index)
      fl = savaglio(xl,a_v,param(2))
      do i = 1, ne
         xh = hc/e(i)
         sh = xh**(-index)
         fh = savaglio(xh,a_v,param(2))
         photar(i) = (sl*fl+sh*fh)/(sl+sh)
         xl = xh
         sl = sh
         fl = fh
      enddo

      return
      end


* ------------------------------------------------------


      function savaglio(rlambda,a_v,gamma)

      implicit none

      real savaglio, lambda, rlambda, xi, a_lambda, a_v, gamma

* convert Angstroms to microns

      lambda = rlambda / 1e4

* build function

      xi = SNGL((5.5d3 / rlambda)**gamma)

* remove a_v normalization on the extinction curve

      a_lambda=a_v*xi
      if (rlambda .lt. 800.) a_lambda = 0.

* linearize extinction factor

      savaglio=10.**(-a_lambda/2.512)

      return
      end
