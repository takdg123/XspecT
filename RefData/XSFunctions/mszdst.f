* -------------------------------------------------------------- 
* IR/optical/UV extinction from Pei et al. (1992, ApJ, 395, 130) 
*
*                                                                     
* number of parameters = 4                                           
* param(1) = model (1=galaxy,2=SMC,3=LMC)
* param(2) = E(B-V)                                                   
* param(3) = Rv                                                
* param(4) = redshift (z)                                                   
*                                                                     
*     created: 22-Mar-2005 by M. Still
*                               
* questions and comments: Martin.Still@gsfc.nasa.gov                  
* --------------------------------------------------------------


      SUBROUTINE mszdst(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(4), photar(ne), photer(ne)

* Redshifted reddening law for the galaxy, SMC and LMC

* Arguments :
*     ear       r        i: the energy ranges on which to calculate the model
*     ne        i        i: the number of energy ranges
*     param     r        i: Model, reddening E(B-V), Rv and redshift
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

      call msdust(ear, ne, param, ifl, photar, photer)

* Now shift the energy array back again

      do ie = 0, ne
         ear(ie) = ear(ie) / zfac
      enddo

      return
      end


* ------------------------------------------------------------------- c


      subroutine msdust(e,ne,param,ifl, photar, photer)

      implicit none
      integer ne,ifl
      real e(0:ne),param(*),photar(ne), photer(ne)

      integer i
      real xl,xh,sl,sh,fl,fh,a_b
      real a_i(3,6), lambda_i(3,6), b_i(3,6), n_i(3,6)
      real pei

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

* extinction at B (a_b)

      a_b = param(2) * (1 + param(3))

* parameter coefficients

      call pei_parameters(a_i,lambda_i,b_i,n_i)

* calculate extinction at lambda (a_lambda)

      xl = hc/e(0)
      sl = xl**(-index)
      fl = pei(int(param(1)),xl,a_b,a_i,lambda_i,b_i,n_i)
      do i = 1, ne
         xh = hc/e(i)
         sh = xh**(-index)
         fh = pei(int(param(1)),xh,a_b,a_i,lambda_i,b_i,n_i)
         photar(i) = (sl*fl+sh*fh)/(sl+sh)
         xl = xh
         sl = sh
         fl = fh
      enddo

      return
      end


* ------------------------------------------------------


      function pei(law,rlambda,a_b,a_i,lambda_i,b_i,n_i)

      implicit none

      integer i, law

      real pei, lambda, rlambda, xi, a_lambda, a_b
      real a_i(3,6), lambda_i(3,6), b_i(3,6), n_i(3,6), term(6)

* convert Angstroms to microns

      lambda = rlambda / 1e4

* build function

      xi = 0.0d0
      do i = 1, 6
         term(i) = a_i(law,i)
     &        / ( ( lambda / lambda_i(law,i) )**n_i(law,i)
     &        + ( lambda_i(law,i) / lambda)**n_i(law,i)
     &        + b_i(law,i) )
         xi = xi + term(i)
      end do

* remove a_b normalization on the extinction curve

      a_lambda=a_b*xi
      if (rlambda .lt. 800.) a_lambda = 0.

* linearize extinction factor

      pei=10.**(-a_lambda/2.512)

      return
      end


* --------------------------------------------------------------


      subroutine pei_parameters(a_i,lambda_i,b_i,n_i) 

* Data from Pei, Y.C., 1992 ApJ, 395, 130 (Table 4).

      implicit none

      real a_i(3,6), lambda_i(3,6), b_i(3,6), n_i(3,6)
      
*     Milky Way Extinction Law
      a_i(1,1) = 165.0e0     ! BKG
      a_i(1,2) =  14.0e0     ! FUV
      a_i(1,3) =   0.045e0   ! 2175 AA
      a_i(1,4) =   0.002e0   ! 9.7 um
      a_i(1,5) =   0.002e0   ! 18 um
      a_i(1,6) =   0.012e0   ! FIR

      lambda_i(1,1) =  0.047e0  ! BKG
      lambda_i(1,2) =  0.08e0   ! FUV
      lambda_i(1,3) =  0.22e0   ! 2175 AA
      lambda_i(1,4) =  9.7e0    ! 9.7 um
      lambda_i(1,5) =  18.0e0   ! 18 um
      lambda_i(1,6) =  25.0e0   ! FIR

      b_i(1,1) = 90.0e0    ! BKG
      b_i(1,2) =  4.0e0    ! FUV
      b_i(1,3) = -1.95e0   ! 2175 AA
      b_i(1,4) = -1.95e0   ! 9.7 um
      b_i(1,5) = -1.80e0   ! 18 um
      b_i(1,6) =  0.0e0    ! FIR

      n_i(1,1) = 2.0e0   ! BKG
      n_i(1,2) = 6.5e0   ! FUV
      n_i(1,3) = 2.0e0   ! 2175 AA
      n_i(1,4) = 2.0e0   ! 9.7 um
      n_i(1,5) = 2.0e0   ! 18 um
      n_i(1,6) = 2.0e0   ! FIR

*     Large Magellanic Cloud Extinction Law

      a_i(2,1) = 175.0e0     ! BKG
      a_i(2,2) =  19.0e0     ! FUV
      a_i(2,3) =   0.023e0   ! 2175 AA
      a_i(2,4) =   0.005e0   ! 9.7 um
      a_i(2,5) =   0.006e0   ! 18 um
      a_i(2,6) =   0.020e0   ! FIR

      lambda_i(2,1) =  0.046e0  ! BKG
      lambda_i(2,2) =  0.08e0   ! FUV
      lambda_i(2,3) =  0.22e0   ! 2175 AA
      lambda_i(2,4) =  9.7e0    ! 9.7 um
      lambda_i(2,5) =  18.0e0   ! 18 um
      lambda_i(2,6) =  25.0e0   ! FIR

      b_i(2,1) = 90.0e0   ! BKG
      b_i(2,2) =  5.5e0   ! FUV
      b_i(2,3) = -1.95e0  ! 2175 AA
      b_i(2,4) = -1.95e0  ! 9.7 um
      b_i(2,5) = -1.80e0  ! 18 um
      b_i(2,6) =  0.0e0   ! FIR

      n_i(2,1) = 2.0e0   ! BKG
      n_i(2,2) = 4.5e0   ! FUV
      n_i(2,3) = 2.0e0   ! 2175 AA
      n_i(2,4) = 2.0e0   ! 9.7 um
      n_i(2,5) = 2.0e0               ! 18 um
      n_i(2,6) = 2.0e0   ! FIR

*     Small Magellanic Extinction Law

      a_i(3,1) = 185.0e0     ! BKG
      a_i(3,2) =  27.0e0     ! FUV
      a_i(3,3) =   0.005e0   ! 2175 AA
      a_i(3,4) =   0.010e0   ! 9.7 um
      a_i(3,5) =   0.012e0   ! 18 um
      a_i(3,6) =   0.030e0   ! FIR

      lambda_i(3,1) =  0.042e0  ! BKG
      lambda_i(3,2) =  0.08e0   ! FUV
      lambda_i(3,3) =  0.22e0   ! 2175 AA
      lambda_i(3,4) =  9.7e0    ! 9.7 um
      lambda_i(3,5) =  18.0e0   ! 18 um
      lambda_i(3,6) =  25.0e0   ! FIR

      b_i(3,1) = 90.0e0   ! BKG
      b_i(3,2) =  5.5e0   ! FUV
      b_i(3,3) = -1.95e0  ! 2175 AA
      b_i(3,4) = -1.95e0  ! 9.7 um
      b_i(3,5) = -1.80e0  ! 18 um
      b_i(3,6) =  0.0e0   ! FIR

      n_i(3,1) = 2.0e0   ! BKG
      n_i(3,2) = 4.0e0   ! FUV
      n_i(3,3) = 2.0e0   ! 2175 AA
      n_i(3,4) = 2.0e0   ! 9.7 um
      n_i(3,5) = 2.0e0   ! 18 um
      n_i(3,6) = 2.0e0   ! FIR

      return
      end
