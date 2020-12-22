
      SUBROUTINE eplogpar(ear, ne, param, ifl, photar)
      IMPLICIT NONE
      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne)
      real      ep,beta,a,b,c,ea,h
      integer      i

      parameter (h=0.5/3.)      ! constant for Simpson's rule

C-----------------------------------------------------
C number of parameters (except normalization): 2
C model form N(E) = 10**(-beta*(log(E/Ep))**2)/E**2
C-----------------------------------------------------
C Ep: peak energy in nuFnu 
C beta: curvature
C the normalization corresponds to the peak flux height
C-----------------------------------------------------
C First Ver. 19 nov 2005 - Tramacere Andrea
C-----------------------------------------------------

c suppress a warning message from the compiler
      i = ifl
        
      ep=param(1)
      beta=param(2)
        
      a=10.0**(-beta*(Alog10(ear(0)/ep))**2)/ear(0)**2        
        
      DO i = 1, ne
         ea=0.5*(ear(i-1)+ear(i))
         b=10.0**(-beta*(Alog10(ea/ep))**2)/ea**2
         c=10.0**(-beta*(Alog10(ear(i)/ep))**2)/ear(i)**2
         photar(i)=h*(ear(i)-ear(i-1))*(a+c+4.0*b) ! integration with Simpson's rule        
         a=c 
      ENDDO
        
      RETURN
      END
      
