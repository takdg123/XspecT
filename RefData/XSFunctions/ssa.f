      SUBROUTINE ssa(e0, ne, param, ifl, photar)

      IMPLICIT NONE

      INTEGER ne, ifl
      REAL e0(0:ne), param(2), photar(ne)
      REAL tmp0, tmp1, tmp2
      real      te,a,b,c,ea,h,ee,eb,y
      integer      i

      parameter (h=0.5/3.)  ! constant for Simpson's rule

c     silly line to suppress a compiler warning
      i = ifl
      
      te = param(1)
      y = param(2)
      
      ee = e0(0)
      tmp0 = 1.0334E-3*(ee)**(2)/(exp(ee/te)-1)
      tmp1 = -0.0039*(ee)**(-3.5)*y*(1-exp(-ee/te))
      tmp2 = 1-exp(tmp1)
      a = tmp0*tmp2
    
      DO i = 1, ne
        ea=0.5*(e0(i-1)+e0(i))
         tmp0 = 1.0334E-3*(ea)**(2)/(exp(ea/te)-1)
         tmp1 = -0.0039*(ea)**(-3.5)*y*(1-exp(-ea/te))
         tmp2 = 1-exp(tmp1)
         b =tmp0*tmp2
         eb = e0(i)
         tmp0 = 1.0334E-3*(eb)**(2)/(exp(eb/te)-1)
         tmp1 = -0.0039*(eb)**(-3.5)*y*(1-exp(-eb/te))
         tmp2 = 1-exp(tmp1)
         c = tmp0*tmp2
        photar(i)=h*(e0(i)-e0(i-1))*(a+c+4*b) ! integration with Simpson's rule
         a = c 
      ENDDO

      RETURN
      END
