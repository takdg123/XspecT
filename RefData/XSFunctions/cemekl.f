C..

      SUBROUTINE cemekl(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(6), photar(ne), photer(ne)

c
c XSPEC model subroutine to calculate Continous Emission
c Measure of the form:
c   Q(T) = Norm*(T/Tmax)^alpha  
c See, for example, Schmitt et al. ApJ 365, 704 (1990), 
c but note that this program yields Schmitt's alpha - 1.0.
c 
c Parameters:
c    param(1) = slope of CEM, alpha
c    param(2) = maximum temperature, tmax
c    param(3) = nH (cm^-3)  Fixed at 1 for most applications
c    param(4) = metallicity used in the Mewe-Kaastra plasma model
c    param(5) = redshift used in the Mewe-Kaastra plasma model
c    param(6) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model)
c
c Robert Dempsey Nov 1, 1991
c K. P. Singh    Feb 25, 1994
c
c Disclaimer: Any resemblance to a real program is purely
c             coincidental
c
c
c Declare variables:
c
 
      REAL pparam(19)

      INTEGER i


      DO i = 1, 3
         pparam(i) = param(i)
      ENDDO
      pparam(4) = 1.
      DO i = 5, 17
         pparam(i) = param(4)
      ENDDO
      pparam(18) = param(5)
      pparam(19) = param(6)

      CALL cevmkl(ear, ne, pparam, ifl, photar, photer)

      RETURN
      END
