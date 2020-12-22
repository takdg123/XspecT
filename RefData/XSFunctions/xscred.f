c ------------------------------------------------------------------- c
c IR/optical/UV extinction from Cardelli et al. (1989, ApJ, 345, 245) c
c                                                                     c
c number of parameters = 1                                            c
c param(1) = E(B-V)                                                   c
c                                                                     c
c created: 02-Apr-1997 by O.Blaes & P.Magdziarz                       c
c questions and comments: pavel@camk.edu.pl                           c
c ------------------------------------------------------------------- c

      subroutine xscred(e,ne,param,ifl, photar, photer)

      implicit none
      integer ne,ifl
      real e(0:ne),param(*),photar(ne), photer(ne)

      integer i
      real xl,xh,sl,sh,fl,fh
      real cardelli

c conversion constant in keV*\AA
      real hc
      parameter (hc=12.3963)

c photon index of the base power law function
      real index
      parameter (index=2)

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      xl=hc/e(0)
      sl=xl**(-index)
      fl=cardelli(xl,param(1))
      do i=1,ne
         xh=hc/e(i)
         sh=xh**(-index)
         fh=cardelli(xh,param(1))
         photar(i)=(sl*fl+sh*fh)/(sl+sh)
         xl=xh
         sl=sh
         fl=fh
      enddo

      return
      end

c *****************************************************************
c Calculates the extinction for a given E(B-V) and wavelength
c according to Cardelli et al. (1989, ApJ, 345, 245).  
c Units of wavelength are angstroms.  Note that
c the Cardelli fits are nominally valid only between 3.3 microns
c and 1250 angstroms.  It is therefore dangerous to use this
c code outside this range.  (Note however that Cardelli et al.
c argue persuasively that their extinction curve is better than
c Seaton's at the shortest UV wavelengths, so even though Seaton's
c extends further, it is probably not as good.
c
c Extension to 900 Angstroms by Martin Still 3/26/05
c *****************************************************************

      function cardelli(rlambda,ebmv)

      implicit none

      real cardelli, rlambda, ebmv

      real rv, av, x, y, ax, bx, al

c note that we assume here that R_V=3.1. Seaton (1979) used R_V=3.2
      rv=3.1
      av=ebmv*rv
      x=1.e4/rlambda
      y=x-1.82

c Warning - it is dangerous to use this for x<0.3 or x>11
      ax = 0
      bx = 0
      if(x.lt.1.1) then
        ax=0.574*x**1.61
        bx=-0.527*x**1.61
      else if(x.le.3.3.and.x.ge.1.1) then
        ax=1.+.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4
     %     +0.01979*y**5-0.77530*y**6+0.32999*y**7
        bx=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4
     %     -0.62251*y**5+5.30260*y**6-2.09002*y**7
      else if(x.le.8.and.x.gt.3.3) then
        if(x.ge.5.9) then
          ax=1.752-0.316*x-0.104/((x-4.67)**2+0.341)
     %       -0.04473*(x-5.9)**2-0.009779*(x-5.9)**3
          bx=-3.090+1.825*x+1.206/((x-4.62)**2+0.263)
     %       +0.2130*(x-5.9)**2+0.1207*(x-5.9)**3
        else
          ax=1.752-0.316*x-0.104/((x-4.67)**2+0.341)
          bx=-3.090+1.825*x+1.206/((x-4.62)**2+0.263)
        endif
      else if(x.le.11.and.x.gt.8.) then
        ax = -1.073-0.628*(x-8.)+0.137*(x-8.)**2-0.07*(x-8.)**3
        bx = 13.670+4.257*(x-8.)-0.42*(x-8.)**2+0.374*(x-8.)**3
      endif
      al=av*(ax+bx/rv)
      cardelli=10.**(-al/2.512)
      return
      end
