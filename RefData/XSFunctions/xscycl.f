
      SUBROUTINE xscycl(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(5), photar(ne), photer(ne)

C---
C XSPEC multiplicative model subroutine:
C interpretes simple exponential high energy cutoff as part of cyclotron
C feature in pulsar spectra
C
C see Mihara etal. submitted to Nature, April 1990
C or Makishima etal., Publ. Astr. Soc. Japan, in press
C and references therein
C---
C see MULMOD for parameter descriptions
C number of model parameters: 5
C       1       A0: Depth of fundamental
C       2       E0: Cyclotron energy
C       3       W0: Wiflh of the fundamental
C       4       A2: Depth 2nd harmonic
C       5       W2: Wiflh of the 2nd harmonic
C model form A(E) = exp{-A0*(W0*E/E0)**2/[(E-E0)**2+W0**2]
C                       +A2*(W2*E/(2*E0))**2/((E-2*E0)**2+W2**2))
C---
C 6 june 1990 - fwj haberl
C---
      REAL a0, e0, w0, a2, w2, w0s, w2s, e0s, ep, eps, a, b, h
      REAL x1, x2, y1, y2
      INTEGER i

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

C---
      a0 = param(1)
      e0 = param(2)
      w0 = param(3)
      a2 = param(4)
      w2 = param(5)
      w0s = w0*w0
      w2s = w2*w2
      e0s = e0*e0
      ep = ear(0)
      eps = ep*ep
      x1 = a0*w0s/e0s
      x2 = a2*w2s/e0s/4.
      y1 = ep - e0
      y2 = y1 - e0
      h = x1*eps/(y1*y1+w0s) + x2*eps/(y2*y2+w2s)
      a = exp(-h)
      DO i = 1, ne
         ep = ear(i)
         eps = ep*ep
         y1 = ep - e0
         y2 = y1 - e0
         h = x1*eps/(y1*y1+w0s) + x2*eps/(y2*y2+w2s)
         b = exp(-h)
         photar(i) = 0.5*(a+b)
         a = b
      ENDDO
      RETURN
      END
