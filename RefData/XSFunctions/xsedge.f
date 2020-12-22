      SUBROUTINE xsedge(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

c	xsedge	rashafer 27 nov 1983
c		XSPEC model subroutine:
c		general "edge" multiplicative absorption program, using
c		an approximation suitable for photoionization.
c	1.1: A zero value for the threshold energy is equivalent to no edge at
c		all.
c		see MULMOD for parameter descriptions
c	number of model parameters: 2
c		1	edgeE	the threshold energy
c		2	maxAbs	the maximum absorption factor at threshold
c	intrinsic energy range.
c		from epsilon to infinity
c	algorithm
c		m(e) = 1. for E<edgeE
c		m(e) = exp(-maxabs*((E/edgeE)**-3))

      REAL cedgei, abs, alow, ahi, xold, x, amax
      INTEGER ie
      LOGICAL qedge

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO


      IF (param(1).LE.0.) THEN
         DO ie = 1, ne
            photar(ie) = 1.
         ENDDO
         RETURN
      ELSE
         cedgei = 1./param(1)
      ENDIF
      abs = -param(2)
      x = ear(0)*cedgei
      IF (x.GE.1.) THEN
         qedge = .TRUE.
         alow = exp(abs/(x*x*x))
      ELSE
         qedge = .FALSE.
      ENDIF
      DO ie = 1, ne
         x = ear(ie)*cedgei
         IF (qedge) THEN
            ahi = exp(abs/(x*x*x))
            photar(ie) = 0.5*(alow+ahi)
            alow = ahi
         ELSEIF (x.GE.1.) THEN
            qedge = .TRUE.
            amax = exp(abs)
            alow = exp(abs/(x*x*x))
            xold = ear(ie-1)*cedgei
            photar(ie) = ((1.*(1.-xold))+(0.5*(amax+alow)*(x-1.)))
     &                   /(x-xold)
         ELSE
            photar(ie) = 1.
         ENDIF
      ENDDO
      RETURN
      END
