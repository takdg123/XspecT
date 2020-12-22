      SUBROUTINE xssmdg(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(4), photar(ne), photer(ne)

c	xssmdg	adopted from xsedge by rashafer 27 nov 1983
c		   by femarshall 25 spt 1992
c		XSPEC model subroutine:
c		reproduces "smeared edge" described in Appendix C
c			of k.ebisawa's thesis
c		general "edge" multiplicative absorption program, using
c		an approximation suitable for photoionization.
c	1.1: A zero value for the threshold energy is equivalent to no edge at
c		all.
c		see MULMOD for parameter descriptions
c	number of model parameters: 4
c		1	edgeE	the threshold energy
c		2	maxAbs	the maximum absorption factor at threshold
c		3	index for photoelectric cross-section
c			    normally (-2.67)
c		4	wiflh for smearing
c	intrinsic energy range.
c		from epsilon to infinity
c	algorithm
c		m(e) = 1. for E<edgeE
c		m(e) = exp(-maxabs*((E/edgeE)**index)(1-exp((edgeE-E))/wiflh))

      REAL cedgei, abs, alow, ahi, xold, x, amax, aindex, wiflh
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
      aindex = param(3)
      wiflh = param(4)

      x = ear(0)*cedgei
      IF (x.GE.1.) THEN
         qedge = .TRUE.
         alow = exp( abs*(x**aindex) *
     &           (1.-exp( (param(1)-ear(0))/wiflh )))
      ELSE
         qedge = .FALSE.
      ENDIF
      DO ie = 1, ne
         x = ear(ie)*cedgei
         IF (qedge) THEN
            ahi = exp(abs*(x**aindex) *
     &           (1.-exp( (param(1)-ear(ie))/wiflh )))
            photar(ie) = 0.5*(alow+ahi)
            alow = ahi
         ELSEIF (x.GE.1.) THEN
            qedge = .TRUE.
            amax = 1.
            alow = exp(abs*(x**aindex) *
     &           (1.-exp( (param(1)-ear(ie))/wiflh )))
            xold = ear(ie-1)*cedgei
            photar(ie) = ((1.*(1.-xold))+(0.5*(amax+alow)*(x-1.)))
     &                   /(x-xold)
         ELSE
            photar(ie) = 1.
         ENDIF
      ENDDO
      RETURN
      END
