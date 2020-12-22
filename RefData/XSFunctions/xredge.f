
      SUBROUTINE xredge(ear, ne, param, ifl, photar, photer)

c Input variable
      INTEGER ne, ifl
      REAL ear(0:ne), param(2), photar(ne), photer(ne)

c	recedge	
c		see ADDITIVE for parameter descriptions
c	number of model parameters: 2
c		1	edgeE	the threshold energy
c		2	KT      plasma temperature
c	algorithm
c		A(e) = 0. for E<edgeE
c		A(e) = (1/KT)*exp(-(E-edgeE)/KT)) dE
c the norm is equal to the integrated number of counts under the line

      REAL kT, edge
      REAL elow, ehigh, flow, fhigh
      INTEGER ie

c suppress a warning message from the compiler
      ie = ifl

c this model does not calculate errors
      DO ie = 1, ne
         photer(ie) = 0.0
      ENDDO



      edge = param(1)
      kT   = param(2)


      elow = MAX(0.,(ear(0)-edge)/kT)
      flow = exp(-elow)

      DO ie = 1, ne

         ehigh = MAX(0.,(ear(ie)-edge)/kT)
         fhigh = exp(-ehigh)

         photar(ie) = flow - fhigh

         elow = ehigh
         flow = fhigh

      ENDDO

      RETURN
      END

