
      subroutine xsexp(ear,ne,param,ifl, photar, photer)

      integer ne, ifl
      real ear(0:ne),param(3),photar(ne), photer(ne)

C---
C multiplicative model for XSPEC  1-param(1)*exp(-param(2)*E)
C---
C see MULMOD for parameter descriptions
C number of parameters: 2
C      1      amplitude
C      2      exponential factor
C      3      beginning energy
C no intrinsic energy range limit
C Model form M(E) = 1.-param(1)*exp(-param(2)*Energy(keV))  for E > param(3)
C            M(E) = 1.                                          E < param(3)
C---
C  07 aug 1991  kj
C---
      REAL el, eh
      integer i
C---

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      do i = 1, ne
          photar(i) = 1.
      end do

c  loop over energies

      do i = 1, ne
         el = ear(i-1)
         eh = ear(i)

         if (el .ne. eh) then

c  if energy range lies partly in the valid interval then apply average over
c  that section

            if (el .lt. param(3) .AND. eh .ge. param(3)) then

               photar(i) = (param(3)-el) +
     &              eh-param(3) - (param(1)/param(2)) *
     &              (exp(-param(2)*eh) - exp(-param(2)*param(3)))
               photar(i) = photar(i)/(eh-el)

c  if energy range lies completely in valid interval then average over entire
c  energy range

            else if (el .gt. param(3))  then

               photar(i) = 1 - ( (param(1)/param(2)) *
     &              (exp(-param(2)*eh)-exp(-param(2)*el)))
     &              / (eh-el)

            end if

         end if

      end do

      return
      end
