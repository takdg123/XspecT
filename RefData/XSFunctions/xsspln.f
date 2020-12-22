      subroutine xsspln(ear,ne,param,ifl, photar, photer)

      integer ne, ifl
      real ear(0:ne),param(6),photar(ne), photer(ne)

C---
C multiplicative spline model for XSPEC
C---
C see MULMOD for parameter descriptions
C number of parameters: 6
C      1      Start x-value
C      2      Start y-value
C      3      End y-value
C      4      Start dy/dx
C      5      End dy/dx
C      6      End x-value
C no intrinsic energy range limit
C Model form M(E) = ax**3 + bx**2 + cx + d where a, b, c, d are
C functions of passed parameters
C---
C 21 may 1990 - kaa
C---
      DOUBLE PRECISION a, b, c, d, el, eh
      DOUBLE PRECISION p, q, w, x, y, z, f
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

c  trap case of no valid energy range

      if (param(6).le.param(1)) return

c  calculate cubic coefficients from passed parameters

      p = param(1)
      w = param(2)
      x = param(3)
      y = param(4)
      z = param(5)
      q = param(6)

      f = 1./(p-q)**3
      a = -(2*w-2*x-p*y+q*y-p*z+q*z)*f
      b = -(-3*p*w-3*q*w+3*p*x+3*q*x+p**2*y+p*q*y-2*q**2*y+
     &            2*p**2*z-p*q*z-q**2*z)*f
      c = (-6*p*q*w+6*p*q*x+2*p**2*q*y-p*q**2*y-q**3*y+p**3*z+
     &            p**2*q*z-2*p*q**2*z)*f
      d = -(-3*p*q**2*w+q**3*w-p**3*x+3*p**2*q*x+p**2*q**2*y-
     &            p*q**3*y+p**3*q*z-p**2*q**2*z)*f

cd
cd      write (*,*) p, w, x, y, z, q
cd      write (*,*) a, b, c, d
cd

c  loop over energies

      do i = 1, ne
         el = ear(i-1)
         eh = ear(i)

         if (el .ne. eh) then

c  if energy range lies partly in the valid interval then apply average over
c  that section

            if (el .lt. p .AND. eh .ge. p) then

               photar(i) = SNGL(((p-el) +
     &                  (d*eh+c*eh**2/2+b*eh**3/3+a*eh**4/4) -
     &                  (d*p+c*p**2/2+b*p**3/3+a*p**4/4))/(eh-el))

c  if energy range lies completely in valid interval then average over entire
c  energy range

            else if (el .gt. p .AND. eh .lt. q) then

               photar(i) = SNGL(((d*eh+c*eh**2/2+b*eh**3/3+a*eh**4/4) -
     &                 (d*el+c*el**2/2+b*el**3/3+a*el**4/4)) / (eh-el))

c  if energy range lies partly in the valid interval then apply average over
c  that section

            else if (eh .gt. q .AND. el .le. q) then

               photar(i) = SNGL(((eh-q) +
     &                  (d*q+c*q**2/2+b*q**3/3+a*q**4/4) -
     &                  (d*el+c*el**2/2+b*el**3/3+a*el**4/4))/(eh-el))

            end if

         end if

      end do

      return
      end
