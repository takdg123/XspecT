      Subroutine xsdiskpn(ear, ne, param, ifl, photar, photer)
      Implicit none
      Integer ne, ifl
      Real ear(0:ne), param(2), photar(ne), photer(ne)
c
c Emission from Shakura-Sunyaev disk. This is an extention of diskbb model,
c including corrections for temperature distribution near the black hole.
c The temperature distribution was calculated in Paczynski-Wiita pseudo-
c Newtonian potential. See Gierlinski et al., 1998, MNRAS, in preparation.
c
c
c Parameters:
c
c   par1 -- maximum temperature in the disk (keV)
c   par2 -- innner disk radius in R_g = GM / c^2 (6 <= par2 <= 1000)
c   K = (M / D)^2 (cos(i) / beta^4) -- normalization, where
c       M - central mass in solar masses
c       D - distance to the source in kpc
c       i - inclination of the disk
c       beta - color/effective temperature ratio
c
c Algorithm:
c
c   n(E) = 2.076e-3 218.02 r_in^2 K kT0^3 / E dE Q(E / par1, Rin)
c
c   Q(E/kT0, Rin) -- tabulated integral
c     = (E/kT0)^3 Integrate[r / (Exp[E/(kT0 t(r, Rin))] - 1), {r,1,Infinity}]
c   t(r) -- temperature distribution over disk, r = R / Rin
c
c Version 3.0, M. Gierlinski, 18 Jan 1999
c
      Real kT0, Rin, Ring
      Real E1, E2, F1, F2, Norm, Tcorr
      Integer i
      Real DiskPNInt

c suppress a warning message from the compiler
      i = ifl

c this model has no errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO

c Er... this should be somhow linked to RMin and RMax in function DiskPNInt

      Ring = param(2)
      if (Ring .lt. 6) then
        call xwrite('Error: R_in < 6 R_g, R_in reset to 6.',5)
        Ring = 6
      elseif (Ring .gt. 1000) then
        call xwrite('Error: R_in > 1000, R_in reset to 1000.',5)
        Ring = 1000
      endif

c Now Rin is in units of R_ms = 6 GM / c^2

      Rin = Ring / 6

c Find maximum temperature in the disk

      if (Rin .lt. 1.5842094) then
        Tcorr = 0.4091842 * Rin**0.75
      else
        Tcorr = 3 * Rin * Rin * (9 * Rin - 1) / (3 * Rin - 1)**3 
     &        * (1 - (3 * Rin - 1) / (2 * Rin**1.5))
        Tcorr = Tcorr**0.25
      endif

c kT0 in our integral is: T0^4 = 3 G M M_dot / 8 Pi R_in^3
c the maximum temperature is: Tcorr * T0

      kT0 = param(1) / Tcorr
      Norm = 2.076e-3 * kT0**3
      E1 = ear(0)

c Change in normalization from diskbb-like to mass:
c     -----------------------------------------
      Norm = Norm * 218.02 * Ring * Ring
c     -----------------------------------------
      F1 = DiskPNInt(E1 / kT0, Rin) / E1
      do i = 1, ne
        E2 = ear(i)
        F2 = DiskPNInt(E2 / kT0, Rin) / E2
        photar(i) = 0.5 * Norm * (F1 + F2) * (E2 - E1)
        E1 = E2
        F1 = F2
      enddo
      return
      end
      
      
      Real Function DiskPNInt(EkT, Rin)
      Implicit none
      Real EkT, Rin
c
c Calculates Integrate[r / (Exp[E/(kT0 t(r, Rin))] - 1), {r, 1, Infinity}]
c interpolating from the tabulated values.
c
c Parameters:
c
c   EkT = E / kT0
c   Rin - inner radius in units of R_ms = 6 GM/c^2
c
      Integer NMax, MMax
      parameter (NMax = 100)
      parameter (MMax = 500)

      Real DiskTab(NMax, MMax)
      Real RMin, RMax, RStep
      Real EMin, EMax, EStep
      Integer N, M
      Real E, R, Pos, Remi, Remj, V, V1, V2
      Integer i, j, ilun, ios
      Logical FirstCall
      Save FirstCall, DiskTab, RMin, RMax, N, EStep, RStep
      Save EMin, EMax, M

      character(256) FileName
      character(128) DataDir

      integer lenact
      character(128) fgmodf
      external lenact, fgmodf

      Data FirstCall /.True./

      if (FirstCall) then
        FirstCall = .False.

        call getlun(ilun)
        DataDir = fgmodf()
        FileName = DataDir(:lenact(DataDir))//'pnint.dat'
        call openwr(ilun, FileName, 'old', ' ', ' ', 0, 0, ios)

        Read(ilun, *, iostat=ios) RMin, RMax, N
        Read(ilun, *, iostat=ios) EMin, EMax, M
        if ((N .gt. NMax) .or. (M .gt. MMax)) then
          call xwrite('Fatal error: internal array overflow',5)
          DiskPNInt = 0.
          return
        endif
        do i = 1, N
          do j = 1, M
            Read(ilun, *) DiskTab(i, j)
          enddo
        enddo

        Close(unit=ilun)
        call frelun(ilun)

        RStep = (RMax - RMin) / (N - 1)
        EStep = (EMax - EMin) / (M - 1)

      endif
      
      E = log10(EkT)
      R = log10(Rin)

      if ((E .lt. EMin) .or. (E .gt. EMax)) then
        DiskPNInt = 0.
      else
        Pos = (E - EMin) / EStep + 1
        j = int(Pos)
        Remj = Pos - j
        Pos = (R - RMin) / RStep + 1
        i = int(Pos)
        Remi = Pos - i
        V1 = DiskTab(i, j) + Remj * (DiskTab(i, j+1) - DiskTab(i, j))
        V2 = DiskTab(i+1,j) + Remj*(DiskTab(i+1,j+1) - DiskTab(i+1,j))
        V = V1 + Remi * (V2 - V1)
        DiskPNInt = 10**V
      endif
      
      return
      end
