
      SUBROUTINE xszvph(ear, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL ear(0:ne), param(19), photar(ne), photer(ne)

c calculates the photoelectric absorption using the cross-sections of
c Balucinska-Church and McCammon, 1992, ApJ 400, 599. The elemental
c abundances are specified relative to the ratios set using the abund
c command.
c Parameters :
c       1     H column in 10^22 cm^-2
c       2-18  Relative abundances of
c               He,C,N,O,Ne,Na,Mg,Al,Si,S,Cl,Ar,Ca,Cr,Fe,Co,Ni
c       19    Redshift

c Arguments :
c      ear     r        i: energy ranges
c      ne      i        i: number of energies
c      param   r        i: model parameters
c      ifl     i        i: file number
c      photar  r        o: transmitted fraction


      INTEGER NPARM
      PARAMETER (NPARM=19)

      REAL pparam(NPARM)

      INTEGER i

      pparam(1) = param(1)
      DO i = 2, NPARM-1
         pparam(i) = param(i) * param(1)
      ENDDO
      pparam(NPARM) = param(NPARM)

      CALL xszvab(ear, ne, pparam, ifl, photar, photer)

      RETURN
      END


