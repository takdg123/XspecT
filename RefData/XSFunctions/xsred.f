      SUBROUTINE xsred(e, ne, param, ifl, photar, photer)

      INTEGER ne, ifl
      REAL e(0:ne), param(1), photar(ne), photer(ne)

********************************************************************************
*      Seaton reddening law; valid from 1000Angstroms to 3704 Angstroms        *
*                                                                              *
*        see MNRAS 187 75p (1979)                                              *
*                                                                              *
*                     P. Barr  18.12.90                                        *
*                                                                              *
*************************   INPUT **********************************************
*                                                                              *
*                                                                              *
*       param(1)        XSPEC reddening parameter, E(B-V)                      *
*       photar()        incident spectrum                                      *
*       e()             energy(keV)                                            *
*       ne              number of points in spectrum                           *
*                                                                              *
*************************  OUTPUT **********************************************
*                                                                              *
*       photar()        reddening correction factors                           *
*                                                                              *
********************************************************************************

      REAL x1, x, q1
      INTEGER i

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      DO i = 1, ne

         x = (e(i-1)+e(i))*806.54658/2
c                          !Convert keV to 1/lamda(microns)

         q1 = 1.01/((x-4.6)*(x-4.6)+0.28)
         x1 = 0.0
         IF (x.GE.2.7 .AND. x.LT.3.65) THEN
            x1 = 1.56 + 1.048*x + q1
         ELSEIF (x.GE.3.65 .AND. x.LT.7.14) THEN
            x1 = 2.29 + 0.848*x + q1
         ELSEIF (x.GE.7.14 .AND. x.LE.10.0) THEN
            x1 = 16.17 - 3.2*x + 0.2975*x*x
         ENDIF

         photar(i) = 10.**(-0.4*param(1)*x1)

      ENDDO
      RETURN
      END
