**==RECOMB.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
      SUBROUTINE RECOMB(Alfa,T,Ic)
c-- version 1.1 j.s. kaastra 06-04-1990
C*** Start of declarations inserted by SPAG
      REAL Alfa , diel , T , t13 , t15 , t25 , tlo , tsq , y , 
     &       YLIM
      INTEGER Ic , ind , k
C*** End of declarations inserted by SPAG

      INCLUDE 'meka.inc'

      PARAMETER (YLIM=69.0)
      DIMENSION Alfa(NUMION)
      IF ( Ic.EQ.6 ) THEN
c		!fe xvii smith
         arec(1,168) = 1.
         arec(4,168) = 0.5601
         arec(5,168) = 2.1754E-11
      ENDIF
      tlo = LOG(T)
      tsq = SQRT(T)
      t13 = T**0.33333333
      t15 = T**1.5
      t25 = T**2.5
      DO 100 k = 1 , NUMION
         ind = NINT(arec(1,k))
         IF ( ind.GT.0 ) THEN
c---------------------------------------------
c autres ions
c---------------------------------------------
            Alfa(k) = arec(2,k)/T**arec(3,k)
c=======================================================
c recombinaison dielectronique
c========================================================
            IF ( ind.EQ.2 ) THEN
c----------------------------------------------------------------
c correction des fits de shull a haute temperature pour
c nevii, neviii, svi, sxiii, sxiv, caxviii
c----------------------------------------------------------------
               IF ( T.GT.arec(4,k) ) THEN
                  Alfa(k) = Alfa(k) + arec(5,k)/T**arec(6,k)
               ELSE
                  y = arec(7,k)/T
                  IF ( y.LT.YLIM ) THEN
                     diel = arec(8,k)*EXP(-y)/t15
                     y = arec(9,k)/T
                     IF ( y.LT.YLIM ) 
     &                 diel = diel*(1.+arec(10,k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  ENDIF
               ENDIF
            ELSEIF ( ind.EQ.3 ) THEN
c---------------------------------------------
c calcul de na, al par interpolation
c---------------------------------------------
               y = arec(4,k)/T
               IF ( y.LT.YLIM ) THEN
                  diel = arec(5,k)*EXP(-y)/t15
                  y = arec(6,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+arec(7,k)*EXP(-y))
     &                 **arec(8,k)
                  y = arec(9,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+arec(10,k)*EXP(-y))
     &                 **arec(11,k)
                  Alfa(k) = Alfa(k) + diel
               ENDIF
            ELSEIF ( ind.EQ.4 ) THEN
c----------------------------------------------------------------
c correction des fits de shull a basse temperature pour si ii
c----------------------------------------------------------------
               IF ( T.GT.arec(4,k) ) THEN
                  y = arec(5,k)/T
                  IF ( y.LT.YLIM ) THEN
                     diel = arec(6,k)*EXP(-y)/t15
                     y = arec(7,k)/T
                     IF ( y.LT.YLIM ) diel = diel*(1.+arec(8,k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  ENDIF
               ELSE
                  Alfa(k) = Alfa(k) + arec(9,k)/T**arec(10,k)
               ENDIF
            ELSEIF ( ind.EQ.5 ) THEN
c----------------------------------------------------------------
c correction des fits de shull a basse temperature pour si iii
c----------------------------------------------------------------
               IF ( T.GT.arec(4,k) ) THEN
                  y = arec(5,k)/T
                  IF ( y.LT.YLIM ) THEN
                     diel = arec(6,k)*EXP(-y)/t15
                     y = arec(7,k)/T
                     IF ( y.LT.YLIM ) diel = diel*(1.+arec(8,k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  ENDIF
               ELSE
                  y = arec(9,k)/T
                  IF ( y.LT.YLIM ) Alfa(k) = Alfa(k) + arec(10,k)
     &                 /t15*EXP(-y)
               ENDIF
            ELSEIF ( ind.EQ.6 ) THEN
c-----------------------------------------------------------------
c recombinaison de cii,ciii,civ;nii,niii,niv,nv;oii,oiii,oiv,ov,ovi
c a basse temperature nussbaumer et storey
c-----------------------------------------------------------------
               y = arec(4,k)/T
               IF ( y.LT.YLIM ) THEN
                  diel = arec(5,k)*EXP(-y)/t15
                  y = arec(6,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+arec(7,k)*EXP(-y))
                  Alfa(k) = Alfa(k) + diel
               ENDIF
               IF ( T.LT.arec(8,k) ) THEN
                  y = arec(9,k)/T
                  IF ( y.LT.YLIM ) Alfa(k) = Alfa(k) + 
     &                 arec(10,k)*EXP(-y)
     &                 /t25*(((T+arec(11,k))*T+arec(12,k))*T+arec(13,k))
               ENDIF
            ELSEIF ( ind.EQ.7 ) THEN
c---------------------------------------------
c autres ions
c---------------------------------------------
               y = arec(4,k)/T
               IF ( y.LT.YLIM ) THEN
                  diel = arec(5,k)*EXP(-y)/t15
                  y = arec(6,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+arec(7,k)*EXP(-y))
                  Alfa(k) = Alfa(k) + diel
               ENDIF
            ELSE
c------------------------------------------------------
c recombinaison des heliumoides de younger
c calculation de smith et al. pour fe xvii (in4=6)
c------------------------------------------------------
               y = arec(4,k)/T
               IF ( y.LT.YLIM ) Alfa(k) = Alfa(k) + 
     &            arec(5,k)/t15*EXP(-y)
            ENDIF
c========================================================
c recombinaison radiative
c========================================================
c------------------------------------------------------
c ions hydrogenoides ( seaton)
c------------------------------------------------------
         ELSEIF ( ind.EQ.0 ) THEN
            Alfa(k) = (arec(2,k)-arec(3,k)*tlo+arec(4,k)*t13)/tsq
         ELSE
            Alfa(k) = 0.
         ENDIF
 100  CONTINUE
      END
**==IONIS.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
      SUBROUTINE IONIS(Sion,T)
C*** Start of declarations inserted by SPAG
      REAL auto , f1 , f2 , f4 , fb , FU , fuy , sa , sd , Sion , 
     &       sk , T , t15 , tsq , x , y , yay , yfuy , YLIM
      INTEGER ia , ind , j , k , n, nj
C*** End of declarations inserted by SPAG

      INCLUDE 'meka.inc'

      PARAMETER (YLIM=69.0)
      DIMENSION Sion(NUMION)

      t15 = T**1.5
      tsq = SQRT(T)
      ia = 0
      DO 200 k = 1 , NUMION
         Sion(k) = 0.
         ia = ia + 1
         nj = NINT(aion(ia))
         IF ( nj.NE.0 ) THEN
            sa = 0.
            DO 60 j = 1 , nj
c-----------------------------------------------
c parametres de younger corriges par les mesures
c-----------------------------------------------
               y = aion(ia+1)/T
               x = 1./y
               IF ( x.GT.0.055 ) THEN
                  f4 = (-5.725E-4+x*(0.01345+x*(0.8691+0.03404*x)))
     &                 /(1.+x*(2.197+x*(0.2454+2.053E-3*x)))
               ELSE
                  f4 = .7699*x**1.9496
               ENDIF
               fuy = FU(y)
               IF ( y.LT.10.0 ) THEN
                  fb = 1. + y - y*fuy*(2.+y)
               ELSE
c calcul du terme en b par integration separee
                  CALL FACTB(y,fb)
               ENDIF
               sa = sa + EXP(-y)
     &              *((aion(ia+2)*(1.-y*fuy)+aion(ia+3)*fb+
     &              aion(ia+4)*fuy)
     &              /y+aion(ia+5)*f4)
               ia = ia + 5
 60         CONTINUE
            Sion(k) = sa/t15
         ENDIF
c--------------------------------
c excitation autoionisation
c--------------------------------
         ia = ia + 1
         ind = NINT(aion(ia))
         IF ( ind.EQ.1 ) THEN
         ELSEIF ( ind.EQ.3 ) THEN
c----------------------
c cas des lithium-like
c----------------------
            y = aion(ia+1)/T
            IF ( y.LT.YLIM ) THEN
               fuy = FU(y)
               yfuy = y*fuy
               sa = EXP(-y)
     &              *(aion(ia+2)*fuy+aion(ia+3)*(1.0-yfuy)+
     &              aion(ia+4)*yfuy+
     &              aion(ia+5)*y*(1.0-yfuy))
               Sion(k) = Sion(k) + sa/tsq
            ENDIF
            ia = ia + 5
         ELSEIF ( ind.EQ.4 ) THEN
c---------------------------------------
c autoionisation deduites des mesures
c---------------------------------------
            y = aion(ia+1)/T
            IF ( y.LT.YLIM ) THEN
               Sion(k) = Sion(k) + aion(ia+2)*EXP(-y)*(1.0-y*FU(y))/tsq
            ENDIF
            ia = ia + 2
         ELSEIF ( ind.EQ.5 ) THEN
            GOTO 100
         ELSEIF ( ind.EQ.6 ) THEN
c-------------------------------------------
c processus dit 'reda' .(la gattuta et hahn)
c-------------------------------------------
            sk = 0.
            DO 80 j = 1 , 12
               y = aion(ia+1)/T
               IF ( y.LT.YLIM ) THEN
                  f1 = EXP(-y)*(y+1.)
                  y = aion(ia+2)/T
                  IF ( y.LT.YLIM ) THEN
                     f2 = EXP(-y)*(y+1.)
                     sk = sk + (f1-f2)*aion(ia+3)
                  ENDIF
               ENDIF
               ia = ia + 3
 80         CONTINUE
            Sion(k) = Sion(k) + sk*tsq
            GOTO 100
         ELSE
c---------------------------------------------------------
c mesures pour ca i,ca ii et calcul de griffin pour le fer
c---------------------------------------------------------
            y = aion(ia+1)/T
            IF ( y.LT.YLIM ) THEN
               auto = aion(ia+2)*EXP(-y)*(1.+aion(ia+3)*FU(y))/tsq
               Sion(k) = Sion(k) + auto
            ENDIF
            ia = ia + 3
         ENDIF
         GOTO 200
c-----------------------------------
c autionisation de sampson(1982)
c-----------------------------------
 100     sd = 0.
         DO 150 n = 1 , 18
            y = aion(ia+1)/T
            IF ( y.LT.YLIM ) THEN
               yay = y*aion(ia+2)
               fuy = FU(yay)
               sd = sd + EXP(-y)
     &              *(FU(y)*aion(ia+3)+aion(ia+4)+y*(aion(ia+5)*fuy+
     &              aion(ia+6)
     &              *(1.-yay*fuy)))
            ENDIF
            ia = ia + 6
 150     CONTINUE
         Sion(k) = Sion(k) + sd/tsq
 200  CONTINUE
      IF ( ia.GT.NAMAX ) WRITE (*,*) 'incompatible dataset'
      RETURN
      END
**==FU.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      FUNCTION FU(U)
c====================================
c approximation from m-s-(1978)
C*** Start of declarations inserted by SPAG
      REAL FU , U
C*** End of declarations inserted by SPAG
      IF ( U.LT.1.0 ) THEN
         FU = LOG(1.+1./U) - (0.36+0.03/SQRT(U+0.01))/(1.+U)**2
      ELSE
         FU = LOG(1.+1./U) - (0.36+0.03*SQRT(U+0.01))/(1.+U)**2
      ENDIF
      RETURN
      END
**==FACTB.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      SUBROUTINE FACTB(Y,Fb)
c=====================================
C*** Start of declarations inserted by SPAG
      REAL c , Fb , ff , sb , t , Y
      INTEGER i
C*** End of declarations inserted by SPAG
      DIMENSION t(16) , c(16)
      DATA t/.8764941047E-1 , .4626963289 , 1.141057774 , 2.129283645 , 
     &     3.437086633 , 5.078018614 , 7.070338535 , 9.438314336 , 
     &     12.21422336 , 15.44152736 , 19.18015685 , 23.51590569 , 
     &     28.57872974 , 34.58339870 , 41.94045264 , 51.70116033/
      DATA c/.2061517149 , .3310578549 , .2657957776 , .1362969342 , 
     &     .4732892869E-1 , .1129990008E-1 , .1849070943E-2 , 
     &     .2042719153E-3 , .1484458687E-4 , .6828319330E-6 , 
     &     .1881024841E-7 , .2862350242E-9 , .2127079033E-11 , 
     &     .6297967002E-14 , .50504737E-17 , .416146237E-21/
      sb = 0.
      DO 100 i = 1 , 16
         ff = t(i)/(t(i)+Y)
         sb = sb + ff**2*c(i)
 100  CONTINUE
      Fb = sb
      RETURN
      END
**==IONCON.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
      SUBROUTINE IONCON(T,Ic,Xe,Xzin,Alfa,Sion,Ed)
c     ******************************************************************
c     *            calculation of ion equilibrium concentrations       *
c     *              1: arnaud and rothenflug			       *
c     *              4: input from calling routine		       *
c     *								       *
c     *  version 1.00 js kaastra 22-12-1989                            *
c     *  version 1.01 js kaastra 10-01-1989 (ic=3 modified and ic=4 added)
c     *  version 1.02 js kaastra 10-01-1989 (storage of alfa,sion added)
c     *  version 1.10 js kaastra 22-01-1989 (ionstx included)	       *
c     *  version 1.11 js kaastra 14-02-1990			       *
c     *  version 1.12 js kaastra 20-02-1990 (ic=4 modified)            *
c     *  version 1.20 js kaastra 06-04-1990 (ic=5,6,7 added; 1 redef.) *
c     ******************************************************************
C*** Start of declarations inserted by SPAG
      REAL Alfa , Ed , fn , Sion , so , T , Xe , Xzin
      INTEGER i , Ic , iel , j , k , nion , NIONMAX , NOEL
C*** End of declarations inserted by SPAG

      INCLUDE 'meka.inc'

      PARAMETER (NOEL=15,NIONMAX=29)
      DIMENSION Xe(NOEL) , Xzin(NUMION) , Alfa(NUMION) , Sion(NUMION) , 
     &          fn(NIONMAX) , nion(NOEL)
      DATA nion/2 , 3 , 7 , 8 , 9 , 11 , 12 , 13 , 14 , 15 , 17 , 19 , 
     &     21 , 27 , 29/
      IF ( Ic.EQ.2 ) THEN
c-- mewe and gronenschild
         STOP
      ELSEIF ( Ic.EQ.3 ) THEN
c-- jacobs and raymond-smith
         STOP
      ELSEIF ( Ic.EQ.4 ) THEN
         GOTO 400
      ELSEIF ( Ic.EQ.5 ) THEN
         STOP
      ELSEIF ( Ic.EQ.7 ) THEN
         CALL IONIS(Sion,T)
         CALL RECOMB(Alfa,T,Ic)
      ELSE
c-- arnaud and rothenflug
         CALL IONIS(Sion,T)
         CALL RECOMB(Alfa,T,Ic)
      ENDIF
c-- calculate the equilibrium concentrations
      k = 2
      Xzin(1) = 0.0
      Xzin(2) = 1.0
c				! hydrogen fully ionized
      DO 200 iel = 2 , NOEL
         fn(1) = 1.0
         so = 1.0
         DO 50 i = 2 , nion(iel)
            fn(i) = fn(i-1)*Sion(k+i-1)/Alfa(k+i)
c to prevent overflow
            IF ( fn(i).GT.1.E25 ) THEN
               DO 10 j = 1 , i
                  fn(j) = fn(j)/1.E25
 10            CONTINUE
               so = so/1.E25
            ENDIF
            so = so + fn(i)
 50      CONTINUE
         DO 100 i = 1 , nion(iel)
            k = k + 1
            Xzin(k) = fn(i)/so
 100     CONTINUE
 200  CONTINUE
c-- correct with the elemental abundances and calculate electron density
      Ed = 0.
      k = 0
      DO 300 iel = 1 , NOEL
         DO 250 i = 1 , nion(iel)
            k = k + 1
            IF ( Xzin(k).LT.1E-5 ) THEN
               Xzin(k) = 0.
            ELSE
               Xzin(k) = Xzin(k)*Xe(iel)
               Ed = Ed + Xzin(k)*FLOAT(i-1)
            ENDIF
 250     CONTINUE
 300  CONTINUE
      RETURN
c-- input from main program
 400  Ed = 0.
      k = 0
      DO 500 iel = 1 , NOEL
         DO 450 i = 1 , nion(iel)
            k = k + 1
            Ed = Ed + Xzin(k)*FLOAT(i-1)
 450     CONTINUE
 500  CONTINUE
      RETURN
      END
**==FLINE.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
      SUBROUTINE FLINE(E1,E2,Flux,Nbin,Ic,Itype,Iem,Cem,H,T,Factor,
     &                 Dilfac,Radt,Xzin,Sion,Alfa,Ed)
c***********************************************************************
c
c this subroutine calculates the ionisation balance and calculates
c the x-ray spectrum
c
c input: e1(nbin)    - lower boundaries energy bin in kev
c        e2(nbin)    - upper boundaries energy bin in kev
c        nbin        - number of energy bins
c        ic          - type of ion concentrations:
c                      1: arnaud and rothenflug
c                      2: mewe and gronenschild
c                      3: jacobs et al.
c                      4: input from calling routine
c                      5: raymond and smith
c                      6: arnaud and rothenflug; fe xvii smith
c                      7: arnaud and rothenflug; dielect. recomb. mewe
c        itype       - 0: only transition rates and ion concentrations
c                      1: continuum and line emission
c                      2: only continuum emission
c                      3: only line strengths and energies:
c                         in this case, input e1(1) as the lower energy
c                         boundary, e1(2) - as the upper energy boundary;
c                         on output e2(1)--e2(nbin) contain the line
c                         energies, flux(1)--flux(nbin) - the line fluxes
c                         or emissivities as specified by iem, except
c                         not per kev (so: phot/cm**2/s  or  phot/m**3/s).
c                         also nbin is now output!
c        iem         - type of spectrum:
c                       1: spectrum   in 1e-4 phot/m**2/s/kev
c                       2: emissivity in      phot/m**3/s/kev
c        cem         - only for for iem = 1:
c                         emitting volume / distance**2
c                         (units: 1e50 cm**3 / pc**2 = 1.0503e11 m)
c        h           - hydrogen density (in cm**-3)
c        t           - temperature in kev
c        factor(15)  - elemental abundances w.r.t. solar values
c        dilfac      - dilution factor
c        radt        - radiation temperature (k)
c
c output:flux(nbin)  - iem=1: spectrum (phot/cm**2/s/kev)
c                      iem=2: emissivity (phot/m**3/s/kev)
c                             (to be multiplied only by the volume)
c        xzin(207)   - ionic concentrations w.r.t. hydrogen
c                      note: for ic=4 this last array is input!
c        alfa(207)   - recombination rates (cm**3/s)
c        sion(207)   - ionisation rates (cm**3/s)
c        ed          - electron density w.r.t. hydrogen
 
c version 1.0  j.s. kaastra 10-01-1990
c version 1.10 j.s. kaastra 23-02-1990 (doppler broadening included)
c version 1.11.j.s. kaastra 06-04-1990
c version 1.2  j.s. kaastra 21-05-1990 (itype = 3 included)
c version 1.3  piotr majer  24-05-1990 (cafe dr *-satellites included)
c version 1.31 j.s. kaastra 28-02-1991 (continuum flux as extra output)
c
c
c***********************************************************************
C*** Start of declarations inserted by SPAG
      REAL Alfa , c , Cem , cst , de , Dilfac , E1 , E2 , Ed , elden , 
     &       elx , Factor , Flux , flx , H , Radt , Sion , T , xe , Xzin
      INTEGER i , Ic , Iem , init , Itype , Nbin , nl , nlx , 
     &          NL_MAX , NOEL
C*** End of declarations inserted by SPAG

      INCLUDE 'meka.inc'

      PARAMETER (NL_MAX=3000)
      PARAMETER (NOEL=15)
      DIMENSION xe(NOEL) , Xzin(NUMION) , Alfa(NUMION) , Sion(NUMION) , 
     &          E1(Nbin) , E2(Nbin) , Flux(Nbin) , Factor(NOEL) , 
     &          elx(NL_MAX) , flx(NL_MAX) , c(2)
c-- c(1) = [2**0.5 * h**2 * alfa**3] /
c--     [(3 * pi * me)**1.5 * (1000 * e)**0.5 * pi * (1 pc)**2 * 1e4]
c--            where all input quantities are in si
c-- c(1) = c(2) * 1.e40 / 4 / pi / (1 pc)**2
      SAVE init
      DATA init/0/
      DATA c/2.53325E-3 , 3.03103E-9/

c-- load the basic data
      IF ( init.EQ.0 ) THEN
         CALL LDDATA()
         init = 1
      ENDIF
c-- set abundances
      CALL ABUN(Factor,xe)
c-- calculate ion concentrations
      CALL IONCON(T,Ic,xe,Xzin,Alfa,Sion,Ed)
      IF ( Itype.EQ.0 ) RETURN
c-- calculate continuum flux
      IF ( Itype.NE.3 ) CALL CONEMX(E1,E2,Flux,Nbin,T,Xzin)
c-- calculate line flux
      IF ( Itype.EQ.1 ) THEN
         elden = Ed*H
         CALL LINEMX(E1(1),E2(Nbin),T,Xzin,elden,Dilfac,
     &               Radt,elx,flx,nlx)
         nl = nlx
         DO 50 i = 1 , Nbin
            de = E2(i) - E1(i)
            DO WHILE ( nl.GT.0 .AND. elx(nl).LE.E2(i) )
               Flux(i) = Flux(i) + flx(nl)/de
               nl = nl - 1
            ENDDO
 50      CONTINUE
      ENDIF
c-- calculate line fluxes separately
      IF ( Itype.EQ.3 ) THEN
         elden = Ed*H
         CALL LINEMX(E1(1),E1(2),T,Xzin,elden,Dilfac,Radt,
     &               E2,Flux,Nbin)
      ENDIF
c-- calculate photon spectrum
      cst = c(Iem)*H**2*Ed/SQRT(T)
      IF ( Iem.EQ.1 ) cst = cst*Cem
      DO 100 i = 1 , Nbin
         Flux(i) = Flux(i)*cst
 100  CONTINUE
      RETURN
      END
**==MGTNAM.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      SUBROUTINE MGTNAM(Filion,Filrec,Filrem,Filjac,Filras,Fillin,
     &                   Filspe,Filkar,Filcaf)
c**************************************************************************
*
c this subroutine gets the names of the basic data files
c
c version 1.1   06-04-1990 j.s. kaastra
c
c**************************************************************************
      character(255) Filion , Filrec , Filrem , Filjac , Filras , 
     &              Fillin , Filspe , Filkar , Filcaf , datdir ,
     &              fgmodf
      INTEGER lenact, dirlen
      EXTERNAL lenact

      datdir = fgmodf()
      dirlen = lenact(datdir)
      Filion = datdir(:dirlen)//'ionis.fits'
      Filrec = datdir(:dirlen)//'recomb.fits'
      Filrem = datdir(:dirlen)//'recombm.fits'
      Fillin = datdir(:dirlen)//'linefile.fits'
      Filjac = datdir(:dirlen)//'ionjac.fits'
      Filras = datdir(:dirlen)//'ionras.fits'
      Filspe = datdir(:dirlen)//'conemx.fits'
      Filkar = datdir(:dirlen)//'karlat.fits'
      Filcaf = datdir(:dirlen)//'cafelines.fits'
      RETURN
      END
**==ABUN.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      SUBROUTINE ABUN(Factor,Xe)
c****************************************************************
c								*
c this subroutine fills xe with the elemental abundances	*
c solar abundances from 					*
c anders and grevesse, geochim. cosmochim. acta 53, 197 (1989)  *
c								*
c  input: factor(15): abundances relative to solar values	*
c output:     xe(15): abundances w.r.t. hydrogen		*
c								*
c	version 2.0 14-02-1990 j.s. kaastra			*
c	version 3.0 10-09-1991 j.s. kaastra			*
c****************************************************************
C*** Start of declarations inserted by SPAG
      REAL abund , Factor , Xe
      INTEGER i , NOEL
C*** End of declarations inserted by SPAG
      PARAMETER (NOEL=15)
      DIMENSION Xe(NOEL) , abund(NOEL) , Factor(NOEL)
      DATA abund/12.00 , 10.99 , 8.56 , 8.05 , 8.93 , 8.09 , 6.33 , 
     &     7.58 , 6.47 , 7.55 , 7.21 , 6.56 , 6.36 , 7.67 , 6.25/
      DO 100 i = 1 , NOEL
         Xe(i) = 10.**(abund(i)-12.)*Factor(i)
 100  CONTINUE
      RETURN
      END
**==GLINE.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
      FUNCTION GLINE(Y,A,B,C,D,E)
c*******************************************************************
c
c this function calculates the gaunt-factor for line emission
c cf. mewe and schrijver 1978, a.a. 65,99, eqn. (36)
c (valid for g(u) = a + b/u + c/u**2 + 2*d/u**3 + e*log(u) )
c
c new approximations are made; the overal accuracy is better than
c (relative error):
c
c 0.0 < y < 0.3 : .00008
c 0.3 < y < 1.0 : .00009
c 1.0 < y < 10. : .00005
c 10. < y       : .00035
c
c version 1.0 j.s. kaastra 12-04-1990
c
c*******************************************************************
C*** Start of declarations inserted by SPAG
      REAL A , B , C , D , E , e1 , f2 , f3 , f4 , f5 , GLINE , Y
C*** End of declarations inserted by SPAG
      IF ( Y.LT.1.0 ) THEN
         IF ( Y.LT.0.3 ) THEN
            e1 = -.57731566 + Y*(.99999193+Y*(-.24991055+Y*.05519968))
     &           - LOG(Y)
            f5 = e1*EXP(Y)
            f2 = Y*f5
         ELSE
            f2 = (.06225196+Y*(1.646421+Y*1.040425))
     &           /(.85539+Y*(2.754082+Y))
            f5 = f2/Y
         ENDIF
         f3 = Y*(1.-f2)
         f4 = Y*(1.-f3)
      ELSE
         IF ( Y.LT.10. ) THEN
            f2 = (.250621+Y*(2.334733+Y))/(1.681534+Y*(3.330657+Y))
            f3 = Y*(1.-f2)
            f4 = Y*(1.-f3)
         ELSE
            f4 = (.37559+Y*(-.6993774+Y*2.))/(-2.520804+Y*(2.617596+Y))
            f3 = 1. - f4/Y
            f2 = 1. - f3/Y
         ENDIF
         f5 = f2/Y
      ENDIF
      GLINE = A + B*f2 + C*f3 + D*f4 + E*f5
      RETURN
      END
**==LINEMX.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
      SUBROUTINE LINEMX(E1,E2,Xkt,Xzin,Elden,Dilfac,Radt,
     &                  Elx,Flx,Nlx)
c     ******************************************************************
c     *                                                                *
c     *                calculation of line emission                    *
c     *                                                                *
c     *  formulas are based on mewe and gronenschild (1981),           *
c     *  astronomy and astrophysics suppl. and contain many            *
c     *  modifications made in later versions                          *
c     *                                                                *
c     *  version 1 october 10, 1980                                    *
c     *  version 2 june 7    , 1983                                    *
c     *  version 3 january 20, 1984                                    *
c     *  version 4 june 5    , 1984                                    *
c     *  version 5 january 15, 1985                                    *
c     *  version 6 bert van den oord  14-6-1985                        *
c     *  version 7 jelle kaastra       3-1-1990			       *
c     *  version 7.01 jelle kaastra   23-2-1990			       *
c     *  version 7.02 jelle kaastra   11-4-1990			       *
c     *  version 7.10 piotr majer     24-5-1990 (cafe dr *-sat. lines) *
c     *                          (also reads new style line data file) *
c     *  version 7.23 jelle kaastra   13-3-1991			       *
c     ******************************************************************
C*** Start of declarations inserted by SPAG
      REAL a , dcor , Dilfac , E1 , E2 , el , Elden , 
     &       Elx , eta , f , fg , fgc , Flx , g , GLINE , Radt , t
      REAL tau , Xkt , xz , Xzin , y
      INTEGER i , lp3 , Nlx , nzz
C*** End of declarations inserted by SPAG

      INCLUDE 'meka.inc'

      DIMENSION Xzin(*) , Elx(*) , Flx(*)
      INTEGER iaux
 
c-- calculation of linestrengths
      t = Xkt*1.16048E+7
      tau = Xkt*1.16048
      eta = Elden*1E-12
      Nlx = 0
      DO 100 i = 1 , ne
         el = rp(1,i)
         IF ( el.LT.E1 ) RETURN
         IF ( el.LT.E2 ) THEN
            nzz = lp(4,i)
            xz = Xzin(nzz)
            IF ( xz.GT.0. ) THEN
               a = rp(3,i)
               y = el*a/Xkt
               IF ( y.LT.40. ) THEN
c--       cf. mewe and gronenschild eqs. 26, 27
                  g = GLINE(y,rp(4,i),rp(5,i),rp(6,i),rp(7,i),rp(8,i))
                  f = rp(2,i)
                  iaux = nrl(i)
                  CALL CORFG(f,g,a,fg,el,y,t,Xzin,nzz,rp(1,i),lp(1,i),
     &                       trans(i),Elden,Dilfac,Radt,iaux,idnr,
     &                       cdrcafe,ncafe)
                  fgc = EXP(-y)*xz*1.646E+5
                  IF ( eta.GT.1E-10 ) THEN
                     IF ( ABS(rp(11,i)).GT.1E-10 .OR. ABS(rp(14,i))
     &                    .GT.1E-10 ) THEN
                        lp3 = lp(3,i)
                        CALL DENCOR(dcor,eta,tau,rp(1,i),lp3)
                        fgc = fgc*dcor
                     ENDIF
                  ENDIF
                  Nlx = Nlx + 1
                  Elx(Nlx) = el
                  Flx(Nlx) = fg*fgc/el
c--                line photon flux (apart from a few constants)
               ENDIF
            ENDIF
         ENDIF
 100  CONTINUE
      RETURN
      END
**==EXPON1.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      FUNCTION EXPON1(X)
c
c        this function calculates the first exponential integral
c        polynomial approximation from abramowitz and stegun,pg 231 (1970)
c
C*** Start of declarations inserted by SPAG
      REAL a , EXPON1 , X
C*** End of declarations inserted by SPAG
      DIMENSION a(6)
      DATA a/ - .57721566 , .99999193 , -.24991055 , .05519968 , 
     &     -.00976004 , .00107857/
      IF ( X.LE.1. ) THEN
         EXPON1 = -ALOG(X) + ((((a(6)*X+a(5))*X+a(4))*X+a(3))*X+a(2))
     &            *X + a(1)
      ELSE
         EXPON1 = (.250621/X+X+2.334733)*EXP(-X)
     &            /(X*X+3.330657*X+1.681534)
      ENDIF
      RETURN
      END
**==CORFG.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      SUBROUTINE CORFG(F,G,A,Fg,El,Y,T,Xzin,Nzz,Rp,Lp,Trans,Elden,
     &                 Dilfac,Radt,Lnum,Idnr,Cdrcafe,Ncafe)
c
c        this subroutine calculates a correction
c        on the f*g value for direct excitation
c
c     in version 3 arrays cdr,bdr and ciid are modified,
c                  a5d is removed and b6 is added.
c                  also some modifications in formulas and constants.
c
C*** Start of declarations inserted by SPAG
      REAL A , b6 , bdr , bii , brf , bri , c1 , c1f , cdr , Cdrcafe , 
     &       cii , ciid , Dilfac , dr , drs , El , Elden , ex , ex1 , 
     &       ex2
      REAL EXPON1 , F , f5 , f6 , FBRF , Fg , fn , fr , G , HELIUM , 
     &       Radt , rii , Rp , rr , T , t2 , t5 , x , xcdr , Xzin
      REAL xzn , xzn1 , Y , yy , z , z12 , z4
      INTEGER ica , idr , idrs , iel , iex , ife , ihel , ii , isat , 
     &          itr , iz , izion , Ncafe , Nzz
C*** End of declarations inserted by SPAG
      INTEGER Lp(*) , Idnr(Ncafe) , Lnum
      DIMENSION Xzin(*) , Rp(*) , Cdrcafe(Ncafe)
      DIMENSION izion(15) , cdr(19) , bdr(19) , b6(19) , ciid(19) , 
     &          bri(15)
      character(1) aster , dum1str
      character(8) Trans
 
      DATA ica , ife/13 , 14/
c				! iel for ca and fe
      DATA aster/'*'/
      DATA izion/1 , 2 , 6 , 7 , 8 , 10 , 11 , 12 , 13 , 14 , 16 , 18 , 
     &     20 , 26 , 28/
      DATA cdr/21.84 , 47.09 , 87.63 , 31.55 , 76.39 , 114.60 , 21.60 , 
     &     36.63 , 55.61 , 15.29 , 8.74 , 3.79 , 0.284 , 0.728 , 21.84 , 
     &     31.55 , 76.39 , 114.60 , 47.09/
      DATA bdr/0.3 , 0.14 , 0.06 , 0.3 , 0.14 , 0.06 , 0.3 , 0.14 , 
     &     0.06 , 0.27 , 0.26 , 0.24 , 0.22 , 0.21 , 0.3 , 0.3 , 0.14 , 
     &     0.06 , 0.14/
      DATA b6/5E-6 , 3E-5 , 5E-5 , 5E-6 , 3E-5 , 5E-5 , 5E-6 , 3E-5 , 
     &     5E-5 , 5E-6 , 5E-6 , 5E-6 , 5E-6 , 5E-6 , 5E-6 , 5E-6 , 
     &     3E-5 , 5E-5 , 3E-5/
      DATA ciid/0.0 , 0.0 , 0.0 , 0.106 , 0.0 , 0.0 , 0.27 , 0.0 , 0.0 , 
     &     0.28 , 0.28 , 0.46 , 0.33 , 0.3 , 0.0 , 0.106 , 0.0 , 0.0 , 
     &     0.0/
      DATA bri/0.0 , 0.0 , 0.111 , 0.225 , 0.293 , 0.338 , 0.353 , 
     &     0.370 , 0.391 , 0.425 , 0.506 , 0.595 , 0.675 , 0.821 , 
     &     0.838/
 
      iel = Lp(1)
      iz = Lp(2) - 1
      z = FLOAT(izion(iel))
      z4 = z**4
      iex = Lp(5)
      idrs = Lp(6)
      idr = Lp(7)
      ii = Lp(8)
      itr = Lp(9)
      Fg = F*G
      t5 = T*1E-5
      dum1str = Trans(1:1)
 
c
c******* excitation ****************************************************
c
      IF ( iex.NE.0 .AND. iex.NE.2 ) THEN
c--               cf. mewe and gronenschild eqs. 28,29
         yy = (Y+0.6)/(Y+1.15)
         c1f = 1. + 7.*yy*(1.065*(1-bri(iel))+0.4*EXP(-0.21*Y))
         brf = FBRF(z,t5,Elden,Dilfac,Radt,bri(iel),iel)
         IF ( iex.EQ.3 ) THEN
c--               only he5 transitions (iex  =  3)
            ex = (bri(iel)+c1f*(1.-brf)/yy/7.46)*1.065
         ELSE
c--               only he6 transitions (iex  =  1)
            ex = c1f*brf
         ENDIF
      ELSE
         ex = 1.
      ENDIF
      ex = ex*Fg
c
c******* dielectronic recombination of satellites **********************
c
      IF ( idrs.NE.0 ) THEN
c--        cf. mewe and gronenschild eqs. 43-45
         iz = iz + 1
         z12 = (iz+1)**2
         x = El*A/(iz+1)
         c1 = EXP(Y/(1.+z12/iz**3/.019))
         c1 = 8538.1*F*c1/((x*81.160+7.7235)*x+1.)/T
         drs = c1*(A*El)**1.5*SQRT(iz/(iz*iz+13.4))*z12
      ELSE
         drs = 0.
      ENDIF
      xzn1 = Xzin(Nzz+1)
      xzn = Xzin(Nzz)
      dr = 0.
      rr = 0.
      IF ( xzn1.GT.0 ) THEN
c
c******* dielectronic recombination ************************************
c           and
c******* radiative recombination ***************************************
c
         IF ( idr.NE.0 ) THEN
            IF ( idr.EQ.1 ) THEN
c====================================================
c--           hes, lis, bes, bs, cs, ns, os, and fs transitions (idr = 1)
c--           cf. mewe and gronenschild eq. 41
c
               IF ( iel.EQ.ica .OR. iel.EQ.ife ) THEN
 
ccccc         special treatement for fe and ca "*-satellites" !!
                  IF ( dum1str.EQ.aster ) THEN
                     DO 2 isat = 1 , Ncafe
                        IF ( Idnr(isat).EQ.Lnum ) GOTO 4
 2                   CONTINUE
                     PRINT '(a/a,i5,1x,a,f10.4)' , 
     &                     ' error: dr satellite not found' , 
     &                     ' lnum, trans, energy   ' , Lnum , Trans , El
                     STOP
 4                   xcdr = Cdrcafe(isat)
cccc            print '(2(a,i4),a,f10.3)', ' cafe line #',isat,
cccc *                '    master #',lnum,'    cdr =',xcdr
                     GOTO 10
                  ENDIF
               ENDIF
               xcdr = cdr(itr)
 10            dr = xcdr/z4*EXP(bdr(itr)*Y)/(b6(itr)+1./(z-1.)**4)
               rr = Rp(9)*T**Rp(10)
c====================================================
            ELSEIF ( idr.EQ.2 ) THEN
c--                     he4 transitions (idr = 2)
               dr = 11.*EXP(0.284*Y)/(1.+z4*6E-6) + 27.*EXP(0.1*Y)
     &              /(1.+z4*3E-5)
               rr = Rp(9)*T**Rp(10)
            ELSE
c--                     he5 + he6 transitions (idr  =  3, 4)
               ex1 = EXP(0.284*Y)
               ex2 = EXP(0.1*Y)
               f5 = 2.3*ex1 + 215.*ex2/(1.+z4*1.9E-4)
               f6 = 16.*ex1/(1.+z4*7E-5) + 90.*ex2/(1.+z4*8E-5)
               brf = FBRF(z,t5,Elden,Dilfac,Radt,bri(iel),iel)
               dr = HELIUM(bri(iel),brf,f5,f6,idr+2)
               f5 = (iz+1.)**2.62/T**0.31*18.72E-4
               f6 = (iz+1.)**2.08/T**0.04*.575E-4
               rr = HELIUM(bri(iel),brf,f5,f6,idr+2)
            ENDIF
            fr = xzn1/xzn*El/12.3985
            rr = fr*EXP(Y)*rr
            dr = fr*0.47*z4/T*dr
         ELSE
            rr = xzn1/xzn*EXP(Y)*El/12.3985*Rp(9)*T**Rp(10)
         ENDIF
      ENDIF
c
c******* innershell ionization *****************************************
c
      rii = 0.
      IF ( ii.NE.0 ) THEN
         xzn1 = Xzin(Nzz-1)
         xzn = Xzin(Nzz)
         IF ( xzn1.GT.0 ) THEN
c--               cf. mewe and gronenschild eqs. 13-14
            IF ( ii.EQ.2 ) THEN
c--               other transitions
               IF ( itr.NE.0 ) THEN
                  cii = ciid(itr)
                  IF ( itr.GE.7 .AND. itr.LE.14 ) THEN
                     bii = 1.242 + 0.009*(z-26)
c					! be - f satellites (itr = 7-14)
                  ELSE
                     bii = 1.296
c					! li satellites (itr = 4-6 , 16-18)
                  ENDIF
               ELSE
                  bii = 1.164 + 0.008*(z-26)
c					!all other transitions (itr set to 100)
                  cii = 0.3
               ENDIF
               t2 = cii*EXPON1(bii*Y)*1.231/bii/(1.+2.E5*z**(-3.5))
            ELSE
c
c
c
c        cf. mewe and gronenschild note table 2
               yy = (z-2.)**4.3*6E12
               fn = yy + 1.33333333*Elden
               f5 = Elden/fn
               f6 = (yy+Elden/3.)/fn
               brf = FBRF(z,t5,Elden,Dilfac,Radt,bri(iel),iel)
               ihel = 5
               IF ( ii.EQ.1 ) ihel = 6
               t2 = 0.22*EXPON1(1.33*Y)*HELIUM(bri(iel),brf,f5,f6,ihel)
            ENDIF
            rii = xzn1/xzn*EXP(Y)*t2
         ENDIF
      ENDIF
c
c******* total correction on the f*g-value ****************************
c
      Fg = ex + drs + dr + rii + rr
      RETURN
      END
**==FBRF.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      FUNCTION FBRF(Z,T5,Elden,Dilfac,Radt,Bri,Iel)
c
c        this function calculates the branching ratio for the
c        forbidden he-like line.
c        the ratio is dependent on the electron density and
c        the external radiation field
c
C*** Start of declarations inserted by SPAG
      REAL beta1 , beta2 , betac , betar , Bri , Dilfac , Elden , ex , 
     &       FBRF , Radt , T5 , Z
      INTEGER Iel
C*** End of declarations inserted by SPAG
      DIMENSION beta1(15) , beta2(15)
      DATA beta1/0.0 , 0.0 , 1.06E7 , 2.45E6 , 7.04E5 , 8.97E4 , 
     &     3.77E4 , 1.75E4 , 8.55E3 , 4.51E3 , 1.48E3 , 5.71E2 , 
     &     2.57E2 , 4.50E1 , 2.79E1/
      DATA beta2/0.0 , 0.0 , 6.33E4 , 7.57E4 , 8.82E4 , 1.14E5 , 
     &     1.28E5 , 1.42E5 , 1.56E5 , 1.72E5 , 2.05E5 , 2.42E5 , 
     &     2.84E5 , 4.51E5 , 5.24E5/
c--        cf. mewe and gronenschild eq. 20
      IF ( beta2(Iel).GE.69.0*Radt ) THEN
         betar = 0.
      ELSE
         betar = beta1(Iel)*Dilfac/(EXP(beta2(Iel)/Radt)-1.)
      ENDIF
c        cf. mewe and gronenschild eq. 19
      ex = -0.73/Z**0.415
      betac = 330./Z**14.57*(T5/Z**2)**ex
c        cf. mewe and gronenschild eq. 18
      FBRF = 1./(1.+(betac*Elden+betar)*Bri)
      RETURN
      END
**==HELIUM.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      FUNCTION HELIUM(Bri,Brf,F5,F6,Itr)
c
c        this function calculates the contribution of
c        the forbidden and intercombination he-like lines
c
C*** Start of declarations inserted by SPAG
      REAL Brf , Bri , F5 , F6 , HELIUM
      INTEGER Itr
C*** End of declarations inserted by SPAG
      IF ( Itr.NE.6 ) THEN
c--        cf. mewe and gronenschild eq. 17
         HELIUM = F5*Bri + (F6+(1.-Bri)*F5)*(1.-Brf)
      ELSE
c--        cf. mewe and gronenschild eq. 16
         HELIUM = (F6+(1.-Bri)*F5)*Brf
      ENDIF
      RETURN
      END
**==DENCOR.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
 
      SUBROUTINE DENCOR(Dcor,Eta,Tau,rp,Lp3)
c
c        this routine calculates the dependence of
c        the line intensities on the electron density
c
C*** Start of declarations inserted by SPAG

      REAL Dcor , Eta , ex , Tau, Rp
      INTEGER Lp3
C*** End of declarations inserted by SPAG
      DIMENSION Rp(*)
      IF ( Lp3.EQ.5 ) THEN
c--        correction on b transitions cf. mewe and gronenschild eq. 32a
         Dcor = 1. + Rp(11)/(Eta**0.03+120./Tau**0.3/Eta)
      ELSEIF ( ABS(Rp(14)).GT.1E-10 ) THEN
         Dcor = 1. + Rp(14)/(1.+3.5/Eta**0.75)
      ELSE
c--        correction on be,c,n,o and f transitions cf. mewe and g. eq. 32
         ex = -Rp(12)*Eta**Rp(13)
         Dcor = 1. + Rp(11)*(1.-EXP(ex))
      ENDIF
      RETURN
      END
**==CONEMX.spg  processed by SPAG 3.09I  at 14:00 on 18 Nov 1992
      SUBROUTINE CONEMX(E1,E2,Flx,Nemx,T,Xzin)
c************************************************************************
c									*
c  calculation of continuum emission					*
c									*
c  formulas are based on mewe, gronenschild, van den oord (1985)	*
c  free-free gaunt factors by bilinear interpolation of logs of		*
c  the gaunt factor tables of carson (1987)				*
c  relativistic correction on free-free gaunt factor using kylafis	*
c  and lamb (1982)							*
c  free-bound gaunt factor cf. karzas and latter (1961)			*
c  two-photon gaunt factor by interpolation of table in	spitzer and	*
c  greenstein (1951) and graph of dalgarno and drake (1969)		*
c									*
c  version 7.4  j.s. kaastra 12-04-1990 (better two-photon-emission)	*
c  version 7.41 j.s. kaastra 18-06-1990 (small error corrected)		*
c									*
c************************************************************************
C*** Start of declarations inserted by SPAG
      REAL a , a31 , a32 , bf , CMAX , dgfbk , dx , E1 , E2 , e2p , 
     &       ej , f2p , Flx , g , g2a , g2p , gfb , gfbk , gff
      REAL gl , GLINE , phi , psi , T , twz , u , ua , ul , UMAX , 
     &       w , W3PI , x , x50 , xez , xzg , Xzin , y
      INTEGER i , i1 , i2 , iel , ifo , in , in6 ,  
     &          inmax , ip , iu , j , j1 , j2 , j3 , ju , k , k1
      INTEGER n , n0 , n2 , NA , NBF , nel , Nemx , nk , NO2P , 
     &          NOEL , NP
C*** End of declarations inserted by SPAG

      INCLUDE 'meka.inc'

      PARAMETER (NOEL=15,NO2P=2)
      PARAMETER (NP=51,NBF=15,NA=NU+5+6*NBF)
      PARAMETER (UMAX=69.0,CMAX=1.E-3,W3PI=0.551328895)
      DIMENSION E1(Nemx) , E2(Nemx) , Flx(Nemx) , Xzin(NUMION) , 
     &          a(NA,NUMION)
      DIMENSION g2a(NG) , ua(NU) , bf(4,NBF) , phi(NP,NO2P)
      DIMENSION nel(NOEL) , e2p(NOEL,NO2P) , f2p(NOEL,NO2P)

c-- bound-free parameters
      DATA bf/6.9248 , -2.1368 , 1.4878 , -0.5867 , 7.0540 , -2.3572 , 
     &     2.4469 , -1.1122 , 7.0799 , -2.3877 , 2.5579 , -1.1304 , 
     &     7.0873 , -2.3932 , 2.5715 , -1.1090 , 7.0865 , -2.3896 , 
     &     2.5563 , -1.0801 , 7.0837 , -2.3845 , 2.5371 , -1.0550 , 
     &     7.0769 , -2.3769 , 2.5122 , -1.0299 , 7.0728 , -2.3724 , 
     &     2.4966 , -1.0133 , 7.0699 , -2.3692 , 2.4848 , -1.0007 , 
     &     7.0658 , -2.3638 , 2.4690 , -0.9871 , 7.0572 , -2.3588 , 
     &     2.4556 , -0.9764 , 7.0530 , -2.3551 , 2.4443 , -0.9667 , 
     &     7.0557 , -2.3555 , 2.4422 , -0.9626 , 7.0487 , -2.3492 , 
     &     2.4259 , -0.9510 , 7.0489 , -2.3481 , 2.4206 , -0.9455/
c-- two-photon parameters
      DATA nel/2 , 5 , 12 , 20 , 29 , 40 , 52 , 65 , 79 , 94 , 111 , 
     &     130 , 151 , 178 , 207/
c-- e2p = energy 1s - 2s transition of the 2 photon decay
c-- f2p = oscillator strength
c-- phi = emissivity distribution, normalised to integral = 2.
      DATA e2p/.0102 , .0408 , 0.368 , 0.500 , 0.654 , 1.022 , 1.236 , 
     &     1.473 , 1.728 , 2.006 , 2.621 , 3.324 , 4.105 , 6.965 , 
     &     8.051 , .0000 , .0200 , 0.304 , 0.426 , 0.569 , 0.915 , 
     &     1.119 , 1.343 , 1.575 , 1.853 , 2.450 , 3.124 , 3.887 , 
     &     6.666 , 7.798/
      DATA f2p/15*0.415 , 0.000 , 0.274 , 0.645 , 0.671 , 0.691 , 
     &     0.719 , 0.729 , 0.737 , 0.745 , 0.751 , 0.761 , 0.768 , 
     &     0.775 , 0.787 , 0.790/
      DATA phi/0.000 , 0.426 , 0.768 , 1.049 , 1.281 , 1.476 , 1.642 , 
     &     1.784 , 1.905 , 2.010 , 2.101 , 2.181 , 2.252 , 2.314 , 
     &     2.367 , 2.412 , 2.451 , 2.484 , 2.513 , 2.538 , 2.559 , 
     &     2.576 , 2.589 , 2.597 , 2.602 , 2.603 , 2.602 , 2.597 , 
     &     2.589 , 2.576 , 2.559 , 2.538 , 2.513 , 2.484 , 2.451 , 
     &     2.412 , 2.367 , 2.314 , 2.252 , 2.181 , 2.101 , 2.010 , 
     &     1.905 , 1.784 , 1.642 , 1.476 , 1.281 , 1.049 , 0.768 , 
     &     0.426 , 0.000 , 0.000 , 0.122 , 0.323 , 0.617 , 0.920 , 
     &     1.176 , 1.417 , 1.600 , 1.772 , 1.923 , 2.065 , 2.177 , 
     &     2.282 , 2.378 , 2.462 , 2.541 , 2.604 , 2.661 , 2.708 , 
     &     2.742 , 2.777 , 2.803 , 2.820 , 2.832 , 2.850 , 2.859 , 
     &     2.850 , 2.832 , 2.820 , 2.803 , 2.777 , 2.742 , 2.708 , 
     &     2.661 , 2.604 , 2.541 , 2.462 , 2.378 , 2.282 , 2.177 , 
     &     2.065 , 1.923 , 1.772 , 1.600 , 1.417 , 1.176 , 0.920 , 
     &     0.617 , 0.323 , 0.122 , 0.000/

      DO 100 i = 1 , NG
         g2a(i) = -4.25 + FLOAT(i)*0.25
 100  CONTINUE
      DO 150 i = 1 , NU
         ua(i) = -4.25 + FLOAT(i)*0.25
 150  CONTINUE

      iel = 1
      k = 0
      DO 200 i = 1 , NUMION
         k1 = k + 1
         IF ( k1.LT.NUMION ) THEN
            DO 160 j = 1 , NA
               a(j,k1) = 0.
 160        CONTINUE
         ENDIF
         ifo = 0
         IF ( Xzin(i).GT.0. ) THEN
c-- preset help arrays for free-free emission
            xez = Xzin(i)*p(1,i)
            IF ( xez.GT.CMAX ) THEN
               ifo = 1
               k = k + 1
               a(1,k) = xez
               gl = ALOG10(.0136*p(1,i)/T)
               i1 = INT(4.*(gl+4.)) + 1
               IF ( i1.LT.1 ) i1 = 1
               IF ( i1.GE.NG ) i1 = NG - 1
               i2 = i1 + 1
               w = (g2a(i2)-gl)/(g2a(i1)-g2a(i2))
               DO 170 ju = 1 , NU
                  a(ju+1,k) = ga(i2,ju) + w*(ga(i2,ju)-ga(i1,ju))
 170           CONTINUE
            ENDIF
c-- preset help arrays for two photon process
            DO 180 i2 = 1 , NO2P
               IF ( i.EQ.nel(iel)-i2 ) THEN
c			!i2=1 h-like, i2=2 he-like
c       determine g(y) - factor
                  IF ( i2.NE.1 ) THEN
                     g = 0.05
                  ELSEIF ( iel.EQ.1 ) THEN
c					!h i
                     g = GLINE(e2p(iel,i2)/T,0.08,-0.16,0.11,0.0,0.0)
                  ELSE
c					!hydrogen-like ions except hydrogen
                     g = 0.055
                  ENDIF
                  twz = 164995.*f2p(iel,i2)*g*Xzin(i)/e2p(iel,i2)
                  IF ( twz.GT.CMAX ) THEN
c				!only for sufficient amplitude
                     u = e2p(iel,i2)/T
                     IF ( u.LT.UMAX ) THEN
c				!only for significant exponent
                        IF ( ifo.EQ.0 ) THEN
                           ifo = 1
                           k = k + 1
                        ENDIF
                        a(27,k) = i2
                        a(28,k) = e2p(iel,i2)
                        a(29,k) = twz*EXP(-u)
                     ENDIF
                  ENDIF
               ENDIF
 180        CONTINUE
            IF ( i.EQ.nel(iel) ) iel = iel + 1
c-- preset help arrays for free-bound emission
            n0 = NINT(p(2,i))
            in = 0
            IF ( n0.GT.0 ) THEN
               DO 190 n = n0 , NBF
                  n2 = n*n
                  IF ( n.EQ.n0 ) THEN
                     a31 = bf(1,n)*p(3,i)*Xzin(i)/T
c							!prefactor
                     a32 = p(4,i)/T
c							!exponential factor
                  ELSE
                     a31 = bf(1,n)*p(5,i)*Xzin(i)/T/FLOAT(n2*n)
                     a32 = p(6,i)/T/FLOAT(n2)
                  ENDIF
                  IF ( a32.LT.UMAX ) THEN
                     IF ( a31*EXP(a32).GT.CMAX ) THEN
                        IF ( ifo.EQ.0 ) THEN
                           ifo = 1
                           k = k + 1
                        ENDIF
                        in6 = 6*in
                        a(31+in6,k) = a31
c							!prefactor
                        a(32+in6,k) = a32
c							!edge energy / kt
                        a(33+in6,k) = a32*FLOAT(n2)
c						!edge energy / kt for n=1
                        a(34+in6,k) = bf(2,n)
                        a(35+in6,k) = bf(3,n)
                        a(36+in6,k) = bf(4,n)
                        in = in + 1
                     ENDIF
                  ENDIF
 190           CONTINUE
            ENDIF
            IF ( ifo.EQ.1 ) a(30,k) = in
         ENDIF
 200  CONTINUE
      nk = k
c-- loop over energies
      DO 300 j = 1 , Nemx
         gff = 0.
         gfb = 0.
         g2p = 0.
         ej = 0.5*(E1(j)+E2(j))
         u = ej/T
         IF ( u.LT.UMAX ) THEN
            ul = ALOG10(u)
            iu = INT(4.*(ul+4.)) + 1
            IF ( iu.LT.1 ) iu = 1
            IF ( iu.GE.NU ) iu = NU - 1
            j1 = iu
            j2 = j1 + 1
            j3 = j2 + 1
            w = (ua(j2)-ul)/(ua(j1)-ua(j2))
            DO 220 k = 1 , nk
c------------------------------------------------------------------
c free-free radiation
c------------------------------------------------------------------
               IF ( a(1,k).GT.0. ) THEN
                  IF ( u.LT.1.E-4 ) THEN
                     xzg = a(1,k)*(EXP(a(2,k))-W3PI*LOG(u/1.E-4))
                  ELSE
                     xzg = a(1,k)*EXP(a(j3,k)+w*(a(j3,k)-a(j2,k)))
                  ENDIF
                  gff = gff + xzg
               ENDIF
c------------------------------------------------------------------
c free-bound radiation
c------------------------------------------------------------------
               gfbk = 0.
               in = 0
               inmax = NINT(a(30,k))
               dgfbk = 1.
               DO WHILE ( in.LT.inmax .AND. dgfbk.GT.CMAX*gfbk )
c			!stop if n > nfb or relative contribution too small
                  in6 = 6*in
                  IF ( a(31+in6,k).GT.0 ) THEN
c						!only for positive prefactor
                     y = u - a(32+in6,k)
                     IF ( y.GE.0.0 .AND. y.LT.UMAX ) THEN
c						!only for e > e_edge
                        x = 1./(1.+SQRT(u/a(33+in6,k)))
                        dgfbk = (((a(36+in6,k)*x+a(35+in6,k))*x+a(34+in6
     &                          ,k))*x+1.)*x*a(31+in6,k)*EXP(-y)
                        gfbk = gfbk + dgfbk
                     ENDIF
                  ENDIF
                  in = in + 1
               ENDDO
               gfb = gfb + gfbk
c------------------------------------------------------------------
c two photon emission
c------------------------------------------------------------------
               IF ( a(29,k).GT.0. ) THEN
                  IF ( ej.LT.a(28,k) ) THEN
                     i2 = NINT(a(27,k))
                     x = ej/a(28,k)
                     x50 = 50.*x
                     ip = INT(x50)
                     dx = x50 - FLOAT(ip)
                     ip = ip + 1
                     psi = phi(ip,i2) + dx*(phi(ip+1,i2)-phi(ip,i2))
                     g2p = g2p + a(29,k)*x*psi
                  ENDIF
               ENDIF
 220        CONTINUE
            gff = gff*EXP(-u)
         ENDIF
c------------------------------------------------------------------
c total continuum photon flux (to be multiplied by a few constants)
c------------------------------------------------------------------
         Flx(j) = (gff+gfb+g2p)/ej
 300  CONTINUE
      RETURN
      END

c*****************************************************************
      SUBROUTINE lddata()

c Subroutine to load the common block arrays for the MEKA code
c Reads the FITS versions of the input files.

      INCLUDE 'meka.inc'

      INTEGER lun, ierr, i, index

      character(255) Filion , Filrec , Filrem , Filjac , Filras , 
     &              Fillin , Filspe , Filkar , Filcaf
      character(72) comment, context

      LOGICAL qanyf

      INTEGER lenact
      EXTERNAL lenact

c Get the names of the input files

      CALL MGTNAM(Filion,Filrec,Filrem,Filjac,Filras,Fillin,
     &            Filspe,Filkar,Filcaf)

c Read the files

c ----------------------------------------------------------------------

      CALL GETLUN(lun)
      ierr = 0
      CALL ftopen(lun, Filrec, 0, index, ierr)
      IF (ierr .NE. 0) THEN
         context = 'Failed to open '//Filrec(:lenact(Filrec))
         GOTO 10
      ENDIF
      CALL ftmrhd(lun, 1, index, ierr)
      IF (ierr .NE. 0) THEN
         context = 'Failed to go to next extension'
         GOTO 10
      ENDIF

      DO i = 1 , NUMION
         CALL ftgcve(lun, 1, i, 1, 13, 0., arec(1,i), qanyf, ierr)
         IF (ierr .NE. 0) THEN
            WRITE(context,'(a,i5)') 'Failed to read row ', i
            GOTO 10
         ENDIF
      ENDDO
      CALL ftclos(lun, ierr)
      CALL FRELUN(lun)

 10   CONTINUE
      IF (ierr .NE. 0 ) THEN
         CALL xwrite('Fatal error : Failed to read Filrec data', 2)
         CALL xwrite(context, 2)
         WRITE(context, '(a,i5)') 'FITSIO error = ', ierr
         CALL xwrite(context, 2)
         STOP
      ENDIF

c ----------------------------------------------------------------------

      CALL GETLUN(lun)
      ierr = 0
      CALL ftopen(lun, Filion, 0, index, ierr)
      CALL ftmrhd(lun, 1, index, ierr)

      DO i = 1 , NREC
         CALL ftgcve(lun, 1, i, 1, 256, 0., aion((i-1)*256+1), 
     &               qanyf, ierr)
      ENDDO
      CALL ftclos(lun, ierr)

      IF (ierr .NE. 0 ) CALL xwrite('Failed to read Filion data', 2)

c ----------------------------------------------------------------------

      CALL GETLUN(lun)
      ierr = 0
      CALL ftopen(lun, Fillin, 0, index, ierr)
      context = 'Failed to open '//Fillin(:lenact(Fillin))
      IF ( ierr .NE. 0 ) GOTO 30
      CALL ftmrhd(lun, 1, index, ierr)
      context = 'Failed to go to next extension'
      IF ( ierr .NE. 0 ) GOTO 30


      CALL ftgkyj(lun, 'NAXIS2', ne, comment, ierr)
      context = 'Failed to get NAXIS2'
      IF ( ierr .NE. 0 ) GOTO 30

      DO i = 1 , ne
         CALL ftgcvj(lun, 1, i, 1, 1, 0, nrl(i), qanyf, ierr)
         CALL ftgcvj(lun, 2, i, 1, 9, 0, lp(1,i), qanyf, ierr)
         CALL ftgcvs(lun, 3, i, 1, 1, ' ', trans(i), qanyf, ierr)
         CALL ftgcve(lun, 4, i, 1, 14, 0., rp(1,i), qanyf, ierr)
         WRITE(context,'(a,i5)') 'Failed to read row ', i
         IF ( ierr .NE. 0 ) GOTO 30
      ENDDO
      CALL ftclos(lun, ierr)
      CALL FRELUN(lun)

 30   CONTINUE
      IF (ierr .NE. 0 ) THEN
         CALL xwrite('Failed to read Fillin data', 2)
         CALL xwrite(context, 2)
         WRITE(context,'(a,i4)') 'CFITSIO error status = ', ierr
         CALL xwrite(context, 2)
      ENDIF

c ----------------------------------------------------------------------

      CALL GETLUN(lun)
      ierr = 0
      CALL ftopen(lun, Filcaf, 0, index, ierr)
      CALL ftmrhd(lun, 1, index, ierr)

      CALL ftgkyj(lun, 'NAXIS2', ncafe, comment, ierr)

      DO i = 1 , ncafe
         CALL ftgcvj(lun, 1, i, 1, 1, 0, idnr(i), qanyf, ierr)
         CALL ftgcve(lun, 2, i, 1, 1, 0., cdrcafe(i), qanyf, ierr)
      ENDDO
      CALL ftclos(lun, ierr)
      CALL FRELUN(lun)

      IF (ierr .NE. 0 ) CALL xwrite('Failed to read Filcaf data', 2)

c ----------------------------------------------------------------------

      CALL GETLUN(lun)
      ierr = 0
      CALL ftopen(lun, Filspe, 0, index, ierr)
      CALL ftmrhd(lun, 1, index, ierr)

      DO i = 1 , NUMION
         CALL ftgcve(lun, 1, i, 1, 6, 0., p(1,i), qanyf, ierr)
      ENDDO
      CALL ftclos(lun, ierr)
      CALL FRELUN(lun)

      IF (ierr .NE. 0 ) CALL xwrite('Failed to read Filspe data', 2)

c ----------------------------------------------------------------------

      CALL GETLUN(lun)
      ierr = 0
      CALL ftopen(lun, Filkar, 0, index, ierr)
      CALL ftmrhd(lun, 1, index, ierr)

      DO i = 1 , NU
         CALL ftgcve(lun, 1, i, 1, NG, 0., ga(1,i), qanyf, ierr)
      ENDDO
      CALL ftclos(lun, ierr)
      CALL FRELUN(lun)

      IF (ierr .NE. 0 ) CALL xwrite('Failed to read Filkar data', 2)

      RETURN
      END
