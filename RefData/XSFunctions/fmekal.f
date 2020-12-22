**==fmekal.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
      SUBROUTINE FMEKAL(E1,E2,Flux,Nbin,Cem,H,T,Factor,Ed)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL abund , Cem , cst , de , E1 , E2 , Ed , elden , elx , 
     &     Factor , Flux , flx , H , T , xe , xzin
      INTEGER i , ierr , index , init , lun , Nbin , nl , nlx , NL_MAX , 
     &        NOEL
C*** End of declarations inserted by SPAG
c***********************************************************************
c
c input: e1(nbin)    - lower boundaries energy bin in keV
c        e2(nbin)    - upper boundaries energy bin in keV
c        nbin        - number of energy bins
c        cem         -
c                         emitting volume / distance**2
c                         (units: 1e50 cm**3 / pc**2 = 1.0503e11 m)
c        h           - hydrogen density (in cm**-3)
c        t           - temperature in keV
c        factor(15)  - elemental abundances w.r.t. solar values
c
c output:flux(nbin)  - spectrum (phot/cm**2/s/keV)
c        ed          - electron density w.r.t. hydrogen
c
c***********************************************************************
      INCLUDE 'mekal.inc'
      PARAMETER (NL_MAX=5500,NOEL=15)
      DIMENSION xe(NOEL) , xzin(NUMION) , E1(Nbin) , E2(Nbin) , 
     &          Flux(Nbin) , abund(NOEL) , Factor(NOEL) , elx(NL_MAX) , 
     &          flx(NL_MAX)
      character(256) filion , filrec , filkar , filcon , fillin , filcaf
      character(72) comment
      character(256) wrtstr
 
      LOGICAL qanyf

      INTEGER lenact, dirlen
      CHARACTER(256) fgmodf, datdir
      EXTERNAL lenact, fgmodf
 
      DATA abund/12.00 , 10.99 , 8.56 , 8.05 , 8.93 , 8.09 , 6.33 , 
     &     7.58 , 6.47 , 7.55 , 7.21 , 6.56 , 6.36 , 7.67 , 6.25/
      DATA init /0/

c Initialize outputs

      ed = 0.0
      DO i = 1, nbin
         Flux(i) = 0.
      ENDDO

c
c-- get basic data file names and initialise
c
      IF ( init.EQ.0 ) THEN
         init = 1
         datdir = fgmodf()
         dirlen = lenact(datdir)
         filion = datdir(:dirlen)//'mekal1.dat'
         filrec = datdir(:dirlen)//'mekal2.dat'
         filkar = datdir(:dirlen)//'mekal3.dat'
         filcon = datdir(:dirlen)//'mekal4.dat'
         fillin = datdir(:dirlen)//'mekal5.dat'
         filcaf = datdir(:dirlen)//'mekal6.dat'
 
         lun = 80
         ierr = 0
 
c
c	open (unit=80,access='direct',file=filion,
c     *  recl=256*4,status='old',form='unformatted')
c	do i = 1,nrec
c	  read (80,rec=i) (aion(j),j=(i-1)*256+1,i*256)
c	enddo
c	close(unit=80)

         CALL ftopen(lun,filion,0,index,ierr)
         CALL FTMRHD(lun,1,index,ierr)
         DO 50 i = 1 , NREC
            CALL FTGCVE(lun,1,i,1,256,0.,AIOn((i-1)*256+1),qanyf,ierr)
 50      CONTINUE
         CALL FTCLOS(lun,ierr)

         IF ( ierr .NE. 0 ) THEN
            wrtstr = 'Failed to get data from '//filion(:lenact(filion))
            CALL xwrite(wrtstr, 10)
            RETURN
         ENDIF
c
c	open (unit=80,access='direct',file=filrec,
c     *    recl=13*4,status='old',form='unformatted')
c	do i = 1,numion
c	  read (unit=80,rec=i) (arec(j,i),j=1,13)
c	enddo
c	close (unit=80)
 
         CALL ftopen(lun,filrec,0,index,ierr)
         CALL FTMRHD(lun,1,index,ierr)
 
         DO 100 i = 1 , NUMION
            CALL FTGCVE(lun,1,i,1,13,0.,AREc(1,i),qanyf,ierr)
 100     CONTINUE
         CALL FTCLOS(lun,ierr)

         IF ( ierr .NE. 0 ) THEN
            wrtstr = 'Failed to get data from '//filrec(:lenact(filrec))
            CALL xwrite(wrtstr, 10)
            RETURN
         ENDIF

c
c	open (unit=80,access='direct',file=filkar,status='old',
c     *    form='unformatted',recl=625*4)
c	read (80,rec=1) ((ga(i,j),j=1,nu),i=1,ng)
c	close (unit=80,status='keep')
 
         CALL ftopen(lun,filkar,0,index,ierr)
         CALL FTMRHD(lun,1,index,ierr)
 
         DO 150 i = 1 , NU
            CALL FTGCVE(lun,1,i,1,NG,0.,GA(1,i),qanyf,ierr)
 150     CONTINUE
         CALL FTCLOS(lun,ierr)

         IF ( ierr .NE. 0 ) THEN
            wrtstr = 'Failed to get data from '//filkar(:lenact(filkar))
            CALL xwrite(wrtstr, 10)
            RETURN
         ENDIF
 
c	open (unit=80,access='direct',file=filcon,status='old',
c     *    form='unformatted',recl=6*4)
c	do i = 1,numion
c	  read (80,rec=i) (p(j,i),j=1,6)
c	enddo
c	close (unit=80,status='keep')
 
         CALL ftopen(lun,filcon,0,index,ierr)
         CALL FTMRHD(lun,1,index,ierr)
 
         DO 200 i = 1 , NUMION
            CALL FTGCVE(lun,1,i,1,6,0.,P(1,i),qanyf,ierr)
 200     CONTINUE
         CALL FTCLOS(lun,ierr)

         IF ( ierr .NE. 0 ) THEN
            wrtstr = 'Failed to get data from '//filcon(:lenact(filcon))
            CALL xwrite(wrtstr, 10)
            RETURN
         ENDIF

c
c        open(unit=80,access='direct',file=fillin,recl=21*4,
c     *    status='old',form='unformatted')
c	read (80,rec=1) ne
c        do i = 1,ne
c          read (80,rec=i+1) nrl(i),
c     *        (lp(k,i),k=1,9),trans(i),(rp(k,i),k=1,14)
c        enddo
c        close(unit=80,status='keep')
 
         CALL ftopen(lun,fillin,0,index,ierr)
         CALL FTMRHD(lun,1,index,ierr)
 
         CALL FTGKYJ(lun,'NAXIS2',NE,comment,ierr)
 
         DO 250 i = 1 , NE
            CALL FTGCVJ(lun,1,i,1,1,0,NRL(i),qanyf,ierr)
            CALL FTGCVJ(lun,2,i,1,9,0,LP(1,i),qanyf,ierr)
            CALL FTGCVS(lun,3,i,1,1,' ',TRAns(i),qanyf,ierr)
            CALL FTGCVE(lun,4,i,1,14,0.,RP(1,i),qanyf,ierr)
 250     CONTINUE
         CALL FTCLOS(lun,ierr)

         IF ( ierr .NE. 0 ) THEN
            wrtstr = 'Failed to get data from '//fillin(:lenact(fillin))
            CALL xwrite(wrtstr, 10)
            RETURN
         ENDIF
 
c
c        open(unit=80,access='direct',file=filcaf,recl=2*4,
c     *    status='old',form='unformatted')
c	read (unit=80,rec=1) ncafe
c        do i = 1,ncafe
c          read (unit=80,rec=i+1) idnr(i),cdrcafe(i)
c        enddo
c        close(unit=80,status='keep')
 
         CALL ftopen(lun,filcaf,0,index,ierr)
         CALL FTMRHD(lun,1,index,ierr)
 
         CALL FTGKYJ(lun,'NAXIS2',NCAfe,comment,ierr)
 
         DO 300 i = 1 , NCAfe
            CALL FTGCVJ(lun,1,i,1,1,0,IDNr(i),qanyf,ierr)
            CALL FTGCVE(lun,2,i,1,1,0.,CDRcafe(i),qanyf,ierr)
 300     CONTINUE
         CALL FTCLOS(lun,ierr)

         IF ( ierr .NE. 0 ) THEN
            wrtstr = 'Failed to get data from '//filcaf(:lenact(filcaf))
            CALL xwrite(wrtstr, 10)
            RETURN
         ENDIF
c
      ENDIF
C
C-- END OF INITIALISATION
C
      DO 400 i = 1 , NOEL
         xe(i) = 10.**(abund(i)-12.)*Factor(i)
 400  CONTINUE
      CALL MEKAL1(T,xe,xzin,Ed)
      CALL MEKAL4(E1,E2,Flux,Nbin,T,xzin)
      elden = Ed*H
      CALL MEKAL5(E1(1),E2(Nbin),T,xzin,elden,elx,flx,nlx)
      nl = nlx
      DO 500 i = 1 , Nbin
         de = E2(i) - E1(i)
         DO WHILE ( nl.GT.0 )
            IF ( elx(nl).GT.E2(i) ) GOTO 500
            Flux(i) = Flux(i) + flx(nl)/de
            nl = nl - 1
         ENDDO
 500  CONTINUE
      cst = 2.53325E-3*H**2*Ed/SQRT(T)*Cem
      DO 600 i = 1 , Nbin
         Flux(i) = Flux(i)*cst
 600  CONTINUE
      RETURN
      END
**==mekal1.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL1(T,Xe,Xzin,Ed)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL alfa , Ed , fn , sion , so , T , Xe , Xzin
      INTEGER i , iel , j , k , nion , NIONMAX , NOEL
C*** End of declarations inserted by SPAG
      INCLUDE 'mekal.inc'
      PARAMETER (NOEL=15,NIONMAX=29)
      DIMENSION Xe(NOEL) , Xzin(NUMION) , alfa(NUMION) , sion(NUMION) , 
     &          fn(NIONMAX) , nion(NOEL)
      DATA nion/2 , 3 , 7 , 8 , 9 , 11 , 12 , 13 , 14 , 15 , 17 , 19 , 
     &     21 , 27 , 29/
      CALL MEKAL2(sion,T)
      CALL MEKAL3(alfa,T)
      k = 2
      Xzin(1) = 0.0
      Xzin(2) = 1.0
      DO 200 iel = 2 , NOEL
         fn(1) = 1.0
         so = 1.0
         DO 50 i = 2 , nion(iel)
            fn(i) = fn(i-1)*sion(k+i-1)/alfa(k+i)
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
      END
**==mekal2.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL2(Sion,T)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL auto , f1 , f2 , f4 , fb , FMEKAL10 , FMEKAL13 , fuy , g , 
     &     sa , sd , Sion , sk , T , t15 , tsq , x , y , yay , yfuy
      REAL YLIM
      INTEGER ia , ind , j , k , n , nj
C*** End of declarations inserted by SPAG
      INCLUDE 'mekal.inc'
      PARAMETER (YLIM=69.0)
      DIMENSION Sion(NUMION)
      t15 = T**1.5
      tsq = SQRT(T)
      ia = 0
      DO 200 k = 1 , NUMION
         Sion(k) = 0.
         ia = ia + 1
         nj = NINT(AIOn(ia))
         IF ( nj.NE.0 ) THEN
            sa = 0.
            DO 20 j = 1 , nj
               y = AIOn(ia+1)/T
               x = 1./y
               IF ( x.GT.0.055 ) THEN
                  f4 = (-5.725E-4+x*(0.01345+x*(0.8691+0.03404*x)))
     &                 /(1.+x*(2.197+x*(0.2454+2.053E-3*x)))
               ELSE
                  f4 = .7699*x**1.9496
               ENDIF
               fuy = FMEKAL13(y)
               IF ( y.LT.10.0 ) THEN
                  fb = 1. + y - y*fuy*(2.+y)
               ELSE
                  CALL MEKAL12(y,fb)
               ENDIF
               IF ( y.LT.YLIM ) sa = sa + EXP(-y)
     &                               *((AIOn(ia+2)*(1.-y*fuy)+AIOn(ia+3)
     &                               *fb+AIOn(ia+4)*fuy)/y+AIOn(ia+5)
     &                               *f4)
               ia = ia + 5
 20         CONTINUE
            Sion(k) = sa/t15
         ENDIF
         ia = ia + 1
         ind = NINT(AIOn(ia))
         IF ( ind.EQ.1 ) THEN
         ELSEIF ( ind.EQ.3 ) THEN
            y = AIOn(ia+1)/T
            IF ( y.LT.YLIM ) THEN
               fuy = FMEKAL13(y)
               yfuy = y*fuy
               sa = EXP(-y)
     &              *(AIOn(ia+2)*fuy+AIOn(ia+3)*(1.0-yfuy)+AIOn(ia+4)
     &              *yfuy+AIOn(ia+5)*y*(1.0-yfuy))
               Sion(k) = Sion(k) + sa/tsq
            ENDIF
            ia = ia + 5
         ELSEIF ( ind.EQ.4 ) THEN
            y = AIOn(ia+1)/T
            IF ( y.LT.YLIM ) Sion(k) = Sion(k) + AIOn(ia+2)*EXP(-y)
     &                                 *(1.0-y*FMEKAL13(y))/tsq
            ia = ia + 2
         ELSEIF ( ind.EQ.5 ) THEN
            sd = 0.
            DO 40 n = 1 , 18
               y = AIOn(ia+1)/T
               IF ( y.LT.YLIM ) THEN
                  yay = y*AIOn(ia+2)
                  fuy = FMEKAL13(yay)
                  sd = sd + EXP(-y)
     &                 *(FMEKAL13(y)*AIOn(ia+3)+AIOn(ia+4)+y*(AIOn(ia+5)
     &                 *fuy+AIOn(ia+6)*(1.-yay*fuy)))
               ENDIF
               ia = ia + 6
 40         CONTINUE
            Sion(k) = Sion(k) + sd/tsq
         ELSEIF ( ind.EQ.6 ) THEN
            sk = 0.
            DO 60 j = 1 , 12
               y = AIOn(ia+1)/T
               IF ( y.LT.YLIM ) THEN
                  f1 = EXP(-y)*(y+1.)
                  y = AIOn(ia+2)/T
                  IF ( y.LT.YLIM ) THEN
                     f2 = EXP(-y)*(y+1.)
                     sk = sk + (f1-f2)*AIOn(ia+3)
                  ENDIF
               ENDIF
               ia = ia + 3
 60         CONTINUE
            Sion(k) = Sion(k) + sk*tsq
            GOTO 100
         ELSEIF ( ind.EQ.7 ) THEN
            GOTO 100
         ELSE
            y = AIOn(ia+1)/T
            IF ( y.LT.YLIM ) THEN
               auto = AIOn(ia+2)*EXP(-y)*(1.+AIOn(ia+3)*FMEKAL13(y))/tsq
               Sion(k) = Sion(k) + auto
            ENDIF
            ia = ia + 3
         ENDIF
         GOTO 200
 100     y = AIOn(ia+1)/T
         g = FMEKAL10(y,AIOn(ia+2),AIOn(ia+3),AIOn(ia+4),AIOn(ia+5),
     &       AIOn(ia+6))
         IF ( y.LT.YLIM ) Sion(k) = Sion(k) + g*EXP(-y)/tsq
         ia = ia + 6
 200  CONTINUE
      IF ( ia.GT.NAMAX ) WRITE (*,*) 'INCOMPATIBLE DATASET'
      RETURN
      END
**==mekal3.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL3(Alfa,T)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL Alfa , diel , T , t13 , t15 , t25 , tlo , tsq , y , YLIM
      INTEGER ind , ip , k , l , m
C*** End of declarations inserted by SPAG
      INCLUDE 'mekal.inc'
      PARAMETER (YLIM=69.0)
      DIMENSION Alfa(NUMION)
      tlo = LOG(T)
      tsq = SQRT(T)
      t13 = T**0.33333333
      t15 = T**1.5
      t25 = T**2.5
      DO 100 k = 1 , NUMION
         ind = NINT(AREc(1,k))
         IF ( ind.GT.0 ) THEN
            IF ( ind.LE.8 ) THEN
               Alfa(k) = AREc(2,k)/T**AREc(3,k)
            ELSE
               Alfa(k) = AREc(2,k)/T**(AREc(3,k)+AREc(4,k)*T)
            ENDIF
            IF ( ind.EQ.2 ) THEN
               IF ( T.GT.AREc(4,k) ) THEN
                  Alfa(k) = Alfa(k) + AREc(5,k)/T**AREc(6,k)
               ELSE
                  y = AREc(7,k)/T
                  IF ( y.LT.YLIM ) THEN
                     diel = AREc(8,k)*EXP(-y)/t15
                     y = AREc(9,k)/T
                     IF ( y.LT.YLIM )
     &                    diel = diel*(1.+AREc(10,k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  ENDIF
               ENDIF
            ELSEIF ( ind.EQ.3 ) THEN
               y = AREc(4,k)/T
               IF ( y.LT.YLIM ) THEN
                  diel = AREc(5,k)*EXP(-y)/t15
                  y = AREc(6,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+AREc(7,k)*EXP(-y))
     &                 **AREc(8,k)
                  y = AREc(9,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+AREc(10,k)*EXP(-y))
     &                 **AREc(11,k)
                  Alfa(k) = Alfa(k) + diel
               ENDIF
            ELSEIF ( ind.EQ.4 ) THEN
               IF ( T.GT.AREc(4,k) ) THEN
                  y = AREc(5,k)/T
                  IF ( y.LT.YLIM ) THEN
                     diel = AREc(6,k)*EXP(-y)/t15
                     y = AREc(7,k)/T
                     IF ( y.LT.YLIM ) diel = diel*(1.+AREc(8,k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  ENDIF
               ELSE
                  Alfa(k) = Alfa(k) + AREc(9,k)/T**AREc(10,k)
               ENDIF
            ELSEIF ( ind.EQ.5 ) THEN
               IF ( T.GT.AREc(4,k) ) THEN
                  y = AREc(5,k)/T
                  IF ( y.LT.YLIM ) THEN
                     diel = AREc(6,k)*EXP(-y)/t15
                     y = AREc(7,k)/T
                     IF ( y.LT.YLIM ) diel = diel*(1.+AREc(8,k)*EXP(-y))
                     Alfa(k) = Alfa(k) + diel
                  ENDIF
               ELSE
                  y = AREc(9,k)/T
                  IF ( y.LT.YLIM ) Alfa(k) = Alfa(k) + AREc(10,k)
     &                 /t15*EXP(-y)
               ENDIF
            ELSEIF ( ind.EQ.6 ) THEN
               y = AREc(4,k)/T
               IF ( y.LT.YLIM ) THEN
                  diel = AREc(5,k)*EXP(-y)/t15
                  y = AREc(6,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+AREc(7,k)*EXP(-y))
                  Alfa(k) = Alfa(k) + diel
               ENDIF
               IF ( T.LT.AREc(8,k) ) THEN
                  y = AREc(9,k)/T
                  IF ( y.LT.YLIM ) Alfa(k) = Alfa(k) + AREc(10,k)
     &                 *EXP(-y)
     &                 /t25*(((T+AREc(11,k))*T+AREc(12,k))*T+AREc(13,k))
               ENDIF
            ELSEIF ( ind.EQ.7 ) THEN
               y = AREc(4,k)/T
               IF ( y.LT.YLIM ) THEN
                  diel = AREc(5,k)*EXP(-y)/t15
                  y = AREc(6,k)/T
                  IF ( y.LT.YLIM ) diel = diel*(1.+AREc(7,k)*EXP(-y))
                  Alfa(k) = Alfa(k) + diel
               ENDIF
            ELSEIF ( ind.EQ.8 .OR. ind.EQ.9 ) THEN
               ip = 4
               IF ( ind.EQ.9 ) ip = ip + 1
               DO 10 l = 1 , 4
                  m = ip + 2*l - 2
                  IF ( AREc(m+1,k).NE.0. ) THEN
                     y = AREc(m,k)/T
                     IF ( y.LT.YLIM ) Alfa(k) = Alfa(k) + AREc(m+1,k)
     &                    /t15*EXP(-y)
                  ENDIF
 10            CONTINUE
            ELSE
               y = AREc(4,k)/T
               IF ( y.LT.YLIM ) Alfa(k) = Alfa(k) + AREc(5,k)
     &              /t15*EXP(-y)
            ENDIF
         ELSEIF ( ind.EQ.0 ) THEN
            Alfa(k) = (AREc(2,k)-AREc(3,k)*tlo+AREc(4,k)*t13)/tsq
         ELSE
            Alfa(k) = 0.
         ENDIF
 100  CONTINUE
      END
**==mekal4.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL4(E1,E2,Flx,Nemx,T,Xzin)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL a , a31 , a32 , bf , CMAX , dgfbk , dx , E1 , E2 , e2p , ej , 
     &     f2p , Flx , FMEKAL10 , g , g2a , g2p , gfb , gfbk , gff
      REAL gl , phi , psi , T , twz , u , ua , ul , UMAX , w , W3PI , 
     &     x , x50 , xez , xzg , Xzin , y
      INTEGER i , i1 , i2 , iel , ifo , in , in6 , inmax , ip , iu , j , 
     &        j1 , j2 , j3 , ju , k , k1 , n , n0 , n2
      INTEGER NA , NBF , nel , Nemx , nk , NO2P , NOEL , NP
C*** End of declarations inserted by SPAG
      INCLUDE 'mekal.inc'
      PARAMETER (NOEL=15,NO2P=2)
      PARAMETER (NP=51,NBF=15,NA=NU+5+6*NBF)
      PARAMETER (UMAX=69.0,CMAX=1.E-3,W3PI=0.551328895)
      DIMENSION E1(Nemx) , E2(Nemx) , Flx(Nemx) , Xzin(NUMION) , 
     &          a(NA,NUMION)
      DIMENSION g2a(NG) , ua(NU) , bf(4,NBF) , phi(NP,NO2P)
      DIMENSION nel(NOEL) , e2p(NOEL,NO2P) , f2p(NOEL,NO2P)
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
      DATA nel/2 , 5 , 12 , 20 , 29 , 40 , 52 , 65 , 79 , 94 , 111 , 
     &     130 , 151 , 178 , 207/
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
      DO 200 i = 1 , NU
         ua(i) = -4.25 + FLOAT(i)*0.25
 200  CONTINUE
      iel = 1
      k = 0
      DO 300 i = 1 , NUMION
         k1 = k + 1
         IF ( k1.LT.NUMION ) THEN
            DO 220 j = 1 , NA
               a(j,k1) = 0.
 220        CONTINUE
         ENDIF
         ifo = 0
         IF ( Xzin(i).GT.0. ) THEN
            xez = Xzin(i)*P(1,i)
            IF ( xez.GT.CMAX ) THEN
               ifo = 1
               k = k + 1
               a(1,k) = xez
               gl = ALOG10(.0136*P(1,i)/T)
               i1 = INT(4.*(gl+4.)) + 1
               IF ( i1.LT.1 ) i1 = 1
               IF ( i1.GE.NG ) i1 = NG - 1
               i2 = i1 + 1
               w = (g2a(i2)-gl)/(g2a(i1)-g2a(i2))
               DO 230 ju = 1 , NU
                  a(ju+1,k) = GA(i2,ju) + w*(GA(i2,ju)-GA(i1,ju))
 230           CONTINUE
            ENDIF
            DO 240 i2 = 1 , NO2P
               IF ( i.EQ.nel(iel)-i2 ) THEN
                  IF ( i2.NE.1 ) THEN
                     g = 0.05
                  ELSEIF ( iel.EQ.1 ) THEN
                     g = FMEKAL10(e2p(iel,i2)/T,0.08,-0.16,0.11,0.0,0.0)
                  ELSE
                     g = 0.055
                  ENDIF
                  twz = 164995.*f2p(iel,i2)*g*Xzin(i)/e2p(iel,i2)
                  IF ( twz.GT.CMAX ) THEN
                     u = e2p(iel,i2)/T
                     IF ( u.LT.UMAX ) THEN
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
 240        CONTINUE
            IF ( i.EQ.nel(iel) ) iel = iel + 1
            n0 = NINT(P(2,i))
            in = 0
            IF ( n0.GT.0 ) THEN
               DO 250 n = n0 , NBF
                  n2 = n*n
                  IF ( n.EQ.n0 ) THEN
                     a31 = bf(1,n)*P(3,i)*Xzin(i)/T
                     a32 = P(4,i)/T
                  ELSE
                     a31 = bf(1,n)*P(5,i)*Xzin(i)/T/FLOAT(n2*n)
                     a32 = P(6,i)/T/FLOAT(n2)
                  ENDIF
                  IF ( a32.LT.UMAX ) THEN
                     IF ( a31*EXP(a32).GT.CMAX ) THEN
                        IF ( ifo.EQ.0 ) THEN
                           ifo = 1
                           k = k + 1
                        ENDIF
                        in6 = 6*in
                        a(31+in6,k) = a31
                        a(32+in6,k) = a32
                        a(33+in6,k) = a32*FLOAT(n2)
                        a(34+in6,k) = bf(2,n)
                        a(35+in6,k) = bf(3,n)
                        a(36+in6,k) = bf(4,n)
                        in = in + 1
                     ENDIF
                  ENDIF
 250           CONTINUE
            ENDIF
            IF ( ifo.EQ.1 ) a(30,k) = in
         ENDIF
 300  CONTINUE
      nk = k
      DO 400 j = 1 , Nemx
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
            DO 320 k = 1 , nk
               IF ( a(1,k).GT.0. ) THEN
                  IF ( u.LT.1.E-4 ) THEN
                     xzg = a(1,k)*(EXP(a(2,k))-W3PI*LOG(u/1.E-4))
                  ELSE
                     xzg = a(1,k)*EXP(a(j3,k)+w*(a(j3,k)-a(j2,k)))
                  ENDIF
                  gff = gff + xzg
               ENDIF
               gfbk = 0.
               in = 0
               inmax = NINT(a(30,k))
               dgfbk = 1.
               DO WHILE ( in.LT.inmax .AND. dgfbk.GT.CMAX*gfbk )
                  in6 = 6*in
                  IF ( a(31+in6,k).GT.0 ) THEN
                     y = u - a(32+in6,k)
                     IF ( y.GE.0.0 .AND. y.LT.UMAX ) THEN
                        x = 1./(1.+SQRT(u/a(33+in6,k)))
                        dgfbk = (((a(36+in6,k)*x+a(35+in6,k))*x+a(34+in6
     &                          ,k))*x+1.)*x*a(31+in6,k)*EXP(-y)
                        gfbk = gfbk + dgfbk
                     ENDIF
                  ENDIF
                  in = in + 1
               ENDDO
               gfb = gfb + gfbk
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
 320        CONTINUE
            gff = gff*EXP(-u)
         ENDIF
         Flx(j) = (gff+gfb+g2p)/ej
 400  CONTINUE
      RETURN
      END
**==mekal5.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL5(E1,E2,Xkt,Xzin,Elden,Elx,Flx,Nlx)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL a , dcor , E1 , E2 , el , Elden , Elx , eta , f , fg , fgc , 
     &     Flx , FMEKAL10 , g , t , tau , Xkt , xz , Xzin , y
      INTEGER i , lp3 , Nlx , nzz
C*** End of declarations inserted by SPAG
      INCLUDE 'mekal.inc'
      DIMENSION Xzin(*) , Elx(*) , Flx(*)
      INTEGER iaux
      t = Xkt*1.16048E+7
      tau = Xkt*1.16048
      eta = Elden*1E-12
      Nlx = 0
      DO 100 i = 1 , NE
         el = RP(1,i)
         IF ( el.LT.E1 ) RETURN
         IF ( el.LT.E2 ) THEN
            nzz = LP(4,i)
            xz = Xzin(nzz)
            IF ( xz.GT.0. ) THEN
               a = RP(3,i)
               y = el*a/Xkt
               IF ( y.LT.40. ) THEN
                  g = FMEKAL10(y,RP(4,i),RP(5,i),RP(6,i),RP(7,i),RP(8,i)
     &                )
                  f = RP(2,i)
                  iaux = NRL(i)
                  CALL MEKAL6(f,g,a,fg,el,y,t,Xzin,nzz,RP(1,i),LP(1,i),
     &                        TRAns(i),Elden,iaux,IDNr,CDRcafe,NCAfe)
                  fgc = EXP(-y)*xz*1.646E+5
                  IF ( eta.GT.1E-10 ) THEN
                     IF ( ABS(RP(11,i)).GT.1E-10 .OR. ABS(RP(14,i))
     &                    .GT.1E-10 ) THEN
                        lp3 = LP(3,i)
                        CALL MEKAL9(dcor,eta,tau,RP(1,i),lp3)
                        fgc = fgc*dcor
                     ENDIF
                  ENDIF
                  Nlx = Nlx + 1
                  Elx(Nlx) = el
                  Flx(Nlx) = fg*fgc/el
               ENDIF
            ENDIF
         ENDIF
 100  CONTINUE
      RETURN
      END
**==mekal6.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL6(F,G,A,Fg,El,Y,T,Xzin,Nzz,Rp,Lp,Trans,Elden,Lnum,
     &                  Idnr,Cdrcafe,Ncafe)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL A , b6 , bdr , bii , brf , bri , c1 , c1f , cdr , Cdrcafe , 
     &     cii , ciid , dr , drs , El , Elden , ex , ex1 , ex2 , F
      REAL f5 , f6 , Fg , FMEKAL11 , FMEKAL7 , FMEKAL8 , fn , fr , G , 
     &     rii , Rp , rr , T , t2 , t5 , x , xcdr , Xzin , xzn , xzn1
      REAL Y , yy , z , z12 , z4
      INTEGER ica , idr , idrs , iel , iex , ife , ihel , ii , img ,
     &        isat , itr , iz , izion , Lnum , Ncafe , Nzz
C*** End of declarations inserted by SPAG
      INTEGER Lp(*)
      INTEGER Idnr(Ncafe)
      DIMENSION Xzin(*) , Rp(*) , Cdrcafe(Ncafe)
      DIMENSION izion(15) , cdr(19) , bdr(19) , b6(19) , ciid(19) , 
     &          bri(15)
      character(1) aster , dum1str
      character(8) Trans
      DATA img , ica , ife/8 , 13 , 14/
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
      IF ( iex.NE.0 .AND. iex.NE.2 ) THEN
         yy = (Y+0.6)/(Y+1.15)
         c1f = 1. + 7.*yy*(1.065*(1-bri(iel))+0.4*EXP(-0.21*Y))
         brf = FMEKAL7(z,t5,Elden,bri(iel))
         IF ( iex.EQ.3 ) THEN
            ex = (bri(iel)+c1f*(1.-brf)/yy/7.46)*1.065
         ELSE
            ex = c1f*brf
         ENDIF
      ELSE
         ex = 1.
      ENDIF
      ex = ex*Fg
      IF ( idrs.NE.0 ) THEN
         iz = iz + 1
         z12 = (iz+1)**2
         x = (El*A)/(iz+1)
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
         IF ( idr.NE.0 ) THEN
            IF ( idr.EQ.1 ) THEN
               IF ( iel.EQ.img .OR. iel.EQ.ica .OR. iel.EQ.ife ) THEN
                  IF ( dum1str.EQ.aster ) THEN
                     DO 2 isat = 1 , Ncafe
                        IF ( Idnr(isat).EQ.Lnum ) GOTO 4
 2                   CONTINUE
 4                   xcdr = Cdrcafe(isat)
                     GOTO 10
                  ENDIF
               ENDIF
               xcdr = cdr(itr)
 10            dr = xcdr/z4*EXP(bdr(itr)*Y)/(b6(itr)+1./(z-1.)**4)
               rr = Rp(9)*T**Rp(10)
            ELSEIF ( idr.EQ.2 ) THEN
               dr = 11.*EXP(0.284*Y)/(1.+z4*6E-6) + 27.*EXP(0.1*Y)
     &              /(1.+z4*3E-5)
               rr = Rp(9)*T**Rp(10)
            ELSE
               ex1 = EXP(0.284*Y)
               ex2 = EXP(0.1*Y)
               f5 = 2.3*ex1 + 215.*ex2/(1.+z4*1.9E-4)
               f6 = 16.*ex1/(1.+z4*7E-5) + 90.*ex2/(1.+z4*8E-5)
               brf = FMEKAL7(z,t5,Elden,bri(iel))
               dr = FMEKAL8(bri(iel),brf,f5,f6,idr+2)
               f5 = (iz+1.)**2.62/T**0.31*18.72E-4
               f6 = (iz+1.)**2.08/T**0.04*.575E-4
               rr = FMEKAL8(bri(iel),brf,f5,f6,idr+2)
            ENDIF
            fr = xzn1/xzn*El/12.3985
            rr = fr*EXP(Y)*rr
            dr = fr*0.47*z4/T*dr
         ELSE
            rr = xzn1/xzn*EXP(Y)*El/12.3985*Rp(9)*T**Rp(10)
         ENDIF
      ENDIF
      rii = 0.
      IF ( ii.NE.0 ) THEN
         xzn1 = Xzin(Nzz-1)
         xzn = Xzin(Nzz)
         IF ( xzn1.GT.0 ) THEN
            IF ( ii.EQ.2 ) THEN
               IF ( itr.NE.0 ) THEN
                  cii = ciid(itr)
                  IF ( itr.GE.7 .AND. itr.LE.14 ) THEN
                     bii = 1.242 + 0.009*(z-26)
                  ELSE
                     bii = 1.296
                  ENDIF
               ELSE
                  bii = 1.164 + 0.008*(z-26)
                  cii = 0.3
               ENDIF
               t2 = cii*FMEKAL11(bii*Y)*1.231/bii/(1.+2.E5*z**(-3.5))
            ELSE
               yy = (z-2.)**4.3*6E12
               fn = yy + 1.33333333*Elden
               f5 = Elden/fn
               f6 = (yy+Elden/3.)/fn
               brf = FMEKAL7(z,t5,Elden,bri(iel))
               ihel = 5
               IF ( ii.EQ.1 ) ihel = 6
               t2 = 0.22*FMEKAL11(1.33*Y)
     &              *FMEKAL8(bri(iel),brf,f5,f6,ihel)
            ENDIF
            rii = xzn1/xzn*EXP(Y)*t2
         ENDIF
      ENDIF
      Fg = ex + drs + dr + rii + rr
      RETURN
      END
**==fmekal7.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      FUNCTION FMEKAL7(Z,T5,Elden,Bri)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL betac , Bri , Elden , ex , FMEKAL7 , T5 , Z
C*** End of declarations inserted by SPAG
      ex = -0.73/Z**0.415
      betac = 330./Z**14.57*(T5/Z**2)**ex
      FMEKAL7 = 1./(1.+(betac*Elden)*Bri)
      RETURN
      END
**==fmekal8.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      FUNCTION FMEKAL8(Bri,Brf,F5,F6,Itr)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL Brf , Bri , F5 , F6 , FMEKAL8
      INTEGER Itr
C*** End of declarations inserted by SPAG
      IF ( Itr.NE.6 ) THEN
         FMEKAL8 = F5*Bri + (F6+(1.-Bri)*F5)*(1.-Brf)
      ELSE
         FMEKAL8 = (F6+(1.-Bri)*F5)*Brf
      ENDIF
      RETURN
      END
**==mekal9.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL9(Dcor,Eta,Tau,Rp,Lp3)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL Dcor , Eta , ex , Rp , Tau
      INTEGER Lp3
C*** End of declarations inserted by SPAG
      DIMENSION Rp(*)
      IF ( Lp3.EQ.5 ) THEN
         Dcor = 1. + Rp(11)/(Eta**0.03+120./Tau**0.3/Eta)
      ELSEIF ( ABS(Rp(14)).GT.1E-10 ) THEN
         Dcor = 1. + Rp(14)/(1.+3.5/Eta**0.75)
      ELSE
         ex = -Rp(12)*Eta**Rp(13)
         Dcor = 1. + Rp(11)*(1.-EXP(ex))
      ENDIF
      RETURN
      END
**==fmekal10.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      FUNCTION FMEKAL10(Y,A,B,C,D,E)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL A , B , C , D , E , e1 , f2 , f3 , f4 , f5 , FMEKAL10 , Y
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
      FMEKAL10 = A + B*f2 + C*f3 + D*f4 + E*f5
      RETURN
      END
**==fmekal11.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      FUNCTION FMEKAL11(X)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL a , FMEKAL11 , X
C*** End of declarations inserted by SPAG
      DIMENSION a(6)
      DATA a/ - .57721566 , .99999193 , -.24991055 , .05519968 , 
     &     -.00976004 , .00107857/
      IF ( X.LE.1. ) THEN
         FMEKAL11 = -ALOG(X) + ((((a(6)*X+a(5))*X+a(4))*X+a(3))*X+a(2))
     &              *X + a(1)
      ELSE
         FMEKAL11 = (.250621/X+X+2.334733)*EXP(-X)
     &              /(X*X+3.330657*X+1.681534)
      ENDIF
      RETURN
      END
**==mekal12.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      SUBROUTINE MEKAL12(Y,Fb)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL c , Fb , ff , sb , t , Y
      INTEGER i
C*** End of declarations inserted by SPAG
C=====================================
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
**==fmekal13.spg  processed by SPAG 4.50J  at 14:35 on 16 May 1995
 
      FUNCTION FMEKAL13(U)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      REAL FMEKAL13 , U
C*** End of declarations inserted by SPAG
      IF ( U.LT.1.0 ) THEN
         FMEKAL13 = LOG(1.+1./U) - (0.36+0.03/SQRT(U+0.01))/(1.+U)**2
      ELSE
         FMEKAL13 = LOG(1.+1./U) - (0.36+0.03*SQRT(U+0.01))/(1.+U)**2
      ENDIF
      RETURN
      END
 
