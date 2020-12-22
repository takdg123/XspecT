!*==compbb.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
      SUBROUTINE COMPBB(Ear,Ne,Param,Ifl,Photar,Photer)
      IMPLICIT NONE
 
      INTEGER Ne , Ifl
      REAL Ear(0:Ne) , Param(3) , Photar(Ne) , Photer(Ne)
 
!     Param(1) = Tbb
!     Param(2) = Te
!     Param(3) = tau
 
!     Originally coded by Kazu Mitsuda for
!     ISAS SPFD on FACOM mainframe  Ported to xspec format
!     by Ken Ebisawa. The architecture specific
!     numerical library (FUJITSU SSL) was removed.
!
!     Last modified on 1999-07-21.

!     Modification on 2004-02-24 by Ken Ebisawa
!     The energy range of the model calculation were hard coded in
!     the previous version as between 0.1 keV and 70 keV. 
!     (This was adequate for Tenma and Ginga!)
!     Now, these values are taken from Ear(0) and Ear(Ne), namely,
!     lower and upper boundaries of the input energy array.
 
      DOUBLE PRECISION fl_low , fl_up, E0, Em
      INTEGER i

! included to suppress compilation warning
      i = ifl

! this model has no errors
      DO i = 1, Ne
         Photer(i) = 0.0
      ENDDO

      E0 = dble(ear(0))
      Em = dble(ear(ne))

      CALL CMPBBK(DBLE(Param(1)),DBLE(Param(2)),DBLE(Param(3)), &
     &            DBLE(Ear(0)),fl_low,E0,Em)
      DO i = 1 , Ne
         CALL CMPBBK(DBLE(Param(1)),DBLE(Param(2)),DBLE(Param(3)), &
     &               DBLE(Ear(i)),fl_up,E0,Em)
         Photar(i) = SNGL((fl_low+fl_up)*(Ear(i)-Ear(i-1))/2.0d0)
         fl_low = fl_up
      ENDDO
 
      END
!*==cmpbbk.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      SUBROUTINE CMPBBK(Ts,Te,Tau,E,Fl,E0,Em)
 
      IMPLICIT NONE
      DOUBLE PRECISION Ts , Te , Tau , E , Fl
      DOUBLE PRECISION FLUxs(0:511) 

      INTEGER N2
      DOUBLE PRECISION E0 , DX, Em
 
      INTEGER j
      INTEGER icon , i
      DOUBLE PRECISION x(0:511) , c(512) , v
      DOUBLE PRECISION ts1 , te1 , tau1, e01, em1

      LOGICAL qnewe

      SAVE ts1, te1, tau1
      SAVE x, c, FLUxs
      SAVE DX, N2

      DATA ts1/0.0/ , te1/0.0/ , tau1/ - 0.01/
      DATA e01/-1.0/, em1/-1.0/


      IF ( Ts.NE.ts1 .OR. Te.NE.te1 .OR. Tau.NE.tau1 .OR. &
     &     E0.NE.e01 .OR. Em.NE.em1 ) THEN
         IF ( Tau.NE.tau1 ) CALL SETSPB(Tau)
         qnewe = .FALSE.
         IF ( E0 .NE. e01 .OR. Em .NE. em1 ) qnewe = .TRUE.
         CALL SCATS(Ts,Te,ts1,te1,FLUxs,E0,DX,N2,Em,qnewe)

         ts1 = Ts
         te1 = Te
         tau1 = Tau
         e01 = E0
         em1 = Em
         DO j = 0 , N2 - 1
            x(j) = DX*j
         ENDDO
         CALL DBIC3(FLUxs,N2,c,icon)
         IF ( icon.NE.0 ) WRITE (6,*) 'DBIC3 RETURN CODE:' , icon
      ENDIF

      v = LOG(E/E0)
      i = INT(v/DX) + 1
      IF ( i.LT.1 .OR. i.GE.N2 ) THEN
         Fl = 0.0
      ELSE
         CALL DBIF3(x,N2,c,v,i,Fl,icon)
         IF ( icon.NE.0 ) WRITE (6,*) 'DBIF3 RETURN CODE:' , icon
      ENDIF

!d
!d      WRITE(6,*) 'E,V,I,FL=',E,V,I,FL
!d
 
      RETURN
      END
!*==scats.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
 
      SUBROUTINE SCATS(Ts,Te,Ts1,Te1,FLUxs,E0,DX,N2,Em,Qnewe)
      IMPLICIT NONE
      DOUBLE PRECISION Ts , Te , Ts1 , Te1, E0, DX
      DOUBLE PRECISION FLUxs(0:511)
      INTEGER N2
      LOGICAL Qnewe
 
!     SUBROUTINE TO CALCULTE COMPTONIZED BLACK BODY SPECTRUM
!     BY FFT.  THE FORMULA FOR THE ENERGY DISTRIBUTION OF THE ONCE-
!     SCATTERED PHOTON IS GIVEN BY NEGLECTING THE DEPENDENCE ON THE
!     SCATTERING ANGLE ON THE ELECTRON REST FRAME.
!     THE DITRIBUTION IS AVERAGED FOR ALL INCIDENT AND SCATTERING ANGLE
!     OF THE PHOTON.  FOR THE ONCE SCATTERED PHOTON THE DISTRIBUTION
!     FOR THE FRONT SCATTERING IS USED. THIS IS VALID WHEN THE ORIGINAL
!     PHOTON COMES FROM THE BOTTOM OF THE HOT GAS, SO THE MOST OF ONCE
!     SCATTERED PHOTON EXPERINCED FRONT SCATTERING.
!     ESCAPE PROBABILTY OF N-TIME SCATTERED PHOTON
!     SHOULD BE GIVEN BY THE SUBROUTINE SETSPB STORED IN
!     'SBSG010.SCAT.FORT77(PROB)'
!     THE POWER LAW DEPENDENCE OF PROBABILTY ON TIMES OF SCATTERING
!     N IS ASSUMED WHEN N IS LARGER THAN AN CERTAIN VALUE WHICH DEPENDS
!     ON TAU.
!     BY K.MITSUDA
 
!     <<< CAUTION >>>
!     THE NAMES OF THE SUBROUTINES AND COMMONS ARE ALMOST THE SAME AS THOSE
!     OF THE MSCATFFT, SO THIS PROGRAM CANNOT BE USED AT THE SAME TIME WITH
!     THE PROGRAM 'MSCATFFT'.
 
!     OUTPUT SPECTRUM IS STORED IN FLUXS(0:511)
!     FOR THE NORMALIZATION OF 'FLUXS', PLEASE SEE SUBROUTINE RESFFT.
 
      COMMON /PROBAB/ SPB
      DOUBLE PRECISION SPB(0:50)
      COMMON /PROBAA/ RATspb
      DOUBLE PRECISION RATspb
      COMMON /MAXTIM/ FMAx
      INTEGER FMAx
      DOUBLE PRECISION GN , FLUx0(0:2047) , F0R(0:2047) , F0I(0:2047) , &
     &       GFR(0:2047) , GFI(0:2047) , FFR(0:2047) , &
     &       FFI(0:2047) , FFR1(0:2047) , FFI1(0:2047)
      DOUBLE PRECISION rr , ii
      INTEGER N , N1
 
!      DOUBLE PRECISION    FNR(0:2047),FNI(0:2047),W,EPS/0.003/,NORM
      DOUBLE PRECISION fnr(0:2047) , fni(0:2047) , w , norm
      DOUBLE PRECISION xm , em

!      INTEGER J,F,CONVFG
      INTEGER j , f
      INTEGER icon

      SAVE GN , FLUx0 , F0R , F0I , GFR , GFI , FFR , FFI , FFR1 , FFI1
      SAVE N
 
!      N=2048
      N = 512
!      N2=512
      N2 = 128

      xm = LOG(em/E0)
      DX = xm/N2
      N1 = N - 1
 
      CALL RESFFT(Ts,Te,Ts1,Te1,GN,FLUx0,F0R,F0I,GFR,GFI,FFR,FFI,FFR1, &
     &            FFI1,E0,DX,Qnewe,N2,N1,N)
 
 
      norm = GN*SPB(1)
      DO j = 0 , N1
         fnr(j) = FFR1(j)*norm
         fni(j) = FFI1(j)*norm
      ENDDO

      DO f = 2 , FMAx
         DO j = 0 , N1
            w = GFR(j)*FFR(j) - GFI(j)*FFI(j)
            FFI(j) = GFI(j)*FFR(j) + GFR(j)*FFI(j)
            FFR(j) = w
         ENDDO
         norm = GN**f*SPB(f)
         IF ( f.LT.FMAx ) THEN
            DO j = 0 , N1
               fnr(j) = fnr(j) + FFR(j)*norm
               fni(j) = fni(j) + FFI(j)*norm
            ENDDO
         ELSE
            DO j = 0 , N1
               rr = 1.0D0 - GFR(j)*RATspb*GN
               ii = -GFI(j)*RATspb*GN
               norm = GN**f*SPB(f)/(rr*rr+ii*ii)
               fnr(j) = (fnr(j)+(FFR(j)*rr+FFI(j)*ii)*norm)/N
               fni(j) = (fni(j)+(FFI(j)*rr-FFR(j)*ii)*norm)/N
            ENDDO
         ENDIF
      ENDDO
 
!     ... FFT:  FNR,FNI (input) => FNR,FNI (output)
!     real imag.         real imag.
      CALL DCFTN(fnr,fni,N,-1,icon)
      IF ( icon.NE.0 ) THEN
         WRITE (6,*) 'DCFTN (FNR,FNI,)  RETURN CODE:' , icon
         STOP
      ENDIF
      CALL DPNR(fnr,fni,N,-1,icon)
!     ... FFT output.  with FNR and FNI in order of increasing wave number

      DO j = 0 , N2 - 1
         FLUxs(j) = FLUx0(j)*SPB(0) + fnr(j)
      ENDDO
      RETURN
!      DEBUG SUBCHK
      END
!*==resfft.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      SUBROUTINE RESFFT(Ts,Te,Ts1,Te1,GN,FLUx0,F0R,F0I,GFR,GFI,FFR,FFI, &
     &                  FFR1,FFI1,E0,DX,Qnewe,N2,N1,N)
      IMPLICIT NONE
      DOUBLE PRECISION Ts , Te , Ts1 , Te1
      DOUBLE PRECISION GN , FLUx0(0:2047) , F0R(0:2047) , F0I(0:2047) , &
     &       GFR(0:2047) , GFI(0:2047) , FFR(0:2047) , &
     &       FFI(0:2047) , FFR1(0:2047) , FFI1(0:2047)
      INTEGER N , N1 , N2
      DOUBLE PRECISION E0 , DX
      LOGICAL Qnewe

      DOUBLE PRECISION mcc, pi, ONCEF
      PARAMETER(mcc=511.0D0, pi=3.14159265358979D0)
      PARAMETER (ONCEF=7.0D0/16.0D0)
 
      INTEGER j , icon
      DOUBLE PRECISION e , x , tm, w
      DOUBLE PRECISION gfr1(0:2047) , gfi1(0:2047)
      DOUBLE PRECISION tm2 , stm2 , spi , z , ae , norm , itm2 , istm2
      DOUBLE PRECISION ERF_CMPBB , ERFC_COMPBB

      SAVE gfr1, gfi1

 
      IF ( Te1.NE.Te .OR. Qnewe ) THEN
         GN = DX
         Te1 = Te
         tm = Te/mcc
         tm2 = tm*2.0D0
         stm2 = SQRT(tm2)
         spi = SQRT(pi)
         istm2 = 1.0D0/stm2
         itm2 = 1.0D0/tm2
         norm = ERF_CMPBB(istm2) - 2.0D0*EXP(-itm2)/stm2/spi
         DO j = 0 , N1
            IF ( j.LT.N/2 ) THEN
               x = DX*j
            ELSE
               x = DX*(j-N)
            ENDIF
            e = EXP(x)
            ae = ABS(e-1)/(e+1)
            z = ae/stm2
!     ... BETA = 0->1
            GFR(j) = (e+1.0D0)/(2.0D0*stm2) &
     &               *((EXP(-z*z)*(1.0D0-tm2)-(ae-tm2)*EXP(-itm2)) &
     &               /spi-z*(ERFC_COMPBB(z)-ERFC_COMPBB(istm2)) &
     &               *(1.0D0-tm2/2.0D0))
            GFR(j) = GFR(j)/norm
            GFI(j) = 0.0
         ENDDO
         CALL DCFTN(GFR,GFI,N,1,icon)
         CALL DPNR(GFR,GFI,N,1,icon)
!      WRITE(6,*) 'FFT OF RESPONSE FUNCTION OF SCATTERING ENDED'
!      WRITE(6,*) '---GFR---'
!      WRITE(6,'(10D13.5)') (GFR(J),J=0,N1,10)
!      WRITE(6,*) '---GFI---'
!      WRITE(6,'(10D13.5)') (GFI(J),J=0,N1,10)
!     ... RESPONCE FOR THE FRONT SCATTERING
         tm = Te/mcc*ONCEF
         tm2 = tm*2.0D0
         stm2 = SQRT(tm2)
         istm2 = 1.0D0/stm2
         itm2 = 1.0D0/tm2
         norm = ERF_CMPBB(istm2) - 2.0D0*EXP(-itm2)/stm2/spi
         DO j = 0 , N1
            IF ( j.LT.N/2 ) THEN
               x = DX*j
            ELSE
               x = DX*(j-N)
            ENDIF
            e = EXP(x)
            ae = ABS(e-1)/(e+1)
            z = ae/stm2
!     ... BETA = 0->1
            gfr1(j) = (e+1.0D0)/(2.0D0*stm2) &
     &                *((EXP(-z*z)*(1.0D0-tm2)-(ae-tm2)*EXP(-itm2)) &
     &                /spi-z*(ERFC_COMPBB(z)-ERFC_COMPBB(istm2)) &
     &                *(1.0D0-tm2/2.0D0))
            gfr1(j) = gfr1(j)/norm
            gfi1(j) = 0.0
         ENDDO
         CALL DCFTN(gfr1,gfi1,N,1,icon)
         CALL DPNR(gfr1,gfi1,N,1,icon)
      ENDIF
 
      IF ( Ts1.NE.Ts .OR. Qnewe ) THEN
         Ts1 = Ts
         DO j = 0 , N2 - 1
            x = j*DX
            e = E0*EXP(x)
            w = e/Ts
            IF ( w.LT.170.0 ) THEN
               FLUx0(j) = e*e/(EXP(w)-1.0)
            ELSE
               FLUx0(j) = e*e*EXP(-w)
            ENDIF
!     Normalization photons/cm2/s/keV for R=1km at D=10kpc
            FLUx0(j) = FLUx0(j)*1.04E-3
            F0R(j) = FLUx0(j)
            F0I(j) = 0.0
         ENDDO
         DO j = N2 , N1
            FLUx0(j) = 0.0
            F0R(j) = FLUx0(j)
            F0I(j) = 0.0
         ENDDO
         CALL DCFTN(F0R,F0I,N,1,icon)
         IF ( icon.NE.0 ) THEN
            WRITE (6,*) 'DCFTN (F0R,F0I,)  RETURN CODE:' , icon
            STOP
         ENDIF
         CALL DPNR(F0R,F0I,N,1,icon)
!     WRITE(6,*) 'FFT OF B-B SPECTRUM ENDED'
!     WRITE(6,'(10D13.5)') (F0R(J),J=0,N1,50)
!     WRITE(6,'(10D13.5)') (F0I(J),J=0,N1,50)
 
      ENDIF
      DO j = 0 , N1
         FFR(j) = GFR(j)*F0R(j) - GFI(j)*F0I(j)
         FFI(j) = GFI(j)*F0R(j) + GFR(j)*F0I(j)
         FFR1(j) = gfr1(j)*F0R(j) - gfi1(j)*F0I(j)
         FFI1(j) = gfi1(j)*F0R(j) + gfr1(j)*F0I(j)
      ENDDO
      RETURN
!      DEBUG SUBCHK
      END
!*==setspb.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
!     ESCAPE PROBABILTY FOR N-TIMES SCATTERED PHOTON  (VER 2.2)
!     ASSUMING ISOTOROPIC SCATTERING (CROSS SECTION= SOLID ANGLE/4PI)
!     AND PLANE PARARELL GEOMETRY AND ISOTOROPIC PHOTON SOURCE
!     AT THE BOTTOM SURFACE OF THE HOT PLASMA CLOUD.
 
!     COPY RIGHT  K.MITSUDA  (ISAS)    1985  AUG.
 
!     SETSPB SETS THE ESCAPE PROBABILY IN /PROBAB/ SPB(0:50),
!     WHERE SPB(J) IS THE PROBABILY FOR J-TIMES SCATTERED PHOTON.
!     INPUT PARAMETER TAU IS THE OPTICAL DEPTH, SHOULD BE < 90.0 .
!     RATSPB= LIMIT(N->INFINITE) SPB(N+1)/SPB(N)
 
      SUBROUTINE SETSPB(Tau)
      IMPLICIT NONE
      DOUBLE PRECISION Tau
      COMMON /PROBAB/ SPB
      DOUBLE PRECISION SPB(0:50)
      COMMON /PROBAA/ RATspb
      DOUBLE PRECISION RATspb
      COMMON /MAXTIM/ FMAx
      INTEGER FMAx
 
      COMMON /OPTDEP/ TAU1
      DOUBLE PRECISION TAU1
      COMMON /CURDEP/ XX
      DOUBLE PRECISION XX
      EXTERNAL SPB1 , J2A , J2B , SPBN , JNA , JNB , E3
      DOUBLE PRECISION SPB1 , J2A , J2B , SPBN , JNA , JNB , E3
      DOUBLE PRECISION x1, x2 , xd
      DOUBLE PRECISION j(0:99) , w1 , w2
      INTEGER icon , nbin , n
      COMMON /SPLINE_COMPBB/ X , C , NBIn1 , MB
      DOUBLE PRECISION X(0:99) , C(100)
      INTEGER NBIn1 , MB , i
      DOUBLE PRECISION ratc0 , ratc1 , rata
      DATA nbin/20/
      DATA x1/0.0D0/ 
      DATA ratc0/ - 4.0711636/ , ratc1/ - 0.089848064/ , rata/0.6524/
 
!     I cannot find this subroutine (Ken Ebisawa 97/05/21)
!      CALL SETEPC
 
      IF ( Tau.LE.0.0 ) THEN
         SPB(0) = 1D0
         DO i = 1 , 50
            SPB(i) = 0D0
         ENDDO
      ELSE
         RATspb = (1.0D0+Tau**rata*EXP(ratc0+ratc1*Tau)) &
     &            *(1.0D0-(0.5D0-E3(Tau))/Tau)
 
         IF ( Tau.LE.0.07 ) THEN
            FMAx = 1
         ELSEIF ( Tau.LE.1.0 ) THEN
            FMAx = INT(9.3*Tau**0.499)
         ELSEIF ( Tau.LE.3.0 ) THEN
            FMAx = INT(9.3*Tau**0.855)
         ELSE
            FMAx = INT(7.586*Tau**1.041)
            IF ( FMAx.GT.50 ) FMAx = 50
         ENDIF
 
!     ... FMAX>=2 FOR THE CURRENT VERSION OF THE SUBROUTINE FOR FFT
         IF ( FMAx.LT.2 ) FMAx = 2
!     FMAX=50
         TAU1 = Tau
         MB = 5
         NBIn1 = nbin + 1
 
         SPB(0) = 2.0D0*E3(Tau)
 
         CALL DAQE(x1,Tau,SPB1,SPB(1),icon)
         IF ( icon.GT.10000 ) WRITE (6,*) ' DAQE FOR SPB1,  ICON = ' ,  &
     &                               icon
 
         X(0) = 0.0D0
         xd = Tau/nbin
         DO i = 1 , nbin - 1
            X(i) = X(i-1) + xd
         ENDDO
         X(nbin) = Tau
 
         DO i = 0 , nbin
            XX = X(i)
            CALL DAQE(x1,X(i),J2A,w1,icon)
            IF ( icon.GT.10000 ) WRITE (6,*) '  EAQE FOR J2A, ICON = ' ,  &
     &                                  icon
            x2 = Tau - X(i)
            CALL DAQE(x1,x2,J2B,w2,icon)
            IF ( icon.GT.10000 ) WRITE (6,*) '  EAQE FOR J2B, ICON = ' ,  &
     &                                  icon
            j(i) = (w1+w2)/4.0D0
         ENDDO
 
         CALL DBIC3(j,NBIn1,C,icon)
         IF ( icon.GT.0 ) WRITE (6,*) '  DBIC3 FOR J2, ICON = ' , icon
!     WRITE(6,*) 'CHKSPL: J2 '
!     CALL CHKSPL(J)
         CALL DAQE(x1,Tau,SPBN,SPB(2),icon)
         IF ( icon.GT.10000 ) WRITE (6,*) ' DAQE FOR SPB2,  ICON = ' ,  &
     &                               icon
         SPB(2) = SPB(2)*2.0D0
 
         DO n = 3 , FMAx
            DO i = 0 , nbin
               XX = X(i)
               CALL DAQE(x1,X(i),JNA,w1,icon)
               IF ( icon.GT.10000 ) WRITE (6,*) &
     &               '  EAQE FOR JNA, ICON = ' , icon
               x2 = Tau - X(i)
               CALL DAQE(x1,x2,JNB,w2,icon)
               IF ( icon.GT.10000 ) WRITE (6,*) &
     &               '  EAQE FOR JNB, ICON = ' , icon
               j(i) = (w1+w2)/2.0D0
            ENDDO
            CALL DBIC3(j,NBIn1,C,icon)
            IF ( icon.GT.0 ) WRITE (6,*) '  DBIC3 FOR JN, ICON = ' ,  &
     &                              icon
!     WRITE(6,*) 'CHKSPL: JN , N= ',N
!     CALL CHKSPL(J)
            CALL DAQE(x1,Tau,SPBN,SPB(n),icon)
            IF ( icon.GT.10000 ) WRITE (6,*) &
     &                                  ' DAQE FOR SPBN,  ICON = ' ,  &
     &                                  icon
            SPB(n) = SPB(n)*2.0D0
         ENDDO
      ENDIF
 
!$$$      WRITE(6,*) '!PROBABILITIES:'
!$$$      WRITE(6,*) '!tau =', tau
!$$$      do i = 0, 50
!$$$         WRITE(6,'(i2,D18.9)') i, max(REAL(SPB(i)),1E-10)
!$$$      end do
 
      RETURN
      END
!*==spb1.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      FUNCTION SPB1(X)
      IMPLICIT NONE
      DOUBLE PRECISION SPB1 , X
 
!     FUNCTION SUBROUTINE FOR THE ONCE SCATTERED PHOTON
 
      COMMON /OPTDEP/ TAU1
      DOUBLE PRECISION TAU1
      EXTERNAL E2_func
      DOUBLE PRECISION E2_func
      DOUBLE PRECISION x1
 
      x1 = TAU1 - X
      SPB1 = E2_func(x1)*E2_func(X)
 
      RETURN
      END
!*==spbn.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
      FUNCTION SPBN(Z)
      IMPLICIT NONE
      DOUBLE PRECISION SPBN , Z
 
!     FUNCTION SUBROUTINE FOR THE N-TIMES SCATTERED PHOTON
 
      COMMON /OPTDEP/ TAU1
      DOUBLE PRECISION TAU1
      DOUBLE PRECISION x1 , f
      INTEGER i , isw
 
      COMMON /SPLINE_COMPBB/ X , C , NBIn1 , MB
      DOUBLE PRECISION X(0:99) , C(100)
      INTEGER NBIn1 , MB
      EXTERNAL E2_func
      DOUBLE PRECISION E2_func
      INTEGER icon
      DATA i/1/ , isw/0/
 
      x1 = TAU1 - Z
      CALL DBIF3(X,NBIn1,C,x1,i,f,icon)
      IF ( icon.GT.10000 ) WRITE (6,*) '  DBIF3 IN SPBN,  ICON = ' ,  &
     &                                 icon
      SPBN = f*E2_func(Z)
 
      RETURN
      END
!*==j2a.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
      FUNCTION J2A(X)
      IMPLICIT NONE
      DOUBLE PRECISION J2A , X
 
!     FUNCTION SUBROUTINE FOR THE SOURCE FUNCTION OF
!     THE TWICE SCATTERED PHOTON, PART ONE
 
      COMMON /CURDEP/ XX
      DOUBLE PRECISION XX
      EXTERNAL E1 , E2_func
      DOUBLE PRECISION E1 , E2_func
      DOUBLE PRECISION x1
 
      x1 = XX - X
      J2A = E2_func(x1)*E1(X)
 
      RETURN
      END
!*==j2b.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
      FUNCTION J2B(X)
      IMPLICIT NONE
      DOUBLE PRECISION J2B , X
 
!     FUNCTION SUBROUTINE FOR THE SOURCE FUNCTION OF
!     THE TWICE SCATTERED PHOTON, PART TWO
 
      COMMON /CURDEP/ XX
      DOUBLE PRECISION XX
      EXTERNAL E1 , E2_func
      DOUBLE PRECISION E1 , E2_func
      DOUBLE PRECISION x1
 
      x1 = XX + X
      J2B = E2_func(x1)*E1(X)
 
      RETURN
      END
!*==jna.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
      FUNCTION JNA(Z)
      IMPLICIT NONE
      DOUBLE PRECISION JNA , Z
!     FUNCTION SUBROUTINE FOR THE SOURCE FUNCTION OF
!     THE N-TIMES SCATTERED PHOTON, PART ONE
 
      COMMON /SPLINE_COMPBB/ X , C , NBIn1 , MB
      DOUBLE PRECISION X(0:99) , C(100)
      INTEGER NBIn1 , MB
      COMMON /CURDEP/ XX
      DOUBLE PRECISION XX
      EXTERNAL E1
      DOUBLE PRECISION E1
 
      DOUBLE PRECISION x1 , f
      INTEGER i , isw
      INTEGER icon
      DATA i/1/ , isw/0/
 
      x1 = XX - Z
      CALL DBIF3(X,NBIn1,C,x1,i,f,icon)
      IF ( icon.GT.10000 ) WRITE (6,*) '  DBIF3 IN JNA,  ICON = ' , icon
      JNA = f*E1(Z)
 
      RETURN
      END
!*==jnb.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
      FUNCTION JNB(Z)
      IMPLICIT NONE
      DOUBLE PRECISION JNB , Z
!     FUNCTION SUBROUTINE FOR THE SOURCE FUNCTION OF
!     THE N-TIMES SCATTERED PHOTON,  PART TWO
 
      COMMON /SPLINE_COMPBB/ X , C , NBIn1 , MB
      DOUBLE PRECISION X(0:99) , C(100)
      INTEGER NBIn1 , MB
      COMMON /CURDEP/ XX
      DOUBLE PRECISION XX
      EXTERNAL E1
      DOUBLE PRECISION E1
 
      DOUBLE PRECISION x1 , f
      INTEGER i , isw
      INTEGER icon
      DATA i/1/ , isw/0/
 
      x1 = XX + Z
      CALL DBIF3(X,NBIn1,C,x1,i,f,icon)
      IF ( icon.GT.10000 ) WRITE (6,*) '  DBIF3 IN JNB,  ICON = ' , icon
      JNB = f*E1(Z)
 
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION E1(X)
!     E_n(x) = \int_1^\infinity{exp(-xs)s^{-n}}\;ds
      IMPLICIT NONE
      DOUBLE PRECISION X
      INTEGER i , n
!      double precision t
      DOUBLE PRECISION s_low , s_high , delta , step
 
      E1 = 0.0D0
      IF ( X.GE.100.0 ) THEN
         RETURN
      ELSEIF ( X.LT.0.3 ) THEN
         s_low = LOG(X)
         n = 10
         DO i = 1 , n
            s_high = s_low + 1.0D0
            E1 = E1 + (s_high-s_low) &
     &           *0.5D0*(EXP(-(EXP(s_low)))+EXP(-(EXP(s_high))))
            s_low = s_high
         ENDDO
      ELSE
!$$$         n = 100
!$$$         step=1.1D0
!$$$         delta = 0.001D0
!$$$         do i = 1, n
!$$$            s_low  = 1.0D0+step**(i-1)*delta
!$$$            s_high = 1.0D0+step**(i  )*delta
!$$$            E1=E1+(s_high-s_low)*(exp(-s_low*x)
!$$$     $           / s_low+exp(-s_high*x)/s_high)/2.0D0
!$$$         end do
         n = 20
         step = 1.2D0
         s_high = 1.0D0
         DO i = 1 , n
            s_low = s_high
            delta = 0.07D0*step**(i-1)
            s_high = s_low + delta
            E1 = E1 + (s_high-s_low) &
     &           *(EXP(-s_low*X)/s_low+EXP(-s_high*X)/s_high)/2.0D0
!            write(*,*) i, delta, s_low, s_high, E1
         ENDDO
      ENDIF
!      write(*,*) 'debug E1', X, E1
      END
!*==e2.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      DOUBLE PRECISION FUNCTION E2_func(X)
      IMPLICIT NONE
      DOUBLE PRECISION X , E1
      EXTERNAL E1
 
      E2_func = EXP(-X) - X*E1(X)
 
      END
!*==e3.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      DOUBLE PRECISION FUNCTION E3(X)
      IMPLICIT NONE
      DOUBLE PRECISION X , E2_func
      EXTERNAL E2_func
 
      E3 = 1.0D0/2.0D0*(EXP(-X)-X*E2_func(X))
 
      END
!*==daqe.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      SUBROUTINE DAQE(Low_b,High_b,FUNCTION,Output,Icon)
      IMPLICIT NONE
      DOUBLE PRECISION Low_b , High_b , FUNCTION , Output
      INTEGER Icon
      EXTERNAL FUNCTION
 
!     Numerical integral.  Original routine is called from FUJUTSU SSL
!     library (source code is not free).
!
!     low_b  (in): lower bounary
!     high_b (in): upper bounary
!     function (in): function to integrate (external)
!     output (out): output integral value
!     icon (out):
 
      DOUBLE PRECISION dx , xm , xr , w(5) , x(5) , low_b_d , high_b_d ,  &
     &                 h , output_d
      INTEGER j , jj , n
      SAVE w , x
      DATA w/.2955242247D0 , .2692667193D0 , .2190863625D0 ,  &
     &     .1494513491D0 , .0666713443D0/
      DATA x/.1488743389D0 , .4333953941D0 , .6794095682D0 ,  &
     &     .8650633666D0 , .9739065285D0/
 
      n = 10
      IF ( High_b.NE.Low_b ) THEN
         h = (High_b-Low_b)/DBLE(n)
         Output = 0.0D0
         DO jj = 1 , n
            low_b_d = Low_b + h*(jj-1)
            high_b_d = Low_b + h*jj
            xm = 0.5D0*(high_b_d+low_b_d)
            xr = 0.5D0*(high_b_d-low_b_d)
            output_d = 0.0D0
            DO j = 1 , 5
               dx = xr*x(j)
               output_d = output_d + w(j) &
     &                    *(FUNCTION(xm+dx)+FUNCTION(xm-dx))
            ENDDO
            Output = Output + output_d
         ENDDO
         Output = Output*xr
      ELSE
         Output = 0.0D0
      ENDIF
      Icon = 0
      END
!*==erf_cmpbb.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      DOUBLE PRECISION FUNCTION ERF_CMPBB(X)
      IMPLICIT NONE
      EXTERNAL ERFC_COMPBB
      DOUBLE PRECISION X , ERFC_COMPBB
      ERF_CMPBB = 1.0D0 - ERFC_COMPBB(X)
      END
!*==erfc_compbb.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      DOUBLE PRECISION FUNCTION ERFC_COMPBB(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      DOUBLE PRECISION a , b , PI
      PARAMETER (PI=3.14159265358979D0)
 
      INTEGER icon
      EXTERNAL FUNC_ERFC
      DOUBLE PRECISION FUNC_ERFC
 
      a = ATAN(X)
      b = PI/2.0D0
 
      CALL DAQE(a,b,FUNC_ERFC,ERFC_COMPBB,icon)
      ERFC_COMPBB = ERFC_COMPBB*2.0D0/SQRT(PI)
 
      END
!*==func_erfc.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      DOUBLE PRECISION FUNCTION FUNC_ERFC(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      FUNC_ERFC = EXP(-(TAN(X))**2)/(COS(X))**2
 
      END
!*==dcftn.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      SUBROUTINE DCFTN(Data_r,Data_i,N,Isign,Icon)
      
      ! to call the fourier transform subroutine
      USE FFTW3_wrap
      
      IMPLICIT NONE
      INTEGER N , Isign , Icon
      DOUBLE PRECISION Data_r(N) , Data_i(N)
 
      INTEGER NMAX , i
      PARAMETER (NMAX=32768)
      DOUBLE PRECISION data(NMAX)
 
      DO i = 1 , N
         data(2*i-1) = Data_r(i)
         data(2*i) = Data_i(i)
      ENDDO
 
      CALL fftw_wrap_four1_dbl(data,N,Isign)
 
      DO i = 1 , N
         Data_r(i) = data(2*i-1)
         Data_i(i) = data(2*i)
      ENDDO
      Icon = 0
      END
!*==dpnr.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      SUBROUTINE DPNR(Data_r,Data_i,N,Isign,Icon)
      IMPLICIT NONE
      INTEGER N , Isign , Icon
      DOUBLE PRECISION Data_r(*) , Data_i(*)
 
!     dummy routine - most of this to suppress compiler warnings
      Icon = INT(Data_r(1))
      Icon = INT(Data_i(1))
      Icon = N
      Icon = Isign
      Icon = 0
      END
!*==dbic3.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      SUBROUTINE DBIC3(Indata,N,C,Icon)
      IMPLICIT NONE
      INTEGER N , Icon
      DOUBLE PRECISION Indata(N)
      DOUBLE PRECISION C(N)
!     dummy routine for interpolation no.1
!     C and XT are created for interpolation parameters, which are
!     used in DBIF3 to obtain interpolated values.
!     actually, just indata are copied to C.
      INTEGER i
 
      DO i = 1 , N
         C(i) = Indata(i)
      ENDDO
      Icon = 0
      END
!*==dbif3.spg  processed by SPAG 4.50J  at 18:42 on 21 Jul 1999
 
      SUBROUTINE DBIF3(X,N,C,V,I_out,Fl,Icon)
      IMPLICIT NONE
      INTEGER N , I_out , Icon
      DOUBLE PRECISION X(N) , C(N), V , Fl
!     dummy routine for interpolation no.2
!     C created in DBIC3 is used in DBIF3 to obtain interpolated values.
      INTEGER i

      DO i = 1 , N - 2
         IF ( V.GE.X(i) .AND. V.LT.X(i+1) ) THEN
            I_out = i
            GOTO 100
         ENDIF
      ENDDO
      i = N - 1
      IF ( V.GE.X(i) .AND. V.LE.X(i+1) ) THEN
         I_out = i
         GOTO 100
      ENDIF
 
 
 
!     This is an error
      WRITE (*,*) 'DBIF3 error. v=' , V
      WRITE (*,*) 'x=' , (X(i),i=1,N)
      WRITE (*,*) 'c=' , (C(i),i=1,N)
      Icon = 9999
      RETURN
 
 100  Fl = (C(I_out+1)*(V-X(I_out))-C(I_out)*(V-X(I_out+1))) &
     &     /(X(I_out+1)-X(I_out))
 
      Icon = 0
      END
 
 
 
 
