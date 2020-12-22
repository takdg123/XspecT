      SUBROUTINE GRAD(Ear,Ne,Param,Ifl,Photar,Photer)
      IMPLICIT NONE
      INTEGER Ne , Ifl
      REAL Ear(0:Ne) , Param(6) , Photar(Ne), Photer(Ne)
 
C     This program was originally written by T. Hanawa in
C     1989 for ISAS GINGA SPFD spectral fitting program.  
C     Ported to xspec by Ken Ebisawa.
C     See Hanawa 1989, ApJ, 341, 948, and Ebisawa, Mitsuda
C     and Hanawa 1991, ApJ, 367, 213
 
C     Several bugs were pointed out by Dr. Dipankar Bhattacharya in 
C     late 2000, and actually three bugs were confirmed and fixed by 
C     T. Hanawa and Ken Ebisawa (search 'bug' in this file) .
C     In short, effect of these bugs was such that the
C     old GRAD spectra with the mass M are identical to
C     the new GRAD spectra (current code) with M/1.4.
C     In other words, the mass obtained by fitting the old GRAD model
C     to the observation was 1.4 times over-estimated.
C     See the report at the following location for detail.
C     http://lheawww.gsfc.nasa.gov/users/ebisawa/report.ps.gz

C     A new flag (sixth parameter) was introduced to switch on/off 
C     the relativistic effects.
C     The disk flux is calculated to more outer radius, 
C     so that the reliable lower energy of the model is ~0.01 keV 
C     (it used to be 0.02 keV).
C     Old legacy of the ISAS SPFD coding style was removed.
C     Last modified by Ken Ebisawa 2001-02-01 (ebisawa@subaru.gsfc.nasa.gov).
 
C     Param(1): Distance in kpc
C     Param(2): i (degree, face-on=0.0)
c     Param(3): Mass (in solar mass unit)
c     Param(4): Mdot (in 10^18 g/s)
C     Param(5): Tcol/Teff (nerver allowed to be free)
C     Param(6): A new flag to switch on/off the relativistic effects.
C     If positive, relativistic calculation.  If negative or zero,
C     Newtonian calculation (nerver allowed to be free).
 
      DOUBLE PRECISION oldinc, flux_high , flux_low
      INTEGER ibin,  n1max, n2max
      real oldflag

      SAVE oldflag, oldinc, n1max, n2max

      DATA oldflag, oldinc /0, 9999D0/
      DATA n1max, n2max / 100, 12/

c suppress a warning message from the compiler
      ibin = ifl

c this model has no errors
      DO ibin = 1, ne
         photer(ibin) = 0.0
      ENDDO

C     n1max is the disk radial bin (<=800), n2max is the
C     phase bin (<=72)
C     It seems like n1max=100 and n2max=12 is precise and
C     fast enough to calculate.

      IF ( DABS(oldinc-Param(2)).GT.1D-5.or.
     $     oldflag*param(6).le.0.0) THEN
         CALL GRAD_DISK(DBLE(SIN(Param(2)*1.745329D-2)), n1max,n2max,
     $        INT(param(6)))
         oldinc  = Param(2)
         oldflag = param(6)
      ENDIF

      DO ibin = 1 , Ne
         CALL FLUX(DBLE(Ear(ibin)),flux_high,DBLE(Param(1)),
     $        DBLE(Param(3)),DBLE(Param(4)),
     &        DBLE(SIN(Param(2)*1.745329D-2)),DBLE(Param(5)),
     $        n2max)
         CALL FLUX(DBLE(Ear(ibin-1)),flux_low,DBLE(Param(1)),
     $        DBLE(Param(3)),DBLE(Param(4)),
     &        DBLE(SIN(Param(2)*1.745329D-2)),DBLE(Param(5)),
     $        n2max)
         photar(ibin) = SNGL((flux_high+flux_low)
     &                       *(Ear(ibin)-Ear(ibin-1))/2.0)
      ENDDO
      END
 
      SUBROUTINE FLUX(E,F,D,Em,Dmdt,Sini,Ratio,N2max)
      IMPLICIT NONE
      double precision E, F, D, Em, Dmdt, Sini, Ratio
      integer  n2max

      DOUBLE PRECISION t

      INTEGER NINT
      PARAMETER (NINT=500)
      double precision GA(0:NINT) , DBS(0:NINT)
      COMMON /DISKR /GA , DBS

      INTEGER k

      DOUBLE PRECISION BCALC
      EXTERNAL BCALC

      F = 0.0D0
      DO k = 1 , NINT
         IF ( DBS(k).GT.0.0 ) THEN
c            t = 1.11267D0*(Ratio/1.5D0)*Dmdt**0.25D0/SQRT(Em/1.4D0)
C     The above is an original bug (2001-2-1).
            t = 1.11267D0*(Ratio/1.5D0)*Dmdt**0.25D0/SQRT(Em)
     &           *GA(k)/0.15D0
            F = F + BCALC(E,t)*DBS(k)                                         
         ENDIF
      end do
c      F = 8.9038633D-3*F/Ratio**4/(D/10.0)**2*(Em/1.4)**2
C     The above is an original bug (2001-2-1).
      F = 8.9038633D-3*F/Ratio**4/(D/10.0)**2*Em**2
      F = F*SQRT(1.0D0-Sini**2)/DBLE(N2max)
      RETURN
      END

C     PROGRAM FOR CALCULATING RADIATION FROM ACCRETION DISK
C
C     (GENERAL RELATIVISTIC COMPUTAITION)   1988.05.23
C
      SUBROUTINE GRAD_DISK(Sini,N1max,N2max, relflag)

      IMPLICIT NONE

      double precision sini
      integer n1max, n2max
C     relflag is newly added flag to switch on(relflag>0)/off(relflag<=0) 
C     the relativistic calculation.  2001-1-17 Ken Ebisawa
      integer relflag

      DOUBLE PRECISION b, c, dbh, f ,phid, pi, r, red
      INTEGER i, j, k, dummy

      double precision RR(0:800,72),G0(0:800,72),BB(0:800),DB(0:800)           

      INTEGER NINT
      PARAMETER (NINT=500)
      double precision GA(0:NINT) , DBS(0:NINT)
      COMMON /DISKR /GA, DBS

      dbh = 8.0D0/DBLE(N1max)
      pi = 4.0D0*ATAN(1.0D0)
      c = (1.0D0-SQRT(0.5D0))/(1.0D0+SQRT(0.5D0))
      DO i = 0 , N1max
         bb(i) = 9.0D0*10**(dbh*DBLE(i))
         DO j = 1 , N2max
            phid = 2.0*pi*(DBLE(j-1))/DBLE(N2max)
            b = SQRT(bb(i)*(1.0D0-(Sini*SIN(phid))**2))
            CALL RAY(r,b,Sini,phid,red, relflag)
            rr(i,j) = r
            IF ( r.LE.3.0D0 ) THEN
               g0(i,j) = 0.0
            ELSE
               if(relflag.le.0) then
                  G0(I,J)=R**(-0.75D0)*(1.0D0-SQRT(3.0D0/R))**0.25D0
               else
                  IF ( r.LE.3.0003D0 ) THEN
                     f = 2.0*(1.0D0-SQRT(3.0D0/r))
     &                    **2*(1.0D0-4.0D0/3.0D0*(1.0D0-SQRT(3.0D0/r)))
                  ELSE
                     f = 1.0D0 - SQRT(3.0D0/r) + SQRT(0.375/r)
     &            *LOG((1.0D0+SQRT(1.5D0/r))/(1.0D0-SQRT(1.5D0/r))*c)
                  ENDIF
                  f = r**(-0.75D0)*SQRT(SQRT(f))
                  g0(i,j) = f*(1.-1.5D0/r)**0.25D0*red
               endif
            ENDIF
         ENDDO
      ENDDO
      db(0) = (bb(1)-bb(0))/2.0
      DO i = 1 , N1max - 1
         db(i) = (bb(i+1)-bb(i-1))/2.0
      ENDDO
      db(N1max) = (bb(N1max)-bb(N1max-1))/2.0
      DO k = 0,  nint
         DBS(k) = 0.0
         GA(k) = 0.0
      ENDDO
      DO j = 1 , N2max
         DO i = 0 , N1max
            dummy = INT(g0(i,j)*DBLE(4*NINT))
            k = dummy
            GA(k) = GA(k) + g0(i,j)*db(i)
            DBS(k) = DBS(k) + db(i)
         ENDDO
      ENDDO
      DO k = 0 , NINT
         IF ( DBS(k).GT.0.0 ) THEN
            GA(k) = GA(k)/DBS(k)
         ELSE
            GA(k) = (DBLE(k)+5.0D-1)/DBLE(4*NINT)
         ENDIF
      ENDDO

c$$$            WRITE(6,600) SINI
c$$$ 600        FORMAT(1H1,//,30X,'SUMMARY OF PHOTON TRAJECTORIES',10X,
c$$$     1           'SIN I = ',F8.4,//)
c$$$            WRITE(6,601)
c$$$ 601        FORMAT(1H0,6X,'B',4X,'(R,G)',3X,'PHI = 0',14X,'90',17X,
c$$$     1           '180',16X,'270',/)
c$$$            WRITE(6,602)(SQRT(BB(I)),(RR(I,J),G0(I,J),J=1,3*N2MAX/4+1,
c$$$     1           N2MAX/4),I=0,N1MAX,10)
c$$$ 602        FORMAT(1H ,0P,2F11.6,F8.5,F11.6,F8.5,F11.6,F8.5,F11.6,F8.5)
c$$$            WRITE(6,601)
      RETURN
      END
C
C   COMPUTING A TRAJECTORY OF A PHOTON IN THE SCHWARZSCHILD METRIC
C
      SUBROUTINE RAY(R0,B,Sini,Phid,Red, relflag)
      IMPLICIT NONE
      DOUBLE PRECISION ACOS , B , bor , COS , cosi ,  dbordx , 
     &                 dphdx , drdx , dtdx , dtdx2 , dthdx , dx1 , dx2 , 
     &                 dx3 , dx4 , dy1 , dy2 , dy3 , dy4
      DOUBLE PRECISION h , phi , Phid , R0 , Red , SIN , Sini , 
     &                 sinth , SQRT , th
C     relflag is newly added flag to switch on(relflag>0)/off(relflag<=0) 
C     the relativistic calculation.  2001-1-17 Ken Ebisawa
      integer relflag

      INTEGER n , i

      DOUBLE PRECISION FCALC
      EXTERNAL FCALC
 
C   B  : IMPACT PARAMETER
C   R0 : THE RADIUS AT THE DISK
C  BOR : = B / R
 
      COMPLEX z , arg

      if(relflag.le.0) then
         R0 = b/sqrt(1.0-sini**2*sin(phid)**2)
         Red = 1.0
         return
      endif

      bor = 0.0D0
      dbordx = 1.0D0
      n = 80
      cosi = SQRT(1.0D0-Sini*Sini)
      z = cmplx(real(COS(PHid)), real(SIN(Phid)*cosi))
      arg = LOG(z)
      phi = AIMAG(arg)
      th = ACOS(Sini*SIN(Phid))
      sinth = SIN(th)
C  INTEGRATION
      h = th/DBLE(n)
      R0 = 1.0D+2
      DO i = 1 , n
         IF ( R0.GT.2.0 ) THEN
            dx1 = h*dbordx
            dy1 = h*FCALC(bor,B)
            dx2 = h*(dbordx+0.5*dy1)
            dy2 = h*FCALC(bor+0.5*dx1,B)
            dx3 = h*(dbordx+0.5*dy2)
            dy3 = h*FCALC(bor+0.5*dx2,B)
            dx4 = h*(dbordx+dy3)
            dy4 = h*FCALC(bor+dx3,B)
            bor = bor + (dx1+2.0D0*dx2+2.0D0*dx3+dx4)/6.0D0
            dbordx = dbordx + (dy1+2.0D0*dy2+2.0D0*dy3+dy4)/6.0D0
            R0 = B/bor
         ENDIF
      ENDDO
      drdx = -R0*R0*dbordx/B
      dthdx = cosi/sinth
      IF ( SIN(phi).GE.0.5D0 ) THEN
         dphdx = -COS(Phid)*SIN(Phid)*Sini*cosi/(SIN(phi)*sinth**2)
      ELSE
         dphdx = -(COS(Phid))**2*Sini/COS(phi)/sinth**2
      ENDIF
      dtdx2 = (drdx**2/(1.0D0-1.0D0/R0)+(R0*dthdx)**2+(R0*dphdx)**2)
     &        /(1.0D0-1.0D0/R0)
      dtdx = -SQRT(dtdx2)
c      Red = 1.0D0/(1.0D0+SQRT(R0/2.0D0)/(1.0D0-1.0D0/R0)*dphdx/dtdx)
C     The above is an original bug (2001-2-1).
      Red =
     $1.0D0/(1.0D0+SQRT(R0/2.0D0)/(1.0D0-1.0D0/R0)**0.5D0*dphdx/dtdx)
      RETURN
      END

      DOUBLE PRECISION FUNCTION BCALC(ep, tb)
      DOUBLE PRECISION ep, tb

      BCALC = ep**2D0*EXP(-ep/tb)/(1.0D0-EXP(-ep/tb))
      RETURN
      END

      DOUBLE PRECISION FUNCTION FCALC(x, B)
      DOUBLE PRECISION x, B

      FCALC = -x*(1.0D0-1.5D0*x/B)
      RETURN
      END
