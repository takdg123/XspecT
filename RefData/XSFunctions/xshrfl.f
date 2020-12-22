      subroutine xshrfl(ear,ne,param,ifl, photar, photer)

c simple reflection model (see below) - valid for E<15 keV in the rest frame
c T. Yaqoob sometime between 1991-1993; cleaned up a lot by C. Day.

      integer ne, ifl
      real ear(0:ne), param(8), photar(ne), photer(ne)

      REAL a,b,hreflutl,zfac,e1,e2
      INTEGER i 

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO


      zfac=1.0+param(8)
      DO i = 1, ne
         e1=zfac*ear(i-1)
         e2=zfac*ear(i)
         a=hreflutl(param(1),param(2),param(3),e1,param(4),
     1        param(5))
         b=hreflutl(param(1),param(2),param(3),e2,param(4),
     1        param(5))
         photar(i) = 0.5*(2.0*param(6)+param(7)*(a+b))
      enddo
      return
      end
      REAL FUNCTION HREFLUTL(D1,D2,TOBS,E,AFE,EKEDGE)
C
C     Function computes the reflected spectrum from a cold semi-infinite
C     slab in the approximation of elastic electron scattering, i.e.
C     incident photon energy << mc^2 in the electron rest-frame. Should
C     be good for photon energies up to ~15 keV.
C     
C     USEAGE: Suppose direct spectrum is: spec (photons/cm/cm/s/keV)
C     then the TOTAL observed spectrum is spec*(1+hrefl) in the same
C     units.
C     
C     PARAMETERS:
C        D1 ....... Minimum angle (degrees) between source photons
C                      incident on the slab and the slab normal.
C        D2 ....... Maximum angle...
C        TOBS ..... Angle (degrees) between the observer's line of sight
C                      and the slab normal.
C        E ........ Energy of incident photon.
C        AFE ...... Iron abundance relative to Solar.
C        EKEDGE ... Iron K-edge energy.
C     Recommended that D1,D2,FEBUN and EKEDGE be fixed.
C
C     The routine uses the solutions to the transfer equations utilizing
C     the H-functions [see Basko, 1978]. Analytic approxiamtions to the
C     H-functions are used, as well as an average of the H-functions for
C     the incident rays. Errors from these approximations are typically
C     less than a few percent, snf the results agree well with Monte-
C     Carlo simulations.
C
      REAL D1,D2,TOBS,E,AFE,EKEDGE,C1,C2,XMU,ALBEDO,ALBD,HMEAN,HFUNA,
     1   F2
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592654D0)

      C1=SNGL(DCOS(D2*PI/180.0d0))
      C2=SNGL(DCOS(D1*PI/180.0d0))
      XMU=SNGL(DCOS(TOBS*PI/180.0d0))
      IF(XMU.LT.-1.0) XMU=-1.0
      IF(XMU.GT.1.0) XMU=1.0
      ALBEDO=ALBD(AFE,EKEDGE,E)
      IF(XMU.EQ.0.0) THEN
         F2=XMU
      ELSE
         F2=XMU*LOG((C2+XMU)/(C1+XMU))
      ENDIF
      HREFLUTL=ALBEDO*HMEAN(ALBEDO)*HFUNA(ALBEDO,XMU)*F2*0.5
      END
C
C -------
C
      REAL FUNCTION ALBD(AFE,EKEDGE,E)
C     Computes the ratio of the Thompson Cross-section to the total
C     (scattering + absorption) cross-section.
      REAL AFE,EKEDGE,E,THOM,ABSCRS
      DATA THOM/0.665E-3/
      ALBD=THOM/(THOM+(ABSCRS(E,1.0,AFE,EKEDGE)/1.2))
      END
C
C -------
C
      REAL FUNCTION HFUNA(ALBEDO,XMU)
      REAL ALBEDO,XMU
      HFUNA=(1.0+SQRT(3.0)*XMU)/(1.0+SQRT(3.0*(1.0-1.0*ALBEDO))*XMU)
      END
C
C -------
C
      REAL FUNCTION HMEAN(ALBEDO)
      REAL ALBEDO,A,B,R
      A=SQRT(3.0)
      IF ( ALBEDO .EQ. 1. ) THEN
         HMEAN = 1.0 + A/2.
      ELSE
         B=SQRT(3.0*(1.0-1.0*ALBEDO))
         R=A/B
         HMEAN=1.0+(R-1.0)*(1.0-1.0*(LOG(1.0+B)/B))
      ENDIF
      END
C
C -------
C
      REAL FUNCTION ABSCRS(E,CD,AFE,EKEDGE)
C
C     Computes photoelectric absorption.
C
C     Uses piecewise polynomial fit of Morrison and McCammon Ap.J. 270,
C     119 for range 0.03 to 10 keV. Below 0.03 keV uses power law fit to
C     hydrogen and helium edge profiles interpolated/extrapolated from
C     Henke data (1982) also used by Morrisom and McMammon. Above 10 keV
C     crude "eyeball" fit provided by Gordon Stewart!? Corrected error
C     in low energy stuff (CD missing!) RW 1988-May-27. Author Dick
C     Willingale 1986-Sep-4
C
C     PARAMETERS
C        E ........ Energy (keV)
C        CD ....... Column density 10^21 cm-2
C        EKEDGE ... Energy (keV)
C
      REAL AFE,E,CD,EKEDGE
      REAL C0(14),C1(14),C2(14),ET(15),TAU,FEABS,AEFF
      INTEGER K
      DATA C0/17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3, 202.7,342.7,
     1   352.2,433.9,629.0,701.2/
      DATA C1/608.1,267.9,18.8,66.8,145.8,-380.6,169.3,146.8, 104.7,
     1   18.7,18.7,-2.4,30.9,25.2/
      DATA C2/-2150.,-476.1,4.3,-51.4,-61.1,294.,-47.7,-31.5,-17.,0.,0.,
     1   0.75,0.,0./
      DATA ET/.03,.1,.284,.4,.532,.707,.867,1.303,1.84,2.471,3.21,4.038,
     1   7.111,8.331,10./
C
      TAU=0.0
      IF(E.LT.0.0136) THEN
         TAU=0.0
      ELSE
         IF((E.GE.0.0136).AND.(E.LT.0.0246)) THEN
            TAU=CD*9.49E-3/(E**3.26)
         ELSE
            IF((E.GE.0.0246).AND.(E.LT.0.03)) THEN
               TAU=CD*8.63E-3/(E**3.26)+289.E-3/(E**2.11)
            ELSE
               IF(E.GE.10.) THEN
                  TAU=(0.30*CD)/(E*E*SQRT(E))
               ELSE
                  DO 10 K=1,14
                     IF((E.GE.ET(K)).AND.(E.LT.ET(K+1))) THEN
                        TAU=(C0(K)+C1(K)*E+C2(K)*E*E)/E/E/E
                        TAU=0.001*TAU*CD
                     ENDIF
10                CONTINUE
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      AEFF=1.0-AFE
      ABSCRS=TAU-FEABS(E,CD,AEFF,EKEDGE)
      END
C
C -------
C
      REAL FUNCTION FEABS(E,CD,AFE,EKEDGE)
C
C     Computes iron absorption
C     PARAMETERS
C        EKEDGE ... Iron K-edge energy (keV)
C        E ........ Energy (keV)
C        AFE ...... Iron abundance relative to Solar
C        CD ....... Column density in units of 10^16.52 cm-2 so that
C                      values are equivalent to hydrogen column in units
C                      10^21 cm-2, assuming abundance ratio Fe/H of
C                      10^(7.52-12).
C
      REAL E,CD,AFE,EKEDGE
      IF((E.GE.0.73).AND.(E.LT.EKEDGE)) THEN
         FEABS=0.043*AFE*CD*(0.73/E)**2.2
      ELSE
         IF(E.GE.EKEDGE) THEN
            FEABS=0.0012*AFE*CD*(EKEDGE/E)**3.11
         ELSE
            FEABS=0.0
         ENDIF
      ENDIF
      END

