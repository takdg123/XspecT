C- This code was written by Jonathan Koehane, and renamed at SPR's request
C- srcut -- Supernova Remnant Cut-off model
C- -- KD 
C- 11/28/1998
C- 
C -       SUBROUTINE syncsr(ear, ne, param, ifl, photar)
       SUBROUTINE srcut(ear, ne, param, ifl, photar, photer)
C---
C Steve Reynolds simple models -- put in X-spec format
C---
C see ADDMOD for parameter descriptions
C number of parameters: 2
C      norm      1 GHz Flux in Jy 
C       1        Radio Energy spectral index (around 0.5)  
C       2        Break Frequency in Hz
C no intrinsic energy range limit
C Coeficient is the 1 GHz flux in Jy
C---
C 08 Apr 1997 - Jonathan Keohane
C---

       INTEGER ne, ifl
       REAL ear(0:ne), param(2), photar(ne), photer(ne)
       INTEGER NMOD
       REAL kevhz, kevfljy, logmin,  logdelta
       PARAMETER (kevhz = 2.417988366E17, kevfljy = 1509.188961)
C  1kev = kevhz Hz,  1 keV/(s kev cm2) = kevfljy Jy  
C              number of bins in model
       Parameter(NMOD = 80)
       Parameter(logmin = -4.1, logdelta = 0.1) 
       REAL a, b, Ebr, Ebr2, recGhz, alphan, gamma
       REAL logR(NMOD)
       INTEGER i, j
       REAL  terp, t, t0

C     10**logR is the factor by which the flux is reduced
      DATA logR/
     &  0.0000E+00, -.1106E-02, -.2339E-02, -.3715E-02, -.5248E-02, 
     &  -.6956E-02, -.8859E-02, -.1098E-01, -.1333E-01, -.1595E-01, 
     &  -.1886E-01, -.2209E-01, -.2567E-01, -.2965E-01, -.3406E-01, 
     &  -.3894E-01, -.4435E-01, -.5033E-01, -.5695E-01, -.6426E-01, 
     &  -.7234E-01, -.8126E-01, -.9109E-01, -.1019E+00, -.1139E+00, 
     &  -.1270E+00, -.1415E+00, -.1574E+00, -.1749E+00, -.1941E+00, 
     &  -.2152E+00, -.2383E+00, -.2636E+00, -.2914E+00, -.3217E+00, 
     &  -.3550E+00, -.3913E+00, -.4309E+00, -.4742E+00, -.5215E+00, 
     &  -.5730E+00, -.6292E+00, -.6904E+00, -.7570E+00, -.8295E+00, 
     &  -.9084E+00, -.9942E+00, -.1087E+01, -.1189E+01, -.1299E+01, 
     &  -.1418E+01, -.1548E+01, -.1689E+01, -.1841E+01, -.2007E+01, 
     &  -.2186E+01, -.2381E+01, -.2591E+01, -.2820E+01, -.3067E+01, 
     &  -.3335E+01, -.3625E+01, -.3940E+01, -.4280E+01, -.4649E+01, 
     &  -.5048E+01, -.5480E+01, -.5947E+01, -.6453E+01, -.7000E+01, 
     &  -.7591E+01, -.8230E+01, -.8922E+01, -.9670E+01, -.1048E+02, 
     &  -.1135E+02, -.1230E+02, -.1332E+02, -.1442E+02, -.1561E+02 /

C---

      recGhz  =  1E-9 * kevhz
      alphan = -1.0*param(1)
      Ebr = param(2)/kevhz

      Ebr2 = log10(2*Ebr)
   
C -- log10(F) = logmin + i * logdelta  where i = 1, 2, 3 etc ...

      gamma = param(1) + 1 
      a = kevfljy * recGhz**alphan
      b = 0.0

      CALL xspwlw(ear, ne, gamma, ifl, photar, photer)


C     PRINT *, 'a = ', a
C     PRINT *, 'Ebr = ', Ebr
C     PRINT *, 'alphan = ', alphan

      t0 = (Ebr2 + logmin) / logdelta
        
        DO j = 1, ne
          t = (log10(Ear(j-1)+Ear(j)))/logdelta - t0
C     PRINT *,'t = ', t
          IF (t .gt. float(NMOD)) then 
            b = 0.0
          ELSEIF (t .lt. 1.0) then 
            b = 1.00
          ELSE
            i = int(t)
            terp = t - float(i)  
            b  = 10.0**( (1.0-terp)*logR(i) + terp*logR(i+1) )
          ENDIF

          photar(j) = a*b*photar(j)
         
        ENDDO
           
        RETURN
        END
