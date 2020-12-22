C- Written by Jonathan Koehane, modified by Kristy Dyer
C- sresc -- Supernova Remnant EScape Model
C- -- KD 
C- 3/8/1999
C- 
C -       SUBROUTINE syncsr(ear, ne, param, ifl, photar)
       SUBROUTINE sresc(ear, ne, param, ifl, photar, photer)
C---
C Steve Reynolds simple models -- put in X-spec format
C---
C see ADDMOD for parameter descriptions
C number of parameters: 2
C      norm      1 GHz Flux in Jy 
C       1        Radio Energy spectral index (around 0.5)  
C       2        Break Frequency in Hz
C (The break frequency should be recast in units of 10e16 Hz or KeV 
C XSPEC is happiest if default values are close to 1)
C no intrinsic energy range limit
C Coeficient is the 1 GHz flux in Jy
C---
C 3/8/1999 K.K. Dyer
C---

! The following five variables are used to call powerLaw
       INTEGER ne ! number of energy bins 
       INTEGER ifl ! never used -- must be something standard to xspec
       REAL ear(0:ne)  ! energy array 
       REAL param(2)   ! param(1)=SI, param(2)=break frequency
       REAL photar(ne) ! photar=photon array (the spectrum) 
!                      [photons/bin] not [keV/bin]
       REAL photer(ne) ! unused, needed for call to powerLaw

! Varibles unique to this model
       INTEGER NMOD ! number of points in the model spectrum

       REAL kevhz, kevfljy ! kevhz [keV]*kevhz=[Hz]
       PARAMETER (kevhz = 2.417988366E17, kevfljy = 1509.188961)
!  1kev = kevhz Hz,  1 keV/(s kev cm2) = kevfljy Jy (keV flux per Jy)
!  I think Jonathan's comments are wrong -- these conversion are:
!  1 Jy = kevfljy * 1 keV/(s kev cm2)
!  1 Hz = kevhz keV

       REAL logmin, logdelta 
! Parameters are constants -- these cannot change within the program now
!       Parameter(NMOD = 80) ! Number of bins in model
       Parameter(NMOD = 121) ! Number of bins in model

!       Parameter(logmin = -4.1, logdelta = 0.1) ! energy to start, increment
       Parameter(logmin = -9.28, logdelta = 0.1) ! energy to start, increment

       REAL a, b                ! [a]=[kevfljy * recGhz**alphan] , b is dimentionless and is related to the reduction
       REAL Ebr, Ebr2 ! break energy in keV, 
       REAL recGhz ! conversion from [keV] -> [GHz]
       REAL alphan ! radio spectral index, ALPHA Negative
       REAL gamma  ! photon spectral index
! gamma is used generally for both the energy spectral index and the photon spectral index -- jonathan's policy is to use lowercase alpha for energy spectral index and uppercase gamma for photon spectral index with gamma = alpha + 1

       REAL logR(NMOD)
       INTEGER i, j
       REAL  terp ! dummy interpolation variable for linear interpolation
       REAL t

!   10**logR is the factor by which the flux is reduced from a straight powerlaw
      DATA logR/
     & 0,0,0,0,0,0, 
     & -4.34316e-05, -4.34316e-05, -4.34316e-05, -4.34316e-05, 
     & -4.34316e-05, -4.34316e-05, -8.68676e-05, -8.68676e-05, 
     & -8.68676e-05, -0.000130308, -0.000130308, -0.000130308, 
     & -0.000173753, -0.000173753, -0.000217202, -0.000260655, 
     & -0.000260655, -0.000304113, -0.000347575, -0.000391041, 
     & -0.000434512, -0.000477987, -0.000521466, -0.000608438, 
     & -0.000651931, -0.000738929, -0.000825944, -0.0009565, 
     & -0.00104356, -0.00117418, -0.00130484, -0.00143554, 
     & -0.00165346, -0.00182788, -0.002046, -0.00226422, 
     & -0.00256993, -0.00287584, -0.00322573, -0.0035759, 
     & -0.00405784, -0.00449643, -0.00506726, -0.00568285, 
     & -0.00638744, -0.00713733, -0.00802121, -0.00899556, 
     & -0.010061, -0.0113519, -0.012736, -0.0140798, 
     & -0.015923, -0.017955, -0.0202241, -0.0223678, 
     & -0.0253041, -0.0285386, -0.0314704, -0.0357404, 
     & -0.040434, -0.0447453, -0.0508052, -0.0567529, 
     & -0.0635868, -0.0714505, -0.0804512, -0.0899629, 
     & -0.101001, -0.113735, -0.126912, -0.142245, 
     & -0.159392, -0.178028, -0.198734, -0.222791, 
     & -0.248567, -0.278189, -0.310958, -0.348528, 
     & -0.390086, -0.437826, -0.492414, -0.553618, 
     & -0.622694, -0.700057, -0.784891, -0.880414, 
     & -0.9897, -1.09903, -1.21724, -1.34496, 
     & -1.48294, -1.6319, -1.79237, -1.99055, 
     & -2.17431, -2.37161, -2.58386, -2.81163, 
     & -3.05631, -3.31867, -3.59998, -3.90136, 
     & -4.22417, -4.56944, -4.9393, -5.3347, 
     & -5.75796, -6.21049, -6.69443, -7.21183, 
     & -7.76574, -8.35803, -8.99225/ 
 
!     &  0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01, 0.1000E+01,
!     &  0.1000E+01, 0.9999E+00, 0.9999E+00, 0.9999E+00, 0.9999E+00,
!     &  0.9999E+00, 0.9999E+00, 0.9998E+00, 0.9998E+00, 0.9998E+00,
!     &  0.9997E+00, 0.9997E+00, 0.9997E+00, 0.9996E+00, 0.9996E+00,
!     &  0.9995E+00, 0.9994E+00, 0.9994E+00, 0.9993E+00, 0.9992E+00,
!     &  0.9991E+00, 0.9990E+00, 0.9989E+00, 0.9988E+00, 0.9986E+00,
!     &  0.9985E+00, 0.9983E+00, 0.9981E+00, 0.9978E+00, 0.9976E+00,
!     &  0.9973E+00, 0.9970E+00, 0.9967E+00, 0.9962E+00, 0.9958E+00,
!     &  0.9953E+00, 0.9948E+00, 0.9941E+00, 0.9934E+00, 0.9926E+00,
!     &  0.9918E+00, 0.9907E+00, 0.9897E+00, 0.9884E+00, 0.9870E+00,
!     &  0.9854E+00, 0.9837E+00, 0.9817E+00, 0.9795E+00, 0.9771E+00,
!     &  0.9742E+00, 0.9711E+00, 0.9681E+00, 0.9640E+00, 0.9595E+00,
!     &  0.9545E+00, 0.9498E+00, 0.9434E+00, 0.9364E+00, 0.9301E+00,
!     &  0.9210E+00, 0.9111E+00, 0.9021E+00, 0.8896E+00, 0.8775E+00,
!     &  0.8638E+00, 0.8483E+00, 0.8309E+00, 0.8129E+00, 0.7925E+00,
!     &  0.7696E+00, 0.7466E+00, 0.7207E+00, 0.6928E+00, 0.6637E+00,
!     &  0.6328E+00, 0.5987E+00, 0.5642E+00, 0.5270E+00, 0.4887E+00,
!     &  0.4482E+00, 0.4073E+00, 0.3649E+00, 0.3218E+00, 0.2795E+00,
!     &  0.2384E+00, 0.1995E+00, 0.1641E+00, 0.1317E+00, 0.1024E+00,
!     &  0.7961E-01, 0.6064E-01, 0.4519E-01, 0.3289E-01, 0.2334E-01,
!     &  0.1613E-01, 0.1022E-01, 0.6694E-02, 0.4250E-02, 0.2607E-02,
!     &  0.1543E-02, 0.8784E-03, 0.4801E-03, 0.2512E-03, 0.1255E-03,
!     &  0.5968E-04, 0.2695E-04, 0.1150E-04, 0.4627E-05, 0.1746E-05,
!     &  0.6159E-06, 0.2021E-06, 0.6140E-07, 0.1715E-07, 0.4385E-08,
!     &  0.1018E-08/        
 




      recGhz  =  1E-9 * kevhz
      alphan = -1.0*param(1) !alpha, negative, spectral index
      Ebr = param(2)/kevhz
      Ebr2 = 2*Ebr
   
C -- log10(F) = logmin + i * logdelta  where i = 1, 2, 3 etc ...

      gamma = param(1) + 1 ! photon spectral index
      a = kevfljy * recGhz**alphan ! normalizing
      b = 0.0

      CALL xspwlw(ear, ne, gamma, ifl, photar, photer) ! call to power law spectra


!      PRINT *, 'a = ', a
!      PRINT *, 'Ebr = ', Ebr
!      PRINT *, 'alphan = ', alphan
        
        DO j = 1, ne
          t = (log10( (ear(j-1)+ear(j))/Ebr2 ) - logmin)/logdelta
!          write (*,*) "---->", ear(j-1)+ear(j), logmin, 
!     &         (log10( (ear(j-1)+ear(j))/Ebr2 ) - logmin),
!     &         logdelta, 
!     &         (log10( (ear(j-1)+ear(j))/Ebr2 ) - logmin)/logdelta
!      PRINT *,'t = ', t
          IF (t .gt. float(NMOD)) then 
            b = 0.0
!          PRINT *, 'b = ', b
          ELSEIF (t .lt. 1.0) then 
            b = 1.00
!          PRINT *, 'b = ', b
          ELSE
            i = int(t)
            terp = t - float(i)  
            b  = 10.0**( (1.0-terp)*logR(i) + terp*logR(i+1) )
!            PRINT *, 'b = ', b
          ENDIF

!          PRINT *, 'a = ', a

          photar(j) = a*b*photar(j)
!         PRINT *, 'photar(', j, ')', photar(j)
         
        ENDDO
           
        RETURN
        END
