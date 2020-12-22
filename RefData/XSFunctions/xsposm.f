
c******************************************************************************
c
c     posm: modified version of the orginal POSM function. The orginal
c           version computed a single-valued differential photon flux at
c           each input energy ear(i). this version does a simple trapezoid
c           integration over the ear(i-1)->ear(i) interval. Also, the 
c           E<=511 was replaced w/E<511 to avoid log(0) overflows
c
c
c   C.Shrader, Code 661  NASA/GSFC  05/2005
c
c
c
c******************************************************************************
c

ch1    ROUTINE NAME:  xsposm
ch1
ch1    VERSION:  I              REVISION:
ch1
ch1    PROGRAMMER(S) AND COMPLETION DATE: Sandhia Bansal -   02/20/97
ch1
ch1    FUNCTION: Returns expected value of a given channel for the
ch1              POSITRONIUM_CONTINUUM model.  This module is taken from
ch1              TGRS Positronium_Continuum and modified for use with
ch1              XSPEC.
ch1
ch1    SOFTWARE SYSTEM AND SPACECRAFT: TGRS Data Analysis Software/WIND.
ch1
ch1    COMPUTER AND LANGUAGE: MicroVax/VAX FORTRAN
ch1
ch2    CALLING SEQUENCE: Xsps (Ear, Ne, Param, Ifl, Photar)
ch2
ch2    ARGUMENT        TYPE    I/O                  DESCRIPTION
ch2    ________        ____    ___      ________________________________________
ch2    bin_energy       real    I       Bin-energy.
ch2    n                real    I       Fit parameter.
ch2    ymodel           real    I       Calculated bin-value.
ch2    nchnl            I       I       Number of channels to be processed.
ch2
ch2    CALLED BY:  qfit, modfit, conf, spplot (process_model);
ch2                lsrch (create_lsrch_model).
ch2
ch2    CALLS:  --
ch3
ch3    COMMON USE:  --
ch3
ch3    SIGNIFICANT LOCAL VARIABLES:
ch3    VARIABLE    TYPE    INI. VAL.            DESCRIPTION
ch3    ________    ____    _________  __________________________________________
ch3
ch3
ch4    LOGICAL UNITS USED:    UNIT #            DESCRIPTION
ch4                           ______  __________________________________________
ch4
ch4
ch4    METHOD:
ch4    PDL
ch4      PROC POSITRONIUM_CONTINUUM  {compute the expected value of a given
ch4                                   channel}
ch4      set E0 to 511
ch4      set NORM to 2.0/((3.14159**2 - 9.0)*E0)
ch4
ch4      DO FOR those channels for which X < E0
ch4         substitute the bin energy and fit parameters in the following
ch4                    equation:
ch4         OUTPUT = N*NORM*(X*(E0-X)/(2*E0-X)**2 +
ch4                          (2*E0*(E0-X)/X**2)*ALOG((E0-X)/E0) -
ch4                          (2*E0*(E0-X)**2/(2*E0-X)**3)*ALOG((E0-X)/E0) +
ch4                          (2*E0-X)/X)
ch4         compute the expected value of that bin
ch4      ENDDO
ch4
ch4      DO FOR rest of the channels
ch4         set OUTPUT to 0
ch4      ENDDO
ch4
ch4      END POSITRONIUM_CONTIMUUM
ch4
ch4    MODIFICATIONS BETWEEN VERSIONS:
ch5    MOD. #   MODIFIER     DATE                   DESCRIPTION
ch5    ______  __________  ________  ___________________________________________
c
c
c_______________________________________________________________________________
c
c
        subroutine xsposm (ear, ne, param, ifl, photar, photer)
c
c  ****************************************************************************
c  *  Given the bin energy and the fit parameters, this module returns        *
c  *  the expected value of that bin for the POSITRONIUM_CONTINUUM model.     *
c  ****************************************************************************
c
        implicit none
c
        integer    ne, ifl, ie
        real       ear(0:ne), param, photar(ne), photer(ne)
        real       e0, norm, e0mx, te0mx, elog, photar0
c
        parameter (e0 = 511.)
        parameter (norm = 2.0/((3.14159**2 - 9.0)*e0))
c
c..............................................................................
c  Compute the bin-energy for the POSITRONIUM_CONTINUUM model.
c..............................................................................
c
c... The equation is valid only for those channels for which the
c...     ear is less than e0 (511.0)
c

c suppress a warning message from the compiler
        ie = ifl
        ie = INT(param)

c this model does not calculate errors
        DO ie = 1, ne
           photer(ie) = 0.0
        ENDDO


c handle zeroth case 

           if ( ear(0) .lt. e0 ) then
              e0mx  = e0 - ear(0)
              te0mx = 2.0 * e0 - ear(0)
              elog  = log(e0mx/e0)
              photar0  = norm * (ear(0)*e0mx/
     *                              te0mx**2 +
     *                              (2.*e0*e0mx/ear(0)**2) *
     *                              elog -
     *                              (2.*e0*e0mx**2/te0mx**3) *
     *                              elog +
     *                              te0mx/ear(0))
           else
              photar0 = 0.0
           endif

c then, same as before for ear(1)->ear(ne):
c
       do ie = 1, ne
           if ( ear(ie) .lt. e0 ) then
              e0mx  = e0 - ear(ie)
              te0mx = 2.0 * e0 - ear(ie)
              elog  = log(e0mx/e0)
              photar (ie) = norm * (ear(ie)*e0mx/
     *                              te0mx**2 +
     *                              (2.*e0*e0mx/ear(ie)**2) *
     *                              elog -
     *                              (2.*e0*e0mx**2/te0mx**3) *
     *                              elog +
     *                              te0mx/ear(ie))
           else
              photar(ie) = 0.0
           endif
        enddo
c
c now do trapezoid integration:

        do ie = ne, 2, -1 
            photar(ie) = photar(ie) + photar(ie-1)
            photar(ie) = photar(ie)/2. * (ear(ie) - ear(ie-1))
            if (ear(ie) .ge. e0) photar(ie) = 0.
        enddo
c
        photar(1) = photar(1) + photar0
        photar(1) = photar(1)/2. * (ear(1)-ear(0))
c
        return
        end
