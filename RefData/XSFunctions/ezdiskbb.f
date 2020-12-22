
      SUBROUTINE EZDISKBB (EAR, NE, PARAM, IFL, PHOTAR, PHOTER)

        IMPLICIT NONE

      Integer IFL, NE
      Real EAR(0:NE), PARAM(*), PHOTAR(NE), PHOTER(NE)      


c      PROGRAM "ezdiskbb.f"
c
c      Written by:       Erik R. Zimmerman
c                  Harvard-Smithsonian Center for Astrophysics
c                  60 Garden St., Cambridge, MA 02138
c                      
c      Version:       July, 2004
c
c      
c      This program computes the Newtonian multi-temperature blackbody 
c      spectrum of a thin, Keplerian accretion disk around a black hole or
c      neutron star, assuming zero torque at the inner boundary of the
c      disk.  The temperature profile of the disk that is assumed here is:
c
c      T(R/R_in) = 2.05 T_max (R/R_in)^(-3/4) (1 - (R/R_in)^(-1/2))^(-3/4)
c
c      This model has one parameter: T_max, the maximum temperature in the
c      disk.  
c      
c      This model is compatible with XSPEC versions 11 and 12
c      and is meant to serve as an alternative version to the XSPEC model
c      "diskbb," which assumes a nonzero torque at the inner boundary
c      of the disk.  The name of the model in XSPEC corresponding to this
c      code is "ezdiskbb."
c
c      INPUTS:
c
c      (For more details on these inputs, see Appendix C of the XSPEC12 
c      manual at http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/xspec12
c      
c      -EAR(0:NE): Array of observed energy bins in units of keV
c      -NE: Number of observed energy bins
c      -PARAM(*): Array of input parameters from XSPEC.  Ezdiskbb has
c            only one parameter: T_max, the maximum temperature of
c            the accretion disk in units of keV
c      -IFL: The number of the spectrum being calculated
c      
c      OUTPUTS:
c
c      -PHOTAR(NE): Array containing the photon flux for each energy bin
c                  in units of photons/cm^2/second/bin.  PHOTAR(i),
c                  contains the flux for the energy bin
c                  between EAR(i-1) and EAR(i)
c      -PHOTER(NE): Optional output containing the flux errors.  This
c                  array is not calculated by this program.


      Integer ILO(0:25000), IHI(0:25000), I, J, M, IREAD, ilun
      Integer ierr, index

      Real FLUX(0:24999), E(0:24999), SLOPE, FIN(0:25000),
     +            F1, F2, F, FNEW(0:24999), ENEW(0:24999), Z

      character(256) datdir, filenm, contxt

      logical qanyf

      integer lenact
      character(256) fgmodf
      external lenact, fgmodf

c      A pre-calculated spectrum for T_max = 1 keV is contained in the
c      file "ezdiskbbdata.fits" with the corresponding energy bins.
c      The data from this files are read by the code once, and the values for 
c      the fluxes and energies are then saved in the variables
c      "flux" and "e," respectively.  The value of "iread" is then set to 99
c       so that the data file will not be read again. 

      save iread, e, flux
      data iread/0/

c suppress a warning message from the compiler
      i = ifl

c this model does not calculate errors
      DO i = 1, ne
         photer(i) = 0.0
      ENDDO
      
      if (iread .ne. 99) then

         datdir = fgmodf()
         filenm = datdir(:lenact(datdir))//'ezdiskbbdata.fits'
         CALL xwrite (' opening and reading ezdiskbb data file ', 10)

         CALL getlun(ilun)
         call ftopen(ilun, filenm, 0, index, ierr)
         call ftmnhd(ilun, 2, 'FLUX', 1, ierr)
         call ftgcve(ilun, 1, 1, 1, 25000, " ", e, qanyf, ierr)
         call ftgcve(ilun, 2, 1, 1, 25000, " ", flux, qanyf, ierr)
         call ftclos(ilun, ierr)
         CALL frelun(ilun)

         if ( ierr .NE. 0 ) then
            contxt = 'EZDISKBB: Failed to read '//
     &               filenm(:MAX(lenact(filenm),225))
            CALL xwrite(contxt, 5)
            RETURN
         else
            iread = 99
         endif
      
      endif
      
c      Using the input parameter, T_max, the new spectrum is calculated from
c      the one in which T_max = 1 keV.  This is done by multiplying all of the
c      values for the fluxes by (T_max)^2 and the values of the energies by
c      T_max.  The new values for the fluxes are saved in "fnew," and the new
c      energies are saved in "enew."
      
      do 30 i = 0,24999
      
            fnew(i) = (param(1)*param(1))*flux(i)
            enew(i) = (param(1))*e(i)
            
30       enddo
    

c      The energy bin in enew that surrounds each value of ear(i) is found 
c      algebraically.  "z" is a real number that corresponds to the "index" 
c      in the enew array of each value of ear(i).  z is then converted to 
c      an integer, "j," and then for each value of i, two indices are assigned: 
c      "ilo" and "ihi."  enew(ilo(i)) will be the largest energy in enew that 
c      is less than ear(i), and enew(ihi(i)) will therefore be the smallest 
c      energy in enew that is greater than ear(i).
c
c      If ear(i)/param(1) is greater than 10^4 or less than 10^-4, then the
c      index of j will be outside of the range of enew and fnew.  This should 
c      not occur in practice, because this model
c      is intended for energies of 0.01 - 100 keV and values of param(1) of
c      0.01 - 100 keV, meaning that ear(i)/param(1) should never exceed
c      10^4 or 10^-4 if the model is used properly.  If this error does occur,
c      then the program assigns a value to ilo(i) of -1.
c      For the rest of the calculations, this value of i is ignored,
c      until the end of the program, when the flux of bin ear(i) is 
c      assigned a value of zero.
      
      do 40 i = 0,NE
            
            z = 3125*(4 + log10(EAR(i)/param(1)))      

            j = int(z)
            
            if (j .LT. 0) then
                  j = -1
            endif
            
            if (j .GT. 24999) then
                  j = -1
            endif
            
            ilo(i) = j
            ihi(i) = j+1
            
40       enddo
           
c      The slope of the spectrum between enew(ilo(i)) and enew(ihi(i)) is 
c      calculated, and then an estimate of the flux at each value of ear(i)
c      is placed in the array fin(i).      
      
      do 60 i = 0 , NE
            
            if (ilo(i) .NE. -1) then
            
            slope = (fnew(ihi(i)) - fnew(ilo(i)))/
     +                  (enew(ihi(i)) - enew(ilo(i)))
     
                 fin(i) = fnew(ilo(i)) + slope*
     +                  (ear(i)-enew(ilo(i)))
     
                 endif
                                                            
60      enddo


c      The fluxes are integrated over the energy bins using the trapezoid
c      rule.  There are four possible cases that are taken into account.  
c      The first case is if ihi(i) is less than ilo(i+1), meaning that 
c      there are more than zero bins in enew between ear(i) and ear(i+1).
c      In this case, the flux is integrated between ear(i) and ihi(i), then
c      between ilo(i+1) and ear(i+1), and then between ihi(i) and ilo(i+1).
c      These fluxes are then summed and placed in photar(i+1).
c
c      The second case is if ihi(i) equals ilo(i+1), meaning that there are
c      zero bins in enew between ear(i) and ear(i+1).  In this case, the flux
c      is integrated between ear(i) and ihi(i), and then between ilo(i+1) and
c      ear(i+1).  These fluxes are then summed and placed in photar(i+1).
c
c      The fourth case is if an enew bin completely surrounds an ear bin,
c      meaning that ilo(i) would equal ilo(i+1).  In this case, the flux is
c      integrated between ear(i) and ear(i+1) and placed in photar(i+1).
c
c      The final case is if, as mentioned above, ear(i) falls outside of the
c      range of enew.  In this case, the value of photar(i) is set to zero.
c
c      In the first three cases, the difference in the energies is multiplied
c      by 1.602E-9 to convert the energies from keV to ergs, because fin
c      and fnew are in units of photons/cm^2/s/erg.  The result is the flux
c       in units of photons/cm^2/s/bin.

      do 70 i = 0, NE - 1
      
            if (ilo(i) .NE. -1) then
            
            if (ihi(i) .LT. ilo(i+1)) then
            
                  f1 = 0.5*(fin(i)+fnew(ihi(i)))*
     +                  (enew(ihi(i))-ear(i))*1.602177E-9
                  f2 = 0.5*(fnew(ilo(i+1))+fin(i+1))*
     +                  (ear(i+1)-enew(ilo(i+1)))*1.602177E-9
                  
                  f = 0.0
                  
                  do 80 m = ihi(i),ilo(i+1)-1
                        f = f + 0.5*(fnew(m)+fnew(m+1))*
     +                        (enew(m+1)-enew(m))*1.602177E-9
80                  enddo

                  photar(i+1) = f1 + f2 + f
            
            elseif (ihi(i) .EQ. ilo(i+1)) then
                  
                  f1 = 0.5*(fin(i)+fnew(ihi(i)))*
     +                  (enew(ihi(i))-ear(i))*1.602177E-9
                  f2 = 0.5*(fnew(ilo(i+1))+fin(i+1))*
     +                  (ear(i+1)-enew(ilo(i+1)))*1.602177E-9
                  
                  photar(i+1) = f1 + f2
            
            else
            
                  photar(i+1) = 0.5*(fin(i)+fin(i+1))*
     *                  (ear(i+1)-ear(i))*1.602177E-9

            endif
            
            else
            
                  photar(i+1) = 0
            
            endif
      
70      enddo
      
      RETURN
      END
      

