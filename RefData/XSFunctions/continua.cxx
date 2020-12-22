#include <xsTypes.h>
#include <memory>
#include <cfortran.h>

// Functions called FROM here to Fortran:

PROTOCCALLSFSUB4(EMISC,emisc,FLOAT,FLOATV,FLOATV,FLOATV)
#define EMISC(t,ei,c,freq) \
            CCALLSFSUB4(EMISC,emisc,FLOAT,FLOATV,FLOATV,FLOATV,t,ei,c,freq)

PROTOCCALLSFSUB11(ZNIBIN,znibin,INT,FLOATV,INT,FLOATV,FLOAT,PINT,INTV,INTV,FLOATV,FLOATV,FLOATV)
#define ZNIBIN(numinbins,ein,numoutbins,eout,z,numcmbins,indexin,zstart,wl,wr,ecmb) \
            CCALLSFSUB11(ZNIBIN,znibin,INT,FLOATV,INT,FLOATV,FLOAT,PINT,INTV,INTV,FLOATV,FLOATV,FLOATV, \
            numinbins,ein,numoutbins,eout,z,numcmbins,indexin,zstart,wl,wr,ecmb)

PROTOCCALLSFSUB10(ZREBIN,zrebin,FLOATV,FLOATV,INT,INT,INTV,FLOATV,FLOATV,INTV,FLOATV,FLOATV)
#define ZREBIN(specin,specout,numoutbins,numcmbins,indexin,wl,wr,zstart,ein,area) \
            CCALLSFSUB10(ZREBIN,zrebin,FLOATV,FLOATV,INT,INT,INTV,FLOATV,FLOATV,INTV,FLOATV,FLOATV, \
            specin,specout,numoutbins,numcmbins,indexin,wl,wr,zstart,ein,area)

// Functions called from Fortran TO here:

void continua(const float t, const float* ei, const float* ear, const int ne, const float z, float* workphot);
FCALLSCSUB6(continua,CONTINUA,continua,FLOAT,FLOATV,FLOATV,INT,FLOAT,FLOATV)


void continua(const float t, const float* ei, const float* ear, const int ne, const float z, float* workphot)
{
   //    This subroutine calculates X-ray continua.

   //    Input:  t          -   electron temperature (K)
   //            ei         -   ionic concentrations by number 
   //                           (but could also be
   //                           emission integral of each ion)
   //                           Hydrogen and helium are completely ionized.
   //                           1:       H^+
   //                           2:       He^++
   //                           3-9:     C ions (including neutrals)
   //                           10-17:   N
   //                           18-26:   O
   //                           27-37:   Ne
   //                           38-50:   Mg
   //                           51-65:   Si
   //                           66-82:   S
   //                           83-103:  Ca
   //                           104-130: Fe
   //                           131-159: Ni
   //            ear         -  energy bins
   //            ne          -  number of energy bins
   //            z           -  redshift
   //    Output: workphot(ibin,ielem)    -  bin continuum fluxes for each 
   //                            element ielem (-1 for H^+, 0 for He^++,
   //                            1-10 for heavy elements)

   //            Translated from cont.f, Feb 2009 (C.Gordon)


   // nfreq = total number of continuum frequency points
   // nelem = 2 + number of heavy elements (for H^+, He^++, and heavy elements)
   // freq[nfreq] = grid of photon energies (keV)c
   const int NFREQ = 69;
   const int NELEM = 12;
   const float freq[] = {.0544,.069,.086,.1,.122,.123,.166,.167,.217,.218,
       .3,.33,.34,.36,.37,.488,.49,.522,.523,.666,.668,.706,.708,
       .87,.872,1.,1.156,1.158,1.359,1.365,1.582,1.75,1.77,
       1.949,1.97,2.04,2.05,2.39,2.44,2.66,2.68,3.223,3.225,
       3.48,3.5,5.128,5.13,5.43,5.5,6.1,6.9,7.7,8.82,8.83,9.27,9.28,
       10.27,10.29,10.78,10.8,13.6,17.1,21.5,27.7,34.2,43.0,54.1,
       68.1,85.8};

   // cont[ielem][ifreq] = continuum only xray spectrum (ergs/sec/keV)
   std::unique_ptr<float[]> apCont(new float[NELEM*NFREQ]);
   float* cont = apCont.get();
   for (int ielem=0; ielem<NELEM; ++ielem)
   {
      int offset = ielem*NFREQ;
      for (int ifreq=0; ifreq<NFREQ; ++ifreq)
         cont[ifreq+offset] = 0.0;
   }

   // Continuum emission
   const float keV = 1.60219e-9;
   EMISC(t, const_cast<float*>(ei), cont, const_cast<float*>(freq));
   for (int ielem=0; ielem<NELEM; ++ielem)
   {
      int offset = ielem*NFREQ;
      for (int ifreq=0; ifreq<NFREQ; ++ifreq)
         cont[ifreq+offset] /= (keV*freq[ifreq]);
   }

   // get memory for rebinning arrays - maximum size required is the sum of the
   // sizes of the input and output arrays.

   const int nmax = ne + NFREQ;
   std::unique_ptr<int[]> apZstart(new int[nmax]);
   std::unique_ptr<int[]> apIndexin(new int[nmax]);
   std::unique_ptr<float[]> apWl(new float[nmax]);
   std::unique_ptr<float[]> apWr(new float[nmax]);
   std::unique_ptr<float[]> apEcmb(new float[nmax]);
   std::unique_ptr<float[]> apArea(new float[nmax]);

   // initialize rebinning
   int numcmbins=0;
   ZNIBIN(NFREQ,const_cast<float*>(freq),ne+1,const_cast<float*>(ear),z,
        numcmbins,apIndexin.get(),apZstart.get(),apWl.get(),apWr.get(),
        apEcmb.get());

   for (int ielem=0; ielem<NELEM;++ielem)
      ZREBIN(&cont[ielem*NFREQ],&workphot[ielem*ne],ne,numcmbins,apIndexin.get(),
           apWl.get(),apWr.get(),apZstart.get(),const_cast<float*>(freq),apArea.get());
}
