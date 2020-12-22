// c
// c     driver for the ntee nonthermal/thermal/reflection program
// c            see MULMOD for parameter descriptions
// c
// c     number of model parameters: 15
// c     parameters are given below in the order of decreasing importance:
// c     1: nonthermal electron compactness
// c     2: blackbody compactness
// c     3: scaling factor for reflection (1 for isotropic source above disk)
// c     4: blackbody temperature in eV
// c     5: the maximum Lorentz factor
// c     6: thermal compactness (0 for pure nonthermal plasma)
// c     7: Thomson optical depth of ionization electrons (e.g., 0)
// c     8: electron injection index (0 for monoenergetic injection)
// c     9: minimum Lorentz factor of the power law injection (not used for mono)
// c    10: minimum Lorentz factor for nonthermal reprocessing (>1; <= par. 9)
// c    11: radius in cm (for Coulomb/bremsstrahlung only)
// c    12: pair escape rate (0-1)
// c    13: cosine of inclination angle
// c    14: iron abundance relative to the solar iron abundance
// c    15: redshift
// c
// c     algorithm:
// c           a(x)=(non-thermal/thermal pair spectrum)+reflection
// c

#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <cmath>
#include <cfortran.h>

// declarations for other routines

void doreflect(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString, Real Emax);

// Called from C++ to F

PROTOCCALLSFSUB5(NONTH,nonth,FLOATV,FLOATV,PINT,FLOATV,FLOATV)
#define NONTH(pars, xnonth, Nnonth, spnth, sphth) \
   CCALLSFSUB5(NONTH,nonth,FLOATV,FLOATV,PINT,FLOATV,FLOATV, \
               pars, xnonth, Nnonth, spnth, sphth)

void xsnteea(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   using namespace XSutility;

   const size_t nE1 = energyArray.size();
   size_t nE = nE1-1;
   flux.resize(nE);

   string msg;
   static bool firstcall(true);
   if ( firstcall ) {
     msg = "NTEEA v1.1. WARNING: This is an experimental model";
     xs_write(const_cast<char*>(msg.c_str()),5);
     msg = "   Consult Andrzej Zdziarski for advice";
     xs_write(const_cast<char*>(msg.c_str()),5);
     firstcall = false;
   }
   
   static RealArray saveParams(0.0,15);
   bool recalc(false);
   for (size_t i=0; i<15; i++) {
     if ( params[i] != saveParams[i] ) recalc = true;
   }
   saveParams = params;

   // In order for the absorption and reflection calculations to be performed
   // correctly the energies on which the input spectrum are calculated should
   // run from 0.005 keV to 511.0/YMIN keV in the rest frame of the emitter.

   const Real YMIN = 0.03333;

   size_t NewBinsHigh = 0;
   size_t NewBinsLow = 0;
   Real dellogEHigh = 0.01;
   Real dellogELow = 0.01;

   Real zshift = 1.0 + params[14];
   Real RestEmin = 0.005/zshift;
   Real RestEmax = 511.0/YMIN/zshift;

   if ( energyArray[0] > RestEmin ) NewBinsLow = (size_t) ((log10(energyArray[0]) - log10(RestEmin)) * 100.0);

   if ( energyArray[nE] < RestEmax ) NewBinsHigh = (size_t) ((log10(RestEmax)-log10(energyArray[nE])) * 100.0);

   size_t NewBins = NewBinsLow + NewBinsHigh;
   size_t nEtot = nE + NewBins;
   size_t nEtot1 = nEtot + 1;

   // instantiate extended arrays

   RealArray ExtendEnergyArray(nEtot1);
   RealArray ExtendFlux(nEtot);
   RealArray ExtendFluxErr(nEtot);

   // and set extended energy array

   for (size_t i=0; i<NewBinsLow; i++) ExtendEnergyArray[i] = pow(10.0,log10(RestEmin)+i*dellogELow);
   for (size_t i=NewBinsLow; i<=nE+NewBinsLow; i++) ExtendEnergyArray[i] = energyArray[i-NewBinsLow];
   for (size_t i=nE+NewBinsLow+1; i<ExtendEnergyArray.size(); i++) 
     ExtendEnergyArray[i] = pow(10.0,log10(ExtendEnergyArray[nE+NewBinsLow])+(i-nE-NewBinsLow)*dellogEHigh);

   // set up internal arrays for NONTH

   const int MAXENE = 900;
   static float xnonth[MAXENE];
   static float spnth[MAXENE];
   static float sphth[MAXENE];
   float normfac(1.0);
   static int   Nnonth(1);

   // if necessary recalculate the thermal and non-thermal spectra

   if ( recalc ) {
     float pars[15];
     for (size_t i=0; i<15; i++) pars[i] = params[i];
     NONTH(pars, xnonth, Nnonth, spnth, sphth);
   }

   // set up array of energies in m_e c^2

   float *xtot=0;
   xtot = new float[nEtot];
   for (size_t i=0; i<nEtot; i++) xtot[i] = (1+params[14])*(ExtendEnergyArray[i]+ExtendEnergyArray[i+1])/2.0/511.0;

   // check for grid non-overlaps

   if ( xtot[0] < xnonth[0] ) {
     msg = "ntee: input energy is below valid energy for model";
     xs_write(const_cast<char*>(msg.c_str()),5);
     return;
   }
   if ( xtot[nEtot-1] > xnonth[Nnonth-1] ) {
     msg = "ntee: input energy is above valid energy for model";
     xs_write(const_cast<char*>(msg.c_str()),5);
     return;
   }

   // Map thermal and non-thermal spectra back onto the input energy array

   RealArray nonthSpectrum(0.0,nEtot);
   RealArray thSpectrum(0.0,nEtot);

   for (size_t j=0; j<nEtot; j++) {
     size_t i = 0;
     while ( xnonth[i] < xtot[j] ) i++;
     size_t il = i - 1;
     Real delta = xtot[j] - xnonth[il];
     Real spo = spnth[il]/xnonth[il]/xnonth[il];
     Real grad = (spnth[i]/xnonth[i]/xnonth[i]-spo)/(xnonth[i]-xnonth[il]);
     nonthSpectrum[j] = (spo+grad*delta)*xtot[j]*xtot[j];
     spo = sphth[il]/xnonth[il]/xnonth[il];
     grad = (sphth[i]/xnonth[i]/xnonth[i]-spo)/(xnonth[i]-xnonth[il]);
     thSpectrum[j] = (spo+grad*delta)*xtot[j]*xtot[j];
   }
  
   // convert to photons and integrate over the bin size,

   for (size_t j=0; j<nEtot; j++) {
     Real emid = (ExtendEnergyArray[j] + ExtendEnergyArray[j+1])/2.0;
     Real esize = ExtendEnergyArray[j+1]-ExtendEnergyArray[j];
     nonthSpectrum[j] *= esize/emid/emid;
     thSpectrum[j] *= esize/emid/emid;
   }

   // calculate reflection if required

   ExtendFlux = nonthSpectrum;

   Real scdir(1.0);
   Real scref(params[2]);
   if ( scref > 0 ) {

     RealArray par_refl(7);

     par_refl[0] = -1.0;              // ensures only reflection returned
     par_refl[1] = params[14];        // Redshift
     par_refl[2] = 1.0;               // Abundances heavier than He
     par_refl[3] = params[13];        // Fe abundance
     par_refl[4] = params[12];        // cosIncl
     par_refl[5] = 0.0;
     par_refl[6] = 0.0;

     msg = "Recalculating reflection...";
     xs_write(const_cast<char*>(msg.c_str()),20);

     //  scale=params[2]=1 for 2pi solid angle
     //  the scaling factor; 1 corresponds to seeing equal
     //  contributions from the reflected and direct spectra
     //  the scaling for the direct component is 1/(1+scale)
     //  the scaling for the reflected component is scale/(1+scale)
     //  the sum is unity: conservation of energy
     //    scdir=1/(1+scale)
     //    scref=parascale/(1+scale)
     //
     //  The above assumptions cause problems in xspec. The normalization cannot
     //  change from one iteration to another. Therefore we keep constant the
     //  normalization of the direct component as below:    (29.10.1994)
     //   scdir = 1
     //   scref = Scale

     // calculate the power in the scaled incident spectrum

     Real clnth = 0.0;
     for (size_t i=0; i<ExtendFlux.size(); i++) {
       clnth += ExtendFlux[i] * (ExtendEnergyArray[i]+ExtendEnergyArray[i+1])/2.0;
     }
     clnth *= scdir;

     // calculate the reflection spectrum (only)

     string modelName = "NTEEA";

     doreflect(ExtendEnergyArray, par_refl, 1, ExtendFlux, fluxErr, modelName, ExtendEnergyArray[nE+NewBinsLow]);

     // calculate the scaled reflected power

     Real clref = 0.0;
     for (size_t i=0; i<ExtendFlux.size(); i++) {
       clref += ExtendFlux[i] * (ExtendEnergyArray[i]+ExtendEnergyArray[i+1])/2.0;
     }
     clref *= scref;

     // the integrated albedo and the scaling for the thermal component

     Real albedo = clref/clnth/scref;
     Real cadd = (1-albedo)*(clnth/params[1])*scref*(4.0*M_PI/3.0);

     // scale the thermal and reflection components

     thSpectrum *= cadd;
     ExtendFlux *= scref/scdir;

   } else {
     ExtendFlux = 0.0;
     thSpectrum = 0.0;
   }

   // calculate the normalization factor - the thermal + non-thermal spectra
   // at 1 keV in the Earth frame.

   size_t i1 = 0;
   while (ExtendEnergyArray[i1] < 1.0) i1++;
   i1--;
   normfac = (ExtendEnergyArray[i1]-ExtendEnergyArray[i1-1])/(nonthSpectrum[i1]+thSpectrum[i1])/((ExtendEnergyArray[i1]+ExtendEnergyArray[i1-1])/2.0);

   // sum all the components and normalize the thermal and non-thermal spectra

   flux.resize(nE);
   for (size_t i=0; i<nE; i++) {
     i1 = i + NewBinsLow;
     flux[i] = ExtendFlux[i1] + nonthSpectrum[i1] + thSpectrum[i1];
   }
   flux *= normfac;
   return;

}



