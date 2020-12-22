// C++ wrapper for compps model
// c*********************************************************************
// c     number of model parameters: 19 
// c   1: Te, electron temperature in keV
// c   2: p,  electron power-law index [ N(gamma)=gamma^-p ] 
// c   3: gmin,  minimum Lorentz factor gamma 
// c   4: gmax,  maximum Lorentz factor gamma
// c      (a) if any of gmin or gmax < 1 then 
// c          Maxwellian electron distribution (IDISTF=0) 
// c          with parameter Te 
// c      (b) if Te=0.  then 
// c         power law electrons (IDISTF=2) 
// c          with parameters p, gmin, gmax
// c       (c) if both gmin,gmax>=1 but gmax<gmin then
// c           cutoff Maxwellian (IDISTF=3) 
// c          with Te, p, gmin (cutoff Lorentz factor) as parameters 
// c      (d) if Te.ne.0, gmin, gmax >=1  then 
// c         hybrid electron distribution (IDISTF=4) 
// c          with parameters Te, p, gmin, gmax 
// c 
// c   5: Tbb, temperature of soft photons 
// c         Tbb>0 blackbody 
// c         Tbb<0 multicolor disk with inner disk temperature Tbb 
// c   6: if > 0 tau, vertical optical depth of the corona 
// c      if < 0 y=4*Theta*tau 
// c       limits: for the slab geometry - tau<1 
// c             if say tau~2 increase MAXTAU to 50
// c            for sphere tau<3 
// c   7: geom,  IGEOM=INT(geom) 
// c      IGEOM=0 approximate treatment of radiative transfer using 
// c      escape probability for a sphere (very fast method) 
// c      IGEOM=1 slab; IGEOM=2 cylinder; IGEOM=3 hemisphere; IGEOM=4,5 sphere
// c      input photons at the bottom of the slab, cylinder, hemisphere 
// c      or center of the sphere 
// c      (or from the central plane of the slab if cov_fact ne 1) 
// c 
// c      if IGEOM<0 then geometry defined by ABS(IGEOM) and 
// c      sources of incident photons are isotropic and homogeneous 
// c      IGEOM=-5 - sphere with the source of photons distributed 
// c      according to the eigen function of the diffusion equation
// c      f(tau')=sim(pi*tau'/tau)/(pi*tau'/tau) 
// c      where tau' varies between 0 and tau. 
// c
// c   8: H/R for cylinder geometry only
// c   9: cosIncl, cosine of inclination angle 
// c      (if < 0 then only black body) 
// c  10: cov_fac, covering factor of cold clouds 
// c               if IGEOM=\pm 4,5 then cov_fac is dummy
// c  11: R, amount of reflection Omega/(2*pi) 
// c      (if R < 0 then only reflection component) 
// c  12: FeAb, iron abundance in units of solar 
// c  13: MeAb, abundance of heavy elements in units of solar
// c  14: xi, disk ionization parameter L/(nR^2)
// c  15: temp, disk temperature for reflection in K
// c  16: beta, reflection emissivity law (r^beta)
// c      if beta=-10 then non-rotating disk
// c      if beta=10  then 1.-sqrt(6./rg))/rg**3 
// c  17: Rin/Rg,  inner radius of the disk (Schwarzschild units)
// c  18: Rout/Rg, outer radius of the disk
// c  19: redshift 
// c
// c*********************************************************************
// c     algorithm:
// c       exact solution of the radiative transfer equation
// c      for Compton scattering 
// c       in a HOT CORONA using iterative scattering method 
// c       (Poutanen J., \&  Svensson R.: 1996, 
// c       The Two-Phase Pair Corona Model for  Active Galactic Nuclei and
// c       X-ray Binaries: How to Obtain Exact Solutions, ApJ, 470, 249-268)
// c 
// c       Compton reflection is NOT scattered 
// c       Energy balance is neglected 
// C normalization (version 3.6) using SEED (INTRINSIC) black body emission
// C       n(E) = K 1.0344E-3  E**2 dE / (exp(E/kt)-1)
// C       The above normalization is such that K = (RKM)**2 /(D10)**2
// C       that is, the norm is the source radius squared (assuming a spherical
// C       surface in units of km, i.e. S=4pi R**2)
// C       divided by the distance to the source (in units of 10 kpc) squared.
// c********************************************************************* 
// c Version 4.00 :  October 18, 2011
// c********************************************************************* 


#include <xsTypes.h>
#include <XSstreams.h>
#include <stlToCArrays.h>
#include <functionMap.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <cmath>
#include <cfortran.h>

// declarations for other routines

void doreflect(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString, Real Emax);

// Functions called FROM here to Fortran:

PROTOCCALLSFSUB8(ISMCO,ismco,DOUBLEV,DOUBLEV,DOUBLEV,DOUBLEV,DOUBLEV,INT,PDOUBLE,LOGICAL)
#define ISMCO(parm,obj,phcon,phblb,phref,ne,phnorm,flagfre)		\
  CCALLSFSUB8(ISMCO,ismco,DOUBLEV,DOUBLEV,DOUBLEV,DOUBLEV,DOUBLEV,INT,PDOUBLE,LOGICAL, \
	      parm,obj,phcon,phblb,phref,ne,phnorm,flagfre)


void xscompps(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   using namespace XSutility;

   // static variables which need to persist between invocations

   static RealArray saveEnergyArray(0.0,0);

   const int nE = static_cast<int>(energyArray.size()) - 1;

//write the start-up message

   static bool firstcall(true);
   string msg;

   if ( firstcall ) {
     msg = "CompPS Version 4.00";
     xs_write(const_cast<char*>(msg.c_str()),5);
     msg = "Comptonization by Iterative Scattering Method";
     xs_write(const_cast<char*>(msg.c_str()),5);
     msg = "Poutanen & Svensson 1996";
     xs_write(const_cast<char*>(msg.c_str()),5);
     msg = "Questions: Juri Poutanen (juri.poutanen@oulu.fi)";
     xs_write(const_cast<char*>(msg.c_str()),5);
     firstcall = false;
   }

// check whether energies have changed

  bool energiesChanged(false);
  if ( (size_t)nE+1 != saveEnergyArray.size() ) {
    energiesChanged = true;
    saveEnergyArray.resize(nE+1);
  } else {
    size_t i = 0;
    while ( energyArray[i] == saveEnergyArray[i] && i <(size_t)nE ) i++;
    if ( i != (size_t)nE ) energiesChanged = true;
  }
  saveEnergyArray = energyArray;

  // In order for the absorption and reflection calculations to be performed
  // correctly the energies on which the input spectrum are calculated should
  // run from 0.005 keV to 511.0/YMIN keV in the rest frame of the emitter.

  const Real YMIN = 0.03333;

  int NewBinsHigh = 0;
  int NewBinsLow = 0;
  Real dellogEHigh = 0.01;
  Real dellogELow = 0.01;

  Real zshift = 1.0 + params[18];
  Real RestEmin = 0.005/zshift;
  Real RestEmax = 511.0/YMIN/zshift;

  if ( energyArray[0] > RestEmin ) NewBinsLow = (int) ((log10(energyArray[0]) - log10(RestEmin)) * 100.0);

  if ( energyArray[nE] < RestEmax ) NewBinsHigh = (int) ((log10(RestEmax)-log10(energyArray[nE])) * 100.0);

  int NewBins = NewBinsLow + NewBinsHigh;
  int nEtot = nE + NewBins;
  int nEtot1 = nEtot + 1;

  // instantiate extended arrays

  RealArray ExtendEnergyArray(nEtot1);
  RealArray ExtendFlux(nEtot);
  RealArray ExtendFluxErr(nEtot);

  // and set extended energy array

  for (int i=0; i<NewBinsLow; i++) ExtendEnergyArray[i] = pow(10.0,log10(RestEmin)+i*dellogELow);
  for (int i=NewBinsLow; i<=nE+NewBinsLow; i++) ExtendEnergyArray[i] = energyArray[i-NewBinsLow];
  for (int i=nE+NewBinsLow+1; i<(int)ExtendEnergyArray.size(); i++) 
    ExtendEnergyArray[i] = pow(10.0,log10(ExtendEnergyArray[nE+NewBinsLow])+(i-nE-NewBinsLow)*dellogEHigh);


   // set parameters required by ismco. these are the first 10 model parameters
   // plus the redshift
   double pac[11];

   for (size_t i=0; i<10; i++) {
     pac[i]=params[i];
   }
   pac[10] = params[18];

   static double *phcon=0, *phblb=0, *phref=0, *obj=0;
   if ( energiesChanged ) {
     delete [] phcon;
     phcon = new double[nEtot];
     delete [] phblb;
     phblb = new double[nEtot];
     delete [] phref;
     phref = new double[nEtot];
     delete [] obj;
     obj   = new double[nEtot];
   }

   // object's frame energies (in mc2)

   for (int i=0; i<nEtot; i++) obj[i] = (1.+params[18])*(ExtendEnergyArray[i]+ExtendEnergyArray[i+1])/2.0/511.0;

   double phnorm = 0.0;

   ISMCO(pac, obj, phcon, phblb, phref, nEtot, phnorm, energiesChanged);

   // convert output spectra back into xspec units

   RealArray spcon(0.0,nEtot);
   RealArray spblb(0.0,nEtot);
   RealArray spref(0.0,nEtot);
   for (int i=0; i<nEtot; i++) {
     Real factor = (ExtendEnergyArray[i+1]-ExtendEnergyArray[i])/phnorm;
     spcon[i] = phcon[i] * factor;
     spblb[i] = phblb[i] * factor;
     spref[i] = phref[i] * factor;
   }

   // set the parameters for the reflection model

   Real ReflScale = fabs(params[10]);

   if ( ReflScale > 0.0 ) {

     RealArray par_refl(7);

     par_refl[0] = -1.0;              // ensures only reflection returned
     par_refl[1] = params[18];        // Redshift
     par_refl[2] = params[12];        // Abundances heavier than He
     par_refl[3] = params[11];        // Fe abundance
     par_refl[4] = params[8];         // cosIncl
     par_refl[5] = params[14];        // Disk temperature
     par_refl[6] = params[13];        // Ionization parameter

     msg = "Recalculating reflection...";
     xs_write(const_cast<char*>(msg.c_str()),20);

     // calculate the reflection spectrum (only)

     string modelName = "COMPPS";

     doreflect(ExtendEnergyArray, par_refl, 1, spref, ExtendFluxErr, modelName, ExtendEnergyArray[nE+NewBinsLow]);

     // now follows rotational broadening in case param[16] is not -10. Uses the rdblur
     // convolution model which is based on diskline

     if ( params[15] != -10 ) {

       msg = "Recalculating blurring...";
       xs_write(const_cast<char*>(msg.c_str()),20);
     
       // set up the rdblur parameters
     
       RealArray par_blur(4);
       par_blur[0] = params[15];                       // Beta
       par_blur[1] = params[16];                       // Rin
       par_blur[2] = params[17];                       // Rout
       par_blur[3] = acos(params[8])*180.0/3.14159265;// Inclination

       rdblur(ExtendEnergyArray, par_blur, 1, spref, ExtendFluxErr, modelName);

     }

   }

   // Sum the total spectrum as required and map back onto original array

   flux.resize(nE);
   for (int i=0; i<nE; i++) {
     int i1 = i + NewBinsLow;
     if (params[10] < 0.) flux[i] = spref[i1];
     if (params[8] < 0.) flux[i] = spblb[i1];
     if (params[8] < 0. && params[10] < 0.) flux[i]= spblb[i1]+spref[i1]*ReflScale; 
     if (params[10] >= 0. && params[8] >= 0.) flux[i] = spcon[i1]+spref[i1]*ReflScale+spblb[i1];
   }

}

