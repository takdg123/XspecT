// c    written by Paolo  Coppi 
// c    changed by Juri to increase the speed (Oct 28, 1997) 
// c    relativistic reflection is added Oct 29, 1997
// c    parameters changed by MG on Nov 5, 1997
// c    diskbb replaced by diskpn by MG, Jun 28, 1998
// c    fits data load routines added by PC, Jul 23, 2001
// 
//       SUBROUTINE XSEQPAIR(Ear,Ne,Param,Ifl,Photar,Photer)
//       
//       IMPLICIT NONE
// 
//       INTEGER Ne , Ifl
//       REAL Param(20) , Ear(0:Ne) , Photar(Ne), Photer(Ne)
// 
// c
// c     driver for the nonthermal/thermal pair program
// c
// c     number of model parameters: 18
// c     parameters are given below in the order of decreasing importance:
// c     1: hard-to-soft compactness ratio
// c     2: blackbody compactness (if <0, uncomptonized black body component only)
// c     3: blackbody temperature in eV (if <0, use diskpn, T_{max}=abs(Par(3)) )
// c     4: non-thermal-to-hard compactness ratio (0 for pure thermal plasma)
// c     5: Thomson optical depth of ionization electrons (e.g., 0)
// c     6: source region radius in cm (for Coulomb/bremsstrahlung only)
// c     7: minimum Lorentz factor for power law injection (not used for mono)
// c     8: maximum Lorentz factor for power law injection 
// c     9: power law injection index (if <0, inject monoenergetically at #8)
// c     10 choose primary pair or electron injection (0=electron, 1=pair)
// c     11: cosIncl, cosine of inclination angle
// c     12, R, amount of reflection (if < 0 then only reflection component)
// c     13: FeAb, iron abundance in units of solar
// c     14: abundance of elements heavier than He relative to solar
// c     15: disk temperature (used for reflection)
// c     16: xi, ionization parameter
// c     17: beta, reflection emissivity law (r^beta) IF=-10 then 
// c         non-rotating disk, IF=10 then  (1.-sqrt(6./rg))/rg**3
// c     18: radin, inner radius of the disk in M
// c     19: radout, outer radius of the disk in M
// c     20: redshift
// c
// c    Note: if #3 OR #12 < 0, only show BB and/or Reflected components
// c
// c     algorithm:
// c           a(x)=(non-thermal/thermal pair spectrum)+reflection
// c

#include <xsTypes.h>
#include <XSstreams.h>
#include <stlToCArrays.h>
#include <functionMap.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <cmath>
#include <cfortran.h>
#include <gsl/gsl_integration.h>

// declaration for main wrap routine for eqpair, eqtherm and compth

void eqpwrap(const RealArray& energyArray, const RealArray& params,
             RealArray& flux, RealArray& fluxerr, int modelNumber);


void xseqpair(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
  eqpwrap(energyArray, params, flux, fluxErr, 1);
}




// declarations for other routines

void doreflect(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString, Real Emax);

void GLrebin(const RealArray& wkev, const RealArray& inSpectrum, const RealArray& outEnergy, RealArray& outSpectrum, Real Redshift, int NewBinsLow, int NewBinsHigh);

// Function called from here to Fortran:

PROTOCCALLSFSUB6(PAIRITER,pairiter,FLOATV,INT,INT,DOUBLEV,DOUBLEV,DOUBLEV)
#define PAIRITER(param,np,imodel,spcomp,spblbd,energy) \
  CCALLSFSUB6(PAIRITER,pairiter,FLOATV,INT,INT,DOUBLEV,DOUBLEV,DOUBLEV,\
	      param,np,imodel,spcomp,spblbd,energy)


const int NGAMBIN = 221;
const int NPARAM = 20;


void eqpwrap(const RealArray& energyArray, const RealArray& params,
             RealArray& flux, RealArray& fluxerror, int modelNumber)
{

// c    written by Paolo  Coppi 
// c    changed by Juri to increase the speed (Oct 28, 1997) 
// c    relativistic reflection is added Oct 29, 1997
// c    parameters changed by MG on Nov 5, 1997
// c    diskbb replaced by diskpn by MG, Jun 28, 1998
// c    fits data load routines added by PC, Jul 23, 2001
// c
// c     driver for the nonthermal/thermal pair program
// c     algorithm:
// c           a(x)=(non-thermal/thermal pair spectrum)+reflection
// c
// 
// c modelNumber input argument determines which model is calculated
// c     1:   eqpair  
// c     2:   eqtherm 
// c     3:   compth  

// c Note that compth model has parameters 1-4 different from other two but this
// c difference is mainly hidden in lower level routines

   using namespace XSutility;

   // static variables which need to persist between invocations

  static RealArray saveParams(0.0,NPARAM);
  static RealArray saveEnergyArray(0.0,0);

  static RealArray wkev(NGAMBIN);
  static RealArray compspect(NGAMBIN);
  static RealArray blbdspect(NGAMBIN);

  static int saveModelNumber(0);
  static bool firstcall(true);

  //  useful variables for output

  std::ostringstream oss;
  string msg, modelName;

  // check for bad value of imodel

  if ( modelNumber < 1 || modelNumber > 3 ) {
    msg = "Unknown value of modelNumber in eqpair code !";
    xs_write(const_cast<char*>(msg.c_str()),5);
    return;
  }

  if ( firstcall ) {
    msg = "      ";
    xs_write(const_cast<char*>(msg.c_str()),5);
    if ( modelNumber == 1 ) {
      modelName = "EQPAIR";
    } else if ( modelNumber == 2 ) {
      modelName = "EQTHERM";
    } else if ( modelNumber == 3 ) {
      modelName = "COMPTH";
    }
    msg = "          " + modelName + " V1.10";
    xs_write(const_cast<char*>(msg.c_str()),5);
    msg = "Uses ireflct for Compton reflection and rdblur for rotational blurring";
    xs_write(const_cast<char*>(msg.c_str()),5);
  }

  // check whether energies have changed

  const int nE = static_cast<int>(energyArray.size()) - 1;

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
  
  // This will be needed to determine whether to perform reflcalc below.
  // Notice that this is ALWAYS true when energiesChanged is true, and
  //   MAY remain true when energiesChanged is false.
  static bool engChangedSinceReflCalc = false;
  if (!engChangedSinceReflCalc)
     engChangedSinceReflCalc = energiesChanged; 

  // In order for the absorption and reflection calculations to be performed
  // correctly the energies on which the input spectrum are calculated should
  // run from 0.005 keV to 511.0/YMIN keV in the rest frame of the emitter.

  const Real YMIN = 0.03333;

  int NewBinsHigh = 0;
  int NewBinsLow = 0;
  Real dellogEHigh = 0.01;
  Real dellogELow = 0.01;

  Real zshift = 1.0 + params[19];
  Real RestEmin = 0.005/zshift;
  Real RestEmax = 511.0/YMIN/zshift;

  if ( energyArray[0] > RestEmin ) NewBinsLow = (int) ((log10(energyArray[0]) - log10(RestEmin)) * 100.0);

  if ( energyArray[nE] < RestEmax ) NewBinsHigh = (int) ((log10(RestEmax)-log10(energyArray[nE])) * 100.0);

  int NewBins = NewBinsLow + NewBinsHigh;
  int nEtot = nE + NewBins;

  // instantiate extended arrays

  RealArray ExtendEnergyArray(nEtot+1);
  RealArray ExtendFlux(nEtot);
  RealArray ExtendFluxErr(nEtot);

  // and set extended energy array

  for (int i=0; i<NewBinsLow; i++) ExtendEnergyArray[i] = pow(10.0,log10(RestEmin)+i*dellogELow);
  for (int i=NewBinsLow; i<=nE+NewBinsLow; i++) ExtendEnergyArray[i] = energyArray[i-NewBinsLow];
  for (int i=nE+NewBinsLow+1; i<(int)ExtendEnergyArray.size(); i++) 
    ExtendEnergyArray[i] = pow(10.0,log10(ExtendEnergyArray[nE+NewBinsLow])+(i-nE-NewBinsLow)*dellogEHigh);

// check for negative second parameter to indicate that only the black body
// component should be output

  bool showbb(false);
  if ( params[1] < 0.0 ) showbb  = true;


// find out whether to recalculate the source spectrum 

  bool redopaircalc(false);
  // jp checking only for the first 10 parameters (forget about accuracy)
  for (size_t i=0; i<10; i++) {
    if ( params[i] != saveParams[i] ) redopaircalc = true;
  }
  // Check also R_in, whether the actual model or the redshift have changed
  if (params[17] != saveParams[17] || modelNumber != saveModelNumber || 
      params[19] != saveParams[19] ) redopaircalc = true;

  if (redopaircalc) {

    msg = "Recalculating pairs...";
    xs_write(const_cast<char*>(msg.c_str()),20);
        
    // this is the routine that does most of the work then passes results back
    // need to convert params, compspect and blbdspect to float arrays then
    // copy the results back into the RealArrays.

    float pars[NPARAM];
    double spcomp[NGAMBIN], spblbd[NGAMBIN], energy[NGAMBIN];
    for (size_t i=0; i<(size_t)NPARAM; i++) pars[i] = params[i];
    for (size_t i=0; i<(size_t)NGAMBIN; i++) {
      spcomp[i] = compspect[i];
      spblbd[i] = blbdspect[i];
      energy[i] = wkev[i];
    }

    // pairiter needs slightly modified parameters

    if ( pars[1] < 0.0 ) pars[1] *= -1;

    PAIRITER(pars, NPARAM, modelNumber, spcomp, spblbd, energy);

    for (size_t i=0; i<(size_t)NGAMBIN; i++) {
      wkev[i] = energy[i];
      compspect[i] = spcomp[i];
      blbdspect[i] = spblbd[i];

    }

  }

  // now decide whether to recalculate the reflection

  bool redoreflcalc = redopaircalc;
  if ( engChangedSinceReflCalc || params[10] != saveParams[10] ) 
     redoreflcalc = true;
  for (size_t i=12; i<19; i++) {
    if (params[i] != saveParams[i]) redoreflcalc = true;
  }

  // if necessary calculate the reflection. The static array ReflSpect holds
  // the output so we can reuse on subsequent invocations if possible

  static RealArray ReflSpect;
  Real RefScale(fabs(params[11]));

  redoreflcalc &= ((RefScale > 0.0 || firstcall) && !showbb);
  if (redoreflcalc) {
    ReflSpect.resize(nEtot);

    // set the parameters for the reflection model

    RealArray par_refl(7);

    par_refl[0] = -1.0;             // ensures only reflection returned
    par_refl[1] = params[19];        // Redshift
    par_refl[2] = params[13];        // Abundances heavier than He
    par_refl[3] = params[12];        // Fe abundance
    par_refl[4] = params[10];        // cosIncl
    par_refl[5] = params[14];        // Disk temperature
    par_refl[6] = params[15];        // Ionization parameter

    msg = "Recalculating reflection...";
    xs_write(const_cast<char*>(msg.c_str()),20);

    // place compspect into the ReflSpect array

    GLrebin(wkev, compspect, ExtendEnergyArray, ReflSpect, params[19], 0, 0);

    // calculate the reflection spectrum (only)

    doreflect(ExtendEnergyArray, par_refl, 1, ReflSpect, ExtendFluxErr, modelName, ExtendEnergyArray[nE+NewBinsLow]);

    // now follows rotational broadening in case param[16] is not -10. Uses the rdblur
    // convolution model which is based on diskline

    if ( params[16] != -10 ) {

      msg = "Recalculating blurring...";
      xs_write(const_cast<char*>(msg.c_str()),20);

      // set up the rdblur parameters

      RealArray par_blur(4);
      par_blur[0] = params[16];                       // Beta
      par_blur[1] = params[17];                       // Rin
      par_blur[2] = params[18];                       // Rout
      par_blur[3] = acos(params[10])*180.0/3.14159265;// Inclination

      rdblur(ExtendEnergyArray, par_blur, 1, ReflSpect, ExtendFluxErr, modelName);

    }
    engChangedSinceReflCalc = false;

  }

  // initialize the output flux array to the reflection array if appropriate
  // otherwise to zero

  if ( RefScale > 0.0 && !showbb ) {
    ExtendFlux = ReflSpect * RefScale;
  } else {
    ExtendFlux = 0.0;
  }

  // rebin the bb and original comptonized spectra and add to the output array as required
  // only need to map back into the section of the flux array which will be preserved

  if ( params[11] >= 0 && !showbb ) {
    GLrebin(wkev, compspect, ExtendEnergyArray, ExtendFlux, params[19], NewBinsLow, NewBinsHigh);
  }

  if ( params[11] >= 0 ) {
    GLrebin(wkev, blbdspect, ExtendEnergyArray, ExtendFlux, params[19], NewBinsLow, NewBinsHigh);
  }

  // load flux back into original array

  flux.resize(nE);
  for (int i=0; i<nE; i++) flux[i] = ExtendFlux[i+NewBinsLow];

  // set all the saved variables used to determine whether the heavy
  // duty calculations are required on the next invocation

  saveParams = params;
  saveModelNumber = modelNumber;
  firstcall = false;

}

// Routine to rebin from the 221 energy bins used internally to another array
// Uses Gauss-Legendre integration over the output energy bins
void GLrebin(const RealArray& wkev, const RealArray& inSpectrum, const RealArray& outEnergy, RealArray& outSpectrum, Real Redshift, int NewBinsLow, int NewBinsHigh)
{

  static bool firstcall(true);
  static RealArray xg(10);
  static RealArray wg(10);

  int Ne = outSpectrum.size();

  // set up the Gauss-Legendre integration coefficients

  if ( firstcall ) {
    gsl_integration_glfixed_table* gltable = gsl_integration_glfixed_table_alloc(10);
    
    for (size_t i=0; i<10; ++i)
    {
       gsl_integration_glfixed_point(0.0, 1.0, i, &xg[i], &wg[i], gltable); 
    }
    gsl_integration_glfixed_table_free(gltable);
    firstcall = false;
  }

  Real zshift(1.0+Redshift);

  Real temp;
  for (size_t it=NewBinsLow; it<(size_t)(Ne-NewBinsHigh); it++) {

    Real xlo = outEnergy[it];
    Real xhi = outEnergy[it+1];
    Real dx = xhi-xlo;

    // Gauss-Legendre integration on this energy bin

    Real sum = 0.0;
    for (size_t j=0; j<10; j++) {
      Real targ = (xg[j]*dx + xlo)*zshift;
      //--- inline expansion
      if ((targ<wkev[0])||(targ>wkev[NGAMBIN-1])) {
	temp  = 1e-30;
      } else {
	int i = int(log10(targ/511.0)*20.0+141.5) - 1;
	if ( i > NGAMBIN - 1 ) i = NGAMBIN - 1;
	if ( i < 0 ) i = 0;
	if (targ < wkev[i]) i -= 1;
	if ( i < 0 ) i = 0;
	int i2 = i + 1;
	if ( i2 > NGAMBIN-1 ) i2 = NGAMBIN-1;
	Real x0 = wkev[i];
	Real x1 = wkev[i2];

	Real y1 = inSpectrum[i2];
	Real y0 = inSpectrum[i];
	if (i != i2) {
	  temp= y0*exp(log(y1/y0)/log(x1/x0)*log(targ/x0));
	} else {
	  temp = y0;
	}
      }

      // sum onto current accumulation

      sum += temp*wg[j]*dx;

    }

    // add the comptonized spectrum into the outSpectrum array
            
    outSpectrum[it] += sum;
            
  }
}

