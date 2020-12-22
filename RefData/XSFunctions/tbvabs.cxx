#include <xsTypes.h>
#include <functionMap.h>
#include <stlToCArrays.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <memory>
#include <cfortran.h>
#include <cstring>

// Functions called FROM here to Fortran:

PROTOCCALLSFSUB10(DOTBVABS,dotbvabs,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV)
#define DOTBVABS(ear,ne,param,ifl,photar,photer,siggas,siggrains,sigmol,sigavg) \
           CCALLSFSUB10(DOTBVABS,dotbvabs,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV,FLOATV, \
             ear,ne,param,ifl,photar,photer,siggas,siggrains,sigmol,sigavg)

// definition of the tbnew function from tbnew.c

extern "C" {
void tbnew(const double *ear, int ne, const double *param, int ifl, double *photar,
		 double *photer, const char *init);
}

// other routines defined in this file

void tbdefaults(RealArray& params);
double fggrain(const char* element);


void tbvabs(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // Wrapper routine for the Tuebingen-Boulder ISM code. Grabs the dynamic 
   // memory required.

   // Translated from tbvabs.f 03/09 (CG)
   // Now either uses the old version (1.0) or the new version (2.0)
   // depending on the xset TBABSVERSION value

   // check for any TBABSVERSION number set
   string version = "2.0";
   string pname = "TBABSVERSION";
   string pvalue = FunctionUtility::getModelString(pname);
   if ( pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY() ) {
     version = pvalue;
   }

   if ( version.substr(0,1) == "1" ) {

     // old version of tbabs

     const int ne = static_cast<int>(energyArray.size() - 1);

     float *ear=0, *pars=0, *photar=0, *photer=0;
     XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
					  ear, pars, photar, photer);
     std::unique_ptr<float[]> apEar(ear);
     std::unique_ptr<float[]> apPars(pars);
     std::unique_ptr<float[]> apPhotar(photar);
     std::unique_ptr<float[]> apPhoter(photer);

     std::unique_ptr<float[]> apSiggas(new float[ne]);
     std::unique_ptr<float[]> apSiggrains(new float[ne]);
     std::unique_ptr<float[]> apSigmol(new float[ne]);
     std::unique_ptr<float[]> apSigavg(new float[ne]);

     DOTBVABS(ear, ne, pars, spectrumNumber, photar, photer, apSiggas.get(),
	      apSiggrains.get(), apSigmol.get(), apSigavg.get());

     XSFunctions::floatFluxToStl<float>(photar, photer, ne, false, flux, fluxErr);       

   } else if ( version.substr(0,1) == "2" ) {

     // new version of tbabs

     const int ne = static_cast<int>(energyArray.size() - 1);

     double *ear=0, *pars=0, *photar=0, *photer=0;
     XSFunctions::stlToFloatArrays<double>(energyArray, params, flux, fluxErr,
					  ear, pars, photar, photer);
     std::unique_ptr<double[]> apEar(ear);
     std::unique_ptr<double[]> apPars(pars);
     std::unique_ptr<double[]> apPhotar(photar);
     std::unique_ptr<double[]> apPhoter(photer);

     tbnew(ear, ne, pars, spectrumNumber, photar, photer, 
	   initString.c_str());

     XSFunctions::floatFluxToStl<double>(photar, photer, ne, false, flux, fluxErr);       

   } else {

     string msg = version + " is not a valid tbabs version number";
     FunctionUtility::xsWrite(msg,5);

   }

   // convert arrays back to C++ like


}

void tbdefaults(RealArray& params)
{
  // Sets the 42 parameters used by tbvabs to their default values. This is useful
  // for models which call tbvabs but only use a subset of the parameters.

   const size_t NELE = 18;
   const char *celts[] = {"H ", "He", "C ", "N ", "O ", "Ne", "Na", 
			  "Mg", "Al", "Si", "S ", "Cl", "Ar", "Ca", 
			  "Cr", "Fe", "Co", "Ni"};


  params.resize(42);

  // nH in units of 1e22 then the abundances wrt Solar of He, C, N, O, Ne, Na,
  // Mg, Al, Si, S, Cl, Ar, Ca, Cr, Fe, Co, Ni

  for (size_t i=0; i<18; i++) params[i] = 1.0;

  // fraction of H existing in molecular form
  params[18] = 0.2;

  // density of dust in g/cm^3
  params[19] = 1.0;

  // minimum and maximum dust grain sizes in mum
  params[20] = 0.025;
  params[21] = 0.25;

  // power-law index for dust size distribution (implicit minus-sign!)
  params[22] = 3.5;

  // depletion factors (ratio of GAS to total ISM abundance)
  for (size_t i=0; i<NELE; i++) params[23+i] = fggrain(celts[i]);

  // redshift
  params[41] = 0.0;
 
  return;
}


double fggrain(const char* element)
{
   //     Returns the grain composition for the grain models summarized
   //     by Wilms, Allen, McCray, 2000, ApJ, in press
   //     The value returned is the ratio of the GAS abundance to the total
   //     ISM abundance (the "depletion factor" commonly used by the ISM
   //     people).
   //
   // Arguments :
   //      element    C    i: The standard chemical element abbreviation
   //      fggrain    R    r: the deplection factor
   //
   // NOTE: This subroutine knows about more elements than the standard
   //     ones used in XSPEC. Therefore, NELTS=21!
   //
   // Version 1.0: Joern Wilms, IAA Tuebingen, Astronomie
   //    wilms@astro.uni-tuebingen.de

   // Translated from tbvabs.f 03/09 (CG)

   const int NELTS = 21;
   const char *elts[] =   {"H ", "He", "C ", "N ", "O ", "Ne", "Na", 
               "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "Ca", "Ti",
               "Cr", "Mn", "Fe", "Co", "Ni"};
   const double fdep[] =  {1.0,  1.0,  0.5,  1.0,  0.6,   1.0,  0.25, 
               0.2,  0.02, 0.1,  0.6,  0.6,  0.5,  1.0,  0.003, 0.002,
               0.03, 0.07, 0.3,  0.05, 0.04};

   double val = 0.0;            
   for (int i=0; i<NELTS; ++i)
   {
      if (strcmp(elts[i], element) == 0)
      {
         val = fdep[i];
         break;
      }
   }
   return val;
}
