// Wrapper for dospin.f

#include <xsTypes.h>
#include <functionMap.h>
#include <stlToCArrays.h>
#include <memory>
#include <cfortran.h>

// Functions called FROM here to Fortran:

PROTOCCALLSFSUB6(DOSPIN,dospin,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
#define DOSPIN(ear,ne,param,ifl,photar,photer) \
           CCALLSFSUB6(DOSPIN,dospin,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV, \
             ear,ne,param,ifl,photar,photer)



void spin(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // Memory allocation wrapper for dospin

   int ne = static_cast<int>(energyArray.size()) - 1;

   float *ear=0, *pars=0, *photar=0, *photer=0;
   XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
           ear, pars, photar, photer);
   std::unique_ptr<float[]> apEar(ear);
   std::unique_ptr<float[]> apPars(pars);
   std::unique_ptr<float[]> apPhotar(photar);
   std::unique_ptr<float[]> apPhoter(photer);

   DOSPIN(ear, ne, pars, spectrumNumber, photar, photer);   

   XSFunctions::floatFluxToStl<float>(photar, photer, ne, false, flux, fluxErr);       
}
