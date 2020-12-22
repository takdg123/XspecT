#include <xsTypes.h>
#include <functionMap.h>
#include <stlToCArrays.h>
#include <memory>
#include <cfortran.h>

// Functions called FROM here to Fortran:

PROTOCCALLSFSUB6(DONTHCOMP,donthcomp,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
#define DONTHCOMP(ear,ne,param,ifl,photar,prim) \
           CCALLSFSUB6(DONTHCOMP,donthcomp,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV,\
             ear,ne,param,ifl,photar,prim)


void nthcomp(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // Actual work is done in donthcomp.f

   int ne = static_cast<int>(energyArray.size()) - 1;

   float *ear=0, *pars=0, *photar=0, *photer=0;
   XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
           ear, pars, photar, photer);
   std::unique_ptr<float[]> apEar(ear);
   std::unique_ptr<float[]> apPars(pars);
   std::unique_ptr<float[]> apPhotar(photar);
   std::unique_ptr<float[]> apPhoter(photer);

   DONTHCOMP(ear, ne, pars, spectrumNumber, photar, photer);

   XSFunctions::floatFluxToStl<float>(photar, photer, ne, false, flux, fluxErr);       

}
