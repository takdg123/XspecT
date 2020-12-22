#include <xsTypes.h>
#include <stlToCArrays.h>
#include <functionMap.h>
#include <memory>

extern "C" void runkbb_(float* ear, int& nE, float& eta, float& astar, 
		float& theta, float& mbh, float& mdd, float& dbh, 
		float& fcol, int& rflag, int& lflag, float& zbh, 
                float* photar, float* fluxE, int* nex1, int* nex2);

void kerrbb(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // Memory allocation wrapper function for Li-Xin Li's runkbb routine.

   int nE = static_cast<int>(energyArray.size()) - 1;

   // The memory for temporary arrays

   std::unique_ptr<float[]> apFluxE(new float[nE*2]);
   std::unique_ptr<int[]> apNex1(new int[nE*2]);
   std::unique_ptr<int[]> apNex2(new int[nE*2]);

   float *ear=0, *pars=0, *photar=0, *photer=0;
   XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
           ear, pars, photar, photer);
   std::unique_ptr<float[]> apEar(ear);
   std::unique_ptr<float[]> apPars(pars);
   std::unique_ptr<float[]> apPhotar(photar);
   std::unique_ptr<float[]> apPhoter(photer);

   // Call the main routine

   float eta = pars[0];
   float astar = pars[1];
   float theta = pars[2];
   float mbh = pars[3];
   float mdd = pars[4];
   float dbh = pars[5];
   float fcol = pars[6];
   int rflag = (int)round(pars[7]);
   int lflag = (int)round(pars[8]);
   float zbh = 0.0;

   runkbb_(ear, nE, eta, astar, theta, mbh, mdd, dbh, fcol, rflag, lflag,
	   zbh, photar, apFluxE.get(), apNex1.get(), apNex2.get());
   XSFunctions::floatFluxToStl<float>(photar, photer, nE, false, flux, fluxErr);

   // no flux errors associated with this model
   fluxErr = 0.0;
}
