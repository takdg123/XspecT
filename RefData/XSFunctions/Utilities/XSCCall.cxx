

#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/XSModelFunction.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <memory>

template <> 
void XSCall<xsccCall>::operator() (const RealArray& energyArray, const RealArray& params, int spectrumNumber, 
                         RealArray& fluxArray, RealArray& fluxErrArray, const string& initString) const
{
        const size_t nArray(energyArray.size());
        XSutility::Carray<Real> convert;

        std::unique_ptr<Real[]> pEnergy(convert(energyArray));
        std::unique_ptr<Real[]> pParam(convert(params));
        Real* energy  = pEnergy.get();
        Real* param   = pParam.get();

        int nEar = nArray - 1; 

        std::unique_ptr<Real[]> pFlux(new Real[nEar]);
        std::unique_ptr<Real[]> pFluxErr(new Real[nEar]);
        Real* flux = pFlux.get();
        Real* fluxErr = pFluxErr.get();    

        if (fluxArray.size() == static_cast<size_t>(nEar))
        {
                std::copy(&fluxArray[0],&fluxArray[0] + nEar,&flux[0]);
        }
        else
        {
                memset(flux,0,nEar*sizeof(Real));
                fluxArray.resize(nEar);
        }
        if ( fluxErrArray.size() == static_cast<size_t>(nEar))
        {
                std::copy(&fluxErrArray[0],&fluxErrArray[0] + nEar,&fluxErr[0]);
        }
        else
        {
                memset(fluxErr,0,nEar*sizeof(Real));
        }

        (*m_generator)(energy,nEar,param, spectrumNumber, flux,fluxErr, initString.c_str());

        fluxArray = RealArray(flux,nEar);

        // keep testing while error array is zero. If it drops out
        // of the loop before the end, the error array is present and
        // needs to be set.
        int testErr(0);
        while ( testErr < nEar && fluxErr[testErr] == 0 ) ++testErr;  

        if (testErr < nEar)
        {
                fluxErrArray.resize(nEar);
                fluxErrArray = RealArray(fluxErr,nEar);
        }
        else 
        {
                fluxErrArray.resize(0); 
        }


}

template <> 
void XSCall<XSCCall>::operator() (const RealArray& energyArray, const RealArray& param, int spectrumNumber,
                         RealArray& fluxArray, RealArray& fluxErrArray, const string& initString) const
{

        (*m_generator)(energyArray,param, spectrumNumber, fluxArray,fluxErrArray,initString);       

}

template  <> 
void XSCall<XSMixCCall>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator,
                        const string& modelName) const
{
        (*m_generator)(energyArray, parameterValues, flux, fluxError, mixGenerator,
                        modelName);       

}

template  <> 
void XSCall<xsmixcall>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator,
                        const string& modelName) const
{

}

template <> 
void XSCall<MdefExpression>::operator() (const RealArray& energyArray, const RealArray& param, int spectrumNumber,
                         RealArray& fluxArray, RealArray& fluxErrArray, const string& initString) const
{

        (*m_generator).evaluate(energyArray, param, fluxArray,fluxErrArray);       

}

template <>
XSCall<MdefExpression>::~XSCall()
{
   delete m_generator;
}

template <>
XSCall<MdefExpression>::XSCall(const XSCall<MdefExpression> &right)
   : XSModelFunction(right)
{
   m_generator = right.m_generator->clone();
}
