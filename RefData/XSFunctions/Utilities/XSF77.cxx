

#include <XSFunctions/Utilities/XSModelFunction.h>
#include <XSUtil/Utils/XSutility.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <memory>


template <> 
void XSCall<xsf77Call>::operator() (const RealArray& energyArray, const RealArray& params, int spectrumNumber, 
                         RealArray& fluxArray, RealArray& fluxErrArray, const string& initString) const
{
  const size_t nArray(energyArray.size());
  const size_t nParams = params.size();
  XSutility::Carray<Real> convert;

  std::unique_ptr<Real[]> enp(convert(energyArray));
  std::unique_ptr<Real[]> parp(convert(params));
  Real* en = enp.get();
  Real* par = parp.get();

  // use auto_ptr rather than plain C arrays to protect from memory
  // leaks. Also, delegate array copying to std::copy to take advantage
  // of potential compiler optimizations.

  int nEar = nArray - 1; 

  std::unique_ptr<float[]> enf(new float[nArray]);
  std::unique_ptr<float[]> parf(new float[params.size()]);
  std::unique_ptr<float[]> flux(new float[nEar]);
  std::unique_ptr<float[]> fluxErr(new float[nEar]);



  float* ear = enf.get();
  float* parr = parf.get();
  float* fluxg = flux.get();
  float* fluxErrg = fluxErr.get();

  std::copy(&en[0],&en[nArray],&ear[0]);
  std::copy(&par[0],&par[nParams],&parr[0]);

  if (fluxArray.size() == static_cast<size_t>(nEar))
  {
          std::copy(&fluxArray[0],&fluxArray[0] + nEar,&fluxg[0]);
  }
  else
  {
          memset(fluxg,0,nEar*sizeof(float));
  }
  if ( fluxErrArray.size() == static_cast<size_t>(nEar))
  {
        std::copy(&fluxErrArray[0],&fluxErrArray[0] + nEar,&fluxErrg[0]);
  }
  else
  {
        memset(fluxErrg,0,nEar*sizeof(float));
  }  
  //memset(fluxg,0,nEar*sizeof(float));

  (*m_generator)(ear,nEar,parr,spectrumNumber,fluxg,fluxErrg);

  std::unique_ptr<Real[]> fluxd(new Real[nEar]);
  Real* fluxrd = fluxd.get();
  std::copy(&fluxg[0],&fluxg[nEar],&fluxrd[0]);
  fluxArray.resize(nEar,0.0);
  fluxArray = RealArray(fluxrd,nEar);

  // keep testing while error array is zero. If it drops out
  // of the loop before the end, the error array is present and
  // needs to be set.
  int testErr(0);
  while ( testErr < nEar && fluxErrg[testErr] == 0 ) ++testErr;

  if (testErr < nEar)
  {
        std::unique_ptr<Real[]> fluxErrd(new Real[nEar]);
        Real* fluxErrdr = fluxErrd.get();
        std::copy(&fluxErrg[0],&fluxErrg[nEar],&fluxErrdr[0]);
        fluxErrArray.resize(nEar,0.0);
        fluxErrArray = RealArray(fluxErrdr,nEar);
  }
  else 
  {
        fluxErrArray.resize(0); 
  }


}

template <> 
void XSCall<xsF77Call>::operator() (const RealArray& energyArray, const RealArray& params, int spectrumNumber,
                         RealArray& fluxArray, RealArray& fluxErrArray, const string& initString) const
{
  const size_t nArray(energyArray.size());
  int nEar = nArray - 1;   
  XSutility::Carray<Real> convert;

  // use auto_ptr rather than plain C arrays to protect from memory
  // leaks. 
  std::unique_ptr<Real[]> enp(convert(energyArray));
  std::unique_ptr<Real[]> parp(convert(params));
  std::unique_ptr<Real[]> flux(new Real[nEar]);
  std::unique_ptr<Real[]> fluxErr(new Real[nEar]);

  Real* ear = enp.get();
  Real* parr = parp.get();
  Real* fluxg = flux.get();
  Real* fluxErrg = fluxErr.get();

  if (fluxArray.size() == static_cast<size_t>(nEar))
  {
          std::copy(&fluxArray[0],&fluxArray[0] + nEar,&fluxg[0]);
  }
  else
  {
          memset(fluxg,0,nEar*sizeof(Real));
  }
  if ( fluxErrArray.size() == static_cast<size_t>(nEar))
  {
        std::copy(&fluxErrArray[0],&fluxErrArray[0] + nEar,&fluxErrg[0]);
  }
  else
  {
        memset(fluxErrg,0,nEar*sizeof(Real));
  }
  //memset(fluxg,0,nEar*sizeof(float));

  // double precision fortran call.
  (*m_generator)(ear,nEar,parr,spectrumNumber,fluxg,fluxErrg);

  fluxArray.resize(nEar,0.0);
  fluxArray = RealArray(fluxg,nEar);

  // keep testing while error array is zero. If it drops out
  // of the loop before the end, the error array is present and
  // needs to be set.
  int testErr(0);
  while ( testErr < nEar && fluxErrg[testErr] == 0 ) ++testErr;

  if (testErr < nEar)
  {
        fluxErrArray.resize(nEar,0.0);
        fluxErrArray = RealArray(fluxErrg,nEar);
  }
  else 
  {
        fluxErrArray.resize(0); 
  }


}

