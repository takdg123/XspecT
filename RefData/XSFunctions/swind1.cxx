#include <xsTypes.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
// debug
#include <iostream>

// prototype from tableInterpolate.cxx

void tableInterpolate(const RealArray& energyArray, const RealArray& params, 
		      string fileName, int spectrumNumber, RealArray& fluxArray, 
		      RealArray& fluxErrArray, const string& initString, const string& tableType, const bool readFull);


void swind1(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& trans, RealArray& transErr, 
	    const string& initString)
{
  // this model does not calculate errors
  transErr.resize(0);

  string DirName = FunctionUtility::modelDataPath();
  string fileName = DirName + "swind1_mtable.fits";

  RealArray eparams(3);
  eparams[0] = params[0]*1.e22;      // nh in units of 1e22 for param(1)      
  eparams[1] = params[1];            // log xi
  eparams[2] = params[3];            // redshift

  Real zfact(1.0+params[3]);

  // Gaussian sigma smearing - renormalised to 1 KeV from 6. See gsmooth.
  RealArray sparams(2);
  sparams[0] = 6.0 * params[2];
  // Sigma Variation law with eneregy - fixed to produce constant deltaE/E
  sparams[1] = 1.0;

  //Rebinning Grid size - fix at 1000
  size_t nex = 1000;

  // Extend energy array

  Real emin = log10(0.1 * energyArray[0]);
  Real emax = log10(10.0 * energyArray[energyArray.size()-1]);
  Real de = (emax - emin) / (Real)nex;
      
  RealArray e(nex+1);
  for (size_t i=0; i<nex+1; i++) e[i] = pow(10.0,emin+i*de);

  RealArray p, perr;
  tableInterpolate(e, eparams, fileName, spectrumNumber, p, perr, " ", "mul", true);
  if ( p.size() == 0 ) {
    FunctionUtility::xsWrite("Failed to read from "+fileName, 5);
    return;
  }

  // fill in below/above the energy limits of 0.1 and 18.0 keV
  // this is zero emission for the emitter, or unity (zero absorption) for
  // the absorber
  // and scale column linearly from 6e22

  Real elo = 0.1 / zfact;
  Real ehi = 18.0 / zfact;

  for (size_t i=0; i<nex; i++) {
    if ( e[i] >= ehi || e[i] <= elo) p[i] = 1.0;
  }

  // gaussian smoothing

  gsmooth(e, sparams, spectrumNumber, p, perr, " ");

  // rebinning
    
  size_t ix = 0;
  trans.resize(energyArray.size()-1);
  for (size_t i=0; i<trans.size(); i++) {
    trans[i] = 0.0;
    Real emin = energyArray[i];
    Real emax = energyArray[i+1];
    while( e[ix+1] < emin) ix++;
    while( e[ix] < emax) {
      Real e1(e[ix]);
      if ( emin > e1 ) e1 = emin;
      Real e2(e[ix+1]);
      if ( emax < e2 ) e2 = emax;
      trans[i] += p[ix] * (e2 - e1);
      ix++;
    }
    ix--;
    trans[i] /= (emax - emin);
  }

  return;
}



