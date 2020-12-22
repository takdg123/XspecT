 // Code for power-law model. Includes calcPowerLaw function that can be called by
// other functions to evaluate the power-law.
//
// Parameters :
//      alpha          power-law index (E^-alpha)
//      redshift       for zpow model
//
// Also checks for xset strings POW_EMIN and POW_EMAX. If these are set and not equal
// then the normalization is the flux over(POW_EMIN, POW_EMAX) keV in units of 
// 10^-12 ergs/cm^2/s. If POW_EMIN = POW_EMAX then the normalization is the flux
// density at POW_EMIN in units of micro-Janskys. If POW_EMIN and POW_EMAX are not
// then the normalization is the flux density at 1 keV in photons/cm^2/s/keV. For the
// zpow model all normalizations are in the source frame.
//

#include <XSFunctions/Utilities/funcType.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"
#include <cmath>
#include <iomanip>

// prototypes

Real calcPowerLaw(const RealArray& energyArray, const Real& index, bool isRenorm, RealArray& fluxArray); 
Real pegRenormalize(const RealArray& energyArray, RealArray& fluxArray);


//**************************************************************************
// Power-law model with no redshift parameter

void 
powerLaw (const RealArray& energyArray, 
                const RealArray& params, 
                int spectrumNumber,
                RealArray& fluxArray, 
                RealArray& fluxErrArray,
                const string& initString)
{

  calcPowerLaw(energyArray, params[0], true, fluxArray);
  fluxErrArray.resize(0);

  return;

}

//**************************************************************************
// Redshifted power-law. Note that normalization is in source frame (even if
// POW_EMIN and POW_EMAX have been applied).

void 
zpowerLaw (const RealArray& energyArray, 
                const RealArray& params, 
                int spectrumNumber,
                RealArray& fluxArray, 
                RealArray& fluxErrArray,
                const string& initString)
{

  const Real z(params[1]);
  const Real fac(1 + z);       

  const RealArray energy(energyArray*fac);

  calcPowerLaw(energy, params[0], true, fluxArray);
  fluxArray /= fac;

  fluxErrArray.resize(0);

}

//**************************************************************************

// Routine to actually calculate power-law. if isRenorm=true apply any renormalization factor
// using the POW_EMIN or POW_EMAX settings and return the factor.

Real calcPowerLaw(const RealArray& energyArray, const Real& index, bool isRenorm,
		  RealArray& fluxArray)
{
  using namespace std;
  const Real alpha ( 1 + index );
  const Real ep=1e2;

  size_t N(energyArray.size());
  fluxArray.resize(N-1);

  // note tolerance to avoid numerical problems

  if ( abs(alpha) < 1e-10 ) {
    RealArray logE(std::log(energyArray));
    Real first (logE[0]);
    for (size_t i = 1; i < N; ++i) {
      Real second(logE[i]);
      fluxArray[i - 1] = second - first;
      first = second;
    }               
  } else {
    const Real alphani (1./alpha);
    Real first (alphani*pow(energyArray[0],alpha));
    for (size_t i = 1; i < N; ++i) {
      Real second(alphani*pow(energyArray[i],alpha));
      fluxArray[i - 1] = (second - first)*pow(ep, -index);
      first =  second;
    }               
  }

  if ( isRenorm ) {

    Real RenormFactor = pegRenormalize(energyArray, fluxArray);
    return(RenormFactor);

  } else {

    return(1.0);
  }

}

//**************************************************************************

// Renormalize power-law based on any POW_EMIN and POW_EMAX settings
// Returns renormalization factor.

Real pegRenormalize(const RealArray& energyArray, RealArray& fluxArray)
{
  using namespace std;

  // check for xset POW_EMIN and POW_EMAX values

  string pname("POW_EMIN");
  string eminstr(FunctionUtility::getModelString(pname));
  pname = "POW_EMAX";
  string emaxstr(FunctionUtility::getModelString(pname));

  if ( !eminstr.length() || !emaxstr.length() || 
       eminstr == FunctionUtility::NOT_A_KEY() ||
       emaxstr == FunctionUtility::NOT_A_KEY() ) return(1.0);

  istringstream eminstream(eminstr);
  Real Emin = -1.0;
  if (!(eminstream >> Emin) || !eminstream.eof()) {
    string message("Failed to read value from POW_EMIN");
    xs_write(const_cast<char*>(message.c_str()),10);
  }

  istringstream emaxstream(emaxstr);
  Real Emax = -1.0;
  if (!(emaxstream >> Emax) || !emaxstream.eof())  {
    string message("Failed to read value from POW_EMAX");
    xs_write(const_cast<char*>(message.c_str()),10);
  }

  // if either Emin or Emax <= 0 then use standard norm

  if ( Emin <= 0.0 || Emax <= 0.0 ) return(1.0);

 // if Emin != Emax then normalization is integrated flux

  if ( Emin != Emax ) {

    // use the integration routine in Numerics to calculate the flux

    pair<Real,Real> flux(0.,0.);
    flux = Numerics::integrationKernel(energyArray, fluxArray, Emin, Emax);

    // the second member of the pair is now the flux in units of erg/cm^2/s
    // so convert to 1e-12 erg/cm^2/s

    flux.second *= 1.0e12;

    // renormalize the power-law model

    if ( flux.second > 0 ) {

      fluxArray /= flux.second;
      return(1.0/flux.second);

    } else {

      return(1.0);

    }

  }

  // if Emin == Emax then normalization is the flux density in micro-Jy at Emin

  else {

    int i;
    XSutility::find(energyArray, Emin, i);

    if ( i == -1 || i == (int)(energyArray.size()-1) ) return(1.0);

    // flux density in keV (ph/cm^2/s/keV)

    Real fluxdens(fluxArray[i]*sqrt(energyArray[i]*energyArray[i+1])/(energyArray[i+1]-energyArray[i]));

    // convert to micro-Jy 

    fluxdens *= Numerics::KEVTOJY * 1e6 / Numerics::KEVTOHZ;

    // renormalize the power-law model

    fluxArray /= fluxdens;

    return(1.0/fluxdens);

  }

}

