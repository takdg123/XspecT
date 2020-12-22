#include "IonizedOpacity.h"

// prototype for routine from powerLaw.cxx

Real calcPowerLaw(const RealArray& energyArray, const Real& gamma, bool isRenorm,
		  RealArray& flux);

// Main model routine

void xsabsori(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{
  // Absorption by an ionized medium.
  //     PARAMETERS
  //               1       gamma - photon index of power spectrum
  //               2       nh - Hydrogen column of the absorper in units of 10^22
  //               3       temp  - absorber temperature in K
  //               4       xi    - absorber ionization state = L/nR^2
  //               5       z     - redshift
  //               6       ab0   - iron abundance relative to the 
  //                               solar iron abundance
  //       algorithm:
  //               a(x) = a(x)*exp(-Nh*sigma(x))
  //

  Real Gamma(params[0]);
  Real nH(params[1]);
  Real Temp(params[2]);
  Real Xi(params[3]);
  Real zshift(1.0+params[4]);
  Real FeAbund(log10(params[5]));

  size_t N(energyArray.size()-1);

  static IonizedOpacity warmabs;

  // set up an input power-law spectrum over the energy range on which the
  // opacity data is tabulated

  size_t powN(5000);
  RealArray powE(powN+1);
  RealArray powF(powN);
  Real dellog((log(20.0)-log(0.005))/(powN+1));

  for (size_t i=0; i<powN+1; i++) powE[i] = pow(10.0,log(0.005)+i*dellog);

  calcPowerLaw(powE, Gamma, false, powF);

  // convert to units required for IonizedOpacity object

  RealArray X(powN);
  RealArray Spinc(powN);

  Real factor = 1.0/(2.0*511.0);
  for (size_t i=0; i<powN; i++) X[i] = (powE[i+1]+powE[i]) * factor;

  for (size_t i=0; i<powN; ++i) {
    if (powE[i+1] != powE[i]) {
      Real avgE = (powE[i+1]+powE[i])/2.;
      Spinc[i] = powF[i]*avgE*avgE/(powE[i+1]-powE[i]);
    }
  }

  // set up opacities

  warmabs.Setup(Xi, Temp, X, Spinc);

  // now make X the input redshifted energies
  // and calculate the opacities

  X.resize(N);
  factor = zshift/(2.0*511.0);
  for (size_t i=0; i<N; i++) X[i] = (energyArray[i+1]+energyArray[i]) * factor;

  RealArray opacity(N);
  Real Abund(0.0);
  warmabs.Get(X, Abund, FeAbund, true, opacity);

  // calculate absorption fraction

  flux.resize(N);
  flux = exp(-nH * opacity);

  return;

}
