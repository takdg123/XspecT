// Broken power law model
//    Parameters are :
//       1      photon power-law index ( E < BreakE )
//       2      BreakE  (keV)
//       3      photon power-law index ( E > BreakE )
// model form N(E) = (E**-par1) [E<par2]
//                 = (par2)**(par3-par1) * E**-par3 [E>par2]
//
// and redshifted version
//    Parameters are :
//       1      photon power-law index ( E < BreakE )
//       2      BreakE  (emitted frame: keV)
//       3      photon power-law index ( E > BreakE )
//       4      redshift
// model form N(E) = ([E(1+z)]**-par1) [E(1+z)<par2]
//                 = (par2)**(par3-par1) * [E(1+z)]**-par3 [E(1+z)>par2]

#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

// prototype for routine in this file

Real calcBrokenPowerLaw(const RealArray& energyArray, const RealArray& params, 
			bool isRenorm, RealArray& fluxArray); 

// prototypes for routines from powerLaw.cxx

Real calcPowerLaw(const RealArray& energyArray, const Real& index, bool isRenorm, 
		  RealArray& fluxArray); 
Real pegRenormalize(const RealArray& energyArray, RealArray& fluxArray);

//****************************************************************************

void 
brokenPowerLaw (const RealArray& energyArray, 
                const RealArray& params, 
                int spectrumNumber,
                RealArray& fluxArray, 
                RealArray& fluxErrArray,
                const string& initString)
{
  fluxErrArray.resize(0);

  calcBrokenPowerLaw(energyArray, params, true, fluxArray); 

  return;
}

//****************************************************************************
void zBrokenPowerLaw (const RealArray& energyArray, 
                      const RealArray& params, 
                      int spectrumNumber,
                      RealArray& fluxArray, 
                      RealArray& fluxErrArray,
                      const string& initString)
{
  fluxErrArray.resize(0);

  Real zfactor = (1.0+params[3]);
  RealArray bparams(3);
  for (size_t i=0; i<3; i++) bparams[i] = params[i];

  const RealArray energy(energyArray*zfactor);
  
  calcBrokenPowerLaw(energy, bparams, true, fluxArray);
  fluxArray /= zfactor;

  return;
}

//****************************************************************************

Real calcBrokenPowerLaw(const RealArray& energyArray, const RealArray& params, bool isRenorm, RealArray& fluxArray)
{
  const Real breakEnergy(params[1]);
  const int Ne(energyArray.size()-1);
  fluxArray.resize(Ne);

  // normalization factor to apply to spectrum above the break.

  const Real normFactor(pow(breakEnergy,params[2]-params[0]));

  // handle the two special cases of break energy below first energy bin or
  // break energy above last energy bin

  if ( breakEnergy < energyArray[0] ) {

    calcPowerLaw(energyArray, params[2], false, fluxArray);
    fluxArray *= normFactor;

  } else if ( breakEnergy > energyArray[Ne] ) {

    calcPowerLaw(energyArray, params[0], false, fluxArray);

  } else {

    // find the energy bin containing the break energy using a binary search

    int breakBin;
    int low = 0;
    int high = Ne;
    while ( (high-low) > 1 ) {
      breakBin = ( low + high) / 2;
      if ( breakEnergy > energyArray[breakBin] ) {
	low = breakBin;
      } else {
	high = breakBin;
      }
    }
    breakBin = low;

    // evaluate power-law below the break energy

    RealArray eArray(breakBin+2);
    RealArray fArray(breakBin+1);

    for (int i=0; i<=breakBin; i++) eArray[i] = energyArray[i];
    eArray[breakBin+1] = breakEnergy;

    calcPowerLaw(eArray, params[0], false, fArray);

    for (int i=0; i<=breakBin; i++) fluxArray[i] = fArray[i];

    // evaluate power-law above the break energy

    const int Neup(Ne - breakBin);
    eArray.resize(Neup+1);
    fArray.resize(Neup);

    for (int i=1; i<=Neup; i++) eArray[i] = energyArray[i+breakBin];
    eArray[0] = breakEnergy;

    calcPowerLaw(eArray, params[2], false, fArray);
    fArray *= normFactor;

    for (int i=breakBin+1; i<Ne; i++) fluxArray[i] = fArray[i-breakBin];
    fluxArray[breakBin] += fArray[0];

  }

  // Perform any required final renormalization of the entire spectrum

  if ( isRenorm ) {

    Real renorm = pegRenormalize(energyArray, fluxArray);
    return(renorm);

  } else {

    return(1.0);

  }

}
