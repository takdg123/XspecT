// Broken power law model with 2 pivot points
//    Parameters are :
//       1      photon power-law index ( E < BreakE1 )
//       2      BreakE1  (keV)
//       3      photon power-law index ( BreakE1 < E < BreakE2 )
//       4      BreakE2  (keV)
//       5      photon power-law index ( E > BreakE2 )
// model form N(E) = (E**-par1) [E<par2]
//                 = (par2)**(par3-par1) * E**-par3 [par2<E<par4]
//                 = (par2)**(par3-par1) * (par4)**(par5-par3) * E**-par5 [E>par4]


#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

// prototypes for routines from brokenPowerLaw.cxx and powerLaw.cxx

Real calcBrokenPowerLaw(const RealArray& energyArray, const RealArray& params, 
			bool isRenorm, RealArray& fluxArray); 

Real calcPowerLaw(const RealArray& energyArray, const Real& index, bool isRenorm, 
		  RealArray& fluxArray); 
Real pegRenormalize(const RealArray& energyArray, RealArray& fluxArray);



void 
broken2PowerLaw (const RealArray& energyArray, 
                const RealArray& params, 
                int spectrumNumber,
                RealArray& fluxArray, 
                RealArray& fluxErrArray,
                const string& initString)
{
  const Real breakEnergy(params[1]);
  const int Ne(energyArray.size()-1);
  fluxArray.resize(Ne);
  fluxErrArray.resize(0);

  // normalization factor to apply to spectrum above the first break.

  const Real normFactor(pow(breakEnergy,params[2]-params[0]));

  // find the energy bin containing the first break energy using a binary search
  // handle the two special cases of break energy below first energy bin or
  // break energy above last energy bin

  int breakBin;
  if ( breakEnergy < energyArray[0] ) {
    RealArray pArray(3);
    for (int i=0; i<3; i++) pArray[i] = params[i+2];
    Real renorm = calcBrokenPowerLaw(energyArray, pArray, true, fluxArray);
    if ( renorm == 1.0 ) fluxArray *= normFactor;
    return;
  } else if ( breakEnergy > energyArray[Ne] ) {
    calcPowerLaw(energyArray, params[0], true, fluxArray);
    return;
  } else {
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
  }

  // evaluate power-law below the first break energy

  RealArray eArray(breakBin+2);
  RealArray fArray(breakBin+1);

  for (int i=0; i<=breakBin; i++) eArray[i] = energyArray[i];
  eArray[breakBin+1] = breakEnergy;

  calcPowerLaw(eArray, params[0], false, fArray);

  for (int i=0; i<=breakBin; i++) fluxArray[i] = fArray[i];

  // evaluate above the first break energy treating this part of the spectrum
  // as a broken power-law with a single pivot point

  const int Neup(Ne - breakBin);
  eArray.resize(Neup+1);
  fArray.resize(Neup);
  RealArray pArray(3);

  for (int i=1; i<=Neup; i++) eArray[i] = energyArray[i+breakBin];
  eArray[0] = breakEnergy;

  for (int i=0; i<3; i++) pArray[i] = params[i+2];

  calcBrokenPowerLaw(eArray, pArray, false, fArray);
  fArray *= normFactor;

  for (int i=breakBin+1; i<Ne; i++) fluxArray[i] = fArray[i-breakBin];
  fluxArray[breakBin] += fArray[0];

  pegRenormalize(energyArray, fluxArray);

  return;
}
